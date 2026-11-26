#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""CRC/irinotecan progressive-disease driver for the consensus nested Cox evaluator.

Compatibility target
--------------------
Python          3.6.8
lifelines       0.26.4
scikit-learn    0.24.2
numpy           1.19.5
pandas          1.1.5
scipy           1.5.4
matplotlib       3.3.4 (recommended for Python 3.6.8)

This driver can load a CRC/irinotecan progressive-disease cohort (or generate
the existing synthetic fallback) and runs two analyses against the *same saved
outer patient-pair plan*:

1. clinical-only formulas;
2. clinical + PDR-index formulas.

It exercises unique inner leave-pair-out or repeated K-fold validation,
leakage-safe training-fold standardization, stability-aware model selection,
Verweij-van Houwelingen cross-validated partial likelihood, two Parvandeh-inspired
consensus modes, selectable cox-filter PLL/C-index group scoring, batched
unique-pair consensus resampling, outer multiprocessing, pickle persistence,
and PGFPlots/PNG/CSV graphical export.
"""

import argparse
import difflib
import os
import re
import warnings

import numpy as np
import pandas as pd

from cox_batch import (
    export_cindex_comparison_graphics,
    export_profile_cindex_summary_graphics,
    load_evaluation_result_pickle,
    load_outer_pair_plan,
    make_outer_pair_plan,
    nested_pair_cox_evaluation,
    save_evaluation_result_pickle,
    save_outer_pair_plan,
)


# -----------------------------------------------------------------------------
# Command-line selectable formula/group profiles
# CRC / irinotecan progressive-disease focused shortlist
# -----------------------------------------------------------------------------

META = "American_Joint_Committee_on_Cancer_Metastasis_Stage_Code_Reduced"
OVERALL_STAGE = "Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code_Encoded"
NODE_STAGE = "Neoplasm_Disease_Lymph_Node_Stage_AJCCC_Encoded"
TUMOR_STAGE = "American_Joint_Committee_on_Cancer_Tumor_Stage_Code_Encoded"
DRUG_PREDICTION = "P_DrugMerge_Prediction"
SCALAR_PREDICTION = "Scalar_prediction"
SCALAR_CRC = "Scalar_crc_prediction"
SCALAR_SCORE = "Scalar_score"

META_FORMULA = "C({})".format(META)
DRUG_FORMULA = "C({})".format(DRUG_PREDICTION)
NT_STAGE_FORMULA = "{} + {}".format(NODE_STAGE, TUMOR_STAGE)


def _profile(formulas, groups, expressions, description, consensus_note=""):
    return {
        "formulas": tuple(formulas),
        "groups": {str(k): tuple(v) for k, v in groups.items()},
        "expressions": {str(k): str(v) for k, v in expressions.items()},
        "description": str(description),
        "consensus_note": str(consensus_note),
    }


STAGE_ONLY_GROUPS = {
    "overall_ajcc_stage": (OVERALL_STAGE,),
}
STAGE_ONLY_EXPRESSIONS = {
    "overall_ajcc_stage": OVERALL_STAGE,
}
STAGE_ONLY_FORMULAS = (OVERALL_STAGE,)

CLINICAL_MODEL_PROFILES = {
    "stage-only": _profile(
        STAGE_ONLY_FORMULAS,
        STAGE_ONLY_GROUPS,
        STAGE_ONLY_EXPRESSIONS,
        "Fixed overall AJCC stage clinical benchmark for the CRC/irinotecan progressive-disease cohort.",
    ),
}


def _groups(*items):
    mapping = {
        "overall_ajcc_stage": (OVERALL_STAGE,),
        "Scalar_prediction": (SCALAR_PREDICTION,),
        "Scalar_crc_prediction": (SCALAR_CRC,),
        "Scalar_score": (SCALAR_SCORE,),
    }
    return {name: mapping[name] for name in items}


def _expressions(*items):
    mapping = {
        "overall_ajcc_stage": OVERALL_STAGE,
        "Scalar_prediction": SCALAR_PREDICTION,
        "Scalar_crc_prediction": SCALAR_CRC,
        "Scalar_score": SCALAR_SCORE,
    }
    return {name: mapping[name] for name in items}


CRC_ONLY_GROUPS = _groups("Scalar_crc_prediction")
CRC_ONLY_EXPRESSIONS = _expressions("Scalar_crc_prediction")
CRC_ONLY_FORMULAS = (SCALAR_CRC,)

STAGE_CRC_FIXED_GROUPS = _groups("overall_ajcc_stage", "Scalar_crc_prediction")
STAGE_CRC_FIXED_EXPRESSIONS = _expressions("overall_ajcc_stage", "Scalar_crc_prediction")
STAGE_CRC_FIXED_FORMULAS = (
    OVERALL_STAGE + " + " + SCALAR_CRC,
)

STAGE_VS_CRC_GROUPS = _groups("overall_ajcc_stage", "Scalar_crc_prediction")
STAGE_VS_CRC_EXPRESSIONS = _expressions("overall_ajcc_stage", "Scalar_crc_prediction")
STAGE_VS_CRC_FORMULAS = (
    OVERALL_STAGE,
    SCALAR_CRC,
    OVERALL_STAGE + " + " + SCALAR_CRC,
)

STAGE_VS_SCALAR_GROUPS = _groups("overall_ajcc_stage", "Scalar_prediction")
STAGE_VS_SCALAR_EXPRESSIONS = _expressions("overall_ajcc_stage", "Scalar_prediction")
STAGE_VS_SCALAR_FORMULAS = (
    OVERALL_STAGE,
    SCALAR_PREDICTION,
    OVERALL_STAGE + " + " + SCALAR_PREDICTION,
)

STAGE_PDR_COMPACT_GROUPS = _groups(
    "overall_ajcc_stage",
    "Scalar_crc_prediction",
    "Scalar_prediction",
    "Scalar_score",
)
STAGE_PDR_COMPACT_EXPRESSIONS = _expressions(
    "overall_ajcc_stage",
    "Scalar_crc_prediction",
    "Scalar_prediction",
    "Scalar_score",
)
STAGE_PDR_COMPACT_FORMULAS = (
    OVERALL_STAGE,
    OVERALL_STAGE + " + " + SCALAR_CRC,
    OVERALL_STAGE + " + " + SCALAR_PREDICTION,
    OVERALL_STAGE + " + " + SCALAR_SCORE,
)

CLINICAL_PDR_MODEL_PROFILES = {
    "crc-only": _profile(
        CRC_ONLY_FORMULAS,
        CRC_ONLY_GROUPS,
        CRC_ONLY_EXPRESSIONS,
        "Fixed Scalar_crc_prediction PDR benchmark for the CRC/irinotecan progressive-disease cohort.",
    ),
    "stage-crc-fixed": _profile(
        STAGE_CRC_FIXED_FORMULAS,
        STAGE_CRC_FIXED_GROUPS,
        STAGE_CRC_FIXED_EXPRESSIONS,
        "Fixed two-predictor benchmark: overall AJCC stage plus Scalar_crc_prediction.",
        "For a completely fixed formula, use --consensus-method none; consensus methods may regenerate formulas from the two groups.",
    ),
    "stage-vs-crc": _profile(
        STAGE_VS_CRC_FORMULAS,
        STAGE_VS_CRC_GROUPS,
        STAGE_VS_CRC_EXPRESSIONS,
        "Three-formula comparison of overall AJCC stage, Scalar_crc_prediction, and their additive combination.",
    ),
    "stage-vs-scalar": _profile(
        STAGE_VS_SCALAR_FORMULAS,
        STAGE_VS_SCALAR_GROUPS,
        STAGE_VS_SCALAR_EXPRESSIONS,
        "Three-formula comparison of overall AJCC stage, Scalar_prediction, and their additive combination.",
    ),
    "stage-pdr-compact": _profile(
        STAGE_PDR_COMPACT_FORMULAS,
        STAGE_PDR_COMPACT_GROUPS,
        STAGE_PDR_COMPACT_EXPRESSIONS,
        "Compact irinotecan progressive-disease PDR library: stage alone and stage augmented separately by Scalar_crc_prediction, Scalar_prediction, or Scalar_score.",
        "For small-cohort consensus runs, --consensus-max-formula-groups 2 is recommended so regenerated models remain low-dimensional.",
    ),
}

# Backward-compatible constant aliases for users importing the old names.
# They now resolve to the irinotecan progressive-disease focused defaults rather
# than the retired cross-cohort formula library.
CLINICAL_FORMULAS = STAGE_ONLY_FORMULAS
CLINICAL_PLUS_PDR_FORMULAS = STAGE_PDR_COMPACT_FORMULAS
CLINICAL_CONSENSUS_GROUPS = STAGE_ONLY_GROUPS
CLINICAL_GROUP_EXPRESSIONS = STAGE_ONLY_EXPRESSIONS
CLINICAL_PDR_CONSENSUS_GROUPS = STAGE_PDR_COMPACT_GROUPS
CLINICAL_PDR_GROUP_EXPRESSIONS = STAGE_PDR_COMPACT_EXPRESSIONS

def get_clinical_model_profile(name):
    return CLINICAL_MODEL_PROFILES[str(name)]


def get_clinical_pdr_model_profile(name):
    return CLINICAL_PDR_MODEL_PROFILES[str(name)]


def _profile_batch_mode(args):
    if getattr(args, "run_all_clinical_and_pdr_profiles", False):
        return "all"
    if getattr(args, "run_all_clinical_profiles", False):
        return "clinical_only"
    if getattr(args, "run_all_clinical_plus_pdr_profiles", False):
        return "clinical_plus_pdr"
    return None


def _batch_profile_specs(args):
    mode = _profile_batch_mode(args)
    specs = []
    if mode in ("clinical_only", "all"):
        for name in sorted(CLINICAL_MODEL_PROFILES):
            specs.append(("clinical_only", name))
    if mode in ("clinical_plus_pdr", "all"):
        for name in sorted(CLINICAL_PDR_MODEL_PROFILES):
            specs.append(("clinical_plus_pdr", name))
    return specs


def _copy_args(args):
    return argparse.Namespace(**dict(vars(args)))


def _args_for_profile(args, experiment, profile_name):
    profile_args = _copy_args(args)
    if experiment == "clinical_only":
        profile_args.clinical_only_model_profile = str(profile_name)
    elif experiment == "clinical_plus_pdr":
        profile_args.clinical_plus_pdr_model_profile = str(profile_name)
    else:
        raise ValueError("Unknown experiment: {}".format(experiment))
    return profile_args



_FORMULA_TOKEN_PATTERN = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_FORMULA_NON_FIELDS = set((
    "C", "I", "Q", "bs", "cr", "strata", "np", "log", "exp", "sqrt",
    "True", "False",
))


def _formula_field_names(formula):
    """Extract identifier-like field names from the simple lifelines formulas used here."""
    return sorted(set(
        token for token in _FORMULA_TOKEN_PATTERN.findall(str(formula))
        if token not in _FORMULA_NON_FIELDS
    ))


def validate_formula_fields(data, formulas, label):
    """Warn about formula fields absent from the supplied DataFrame.

    The formulas in this driver use ordinary identifier-style pandas column
    names and ``C(name)`` categorical wrappers, so identifier extraction is
    sufficient and remains independent of patsy/formulaic internals.
    """
    columns = set(str(value) for value in data.columns)
    missing_any = False
    for formula_index, formula in enumerate(formulas, start=1):
        missing = [name for name in _formula_field_names(formula) if name not in columns]
        if not missing:
            continue
        missing_any = True
        details = []
        for name in missing:
            close = difflib.get_close_matches(name, sorted(columns), n=3, cutoff=0.55)
            if close:
                details.append("{} (possible match: {})".format(name, ", ".join(close)))
            else:
                details.append(name)
        warnings.warn(
            "{} formula {} references missing DataFrame field(s): {}. Formula: {}".format(
                label, formula_index, "; ".join(details), formula
            ),
            RuntimeWarning,
            stacklevel=2,
        )
    if not missing_any:
        print("{} field-name check: all formula fields are present.".format(label))
    return not missing_any


def validate_group_fields(data, groups, label):
    columns = set(str(value) for value in data.columns)
    missing_any = False
    for group_name in sorted(groups):
        missing = [str(name) for name in groups[group_name] if str(name) not in columns]
        if not missing:
            continue
        missing_any = True
        details = []
        for name in missing:
            close = difflib.get_close_matches(name, sorted(columns), n=3, cutoff=0.55)
            if close:
                details.append("{} (possible match: {})".format(name, ", ".join(close)))
            else:
                details.append(name)
        warnings.warn(
            "{} group {!r} references missing DataFrame field(s): {}".format(
                label, group_name, "; ".join(details)
            ),
            RuntimeWarning,
            stacklevel=2,
        )
    if not missing_any:
        print("{} field-name check: all group fields are present.".format(label))
    return not missing_any


def _validate_profile_fields(data, profile, label):
    validate_formula_fields(data, profile["formulas"], label + " seed formulas")
    validate_group_fields(data, profile["groups"], label + " consensus groups")
    validate_formula_fields(
        data, tuple(profile["expressions"].values()), label + " group expressions"
    )


def validate_input_frame(data, clinical_profile=None, pdr_profile=None):
    if not isinstance(data, pd.DataFrame):
        raise TypeError("Input pickle must contain a pandas DataFrame.")
    required = ("Response_duration", "Response_coded")
    missing = [name for name in required if name not in data.columns]
    if missing:
        raise ValueError(
            "Input DataFrame is missing required survival field(s): {}".format(
                ", ".join(missing)
            )
        )
    if clinical_profile is None:
        clinical_profile = get_clinical_model_profile("stage-only")
    if pdr_profile is None:
        pdr_profile = get_clinical_pdr_model_profile("stage-pdr-compact")
    _validate_profile_fields(data, clinical_profile, "Clinical-only profile")
    _validate_profile_fields(data, pdr_profile, "Clinical+PDR profile")
    return data


def load_main_input_data(args):
    if args.data_pickle:
        path = os.path.abspath(args.data_pickle)
        warnings.warn(
            "Loading a pandas pickle executes Python pickle deserialization. "
            "Only load files from trusted sources: {}".format(path),
            RuntimeWarning,
            stacklevel=2,
        )
        data = pd.read_pickle(path)
        print("Loaded input pandas DataFrame from: {}".format(path))
    else:
        data = generate_brca_like_22_patient_dataset(args.dataset_seed)
        print("Using generated 22-patient TCGA-BRCA-like synthetic DataFrame.")
    batch_specs = _batch_profile_specs(args)
    if batch_specs:
        validate_input_frame(data)
        for experiment, profile_name in batch_specs:
            if experiment == "clinical_only":
                profile = get_clinical_model_profile(profile_name)
                _validate_profile_fields(
                    data, profile, "Clinical-only profile {}".format(profile_name)
                )
            else:
                profile = get_clinical_pdr_model_profile(profile_name)
                _validate_profile_fields(
                    data, profile, "Clinical+PDR profile {}".format(profile_name)
                )
        return data
    clinical_profile = get_clinical_model_profile(args.clinical_only_model_profile)
    pdr_profile = get_clinical_pdr_model_profile(args.clinical_plus_pdr_model_profile)
    return validate_input_frame(data, clinical_profile, pdr_profile)

def generate_brca_like_22_patient_dataset(random_seed=20260902):
    """Generate a small survival cohort with TCGA-BRCA-like variable names.

    The synthetic cohort deliberately resembles the 22-patient real-data subset
    discussed in the accompanying analysis: all patients experience the event,
    observed durations are integer-like and include ties, clinical staging is
    low-dimensional, and three continuous PDR-like indices contain a moderate
    but noisy prognostic signal.
    """

    rng = np.random.RandomState(int(random_seed))
    n = 22

    tumor_stage = rng.choice([1, 2, 3, 4], size=n, p=[0.10, 0.40, 0.35, 0.15])
    node_stage = rng.choice([0, 1, 2, 3], size=n, p=[0.25, 0.35, 0.25, 0.15])
    metastasis_binary = rng.binomial(1, 0.14, size=n)

    overall_stage = np.maximum(tumor_stage, np.minimum(4, node_stage + 1))
    overall_stage = np.maximum(overall_stage, 2 + 2 * metastasis_binary)
    overall_stage = np.minimum(overall_stage, 4).astype(int)

    pdr_latent = rng.normal(size=n)
    scalar_prediction = 0.85 * pdr_latent + 0.35 * rng.normal(size=n)
    scalar_crc_prediction = 0.55 * pdr_latent + 0.65 * rng.normal(size=n)
    scalar_score = 0.40 * pdr_latent + 0.80 * rng.normal(size=n)

    drug_class = np.where(
        pdr_latent > 0.65,
        "PDR_high",
        np.where(pdr_latent < -0.65, "PDR_low", "PDR_mid"),
    )
    metastasis_label = np.where(metastasis_binary == 1, "M1", "M0")

    # Moderate prognostic signal with substantial noise, matching the difficult
    # small-sample setting rather than creating an unrealistically easy cohort.
    linear_predictor = (
        0.20 * (tumor_stage - 2.0)
        + 0.22 * node_stage
        + 0.45 * metastasis_binary
        + 0.45 * scalar_prediction
        + 0.18 * scalar_crc_prediction
        + 0.10 * scalar_score
    )
    # Integer-like observed durations in a range similar to the small BRCA
    # subset used in the real-data experiment. Higher latent risk tends to
    # shorten the duration, but substantial noise keeps discrimination modest.
    event_time = 82.0 - 11.0 * linear_predictor + rng.normal(0.0, 16.0, size=n)
    event_time = np.clip(event_time, 38.0, 125.0)
    observed_time = np.round(event_time).astype(float)

    frame = pd.DataFrame(
        {
            "Response_duration": observed_time,
            "Response_coded": np.ones(n, dtype=int),
            "American_Joint_Committee_on_Cancer_Metastasis_Stage_Code_Reduced": metastasis_label,
            "P_DrugMerge_Prediction": drug_class,
            "Neoplasm_Disease_Lymph_Node_Stage_AJCCC_Encoded": node_stage.astype(float),
            "Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code_Encoded": overall_stage.astype(float),
            "American_Joint_Committee_on_Cancer_Tumor_Stage_Code_Encoded": tumor_stage.astype(float),
            "Scalar_prediction": scalar_prediction,
            "Scalar_crc_prediction": scalar_crc_prediction,
            "Scalar_score": scalar_score,
        },
        index=["TCGA-SYN-{:04d}".format(i + 1) for i in range(n)],
    )
    return frame


def _resolved_consensus_top_features(args, experiment):
    if experiment == "clinical_only":
        value = args.clinical_only_consensus_top_features
    elif experiment == "clinical_plus_pdr":
        value = args.clinical_plus_pdr_consensus_top_features
    else:
        raise ValueError("Unknown experiment: {}".format(experiment))
    if value is None:
        value = args.consensus_top_features
    value = int(value)
    if value < 1:
        raise ValueError("Consensus top-features values must be >= 1.")
    return value


def print_dataset_summary(data, args):
    clinical_profile = get_clinical_model_profile(args.clinical_only_model_profile)
    pdr_profile = get_clinical_pdr_model_profile(args.clinical_plus_pdr_model_profile)
    print("\n{}".format("=" * 78))
    print("MAIN INPUT SURVIVAL COHORT")
    print("{}".format("=" * 78))
    print("Patients : {}".format(len(data)))
    print("Events   : {}".format(int(data["Response_coded"].sum())))
    print("Censored : {}".format(int((1 - data["Response_coded"]).sum())))
    print("Durations: min={:.1f}, median={:.1f}, max={:.1f}".format(
        float(data["Response_duration"].min()),
        float(data["Response_duration"].median()),
        float(data["Response_duration"].max()),
    ))
    print("\nClinical model profile: {}".format(args.clinical_only_model_profile))
    print("  {}".format(clinical_profile["description"]))
    print("  seed formulas: {}; groups: {}; consensus_top_features: {}".format(
        len(clinical_profile["formulas"]), len(clinical_profile["groups"]),
        _resolved_consensus_top_features(args, "clinical_only")
    ))
    print("Clinical + PDR model profile: {}".format(args.clinical_plus_pdr_model_profile))
    print("  {}".format(pdr_profile["description"]))
    print("  seed formulas: {}; groups: {}; consensus_top_features: {}".format(
        len(pdr_profile["formulas"]), len(pdr_profile["groups"]),
        _resolved_consensus_top_features(args, "clinical_plus_pdr")
    ))
    if pdr_profile.get("consensus_note"):
        warnings.warn(pdr_profile["consensus_note"], RuntimeWarning, stacklevel=2)
    print(data.head())


def _common_evaluation_kwargs(args, outer_pair_plan, consensus_top_features):
    return dict(
        duration_col="Response_duration",
        event_col="Response_coded",
#        base_models=("breslow", "spline", "piecewise"),
        base_models=("breslow",),
        knot_numbers=(2, 3),
        piecewise_break_numbers=(1, 2),
#        penalizers=(0.03, 0.06, 0.1, 0.15, 0.2),
        penalizers=(0.03, ),
        l1_ratios=(0.0,),
        inner_splits=args.inner_splits,
        inner_repeats=args.inner_repeats,
        inner_validation_method=args.inner_method,
        inner_max_pairs=args.inner_max_pairs,
        inner_scoring=args.inner_scoring,
        inner_selection_rule=args.inner_selection_rule,
        inner_stability_penalty=args.inner_stability_penalty,
        consensus_method=args.consensus_method,
        consensus_top_models=args.consensus_top_models,
        consensus_top_features=int(consensus_top_features),
        consensus_threshold=args.consensus_threshold,
        consensus_max_formula_groups=args.consensus_max_formula_groups,
        consensus_pair_batch_size=args.consensus_pair_batch_size,
        consensus_filter_scoring=args.consensus_filter_scoring,
        prediction_aggregation=args.prediction_aggregation,
        ensemble_size=args.ensemble_size,
        ensemble_weighting=args.ensemble_weighting,
        ensemble_temperature=args.ensemble_temperature,
        ensemble_distinct_formulas=args.ensemble_distinct_formulas,
        ensemble_rank_level=args.ensemble_rank_level,
        standardize_predictors=args.standardize_predictors,
        standardize_columns=None,
        max_outer_pairs=args.max_outer_pairs,
        event_policy=args.event_policy,
        outer_pairs=outer_pair_plan,
        n_jobs=args.n_jobs,
        multiprocessing_start_method="fork",
        n_bootstrap=args.bootstrap,
        confidence_level=0.95,
        random_state=args.analysis_seed,
        fit_final_model=(not args.no_final_model),
        verbose=True,
    )


def _tag_result_profile(result, profile_name, experiment_name):
    result.model_profile = str(profile_name)
    result.experiment_name = str(experiment_name)
    result.outer_results["model_profile"] = str(profile_name)
    result.outer_results["experiment_name"] = str(experiment_name)
    return result


def run_clinical_only(data, args, outer_pair_plan):
    profile = get_clinical_model_profile(args.clinical_only_model_profile)
    kwargs = _common_evaluation_kwargs(
        args, outer_pair_plan,
        _resolved_consensus_top_features(args, "clinical_only")
    )
    kwargs["formulas"] = profile["formulas"]
    kwargs["consensus_feature_groups"] = profile["groups"]
    kwargs["consensus_group_expressions"] = profile["expressions"]
    result = nested_pair_cox_evaluation(data, **kwargs)
    return _tag_result_profile(result, args.clinical_only_model_profile, "clinical_only")


def run_clinical_plus_pdr(data, args, outer_pair_plan):
    profile = get_clinical_pdr_model_profile(args.clinical_plus_pdr_model_profile)
    kwargs = _common_evaluation_kwargs(
        args, outer_pair_plan,
        _resolved_consensus_top_features(args, "clinical_plus_pdr")
    )
    kwargs["formulas"] = profile["formulas"]
    kwargs["consensus_feature_groups"] = profile["groups"]
    kwargs["consensus_group_expressions"] = profile["expressions"]
    result = nested_pair_cox_evaluation(data, **kwargs)
    return _tag_result_profile(
        result, args.clinical_plus_pdr_model_profile, "clinical_plus_pdr"
    )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Run two 22-patient BRCA-like nested Cox analyses using one reusable "
            "outer-pair plan and export pickle/PGFPlots/PNG/CSV results."
        )
    )
    parser.add_argument(
        "--mode",
        choices=("run", "graphics-only"),
        default="run",
        help="Fresh analysis or reload existing pickle files (default: run).",
    )
    parser.add_argument("--n-jobs", type=int, default=8, help="Outer-pair workers (1-32).")
    parser.add_argument("--bootstrap", type=int, default=500, help="Bootstrap replicates.")
    parser.add_argument(
        "--event-policy",
        choices=("all", "both_events", "valid_pairs", "harrell_comparable"),
        default="harrell_comparable",
        help="Outer pair policy (default: harrell_comparable).",
    )
    parser.add_argument(
        "--inner-method",
        choices=("kfold", "unique_pairs"),
        default="unique_pairs",
        help="Inner validation design (default: unique_pairs).",
    )
    parser.add_argument(
        "--inner-max-pairs",
        type=int,
        default=60,
        help="Unique Harrell-comparable inner pairs when --inner-method=unique_pairs.",
    )
    parser.add_argument(
        "--inner-splits",
        type=int,
        default=5,
        help="Folds for repeated K-fold mode (default: 5).",
    )
    parser.add_argument(
        "--inner-repeats",
        type=int,
        default=1,
        help="Repeated randomizations for K-fold mode; ignored by unique_pairs.",
    )
    parser.add_argument(
        "--inner-scoring",
        choices=("c_index", "cv_partial_likelihood"),
        default="cv_partial_likelihood",
        help=(
            "Inner model-selection metric. cv_partial_likelihood uses the "
            "Verweij-van Houwelingen grouped partial log-likelihood; "
            "default: cv_partial_likelihood."
        ),
    )
    parser.add_argument(
        "--inner-selection-rule",
        choices=("mean", "lower_confidence"),
        default="lower_confidence",
        help=(
            "Rank by the raw mean selected metric or mean minus lambda*SE. "
            "The rule applies to either C-index or CV partial likelihood."
        ),
    )
    parser.add_argument(
        "--inner-stability-penalty",
        type=float,
        default=1.0,
        help="Lambda in the lower-confidence stability score (default: 1.0).",
    )
    parser.add_argument(
        "--consensus-method",
        choices=("none", "top_models", "cox_filter"),
        default="top_models",
        help=(
            "Parvandeh-inspired consensus mode: none keeps ordinary nested "
            "selection; top_models builds grouped consensus from the best "
            "models in each inner resample; cox_filter ranks groups by "
            "training-only Cox partial-likelihood improvement. Default: top_models."
        ),
    )
    parser.add_argument(
        "--consensus-top-models",
        type=int,
        default=5,
        help="Top models retained per inner resample in top_models mode (default: 5).",
    )
    parser.add_argument(
        "--consensus-top-features",
        type=int,
        default=3,
        help="Top feature/groups retained per training resample in cox_filter mode (default: 3).",
    )
    parser.add_argument(
        "--clinical-only-consensus-top-features",
        type=int,
        default=None,
        help=(
            "Clinical-only override for --consensus-top-features. If omitted, "
            "the shared --consensus-top-features value is used."
        ),
    )
    parser.add_argument(
        "--clinical-plus-pdr-consensus-top-features",
        type=int,
        default=None,
        help=(
            "Clinical+PDR override for --consensus-top-features. If omitted, "
            "the shared --consensus-top-features value is used."
        ),
    )
    parser.add_argument(
        "--consensus-threshold",
        type=float,
        default=0.60,
        help="Minimum grouped selection frequency for consensus inclusion (default: 0.60).",
    )
    parser.add_argument(
        "--consensus-max-formula-groups",
        type=int,
        default=3,
        help="Maximum consensus groups in an internally generated formula (default: 3).",
    )
    parser.add_argument(
        "--consensus-pair-batch-size",
        type=int,
        default=10,
        help=(
            "Unique inner pairs pooled into each consensus resample when "
            "--inner-method=unique_pairs (default: 10). Use 1 for the previous "
            "single-pair consensus behavior; ignored by K-fold consensus."
        ),
    )
    parser.add_argument(
        "--consensus-filter-scoring",
        choices=("partial_likelihood", "c_index"),
        default="partial_likelihood",
        help=(
            "Group-scoring rule for --consensus-method=cox_filter. "
            "partial_likelihood preserves training-only Cox PLL improvement; "
            "c_index fits each group on inner training patients and scores "
            "held-out Harrell concordance within each fold/pair batch. "
            "Ignored by other consensus methods. Default: partial_likelihood."
        ),
    )
    parser.add_argument(
        "--prediction-aggregation",
        choices=("winner", "stacked"),
        default="stacked",
        help=(
            "Outer prediction rule after inner model ranking. winner keeps the "
            "single best model; stacked averages standardized log-risk scores "
            "from several top eligible models (default: stacked)."
        ),
    )
    parser.add_argument(
        "--ensemble-size",
        type=int,
        default=5,
        help="Maximum fitted components in stacked prediction mode (default: 5).",
    )
    parser.add_argument(
        "--ensemble-weighting",
        choices=("equal", "rank", "softmax"),
        default="equal",
        help=(
            "Stacked prediction weights: equal, inverse inner-rank, or softmax "
            "of standardized inner selection scores (default: equal)."
        ),
    )
    parser.add_argument(
        "--ensemble-rank-level",
        choices=("specification", "formula"),
        default="specification",
        help=(
            "Ranks used by --ensemble-weighting=rank. specification uses raw "
            "hyperparameter-specification ranks; formula uses dense ranks of "
            "distinct formulas (default: specification for backward compatibility)."
        ),
    )
    parser.add_argument(
        "--ensemble-temperature",
        type=float,
        default=1.0,
        help="Softmax temperature; ignored by equal/rank weighting (default: 1.0).",
    )
    parser.add_argument(
        "--ensemble-allow-duplicate-formulas",
        action="store_false",
        dest="ensemble_distinct_formulas",
        default=True,
        help=(
            "Allow several hyperparameter/baseline variants of the same formula "
            "to enter one ensemble. By default only the best fitted specification "
            "per distinct formula is used."
        ),
    )
    parser.add_argument(
        "--no-standardize-predictors",
        action="store_false",
        dest="standardize_predictors",
        default=True,
        help="Disable training-fold z-score standardization of numeric predictors.",
    )
    parser.add_argument(
        "--max-outer-pairs",
        type=int,
        default=120,
        help="Maximum outer pairs when generating the reusable plan (default: 120).",
    )
    parser.add_argument(
        "--data-pickle",
        default=None,
        help=(
            "Trusted pickle containing the pandas DataFrame to analyze. If omitted, "
            "the 22-patient TCGA-BRCA-like synthetic generator is used."
        ),
    )
    parser.add_argument(
        "--clinical-only-model-profile",
        choices=tuple(sorted(CLINICAL_MODEL_PROFILES)),
        default="stage-only",
        help=(
            "Seed-formula + consensus-group profile for the clinical-only experiment. "
            "Use --list-model-profiles to print detailed definitions. "
            "Default: stage-only."
        ),
    )
    parser.add_argument(
        "--clinical-plus-pdr-model-profile",
        choices=tuple(sorted(CLINICAL_PDR_MODEL_PROFILES)),
        default="stage-pdr-compact",
        help=(
            "Seed-formula + consensus-group profile for the clinical+PDR experiment. "
            "Use --list-model-profiles to print detailed definitions. "
            "Default: stage-pdr-compact."
        ),
    )
    batch_group = parser.add_mutually_exclusive_group()
    batch_group.add_argument(
        "--run-all-clinical-profiles",
        action="store_true",
        help=(
            "Run every available clinical-only model profile using the same "
            "outer-pair plan. Each profile is written to its own result/pickle "
            "subfolder, with a common C-index summary figure at the root."
        ),
    )
    batch_group.add_argument(
        "--run-all-clinical-plus-pdr-profiles",
        action="store_true",
        help=(
            "Run every available clinical+PDR model profile using the same "
            "outer-pair plan. Each profile is written to its own result/pickle "
            "subfolder, with a common C-index summary figure at the root."
        ),
    )
    batch_group.add_argument(
        "--run-all-clinical-and-pdr-profiles",
        action="store_true",
        help=(
            "Run every available clinical-only and clinical+PDR model profile "
            "using one common outer-pair plan. Each profile is written to its "
            "own result/pickle subfolder, with one combined C-index summary "
            "figure plus family-specific summaries."
        ),
    )
    parser.add_argument(
        "--list-model-profiles",
        action="store_true",
        help="Print all selectable clinical and PDR model profiles and exit.",
    )
    parser.add_argument(
        "--dataset-seed", type=int, default=20260902, help="Synthetic-data seed."
    )
    parser.add_argument(
        "--no-final-model",
        action="store_true",
        help=(
            "Skip the expensive final full-data inner ranking/final-model fit. "
            "Outer C-index estimation, bootstrap, pickle, and graphics are retained."
        ),
    )
    parser.add_argument(
        "--pair-seed", type=int, default=1001, help="Reusable outer-pair-plan seed."
    )
    parser.add_argument(
        "--analysis-seed", type=int, default=2001, help="Inner/bootstrap analysis seed."
    )
    parser.add_argument(
        "--output-dir", default="brca_like_fast_figures", help="Graphics/CSV directory."
    )
    parser.add_argument(
        "--pickle-dir", default="brca_like_fast_pickles", help="Pickle directory."
    )
    parser.add_argument(
        "--outer-pair-file",
        default="brca_like_fast_fixed_outer_pairs.csv",
        help="Reusable outer-pair CSV. Existing file is validated and reused.",
    )
    parser.add_argument("--png-dpi", type=int, default=300)
    parser.add_argument("--histogram-bins", type=int, default=20)
    parser.add_argument("--top-n-models", type=int, default=10)
    parser.add_argument(
        "--include-final-model-in-pickle",
        action="store_true",
        help="Also serialize the fitted final lifelines Cox model.",
    )
    return parser.parse_args()


def _pickle_paths(args):
    stem = "{}_{}_{}_{}_{}_{}_batch{}_{}_{}".format(
        args.event_policy,
        args.inner_method,
        args.inner_scoring,
        args.inner_selection_rule,
        args.consensus_method,
        args.consensus_filter_scoring,
        args.consensus_pair_batch_size,
        args.prediction_aggregation,
        "{}_{}".format(args.ensemble_weighting, args.ensemble_rank_level),
    )
    profile_suffix = []
    if args.clinical_only_model_profile != "stage-only":
        profile_suffix.append("clin-{}".format(args.clinical_only_model_profile))
    if args.clinical_plus_pdr_model_profile != "stage-pdr-compact":
        profile_suffix.append("pdr-{}".format(args.clinical_plus_pdr_model_profile))
    clinical_tf = _resolved_consensus_top_features(args, "clinical_only")
    pdr_tf = _resolved_consensus_top_features(args, "clinical_plus_pdr")
    if (args.clinical_only_consensus_top_features is not None or
            args.clinical_plus_pdr_consensus_top_features is not None):
        profile_suffix.append("ctf{}-ptf{}".format(clinical_tf, pdr_tf))
    if profile_suffix:
        stem = stem + "__" + "__".join(profile_suffix)
    return (
        os.path.join(args.pickle_dir, "clinical_only_{}.pkl".format(stem)),
        os.path.join(args.pickle_dir, "clinical_plus_pdr_{}.pkl".format(stem)),
    )


def obtain_fixed_outer_pair_plan(data, args):
    path = os.path.abspath(args.outer_pair_file)
    if os.path.isfile(path):
        print("Reusing fixed outer-pair plan: {}".format(path))
        plan = load_outer_pair_plan(path)
    else:
        plan = make_outer_pair_plan(
            data,
            duration_col="Response_duration",
            event_col="Response_coded",
            max_outer_pairs=args.max_outer_pairs,
            event_policy=args.event_policy,
            random_state=args.pair_seed,
        )
        save_outer_pair_plan(plan, path)
        print("Created fixed outer-pair plan: {}".format(path))
    print("Outer pairs in plan: {}".format(len(plan)))
    return plan


def run_and_pickle_results(args):
    data = load_main_input_data(args)
    print_dataset_summary(data, args)
    pair_plan = obtain_fixed_outer_pair_plan(data, args)

    print("\nRunning CLINICAL-ONLY analysis")
    clinical = run_clinical_only(data, args, pair_plan)
    print("\nRunning CLINICAL + PDR analysis")
    clinical_pdr = run_clinical_plus_pdr(data, args, pair_plan)

    # Verify the paired comparison really used the exact same outer pairs.
    clinical_pairs = list(zip(
        clinical.outer_results["patient_position_i"],
        clinical.outer_results["patient_position_j"],
    ))
    pdr_pairs = list(zip(
        clinical_pdr.outer_results["patient_position_i"],
        clinical_pdr.outer_results["patient_position_j"],
    ))
    if clinical_pairs != pdr_pairs:
        raise RuntimeError("The two analyses did not use identical outer pairs.")

    clinical_pickle, pdr_pickle = _pickle_paths(args)
    include_model = bool(args.include_final_model_in_pickle)
    save_evaluation_result_pickle(clinical, clinical_pickle, include_model, protocol=4)
    save_evaluation_result_pickle(clinical_pdr, pdr_pickle, include_model, protocol=4)
    print("\nSaved pickle files:")
    print("  {}".format(os.path.abspath(clinical_pickle)))
    print("  {}".format(os.path.abspath(pdr_pickle)))
    return clinical, clinical_pdr


def load_pickled_results(args):
    clinical_pickle, pdr_pickle = _pickle_paths(args)
    print("Loading:\n  {}\n  {}".format(clinical_pickle, pdr_pickle))
    return (
        load_evaluation_result_pickle(clinical_pickle),
        load_evaluation_result_pickle(pdr_pickle),
    )


def write_latex_include_example(output_dir):
    path = os.path.join(output_dir, "latex_include_example.tex")
    content = r"""% Preamble:
% \usepackage{pgfplots}
% \pgfplotsset{compat=1.15}
% \usepackage{graphicx}

\begin{figure}[htbp]
\centering
\input{clinical_only/clinical_only_bootstrap_histogram.tex}
\caption{Bootstrap distribution of outer concordance for the clinical-only Cox analysis.}
\end{figure}

\begin{figure}[htbp]
\centering
\input{clinical_plus_pdr/clinical_plus_pdr_bootstrap_histogram.tex}
\caption{Bootstrap distribution of outer concordance for the clinical plus PDR Cox analysis.}
\end{figure}

\begin{figure}[htbp]
\centering
\input{clinical_vs_pdr_cindex.tex}
\caption{Outer concordance estimates and confidence intervals using the same fixed outer pairs.}
\end{figure}
"""
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(content)
    return path


def export_all_graphics(clinical, clinical_pdr, args):
    root = os.path.abspath(args.output_dir)
    clinical_dir = os.path.join(root, "clinical_only")
    pdr_dir = os.path.join(root, "clinical_plus_pdr")
    clinical_paths = clinical.export_graphics(
        clinical_dir,
        prefix="clinical_only",
        top_n_models=args.top_n_models,
        histogram_bins=args.histogram_bins,
        png_dpi=args.png_dpi,
    )
    pdr_paths = clinical_pdr.export_graphics(
        pdr_dir,
        prefix="clinical_plus_pdr",
        top_n_models=args.top_n_models,
        histogram_bins=args.histogram_bins,
        png_dpi=args.png_dpi,
    )
    comparison = export_cindex_comparison_graphics(
        [clinical, clinical_pdr],
        ["Clinical only", "Clinical + PDR"],
        os.path.join(root, "clinical_vs_pdr_cindex"),
        png_dpi=args.png_dpi,
    )
    if not os.path.isdir(clinical_dir):
        os.makedirs(clinical_dir)
    if not os.path.isdir(pdr_dir):
        os.makedirs(pdr_dir)
    clinical.outer_results.to_csv(
        os.path.join(clinical_dir, "clinical_only_outer_results.csv"), index=False
    )
    clinical_pdr.outer_results.to_csv(
        os.path.join(pdr_dir, "clinical_plus_pdr_outer_results.csv"), index=False
    )
    clinical.outer_pair_plan.to_csv(
        os.path.join(root, "outer_pair_plan_used.csv"), index=False
    )
    include_file = write_latex_include_example(root)

    print("\nGenerated graphical outputs:")
    for key in sorted(clinical_paths):
        print("  clinical {:24s} {}".format(key + ":", clinical_paths[key]))
    for key in sorted(pdr_paths):
        print("  clinical+PDR {:20s} {}".format(key + ":", pdr_paths[key]))
    print("  comparison tex: {}".format(comparison["tex"]))
    print("  comparison png: {}".format(comparison["png"]))
    print("  LaTeX example : {}".format(include_file))


def _batch_profile_output_dir(args, experiment, profile_name):
    return os.path.join(
        os.path.abspath(args.output_dir), str(experiment), str(profile_name)
    )


def _batch_profile_pickle_dir(args, experiment, profile_name):
    return os.path.join(
        os.path.abspath(args.pickle_dir), str(experiment), str(profile_name)
    )


def _batch_profile_pickle_path(args, experiment, profile_name):
    profile_args = _copy_args(args)
    # Batch pickle names depend only on the profile being evaluated, not on an
    # irrelevant single-profile selector that the user may also have supplied.
    profile_args.clinical_only_model_profile = "stage-only"
    profile_args.clinical_plus_pdr_model_profile = "stage-pdr-compact"
    profile_args = _args_for_profile(profile_args, experiment, profile_name)
    profile_args.pickle_dir = _batch_profile_pickle_dir(
        args, experiment, profile_name
    )
    clinical_path, pdr_path = _pickle_paths(profile_args)
    if experiment == "clinical_only":
        return clinical_path
    if experiment == "clinical_plus_pdr":
        return pdr_path
    raise ValueError("Unknown experiment: {}".format(experiment))


def _print_batch_dataset_summary(data, args, specs):
    print("\n{}".format("=" * 78))
    print("MAIN INPUT SURVIVAL COHORT -- PROFILE BATCH")
    print("{}".format("=" * 78))
    print("Patients : {}".format(len(data)))
    print("Events   : {}".format(int(data["Response_coded"].sum())))
    print("Censored : {}".format(int((1 - data["Response_coded"]).sum())))
    print("Durations: min={:.1f}, median={:.1f}, max={:.1f}".format(
        float(data["Response_duration"].min()),
        float(data["Response_duration"].median()),
        float(data["Response_duration"].max()),
    ))
    print("Batch mode: {}".format(_profile_batch_mode(args)))
    print("Profiles to evaluate: {}".format(len(specs)))
    for experiment, profile_name in specs:
        if experiment == "clinical_only":
            profile = get_clinical_model_profile(profile_name)
            top_features = _resolved_consensus_top_features(args, "clinical_only")
            label = "clinical-only"
        else:
            profile = get_clinical_pdr_model_profile(profile_name)
            top_features = _resolved_consensus_top_features(args, "clinical_plus_pdr")
            label = "clinical+PDR"
        print("  {:13s} {:30s} formulas={} groups={} top_features={}".format(
            label,
            profile_name,
            len(profile["formulas"]),
            len(profile["groups"]),
            top_features,
        ))
        if profile.get("consensus_note"):
            warnings.warn(
                "{} / {}: {}".format(label, profile_name, profile["consensus_note"]),
                RuntimeWarning,
                stacklevel=2,
            )


def _run_one_batch_profile(data, args, outer_pair_plan, experiment, profile_name):
    profile_args = _args_for_profile(args, experiment, profile_name)
    print("\n{}".format("-" * 78))
    print("Running {} profile: {}".format(experiment, profile_name))
    print("{}".format("-" * 78))
    if experiment == "clinical_only":
        return run_clinical_only(data, profile_args, outer_pair_plan)
    if experiment == "clinical_plus_pdr":
        return run_clinical_plus_pdr(data, profile_args, outer_pair_plan)
    raise ValueError("Unknown experiment: {}".format(experiment))


def _save_one_batch_profile(result, args, experiment, profile_name):
    pickle_path = _batch_profile_pickle_path(
        args, experiment, profile_name
    )
    pickle_parent = os.path.dirname(pickle_path)
    if pickle_parent and not os.path.isdir(pickle_parent):
        os.makedirs(pickle_parent)
    save_evaluation_result_pickle(
        result,
        pickle_path,
        bool(args.include_final_model_in_pickle),
        protocol=4,
    )
    return pickle_path


def _load_one_batch_profile(args, experiment, profile_name):
    pickle_path = _batch_profile_pickle_path(
        args, experiment, profile_name
    )
    print("Loading {} / {}: {}".format(
        experiment, profile_name, pickle_path
    ))
    return load_evaluation_result_pickle(pickle_path)


def _export_one_batch_profile(result, args, experiment, profile_name):
    output_dir = _batch_profile_output_dir(
        args, experiment, profile_name
    )
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    prefix = "{}_{}".format(experiment, profile_name.replace("-", "_"))
    paths = result.export_graphics(
        output_dir,
        prefix=prefix,
        top_n_models=args.top_n_models,
        histogram_bins=args.histogram_bins,
        png_dpi=args.png_dpi,
    )
    result.outer_results.to_csv(
        os.path.join(output_dir, "outer_results.csv"), index=False
    )
    result.outer_pair_plan.to_csv(
        os.path.join(output_dir, "outer_pair_plan_used.csv"), index=False
    )
    return paths


def _summary_for_profile_results(results_by_key, args, experiment=None):
    selected = []
    if experiment is None:
        for key in sorted(results_by_key):
            selected.append((key, results_by_key[key]))
        stem_name = "all_profiles_cindex"
    else:
        for key in sorted(results_by_key):
            if key[0] == experiment:
                selected.append((key, results_by_key[key]))
        if experiment == "clinical_only":
            stem_name = "clinical_profiles_cindex"
        else:
            stem_name = "clinical_plus_pdr_profiles_cindex"
    if not selected:
        return None

    results = []
    labels = []
    experiment_labels = []
    profile_names = []
    for (experiment_name, profile_name), result in selected:
        results.append(result)
        experiment_labels.append(experiment_name)
        profile_names.append(profile_name)
        if experiment is None:
            family = "Clinical" if experiment_name == "clinical_only" else "Clinical+PDR"
            labels.append("{} | {}".format(family, profile_name))
        else:
            labels.append(profile_name)

    root = os.path.abspath(args.output_dir)
    if not os.path.isdir(root):
        os.makedirs(root)
    paths = export_profile_cindex_summary_graphics(
        results=results,
        labels=labels,
        output_stem=os.path.join(root, stem_name),
        experiment_labels=experiment_labels,
        profile_names=profile_names,
        png_dpi=args.png_dpi,
    )
    return paths


def _verify_batch_common_outer_pairs(results_by_key):
    reference_key = None
    reference_pairs = None
    for key in sorted(results_by_key):
        result = results_by_key[key]
        columns = result.outer_results.columns
        if "patient_position_i" in columns and "patient_position_j" in columns:
            pairs = list(zip(
                result.outer_results["patient_position_i"],
                result.outer_results["patient_position_j"],
            ))
        elif "patient_id_i" in columns and "patient_id_j" in columns:
            pairs = list(zip(
                result.outer_results["patient_id_i"],
                result.outer_results["patient_id_j"],
            ))
        else:
            raise RuntimeError(
                "Cannot verify common outer pairs for profile {}: pair identifier "
                "columns are missing.".format(key)
            )
        if reference_pairs is None:
            reference_key = key
            reference_pairs = pairs
        elif pairs != reference_pairs:
            raise RuntimeError(
                "Profile {} does not use the same ordered outer-pair sequence as {}. "
                "Refusing to create a common profile comparison figure.".format(
                    key, reference_key
                )
            )
    return True


def _print_profile_batch_summary(results_by_key):
    print("\n{}".format("=" * 78))
    print("PROFILE C-INDEX SUMMARY")
    print("{}".format("=" * 78))
    print("{:<18s} {:<32s} {:>9s} {:>11s} {:>11s}".format(
        "Experiment", "Profile", "C-index", "CI lower", "CI upper"
    ))
    for (experiment, profile_name) in sorted(results_by_key):
        result = results_by_key[(experiment, profile_name)]
        print("{:<18s} {:<32s} {:>9.4f} {:>11.4f} {:>11.4f}".format(
            experiment,
            profile_name,
            float(result.c_index),
            float(result.confidence_interval[0]),
            float(result.confidence_interval[1]),
        ))


def run_profile_batch(args):
    specs = _batch_profile_specs(args)
    if not specs:
        raise ValueError("No batch profiles were requested.")
    data = load_main_input_data(args)
    _print_batch_dataset_summary(data, args, specs)
    outer_pair_plan = obtain_fixed_outer_pair_plan(data, args)
    results_by_key = {}
    for experiment, profile_name in specs:
        result = _run_one_batch_profile(
            data, args, outer_pair_plan, experiment, profile_name
        )
        results_by_key[(experiment, profile_name)] = result
        pickle_path = _save_one_batch_profile(
            result, args, experiment, profile_name
        )
        _export_one_batch_profile(
            result, args, experiment, profile_name
        )
        print("Saved profile pickle: {}".format(os.path.abspath(pickle_path)))

    _verify_batch_common_outer_pairs(results_by_key)
    root = os.path.abspath(args.output_dir)
    if not os.path.isdir(root):
        os.makedirs(root)
    outer_pair_plan.to_csv(
        os.path.join(root, "outer_pair_plan_used.csv"), index=False
    )
    _print_profile_batch_summary(results_by_key)
    mode = _profile_batch_mode(args)
    if mode in ("clinical_only", "all"):
        paths = _summary_for_profile_results(
            results_by_key, args, experiment="clinical_only"
        )
        if paths:
            print("Clinical profile summary PNG: {}".format(paths["png"]))
            print("Clinical profile summary TEX: {}".format(paths["tex"]))
            print("Clinical profile summary CSV: {}".format(paths["summary_csv"]))
    if mode in ("clinical_plus_pdr", "all"):
        paths = _summary_for_profile_results(
            results_by_key, args, experiment="clinical_plus_pdr"
        )
        if paths:
            print("PDR profile summary PNG: {}".format(paths["png"]))
            print("PDR profile summary TEX: {}".format(paths["tex"]))
            print("PDR profile summary CSV: {}".format(paths["summary_csv"]))
    if mode == "all":
        paths = _summary_for_profile_results(results_by_key, args, experiment=None)
        if paths:
            print("Combined profile summary PNG: {}".format(paths["png"]))
            print("Combined profile summary TEX: {}".format(paths["tex"]))
            print("Combined profile summary CSV: {}".format(paths["summary_csv"]))
    return results_by_key


def load_profile_batch(args):
    specs = _batch_profile_specs(args)
    if not specs:
        raise ValueError("No batch profiles were requested.")
    results_by_key = {}
    for experiment, profile_name in specs:
        result = _load_one_batch_profile(args, experiment, profile_name)
        results_by_key[(experiment, profile_name)] = result
        _export_one_batch_profile(result, args, experiment, profile_name)
    _verify_batch_common_outer_pairs(results_by_key)
    _print_profile_batch_summary(results_by_key)
    mode = _profile_batch_mode(args)
    if mode in ("clinical_only", "all"):
        _summary_for_profile_results(results_by_key, args, experiment="clinical_only")
    if mode in ("clinical_plus_pdr", "all"):
        _summary_for_profile_results(results_by_key, args, experiment="clinical_plus_pdr")
    if mode == "all":
        _summary_for_profile_results(results_by_key, args, experiment=None)
    return results_by_key


def print_model_profiles():
    print("CLINICAL-ONLY MODEL PROFILES")
    for name in sorted(CLINICAL_MODEL_PROFILES):
        profile = CLINICAL_MODEL_PROFILES[name]
        print("\n{}".format(name))
        print("  {}".format(profile["description"]))
        print("  groups: {}".format(", ".join(sorted(profile["groups"]))))
        for formula in profile["formulas"]:
            print("    - {}".format(formula))
    print("\nCLINICAL+PDR MODEL PROFILES")
    for name in sorted(CLINICAL_PDR_MODEL_PROFILES):
        profile = CLINICAL_PDR_MODEL_PROFILES[name]
        print("\n{}".format(name))
        print("  {}".format(profile["description"]))
        print("  groups: {}".format(", ".join(sorted(profile["groups"]))))
        if profile.get("consensus_note"):
            print("  note: {}".format(profile["consensus_note"]))
        for formula in profile["formulas"]:
            print("    - {}".format(formula))


def main():
    args = parse_arguments()
    if args.list_model_profiles:
        print_model_profiles()
        return
    if _profile_batch_mode(args) is not None:
        if args.mode == "run":
            run_profile_batch(args)
        else:
            load_profile_batch(args)
        return

    if args.mode == "run":
        clinical, clinical_pdr = run_and_pickle_results(args)
    else:
        clinical, clinical_pdr = load_pickled_results(args)

    print("\nCLINICAL-ONLY PERFORMANCE")
    print(clinical.summary_frame())
    print("Prediction aggregation: {} (weighting={}, rank_level={}, requested size={})".format(
        clinical.prediction_aggregation, clinical.ensemble_weighting,
        getattr(clinical, "ensemble_rank_level", "specification"), clinical.ensemble_size
    ))
    print("Bootstrap workers used: {}".format(getattr(clinical, "bootstrap_n_jobs_used", 1)))
    print("Consensus pair batch size: {}".format(clinical.consensus_pair_batch_size))
    print("Consensus filter scoring: {}".format(clinical.consensus_filter_scoring))
    print("Model profile: {}".format(getattr(clinical, "model_profile", args.clinical_only_model_profile)))
    print("Consensus top features: {}".format(clinical.consensus_top_features))
    print("\nCLINICAL + PDR PERFORMANCE")
    print(clinical_pdr.summary_frame())
    print("Prediction aggregation: {} (weighting={}, rank_level={}, requested size={})".format(
        clinical_pdr.prediction_aggregation, clinical_pdr.ensemble_weighting,
        getattr(clinical_pdr, "ensemble_rank_level", "specification"), clinical_pdr.ensemble_size
    ))
    print("Bootstrap workers used: {}".format(getattr(clinical_pdr, "bootstrap_n_jobs_used", 1)))
    print("Consensus pair batch size: {}".format(clinical_pdr.consensus_pair_batch_size))
    print("Consensus filter scoring: {}".format(clinical_pdr.consensus_filter_scoring))
    print("Model profile: {}".format(getattr(clinical_pdr, "model_profile", args.clinical_plus_pdr_model_profile)))
    print("Consensus top features: {}".format(clinical_pdr.consensus_top_features))
    export_all_graphics(clinical, clinical_pdr, args)


if __name__ == "__main__":
    main()
