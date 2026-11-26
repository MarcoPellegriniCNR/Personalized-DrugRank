"""Nested leave-pair-out evaluation for Cox proportional-hazards models.

The outer resampling requested here withholds two patients at a time.  This is
leave-pair-out (or sampled leave-pair-out) cross-validation rather than a
classical bootstrap sample.  Model selection is repeated *inside* every outer
training set, so the held-out patients are never used to choose the formula,
baseline model, number of spline knots, or penalization.

Independent outer pairs can be evaluated concurrently with Python's
multiprocessing library. The user selects ``n_jobs`` from 1 through 32; the
implementation never creates more than 32 worker processes. On Linux, CPU
affinity is also respected, and native BLAS/OpenMP libraries are limited to one
thread per worker to avoid oversubscription.

The pooled performance estimate is Harrell's concordance index computed from
within-pair rankings.  This is important because independently fitted Cox
models can put their risk scores on different numerical scales.

Inner model selection can use repeated shuffled cross-validation through
``inner_repeats`` or unique Harrell-comparable leave-pair-out validation.
Fold/pair-specific predictor standardization, reusable outer-pair plans, and
stability-aware candidate ranking are supported. Inner model selection can be
scored either by Harrell concordance or by Verweij-van Houwelingen grouped
cross-validated partial log-likelihood. Parvandeh-inspired grouped consensus
can be formed either from top-performing inner models or from training-only Cox
partial-likelihood feature/group rankings. Prediction can use the single best
configuration or an outer-training-only stacked/averaged ensemble of top
configurations. Repeat-to-repeat variability is retained as a stability diagnostic. ``inner_repeats=1`` reproduces the original
single-randomization inner-CV behavior.

Baseline estimation may use Breslow, spline, or piecewise-constant hazards.
For piecewise models, breakpoint locations are recomputed from training-set
duration quantiles inside every inner/outer fit to avoid information leakage.

Dependencies
------------
pandas, numpy, lifelines, scikit-learn, matplotlib (for PNG export)
Python standard-library pickle (for result persistence)

Outer-pair event policies
-------------------------
``event_policy`` controls which patients can enter an outer test pair:
``"all"`` uses all patient pairs; ``"both_events"`` uses only pairs in which
both patients experienced the event; ``"valid_pairs"`` uses pairs in which at
least one patient experienced the event; ``"harrell_comparable"`` uses only
pairs whose observed duration/event pattern is comparable under Harrell's
concordance rule. For the first three policies Harrell comparability is still
checked after model fitting.

Compatibility target
--------------------
This build targets exactly:
Python 3.6.8; lifelines 0.26.4; scikit-learn 0.24.2;
numpy 1.19.5; pandas 1.1.5; scipy 1.5.4; matplotlib 3.3.4 for PNG export.

It intentionally avoids dataclasses, postponed annotations, PEP 585 built-in
generics, and math.comb so that Python 3.6.8 is sufficient. Multiprocessing is
implemented with module-level worker functions and ``multiprocessing.Pool``.

Selected references
-------------------
Cox DR (1972), Journal of the Royal Statistical Society B, 34:187-202.
Harrell FE Jr, Lee KL, Mark DB (1996), Statistics in Medicine, 15:361-387.
Royston P, Parmar MKB (2002), Statistics in Medicine, 21:2175-2197.
Varma S, Simon R (2006), BMC Bioinformatics, 7:91.
Efron B (1979), Annals of Statistics, 7:1-26.
Uno H et al. (2011), Statistics in Medicine, 30:1105-1117.
Verweij PJM, van Houwelingen HC (1993), Statistics in Medicine, 12:2305-2314.
Parvandeh S et al. (2020), Bioinformatics, 36:3093-3098.
Wey A, Connett J, Rudser K (2015), Biostatistics, 16:537-549.
"""


from functools import wraps
from itertools import combinations
import multiprocessing as mp
from numbers import Integral
import os
import pickle
import re
import time
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple
import warnings

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, ParameterGrid, StratifiedKFold

try:
    import matplotlib
    # Use a non-interactive backend so PNG export works on headless CentOS nodes.
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib = None
    plt = None

try:
    from threadpoolctl import threadpool_limits
except ImportError:
    # scikit-learn 0.24.2 normally installs threadpoolctl as a dependency.
    # The module still works without it, but native BLAS/OpenMP thread counts
    # cannot then be changed dynamically after those libraries are loaded.
    threadpool_limits = None

try:
    from lifelines import CoxPHFitter
except ImportError as exc:  # Allows importing/testing non-fitting helpers.
    CoxPHFitter = None  # type: ignore[assignment]
    _LIFELINES_IMPORT_ERROR: Optional[ImportError] = exc
else:
    _LIFELINES_IMPORT_ERROR = None


MAX_MULTIPROCESSING_WORKERS = 32
_OUTER_WORKER_STATE = None
_WORKER_THREADPOOL_LIMITER = None


class CoxModelSpec(object):
    """One candidate in the inner model-selection grid.

    Implemented as a regular class rather than a dataclass so this module has
    no dependency on the Python 3.7+ :mod:`dataclasses` standard-library module.
    """

    __slots__ = (
        "formula",
        "base_model",
        "knot_number",
        "piecewise_break_number",
        "penalizer",
        "l1_ratio",
        "order",
    )

    def __init__(
        self,
        formula,
        base_model,
        knot_number,
        piecewise_break_number,
        penalizer,
        l1_ratio,
        order,
    ):
        self.formula = formula
        self.base_model = base_model
        self.knot_number = knot_number
        self.piecewise_break_number = piecewise_break_number
        self.penalizer = penalizer
        self.l1_ratio = l1_ratio
        self.order = order

    @property
    def formula_label(self):
        return "<all covariates>" if self.formula is None else self.formula

    @property
    def label(self):
        knot = "-" if self.knot_number is None else str(self.knot_number)
        piecewise_breaks = (
            "-"
            if self.piecewise_break_number is None
            else str(self.piecewise_break_number)
        )
        return (
            "formula={}; base_model={}; knots={}; piecewise_breaks={}; "
            "penalizer={:g}; l1_ratio={:g}"
        ).format(
            self.formula_label,
            self.base_model,
            knot,
            piecewise_breaks,
            self.penalizer,
            self.l1_ratio,
        )

    def __repr__(self):
        return "CoxModelSpec({})".format(self.label)

    def __reduce__(self):
        """Provide an explicit Python-3.6-compatible pickle representation."""

        return (
            self.__class__,
            (
                self.formula,
                self.base_model,
                self.knot_number,
                self.piecewise_break_number,
                self.penalizer,
                self.l1_ratio,
                self.order,
            ),
        )


class CoxPredictionEnsemble(object):
    """Leakage-safe prediction ensemble fitted on one training set.

    Each component Cox model is put on a common scale by centering and scaling
    its log-relative-risk score using *training-set predictions only*. The
    transformed component scores are then averaged with fixed non-negative
    weights. This deliberately implements prediction averaging rather than
    estimating another set of stacking coefficients from the tiny outer
    training cohort.
    """

    def __init__(
        self,
        models,
        specs,
        weights,
        score_centers,
        score_scales,
        covariates,
        duration_col,
        event_col,
        training_data,
        weighting,
        rank_level="specification",
    ):
        self.models = list(models)
        self.specs = list(specs)
        self.weights = np.asarray(weights, dtype=float)
        self.score_centers = np.asarray(score_centers, dtype=float)
        self.score_scales = np.asarray(score_scales, dtype=float)
        self.covariates = list(covariates)
        self.duration_col = str(duration_col)
        self.event_col = str(event_col)
        self.training_data = training_data.copy()
        self.weighting = str(weighting)
        self.rank_level = str(rank_level)

    def predict_log_risk(self, data):
        component_scores = []
        for model, center, scale in zip(
            self.models, self.score_centers, self.score_scales
        ):
            raw = _predict_log_relative_risk(
                model,
                train_data=self.training_data,
                evaluation_data=data,
                covariates=self.covariates,
                duration_col=self.duration_col,
                event_col=self.event_col,
            )
            component_scores.append((np.asarray(raw, dtype=float) - center) / scale)
        matrix = np.vstack(component_scores)
        return np.dot(self.weights, matrix)

    def predict_risk(self, data):
        """Alias returning the ensemble score; larger values mean higher risk."""
        return self.predict_log_risk(data)


class CoxEvaluationResult(object):
    """Outputs from :func:`nested_pair_cox_evaluation`."""

    def __init__(
        self,
        c_index,
        confidence_interval,
        confidence_level,
        bootstrap_standard_error,
        n_patients,
        total_possible_pairs,
        total_policy_eligible_pairs,
        event_policy,
        n_outer_pairs_requested,
        n_outer_pairs_fitted,
        n_comparable_pairs,
        n_failed_outer_pairs,
        n_jobs_requested,
        n_jobs_used,
        multiprocessing_start_method,
        inner_repeats,
        inner_validation_method,
        inner_scoring,
        inner_selection_rule,
        inner_stability_penalty,
        standardize_predictors,
        standardized_columns,
        consensus_method,
        consensus_top_models,
        consensus_top_features,
        consensus_threshold,
        consensus_max_formula_groups,
        consensus_pair_batch_size,
        consensus_filter_scoring,
        consensus_feature_groups,
        consensus_group_expressions,
        outer_pairs_reused,
        outer_pair_plan,
        outer_results,
        selection_frequency,
        bootstrap_c_indices,
        final_model,
        final_spec,
        final_inner_ranking,
        prediction_aggregation="winner",
        ensemble_size=1,
        ensemble_weighting="equal",
        ensemble_temperature=1.0,
        ensemble_distinct_formulas=True,
        ensemble_rank_level="specification",
        bootstrap_n_jobs_used=1,
    ):
        self.c_index = c_index
        self.confidence_interval = confidence_interval
        self.confidence_level = confidence_level
        self.bootstrap_standard_error = bootstrap_standard_error
        self.n_patients = n_patients
        self.total_possible_pairs = total_possible_pairs
        self.total_policy_eligible_pairs = total_policy_eligible_pairs
        self.event_policy = event_policy
        self.n_outer_pairs_requested = n_outer_pairs_requested
        self.n_outer_pairs_fitted = n_outer_pairs_fitted
        self.n_comparable_pairs = n_comparable_pairs
        self.n_failed_outer_pairs = n_failed_outer_pairs
        self.n_jobs_requested = n_jobs_requested
        self.n_jobs_used = n_jobs_used
        self.multiprocessing_start_method = multiprocessing_start_method
        self.inner_repeats = inner_repeats
        self.inner_validation_method = inner_validation_method
        self.inner_scoring = inner_scoring
        self.inner_selection_rule = inner_selection_rule
        self.inner_stability_penalty = inner_stability_penalty
        self.standardize_predictors = standardize_predictors
        self.standardized_columns = list(standardized_columns)
        self.consensus_method = consensus_method
        self.consensus_top_models = int(consensus_top_models)
        self.consensus_top_features = int(consensus_top_features)
        self.consensus_threshold = float(consensus_threshold)
        self.consensus_max_formula_groups = int(consensus_max_formula_groups)
        self.consensus_pair_batch_size = int(consensus_pair_batch_size)
        self.consensus_filter_scoring = str(consensus_filter_scoring)
        self.consensus_feature_groups = (
            None if consensus_feature_groups is None
            else {str(k): tuple(v) for k, v in consensus_feature_groups.items()}
        )
        self.consensus_group_expressions = (
            None if consensus_group_expressions is None
            else {str(k): str(v) for k, v in consensus_group_expressions.items()}
        )
        self.outer_pairs_reused = bool(outer_pairs_reused)
        self.outer_pair_plan = outer_pair_plan
        self.outer_results = outer_results
        self.selection_frequency = selection_frequency
        self.bootstrap_c_indices = bootstrap_c_indices
        self.final_model = final_model
        self.final_spec = final_spec
        self.final_inner_ranking = final_inner_ranking
        self.prediction_aggregation = str(prediction_aggregation)
        self.ensemble_size = int(ensemble_size)
        self.ensemble_weighting = str(ensemble_weighting)
        self.ensemble_temperature = float(ensemble_temperature)
        self.ensemble_distinct_formulas = bool(ensemble_distinct_formulas)
        self.ensemble_rank_level = str(ensemble_rank_level)
        self.bootstrap_n_jobs_used = int(bootstrap_n_jobs_used)

    def summary_frame(self):
        """Return the main result as a one-row DataFrame."""

        return pd.DataFrame(
            {
                "c_index": [self.c_index],
                "ci_lower": [self.confidence_interval[0]],
                "ci_upper": [self.confidence_interval[1]],
                "confidence_level": [self.confidence_level],
                "bootstrap_standard_error": [self.bootstrap_standard_error],
                "n_patients": [self.n_patients],
                "total_possible_pairs": [self.total_possible_pairs],
                "total_policy_eligible_pairs": [self.total_policy_eligible_pairs],
                "event_policy": [self.event_policy],
                "n_outer_pairs_requested": [self.n_outer_pairs_requested],
                "n_outer_pairs_fitted": [self.n_outer_pairs_fitted],
                "n_comparable_pairs": [self.n_comparable_pairs],
                "n_failed_outer_pairs": [self.n_failed_outer_pairs],
                "n_jobs_requested": [self.n_jobs_requested],
                "n_jobs_used": [self.n_jobs_used],
                "multiprocessing_start_method": [self.multiprocessing_start_method],
                "inner_repeats": [self.inner_repeats],
                "inner_validation_method": [self.inner_validation_method],
                "inner_scoring": [self.inner_scoring],
                "inner_selection_rule": [self.inner_selection_rule],
                "inner_stability_penalty": [self.inner_stability_penalty],
                "standardize_predictors": [self.standardize_predictors],
                "standardized_columns": [",".join(self.standardized_columns)],
                "consensus_method": [self.consensus_method],
                "consensus_top_models": [self.consensus_top_models],
                "consensus_top_features": [self.consensus_top_features],
                "consensus_threshold": [self.consensus_threshold],
                "consensus_max_formula_groups": [self.consensus_max_formula_groups],
                "consensus_pair_batch_size": [self.consensus_pair_batch_size],
                "consensus_filter_scoring": [self.consensus_filter_scoring],
                "consensus_groups_defined": [
                    "" if self.consensus_feature_groups is None
                    else ",".join(sorted(self.consensus_feature_groups))
                ],
                "outer_pairs_reused": [self.outer_pairs_reused],
                "prediction_aggregation": [self.prediction_aggregation],
                "ensemble_size": [self.ensemble_size],
                "ensemble_weighting": [self.ensemble_weighting],
                "ensemble_temperature": [self.ensemble_temperature],
                "ensemble_distinct_formulas": [self.ensemble_distinct_formulas],
                "ensemble_rank_level": [self.ensemble_rank_level],
                "bootstrap_n_jobs_used": [self.bootstrap_n_jobs_used],
                "model_profile": [getattr(self, "model_profile", "")],
                "experiment_name": [getattr(self, "experiment_name", "")],
            },
            index=["nested_leave_pair_out_cox"],
        )

    def to_pickle(self, output_path, include_final_model=True, protocol=4):
        """Serialize this evaluation result to a pickle file."""
        return save_evaluation_result_pickle(
            self,
            output_path=output_path,
            include_final_model=include_final_model,
            protocol=protocol,
        )

    @classmethod
    def from_pickle(cls, input_path):
        """Load a :class:`CoxEvaluationResult` from a trusted pickle file."""
        result = load_evaluation_result_pickle(input_path)
        if not isinstance(result, cls):
            raise TypeError(
                "Pickle did not contain a CoxEvaluationResult compatible with this module."
            )
        return result

    def export_graphics(
        self,
        output_dir,
        prefix="cox_performance",
        top_n_models=10,
        histogram_bins=20,
        png_dpi=300,
    ):
        """Export each diagnostic as PGFPlots ``.tex``, PNG, and CSV data."""
        return export_performance_graphics(
            self,
            output_dir=output_dir,
            prefix=prefix,
            top_n_models=top_n_models,
            histogram_bins=histogram_bins,
            png_dpi=png_dpi,
        )

    def export_latex_graphics(
        self,
        output_dir,
        prefix="cox_performance",
        top_n_models=10,
        histogram_bins=20,
        png_dpi=300,
    ):
        """Backward-compatible alias that now exports both ``.tex`` and PNG."""
        return self.export_graphics(
            output_dir=output_dir,
            prefix=prefix,
            top_n_models=top_n_models,
            histogram_bins=histogram_bins,
            png_dpi=png_dpi,
        )

_PICKLE_FORMAT_NAME = "nested_pair_cox_evaluation_result"
_PICKLE_FORMAT_VERSION = 9


def _result_copy_for_pickle(result, include_final_model):
    """Return a result copy suitable for serialization."""

    if include_final_model:
        return result

    stored = CoxEvaluationResult(
        c_index=result.c_index,
        confidence_interval=result.confidence_interval,
        confidence_level=result.confidence_level,
        bootstrap_standard_error=result.bootstrap_standard_error,
        n_patients=result.n_patients,
        total_possible_pairs=result.total_possible_pairs,
        total_policy_eligible_pairs=result.total_policy_eligible_pairs,
        event_policy=result.event_policy,
        n_outer_pairs_requested=result.n_outer_pairs_requested,
        n_outer_pairs_fitted=result.n_outer_pairs_fitted,
        n_comparable_pairs=result.n_comparable_pairs,
        n_failed_outer_pairs=result.n_failed_outer_pairs,
        n_jobs_requested=result.n_jobs_requested,
        n_jobs_used=result.n_jobs_used,
        multiprocessing_start_method=result.multiprocessing_start_method,
        inner_repeats=result.inner_repeats,
        inner_validation_method=result.inner_validation_method,
        inner_scoring=result.inner_scoring,
        inner_selection_rule=result.inner_selection_rule,
        inner_stability_penalty=result.inner_stability_penalty,
        standardize_predictors=result.standardize_predictors,
        standardized_columns=result.standardized_columns,
        consensus_method=result.consensus_method,
        consensus_top_models=result.consensus_top_models,
        consensus_top_features=result.consensus_top_features,
        consensus_threshold=result.consensus_threshold,
        consensus_max_formula_groups=result.consensus_max_formula_groups,
        consensus_pair_batch_size=getattr(result, "consensus_pair_batch_size", 1),
        consensus_filter_scoring=getattr(
            result, "consensus_filter_scoring", "partial_likelihood"
        ),
        consensus_feature_groups=result.consensus_feature_groups,
        consensus_group_expressions=result.consensus_group_expressions,
        outer_pairs_reused=result.outer_pairs_reused,
        outer_pair_plan=result.outer_pair_plan,
        outer_results=result.outer_results,
        selection_frequency=result.selection_frequency,
        bootstrap_c_indices=result.bootstrap_c_indices,
        final_model=None,
        final_spec=result.final_spec,
        final_inner_ranking=result.final_inner_ranking,
        prediction_aggregation=getattr(result, "prediction_aggregation", "winner"),
        ensemble_size=getattr(result, "ensemble_size", 1),
        ensemble_weighting=getattr(result, "ensemble_weighting", "equal"),
        ensemble_temperature=getattr(result, "ensemble_temperature", 1.0),
        ensemble_distinct_formulas=getattr(result, "ensemble_distinct_formulas", True),
        ensemble_rank_level=getattr(result, "ensemble_rank_level", "specification"),
        bootstrap_n_jobs_used=getattr(result, "bootstrap_n_jobs_used", 1),
    )
    for attribute_name in ("model_profile", "experiment_name"):
        if hasattr(result, attribute_name):
            setattr(stored, attribute_name, getattr(result, attribute_name))
    return stored

def save_evaluation_result_pickle(
    result,
    output_path,
    include_final_model=True,
    protocol=4,
):
    """Save a :class:`CoxEvaluationResult` to a versioned pickle payload.

    Protocol 4 is used by default because it is the highest protocol supported
    by Python 3.6.8. For later graphical regeneration the fitted final Cox model
    is not needed, so callers may set ``include_final_model=False``.

    Notes
    -----
    Pickle is Python-specific and should be treated as a trusted-data format.
    Do not load pickle files obtained from untrusted sources.
    """

    if not isinstance(result, CoxEvaluationResult):
        raise TypeError("result must be a CoxEvaluationResult instance.")

    if isinstance(protocol, bool) or not isinstance(protocol, Integral):
        raise TypeError("protocol must be an integer pickle protocol number.")
    protocol = int(protocol)
    if protocol < 0 or protocol > pickle.HIGHEST_PROTOCOL:
        raise ValueError(
            "protocol must be between 0 and {} for this Python runtime.".format(
                pickle.HIGHEST_PROTOCOL
            )
        )

    output_path = os.path.abspath(str(output_path))
    parent = os.path.dirname(output_path)
    if parent and not os.path.isdir(parent):
        os.makedirs(parent)

    stored_result = _result_copy_for_pickle(
        result,
        include_final_model=bool(include_final_model),
    )
    payload = {
        "format": _PICKLE_FORMAT_NAME,
        "format_version": _PICKLE_FORMAT_VERSION,
        "pickle_protocol": protocol,
        "include_final_model": bool(include_final_model),
        "result": stored_result,
    }

    try:
        with open(output_path, "wb") as handle:
            pickle.dump(payload, handle, protocol=protocol)
    except Exception as exc:
        # Remove a partially written file so it is not mistaken for a valid
        # result dump on a later graphics-only run.
        try:
            if os.path.isfile(output_path):
                os.remove(output_path)
        except OSError:
            pass
        if include_final_model:
            raise RuntimeError(
                "Could not pickle the complete result. The fitted final_model "
                "may contain state that is not serializable in this environment. "
                "For graphics-only persistence, retry with "
                "include_final_model=False. Original error: {}".format(exc)
            )
        raise

    return output_path


def load_evaluation_result_pickle(input_path):
    """Load a previously saved evaluation result from a trusted pickle file.

    The returned object can immediately be passed to ``export_graphics`` or to
    ``export_cindex_comparison_graphics`` without rerunning Cox fitting.

    Security
    --------
    Python pickle files can execute arbitrary code while loading. Only load
    files that you created yourself or otherwise trust completely.
    """

    input_path = os.path.abspath(str(input_path))
    if not os.path.isfile(input_path):
        raise IOError("Pickle file does not exist: {}".format(input_path))

    with open(input_path, "rb") as handle:
        payload = pickle.load(handle)

    # Accept a direct CoxEvaluationResult for convenience with very early
    # experimental dumps, while preferring the versioned payload above.
    if isinstance(payload, CoxEvaluationResult):
        return payload

    if not isinstance(payload, dict):
        raise ValueError("Unrecognized Cox evaluation pickle payload.")
    if payload.get("format") != _PICKLE_FORMAT_NAME:
        raise ValueError(
            "Unrecognized pickle format: {!r}.".format(payload.get("format"))
        )
    if payload.get("format_version") != _PICKLE_FORMAT_VERSION:
        raise ValueError(
            "Unsupported pickle format version: {!r}.".format(
                payload.get("format_version")
            )
        )

    result = payload.get("result")
    if not isinstance(result, CoxEvaluationResult):
        raise TypeError(
            "Pickle payload does not contain a CoxEvaluationResult compatible "
            "with this module."
        )
    if not hasattr(result, "consensus_filter_scoring"):
        result.consensus_filter_scoring = "partial_likelihood"
    return result


def _require_lifelines() -> None:
    if CoxPHFitter is None:
        raise ImportError(
            "lifelines is required. This compatibility build targets Python 3.6.8, "
            "lifelines 0.26.4, scikit-learn 0.24.2, numpy 1.19.5, "
            "pandas 1.1.5, and scipy 1.5.4."
        ) from _LIFELINES_IMPORT_ERROR


def _validate_and_prepare_data(
    data: pd.DataFrame,
    duration_col: str,
    event_col: str,
    exclude_covariates: Sequence[str],
    stabilize_categorical_levels: bool,
) -> Tuple[pd.DataFrame, List[str], np.ndarray]:
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame.")
    if duration_col not in data.columns or event_col not in data.columns:
        raise KeyError(
            f"data must contain {duration_col!r} and {event_col!r}."
        )
    if duration_col == event_col:
        raise ValueError("duration_col and event_col must be different columns.")

    excluded = list(dict.fromkeys(exclude_covariates))
    forbidden = {duration_col, event_col}.intersection(excluded)
    if forbidden:
        raise ValueError(
            "Response columns cannot be excluded: " + ", ".join(sorted(forbidden))
        )
    unknown_excluded = sorted(set(excluded).difference(data.columns))
    if unknown_excluded:
        raise KeyError(f"Unknown excluded covariates: {unknown_excluded}")

    original_index = data.index.to_numpy(copy=True)
    covariates = [
        c
        for c in data.columns
        if c not in {duration_col, event_col} and c not in excluded
    ]
    if not covariates:
        raise ValueError("No covariate columns remain after exclusions.")

    work = data[[duration_col, event_col, *covariates]].copy().reset_index(drop=True)
    if len(work) < 6:
        raise ValueError(
            "At least 6 patients are required for an outer test pair plus "
            "a minimally useful inner cross-validation."
        )

    try:
        work[duration_col] = pd.to_numeric(work[duration_col], errors="raise").astype(float)
        event_numeric = pd.to_numeric(work[event_col], errors="raise")
    except (TypeError, ValueError) as exc:
        raise ValueError("Duration and event columns must be numeric.") from exc

    if work[duration_col].isna().any() or event_numeric.isna().any():
        raise ValueError("Duration and event columns cannot contain missing values.")
    if not np.isfinite(work[duration_col].to_numpy(dtype=float)).all():
        raise ValueError("Durations must be finite.")
    if (work[duration_col] <= 0).any():
        raise ValueError("All durations must be strictly positive.")

    unique_events = set(pd.unique(event_numeric.dropna()))
    if not unique_events.issubset({0, 1, False, True}):
        raise ValueError("The event indicator must contain only 0/1 values.")
    work[event_col] = event_numeric.astype(int)
    if work[event_col].sum() < 2:
        raise ValueError("At least two observed events are required.")

    if work[covariates].isna().any().any():
        missing = work[covariates].columns[work[covariates].isna().any()].tolist()
        raise ValueError(
            "Covariates contain missing values. Impute them inside the nested "
            f"resampling procedure before fitting. Columns: {missing}"
        )

    numeric_covariates = work[covariates].select_dtypes(include=[np.number])
    if not numeric_covariates.empty:
        numeric_values = numeric_covariates.to_numpy(dtype=float)
        if not np.isfinite(numeric_values).all():
            raise ValueError("Numeric covariates must be finite.")

    if stabilize_categorical_levels:
        for column in covariates:
            dtype = work[column].dtype
            if (
                pd.api.types.is_object_dtype(dtype)
                or pd.api.types.is_string_dtype(dtype)
                or pd.api.types.is_bool_dtype(dtype)
            ):
                work[column] = pd.Categorical(work[column])

    constant_columns = [c for c in covariates if work[c].nunique(dropna=False) <= 1]
    if constant_columns:
        warnings.warn(
            "Constant covariates can make Cox fits singular and should normally "
            f"be excluded: {constant_columns}",
            RuntimeWarning,
            stacklevel=2,
        )

    return work, covariates, original_index


def _build_candidate_specs(
    formulas: Sequence[Optional[str]],
    base_models: Sequence[str],
    knot_numbers: Sequence[int],
    piecewise_break_numbers: Sequence[int],
    penalizers: Sequence[float],
    l1_ratios: Sequence[float],
) -> List[CoxModelSpec]:
    if not formulas:
        raise ValueError("formulas cannot be empty. Use [None] for all covariates.")
    if not base_models:
        raise ValueError("base_models cannot be empty.")

    allowed = {"breslow", "spline", "piecewise"}
    invalid_models = sorted(set(base_models).difference(allowed))
    if invalid_models:
        raise ValueError(
            f"Unsupported base_models {invalid_models}; allowed values are {sorted(allowed)}."
        )
    if "spline" in base_models:
        if not knot_numbers:
            raise ValueError("knot_numbers cannot be empty when using a spline model.")
        if any(int(k) != k or int(k) < 2 for k in knot_numbers):
            raise ValueError("Every spline knot number must be an integer >= 2.")
    if "piecewise" in base_models:
        if not piecewise_break_numbers:
            raise ValueError(
                "piecewise_break_numbers cannot be empty when using a piecewise model."
            )
        if any(int(k) != k or int(k) < 1 for k in piecewise_break_numbers):
            raise ValueError(
                "Every piecewise break number must be an integer >= 1."
            )
    if any(float(p) < 0 for p in penalizers):
        raise ValueError("penalizers must be non-negative.")
    if any(not 0 <= float(x) <= 1 for x in l1_ratios):
        raise ValueError("l1_ratios must be between 0 and 1.")

    specs: List[CoxModelSpec] = []
    order = 0
    for base_model in base_models:
        if base_model == "breslow":
            knot_values = [None]
            piecewise_break_values = [None]
        elif base_model == "spline":
            knot_values = [int(k) for k in knot_numbers]
            piecewise_break_values = [None]
        else:
            knot_values = [None]
            piecewise_break_values = [int(k) for k in piecewise_break_numbers]

        conditional_grid: Dict[str, Sequence[Any]] = {
            "formula": list(formulas),
            "penalizer": [float(p) for p in penalizers],
            "l1_ratio": [float(x) for x in l1_ratios],
            "knot_number": knot_values,
            "piecewise_break_number": piecewise_break_values,
        }
        for params in ParameterGrid(conditional_grid):
            specs.append(
                CoxModelSpec(
                    formula=params["formula"],
                    base_model=base_model,
                    knot_number=params["knot_number"],
                    piecewise_break_number=params["piecewise_break_number"],
                    penalizer=params["penalizer"],
                    l1_ratio=params["l1_ratio"],
                    order=order,
                )
            )
            order += 1
    return specs


def _number_of_pairs(n):
    """Return n choose 2 without requiring math.comb (Python >=3.8)."""
    n = int(n)
    if n < 0:
        raise ValueError("n must be non-negative.")
    return n * (n - 1) // 2


def _validate_event_policy(event_policy):
    """Validate and normalize the outer-pair event-selection policy."""

    if not isinstance(event_policy, str):
        raise TypeError("event_policy must be a string.")
    allowed = ("all", "both_events", "valid_pairs", "harrell_comparable")
    if event_policy not in allowed:
        raise ValueError(
            "event_policy must be one of {}.".format(", ".join(repr(v) for v in allowed))
        )
    return event_policy


def _validate_event_values(event_values, n_patients):
    """Return a validated 0/1 integer event array of length n_patients."""

    if event_values is None:
        raise ValueError(
            "event_values must be supplied when event_policy is not 'all'."
        )
    values = np.asarray(event_values, dtype=int).reshape(-1)
    if len(values) != int(n_patients):
        raise ValueError("event_values length does not match n_patients.")
    unique = set(values.tolist())
    if not unique.issubset({0, 1}):
        raise ValueError("event_values must contain only 0/1 values.")
    return values


def _validate_duration_values(duration_values, n_patients):
    """Return a validated finite duration array of length n_patients."""

    if duration_values is None:
        raise ValueError(
            "duration_values must be supplied for event_policy='harrell_comparable'."
        )
    values = np.asarray(duration_values, dtype=float).reshape(-1)
    if len(values) != int(n_patients):
        raise ValueError("duration_values length does not match n_patients.")
    if not np.isfinite(values).all():
        raise ValueError("duration_values must be finite.")
    return values


def _is_harrell_comparable_pair(duration_i, event_i, duration_j, event_j):
    """Return True exactly when a pair is comparable for Harrell's C-index.

    If observed times differ, the patient with the shorter observed time must
    have experienced the event. If observed times are equal, the current
    evaluator treats only event/censor pairs as informative; event/event and
    censor/censor ties are not comparable. This is deliberately identical to
    the comparability rule used later by ``_pair_contribution``.
    """

    if duration_i < duration_j:
        return int(event_i) == 1
    if duration_j < duration_i:
        return int(event_j) == 1
    return int(event_i) != int(event_j)


def _count_policy_eligible_pairs(
    n_patients,
    event_values,
    event_policy,
    duration_values=None,
):
    """Count all pairs eligible under the requested event policy."""

    event_policy = _validate_event_policy(event_policy)
    n_patients = int(n_patients)
    if event_policy == "all":
        return _number_of_pairs(n_patients)

    values = _validate_event_values(event_values, n_patients)
    n_events = int(values.sum())
    n_censored = n_patients - n_events

    if event_policy == "both_events":
        return _number_of_pairs(n_events)

    if event_policy == "valid_pairs":
        # At least one member of the pair experienced the event.
        return _number_of_pairs(n_patients) - _number_of_pairs(n_censored)

    durations = _validate_duration_values(duration_values, n_patients)
    comparable_count = 0
    for i, j in combinations(range(n_patients), 2):
        if _is_harrell_comparable_pair(
            float(durations[i]),
            int(values[i]),
            float(durations[j]),
            int(values[j]),
        ):
            comparable_count += 1
    return comparable_count


def _all_policy_eligible_pairs(
    n_patients,
    event_values,
    event_policy,
    duration_values=None,
):
    """Enumerate every eligible pair under an event-selection policy."""

    event_policy = _validate_event_policy(event_policy)
    all_indices = list(range(int(n_patients)))
    if event_policy == "all":
        return list(combinations(all_indices, 2))

    values = _validate_event_values(event_values, n_patients)
    event_indices = np.flatnonzero(values == 1).astype(int).tolist()
    if event_policy == "both_events":
        return list(combinations(event_indices, 2))

    if event_policy == "valid_pairs":
        # Exclude only censored-censored pairs at generation time.
        return [
            (i, j)
            for i, j in combinations(all_indices, 2)
            if values[i] == 1 or values[j] == 1
        ]

    durations = _validate_duration_values(duration_values, n_patients)
    return [
        (i, j)
        for i, j in combinations(all_indices, 2)
        if _is_harrell_comparable_pair(
            float(durations[i]),
            int(values[i]),
            float(durations[j]),
            int(values[j]),
        )
    ]


def _pair_is_eligible_under_policy(
    i,
    j,
    event_values,
    duration_values,
    event_policy,
):
    """Return whether one concrete pair satisfies an outer-pair policy."""

    event_policy = _validate_event_policy(event_policy)
    if event_policy == "all":
        return True
    event_i = int(event_values[i])
    event_j = int(event_values[j])
    if event_policy == "both_events":
        return event_i == 1 and event_j == 1
    if event_policy == "valid_pairs":
        return event_i == 1 or event_j == 1
    return _is_harrell_comparable_pair(
        float(duration_values[i]),
        event_i,
        float(duration_values[j]),
        event_j,
    )


def _canonical_patient_id(value):
    """Return a stable text representation for patient identifiers.

    CSV round-trips can otherwise turn an integer identifier such as ``0``
    into the floating-point value ``0.0`` when a heterogeneous pandas row is
    materialized.  Identifiers are therefore compared in canonical text form.
    """
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass

    if isinstance(value, Integral):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        numeric = float(value)
        if np.isfinite(numeric) and numeric.is_integer():
            return str(int(numeric))
        if np.isfinite(numeric):
            return "{:.17g}".format(numeric)
        return str(numeric)

    text = str(value).strip()
    # Defensive support for legacy CSVs whose identifier column was parsed as
    # text after first becoming a floating-point integer, e.g. "12.0".
    if re.match(r"^[+-]?\d+\.0+$", text):
        try:
            return str(int(float(text)))
        except (TypeError, ValueError, OverflowError):
            pass
    return text


def _identifier_values_are_usable(values):
    canonical = [_canonical_patient_id(value) for value in values]
    if any(value is None or value == "" for value in canonical):
        return False, canonical
    if len(set(canonical)) != len(canonical):
        return False, canonical
    return True, canonical


def _patient_identifier_sources(work, original_index):
    """Return unique patient-ID sources available in the current data.

    The original DataFrame index is always considered.  When a conventional
    ``Patient_ID`` column (or close spelling variant) is present, unique and
    non-missing, it is also available and is preferred for newly written pair
    plans.  This lets TCGA barcodes remain stable if rows are reordered.
    """
    sources = []

    index_ok, index_values = _identifier_values_are_usable(list(original_index))
    if index_ok:
        sources.append(("__index__", index_values))

    candidate_columns = (
        "Patient_ID",
        "patient_id",
        "PatientID",
        "patientID",
    )
    for column in candidate_columns:
        if column not in work.columns:
            continue
        ok, values = _identifier_values_are_usable(work[column].tolist())
        if ok:
            # Put explicit patient-ID columns before the DataFrame index so
            # they are preferred for newly generated plans.
            sources.insert(0, (column, values))
            break
        warnings.warn(
            "Column {!r} exists but cannot be used as a reusable outer-pair "
            "identifier because it contains missing or duplicate values. "
            "Falling back to another identifier source.".format(column),
            RuntimeWarning,
            stacklevel=3,
        )
        break

    if not sources:
        raise ValueError(
            "No unique, non-missing patient identifier source is available. "
            "Use a unique DataFrame index or a unique Patient_ID column."
        )
    return sources


def _preferred_patient_identifiers(work, original_index):
    sources = _patient_identifier_sources(work, original_index)
    return sources[0][0], sources[0][1]


def _match_pair_plan_identifier_source(plan_frame, work, original_index):
    """Match stored pair-plan IDs to the current DataFrame and return positions."""
    if "patient_id_i" not in plan_frame.columns or "patient_id_j" not in plan_frame.columns:
        return None, None

    plan_i = [_canonical_patient_id(value) for value in plan_frame["patient_id_i"].tolist()]
    plan_j = [_canonical_patient_id(value) for value in plan_frame["patient_id_j"].tolist()]
    if any(value is None or value == "" for value in plan_i + plan_j):
        raise ValueError("Reusable outer-pair plan contains missing patient identifiers.")

    sources = _patient_identifier_sources(work, original_index)
    source_by_name = dict(sources)

    hinted_source = None
    if "patient_id_source" in plan_frame.columns:
        hints = [
            str(value).strip()
            for value in plan_frame["patient_id_source"].dropna().tolist()
            if str(value).strip()
        ]
        unique_hints = sorted(set(hints))
        if len(unique_hints) > 1:
            raise ValueError(
                "Reusable outer-pair plan contains multiple patient_id_source values: {}".format(
                    unique_hints
                )
            )
        if unique_hints:
            hinted_source = unique_hints[0]

    ordered_sources = list(sources)
    if hinted_source is not None:
        if hinted_source not in source_by_name:
            raise ValueError(
                "Reusable outer-pair plan expects patient identifier source {!r}, "
                "but the current DataFrame does not provide a unique usable source "
                "with that name.".format(hinted_source)
            )
        ordered_sources = [(hinted_source, source_by_name[hinted_source])] + [
            item for item in sources if item[0] != hinted_source
        ]

    all_plan_ids = set(plan_i + plan_j)
    matched = None
    for source_name, current_ids in ordered_sources:
        current_set = set(current_ids)
        if all_plan_ids.issubset(current_set):
            matched = (source_name, current_ids)
            break

    if matched is None:
        sample_ids = sorted(all_plan_ids)[:5]
        available = [name for name, _values in sources]
        raise ValueError(
            "Reusable outer-pair plan patient IDs cannot be matched to the current "
            "DataFrame. Example stored IDs: {}. Available identifier sources: {}. "
            "This usually means the pair plan belongs to a different cohort, or the "
            "patient identifier column/index changed.".format(sample_ids, available)
        )

    source_name, current_ids = matched
    position_by_id = {patient_id: position for position, patient_id in enumerate(current_ids)}
    mapped_pairs = [
        (position_by_id[patient_i], position_by_id[patient_j])
        for patient_i, patient_j in zip(plan_i, plan_j)
    ]
    return source_name, mapped_pairs


def _outer_pair_plan_from_pairs(
    pairs,
    work,
    original_index,
    duration_col,
    event_col,
):
    patient_id_source, patient_ids = _preferred_patient_identifiers(work, original_index)
    rows = []
    for pair_number, pair in enumerate(pairs, start=1):
        i, j = int(pair[0]), int(pair[1])
        rows.append(
            {
                "pair_number": pair_number,
                "patient_position_i": i,
                "patient_position_j": j,
                "patient_id_source": patient_id_source,
                "patient_id_i": patient_ids[i],
                "patient_id_j": patient_ids[j],
                "duration_i": float(work.at[i, duration_col]),
                "duration_j": float(work.at[j, duration_col]),
                "event_i": int(work.at[i, event_col]),
                "event_j": int(work.at[j, event_col]),
            }
        )
    return pd.DataFrame(rows)


def _normalize_reusable_outer_pairs(
    outer_pairs,
    work,
    original_index,
    duration_col,
    event_col,
    event_policy,
):
    """Validate and safely remap a caller-supplied fixed outer-pair plan.

    When patient identifiers are present they are authoritative. Stored row
    positions are treated only as a cache and may change between runs. This
    makes pair plans invariant to harmless DataFrame row reordering while still
    rejecting genuinely stale or incompatible plans.
    """

    plan_frame = outer_pairs if isinstance(outer_pairs, pd.DataFrame) else None
    matched_id_source = None
    if plan_frame is not None:
        required = {"patient_position_i", "patient_position_j"}
        missing = required.difference(plan_frame.columns)
        if missing:
            raise ValueError(
                "outer_pairs DataFrame is missing columns: {}".format(sorted(missing))
            )

        matched_id_source, id_mapped_pairs = _match_pair_plan_identifier_source(
            plan_frame, work, original_index
        )
        if id_mapped_pairs is not None:
            raw_pairs = id_mapped_pairs
            stored_pairs = list(
                zip(
                    plan_frame["patient_position_i"].tolist(),
                    plan_frame["patient_position_j"].tolist(),
                )
            )
            remapped_rows = 0
            for stored, mapped in zip(stored_pairs, raw_pairs):
                stored_pair = (int(stored[0]), int(stored[1]))
                mapped_pair = (int(mapped[0]), int(mapped[1]))
                if stored_pair != mapped_pair:
                    remapped_rows += 1
            if remapped_rows:
                warnings.warn(
                    "Reusable outer-pair plan was remapped by patient ID using {!r}: "
                    "{} of {} rows had different cached patient positions. This is "
                    "expected after harmless DataFrame row reordering.".format(
                        matched_id_source, remapped_rows, len(raw_pairs)
                    ),
                    RuntimeWarning,
                    stacklevel=3,
                )
        else:
            raw_pairs = list(
                zip(
                    plan_frame["patient_position_i"].tolist(),
                    plan_frame["patient_position_j"].tolist(),
                )
            )
    else:
        raw_pairs = list(outer_pairs)

    n_patients = len(work)
    event_values = work[event_col].to_numpy(dtype=int)
    duration_values = work[duration_col].to_numpy(dtype=float)
    normalized = []
    seen = set()
    for row_index, value in enumerate(raw_pairs):
        if not isinstance(value, (tuple, list, np.ndarray)) or len(value) != 2:
            raise ValueError("Each outer pair must contain exactly two patient positions.")

        # Keep endpoint orientation for row-level event/duration validation.
        raw_i, raw_j = int(value[0]), int(value[1])
        if raw_i == raw_j:
            raise ValueError("An outer pair cannot contain the same patient twice.")
        if raw_i < 0 or raw_i >= n_patients or raw_j < 0 or raw_j >= n_patients:
            raise IndexError(
                "Outer pair ({}, {}) is outside the patient range.".format(raw_i, raw_j)
            )

        i, j = sorted((raw_i, raw_j))
        pair = (i, j)
        if pair in seen:
            raise ValueError("Duplicate outer pair {} in reusable plan.".format(pair))
        if not _pair_is_eligible_under_policy(
            i,
            j,
            event_values,
            duration_values,
            event_policy,
        ):
            raise ValueError(
                "Outer pair {} is not eligible under event_policy={!r}.".format(
                    pair, event_policy
                )
            )

        if plan_frame is not None:
            # Access each column directly rather than through a heterogeneous
            # pandas Series row. This prevents integer IDs from becoming 0.0,
            # 1.0, ... solely due to dtype promotion.
            if "event_i" in plan_frame.columns and int(plan_frame["event_i"].iloc[row_index]) != int(event_values[raw_i]):
                raise ValueError("Reusable outer-pair plan event_i does not match current data after patient-ID mapping.")
            if "event_j" in plan_frame.columns and int(plan_frame["event_j"].iloc[row_index]) != int(event_values[raw_j]):
                raise ValueError("Reusable outer-pair plan event_j does not match current data after patient-ID mapping.")
            if "duration_i" in plan_frame.columns and not np.isclose(float(plan_frame["duration_i"].iloc[row_index]), float(duration_values[raw_i])):
                raise ValueError("Reusable outer-pair plan duration_i does not match current data after patient-ID mapping.")
            if "duration_j" in plan_frame.columns and not np.isclose(float(plan_frame["duration_j"].iloc[row_index]), float(duration_values[raw_j])):
                raise ValueError("Reusable outer-pair plan duration_j does not match current data after patient-ID mapping.")

        seen.add(pair)
        normalized.append(pair)

    if not normalized:
        raise ValueError("outer_pairs must contain at least one pair.")

    return normalized

def make_outer_pair_plan(
    data,
    duration_col="Response_duration",
    event_col="Response_coded",
    max_outer_pairs=120,
    event_policy="all",
    random_state=42,
):
    """Create a deterministic reusable outer-pair plan as a DataFrame.

    Pass the returned DataFrame as ``outer_pairs=...`` to multiple calls of
    :func:`nested_pair_cox_evaluation` to guarantee identical outer test pairs.
    """

    work, _, original_index = _validate_and_prepare_data(
        data=data,
        duration_col=duration_col,
        event_col=event_col,
        exclude_covariates=(),
        stabilize_categorical_levels=True,
    )
    event_policy = _validate_event_policy(event_policy)
    pair_seed = np.random.SeedSequence(random_state).spawn(1)[0]
    rng = np.random.default_rng(pair_seed)
    pairs = _choose_outer_pairs(
        n_patients=len(work),
        max_pairs=max_outer_pairs,
        rng=rng,
        event_values=work[event_col].to_numpy(dtype=int),
        duration_values=work[duration_col].to_numpy(dtype=float),
        event_policy=event_policy,
    )
    return _outer_pair_plan_from_pairs(
        pairs,
        work,
        original_index,
        duration_col,
        event_col,
    )


def save_outer_pair_plan(pair_plan, output_path):
    """Save a reusable outer-pair plan to CSV."""
    if not isinstance(pair_plan, pd.DataFrame):
        raise TypeError("pair_plan must be a pandas DataFrame.")
    output_path = os.path.abspath(str(output_path))
    parent = os.path.dirname(output_path)
    if parent and not os.path.isdir(parent):
        os.makedirs(parent)
    pair_plan.to_csv(output_path, index=False)
    return output_path


def load_outer_pair_plan(input_path):
    """Load a reusable outer-pair plan previously saved as CSV.

    Patient identifiers are read as text so a CSV round-trip cannot silently
    change integer-like IDs (for example 0 -> 0.0). Legacy plans without a
    ``patient_id_source`` column remain supported.
    """
    input_path = os.path.abspath(str(input_path))
    if not os.path.isfile(input_path):
        raise IOError("Outer-pair plan does not exist: {}".format(input_path))
    dtype = {
        "patient_id_i": str,
        "patient_id_j": str,
        "patient_id_source": str,
    }
    return pd.read_csv(input_path, dtype=dtype)


def _choose_outer_pairs(
    n_patients,
    max_pairs,
    rng,
    event_values=None,
    duration_values=None,
    event_policy="all",
):
    """Choose unique outer pairs according to an event-selection policy.

    Policies
    --------
    all
        Generate pairs from all patients.
    both_events
        Generate pairs only when both patients have event == 1.
    valid_pairs
        Generate pairs with at least one patient having event == 1.
    harrell_comparable
        Generate only pairs that are comparable under Harrell's C-index rule,
        based solely on observed durations and event indicators.

    For n_patients <= 16 every policy-eligible pair is returned. For larger
    cohorts, a uniform unique sample of at most max_pairs eligible pairs is
    returned. The usual scoring-time comparability check is retained as a
    consistency safeguard.
    """

    if int(max_pairs) != max_pairs or not 1 <= int(max_pairs) <= 120:
        raise ValueError("max_pairs must be an integer from 1 through 120.")
    max_pairs = int(max_pairs)
    n_patients = int(n_patients)
    event_policy = _validate_event_policy(event_policy)

    if event_policy == "all":
        values = None
    else:
        values = _validate_event_values(event_values, n_patients)

    if event_policy == "harrell_comparable":
        durations = _validate_duration_values(duration_values, n_patients)
    else:
        durations = None

    total_eligible = _count_policy_eligible_pairs(
        n_patients=n_patients,
        event_values=values,
        event_policy=event_policy,
        duration_values=durations,
    )
    if total_eligible < 1:
        raise ValueError(
            "No outer patient pairs are eligible under event_policy={!r}.".format(
                event_policy
            )
        )

    if n_patients <= 16 or total_eligible <= max_pairs:
        return _all_policy_eligible_pairs(
            n_patients=n_patients,
            event_values=values,
            event_policy=event_policy,
            duration_values=durations,
        )

    target = min(max_pairs, total_eligible)
    sampled = set()

    if event_policy == "all":
        while len(sampled) < target:
            first, second = rng.choice(n_patients, size=2, replace=False)
            i, j = sorted((int(first), int(second)))
            sampled.add((i, j))
        return sorted(sampled)

    event_indices = np.flatnonzero(values == 1).astype(int)

    if event_policy == "both_events":
        if len(event_indices) < 2:
            raise ValueError(
                "At least two observed-event patients are required for "
                "event_policy='both_events'."
            )
        while len(sampled) < target:
            selected = rng.choice(event_indices, size=2, replace=False)
            i, j = sorted((int(selected[0]), int(selected[1])))
            sampled.add((i, j))
        return sorted(sampled)

    if event_policy == "harrell_comparable":
        # Reservoir sampling gives an exactly uniform sample of up to max_pairs
        # from the stream of all Harrell-comparable unordered pairs without
        # storing the complete eligible-pair set. This is preferable to rejection
        # sampling when comparable pairs are rare.
        reservoir = []
        seen = 0
        for i, j in combinations(range(n_patients), 2):
            if not _is_harrell_comparable_pair(
                float(durations[i]),
                int(values[i]),
                float(durations[j]),
                int(values[j]),
            ):
                continue
            seen += 1
            pair = (int(i), int(j))
            if len(reservoir) < target:
                reservoir.append(pair)
            else:
                replacement = int(rng.integers(0, seen))
                if replacement < target:
                    reservoir[replacement] = pair
        return sorted(reservoir)

    # valid_pairs. Sample uniformly over all eligible unordered pairs by first
    # choosing between event-event and event-censored categories in proportion
    # to the number of pairs in each category, then sampling uniformly within
    # that category.
    censored_indices = np.flatnonzero(values == 0).astype(int)
    n_event_event = _number_of_pairs(len(event_indices))
    n_event_censored = len(event_indices) * len(censored_indices)
    valid_total = n_event_event + n_event_censored

    while len(sampled) < target:
        choose_event_event = (
            n_event_event > 0
            and (
                n_event_censored == 0
                or float(rng.random()) < (n_event_event / float(valid_total))
            )
        )
        if choose_event_event:
            selected = rng.choice(event_indices, size=2, replace=False)
            i, j = sorted((int(selected[0]), int(selected[1])))
        else:
            event_index = int(rng.choice(event_indices))
            censored_index = int(rng.choice(censored_indices))
            i, j = sorted((event_index, censored_index))
        sampled.add((i, j))

    return sorted(sampled)


def _pair_contribution(
    duration_i: float,
    event_i: int,
    risk_i: float,
    duration_j: float,
    event_j: int,
    risk_j: float,
    risk_tie_tolerance: float,
) -> Tuple[bool, float]:
    """Return whether a pair is comparable and its Harrell-C contribution."""

    if not _is_harrell_comparable_pair(
        duration_i, event_i, duration_j, event_j
    ):
        return False, np.nan

    if duration_i < duration_j:
        earlier_risk, later_risk = risk_i, risk_j
    elif duration_j < duration_i:
        earlier_risk, later_risk = risk_j, risk_i
    else:
        # Equal observed times are comparable here only for event/censor pairs.
        if event_i == 1:
            earlier_risk, later_risk = risk_i, risk_j
        else:
            earlier_risk, later_risk = risk_j, risk_i

    difference = float(earlier_risk - later_risk)
    tolerance = risk_tie_tolerance * max(
        1.0, abs(float(earlier_risk)), abs(float(later_risk))
    )
    if abs(difference) <= tolerance:
        return True, 0.5
    return True, 1.0 if difference > 0 else 0.0


def _concordance_counts(
    durations: np.ndarray,
    events: np.ndarray,
    risks: np.ndarray,
    risk_tie_tolerance: float,
) -> Tuple[float, int]:
    contribution_sum = 0.0
    comparable = 0
    for i, j in combinations(range(len(durations)), 2):
        is_comparable, contribution = _pair_contribution(
            float(durations[i]),
            int(events[i]),
            float(risks[i]),
            float(durations[j]),
            int(events[j]),
            float(risks[j]),
            risk_tie_tolerance,
        )
        if is_comparable:
            comparable += 1
            contribution_sum += contribution
    return contribution_sum, comparable


def _make_inner_splits(
    train_data: pd.DataFrame,
    event_col: str,
    requested_splits: int,
    random_state: int,
) -> List[Tuple[np.ndarray, np.ndarray]]:
    n = len(train_data)
    if requested_splits < 2:
        raise ValueError("inner_splits must be at least 2.")

    # Keep validation folds at roughly two or more patients where possible.
    n_splits = min(int(requested_splits), max(2, n // 2))
    n_splits = min(n_splits, n)
    if n_splits < 2:
        raise ValueError("Not enough outer-training patients for inner CV.")

    y = train_data[event_col].to_numpy(dtype=int)
    class_counts = np.bincount(y, minlength=2)
    can_stratify = class_counts.min() >= 2

    if can_stratify:
        n_splits = min(n_splits, int(class_counts.min()))
        splitter = StratifiedKFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=random_state,
        )
        return list(splitter.split(np.zeros(n), y))

    splitter = KFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=random_state,
    )
    return list(splitter.split(np.zeros(n)))



def _quantile_breakpoints(train_durations, n_breaks=3):
    """Deterministically derive piecewise baseline breakpoints from training data.

    Interior duration quantiles are used, so held-out/validation durations never
    influence breakpoint placement. Duplicate quantiles are collapsed. Values
    equal to the observed minimum or maximum are removed because lifelines
    piecewise breakpoints should define genuine interior changes in baseline
    hazard.
    """

    if isinstance(n_breaks, bool) or int(n_breaks) != n_breaks or int(n_breaks) < 1:
        raise ValueError("n_breaks must be an integer >= 1.")
    n_breaks = int(n_breaks)

    durations = np.asarray(train_durations, dtype=float).reshape(-1)
    durations = durations[np.isfinite(durations)]
    if len(durations) < 2:
        raise ValueError(
            "At least two finite training durations are required for piecewise breakpoints."
        )

    minimum = float(np.min(durations))
    maximum = float(np.max(durations))
    if not maximum > minimum:
        raise ValueError(
            "Piecewise baseline estimation requires at least two distinct training durations."
        )

    qs = np.linspace(0.0, 1.0, n_breaks + 2)[1:-1]
    breakpoints = np.quantile(durations, qs).tolist()
    breakpoints = [
        float(value)
        for value in breakpoints
        if np.isfinite(value) and minimum < float(value) < maximum
    ]
    breakpoints = sorted(set(breakpoints))

    if not breakpoints:
        raise ValueError(
            "Quantile breakpoint calculation produced no distinct interior breakpoints."
        )
    return breakpoints

def _formula_categorical_columns(formula):
    """Return simple column names explicitly wrapped by Formulaic ``C(...)``."""
    if formula is None:
        return set()
    pattern = r"C\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)"
    return set(re.findall(pattern, str(formula)))


def _resolve_standardize_columns(
    work,
    covariates,
    standardize_predictors,
    standardize_columns,
):
    """Resolve numeric predictors eligible for training-fold standardization."""
    if not standardize_predictors:
        return []

    if standardize_columns is None:
        candidates = list(covariates)
    else:
        candidates = list(dict.fromkeys(standardize_columns))
        unknown = sorted(set(candidates).difference(covariates))
        if unknown:
            raise KeyError("Unknown standardize_columns: {}".format(unknown))

    resolved = []
    for column in candidates:
        if not pd.api.types.is_numeric_dtype(work[column].dtype):
            if standardize_columns is not None:
                raise TypeError(
                    "standardize_columns contains non-numeric predictor {!r}.".format(column)
                )
            continue
        resolved.append(column)
    return resolved


def _standardize_training_frame(train_data, spec, candidate_columns):
    """Fit z-score parameters on training data and return transformed training data."""
    transformed = train_data.copy()
    stats = {}
    categorical_columns = _formula_categorical_columns(spec.formula)
    for column in candidate_columns:
        if column in categorical_columns:
            continue
        values = pd.to_numeric(transformed[column], errors="raise").to_numpy(dtype=float)
        mean = float(np.mean(values))
        scale = float(np.std(values, ddof=0))
        if not np.isfinite(mean) or not np.isfinite(scale):
            raise ValueError("Non-finite standardization statistic for {!r}.".format(column))
        if scale <= 1.0e-12:
            continue
        transformed[column] = (values - mean) / scale
        stats[column] = (mean, scale)
    return transformed, stats


def _apply_model_standardization(model, data):
    transformed = data.copy()
    stats = getattr(model, "_nested_standardization_stats", {})
    for column, values in stats.items():
        mean, scale = values
        transformed[column] = (
            pd.to_numeric(transformed[column], errors="raise").astype(float) - mean
        ) / scale
    return transformed


def transform_for_fitted_cox_model(model, data):
    """Apply training-derived standardization before direct final-model prediction."""
    return _apply_model_standardization(model, data)


def _fit_cox_model(
    train_data: pd.DataFrame,
    duration_col: str,
    event_col: str,
    spec: CoxModelSpec,
    covariates=None,
    standardize_predictors=False,
    standardize_columns=None,
) -> Tuple[Any, List[str]]:
    _require_lifelines()

    fit_data = train_data
    stats = {}
    if standardize_predictors:
        fit_data, stats = _standardize_training_frame(
            train_data,
            spec=spec,
            candidate_columns=list(standardize_columns or []),
        )

    kwargs: Dict[str, Any] = {
        "baseline_estimation_method": spec.base_model,
        "penalizer": spec.penalizer,
        "l1_ratio": spec.l1_ratio,
    }
    if spec.base_model == "spline":
        kwargs["n_baseline_knots"] = spec.knot_number
    elif spec.base_model == "piecewise":
        kwargs["breakpoints"] = _quantile_breakpoints(
            fit_data[duration_col].to_numpy(dtype=float),
            n_breaks=spec.piecewise_break_number,
        )

    model = CoxPHFitter(**kwargs)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        model.fit(
            fit_data,
            duration_col=duration_col,
            event_col=event_col,
            formula=spec.formula,
            show_progress=False,
        )
    model._nested_standardization_stats = stats
    model._nested_standardized_columns = tuple(sorted(stats))
    warning_messages = [str(item.message) for item in captured]
    return model, warning_messages


def _reference_prediction_time(
    train_data: pd.DataFrame,
    duration_col: str,
    event_col: str,
) -> float:
    event_times = train_data.loc[train_data[event_col] == 1, duration_col]
    if event_times.empty:
        raise ValueError("The training set has no observed event.")
    return float(event_times.max())


def _predict_risk(
    model: Any,
    train_data: pd.DataFrame,
    test_data: pd.DataFrame,
    covariates: Sequence[str],
    duration_col: str,
    event_col: str,
) -> Tuple[np.ndarray, float]:
    """Predict cumulative hazard at a common time; larger means higher risk."""

    prediction_time = _reference_prediction_time(train_data, duration_col, event_col)
    transformed_test = _apply_model_standardization(model, test_data)
    prediction = model.predict_cumulative_hazard(
        transformed_test[list(covariates)],
        times=[prediction_time],
    )
    values = np.asarray(prediction, dtype=float)
    if values.ndim == 1:
        risks = values
    elif values.shape[0] == 1:
        risks = values[0, :]
    elif values.shape[1] == 1 and len(test_data) == values.shape[0]:
        risks = values[:, 0]
    else:
        raise RuntimeError(
            "Unexpected shape from predict_cumulative_hazard: {}.".format(values.shape)
        )
    risks = np.asarray(risks, dtype=float).reshape(-1)
    if len(risks) != len(test_data):
        raise RuntimeError(
            "Risk prediction length does not match the number of test patients."
        )
    if not np.isfinite(risks).all():
        raise FloatingPointError("The fitted model produced non-finite risks.")
    return risks, prediction_time


def _predict_log_relative_risk(
    model,
    train_data,
    evaluation_data,
    covariates,
    duration_col,
    event_col,
):
    """Return a log relative-risk score suitable for Cox partial likelihood.

    ``predict_partial_hazard`` is preferred when the fitted lifelines backend
    exposes it. For parametric Cox baseline backends where that method is not
    available, cumulative hazard at one common time is used instead. Under a
    proportional-hazards model this differs from ``exp(x beta)`` only by a
    positive baseline factor common to all evaluated patients, which cancels
    exactly from the Cox partial likelihood.
    """
    transformed = _apply_model_standardization(model, evaluation_data)
    try:
        prediction = model.predict_partial_hazard(transformed[list(covariates)])
        values = np.asarray(prediction, dtype=float).reshape(-1)
        if len(values) != len(evaluation_data):
            raise RuntimeError("Unexpected predict_partial_hazard length.")
        if (values <= 0).any() or not np.isfinite(values).all():
            raise FloatingPointError("Non-positive/non-finite partial hazards.")
        return np.log(values)
    except (AttributeError, NotImplementedError):
        pass

    risks, _ = _predict_risk(
        model,
        train_data=train_data,
        test_data=evaluation_data,
        covariates=covariates,
        duration_col=duration_col,
        event_col=event_col,
    )
    if (risks <= 0).any():
        raise FloatingPointError(
            "Cumulative-hazard fallback produced a non-positive relative risk."
        )
    return np.log(risks)



def _fit_prediction_ensemble(
    ranking,
    train_data,
    evaluation_data,
    covariates,
    duration_col,
    event_col,
    standardize_predictors,
    standardize_columns,
    ensemble_size,
    ensemble_weighting,
    ensemble_temperature,
    ensemble_distinct_formulas,
    ensemble_rank_level,
):
    """Fit top inner-ranked models and average standardized log-risk scores.

    Selection and all score normalization use only ``train_data``. Candidate
    models that fail on the full training set are skipped and lower-ranked
    eligible candidates are tried until the requested ensemble size is reached.
    """
    fitted = []
    errors = []
    seen_formulas = set()
    eligible = ranking.loc[ranking["eligible"]]

    for _, candidate_row in eligible.iterrows():
        spec = candidate_row["spec"]
        formula_key = spec.formula_label
        if ensemble_distinct_formulas and formula_key in seen_formulas:
            continue
        try:
            model, fit_warnings = _fit_cox_model(
                train_data,
                duration_col=duration_col,
                event_col=event_col,
                spec=spec,
                covariates=covariates,
                standardize_predictors=standardize_predictors,
                standardize_columns=standardize_columns,
            )
            train_scores = _predict_log_relative_risk(
                model,
                train_data=train_data,
                evaluation_data=train_data,
                covariates=covariates,
                duration_col=duration_col,
                event_col=event_col,
            )
            evaluation_scores = _predict_log_relative_risk(
                model,
                train_data=train_data,
                evaluation_data=evaluation_data,
                covariates=covariates,
                duration_col=duration_col,
                event_col=event_col,
            )
            train_scores = np.asarray(train_scores, dtype=float).reshape(-1)
            evaluation_scores = np.asarray(evaluation_scores, dtype=float).reshape(-1)
            center = float(np.mean(train_scores))
            scale = float(np.std(train_scores, ddof=0))
            if not np.isfinite(center) or not np.isfinite(scale) or scale <= 1e-12:
                raise FloatingPointError(
                    "Component log-risk score has zero/non-finite training variance."
                )
            z_scores = (evaluation_scores - center) / scale
            if not np.isfinite(z_scores).all():
                raise FloatingPointError("Component produced non-finite standardized scores.")
            fitted.append(
                {
                    "candidate_row": candidate_row,
                    "spec": spec,
                    "model": model,
                    "warnings": fit_warnings,
                    "center": center,
                    "scale": scale,
                    "evaluation_scores": z_scores,
                }
            )
            seen_formulas.add(formula_key)
            if len(fitted) >= int(ensemble_size):
                break
        except Exception as exc:
            errors.append("{}: {}: {}".format(spec.label, type(exc).__name__, exc))

    if not fitted:
        raise RuntimeError(
            "No eligible candidate could be fitted as an ensemble component. "
            + " | ".join(errors[:3])
        )

    candidate_rows = [item["candidate_row"] for item in fitted]
    weights, weight_ranks = _ensemble_weights(
        candidate_rows, ensemble_weighting, ensemble_temperature, ensemble_rank_level
    )
    for item, weight_rank in zip(fitted, weight_ranks):
        item["ensemble_weight_rank"] = float(weight_rank)
    matrix = np.vstack([item["evaluation_scores"] for item in fitted])
    ensemble_scores = np.dot(weights, matrix)
    if not np.isfinite(ensemble_scores).all():
        raise FloatingPointError("Ensemble produced non-finite risk scores.")

    ensemble_model = CoxPredictionEnsemble(
        models=[item["model"] for item in fitted],
        specs=[item["spec"] for item in fitted],
        weights=weights,
        score_centers=[item["center"] for item in fitted],
        score_scales=[item["scale"] for item in fitted],
        covariates=covariates,
        duration_col=duration_col,
        event_col=event_col,
        training_data=train_data,
        weighting=ensemble_weighting,
        rank_level=ensemble_rank_level,
    )
    return ensemble_model, fitted, np.asarray(weights, dtype=float), ensemble_scores, errors


def _efron_partial_log_likelihood(durations, events, log_relative_risks):
    """Unpenalized Cox partial log-likelihood with Efron handling of ties.

    This scoring implementation is intentionally independent of the baseline
    hazard parameterization. It evaluates only the relative-risk component of
    the fitted Cox model. Lifelines' traditional Breslow-baseline Cox fitter
    uses Efron's approximation for tied event times, so Efron is used here too.
    """
    durations = np.asarray(durations, dtype=float).reshape(-1)
    events = np.asarray(events, dtype=int).reshape(-1)
    eta = np.asarray(log_relative_risks, dtype=float).reshape(-1)
    if not (len(durations) == len(events) == len(eta)):
        raise ValueError("durations, events, and log_relative_risks must align.")
    if len(durations) == 0:
        raise ValueError("Partial likelihood requires at least one patient.")
    if not np.isfinite(durations).all() or not np.isfinite(eta).all():
        raise ValueError("Partial-likelihood inputs must be finite.")
    if not set(np.unique(events)).issubset({0, 1}):
        raise ValueError("events must contain only 0/1 values.")
    if int(events.sum()) == 0:
        return 0.0

    # A global shift prevents overflow; it cancels algebraically from every
    # numerator/denominator contribution to the partial log-likelihood.
    shift = float(np.max(eta))
    eta_shifted = eta - shift
    relative = np.exp(eta_shifted)
    total = 0.0

    event_times = np.unique(durations[events == 1])
    for event_time in event_times:
        deaths = (durations == event_time) & (events == 1)
        risk_set = durations >= event_time
        d = int(deaths.sum())
        if d <= 0:
            continue
        risk_sum = float(relative[risk_set].sum())
        death_sum = float(relative[deaths].sum())
        total += float(eta_shifted[deaths].sum())
        for tied_index in range(d):
            denominator = risk_sum - (float(tied_index) / float(d)) * death_sum
            if denominator <= 0 or not np.isfinite(denominator):
                raise FloatingPointError(
                    "Invalid Efron partial-likelihood denominator."
                )
            total -= float(np.log(denominator))
    return float(total)


def _verweij_van_houwelingen_fold_score(
    model,
    full_data,
    training_positions,
    covariates,
    duration_col,
    event_col,
):
    """Grouped cross-validated partial log-likelihood contribution.

    For coefficients fitted without the held-out fold, compute

        PLL(full data; beta_minus_fold) - PLL(training data; beta_minus_fold).

    This is the Verweij-van Houwelingen construction and avoids forming Cox
    risk sets from a tiny validation fold alone.
    """
    training_positions = np.asarray(training_positions, dtype=int)
    fit_data = full_data.iloc[training_positions]
    eta_full = _predict_log_relative_risk(
        model,
        train_data=fit_data,
        evaluation_data=full_data,
        covariates=covariates,
        duration_col=duration_col,
        event_col=event_col,
    )
    pll_full = _efron_partial_log_likelihood(
        full_data[duration_col].to_numpy(dtype=float),
        full_data[event_col].to_numpy(dtype=int),
        eta_full,
    )
    pll_training = _efron_partial_log_likelihood(
        fit_data[duration_col].to_numpy(dtype=float),
        fit_data[event_col].to_numpy(dtype=int),
        eta_full[training_positions],
    )
    return float(pll_full - pll_training)


def _inner_repeat_seeds(random_state, inner_repeats):
    """Return deterministic shuffle seeds for repeated inner CV.

    The first seed is exactly ``random_state`` so ``inner_repeats=1`` preserves
    the historical single-randomization behavior. Additional seeds are drawn
    deterministically from a local RandomState seeded by the same value.
    """

    if isinstance(inner_repeats, bool) or int(inner_repeats) != inner_repeats:
        raise TypeError("inner_repeats must be a positive integer.")
    inner_repeats = int(inner_repeats)
    if inner_repeats < 1:
        raise ValueError("inner_repeats must be >= 1.")

    base_seed = int(random_state)
    seeds = [base_seed]
    if inner_repeats == 1:
        return seeds

    local_rng = np.random.RandomState(base_seed)
    maximum = int(np.iinfo(np.int32).max)
    used = set(seeds)
    while len(seeds) < inner_repeats:
        candidate = int(local_rng.randint(0, maximum))
        if candidate not in used:
            used.add(candidate)
            seeds.append(candidate)
    return seeds


def _validate_inner_validation_method(value):
    allowed = ("kfold", "unique_pairs")
    if value not in allowed:
        raise ValueError("inner_validation_method must be one of {}.".format(allowed))
    return value


def _validate_inner_scoring(value):
    allowed = ("c_index", "cv_partial_likelihood")
    if value not in allowed:
        raise ValueError("inner_scoring must be one of {}.".format(allowed))
    return value


def _validate_inner_selection_rule(value):
    allowed = ("mean", "lower_confidence")
    if value not in allowed:
        raise ValueError("inner_selection_rule must be one of {}.".format(allowed))
    return value



def _validate_prediction_aggregation(value):
    allowed = ("winner", "stacked")
    if value not in allowed:
        raise ValueError("prediction_aggregation must be one of {}.".format(allowed))
    return value


def _validate_ensemble_weighting(value):
    allowed = ("equal", "rank", "softmax")
    if value not in allowed:
        raise ValueError("ensemble_weighting must be one of {}.".format(allowed))
    return value


def _validate_ensemble_rank_level(value):
    allowed = ("specification", "formula")
    if value not in allowed:
        raise ValueError("ensemble_rank_level must be one of {}.".format(allowed))
    return value


def _validate_ensemble_controls(size, temperature):
    if isinstance(size, bool) or int(size) != size or int(size) < 1:
        raise ValueError("ensemble_size must be a positive integer.")
    temperature = float(temperature)
    if not np.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("ensemble_temperature must be a positive finite number.")
    return int(size), temperature


def _ensemble_weight_ranks(candidate_rows, rank_level):
    """Return the ranks used specifically for ensemble rank weighting.

    ``specification`` uses the raw inner rank of each hyperparameter/baseline
    specification. ``formula`` performs a dense rank over distinct formulas in
    the order of their best successfully fitted specification. This prevents
    a formula whose best variant appears at raw rank 26 from being treated as
    twenty-five times worse merely because ranks 2--25 are alternate variants
    of the first formula.
    """
    if rank_level == "specification":
        return np.asarray(
            [max(1.0, float(row["inner_rank"])) for row in candidate_rows],
            dtype=float,
        )

    formula_to_rank = {}
    next_rank = 1
    ranks = []
    for row in candidate_rows:
        spec = row["spec"]
        key = spec.formula_label
        if key not in formula_to_rank:
            formula_to_rank[key] = next_rank
            next_rank += 1
        ranks.append(float(formula_to_rank[key]))
    return np.asarray(ranks, dtype=float)


def _ensemble_weights(candidate_rows, weighting, temperature, rank_level):
    n = len(candidate_rows)
    if n < 1:
        raise ValueError("At least one fitted ensemble component is required.")
    weight_ranks = _ensemble_weight_ranks(candidate_rows, rank_level)
    if weighting == "equal" or n == 1:
        return np.repeat(1.0 / float(n), n), weight_ranks
    if weighting == "rank":
        raw = 1.0 / weight_ranks
        return raw / raw.sum(), weight_ranks

    values = np.asarray(
        [float(row["selection_score"]) for row in candidate_rows], dtype=float
    )
    if not np.isfinite(values).all():
        return np.repeat(1.0 / float(n), n), weight_ranks
    scale = float(np.std(values, ddof=0))
    if not np.isfinite(scale) or scale <= 1e-12:
        return np.repeat(1.0 / float(n), n), weight_ranks
    z = (values - float(np.mean(values))) / scale
    z = z / float(temperature)
    z = z - float(np.max(z))
    raw = np.exp(z)
    if not np.isfinite(raw).all() or float(raw.sum()) <= 0.0:
        return np.repeat(1.0 / float(n), n), weight_ranks
    return raw / raw.sum(), weight_ranks


def _selection_score(mean_value, standard_error, rule, penalty):
    if not np.isfinite(mean_value):
        return np.nan
    if rule == "mean":
        return float(mean_value)
    if not np.isfinite(standard_error):
        return np.nan
    return float(mean_value) - float(penalty) * float(standard_error)


def _rank_candidate_rows(rows):
    ranking = pd.DataFrame(rows)
    ranking = ranking.sort_values(
        ["eligible", "selection_score", "inner_primary_score", "inner_primary_score_std", "candidate_order"],
        ascending=[False, False, False, True, True],
        na_position="last",
        kind="mergesort",
    ).reset_index(drop=True)
    ranking["inner_rank"] = np.arange(1, len(ranking) + 1)
    return ranking


def _rank_candidates_repeated_kfold(
    train_data,
    covariates,
    duration_col,
    event_col,
    specs,
    inner_splits,
    inner_repeats,
    random_state,
    risk_tie_tolerance,
    inner_scoring,
    inner_selection_rule,
    inner_stability_penalty,
    standardize_predictors,
    standardize_columns,
):
    repeat_seeds = _inner_repeat_seeds(random_state, inner_repeats)
    repeated_splits = []
    for repeat_index, repeat_seed in enumerate(repeat_seeds, start=1):
        splits = _make_inner_splits(
            train_data,
            event_col=event_col,
            requested_splits=inner_splits,
            random_state=repeat_seed,
        )
        repeated_splits.append((repeat_index, repeat_seed, splits))

    rows = []
    for spec in specs:
        total_comparable = 0
        failed_folds = 0
        failed_repeats = 0
        warning_count = 0
        first_error = ""
        repeat_c_indices = []
        repeat_comparable_counts = []
        repeat_cvpl_values = []

        for repeat_index, repeat_seed, splits in repeated_splits:
            repeat_contribution = 0.0
            repeat_comparable = 0
            repeat_cvpl = 0.0
            repeat_failed_folds = 0
            for fit_positions, validation_positions in splits:
                fit_data = train_data.iloc[fit_positions]
                validation_data = train_data.iloc[validation_positions]
                try:
                    if int(fit_data[event_col].sum()) == 0:
                        raise ValueError("An inner training fold contains no observed event.")
                    model, fit_warnings = _fit_cox_model(
                        fit_data,
                        duration_col=duration_col,
                        event_col=event_col,
                        spec=spec,
                        covariates=covariates,
                        standardize_predictors=standardize_predictors,
                        standardize_columns=standardize_columns,
                    )
                    warning_count += len(fit_warnings)

                    # Always retain concordance as a diagnostic, even when
                    # partial likelihood is the model-selection objective.
                    risks, _ = _predict_risk(
                        model,
                        train_data=fit_data,
                        test_data=validation_data,
                        covariates=covariates,
                        duration_col=duration_col,
                        event_col=event_col,
                    )
                    contribution, comparable = _concordance_counts(
                        validation_data[duration_col].to_numpy(dtype=float),
                        validation_data[event_col].to_numpy(dtype=int),
                        risks,
                        risk_tie_tolerance=risk_tie_tolerance,
                    )
                    repeat_contribution += contribution
                    repeat_comparable += comparable

                    repeat_cvpl += _verweij_van_houwelingen_fold_score(
                        model,
                        full_data=train_data,
                        training_positions=fit_positions,
                        covariates=covariates,
                        duration_col=duration_col,
                        event_col=event_col,
                    )
                except Exception as exc:
                    failed_folds += 1
                    repeat_failed_folds += 1
                    if not first_error:
                        first_error = "repeat {} seed {}: {}: {}".format(
                            repeat_index, repeat_seed, type(exc).__name__, exc
                        )

            if repeat_failed_folds > 0:
                failed_repeats += 1
            else:
                repeat_cvpl_values.append(float(repeat_cvpl))
                if repeat_comparable > 0:
                    repeat_c_indices.append(
                        float(repeat_contribution) / float(repeat_comparable)
                    )
                    repeat_comparable_counts.append(int(repeat_comparable))
            total_comparable += repeat_comparable

        if len(repeat_c_indices) > 0:
            c_values = np.asarray(repeat_c_indices, dtype=float)
            mean_c = float(np.mean(c_values))
            sd_c = float(np.std(c_values, ddof=1)) if len(c_values) > 1 else 0.0
            se_c = sd_c / np.sqrt(float(len(c_values)))
            min_c = float(np.min(c_values))
            max_c = float(np.max(c_values))
        else:
            mean_c = sd_c = se_c = min_c = max_c = np.nan

        if len(repeat_cvpl_values) > 0:
            cvpl_values = np.asarray(repeat_cvpl_values, dtype=float)
            mean_cvpl = float(np.mean(cvpl_values))
            sd_cvpl = float(np.std(cvpl_values, ddof=1)) if len(cvpl_values) > 1 else 0.0
            se_cvpl = sd_cvpl / np.sqrt(float(len(cvpl_values)))
            min_cvpl = float(np.min(cvpl_values))
            max_cvpl = float(np.max(cvpl_values))
        else:
            mean_cvpl = sd_cvpl = se_cvpl = min_cvpl = max_cvpl = np.nan

        if inner_scoring == "c_index":
            primary_mean, primary_sd, primary_se = mean_c, sd_c, se_c
            primary_count = len(repeat_c_indices)
        else:
            primary_mean, primary_sd, primary_se = mean_cvpl, sd_cvpl, se_cvpl
            primary_count = len(repeat_cvpl_values)
        eligible = failed_folds == 0 and np.isfinite(primary_mean)
        score = _selection_score(
            primary_mean,
            primary_se,
            inner_selection_rule,
            inner_stability_penalty,
        )
        rows.append(
            {
                "spec": spec,
                "candidate_order": spec.order,
                "candidate_label": spec.label,
                "formula": spec.formula_label,
                "base_model": spec.base_model,
                "knot_number": spec.knot_number,
                "piecewise_break_number": spec.piecewise_break_number,
                "penalizer": spec.penalizer,
                "l1_ratio": spec.l1_ratio,
                "inner_validation_method": "kfold",
                "inner_scoring": inner_scoring,
                "inner_primary_score": primary_mean,
                "inner_primary_score_std": primary_sd,
                "inner_primary_score_se": primary_se,
                "inner_c_index": mean_c,
                "inner_c_index_std": sd_c,
                "inner_c_index_se": se_c,
                "inner_cv_partial_log_likelihood": mean_cvpl,
                "inner_cv_partial_log_likelihood_std": sd_cvpl,
                "inner_cv_partial_log_likelihood_se": se_cvpl,
                "selection_score": score,
                "inner_c_index_min": min_c,
                "inner_c_index_max": max_c,
                "inner_cv_partial_log_likelihood_min": min_cvpl,
                "inner_cv_partial_log_likelihood_max": max_cvpl,
                "inner_repeats_requested": int(inner_repeats),
                "inner_repeats_scored": primary_count,
                "inner_repeat_c_indices": ",".join(
                    "{:.12g}".format(v) for v in repeat_c_indices
                ),
                "inner_repeat_cv_partial_log_likelihood": ",".join(
                    "{:.12g}".format(v) for v in repeat_cvpl_values
                ),
                "inner_repeat_comparable_pairs": ",".join(
                    str(v) for v in repeat_comparable_counts
                ),
                "inner_comparable_pairs": total_comparable,
                "inner_unique_pairs_requested": 0,
                "inner_unique_pairs_scored": 0,
                "inner_failed_repeats": failed_repeats,
                "inner_failed_folds": failed_folds,
                "inner_warning_count": warning_count,
                "eligible": eligible,
                "first_error": first_error,
            }
        )
    return _rank_candidate_rows(rows)


def _unique_harrell_inner_pairs(train_data, duration_col, event_col, max_pairs, random_state):
    if isinstance(max_pairs, bool) or int(max_pairs) != max_pairs or int(max_pairs) < 1:
        raise ValueError("inner_max_pairs must be a positive integer.")
    max_pairs = int(max_pairs)
    durations = train_data[duration_col].to_numpy(dtype=float)
    events = train_data[event_col].to_numpy(dtype=int)
    eligible = [
        (i, j)
        for i, j in combinations(range(len(train_data)), 2)
        if _is_harrell_comparable_pair(
            float(durations[i]), int(events[i]), float(durations[j]), int(events[j])
        )
    ]
    if not eligible:
        raise ValueError("No Harrell-comparable inner patient pairs are available.")
    if len(eligible) <= max_pairs:
        return eligible
    rng = np.random.default_rng(np.random.SeedSequence(int(random_state)))
    chosen = rng.choice(len(eligible), size=max_pairs, replace=False)
    return sorted([eligible[int(index)] for index in chosen])


def _rank_candidates_unique_pairs(
    train_data,
    covariates,
    duration_col,
    event_col,
    specs,
    inner_max_pairs,
    random_state,
    risk_tie_tolerance,
    inner_scoring,
    inner_selection_rule,
    inner_stability_penalty,
    standardize_predictors,
    standardize_columns,
):
    pairs = _unique_harrell_inner_pairs(
        train_data,
        duration_col=duration_col,
        event_col=event_col,
        max_pairs=inner_max_pairs,
        random_state=random_state,
    )
    n_patients = len(train_data)
    rows = []
    for spec in specs:
        contributions = []
        cvpl_values = []
        failed_pairs = 0
        warning_count = 0
        first_error = ""
        for i, j in pairs:
            excluded = {i, j}
            fit_positions = [p for p in range(n_patients) if p not in excluded]
            fit_data = train_data.iloc[fit_positions]
            validation_data = train_data.iloc[[i, j]]
            try:
                if int(fit_data[event_col].sum()) == 0:
                    raise ValueError("An inner pair training set contains no observed event.")
                model, fit_warnings = _fit_cox_model(
                    fit_data,
                    duration_col=duration_col,
                    event_col=event_col,
                    spec=spec,
                    covariates=covariates,
                    standardize_predictors=standardize_predictors,
                    standardize_columns=standardize_columns,
                )
                warning_count += len(fit_warnings)

                risks, _ = _predict_risk(
                    model,
                    train_data=fit_data,
                    test_data=validation_data,
                    covariates=covariates,
                    duration_col=duration_col,
                    event_col=event_col,
                )
                comparable, contribution = _pair_contribution(
                    float(validation_data.iloc[0][duration_col]),
                    int(validation_data.iloc[0][event_col]),
                    float(risks[0]),
                    float(validation_data.iloc[1][duration_col]),
                    int(validation_data.iloc[1][event_col]),
                    float(risks[1]),
                    risk_tie_tolerance,
                )
                if comparable:
                    contributions.append(float(contribution))

                cvpl_values.append(
                    _verweij_van_houwelingen_fold_score(
                        model,
                        full_data=train_data,
                        training_positions=fit_positions,
                        covariates=covariates,
                        duration_col=duration_col,
                        event_col=event_col,
                    )
                )
            except Exception as exc:
                failed_pairs += 1
                if not first_error:
                    first_error = "inner pair ({}, {}): {}: {}".format(
                        i, j, type(exc).__name__, exc
                    )

        if len(contributions) > 0:
            c_values = np.asarray(contributions, dtype=float)
            mean_c = float(np.mean(c_values))
            sd_c = float(np.std(c_values, ddof=1)) if len(c_values) > 1 else 0.0
            se_c = sd_c / np.sqrt(float(len(c_values)))
            min_c = float(np.min(c_values))
            max_c = float(np.max(c_values))
        else:
            mean_c = sd_c = se_c = min_c = max_c = np.nan

        if len(cvpl_values) > 0:
            pll_values = np.asarray(cvpl_values, dtype=float)
            mean_cvpl = float(np.mean(pll_values))
            sd_cvpl = float(np.std(pll_values, ddof=1)) if len(pll_values) > 1 else 0.0
            se_cvpl = sd_cvpl / np.sqrt(float(len(pll_values)))
            min_cvpl = float(np.min(pll_values))
            max_cvpl = float(np.max(pll_values))
        else:
            mean_cvpl = sd_cvpl = se_cvpl = min_cvpl = max_cvpl = np.nan

        if inner_scoring == "c_index":
            primary_mean, primary_sd, primary_se = mean_c, sd_c, se_c
            primary_count = len(contributions)
        else:
            primary_mean, primary_sd, primary_se = mean_cvpl, sd_cvpl, se_cvpl
            primary_count = len(cvpl_values)
        eligible = failed_pairs == 0 and np.isfinite(primary_mean)
        score = _selection_score(
            primary_mean,
            primary_se,
            inner_selection_rule,
            inner_stability_penalty,
        )
        rows.append(
            {
                "spec": spec,
                "candidate_order": spec.order,
                "candidate_label": spec.label,
                "formula": spec.formula_label,
                "base_model": spec.base_model,
                "knot_number": spec.knot_number,
                "piecewise_break_number": spec.piecewise_break_number,
                "penalizer": spec.penalizer,
                "l1_ratio": spec.l1_ratio,
                "inner_validation_method": "unique_pairs",
                "inner_scoring": inner_scoring,
                "inner_primary_score": primary_mean,
                "inner_primary_score_std": primary_sd,
                "inner_primary_score_se": primary_se,
                "inner_c_index": mean_c,
                "inner_c_index_std": sd_c,
                "inner_c_index_se": se_c,
                "inner_cv_partial_log_likelihood": mean_cvpl,
                "inner_cv_partial_log_likelihood_std": sd_cvpl,
                "inner_cv_partial_log_likelihood_se": se_cvpl,
                "selection_score": score,
                "inner_c_index_min": min_c,
                "inner_c_index_max": max_c,
                "inner_cv_partial_log_likelihood_min": min_cvpl,
                "inner_cv_partial_log_likelihood_max": max_cvpl,
                "inner_repeats_requested": 1,
                "inner_repeats_scored": 1 if eligible else 0,
                "inner_repeat_c_indices": "",
                "inner_repeat_cv_partial_log_likelihood": "",
                "inner_repeat_comparable_pairs": "",
                "inner_comparable_pairs": len(contributions),
                "inner_unique_pairs_requested": len(pairs),
                "inner_unique_pairs_scored": primary_count,
                "inner_failed_repeats": 0,
                "inner_failed_folds": failed_pairs,
                "inner_warning_count": warning_count,
                "eligible": eligible,
                "first_error": first_error,
            }
        )
    return _rank_candidate_rows(rows)


def _rank_candidates_inner_cv_base(
    train_data: pd.DataFrame,
    covariates: Sequence[str],
    duration_col: str,
    event_col: str,
    specs: Sequence[CoxModelSpec],
    inner_splits: int,
    inner_repeats: int,
    random_state: int,
    risk_tie_tolerance: float,
    inner_validation_method="kfold",
    inner_max_pairs=100,
    inner_scoring="c_index",
    inner_selection_rule="mean",
    inner_stability_penalty=1.0,
    standardize_predictors=False,
    standardize_columns=None,
) -> pd.DataFrame:
    """Rank candidates by C-index or grouped CV partial log-likelihood."""
    method = _validate_inner_validation_method(inner_validation_method)
    scoring = _validate_inner_scoring(inner_scoring)
    rule = _validate_inner_selection_rule(inner_selection_rule)
    if float(inner_stability_penalty) < 0:
        raise ValueError("inner_stability_penalty must be non-negative.")
    if method == "unique_pairs":
        return _rank_candidates_unique_pairs(
            train_data,
            covariates,
            duration_col,
            event_col,
            specs,
            inner_max_pairs,
            random_state,
            risk_tie_tolerance,
            scoring,
            rule,
            inner_stability_penalty,
            standardize_predictors,
            standardize_columns,
        )
    return _rank_candidates_repeated_kfold(
        train_data,
        covariates,
        duration_col,
        event_col,
        specs,
        inner_splits,
        inner_repeats,
        random_state,
        risk_tie_tolerance,
        scoring,
        rule,
        inner_stability_penalty,
        standardize_predictors,
        standardize_columns,
    )



def _validate_consensus_method(value):
    allowed = ("none", "top_models", "cox_filter")
    if value not in allowed:
        raise ValueError("consensus_method must be one of {}.".format(allowed))
    return value


def _validate_consensus_filter_scoring(value):
    """Validate the feature-group scoring rule used by cox_filter consensus."""
    allowed = ("partial_likelihood", "c_index")
    if value not in allowed:
        raise ValueError(
            "consensus_filter_scoring must be one of {}.".format(allowed)
        )
    return value


def _validate_consensus_controls(
    consensus_top_models,
    consensus_top_features,
    consensus_threshold,
    consensus_max_formula_groups,
):
    for name, value in (
        ("consensus_top_models", consensus_top_models),
        ("consensus_top_features", consensus_top_features),
        ("consensus_max_formula_groups", consensus_max_formula_groups),
    ):
        if isinstance(value, bool) or int(value) != value or int(value) < 1:
            raise ValueError("{} must be a positive integer.".format(name))
    threshold = float(consensus_threshold)
    if not (0.0 < threshold <= 1.0):
        raise ValueError("consensus_threshold must be in (0, 1].")
    return (
        int(consensus_top_models),
        int(consensus_top_features),
        threshold,
        int(consensus_max_formula_groups),
    )


def _validate_consensus_pair_batch_size(value):
    """Validate the number of unique inner pairs pooled per consensus resample."""
    if isinstance(value, bool) or int(value) != value:
        raise TypeError("consensus_pair_batch_size must be a positive integer.")
    value = int(value)
    if value < 1:
        raise ValueError("consensus_pair_batch_size must be >= 1.")
    return value


def _chunk_sequence(items, batch_size):
    """Return deterministic consecutive batches from an already randomized sequence."""
    size = _validate_consensus_pair_batch_size(batch_size)
    return [list(items[start:start + size]) for start in range(0, len(items), size)]


def _formula_feature_columns(formula, covariates):
    """Return covariate names referenced by a simple Formulaic formula.

    The evaluator's supported formula lists consist of ordinary column names,
    optional ``C(column)`` wrappers, and ``+`` combinations. Matching against
    the known covariate names is deliberately conservative and avoids trying to
    implement a second formula parser.
    """
    if formula is None:
        return set(covariates)
    text = str(formula)
    selected = set()
    for column in covariates:
        pattern = r"(?<![A-Za-z0-9_]){}(?![A-Za-z0-9_])".format(
            re.escape(str(column))
        )
        if re.search(pattern, text):
            selected.add(column)
    return selected


def _normalize_consensus_groups(covariates, feature_groups, group_expressions):
    """Validate grouped-consensus definitions and return deterministic maps.

    ``feature_groups`` is a mapping ``group_name -> sequence of data columns``.
    When it is omitted, every covariate is treated as its own group. When it is
    supplied, only the listed columns belong to the consensus universe. This is
    useful for a clinical-only analysis performed on a DataFrame that also
    contains PDR variables.
    """
    covariate_set = set(covariates)
    if feature_groups is None:
        names = [str(column) for column in covariates]
        members = {str(column): (column,) for column in covariates}
    else:
        if not hasattr(feature_groups, "items"):
            raise TypeError("consensus_feature_groups must be a mapping or None.")
        names = sorted([str(name) for name in feature_groups.keys()])
        members = {}
        used_columns = {}
        for name in names:
            raw = feature_groups[name]
            values = tuple(dict.fromkeys(list(raw)))
            if not values:
                raise ValueError(
                    "Consensus group {!r} contains no columns.".format(name)
                )
            unknown = sorted(set(values).difference(covariate_set))
            if unknown:
                raise KeyError(
                    "Consensus group {!r} contains unknown covariates: {}".format(
                        name, unknown
                    )
                )
            for column in values:
                if column in used_columns:
                    raise ValueError(
                        "Covariate {!r} occurs in consensus groups {!r} and {!r}. "
                        "Groups must be non-overlapping.".format(
                            column, used_columns[column], name
                        )
                    )
                used_columns[column] = name
            members[name] = values

    expressions = {}
    supplied = {} if group_expressions is None else dict(group_expressions)
    unknown_expression_groups = sorted(set(supplied).difference(names))
    if unknown_expression_groups:
        raise KeyError(
            "consensus_group_expressions contains unknown groups: {}".format(
                unknown_expression_groups
            )
        )
    for name in names:
        if name in supplied:
            expression = str(supplied[name]).strip()
        else:
            expression = " + ".join(str(column) for column in members[name])
        if not expression:
            raise ValueError(
                "Consensus group {!r} has an empty formula expression.".format(name)
            )
        expressions[name] = expression
    return names, members, expressions


def _groups_in_formula(formula, covariates, group_names, group_members):
    selected_columns = _formula_feature_columns(formula, covariates)
    result = []
    for name in group_names:
        if selected_columns.intersection(group_members[name]):
            result.append(name)
    return result


def _format_consensus_frequencies(group_names, frequencies):
    ordered = sorted(
        group_names,
        key=lambda name: (-float(frequencies.get(name, 0.0)), str(name)),
    )
    return ";".join(
        "{}={:.6g}".format(name, float(frequencies.get(name, 0.0)))
        for name in ordered
    )


def _select_consensus_groups(group_names, frequencies, threshold):
    selected = [
        name for name in group_names
        if float(frequencies.get(name, 0.0)) >= float(threshold)
    ]
    if selected:
        return sorted(
            selected,
            key=lambda name: (-float(frequencies.get(name, 0.0)), str(name)),
        )
    if not group_names:
        raise ValueError("No consensus feature groups are available.")
    # Small cohorts may not place any feature above a stringent threshold. A
    # deterministic one-group fallback keeps the outer fold estimable while
    # making the low stability visible in the reported frequencies.
    best = sorted(
        group_names,
        key=lambda name: (-float(frequencies.get(name, 0.0)), str(name)),
    )[0]
    return [best]


def _generate_consensus_formulas(
    selected_groups,
    frequencies,
    group_expressions,
    max_formula_groups,
):
    ordered = sorted(
        selected_groups,
        key=lambda name: (-float(frequencies.get(name, 0.0)), str(name)),
    )
    maximum = min(int(max_formula_groups), len(ordered))
    formulas = []
    seen = set()
    for size in range(1, maximum + 1):
        for combination in combinations(ordered, size):
            formula = " + ".join(group_expressions[name] for name in combination)
            if formula not in seen:
                seen.add(formula)
                formulas.append(formula)
    if not formulas:
        raise RuntimeError("Consensus formula generation produced no formulas.")
    return formulas


def _specs_for_generated_formulas(formulas, template_specs):
    """Apply the existing baseline/penalty grid to generated formulas."""
    configurations = []
    seen = set()
    for spec in template_specs:
        key = (
            spec.base_model,
            spec.knot_number,
            spec.piecewise_break_number,
            float(spec.penalizer),
            float(spec.l1_ratio),
        )
        if key not in seen:
            seen.add(key)
            configurations.append(key)
    result = []
    order = 0
    for config in configurations:
        base_model, knot_number, piecewise_break_number, penalizer, l1_ratio = config
        for formula in formulas:
            result.append(
                CoxModelSpec(
                    formula=formula,
                    base_model=base_model,
                    knot_number=knot_number,
                    piecewise_break_number=piecewise_break_number,
                    penalizer=penalizer,
                    l1_ratio=l1_ratio,
                    order=order,
                )
            )
            order += 1
    return result


def _rank_candidates_given_unique_pairs(
    train_data,
    covariates,
    duration_col,
    event_col,
    specs,
    pairs,
    risk_tie_tolerance,
    inner_scoring,
    inner_selection_rule,
    inner_stability_penalty,
    standardize_predictors,
    standardize_columns,
):
    """Rank specs on an explicitly supplied list of unique held-out pairs."""
    n_patients = len(train_data)
    rows = []
    for spec in specs:
        contributions = []
        cvpl_values = []
        failed_pairs = 0
        warning_count = 0
        first_error = ""
        for i, j in pairs:
            excluded = {int(i), int(j)}
            fit_positions = [p for p in range(n_patients) if p not in excluded]
            fit_data = train_data.iloc[fit_positions]
            validation_data = train_data.iloc[[int(i), int(j)]]
            try:
                if int(fit_data[event_col].sum()) == 0:
                    raise ValueError(
                        "An inner pair training set contains no observed event."
                    )
                model, fit_warnings = _fit_cox_model(
                    fit_data,
                    duration_col=duration_col,
                    event_col=event_col,
                    spec=spec,
                    covariates=covariates,
                    standardize_predictors=standardize_predictors,
                    standardize_columns=standardize_columns,
                )
                warning_count += len(fit_warnings)
                risks, _ = _predict_risk(
                    model,
                    train_data=fit_data,
                    test_data=validation_data,
                    covariates=covariates,
                    duration_col=duration_col,
                    event_col=event_col,
                )
                comparable, contribution = _pair_contribution(
                    float(validation_data.iloc[0][duration_col]),
                    int(validation_data.iloc[0][event_col]),
                    float(risks[0]),
                    float(validation_data.iloc[1][duration_col]),
                    int(validation_data.iloc[1][event_col]),
                    float(risks[1]),
                    risk_tie_tolerance,
                )
                if comparable:
                    contributions.append(float(contribution))
                cvpl_values.append(
                    _verweij_van_houwelingen_fold_score(
                        model,
                        full_data=train_data,
                        training_positions=fit_positions,
                        covariates=covariates,
                        duration_col=duration_col,
                        event_col=event_col,
                    )
                )
            except Exception as exc:
                failed_pairs += 1
                if not first_error:
                    first_error = "inner pair ({}, {}): {}: {}".format(
                        i, j, type(exc).__name__, exc
                    )

        if contributions:
            c_values = np.asarray(contributions, dtype=float)
            mean_c = float(np.mean(c_values))
            sd_c = float(np.std(c_values, ddof=1)) if len(c_values) > 1 else 0.0
            se_c = sd_c / np.sqrt(float(len(c_values)))
            min_c = float(np.min(c_values))
            max_c = float(np.max(c_values))
        else:
            mean_c = sd_c = se_c = min_c = max_c = np.nan

        if cvpl_values:
            pll_values = np.asarray(cvpl_values, dtype=float)
            mean_cvpl = float(np.mean(pll_values))
            sd_cvpl = (
                float(np.std(pll_values, ddof=1)) if len(pll_values) > 1 else 0.0
            )
            se_cvpl = sd_cvpl / np.sqrt(float(len(pll_values)))
            min_cvpl = float(np.min(pll_values))
            max_cvpl = float(np.max(pll_values))
        else:
            mean_cvpl = sd_cvpl = se_cvpl = min_cvpl = max_cvpl = np.nan

        if inner_scoring == "c_index":
            primary_mean, primary_sd, primary_se = mean_c, sd_c, se_c
            primary_count = len(contributions)
        else:
            primary_mean, primary_sd, primary_se = mean_cvpl, sd_cvpl, se_cvpl
            primary_count = len(cvpl_values)
        eligible = failed_pairs == 0 and np.isfinite(primary_mean)
        score = _selection_score(
            primary_mean,
            primary_se,
            inner_selection_rule,
            inner_stability_penalty,
        )
        rows.append(
            {
                "spec": spec,
                "candidate_order": spec.order,
                "candidate_label": spec.label,
                "formula": spec.formula_label,
                "base_model": spec.base_model,
                "knot_number": spec.knot_number,
                "piecewise_break_number": spec.piecewise_break_number,
                "penalizer": spec.penalizer,
                "l1_ratio": spec.l1_ratio,
                "inner_validation_method": "unique_pairs",
                "inner_scoring": inner_scoring,
                "inner_primary_score": primary_mean,
                "inner_primary_score_std": primary_sd,
                "inner_primary_score_se": primary_se,
                "inner_c_index": mean_c,
                "inner_c_index_std": sd_c,
                "inner_c_index_se": se_c,
                "inner_cv_partial_log_likelihood": mean_cvpl,
                "inner_cv_partial_log_likelihood_std": sd_cvpl,
                "inner_cv_partial_log_likelihood_se": se_cvpl,
                "selection_score": score,
                "inner_c_index_min": min_c,
                "inner_c_index_max": max_c,
                "inner_cv_partial_log_likelihood_min": min_cvpl,
                "inner_cv_partial_log_likelihood_max": max_cvpl,
                "inner_repeats_requested": 1,
                "inner_repeats_scored": 1 if eligible else 0,
                "inner_repeat_c_indices": "",
                "inner_repeat_cv_partial_log_likelihood": "",
                "inner_repeat_comparable_pairs": "",
                "inner_comparable_pairs": len(contributions),
                "inner_unique_pairs_requested": len(pairs),
                "inner_unique_pairs_scored": primary_count,
                "inner_failed_repeats": 0,
                "inner_failed_folds": failed_pairs,
                "inner_warning_count": warning_count,
                "eligible": eligible,
                "first_error": first_error,
            }
        )
    return _rank_candidate_rows(rows)



def _top_distinct_formula_rows(ranking, top_models):
    """Return the best eligible rows with distinct formula labels.

    Hyperparameter variants of one formula should not occupy multiple places in
    the top-model consensus set, otherwise a formula represented by many
    baseline/penalty combinations would receive artificial extra weight.
    """
    rows = []
    seen = set()
    for _, row in ranking.loc[ranking["eligible"]].iterrows():
        label = str(row["formula"])
        if label in seen:
            continue
        seen.add(label)
        rows.append(row)
        if len(rows) >= int(top_models):
            break
    return rows

def _consensus_from_top_models(
    train_data,
    covariates,
    duration_col,
    event_col,
    specs,
    inner_splits,
    inner_repeats,
    random_state,
    risk_tie_tolerance,
    inner_validation_method,
    inner_max_pairs,
    inner_scoring,
    inner_selection_rule,
    inner_stability_penalty,
    standardize_predictors,
    standardize_columns,
    group_names,
    group_members,
    top_models,
    pair_batch_size,
):
    counts = {name: 0 for name in group_names}
    retained_model_sets = 0
    resamples = 0

    if inner_validation_method == "kfold":
        seeds = _inner_repeat_seeds(random_state, inner_repeats)
        for seed in seeds:
            ranking = _rank_candidates_repeated_kfold(
                train_data,
                covariates,
                duration_col,
                event_col,
                specs,
                inner_splits,
                1,
                seed,
                risk_tie_tolerance,
                inner_scoring,
                inner_selection_rule,
                inner_stability_penalty,
                standardize_predictors,
                standardize_columns,
            )
            retained_rows = _top_distinct_formula_rows(ranking, top_models)
            if not retained_rows:
                continue
            resamples += 1
            for row in retained_rows:
                retained_model_sets += 1
                for name in _groups_in_formula(
                    row["spec"].formula,
                    covariates,
                    group_names,
                    group_members,
                ):
                    counts[name] += 1
    else:
        pairs = _unique_harrell_inner_pairs(
            train_data,
            duration_col=duration_col,
            event_col=event_col,
            max_pairs=inner_max_pairs,
            random_state=random_state,
        )
        for pair_batch in _chunk_sequence(pairs, pair_batch_size):
            ranking = _rank_candidates_given_unique_pairs(
                train_data,
                covariates,
                duration_col,
                event_col,
                specs,
                pair_batch,
                risk_tie_tolerance,
                inner_scoring,
                inner_selection_rule,
                inner_stability_penalty,
                standardize_predictors,
                standardize_columns,
            )
            retained_rows = _top_distinct_formula_rows(ranking, top_models)
            if not retained_rows:
                continue
            resamples += 1
            for row in retained_rows:
                retained_model_sets += 1
                for name in _groups_in_formula(
                    row["spec"].formula,
                    covariates,
                    group_names,
                    group_members,
                ):
                    counts[name] += 1

    if retained_model_sets <= 0:
        raise RuntimeError(
            "Consensus top-model stage produced no eligible retained models."
        )
    frequencies = {
        name: float(counts[name]) / float(retained_model_sets)
        for name in group_names
    }
    return {
        "frequencies": frequencies,
        "resamples": int(resamples),
        "source_items": int(retained_model_sets),
    }


def _cox_filter_group_improvement(
    fit_data,
    covariates,
    duration_col,
    event_col,
    expression,
    penalizer,
    standardize_predictors,
    standardize_columns,
):
    """Training-only Cox partial-likelihood improvement over the null model."""
    spec = CoxModelSpec(
        formula=expression,
        base_model="breslow",
        knot_number=None,
        piecewise_break_number=None,
        penalizer=float(penalizer),
        l1_ratio=0.0,
        order=0,
    )
    model, _ = _fit_cox_model(
        fit_data,
        duration_col=duration_col,
        event_col=event_col,
        spec=spec,
        covariates=covariates,
        standardize_predictors=standardize_predictors,
        standardize_columns=standardize_columns,
    )
    eta = _predict_log_relative_risk(
        model,
        train_data=fit_data,
        evaluation_data=fit_data,
        covariates=covariates,
        duration_col=duration_col,
        event_col=event_col,
    )
    durations = fit_data[duration_col].to_numpy(dtype=float)
    events = fit_data[event_col].to_numpy(dtype=int)
    model_pll = _efron_partial_log_likelihood(durations, events, eta)
    null_pll = _efron_partial_log_likelihood(
        durations, events, np.zeros(len(fit_data), dtype=float)
    )
    return float(model_pll - null_pll)


def _cox_filter_group_c_index_counts(
    fit_data,
    validation_data,
    covariates,
    duration_col,
    event_col,
    expression,
    penalizer,
    standardize_predictors,
    standardize_columns,
    risk_tie_tolerance,
):
    """Return held-out concordance contribution/count for one group-only fit."""
    spec = CoxModelSpec(
        formula=expression,
        base_model="breslow",
        knot_number=None,
        piecewise_break_number=None,
        penalizer=float(penalizer),
        l1_ratio=0.0,
        order=0,
    )
    model, _ = _fit_cox_model(
        fit_data,
        duration_col=duration_col,
        event_col=event_col,
        spec=spec,
        covariates=covariates,
        standardize_predictors=standardize_predictors,
        standardize_columns=standardize_columns,
    )
    risks = _predict_log_relative_risk(
        model,
        train_data=fit_data,
        evaluation_data=validation_data,
        covariates=covariates,
        duration_col=duration_col,
        event_col=event_col,
    )
    contribution, comparable = _concordance_counts(
        validation_data[duration_col].to_numpy(dtype=float),
        validation_data[event_col].to_numpy(dtype=int),
        np.asarray(risks, dtype=float).reshape(-1),
        risk_tie_tolerance,
    )
    return float(contribution), int(comparable)


def _consensus_from_cox_filter(
    train_data,
    covariates,
    duration_col,
    event_col,
    specs,
    inner_splits,
    inner_repeats,
    random_state,
    inner_validation_method,
    inner_max_pairs,
    standardize_predictors,
    standardize_columns,
    group_names,
    group_expressions,
    top_features,
    pair_batch_size,
    filter_scoring,
    risk_tie_tolerance,
):
    """Build grouped consensus using PLL improvement or held-out C-index.

    ``partial_likelihood`` preserves the previous behavior: each group's
    training-only Cox partial-log-likelihood improvement over the null is used.
    ``c_index`` fits the same group-only Cox model on inner training patients
    and scores only held-out validation patients. For unique-pair batches,
    concordance contributions are pooled across the held-out pairs before group
    ranking.
    """
    filter_scoring = _validate_consensus_filter_scoring(filter_scoring)
    counts = {name: 0 for name in group_names}
    scored_resamples = 0
    source_items = 0
    if specs:
        filter_penalizer = min(float(spec.penalizer) for spec in specs)
    else:
        filter_penalizer = 0.01

    # Each consensus resample contains one or more
    # (fit_positions, validation_positions) units. K-fold has one unit per
    # resample; unique-pair batching may have several.
    resample_units = []
    n_patients = len(train_data)
    if inner_validation_method == "kfold":
        for seed in _inner_repeat_seeds(random_state, inner_repeats):
            splits = _make_inner_splits(
                train_data,
                event_col=event_col,
                requested_splits=inner_splits,
                random_state=seed,
            )
            for fit_positions, validation_positions in splits:
                resample_units.append([
                    (
                        list(np.asarray(fit_positions, dtype=int)),
                        list(np.asarray(validation_positions, dtype=int)),
                    )
                ])
    else:
        pairs = _unique_harrell_inner_pairs(
            train_data,
            duration_col=duration_col,
            event_col=event_col,
            max_pairs=inner_max_pairs,
            random_state=random_state,
        )
        for pair_batch in _chunk_sequence(pairs, pair_batch_size):
            batch_units = []
            for i, j in pair_batch:
                excluded = {int(i), int(j)}
                fit_positions = [
                    position for position in range(n_patients)
                    if position not in excluded
                ]
                batch_units.append((fit_positions, [int(i), int(j)]))
            if batch_units:
                resample_units.append(batch_units)

    for units in resample_units:
        scores = []
        for name in group_names:
            if filter_scoring == "partial_likelihood":
                values = []
                for fit_positions, _validation_positions in units:
                    fit_data = train_data.iloc[fit_positions]
                    try:
                        value = _cox_filter_group_improvement(
                            fit_data,
                            covariates,
                            duration_col,
                            event_col,
                            group_expressions[name],
                            filter_penalizer,
                            standardize_predictors,
                            standardize_columns,
                        )
                        if np.isfinite(value):
                            values.append(float(value))
                    except Exception:
                        continue
                if values:
                    scores.append((float(np.mean(values)), str(name)))
            else:
                contribution_sum = 0.0
                comparable_sum = 0
                for fit_positions, validation_positions in units:
                    fit_data = train_data.iloc[fit_positions]
                    validation_data = train_data.iloc[validation_positions]
                    try:
                        contribution, comparable = _cox_filter_group_c_index_counts(
                            fit_data,
                            validation_data,
                            covariates,
                            duration_col,
                            event_col,
                            group_expressions[name],
                            filter_penalizer,
                            standardize_predictors,
                            standardize_columns,
                            risk_tie_tolerance,
                        )
                        if comparable > 0 and np.isfinite(contribution):
                            contribution_sum += float(contribution)
                            comparable_sum += int(comparable)
                    except Exception:
                        continue
                if comparable_sum > 0:
                    scores.append((
                        float(contribution_sum) / float(comparable_sum),
                        str(name),
                    ))
        if not scores:
            continue
        scores.sort(key=lambda item: (-item[0], item[1]))
        selected = scores[: min(int(top_features), len(scores))]
        scored_resamples += 1
        source_items += len(selected)
        for item in selected:
            counts[item[1]] += 1

    if scored_resamples <= 0:
        metric_label = (
            "held-out concordance" if filter_scoring == "c_index"
            else "training partial likelihood"
        )
        raise RuntimeError(
            "Cox-filter consensus stage could not score any inner resample "
            "using {}.".format(metric_label)
        )
    frequencies = {
        name: float(counts[name]) / float(scored_resamples)
        for name in group_names
    }
    return {
        "frequencies": frequencies,
        "resamples": int(scored_resamples),
        "source_items": int(source_items),
    }


def _annotate_consensus_ranking(
    ranking,
    consensus_method,
    selected_groups,
    group_names,
    frequencies,
    generated_formulas,
    resamples,
    source_items,
    pair_batch_size,
    filter_scoring,
):
    ranking = ranking.copy()
    ranking["consensus_method"] = consensus_method
    ranking["consensus_groups"] = ",".join(selected_groups)
    ranking["consensus_group_frequencies"] = _format_consensus_frequencies(
        group_names, frequencies
    )
    ranking["consensus_generated_formula_count"] = int(len(generated_formulas))
    ranking["consensus_resamples"] = int(resamples)
    ranking["consensus_source_items"] = int(source_items)
    ranking["consensus_pair_batch_size"] = int(pair_batch_size)
    ranking["consensus_filter_scoring"] = str(filter_scoring)
    return ranking


def _rank_candidates_inner_cv(
    train_data: pd.DataFrame,
    covariates: Sequence[str],
    duration_col: str,
    event_col: str,
    specs: Sequence[CoxModelSpec],
    inner_splits: int,
    inner_repeats: int,
    random_state: int,
    risk_tie_tolerance: float,
    inner_validation_method="kfold",
    inner_max_pairs=100,
    inner_scoring="c_index",
    inner_selection_rule="mean",
    inner_stability_penalty=1.0,
    standardize_predictors=False,
    standardize_columns=None,
    consensus_method="none",
    consensus_top_models=5,
    consensus_top_features=3,
    consensus_threshold=0.60,
    consensus_max_formula_groups=3,
    consensus_pair_batch_size=1,
    consensus_filter_scoring="partial_likelihood",
    consensus_feature_groups=None,
    consensus_group_expressions=None,
) -> pd.DataFrame:
    """Rank candidates, optionally after a Parvandeh-inspired consensus stage."""
    method = _validate_inner_validation_method(inner_validation_method)
    scoring = _validate_inner_scoring(inner_scoring)
    rule = _validate_inner_selection_rule(inner_selection_rule)
    consensus = _validate_consensus_method(consensus_method)
    (
        top_models,
        top_features,
        threshold,
        max_formula_groups,
    ) = _validate_consensus_controls(
        consensus_top_models,
        consensus_top_features,
        consensus_threshold,
        consensus_max_formula_groups,
    )
    pair_batch_size = _validate_consensus_pair_batch_size(consensus_pair_batch_size)
    filter_scoring = _validate_consensus_filter_scoring(consensus_filter_scoring)

    if consensus == "none":
        ranking = _rank_candidates_inner_cv_base(
            train_data=train_data,
            covariates=covariates,
            duration_col=duration_col,
            event_col=event_col,
            specs=specs,
            inner_splits=inner_splits,
            inner_repeats=inner_repeats,
            random_state=random_state,
            risk_tie_tolerance=risk_tie_tolerance,
            inner_validation_method=method,
            inner_max_pairs=inner_max_pairs,
            inner_scoring=scoring,
            inner_selection_rule=rule,
            inner_stability_penalty=inner_stability_penalty,
            standardize_predictors=standardize_predictors,
            standardize_columns=standardize_columns,
        )
        frequencies = {}
        return _annotate_consensus_ranking(
            ranking,
            "none",
            [],
            [],
            frequencies,
            [],
            0,
            0,
            pair_batch_size,
            filter_scoring,
        )

    group_names, group_members, group_expressions = _normalize_consensus_groups(
        covariates,
        consensus_feature_groups,
        consensus_group_expressions,
    )
    if not group_names:
        raise ValueError("Consensus selection requires at least one feature group.")

    if consensus == "top_models":
        information = _consensus_from_top_models(
            train_data,
            covariates,
            duration_col,
            event_col,
            specs,
            inner_splits,
            inner_repeats,
            random_state,
            risk_tie_tolerance,
            method,
            inner_max_pairs,
            scoring,
            rule,
            inner_stability_penalty,
            standardize_predictors,
            standardize_columns,
            group_names,
            group_members,
            top_models,
            pair_batch_size,
        )
    else:
        information = _consensus_from_cox_filter(
            train_data,
            covariates,
            duration_col,
            event_col,
            specs,
            inner_splits,
            inner_repeats,
            random_state,
            method,
            inner_max_pairs,
            standardize_predictors,
            standardize_columns,
            group_names,
            group_expressions,
            top_features,
            pair_batch_size,
            filter_scoring,
            risk_tie_tolerance,
        )

    frequencies = information["frequencies"]
    selected_groups = _select_consensus_groups(
        group_names,
        frequencies,
        threshold,
    )
    generated_formulas = _generate_consensus_formulas(
        selected_groups,
        frequencies,
        group_expressions,
        max_formula_groups,
    )
    generated_specs = _specs_for_generated_formulas(generated_formulas, specs)

    ranking = _rank_candidates_inner_cv_base(
        train_data=train_data,
        covariates=covariates,
        duration_col=duration_col,
        event_col=event_col,
        specs=generated_specs,
        inner_splits=inner_splits,
        inner_repeats=inner_repeats,
        random_state=random_state,
        risk_tie_tolerance=risk_tie_tolerance,
        inner_validation_method=method,
        inner_max_pairs=inner_max_pairs,
        inner_scoring=scoring,
        inner_selection_rule=rule,
        inner_stability_penalty=inner_stability_penalty,
        standardize_predictors=standardize_predictors,
        standardize_columns=standardize_columns,
    )
    return _annotate_consensus_ranking(
        ranking,
        consensus,
        selected_groups,
        group_names,
        frequencies,
        generated_formulas,
        information["resamples"],
        information["source_items"],
        pair_batch_size,
        filter_scoring,
    )

def _patient_bootstrap_worker(payload):
    """Module-level worker for deterministic patient-bootstrap replicates."""
    (seeds, edge_i, edge_j, contributions, n_patients) = payload
    probabilities = np.full(int(n_patients), 1.0 / float(n_patients))
    values = np.full(len(seeds), np.nan, dtype=float)
    for local_index, seed in enumerate(seeds):
        local_rng = np.random.default_rng(int(seed))
        multiplicity = local_rng.multinomial(int(n_patients), probabilities)
        weights = multiplicity[edge_i] * multiplicity[edge_j]
        denominator = float(weights.sum())
        if denominator > 0:
            values[local_index] = float(
                np.dot(weights, contributions) / denominator
            )
    return values


def _patient_bootstrap_ci(
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    contributions: np.ndarray,
    n_patients: int,
    n_bootstrap: int,
    confidence_level: float,
    rng: np.random.Generator,
    n_jobs: int = 1,
    multiprocessing_start_method: str = "fork",
) -> Tuple[Tuple[float, float], float, np.ndarray]:
    """Patient-resampling percentile CI for fixed out-of-sample pair scores.

    Replicate seeds are generated once in the parent process, so results are
    deterministic for a fixed random state and do not depend on the number of
    bootstrap worker processes.
    """
    if not 0 < confidence_level < 1:
        raise ValueError("confidence_level must be strictly between 0 and 1.")
    if n_bootstrap < 0:
        raise ValueError("n_bootstrap cannot be negative.")
    if n_bootstrap == 0:
        return (np.nan, np.nan), np.nan, np.empty(0, dtype=float)

    _, n_jobs_used = _resolve_n_jobs(n_jobs=n_jobs, n_tasks=n_bootstrap)
    maximum_seed = np.iinfo(np.int64).max
    seeds = rng.integers(0, maximum_seed, size=int(n_bootstrap), dtype=np.int64)

    if n_jobs_used == 1:
        values = _patient_bootstrap_worker(
            (seeds, edge_i, edge_j, contributions, n_patients)
        )
    else:
        if mp.current_process().daemon:
            raise RuntimeError(
                "Parallel bootstrap cannot start from a daemonic process."
            )
        # A few chunks per worker amortize Pool communication while retaining
        # enough work units for balanced scheduling.
        n_chunks = min(int(n_bootstrap), max(n_jobs_used, n_jobs_used * 4))
        seed_chunks = [chunk for chunk in np.array_split(seeds, n_chunks) if len(chunk)]
        payloads = [
            (chunk, edge_i, edge_j, contributions, n_patients)
            for chunk in seed_chunks
        ]
        context = mp.get_context(multiprocessing_start_method)
        pool = None
        try:
            pool = context.Pool(processes=n_jobs_used)
            pieces = pool.map(_patient_bootstrap_worker, payloads, chunksize=1)
            pool.close()
        except BaseException:
            if pool is not None:
                pool.terminate()
            raise
        finally:
            if pool is not None:
                pool.join()
        values = np.concatenate(pieces) if pieces else np.empty(0, dtype=float)

    valid = values[np.isfinite(values)]
    if len(valid) < max(20, int(0.8 * n_bootstrap)):
        warnings.warn(
            "Many patient-bootstrap replicates contained no sampled comparable "
            "edge; the confidence interval may be unstable.",
            RuntimeWarning,
            stacklevel=2,
        )
    if len(valid) < 2:
        return (np.nan, np.nan), np.nan, valid

    alpha = 1.0 - confidence_level
    lower, upper = np.quantile(valid, [alpha / 2.0, 1.0 - alpha / 2.0])
    standard_error = float(np.std(valid, ddof=1))
    return (float(lower), float(upper)), standard_error, valid



def _selection_frequency(outer_results: pd.DataFrame) -> pd.DataFrame:
    fitted = outer_results.loc[outer_results["fit_succeeded"]].copy()
    if fitted.empty:
        return pd.DataFrame(
            columns=[
                "formula",
                "base_model",
                "knot_number",
                "piecewise_break_number",
                "penalizer",
                "l1_ratio",
                "count",
                "proportion",
            ]
        )

    fitted["knot_number"] = fitted["knot_number"].fillna("-")
    fitted["piecewise_break_number"] = fitted["piecewise_break_number"].fillna("-")
    fitted["penalizer"] = fitted["penalizer"].fillna("-")
    fitted["l1_ratio"] = fitted["l1_ratio"].fillna("-")
    frequency = (
        fitted.groupby(
            [
                "formula",
                "base_model",
                "knot_number",
                "piecewise_break_number",
                "penalizer",
                "l1_ratio",
            ],
            observed=False,
        )
        .size()
        .rename("count")
        .reset_index()
    )
    frequency["proportion"] = frequency["count"] / len(fitted)
    return frequency.sort_values(
        ["count", "formula", "base_model"],
        ascending=[False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)


class _NativeThreadLimiter(object):
    """Temporarily limit native BLAS/OpenMP libraries to one thread.

    Process-level parallelism is used for the outer patient pairs. Limiting
    native numerical libraries to one thread per process prevents a pool of,
    for example, 32 Python workers from each creating many additional BLAS or
    OpenMP threads.
    """

    _ENVIRONMENT_VARIABLES = (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    )

    def __init__(self, enabled):
        self.enabled = bool(enabled)
        self._saved_environment = {}
        self._limiter = None

    def __enter__(self):
        if not self.enabled:
            return self

        for variable in self._ENVIRONMENT_VARIABLES:
            self._saved_environment[variable] = os.environ.get(variable)
            os.environ[variable] = "1"

        if threadpool_limits is None:
            warnings.warn(
                "threadpoolctl is unavailable. The multiprocessing worker count "
                "is capped, but native BLAS/OpenMP threads cannot be limited "
                "dynamically. scikit-learn 0.24.2 normally supplies "
                "threadpoolctl>=2.0.0.",
                RuntimeWarning,
                stacklevel=3,
            )
        else:
            try:
                # Construction applies the limit immediately. Keeping this object
                # alive preserves the limit until __exit__ restores it.
                self._limiter = threadpool_limits(limits=1)
            except Exception as exc:
                warnings.warn(
                    "Unable to apply threadpoolctl limits: {}: {}".format(
                        type(exc).__name__, exc
                    ),
                    RuntimeWarning,
                    stacklevel=3,
                )
                self._limiter = None
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._limiter is not None:
            restore = getattr(self._limiter, "restore_original_limits", None)
            if restore is None:
                restore = getattr(self._limiter, "unregister", None)
            if restore is not None:
                restore()

        if self.enabled:
            for variable, previous_value in self._saved_environment.items():
                if previous_value is None:
                    os.environ.pop(variable, None)
                else:
                    os.environ[variable] = previous_value
        return False


def _limit_native_threads_during_call(function):
    """Decorator used by the public evaluator without changing its signature."""

    @wraps(function)
    def wrapped(*args, **kwargs):
        with _NativeThreadLimiter(True):
            return function(*args, **kwargs)

    return wrapped


def _available_cpu_count():
    """Return the CPUs available to this process, respecting Linux affinity."""

    try:
        count = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        try:
            count = mp.cpu_count()
        except NotImplementedError:
            count = 1
    return max(1, int(count))


def _resolve_n_jobs(n_jobs, n_tasks):
    """Validate the requested process count and return requested/effective values."""

    if isinstance(n_jobs, bool) or not isinstance(n_jobs, Integral):
        raise TypeError("n_jobs must be an integer from 1 through 32.")
    requested = int(n_jobs)
    if requested < 1 or requested > MAX_MULTIPROCESSING_WORKERS:
        raise ValueError(
            "n_jobs must be between 1 and {} inclusive.".format(
                MAX_MULTIPROCESSING_WORKERS
            )
        )
    if n_tasks < 1:
        raise ValueError("At least one outer-pair task is required.")

    effective = min(
        requested,
        MAX_MULTIPROCESSING_WORKERS,
        _available_cpu_count(),
        int(n_tasks),
    )
    return requested, max(1, effective)


def _resolve_multiprocessing_start_method(start_method):
    """Validate and return a multiprocessing start method."""

    if not isinstance(start_method, str):
        raise TypeError("multiprocessing_start_method must be a string.")
    available = mp.get_all_start_methods()
    if start_method not in available:
        raise ValueError(
            "Unsupported multiprocessing_start_method {!r}. Available methods: "
            "{}".format(start_method, available)
        )
    return start_method


def _set_worker_native_thread_limit(enabled):
    """Apply one-native-thread limits inside a worker process."""

    global _WORKER_THREADPOOL_LIMITER

    if not enabled:
        _WORKER_THREADPOOL_LIMITER = None
        return

    for variable in _NativeThreadLimiter._ENVIRONMENT_VARIABLES:
        os.environ[variable] = "1"

    if threadpool_limits is not None:
        try:
            _WORKER_THREADPOOL_LIMITER = threadpool_limits(limits=1)
        except Exception:
            _WORKER_THREADPOOL_LIMITER = None


def _initialize_outer_pair_worker(
    work,
    covariates,
    original_index,
    duration_col,
    event_col,
    specs,
    inner_splits,
    inner_repeats,
    inner_validation_method,
    inner_max_pairs,
    inner_scoring,
    inner_selection_rule,
    inner_stability_penalty,
    standardize_predictors,
    standardize_columns,
    consensus_method,
    consensus_top_models,
    consensus_top_features,
    consensus_threshold,
    consensus_max_formula_groups,
    consensus_pair_batch_size,
    consensus_filter_scoring,
    consensus_feature_groups,
    consensus_group_expressions,
    prediction_aggregation,
    ensemble_size,
    ensemble_weighting,
    ensemble_temperature,
    ensemble_distinct_formulas,
    ensemble_rank_level,
    risk_tie_tolerance,
):
    """Initialize immutable state once in each multiprocessing worker."""
    global _OUTER_WORKER_STATE
    _set_worker_native_thread_limit(True)
    _OUTER_WORKER_STATE = {
        "work": work,
        "covariates": covariates,
        "original_index": original_index,
        "duration_col": duration_col,
        "event_col": event_col,
        "specs": specs,
        "inner_splits": inner_splits,
        "inner_repeats": inner_repeats,
        "inner_validation_method": inner_validation_method,
        "inner_max_pairs": inner_max_pairs,
        "inner_scoring": inner_scoring,
        "inner_selection_rule": inner_selection_rule,
        "inner_stability_penalty": inner_stability_penalty,
        "standardize_predictors": standardize_predictors,
        "standardize_columns": standardize_columns,
        "consensus_method": consensus_method,
        "consensus_top_models": consensus_top_models,
        "consensus_top_features": consensus_top_features,
        "consensus_threshold": consensus_threshold,
        "consensus_max_formula_groups": consensus_max_formula_groups,
        "consensus_pair_batch_size": consensus_pair_batch_size,
        "consensus_filter_scoring": consensus_filter_scoring,
        "consensus_feature_groups": consensus_feature_groups,
        "consensus_group_expressions": consensus_group_expressions,
        "prediction_aggregation": prediction_aggregation,
        "ensemble_size": ensemble_size,
        "ensemble_weighting": ensemble_weighting,
        "ensemble_temperature": ensemble_temperature,
        "ensemble_distinct_formulas": ensemble_distinct_formulas,
        "ensemble_rank_level": ensemble_rank_level,
        "risk_tie_tolerance": risk_tie_tolerance,
    }


def _make_outer_result_row(
    pair_number,
    i,
    j,
    inner_seed,
    work,
    original_index,
    duration_col,
    event_col,
):
    """Create a complete default result record for one held-out pair."""
    return {
        "pair_number": pair_number,
        "patient_position_i": i,
        "patient_position_j": j,
        "patient_id_i": original_index[i],
        "patient_id_j": original_index[j],
        "duration_i": float(work.at[i, duration_col]),
        "duration_j": float(work.at[j, duration_col]),
        "event_i": int(work.at[i, event_col]),
        "event_j": int(work.at[j, event_col]),
        "inner_random_state": int(inner_seed),
        "worker_pid": int(os.getpid()),
        "risk_i": np.nan,
        "risk_j": np.nan,
        "prediction_time": np.nan,
        "comparable": False,
        "concordance_contribution": np.nan,
        "fit_succeeded": False,
        "selected_inner_rank": np.nan,
        "formula": "",
        "base_model": "",
        "knot_number": np.nan,
        "piecewise_break_number": np.nan,
        "piecewise_breakpoints": "",
        "penalizer": np.nan,
        "l1_ratio": np.nan,
        "inner_validation_method": "",
        "inner_scoring": "",
        "inner_primary_score": np.nan,
        "inner_primary_score_std": np.nan,
        "inner_primary_score_se": np.nan,
        "inner_c_index": np.nan,
        "inner_c_index_std": np.nan,
        "inner_c_index_se": np.nan,
        "inner_cv_partial_log_likelihood": np.nan,
        "inner_cv_partial_log_likelihood_std": np.nan,
        "inner_cv_partial_log_likelihood_se": np.nan,
        "selection_score": np.nan,
        "consensus_method": "",
        "consensus_groups": "",
        "consensus_group_frequencies": "",
        "consensus_generated_formula_count": 0,
        "consensus_resamples": 0,
        "consensus_source_items": 0,
        "consensus_pair_batch_size": 1,
        "consensus_filter_scoring": "",
        "inner_repeats_requested": np.nan,
        "inner_repeats_scored": np.nan,
        "inner_comparable_pairs": 0,
        "inner_unique_pairs_requested": 0,
        "inner_unique_pairs_scored": 0,
        "standardized_predictors": "",
        "prediction_aggregation": "",
        "ensemble_weighting": "",
        "ensemble_requested_size": 0,
        "ensemble_fitted_size": 0,
        "ensemble_temperature": np.nan,
        "ensemble_distinct_formulas": False,
        "ensemble_rank_level": "",
        "ensemble_component_inner_ranks": "",
        "ensemble_component_weight_ranks": "",
        "ensemble_component_weights": "",
        "ensemble_component_labels": "",
        "ensemble_component_formulas": "",
        "ensemble_component_base_models": "",
        "ensemble_component_penalizers": "",
        "fit_warning_count": 0,
        "status": "failed",
        "error": "",
    }


def _evaluate_outer_pair(
    task,
    work,
    covariates,
    original_index,
    duration_col,
    event_col,
    specs,
    inner_splits,
    inner_repeats,
    inner_validation_method,
    inner_max_pairs,
    inner_scoring,
    inner_selection_rule,
    inner_stability_penalty,
    standardize_predictors,
    standardize_columns,
    consensus_method,
    consensus_top_models,
    consensus_top_features,
    consensus_threshold,
    consensus_max_formula_groups,
    consensus_pair_batch_size,
    consensus_filter_scoring,
    consensus_feature_groups,
    consensus_group_expressions,
    prediction_aggregation,
    ensemble_size,
    ensemble_weighting,
    ensemble_temperature,
    ensemble_distinct_formulas,
    ensemble_rank_level,
    risk_tie_tolerance,
):
    """Evaluate one outer pair; suitable for serial or worker-process use."""
    pair_number, i, j, inner_seed = task
    n_patients = len(work)
    excluded_positions = {i, j}
    train_positions = [p for p in range(n_patients) if p not in excluded_positions]
    outer_train = work.iloc[train_positions].copy()
    outer_test = work.iloc[[i, j]].copy()
    row = _make_outer_result_row(
        pair_number, i, j, inner_seed, work, original_index, duration_col, event_col
    )

    try:
        ranking = _rank_candidates_inner_cv(
            train_data=outer_train,
            covariates=covariates,
            duration_col=duration_col,
            event_col=event_col,
            specs=specs,
            inner_splits=inner_splits,
            inner_repeats=inner_repeats,
            random_state=inner_seed,
            risk_tie_tolerance=risk_tie_tolerance,
            inner_validation_method=inner_validation_method,
            inner_max_pairs=inner_max_pairs,
            inner_scoring=inner_scoring,
            inner_selection_rule=inner_selection_rule,
            inner_stability_penalty=inner_stability_penalty,
            standardize_predictors=standardize_predictors,
            standardize_columns=standardize_columns,
            consensus_method=consensus_method,
            consensus_top_models=consensus_top_models,
            consensus_top_features=consensus_top_features,
            consensus_threshold=consensus_threshold,
            consensus_max_formula_groups=consensus_max_formula_groups,
            consensus_pair_batch_size=consensus_pair_batch_size,
            consensus_filter_scoring=consensus_filter_scoring,
            consensus_feature_groups=consensus_feature_groups,
            consensus_group_expressions=consensus_group_expressions,
        )
        eligible = ranking.loc[ranking["eligible"]]
        if eligible.empty:
            errors = ranking.loc[ranking["first_error"] != "", "first_error"]
            detail = errors.iloc[0] if not errors.empty else "no comparable inner pairs"
            raise RuntimeError("No candidate passed inner validation ({}).".format(detail))

        if prediction_aggregation == "stacked":
            ensemble_model, components, weights, risks, outer_errors = _fit_prediction_ensemble(
                ranking=ranking,
                train_data=outer_train,
                evaluation_data=outer_test,
                covariates=covariates,
                duration_col=duration_col,
                event_col=event_col,
                standardize_predictors=standardize_predictors,
                standardize_columns=standardize_columns,
                ensemble_size=ensemble_size,
                ensemble_weighting=ensemble_weighting,
                ensemble_temperature=ensemble_temperature,
                ensemble_distinct_formulas=ensemble_distinct_formulas,
                ensemble_rank_level=ensemble_rank_level,
            )
            lead = components[0]["candidate_row"]
            is_comparable, contribution = _pair_contribution(
                row["duration_i"], row["event_i"], float(risks[0]),
                row["duration_j"], row["event_j"], float(risks[1]),
                risk_tie_tolerance,
            )
            standardized_union = set()
            warning_count = 0
            for item in components:
                standardized_union.update(
                    getattr(item["model"], "_nested_standardized_columns", ())
                )
                warning_count += len(item["warnings"])
            row.update(
                {
                    "risk_i": float(risks[0]),
                    "risk_j": float(risks[1]),
                    "prediction_time": _reference_prediction_time(
                        outer_train, duration_col, event_col
                    ),
                    "comparable": bool(is_comparable),
                    "concordance_contribution": contribution,
                    "fit_succeeded": True,
                    "selected_inner_rank": int(lead["inner_rank"]),
                    "formula": "STACKED[{}]".format(len(components)),
                    "base_model": "ensemble",
                    "knot_number": np.nan,
                    "piecewise_break_number": np.nan,
                    "piecewise_breakpoints": "",
                    "penalizer": np.nan,
                    "l1_ratio": np.nan,
                    "inner_validation_method": lead["inner_validation_method"],
                    "inner_scoring": lead["inner_scoring"],
                    "inner_primary_score": float(lead["inner_primary_score"]),
                    "inner_primary_score_std": float(lead["inner_primary_score_std"]),
                    "inner_primary_score_se": float(lead["inner_primary_score_se"]),
                    "inner_c_index": float(lead["inner_c_index"]),
                    "inner_c_index_std": float(lead["inner_c_index_std"]),
                    "inner_c_index_se": float(lead["inner_c_index_se"]),
                    "inner_cv_partial_log_likelihood": float(
                        lead["inner_cv_partial_log_likelihood"]
                    ),
                    "inner_cv_partial_log_likelihood_std": float(
                        lead["inner_cv_partial_log_likelihood_std"]
                    ),
                    "inner_cv_partial_log_likelihood_se": float(
                        lead["inner_cv_partial_log_likelihood_se"]
                    ),
                    "selection_score": float(lead["selection_score"]),
                    "consensus_method": str(lead["consensus_method"]),
                    "consensus_groups": str(lead["consensus_groups"]),
                    "consensus_group_frequencies": str(lead["consensus_group_frequencies"]),
                    "consensus_generated_formula_count": int(
                        lead["consensus_generated_formula_count"]
                    ),
                    "consensus_resamples": int(lead["consensus_resamples"]),
                    "consensus_source_items": int(lead["consensus_source_items"]),
                    "consensus_pair_batch_size": int(lead["consensus_pair_batch_size"]),
                    "consensus_filter_scoring": str(lead["consensus_filter_scoring"]),
                    "inner_repeats_requested": int(lead["inner_repeats_requested"]),
                    "inner_repeats_scored": int(lead["inner_repeats_scored"]),
                    "inner_comparable_pairs": int(lead["inner_comparable_pairs"]),
                    "inner_unique_pairs_requested": int(lead["inner_unique_pairs_requested"]),
                    "inner_unique_pairs_scored": int(lead["inner_unique_pairs_scored"]),
                    "standardized_predictors": ",".join(sorted(standardized_union)),
                    "prediction_aggregation": "stacked",
                    "ensemble_weighting": ensemble_weighting,
                    "ensemble_requested_size": int(ensemble_size),
                    "ensemble_fitted_size": int(len(components)),
                    "ensemble_temperature": float(ensemble_temperature),
                    "ensemble_distinct_formulas": bool(ensemble_distinct_formulas),
                    "ensemble_rank_level": str(ensemble_rank_level),
                    "ensemble_component_inner_ranks": ";".join(
                        str(int(item["candidate_row"]["inner_rank"])) for item in components
                    ),
                    "ensemble_component_weight_ranks": ";".join(
                        "{:.12g}".format(float(item.get("ensemble_weight_rank", np.nan)))
                        for item in components
                    ),
                    "ensemble_component_weights": ";".join(
                        "{:.12g}".format(float(value)) for value in weights
                    ),
                    "ensemble_component_labels": " || ".join(
                        item["spec"].label for item in components
                    ),
                    "ensemble_component_formulas": " || ".join(
                        item["spec"].formula_label for item in components
                    ),
                    "ensemble_component_base_models": ";".join(
                        item["spec"].base_model for item in components
                    ),
                    "ensemble_component_penalizers": ";".join(
                        "{:.12g}".format(float(item["spec"].penalizer)) for item in components
                    ),
                    "fit_warning_count": int(warning_count),
                    "status": "ok" if is_comparable else "incomparable",
                }
            )
            return row

        outer_errors = []
        for _, candidate_row in eligible.iterrows():
            spec = candidate_row["spec"]
            try:
                model, fit_warnings = _fit_cox_model(
                    outer_train,
                    duration_col=duration_col,
                    event_col=event_col,
                    spec=spec,
                    covariates=covariates,
                    standardize_predictors=standardize_predictors,
                    standardize_columns=standardize_columns,
                )
                risks, prediction_time = _predict_risk(
                    model,
                    train_data=outer_train,
                    test_data=outer_test,
                    covariates=covariates,
                    duration_col=duration_col,
                    event_col=event_col,
                )
                is_comparable, contribution = _pair_contribution(
                    row["duration_i"], row["event_i"], float(risks[0]),
                    row["duration_j"], row["event_j"], float(risks[1]),
                    risk_tie_tolerance,
                )
                row.update(
                    {
                        "risk_i": float(risks[0]),
                        "risk_j": float(risks[1]),
                        "prediction_time": prediction_time,
                        "comparable": bool(is_comparable),
                        "concordance_contribution": contribution,
                        "fit_succeeded": True,
                        "selected_inner_rank": int(candidate_row["inner_rank"]),
                        "formula": spec.formula_label,
                        "base_model": spec.base_model,
                        "knot_number": spec.knot_number,
                        "piecewise_break_number": spec.piecewise_break_number,
                        "piecewise_breakpoints": (
                            ",".join(
                                "{:.12g}".format(float(value))
                                for value in getattr(model, "breakpoints", [])
                            ) if spec.base_model == "piecewise" else ""
                        ),
                        "penalizer": spec.penalizer,
                        "l1_ratio": spec.l1_ratio,
                        "inner_validation_method": candidate_row["inner_validation_method"],
                        "inner_scoring": candidate_row["inner_scoring"],
                        "inner_primary_score": float(candidate_row["inner_primary_score"]),
                        "inner_primary_score_std": float(candidate_row["inner_primary_score_std"]),
                        "inner_primary_score_se": float(candidate_row["inner_primary_score_se"]),
                        "inner_c_index": float(candidate_row["inner_c_index"]),
                        "inner_c_index_std": float(candidate_row["inner_c_index_std"]),
                        "inner_c_index_se": float(candidate_row["inner_c_index_se"]),
                        "inner_cv_partial_log_likelihood": float(
                            candidate_row["inner_cv_partial_log_likelihood"]
                        ),
                        "inner_cv_partial_log_likelihood_std": float(
                            candidate_row["inner_cv_partial_log_likelihood_std"]
                        ),
                        "inner_cv_partial_log_likelihood_se": float(
                            candidate_row["inner_cv_partial_log_likelihood_se"]
                        ),
                        "selection_score": float(candidate_row["selection_score"]),
                        "consensus_method": str(candidate_row["consensus_method"]),
                        "consensus_groups": str(candidate_row["consensus_groups"]),
                        "consensus_group_frequencies": str(
                            candidate_row["consensus_group_frequencies"]
                        ),
                        "consensus_generated_formula_count": int(
                            candidate_row["consensus_generated_formula_count"]
                        ),
                        "consensus_resamples": int(candidate_row["consensus_resamples"]),
                        "consensus_source_items": int(candidate_row["consensus_source_items"]),
                        "consensus_pair_batch_size": int(candidate_row["consensus_pair_batch_size"]),
                        "consensus_filter_scoring": str(candidate_row["consensus_filter_scoring"]),
                        "inner_repeats_requested": int(candidate_row["inner_repeats_requested"]),
                        "inner_repeats_scored": int(candidate_row["inner_repeats_scored"]),
                        "inner_comparable_pairs": int(candidate_row["inner_comparable_pairs"]),
                        "inner_unique_pairs_requested": int(candidate_row["inner_unique_pairs_requested"]),
                        "inner_unique_pairs_scored": int(candidate_row["inner_unique_pairs_scored"]),
                        "standardized_predictors": ",".join(
                            getattr(model, "_nested_standardized_columns", ())
                        ),
                        "prediction_aggregation": "winner",
                        "ensemble_weighting": "",
                        "ensemble_requested_size": 1,
                        "ensemble_fitted_size": 1,
                        "ensemble_temperature": np.nan,
                        "ensemble_distinct_formulas": False,
                        "ensemble_rank_level": str(ensemble_rank_level),
                        "ensemble_component_inner_ranks": str(int(candidate_row["inner_rank"])),
                        "ensemble_component_weight_ranks": "1",
                        "ensemble_component_weights": "1",
                        "ensemble_component_labels": spec.label,
                        "ensemble_component_formulas": spec.formula_label,
                        "ensemble_component_base_models": spec.base_model,
                        "ensemble_component_penalizers": "{:.12g}".format(float(spec.penalizer)),
                        "fit_warning_count": len(fit_warnings),
                        "status": "ok" if is_comparable else "incomparable",
                    }
                )
                break
            except Exception as exc:
                outer_errors.append(
                    "{}: {}: {}".format(spec.label, type(exc).__name__, exc)
                )
        if not row["fit_succeeded"]:
            raise RuntimeError(
                "Every inner-eligible candidate failed on the full outer training set. "
                + " | ".join(outer_errors[:3])
            )
    except Exception as exc:
        row["error"] = "{}: {}".format(type(exc).__name__, exc)
    return row


def _outer_pair_worker(task):
    """Multiprocessing entry point. It must remain defined at module level."""
    if _OUTER_WORKER_STATE is None:
        raise RuntimeError("Outer-pair worker state was not initialized.")
    state = _OUTER_WORKER_STATE
    return _evaluate_outer_pair(
        task=task,
        work=state["work"],
        covariates=state["covariates"],
        original_index=state["original_index"],
        duration_col=state["duration_col"],
        event_col=state["event_col"],
        specs=state["specs"],
        inner_splits=state["inner_splits"],
        inner_repeats=state["inner_repeats"],
        inner_validation_method=state["inner_validation_method"],
        inner_max_pairs=state["inner_max_pairs"],
        inner_scoring=state["inner_scoring"],
        inner_selection_rule=state["inner_selection_rule"],
        inner_stability_penalty=state["inner_stability_penalty"],
        standardize_predictors=state["standardize_predictors"],
        standardize_columns=state["standardize_columns"],
        consensus_method=state["consensus_method"],
        consensus_top_models=state["consensus_top_models"],
        consensus_top_features=state["consensus_top_features"],
        consensus_threshold=state["consensus_threshold"],
        consensus_max_formula_groups=state["consensus_max_formula_groups"],
        consensus_pair_batch_size=state["consensus_pair_batch_size"],
        consensus_filter_scoring=state["consensus_filter_scoring"],
        consensus_feature_groups=state["consensus_feature_groups"],
        consensus_group_expressions=state["consensus_group_expressions"],
        prediction_aggregation=state["prediction_aggregation"],
        ensemble_size=state["ensemble_size"],
        ensemble_weighting=state["ensemble_weighting"],
        ensemble_temperature=state["ensemble_temperature"],
        ensemble_distinct_formulas=state["ensemble_distinct_formulas"],
        ensemble_rank_level=state["ensemble_rank_level"],
        risk_tie_tolerance=state["risk_tie_tolerance"],
    )


def _should_report_progress(completed, total):
    return completed == 1 or completed % 10 == 0 or completed == total


def _run_outer_pair_evaluations(
    tasks,
    work,
    covariates,
    original_index,
    duration_col,
    event_col,
    specs,
    inner_splits,
    inner_repeats,
    inner_validation_method,
    inner_max_pairs,
    inner_scoring,
    inner_selection_rule,
    inner_stability_penalty,
    standardize_predictors,
    standardize_columns,
    consensus_method,
    consensus_top_models,
    consensus_top_features,
    consensus_threshold,
    consensus_max_formula_groups,
    consensus_pair_batch_size,
    consensus_filter_scoring,
    consensus_feature_groups,
    consensus_group_expressions,
    prediction_aggregation,
    ensemble_size,
    ensemble_weighting,
    ensemble_temperature,
    ensemble_distinct_formulas,
    ensemble_rank_level,
    risk_tie_tolerance,
    n_jobs_used,
    multiprocessing_start_method,
    verbose,
):
    """Run pair tasks serially or through a multiprocessing.Pool."""
    total = len(tasks)
    rows = []
    common = dict(
        work=work,
        covariates=covariates,
        original_index=original_index,
        duration_col=duration_col,
        event_col=event_col,
        specs=specs,
        inner_splits=inner_splits,
        inner_repeats=inner_repeats,
        inner_validation_method=inner_validation_method,
        inner_max_pairs=inner_max_pairs,
        inner_scoring=inner_scoring,
        inner_selection_rule=inner_selection_rule,
        inner_stability_penalty=inner_stability_penalty,
        standardize_predictors=standardize_predictors,
        standardize_columns=standardize_columns,
        consensus_method=consensus_method,
        consensus_top_models=consensus_top_models,
        consensus_top_features=consensus_top_features,
        consensus_threshold=consensus_threshold,
        consensus_max_formula_groups=consensus_max_formula_groups,
        consensus_pair_batch_size=consensus_pair_batch_size,
        consensus_filter_scoring=consensus_filter_scoring,
        consensus_feature_groups=consensus_feature_groups,
        consensus_group_expressions=consensus_group_expressions,
        prediction_aggregation=prediction_aggregation,
        ensemble_size=ensemble_size,
        ensemble_weighting=ensemble_weighting,
        ensemble_temperature=ensemble_temperature,
        ensemble_distinct_formulas=ensemble_distinct_formulas,
        ensemble_rank_level=ensemble_rank_level,
        risk_tie_tolerance=risk_tie_tolerance,
    )
    if n_jobs_used == 1:
        for completed, task in enumerate(tasks, start=1):
            row = _evaluate_outer_pair(task=task, **common)
            rows.append(row)
            if verbose and _should_report_progress(completed, total):
                print(
                    "Completed outer pair {}/{}: pair_number={}, status={}".format(
                        completed, total, row["pair_number"], row["status"]
                    ),
                    flush=True,
                )
        return rows

    if mp.current_process().daemon:
        raise RuntimeError(
            "n_jobs>1 cannot be used from a daemonic process because daemonic processes cannot create multiprocessing children."
        )
    context = mp.get_context(multiprocessing_start_method)
    pool = None
    try:
        pool = context.Pool(
            processes=n_jobs_used,
            initializer=_initialize_outer_pair_worker,
            initargs=(
                work, list(covariates), original_index, duration_col, event_col,
                list(specs), inner_splits, inner_repeats, inner_validation_method,
                inner_max_pairs, inner_scoring, inner_selection_rule, inner_stability_penalty,
                standardize_predictors, list(standardize_columns),
                consensus_method, consensus_top_models, consensus_top_features,
                consensus_threshold, consensus_max_formula_groups,
                consensus_pair_batch_size, consensus_filter_scoring,
                consensus_feature_groups, consensus_group_expressions,
                prediction_aggregation, ensemble_size, ensemble_weighting,
                ensemble_temperature, ensemble_distinct_formulas, ensemble_rank_level,
                risk_tie_tolerance,
            ),
        )
        iterator = pool.imap_unordered(_outer_pair_worker, tasks, chunksize=1)
        for completed, row in enumerate(iterator, start=1):
            rows.append(row)
            if verbose and _should_report_progress(completed, total):
                print(
                    "Completed outer pair {}/{}: pair_number={}, status={}, worker_pid={}".format(
                        completed, total, row["pair_number"], row["status"], row["worker_pid"]
                    ),
                    flush=True,
                )
        pool.close()
    except BaseException:
        if pool is not None:
            pool.terminate()
        raise
    finally:
        if pool is not None:
            pool.join()
    rows.sort(key=lambda row: row["pair_number"])
    return rows


@_limit_native_threads_during_call
def nested_pair_cox_evaluation(
    data: pd.DataFrame,
    *,
    duration_col: str = "Response_duration",
    event_col: str = "Response_coded",
    formulas: Sequence[Optional[str]] = (None,),
    base_models: Sequence[str] = ("breslow", "spline"),
    knot_numbers: Sequence[int] = (2, 3, 4),
    piecewise_break_numbers: Sequence[int] = (3,),
    penalizers: Sequence[float] = (0.01,),
    l1_ratios: Sequence[float] = (0.0,),
    exclude_covariates: Sequence[str] = (),
    inner_splits: int = 4,
    inner_repeats: int = 1,
    inner_validation_method: str = "kfold",
    inner_max_pairs: int = 100,
    inner_scoring: str = "c_index",
    inner_selection_rule: str = "mean",
    inner_stability_penalty: float = 1.0,
    standardize_predictors: bool = False,
    standardize_columns: Optional[Sequence[str]] = None,
    consensus_method: str = "none",
    consensus_top_models: int = 5,
    consensus_top_features: int = 3,
    consensus_threshold: float = 0.60,
    consensus_max_formula_groups: int = 3,
    consensus_pair_batch_size: int = 1,
    consensus_filter_scoring: str = "partial_likelihood",
    consensus_feature_groups=None,
    consensus_group_expressions=None,
    prediction_aggregation: str = "winner",
    ensemble_size: int = 5,
    ensemble_weighting: str = "equal",
    ensemble_temperature: float = 1.0,
    ensemble_distinct_formulas: bool = True,
    ensemble_rank_level: str = "specification",
    max_outer_pairs: int = 120,
    event_policy: str = "all",
    outer_pairs=None,
    n_jobs: int = 1,
    multiprocessing_start_method: str = "fork",
    n_bootstrap: int = 2000,
    confidence_level: float = 0.95,
    random_state: Optional[int] = 42,
    risk_tie_tolerance: float = 1e-12,
    stabilize_categorical_levels: bool = True,
    fit_final_model: bool = True,
    verbose: bool = False,
) -> CoxEvaluationResult:
    """Estimate nested out-of-sample concordance for Cox models.

    New stability controls
    ----------------------
    ``outer_pairs`` may contain a previously generated pair plan, guaranteeing
    identical outer test pairs across competing analyses. ``standardize_predictors``
    performs z-scoring using each fit's training patients only. Set
    ``inner_validation_method='unique_pairs'`` to tune on a deterministic unique
    sample of Harrell-comparable inner leave-pair-out tests. With
    ``inner_scoring='cv_partial_likelihood'`` selects models using grouped
    Verweij-van Houwelingen cross-validated partial log-likelihood. With
    ``inner_selection_rule='lower_confidence'`` candidates are ranked by the
    selected metric minus ``inner_stability_penalty * SE``.

    Parvandeh-inspired consensus controls
    --------------------------------------
    ``consensus_method='top_models'`` retains the top models within each inner
    resample, estimates grouped feature-selection frequencies, generates a
    reduced formula family from stable groups, and then reruns the ordinary
    inner model-selection procedure on that family.

    ``consensus_method='cox_filter'`` supports two group-scoring rules through
    ``consensus_filter_scoring``. ``partial_likelihood`` preserves training-only
    Cox partial-log-likelihood improvement over the null. ``c_index`` fits each
    group-only Cox model on inner training patients and ranks groups by held-out
    Harrell concordance. The most frequently retained groups define the reduced
    formula family. All consensus operations remain inside the outer training data.

    With ``inner_validation_method='unique_pairs'``,
    ``consensus_pair_batch_size`` pools several distinct inner held-out pairs into
    each consensus resample. In ``top_models`` mode candidate models are ranked
    on the pooled batch. In ``cox_filter`` mode, partial-likelihood scores are
    averaged across pair-specific training sets, whereas C-index scores pool
    held-out concordance contributions across all pairs in the batch. A batch
    size of 1 preserves the previous resampling granularity.

    Prediction averaging
    --------------------
    ``prediction_aggregation='stacked'`` fits the top eligible inner-ranked
    Cox models on each full outer training set, converts every component to a
    training-standardized log-risk score, and averages those scores. ``equal``
    weights are the most conservative default; ``rank`` uses inverse inner rank
    and ``softmax`` uses standardized inner selection scores. With rank weighting, ``ensemble_rank_level='formula'`` uses dense formula ranks rather than raw hyperparameter-specification ranks. No ensemble weight
    is estimated from the held-out outer pair. ``winner`` preserves the prior
    winner-takes-all behavior.
    """
    _require_lifelines()
    event_policy = _validate_event_policy(event_policy)
    method = _validate_inner_validation_method(inner_validation_method)
    scoring = _validate_inner_scoring(inner_scoring)
    selection_rule = _validate_inner_selection_rule(inner_selection_rule)
    consensus_mode = _validate_consensus_method(consensus_method)
    (
        consensus_top_models,
        consensus_top_features,
        consensus_threshold,
        consensus_max_formula_groups,
    ) = _validate_consensus_controls(
        consensus_top_models,
        consensus_top_features,
        consensus_threshold,
        consensus_max_formula_groups,
    )
    consensus_pair_batch_size = _validate_consensus_pair_batch_size(
        consensus_pair_batch_size
    )
    consensus_filter_scoring = _validate_consensus_filter_scoring(
        consensus_filter_scoring
    )
    prediction_mode = _validate_prediction_aggregation(prediction_aggregation)
    ensemble_weighting = _validate_ensemble_weighting(ensemble_weighting)
    ensemble_rank_level = _validate_ensemble_rank_level(ensemble_rank_level)
    ensemble_size, ensemble_temperature = _validate_ensemble_controls(
        ensemble_size, ensemble_temperature
    )
    ensemble_distinct_formulas = bool(ensemble_distinct_formulas)
    if isinstance(inner_repeats, bool) or int(inner_repeats) != inner_repeats:
        raise TypeError("inner_repeats must be a positive integer.")
    inner_repeats = int(inner_repeats)
    if inner_repeats < 1:
        raise ValueError("inner_repeats must be >= 1.")
    if isinstance(inner_max_pairs, bool) or int(inner_max_pairs) != inner_max_pairs:
        raise TypeError("inner_max_pairs must be a positive integer.")
    inner_max_pairs = int(inner_max_pairs)
    if inner_max_pairs < 1:
        raise ValueError("inner_max_pairs must be >= 1.")
    if float(inner_stability_penalty) < 0:
        raise ValueError("inner_stability_penalty must be non-negative.")
    if risk_tie_tolerance < 0:
        raise ValueError("risk_tie_tolerance must be non-negative.")
    if method == "unique_pairs" and inner_repeats != 1:
        warnings.warn(
            "inner_repeats is ignored when inner_validation_method='unique_pairs'; "
            "each selected inner pair is evaluated exactly once.",
            RuntimeWarning,
            stacklevel=2,
        )
    multiprocessing_start_method = _resolve_multiprocessing_start_method(
        multiprocessing_start_method
    )

    work, covariates, original_index = _validate_and_prepare_data(
        data=data,
        duration_col=duration_col,
        event_col=event_col,
        exclude_covariates=exclude_covariates,
        stabilize_categorical_levels=stabilize_categorical_levels,
    )
    standardized_columns = _resolve_standardize_columns(
        work,
        covariates,
        standardize_predictors=standardize_predictors,
        standardize_columns=standardize_columns,
    )
    normalized_consensus_groups = consensus_feature_groups
    normalized_consensus_expressions = consensus_group_expressions
    if consensus_mode != "none":
        _group_names, _group_members, _group_expressions = _normalize_consensus_groups(
            covariates, consensus_feature_groups, consensus_group_expressions
        )
        normalized_consensus_groups = {
            name: tuple(_group_members[name]) for name in _group_names
        }
        normalized_consensus_expressions = {
            name: _group_expressions[name] for name in _group_names
        }
    specs = _build_candidate_specs(
        formulas=formulas,
        base_models=base_models,
        knot_numbers=knot_numbers,
        piecewise_break_numbers=piecewise_break_numbers,
        penalizers=penalizers,
        l1_ratios=l1_ratios,
    )

    seed_sequence = np.random.SeedSequence(random_state)
    pair_seed, inner_seed_seed, ci_seed, final_seed = seed_sequence.spawn(4)
    pair_rng = np.random.default_rng(pair_seed)
    inner_seed_rng = np.random.default_rng(inner_seed_seed)
    ci_rng = np.random.default_rng(ci_seed)
    final_rng = np.random.default_rng(final_seed)

    n_patients = len(work)
    total_possible_pairs = _number_of_pairs(n_patients)
    event_values = work[event_col].to_numpy(dtype=int)
    duration_values = work[duration_col].to_numpy(dtype=float)
    total_policy_eligible_pairs = _count_policy_eligible_pairs(
        n_patients=n_patients,
        event_values=event_values,
        event_policy=event_policy,
        duration_values=duration_values,
    )
    if outer_pairs is None:
        selected_outer_pairs = _choose_outer_pairs(
            n_patients=n_patients,
            max_pairs=max_outer_pairs,
            rng=pair_rng,
            event_values=event_values,
            duration_values=duration_values,
            event_policy=event_policy,
        )
        outer_pairs_reused = False
    else:
        selected_outer_pairs = _normalize_reusable_outer_pairs(
            outer_pairs,
            work,
            original_index,
            duration_col,
            event_col,
            event_policy,
        )
        outer_pairs_reused = True

    pair_plan = _outer_pair_plan_from_pairs(
        selected_outer_pairs,
        work,
        original_index,
        duration_col,
        event_col,
    )
    n_jobs_requested, n_jobs_used = _resolve_n_jobs(
        n_jobs=n_jobs,
        n_tasks=len(selected_outer_pairs),
    )
    inner_seeds = [
        int(inner_seed_rng.integers(0, np.iinfo(np.int32).max))
        for _ in selected_outer_pairs
    ]
    outer_tasks = [
        (pair_number, pair[0], pair[1], inner_seeds[pair_number - 1])
        for pair_number, pair in enumerate(selected_outer_pairs, start=1)
    ]
    outer_rows = _run_outer_pair_evaluations(
        tasks=outer_tasks,
        work=work,
        covariates=covariates,
        original_index=original_index,
        duration_col=duration_col,
        event_col=event_col,
        specs=specs,
        inner_splits=inner_splits,
        inner_repeats=inner_repeats,
        inner_validation_method=method,
        inner_max_pairs=inner_max_pairs,
        inner_scoring=scoring,
        inner_selection_rule=selection_rule,
        inner_stability_penalty=float(inner_stability_penalty),
        standardize_predictors=bool(standardize_predictors),
        standardize_columns=standardized_columns,
        consensus_method=consensus_mode,
        consensus_top_models=consensus_top_models,
        consensus_top_features=consensus_top_features,
        consensus_threshold=consensus_threshold,
        consensus_max_formula_groups=consensus_max_formula_groups,
        consensus_pair_batch_size=consensus_pair_batch_size,
        consensus_filter_scoring=consensus_filter_scoring,
        consensus_feature_groups=normalized_consensus_groups,
        consensus_group_expressions=normalized_consensus_expressions,
        prediction_aggregation=prediction_mode,
        ensemble_size=ensemble_size,
        ensemble_weighting=ensemble_weighting,
        ensemble_temperature=ensemble_temperature,
        ensemble_distinct_formulas=ensemble_distinct_formulas,
        ensemble_rank_level=ensemble_rank_level,
        risk_tie_tolerance=risk_tie_tolerance,
        n_jobs_used=n_jobs_used,
        multiprocessing_start_method=multiprocessing_start_method,
        verbose=verbose,
    )

    outer_results = pd.DataFrame(outer_rows)
    fitted_mask = outer_results["fit_succeeded"].astype(bool)
    comparable_mask = fitted_mask & outer_results["comparable"].astype(bool)
    comparable_results = outer_results.loc[comparable_mask].copy()
    if comparable_results.empty:
        representative_errors = outer_results.loc[
            outer_results["error"] != "", "error"
        ].head(3)
        detail = " | ".join(representative_errors.tolist())
        raise RuntimeError(
            "No outer test pair was both successfully fitted and comparable. " + detail
        )

    c_index = float(comparable_results["concordance_contribution"].mean())
    edge_i = comparable_results["patient_position_i"].to_numpy(dtype=int)
    edge_j = comparable_results["patient_position_j"].to_numpy(dtype=int)
    contributions = comparable_results["concordance_contribution"].to_numpy(dtype=float)
    _bootstrap_requested, bootstrap_n_jobs_used = _resolve_n_jobs(
        n_jobs=n_jobs, n_tasks=max(1, int(n_bootstrap))
    )
    if verbose and n_bootstrap > 0:
        print(
            "Starting patient bootstrap: {} replicates with {} worker(s).".format(
                int(n_bootstrap), int(bootstrap_n_jobs_used)
            ),
            flush=True,
        )
    bootstrap_started_at = time.time()
    confidence_interval, bootstrap_se, bootstrap_values = _patient_bootstrap_ci(
        edge_i=edge_i,
        edge_j=edge_j,
        contributions=contributions,
        n_patients=n_patients,
        n_bootstrap=n_bootstrap,
        confidence_level=confidence_level,
        rng=ci_rng,
        n_jobs=bootstrap_n_jobs_used,
        multiprocessing_start_method=multiprocessing_start_method,
    )
    if verbose and n_bootstrap > 0:
        print(
            "Completed patient bootstrap in {:.3f} seconds.".format(
                time.time() - bootstrap_started_at
            ),
            flush=True,
        )

    failed_outer_pairs = int((~fitted_mask).sum())
    if failed_outer_pairs:
        warnings.warn(
            "{} of {} outer pairs failed and were excluded. Inspect outer_results['error'].".format(
                failed_outer_pairs, len(outer_results)
            ),
            RuntimeWarning,
            stacklevel=2,
        )

    final_model = None
    final_spec = None
    final_ranking = pd.DataFrame()
    if fit_final_model:
        if verbose:
            print(
                "Starting final full-data inner ranking/final-model fit (serial stage).",
                flush=True,
            )
        final_fit_started_at = time.time()
        final_ranking = _rank_candidates_inner_cv(
            train_data=work,
            covariates=covariates,
            duration_col=duration_col,
            event_col=event_col,
            specs=specs,
            inner_splits=inner_splits,
            inner_repeats=inner_repeats,
            random_state=int(final_rng.integers(0, np.iinfo(np.int32).max)),
            risk_tie_tolerance=risk_tie_tolerance,
            inner_validation_method=method,
            inner_max_pairs=inner_max_pairs,
            inner_scoring=scoring,
            inner_selection_rule=selection_rule,
            inner_stability_penalty=float(inner_stability_penalty),
            standardize_predictors=bool(standardize_predictors),
            standardize_columns=standardized_columns,
            consensus_method=consensus_mode,
            consensus_top_models=consensus_top_models,
            consensus_top_features=consensus_top_features,
            consensus_threshold=consensus_threshold,
            consensus_max_formula_groups=consensus_max_formula_groups,
            consensus_pair_batch_size=consensus_pair_batch_size,
            consensus_filter_scoring=consensus_filter_scoring,
            consensus_feature_groups=normalized_consensus_groups,
            consensus_group_expressions=normalized_consensus_expressions,
        )
        if prediction_mode == "stacked":
            try:
                final_model, final_components, _final_weights, _full_scores, _ = (
                    _fit_prediction_ensemble(
                        ranking=final_ranking,
                        train_data=work,
                        evaluation_data=work,
                        covariates=covariates,
                        duration_col=duration_col,
                        event_col=event_col,
                        standardize_predictors=bool(standardize_predictors),
                        standardize_columns=standardized_columns,
                        ensemble_size=ensemble_size,
                        ensemble_weighting=ensemble_weighting,
                        ensemble_temperature=ensemble_temperature,
                        ensemble_distinct_formulas=ensemble_distinct_formulas,
                        ensemble_rank_level=ensemble_rank_level,
                    )
                )
                final_spec = tuple(item["spec"] for item in final_components)
            except Exception:
                final_model = None
                final_spec = None
        else:
            for _, candidate_row in final_ranking.loc[final_ranking["eligible"]].iterrows():
                candidate_spec = candidate_row["spec"]
                try:
                    final_model, _ = _fit_cox_model(
                        work,
                        duration_col=duration_col,
                        event_col=event_col,
                        spec=candidate_spec,
                        covariates=covariates,
                        standardize_predictors=bool(standardize_predictors),
                        standardize_columns=standardized_columns,
                    )
                    final_spec = candidate_spec
                    break
                except Exception:
                    continue
        if verbose:
            print(
                "Completed final full-data ranking/final-model stage in {:.3f} seconds.".format(
                    time.time() - final_fit_started_at
                ),
                flush=True,
            )
        if final_model is None:
            warnings.warn(
                "Performance estimation succeeded, but no final model could be fitted on all patients.",
                RuntimeWarning,
                stacklevel=2,
            )

    return CoxEvaluationResult(
        c_index=c_index,
        confidence_interval=confidence_interval,
        confidence_level=confidence_level,
        bootstrap_standard_error=bootstrap_se,
        n_patients=n_patients,
        total_possible_pairs=total_possible_pairs,
        total_policy_eligible_pairs=total_policy_eligible_pairs,
        event_policy=event_policy,
        n_outer_pairs_requested=len(selected_outer_pairs),
        n_outer_pairs_fitted=int(fitted_mask.sum()),
        n_comparable_pairs=len(comparable_results),
        n_failed_outer_pairs=failed_outer_pairs,
        n_jobs_requested=n_jobs_requested,
        n_jobs_used=n_jobs_used,
        multiprocessing_start_method=(
            multiprocessing_start_method if n_jobs_used > 1 else "serial"
        ),
        inner_repeats=(inner_repeats if method == "kfold" else 1),
        inner_validation_method=method,
        inner_scoring=scoring,
        inner_selection_rule=selection_rule,
        inner_stability_penalty=float(inner_stability_penalty),
        standardize_predictors=bool(standardize_predictors),
        standardized_columns=standardized_columns,
        consensus_method=consensus_mode,
        consensus_top_models=consensus_top_models,
        consensus_top_features=consensus_top_features,
        consensus_threshold=consensus_threshold,
        consensus_max_formula_groups=consensus_max_formula_groups,
        consensus_pair_batch_size=consensus_pair_batch_size,
        consensus_filter_scoring=consensus_filter_scoring,
        consensus_feature_groups=normalized_consensus_groups,
        consensus_group_expressions=normalized_consensus_expressions,
        outer_pairs_reused=outer_pairs_reused,
        outer_pair_plan=pair_plan,
        outer_results=outer_results,
        selection_frequency=_selection_frequency(outer_results),
        bootstrap_c_indices=bootstrap_values,
        final_model=final_model,
        final_spec=final_spec,
        final_inner_ranking=final_ranking,
        prediction_aggregation=prediction_mode,
        ensemble_size=ensemble_size,
        ensemble_weighting=ensemble_weighting,
        ensemble_temperature=ensemble_temperature,
        ensemble_distinct_formulas=ensemble_distinct_formulas,
        ensemble_rank_level=ensemble_rank_level,
        bootstrap_n_jobs_used=bootstrap_n_jobs_used,
    )

# -----------------------------------------------------------------------------
# LaTeX / PGFPlots export helpers
# -----------------------------------------------------------------------------


def _latex_escape(value):
    """Escape a value for use as ordinary LaTeX text."""

    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    text = str(value)
    return "".join(replacements.get(char, char) for char in text)


def _validate_graphics_prefix(prefix):
    """Validate a filename prefix without requiring pathlib or extra packages."""

    prefix = str(prefix)
    if not prefix:
        raise ValueError("prefix cannot be empty.")
    if os.path.basename(prefix) != prefix:
        raise ValueError("prefix must be a filename prefix, not a path.")
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")
    if any(char not in allowed for char in prefix):
        raise ValueError(
            "prefix can contain only letters, numbers, underscores, and hyphens."
        )
    return prefix


def _write_utf8_text(path, text):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def _pgf_float(value):
    value = float(value)
    if not np.isfinite(value):
        raise ValueError("Non-finite value cannot be written as a PGFPlots coordinate.")
    return "{:.10g}".format(value)


def _bootstrap_histogram_data(result, histogram_bins):
    if int(histogram_bins) != histogram_bins or int(histogram_bins) < 5:
        raise ValueError("histogram_bins must be an integer >= 5.")
    histogram_bins = int(histogram_bins)

    values = np.asarray(result.bootstrap_c_indices, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return pd.DataFrame(columns=["bin_left", "bin_right", "bin_center", "proportion"])

    counts, edges = np.histogram(values, bins=histogram_bins, range=(0.0, 1.0))
    proportions = counts.astype(float) / float(counts.sum())
    centers = (edges[:-1] + edges[1:]) / 2.0
    return pd.DataFrame(
        {
            "bin_left": edges[:-1],
            "bin_right": edges[1:],
            "bin_center": centers,
            "proportion": proportions,
        }
    )


def _pair_outcome_data(result):
    rows = result.outer_results
    mask = rows["fit_succeeded"].astype(bool) & rows["comparable"].astype(bool)
    values = rows.loc[mask, "concordance_contribution"].to_numpy(dtype=float)

    discordant = int(np.isclose(values, 0.0).sum())
    tied = int(np.isclose(values, 0.5).sum())
    concordant = int(np.isclose(values, 1.0).sum())
    total = discordant + tied + concordant

    if total == 0:
        proportions = [np.nan, np.nan, np.nan]
    else:
        proportions = [
            discordant / float(total),
            tied / float(total),
            concordant / float(total),
        ]

    return pd.DataFrame(
        {
            "outcome": ["Discordant", "Tied", "Concordant"],
            "count": [discordant, tied, concordant],
            "proportion": proportions,
        }
    )


def _selection_plot_data(result, top_n_models):
    if int(top_n_models) != top_n_models or int(top_n_models) < 1:
        raise ValueError("top_n_models must be a positive integer.")
    top_n_models = int(top_n_models)

    frequency = result.selection_frequency.copy().head(top_n_models).reset_index(drop=True)
    if frequency.empty:
        frequency["model_id"] = []
        return frequency

    frequency.insert(
        0,
        "model_id",
        ["M{}".format(index + 1) for index in range(len(frequency))],
    )
    return frequency


def _make_bootstrap_tex(result, histogram_data):
    if histogram_data.empty:
        return "% No finite bootstrap C-index values were available.\n"

    ymax = float(histogram_data["proportion"].max())
    ymax = max(0.05, ymax * 1.15)
    width = float(histogram_data.iloc[0]["bin_right"] - histogram_data.iloc[0]["bin_left"])
    bar_width = width * 0.85

    coords = " ".join(
        "({},{})".format(_pgf_float(row.bin_center), _pgf_float(row.proportion))
        for row in histogram_data.itertuples(index=False)
    )

    lower = float(result.confidence_interval[0])
    upper = float(result.confidence_interval[1])
    c_index = float(result.c_index)

    return r"""%% Requires: \usepackage{pgfplots}
%% Recommended: \pgfplotsset{compat=1.15}
\begin{tikzpicture}
\begin{axis}[
    width=0.92\linewidth,
    height=0.55\linewidth,
    xmin=0, xmax=1,
    ymin=0, ymax=%(ymax)s,
    xlabel={Bootstrap C-index},
    ylabel={Bootstrap proportion},
    ybar,
    bar width=%(barwidth)s,
    grid=major,
    legend style={at={(0.5,-0.22)},anchor=north,legend columns=2},
]
\addplot coordinates {%(coords)s};
\addplot+[mark=none,dashed] coordinates {(%(c)s,0) (%(c)s,%(ymax)s)};
\addplot+[mark=none,dotted] coordinates {(%(lo)s,0) (%(lo)s,%(ymax)s)};
\addplot+[mark=none,dotted] coordinates {(%(hi)s,0) (%(hi)s,%(ymax)s)};
\legend{Bootstrap distribution,C-index,CI lower,CI upper}
\end{axis}
\end{tikzpicture}
""" % {
        "ymax": _pgf_float(ymax),
        "barwidth": _pgf_float(bar_width),
        "coords": coords,
        "c": _pgf_float(c_index),
        "lo": _pgf_float(lower),
        "hi": _pgf_float(upper),
    }


def _make_pair_outcomes_tex(pair_data):
    if pair_data["proportion"].isna().all():
        return "% No comparable outer pairs were available.\n"

    coords = " ".join(
        "({},{})".format(_latex_escape(row.outcome), _pgf_float(row.proportion))
        for row in pair_data.itertuples(index=False)
    )

    return r"""%% Requires: \usepackage{pgfplots}
%% Recommended: \pgfplotsset{compat=1.15}
\begin{tikzpicture}
\begin{axis}[
    width=0.82\linewidth,
    height=0.52\linewidth,
    ymin=0, ymax=1,
    ylabel={Proportion of comparable outer pairs},
    symbolic x coords={Discordant,Tied,Concordant},
    xtick=data,
    ybar,
    bar width=18pt,
    nodes near coords,
    nodes near coords align={vertical},
    grid=major,
]
\addplot coordinates {%(coords)s};
\end{axis}
\end{tikzpicture}
""" % {"coords": coords}


def _make_selection_frequency_tex(selection_data):
    if selection_data.empty:
        return "% No successful outer-model selections were available.\n"

    model_ids = selection_data["model_id"].tolist()
    coords = " ".join(
        "({},{})".format(row.model_id, _pgf_float(row.proportion))
        for row in selection_data.itertuples(index=False)
    )
    max_prop = max(0.05, float(selection_data["proportion"].max()) * 1.15)

    key_comments = []
    for row in selection_data.itertuples(index=False):
        key_comments.append(
            "%% %s: formula=%s; base_model=%s; knots=%s; piecewise_breaks=%s; "
            "penalizer=%s; l1_ratio=%s"
            % (
                row.model_id,
                str(row.formula),
                str(row.base_model),
                str(row.knot_number),
                str(row.piecewise_break_number),
                str(row.penalizer),
                str(row.l1_ratio),
            )
        )

    return r"""%% Requires: \usepackage{pgfplots}
%% Recommended: \pgfplotsset{compat=1.15}
%(key_comments)s
\begin{tikzpicture}
\begin{axis}[
    width=0.92\linewidth,
    height=0.55\linewidth,
    ymin=0, ymax=%(ymax)s,
    ylabel={Selection proportion},
    xlabel={Model specification ID},
    symbolic x coords={%(ids)s},
    xtick=data,
    x tick label style={rotate=45,anchor=east},
    ybar,
    bar width=12pt,
    grid=major,
]
\addplot coordinates {%(coords)s};
\end{axis}
\end{tikzpicture}
""" % {
        "key_comments": "\n".join(key_comments),
        "ymax": _pgf_float(max_prop),
        "ids": ",".join(model_ids),
        "coords": coords,
    }



def _require_matplotlib_for_png():
    if plt is None:
        raise ImportError(
            "PNG export requires matplotlib. For Python 3.6.8, install a "
            "Python-3.6-compatible release such as matplotlib==3.3.4."
        )


def _validate_png_dpi(png_dpi):
    if isinstance(png_dpi, bool) or int(png_dpi) != png_dpi or int(png_dpi) < 72:
        raise ValueError("png_dpi must be an integer >= 72.")
    return int(png_dpi)


def _save_current_figure_png(path, png_dpi):
    plt.tight_layout()
    plt.savefig(path, dpi=png_dpi, bbox_inches="tight")
    plt.close()


def _make_bootstrap_png(result, histogram_data, output_path, png_dpi):
    _require_matplotlib_for_png()
    if histogram_data.empty:
        return False

    centers = histogram_data["bin_center"].to_numpy(dtype=float)
    proportions = histogram_data["proportion"].to_numpy(dtype=float)
    width = float(
        histogram_data.iloc[0]["bin_right"] - histogram_data.iloc[0]["bin_left"]
    ) * 0.85

    plt.figure(figsize=(7.2, 4.4))
    plt.bar(centers, proportions, width=width, label="Bootstrap distribution")
    plt.axvline(float(result.c_index), linestyle="--", label="C-index")
    plt.axvline(float(result.confidence_interval[0]), linestyle=":", label="CI lower")
    plt.axvline(float(result.confidence_interval[1]), linestyle=":", label="CI upper")
    plt.xlim(0.0, 1.0)
    plt.xlabel("Bootstrap C-index")
    plt.ylabel("Bootstrap proportion")
    plt.grid(True)
    plt.legend()
    _save_current_figure_png(output_path, png_dpi)
    return True


def _make_pair_outcomes_png(pair_data, output_path, png_dpi):
    _require_matplotlib_for_png()
    if pair_data["proportion"].isna().all():
        return False

    labels = pair_data["outcome"].astype(str).tolist()
    proportions = pair_data["proportion"].to_numpy(dtype=float)
    x_values = np.arange(len(labels))

    plt.figure(figsize=(6.6, 4.4))
    bars = plt.bar(x_values, proportions)
    plt.xticks(x_values, labels)
    plt.ylim(0.0, 1.0)
    plt.ylabel("Proportion of comparable outer pairs")
    plt.grid(True, axis="y")
    for bar, value in zip(bars, proportions):
        if np.isfinite(value):
            plt.text(
                bar.get_x() + bar.get_width() / 2.0,
                value,
                "{:.3f}".format(value),
                ha="center",
                va="bottom",
            )
    _save_current_figure_png(output_path, png_dpi)
    return True


def _make_selection_frequency_png(selection_data, output_path, png_dpi):
    _require_matplotlib_for_png()
    if selection_data.empty:
        return False

    labels = selection_data["model_id"].astype(str).tolist()
    proportions = selection_data["proportion"].to_numpy(dtype=float)
    x_values = np.arange(len(labels))

    plt.figure(figsize=(7.2, 4.6))
    plt.bar(x_values, proportions)
    plt.xticks(x_values, labels, rotation=45, ha="right")
    plt.xlabel("Model specification ID")
    plt.ylabel("Selection proportion")
    plt.ylim(0.0, max(0.05, float(np.max(proportions)) * 1.15))
    plt.grid(True, axis="y")
    _save_current_figure_png(output_path, png_dpi)
    return True


def _make_cindex_comparison_png(results, labels, output_path, png_dpi):
    _require_matplotlib_for_png()

    c_values = np.asarray([result.c_index for result in results], dtype=float)
    lower = np.asarray([result.confidence_interval[0] for result in results], dtype=float)
    upper = np.asarray([result.confidence_interval[1] for result in results], dtype=float)
    positions = np.arange(len(results), dtype=float)

    plt.figure(figsize=(7.0, max(3.2, 0.7 * len(results) + 1.8)))
    for y_value, c_value, lo_value, hi_value in zip(
        positions, c_values, lower, upper
    ):
        plt.plot([lo_value, hi_value], [y_value, y_value])
        plt.plot([c_value], [y_value], marker="o")
    plt.yticks(positions, [str(label) for label in labels])
    plt.xlim(0.0, 1.0)
    plt.xlabel("Concordance index")
    plt.grid(True)
    plt.gca().invert_yaxis()
    _save_current_figure_png(output_path, png_dpi)
    return True


def export_performance_graphics(
    result,
    output_dir,
    prefix="cox_performance",
    top_n_models=10,
    histogram_bins=20,
    png_dpi=300,
):
    """Export performance figures as PGFPlots ``.tex`` and equivalent PNG.

    Three diagnostics are exported in both formats:

    1. bootstrap C-index distribution with point estimate and confidence limits;
    2. concordant/tied/discordant proportions among comparable outer test pairs;
    3. selection frequency of the most frequently chosen nested-CV model specs.

    CSV source tables and the scalar summary are also written. The returned
    dictionary contains paths for all generated files.
    """

    prefix = _validate_graphics_prefix(prefix)
    png_dpi = _validate_png_dpi(png_dpi)
    output_dir = os.path.abspath(str(output_dir))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    bootstrap_data = _bootstrap_histogram_data(result, histogram_bins)
    pair_data = _pair_outcome_data(result)
    selection_data = _selection_plot_data(result, top_n_models)

    paths = {
        "summary_csv": os.path.join(output_dir, prefix + "_summary.csv"),
        "bootstrap_csv": os.path.join(output_dir, prefix + "_bootstrap_histogram.csv"),
        "bootstrap_tex": os.path.join(output_dir, prefix + "_bootstrap_histogram.tex"),
        "bootstrap_png": os.path.join(output_dir, prefix + "_bootstrap_histogram.png"),
        "pair_outcomes_csv": os.path.join(output_dir, prefix + "_pair_outcomes.csv"),
        "pair_outcomes_tex": os.path.join(output_dir, prefix + "_pair_outcomes.tex"),
        "pair_outcomes_png": os.path.join(output_dir, prefix + "_pair_outcomes.png"),
        "selection_csv": os.path.join(output_dir, prefix + "_selection_frequency.csv"),
        "selection_tex": os.path.join(output_dir, prefix + "_selection_frequency.tex"),
        "selection_png": os.path.join(output_dir, prefix + "_selection_frequency.png"),
    }

    result.summary_frame().to_csv(paths["summary_csv"], index=True)
    bootstrap_data.to_csv(paths["bootstrap_csv"], index=False)
    pair_data.to_csv(paths["pair_outcomes_csv"], index=False)
    selection_data.to_csv(paths["selection_csv"], index=False)

    _write_utf8_text(paths["bootstrap_tex"], _make_bootstrap_tex(result, bootstrap_data))
    _write_utf8_text(paths["pair_outcomes_tex"], _make_pair_outcomes_tex(pair_data))
    _write_utf8_text(paths["selection_tex"], _make_selection_frequency_tex(selection_data))

    _make_bootstrap_png(result, bootstrap_data, paths["bootstrap_png"], png_dpi)
    _make_pair_outcomes_png(pair_data, paths["pair_outcomes_png"], png_dpi)
    _make_selection_frequency_png(selection_data, paths["selection_png"], png_dpi)

    return paths


def export_performance_graphics_latex(
    result,
    output_dir,
    prefix="cox_performance",
    top_n_models=10,
    histogram_bins=20,
    png_dpi=300,
):
    """Backward-compatible alias; now exports both PGFPlots and PNG files."""

    return export_performance_graphics(
        result=result,
        output_dir=output_dir,
        prefix=prefix,
        top_n_models=top_n_models,
        histogram_bins=histogram_bins,
        png_dpi=png_dpi,
    )


def export_cindex_comparison_graphics(
    results,
    labels,
    output_stem,
    png_dpi=300,
):
    """Export a forest-style C-index comparison as both ``.tex`` and PNG.

    ``output_stem`` is the path without an extension. A ``.tex`` PGFPlots
    snippet and a ``.png`` file are created next to each other.
    """

    if len(results) != len(labels):
        raise ValueError("results and labels must have the same length.")
    if len(results) == 0:
        raise ValueError("At least one result is required.")
    png_dpi = _validate_png_dpi(png_dpi)

    output_stem = os.path.abspath(str(output_stem))
    root, extension = os.path.splitext(output_stem)
    if extension.lower() in (".tex", ".png"):
        output_stem = root
    parent = os.path.dirname(output_stem)
    if parent and not os.path.isdir(parent):
        os.makedirs(parent)

    tex_path = output_stem + ".tex"
    png_path = output_stem + ".png"

    y_positions = list(range(len(results), 0, -1))
    yticklabels = [_latex_escape(label) for label in labels]

    plots = []
    for result, y_value in zip(results, y_positions):
        c_value = _pgf_float(result.c_index)
        lower = _pgf_float(result.confidence_interval[0])
        upper = _pgf_float(result.confidence_interval[1])
        y_text = str(y_value)
        plots.append(
            "\\addplot+[mark=none] coordinates {(%s,%s) (%s,%s)};\n"
            "\\addplot+[only marks,mark=*] coordinates {(%s,%s)};"
            % (lower, y_text, upper, y_text, c_value, y_text)
        )

    tex_height = max(0.48, min(1.0, 0.08 * len(results) + 0.18))
    tex = r"""%% Requires: \usepackage{pgfplots}
%% Recommended: \pgfplotsset{compat=1.15}
\begin{tikzpicture}
\begin{axis}[
    width=0.88\linewidth,
    height=%(height)s\linewidth,
    xmin=0, xmax=1,
    ymin=0.5, ymax=%(ymax)s,
    xlabel={Concordance index},
    ytick={%(yticks)s},
    yticklabels={%(yticklabels)s},
    grid=major,
]
%(plots)s
\end{axis}
\end{tikzpicture}
""" % {
        "height": _pgf_float(tex_height),
        "ymax": _pgf_float(len(results) + 0.5),
        "yticks": ",".join(str(v) for v in y_positions),
        "yticklabels": ",".join("{%s}" % value for value in yticklabels),
        "plots": "\n".join(plots),
    }
    _write_utf8_text(tex_path, tex)
    _make_cindex_comparison_png(results, labels, png_path, png_dpi)
    return {"tex": tex_path, "png": png_path}


def export_profile_cindex_summary_graphics(
    results,
    labels,
    output_stem,
    experiment_labels=None,
    profile_names=None,
    png_dpi=300,
):
    """Export a multi-profile C-index summary as PNG, PGFPlots, and CSV.

    This is a thin profile-aware wrapper around
    :func:`export_cindex_comparison_graphics`.  All supplied results should be
    based on the same outer-pair plan when the figure is intended for a paired
    comparison.  ``experiment_labels`` and ``profile_names`` are optional
    metadata columns for the summary CSV.
    """

    if len(results) != len(labels):
        raise ValueError("results and labels must have the same length.")
    if len(results) == 0:
        raise ValueError("At least one result is required.")
    if experiment_labels is not None and len(experiment_labels) != len(results):
        raise ValueError("experiment_labels must match results length.")
    if profile_names is not None and len(profile_names) != len(results):
        raise ValueError("profile_names must match results length.")

    paths = export_cindex_comparison_graphics(
        results=results,
        labels=labels,
        output_stem=output_stem,
        png_dpi=png_dpi,
    )
    rows = []
    for index, (result, label) in enumerate(zip(results, labels)):
        experiment = "" if experiment_labels is None else str(experiment_labels[index])
        profile = (
            str(getattr(result, "model_profile", ""))
            if profile_names is None
            else str(profile_names[index])
        )
        rows.append({
            "label": str(label),
            "experiment": experiment,
            "model_profile": profile,
            "c_index": float(result.c_index),
            "confidence_interval_lower": float(result.confidence_interval[0]),
            "confidence_interval_upper": float(result.confidence_interval[1]),
            "standard_error": float(getattr(result, "bootstrap_standard_error", np.nan)),
            "outer_pairs": int(len(result.outer_results)),
            "comparable_outer_pairs": int(np.sum(result.outer_results["comparable"].astype(bool))),
            "fit_successes": int(np.sum(result.outer_results["fit_succeeded"].astype(bool))),
        })
    csv_path = os.path.abspath(str(output_stem)) + "_summary.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    paths["summary_csv"] = csv_path
    return paths


def export_cindex_comparison_latex(results, labels, output_path, png_dpi=300):
    """Backward-compatible wrapper that also creates a sibling PNG file.

    The return value remains the ``.tex`` path for compatibility. If
    ``output_path`` ends in ``.tex``, the PNG is created using the same stem.
    """

    output_path = os.path.abspath(str(output_path))
    stem, extension = os.path.splitext(output_path)
    if extension.lower() != ".tex":
        stem = output_path
    paths = export_cindex_comparison_graphics(
        results=results,
        labels=labels,
        output_stem=stem,
        png_dpi=png_dpi,
    )
    return paths["tex"]


__all__ = [
    "MAX_MULTIPROCESSING_WORKERS",
    "CoxModelSpec",
    "CoxEvaluationResult",
    "nested_pair_cox_evaluation",
    "make_outer_pair_plan",
    "save_outer_pair_plan",
    "load_outer_pair_plan",
    "transform_for_fitted_cox_model",
    "save_evaluation_result_pickle",
    "load_evaluation_result_pickle",
    "export_performance_graphics",
    "export_performance_graphics_latex",
    "export_cindex_comparison_graphics",
    "export_cindex_comparison_latex",
    "export_profile_cindex_summary_graphics",
]
