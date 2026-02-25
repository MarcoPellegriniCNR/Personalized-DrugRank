"""
Reproducible k-fold CV for lifelines CoxPHFitter + sklearn preprocessing
Input: a single DataFrame containing:
  - Response_duration (duration)
  - Response_coded    (event indicator)
All other columns are treated as covariates by default.

Dependencies: numpy, pandas, scikit-learn, lifelines
"""

import time
import traceback
import random
from typing import Dict, Any, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold

from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

#print("lifelines version:", lifelines.__version__)

############
import os
import sys

######################
# begin handling  command parameters
import argparse
activate_command_line = True
if activate_command_line :
    parser = argparse.ArgumentParser(description='Command line to run Cox regresion analysis.')

    parser.add_argument('experiment_labled',
                        action='store',
                        type = str,
                        choices=['STAD_fluorouracil_only',
                                 'STAD_fluorouracil_only_subtype',
                                 'CRC_irinotecan',
                                 'CRC_fluorouracil',
                                 'CRC_oxaliplatin',
                                 'CRC_leucovorin',
                                 'BRCA_A_paclitaxel',
                                 'BRCA_A_paclitaxel_subtype',
                                 'BRCA_B_cyclophosphamide',
                                 'BRCA_B_cyclophosphamide_subtype',
                                 'BRCA_C_fluorouracil',
                                 'BRCA_C_fluorouracil_subtype',
                                 'BRCA_D_docetaxel',
                                 'BRCA_D_docetaxel_subtype',
                                 'BRCA_E_doxorubicin',
                                 'BRCA_E_doxorubicin_subtype',
                                 ],
                        help='Defines the  experiment type: choose amnong 16 main experiments')

    parser.add_argument('local_event_coding_lable',
                        action='store',
                        type = str,
                        choices=[ 'PD1_all0', 'PDSD1_all0', 'CRPR1_all0', 'CR1_all0'],
                        help='Defines the  experiment coding for the mian events')

    parser.add_argument('baseline_method',
                        action='store',
                        type = str,
                        choices=[  "breslow",  "piecewise",  "spline"],
                        help='Defines the  function to be use in the Cox PH model')

    parser.add_argument( "--clinical_only_flag",
                    action="store_true",   # if the flag is present => True
                    default=False,         # optional; store_true already defaults to False
                    help="Enable only the clinical only parameters",)

    parser.add_argument( "--clinical_and_PDR_flag",
                    action="store_true",   # if the flag is present => True
                    default=False,         # optional; store_true already defaults to False
                    help="Enable only the clinical  plus PDR parameters",)

    args = parser.parse_args()

    print('Experiment {}'.format(args.experiment_labled))
    print('Event coding {}'.format(args.local_event_coding_lable))
    print('Function for Cox model {}'.format(args.baseline_method))
    print('Enable clinical only formuale {}'.format(args.clinical_only_flag))
    print('Enable clinical and PDR formuale {}'.format(args.clinical_and_PDR_flag))

    #assigments
    experiment_lable = args.experiment_labled
    local_event_coding_lable = args.local_event_coding_lable
    baseline_method = args.baseline_method
    clinical_only_flag = args.clinical_only_flag
    clinical_and_PDR_flag = args.clinical_and_PDR_flag

##############
# end  handling  command parameters   

          

##path = os.path.join("..", "Personalized-Drug-Ranking", "P_Drugmerge_Project_code_2025_02_27", "Personalized_DrugMerge" , "Cox_regression_analsysis")
path = ".."
print("Relative path exists", os.path.isdir(path))
abspath = os.path.abspath(path)
print("Absolute path exists", os.path.isdir(abspath))
sys.path.append(abspath)
print("Including in sys path to ", abspath)
from generic_util import printlist

####################
# select one

#STAD
#experiment_lable = 'STAD_fluorouracil_only'
#experiment_lable = 'STAD_fluorouracil_only_subtype'

#CRC
#experiment_lable = 'CRC_irinotecan'
#experiment_lable = 'CRC_fluorouracil'
#experiment_lable = 'CRC_oxaliplatin'
#experiment_lable = 'CRC_leucovorin'

#BRCA
#experiment_lable = 'BRCA_A_paclitaxel'
#experiment_lable = 'BRCA_A_paclitaxel_subtype'

#experiment_lable = 'BRCA_B_cyclophosphamide'
#experiment_lable = 'BRCA_B_cyclophosphamide_subtype'

#experiment_lable = 'BRCA_C_fluorouracil'
#experiment_lable = 'BRCA_C_fluorouracil_subtype'

#experiment_lable = 'BRCA_D_docetaxel'
#experiment_lable = 'BRCA_D_docetaxel_subtype'

#experiment_lable = 'BRCA_E_doxorubicin'
#experiment_lable = 'BRCA_E_doxorubicin_subtype'


f0 = 'C(American_Joint_Committee_on_Cancer_Metastasis_Stage_Code_Reduced)'
f1 = 'American_Joint_Committee_on_Cancer_Metastasis_Stage_Code_Reduced'


def remove_forbidden(in_list, forbidden) :
    out_list = []
    for formula in in_list :
        found = False
        for x in forbidden :
            if x in formula :
                found = True
                formula1  =formula.replace(x, "")
                formula1 = formula1.strip(" +")
                if formula1 != "" :
                    out_list.append(formula1)
                    print('substitution', formula, '====' ,  formula1)
                break
        if not found :        
            out_list.append(formula)
    printlist(out_list, 'Reduced List of fingerprints')
    return out_list
            
  


#STAD
if experiment_lable == 'STAD_fluorouracil_only' :
    from Cox_STAD_fluorouracil_only_analysis_iterated import  df, formulas_clinical_only, formulas_clinical_plus_PDR_indices, event_coding_lable
    
elif experiment_lable ==  'STAD_fluorouracil_only_subtype' :
    from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import  df, formulas_clinical_only, formulas_clinical_plus_PDR_indices, event_coding_lable
    
# CRC
elif experiment_lable ==  'CRC_irinotecan' :
    from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_CRC_analysis_irinotecan_patients  import  df_cleaned, event_coding_lable  
    df = df_cleaned
    
elif  experiment_lable ==  'CRC_fluorouracil' :
    from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_CRC_analysis_fluorouracil_patients import df_cleaned, event_coding_lable  # file da duplicare ed editare
    from clinical_data_analysis_CRC import sub_population_flag
    assert  sub_population_flag == 'fl_only',  "Error: Make sure the subpopluation flag  is aligned with the local  experiment lable CRC_fluorouraci and fl_only"
    df = df_cleaned
    formulas_clinical_only = remove_forbidden(formulas_clinical_only, [f0, f1])
    formulas_clinical_plus_PDR_indices = remove_forbidden(formulas_clinical_plus_PDR_indices, [f0, f1])
    
elif experiment_lable == 'CRC_oxaliplatin' :
    from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_CRC_analysis_fluorouracil_patients import df_cleaned,  event_coding_lable  # file da duplicare ed editare
    from clinical_data_analysis_CRC import sub_population_flag
    assert  sub_population_flag == 'fl_ox',  "Error: Make sure the subpopluation flag  is aligned with the local  experiment lable CRC_oxaliplatin and fl_ox "
    df = df_cleaned
    formulas_clinical_only = remove_forbidden(formulas_clinical_only, [f0, f1])
    formulas_clinical_plus_PDR_indices = remove_forbidden(formulas_clinical_plus_PDR_indices, [f0, f1])   

elif experiment_lable == 'CRC_leucovorin' :
    from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_CRC_analysis_fluorouracil_patients import df_cleaned, event_coding_lable  # file da duplicare ed editare
    from clinical_data_analysis_CRC import sub_population_flag
    assert  sub_population_flag == 'fl_leu',  "Error: Make sure the subpopluation flag  is aligned with the local  experiment lable CRC_leucovrin and fl_leu "
    df = df_cleaned
    formulas_clinical_only = remove_forbidden(formulas_clinical_only, [f0, f1])
    formulas_clinical_plus_PDR_indices = remove_forbidden(formulas_clinical_plus_PDR_indices, [f0, f1])   

# BRCA
elif experiment_lable ==  'BRCA_A_paclitaxel' :
    from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from   Cox_BRCA_paclitaxel_analysis import df_cleaned,  event_coding_lable  
    df = df_cleaned
    
elif experiment_lable == 'BRCA_A_paclitaxel_subtype' :
    from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from   Cox_BRCA_paclitaxel_analysis_with_subtypes import df_cleaned,  event_coding_lable 
    df = df_cleaned
    
elif experiment_lable ==  'BRCA_B_cyclophosphamide' :
    from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_BRCA_docetaxel_analysis import   df_cleaned,  event_coding_lable 
    df = df_cleaned
     
elif experiment_lable ==  'BRCA_B_cyclophosphamide_subtype':
    from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_BRCA_docetaxel_analysis_with_subtypes import  df_cleaned,  event_coding_lable 
    df = df_cleaned
     
elif experiment_lable == 'BRCA_C_fluorouracil':
    from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_BRCA_fluorouracil_analysis  import   df_cleaned,  event_coding_lable   
    df = df_cleaned
     
elif experiment_lable ==  'BRCA_C_fluorouracil_subtype' :
    from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_BRCA_fluorouracil_analysis_with_subtypes import   df_cleaned,  event_coding_lable
    df = df_cleaned

elif experiment_lable == 'BRCA_D_docetaxel' :
    from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_BRCA_docetaxel_D_2_analysis import  df_cleaned,  event_coding_lable 
    df = df_cleaned
     
elif experiment_lable == 'BRCA_D_docetaxel_subtype' :
    from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_BRCA_docetaxel_D_2_analysis_with_subtypes import  df_cleaned,  event_coding_lable 
    df = df_cleaned
     
elif experiment_lable == 'BRCA_E_doxorubicin':
    from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_BRCA_doxorubicin_analysis  import  df_cleaned,  event_coding_lable 
    df = df_cleaned
     
elif experiment_lable == 'BRCA_E_doxorubicin_subtype' :
    from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
    from Cox_BRCA_doxorubicin_analysis_with_subtypes  import df_cleaned,  event_coding_lable 
    df = df_cleaned

else :
    print("Wrong experiment_lable : ", experiment_lable )
    assert False


# checking consistency of the imported data with the local setting
#local_event_coding_lable = 'PD1_all0'
#local_event_coding_lable = 'PD2_CRPR1_all0'
#local_event_coding_lable ='CRPR1_all0'
#local_event_coding_lable ='CR1_all0'

assert local_event_coding_lable == event_coding_lable , "Error: Make sure the imported event_coding_lable is aligned with the local value"


def df_fingerprint(df: pd.DataFrame) -> int:
    # Stable hash across runs for the *content* (including index)
    return int(pd.util.hash_pandas_object(df, index=True).sum())


print("fingerprint:", df_fingerprint(df))
print(df.dtypes)
df = df.sort_values('Patient_ID').reset_index(drop=True)
print("fingerprint:", df_fingerprint(df))
print(df.dtypes)

###############


# -------------------------
# Deterministic seeding
# -------------------------
def set_global_seed(seed: int) -> None:
    np.random.seed(seed)
    random.seed(seed)
#    os.environ["PYTHONHASHSEED"] = str(seed)


# -------------------------
# Synthetic data generator (for demo) -> includes Response_* columns
# -------------------------
def generate_synthetic_df(
    n_samples: int = 800,
    n_features: int = 6,
    seed: int = 2026,
) -> pd.DataFrame:
    set_global_seed(seed)

    X = np.random.normal(0, 1, size=(n_samples, n_features))
    betas = np.arange(1, n_features + 1) * 0.5
    lin_pred = X.dot(betas)

    u = np.random.uniform(size=n_samples)
    baseline_scale = 10.0
    times = -np.log(u) * np.exp(-lin_pred) * baseline_scale

    censoring_threshold = np.quantile(times, 0.7)
    events = (times <= censoring_threshold).astype(int)
    durations = np.minimum(times, censoring_threshold)

    df = pd.DataFrame(X, columns=[f"x{i+1}" for i in range(n_features)])
    df["Response_duration"] = durations
    df["Response_coded"] = events
    return df


# -------------------------
# Helpers for piecewise baseline
# -------------------------
def _quantile_breakpoints(train_durations: np.ndarray, n_breaks: int = 3) -> List[float]:
    """
    Deterministically pick piecewise breakpoints using training durations only.
    Uses interior quantiles to avoid 0%/100% extremes.
    """
    if n_breaks < 1:
        raise ValueError("n_breaks must be >= 1")
    qs = np.linspace(0, 1, n_breaks + 2)[1:-1]  # interior quantiles
    bps = np.quantile(train_durations, qs).tolist()
    bps = [float(x) for x in bps if np.isfinite(x)]
    bps = sorted(set(bps))
    return bps


# -------------------------
# Main CV function (single DataFrame input)
# -------------------------
def k_fold_cox_cv_from_df(
    df: pd.DataFrame,
    duration_col: str = "Response_duration",
    event_col: str = "Response_coded",
    feature_cols: Optional[List[str]] = None,     # if None, use all columns except duration/event
    formula: Optional[str] = None,                # lifelines formula string (optional)
    n_splits: int = 5,
    seed: int = 42,
    penalizer: float = 0.0,
    fallback_penalizer: float = 1.0,
    # baseline switch
    baseline_estimation_method: str = "breslow",  # "breslow" | "spline" | "piecewise"
    n_baseline_knots: Optional[int] = None,       # used if method="spline"
    knots: Optional[List[float]] = None,          # used if method="spline"
    breakpoints: Optional[List[float]] = None,    # used if method="piecewise"
    n_piecewise_breaks: int = 3,                  # used if breakpoints is None & method="piecewise"
) -> pd.DataFrame:
    """
    K-fold CV for CoxPHFitter where df already contains duration/event columns.

    - Scales numeric covariates fold-wise (fit on train, apply to test)
    - Leaves non-numeric covariates untouched (useful for formula terms like C(cat))
    - Uses deterministic KFold shuffle with random_state=seed
    - Fold-wise exception handling; tries fallback penalizer on failure
    """

    if duration_col not in df.columns:
        raise ValueError(f"duration_col='{duration_col}' not found in df.columns")
    if event_col not in df.columns:
        raise ValueError(f"event_col='{event_col}' not found in df.columns")

    allowed = {"breslow", "spline", "piecewise"}
    if baseline_estimation_method not in allowed:
        raise ValueError(f"baseline_estimation_method must be one of {sorted(allowed)}")

#    set_global_seed(seed)  # do he call externally 

    # Decide covariate columns
    if feature_cols is None:
        feature_cols = [c for c in df.columns if c not in {duration_col, event_col}]
    else:
        missing = [c for c in feature_cols if c not in df.columns]
        if missing:
            raise ValueError(f"feature_cols contains columns not in df: {missing}")
        # Prevent accidental inclusion of response columns
        for bad in (duration_col, event_col):
            if bad in feature_cols:
                raise ValueError(f"feature_cols must NOT include '{bad}'")

    # Basic sanity on response types
    work_df = df.copy()
    work_df[duration_col] = pd.to_numeric(work_df[duration_col], errors="raise")
    # event can be bool/int; coerce to int 0/1
    work_df[event_col] = pd.to_numeric(work_df[event_col], errors="raise").astype(int)

#    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed) # unbalanced
# blancing code
    y = work_df[event_col].astype(int).to_numpy()
    min_events_per_fold = 2 
    n_events = int(y.sum())

    if n_events < min_events_per_fold * n_splits:
        raise ValueError(
            f"Not enough events to guarantee {min_events_per_fold} events per fold: "
            f"n_events={n_events}, n_splits={n_splits} "
            f"(need >= {min_events_per_fold * n_splits})."
        )
    kf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)

    results: List[Dict[str, Any]] = []

    # Identify which feature columns are numeric vs non-numeric (global, but scaling is fit fold-wise)
#    numeric_features = [c for c in feature_cols if pd.api.types.is_numeric_dtype(work_df[c])]  # prduces typing probems 
    numeric_features = [c for c in feature_cols if pd.api.types.is_float_dtype(work_df[c])]    
    non_numeric_features = [c for c in feature_cols if c not in numeric_features]

    for fold, (train_idx, test_idx) in enumerate(kf.split(work_df, y), start=1):

        assert y[test_idx].sum() >= min_events_per_fold
        
        record: Dict[str, Any] = {
            "fold": fold,
            "train_size": int(len(train_idx)),
            "test_size": int(len(test_idx)),
            "baseline_estimation_method": baseline_estimation_method,
            "c_index": np.nan,
            "fit_time_sec": np.nan,
            "fallback_used": False,
            "exception": None,
        }

        train_raw = work_df.iloc[train_idx].copy()
        test_raw = work_df.iloc[test_idx].copy()

        # Fold-wise scaling for numeric covariates only
        scaler = StandardScaler()
        train_scaled = train_raw.copy()
        test_scaled = test_raw.copy()

        if numeric_features:
            scaler.fit(train_raw[numeric_features].values)
            train_scaled.loc[:, numeric_features] = scaler.transform(train_raw[numeric_features].values)
            test_scaled.loc[:, numeric_features] = scaler.transform(test_raw[numeric_features].values)

        # (non-numeric features remain untouched)
        # Ensure we keep only requested covariates + response columns
        keep_cols = feature_cols + [duration_col, event_col]
        train_fold_df = train_scaled[keep_cols].copy()
        test_fold_df = test_scaled[keep_cols].copy()

        # Baseline kwargs, some depend on fold's training durations
        baseline_kwargs: Dict[str, Any] = {"baseline_estimation_method": baseline_estimation_method}

        if baseline_estimation_method == "spline":
            if knots is not None:
                baseline_kwargs["knots"] = knots
            else:
                baseline_kwargs["n_baseline_knots"] = int(n_baseline_knots) if n_baseline_knots is not None else 4

        elif baseline_estimation_method == "piecewise":
            if breakpoints is not None:
                baseline_kwargs["breakpoints"] = sorted(set(float(x) for x in breakpoints))
            else:
                baseline_kwargs["breakpoints"] = _quantile_breakpoints(
                    train_fold_df[duration_col].to_numpy(),
                    n_breaks=n_piecewise_breaks,
                )
            record["breakpoints_used"] = baseline_kwargs["breakpoints"]

        ######### begin internal function  called twice ######
        def _fit_and_score(current_penalizer: float) -> Tuple[float, float]:
            t0 = time.time()

            cph = CoxPHFitter(penalizer=current_penalizer, **baseline_kwargs)
            cph.fit(
                train_fold_df,
                duration_col=duration_col,
                event_col=event_col,
                formula=formula,          # can be None
                show_progress=False,
            )

            testing_attributes = False
            if testing_attributes :
                print("has predict_partial_hazard:", hasattr(cph, "predict_partial_hazard"))
                print("has predict_log_partial_hazard:", hasattr(cph, "predict_log_partial_hazard"))
                print("has score:", hasattr(cph, "score"))

            if hasattr(cph, "predict_partial_hazard") :
                # Predict on covariates only (safer than passing response columns)
                X_test = test_fold_df[feature_cols]
                risk = cph.predict_partial_hazard(X_test).astype(float)

                T_test = test_fold_df[duration_col].to_numpy()
                E_test = test_fold_df[event_col].to_numpy().astype(int)

                # Higher hazard => shorter survival, so negate for concordance_index
                c_idx = concordance_index(T_test, -risk.values, event_observed=E_test)

            if hasattr(cph, "score") :
                c_idx2 = cph.score(test_fold_df, scoring_method="concordance_index")

            if hasattr(cph, "predict_partial_hazard") and  hasattr(cph, "score") :
                if c_idx != c_idx2 :
                    print('Two c_idx values', c_idx, c_idx2)

            try:
               c_idx # does c_idx exist in the current namespace
            except NameError:
                c_idx  = c_idx2

            return float(c_idx), float(time.time() - t0)
        ######### end internal function  ######
        
        try:
            c_idx, t_sec = _fit_and_score(penalizer)
            record["c_index"] = c_idx
            record["fit_time_sec"] = t_sec

        except Exception as e:
            record["exception"] = f"{type(e).__name__}: {str(e)}"
            record["traceback"] = traceback.format_exc()

            # Fallback attempt with stronger penalizer
            try:
                record["fallback_used"] = True
                c_idx, t_sec = _fit_and_score(fallback_penalizer)
                record["c_index"] = c_idx
                record["fit_time_sec"] = t_sec
                record["exception"] += " | fallback_success"
            except Exception as e2:
                record["exception"] += f" || fallback_failed: {type(e2).__name__}: {str(e2)}"
                record["traceback_fallback"] = traceback.format_exc()

        results.append(record)

    return pd.DataFrame(results), baseline_kwargs


################################
# wrpper for one call

def one_call(df, formula, baseline_method, seed, k_fold, n_baseline_knot) :
    

    res, baseline_kwargs = k_fold_cox_cv_from_df(
        df=df,
        duration_col="Response_duration",
        event_col="Response_coded",
        feature_cols=None,              # all columns except Response_* are covariates
        formula=formula,
        n_splits=k_fold,
        seed=seed,
        penalizer=0.01,
        fallback_penalizer=1.0,
        baseline_estimation_method=baseline_method,
        # spline params (used only if baseline_method == "spline")
        n_baseline_knots=n_baseline_knot,
        knots=None,
        # piecewise params (used only if baseline_method == "piecewise")
        breakpoints=None,
        n_piecewise_breaks=3,
    )

    print("\nPer-fold results:")
    show_cols = ["fold", "train_size", "test_size", "baseline_estimation_method", "c_index", "fit_time_sec", "fallback_used", "exception"]
    if "breakpoints_used" in res.columns:
        show_cols.append("breakpoints_used")
#    print(res[show_cols])
    with pd.option_context(
        "display.max_rows", None,
        "display.max_columns", None,
        "display.width", None,
        "display.max_colwidth", None,   # show full text (e.g., exceptions)
    ):
        print(res[show_cols])


    cvals = res["c_index"].dropna().to_numpy()
    if cvals.size:
        mean = float(np.mean(cvals))
        std = float(np.std(cvals, ddof=1)) if cvals.size > 1 else 0.0
        max_c = float(max(cvals)) 
        min_c = float(min(cvals))
        print(f"\nAggregate C-index over {cvals.size} successful folds: mean={mean:.4f}, std={std:.4f} , min={min_c:.4f}, max={max_c:.4f}")
        return (mean, std, min_c, max_c) , baseline_kwargs
    else:
        print("\nNo successful folds to aggregate.")
        return None , baseline_kwargs



##############
# proces one formula

def process_one_formula(df, formula, baseline_method, n_baseline_knot) :

    repeats_number = NUM_REPEATS 
    rng = np.random.default_rng(SEED)
    seeds = rng.integers(0, 2**32 - 1, size=repeats_number)

    y = df[event_col].astype(int).to_numpy()     
    n_events = int(y.sum())
    max_k_fold = n_events // 2
    print('max number of folds ', max_k_fold)

    accum_results ={}
    
    for k_fold in range(2, max_k_fold +1) :
# beign  k_fold loop        
        repeat_result = []
        
        for iteration in range(repeats_number) :
            res, baseline_kwargs = one_call(df, formula, baseline_method, seed=seeds[iteration] , k_fold=k_fold,  n_baseline_knot=n_baseline_knot)
            if res != None :
                repeat_result.append(res)

        for x in  range(len(repeat_result)):
            print(x, repeat_result[x])

        print("===== Results for formula {}  k-fold {} , n_baseline_knot {} ====".format(formula, k_fold, n_baseline_knot))
        num_res = len(repeat_result)
        if num_res > 0 :
            global_mean = float(np.mean([x[0] for x in repeat_result]))
            global_std =  float(np.mean([x[1] for x in repeat_result]))
            print('global_mean = {}, global std = {} over {} successful repetitions '.format(global_mean,global_std , num_res))
            accum_results[(formula, k_fold, n_baseline_knot)] = [global_mean, global_std, num_res]
        else :
            print('No global mean')
# end  k_fold loop

    work_df = df.copy()
#    cph = CoxPHFitter(baseline_method, penalizer=current_penalizer)
    cph = CoxPHFitter(penalizer=current_penalizer, **baseline_kwargs)
    try :
        cph.fit(
                work_df,
                duration_col="Response_duration",
                event_col="Response_coded",
                formula=formula,          # can be None
                show_progress=False,
#                baseline_estimation_method  =  baseline_method,  # "breslow" | "spline" | "piecewise"
                 # other aprameters to be used in other options besides breslow"
            )
        print("===== Results for formula {}  model and LL p val ====".format(formula))
        cph.print_summary()
        lrt_result = cph.log_likelihood_ratio_test()
        p_value = lrt_result.p_value
        log2_p_value = -np.log2(p_value)
        p_value_f = float(p_value)
        print('recomputed pval and  - log 2 pval ', p_value,  log2_p_value)

    except Exception as e:
        print( f"{type(e).__name__}: {str(e)}")
        print( traceback.format_exc())
        p_value_f = None

    for key in accum_results.keys() :
        accum_results[key].append(p_value_f)
        
    return accum_results

# -------------------------
# Example run
# -------------------------
if __name__ == "__main__":
# Choose a seed
#    SEED = 2026
#    SEED = 7856432
# Choose a number of repeats
    set_global_seed(SEED)
#    NUM_REPEATS = 10
#    NUM_REPEATS = 20

#    df = generate_synthetic_df(n_samples=800, n_features=6, seed=SEED)

    # Pick baseline mode:
    #   "breslow"   -> classic Cox (non-param baseline)
    #   "spline"    -> parametric-ish baseline via splines
    #   "piecewise" -> parametric-ish baseline via piecewise hazard
    
#    baseline_method = "breslow"
#    baseline_method = "piecewise"
#    baseline_method = "spline"
    
    current_penalizer = 0.01
    event_col="Response_coded"
    
    PDR_start_index = 0
#    PDR_start_index = 31
    clinical_start_index = 0

    # Optional: explicit formula (otherwise uses all covariates in feature_cols)
    # formula = "x1 + x2 + x3 + x4"
#    formula = None
#    formula = formulas_clinical_only[0]

    if baseline_method ==  "spline" :
        n_baseline_knots_list = [2,3,4,5,6,7,8]
    else :
        n_baseline_knots_list = [0]


#    clinical_only_flag = True
#    clinical_and_PDR_flag =False

    processed_formulas = set()

# fixing long names in columns, neded for unix runs
    df = df.rename(columns={
    "Neoplasm_Disease_Lymph_Node_Stage_American_Joint_Committee_on_Cancer_Code_Encoded": "Neoplasm_Disease_Lymph_Node_Stage_AJCCC_Encoded"})
#    df = df.rename(columns=lambda c: c.strip().strip('"').strip("'").replace("\r","").replace("\n",""))


    print("=============== Main parametrs setting ===============")
    print("Experiment_lable {}, event_coding {}, baseline_method {}, current_penalizer {}, clinical_only_flag {}, clinical_and_PDR_flag {}, ".format(experiment_lable,  event_coding_lable, baseline_method, current_penalizer, clinical_only_flag, clinical_and_PDR_flag))

    if clinical_only_flag :
        print("=============== Begin analysis of  formulas_clinical_only ===============")
        formulas_clinical_only_results = []
        for formula in formulas_clinical_only[clinical_start_index:] :
            if formula not  in  processed_formulas :
                for n_baseline_knots  in n_baseline_knots_list :
                    accum_results = process_one_formula(df, formula, baseline_method, n_baseline_knots)
                    formulas_clinical_only_results.extend(accum_results.items())
                processed_formulas.add(formula)                          
        sorted_formulas_clinical_only_results = sorted(formulas_clinical_only_results, key=lambda tup: tup[1][0], reverse=True) # sort by Concordance Index
        printlist( sorted_formulas_clinical_only_results, ' sorted_formulas_clinical_only_results - sort by Concordance Index')

    if clinical_and_PDR_flag :
        formulas_clinical_plus_PDR_results = []
        print("=============== Begin analysis of formulas_clinical_plus_PDR_indices  ===============")
        for formula in formulas_clinical_plus_PDR_indices[PDR_start_index:] :
            if formula not  in  processed_formulas :
                for n_baseline_knots  in n_baseline_knots_list :
                    accum_results = process_one_formula(df, formula, baseline_method, n_baseline_knots)
                    formulas_clinical_plus_PDR_results.extend(accum_results.items())                
                processed_formulas.add(formula)
               
        sorted_formulas_clinical_plus_PDR_results = sorted(formulas_clinical_plus_PDR_results, key=lambda tup: tup[1][0], reverse=True)  # sort by Concordance Index
        printlist( sorted_formulas_clinical_plus_PDR_results, ' sorted_formulas_clinical_plus_PDR_results - sort by Concordance Index')

    print ('***Analysed {} formulae'.format(len(processed_formulas)))
            
        
