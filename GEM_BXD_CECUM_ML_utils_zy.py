import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.preprocessing import LabelEncoder
from pathlib import Path


def set_random_seed(seed=0):
    np.random.seed(seed)
    tf.random.set_seed(seed)


def load_BXD_micro_data(data_path, tax_level="species"):

    omicsData = pd.read_csv(data_path, sep="\t", na_values="",
                            dtype={'taxonomy': str})
    omicsData = omicsData.dropna(subset=['taxonomy'])
    omicsData.iloc[:, 1:] = omicsData.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')

    tax_levels = ["phylum", "class", "order", "family", "genus", "species"]
    tax_patterns = ["\\|p__", "\\|c__", "\\|o__", "\\|f__", "\\|g__", "\\|s__", "\\|t__"]

    matching_index = tax_levels.index(tax_level)
    pattern = tax_patterns[matching_index]
    pattern_exclude = tax_patterns[matching_index + 1] if matching_index + 1 < len(tax_patterns) else None

    if pattern_exclude:
        omicsData = omicsData[
            omicsData['taxonomy'].str.contains(pattern) &
            ~omicsData['taxonomy'].str.contains(pattern_exclude)
            ]
    else:
        omicsData = omicsData[omicsData['taxonomy'].str.contains(pattern)]

    pattern_remove = pattern.lstrip("\\|")
    omicsData['taxonomy'] = omicsData['taxonomy'].str.replace(
        f'.*{pattern_remove}', pattern_remove, regex=True
    )

    omicsData = omicsData.transpose()
    omicsData.columns = omicsData.iloc[0]
    omicsData = omicsData.iloc[1:]

    return omicsData


def create_binary_groups(df, target_col, low_threshold, high_threshold,
                         group_col_name, low_label, high_label,
                         filter_col=None, filter_value=None):

    df_copy = df.copy()

    if filter_col and filter_value:
        mask = (df_copy[filter_col] == filter_value) & \
               ((df_copy[target_col] <= low_threshold) |
                (df_copy[target_col] >= high_threshold))
    else:
        mask = (df_copy[target_col] <= low_threshold) | \
               (df_copy[target_col] >= high_threshold)

    df_filtered = df_copy[mask].copy()

    df_filtered.loc[df_filtered[target_col] <= low_threshold, group_col_name] = low_label
    df_filtered.loc[df_filtered[target_col] >= high_threshold, group_col_name] = high_label

    return df_filtered


def prepare_classification_data(df, target_col, drop_cols):

    X = df.drop(drop_cols, axis=1)
    y = df[target_col]

    le = LabelEncoder()
    y_encoded = le.fit_transform(y)

    X = X.fillna(X.mean())

    return X, y_encoded, le


def calculate_tpr_fpr(y_true, y_pred, threshold):

    tp = np.sum((y_true == 1) & (y_pred >= threshold))
    fn = np.sum((y_true == 1) & (y_pred < threshold))
    fp = np.sum((y_true == 0) & (y_pred >= threshold))
    tn = np.sum((y_true == 0) & (y_pred < threshold))

    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0

    return tpr, fpr


def extract_shap_values_for_class(shap_values, model_name, class_idx=1):

    if isinstance(shap_values, list):
        return shap_values[class_idx]
    elif isinstance(shap_values, np.ndarray):
        if shap_values.ndim == 3:
            return shap_values[:, :, class_idx]
        return shap_values
    return shap_values


def save_feature_importance(feature_names, importance_scores,
                            output_path, model_name):

    importance_df = pd.DataFrame({
        'feature': feature_names,
        'importance': importance_scores
    }).sort_values('importance', ascending=False)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    importance_df.to_csv(output_path, sep='\t', index=False)
    print(f"[{model_name}] Feature importance saved to {output_path}")


def compute_mean_shap_importance(shap_values, feature_names):

    mean_abs_shap = np.abs(shap_values).mean(axis=0)

    importance_df = pd.DataFrame({
        'feature': feature_names,
        'mean_abs_shap': mean_abs_shap
    }).sort_values('mean_abs_shap', ascending=False)

    return importance_df