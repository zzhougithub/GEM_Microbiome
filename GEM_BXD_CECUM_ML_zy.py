import pandas as pd
import numpy as np
import tensorflow as tf
import datetime
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.metrics import roc_auc_score, auc
import shap


meta_name = "/Users/path/0_info_summary.txt"
metaData = pd.read_csv(meta_name, sep="\t", index_col=0)

data_name = "/Users/path/0_metaDNA_metaphlan4_abundance_all_zy.txt"
tax_level = "species"
omicsData = pd.read_csv(data_name, sep="\t", na_values="", dtype={'taxonomy': str})
omicsData = omicsData.dropna(subset=['taxonomy'])
omicsData.iloc[:, 1:] = omicsData.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')

tax = ["phylum", "class", "order", "family", "genus", "species"]
tax_pattern = ["\\|p__", "\\|c__", "\\|o__", "\\|f__", "\\|g__", "\\|s__", "\\|t__"]
matching_index = tax.index(tax_level)
pattern = tax_pattern[matching_index]
pattern_exclude = tax_pattern[matching_index + 1] if matching_index + 1 < len(tax_pattern) else None

if pattern_exclude:
    omicsData = omicsData[omicsData['taxonomy'].str.contains(pattern) & ~omicsData['taxonomy'].str.contains(pattern_exclude)]
else:
    omicsData = omicsData[omicsData['taxonomy'].str.contains(pattern)]
pattern_remove = pattern.lstrip("\\|")
omicsData['taxonomy'] = omicsData['taxonomy'].str.replace(f'.*{pattern_remove}', pattern_remove, regex=True)
omicsData = omicsData.transpose()
omicsData.columns = omicsData.iloc[0]
omicsData = omicsData.iloc[1:]

df1 = pd.merge(omicsData, metaData[['Diet', 'Age']], left_index=True, right_index=True, how='left')

Age1, Age2 = 350, 450

df2 = df1[(df1['Diet'] == 'CD') & ((df1['Age'] <= Age1) | (df1['Age'] >= Age2))].copy()
df2.loc[df2['Age'] <= Age1, 'AgeGroup'] = 'Young'
df2.loc[df2['Age'] >= Age2, 'AgeGroup'] = 'Old'
df2 = df2[['Diet', 'Age', 'AgeGroup'] + [col for col in df2.columns if col not in ['Diet', 'Age', 'AgeGroup']]].copy()

def calculate_tpr_fpr(y_test, y_pred, threshold):
    tp = np.sum((y_test == 1) & (y_pred >= threshold))
    fn = np.sum((y_test == 1) & (y_pred < threshold))
    tpr = tp / (tp + fn)

    fp = np.sum((y_test == 0) & (y_pred >= threshold))
    tn = np.sum((y_test == 0) & (y_pred < threshold))
    fpr = fp / (fp + tn)
    return tpr, fpr

random_state_value = 0
def set_seed(seed = random_state_value):
    np.random.seed(seed)
    tf.random.set_seed(seed)
set_seed()

def age_classification():

    data = df2

    X = data.drop(['Diet', 'Age', 'AgeGroup'], axis = 1)

    y = data['AgeGroup']

    le = LabelEncoder()
    y_encoded = le.fit_transform(y)

    means = X.mean()
    X = X.fillna(means)

    kfold = KFold(n_splits=5, shuffle=True, random_state=random_state_value)

    models = [
        ('Random Forest', RandomForestClassifier(random_state=random_state_value)),
        ('XGBoost', XGBClassifier(random_state=random_state_value))
    ]

    current_date = datetime.datetime.now().strftime("%Y%m%d")

    auc_results = {}
    roc_info = {'Random Forest': [], 'XGBoost': []}
    auc_by_num_features = {'Random Forest': {'mean_auc': [], 'std_auc': []},
                           'XGBoost': {'mean_auc': [], 'std_auc': []}}
    feature_weights = {'Random Forest': [], 'XGBoost': []}

    num_thresholds = 500

    def save_feature_weights(model_name, feature_names, avg_weights):
        feature_importances = pd.DataFrame({
            'feature': feature_names,
            'weight': avg_weights
        })
        feature_importances = feature_importances.sort_values(by='feature', ascending=False)
        filename = f'1_{model_name}_feature_importances.txt'
        feature_importances.to_csv(filename, sep='\t', index=False)

    for name, model in models:
        print('')
        print(f'{name}:')

        auc_scores = []
        weights = []
        shap_values_folds = []
        X_test_folds = []

        for fold, (train_index, test_index) in enumerate(kfold.split(X), ):
            #print(f"Fold {fold + 1}:")

            X_train, X_test = X.iloc[train_index], X.iloc[test_index]
            y_train, y_test = y_encoded[train_index], y_encoded[test_index]

            model.fit(X_train, y_train)
            y_pred = model.predict_proba(X_test)[:, 1]
            auc_scores.append(roc_auc_score(y_test, y_pred))

            thresholds = np.linspace(0, 1, num_thresholds)
            tpr = np.zeros(num_thresholds)
            fpr = np.zeros(num_thresholds)
            for i, threshold in enumerate(thresholds):
                tpr[i], fpr[i] = calculate_tpr_fpr(y_test, y_pred, threshold)

            roc_auc = auc(fpr, tpr)
            roc_info[name].append((fpr, tpr, roc_auc))

            w_model = model.feature_importances_
            weights.append(w_model)

            explainer = shap.TreeExplainer(model, )
            x_test = pd.DataFrame(X_test, )
            X_test_folds.append(x_test)
            shap_values = explainer.shap_values(x_test, )
            if name == 'Random Forest':
                shap_values_class = shap_values[:, :, 1]  # always 3D
            elif name == 'XGBoost':
                if isinstance(shap_values, list):
                    shap_values_class = shap_values[1]  # list of (samples, features)
                elif isinstance(shap_values, np.ndarray) and shap_values.ndim == 3:
                    shap_values_class = shap_values[:, :, 1]  # 3D array
                else:
                    shap_values_class = shap_values
            shap_values_folds.append(shap_values_class)


        mean_auc = np.mean(auc_scores)
        std_auc = np.std(auc_scores)
        auc_results[name] = mean_auc
        print(f'{name} - Mean AUC : {mean_auc:.4f}, Standard Deviation: {std_auc:.4f}')

        avg_weights = np.mean(weights, axis=0)
        save_feature_weights(name, X.columns, avg_weights)
        feature_weights[name] = avg_weights

        shap_fig = f"../images/1_age_Y&O_top20_features_shap_seed{random_state_value}_{name}_zy_{current_date}.pdf"
        X_shap_total = pd.concat(X_test_folds, axis=0)
        shap_values_total = np.vstack(shap_values_folds)
        shap.summary_plot(shap_values_total, X_shap_total, max_display=20, show=False)
        #plt.savefig(shap_fig, format='pdf')
        plt.show()

        mean_abs_shap = np.abs(shap_values_total).mean(axis=0)
        shap_importance_df = pd.DataFrame({
            'feature': X_shap_total.columns,
            'mean_abs_shap': mean_abs_shap
        })
        shap_importance_df = shap_importance_df.sort_values(by='mean_abs_shap', ascending=False)
        shap_data = f"../1_age_Y&O_top20_features_shap_seed{random_state_value}_{name}_zy_{current_date}.txt"

        top_features = shap_importance_df['feature'].head(50).tolist()
        X_selected = X.loc[:, top_features]

        for num_features in range(5, 51, 5):
            X_selected_subset = X_selected.iloc[:, :num_features]
            auc_scores = []
            for train_index, test_index in kfold.split(X_selected_subset):
                X_train, X_test = X_selected_subset.iloc[train_index], X_selected_subset.iloc[test_index]
                y_train, y_test = y_encoded[train_index], y_encoded[test_index]
                model.fit(X_train, y_train)
                y_pred = model.predict_proba(X_test)[:, 1]
                auc_scores.append(roc_auc_score(y_test, y_pred))
            mean_auc = np.mean(auc_scores)
            std_auc = np.std(auc_scores)
            auc_by_num_features[name]['mean_auc'].append(mean_auc)
            auc_by_num_features[name]['std_auc'].append(std_auc)
            print(f'{name} - {num_features} Features - Mean AUC: {mean_auc:.4f}, Standard Deviation: {std_auc:.4f}')

        Optimal_features = 50
        X_features = X.loc[:, top_features[:Optimal_features]]


    file_name_2 = f"../images/1_age_Y&O_mean_auc_seed{random_state_value}_features_select_errorbar_zy_{current_date}.pdf"
    file_name_3 = f"../images/1_age_Y&O_mean_auc_seed{random_state_value}_features_select_SHAP_errorbar_zy_{current_date}.pdf"
    colors = {'Random Forest': 'skyblue',
              'XGBoost': 'salmon'}
    plt.figure(figsize=(7.2, 5.4))
    width = 2  
    positions = [-0.5, 0.5]

    adjusted_positions = [i * (width + 0.1) for i in positions]
    for i, (model_name, auc_data) in enumerate(auc_by_num_features.items()):
        color = colors[model_name]
        auc_list = auc_data['mean_auc']
        std_list = auc_data['std_auc']
        pos_offset = adjusted_positions[i]

        plt.bar(np.arange(5, 51, 5) + pos_offset, auc_list, width=width, label=f'{model_name}',
                color=color, alpha=0.75, zorder=3)
        plt.errorbar(np.arange(5, 51, 5) + pos_offset, auc_list, yerr=std_list, fmt='none', color='darkgray',
                     capsize=2)
    plt.legend()

    plt.grid(axis='y', linestyle='--', zorder=2)
    plt.title('AUC Prediction of Age for MG Data Species by Number of Features')
    plt.xlabel('Number of Selected Species', fontsize = 14)
    plt.ylabel('Mean AUC', fontsize = 14)
    plt.legend(loc='upper right', ncol=len(auc_by_num_features))
    plt.xticks(np.arange(5, 51, 5), fontsize = 14)

    plt.ylim(0.5, 1)
    plt.yticks([i * 0.05 + 0.5 for i in range(11)], fontsize = 14)

    for spine in plt.gca().spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    plt.show()


    for name, weights in feature_weights.items():
        file_name_4 = f"../images/1_age_Y&O_top20_features_mean_weight_seed{random_state_value}_{name}_zy_{current_date}.pdf"

        feature_importances = pd.Series(weights, index=X.columns)

        top_features = feature_importances.nlargest(20).sort_values(ascending=True)

        plt.figure(figsize=(10, 8))
        top_features.plot(kind='barh')
        plt.title(f'Top 20 Feature Importances for {name}')
        plt.xlabel('Importance')
        plt.ylabel('Features')
        #plt.savefig(file_name_4, format='pdf')
        plt.show()


if __name__ == '__main__':
    age_classification()

