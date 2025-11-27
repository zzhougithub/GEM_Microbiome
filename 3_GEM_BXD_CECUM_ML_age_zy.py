import os
os.environ['PYTHONHASHSEED'] = '0'
os.environ['TF_DETERMINISTIC_OPS'] = '1'
os.environ['TF_CUDNN_DETERMINISTIC'] = '1'

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from sklearn.metrics import roc_auc_score, auc
import shap

from GEM_BXD_CECUM_ML_utils_zy import (
    set_random_seed, load_BXD_micro_data, create_binary_groups,
    prepare_classification_data, calculate_tpr_fpr,
    extract_shap_values_for_class, save_feature_importance,
    compute_mean_shap_importance
)


class AgeClassificationPipeline:

    def __init__(self, data_dir, output_dir, random_state=0):
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.random_state = random_state

        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / 'images').mkdir(exist_ok=True)

        set_random_seed(random_state)

    def load_data(self):
        print("Loading data...")

        meta_path = self.data_dir / "S1.2_metadata_summary_zy.txt"
        data_path = self.data_dir / "S4_metaDNA_metaphlan4_abundance_zy.txt"

        metaData = pd.read_csv(meta_path, sep="\t", index_col=0)
        omicsData = load_BXD_micro_data(data_path, tax_level="species")

        df = pd.merge(omicsData, metaData[['Diet', 'Age']],
                      left_index=True, right_index=True, how='left')

        self.df = create_binary_groups(
            df, target_col='Age',
            low_threshold=350, high_threshold=450,
            group_col_name='AgeGroup',
            low_label='Young', high_label='Old',
            filter_col='Diet', filter_value='CD'
        )

        print(f"Data loaded: {self.df.shape}")
        print(f"Age group distribution:\n{self.df['AgeGroup'].value_counts()}")

    def step1_compare_models(self, n_splits=5):
        print("\n" + "=" * 50)
        print("STEP 1: Model Comparison")
        print("=" * 50)

        X, y_encoded, le = prepare_classification_data(
            self.df, target_col='AgeGroup',
            drop_cols=['Diet', 'Age', 'AgeGroup']
        )

        def build_dnn():
            import random
            import tensorflow as tf
            random.seed(self.random_state)
            np.random.seed(self.random_state)
            tf.random.set_seed(self.random_state)

            model = Sequential([
                Dense(64, input_dim=X.shape[1], activation='relu',
                      kernel_initializer=tf.keras.initializers.GlorotUniform(seed=self.random_state)),
                Dense(32, activation='relu',
                      kernel_initializer=tf.keras.initializers.GlorotUniform(seed=self.random_state)),
                Dense(1, activation='sigmoid',
                      kernel_initializer=tf.keras.initializers.GlorotUniform(seed=self.random_state))
            ])
            model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.01),
                          loss='binary_crossentropy',
                          metrics=['accuracy'])
            return model

        models = [
            ('Random Forest', RandomForestClassifier(random_state=self.random_state)),
            ('XGBoost', XGBClassifier(random_state=self.random_state)),
            ('Logistic Regression', LogisticRegression(random_state=self.random_state)),
            ('DNN', build_dnn)
        ]



        colors = {
            'Random Forest': '#1F77B4',
            'XGBoost': '#FF7F0E',
            'Logistic Regression': '#2CA02C',
            'DNN': '#D62728'
        }

        kfold = KFold(n_splits=n_splits, shuffle=True, random_state=self.random_state)
        roc_curves = {name: [] for name, _ in models}
        num_thresholds = 500

        for name, model_fn in models:
            print(f"\n{name}:")
            auc_scores = []

            for fold, (train_idx, test_idx) in enumerate(kfold.split(X), 1):
                X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
                y_train, y_test = y_encoded[train_idx], y_encoded[test_idx]

                if name == 'DNN':
                    model = model_fn()
                    model.fit(X_train, y_train, epochs=10, batch_size=32, verbose=0)
                    y_pred = model.predict(X_test).ravel()
                else:
                    model = model_fn if not callable(model_fn) else model_fn()
                    model.fit(X_train, y_train)
                    y_pred = model.predict_proba(X_test)[:, 1]

                thresholds = np.linspace(0, 1, num_thresholds)
                tpr = np.zeros(num_thresholds)
                fpr = np.zeros(num_thresholds)
                for i, thresh in enumerate(thresholds):
                    tpr[i], fpr[i] = calculate_tpr_fpr(y_test, y_pred, thresh)

                auc_val = auc(fpr, tpr)
                roc_curves[name].append((fpr, tpr, auc_val))
                auc_scores.append(auc_val)

                print(f"  Fold {fold}: AUC = {auc_val:.4f}")

            print(f"  Mean AUC: {np.mean(auc_scores):.4f} ± {np.std(auc_scores):.4f}")

        self._plot_roc_curves(roc_curves, colors)

    def step2_feature_selection_and_shap(self, n_splits=5, n_top_features=50):
        print("\n" + "=" * 50)
        print("STEP 2: Feature Selection & SHAP Analysis")
        print("=" * 50)

        X, y_encoded, le = prepare_classification_data(
            self.df, target_col='AgeGroup',
            drop_cols=['Diet', 'Age', 'AgeGroup']
        )

        models = {
            'Random Forest': RandomForestClassifier(random_state=self.random_state),
            'XGBoost': XGBClassifier(random_state=self.random_state)
        }

        kfold = KFold(n_splits=n_splits, shuffle=True, random_state=self.random_state)

        all_models_auc_results = {}

        for name, model in models.items():
            print(f"\n{name}:")

            shap_values_folds = []
            X_test_folds = []
            feature_importances = []

            for fold, (train_idx, test_idx) in enumerate(kfold.split(X), 1):
                X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
                y_train, y_test = y_encoded[train_idx], y_encoded[test_idx]

                model.fit(X_train, y_train)

                feature_importances.append(model.feature_importances_)

                explainer = shap.TreeExplainer(model)
                shap_vals = explainer.shap_values(X_test)
                shap_vals_class = extract_shap_values_for_class(shap_vals, name)

                shap_values_folds.append(shap_vals_class)
                X_test_folds.append(X_test)

            avg_importance = np.mean(feature_importances, axis=0)
            save_feature_importance(
                X.columns, avg_importance,
                self.output_dir / f"{name}_feature_importances.txt",
                name
            )

            X_shap_total = pd.concat(X_test_folds, axis=0)
            shap_total = np.vstack(shap_values_folds)

            shap.summary_plot(shap_total, X_shap_total, max_display=20, show=False)
            plt.tight_layout()
            plt.savefig(self.output_dir / 'images' / f'{name}_shap_summary.pdf')
            plt.show()

            shap_importance = compute_mean_shap_importance(shap_total, X_shap_total.columns)
            shap_importance.to_csv(
                self.output_dir / f"{name}_shap_importance.txt",
                sep='\t', index=False
            )

            auc_results = self._evaluate_feature_selection(
                X, y_encoded, shap_importance, name, kfold, n_top_features
            )

            all_models_auc_results[name] = auc_results

        self._plot_feature_selection_results(all_models_auc_results, n_top_features)

    def _evaluate_feature_selection(self, X, y, importance_df, model_name,
                                    kfold, max_features=50):

        top_features = importance_df['feature'].head(max_features).tolist()

        model = RandomForestClassifier(random_state=self.random_state) if model_name == 'Random Forest' \
            else XGBClassifier(random_state=self.random_state)

        auc_results = {'mean': [], 'std': []}

        for n_feat in range(5, max_features + 1, 5):
            X_selected = X.loc[:, top_features[:n_feat]]
            auc_scores = []

            for train_idx, test_idx in kfold.split(X_selected):
                X_train, X_test = X_selected.iloc[train_idx], X_selected.iloc[test_idx]
                y_train, y_test = y[train_idx], y[test_idx]

                model.fit(X_train, y_train)
                y_pred = model.predict_proba(X_test)[:, 1]
                auc_scores.append(roc_auc_score(y_test, y_pred))

            auc_results['mean'].append(np.mean(auc_scores))
            auc_results['std'].append(np.std(auc_scores))

            print(f"  {n_feat} features: AUC = {np.mean(auc_scores):.4f} ± {np.std(auc_scores):.4f}")

        return auc_results

    def _plot_roc_curves(self, roc_curves, colors):
        plt.figure(figsize=(8, 6), facecolor='white')
        ax = plt.gca()
        ax.set_facecolor('white')

        for model_name, curves in roc_curves.items():
            fpr_avg = np.mean([c[0] for c in curves], axis=0)
            tpr_avg = np.mean([c[1] for c in curves], axis=0)
            auc_avg = np.mean([c[2] for c in curves])

            plt.plot(fpr_avg, tpr_avg,
                     label=f'{model_name} (AUC = {auc_avg:.2f})',
                     color=colors[model_name], linewidth=2)

            tpr_std = np.std([c[1] for c in curves], axis=0)
            plt.fill_between(fpr_avg,
                             np.maximum(tpr_avg - tpr_std, 0),
                             np.minimum(tpr_avg + tpr_std, 1),
                             color=colors[model_name], alpha=0.1)

        plt.plot([0, 1], [0, 1], 'k--', linewidth=1)
        plt.xlabel('False Positive Rate', fontsize=12)
        plt.ylabel('True Positive Rate', fontsize=12)
        plt.title('ROC Curve for Age Classification', fontsize=14)
        plt.legend(loc='lower right')
        plt.grid(False)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'images' / 'roc_comparison.pdf')
        plt.show()

    def _plot_feature_selection_results(self, all_models_auc_results, max_features=50):
        print("\nPlotting feature selection results...")

        colors = {
            'Random Forest': 'skyblue',
            'XGBoost': 'salmon'
        }

        plt.figure(figsize=(7.2, 5.4), facecolor='white')
        ax = plt.gca()
        ax.set_facecolor('white')

        width = 2
        positions = [-0.5, 0.5]
        adjusted_positions = [i * (width + 0.1) for i in positions]

        for i, (model_name, auc_data) in enumerate(all_models_auc_results.items()):
            color = colors[model_name]
            auc_list = auc_data['mean']
            std_list = auc_data['std']
            pos_offset = adjusted_positions[i]

            plt.bar(np.arange(5, max_features + 1, 5) + pos_offset, auc_list,
                    width=width, label=f'{model_name}',
                    color=color, alpha=0.75, zorder=3)

            plt.errorbar(np.arange(5, max_features + 1, 5) + pos_offset, auc_list,
                         yerr=std_list, fmt='none', color='darkgray', capsize=2)

        plt.legend(loc='upper right', ncol=len(all_models_auc_results))
        plt.grid(axis='y', linestyle='--', zorder=2, alpha=0.3)
        plt.title('AUC Prediction of Age for MG Data Species by Number of Features', fontsize=12)
        plt.xlabel('Number of Selected Species', fontsize=14)
        plt.ylabel('Mean AUC', fontsize=14)
        plt.xticks(np.arange(5, max_features + 1, 5), fontsize=14)
        plt.ylim(0.5, 1)
        plt.yticks([i * 0.05 + 0.5 for i in range(11)], fontsize=14)

        for spine in ax.spines.values():
            spine.set_edgecolor('black')
            spine.set_linewidth(1.5)

        plt.tight_layout()
        plt.savefig(self.output_dir / 'images' / 'feature_selection_auc_barplot.pdf',
                    facecolor='white', edgecolor='none')
        plt.show()


    def run_full_pipeline(self):

        self.load_data()
        self.step1_compare_models()
        self.step2_feature_selection_and_shap()


if __name__ == '__main__':
    DATA_DIR = "/Users/zzhou/Documents/LCSB/Work/MetaDNA_ML"
    OUTPUT_DIR = "./results/age_classification"

    pipeline = AgeClassificationPipeline(DATA_DIR, OUTPUT_DIR, random_state=0)
    pipeline.run_full_pipeline()