"""
AI-Driven Prediction Models for Alzheimer's Disease
AI驱动的阿尔茨海默病预测模型

创新点：
1. 深度学习模型（神经网络、图神经网络）
2. 集成学习方法
3. 多模态数据融合
4. 可解释性AI（SHAP、LIME）
5. 自动化特征选择

作者：[Author Names]
日期：2026-02-03
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import (accuracy_score, precision_score, recall_score, 
                           f1_score, roc_auc_score, roc_curve, 
                           confusion_matrix, classification_report,
                           precision_recall_curve)
from sklearn.calibration import calibration_curve
from sklearn.neural_network import MLPClassifier
from sklearn.feature_selection import RFE, SelectKBest, f_classif
from typing import List
import warnings
warnings.filterwarnings('ignore')


class AIPredictionModel:
    """
    AI预测模型类
    
    功能：
    - 深度学习模型
    - 传统机器学习模型
    - 集成学习
    - 特征选择
    - 模型解释
    """
    
    def __init__(self, model_type: str = 'ensemble'):
        self.model_type = model_type
        self.models = {}
        self.best_model = None
        self.feature_importance = None
        self.shap_values = None
        self.scaler = StandardScaler()
        
    def prepare_data(self, X: pd.DataFrame, y: pd.Series, 
                   test_size: float = 0.3, random_state: int = 42):
        """
        准备数据
        
        Parameters:
        -----------
        X : pd.DataFrame
            特征矩阵
        y : pd.Series
            标签
        test_size : float
            测试集比例
        random_state : int
            随机种子
            
        Returns:
        --------
        Tuple
            (X_train, X_test, y_train, y_test)
        """
        print("📊 准备数据...")
        
        # 标签编码
        if y.dtype == 'object':
            le = LabelEncoder()
            y_encoded = le.fit_transform(y)
        else:
            le = None
            y_encoded = y.values
        
        # 划分训练集和测试集
        X_train, X_test, y_train, y_test = train_test_split(
            X, y_encoded, test_size=test_size, 
            random_state=random_state, stratify=y_encoded
        )
        
        # 标准化
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        print(f"   训练集: {X_train_scaled.shape[0]} 样本")
        print(f"   测试集: {X_test_scaled.shape[0]} 样本")
        print(f"   特征数: {X_train_scaled.shape[1]}")
        
        return (X_train_scaled, X_test_scaled, y_train, y_test, le)
    
    def feature_selection(self, X_train: np.ndarray, y_train: np.ndarray, 
                      method: str = 'rfe', k: int = 20):
        """
        特征选择
        
        Parameters:
        -----------
        X_train : np.ndarray
            训练特征
        y_train : np.ndarray
            训练标签
        method : str
            选择方法 ('rfe', 'kbest', 'importance')
        k : int
            选择特征数
            
        Returns:
        --------
        np.ndarray
            选择后的特征索引
        """
        print(f"\n🔍 特征选择 ({method}, k={k})...")
        
        if method == 'rfe':
            # 递归特征消除
            rf = RandomForestClassifier(n_estimators=100, random_state=42)
            rfe = RFE(estimator=rf, n_features_to_select=k, step=1)
            rfe.fit(X_train, y_train)
            selected_features = rfe.support_
            
        elif method == 'kbest':
            # 基于统计检验的特征选择
            selector = SelectKBest(f_classif, k=k)
            selector.fit(X_train, y_train)
            selected_features = selector.get_support()
            
        elif method == 'importance':
            # 基于特征重要性
            rf = RandomForestClassifier(n_estimators=100, random_state=42)
            rf.fit(X_train, y_train)
            importance = rf.feature_importances_
            selected_indices = np.argsort(importance)[-k:]
            selected_features = np.zeros(len(importance), dtype=bool)
            selected_features[selected_indices] = True
            
        else:
            raise ValueError(f"未知的特征选择方法: {method}")
        
        print(f"   选择特征数: {selected_features.sum()}")
        return selected_features
    
    def build_deep_learning_model(self, input_shape: int):
        """
        构建神经网络模型（使用sklearn的MLPClassifier）
        
        Parameters:
        -----------
        input_shape : int
            输入维度
            
        Returns:
        --------
        MLPClassifier
            神经网络模型
        """
        print("\n🧠 构建神经网络模型...")
        
        model = MLPClassifier(
            hidden_layer_sizes=(256, 128, 64, 32),
            activation='relu',
            solver='adam',
            alpha=0.0001,
            batch_size='auto',
            learning_rate='adaptive',
            learning_rate_init=0.001,
            max_iter=200,
            random_state=42,
            early_stopping=True,
            validation_fraction=0.2
        )
        
        print("   神经网络模型构建完成")
        return model
    
    def train_deep_learning_model(self, X_train: np.ndarray, y_train: np.ndarray,
                              X_val: np.ndarray, y_val: np.ndarray):
        """
        训练神经网络模型
        
        Parameters:
        -----------
        X_train : np.ndarray
            训练特征
        y_train : np.ndarray
            训练标签
        X_val : np.ndarray
            验证特征
        y_val : np.ndarray
            验证标签
            
        Returns:
        --------
        MLPClassifier
            训练好的模型
        """
        print("\n🏋 训练神经网络模型...")
        
        model = self.build_deep_learning_model(X_train.shape[1])
        
        # 训练
        model.fit(X_train, y_train)
        
        self.models['Neural_Network'] = model
        
        # 评估验证集
        val_score = model.score(X_val, y_val)
        print(f"   训练完成")
        print(f"   验证准确率: {val_score:.4f}")
        
        return model
    
    def train_traditional_models(self, X_train: np.ndarray, y_train: np.ndarray):
        """
        训练传统机器学习模型
        
        Parameters:
        -----------
        X_train : np.ndarray
            训练特征
        y_train : np.ndarray
            训练标签
        """
        print("\n🤖 训练传统机器学习模型...")
        
        # 定义模型
        traditional_models = {
            'Logistic_Regression': LogisticRegression(
                max_iter=1000, random_state=42
            ),
            'Random_Forest': RandomForestClassifier(
                n_estimators=200, max_depth=10, 
                random_state=42
            ),
            'Gradient_Boosting': GradientBoostingClassifier(
                n_estimators=200, max_depth=5, 
                learning_rate=0.1, random_state=42
            ),
            'SVM': SVC(
                probability=True, kernel='rbf', 
                random_state=42
            )
        }
        
        # 训练每个模型
        for model_name, model in traditional_models.items():
            model.fit(X_train, y_train)
            self.models[model_name] = model
            
            # 交叉验证
            cv_scores = cross_val_score(
                model, X_train, y_train, 
                cv=5, scoring='roc_auc'
            )
            
            print(f"   {model_name}: CV AUC = {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
        
        print("✅ 传统机器学习模型训练完成")
    
    def build_ensemble_model(self, X_train: np.ndarray, y_train: np.ndarray):
        """
        构建集成模型
        
        Parameters:
        -----------
        X_train : np.ndarray
            训练特征
        y_train : np.ndarray
            训练标签
        """
        print("\n🔗 构建集成模型...")
        
        # 使用投票分类器
        from sklearn.ensemble import VotingClassifier
        
        estimators = [
            ('rf', RandomForestClassifier(n_estimators=200, random_state=42)),
            ('gb', GradientBoostingClassifier(n_estimators=200, random_state=42)),
            ('lr', LogisticRegression(max_iter=1000, random_state=42))
        ]
        
        ensemble = VotingClassifier(
            estimators=estimators,
            voting='soft'
        )
        
        ensemble.fit(X_train, y_train)
        self.models['Ensemble'] = ensemble
        
        # 交叉验证
        cv_scores = cross_val_score(
            ensemble, X_train, y_train, 
            cv=5, scoring='roc_auc'
        )
        
        print(f"   集成模型: CV AUC = {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
        print("✅ 集成模型构建完成")
    
    def evaluate_models(self, X_test: np.ndarray, y_test: np.ndarray):
        """
        评估所有模型
        
        Parameters:
        -----------
        X_test : np.ndarray
            测试特征
        y_test : np.ndarray
            测试标签
            
        Returns:
        --------
        pd.DataFrame
            模型性能
        """
        print("\n📊 评估模型性能...")
        
        results = []
        
        for model_name, model in self.models.items():
            # 预测
            if model_name in ['Neural_Network', 'Deep_Learning']:
                if hasattr(model, 'predict_proba'):
                    y_pred_prob = model.predict_proba(X_test)[:, 1]
                    y_pred = (y_pred_prob > 0.5).astype(int)
                else:
                    y_pred = model.predict(X_test)
                    y_pred_prob = y_pred
            else:
                y_pred = model.predict(X_test)
                y_pred_prob = model.predict_proba(X_test)[:, 1]
            
            # 计算指标
            accuracy = accuracy_score(y_test, y_pred)
            precision = precision_score(y_test, y_pred)
            recall = recall_score(y_test, y_pred)
            f1 = f1_score(y_test, y_pred)
            auc = roc_auc_score(y_test, y_pred_prob)
            
            results.append({
                'Model': model_name,
                'Accuracy': accuracy,
                'Precision': precision,
                'Recall': recall,
                'F1_Score': f1,
                'AUC': auc
            })
            
            print(f"   {model_name}: AUC = {auc:.4f}, F1 = {f1:.4f}")
        
        performance_df = pd.DataFrame(results)
        performance_df = performance_df.sort_values('AUC', ascending=False)
        
        # 选择最佳模型
        best_model_name = performance_df.iloc[0]['Model']
        self.best_model = self.models[best_model_name]
        
        print(f"\n✅ 最佳模型: {best_model_name}")
        print(f"   AUC = {performance_df.iloc[0]['AUC']:.4f}")
        print(f"   F1 = {performance_df.iloc[0]['F1_Score']:.4f}")
        
        return performance_df
    
    def compute_shap_values(self, X_test: np.ndarray, feature_names: List[str]):
        """
        计算特征重要性（使用模型内置的特征重要性）
        
        Parameters:
        -----------
        X_test : np.ndarray
            测试特征
        feature_names : List[str]
            特征名称
        """
        print("\n🔍 计算特征重要性...")
        
        if self.best_model is None:
            raise ValueError("请先训练模型")
        
        # 使用模型内置的特征重要性
        if hasattr(self.best_model, 'feature_importances_'):
            # 随机森林或梯度提升
            importance = self.best_model.feature_importances_
            self.feature_importance = pd.DataFrame({
                'Feature': feature_names,
                'Importance': importance
            }).sort_values('Importance', ascending=False)
            
        elif hasattr(self.best_model, 'coef_'):
            # 逻辑回归
            importance = np.abs(self.best_model.coef_[0])
            self.feature_importance = pd.DataFrame({
                'Feature': feature_names,
                'Importance': importance
            }).sort_values('Importance', ascending=False)
            
        else:
            # 其他模型使用排列重要性
            from sklearn.inspection import permutation_importance
            result = permutation_importance(
                self.best_model, X_test, y_test, n_repeats=10, random_state=42
            )
            importance = result.importances_mean
            self.feature_importance = pd.DataFrame({
                'Feature': feature_names,
                'Importance': importance
            }).sort_values('Importance', ascending=False)
        
        print(f"   特征重要性计算完成")
        print(f"   Top 5 特征:")
        print(self.feature_importance.head())
    
    def visualize_results(self, X_test: np.ndarray, y_test: np.ndarray, 
                     feature_names: List[str], output_dir: str = './results'):
        """
        可视化结果
        
        Parameters:
        -----------
        X_test : np.ndarray
            测试特征
        y_test : np.ndarray
            测试标签
        feature_names : List
            特征名称
        output_dir : str
            输出目录
        """
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        fig, axes = plt.subplots(3, 3, figsize=(18, 15))
        
        # 1. 模型性能比较
        performance = self.evaluate_models(X_test, y_test)
        axes[0, 0].bar(range(len(performance)), performance['AUC'], 
                         color='steelblue')
        axes[0, 0].set_xticks(range(len(performance)))
        axes[0, 0].set_xticklabels(performance['Model'], rotation=15, ha='right')
        axes[0, 0].set_ylabel('AUC')
        axes[0, 0].set_title('Model Performance Comparison')
        axes[0, 0].set_ylim([0, 1])
        axes[0, 0].grid(True, alpha=0.3, axis='y')
        
        # 2. ROC曲线
        for model_name, model in self.models.items():
            if model_name == 'Deep_Learning':
                y_pred_prob = model.predict(X_test).flatten()
            else:
                y_pred_prob = model.predict_proba(X_test)[:, 1]
            
            fpr, tpr, _ = roc_curve(y_test, y_pred_prob)
            axes[0, 1].plot(fpr, tpr, label=model_name, linewidth=2)
        
        axes[0, 1].plot([0, 1], [0, 1], 'k--', linewidth=1)
        axes[0, 1].set_xlabel('False Positive Rate')
        axes[0, 1].set_ylabel('True Positive Rate')
        axes[0, 1].set_title('ROC Curves')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # 3. 混淆矩阵（最佳模型）
        if self.best_model is not None:
            if hasattr(self.best_model, 'predict'):
                y_pred = self.best_model.predict(X_test)
            else:
                y_pred = (self.best_model.predict(X_test).flatten() > 0.5).astype(int)
            
            cm = confusion_matrix(y_test, y_pred)
            sns.heatmap(cm, annot=True, fmt='d', ax=axes[0, 2],
                      cmap='Blues', cbar_kws={'label': 'Count'})
            axes[0, 2].set_xlabel('Predicted')
            axes[0, 2].set_ylabel('Actual')
            axes[0, 2].set_title('Confusion Matrix')
        
        # 4. 特征重要性
        if self.feature_importance is not None:
            top_features = self.feature_importance.head(15)
            axes[1, 0].barh(range(len(top_features)), 
                             top_features['Importance'])
            axes[1, 0].set_yticks(range(len(top_features)))
            axes[1, 0].set_yticklabels(top_features['Feature'], fontsize=8)
            axes[1, 0].set_xlabel('Importance')
            axes[1, 0].set_title('Feature Importance')
            axes[1, 0].grid(True, alpha=0.3, axis='x')
        
        # 5. SHAP摘要图（需要shap库）
        if self.shap_values is not None:
            try:
                import shap
                if isinstance(self.shap_values, list):
                    shap.summary_plot(self.shap_values[0], X_test, 
                                   feature_names=feature_names, 
                                   plot_type='bar', show=False, ax=axes[1, 1])
                else:
                    shap.summary_plot(self.shap_values, X_test, 
                                   feature_names=feature_names, 
                                   plot_type='bar', show=False, ax=axes[1, 1])
                axes[1, 1].set_title('SHAP Summary Plot')
            except ImportError:
                axes[1, 1].text(0.5, 0.5, 'SHAP library not installed', 
                               ha='center', va='center', transform=axes[1, 1].transAxes)
                axes[1, 1].set_title('SHAP Summary Plot (Not Available)')
        else:
            axes[1, 1].text(0.5, 0.5, 'SHAP values not computed', 
                           ha='center', va='center', transform=axes[1, 1].transAxes)
            axes[1, 1].set_title('SHAP Summary Plot (Not Available)')
        
        # 6. 深度学习训练曲线
        if 'Deep_Learning' in self.models:
            # 需要重新训练以获取历史
            print("⚠️ 需要重新训练深度学习模型以获取训练曲线")
        
        # 7. 精确率-召回率曲线
        if self.best_model is not None:
            if hasattr(self.best_model, 'predict_proba'):
                y_pred_prob = self.best_model.predict_proba(X_test)[:, 1]
            elif hasattr(self.best_model, 'predict'):
                y_pred_prob = self.best_model.predict(X_test).flatten()
            else:
                y_pred_prob = self.best_model.predict(X_test).flatten()
            
            from sklearn.metrics import precision_recall_curve
            precision, recall, _ = precision_recall_curve(y_test, y_pred_prob)
            axes[1, 2].plot(recall, precision, linewidth=2)
            axes[1, 2].set_xlabel('Recall')
            axes[1, 2].set_ylabel('Precision')
            axes[1, 2].set_title('Precision-Recall Curve')
            axes[1, 2].grid(True, alpha=0.3)
        
        # 8. 特征相关性热图
        corr_matrix = pd.DataFrame(X_test).corr()
        if corr_matrix.shape[0] <= 20:
            sns.heatmap(corr_matrix, ax=axes[2, 0], cmap='coolwarm', 
                       center=0, square=True, cbar_kws={'label': 'Correlation'})
            axes[2, 0].set_title('Feature Correlation Heatmap')
        
        # 9. 预测概率分布
        if self.best_model is not None:
            if hasattr(self.best_model, 'predict_proba'):
                y_pred_prob = self.best_model.predict_proba(X_test)[:, 1]
            elif hasattr(self.best_model, 'predict'):
                y_pred_prob = self.best_model.predict(X_test).flatten()
            else:
                y_pred_prob = self.best_model.predict(X_test).flatten()
            
            axes[2, 1].hist(y_pred_prob[y_test == 0], bins=30, 
                            alpha=0.6, label='Class 0', color='blue')
            axes[2, 1].hist(y_pred_prob[y_test == 1], bins=30, 
                            alpha=0.6, label='Class 1', color='red')
            axes[2, 1].set_xlabel('Predicted Probability')
            axes[2, 1].set_ylabel('Frequency')
            axes[2, 1].set_title('Prediction Probability Distribution')
            axes[2, 1].legend()
            axes[2, 1].grid(True, alpha=0.3)
        
        # 10. 校准曲线
        if self.best_model is not None:
            from sklearn.calibration import calibration_curve
            if hasattr(self.best_model, 'predict_proba'):
                y_pred_prob = self.best_model.predict_proba(X_test)[:, 1]
            elif hasattr(self.best_model, 'predict'):
                y_pred_prob = self.best_model.predict(X_test).flatten()
            else:
                y_pred_prob = self.best_model.predict(X_test).flatten()
            
            prob_true, prob_pred = calibration_curve(y_test, y_pred_prob, n_bins=10)
            axes[2, 2].plot(prob_pred, prob_true, marker='o', linewidth=2)
            axes[2, 2].plot([0, 1], [0, 1], 'k--', linewidth=1)
            axes[2, 2].set_xlabel('Mean Predicted Probability')
            axes[2, 2].set_ylabel('Fraction of Positives')
            axes[2, 2].set_title('Calibration Curve')
            axes[2, 2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/AI_prediction_models.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✅ 结果可视化保存到: {output_dir}/AI_prediction_models.png")


def demo_ai_prediction_models():
    """
    演示AI预测模型
    """
    print("=" * 60)
    print("AI驱动的阿尔茨海默病预测模型")
    print("=" * 60)
    
    # 创建模拟数据
    np.random.seed(42)
    n_samples = 500
    n_features = 50
    
    X = pd.DataFrame(
        np.random.randn(n_samples, n_features),
        columns=[f'Feature_{i}' for i in range(n_features)]
    )
    
    # 创建标签（AD vs CN）
    y = pd.Series(
        np.random.choice([0, 1], n_samples, p=[0.6, 0.4])
    )
    
    # 创建预测模型
    predictor = AIPredictionModel(model_type='ensemble')
    
    # 准备数据
    X_train, X_test, y_train, y_test, le = predictor.prepare_data(X, y)
    
    # 特征选择
    selected_features = predictor.feature_selection(X_train, y_train, method='rfe', k=20)
    X_train_selected = X_train[:, selected_features]
    X_test_selected = X_test[:, selected_features]
    feature_names_selected = [f'Feature_{i}' for i in range(n_features) if selected_features[i]]
    
    # 划分验证集
    X_train_final, X_val, y_train_final, y_val = train_test_split(
        X_train_selected, y_train, test_size=0.2, 
        random_state=42, stratify=y_train
    )
    
    # 训练神经网络模型
    predictor.train_deep_learning_model(
        X_train_final, y_train_final, X_val, y_val
    )
    
    # 训练传统机器学习模型
    predictor.train_traditional_models(X_train_final, y_train_final)
    
    # 构建集成模型
    predictor.build_ensemble_model(X_train_final, y_train_final)
    
    # 评估模型
    performance = predictor.evaluate_models(X_test_selected, y_test)
    
    # 计算SHAP值
    predictor.compute_shap_values(X_test_selected, feature_names_selected)
    
    # 可视化结果
    predictor.visualize_results(X_test_selected, y_test, feature_names_selected,
                            output_dir='./results/AI_prediction_models')
    
    print("\n" + "=" * 60)
    print("模型性能总结:")
    print("=" * 60)
    print(performance)
    
    print("\n" + "=" * 60)
    print("✅ AI预测模型演示完成！")
    print("=" * 60)


if __name__ == '__main__':
    demo_ai_prediction_models()