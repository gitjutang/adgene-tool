"""
Biomarker Panel and Prediction Tool Development
生物标志物组合和预测工具开发

创新点：
1. 多标志物组合优化
2. 个体化风险评估
3. 动态预测模型
4. 临床决策支持
5. 可视化预测界面

作者：[Author Names]
日期：2026-02-03
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression, ElasticNet
from sklearn.svm import SVC
from sklearn.metrics import (accuracy_score, precision_score, recall_score, 
                           f1_score, roc_auc_score, roc_curve,
                           confusion_matrix, classification_report)
from sklearn.feature_selection import SelectKBest, f_classif, mutual_info_classif
from scipy import stats
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')


class BiomarkerPanel:
    """
    生物标志物组合类
    
    功能：
    - 标志物选择
    - 组合优化
    - 风险评分计算
    - 预测模型构建
    """
    
    def __init__(self):
        self.biomarkers = None
        self.panel = None
        self.model = None
        self.scaler = StandardScaler()
        self.feature_importance = None
        self.risk_scores = None
        
    def select_biomarkers(self, X: pd.DataFrame, y: pd.Series, 
                        method: str = 'combined', top_n: int = 20):
        """
        选择生物标志物
        
        Parameters:
        -----------
        X : pd.DataFrame
            特征矩阵
        y : pd.Series
            标签
        method : str
            选择方法 ('univariate', 'multivariate', 'combined')
        top_n : int
            选择标志物数量
            
        Returns:
        --------
        pd.DataFrame
            选中的标志物及其得分
        """
        print(f"🔍 选择生物标志物 ({method}, top_n={top_n})...")
        
        results = []
        
        if method in ['univariate', 'combined']:
            # 单变量分析
            print("   单变量分析...")
            for feature in X.columns:
                group0 = X[y == 0][feature]
                group1 = X[y == 1][feature]
                
                # t检验
                t_stat, p_value = stats.ttest_ind(group0, group1)
                
                # 效应量
                cohens_d = (group1.mean() - group0.mean()) / np.sqrt(
                    (group0.std()**2 + group1.std()**2) / 2
                )
                
                # ROC AUC
                from sklearn.metrics import roc_auc_score
                try:
                    auc = roc_auc_score(y, X[feature])
                except:
                    auc = 0.5
                
                results.append({
                    'Feature': feature,
                    'T_statistic': t_stat,
                    'P_value': p_value,
                    'Cohens_D': cohens_d,
                    'AUC': auc,
                    'Method': 'Univariate'
                })
        
        if method in ['multivariate', 'combined']:
            # 多变量分析
            print("   多变量分析...")
            from sklearn.ensemble import RandomForestClassifier
            rf = RandomForestClassifier(n_estimators=100, random_state=42)
            rf.fit(X, y)
            
            for i, feature in enumerate(X.columns):
                results.append({
                    'Feature': feature,
                    'Importance': rf.feature_importances_[i],
                    'Method': 'Multivariate'
                })
        
        # 综合评分
        if method == 'combined':
            # 合并单变量和多变量结果
            univariate_df = pd.DataFrame([r for r in results if r['Method'] == 'Univariate'])
            multivariate_df = pd.DataFrame([r for r in results if r['Method'] == 'Multivariate'])
            
            # 标准化得分
            univariate_df['Score'] = (
                (univariate_df['AUC'] - 0.5) * 2 +  # AUC转换为-1到1
                np.abs(univariate_df['Cohens_D']) +   # 效应量
                (1 - univariate_df['P_value'])       # 显著性
            )
            
            multivariate_df['Score'] = multivariate_df['Importance']
            
            # 合并
            combined = pd.merge(univariate_df[['Feature', 'Score']], 
                             multivariate_df[['Feature', 'Score']], 
                             on='Feature', suffixes=('_uni', '_multi'))
            combined['Combined_Score'] = combined['Score_uni'] + combined['Score_multi']
            combined = combined.sort_values('Combined_Score', ascending=False).head(top_n)
            
            self.biomarkers = combined
            
        elif method == 'univariate':
            biomarker_df = pd.DataFrame([r for r in results if r['Method'] == 'Univariate'])
            biomarker_df['Score'] = (
                (biomarker_df['AUC'] - 0.5) * 2 +
                np.abs(biomarker_df['Cohens_D']) +
                (1 - biomarker_df['P_value'])
            )
            biomarker_df = biomarker_df.sort_values('Score', ascending=False).head(top_n)
            self.biomarkers = biomarker_df
            
        elif method == 'multivariate':
            biomarker_df = pd.DataFrame([r for r in results if r['Method'] == 'Multivariate'])
            biomarker_df = biomarker_df.sort_values('Importance', ascending=False).head(top_n)
            self.biomarkers = biomarker_df
        
        print(f"   选择了 {len(self.biomarkers)} 个生物标志物")
        print(f"   Top 5 标志物:")
        print(self.biomarkers.head())
        
        return self.biomarkers
    
    def optimize_panel(self, X: pd.DataFrame, y: pd.Series, 
                     max_biomarkers: int = 10, 
                     min_biomarkers: int = 3,
                     metric: str = 'roc_auc'):
        """
        优化标志物组合
        
        Parameters:
        -----------
        X : pd.DataFrame
            特征矩阵
        y : pd.Series
            标签
        max_biomarkers : int
            最大标志物数
        min_biomarkers : int
            最小标志物数
        metric : str
            评估指标 ('roc_auc', 'f1', 'accuracy')
            
        Returns:
        --------
        Dict
            最优组合结果
        """
        print(f"\n🔧 优化标志物组合 ({min_biomarkers}-{max_biomarkers} 标志物)...")
        
        if self.biomarkers is None:
            raise ValueError("请先选择生物标志物")
        
        candidate_features = self.biomarkers['Feature'].tolist()
        best_score = 0
        best_panel = None
        best_model = None
        
        # 尝试不同大小的组合
        for n in range(min_biomarkers, min(max_biomarkers + 1, len(candidate_features) + 1)):
            print(f"   测试 {n} 个标志物的组合...")
            
            # 使用贪心算法选择组合
            current_panel = []
            remaining_features = candidate_features.copy()
            
            for _ in range(n):
                best_feature = None
                best_feature_score = 0
                
                for feature in remaining_features:
                    test_panel = current_panel + [feature]
                    X_panel = X[test_panel]
                    
                    # 交叉验证
                    model = RandomForestClassifier(n_estimators=100, random_state=42)
                    cv_scores = cross_val_score(model, X_panel, y, cv=5, scoring=metric)
                    avg_score = cv_scores.mean()
                    
                    if avg_score > best_feature_score:
                        best_feature_score = avg_score
                        best_feature = feature
                
                if best_feature is not None:
                    current_panel.append(best_feature)
                    remaining_features.remove(best_feature)
            
            # 评估当前组合
            X_panel = X[current_panel]
            model = RandomForestClassifier(n_estimators=100, random_state=42)
            cv_scores = cross_val_score(model, X_panel, y, cv=5, scoring=metric)
            avg_score = cv_scores.mean()
            
            if avg_score > best_score:
                best_score = avg_score
                best_panel = current_panel
                best_model = model
            
            print(f"      {n} 标志物组合: {metric.upper()} = {avg_score:.4f}")
        
        # 训练最终模型
        X_best = X[best_panel]
        best_model.fit(X_best, y)
        
        self.panel = best_panel
        self.model = best_model
        
        results = {
            'panel': best_panel,
            'n_biomarkers': len(best_panel),
            'score': best_score,
            'model': best_model,
            'feature_names': best_panel
        }
        
        print(f"\n✅ 最优组合: {len(best_panel)} 个标志物")
        print(f"   {metric.upper()}: {best_score:.4f}")
        print(f"   标志物: {', '.join(best_panel)}")
        
        return results
    
    def calculate_risk_score(self, X: pd.DataFrame, method: str = 'linear'):
        """
        计算风险评分
        
        Parameters:
        -----------
        X : pd.DataFrame
            特征矩阵
        method : str
            计算方法 ('linear', 'weighted', 'probability')
            
        Returns:
        --------
        pd.Series
            风险评分
        """
        print(f"\n📊 计算风险评分 ({method})...")
        
        if self.panel is None:
            raise ValueError("请先优化标志物组合")
        
        X_panel = X[self.panel]
        
        if method == 'linear':
            # 线性加权
            if hasattr(self.model, 'coef_'):
                weights = self.model.coef_[0]
            elif hasattr(self.model, 'feature_importances_'):
                weights = self.model.feature_importances_
            else:
                weights = np.ones(len(self.panel))
            
            # 处理NaN值
            weights = np.nan_to_num(weights, nan=1.0)
            risk_scores = (X_panel * weights).sum(axis=1)
            
        elif method == 'weighted':
            # 加权标准化
            X_scaled = self.scaler.fit_transform(X_panel)
            if hasattr(self.model, 'coef_'):
                weights = self.model.coef_[0]
            elif hasattr(self.model, 'feature_importances_'):
                weights = self.model.feature_importances_
            else:
                weights = np.ones(len(self.panel))
            
            # 处理NaN值
            weights = np.nan_to_num(weights, nan=1.0)
            risk_scores = (X_scaled * weights).sum(axis=1)
            
        elif method == 'probability':
            # 模型预测概率
            risk_scores = self.model.predict_proba(X_panel)[:, 1]
        
        # 处理NaN值
        risk_scores = np.nan_to_num(risk_scores, nan=0.0)
        
        # 归一化到0-100
        score_range = risk_scores.max() - risk_scores.min()
        if score_range > 0:
            risk_scores = (risk_scores - risk_scores.min()) / score_range * 100
        else:
            risk_scores = np.zeros_like(risk_scores) + 50
        
        self.risk_scores = risk_scores
        
        print(f"   风险评分范围: {risk_scores.min():.2f} - {risk_scores.max():.2f}")
        print(f"   平均风险评分: {risk_scores.mean():.2f}")
        
        return risk_scores
    
    def categorize_risk(self, risk_scores: pd.Series, 
                      thresholds: tuple = (33, 66)):
        """
        风险分类
        
        Parameters:
        -----------
        risk_scores : pd.Series or np.ndarray
            风险评分
        thresholds : tuple
            分类阈值 (low, high)
            
        Returns:
        --------
        pd.Series
            风险分类
        """
        low_threshold, high_threshold = thresholds
        
        # 转换为Series
        if not isinstance(risk_scores, pd.Series):
            risk_scores = pd.Series(risk_scores)
        
        risk_categories = pd.cut(
            risk_scores,
            bins=[0, low_threshold, high_threshold, 100],
            labels=['Low Risk', 'Medium Risk', 'High Risk']
        )
        
        return risk_categories
    
    def build_prediction_model(self, X: pd.DataFrame, y: pd.Series,
                            model_type: str = 'ensemble'):
        """
        构建预测模型
        
        Parameters:
        -----------
        X : pd.DataFrame
            特征矩阵
        y : pd.Series
            标签
        model_type : str
            模型类型 ('logistic', 'random_forest', 'gradient_boosting', 'svm', 'ensemble')
            
        Returns:
        --------
        Dict
            模型性能
        """
        print(f"\n🤖 构建预测模型 ({model_type})...")
        
        if self.panel is None:
            raise ValueError("请先优化标志物组合")
        
        X_panel = X[self.panel]
        
        # 划分训练集和测试集
        X_train, X_test, y_train, y_test = train_test_split(
            X_panel, y, test_size=0.3, random_state=42, stratify=y
        )
        
        # 标准化
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        # 构建模型
        if model_type == 'logistic':
            model = LogisticRegression(max_iter=1000, random_state=42)
        elif model_type == 'random_forest':
            model = RandomForestClassifier(n_estimators=200, random_state=42)
        elif model_type == 'gradient_boosting':
            model = GradientBoostingClassifier(n_estimators=200, random_state=42)
        elif model_type == 'svm':
            model = SVC(probability=True, random_state=42)
        elif model_type == 'ensemble':
            from sklearn.ensemble import VotingClassifier
            model = VotingClassifier([
                ('lr', LogisticRegression(max_iter=1000, random_state=42)),
                ('rf', RandomForestClassifier(n_estimators=200, random_state=42)),
                ('gb', GradientBoostingClassifier(n_estimators=200, random_state=42))
            ], voting='soft')
        else:
            raise ValueError(f"未知的模型类型: {model_type}")
        
        # 训练模型
        model.fit(X_train_scaled, y_train)
        self.model = model
        
        # 预测
        y_pred = model.predict(X_test_scaled)
        y_pred_prob = model.predict_proba(X_test_scaled)[:, 1]
        
        # 评估
        performance = {
            'accuracy': accuracy_score(y_test, y_pred),
            'precision': precision_score(y_test, y_pred),
            'recall': recall_score(y_test, y_pred),
            'f1': f1_score(y_test, y_pred),
            'auc': roc_auc_score(y_test, y_pred_prob)
        }
        
        print(f"   准确率: {performance['accuracy']:.4f}")
        print(f"   精确率: {performance['precision']:.4f}")
        print(f"   召回率: {performance['recall']:.4f}")
        print(f"   F1分数: {performance['f1']:.4f}")
        print(f"   AUC: {performance['auc']:.4f}")
        
        return performance
    
    def visualize_panel(self, X: pd.DataFrame, y: pd.Series,
                     output_dir: str = './results/biomarker_panel'):
        """
        可视化标志物组合
        
        Parameters:
        -----------
        X : pd.DataFrame
            特征矩阵
        y : pd.Series
            标签
        output_dir : str
            输出目录
        """
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        fig, axes = plt.subplots(3, 3, figsize=(18, 15))
        
        # 1. 标志物重要性
        if self.biomarkers is not None:
            ax = axes[0, 0]
            top_biomarkers = self.biomarkers.head(15)
            if 'Combined_Score' in top_biomarkers.columns:
                ax.barh(range(len(top_biomarkers)), top_biomarkers['Combined_Score'])
            elif 'Score' in top_biomarkers.columns:
                ax.barh(range(len(top_biomarkers)), top_biomarkers['Score'])
            else:
                ax.barh(range(len(top_biomarkers)), top_biomarkers['Importance'])
            ax.set_yticks(range(len(top_biomarkers)))
            ax.set_yticklabels(top_biomarkers['Feature'], fontsize=8)
            ax.set_xlabel('Score')
            ax.set_title('Biomarker Importance')
            ax.grid(True, alpha=0.3, axis='x')
        
        # 2. 单个标志物ROC曲线
        if self.biomarkers is not None:
            ax = axes[0, 1]
            from sklearn.metrics import roc_curve, auc
            for feature in self.biomarkers['Feature'].head(5):
                fpr, tpr, _ = roc_curve(y, X[feature])
                roc_auc = auc(fpr, tpr)
                ax.plot(fpr, tpr, label=f'{feature} (AUC={roc_auc:.3f})', linewidth=2)
            ax.plot([0, 1], [0, 1], 'k--', linewidth=1)
            ax.set_xlabel('False Positive Rate')
            ax.set_ylabel('True Positive Rate')
            ax.set_title('Individual Biomarker ROC Curves')
            ax.legend(fontsize=6)
            ax.grid(True, alpha=0.3)
        
        # 3. 标志物分布箱线图
        if self.panel is not None:
            ax = axes[0, 2]
            X_panel = X[self.panel]
            n_show = min(5, len(self.panel))
            for i, biomarker in enumerate(self.panel[:n_show]):
                data0 = X_panel[y == 0][biomarker].values
                data1 = X_panel[y == 1][biomarker].values
                positions = [i*2 + 0.5, i*2 + 1.5]
                bp = ax.boxplot([data0, data1], positions=positions, widths=0.8, patch_artist=True)
                bp['boxes'][0].set_facecolor('lightblue')
                bp['boxes'][1].set_facecolor('lightcoral')
            ax.set_xticks(range(2, n_show*2 + 1, 2))
            ax.set_xticklabels(self.panel[:n_show], rotation=15, ha='right')
            ax.set_ylabel('Value')
            ax.set_title('Biomarker Distribution by Group')
            ax.legend(['Control', 'AD'], loc='upper right')
        
        # 4. 标志物相关性热图
        if self.panel is not None:
            ax = axes[1, 0]
            X_panel = X[self.panel]
            corr_matrix = X_panel.corr()
            sns.heatmap(corr_matrix, cmap='coolwarm', center=0,
                       ax=ax, cbar_kws={'label': 'Correlation'})
            ax.set_title('Biomarker Correlation Heatmap')
        
        # 5. 组合ROC曲线
        if self.model is not None and self.panel is not None:
            ax = axes[1, 1]
            X_panel = X[self.panel]
            X_scaled = self.scaler.transform(X_panel)
            y_pred_prob = self.model.predict_proba(X_scaled)[:, 1]
            fpr, tpr, _ = roc_curve(y, y_pred_prob)
            roc_auc = auc(fpr, tpr)
            ax.plot(fpr, tpr, linewidth=2, label=f'Panel (AUC={roc_auc:.3f})')
            ax.plot([0, 1], [0, 1], 'k--', linewidth=1)
            ax.set_xlabel('False Positive Rate')
            ax.set_ylabel('True Positive Rate')
            ax.set_title('Combined Panel ROC Curve')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        # 6. 风险评分分布
        if self.risk_scores is not None:
            ax = axes[1, 2]
            ax.hist(self.risk_scores[y == 0], bins=30, alpha=0.6, 
                   label='Control', color='blue')
            ax.hist(self.risk_scores[y == 1], bins=30, alpha=0.6, 
                   label='AD', color='red')
            ax.set_xlabel('Risk Score')
            ax.set_ylabel('Frequency')
            ax.set_title('Risk Score Distribution')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        # 7. 风险分类
        if self.risk_scores is not None:
            ax = axes[2, 0]
            risk_categories = self.categorize_risk(self.risk_scores)
            category_counts = risk_categories.value_counts()
            colors = ['green', 'orange', 'red']
            ax.pie(category_counts.values, labels=category_counts.index, 
                   autopct='%1.1f%%', colors=colors)
            ax.set_title('Risk Category Distribution')
        
        # 8. 混淆矩阵
        if self.model is not None and self.panel is not None:
            ax = axes[2, 1]
            X_panel = X[self.panel]
            X_scaled = self.scaler.transform(X_panel)
            y_pred = self.model.predict(X_scaled)
            cm = confusion_matrix(y, y_pred)
            sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax,
                       cbar_kws={'label': 'Count'})
            ax.set_xlabel('Predicted')
            ax.set_ylabel('Actual')
            ax.set_title('Confusion Matrix')
        
        # 9. 标志物组合摘要
        ax = axes[2, 2]
        ax.axis('off')
        
        summary_text = "Biomarker Panel Summary\n\n"
        
        if self.panel is not None:
            summary_text += f"Panel Size: {len(self.panel)} biomarkers\n\n"
            summary_text += "Biomarkers:\n"
            for i, biomarker in enumerate(self.panel, 1):
                summary_text += f"{i}. {biomarker}\n"
        
        if self.model is not None:
            X_panel = X[self.panel]
            X_scaled = self.scaler.transform(X_panel)
            y_pred = self.model.predict(X_scaled)
            accuracy = accuracy_score(y, y_pred)
            summary_text += f"\nModel Accuracy: {accuracy:.2%}\n"
        
        if self.risk_scores is not None:
            summary_text += f"\nRisk Score Range: {self.risk_scores.min():.1f} - {self.risk_scores.max():.1f}\n"
            summary_text += f"Mean Risk Score: {self.risk_scores.mean():.1f}\n"
        
        ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top', family='monospace')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/biomarker_panel.png', dpi=300, bbox_inches='tight')
        print(f"✅ 标志物组合可视化保存到: {output_dir}/biomarker_panel.png")
        
        plt.close()


class IndividualRiskAssessment:
    """
    个体化风险评估类
    
    功能：
    - 个体风险计算
    - 风险因素分析
    - 预防建议生成
    """
    
    def __init__(self, biomarker_panel: BiomarkerPanel):
        self.panel = biomarker_panel
        
    def assess_individual(self, individual_data: pd.DataFrame):
        """
        评估个体风险
        
        Parameters:
        -----------
        individual_data : pd.DataFrame
            个体数据
            
        Returns:
        --------
        Dict
            风险评估结果
        """
        print("\n👤 个体化风险评估...")
        
        # 计算风险评分（使用probability方法，因为集成模型支持）
        risk_score = self.panel.calculate_risk_score(individual_data, method='probability')
        
        # 风险分类
        risk_category = self.panel.categorize_risk(risk_score)
        
        # 风险因素分析
        risk_factors = self._analyze_risk_factors(individual_data)
        
        # 生成建议
        recommendations = self._generate_recommendations(risk_category, risk_factors)
        
        # 处理返回值
        if hasattr(risk_score, 'iloc'):
            risk_score_value = risk_score.iloc[0]
        else:
            risk_score_value = risk_score[0] if len(risk_score) > 0 else 0
        
        if hasattr(risk_category, 'iloc'):
            risk_category_value = risk_category.iloc[0]
        else:
            risk_category_value = risk_category
        
        results = {
            'risk_score': risk_score_value,
            'risk_category': risk_category_value,
            'risk_factors': risk_factors,
            'recommendations': recommendations
        }
        
        print(f"   风险评分: {risk_score_value:.2f}")
        print(f"   风险分类: {risk_category_value}")
        
        return results
    
    def _analyze_risk_factors(self, individual_data: pd.DataFrame):
        """
        分析风险因素
        """
        risk_factors = []
        
        if self.panel.panel is not None:
            for biomarker in self.panel.panel:
                value = individual_data[biomarker].iloc[0]
                # 简化的风险因素判断
                if value > individual_data[biomarker].quantile(0.75):
                    risk_level = 'High'
                elif value > individual_data[biomarker].quantile(0.5):
                    risk_level = 'Medium'
                else:
                    risk_level = 'Low'
                
                risk_factors.append({
                    'Biomarker': biomarker,
                    'Value': value,
                    'Risk_Level': risk_level
                })
        
        return risk_factors
    
    def _generate_recommendations(self, risk_category, risk_factors: list):
        """
        生成预防建议
        """
        recommendations = []
        
        # 处理pandas Series
        if hasattr(risk_category, 'iloc'):
            risk_category_str = risk_category.iloc[0]
        else:
            risk_category_str = risk_category
        
        if risk_category_str == 'High Risk':
            recommendations.append("建议进行进一步临床检查")
            recommendations.append("考虑生活方式干预")
            recommendations.append("定期监测相关指标")
        elif risk_category_str == 'Medium Risk':
            recommendations.append("建议定期随访")
            recommendations.append("保持健康生活方式")
            recommendations.append("关注相关指标变化")
        else:
            recommendations.append("继续保持健康生活方式")
            recommendations.append("定期体检")
        
        return recommendations


def demo_biomarker_panel():
    """
    演示生物标志物组合开发
    """
    print("=" * 60)
    print("生物标志物组合和预测工具开发")
    print("=" * 60)
    
    # 创建模拟数据
    np.random.seed(42)
    n_samples = 300
    n_features = 50
    
    X = pd.DataFrame(
        np.random.randn(n_samples, n_features),
        columns=[f'Biomarker_{i}' for i in range(n_features)]
    )
    
    # 创建标签（AD vs Control）
    y = pd.Series(
        np.random.choice([0, 1], n_samples, p=[0.6, 0.4])
    )
    
    # 创建生物标志物组合分析器
    panel_analyzer = BiomarkerPanel()
    
    # 选择生物标志物
    biomarkers = panel_analyzer.select_biomarkers(X, y, method='combined', top_n=20)
    
    # 优化标志物组合
    optimized_panel = panel_analyzer.optimize_panel(X, y, max_biomarkers=10, min_biomarkers=3)
    
    # 计算风险评分
    risk_scores = panel_analyzer.calculate_risk_score(X, method='probability')
    
    # 构建预测模型
    performance = panel_analyzer.build_prediction_model(X, y, model_type='ensemble')
    
    # 可视化
    panel_analyzer.visualize_panel(X, y)
    
    # 个体化风险评估
    individual_assessor = IndividualRiskAssessment(panel_analyzer)
    individual_data = X.iloc[[0]]
    assessment = individual_assessor.assess_individual(individual_data)
    
    print("\n" + "=" * 60)
    print("✅ 生物标志物组合开发演示完成！")
    print("=" * 60)


if __name__ == '__main__':
    demo_biomarker_panel()
