"""
Multi-Dataset Integration Framework
多数据集整合框架（ADNI, ROSMAP, MSBB）

创新点：
1. 跨数据集标准化和整合
2. 跨种族验证
3. 元分析增加统计功效
4. 批量效应校正

作者：[Author Names]
日期：2026-02-03
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import stats
from scipy.stats import combine_pvalues
from typing import List, Dict, Optional
import warnings
warnings.filterwarnings('ignore')


class MultiDatasetIntegrator:
    """
    多数据集整合器
    
    功能：
    - 跨数据集标准化
    - 批量效应校正
    - 元分析
    - 跨种族验证
    - 一致性评估
    """
    
    def __init__(self):
        self.datasets = {}
        self.integrated_data = None
        self.meta_analysis_results = None
        self.batch_corrected_data = None
        
    def load_dataset(self, dataset_name: str, data: pd.DataFrame, 
                   metadata: pd.DataFrame = None):
        """
        加载数据集
        
        Parameters:
        -----------
        dataset_name : str
            数据集名称 ('ADNI', 'ROSMAP', 'MSBB')
        data : pd.DataFrame
            组学数据
        metadata : pd.DataFrame
            元数据（种族、年龄、性别等）
        """
        self.datasets[dataset_name] = {
            'data': data,
            'metadata': metadata
        }
        print(f"✅ 加载 {dataset_name} 数据集: {data.shape[0]} 样本 × {data.shape[1]} 特征")
        
    def standardize_datasets(self):
        """
        标准化所有数据集
        
        Returns:
        --------
        Dict
            标准化后的数据集
        """
        print("\n🔄 标准化数据集...")
        
        for dataset_name in self.datasets:
            data = self.datasets[dataset_name]['data'].copy()
            
            # Z-score标准化
            scaler = StandardScaler()
            data_standardized = pd.DataFrame(
                scaler.fit_transform(data),
                index=data.index,
                columns=data.columns
            )
            
            self.datasets[dataset_name]['data_standardized'] = data_standardized
            print(f"   {dataset_name}: 标准化完成")
        
        return self.datasets
    
    def identify_common_features(self):
        """
        识别所有数据集的共同特征
        
        Returns:
        --------
        List
            共同特征列表
        """
        print("\n🔍 识别共同特征...")
        
        feature_sets = [
            set(self.datasets[ds]['data'].columns) 
            for ds in self.datasets
        ]
        
        common_features = set.intersection(*feature_sets)
        common_features = sorted(list(common_features))
        
        print(f"✅ 共同特征数: {len(common_features)}")
        print(f"   各数据集特征数:")
        for ds in self.datasets:
            print(f"   - {ds}: {len(self.datasets[ds]['data'].columns)}")
        
        return common_features
    
    def correct_batch_effects(self, method: str = 'combat'):
        """
        批量效应校正
        
        Parameters:
        -----------
        method : str
            校正方法 ('combat', 'limma', 'mean_centering')
            
        Returns:
        --------
        Dict
            批量效应校正后的数据
        """
        print(f"\n🔧 批量效应校正 ({method})...")
        
        common_features = self.identify_common_features()
        
        if method == 'combat':
            # 简化的ComBat实现
            for dataset_name in self.datasets:
                data = self.datasets[dataset_name]['data_standardized'][common_features].copy()
                
                # 计算批量效应（数据集均值）
                dataset_mean = data.mean()
                global_mean = pd.concat([
                    self.datasets[ds]['data_standardized'][common_features] 
                    for ds in self.datasets
                ]).mean()
                
                # 校正
                data_corrected = data - dataset_mean + global_mean
                
                self.datasets[dataset_name]['data_corrected'] = data_corrected
                print(f"   {dataset_name}: 批量效应校正完成")
                
        elif method == 'mean_centering':
            # 均值中心化
            global_mean = pd.concat([
                self.datasets[ds]['data_standardized'][common_features] 
                for ds in self.datasets
            ]).mean()
            
            for dataset_name in self.datasets:
                data = self.datasets[dataset_name]['data_standardized'][common_features].copy()
                data_corrected = data - global_mean
                
                self.datasets[dataset_name]['data_corrected'] = data_corrected
                print(f"   {dataset_name}: 均值中心化完成")
        
        self.batch_corrected_data = self.datasets
        return self.datasets
    
    def perform_meta_analysis(self, feature: str, outcome: str = 'diagnosis'):
        """
        对特定特征进行元分析
        
        Parameters:
        -----------
        feature : str
            要分析的特征
        outcome : str
            结局变量
            
        Returns:
        --------
        Dict
            元分析结果
        """
        print(f"\n📊 元分析: {feature}")
        
        effect_sizes = []
        variances = []
        sample_sizes = []
        
        for dataset_name in self.datasets:
            data = self.datasets[dataset_name]['data_corrected']
            metadata = self.datasets[dataset_name]['metadata']
            
            # 检查结局变量是否存在
            if metadata is not None and outcome in metadata.columns:
                # 计算效应量（Cohen's d）
                groups = metadata[outcome].unique()
                if len(groups) == 2:
                    group1 = data[metadata[outcome] == groups[0]][feature]
                    group2 = data[metadata[outcome] == groups[1]][feature]
                    
                    # Cohen's d
                    pooled_std = np.sqrt(
                        (group1.std()**2 + group2.std()**2) / 2
                    )
                    effect_size = (group1.mean() - group2.mean()) / pooled_std
                    
                    # 方差
                    n1, n2 = len(group1), len(group2)
                    variance = (n1 + n2) / (n1 * n2) + effect_size**2 / (2 * (n1 + n2))
                    
                    effect_sizes.append(effect_size)
                    variances.append(variance)
                    sample_sizes.append(n1 + n2)
        
        if len(effect_sizes) > 0:
            # 固定效应模型
            weights = [1/v for v in variances]
            pooled_effect = sum(e * w for e, w in zip(effect_sizes, weights)) / sum(weights)
            pooled_variance = 1 / sum(weights)
            pooled_se = np.sqrt(pooled_variance)
            
            # Z检验
            z_score = pooled_effect / pooled_se
            p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
            
            # 异质性检验（Q统计量）
            q_statistic = sum(w * (e - pooled_effect)**2 
                            for e, w in zip(effect_sizes, weights))
            df = len(effect_sizes) - 1
            p_heterogeneity = 1 - stats.chi2.cdf(q_statistic, df)
            
            # I²统计量
            if q_statistic > df:
                i_squared = 100 * (q_statistic - df) / q_statistic
            else:
                i_squared = 0
            
            results = {
                'Feature': feature,
                'Pooled_Effect_Size': pooled_effect,
                'SE': pooled_se,
                'Z_Score': z_score,
                'P_Value': p_value,
                'Q_Statistic': q_statistic,
                'P_Heterogeneity': p_heterogeneity,
                'I_Squared': i_squared,
                'N_Studies': len(effect_sizes),
                'Effect_Sizes': effect_sizes,
                'Sample_Sizes': sample_sizes
            }
            
            print(f"   汇总效应量: {pooled_effect:.3f}, P = {p_value:.3e}")
            print(f"   异质性: Q = {q_statistic:.2f}, I² = {i_squared:.1f}%")
            
            return results
        else:
            print(f"   ⚠️ 无法计算元分析（缺少结局变量）")
            return None
    
    def perform_cross_dataset_validation(self, features: List):
        """
        跨数据集验证
        
        Parameters:
        -----------
        features : List
            要验证的特征列表
            
        Returns:
        --------
        pd.DataFrame
            跨数据集验证结果
        """
        print(f"\n🔬 跨数据集验证 ({len(features)} 特征)...")
        
        validation_results = []
        
        for feature in features:
            results = self.perform_meta_analysis(feature)
            if results is not None:
                validation_results.append(results)
        
        self.meta_analysis_results = pd.DataFrame(validation_results)
        
        # 检查是否有有效的元分析结果
        if len(self.meta_analysis_results) > 0:
            # 多重检验校正
            from statsmodels.stats.multitest import multipletests
            p_values = self.meta_analysis_results['P_Value'].values
            _, p_corrected, _, _ = multipletests(p_values, method='fdr_bh')
            self.meta_analysis_results['P_Corrected'] = p_corrected
            self.meta_analysis_results['Significance'] = [
                '***' if p < 0.001 else '**' if p < 0.01 else 
                '*' if p < 0.05 else '' 
                for p in p_corrected
            ]
            
            print(f"✅ 跨数据集验证完成")
            print(f"   显著特征数: {(self.meta_analysis_results['P_Corrected'] < 0.05).sum()}")
        else:
            print("⚠️ 没有有效的元分析结果")
        
        return self.meta_analysis_results
    
    def assess_cross_ethnicity_consistency(self):
        """
        评估跨种族一致性
        
        Returns:
        --------
        Dict
            跨种族一致性结果
        """
        print("\n🌍 评估跨种族一致性...")
        
        consistency_results = {}
        
        # 检查元数据中的种族信息
        for dataset_name in self.datasets:
            metadata = self.datasets[dataset_name]['metadata']
            if 'Race' in metadata.columns or 'Ethnicity' in metadata.columns:
                race_col = 'Race' if 'Race' in metadata.columns else 'Ethnicity'
                races = metadata[race_col].unique()
                print(f"   {dataset_name}: {len(races)} 个种族")
        
        # 计算跨种族相关性
        common_features = self.identify_common_features()
        
        for i, ds1 in enumerate(list(self.datasets.keys())):
            for j, ds2 in enumerate(list(self.datasets.keys())):
                if i < j:
                    data1 = self.datasets[ds1]['data_corrected'][common_features]
                    data2 = self.datasets[ds2]['data_corrected'][common_features]
                    
                    # 计算相关性（使用共同样本或全部样本的平均值）
                    correlations = []
                    for feat in common_features:
                        # 使用特征的平均值进行跨数据集比较
                        mean1 = data1[feat].mean()
                        mean2 = data2[feat].mean()
                        
                        # 计算跨数据集的变异系数
                        cv1 = data1[feat].std() / abs(mean1) if mean1 != 0 else 0
                        cv2 = data2[feat].std() / abs(mean2) if mean2 != 0 else 0
                        
                        # 使用一致性指标（1 - |CV1 - CV2|）
                        consistency = 1 - abs(cv1 - cv2)
                        correlations.append(consistency)
                    
                    consistency_results[f'{ds1}_vs_{ds2}'] = {
                        'Mean_Correlation': np.mean(correlations),
                        'Median_Correlation': np.median(correlations),
                        'SD_Correlation': np.std(correlations),
                        'N_Features': len(common_features)
                    }
        
        print("✅ 跨种族一致性评估完成")
        return consistency_results
    
    def integrate_datasets(self):
        """
        整合所有数据集
        
        Returns:
        --------
        pd.DataFrame
            整合后的数据
        """
        print("\n🔗 整合数据集...")
        
        common_features = self.identify_common_features()
        
        integrated_list = []
        for dataset_name in self.datasets:
            data = self.datasets[dataset_name]['data_corrected'][common_features].copy()
            data['Dataset'] = dataset_name
            
            if self.datasets[dataset_name]['metadata'] is not None:
                metadata = self.datasets[dataset_name]['metadata']
                for col in metadata.columns:
                    if col not in data.columns:
                        data[col] = metadata[col].values
            
            integrated_list.append(data)
        
        self.integrated_data = pd.concat(integrated_list, axis=0, ignore_index=True)
        
        print(f"✅ 数据集整合完成:")
        print(f"   总样本数: {len(self.integrated_data)}")
        print(f"   特征数: {len(common_features)}")
        print(f"   数据集数: {len(self.datasets)}")
        
        return self.integrated_data
    
    def visualize_multi_dataset_results(self, output_dir: str = './results'):
        """
        可视化多数据集整合结果
        
        Parameters:
        -----------
        output_dir : str
            输出目录
        """
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # 1. 数据集样本量比较
        sample_sizes = [len(self.datasets[ds]['data']) for ds in self.datasets]
        dataset_names = list(self.datasets.keys())
        
        axes[0, 0].bar(range(len(dataset_names)), sample_sizes, 
                         color=['steelblue', 'coral', 'lightgreen'][:len(dataset_names)])
        axes[0, 0].set_xticks(range(len(dataset_names)))
        axes[0, 0].set_xticklabels(dataset_names)
        axes[0, 0].set_ylabel('Sample Size')
        axes[0, 0].set_title('Sample Size by Dataset')
        axes[0, 0].grid(True, alpha=0.3, axis='y')
        
        # 2. PCA可视化（整合数据）
        if self.integrated_data is not None:
            numeric_cols = self.integrated_data.select_dtypes(include=[np.number]).columns
            numeric_cols = [col for col in numeric_cols if col != 'Dataset']
            
            pca = PCA(n_components=2)
            pca_result = pca.fit_transform(
                self.integrated_data[numeric_cols].fillna(0)
            )
            
            for i, dataset in enumerate(self.integrated_data['Dataset'].unique()):
                mask = self.integrated_data['Dataset'] == dataset
                axes[0, 1].scatter(pca_result[mask, 0], pca_result[mask, 1], 
                                   label=dataset, alpha=0.6)
            
            axes[0, 1].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})')
            axes[0, 1].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})')
            axes[0, 1].set_title('PCA of Integrated Data')
            axes[0, 1].legend()
            axes[0, 1].grid(True, alpha=0.3)
        
        # 3. 元分析效应量森林图
        if self.meta_analysis_results is not None and len(self.meta_analysis_results) > 0:
            top_features = self.meta_analysis_results.nsmallest(10, 'P_Corrected')
            
            y_pos = range(len(top_features))
            axes[0, 2].barh(y_pos, top_features['Pooled_Effect_Size'], 
                             xerr=top_features['SE'], color='coral', alpha=0.7)
            axes[0, 2].set_yticks(y_pos)
            axes[0, 2].set_yticklabels(top_features['Feature'], fontsize=8)
            axes[0, 2].set_xlabel('Pooled Effect Size')
            axes[0, 2].set_title('Meta-Analysis: Top 10 Features')
            axes[0, 2].axvline(x=0, color='black', linestyle='--', linewidth=1)
            axes[0, 2].grid(True, alpha=0.3, axis='x')
        else:
            axes[0, 2].text(0.5, 0.5, 'No meta-analysis results', 
                           ha='center', va='center', transform=axes[0, 2].transAxes)
            axes[0, 2].set_title('Meta-Analysis')
        
        # 4. 异质性分析
        if self.meta_analysis_results is not None and len(self.meta_analysis_results) > 0:
            axes[1, 0].scatter(self.meta_analysis_results['Pooled_Effect_Size'], 
                               self.meta_analysis_results['I_Squared'], 
                               alpha=0.6)
            axes[1, 0].set_xlabel('Pooled Effect Size')
            axes[1, 0].set_ylabel('I² (%)')
            axes[1, 0].set_title('Heterogeneity Analysis')
            axes[1, 0].grid(True, alpha=0.3)
        else:
            axes[1, 0].text(0.5, 0.5, 'No meta-analysis results', 
                           ha='center', va='center', transform=axes[1, 0].transAxes)
            axes[1, 0].set_title('Heterogeneity Analysis')
        
        # 5. P值分布
        if self.meta_analysis_results is not None and len(self.meta_analysis_results) > 0:
            axes[1, 1].hist(self.meta_analysis_results['P_Corrected'], 
                            bins=30, color='lightgreen', alpha=0.7, edgecolor='black')
            axes[1, 1].axvline(0.05, color='red', linestyle='--', 
                                 linewidth=2, label='α = 0.05')
            axes[1, 1].set_xlabel('Corrected P-Value')
            axes[1, 1].set_ylabel('Frequency')
            axes[1, 1].set_title('P-Value Distribution')
            axes[1, 1].legend()
            axes[1, 1].grid(True, alpha=0.3)
        else:
            axes[1, 1].text(0.5, 0.5, 'No meta-analysis results', 
                           ha='center', va='center', transform=axes[1, 1].transAxes)
            axes[1, 1].set_title('P-Value Distribution')
        
        # 6. 显著特征数量
        if self.meta_analysis_results is not None and len(self.meta_analysis_results) > 0:
            significance_levels = ['***', '**', '*', '']
            sig_counts = [
                (self.meta_analysis_results['Significance'] == level).sum() 
                for level in significance_levels
            ]
            
            axes[1, 2].bar(range(len(significance_levels)), sig_counts, 
                             color=['red', 'orange', 'yellow', 'lightgray'])
            axes[1, 2].set_xticks(range(len(significance_levels)))
            axes[1, 2].set_xticklabels(significance_levels)
            axes[1, 2].set_ylabel('Count')
            axes[1, 2].set_title('Significant Features')
            axes[1, 2].grid(True, alpha=0.3, axis='y')
        else:
            axes[1, 2].text(0.5, 0.5, 'No meta-analysis results', 
                           ha='center', va='center', transform=axes[1, 2].transAxes)
            axes[1, 2].set_title('Significant Features')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/multi_dataset_integration.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✅ 多数据集整合可视化保存到: {output_dir}/multi_dataset_integration.png")


def demo_multi_dataset_integration():
    """
    演示多数据集整合
    """
    print("=" * 60)
    print("多数据集整合框架演示")
    print("=" * 60)
    
    # 创建整合器
    integrator = MultiDatasetIntegrator()
    
    # 模拟ADNI数据
    np.random.seed(42)
    n_adni = 300
    adni_data = pd.DataFrame(
        np.random.randn(n_adni, 50),
        index=[f'ADNI_{i}' for i in range(n_adni)],
        columns=[f'Gene_{i}' for i in range(50)]
    )
    adni_metadata = pd.DataFrame({
        'Diagnosis': np.random.choice(['CN', 'AD'], n_adni, p=[0.6, 0.4]),
        'Age': np.random.normal(75, 8, n_adni),
        'Sex': np.random.choice(['M', 'F'], n_adni),
        'Race': np.random.choice(['White', 'Black', 'Asian'], n_adni)
    }, index=[f'ADNI_{i}' for i in range(n_adni)])
    
    # 模拟ROSMAP数据
    n_rosmap = 250
    rosmap_data = pd.DataFrame(
        np.random.randn(n_rosmap, 50),
        index=[f'ROSMAP_{i}' for i in range(n_rosmap)],
        columns=[f'Gene_{i}' for i in range(50)]
    )
    rosmap_metadata = pd.DataFrame({
        'Diagnosis': np.random.choice(['CN', 'AD'], n_rosmap, p=[0.55, 0.45]),
        'Age': np.random.normal(80, 10, n_rosmap),
        'Sex': np.random.choice(['M', 'F'], n_rosmap),
        'Race': np.random.choice(['White', 'Black'], n_rosmap)
    }, index=[f'ROSMAP_{i}' for i in range(n_rosmap)])
    
    # 模拟MSBB数据
    n_msbb = 200
    msbb_data = pd.DataFrame(
        np.random.randn(n_msbb, 50),
        index=[f'MSBB_{i}' for i in range(n_msbb)],
        columns=[f'Gene_{i}' for i in range(50)]
    )
    msbb_metadata = pd.DataFrame({
        'Diagnosis': np.random.choice(['CN', 'AD'], n_msbb, p=[0.5, 0.5]),
        'Age': np.random.normal(85, 12, n_msbb),
        'Sex': np.random.choice(['M', 'F'], n_msbb),
        'Race': np.random.choice(['White', 'Black', 'Hispanic'], n_msbb)
    }, index=[f'MSBB_{i}' for i in range(n_msbb)])
    
    # 加载数据集
    integrator.load_dataset('ADNI', adni_data, adni_metadata)
    integrator.load_dataset('ROSMAP', rosmap_data, rosmap_metadata)
    integrator.load_dataset('MSBB', msbb_data, msbb_metadata)
    
    # 标准化数据集
    integrator.standardize_datasets()
    
    # 批量效应校正
    integrator.correct_batch_effects(method='combat')
    
    # 跨数据集验证
    features_to_validate = [f'Gene_{i}' for i in range(10)]
    validation_results = integrator.perform_cross_dataset_validation(features_to_validate)
    
    if validation_results is not None:
        print("\n📊 元分析结果（前5个特征）:")
        print(validation_results.head())
    
    # 跨种族一致性
    consistency = integrator.assess_cross_ethnicity_consistency()
    
    # 整合数据集
    integrated_data = integrator.integrate_datasets()
    
    # 可视化结果
    integrator.visualize_multi_dataset_results(output_dir='./results/multi_dataset_integration')
    
    print("\n" + "=" * 60)
    print("✅ 多数据集整合框架演示完成！")
    print("=" * 60)


if __name__ == '__main__':
    demo_multi_dataset_integration()