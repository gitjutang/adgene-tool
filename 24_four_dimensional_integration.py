"""
Four-Dimensional Integration Analysis (Spatial-Temporal-Cellular-Molecular)
四维整合分析（空间-时间-细胞-分子）

创新点：
1. 空间维度：空间转录组学数据整合
2. 时间维度：纵向时间序列分析
3. 细胞维度：单细胞分辨率分析
4. 分子维度：多组学数据整合
5. 四维联合建模和可视化

作者：[Author Names]
日期：2026-02-03
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA, NMF
from sklearn.cluster import KMeans, DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')


class FourDimensionalIntegration:
    """
    四维整合分析类
    
    功能：
    - 空间数据分析
    - 时间序列分析
    - 单细胞分析
    - 多组学整合
    - 四维联合建模
    """
    
    def __init__(self):
        self.spatial_data = None
        self.temporal_data = None
        self.cellular_data = None
        self.molecular_data = None
        self.integrated_data = None
        self.spatial_clusters = None
        self.temporal_patterns = None
        self.cell_types = None
        self.molecular_modules = None
        
    def load_spatial_data(self, data: pd.DataFrame, coordinates: pd.DataFrame):
        """
        加载空间转录组学数据
        
        Parameters:
        -----------
        data : pd.DataFrame
            基因表达矩阵 (spots × genes)
        coordinates : pd.DataFrame
            空间坐标 (x, y)
        """
        print("📍 加载空间转录组学数据...")
        self.spatial_data = {
            'expression': data,
            'coordinates': coordinates
        }
        print(f"   空间点数: {data.shape[0]}")
        print(f"   基因数: {data.shape[1]}")
        
    def load_temporal_data(self, data: pd.DataFrame, time_points: list):
        """
        加载时间序列数据
        
        Parameters:
        -----------
        data : pd.DataFrame
            时间序列数据 (samples × features)
        time_points : list
            时间点
        """
        print("⏰ 加载时间序列数据...")
        self.temporal_data = {
            'data': data,
            'time_points': time_points
        }
        print(f"   样本数: {data.shape[0]}")
        print(f"   特征数: {data.shape[1]}")
        print(f"   时间点数: {len(time_points)}")
        
    def load_cellular_data(self, data: pd.DataFrame, cell_metadata: pd.DataFrame):
        """
        加载单细胞数据
        
        Parameters:
        -----------
        data : pd.DataFrame
            单细胞表达矩阵 (cells × genes)
        cell_metadata : pd.DataFrame
            细胞元数据
        """
        print("🔬 加载单细胞数据...")
        self.cellular_data = {
            'expression': data,
            'metadata': cell_metadata
        }
        print(f"   细胞数: {data.shape[0]}")
        print(f"   基因数: {data.shape[1]}")
        
    def load_molecular_data(self, transcriptomics: pd.DataFrame = None, 
                           proteomics: pd.DataFrame = None,
                           metabolomics: pd.DataFrame = None):
        """
        加载多组学数据
        
        Parameters:
        -----------
        transcriptomics : pd.DataFrame
            转录组数据
        proteomics : pd.DataFrame
            蛋白质组数据
        metabolomics : pd.DataFrame
            代谢组数据
        """
        print("🧬 加载多组学数据...")
        self.molecular_data = {}
        
        if transcriptomics is not None:
            self.molecular_data['transcriptomics'] = transcriptomics
            print(f"   转录组: {transcriptomics.shape[0]} 样本 × {transcriptomics.shape[1]} 基因")
            
        if proteomics is not None:
            self.molecular_data['proteomics'] = proteomics
            print(f"   蛋白质组: {proteomics.shape[0]} 样本 × {proteomics.shape[1]} 蛋白质")
            
        if metabolomics is not None:
            self.molecular_data['metabolomics'] = metabolomics
            print(f"   代谢组: {metabolomics.shape[0]} 样本 × {metabolomics.shape[1]} 代谢物")
    
    def analyze_spatial_patterns(self, n_clusters: int = 5):
        """
        分析空间模式
        
        Parameters:
        -----------
        n_clusters : int
            空间聚类数
            
        Returns:
        --------
        Dict
            空间分析结果
        """
        print("\n📍 分析空间模式...")
        
        if self.spatial_data is None:
            raise ValueError("请先加载空间数据")
        
        expression = self.spatial_data['expression']
        coordinates = self.spatial_data['coordinates']
        
        # 空间聚类
        scaler = StandardScaler()
        expression_scaled = scaler.fit_transform(expression)
        
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        self.spatial_clusters = kmeans.fit_predict(expression_scaled)
        
        # 计算空间自相关
        distances = pdist(coordinates.values)
        distance_matrix = squareform(distances)
        
        # 计算Moran's I（简化版）
        moran_i = self._calculate_morans_i(expression_scaled, distance_matrix)
        
        # 识别空间差异基因
        spatial_genes = self._identify_spatial_genes(expression, self.spatial_clusters)
        
        results = {
            'clusters': self.spatial_clusters,
            'cluster_centers': kmeans.cluster_centers_,
            'morans_i': moran_i,
            'spatial_genes': spatial_genes
        }
        
        print(f"   识别了 {len(spatial_genes)} 个空间差异基因")
        print(f"   平均Moran's I: {np.mean(moran_i):.4f}")
        
        return results
    
    def _calculate_morans_i(self, data: np.ndarray, distance_matrix: np.ndarray):
        """
        计算Moran's I（空间自相关）
        """
        moran_i_values = []
        
        for i in range(data.shape[1]):
            z = data[:, i] - np.mean(data[:, i])
            w = 1 / (distance_matrix + 1e-6)
            np.fill_diagonal(w, 0)
            w = w / np.sum(w)
            
            numerator = np.sum(w * np.outer(z, z))
            denominator = np.sum(z ** 2)
            moran_i = numerator / denominator if denominator != 0 else 0
            moran_i_values.append(moran_i)
        
        return moran_i_values
    
    def _identify_spatial_genes(self, expression: pd.DataFrame, clusters: np.ndarray, 
                              top_n: int = 50):
        """
        识别空间差异基因
        """
        spatial_genes = []
        
        for gene in expression.columns:
            f_stat, p_value = stats.f_oneway(
                *[expression[clusters == i][gene] for i in np.unique(clusters)]
            )
            spatial_genes.append({
                'Gene': gene,
                'F_statistic': f_stat,
                'P_value': p_value
            })
        
        spatial_df = pd.DataFrame(spatial_genes)
        spatial_df = spatial_df.sort_values('F_statistic', ascending=False).head(top_n)
        
        return spatial_df
    
    def analyze_temporal_patterns(self, method: str = 'trend'):
        """
        分析时间模式
        
        Parameters:
        -----------
        method : str
            分析方法 ('trend', 'periodicity', 'change_point')
            
        Returns:
        --------
        Dict
            时间分析结果
        """
        print("\n⏰ 分析时间模式...")
        
        if self.temporal_data is None:
            raise ValueError("请先加载时间序列数据")
        
        data = self.temporal_data['data']
        time_points = self.temporal_data['time_points']
        
        results = {}
        
        if method == 'trend':
            # 趋势分析
            trends = {}
            for feature in data.columns:
                values = data[feature].values
                slope, intercept, r_value, p_value, std_err = stats.linregress(
                    range(len(values)), values
                )
                trends[feature] = {
                    'slope': slope,
                    'intercept': intercept,
                    'r_squared': r_value ** 2,
                    'p_value': p_value
                }
            
            results['trends'] = pd.DataFrame(trends).T
            
        elif method == 'periodicity':
            # 周期性分析（简化版）
            periodicities = {}
            for feature in data.columns:
                values = data[feature].values
                fft = np.fft.fft(values)
                power = np.abs(fft) ** 2
                freq = np.fft.fftfreq(len(values))
                
                dominant_freq = freq[np.argmax(power[1:len(power)//2]) + 1]
                periodicities[feature] = dominant_freq
            
            results['periodicities'] = periodicities
            
        elif method == 'change_point':
            # 变点检测（简化版）
            change_points = {}
            for feature in data.columns:
                values = data[feature].values
                max_change = 0
                change_point = 0
                
                for i in range(1, len(values) - 1):
                    before = values[:i]
                    after = values[i:]
                    
                    if len(before) > 0 and len(after) > 0:
                        change = abs(np.mean(before) - np.mean(after))
                        if change > max_change:
                            max_change = change
                            change_point = i
                
                change_points[feature] = {
                    'change_point': change_point,
                    'magnitude': max_change
                }
            
            results['change_points'] = pd.DataFrame(change_points).T
        
        self.temporal_patterns = results
        print(f"   时间模式分析完成")
        
        return results
    
    def analyze_cellular_types(self, n_clusters: int = 10):
        """
        分析细胞类型
        
        Parameters:
        -----------
        n_clusters : int
            细胞类型数
            
        Returns:
        --------
        Dict
            细胞类型分析结果
        """
        print("\n🔬 分析细胞类型...")
        
        if self.cellular_data is None:
            raise ValueError("请先加载单细胞数据")
        
        expression = self.cellular_data['expression']
        
        # 标准化
        scaler = StandardScaler()
        expression_scaled = scaler.fit_transform(expression)
        
        # 降维
        n_components_pca = min(50, min(expression_scaled.shape[0], expression_scaled.shape[1]) - 1)
        pca = PCA(n_components=n_components_pca, random_state=42)
        expression_pca = pca.fit_transform(expression_scaled)
        
        # 聚类
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        self.cell_types = kmeans.fit_predict(expression_pca)
        
        # 计算轮廓系数
        silhouette = silhouette_score(expression_pca, self.cell_types)
        
        # 识别细胞类型标记基因
        marker_genes = self._identify_marker_genes(expression, self.cell_types)
        
        results = {
            'cell_types': self.cell_types,
            'cell_type_centers': kmeans.cluster_centers_,
            'silhouette_score': silhouette,
            'marker_genes': marker_genes,
            'pca_variance': pca.explained_variance_ratio_
        }
        
        print(f"   识别了 {n_clusters} 个细胞类型")
        print(f"   轮廓系数: {silhouette:.4f}")
        print(f"   PCA解释方差: {pca.explained_variance_ratio_.sum():.2%}")
        
        return results
    
    def _identify_marker_genes(self, expression: pd.DataFrame, cell_types: np.ndarray, 
                              top_n: int = 10):
        """
        识别细胞类型标记基因
        """
        marker_genes = {}
        
        for cell_type in np.unique(cell_types):
            type_expression = expression[cell_types == cell_type]
            other_expression = expression[cell_types != cell_type]
            
            fold_changes = []
            p_values = []
            
            for gene in expression.columns:
                type_mean = type_expression[gene].mean()
                other_mean = other_expression[gene].mean()
                
                if other_mean == 0:
                    fold_change = np.inf
                else:
                    fold_change = type_mean / other_mean
                
                # t检验
                t_stat, p_value = stats.ttest_ind(
                    type_expression[gene],
                    other_expression[gene]
                )
                
                fold_changes.append(fold_change)
                p_values.append(p_value)
            
            marker_df = pd.DataFrame({
                'Gene': expression.columns,
                'Fold_Change': fold_changes,
                'P_Value': p_values
            })
            
            marker_df = marker_df.sort_values('Fold_Change', ascending=False).head(top_n)
            marker_genes[f'CellType_{cell_type}'] = marker_df
        
        return marker_genes
    
    def analyze_molecular_modules(self, n_modules: int = 10):
        """
        分析分子模块
        
        Parameters:
        -----------
        n_modules : int
            模块数
            
        Returns:
        --------
        Dict
            分子模块分析结果
        """
        print("\n🧬 分析分子模块...")
        
        if self.molecular_data is None:
            raise ValueError("请先加载多组学数据")
        
        results = {}
        
        for omics_type, data in self.molecular_data.items():
            print(f"   分析 {omics_type}...")
            
            # 标准化
            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(data)
            
            # 确保数据非负（NMF要求）
            data_nonneg = data_scaled - data_scaled.min() + 1e-6
            
            # NMF分解
            actual_n_modules = min(n_modules, min(data_nonneg.shape[0], data_nonneg.shape[1]))
            nmf = NMF(n_components=actual_n_modules, random_state=42, max_iter=200)
            W = nmf.fit_transform(data_nonneg)
            H = nmf.components_
            
            # 计算模块相关性
            module_corr = np.corrcoef(W.T)
            
            results[omics_type] = {
                'sample_scores': W,
                'module_features': H,
                'module_correlation': module_corr,
                'reconstruction_error': nmf.reconstruction_err_
            }
            
            print(f"      重构误差: {nmf.reconstruction_err_:.4f}")
        
        self.molecular_modules = results
        print("   分子模块分析完成")
        
        return results
    
    def integrate_four_dimensions(self, method: str = 'concatenation'):
        """
        四维整合
        
        Parameters:
        -----------
        method : str
            整合方法 ('concatenation', 'pca', 'nmf')
            
        Returns:
        --------
        pd.DataFrame
            整合后的数据
        """
        print("\n🔗 四维整合...")
        
        if self.spatial_data is None and self.temporal_data is None and \
           self.cellular_data is None and self.molecular_data is None:
            raise ValueError("请至少加载一种数据")
        
        # 收集所有数据
        all_data = []
        
        if self.spatial_data is not None:
            spatial_features = self.spatial_data['expression'].values
            all_data.append(('Spatial', spatial_features))
            
        if self.temporal_data is not None:
            temporal_features = self.temporal_data['data'].values
            all_data.append(('Temporal', temporal_features))
            
        if self.cellular_data is not None:
            cellular_features = self.cellular_data['expression'].values
            all_data.append(('Cellular', cellular_features))
            
        if self.molecular_data is not None:
            for omics_type, data in self.molecular_data.items():
                all_data.append((omics_type.capitalize(), data.values))
        
        # 整合方法
        if method == 'concatenation':
            # 简单拼接
            max_samples = max([d.shape[0] for _, d in all_data])
            integrated = np.zeros((max_samples, sum([d.shape[1] for _, d in all_data])))
            
            idx = 0
            for data_type, data in all_data:
                n_samples = data.shape[0]
                n_features = data.shape[1]
                integrated[:n_samples, idx:idx+n_features] = data
                idx += n_features
            
            feature_names = []
            for data_type, data in all_data:
                feature_names.extend([f'{data_type}_Feature_{i}' for i in range(data.shape[1])])
            
            self.integrated_data = pd.DataFrame(integrated, columns=feature_names)
            
        elif method == 'pca':
            # PCA降维后整合
            max_samples = max([d.shape[0] for _, d in all_data])
            n_components = 20
            
            integrated = np.zeros((max_samples, len(all_data) * n_components))
            
            for i, (data_type, data) in enumerate(all_data):
                actual_n_components = min(n_components, min(data.shape[0], data.shape[1]))
                pca = PCA(n_components=actual_n_components, random_state=42)
                data_pca = pca.fit_transform(data)
                integrated[:data_pca.shape[0], i*n_components:(i+1)*n_components] = np.pad(
                    data_pca, ((0, 0), (0, n_components - actual_n_components)), mode='constant'
                )
            
            feature_names = []
            for data_type, _ in all_data:
                feature_names.extend([f'{data_type}_PC_{i}' for i in range(n_components)])
            
            self.integrated_data = pd.DataFrame(integrated, columns=feature_names)
            
        elif method == 'nmf':
            # NMF分解后整合
            max_samples = max([d.shape[0] for _, d in all_data])
            n_components = 10
            
            integrated = np.zeros((max_samples, len(all_data) * n_components))
            
            for i, (data_type, data) in enumerate(all_data):
                # 确保数据非负
                data_nonneg = data - data.min() + 1e-6
                actual_n_components = min(n_components, min(data_nonneg.shape[0], data_nonneg.shape[1]))
                nmf = NMF(n_components=actual_n_components, random_state=42, max_iter=200)
                data_nmf = nmf.fit_transform(data_nonneg)
                integrated[:data_nmf.shape[0], i*n_components:(i+1)*n_components] = np.pad(
                    data_nmf, ((0, 0), (0, n_components - actual_n_components)), mode='constant'
                )
            
            feature_names = []
            for data_type, _ in all_data:
                feature_names.extend([f'{data_type}_Module_{i}' for i in range(n_components)])
            
            self.integrated_data = pd.DataFrame(integrated, columns=feature_names)
        
        print(f"   整合完成: {self.integrated_data.shape[0]} 样本 × {self.integrated_data.shape[1]} 特征")
        
        return self.integrated_data
    
    def visualize_four_dimensions(self, output_dir: str = './results/four_dimensional_integration'):
        """
        可视化四维整合结果
        
        Parameters:
        -----------
        output_dir : str
            输出目录
        """
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        fig = plt.figure(figsize=(24, 16))
        gs = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)
        
        # 1. 空间聚类可视化
        if self.spatial_data is not None and self.spatial_clusters is not None:
            ax = fig.add_subplot(gs[0, 0])
            coordinates = self.spatial_data['coordinates']
            scatter = ax.scatter(coordinates.iloc[:, 0], coordinates.iloc[:, 1], 
                               c=self.spatial_clusters, cmap='tab10', s=50, alpha=0.7)
            ax.set_xlabel('X Coordinate')
            ax.set_ylabel('Y Coordinate')
            ax.set_title('Spatial Clustering')
            plt.colorbar(scatter, ax=ax, label='Cluster')
        
        # 2. 空间差异基因热图
        if self.spatial_data is not None and self.spatial_clusters is not None:
            ax = fig.add_subplot(gs[0, 1])
            spatial_genes = self._identify_spatial_genes(
                self.spatial_data['expression'], self.spatial_clusters, top_n=20
            )
            sns.heatmap(self.spatial_data['expression'][spatial_genes['Gene'].head(10)].T,
                       cmap='viridis', ax=ax, cbar_kws={'label': 'Expression'})
            ax.set_title('Top Spatial Variable Genes')
        
        # 3. 时间趋势
        if self.temporal_data is not None and 'trends' in self.temporal_patterns:
            ax = fig.add_subplot(gs[0, 2])
            trends = self.temporal_patterns['trends']
            top_trends = trends.sort_values('r_squared', ascending=False).head(10)
            for feature in top_trends.index:
                slope = top_trends.loc[feature, 'slope']
                intercept = top_trends.loc[feature, 'intercept']
                x = np.arange(len(self.temporal_data['time_points']))
                y = slope * x + intercept
                ax.plot(x, y, label=feature, alpha=0.7)
            ax.set_xlabel('Time Point')
            ax.set_ylabel('Value')
            ax.set_title('Top Temporal Trends')
            ax.legend(fontsize=6)
        
        # 4. 时间序列热图
        if self.temporal_data is not None:
            ax = fig.add_subplot(gs[0, 3])
            top_features = self.temporal_data['data'].var().nlargest(10).index
            sns.heatmap(self.temporal_data['data'][top_features].T,
                       cmap='viridis', ax=ax, cbar_kws={'label': 'Value'})
            ax.set_title('Temporal Expression Heatmap')
        
        # 5. 细胞类型PCA
        if self.cellular_data is not None and self.cell_types is not None:
            ax = fig.add_subplot(gs[1, 0])
            expression = self.cellular_data['expression']
            scaler = StandardScaler()
            expression_scaled = scaler.fit_transform(expression)
            pca = PCA(n_components=2, random_state=42)
            expression_pca = pca.fit_transform(expression_scaled)
            scatter = ax.scatter(expression_pca[:, 0], expression_pca[:, 1],
                               c=self.cell_types, cmap='tab10', s=20, alpha=0.7)
            ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})')
            ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})')
            ax.set_title('Cell Type PCA')
            plt.colorbar(scatter, ax=ax, label='Cell Type')
        
        # 6. 细胞类型标记基因
        if self.cellular_data is not None and self.cell_types is not None:
            ax = fig.add_subplot(gs[1, 1])
            marker_genes = self._identify_marker_genes(
                self.cellular_data['expression'], self.cell_types, top_n=5
            )
            for i, (cell_type, markers) in enumerate(marker_genes.items()):
                if i < 3:
                    ax.barh([f'{cell_type}_{gene}' for gene in markers['Gene']],
                           markers['Fold_Change'], label=cell_type, alpha=0.7)
            ax.set_xlabel('Fold Change')
            ax.set_title('Cell Type Marker Genes')
            ax.legend()
        
        # 7. 细胞类型分布
        if self.cell_types is not None:
            ax = fig.add_subplot(gs[1, 2])
            unique, counts = np.unique(self.cell_types, return_counts=True)
            ax.bar(unique, counts, color='steelblue', alpha=0.7)
            ax.set_xlabel('Cell Type')
            ax.set_ylabel('Count')
            ax.set_title('Cell Type Distribution')
        
        # 8. 分子模块热图
        if self.molecular_modules is not None:
            ax = fig.add_subplot(gs[1, 3])
            for omics_type, modules in self.molecular_modules.items():
                if 'module_correlation' in modules:
                    sns.heatmap(modules['module_correlation'], 
                               cmap='coolwarm', center=0,
                               ax=ax, cbar_kws={'label': 'Correlation'})
                    ax.set_title(f'{omics_type} Module Correlation')
                    break
        
        # 9. 整合数据PCA
        if self.integrated_data is not None:
            ax = fig.add_subplot(gs[2, 0])
            scaler = StandardScaler()
            integrated_scaled = scaler.fit_transform(self.integrated_data)
            pca = PCA(n_components=2, random_state=42)
            integrated_pca = pca.fit_transform(integrated_scaled)
            ax.scatter(integrated_pca[:, 0], integrated_pca[:, 1], 
                      alpha=0.5, s=30)
            ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})')
            ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})')
            ax.set_title('Integrated Data PCA')
        
        # 10. 整合数据相关性热图
        if self.integrated_data is not None:
            ax = fig.add_subplot(gs[2, 1])
            corr_matrix = self.integrated_data.corr()
            # 只显示部分相关性以避免过于密集
            n_show = min(20, corr_matrix.shape[0])
            sns.heatmap(corr_matrix.iloc[:n_show, :n_show],
                       cmap='coolwarm', center=0, ax=ax,
                       cbar_kws={'label': 'Correlation'})
            ax.set_title('Integrated Data Correlation')
        
        # 11. 四维特征重要性
        if self.integrated_data is not None:
            ax = fig.add_subplot(gs[2, 2])
            feature_importance = self.integrated_data.var().nlargest(15)
            feature_importance.plot(kind='barh', ax=ax, color='steelblue')
            ax.set_xlabel('Variance')
            ax.set_title('Top Feature Importance')
        
        # 12. 四维聚类
        if self.integrated_data is not None:
            ax = fig.add_subplot(gs[2, 3])
            scaler = StandardScaler()
            integrated_scaled = scaler.fit_transform(self.integrated_data)
            kmeans = KMeans(n_clusters=5, random_state=42)
            clusters = kmeans.fit_predict(integrated_scaled)
            pca = PCA(n_components=2, random_state=42)
            integrated_pca = pca.fit_transform(integrated_scaled)
            scatter = ax.scatter(integrated_pca[:, 0], integrated_pca[:, 1],
                               c=clusters, cmap='tab10', s=30, alpha=0.7)
            ax.set_xlabel('PC1')
            ax.set_ylabel('PC2')
            ax.set_title('Four-Dimensional Clustering')
            plt.colorbar(scatter, ax=ax, label='Cluster')
        
        # 13. 空间-时间整合
        if self.spatial_data is not None and self.temporal_data is not None:
            ax = fig.add_subplot(gs[3, 0])
            spatial_mean = self.spatial_data['expression'].mean(axis=1)
            temporal_mean = self.temporal_data['data'].mean(axis=1)
            min_len = min(len(spatial_mean), len(temporal_mean))
            ax.scatter(spatial_mean[:min_len], temporal_mean[:min_len], 
                      alpha=0.5, s=30)
            ax.set_xlabel('Spatial Mean Expression')
            ax.set_ylabel('Temporal Mean Expression')
            ax.set_title('Spatial-Temporal Integration')
        
        # 14. 细胞-分子整合
        if self.cellular_data is not None and self.molecular_data is not None:
            ax = fig.add_subplot(gs[3, 1])
            cellular_mean = self.cellular_data['expression'].mean(axis=1)
            molecular_mean = list(self.molecular_data.values())[0].mean(axis=1)
            min_len = min(len(cellular_mean), len(molecular_mean))
            ax.scatter(cellular_mean[:min_len], molecular_mean[:min_len],
                      alpha=0.5, s=30)
            ax.set_xlabel('Cellular Mean Expression')
            ax.set_ylabel('Molecular Mean Expression')
            ax.set_title('Cellular-Molecular Integration')
        
        # 15. 四维网络图
        if self.integrated_data is not None:
            ax = fig.add_subplot(gs[3, 2])
            corr_matrix = self.integrated_data.corr()
            n_show = min(10, corr_matrix.shape[0])
            threshold = 0.5
            for i in range(n_show):
                for j in range(i+1, n_show):
                    if abs(corr_matrix.iloc[i, j]) > threshold:
                        ax.plot([i, j], [0, 0], 'k-', alpha=0.3, linewidth=abs(corr_matrix.iloc[i, j])*2)
            ax.scatter(range(n_show), [0]*n_show, s=100, c='red', zorder=5)
            ax.set_xlim(-0.5, n_show-0.5)
            ax.set_ylim(-1, 1)
            ax.set_yticks([])
            ax.set_title('Four-Dimensional Network')
        
        # 16. 四维统计摘要
        ax = fig.add_subplot(gs[3, 3])
        ax.axis('off')
        summary_text = "Four-Dimensional Integration Summary\n\n"
        
        if self.spatial_data is not None:
            summary_text += f"Spatial: {self.spatial_data['expression'].shape[0]} spots × {self.spatial_data['expression'].shape[1]} genes\n"
        
        if self.temporal_data is not None:
            summary_text += f"Temporal: {self.temporal_data['data'].shape[0]} samples × {self.temporal_data['data'].shape[1]} features\n"
        
        if self.cellular_data is not None:
            summary_text += f"Cellular: {self.cellular_data['expression'].shape[0]} cells × {self.cellular_data['expression'].shape[1]} genes\n"
        
        if self.molecular_data is not None:
            for omics_type, data in self.molecular_data.items():
                summary_text += f"{omics_type.capitalize()}: {data.shape[0]} samples × {data.shape[1]} features\n"
        
        if self.integrated_data is not None:
            summary_text += f"\nIntegrated: {self.integrated_data.shape[0]} samples × {self.integrated_data.shape[1]} features"
        
        ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top', family='monospace')
        
        plt.savefig(f'{output_dir}/four_dimensional_integration.png', 
                   dpi=300, bbox_inches='tight')
        print(f"✅ 四维整合可视化保存到: {output_dir}/four_dimensional_integration.png")
        
        plt.close()


def demo_four_dimensional_integration():
    """
    演示四维整合分析
    """
    print("=" * 60)
    print("四维整合分析（空间-时间-细胞-分子）")
    print("=" * 60)
    
    # 创建四维整合分析器
    integrator = FourDimensionalIntegration()
    
    # 创建模拟空间数据
    print("\n📍 创建模拟空间数据...")
    n_spots = 100
    n_genes = 50
    spatial_expression = pd.DataFrame(
        np.random.randn(n_spots, n_genes),
        columns=[f'Gene_{i}' for i in range(n_genes)]
    )
    spatial_coordinates = pd.DataFrame({
        'x': np.random.rand(n_spots) * 10,
        'y': np.random.rand(n_spots) * 10
    })
    integrator.load_spatial_data(spatial_expression, spatial_coordinates)
    
    # 创建模拟时间数据
    print("\n⏰ 创建模拟时间数据...")
    n_timepoints = 10
    n_features = 30
    temporal_data = pd.DataFrame(
        np.random.randn(n_timepoints, n_features),
        columns=[f'Feature_{i}' for i in range(n_features)]
    )
    time_points = list(range(n_timepoints))
    integrator.load_temporal_data(temporal_data, time_points)
    
    # 创建模拟单细胞数据
    print("\n🔬 创建模拟单细胞数据...")
    n_cells = 200
    sc_genes = 40
    cellular_expression = pd.DataFrame(
        np.random.randn(n_cells, sc_genes),
        columns=[f'SCGene_{i}' for i in range(sc_genes)]
    )
    cell_metadata = pd.DataFrame({
        'Cell_ID': [f'Cell_{i}' for i in range(n_cells)]
    })
    integrator.load_cellular_data(cellular_expression, cell_metadata)
    
    # 创建模拟多组学数据
    print("\n🧬 创建模拟多组学数据...")
    n_samples = 50
    transcriptomics = pd.DataFrame(
        np.random.randn(n_samples, 35),
        columns=[f'Transcript_{i}' for i in range(35)]
    )
    proteomics = pd.DataFrame(
        np.random.randn(n_samples, 25),
        columns=[f'Protein_{i}' for i in range(25)]
    )
    metabolomics = pd.DataFrame(
        np.random.randn(n_samples, 20),
        columns=[f'Metabolite_{i}' for i in range(20)]
    )
    integrator.load_molecular_data(transcriptomics, proteomics, metabolomics)
    
    # 分析空间模式
    spatial_results = integrator.analyze_spatial_patterns(n_clusters=5)
    
    # 分析时间模式
    temporal_results = integrator.analyze_temporal_patterns(method='trend')
    
    # 分析细胞类型
    cellular_results = integrator.analyze_cellular_types(n_clusters=8)
    
    # 分析分子模块
    molecular_results = integrator.analyze_molecular_modules(n_modules=8)
    
    # 四维整合
    integrated_data = integrator.integrate_four_dimensions(method='concatenation')
    
    # 可视化
    integrator.visualize_four_dimensions()
    
    print("\n" + "=" * 60)
    print("✅ 四维整合分析演示完成！")
    print("=" * 60)


if __name__ == '__main__':
    demo_four_dimensional_integration()
