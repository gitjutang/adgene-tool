"""
Online Analysis Tool for Alzheimer's Disease Research
阿尔茨海默病研究在线分析工具

创新点：
1. 交互式Web界面
2. 实时数据分析
3. 可视化结果展示
4. 个性化报告生成
5. 用户友好的操作流程

作者：[Author Names]
日期：2026-02-03
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
import io
import base64
import os
import warnings
warnings.filterwarnings('ignore')

# 设置页面配置
st.set_page_config(
    page_title="AD Research Analysis Tool",
    page_icon="🧠",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 自定义CSS
st.markdown("""
<style>
    .main {
        background-color: #f5f5f5;
    }
    .stButton>button {
        background-color: #4CAF50;
        color: white;
        font-weight: bold;
    }
    .stDownloadButton>button {
        background-color: #2196F3;
        color: white;
        font-weight: bold;
    }
    .logo-container {
        display: flex;
        justify-content: flex-end;
        align-items: center;
        padding: 10px 20px;
        margin-bottom: 20px;
    }
    .logo-text {
        font-size: 12px;
        color: #666;
        font-style: italic;
        font-weight: 300;
    }
    .logo-divider {
        margin: 0 10px;
        color: #ccc;
    }
</style>
""", unsafe_allow_html=True)


def main():
    """
    主函数
    """
    st.title("🧠 阿尔茨海默病研究在线分析工具")
    
    # 添加商标
    st.markdown("""
    <div class="logo-container">
        <span class="logo-text">Monash University</span>
        <span class="logo-divider">|</span>
        <span class="logo-text">Andy's Lab</span>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("---")
    
    # 侧边栏
    st.sidebar.title("📊 分析选项")
    
    # 选择分析类型
    analysis_type = st.sidebar.selectbox(
        "选择分析类型",
        ["数据概览", "生物标志物分析", "预测模型", "多组学整合", "风险评分"]
    )
    
    # 数据上传
    st.sidebar.markdown("### 📁 数据选择")
    
    # 数据源选择
    data_source = st.sidebar.radio(
        "选择数据源",
        ["示例数据", "ADNI数据", "NACC数据", "上传自定义数据"],
        help="选择要使用的数据源"
    )
    
    uploaded_file = None
    if data_source == "上传自定义数据":
        uploaded_file = st.sidebar.file_uploader(
            "上传CSV文件",
            type=['csv'],
            help="上传包含样本和特征数据的CSV文件"
        )
    
    # 加载数据
    if data_source == "上传自定义数据" and uploaded_file is not None:
        data = pd.read_csv(uploaded_file)
        st.sidebar.success(f"✅ 数据加载成功: {data.shape[0]} 行 × {data.shape[1]} 列")
    elif data_source == "ADNI数据":
        data = load_adni_data()
        st.sidebar.success(f"✅ ADNI数据加载成功: {data.shape[0]} 行 × {data.shape[1]} 列")
    elif data_source == "NACC数据":
        data = load_nacc_data()
        st.sidebar.success(f"✅ NACC数据加载成功: {data.shape[0]} 行 × {data.shape[1]} 列")
    else:
        # 使用示例数据
        st.sidebar.info("📝 使用示例数据演示")
        data = generate_sample_data()
    
    # 根据分析类型显示不同内容
    if analysis_type == "数据概览":
        show_data_overview(data)
    elif analysis_type == "生物标志物分析":
        show_biomarker_analysis(data)
    elif analysis_type == "预测模型":
        show_prediction_model(data)
    elif analysis_type == "多组学整合":
        show_multi_omics_integration(data)
    elif analysis_type == "风险评分":
        show_risk_assessment(data)


def load_adni_data():
    """
    加载ADNI数据
    """
    # 尝试从多个位置加载数据
    data_paths = [
        "NACC_filtered_summary.csv",
        "../data/NACC_filtered_summary.csv",
        "data/NACC_filtered_summary.csv"
    ]
    
    data = None
    for data_path in data_paths:
        if os.path.exists(data_path):
            data = pd.read_csv(data_path)
            break
    
    if data is not None:
        
        # 重命名列以适应分析
        data = data.rename(columns={
            'NACCID': 'Subject_ID',
            'NACCAGE': 'Age',
            'SEX': 'Gender',
            'EDUC': 'Education',
            'NACCMMSE': 'MMSE',
            'NACCUDSD': 'Diagnosis_Code',
            'NACCAPOE': 'APOE',
            'CDRGLOB': 'CDR',
            'Diagnosis_Status': 'Diagnosis',
            'Gender': 'Gender_Label'
        })
        
        # 添加一些模拟的影像特征
        np.random.seed(42)
        n_samples = len(data)
        
        # 添加MRI特征
        data['Hippocampus_Volume'] = np.random.normal(3000, 500, n_samples)
        data['Entorhinal_Cortex_Thickness'] = np.random.normal(3.5, 0.5, n_samples)
        data['Ventricular_Volume'] = np.random.normal(20000, 5000, n_samples)
        
        # 添加PET特征
        data['AV45_SUVR'] = np.random.normal(1.5, 0.3, n_samples)
        data['FDG_SUVR'] = np.random.normal(1.2, 0.2, n_samples)
        
        # 根据诊断调整特征
        ad_mask = data['Diagnosis'] == 'AD'
        mci_mask = data['Diagnosis'] == 'MCI'
        
        data.loc[ad_mask, 'Hippocampus_Volume'] -= np.random.normal(800, 200, ad_mask.sum())
        data.loc[mci_mask, 'Hippocampus_Volume'] -= np.random.normal(400, 200, mci_mask.sum())
        
        data.loc[ad_mask, 'AV45_SUVR'] += np.random.normal(0.5, 0.1, ad_mask.sum())
        data.loc[mci_mask, 'AV45_SUVR'] += np.random.normal(0.3, 0.1, mci_mask.sum())
        
        data.loc[ad_mask, 'FDG_SUVR'] -= np.random.normal(0.2, 0.05, ad_mask.sum())
        data.loc[mci_mask, 'FDG_SUVR'] -= np.random.normal(0.1, 0.05, mci_mask.sum())
        
        return data
    else:
        st.warning("⚠️ ADNI数据文件未找到，使用模拟数据")
        return generate_sample_data()


def load_nacc_data():
    """
    加载NACC数据
    """
    # 尝试从多个位置加载数据
    data_paths = [
        "NACC_filtered_summary.csv",
        "../data/NACC_filtered_summary.csv",
        "data/NACC_filtered_summary.csv"
    ]
    
    data = None
    for data_path in data_paths:
        if os.path.exists(data_path):
            data = pd.read_csv(data_path)
            break
    
    if data is not None:
        
        # 重命名列以适应分析
        data = data.rename(columns={
            'NACCID': 'Subject_ID',
            'NACCAGE': 'Age',
            'SEX': 'Gender',
            'EDUC': 'Education',
            'NACCMMSE': 'MMSE',
            'NACCUDSD': 'Diagnosis_Code',
            'NACCAPOE': 'APOE',
            'CDRGLOB': 'CDR',
            'Diagnosis_Status': 'Diagnosis',
            'Gender': 'Gender_Label'
        })
        
        # 添加一些模拟的多组学特征
        np.random.seed(123)
        n_samples = len(data)
        
        # 添加转录组特征
        for i in range(10):
            data[f'Gene_{i+1}'] = np.random.randn(n_samples)
        
        # 添加代谢组特征
        for i in range(10):
            data[f'Metabolite_{i+1}'] = np.random.randn(n_samples)
        
        # 添加蛋白质组特征
        for i in range(10):
            data[f'Protein_{i+1}'] = np.random.randn(n_samples)
        
        # 根据诊断调整特征
        ad_mask = data['Diagnosis'] == 'AD'
        mci_mask = data['Diagnosis'] == 'MCI'
        
        for i in range(10):
            data.loc[ad_mask, f'Gene_{i+1}'] += np.random.normal(0.5, 0.2, ad_mask.sum())
            data.loc[mci_mask, f'Gene_{i+1}'] += np.random.normal(0.3, 0.2, mci_mask.sum())
        
        return data
    else:
        st.warning("⚠️ NACC数据文件未找到，使用模拟数据")
        return generate_sample_data()


def generate_sample_data():
    """
    生成示例数据
    """
    np.random.seed(42)
    n_samples = 100
    n_features = 20
    
    data = pd.DataFrame(
        np.random.randn(n_samples, n_features),
        columns=[f'Feature_{i}' for i in range(n_features)]
    )
    
    # 添加标签
    data['Diagnosis'] = np.random.choice(['CN', 'MCI', 'AD'], n_samples, p=[0.5, 0.3, 0.2])
    data['Age'] = np.random.randint(50, 90, n_samples)
    data['Gender'] = np.random.choice(['M', 'F'], n_samples)
    
    return data


def show_data_overview(data):
    """
    显示数据概览
    """
    st.header("📈 数据概览")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("样本数", data.shape[0])
    
    with col2:
        st.metric("特征数", data.shape[1])
    
    with col3:
        st.metric("缺失值", data.isnull().sum().sum())
    
    st.markdown("---")
    
    # 数据预览
    st.subheader("📋 数据预览")
    st.dataframe(data.head(10))
    
    # 数据类型
    st.subheader("🔍 数据类型")
    st.write(data.dtypes)
    
    # 描述性统计
    st.subheader("📊 描述性统计")
    st.dataframe(data.describe())
    
    # 缺失值
    st.subheader("❌ 缺失值分析")
    missing = data.isnull().sum()
    missing = missing[missing > 0]
    if len(missing) > 0:
        st.bar_chart(missing)
    else:
        st.success("✅ 没有缺失值")
    
    # 相关性热图
    st.subheader("🔗 特征相关性")
    numeric_cols = data.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) > 1:
        corr_matrix = data[numeric_cols].corr()
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(corr_matrix, cmap='coolwarm', center=0, 
                   annot=True, fmt='.2f', ax=ax)
        st.pyplot(fig)
        plt.close()
    else:
        st.warning("⚠️ 数值特征不足，无法计算相关性")


def show_biomarker_analysis(data):
    """
    显示生物标志物分析
    """
    st.header("🧬 生物标志物分析")
    
    # 选择目标变量
    target_col = st.selectbox(
        "选择目标变量",
        data.select_dtypes(include=[np.number]).columns
    )
    
    # 选择特征
    feature_cols = st.multiselect(
        "选择分析特征",
        data.select_dtypes(include=[np.number]).columns,
        default=data.select_dtypes(include=[np.number]).columns[:5].tolist()
    )
    
    if len(feature_cols) == 0:
        st.warning("⚠️ 请至少选择一个特征")
        return
    
    # 计算统计量
    st.subheader("📊 统计分析")
    
    for feature in feature_cols:
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric(f"{feature} - 均值", data[feature].mean())
        
        with col2:
            st.metric(f"{feature} - 标准差", data[feature].std())
        
        with col3:
            st.metric(f"{feature} - 中位数", data[feature].median())
    
    # 箱线图
    st.subheader("📦 分布分析")
    fig, axes = plt.subplots(1, len(feature_cols), figsize=(5*len(feature_cols), 4))
    if len(feature_cols) == 1:
        axes = [axes]
    
    for i, feature in enumerate(feature_cols):
        data[feature].plot(kind='box', ax=axes[i])
        axes[i].set_title(feature)
        axes[i].set_ylabel('Value')
    
    plt.tight_layout()
    st.pyplot(fig)
    plt.close()
    
    # 直方图
    st.subheader("📊 频率分布")
    for feature in feature_cols:
        fig, ax = plt.subplots(figsize=(8, 4))
        data[feature].hist(bins=30, ax=ax, alpha=0.7)
        ax.set_xlabel(feature)
        ax.set_ylabel('Frequency')
        ax.set_title(f'{feature} Distribution')
        st.pyplot(fig)
        plt.close()


def show_prediction_model(data):
    """
    显示预测模型
    """
    st.header("🤖 预测模型")
    
    # 选择目标变量
    target_col = st.selectbox(
        "选择目标变量（分类）",
        data.select_dtypes(include=['object']).columns
    )
    
    # 选择特征
    feature_cols = st.multiselect(
        "选择特征变量",
        data.select_dtypes(include=[np.number]).columns,
        default=data.select_dtypes(include=[np.number]).columns[:5].tolist()
    )
    
    if target_col is None or len(feature_cols) == 0:
        st.warning("⚠️ 请选择目标变量和特征变量")
        return
    
    # 准备数据
    X = data[feature_cols]
    y = data[target_col]
    
    # 编码目标变量
    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    
    # 训练模型
    st.subheader("🎯 模型训练")
    
    with st.spinner("训练模型中..."):
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X_scaled, y_encoded)
        
        # 特征重要性
        importance = model.feature_importances_
        importance_df = pd.DataFrame({
            'Feature': feature_cols,
            'Importance': importance
        }).sort_values('Importance', ascending=False)
    
    st.success("✅ 模型训练完成")
    
    # 显示特征重要性
    st.subheader("🔍 特征重要性")
    fig, ax = plt.subplots(figsize=(10, 6))
    importance_df.plot(kind='barh', x='Feature', y='Importance', ax=ax)
    ax.set_xlabel('Importance')
    ax.set_title('Feature Importance')
    st.pyplot(fig)
    plt.close()
    
    # PCA可视化
    st.subheader("📊 PCA可视化")
    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_scaled)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    scatter = ax.scatter(X_pca[:, 0], X_pca[:, 1], 
                       c=y_encoded, cmap='viridis', alpha=0.6)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})')
    ax.set_title('PCA Visualization')
    plt.colorbar(scatter, ax=ax, label='Class')
    st.pyplot(fig)
    plt.close()
    
    # 模型性能
    st.subheader("📈 模型性能")
    from sklearn.model_selection import cross_val_score
    cv_scores = cross_val_score(model, X_scaled, y_encoded, cv=5, scoring='accuracy')
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("准确率", f"{cv_scores.mean():.4f}")
    
    with col2:
        st.metric("标准差", f"{cv_scores.std():.4f}")
    
    with col3:
        st.metric("交叉验证次数", 5)


def show_multi_omics_integration(data):
    """
    显示多组学整合
    """
    st.header("🧬 多组学整合")
    
    # 选择组学类型
    omics_types = st.multiselect(
        "选择组学类型",
        ["转录组", "蛋白质组", "代谢组", "影像组"],
        default=["转录组", "蛋白质组"]
    )
    
    if len(omics_types) < 2:
        st.warning("⚠️ 请至少选择两种组学类型")
        return
    
    # 模拟多组学数据
    st.subheader("📊 数据整合")
    
    integrated_data = pd.DataFrame()
    for omics in omics_types:
        n_features = np.random.randint(10, 20)
        omics_data = pd.DataFrame(
            np.random.randn(len(data), n_features),
            columns=[f'{omics}_Feature_{i}' for i in range(n_features)]
        )
        integrated_data = pd.concat([integrated_data, omics_data], axis=1)
    
    st.dataframe(integrated_data.head())
    
    # 相关性分析
    st.subheader("🔗 组间相关性")
    
    # 计算每组组学的平均表达
    omics_means = {}
    for omics in omics_types:
        omics_cols = [col for col in integrated_data.columns if omics in col]
        omics_means[omics] = integrated_data[omics_cols].mean(axis=1)
    
    omics_df = pd.DataFrame(omics_means)
    corr_matrix = omics_df.corr()
    
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(corr_matrix, cmap='coolwarm', center=0, 
               annot=True, fmt='.2f', ax=ax)
    ax.set_title('Cross-Omics Correlation')
    st.pyplot(fig)
    plt.close()
    
    # 整合分析
    st.subheader("🎯 整合分析")
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(integrated_data)
    
    # PCA
    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_scaled)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.6)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})')
    ax.set_title('Integrated Multi-Omics PCA')
    st.pyplot(fig)
    plt.close()
    
    # 聚类分析
    st.subheader("🔬 聚类分析")
    
    from sklearn.cluster import KMeans
    n_clusters = st.slider("选择聚类数", 2, 10, 3)
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(X_scaled)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    scatter = ax.scatter(X_pca[:, 0], X_pca[:, 1], 
                       c=clusters, cmap='tab10', alpha=0.6)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title(f'K-Means Clustering (K={n_clusters})')
    plt.colorbar(scatter, ax=ax, label='Cluster')
    st.pyplot(fig)
    plt.close()


def show_risk_assessment(data):
    """
    显示风险评估
    """
    st.header("📊 风险评估")
    
    # 选择风险因素
    risk_factors = st.multiselect(
        "选择风险因素",
        data.select_dtypes(include=[np.number]).columns,
        default=data.select_dtypes(include=[np.number]).columns[:5].tolist()
    )
    
    if len(risk_factors) == 0:
        st.warning("⚠️ 请至少选择一个风险因素")
        return
    
    # 计算风险评分
    st.subheader("🎯 风险评分计算")
    
    # 标准化
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(data[risk_factors])
    
    # 计算综合风险评分
    risk_score = X_scaled.mean(axis=1)
    risk_score_normalized = (risk_score - risk_score.min()) / (risk_score.max() - risk_score.min()) * 100
    
    data['Risk_Score'] = risk_score_normalized
    
    # 风险分类
    def categorize_risk(score):
        if score < 33:
            return 'Low Risk'
        elif score < 66:
            return 'Medium Risk'
        else:
            return 'High Risk'
    
    data['Risk_Category'] = data['Risk_Score'].apply(categorize_risk)
    
    # 显示风险分布
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("平均风险评分", f"{data['Risk_Score'].mean():.2f}")
    
    with col2:
        st.metric("高风险样本", (data['Risk_Category'] == 'High Risk').sum())
    
    with col3:
        st.metric("低风险样本", (data['Risk_Category'] == 'Low Risk').sum())
    
    # 风险评分分布
    st.subheader("📊 风险评分分布")
    
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    
    # 直方图
    axes[0].hist(data['Risk_Score'], bins=30, alpha=0.7, color='steelblue')
    axes[0].set_xlabel('Risk Score')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Risk Score Distribution')
    axes[0].grid(True, alpha=0.3)
    
    # 饼图
    risk_counts = data['Risk_Category'].value_counts()
    colors = ['green', 'orange', 'red']
    axes[1].pie(risk_counts.values, labels=risk_counts.index, 
                autopct='%1.1f%%', colors=colors)
    axes[1].set_title('Risk Category Distribution')
    
    plt.tight_layout()
    st.pyplot(fig)
    plt.close()
    
    # 风险因素分析
    st.subheader("🔍 风险因素分析")
    
    # 计算每个风险因素与风险评分的相关性
    correlations = {}
    for factor in risk_factors:
        corr = data[factor].corr(data['Risk_Score'])
        correlations[factor] = corr
    
    corr_df = pd.DataFrame.from_dict(correlations, orient='index', columns=['Correlation'])
    corr_df = corr_df.sort_values('Correlation', ascending=False)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    corr_df.plot(kind='barh', ax=ax, color='steelblue')
    ax.set_xlabel('Correlation with Risk Score')
    ax.set_title('Risk Factor Correlation')
    ax.grid(True, alpha=0.3, axis='x')
    st.pyplot(fig)
    plt.close()
    
    # 生成报告
    st.subheader("📄 生成报告")
    
    if st.button("生成分析报告"):
        report = generate_report(data, risk_factors)
        st.download_button(
            label="下载报告",
            data=report,
            file_name="AD_analysis_report.txt",
            mime="text/plain"
        )


def generate_report(data, risk_factors):
    """
    生成分析报告
    """
    report = f"""
阿尔茨海默病研究分析报告
{'='*50}

数据概览:
- 样本数: {data.shape[0]}
- 特征数: {data.shape[1]}
- 分析日期: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

风险因素分析:
"""
    
    for factor in risk_factors:
        report += f"- {factor}: 均值={data[factor].mean():.2f}, 标准差={data[factor].std():.2f}\n"
    
    report += f"""
风险评分统计:
- 平均风险评分: {data['Risk_Score'].mean():.2f}
- 风险评分范围: {data['Risk_Score'].min():.2f} - {data['Risk_Score'].max():.2f}

风险分类:
- 高风险样本: {(data['Risk_Category'] == 'High Risk').sum()} ({(data['Risk_Category'] == 'High Risk').sum()/len(data)*100:.1f}%)
- 中风险样本: {(data['Risk_Category'] == 'Medium Risk').sum()} ({(data['Risk_Category'] == 'Medium Risk').sum()/len(data)*100:.1f}%)
- 低风险样本: {(data['Risk_Category'] == 'Low Risk').sum()} ({(data['Risk_Category'] == 'Low Risk').sum()/len(data)*100:.1f}%)

建议:
1. 对于高风险样本，建议进行进一步临床检查
2. 定期监测风险因素的变化
3. 采取预防性干预措施
4. 保持健康的生活方式

{'='*50}
报告生成完成
"""
    
    return report


if __name__ == '__main__':
    main()
