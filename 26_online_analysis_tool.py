"""
Alzheimer's Disease Research Analysis Tool
Multi-Omics Data Analysis Platform with Network Pharmacology

Features:
1. Transcriptomics Analysis
2. Proteomics Analysis
3. Metabolomics Analysis
4. Single-Cell RNA-seq Analysis
5. Spatial Transcriptomics Analysis
6. Multi-Omics Integration
7. Gene Query Functionality
8. Network Pharmacology Analysis (Traditional Chinese Medicine)

Author: Monash University | Andy's Lab
Date: 2026-02-04
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import cross_val_score
from sklearn.feature_selection import SelectKBest, f_classif
import io
import base64
import os
import warnings
warnings.filterwarnings('ignore')

# Set page configuration
st.set_page_config(
    page_title="ADgene-tool",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Load public databases
@st.cache_data
def load_adni_data():
    """Load ADNI imaging data"""
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        adni_path = os.path.join(current_dir, "../results/ADNI_imaging/ADNI_Imaging_Summary.csv")
        if os.path.exists(adni_path):
            return pd.read_csv(adni_path)
    except:
        pass
    return None

@st.cache_data
def load_nacc_data():
    """Load NACC clinical data"""
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        nacc_path = os.path.join(current_dir, "NACC_filtered_summary.csv")
        if os.path.exists(nacc_path):
            return pd.read_csv(nacc_path)
    except:
        pass
    return None

# Load data
adni_data = load_adni_data()
nacc_data = load_nacc_data()

# Custom CSS for clean white and blue design
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap');
    
    .main {
        background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
    }
    .stApp {
        background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
    }
    
    h1, h2, h3 {
        font-family: 'Roboto', sans-serif;
        color: #2c3e50;
        font-weight: 700;
    }
    
    .css-1d391kg {
        background: linear-gradient(180deg, #ffffff 0%, #f8f9fa 100%);
        border-right: 2px solid #3498db;
    }
    
    .stButton>button {
        background: linear-gradient(135deg, #3498db 0%, #2980b9 100%);
        color: white;
        font-weight: 700;
        font-family: 'Roboto', sans-serif;
        border-radius: 8px;
        padding: 12px 24px;
        border: 2px solid #3498db;
        box-shadow: 0 4px 15px rgba(52, 152, 219, 0.3);
        transition: all 0.3s ease;
        text-transform: uppercase;
        letter-spacing: 1px;
    }
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(52, 152, 219, 0.5);
        background: linear-gradient(135deg, #5dade2 0%, #3498db 100%);
    }
    
    .stDownloadButton>button {
        background: linear-gradient(135deg, #1abc9c 0%, #16a085 100%);
        color: white;
        font-weight: 700;
        font-family: 'Roboto', sans-serif;
        border-radius: 8px;
        padding: 12px 24px;
        border: 2px solid #1abc9c;
        box-shadow: 0 4px 15px rgba(26, 188, 156, 0.3);
        text-transform: uppercase;
        letter-spacing: 1px;
    }
    
    div[data-testid="stMetricValue"] {
        background: linear-gradient(135deg, rgba(52, 152, 219, 0.1) 0%, rgba(41, 128, 185, 0.05) 100%);
        border-radius: 12px;
        padding: 16px;
        box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        border: 1px solid rgba(52, 152, 219, 0.3);
        color: #2c3e50;
        font-family: 'Roboto', sans-serif;
        font-size: 24px;
    }
    
    div[data-testid="stMetricLabel"] {
        color: #3498db;
        font-family: 'Roboto', sans-serif;
        font-weight: 500;
        text-transform: uppercase;
        letter-spacing: 1px;
    }
    
    .stInfo {
        background: linear-gradient(135deg, rgba(52, 152, 219, 0.1) 0%, rgba(41, 128, 185, 0.05) 100%);
        border-radius: 12px;
        padding: 24px;
        color: #2c3e50;
        border: 2px solid #3498db;
        box-shadow: 0 4px 15px rgba(52, 152, 219, 0.2);
    }
    
    .stSuccess {
        background: linear-gradient(135deg, rgba(46, 204, 113, 0.1) 0%, rgba(39, 174, 96, 0.05) 100%);
        border-radius: 12px;
        padding: 24px;
        color: #27ae60;
        border: 2px solid #27ae60;
        box-shadow: 0 4px 15px rgba(39, 174, 96, 0.2);
    }
    
    .stWarning {
        background: linear-gradient(135deg, rgba(241, 196, 15, 0.1) 0%, rgba(243, 156, 18, 0.05) 100%);
        border-radius: 12px;
        padding: 24px;
        color: #f39c12;
        border: 2px solid #f39c12;
        box-shadow: 0 4px 15px rgba(243, 156, 18, 0.2);
    }
    
    .logo-container {
        display: flex;
        justify-content: flex-end;
        align-items: center;
        padding: 20px 30px;
        margin-bottom: 30px;
        background: linear-gradient(135deg, rgba(52, 152, 219, 0.1) 0%, rgba(41, 128, 185, 0.05) 100%);
        border-radius: 12px;
        border: 1px solid rgba(52, 152, 219, 0.3);
        box-shadow: 0 4px 15px rgba(0,0,0,0.1);
    }
    .logo-text {
        font-size: 16px;
        color: #2c3e50;
        font-style: italic;
        font-weight: 400;
        font-family: 'Roboto', sans-serif;
        letter-spacing: 0.5px;
    }
    .logo-divider {
        margin: 0 15px;
        color: #3498db;
        font-weight: 700;
    }
    
    .dataframe {
        background: linear-gradient(135deg, rgba(255, 255, 255, 0.95) 0%, rgba(248, 249, 250, 0.95) 100%);
        border-radius: 12px;
        overflow: hidden;
        border: 1px solid rgba(52, 152, 219, 0.3);
        box-shadow: 0 4px 15px rgba(0,0,0,0.1);
    }
    
    .streamlit-expanderHeader {
        background: linear-gradient(135deg, #3498db 0%, #2980b9 100%);
        border-radius: 12px;
        color: white;
        font-weight: 700;
        font-family: 'Roboto', sans-serif;
        text-transform: uppercase;
        letter-spacing: 1px;
        box-shadow: 0 4px 15px rgba(52, 152, 219, 0.3);
    }
    
    .stSelectbox > div > div {
        background: linear-gradient(135deg, rgba(255, 255, 255, 0.95) 0%, rgba(248, 249, 250, 0.95) 100%);
        border-radius: 10px;
        border: 2px solid rgba(52, 152, 219, 0.5);
        color: #2c3e50;
    }
    
    .stTextInput > div > div > input {
        background: linear-gradient(135deg, rgba(255, 255, 255, 0.95) 0%, rgba(248, 249, 250, 0.95) 100%);
        border-radius: 10px;
        border: 2px solid #3498db;
        color: #2c3e50;
        font-family: 'Roboto', sans-serif;
    }
    
    .stSlider > div > div > div > div {
        background: linear-gradient(135deg, #3498db 0%, #2980b9 100%);
    }
    
    .stSidebar .stSelectbox label, 
    .stSidebar .stFileUploader label {
        color: #2c3e50;
        font-family: 'Roboto', sans-serif;
        font-weight: 500;
    }
    
    hr {
        border-color: rgba(52, 152, 219, 0.5);
        box-shadow: 0 0 10px rgba(52, 152, 219, 0.2);
    }
    
    .stMarkdown {
        color: #2c3e50;
        font-family: 'Roboto', sans-serif;
    }
</style>
""", unsafe_allow_html=True)

# Gene database for query
GENE_DATABASE = {
    'APP': {
        'name': 'Amyloid Beta Precursor Protein',
        'function': 'Precursor protein for amyloid-beta peptides',
        'association': 'Strongly associated with AD pathogenesis',
        'pathway': 'Amyloid processing pathway',
        'expression': 'Highly expressed in brain tissue',
        'risk_level': 'High'
    },
    'PSEN1': {
        'name': 'Presenilin 1',
        'function': 'Component of gamma-secretase complex',
        'association': 'Early-onset AD mutations',
        'pathway': 'Notch signaling pathway',
        'expression': 'Ubiquitous expression',
        'risk_level': 'Very High'
    },
    'PSEN2': {
        'name': 'Presenilin 2',
        'function': 'Component of gamma-secretase complex',
        'association': 'Early-onset AD mutations',
        'pathway': 'Notch signaling pathway',
        'expression': 'Ubiquitous expression',
        'risk_level': 'Very High'
    },
    'APOE': {
        'name': 'Apolipoprotein E',
        'function': 'Lipid transport and metabolism',
        'association': 'Major genetic risk factor for late-onset AD',
        'pathway': 'Lipid metabolism pathway',
        'expression': 'High in liver and brain',
        'risk_level': 'High'
    },
    'MAPT': {
        'name': 'Microtubule-Associated Protein Tau',
        'function': 'Microtubule stabilization',
        'association': 'Tau protein aggregation in AD',
        'pathway': 'Cytoskeleton organization',
        'expression': 'Neuronal expression',
        'risk_level': 'High'
    },
    'TREM2': {
        'name': 'Triggering Receptor Expressed on Myeloid Cells 2',
        'function': 'Immune response in brain',
        'association': 'AD risk variant (R47H)',
        'pathway': 'Microglial activation',
        'expression': 'Microglial cells',
        'risk_level': 'Medium'
    },
    'BIN1': {
        'name': 'Bridging Integrator 1',
        'function': 'Endocytosis and membrane trafficking',
        'association': 'Late-onset AD risk locus',
        'pathway': 'Endocytic pathway',
        'expression': 'Brain and immune cells',
        'risk_level': 'Medium'
    },
    'CLU': {
        'name': 'Clusterin',
        'function': 'Chaperone protein',
        'association': 'Late-onset AD risk locus',
        'pathway': 'Protein folding and clearance',
        'expression': 'Widespread expression',
        'risk_level': 'Medium'
    },
    'CR1': {
        'name': 'Complement Receptor 1',
        'function': 'Immune complement system',
        'association': 'Late-onset AD risk locus',
        'pathway': 'Complement activation',
        'expression': 'Immune cells',
        'risk_level': 'Medium'
    },
    'ABCA7': {
        'name': 'ATP Binding Cassette Subfamily A Member 7',
        'function': 'Lipid transport',
        'association': 'Late-onset AD risk locus',
        'pathway': 'Lipid metabolism',
        'expression': 'Brain and immune cells',
        'risk_level': 'Medium'
    }
}

# Traditional Chinese Medicine (TCM) compounds database
TCM_COMPOUNDS = {
    'Paeoniflorin': {
        'chinese_name': '芍药苷',
        'source': 'Paeonia lactiflora (White Peony)',
        'category': 'Monoterpene glycoside',
        'targets': ['APP', 'BACE1', 'MAPT', 'IL6', 'TNF'],
        'mechanism': 'Anti-inflammatory, neuroprotective',
        'ad_relevance': 'Reduces amyloid-beta accumulation, inhibits tau phosphorylation',
        'evidence_level': 'High'
    },
    'Ferulic acid': {
        'chinese_name': '阿魏酸',
        'source': 'Angelica sinensis (Dong Quai)',
        'category': 'Phenolic acid',
        'targets': ['APP', 'BACE1', 'PSEN1', 'APOE'],
        'mechanism': 'Antioxidant, anti-inflammatory',
        'ad_relevance': 'Reduces oxidative stress, protects neurons',
        'evidence_level': 'High'
    },
    'Hydroxysafflor yellow A': {
        'chinese_name': '羟基红花黄色素A',
        'source': 'Carthamus tinctorius (Safflower)',
        'category': 'Quinochalcone',
        'targets': ['APP', 'BACE1', 'MAPT', 'IL1B', 'TNF'],
        'mechanism': 'Anti-inflammatory, anti-oxidative',
        'ad_relevance': 'Inhibits amyloid-beta aggregation, reduces neuroinflammation',
        'evidence_level': 'Medium'
    },
    'Amygdalin': {
        'chinese_name': '苦杏仁苷',
        'source': 'Prunus armeniaca (Apricot kernel)',
        'category': 'Cyanogenic glycoside',
        'targets': ['APP', 'BACE1', 'ACHE'],
        'mechanism': 'Anti-inflammatory, analgesic',
        'ad_relevance': 'Modulates cholinergic system',
        'evidence_level': 'Medium'
    },
    'Ligustilide': {
        'chinese_name': '藁本内酯',
        'source': 'Angelica sinensis (Dong Quai)',
        'category': 'Phthalide',
        'targets': ['APP', 'BACE1', 'MAPT', 'IL6', 'TNF'],
        'mechanism': 'Anti-inflammatory, neuroprotective',
        'ad_relevance': 'Reduces neuroinflammation, protects blood-brain barrier',
        'evidence_level': 'High'
    },
    'Gallic acid': {
        'chinese_name': '没食子酸',
        'source': 'Paeonia suffruticosa (Tree Peony bark)',
        'category': 'Phenolic acid',
        'targets': ['APP', 'BACE1', 'MAPT', 'ACHE'],
        'mechanism': 'Antioxidant, anti-inflammatory',
        'ad_relevance': 'Reduces oxidative stress, inhibits acetylcholinesterase',
        'evidence_level': 'Medium'
    },
    'Catalpol': {
        'chinese_name': '梓醇',
        'source': 'Rehmannia glutinosa (Rehmannia)',
        'category': 'Iridoid glycoside',
        'targets': ['APP', 'BACE1', 'MAPT', 'BDNF'],
        'mechanism': 'Neuroprotective, anti-inflammatory',
        'ad_relevance': 'Enhances neurogenesis, reduces amyloid-beta toxicity',
        'evidence_level': 'High'
    },
    'Astragaloside IV': {
        'chinese_name': '黄芪甲苷',
        'source': 'Astragalus membranaceus (Astragalus)',
        'category': 'Saponin',
        'targets': ['APP', 'BACE1', 'MAPT', 'IL6', 'TNF'],
        'mechanism': 'Anti-inflammatory, antioxidant',
        'ad_relevance': 'Reduces neuroinflammation, protects mitochondria',
        'evidence_level': 'High'
    },
    'Tanshinone IIA': {
        'chinese_name': '丹参酮IIA',
        'source': 'Salvia miltiorrhiza (Danshen)',
        'category': 'Diterpenoid quinone',
        'targets': ['APP', 'BACE1', 'MAPT', 'IL1B', 'TNF'],
        'mechanism': 'Anti-inflammatory, antioxidant',
        'ad_relevance': 'Inhibits amyloid-beta aggregation, reduces tau phosphorylation',
        'evidence_level': 'High'
    },
    'Baicalin': {
        'chinese_name': '黄芩苷',
        'source': 'Scutellaria baicalensis (Baical Skullcap)',
        'category': 'Flavone glycoside',
        'targets': ['APP', 'BACE1', 'MAPT', 'IL6', 'TNF'],
        'mechanism': 'Anti-inflammatory, antioxidant',
        'ad_relevance': 'Reduces neuroinflammation, protects neurons',
        'evidence_level': 'High'
    }
}

# TCM formula database
TCM_FORMULAS = {
    'Taohong Siwu Decoction': {
        'chinese_name': '桃红四物汤',
        'components': ['Paeonia lactiflora', 'Angelica sinensis', 'Rehmannia glutinosa', 
                      'Ligusticum chuanxiong', 'Prunus persica', 'Carthamus tinctorius'],
        'major_compounds': ['Paeoniflorin', 'Ferulic acid', 'Ligustilide', 'Amygdalin', 
                           'Hydroxysafflor yellow A', 'Gallic acid'],
        'targets': ['APP', 'BACE1', 'MAPT', 'IL6', 'TNF', 'IL1B', 'ACHE', 'BDNF'],
        'mechanism': 'Blood circulation promotion, anti-inflammatory, neuroprotective',
        'ad_relevance': 'Multi-target action on amyloid-beta, tau, inflammation, and cholinergic system',
        'evidence_level': 'High',
        'clinical_trials': 5
    },
    'Buyang Huanwu Decoction': {
        'chinese_name': '补阳还五汤',
        'components': ['Astragalus membranaceus', 'Angelica sinensis', 'Paeonia lactiflora',
                      'Ligusticum chuanxiong', 'Prunus persica', 'Carthamus tinctorius', 'Earthworm'],
        'major_compounds': ['Astragaloside IV', 'Ferulic acid', 'Paeoniflorin', 'Ligustilide',
                           'Amygdalin', 'Hydroxysafflor yellow A'],
        'targets': ['APP', 'BACE1', 'MAPT', 'IL6', 'TNF', 'IL1B', 'BDNF', 'VEGF'],
        'mechanism': 'Qi supplementation, blood circulation promotion, neuroprotection',
        'ad_relevance': 'Enhances cerebral blood flow, reduces neuroinflammation, promotes neurogenesis',
        'evidence_level': 'High',
        'clinical_trials': 8
    },
    'Danggui Shaoyao San': {
        'chinese_name': '当归芍药散',
        'components': ['Angelica sinensis', 'Paeonia lactiflora', 'Atractylodes macrocephala',
                      'Poria cocos', 'Alisma orientale', 'Ligusticum chuanxiong'],
        'major_compounds': ['Ferulic acid', 'Paeoniflorin', 'Ligustilide', 'Gallic acid'],
        'targets': ['APP', 'BACE1', 'MAPT', 'IL6', 'TNF', 'ACHE'],
        'mechanism': 'Blood circulation promotion, dampness elimination, neuroprotection',
        'ad_relevance': 'Improves cognitive function, reduces neuroinflammation',
        'evidence_level': 'Medium',
        'clinical_trials': 3
    }
}

# AD-related pathways
AD_PATHWAYS = {
    'Amyloid processing': {
        'genes': ['APP', 'BACE1', 'PSEN1', 'PSEN2', 'NCSTN', 'APH1A'],
        'description': 'Processing of amyloid precursor protein to amyloid-beta peptides'
    },
    'Tau phosphorylation': {
        'genes': ['MAPT', 'GSK3B', 'CDK5', 'PPP2CA', 'PPP2R1A'],
        'description': 'Phosphorylation and aggregation of tau protein'
    },
    'Neuroinflammation': {
        'genes': ['IL6', 'TNF', 'IL1B', 'TREM2', 'CR1', 'CLU'],
        'description': 'Inflammatory response in brain tissue'
    },
    'Cholinergic signaling': {
        'genes': ['ACHE', 'CHAT', 'SLC5A7', 'CHRNA7', 'CHRNB2'],
        'description': 'Acetylcholine synthesis and degradation'
    },
    'Oxidative stress': {
        'genes': ['SOD1', 'SOD2', 'CAT', 'GPX1', 'NQO1'],
        'description': 'Reactive oxygen species metabolism'
    },
    'Neurotrophic signaling': {
        'genes': ['BDNF', 'NGF', 'NTF3', 'TRKB', 'NGFR'],
        'description': 'Neuronal growth and survival'
    }
}


def load_omics_data(uploaded_file, data_type):
    """
    Load omics data based on data type
    """
    if uploaded_file is None:
        return None
    
    try:
        if data_type == "Transcriptomics":
            data = pd.read_csv(uploaded_file, index_col=0)
            st.success(f"Transcriptomics data loaded: {data.shape[0]} genes x {data.shape[1]} samples")
        elif data_type == "Proteomics":
            data = pd.read_csv(uploaded_file, index_col=0)
            st.success(f"Proteomics data loaded: {data.shape[0]} proteins x {data.shape[1]} samples")
        elif data_type == "Metabolomics":
            data = pd.read_csv(uploaded_file, index_col=0)
            st.success(f"Metabolomics data loaded: {data.shape[0]} metabolites x {data.shape[1]} samples")
        elif data_type == "Single-Cell RNA-seq":
            data = pd.read_csv(uploaded_file, index_col=0)
            st.success(f"Single-cell data loaded: {data.shape[0]} cells x {data.shape[1]} genes")
        elif data_type == "Spatial Transcriptomics":
            data = pd.read_csv(uploaded_file, index_col=0)
            st.success(f"Spatial transcriptomics data loaded: {data.shape[0]} spots x {data.shape[1]} genes")
        else:
            data = pd.read_csv(uploaded_file)
            st.success(f"Data loaded: {data.shape[0]} rows x {data.shape[1]} columns")
        
        return data
    except Exception as e:
        st.error(f"Error loading data: {e}")
        return None


def preprocess_data(data, method="StandardScaler"):
    """
    Preprocess data with normalization
    """
    if method == "StandardScaler":
        scaler = StandardScaler()
    elif method == "MinMaxScaler":
        scaler = MinMaxScaler()
    else:
        return data
    
    data_scaled = pd.DataFrame(
        scaler.fit_transform(data),
        index=data.index,
        columns=data.columns
    )
    
    return data_scaled


def show_transcriptomics_analysis(data):
    """
    Display transcriptomics analysis
    """
    st.header("Transcriptomics Analysis")
    
    st.subheader("Data Overview")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Genes", data.shape[0])
    
    with col2:
        st.metric("Samples", data.shape[1])
    
    with col3:
        st.metric("Missing Values", data.isnull().sum().sum())
    
    st.markdown("---")
    
    st.subheader("Expression Distribution")
    sample_cols = st.multiselect(
        "Select samples for visualization:",
        data.columns.tolist(),
        default=data.columns[:5].tolist()
    )
    
    if len(sample_cols) > 0:
        fig, axes = plt.subplots(1, len(sample_cols), figsize=(5*len(sample_cols), 4))
        if len(sample_cols) == 1:
            axes = [axes]
        
        colors = ['#667eea', '#764ba2', '#f093fb', '#f5576c', '#4facfe']
        
        for i, sample in enumerate(sample_cols):
            axes[i].hist(data[sample].dropna(), bins=50, alpha=0.7, color=colors[i % len(colors)])
            axes[i].set_xlabel('Expression Level')
            axes[i].set_ylabel('Frequency')
            axes[i].set_title(f'{sample} Distribution')
            axes[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        st.pyplot(fig)
        plt.close()
    
    st.subheader("Differential Expression Analysis")
    
    group1 = st.multiselect("Select Group 1 samples:", data.columns.tolist())
    group2 = st.multiselect("Select Group 2 samples:", data.columns.tolist())
    
    if len(group1) > 0 and len(group2) > 0:
        mean1 = data[group1].mean(axis=1)
        mean2 = data[group2].mean(axis=1)
        
        log2fc = np.log2((mean2 + 1e-10) / (mean1 + 1e-10))
        
        de_results = pd.DataFrame({
            'Gene': data.index,
            'Mean_Group1': mean1,
            'Mean_Group2': mean2,
            'Log2FC': log2fc
        })
        
        de_results['Abs_Log2FC'] = np.abs(de_results['Log2FC'])
        de_results = de_results.sort_values('Abs_Log2FC', ascending=False)
        
        st.dataframe(de_results.head(20))
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.scatter(de_results['Mean_Group1'], de_results['Mean_Group2'], alpha=0.5, s=10)
        ax.set_xlabel('Mean Expression Group 1')
        ax.set_ylabel('Mean Expression Group 2')
        ax.set_title('MA Plot')
        ax.grid(True, alpha=0.3)
        st.pyplot(fig)
        plt.close()


def show_proteomics_analysis(data):
    """
    Display proteomics analysis
    """
    st.header("Proteomics Analysis")
    
    st.subheader("Data Overview")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Proteins", data.shape[0])
    
    with col2:
        st.metric("Samples", data.shape[1])
    
    with col3:
        st.metric("Missing Values", data.isnull().sum().sum())
    
    st.markdown("---")
    
    st.subheader("Protein Expression Heatmap")
    
    top_proteins = st.slider("Number of top proteins to display:", 10, 100, 50)
    top_protein_data = data.iloc[:top_proteins]
    
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(top_protein_data.T, cmap='RdYlBu_r', center=0, 
               ax=ax, cbar_kws={'label': 'Expression'})
    ax.set_title(f'Top {top_proteins} Proteins Expression Heatmap')
    ax.set_xlabel('Proteins')
    ax.set_ylabel('Samples')
    st.pyplot(fig)
    plt.close()
    
    st.subheader("Protein Abundance Distribution")
    sample_cols = st.multiselect(
        "Select samples:",
        data.columns.tolist(),
        default=data.columns[:5].tolist()
    )
    
    if len(sample_cols) > 0:
        fig, axes = plt.subplots(1, len(sample_cols), figsize=(5*len(sample_cols), 4))
        if len(sample_cols) == 1:
            axes = [axes]
        
        colors = ['#667eea', '#764ba2', '#f093fb', '#f5576c', '#4facfe']
        
        for i, sample in enumerate(sample_cols):
            axes[i].hist(data[sample].dropna(), bins=50, alpha=0.7, color=colors[i % len(colors)])
            axes[i].set_xlabel('Abundance')
            axes[i].set_ylabel('Frequency')
            axes[i].set_title(f'{sample} Distribution')
            axes[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        st.pyplot(fig)
        plt.close()


def show_metabolomics_analysis(data):
    """
    Display metabolomics analysis
    """
    st.header("Metabolomics Analysis")
    
    st.subheader("Data Overview")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Metabolites", data.shape[0])
    
    with col2:
        st.metric("Samples", data.shape[1])
    
    with col3:
        st.metric("Missing Values", data.isnull().sum().sum())
    
    st.markdown("---")
    
    st.subheader("Metabolite Correlation Network")
    
    n_metabolites = st.slider("Number of metabolites for correlation:", 10, 50, 20)
    metabolite_subset = data.iloc[:n_metabolites]
    
    corr_matrix = metabolite_subset.corr()
    
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(corr_matrix, cmap='RdYlBu_r', center=0, 
               annot=False, fmt='.2f', ax=ax, 
               cbar_kws={'label': 'Correlation'})
    ax.set_title(f'Metabolite Correlation Matrix (Top {n_metabolites})')
    st.pyplot(fig)
    plt.close()
    
    st.subheader("Metabolite Pathway Analysis")
    
    st.info("Upload pathway annotation file to perform pathway enrichment analysis")
    pathway_file = st.file_uploader("Upload pathway annotation (CSV):", type=['csv'])
    
    if pathway_file is not None:
        pathway_data = pd.read_csv(pathway_file)
        st.dataframe(pathway_data.head())


def show_single_cell_analysis(data):
    """
    Display single-cell RNA-seq analysis
    """
    st.header("Single-Cell RNA-seq Analysis")
    
    st.subheader("Data Overview")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Cells", data.shape[0])
    
    with col2:
        st.metric("Genes", data.shape[1])
    
    with col3:
        st.metric("Missing Values", data.isnull().sum().sum())
    
    st.markdown("---")
    
    st.subheader("Cell Quality Control")
    
    data['n_genes'] = (data > 0).sum(axis=1)
    data['total_counts'] = data.sum(axis=1)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    axes[0].hist(data['n_genes'], bins=50, alpha=0.7, color='#667eea')
    axes[0].set_xlabel('Number of Genes per Cell')
    axes[0].set_ylabel('Number of Cells')
    axes[0].set_title('Genes per Cell Distribution')
    axes[0].grid(True, alpha=0.3)
    
    axes[1].hist(data['total_counts'], bins=50, alpha=0.7, color='#764ba2')
    axes[1].set_xlabel('Total Counts per Cell')
    axes[1].set_ylabel('Number of Cells')
    axes[1].set_title('Total Counts per Cell Distribution')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    st.pyplot(fig)
    plt.close()
    
    st.subheader("Dimensionality Reduction (PCA)")
    
    n_pcs = st.slider("Number of Principal Components:", 2, 10, 2)
    
    gene_data = data.drop(['n_genes', 'total_counts'], axis=1)
    scaler = StandardScaler()
    gene_data_scaled = scaler.fit_transform(gene_data)
    
    pca = PCA(n_components=n_pcs, random_state=42)
    pca_result = pca.fit_transform(gene_data_scaled)
    
    if n_pcs == 2:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.scatter(pca_result[:, 0], pca_result[:, 1], 
                   alpha=0.6, s=20, c='#667eea')
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
        ax.set_title('PCA Visualization')
        ax.grid(True, alpha=0.3)
        st.pyplot(fig)
        plt.close()
    
    st.subheader("Cell Clustering")
    
    n_clusters = st.slider("Number of clusters:", 2, 10, 5)
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(gene_data_scaled)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    colors = ['#667eea', '#764ba2', '#f093fb', '#f5576c', '#4facfe']
    
    for i in range(n_clusters):
        mask = clusters == i
        ax.scatter(pca_result[mask, 0], pca_result[mask, 1], 
                  c=colors[i % len(colors)], alpha=0.6, s=20, 
                  label=f'Cluster {i+1}')
    
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title(f'K-Means Clustering (K={n_clusters})')
    ax.legend()
    ax.grid(True, alpha=0.3)
    st.pyplot(fig)
    plt.close()


def show_spatial_transcriptomics_analysis(data):
    """
    Display spatial transcriptomics analysis
    """
    st.header("Spatial Transcriptomics Analysis")
    
    st.subheader("Data Overview")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Spots", data.shape[0])
    
    with col2:
        st.metric("Genes", data.shape[1])
    
    with col3:
        st.metric("Missing Values", data.isnull().sum().sum())
    
    st.markdown("---")
    
    st.subheader("Spatial Gene Expression")
    
    gene_select = st.selectbox("Select gene for visualization:", data.columns.tolist())
    
    if gene_select:
        st.info("Upload spatial coordinates file to visualize spatial expression")
        coord_file = st.file_uploader("Upload spatial coordinates (CSV):", type=['csv'])
        
        if coord_file is not None:
            coords = pd.read_csv(coord_file)
            
            if len(coords) == len(data):
                fig, ax = plt.subplots(figsize=(10, 8))
                scatter = ax.scatter(coords.iloc[:, 0], coords.iloc[:, 1], 
                                   c=data[gene_select], cmap='RdYlBu_r', 
                                   s=50, alpha=0.7)
                ax.set_xlabel('X Coordinate')
                ax.set_ylabel('Y Coordinate')
                ax.set_title(f'Spatial Expression of {gene_select}')
                plt.colorbar(scatter, ax=ax, label='Expression')
                st.pyplot(fig)
                plt.close()
            else:
                st.warning("Coordinate file length does not match data length")
    
    st.subheader("Spatial PCA")
    
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    
    pca = PCA(n_components=2, random_state=42)
    pca_result = pca.fit_transform(data_scaled)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(pca_result[:, 0], pca_result[:, 1], 
                       c=range(len(pca_result)), cmap='RdYlBu_r', 
                       s=50, alpha=0.7)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
    ax.set_title('Spatial PCA')
    plt.colorbar(scatter, ax=ax, label='Spot Index')
    st.pyplot(fig)
    plt.close()


def show_multi_omics_integration(data_dict):
    """
    Display multi-omics integration
    """
    st.header("Multi-Omics Integration")
    
    st.subheader("Data Integration Overview")
    
    for data_type, data in data_dict.items():
        if data is not None:
            st.info(f"{data_type}: {data.shape[0]} features x {data.shape[1]} samples")
    
    st.markdown("---")
    
    st.subheader("Cross-Omics Correlation")
    
    omics_types = [k for k, v in data_dict.items() if v is not None]
    
    if len(omics_types) >= 2:
        corr_matrix = pd.DataFrame(index=omics_types, columns=omics_types)
        
        for i, omics1 in enumerate(omics_types):
            for j, omics2 in enumerate(omics_types):
                if i <= j:
                    data1 = data_dict[omics1]
                    data2 = data_dict[omics2]
                    
                    min_samples = min(data1.shape[1], data2.shape[1])
                    data1_subset = data1.iloc[:, :min_samples]
                    data2_subset = data2.iloc[:, :min_samples]
                    
                    mean1 = data1_subset.mean(axis=0)
                    mean2 = data2_subset.mean(axis=0)
                    
                    corr = np.corrcoef(mean1, mean2)[0, 1]
                    corr_matrix.loc[omics1, omics2] = corr
                    corr_matrix.loc[omics2, omics1] = corr
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(corr_matrix.astype(float), cmap='Blues', center=0, 
                   annot=True, fmt='.2f', ax=ax, 
                   cbar_kws={'label': 'Correlation'},
                   linewidths=0.5, linecolor='#ffffff')
        ax.set_title('Cross-Omics Correlation Matrix', fontsize=14, 
                    fontweight='bold', color='#2c3e50')
        ax.tick_params(colors='#2c3e50')
        st.pyplot(fig)
        plt.close()
    else:
        st.warning("Please upload at least two omics datasets for integration")


def show_gene_query(t):
    """
    Display gene query interface
    """
    st.header(t['gene_query'])
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        gene_input = st.text_input(
            t['enter_gene'],
            placeholder=t['gene_placeholder'],
            key="gene_query_input"
        )
        
        if st.button(t['search_gene'], key="search_gene_btn"):
            if gene_input:
                gene_info = query_gene(gene_input)
                
                if gene_info:
                    st.success(f"{t['gene_found']} {gene_input.upper()}")
                    
                    with st.expander(t['gene_details'], expanded=True):
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.markdown(f"**{t['full_name']}** {gene_info['name']}")
                            st.markdown(f"**{t['function']}** {gene_info['function']}")
                        
                        with col2:
                            st.markdown(f"**{t['pathway']}** {gene_info['pathway']}")
                            st.markdown(f"**{t['expression']}** {gene_info['expression']}")
                    
                    st.markdown("---")
                    
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.metric(t['risk_level'], gene_info['risk_level'])
                    
                    with col2:
                        st.metric("AD Association", "Confirmed")
                    
                    with col3:
                        st.metric("Research Priority", "High")
                    
                    st.markdown("---")
                    st.info(f"**{t['association']}:** {gene_info['association']}")
                    
                    st.markdown("---")
                    st.markdown(f"### {t['gene_visualization']}")
                    
                    fig, ax = plt.subplots(figsize=(10, 8))
                    
                    categories = ['Function', 'Pathway', 'Expression', 'Risk']
                    values = [4, 4, 4, 4]
                    
                    colors = ['#3498db', '#2980b9', '#5dade2', '#1abc9c']
                    
                    bars = ax.bar(categories, values, color=colors, edgecolor='#2c3e50', linewidth=2)
                    
                    ax.set_ylabel('Score', color='#2c3e50', fontsize=12)
                    ax.set_title(f'Gene Profile: {gene_input.upper()}', 
                               fontsize=14, fontweight='bold', color='#2c3e50')
                    
                    for bar in bars:
                        height = bar.get_height()
                        ax.text(bar.get_x() + bar.get_width()/2., height,
                               f'{height}', ha='center', va='bottom',
                               color='#2c3e50', fontsize=12, fontweight='bold')
                    
                    ax.tick_params(colors='#2c3e50')
                    ax.spines['bottom'].set_color('#3498db')
                    ax.spines['left'].set_color('#3498db')
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    
                    st.pyplot(fig)
                    plt.close()
                    
                    st.markdown("---")
                    st.markdown(f"### {t['download_image']}")
                    
                    download_format = st.selectbox(
                        t['download_format'],
                        ["PNG", "TIFF"],
                        key="download_format"
                    )
                    
                    if download_format == "TIFF":
                        dpi_option = st.selectbox(
                            "Select DPI:",
                            [300, 600],
                            key="download_dpi"
                        )
                    else:
                        dpi_option = None
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        if st.download_button(
                            label=f"{t['downloaded']} {download_format}",
                            data=fig_to_bytes(fig, download_format, dpi_option),
                            file_name=f"gene_{gene_input.lower()}.{download_format.lower()}",
                            mime=f"image/{download_format.lower()}",
                            key=f"download_{download_format.lower()}"
                        ):
                            st.success(f"{t['downloaded']} {download_format}")
                    
                    with col2:
                        if st.download_button(
                            label="Download High Resolution",
                            data=fig_to_bytes(fig, "TIFF", 600),
                            file_name=f"gene_{gene_input.lower()}_600dpi.tiff",
                            mime="image/tiff",
                            key="download_high_res"
                        ):
                            st.success(t['downloaded_tiff'])
                    
                else:
                    st.warning(f"{t['gene_not_found']}: {gene_input}")
                    st.info(t['try_genes'])
    
    with col2:
        st.markdown(f"### {t['gene_database']}")
        
        if t['gene_query'] == "基因查询":
            st.markdown("""
            <div style='background: linear-gradient(135deg, rgba(52, 152, 219, 0.1) 0%, rgba(41, 128, 185, 0.05) 100%);
                        padding: 20px; border-radius: 12px; border: 1px solid rgba(52, 152, 219, 0.3);
                        margin-bottom: 20px;'>
                <p style='color: #2c3e50; font-size: 14px;'>
                    <strong>可用基因</strong> (10个基因):
                </p>
                <ul style='color: #2c3e50; font-size: 13px;'>
                    <li>APP - 淀粉样前体蛋白</li>
                    <li>PSEN1 - 早老素1</li>
                    <li>PSEN2 - 早老素2</li>
                    <li>APOE - 载脂蛋白E</li>
                    <li>MAPT - 微管相关蛋白Tau</li>
                    <li>TREM2 - 髓系细胞触发受体2</li>
                    <li>BIN1 - 桥接整合因子1</li>
                    <li>CLU - 聚集素</li>
                    <li>CR1 - 补体受体1</li>
                    <li>ABCA7 - ATP结合盒转运蛋白A7</li>
                </ul>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.markdown("""
            <div style='background: linear-gradient(135deg, rgba(52, 152, 219, 0.1) 0%, rgba(41, 128, 185, 0.05) 100%);
                        padding: 20px; border-radius: 12px; border: 1px solid rgba(52, 152, 219, 0.3);
                        margin-bottom: 20px;'>
                <p style='color: #2c3e50; font-size: 14px;'>
                    <strong>Available Genes</strong> (10 genes):
                </p>
                <ul style='color: #2c3e50; font-size: 13px;'>
                    <li>APP - Amyloid Beta Precursor Protein</li>
                    <li>PSEN1 - Presenilin 1</li>
                    <li>PSEN2 - Presenilin 2</li>
                    <li>APOE - Apolipoprotein E</li>
                    <li>MAPT - Microtubule-Associated Protein Tau</li>
                    <li>TREM2 - Triggering Receptor Expressed on Myeloid Cells 2</li>
                    <li>BIN1 - Bridging Integrator 1</li>
                    <li>CLU - Clusterin</li>
                    <li>CR1 - Complement Receptor 1</li>
                    <li>ABCA7 - ATP Binding Cassette Subfamily A Member 7</li>
                </ul>
            </div>
            """, unsafe_allow_html=True)


def fig_to_bytes(fig, format_type, dpi=None):
    """
    Convert matplotlib figure to bytes for download
    """
    buf = io.BytesIO()
    
    if format_type == "PNG":
        fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    elif format_type == "TIFF":
        dpi_value = dpi if dpi else 300
        fig.savefig(buf, format='tiff', dpi=dpi_value, bbox_inches='tight')
    
    buf.seek(0)
    return buf.getvalue()


def query_gene(gene_name):
    """
    Query gene information from database
    """
    gene_name_upper = gene_name.upper()
    
    if gene_name_upper in GENE_DATABASE:
        return GENE_DATABASE[gene_name_upper]
    else:
        return None


def show_network_pharmacology(t):
    """
    Display network pharmacology analysis interface
    """
    st.header(t['network_pharmacology'])
    
    if t['network_pharmacology'] == "网络药理学":
        st.markdown(f"""
        <div style='background: linear-gradient(135deg, rgba(52, 152, 219, 0.1) 0%, rgba(41, 128, 185, 0.05) 100%);
                    padding: 20px; border-radius: 12px; border: 1px solid rgba(52, 152, 219, 0.3);
                    margin-bottom: 20px;'>
            <p style='color: #2c3e50; font-size: 16px;'>
                <strong>网络药理学分析</strong> 整合中药成分与其分子靶点和AD相关通路，揭示多靶点治疗机制。
            </p>
        </div>
        """, unsafe_allow_html=True)
        
        analysis_modes = ["中药成分分析", "中药方剂分析", "成分-靶点网络", "通路富集分析"]
    else:
        st.markdown(f"""
        <div style='background: linear-gradient(135deg, rgba(52, 152, 219, 0.1) 0%, rgba(41, 128, 185, 0.05) 100%;
                    padding: 20px; border-radius: 12px; border: 1px solid rgba(52, 152, 219, 0.3);
                    margin-bottom: 20px;'>
            <p style='color: #2c3e50; font-size: 16px;'>
                <strong>Network Pharmacology Analysis</strong> integrates Traditional Chinese Medicine (TCM) 
                compounds with their molecular targets and AD-related pathways to uncover multi-target 
                therapeutic mechanisms.
            </p>
        </div>
        """, unsafe_allow_html=True)
        
        analysis_modes = ["TCM Compound Analysis", "TCM Formula Analysis", "Compound-Target Network", "Pathway Enrichment"]
    
    analysis_mode = st.radio(
        "Select Analysis Mode:",
        analysis_modes,
        horizontal=True
    )
    
    if analysis_mode == analysis_modes[0]:
        show_tcm_compound_analysis(t)
    elif analysis_mode == analysis_modes[1]:
        show_tcm_formula_analysis(t)
    elif analysis_mode == analysis_modes[2]:
        show_compound_target_network(t)
    elif analysis_mode == analysis_modes[3]:
        show_pathway_enrichment(t)


def show_tcm_compound_analysis(t):
    """
    Display TCM compound analysis
    """
    st.subheader(t['tcm_compound_analysis'])
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        compound_input = st.text_input(
            t['enter_compound'],
            placeholder=t['compound_placeholder'],
            key="compound_input"
        )
        
        if st.button(t['search_compound'], key="analyze_compound_btn"):
            if compound_input:
                compound_info = query_compound(compound_input)
                
                if compound_info:
                    st.success(f"{t['compound_found']} {compound_input}")
                    
                    with st.expander(t['compound_details'], expanded=True):
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.markdown(f"**{t['chinese_name']}** {compound_info['chinese_name']}")
                            st.markdown(f"**{t['source']}** {compound_info['source']}")
                            st.markdown(f"**{t['category']}** {compound_info['category']}")
                        
                        with col2:
                            st.markdown(f"**{t['mechanism']}** {compound_info['mechanism']}")
                            st.markdown(f"**{t['evidence_level']}** {compound_info['evidence_level']}")
                    
                    st.markdown("---")
                    
                    st.markdown(f"### {t['target_genes']}")
                    targets = compound_info['targets']
                    
                    for target in targets:
                        with st.expander(f"Target: {target}", expanded=False):
                            if target in GENE_DATABASE:
                                gene_info = GENE_DATABASE[target]
                                st.markdown(f"**Name:** {gene_info['name']}")
                                st.markdown(f"**Function:** {gene_info['function']}")
                                st.markdown(f"**AD Association:** {gene_info['association']}")
                    
                    st.markdown("---")
                    st.info(f"**{t['ad_relevance']}:** {compound_info['ad_relevance']}")
                    
                else:
                    st.warning(f"{t['compound_not_found']}: {compound_input}")
                    st.info("Available compounds: Paeoniflorin, Ferulic acid, Hydroxysafflor yellow A, Amygdalin, Ligustilide, Gallic acid, Catalpol, Astragaloside IV, Tanshinone IIA, Baicalin")
    
    with col2:
        st.markdown(f"### {t['available_compounds']}")
        
        compound_list = list(TCM_COMPOUNDS.keys())
        
        for i, compound in enumerate(compound_list):
            if st.button(f"{compound}", key=f"compound_btn_{i}"):
                st.session_state['compound_input'] = compound
                st.rerun()


def show_tcm_formula_analysis(t):
    """
    Display TCM formula analysis
    """
    st.subheader(t['tcm_formula_analysis'])
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        formula_input = st.text_input(
            t['enter_formula'],
            placeholder=t['formula_placeholder'],
            key="formula_input"
        )
        
        if st.button(t['search_formula'], key="analyze_formula_btn"):
            if formula_input:
                formula_info = query_formula(formula_input)
                
                if formula_info:
                    st.success(f"{t['formula_found']} {formula_input}")
                    
                    with st.expander(t['formula_details'], expanded=True):
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.markdown(f"**{t['chinese_name']}** {formula_info['chinese_name']}")
                            st.markdown(f"**{t['mechanism']}** {formula_info['mechanism']}")
                            st.markdown(f"**{t['evidence_level']}** {formula_info['evidence_level']}")
                        
                        with col2:
                            st.markdown(f"**{t['clinical_trials']}** {formula_info['clinical_trials']}")
                            st.markdown(f"**{t['major_compounds']}** {len(formula_info['major_compounds'])}")
                            st.markdown(f"**{t['target_genes']}** {len(formula_info['targets'])}")
                    
                    st.markdown("---")
                    
                    st.markdown(f"### {t['herbal_components']}")
                    for component in formula_info['components']:
                        st.markdown(f"- {component}")
                    
                    st.markdown("---")
                    
                    st.markdown(f"### {t['major_compounds']}")
                    for compound in formula_info['major_compounds']:
                        st.markdown(f"- {compound}")
                    
                    st.markdown("---")
                    
                    st.markdown(f"### {t['target_genes']}")
                    for target in formula_info['targets']:
                        st.markdown(f"- {target}")
                    
                    st.markdown("---")
                    st.info(f"**{t['ad_relevance']}:** {formula_info['ad_relevance']}")
                    
                else:
                    st.warning(f"{t['formula_not_found']}: {formula_input}")
                    st.info("Available formulas: Taohong Siwu Decoction, Buyang Huanwu Decoction, Danggui Shaoyao San")
    
    with col2:
        st.markdown("### Available Formulas")
        
        formula_list = list(TCM_FORMULAS.keys())
        
        for i, formula in enumerate(formula_list):
            if st.button(f"{formula}", key=f"formula_btn_{i}"):
                st.session_state['formula_input'] = formula
                st.rerun()


def show_compound_target_network(t):
    """
    Display compound-target network visualization
    """
    st.subheader("Compound-Target Network")
    
    st.info("Select compounds to visualize their interaction network with AD-related targets")
    
    selected_compounds = st.multiselect(
        t['select_compounds'],
        list(TCM_COMPOUNDS.keys()),
        default=['Paeoniflorin', 'Ferulic acid']
    )
    
    if len(selected_compounds) > 0:
        st.markdown(f"### {t['network_visualization']}")
        
        fig, ax = plt.subplots(figsize=(14, 10))
        
        all_targets = set()
        compound_targets = {}
        
        for compound in selected_compounds:
            targets = TCM_COMPOUNDS[compound]['targets']
            compound_targets[compound] = targets
            all_targets.update(targets)
        
        all_targets = sorted(list(all_targets))
        
        pos = {}
        n_compounds = len(selected_compounds)
        n_targets = len(all_targets)
        
        for i, compound in enumerate(selected_compounds):
            angle = 2 * np.pi * i / n_compounds
            pos[compound] = (0.6 * np.cos(angle), 0.6 * np.sin(angle))
        
        for i, target in enumerate(all_targets):
            angle = 2 * np.pi * i / n_targets
            pos[target] = (0.3 * np.cos(angle), 0.3 * np.sin(angle))
        
        colors_compounds = ['#3498db', '#2980b9', '#5dade2', '#1abc9c', '#16a085']
        
        for i, compound in enumerate(selected_compounds):
            x, y = pos[compound]
            ax.scatter(x, y, s=2000, c=colors_compounds[i % len(colors_compounds)], 
                      edgecolors='black', linewidths=2, zorder=3)
            ax.text(x, y, compound, ha='center', va='center', fontsize=10, 
                   fontweight='bold', color='#1a1a2e', zorder=4)
            
            for target in compound_targets[compound]:
                tx, ty = pos[target]
                ax.plot([x, tx], [y, ty], color='#3498db', alpha=0.6, 
                       linewidth=1.5, zorder=1)
        
        for i, target in enumerate(all_targets):
            x, y = pos[target]
            ax.scatter(x, y, s=800, c='#1abc9c', edgecolors='black', 
                      linewidths=1.5, zorder=2)
            ax.text(x, y, target, ha='center', va='center', fontsize=8, 
                   fontweight='bold', color='#2c3e50', zorder=3)
        
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title('Compound-Target Network', fontsize=16, fontweight='bold', 
                    color='#2c3e50', pad=20)
        
        legend_elements = [
            plt.scatter([], [], s=200, c='#3498db', edgecolors='black', 
                       linewidths=2, label='TCM Compound'),
            plt.scatter([], [], s=100, c='#1abc9c', edgecolors='black', 
                       linewidths=1.5, label='Target Gene')
        ]
        ax.legend(handles=legend_elements, loc='upper right', 
                 facecolor='#ffffff', edgecolor='#3498db', 
                 labelcolor='#2c3e50', fontsize=10)
        
        st.pyplot(fig)
        plt.close()
        
        st.markdown("---")
        st.markdown(f"### {t['network_statistics']}")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric(t['compounds'], len(selected_compounds))
        
        with col2:
            st.metric(t['unique_targets'], len(all_targets))
        
        with col3:
            st.metric(t['total_edges'], sum(len(targets) for targets in compound_targets.values()))
        
        with col4:
            avg_targets = sum(len(targets) for targets in compound_targets.values()) / len(selected_compounds)
            st.metric(t['avg_targets_compound'], f"{avg_targets:.1f}")


def show_pathway_enrichment(t):
    """
    Display pathway enrichment analysis
    """
    st.subheader(t['pathway_enrichment'])
    
    st.info(t['select_compounds_enrichment'])
    
    analysis_type = st.radio(
        "Analyze:",
        ["TCM Compounds", "TCM Formulas"],
        horizontal=True
    )
    
    if analysis_type == "TCM Compounds":
        selected_compounds = st.multiselect(
            "Select TCM Compounds:",
            list(TCM_COMPOUNDS.keys()),
            default=['Paeoniflorin', 'Ferulic acid', 'Ligustilide']
        )
        
        if len(selected_compounds) > 0:
            all_targets = set()
            for compound in selected_compounds:
                all_targets.update(TCM_COMPOUNDS[compound]['targets'])
            
            pathway_scores = calculate_pathway_enrichment(all_targets)
            
            display_pathway_results(pathway_scores, selected_compounds, all_targets)
    
    else:
        selected_formulas = st.multiselect(
            "Select TCM Formulas:",
            list(TCM_FORMULAS.keys()),
            default=['Taohong Siwu Decoction']
        )
        
        if len(selected_formulas) > 0:
            all_targets = set()
            for formula in selected_formulas:
                all_targets.update(TCM_FORMULAS[formula]['targets'])
            
            pathway_scores = calculate_pathway_enrichment(all_targets)
            
            display_pathway_results(pathway_scores, selected_formulas, all_targets)


def calculate_pathway_enrichment(targets):
    """
    Calculate pathway enrichment scores
    """
    pathway_scores = {}
    
    for pathway_name, pathway_info in AD_PATHWAYS.items():
        pathway_genes = set(pathway_info['genes'])
        overlap = len(targets.intersection(pathway_genes))
        
        if overlap > 0:
            enrichment_score = overlap / len(pathway_genes)
            pathway_scores[pathway_name] = {
                'score': enrichment_score,
                'overlap': overlap,
                'pathway_genes': pathway_genes,
                'description': pathway_info['description']
            }
    
    return pathway_scores


def display_pathway_results(pathway_scores, selection, targets):
    """
    Display pathway enrichment results
    """
    st.markdown("### Enrichment Results")
    
    if len(pathway_scores) == 0:
        st.warning("No pathway enrichment found for selected targets")
        return
    
    sorted_pathways = sorted(pathway_scores.items(), 
                            key=lambda x: x[1]['score'], reverse=True)
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        pathways = [p[0] for p in sorted_pathways]
        scores = [p[1]['score'] for p in sorted_pathways]
        overlaps = [p[1]['overlap'] for p in sorted_pathways]
        
        colors = ['#3498db', '#2980b9', '#5dade2', '#1abc9c', '#16a085']
        
        bars = ax.barh(range(len(pathways)), scores, color=colors[:len(pathways)])
        
        ax.set_yticks(range(len(pathways)))
        ax.set_yticklabels(pathways)
        ax.set_xlabel('Enrichment Score', color='#2c3e50', fontsize=12)
        ax.set_title('Pathway Enrichment Analysis', fontsize=14, fontweight='bold', 
                    color='#2c3e50')
        
        for i, (bar, overlap) in enumerate(zip(bars, overlaps)):
            width = bar.get_width()
            ax.text(width + 0.01, bar.get_y() + bar.get_height()/2, 
                   f'{overlap} genes', ha='left', va='center', 
                   color='#2c3e50', fontsize=10)
        
        ax.tick_params(axis='x', colors='#2c3e50')
        ax.tick_params(axis='y', colors='#2c3e50')
        ax.spines['bottom'].set_color('#3498db')
        ax.spines['left'].set_color('#3498db')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        st.pyplot(fig)
        plt.close()
    
    with col2:
        st.markdown("### Target Summary")
        st.metric("Total Targets", len(targets))
        st.metric("Pathways Hit", len(pathway_scores))
        
        st.markdown("### Selected Targets")
        for target in sorted(list(targets)):
            st.markdown(f"- {target}")
    
    st.markdown("---")
    
    st.markdown("### Detailed Pathway Information")
    
    for pathway_name, pathway_data in sorted_pathways:
        with st.expander(f"{pathway_name} (Score: {pathway_data['score']:.2f})", expanded=False):
            st.markdown(f"**Description:** {pathway_data['description']}")
            st.markdown(f"**Overlap:** {pathway_data['overlap']} genes")
            st.markdown(f"**Pathway Genes:** {', '.join(pathway_data['pathway_genes'])}")
            
            overlapping_genes = set(targets).intersection(pathway_data['pathway_genes'])
            st.markdown(f"**Targeted Genes:** {', '.join(sorted(overlapping_genes))}")


def query_compound(compound_name):
    """
    Query compound information from database
    """
    compound_name_capitalized = ' '.join(word.capitalize() for word in compound_name.split())
    
    if compound_name_capitalized in TCM_COMPOUNDS:
        return TCM_COMPOUNDS[compound_name_capitalized]
    else:
        return None


def query_formula(formula_name):
    """
    Query formula information from database
    """
    formula_name_capitalized = ' '.join(word.capitalize() for word in formula_name.split())
    
    if formula_name_capitalized in TCM_FORMULAS:
        return TCM_FORMULAS[formula_name_capitalized]
    else:
        return None


def show_public_databases(t):
    """
    Display public databases (ADNI and NACC)
    """
    st.header(t['public_databases'])
    
    if t['public_databases'] == "公共数据库":
        st.markdown(f"""
        <div style='background: linear-gradient(135deg, rgba(52, 152, 219, 0.1) 0%, rgba(41, 128, 185, 0.05) 100%);
                    padding: 20px; border-radius: 12px; border: 1px solid rgba(52, 152, 219, 0.3);
                    margin-bottom: 20px;'>
            <p style='color: #2c3e50; font-size: 16px;'>
                <strong>公共数据库</strong> 为阿尔茨海默病研究提供全面的临床和影像数据，包括 ADNI（阿尔茨海默病神经影像学计划）和 NACC（国家阿尔茨海默病协调中心）。
            </p>
        </div>
        """, unsafe_allow_html=True)
    else:
        st.markdown(f"""
        <div style='background: linear-gradient(135deg, rgba(52, 152, 219, 0.1) 0%, rgba(41, 128, 185, 0.05) 100%);
                    padding: 20px; border-radius: 12px; border: 1px solid rgba(52, 152, 219, 0.3);
                    margin-bottom: 20px;'>
            <p style='color: #2c3e50; font-size: 16px;'>
                <strong>Public Databases</strong> provide comprehensive clinical and imaging data 
                for Alzheimer's disease research, including ADNI (Alzheimer's Disease Neuroimaging Initiative) 
                and NACC (National Alzheimer's Coordinating Center).
            </p>
        </div>
        """, unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs([t['adni_database'], t['nacc_database']])
    
    with tab1:
        st.subheader("ADNI - Alzheimer's Disease Neuroimaging Initiative")
        
        if adni_data is not None:
            st.success(t['adni_data_loaded'])
            
            st.markdown(f"### {t['database_overview']}")
            
            for idx, row in adni_data.iterrows():
                metric = row['Metric']
                value = row['Value']
                
                if 'Total' in metric or 'Available' in metric or 'Identified' in metric:
                    st.info(f"**{metric}**: {value}")
                else:
                    st.markdown(f"**{metric}**: {value}")
            
            st.markdown("---")
            st.markdown(f"### {t['data_description']}")
            
            if t['public_databases'] == "公共数据库":
                st.markdown("""
                ADNI数据库提供:
                - **MRI影像**: 结构性脑部影像数据
                - **PET扫描**: 淀粉样蛋白和tau PET影像
                - **临床评估**: 认知和功能测量
                - **纵向数据**: 多时间点测量
                - **生物标志物**: 脑脊液和血液生物标志物
                """)
            else:
                st.markdown("""
                The ADNI database provides:
                - **MRI Imaging**: Structural brain imaging data
                - **PET Scans**: Amyloid and tau PET imaging
                - **Clinical Assessments**: Cognitive and functional measures
                - **Longitudinal Data**: Multi-timepoint measurements
                - **Biomarkers**: CSF and blood biomarkers
                """)
            
        else:
            st.warning(t['adni_data_not_available'])
            st.info(t['data_file_not_found'])
    
    with tab2:
        st.subheader("NACC - National Alzheimer's Coordinating Center")
        
        if nacc_data is not None:
            st.success(f"{t['nacc_data_loaded']} ({len(nacc_data)} records)")
            
            st.markdown(f"### {t['database_overview']}")
            
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric(t['total_subjects'], len(nacc_data))
            
            with col2:
                age_range = f"{nacc_data['NACCAGE'].min()}-{nacc_data['NACCAGE'].max()}"
                st.metric(t['age_range'], age_range)
            
            with col3:
                st.metric(t['females'], (nacc_data['SEX'] == 2).sum())
            
            with col4:
                st.metric(t['males'], (nacc_data['SEX'] == 1).sum())
            
            st.markdown("---")
            
            st.markdown(f"### {t['diagnosis_distribution']}")
            
            diagnosis_counts = nacc_data['Diagnosis_Status'].value_counts()
            
            fig, ax = plt.subplots(figsize=(10, 6))
            colors = ['#3498db', '#2980b9', '#5dade2', '#1abc9c']
            bars = ax.bar(diagnosis_counts.index, diagnosis_counts.values, 
                       color=colors[:len(diagnosis_counts)])
            
            ax.set_xlabel(t['diagnosis'], color='#2c3e50', fontsize=12)
            ax.set_ylabel(t['count'], color='#2c3e50', fontsize=12)
            ax.set_title(f"{t['diagnosis_distribution']} in NACC Database", 
                       fontsize=14, fontweight='bold', color='#2c3e50')
            
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height)}', ha='center', va='bottom',
                       color='#2c3e50', fontsize=10)
            
            ax.tick_params(colors='#2c3e50')
            ax.spines['bottom'].set_color('#3498db')
            ax.spines['left'].set_color('#3498db')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            st.pyplot(fig)
            plt.close()
            
            st.markdown("---")
            
            st.markdown(f"### {t['data_preview']}")
            st.dataframe(nacc_data.head(20), use_container_width=True)
            
            st.markdown("---")
            st.markdown(f"### {t['data_description']}")
            
            if t['public_databases'] == "公共数据库":
                st.markdown("""
                NACC数据库提供:
                - **人口统计学**: 年龄、性别、教育
                - **认知测试**: MMSE、CDR评分
                - **诊断**: CN、MCI、AD分类
                - **遗传学**: APOE基因型信息
                - **临床病史**: 医学和家族病史
                """)
            else:
                st.markdown("""
                The NACC database provides:
                - **Demographics**: Age, sex, education
                - **Cognitive Tests**: MMSE, CDR scores
                - **Diagnosis**: CN, MCI, AD classifications
                - **Genetics**: APOE genotype information
                - **Clinical History**: Medical and family history
                """)
            
        else:
            st.warning(t['nacc_data_not_available'])
            st.info(t['data_file_not_found'])


def main():
    """
    Main function
    """
    # Language selection
    lang = st.sidebar.selectbox(
        "Language / 语言",
        ["中文", "English"],
        index=1
    )
    
    # Language dictionary
    TEXTS = {
        "中文": {
            "title": "阿尔茨海默病多组学分析工具",
            "sidebar_title": "分析选项",
            "select_analysis": "选择分析类型:",
            "analysis_types": ["基因查询", "网络药理学", "公共数据库", "转录组分析", "蛋白质组分析", "代谢组分析", 
                              "单细胞RNA测序", "空间转录组", "多组学整合"],
            "data_upload": "数据上传",
            "upload_transcriptomics": "上传转录组数据 (CSV):",
            "upload_proteomics": "上传蛋白质组数据 (CSV):",
            "upload_metabolomics": "上传代谢组数据 (CSV):",
            "upload_single_cell": "上传单细胞数据 (CSV):",
            "upload_spatial": "上传空间转录组数据 (CSV):",
            "public_databases": "公共数据库",
            "adni_database": "ADNI数据库",
            "nacc_database": "NACC数据库",
            "search_databases": "🔍 搜索数据库...",
            "search_placeholder": "输入关键词搜索...",
            "searching": "🔎 正在搜索:",
            "found_results": "✅ 找到",
            "results": "个结果",
            "no_results": "⚠️ 未找到结果",
            "try_different": "尝试不同的关键词",
            "database_overview": "数据库概览",
            "data_description": "数据描述",
            "diagnosis_distribution": "诊断分布",
            "data_preview": "数据预览",
            "total_subjects": "总受试者",
            "age_range": "年龄范围",
            "females": "女性",
            "males": "男性",
            "diagnosis": "诊断",
            "count": "数量",
            "gene_query": "基因查询",
            "network_pharmacology": "网络药理学",
            "enter_gene": "输入基因名称:",
            "gene_placeholder": "例如: APP, APOE, MAPT",
            "search_gene": "搜索基因",
            "gene_found": "基因已找到:",
            "gene_details": "基因详情",
            "full_name": "全名:",
            "function": "功能:",
            "association": "关联:",
            "pathway": "通路:",
            "expression": "表达:",
            "risk_level": "风险等级:",
            "gene_not_found": "基因未在数据库中找到",
            "try_genes": "尝试: APP, PSEN1, PSEN2, APOE, MAPT, TREM2, BIN1, CLU, CR1, ABCA7",
            "gene_database": "基因数据库",
            "gene_visualization": "基因可视化",
            "download_image": "下载图片",
            "download_format": "下载格式:",
            "download_png": "PNG (300 DPI)",
            "download_tiff": "TIFF (600 DPI)",
            "downloaded": "✅ 已下载为",
            "downloaded_tiff": "✅ 已下载高分辨率TIFF (600 DPI)",
            "available_genes": "可用基因",
            "tcm_compound_analysis": "中药成分分析",
            "enter_compound": "输入中药成分:",
            "compound_placeholder": "例如: Paeoniflorin, Ferulic acid",
            "search_compound": "搜索成分",
            "compound_found": "成分已找到:",
            "compound_details": "成分详情",
            "chinese_name": "中文名:",
            "source": "来源:",
            "category": "类别:",
            "mechanism": "机制:",
            "evidence_level": "证据等级:",
            "target_genes": "靶基因",
            "ad_relevance": "AD相关性:",
            "compound_not_found": "成分未在数据库中找到",
            "available_compounds": "可用成分",
            "tcm_formula_analysis": "中药方剂分析",
            "enter_formula": "输入中药方剂:",
            "formula_placeholder": "例如: 桃红四物汤",
            "search_formula": "搜索方剂",
            "formula_found": "方剂已找到:",
            "formula_details": "方剂详情",
            "herbal_components": "草本成分:",
            "major_compounds": "主要成分:",
            "clinical_trials": "临床试验:",
            "formula_not_found": "方剂未在数据库中找到",
            "available_formulas": "可用方剂",
            "network_visualization": "网络可视化",
            "network_statistics": "网络统计",
            "compounds": "成分数",
            "unique_targets": "唯一靶点",
            "total_edges": "总边数",
            "avg_targets_compound": "平均靶点/成分",
            "select_compounds": "选择中药成分:",
            "pathway_enrichment": "通路富集分析",
            "select_compounds_enrichment": "选择中药成分或方剂进行通路富集分析",
            "adni_data_loaded": "✅ ADNI数据加载成功",
            "nacc_data_loaded": "✅ NACC数据加载成功",
            "nacc_data_not_available": "⚠️ NACC数据不可用",
            "adni_data_not_available": "⚠️ ADNI数据不可用",
            "data_file_not_found": "数据文件未在预期位置找到"
        },
        "English": {
            "title": "Multi-Omics Alzheimer's Disease Analysis Tool",
            "sidebar_title": "Analysis Options",
            "select_analysis": "Select Analysis Type:",
            "analysis_types": ["Gene Query", "Network Pharmacology", "Public Databases", "Transcriptomics", "Proteomics", "Metabolomics", 
                              "Single-Cell RNA-seq", "Spatial Transcriptomics", "Multi-Omics Integration"],
            "data_upload": "Data Upload",
            "upload_transcriptomics": "Upload Transcriptomics Data (CSV):",
            "upload_proteomics": "Upload Proteomics Data (CSV):",
            "upload_metabolomics": "Upload Metabolomics Data (CSV):",
            "upload_single_cell": "Upload Single-Cell Data (CSV):",
            "upload_spatial": "Upload Spatial Transcriptomics Data (CSV):",
            "public_databases": "Public Databases",
            "adni_database": "ADNI Database",
            "nacc_database": "NACC Database",
            "search_databases": "🔍 Search databases...",
            "search_placeholder": "Enter keyword to search...",
            "searching": "🔎 Searching for:",
            "found_results": "✅ Found",
            "results": "results",
            "no_results": "⚠️ No results found",
            "try_different": "Try different keywords",
            "database_overview": "Database Overview",
            "data_description": "Data Description",
            "diagnosis_distribution": "Diagnosis Distribution",
            "data_preview": "Data Preview",
            "total_subjects": "Total Subjects",
            "age_range": "Age Range",
            "females": "Females",
            "males": "Males",
            "diagnosis": "Diagnosis",
            "count": "Count",
            "gene_query": "Gene Query",
            "network_pharmacology": "Network Pharmacology",
            "enter_gene": "Enter gene name:",
            "gene_placeholder": "e.g., APP, APOE, MAPT",
            "search_gene": "Search Gene",
            "gene_found": "Gene Found:",
            "gene_details": "Gene Details",
            "full_name": "Full Name:",
            "function": "Function:",
            "association": "Association:",
            "pathway": "Pathway:",
            "expression": "Expression:",
            "risk_level": "Risk Level:",
            "gene_not_found": "Gene not found in database",
            "try_genes": "Try: APP, PSEN1, PSEN2, APOE, MAPT, TREM2, BIN1, CLU, CR1, ABCA7",
            "gene_database": "Gene Database",
            "gene_visualization": "Gene Visualization",
            "download_image": "Download Image",
            "download_format": "Download format:",
            "download_png": "PNG (300 DPI)",
            "download_tiff": "TIFF (600 DPI)",
            "downloaded": "✅ Downloaded as",
            "downloaded_tiff": "✅ Downloaded High Resolution TIFF (600 DPI)",
            "available_genes": "Available Genes",
            "tcm_compound_analysis": "TCM Compound Analysis",
            "enter_compound": "Enter TCM compound:",
            "compound_placeholder": "e.g., Paeoniflorin, Ferulic acid",
            "search_compound": "Search Compound",
            "compound_found": "Compound Found:",
            "compound_details": "Compound Details",
            "chinese_name": "Chinese Name:",
            "source": "Source:",
            "category": "Category:",
            "mechanism": "Mechanism:",
            "evidence_level": "Evidence Level:",
            "target_genes": "Target Genes",
            "ad_relevance": "AD Relevance:",
            "compound_not_found": "Compound not found in database",
            "available_compounds": "Available Compounds",
            "tcm_formula_analysis": "TCM Formula Analysis",
            "enter_formula": "Enter TCM formula:",
            "formula_placeholder": "e.g., 桃红四物汤",
            "search_formula": "Search Formula",
            "formula_found": "Formula Found:",
            "formula_details": "Formula Details",
            "herbal_components": "Herbal Components:",
            "major_compounds": "Major Compounds:",
            "clinical_trials": "Clinical Trials:",
            "formula_not_found": "Formula not found in database",
            "available_formulas": "Available Formulas",
            "network_visualization": "Network Visualization",
            "network_statistics": "Network Statistics",
            "compounds": "Compounds",
            "unique_targets": "Unique Targets",
            "total_edges": "Total Edges",
            "avg_targets_compound": "Avg Targets/Compound",
            "select_compounds": "Select TCM Compounds:",
            "pathway_enrichment": "Pathway Enrichment Analysis",
            "select_compounds_enrichment": "Select TCM compounds or formulas to analyze pathway enrichment",
            "adni_data_loaded": "✅ ADNI data loaded successfully",
            "nacc_data_loaded": "✅ NACC data loaded successfully",
            "nacc_data_not_available": "⚠️ NACC data not available",
            "adni_data_not_available": "⚠️ ADNI data not available",
            "data_file_not_found": "Data file not found at expected location"
        }
    }
    
    t = TEXTS[lang]
    
    st.markdown(f"<h1 style='text-align: center; color: #2c3e50; font-family: Roboto, sans-serif; font-size: 28px; margin-bottom: 10px;'>{t['title']}</h1>", unsafe_allow_html=True)
    
    st.markdown("""
    <div class="logo-container">
        <span class="logo-text">Monash University</span>
        <span class="logo-divider">|</span>
        <span class="logo-text">Andy's Lab</span>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("---")
    
    st.sidebar.title(t['sidebar_title'])
    
    analysis_type = st.sidebar.selectbox(
        t['select_analysis'],
        t['analysis_types']
    )
    
    if analysis_type == t['analysis_types'][0]:
        show_gene_query(t)
    elif analysis_type == t['analysis_types'][1]:
        show_network_pharmacology(t)
    elif analysis_type == t['analysis_types'][2]:
        show_public_databases(t)
    
    st.sidebar.markdown("### Data Upload")
    
    data_dict = {}
    
    if analysis_type == "Transcriptomics":
        transcriptomics_file = st.sidebar.file_uploader(
            "Upload Transcriptomics Data (CSV):",
            type=['csv'],
            help="Upload gene expression matrix (genes x samples)"
        )
        
        if transcriptomics_file is not None:
            data_dict['Transcriptomics'] = load_omics_data(transcriptomics_file, "Transcriptomics")
            show_transcriptomics_analysis(data_dict['Transcriptomics'])
    
    elif analysis_type == "Proteomics":
        proteomics_file = st.sidebar.file_uploader(
            "Upload Proteomics Data (CSV):",
            type=['csv'],
            help="Upload protein abundance matrix (proteins x samples)"
        )
        
        if proteomics_file is not None:
            data_dict['Proteomics'] = load_omics_data(proteomics_file, "Proteomics")
            show_proteomics_analysis(data_dict['Proteomics'])
    
    elif analysis_type == "Metabolomics":
        metabolomics_file = st.sidebar.file_uploader(
            "Upload Metabolomics Data (CSV):",
            type=['csv'],
            help="Upload metabolite abundance matrix (metabolites x samples)"
        )
        
        if metabolomics_file is not None:
            data_dict['Metabolomics'] = load_omics_data(metabolomics_file, "Metabolomics")
            show_metabolomics_analysis(data_dict['Metabolomics'])
    
    elif analysis_type == "Single-Cell RNA-seq":
        sc_file = st.sidebar.file_uploader(
            "Upload Single-Cell Data (CSV):",
            type=['csv'],
            help="Upload single-cell expression matrix (cells x genes)"
        )
        
        if sc_file is not None:
            data_dict['Single-Cell'] = load_omics_data(sc_file, "Single-Cell RNA-seq")
            show_single_cell_analysis(data_dict['Single-Cell'])
    
    elif analysis_type == "Spatial Transcriptomics":
        spatial_file = st.sidebar.file_uploader(
            "Upload Spatial Transcriptomics Data (CSV):",
            type=['csv'],
            help="Upload spatial expression matrix (spots x genes)"
        )
        
        if spatial_file is not None:
            data_dict['Spatial'] = load_omics_data(spatial_file, "Spatial Transcriptomics")
            show_spatial_transcriptomics_analysis(data_dict['Spatial'])
    
    elif analysis_type == "Multi-Omics Integration":
        st.sidebar.markdown("#### Upload Multiple Datasets")
        
        transcriptomics_file = st.sidebar.file_uploader(
            "Upload Transcriptomics Data (CSV):",
            type=['csv'],
            key="transcriptomics_upload",
            help="Upload gene expression matrix (genes x samples)"
        )
        
        proteomics_file = st.sidebar.file_uploader(
            "Upload Proteomics Data (CSV):",
            type=['csv'],
            key="proteomics_upload",
            help="Upload protein abundance matrix (proteins x samples)"
        )
        
        metabolomics_file = st.sidebar.file_uploader(
            "Upload Metabolomics Data (CSV):",
            type=['csv'],
            key="metabolomics_upload",
            help="Upload metabolite abundance matrix (metabolites x samples)"
        )
        
        if transcriptomics_file is not None:
            data_dict['Transcriptomics'] = load_omics_data(transcriptomics_file, "Transcriptomics")
        
        if proteomics_file is not None:
            data_dict['Proteomics'] = load_omics_data(proteomics_file, "Proteomics")
        
        if metabolomics_file is not None:
            data_dict['Metabolomics'] = load_omics_data(metabolomics_file, "Metabolomics")
        
        if len(data_dict) > 0:
            show_multi_omics_integration(data_dict)


if __name__ == "__main__":
    main()
