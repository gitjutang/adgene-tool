# 阿尔茨海默病研究在线分析工具

## 📋 项目简介

这是一个基于Streamlit的阿尔茨海默病（AD）研究在线分析工具，集成了ADNI和NACC的真实数据，提供多种分析功能。

## 🎯 主要功能

- **数据概览**：查看数据统计、分布和相关性
- **生物标志物分析**：识别和验证AD相关生物标志物
- **预测模型**：使用机器学习模型进行AD预测
- **多组学整合**：整合转录组、代谢组和蛋白质组数据
- **风险评估**：计算个体的AD风险评分

## 📊 数据源

1. **示例数据**：用于快速演示工具功能
2. **ADNI数据**：包含临床和影像特征
3. **NACC数据**：包含临床和多组学特征
4. **自定义数据**：支持上传CSV文件

## 🚀 部署到Streamlit Cloud

### 步骤1：准备代码

确保以下文件在项目目录中：
- `26_online_analysis_tool.py` - 主应用文件
- `requirements.txt` - Python依赖包
- `NACC_filtered_summary.csv` - 数据文件
- `.gitignore` - Git忽略文件

### 步骤2：创建GitHub仓库

1. 访问 [GitHub](https://github.com) 并登录
2. 点击右上角的 "+" 按钮，选择 "New repository"
3. 填写仓库信息：
   - Repository name: `ad-analysis-tool`
   - Description: `阿尔茨海默病研究在线分析工具`
   - 选择 "Public" 或 "Private"
   - 勾选 "Add a README file"
4. 点击 "Create repository"

### 步骤3：上传代码到GitHub

**方法A：使用GitHub网页界面（推荐新手）**

1. 在GitHub仓库页面，点击 "uploading an existing file"
2. 拖拽以下文件到上传区域：
   - `26_online_analysis_tool.py`
   - `requirements.txt`
   - `NACC_filtered_summary.csv`
   - `.gitignore`
3. 在 "Commit changes" 中填写提交信息
4. 点击 "Commit changes"

**方法B：使用Git命令行**

```bash
# 初始化Git仓库
git init

# 添加所有文件
git add .

# 提交更改
git commit -m "Initial commit"

# 添加远程仓库
git remote add origin https://github.com/YOUR_USERNAME/ad-analysis-tool.git

# 推送到GitHub
git branch -M main
git push -u origin main
```

### 步骤4：部署到Streamlit Cloud

1. 访问 [Streamlit Cloud](https://share.streamlit.io)
2. 点击 "Sign up" 或 "Log in"
3. 选择使用GitHub账号登录
4. 授权Streamlit Cloud访问您的GitHub仓库
5. 点击 "New app"
6. 填写应用信息：
   - Repository: 选择 `ad-analysis-tool`
   - Branch: 选择 `main`
   - Main file path: 输入 `26_online_analysis_tool.py`
7. 点击 "Deploy"

### 步骤5：等待部署完成

- Streamlit Cloud会自动构建和部署应用
- 通常需要2-5分钟
- 部署完成后，您会获得一个公网URL，例如：
  - `https://your-app-name.streamlit.app`

## 📝 自定义配置

### 修改应用标题

编辑 `26_online_analysis_tool.py` 文件：

```python
st.title("🧠 阿尔茨海默病研究在线分析工具")
```

### 修改商标

编辑 `26_online_analysis_tool.py` 文件中的商标部分：

```python
st.markdown("""
<div class="logo-container">
    <span class="logo-text">Monash University</span>
    <span class="logo-divider">|</span>
    <span class="logo-text">Andy's Lab</span>
</div>
""", unsafe_allow_html=True)
```

### 添加更多数据

将CSV文件放在项目目录中，然后修改 `load_adni_data()` 或 `load_nacc_data()` 函数。

## 🔧 本地运行

### 安装依赖

```bash
pip install -r requirements.txt
```

### 运行应用

```bash
streamlit run 26_online_analysis_tool.py
```

应用将在浏览器中自动打开：http://localhost:8501

## 📧 联系方式

如有问题，请联系：
- Monash University
- Andy's Lab

## 📄 许可证

本项目仅供学术研究使用。

---

**注意**：部署到Streamlit Cloud后，任何人都可以通过公网URL访问您的应用。如果需要限制访问，请在Streamlit Cloud中设置密码保护。
