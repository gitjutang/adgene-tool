# ADgene-tool Streamlit Cloud 部署指南

## 📦 部署步骤（5分钟完成）

### 第1步：创建GitHub仓库

1. 访问：https://github.com/new
2. 填写：
   - Repository name: `adgene-tool`
   - Description: `阿尔茨海默病多组学分析工具 - Monash University Andy's Lab`
   - 选择 Public
3. 点击 Create repository

### 第2步：上传文件

将以下3个文件上传到GitHub仓库：
1. `26_online_analysis_tool.py`
2. `requirements.txt`
3. `NACC_filtered_summary.csv`

上传方法：
- 在GitHub仓库页面点击 "uploading an existing file"
- 拖拽文件到上传区域
- 点击 "Commit changes"

### 第3步：部署到Streamlit Cloud

1. 访问：https://share.streamlit.io
2. 点击 "Sign up" 或 "Log in"
3. 选择 "Continue with GitHub"
4. 授权Streamlit访问GitHub
5. 点击 "New app"
6. 填写：
   - Repository: 选择 `adgene-tool`
   - Branch: `main`
   - Main file path: `26_online_analysis_tool.py`
7. 点击 "Deploy"

### 第4步：获取公网地址

等待3-5分钟，部署完成后会显示：
```
https://your-app-name.streamlit.app
```

## ✅ 完成！

现在任何人都可以通过这个URL访问您的ADgene-tool了！

## 🔄 更新应用

修改代码后，在GitHub更新文件，Streamlit Cloud会自动重新部署。

## 📞 需要帮助？

- Streamlit Cloud文档：https://docs.streamlit.io/streamlit-cloud
- Streamlit社区：https://discuss.streamlit.io
