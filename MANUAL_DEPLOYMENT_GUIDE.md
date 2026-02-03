# 手动部署指南

由于网络连接问题，请按照以下步骤手动部署到 GitHub：

## 📋 步骤 1：访问 GitHub 仓库

1. 打开浏览器，访问：https://github.com/gitjutang/ad-analysis-tool
2. 如果仓库不存在，请先创建：
   - 访问：https://github.com/new
   - Repository name: `ad-analysis-tool`
   - Description: `阿尔茨海默病研究在线分析工具 - Monash University Andy's Lab`
   - 选择 Public 或 Private
   - 点击 Create repository

## 📋 步骤 2：上传文件到 GitHub

### 方法 A：使用 GitHub 网页界面上传（推荐）

1. 在 GitHub 仓库页面，点击 **uploading an existing file** 链接
2. 打开 Finder，导航到：
   ```
   /Users/tomli/mac-rworking/01-梁老师AD疾病-（MR中介分析转录组分析机器学习免疫浸润）/001-AD/scripts/
   ```
3. 选择以下文件，拖拽到 GitHub 上传区域：

**必需文件：**
- ✅ `26_online_analysis_tool.py` - 主应用文件
- ✅ `requirements.txt` - Python 依赖包
- ✅ `NACC_filtered_summary.csv` - ADNI/NACC 数据
- ✅ `.gitignore` - Git 配置文件

**文档文件：**
- ✅ `README.md` - 项目说明
- ✅ `DEPLOYMENT_GUIDE.md` - 部署指南
- ✅ `DEPLOYMENT_CHECKLIST.md` - 文件清单
- ✅ `DEPLOYMENT_SUMMARY.md` - 部署总结
- ✅ `GITHUB_DEPLOYMENT_GUIDE.md` - GitHub 部署指南
- ✅ `DEPLOYMENT_OVERVIEW.md` - 部署概览
- ✅ `deploy_to_github.sh` - 部署脚本

**结果图片（可选）：**
- ✅ `results/ADNI_longitudinal/ADNI_longitudinal_analysis.png`
- ✅ `results/AI_prediction_models/AI_prediction_models.png`
- ✅ `results/biomarker_panel/biomarker_panel.png`
- ✅ `results/four_dimensional_integration/four_dimensional_integration.png`
- ✅ `results/multi_dataset_integration/multi_dataset_integration.png`
- ✅ `results/multi_omics_integration/multi_omics_integration.png`

4. 在 "Commit changes" 框中输入：`Initial commit: 阿尔茨海默病研究在线分析工具 - Monash University Andy's Lab`
5. 点击 **Commit changes**

### 方法 B：使用 GitHub Desktop（如果已安装）

1. 打开 GitHub Desktop 应用
2. 选择 **File > Clone Repository**
3. 输入仓库 URL：`https://github.com/gitjutang/ad-analysis-tool`
4. 选择本地路径，点击 Clone
5. 将文件复制到克隆的仓库中
6. 在 GitHub Desktop 中提交并推送

## 📋 步骤 3：验证文件上传成功

1. 访问：https://github.com/gitjutang/ad-analysis-tool
2. 检查以下文件是否存在：
   - `26_online_analysis_tool.py`
   - `requirements.txt`
   - `NACC_filtered_summary.csv`
   - `.gitignore`
   - `README.md`

## 📋 步骤 4：部署到 Streamlit Cloud

### 4.1 访问 Streamlit Cloud

1. 打开浏览器，访问：https://share.streamlit.io
2. 点击 **Sign up** 或 **Log in**
3. 选择 **Continue with GitHub**
4. 授权 Streamlit Cloud 访问您的 GitHub 账号

### 4.2 创建新应用

1. 点击 **New app** 按钮
2. 填写应用信息：

   **Repository（仓库）：**
   - 选择：`ad-analysis-tool`

   **Branch（分支）：**
   - 选择：`main`

   **Main file path（主文件路径）：**
   - 输入：`26_online_analysis_tool.py`

   **App URL（应用 URL，可选）：**
   - 输入：`monash-andy-ad-tool`
   - 最终 URL 将是：`https://monash-andy-ad-tool.streamlit.app`

3. 点击 **Deploy** 按钮

### 4.3 等待部署完成

- Streamlit Cloud 会自动：
  1. 克隆您的 GitHub 仓库
  2. 安装 requirements.txt 中的依赖包
  3. 启动 Streamlit 应用
  4. 生成公网访问 URL

- 部署时间：3-5 分钟

### 4.4 访问您的应用

部署成功后，您会看到：
- ✅ 应用状态：Running
- ✅ 应用 URL：`https://monash-andy-ad-tool.streamlit.app`

点击 URL 即可在浏览器中访问您的在线分析工具！

## 📋 步骤 5：测试应用功能

访问应用后，测试以下功能：

1. **数据源选择**
   - ✅ 示例数据
   - ✅ ADNI 数据
   - ✅ NACC 数据
   - ✅ 自定义数据上传

2. **分析类型**
   - ✅ 数据概览
   - ✅ 生物标志物分析
   - ✅ 预测模型
   - ✅ 多组学整合
   - ✅ 风险评估

3. **商标显示**
   - ✅ Monash University | Andy's Lab

## 📋 步骤 6：分享应用

部署成功后，您可以：

1. **分享 URL**
   - 将 `https://monash-andy-ad-tool.streamlit.app` 分享给其他人
   - 任何人都可以通过浏览器访问

2. **更新应用**
   - 在 GitHub 上更新代码
   - Streamlit Cloud 会自动重新部署
   - 无需手动操作

## ❓ 常见问题

### Q1: 文件上传失败怎么办？

**A:** 
- 检查文件大小（GitHub 单文件限制 100MB）
- 使用 GitHub Desktop 或 Git 命令行工具
- 将大文件压缩后上传

### Q2: Streamlit Cloud 部署失败怎么办？

**A:**
- 检查 requirements.txt 是否正确
- 查看 Streamlit Cloud 的部署日志
- 确保主文件路径正确：`26_online_analysis_tool.py`

### Q3: 如何更新应用？

**A:**
1. 在 GitHub 上修改文件
2. 提交更改
3. Streamlit Cloud 会自动检测并重新部署

### Q4: 应用可以私有吗？

**A:**
- Streamlit Cloud 支持私有仓库
- 需要付费订阅（免费版仅支持公开仓库）
- 私有应用需要登录才能访问

## 🎉 完成！

恭喜！您的阿尔茨海默病在线分析工具已经部署成功！

**应用地址：** https://monash-andy-ad-tool.streamlit.app

**GitHub 仓库：** https://github.com/gitjutang/ad-analysis-tool

---

**Monash University | Andy's Lab**
