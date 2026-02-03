# GitHub仓库部署指南

## 📋 准备工作

### 已完成的文件
在 `/Users/tomli/mac-rworking/01-梁老师AD疾病-（MR中介分析转录组分析机器学习免疫浸润）/001-AD/scripts/` 目录下，以下文件已准备好：

✅ `26_online_analysis_tool.py` - 主应用文件
✅ `requirements.txt` - Python依赖包
✅ `NACC_filtered_summary.csv` - ADNI/NACC数据
✅ `.gitignore` - Git配置文件
✅ `README.md` - 项目说明
✅ `DEPLOYMENT_GUIDE.md` - 部署指南
✅ `DEPLOYMENT_CHECKLIST.md` - 文件清单
✅ `DEPLOYMENT_SUMMARY.md` - 部署总结
✅ `deploy_to_github.sh` - 自动部署脚本

## 🚀 两种部署方法

### 方法1：使用GitHub网页界面（推荐，最简单）

#### 步骤1：创建GitHub仓库（2分钟）

1. 打开浏览器，访问：https://github.com/new
2. 填写仓库信息：
   - **Repository name**: `ad-analysis-tool`
   - **Description**: `阿尔茨海默病研究在线分析工具 - Monash University Andy's Lab`
   - **Public/Private**: 选择 `Public`（公开）或 `Private`（私有）
3. 点击 **Create repository**

#### 步骤2：上传文件到GitHub（3分钟）

1. 在新创建的GitHub仓库页面，点击 **uploading an existing file**
2. 打开Finder，导航到：
   ```
   /Users/tomli/mac-rworking/01-梁老师AD疾病-（MR中介分析转录组分析机器学习免疫浸润）/001-AD/scripts/
   ```
3. 选择以下文件，拖拽到GitHub上传区域：
   ```
   26_online_analysis_tool.py
   requirements.txt
   NACC_filtered_summary.csv
   .gitignore
   README.md
   DEPLOYMENT_GUIDE.md
   DEPLOYMENT_CHECKLIST.md
   DEPLOYMENT_SUMMARY.md
   ```
4. 在 "Commit changes" 框中输入：`Initial commit`
5. 点击 **Commit changes**

#### 步骤3：部署到Streamlit Cloud（5分钟）

1. 打开浏览器，访问：https://share.streamlit.io
2. 点击 **Sign up** 或 **Log in**
3. 选择 **Continue with GitHub**
4. 授权Streamlit Cloud访问您的GitHub账号
5. 点击 **New app**
6. 填写应用信息：
   - **Repository**: 选择 `ad-analysis-tool`
   - **Branch**: 选择 `main`
   - **Main file path**: 输入 `26_online_analysis_tool.py`
   - **App URL**: 输入 `monash-andy-ad-tool`（可选）
7. 点击 **Deploy**

#### 步骤4：等待部署完成（3-5分钟）

- Streamlit Cloud会自动构建和部署应用
- 部署完成后，您会看到类似这样的URL：
  ```
  https://monash-andy-ad-tool.streamlit.app
  ```

### 方法2：使用命令行（适合有Git经验的用户）

#### 步骤1：创建GitHub仓库

1. 访问：https://github.com/new
2. 创建名为 `ad-analysis-tool` 的仓库

#### 步骤2：使用自动部署脚本

```bash
# 进入scripts目录
cd /Users/tomli/mac-rworking/01-梁老师AD疾病-（MR中介分析转录组分析机器学习免疫浸润）/001-AD/scripts

# 给脚本添加执行权限
chmod +x deploy_to_github.sh

# 运行部署脚本
./deploy_to_github.sh
```

#### 步骤3：按照脚本提示操作

脚本会自动完成以下步骤：
1. 检查文件是否齐全
2. 初始化Git仓库
3. 添加文件到Git
4. 提交更改
5. 推送到GitHub

#### 步骤4：部署到Streamlit Cloud

按照方法1的步骤3和步骤4操作。

## ✅ 部署成功后的结果

### 公网访问URL
```
https://your-app-name.streamlit.app
```

### 应用功能
- ✅ 4种数据源选择（示例数据、ADNI、NACC、自定义）
- ✅ 5种分析类型（数据概览、生物标志物、预测模型、多组学整合、风险评估）
- ✅ 商标展示（Monash University | Andy's Lab）
- ✅ 任何人都可以通过浏览器访问

### 自动更新
- 更新GitHub代码后，Streamlit Cloud会自动重新部署
- 无需手动操作

## 📞 获取帮助

如果遇到问题：

1. 查看 `DEPLOYMENT_GUIDE.md` 中的故障排除部分
2. 访问 Streamlit 文档：https://docs.streamlit.io
3. 在 Streamlit 社区论坛提问：https://discuss.streamlit.io

## 🎉 开始部署

**推荐使用方法1（GitHub网页界面）**，因为：
- 操作简单直观
- 不需要Git命令
- 适合所有用户

**开始步骤**：
1. 打开 https://github.com/new
2. 创建仓库
3. 上传文件
4. 部署到Streamlit Cloud

10-15分钟后，您的应用就可以在公网访问了！

---

**Monash University | Andy's Lab**
