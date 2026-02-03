# ADgene-tool Streamlit Cloud 部署总结

## 📋 准备好的文件

✅ `26_online_analysis_tool.py` - 主应用文件（已更新为相对路径）
✅ `requirements.txt` - Python依赖包
✅ `NACC_filtered_summary.csv` - NACC数据文件
✅ `.gitignore` - Git忽略文件配置
✅ `QUICK_DEPLOY.md` - 快速部署指南

## 🚀 部署步骤（5分钟）

### 第1步：创建GitHub仓库

1. 访问：https://github.com/new
2. 填写信息：
   - **Repository name**: `adgene-tool`
   - **Description**: `阿尔茨海默病多组学分析工具 - Monash University Andy's Lab`
   - **Public/Private**: 选择 `Public`
3. 点击 **Create repository**

### 第2步：上传文件到GitHub

**方法A：网页上传（推荐）**

1. 在GitHub仓库页面，点击 **uploading an existing file**
2. 拖拽以下3个文件到上传区域：
   - `26_online_analysis_tool.py`
   - `requirements.txt`
   - `NACC_filtered_summary.csv`
3. 在 "Commit changes" 框中输入：`Initial commit`
4. 点击 **Commit changes**

**方法B：命令行上传**

```bash
cd /Users/tomli/mac-rworking/01-梁老师AD疾病-（MR中介分析转录组分析机器学习免疫浸润）/001-AD/scripts

# 初始化Git
git init

# 添加文件
git add 26_online_analysis_tool.py
git add requirements.txt
git add NACC_filtered_summary.csv
git add .gitignore

# 提交
git commit -m "Initial commit"

# 添加远程仓库（替换YOUR_USERNAME）
git remote add origin https://github.com/YOUR_USERNAME/adgene-tool.git

# 推送到GitHub
git branch -M main
git push -u origin main
```

### 第3步：部署到Streamlit Cloud

1. 访问：https://share.streamlit.io
2. 点击 **Sign up** 或 **Log in**
3. 选择 **Continue with GitHub**
4. 授权Streamlit Cloud访问您的GitHub仓库
5. 点击 **New app**
6. 填写应用信息：
   - **Repository**: 选择 `adgene-tool`
   - **Branch**: 选择 `main`
   - **Main file path**: 输入 `26_online_analysis_tool.py`
   - **App URL**: 输入 `adgene-tool`（可选，自定义URL）
7. 点击 **Deploy**

### 第4步：等待部署完成

- Streamlit Cloud会自动构建和部署
- 通常需要3-5分钟
- 部署完成后，您会看到类似这样的URL：
  ```
  https://adgene-tool.streamlit.app
  ```

## ✅ 部署成功！

现在任何人都可以通过这个URL访问您的ADgene-tool了！

## 🎯 功能特性

- ✅ 中英文语言切换
- ✅ 基因查询系统
- ✅ 网络药理学分析
- ✅ 公共数据库（ADNI、NACC）
- ✅ 多组学数据整合
- ✅ 可下载的高质量图片（PNG、TIFF 300/600 DPI）

## 🔄 更新应用

如果需要更新应用：

1. 在本地修改 `26_online_analysis_tool.py`
2. 在GitHub仓库中更新文件
3. Streamlit Cloud会自动检测到更改并重新部署

## 📞 常见问题

### Q1: 部署失败怎么办？

**A**: 检查以下几点：
- `requirements.txt` 格式是否正确
- 代码是否有语法错误
- 在本地运行 `streamlit run 26_online_analysis_tool.py` 测试

### Q2: 数据文件未找到？

**A**: 确保 `NACC_filtered_summary.csv` 已上传到GitHub仓库根目录

### Q3: 应用运行缓慢？

**A**: 这是正常的，Streamlit Cloud免费版有性能限制。可以考虑：
- 优化数据处理逻辑
- 使用数据缓存
- 升级到付费计划

### Q4: 如何限制访问？

**A**: 在Streamlit Cloud中可以设置密码保护：
1. 进入应用设置
2. 启用 "Require password"
3. 设置访问密码

## 📞 获取帮助

- Streamlit Cloud文档：https://docs.streamlit.io/streamlit-cloud
- Streamlit社区论坛：https://discuss.streamlit.io
- GitHub仓库：https://github.com/YOUR_USERNAME/adgene-tool

## 🎉 完成！

恭喜您成功部署了ADgene-tool！

---

**注意**：
- Streamlit Cloud免费版提供无限公开应用
- 每月750小时的运行时间
- 基本的技术支持
- 如需更多功能，可升级到付费计划
