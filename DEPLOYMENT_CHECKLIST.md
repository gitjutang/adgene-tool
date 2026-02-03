# Streamlit Cloud 部署文件清单

## 📁 需要上传到GitHub的文件

在 `/Users/tomli/mac-rworking/01-梁老师AD疾病-（MR中介分析转录组分析机器学习免疫浸润）/001-AD/scripts/` 目录下，以下文件需要上传：

### 必需文件

1. **26_online_analysis_tool.py** ✅
   - 主应用文件
   - 包含所有分析功能
   - 已集成ADNI和NACC数据

2. **requirements.txt** ✅
   - Python依赖包列表
   - Streamlit Cloud自动安装这些包

3. **NACC_filtered_summary.csv** ✅
   - ADNI/NACC真实数据
   - 包含临床特征和诊断信息

4. **.gitignore** ✅
   - Git忽略文件配置
   - 避免上传不必要的文件

### 可选文件

5. **README.md** ✅
   - 项目说明文档
   - 包含使用说明和部署指南

6. **DEPLOYMENT_GUIDE.md** ✅
   - 快速部署指南
   - 详细的步骤说明

## 📋 文件检查清单

在部署前，请确认：

- [ ] `26_online_analysis_tool.py` 文件存在
- [ ] `requirements.txt` 文件存在
- [ ] `NACC_filtered_summary.csv` 文件存在
- [ ] `.gitignore` 文件存在
- [ ] `README.md` 文件存在（可选）
- [ ] `DEPLOYMENT_GUIDE.md` 文件存在（可选）

## 🔍 文件内容检查

### requirements.txt 内容
```
streamlit==1.53.1
pandas==2.2.0
numpy==1.26.3
matplotlib==3.8.2
seaborn==0.13.1
scikit-learn==1.4.0
```

### .gitignore 关键配置
```
# 忽略Python缓存
__pycache__/
*.py[cod]

# 忽略虚拟环境
venv/
ENV/

# 忽略IDE配置
.vscode/
.idea/

# 忽略Streamlit配置
.streamlit/

# 保留数据文件
!NACC_filtered_summary.csv
```

## 📦 部署步骤总结

1. **创建GitHub仓库**
   - 访问 https://github.com/new
   - 创建名为 `ad-analysis-tool` 的仓库

2. **上传文件**
   - 上传上述6个文件到GitHub仓库

3. **部署到Streamlit Cloud**
   - 访问 https://share.streamlit.io
   - 连接GitHub仓库
   - 选择 `26_online_analysis_tool.py` 作为主文件
   - 点击部署

4. **获取公网URL**
   - 部署完成后获得类似：`https://your-app-name.streamlit.app`

## 🎯 预期结果

部署成功后，您将获得：

- ✅ 一个公网可访问的URL
- ✅ 任何人都可以通过浏览器访问
- ✅ 自动更新功能（更新GitHub代码后自动重新部署）
- ✅ 免费的托管服务

## 📊 应用功能

部署后的应用将提供：

1. **4种数据源选择**
   - 示例数据
   - ADNI数据（真实数据）
   - NACC数据（真实数据）
   - 自定义数据上传

2. **5种分析类型**
   - 数据概览
   - 生物标志物分析
   - 预测模型
   - 多组学整合
   - 风险评估

3. **商标展示**
   - Monash University
   - Andy's Lab

## 🚀 开始部署

现在您可以按照 `DEPLOYMENT_GUIDE.md` 中的步骤开始部署了！

---

**提示**：整个部署过程大约需要10-15分钟。
