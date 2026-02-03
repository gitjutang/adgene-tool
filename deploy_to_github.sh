#!/bin/bash

# GitHub仓库克隆和部署脚本
# 作者：Monash University Andy's Lab

echo "=========================================="
echo "  阿尔茨海默病在线分析工具 - 部署脚本"
echo "  Monash University | Andy's Lab"
echo "=========================================="
echo ""

# 设置变量
GITHUB_USERNAME="gitjutang"
REPO_NAME="ad-analysis-tool"
LOCAL_DIR="/Users/tomli/mac-rworking/01-梁老师AD疾病-（MR中介分析转录组分析机器学习免疫浸润）/github_repos"
SCRIPT_DIR="/Users/tomli/mac-rworking/01-梁老师AD疾病-（MR中介分析转录组分析机器学习免疫浸润）/001-AD/scripts"

# 进入脚本目录
cd "$SCRIPT_DIR"

echo "📋 步骤1：准备文件"
echo "----------------------------------------"
echo "检查必需文件..."

# 检查文件是否存在
FILES=(
    "26_online_analysis_tool.py"
    "requirements.txt"
    "NACC_filtered_summary.csv"
    ".gitignore"
    "README.md"
    "DEPLOYMENT_GUIDE.md"
    "DEPLOYMENT_CHECKLIST.md"
    "DEPLOYMENT_SUMMARY.md"
)

ALL_FILES_EXIST=true
for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "✅ $file"
    else
        echo "❌ $file (未找到)"
        ALL_FILES_EXIST=false
    fi
done

if [ "$ALL_FILES_EXIST" = false ]; then
    echo ""
    echo "❌ 错误：部分文件缺失，请检查！"
    exit 1
fi

echo ""
echo "📋 步骤2：创建GitHub仓库"
echo "----------------------------------------"
echo "请按照以下步骤操作："
echo ""
echo "1. 访问：https://github.com/new"
echo "2. 填写仓库信息："
echo "   - Repository name: $REPO_NAME"
echo "   - Description: 阿尔茨海默病研究在线分析工具 - Monash University Andy's Lab"
echo "   - Public/Private: 选择 Public（公开）或 Private（私有）"
echo "3. 点击 'Create repository'"
echo ""
read -p "按 Enter 键继续，确认您已创建GitHub仓库..."

echo ""
echo "📋 步骤3：初始化Git仓库"
echo "----------------------------------------"

# 初始化Git仓库
if [ -d ".git" ]; then
    echo "⚠️  Git仓库已存在，将重新初始化"
    rm -rf .git
fi

git init
echo "✅ Git仓库初始化完成"

# 添加远程仓库
git remote add origin "https://github.com/$GITHUB_USERNAME/$REPO_NAME.git"
echo "✅ 远程仓库已添加"

echo ""
echo "📋 步骤4：添加文件到Git"
echo "----------------------------------------"

# 添加所有文件
git add .
echo "✅ 文件已添加到Git"

# 提交更改
git commit -m "Initial commit: 阿尔茨海默病研究在线分析工具"
echo "✅ 文件已提交"

echo ""
echo "📋 步骤5：推送到GitHub"
echo "----------------------------------------"
echo "正在推送到GitHub..."
echo ""

# 推送到GitHub
git branch -M main
git push -u origin main

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "  ✅ 部署成功！"
    echo "=========================================="
    echo ""
    echo "📊 仓库地址："
    echo "   https://github.com/$GITHUB_USERNAME/$REPO_NAME"
    echo ""
    echo "🚀 下一步：部署到Streamlit Cloud"
    echo "----------------------------------------"
    echo "1. 访问：https://share.streamlit.io"
    echo "2. 使用GitHub账号登录"
    echo "3. 点击 'New app'"
    echo "4. 选择仓库：$REPO_NAME"
    echo "5. 选择分支：main"
    echo "6. 主文件路径：26_online_analysis_tool.py"
    echo "7. 点击 'Deploy'"
    echo ""
    echo "等待3-5分钟，部署完成后您将获得公网URL！"
    echo ""
else
    echo ""
    echo "❌ 推送失败，请检查："
    echo "   1. GitHub仓库是否已创建"
    echo "   2. GitHub账号是否已登录"
    echo "   3. 网络连接是否正常"
    echo ""
    echo "手动推送命令："
    echo "   git push -u origin main"
    echo ""
fi
