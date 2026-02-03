#!/bin/bash

# Streamlit应用启动脚本

echo "🚀 启动阿尔茨海默病研究在线分析工具..."

# 设置环境变量
export STREAMLIT_SERVER_PORT=8502
export STREAMLIT_SERVER_ADDRESS=localhost
export STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

# 启动Streamlit应用
echo "📱 应用将在浏览器中自动打开..."
echo "📍 访问地址: http://localhost:8502"
echo ""
echo "💡 提示："
echo "  - 如果浏览器没有自动打开，请手动访问上面的地址"
echo "  - 要停止应用，请按 Ctrl+C"
echo ""

# 切换到脚本目录
cd "$(dirname "$0")"

# 启动应用
streamlit run 26_online_analysis_tool.py
