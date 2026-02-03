"""
Streamlit应用启动器
解决权限问题的简化版本
"""

import subprocess
import sys
import os

def main():
    """
    启动Streamlit应用
    """
    print("=" * 60)
    print("🚀 启动阿尔茨海默病研究在线分析工具")
    print("=" * 60)
    print()
    
    # 获取脚本目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    app_script = os.path.join(script_dir, "26_online_analysis_tool.py")
    
    # 检查脚本是否存在
    if not os.path.exists(app_script):
        print(f"❌ 错误: 找不到应用脚本 {app_script}")
        return
    
    print(f"📁 应用脚本: {app_script}")
    print(f"📍 访问地址: http://localhost:8501")
    print()
    print("💡 提示：")
    print("  - 应用将在浏览器中自动打开")
    print("  - 如果浏览器没有自动打开，请手动访问上面的地址")
    print("  - 要停止应用，请按 Ctrl+C")
    print()
    print("=" * 60)
    print()
    
    # 启动Streamlit应用
    try:
        subprocess.run([
            sys.executable,
            "-m", "streamlit", "run", app_script,
            "--server.port", "8501",
            "--browser.gatherUsageStats", "false"
        ])
    except KeyboardInterrupt:
        print("\n\n✅ 应用已停止")
    except Exception as e:
        print(f"\n\n❌ 启动失败: {e}")

if __name__ == "__main__":
    main()
