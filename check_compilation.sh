#!/bin/bash
# 切换到项目目录
cd /home/luck/git/poly-commit

# 运行 cargo check 来验证修复
echo "正在检查编译错误..."
cargo check

echo "检查完成！"
