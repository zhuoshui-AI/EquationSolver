#!/usr/bin/env python3
"""
非线性方程求解器演示脚本
展示C#方程求解程序中非线性方程求解功能的用法和效果
"""

def print_header(title):
    """打印标题分隔符"""
    print("\n" + "="*60)
    print(f"🚀 {title}")
    print("="*60)

def main():
    print_header("非线性方程求解器功能演示")
    
    print("本演示展示了C#方程求解器中非线性方程求解组件的功能:")
    print("- 牛顿-拉弗森法 (Newton-Raphson)")
    print("- 二分法 (Bisection Method)") 
    print("- 割线法 (Secant Method)")
    print("- 布伦特法 (Brent's Method)")
    print("- 不动点迭代法 (Fixed-Point Iteration)")
    print("- 多项式方程求解器 (Polynomial Solver)")
    print("- 自然语言非线性方程求解器")
    
    print_header("1. 牛顿-拉弗森法示例")
    print("""
用途: 求解具有连续导数的光滑函数的根
优势: 二阶收敛速度，效率高
限制: 需要导数信息，对初值敏感

示例方程: x² - 4 = 0
初始猜测: x₀ = 3.0
期望解: x = 2.0

C#代码示例:
var solver = new NewtonRaphsonSolver()
    .SetFunction(x => x * x - 4, x => 2 * x)
    .WithInitialGuess(3.0)
    .WithTolerance(1e-6)
    .WithMaxIterations(100);
""")

    print_header("2. 二分法示例")
    print("""
用途: 在已知有根区间内求解连续函数的根
优势: 绝对可靠，总能收敛
限制: 收敛速度较慢，需要知道有根区间

示例方程: x² - 2 = 0
搜索区间: [1, 2]
期望解: x ≈ 1.41421356

C#代码示例:
var solver = new BisectionSolver()
    .SetFunction(x => x * x - 2)
    .WithBounds(1, 2)
    .WithTolerance(1e-6)
    .WithMaxIterations(100);
""")

    print_header("3. 割线法示例")
    print("""
用途: 无需导数信息的拟牛顿法
优势: 比二分法快，不需要导数
限制: 收敛性不如牛顿法稳定

示例方程: eˣ - 10 = 0
初始猜测: x₀ = 1.0, x₁ = 3.0
期望解: x ≈ 2.302585

C#代码示例:
var solver = new SecantSolver()
    .SetFunction(x => Math.Exp(x) - 10)
    .WithGuesses(1.0, 3.0)
    .WithTolerance(1e-6)
    .WithMaxIterations(100);
""")

    print_header("4. 布伦特法示例")
    print("""
用途: 结合二分法、割线法和反二次插值的混合方法
优势: 兼具可靠性和高效性
限制: 实现相对复杂

示例方程: sin(x) - x/2 = 0
搜索区间: [-π, π]
期望解: 非平凡解（除了x=0）

C#代码示例:
var solver = new BrentSolver()
    .SetFunction(x => Math.Sin(x) - x/2)
    .WithBounds(-Math.PI, Math.PI)
    .WithTolerance(1e-6)
    .WithMaxIterations(100);
""")

    print_header("5. 多项式方程求解器")
    print("""
支持方法:
- 伴随矩阵法 (Companion Matrix)
- Durand-Kerner方法
- Bairstow方法（实系数多项式）

示例方程: x³ - 6x² + 11x - 6 = 0
期望根: x = 1, 2, 3

C#代码示例:
var poly = new PolynomialSolver(1, -6, 11, -6);
var roots = poly.SolveUsingDurandKerner();
// 或: var roots = poly.SolveUsingCompanionMatrix();
""")

    print_header("6. 自然语言非线性方程求解器")
    print("""
支持的中文输入示例:
- "求解方程：x的平方减四等于零"
- "求sin(x)=0.5的解"
- "解方程e的x次方等于10"

支持的英文输入示例:
- "Solve the equation: x squared minus four equals zero"
- "Find roots of sin(x) = 0.5"
- "Calculate solution for exp(x) = 10"

C#代码示例:
var solver = new NonlinEqNLPSolver()
    .WithParameter("pi", Math.PI)
    .WithParameter("e", Math.E);
    
var result = solver.Solve("求解方程：x^2 - 4 = 0");
""")

    print_header("7. 与MATLAB的比较")
    print("""
MATLAB等效命令 vs C#实现:

MATLAB:
>> f = @(x) x^2 - 4;
>> x = fzero(f, 3)

C#:
var solver = new NewtonRaphsonSolver()
    .SetFunction(x => x*x - 4)
    .WithInitialGuess(3.0);

MATLAB多项式求根:
>> roots([1, -6, 11, -6])

C#多项式求根:
var poly = new PolynomialSolver(1, -6, 11, -6);
var roots = poly.SolveUsingCompanionMatrix();
""")

    print_header("8. 高级特性")
    print("""
🔄 自适应求解策略:
- 系统会自动选择最适合的求解方法
- 牛顿法失败时自动切换到二分法或割线法
- 针对不同类型方程采用专用算法

📊 收敛监控:
- 实时监测迭代过程和收敛状态
- 自动调整步长和容忍度
- 防止无限循环和数值不稳定

🔢 数值稳定性:
- 精心设计的数值微分算法
- 异常情况和边界条件处理
- 浮点数精度优化

🌐 多语言支持:
- 中英文自然语言理解
- 国际化数学符号识别
- 文化敏感的数值格式化
""")

    print_header("9. 实用示例场景")
    
    examples = [
        ("工程计算", "结构力学中的挠度方程", "EI·d²y/dx² = -M(x)", "数值积分和微分方程"),
        ("金融建模", "Black-Scholes期权定价", "偏微分方程求解", "有限差分法"),
        ("物理仿真", "弹簧质点系统", "mx'' + cx' + kx = F(t)", "龙格-库塔法"),
        ("化学平衡", "反应动力学方程", "速率方程和平衡常数", "非线性最小二乘拟合"),
        ("机器学习", "神经网络训练", "梯度下降优化", "自动微分和反向传播")
    ]
    
    for category, application, equation, method in examples:
        print(f"🏷️ {category}: {application}")
        print(f"   方程: {equation}")
        print(f"   方法: {method}\n")

    print_header("10. 性能优化建议")
    print("""
💡 提高求解效率的技巧:

1. 选择合适的初值:
   - 图形分析法确定大致范围
   - 使用粗略扫描定位根的大致位置

2. 调整容忍度设置:
   - 工程应用: 1e-4 ~ 1e-6
   - 科学研究: 1e-8 ~ 1e-12
   - 根据实际需求平衡精度和速度

3. 利用方程特性:
   - 对称性简化计算
   - 单调性指导搜索方向
   - 周期性减少搜索范围

4. 并行计算:
   - 同时求解多个初值
   - 分布式根查找算法
   - GPU加速数值计算
""")

    print_header("演示结束")
    print("🎉 非线性方程求解器已成功集成到C#方程求解程序中！")
    print("📁 相关文件:")
    print("   - NonlinearEquationSolver.cs (核心求解器)")
    print("   - PolynomialSolver.cs (多项式求解器)") 
    print("   - NonlinEqNLPSolver.cs (自然语言求解器)")
    print("   - NonlinearEquationTests.cs (单元测试)")
    print("   - Program.cs (主程序集成)")
    print("\n🚀 下一步: 实现高级矩阵操作功能")

if __name__ == "__main__":
    main()