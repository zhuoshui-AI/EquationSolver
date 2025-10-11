# EquationSolver - 智能方程求解器

## 项目概述

这是一个用C#开发的智能方程求解程序，支持多种方程形式的自动识别与求解，具备类似于MATLAB的高度优化矩阵操作能力，并且支持自然语言输入。

## 功能特性

### 🧮 支持的方程类型
- **一元一次方程**: ax + b = 0
- **一元二次方程**: ax² + bx + c = 0  
- **多项式方程**: 高次多项式求解
- **非线性方程**: 使用数值方法求解
- **三角函数方程**
- **指数和对数方程**

### 🤖 自然语言处理
- 支持中文自然语言输入："求方程 x² + 2x - 3 = 0 的解"
- 自动识别方程类型和变量
- 智能推荐合适的求解算法

### ⚡ 高性能计算
- 优化的矩阵运算库
- 多重精度数值计算
- 并行化算法实现

## 项目结构

```
EquationSolver/
├── Interfaces/           # 接口定义
│   ├── IEquationSolver.cs          # 方程求解器接口
│   └── INaturalLanguageProcessor.cs # 自然语言处理接口
├── Models/              # 数据模型
│   └── MathExpressions.cs          # 数学表达式模型
├── Parsers/             # 解析器
│   ├── SimpleNaturalLanguageProcessor.cs # 自然语言处理器
│   └── ShuntingYardParser.cs       # 数学表达式解析器
├── EquationSolvers/     # 求解器实现
│   ├── BaseEquationSolver.cs       # 抽象基类
│   └── BasicEquationSolvers.cs     # 基础求解器
├── MatrixOperations/    # 矩阵运算模块
└── Tests/               # 测试用例
```

## 快速开始

### 使用方法

```bash
# 克隆项目后进入目录
cd EquationSolver

# 编译项目
dotnet build

# 运行程序
dotnet run
```

### 示例输入

**自然语言输入:**
```
请输入方程或指令 ('quit'退出, 'help'帮助): 
求方程 x² + 2x - 3 = 0 的解
```

**数学表达式输入:**
```
请输入方程或指令 ('quit'退出, 'help'帮助): 
solve(x^2 + 2*x - 3 = 0)
```

**简单方程输入:**
```
请输入方程或指令 ('quit'退出, 'help'帮助): 
2x + 5 = 16
```

## 技术架构

### 核心组件

1. **自然语言处理器 (NLP)** - 将自然语言转换为数学表达式
2. **数学表达式解析器** - 使用调度场算法解析数学表达式
3. **方程路由器** - 根据方程类型分配合适的求解器
4. **专用求解器** - 针对不同类型方程的优化算法
5. **矩阵运算引擎** - 高效的线性代数计算

### 算法特点

- **自适应求解策略**: 根据方程复杂度自动选择最优算法
- **数值稳定性**: 内置误差控制和边界检查
- **性能优化**: 利用现代CPU特性和并行计算
- **扩展性设计**: 易于添加新的方程类型和求解算法

## 依赖说明

本项目使用以下技术和库：

- .NET 8.0 Runtime
- 标准数学库 (System.Math)
- 正则表达式引擎 (System.Text.RegularExpressions)

## 开发指南

### 添加新方程类型

1. 在 `Interfaces/IEquationSolver.cs` 中添加新接口
2. 在 `EquationSolvers/` 下实现对应的求解器类
3. 更新 `MasterEquationSolver` 的路由逻辑

### 扩展自然语言理解

修改 `SimpleNaturalLanguageProcessor.cs` 中的模式识别规则和词汇表。

## 许可证

MIT License

## 贡献

欢迎提交Issue和Pull Request来改进这个项目！

---

*最后更新时间: 2025年10月*
