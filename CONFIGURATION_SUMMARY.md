# EquationSolver 项目配置和命名空间说明

## 项目概述

EquationSolver 是一个功能强大的方程求解系统，支持线性和非线性方程的求解，提供自然语言处理接口，并包含高级矩阵运算功能。

## 命名空间结构说明

### 根命名空间
- `EquationSolver` - 项目根命名空间，包含程序入口点和核心协调类

### 核心功能命名空间

#### `EquationSolver.Interfaces`
定义了项目中使用的所有接口，包括：
- `IEquationSolver` - 方程求解器基础接口
- `ILinearSystemSolver` - 线性方程组求解器接口
- `INonlinearEquationSolver` - 非线性方程求解器接口
- `INaturalLanguageProcessor` - 自然语言处理器接口
- `IMathExpressionParser` - 数学表达式解析器接口

#### `EquationSolver.Models`
包含数据模型和结果封装类：
- `SolveResult` - 求解结果封装
- `ExpressionTree` - 表达式树结构
- `NaturalLanguagePattern` - 自然语言模式识别结果

#### `EquationSolver.Parsers`
包含各种解析器实现：
- `SimpleNaturalLanguageProcessor` - 简单自然语言处理器
- `MathExpressionTokenizer` - 数学表达式词法分析器
- `ShuntingYardParser` - 调度场算法解析器
- `RPNExpressionParser` - 逆波兰表达式解析器

#### `EquationSolver.EquationSolvers`
包含各种方程求解器实现：
- `BaseEquationSolver` - 方程求解器抽象基类
- `MasterEquationSolver` - 主方程求解器（路由不同类型的方程到相应求解器）
- `LinearEquationSolver` - 线性方程求解器
- `QuadraticEquationSolver` - 二次方程求解器
- `GenericEquationSolver` - 通用方程求解器
- `NonlinEqNLPSolver` - 非线性方程自然语言求解器

#### `EquationSolver.EquationSolvers.LinearEquations`
线性方程专用求解器：
- `LinEqNLPSolver` - 线性方程自然语言求解器

#### `EquationSolver.EquationSolvers.NonlinearEquations`
非线性方程专用求解器：
- `NonlinearEquationSolver` - 非线性方程求解器
- `NonlinEqNLPSolver` - 非线性方程自然语言求解器
- `PolynomialSolver` - 多项式求解器

#### `EquationSolver.MatrixOperations`
基础矩阵运算：
- `Matrix` - 矩阵类，提供基本矩阵操作
- 各种矩阵分解结果类（LUDecompositionResult, QRDecompositionResult等）

#### `EquationSolver.AdvancedMatrixOperations`
高级矩阵运算：
- `EigenvalueSolver` - 特征值求解器
- `SvdSolver` - 奇异值分解求解器
- `SparseMatrix` - 稀疏矩阵实现
- `ComprehensiveMatrixAnalyzer` - 综合矩阵分析器

#### `EquationSolver.LinearAlgebra`
线性代数相关功能：
- `LinearSystemSolver` - 线性系统求解器

#### `EquationSolver.Preprocessing`
预处理功能：
- `EquationPreprocessor` - 方程预处理器

## 项目结构

```
EquationSolver/
├── AdvancedMatrixOperations/     高级矩阵运算模块
├── Demos/                        演示代码
├── EquationEngines/              方程引擎
├── EquationSolvers/              方程求解器
│   ├── LinearEquations/         线性方程求解器
│   └── NonlinearEquations/      非线性方程求解器
├── Interfaces/                   接口定义
├── LinearAlgebra/               线性代数模块
├── MatrixOperations/            基础矩阵运算
├── Models/                      数据模型
├── Parsers/                     解析器
├── Preprocessing/               预处理模块
├── Tests/                       测试代码
├── Program.cs                   程序入口点
├── InteractiveConsoleInterface.cs 交互式控制台接口
└── ...
```

## 解决的二义性问题

1. **项目名称与命名空间**：项目名称 "EquationSolver" 与根命名空间相同是 .NET 项目的常见做法，不会造成实际的二义性问题。

2. **命名空间职责明确**：通过上述说明，每个命名空间的职责已经明确定义，避免了功能重叠和未定义的问题。

3. **组件关系清晰**：通过接口和基类的定义，各组件之间的依赖关系和职责分工已经明确。

## 使用建议

1. 对于方程求解功能，建议使用 `MasterEquationSolver` 作为入口点
2. 对于矩阵运算，基础操作使用 `MatrixOperations.Matrix`，高级功能使用 `AdvancedMatrixOperations` 中的专门类
3. 对于自然语言处理，使用 `Parsers.SimpleNaturalLanguageProcessor`
4. 对于特定类型的方程求解，可以直接使用相应的专用求解器