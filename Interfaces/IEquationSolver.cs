using EquationSolver.Models;

namespace EquationSolver.Interfaces
{
    /// <summary>
    /// 方程求解器基础接口
    /// </summary>
    public interface IEquationSolver
    {
        /// <summary>
        /// 验证方程格式是否正确
        /// </summary>
        bool ValidateEquation(string equation);
        
        /// <summary>
        /// 解析方程字符串为内部表示
        /// </summary>
        void ParseEquation(string equation);
        
        /// <summary>
        /// 解方程并返回结果
        /// </summary>
        SolveResult Solve();
        
        /// <summary>
        /// 获取支持的方程类型描述
        /// </summary>
        string GetSupportedTypes();
    }

    /// <summary>
    /// 线性方程组求解器接口
    /// </summary>
    public interface ILinearSystemSolver : IEquationSolver
    {
        /// <summary>
        /// 设置系数矩阵
        /// </summary>
        void SetCoefficientMatrix(double[,] matrix);
        
        /// <summary>
        /// 设置常数项向量
        /// </summary>
        void SetConstantVector(double[] vector);
        
        /// <summary>
        /// 高斯消元法求解
        /// </summary>
        double[] SolveByGaussianElimination();
        
        /// <summary>
        /// LU分解法求解
        /// </summary>
        double[] SolveByLUDecomposition();
    }

    /// <summary>
    /// 非线性方程求解器接口
    /// </summary>
    public interface INonlinearEquationSolver : IEquationSolver
    {
        /// <summary>
        /// 牛顿迭代法求解
        /// </summary>
        double SolveByNewtonRaphson(double initialGuess, double tolerance = 1e-15, int maxIterations = 1000);
        
        /// <summary>
        /// 二分法求解
        /// </summary>
        double SolveByBisection(double leftBound, double rightBound, double tolerance = 1e-14);
        
        /// <summary>
        /// 割线法求解
        /// </summary>
        double SolveBySecantMethod(double x0, double x1, double tolerance = 1e-13);
    }
}

// 数学表达式解析器接口
public interface IMathExpressionParser
{
    EquationSolver.Models.ExpressionTree Parse(string expression);
    Dictionary<string, double> ExtractVariables(string expression);
    bool ValidateSyntax(string expression);
}
