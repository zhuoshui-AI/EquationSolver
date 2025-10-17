using System;
using System.Diagnostics;
using System.Text.RegularExpressions;
using EquationSolver.Interfaces;
using EquationSolver.Models;

namespace EquationSolver.EquationSolvers
{
    /// <summary>
    /// 方程求解器抽象基类
    /// </summary>
    public abstract class BaseEquationSolver : IEquationSolver
    {
        protected string CurrentEquation { get; private set; } = string.Empty;
        protected ExpressionTree ParsedExpression { get; private set; }
        protected readonly Regex ValidEquationRegex = new Regex(@"^[\dxX\+\-\*\/\^\=\(\)\.\,\s]+$", RegexOptions.Compiled);

        /// <summary>
        /// 验证方程格式
        /// </summary>
        public virtual bool ValidateEquation(string equation)
        {
            if (string.IsNullOrWhiteSpace(equation))
                return false;

            // 基本字符集验证
            if (!ValidEquationRegex.IsMatch(equation.Replace(" ", "")))
                return false;

            // 等式符号数量验证
            var equalsCount = CountOccurrences(equation, '=');
            if (equalsCount != 1)
                return false;

            return true;
        }

        /// <summary>
        /// 解析方程
        /// </summary>
        public virtual void ParseEquation(string equation)
        {
            if (!ValidateEquation(equation))
                throw new ArgumentException("无效的方程格式");

            CurrentEquation = equation.Trim();
            
            // 分离左右两边
            var parts = equation.Split('=', 2);
            if (parts.Length != 2)
                throw new ArgumentException("方程必须包含且仅包含一个等号");

            var leftSide = parts[0].Trim();
            var rightSide = parts[1].Trim();

            // 移项到左边形成 f(x) = 0 的形式
            var normalizedEquation = NormalizeEquation(leftSide, rightSide);
            
            // 具体的解析逻辑由子类实现
            PerformSpecificParsing(normalizedEquation);
        }

        /// <summary>
        /// 解方程
        /// </summary>
        public abstract SolveResult Solve();

        /// <summary>
        /// 获取支持的方程类型
        /// </summary>
        public abstract string GetSupportedTypes();

        /// <summary>
        /// 执行具体解析逻辑
        /// </summary>
        protected abstract void PerformSpecificParsing(string normalizedEquation);

        /// <summary>
        /// 规范化方程到 f(x) = 0 形式
        /// </summary>
        protected virtual string NormalizeEquation(string leftSide, string rightSide)
        {
            if (rightSide == "0")
                return leftSide;

            // 如果右边不是零，则移到左边
            return $"{leftSide} - ({rightSide})";
        }

        /// <summary>
        /// 统计字符出现次数
        /// </summary>
        protected int CountOccurrences(string text, char character)
        {
            int count = 0;
            foreach (char c in text)
            {
                if (c == character) count++;
            }
            return count;
        }

        /// <summary>
        /// 计时执行方法
        /// </summary>
        protected T MeasureExecution<T>(Func<T> func, out TimeSpan elapsedTime)
        {
            var stopwatch = Stopwatch.StartNew();
            try
            {
                return func();
            }
            finally
            {
                stopwatch.Stop();
                elapsedTime = stopwatch.Elapsed;
            }
        }

        /// <summary>
        /// 检查数值稳定性
        /// </summary>
        protected virtual bool CheckNumericalStability(double[] values, double tolerance = 1e-15)
        {
            foreach (var value in values)
            {
                if (double.IsNaN(value) || double.IsInfinity(value))
                    return false;
                    
                // 检查数值是否超出双精度浮点数的表示范围
                const double MAX_DOUBLE_VALUE = double.MaxValue;
                if (Math.Abs(value) > MAX_DOUBLE_VALUE)
                    return false;
            }
            return true;
        }
    }

    /// <summary>
    /// 主方程求解器 - 负责路由不同类型的方程到相应的专门求解器
    /// </summary>
    public class MasterEquationSolver
    {
        private readonly INaturalLanguageProcessor _nlpProcessor;
        private readonly IMathExpressionParser _mathParser;

        public MasterEquationSolver(INaturalLanguageProcessor nlpProcessor, IMathExpressionParser mathParser)
        {
            _nlpProcessor = nlpProcessor ?? throw new ArgumentNullException(nameof(nlpProcessor));
            _mathParser = mathParser ?? throw new ArgumentNullException(nameof(mathParser));
        }

        /// <summary>
        /// 异步求解方程
        /// </summary>
        public async Task<SolveResult> SolveAsync(string input)
        {
            // 判断输入是自然语言还是纯数学表达式
            var equationType = DetectInputType(input);
            
            string mathematicalExpression;
            Dictionary<string, double> parameters = new Dictionary<string, double>();

            if (equationType == InputType.NaturalLanguage)
            {
                // 使用NLP处理器转换自然语言
                var nlpResult = _nlpProcessor.ConvertToMathematicalNotation(
                    _nlpProcessor.RecognizeEquationPattern(_nlpProcessor.PreprocessText(input)));
                
                mathematicalExpression = nlpResult.MathematicalExpression;
                parameters = nlpResult.InitialValues;
            }
            else
            {
                mathematicalExpression = input;
            }

            // 解析数学表达式
            var expressionTree = _mathParser.Parse(mathematicalExpression);
            
            // 根据方程复杂度选择适当的求解器
            var solver = SelectAppropriateSolver(expressionTree, mathematicalExpression);
            
            // 执行求解
            return await Task.Run(() => solver.Solve());
        }

        /// <summary>
        /// 检测输入类型
        /// </summary>
        private InputType DetectInputType(string input)
        {
            // 简单的启发式规则来判断是否是自然语言
            var naturalLanguageIndicators = new[]
            {
                "求", "解", "方程", "等于", "多少", "什么", "如何", "怎么",
                "calculate", "solve", "equation", "equal", "what", "how"
            };

            var containsChinese = ContainsChineseCharacters(input);
            var containsKeywords = naturalLanguageIndicators.Any(keyword => 
                input.Contains(keyword, StringComparison.OrdinalIgnoreCase));

            return (containsChinese || containsKeywords) ? 
                   InputType.NaturalLanguage : InputType.MathematicalExpression;
        }

        /// <summary>
        /// 检查是否包含中文字符
        /// </summary>
        private bool ContainsChineseCharacters(string text)
        {
            return text.Any(c => c >= 0x4E00 && c <= 0x9FFF);
        }

        /// <summary>
        /// 选择合适的求解器
        /// </summary>
        private IEquationSolver SelectAppropriateSolver(ExpressionTree expressionTree, string expression)
        {
            var variables = expressionTree.Variables;
            var complexity = EstimateComplexity(expressionTree);
            var equationType = DetectEquationCategory(expression);

            // 基于方程类别、复杂度和变量数量的智能路由逻辑
            if (variables.Count == 1)
            {
                switch (complexity)
                {
                    case Complexity.Simple:
                    case Complexity.Linear:
                        return new LinearEquationSolver();
                        
                    case Complexity.Quadratic:
                        return new QuadraticEquationSolver();
                        
                    case Complexity.Polynomial:
                        // 对于多项式，优先使用专门的求解器
                        return CreatePolynomialSolverIfPossible(expression) ?? 
                               new GenericEquationSolver(_mathParser);
                               
                    case Complexity.Complex:
                        // 复杂方程使用非线性求解器
                        return new NonlinEqNLPSolver().WithParameters(GetDefaultParameters());
                }
            }
            
            // 多变量情况
            // 当变量数量较多时（这里设定阈值为100），考虑使用线性方程组求解器
            const int MULTI_VARIABLE_THRESHOLD = 100;
            if (variables.Count > MULTI_VARIABLE_THRESHOLD)
            {
                // 如果是线性方程组，使用线性求解器
                if (IsLinearSystem(expressionTree))
                {
                    return new LinEqNLPSolver();
                }
            }
            
            // 默认情况下使用非线性自然语言求解器
            return new NonlinEqNLPSolver().WithParameters(GetDefaultParameters());
        }

        /// <summary>
        /// 估计方程复杂度
        /// </summary>
        private Complexity EstimateComplexity(ExpressionTree expressionTree)
        {
            // 改进的复杂度估算逻辑
            var expressionStr = expressionTree.ToString().ToLower();
            
            // 检测复杂函数
            var complexFunctions = new[] { "sin(", "cos(", "tan(", "exp(", "log(", "ln(", "sqrt(" };
            if (complexFunctions.Any(func => expressionStr.Contains(func)))
                return Complexity.Complex;
            
            // 检测多项式
            if (expressionStr.Contains("^") || Regex.IsMatch(expressionStr, @"x\^?\d+"))
                return Complexity.Polynomial;
                
            // 检测二次项
            if (expressionStr.Contains("x²") || expressionStr.Contains("x*x") || 
                Regex.IsMatch(expressionStr, @"x\*\s*x"))
                return Complexity.Quadratic;
                
            // 检测线性关系
            if (expressionStr.Contains("*") || expressionStr.Contains("/"))
                return Complexity.Linear;
                
            return Complexity.Simple;
        }

        /// <summary>
        /// 检测方程类别
        /// </summary>
        private EquationCategory DetectEquationCategory(string expression)
        {
            var lowerExpr = expression.ToLower();
            
            if (lowerExpr.Contains("sin") || lowerExpr.Contains("cos") || lowerExpr.Contains("tan"))
                return EquationCategory.Trigonometric;
                
            if (lowerExpr.Contains("exp") || lowerExpr.Contains("log") || lowerExpr.Contains("ln"))
                return EquationCategory.ExponentialLogarithmic;
                
            if (Regex.IsMatch(lowerExpr, @"x\^\d+") || lowerExpr.Contains("x²") || lowerExpr.Contains("x³"))
                return EquationCategory.Polynomial;
                
            return EquationCategory.General;
        }

        /// <summary>
        /// 如果可能的话创建多项式求解器
        /// </summary>
        private IEquationSolver CreatePolynomialSolverIfPossible(string expression)
        {
            try
            {
                var coefficients = ExtractPolynomialCoefficients(expression);
                // 当多项式阶数较高时（这里设定阈值为100），使用通用求解器
                const int HIGH_ORDER_POLYNOMIAL_THRESHOLD = 100;
                if (coefficients != null && coefficients.Length > HIGH_ORDER_POLYNOMIAL_THRESHOLD)
                {
                    // 这里可以进一步封装成适配器模式
                    return new SimpleGenericEquationSolver(_mathParser); // 临时方案
                }
            }
            catch
            {
                // 如果不能提取系数，返回null使用默认求解器
            }
            return null;
        }

        /// <summary>
        /// 检查是否为线性系统
        /// </summary>
        private bool IsLinearSystem(ExpressionTree expressionTree)
        {
            // 简化的线性系统检测
            var expressionStr = expressionTree.ToString();
            return !expressionStr.Contains("^") && 
                   !expressionStr.Contains("sin") &&
                   !expressionStr.Contains("cos") &&
                   !expressionStr.Contains("exp");
        }

        /// <summary>
        /// 获取默认参数
        /// </summary>
        private Dictionary<string, double> GetDefaultParameters()
        {
            return new Dictionary<string, double>
            {
                { "pi", Math.PI },
                { "e", Math.E }
            };
        }

        /// <summary>
        /// 从表达式中提取多项式系数
        /// </summary>
        private double[] ExtractPolynomialCoefficients(string expression)
        {
            // 简化的系数提取逻辑
            // 在实际实现中应该使用更复杂的解析器
            try
            {
                var matches = Regex.Matches(expression, @"([+-]?\d*\.?\d*)\*?x\^?(\d*)");
                // 这里应该是检查匹配数量是否合理，而不是与一个魔数比较
                if (matches.Count == 0)
                    return null;
                    
                var maxPower = 0;
                foreach (Match match in matches)
                {
                    // 正确访问正则表达式的组
                    if (match.Groups.Count > 2 && !string.IsNullOrEmpty(match.Groups[2].Value))
                    {
                        maxPower = Math.Max(maxPower, int.Parse(match.Groups[2].Value));
                    }
                }
                
                var coefficients = new double[maxPower + 1];
                
                foreach (Match match in matches)
                {
                    // 正确访问正则表达式的组
                    var coefStr = match.Groups[1].Value;
                    var powerStr = match.Groups[2].Value;
                    
                    // 修正系数解析逻辑
                    double coefficient = 1.0;
                    if (!string.IsNullOrEmpty(coefStr))
                    {
                        if (coefStr == "+")
                            coefficient = 1.0;
                        else if (coefStr == "-")
                            coefficient = -1.0;
                        else
                            coefficient = double.Parse(coefStr);
                    }
                                       
                    int power = string.IsNullOrEmpty(powerStr) ? 0 : int.Parse(powerStr);
                    
                    coefficients[power] += coefficient;
                }
                
                return coefficients;
            }
            catch
            {
                return null;
            }
        }
    }

    /// <summary>
    /// 输入类型枚举
    /// </summary>
    public enum InputType
    {
        NaturalLanguage,
        MathematicalExpression
    }

    /// <summary>
    /// 复杂度级别
    /// </summary>
    public enum Complexity
    {
        Simple,
        Linear,
        Quadratic,
        Polynomial,
        Complex
    }

    /// <summary>
    /// 方程类别枚举
    /// </summary>
    public enum EquationCategory
    {
        General,
        Polynomial,
        Trigonometric,
        ExponentialLogarithmic
    }
}

// 为非线方程求解器添加扩展方法
public static class NonlinearSolverExtensions
{
    /// <summary>
    /// 为非线性求解器设置参数
    /// </summary>
    public static NonlinEqNLPSolver WithParameters(this NonlinEqNLPSolver solver, Dictionary<string, double> parameters)
    {
        foreach (var param in parameters)
        {
            solver.WithParameter(param.Key, param.Value);
        }
        return solver;
    }
}
