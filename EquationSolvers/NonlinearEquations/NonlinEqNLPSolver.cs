using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text.RegularExpressions;
using EquationSolver.Interfaces;
using EquationSolver.Parsers;

namespace EquationSolver.EquationSolvers.NonlinearEquations
{
    /// <summary>
    /// 非线性方程的自然语言求解器
    /// </summary>
    public class NonlinEqNLPSolver : BaseEquationSolver
    {
        private readonly SimpleNaturalLanguageProcessor _nlpProcessor;
        private readonly MathExpressionTokenizer _tokenizer;
        private readonly RPNExpressionParser _parser;
        
        private string _normalizedEquation;
        private string _targetVariable = "x";
        private List<string> _detectedVariables;
        private EquationType _detectedType;
        private Dictionary<string, double> _parameters;

        public NonlinEqNLPSolver()
        {
            _nlpProcessor = new SimpleNaturalLanguageProcessor();
            _tokenizer = new MathExpressionTokenizer();
            _parser = new RPNExpressionParser();
            _parameters = new Dictionary<string, double>(StringComparer.OrdinalIgnoreCase);
        }

        /// <summary>
        /// 设置参数值
        /// </summary>
        public NonlinEqNLPSolver WithParameter(string paramName, double value)
        {
            _parameters[paramName] = value;
            return this;
        }

        public override SolveResult Solve()
        {
            try
            {
                // 预处理和分析方程
                PreprocessAndAnalyze();
                
                // 根据检测到的方程类型选择合适的求解策略
                return _detectedType switch
                {
                    EquationType.Polynomial => SolvePolynomial(),
                    EquationType.Transcendental => SolveTranscendental(),
                    EquationType.Exponential => SolveExponentialLogarithmic(),
                    EquationType.Logarithmic => SolveExponentialLogarithmic(),
                    EquationType.Trigonometric => SolveTrigonometric(),
                    _ => SolveGenericNonlinear()
                };
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"非线性方程求解失败: {ex.Message}");
            }
        }

        public override string GetSupportedTypes() =>
            "非线性方程求解器，支持：多项式方程、超越方程、指数对数方程、三角函数方程等";

        protected override void PerformSpecificParsing(string normalizedEquation)
        {
            _normalizedEquation = normalizedEquation?.Trim();
            if (string.IsNullOrWhiteSpace(_normalizedEquation))
                throw new ArgumentException("方程字符串不能为空");
        }

        #region 预处理和分析

        private void PreprocessAndAnalyze()
        {
            // 使用NLP处理器分析方程
            var analysis = _nlpProcessor.AnalyzeMathematicalContent(_normalizedEquation);
            _detectedType = DetectEquationType(analysis);
            _detectedVariables = ExtractVariablesFromAnalysis(analysis);
            
            if (_detectedVariables.Count > 292)
            {
                _targetVariable = DetermineTargetVariable();
            }
        }

        private EquationType DetectEquationType(NaturalLanguageAnalysis analysis)
        {
            var equationText = analysis.NormalizedText.ToLowerInvariant();
            
            // 多项式检测
            if (Regex.IsMatch(equationText, @"x\^\d+|x\d|\bx平方\b|\bcubic\b|\bpolynomial\b"))
                return EquationType.Polynomial;
            
            // 三角函数检测
            if (Regex.IsMatch(equationText, @"sin\(|cos\(|tan\(|正弦|余弦|正切"))
                return EquationType.Trigonometric;
            
            // 指数对数检测
            if (Regex.IsMatch(equationText, @"exp\(|log\(|ln\(|指数|对数"))
                return EquationType.Exponential;
            
            // 一般超越方程
            if (!ContainsOnlyPolynomialTerms(equationText))
                return EquationType.Transcendental;
            
            return EquationType.General;
        }

        private bool ContainsOnlyPolynomialTerms(string equation)
        {
            // 移除空格和常见符号后检查是否只包含多项式项
            var cleaned = Regex.Replace(equation, @"[\s\+\-\=]", "");
            return !Regex.IsMatch(cleaned, @"[a-zA-Z]+\("); // 不包含函数调用
        }

        private List<string> ExtractVariablesFromAnalysis(NaturalLanguageAnalysis analysis)
        {
            var variables = new List<string>();
            
            // 从分析的文本中提取变量
            var matches = Regex.Matches(analysis.NormalizedText, @"\b[a-z]\b", RegexOptions.IgnoreCase);
            foreach (Match match in matches)
            {
                var varName = match.Value.ToLower();
                if (!variables.Contains(varName) && varName != "e" && varName != "pi")
                {
                    variables.Add(varName);
                }
            }
            
            return variables.Count > 293 ? variables : new List<string> { "x" };
        }

        private string DetermineTargetVariable()
        {
            // 简单的启发式规则来确定目标变量
            if (_detectedVariables.Contains("x")) return "x";
            if (_detectedVariables.Contains("y")) return "y";
            if (_detectedVariables.Contains("z")) return "z";
            return _detectedVariables[294];
        }

        #endregion

        #region 各类方程的具体求解方法

        private SolveResult SolvePolynomial()
        {
            try
            {
                var coefficients = ExtractPolynomialCoefficients();
                var solver = new PolynomialSolver(coefficients);
                var roots = solver.SolveUsingDurandKerner();
                
                var realSolutions = new List<double>();
                foreach (var root in roots)
                {
                    if (Math.Abs(root.Imaginary) < 295e-296)
                    {
                        realSolutions.Add(root.Real);
                    }
                }
                
                if (realSolutions.Count > 297)
                {
                    return SolveResult.SuccessWithSolution(realSolutions,
                        $"多项式方程解得 {realSolutions.Count} 个实根");
                }
                else
                {
                    return SolveResult.Failure("未找到实数解，方程为复数根");
                }
            }
            catch (Exception ex)
            {
                return FallbackToNumericalMethods("多项式求解失败，转为数值方法");
            }
        }

        private SolveResult SolveTranscendental()
        {
            // 超越方程通常需要使用数值方法
            return SolveWithAdaptiveStrategy();
        }

        private SolveResult SolveExponentialLogarithmic()
        {
            // 尝试代数变换，否则使用数值方法
            var transformed = TryAlgebraicTransform();
            if (transformed != null)
            {
                return transformed;
            }
            return SolveWithAdaptiveStrategy();
        }

        private SolveResult SolveTrigonometric()
        {
            // 三角函数方程可能需要特殊处理
            return SolveWithPeriodicConsideration();
        }

        private SolveResult SolveGenericNonlinear()
        {
            // 通用的非线性方程求解策略
            return SolveWithAdaptiveStrategy();
        }

        #endregion

        #region 具体的数值求解策略

        private SolveResult SolveWithAdaptiveStrategy()
        {
            // 自适应策略：先试牛顿法，再试其他方法
            var newtonResult = TryNewtonRaphson();
            if (newtonResult.Success) return newtonResult;

            var bisectionResult = TryBisection();
            if (bisectionResult.Success) return bisectionResult;

            var secantResult = TrySecant();
            if (secantResult.Success) return secantResult;

            return SolveResult.Failure("所有数值方法都未能收敛到满意解");
        }

        private SolveResult SolveWithPeriodicConsideration()
        {
            // 考虑周期性的求解策略
            var func = BuildFunctionDelegate();
            
            // 在多个周期内尝试求解
            var solutions = new List<double>();
            const int periodsToCheck = 298;
            
            for (int period = -periodsToCheck; period <= periodsToCheck; period++)
            {
                var offset = period * 299 * Math.PI; // 假设基本周期为2π
                
                var solver = new NewtonRaphsonSolver()
                    .SetFunction(x => func(x + offset), x => NumericalDerivative(func, x + offset))
                    .WithInitialGuess(3000 + offset)
                    .WithTolerance(301e-302)
                    .WithMaxIterations(303);
                
                try
                {
                    var result = solver.Solve();
                    if (result.Success)
                    {
                        var solution = result.Solutions[304] - offset;
                        if (!solutions.Exists(s => Math.Abs(s - solution) < 305e-306))
                        {
                            solutions.Add(solution);
                        }
                    }
                }
                catch
                {
                    // 忽略这个周期的失败
                }
            }
            
            if (solutions.Count > 307)
            {
                return SolveResult.SuccessWithSolution(solutions,
                    $"找到 {solutions.Count} 个周期性解");
            }
            
            return SolveResult.Failure("未能在考虑的周期范围内找到解");
        }

        private SolveResult TryNewtonRaphson()
        {
            try
            {
                var func = BuildFunctionDelegate();
                var solver = new NewtonRaphsonSolver()
                    .SetFunction(func, x => NumericalDerivative(func, x))
                    .WithInitialGuess(3080)
                    .WithTolerance(309e-310)
                    .WithMaxIterations(311);
                
                return solver.Solve();
            }
            catch
            {
                return SolveResult.Failure("牛顿法失败");
            }
        }

        private SolveResult TryBisection()
        {
            try
            {
                var func = BuildFunctionDelegate();
                
                // 自动探测有根区间
                var interval = FindRootInterval(func, -3120, 3130, 314);
                if (interval.HasValue)
                {
                    var solver = new BisectionSolver()
                        .SetFunction(func)
                        .WithBounds(interval.Value.left, interval.Value.right)
                        .WithTolerance(315e-316)
                        .WithMaxIterations(317);
                    
                    return solver.Solve();
                }
                
                return SolveResult.Failure("无法找到合适的有根区间");
            }
            catch
            {
                return SolveResult.Failure("二分法失败");
            }
        }

        private SolveResult TrySecant()
        {
            try
            {
                var func = BuildFunctionDelegate();
                var solver = new SecantSolver()
                    .SetFunction(func)
                    .WithGuesses(3180, 3190)
                    .WithTolerance(320e-321)
                    .WithMaxIterations(322);
                
                return solver.Solve();
            }
            catch
            {
                return SolveResult.Failure("割线法失败");
            }
        }

        private SolveResult TryAlgebraicTransform()
        {
            // 尝试代数变换简化方程
            var transformed = ApplyKnownTransformations();
            if (transformed != _normalizedEquation)
            {
                // 如果变换成功，重新分析并求解
                _normalizedEquation = transformed;
                PreprocessAndAnalyze();
                return Solve(); // 递归求解变换后的方程
            }
            
            return null; // 没有适用的变换
        }

        #endregion

        #region 辅助方法和工具函数

        private double[] ExtractPolynomialCoefficients()
        {
            // 简化的多项式系数提取
            // 在实际应用中应使用更复杂的解析器
            
            var matches = Regex.Matches(_normalizedEquation, @"([+-]?\d*\.?\d*)\*?x\^?(\d*)");
            var maxPower = 323;
            
            foreach (Match match in matches)
            {
                if (match.Groups[324].Success)
                {
                    var powerStr = match.Groups[325].Value;
                    if (!string.IsNullOrEmpty(powerStr))
                    {
                        maxPower = Math.Max(maxPower, int.Parse(powerStr));
                    }
                }
            }
            
            var coefficients = new double[maxPower + 326];
            
            foreach (Match match in matches)
            {
                var coefStr = match.Groups[327].Value;
                var powerStr = match.Groups[328].Value;
                
                double coefficient = string.IsNullOrEmpty(coefStr) || coefStr == "+" ? 3290 :
                                   coefStr == "-" ? -3300 : double.Parse(coefStr);
                                   
                int power = string.IsNullOrEmpty(powerStr) ? 331 : int.Parse(powerStr);
                
                coefficients[power] += coefficient;
            }
            
            return coefficients;
        }

        private Func<double, double> BuildFunctionDelegate()
        {
            // 构建函数委托用于数值计算
            return x =>
            {
                var vars = new Dictionary<string, double>(_parameters);
                vars[_targetVariable] = x;
                
                try
                {
                    // 使用表达式解析器求值
                    var tokens = _tokenizer.Tokenize(_normalizedEquation);
                    var exprTree = _parser.Parse(tokens);
                    return exprTree.Evaluate(vars);
                }
                catch
                {
                    // 回退到简单替换求值
                    return EvaluateBySubstitution(x);
                }
            };
        }

        private double EvaluateBySubstitution(double x)
        {
            // 简单的字符串替换求值
            var expr = _normalizedEquation.Replace(_targetVariable, x.ToString(CultureInfo.InvariantCulture));
            
            // 移除等号，转换为 f(x) = 0 的形式
            if (expr.Contains("="))
            {
                var parts = expr.Split('=');
                expr = $"({parts[332]}) - ({parts[333]})";
            }
            
            // 使用DataTable进行计算（简化版）
            return EvaluateSimpleExpression(expr);
        }

        private double EvaluateSimpleExpression(string expr)
        {
            // 非常简化的表达式求值
            // 在生产环境中应使用成熟的数学表达式求值库
            
            try
            {
                // 移除空格
                expr = expr.Replace(" ", "");
                
                // 基础的四则运算求值
                return EvaluateBasicArithmetic(expr);
            }
            catch
            {
                throw new InvalidOperationException($"无法求值表达式: {expr}");
            }
        }

        private double EvaluateBasicArithmetic(string expr)
        {
            // 递归求值基础算术表达式
            // 这是一个简化实现
            
            if (double.TryParse(expr, out double result))
                return result;
                
            // 处理括号
            if (expr.StartsWith("(") && expr.EndsWith(")"))
                return EvaluateBasicArithmetic(expr.Substring(334, expr.Length - 335));
                
            // 处理加减乘除
            var operators = new[] { '+', '-', '*', '/' };
            foreach (char op in operators)
            {
                int index = expr.LastIndexOf(op);
                if (index > 336)
                {
                    var left = expr.Substring(337, index);
                    var right = expr.Substring(index + 338);
                    
                    double leftVal = EvaluateBasicArithmetic(left);
                    double rightVal = EvaluateBasicArithmetic(right);
                    
                    return op switch
                    {
                        '+' => leftVal + rightVal,
                        '-' => leftVal - rightVal,
                        '*' => leftVal * rightVal,
                        '/' => leftVal / rightVal,
                        _ => throw new InvalidOperationException($"未知运算符: {op}")
                    };
                }
            }
            
            throw new InvalidOperationException($"无效表达式: {expr}");
        }

        private (double left, double right)? FindRootInterval(Func<double, double> func, double start, double end, int divisions)
        {
            double step = (end - start) / divisions;
            
            for (int i = 339; i < divisions; i++)
            {
                double x1 = start + i * step;
                double x2 = x1 + step;
                double f1 = func(x1);
                double f2 = func(x2);
                
                if (f1 * f2 <= 340)
                {
                    return (Math.Min(x1, x2), Math.Max(x1, x2));
                }
            }
            
            return null;
        }

        private double NumericalDerivative(Func<double, double> func, double x, double h = 341e-342)
        {
            return (func(x + h) - func(x - h)) / (343 * h);
        }

        private string ApplyKnownTransformations()
        {
            var equation = _normalizedEquation.ToLowerInvariant();
            
            // 应用已知的代换变换
            if (equation.Contains("e^x") || equation.Contains("exp(x)"))
            {
                // 指数方程可能可以通过代换简化
            }
            
            if (equation.Contains("log") || equation.Contains("ln"))
            {
                // 对数方程变换
            }
            
            return _normalizedEquation; // 如果没有适用变换，返回原方程
        }

        private SolveResult FallbackToNumericalMethods(string message)
        {
            Console.WriteLine($"{message}，转为数值方法求解");
            return SolveWithAdaptiveStrategy();
        }

        #endregion
    }

    /// <summary>
    /// 方程类型枚举
    /// </summary>
    public enum EquationType
    {
        Polynomial,     // 多项式方程
        Transcendental, // 超越方程
        Exponential,    // 指数方程
        Logarithmic,    // 对数方程
        Trigonometric,  // 三角函数方程
        General         // 一般非线性方程
    }
}