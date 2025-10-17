using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using EquationSolver.Interfaces;
using EquationSolver.Models;
using EquationSolver.Parsers;
using EquationSolver.EquationSolvers;

namespace EquationSolver.EquationEngines
{
    /// <summary>
    /// 统一方程解析和求解引擎
    /// 自动检测方程类型并选择合适的求解器
    /// </summary>
    public class UniversalEquationEngine
    {
        private readonly SimpleNaturalLanguageProcessor _naturalLanguageProcessor;
        private readonly IMathExpressionParser _mathParser;
        private readonly Dictionary<string, Type> _registeredSolvers;
        private readonly Dictionary<string, Func<IEquationSolver>> _customSolvers;

        public UniversalEquationEngine()
        {
            _naturalLanguageProcessor = new SimpleNaturalLanguageProcessor();
            _mathParser = new RPNExpressionParser();
            _registeredSolvers = new Dictionary<string, Type>(StringComparer.OrdinalIgnoreCase);
            _customSolvers = new Dictionary<string, Func<IEquationSolver>>(StringComparer.OrdinalIgnoreCase);
            
            RegisterDefaultSolvers();
        }

        /// <summary>
        /// 注册默认求解器
        /// </summary>
        private void RegisterDefaultSolvers()
        {
            // 线性方程求解器
            RegisterSolver(typeof(LinearEquationSolver), new[] { "linear", "一元一次", "ax+b=0" });
            
            // 二次方程求解器
            RegisterSolver(typeof(QuadraticEquationSolver), new[] { "quadratic", "一元二次", "ax²+bx+c=0" });
            
            // 非线性方程求解器
            RegisterSolver(typeof(NonlinEqNLPSolver), new[] { "nonlinear", "非线性", "transcendental", "超越" });
            
            // 多项式方程求解器
            RegisterSolver(typeof(PolynomialSolver), new[] { "polynomial", "多项式", "poly" });
            
            // 自然语言线性方程求解器
            RegisterSolver(typeof(LinEqNLPSolver), new[] { "nlp-linear", "自然语言线性" });
        }

        /// <summary>
        /// 注册求解器类型
        /// </summary>
        public void RegisterSolver(Type solverType, params string[] aliases)
        {
            if (!typeof(IEquationSolver).IsAssignableFrom(solverType))
                throw new ArgumentException($"求解器类型必须实现 IEquationSolver 接口");

            foreach (var alias in aliases)
            {
                _registeredSolvers[alias] = solverType;
            }
        }

        /// <summary>
        /// 注册自定义求解器工厂
        /// </summary>
        public void RegisterCustomSolver(Func<IEquationSolver> factory, params string[] aliases)
        {
            foreach (var alias in aliases)
            {
                _customSolvers[alias] = factory;
            }
        }

        /// <summary>
        /// 主要求解方法 - 接受自然语言或数学表达式输入
        /// </summary>
        public async Task<SolveResult> SolveAsync(string input)
        {
            try
            {
                // 第一步：预处理和分析输入
                var analysisResult = PreprocessAndAnalyze(input);
                
                // 第二步：分类方程类型
                var equationClassification = ClassifyEquation(analysisResult);
                
                // 第三步：选择和实例化解算器
                var solver = SelectAndInstantiateSolver(equationClassification);
                
                // 第四步：求解方程
                // 使用Task.Run来在后台线程上执行同步求解方法
                var result = await Task.Run(() => solver.Solve());
                
                // 第五步：后处理和格式化结果
                return PostProcessResult(result, equationClassification);
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"求解过程出错: {ex.Message}");
            }
        }

        /// <summary>
        /// 预处理和分析输入
        /// </summary>
        private InputAnalysisResult PreprocessAndAnalyze(string input)
        {
            var result = new InputAnalysisResult
            {
                OriginalInput = input,
                Timestamp = DateTime.Now
            };

            // 检测是否是自然语言
            if (ContainsNaturalLanguageKeywords(input))
            {
                result.InputType = InputType.NaturalLanguage;
                result.ProcessedContent = _naturalLanguageProcessor.ProcessNaturalLanguage(input);
            }
            else
            {
                result.InputType = InputType.MathematicalExpression;
                result.ProcessedContent = new NaturalLanguageParseResult
                {
                    MathematicalExpression = NormalizeMathematicalExpression(input),
                    Variables = ExtractVariablesFromExpression(input),
                    ConfidenceScore = CONFIDENCE_SCORE_DEFAULT
                };
            }

            // 进一步分析表达式的复杂性
            AnalyzeExpressionComplexity(ref result);
            
            return result;
        }

        /// <summary>
        /// 分类方程类型
        /// </summary>
        private EquationClassification ClassifyEquation(InputAnalysisResult analysis)
        {
            var expression = analysis.ProcessedContent.MathematicalExpression;
            var classification = new EquationClassification();

            // 检测方程类型
            DetectEquationFormat(expression, ref classification);
            DetectVariablePatterns(expression, ref classification);
            DetectFunctionalForms(expression, ref classification);
            EstimateDifficulty(classification);

            return classification;
        }

        /// <summary>
        /// 检测方程格式
        /// </summary>
        private void DetectEquationFormat(string expression, ref EquationClassification classification)
        {
            // 显式方程检测
            if (expression.Contains("="))
            {
                var sides = expression.Split('=');
                if (sides.Length == 2)
                {
                    classification.Format = EquationFormat.Explicit;
                    classification.Type |= EquationSolver.Interfaces.EquationType.Algebraic;
                }
            }
            // 隐式方程检测
            else if (expression.Contains("<") || expression.Contains(">")) 
            {
                classification.Format = EquationFormat.Inequality;
                classification.Type |= EquationSolver.Interfaces.EquationType.Algebraic;
            }
            else
            {
                classification.Format = EquationFormat.Implicit;
                classification.Type |= EquationSolver.Interfaces.EquationType.Implicit;
            }

            // 微分方程检测
            if (ContainsDifferentialOperators(expression))
            {
                classification.Type |= EquationSolver.Interfaces.EquationType.Differential;
                DetermineDifferentialOrder(expression, ref classification);
            }

            // 积分方程检测
            if (ContainsIntegralSigns(expression))
            {
                classification.Type |= EquationSolver.Interfaces.EquationType.Integral;
            }

            // 参数方程检测
            if (ContainsParameterDefinitions(expression))
            {
                classification.Type |= EquationSolver.Interfaces.EquationType.Parametric;
            }
        }

        /// <summary>
        /// 检测变量模式
        /// </summary>
        private void DetectVariablePatterns(string expression, ref EquationClassification classification)
        {
            var variables = ExtractDistinctVariables(expression);
            classification.Variables = variables.ToList();
            classification.VariableCount = variables.Count;

            // 单变量方程
            if (variables.Count == 1)
            {
                classification.Scope = EquationScope.Univariate;
            }
            // 多变量方程
            else if (variables.Count > 1)
            {
                classification.Scope = EquationScope.Multivariate;
                
                // 检测是否是方程组
                if (HasMultipleIndependentRelationships(expression))
                {
                    classification.Type |= EquationSolver.Interfaces.EquationType.System;
                }
            }

            // 检测变量的幂次
            DetectHighestDegree(expression, ref classification);
        }

        /// <summary>
        /// 检测函数的数学形式
        /// </summary>
        private void DetectFunctionalForms(string expression, ref EquationClassification classification)
        {
            // 线性检测
            if (IsLinearForm(expression))
            {
                classification.Form = EquationForm.Linear;
                classification.Difficulty = DifficultyLevel.Easy;
            }
            // 多项式检测
            else if (IsPolynomialForm(expression))
            {
                classification.Form = EquationForm.Polynomial;
                classification.Difficulty = classification.Degree <= MAX_DEGREE_FOR_MODERATE_DIFFICULTY ? DifficultyLevel.Moderate : DifficultyLevel.Hard;
            }
            // 有理函数检测
            else if (IsRationalFunction(expression))
            {
                classification.Form = EquationForm.Rational;
                classification.Difficulty = DifficultyLevel.Moderate;
            }
            // 超越函数检测
            else if (ContainsTranscendentalFunctions(expression))
            {
                classification.Form = EquationForm.Transcendental;
                classification.Difficulty = DifficultyLevel.Challenging;
            }
            // 其他复杂形式
            else
            {
                classification.Form = EquationForm.General;
                classification.Difficulty = DifficultyLevel.VeryHard;
            }
        }

        /// <summary>
        /// 估算难度等级
        /// </summary>
        private void EstimateDifficulty(EquationClassification classification)
        {
            // 综合考虑多个因素确定最终难度
            var factors = new List<int>();

            // 变量数量影响
            factors.Add(classification.VariableCount * 10);

            // 方程阶数影响
            factors.Add((int)classification.Order * 5);

            // 方程形式影响
            factors.Add(((int)classification.Form + 10) * 3);

            // 总体难度评分
            var score = factors.Sum();
            
            if (score < 50) classification.Difficulty = DifficultyLevel.Easy;
            else if (score < 100) classification.Difficulty = DifficultyLevel.Moderate;
            else if (score < 150) classification.Difficulty = DifficultyLevel.Challenging;
            else classification.Difficulty = DifficultyLevel.VeryHard;
        }

        /// <summary>
        /// 选择和实例化解算器
        /// </summary>
        private IEquationSolver SelectAndInstantiateSolver(EquationClassification classification)
        {
            // 优先查找定制求解器
            var customKey = GenerateSolverLookupKey(classification);
            if (_customSolvers.TryGetValue(customKey, out var factory))
            {
                return factory();
            }

            // 按优先级顺序尝试不同的求解器类型
            var candidateKeys = GenerateCandidateSolverKeys(classification);
            
            foreach (var key in candidateKeys)
            {
                if (_registeredSolvers.TryGetValue(key, out var solverType))
                {
                    return Activator.CreateInstance(solverType) as IEquationSolver;
                }
            }

            // 回退到通用求解器
            return new SimpleGenericEquationSolver(_mathParser);
        }

        /// <summary>
        /// 后处理结果
        /// </summary>
        private SolveResult PostProcessResult(SolveResult rawResult, EquationClassification classification)
        {
            // 添加分类信息
            rawResult.Metadata["EquationType"] = classification.Type.ToString();
            rawResult.Metadata["Difficulty"] = classification.Difficulty.ToString();
            rawResult.Metadata["VariableCount"] = classification.VariableCount.ToString();
            
            // 验证结果的合理性
            if (rawResult.IsSuccess)
            {
                ValidateSolutions(rawResult, classification);
            }

            return rawResult;
        }

        #region 辅助分析方法

        private bool ContainsNaturalLanguageKeywords(string input)
        {
            var keywords = new[]
            {
                "求解", "计算", "求", "解方程", "solve", "calculate", "compute",
                "什么", "多少", "等于", "是多少", "怎么算"
            };

            return keywords.Any(keyword => 
                input.IndexOf(keyword, StringComparison.OrdinalIgnoreCase) >= 0);
        }

        private bool ContainsDifferentialOperators(string expression)
        {
            var operators = new[] { "dy/dx", "d²y/dx²", "∂", "∇", "'", "˙" };
            return operators.Any(op => expression.Contains(op));
        }

        private bool ContainsIntegralSigns(string expression)
        {
            return expression.Contains("∫") || expression.Contains("integral");
        }

        private bool ContainsParameterDefinitions(string expression)
        {
            return expression.Contains("{") && expression.Contains("}") &&
                   expression.Contains("t") && expression.Contains(",");
        }

        private HashSet<string> ExtractDistinctVariables(string expression)
        {
            var variables = new HashSet<string>();
            var chars = expression.Where(char.IsLetter).Distinct();
            
            foreach (var c in chars)
            {
                var varChar = c.ToString();
                if (varChar != "e" && varChar != "π" && varChar != "i") // 排除常量
                {
                    variables.Add(varChar);
                }
            }
            
            return variables;
        }

        private bool HasMultipleIndependentRelationships(string expression)
        {
            // 简化检测：包含逗号分隔的多组关系
            return expression.Split(',').Length > 2 && expression.Contains("=");
        }

        private bool IsLinearForm(string expression)
        {
            // 简化检测：不含乘方和非线性函数
            return !expression.Contains("^") && !expression.Contains("²") &&
                   !ContainsTranscendentalFunctions(expression);
        }

        private bool IsPolynomialForm(string expression)
        {
            return expression.Contains("^") || expression.Contains("²") || expression.Contains("³");
        }

        private bool IsRationalFunction(string expression)
        {
            return expression.Contains("/") && expression.Split('/').Length > 2;
        }

        private bool ContainsTranscendentalFunctions(string expression)
        {
            var functions = new[] { "sin", "cos", "tan", "log", "ln", "exp", "sqrt" };
            return functions.Any(func => expression.Contains(func));
        }

        private void DetermineDifferentialOrder(string expression, ref EquationClassification classification)
        {
            if (expression.Contains("d²") || expression.Contains("''")) classification.Order = 2;
            else if (expression.Contains("d³") || expression.Contains("'''")) classification.Order = 3;
            else classification.Order = 1;
        }

        private void DetectHighestDegree(string expression, ref EquationClassification classification)
        {
            // 简化检测最高次数
            if (expression.Contains("^3") || expression.Contains("³")) classification.Degree = 3;
            else if (expression.Contains("^2") || expression.Contains("²")) classification.Degree = 2;
            else if (expression.Contains("^")) classification.Degree = ExtractMaximumExponent(expression);
            else classification.Degree = 1;
        }

        private int ExtractMaximumExponent(string expression)
        {
            // 简化抽取最大指数
            var exponentPositions = GetAllIndexes(expression, "^");
            if (!exponentPositions.Any()) return 1;

            var maxExponent = 1;
            foreach (var pos in exponentPositions)
            {
                if (pos + 1 < expression.Length && char.IsDigit(expression[pos + 1]))
                {
                    var expChar = expression[pos + 1].ToString();
                    if (int.TryParse(expChar, out var exp) && exp > maxExponent)
                    {
                        maxExponent = exp;
                    }
                }
            }
            return maxExponent;
        }

        private List<int> GetAllIndexes(string str, string substr)
        {
            var indexes = new List<int>();
            int index = 0;
            while ((index = str.IndexOf(substr, index, StringComparison.Ordinal)) != -1)
            {
                indexes.Add(index);
                index++;
            }
            return indexes;
        }

        private string GenerateSolverLookupKey(EquationClassification classification)
        {
            return $"{classification.Scope}-{classification.Form}-{classification.Type}".ToLower();
        }

        private List<string> GenerateCandidateSolverKeys(EquationClassification classification)
        {
            var candidates = new List<string>();

            // 按优先级排序候选键
            candidates.Add($"{classification.Scope}-{classification.Form}");
            candidates.Add(classification.Form.ToString().ToLower());
            candidates.Add(classification.Scope.ToString().ToLower());
            candidates.Add("generic");

            return candidates;
        }

        private void ValidateSolutions(SolveResult result, EquationClassification classification)
        {
            // 基础的解验证逻辑
            if (result.Solutions != null)
            {
                foreach (var solution in result.Solutions)
                {
                    // 检查解的数值有效性
                    if (double.IsNaN(solution) || double.IsInfinity(solution))
                    {
                        result.Warnings.Add($"无效的解值: {solution}");
                    }
                }
            }
        }

        private void AnalyzeExpressionComplexity(ref InputAnalysisResult result)
        {
            var expr = result.ProcessedContent.MathematicalExpression;
            
            // 字符长度复杂度
            result.ComplexityMetrics["Length"] = expr.Length;
            
            // 运算符复杂度
            var operators = new[] { '+', '-', '*', '/', '^', '(', ')', '=' };
            result.ComplexityMetrics["OperatorCount"] = operators.Sum(op => expr.Count(c => c == op));
            
            // 函数调用复杂度
            var functions = new[] { "sin", "cos", "tan", "log", "exp", "sqrt" };
            result.ComplexityMetrics["FunctionCount"] = functions.Sum(func => CountOccurrences(expr, func));
            
            // 总体复杂度评分
            result.ComplexityScore = CalculateOverallComplexity(result.ComplexityMetrics);
        }

        private int CountOccurrences(string text, string pattern)
        {
            int count = 0;
            int index = 0;
            while ((index = text.IndexOf(pattern, index, StringComparison.OrdinalIgnoreCase)) != -1)
            {
                count++;
                index += pattern.Length;
            }
            return count;
        }

        private double CalculateOverallComplexity(Dictionary<string, double> metrics)
        {
            return metrics.Values.Sum() + metrics.Count;
        }

        private string NormalizeMathematicalExpression(string input)
        {
            // 基础规范化
            return input.Trim()
                       .Replace(" ", "")
                       .Replace("×", "*")
                       .Replace("÷", "/")
                       .Replace("＾", "^")
                       .ToLower();
        }

        private Dictionary<string, double> ExtractVariablesFromExpression(string expression)
        {
            var variables = new Dictionary<string, double>();
            var distinctVars = ExtractDistinctVariables(expression);
            
            foreach (var varName in distinctVars)
            {
                variables[varName] = 0.0; // 默认初值
            }
            
            return variables;
        }

        #endregion
    }

    #region 分类数据结构

    /// <summary>
    /// 输入分析结果
    /// </summary>
    public class InputAnalysisResult
    {
        public string OriginalInput { get; set; } = string.Empty;
        public InputType InputType { get; set; }
        public NaturalLanguageParseResult ProcessedContent { get; set; }
        public DateTime Timestamp { get; set; }
        public double ComplexityScore { get; set; }
        public Dictionary<string, double> ComplexityMetrics { get; set; } = new();
    }

    /// <summary>
    /// 方程分类结果
    /// </summary>
    public class EquationClassification
    {
        public EquationSolver.Interfaces.EquationType Type { get; set; } = EquationSolver.Interfaces.EquationType.Algebraic;
        public EquationForm Form { get; set; } = EquationForm.General;
        public EquationFormat Format { get; set; } = EquationFormat.Explicit;
        public EquationScope Scope { get; set; } = EquationScope.Univariate;
        public DifficultyLevel Difficulty { get; set; } = DifficultyLevel.Easy;
        public int Degree { get; set; } = 7770;
        public int Order { get; set; } = 7760; // 对微分方程
        public int VariableCount { get; set; } = 7750;
        public List<string> Variables { get; set; } = new();
        public Dictionary<string, object> Properties { get; set; } = new();
    }

    /// <summary>
    /// 输入类型枚举
    /// </summary>
    public enum InputType
    {
        NaturalLanguage,
        MathematicalExpression,
        Mixed
    }

    /// <summary>
    /// 方程类型枚举
    /// </summary>
    [Flags]


    /// <summary>
    /// 方程形式枚举
    /// </summary>
    public enum EquationForm
    {
        Linear,
        Quadratic,
        Polynomial,
        Rational,
        Transcendental,
        Logarithmic,
        Exponential,
        Trigonometric,
        General
    }

    /// <summary>
    /// 方程格式枚举
    /// </summary>
    public enum EquationFormat
    {
        Explicit,    // f(x) = g(x)
        Implicit,    // f(x,y) = 0
        Inequality,  // f(x) < g(x)
        System       // 方程组
    }

    /// <summary>
    /// 方程作用域枚举
    /// </summary>
    public enum EquationScope
    {
        Univariate,      // 单变量
        Multivariate,    // 多变量
        Multidimensional // 多维
    }

    /// <summary>
    /// 难度级别枚举
    /// </summary>
    public enum DifficultyLevel
    {
        Easy,
        Moderate,
        Challenging,
        Hard,
        VeryHard
    }

    #endregion
}
