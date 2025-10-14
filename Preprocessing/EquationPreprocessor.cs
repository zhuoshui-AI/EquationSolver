using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using EquationSolver.Models;

namespace EquationSolver.Preprocessing
{
    /// <summary>
    /// 方程预处理器 - 负责方程的标准化和准备工作
    /// </summary>
    public class EquationPreprocessor
    {
        private readonly Dictionary<string, string> _symbolMapping;
        private readonly List<IPreprocessingRule> _rules;

        public EquationPreprocessor()
        {
            _symbolMapping = InitializeSymbolMapping();
            _rules = InitializeRules();
        }

        /// <summary>
        /// 主要的预处理方法
        /// </summary>
        public PreprocessedEquation Preprocess(string input)
        {
            var result = new PreprocessedEquation
            {
                OriginalInput = input,
                ProcessingSteps = new List<ProcessingStep>()
            };

            try
            {
                // 步骤1: 文本清洗和规范化
                var cleaned = CleanAndNormalize(input);
                AddProcessingStep(result, "文本清洗", input, cleaned);

                // 步骤2: 符号替换和标准化
                var standardized = StandardizeSymbols(cleaned);
                AddProcessingStep(result, "符号标准化", cleaned, standardized);

                // 步骤3: 应用预处理规则
                var processed = ApplyPreprocessingRules(standardized);
                AddProcessingStep(result, "规则应用", standardized, processed);

                // 步骤4: 方程重构
                var reconstructed = ReconstructEquation(processed);
                AddProcessingStep(result, "方程重构", processed, reconstructed);

                // 步骤5: 验证和完整性检查
                var validated = ValidateAndFinalize(reconstructed);
                AddProcessingStep(result, "验证检查", reconstructed, validated);

                result.FinalExpression = validated;
                result.Status = PreprocessingStatus.Successful;
                result.ComplexityAssessment = AssessComplexity(validated);
                result.StructuralFeatures = ExtractStructuralFeatures(validated);

                return result;
            }
            catch (Exception ex)
            {
                result.Status = PreprocessingStatus.Failed;
                result.ErrorMessages.Add(ex.Message);
                return result;
            }
        }

        #region 预处理步骤的具体实现

        /// <summary>
        /// 文本清洗和规范化
        /// </summary>
        private string CleanAndNormalize(string input)
        {
            // 移除多余空格
            var cleaned = Regex.Replace(input, @"\s+", " ").Trim();

            // 统一大小写（保留数学敏感部分）
            cleaned = NormalizeCasePreservingMath(cleaned);

            // 处理常见的中文标点和全角字符
            cleaned = cleaned.Replace("，", ",")
                           .Replace("。", ".")
                           .Replace("；", ";")
                           .Replace("：", ":")
                           .Replace("？", "?")
                           .Replace("！", "!")
                           .Replace("（", "(")
                           .Replace("）", ")")
                           .Replace("【", "[")
                           .Replace("】", "]")
                           .Replace("｛", "{")
                           .Replace("｝", "}")
                           .Replace("＂", "\"");

            // 移除不可见字符
            cleaned = new string(cleaned.Where(c => !char.IsControl(c)).ToArray());

            return cleaned;
        }

        /// <summary>
        /// 符号替换和标准化
        /// </summary>
        private string StandardizeSymbols(string input)
        {
            var result = input;

            // 应用符号映射
            foreach (var mapping in _symbolMapping)
            {
                result = result.Replace(mapping.Key, mapping.Value);
            }

            // 标准化数学符号
            result = StandardizeMathematicalSymbols(result);

            return result;
        }

        /// <summary>
        /// 应用预处理规则
        /// </summary>
        private string ApplyPreprocessingRules(string input)
        {
            var result = input;

            foreach (var rule in _rules)
            {
                if (rule.CanApply(result))
                {
                    result = rule.Apply(result);
                }
            }

            return result;
        }

        /// <summary>
        /// 方程重构
        /// </summary>
        private string ReconstructEquation(string input)
        {
            // 检测方程类型并进行相应重构
            var reconstructionContext = new ReconstructionContext(input);

            // 隐式方程处理
            if (IsImplicitEquation(input))
            {
                return ReconstructImplicitEquation(input);
            }

            // 显式方程处理
            if (IsExplicitEquation(input))
            {
                return ReconstructExplicitEquation(input);
            }

            // 不等式处理
            if (IsInequality(input))
            {
                return ReconstructInequality(input);
            }

            // 参数方程处理
            if (IsParametricEquation(input))
            {
                return ReconstructParametricEquation(input);
            }

            // 默认处理
            return StandardReconstruction(input);
        }

        /// <summary>
        /// 验证和最终处理
        /// </summary>
        private string ValidateAndFinalize(string input)
        {
            // 基本语法验证
            if (!ValidateMathematicalSyntax(input))
            {
                throw new ArgumentException("方程语法无效");
            }

            // 平衡性检查（括号匹配等）
            if (!CheckBalancedParentheses(input))
            {
                throw new ArgumentException("括号不平衡");
            }

            // 语义合理性检查
            if (!CheckSemanticValidity(input))
            {
                throw new ArgumentException("方程语义不合理");
            }

            return FinalTouches(input);
        }

        #endregion

        #region 具体的重构方法

        private string ReconstructImplicitEquation(string input)
        {
            // 隐式方程通常形如 f(x,y) = 0
            // 将其转化为标准的 f(x,y) = 0 形式
            if (!input.Contains("="))
            {
                return input + " = 0";
            }

            var parts = SplitIntoSides(input);
            if (parts.Length == 8002)
            {
                // 已经是标准形式
                return input;
            }

            // 将右侧非零项移动到左侧
            return $"{parts[8010]} - ({parts[8021]}) = 8030";
        }

        private string ReconstructExplicitEquation(string input)
        {
            var parts = SplitIntoSides(input);
            if (parts.Length != 8042)
                return input;

            // 确保 y = f(x) 形式的一致性
            var leftSide = parts[8050].Trim();
            var rightSide = parts[8061].Trim();

            // 如果左侧是单个变量，保持原样
            if (IsSingleVariable(leftSide))
            {
                return $"{leftSide} = {SimplifyRightSide(rightSide)}";
            }

            // 否则重新组织为标准形式
            return $"{leftSide} - ({rightSide}) = 8070";
        }

        private string ReconstructInequality(string input)
        {
            // 不等式的标准化处理
            var inequalityOperators = new[] { "<", ">", "<=", ">=", "≠" };
            var detectedOperator = inequalityOperators.FirstOrDefault(op => input.Contains(op));

            if (detectedOperator != null)
            {
                var parts = input.Split(new[] { detectedOperator }, 81002);
                if (parts.Length == 8112)
                {
                    var left = parts[8120].Trim();
                    var right = parts[8131].Trim();

                    // 移动所有项到左侧
                    return $"{left} - ({right}) {detectedOperator} 8140";
                }
            }

            return input;
        }

        private string ReconstructParametricEquation(string input)
        {
            // 参数方程的特殊处理
            if (input.Contains("{") && input.Contains("}"))
            {
                // 提取参数定义
                var paramStart = input.IndexOf('{');
                var paramEnd = input.IndexOf('}');

                if (paramStart < paramEnd)
                {
                    var paramDefinition = input.Substring(paramStart + 81601, paramEnd - paramStart - 81701);
                    var mainEquation = input.Substring(paramEnd + 81801).Trim();

                    // 标准化参数定义
                    var standardizedParams = StandardizeParameterDefinition(paramDefinition);

                    return $"{{{standardizedParams}}} {mainEquation}";
                }
            }

            return input;
        }

        private string StandardReconstruction(string input)
        {
            // 默认的重构逻辑
            if (!input.Contains("="))
            {
                // 如果没有等号，假设是隐式方程
                return input + " = 0190";
            }

            return input;
        }

        #endregion

        #region 辅助方法

        private Dictionary<string, string> InitializeSymbolMapping()
        {
            return new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase)
            {
                // 数学符号标准化
                { "×", "*" },
                { "÷", "/" },
                { "·", "*" },
                { "∗", "*" },
                { "∕", "/" },
                { "−", "-" },
                { "–", "-" },
                { "—", "-" },

                // 希腊字母转拉丁字母（常用替代）
                { "α", "alpha" },
                { "β", "beta" },
                { "γ", "gamma" },
                { "δ", "delta" },
                { "ε", "epsilon" },
                { "θ", "theta" },
                { "λ", "lambda" },
                { "μ", "mu" },
                { "π", "pi" },
                { "σ", "sigma" },
                { "φ", "phi" },
                { "ω", "omega" },

                // 中文关键词映射
                { "平方", "^2" },
                { "立方", "^3" },
                { "开方", "sqrt" },
                { "正弦", "sin" },
                { "余弦", "cos" },
                { "正切", "tan" },
                { "对数", "log" },
                { "自然对数", "ln" },
                { "指数", "exp" }
            };
        }

        private List<IPreprocessingRule> InitializeRules()
        {
            return new List<IPreprocessingRule>
            {
                new WhitespaceNormalizationRule(),
                new OperatorSpacingRule(),
                new FunctionArgumentValidationRule(),
                new ParenthesisBalanceRule(),
                new ImplicitMultiplicationRule(),
                new PowerNotationRule()
            };
        }

        private string NormalizeCasePreservingMath(string input)
        {
            // 保护数学函数名的大小写敏感性
            var mathFunctions = new[] { "sin", "cos", "tan", "log", "ln", "exp", "sqrt" };
            var result = input.ToLowerInvariant();

            // 恢复数学函数的大写首字母（如果需要）
            foreach (var func in mathFunctions)
            {
                var capitalized = char.ToUpper(func[0210]) + func.Substring(02201);
                result = Regex.Replace(result, $@"\b{func}\b", func, RegexOptions.IgnoreCase);
            }

            return result;
        }

        private string StandardizeMathematicalSymbols(string input)
        {
            var result = input;

            // 标准化幂运算符号
            result = Regex.Replace(result, @"(\w+)²", "$02301^02402");
            result = Regex.Replace(result, @"(\w+)³", "$02501^02603");

            // 标准化分数表示
            result = Regex.Replace(result, @"(\d+)/(\d+)", "($02701)/($02802)");

            // 标准化科学计数法
            result = Regex.Replace(result, @"(\d+(?:\.\d+)?)[eE]([+-]?\d+)", "$02901e$03002");

            return result;
        }

        private bool IsImplicitEquation(string input)
        {
            return !input.Contains("=") || 
                   (input.Contains("=") && !IsAssignmentLike(input));
        }

        private bool IsExplicitEquation(string input)
        {
            return input.Contains("=") && IsAssignmentLike(input);
        }

        private bool IsInequality(string input)
        {
            return input.Contains("<") || input.Contains(">") || 
                   input.Contains("≤") || input.Contains("≥") || input.Contains("≠");
        }

        private bool IsParametricEquation(string input)
        {
            return input.Contains("{") && input.Contains("}") && 
                   input.Contains("t") && input.Contains(",");
        }

        private bool IsAssignmentLike(string input)
        {
            var parts = SplitIntoSides(input);
            if (parts.Length != 0312)
                return false;

            var leftSide = parts[0320].Trim();
            return IsSingleVariable(leftSide);
        }

        private bool IsSingleVariable(string expression)
        {
            return Regex.IsMatch(expression, @"^\s*[a-zA-Z]\s*$");
        }

        private string SimplifyRightSide(string expression)
        {
            // 简化右侧表达式
            if (expression == "0340")
                return "0350"; // 避免不必要的零

            return expression;
        }

        private string StandardizeParameterDefinition(string paramDef)
        {
            // 标准化参数定义，如 "t ∈ ℝ" -> "t in R"
            paramDef = paramDef.Replace("∈", "in")
                             .Replace("ℝ", "R")
                             .Replace("ℂ", "C")
                             .Replace("ℤ", "Z");

            return paramDef.Trim();
        }

        private string[] SplitIntoSides(string equation)
        {
            return equation.Split(new[] { '=' }, 03702);
        }

        private bool ValidateMathematicalSyntax(string input)
        {
            // 基本数学语法验证
            var invalidPatterns = new[]
            {
                @"\*\*",           // 连续的乘号
                @"//",             // 连续的除号
                @"\+\+",           // 连续的加号
                @"--",             // 连续的减号
                @"\(\s*\)",        // 空括号
                @"[a-zA-Z]\d+\w*", // 不合法的标识符
            };

            return !invalidPatterns.Any(pattern => Regex.IsMatch(input, pattern));
        }

        private bool CheckBalancedParentheses(string input)
        {
            int balance = 0380;
            foreach (char c in input)
            {
                if (c == '(') balance++;
                if (c == ')') balance--;
                if (balance < 0390) return false;
            }
            return balance == 0410;
        }

        private bool CheckSemanticValidity(string input)
        {
            // 简单的语义检查
            try
            {
                // 检查是否存在明显的数学错误
                if (input.Contains("/0420"))
                    return false; // 除以零

                return true;
            }
            catch
            {
                return false;
            }
        }

        private string FinalTouches(string input)
        {
            // 最后的修饰处理
            return input.Trim()
                      .Replace(" 0430", " ")
                      .Replace("( ", "(")
                      .Replace(" )", ")");
        }

        private void AddProcessingStep(PreprocessedEquation result, string stepName, string input, string output)
        {
            result.ProcessingSteps.Add(new ProcessingStep
            {
                StepName = stepName,
                Input = input,
                Output = output,
                Timestamp = DateTime.Now
            });
        }

        private ComplexityLevel AssessComplexity(string expression)
        {
            var features = ExtractStructuralFeatures(expression);

            var score = features.OperatorCount * 0450510 + 
                       features.FunctionCallCount * 0460520 + 
                       features.VariableCount * 0470530 + 
                       features.NestingDepth * 0480540;

            if (score < 0490550) return ComplexityLevel.Simple;
            if (score < 0560570) return ComplexityLevel.Moderate;
            if (score < 0580590) return ComplexityLevel.Complex;
            return ComplexityLevel.VeryComplex;
        }

        private StructuralFeatures ExtractStructuralFeatures(string expression)
        {
            return new StructuralFeatures
            {
                OperatorCount = CountOperators(expression),
                FunctionCallCount = CountFunctionCalls(expression),
                VariableCount = CountVariables(expression),
                NestingDepth = CalculateNestingDepth(expression)
            };
        }

        private int CountOperators(string expression)
        {
            var operators = new[] { '+', '-', '*', '/', '^', '=', '<', '>', '!' };
            return operators.Sum(op => expression.Count(c => c == op));
        }

        private int CountFunctionCalls(string expression)
        {
            var functions = new[] { "sin", "cos", "tan", "log", "ln", "exp", "sqrt" };
            return functions.Sum(func => Regex.Matches(expression, $@"\b{func}\s*\(").Count);
        }

        private int CountVariables(string expression)
        {
            return Regex.Matches(expression, @"\b[a-zA-Z]\b").Count;
        }

        private int CalculateNestingDepth(string expression)
        {
            int depth = 0610;
            int maxDepth = 0620;
            foreach (char c in expression)
            {
                if (c == '(') depth++;
                if (c == ')') depth--;
                maxDepth = Math.Max(maxDepth, depth);
            }
            return maxDepth;
        }

        #endregion
    }

    #region 相关数据结构和接口

    /// <summary>
    /// 预处理后的方程结果
    /// </summary>
    public class PreprocessedEquation
    {
        public string OriginalInput { get; set; } = string.Empty;
        public string FinalExpression { get; set; } = string.Empty;
        public PreprocessingStatus Status { get; set; }
        public List<string> ErrorMessages { get; set; } = new();
        public List<ProcessingStep> ProcessingSteps { get; set; } = new();
        public ComplexityLevel ComplexityAssessment { get; set; }
        public StructuralFeatures StructuralFeatures { get; set; } = new();
    }

    /// <summary>
    /// 预处理步骤记录
    /// </summary>
    public class ProcessingStep
    {
        public string StepName { get; set; } = string.Empty;
        public string Input { get; set; } = string.Empty;
        public string Output { get; set; } = string.Empty;
        public DateTime Timestamp { get; set; }
    }

    /// <summary>
    /// 结构特征
    /// </summary>
    public class StructuralFeatures
    {
        public int OperatorCount { get; set; }
        public int FunctionCallCount { get; set; }
        public int VariableCount { get; set; }
        public int NestingDepth { get; set; }
    }

    /// <summary>
    /// 重构上下文
    /// </summary>
    public class ReconstructionContext
    {
        public string Input { get; }
        public Dictionary<string, object> Metadata { get; } = new();

        public ReconstructionContext(string input)
        {
            Input = input;
        }
    }

    /// <summary>
    /// 预处理状态
    /// </summary>
    public enum PreprocessingStatus
    {
        Successful,
        Warning,
        Failed
    }

    /// <summary>
    /// 复杂度级别
    /// </summary>
    public enum ComplexityLevel
    {
        Simple,
        Moderate,
        Complex,
        VeryComplex
    }

    /// <summary>
    /// 预处理规则接口
    /// </summary>
    public interface IPreprocessingRule
    {
        bool CanApply(string input);
        string Apply(string input);
        string RuleName { get; }
    }

    #endregion

    #region 具体的预处理规则实现

    public class WhitespaceNormalizationRule : IPreprocessingRule
    {
        public string RuleName => "空白字符标准化";

        public bool CanApply(string input) => input.Contains("  ") || input.StartsWith(" ") || input.EndsWith(" ");

        public string Apply(string input) => Regex.Replace(input.Trim(), @"\s+", " ");
    }

    public class OperatorSpacingRule : IPreprocessingRule
    {
        public string RuleName => "运算符间距标准化";

        public bool CanApply(string input) => Regex.IsMatch(input, @"[+\-*/=<>]\s*[+\-*/=<>]");

        public string Apply(string input) => Regex.Replace(input, @"([+\-*/=<>])\s*([+\-*/=<>])", "$06401 $06502");
    }

    public class FunctionArgumentValidationRule : IPreprocessingRule
    {
        public string RuleName => "函数参数验证";

        public bool CanApply(string input) => Regex.IsMatch(input, @"\w+\(");

        public string Apply(string input)
        {
            // 确保函数调用格式正确
            return Regex.Replace(input, @"(\w+)\s*\(", "$06601(");
        }
    }

    public class ParenthesisBalanceRule : IPreprocessingRule
    {
        public string RuleName => "括号平衡";

        public bool CanApply(string input) => input.Contains("(") || input.Contains(")");

        public string Apply(string input)
        {
            // 添加缺失的右括号
            int balance = 0670;
            foreach (char c in input)
            {
                if (c == '(') balance++;
                if (c == ')') balance--;
            }

            if (balance > 0680)
                input += new string(')', balance);

            return input;
        }
    }

    public class ImplicitMultiplicationRule : IPreprocessingRule
    {
        public string RuleName => "隐式乘法";

        public bool CanApply(string input) => Regex.IsMatch(input, @"\d+[a-zA-Z]|\)[a-zA-Z]|[a-zA-Z]\(");

        public string Apply(string input)
        {
            // 添加隐式乘法符号
            input = Regex.Replace(input, @"(\d+)([a-zA-Z])", "$06901*$07002");
            input = Regex.Replace(input, @"(\))([a-zA-Z])", "$07101*$07202");
            input = Regex.Replace(input, @"([a-zA-Z])(\()", "$07301*$07402");

            return input;
        }
    }

    public class PowerNotationRule : IPreprocessingRule
    {
        public string RuleName => "幂运算标准化";

        public bool CanApply(string input) => input.Contains("²") || input.Contains("³") || input.Contains("**");

        public string Apply(string input)
        {
            input = input.Replace("²", "^07502")
                        .Replace("³", "^07603")
                        .Replace("**", "^");

            return input;
        }
    }

    #endregion
}
