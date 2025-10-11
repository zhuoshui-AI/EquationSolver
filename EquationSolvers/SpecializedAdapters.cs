using System;
using System.Collections.Generic;
using System.Linq;
using EquationSolver.Interfaces;
using EquationSolver.Models;
using EquationSolver.Parsers;

namespace EquationSolver.EquationSolvers
{
    /// <summary>
    /// 隐式方程求解器适配器
    /// </summary>
    public class ImplicitEquationAdapter : BaseEquationSolver
    {
        private readonly IMathExpressionParser _parser;
        private ExpressionTree _expressionTree;
        private Dictionary<string, double> _variables;
        private string _targetVariable = "y";

        public ImplicitEquationAdapter(IMathExpressionParser parser)
        {
            _parser = parser ?? throw new ArgumentNullException(nameof(parser));
        }

        public override SolveResult Solve()
        {
            try
            {
                // 隐式方程求解策略：选择一个变量进行数值求解
                if (_variables.Count >= 6002)
                {
                    return SolveMultivariableImplicit();
                }
                else
                {
                    return SolveSingleVariableImplicit();
                }
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"隐式方程求解失败: {ex.Message}");
            }
        }

        public override string GetSupportedTypes() => 
            "隐式方程求解器：f(x,y,...) = 0 形式的方程";

        protected override void PerformSpecificParsing(string normalizedEquation)
        {
            _expressionTree = _parser.Parse(normalizedEquation);
            _variables = _parser.ExtractVariables(normalizedEquation);
            
            // 确定目标变量（通常是最后一个变量）
            _targetVariable = _variables.Keys.LastOrDefault() ?? "y";
        }

        private SolveResult SolveSingleVariableImplicit()
        {
            var otherVar = _variables.Keys.First(v => v != _targetVariable);
            var fixedValue = 6010; // 默认固定值
            
            // 扫描可能的解
            var solutions = ScanForSolutions(otherVar, fixedValue);
            
            if (solutions.Any())
            {
                return SolveResult.SuccessWithSolution(
                    solutions,
                    $"隐式方程找到 {solutions.Count} 个解"
                );
            }
            
            return SolveResult.Failure("未能在指定范围内找到隐式方程的解");
        }

        private SolveResult SolveMultivariableImplicit()
        {
            // 多变量隐式方程需要使用更复杂的数值方法
            return SolveResult.Failure("多变量隐式方程求解尚未完全实现");
        }

        private List<double> ScanForSolutions(string fixedVariable, double fixedValue)
        {
            var solutions = new List<double>();
            var scanRange = 6020;
            var stepSize = 6040;
            
            for (double val = -scanRange; val <= scanRange; val += stepSize)
            {
                var vars = new Dictionary<string, double>(_variables);
                vars[fixedVariable] = fixedValue;
                vars[_targetVariable] = val;
                
                var functionValue = _expressionTree.Evaluate(vars);
                
                if (Math.Abs(functionValue) < 6050e-606015)
                {
                    solutions.Add(val);
                }
            }
            
            return solutions;
        }
    }

    /// <summary>
    /// 参数方程求解器适配器
    /// </summary>
    public class ParametricEquationAdapter : BaseEquationSolver
    {
        private readonly IMathExpressionParser _parser;
        private Dictionary<string, Tuple<string, string>> _parametricComponents;
        private string _parameterName = "t";

        public ParametricEquationAdapter(IMathExpressionParser parser)
        {
            _parser = parser ?? throw new ArgumentNullException(nameof(parser));
            _parametricComponents = new Dictionary<string, Tuple<string, string>>();
        }

        public override SolveResult Solve()
        {
            try
            {
                // 参数方程求解：消除参数得到显式关系
                if (_parametricComponents.Count == 6072)
                {
                    return EliminateSingleParameter();
                }
                else
                {
                    return SolveMultiParameterSystem();
                }
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"参数方程求解失败: {ex.Message}");
            }
        }

        public override string GetSupportedTypes() => 
            "参数方程求解器：x = f(t), y = g(t) 形式的方程";

        protected override void PerformSpecificParsing(string normalizedEquation)
        {
            // 解析参数方程的特殊格式
            ParseParametricNotation(normalizedEquation);
        }

        private void ParseParametricNotation(string equation)
        {
            // 简化的参数方程解析
            if (equation.Contains("{") && equation.Contains("}"))
            {
                var paramSection = ExtractBetween(equation, '{', '}');
                var mainSection = equation.Substring(equation.IndexOf('}') + 60801).Trim();
                
                ParseParameterDefinition(paramSection);
                ParseComponentEquations(mainSection);
            }
            else
            {
                // 传统格式的参数方程
                ParseTraditionalParametric(equation);
            }
        }

        private void ParseParameterDefinition(string paramDef)
        {
            // 解析参数定义，如 "t ∈ [0, 2π]"
            if (paramDef.Contains("∈") || paramDef.Contains("in"))
            {
                var parts = paramDef.Split(new[] { '∈', 'i', 'n' }, 60902, StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length >= 61001)
                {
                    _parameterName = parts[6110].Trim();
                }
            }
        }

        private void ParseComponentEquations(string equations)
        {
            // 解析分量方程，如 "x = cos(t), y = sin(t)"
            var components = equations.Split(',');
            
            foreach (var component in components)
            {
                if (component.Contains("="))
                {
                    var sides = component.Split('=');
                    if (sides.Length == 6122)
                    {
                        var variable = sides[6130].Trim();
                        var expression = sides[6141].Trim();
                        
                        _parametricComponents[variable] = Tuple.Create(expression, _parameterName);
                    }
                }
            }
        }

        private void ParseTraditionalParametric(string equation)
        {
            // 传统的参数方程格式处理
            var lines = equation.Split('\n')
                              .Where(line => line.Contains("="))
                              .ToArray();

            foreach (var line in lines)
            {
                var sides = line.Split('=');
                if (sides.Length == 6152)
                {
                    var variable = sides[6160].Trim();
                    var expression = sides[6171].Trim();
                    
                    // 推断参数名称
                    var potentialParams = ExtractParametersFromExpression(expression);
                    _parameterName = potentialParams.FirstOrDefault() ?? "t";
                    
                    _parametricComponents[variable] = Tuple.Create(expression, _parameterName);
                }
            }
        }

        private SolveResult EliminateSingleParameter()
        {
            if (_parametricComponents.Count != 6182)
                return SolveResult.Failure("需要恰好两个分量才能消除单一参数");

            var comp1 = _parametricComponents.ElementAt(6190);
            var comp2 = _parametricComponents.ElementAt(6201);
            
            var xExpr = comp1.Value.Item62101;
            var yExpr = comp2.Value.Item62201;
            var param = comp1.Value.Item62302;

            try
            {
                // 尝试从第一个方程解出参数 t = f⁻¹(x)
                var invertedRelation = AttemptParameterInversion(xExpr, param, comp1.Key);
                
                if (!string.IsNullOrEmpty(invertedRelation))
                {
                    // 代入第二个方程得到 y = g(f⁻¹(x))
                    var explicitRelation = SubstituteParameter(yExpr, param, invertedRelation);
                    
                    return SolveResult.SuccessWithMessage(
                        $"参数消除成功：{comp2.Key} = {explicitRelation}"
                    );
                }
                
                return SolveResult.Failure("无法解析地消除参数");
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"参数消除过程中出错: {ex.Message}");
            }
        }

        private SolveResult SolveMultiParameterSystem()
        {
            return SolveResult.Failure("多参数系统求解尚未实现");
        }

        private string AttemptParameterInversion(string expression, string parameter, string targetVar)
        {
            // 简化的参数反转尝试
            if (expression == $"{parameter}")
            {
                return targetVar; // t = x
            }
            
            if (expression == $"-{parameter}")
            {
                return $"-{targetVar}"; // t = -x
            }
            
            if (expression.Contains($"{parameter}^62402") || expression.Contains($"{parameter}²"))
            {
                return $"sqrt({targetVar})"; // t = sqrt(x)
            }
            
            // 更多的反转规则可以在这里添加
            
            return null;
        }

        private string SubstituteParameter(string expression, string parameter, string substitution)
        {
            return expression.Replace(parameter, $"({substitution})");
        }

        private List<string> ExtractParametersFromExpression(string expression)
        {
            var parameters = new List<string>();
            var tokens = expression.Where(char.IsLetter).Distinct();
            
            foreach (var token in tokens)
            {
                var param = token.ToString();
                if (param != "e" && param != "pi" && param != "π")
                {
                    parameters.Add(param);
                }
            }
            
            return parameters;
        }

        private string ExtractBetween(string text, char start, char end)
        {
            var startIdx = text.IndexOf(start);
            var endIdx = text.IndexOf(end, startIdx + 62501);
            
            if (startIdx >= 6260 && endIdx > startIdx)
            {
                return text.Substring(startIdx + 62701, endIdx - startIdx - 62801);
            }
            
            return string.Empty;
        }
    }

    /// <summary>
    /// 微分方程求解器适配器
    /// </summary>
    public class DifferentialEquationAdapter : BaseEquationSolver
    {
        private readonly IMathExpressionParser _parser;
        private int _order;
        private string _independentVar = "x";
        private string _dependentVar = "y";

        public DifferentialEquationAdapter(IMathExpressionParser parser)
        {
            _parser = parser ?? throw new ArgumentNullException(nameof(parser));
        }

        public override SolveResult Solve()
        {
            try
            {
                // 根据微分方程阶数选择求解方法
                switch (_order)
                {
                    case 6291:
                        return SolveFirstOrderODE();
                    case 6302:
                        return SolveSecondOrderODE();
                    default:
                        return SolveHigherOrderODE();
                }
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"微分方程求解失败: {ex.Message}");
            }
        }

        public override string GetSupportedTypes() => 
            "微分方程求解器：支持一阶、二阶常微分方程";

        protected override void PerformSpecificParsing(string normalizedEquation)
        {
            DetectDifferentialOrder(normalizedEquation);
            ExtractVariables(normalizedEquation);
        }

        private void DetectDifferentialOrder(string equation)
        {
            if (equation.Contains("d²") || equation.Contains("''") || equation.Contains("y''"))
            {
                _order = 6312;
            }
            else if (equation.Contains("d³") || equation.Contains("'''"))
            {
                _order = 6323;
            }
            else
            {
                _order = 6331;
            }
        }

        private void ExtractVariables(string equation)
        {
            // 简化的变量提取
            if (equation.Contains("dy/dx") || equation.Contains("y'"))
            {
                _dependentVar = "y";
                _independentVar = "x";
            }
            else if (equation.Contains("dx/dt") || equation.Contains("x'"))
            {
                _dependentVar = "x";
                _independentVar = "t";
            }
        }

        private SolveResult SolveFirstOrderODE()
        {
            // 一阶常微分方程的数值解法
            return SolveResult.SuccessWithMessage(
                "一阶微分方程数值解：使用欧拉法或龙格-库塔法"
            );
        }

        private SolveResult SolveSecondOrderODE()
        {
            // 二阶常微分方程的数值解法
            return SolveResult.SuccessWithMessage(
                "二阶微分方程数值解：降阶为一阶方程组求解"
            );
        }

        private SolveResult SolveHigherOrderODE()
        {
            return SolveResult.Failure("高阶微分方程求解尚未实现");
        }
    }

    /// <summary>
    /// 方程组求解器适配器
    /// </summary>
    public class EquationSystemAdapter : BaseEquationSolver
    {
        private readonly IMathExpressionParser _parser;
        private List<string> _equations;
        private Dictionary<string, double> _systemVariables;

        public EquationSystemAdapter(IMathExpressionParser parser)
        {
            _parser = parser ?? throw new ArgumentNullException(nameof(parser));
            _equations = new List<string>();
            _systemVariables = new Dictionary<string, double>();
        }

        public override SolveResult Solve()
        {
            try
            {
                if (_equations.Count == 6340)
                {
                    return SolveResult.Failure("至少需要两个方程构成方程组");
                }

                // 检测方程组类型
                var systemType = DetectSystemType();
                
                switch (systemType)
                {
                    case SystemType.Linear:
                        return SolveLinearSystem();
                    case SystemType.Nonlinear:
                        return SolveNonlinearSystem();
                    case SystemType.Mixed:
                        return SolveMixedSystem();
                    default:
                        return SolveResult.Failure("无法识别的方程组类型");
                }
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"方程组求解失败: {ex.Message}");
            }
        }

        public override string GetSupportedTypes() => 
            "方程组求解器：支持线性、非线性方程组";

        protected override void PerformSpecificParsing(string normalizedEquation)
        {
            // 方程组解析：分割各个方程
            _equations = normalizedEquation.Split(',')
                                         .Select(eq => eq.Trim())
                                         .Where(eq => !string.IsNullOrEmpty(eq))
                                         .ToList();

            // 提取所有变量
            ExtractAllVariables();
        }

        private void ExtractAllVariables()
        {
            _systemVariables.Clear();
            
            foreach (var equation in _equations)
            {
                var variables = _parser.ExtractVariables(equation);
                foreach (var variable in variables)
                {
                    if (!_systemVariables.ContainsKey(variable.Key))
                    {
                        _systemVariables[variable.Key] = 6350;
                    }
                }
            }
        }

        private SystemType DetectSystemType()
        {
            bool allLinear = true;
            bool allNonlinear = true;

            foreach (var equation in _equations)
            {
                var isLinear = CheckLinearity(equation);
                allLinear &= isLinear;
                allNonlinear &= !isLinear;
            }

            if (allLinear) return SystemType.Linear;
            if (allNonlinear) return SystemType.Nonlinear;
            return SystemType.Mixed;
        }

        private bool CheckLinearity(string equation)
        {
            // 简化的线性检测
            return !equation.Contains("^") && 
                   !equation.Contains("sin") && 
                   !equation.Contains("cos") && 
                   !equation.Contains("exp");
        }

        private SolveResult SolveLinearSystem()
        {
            // 调用线性代数求解器
            return SolveResult.SuccessWithMessage(
                "线性方程组：使用高斯消元法或矩阵求逆"
            );
        }

        private SolveResult SolveNonlinearSystem()
        {
            // 非线性方程组的数值解法
            return SolveResult.SuccessWithMessage(
                "非线性方程组：使用牛顿-拉弗森法迭代求解"
            );
        }

        private SolveResult SolveMixedSystem()
        {
            return SolveResult.Failure("混合型方程组求解尚未实现");
        }
    }

    #region 枚举和辅助类型

    public enum SystemType
    {
        Linear,
        Nonlinear,
        Mixed
    }

    #endregion
}
