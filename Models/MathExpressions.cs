using System;
using System.Collections.Generic;
using System.Linq;

namespace EquationSolver.Models
{
    /// <summary>
    /// 数学表达式树节点
    /// </summary>
    public abstract class ExpressionNode
    {
        public abstract double Evaluate(Dictionary<string, double> variables);
        public abstract ExpressionNode Simplify();
        public abstract IEnumerable<string> GetVariables();
        public abstract override string ToString();
    }

    /// <summary>
    /// 数值常量节点
    /// </summary>
    public class NumberNode : ExpressionNode
    {
        public double Value { get; }

        public NumberNode(double value) => Value = value;

        public override double Evaluate(Dictionary<string, double> variables) => Value;
        
        public override ExpressionNode Simplify() => this;
        
        public override IEnumerable<string> GetVariables() => Enumerable.Empty<string>();
        
        public override string ToString() => Value.ToString();
    }

    /// <summary>
    /// 变量节点
    /// </summary>
    public class VariableNode : ExpressionNode
    {
        public string Name { get; }

        public VariableNode(string name) => Name = name ?? throw new ArgumentNullException(nameof(name));

        public override double Evaluate(Dictionary<string, double> variables)
            => variables.TryGetValue(Name, out var value) ? value : throw new KeyNotFoundException($"变量 '{Name}' 未定义");

        public override ExpressionNode Simplify() => this;
        
        public override IEnumerable<string> GetVariables() => new[] { Name };
        
        public override string ToString() => Name;
    }

    /// <summary>
    /// 二元运算节点
    /// </summary>
    public class BinaryOperatorNode : ExpressionNode
    {
        public ExpressionNode Left { get; }
        public ExpressionNode Right { get; }
        public OperatorType Operator { get; }

        public BinaryOperatorNode(ExpressionNode left, ExpressionNode right, OperatorType op)
        {
            Left = left ?? throw new ArgumentNullException(nameof(left));
            Right = right ?? throw new ArgumentNullException(nameof(right));
            Operator = op;
        }

        public override double Evaluate(Dictionary<string, double> variables)
        {
            var leftVal = Left.Evaluate(variables);
            var rightVal = Right.Evaluate(variables);

            return Operator switch
            {
                OperatorType.Add => leftVal + rightVal,
                OperatorType.Subtract => leftVal - rightVal,
                OperatorType.Multiply => leftVal * rightVal,
                OperatorType.Divide => rightVal != 0 ? leftVal / rightVal : throw new DivideByZeroException(),
                OperatorType.Power => Math.Pow(leftVal, rightVal),
                _ => throw new InvalidOperationException($"未知运算符: {Operator}")
            };
        }

        public override ExpressionNode Simplify()
        {
            var simplifiedLeft = Left.Simplify();
            var simplifiedRight = Right.Simplify();

            // 简化规则实现
            return ApplySimplificationRules(simplifiedLeft, simplifiedRight);
        }

        private ExpressionNode ApplySimplificationRules(ExpressionNode left, ExpressionNode right)
        {
            // 这里可以实现各种代数简化规则
            return new BinaryOperatorNode(left, right, Operator);
        }

        public override IEnumerable<string> GetVariables() =>
            Left.GetVariables().Concat(Right.GetVariables()).Distinct();

        public override string ToString() => $"({Left} {GetOperatorSymbol()} {Right})";

        private string GetOperatorSymbol() => Operator switch
        {
            OperatorType.Add => "+",
            OperatorType.Subtract => "-",
            OperatorType.Multiply => "*",
            OperatorType.Divide => "/",
            OperatorType.Power => "^",
            _ => "?"
        };
    }

    /// <summary>
    /// 函数调用节点
    /// </summary>
    public class FunctionCallNode : ExpressionNode
    {
        public string FunctionName { get; }
        public List<ExpressionNode> Arguments { get; }

        public FunctionCallNode(string functionName, params ExpressionNode[] arguments)
        {
            FunctionName = functionName ?? throw new ArgumentNullException(nameof(functionName));
            Arguments = arguments?.ToList() ?? new List<ExpressionNode>();
        }

        public override double Evaluate(Dictionary<string, double> variables)
        {
            var argValues = Arguments.Select(arg => arg.Evaluate(variables)).ToArray();

            return FunctionName.ToLowerInvariant() switch
            {
                "sin" => Math.Sin(argValues[0]),
                "cos" => Math.Cos(argValues[0]),
                "tan" => Math.Tan(argValues[0]),
                "log" => Math.Log(argValues[0]),
                "ln" => Math.Log(argValues[0]),
                "exp" => Math.Exp(argValues[0]),
                "sqrt" => Math.Sqrt(argValues[0]),
                "abs" => Math.Abs(argValues[0]),
                "max" => Math.Max(argValues[0], argValues.Length > 1 ? argValues[1] : double.MinValue),
                "min" => Math.Min(argValues[0], argValues.Length > 1 ? argValues[1] : double.MinValue),
                _ => throw new NotImplementedException($"不支持函数: {FunctionName}")
            };
        }

        public override ExpressionNode Simplify()
        {
            var simplifiedArgs = Arguments.Select(arg => arg.Simplify()).ToList();
            return new FunctionCallNode(FunctionName, simplifiedArgs.ToArray());
        }

        public override IEnumerable<string> GetVariables() =>
            Arguments.SelectMany(arg => arg.GetVariables()).Distinct();

        public override string ToString() => $"{FunctionName}({string.Join(", ", Arguments)})";
    }

    /// <summary>
    /// 运算符类型枚举
    /// </summary>
    public enum OperatorType
    {
        Add,
        Subtract,
        Multiply,
        Divide,
        Power
    }

    /// <summary>
    /// 表达式树封装类
    /// </summary>
    public class ExpressionTree
    {
        public ExpressionNode Root { get; }
        public HashSet<string> Variables => new HashSet<string>(Root.GetVariables());

        public ExpressionTree(ExpressionNode root) => Root = root ?? throw new ArgumentNullException(nameof(root));

        public double Evaluate(Dictionary<string, double> variableValues) => Root.Evaluate(variableValues);
        
        public ExpressionTree Simplify() => new ExpressionTree(Root.Simplify());
        
        public override string ToString() => Root.ToString();
    }

    /// <summary>
    /// 方程求解结果
    /// </summary>
    public class SolveResult
    {
        public bool Success { get; set; }
        public List<double> Solutions { get; set; } = new List<double>();
        public string Message { get; set; } = "";
        public int Iterations { get; set; }
        public double Residual { get; set; }
        public TimeSpan ElapsedTime { get; set; }
        public ConvergenceStatus Convergence { get; set; }

        public static SolveResult Failure(string message) => new SolveResult 
        { 
            Success = false, 
            Message = message 
        };

        public static SolveResult SuccessWithSolution(List<double> solutions, string message = "") => new SolveResult
        {
            Success = true,
            Solutions = solutions,
            Message = message
        };
    }

    /// <summary>
    /// 收敛状态
    /// </summary>
    public enum ConvergenceStatus
    {
        Converged,
        Diverged,
        MaximumIterationsReached,
        ToleranceAchieved,
        NumericalInstabilityDetected
    }
}
