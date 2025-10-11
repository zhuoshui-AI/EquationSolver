using System;
using System.Collections.Generic;

namespace EquationSolver.Parsers
{
    /// <summary>
    /// 语法树节点基类
    /// </summary>
    public abstract class SyntaxTreeNode
    {
        public abstract double Evaluate(Dictionary<string, double> variables);
        
        /// <summary>
        /// 返回节点的字符串表示
        /// </summary>
        public abstract override string ToString();
        
        /// <summary>
        /// 递归复制节点及其子节点
        /// </summary>
        public abstract SyntaxTreeNode Clone();
    }

    /// <summary>
    /// 数字常量节点
    /// </summary>
    public class ConstantNode : SyntaxTreeNode
    {
        public double Value { get; set; }

        public ConstantNode(double value)
        {
            Value = value;
        }

        public override double Evaluate(Dictionary<string, double> variables)
        {
            return Value;
        }

        public override string ToString()
        {
            return Value.ToString(System.Globalization.CultureInfo.InvariantCulture);
        }

        public override SyntaxTreeNode Clone()
        {
            return new ConstantNode(Value);
        }
    }

    /// <summary>
    /// 变量节点
    /// </summary>
    public class VariableNode : SyntaxTreeNode
    {
        public string Name { get; set; }

        public VariableNode(string name)
        {
            Name = name ?? throw new ArgumentNullException(nameof(name));
        }

        public override double Evaluate(Dictionary<string, double> variables)
        {
            if (variables == null || !variables.ContainsKey(Name))
                throw new KeyNotFoundException($"未定义的变量: {Name}");
            
            return variables[Name];
        }

        public override string ToString()
        {
            return Name;
        }

        public override SyntaxTreeNode Clone()
        {
            return new VariableNode(Name);
        }
    }

    /// <summary>
    /// 二元运算符节点
    /// </summary>
    public class BinaryOperatorNode : SyntaxTreeNode
    {
        public string Operator { get; set; }
        public SyntaxTreeNode Left { get; set; }
        public SyntaxTreeNode Right { get; set; }

        public BinaryOperatorNode(string @operator, SyntaxTreeNode left, SyntaxTreeNode right)
        {
            Operator = @operator ?? throw new ArgumentNullException(nameof(@operator));
            Left = left ?? throw new ArgumentNullException(nameof(left));
            Right = right ?? throw new ArgumentNullException(nameof(right));
        }

        public override double Evaluate(Dictionary<string, double> variables)
        {
            var leftVal = Left.Evaluate(variables);
            var rightVal = Right.Evaluate(variables);

            switch (Operator)
            {
                case "+": return leftVal + rightVal;
                case "-": return leftVal - rightVal;
                case "*": return leftVal * rightVal;
                case "/":
                    if (Math.Abs(rightVal) < double.Epsilon)
                        throw new DivideByZeroException("除以零错误");
                    return leftVal / rightVal;
                case "^": return Math.Pow(leftVal, rightVal);
                case "%": return leftVal % rightVal;
                default:
                    throw new NotSupportedException($"不支持的运算符: {Operator}");
            }
        }

        public override string ToString()
        {
            return $"({Left} {Operator} {Right})";
        }

        public override SyntaxTreeNode Clone()
        {
            return new BinaryOperatorNode(Operator, Left.Clone(), Right.Clone());
        }
    }

    /// <summary>
    /// 一元运算符节点
    /// </summary>
    public class UnaryOperatorNode : SyntaxTreeNode
    {
        public string Operator { get; set; }
        public SyntaxTreeNode Operand { get; set; }

        public UnaryOperatorNode(string @operator, SyntaxTreeNode operand)
        {
            Operator = @operator ?? throw new ArgumentNullException(nameof(@operator));
            Operand = operand ?? throw new ArgumentNullException(nameof(operand));
        }

        public override double Evaluate(Dictionary<string, double> variables)
        {
            var operandVal = Operand.Evaluate(variables);

            switch (Operator)
            {
                case "+": return operandVal;
                case "-": return -operandVal;
                case "!":
                    if (operandVal < 0)
                        throw new ArgumentOutOfRangeException("阶乘不能应用于负数");
                    return Factorial((int)operandVal);
                default:
                    throw new NotSupportedException($"不支持的一元运算符: {Operator}");
            }
        }

        private double Factorial(int n)
        {
            if (n <= 1)
                return 1;
            
            double result = 1;
            for (int i = 2; i <= n; i++)
                result *= i;
            return result;
        }

        public override string ToString()
        {
            return $"{Operator}({Operand})";
        }

        public override SyntaxTreeNode Clone()
        {
            return new UnaryOperatorNode(Operator, Operand.Clone());
        }
    }

    /// <summary>
    /// 函数调用节点
    /// </summary>
    public class FunctionCallNode : SyntaxTreeNode
    {
        public string FunctionName { get; set; }
        public List<SyntaxTreeNode> Arguments { get; set; }

        public FunctionCallNode(string functionName, params SyntaxTreeNode[] arguments)
        {
            FunctionName = functionName ?? throw new ArgumentNullException(nameof(functionName));
            Arguments = new List<SyntaxTreeNode>(arguments ?? Array.Empty<SyntaxTreeNode>());
        }

        public override double Evaluate(Dictionary<string, double> variables)
        {
            var argValues = new List<double>();
            foreach (var arg in Arguments)
            {
                argValues.Add(arg.Evaluate(variables));
            }

            return CallFunction(FunctionName.ToLowerInvariant(), argValues.ToArray());
        }

        private double CallFunction(string functionName, double[] args)
        {
            switch (functionName)
            {
                case "sin": return Math.Sin(args[0]);
                case "cos": return Math.Cos(args[0]);
                case "tan": return Math.Tan(args[0]);
                case "asin": return Math.Asin(args[0]);
                case "acos": return Math.Acos(args[0]);
                case "atan": return Math.Atan(args[0]);
                case "log": return Math.Log(args[0]); // 自然对数
                case "ln": return Math.Log(args[0]); // 自然对数别名
                case "lg": return Math.Log10(args[0]); // 常用对数
                case "exp": return Math.Exp(args[0]);
                case "sqrt": return Math.Sqrt(args[0]);
                case "pow": return Math.Pow(args[0], args[1]);
                case "max": return Math.Max(args[0], args[1]);
                case "min": return Math.Min(args[0], args[1]);
                case "abs": return Math.Abs(args[0]);
                case "floor": return Math.Floor(args[0]);
                case "ceil": return Math.Ceiling(args[0]);
                case "round": return Math.Round(args[0]);

                // 自定义函数
                case "factorial":
                    if (args[0] < 0)
                        throw new ArgumentOutOfRangeException("阶乘不能应用于负数");
                    return Factorial((int)args[0]);
                
                default:
                    throw new NotSupportedException($"未知的函数: {functionName}");
            }
        }

        private double Factorial(int n)
        {
            if (n <= 1)
                return 1;
            
            double result = 1;
            for (int i = 2; i <= n; i++)
                result *= i;
            return result;
        }

        public override string ToString()
        {
            var args = string.Join(", ", Arguments.Select(a => a.ToString()));
            return $"{FunctionName}({args})";
        }

        public override SyntaxTreeNode Clone()
        {
            var clonedArgs = Arguments.Select(arg => arg.Clone()).ToArray();
            return new FunctionCallNode(FunctionName, clonedArgs);
        }
    }
}

// 为了简化LINQ查询，我们需要添加必要的using语句
public static partial class Extensions
{
    public static IEnumerable<T> SelectMany<T>(this T[] items, Func<T, IEnumerable<T>> selector)
    {
        foreach (T item in items)
        {
            foreach (T result in selector(item))
            {
                yield return result;
            }
        }
    }
}