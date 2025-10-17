using System;

namespace EquationSolver.Parsers
{
    /// <summary>
    /// 数学表达式标记类型枚举
    /// </summary>
    public enum TokenType
    {
        Number,         // 数字
        Variable,       // 变量
        Operator,       // 运算符
        Function,       // 函数
        LeftParenthesis,// 左括号
        RightParenthesis,// 右括号
        Comma,          // 逗号（函数参数分隔符）
        Constant,       // 数学常数
        EndOfExpression // 表达式结束
    }

    /// <summary>
    /// 运算符结合性
    /// </summary>
    public enum Associativity
    {
        Left,
        Right
    }

    /// <summary>
    /// 数学表达式标记类
    /// </summary>
    public class Token
    {
        public TokenType Type { get; set; }
        public string Value { get; set; }
        public int Position { get; set; }
        
        // 仅适用于运算符
        public int Precedence { get; set; }
        public Associativity Associtivity { get; set; }

        public Token(TokenType type, string value, int position)
        {
            Type = type;
            Value = value;
            Position = position;
            
            // 设置默认优先级和结合性
            SetOperatorProperties(value);
        }

        private void SetOperatorProperties(string opValue)
        {
            switch (opValue)
            {
                case "+":
                case "-":
                    Precedence = 2;
                    Associtivity = Associativity.Left;
                    break;
                case "*":
                case "/":
                    Precedence = 3;
                    Associtivity = Associativity.Left;
                    break;
                case "^":
                    Precedence = 4;
                    Associtivity = Associativity.Right;
                    break;
                default:
                    Precedence = 0;
                    Associtivity = Associativity.Left;
                    break;
            }
        }

        public override string ToString()
        {
            return $"Token({Type}: '{Value}' at pos {Position})";
        }

        public bool IsBinaryOperator()
        {
            return Type == TokenType.Operator && 
                   (Value == "+" || Value == "-" || Value == "*" || Value == "/" || Value == "^");
        }

        public bool IsUnaryOperator()
        {
            return Type == TokenType.Operator && 
                   (Value == "+" || Value == "-"); // 一元加减号
        }

        public bool IsRightAssociative()
        {
            return Associtivity == Associativity.Right;
        }

        public bool HasHigherPrecedenceThan(Token other)
        {
            return Precedence > other.Precedence;
        }

        public bool HasSamePrecedenceAs(Token other)
        {
            return Precedence == other.Precedence;
        }
    }

    /// <summary>
    /// 支持的数学函数列表
    /// </summary>
    public static class SupportedFunctions
    {
        public static readonly string[] Functions = 
        {
            "sin", "cos", "tan", "cot", "sec", "csc",
            "asin", "acos", "atan", "acot", "asec", "acsc",
            "sinh", "cosh", "tanh", "coth", "sech", "csch",
            "asinh", "acosh", "atanh", "acoth", "asech", "acsch",
            "log", "lg", "ln", "exp", "sqrt", "cbrt",
            "abs", "floor", "ceil", "round", "sign",
            "factorial", "gamma", "beta", "erf", "zeta"
        };

        public static bool IsFunction(string name)
        {
            return Array.Exists(Functions, func => func.Equals(name, StringComparison.OrdinalIgnoreCase));
        }

        public static int GetArgumentCount(string functionName)
        {
            switch (functionName.ToLower())
            {
                case "sin": case "cos": case "tan": case "cot": case "sec": case "csc":
                case "asin": case "acos": case "atan": case "acot": case "asec": case "acsc":
                case "sinh": case "cosh": case "tanh": case "coth": case "sech": case "csch":
                case "asinh": case "acosh": case "atanh": case "acoth": case "asech": case "acsch":
                case "log": case "lg": case "ln": case "exp": case "sqrt": case "cbrt":
                case "abs": case "floor": case "ceil": case "round": case "sign":
                case "factorial": case "gamma": case "erf": case "zeta":
                    return 1;
                case "pow": case "mod": case "hypot": case "gcd": case "lcm":
                case "beta": case "polyval":
                    return 2;
                case "if": case "choose": case "piecewise":
                    return -1; // 可变参数
                default:
                    return 0;
            }
        }
    }

    /// <summary>
    /// 数学常数定义
    /// </summary>
    public static class MathConstants
    {
        public static readonly (string name, double value)[] Constants =
        {
            ("pi", Math.PI),
            ("π", Math.PI),
            ("e", Math.E),
            ("infinity", double.PositiveInfinity),
            ("∞", double.PositiveInfinity),
            ("nan", double.NaN),
            ("phi", 1.618033988749895), // 黄金比例
            ("γ", 0.577215664901532),   // Euler-Mascheroni常数
            ("deg", Math.PI / 180.0), // 角度转弧度因子
            ("grad", Math.PI / 200.0) // 梯度转弧度因子
        };

        public static bool IsConstant(string name)
        {
            return Array.Exists(Constants, constant => 
                constant.name.Equals(name, StringComparison.OrdinalIgnoreCase));
        }

        public static double GetConstantValue(string name)
        {
            var constant = Array.Find(Constants, c => 
                c.name.Equals(name, StringComparison.OrdinalIgnoreCase));
            return constant.value;
        }
    }
}