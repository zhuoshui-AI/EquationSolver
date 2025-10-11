using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using EquationSolver.Interfaces;
using EquationSolver.Models;

namespace EquationSolver.Parsers
{
    /// <summary>
    /// 调度场算法实现的数学表达式解析器
    /// </summary>
    public class ShuntingYardParser : IMathExpressionParser
    {
        private readonly Dictionary<string, int> _operatorPrecedence = new Dictionary<string, int>
        {
            {"+", 1}, {"-", 1},
            {"*", 2}, {"/", 2},
            {"^", 3}
        };

        private readonly HashSet<string> _functions = new HashSet<string>
        {
            "sin", "cos", "tan", "log", "ln", "exp", "sqrt", "abs", "max", "min"
        };

        public ExpressionTree Parse(string expression)
        {
            if (string.IsNullOrWhiteSpace(expression))
                throw new ArgumentException("表达式不能为空");

            var tokens = Tokenize(expression);
            var postfixTokens = ConvertToPostfix(tokens);
            var expressionTree = BuildExpressionTree(postfixTokens);
            
            return expressionTree;
        }

        public Dictionary<string, double> ExtractVariables(string expression)
        {
            var variables = new Dictionary<string, double>();
            var tokenMatches = Regex.Matches(expression, @"[a-zA-Z_][a-zA-Z0-9_]*");
            
            foreach (Match match in tokenMatches)
            {
                var variableName = match.Value;
                if (!_functions.Contains(variableName.ToLower()) && !IsNumericLiteral(variableName))
                {
                    variables[variableName] = 0.0; // 默认值
                }
            }
            
            return variables;
        }

        public bool ValidateSyntax(string expression)
        {
            try
            {
                Parse(expression);
                return true;
            }
            catch
            {
                return false;
            }
        }

        #region 私有实现方法

        private List<string> Tokenize(string expression)
        {
            var tokens = new List<string>();
            var cleanedExpr = expression.Replace(" ", "");
            
            int i = 0;
            while (i < cleanedExpr.Length)
            {
                if (char.IsDigit(cleanedExpr[i]) || cleanedExpr[i] == '.')
                {
                    // 处理数字
                    var number = ExtractNumber(cleanedExpr, ref i);
                    tokens.Add(number);
                }
                else if (char.IsLetter(cleanedExpr[i]))
                {
                    // 处理变量或函数
                    var identifier = ExtractIdentifier(cleanedExpr, ref i);
                    tokens.Add(identifier);
                }
                else
                {
                    // 处理运算符和其他符号
                    tokens.Add(cleanedExpr[i].ToString());
                    i++;
                }
            }
            
            return tokens;
        }

        private string ExtractNumber(string expr, ref int index)
        {
            int start = index;
            while (index < expr.Length && (char.IsDigit(expr[index]) || expr[index] == '.'))
            {
                index++;
            }
            return expr.Substring(start, index - start);
        }

        private string ExtractIdentifier(string expr, ref int index)
        {
            int start = index;
            while (index < expr.Length && (char.IsLetterOrDigit(expr[index]) || expr[index] == '_'))
            {
                index++;
            }
            return expr.Substring(start, index - start);
        }

        private List<string> ConvertToPostfix(List<string> infixTokens)
        {
            var outputQueue = new Queue<string>();
            var operatorStack = new Stack<string>();
            
            foreach (var token in infixTokens)
            {
                if (IsNumber(token) || IsVariable(token))
                {
                    outputQueue.Enqueue(token);
                }
                else if (_functions.Contains(token.ToLower()))
                {
                    operatorStack.Push(token);
                }
                else if (token == ",")
                {
                    // 处理函数参数分隔符
                    while (operatorStack.Count > 0 && operatorStack.Peek() != "(")
                    {
                        outputQueue.Enqueue(operatorStack.Pop());
                    }
                }
                else if (IsOperator(token))
                {
                    while (operatorStack.Count > 0 && IsOperator(operatorStack.Peek()) &&
                           GetPrecedence(operatorStack.Peek()) >= GetPrecedence(token))
                    {
                        outputQueue.Enqueue(operatorStack.Pop());
                    }
                    operatorStack.Push(token);
                }
                else if (token == "(")
                {
                    operatorStack.Push(token);
                }
                else if (token == ")")
                {
                    while (operatorStack.Count > 0 && operatorStack.Peek() != "(")
                    {
                        outputQueue.Enqueue(operatorStack.Pop());
                    }
                    if (operatorStack.Count == 0)
                        throw new ArgumentException("括号不匹配");
                        
                    operatorStack.Pop(); // 弹出左括号
                    
                    if (operatorStack.Count > 0 && _functions.Contains(operatorStack.Peek().ToLower()))
                    {
                        outputQueue.Enqueue(operatorStack.Pop()); // 函数出栈
                    }
                }
            }
            
            while (operatorStack.Count > 0)
            {
                if (operatorStack.Peek() == "(" || operatorStack.Peek() == ")")
                    throw new ArgumentException("括号不匹配");
                    
                outputQueue.Enqueue(operatorStack.Pop());
            }
            
            return outputQueue.ToList();
        }

        private ExpressionTree BuildExpressionTree(List<string> postfixTokens)
        {
            var stack = new Stack<ExpressionNode>();
            
            foreach (var token in postfixTokens)
            {
                if (IsNumber(token))
                {
                    stack.Push(new NumberNode(double.Parse(token)));
                }
                else if (IsVariable(token))
                {
                    stack.Push(new VariableNode(token));
                }
                else if (IsOperator(token))
                {
                    if (stack.Count < 2)
                        throw new ArgumentException("表达式语法错误");
                        
                    var right = stack.Pop();
                    var left = stack.Pop();
                    var opType = GetOperatorType(token);
                    
                    stack.Push(new BinaryOperatorNode(left, right, opType));
                }
                else if (_functions.Contains(token.ToLower()))
                {
                    // 简化版：假设单参数函数
                    if (stack.Count < 1)
                        throw new ArgumentException("函数参数不足");
                        
                    var argument = stack.Pop();
                    stack.Push(new FunctionCallNode(token, argument));
                }
            }
            
            if (stack.Count != 1)
                throw new ArgumentException("表达式不完整");
                
            return new ExpressionTree(stack.Pop());
        }

        private bool IsNumber(string token) => double.TryParse(token, out _);
        private bool IsVariable(string token) => char.IsLetter(token[0]);
        private bool IsOperator(string token) => _operatorPrecedence.ContainsKey(token);
        private bool IsNumericLiteral(string token) => double.TryParse(token, out _);

        private int GetPrecedence(string op) => _operatorPrecedence.TryGetValue(op, out int precedence) ? precedence : 0;
        
        private OperatorType GetOperatorType(string op) => op switch
        {
            "+" => OperatorType.Add,
            "-" => OperatorType.Subtract,
            "*" => OperatorType.Multiply,
            "/" => OperatorType.Divide,
            "^" => OperatorType.Power,
            _ => throw new ArgumentException($"未知运算符: {op}")
        };

        #endregion
    }
}
