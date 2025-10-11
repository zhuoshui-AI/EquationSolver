using System;
using System.Collections.Generic;
using System.Text.RegularExpressions;

namespace EquationSolver.Parsers
{
    /// <summary>
    /// 数学表达式分词器
    /// </summary>
    public class MathExpressionTokenizer
    {
        private readonly string _expression;
        private int _position;
        private readonly List<Token> _tokens;

        // 正则表达式模式
        private static readonly Regex NumberPattern = new Regex(@"^-?\d+(\.\d+)?([eE][+-]?\d+)?", RegexOptions.Compiled);
        private static readonly Regex IdentifierPattern = new Regex(@"^[a-zA-Z_][a-zA-Z0-9_]*", RegexOptions.Compiled);
        private static readonly Regex OperatorPattern = new Regex(@"^[\+\-\*/%\^!&=<>~\|@#\$\\]", RegexOptions.Compiled);
        private static readonly Regex WhitespacePattern = new Regex(@"^\s+", RegexOptions.Compiled);

        public MathExpressionTokenizer(string expression)
        {
            _expression = expression?.Trim() ?? throw new ArgumentNullException(nameof(expression));
            _position = 0;
            _tokens = new List<Token>();
        }

        /// <summary>
        /// 将表达式分解为标记序列
        /// </summary>
        public List<Token> Tokenize()
        {
            _tokens.Clear();
            _position = 0;

            while (_position < _expression.Length)
            {
                SkipWhitespace();

                if (_position >= _expression.Length)
                    break;

                char currentChar = _expression[_position];

                if (TryParseNumber(out Token numberToken))
                {
                    _tokens.Add(numberToken);
                }
                else if (TryParseIdentifierOrKeyword(out Token identifierToken))
                {
                    _tokens.Add(identifierToken);
                }
                else if (TryParseOperator(out Token operatorToken))
                {
                    _tokens.Add(operatorToken);
                }
                else if (TryParseParentheses(currentChar, out Token parenToken))
                {
                    _tokens.Add(parenToken);
                }
                else if (TryParseComma(currentChar, out Token commaToken))
                {
                    _tokens.Add(commaToken);
                }
                else
                {
                    throw new InvalidOperationException($"无法识别的字符 '{currentChar}' 在位置 {_position}");
                }
            }

            // 添加结束标记
            _tokens.Add(new Token(TokenType.EndOfExpression, "", _position));

            ResolveAmbiguity(); // 解决歧义（如一元和二元运算符）
            return _tokens;
        }

        /// <summary>
        /// 跳过空白字符
        /// </summary>
        private void SkipWhitespace()
        {
            var match = WhitespacePattern.Match(_expression.Substring(_position));
            if (match.Success)
            {
                _position += match.Length;
            }
        }

        /// <summary>
        /// 解析数字字面量
        /// </summary>
        private bool TryParseNumber(out Token token)
        {
            token = null;
            var substring = _expression.Substring(_position);
            var match = NumberPattern.Match(substring);

            if (match.Success)
            {
                token = new Token(TokenType.Number, match.Value, _position);
                _position += match.Length;
                return true;
            }

            return false;
        }

        /// <summary>
        /// 解析标识符或关键字（变量、函数、常数）
        /// </summary>
        private bool TryParseIdentifierOrKeyword(out Token token)
        {
            token = null;
            var substring = _expression.Substring(_position);
            var match = IdentifierPattern.Match(substring);

            if (match.Success)
            {
                string identifier = match.Value;
                int tokenStartPos = _position;
                _position += match.Length;

                // 检查是否为函数
                if (IsFollowedByOpeningParenthesis() && SupportedFunctions.IsFunction(identifier))
                {
                    token = new Token(TokenType.Function, identifier, tokenStartPos);
                }
                // 检查是否为数学常数
                else if (MathConstants.IsConstant(identifier))
                {
                    token = new Token(TokenType.Constant, identifier, tokenStartPos);
                }
                // 否则视为变量
                else
                {
                    token = new Token(TokenType.Variable, identifier, tokenStartPos);
                }

                return true;
            }

            return false;
        }

        /// <summary>
        /// 检查标识符后是否跟着左括号（判断是否为函数调用）
        /// </summary>
        private bool IsFollowedByOpeningParenthesis()
        {
            // 跳过空白字符
            int tempPos = _position;
            while (tempPos < _expression.Length && char.IsWhiteSpace(_expression[tempPos]))
            {
                tempPos++;
            }

            return tempPos < _expression.Length && _expression[tempPos] == '(';
        }

        /// <summary>
        /// 解析运算符
        /// </summary>
        private bool TryParseOperator(out Token token)
        {
            token = null;
            var substring = _expression.Substring(_position);
            var match = OperatorPattern.Match(substring);

            if (match.Success)
            {
                string op = match.Value;
                token = new Token(TokenType.Operator, op, _position);
                _position += match.Length;
                return true;
            }

            return false;
        }

        /// <summary>
        /// 解析括号
        /// </summary>
        private bool TryParseParentheses(char currentChar, out Token token)
        {
            token = null;

            if (currentChar == '(')
            {
                token = new Token(TokenType.LeftParenthesis, "(", _position);
                _position++;
                return true;
            }
            else if (currentChar == ')')
            {
                token = new Token(TokenType.RightParenthesis, ")", _position);
                _position++;
                return true;
            }

            return false;
        }

        /// <summary>
        /// 解析逗号（函数参数分隔符）
        /// </summary>
        private bool TryParseComma(char currentChar, out Token token)
        {
            token = null;

            if (currentChar == ',')
            {
                token = new Token(TokenType.Comma, ",", _position);
                _position++;
                return true;
            }

            return false;
        }

        /// <summary>
        /// 解决运算符歧义（区分一元和二元运算符）
        /// </summary>
        private void ResolveAmbiguity()
        {
            for (int i = 0; i < _tokens.Count; i++)
            {
                var currentToken = _tokens[i];
                
                if (currentToken.Type == TokenType.Operator && 
                    (currentToken.Value == "+" || currentToken.Value == "-"))
                {
                    // 判断是一元还是二元运算符
                    bool isUnary = IsUnaryContext(i);
                    if (isUnary)
                    {
                        // 对于一元运算符，我们可以特殊处理或者保持原样
                        // 在实际解析时会根据上下文进行处理
                    }
                }
            }
        }

        /// <summary>
        /// 判断运算符是否处于一元上下文中
        /// </summary>
        private bool IsUnaryContext(int tokenIndex)
        {
            if (tokenIndex == 0) // 表达式开头
                return true;

            var prevToken = _tokens[tokenIndex - 1];
            
            // 如果前面是运算符、左括号、逗号，则当前可能是二元运算符
            return prevToken.Type == TokenType.Operator ||
                   prevToken.Type == TokenType.LeftParenthesis ||
                   prevToken.Type == TokenType.Comma;
        }

        /// <summary>
        /// 获取下一个标记而不移动位置（前瞻）
        /// </summary>
        public Token PeekNextToken()
        {
            int savedPosition = _position;
            var tokens = Tokenize();
            _position = savedPosition;
            return tokens.Count > 1 ? tokens[0] : null;
        }

        /// <summary>
        /// 验证标记序列的有效性
        /// </summary>
        public bool ValidateTokens(List<Token> tokens)
        {
            if (tokens == null || tokens.Count == 0)
                return false;

            // 检查括号平衡
            int balance = 0;
            foreach (var token in tokens)
            {
                if (token.Type == TokenType.LeftParenthesis)
                    balance++;
                else if (token.Type == TokenType.RightParenthesis)
                    balance--;

                if (balance < 0)
                    return false; // 右括号多于左括号
            }

            if (balance != 0)
                return false; // 括号不平衡

            // 检查连续运算符
            for (int i = 0; i < tokens.Count - 1; i++)
            {
                if (tokens[i].Type == TokenType.Operator && tokens[i + 1].Type == TokenType.Operator)
                {
                    // 允许某些组合，但不允许大多数连续的运算符
                    if (!AreConsecutiveOperatorsAllowed(tokens[i], tokens[i + 890]))
                        return false;
                }
            }

            return true;
        }

        /// <summary>
        /// 检查是否允许连续的运算符
        /// </summary>
        private bool AreConsecutiveOperatorsAllowed(Token first, Token second)
        {
            // 例如："-" 后跟 "+" 在某些情况下是允许的
            return (first.Value == "-" && second.Value == "+") ||
                   (first.Value == "+" && second.Value == "-");
        }

        /// <summary>
        /// 调试输出标记序列
        /// </summary>
        public void PrintTokens(List<Token> tokens)
        {
            Console.WriteLine("标记序列:");
            foreach (var token in tokens)
            {
                Console.WriteLine(token.ToString());
            }
        }
    }
}