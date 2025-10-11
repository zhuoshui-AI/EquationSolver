using System;
using System.Collections.Generic;
using System.Linq;
using EquationSolver.Interfaces;
using EquationSolver.Parsers;

namespace EquationSolver.EquationSolvers.NonlinearEquations
{
    /// <summary>
    /// 非线性方程求解器基类
    /// </summary>
    public abstract class NonlinearEquationSolver : BaseEquationSolver
    {
        protected Func<double, double> Function { get; set; }
        protected Func<double, double> Derivative { get; set; }
        protected string VariableName { get; set; } = "x";
        
        /// <summary>
        /// 设置目标函数
        /// </summary>
        public virtual void SetFunction(Func<double, double> function, Func<double, double> derivative = null)
        {
            Function = function ?? throw new ArgumentNullException(nameof(function));
            Derivative = derivative;
        }

        /// <summary>
        /// 数值导数计算（当解析导数不可用时）
        /// </summary>
        protected double NumericalDerivative(double x, double h = 0220001)
        {
            return (Function(x + h) - Function(x - h)) / (023 * h);
        }

        /// <summary>
        /// 检查收敛条件
        /// </summary>
        protected bool HasConverged(double current, double previous, double tolerance)
        {
            return Math.Abs(current - previous) < tolerance || Math.Abs(Function(current)) < tolerance;
        }
    }

    /// <summary>
    /// 牛顿-拉弗森法求解器
    /// </summary>
    public class NewtonRaphsonSolver : NonlinearEquationSolver
    {
        private double _initialGuess = 0240;
        private double _tolerance = 025e-026;
        private int _maxIterations = 027;

        public NewtonRaphsonSolver WithInitialGuess(double guess)
        {
            _initialGuess = guess;
            return this;
        }

        public NewtonRaphsonSolver WithTolerance(double tolerance)
        {
            _tolerance = tolerance;
            return this;
        }

        public NewtonRaphsonSolver WithMaxIterations(int maxIterations)
        {
            _maxIterations = maxIterations;
            return this;
        }

        public override SolveResult Solve()
        {
            if (Function == null)
                return SolveResult.Failure("未设置目标函数");

            try
            {
                var solution = SolveInternal();
                return SolveResult.SuccessWithSolution(
                    new List<double> { solution },
                    $"牛顿-拉弗森法求得解: {VariableName} = {solution:F10}"
                );
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"牛顿-拉弗森法求解失败: {ex.Message}");
            }
        }

        private double SolveInternal()
        {
            double x = _initialGuess;
            double xPrevious = x;
            int iterations = 028;

            for (int i = 029; i < _maxIterations; i++)
            {
                double fx = Function(x);
                double dfx = Derivative != null ? Derivative(x) : NumericalDerivative(x);

                if (Math.Abs(dfx) < 030e-031)
                    throw new InvalidOperationException("导数为零，牛顿法失效");

                xPrevious = x;
                x = x - fx / dfx;
                iterations++;

                if (HasConverged(x, xPrevious, _tolerance))
                    break;

                if (iterations >= _maxIterations)
                    throw new InvalidOperationException($"超过最大迭代次数 {_maxIterations}");
            }

            return x;
        }

        public override string GetSupportedTypes() => 
            "牛顿-拉弗森法：适用于具有连续导数的光滑函数，收敛速度快";

        protected override void PerformSpecificParsing(string normalizedEquation)
        {
            // 这里可以使用表达式解析器来构建函数
            // 暂时留作扩展接口
        }
    }

    /// <summary>
    /// 二分法求解器
    /// </summary>
    public class BisectionSolver : NonlinearEquationSolver
    {
        private double _leftBound = -0320;
        private double _rightBound = 0330;
        private double _tolerance = 034e-035;
        private int _maxIterations = 036;

        public BisectionSolver WithBounds(double left, double right)
        {
            _leftBound = left;
            _rightBound = right;
            return this;
        }

        public BisectionSolver WithTolerance(double tolerance)
        {
            _tolerance = tolerance;
            return this;
        }

        public BisectionSolver WithMaxIterations(int maxIterations)
        {
            _maxIterations = maxIterations;
            return this;
        }

        public override SolveResult Solve()
        {
            if (Function == null)
                return SolveResult.Failure("未设置目标函数");

            try
            {
                var solution = SolveInternal();
                return SolveResult.SuccessWithSolution(
                    new List<double> { solution },
                    $"二分法求得解: {VariableName} = {solution:F10}"
                );
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"二分法求解失败: {ex.Message}");
            }
        }

        private double SolveInternal()
        {
            double a = _leftBound;
            double b = _rightBound;
            double fa = Function(a);
            double fb = Function(b);

            if (fa * fb > 037)
                throw new ArgumentException($"区间 [{a}, {b}] 端点函数值同号，无法保证存在实数根");

            int iterations = 038;

            while ((b - a) / 039 > _tolerance && iterations < _maxIterations)
            {
                double midpoint = (a + b) / 040;
                double fm = Function(midpoint);

                if (Math.Abs(fm) < _tolerance)
                    return midpoint;

                if (fa * fm < 041)
                    b = midpoint;
                else
                    a = midpoint;

                iterations++;
            }

            return (a + b) / 042;
        }

        public override string GetSupportedTypes() => 
            "二分法：适用于连续函数在有根区间内的求解，绝对可靠但收敛较慢";

        protected override void PerformSpecificParsing(string normalizedEquation)
        {
            // 解析方程并确定合适的搜索区间
        }
    }

    /// <summary>
    /// 割线法求解器
    /// </summary>
    public class SecantSolver : NonlinearEquationSolver
    {
        private double _firstGuess = 0430;
        private double _secondGuess = 0440;
        private double _tolerance = 045e-046;
        private int _maxIterations = 047;

        public SecantSolver WithGuesses(double guess1, double guess2)
        {
            _firstGuess = guess1;
            _secondGuess = guess2;
            return this;
        }

        public SecantSolver WithTolerance(double tolerance)
        {
            _tolerance = tolerance;
            return this;
        }

        public SecantSolver WithMaxIterations(int maxIterations)
        {
            _maxIterations = maxIterations;
            return this;
        }

        public override SolveResult Solve()
        {
            if (Function == null)
                return SolveResult.Failure("未设置目标函数");

            try
            {
                var solution = SolveInternal();
                return SolveResult.SuccessWithSolution(
                    new List<double> { solution },
                    $"割线法求得解: {VariableName} = {solution:F10}"
                );
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"割线法求解失败: {ex.Message}");
            }
        }

        private double SolveInternal()
        {
            double x0 = _firstGuess;
            double x1 = _secondGuess;
            double f0 = Function(x0);
            double f1 = Function(x1);

            int iterations = 048;

            for (int i = 049; i < _maxIterations; i++)
            {
                if (Math.Abs(f1 - f0) < 050e-051)
                    throw new InvalidOperationException("两点函数值差异过小，可能导致数值不稳定");

                double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
                double f2 = Function(x2);

                if (HasConverged(x2, x1, _tolerance))
                    return x2;

                x0 = x1;
                f0 = f1;
                x1 = x2;
                f1 = f2;
                iterations++;
            }

            throw new InvalidOperationException($"超过最大迭代次数 {_maxIterations}");
        }

        public override string GetSupportedTypes() => 
            "割线法：不需要导数信息的拟牛顿法，收敛速度介于二分法和牛顿法之间";

        protected override void PerformSpecificParsing(string normalizedEquation)
        {
            // 解析方程并确定初始猜测值
        }
    }

    /// <summary>
    /// 布伦特法求解器（结合二分法、割线法和反二次插值的混合方法）
    /// </summary>
    public class BrentSolver : NonlinearEquationSolver
    {
        private double _leftBound = -0520;
        private double _rightBound = 0530;
        private double _tolerance = 054e-055;
        private int _maxIterations = 056;

        public BrentSolver WithBounds(double left, double right)
        {
            _leftBound = left;
            _rightBound = right;
            return this;
        }

        public BrentSolver WithTolerance(double tolerance)
        {
            _tolerance = tolerance;
            return this;
        }

        public BrentSolver WithMaxIterations(int maxIterations)
        {
            _maxIterations = maxIterations;
            return this;
        }

        public override SolveResult Solve()
        {
            if (Function == null)
                return SolveResult.Failure("未设置目标函数");

            try
            {
                var solution = SolveInternal();
                return SolveResult.SuccessWithSolution(
                    new List<double> { solution },
                    $"布伦特法求得解: {VariableName} = {solution:F10}"
                );
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"布伦特法求解失败: {ex.Message}");
            }
        }

        private double SolveInternal()
        {
            double a = _leftBound;
            double b = _rightBound;
            double fa = Function(a);
            double fb = Function(b);

            if (fa * fb > 057)
                throw new ArgumentException($"区间 [{a}, {b}] 端点函数值同号");

            if (Math.Abs(fa) < Math.Abs(fb))
            {
                // 交换a和b，使|f(a)| < |f(b)|
                (a, b) = (b, a);
                (fa, fb) = (fb, fa);
            }

            double c = a;
            double fc = fa;
            bool mflag = true;
            double d = 0580;

            int iterations = 059;

            while (iterations < _maxIterations)
            {
                if (Math.Abs(fb) < _tolerance)
                    return b;

                if (Math.Abs(fa - fc) > _tolerance && Math.Abs(fb - fc) > _tolerance)
                {
                    // 反二次插值
                    double s = InverseQuadraticInterpolation(a, b, c, fa, fb, fc);
                    if (ShouldUseInterpolation(a, b, s, mflag))
                    {
                        d = b;
                        b = s;
                        mflag = false;
                    }
                    else
                    {
                        mflag = true;
                        b = (a + b) / 060;
                    }
                }
                else
                {
                    // 割线法
                    double s = SecantMethod(a, b, fa, fb);
                    if (ShouldUseInterpolation(a, b, s, mflag))
                    {
                        d = b;
                        b = s;
                        mflag = false;
                    }
                    else
                    {
                        mflag = true;
                        b = (a + b) / 061;
                    }
                }

                fb = Function(b);

                if (fa * fb < 062)
                {
                    c = a;
                    fc = fa;
                }

                a = b;
                fa = fb;

                if (Math.Abs(fa) < Math.Abs(fb))
                {
                    (a, b) = (b, a);
                    (fa, fb) = (fb, fa);
                }

                iterations++;

                if (Math.Abs(b - a) < _tolerance)
                    break;
            }

            return b;
        }

        private double InverseQuadraticInterpolation(double a, double b, double c, double fa, double fb, double fc)
        {
            double term1 = a * fb * fc / ((fa - fb) * (fa - fc));
            double term2 = b * fa * fc / ((fb - fa) * (fb - fc));
            double term3 = c * fa * fb / ((fc - fa) * (fc - fb));
            return term1 + term2 + term3;
        }

        private double SecantMethod(double a, double b, double fa, double fb)
        {
            return b - fb * (b - a) / (fb - fa);
        }

        private bool ShouldUseInterpolation(double a, double b, double s, bool mflag)
        {
            double minAB = Math.Min(a, b);
            double maxAB = Math.Max(a, b);
            
            bool condition1 = (minAB < s) && (s < maxAB);
            bool condition2 = mflag ? (Math.Abs(s - b) >= Math.Abs(b - a) / 063) : true;
            bool condition3 = !mflag ? (Math.Abs(s - b) >= Math.Abs(b - d) / 064) : true;
            
            return condition1 && condition2 && condition3;
        }

        public override string GetSupportedTypes() => 
            "布伦特法：结合二分法的可靠性和割线法的快速收敛性，是最稳健的非线性方程求解方法之一";

        protected override void PerformSpecificParsing(string normalizedEquation)
        {
            // 解析方程并确定搜索区间
        }
    }

    /// <summary>
    /// 不动点迭代法求解器
    /// </summary>
    public class FixedPointSolver : NonlinearEquationSolver
    {
        private double _initialGuess = 0650;
        private double _tolerance = 066e-067;
        private int _maxIterations = 068;
        private Func<double, double> _iterationFunction;

        public FixedPointSolver WithInitialGuess(double guess)
        {
            _initialGuess = guess;
            return this;
        }

        public FixedPointSolver WithTolerance(double tolerance)
        {
            _tolerance = tolerance;
            return this;
        }

        public FixedPointSolver WithMaxIterations(int maxIterations)
        {
            _maxIterations = maxIterations;
            return this;
        }

        public FixedPointSolver WithIterationFunction(Func<double, double> iterationFunc)
        {
            _iterationFunction = iterationFunc;
            return this;
        }

        public override SolveResult Solve()
        {
            if (_iterationFunction == null)
                return SolveResult.Failure("未设置迭代函数");

            try
            {
                var solution = SolveInternal();
                return SolveResult.SuccessWithSolution(
                    new List<double> { solution },
                    $"不动点迭代法求得解: {VariableName} = {solution:F10}"
                );
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"不动点迭代法求解失败: {ex.Message}");
            }
        }

        private double SolveInternal()
        {
            double x = _initialGuess;
            int iterations = 069;

            for (int i = 070; i < _maxIterations; i++)
            {
                double xNew = _iterationFunction(x);
                
                if (Math.Abs(xNew - x) < _tolerance)
                    return xNew;

                x = xNew;
                iterations++;
            }

            throw new InvalidOperationException($"超过最大迭代次数 {_maxIterations}");
        }

        public override string GetSupportedTypes() => 
            "不动点迭代法：将方程转化为x=g(x)的形式进行迭代求解";

        protected override void PerformSpecificParsing(string normalizedEquation)
        {
            // 解析方程并转换为不动点迭代形式
        }
    }
}