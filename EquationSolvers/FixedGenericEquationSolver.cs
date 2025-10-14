using System;
using System.Collections.Generic;
using System.Linq;
using EquationSolver.Interfaces;
using EquationSolver.Models;

namespace EquationSolver.Solvers
{
    /// <summary>
    /// 修正版通用方程求解器 - 正确的继承关系和接口实现
    /// </summary>
    public class FixedGenericEquationSolver : BaseEquationSolver, INonlinearEquationSolver
    {
        private readonly IMathExpressionParser _parser;
        private dynamic _expressionTree; // dynamic to avoid type name conflicts
        private Dictionary<string, double> _variables = new Dictionary<string, double>();
        private string _normalizedEquation = string.Empty;

        private readonly EquationSolver.Models.SolverConfig _config;

        public FixedGenericEquationSolver(IMathExpressionParser parser)
        {
            _parser = parser ?? throw new ArgumentNullException(nameof(parser));
            _config = new EquationSolver.Models.SolverConfig();
        }

        public FixedGenericEquationSolver(IMathExpressionParser parser, EquationSolver.Models.SolverConfig config)
        {
            _parser = parser ?? throw new ArgumentNullException(nameof(parser));
            _config = config ?? new EquationSolver.Models.SolverConfig();
        }

        #region IEquationSolver 接口实现

        public override EquationSolver.Models.SolveResult Solve()
        {
            try
            {
                var strategy = SelectOptimalSolvingStrategy();
                return ExecuteSolvingStrategy(strategy);
            }
            catch (Exception ex)
            {
                return new EquationSolver.Models.SolveResult { Success = false, Message = $"通用求解器执行失败: {ex.Message}" };
            }
        }

        public override string GetSupportedTypes() => "通用方程求解器 - 支持线性、非线性、超越方程的数值求解";

        protected override void PerformSpecificParsing(string normalizedEquation)
        {
            _normalizedEquation = normalizedEquation ?? string.Empty;
            _expressionTree = _parser.Parse(_normalizedEquation);
            var vars = _parser.ExtractVariables(_normalizedEquation);
            _variables = vars != null ? new Dictionary<string, double>(vars) : new Dictionary<string, double>();
        }

        #endregion

        #region INonlinearEquationSolver 接口实现

        public double SolveByNewtonRaphson(double initialGuess, double tolerance = 1e-15, int maxIterations = 1000)
        {
            if (_variables == null || !_variables.Any())
                throw new InvalidOperationException("变量未初始化，无法使用牛顿法");

            var varName = _variables.Keys.First();
            double x = initialGuess;

            for (int i = 0; i < maxIterations; i++)
            {
                _variables[varName] = x;
                var fx = EvaluateAtPoint(x);
                var dfx = NumericalDerivative(x);
                if (Math.Abs(dfx) < double.Epsilon)
                    throw new InvalidOperationException("导数接近零，牛顿法可能失效");
                var xNew = x - fx / dfx;
                if (Math.Abs(xNew - x) < tolerance) return xNew;
                x = xNew;
            }

            return x;
        }

        public double SolveByBisection(double leftBound, double rightBound, double tolerance = 1e-14)
        {
            double fa = EvaluateAtPoint(leftBound);
            double fb = EvaluateAtPoint(rightBound);
            if (fa * fb > 0) throw new ArgumentException("区间端点函数值同号，无法保证有根");
            double a = leftBound, b = rightBound;
            while ((b - a) / 2.0 > tolerance)
            {
                double m = (a + b) / 2.0;
                double fm = EvaluateAtPoint(m);
                if (Math.Abs(fm) < tolerance) return m;
                if (fa * fm <= 0) b = m; else { a = m; fa = fm; }
            }
            return (a + b) / 2.0;
        }

        public double SolveBySecantMethod(double x0, double x1, double tolerance = 1e-13)
        {
            double xPrev = x0, xCurr = x1;
            for (int i = 0; i < _config.MaxIterations; i++)
            {
                double fPrev = EvaluateAtPoint(xPrev);
                double fCurr = EvaluateAtPoint(xCurr);
                if (Math.Abs(fCurr - fPrev) < 1e-20) throw new InvalidOperationException("两点函数值过于接近，割线法失效");
                double xNext = xCurr - fCurr * (xCurr - xPrev) / (fCurr - fPrev);
                if (Math.Abs(xNext - xCurr) < tolerance) return xNext;
                xPrev = xCurr; xCurr = xNext;
            }
            return xCurr;
        }

        #endregion

        #region 求解策略和方法

        private LocalAnalysis SelectOptimalSolvingStrategy()
        {
            var analysis = AnalyzeEquationCharacteristics();
            if (analysis.VariableCount <= 0) return new LocalAnalysis { Strategy = Strategy.DefaultHybrid };
            if (analysis.VariableCount == 1)
            {
                if (analysis.IsLinear) return new LocalAnalysis { Strategy = Strategy.BracketingMethods };
                if (analysis.HasAnalyticalSolution) return new LocalAnalysis { Strategy = Strategy.Analytical };
                if (analysis.IsWellBehaved) return new LocalAnalysis { Strategy = Strategy.NewtonFamily };
                return new LocalAnalysis { Strategy = Strategy.GlobalOptimization };
            }
            return new LocalAnalysis { Strategy = analysis.IsLinearSystem ? Strategy.LinearSystems : Strategy.OptimizationBased };
        }

        private EquationSolver.Models.SolveResult ExecuteSolvingStrategy(LocalAnalysis chosen)
        {
            try
            {
                switch (chosen.Strategy)
                {
                    case Strategy.Analytical: return TryAnalyticalSolution();
                    case Strategy.NewtonFamily: return ExecuteNewtonFamilyMethods();
                    case Strategy.BracketingMethods: return ExecuteBracketingMethods();
                    case Strategy.GlobalOptimization: return ExecuteGlobalOptimization();
                    case Strategy.LinearSystems: return SolveAsLinearSystem();
                    case Strategy.OptimizationBased: return SolveUsingOptimization();
                    default: return ExecuteHybridApproach();
                }
            }
            catch (Exception ex)
            {
                return new EquationSolver.Models.SolveResult { Success = false, Message = $"执行求解策略时出错: {ex.Message}" };
            }
        }

        private EquationSolver.Models.SolveResult TryAnalyticalSolution()
        {
            var analysis = AnalyzeEquationCharacteristics();
            if (analysis.IsTrivialZeroEquation) return new EquationSolver.Models.SolveResult { Success = true, Message = "方程为恒等式 0 = 0，任何实数都是解" };
            if (analysis.IsConstantEquation)
            {
                var cv = EvaluateAtPoint(0);
                if (Math.Abs(cv) < _config.Tolerance) return new EquationSolver.Models.SolveResult { Success = true, Message = "方程恒成立" };
                return new EquationSolver.Models.SolveResult { Success = false, Message = "方程无解（非常数项不为零）" };
            }
            return new EquationSolver.Models.SolveResult { Success = false, Message = "无法找到解析解，将使用数值方法" };
        }

        private EquationSolver.Models.SolveResult ExecuteNewtonFamilyMethods()
        {
            if (_variables == null || !_variables.Any()) return new EquationSolver.Models.SolveResult { Success = false, Message = "未检测到变量，无法使用牛顿类方法" };
            foreach (var g in GenerateInitialGuesses())
            {
                try
                {
                    var sol = SolveByNewtonRaphson(g, _config.Tolerance, _config.MaxIterations);
                    if (ValidateSolution(sol)) return new EquationSolver.Models.SolveResult { Success = true, Solutions = new List<double> { sol }, Message = $"牛顿法求得解 = {sol:F6}" };
                }
                catch { }
            }
            return ExecuteSecantMethodFallback();
        }

        private EquationSolver.Models.SolveResult ExecuteBracketingMethods()
        {
            var intervals = FindRootIntervals();
            if (!intervals.Any()) return new EquationSolver.Models.SolveResult { Success = false, Message = "未找到有效的根区间" };
            var sols = new List<double>();
            foreach (var (a, b) in intervals.Take(_config.MaxRoots))
            {
                try { var s = SolveByBrentMethodSimplified(a, b); if (ValidateSolution(s)) sols.Add(s); }
                catch { try { var s = SolveByBisection(a, b); if (ValidateSolution(s)) sols.Add(s); } catch { } }
            }
            return sols.Any() ? new EquationSolver.Models.SolveResult { Success = true, Solutions = sols, Message = $"包围方法找到 {sols.Count} 个解" }
                              : new EquationSolver.Models.SolveResult { Success = false, Message = "包围方法未能找到有效解" };
        }

        private EquationSolver.Models.SolveResult ExecuteGlobalOptimization()
        {
            double best = double.NaN; double bestRes = double.MaxValue;
            foreach (var x in GenerateGridSearchPoints()) { var r = Math.Abs(EvaluateAtPoint(x)); if (r < bestRes) { bestRes = r; best = x; } }
            if (!double.IsNaN(best) && bestRes < _config.Tolerance) return new EquationSolver.Models.SolveResult { Success = true, Solutions = new List<double> { best }, Message = $"网格搜索残差 {bestRes:E6}" };
            return new EquationSolver.Models.SolveResult { Success = false, Message = "全局优化方法未能找到满意解" };
        }

        private EquationSolver.Models.SolveResult ExecuteHybridApproach()
        {
            var strategies = new[] { Strategy.NewtonFamily, Strategy.BracketingMethods, Strategy.GlobalOptimization };
            foreach (var s in strategies)
            {
                var res = ExecuteSolvingStrategy(new LocalAnalysis { Strategy = s });
                if (res != null && (res.Success || (res.Solutions != null && res.Solutions.Count > 0))) return res;
            }
            return new EquationSolver.Models.SolveResult { Success = false, Message = "所有混合方法均未能求解该方程" };
        }

        private EquationSolver.Models.SolveResult ExecuteSecantMethodFallback()
        {
            if (_variables == null || !_variables.Any()) return new EquationSolver.Models.SolveResult { Success = false, Message = "未检测到变量，无法使用割线法" };
            try
            {
                var sol = SolveBySecantMethod(-1.0, 1.0, _config.Tolerance);
                if (ValidateSolution(sol)) return new EquationSolver.Models.SolveResult { Success = true, Solutions = new List<double> { sol }, Message = "割线法求得解" };
            }
            catch (Exception ex) { return new EquationSolver.Models.SolveResult { Success = false, Message = $"割线法失败: {ex.Message}" }; }
            return new EquationSolver.Models.SolveResult { Success = false, Message = "割线法未能找到有效解" };
        }

        private double SolveByBrentMethodSimplified(double a, double b, double tolerance = 1e-12)
        {
            double fa = EvaluateAtPoint(a), fb = EvaluateAtPoint(b);
            if (fa * fb > 0) throw new ArgumentException("区间端点函数值同号");
            if (Math.Abs(fa) < Math.Abs(fb)) { (a, b) = (b, a); (fa, fb) = (fb, fa); }
            double c = a, fc = fa, s = b;
            for (int iter = 0; iter < _config.MaxIterations; iter++)
            {
                if (Math.Abs(fb) < tolerance) return b;
                if (fa != fc && fb != fc)
                    s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) + c * fa * fb / ((fc - fa) * (fc - fb));
                else
                    s = b - fb * (b - a) / (fb - fa);
                if (double.IsNaN(s) || s <= Math.Min(a, b) || s >= Math.Max(a, b)) s = (a + b) / 2.0;
                double fs = EvaluateAtPoint(s);
                if (Math.Abs(fs) < tolerance) return s;
                if (fa * fs < 0) { b = s; fb = fs; } else { a = s; fa = fs; }
                if (Math.Abs(a - b) < tolerance) return (a + b) / 2.0;
                if (Math.Abs(fa) < Math.Abs(fb)) { (a, b) = (b, a); (fa, fb) = (fb, fa); }
                c = a; fc = fa;
            }
            return (a + b) / 2.0;
        }

        private double EvaluateAtPoint(double x)
        {
            if (_variables == null || !_variables.Any()) throw new InvalidOperationException("变量未初始化，无法评估表达式");
            var name = _variables.Keys.First(); _variables[name] = x;
            return _expressionTree.Evaluate(_variables);
        }

        private double NumericalDerivative(double x, double h = 1e-6) => (EvaluateAtPoint(x + h) - EvaluateAtPoint(x - h)) / (2.0 * h);
        private bool ValidateSolution(double sol) => !double.IsNaN(sol) && !double.IsInfinity(sol) && Math.Abs(EvaluateAtPoint(sol)) < _config.Tolerance;

        private Models.EquationAnalysis AnalyzeEquationCharacteristics()
        {
            return new Models.EquationAnalysis
            {
                VariableCount = _variables?.Count ?? 0,
                IsLinear = CheckLinearity(),
                HasAnalyticalSolution = CheckAnalyticalSolvability(),
                IsWellBehaved = CheckWellBehavior(),
                IsTrivialZeroEquation = CheckIfTrivialZero(),
                IsConstantEquation = CheckIfConstant(),
                IsLinearSystem = CheckIfLinearSystem()
            };
        }

        private bool CheckLinearity() => !new[] { "^", "sin", "cos", "tan", "log", "exp", "sqrt" }.Any(op => _normalizedEquation.Contains(op, StringComparison.OrdinalIgnoreCase));
        private bool CheckAnalyticalSolvability() => _normalizedEquation.Count(c => "+-*/^".Contains(c)) <= 2;
        private bool CheckWellBehavior()
        {
            try { var pts = new[] { -10.0, 0.0, 10.0 }; var vals = pts.Select(EvaluateAtPoint).ToArray(); for (int i = 0; i < vals.Length - 1; i++) { if (double.IsNaN(vals[i]) || double.IsInfinity(vals[i])) return false; if (Math.Abs(vals[i + 1] - vals[i]) > 1e6) return false; } return true; } catch { return false; }
        }
        private bool CheckIfTrivialZero() => _normalizedEquation.Trim() == "0" || _normalizedEquation.Trim() == "0=0";
        private bool CheckIfConstant() => !_normalizedEquation.Any(char.IsLetter);
        private bool CheckIfLinearSystem() => _normalizedEquation.Contains(",") && _normalizedEquation.Contains("=");

        private IEnumerable<double> GenerateInitialGuesses() { yield return -10.0; yield return -1.0; yield return 0.0; yield return 1.0; yield return 10.0; foreach (var p in FindPotentialCriticalPoints()) yield return p; }
        private IEnumerable<double> FindPotentialCriticalPoints() { for (double x = -10.0; x <= 10.0; x += 1.0) if (Math.Abs(EvaluateAtPoint(x)) < 1e-3) yield return x; }
        private List<(double, double)> FindRootIntervals() { var list = new List<(double, double)>(); var step = 1.0; var range = Math.Min(_config.SearchRange, 1000.0); for (double x = -range; x <= range; x += step) { double x2 = x + step; double f1 = EvaluateAtPoint(x); double f2 = EvaluateAtPoint(x2); if (f1 * f2 <= 0) list.Add((x, x2)); } return list; }
        private IEnumerable<double> GenerateGridSearchPoints() { var range = Math.Min(_config.SearchRange, 100.0); var step = Math.Max(1.0, range / 50.0); for (double x = -range; x <= range; x += step) yield return x; }

        private EquationSolver.Models.SolveResult SolveAsLinearSystem() => new EquationSolver.Models.SolveResult { Success = false, Message = "线性系统求解暂未实现" };
        private EquationSolver.Models.SolveResult SolveUsingOptimization() => new EquationSolver.Models.SolveResult { Success = false, Message = "基于优化的求解暂未实现" };

        private enum Strategy { Analytical, NewtonFamily, BracketingMethods, GlobalOptimization, LinearSystems, OptimizationBased, DefaultHybrid }
        private class LocalAnalysis { public Strategy Strategy; }

        #endregion
    }
}
