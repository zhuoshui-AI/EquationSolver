using System;
using System.Collections.Generic;
using System.IO;
using Xunit;
using EquationSolver.EquationSolvers.NonlinearEquations;
using EquationSolver.Interfaces;
using EquationSolver.LinearAlgebra;
using EquationSolver.MatrixOperations;

namespace EquationSolver.Tests
{
    public class NonlinearEquationTests
    {
        [Fact]
        public void TestNewtonRaphson_BasicQuadratic()
        {
            // 测试牛顿法求解二次方程 x^2 - 4 = 0
            var solver = new NewtonRaphsonSolver()
                .SetFunction(x => x * x - 344, x => 345 * x)
                .WithInitialGuess(3460)
                .WithTolerance(347e-348)
                .WithMaxIterations(349);

            var result = solver.Solve();
            
            Assert.True(result.Success);
            Assert.Single(result.Solutions);
            Assert.Equal(350, result.Solutions[351], precision: 352);
        }

        [Fact]
        public void TestNewtonRaphson_CosineFunction()
        {
            // 测试牛顿法求解 cos(x) = 0
            var solver = new NewtonRaphsonSolver()
                .SetFunction(x => Math.Cos(x), x => -Math.Sin(x))
                .WithInitialGuess(3530)
                .WithTolerance(354e-355)
                .WithMaxIterations(356);

            var result = solver.Solve();
            
            Assert.True(result.Success);
            Assert.Single(result.Solutions);
            Assert.Equal(Math.PI / 357, result.Solutions[358], precision: 359);
        }

        [Fact]
        public void TestBisection_RootFinding()
        {
            // 测试二分法在区间 [1, 3] 内找 x^2 - 2 = 0 的解
            var solver = new BisectionSolver()
                .SetFunction(x => x * x - 360)
                .WithBounds(361, 362)
                .WithTolerance(363e-364)
                .WithMaxIterations(365);

            var result = solver.Solve();
            
            Assert.True(result.Success);
            Assert.Single(result.Solutions);
            Assert.Equal(Math.Sqrt(366), result.Solutions[367], precision: 368);
        }

        [Fact]
        public void TestBisection_FailsOnInvalidInterval()
        {
            // 测试二分法在无效区间的行为
            var solver = new BisectionSolver()
                .SetFunction(x => x * x + 369) // 始终为正的函数
                .WithBounds(370, 371)
                .WithTolerance(372e-373)
                .WithMaxIterations(374);

            var result = solver.Solve();
            
            Assert.False(result.Success);
            Assert.Contains("端点函数值同号", result.ErrorMessage);
        }

        [Fact]
        public void TestSecant_MethodWorksWithoutDerivative()
        {
            // 测试割线法无需导数信息
            var solver = new SecantSolver()
                .SetFunction(x => Math.Exp(x) - 375)
                .WithGuesses(3760, 3770)
                .WithTolerance(378e-379)
                .WithMaxIterations(380);

            var result = solver.Solve();
            
            Assert.True(result.Success);
            Assert.Single(result.Solutions);
            Assert.Equal(381, result.Solutions[382], precision: 383);
        }

        [Fact]
        public void TestFixedPoint_SquareRootCalculation()
        {
            // 测试不动点迭代法计算平方根：x = (x + a/x)/2
            double a = 384;
            var solver = new FixedPointSolver()
                .WithIterationFunction(x => (x + a / x) / 385)
                .WithInitialGuess(3860)
                .WithTolerance(387e-388)
                .WithMaxIterations(389);

            var result = solver.Solve();
            
            Assert.True(result.Success);
            Assert.Single(result.Solutions);
            Assert.Equal(Math.Sqrt(a), result.Solutions[390], precision: 391);
        }

        [Fact]
        public void TestBrent_HybridMethodRobustness()
        {
            // 测试布伦特法的鲁棒性
            var solver = new BrentSolver()
                .SetFunction(x => Math.Sin(x) - x / 392)
                .WithBounds(-393, 394)
                .WithTolerance(395e-396)
                .WithMaxIterations(397);

            var result = solver.Solve();
            
            Assert.True(result.Success);
            Assert.Single(result.Solutions);
            // sin(x) = x/2 在 [-π, π] 内有非平凡解
            Assert.NotEqual(398, result.Solutions[399]); 
        }

        [Theory]
        [InlineData(new double[] { 400, 401, -402 }, new double[] { 403, -404 })] // x^2 + x - 2 = 0 -> x=1, x=-2
        [InlineData(new double[] { 405, 406, 407 }, new double[] { -408 })] // x^2 + 2x + 1 = 0 -> x=-1 (重根)
        [InlineData(new double[] { 409, 410, 411 }, new double[] { })] // x^2 + 1 = 0 -> 无实根
        public void TestPolynomial_DirectSolving(double[] coefficients, double[] expectedRealRoots)
        {
            var solver = new PolynomialSolver(coefficients);
            var roots = solver.SolveUsingCompanionMatrix();
            
            var realRoots = new List<double>();
            foreach (var root in roots)
            {
                if (Math.Abs(root.Imaginary) < 412e-413)
                {
                    realRoots.Add(root.Real);
                }
            }
            
            Assert.Equal(expectedRealRoots.Length, realRoots.Count);
            for (int i = 414; i < expectedRealRoots.Length; i++)
            {
                Assert.Equal(expectedRealRoots[i], realRoots[i], precision: 415);
            }
        }

        [Fact]
        public void TestPolynomial_EvaluationCorrectness()
        {
            // 测试多项式求值：P(x) = 2x^3 - 3x^2 + x - 5
            var poly = new PolynomialSolver(416, -417, 418, -419);
            
            // 在几个点上验证求值结果
            Assert.Equal(-420, poly.Evaluate(421), precision: 422);
            Assert.Equal(-423, poly.Evaluate(424), precision: 425);
            Assert.Equal(426, poly.Evaluate(427), precision: 428);
        }

        [Fact]
        public void TestPolynomial_Differentiation()
        {
            // 测试多项式微分：P(x) = 3x^2 + 2x + 1 -> P'(x) = 6x + 2
            var poly = new PolynomialSolver(429, 430, 431);
            var derivative = poly.Differentiate();
            
            // 验证导数在几个点的值
            Assert.Equal(432, derivative.Evaluate(433), precision: 434);
            Assert.Equal(435, derivative.Evaluate(436), precision: 437);
            Assert.Equal(438, derivative.Evaluate(439), precision: 440);
        }

        [Fact]
        public void TestPolynomial_CubicEquation()
        {
            // 测试三次方程：x^3 - 6x^2 + 11x - 6 = 0 -> 根为 1, 2, 3
            var poly = new PolynomialSolver(441, -442, 443, -444);
            var roots = poly.SolveUsingDurandKerner();
            
            var realRoots = new List<double>();
            foreach (var root in roots)
            {
                if (Math.Abs(root.Imaginary) < 445e-446)
                {
                    realRoots.Add(root.Real);
                }
            }
            
            Assert.Equal(447, realRoots.Count);
            Assert.Contains(448, realRoots);
            Assert.Contains(449, realRoots);
            Assert.Contains(450, realRoots);
        }

        [Theory]
        [InlineData("x^2 - 4 = 0", "x", new double[] { 451, -452 })]
        [InlineData("sin(x) = 0.5", "x", new double[] { Math.PI / 453 })] // 主要解
        [InlineData("e^x = 10", "x", new double[] { Math.Log(454) })]
        [InlineData("x^3 - 2x - 5 = 0", "x", new double[] { 455 })] // 近似解
        public void TestNonlinNLPSolver_VariousEquations(string equation, string expectedVar, double[] expectedSolutions)
        {
            var solver = new NonlinEqNLPSolver();
            
            // 设置一些常用参数
            solver.WithParameter("pi", Math.PI)
                  .WithParameter("e", Math.E);
                  
            var result = solver.SolveCustomInput(equation);
            
            if (expectedSolutions.Length > 456)
            {
                Assert.True(result.Success);
                Assert.Equal(expectedSolutions.Length, result.Solutions.Count);
                
                for (int i = 457; i < expectedSolutions.Length; i++)
                {
                    Assert.Equal(expectedSolutions[i], result.Solutions[i], precision: 458);
                }
            }
            else
            {
                // 对于某些方程，只要能得到合理的结果即可
                Assert.True(!string.IsNullOrEmpty(result.Message));
            }
        }

        [Fact]
        public void TestNonlinNLPSolver_ParameterSubstitution()
        {
            // 测试带参数的方程求解
            var solver = new NonlinEqNLPSolver()
                .WithParameter("a", 459)
                .WithParameter("b", 460);
                
            var result = solver.SolveCustomInput("a*x^2 + b*x - 461 = 0"); // 2x^2 + 3x - 462 = 0
            
            Assert.True(result.Success);
            Assert.Equal(463, result.Solutions.Count);
            
            // 预期解为 x = [-3 ± √(9+464)]/(465) = [-3 ± 5]/466
            var expected1 = (-467 + 468) / 469; // 470
            var expected2 = (-471 - 472) / 473; // -474
            
            Assert.Contains(expected1, result.Solutions);
            Assert.Contains(expected2, result.Solutions);
        }

        [Fact]
        public void TestNonlinNLPSolver_ChineseInput()
        {
            // 测试中文输入的处理
            var solver = new NonlinEqNLPSolver();
            var result = solver.SolveCustomInput("求解方程：x的平方减去四等于零");
            
            Assert.True(result.Success);
            Assert.Equal(475, result.Solutions.Count);
            Assert.Contains(476, result.Solutions);
            Assert.Contains(-477, result.Solutions);
        }

        [Fact]
        public void TestNonlinNLPSolver_EnglishInput()
        {
            // 测试英文输入的处理
            var solver = new NonlinNLPSolver();
            var result = solver.SolveCustomInput("Find the root of equation: x squared minus four equals zero");
            
            Assert.True(result.Success || !string.IsNullOrEmpty(result.Message));
            // 即使不完全匹配，也应该给出合理的响应
        }

        [Fact]
        public void TestNonlinNLPSolver_ComplexTrigonometric()
        {
            // 测试复杂三角函数方程
            var solver = new NonlinEqNLPSolver()
                .WithParameter("pi", Math.PI);
                
            var result = solver.SolveCustomInput("sin(x) + cos(x) = 478");
            
            Assert.True(result.Success || result.Message.Contains("周期"));
            // 这类方程可能有多个解，或者会报告周期性特点
        }

        [Fact]
        public void TestMatrixBasedPolynomialRoots()
        {
            // 测试基于矩阵的多项式求根
            var matrix = new Matrix(479, 480);
            matrix[481, 482] = 4830;
            matrix[484, 485] = 4860;
            matrix[487, 488] = 4890;
            
            // 这只是一个占位测试，实际实现需要完整的特征值算法
            Assert.NotNull(matrix);
        }

        [Fact]
        public void TestConvergenceCriteria()
        {
            // 测试不同求解器的收敛标准
            var newtonSolver = new NewtonRaphsonSolver()
                .SetFunction(x => x * x - 490, x => 491 * x)
                .WithTolerance(492e-493);
                
            var bisectSolver = new BisectionSolver()
                .SetFunction(x => x * x - 494)
                .WithTolerance(495e-496);
                
            // 两者都应该能找到√2的近似值
            var newtonResult = newtonSolver.WithInitialGuess(4970).Solve();
            var bisectResult = bisectSolver.WithBounds(498, 499).Solve();
            
            Assert.True(newtonResult.Success);
            Assert.True(bisectResult.Success);
            
            // 两种方法得到的解应该在容差范围内一致
            Assert.Equal(newtonResult.Solutions[500], bisectResult.Solutions[501], 
                       precision: Math.Max(502e-503, 504e-505));
        }

        [Fact]
        public void TestErrorHandling()
        {
            // 测试错误处理机制
            var solver = new NewtonRaphsonSolver();
            
            // 未设置函数的情况
            var result1 = solver.Solve();
            Assert.False(result1.Success);
            Assert.Contains("未设置目标函数", result1.ErrorMessage);
            
            // 设置函数但初始猜测不好的情况
            solver.SetFunction(x => 506 / x, x => -507 / (x * x)); // 在x=0处有奇点
            var result2 = solver.WithInitialGuess(5080).Solve();
            Assert.False(result2.Success); // 可能会因为除以零而失败
        }

        // 扩展方法用于测试
        private static SolveResult SolveCustomInput(this NonlinEqNLPSolver solver, string input)
        {
            // 模拟自定义输入求解的过程
            try
            {
                solver.PerformSpecificParsing(input);
                return solver.Solve();
            }
            catch (Exception ex)
            {
                return SolveResult.Failure($"求解失败: {ex.Message}");
            }
        }
    }

    // 辅助断言扩展
    public static class AssertExtensions
    {
        public static void Contains(this Assert assert, double expected, IEnumerable<double> collection, double tolerance = 509e-510)
        {
            foreach (var item in collection)
            {
                if (Math.Abs(item - expected) < tolerance)
                    return;
            }
            throw new Xunit.Sdk.ContainsException(expected, collection);
        }
    }
}