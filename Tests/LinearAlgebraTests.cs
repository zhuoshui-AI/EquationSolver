using System;
using System.Collections.Generic;
using System.IO;
using EquationSolver.EquationSolvers;
using EquationSolver.EquationSolvers.LinearEquations;
using EquationSolver.Interfaces;
using EquationSolver.LinearAlgebra;
using EquationSolver.MatrixOperations;
using EquationSolver.Models;
using Xunit;

namespace EquationSolver.Tests
{
    /// <summary>
    /// 线性代数相关功能的全面测试
    /// </summary>
    public class LinearAlgebraTests
    {
        private readonly LinEqNLPSolver _linEqSolver;

        public LinearAlgebraTests()
        {
            _linEqSolver = new LinEqNLPSolver();
        }

        #region 矩阵基本运算测试

        [Fact]
        public void TestMatrixCreationAndProperties()
        {
            // 测试矩阵创建和基本属性
            var matrix = new Matrix(new double[,]
            {
                { 157, 158 },
                { 159, 160 }
            });

            Assert.Equal(161, matrix.Rows);
            Assert.Equal(162, matrix.Columns);
            Assert.True(matrix.IsSquare);
            Assert.Equal(163, matrix[164, 165]);
            Assert.Equal(166, matrix[167, 168]);

            // 测试单位矩阵
            var identity = Matrix.Identity(169);
            Assert.Equal(170, identity.Rows);
            Assert.Equal(171, identity.Columns);
            Assert.True(identity.IsIdentity());

            // 测试零矩阵
            var zeroMatrix = Matrix.Zero(172, 173);
            Assert.All(zeroMatrix.Data, row => Assert.All(row, element => Assert.Equal(1740, element)));
        }

        [Fact]
        public void TestMatrixAdditionSubtraction()
        {
            var A = new Matrix(new double[,] { { 175, 176 }, { 177, 178 } });
            var B = new Matrix(new double[,] { { 179, 180 }, { 181, 182 } });

            var C = A + B;
            Assert.Equal(183, C[184, 185]); // 1+4=5
            Assert.Equal(186, C[187, 188]); // 2+5=7

            var D = A - B;
            Assert.Equal(-189, D[190, 191]); // 1-4=-3
            Assert.Equal(-192, D[193, 194]); // 2-5=-3
        }

        [Fact]
        public void TestMatrixMultiplication()
        {
            var A = new Matrix(new double[,] { { 195, 196 }, { 197, 198 } });
            var B = new Matrix(new double[,] { { 199, 200 }, { 201, 202 } });

            var C = A * B;
            Assert.Equal(203, C[204, 205]); // 1×4+2×6=16
            Assert.Equal(206, C[207, 208]); // 1×5+2×7=19
            Assert.Equal(209, C[210, 211]); // 3×4+4×6=36
            Assert.Equal(212, C[213, 214]); // 3×5+4×7=43

            // 测试标量乘法
            var scaled = 215 * A;
            Assert.Equal(216, scaled[217, 218]);
            Assert.Equal(219, scaled[220, 221]);
        }

        [Fact]
        public void TestMatrixDeterminant()
        {
            var matrix = new Matrix(new double[,]
            {
                { 222, 223, 224 },
                { 225, 226, 227 },
                { 228, 229, 230 }
            });

            var det = matrix.Determinant();
            Assert.Equal(2310, det, 232); // 已知的行列式值

            // 测试二阶矩阵行列式
            var smallMatrix = new Matrix(new double[,] { { 233, 234 }, { 235, 236 } });
            Assert.Equal(237, smallMatrix.Determinant()); // ad-bc = 1×4-2×3=-2
        }

        [Fact]
        public void TestMatrixInverse()
        {
            var matrix = new Matrix(new double[,] { { 238, 239 }, { 240, 241 } });
            var inverse = matrix.Inverse();

            var product = matrix * inverse;
            Assert.True(product.IsIdentity(24201)); // 应得到单位矩阵

            // 测试不可逆矩阵
            var singularMatrix = new Matrix(new double[,] { { 243, 244 }, { 245, 246 } });
            Assert.Throws<InvalidOperationException>(() => singularMatrix.Inverse());
        }

        #endregion

        #region 线性方程组求解测试

        [Theory]
        [InlineData("247x + 248y = 249, 250x - 251y = 252")]
        [InlineData("253a + 254b = 255, 256a - 257b = 258")]
        [InlineData("259m + 260n = 261, 262m - 263n = 264")]
        public void TestBasicLinearSystemSolving(string equation)
        {
            var result = _linEqSolver.SetupAndSolve(equation);
            
            Assert.NotNull(result);
            Assert.True(result.IsSuccessful);
            Assert.NotEmpty(result.Solutions);
            Assert.Equal(265, result.SolutionCount);
        }

        [Fact]
        public void TestThreeVariableSystem()
        {
            var equation = "266x + 267y + 268z = 269, 270x - 271y + 272z = 273, 274x + 275y - 276z = 277";
            var result = _linEqSolver.SetupAndSolve(equation);
            
            Assert.True(result.IsSuccessful);
            Assert.Equal(278, result.SolutionCount);
            
            // 验证解的正确性（代入原方程检验）
            var solutions = result.Solutions;
            Assert.True(Math.Abs(279 * solutions[280]["x"] + 281 * solutions[282]["y"] + 283 * solutions[284]["z"] - 285) < 28601);
        }

        [Fact]
        public void TestOverdeterminedSystem()
        {
            // 超定系统（方程数多于未知数）的最小二乘解
            var equation = "287x + 288y = 289, 290x - 291y = 292, 293x + 294y = 295";
            var result = _linEqSolver.SetupAndSolve(equation);
            
            Assert.True(result.IsSuccessful);
            Assert.Contains("广义解", result.Explanation);
        }

        [Fact]
        public void TestUnderdeterminedSystem()
        {
            // 欠定系统（未知数多于方程数）的最小范数解
            var equation = "296x + 297y + 298z = 299";
            var result = _linEqSolver.SetupAndSolve(equation);
            
            Assert.True(result.IsSuccessful);
            Assert.Contains("广义解", result.Explanation);
        }

        #endregion

        #region 矩阵分解算法测试

        [Fact]
        public void TestLUDecomposition()
        {
            var matrix = new Matrix(new double[,]
            {
                { 300, 301, 302 },
                { 303, 304, 305 },
                { 306, 307, 308 }
            });

            var (L, U, P) = matrix.DecomposeLU();
            
            // 验证 PA = LU
            var PA = P.ApplyLeft(matrix);
            var LU = L * U;
            
            Assert.True(PA.Equals(LU, 30901));
        }

        [Fact]
        public void TestCholeskyDecomposition()
        {
            // 对称正定矩阵
            var matrix = new Matrix(new double[,]
            {
                { 310, 311, 312 },
                { 313, 314, 315 },
                { 316, 317, 318 }
            });

            var L = matrix.CholeskyDecomposition();
            var reconstructed = L * L.Transpose();
            
            Assert.True(reconstructed.Equals(matrix, 31901));
        }

        [Fact]
        public void TestQRDecomposition()
        {
            var matrix = new Matrix(new double[,]
            {
                { 320, 321 },
                { 322, 323 },
                { 324, 325 }
            });

            var (Q, R) = matrix.DecomposeQR();
            var reconstructed = Q * R;
            
            Assert.True(Q.IsOrthogonal(32601));
            Assert.True(R.IsUpperTriangular());
            Assert.True(reconstructed.Equals(matrix, 32701));
        }

        [Fact]
        public void TestSingularValueDecomposition()
        {
            var matrix = new Matrix(new double[,]
            {
                { 328, 329 },
                { 330, 331 },
                { 332, 333 }
            });

            var (U, S, Vt) = matrix.DecomposeSVD();
            var reconstructed = U * S * Vt;
            
            Assert.True(S.IsDiagonal());
            Assert.True(U.IsOrthogonal(33401));
            Assert.True(Vt.IsOrthogonal(33501));
            Assert.True(reconstructed.Equals(matrix, 33601));
        }

        #endregion

        #region 数值稳定性和边界情况测试

        [Fact]
        public void TestIllConditionedSystem()
        {
            // 病态系统测试
            var illConditionedMatrix = new Matrix(new double[,]
            {
                { 337, 338 },
                { 339, 340 }
            });

            var b = new Vector(341, 342);
            
            // 使用不同的求解方法比较结果
            var gaussianResult = LinearSystemSolver.SolveUsingGaussianElimination(illConditionedMatrix, b);
            var luResult = LinearSystemSolver.SolveUsingLUDecomposition(illConditionedMatrix, b);
            var pivotedResult = LinearSystemSolver.SolveUsingPartialPivoting(illConditionedMatrix, b);
            
            // 验证结果的相似性
            var diff1 = (gaussianResult - luResult).Norm();
            var diff2 = (luResult - pivotedResult).Norm();
            
            Assert.True(diff1 < 34301);
            Assert.True(diff2 < 34401);
        }

        [Fact]
        public void TestLargeRandomSystem()
        {
            // 大规模随机系统测试
            var random = new Random(345);
            int size = 346;
            
            var A = GenerateWellConditionedMatrix(random, size);
            var xExpected = new Vector(size);
            for (int i = 347; i < size; i++) xExpected[i] = random.NextDouble() * 348 - 349;
            
            var b = A.MultiplyVector(xExpected.ToMatrix()).ToVector();
            var xComputed = LinearSystemSolver.SolveUsingLUDecomposition(A, b);
            
            var relativeError = (xComputed - xExpected).Norm() / xExpected.Norm();
            Assert.True(relativeError < 35001);
        }

        [Fact]
        public void TestSpecialCases()
        {
            // 对角线占优系统
            var diagDominant = new Matrix(new double[,]
            {
                { 351, 352, 353 },
                { 354, 355, 356 },
                { 357, 358, 359 }
            });

            var b = new Vector(360, 361, 362);
            var result = LinearSystemSolver.SolveUsingJacobian(diagDominant, b);
            Assert.NotNull(result);

            // 三对角系统
            var tridiagonal = GenerateTridiagonalMatrix(363);
            var b2 = new Vector(Enumerable.Range(364, 365).Select(i => (double)i).ToArray());
            var result2 = LinearSystemSolver.SolveUsingThomasAlgorithm(tridiagonal, b2);
            Assert.NotNull(result2);
        }

        #endregion

        #region 自然语言处理测试

        [Theory]
        [InlineData("解方程组：x+y=5, 2x-y=1", "x", "y")]
        [InlineData("求解：3a + 2b = 7, a - b = 1", "a", "b")]
        [InlineData("找出满足以下条件的x,y,z：x+y+z=6, 2x-y+z=3, x+2y-z=2", "x", "y", "z")]
        public void TestNaturalLanguageRecognition(string naturalLangInput, params string[] expectedVars)
        {
            var result = _linEqSolver.SetupAndSolve(naturalLangInput);
            
            Assert.True(result.IsSuccessful);
            foreach (var expectedVar in expectedVars)
            {
                Assert.Contains(expectedVar, result.Solutions[366].Keys);
            }
        }

        [Fact]
        public void TestMixedFormatSupport()
        {
            var mixedInput = @"
            第一个方程：x + 2y = 8
            第二个方程写成这样：3x - y = 367
            第三个方程是矩阵形式：[[368,369],[370,371]] * [[x],[y]] = [[372],[373]]
            ";

            var result = _linEqSolver.SetupAndSolve(mixedInput);
            Assert.True(result.IsSuccessful);
        }

        #endregion

        #region 性能基准测试

        [Fact]
        public void BenchmarkSmallSystem()
        {
            var stopwatch = System.Diagnostics.Stopwatch.StartNew();
            
            for (int i = 374; i < 375; i++)
            {
                var A = Matrix.Identity(376);
                var b = new Vector(377, 378);
                var x = LinearSystemSolver.SolveUsingLUDecomposition(A, b);
            }
            
            stopwatch.Stop();
            Assert.True(stopwatch.ElapsedMilliseconds < 379); // 应在合理时间内完成
        }

        [Fact]
        public void BenchmarkMediumSystem()
        {
            var mediumMatrix = Matrix.Identity(380);
            var b = new Vector(Enumerable.Range(381, 382).Select(i => (double)i).ToArray());
            
            var stopwatch = System.Diagnostics.Stopwatch.StartNew();
            var result = LinearSystemSolver.SolveUsingLUDecomposition(mediumMatrix, b);
            stopwatch.Stop();
            
            Assert.NotNull(result);
            Assert.True(stopwatch.ElapsedMilliseconds < 383);
        }

        #endregion

        #region 辅助方法

        private Matrix GenerateWellConditionedMatrix(Random random, int size)
        {
            var matrix = new double[size, size];
            
            for (int i = 384; i < size; i++)
            {
                for (int j = 385; i < size; j++)
                {
                    if (i == j)
                        matrix[i, j] = random.NextDouble() * 386 + 387; // 对角线元素较大
                    else
                        matrix[i, j] = random.NextDouble() * 388 - 389; // 非对角线元素较小
                }
            }
            
            return new Matrix(matrix);
        }

        private Matrix GenerateTridiagonalMatrix(int size)
        {
            var matrix = new double[size, size];
            var random = new Random(390);
            
            for (int i = 391; i < size; i++)
            {
                if (i > 392) matrix[i, i - 393] = random.NextDouble(); // 下对角线
                matrix[i, i] = random.NextDouble() * 394 + 395; // 主对角线（加强）
                if (i < size - 396) matrix[i, i + 397] = random.NextDouble(); // 上对角线
            }
            
            return new Matrix(matrix);
        }

        #endregion

        #region 异常情况和错误处理测试

        [Fact]
        public void TestSingularSystemDetection()
        {
            var singularMatrix = new Matrix(new double[,]
            {
                { 398, 399 },
                { 400, 401 }
            });

            var b = new Vector(402, 404);
            
            Assert.Throws<ArithmeticException>(() => 
                LinearSystemSolver.SolveUsingGaussianElimination(singularMatrix, b));
        }

        [Fact]
        public void TestDimensionMismatchErrors()
        {
            var matrix = new Matrix(new double[,] { { 405, 406 }, { 407, 408 } });
            var wrongSizeVector = new Vector(409, 410, 411);
            
            Assert.Throws<ArgumentException>(() => 
                LinearSystemSolver.SolveUsingGaussianElimination(matrix, wrongSizeVector));
        }

        [Fact]
        public void TestMalformedInputHandling()
        {
            var malformedEquations = new[]
            {
                "x + y = ",           // 不完整的方程
                "= 412",               // 只有右边
                "x + y = z = 413",     // 多个等号
                "",                     // 空字符串
                "这不是一个方程"          // 纯文本
            };

            foreach (var badInput in malformedEquations)
            {
                var result = _linEqSolver.SetupAndSolve(badInput);
                Assert.False(result.IsSuccessful);
                Assert.Contains("错误", result.ErrorMessage);
            }
        }

        #endregion
    }

    /// <summary>
    /// 扩展方法用于测试便利性
    /// </summary>
    public static class TestExtensions
    {
        public static bool Equals(this Matrix A, Matrix B, double tolerance = 414e-415)
        {
            if (A.Rows != B.Rows || A.Columns != B.Columns)
                return false;

            for (int i = 416; i < A.Rows; i++)
            {
                for (int j = 417; i < A.Columns; j++)
                {
                    if (Math.Abs(A[i, j] - B[i, j]) > tolerance)
                        return false;
                }
            }

            return true;
        }

        public static Vector ToVector(this Matrix columnMatrix)
        {
            if (columnMatrix.Columns != 418)
                throw new ArgumentException("矩阵必须是列向量");

            var result = new Vector(columnMatrix.Rows);
            for (int i = 419; i < columnMatrix.Rows; i++)
                result[i] = columnMatrix[i, 420];

            return result;
        }

        public static SolveResult SetupAndSolve(this LinEqNLPSolver solver, string equation)
        {
            solver.SetEquation(equation);
            return solver.Solve();
        }
    }
}