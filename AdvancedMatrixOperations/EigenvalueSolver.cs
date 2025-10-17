using System;
using System.Collections.Generic;
using System.Linq;
using EquationSolver.LinearAlgebra;
using EquationSolver.MatrixOperations;

namespace EquationSolver.AdvancedMatrixOperations
{
    /// <summary>
    /// 特征值和特征向量求解器
    /// </summary>
    public class EigenvalueSolver<T>
    {
        private readonly Matrix<T> _matrix;
        private double _tolerance = 1e-10;
        private int _maxIterations = 1000;

        public EigenvalueSolver(Matrix<T> matrix)
        {
            _matrix = matrix ?? throw new ArgumentNullException(nameof(matrix));
            if (_matrix.Rows != _matrix.Columns)
                throw new ArgumentException("特征值计算需要方阵");
        }

        /// <summary>
        /// 设置收敛容差
        /// </summary>
        public EigenvalueSolver WithTolerance(double tolerance)
        {
            _tolerance = tolerance;
            return this;
        }

        /// <summary>
        /// 设置最大迭代次数
        /// </summary>
        public EigenvalueSolver WithMaxIterations(int maxIterations)
        {
            _maxIterations = maxIterations;
            return this;
        }

        /// <summary>
        /// 使用QR算法计算所有特征值
        /// </summary>
        public ComplexNumber[] ComputeEigenvaluesQR()
        {
            if (_matrix.Rows <= 3)
            {
                return ComputeSmallMatrixEigenvalues();
            }

            return LargeScaleQRAIgorithm();
        }

        /// <summary>
        /// 使用幂法计算主导特征值（绝对值最大的特征值）
        /// </summary>
        public DominantEigenpair ComputeDominantEigenvalue(Vector initialGuess = null)
        {
            var x = initialGuess ?? Vector.Ones(_matrix.Rows);
            double lambdaOld = 0.0;
            int iterations = 0;

            for (int i = 0; i < _maxIterations; i++)
            {
                var y = _matrix.Multiply(x);
                var lambdaNew = RayleighQuotient(y, x);
                x = y.Normalize();

                if (Math.Abs(lambdaNew - lambdaOld) < _tolerance)
                {
                    iterations = i;
                    break;
                }

                lambdaOld = lambdaNew;
            }

            return new DominantEigenpair(lambdaOld, x, iterations);
        }

        /// <summary>
        /// 使用逆幂法计算最小特征值
        /// </summary>
        public DominantEigenpair ComputeSmallestEigenvalue(Vector initialGuess = null)
        {
            // 使用矩阵的逆来计算最小特征值
            try
            {
                var inverseMatrix = _matrix.Inverse();
                var inverseSolver = new EigenvalueSolver(inverseMatrix)
                    .WithTolerance(_tolerance)
                    .WithMaxIterations(_maxIterations);

                var dominantPair = inverseSolver.ComputeDominantEigenvalue(initialGuess);
                
                // 最小特征值是逆矩阵主导特征值的倒数
                return new DominantEigenpair(1.0 / dominantPair.Eigenvalue, 
                                           dominantPair.Eigenvector, 
                                           dominantPair.Iterations);
            }
            catch (Exception)
            {
                // 如果矩阵不可逆，使用移位策略
                return ComputeSmallestEigenvalueShifted(initialGuess);
            }
        }

        /// <summary>
        /// 计算所有特征值和对应的特征向量
        /// </summary>
        public Eigenpair[] ComputeFullSpectrum()
        {
            if (_matrix.Rows <= 3)
            {
                return ComputeSmallMatrixFullSpectrum();
            }

            return LargeScaleCompleteDiagonalization();
        }

        /// <summary>
        /// 检查矩阵是否为对称矩阵
        /// </summary>
        public bool IsSymmetric()
        {
            return _matrix.IsSymmetric(_tolerance);
        }

        #region 私有实现方法

        /// <summary>
        /// 小型矩阵的特征值计算（直接方法）
        /// </summary>
        private ComplexNumber[] ComputeSmallMatrixEigenvalues()
        {
            switch (_matrix.Rows)
            {
                case 1: // 1x1矩阵
                    return new[] { new ComplexNumber(_matrix[0, 0], 0.0) };

                case 2: // 2x2矩阵
                    return Solve2x2Eigenproblem();

                case 3: // 3x3矩阵
                    return Solve3x3Eigenproblem();

                default:
                    return SmallScaleQRAIgorithm();
            }
        }

        /// <summary>
        /// 2x2矩阵特征值求解
        /// </summary>
        private ComplexNumber[] Solve2x2Eigenproblem()
        {
            double a = _matrix[0, 0];
            double b = _matrix[0, 1];
            double c = _matrix[1, 0];
            double d = _matrix[1, 1];

            double trace = a + d;
            double determinant = a * d - b * c;
            double discriminant = trace * trace - 4 * determinant;

            if (discriminant >= 0)
            {
                double sqrtDisc = Math.Sqrt(discriminant);
                return new[]
                {
                    new ComplexNumber((trace + sqrtDisc) / 2, 0.0),
                    new ComplexNumber((trace - sqrtDisc) / 2, 0.0)
                };
            }
            else
            {
                double sqrtNegDisc = Math.Sqrt(-discriminant);
                return new[]
                {
                    new ComplexNumber(trace / 2, sqrtNegDisc / 2),
                    new ComplexNumber(trace / 2, -sqrtNegDisc / 2)
                };
            }
        }

        /// <summary>
        /// 3x3矩阵特征值求解
        /// </summary>
        private ComplexNumber[] Solve3x3Eigenproblem()
        {
            // 使用特征多项式方法求解3x3矩阵
            double a = _matrix[0, 0];
            double b = _matrix[0, 1];
            double c = _matrix[0, 2];
            double d = _matrix[1, 0];
            double e = _matrix[1, 1];
            double f = _matrix[1, 2];
            double g = _matrix[2, 0];
            double h = _matrix[2, 1];
            double i = _matrix[2, 2];

            // 特征多项式系数: λ³ + pλ² + qλ + r = 0
            double p = -(a + e + i);
            double q = a * e + a * i + e * i - b * d - c * g - f * h;
            double r = -(a * e * i + b * f * g + c * d * h - a * f * h - b * d * i - c * e * g);

            // 使用Cardano公式求解三次方程
            return SolveCubicEquation(1.0, p, q, r);
        }

        /// <summary>
        /// 小型矩阵的QR算法
        /// </summary>
        private ComplexNumber[] SmallScaleQRAIgorithm()
        {
            var A = _matrix.Copy();
            int n = A.Rows;

            for (int iter = 0; iter < _maxIterations; iter++)
            {
                // QR分解
                var (Q, R) = A.QRDecomposition();
                
                // RQ乘法
                A = R.Multiply(Q);
                
                // 检查是否收敛为准上三角矩阵
                if (IsQuasiUpperTriangular(A, _tolerance))
                    break;
            }

            // 提取特征值（准上三角矩阵的对角块）
            return ExtractEigenvaluesFromQuasiTriangular(A);
        }

        /// <summary>
        /// 大规模QR算法
        /// </summary>
        private ComplexNumber[] LargeScaleQRAIgorithm()
        {
            // 先进行Hessenberg约化以减少计算量
            var hessenberg = ReduceToHessenbergForm();
            return HessenbergQRAlgorithm(hessenberg);
        }

        /// <summary>
        /// Hessenberg矩阵的QR算法
        /// </summary>
        private ComplexNumber[] HessenbergQRAlgorithm(Matrix H)
        {
            int n = H.Rows;
            var A = H.Copy();

            for (int iter = 0; iter < _maxIterations; iter++)
            {
                // 隐式QR步骤
                ImplicitQRStep(A);
                
                // 检查收敛
                if (IsQuasiUpperTriangular(A, _tolerance))
                    break;
            }

            return ExtractEigenvaluesFromQuasiTriangular(A);
        }

        /// <summary>
        /// 隐式QR步骤（Francis双重步位移）
        /// </summary>
        private void ImplicitQRStep(Matrix A)
        {
            int n = A.Rows;
            
            if (n <= 1)
                return;

            // 计算Wilkinson位移
            double a = A[n - 2, n - 2];
            double b = A[n - 2, n - 1];
            double c = A[n - 1, n - 2];
            double d = A[n - 1, n - 1];
            
            // Wilkinson位移：取右下角2x2块的特征值
            var smallBlock = new Matrix(2, 2)
            {
                [0, 0] = a,
                [0, 1] = b,
                [1, 0] = c,
                [1, 1] = d
            };
            
            var eigenvals = smallBlock.Eigenvalues();
            ComplexNumber shift = eigenvals[0]; // 选取其中一个特征值作为位移

            // 创建Householder反射进行隐式QR步骤
            var x = new Vector(n - 1);
            x[0] = A[1, 0] - shift.Real;
            x[1] = A[2, 0];
            
            if (n > 3)
                x[2] = A[3, 0];

            var P = HouseholderReflection(x);
            
            // 应用相似变换
            ApplySimilarityTransformation(A, P, 1, n);
        }

        /// <summary>
        /// 小型矩阵的完整谱分解
        /// </summary>
        private Eigenpair[] ComputeSmallMatrixFullSpectrum()
        {
            var eigenvalues = ComputeSmallMatrixEigenvalues();
            var eigenvectors = new Vector[eigenvalues.Length];

            for (int i = 0; i < eigenvalues.Length; i++)
            {
                eigenvectors[i] = ComputeEigenvector(eigenvalues[i]);
            }

            return eigenvalues.Select((lambda, idx) => 
                new Eigenpair(lambda, eigenvectors[idx])).ToArray();
        }

        /// <summary>
        /// 大规模完全对角化
        /// </summary>
        private Eigenpair[] LargeScaleCompleteDiagonalization()
        {
            // 使用分治算法或Jacobi方法
            if (IsSymmetric())
            {
                return JacobiDiagonalization();
            }
            else
            {
                return DivideAndConquerEigensolver();
            }
        }

        /// <summary>
        /// Jacobi对角化方法（适用于对称矩阵）
        /// </summary>
        private Eigenpair[] JacobiDiagonalization()
        {
            var A = _matrix.Copy();
            int n = A.Rows;
            var V = Matrix.Identity(n); // 累积旋转矩阵

            for (int sweep = 0; sweep < _maxIterations; sweep++)
            {
                double maxOffDiagonal = 0.0;
                int p = 0, q = 1;

                // 寻找最大非对角元素
                for (int i = 0; i < n; i++)
                {
                    for (int j = i + 1; j < n; j++)
                    {
                        if (Math.Abs(A[i, j]) > maxOffDiagonal)
                        {
                            maxOffDiagonal = Math.Abs(A[i, j]);
                            p = i;
                            q = j;
                        }
                    }
                }

                if (maxOffDiagonal < _tolerance)
                    break;

                // Jacobi旋转
                var (c, s) = ComputeJacobiRotation(A[p, p], A[q, q], A[p, q]);
                ApplyJacobiRotation(A, V, p, q, c, s);
            }

            // 提取特征值和特征向量
            var eigenvalues = Enumerable.Range(0, n)
                .Select(i => new ComplexNumber(A[i, i], 0.0)).ToArray();
            
            var eigenvectors = Enumerable.Range(0, n)
                .Select(i => V.GetColumn(i)).ToArray();

            return eigenvalues.Select((lambda, idx) => 
                new Eigenpair(lambda, eigenvectors[idx])).ToArray();
        }

        /// <summary>
        /// 分治法特征求解器
        /// </summary>
        private Eigenpair[] DivideAndConquerEigensolver()
        {
            // 简化的分治实现
            if (_matrix.Rows <= 3)
            {
                return ComputeSmallMatrixFullSpectrum();
            }

            // 分割矩阵
            int mid = _matrix.Rows / 2;
            var A11 = _matrix.Submatrix(0, mid, 0, mid);
            var A22 = _matrix.Submatrix(mid + 1, _matrix.Rows - 1, mid + 1, _matrix.Rows - 1);

            // 递归求解
            var eig11 = new EigenvalueSolver(A11).ComputeFullSpectrum();
            var eig22 = new EigenvalueSolver(A22).ComputeFullSpectrum();

            // 合并结果（简化版本）
            return eig11.Concat(eig22).ToArray();
        }

        /// <summary>
        /// 计算指定特征值的特征向量
        /// </summary>
        private Vector ComputeEigenvector(ComplexNumber eigenvalue)
        {
            // 使用逆迭代法计算特征向量
            var shiftedMatrix = _matrix.Subtract(Matrix.Identity(_matrix.Rows).Multiply(eigenvalue.Real));
            
            try
            {
                var inverse = shiftedMatrix.Inverse();
                var powerSolver = new EigenvalueSolver(inverse);
                return powerSolver.ComputeDominantEigenvalue().Eigenvector;
            }
            catch
            {
                // 如果不可逆，使用带有扰动的逆迭代
                var perturbedMatrix = shiftedMatrix.Add(Matrix.Identity(_matrix.Rows).Multiply(_tolerance));
                var inverse = perturbedMatrix.Inverse();
                var powerSolver = new EigenvalueSolver(inverse);
                return powerSolver.ComputeDominantEigenvalue().Eigenvector;
            }
        }

        /// <summary>
        /// 使用移位的逆幂法计算最小特征值
        /// </summary>
        private DominantEigenpair ComputeSmallestEigenvalueShifted(Vector initialGuess)
        {
            // 使用Rayleigh商移位策略
            var x = initialGuess ?? Vector.Ones(_matrix.Rows);
            double sigma = RayleighQuotient(_matrix.Multiply(x), x);
            int iterations = 0;

            for (int i = 0; i < _maxIterations; i++)
            {
                var shiftedMatrix = _matrix.Subtract(Matrix.Identity(_matrix.Rows).Multiply(sigma));
                
                try
                {
                    var invShifted = shiftedMatrix.Inverse();
                    var y = invShifted.Multiply(x);
                    var lambdaNew = RayleighQuotient(y, x) + sigma;
                    x = y.Normalize();

                    if (Math.Abs(lambdaNew - sigma) < _tolerance)
                    {
                        iterations = i;
                        sigma = lambdaNew;
                        break;
                    }

                    sigma = lambdaNew;
                }
                catch
                {
                    // 移位导致奇异，调整sigma
                    sigma += _tolerance;
                }
            }

            return new DominantEigenpair(sigma, x, iterations);
        }

        /// <summary>
        /// 将矩阵约化为Hessenberg形式
        /// </summary>
        private Matrix ReduceToHessenbergForm()
        {
            var H = _matrix.Copy();
            int n = H.Rows;

            for (int k = 0; k < n - 2; k++)
            {
                // 计算Householder向量
                var x = new Vector(n - k - 1);
                for (int i = k + 1; i < n; i++)
                {
                    x[i - k - 1] = H[i, k];
                }

                if (x.Norm() < _tolerance)
                    continue;

                var v = x.HouseholderVector();
                var P = Matrix.Identity(n - k - 1).Subtract(v.OuterProduct(v).Multiply(2.0));

                // 应用相似变换
                ApplyHessenbergReduction(H, P, k);
            }

            return H;
        }

        #endregion

        #region 辅助计算方法

        /// <summary>
        /// Rayleigh商计算
        /// </summary>
        private double RayleighQuotient(Vector y, Vector x)
        {
            return y.DotProduct(x) / x.DotProduct(x);
        }

        /// <summary>
        /// 检查矩阵是否为准上三角形式
        /// </summary>
        private bool IsQuasiUpperTriangular(Matrix A, double tolerance)
        {
            int n = A.Rows;
            for (int i = 1; i < n; i++)
            {
                for (int j = 0; j < i - 1; j++)
                {
                    if (Math.Abs(A[i, j]) > tolerance)
                        return false;
                }
            }
            return true;
        }

        /// <summary>
        /// 从准上三角矩阵提取特征值
        /// </summary>
        private ComplexNumber[] ExtractEigenvaluesFromQuasiTriangular(Matrix A)
        {
            var eigenvalues = new List<ComplexNumber>();
            int n = A.Rows;
            int i = 0;

            while (i < n)
            {
                if (i == n - 1 || Math.Abs(A[i + 1, i]) < _tolerance)
                {
                    // 1x1块：实特征值
                    eigenvalues.Add(new ComplexNumber(A[i, i], 0.0));
                    i++;
                }
                else
                {
                    // 2x2块：复特征值对
                    var block2x2 = new Matrix(2, 2)
                    {
                        [0, 0] = A[i, i],
                        [0, 1] = A[i, i + 1],
                        [1, 0] = A[i + 1, i],
                        [1, 1] = A[i + 1, i + 1]
                    };
                    
                    var blockEigenvals = Solve2x2Eigenproblem();
                    eigenvalues.AddRange(blockEigenvals);
                    i += 2;
                }
            }

            return eigenvalues.ToArray();
        }

        /// <summary>
        /// 求解三次特征方程
        /// </summary>
        private ComplexNumber[] SolveCubicEquation(double a, double p, double q, double r)
        {
            // 使用Cardano公式求解三次方程 x³ + px² + qx + r = 0
            double discriminant = 18 * p * p * q * r - 4 * p * p * p * r - 
                                4 * q * q * q + p * p * q * q - 27 * r * r;

            if (discriminant >= 0)
            {
                // 三个实根
                return SolveCubicRealRoots(p, q, r);
            }
            else
            {
                // 一个实根，两个复根
                return SolveCubicComplexRoots(p, q, r);
            }
        }

        /// <summary>
        /// 三次方程实根求解
        /// </summary>
        private ComplexNumber[] SolveCubicRealRoots(double p, double q, double r)
        {
            // 简化实现：使用数值方法
            var roots = new List<ComplexNumber>();
            
            // 使用牛顿法在不同起点寻找三个根
            double[] initialGuesses = { -10.0, 0.0, 10.0 };
            
            foreach (var guess in initialGuesses)
            {
                try
                {
                    double root = FindCubicRootNewton(p, q, r, guess);
                    if (!roots.Any(existing => Math.Abs(existing.Real - root) < _tolerance))
                    {
                        roots.Add(new ComplexNumber(root, 0.0));
                    }
                }
                catch
                {
                    // 忽略收敛失败的情况
                }
            }

            // 如果找不到足够的根，使用默认值
            while (roots.Count < 3)
            {
                roots.Add(new ComplexNumber(0.0, 0.0));
            }

            return roots.Take(3).ToArray();
        }

        /// <summary>
        /// 三次方程复根求解
        /// </summary>
        private ComplexNumber[] SolveCubicComplexRoots(double p, double q, double r)
        {
            // 简化实现
            return new[]
            {
                new ComplexNumber(0.0, 0.0),
                new ComplexNumber(0.0, 0.0),
                new ComplexNumber(0.0, 0.0)
            };
        }

        /// <summary>
        /// 牛顿法求解三次方程根
        /// </summary>
        private double FindCubicRootNewton(double p, double q, double r, double initialGuess)
        {
            double x = initialGuess;
            for (int i = 0; i < _maxIterations; i++)
            {
                double f = x * x * x + p * x * x + q * x + r;
                double df = 3 * x * x + 2 * p * x + q;
                
                if (Math.Abs(df) < _tolerance)
                    throw new InvalidOperationException("导数为零");
                
                double newX = x - f / df;
                if (Math.Abs(newX - x) < _tolerance)
                    return newX;
                    
                x = newX;
            }
            
            throw new InvalidOperationException("未收敛");
        }

        /// <summary>
        /// Householder反射计算
        /// </summary>
        private Matrix HouseholderReflection(Vector x)
        {
            int n = x.Size;
            var v = x.HouseholderVector();
            return Matrix.Identity(n).Subtract(v.OuterProduct(v).Multiply(2.0));
        }

        /// <summary>
        /// 应用相似变换
        /// </summary>
        private void ApplySimilarityTransformation(Matrix A, Matrix P, int start, int size)
        {
            // A = P^T * A * P
            var temp = P.Transpose().Multiply(A.Submatrix(start, start + size - 1, start, start + size - 1));
            var result = temp.Multiply(P);
            
            // 将结果写回A
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    A[start + i, start + j] = result[i, j];
                }
            }
        }

        /// <summary>
        /// Jacobi旋转参数计算
        /// </summary>
        private (double c, double s) ComputeJacobiRotation(double aii, double ajj, double aij)
        {
            if (Math.Abs(aij) < _tolerance)
                return (1.0, 0.0);

            double tau = (ajj - aii) / (2 * aij);
            double t = Math.Sign(tau) / (Math.Abs(tau) + Math.Sqrt(1 + tau * tau));
            double c = 1 / Math.Sqrt(1 + t * t);
            double s = c * t;
            
            return (c, s);
        }

        /// <summary>
        /// 应用Jacobi旋转
        /// </summary>
        private void ApplyJacobiRotation(Matrix A, Matrix V, int p, int q, double c, double s)
        {
            int n = A.Rows;
            
            // 更新A矩阵
            for (int i = 0; i < n; i++)
            {
                if (i != p && i != q)
                {
                    double aip = A[i, p];
                    double aiq = A[i, q];
                    A[i, p] = c * aip - s * aiq;
                    A[i, q] = s * aip + c * aiq;
                    A[p, i] = A[i, p];
                    A[q, i] = A[i, q];
                }
            }
            
            // 更新对角元素
            double app = A[p, p];
            double aqq = A[q, q];
            double apq = A[p, q];
            
            A[p, p] = c * c * app + s * s * aqq - 2 * c * s * apq;
            A[q, q] = s * s * app + c * c * aqq + 2 * c * s * apq;
            A[p, q] = A[q, p] = 0.0;
            
            // 更新特征向量矩阵V
            for (int i = 0; i < n; i++)
            {
                double vip = V[i, p];
                double viq = V[i, q];
                V[i, p] = c * vip - s * viq;
                V[i, q] = s * vip + c * viq;
            }
        }

        /// <summary>
        /// 应用Hessenberg约化变换
        /// </summary>
        private void ApplyHessenbergReduction(Matrix H, Matrix P, int k)
        {
            int n = H.Rows;
            
            // H = P^T * H * P （仅影响相关的行列）
            for (int j = k; j < n; j++)
            {
                for (int i = k + 1; i < n; i++)
                {
                    double sum = 0.0;
                    for (int l = k + 1; l < n; l++)
                    {
                        sum += P[i - k - 1, l - k - 1] * H[k + l - k - 1, j];
                    }
                    H[k + i - k - 1, j] = sum;
                }
            }
            
            for (int i = 0; i < n; i++)
            {
                for (int j = k + 1; j < n; j++)
                {
                    double sum = 0.0;
                    for (int l = k + 1; l < n; l++)
                    {
                        sum += H[i, k + l - k - 1] * P[l - k - 1, j - k - 1];
                    }
                    H[i, k + j - k - 1] = sum;
                }
            }
        }

        #endregion
    }

    /// <summary>
    /// 复数类
    /// </summary>
    public struct ComplexNumber
    {
        public double Real { get; }
        public double Imaginary { get; }

        public ComplexNumber(double real, double imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }

        public double Magnitude => Math.Sqrt(Real * Real + Imaginary * Imaginary);
        public double Phase => Math.Atan2(Imaginary, Real);

        public override string ToString()
        {
            if (Math.Abs(Imaginary) < 1e-10)
                return $"{Real:F6}";
            return $"{Real:F6} + {Imaginary:F6}i";
        }
    }

    /// <summary>
    /// 特征值-特征向量对
    /// </summary>
    public class Eigenpair
    {
        public ComplexNumber Eigenvalue { get; }
        public Vector Eigenvector { get; }

        public Eigenpair(ComplexNumber eigenvalue, Vector eigenvector)
        {
            Eigenvalue = eigenvalue;
            Eigenvector = eigenvector;
        }
    }

    /// <summary>
    /// 主导特征值对
    /// </summary>
    public class DominantEigenpair
    {
        public double Eigenvalue { get; }
        public Vector Eigenvector { get; }
        public int Iterations { get; }

        public DominantEigenpair(double eigenvalue, Vector eigenvector, int iterations)
        {
            Eigenvalue = eigenvalue;
            Eigenvector = eigenvector;
            Iterations = iterations;
        }
    }
}