using System;
using System.Collections.Generic;
using System.Linq;
using EquationSolver.LinearAlgebra;
using EquationSolver.MatrixOperations;

namespace EquationSolver.AdvancedMatrixOperations
{
    /// <summary>
    /// 奇异值分解(SVD)求解器
    /// </summary>
    public class SvdSolver
    {
        private readonly Matrix<double> _matrix;
        private double _tolerance = 1e-10;
        private int _maxIterations = 1000;
        private bool _computeThinSvd = true;

        public SvdSolver(Matrix<double> matrix)
        {
            _matrix = matrix ?? throw new ArgumentNullException(nameof(matrix));
        }

        /// <summary>
        /// 设置收敛容差
        /// </summary>
        public SvdSolver WithTolerance(double tolerance)
        {
            _tolerance = tolerance;
            return this;
        }

        /// <summary>
        /// 设置最大迭代次数
        /// </summary>
        public SvdSolver WithMaxIterations(int maxIterations)
        {
            _maxIterations = maxIterations;
            return this;
        }

        /// <summary>
        /// 设置是否计算精简SVD
        /// </summary>
        public SvdSolver WithThinSvd(bool thinSvd)
        {
            _computeThinSvd = thinSvd;
            return this;
        }

        /// <summary>
        /// 执行奇异值分解
        /// </summary>
        public SvdResult Compute()
        {
            if (_matrix.Rows <= 3 || _matrix.Columns <= 3)
            {
                return ComputeSmallMatrixSvd();
            }

            return LargeScaleSvdUsingBidiagonalization();
        }

        /// <summary>
        /// 获取前k个奇异值和相应的奇异向量
        /// </summary>
        public PartialSvdResult ComputePartial(int k)
        {
            if (k <= 0)
                throw new ArgumentException("k必须为正整数");

            if (k >= Math.Min(_matrix.Rows, _matrix.Columns))
                return ConvertToPartialResult(Compute());

            return RandomizedSvd(k);
        }

        /// <summary>
        /// 计算矩阵的条件数（最大奇异值与最小奇异值的比值）
        /// </summary>
        public double ConditionNumber()
        {
            var svd = Compute();
            return svd.SingularValues.Max() / svd.SingularValues.Where(s => s > _tolerance).Min();
        }

        /// <summary>
        /// 计算矩阵的有效秩（大于容差的奇异值数量）
        /// </summary>
        public int EffectiveRank()
        {
            var svd = Compute();
            double threshold = _tolerance * Math.Max(_matrix.Rows, _matrix.Columns) * svd.SingularValues.Max();
            return svd.SingularValues.Count(s => s > threshold);
        }

        /// <summary>
        /// Moore-Penrose伪逆计算
        /// </summary>
        public Matrix Pseudoinverse()
        {
            var svd = Compute();
            var m = _matrix.Rows;
            var n = _matrix.Columns;
            var k = Math.Min(m, n);

            // Σ⁺ (Sigma plus)
            var sigmaPlusData = new double[m, n];
            for (int i = 0; i < k; i++)
            {
                if (svd.SingularValues[i] > _tolerance)
                {
                    sigmaPlusData[i, i] = 1.0 / svd.SingularValues[i];
                }
            }

            var sigmaPlus = new Matrix(sigmaPlusData);
            
            // A⁺ = VΣ⁺Uᵀ
            return svd.V.Multiply(sigmaPlus).Multiply(svd.U.Transpose());
        }

        #region 私有实现方法

        /// <summary>
        /// 小型矩阵的直接SVD计算
        /// </summary>
        private SvdResult ComputeSmallMatrixSvd()
        {
            // 对于小矩阵，可以通过特征值分解计算SVD
            var m = _matrix.Rows;
            var n = _matrix.Columns;

            if (m >= n)
            {
                // AᵀA的特征值分解
                var ata = _matrix.Transpose().Multiply(_matrix);
                var eigenSolver = new EigenvalueSolver(ata);
                var eigenPairs = eigenSolver.ComputeFullSpectrum();

                // 排序特征值（降序）
                var sortedPairs = eigenPairs.OrderByDescending(pair => pair.Eigenvalue.Magnitude).ToArray();

                var singularValues = sortedPairs.Select(pair => Math.Sqrt(Math.Max(0.0, pair.Eigenvalue.Real))).ToArray();
                var vColumns = sortedPairs.Select(pair => pair.Eigenvector).ToArray();
                var v = Matrix.FromColumns(vColumns);

                // 计算U矩阵
                var uColumns = new List<Vector>();
                for (int i = 0; i < n; i++)
                {
                    if (singularValues[i] > _tolerance)
                    {
                        var ui = _matrix.Multiply(vColumns[i]).Divide(singularValues[i]);
                        uColumns.Add(ui);
                    }
                    else
                    {
                        uColumns.Add(Vector.Zeros(m));
                    }
                }

                // 补充U的正交基
                while (uColumns.Count < m)
                {
                    uColumns.Add(CreateOrthogonalVector(uColumns.ToArray(), m));
                }

                var u = Matrix.FromColumns(uColumns.Take(m).ToArray());

                return new SvdResult(u, singularValues, v, _computeThinSvd);
            }
            else
            {
                // AAᵀ的特征值分解
                var aat = _matrix.Multiply(_matrix.Transpose());
                var eigenSolver = new EigenvalueSolver(aat);
                var eigenPairs = eigenSolver.ComputeFullSpectrum();

                var sortedPairs = eigenPairs.OrderByDescending(pair => pair.Eigenvalue.Magnitude).ToArray();
                var singularValues = sortedPairs.Select(pair => Math.Sqrt(Math.Max(0.0, pair.Eigenvalue.Real))).ToArray();
                var uColumns = sortedPairs.Select(pair => pair.Eigenvector).ToArray();
                var u = Matrix.FromColumns(uColumns);

                // 计算V矩阵
                var vColumns = new List<Vector>();
                for (int i = 0; i < m; i++)
                {
                    if (singularValues[i] > _tolerance)
                    {
                        var vi = _matrix.Transpose().Multiply(uColumns[i]).Divide(singularValues[i]);
                        vColumns.Add(vi);
                    }
                    else
                    {
                        vColumns.Add(Vector.Zeros(n));
                    }
                }

                while (vColumns.Count < n)
                {
                    vColumns.Add(CreateOrthogonalVector(vColumns.ToArray(), n));
                }

                var v = Matrix.FromColumns(vColumns.Take(n).ToArray());

                return new SvdResult(u, singularValues, v, _computeThinSvd);
            }
        }

        /// <summary>
        /// 大型矩阵的双对角化和隐式QR算法
        /// </summary>
        private SvdResult LargeScaleSvdUsingBidiagonalization()
        {
            // Golub-Reinsch算法
            var (u, b, v) = BidiagonalizeMatrix();
            var (uFinal, sigma, vFinal) = DiagonalizeBidiagonal(b, u, v);

            return new SvdResult(uFinal, sigma, vFinal, _computeThinSvd);
        }

        /// <summary>
        /// 矩阵双对角化
        /// </summary>
        private (Matrix U, Matrix B, Matrix V) BidiagonalizeMatrix()
        {
            var a = _matrix.Copy();
            int m = a.Rows;
            int n = a.Columns;
            int k = Math.Min(m, n);

            var u = Matrix.Identity(m);
            var v = Matrix.Identity(n);

            for (int i = 0; i < k; i++)
            {
                // 列方向Householder变换（消除下三角部分）
                if (i < m)
                {
                    var xCol = new Vector(m - i);
                    for (int j = i; j < m; j++)
                    {
                        xCol[j - i] = a[j, i];
                    }

                    if (xCol.Norm() > _tolerance)
                    {
                        var householderCol = HouseholderReflection(xCol);
                        ApplyLeftTransformation(a, householderCol, i, m);
                        ApplyLeftTransformation(u, householderCol, i, m);
                    }
                }

                // 行方向Householder变换（消除右上三角部分）
                if (i < n - 1)
                {
                    var xRow = new Vector(n - i - 1);
                    for (int j = i + 1; j < n; j++)
                    {
                        xRow[j - i - 1] = a[i, j];
                    }

                    if (xRow.Norm() > _tolerance)
                    {
                        var householderRow = HouseholderReflection(xRow);
                        ApplyRightTransformation(a, householderRow, i + 1, n);
                        ApplyRightTransformation(v, householderRow, i + 1, n);
                    }
                }
            }

            return (u, a, v);
        }

        /// <summary>
        /// 双对角矩阵的对角化
        /// </summary>
        private (Matrix U, double[] Sigma, Matrix V) DiagonalizeBidiagonal(Matrix b, Matrix u, Matrix v)
        {
            int m = b.Rows;
            int n = b.Columns;
            var sigma = new double[Math.Min(m, n)];

            for (int iter = 0; iter < _maxIterations; iter++)
            {
                // 检查对角线上的零元素
                for (int i = 0; i < Math.Min(m, n) - 1; i++)
                {
                    if (Math.Abs(b[i + 1, i]) < _tolerance)
                    {
                        b[i + 1, i] = 0.0;
                    }
                }

                // 查找最后一个显著的超对角线元素
                int q = Math.Min(m, n) - 1;
                while (q > 0 && Math.Abs(b[q, q - 1]) < _tolerance)
                {
                    q--;
                }

                if (q == 0)
                    break;

                // 查找第一个显著的超对角线元素
                int p = q - 1;
                while (p > 0 && Math.Abs(b[p, p - 1]) >= _tolerance)
                {
                    p--;
                }

                // 对子矩阵B[p:q, p:q]应用隐式QR步骤
                ImplicitQRStepForBidiagonal(b, u, v, p, q);
            }

            // 提取奇异值
            for (int i = 0; i < Math.Min(m, n); i++)
            {
                sigma[i] = Math.Abs(b[i, i]);
            }

            Array.Sort(sigma);
            Array.Reverse(sigma);

            return (u, sigma, v);
        }

        /// <summary>
        /// 随机化SVD算法（适用于大型稀疏矩阵）
        /// </summary>
        private PartialSvdResult RandomizedSvd(int k)
        {
            int m = _matrix.Rows;
            int n = _matrix.Columns;
            int oversampling = Math.Min(10, Math.Min(m, n) - k);

            // 步骤1: 生成随机高斯矩阵
            var random = new Random();
            var omega = new Matrix(n, k + oversampling);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < k + oversampling; j++)
                {
                    omega[i, j] = random.NextDouble() - 0.5;
                }
            }

            // 步骤2: 形成采样矩阵 Y = AΩ
            var y = _matrix.Multiply(omega);

            // 步骤3: 正交化Y得到Q
            var q = OrthogonalizeColumns(y);

            // 步骤4: 形成小矩阵 B = QᵀA
            var b = q.Transpose().Multiply(_matrix);

            // 步骤5: 对小矩阵B进行SVD
            var smallSvd = new SvdSolver(b).WithThinSvd(true).Compute();

            // 步骤6: 恢复原始空间的奇异向量
            var u = q.Multiply(smallSvd.U);
            var sigma = smallSvd.SingularValues.Take(k).ToArray();
            var v = smallSvd.V;

            return new PartialSvdResult(u, sigma, v, k);
        }

        /// <summary>
        /// 双对角矩阵的隐式QR步骤
        /// </summary>
        private void ImplicitQRStepForBidiagonal(Matrix b, Matrix u, Matrix v, int p, int q)
        {
            // 计算Wilkinson位移（取自右下角2×2块）
            int n = q - p + 1;
            if (n < 2)
                return;

            double d1 = b[p, p] * b[p, p];
            double d2 = b[q, q] * b[q, q];
            double f = b[p, p + 1] * b[p, p + 1];
            double g = b[q - 1, q] * b[q - 1, q];

            // 计算位移（Francis双重步位移的简化版本）
            double mu = ComputeWilkinsonShift(d1, d2, f, g);

            // 初始化Givens旋转
            double x = b[p, p] * b[p, p] - mu;
            double y = b[p, p] * b[p, p + 1];

            for (int j = p; j < q; j++)
            {
                // 计算Givens旋转参数
                var (c1, s1) = ComputeGivensRotation(x, y);
                
                // 从右侧应用旋转到B
                ApplyGivensRotationRight(b, j, j + 1, c1, s1);
                ApplyGivensRotationRight(v, j, j + 1, c1, s1);

                // 更新x,y
                x = b[j, j];
                y = b[j + 1, j];

                // 计算第二个Givens旋转
                var (c2, s2) = ComputeGivensRotation(x, y);
                
                // 从左侧应用旋转到B
                ApplyGivensRotationLeft(b, j, j + 1, c2, s2);
                ApplyGivensRotationLeft(u, j, j + 1, c2, s2);

                if (j < q - 1)
                {
                    x = b[j, j + 1];
                    y = b[j, j + 2];
                }
            }
        }

        #endregion

        #region 辅助计算方法

        /// <summary>
        /// 创建与给定向量组正交的新向量
        /// </summary>
        private Vector CreateOrthogonalVector(Vector[] existingVectors, int dimension)
        {
            var random = new Random();
            var candidate = new Vector(dimension);
            
            for (int attempt = 0; attempt < 10; attempt++)
            {
                // 生成随机向量
                for (int i = 0; i < dimension; i++)
                {
                    candidate[i] = random.NextDouble() - 0.5;
                }

                // Gram-Schmidt正交化
                foreach (var vec in existingVectors)
                {
                    if (vec.Norm() > _tolerance)
                    {
                        var projection = candidate.DotProduct(vec) / vec.DotProduct(vec);
                        candidate = candidate.Subtract(vec.Multiply(projection));
                    }
                }

                // 归一化
                var norm = candidate.Norm();
                if (norm > _tolerance)
                {
                    return candidate.Divide(norm);
                }
            }

            throw new InvalidOperationException("无法找到合适的正交向量");
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
        /// 左乘变换应用
        /// </summary>
        private void ApplyLeftTransformation(Matrix a, Matrix transformation, int start, int end)
        {
            int rows = end - start;
            var subMatrix = a.Submatrix(start, end - 1, 0, a.Columns - 1);
            var transformed = transformation.Multiply(subMatrix);
            
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < a.Columns; j++)
                {
                    a[start + i, j] = transformed[i, j];
                }
            }
        }

        /// <summary>
        /// 右乘变换应用
        /// </summary>
        private void ApplyRightTransformation(Matrix a, Matrix transformation, int start, int end)
        {
            int cols = end - start;
            var subMatrix = a.Submatrix(0, a.Rows - 1, start, end - 1);
            var transformed = subMatrix.Multiply(transformation);
            
            for (int i = 0; i < a.Rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    a[i, start + j] = transformed[i, j];
                }
            }
        }

        /// <summary>
        /// 列向量的正交化（Gram-Schmidt过程）
        /// </summary>
        private Matrix OrthogonalizeColumns(Matrix a)
        {
            int n = a.Columns;
            var q = a.Copy();

            for (int j = 0; j < n; j++)
            {
                var columnJ = q.GetColumn(j);
                
                // 正交化过程
                for (int k = 0; k < j; k++)
                {
                    var columnK = q.GetColumn(k);
                    var projection = columnJ.DotProduct(columnK) / columnK.DotProduct(columnK);
                    columnJ = columnJ.Subtract(columnK.Multiply(projection));
                }

                // 归一化
                var norm = columnJ.Norm();
                if (norm > _tolerance)
                {
                    columnJ = columnJ.Divide(norm);
                }

                // 设置列
                q.SetColumn(j, columnJ);
            }

            return q;
        }

        /// <summary>
        /// 计算Wilkinson位移
        /// </summary>
        private double ComputeWilkinsonShift(double d1, double d2, double f, double g)
        {
            double delta = (d1 - d2) / 2.0;
            return d2 - (f * g) / (delta + Math.Sign(delta) * Math.Sqrt(delta * delta + f * g));
        }

        /// <summary>
        /// 计算Givens旋转参数
        /// </summary>
        private (double c, double s) ComputeGivensRotation(double x, double y)
        {
            if (Math.Abs(y) < _tolerance)
                return (1.0, 0.0);

            if (Math.Abs(y) > Math.Abs(x))
            {
                double tau = -x / y;
                double s = 1.0 / Math.Sqrt(1.0 + tau * tau);
                double c = s * tau;
                return (c, s);
            }
            else
            {
                double tau = -y / x;
                double c = 1.0 / Math.Sqrt(1.0 + tau * tau);
                double s = c * tau;
                return (c, s);
            }
        }

        /// <summary>
        /// 从右侧应用Givens旋转
        /// </summary>
        private void ApplyGivensRotationRight(Matrix a, int i, int j, double c, double s)
        {
            for (int k = 0; k < a.Rows; k++)
            {
                double temp1 = c * a[k, i] - s * a[k, j];
                double temp2 = s * a[k, i] + c * a[k, j];
                a[k, i] = temp1;
                a[k, j] = temp2;
            }
        }

        /// <summary>
        /// 从左侧应用Givens旋转
        /// </summary>
        private void ApplyGivensRotationLeft(Matrix a, int i, int j, double c, double s)
        {
            for (int k = 0; k < a.Columns; k++)
            {
                double temp1 = c * a[i, k] - s * a[j, k];
                double temp2 = s * a[i, k] + c * a[j, k];
                a[i, k] = temp1;
                a[j, k] = temp2;
            }
        }

        /// <summary>
        /// 转换完整SVD结果为部分SVD结果
        /// </summary>
        private PartialSvdResult ConvertToPartialResult(SvdResult fullSvd)
        {
            return new PartialSvdResult(fullSvd.U, fullSvd.SingularValues, fullSvd.V, Math.Min(fullSvd.U.Rows, fullSvd.V.Rows));
        }

        #endregion
    }

    /// <summary>
    /// SVD分解结果
    /// </summary>
    public class SvdResult
    {
        public Matrix U { get; }
        public double[] SingularValues { get; }
        public Matrix V { get; }
        public bool IsThinSvd { get; }

        public SvdResult(Matrix u, double[] singularValues, Matrix v, bool isThinSvd)
        {
            U = u;
            SingularValues = singularValues;
            V = v;
            IsThinSvd = isThinSvd;
        }

        /// <summary>
        /// 重构原始矩阵
        /// </summary>
        public Matrix Reconstruct()
        {
            int m = U.Rows;
            int n = V.Rows;
            int k = Math.Min(m, n);

            var sigma = Matrix.Zeros(m, n);
            for (int i = 0; i < Math.Min(k, SingularValues.Length); i++)
            {
                sigma[i, i] = SingularValues[i];
            }

            return U.Multiply(sigma).Multiply(V.Transpose());
        }

        /// <summary>
        /// 低秩近似
        /// </summary>
        public Matrix LowRankApproximation(int rank)
        {
            if (rank <= 0 || rank > SingularValues.Length)
                throw new ArgumentException("秩数无效");

            var uTruncated = U.Submatrix(0, U.Rows - 1, 0, rank - 1);
            var vTruncated = V.Submatrix(0, V.Rows - 1, 0, rank - 1);
            var sigmaTruncated = Matrix.Zeros(rank, rank);
            
            for (int i = 0; i < rank; i++)
            {
                sigmaTruncated[i, i] = SingularValues[i];
            }

            return uTruncated.Multiply(sigmaTruncated).Multiply(vTruncated.Transpose());
        }
    }

    /// <summary>
    /// 部分SVD分解结果
    /// </summary>
    public class PartialSvdResult
    {
        public Matrix U { get; }
        public double[] SingularValues { get; }
        public Matrix V { get; }
        public int K { get; }

        public PartialSvdResult(Matrix u, double[] singularValues, Matrix v, int k)
        {
            U = u;
            SingularValues = singularValues;
            V = v;
            K = k;
        }
    }
}