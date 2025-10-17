using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EquationSolver.MatrixOperations
{
    /// <summary>
    /// 矩阵类 - 提供类似MATLAB的高级矩阵操作
    /// </summary>
    public class Matrix<T>
    {
        private readonly T[,] _data;
        public int Rows { get; }
        public int Columns { get; }

        /// <summary>
        /// 构造指定大小的零矩阵
        /// </summary>
        public Matrix(int rows, int columns)
        {
            if (rows <= 0 || columns <= 0)
                throw new ArgumentException("矩阵行列数必须为正整数");

            Rows = rows;
            Columns = columns;
            _data = new T[rows, columns];
        }

        /// <summary>
        /// 从二维数组构造矩阵
        /// </summary>
        public Matrix(T[,] data)
        {
            if (data == null)
                throw new ArgumentNullException(nameof(data));

            Rows = data.GetLength(0);
            Columns = data.GetLength(1);
            _data = (T[,])data.Clone();
        }

        /// <summary>
        /// 从向量构造列矩阵
        /// </summary>
        public Matrix(IEnumerable<T> vector)
        {
            if (vector == null)
                throw new ArgumentNullException(nameof(vector));

            var vecArray = vector.ToArray();
            Rows = vecArray.Length;
            Columns = 1;  // 列矩阵只有一列
            _data = new T[Rows, Columns];

            for (int i = 0; i < Rows; i++)
            {
                _data[i, 0] = vecArray[i];
            }
        }

        /// <summary>
        /// 索引器
        /// </summary>
        public T this[int row, int col]
        {
            get
            {
                ValidateIndices(row, col);
                return _data[row, col];
            }
            set
            {
                ValidateIndices(row, col);
                _data[row, col] = value;
            }
        }

        /// <summary>
        /// 验证行列索引有效性
        /// </summary>
        private void ValidateIndices(int row, int col)
        {
            if (row < 0 || row >= Rows)
                throw new IndexOutOfRangeException($"行索引 {row} 超出范围 [0, {Rows - 1}]");
            if (col < 0 || col >= Columns)
                throw new IndexOutOfRangeException($"列索引 {col} 超出范围 [0, {Columns - 1}]");
        }

        #region 静态工厂方法

        /// <summary>
        /// 创建单位矩阵
        /// </summary>
        public static Matrix<double> Identity(int size)
        {
            var identity = new Matrix<double>(size, size);
            for (int i = 0; i < size; i++)
            {
                identity[i, i] = 1.0f;
            }
            return identity;
        }

        /// <summary>
        /// 创建全零矩阵
        /// </summary>
        public static Matrix<double> Zeros(int rows, int columns)
        {
            return new Matrix<double>(rows, columns);
        }

        /// <summary>
        /// 创建全一矩阵
        /// </summary>
        public static Matrix<double> Ones(int rows, int columns)
        {
            var ones = new Matrix<double>(rows, columns);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    ones[i, j] = 1.0f;
                }
            }
            return ones;
        }

        /// <summary>
        /// 创建对角矩阵
        /// </summary>
        public static Matrix<double> Diagonal(params T[] diagonalElements)
        {
            int size = diagonalElements.Length;
            var diag = new Matrix<double>(size, size);
            for (int i = 0; i < size; i++)
            {
                diag[i, i] = diagonalElements[i];
            }
            return diag;
        }

        /// <summary>
        /// 创建随机矩阵
        /// </summary>
        public static Matrix<double> Random(int rows, int columns, T minValue = 0.0, T maxValue = 1.0)
        {
            var random = new Random();
            var randMatrix = new Matrix<double>(rows, columns);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    randMatrix[i, j] = random.NextDouble() * (maxValue - minValue) + minValue;
                }
            }
            return randMatrix;
        }

        #endregion

        #region 基本矩阵操作

        /// <summary>
        /// 矩阵相加
        /// </summary>
        public static Matrix<double> operator +(Matrix<double> a, Matrix<double> b)
        {
            if (a.Rows != b.Rows || a.Columns != b.Columns)
                throw new ArgumentException("矩阵尺寸不匹配，无法相加");

            var result = new Matrix<double>(a.Rows, a.Columns);
            for (int i = 0; i < a.Rows; i++)
            {
                for (int j = 0; j < a.Columns; j++)
                {
                    result[i, j] = a[i, j] + b[i, j];
                }
            }
            return result;
        }

        /// <summary>
        /// 矩阵相减
        /// </summary>
        public static Matrix<double> operator -(Matrix<double> a, Matrix<double> b)
        {
            if (a.Rows != b.Rows || a.Columns != b.Columns)
                throw new ArgumentException("矩阵尺寸不匹配，无法相减");

            var result = new Matrix<double>(a.Rows, a.Columns);
            for (int i = 0; i < a.Rows; i++)
            {
                for (int j = 0; j < a.Columns; j++)
                {
                    result[i, j] = a[i, j] - b[i, j];
                }
            }
            return result;
        }

        /// <summary>
        /// 标量与矩阵相乘
        /// </summary>
        public static Matrix<double> operator *(T scalar, Matrix<double> Matrix<double>)
        {
            var result = new Matrix<double>(Matrix<double>.Rows, Matrix<double>.Columns);
            for (int i = 0; i < Matrix<double>.Rows; i++)
            {
                for (int j = 0; j < Matrix<double>.Columns; j++)
                {
                    result[i, j] = scalar * Matrix<double>[i, j];
                }
            }
            return result;
        }

        /// <summary>
        /// 矩阵与标量相乘
        /// </summary>
        public static Matrix<double> operator *(Matrix<double> Matrix<double>, T scalar)
        {
            return scalar * Matrix<double>;
        }

        /// <summary>
        /// 矩阵相乘
        /// </summary>
        public static Matrix<double> operator *(Matrix<double> a, Matrix<double> b)
        {
            if (a.Columns != b.Rows)
                throw new ArgumentException($"矩阵尺寸不匹配: ({a.Rows}x{a.Columns}) × ({b.Rows}x{b.Columns})");

            var result = new Matrix<double>(a.Rows, b.Columns);
            for (int i = 0; i < a.Rows; i++)
            {
                for (int j = 0; j < b.Columns; j++)
                {
                    T sum = 0.0;
                    for (int k = 0; k < a.Columns; k++)
                    {
                        sum += a[i, k] * b[k, j];
                    }
                    result[i, j] = sum;
                }
            }
            return result;
        }

        /// <summary>
        /// 矩阵转置
        /// </summary>
        public Matrix<double> Transpose()
        {
            var transposed = new Matrix<double>(Columns, Rows);
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns; j++)
                {
                    transposed[j, i] = _data[i, j];
                }
            }
            return transposed;
        }

        /// <summary>
        /// 获取子矩阵
        /// </summary>
        public Matrix<double> Submatrix(int startRow, int endRow, int startCol, int endCol)
        {
            ValidateSubmatrixBounds(startRow, endRow, startCol, endCol);

            int subRows = endRow - startRow + 1;
            int subCols = endCol - startCol + 1;
            var submatrix = new Matrix<double>(subRows, subCols);

            for (int i = 0; i < subRows; i++)
            {
                for (int j = 0; j < subCols; j++)
                {
                    submatrix[i, j] = _data[startRow + i, startCol + j];
                }
            }
            return submatrix;
        }

        /// <summary>
        /// 验证子矩阵边界
        /// </summary>
        private void ValidateSubmatrixBounds(int startRow, int endRow, int startCol, int endCol)
        {
            if (startRow < 0 || endRow >= Rows || startRow > endRow)
                throw new ArgumentException($"无效的行范围 [{startRow}, {endRow}]");
            if (startCol < 0 || endCol >= Columns || startCol > endCol)
                throw new ArgumentException($"无效的列范围 [{startCol}, {endCol}]");
        }

        #endregion

        #region 高级矩阵操作

        /// <summary>
        /// 计算矩阵行列式
        /// </summary>
        public T Determinant()
        {
            if (Rows != Columns)
                throw new InvalidOperationException("只能计算方阵的行列式");

            if (Rows == 1)
                return _data[0, 0];

            if (Rows == 2)
                return _data[0, 0] * _data[1, 1] - _data[0, 1] * _data[1, 0];

            // 使用LU分解计算行列式
            var luDecomposition = LUDecomposition();
            T det = 1.0;
            for (int i = 0; i < Rows; i++)
            {
                det *= luDecomposition.LU[i, i];
            }
            return det * luDecomposition.Sign;
        }

        /// <summary>
        /// 计算矩阵的迹
        /// </summary>
        public T Trace()
        {
            if (Rows != Columns)
                throw new InvalidOperationException("只能计算方阵的迹");

            T trace = 0.0;
            for (int i = 0; i < Rows; i++)
            {
                trace += _data[i, i];
            }
            return trace;
        }

        /// <summary>
        /// 计算矩阵的秩
        /// </summary>
        public int Rank()
        {
            // 使用奇异值分解计算矩阵的秩
            var svd = SingularValueDecomposition();
            T epsilon = 1e-15 * Math.Max(Rows, Columns) * svd.SigmaMax;

            int rank = 0;
            for (int i = 0; i < Math.Min(Rows, Columns); i++)
            {
                if (svd.SingularValues[i] > epsilon)
                    rank++;
            }
            return rank;
        }

        /// <summary>
        /// 计算矩阵范数
        /// </summary>
        public T Norm(NormType normType = NormType.Frobenius)
        {
            switch (normType)
            {
                case NormType.Frobenius:
                    return FrobeniusNorm();
                case NormType.OneNorm:
                    return OneNorm();
                case NormType.InfinityNorm:
                    return InfinityNorm();
                case NormType.TwoNorm:
                    return TwoNorm();
                default:
                    throw new ArgumentException("不支持的范数类型");
            }
        }

        /// <summary>
        /// F-范数（弗罗贝尼乌斯范数）
        /// </summary>
        private T FrobeniusNorm()
        {
            T sum = 0.0;
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns; j++)
                {
                    sum += _data[i, j] * _data[i, j];
                }
            }
            return Math.Sqrt(sum);
        }

        /// <summary>
        /// 1-范数（列和范数）
        /// </summary>
        private T OneNorm()
        {
            T maxSum = 0.0;
            for (int j = 0; j < Columns; j++)
            {
                T columnSum = 0.0;
                for (int i = 0; i < Rows; i++)
                {
                    columnSum += Math.Abs(_data[i, j]);
                }
                maxSum = Math.Max(maxSum, columnSum);
            }
            return maxSum;
        }

        /// <summary>
        /// ∞-范数（行和范数）
        /// </summary>
        private T InfinityNorm()
        {
            T maxSum = 0.0;
            for (int i = 0; i < Rows; i++)
            {
                T rowSum = 0.0;
                for (int j = 0; j < Columns; j++)
                {
                    rowSum += Math.Abs(_data[i, j]);
                }
                maxSum = Math.Max(maxSum, rowSum);
            }
            return maxSum;
        }

        /// <summary>
        /// 2-范数（谱范数）
        /// </summary>
        private T TwoNorm()
        {
            var svd = SingularValueDecomposition();
            return svd.SingularValues[0];
        }

        /// <summary>
        /// 计算矩阵的条件数
        /// </summary>
        public T ConditionNumber(NormType normType = NormType.TwoNorm)
        {
            if (Rows != Columns)
                throw new InvalidOperationException("只能计算方阵的条件数");

            var inv = Inverse();
            return Norm(normType) * inv.Norm(normType);
        }

        #endregion

        #region 矩阵分解

        /// <summary>
        /// LU分解
        /// </summary>
        public LUDecompositionResult LUDecomposition()
        {
            if (Rows != Columns)
                throw new InvalidOperationException("LU分解仅适用于方阵");

            int n = Rows;
            var lu = (T[,])_data.Clone();
            int[] pivot = new int[n];
            int sign = 1;

            for (int i = 0; i < n; i++)
            {
                pivot[i] = i;
            }

            for (int k = 0; k < n; k++)
            {
                // 寻找主元
                int p = k;
                T max = Math.Abs(lu[k, k]);
                for (int i = k + 1; i < n; i++)
                {
                    T absValue = Math.Abs(lu[i, k]);
                    if (absValue > max)
                    {
                        max = absValue;
                        p = i;
                    }
                }

                // 交换行
                if (p != k)
                {
                    for (int j = 0; j < n; j++)
                    {
                        (lu[k, j], lu[p, j]) = (lu[p, j], lu[k, j]);
                    }
                    (pivot[k], pivot[p]) = (pivot[p], pivot[k]);
                    sign = -sign;
                }

                // 消元
                for (int i = k + 1; i < n; i++)
                {
                    lu[i, k] /= lu[k, k];
                    for (int j = k + 1; j < n; j++)
                    {
                        lu[i, j] -= lu[i, k] * lu[k, j];
                    }
                }
            }

            return new LUDecompositionResult(new Matrix<double>(lu), pivot, sign);
        }

        /// <summary>
        /// Cholesky分解
        /// </summary>
        public Matrix<double> CholeskyDecomposition()
        {
            if (Rows != Columns)
                throw new InvalidOperationException("Cholesky分解仅适用于对称正定方阵");

            int n = Rows;
            var L = new Matrix<double>(n, n);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    T sum = default;
                    for (int k = 0; k < j; k++)
                    {
                        sum += L[i, k] * L[j, k];
                    }

                    if (i == j)
                    {
                        L[i, j] = Math.Sqrt(_data[i, i] - sum);
                    }
                    else
                    {
                        L[i, j] = (_data[i, j] - sum) / L[j, j];
                    }
                }
            }
            return L;
        }

        /// <summary>
        /// QR分解
        /// </summary>
        public QRDecompositionResult<T> QRDecomposition()
        {
            int m = Rows;
            int n = Columns;
            var Q = new Matrix<T>(m, m);
            var R = new Matrix<T>(m, n);

            // Householder变换实现QR分解
            var A = (T[,])_data.Clone();
            var householderVectors = new List<T[]>();

            for (int k = 0; k < Math.Min(m, n); k++)
            {
                T[] x = new T[m - k];
                for (int i = k; i < m; i++)
                {
                    x[i - k] = A[i, k];
                }

                double normX = EuclideanNorm(x);
                double alpha = -Math.Sign(normX) * normX;

                T[] v = new T[x.Length];
                Array.Copy(x, v, x.Length);
                v[0] -= (T)Convert.ChangeType(alpha, typeof(T));

                double normV = EuclideanNorm(v);
                if (normV > 1e-15)
                {
                    for (int i = 0; i < v.Length; i++)
                    {
                        v[i] = (T)Convert.ChangeType(Convert.ToDouble(v[i]) / normV, typeof(T));
                    }
                }

                // 应用Householder变换
                for (int j = k; j < n; j++)
                {
                    double dot = 0.0;
                    for (int i = k; i < m; i++)
                    {
                        dot += Convert.ToDouble(v[i - k]) * Convert.ToDouble(A[i, j]);
                    }

                    for (int i = k; i < m; i++)
                    {
                        A[i, j] = (T)Convert.ChangeType(Convert.ToDouble(A[i, j]) - 2.0 * dot * Convert.ToDouble(v[i - k]), typeof(T));
                    }
                }

                householderVectors.Add(v);
            }

            // 构建R矩阵
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    R[i, j] = (i <= j) ? A[i, j] : (T)Convert.ChangeType(0, typeof(T));
                }
            }

            // 构建Q矩阵
            Q = new Matrix<T>(m, m);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    Q[i, j] = (i == j) ? (T)Convert.ChangeType(1, typeof(T)) : (T)Convert.ChangeType(0, typeof(T));
                }
            }
            
            for (int k = Math.Min(m, n) - 1; k >= 0; k--)
            {
                var v = householderVectors[k];
                for (int j = 0; j < m; j++)
                {
                    double dot = 0.0;
                    for (int i = k; i < m; i++)
                    {
                        dot += Convert.ToDouble(v[i - k]) * Convert.ToDouble(Q[i, j]);
                    }

                    for (int i = k; i < m; i++)
                    {
                        Q[i, j] = (T)Convert.ChangeType(Convert.ToDouble(Q[i, j]) - 2.0 * dot * Convert.ToDouble(v[i - k]), typeof(T));
                    }
                }
            }

            return new QRDecompositionResult<T>(Q, R);
        }

        /// <summary>
        /// 计算向量的欧几里得范数
        /// </summary>
        private T EuclideanNorm(T[] vector)
        {
            T sum = 0.0;
            foreach (T val in vector)
            {
                sum += val * val;
            }
            return Math.Sqrt(sum);
        }

        /// <summary>
        /// 奇异值分解 (简化版本)
        /// </summary>
        public SVDResult SingularValueDecomposition()
        {
            // 这里使用简化的雅克比旋转方法来计算奇异值分解
            // 在实际应用中，应该使用更先进的算法如LAPACK或Householder双对角化方法

            var ATA = Transpose() * this;
            var eigenDecomposition = ATA.EigenDecomposition();
            
            int singularValueCount = Math.Min(Rows, Columns);
            var singularValues = new T[singularValueCount];
            
            for (int i = 0; i < singularValueCount; i++)
            {
                singularValues[i] = Math.Sqrt(Math.Max(0.0, eigenDecomposition.Eigenvalues[i]));
            }

            // 简化版本：这里不计算完整的U和V矩阵
            return new SVDResult(null, null, singularValues);
        }

        /// <summary>
        /// 特征值分解（幂迭代法）
        /// </summary>
        public EigenDecompositionResult<T> EigenDecomposition()
        {
            if (Rows != Columns)
                throw new InvalidOperationException("特征值分解仅适用于方阵");

            int n = Rows;
            var eigenvalues = new T[n];
            var eigenvectors = new Matrix<T>(n, n);

            // 使用QR算法计算特征值
            var A = (Matrix<T>)Clone();
            const int maxIterations = 1000;
            const double tolerance = 1e-15;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                var qr = A.QRDecomposition();
                A = qr.R * qr.Q;

                // 检查收敛性
                double offDiagonalSum = 0.0;
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (i != j)
                        {
                            offDiagonalSum += Math.Abs(Convert.ToDouble(A[i, j]));
                        }
                    }
                }

                if (offDiagonalSum < tolerance)
                    break;
            }

            // 提取特征值
            for (int i = 0; i < n; i++)
            {
                eigenvalues[i] = A[i, i];
            }

            return new EigenDecompositionResult<T>(eigenvalues, eigenvectors);
        }

        /// <summary>
        /// 奇异值分解 (简化版本)
        /// </summary>
        public SVDResult<T> SingularValueDecomposition()
        {
            // 这里使用简化的雅克比旋转方法来计算奇异值分解
            // 在实际应用中，应该使用更先进的算法如LAPACK或Householder双对角化方法

            var ATA = Transpose() * this;
            var eigenDecomposition = ATA.EigenDecomposition();
            
            int singularValueCount = Math.Min(Rows, Columns);
            var singularValues = new T[singularValueCount];
            
            for (int i = 0; i < singularValueCount; i++)
            {
                singularValues[i] = (T)Convert.ChangeType(Math.Sqrt(Math.Max(0.0, Convert.ToDouble(eigenDecomposition.Eigenvalues[i]))), typeof(T));
            }

            // 简化版本：这里不计算完整的U和V矩阵
            return new SVDResult<T>(null, null, singularValues);
        }

        #endregion

        #region 线性系统求解

        /// <summary>
        /// 求解线性系统 Ax = b
        /// </summary>
        public T[] Solve(T[] b)
        {
            if (b.Length != Rows)
                throw new ArgumentException("右侧向量长度与矩阵行数不匹配");

            if (Rows == Columns)
            {
                // 方阵系统，使用LU分解
                return SolveLU(b);
            }
            else if (Rows > Columns)
            {
                // 超定系统，使用最小二乘法
                return SolveLeastSquares(b);
            }
            else
            {
                // 欠定系统，使用最小范数解
                return SolveMinimumNorm(b);
            }
        }

        /// <summary>
        /// 使用LU分解求解方程组
        /// </summary>
        private T[] SolveLU(T[] b)
        {
            var lu = LUDecomposition();
            int n = Rows;
            var x = new T[n];
            var y = new T[n];

            // 前向替换求解Ly = Pb
            for (int i = 0; i < n; i++)
            {
                y[i] = b[lu.Pivot[i]];
                for (int j = 0; j < i; j++)
                {
                    y[i] -= lu.LU[i, j] * y[j];
                }
            }

            // 后向替换求解Ux = y
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = y[i];
                for (int j = i + 1; j < n; j++)
                {
                    x[i] -= lu.LU[i, j] * x[j];
                }
                x[i] /= lu.LU[i, i];
            }

            return x;
        }

        /// <summary>
        /// 最小二乘法求解超定系统
        /// </summary>
        private T[] SolveLeastSquares(T[] b)
        {
            // 求解 A^T A x = A^T b
            var AT = Transpose();
            var ATA = AT * this;
            var ATb = AT * new Matrix(b);

            return ATA.SolveLU(ATb.ToArray());
        }

        /// <summary>
        /// 求解欠定系统的最小范数解
        /// </summary>
        private T[] SolveMinimumNorm(T[] b)
        {
            // 使用伪逆: x = A^T (A A^T)^(-1) b
            var AT = Transpose();
            var AAT = this * AT;
            var AATInv = AAT.Inverse();
            var result = AT * AATInv * new Matrix(b);

            return result.ToArray();
        }

        /// <summary>
        /// 计算矩阵的逆矩阵
        /// </summary>
        public Matrix Inverse()
        {
            if (Rows != Columns)
                throw new InvalidOperationException("只有方阵才能计算逆矩阵");

            if (Math.Abs(Determinant()) < 1e-15)
                throw new InvalidOperationException("奇异矩阵无法求逆");

            int n = Rows;
            var augmented = new Matrix(n, 2 * n);

            // 构造增广矩阵 [A | I]
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    augmented[i, j] = _data[i, j];
                    augmented[i, j + n] = (i == j) ? 1.0 : 0.0;
                }
            }

            // 高斯-约旦消元
            for (int i = 0; i < n; i++)
            {
                // 选择主元
                int maxRow = i;
                for (int k = i + 1; k < n; k++)
                {
                    if (Math.Abs(augmented[k, i]) > Math.Abs(augmented[maxRow, i]))
                    {
                        maxRow = k;
                    }
                }

                // 交换行
                if (maxRow != i)
                {
                    for (int k = 0; k < 2 * n; k++)
                    {
                        (augmented[i, k], augmented[maxRow, k]) = (augmented[maxRow, k], augmented[i, k]);
                    }
                }

                // 消元
                T pivot = augmented[i, i];
                for (int k = 0; k < 2 * n; k++)
                {
                    augmented[i, k] /= pivot;
                }

                for (int k = 0; k < n; k++)
                {
                    if (k != i)
                    {
                        T factor = augmented[k, i];
                        for (int j = 0; j < 2 * n; j++)
                        {
                            augmented[k, j] -= factor * augmented[i, j];
                        }
                    }
                }
            }

            // 提取逆矩阵
            var inverse = new Matrix(n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    inverse[i, j] = augmented[i, j + n];
                }
            }

            return inverse;
        }

        #endregion

        #region 辅助方法

        /// <summary>
        /// 转换为一维数组
        /// </summary>
        public T[] ToArray()
        {
            var array = new T[Rows * Columns];
            int index = 0;
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns; j++)
                {
                    array[index++] = _data[i, j];
                }
            }
            return array;
        }

        /// <summary>
        /// 克隆矩阵
        /// </summary>
        public object Clone()
        {
            return new Matrix(_data);
        }

        /// <summary>
        /// 判断是否为对称矩阵
        /// </summary>
        public bool IsSymmetric(T tolerance = 1e-15)
        {
            if (Rows != Columns)
                return false;

            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns; j++)
                {
                    if (Math.Abs(_data[i, j] - _data[j, i]) > tolerance)
                        return false;
                }
            }
            return true;
        }

        /// <summary>
        /// 判断是否为正定矩阵
        /// </summary>
        public bool IsPositiveDefinite()
        {
            if (!IsSymmetric())
                return false;

            try
            {
                CholeskyDecomposition();
                return true;
            }
            catch
            {
                return false;
            }
        }

        /// <summary>
        /// 格式化输出矩阵
        /// </summary>
        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine($"Matrix ({Rows}x{Columns}):");

            for (int i = 0; i < Rows; i++)
            {
                sb.Append("[");
                for (int j = 0; j < Columns; j++)
                {
                    sb.Append($"{_data[i, j]:F6}");
                    if (j < Columns - 1)
                        sb.Append(", ");
                }
                sb.AppendLine("]");
            }

            return sb.ToString();
        }

        #endregion
    }

    #region 结果类型定义

    /// <summary>
    /// 范数类型枚举
    /// </summary>
    public enum NormType
    {
        Frobenius,   // F-范数
        OneNorm,     // 1-范数
        InfinityNorm, // ∞-范数
        TwoNorm      // 2-范数（谱范数）
    }

    /// <summary>
    /// LU分解结果
    /// </summary>
    public class LUDecompositionResult
    {
        public Matrix<T> LU { get; }
        public int[] Pivot { get; }
        public int Sign { get; }

        public LUDecompositionResult(Matrix<T> lu, int[] pivot, int sign)
        {
            LU = lu;
            Pivot = pivot;
            Sign = sign;
        }
    }

    /// <summary>
    /// QR分解结果
    /// </summary>
    public class QRDecompositionResult<T>
    {
        public Matrix<T> Q { get; }
        public Matrix<T> R { get; }

        public QRDecompositionResult(Matrix<T> q, Matrix<T> r)
        {
            Q = q;
            R = r;
        }
    }

    /// <summary>
    /// SVD分解结果
    /// </summary>
    public class SVDResult<T>
    {
        public Matrix<T> U { get; }
        public Matrix<T> V { get; }
        public T[] SingularValues { get; }
        public T SigmaMax => SingularValues.Max();

        public SVDResult(Matrix<T> u, Matrix<T> v, T[] singularValues)
        {
            U = u;
            V = v;
            SingularValues = singularValues;
        }
    }

    /// <summary>
    /// 特征值分解结果
    /// </summary>
    public class EigenDecompositionResult<T>
    {
        public T[] Eigenvalues { get; }
        public Matrix<T> Eigenvectors { get; }

        public EigenDecompositionResult(T[] eigenvalues, Matrix<T> eigenvectors)
        {
            Eigenvalues = eigenvalues;
            Eigenvectors = eigenvectors;
        }
    }

    #endregion
}