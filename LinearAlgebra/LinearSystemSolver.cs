using System;
using System.Collections.Generic;
using System.Linq;
using EquationSolver.MatrixOperations;

namespace EquationSolver.LinearAlgebra
{
    /// <summary>
    /// 线性方程组求解器 - 支持多种高效的求解算法
    /// </summary>
    public static class LinearSystemSolver
    {
        /// <summary>
        /// 使用高斯消元法求解线性方程组 Ax = b
        /// </summary>
        /// <param name="A">系数矩阵</param>
        /// <param name="b">右侧向量</param>
        /// <returns>解向量 x</returns>
        public static Vector SolveUsingGaussianElimination(Matrix A, Vector b)
        {
            ValidateInput(A, b);
            
            // 增广矩阵 [A|b]
            var augmentedMatrix = AugmentMatrix(A, b);
            
            // 前向消去过程
            ForwardElimination(augmentedMatrix);
            
            // 回代过程
            return BackSubstitution(augmentedMatrix);
        }

        /// <summary>
        /// 使用部分主元高斯消元法（提高数值稳定性）
        /// </summary>
        public static Vector SolveUsingPartialPivoting(Matrix A, Vector b)
        {
            ValidateInput(A, b);
            
            var n = A.Rows;
            var augMat = AugmentMatrix(A, b);
            
            for (int col = 0; col < n - 1; col++)
            {
                // 寻找列主元
                int pivotRow = FindPivotRow(augMat, col);
                
                if (pivotRow != col)
                {
                    SwapRows(augMat, col, pivotRow);
                }
                
                // 消去下方元素
                EliminateColumn(augMat, col);
            }
            
            return BackSubstitution(augMat);
        }

        /// <summary>
        /// 使用LU分解求解线性方程组
        /// </summary>
        public static Vector SolveUsingLUDecomposition(Matrix A, Vector b)
        {
            ValidateInput(A, b);
            
            var luResult = A.DecomposeLU();
            var L = luResult.Item1;
            var U = luResult.Item2;
            var permutation = luResult.Item3;
            
            // 应用置换矩阵到b向量
            var permutedB = ApplyPermutation(b, permutation);
            
            // Ly = Pb （前向替换）
            var y = ForwardSubstitution(L, permutedB);
            
            // Ux = y （后向替换）
            return BackwardSubstitution(U, y);
        }

        /// <summary>
        /// 使用乔莱斯基分解求解对称正定方程组
        /// </summary>
        public static Vector SolveUsingCholesky(Matrix<T> A, Vector b)
        {
            ValidateInput(A, b);
            
            if (!A.IsSquare)
                throw new ArgumentException("矩阵必须是方阵");
            
            if (!A.IsSymmetric())
                throw new ArgumentException("矩阵必须是对称的");
            
            var L = A.CholeskyDecomposition();
            
            // Ly = b （前向替换）
            var y = ForwardSubstitutionLowerTriangular(L, b);
            
            // Lᵀx = y （后向替换）
            return BackwardSubstitutionUpperTriangular(L.Transpose(), y);
        }

        /// <summary>
        /// 使用雅可比迭代法求解大型稀疏方程组
        /// </summary>
        public static Vector SolveUsingJacobi(Matrix A, Vector b, double tolerance = 1e-15, int maxIterations = 1000)
        {
            ValidateInput(A, b);
            
            var n = A.Rows;
            var x = new Vector(n); // 初始猜测为零向量
            
            for (int iter = 0; iter < maxIterations; iter++)
            {
                var xNew = new Vector(n);
                
                for (int i = 0; i < n; i++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < n; j++)
                    {
                        if (j != i)
                            sum += A[i, j] * x[j];
                    }
                    xNew[i] = (b[i] - sum) / A[i, i];
                }
                
                // 检查收敛性
                if ((xNew - x).Norm(NormType.Two) < tolerance)
                    return xNew;
                
                x = xNew;
            }
            
            throw new Exception($"雅可比迭代未能在 {maxIterations} 次内收敛");
        }

        /// <summary>
        /// 使用高斯-赛德尔迭代法（比雅可比更快收敛）
        /// </summary>
        public static Vector SolveUsingGaussSeidel(Matrix A, Vector b, double tolerance = 1e-15, int maxIterations = 800)
        {
            ValidateInput(A, b);
            
            var n = A.Rows;
            var x = new Vector(n);
            
            for (int iter = 0; iter < maxIterations; iter++)
            {
                var xOld = x.Copy();
                
                for (int i = 0; i < n; i++)
                {
                    double sum1 = 0.0;
                    for (int j = 0; j < i; j++)
                        sum1 += A[i, j] * x[j];
                        
                    double sum2 = 0.0;
                    for (int j = i + 1; j < n; j++)
                        sum2 += A[i, j] * xOld[j];
                    
                    x[i] = (b[i] - sum1 - sum2) / A[i, i];
                }
                
                if ((x - xOld).Norm(NormType.Two) < tolerance)
                    return x;
            }
            
            throw new Exception($"高斯-赛德尔迭代未能在 {maxIterations} 次内收敛");
        }

        /// <summary>
        /// 求解超定方程组的最小二乘解（Ax ≈ b）
        /// </summary>
        public static Vector SolveLeastSquares(Matrix A, Vector b)
        {
            // AᵀAx = Aᵀb
            var ATA = A.Transpose().Multiply(A);
            var ATb = A.Transpose().MultiplyVector(b.ToMatrix()).ToVector();
            
            return SolveUsingLUDecomposition(ATA, ATb);
        }

        /// <summary>
        /// 求解欠定方程组的最小范数解
        /// </summary>
        public static Vector SolveMinimumNorm(Matrix A, Vector b)
        {
            // x = Aᵀ(AAᵀ)^{-1}b
            var AAT = A.Multiply(A.Transpose());
            var pseudoInvTimesB = SolveUsingLUDecomposition(AAT, b);
            return A.Transpose().MultiplyVector(pseudoInvTimesB.ToMatrix()).ToVector();
        }

        #region 私有辅助方法

        private static void ValidateInput(Matrix A, Vector b)
        {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (!A.IsSquare) throw new ArgumentException("系数矩阵必须是方阵");
            if (A.Rows != b.Size) throw new ArgumentException("矩阵行数与向量维度不匹配");
        }

        private static Matrix AugmentMatrix(Matrix A, Vector b)
        {
            var augmentedData = new double[A.Rows, A.Columns + 1];
            
            for (int i = 0; i < A.Rows; i++)
            {
                for (int j = 0; j < A.Columns; j++)
                    augmentedData[i, j] = A[i, j];
                
                augmentedData[i, A.Columns] = b[i];
            }
            
            return new Matrix(augmentedData);
        }

        private static void ForwardElimination(Matrix augmentedMatrix)
        {
            var rows = augmentedMatrix.Rows;
            var cols = augmentedMatrix.Columns - 1; // 排除最后一列（b向量）
            
            for (int pivotRow = 0; pivotRow < rows - 1; pivotRow++)
            {
                // 归一化主元行
                double pivotValue = augmentedMatrix[pivotRow, pivotRow];
                if (Math.Abs(pivotValue) < 1e-16)
                    throw new ArithmeticException("矩阵是奇异的或近似奇异的");
                
                // 消去下方行的对应列
                for (int rowBelow = pivotRow + 1; rowBelow < rows; rowBelow++)
                {
                    double factor = augmentedMatrix[rowBelow, pivotRow] / pivotValue;
                    
                    for (int col = pivotRow; col < cols + 1; col++) // 包括b向量列
                    {
                        augmentedMatrix[rowBelow, col] -= factor * augmentedMatrix[pivotRow, col];
                    }
                }
            }
        }

        private static Vector BackSubstitution(Matrix augmentedMatrix)
        {
            var n = augmentedMatrix.Rows;
            var solution = new Vector(n);
            
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = 0.0;
                for (int j = i + 1; j < n; j++)
                    sum += augmentedMatrix[i, j] * solution[j];
                
                solution[i] = (augmentedMatrix[i, n] - sum) / augmentedMatrix[i, i];
            }
            
            return solution;
        }

        private static int FindPivotRow(Matrix matrix, int column)
        {
            int maxRow = column;
            double maxVal = Math.Abs(matrix[column, column]);
            
            for (int row = column + 1; row < matrix.Rows; row++)
            {
                double absVal = Math.Abs(matrix[row, column]);
                if (absVal > maxVal)
                {
                    maxVal = absVal;
                    maxRow = row;
                }
            }
            
            return maxRow;
        }

        private static void SwapRows(Matrix matrix, int row1, int row2)
        {
            for (int col = 0; col < matrix.Columns; col++)
            {
                double temp = matrix[row1, col];
                matrix[row1, col] = matrix[row2, col];
                matrix[row2, col] = temp;
            }
        }

        private static void EliminateColumn(Matrix matrix, int pivotCol)
        {
            double pivotElement = matrix[pivotCol, pivotCol];
            
            for (int row = pivotCol + 1; row < matrix.Rows; row++)
            {
                double multiplier = matrix[row, pivotCol] / pivotElement;
                
                for (int col = pivotCol; col < matrix.Columns; col++)
                {
                    matrix[row, col] -= multiplier * matrix[pivotCol, col];
                }
            }
        }

        private static Vector ApplyPermutation(Vector v, int[] permutation)
        {
            var result = new Vector(v.Size);
            for (int i = 0; i < v.Size; i++)
                result[i] = v[permutation[i]];
            return result;
        }

        private static Vector ForwardSubstitution(Matrix L, Vector b)
        {
            var n = L.Rows;
            var y = new Vector(n);
            
            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < i; j++)
                    sum += L[i, j] * y[j];
                
                y[i] = (b[i] - sum) / L[i, i];
            }
            
            return y;
        }

        private static Vector BackwardSubstitution(Matrix U, Vector y)
        {
            var n = U.Rows;
            var x = new Vector(n);
            
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = 0.0;
                for (int j = i + 1; j < n; j++)
                    sum += U[i, j] * x[j];
                
                x[i] = (y[i] - sum) / U[i, i];
            }
            
            return x;
        }

        private static Vector ForwardSubstitutionLowerTriangular(Matrix L, Vector b)
        {
            return ForwardSubstitution(L, b);
        }

        private static Vector BackwardSubstitutionUpperTriangular(Matrix U, Vector y)
        {
            return BackwardSubstitution(U, y);
        }

        #endregion
    }

    /// <summary>
    /// 向量类 - 简化的一维数组包装
    /// </summary>
    public class Vector
    {
        private readonly double[] _data;

        public int Size => _data.Length;
        
        public double this[int index]
        {
            get => _data[index];
            set => _data[index] = value;
        }

        public Vector(int size)
        {
            _data = new double[size];
        }

        public Vector(params double[] elements)
        {
            _data = (double[])elements.Clone();
        }

        public Vector Copy()
        {
            return new Vector((double[])_data.Clone());
        }

        public Matrix ToMatrix(bool asColumnVector = true)
        {
            if (asColumnVector)
            {
                var matrixData = new double[_data.Length, 1];
                for (int i = 0; i < _data.Length; i++)
                    matrixData[i, 0] = _data[i];
                return new Matrix(matrixData);
            }
            else
            {
                var matrixData = new double[1, _data.Length];
                for (int i = 0; i < _data.Length; i++)
                    matrixData[0, i] = _data[i];
                return new Matrix(matrixData);
            }
        }

        public static Vector operator +(Vector a, Vector b)
        {
            if (a.Size != b.Size)
                throw new ArgumentException("向量维度不匹配");
            
            var result = new Vector(a.Size);
            for (int i = 0; i < a.Size; i++)
                result[i] = a[i] + b[i];
            return result;
        }

        public static Vector operator -(Vector a, Vector b)
        {
            if (a.Size != b.Size)
                throw new ArgumentException("向量维度不匹配");
            
            var result = new Vector(a.Size);
            for (int i = 0; i < a.Size; i++)
                result[i] = a[i] - b[i];
            return result;
        }

        public static Vector operator *(double scalar, Vector vector)
        {
            var result = new Vector(vector.Size);
            for (int i = 0; i < vector.Size; i++)
                result[i] = scalar * vector[i];
            return result;
        }

        public double DotProduct(Vector other)
        {
            if (Size != other.Size)
                throw new ArgumentException("向量维度不匹配");
            
            double result = 0.0;
            for (int i = 0; i < Size; i++)
                result += _data[i] * other._data[i];
            return result;
        }

        public double Norm(NormType normType = NormType.Two)
        {
            switch (normType)
            {
                case NormType.One:
                    return _data.Select(Math.Abs).Sum();
                case NormType.Two:
                    return Math.Sqrt(DotProduct(this));
                case NormType.Infinity:
                    return _data.Max(Math.Abs);
                default:
                    throw new ArgumentOutOfRangeException(nameof(normType));
            }
        }

        public override string ToString()
        {
            return $"Vector[{Size}] {{ {string.Join(", ", _data)} }}";
        }
    }

    /// <summary>
    /// 范数类型枚举
    /// </summary>
    public enum NormType
    {
        One,
        Two,
        Infinity
    }
}