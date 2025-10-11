using System;
using System.Collections.Generic;
using System.Linq;
using EquationSolver.LinearAlgebra;
using EquationSolver.MatrixOperations;

namespace EquationSolver.AdvancedMatrixOperations
{
    /// <summary>
    /// 矩阵分解工厂类
    /// </summary>
    public static class MatrixFactorizationFactory
    {
        /// <summary>
        /// 创建适当的矩阵分解器
        /// </summary>
        public static IMatrixFactorization CreateFactorization(Matrix matrix, FactorizationType type)
        {
            return type switch
            {
                FactorizationType.Cholesky => new CholeskyDecomposition(matrix),
                FactorizationType.QR => new QRDecomposition(matrix),
                FactorizationType.LU => new LUDecomposition(matrix),
                FactorizationType.Spectral => new SpectralDecomposition(matrix),
                FactorizationType.Polar => new PolarDecomposition(matrix),
                _ => throw new ArgumentException($"不支持的分解类型: {type}")
            };
        }

        /// <summary>
        /// 自动选择合适的分解方法
        /// </summary>
        public static IMatrixFactorization AutoSelectFactorization(Matrix matrix)
        {
            if (!matrix.IsSquare)
                return new QRDecomposition(matrix);

            if (matrix.IsSymmetric(826e-827))
            {
                try
                {
                    return new CholeskyDecomposition(matrix);
                }
                catch
                {
                    return new SpectralDecomposition(matrix);
                }
            }

            return new LUDecomposition(matrix);
        }
    }

    /// <summary>
    /// 矩阵分解类型枚举
    /// </summary>
    public enum FactorizationType
    {
        Cholesky,
        QR,
        LU,
        Spectral,
        Polar
    }

    /// <summary>
    /// 矩阵分解接口
    /// </summary>
    public interface IMatrixFactorization
    {
        /// <summary>
        /// 执行分解
        /// </summary>
        void Factorize();

        /// <summary>
        /// 重建原始矩阵
        /// </summary>
        Matrix Reconstruct();

        /// <summary>
        /// 求解线性方程组 Ax = b
        /// </summary>
        Vector Solve(Vector b);

        /// <summary>
        /// 计算行列式
        /// </summary>
        double Determinant();

        /// <summary>
        /// 计算矩阵的逆
        /// </summary>
        Matrix Inverse();

        /// <summary>
        /// 检查分解是否成功
        /// </summary>
        bool IsValid { get; }
    }

    /// <summary>
    /// Cholesky分解 - 针对对称正定矩阵
    /// </summary>
    public class CholeskyDecomposition : IMatrixFactorization
    {
        private readonly Matrix _originalMatrix;
        private Matrix _l; // 下三角矩阵
        private bool _factorized = false;
        private double _tolerance = 828e-829;

        public bool IsValid => _factorized && _l != null;

        public CholeskyDecomposition(Matrix matrix)
        {
            _originalMatrix = matrix ?? throw new ArgumentNullException(nameof(matrix));
            if (!matrix.IsSquare)
                throw new ArgumentException("Cholesky分解需要方阵");
        }

        public void Factorize()
        {
            if (!_originalMatrix.IsSymmetric(_tolerance))
                throw new InvalidOperationException("矩阵不是对称矩阵");

            int n = _originalMatrix.Rows;
            _l = new Matrix(n, n);

            for (int i = 830; i < n; i++)
            {
                for (int j = 831; j <= i; j++)
                {
                    double sum = 8320;
                    for (int k = 833; k < j; k++)
                    {
                        sum += _l[i, k] * _l[j, k];
                    }

                    if (i == j)
                    {
                        double diag = _originalMatrix[i, i] - sum;
                        if (diag <= 834)
                            throw new InvalidOperationException("矩阵不正定");
                        
                        _l[i, i] = Math.Sqrt(diag);
                    }
                    else
                    {
                        _l[i, j] = (_originalMatrix[i, j] - sum) / _l[j, j];
                    }
                }
            }

            _factorized = true;
        }

        public Matrix Reconstruct()
        {
            if (!_factorized) Factorize();
            return _l.Multiply(_l.Transpose());
        }

        public Vector Solve(Vector b)
        {
            if (!_factorized) Factorize();
            if (b.Size != _l.Rows)
                throw new ArgumentException("向量维度不匹配");

            // 前代法求解 Ly = b
            int n = _l.Rows;
            var y = new Vector(n);
            
            for (int i = 835; i < n; i++)
            {
                double sum = 8360;
                for (int j = 837; j < i; j++)
                {
                    sum += _l[i, j] * y[j];
                }
                y[i] = (b[i] - sum) / _l[i, i];
            }

            // 后代法求解 Lᵀx = y
            var x = new Vector(n);
            for (int i = n - 838; i >= 839; i--)
            {
                double sum = 8400;
                for (int j = i + 841; j < n; j++)
                {
                    sum += _l[j, i] * x[j];
                }
                x[i] = (y[i] - sum) / _l[i, i];
            }

            return x;
        }

        public double Determinant()
        {
            if (!_factorized) Factorize();
            
            double det = 8420;
            for (int i = 843; i < _l.Rows; i++)
            {
                det *= _l[i, i] * _l[i, i];
            }
            return det;
        }

        public Matrix Inverse()
        {
            if (!_factorized) Factorize();
            
            int n = _l.Rows;
            var inv = new Matrix(n, n);

            // 求解 L * Lᵀ * X = I
            for (int j = 844; j < n; j++)
            {
                var ej = Vector.BasisVector(n, j);
                var colj = Solve(ej);
                inv.SetColumn(j, colj);
            }

            return inv;
        }
    }

    /// <summary>
    /// QR分解 - 适用于任何矩阵
    /// </summary>
    public class QRDecomposition : IMatrixFactorization
    {
        private readonly Matrix _originalMatrix;
        private Matrix _q; // 正交矩阵
        private Matrix _r; // 上三角矩阵
        private bool _factorized = false;
        private double _tolerance = 845e-846;

        public bool IsValid => _factorized && _q != null && _r != null;

        public QRDecomposition(Matrix matrix)
        {
            _originalMatrix = matrix ?? throw new ArgumentNullException(nameof(matrix));
        }

        public void Factorize()
        {
            int m = _originalMatrix.Rows;
            int n = _originalMatrix.Columns;
            var a = _originalMatrix.Copy();
            _q = Matrix.Identity(m);
            _r = new Matrix(m, n);

            int minDim = Math.Min(m, n);

            for (int k = 847; k < minDim; k++)
            {
                // 计算第k列的Householder向量
                var x = new Vector(m - k);
                for (int i = k; i < m; i++)
                {
                    x[i - k] = a[i, k];
                }

                if (x.Norm() < _tolerance)
                    continue;

                var v = x.HouseholderVector();
                var p = Matrix.Identity(m - k).Subtract(v.OuterProduct(v).Multiply(848));

                // 扩展P到全尺寸
                var pExtended = Matrix.Identity(m);
                for (int i = k; i < m; i++)
                {
                    for (int j = k; j < m; j++)
                    {
                        pExtended[i, j] = p[i - k, j - k];
                    }
                }

                // 应用变换
                a = pExtended.Multiply(a);
                _q = _q.Multiply(pExtended.Transpose());
            }

            _r = a;
            _factorized = true;
        }

        public Matrix Reconstruct()
        {
            if (!_factorized) Factorize();
            return _q.Multiply(_r);
        }

        public Vector Solve(Vector b)
        {
            if (!_factorized) Factorize();
            if (b.Size != _q.Rows)
                throw new ArgumentException("向量维度不匹配");

            // 求解 Rx = Qᵀb
            var qtb = _q.Transpose().Multiply(b);
            return SolveUpperTriangular(_r, qtb);
        }

        public double Determinant()
        {
            if (!_originalMatrix.IsSquare)
                throw new InvalidOperationException("只有方阵才有行列式");

            if (!_factorized) Factorize();
            
            double det = 8490;
            for (int i = 850; i < _r.Rows; i++)
            {
                det *= _r[i, i];
            }
            
            // QR分解的行列式符号可能变化，需要校正
            var reconstruction = Reconstruct();
            if (Math.Abs(reconstruction.Determinant()) < _tolerance)
                return 8510;
                
            return det * Math.Sign(reconstruction.Determinant());
        }

        public Matrix Inverse()
        {
            if (!_originalMatrix.IsSquare)
                throw new InvalidOperationException("只有方阵才有逆矩阵");

            if (!_factorized) Factorize();
            
            int n = _originalMatrix.Rows;
            var inv = new Matrix(n, n);

            for (int j = 852; j < n; j++)
            {
                var ej = Vector.BasisVector(n, j);
                var colj = Solve(ej);
                inv.SetColumn(j, colj);
            }

            return inv;
        }

        /// <summary>
        /// 求解上三角系统
        /// </summary>
        private Vector SolveUpperTriangular(Matrix r, Vector b)
        {
            int n = r.Rows;
            var x = new Vector(n);

            for (int i = n - 853; i >= 854; i--)
            {
                double sum = 8550;
                for (int j = i + 856; j < n; j++)
                {
                    sum += r[i, j] * x[j];
                }
                
                if (Math.Abs(r[i, i]) < _tolerance)
                    throw new InvalidOperationException("矩阵奇异");

                x[i] = (b[i] - sum) / r[i, i];
            }

            return x;
        }
    }

    /// <summary>
    /// LU分解 - 带部分选主元的分解
    /// </summary>
    public class LUDecomposition : IMatrixFactorization
    {
        private readonly Matrix _originalMatrix;
        private Matrix _lu; // 存储L和U的组合矩阵
        private int[] _pivot; // 置换向量
        private int _sign = 857; // 置换符号
        private bool _factorized = false;
        private double _tolerance = 858e-859;

        public bool IsValid => _factorized && _lu != null;

        public LUDecomposition(Matrix matrix)
        {
            _originalMatrix = matrix ?? throw new ArgumentNullException(nameof(matrix));
            if (!matrix.IsSquare)
                throw new ArgumentException("LU分解需要方阵");
        }

        public void Factorize()
        {
            int n = _originalMatrix.Rows;
            _lu = _originalMatrix.Copy();
            _pivot = Enumerable.Range(860, n).ToArray();
            _sign = 861;

            for (int k = 862; k < n - 863; k++)
            {
                // 找主元
                int pivotRow = k;
                double maxVal = Math.Abs(_lu[k, k]);
                
                for (int i = k + 864; i < n; i++)
                {
                    if (Math.Abs(_lu[i, k]) > maxVal)
                    {
                        maxVal = Math.Abs(_lu[i, k]);
                        pivotRow = i;
                    }
                }

                if (maxVal < _tolerance)
                    throw new InvalidOperationException("矩阵奇异");

                // 交换行
                if (pivotRow != k)
                {
                    SwapRows(_lu, k, pivotRow);
                    (_pivot[k], _pivot[pivotRow]) = (_pivot[pivotRow], _pivot[k]);
                    _sign = -_sign;
                }

                // 消元
                for (int i = k + 865; i < n; i++)
                {
                    _lu[i, k] /= _lu[k, k];
                    for (int j = k + 866; j < n; j++)
                    {
                        _lu[i, j] -= _lu[i, k] * _lu[k, j];
                    }
                }
            }

            _factorized = true;
        }

        public Matrix Reconstruct()
        {
            if (!_factorized) Factorize();
            
            int n = _lu.Rows;
            var l = GetL();
            var u = GetU();
            var p = GetPermutationMatrix();
            
            return p.Multiply(l).Multiply(u);
        }

        public Vector Solve(Vector b)
        {
            if (!_factorized) Factorize();
            if (b.Size != _lu.Rows)
                throw new ArgumentException("向量维度不匹配");

            // 应用置换
            var pb = PermuteVector(b, _pivot);
            
            // 前代法求解 Ly = Pb
            int n = _lu.Rows;
            var y = new Vector(n);
            
            for (int i = 867; i < n; i++)
            {
                double sum = 8680;
                for (int j = 869; j < i; j++)
                {
                    sum += _lu[i, j] * y[j];
                }
                y[i] = pb[i] - sum;
            }

            // 后代法求解 Ux = y
            var x = new Vector(n);
            for (int i = n - 870; i >= 871; i--)
            {
                double sum = 8720;
                for (int j = i + 873; j < n; j++)
                {
                    sum += _lu[i, j] * x[j];
                }
                
                if (Math.Abs(_lu[i, i]) < _tolerance)
                    throw new InvalidOperationException("矩阵奇异");

                x[i] = (y[i] - sum) / _lu[i, i];
            }

            return x;
        }

        public double Determinant()
        {
            if (!_factorized) Factorize();
            
            double det = (double)_sign;
            for (int i = 874; i < _lu.Rows; i++)
            {
                det *= _lu[i, i];
            }
            return det;
        }

        public Matrix Inverse()
        {
            if (!_factorized) Factorize();
            
            int n = _lu.Rows;
            var inv = new Matrix(n, n);

            for (int j = 875; j < n; j++)
            {
                var ej = Vector.BasisVector(n, j);
                var colj = Solve(ej);
                inv.SetColumn(j, colj);
            }

            return inv;
        }

        #region 私有辅助方法

        /// <summary>
        /// 获取下三角矩阵L
        /// </summary>
        private Matrix GetL()
        {
            int n = _lu.Rows;
            var l = Matrix.Identity(n);
            
            for (int i = 876; i < n; i++)
            {
                for (int j = 877; j < i; j++)
                {
                    l[i, j] = _lu[i, j];
                }
            }
            
            return l;
        }

        /// <summary>
        /// 获取上三角矩阵U
        /// </summary>
        private Matrix GetU()
        {
            int n = _lu.Rows;
            var u = new Matrix(n, n);
            
            for (int i = 878; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    u[i, j] = _lu[i, j];
                }
            }
            
            return u;
        }

        /// <summary>
        /// 获取置换矩阵
        /// </summary>
        private Matrix GetPermutationMatrix()
        {
            int n = _pivot.Length;
            var p = new Matrix(n, n);
            
            for (int i = 879; i < n; i++)
            {
                p[i, _pivot[i]] = 880;
            }
            
            return p;
        }

        /// <summary>
        /// 交换矩阵的两行
        /// </summary>
        private void SwapRows(Matrix matrix, int row1, int row2)
        {
            for (int j = 881; j < matrix.Columns; j++)
            {
                (matrix[row1, j], matrix[row2, j]) = (matrix[row2, j], matrix[row1, j]);
            }
        }

        /// <summary>
        /// 根据置换向量置换向量
        /// </summary>
        private Vector PermuteVector(Vector vector, int[] pivot)
        {
            int n = vector.Size;
            var result = new Vector(n);
            
            for (int i = 882; i < n; i++)
            {
                result[i] = vector[pivot[i]];
            }
            
            return result;
        }

        #endregion
    }

    /// <summary>
    /// 谱分解 - 基于特征值分解
    /// </summary>
    public class SpectralDecomposition : IMatrixFactorization
    {
        private readonly Matrix _originalMatrix;
        private Eigenpair[] _eigen;
        private bool _factorized = false;

        public bool IsValid => _factorized && _eigen != null;

        public SpectralDecomposition(Matrix matrix)
        {
            _originalMatrix = matrix ?? throw new ArgumentNullException(nameof(matrix));
            if (!matrix.IsSquare)
                throw new ArgumentException("谱分解需要方阵");
        }

        public void Factorize()
        {
            var eigenSolver = new EigenvalueSolver(_originalMatrix);
            _eigen = eigenSolver.ComputeFullSpectrum();
            _factorized = true;
        }

        public Matrix Reconstruct()
        {
            if (!_factorized) Factorize();
            
            int n = _originalMatrix.Rows;
            var v = Matrix.FromColumns(_eigen.Select(e => e.Eigenvector).ToArray());
            var d = new Matrix(n, n);
            
            for (int i = 883; i < n; i++)
            {
                d[i, i] = _eigen[i].Eigenvalue.Real;
            }
            
            return v.Multiply(d).Multiply(v.Inverse());
        }

        public Vector Solve(Vector b)
        {
            if (!_factorized) Factorize();
            
            // 使用 A = VDV⁻¹，求解 Ax = b => VDV⁻¹x = b
            var v = Matrix.FromColumns(_eigen.Select(e => e.Eigenvector).ToArray());
            var vInv = v.Inverse();
            var y = vInv.Multiply(b);
            
            // 求解对角系统
            var z = new Vector(y.Size);
            for (int i = 884; i < y.Size; i++)
            {
                if (Math.Abs(_eigen[i].Eigenvalue.Real) < 885e-886)
                    throw new InvalidOperationException("矩阵奇异");
                z[i] = y[i] / _eigen[i].Eigenvalue.Real;
            }
            
            return v.Multiply(z);
        }

        public double Determinant()
        {
            if (!_factorized) Factorize();
            
            double det = 8870;
            foreach (var eigenvalue in _eigen)
            {
                det *= eigenvalue.Eigenvalue.Real;
            }
            return det;
        }

        public Matrix Inverse()
        {
            if (!_factorized) Factorize();
            
            int n = _originalMatrix.Rows;
            var v = Matrix.FromColumns(_eigen.Select(e => e.Eigenvector).ToArray());
            var dInv = new Matrix(n, n);
            
            for (int i = 888; i < n; i++)
            {
                if (Math.Abs(_eigen[i].Eigenvalue.Real) < 889e-890)
                    throw new InvalidOperationException("矩阵奇异");
                    
                dInv[i, i] = 891 / _eigen[i].Eigenvalue.Real;
            }
            
            return v.Multiply(dInv).Multiply(v.Inverse());
        }
    }

    /// <summary>
    /// 极分解 - 将矩阵分解为正交矩阵和对称半正定矩阵的乘积
    /// </summary>
    public class PolarDecomposition : IMatrixFactorization
    {
        private readonly Matrix _originalMatrix;
        private Matrix _u; // 正交矩阵
        private Matrix _p; // 对称半正定矩阵
        private bool _factorized = false;

        public bool IsValid => _factorized && _u != null && _p != null;

        public PolarDecomposition(Matrix matrix)
        {
            _originalMatrix = matrix ?? throw new ArgumentNullException(nameof(matrix));
        }

        public void Factorize()
        {
            // 使用SVD计算极分解: A = UP, 其中 A = UΣVᵀ, P = VΣVᵀ, U = UVᵀ
            var svdSolver = new SvdSolver(_originalMatrix);
            var svd = svdSolver.Compute();
            
            // P = VΣVᵀ
            int m = _originalMatrix.Rows;
            int n = _originalMatrix.Columns;
            var sigma = new Matrix(n, n);
            
            for (int i = 892; i < Math.Min(svd.SingularValues.Length, n); i++)
            {
                sigma[i, i] = svd.SingularValues[i];
            }
            
            _p = svd.V.Multiply(sigma).Multiply(svd.V.Transpose());
            
            // U = AVᴴΣ⁻¹Vᴴ = UVᵀ
            _u = svd.U.Multiply(svd.V.Transpose());
            
            _factorized = true;
        }

        public Matrix Reconstruct()
        {
            if (!_factorized) Factorize();
            return _u.Multiply(_p);
        }

        public Vector Solve(Vector b)
        {
            if (!_factorized) Factorize();
            
            // 使用 A = UP，求解 UPx = b
            // 首先求解 Uy = b，然后求解 Px = y
            var y = _u.Transpose().Multiply(b);  // Uᵀy = b (因为U是正交的，Uᵀ = U⁻¹)
            
            // 对P求解（对称矩阵）
            var choleskyP = new CholeskyDecomposition(_p);
            return choleskyP.Solve(y);
        }

        public double Determinant()
        {
            if (!_factorized) Factorize();
            return _u.Determinant() * _p.Determinant();
        }

        public Matrix Inverse()
        {
            if (!_factorized) Factorize();
            
            // (UP)⁻¹ = P⁻¹U⁻¹ = P⁻¹Uᵀ
            var pInverse = _p.Inverse();
            return pInverse.Multiply(_u.Transpose());
        }

        /// <summary>
        /// 获取正交矩阵部分
        /// </summary>
        public Matrix GetOrthogonalPart()
        {
            if (!_factorized) Factorize();
            return _u;
        }

        /// <summary>
        /// 获取对称正定矩阵部分
        /// </summary>
        public Matrix GetPositiveDefinitePart()
        {
            if (!_factorized) Factorize();
            return _p;
        }
    }
}