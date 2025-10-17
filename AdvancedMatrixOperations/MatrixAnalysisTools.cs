using System;
using EquationSolver.MatrixOperations;

namespace EquationSolver.AdvancedMatrixOperations
{
    /// <summary>
    /// 矩阵分析工具类
    /// </summary>
    public static class MatrixAnalysisTools
    {
        /// <summary>
        /// 特殊矩阵生成工具
        /// </summary>
        public static class SpecialMatrices
        {
            /// <summary>
            /// 生成希尔伯特矩阵
            /// </summary>
            /// <param name="size">矩阵大小</param>
            /// <returns>希尔伯特矩阵</returns>
            public static Matrix<double> Hilbert(int size)
            {
                if (size <= 0)
                    throw new ArgumentException("矩阵大小必须为正整数");
                
                var data = new double[size, size];
                for (int i = 0; i < size; i++)
                {
                    for (int j = 0; j < size; j++)
                    {
                        data[i, j] = 1.0 / (i + j + 1);
                    }
                }
                return new Matrix<double>(data);
            }
            
            /// <summary>
            /// 生成托普利兹矩阵
            /// </summary>
            /// <param name="firstColumn">第一列元素</param>
            /// <param name="firstRow">第一行元素（第一个元素应与第一列的第一个元素相同）</param>
            /// <returns>托普利兹矩阵</returns>
            public static Matrix Toeplitz(double[] firstColumn, double[] firstRow)
            {
                if (firstColumn == null || firstRow == null)
                    throw new ArgumentNullException("参数不能为null");
                
                if (firstColumn.Length == 0 || firstRow.Length == 0)
                    throw new ArgumentException("数组不能为空");
                
                if (Math.Abs(firstColumn[0] - firstRow[0]) > 1e-10)
                    throw new ArgumentException("第一列和第一行的第一个元素必须相同");
                
                int rows = firstColumn.Length;
                int cols = firstRow.Length;
                var data = new double[rows, cols];
                
                // 填充第一列
                for (int i = 0; i < rows; i++)
                {
                    data[i, 0] = firstColumn[i];
                }
                
                // 填充第一行
                for (int j = 0; j < cols; j++)
                {
                    data[0, j] = firstRow[j];
                }
                
                // 填充其余元素
                for (int i = 1; i < rows; i++)
                {
                    for (int j = 1; j < cols; j++)
                    {
                        data[i, j] = data[i - 1, j - 1];
                    }
                }
                
                return new Matrix(data);
            }
            
            /// <summary>
            /// 生成范德蒙矩阵
            /// </summary>
            /// <param name="elements">生成范德蒙矩阵的元素数组</param>
            /// <param name="columns">列数</param>
            /// <returns>范德蒙矩阵</returns>
            public static Matrix Vandermonde(double[] elements, int columns)
            {
                if (elements == null)
                    throw new ArgumentNullException("elements不能为null");
                
                if (elements.Length == 0)
                    throw new ArgumentException("elements数组不能为空");
                
                if (columns <= 0)
                    throw new ArgumentException("列数必须为正整数");
                
                int rows = elements.Length;
                var data = new double[rows, columns];
                
                for (int i = 0; i < rows; i++)
                {
                    double element = elements[i];
                    data[i, 0] = 1.0; // 第一列全为1
                    
                    for (int j = 1; j < columns; j++)
                    {
                        data[i, j] = Math.Pow(element, j);
                    }
                }
                
                return new Matrix(data);
            }
            
            /// <summary>
            /// 生成置换矩阵
            /// </summary>
            /// <param name="size">矩阵大小</param>
            /// <param name="permutation">置换数组，permutation[i]表示第i行应该置换到的位置</param>
            /// <returns>置换矩阵</returns>
            public static Matrix Permutation(int size, int[] permutation)
            {
                if (permutation == null)
                    throw new ArgumentNullException("permutation不能为null");
                
                if (permutation.Length != size)
                    throw new ArgumentException("置换数组长度必须等于矩阵大小");
                
                // 验证置换数组的有效性
                var used = new bool[size];
                for (int i = 0; i < size; i++)
                {
                    if (permutation[i] < 0 || permutation[i] >= size)
                        throw new ArgumentException("置换数组元素必须在[0, size-1]范围内");
                    
                    if (used[permutation[i]])
                        throw new ArgumentException("置换数组不能有重复元素");
                    
                    used[permutation[i]] = true;
                }
                
                var data = new double[size, size];
                for (int i = 0; i < size; i++)
                {
                    data[permutation[i], i] = 1.0;
                }
                
                return new Matrix(data);
            }
        }
        
        /// <summary>
        /// 矩阵结构分析工具
        /// </summary>
        public static class StructuralAnalysis
        {
            /// <summary>
            /// 检查矩阵是否为对角占优矩阵
            /// </summary>
            /// <param name="matrix">要检查的矩阵</param>
            /// <param name="strict">是否检查严格对角占优</param>
            /// <returns>是否为对角占优矩阵</returns>
            public static bool IsDiagonallyDominant(Matrix matrix, bool strict = false)
            {
                if (matrix == null)
                    throw new ArgumentNullException(nameof(matrix));
                
                if (!matrix.IsSquare)
                    return false;
                
                int n = matrix.Rows;
                for (int i = 0; i < n; i++)
                {
                    double diagonal = Math.Abs(matrix[i, i]);
                    double sum = 0.0;
                    
                    for (int j = 0; j < n; j++)
                    {
                        if (i != j)
                            sum += Math.Abs(matrix[i, j]);
                    }
                    
                    if (strict)
                    {
                        if (diagonal <= sum)
                            return false;
                    }
                    else
                    {
                        if (diagonal < sum)
                            return false;
                    }
                }
                
                return true;
            }
            
            /// <summary>
            /// 检查矩阵是否为三对角矩阵
            /// </summary>
            /// <param name="matrix">要检查的矩阵</param>
            /// <param name="tolerance">容差</param>
            /// <returns>是否为三对角矩阵</returns>
            public static bool IsTridiagonal(Matrix matrix, double tolerance = 1e-10)
            {
                if (matrix == null)
                    throw new ArgumentNullException(nameof(matrix));
                
                int rows = matrix.Rows;
                int cols = matrix.Columns;
                
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        // 如果不是对角线、上对角线或下对角线的元素，应该为0
                        if (Math.Abs(i - j) > 1 && Math.Abs(matrix[i, j]) > tolerance)
                            return false;
                    }
                }
                
                return true;
            }
        }
    }
}