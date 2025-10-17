using System;
using System.Text;
using EquationSolver.MatrixOperations;

namespace EquationSolver.AdvancedMatrixOperations
{
    /// <summary>
    /// 综合矩阵分析器
    /// </summary>
    public static class ComprehensiveMatrixAnalyzer
    {
        /// <summary>
        /// 对矩阵进行全面分析
        /// </summary>
        public static MatrixAnalysisResult AnalyzeMatrix<T>(Matrix<T> matrix)
        {
            var result = new MatrixAnalysisResult();
            result.Rows = matrix.Rows;
            result.Columns = matrix.Columns;
            result.IsSquare = matrix.Rows == matrix.Columns;
            
            // 计算行列式（仅对方阵）
            if (result.IsSquare)
            {
                try
                {
                    // 转换为double矩阵进行计算
                    var doubleMatrix = ConvertToDoubleMatrix(matrix);
                    result.Determinant = doubleMatrix.Determinant();
                }
                catch
                {
                    result.Determinant = double.NaN;
                }
            }
            
            // 计算条件数
            try
            {
                var doubleMatrix = ConvertToDoubleMatrix(matrix);
                result.ConditionNumber = doubleMatrix.ConditionNumber();
            }
            catch
            {
                result.ConditionNumber = double.NaN;
            }
            
            // 计算秩
            try
            {
                var doubleMatrix = ConvertToDoubleMatrix(matrix);
                result.Rank = doubleMatrix.Rank();
            }
            catch
            {
                result.Rank = -1; // 表示计算失败
            }
            
            // 计算范数
            try
            {
                var doubleMatrix = ConvertToDoubleMatrix(matrix);
                result.FrobeniusNorm = doubleMatrix.Norm(NormType.Frobenius);
                result.OneNorm = doubleMatrix.Norm(NormType.OneNorm);
                result.InfinityNorm = doubleMatrix.Norm(NormType.InfinityNorm);
            }
            catch
            {
                result.FrobeniusNorm = double.NaN;
                result.OneNorm = double.NaN;
                result.InfinityNorm = double.NaN;
            }
            
            // 检查矩阵特性
            if (result.IsSquare)
            {
                try
                {
                    var doubleMatrix = ConvertToDoubleMatrix(matrix);
                    result.IsSymmetric = doubleMatrix.IsSymmetric(1e-10);
                    result.IsPositiveDefinite = doubleMatrix.IsPositiveDefinite();
                }
                catch
                {
                    result.IsSymmetric = false;
                    result.IsPositiveDefinite = false;
                }
            }
            
            return result;
        }
        
        /// <summary>
        /// 将泛型矩阵转换为double矩阵
        /// </summary>
        private static Matrix ConvertToDoubleMatrix<T>(Matrix<T> matrix)
        {
            var doubleData = new double[matrix.Rows, matrix.Columns];
            for (int i = 0; i < matrix.Rows; i++)
            {
                for (int j = 0; j < matrix.Columns; j++)
                {
                    doubleData[i, j] = Convert.ToDouble(matrix[i, j]);
                }
            }
            return new Matrix(doubleData);
        }
    }
    
    /// <summary>
    /// 矩阵分析结果
    /// </summary>
    public class MatrixAnalysisResult
    {
        public int Rows { get; set; }
        public int Columns { get; set; }
        public bool IsSquare { get; set; }
        public double Determinant { get; set; }
        public double ConditionNumber { get; set; }
        public int Rank { get; set; }
        public double FrobeniusNorm { get; set; }
        public double OneNorm { get; set; }
        public double InfinityNorm { get; set; }
        public bool IsSymmetric { get; set; }
        public bool IsPositiveDefinite { get; set; }
        
        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("=== 矩阵分析结果 ===");
            sb.AppendLine($"维度: {Rows} × {Columns}");
            sb.AppendLine($"是否为方阵: {(IsSquare ? "是" : "否")}");
            
            if (IsSquare)
            {
                sb.AppendLine($"行列式: {Determinant}");
                sb.AppendLine($"条件数: {ConditionNumber}");
                sb.AppendLine($"秩: {Rank}");
                sb.AppendLine($"是否对称: {(IsSymmetric ? "是" : "否")}");
                sb.AppendLine($"是否正定: {(IsPositiveDefinite ? "是" : "否")}");
            }
            else
            {
                sb.AppendLine($"条件数: {ConditionNumber}");
                sb.AppendLine($"秩: {Rank}");
            }
            
            sb.AppendLine($"F-范数: {FrobeniusNorm}");
            sb.AppendLine($"1-范数: {OneNorm}");
            sb.AppendLine($"∞-范数: {InfinityNorm}");
            
            return sb.ToString();
        }
    }
}