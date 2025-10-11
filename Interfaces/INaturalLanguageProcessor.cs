using System.Collections.Generic;

namespace EquationSolver.Interfaces
{
    /// <summary>
    /// 自然语言处理器接口 - 将自然语言转换为数学表达式
    /// </summary>
    public interface INaturalLanguageProcessor
    {
        /// <summary>
        /// 预处理文本，移除无关信息
        /// </summary>
        string PreprocessText(string inputText);
        
        /// <summary>
        /// 识别方程关键词和模式
        /// </summary>
        NaturalLanguagePattern RecognizeEquationPattern(string processedText);
        
        /// <summary>
        /// 提取变量名称和约束条件
        /// </summary>
        VariableExtractionResult ExtractVariablesAndConstraints(string text);
        
        /// <summary>
        /// 转换自然语言方程为标准数学表达式
        /// </summary>
        StandardizedEquation ConvertToMathematicalNotation(NaturalLanguagePattern pattern);
        
        /// <summary>
        /// 验证自然语言输入的完整性
        /// </summary>
        ValidationResult ValidateInputCompleteness(string inputText);
    }

    /// <summary>
    /// 自然语言模式识别结果
    /// </summary>
    public struct NaturalLanguagePattern
    {
        public string OriginalText { get; set; }
        public EquationType Type { get; set; }
        public List<string> Variables { get; set; }
        public List<string> MathematicalTerms { get; set; }
        public List<string> Constraints { get; set; }
        public ConfidenceScore Confidence { get; set; }
    }

    /// <summary>
    /// 标准化方程结果
    /// </summary>
    public struct StandardizedEquation
    {
        public string MathematicalExpression { get; set; }
        public Dictionary<string, double> InitialValues { get; set; }
        public List<string> DomainRestrictions { get; set; }
        public SolverRecommendation RecommendedSolver { get; set; }
    }

    /// <summary>
    /// 变量提取结果
    /// </summary>
    public struct VariableExtractionResult
    {
        public Dictionary<string, VariableInfo> Variables { get; set; }
        public List<ConstraintInfo> Constraints { get; set; }
        public List<ParameterRange> ParameterRanges { get; set; }
    }

    // 枚举和辅助结构定义
    public enum EquationType
    {
        Linear,
        Quadratic,
        Polynomial,
        Exponential,
        Logarithmic,
        Trigonometric,
        Differential,
        Integral,
        Matrix,
        SystemOfEquations,
        Inequality,
        Undetermined
    }

    public struct ConfidenceScore
    {
        public double Score { get; set; }
        public string Reason { get; set; }
    }

    public struct VariableInfo
    {
        public string Name { get; set; }
        public VariableType Type { get; set; }
        public double? SuggestedValue { get; set; }
        public List<string> Synonyms { get; set; }
    }

    public enum VariableType
    {
        IndependentVariable,
        DependentVariable,
        Constant,
        Parameter,
        Coefficient
    }

    public struct ConstraintInfo
    {
        public string Condition { get; set; }
        public ConstraintType Type { get; set; }
        public string Description { get; set; }
    }

    public enum ConstraintType
    {
        Equality,
        Inequality,
        Range,
        BoundaryCondition,
        InitialCondition
    }

    public struct ParameterRange
    {
        public string ParameterName { get; set; }
        public double MinValue { get; set; }
        public double MaxValue { get; set; }
        public bool InclusiveMin { get; set; }
        public bool InclusiveMax { get; set; }
    }

    public struct SolverRecommendation
    {
        public string SolverType { get; set; }
        public string Algorithm { get; set; }
        public string Justification { get; set; }
        public double SuitabilityScore { get; set; }
    }

    public struct ValidationResult
    {
        public bool IsValid { get; set; }
        public List<string> MissingInformation { get; set; }
        public List<string> Suggestions { get; set; }
        public List<string> Warnings { get; set; }
    }
}
