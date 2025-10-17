using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EquationSolver.MatrixOperations
{
    /// <summary>
    /// 向量类 - 提供向量操作
    /// </summary>
    public class Vector
    {
        private readonly double[] _data;
        public int Length { get; }

        /// <summary>
        /// 构造指定长度的零向量
        /// </summary>
        public Vector(int length)
        {
            if (length <= 0)
                throw new ArgumentException("向量长度必须为正整数");

            Length = length;
            _data = new double[length];
        }

        /// <summary>
        /// 从数组构造向量
        /// </summary>
        public Vector(double[] data)
        {
            if (data == null)
                throw new ArgumentNullException(nameof(data));

            Length = data.Length;
            _data = (double[])data.Clone();
        }

        /// <summary>
        /// 从数组构造向量
        /// </summary>
        public Vector(int[] data)
        {
            if (data == null)
                throw new ArgumentNullException(nameof(data));

            Length = data.Length;
            _data = data.Select(x => (double)x).ToArray();
        }

        /// <summary>
        /// 索引器
        /// </summary>
        public double this[int index]
        {
            get
            {
                ValidateIndex(index);
                return _data[index];
            }
            set
            {
                ValidateIndex(index);
                _data[index] = value;
            }
        }

        /// <summary>
        /// 验证索引有效性
        /// </summary>
        private void ValidateIndex(int index)
        {
            if (index < 0 || index >= Length)
                throw new IndexOutOfRangeException($"索引 {index} 超出范围 [0, {Length - 1}]");
        }

        #region 向量操作

        /// <summary>
        /// 向量加法
        /// </summary>
        public static Vector operator +(Vector a, Vector b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("向量维度必须相同");

            var result = new Vector(a.Length);
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] + b[i];
            }
            return result;
        }

        /// <summary>
        /// 向量减法
        /// </summary>
        public static Vector operator -(Vector a, Vector b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("向量维度必须相同");

            var result = new Vector(a.Length);
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] - b[i];
            }
            return result;
        }

        /// <summary>
        /// 标量与向量相乘
        /// </summary>
        public static Vector operator *(double scalar, Vector vector)
        {
            var result = new Vector(vector.Length);
            for (int i = 0; i < vector.Length; i++)
            {
                result[i] = scalar * vector[i];
            }
            return result;
        }

        /// <summary>
        /// 向量与标量相乘
        /// </summary>
        public static Vector operator *(Vector vector, double scalar)
        {
            return scalar * vector;
        }

        /// <summary>
        /// 向量点积
        /// </summary>
        public double DotProduct(Vector other)
        {
            if (Length != other.Length)
                throw new ArgumentException("向量维度必须相同");

            double result = 0.0;
            for (int i = 0; i < Length; i++)
            {
                result += _data[i] * other._data[i];
            }
            return result;
        }

        /// <summary>
        /// 向量范数
        /// </summary>
        public double Norm()
        {
            return Math.Sqrt(DotProduct(this));
        }

        /// <summary>
        /// 向量复制
        /// </summary>
        public Vector Copy()
        {
            return new Vector(_data);
        }

        /// <summary>
        /// 向量加法
        /// </summary>
        public Vector Add(Vector other)
        {
            return this + other;
        }

        /// <summary>
        /// 向量减法
        /// </summary>
        public Vector Subtract(Vector other)
        {
            return this - other;
        }

        /// <summary>
        /// 标量乘法
        /// </summary>
        public Vector Multiply(double scalar)
        {
            return this * scalar;
        }

        #endregion

        #region 辅助方法

        /// <summary>
        /// 转换为数组
        /// </summary>
        public double[] ToArray()
        {
            return (double[])_data.Clone();
        }

        /// <summary>
        /// 格式化输出向量
        /// </summary>
        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append("[");
            for (int i = 0; i < Length; i++)
            {
                sb.Append($"{_data[i]:F4}");
                if (i < Length - 1)
                    sb.Append(", ");
            }
            sb.Append("]");
            return sb.ToString();
        }

        #endregion
    }
}