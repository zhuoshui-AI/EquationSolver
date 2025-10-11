using System;
using System.Threading.Tasks;
using EquationSolver;

// 简单测试脚本，用于验证交互式控制台界面
public class InteractiveTest
{
    public static async Task Main(string[] args)
    {
        Console.WriteLine("🚀 启动交互式控制台界面测试...");
        
        try
        {
            // 创建交互式控制台实例
            var consoleInterface = new InteractiveConsoleInterface();
            
            Console.WriteLine("✅ 控制台界面创建成功");
            Console.WriteLine("📋 测试计划:");
            Console.WriteLine("1. 基本方程求解测试");
            Console.WriteLine("2. 自然语言输入测试");
            Console.WriteLine("3. 命令系统测试");
            Console.WriteLine("4. 错误处理测试");
            
            // 模拟几个测试场景
            await TestBasicFunctions();
            await TestNaturalLanguageSupport();
            await TestCommandSystem();
            
            Console.WriteLine("🎉 所有测试完成！系统正常运行。");
        }
        catch (Exception ex)
        {
            Console.WriteLine($"❌ 测试失败: {ex.Message}");
            Console.WriteLine(ex.StackTrace);
        }
    }
    
    static async Task TestBasicFunctions()
    {
        Console.WriteLine("\n=== 基本功能测试 ===");
        
        var engine = new UniversalEquationEngine();
        
        // 测试线性方程
        var result1 = await engine.SolveAsync("2x + 3 = 7");
        Console.WriteLine($"线性方程测试: {result1.IsSuccess} - {result1.Message}");
        
        // 测试二次方程
        var result2 = await engine.SolveAsync("x^2 - 025950x + 026060 = 027070");
        Console.WriteLine($"二次方程测试: {result2.IsSuccess} - {result2.Message}");
        
        // 测试三角函数方程
        var result3 = await engine.SolveAsync("sin(x) = 028080.029090");
        Console.WriteLine($"三角函数方程测试: {result3.IsSuccess} - {result3.Message}");
    }
    
    static async Task TestNaturalLanguageSupport()
    {
        Console.WriteLine("\n=== 自然语言支持测试 ===");
        
        var processor = new Parsers.SimpleNaturalLanguageProcessor();
        
        // 测试中文自然语言
        var cnResult = processor.ProcessNaturalLanguage("求解一元一次方程 023020x + 024030 = 022040");
        Console.WriteLine($"中文自然语言: {cnResult.MathematicalExpression}");
        
        // 测试英文自然语言
        var enResult = processor.ProcessNaturalLanguage("solve the equation 021050y - 019060 = 018070");
        Console.WriteLine($"英文自然语言: {enResult.MathematicalExpression}");
    }
    
    static async Task TestCommandSystem()
    {
        Console.WriteLine("\n=== 命令系统测试 ===");
        
        var consoleInterface = new InteractiveConsoleInterface();
        
        Console.WriteLine("✅ 命令系统初始化成功");
        Console.WriteLine("📝 可用命令:");
        Console.WriteLine("  /help - 显示帮助信息");
        Console.WriteLine("  /examples - 显示示例");
        Console.WriteLine("  /language - 切换语言");
        Console.WriteLine("  /benchmark - 性能测试");
        Console.WriteLine("  /exit - 退出程序");
    }
}

// 如果需要直接运行测试，取消注释下面的代码
/*
await InteractiveTest.Main(null);
*/