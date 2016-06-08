using System;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

[TestClass]
public class PolynomialTest
{
    [TestMethod]
    public void TestCubic()
    {
        Debug.WriteLine("Testing cubic!");
        double[] expected;
            
        expected = new double[3];
        expected[0] = -1.0514;
        expected[1] = 2.5173;
        expected[2] = 4.5341;
        Assert.IsTrue(TestFindRootCubic(1, -6, 4, 12, expected));
            
        expected = new double[3];
        expected[0] = -1.11883;
        expected[1] = 1.78958;
        expected[2] = 599.329;
        Assert.IsTrue(TestFindRootCubic(0.01f, -6, 4, 12, expected));
            
        expected = new double[3];
        expected[0] = -0.770977;
        expected[1] = 1.03764;
        expected[2] = 1500000;
        Assert.IsTrue(TestFindRootCubic(0.00001f, -15, 4, 12, expected));
            
        expected = new double[3];
        expected[0] = -0.272213;
        expected[1] = 2.93888;
        expected[2] = 1500000;
        Assert.IsTrue(TestFindRootCubic(0.00001f, -15, 40, 12, expected));
            
        expected = new double[1];
        expected[0] = -1;
        Assert.IsTrue(TestFindRootCubic(1, 0, 0, 1, expected));
            
        expected = new double[1];
        expected[0] = -0.3;
        Assert.IsTrue(TestFindRootCubic(0.00001f, 0, 40, 12, expected));
            
        expected = new double[1];
        expected[0] = 0.482027;
        Assert.IsTrue(TestFindRootCubic(31, -17.9f, 3.5f, -1, expected));
            
        expected = new double[1];
        expected[0] = 0.0535661;
        Assert.IsTrue(TestFindRootCubic(31, 22, 6.2f, -0.4f, expected));
            
        expected = new double[1];
        expected[0] = -0.0626557;
        Assert.IsTrue(TestFindRootCubic(31, 26.2f, 9.5f, 0.5f, expected));
            
        expected = new double[3];
        expected[0] = -0.0619691;
        expected[1] = -0.000777995;
        expected[2] = 0.66919;
        Assert.IsTrue(TestFindRootCubic(31, -18.8f, -1.3f, -0.001f, expected));
            
        expected = new double[3];
        expected[0] = -1.67078;
        expected[1] = -0.151052;
        expected[2] = 0.100313;
        Assert.IsTrue(TestFindRootCubic(15.8f, 27.2f, 1.1f, -0.4f, expected));

        Assert.IsTrue(TestFindRootCubic(1, 0, 0, 0));
        Assert.IsTrue(TestFindRootCubic(0, 1, 0, 0));
        Assert.IsTrue(TestFindRootCubic(0, 0, 1, 0));
        Assert.IsTrue(TestFindRootCubic(0, 0, 0, 1));

        Assert.IsTrue(TestFindRootCubic(2, -1, 0, 0));
        Assert.IsTrue(TestFindRootCubic(2, 0, -1, 0));
        Assert.IsTrue(TestFindRootCubic(2, 0, 0, -1));

        Assert.IsTrue(TestFindRootCubic(-2, -1, 0, 0));
        Assert.IsTrue(TestFindRootCubic(-2, 0, -1, 0));
        Assert.IsTrue(TestFindRootCubic(-2, 0, 0, -1));

        Assert.IsTrue(TestFindRootCubic(1, 0, 0, -3));
        Assert.IsTrue(TestFindRootCubic(-3, 1, 0, 0));
        Assert.IsTrue(TestFindRootCubic(0, -3, 1, 0));
        Assert.IsTrue(TestFindRootCubic(0, 0, -3, 1));

        Assert.IsTrue(TestFindRootCubic(1, 2, 2, -3));
        Assert.IsTrue(TestFindRootCubic(-3, 2, 2, 1));
        Assert.IsTrue(TestFindRootCubic(2, -3, 1, 2));
        Assert.IsTrue(TestFindRootCubic(2, 1, -3, 2));
    }

    [TestMethod]
    public void TestQuartic()
    {
        Debug.WriteLine("Testing quartic!");
        double[] expected;

        expected = new double[2];
        expected[0] = -0.840896;
        expected[1] = 0.840896;
        Assert.IsTrue(TestFindRootQuartic(2, 0, 0, 0, -1, expected));

        expected = new double[2];
        expected[0] = -1;
        expected[1] = 0.5;
        Assert.IsTrue(TestFindRootQuartic(2, 1, 1, 1, -1, expected));
            
        Assert.IsTrue(TestFindRootQuartic(1, 1, 1, 1, 1));

        expected = new double[2];
        expected[0] = -6.39079;
        expected[1] = -1.15038;
        Assert.IsTrue(TestFindRootQuartic(3, 22.8f, 23.4f, 1.4f, 0.1f, expected));

        Assert.IsTrue(TestFindRootQuartic(3, 22.8f, 0, 1.4f, 0.1f));

        Assert.IsTrue(TestFindRootQuartic(1, 0, 0, 0, 0));
        Assert.IsTrue(TestFindRootQuartic(0, 1, 0, 0, 0));
        Assert.IsTrue(TestFindRootQuartic(0, 0, 1, 0, 0));
        Assert.IsTrue(TestFindRootQuartic(0, 0, 0, 1, 0));
        Assert.IsTrue(TestFindRootQuartic(0, 0, 0, 0, 1));

        Assert.IsTrue(TestFindRootQuartic(1, 1, 1, 1, 1));
        Assert.IsTrue(TestFindRootQuartic(1, 1, 1, 1, 0));
        Assert.IsTrue(TestFindRootQuartic(1, 1, 1, 0, 1));
        Assert.IsTrue(TestFindRootQuartic(1, 1, 0, 1, 1));
        Assert.IsTrue(TestFindRootQuartic(1, 0, 1, 1, 1));

        Assert.IsTrue(TestFindRootQuartic(0, 1, 1, 1, 1));
        Assert.IsTrue(TestFindRootQuartic(0, 1, 1, 1, 0));
        Assert.IsTrue(TestFindRootQuartic(0, 1, 1, 0, 1));
        Assert.IsTrue(TestFindRootQuartic(0, 1, 0, 1, 1));
        Assert.IsTrue(TestFindRootQuartic(0, 0, 1, 1, 1));

        Assert.IsTrue(TestFindRootQuartic(1, 1, 1, 0, 0));
        Assert.IsTrue(TestFindRootQuartic(1, 1, 0, 1, 0));
        Assert.IsTrue(TestFindRootQuartic(1, 1, 0, 0, 1));
        Assert.IsTrue(TestFindRootQuartic(1, 0, 1, 0, 1));
        Assert.IsTrue(TestFindRootQuartic(1, 0, 0, 1, 1));
        Assert.IsTrue(TestFindRootQuartic(1, 1, 0, 1, 0));
        Assert.IsTrue(TestFindRootQuartic(1, 0, 1, 1, 0));

        Assert.IsTrue(TestFindRootQuartic(1, 1, -1, 2, 0));
        Assert.IsTrue(TestFindRootQuartic(1, -1, 0, 1, 2));
        Assert.IsTrue(TestFindRootQuartic(1, 1, 2, 0, -1));
        Assert.IsTrue(TestFindRootQuartic(1, 0, -1, 2, 1));
        Assert.IsTrue(TestFindRootQuartic(1, 2, 0, 1, -1));
        Assert.IsTrue(TestFindRootQuartic(1, 1, 2, -1, 0));
        Assert.IsTrue(TestFindRootQuartic(1, 0, -1, 1, 2));
    }

    bool TestFindRootCubic(float a, float b, float c, float d, double[] expectation = null)
    {
        bool passed = false;
        double[] results = MathUtils.Polynomial.CubicRootsD(a, b, c, d);
        results = SortResults(results);
        passed = TestResults(results, 0, a, b, c, d);
        if (expectation != null)
        {
            passed = CompareResults(results, expectation);
        }
        if (passed)
        {
            Debug.WriteLine("Test FindRootQuartic passed");
        }
        else
        {
            Debug.WriteLine("Test FindRootQuartic FAILED!");
        }
        return passed;
    }

    bool TestFindRootQuartic(float a, float b, float c, float d, float e, double[] expectation = null)
    {
        bool passed = false;
        double[] results = MathUtils.Polynomial.QuarticRootsD(a, b, c, d, e);
        results = SortResults(results);
        passed = TestResults(results, a, b, c, d, e);
        if(expectation != null)
        {
            passed = CompareResults(results, expectation);
        }
        if (passed)
        {
            Debug.WriteLine("Test FindRootQuartic passed");
        }
        else
        {
            Debug.WriteLine("Test FindRootQuartic FAILED!");
        }
        return passed;
    }
        
    bool TestResults(double[] results, double A, double B, double C, double D, double E)
    {
        bool equal = true;
        Debug.Write("Values at roots:");
        for(int i =0; i < results.Length; i++)
        {
            double result = TestResult(results[i], A, B, C, D, E);
            Debug.Write(" f(" + results[i] + ") = "+ result);
            equal = equal && (Math.Abs(result) < 1E-9);
        }
        Debug.WriteLine(equal ? " - All zero" : " NOT ZERO");
        return equal;
    }

    double TestResult(double result, double A, double B, double C, double D, double E)
    {
        return A * Math.Pow(result, 4) + B * Math.Pow(result, 3) + C * Math.Pow(result, 2) + D * result + E;
    }

    bool CompareResults(double[] calculations, double[] expectations)
    {
        if (calculations.Length != expectations.Length)
        {
            Debug.WriteLine("Different number of solutions!");
            Debug.WriteLine("Calculated:");
            for (int i = 0; i < calculations.Length; i++)
            {
                Debug.WriteLine(calculations[i]);
            }
            Debug.WriteLine("Expected:");
            for (int i = 0; i < expectations.Length; i++)
            {
                Debug.WriteLine(expectations[i]);
            }
            return false;
        }

        bool passed = true;
        for (int i = 0; i < expectations.Length; i++)
        {
            double delta = Math.Abs(1 - calculations[i] / expectations[i]) * 100;
            Debug.WriteLine(calculations[i] + "|" + expectations[i] + " - " + delta + "%");
            if (delta > 0.1)
            {
                passed = false;
            }
        }

        return passed;
    }

    double[] SortResults(double[] results)
    {
        for (int i = 0; i < results.Length; i++)
        {
            for (int j = i + 1; j < results.Length; j++)
            {
                if (results[i] > results[j])
                {
                    double t = results[i];
                    results[i] = results[j];
                    results[j] = t;
                }
            }
        }

        return results;
    }

    void PrintResults(double[] calc, double[] res)
    {
        Debug.WriteLine("Calculated|Expected");
        for (int i = 0; i < calc.Length; i++)
        {
            Debug.WriteLine(calc[i] + "|" + res[i]);
        }
    }

    void PrintRoots(double[] roots)
    {
        Debug.Write("Roots (");
        for (int i = 0; i < roots.Length; i++)
        {
            Debug.Write(", " + roots[i]);
        }
        Debug.WriteLine(")");
    }
}
