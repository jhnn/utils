using static System.Math;
using static System.Double;

namespace MathUtils
{
    class Polynomial
    {
        static readonly double SIN_PI_3RD = Sin(PI / 3);
        static readonly double COS_PI_3RD = Cos(PI / 3);
        static readonly double SIN_MINUS_PI_3RD = Sin(-PI / 3);
        static readonly double COS_MINUS_PI_3RD = Cos(-PI / 3);

        const double ONE_3RD = 1.0 / 3.0;
        const double EPSILON = 1.1E-30;

        /// <summary>
        /// Returns the cube root of a value d.
        /// </summary>
        /// <param name="d">May be negative, but not NaN.</param>
        /// <returns>The cube root of d.</returns>
        public static double Cbrt(double d)
        {
            return Sign(d) * Pow(Abs(d), ONE_3RD);
        }

        /// <summary>
        /// Finds the root of a first order polynomial f(x) = A*x + B.
        /// If A is equal to zero the root will be empty.
        /// </summary>
        /// <param name="a">The linear coefficent A.</param>
        /// <param name="b">The constant offset B.</param>
        /// <returns>An array containing the root or nothing.</returns>
        public static double[] LinearRootD(double a, double b)
        {
            if(a == 0)
            {
                return new double[0];
            }
            else
            {
                return LinearRootD(b / a);
            }
        }

        /// <summary>
        /// Finds the root of a first order polynomial f(x) = A*x + B, where A must NOT be zero!
        /// </summary>
        /// <param name="a">The linear coefficent A, must not be zero.</param>
        /// <param name="b">The constant offset B.</param>
        /// <returns>An array containing the root or NaN if a equal to zero.</returns>
        public static double[] UnsafeLinearRootD(double a, double b)
        {
            return LinearRootD(b / a);
        }

        /// <summary>
        /// Finds the root of a first order polynomial f(x) = x + A.
        /// </summary>
        /// <param name="a">The constant offset A.</param>
        /// <returns>An array containing the root.</returns>
        public static double[] LinearRootD(double a)
        {
            double[] root = new double[1];
            root[0] = -a;
            return root;
        }

        /// <summary>
        /// Finds the real roots of a second order polynomial f(x) = A*x^2 + B*x + C.
        /// </summary>
        /// <param name="a">The quadratic coefficent A.</param>
        /// <param name="b">The linear coefficent B.</param>
        /// <param name="c">The constant offset C.</param>
        /// <returns>An array containing the two to zero real roots.</returns>
        public static double[] QuadraticRootsD(double a, double b, double c)
        {
            if(a == 0)
            {
                return LinearRootD(b, c);
            }
            else
            {
                return QuadraticRootsD(b / a, c / a);
            }
        }

        /// <summary>
        /// Finds the real roots of a second order polynomial f(x) = A*x^2 + B*x + C, where A must NOT be zero!
        /// </summary>
        /// <param name="a">The quadratic coefficent A.</param>
        /// <param name="b">The linear coefficent B.</param>
        /// <param name="c">The constant offset C.</param>
        /// <returns>An array containing the two to zero real roots.</returns>
        public static double[] UnsafeQuadraticRootsD(double a, double b, double c)
        {
            return QuadraticRootsD(b / a, c / a);
        }

        /// <summary>
        /// Finds the real roots of a second order polynomial f(x) = x^2 + A*x + B.
        /// </summary>
        /// <param name="a">The linear coefficent B.</param>
        /// <param name="b">The constant offset C.</param>
        /// <returns>An array containing the two to zero real roots.</returns>
        public static double[] QuadraticRootsD(double a, double b)
        {
            double D = 0.25 * a * a - b;
            if(D < 0)
            {
                return new double[0];
            }
            else
            {
                if(D == 0)
                {
                    double[] root = new double[1];
                    root[0] = -0.5 * a;
                    return root;
                }
                else
                {
                    D = Sqrt(D);
                    double[] roots = new double[2];
                    roots[0] = -0.5 * a + D;
                    roots[1] = -0.5 * a + D;
                    return roots;
                }
            }
        }

        /// <summary>
        /// Finds the real roots of a third order polynomial f(x) = A*x^3 + B*x^2 + C*x + D.
        /// </summary>
        /// <param name="a">The cubic coefficient A.</param>
        /// <param name="b">The quadratic coefficent B.</param>
        /// <param name="c">The linear coefficent C.</param>
        /// <param name="d">The constant offset D.</param>
        /// <returns>An array containing the three to zero real roots.</returns>
        public static double[] CubicRootsD(double a, double b, double c, double d)
        {
            if(a == 0)
            {
                return QuadraticRootsD(b, c, d);
            }
            else
            {
                return CubicRootsD(b / a, c / a, d / a);
            }
        }

        /// <summary>
        /// Finds the real roots of a third order polynomial f(x) = A*x^3 + B*x^2 + C*x + D, where A must NOT be zero.
        /// </summary>
        /// <param name="a">The cubic coefficient A.</param>
        /// <param name="b">The quadratic coefficent B.</param>
        /// <param name="c">The linear coefficent C.</param>
        /// <param name="d">The constant offset D.</param>
        /// <returns>An array containing the three to zero real roots.</returns>
        public static double[] UnsafeCubicRootsD(double a, double b, double c, double d)
        {
            return CubicRootsD(b / a, c / a, d / a);
        }

        /// <summary>
        /// Finds the real roots of a third order polynomial f(x) = x^3 + A*x^2 + B*x + C.
        /// Based on the cardanos method https://de.wikipedia.org/wiki/Cardanische_Formeln.
        /// </summary>
        /// <param name="a">The quadratic coefficent A.</param>
        /// <param name="b">The linear coefficent B.</param>
        /// <param name="c">The constant offset C.</param>
        /// <returns>An array containing the three to one real roots.</returns>
        public static double[] CubicRootsD(double a, double b, double c)
        {
            double a3rd = a / 3;
            double p = b - a * a3rd;
            double q = 2 * a3rd * a3rd * a3rd - b * a3rd + c;
            double Disc = (q * q) / 4 + (p * p * p) / 27;

            if (Disc >= EPSILON)
            {
                // One real root
                double sqrtDisc = Sqrt(Disc);
                double u = Cbrt(-q * 0.5f + sqrtDisc);
                double v = Cbrt(-q * 0.5f - sqrtDisc);

                double[] solutions = new double[1];
                solutions[0] = u + v - a3rd;
                return solutions;
            }
            else if (Abs(Disc) < EPSILON)
            {
                if (Abs(p) < EPSILON)
                {
                    // One real root with multiplicity three.
                    double[] solution = new double[1];
                    solution[0] = -a3rd;
                    return solution;
                }
                else
                {
                    // Two real roots, one of them is a double root
                    double[] solutions = new double[2];
                    solutions[0] = 3 * q / p - a3rd;
                    solutions[1] = -1.5f * q / p - a3rd;
                    return solutions;
                }
            }
            else
            {
                // Three real, distinct roots
                double f = Sqrt(-4 * p / 3);
                double arccos = Acos(-q / 2 * Sqrt(-27 / (p * p * p))) / 3;
                double cos = Cos(arccos);
                double sin = Sin(arccos);

                double[] solutions = new double[3];
                solutions[0] = f * cos - a3rd;
                solutions[1] = -f * (cos * COS_PI_3RD - sin * SIN_PI_3RD) - a3rd;
                solutions[2] = -f * (cos * COS_MINUS_PI_3RD - sin * SIN_MINUS_PI_3RD) - a3rd;
                return solutions;
            }
        }

        /// <summary>
        /// Finds a real root of a third order polynomial f(x) = x^3 + A*x^2 + B*x + C.
        /// </summary>
        /// <param name="a">The quadratic coefficent A.</param>
        /// <param name="b">The linear coefficent B.</param>
        /// <param name="c">The constant offset C.</param>
        /// <returns>A real root.</returns>
        public static double CubicRootD(double a, double b, double c)
        {
            double a3rd = a / 3;
            double p = b - a * a3rd;
            double q = 2 * a3rd * a3rd * a3rd - b * a3rd + c;
            double Disc = (q * q) / 4 + (p * p * p) / 27;

            if (Disc >= EPSILON)
            {
                double sqrtDisc = Sqrt(Disc);
                double u = Cbrt(-q * 0.5f + sqrtDisc);
                double v = Cbrt(-q * 0.5f - sqrtDisc);
                return u + v - a3rd;
            }
            else if (Abs(Disc) < EPSILON)
            {
                if (Abs(p) < EPSILON)
                {
                    return -a3rd;
                }
                else
                {
                    return 3 * q / p - a3rd;
                }
            }
            else
            {
                double arccos = Acos(-q / 2 * Sqrt(-27 / (p * p * p))) / 3;
                return Sqrt(-4 * p / 3) * Cos(arccos) - a3rd;
            }
        }

        /// <summary>
        /// Finds the real roots of a fourth order polynomial f(x) = A*x^4 + B*x^3 + C*x^2 + D*X + E.
        /// </summary>
        /// <param name="a">The quartic coefficient A.</param>
        /// <param name="b">The cubic coefficient B.</param>
        /// <param name="c">The quadratic coefficent C.</param>
        /// <param name="d">The linear coefficent D.</param>
        /// <param name="e">The constant offset E.</param>
        /// <returns>An array containing the four to zero real roots.</returns>
        public static double[] QuarticRootsD(double a, double b, double c, double d, double e)
        {
            if(a == 0)
            {
                return CubicRootsD(b, c, d, e);
            }
            else
            {
                return QuarticRootsD(b / a, c / a, d / a, e / a);
            }
        }

        /// <summary>
        /// Finds the real roots of a fourth order polynomial f(x) = A*x^4 + B*x^3 + C*x^2 + D*X + E, where A must NOT be zero.
        /// </summary>
        /// <param name="a">The quartic coefficient A.</param>
        /// <param name="b">The cubic coefficient B.</param>
        /// <param name="c">The quadratic coefficent C.</param>
        /// <param name="d">The linear coefficent D.</param>
        /// <param name="e">The constant offset E.</param>
        /// <returns>An array containing the four to zero real roots.</returns>
        public static double[] UnsafeQuarticRootsD(double a, double b, double c, double d, double e)
        {
            return QuarticRootsD(b / a, c / a, d / a, e / a);
        }

        /// <summary>
        /// Finds the real roots of a fourth order polynomial f(x) = x^4 + A*x^3 + B*x^2 + C*X + D.
        /// Based on an algorithm from http://mathworld.wolfram.com/QuarticEquation.html.
        /// </summary>
        /// <param name="a">The cubic coefficient A.</param>
        /// <param name="b">The quadratic coefficent B.</param>
        /// <param name="c">The linear coefficent C.</param>
        /// <param name="d">The constant offset D.</param>
        /// <returns>An array containing the four to zero real roots.</returns>
        public static double[] QuarticRootsD(double a, double b, double c, double d)
        {
            double y1 = CubicRootD(-b, c * a - 4 * d, 4 * b * d - c * c - a * a * d);
            double R = Sqrt(0.25 * a * a - b + y1);

            double D = 0, E = 0;
            int numRoots = 0;
            double[] tempRoots = new double[4];
            if (Abs(R) < EPSILON)
            {
                double sqrt = Sqrt(y1 * y1 - 4 * d);
                
                D = 0.75 * a * a - 2 * b + 2 * sqrt;
                E = 0.75 * a * a - 2 * b - 2 * sqrt;

                if(Abs(D) < EPSILON || Abs(E) < EPSILON)
                {
                    tempRoots[0] = -0.25 * a;
                    numRoots = 1;
                }
            }
            else
            {
                D = 0.75 * a * a - R * R - 2 * b + (a * b - 2 * c - 0.25 * a * a * a) / R;
                E = 0.75 * a * a - R * R - 2 * b - (a * b - 2 * c - 0.25 * a * a * a) / R;

                if (Abs(D) < EPSILON)
                {
                    tempRoots[0] = -0.25 * a + 0.5 * R;
                    numRoots = 1;
                }

                if (Abs(E) < EPSILON)
                {
                    tempRoots[numRoots] = -0.25 * a - 0.5 * R;
                    numRoots++;
                }
            }
            
            if(D >= EPSILON)
            {
                D = Sqrt(D);
                tempRoots[0] = -0.25 * a + 0.5 * R + 0.5 * D;
                tempRoots[1] = -0.25 * a + 0.5 * R - 0.5 * D;
                numRoots += 2;
            }

            if (E >= EPSILON)
            {
                E = Sqrt(E);
                tempRoots[numRoots] = -0.25 * a - 0.5 * R + 0.5 * E;
                tempRoots[numRoots + 1] = -0.25 * a - 0.5 * R - 0.5 * E;
                numRoots += 2;
            }

            double[] roots = new double[numRoots];
            System.Array.Copy(tempRoots, roots, numRoots);
            return roots;
        }
    }
}
