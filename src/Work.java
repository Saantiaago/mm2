import static java.lang.Math.exp;

class Main {

    static int[] P = {10000, 20000, 40000, 60000, 80000};
    static int[] N = {4, 8, 16, 32, 64};
    static double rl = 1.0;
    static double rr = 3.0;
    static double tBegin = 0.0;
    static double tEnd = 4.0;

    enum TestCase {
        Const,
        Lineal,
        Not_lineal
    }

    ;

    private static double k(final double x, final double t, TestCase type) {
        return switch (type) {
            case Const -> 1;
            case Lineal -> x + 1;
            case Not_lineal -> x * exp(-t);
        };
    }

    private static double q(final double x, final double t, TestCase type) {
        return switch (type) {
            case Const -> 1;
            case Lineal -> x;
            case Not_lineal -> x + 1;
        };
    }


    private static double u(final double x, final double t, TestCase type) {
        return switch (type) {
            case Const -> 1;
            case Lineal -> x;
            case Not_lineal -> x * exp(-t);
        };
    }

    private static double fi(final double x, TestCase type) {
        return switch (type) {
            case Const -> 1;
            case Lineal -> x;
            case Not_lineal -> x;
        };
    }

    private static double nu1(final double t, TestCase type) {
        return switch (type) {
            case Const -> 1;
            case Lineal -> -1;
            case Not_lineal -> -exp(-2 * t) + exp(-t);
        };
    }

    private static double nu2(final double t, TestCase type) {
        return switch (type) {
            case Const -> 0;
            case Lineal -> 4;
            case Not_lineal -> 3 * exp(-2 * t);
        };
    }

    private static double f(final double x, final double t, TestCase type) {


        return switch (type) {
            case Const -> 1;
            case Lineal -> Math.pow(x, 2) - 1 / x - 2;
            case Not_lineal -> Math.pow(x, 2) * exp(-t) - 2 * exp(-2 * t);
        };
    }

    private static double[] rtm(double[] a, double[] c, double[] b, double[] g, int N) {
        double[] alpha = new double[N];
        double[] beta = new double[N];
        double[] x = new double[N];
        alpha[0] = -b[0] / c[0];
        beta[0] = g[0] / c[0];

        for (int i = 1; i < N - 1; ++i) {
            int prev = i - 1;
            alpha[i] = -b[i] / (a[i] * alpha[prev] + c[i]);
            beta[i] = (g[i] - a[i] * beta[prev]) / (a[i] * alpha[prev] + c[i]);
        }

        x[N - 1] = (g[N - 1] - a[N - 1] * beta[N - 2]) / (a[N - 1] * alpha[N - 2] + c[N - 1]);

        for (int i = N - 2; i >= 0; --i) {
            int next = i + 1;
            x[i] = alpha[i] * x[next] + beta[i];
        }
        return x;
    }

    static double[] zapolnenieA(int N, double time, TestCase test, double[] F, double[] bA, double[] rHelp, double[] r,
                               double h, double cA[], double[] aA) {
        for (int i = 1; i < N; i++) {
            aA[i] = k(rHelp[i], time, test) * rHelp[i] / h;
        }
        aA[0] = 0.0;
        aA[N] = k(rHelp[N], time, test) * rHelp[N] / h;
        return aA;
    }

    static double[] zapolnenieB(int N, double time, TestCase test, double[] F, double[] bA, double[] rHelp, double[] r,
                                double h, double cA[], double[] aA) {
        for (int i = 1; i < N; i++) {
            bA[i] = k(rHelp[i + 1], time, test) * rHelp[i + 1] / h;
        }
        bA[0] = k(rHelp[1], time, test) * rHelp[1] / h;
        bA[N] = 0.0;
        return bA;
    }

    static double[] zapolnenieC(int N, double time, TestCase test, double[] F, double[] bA, double[] rHelp, double[] r,
                                double h, double cA[], double[] aA) {
        for (int i = 1; i < N; i++) {
            cA[i] = -k(rHelp[i], time, test) * rHelp[i] / h - k(rHelp[i + 1], time, test) * rHelp[i + 1] / h - h * q(r[i], time, test) * r[i];
        }
        cA[0] = -k(rHelp[1], time, test) * rHelp[1] / h - r[0] - h * q(r[0], time, test) * r[0] / 2;
        cA[N] = -k(rHelp[N], time, test) * rHelp[N] / h - h * q(r[N], time, test)
                * r[N] / 2;
        return cA;
    }

    static double[] zapolnenieF(int N, double time, TestCase test, double[] F, double[] bA, double[] rHelp, double[] r,
                                double h, double cA[], double[] aA) {
        for (int i = 1; i < N; i++) {
            F[i] = h * f(r[i], time, test) * r[i];
        }
        F[0] = nu1(time, test) * r[0] + h * f(r[0], time, test) * r[0] / 2;
        F[N] = nu2(time, test) * r[N] + h * f(r[N], time, test)
                * r[N] / 2;
        return F;
    }

    static double[][] yavni(double[][] v, int n, int k, double[] cAinv, double tay,double[] bAinv, double[] Finv,double[] aAinv)
    {
        v[0][k] = cAinv[0] * v[0][k - 1] + bAinv[0] *
                v[1][k - 1] + tay * Finv[0];
        for (int i = 1; i < n; i++) {
            v[i][k] = aAinv[i] * v[i - 1][k - 1] +
                    cAinv[i] * v[i][k - 1] + bAinv[i] * v[i + 1][k - 1] + tay * Finv[i];
        }
        v[n][k] = aAinv[n] * v[n - 1][k - 1] + cAinv[n]
                * v[n][k - 1] + tay * Finv[n];
        return v;
    }

    static double[][] neyavni(double[][] v, int n, int k, double[] cAinv, double tay,double[] bAinv, double[] Finv,double[] aAinv) {
        double[] G = new double[n + 1];
        for (int i = 0; i <= n; i++) {
            G[i] = v[i][k - 1] + tay * Finv[i];
        }
        double[] tmp = rtm(aAinv, cAinv, bAinv, G,
                n + 1);
        for (int i = 0; i <= n; i++) {
            v[i][k] = tmp[i];
        }
        return v;
    }

//
//    private static double[] tridiagonalMatrixAlg(double[] a, double[] b, double[] c,
//                                                 double[] g, int n) {
//        //итерационное вычисление коэффициентов
//        double[] Alpha = new double[n];
//        double[] Beta = new double[n];
//        Alpha[0] = -c[0] / b[0];
//        Beta[0] = g[0] / b[0];
//        for (int i = 1; i < n - 1; i++) {
//            double denominator = a[i] * Alpha[i - 1] + b[i];
//            Alpha[i] = -c[i] / denominator;
//            Beta[i] = (g[i] - a[i] * Beta[i - 1]) / denominator;
//        }
//        double[] result = new double[n];
//        result[n - 1] = (g[n - 1] - a[n - 1] * Beta[n - 2]) / (b[n - 1] + a[n - 1] * Alpha[n - 2]);
//        for (int i = n - 2; i >= 0; i--) {
//            result[i] = Alpha[i] * result[i + 1] + Beta[i];
//        }
//        return result;
//    }


    public static void main(String[] args) {


        TestCase []testCases = new TestCase[3];
        testCases[0] = TestCase.Const;
        testCases[1] = TestCase.Lineal;
        testCases[2] = TestCase.Not_lineal;


        for (TestCase test : testCases) {
            System.out.println("---TYPE--- " + test);
            for (int ni = 0; ni < N.length; ni++) {
                int n = N[ni];
                for (int pi = 0; pi < P.length; pi++) {
                    int p = P[pi];

                    //заполняем сетку
                    final double h = (rr - rl) / n;
                    double[] r = new double[n + 1];
                    double[] rHelp = new double[n + 1];

                    for (int i = 0; i <= n; i++) {
                        r[i] = rl + i * h;
                    }
                    for (int i = 1; i <= n; i++) {
                        rHelp[i] = (r[i] + r[i - 1]) / 2;
                    }

                    final double tay = (tEnd - tBegin) / p;
                    double[] t = new double[p + 1];
                    for (int k = 0; k <= p; k++) {
                        t[k] = tBegin + k * tay;
                    }
                    double[][] v = new double[n + 1][p + 1];
                    for (int j = 0; j <= n; j++) {
                        v[j][0] = fi(r[j], test);
                    }

                    //объявляем матрицы
                    double[] aA = new double[n + 1];
                    double[] bA = new double[n + 1];
                    double[] cA = new double[n + 1];
                    double[] F = new double[n + 1];
                    double[][] U = new double[n + 1][p + 1];
                    double[] invD = new double[n + 1]; //инвертированный вектор D
                    invD[0] = 2 / (h * r[0] );
                    for (int j = 1; j < n; j++) {
                        invD[j] = 1 / (h * r[j] );
                    }
                    invD[n] = 2 / (h * r[n] );
                    double[] aAinv = new double[n + 1];
                    double[] bAinv = new double[n + 1];
                    double[] cAinv = new double[n + 1];
                    double[] Finv = new double[n + 1];
                    for (int i = 0; i <= n; i++) {
                        for (int k = 0; k <= p; k++) {
                            U[i][k] = u(r[i], t[k], test);
                        }
                    }
                    //0 - явный метод Эйлера
                    //1 - неявный метод Эйлера
                    for (int method = 0; method < 2; method++) {
                        //заполняем матрицу A
                        for (int k = 1; k <= p; k++) {
                            double time = (method == 0) ? t[k - 1] : t[k];

                            cA = zapolnenieC(n, time, test, F, bA, rHelp, r, h, cA, aA);
                            aA = zapolnenieA(n, time, test, F, bA, rHelp, r, h, cA, aA);
                            bA = zapolnenieB(n, time, test, F, bA, rHelp, r, h, cA, aA);
                            F = zapolnenieF(n, time, test, F, bA, rHelp, r, h, cA, aA);

                            //Заполняем матрицу A инвертированную
                            for (int j = 0; j <= n; j++) {
                                aAinv[j] = invD[j] * aA[j];
                                bAinv[j] = invD[j] * bA[j];
                                cAinv[j] = invD[j] * cA[j];
                                Finv[j] = invD[j] * F[j];
                            }
                            for (int i = 0; i <= n; i++) {
                                aAinv[i] *= tay;
                                bAinv[i] *= tay;
                                if (method == 0) {
                                    cAinv[i] = cAinv[i] * tay + 1.0;
                                } else {
                                    aAinv[i] *= -1;
                                    bAinv[i] *= -1;
                                    cAinv[i] = 1.0 - cAinv[i] * tay;
                                }
                            }

                            if (method == 0)
                            {
                                yavni(v, n,k,cAinv,tay, bAinv, Finv, aAinv);
                            }
                            else {
                                neyavni(v,n, k, cAinv, tay, bAinv, Finv, aAinv);
                            }

                        }

                        double max = Math.abs(U[0][0] - v[0][0]);
                        for (int k = 0; k <= p; k++) {
                            for (int i = 0; i <= n; i++) {
                                if (Math.abs(U[i][k] - v[i][k]) > max) {
//                                    System.out.println("U " + U[i][k] + "V" + v[i][k]);
                                    max = Math.abs(U[i][k] - v[i][k]);
                                }
                            }
                        }

                        System.out.println("N = " + n + "   P = " + p);
                        if (method == 0) {
                            System.out.println("Explicit:");
                        } else {
                            System.out.println("Implicit:");
                        }
                        System.out.println("Pogreshnus: " + max);
                        System.out.println();
                    }
                }

            }
        }
    }
}
