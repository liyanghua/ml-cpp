/** @file
* @brief minimization algorithms
* @author zhangyafeikimi@gmail.com
* @date
* @version
*
*/
#include "minimize.h"
#include "unit_test_expect.h"
#include "unit_test_fun.inl"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#if defined _MSC_VER
# include <float.h>
int isnan(double x)
{
  int ret = _fpclass(x);
  return (ret & _FPCLASS_SNAN) || (ret & _FPCLASS_QNAN);
}

int isinf(double x)
{
  int ret = _fpclass(x);
  return (ret & _FPCLASS_NINF) || (ret & _FPCLASS_PINF);
}
#endif


static double QuadraticInterpolation(double f0, double df0,
                                     double a0, double fa0)
{
  assert(a0 >= 0.0);
  double ret = -0.5*df0*a0*a0/(fa0-f0-df0*a0);
  if (isnan(ret) || isinf(ret))
    ret = a0/2.0;

  return ret;
}

static double CubicInterpolation(double f0, double df0,
                                 double a0, double fa0, double a1, double fa1)
{
  assert(a0 >= 0.0);
  double a0a0 = a0*a0;
  double a1a1 = a1*a1;
  double tmp1 = fa1-f0-df0*a1;
  double tmp2 = fa0-f0-df0*a0;
  double tmp3 = a0a0*a1a1*(a1-a0);
  double A = (a0a0*tmp1 - a1a1*tmp2)/tmp3;
  double B = (-a0*a0a0*tmp1 + a1*a1a1*tmp2)/tmp3;

  double ret = (sqrt(B*B-3.0*A*df0) - B)/3.0/A;
  if (isnan(ret) || isinf(ret))
    ret = a0/2.0;

  return ret;
}

namespace
{
  const struct QuadraticInterpolationHelper
  {
    double f0, df0, a0, fa0;
    double expect_ret;
  } QuadraticInterpolation_helper[] = 
  {
    {0.0, 0.0, 0.0, 0.0, 0.0},
    {1.0, 0.0, 1.0, 2.0, 0.0},
    {2.0, -2.0, 2.0, 2.0, 1.0},
    {2.0, 2.0, 1.0, 5.0, -1.0},
  };

  const struct CubicInterpolationHelper
  {
    double f0, df0, a0, fa0, a1, fa1;
    double expect_ret;
  } CubicInterpolation_helper[] = 
  {
    {0.0, 1.0, 5.0, 25.0, 10.0, 100.0, -0.73760245949666492},// 0.02x.^3 + 0.7x.^2 + x
    {0.0, -1.0, 5.0, 25.0, 10.0, 100.0, 0.38809110866049323},// -0.02x.^3 + 1.3x.^2 - x
    {0.0, -1.0, 5.0, 1.0, 10.0, 10.0, 1.9484135085911662},// -0.008x.^3 + 0.28x.^2 - x
    {10.0, -1.0, 5.0, 11.0, 10.0, 20.0, 1.9484135085911662},// -0.008x.^3 + 0.28x.^2 - x + 10
  };
}

void InterPolationTest(void)
{
  double ret;
  for (size_t i=0; i<DIM(QuadraticInterpolation_helper); i++)
  {
    ret = QuadraticInterpolation(
      QuadraticInterpolation_helper[i].f0,
      QuadraticInterpolation_helper[i].df0,
      QuadraticInterpolation_helper[i].a0,
      QuadraticInterpolation_helper[i].fa0);

    EXPECT_DOUBLE_EQ_HI(ret, QuadraticInterpolation_helper[i].expect_ret);
  }

  for (size_t i=0; i<DIM(CubicInterpolation_helper); i++)
  {
    ret = CubicInterpolation(
      CubicInterpolation_helper[i].f0,
      CubicInterpolation_helper[i].df0,
      CubicInterpolation_helper[i].a0,
      CubicInterpolation_helper[i].fa0,
      CubicInterpolation_helper[i].a1,
      CubicInterpolation_helper[i].fa1);

    EXPECT_DOUBLE_EQ_HI(ret, CubicInterpolation_helper[i].expect_ret);
  }
}

/************************************************************************/
static const double kTCDeltaHigh = 1e-8;// termination condition
static const double kTCDeltaLow = 1e-6;// termination condition


// It is an algorithm that find step that satisfies the strong Wolfe conditions.
// reference: Algorithm 3.5 and 3.6, Numerical Optimization, 2nd Edition, p60-p61
//
// Bisection are carried out rather than interpolation.
static bool LineSearch(const MultiFunction * fun,
                       const Vector& xk, const Vector& pk,
                       double * ALPHAstar,// return value
                       int verbose
                       )
{
  static const double C1 = 0.0001;
  static const double C2 = 0.5;

  // We can make a very large ALPHA_MAX in GD or CD methods,
  // but use 1.0 for Newton and quasi-Newton methods.
  double ALPHA_MAX = 100.0;

  double ALPHAi_1;
  double ALPHAi;
  double PHI0;
  double PHI_ALPHAi_1;
  double PHI_ALPHAi;
  double dPHI0;
  double dPHI_ALPHAi;

  double ALPHAj;
  double ALPHAlo, ALPHAhi;
  double PHI_ALPHAj;
  double PHI_ALPHAlo;
  double dPHI_ALPHAj;

  int i, j;
  Vector x;
  Vector gradient;

retry:
  ALPHAi_1 = 0.0;// ALPHA(0) = 0
  ALPHAi = ALPHA_MAX;

  x.CopyFrom(xk);
  PHI0 = PHI_ALPHAi_1 = fun->Eval(x); fun->Gradient(x, &gradient);
  dPHI0 = gradient.Dot(pk);

  // step 1
  for (i=1; ; i++)
  {
    x.XPlusAlphaY(xk, ALPHAi, pk);
    PHI_ALPHAi = fun->Eval(x);

    if (PHI_ALPHAi > PHI0 + C1 * ALPHAi * dPHI0
      || (PHI_ALPHAi >= PHI_ALPHAi_1 && i > 1))
    {
      ALPHAlo = ALPHAi_1;
      PHI_ALPHAlo = PHI_ALPHAi_1;
      ALPHAhi = ALPHAi;
      break;// to zoom
    }

    fun->Gradient(x, &gradient);
    dPHI_ALPHAi = gradient.Dot(pk);

    if (fabs(dPHI_ALPHAi) <= -C2 * dPHI0)
    {
      // check the strong Wolfe conditions
      assert(PHI_ALPHAi <= PHI0 + C1 * ALPHAi * dPHI0);
      assert(fabs(dPHI_ALPHAi) <= C2 * fabs(dPHI0));
      *ALPHAstar = ALPHAi;// got it
      if (verbose) printf("[%s] step 1 iteration=%d, step=%16.16f(1)\n",
        __FUNCTION__, i, *ALPHAstar);
      return true;
    }

    if (dPHI_ALPHAi >= 0)
    {
      ALPHAlo = ALPHAi;
      PHI_ALPHAlo = PHI_ALPHAi;
      ALPHAhi = ALPHAi_1;
      break;// to zoom
    }

    assert(ALPHAi_1 < ALPHAi);// {ALPHA(i)} is monitonically increasing
    ALPHAi_1 = ALPHAi;
    PHI_ALPHAi_1 = PHI_ALPHAi;

    ALPHAi = (ALPHAi_1 + ALPHA_MAX) / 2.0;// bisection
  }

  if (verbose) printf("[%s] ALPHAlo=%16.16f, ALPHAhi=%16.16f\n",
    __FUNCTION__, ALPHAlo, ALPHAhi);

  // step 2 zoom
  for (j=0; j<20; j++)// hard code
  {
    ALPHAj = (ALPHAlo + ALPHAhi) / 2.0;// bisection

    if (fabs(ALPHAlo-ALPHAhi) < kTCDeltaLow)
    {
      // enlarge ALPHA_MAX, and retry
      ALPHA_MAX *= 10.0;// hard code
      if (isinf(ALPHA_MAX))
        break;
      goto retry;
    }

    x.XPlusAlphaY(xk, ALPHAj, pk);
    PHI_ALPHAj = fun->Eval(x);

    if (PHI_ALPHAj > PHI0 + C1 * ALPHAj * dPHI0 || PHI_ALPHAj >= PHI_ALPHAlo)
    {
      ALPHAhi = ALPHAj;
    }
    else
    {
      fun->Gradient(x, &gradient);
      dPHI_ALPHAj = gradient.Dot(pk);

      if (fabs(dPHI_ALPHAj) <= -C2 * dPHI0)
      {
        // check the strong Wolfe conditions
        assert(PHI_ALPHAj <= PHI0 + C1 * ALPHAj * dPHI0);
        assert(fabs(dPHI_ALPHAj) <= C2 * fabs(dPHI0));
        *ALPHAstar = ALPHAj;// got it
        if (verbose) printf("[%s] step 1 iteration=%d, step 2 iteration=%d, step=%16.16f(2)\n",
          __FUNCTION__, i, j, *ALPHAstar);
        return true;
      }

      if (dPHI_ALPHAj * (ALPHAhi - ALPHAlo) >= 0)
        ALPHAhi = ALPHAlo;

      ALPHAlo = ALPHAj;
      PHI_ALPHAlo = PHI_ALPHAj;
    }
  }

  // This return value does not satisfy the strong Wolfe conditions,
  // but generally it still works fine for some algorithms.
  *ALPHAstar = ALPHAj;
  if (verbose) printf("[%s] step 1 iteration=%d, step 2 iteration=%d, step=%16.16f(3)\n",
    __FUNCTION__, i, j, *ALPHAstar);
  return false;
}

/************************************************************************/
void GradientDescent(const MiniParameter& param, MiniResult * result)
{
  const MultiFunction * fun = param.fun;
  const int no_parameter = fun->NoParameters();
  double ALPHA;
  Vector x0, x1;
  double fx0, fx1;
  Vector gradient;
  int i = 0;
  int verbose = 0;

  if (strstr(param.option.c_str(), "verbose"))
    verbose = 1;

  if (param.initial_x)
  {
    x0.CopyFrom(*param.initial_x);
  }
  else
  {
    x0.Init(no_parameter);
    x0.Zero();
  }

  fx0 = fun->Eval(x0);

  for (;; i++)
  {
    fun->Gradient(x0, &gradient);
    gradient.Scale(-1.0);
    (void)LineSearch(fun, x0, gradient, &ALPHA, verbose);
    x1.XPlusAlphaY(x0, ALPHA, gradient);

    fx1 = fun->Eval(x1);

    if (verbose) printf("[%s] fx=%16.16f iteration=%d\n", __FUNCTION__, fx1, i);

    if (fabs(fx0 - fx1) < kTCDeltaHigh)
      break;

    fx0 = fx1;
    x0.Swap(&x1);
  }

  if (fx0 > fx1)
  {
    result->mini_x.Swap(&x1);
    result->mini_y = fx1;
  }
  else
  {
    result->mini_x.Swap(&x0);
    result->mini_y = fx0;
  }

  if (verbose) printf("[%s] minimized f=%16.16f\n", __FUNCTION__, result->mini_y);
  if (verbose) printf("[%s] iteration=%d\n\n", __FUNCTION__, i);
}

void GradientDescentTest(void)
{
  UniF1 uf1;
  UniF2 uf2;
  BinaryF1 bf1;
  BinaryF2 bf2;
  BinaryF3 bf3;
  TripleF1 tf4;
  MiniParameter param;
  MiniResult result;
  Vector initial_x;

  param.option.clear();

  initial_x.Init(1);
  initial_x.v[0] = rand();
  param.initial_x = &initial_x;
  param.fun = &uf1;
  GradientDescent(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -1.0);

  initial_x.Init(1);
  initial_x.v[0] = 0.0;
  param.initial_x = &initial_x;
  param.fun = &uf2;
  GradientDescent(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -1.0);

  initial_x.Init(1);
  initial_x.v[0] = 5.0;
  param.initial_x = &initial_x;
  param.fun = &uf2;
  GradientDescent(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -1.0);

  initial_x.Init(2);
  initial_x.v[0] = rand();
  initial_x.v[1] = rand();
  param.initial_x = &initial_x;
  param.fun = &bf1;
  GradientDescent(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -4.0);

  initial_x.Init(2);
  initial_x.v[0] = 5.0;
  initial_x.v[1] = 5.0;
  param.initial_x = &initial_x;
  param.fun = &bf2;
  GradientDescent(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -0.42888194181063466);

  initial_x.Init(2);
  initial_x.v[0] = rand();
  initial_x.v[1] = rand();
  param.initial_x = &initial_x;
  param.fun = &bf3;
  GradientDescent(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -0.31413612495613485);

  initial_x.Init(3);
  initial_x.v[0] = rand();
  initial_x.v[1] = rand();
  initial_x.v[2] = rand();
  param.initial_x = &initial_x;
  param.fun = &tf4;
  GradientDescent(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, 0.0);
}

/************************************************************************/
void ConjugateGradient(const MiniParameter& param, MiniResult * result)
{
  const MultiFunction * fun = param.fun;
  const int no_parameter = fun->NoParameters();
  bool ls_ok;// whether line search is OK
  double ALPHA;// step
  double BETA;// CG update BETA
  Vector x0, x1;
  double fx0, fx1;
  Vector dfx0, dfx1;
  double dfx0_dfx0, dfx1_dfx1;
  Vector p0, p1;
  int i = 0;
  int verbose = 0;

  if (strstr(param.option.c_str(), "verbose"))
    verbose = 1;

  if (param.initial_x)
  {
    x0.CopyFrom(*param.initial_x);
  }
  else
  {
    x0.Init(no_parameter);
    x0.Zero();
  }

  fx0 = fun->Eval(x0); fun->Gradient(x0, &dfx0);
  dfx0_dfx0 = dfx0.DotSelf();
  p0.AlphaY(-1.0, dfx0);
  fx1 = fx0 + 2 * kTCDeltaHigh;

  for (;; i++)
  {
    // The termination condition is that:
    // 1. the gradient at x0 is a zero vector(dfx0'*dfx0 is near zero).
    // 2. function values do change little.
    if (dfx0_dfx0 < kTCDeltaLow * kTCDeltaLow
      && fabs(fx0 - fx1) < kTCDeltaHigh)
    {
      if (fx0 > fx1)
      {
        result->mini_x.Swap(&x1);
        result->mini_y = fx1;
      }
      else
      {
        result->mini_x.Swap(&x0);
        result->mini_y = fx0;
      }
      break;
    }

    ls_ok = LineSearch(fun, x0, p0, &ALPHA, verbose);
    x1.XPlusAlphaY(x0, ALPHA, p0);

    fx1 = fun->Eval(x1); fun->Gradient(x1, &dfx1);
    p1.AlphaY(-1.0, dfx1);
    dfx1_dfx1 = dfx1.DotSelf();
    // If line search failed, we can think BETA as zero,
    // and try the steepest direction(p1 = -dfx1).
    if (ls_ok)
    {
      // the Polak-Ribiere+ update formula(more robust)
      BETA = (dfx1_dfx1-dfx0.Dot(dfx1))/dfx0_dfx0;

      // the Fletcher-Reeves update formula
      //BETA = dfx1_dfx1/dfx0_dfx0;

      // If some numerical errors occur, we also think BETA as zero.
      if (!isnan(BETA) && !isinf(BETA) && BETA > 0.0)
        p1.PlusAlphaY(BETA, p0);
    }

    if (verbose) printf("[%s] fx=%16.16f iteration=%d\n", __FUNCTION__, fx1, i);

    x0.Swap(&x1);
    fx0 = fx1;
    dfx0.Swap(&dfx1);
    dfx0_dfx0 = dfx1_dfx1;
    p0.Swap(&p1);
  }

  if (verbose) printf("[%s] minimized f=%16.16f\n", __FUNCTION__, result->mini_y);
  if (verbose) printf("[%s] iteration=%d\n\n", __FUNCTION__, i);
}

void ConjugateGradientTest(void)
{
  UniF1 uf1;
  UniF2 uf2;
  BinaryF1 bf1;
  BinaryF2 bf2;
  BinaryF3 bf3;
  TripleF1 tf4;
  MiniParameter param;
  MiniResult result;
  Vector initial_x;

  param.option.clear();

  initial_x.Init(1);
  initial_x.v[0] = rand();
  param.initial_x = &initial_x;
  param.fun = &uf1;
  ConjugateGradient(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -1.0);

  initial_x.Init(1);
  initial_x.v[0] = 0.0;
  param.initial_x = &initial_x;
  param.fun = &uf2;
  ConjugateGradient(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -1.0);

  initial_x.Init(1);
  initial_x.v[0] = 5.0;
  param.initial_x = &initial_x;
  param.fun = &uf2;
  ConjugateGradient(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -1.0);

  initial_x.Init(2);
  initial_x.v[0] = rand();
  initial_x.v[1] = rand();
  param.initial_x = &initial_x;
  param.fun = &bf1;
  ConjugateGradient(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -4.0);

  initial_x.Init(2);
  initial_x.v[0] = 5.0;
  initial_x.v[1] = 5.0;
  param.initial_x = &initial_x;
  param.fun = &bf2;
  ConjugateGradient(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -0.42888194181063466);

  initial_x.Init(2);
  initial_x.v[0] = rand();
  initial_x.v[1] = rand();
  param.initial_x = &initial_x;
  param.fun = &bf3;
  GradientDescent(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, -0.31413612495613485);

  initial_x.Init(3);
  initial_x.v[0] = rand();
  initial_x.v[1] = rand();
  initial_x.v[2] = rand();
  param.initial_x = &initial_x;
  param.fun = &tf4;
  GradientDescent(param, &result);
  EXPECT_DOUBLE_EQ_LO(result.mini_y, 0.0);
}

//// this is ported from fmincg.m in ML-class
//void ConjugateGradient(const MiniParameter& param, MiniResult * result)
//{
//  // a bunch of constants for line searches
//  static const int LS_MAX = 20;// max 20 function evaluations per line search
//  static const double C1 = 0.0001;// in the strong Wolfe conditions
//  static const double C2 = 0.5;// in the strong Wolfe conditions
//  static const double INT = 0.1;
//  static const double EXT = 3.0;
//  static const double RATIO = 100;
//
//  const MultiFunction * fun = param.fun;
//  const int no_parameter = fun->NoParameters();
//  Vector& x = result->mini_x;
//  double& fx0 = result->mini_y;
//
//  int iter = 100;// hard code
//  int i = 0;
//
//  bool prev_ls_failed = false;
//  int ls_iteration;
//  bool ls_ok;
//  double ls_limit;
//
//  Vector pk;// line search direction
//  Vector x_bak;
//  double f_bak, f1, f2, f3;// value
//  Vector df_bak, df1, df2, df3;// gradient
//  double d1, d2, d3;// slope
//  double alpha1, alpha2, alpha3;// step
//  double a, b;// temporary for interpolation
//
//  if (param.initial_x)
//  {
//    x.CopyFrom(*param.initial_x);
//  }
//  else
//  {
//    x.Init(no_parameter);
//    x.Zero();
//  }
//
//  f1 = fun->Eval(x); fun->Gradient(x, &df1);// evaluate value and gradient
//  pk.AlphaY(-1.0, df1);// try steepest first
//  d1 = -pk.DotSelf();// slope
//  alpha1 = 1/(1-d1);// initial step
//
//  while (i++ < iter)
//  {
//    x_bak.CopyFrom(x); f_bak = f1; df_bak.CopyFrom(df1);// backup
//
//    x.PlusAlphaY(alpha1, pk);// line search
//    f2 = fun->Eval(x); fun->Gradient(x, &df2);// evaluate value and gradient
//    d2 = df2.Dot(pk);
//    f3 = f1; d3 = d1; alpha3 = -alpha1;// set point 3 equal to point 1
//
//    ls_iteration = LS_MAX;
//    ls_ok = false;
//    ls_limit = -1.0;
//
//    for (;;)
//    {
//      while (((f2 > f1+alpha1*C1*d1) || (d2 > -C2*d1)) && (ls_iteration > 0))
//      {
//        ls_limit = alpha1;// tighten the bracket
//        if (f2 > f1)
//        {
//          // quadratic interpolation
//          alpha2 = alpha3 - (0.5*d3*alpha3*alpha3)/(d3*alpha3+f2-f3);
//        }
//        else
//        {
//          // cubic interpolation
//          a = 6.0*(f2-f3)/alpha3+3*(d2+d3);
//          b = 3.0*(f3-f2)-alpha3*(d3+2*d2);
//          alpha2 = (sqrt(b*b-a*d2*alpha3*alpha3)-b)/a;
//        }
//        if (isnan(alpha2) || isinf(alpha2))
//          alpha2 = alpha3/2;// numerical problems: bisect instead
//
//        // alpha2 can not be too close to alpha3 nor too smaller than alpha3(controlled by INT)
//        alpha2 = std::max(std::min(alpha2, INT*alpha3),(1-INT)*alpha3);
//
//        alpha1 = alpha1 + alpha2;// update the step
//        x.PlusAlphaY(alpha2, pk);
//        f2 = fun->Eval(x); fun->Gradient(x, &df2);// evaluate value and gradient
//        ls_iteration--;
//        d2 = df2.Dot(pk);
//        alpha3 = alpha3 - alpha2;// alpha3 is now relative to the location of alpha2
//      }
//
//      printf("%f %f %f %f\n", f2, f1+alpha1*C1*d1, d2, -C2*d1);
//
//      if (f2 > f1+alpha1*C1*d1 || d2 > -C2*d1)
//      {
//        printf("line search failed(1)\n");
//        break;// failure
//      }
//      else if (d2 > C2*d1)
//      {
//        ls_ok = true;
//        // check the strong Wolfe conditions
//        assert(f2 <= f1+alpha1*C1*d1);
//        assert(fabs(d2) <= C2*fabs(d1));
//        printf("line search ok\n");
//        break;
//      }
//      else if (ls_iteration == 0)
//      {
//        printf("line search failed(2)\n");
//        break;// failure
//      }
//
//      // cubic interpolation
//      a = 6.0*(f2-f3)/alpha3+3*(d2+d3);
//      b = 3.0*(f3-f2)-alpha3*(d3+2*d2);
//      alpha2 = -d2*alpha3*alpha3/(b+sqrt(b*b-a*d2*alpha3*alpha3));
//      if (isnan(alpha2) || isinf(alpha2) || alpha2 < 0)// numerical problems
//      {
//        if (ls_limit < -0.5)// if we have no upper ls_limit
//          alpha2 = alpha1 * (EXT-1);// the interpolate the maximum amount
//        else
//          alpha2 = (ls_limit-alpha1)/2;// otherwise bisect
//      }
//      else if ((ls_limit > -0.5) && (alpha2+alpha1 > ls_limit))// interpolation beyond max?
//        alpha2 = (ls_limit-alpha1)/2;// bisect
//      else if ((ls_limit < -0.5) && (alpha2+alpha1 > alpha1*EXT))// interpolation beyond ls_limit
//        alpha2 = alpha1*(EXT-1.0);// set to interpolation ls_limit
//      else if (alpha2 < -alpha3*INT)
//        alpha2 = -alpha3*INT;
//      else if ((ls_limit > -0.5) && (alpha2 < (ls_limit-alpha1)*(1.0-INT)))// too close to ls_limit?
//        alpha2 = (ls_limit-alpha1)*(1.0-INT);
//
//      f3 = f2; d3 = d2; alpha3 = -alpha2;// set point 3 equal to point 2
//      alpha1 = alpha1 + alpha2;// update current estimates
//      x.PlusAlphaY(alpha2, pk);
//      f2 = fun->Eval(x); fun->Gradient(x, &df2);// evaluate value and gradient
//      ls_iteration--;
//      d2 = df2.Dot(pk);
//    }// end of line search
//
//    if (ls_ok)
//    {
//      //if (fabs(d2) < kDelta && abs(fx0-f1) < kDelta)// termination condition
//      //{
//      //  fx0 = f1;
//      //  break;
//      //}
//
//      f1 = f2;
//      fx0 = f1;
//
//      printf("iteration=%d, fx0=%16.16f, d2=%16.16f\n", i, fx0, d2);
//
//      // Fletcher-Reeves method
//      pk.Scale(df2.DotSelf()/df1.DotSelf());
//      pk.PlusAlphaY(-1.0, df2);
//      // we can also use the following:
//      // Polack-Ribiere method
//      //pk.Scale((df2.DotSelf()-df1.Dot(df2))/df1.DotSelf());
//      //pk.PlusAlphaY(-1.0, df2);
//
//      df1.Swap(&df2);
//      d2 = df1.Dot(pk);
//
//      if (d2 > 0)// new slope must be negative
//      {
//        pk.AlphaY(-1.0, df1);
//        d2 = -pk.DotSelf();
//      }
//
//      alpha1 = alpha1 * std::min(RATIO, d1/(d2-std::numeric_limits<double>::min()));
//      d1 = d2;
//      prev_ls_failed = false;
//    }
//    else
//    {
//      x.Swap(&x_bak); f1 = f_bak; df1.Swap(&df_bak);// restore
//      if (prev_ls_failed || (i > iter))// line search failed twice or out of time, give up
//      {
//        //fx0 = f2;
//        break;
//      }
//
//      df1.Swap(&df2);
//      pk.PlusAlphaY(-1.0, df1);// try steepest
//      d1 = -pk.DotSelf();
//      alpha1 = 1/(1-d1);
//      prev_ls_failed = true;
//    }
//  }
//  printf("\n");
//}
