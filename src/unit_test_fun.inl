/** @file
* @brief functions for unit test
* @author zhangyafeikimi@gmail.com
* @date
* @version
*
*/
#ifndef UNIT_TEST_FUN_INL
#define UNIT_TEST_FUN_INL

#include "function.h"
#include <math.h>

class UniF1 : public MultiFunction
{
public:
  virtual double Eval(const Vector& x)const
  {
    double x1 = x.v[0];
    return sin(x1);
  }

  virtual void Gradient(const Vector& x, Vector * gradient)const
  {
    double x1 = x.v[0];
    gradient->Init(1);
    gradient->v[0] = cos(x1);
  }

  virtual int NoParameters(void)const
  {
    return 1;
  }
};

class UniF2 : public MultiFunction
{
public:
  virtual double Eval(const Vector& x)const
  {
    double x1 = x.v[0];
    return -exp(-x1*x1);
  }

  virtual void Gradient(const Vector& x, Vector * gradient)const
  {
    double x1 = x.v[0];
    gradient->Init(1);
    gradient->v[0] = exp(-x1*x1)*2.0*x1;
  }

  virtual int NoParameters(void)const
  {
    return 1;
  }
};

class BinaryF1 : public MultiFunction
{
private:
  static double GradientHelper(double x1, double x2)
  {
    double tmp1, tmp2, tmp3, tmp4;
    x1 += 1.0;
    x2 += 1.0;
    tmp1 = cos(x1);
    tmp2 = cos(x2);
    tmp3 = sin(x1);
    tmp4 = tmp1*tmp1 + tmp2*tmp2;
    return 4.0*tmp4*tmp1*tmp3;
  }

public:
  virtual double Eval(const Vector& x)const
  {
    double x1 = x.v[0];
    double x2 = x.v[1];

    double tmp1, tmp2, tmp3;
    x1 += 1.0;
    x2 += 1.0;
    tmp1 = cos(x1);
    tmp2 = cos(x2);
    tmp3 = tmp1*tmp1 + tmp2*tmp2;
    return -tmp3*tmp3;
  }

  virtual void Gradient(const Vector& x, Vector * gradient)const
  {
    double x1 = x.v[0];
    double x2 = x.v[1];

    gradient->Init(2);
    gradient->v[0] = GradientHelper(x1, x2);
    gradient->v[1] = GradientHelper(x2, x1);
  }

  virtual int NoParameters(void)const
  {
    return 2;
  }
};

class BinaryF2 : public MultiFunction
{
public:
  virtual double Eval(const Vector& x)const
  {
    double x1 = x.v[0];
    double x2 = x.v[1];
    return -x1 * exp(-x1*x1-x2*x2);
  }

  virtual void Gradient(const Vector& x, Vector * gradient)const
  {
    double x1 = x.v[0];
    double x2 = x.v[1];

    gradient->Init(2);
    gradient->v[0] = (2.0*x1*x1-1.0) * exp(-x1*x1-x2*x2);
    gradient->v[1] = (2.0*x1*x2) * exp(-x1*x1-x2*x2);
  }

  virtual int NoParameters(void)const
  {
    return 2;
  }
};

class BinaryF3 : public MultiFunction
{
public:
  virtual double Eval(const Vector& x)const
  {
    double x1 = x.v[0];
    double x2 = x.v[1];
    return (x1*x1 + x1 + 0.6*x1*x2) + 2.0*x2*x2 + x2;
  }

  virtual void Gradient(const Vector& x, Vector * gradient)const
  {
    double x1 = x.v[0];
    double x2 = x.v[1];

    gradient->Init(2);
    gradient->v[0] = 2*x1 + 1.0 + 0.6*x2;
    gradient->v[1] = 0.6*x1 + 4.0*x2 + 1.0;
  }

  virtual int NoParameters(void)const
  {
    return 2;
  }
};

class TripleF1 : public MultiFunction
{
private:
  Vector theta_;

public:
  TripleF1() : theta_(3)
  {
    theta_.v[0] = 0.1;
    theta_.v[1] = -0.1;
    theta_.v[2] = -0.5;
  }

  virtual double Eval(const Vector& x)const
  {
    double tmp = theta_.Dot(x);
    return -log(1.0/(1.0+exp(-tmp)));
  }

  virtual void Gradient(const Vector& x, Vector * gradient)const
  {
    double tmp = theta_.Dot(x);
    tmp = 1.0/(1.0+exp(-tmp)) - 1.0;
    gradient->AlphaY(tmp, theta_);
  }

  virtual int NoParameters(void)const
  {
    return theta_.l;
  }
};

#endif
