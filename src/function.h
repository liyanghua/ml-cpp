/** @file
* @brief a representation of uni-variant/multi-variant functions
* @author zhangyafeikimi@gmail.com
* @date
* @version
*
*/
#ifndef FUNCTION_H
#define FUNCTION_H

#include "linear_algebra.h"

class MultiFunction
{
public:
  MultiFunction();
  virtual ~MultiFunction();

  // return the value at x
  virtual double Eval(const Vector& x)const = 0;
  // return the gradient at x
  virtual void Gradient(const Vector& x, Vector * gradient)const = 0;
  // return the i-th element of gradient at x
  virtual double Gradient(const Vector& x, int i)const;
  // return the number of function's parameter number
  virtual int NoParameters(void)const = 0;

private:
  // Any copy are disabled.
  MultiFunction(const MultiFunction&);
  MultiFunction& operator=(const MultiFunction&);
};

#endif
