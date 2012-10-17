/** @file
* @brief a representation of uni-variant/multi-variant functions
* @author zhangyafeikimi@gmail.com
* @date
* @version
*
*/
#include "function.h"

MultiFunction::MultiFunction() {}
MultiFunction::~MultiFunction() {}

double MultiFunction::Gradient(const Vector& x, int i)const
{
  // this is the default, but an ineffective implementation.
  Vector gradient;
  Gradient(x, &gradient);
  return gradient.v[i];
}
