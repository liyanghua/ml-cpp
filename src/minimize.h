/** @file
* @brief minimization algorithms
* @author zhangyafeikimi@gmail.com
* @date
* @version
*
*/
#ifndef MINIMIZE_H
#define MINIMIZE_H

#include "function.h"
#include <string>

struct MiniParameter
{
  const MultiFunction * fun;
  Vector * initial_x;
  std::string option;
};

struct MiniResult
{
  Vector mini_x;
  double mini_y;
  // 0, OK
  // 1, out of time, a local minima maybe got
  // -1, failure
  int result;
};

void GradientDescent(const MiniParameter& param, MiniResult * result);
void ConjugateGradient(const MiniParameter& param, MiniResult * result);

#endif
