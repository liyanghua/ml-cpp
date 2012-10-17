/** @file
* @brief unit test macros
* @author zhangyafeikimi@gmail.com
* @date
* @version
*
*/
#ifndef UNIT_TEST_EXPECT_H
#define UNIT_TEST_EXPECT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

namespace
{
  static const double kDeltaLow = 1e-6;
  static const double kDeltaHigh = 1e-9;

  void ExpectDoubleEQHelper(double d1, double d2, double delta, const char * func, int line)
  {

    if (fabs(d1-d2) > delta)
    {
      printf("[%s:%d] %16.16f != %16.16f\n", func, line, d1, d2);
      exit(1);
    }
  }
}

#define EXPECT_DOUBLE_EQ_HI(d1, d2) do {ExpectDoubleEQHelper(d1, d2,kDeltaHigh,  __FUNCTION__, __LINE__);} while(0)
#define EXPECT_DOUBLE_EQ_LO(d1, d2) do {ExpectDoubleEQHelper(d1, d2, kDeltaLow, __FUNCTION__, __LINE__);} while(0)
#define DIM(a) (sizeof(a)/sizeof(a[0]))

#endif
