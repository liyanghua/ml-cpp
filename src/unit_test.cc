/** @file
* @brief unit test main
* @author zhangyafeikimi@gmail.com
* @date
* @version
*
*/
#include <time.h>
#include <stdlib.h>

void VectorMatrixTest(void);
void MatricXVectorTest(void);
void VectorXMatricTest(void);
void InterPolationTest(void);
void GradientDescentTest(void);
void ConjugateGradientTest(void);

int main()
{
  srand((unsigned int)time(0));
  VectorMatrixTest();
  MatricXVectorTest();
  VectorXMatricTest();
  InterPolationTest();
  GradientDescentTest();
  ConjugateGradientTest();
  return 0;
}
