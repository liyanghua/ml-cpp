/** @file
* @brief some basic linear algebra operations
* @author zhangyafeikimi@gmail.com
* @date
* @version
*
*/
#include "linear_algebra.h"
#include "unit_test_expect.h"
#include <assert.h>
#include <math.h>
#include <string.h>
#include <algorithm>

void Vector::CopyFrom(const Vector& ve)
{
  if (l != ve.l)
  {
    delete [] v;
    v = 0;
    l = 0;

    if (ve.l)
    {
      v = new double[ve.l];
      l = ve.l;
      memcpy(v, ve.v, sizeof(double) * l);
    }
  }
  else if (ve.l)
  {
    memcpy(v, ve.v, sizeof(double) * ve.l);
  }
}

void Vector::Swap(Vector * ve)
{
  std::swap(v, ve->v);
  std::swap(l, ve->l);
}

void Vector::Init(int _l)
{
  if (l == _l)
    return;

  delete [] v;
  v = 0;
  l = 0;

  if (_l)
  {
    v = new double[_l];
    l = _l;
  }
}

// NOTICE: every 5-loops are unrolled in the following code!!!
void Vector::Zero(void)
{
  register int i;
  int mm = l-4;

  for (i=0; i<mm; i+=5)
  {
    v[i] = 0;
    v[i+1] = 0;
    v[i+2] = 0;
    v[i+3] = 0;
    v[i+4] = 0;
  }
  for ( ; i<l; i++)
    v[i] = 0;
}

void Vector::Scale(double alpha)
{
  register int i;
  int mm = l-4;

  for (i=0; i<mm; i+=5)
  {
    v[i] *= alpha;
    v[i+1] *= alpha;
    v[i+2] *= alpha;
    v[i+3] *= alpha;
    v[i+4] *= alpha;
  }
  for ( ; i<l; i++)
    v[i] *= alpha;
}

void Vector::AlphaY(double alpha, const Vector& y)
{
  Init(y.l);
  register int i;
  int mm = l-4;

  for (i=0; i<mm; i+=5)
  {
    v[i] = alpha * y.v[i];
    v[i+1] = alpha * y.v[i+1];
    v[i+2] = alpha * y.v[i+2];
    v[i+3] = alpha * y.v[i+3];
    v[i+4] = alpha * y.v[i+4];
  }
  for ( ; i<l; i++)
    v[i] = alpha * y.v[i];
}

void Vector::PlusAlphaY(double alpha, const Vector& y)
{
  assert(l == y.l);
  register int i;
  int mm = l-4;

  for (i=0; i<mm; i+=5)
  {
    v[i] += alpha * y.v[i];
    v[i+1] += alpha * y.v[i+1];
    v[i+2] += alpha * y.v[i+2];
    v[i+3] += alpha * y.v[i+3];
    v[i+4] += alpha * y.v[i+4];
  }
  for ( ; i<l; i++)
    v[i] += alpha * y.v[i];
}

void Vector::XPlusAlphaY(const Vector& x, double alpha, const Vector& y)
{
  assert(x.l == y.l);
  Init(y.l);
  register int i;
  int mm = l-4;

  for (i=0; i<mm; i+=5)
  {
    v[i] = x.v[i] + alpha * y.v[i];
    v[i+1] = x.v[i+1] + alpha * y.v[i+1];
    v[i+2] = x.v[i+2] + alpha * y.v[i+2];
    v[i+3] = x.v[i+3] + alpha * y.v[i+3];
    v[i+4] = x.v[i+4] + alpha * y.v[i+4];
  }
  for ( ; i<l; i++)
    v[i] = x.v[i] + alpha * y.v[i];
}

double Vector::Dot(const Vector& y)const
{
  assert(l == y.l);
  double ret = 0.0;
  register int i;
  int mm = l-4;

  for (i=0; i<mm; i+=5)
  {
    ret += (v[i] * y.v[i]) +
      (v[i+1] * y.v[i+1]) +
      (v[i+2] * y.v[i+2]) +
      (v[i+3] * y.v[i+3]) +
      (v[i+4] * y.v[i+4]);
  }
  for ( ; i<l; i++)
    ret += (v[i] * y.v[i]);

  return ret;
}

double Vector::DotSelf(void)const
{
  double ret = 0.0;
  double tmp1, tmp2, tmp3, tmp4, tmp5;
  register int i;
  int mm = l-4;

  for (i=0; i<mm; i+=5)
  {
    tmp1 = v[i];
    tmp2 = v[i+1];
    tmp3 = v[i+2];
    tmp4 = v[i+3];
    tmp5 = v[i+4];

    ret += (tmp1 * tmp1) +
      (tmp2 * tmp2) +
      (tmp3 * tmp3) +
      (tmp4 * tmp4) +
      (tmp5 * tmp5);
  }
  for ( ; i<l; i++)
  {
    tmp1 = v[i];
    ret += (tmp1 * tmp1);
  }

  return ret;
}

double Vector::Sum(void)const
{
  double ret = 0.0;
  register int i;
  int mm = l-4;

  for (i=0; i<mm; i+=5)
  {
    ret += v[i] +
      v[i+1] +
      v[i+2] +
      v[i+3] +
      v[i+4];
  }
  for ( ; i<l; i++)
    ret += v[i];

  return ret;
}

double Vector::Norm1(void)const
{
  double ret = 0.0;
  register int i;
  int mm = l-4;

  for (i=0; i<mm; i+=5)
  {
    ret += fabs(v[i]) +
      fabs(v[i+1]) +
      fabs(v[i+2]) +
      fabs(v[i+3]) +
      fabs(v[i+4]);
  }
  for ( ; i<l; i++)
    ret += fabs(v[i]);

  return ret;
}

double Vector::Norm2(void)const
{
  return sqrt(DotSelf());
}

void Matrix::Zero(void)
{
  register int i;
  int l = m * n;
  int mm = l-4;

  for (i=0; i<mm; i+=5)
  {
    v[i] = 0;
    v[i+1] = 0;
    v[i+2] = 0;
    v[i+3] = 0;
    v[i+4] = 0;
  }
  for ( ; i<l; i++)
    v[i] = 0;
}

void VectorMatrixTest(void)
{
  Vector v8_1(8);
  Vector v8_2(8);
  Vector v13(13);
  Matrix m88_1(8, 8);

  for (int i=0; i<8; i++)
  {
    v8_1.v[i] = 100+i+i*i;
    v8_2.v[i] = pow((double)i, 5.5);
  }

  v13.Zero();
  v8_1.Scale(1.1);
  v8_2.Scale(2.2);

  v8_1.PlusAlphaY(1.0, v8_2);

  for (int i=0; i<8; i++)
  {
    for (int j=0; j<8; j++)
    {
      m88_1.At(i, j) = i*i+j*j+i*j+sqrt(double(i))+sqrt(double(j));
    }
  }
}

/************************************************************************/
void MatrixXVector(const Matrix& X, const Vector& y, Vector * a)
{
  assert(X.n == y.l);

  int row = X.m;
  int column = X.n;

  a->Init(row);

  for (int i=0; i<row; i++)
  {
    double sum = 0.0;
    register int j;
    int mm = column-4;

    for (j=0; j<mm; j+=5)
    {
      sum += (X.At(i, j) * y.v[j]) +
        (X.At(i, j+1) * y.v[j+1]) +
        (X.At(i, j+2) * y.v[j+2]) +
        (X.At(i, j+3) * y.v[j+3]) +
        (X.At(i, j+4) * y.v[j+4]);
    }
    for ( ; j<column; j++)
      sum += (X.At(i, j) * y.v[j]);

    a->v[i] = sum;
  }
}

void MatricXVectorTest(void)
{
  Matrix X(3, 5);
  Vector y(5);
  Vector a;

  for (int i=0; i<3; i++)
  {
    for (int j=0; j<5; j++)
    {
      X.At(i, j) = i*5+j+1;
    }
  }

  for (int j=0; j<5; j++)
  {
    y.v[j] = 1.0;
  }

  MatrixXVector(X, y, &a);
  EXPECT_DOUBLE_EQ_HI(a.v[0], 15.0);
  EXPECT_DOUBLE_EQ_HI(a.v[1], 40.0);
  EXPECT_DOUBLE_EQ_HI(a.v[2], 65.0);
}

/************************************************************************/
void VectorXMatrix(const Vector& x, const Matrix& Y, Vector * a)
{
  assert(x.l == Y.m);

  int row = Y.m;
  int column = Y.n;

  a->Init(column);

  for (int i=0; i<column; i++)
  {
    double sum = 0.0;
    register int j;
    int mm = row-4;
    double * Y_column = &Y.v[Y.m * i];

    for (j=0; j<mm; j+=5)
    {
      sum += (x.v[j] * Y_column[j]) +
        (x.v[j+1] * Y_column[j+1]) +
        (x.v[j+2] * Y_column[j+2]) +
        (x.v[j+3] * Y_column[j+3]) +
        (x.v[j+4] * Y_column[j+4]);
    }
    for ( ; j<row; j++)
      sum += (x.v[j] * Y_column[j]);

    a->v[i] = sum;
  }
}

void VectorXMatricTest(void)
{
  Vector x(3);
  Matrix Y(3, 5);
  Vector a;

  for (int i=0; i<3; i++)
  {
    for (int j=0; j<5; j++)
    {
      Y.At(i, j) = i*5+j+1;
    }
  }

  for (int j=0; j<3; j++)
  {
    x.v[j] = 1.0;
  }

  VectorXMatrix(x, Y, &a);
  EXPECT_DOUBLE_EQ_HI(a.v[0], 18.0);
  EXPECT_DOUBLE_EQ_HI(a.v[1], 21.0);
  EXPECT_DOUBLE_EQ_HI(a.v[2], 24.0);
  EXPECT_DOUBLE_EQ_HI(a.v[3], 27.0);
  EXPECT_DOUBLE_EQ_HI(a.v[4], 30.0);
}
