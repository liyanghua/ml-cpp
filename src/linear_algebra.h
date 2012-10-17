/** @file
* @brief some basic linear algebra operations
* @author zhangyafeikimi@gmail.com
* @date
* @version
*
*/
#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

struct Vector;
struct Matrix;

// a column-based vector
struct Vector
{
  double * v;
  int l;

  Vector() :v(0), l(0) {}
  explicit Vector(int _l) :l(_l)
  {
    v = new double[l];
  }
  ~Vector()
  {
    delete [] v;
  }
  void CopyFrom(const Vector& ve);
  void Swap(Vector * ve);


  void Init(int _l);
  // zero all elements
  void Zero(void);
  // a = alpha * a
  void Scale(double alpha);
  // a = alpha * y
  void AlphaY(double alpha, const Vector& y);
  // a = a + alpha * y
  void PlusAlphaY(double alpha, const Vector& y);
  // a = x + alpha * y
  void XPlusAlphaY(const Vector& x, double alpha, const Vector& y);


  // inner product
  // return a^T * y
  double Dot(const Vector& y)const;
  // return a^T * a
  double DotSelf(void)const;
  // return the sum of all elements
  double Sum(void)const;
  // return 1-norm
  double Norm1(void)const;
  // return 2-norm
  double Norm2(void)const;

private:
  // we can copy Vectors by CopyFrom, but any implicit copy are disabled.
  Vector(const Vector&);
  Vector& operator=(const Vector&);
};


// a m*n matrix
struct Matrix
{
  double * v;
  int m, n;

  Matrix() :v(0), m(0), n(0) {}
  Matrix(int _m, int _n) :m(_m), n(_n)
  {
    v = new double[m*n];
  }
  ~Matrix()
  {
    delete [] v;
  }

  // accessors
  inline double& At(int i, int j)
  {
    return v[i+m*j];
  }
  // const accessors
  inline const double& At(int i, int j)const
  {
    return v[i+m*j];
  }


  // zero all elements
  void Zero(void);

private:
  // Any implicit copy are disabled.
  Matrix(const Matrix&);
  Matrix& operator=(const Matrix&);
};


// TODO sparse matrix


// a = X * y
// X: m*n
// y: n*1
// a: m*1
void MatrixXVector(const Matrix& X, const Vector& y, Vector * a);

// a^T = x^T * Y
// x: m*1
// Y: m*n
// a: 1*n
void VectorXMatrix(const Vector& x, const Matrix& Y, Vector * a);

#endif
