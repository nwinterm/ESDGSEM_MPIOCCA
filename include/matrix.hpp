#ifndef RASCALS_MATRIX
#define RASCALS_MATRIX


template <class T> class matrix;
#include "headers3d.hpp"
#include "blaslapack.hpp"

template <class T>
class matrix {
private:
  int Nrows;
  int Ncolumns;
  T *data;

  void allocate();
public:

  matrix();
  matrix(int n);
  matrix(int nr, int nc);

  // Copy constructor
  matrix(const matrix <T> &A);

  // Assignment operator
  matrix <T> operator = (const matrix <T> &A);

  // Assignment operator from double **array (assumes initialized)
  matrix <T> operator = (const double *A);

  // Assignment operator (all values set to d)
  matrix <T> operator = (const double &d);

  // Assignment operator (all values set to i)
  matrix <T> operator = (const int &i);

  matrix<T> operator += (const int &d);
  matrix<T> operator -= (const int &d);

  ~matrix();

  int nrows() const;
  int ncolumns() const;
  int size() const;

  T *c_array();

  void resize(int n);
  void resize(int nr, int nc);

  void reshape(int nr, int nc);

  matrix <T> transpose();

  // Column major serial access - 1-indexed
  T operator [] (int r) const;
  T & operator [] (int r);

  T operator () (int r, int c) const;
  T & operator () (int r, int c);
  T & operator () (int r);

  // Permute
  matrix<T> operator [] (const matrix <int> &ind);

  // Member function to randomize entries
  void randomize();

  // Sort columns using user supplied comparison function
  void sort( int (*compare)(const void *, const void *) );

  // Not implemented for general T
  void symeig(matrix <T> &d, matrix <T> &v);

  // Not implemented for general T
  // Beware assume eigenvecs are real
  void eig(matrix <T> &WR, matrix <T> &WI, matrix <T> &VL, matrix <T> &VR);

  int byteCount();

  long double frobenius();

  T maxEntry();
  T minEntry();
};

template <class T>
ostream & operator  << (ostream &os, matrix <T> & A);

template <class T>
matrix <T> operator * (const matrix <T> & A, const matrix <T> &B);

template <class T>
matrix <T> operator + (const matrix <T> & A, const matrix <T> &B);

template <class T>
matrix <T> operator - (const matrix <T> & A, const matrix <T> &B);

// General left matrix inverse not implemented
matrix <double> operator | (const matrix <double> & A, const matrix <double> &B);

// General left matrix inverse not implemented
matrix <float> operator | (const matrix <float> & A, const matrix <float> &B);

#include "matrix.tpp"

#endif
