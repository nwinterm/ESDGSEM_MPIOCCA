template <class T>
matrix<T>::matrix(){
  data = NULL;
  Nrows = Ncolumns = 0;
}

template <class T>
matrix<T>::matrix(int n){
  Nrows    = n;
  Ncolumns = 1;

  allocate();
}

template <class T>
matrix<T>::matrix(int nr, int nc){
  Nrows    = nr;
  Ncolumns = nc;

  allocate();
}

// Copy constructor
template <class T>
matrix<T>::matrix(const matrix<T> &A){
  Nrows    = A.Nrows;
  Ncolumns = A.Ncolumns;

  allocate();

  for(int r=1; r<=Nrows; r=r+1)
    for(int c=1; c<=Ncolumns; c=c+1)
      (*this)(r,c) = A(r,c);
}

// Assignment operator
template <class T>
matrix<T> matrix<T>::operator = (const matrix<T> &A){
  if(this != &A){
    if(Nrows != A.Nrows && Ncolumns != A.Ncolumns){

      if(data)
	delete [] data;

      Nrows = A.Nrows;
      Ncolumns = A.Ncolumns;

      allocate();
    }

    for(int r=1; r<=Nrows; r=r+1)
      for(int c=1; c<=Ncolumns; c=c+1)
	(*this)(r,c) = A(r,c);
  }

  return *this;
}


// Assignment operator from double **array (assumes initialized)
template <class T>
matrix<T> matrix<T>::operator = (const double *A){

  // Load from C style matrix
  for(int r=1; r<=Nrows; r=r+1)
    for(int c=1; c<=Ncolumns; c=c+1)
      (*this)(r,c) = A[Ncolumns*(r-1) + c-1];

  return *this;
}

// Assignment operator (all values set to d)
template <class T>
matrix<T> matrix<T>::operator = (const double &d){
  for(int r=1; r<=Nrows; r=r+1)
    for(int c=1; c<=Ncolumns; c=c+1)
      (*this)(r,c) = d;

  return *this;
}

// Assignment operator (all values set to i)
template <class T>
matrix<T> matrix<T>::operator = (const int &i){
  for(int r=1; r<=Nrows; r=r+1)
    for(int c=1; c<=Ncolumns; c=c+1)
      (*this)(r,c) = i;

  return *this;
}

// Increment operator 
template <class T>
matrix<T> matrix<T>::operator += (const int &d){
  for(int r=1; r<=Nrows; r=r+1)
    for(int c=1; c<=Ncolumns; c=c+1)
      (*this)(r,c) += d;

  return *this;
}

// Increment operator 
template <class T>
matrix<T> matrix<T>::operator -= (const int &d){
  for(int r=1; r<=Nrows; r=r+1)
    for(int c=1; c<=Ncolumns; c=c+1)
      (*this)(r,c) -= d;

  return *this;
}


template <class T>
matrix<T>::~matrix(){
  if(data)
    delete [] data;
  data = NULL;   
}

template <class T>
int matrix<T>::nrows() const{
  return Nrows;
}

template <class T>
int matrix<T>::ncolumns() const{
  return Ncolumns;
}

template <class T>
int matrix<T>::size() const{
  return Nrows*Ncolumns;
}

template <class T>
T* matrix<T>::c_array(){
  return data;
}

template <class T>
void matrix<T>::allocate(){
  assert( Nrows*Ncolumns );

  data = new T[Nrows*Ncolumns]();

}

template <class T>
void matrix<T>::resize(int n){
  resize(n,1);
}

template <class T>
void matrix<T>::resize(int nr, int nc){
  // Temporary
  matrix<T> dupe(nr,nc);
  dupe = 0.0;

  if(data){
    // Entry by entry copy to dupe
    for(int r=1; r<=nr; r=r+1)
      for(int c=1; c<=nc; c=c+1)
	if(r<=Nrows && c<=Ncolumns)
	  dupe(r,c) = (*this)(r,c);

    delete [] data;
  }

  Nrows    = nr;
  Ncolumns = nc;

  allocate();

  for(int r=1; r<=nr; r=r+1)
    for(int c=1; c<=nc; c=c+1)
      (*this)(r,c) = dupe(r,c);
}

template <class T>
void matrix<T>::reshape(int nr, int nc){
  assert( nr*nc == Nrows*Ncolumns );

  Nrows    = nr;
  Ncolumns = nc;
}

template <class T>
matrix<T> matrix<T>::transpose(){
  matrix<T> AT(Ncolumns, Nrows);

  for(int r=1; r<=Nrows; r=r+1)
    for(int c=1; c<=Ncolumns; c=c+1)
      AT(c,r) = (*this)(r,c);

  return AT;
}

// Column major serial access - 1-indexed
template <class T>
T matrix<T>::operator [] (int r) const {
  assert( (0 < r) && (r <= Nrows*Ncolumns) );

  return data[r-1];
}

// Column major serial access - 1-indexed
template <class T>
T& matrix<T>::operator [] (int r) {
  assert( (0 < r) && (r <= Nrows*Ncolumns) );

  return data[r-1];
}

template <class T>
T matrix<T>::operator () (int r, int c) const {
  // Bounds check here [ could also throw an exception ]
  assert( (0 < r) && (r <= Nrows) && (0 < c) && (c <= Ncolumns) );

  return data[ (r-1) + (c-1)*Nrows ];
}

template <class T>
T& matrix<T>::operator () (int r, int c) {
  assert( (0 < r) && (r <= Nrows) && (0 < c) && (c <= Ncolumns) );

  return data[ (r-1) + (c-1)*Nrows ];
}

template <class T>
T& matrix<T>::operator () (int r) {
  assert( (0 < r) && (r <= Nrows*Ncolumns) );

  return data[r-1];
}

// Permute
template <class T>
matrix<T> matrix<T>::operator [] (const matrix <int> &ind)  {
  matrix<T> C(ind.nrows(), ind.ncolumns());

  for(int r=1;r<=ind.nrows();r=r+1){
    for(int c=1;c<=ind.ncolumns();c=c+1){
      C(r,c) = (*this)[ind(r,c)];
    }
  }

  return C;
}

// Member function to randomize entries
template <class T>
void matrix<T>::randomize(){
  int r, c;
  for(c=1;c<=Ncolumns;++c){
    for(r=1;r<=Nrows;++r){
      (*this)(r,c) = drand48();
    }
  }
}

// Sort columns using user supplied comparison function
template <class T>
void matrix<T>::sort( int (*compare)(const void *, const void *) ){
  qsort(data, Ncolumns, Nrows*sizeof(T), compare);
}

// Not implemented for general T
template <class T>
void matrix<T>::symeig(matrix<T> &d, matrix<T> &v){}

// Not implemented for general T
// Beware assume eigenvecs are real
template <class T>
void matrix<T>::eig(matrix<T> &WR, matrix<T> &WI, matrix<T> &VL, matrix<T> &VR){}

template <>
void matrix<float>::eig(matrix<float> &WR, matrix<float> &WI, matrix<float> &VL, matrix<float> &VR);

template <>
void matrix<double>::eig(matrix<double> &WR, matrix<double> &WI, matrix<double> &VL, matrix<double> &VR);

template <class T>
int matrix<T>::byteCount(){
  return sizeof(T)*Nrows*Ncolumns;
}

template <class T>
long double matrix<T>::frobenius(){
  long double f(0);

  for(int c=1;c<=Ncolumns;++c)
    for(int r=1;r<=Nrows;++r)
      f += (*this)(r,c)*(*this)(r,c);

  return f;
}

template <class T>
T matrix<T>::maxEntry(){
  T maxval = (*this)(1,1);

  for(int c=1;c<=Ncolumns;++c){
    for(int r=1;r<=Nrows;++r){
      T val = (*this)(r,c);
      maxval = max(maxval, val);
    }
  }

  return maxval;
}

template <class T>
T matrix<T>::minEntry(){
  T minval = (*this)(1,1);

  for(int c=1;c<=Ncolumns;++c){
    for(int r=1;r<=Nrows;++r){
      T val = (*this)(r,c);
      minval = min(minval, val);
    }
  }

  return minval;
}

template <class T>
ostream & operator << (ostream &os, matrix<T> & A){
  int r,c;

  os << "[";
  for(r=1;r<=A.nrows();r=r+1){

    os << "[";

    for(c=1;c<A.ncolumns();c=c+1)
      os << A(r,c) << ",";

    os << A(r,c) << "]";

    if(r<A.nrows())
      os << endl;
  }

  os << "];" << endl;

  return os;
}

template <class T>
matrix<T> operator * (const matrix<T> & A, const matrix<T> &B){
  matrix<T> C(A.nrows(), B.ncolumns());

  for(int r=1;r<=A.nrows();r=r+1){
    for(int c=1;c<=B.ncolumns();c=c+1){
      // Carefully use this to start row-column product
      T s = A(r,1)*B(1,c);
      for(int i=2;i<=A.ncolumns();i=i+1)
	s = s + A(r,i)*B(i,c);
      C(r,c) = s;
    }
  }

  return C;
}

template <class T>
matrix<T> operator + (const matrix<T> & A, const matrix<T> &B){
  matrix<T> C = A;

  for(int r=1;r<=A.nrows();r=r+1)
    for(int c=1;c<=B.ncolumns();c=c+1)
      C(r,c) = C(r,c)+B(r,c);

  return C;
}

template <class T>
matrix<T> operator - (const matrix<T> & A, const matrix<T> &B){
  matrix<T> C = A;

  for(int r=1;r<=A.nrows();r=r+1)
    for(int c=1;c<=B.ncolumns();c=c+1)
      C(r,c) = C(r,c)-B(r,c);

  return C;
}
