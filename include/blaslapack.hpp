extern "C" {
  // General linear algebra routines
  void sgemm_ (char *TRANSA, char *TRANSB, int *M, int *N, int *K, float *ALPHA , float *A , int *LDA, float *B , int *LDB, float *BETA , float *C , int *LDC);
  void dgemm_ (char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *A, int *LDA, double *B, int *LDB, double *BETA, double *C, int *LDC);

  // Solves AX = B (for X)
  void sgesv_ (int *N, int *NRHS, float *A , int *LDA, int *IPIV, float *B , int *LDB, int *INFO );
  void dgesv_ (int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );

  // Calculates eigenvalues/eigenvectors of A (symmetric)
  void ssyev_(char *JOBZ, char *UPLO, int *N, float *A , int *LDA, float *W , float *WORK , int *LWORK, int *INFO );
  void dsyev_(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO );

  // Calculates eigenvalues/eigenvectors of A (symmetric)
  void sgeev_(char *JOBVL, char *JOBVR, int *N, float *A , int *LDA, float *WR , float *WI ,
	      float *VL , int *LDVL, float *VR , int *LDVR, float *WORK , int *LWORK, int *INFO );
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
	      double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
}
