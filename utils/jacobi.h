int jacobi(int n, float **a, float *d, float **v);
  /*  computes eigenvalues and -vectors from matrix a of dimension
   *  N. Elements of A above diagonal are destroyed. Eigenvectors are
   *  returned in array D, and eigenvectors are the columns of
   *  V. Eigenvalues and -vectors are sorted to descending order. The
   *  number of rotations needed is returned 
   */ 
