#include "system_includes.h"
#include "matrix_transpose.h"
void matrix_transpose(int nrows, int ncols, double *matrix, double *matrix_t) {
  /*
    transpose a nrows x ncols matrix in row major order into a column major
    order.
    Called by: prep
  */
  double *column;
  double *matrix_p;
  int i;
  int j;
  column = matrix;
  matrix_p = matrix_t;
  /*
    Loop accross columns of matrix.
    Extracting one column at a time.
  */
  for (j=0;j<ncols;j++) {
    column = &matrix[j];
    /*
      Loop over rows; Note the writing to consecutive memory locations
      while reading from strided locations.
    */
    for (i=0;i<nrows;i++) {
      *matrix_p = *column;
      matrix_p += 1; /* Caution, address arithmetic */
      column   += ncols; /* Caution, address arithmetic */
    }
  }
}
 
