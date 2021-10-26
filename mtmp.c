#include "system_includes.h"
#include "mtmp.h"
int mtmp(int num_a_rows, int num_a_cols, double *a,
	 int num_b_rows, int num_b_cols, double *b,
	 double *atb, FILE *lfp) {
  /*
    inner product formulation of matrix transpose * matrix product: a^t * b
    num_a_rows must equal num_b_rows, all matrices stored in column
    major order.
    Called by: prep
  */
 
  double sum;
  double *atb_pos;
  double *a_col;
  double *b_col;
  
  int i;
  int j;
  int k;
  int success;
  success = 1;
  if (num_a_rows != num_b_rows) {
    if (lfp) {
      fprintf(lfp,"mtm error: num_a_rows must match num_b_rows\n");
      fflush(lfp);
      success =0;
    }
  }
  if (success) {
    atb_pos = atb;
    b_col   = b;
    for (j=0;j<num_b_cols;j++) {
      a_col   = a;
      for (i=0;i<num_a_cols;i++) {
	sum = 0.0;
	for (k=0;k<num_a_rows;k++) {
	  sum = sum + (a_col[k] * b_col[k]);
	}
	*atb_pos = sum;
	atb_pos  = atb_pos + 1; /* Caution address arithmetic */
	a_col    = a_col + num_a_rows; /* Caution address arithmetic */
      }
      b_col = b_col + num_b_rows; /* Caution address arithmetic */
    }
  }
  return(success);
}
    
