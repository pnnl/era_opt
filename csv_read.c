#include "system_includes.h"
#include "field_width.h"
#include "csv_read.h"
int csv_read(char *filename, char *io_buffer,
	     int buffer_len,
	     int nrows, int ncols, double *values,
	     FILE *lfp) {
  /*
    Read a comma separated values file given by filename, with 
    nrows rows and ncols columns. Return 1 on success 0 on failure.
    values is assumed to be an array of length nrows*ncols.
    io_buffer is of lengthe buffer_len and
    assumed to be long enough to hold any row.
    Called by: thous_opt_2
    Calls:     fopen, fclose, fgets, sscanf, field_width
  */
  FILE *csv_fp;
  char *next_val;
  double *value;
  int success;
  int ns;

  int left;
  int ifw;

  int line_no;
  int j;
  
  success = 1;
  csv_fp = fopen(filename,"r");
  if (csv_fp == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"csv_read: could not open %s for reading\n",
	      filename);
      fflush(lfp);
    }
  }
  if (success) {
    value = values;
    for (line_no = 1;((line_no <= nrows) && success);line_no++) {
      next_val = fgets(io_buffer,buffer_len,csv_fp);
      if (next_val != io_buffer) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"csv_read: Error reading line %d of %s\n",
		  line_no,filename);
	  fflush(lfp);
	}
      } else {
	left = strlen(next_val);
	for (j = 1;((j<=ncols) && (left > 0) && success);j++) {
	  ifw = field_width(next_val,left);
	  /* replace field ending , or lf with null for sscan f */
	  next_val[ifw] = '\0';
	  if (ifw > 0) {
	    ns = sscanf(next_val,"%le",value);
	    if (ns != 1) { 
	      success = 0;
	      if (lfp) {
		fprintf(lfp,"csv_read error reading column %d in row %d\n%s\n",
			j,line_no,io_buffer);
	      }
	    }
	  } else {
	    *value = 0.0;
	  }
	  left -= (ifw + 1);
	  next_val += (ifw + 1); /* Caution address arithmetic. */
	  value += 1; /* Caution address arithmetic. */
	} /* end for (j...) */
      } /* else succeeded in reading line */
    } /* for line_no */
    fclose(csv_fp);
  } /* end if file open succeeded */
  return (success);
}
