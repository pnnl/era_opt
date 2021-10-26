#ifndef _CSV_READ_
#define _CSV_READ_ 1
extern int csv_read(char *filename, char *io_buffer,
		    int buffer_len,
		    int nrows, int ncols, double *values,
		    FILE *lfp);
#endif
