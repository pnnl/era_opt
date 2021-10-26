#include "system_includes.h"
#include "size_file.h"
int64_t size_file(char *filename, FILE *lfp) {
  /*
    Determine the size of a file 
    Called by: 
    Calls:     stat, fprintf, fflush
  */
  struct stat stat_buff;
  int64_t file_len;
  int     stat_ret;
  int     save_err;
  
  stat_ret = stat(filename, &stat_buff);
  save_err = errno;
  if (stat_ret != 0) {
    file_len = (int64_t)0;
    if (lfp) {
      fprintf(lfp,
	      "size_file: Error from stat call on %s was %s\n",
	      filename,strerror(save_err));
      fflush(lfp);
    }
  } else {
    file_len = stat_buff.st_size;
  }
  return(file_len);
}

