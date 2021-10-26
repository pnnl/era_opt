#include "system_includes.h"
#include "state_struct.h"
#include "size_file.h"
#include "mmap_coord_files.h"
int mmap_coord_files(struct state_struct *state) {
  /*
    memory map the input coordinate files.
    Called by: era_opt
    Calls:     size_file,fopen,fileno,mmap,strerror,fprintf,fflush
  */
  char *z_dom_file;
  char *y_dom_file;
  char *x_dom_file;
  char *combo_dir;
  char *dir_plus_fn;
  char *file_pos;
  
  double *z_dom_data;
  double *y_dom_data;
  double *x_dom_data;
  double *d_z_dom_base;
  double *d_y_dom_base;
  double *d_x_dom_base;
  double *r_ptr;

  int64_t *iy_dom_bins;
  int64_t *ix_dom_bins;
  int64_t *j_ptr;

  void *null_addr;
  void *x_dom_addr;
  void *y_dom_addr;
  void *z_dom_addr;

  int  *iz_dom_meta_data;
  int  *iy_dom_meta_data;
  int  *ix_dom_meta_data;

  double y_dom_combo_bin_width;
  double x_dom_combo_bin_width;

  double y_dom_vb_min;
  double y_dom_vb_max;
  double y_dom_ub_min;
  double y_dom_ub_max;

  double x_dom_vb_min;
  double x_dom_vb_max;
  double x_dom_ub_min;
  double x_dom_ub_max;

  int64_t offset;
  int64_t z_dom_size;
  int64_t y_dom_size;
  int64_t x_dom_size;
  int64_t j0;

  int64_t izc;
  int64_t iyc;
  int64_t ixc;

  int nz_dom_fields;
  int iz_dom_tsv_offset;
  
  int iz_dom_sbpi;
  int ny_dom_fields;
  
  int iy_dom_tsv_offset;
  int iy_dom_sbpi;

  int nx_dom_fields;
  int ix_dom_tsv_offset;

  int ix_dom_sbpi;
  int success;

  int prot;
  int flags;

  int xfd;
  int yfd;

  int zfd;
  int mask;

  int iz_dom_ner;
  int iy_dom_ner;
  
  int ix_dom_ner;
  int iz_dom_data_offset;

  int iy_dom_data_offset;
  int ix_dom_data_offset;

  int save_err;
  int extra_log_info;

  int nx_dom_bins;
  int ny_dom_bins;

  int icombo_dir_len;
  int ipad;

  FILE *lfp;
  FILE *xfp;
  FILE *yfp;
  FILE *zfp;

  success        = 1;
  lfp            = state->lfp;
  z_dom_file     = state->z_dom_file;
  y_dom_file     = state->y_dom_file;
  x_dom_file     = state->x_dom_file;
  combo_dir      = state->combo_dir;
  dir_plus_fn    = state->dir_plus_fn;
  extra_log_info = state->extra_log_info;
  icombo_dir_len = strlen(combo_dir);
  /*
    prepend the combo_dir to the z_dom_file name for
    passing to size_file.
  */
  file_pos       = dir_plus_fn;
  if (icombo_dir_len > 0) {
    strcpy(dir_plus_fn,combo_dir);
    file_pos += icombo_dir_len; /* Caution address arithmetic */
  }
  strcpy(file_pos,z_dom_file);
  /*
  z_dom_size     = size_file(z_dom_file,lfp);
  */
  z_dom_size     = size_file(dir_plus_fn,lfp);
  j0 = (int64_t)0;
  if (z_dom_size == j0) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"mmap_color_files: Error %s was of zero length \n",
	      z_dom_file);
      fflush(lfp);
    }
  }
  if (success) {
    /*
      Open the z-dominant(blue) file.
    zfp = fopen(z_dom_file,"r");
    */
    zfp = fopen(dir_plus_fn,"r");
    if (zfp == NULL) {
      success = 0;
      if (lfp) {
	/*
	fprintf(lfp,"mmap_color_files: Unable to open %s\n",z_dom_file);
	*/
	fprintf(lfp,"mmap_color_files: Unable to open %s\n",dir_plus_fn);
	fflush(lfp);
      }
    }
  }
  if (success) {
    /*
      memory map the z-dominant(blue) file.
    */
    zfd       = fileno(zfp);
    null_addr = NULL;
    offset    = (int64_t)0;
    prot      = PROT_READ;
    flags     = MAP_PRIVATE;
    z_dom_addr = mmap(null_addr,z_dom_size,prot,flags,zfd,offset);
    save_err  = errno;
    if (z_dom_addr == MAP_FAILED) {
      success = 0;
      if (lfp) {
	/*
	fprintf(lfp,
	"mmap_color_files: Error, unable to mmap %s, error = %d, %s \n",
		z_dom_file,save_err,strerror(save_err));
	*/
	fprintf(lfp,
	"mmap_color_files: Error, unable to mmap %s, error = %d, %s \n",
		dir_plus_fn,save_err,strerror(save_err));
	fflush(lfp);
      }
    }
  }
  if (success) {
    /*
      Here we reuse that fact that we have already copied
      combo_dir to the start of dir_plus_fn variable
      and file_pos still points to the next character after the directory.
    */
    strcpy(file_pos,y_dom_file);
    /*
    y_dom_size  = size_file(y_dom_file,lfp);
    */
    y_dom_size  = size_file(dir_plus_fn,lfp);
    if (y_dom_size == j0) {
      if (lfp) {
	/*
	fprintf(lfp,"map_color_files: Error %s was of zero length \n",
		y_dom_file);
	*/
	fprintf(lfp,"map_color_files: Error %s was of zero length \n",
		dir_plus_fn);
	fflush(lfp);
      }
    }
  }
  if (success) {
    /*
      Open the y-dominant file.
    yfp = fopen(y_dom_file,"r");
    */
    yfp = fopen(dir_plus_fn,"r");
    if (yfp == NULL) {
      success = 0;
      if (lfp) {
	/*
	fprintf(lfp,"mmap_color_files: Unable to open %s\n",y_dom_file);
	*/
	fprintf(lfp,"mmap_color_files: Unable to open %s\n",dir_plus_fn);
	fflush(lfp);
      }
    }
  }
  if (success) {
    /*
      memory map the y-dominant file.
    */
    yfd        = fileno(yfp);
    null_addr  = NULL;
    offset     = (int64_t)0;
    prot       = PROT_READ;
    flags      = MAP_PRIVATE;
    y_dom_addr = mmap(null_addr,y_dom_size,prot,flags,yfd,offset);
    save_err   = errno;
    if (y_dom_addr == MAP_FAILED) {
      success = 0;
      if (lfp) {
	/*
	fprintf(lfp,
	    "mmap_color_files: Error, unable to mmap %s, error = %d, %s \n",
		y_dom_file,save_err,strerror(save_err));
	*/
	fprintf(lfp,
	    "mmap_color_files: Error, unable to mmap %s, error = %d, %s \n",
		dir_plus_fn,save_err,strerror(save_err));
	fflush(lfp);
      }
    }
  }
  if (success) {
    /*
      Here we reuse that fact that we have already copied
      combo_dir to the start of dir_plus_fn variable
      and file_pos still points to the next character after the directory.
    */
    strcpy(file_pos,x_dom_file);
    /*
    x_dom_size  = size_file(x_dom_file,lfp);
    */
    x_dom_size  = size_file(dir_plus_fn,lfp);
    if (x_dom_size == j0) {
      if (lfp) {
	/*
	fprintf(lfp,"map_color_files: Error %s was of zero length \n",
		x_dom_file);
	*/
	fprintf(lfp,"map_color_files: Error %s was of zero length \n",
		dir_plus_fn);
	fflush(lfp);
      }
    }
  }
  if (success) {
    /*
      Open the x-dominant file.
    */
    /*
    xfp = fopen(x_dom_file,"r");
    */
    xfp = fopen(dir_plus_fn,"r");
    if (xfp == NULL) {
      success = 0;
      if (lfp) {
	/*
	fprintf(lfp,"mmap_color_files: Unable to open %s\n",x_dom_file);
	*/
	fprintf(lfp,"mmap_color_files: Unable to open %s\n",dir_plus_fn);
	fflush(lfp);
      }
    }
  }
  if (success) {
    /*
      memory map the x-dominant file.
    */
    xfd       = fileno(xfp);
    null_addr = NULL;
    offset    = (int64_t)0;
    prot      = PROT_READ;
    flags     = MAP_PRIVATE;
    x_dom_addr  = mmap(null_addr,x_dom_size,prot,flags,xfd,offset);
    save_err  = errno;
    if (z_dom_addr == MAP_FAILED) {
      success = 0;
      if (lfp) {
	/*
	fprintf(lfp,
	"mmap_color_files: Error, unable to mmap %s, error = %d, %s \n",
		x_dom_file,save_err,strerror(save_err));
	*/
	fprintf(lfp,
	"mmap_color_files: Error, unable to mmap %s, error = %d, %s \n",
		dir_plus_fn,save_err,strerror(save_err));
	fflush(lfp);
      }
    }
  }
  if (success) {
    /* Now extract information from the meta data fields and 
       set the data and bin array info.
    */
    mask               = 255;
    iz_dom_meta_data    = (int*)z_dom_addr;
    d_z_dom_base        = (double*)z_dom_addr;
    iz_dom_sbpi         = iz_dom_meta_data[0] & mask;
    iz_dom_tsv_offset   = (iz_dom_meta_data[1] >> 24) & mask;
    iz_dom_tsv_offset   = iz_dom_tsv_offset >> 3; 
    nz_dom_fields       = iz_dom_tsv_offset + 3;
    iz_dom_ner          = iz_dom_meta_data[4];
    j_ptr              = (int64_t*)&iz_dom_meta_data[8];
    izc                = *j_ptr;
    iz_dom_data_offset  = (iz_dom_ner + 3) * nz_dom_fields;
    z_dom_data          = (double *)&d_z_dom_base[iz_dom_data_offset];

    iy_dom_meta_data   =  (int*)y_dom_addr;
    d_y_dom_base       =  (double*)y_dom_addr;
    iy_dom_sbpi        =  iy_dom_meta_data[0] & mask;
    iy_dom_tsv_offset  =  (iy_dom_meta_data[1] >> 24) & mask; 
    iy_dom_tsv_offset  =  iy_dom_tsv_offset >> 3;
    ny_dom_fields      =  iy_dom_tsv_offset + 5;
    iy_dom_ner         =  iy_dom_meta_data[4];
    ny_dom_bins        =  iy_dom_meta_data[5];
    j_ptr              =  (int64_t*)&iy_dom_meta_data[8];
    iyc                =  *j_ptr;
    r_ptr              =  (double*)&iy_dom_meta_data[14];
    y_dom_combo_bin_width = *r_ptr;
    r_ptr              =  (double*)&iy_dom_meta_data[16];
    y_dom_vb_min       =  *r_ptr;
    r_ptr              =  (double*)&iy_dom_meta_data[18];
    y_dom_vb_max       =  *r_ptr;
    r_ptr              =  (double*)&iy_dom_meta_data[20];
    y_dom_ub_min       =  *r_ptr;
    r_ptr              =  (double*)&iy_dom_meta_data[22];
    y_dom_ub_max       =  *r_ptr;
    iy_dom_data_offset =  (iy_dom_ner + 3) * ny_dom_fields;
    y_dom_data         =  (double*)&d_y_dom_base[iy_dom_data_offset];
    iy_dom_bins        =  (int64_t*)&y_dom_data[ny_dom_fields*iyc];

    ix_dom_meta_data     =  (int*)x_dom_addr;
    d_x_dom_base         =  (double*)x_dom_addr;
    ix_dom_sbpi          =  ix_dom_meta_data[0] & mask;
    ix_dom_tsv_offset    =  (ix_dom_meta_data[1] >> 24) & mask; 
    ix_dom_tsv_offset    =  ix_dom_tsv_offset >> 3;
    nx_dom_fields        =  ix_dom_tsv_offset + 5;
    ix_dom_ner           =  ix_dom_meta_data[4];
    nx_dom_bins          =  ix_dom_meta_data[5];
    j_ptr              =  (int64_t*)&ix_dom_meta_data[8];
    ixc                =  *j_ptr;
    r_ptr              =  (double*)&ix_dom_meta_data[14];
    x_dom_combo_bin_width = *r_ptr;
    r_ptr              =  (double*)&ix_dom_meta_data[16];
    x_dom_vb_min         = *r_ptr;
    r_ptr              =  (double*)&ix_dom_meta_data[18];
    x_dom_vb_max         = *r_ptr;
    r_ptr              =  (double*)&ix_dom_meta_data[20];
    x_dom_ub_min         = *r_ptr;
    r_ptr              =  (double*)&ix_dom_meta_data[22];
    x_dom_ub_max         = *r_ptr;
    ix_dom_data_offset   =  (ix_dom_ner + 3) * nx_dom_fields;
    x_dom_data           =  (double*)&d_x_dom_base[ix_dom_data_offset];
    ix_dom_bins          =  (int64_t*)&x_dom_data[nx_dom_fields*ixc];

    state->iz_dom_sbpi        = iz_dom_sbpi;
    state->iz_dom_tsv_offset  = iz_dom_tsv_offset;
    state->nz_dom_fields      = nz_dom_fields;
    state->izc               = izc;
    state->z_dom_data         = z_dom_data;

    state->iy_dom_sbpi       = iy_dom_sbpi;
    state->iy_dom_tsv_offset = iy_dom_tsv_offset;
    state->ny_dom_fields     = ny_dom_fields;
    state->ny_dom_bins       = ny_dom_bins;
    state->iyc               = iyc;
    state->y_dom_data        = y_dom_data;
    state->iy_dom_bins       = iy_dom_bins;
    state->y_dom_combo_bin_width  = y_dom_combo_bin_width;
    state->y_dom_vb_min      = y_dom_vb_min;
    state->y_dom_vb_max      = y_dom_vb_max;
    state->y_dom_ub_min      = y_dom_ub_min;
    state->y_dom_ub_max      = y_dom_ub_max;

    state->ix_dom_sbpi         = ix_dom_sbpi;
    state->ix_dom_tsv_offset   = ix_dom_tsv_offset;
    state->nx_dom_fields       = nx_dom_fields;
    state->ixc               = ixc;
    state->nx_dom_bins         = nx_dom_bins;
    state->x_dom_data          = x_dom_data;
    state->ix_dom_bins         = ix_dom_bins;
    state->x_dom_combo_bin_width = x_dom_combo_bin_width;
    state->x_dom_vb_min        = x_dom_vb_min;
    state->x_dom_vb_max        = x_dom_vb_max;
    state->x_dom_ub_min        = x_dom_ub_min;
    state->x_dom_ub_max        = x_dom_ub_max;
    if (extra_log_info) {
      if (lfp) {
	fprintf (lfp,"iz_dom_sbpi       	= %d\n",iz_dom_sbpi);
	fprintf (lfp,"iz_dom_tsv_offset 	= %d\n",iz_dom_tsv_offset);
	fprintf (lfp,"nz_dom_fields     	= %d\n",nz_dom_fields);
	fprintf (lfp,"izc              	        = %ld\n",izc);
	fprintf (lfp,"z_dom_data        	= %p\n",z_dom_data);
	fprintf (lfp,"iy_dom_sbpi      	= %d\n",iy_dom_sbpi);
	fprintf (lfp,"iy_dom_tsv_offset	= %d\n",iy_dom_tsv_offset);
	fprintf (lfp,"ny_dom_fields   	= %d\n",ny_dom_fields);
	fprintf (lfp,"iyc              	= %ld\n",iyc);
	fprintf (lfp,"y_dom_combo_bin_width = %le\n",y_dom_combo_bin_width);
	fprintf (lfp,"y_dom_data       	= %p\n",y_dom_data);
	fprintf (lfp,"iy_dom_bins      	= %p\n",iy_dom_bins);
	fprintf (lfp,"ix_dom_sbpi      	= %d\n",ix_dom_sbpi);
	fprintf (lfp,"ix_dom_tsv_offset	= %d\n",ix_dom_tsv_offset);
	fprintf (lfp,"nx_dom_fields   	= %d\n",nx_dom_fields);
	fprintf (lfp,"ixc              	= %ld\n",ixc);
	fprintf (lfp,"x_dom_combo_bin_width = %le\n",x_dom_combo_bin_width);
	fprintf (lfp,"x_dom_data       	= %p\n",x_dom_data);
	fprintf (lfp,"ix_dom_bins      	= %p\n",ix_dom_bins);
      }
    }
  }
  return (success);
}
