#include "system_includes.h"
#include "state_struct.h"
#include "pack_imap.h"
#include "init_z_dom_combos.h"
int init_z_dom_combos(struct state_struct *state, 
		      double **ddata_p, 
		      int **imeta_data_p,
		      int *zfd_p,
		      int *tsv_offset_p) {
  
  /*
    Initialize the *.z.dat file and meta data.
    Called by: gen_z_dom_combos
    Calls:     mmap,pack_imap,fopen,fprintf,fflush,fileno,ftruncate,strerror

    The meta-data header for a combos file will have the
    following format:

byte  Description
      First 4 byte integer
0     Length in bytes of this meta data ((ner+3)*bpe) <= 255
1     Number of channels, nc, of this color <= 255
2     Number of active bits/channel_index <= 255, abpi
      (7 for .01 power levels 10 for for .001 power levels. - number
       of bits it takes to represent max_index)

3     Number of stored bits/channel_index, sbpi 
      (8 for .01 power levels 16 for .001 power levels).
      sbpi = (abpi + 7) & 248
      					   
      Second 4 byte integer
4     Total bytes per index, bpi (8) (#3 * #1)/8 rounded up to a multiple of 8.
      tbpi  = ((sbpi/8) * nc) + 7) & 248

5     Bytes per ts value         (8) bptsv
6     Total bytes per entry. bpe = 3*bptsv + tbpi (>= 32)
7     dominant coordinte indicator (0 = x, 1 = y, 2 = z)

8     nx <= 255
9     ny <= 255
10    nz <= 255
11    nx + ny + nz <= 255

      Fourth four byte integer, max_index
12:15 max_index

      Fifth four byte integer, ner
16:19 Number of extra records to store mapping,ner: = 0 if 3bpe-64 <= nc
      else ner = ceil (nc/bpe)

20:23 Sixth 4 byte integer, 
      Number of bins for sorted tristimulus values, num_bins      

24:31 total Number of bytes as an 8 byte integer including these first ner+3
      meta records.

32:39 Number of combinations stored as an 8 byte integer 

40:47 Position in 8 byte words of first index element.

48:55 8byte double, min_fract = minimum power step value (1/max_index, e.g. .01)

56:64 8byte double combo_bin_width.

64:3bpe-1 imap for channels represented if ner = 0, else unused.
3bpe:(3+ner)*bpe imap characters for channels if ner > 0

*/
  struct loop_state_struct *loop_state;
  void 	  *null_addr;
  void 	  *md_addr;
  double  *dmd_p;
  double  *ddata;
  int64_t *jmd_p;
  int  	  *imeta_data;
  int     *imap_sort;
  int     *cimap;
  char 	  *z_dom_file;
  double  dmd;
  double  min_fract;
  int64_t max_elems;
  int64_t zfd_size;
  int64_t offset;
  int64_t jmd;
  

  int 	  zfd;
  int 	  nc;

  int 	  max_index;
  int 	  abpi;

  int 	  sbpi;
  int 	  tbpi;

  int 	  bptsv;
  int 	  bytes_pi;

  int 	  limp1;
  int 	  mask;

  int 	  bpe;
  int 	  step;

  int 	  flags;
  int 	  prot;

  int 	  md_len;
  int 	  x_dom;

  int 	  y_dom;
  int 	  z_dom;

  int     success;
  int     nx;
  
  int     ny;
  int     nz;
  
  int     sum_nxyz;
  int     imd0;

  int     imd1;
  int     imd2;

  int     ner;
  int     i;

  int     mdl;
  int     color;

  int     save_err;
  int     ierr;

  FILE    *zfp;
  FILE    *lfp;

  success       = 1;
  x_dom         = 0;
  y_dom         = 1;
  z_dom         = 2;
  loop_state    = state->loop_state;
  z_dom_file    = state->z_dom_file;
  nx            = loop_state->nx;
  ny            = loop_state->ny;
  nz            = loop_state->nz;
  imap_sort     = loop_state->imap_sort;
  min_fract     = loop_state->min_fract;
  lfp           = state->lfp;

  zfp = fopen(z_dom_file,"w+");
  if (zfp == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"init_z_dom_combos: Error, could not open %s\n",
	      z_dom_file);
      fflush(lfp);
    }
  }
  if (success) {
    state->zfp = zfp;
    zfd = fileno(zfp);
    state->zfd = zfd;
    *zfd_p = zfd;
    /*
      Number of z dominiant channels.
    */
    nc = loop_state->nz;
    max_index = loop_state->max_index;
    /*
      Want ceil(log_2(max_index) rounded up to a multiple of 8.
      = bits/index. abpi is the number of bits needed to represent max_index.
    */
    abpi = 0;
    step = 1;
    while (step <= max_index) {
      step = step + step;
      abpi += 1;
    }
    mask = 248;
    /*
      spbpi is abpi rounded up to a multiple of 8 bits.
    */
    sbpi = (abpi+7) & mask;
    /*
      tbpi is the total bytes per nc indices (1 combination)
    */
    tbpi = (((sbpi/8) * nc) + 7) & mask;
    *tsv_offset_p = tbpi >> 3;
    /*
      Bytes per tristimulus value.
    */
    bptsv = 8;
    /*
      Bytes per entree, bpe = 3 * size of tristimulus values plus index size.
    */
    bpe   = (3 * bptsv) + tbpi;
    /*
      Check to see if there is sufficient space to store imap in remaining
      bytes in records 3 and 4 of the metadata file.
    */
    if (((3*bpe) - 64) >= nc) {
      ner = 0;
    } else {
      ner = (nc + bpe - 1)/bpe;
    }
    md_len = (ner + 3) * bpe;
    
    max_elems = 1;
    for (i=0;i<nz;i++) {
      max_elems = max_elems * (max_index + 1);
    }
    /*
      Actually allocate scratch space for the sorting of the combinations
      and the meta data.
      memory map the file.
    */
    zfd_size  = 2*(max_elems * bpe) + md_len;
    ierr      = ftruncate(zfd,zfd_size);
    save_err  = errno;
    if (ierr < 0) {
      success = 0;
      if (lfp) {
	fprintf(lfp,
		"init_z_dom_combos: Error unable to truncate"
		" %s to %ld bytes, err = %d, %s\n",
		z_dom_file,zfd_size,save_err,strerror(save_err));
	fflush(lfp);
      }
    }
  }
  if (success) {
    null_addr = NULL;
    offset    = (int64_t)0;
    prot      = PROT_READ | PROT_WRITE;
    flags     = MAP_SHARED;
    md_addr   = mmap(null_addr,zfd_size,prot,flags,zfd,offset);
    save_err  = errno;
    if (md_addr == MAP_FAILED) {
      success = 0;
      if (lfp) {
	fprintf(lfp,
		"init_z_dom_combos: Error unable to mmap %s, error = %d, %s\n",
		z_dom_file,save_err,strerror(save_err));
	fflush(lfp);
      }
    }
  }
  if (success) {
    imeta_data = (int*)md_addr;
    mdl = (ner+3)*bpe;
    imd0 = 0;
    mask = 255;
    if (mdl > 255) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"init_z_dom_combos: Error, meta data length = %d > 255\n",
		mdl);
      }
    }
  }
  if (success) {
    /* mdl, nc, abpi,sbpi */
    mask = 255;
    imd0 = (mdl & mask) << 24;
    imd0 = imd0 | ((nc & mask) << 16);
    imd0 = imd0 | ((abpi & mask) << 8);
    imd0 = imd0 | (sbpi & mask);

    imeta_data[0] = imd0;
    /* bpi,bptsv,bpi,color */
    color = z_dom;
    imd1 = tbpi & mask;
    imd1 = (imd1 << 8) | (bptsv & mask);
    imd1 = (imd1 << 8) | (bpe & mask);
    imd1 = (imd1 << 8) | (color & mask);
    imeta_data[1] = imd1;

    imd2 = nx & mask;
    imd2 = (imd2 << 8) | (ny & mask);
    imd2 = (imd2 << 8) | (nz & mask);
    sum_nxyz = nx + ny + nz;
    imd2 = (imd2 << 8) | (sum_nxyz & mask);
    imeta_data[2] = imd2;

    imeta_data[3] = max_index;
    imeta_data[4] = ner;
    /*
      We don't set imeta_data[5] for the z_dominant combinations,
      as we don't bin them.
    */
    /* imd6 and imd7  get zfd_size */
    jmd_p = (int64_t*)&imeta_data[6];
    *jmd_p = zfd_size;
    /* 
       imd8 and imd9  will store the number of combinations <= max_elem 
    */
    /* 
       imd10 and imd11 will store the position in 8 byte words of the
       position of the first index element, only used for y_dom and x_dom 
       channels.
    */
    dmd_p = (double*)&imeta_data[12];
    *dmd_p = min_fract;
    /*
      We also don't need to store the combo_bin_width in imeta_data[14:15]
    */
    cimap = (int*)&imeta_data[16];
    if (ner > 0) {
      cimap = (int*)&imeta_data[3*bpe/4];
    }
    pack_imap(nz,&imap_sort[0],cimap);
    ddata = (double*)&imeta_data[mdl/4];

    *ddata_p  = ddata;
    *imeta_data_p = imeta_data;
  }
  return (success);
}
