#include "system_includes.h"
#include "state_struct.h"
#include "size_file.h"
#include "unpack_isteps.h"
#include "unpack_imap.h"
int main (int argc, char **argv) {
  int64_t *j_ptr;
  int64_t *ibin_data;
  double  *d_ptr;
  double  *d_data;
  char *coord_file;
  char *output_file;
  char *coord_words[3];
  char coord_words_c[24];
  char coord_file_c[128];
  char output_file_c[128];
  void *null_addr;
  void *m_addr;
  int *imeta_data;
  int isteps_index[64];
  int imap[64];
  int *isteps_packed;
  int *imap_packed;
  double min_fract;
  double combo_bin_width;
  double vb_min;
  double vb_max;
  double ub_min;
  double ub_max;
  int64_t file_size;
  int64_t offset;
  int64_t num_combos;
  int64_t n_combo_words;
  int64_t coord_size;
  int64_t j0;
  int64_t ibin_size;
  int cfd;
  int prot;

  int flags;
  int mask;

  int sbpi;
  int abpi;

  int nc;
  int mdl;

  int imde;
  int tbpi;

  int bptsv;
  int bpe;

  int coord;
  int nx;

  int ny;
  int nz;

  int sum_nxyz;
  int max_index;

  int ner;
  int num_bins;

  int nfields;
  int tsv_offset;

  int id_offset;
  int i;

  int j;
  int save_err;

  FILE *cfp;
  FILE *ofp;

  int success;
  success = 1;
  coord_words[0] = (char*)&coord_words_c[0];
  coord_words[1] = (char*)&coord_words_c[8];
  coord_words[2] = (char*)&coord_words_c[16];
  strcpy(coord_words[0],"x_dom");
  strcpy(coord_words[1],"y_dom");
  strcpy(coord_words[2],"z_dom");
  coord_file = &coord_file_c[0];
  output_file = &output_file_c[0];
  if (argc < 2) {
    fprintf(stderr,"view_coord_file error: no coord file specified\n");
    success = 0;
  } 
  if (success) {
    j0 = (int64_t)0;
    strcpy(coord_file,argv[1]);
    file_size = size_file(coord_file,stderr);
    if (file_size == j0) {
      success = 0;
      fprintf(stderr,"view_coord_file error: coord file %s is empty\n",
	      coord_file);
      fflush(stderr);
    }
  }
  if (success) {
    cfp = fopen(coord_file,"r");
    if (cfp == NULL) {
      success = 0;
      fprintf(stderr,"view_coord_file: Unable to open %s\n",coord_file);
      fflush(stderr);
    }
  }
  if (success) {
    strcpy(output_file,coord_file);
    strcat(output_file,".ascii");
    ofp = fopen(output_file,"w");
    if (ofp == NULL) {
      success = 0;
      fprintf(stderr,"view_coord_file: Unable to open %s\n",output_file);
      fflush(stderr);
    }
  }
  if (success) {
    /*
      mmap coord file.
    */
    cfd       = fileno(cfp);
    null_addr = NULL;
    offset    = (int64_t)0;
    prot      = PROT_READ;
    flags     = MAP_PRIVATE;
    m_addr    = mmap(null_addr,file_size,prot,flags,cfd,offset);
    save_err  = errno;
    if (m_addr == MAP_FAILED) {
      success = 0;
      fprintf(stderr,
	      "view_coord_file: Error, unable to mmap %s, error = %d, %s \n",
	      coord_file,save_err,strerror(save_err));
      fflush(stderr);
    }
  }
  if (success) {
    mask = 255;
    imeta_data = (int*)m_addr;
    imde       = imeta_data[0];
    sbpi = imde & mask;
    imde = imde >> 8;
    abpi = imde & mask;
    imde = imde >> 8;
    nc   = imde & mask;
    imde = imde >> 8;
    mdl  = imde & mask;

    fprintf(ofp,"Meta_data length, mdl              = %d\n",mdl);
    fprintf(ofp,"Number of channels, nc             = %d\n",nc);
    fprintf(ofp,"Active bytes per istep index, abip = %d\n",abpi);
    fprintf(ofp,"Stored bytes per istep index, sbpi = %d\n",sbpi);
    
    imde  = imeta_data[1];
    coord = imde & mask;
    imde  = imde >> 8;
    bpe   = imde & mask;
    imde  = imde >> 8;
    bptsv = imde & mask;
    imde  = imde >> 8;
    tbpi  = imde & mask;
    fprintf(ofp,"Total bytes per index set,tbpi     = %d\n",tbpi);
    fprintf(ofp,"bytes per tristimulus value,bptsv  = %d\n",bptsv);
    fprintf(ofp,"bytes per combination element      = %d\n",bpe);
    fprintf(ofp,"dominant_coord                     = %d, %s\n",
	    coord,coord_words[coord]);

    imde  = imeta_data[2];
    sum_nxyz = imde & mask;
    imde     = imde >> 8;
    nz       = imde & mask;
    imde     = imde >> 8;
    ny       = imde & mask;
    imde     = imde >> 8;
    nx       = imde & mask;
    fprintf(ofp,"nx (number x_dom channels)         = %d\n",nx);
    fprintf(ofp,"ny (number y_dom channels)         = %d\n",ny);
    fprintf(ofp,"nz (number z_dom channels)         = %d\n",nz);
    fprintf(ofp,"nx+ny+nz (total number channels)   = %d\n",sum_nxyz);

    max_index = imeta_data[3];
    fprintf(ofp,"max_index                          = %d\n",max_index);
    ner       = imeta_data[4];
    fprintf(ofp,"Number of extra records for imap   = %d\n",ner);
    num_bins  = imeta_data[5];
    fprintf(ofp,"Number of bins                     = %d\n",num_bins);

    j_ptr     = (int64_t*)&imeta_data[6];
    coord_size = *j_ptr;
    fprintf(ofp,"stored coord file size             = %ld\n",coord_size);
    fprintf(ofp,"measured coord file size           = %ld\n",file_size);
    
    j_ptr     = (int64_t*)&imeta_data[8];
    num_combos = *j_ptr;
    fprintf(ofp,"Number of combinations             = %ld\n",num_combos);

    j_ptr     = (int64_t*)&imeta_data[10];
    n_combo_words = *j_ptr;
    fprintf(ofp,"Number of combo words              = %ld\n",n_combo_words);
    
    d_ptr     = (double*)&imeta_data[12];
    min_fract = *d_ptr;
    fprintf(ofp,"min_fract, (smmalest power_step    = %le\n",min_fract);
    d_ptr     = (double*)&imeta_data[14];
    combo_bin_width = *d_ptr;
    fprintf(ofp,"combo_bin_width                    = %le\n",combo_bin_width);

    d_ptr     = (double*)&imeta_data[16];
    vb_min    = *d_ptr;
    fprintf(ofp,"vb_min                             = %le\n",vb_min);
    d_ptr     = (double*)&imeta_data[18];
    vb_max    = *d_ptr;
    fprintf(ofp,"vb_max                             = %le\n",vb_max);
    d_ptr     = (double*)&imeta_data[20];
    ub_min    = *d_ptr;
    fprintf(ofp,"ub_min                             = %le\n",ub_min);
    d_ptr     = (double*)&imeta_data[22];
    ub_max    = *d_ptr;
    fprintf(ofp,"ub_max                             = %le\n",ub_max);

    tsv_offset = tbpi >> 3;
    if (coord != 2) {
      nfields    = tsv_offset + 5;
    } else {
      nfields    = tsv_offset + 3;
    }

    fprintf(ofp,"tsv_offset                        = %d\n",tsv_offset);
    fprintf(ofp,"nfields                           = %d\n",nfields);

    d_ptr = (double *)m_addr;

    if (ner > 0) {
      imap_packed = (int*)&d_ptr[3*nfields];
    } else {
      imap_packed = &imeta_data[24];
    }
    unpack_imap(nc,imap,imap_packed);
    fprintf(ofp,"imap:\n");
    for (i=0;i<nc;i++) {
      fprintf(ofp,"%d\t",imap[i]);
    }
    fprintf(ofp,"\n");

    id_offset = (ner + 3)*nfields;
    d_data = &d_ptr[id_offset];
    
    if (coord != 2) {
      fprintf(ofp,"\n\t\tsteps\t\t\ttsv_x\ttsv_y\ttsv_z\tvb\tub\n");
    } else {
      fprintf(ofp,"\n\t\tsteps\t\t\ttsv_x\ttsv_y\ttsv_z\n");
    }
    fflush(ofp);
    for (i=0;i<num_combos;i++) {
      isteps_packed = (int*)d_data;
      unpack_isteps(nc,sbpi,isteps_index,isteps_packed);
      for (j=0;j<nc;j++) {
	fprintf(ofp,"%d\t",isteps_index[j]);
      }
      if (coord != 2) {
	fprintf(ofp,"\t%le\t%le\t%le\t%le\t%le\n",d_data[tsv_offset],
		d_data[tsv_offset+1],d_data[tsv_offset+2],
		d_data[tsv_offset+3],d_data[tsv_offset+4]);
      } else {
	fprintf(ofp,"\t%le\t%le\t%le\n",d_data[tsv_offset],
		d_data[tsv_offset+1],d_data[tsv_offset+2]);
      }
      d_data += nfields; /* Caution address arithmetic */
    }
    if (coord != 2) {
      fprintf(ofp,"\nbin_number\tcombo_pos\n");
      if (num_bins > 0) {
	ibin_data = (int64_t*)d_data;
	for (i=1;i<num_bins;i++) {
	  ibin_size = ibin_data[i]-ibin_data[i-1];
	  fprintf(ofp,"%d\t%ld\t%ld\n",i,ibin_data[i],ibin_size);
	}
      }
    }
    fclose(ofp);
    munmap(m_addr,file_size);
    fclose(cfp);
  }
  return(0);
}

	      
      
