#include "system_includes.h"
#include "state_struct.h"
#include "count_ws.h"
#include "count_nws.h"
#include "upcase.h"
#include "record_steps_per_channel.h"
#include "record_scales.h"
#include "record_subset.h"
#include "record_vrb.h"
#include "record_dbg_steps.h"
#include "set_test_radii_sq.h"
#include "read_params.h"
int read_params(struct state_struct *state) {
  /*
    Set the defaults for scalar fields of the era_opt state structure.
    Open and read a parameter file if it exists overwriting the defaults.
    Called by: era_opt,gen_combos
    Calls:     fopen, strncmp, strlen, strcpy, fopen, fgets, sscanf,
               fprintf, fflush
  */
  double *scales;
  double *vrb_radii;
  double *test_radii_sq;
  double *dbg_steps;
  char *parameter_file;      
  char *log_file;            
  char *source_file;         
  char *channels_file;       
  char *z_dom_file;
  char *y_dom_file;
  char *x_dom_file;
  char *input_dir;
  char *combo_dir;
  char *output_dir;
  /*
  char *quadrangles_file;    
  */
  char *huebins_file;      
  char *output_file;         
  char *summary_file;
  char *cct_file;
  char *s026_file;
  char *subset_file;
  char *combos_in_file;
  char *combos_out_file;
  char *io_buffer;
  char *key;
  char *value;
  char *fgp;
  char *dir_sep;
  char *dir_end;
  char *mcat02inv_file;
  int  *isteps_per_channel;
  int  *vrb_first_index;
  int  *vrb_last_index;
  int  *subset_indices;

  int64_t buffer_length;

  double u_target;
  double v_target;
  double test_rad;
  double rad_skosh;
  double rf_thresh;
  double neutral_rcs_radius;  /* if all rcsh huebin lengths are < neutral_rcs_radius, the LCS group is neutral */
  double neutral_skew_radius; /* if major axis length/minor axis_length < 1 + neutral_skew_radius ellipse 
				 group is neutral */

  double watchme;
  int source_rows;
  int source_columns;

  int channels_rows;
  int channels_columns;

  int first_channel_index;
  int last_channel_index;

  int num_radiator_temps;
  int num_samples;

  /*
  int num_quadrangles;
  */
  int num_huebins;
  int wavelength_index;

  int xbar_index;
  int ybar_index;

  int zbar_index;
  int xbar_10_index;

  int ybar_10_index;
  int zbar_10_index;

  int s_0_index;
  int s_1_index;

  int s_2_index;
  int ns;

  int nskip;
  int key_len;

  int line_no;
  int success;

  int num_channels;
  int nc_steps_set;

  int found;
  int bad_scan;

  int num_cct_rows;
  int max_scales;

  int number_scales;
  int use_multiscale;

  int extra_log_info;
  int use_vrb;

  int number_vrb;
  int max_radii_bins;

  int sz_index;
  int s026_columns;

  int max_num_bins;
  int nz_cpcg_max;

  int subset_size;
  int max_subset_size;

  int full_num_channels;
  int use_bruteforce;

  int full_output;
  int select_only;

  int discount_illuminant;
  int new_cam02ucs;

  int dir_len;
  int len_dir_sep;

  int num_dbg_steps;
  int hp_out;

  int save_steps;
  int read_combos;

  int use_alt_cct;
  int read_mcat02inv;

  int cam02ucs_1mvp;
  int cam02ucs_2mvp;
    
  int sref_mode;    /* 0 for compute sref from t_t (original)
		       1 for lookup sref
		       2 for interpolate sref
		    */
  int t_t_inc;      /* increment value to use in building tables
		       for sref_mode = 1 or 2 */
  int sref_t_t_min; /* minimum t_t value for sref_mode =1 or 2*/
  int sref_t_t_max; /* maximum t_t value for sref_mode =1 or 2*/

  int dump_combos;
  int ipad;

  size_t ival_len;
  size_t longest_dir_len;
  size_t filename_length;

  FILE *pfp;
  FILE *lfp;
  parameter_file   = state->parameter_file;
  /*
    Set default values.
  */
  log_file         = state->log_file;
  source_file      = state->source_file;
  channels_file    = state->channels_file;
  output_file      = state->output_file;
  summary_file     = state->summary_file;
  combos_in_file   = state->combos_in_file;
  combos_out_file  = state->combos_out_file;
  cct_file         = state->cct_file;
  s026_file        = state->s026_file;
  subset_file      = state->subset_file;
  combos_in_file   = state->combos_in_file;
  combos_out_file  = state->combos_out_file;
  z_dom_file       = state->z_dom_file;
  y_dom_file       = state->y_dom_file;
  x_dom_file       = state->x_dom_file;
  input_dir        = state->input_dir;
  combo_dir        = state->combo_dir;
  output_dir       = state->output_dir;
  mcat02inv_file   = state->mcat02inv_file;
  isteps_per_channel = state->isteps_per_channel;
  scales           = state->scales;
  dbg_steps        = state->dbg_steps;
  max_scales       = state->max_scales;
  max_subset_size  = state->max_subset_size;
  subset_indices   = state->subset_indices;
  max_radii_bins   = state->max_radii_bins;
  vrb_radii        = state->vrb_radii;
  test_radii_sq    = state->test_radii_sq;
  vrb_first_index  = state->vrb_first_index;
  vrb_last_index   = state->vrb_last_index;
  /*
  quadrangles_file = state->quadrangles_file;
  */
  huebins_file     = state->huebins_file;
  io_buffer        = state->io_buffer;
  dir_sep = "/"; /* string that terminates a directory name. Use \ for broken Windows */
  len_dir_sep = 1;
  filename_length = (size_t)state->filename_length - 1;
  longest_dir_len = filename_length - len_dir_sep;
  buffer_length    = (int64_t)state->buffer_length;

  success          = 1;
  strcpy(log_file,"era_opt.log");
  strcpy(source_file,"source_tm30_new_20210623.csv");
  strcpy(channels_file,"test_spds.csv");
  strcpy(huebins_file,"huebins.csv");
  strcpy(cct_file,"cct_20210708.csv");
  strcpy(output_file,"era_opt.out");
  strcpy(summary_file,"era_opt.sum");
  strcpy(s026_file,"cie_s026.csv");
  strcpy(input_dir,""); /* default is empty == current working directory */
  strcpy(combo_dir,""); /* default is empty == current working directory */
  strcpy(output_dir,""); /* default is empty == current working directory */
  /*
  strcpy(quadrangles_file,"checkansi.csv");
  */
  source_rows = 401;
  source_columns = 109;
  channels_rows = 401;
  channels_columns = 19;
  num_channels = 4;
  full_num_channels = 18;
  max_subset_size = 64;
  num_radiator_temps = 951;
  /*
  num_cct_rows        = 1201;
  changing to 1217 for the new cct_20210708.csv file
  */
  num_cct_rows        = 1217;
  first_channel_index =  2;
  last_channel_index  =  5;
  sz_index            =  4;
  s026_columns        =  8;
  num_samples         =  99;
  wavelength_index    =  100;
  xbar_index          =  101;
  ybar_index          =  102;
  zbar_index          =  103;
  xbar_10_index       =  104;
  ybar_10_index       =  105;
  zbar_10_index       =  106;
  s_0_index           =  107;
  s_1_index           =  108;
  s_2_index           =  109;
  num_huebins         =  16;
  u_target            = 0.2357;
  v_target            = 0.5112;
  test_rad            = 0.0005;
  rad_skosh           = 0.0;
  rf_thresh           = 0.0;
  isteps_per_channel[0] = 6;
  nc_steps_set          = 1;
  use_multiscale        = 0;
  use_bruteforce        = 0;
  select_only           = 0;
  number_scales         = 1;
  subset_size           = 0;
  extra_log_info        = 0;
  full_output           = 0;
  use_vrb               = 0;
  discount_illuminant   = 1;
  new_cam02ucs          = 0;
  use_alt_cct           = 0;
  number_vrb            = 4;
  state->number_vrb     = 4;
  max_num_bins          = 4096;
  sref_mode             = 0; /* Calclate sref from t_t using create_sref */
                             /* = 1 to use look up table determined in prep */
                             /* = 2 to interpolate from table in prep. */
  t_t_inc               = 1;
  sref_t_t_min          = 3480;
  sref_t_t_max          = 3520;
  dump_combos           = 0;
  /*
    By default multiply by mhpe and mcat02^-1 separately
    in cam02ucs
  */
  cam02ucs_1mvp         = 0;
  cam02ucs_2mvp         = 1;
  nz_cpcg_max           = -1;
  vrb_first_index[0]    = 1;
  vrb_last_index[0]     = 5;
  vrb_radii[0]          = 5.0e-4;
  vrb_first_index[1]    = 6;
  vrb_last_index[1]     = 9;
  vrb_radii[1]          = 1.0e-4;
  vrb_first_index[2]    = 10;
  vrb_last_index[2]     = 13;
  vrb_radii[2]          = 2.0e-5;
  vrb_first_index[3]    = 14;
  vrb_last_index[3]     = 18;
  vrb_radii[3]          = 4.0e-6;
  /*
    The r_cs_h values are percentages so change this to a percentage
    as per discussion with Michael Royer, 8/20/21
  neutral_rcs_radius    = 0.03;   
  */
  neutral_rcs_radius    = 3.0;   
  neutral_skew_radius   = 0.04;   
  hp_out                = 0; /* 1 for high precsion output, f16.12 instead of f6.2 */
  save_steps            = 0; /* 1 for save step indices and distances of combos in tca */
  read_combos           = 0; /* 1 =read combos from file instead of using era_opt to 
				generate them */

  read_mcat02inv        = 1;   /* 0 to compute the inverse of mcat02, 1 to read it in */
  watchme               = 5.0; /* Hardwired for now will allow change later 
				  Print progress message every 5% of the way */
  /*
    Initialize the coordinate dominant files to be null.
  */
  z_dom_file[0]         = '\0';
  y_dom_file[0]         = '\0';
  x_dom_file[0]         = '\0';
  strcpy(subset_file,"subset.csv");
  strcpy(combos_in_file,"combos.csv");
  strcpy(combos_out_file,"combos_out.csv");
  strcpy(mcat02inv_file,"mcat02inv.csv");
  /*
    Try to open parameter file.
  */
  if (strlen(parameter_file) > 0) {
    pfp = fopen(parameter_file,"r");
    if (pfp != NULL) {
      /*
	Read keyword value pairs from parameter file.
      */
      fgp = fgets(io_buffer,buffer_length,pfp);
      line_no = 1;
      while (fgp) {
	key = io_buffer;
	nskip = count_ws(key);
	key += nskip; /* Caution address arithmetic */
	if (key[0] != '#') {
	  key_len = count_nws(key);
	  upcase(key_len,key);
	  
	  value = key + key_len + 1; /* Caution address arithmetic */
	  /*
	    Remove leading whitespace from value.
	  */
	  nskip = count_ws(value);
	  value += nskip; /* Caution address arithmetic */
	  ival_len = strlen(value);
	  found = 1;
	  /*
	    fprintf(stdout,"read_params: key   = %s\n"
	    "                    value = %s\n",key,value);
	    fflush(stdout);
	  */
	  bad_scan = 0;
	  if (key_len > 13) {
	    if (key_len > 16) {
	      if (key_len == 19) {
		if (strncmp(key,"FIRST_CHANNEL_INDEX",19) == 0) {
		  ns = sscanf(value,"%d",&first_channel_index);
		  bad_scan = (ns == 0);
		} else {
		  if (strncmp(key,"NUM_CHANNEL_COLUMNS",19) == 0) {
		    ns = sscanf(value,"%d",&channels_columns);
		    bad_scan = (ns == 0);
		  } else {
		    if (strncmp(key,"NEUTRAL_SKEW_RADIUS",19) == 0) {
		      ns = sscanf(value,"%le",&neutral_skew_radius);
		      bad_scan = (ns == 0);
		    } else {
		      if (strncmp(key,"DISCOUNT_ILLUMINANT",19) == 0) {
			ns = sscanf(value,"%d",&discount_illuminant);
			bad_scan = (ns == 0);
		      } else {
			found = 0;
		      }
		    }
		  }
		}
	      } else {
		if (key_len == 18) {
		  if(strncmp(key,"LAST_CHANNEL_INDEX",18) == 0) {
		    ns = sscanf(value,"%d",&last_channel_index);
		    bad_scan = (ns == 0);
		  } else {
		    if(strncmp(key,"NUM_SOURCE_COLUMNS",18) == 0) {
		      ns = sscanf(value,"%d",&source_columns);
		      bad_scan = (ns == 0);
		    } else {
		      if (strncmp(key,"NEUTRAL_RCS_RADIUS",18) == 0) {
			ns = sscanf(value,"%le",&neutral_rcs_radius);
			bad_scan = (ns == 0);
		      } else {
			found = 0;
		      }
		    }
		  }
		} else { 
		  if (strncmp(key,"STEPS_PER_CHANNEL",17) == 0) {
		    record_steps_per_channel(value,&nc_steps_set,
					     isteps_per_channel,stderr);
		    /*
		      fprintf(stderr,"isteps_per_channel[0] = %d\n",
		      isteps_per_channel[0]);
		    */
		    state->nc_steps_set=nc_steps_set;
		  } else { /* key_len was 17 or >19 - no match */
		    found = 0;
		  }		
		}
	      }
	    } else { 
	      /* key_len <= 16 */
	      if (key_len == 16) {
		if (strncmp(key,"NUM_CHANNEL_ROWS",16) == 0) {
		  ns = sscanf(value,"%d",&channels_rows);
		  bad_scan = (ns == 0);
		} else {
		  if (strncmp(key,"WAVELENGTH_INDEX",16) == 0) {
		    ns = sscanf(value,"%d",&wavelength_index);
		    bad_scan = (ns == 0);
		  } else {
		    found = 0;
		  }
		}
	      } else {
		if (key_len == 15) {
		  /*
		    if (strncmp(key,"NUM_QUADRANGLES",15) == 0) {
		    ns = sscanf(value,"%d",&num_quadrangles);
		    bad_scan = (ns == 0);
		    } else {
		    }
		  */               
		  if (strncmp(key,"NUM_SOURCE_ROWS",15) == 0) {
		    ns = sscanf(value,"%d",&source_rows);
		    channels_rows = source_rows;
		    bad_scan = (ns == 0);
		  } else {
		    if (strncmp(key,"COMBOS_OUT_FILE",15) == 0) {
		      if (ival_len <= filename_length) {
			ns = sscanf(value,"%s",combos_out_file);
			bad_scan = (ns == 0);
		      } else {
			strncpy(combos_out_file,value,filename_length);
			combos_in_file[filename_length] = '\0';
			bad_scan = 0;
			fprintf(stdout,"read_params: Warning "
				"combos_out_file truncated to %s\n",
				combos_out_file);
			fflush(stdout);
		      }
		    } else {
		      found = 0;
		    }
		  }
		} else {
		  if (key_len == 14) {
		    if (strncmp(key,"EXTRA_LOG_INFO",14) == 0) {
		      ns = sscanf(value,"%d",&extra_log_info);
		      bad_scan = (ns == 0);
		    } else {
		      if (strncmp(key,"USE_MULTISCALE",14) == 0) {
			ns = sscanf(value,"%d",&use_multiscale);
			bad_scan = (ns == 0);
		      } else {
			if (strncmp(key,"USE_BRUTEFORCE",14) == 0) {
			  ns = sscanf(value,"%d",&use_bruteforce);
			  bad_scan = (ns == 0);
			} else {
			  if (strncmp(key,"READ_MCAT02INV",14) == 0) {
			    ns = sscanf(value,"%d",&read_mcat02inv);
			    bad_scan = (ns == 0);
			  } else {
			    if (strncmp(key,"MCAT02INV_FILE",14) == 0) {
			      if (ival_len <= filename_length) {
				ns = sscanf(value,"%s",mcat02inv_file);
				bad_scan = (ns == 0);
			      } else {
				/* Should probably print a warning to stdout here */
				strncpy(mcat02inv_file,value,filename_length);
				mcat02inv_file[filename_length] = '\0';
				fprintf(stdout,"WARNING: mcat02inv_file truncated to: %s\n"<
					mcat02inv_file);
				fflush(stdout);
				bad_scan = (ns == 0);
			      }
			    } else {
			      if (strncmp(key,"COMBOS_IN_FILE",14) == 0) {
				if (ival_len <= filename_length) {
				  ns = sscanf(value,"%s",combos_in_file);
				  bad_scan = (ns == 0);
				} else {
				  strncpy(combos_in_file,value,filename_length);
				  combos_in_file[filename_length] = '\0';
				  bad_scan = 0;
				  fprintf(stdout,"read_params: Warning "
					  "combos_in_file truncated to %s\n",
					  combos_in_file);
				  fflush(stdout);
				}
			      } else {
				found = 0;
			      }
			    }
			  }
			}
		      }
		    }
		  } else {
		    found = 0;
		  }
		}
	      } /* end else key_len != 16 */
	    } /* end else key_len <= 16 */
	  } else {
	    /* 
	       key_len <= 13
	    */
	    if (key_len > 11) {
	      if (key_len == 13) {
		if (strncmp(key,"CHANNELS_FILE",13) == 0) {
		  if (ival_len <= filename_length) {
		    ns = sscanf(value,"%s",channels_file);
		    bad_scan = (ns == 0);
		  } else {
		    /* Should probably print a warning to stdout here */
		    strncpy(channels_file,value,filename_length);
		    channels_file[filename_length] = '\0';
		    fprintf(stdout,"WARNING: channels_file truncated to: %s\n"<
			    channels_file);
		    fflush(stdout);
		    bad_scan = 0;
		  }
		} else {
		  if (strncmp(key,"XBAR_10_INDEX",13) == 0) {
		    ns = sscanf(value,"%d",&xbar_10_index);
		    bad_scan = (ns == 0);
		  } else {
		    if (strncmp(key,"YBAR_10_INDEX",13) == 0) {
		      ns = sscanf(value,"%d",&ybar_10_index);
		      bad_scan = (ns == 0);
		    } else {
		      if (strncmp(key,"ZBAR_10_INDEX",13) == 0) {
			ns = sscanf(value,"%d",&zbar_10_index);
			bad_scan = (ns == 0);
		      } else {
			if (strncmp(key,"CAM02UCS_1MVP",13) == 0) {
			  ns = sscanf(value,"%d",&cam02ucs_1mvp);
   		          bad_scan = (ns == 0);
			  if (cam02ucs_1mvp != 0) {
			    cam02ucs_2mvp = 0;
			    cam02ucs_1mvp = 1;
			  } else {
			    cam02ucs_2mvp = 1;
			  }
			} else {
			  if (strncmp(key,"CAM02UCS_2MVP",13) == 0) {
			    ns = sscanf(value,"%d",&cam02ucs_2mvp);
			    bad_scan = (ns == 0);
			    if (cam02ucs_2mvp != 0) {
			      cam02ucs_1mvp = 0;
			      cam02ucs_2mvp = 1;
			    } else {
			      cam02ucs_1mvp = 1;
  	   		    }  
			  } else {
			    found = 0;
			  }
			}
		      }
		    }
		  }
		}
	      } else {
		/* 
		   key_len must be 12 
		*/
		if (strncmp(key,"HUEBINS_FILE",12) == 0) {
		  if (ival_len <= filename_length) {
		    ns = sscanf(value,"%s",huebins_file);
		    bad_scan = (ns == 0);
		  } else {
		    /* Should probably print a warning to stdout here */
		    strncpy(huebins_file,value,filename_length);
		    huebins_file[filename_length] = '\0';
		    fprintf(stdout,"WARNING: huebins_file truncated to: %s\n"<
			    huebins_file);
		    fflush(stdout);
		    bad_scan = 0;
		  }
		} else {
		  if (strncmp(key,"NUM_CHANNELS",12) == 0) {
		    ns = sscanf(value,"%d",&num_channels);
		    channels_columns = num_channels + 1;
		    first_channel_index = 2;
		    last_channel_index = channels_columns;
		    bad_scan = (ns == 0);
		  } else {
		    if (strncmp(key,"NUM_CCT_ROWS",12) == 0) {
		      ns = sscanf(value,"%d",&num_cct_rows);
		      bad_scan = (ns == 0);
		    } else { 
		      if (strncmp(key,"S026_COLUMNS",12) == 0) {
			ns = sscanf(value,"%d",&s026_columns);
			bad_scan = (ns == 0);
		      } else {
			if (strncmp(key,"MAX_NUM_BINS",12) == 0) {
			  ns = sscanf(value,"%d",&max_num_bins);
			  bad_scan = (ns == 0);
			} else {
			  if (strncmp(key,"NEW_CAM02UCS",12) == 0) {
			    ns = sscanf(value,"%d",&new_cam02ucs);
			    bad_scan = (ns == 0);
			  } else {
			    if (strncmp(key,"SUMMARY_FILE",12) == 0) {
			      ns = sscanf(value,"%s",summary_file);
			      bad_scan = (ns == 0);
			    } else {
			      if (strncmp(key,"SREF_T_T_MIN",12) == 0) {
				ns = sscanf(value,"%d",&sref_t_t_min);
				bad_scan = (ns == 0);
			      } else {
				if (strncmp(key,"SREF_T_T_MAX",12) == 0) {
				  ns = sscanf(value,"%d",&sref_t_t_max);
				  bad_scan = (ns == 0);
				} else {
				  found = 0;
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    } else {
	      /*
		key_len <= 11
	      */
	      if (key_len == 11) {
		if (strncmp(key,"NUM_HUEBINS",11) == 0) {
		  ns = sscanf(value,"%d",&num_huebins);
		  bad_scan = (ns == 0);
		} else {
		  if (strncmp(key,"NUM_SAMPLES",11) == 0) {
		    ns = sscanf(value,"%d",&num_samples);
		    bad_scan = (ns == 0);
		  } else {
		    if (strncmp(key,"SOURCE_FILE",11) == 0) {
		      if (ival_len <= filename_length) {
			ns = sscanf(value,"%s",source_file);
			bad_scan = (ns == 0);
		      } else {
			strncpy(source_file,value,filename_length);
			source_file[filename_length] = '\0';
			fprintf(stdout,"read_params: WARNING source file truncated to %s\n",
				source_file);
			fflush(stdout);
			bad_scan = 0;
		      }
		    } else {
		      if (strncmp(key,"OUTPUT_FILE",11) == 0) {
			if (ival_len <= filename_length) {
			  ns = sscanf(value,"%s",output_file);
			  bad_scan = (ns == 0);
			} else {
			  strncpy(output_file,value,filename_length);
			  output_file[filename_length] = '\0';
			  bad_scan = 0;
			  fprintf(stdout,"read_params: Warning output file name truncated to %s\n",
				  output_file);
			  fflush(stdout);
			}
		      } else {
			if (strncmp(key,"NZ_CPCG_MAX",11) == 0) {
			  ns = sscanf(value,"%d",&nz_cpcg_max);
			  bad_scan = (ns == 0);
			} else {
			  if (strncmp(key,"SUBSET_FILE",11) == 0) {
			    if (ival_len <= filename_length) {
			      ns = sscanf(value,"%s",subset_file);
			      bad_scan = (ns == 0);
			    } else {
			      strncpy(subset_file,value,filename_length);
			      subset_file[filename_length] = '\0';
			      bad_scan = 0;
			      fprintf(stdout,"read_params: Warning subset_file truncated to %s\n",
				      subset_file);
			      fflush(stdout);
			    }
			  } else {
			    if (strncmp(key,"FULL_OUTPUT",11) == 0) {
			      ns = sscanf(value,"%d",&full_output);
			      bad_scan = (ns == 0);
			    } else {
			      if (strncmp(key,"SELECT_ONLY",11) == 0) {
				ns = sscanf(value,"%d",&select_only);
				bad_scan = (ns == 0);
			      } else {
				if (strncmp(key,"READ_COMBOS",11) == 0) {
				  ns = sscanf(value,"%d",&read_combos);
				  bad_scan = (ns == 0);
				} else {
				  if (strncmp(key,"USE_ALT_CCT",11) == 0) {
				    ns = sscanf(value,"%d",&use_alt_cct);
				    bad_scan = (ns == 0);
				  } else {
				    if (strncmp(key,"DUMP_COMBOS",11) == 0) {
				      ns = sscanf(value,"%d",&dump_combos);
				      bad_scan = (ns == 0);
				    } else {
				      found = 0;
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      } else {
		if (key_len == 10) {
		  if (strncmp(key,"XBAR_INDEX",10) == 0) {
		    ns = sscanf(value,"%d",&xbar_index);
		    bad_scan = (ns == 0);
		  } else {
		    if (strncmp(key,"YBAR_INDEX",10) == 0) {
		      ns = sscanf(value,"%d",&ybar_index);
		      bad_scan = (ns == 0);
		    } else {
		      if (strncmp(key,"ZBAR_INDEX",10) == 0) {
			ns = sscanf(value,"%d",&zbar_index);
			bad_scan = (ns == 0);
		      } else {
			if (strncmp(key,"X_DOM_FILE",10) == 0) {
			  if (ival_len <= filename_length) {
			    ns = sscanf(value,"%s",x_dom_file);
			    bad_scan = (ns == 0);
			  } else {
			    strncpy(x_dom_file,value,filename_length);
			    x_dom_file[filename_length] = '\0';
			    bad_scan = 0;
			    fprintf(stdout,"read_params: WARNING x_dom_file truncated to %s\n",
				    x_dom_file);
			    fflush(stdout);
			  }
			} else {
			  if (strncmp(key,"Y_DOM_FILE",10) == 0) {
			    if (ival_len <= filename_length) {
			      ns = sscanf(value,"%s",y_dom_file);
			      bad_scan = (ns == 0);
			    } else {
			      strncpy(y_dom_file,value,filename_length);
			      y_dom_file[filename_length] = '\0';
			      bad_scan = 0;
			      fprintf(stdout,"read_params: WARNING y_dom_file truncated to %s\n",
				      y_dom_file);
			      fflush(stdout);
			    }
			  } else {
			    if (strncmp(key,"Z_DOM_FILE",10) == 0) {
			      if (ival_len < filename_length) {
				ns = sscanf(value,"%s",z_dom_file);
				bad_scan = (ns == 0);
			      } else {
				strncpy(z_dom_file,value,filename_length);
				z_dom_file[filename_length] = '\0';
				bad_scan = 0;
				fprintf(stdout,"read_params: WARNING z_dom_file truncated to %s\n",
					z_dom_file);
			      }
			    } else {
			      if (strncmp(key,"OUTPUT_DIR",10) == 0) {
				if (ival_len <= longest_dir_len) {
				  ns = sscanf(value,"%s",output_dir);
				  bad_scan = (ns == 0) ;
				} else {
				  strncpy(output_dir,value,longest_dir_len);
				  output_dir[longest_dir_len] = '\0';
				  bad_scan = 0;
				  ns = 1;
				  fprintf(stdout,"read_params: WARNING output_dir truncated to %s\n",
					  output_dir);
				  fflush(stdout);
				}
				/*
				  Ensure that a nonempty directory ends in dir_sep ("/")
				*/
				if (ns != 0) {
				  dir_len = strlen(output_dir);
				  if (dir_len > 0) {
				    dir_end = &(output_dir[dir_len-len_dir_sep]);
				    if (strcmp(dir_sep,dir_end) != 0) {
				      if (dir_len < longest_dir_len) {
					strcat(output_dir,dir_sep);
				      } else {
					dir_end = &output_dir[longest_dir_len];
					strcat(dir_end,dir_sep);
				      }
				    }
				  }
				}
			      } else {
				if (strncmp(key,"SAVE_STEPS",10) == 0) {
				  ns = sscanf(value,"%d",&save_steps);
				  bad_scan = (ns == 0);
				} else {
				  found = 0;
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		} else {
		  if (key_len == 9) {
		    if (strncmp(key,"S_0_INDEX",9) == 0) {
		      ns = sscanf(value,"%d",&s_0_index);
		      bad_scan = (ns == 0);
		    } else {
		      if (strncmp(key,"S_1_INDEX",9) == 0) {
			ns = sscanf(value,"%d",&s_1_index);
			bad_scan = (ns == 0);
		      } else {
			if (strncmp(key,"S_2_INDEX",9) == 0) {
			  ns = sscanf(value,"%d",&s_2_index);
			  bad_scan = (ns == 0);
			} else {
			  if (strncmp(key,"RAD_SKOSH",9) == 0) {
			    ns = sscanf(value,"%le",&rad_skosh);
			    bad_scan = (ns == 0);
			  } else {
			    if (strncmp(key,"RF_THRESH",9) == 0) {
			      ns = sscanf(value,"%le",&rf_thresh);
			      bad_scan = (ns == 0);
			    } else {
			      if (strncmp(key,"S026_FILE",9) == 0) {
				if (ival_len <= filename_length) {
				  ns = sscanf(value,"%s",s026_file);
				  bad_scan = (ns == 0);
				} else {
				  strncpy(s026_file,value,filename_length);
				  s026_file[filename_length] = '\0';
				  bad_scan = 0;
				  ns = 1;
				  fprintf(stdout,"read_params: WARNING s026_file name "
					  "truncated to %s\n",s026_file);
				  fflush(stdout);
				}
			      } else {
				if (strncmp(key,"INPUT_DIR",9) == 0) {
				  if (ival_len < longest_dir_len) {
				    ns = sscanf(value,"%s",input_dir);
				    bad_scan = (ns == 0) ;
				  } else {
				    strncpy(input_dir,value,longest_dir_len);
				    input_dir[longest_dir_len] = '\0';
				    bad_scan = 0;
				    ns = 1;
				    fprintf(stdout,"read_params: WARNING input_dir "
					    "truncated to %s\n",input_dir);
				    fflush(stdout);
				  }
				  /*
				    Ensure the nonempty directory ends in dir_sep ("/")
				  */
				  if (ns != 0) {
				    dir_len = strlen(input_dir);
				    if (dir_len > 0) {
				      dir_end = &(input_dir[dir_len-len_dir_sep]);
				      if (strcmp(dir_sep,dir_end) != 0) {
					if (dir_len < longest_dir_len) {
					  strcat(input_dir,dir_sep);
					} else {
					  dir_end = &input_dir[longest_dir_len];
					  strcat(dir_end,dir_sep);
					}
				      }
				    }
				  }
				} else {
				  if (strncmp(key,"COMBO_DIR",9) == 0) {
				    if (ival_len < longest_dir_len) {
				      ns = sscanf(value,"%s",combo_dir);
				      bad_scan = (ns == 0);
				    } else {
				      strncpy(combo_dir,value,longest_dir_len);
				      combo_dir[longest_dir_len] = '\0';
				      bad_scan = 0;
				      ns = 1;
				      fprintf(stdout,"read_params: WARNING combo_dir "
					      "truncated to %s\n",combo_dir);
				      fflush(stdout);
				    }
				    /*
				      Ensure the nonempty directory ends in dir_sep ("/")
				    */
				    if (ns != 0) {
				      dir_len = strlen(combo_dir);
				      if (dir_len > 0) {
					dir_end = &(combo_dir[dir_len-len_dir_sep]);
					if (strcmp(dir_sep,dir_end) != 0) {
					  if (dir_len < longest_dir_len) {
					    strcat(combo_dir,dir_sep);
					  } else {
					    dir_end = &combo_dir[longest_dir_len];
					    strcat(dir_end,dir_sep);
					  }
					}
				      }
				    }
				  } else {
				    if (strncmp(key,"DBG_STEPS",9) == 0) {
				      record_dbg_steps(value,num_channels,&num_dbg_steps,
						       dbg_steps,stderr);
				      if (subset_size > 0) {
					bad_scan = (num_dbg_steps < subset_size);
				      }
				    } else {
				      if (strncmp(key,"SREF_MODE",9) == 0) {
					ns = sscanf(value,"%d",&sref_mode);
					bad_scan = (ns == 0);
				      } else {
					found = 0;
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  } else {
		    if (key_len == 8) {
		      if (strncmp(key,"LOG_FILE",8) == 0) {
			if (ival_len < filename_length) {
			  ns = sscanf(value,"%s",log_file);
			  bad_scan = (ns == 0);
			} else {
			  strncpy(log_file,value,filename_length);
			  log_file[filename_length] = '\0';
			  bad_scan = 0;
			  ns       = 1;
			  fprintf(stdout,"read_params: WARNING log_file truncated "
				  "to %s\n",log_file);
			  fflush(stdout);
			}
		      } else {
			if (strncmp(key,"U_TARGET",8) == 0) {
			  ns = sscanf(value,"%le",&u_target);
			  bad_scan = (ns == 0);
			} else {
			  if (strncmp(key,"V_TARGET",8) == 0) {
			    ns = sscanf(value,"%le",&v_target);
			    bad_scan = (ns == 0);
			  } else {
			    if (strncmp(key,"TEST_RAD",8) == 0) {
			      ns = sscanf(value,"%le",&test_rad);
			      bad_scan = (ns == 0);
			    } else {
			      if (strncmp(key,"CCT_FILE",8) == 0) {
				if (ival_len < filename_length) {
				  ns = sscanf(value,"%s",cct_file);
				  bad_scan = (ns == 0);
				} else {
				  strncpy(cct_file,value,filename_length);
				  cct_file[filename_length] = '\0';
				  bad_scan = 0;
				  ns = 1;
				  fprintf(stdout,"read_params: WARNING cct_file was "
					  "truncated to %s\n",cct_file);
				  fflush(stdout);
				}
			      } else {
				if (strncmp(key,"SZ_INDEX",8) == 0) {
				  ns = sscanf(value,"%d",&sz_index);
				  bad_scan = (ns == 0);
				} else {
				  found = 0;
				}
			      }
			    }
			  }
			}
		      } 
		    } else {
		      if (strncmp(key,"USE_VRB",7) == 0) {
			ns = sscanf(value,"%d",&use_vrb);
			bad_scan = (ns == 0);
		      } else {
			if (strncmp(key,"T_T_INC",7) == 0) {
			  ns = sscanf(value,"%d",&t_t_inc);
			  bad_scan = (ns == 0);
			} else {
			  if (strncmp(key,"SCALES",6) == 0) {
			    use_multiscale = 1;
			    record_scales(value,max_scales,&number_scales,
					  scales,stderr);
			    bad_scan = (number_scales < 1);
			  } else {
			    if (strncmp(key,"SUBSET",6) == 0) {
			      record_subset(value,max_subset_size,num_channels,&subset_size,subset_indices,stderr);
			      bad_scan = (subset_size < 1);
			    } else {
			      if (strncmp(key,"HP_OUT",6) == 0) {
				ns = sscanf(value,"%d",&hp_out);
				bad_scan = (ns == 0);
			      } else {
				if (strncmp(key,"VRB",3) == 0) {
				  use_vrb = 1;
				  record_vrb(value,state);
				} else {
				  found = 0;
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  } /* end else key_len <= 13 */
	  if (found) {
	    if (bad_scan) {
	      fprintf(stderr,"read params: Invalid value in parameter file "
		      "on line %d\n%s\n",line_no,io_buffer);
	      fflush(stderr);
	    }
	  } else {
	    fprintf(stderr,"read params: Invalid keyword in parameter file "
		    "on line %d\n%s\n",line_no,io_buffer);
	    fflush(stderr);
	  }
	}
	fgp = fgets(io_buffer,buffer_length,pfp);
	line_no += 1;
      } /* end while fgp */
      fclose(pfp);
    } /* end if (pfp) */
  } /* end if (parameter_file was defined */
  /* 
    check consistencey of params.
  */
  if (source_rows != channels_rows) {
    success = 0;
    fprintf(stderr,"read_params: source_rows = %d and "
	    "channels_rows = %d, they must be the same\n",
	    source_rows,channels_rows);
    fflush(stderr);
  } else {
    if ((ybar_index != (xbar_index + 1)) || 
	(zbar_index != (ybar_index + 1))) {
      success = 0;
      fprintf(stderr,"read_params: xbar_index, ybar_index, "
	      "zbar_index = %d,%d,%d are expected to be contiguous\n",
	      xbar_index,ybar_index,zbar_index);
      fflush(stderr);
    } else {
      if ((ybar_10_index != (xbar_10_index + 1)) || 
	  (zbar_10_index != (ybar_10_index + 1))) {
	success = 0;
	fprintf(stderr,"read_params: xbar_10_index, ybar_10_index, "
	      "zbar_10_index = %d,%d,%d are expected to be contiguous\n",
	      xbar_10_index,ybar_10_index,zbar_10_index);
	fflush(stderr);
      } else {
	if ((s_1_index != (s_0_index + 1)) || 
	  (s_2_index != (s_1_index + 1))) {
	  success = 0;
	  fprintf(stderr,"read_params: s_0_index, s_1_index, "
		  "s_2_index = %d,%d,%d are expected to be contiguous\n",
		  s_0_index,s_1_index,s_2_index);
	  fflush(stderr);
	} else {
	  if (num_channels < 3) {
	    success = 0;
	    fprintf(stderr,"read_params: num_channels must be >= 3\n");
	    fflush(stderr);
	  } else {
	    if (channels_rows != source_rows) {
	      success = 0;
	      fprintf(stderr,"read_params: num_channel_rows must = "
		      "num_source_rows = %d\n",
		      source_rows);
	      fflush(stderr);
	    }
	  }
	}
      }
    }
  }
  if (success) {
    number_vrb = state->number_vrb;
    success = set_test_radii_sq(num_channels,use_vrb,number_vrb,vrb_first_index,
				vrb_last_index,vrb_radii,test_rad,
				test_radii_sq);
    if (success == 0) {
      fprintf(stderr,"read_params: vrb_first_index:vrb_last_index arrays did not span 1:num_channels in a contiguous non-overlapping way\n");
      fflush(stderr);
    }
  }
  if (nz_cpcg_max == -1) {
    nz_cpcg_max = num_channels;
  }
  /*
    Pack parameters into state.
  */
  state->source_rows          = source_rows;
  state->source_columns       = source_columns;
  state->num_samples          = num_samples;
  state->channels_rows        = channels_rows;
  state->channels_columns     = channels_columns;
  state->first_channel_index  = first_channel_index;
  state->last_channel_index   = last_channel_index;
  state->num_channels         = num_channels;
  state->full_num_channels    = num_channels;
  if (subset_size == 0) {
    /*
      No SUBSET line in parameter so set subset size to
      be num_channels to indicate a run against all channels.
    */
    state->subset_size = num_channels;
  } else {
    state->subset_size = subset_size;
  }
  state->num_cct_rows         = num_cct_rows;
  state->wavelength_index     = wavelength_index;
  state->xbar_index           = xbar_index;
  state->ybar_index           = ybar_index;
  state->zbar_index           = zbar_index;
  state->xbar_10_index        = xbar_10_index;
  state->ybar_10_index        = ybar_10_index;
  state->zbar_10_index        = zbar_10_index;
  state->s_0_index            = s_0_index;
  state->s_1_index            = s_1_index;
  state->s_2_index            = s_2_index;
  state->num_radiator_temps   = num_radiator_temps;
  state->num_huebins          = num_huebins;
  state->u_target             = u_target;
  state->v_target             = v_target;
  state->test_rad             = test_rad;
  state->rad_skosh            = rad_skosh;
  state->rf_thresh            = rf_thresh;
  state->use_multiscale       = use_multiscale;
  state->use_bruteforce       = use_bruteforce;
  state->select_only          = select_only;
  state->number_scales        = number_scales;
  state->extra_log_info       = extra_log_info;
  state->full_output          = full_output;
  state->number_vrb           = number_vrb;
  state->use_vrb              = use_vrb;
  state->sz_index             = sz_index;
  state->use_alt_cct          = use_alt_cct;
  state->s026_columns         = s026_columns;
  state->max_num_bins         = max_num_bins;
  state->nz_cpcg_max          = nz_cpcg_max;
  state->neutral_rcs_radius   = neutral_rcs_radius;  /* This should be a parameter probably */
  state->neutral_skew_radius  = neutral_skew_radius; /* This should be a parameter probably */
  /*
  state->num_quadrangles      = num_quadrangles;
  */
  state->discount_illuminant  = discount_illuminant;
  state->new_cam02ucs         = new_cam02ucs;
  state->hp_out               = hp_out;
  state->save_steps           = save_steps;
  state->read_combos          = read_combos;
  state->watchme              = watchme;
  state->read_mcat02inv       = read_mcat02inv;
  state->cam02ucs_1mvp        = cam02ucs_1mvp;
  state->cam02ucs_2mvp        = cam02ucs_2mvp;
  state->sref_mode            = sref_mode;
  state->t_t_inc              = t_t_inc;
  state->sref_t_t_min         = sref_t_t_min;
  state->sref_t_t_max         = sref_t_t_max;
  /*
    Make sure dump_combos is 0 if select_only is not set.
  */
  if (select_only == 0) {
    dump_combos = 0;
  }
  state->dump_combos          = dump_combos;
  return(success);
}
