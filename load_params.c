#include "system_includes.h"
#include "state_struct.h"
#include "load_params.h"
void load_params(int success, struct state_struct *state, int *iparams,
		 double *rparams){
  /*
    Extract the iparams and rparams vectors from state and success.
    Called by: rthous_opt
    Uses the following fields of state:
      num_procs           
      num_threads 	     
      source_rows         
      source_columns	     
      channels_rows       
      channels_columns    
      first_channel_index 
      last_channel_index  
      num_channels        
      num_huebins         
      buffer_length       
      num_samples         
      wavelength_index    
      xbar_index	     
      ybar_index          
      zbar_index	     
      xbar_10_index	     
      ybar_10_index	     
      zbar_10_index	     
      s_0_index 	     
      s_1_index	     
      s_2_index	     
      num_radiator_temps  
      max_num_channels    
      nc_stesp_set	     
      num_cct_rows        
      max_scales
      number_scales
      use_multiscale
      extra_log_info
      use_vrb
      number_vrb
      max_index
      sz_index
      s026_columns
      max_num_bins
      nz_cpcg_max
      max_subset_size
      subset_size
      full_num_channels
      full_output
      select_only
      use_bruteforce
      discount_illuminant
      new_cam02ucs
      hp_out
      read_combos
      save_steps

      u_target            
      v_target            
      test_rad            
      rf_thresh           
      rad_skosh                 
  */
  iparams[0]  = success;
  iparams[1]  = state->num_procs;
  iparams[2]  = state->num_threads;
  iparams[3]  = state->source_rows;
  iparams[4]  = state->source_columns;
  iparams[5]  = state->channels_rows;
  iparams[6]  = state->channels_columns;
  iparams[7]  = state->first_channel_index;
  iparams[8]  = state->last_channel_index;
  iparams[9]  = state->num_channels;
  iparams[10] = state->num_huebins;
  iparams[11] = state->buffer_length;
  iparams[12] = state->num_samples;
  iparams[13] = state->wavelength_index;
  iparams[14] = state->xbar_index;
  iparams[15] = state->ybar_index;
  iparams[16] = state->zbar_index;
  iparams[17] = state->xbar_10_index;
  iparams[18] = state->ybar_10_index;
  iparams[19] = state->zbar_10_index;
  iparams[20] = state->s_0_index;
  iparams[21] = state->s_1_index;
  iparams[22] = state->s_2_index;
  iparams[23] = state->num_radiator_temps;
  iparams[24] = state->max_num_channels;
  iparams[25] = state->nc_steps_set;
  iparams[26] = state->num_cct_rows;
  iparams[27] = state->max_scales;
  iparams[28] = state->number_scales;
  iparams[29] = state->use_multiscale;
  iparams[30] = state->extra_log_info;
  iparams[31] = state->use_vrb;
  iparams[32] = state->number_vrb;
  iparams[33] = state->max_index;
  iparams[34] = state->sz_index;
  iparams[35] = state->s026_columns;
  iparams[36] = state->max_num_bins;
  iparams[37] = state->nz_cpcg_max;
  iparams[38] = state->max_subset_size;
  iparams[39] = state->subset_size;
  iparams[40] = state->full_num_channels;
  iparams[41] = state->full_output;
  iparams[42] = state->select_only;
  iparams[43] = state->use_bruteforce;
  iparams[44] = state->discount_illuminant;
  iparams[45] = state->new_cam02ucs;
  iparams[46] = state->hp_out;
  iparams[47] = state->read_combos;
  iparams[48] = state->save_steps;
  iparams[49] = state->use_alt_cct;
  iparams[50] = state->read_mcat02inv;
  iparams[51] = state->cam02ucs_1mvp;
  iparams[52] = state->cam02ucs_2mvp;
  iparams[53] = state->sref_mode;
  iparams[54] = state->t_t_inc;
  iparams[55] = state->sref_t_t_min;
  iparams[56] = state->sref_t_t_max;
  iparams[57] = state->dump_combos;
  
  rparams[0]  = state->u_target;
  rparams[1]  = state->v_target;
  rparams[2]  = state->test_rad;
  rparams[3]  = state->rf_thresh;
  rparams[4]  = state->rad_skosh;
  rparams[5]  = state->watchme;
}
