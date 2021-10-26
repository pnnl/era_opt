#include "system_includes.h"
#include "state_struct.h"
#include "unload_params.h"
void unload_params(int *success_p, struct state_struct *state, int *iparams,
		 double *rparams){
  /*
    Extract success and state scalars from the iparams and rparams vectors.
    Called by: rthous_opt

    Sets the following fields of state:
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
      use_brute_force
      discount_illuminant
      new_cam02ucs

      u_target            
      v_target            
      test_rad            
      rf_thresh           
      rad_skosh                 
  */
  *success_p                 = iparams[0];
  state->num_procs           = iparams[1];
  state->num_threads 	     = iparams[2];
  state->source_rows         = iparams[3];
  state->source_columns	     = iparams[4];
  state->channels_rows       = iparams[5];
  state->channels_columns    = iparams[6];
  state->first_channel_index = iparams[7];
  state->last_channel_index  = iparams[8];
  state->num_channels        = iparams[9];
  state->num_huebins         = iparams[10];
  state->buffer_length       = iparams[11];
  state->num_samples         = iparams[12];
  state->wavelength_index    = iparams[13];
  state->xbar_index	     = iparams[14];
  state->ybar_index          = iparams[15];
  state->zbar_index	     = iparams[16];
  state->xbar_10_index	     = iparams[17];
  state->ybar_10_index	     = iparams[18];
  state->zbar_10_index	     = iparams[19];
  state->s_0_index 	     = iparams[20];
  state->s_1_index	     = iparams[21];
  state->s_2_index	     = iparams[22];
  state->num_radiator_temps  = iparams[23];
  state->max_num_channels    = iparams[24];
  state->nc_steps_set	     = iparams[25];
  state->num_cct_rows        = iparams[26];
  state->max_scales          = iparams[27];
  state->number_scales       = iparams[28];
  state->use_multiscale      = iparams[29];
  state->extra_log_info      = iparams[30];
  state->use_vrb             = iparams[31];
  state->number_vrb          = iparams[32];
  state->max_index           = iparams[33];
  state->sz_index            = iparams[34];
  state->s026_columns        = iparams[35];
  state->max_num_bins        = iparams[36];
  state->nz_cpcg_max         = iparams[37];
  state->max_subset_size     = iparams[38];
  state->subset_size         = iparams[39];
  state->full_num_channels   = iparams[40];
  state->full_output         = iparams[41];
  state->select_only         = iparams[42];
  state->use_bruteforce      = iparams[43];
  state->discount_illuminant = iparams[44];
  state->new_cam02ucs        = iparams[45];
  state->hp_out              = iparams[46];
  state->read_combos         = iparams[47];
  state->save_steps          = iparams[48];
  state->use_alt_cct         = iparams[49];
  state->read_mcat02inv      = iparams[50];
  state->cam02ucs_1mvp       = iparams[51];
  state->cam02ucs_2mvp       = iparams[52];
  state->sref_mode           = iparams[53];  
  state->t_t_inc	     = iparams[54];  
  state->sref_t_t_min	     = iparams[55];  
  state->sref_t_t_max	     = iparams[56];  
  state->dump_combos         = iparams[57];

  state->u_target            = rparams[0];
  state->v_target            = rparams[1];
  state->test_rad            = rparams[2];
  state->rf_thresh           = rparams[3];
  state->rad_skosh           = rparams[4];
  state->watchme             = rparams[5];
}
