#include "system_includes.h"
#include "loop_state_struct.h"
#include "gen_x_combo.h"
int gen_x_combo(struct loop_state_struct *loop_state, double *xcombo) {
  /*
    Get the next combination of x (red) channels.
    Called by: gen_combos.
    Returns 1 if a new red channel combination was found, 0 otherwise.

    Arguments :             TMF     Description
       loop_state           p*b     loop_state structure.
       
    Fields of loop_state used:
       ny                   number of y-dominant channels
       nz                   number of z-dominant channels
       num_channels         total number of channels - 
                            number of x-dominant channels = num_channels-ny-nz

       max_index,           maximum value in istep_index array.

       prev_xyz,            (3 x num_channels),
                            the previous tristimulus values for
                            channels at the power levels specified in
			    istep_index, locally used as prev_xyzs

       stepped_channel_xyz, (3 x total steps) array with the scaled
                            x,y,z tristimulus values at each power level,
			    for each channel. Indexed by istepped_xyz_ptrs.


       istepped_xyz_ptrs,   Index vector of length nun_channels pointing to
                            the first element of the tristimulus power level values for
			    each channel, locally used as is_xyz_p

       istep_index,         array of index values (length is num_channels)
                            the index values for the current z, current y,
			    and previous x combo.

       last_index           maximum number of steps per channel

       lfp                  log file pointer.			    
			    
    Fields of loop_state set:
       prev_xyz,            tristimulus values for
			    channels at the power levels specified in
			    istep_index, locally used as prev_xyzs
   
       curr_xyz,            (x,y,z) values for new combination of red channels
                            combined with current blue and green channel combo.

       istep_index          array of index values (length is num_channels)
                            the index values for current z and 
			    the next y combo 
			    istep_index[ny+nz:num_channels-1] are modified by
			    this routine.

       last_index           maximum number of steps per channel
                            last_index[ny+nz:num_channels-1] are modified by
			    this routine

  */
  double *stepped_channel_xyz;
  double *prev_xyzs;
  double *curr_xyz;

  double *channel_steps;
  double *prev_xyz;
  double *prev_xyzf;
  double *xyz_sort;
  double *vb_plane;
  double *ub_plane;

  double delta_x;
  double delta_y;
  double delta_z;
  double fract;

  int *is_xyz_p;
  int *istep_index;
  int *msteps_per_channel;
  
  int max_index;

  int x_combos;
  int ix;

  int i3x;
  int interval_ix;

  int istep;
  int iix;

  int i3step;

  int ny;
  int nz;

  int num_channels;
  int more_x_combos;

  int step_taken;
  int irix;

  int istepd;
  int istep_min;

  int irl;
  int is3m;

  int inc1;
  int ithree;

  FILE *lfp;
  FILE *efp;

  inc1                 = 1;
  ithree               = 3;
  ny                   = loop_state->ny;
  nz                   = loop_state->nz;
  num_channels         = loop_state->num_channels;
  max_index            = loop_state->max_index;
  stepped_channel_xyz  = loop_state->stepped_channel_xyz;
  prev_xyzs            = loop_state->prev_xyz;
  msteps_per_channel   = loop_state->msteps_per_channel;
  is_xyz_p             = loop_state->istepped_xyz_ptrs;
  istep_index          = loop_state->istep_index;
  xyz_sort             = loop_state->xyz_sort;
  lfp                  = loop_state->lfp;
  vb_plane             = loop_state->vb_plane;
  ub_plane             = loop_state->ub_plane;

  fract                = 1.0/((double)max_index);
  x_combos = 0;
  for (ix = num_channels - 1; ((ix >= ny+nz) && (x_combos == 0)); ix--) {
    i3x = ix + ix + ix;
    if (istep_index[ix] < msteps_per_channel[ix]-1) {
      x_combos        = 1;
      istep_index[ix] += 1;
      istep           = istep_index[ix];
      /*
	((istep > 0) && (istep <= last_index[ix]))
      */
      i3step          = istep + istep + istep;
      prev_xyz        = &prev_xyzs[i3x];
      channel_steps   = &stepped_channel_xyz[is_xyz_p[ix]];
      delta_x         = channel_steps[i3step];
      delta_y         = channel_steps[i3step+1];
      delta_z         = channel_steps[i3step+2];
      xcombo[0]       = prev_xyz[0] + delta_x;
      xcombo[1]       = prev_xyz[1] + delta_y;
      xcombo[2]       = prev_xyz[2] + delta_z;
      /*
      xcombo[3] = ddot_(&ithree,xcombo,&inc1,vb_plane,&inc1);
      xcombo[4] = ddot_(&ithree,xcombo,&inc1,ub_plane,&inc1);
      */
      xcombo[3] = ((xcombo[0] * vb_plane[0]) + (xcombo[1] * vb_plane[1])) +
   	          (xcombo[2] * vb_plane[2]);
      xcombo[4] = ((xcombo[0] * ub_plane[0]) + (xcombo[1] * ub_plane[1])) +
   	          (xcombo[2] * ub_plane[2]);

      prev_xyzf = prev_xyz + 3; /* Caution address arithmetic */
      for (iix = ix+1;iix < num_channels;iix++) {
	prev_xyzf[0] = xcombo[0];
	prev_xyzf[1] = xcombo[1];
	prev_xyzf[2] = xcombo[2];
	prev_xyzf    = prev_xyzf + 3; /* Caution address arithmetic */
      }
    } else {
      istep_index[ix] = 0;
    }
  } /* end for (ix ...) */ 
  return(x_combos);
}

