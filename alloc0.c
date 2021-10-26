#include "system_includes.h"
#include "state_struct.h"
#include "alloc0.h"
int alloc0(struct state_struct **state_p) {
  /*
    Allocate space for the state struct and its filename subfields,
    and the io_buffer fields.
    
    Allocates state and allocates and/or sets the following fields of state:
      
    scales
    vrb_radii
    test_radia_sq
    dbg_steps
    parameter_file
    log_file
    source_file
    channels_file
    cct_file
    huebins_files
    output_file
    z_dom_file
    y_dom_file
    x_dom_file
    s26_file
    subset_file
    summary_file
    combos_in_file
    combos_out_file
    input_dir
    combo_dir
    output_dir
    dir_plus_fn,
    mcat02inv_file,
    io_buffer,
    isteps_per_channel
    vrb_first_index
    vrb_last_index
    subset_indices

    buffer_length
    filename_length
    max_scales
    max_num_channels
    max_radii_bins
    max_subsetsize
    usage
    
    Returns 1 on success 0 on failure.

    If successful, sets the usage field of the structure

    Called by: thous_opt
    Calls:     calloc
  */
  struct state_struct *state;
  struct state_struct state_i;
  double *scales;            /* 16 */
  double *vrb_radii;         /* 64 */
  double *test_radii_sq;     /* 64 */
  double *dbg_steps;         /* 64 */
  char *parameter_file;      /* 128 */
  char *log_file;            /* 128 */
  char *source_file;         /* 128 */
  char *channels_file;       /* 128 */
  char *cct_file;            /* 128 */
  char *huebins_file;        /* 128 */
  char *output_file;         /* 128 */
  char *z_dom_file;           /* 128 */
  char *y_dom_file;          /* 128 */
  char *x_dom_file;            /* 128 */
  char *s026_file;           /* 128 */
  char *subset_file;         /* 128 */
  char *summary_file;        /* 128 */
  char *combos_in_file;      /* 128 */
  char *combos_out_file;     /* 128 */
  char *input_dir;           /* 128 */
  char *combo_dir;           /* 128 */
  char *output_dir;          /* 128 */
  char *mcat02inv_file;      /* 128 */
  char *dir_plus_fn;         /* 256 */
  char *io_buffer;           /* 16384 */
  int  *isteps_per_channel;  /* 64 */
  int  *vrb_first_index;     /* 64 */
  int  *vrb_last_index;      /* 64 */
  int  *subset_indices;      /* 64 */
  
  int64_t ask_for;
  int64_t one_l;
  int64_t usage;
  int64_t filename_length;
  int64_t buffer_length;

  int success;
  int num_files;

  int max_num_channels;
  int max_scales;

  int max_radii_bins;
  int max_radii_bins_p1;

  int max_subset_size;

  success           = 1;
  num_files         = 32; /* >= 18 for files and directories + 2 for dir_plus_fn work */
  one_l             = (int64_t)1;
  max_scales        = 16;
  max_radii_bins    = 64;
  max_num_channels  = 64;
  max_radii_bins_p1 = max_radii_bins+1;
  max_subset_size   = 64;
  buffer_length = (int64_t)16384;
  filename_length = (int64_t)128;
  ask_for = sizeof(state_i);
  usage   = ask_for;
  state = (struct state_struct *)calloc(one_l,ask_for);
  if (state == NULL) {
    success = 0;
  } 
  if (success) {
    *state_p = state;
    state->filename_length = filename_length;
    ask_for = num_files * filename_length;
    usage   += ask_for;
    parameter_file = (char *)calloc(one_l,ask_for);
    if (parameter_file == NULL) {
      success = 0;
    }
  }
  if (success) {
    /* 
       Caution the following 6 statments use address arithmetic 
    */
    log_file = parameter_file + filename_length;    
    source_file      = log_file + filename_length;
    channels_file    = source_file + filename_length;
    huebins_file     = channels_file + filename_length;
    output_file      = huebins_file + filename_length;
    cct_file         = output_file + filename_length;
    z_dom_file       = cct_file + filename_length;
    y_dom_file       = z_dom_file + filename_length;
    x_dom_file       = y_dom_file + filename_length;
    s026_file        = x_dom_file + filename_length;
    subset_file      = s026_file + filename_length;
    summary_file     = subset_file + filename_length;
    combos_in_file   = summary_file + filename_length;
    combos_out_file  = combos_in_file + filename_length;
    input_dir        = combos_out_file + filename_length;
    combo_dir        = input_dir + filename_length;
    output_dir       = combo_dir + filename_length;
    mcat02inv_file   = output_dir + filename_length;
    dir_plus_fn      = mcat02inv_file + filename_length;


    state->filename_length  = (int)filename_length;
    state->parameter_file   = parameter_file;
    state->log_file         = log_file;
    state->source_file      = source_file;
    state->channels_file    = channels_file;
    state->huebins_file     = huebins_file;
    state->output_file      = output_file;
    state->cct_file         = cct_file;
    state->z_dom_file       = z_dom_file;
    state->y_dom_file       = y_dom_file;
    state->x_dom_file       = x_dom_file;
    state->s026_file        = s026_file;
    state->subset_file      = subset_file;
    state->summary_file     = summary_file;
    state->combos_in_file   = combos_in_file;
    state->combos_out_file  = combos_out_file;
    state->input_dir        = input_dir;
    state->combo_dir        = combo_dir;
    state->output_dir       = output_dir;
    state->dir_plus_fn      = dir_plus_fn;
    state->mcat02inv_file   = mcat02inv_file;
    state->buffer_length    = buffer_length;

    ask_for = buffer_length;
    usage += ask_for;
    io_buffer = (char*)calloc(one_l,ask_for);
    if (io_buffer == NULL) {
      success = 0;
    } else {
      state->io_buffer        = io_buffer;
    }
  }
  if (success) {
    state->max_num_channels = max_num_channels;
    ask_for = (int64_t)(max_num_channels * sizeof(int));
    usage += ask_for;
    isteps_per_channel = (int*)calloc(one_l,ask_for);
    if (isteps_per_channel == NULL) {
      success = 0;
    } else {
      state->isteps_per_channel = isteps_per_channel;
    }
  }
  if (success) {
    state->max_scales         = max_scales;
    ask_for = (int64_t)(max_scales*sizeof(double));
    usage += ask_for;
    scales                    = (double*)calloc(one_l,ask_for);
    if (scales == NULL) {
      success = 0;
    } else {
      state->scales = scales;
    }
  }
  if (success) {
    state->max_radii_bins = max_radii_bins;
    ask_for = (int64_t)(max_radii_bins_p1 * sizeof(double));
    usage   += ask_for;
    vrb_radii                 = (double*)calloc(one_l,ask_for);
    if (vrb_radii == NULL) {
      success = 0;
    } else {
      state->vrb_radii = vrb_radii;
    }
  }
  if (success) {
    ask_for = (int64_t)(max_radii_bins_p1 * sizeof(double));
    usage   += ask_for;
    test_radii_sq              = (double*)calloc(one_l,ask_for);
    if (test_radii_sq == NULL) {
      success = 0;
    } else {
      state->test_radii_sq = test_radii_sq;
    }
  }
  if (success) {
    ask_for = (int64_t)(max_radii_bins_p1 * sizeof(int));
    usage   += ask_for;
    vrb_first_index           = (int*)calloc(one_l,ask_for);
    if (vrb_first_index == NULL) {
      success = 0;
    } else {
      state->vrb_first_index    = vrb_first_index;
    }
  }
  if (success) {
    ask_for = (int64_t)(max_radii_bins_p1 * sizeof(int));
    usage   += ask_for;
    vrb_last_index            = (int*)calloc(one_l,ask_for);
    if (vrb_last_index == NULL) {
      success = 0;
    } else {
      state->vrb_last_index     = vrb_last_index;
    }
  }
  if (success) {
    ask_for = (int64_t)(max_subset_size * sizeof(int));
    usage += ask_for;
    subset_indices            = (int*)calloc(one_l,ask_for);
    if (subset_indices == NULL) {
      success = 0;
    } else {
      state->subset_indices     = subset_indices;
    }
  }
  if (success) {
    ask_for = (int64_t)(max_num_channels * sizeof(double));
    usage += ask_for;
    dbg_steps = (double*)calloc(one_l,ask_for);
    if (dbg_steps == NULL) {
      success = 0;
    } else {
      state->dbg_steps = dbg_steps;
    }
  }
  if (success) {
    state->max_subset_size    = max_subset_size;
    state->usage              = usage;
  }
  return(success);
}
