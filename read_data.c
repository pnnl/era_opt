#include "system_includes.h"
#include "state_struct.h"
#include "csv_read.h"
#include "read_data.h"
int read_data(struct state_struct *state) {
  /*
    Read the source data, channels data, 
    huebins coordinates, cct_matrix, all expected to be in csv files
    Called by: era_opt
    Calls:     csv_read

    Uses the following fields of state
      channels_file,
      source_file,
      huebins_file,
      cct_file,
      source_rows,
      source_columns,
      channels_rows,
      channels_columns
      num_huebins
      num_cct_rows
      source_matrix,
      channel_matrix,
      huebins,
      cct_matrix

    Sets the following fields of state
      *channel_matrix,
      *source_matrix,
      *huebins
      *cct_matrix
  */
  double *source_matrix;
  double *channel_matrix;
  double *s026_matrix;
  /*
  double *quadrangles;
  */
  double *huebins;
  double *cct_matrix;
  char *io_buffer;
  char *source_file;
  char *channels_file;
  char *cct_file;
  char *s026_file;
  /*
  char *quadrangles_file;
  */
  char *huebins_file;
  char *input_dir;
  char *dir_plus_fn;
  char *filename_pos;

  int source_rows;
  int source_columns;

  int channels_rows;
  int channels_columns;
  
  /*
  int num_quadrangles;
  */

  int num_huebins;
  int buffer_length;

  int success;
  int i6;

  int i3;
  int num_cct_rows;

  int s026_columns;
  int idir_len;

  int select_only;
  int ipad;

  FILE *lfp;
  /* 
    Extract state fields.
  */
  success          = 1;
  channels_file    = state->channels_file;
  select_only      = state->select_only;
  channels_rows    = state->channels_rows;
  channels_columns = state->channels_columns;
  source_rows      = state->source_rows;
  source_columns   = state->source_columns;
  input_dir        = state->input_dir;
  dir_plus_fn      = state->dir_plus_fn;
  idir_len         = strlen(input_dir);
  filename_pos     = dir_plus_fn;
  lfp              = state->lfp;
  io_buffer        = state->io_buffer;
  buffer_length    = state->buffer_length;
  channel_matrix   = state->channel_matrix;
  source_file      = state->source_file;
  source_matrix    = state->source_matrix;
  if (idir_len > 0) {
    strcpy(dir_plus_fn,input_dir);
    filename_pos += idir_len; /* Caution address arithmetic */
  }
  if (success) {
    /*
      Prefix channels file with input directory
    */
    strcpy(filename_pos,channels_file);
    success = csv_read(dir_plus_fn,
		       io_buffer,
		       buffer_length,
		       channels_rows,
		       channels_columns,
		       channel_matrix,
		       lfp);
  }
  if (success) {
    /*
      Prefix source file with input directory
    */
    strcpy(filename_pos,source_file);
    success = csv_read(dir_plus_fn,
		       io_buffer,
		       buffer_length,
		       source_rows,
		       source_columns,
		       source_matrix,
		       lfp);
  }
  if ((select_only == 0) && success) {
    s026_file        = state->s026_file;
    /*
      quadrangles_file = state->quadrantgles_file;
    */
    huebins_file     = state->huebins_file;
    cct_file         = state->cct_file;
    num_cct_rows     = state->num_cct_rows;
    s026_columns     = state->s026_columns;
    /*
      num_quadrangales = state->num_quadrangales;
    */
    num_huebins      = state->num_huebins;
    cct_matrix       = state->cct_matrix;
    s026_matrix      = state->s026_matrix;
    /*
      quadrangles      = state->quadrantles;
    */
    huebins        = state->huebins;
    /*
    if (success) {
    i8 = 8;
    strcpy(filename_pos,quadrangles_file);
    success = csv_read(dir_plus_fn,
                       io_buffer,
    		         buffer_length,
    		         num_quadrangles,
    		         i8,
    		         quadrangles,
    		         lfp);
    }
    */
    if (success) {
      i6 = 6;
      /*
	Prefix huebins file with input directory
      */
      strcpy(filename_pos,huebins_file);
      success = csv_read(dir_plus_fn,
			 io_buffer,
			 buffer_length,
			 num_huebins,
			 i6,
			 huebins,
			 lfp);
    }
    if (success) {
      i3 = 3;
      /*
	Prefix cct file with input directory
      */
      strcpy(filename_pos,cct_file);
      success = csv_read(dir_plus_fn,
			 io_buffer,
			 buffer_length,
			 num_cct_rows,
			 i3,
			 cct_matrix,
			 lfp);
    }
    if (success) {
      /*
	Prefix s026 file with input directory
      */
      strcpy(filename_pos,s026_file);
      success = csv_read(dir_plus_fn,
			 io_buffer,
			 buffer_length,
			 source_rows,
			 s026_columns,
			 s026_matrix,
			 lfp);
    }
  } /* end if (select_only == 0) && success */
  return(success);
}
 
