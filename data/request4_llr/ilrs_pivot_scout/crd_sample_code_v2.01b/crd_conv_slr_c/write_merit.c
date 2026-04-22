#include <stdio.h>
#include "merit.h"

// check whether number is a standin for "na".
int iisna (int);
int disna (double);
int ldisna (long double);

// Write different types of variables to output.
static int writei_def (FILE *str, char *fmt, int, int);
static int writed_def (FILE *str, char *fmt, double, double);
static int writeld_def (FILE *str, char *fmt, long double, long double);

/*-------------------------------------------------------------------------
 * Subroutine:  write old ILRS (CSTG) fullrate data format records
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   Nov 06, 2007 - Initial version
 *   Dec 15, 2020 - CRD v2 version for conversion of "na" representation to "0".
 *                  NOTE: include ../common_c/write_crd.c to capture *isna
 *                  routines! Be sure to pick the best numberical equivalents
 *                  to "NA" in write*0 for your system!
 *                  rlr.
 *
**-----------------------------------------------------------------------*/
/*-------------------------------------------------------------------------
**
**      write_merit_fr - write old ILRS fullrate format records (previously
**			 known as the Merit II format)
**
**-----------------------------------------------------------------------*/
void
write_merit_fr (FILE *str_out, struct merit_fr fr)
{
  int nstat= 0;

  if (writei_def (str_out, "%7d", fr.ilrs_id, 0) < 0) nstat--;
  if (writei_def (str_out, "%02d", fr.year%100, 0) < 0) nstat--;
  if (writei_def (str_out, "%3d", fr.doy, 0) < 0) nstat--;
  if (writeld_def (str_out, "%12.0Lf", fr.sec_of_day, 0.0e0) < 0) nstat--;
  if (writei_def (str_out, "%4d", fr.cdp_pad_id, 0) < 0) nstat--;
  if (writei_def (str_out, "%2d", fr.cdp_sys_num, 0) < 0) nstat--;
  if (writei_def (str_out, "%2d", fr.cdp_occ_num, 0) < 0) nstat--;

  // Note: the following are doubles that need to be written as ints. 
  //   However, must not cast them to int here, as will not recognize
  //   "na" numerical values, as is different for doubles and ints!
  if (writed_def (str_out, "%7.0f", fr.azimuth, 0.0) < 0) nstat--;
  if (writed_def (str_out, "%6.0f", fr.elevation, 0.0) < 0) nstat--;
  if (writeld_def (str_out, "%12.0Lf", fr.time_of_flight, 0.0e0) < 0) nstat--;
  if (writed_def (str_out, "%7.0f", fr.sess_rms, 0.0) < 0) nstat--;
  if (writed_def (str_out, "%4.0f", fr.xmit_wavelength, 0.0) < 0) nstat--;
  if (writed_def (str_out, "%5.0f", fr.pressure, 0.0) < 0) nstat--;
  if (writed_def (str_out, "%4.0f", fr.temperature, 0.0) < 0) nstat--;
  if (writed_def (str_out, "%3.0f", fr.humidity, 0.0) < 0) nstat--;
  if (writed_def (str_out, "%5.0f", fr.refraction_corr, 0.0) < 0) nstat--;
  if (writed_def (str_out, "%6.0f", fr.target_CofM_corr, 0.0) < 0) nstat--;
  if (writed_def (str_out, "%5.0f", fr.xcv_amp, 0.0) < 0) nstat--;
  if (writed_def (str_out, "%8.0f", fr.cal_sys_delay, 0.0) < 0) nstat--;
  if (writed_def (str_out, "%6.0f", fr.cal_delay_shift, 0.0) < 0) nstat--;
  if (writed_def (str_out, "%4.0f", fr.cal_rms, 0.0) < 0) nstat--;

  if (writei_def (str_out, "%1d", fr.np_window_ind, 0) < 0) nstat--;
  if (writei_def (str_out, "%4d", fr.num_ranges, 0) < 0) nstat--;
  if (writei_def (str_out, "%1d", fr.epoch_event, 0) < 0) nstat--;
  if (writei_def (str_out, "%1d", fr.stn_timescale, 9) < 0) nstat--;
  if (writei_def (str_out, "%1d", fr.angle_origin_ind, 9) < 0) nstat--;
  if (writei_def (str_out, "%1d", fr.refraction_app_ind, 9) < 0) nstat--;
  if (writei_def (str_out, "%1d", fr.CofM_app_ind, 9) < 0) nstat--;
  if (writei_def (str_out, "%1d", fr.xcv_amp_app_ind, 9) < 0) nstat--;
  if (writei_def (str_out, "%1d", fr.cal_type_ind, 9) < 0) nstat--;
  if (writei_def (str_out, "%1d", fr.sys_change_ind, 9) < 0) nstat--;
  if (writei_def (str_out, "%1d", fr.sys_config_ind, 9) < 0) nstat--;
  if (writei_def (str_out, "%1d", fr.format_version, 9) < 0) nstat--;
  if (writei_def (str_out, "%1d\n", fr.data_release, 0) < 0) nstat--;

  if (nstat < 0) printf ("Write error: write_merit_fr.c\n");

// 7603901203556562483324697840350100000000000000533269598190000079532110130279109100000000000000000010445300000000237068704011100020
}

// Write an integer number or a number that substitutes for "na"
// Input:
//   str - stream id for output file
//   fmt - format used to output this field
//   iv  - integer variable to write to file
//   ivna- integer variable to write if iv is the numerical representation
//         of "na" (not available or not applicable)
// Output:
//   a character string is written to the file attached to "str".
//
int
writei_def (FILE *str, char *fmt, int iv, int ivna)
{
  if (!iisna (iv))
    {
      fprintf (str, fmt, iv);
    }
  else
    {
      fprintf (str, fmt, ivna);
    }
}

// Write a double precisions number or a number that substitutes for "na"
// Input:
//   str - stream id for output file
//   fmt - format used to output this field
//   iv  - double precision variable to write to file
//   ivna- double precision variable to write if iv is the numerical 
//         representation
//         of "na" (not available or not applicable)
// Output:
//   a character string is written to the file attached to "str".
int
writed_def (FILE *str, char *fmt, double dv, double dvna)
{
  if (!disna (dv))
    {
      fprintf (str, fmt, dv);
    }
  else
    {
      fprintf (str, fmt, dvna);
    }
}

// Write a long double precisions number or a number that substitutes for "na"
// Input:
//   str - stream id for output file
//   fmt - format used to output this field
//   iv  - long double precision variable to write to file
//   ivna- long double precision variable to write if iv is the numerical 
//         representation
//         of "na" (not available or not applicable)
// Output:
//   a character string is written to the file attached to "str".
int
writeld_def (FILE *str, char *fmt, long double ldv, long double ldvna)
{
  if (!ldisna (ldv))
    {
      fprintf (str, fmt, ldv);
    }
  else
    {
      fprintf (str, fmt, ldvna);
    }
}

