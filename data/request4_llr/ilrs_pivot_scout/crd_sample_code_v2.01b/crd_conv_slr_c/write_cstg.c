#include <stdio.h>
#include <string.h>
#include "cstg.h"


// check whether number is a standin for "na".
int iisna (int);
int disna (double);
int ldisna (long double);

static int first_sed= 1;
int check_sum ();

// Write different types of variables to output.
static int writei_defcs (FILE *str, char *fmt, int, int, int *);
static int writed_defcs (FILE *str, char *fmt, double, double, int *);
static int writeld_defcs (FILE *str, char *fmt, long double, long double, int *);

/*-------------------------------------------------------------------------
 * Subroutines:  write old ILRS (CSTG) normalpoint and sampled
 *               engineering data format records
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   Nov 06, 2007 - Initial version
 *   Sep 18, 2018 - added "void" to routine declarations;
 *                  expanded str in write_cstg_sed from 70 -> 71 long. rlr.
 *   Feb 02, 2021 - Rewrite to handle "na" in CRD v2. rlr.
 *
**-----------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
**
**      write_cstg_hdr - write header record for old ilrs (cstg)
**                       normalpoint and sampled engineering format
**
**-----------------------------------------------------------------------*/
void
write_cstg_hdr (FILE *str_out, struct cstg_hdr hdr, int file_type)
{
  char str[60];
  int cs= 0;
  int nstat= 0;

  if (file_type == 0)
    fprintf(str_out, "88888\n");
  else if (file_type == 1)
    fprintf(str_out, "99999\n");

  // Write fiel-by-field fir CRD v2 to properly handle "na" fields.
  // The final value in each field is what YOU want to be substituted for the 
  // "NA" value. It's normally assumed to be zero here.
  if (writei_defcs (str_out, "%07d", hdr.ilrs_id, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%02d", hdr.year, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%03d", hdr.doy, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%04d", hdr.cdp_pad_id, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%02d", hdr.cdp_sys_num, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%02d", hdr.cdp_occ_num, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%04.0f", hdr.xmit_wavelength, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%08.0f", hdr.cal_sys_delay, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%06.0f", hdr.cal_delay_shift, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%04.0f", hdr.cal_rms, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%01d", hdr.np_window_ind, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%01d", hdr.stn_timescale, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%01d", hdr.cal_type_ind, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%01d", hdr.sys_change_ind, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%01d", hdr.sys_config_ind, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%04.0f", hdr.sess_rms, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%01d", hdr.data_qual_ind, 0, &cs) < 0) nstat--;
  hdr.checksum= cs % 100;
  if (writei_defcs (str_out, "%02d", hdr.checksum, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%01d\n", hdr.format_version, 0, &cs) < 0) nstat--;
}

/*-------------------------------------------------------------------------
**
**      write_cstg_np - write old ilrs (cstg) normalpoint record
**
**-----------------------------------------------------------------------*/
void
write_cstg_np (FILE *str_out, struct cstg_np np)
{
  char str[60];
  int cs= 0;
  int nstat= 0;

  if (writeld_defcs (str_out, "%012.0Lf", np.sec_of_day, 0.0e0, &cs) < 0) nstat--;
  if (writeld_defcs (str_out, "%012.0Lf", np.time_of_flight, 0.0e0, &cs) < 0) nstat--;

  // Note: the following are doubles that need to be written as ints.
  //   However, must not cast them to int here, as will not recognize
  //   "na" numerical values, as is different for doubles and ints!
  if (writed_defcs (str_out, "%07.0f", np.bin_rms, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%05.0f", np.pressure, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%04.0f", np.temperature, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%03.0f", np.humidity, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%04d", np.num_ranges, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%01d", np.data_release, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%01d", np.scale_or_tof_sec, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%01d", np.llr_np_window_ind, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%02.0f", np.snr, 0, &cs) < 0) nstat--;
  np.checksum= cs % 100;
  if (writei_defcs (str_out, "%02d\n", np.checksum, 0, &cs) < 0) nstat--;
}

/*-------------------------------------------------------------------------
**
**      write_cstg_sed - write old ilrs (cstg) sampled engineering record
**
**-----------------------------------------------------------------------*/
void
write_cstg_sed (FILE *str_out, struct cstg_sed sed)
{
  char str[71];
  int cs= 0;
  int nstat= 0;

  if (writeld_defcs (str_out, "%012.0Lf", sed.sec_of_day, 0.0e0, &cs) < 0) nstat--;
  if (writeld_defcs (str_out, "%012.0Lf", sed.time_of_flight, 0.0e0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%05.0f", sed.pressure, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%04.0f", sed.temperature, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%03.0f", sed.humidity, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%08.0f", sed.internal_burst_cal, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%04d", sed.xcv_amp, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%01d", sed.angle_origin_ind, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%07.0f", sed.azimuth, 0, &cs) < 0) nstat--;
  if (writed_defcs (str_out, "%06.0f", sed.elevation, 0, &cs) < 0) nstat--;
  if (writei_defcs (str_out, "%05d", 0, 0, &cs) < 0) nstat--;
  sed.checksum= cs % 100;
  if (writei_defcs (str_out, "%02d\n", sed.checksum, 0, &cs) < 0) nstat--;
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
writei_defcs (FILE *str, char *fmt, int iv, int ivna, int *cs)
{
  char test[20];

  if (!iisna (iv))
    {
      fprintf (str, fmt, iv);
      sprintf (test, fmt, iv);
    }
  else
    {
      fprintf (str, fmt, ivna);
      sprintf (test, fmt, ivna);
    }
  (*cs)+= check_sum(test, strlen(test)); 
}

// Write a double precisions number or a number that substitutes for "na"
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
writed_defcs (FILE *str, char *fmt, double dv, double dvna, int *cs)
{
  char test[20];

  if (!disna (dv))
    {
      fprintf (str, fmt, dv);
      sprintf (test, fmt, dv);
    }
  else
    {
      fprintf (str, fmt, dvna);
      sprintf (test, fmt, dvna);
    }
  (*cs)+= check_sum(test, strlen(test)); 
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
writeld_defcs (FILE *str, char *fmt, long double ldv, long double ldvna, int *cs)
{
  char test[20];

  if (!ldisna (ldv))
    {
      fprintf (str, fmt, ldv);
      sprintf (test, fmt, ldv);
    }
  else
    {
      fprintf (str, fmt, ldvna);
      sprintf (test, fmt, ldvna);
    }
  (*cs)+= check_sum(test, strlen(test)); 
}


/*********************************************************************
*
*PURPOSE
*       Determine and return the checksum for an output field.
*
*********************************************************************/
int
check_sum(char *line,int line_len)
{
        char check[3];
        int i, sum;

        sum = 0;
        for (i=0;i<line_len;i++)
	  {
             if (line[i] == '1') sum+= 1;
             if (line[i] == '2') sum+= 2;
             if (line[i] == '3') sum+= 3;
             if (line[i] == '4') sum+= 4;
             if (line[i] == '5') sum+= 5;
             if (line[i] == '6') sum+= 6;
             if (line[i] == '7') sum+= 7;
             if (line[i] == '8') sum+= 8;
             if (line[i] == '9') sum+= 9;
          }
	return (sum);
}
