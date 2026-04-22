#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/crd.h"

struct rh1 h1;
struct rh2 h2;
struct rh3 h3;
struct rh4 h4;
struct rh5 h5;
struct rc0 c0;
struct rc1 c1;
struct rc2 c2;
struct rc3 c3;
struct rc4 c4;
struct rc5 c5;
struct rc6 c6;
struct rc7 c7;
struct rd10 d10;
struct rd11 d11;
struct rd12 d12;
struct rd20 d20;
struct rd21 d21;
struct rd30 d30;
struct rd40 d40;
struct rd40 d41;
struct rd42 d42;
struct rd50 d50;
struct rd60 d60;
struct rd00 d00;

/*-------------------------------------------------------------------------
 * Program: crd_conv
 *
 * Purpose:
 * Converts ILRS CRD file between version 1 and 2. 
 *   If v1 is input, v2 is output
 *   If v2 is input, v1 is output
 * To prevent "collisions" with real data in testing, the CDP pad ID can
 * be replaced. Suggested values run from 9990 through 9999.
 *
 * This program can be used to test the consistency of the CRD c read/write 
 * routines by converting from one format to the other and then back.
 *
 * Calling sequence:
 *   crd_conv crd_input crd_output [pad_id_override]
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   April 20, 2018 - Initial version
 *   October 16, 2019 - Updated for v2.01 ("NA", etc.).
 *
**-----------------------------------------------------------------------*/

int
main (argc, argv)
     int argc;
     char *argv[];
{
  char str[512];
  FILE *str_in, *str_out;
  int pad_override= 0;

  // Is the input OK?
  if (argc == 3)
    {
      printf("in: %s\nout: %s\n", argv[1], argv[2]);
    }
  else if (argc == 4)
    {
      printf("in: %s\nout: %s\npad_id_override: %s\n", 
             argv[1], argv[2], argv[3]);
      pad_override= atoi (argv[3]);
    }
  else
    {
      printf("Too few or too many arguments: %d\n", argc);
      printf("Usage: input_CRD output_CRD [pad_id_override]\n");
      exit (-1);
    }

  // Try to topen the files
  if ((str_in = fopen (argv[1], "r")) == NULL)
    {
      printf ("Could not open file %s\n", argv[1]);
      exit (1);
    }
  if ((str_out = fopen (argv[2], "w")) == NULL)
    {
      printf ("Could not open file %s\n", argv[2]);
      exit (1);
    }

  /*header.format_version Copy and reformat data */
  while (fgets (str, 512, str_in) != NULL)
    {
      if (strncmp (str, "H1", 2) == 0 ||
          strncmp (str, "h1", 2) == 0)
        {
          read_h1 (str, &h1);

          // convert to the other version!
          if (h1.format_version == 1)
            h1.format_version= 2;
          else if (h1.format_version == 2)
            h1.format_version= 1;

	  write_h1 (str_out, h1);
        }
      else if (strncmp (str, "H2", 2) == 0 ||
               strncmp (str, "h2", 2) == 0)
        {
          read_h2 (str, &h2);
          if (pad_override != 0)
            h2.cdp_pad_id= pad_override;
	  write_h2 (str_out, h2);
        }
      else if (strncmp (str, "H3", 2) == 0 ||
               strncmp (str, "h3", 2) == 0)
        {
          read_h3 (str, &h3);
          if (h1.format_version == 2 && h3.target_loc < 0)
                  h3.target_loc= NA_VALUE;
	  write_h3 (str_out, h3);
        }
      else if (strncmp (str, "H4", 2) == 0 ||
               strncmp (str, "h4", 2) == 0)
        {
          read_h4 (str, &h4);
          if (h1.format_version == 2)
            {
              if (h4.end_year < 0) h4.end_year= NA_VALUE;
              if (h4.end_mon < 0) h4.end_mon= NA_VALUE;
              if (h4.end_day < 0) h4.end_day= NA_VALUE;
              if (h4.end_hour < 0) h4.end_hour= NA_VALUE;
              if (h4.end_min < 0) h4.end_min= NA_VALUE;
              if (h4.end_sec < 0) h4.end_sec= NA_VALUE;
            }
	  write_h4 (str_out, h4);
        }
      else if (strncmp (str, "H5", 2) == 0 ||
               strncmp (str, "h5", 2) == 0)
        {
          read_h5 (str, &h5);
	  write_h5 (str_out, h5);
        }
      else if (strncmp (str, "H8", 2) == 0 ||
               strncmp (str, "h8", 2) == 0)
        {
          read_h8 (str);
	  write_h8 (str_out);
        }
      else if (strncmp (str, "H9", 2) == 0 ||
               strncmp (str, "h9", 2) == 0)
        {
          read_h9 (str);
	  write_h9 (str_out);
        }
      else if (strncmp (str, "C0", 2) == 0 ||
               strncmp (str, "c0", 2) == 0)
        {
          read_c0 (str, &c0);
          write_c0 (str_out, c0);
        }
      else if (strncmp (str, "C1", 2) == 0 ||
               strncmp (str, "c1", 2) == 0)
        {
          read_c1 (str, &c1);
          if (h1.format_version == 2)
            {
              if (c1.nom_fire_rate < 0) c1.nom_fire_rate= NA_VALUEF;
              if (c1.pulse_energy < 0) c1.pulse_energy= NA_VALUEF;
              if (c1.pulse_width < 0) c1.pulse_width= NA_VALUEF;
              if (c1.beam_div < 0) c1.beam_div= NA_VALUEF;
              if (c1.pulses_in_semitrain < 0) c1.pulses_in_semitrain= NA_VALUE;
            }
          write_c1 (str_out, c1);
        }
      else if (strncmp (str, "C2", 2) == 0 ||
               strncmp (str, "c2", 2) == 0)
        {
          read_c2 (str, &c2);
          if (h1.format_version == 2)
            {
             // printf ("C2a, %lf %lf %lf %lf %lf %lf  %lf %lf %lf  \n", c2.qe, c2.voltage, c2.dark_count, c2.output_pulse_width, c2.spectral_filter, c2.spectral_filter_xmission, c2.spatial_filter, c2.amp_gain, c2.amp_bandwidth);
              if (c2.qe < 0) c2.qe= NA_VALUEF;
              if (c2.voltage < 0) c2.voltage= NA_VALUEF;
              if (c2.dark_count < 0) c2.dark_count= NA_VALUEF;
              if (c2.output_pulse_width < 0) c2.output_pulse_width= NA_VALUEF;
              if (c2.spectral_filter < 0) c2.spectral_filter= NA_VALUEF;
              if (c2.spectral_filter_xmission < 0)
                  c2.spectral_filter_xmission= NA_VALUEF;
              if (c2.spatial_filter < 0) c2.spatial_filter= NA_VALUEF;
              if (c2.amp_gain < 0) c2.amp_gain= NA_VALUEF;
              if (c2.amp_bandwidth < 0) c2.amp_bandwidth= NA_VALUEF;
              if (c2.amp_in_use < 0) c2.amp_in_use= NA_VALUE;
              //printf ("C2b, %lf %lf %lf %lf %lf %lf  %lf %lf %lf  \n", c2.qe, c2.voltage, c2.dark_count, c2.output_pulse_width, c2.spectral_filter, c2.spectral_filter_xmission, c2.spatial_filter, c2.amp_gain, c2.amp_bandwidth);
            }
          write_c2 (str_out, c2);
        }
      else if (strncmp (str, "C3", 2) == 0 ||
               strncmp (str, "c3", 2) == 0)
        {
          read_c3 (str, &c3);
          if (h1.format_version == 2 && c3.epoch_delay_corr == -1)
                  c3.epoch_delay_corr= NA_VALUEF;
          write_c3 (str_out, c3);
        }
      else if (strncmp (str, "C4", 2) == 0 ||
               strncmp (str, "c4", 2) == 0)
        {
          read_c4 (str, &c4);
          if (h1.format_version == 2)
            {
              if (c4.est_stn_utc_offset < 0) c4.est_stn_utc_offset= NA_VALUEF;
              if (c4.est_stn_osc_drift < 0) c4.est_stn_osc_drift= NA_VALUEF;
              if (c4.est_xponder_utc_offset < 0)
                  c4.est_xponder_utc_offset= NA_VALUEF;
              if (c4.est_xponder_osc_drift < 0)
                  c4.est_xponder_osc_drift= NA_VALUEF;
              if (c4.xponder_clock_ref_time < 0)
                  c4.xponder_clock_ref_time= NA_VALUEF;
            }
          write_c4 (str_out, c4);
        }
      else if (strncmp (str, "C5", 2) == 0 ||
               strncmp (str, "c5", 2) == 0)
        {
          // Does not exist in V1
          read_c5 (str, &c5);
	  if (h1.format_version == 2)
            write_c5 (str_out, c5);
        }
      else if (strncmp (str, "C6", 2) == 0 ||
               strncmp (str, "c6", 2) == 0)
        {
          // Does not exist in V1
          read_c6 (str, &c6);
	  if (h1.format_version == 2)
            write_c6 (str_out, c6);
        }
      else if (strncmp (str, "C7", 2) == 0 ||
               strncmp (str, "c7", 2) == 0)
        {
          // Does not exist in V1
          read_c7 (str, &c7);
	  if (h1.format_version == 2)
            write_c7 (str_out, c7);
        }
      else if (strncmp (str, "10", 2) == 0)
        {
          read_10 (str, &d10);
          if (h1.format_version == 2)
            {
              d10.xmt_amp= NA_VALUE;
              if (d10.time_of_flight < 0) d10.time_of_flight= NA_VALUEF;
              if (d10.xcv_amp < 0) d10.xcv_amp= NA_VALUE;
              if (d10.xmt_amp < 0) d10.xmt_amp= NA_VALUE;
            }
          write_10 (str_out, d10);
        }
      else if (strncmp (str, "11", 2) == 0)
        {
          read_11 (str, &d11);
          if (h1.format_version == 2)
            {
              d11.signal_to_noise= NA_VALUEF;
              if (d11.time_of_flight < 0) d11.time_of_flight= NA_VALUEF;
              if (d11.bin_skew == -1) d11.bin_skew= NA_VALUEF;
              if (d11.bin_kurtosis == -1) d11.bin_kurtosis= NA_VALUEF;
              if (d11.bin_PmM == -1) d11.bin_PmM= NA_VALUEF;
              if (d11.return_rate == -1) d11.return_rate= NA_VALUEF;
            }
          write_11 (str_out, d11);
        }
      else if (strncmp (str, "12", 2) == 0)
        {
          read_12 (str, &d12);
          if (h1.format_version == 2)
            {
              d12.range_rate= NA_VALUEF;
              if (d12.refraction_corr < 0) d12.refraction_corr= NA_VALUEF;
              if (d12.target_CofM_corr < 0) d12.target_CofM_corr= NA_VALUEF;
              if (d12.nd_value < 0) d12.nd_value= NA_VALUEF;
              if (d12.time_bias < 0) d12.time_bias= NA_VALUEF;
            }
          write_12 (str_out, d12);
        }
      else if (strncmp (str, "20", 2) == 0)
        {
          read_20 (str, &d20);
          if (h1.format_version == 2 && d20.value_origin < 0)
                  d20.value_origin= NA_VALUE;
          write_20 (str_out, d20);
        }
      else if (strncmp (str, "21", 2) == 0)
        {
          read_21 (str, &d21);
          if (h1.format_version == 2)
            {
              d21.sky_temperature= NA_VALUEF;
              if (d21.wind_speed < 0) d21.wind_speed= NA_VALUEF;
              if (d21.wind_direction < 0) d21.wind_direction= NA_VALUEF;
              if (d21.visibility < 0) d21.visibility= NA_VALUE;
              if (d21.sky_clarity < 0) d21.sky_clarity= NA_VALUEF;
              if (d21.atmospheric_seeing < 0) d21.atmospheric_seeing= NA_VALUE;
              if (d21.cloud_cover < 0) d21.cloud_cover= NA_VALUE;
              if (d21.sky_temperature < 0) d21.sky_temperature= NA_VALUEF;
            }
          write_21 (str_out, d21);
        }
      else if (strncmp (str, "30", 2) == 0)
        {
          read_30 (str, &d30);
          if (h1.format_version == 2)
            {
              d30.azimuth_rate= NA_VALUEF;
              d30.elevation_rate= NA_VALUEF;
              if (d30.azimuth == -1) d30.azimuth= NA_VALUEF;
              if (d30.elevation == -1) d30.elevation= NA_VALUEF;
              if (d30.direction_ind < 0) d30.direction_ind= NA_VALUE;
              if (d30.azimuth_rate == -1) d30.azimuth_rate= NA_VALUEF;
              if (d30.elevation_rate == -1) d30.elevation_rate= NA_VALUEF;
            }
          write_30 (str_out, d30);
        }
      else if (strncmp (str, "40", 2) == 0)
        {
          read_40 (str, &d40);
          if (h1.format_version == 2)
            {
              d40.cal_span= NA_VALUE;
              d40.cal_return_rate= NA_VALUEF;
              //printf ("40a %d %d %f %f %f %f %f %f\n", d40.num_points_recorded, d40.num_points_used, d40.one_way_target_dist, d40.cal_rms, d40.cal_skew, d40.cal_kurtosis, d40.cal_PmM, d40.cal_return_rate);
              if (d40.num_points_recorded == -1)
                      d40.num_points_recorded= NA_VALUE;
              if (d40.num_points_used == -1)
                      d40.num_points_used= NA_VALUE;
              if (d40.one_way_target_dist == -1)
                      d40.one_way_target_dist= NA_VALUEF;
              if (d40.cal_delay_shift == -1)
                      d40.cal_delay_shift= NA_VALUEF;
              if (d40.cal_rms == -1)
                      d40.cal_rms= NA_VALUEF;
              if (d40.cal_skew == -1)
                      d40.cal_skew= NA_VALUEF;
              if (d40.cal_kurtosis == -1)
                      d40.cal_kurtosis= NA_VALUEF;
              if (d40.cal_PmM == -1)
                      d40.cal_PmM= NA_VALUEF;
              if (d40.cal_return_rate == -1)
                      d40.cal_return_rate= NA_VALUEF;
              //printf ("40b %d %d %f %f %f %f %f %f\n", d40.num_points_recorded, d40.num_points_used, d40.one_way_target_dist, d40.cal_rms, d40.cal_skew, d40.cal_kurtosis, d40.cal_PmM, d40.cal_return_rate);
            }
          write_40 (str_out, d40);
        }
      else if (strncmp (str, "41", 2) == 0)
        {
          // Does not exist in V1
          read_41 (str, &d41);
	  if (h1.format_version == 2)
            write_41 (str_out, d41);
        }
      else if (strncmp (str, "42", 2) == 0)
        {
          // Does not exist in V1
          read_42 (str, &d42);
	  if (h1.format_version == 2)
            write_42 (str_out, d42);
        }
      else if (strncmp (str, "50", 2) == 0)
        {
          read_50 (str, &d50);
          if (h1.format_version == 2)
            {
              if (d50.sess_rms == -1)
                  d50.sess_rms= NA_VALUEF;
              if (d50.sess_skew == -1)
                  d50.sess_skew= NA_VALUEF;
              if (d50.sess_kurtosis == -1)
                  d50.sess_kurtosis= NA_VALUEF;
              if (d50.sess_PmM == -1)
                  d50.sess_PmM= NA_VALUEF;
            }
          write_50 (str_out, d50);
        }
      else if (strncmp (str, "60", 2) == 0)
        {
          read_60 (str, &d60);
          write_60 (str_out, d60);
        }
      else if (strncmp (str, "00", 2) == 0)
        {
          read_00 (str, &d00);
          write_00 (str_out, d00);
        }

    }
}
