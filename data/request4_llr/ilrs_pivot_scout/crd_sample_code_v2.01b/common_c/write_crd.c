#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "../include/crd.h"

/*-------------------------------------------------------------------------
 * Subroutines: write CRD data records to an output file.
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   July 06, 2007 - Initial version
 *   05/07/08   - Expand all configuration and data section character strings
 *                to allow up to 40 characters (plus NULL). 
 *              - Added detector channel to normalpoint (11) and calibration
 *                (40) records. 
 *              - Added field for 'crd' literal to H1.
 *              - Record '21' sky_clarity is now double rather than int.
 *   06/24/08   - Record '11', np window length is now double rather than int.
 *                (v1.00 rlr)
 *   03/10/09   - Record H2 Epoch Timescale corrected from I1 to I2.
 *                (v1.00a rlr).
 *   03/11/09   - Record H3 changed to print leading zeros rather than spaces
 *                for ilrs_id. (v1.00a rlr)
 *   03/12/18   - Changes for CRD v2.00. rlr.
 *   08/12/18   - Changes for v02 n/a. rlr.
 *
**-----------------------------------------------------------------------*/

char stro[256];
int iisna (int);
int disna (double);
int ldisna (long double);
static int nstat;
static int format_version;
static void remove_blanks (char *, int);
static int writei (FILE *str, char *fmt, int iv);
static int writed (FILE *str, char *fmt, double dv);
static int writeld (FILE *str, char *fmt, long double ldv);

/* Ranging data header/footer records */
  /* H1 - format header */
int
write_h1 (FILE * str_out, struct rh1 header)
{
  format_version= header.format_version;
  if (format_version == 1)
    {
      sprintf (stro,
	       "h1     %2d %4d %2d %2d %2d",
	       header.format_version, header.prod_year, header.prod_mon,
	       header.prod_day, header.prod_hour);
      strncpy (&stro[3], header.crd_literal, 3);
      fprintf (str_out, "%s\n", stro);
      return (0);
    }
  else if (format_version == 2)
    {
      nstat= 0;
      fprintf (str_out,"h1 ");
      fprintf (str_out,"%-s "  , header.crd_literal);
      if (writei (str_out, "%d ", header.format_version) < 0) nstat--;
      if (writei (str_out, "%d ", header.prod_year) < 0) nstat--;
      if (writei (str_out, "%d ", header.prod_mon) < 0) nstat--;
      if (writei (str_out, "%d ", header.prod_day) < 0) nstat--;
      if (writei (str_out, "%d", header.prod_hour) < 0) nstat--;
      fprintf (str_out, "\n");
      return (nstat);
    }
}

  /* H2 - station header */
int
write_h2 (FILE * str_out, struct rh2 header)
{
  char temp[11];

  if (format_version == 1)
    {
      sprintf (stro,
           "h2            %4d %2d %2d %2d",
           header.cdp_pad_id, header.cdp_sys_num, header.cdp_occ_num, 
           header.stn_timescale);
      strncpy (&stro[3], header.stn_name, strlen(header.stn_name));
      fprintf (str_out, "%s\n", stro);
      return (0);
    }
  else if (format_version == 2)
    {
      strncpy (temp, header.stn_name, 11);
      remove_blanks (temp, strlen(temp));

      nstat= 0;
      fprintf (str_out,"h2 ");
      fprintf (str_out,"%-s "  , temp);
      if (writei (str_out, "%04d ", header.cdp_pad_id) < 0) nstat--;
      if (writei (str_out, "%02d ", header.cdp_sys_num) < 0) nstat--;
      if (writei (str_out, "%02d ", header.cdp_occ_num) < 0) nstat--;
      if (writei (str_out, "%d ", header.stn_timescale) < 0) nstat--;
      fprintf (str_out,"%-s"  , header.stn_network);
      fprintf (str_out, "\n");
      //fprintf (str_out, "%s\n", stro);
      return (nstat);
    }
}

  /* H3 - spacecraft header */
int
write_h3 (FILE * str_out, struct rh3 header)
{
  char temp[11];

  // Does code need to convert between target_type and _class/_loc?
  if (format_version == 1)
    {
      sprintf (stro,
           "h3             %07d %4d %8d %1d %1d", header.ilrs_id, header.sic, 
	   header.norad, header.SC_timescale, header.target_type);
      strncpy (&stro[3], header.target_name, strlen(header.target_name));
      fprintf (str_out, "%s\n", stro);
      return (0);
    }
  else if (format_version == 2)
    {
      strncpy (temp, header.target_name, 11);
      remove_blanks (temp, strlen(temp));

      nstat= 0;
      fprintf (str_out,"h3 ");
      fprintf (str_out,"%-s " , temp);
      if (writei (str_out, "%07d ", header.ilrs_id) < 0) nstat--;
      if (writei (str_out, "%04d ", header.sic) < 0) nstat--;
      if (writei (str_out, "%05d ", header.norad) < 0) nstat--;
      if (writei (str_out, "%d ", header.SC_timescale) < 0) nstat--;
      if (writei (str_out, "%d ", header.target_class) < 0) nstat--;
      if (writei (str_out, "%d", header.target_loc) < 0) nstat--;
      fprintf (str_out, "\n");
      return (nstat);
    }
}

  /* H4 - Session header */
int
write_h4 (FILE * str_out, struct rh4 header)
{
  if (format_version == 1)
    {
      fprintf (str_out,
           "h4 %2d %4d %2d %2d %2d %2d %2d %4d %2d %2d %2d %2d %2d %2d %1d %1d %1d %1d %1d %1d %1d\n",
           header.data_type, 
	   header.start_year, header.start_mon, header.start_day, 
	   header.start_hour, header.start_min, header.start_sec,
           header.end_year, header.end_mon, header.end_day, 
	   header.end_hour, header.end_min, header.end_sec, 
	   header.data_release, header.refraction_app_ind, 
	   header.CofM_app_ind, header.xcv_amp_app_ind,
           header.stn_sysdelay_app_ind, 
	   header.SC_sysdelay_app_ind, 
	   header.range_type_ind, header.data_qual_alert_ind);
      return (0);
    }
  else if (format_version == 2)
    {
      nstat= 0;
      fprintf (str_out,"h4 ");
      if (writei (str_out, "%d ", header.data_type) < 0) nstat--;
      if (writei (str_out, "%d ", header.start_year) < 0) nstat--;
      if (writei (str_out, "%d ", header.start_mon) < 0) nstat--;
      if (writei (str_out, "%d ", header.start_day) < 0) nstat--;
      if (writei (str_out, "%d ", header.start_hour) < 0) nstat--;
      if (writei (str_out, "%d ", header.start_min) < 0) nstat--;
      if (writei (str_out, "%d ", header.start_sec) < 0) nstat--;
      if (writei (str_out, "%d ", header.end_year) < 0) nstat--;
      if (writei (str_out, "%d ", header.end_mon) < 0) nstat--;
      if (writei (str_out, "%d ", header.end_day) < 0) nstat--;
      if (writei (str_out, "%d ", header.end_hour) < 0) nstat--;
      if (writei (str_out, "%d ", header.end_min) < 0) nstat--;
      if (writei (str_out, "%d ", header.end_sec) < 0) nstat--;
      if (writei (str_out, "%d ", header.data_release) < 0) nstat--;
      if (writei (str_out, "%d ", header.refraction_app_ind) < 0) nstat--;
      if (writei (str_out, "%d ", header.CofM_app_ind) < 0) nstat--;
      if (writei (str_out, "%d ", header.xcv_amp_app_ind) < 0) nstat--;
      if (writei (str_out, "%d ", header.stn_sysdelay_app_ind) < 0) nstat--;
      if (writei (str_out, "%d ", header.SC_sysdelay_app_ind) < 0) nstat--;
      if (writei (str_out, "%d ", header.range_type_ind) < 0) nstat--;
      if (writei (str_out, "%d", header.data_qual_alert_ind) < 0) nstat--;
      fprintf (str_out, "\n");
      return (nstat);
    }
}

  /* H5 - prediction header */
  /* New for V2 */
int
write_h5 (FILE * str_out, struct rh5 header)
{
  char temp1[13];
  char temp2[11];

  strncpy (temp1, header.date_and_time, 13);
  remove_blanks (temp1, strlen(temp1));
  strncpy (temp2, header.prediction_provider, 11);
  remove_blanks (temp2, strlen(temp2));

  if (format_version == 1) return (0);	// no "h5" in v1
  nstat= 0;
  fprintf (str_out,"h5 ");
  if (writei (str_out, "%d ", header.prediction_type) < 0) nstat--;
  if (writei (str_out, "%0d ", header.year_of_century) < 0) nstat--;
  fprintf (str_out,"%-s "  , temp1);
  fprintf (str_out,"%-s "  , temp2);
  if (writei (str_out, "%d", header.sequence_number) < 0) nstat--;
  fprintf (str_out, "\n");
  return (nstat);
}

  /* Need indicators that these have been read? */
  /* H8 - End of Session footer */
int
write_h8 (FILE * str_out)
{
  fprintf (str_out, "h8\n");
  return (0);
}

  /* H9 - End of File footer */
int
write_h9 (FILE * str_out)
{
  fprintf (str_out, "h9\n");
  return (0);
}

/* Ranging data configuration records (1 of n) */
    /* C0 - System Configuration Record */
int
write_c0 (FILE * str_out, struct rc0 config)
{
  int i;

  fprintf (str_out, 
	   "c0 %1d %.3f %-s", 
	   config.detail_type, config.xmit_wavelength, config.config_ids[0]);

// This presumes that config_ids that are not used are blank.
  for (i=1; i<10; i++)
    {
      if (config.config_ids[i][0] != '\0') 
	fprintf (str_out, " %-s", config.config_ids[i]);
    }
  fprintf (str_out, "\n");
  return (0);
}

    /* C1 - Laser Configuration Record */
int
write_c1 (FILE * str_out, struct rc1 config)
{
  nstat= 0;
  fprintf (str_out,"c1 ");
  if (writei (str_out, "%d ", config.detail_type) < 0) nstat--;
  fprintf (str_out,"%-s "    , config.laser_config_id);
  fprintf (str_out,"%-s "    , config.laser_type);
  if (writed (str_out, "%.2lf ", config.prim_wavelength) < 0) nstat--;
  if (writed (str_out, "%.2lf ", config.nom_fire_rate) < 0) nstat--;
  if (writed (str_out, "%.2lf ", config.pulse_energy) < 0) nstat--;
  if (writed (str_out, "%.1lf ", config.pulse_width) < 0) nstat--;
  if (writed (str_out, "%.2lf ", config.beam_div) < 0) nstat--;
  if (writei (str_out, "%d", config.pulses_in_semitrain) < 0) nstat--;
  fprintf (str_out,"\n");
  return (nstat);
}

    /* C2 - Detector Configuration Record */
int
write_c2 (FILE * str_out, struct rc2 config)
{
  if (format_version == 1)
    {
      fprintf (str_out, 
               "c2 %d %-s %-s %.3f %.2f %.1f %.1f %-s %.1f %.2f %.1f %.1f %-s\n", 
               config.detail_type, config.detector_config_id, 
	       config.detector_type, config.app_wavelength, config.qe, 
	       config.voltage, config.dark_count,
	       config.output_pulse_type, config.output_pulse_width, 
	       config.spectral_filter, config.spectral_filter_xmission,
	       config.spatial_filter, config.signal_proc);
      return (0);
    }
  else if (format_version == 2)
    {
      nstat= 0;
      fprintf (str_out,"c2 ");
      if (writei (str_out, "%d ", config.detail_type) < 0) nstat--;
      fprintf (str_out,"%-s "    , config.detector_config_id);
      fprintf (str_out,"%-s "    , config.detector_type);
      if (writed (str_out, "%.3lf ", config.app_wavelength) < 0) nstat--;
      if (writed (str_out, "%.2lf ", config.qe) < 0) nstat--;
      if (writed (str_out, "%.1lf ", config.voltage) < 0) nstat--;
      if (writed (str_out, "%.1lf ", config.dark_count) < 0) nstat--;
      fprintf (str_out,"%-s "    , config.output_pulse_type);
      if (writed (str_out, "%.1lf ", config.output_pulse_width) < 0) nstat--;
      if (writed (str_out, "%.2lf ", config.spectral_filter) < 0) nstat--;
      if (writed (str_out, "%.1lf ", config.spectral_filter_xmission) < 0) nstat--;
      if (writed (str_out, "%.2lf ", config.spatial_filter) < 0) nstat--;
      fprintf (str_out,"%-s "    , config.signal_proc);
      if (writed (str_out, "%.1lf ", config.amp_gain) < 0) nstat--;
      if (writed (str_out, "%.1lf ", config.amp_bandwidth) < 0) nstat--;
      if (writei (str_out, "%d", config.amp_in_use) < 0) nstat--;
      fprintf (str_out,"\n");
      return (nstat);
    }
}

    /* C3 - Timing Configuration Record */
int
write_c3 (FILE * str_out, struct rc3 config)
{
/**
  fprintf (str_out, 
           "c3 %d %-s %-s %-s %-s %-s %.1f\n",
           config.detail_type, config.timing_config_id, config.time_source,
	   config.freq_source, config.timer, config.timer_serial_num,
	   config.epoch_delay_corr);
	   **/
  fprintf (str_out, 
           "c3 %d %-s %-s %-s %-s %-s ",
           config.detail_type, config.timing_config_id, config.time_source,
	   config.freq_source, config.timer, config.timer_serial_num);
  if (writed (str_out, "%.1lf", config.epoch_delay_corr) < 0) nstat--;
  fprintf (str_out,"\n");
  return (0);
}

    /* C4 - Transponder Configuration Record */
int
write_c4 (FILE * str_out, struct rc4 config)
{
  nstat= 0;

  // Warning: There may be problems with NA_VALUE being too small for some of
  // these fields.
  fprintf (str_out,"c4 ");
  if (writei (str_out, "%d ", config.detail_type) < 0) nstat--;
  fprintf (str_out,"%-s "    , config.xponder_config_id);
  if (writeld (str_out, "%.3Lf ", config.est_stn_utc_offset) < 0) nstat--;
  if (writed (str_out, "%.2lf ", config.est_stn_osc_drift) < 0) nstat--;
  if (writeld (str_out, "%.3Lf ", config.est_xponder_utc_offset) < 0) nstat--;
  if (writed (str_out, "%.2lf ", config.est_xponder_osc_drift) < 0) nstat--;
  if (writeld (str_out, "%.2Lf ", config.xponder_clock_ref_time) < 0) nstat--;
  if (writei (str_out, "%d ", config.stn_off_drift_app_ind) < 0) nstat--;
  if (writei (str_out, "%d ", config.SC_off_drift_app_ind) < 0) nstat--;
  if (writei (str_out, "%d", config.SC_time_simplified_ind) < 0) nstat--;
  fprintf (str_out,"\n");
  return (nstat);
}

    /* C5 - Software Configuration Record */
int
write_c5 (FILE * str_out, struct rc5 config)
{
  if (format_version == 1) return(0);	// no "c5" in v1
  fprintf (str_out, 
           "c5 %d %-s %-s %-s %-s %-s\n",
           config.detail_type, config.software_config_id, 
           config.tracking_software, config.tracking_software_versions,
           config.processing_software, config.processing_software_versions);
  return (0);
}

    /* C6 - Meteorological Instrument Configuration Record */
int
write_c6 (FILE * str_out, struct rc6 config)
{
  if (format_version == 1) return(0);	// no "c6" in v1
  fprintf (str_out, 
           "c6 %d %-s %-s %-s %-s %-s %-s %-s %-s %-s %-s\n",
           config.detail_type, config.met_config_id, 
           config.pressure_sensor_manufacturer, 
           config.pressure_sensor_model, 
           config.pressure_sensor_serial_num,
           config.temperature_sensor_manufacturer, 
           config.temperature_sensor_model, 
           config.temperature_sensor_serial_num,
           config.humidity_sensor_manufacturer, 
           config.humidity_sensor_model, 
           config.humidity_sensor_serial_num);
  return (0);
}

    /* C7 - Calibration Configuration Record */
int
write_c7 (FILE * str_out, struct rc7 config)
{
  if (format_version == 1) return(0);	// no "c7" in v1
  nstat= 0;
  fprintf (str_out,"c7 ");
  if (writei (str_out, "%d ", config.detail_type) < 0) nstat--;
  fprintf (str_out,"%-s "    , config.calconfig_id);
  fprintf (str_out,"%-s "    , config.target_name);
  if (writed (str_out, "%.5lf ", config.surveyed_target_dist) < 0) nstat--;
  if (writed (str_out, "%.2lf ", config.survey_error) < 0) nstat--;
  if (writed (str_out, "%.5lf ", config.other_fixed_delays) < 0) nstat--;
  if (writed (str_out, "%.2lf ", config.pulse_energy) < 0) nstat--;
  fprintf (str_out,"%-s "    , config.processing_software);
  fprintf (str_out,"%-s"     , config.processing_software_version);
  fprintf (str_out,"\n");
  return (nstat);
}

/* Ranging data records */
    /* 10 - Range Record */
int
write_10 (FILE * str_out, struct rd10 data_recd)
{
  if (format_version == 1)
    {
      fprintf (str_out,
           "10 %.12Lf %.12Lf %-s %d %d %d %d %d\n",
           data_recd.sec_of_day, data_recd.time_of_flight, 
	   data_recd.sysconfig_id, data_recd.epoch_event, 
	   data_recd.filter_flag, data_recd.detector_channel, 
	   data_recd.stop_number, data_recd.xcv_amp);
      return (0);
    }
  else if (format_version == 2)
    {
      nstat= 0;
      fprintf (str_out,"10 ");
      if (writeld (str_out, "%.12Lf ", data_recd.sec_of_day) < 0) nstat--;
      if (writeld (str_out, "%.12Lf ", data_recd.time_of_flight) < 0) nstat--;
      fprintf (str_out,"%-s "  , data_recd.sysconfig_id);
      if (writei (str_out, "%d ", data_recd.epoch_event) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.filter_flag) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.detector_channel) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.stop_number) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.xcv_amp) < 0) nstat--;
      if (writei (str_out, "%d", data_recd.xmt_amp) < 0) nstat--;
      fprintf (str_out,"\n");
      return (nstat);
    }
}

    /* 11 - Normal Point Record */
int
write_11 (FILE * str_out, struct rd11 data_recd)
{
  if (format_version == 1)
    {
      fprintf (str_out,
           "11 %.12Lf %.12Lf %-s %d %.1f %d %.1f %.3f %.3f %.1f %.1f %d\n",
           data_recd.sec_of_day, data_recd.time_of_flight, 
	   data_recd.sysconfig_id, data_recd.epoch_event, 
	   data_recd.np_window_length, data_recd.num_ranges, 
	   data_recd.bin_rms, data_recd.bin_skew, data_recd.bin_kurtosis, 
	   data_recd.bin_PmM, data_recd.return_rate, 
           data_recd.detector_channel);
      return (0);
    }
  else if (format_version == 2)
    {
      nstat= 0;
      fprintf (str_out,"11 ");
      if (writeld (str_out, "%.12Lf ", data_recd.sec_of_day) < 0) nstat--;
      if (writeld (str_out, "%.12Lf ", data_recd.time_of_flight) < 0) nstat--;
      fprintf (str_out,"%-s "   ,   data_recd.sysconfig_id);
      if (writei (str_out, "%d ",   data_recd.epoch_event) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.np_window_length) < 0) nstat--;
      if (writei (str_out, "%d ",   data_recd.num_ranges) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.bin_rms) < 0) nstat--;
      if (writed (str_out, "%.3f ", data_recd.bin_skew) < 0) nstat--;
      if (writed (str_out, "%.3f ", data_recd.bin_kurtosis) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.bin_PmM) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.return_rate) < 0) nstat--;
      if (writei (str_out, "%d ",   data_recd.detector_channel) < 0) nstat--;
      if (writed (str_out, "%.1f", data_recd.signal_to_noise) < 0) nstat--;
      fprintf (str_out,"\n");
      return (nstat);
    }
}

    /* 12 - Range Supplement Record */
int
write_12 (FILE * str_out, struct rd12 data_recd)
{
  if (format_version == 1)
    {
      fprintf (str_out,
           "12 %.7Lf %-s %.1f %.4f %.2f %.4f\n",
           data_recd.sec_of_day, data_recd.sysconfig_id,
	   data_recd.refraction_corr, data_recd.target_CofM_corr, 
	   data_recd.nd_value, data_recd.time_bias);
      return (0);
    }
  else if (format_version == 2)
    {
      nstat= 0;
      fprintf (str_out,"12 ");
      if (writeld (str_out, "%.7Lf ", data_recd.sec_of_day) < 0) nstat--;
      fprintf (str_out,"%-s "   , data_recd.sysconfig_id);
      if (writed (str_out, "%.1f ", data_recd.refraction_corr) < 0) nstat--;
      if (writed (str_out, "%.4f ", data_recd.target_CofM_corr) < 0) nstat--;
      if (writed (str_out, "%.2f ", data_recd.nd_value) < 0) nstat--;
      if (writed (str_out, "%.4f ", data_recd.time_bias) < 0) nstat--;
      if (writed (str_out, "%.16f", data_recd.range_rate) < 0) nstat--;
      fprintf (str_out,"\n");
      return (nstat);
    }
}

    /* 20 - Meteorological Record */
int
write_20 (FILE * str_out, struct rd20 data_recd)
{
  nstat= 0;
  fprintf (str_out,"20 ");
  if (writeld (str_out, "%.3Lf ", data_recd.sec_of_day) < 0) nstat--;
  if (writed (str_out, "%.2f ", data_recd.pressure) < 0) nstat--;
  if (writed (str_out, "%.2f ", data_recd.temperature) < 0) nstat--;
  if (writed (str_out, "%.0f ", data_recd.humidity) < 0) nstat--;
  if (writei (str_out, "%d", data_recd.value_origin) < 0) nstat--;
  fprintf (str_out,"\n");
  return (nstat);
}

    /* 21 - Meteorological Supplement Record */
int
write_21 (FILE * str_out, struct rd21 data_recd)
{
  if (format_version == 1)
    {
      fprintf (str_out,
           "21 %.3Lf %.1f %.1f %-s %d %.2f %d %d\n",
           data_recd.sec_of_day, data_recd.wind_speed,data_recd.wind_direction, 
	   data_recd.weather_conditions, data_recd.visibility, 
           data_recd.sky_clarity, data_recd.atmospheric_seeing, 
           data_recd.cloud_cover);
      return (0);
    }
  else if (format_version == 2)
    {
      nstat= 0;
      fprintf (str_out,"21 ");
      if (writeld (str_out, "%.3Lf ", data_recd.sec_of_day) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.wind_speed) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.wind_direction) < 0) nstat--;
      fprintf (str_out,"%-s "   , data_recd.weather_conditions);
      if (writei (str_out, "%d ", data_recd.visibility) < 0) nstat--;
      if (writed (str_out, "%.2f ", data_recd.sky_clarity) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.atmospheric_seeing) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.cloud_cover) < 0) nstat--;
      if (writed (str_out, "%.2f", data_recd.sky_temperature) < 0) nstat--;
      fprintf (str_out,"\n");
      return (nstat);
    }
}

    /* 30 - Pointing Angles Record */
int
write_30 (FILE * str_out, struct rd30 data_recd)
{
  if (format_version == 1)
    {
      fprintf (str_out,
           "30 %.3Lf %.4f %.4f %d %d %d\n",
           data_recd.sec_of_day, data_recd.azimuth, data_recd.elevation, 
	   data_recd.direction_ind, data_recd.angle_origin_ind, 
	   data_recd.refraction_corr_ind);
      return (0);
    }
  else if (format_version == 2)
    {
      nstat= 0;
      fprintf (str_out,"30 ");
      if (writeld (str_out, "%.3Lf ", data_recd.sec_of_day) < 0) nstat--;
      if (writed (str_out, "%.4f ", data_recd.azimuth) < 0) nstat--;
      if (writed (str_out, "%.4f ", data_recd.elevation) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.direction_ind) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.angle_origin_ind) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.refraction_corr_ind) < 0) nstat--;
      if (writed (str_out, "%.7f ", data_recd.azimuth_rate) < 0) nstat--;
      if (writed (str_out, "%.7f", data_recd.elevation_rate) < 0) nstat--;
      fprintf (str_out,"\n");
      return (nstat);
    }
}

    /* 40 - Calibration Record */
int
write_40 (FILE * str_out, struct rd40 data_recd)
{
/* Due to some problem w/ fedora 8 or crd_cal.... */
  if (format_version == 1)
    {
      fprintf (str_out,
               "40 %.7Lf %1d %-s %d %d %.3f",
               data_recd.sec_of_day, data_recd.type_of_data, 
	       data_recd.sysconfig_id, data_recd.num_points_recorded, 
	       data_recd.num_points_used, data_recd.one_way_target_dist);
      fprintf (str_out,
               " %.1f %.1f %.1f %.3f %.3f %.1f %1d %1d %1d\n",
               data_recd.cal_sys_delay, data_recd.cal_delay_shift, 
	       data_recd.cal_rms, data_recd.cal_skew, data_recd.cal_kurtosis, 
	       data_recd.cal_PmM, data_recd.cal_type_ind,
               data_recd.cal_shift_type_ind, data_recd.detector_channel);
      return (0);
    }
  else if (format_version == 2)
    {
      nstat= 0;
      fprintf (str_out,"40 ");
      if (writeld (str_out, "%.12Lf ", data_recd.sec_of_day) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.type_of_data) < 0) nstat--;
      fprintf (str_out,"%-s "   , data_recd.sysconfig_id);
      if (writei (str_out, "%d ", data_recd.num_points_recorded) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.num_points_used) < 0) nstat--;
      if (writed (str_out, "%.3f ", data_recd.one_way_target_dist) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.cal_sys_delay) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.cal_delay_shift) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.cal_rms) < 0) nstat--;
      if (writed (str_out, "%.3f ", data_recd.cal_skew) < 0) nstat--;
      if (writed (str_out, "%.3f ", data_recd.cal_kurtosis) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.cal_PmM) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.cal_type_ind) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.cal_shift_type_ind) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.detector_channel) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.cal_span) < 0) nstat--;
      if (writed (str_out, "%.2f", data_recd.cal_return_rate) < 0) nstat--;
      fprintf (str_out,"\n");
      return (nstat);
    }
}

    /* 41 - Calibration Detail Record */
int
write_41 (FILE * str_out, struct rd40 data_recd)
{
/* Two prints is due to some problem w/ fedora 8 or crd_cal.... */
  if (format_version == 1) return (0);	// No "41" in v1
      nstat= 0;
      fprintf (str_out,"41 ");
      if (writeld (str_out, "%.12Lf ", data_recd.sec_of_day) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.type_of_data) < 0) nstat--;
      fprintf (str_out,"%-s "   , data_recd.sysconfig_id);
      if (writei (str_out, "%d ", data_recd.num_points_recorded) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.num_points_used) < 0) nstat--;
      if (writed (str_out, "%.3f ", data_recd.one_way_target_dist) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.cal_sys_delay) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.cal_delay_shift) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.cal_rms) < 0) nstat--;
      if (writed (str_out, "%.3f ", data_recd.cal_skew) < 0) nstat--;
      if (writed (str_out, "%.3f ", data_recd.cal_kurtosis) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.cal_PmM) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.cal_type_ind) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.cal_shift_type_ind) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.detector_channel) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.cal_span) < 0) nstat--;
      if (writed (str_out, "%.2f", data_recd.cal_return_rate) < 0) nstat--;
      fprintf (str_out,"\n");
      return (nstat);
}

    /* 42 - Calibration "Shot" Record */
int
write_42 (FILE * str_out, struct rd42 data_recd)
{
  if (format_version == 1) return (0);	// No "42" in v1
      nstat= 0;
      fprintf (str_out,"42 ");
      if (writeld (str_out, "%.12Lf ", data_recd.sec_of_day) < 0) nstat--;
      if (writeld (str_out, "%.12Lf ", data_recd.time_of_flight) < 0) nstat--;
      fprintf (str_out,"%-s "   , data_recd.sysconfig_id);
      fprintf (str_out,"%-s "   , data_recd.calconfig_id);
      if (writed (str_out, "%.5f ", data_recd.other_variable_delays) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.type_of_data) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.cal_type_ind) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.filter_flag) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.detector_channel) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.stop_number) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.cal_span) < 0) nstat--;
      if (writei (str_out, "%d ", data_recd.xcv_amp) < 0) nstat--;
      if (writei (str_out, "%d", data_recd.xmt_amp) < 0) nstat--;
      fprintf (str_out,"\n");
      return (nstat);
}

    /* 50 - Session Statistics Record */
int
write_50 (FILE * str_out, struct rd50 data_recd)
{
      nstat= 0;
      fprintf (str_out,"50 ");
      fprintf (str_out,"%-s "   , data_recd.sysconfig_id);
      if (writed (str_out, "%.1f ", data_recd.sess_rms) < 0) nstat--;
      if (writed (str_out, "%.3f ", data_recd.sess_skew) < 0) nstat--;
      if (writed (str_out, "%.3f ", data_recd.sess_kurtosis) < 0) nstat--;
      if (writed (str_out, "%.1f ", data_recd.sess_PmM) < 0) nstat--;
      if (writei (str_out, "%d", data_recd.data_qual_ind) < 0) nstat--;
      fprintf (str_out,"\n");
      return (nstat);
}

    /* 60 - Compatibility Record */
int
write_60 (FILE * str_out, struct rd60 data_recd)
{
  // Not used in version 2.
  // Enable in version 2 ONLY if converting pre-CRD files.
  if (format_version == 1)
    {
      fprintf (str_out, 
	 "60 %-s %d %d\n", 
	 data_recd.sysconfig_id, data_recd.sys_change_ind, 
         data_recd.sys_config_ind); 
      return (0);
    }
}

    /* 9X - User Defined Record */
int
write_9x (FILE * str_out, struct rd9x data_recd)
{
}

    /* 00 - Comment Record */
int
write_00 (FILE * str_out, struct rd00 data_recd)
{
  fprintf (str_out,
	   "00 %-s\n", 
	   data_recd.comment);
  return (0);
}

void
remove_blanks (char *str, int slen)
{
  int i, j, k;

// Handle left- and right-justified fields.
// Assumes no blanks mid-field.
  if (isalnum(str[0]))
    {
      // left-justified
      for (i=slen; i>0; i--)	// be sure to have at least one ' '
        { 
          if (str[i] == ' ') str[i]= '\0';
        }
    }
  else
    {
      // right-justified
      for (i=0; i<slen; i++)
        {
          if (str[i] != ' ') break;
        }
      k= 0;
      for (j=i; j<slen; j++)
        {
          str[k++]= str[j];
        }
      str[k]= '\0';
    }
}

// Write an integer number or "na"
int
writei (FILE *str, char *fmt, int iv)
{
  //if (iv > NA_VALUE)
  if (!iisna (iv))
    { 
      fprintf (str, fmt, iv);
    }    
  else
    {
      if (fmt[strlen(fmt)-1] == ' ')
	fprintf (str,"%-s "   , "na");
      else
        fprintf (str,"%-s"   , "na");
    }
}

// Write a double precisions number or "na"
int
writed (FILE *str, char *fmt, double dv)
{
  //if (dv > NA_VALUE)
  if (!disna (dv))
    { 
      fprintf (str, fmt, dv);
    }    
  else
    {
      if (fmt[strlen(fmt)-1] == ' ')
	fprintf (str,"%-s "   , "na");
      else
        fprintf (str,"%-s"   , "na");
    }
}

// Write a long double precisions number or "na"
int
writeld (FILE *str, char *fmt, long double ldv)
{
  //if (ldv > NA_VALUE)
  if (!ldisna (ldv))
    { 
      fprintf (str, fmt, ldv);
    }    
  else
    {
      if (fmt[strlen(fmt)-1] == ' ')
	fprintf (str,"%-s "   , "na");
      else
        fprintf (str,"%-s"   , "na");
    }
}

// Is the field "na"?
int
iisna (int iv)
{
  if (abs (iv- NA_VALUE) < 1) return (1);
  return (0);
}

// Is the field "na"?
int
disna (double dv)
{
  if (fabs (dv- NA_VALUEF) < 1.0e-3) return (1);
  return (0);
}

// Is the field "na"?
int
ldisna (long double ldv)
{
  if (fabsl (ldv- NA_VALUEF) < 1.0e-3) return (1);
  return (0);
}
