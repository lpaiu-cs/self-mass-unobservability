C**-------------------------------------------------------------------------
C * Subroutines: write CRD data records to an output file
C *
C * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
C *
C * History:
C *   July 06, 2007 - Initial version
C  05/07/08   - Expand configuration and data record character fields to
C               allow up to 40 characters.
C             - Added detector channel to normalpoint (11) and calibration (40)
C               records.
C             - Added field for 'crd' literal to 'h1'.
C             - Record '21' sky_clarity is not double rather than int.
C  06/24/08   - Record '11' np window length is now double rather than
C               int. (v1.00 rlr)
C  03/10/09   - Record H2 Epoch Timescale corrected from I1 to I2.
C               (v1.00a rlr).
C  03/10/09   - Record H3 changed to print leading zeros rather than
C               spaces for ilrs_id. (v1.00a rlr).
C  03/20/18   - Changes for CRD v2.00. rlr.
C  06/26/19   - Add calibartion record c7 and 42 to v2.00. rlr.
C  08/08/19   - Added n/a to output. rlr.
C  10/08/20   - write_60: changed location of 'str=""' so that for a CRD
C               v2 pass, the return string would have length=0. v 2.01a. rlr.
C *
C**-------------------------------------------------------------------------

C Ranging data header/footer records
C H1 - format header
      SUBROUTINE write_h1 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

      write (str,"(a)") "h1"
      write (str,"(A)") trim(str)//" "//trim(crd_literal);

      if (writei (str, "(i2)", format_version) .lt. 0) go to 100
      if (writei (str, "(i4)", prod_year) .lt. 0) go to 100
      if (writei (str, "(i2)", prod_mon) .lt. 0) go to 100
      if (writei (str, "(i2)", prod_day) .lt. 0) go to 100
      if (writei (str, "(i2)", prod_hour) .lt. 0) go to 100
      RETURN

 100  WRITE(*,*) "Error writing CRD record type h1"
      STOP 1

      END

C H2 - station header
      SUBROUTINE write_h2 (str)
      IMPLICIT none
CC      INTEGER trimlen
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        WRITE (str,1000,err=100) stn_name,
     &           cdp_pad_id, cdp_sys_num, cdp_occ_num, stn_timescale
      ELSE IF (format_version .EQ. 2) THEN
        write (str,"(a)") "h2"
        write (str,"(A)") trim(str)//" "//trim(stn_name);

        if (writei (str, "(i4)", cdp_pad_id) .lt. 0) go to 100
        if (writei (str, "(i2)", cdp_sys_num) .lt. 0) go to 100
        if (writei (str, "(i2)", cdp_occ_num) .lt. 0) go to 100
        if (writei (str, "(i2)", stn_timescale) .lt. 0) go to 100
        write (str,"(A)") trim(str)//" "//trim(stn_network);
      ENDIF

      RETURN

 100  WRITE(*,*) "Error writing CRD record type h2"
      STOP 1

 1000 FORMAT ("h2",1x,a10,1x,i4,1x,i2,1x,i2,1x,i2)
      END

C H3 - spacecraft header
      SUBROUTINE write_h3 (str)
      IMPLICIT none
      INTEGER trimlen
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        WRITE (str,1000,err=100) target_name,
     &           ilrs_id, sic, norad, SC_timescale, target_type
      ELSE IF (format_version .EQ. 2) THEN
        write (str,"(a)") "h3"
        write (str,"(A)") trim(str)//" "//trim(target_name);

        if (writei (str, "(i7.7)", ilrs_id) .lt. 0) go to 100
        if (writei (str, "(i4.4)", sic) .lt. 0) go to 100
        if (writei (str, "(i5.5)", norad) .lt. 0) go to 100
        if (writei (str, "(i1)", SC_timescale) .lt. 0) go to 100
        if (writei (str, "(i1)", target_class) .lt. 0) go to 100
        if (writei (str, "(i1)", target_loc) .lt. 0) go to 100
      ENDIF
      RETURN

 100  WRITE(*,*) "Error writing CRD record type h3"
      STOP 1

 1000 FORMAT ("h3",1x,a10,1x,i8.7,1x,i4,1x,i8,1x,i1,1x,i1)
      END

C H4 - Session header
      SUBROUTINE write_h4 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

      write (str,"(a)") "h4"
      if (writei (str, "(i2)", data_type) .lt. 0) go to 100
      if (writei (str, "(i4)", start_year) .lt. 0) go to 100
      if (writei (str, "(i2)", start_mon) .lt. 0) go to 100
      if (writei (str, "(i2)", start_day) .lt. 0) go to 100
      if (writei (str, "(i2)", start_hour) .lt. 0) go to 100
      if (writei (str, "(i2)", start_min) .lt. 0) go to 100
      if (writei (str, "(i2)", start_sec) .lt. 0) go to 100
      if (writei (str, "(i4)", end_year) .lt. 0) go to 100
      if (writei (str, "(i2)", end_mon) .lt. 0) go to 100
      if (writei (str, "(i2)", end_day) .lt. 0) go to 100
      if (writei (str, "(i2)", end_hour) .lt. 0) go to 100
      if (writei (str, "(i2)", end_min) .lt. 0) go to 100
      if (writei (str, "(i2)", end_sec) .lt. 0) go to 100
      if (writei (str, "(i2)", data_release) .lt. 0) go to 100
      if (writei (str, "(i1)", refraction_app_ind) .lt. 0) go to 100
      if (writei (str, "(i1)", CofM_app_ind) .lt. 0) go to 100
      if (writei (str, "(i1)", xcv_amp_app_ind) .lt. 0) go to 100
      if (writei (str, "(i1)", stn_sysdelay_app_ind) .lt. 0) go to 100
      if (writei (str, "(i1)", SC_sysdelay_app_ind) .lt. 0) go to 100
      if (writei (str, "(i1)", range_type_ind) .lt. 0) go to 100
      if (writei (str, "(i1)", data_qual_alert_ind) .lt. 0) go to 100

      RETURN

 100  WRITE(*,*) "Error writing CRD record type h4"
      STOP 1

      END

C H5 - Prediction header
C New in v2
      SUBROUTINE write_h5 (str)
      IMPLICIT none
      INTEGER trimlen
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

C     No "h5" in version 1
      str= ""
      IF (format_version .EQ. 1) RETURN
      write (str,"(a)") "h5"
      if (writei (str, "(i1)", prediction_type) .lt. 0) go to 100
      if (writei (str, "(i2)", year_of_century) .lt. 0) go to 100
      write (str,"(A)",err=100) trim(str)//" "//trim(date_time)
      write (str,"(A)",err=100) 
     &       trim(str)//" "//trim(prediction_provider)
      if (writei (str, "(i5)", sequence_number) .lt. 0) go to 100
      RETURN

 100  WRITE(*,*) "Error writing CRD record type h5"
      STOP 1

      END

C H8 - End of Session footer
      SUBROUTINE write_h8 (str)
      IMPLICIT none
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      WRITE (str,1000,err=100)
      RETURN

 100  WRITE(*,*) "Error writing CRD record type h8"
      STOP 1

 1000 FORMAT ("h8")
      END

C H9 - End of File footer
      SUBROUTINE write_h9 (str)
      IMPLICIT none
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      WRITE (str,1000,err=100)
      RETURN

 100  WRITE(*,*) "Error writing CRD record type h9"
      STOP 1

 1000 FORMAT ("h9")
      END

C Ranging data configuration records (1 of n)
C C0 - System Configuration Record
      SUBROUTINE write_c0 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

      write (str,"(a)") "c0"

      if (writei (str, "(i1)", c0_detail_type) .lt. 0) go to 100
      if (writed (str, "(f8.3)", xmit_wavelength) .lt. 0) go to 100

CC This presumes that config_ids that are not used are blank.
      write (str,"(A)",err=100) trim(str)//" "//trim(config_ids(1))
      write (str,"(A)",err=100) trim(str)//" "//trim(config_ids(2))
      write (str,"(A)",err=100) trim(str)//" "//trim(config_ids(3))
      write (str,"(A)",err=100) trim(str)//" "//trim(config_ids(4))
      write (str,"(A)",err=100) trim(str)//" "//trim(config_ids(5))
      write (str,"(A)",err=100) trim(str)//" "//trim(config_ids(6))
      write (str,"(A)",err=100) trim(str)//" "//trim(config_ids(7))

CC      ENDIF
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c0"
      STOP 1

      END

C C1 - Laser Configuration Record
      SUBROUTINE write_c1 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writed, writei
      INTEGER trimlen
      INTEGER lci_len, lt_len
      INCLUDE '../include/crd.inc'

      write (str,"(a)") "c1"

      if (writei (str, "(i1)", c1_detail_type) .lt. 0) go to 100

      write (str,"(A)",err=100) trim(str)//" "//trim(laser_config_id)

      write (str,"(A)",err=100) trim(str)//" "//trim(laser_type)

      if (writed (str, "(f10.2)", prim_wavelength) .lt. 0) go to 100
      if (writed (str, "(f10.2)", nom_fire_rate) .lt. 0) go to 100
      if (writed (str, "(f10.2)", pulse_energy) .lt. 0) go to 100
      if (writed (str, "(f6.1)", pulse_width) .lt. 0) go to 100
      if (writed (str, "(f5.2)", beam_div) .lt. 0) go to 100
      if (writei (str, "(i4)", pulses_in_semitrain) .lt. 0) go to 100
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c1"
      STOP 1

      END

C C2 - Detector Configdetector_type,uration Record
      SUBROUTINE write_c2 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writed, writei
      INTEGER trimlen
      INTEGER dci_len, dt_len, opt_len, sp_len
      INCLUDE '../include/crd.inc'

      dci_len= trimlen(detector_config_id)
      dt_len= trimlen(detector_type)
      opt_len= trimlen(output_pulse_type)
      sp_len= trimlen(signal_proc)
      IF (format_version .EQ. 1) THEN
        WRITE (str,1000,err=100) 
     &           c2_detail_type, detector_config_id(1:dci_len), 
     &           detector_type(1:dt_len), app_wavelength, qe, voltage, 
     &           dark_count, output_pulse_type(1:opt_len), 
     &           output_pulse_width, spectral_filter,
     &           spectral_filter_xmission, spatial_filter, 
     &           signal_proc(1:sp_len)
      ELSE IF (format_version .EQ. 2) THEN
        write (str,"(a)") "c2"
  
        if (writei (str, "(i1)", c2_detail_type) .lt. 0) go to 100
  
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(detector_config_id)
  
        write (str,"(A)",err=100) trim(str)//" "//trim(detector_type)
  
        if (writed (str, "(f10.3)", app_wavelength) .lt. 0) go to 100
        if (writed (str, "(f5.1)", qe) .lt. 0) go to 100
        if (writed (str, "(f6.1)", voltage) .lt. 0) go to 100
        if (writed (str, "(f5.1)", dark_count) .lt. 0) go to 100
  
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(output_pulse_type)
  
        if (writed (str, "(f5.1)", output_pulse_width) .lt. 0) go to 100
        if (writed (str, "(f5.2)", spectral_filter) .lt. 0) go to 100
        if (writed (str, "(f5.1)", spectral_filter_xmission) .lt. 0) 
     &      go to 100
        if (writed (str, "(f5.2)", spatial_filter) .lt. 0) go to 100
        
        write (str,"(A)",err=100) trim(str)//" "//trim(signal_proc)

        if (writed (str, "(f6.1)", amp_gain) .lt. 0) go to 100
        if (writed (str, "(f6.1)", amp_bandwidth) .lt. 0) go to 100
        if (writei (str, "(i2)", amp_in_use) .lt. 0) go to 100

      ENDIF
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c2"
      STOP 1

 1000 FORMAT ("c2 ",i1,1x,a,1x,a,1x,f10.3,1x,f5.1,1x,f6.1,1x,f5.1,
     &          1x,a,1x,f5.1,1x,f5.2,1x,f5.1,1x,f5.2,1x,a)
      END

C C3 - Timing Configuration Record
      SUBROUTINE write_c3 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

        write (str,"(a)") "c3"

        if (writei (str, "(i1)", c3_detail_type) .lt. 0) go to 100
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(timing_config_id)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(time_source)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(freq_source)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(timer)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(timer_serial_num)
        if (writed (str, "(f6.1)", epoch_delay_corr) .lt. 0) 
     &    go to 100
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c3"
      STOP 1

      END

C C4 - Transponder Configuration Record
      SUBROUTINE write_c4 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

CC Warning: There may be problems with NA_VALUE being too small for some of
CC these fields.
        write (str,"(a)") "c4"

        if (writei (str, "(i1)", c4_detail_type) .lt. 0) go to 100
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(xponder_config_id)
        if (writed (str, "(f20.3)", est_stn_utc_offset) .lt. 0) 
     &    go to 100
        if (writed (str, "(f11.2)", est_stn_osc_drift) .lt. 0)
     &    go to 100
        if (writed (str, "(f20.3)", est_xponder_utc_offset) .lt. 0)
     &    go to 100
        if (writed (str, "(f11.2)", xponder_clock_ref_time) .lt. 0)
     &    go to 100
        if (writed (str, "(f20.12)", stn_off_drift_app_ind) .lt. 0)
     &    go to 100
        if (writei (str, "(i1)", stn_off_drift_app_ind) .lt. 0)
     &    go to 100
        if (writei (str, "(i1)", SC_off_drift_app_ind) .lt. 0)
     &    go to 100
        if (writei (str, "(i1)", SC_time_simplified_ind) .lt. 0)
     &    go to 100
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c4"
      STOP 1

      END

C C5 - Timing Configuratiming_config_id,tion Record
C New in v2
      SUBROUTINE write_c5 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

C     No "c5" in version 1
      str= ""
      IF (format_version .EQ. 1) RETURN

        write (str,"(a)") "c5"

        if (writei (str, "(i1)", c5_detail_type) .lt. 0) go to 100
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(software_config_id)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(track_software)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(track_software_versions)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(process_software)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(process_software_versions)
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c5"
      STOP 1

      END

C C6 - Timing Configuration Record
C New in v2
      SUBROUTINE write_c6 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

C     No "c6" in version 1
      str= ""
      IF (format_version .EQ. 1) RETURN

        write (str,"(a)") "c6"

        if (writei (str, "(i1)", c6_detail_type) .lt. 0) go to 100
        write (str,"(A)",err=100) trim(str)//" "//trim(met_config_id)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(press_sensor_manufacturer)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(press_sensor_model)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(press_sensor_serial_num)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(temp_sensor_manufacturer)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(temp_sensor_model)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(temp_sensor_serial_num)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(humid_sensor_manufacturer)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(humid_sensor_model)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(humid_sensor_serial_num)
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c6"
      STOP 1

      END

C C7 - Calibration configuration record
C New in v2
      SUBROUTINE write_c7 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

C     No "c7" in version 1
      str= ""
      IF (format_version .EQ. 1) RETURN
        write (str,"(a)") "c7"

        if (writei (str, "(i1)", c7_detail_type) .lt. 0) go to 100
        write (str,"(A)",err=100) trim(str)//" "//trim(c7_cal_config_id)
        write (str,"(A)",err=100) trim(str)//" "//trim(cal_target_name)
        if (writed (str, "(f12.5)", surveyed_target_dist) .lt. 0) 
     &    go to 100
        if (writed (str, "(f6.2)", cal_survey_error) .lt. 0) go to 100
        if (writed (str, "(f10.5)", other_fixed_delays) .lt. 0) 
     &    go to 100
        if (writed (str, "(f10.2)", cal_pulse_energy) .lt. 0) go to 100
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(cal_processing_software)
        write (str,"(A)",err=100) 
     &         trim(str)//" "//trim(cal_processing_software_version)
CC        write (*,*) "write_c7: ", str
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c7"
      STOP 1

      END

C Ranging data records
C 10 - Range Record
      SUBROUTINE write_10 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writei, writed
      INTEGER trimlen, sci_len
      INCLUDE '../include/crd.inc'

      sci_len= trimlen(d10_sysconfig_id)
      IF (format_version .EQ. 1) THEN
        WRITE (str,1000,err=100) d10_sec_of_day,d10_time_of_flight,
     &           d10_sysconfig_id(1:sci_len), d10_epoch_event, 
     &           filter_flag, d10_detector_channel, stop_number, xcv_amp
      ELSE IF (format_version .EQ. 2) THEN
        write (str,"(a)") "10"
        if (writed (str, "(f18.12)", d10_sec_of_day) .lt. 0) go to 100
        if (writed (str, "(f18.12)", d10_time_of_flight) .lt. 0) 
     &      go to 100
        write (str,"(A)",err=100) trim(str)//" "//trim(d10_sysconfig_id)
        if (writei (str, "(i1)", d10_epoch_event) .lt. 0) go to 100
        if (writei (str, "(i1)", filter_flag) .lt. 0) go to 100
        if (writei (str, "(i1)", d10_detector_channel) .lt. 0) go to 100
        if (writei (str, "(i1)", stop_number) .lt. 0) go to 100
        if (writei (str, "(i5)", xcv_amp) .lt. 0) go to 100
        if (writei (str, "(i5)", xmt_amp) .lt. 0) go to 100
      ENDIF
      RETURN

 100  WRITE(*,*) "Error writing CRD record type 10"
      STOP 1

 1000 FORMAT ("10 ",f18.12,1x,f18.12,1x,a,1x,4(i1,1x),i5)
      END

C 11 - Normal Point Record
      SUBROUTINE write_11 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writei, writed
      INTEGER trimlen, sci_len
      INCLUDE '../include/crd.inc'

      sci_len= trimlen(d11_sysconfig_id)
      IF (format_version .EQ. 1) THEN
        WRITE (str,1000,err=100)  d11_sec_of_day, d11_time_of_flight,
     &           d11_sysconfig_id(1:sci_len), d11_epoch_event, 
     &           np_window_length, num_ranges, bin_rms, bin_skew, 
     &           bin_kurtosis, bin_PmM, return_rate, 
     &           d11_detector_channel
      ELSE IF (format_version .EQ. 2) THEN
        write (str,"(a)") "11"

        if (writed (str, "(f18.12)", d11_sec_of_day) .lt. 0) go to 100
        if (writed (str, "(f18.12)", d11_time_of_flight) .lt. 0) 
     &      go to 100

        write (str,"(A)",err=100) trim(str)//" "//trim(d11_sysconfig_id)
        if (writei (str, "(i1)", d11_epoch_event) .lt. 0) go to 100
        if (writed (str, "(f6.1)", np_window_length) .lt. 0) go to 100
        if (writei (str, "(i6)", num_ranges) .lt. 0) go to 100
        if (writed (str, "(f6.1)", bin_rms) .lt. 0) go to 100
        if (writed (str, "(f7.3)", bin_skew) .lt. 0) go to 100
        if (writed (str, "(f7.3)", bin_kurtosis) .lt. 0) go to 100
        if (writed (str, "(f9.1)", bin_PmM) .lt. 0) go to 100
        if (writed (str, "(f6.2)", return_rate) .lt. 0) go to 100
        if (writei (str, "(i1)", d11_detector_channel) .lt. 0) go to 100
        if (writed (str, "(f9.1)", signal_to_noise) .lt. 0) go to 100
      ENDIF

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 11"
      STOP 1

 1000 FORMAT ("11 ",f18.12,1x,f18.12,1x,a,1x,i1,1x,f6.1,1x,i6,
     &            1x,f6.1,1x,f7.3,1x,f7.3,1x,f9.1,1x,f6.2,1x,i1)
      END

C 12 - Range Supplement Record
      SUBROUTINE write_12 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writei, writed
      INTEGER trimlen, sci_len
      INCLUDE '../include/crd.inc'

      sci_len= trimlen(d12_sysconfig_id)
      IF (format_version .EQ. 1) THEN
        WRITE (str,1000,err=100) d12_sec_of_day, 
     &           d12_sysconfig_id(1:sci_len), refraction_corr, 
     &           target_CofM_corr, nd_value, time_bias
      ELSE IF (format_version .EQ. 2) THEN
        write (str,"(a)") "12"

        if (writed (str, "(f18.12)", d12_sec_of_day) .lt. 0) go to 100

        write (str,"(A)",err=100) trim(str)//" "//trim(d12_sysconfig_id)

        if (writed (str, "(f6.1)", refraction_corr) .lt. 0) go to 100
        if (writed (str, "(f6.4)", target_CofM_corr) .lt. 0) go to 100
        if (writed (str, "(f5.2)", nd_value) .lt. 0) go to 100
        if (writed (str, "(f8.4)", time_bias) .lt. 0) go to 100
        if (writed (str, "(f20.16)", range_rate) .lt. 0) go to 100

      ENDIF

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 12"
      STOP 1

 1000 FORMAT ("12 ",f18.12,1x,a,1x,f6.1,1x,f6.4,1x,f5.2,1x,f8.4)
      END

C 20 - Meteorological Record
      SUBROUTINE write_20 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writei, writed
      INCLUDE '../include/crd.inc'

      write (str,"(a)") "20"

      if (writed (str, "(f9.3)", d20_sec_of_day) .lt. 0) go to 100
      if (writed (str, "(f7.2)", pressure) .lt. 0) go to 100
      if (writed (str, "(f6.2)", temperature) .lt. 0) go to 100
      if (writed (str, "(f4.0)", humidity) .lt. 0) go to 100
      if (writei (str, "(i1)", value_origin) .lt. 0) go to 100

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 20"
      STOP 1

      END

C 21 - Meteorological Supplement Record
      SUBROUTINE write_21 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writei, writed
      INTEGER trimlen, wx_len
      INCLUDE '../include/crd.inc'

      wx_len= trimlen(wx_conditions)
      IF (format_version .EQ. 1) THEN
        WRITE (str,1000,err=100) d21_sec_of_day, wind_speed, 
     &           wind_direction, wx_conditions(1:wx_len), 
     &           visibility, sky_clarity, atmospheric_seeing, 
     &           cloud_cover
      ELSE IF (format_version .EQ. 2) THEN
        write (str,"(a)") "21"

        if (writed (str, "(f9.3)", d21_sec_of_day) .lt. 0) go to 100
        if (writed (str, "(f5.1)", wind_speed) .lt. 0) go to 100
        if (writed (str, "(f5.1)", wind_direction) .lt. 0) go to 100

        write (str,"(A)",err=100) trim(str)//" "//trim(wx_conditions)

        if (writei (str, "(i3)", visibility) .lt. 0) go to 100
        if (writed (str, "(f4.2)", sky_clarity) .lt. 0) go to 100
        if (writei (str, "(i2)", atmospheric_seeing) .lt. 0) go to 100
        if (writei (str, "(i2)", cloud_cover) .lt. 0) go to 100
        if (writed (str, "(f6.2)", sky_temperature) .lt. 0) go to 100
      ENDIF

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 21"
      STOP 1

 1000 FORMAT ("21 ",f9.3,1x,f5.1,1x,f5.1,1x,a,1x,i3,1x,f4.2,1x,
     &        i2,1x,i2)
      END

C 30 - Pointing Angles Record
      SUBROUTINE write_30 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writei, writed
      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        WRITE (str,1000,err=100) d30_sec_of_day, azimuth, elevation, 
     &           direction_ind, angle_origin_ind, refraction_corr_ind
      ELSE IF (format_version .EQ. 2) THEN

        write (str,"(a)") "30"

        if (writed (str, "(f9.3)", d30_sec_of_day) .lt. 0) go to 100
        if (writed (str, "(f8.4)", azimuth) .lt. 0) go to 100
        if (writed (str, "(f8.4)", elevation) .lt. 0) go to 100
        if (writei (str, "(i1)", direction_ind) .lt. 0) go to 100
        if (writei (str, "(i1)", angle_origin_ind) .lt. 0) go to 100
        if (writei (str, "(i1)", refraction_corr_ind) .lt. 0) go to 100
        if (writed (str, "(f10.7)", azimuth_rate) .lt. 0) go to 100
        if (writed (str, "(f10.7)", elevation_rate) .lt. 0) go to 100
      ENDIF

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 30"
      STOP 1

 1000 FORMAT ("30 ",f9.3,2(1x,f8.4),3(1x,i1))
      END

C 40 - Calibration Record
      SUBROUTINE write_40 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER trimlen, sci_len
      INTEGER writei, writed
      INCLUDE '../include/crd.inc'

      sci_len= trimlen(d40_sysconfig_id)
      IF (format_version .EQ. 1) THEN
        WRITE (str,1000,err=100) d40_sec_of_day, 
     &           type_of_data, d40_sysconfig_id(1:sci_len),
     &           num_points_recorded, num_points_used,
     &           one_way_target_dist, cal_sys_delay, cal_delay_shift,
     &           cal_rms, cal_skew, cal_kurtosis, cal_PmM, cal_type_ind,
     &           cal_shift_type_ind, d40_detector_channel
      ELSE IF (format_version .EQ. 2) THEN
        write (str,"(a)") "40"

        if (writed (str, "(f18.12)", d40_sec_of_day) .lt. 0) go to 100
        if (writei (str, "(i1)", type_of_data) .lt. 0) go to 100

        write (str,"(A)",err=100) trim(str)//" "//trim(d40_sysconfig_id)
        
        if (writei (str, "(i8)", num_points_recorded) .lt. 0) go to 100
        if (writei (str, "(i8)", num_points_used) .lt. 0) go to 100
        if (writed (str, "(f7.3)", one_way_target_dist) .lt. 0) 
     &      go to 100
        if (writed (str, "(f10.1)", cal_sys_delay) .lt. 0) go to 100
        if (writed (str, "(f8.1)", cal_delay_shift) .lt. 0) go to 100
        if (writed (str, "(f6.1)", cal_rms) .lt. 0) go to 100
        if (writed (str, "(f7.3)", cal_skew) .lt. 0) go to 100
        if (writed (str, "(f7.3)", cal_kurtosis) .lt. 0) go to 100
        if (writed (str, "(f6.1)", cal_PmM) .lt. 0) go to 100
        if (writei (str, "(i1)", cal_type_ind) .lt. 0) go to 100
        if (writei (str, "(i1)", cal_shift_type_ind) .lt. 0) go to 100
        if (writei (str, "(i1)", d40_detector_channel) .lt. 0) go to 100
        if (writei (str, "(i1)", cal_span) .lt. 0) go to 100
        if (writed (str, "(f6.2)", cal_return_rate) .lt. 0) go to 100
      ENDIF

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 40"
      STOP 1

 1000 FORMAT ("40 ",f18.12,1x,i1,1x,a,1x,i8,1x,i8,1x,f7.3,1x,f10.1,
     &            1x,f8.1,1x,f6.1,1x,f7.3,1x,f7.3,1x,f6.1,3(1x,i1))
      END

C 41 - Calibration Detail Record
      SUBROUTINE write_41 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writei, writed
      INTEGER trimlen, sci_len
      INCLUDE '../include/crd.inc'

C     No "41" in version 1
      str= ""
      IF (format_version .EQ. 1) RETURN

      write (str,"(a)") "41"

      if (writed (str, "(f18.12)", d41_sec_of_day) .lt. 0) go to 100
      if (writei (str, "(i1)", d41_type_of_data) .lt. 0) go to 100

      write (str,"(A)",err=100) trim(str)//" "//trim(d41_sysconfig_id)
      
      if (writei (str, "(i8)", d41_num_points_recorded) .lt. 0) 
     &      go to 100
      if (writei (str, "(i8)", d41_num_points_used) .lt. 0) go to 100
      if (writed (str, "(f7.3)", d41_one_way_target_dist) .lt. 0) 
     &      go to 100
      if (writed (str, "(f10.1)", d41_cal_sys_delay) .lt. 0) go to 100
      if (writed (str, "(f8.1)", d41_cal_delay_shift) .lt. 0) go to 100
      if (writed (str, "(f6.1)", d41_cal_rms) .lt. 0) go to 100
      if (writed (str, "(f7.3)", d41_cal_skew) .lt. 0) go to 100
      if (writed (str, "(f7.3)", d41_cal_kurtosis) .lt. 0) go to 100
      if (writed (str, "(f6.1)", d41_cal_PmM) .lt. 0) go to 100
      if (writei (str, "(i1)", d41_cal_type_ind) .lt. 0) go to 100
      if (writei (str, "(i1)", d41_cal_shift_type_ind) .lt. 0) go to 100
      if (writei (str, "(i1)", d41_detector_channel) .lt. 0) go to 100
      if (writei (str, "(i1)", d41_cal_span) .lt. 0) go to 100
      if (writed (str, "(f6.2)", d41_cal_return_rate) .lt. 0) go to 100
      RETURN

 100  WRITE(*,*) "Error writing CRD record type 41"
      STOP 1

      END

C 42 - Calibration "Shot" Record
      SUBROUTINE write_42 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER trimlen, sci_len, cci_len
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

C     No "42" in version 1
      str= ""
      IF (format_version .EQ. 1) RETURN

      write (str,"(a)") "42"

      if (writed (str, "(f18.12)", d42_sec_of_day) .lt. 0) go to 100
      if (writed (str, "(f18.12)", d42_time_of_flight) .lt. 0) go to 100

      write (str,"(A)",err=100) trim(str)//" "//trim(d42_sysconfig_id)

      write (str,"(A)",err=100) trim(str)//" "//trim(d42_calconfig_id)

      if (writed (str, "(f9.5)", d42_other_variable_delays) .lt. 0) 
     &  go to 100
      if (writei (str, "(i1)", d41_type_of_data) .lt. 0) go to 100
      if (writei (str, "(i1)", d42_cal_type_ind) .lt. 0) go to 100
      if (writei (str, "(i1)", d42_filter_flag) .lt. 0) go to 100
      if (writei (str, "(i1)", d42_detector_channel) .lt. 0) go to 100
      if (writei (str, "(i1)", d42_stop_number) .lt. 0) go to 100
      if (writei (str, "(i1)", d42_cal_span) .lt. 0) go to 100
      if (writei (str, "(i5)", d42_xcv_amp) .lt. 0) go to 100
      if (writei (str, "(i5)", d42_xmt_amp) .lt. 0) go to 100

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 41"
      STOP 1

      END

C 50 - Session Statistics Record
      SUBROUTINE write_50 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER writed, writei
      INCLUDE '../include/crd.inc'

      write (str,"(a)") "50"
      write (str,"(A)") trim(str)//" "//trim(d50_sysconfig_id);

      if (writed (str, "(f6.1)", sess_rms) .lt. 0) go to 100
      if (writed (str, "(f7.3)", sess_skew) .lt. 0) go to 100
      if (writed (str, "(f7.3)", sess_kurtosis) .lt. 0) go to 100
      if (writed (str, "(f6.1)", sess_PmM) .lt. 0) go to 100
      if (writei (str, "(i1)", data_qual_ind) .lt. 0) go to 100

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 50"
      STOP 1

      END

C 60 - Compatibility Record
C V2 -- OBSOLETE
      SUBROUTINE write_60 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER trimlen, sci_len
      INCLUDE '../include/crd.inc'

      str= ""
      IF (format_version .EQ. 1) THEN
        sci_len= trimlen(d60_sysconfig_id)
        WRITE (str,1000,err=100) d60_sysconfig_id(1:sci_len),
     &           sys_change_ind, sys_config_ind
       ENDIF

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 60"
      STOP 1

 1000 FORMAT ("60 ",a,1x,i1,1x,i1)
      END

C 9X - User Defined Records 90-99
      SUBROUTINE write_9x (str)
      IMPLICIT none
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      RETURN
      END

C 00 - Comment Record
      SUBROUTINE write_00 (str)
      IMPLICIT none
      CHARACTER*512 str
      INTEGER trimlen, c_len
      INCLUDE '../include/crd.inc'

C     Maybe should be 77 rather than 80, for record length to be 80
      c_len= MIN (trimlen(comment), 80)
      WRITE (str,1000,err=100) comment(1:c_len)

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 00"
      STOP 1

 1000 FORMAT ("00 ",a)
      END

CC  Write an integer number or "na"
      INTEGER FUNCTION writei (str, ifmt, iv)
      IMPLICIT none
      CHARACTER*512 str
      CHARACTER*20 ifmt
      CHARACTER*256 temp
      LOGICAL iisna
      INTEGER iv
      INCLUDE '../include/crd.inc'

      writei= 0
CC      IF (iv .GT. NA_VALUE) THEN
      IF (.NOT.iisna (iv)) THEN
        WRITE (temp,ifmt,err=100) iv
        WRITE (str,"(A)",err=100) trim(str)//" "//trim(temp)
      ELSE
        WRITE (str,"(A)",err=100) trim(str)//" na"
      ENDIF
      RETURN

100   writei= -1
      RETURN
      END

CC  Write a double precision number or "na"
      integer FUNCTION writed (str, dfmt, dp)
      IMPLICIT none
      CHARACTER*512 str
      CHARACTER*20 dfmt
      CHARACTER*256 temp
      DOUBLE PRECISION dp
      LOGICAL disna
      INCLUDE '../include/crd.inc'

      writed= 0
CC      IF (dp .GT. NA_VALUE) THEN
      if (.NOT.disna (dp)) THEN
        WRITE (temp,dfmt,err=100) dp
        WRITE (str,"(A)",err=100) trim(str)//" "//trim(temp)
      ELSE
        WRITE (str,"(A)",err=100) trim(str)//" na"
      ENDIF
      RETURN

100   writed= -1
      RETURN
      END

CC Is the field "na"?
      LOGICAL FUNCTION iisna (iv)
      IMPLICIT none
      INTEGER iv
      INCLUDE '../include/crd.inc'

      iisna= .FALSE.
      IF (iabs (iv- NA_VALUE) .LT. 1) then
        iisna= .TRUE.
      ENDIF
      RETURN
      END

CC Is the field "na"?
      LOGICAL FUNCTION disna (dp)
      IMPLICIT none
      DOUBLE PRECISION dp
      INCLUDE '../include/crd.inc'

CC      write (*,*) "dp= ",dp, "NA_V= ",NA_VALUE, " NA_VF= ", NA_VALUEF
      disna= .FALSE.
      IF (dabs (dp- NA_VALUEF) .LT. 1.0d-3) then
        disna= .TRUE.
      ENDIF
CC      write (*,*) "dabs= ",dabs (dp- NA_VALUEF), "disna= ",disna
      RETURN
      END
