C**-------------------------------------------------------------------------
C * Subroutines: read CRD data record from an input string
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
C  09/29/08   - Initialize the variable in_arg in the C0 record. (v1.00a
C               rlr)
C  03/10/09   - Record H2 Epoch Timescale corrected from I1 to I2.
C               (v1.00a rlr).
C  03/13/18   - Changes for Version 2.0. rlr.
C  06/26/19   - Add calibration records c7 and 42 to CRD v2.00. rlr.
C  08/01/19   - Add "na" processing (i.i., "na" rather than "-1". v2.01. rlr.
C *
C**-------------------------------------------------------------------------
C Ranging data header/footer records
C H1 - format header
      SUBROUTINE read_h1 (str)
      IMPLICIT none
      CHARACTER*512 str
      CHARACTER*40 tokens(6)
      INTEGER readd, readi
      INCLUDE '../include/crd.inc'

CC      READ (str(4:512),1000,err=100) crd_literal, format_version,
CC      READ (str(4:512),*,err=100) crd_literal, format_version,
CC     &          prod_year, prod_mon, prod_day, prod_hour
      READ (str(3:512),*,err=100) tokens
  
      READ(tokens(1),*,err=100) crd_literal
      if (readi (tokens(2), format_version) .lt. 0) go to 100
      if (readi (tokens(3), prod_year) .lt. 0) go to 100
      if (readi (tokens(4), prod_mon) .lt. 0) go to 100
      if (readi (tokens(5), prod_day) .lt. 0) go to 100
      if (readi (tokens(6), prod_hour) .lt. 0) go to 100

C     Version 0 was for beta-version CRD files prior to v1.
      if (format_version == 0) format_version= 1;
      RETURN

 100  WRITE(*,*) "Error reading CRD record type h1"
      STOP 1

CC 1000 FORMAT (a3,1x,i2,1x,i4,1x,i2,1x,i2,1x,i2)
      END

C H2 - station header
      SUBROUTINE read_h2 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      CHARACTER*512 str
      CHARACTER*512 temp_stn_name, temp_stn_network
      CHARACTER*40 tokens(6)
      INTEGER readd, readi
      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        READ (str(4:512),*,err=100) temp_stn_name,
     &            cdp_pad_id, cdp_sys_num, cdp_occ_num, stn_timescale
        stn_name= temp_stn_name(1:10)
        stn_network= ""
      ELSE IF (format_version .EQ. 2) THEN
CC      READ (str(4:512),1000,err=100) stn_name,
CC        READ (str(4:512),*,err=100) temp_stn_name,
CC     &            cdp_pad_id, cdp_sys_num, cdp_occ_num, stn_timescale,
CC     &            temp_stn_network
CC        tl= MIN (trimlen(temp_stn_name), 10)
CC        stn_name= temp_stn_name(1:tl)
CC        tl= MIN (trimlen(temp_stn_network), 10)
CC        stn_network= temp_stn_network(1:tl)
        READ (str(3:512),*,err=100) tokens
  
        READ(tokens(1),*,err=100) stn_name
CC        write (stn_name, "(A)") trim (tokens(1))
        if (readi (tokens(2), cdp_pad_id) .lt. 0) go to 100
        if (readi (tokens(3), cdp_sys_num) .lt. 0) go to 100
        if (readi (tokens(4), cdp_occ_num) .lt. 0) go to 100
        if (readi (tokens(5), stn_timescale) .lt. 0) go to 100
        READ(tokens(6),*,err=100) stn_network
CC        write (stn_network, "(A)") trim (tokens(6))
CC        write (*,*) ">[",stn_name,"]"
CC        write (*,*) ">[",stn_network,"]"
      ENDIF
      RETURN

 100  WRITE(*,*) "Error reading CRD record type h2"
      STOP 1

CC 1000 FORMAT (a10,1x,i4,1x,i2,1x,i2,1x,i2)
      END

C H3 - spacecraft header
      SUBROUTINE read_h3 (str)
      IMPLICIT none
      CHARACTER*512 str
      CHARACTER*512 temp_target_name
      CHARACTER*40 tokens(7)
      INTEGER trimlen, slen
      INTEGER readd, readi
      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        READ (str(4:512),*,err=100) temp_target_name,
     &          ilrs_id, sic, norad, SC_timescale, target_type
C       try to recover target ckass and loc in case CRD v2 is being written
        target_class= target_type
        IF (target_type .EQ. 2) THEN
          target_class= 1
          target_loc= 3
        ENDIF
        slen= MIN (trimlen(temp_target_name), 10)
        target_name= temp_target_name(1:slen)
      ELSE IF (format_version .EQ. 2) THEN
        READ (str(3:512),*,err=100) tokens
  
        READ(tokens(1),*,err=100) target_name
        if (readi (tokens(2), ilrs_id) .lt. 0) go to 100
        if (readi (tokens(3), sic) .lt. 0) go to 100
        if (readi (tokens(4), norad) .lt. 0) go to 100
        if (readi (tokens(5), SC_timescale) .lt. 0) go to 100
        if (readi (tokens(6), target_class) .lt. 0) go to 100
        if (readi (tokens(7), target_loc) .lt. 0) go to 100

C       try to recover target_type in case CRD v1 is being written
        target_type= target_class
        IF (target_class .EQ. 1 .AND. target_loc .EQ. 3) THEN
          target_type= 2
        ENDIF
      ENDIF
      RETURN

 100  WRITE(*,*) "Error reading CRD record type h3"
      STOP 1

CC 1000 FORMAT (a10,1x,i8,1x,i4,1x,i8,1x,i1,1x,i1)
      END

C H4 - Session header
      SUBROUTINE read_h4 (str)
      IMPLICIT none
      CHARACTER*512 str
      CHARACTER*40 tokens(21)
      INCLUDE '../include/crd.inc'
      INTEGER readd, readi

      READ (str(3:512),*,err=100) tokens

      if (readi (tokens(1), data_type) .lt. 0) go to 100
      if (readi (tokens(2), start_year) .lt. 0) go to 100
      if (readi (tokens(3), start_mon) .lt. 0) go to 100
      if (readi (tokens(4), start_day) .lt. 0) go to 100
      if (readi (tokens(5), start_hour) .lt. 0) go to 100
      if (readi (tokens(6), start_min) .lt. 0) go to 100
      if (readi (tokens(7), start_sec) .lt. 0) go to 100
      if (readi (tokens(8), end_year) .lt. 0) go to 100
      if (readi (tokens(9), end_mon) .lt. 0) go to 100
      if (readi (tokens(10), end_day) .lt. 0) go to 100
      if (readi (tokens(11), end_hour) .lt. 0) go to 100
      if (readi (tokens(12), end_min) .lt. 0) go to 100
      if (readi (tokens(13), end_sec) .lt. 0) go to 100
      if (readi (tokens(14), data_release) .lt. 0) go to 100
      if (readi (tokens(15), refraction_app_ind) .lt. 0) go to 100
      if (readi (tokens(16), CofM_app_ind) .lt. 0) go to 100
      if (readi (tokens(17), xcv_amp_app_ind) .lt. 0) go to 100
      if (readi (tokens(18), stn_sysdelay_app_ind) .lt. 0) go to 100
      if (readi (tokens(19), SC_sysdelay_app_ind) .lt. 0) go to 100
      if (readi (tokens(20), range_type_ind) .lt. 0) go to 100
      if (readi (tokens(21), data_qual_alert_ind) .lt. 0) go to 100
      RETURN

 100  WRITE(*,*) "Error reading CRD record type h4"
      STOP 1

      END

C H5 - Prediction header
C  New in v2.
      SUBROUTINE read_h5 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      INTEGER readd, readi
      CHARACTER*512 str
      CHARACTER*40 tokens(5)
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) tokens

      if (readi (tokens(1), prediction_type) .lt. 0) go to 100
      if (readi (tokens(2), year_of_century) .lt. 0) go to 100
      READ(tokens(3),*,err=100) date_time
      READ(tokens(4),*,err=100) prediction_provider
      if (readi (tokens(5), sequence_number) .lt. 0) go to 100

      RETURN

 100  WRITE(*,*) "Error reading CRD record type h4"
      STOP 1

      END

C H8 - End of Session footer
      SUBROUTINE read_h8 (str)
      IMPLICIT none
      CHARACTER*512 str
      END

C H9 - End of File footer
      SUBROUTINE read_h9 (str)
      IMPLICIT none
      CHARACTER*512 str
      END

C Ranging data configuration records (1 of n)
C C0 - System Configuration Record
      SUBROUTINE read_c0 (str)
      IMPLICIT none
      character*512 str
      character*256 temp_config_ids(7)
      INCLUDE '../include/crd.inc'

      integer n_arg
      integer i, ii
      logical in_arg

C  See how many parameters are on this line.
C  BE VERY CAREFUL that str is declared length 512 in calling program!
      n_arg= 0
      in_arg= .false.
      do ii= 1, len(str)
        if (str(ii:ii) .ne. " ") then
          if (.not.in_arg) then
             in_arg= .true.
           endif
         else
           if (in_arg .eqv. .true.) then
             in_arg= .false.
             n_arg= n_arg+ 1;
           endif
         endif
       enddo

C  Create default values
      do i=1,7
        temp_config_ids(i)= ""
      enddo

C  Read the correct number of parameters
      if (n_arg .eq. 4) then
        READ (str(3:120),*,err=100)
     &          c0_detail_type, xmit_wavelength, 
     &          temp_config_ids(1)
      elseif (n_arg .eq. 5) then
        READ (str(3:120),*,err=100)
     &          c0_detail_type, xmit_wavelength, 
     &          temp_config_ids(1),temp_config_ids(2)
      elseif (n_arg .eq. 6) then
        READ (str(3:120),*,err=100)
     &          c0_detail_type, xmit_wavelength, 
     &          temp_config_ids(1),temp_config_ids(2),
     &          temp_config_ids(3)
      elseif (n_arg .eq. 7) then
        READ (str(3:120),*,err=100)
     &          c0_detail_type, xmit_wavelength, 
     &          temp_config_ids(1),temp_config_ids(2),
     &          temp_config_ids(3),temp_config_ids(4)
      elseif (n_arg .eq. 8) then
        READ (str(3:120),*,err=100)
     &          c0_detail_type, xmit_wavelength, 
     &          temp_config_ids(1),temp_config_ids(2),
     &          temp_config_ids(3),temp_config_ids(4),
     &          temp_config_ids(5)
      elseif (n_arg .eq. 9) then
        READ (str(3:120),*,err=100)
     &          c0_detail_type, xmit_wavelength, 
     &          temp_config_ids(1),temp_config_ids(2),
     &          temp_config_ids(3),temp_config_ids(4),
     &          temp_config_ids(5),temp_config_ids(6)
      else 
        goto 100
      endif

C In the future, there could be more than 7 config types.
      do i=1,7
        config_ids(i)= temp_config_ids(i)(1:40)
      enddo
      return

 100  write(*,*) "Error reading CRD record type c0"
      stop 1
      end

C C1 - Laser Configuration Record
      SUBROUTINE read_c1 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      CHARACTER*40 tokens(9)
      character*512 str
      INTEGER readd, readi
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) tokens

      if (readi (tokens(1), c1_detail_type) .lt. 0) go to 100
      READ(tokens(2),*,err=100) laser_config_id
      READ(tokens(3),*,err=100) laser_type
      if (readd (tokens(4), prim_wavelength) .lt. 0) go to 100
      if (readd (tokens(5), nom_fire_rate) .lt. 0) go to 100
      if (readd (tokens(6), pulse_energy) .lt. 0) go to 100
      if (readd (tokens(7), pulse_width) .lt. 0) go to 100
      if (readd (tokens(8), beam_div) .lt. 0) go to 100
      if (readi (tokens(9), pulses_in_semitrain) .lt. 0) go to 100

      RETURN

 100      write(*,*) "Error reading CRD record type c1"
      stop 1
      END

C C2 - Detector Configuration Record
      SUBROUTINE read_c2 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      CHARACTER*512 str
      CHARACTER*40 tokens(16)
      INTEGER readd, readi
      character*256 temp_detector_config_id, temp_detector_type
      character*256 temp_output_pulse_type, temp_signal_proc

      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        READ (str(3:512),*,err=100) 
     &          c2_detail_type, temp_detector_config_id, 
     &          temp_detector_type, app_wavelength, qe, voltage, 
     &          dark_count, temp_output_pulse_type, output_pulse_width,
     &          spectral_filter, spectral_filter_xmission, 
     &          spatial_filter, temp_signal_proc
        tl= MIN (trimlen(temp_detector_config_id), 40)
        detector_config_id= temp_detector_config_id(1:tl)
        tl= MIN (trimlen(temp_detector_type), 40)
        detector_type= temp_detector_type(1:tl)
        tl= MIN (trimlen(temp_output_pulse_type), 40)
        output_pulse_type= temp_output_pulse_type(1:tl)
        tl= MIN (trimlen(temp_signal_proc), 40)
        signal_proc= temp_signal_proc(1:tl)
      ELSE IF (format_version .EQ. 2) THEN
        READ (str(3:512),*,err=100,end=100) tokens

        if (readi (tokens(1), c2_detail_type) .lt. 0) go to 100
        READ(tokens(2),*,err=100) detector_config_id
        READ(tokens(3),*,err=100) detector_type
        if (readd (tokens(4), app_wavelength) .lt. 0) go to 100
        if (readd (tokens(5), qe) .lt. 0) go to 100
        if (readd (tokens(6), voltage) .lt. 0) go to 100
        if (readd (tokens(7), dark_count) .lt. 0) go to 100
        READ(tokens(8),*,err=100) output_pulse_type
        if (readd (tokens(9), output_pulse_width) .lt. 0) go to 100
        if (readd (tokens(10), spectral_filter) .lt. 0) go to 100
        if (readd (tokens(11), spectral_filter_xmission) .lt. 0) 
     &      go to 100
        if (readd (tokens(12), spatial_filter) .lt. 0) go to 100
        READ(tokens(13),*,err=100) signal_proc
        if (readd (tokens(14), amp_gain) .lt. 0) go to 100
        if (readd (tokens(15), amp_bandwidth) .lt. 0) go to 100
        if (readi (tokens(16), amp_in_use) .lt. 0) go to 100
      ENDIF

      RETURN

 100      write(*,*) "Error reading CRD record type c2"
      stop 1
      end

C C3 - Timing Configuration Record
      SUBROUTINE read_c3 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      character*512 str
      character*256 temp_timing_config_id, temp_time_source
      character*256 temp_freq_source, temp_timer, temp_timer_serial_num
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) 
     &          c3_detail_type, temp_timing_config_id, temp_time_source,
     &          temp_freq_source, temp_timer, temp_timer_serial_num, 
     &          epoch_delay_corr
      tl= MIN (trimlen(temp_timing_config_id), 40)
      timing_config_id= temp_timing_config_id(1:tl)
      tl= MIN (trimlen(temp_time_source), 40)
      time_source= temp_time_source(1:tl)
      tl= MIN (trimlen(temp_freq_source), 40)
      freq_source= temp_freq_source(1:tl)
      tl= MIN (trimlen(temp_timer), 40)
      timer= temp_timer(1:tl)
      tl= MIN (trimlen(temp_timer_serial_num), 40)
      timer_serial_num= temp_timer_serial_num(1:tl)
      RETURN

 100      write(*,*) "Error reading CRD record type c3"
      stop 1
      end

C C4 - Transponder Configuration Record
      SUBROUTINE read_c4 (str)
      IMPLICIT none
      CHARACTER*40 tokens(10)
      INTEGER readd, readi
      character*512 str
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100,end=100) tokens

      if (readi (tokens(1), c4_detail_type) .lt. 0) go to 100
      READ(tokens(2),*,err=100) xponder_config_id
      if (readd (tokens(3), est_stn_utc_offset) .lt. 0) go to 100
      if (readd (tokens(4), est_stn_osc_drift) .lt. 0) go to 100
      if (readd (tokens(5), est_xponder_utc_offset) .lt. 0) go to 100
      if (readd (tokens(6), est_xponder_osc_drift) .lt. 0) go to 100
      if (readd (tokens(7), xponder_clock_ref_time) .lt. 0) go to 100
      if (readi (tokens(8), stn_off_drift_app_ind) .lt. 0) go to 100
      if (readi (tokens(9), SC_off_drift_app_ind) .lt. 0) go to 100
      if (readi (tokens(10), SC_time_simplified_ind) .lt. 0) go to 100
      RETURN

 100      write(*,*) "Error reading CRD record type c4"
      stop 1
      end

C C5 - Software Configuration Record
C New in v2
C Lots of extra work because READ ends with a comma!
      SUBROUTINE read_c5 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      INTEGER i, ifld, gotit
      character*512 str, temp
      character*256 temp_c5_detail_type
      character*256 temp_software_config_id
      character*256 temp_track_software, temp_track_software_versions
      character*256 temp_process_software,temp_process_software_versions
      INCLUDE '../include/crd.inc'

      temp_c5_detail_type=""
      temp_software_config_id=""
      temp_track_software=""
      temp_track_software_versions=""
      temp_process_software=""
      temp_process_software_versions=""
      temp=""
      ifld= 0
      gotit= 0;
      tl= len(str);
      do i=3,tl,1
        if (ifld .eq. 0) then
          if (str(i:i) .eq. " ") then
            go to 5
          else
            ifld= 1
          endif
        endif
        if (ifld.eq.1) then
          if (str(i:i) .eq." " .AND. gotit.eq.1) then
            ifld= 2
            gotit= 0;
            go to 5
          endif
          if (str(i:i) .eq." ") go to 5
          temp_c5_detail_type=
     &      temp_c5_detail_type(1:trimlen(temp_c5_detail_type))
     &      // str(i:i)
CC          write (*,*) "t1 ",temp_c5_detail_type
          gotit= 1;
        endif
        if (ifld.eq.2) then
          if (str(i:i) .eq.' ' .AND. gotit.eq.1) then
            ifld= 3
            gotit= 0;
            go to 5
          endif
          if (str(i:i) .eq.' ') go to 5
          temp_software_config_id=
     &      temp_software_config_id(1:trimlen(temp_software_config_id))
     &      //str(i:i)
CC          write (*,*) "t2 ",temp_software_config_id
          gotit= 1;
        endif
        if (ifld.eq.3) then
          if (str(i:i) .eq.' ' .AND. gotit.eq.1) then
            ifld= 4
            gotit= 0;
            go to 5
          endif
          if (str(i:i) .eq.' ') go to 5
          temp_track_software=
     &      temp_track_software(1:trimlen(temp_track_software))
     &      //str(i:i)
          gotit= 1;
        endif
        if (ifld.eq.4) then
          if (str(i:i) .eq.' ' .AND. gotit.eq.1) then
            ifld= 5
            gotit= 0;
            go to 5
          endif
          if (str(i:i) .eq.' ') go to 5
          temp_track_software_versions=
     &      temp_track_software_versions(1:
     &        trimlen(temp_track_software_versions))
     &        //str(i:i)
          gotit= 1;
        endif
        if (ifld.eq.5) then
          if (str(i:i) .eq.' ' .AND. gotit.eq.1) then
            ifld= 6
            gotit= 0;
            go to 5
          endif
          if (str(i:i) .eq.' ') go to 5
          temp_process_software=
     &      temp_process_software(1:trimlen(temp_process_software))
     &      //str(i:i)
          gotit= 1;
        endif
        if (ifld.eq.6) then
          if (str(i:i) .eq.' ' .AND. gotit.eq.1) then
            gotit= 0;
            go to 10
          endif
          if (str(i:i) .eq.' ') go to 5
          temp_process_software_versions=
     &      temp_process_software_versions(1:
     &        trimlen(temp_process_software_versions))
     &      //str(i:i)
          gotit= 1;
        endif
 5      continue
        enddo

 10   read (temp_c5_detail_type,*) c5_detail_type
      tl= MIN (trimlen(temp_software_config_id), 40)
      software_config_id= temp_software_config_id(1:tl)
      tl= MIN (trimlen(temp_track_software), 40)
      track_software= temp_track_software(1:tl)
      tl= MIN (trimlen(temp_track_software_versions), 40)
      track_software_versions= temp_track_software_versions(1:tl)
      tl= MIN (trimlen(temp_process_software), 40)
      process_software= temp_process_software(1:tl)
      tl= MIN (trimlen(temp_process_software_versions), 40)
      process_software_versions= temp_process_software_versions(1:tl)

      RETURN

 100      write(*,*) "Error reading CRD record type c5"
      stop 1

      end

C C6 - Meteorology Configuration Record
C New in v2
      SUBROUTINE read_c6 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      character*512 str
      character*256 temp_met_config_id
      character*256 
     &        temp_press_sensor_manufacturer, 
     &        temp_press_sensor_model, 
     &        temp_press_sensor_serial_num,
     &        temp_temp_sensor_manufacturer, 
     &        temp_temp_sensor_model, 
     &        temp_temp_sensor_serial_num,
     &        temp_humid_sensor_manufacturer, 
     &        temp_humid_sensor_model, 
     &        temp_humid_sensor_serial_num
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) 
     &        c6_detail_type, temp_met_config_id,
     &        temp_press_sensor_manufacturer, 
     &        temp_press_sensor_model, 
     &        temp_press_sensor_serial_num,
     &        temp_temp_sensor_manufacturer, 
     &        temp_temp_sensor_model, 
     &        temp_temp_sensor_serial_num,
     &        temp_humid_sensor_manufacturer, 
     &        temp_humid_sensor_model, 
     &        temp_humid_sensor_serial_num

      tl= MIN (trimlen(temp_met_config_id), 40)
      met_config_id= temp_met_config_id(1:tl)
CC      write (*,*) "metci: ",met_config_id
      tl= MIN (trimlen(temp_press_sensor_manufacturer), 40)
      press_sensor_manufacturer= temp_press_sensor_manufacturer(1:tl)
      tl= MIN (trimlen(temp_press_sensor_model), 40)
      press_sensor_model= temp_press_sensor_model(1:tl)
      tl= MIN (trimlen(temp_press_sensor_serial_num), 40)
      press_sensor_serial_num= temp_press_sensor_serial_num(1:tl)
      tl= MIN (trimlen(temp_temp_sensor_manufacturer), 40)
      temp_sensor_manufacturer= temp_temp_sensor_manufacturer(1:tl)
      tl= MIN (trimlen(temp_temp_sensor_model), 40)
      temp_sensor_model= temp_temp_sensor_model(1:tl)
      tl= MIN (trimlen(temp_temp_sensor_serial_num), 40)
      temp_sensor_serial_num= temp_temp_sensor_serial_num(1:tl)
      tl= MIN (trimlen(temp_humid_sensor_manufacturer), 40)
      humid_sensor_manufacturer= temp_humid_sensor_manufacturer(1:tl)
      tl= MIN (trimlen(temp_humid_sensor_model), 40)
      humid_sensor_model= temp_humid_sensor_model(1:tl)
      tl= MIN (trimlen(temp_humid_sensor_serial_num), 40)
      humid_sensor_serial_num= temp_humid_sensor_serial_num(1:tl)

      RETURN

 100      write(*,*) "Error reading CRD record type c6"
      stop 1
      end

C C7 - Calibration Configuration Record
C New in v2
      SUBROUTINE read_c7 (str)
      IMPLICIT none
      INTEGER readd, readi
      character*512 str
      CHARACTER *40 tokens(9)
      INCLUDE '../include/crd.inc'

        READ (str(3:512),*,err=100,end=100) tokens

        if (readi (tokens(1), c7_detail_type) .lt. 0) go to 100
        READ(tokens(2),*,err=100) c7_cal_config_id
        READ(tokens(3),*,err=100) cal_target_name
        if (readd (tokens(4), surveyed_target_dist) .lt. 0) go to 100
        if (readd (tokens(5), cal_survey_error) .lt. 0) go to 100
        if (readd (tokens(6), other_fixed_delays) .lt. 0) go to 100
        if (readd (tokens(7), cal_pulse_energy) .lt. 0) go to 100
        READ(tokens(8),*,err=100) cal_processing_software
        READ(tokens(9),*,err=100) cal_processing_software_version

      RETURN

 100      write(*,*) "Error reading CRD record type c7"
      stop 1
      end

C Ranging data records
C V2: added xmit_amp
C 10 - Range Record
      SUBROUTINE read_10 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      INTEGER readd, readi
      CHARACTER *40 tokens(9)
      character*512 str
      character*256 temp_sysconfig_id
      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        READ (str(3:512),*,err=100) d10_sec_of_day,d10_time_of_flight,
     &          temp_sysconfig_id, d10_epoch_event, filter_flag, 
     &          d10_detector_channel, stop_number, xcv_amp
        tl= MIN (trimlen(temp_sysconfig_id), 40)
        d10_sysconfig_id= temp_sysconfig_id(1:tl)
      ELSE IF (format_version .EQ. 2) THEN
        READ (str(3:512),*,err=100,end=100) tokens

        if (readd (tokens(1), d10_sec_of_day) .lt. 0) go to 100
        if (readd (tokens(2), d10_time_of_flight) .lt. 0) go to 100
        READ(tokens(3),*,err=100) d10_sysconfig_id
        if (readi (tokens(4), d10_epoch_event) .lt. 0) go to 100
        if (readi (tokens(5), filter_flag) .lt. 0) go to 100
        if (readi (tokens(6), d10_detector_channel) .lt. 0) go to 100
        if (readi (tokens(7), stop_number) .lt. 0) go to 100
        if (readi (tokens(8), xcv_amp) .lt. 0) go to 100
        if (readi (tokens(9), xmt_amp) .lt. 0) go to 100
      ENDIF

      return

 100      write(*,*) "Error reading CRD record type 10"
C      stop 1
      end

C 11 - Normal Point Record 
C V2: added signal_to_noise
      SUBROUTINE read_11 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      INTEGER readd, readi
      CHARACTER *40 tokens(13)
      character*512 str
      character*256 temp_sysconfig_id
      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        READ (str(3:512),*,err=100) d11_sec_of_day, d11_time_of_flight,
     &          temp_sysconfig_id, d11_epoch_event, np_window_length,
     &          num_ranges, bin_rms, bin_skew, bin_kurtosis, bin_PmM,
     &          return_rate, d11_detector_channel
        tl= MIN (trimlen(temp_sysconfig_id), 40)
        d11_sysconfig_id= temp_sysconfig_id(1:tl)
      ELSE IF (format_version .EQ. 2) THEN
        READ (str(3:512),*,err=100,end=100) tokens

        if (readd (tokens(1), d11_sec_of_day) .lt. 0) go to 100
        if (readd (tokens(2), d11_time_of_flight) .lt. 0) go to 100
        READ(tokens(3),*,err=100) d11_sysconfig_id
        if (readi (tokens(4), d11_epoch_event) .lt. 0) go to 100
        if (readd (tokens(5), np_window_length) .lt. 0) go to 100
        if (readi (tokens(6), num_ranges) .lt. 0) go to 100
        if (readd (tokens(7), bin_rms) .lt. 0) go to 100
        if (readd (tokens(8), bin_skew) .lt. 0) go to 100
        if (readd (tokens(9), bin_kurtosis) .lt. 0) go to 100
        if (readd (tokens(10), bin_PmM) .lt. 0) go to 100
        if (readd (tokens(11), return_rate) .lt. 0) go to 100
        if (readi (tokens(12), d11_detector_channel) .lt. 0) go to 100
        if (readd (tokens(13), signal_to_noise) .lt. 0) go to 100
      ENDIF

      return

 100      write(*,*) "Error reading CRD record type 11"
      stop 1
      end

C 12 - Range Supplement Record
C V2: added range rate
      SUBROUTINE read_12 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      INTEGER readd, readi
      CHARACTER *40 tokens(7)
      character*512 str
      character*256 temp_sysconfig_id
      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        READ (str(3:512),*,err=100) d12_sec_of_day, temp_sysconfig_id,
     &          refraction_corr, target_CofM_corr, nd_value, 
     &          time_bias
        tl= MIN (trimlen(temp_sysconfig_id), 40)
        d12_sysconfig_id= temp_sysconfig_id(1:tl)
      ELSE IF (format_version .EQ. 2) THEN
        READ (str(3:512),*,err=100,end=100) tokens

        if (readd (tokens(1), d12_sec_of_day) .lt. 0) go to 100
        READ(tokens(2),*,err=100) d12_sysconfig_id
        if (readd (tokens(3), refraction_corr) .lt. 0) go to 100
        if (readd (tokens(4), target_CofM_corr) .lt. 0) go to 100
        if (readd (tokens(5), nd_value) .lt. 0) go to 100
        if (readd (tokens(6), time_bias) .lt. 0) go to 100
        if (readd (tokens(7), range_rate) .lt. 0) go to 100
      ENDIF
      return

 100      write(*,*) "Error reading CRD record type 12"
      stop 1
      end

C 20 - Meteorological Record
      SUBROUTINE read_20 (str)
      IMPLICIT none
      INTEGER readd, readi
      CHARACTER *40 tokens(5)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100,end=100) tokens

      if (readd (tokens(1), d20_sec_of_day) .lt. 0) go to 100
      if (readd (tokens(2), pressure) .lt. 0) go to 100
      if (readd (tokens(3), temperature) .lt. 0) go to 100
      if (readd (tokens(4), humidity) .lt. 0) go to 100
      if (readi (tokens(5), value_origin) .lt. 0) go to 100
      return

 100      write(*,*) "Error reading CRD record type 20"
      stop 1
      end

C 21 - Meteorological Supplement Record
C V2: precip_type changed to wx_conditions;
C     added sky_temperature

      SUBROUTINE read_21 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      INTEGER readd, readi
      CHARACTER *40 tokens(9)
      character*512 str
      character*256 temp_precip_type, temp_wx_conditions
      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        READ (str(3:512),*,err=100) d21_sec_of_day, wind_speed, 
     &          wind_direction, temp_precip_type, visibility, 
     &          sky_clarity, atmospheric_seeing, cloud_cover
        tl= MIN (trimlen(temp_precip_type), 40)
        wx_conditions= temp_precip_type(1:tl);
      ELSE IF (format_version .EQ. 2) THEN
        READ (str(3:512),*,err=100,end=100) tokens

        if (readd (tokens(1), d21_sec_of_day) .lt. 0) go to 100
        if (readd (tokens(2), wind_speed) .lt. 0) go to 100
        if (readd (tokens(3), wind_direction) .lt. 0) go to 100
        READ(tokens(4),*,err=100) wx_conditions
        if (readi (tokens(5), visibility) .lt. 0) go to 100
        if (readd (tokens(6), sky_clarity) .lt. 0) go to 100
        if (readi (tokens(7), atmospheric_seeing) .lt. 0) go to 100
        if (readi (tokens(8), cloud_cover) .lt. 0) go to 100
        if (readd (tokens(9), sky_temperature) .lt. 0) go to 100
      ENDIF
      return

 100      write(*,*) "Error reading CRD record type 21"
      stop 1
      end

C 30 - Pointing Angles Record
C V2: added azimuth_rate, elevation_rate
      SUBROUTINE read_30 (str)
      IMPLICIT none
      INTEGER readd, readi
      CHARACTER *40 tokens(8)
      character*512 str
      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        READ (str(3:512),*,err=100) d30_sec_of_day, azimuth, elevation, 
     &          direction_ind, angle_origin_ind, refraction_corr_ind
      ELSE IF (format_version .EQ. 2) THEN
        READ (str(3:512),*,err=100,end=100) tokens

        if (readd (tokens(1), d30_sec_of_day) .lt. 0) go to 100
        if (readd (tokens(2), azimuth) .lt. 0) go to 100
        if (readd (tokens(3), elevation) .lt. 0) go to 100
        if (readi (tokens(4), direction_ind) .lt. 0) go to 100
        if (readi (tokens(5), angle_origin_ind) .lt. 0) go to 100
        if (readi (tokens(6), refraction_corr_ind) .lt. 0) go to 100
        if (readd (tokens(7), azimuth_rate) .lt. 0) go to 100
        if (readd (tokens(8), elevation_rate) .lt. 0) go to 100
      ENDIF
      return

 100      write(*,*) "Error reading CRD record type 30"
      stop 1
      end

C 40 - Calibration Record
      SUBROUTINE read_40 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      INTEGER readd, readi
      CHARACTER*512 str
      CHARACTER*40 tokens(17)
      CHARACTER*256 temp_sysconfig_id
      INCLUDE '../include/crd.inc'

      IF (format_version .EQ. 1) THEN
        READ (str(3:512),*,err=100) d40_sec_of_day, 
     &          type_of_data, temp_sysconfig_id,
     &          num_points_recorded, num_points_used,
     &          one_way_target_dist, cal_sys_delay, cal_delay_shift,
     &          cal_rms, cal_skew, cal_kurtosis, cal_PmM, cal_type_ind,
     &          cal_shift_type_ind, d40_detector_channel
        tl= MIN (trimlen(temp_sysconfig_id), 40)
        d40_sysconfig_id= temp_sysconfig_id(1:tl)
      ELSE IF (format_version .EQ. 2) THEN
        READ (str(3:512),*,err=100,end=100) tokens

        if (readd (tokens(1), d40_sec_of_day) .lt. 0) go to 100
        if (readi (tokens(2), type_of_data) .lt. 0) go to 100
        READ(tokens(3),*,err=100) d40_sysconfig_id
        if (readi (tokens(4), num_points_recorded) .lt. 0) go to 100
        if (readi (tokens(5), num_points_used) .lt. 0) go to 100
        if (readd (tokens(6), one_way_target_dist) .lt. 0) go to 100
        if (readd (tokens(7), cal_sys_delay) .lt. 0) go to 100
        if (readd (tokens(8), cal_delay_shift) .lt. 0) go to 100
        if (readd (tokens(9), cal_rms) .lt. 0) go to 100
        if (readd (tokens(10), cal_skew) .lt. 0) go to 100
        if (readd (tokens(11), cal_kurtosis) .lt. 0) go to 100
        if (readd (tokens(12), cal_PmM) .lt. 0) go to 100
        if (readi (tokens(13), cal_type_ind) .lt. 0) go to 100
        if (readi (tokens(14), cal_shift_type_ind) .lt. 0) go to 100
        if (readi (tokens(15), d40_detector_channel) .lt. 0) go to 100
        if (readi (tokens(16), cal_span) .lt. 0) go to 100
        if (readd (tokens(17), cal_return_rate) .lt. 0) go to 100
      ENDIF

      return

 100      write(*,*) "Error reading CRD record type 40"
      stop 1
      end

C 41 - Calibration Detail Record
      SUBROUTINE read_41 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      INTEGER readd, readi
      CHARACTER*512 str
      CHARACTER*40 tokens(17)
      INCLUDE '../include/crd.inc'
      READ (str(3:512),*,err=100,end=100) tokens
  
      if (readd (tokens(1), d41_sec_of_day) .lt. 0) go to 100
      if (readi (tokens(2), d41_type_of_data) .lt. 0) go to 100
      READ(tokens(3),*,err=100) d41_sysconfig_id
      if (readi (tokens(4), d41_num_points_recorded) .lt. 0) go to 100
      if (readi (tokens(5), d41_num_points_used) .lt. 0) go to 100
      if (readd (tokens(6), d41_one_way_target_dist) .lt. 0) go to 100
      if (readd (tokens(7), d41_cal_sys_delay) .lt. 0) go to 100
      if (readd (tokens(8), d41_cal_delay_shift) .lt. 0) go to 100
      if (readd (tokens(9), d41_cal_rms) .lt. 0) go to 100
      if (readd (tokens(10), d41_cal_skew) .lt. 0) go to 100
      if (readd (tokens(11), d41_cal_kurtosis) .lt. 0) go to 100
      if (readd (tokens(12), d41_cal_PmM) .lt. 0) go to 100
      if (readi (tokens(13), d41_cal_type_ind) .lt. 0) go to 100
      if (readi (tokens(14), d41_cal_shift_type_ind) .lt. 0) go to 100
      if (readi (tokens(15), d41_detector_channel) .lt. 0) go to 100
      if (readi (tokens(16), d41_cal_span) .lt. 0) go to 100
      if (readd (tokens(17), d41_cal_return_rate) .lt. 0) go to 100
      return

 100      write(*,*) "Error reading CRD record type 41"
      stop 1
      end

C 42 - Calibration Detail Record
      SUBROUTINE read_42 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      INTEGER readd, readi
      CHARACTER*512 str
      CHARACTER*40 tokens(13)
      INCLUDE '../include/crd.inc'

CC      d42_calconfig_id= temp_calconfig_id(1:tl)
      READ (str(3:512),*,err=100,end=100) tokens
  
      if (readd (tokens(1), d42_sec_of_day) .lt. 0) go to 100
      if (readd (tokens(2), d42_time_of_flight) .lt. 0) go to 100
      READ(tokens(3),*,err=100) d42_sysconfig_id
      READ(tokens(4),*,err=100) d42_calconfig_id
      if (readd (tokens(5), d42_other_variable_delays) .lt. 0) go to 100
      if (readi (tokens(6), d42_type_of_data) .lt. 0) go to 100
      if (readi (tokens(7), d42_cal_type_ind) .lt. 0) go to 100
      if (readi (tokens(8), d42_filter_flag) .lt. 0) go to 100
      if (readi (tokens(9), d42_detector_channel) .lt. 0) go to 100
      if (readi (tokens(10), d42_stop_number) .lt. 0) go to 100
      if (readi (tokens(11), d42_cal_span) .lt. 0) go to 100
      if (readi (tokens(12), d42_xcv_amp) .lt. 0) go to 100
      if (readi (tokens(13), d42_xmt_amp) .lt. 0) go to 100
      return

 100      write(*,*) "Error reading CRD record type 42"
      stop 1
      end

C 50 - Session Statistics Record
      SUBROUTINE read_50 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      INTEGER readd, readi
      character*512 str
      character*256 temp_sysconfig_id
      CHARACTER*40 tokens(6)
      INCLUDE '../include/crd.inc'
      READ (str(3:512),*,err=100) tokens

      READ(tokens(1),*,err=100) d50_sysconfig_id
      if (readd (tokens(2), sess_rms) .lt. 0) go to 100
      if (readd (tokens(3), sess_skew) .lt. 0) go to 100
      if (readd (tokens(4), sess_kurtosis) .lt. 0) go to 100
      if (readd (tokens(5), sess_PmM) .lt. 0) go to 100
      if (readi (tokens(6), data_qual_ind) .lt. 0) go to 100
      return

 100  write(*,*) "Error reading CRD record type 50"
      stop 1
      end

C 60 - Compatibility Record
C V2 -- OBSOLETE
      SUBROUTINE read_60 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      character*512 str
      character*256 temp_sysconfig_id
      INCLUDE '../include/crd.inc'

C The 60 record should only be used for reformatting of pre_CRD data
      IF (format_version .EQ. 1) THEN
        READ (str(4:512),*,err=100) 
     &          temp_sysconfig_id, sys_change_ind, sys_config_ind
        tl= MIN (trimlen(temp_sysconfig_id), 40)
        d60_sysconfig_id= temp_sysconfig_id(1:tl)
      ENDIF

      return

 100  write(*,*) "Error reading CRD record type 60"
      stop 1
      end

C 9X - User Defined Records 90-99
      SUBROUTINE read_9x (str)
      IMPLICIT none
      character*512 str
      end

C 00 - Comment Record
      SUBROUTINE read_00 (str)
      IMPLICIT none
      INTEGER tl, trimlen
      character*512 str
      character*256 temp_comment
      INCLUDE '../include/crd.inc'

      tl= MIN (trimlen(str), 80)
      comment= str(4:tl)

      return

 100      write(*,*) "Error reading CRD record type 00"
      stop 1
      end

      INTEGER FUNCTION readd (str, fval)
      IMPLICIT none
      LOGICAL isna
      CHARACTER*512 str
      DOUBLE PRECISION fval
      INCLUDE '../include/crd.inc'

      if (isna (str(1:2))) then
        fval= NA_VALUEF
      else
        READ(str,*,err=100) fval
      endif
      readd= 0;
      RETURN

 100  readd= -1
      RETURN
      END

      INTEGER FUNCTION readi (str, ival)
      IMPLICIT none
      LOGICAL isna
      CHARACTER*512 str
      INTEGER ival
      INCLUDE '../include/crd.inc'

      if (isna (str(1:2))) then
        ival= NA_VALUE
      else
        READ(str,*,err=100) ival
      endif
      readi= 0;
      RETURN

 100  readi= -1
      RETURN
      END

C isna - is the token na?
      LOGICAL FUNCTION isna (str)
      CHARACTER str(2)

      isna= .false.
      IF (str(1) .ne. 'n' .and. str(1) .ne. 'N') RETURN
      IF (str(2) .ne. 'a' .and. str(2) .ne. 'A') RETURN
      isna= .true.
      RETURN
      END
