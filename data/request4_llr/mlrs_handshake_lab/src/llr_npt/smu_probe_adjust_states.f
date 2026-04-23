      SUBROUTINE SMU_PROBE_ADJUST_STATES(TJDINT,SUTCT,EARTHPV,MOONPV,
     1  SUNPV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  Bounded state-seam probe for Request 4.
C  Default is a no-op.
C
C  Optional controls:
C    SMU_STATE_EM_SHIFT_M
C      Constant displacement in meters applied to the Moon position
C      along the instantaneous Earth->Moon line.
C    SMU_STATE_SUN_SHIFT_M
C      Constant displacement in meters applied to the Moon position
C      along the instantaneous Earth->Sun line.
C    SMU_WF_REL_AMP_M
C      Weak-field-inspired relative Earth-Moon displacement amplitude in
C      meters. The helper converts this into a Sun-driven relative
C      state correction:
C        delta r_rel(t) = A * u_ES(t)
C        delta v_rel(t) = A * d(u_ES)/dt
C      and then distributes it between Earth and Moon while preserving the
C      Earth-Moon barycenter.
C
      CHARACTER*64 ENVBUF
      INTEGER LENV, ESTAT, SMU_INITIALIZED
      DOUBLE PRECISION EARTHPV(6), MOONPV(6), SUNPV(6)
      DOUBLE PRECISION SMU_STATE_EM_SHIFT_M, SMU_STATE_SUN_SHIFT_M,
     1  SMU_WF_REL_AMP_M
      DOUBLE PRECISION SMU_AU_M, SMU_SHIFT_AU, VX, VY, VZ, VLEN
      DOUBLE PRECISION UDOTX, UDOTY, UDOTZ, VDOTU, DRELX, DRELY, DRELZ
      DOUBLE PRECISION DVRELX, DVRELY, DVRELZ, MU_E, MU_M
      COMMON /SMU_STATE_CFG/ SMU_STATE_EM_SHIFT_M, SMU_STATE_SUN_SHIFT_M,
     1  SMU_WF_REL_AMP_M, SMU_INITIALIZED
      SAVE /SMU_STATE_CFG/
      DATA SMU_STATE_EM_SHIFT_M /0.D0/, SMU_STATE_SUN_SHIFT_M /0.D0/,
     1  SMU_WF_REL_AMP_M /0.D0/, SMU_INITIALIZED /0/

      IF (SMU_INITIALIZED .EQ. 0) THEN
        ENVBUF = ' '
        LENV = 0
        ESTAT = 1
        CALL GET_ENVIRONMENT_VARIABLE('SMU_STATE_EM_SHIFT_M', ENVBUF,
     1        LENV, ESTAT)
        IF (ESTAT .EQ. 0 .AND. LENV .GT. 0) THEN
          READ(ENVBUF(1:LENV),*,ERR=100) SMU_STATE_EM_SHIFT_M
        ENDIF
  100   CONTINUE

        ENVBUF = ' '
        LENV = 0
        ESTAT = 1
        CALL GET_ENVIRONMENT_VARIABLE('SMU_STATE_SUN_SHIFT_M', ENVBUF,
     1        LENV, ESTAT)
        IF (ESTAT .EQ. 0 .AND. LENV .GT. 0) THEN
          READ(ENVBUF(1:LENV),*,ERR=110) SMU_STATE_SUN_SHIFT_M
        ENDIF
  110   CONTINUE

        ENVBUF = ' '
        LENV = 0
        ESTAT = 1
        CALL GET_ENVIRONMENT_VARIABLE('SMU_WF_REL_AMP_M', ENVBUF,
     1        LENV, ESTAT)
        IF (ESTAT .EQ. 0 .AND. LENV .GT. 0) THEN
          READ(ENVBUF(1:LENV),*,ERR=120) SMU_WF_REL_AMP_M
        ENDIF
  120   CONTINUE

        SMU_INITIALIZED = 1
      ENDIF

      SMU_AU_M = 149597870700.D0
      MU_E = 5.9722D24 / (5.9722D24 + 7.34767309D22)
      MU_M = 7.34767309D22 / (5.9722D24 + 7.34767309D22)

      IF (DABS(SMU_STATE_EM_SHIFT_M) .GT. 0.D0) THEN
        VX = MOONPV(1) - EARTHPV(1)
        VY = MOONPV(2) - EARTHPV(2)
        VZ = MOONPV(3) - EARTHPV(3)
        VLEN = DSQRT(VX*VX + VY*VY + VZ*VZ)
        IF (VLEN .GT. 0.D0) THEN
          SMU_SHIFT_AU = SMU_STATE_EM_SHIFT_M / SMU_AU_M
          MOONPV(1) = MOONPV(1) + SMU_SHIFT_AU * VX / VLEN
          MOONPV(2) = MOONPV(2) + SMU_SHIFT_AU * VY / VLEN
          MOONPV(3) = MOONPV(3) + SMU_SHIFT_AU * VZ / VLEN
        ENDIF
      ENDIF

      IF (DABS(SMU_STATE_SUN_SHIFT_M) .GT. 0.D0) THEN
        VX = SUNPV(1) - EARTHPV(1)
        VY = SUNPV(2) - EARTHPV(2)
        VZ = SUNPV(3) - EARTHPV(3)
        VLEN = DSQRT(VX*VX + VY*VY + VZ*VZ)
        IF (VLEN .GT. 0.D0) THEN
          SMU_SHIFT_AU = SMU_STATE_SUN_SHIFT_M / SMU_AU_M
          MOONPV(1) = MOONPV(1) + SMU_SHIFT_AU * VX / VLEN
          MOONPV(2) = MOONPV(2) + SMU_SHIFT_AU * VY / VLEN
          MOONPV(3) = MOONPV(3) + SMU_SHIFT_AU * VZ / VLEN
        ENDIF
      ENDIF

      IF (DABS(SMU_WF_REL_AMP_M) .GT. 0.D0) THEN
        VX = SUNPV(1) - EARTHPV(1)
        VY = SUNPV(2) - EARTHPV(2)
        VZ = SUNPV(3) - EARTHPV(3)
        VLEN = DSQRT(VX*VX + VY*VY + VZ*VZ)
        IF (VLEN .GT. 0.D0) THEN
          SMU_SHIFT_AU = SMU_WF_REL_AMP_M / SMU_AU_M
          VDOTU = (VX*(SUNPV(4)-EARTHPV(4)) +
     1             VY*(SUNPV(5)-EARTHPV(5)) +
     2             VZ*(SUNPV(6)-EARTHPV(6))) / VLEN
          UDOTX = ((SUNPV(4)-EARTHPV(4)) - (VX/VLEN)*VDOTU) / VLEN
          UDOTY = ((SUNPV(5)-EARTHPV(5)) - (VY/VLEN)*VDOTU) / VLEN
          UDOTZ = ((SUNPV(6)-EARTHPV(6)) - (VZ/VLEN)*VDOTU) / VLEN

          DRELX = SMU_SHIFT_AU * VX / VLEN
          DRELY = SMU_SHIFT_AU * VY / VLEN
          DRELZ = SMU_SHIFT_AU * VZ / VLEN

          DVRELX = SMU_SHIFT_AU * UDOTX
          DVRELY = SMU_SHIFT_AU * UDOTY
          DVRELZ = SMU_SHIFT_AU * UDOTZ

          EARTHPV(1) = EARTHPV(1) - MU_M * DRELX
          EARTHPV(2) = EARTHPV(2) - MU_M * DRELY
          EARTHPV(3) = EARTHPV(3) - MU_M * DRELZ
          MOONPV(1) = MOONPV(1) + MU_E * DRELX
          MOONPV(2) = MOONPV(2) + MU_E * DRELY
          MOONPV(3) = MOONPV(3) + MU_E * DRELZ

          EARTHPV(4) = EARTHPV(4) - MU_M * DVRELX
          EARTHPV(5) = EARTHPV(5) - MU_M * DVRELY
          EARTHPV(6) = EARTHPV(6) - MU_M * DVRELZ
          MOONPV(4) = MOONPV(4) + MU_E * DVRELX
          MOONPV(5) = MOONPV(5) + MU_E * DVRELY
          MOONPV(6) = MOONPV(6) + MU_E * DVRELZ
        ENDIF
      ENDIF

      RETURN
      END
