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
C
      CHARACTER*64 ENVBUF
      INTEGER LENV, ESTAT, SMU_INITIALIZED
      DOUBLE PRECISION EARTHPV(6), MOONPV(6), SUNPV(6)
      DOUBLE PRECISION SMU_STATE_EM_SHIFT_M, SMU_STATE_SUN_SHIFT_M
      DOUBLE PRECISION SMU_AU_M, SMU_SHIFT_AU, VX, VY, VZ, VLEN
      COMMON /SMU_STATE_CFG/ SMU_STATE_EM_SHIFT_M, SMU_STATE_SUN_SHIFT_M,
     1  SMU_INITIALIZED
      SAVE /SMU_STATE_CFG/
      DATA SMU_STATE_EM_SHIFT_M /0.D0/, SMU_STATE_SUN_SHIFT_M /0.D0/,
     1  SMU_INITIALIZED /0/

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

        SMU_INITIALIZED = 1
      ENDIF

      SMU_AU_M = 149597870700.D0

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

      RETURN
      END
