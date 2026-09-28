J/ApJS/216/21       Matrix elements and energies for Hydrogen       (Babb, 2015)
================================================================================
State resolved data for radiative association of H and H^+^ and for
photodissociation of H_2_^+^.
    Babb J.F.
   <Astrophys. J. Suppl. Ser., 216, 21 (2015)>
   =2015ApJS..216...21B
================================================================================
ADC_Keywords: Atomic physics
Keywords: early universe; molecular data; molecular processes; Sun: atmosphere;
          white dwarfs

Abstract:
    The matrix elements and energies needed to calculate
    vibrational-rotational state resolved cross sections and rate
    coefficients for radiative association of H and H^+^ and for
    photodissociation of H_2_^+^ are presented for applications to
    simulations of chemistry in the early universe and to stellar
    atmospheres.

File Summary:
--------------------------------------------------------------------------------
 FileName    Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe          80        .   This file
table1.dat      48   284168  *Calculated values of the squared matrix element
--------------------------------------------------------------------------------
Note on table1.dat: Calculated values of the squared matrix element,
  M^2^_{v,N;k,N}_, for bound vibrational-rotational levels (v,N) of the
  1s{sigma}_g_ state of H^+^_2_ with eigenvalue |E_v,N_| to a continuum level
  of the 2p{sigma}_u_ state of relative kinetic energy E for a photon of
  energy h{nu}
--------------------------------------------------------------------------------

See also:
 J/MNRAS/403/1541 : Lyman & Werner molecular hydrogen transitions (Malec+, 2010)
 J/ApJ/667/1119 : Photoionization of H in strong magnetic fields (Zhao+, 2007)
 J/A+AS/141/297 : H_2_ total transition probability (Abgrall+, 2000)
 J/A+AS/114/363 : Polarization profiles of H lines (Casini 1995)

Byte-by-byte Description of file: table1.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label   Explanations
--------------------------------------------------------------------------------
   1-  2  I2    ---     v       [0/19] The v vibrational level
   4-  5  I2    ---     N       [0/23] The N rotational level
   7- 16  E10.4 ---     E       [0/0.83] Relative kinetic energy in atomic units
  18- 25  F8.2  cm-1    E(v,N)  [0.2/21376] Bound state energy
  27- 37  F11.2 nm      lambda  [55/46094692] Photon wavelength
  39- 48  E10.4 ---     M2      Matrix element squared of the transition dipole
                                  moment in atomic units (1)
--------------------------------------------------------------------------------
Note (1): As in equation 6 (section 2).
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

================================================================================
(End)                 Greg Schwarz [AAS], Emmanuelle Perret [CDS]    17-Feb-2015
