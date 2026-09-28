J/A+A/595/A130   Thermodynamic quantities of molecular hydrogen (Popovas+, 2016)
================================================================================
Partition functions. I. Improved partition functions and thermodynamic
quantities for normal, equilibrium, and ortho and para molecular hydrogen.
    Popovas A., Jorgensen U.G.
    <Astron. Astrophys. 595, A130 (2016)>
    =2016A&A...595A.130P        (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Atomic physics
Keywords: miscellaneous - molecular data - astrochemistry - equation of state

Abstract:
    Hydrogen is the most abundant molecule in the Universe. Its
    thermodynamic quantities dominate the physical conditions in molecular
    clouds, protoplanetary disks, etc. It is also of high interest in
    plasma physics. Therefore thermodynamic data for molecular hydrogen
    have to be as accurate as possible in a wide temperature range.

    We here rigorously show the shortcomings of various simplifications
    that are used to calculate the total internal partition function.
    These shortcomings can lead to errors of up to 40 percent or more in
    the estimated partition function. These errors carry on to
    calculations of thermodynamic quantities. Therefore a more complicated
    approach has to be taken.

    Seven possible simplifications of various complexity are described,
    together with advantages and disadvantages of direct summation of
    experimental values. These were compared to what we consider the most
    accurate and most complete treatment (case 8). Dunham coefficients
    were determined from experimental and theoretical energy levels of a
    number of electronically excited states of H2. Both equilibrium and
    normal hydrogen was taken into consideration.

Description:
    New partition functions for equilibrium, normal, and ortho and para
    hydrogen are calculated and thermodynamic quantities are reported for
    the temperature range 1-20000K. Our results are compared to
    previous estimates in the literature. The calculations are not limited
    to the ground electronic state, but include all bound and quasi-bound
    levels of excited electronic states. Dunham coefficients of these
    states of H_2_ are also reported. Reported internal partition functions
    and thermodynamic quantities in the present work are shown to be more
    accurate than previously available data.

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
h2equil.dat      177    20000   Partition function and thermodynamic quantities
                                 for equilibrium H_2_ in 1K temperature steps,
                                 ranging 1-20000K
h2normal.dat     177    20000   Partition function and thermodynamic quantities
                                 for normal H_2_ in 1K temperature steps,
                                 ranging 1-20000K
h2ortho.dat      177    20000   Partition function and thermodynamic quantities
                                 for ortho-H_2_ in 1K temperature steps,
                                 ranging 1-20000K
h2para.dat       177    20000   Partition function and thermodynamic quantities
                                 for para-H_2_ in 1K temperature steps,
                                 ranging 1-20000K
--------------------------------------------------------------------------------

Byte-by-byte Description of file: h2*.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label        Explanations
--------------------------------------------------------------------------------
   1-  7  F7.1  K       T            Temperature
   9- 29 F21.17 ---     Q            Partition function
  31- 52 E22.17 J       Eint/RT      Internal contribution to the
                                      thermodynamic energy
  54- 74 F21.14 J/mol   H-H(0)       Enthalpy
  76- 95 F20.16 J/K.mol S            Entropy
  97-114 F18.15 J/K.mol Cp           Constant-pressure specific heat
 116-137 E22.17 J/K.mol Cv           Constant-volume specific heat
 139-158 E20.15 J/K.mol [-G-H(0)]/T  Gibbs free energy
 160-177 F18.16 ---     gamma        Adiabatic index
--------------------------------------------------------------------------------

Acknowledgements:
   Andrius Popovas, popovas(at)nbi.ku.dk
   StarPlan & NBI, University of Copenhagen, Denmark

================================================================================
(End)  A. Popovas [Univ. Copenhagen, Denmark], P. Vannier [CDS]    18-Jul-2016
