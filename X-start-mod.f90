! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
!
!     X: The Start Module
!     ===================
!     -- Contains settings and initial conditions.
!
! ----------------------------------------------------------------------------------------------------------------------------------
! --- This code is best viewed with a window at least 135 characters wide. ---------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------



module start

  implicit none


  ! Settings and initial conditions.
  ! --------------------------------

  character, parameter :: TEST_PROBLEM = 'B'
  character, parameter :: OUTPUT_TYPE  = 'C'

  character, parameter :: RECONSTRUCT_TYPE = 'C'
  character, parameter :: WAVESPEED_TYPE   = 'J'
  
  logical, parameter :: VARIABLE_DELTAT = .true.
  logical, parameter :: PRINT_DELTAT    = .true.
  logical, parameter :: PRESSURE_FIX    = .true.

  integer, parameter :: PRINT_FREQ = 100

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! [STA_A] GUIDE: Settings and Initial Conditions
  ! --------------------------------------------------------------------------------------------------------------------------------
  !   Here are all the main options for the program.
  !
  !     TEST_PROBLEM :: This setting determines the initial conditions. The options are:
  !     ------------
  !       = 'A' - Sod's 1D Shock Tube Test (magnetic)
  !       = 'B' - Circularly Polarised Alfven Waves
  !       = 'C' - The Blast Problem
  !       = 'D' - The Orszag-Tang Vortex
  !       = 'E' - The MHD Rotor Problem
  !       = 'F' - The Kelvin-Helmholtz Instability
  !       = 'X' - Load data from restart file
  !     Alternatively, you can create your own subroutine in C-tests-mod.f90, using the existing ones as templates
  !     (and make sure to update the code in G-evolve-mod and H-error-mod so the program recognises any new options...)
  !
  !
  !     OUTPUT_TYPE :: This setting determines the format of the output files.
  !     -----------
  !       = 'A' - 1D output: density; pressure; y-velocity; y-magnetic field.
  !       = 'B' - 2D output: density; thermal pressure; magnetic pressure; velocity squared. 
  !       = 'C' - mixed output (for test E): magnetic field components: z- in 2D; transverse in 1D.
  !
  ! 
  !     RECONSTRUCT_TYPE :: This setting determines the method used to calculate the input states for the Riemann problem.
  !     ----------------
  !       = 'T' - Balsara's TVD Scheme with steepening (second-order)
  !       = 'P' - Piecewise Parabolic Method (third-order) - based on Stone et al., Astrophys. J. Suppl. S., 178, 137 (2008)
  !       = 'C' - Piecewise Parabolic Method (third-order) - based on Colella and Sekora, J. Comput. Phys., 227, 7069 (2008)
  !
  !
  !     WAVESPEED_TYPE :: This setting determines the method used to estimate the SL/SR wavespeeds in the HLLD Riemann solver.
  !     --------------
  !       = 'M' - The method of Miyoshi and Kusano, J. Comput. Phys., 208, 315 (2005) -- eq. (67)
  !       = 'J' - The method of Janhunen, J. Comput. Phys., 160, 649 (2000) -- eq. (27)
  !
  !
  !     Other Settings
  !     --------------
  !     VARIABLE_DELTAT  - if 'true', DELTAT will change with time according to the CFL condition, in <FLU07>.
  !                      - if 'false', then DELTAT will be fixed at the value of MAX_DELTAT, set in the <TESxx> subroutines.
  !     PRINT_DELTAT     - if 'true', then DELTAT will be stored at each timestep in 'output/timestep.txt'
  !     PRESSURE_FIX     - if 'true', then the pressure and energy will be corrected at the end of each timestep, to
  !                      -   compensate for the new, divergence-free maintaining values of the magnetic field. 
  !
  !     PRINT_FREQ       - the frequency (in number of iterations) at which output files are written.
  !
  ! --------------------------------------------------------------------------------------------------------------------------------



end module start
