! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.1 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
!
!     A: The Setup Module
!     ===================
!     -- Contains initial settings and some basic functions.
!
! ----------------------------------------------------------------------------------------------------------------------------------
! --- This code is best viewed with a window at least 135 characters wide. ---------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------



module setup

  implicit none


  ! Settings and initial conditions.
  ! --------------------------------

  character, parameter :: TEST_PROBLEM     = 'B'
  character, parameter :: RECONSTRUCT_TYPE = 'P'
  character, parameter :: OUTPUT_TYPE      = 'C'

  logical, parameter :: VARIABLE_DELTAT  = .false.
  logical, parameter :: PRINT_DELTAT     = .true.
  logical, parameter :: MAG_WAVESPEED    = .true.
  logical, parameter :: PRESSURE_FIX     = .true.

  integer, parameter :: PRINT_FREQ = 100

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! [SET_A] GUIDE: Settings and Initial Conditions
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
  !     RECONSTRUCT_TYPE :: This setting determines the method used to calculate the input states for the Riemann problem.
  !     ----------------
  !       = 'T' - Balsara's TVD Scheme with steepening (second-order)
  !       = 'P' - Piecewise Parabolic Method (third-order) - WARNING: not bug-free
  !
  !
  !     OUTPUT_TYPE :: This setting determines the format of the output files.
  !     -----------
  !       = 'A' - 1D output: density; pressure; y-velocity; y-magnetic field.
  !       = 'B' - 2D output: density; thermal pressure; magnetic pressure; velocity squared. 
  !       = 'C' - mixed output (for test E): magnetic field components: z- in 2D; transverse in 1D.
  !
  !
  !     Other Settings
  !     --------------
  !     VARIABLE_DELTAT  - if 'true', DELTAT will change with time according to the CFL condition, in <FLU07>.
  !                      - if 'false', then DELTAT will be fixed at the value of MAX_DELTAT, set in the <TESxx> subroutines.
  !     PRINT_DELTAT     - if 'true', then DELTAT will be stored at each timestep in 'output/timestep.txt'
  !     MAG_WAVESPEED    - if 'true', then the magnetic field will be included in the SL/SR estimates in the Riemann solver.
  !     PRESSURE_FIX     - if 'true', then the pressure and energy will be corrected at the end of each timestep, to
  !                      -   compensate for the new, divergence-free maintaining values of the magnetic field. 
  !
  !     PRINT_FREQ       - the frequency (in number of iterations) at which output files are written.
  !
  ! --------------------------------------------------------------------------------------------------------------------------------



  ! Defining real number precision.
  ! -------------------------------

  integer, parameter :: PREC = selected_real_kind (15, 307)



  ! Defining time and space.
  ! ------------------------
      ! Note: the values of these variables can be set in the <TESxx> subroutines.

  real (PREC) :: DELTAT, FULLTIME        ! dt, total required time
  real (PREC) :: MAX_DELTAT, COURANT     ! maximum possible dt, Courant number
  
  real (PREC) :: XDOMN, YDOMN, ZDOMN     ! size of domain, i.e. x in [0, XDOMN], etc.
  integer     :: XSIZE, YSIZE, ZSIZE     ! number of points in each direction

  integer     :: DIMENSIONS     ! number of dimensions modelled

  integer     :: BOUNDARY      ! number of boundary points (at each end)
  integer     :: BOUNDTYPE     ! the boundary type

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! [SET_B] Guide: Boundary Conditions             
  ! --------------------------------------------------------------------------------------------------------------------------------
  !   When setting initial conditions, the BOUNDTYPE flag below determines the boundary conditions.
  !
  !     Set BOUNDTYPE = 1 for EXTENDED boundaries, i.e. ---123 ... 789---  >>  111123 ... 789999
  !     Set BOUNDTYPE = 2 for PERIODIC boundaries, i.e. ---123 ... 789---  >>  789123 ... 789123
  ! --------------------------------------------------------------------------------------------------------------------------------


  integer :: MHDUNITS     ! the MHD units

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! [SET_C] Guide: Magnetohydrodynamic Units       
  ! --------------------------------------------------------------------------------------------------------------------------------
  !   When setting initial conditions, the MHDUNITS flag below determines what units are used for the magnetic field.
  !
  !     Set MHDUNITS = 1 so that MHDF(1) returns 1.0D0 / dsqrt(4.0D0 * PI). Used for CGS units.
  !     Set MHDUNITS = 2 so that MHDF(1) returns 1.0D0. Used for SI units.
  !
  !   A more detailed explanation is given in readmes/equations.pdf, S1.1
  ! --------------------------------------------------------------------------------------------------------------------------------  



  ! Other constants (determined by program).
  ! ----------------------------------------

  real (PREC) :: DELTAX, DELTAY, DELTAZ     ! dx, dy, dz
  integer     :: XFULL,  YFULL,  ZFULL      ! XFULL = XSIZE + BOUNDARY, etc.

  real (PREC) :: NTIME      ! current time
  integer     :: NSTEPS     ! total number of timesteps

  integer :: MSIZE     ! = max(XSIZE, YSIZE, ZSIZE)
  integer :: MFULL     ! = max(XFULL, YFULL, ZFULL)

  real (PREC) :: BASE_NTIME       ! used for saving output files
  integer     :: BASE_FILE_NO     ! used for saving output files



  ! Useful constants.
  ! -----------------

  real (PREC), parameter :: PI = 3.141592653589793
  real (PREC) :: alpha     ! Used for circularly-polarised Alfven waves.

contains



! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Calculate the Dot Product ------------------------------------------------------------------------------ [SET01] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  function Dotproduct (a, b)

    real (PREC), intent (in), dimension(3) :: a, b
    real (PREC) :: Dotproduct

    Dotproduct = (a(1) * b(1)) + (a(2) * b(2)) + (a(3) * b(3))

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end function Dotproduct
! ----------------------------------------------------------------------------------------------------------------------------------


! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Create Vector (from three scalar values) --------------------------------------------------------------- [SET02] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  function Create_vector (a, b, c)

    real (PREC), intent (in)  :: a, b, c
    real (PREC), dimension(3) :: Create_vector

    Create_vector(1) = a
    Create_vector(2) = b
    Create_vector(3) = c

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end function Create_vector
! ----------------------------------------------------------------------------------------------------------------------------------


! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Create Vector (from seven scalar values) --------------------------------------------------------------- [SET03] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  function Create_vector7 (a, b, c, d, e, f, g)

    real (PREC), intent (in)  :: a, b, c, d, e, f, g
    real (PREC), dimension(7) :: Create_vector7

    Create_vector7(1) = a
    Create_vector7(2) = b
    Create_vector7(3) = c
    Create_vector7(4) = d
    Create_vector7(5) = e
    Create_vector7(6) = f
    Create_vector7(7) = g

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end function Create_vector7
! ----------------------------------------------------------------------------------------------------------------------------------



end module setup
