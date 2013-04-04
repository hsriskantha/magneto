! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
!
!     A: The Setup Module
!     ===================
!     -- Contains parameters and some basic functions.
!
! ----------------------------------------------------------------------------------------------------------------------------------
! --- This code is best viewed with a window at least 135 characters wide. ---------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------



module setup

  implicit none


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

  integer :: MHDUNITS     ! the MHD units

 

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
