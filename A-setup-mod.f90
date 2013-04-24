! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
!
!     A: The Setup Module
!     ===================
!     -- Contains parameters and some basic functions.
!
! ----------------------------------------------------------------------------------------------------------------------------------
!
!     Copyright 2012, 2013 Hari Sriskantha (School of Mathematics, University of Edinburgh).
!     This file is part of Magneto.
!
!     Magneto is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!     Magneto is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!   
!     A copy of the GNU General Public License can be found in the folder 'readmes', or at <http://www.gnu.org/licenses/>.
!
!     The development of this software was supported by the Engineering and Physical Sciences Research Council (UK). 
!
! ----------------------------------------------------------------------------------------------------------------------------------
! --- This code is best viewed with a window at least 135 characters wide. ---------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------



module setup

  implicit none


  ! Defining real number precision.
  ! -------------------------------

  integer, parameter :: PREC = selected_real_kind (15, 307)



  ! Settings and initial conditions.
  ! --------------------------------

  character :: TEST_PROBLEM
  character :: OUTPUT_TYPE
  character :: RECONSTRUCT_TYPE
  character :: WAVESPEED_TYPE

  logical :: VARIABLE_DELTAT
  logical :: PRINT_DELTAT
  logical :: PRESSURE_FIX
  logical :: DEBUG_MODE

  integer :: PRINT_FREQ



  ! Defining time and space.
  ! ------------------------

    ! Note: the values of these variables can be set in the <TESxx> subroutines.

  real (PREC) :: DELTAT, FULLTIME        ! dt, total required time
  real (PREC) :: MAX_DELTAT, COURANT     ! maximum possible dt, Courant number
  
  real (PREC) :: XDOMN, YDOMN, ZDOMN     ! size of domain, i.e. x in [0, XDOMN], etc.
  integer     :: XSIZE, YSIZE, ZSIZE     ! number of points in each direction

  integer     :: DIMENSIONS     ! number of dimensions modelled

  integer     :: BOUNDARY      ! number of boundary cells (at each end)
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



  ! Variables for Debug Mode.
  ! -------------------------
  
  integer :: CFL_violations = 0

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





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Load Settings (from start.txt) ----------------------------------------------------------------------- [SET04] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Load_settings ()


  ! Namelist of parameters.
  ! -----------------------

    namelist /parameters/ &
         TEST_PROBLEM, OUTPUT_TYPE, &
         RECONSTRUCT_TYPE, WAVESPEED_TYPE, &
         VARIABLE_DELTAT, PRINT_DELTAT, PRESSURE_FIX, DEBUG_MODE, &
         PRINT_FREQ         



  ! Reading file.
  ! -------------

    open (100, file = "start.txt")

    read(100, nml=parameters)

    close (100)

    return

    

! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Load_settings
! ----------------------------------------------------------------------------------------------------------------------------------



end module setup
