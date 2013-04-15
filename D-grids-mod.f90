! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
! 
!     D: The Grids Module    
!     ===================
!     -- Main variables: 1D row versions of the fluid and magnetic field variables.
!     -- Main actions: extract 1D row from 3D grid; return 1D row to 3D grid.
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



module grids

  use setup
  use fluid
  implicit none


  ! The working 1D arrays.
  ! ----------------------

  real (PREC), dimension(:), allocatable :: density_1D
  real (PREC), dimension(:), allocatable :: x_velocity_1D
  real (PREC), dimension(:), allocatable :: y_velocity_1D
  real (PREC), dimension(:), allocatable :: z_velocity_1D
  real (PREC), dimension(:), allocatable :: x_momentum_1D
  real (PREC), dimension(:), allocatable :: y_momentum_1D
  real (PREC), dimension(:), allocatable :: z_momentum_1D
  real (PREC), dimension(:), allocatable :: energy_1D
  real (PREC), dimension(:), allocatable :: pressure_1D
  real (PREC), dimension(:), allocatable :: gamma_1D

  real (PREC), dimension(:), allocatable :: x_magfield_1D
  real (PREC), dimension(:), allocatable :: y_magfield_1D
  real (PREC), dimension(:), allocatable :: z_magfield_1D

  real (PREC), dimension(:), allocatable :: flux_y1
  real (PREC), dimension(:), allocatable :: flux_y2
  real (PREC), dimension(:), allocatable :: flux_z1
  real (PREC), dimension(:), allocatable :: flux_z2



  ! Other variables.
  ! ----------------

  logical, dimension(2,3) :: edgecheck

  real (PREC), dimension(:), allocatable :: dx, dy, dz, dt
  real (PREC), dimension(:), allocatable :: dtdx, dtdy, dtdz

  integer :: rowsize, rowfull
  
contains



! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Initialise the Grids --------------------------------------------------------------------------------- [GRI01] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Initialise_grids ()


  ! Allocating the working 1D arrays.
  ! ---------------------------------

    allocate (density_1D   (MFULL))
    allocate (x_velocity_1D(MFULL))
    allocate (y_velocity_1D(MFULL))
    allocate (z_velocity_1D(MFULL))
    allocate (x_momentum_1D(MFULL))
    allocate (y_momentum_1D(MFULL))
    allocate (z_momentum_1D(MFULL))
    allocate (energy_1D    (MFULL))
    allocate (pressure_1D  (MFULL))
    allocate (gamma_1D     (MFULL))

    allocate (x_magfield_1D(MFULL))
    allocate (y_magfield_1D(MFULL))
    allocate (z_magfield_1D(MFULL))

    allocate (flux_y1(MFULL))
    allocate (flux_y2(MFULL))
    allocate (flux_z1(MFULL))
    allocate (flux_z2(MFULL))



  ! Allocating the other variables.
  ! -------------------------------

    allocate (dx(MFULL))
    allocate (dy(MFULL))
    allocate (dz(MFULL))
    allocate (dt(MFULL))

    allocate (dtdx(MFULL))
    allocate (dtdy(MFULL))
    allocate (dtdz(MFULL))
    
    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Initialise_grids
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Clear the 1D Data ------------------------------------------------------------------------------------ [GRI02] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Clear_1D_data ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: m



  ! Clearing 1D working data.
  ! -------------------------

    do m = 1, MFULL

       density_1D(m)    = 0.00D0
       x_velocity_1D(m) = 0.00D0
       y_velocity_1D(m) = 0.00D0
       z_velocity_1D(m) = 0.00D0
       x_momentum_1D(m) = 0.00D0
       y_momentum_1D(m) = 0.00D0
       z_momentum_1D(m) = 0.00D0
       energy_1D(m)     = 0.00D0
       pressure_1D(m)   = 0.00D0
       gamma_1D(m)      = 0.00D0

       x_magfield_1D(m) = 0.00D0
       y_magfield_1D(m) = 0.00D0
       z_magfield_1D(m) = 0.00D0
       
       flux_y1(m) = 0.00D0
       flux_y2(m) = 0.00D0
       flux_z1(m) = 0.00D0
       flux_z2(m) = 0.00D0

    end do



  ! Resetting the boundary check.
  ! -----------------------------

    do m = 1, 3

       edgecheck(1, m) = .true.
       edgecheck(2, m) = .true.

    end do

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Clear_1D_data
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Initialise the Magnetic Field Arrays ----------------------------------------------------------------- [GRI03] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Initialise_magnetic_field_arrays ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, j, k



  ! Clearing the delta_B arrays.
  ! ----------------------------

    do k = 1, ZSIZE
       do j = 1, YSIZE
          do i = 1, XSIZE

             x_magfield_stored(i, j, k) = x_magfield(i, j, k)
             y_magfield_stored(i, j, k) = y_magfield(i, j, k)
             z_magfield_stored(i, j, k) = z_magfield(i, j, k)

             delta_bx(i, j, k) = 0.00D0
             delta_by(i, j, k) = 0.00D0
             delta_bz(i, j, k) = 0.00D0

          end do
       end do
    end do

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Initialise_magnetic_field_arrays
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Extract Row from the x-Direction --------------------------------------------------------------------- [GRI04] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Extract_row_x (j_fixed, k_fixed)

    integer, intent (in) :: j_fixed
    integer, intent (in) :: k_fixed


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, iB



  ! Extracting the inner row (minus boundaries).
  ! --------------------------------------------

    rowsize = XSIZE
    rowfull = XFULL

    if (j_fixed == 1)     edgecheck(1,2) = .false.
    if (j_fixed == YSIZE) edgecheck(2,2) = .false.

    if (k_fixed == 1)     edgecheck(1,3) = .false.
    if (k_fixed == ZSIZE) edgecheck(2,3) = .false.

    do i = 1, rowsize

       iB = i + BOUNDARY

       density_1D   (iB) = density   (i, j_fixed, k_fixed)
       x_velocity_1D(iB) = x_velocity(i, j_fixed, k_fixed)
       y_velocity_1D(iB) = y_velocity(i, j_fixed, k_fixed)
       z_velocity_1D(iB) = z_velocity(i, j_fixed, k_fixed)
       x_momentum_1D(iB) = x_momentum(i, j_fixed, k_fixed)
       y_momentum_1D(iB) = y_momentum(i, j_fixed, k_fixed)
       z_momentum_1D(iB) = z_momentum(i, j_fixed, k_fixed)
       energy_1D    (iB) = energy    (i, j_fixed, k_fixed)
       pressure_1D  (iB) = pressure  (i, j_fixed, k_fixed)
       gamma_1D     (iB) = gamma     (i, j_fixed, k_fixed)

       x_magfield_1D(iB) = x_magfield(i, j_fixed, k_fixed)
       y_magfield_1D(iB) = y_magfield(i, j_fixed, k_fixed)
       z_magfield_1D(iB) = z_magfield(i, j_fixed, k_fixed)

    end do



  ! Extracting the whole-row values.
  ! --------------------------------

    do i = 1, rowfull

       dx(i) = DELTAX
       dy(i) = DELTAY
       dz(i) = DELTAZ
       dt(i) = DELTAT / 2.0D0

       dtdx(i) = dt(i) / dx(i)
       dtdy(i) = dt(i) / dy(i)
       dtdz(i) = dt(i) / dz(i)

    end do



  ! Adding 'boundaries'.
  ! --------------------

    call Add_boundaries (density_1D)     ! --> <GRI07>
    call Add_boundaries (x_velocity_1D)
    call Add_boundaries (y_velocity_1D)
    call Add_boundaries (z_velocity_1D)
    call Add_boundaries (x_momentum_1D)
    call Add_boundaries (y_momentum_1D)
    call Add_boundaries (z_momentum_1D)
    call Add_boundaries (energy_1D)
    call Add_boundaries (pressure_1D)
    call Add_boundaries (gamma_1D)
    
    call Add_boundaries (x_magfield_1D)
    call Add_boundaries (y_magfield_1D)
    call Add_boundaries (z_magfield_1D)

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Extract_row_x
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE 04: Extract Row from the y-Direction ------------------------------------------------------------------ [GRI05] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Extract_row_y (i_fixed, k_fixed)

    integer, intent (in) :: i_fixed
    integer, intent (in) :: k_fixed


  ! Declaration of local variables.
  ! -------------------------------

    integer :: j, jB



  ! Extracting the inner row (minus boundaries).
  ! --------------------------------------------

    rowsize = YSIZE
    rowfull = YFULL

    if (i_fixed == 1)     edgecheck(1,2) = .false.
    if (i_fixed == XSIZE) edgecheck(2,2) = .false.

    if (k_fixed == 1)     edgecheck(1,3) = .false.
    if (k_fixed == ZSIZE) edgecheck(2,3) = .false.

    do j = 1, rowsize

       jB = j + BOUNDARY

       density_1D   (jB) = density   (i_fixed, j, k_fixed)
       x_velocity_1D(jB) = y_velocity(i_fixed, j, k_fixed)
       y_velocity_1D(jB) = x_velocity(i_fixed, j, k_fixed)
       z_velocity_1D(jB) = z_velocity(i_fixed, j, k_fixed)
       x_momentum_1D(jB) = y_momentum(i_fixed, j, k_fixed)
       y_momentum_1D(jB) = x_momentum(i_fixed, j, k_fixed)
       z_momentum_1D(jB) = z_momentum(i_fixed, j, k_fixed)
       energy_1D    (jB) = energy    (i_fixed, j, k_fixed)
       pressure_1D  (jB) = pressure  (i_fixed, j, k_fixed)
       gamma_1D     (jB) = gamma     (i_fixed, j, k_fixed)

       x_magfield_1D(jB) = y_magfield(i_fixed, j, k_fixed)
       y_magfield_1D(jB) = x_magfield(i_fixed, j, k_fixed)
       z_magfield_1D(jB) = z_magfield(i_fixed, j, k_fixed)

    end do



  ! Extracting the whole-row values.
  ! --------------------------------

    do j = 1, rowfull

       dx(j) = DELTAY
       dy(j) = DELTAX
       dz(j) = DELTAZ
       dt(j) = DELTAT / 2.0D0

       dtdx(j) = dt(j) / dx(j)
       dtdy(j) = dt(j) / dy(j)
       dtdz(j) = dt(j) / dz(j)

    end do



  ! Adding 'boundaries'.
  ! --------------------

    call Add_boundaries (density_1D)     ! --> <GRI07>
    call Add_boundaries (x_velocity_1D)
    call Add_boundaries (y_velocity_1D)
    call Add_boundaries (z_velocity_1D)
    call Add_boundaries (x_momentum_1D)
    call Add_boundaries (y_momentum_1D)
    call Add_boundaries (z_momentum_1D)
    call Add_boundaries (energy_1D)
    call Add_boundaries (pressure_1D)
    call Add_boundaries (gamma_1D)
    
    call Add_boundaries (x_magfield_1D)
    call Add_boundaries (y_magfield_1D)
    call Add_boundaries (z_magfield_1D)

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Extract_row_y
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Extract Row from the z-Direction --------------------------------------------------------------------- [GRI06] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Extract_row_z (i_fixed, j_fixed)

    integer, intent (in) :: i_fixed
    integer, intent (in) :: j_fixed


  ! Declaration of local variables.
  ! -------------------------------

    integer :: k, kB



  ! Extracting the inner row (minus boundaries).
  ! --------------------------------------------
    
    rowsize = ZSIZE
    rowfull = ZFULL

    if (i_fixed == 1)     edgecheck(1,2) = .false.
    if (i_fixed == XSIZE) edgecheck(2,2) = .false.

    if (j_fixed == 1)     edgecheck(1,3) = .false.
    if (j_fixed == YSIZE) edgecheck(2,3) = .false.

    do k = 1, rowsize

       kB = k + BOUNDARY

       density_1D   (kB) = density   (i_fixed, j_fixed, k)
       x_velocity_1D(kB) = z_velocity(i_fixed, j_fixed, k)
       y_velocity_1D(kB) = x_velocity(i_fixed, j_fixed, k)
       z_velocity_1D(kB) = y_velocity(i_fixed, j_fixed, k)
       x_momentum_1D(kB) = z_momentum(i_fixed, j_fixed, k)
       y_momentum_1D(kB) = x_momentum(i_fixed, j_fixed, k)
       z_momentum_1D(kB) = y_momentum(i_fixed, j_fixed, k)
       energy_1D    (kB) = energy    (i_fixed, j_fixed, k)
       pressure_1D  (kB) = pressure  (i_fixed, j_fixed, k)
       gamma_1D     (kB) = gamma     (i_fixed, j_fixed, k)

       x_magfield_1D(kB) = z_magfield(i_fixed, j_fixed, k)
       y_magfield_1D(kB) = x_magfield(i_fixed, j_fixed, k)
       z_magfield_1D(kB) = y_magfield(i_fixed, j_fixed, k)

    end do



  ! Extracting the whole-row values.
  ! --------------------------------

    do k = 1, rowfull

       dx(k) = DELTAZ
       dy(k) = DELTAX
       dz(k) = DELTAY
       dt(k) = DELTAT / 2.0D0

       dtdx(k) = dt(k) / dx(k)
       dtdy(k) = dt(k) / dy(k)
       dtdz(k) = dt(k) / dz(k)

    end do



  ! Adding 'boundaries'.
  ! --------------------

    call Add_boundaries (density_1D)     ! --> <GRI07>
    call Add_boundaries (x_velocity_1D)
    call Add_boundaries (y_velocity_1D)
    call Add_boundaries (z_velocity_1D)
    call Add_boundaries (x_momentum_1D)
    call Add_boundaries (y_momentum_1D)
    call Add_boundaries (z_momentum_1D)
    call Add_boundaries (energy_1D)
    call Add_boundaries (pressure_1D)
    call Add_boundaries (gamma_1D)
    
    call Add_boundaries (x_magfield_1D)
    call Add_boundaries (y_magfield_1D)
    call Add_boundaries (z_magfield_1D)
    
    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Extract_row_z
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Add Boundaries --------------------------------------------------------------------------------------- [GRI07] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Add_boundaries (quant)

    real (PREC), dimension(:), intent (inout) :: quant


  ! Declaration of local variables.
  ! -------------------------------

    integer :: m



  ! Adding boundaries (extended).
  ! -----------------------------

    if (BOUNDTYPE == int(1)) then

       do m = 1, BOUNDARY
          
          quant(m) = quant(BOUNDARY + 1)

       end do

       do m = (BOUNDARY + rowsize + 1), rowfull

          quant(m) = quant(BOUNDARY + rowsize)

       end do

    end if



  ! Adding boundaries (periodic).
  ! -----------------------------
    
    if (BOUNDTYPE == int(2)) then

       do m = 1, BOUNDARY

          quant(m) = quant(rowsize + m)
          quant(BOUNDARY + rowsize + m) = quant(BOUNDARY + m)

       end do

    end if

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Add_boundaries
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Return Row from the x-Direction ---------------------------------------------------------------------- [GRI08] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Return_row_x (j_fixed, k_fixed)

    integer, intent (in) :: j_fixed
    integer, intent (in) :: k_fixed


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, iB



  ! Returning the inner row.
  ! ------------------------

    do i = 1, rowsize

       iB = i + BOUNDARY

       density   (i, j_fixed, k_fixed) = density_1D   (iB)
       x_velocity(i, j_fixed, k_fixed) = x_velocity_1D(iB)
       y_velocity(i, j_fixed, k_fixed) = y_velocity_1D(iB)
       z_velocity(i, j_fixed, k_fixed) = z_velocity_1D(iB)
       x_momentum(i, j_fixed, k_fixed) = x_momentum_1D(iB)
       y_momentum(i, j_fixed, k_fixed) = y_momentum_1D(iB)
       z_momentum(i, j_fixed, k_fixed) = z_momentum_1D(iB)
       energy    (i, j_fixed, k_fixed) = energy_1D    (iB)
       gamma     (i, j_fixed, k_fixed) = gamma_1D     (iB)
       pressure  (i, j_fixed, k_fixed) = pressure_1D  (iB)

       x_magfield(i, j_fixed, k_fixed) = x_magfield_1D(iB)
       y_magfield(i, j_fixed, k_fixed) = y_magfield_1D(iB)
       z_magfield(i, j_fixed, k_fixed) = z_magfield_1D(iB)
          

       if (DIMENSIONS > 1) then

          delta_by(i, j_fixed, k_fixed) = delta_by(i, j_fixed, k_fixed) + (2.0D0 * flux_y2(iB))
          delta_bz(i, j_fixed, k_fixed) = delta_bz(i, j_fixed, k_fixed) + (2.0D0 * flux_z2(iB))

          if (edgecheck(1, 2)) then
             delta_bx(i, j_fixed-1, k_fixed) = delta_bx(i, j_fixed-1, k_fixed) + flux_y1(iB)
             delta_by(i, j_fixed-1, k_fixed) = delta_by(i, j_fixed-1, k_fixed) + flux_y2(iB)
          else
             if (BOUNDTYPE == 2) then
                delta_bx(i, YSIZE, k_fixed) = delta_bx(i, YSIZE, k_fixed) + flux_y1(iB)
                delta_by(i, YSIZE, k_fixed) = delta_by(i, YSIZE, k_fixed) + flux_y2(iB)
             end if
          end if
          
          if (edgecheck(2, 2)) then
             delta_bx(i, j_fixed+1, k_fixed) = delta_bx(i, j_fixed+1, k_fixed) - flux_y1(iB)
             delta_by(i, j_fixed+1, k_fixed) = delta_by(i, j_fixed+1, k_fixed) + flux_y2(iB)
          else
             if (BOUNDTYPE == 2) then
                delta_bx(i, 1, k_fixed) = delta_bx(i, 1, k_fixed) - flux_y1(iB)
                delta_by(i, 1, k_fixed) = delta_by(i, 1, k_fixed) + flux_y2(iB)
             end if
          end if
          
          if (edgecheck(1, 3)) then
             delta_bx(i, j_fixed, k_fixed-1) = delta_bx(i, j_fixed, k_fixed-1) + flux_z1(iB)
             delta_bz(i, j_fixed, k_fixed-1) = delta_bz(i, j_fixed, k_fixed-1) + flux_z2(iB)
          else
             if (BOUNDTYPE == 2) then
                delta_bx(i, j_fixed, ZSIZE) = delta_bx(i, j_fixed, ZSIZE) + flux_z1(iB)
                delta_bz(i, j_fixed, ZSIZE) = delta_bz(i, j_fixed, ZSIZE) + flux_z2(iB)
             end if
          end if
          
          if (edgecheck(2, 3)) then
             delta_bx(i, j_fixed, k_fixed+1) = delta_bx(i, j_fixed, k_fixed+1) - flux_z1(iB)
             delta_bz(i, j_fixed, k_fixed+1) = delta_bz(i, j_fixed, k_fixed+1) + flux_z2(iB)
          else
             if (BOUNDTYPE == 2) then
                delta_bx(i, j_fixed, 1) = delta_bx(i, j_fixed, 1) - flux_z1(iB)
                delta_bz(i, j_fixed, 1) = delta_bz(i, j_fixed, 1) + flux_z2(iB)
             end if
          end if
          
       end if

    end do

    return

       

! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Return_row_x
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Return Row from the y-Direction ---------------------------------------------------------------------- [GRI09] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Return_row_y (i_fixed, k_fixed)

    integer, intent (in) :: i_fixed
    integer, intent (in) :: k_fixed


    ! Declaration of local variables.
    ! -------------------------------

    integer :: j, jB



    ! Returning the inner row.
    ! ------------------------

    do j = 1, rowsize

       jB = j + BOUNDARY

       density   (i_fixed, j, k_fixed) = density_1D   (jB)
       x_velocity(i_fixed, j, k_fixed) = y_velocity_1D(jB)
       y_velocity(i_fixed, j, k_fixed) = x_velocity_1D(jB)
       z_velocity(i_fixed, j, k_fixed) = z_velocity_1D(jB)
       x_momentum(i_fixed, j, k_fixed) = y_momentum_1D(jB)
       y_momentum(i_fixed, j, k_fixed) = x_momentum_1D(jB)
       z_momentum(i_fixed, j, k_fixed) = z_momentum_1D(jB)
       energy    (i_fixed, j, k_fixed) = energy_1D    (jB)
       gamma     (i_fixed, j, k_fixed) = gamma_1D     (jB)
       pressure  (i_fixed, j, k_fixed) = pressure_1D  (jB)

       x_magfield(i_fixed, j, k_fixed) = y_magfield_1D(jB)
       y_magfield(i_fixed, j, k_fixed) = x_magfield_1D(jB)
       z_magfield(i_fixed, j, k_fixed) = z_magfield_1D(jB)
       
       if (DIMENSIONS > 1) then

          delta_bx(i_fixed, j, k_fixed) = delta_bx(i_fixed, j, k_fixed) + (2.0D0 * flux_y2(jB))
          delta_bz(i_fixed, j, k_fixed) = delta_bz(i_fixed, j, k_fixed) + (2.0D0 * flux_z2(jB))
          
          if (edgecheck(1, 2)) then
             delta_by(i_fixed-1, j, k_fixed) = delta_by(i_fixed-1, j, k_fixed) + flux_y1(jB)
             delta_bx(i_fixed-1, j, k_fixed) = delta_bx(i_fixed-1, j, k_fixed) + flux_y2(jB)
          else
             if (BOUNDTYPE == 2) then
                delta_by(XSIZE, j, k_fixed) = delta_by(XSIZE, j, k_fixed) + flux_y1(jB)
                delta_bx(XSIZE, j, k_fixed) = delta_bx(XSIZE, j, k_fixed) + flux_y2(jB)
             end if
          end if
          
          if (edgecheck(2, 2)) then
             delta_by(i_fixed+1, j, k_fixed) = delta_by(i_fixed+1, j, k_fixed) - flux_y1(jB)
             delta_bx(i_fixed+1, j, k_fixed) = delta_bx(i_fixed+1, j, k_fixed) + flux_y2(jB)
          else
             if (BOUNDTYPE == 2) then
                delta_by(1, j, k_fixed) = delta_by(1, j, k_fixed) - flux_y1(jB)
                delta_bx(1, j, k_fixed) = delta_bx(1, j, k_fixed) + flux_y2(jB)
             end if
          end if
          
          if (edgecheck(1, 3)) then
             delta_by(i_fixed, j, k_fixed-1) = delta_by(i_fixed, j, k_fixed-1) + flux_z1(jB)
             delta_bz(i_fixed, j, k_fixed-1) = delta_bz(i_fixed, j, k_fixed-1) + flux_z2(jB)
          else
             if (BOUNDTYPE == 2) then
                delta_by(i_fixed, j, ZSIZE) = delta_by(i_fixed, j, ZSIZE) + flux_z1(jB)
                delta_bz(i_fixed, j, ZSIZE) = delta_bz(i_fixed, j, ZSIZE) + flux_z2(jB)
             end if
          end if
          
          if (edgecheck(2, 3)) then
             delta_by(i_fixed, j, k_fixed+1) = delta_by(i_fixed, j, k_fixed+1) - flux_z1(jB)
             delta_bz(i_fixed, j, k_fixed+1) = delta_bz(i_fixed, j, k_fixed+1) + flux_z2(jB)
          else
             if (BOUNDTYPE == 2) then
                delta_by(i_fixed, j, 1) = delta_by(i_fixed, j, 1) - flux_z1(jB)
                delta_bz(i_fixed, j, 1) = delta_bz(i_fixed, j, 1) + flux_z2(jB)
             end if
          end if
          
       end if

    end do
    
    return
    


! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Return_row_y
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Return Row from the z-Direction ---------------------------------------------------------------------- [GRI10] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Return_row_z (i_fixed, j_fixed)

    integer, intent (in) :: i_fixed
    integer, intent (in) :: j_fixed


  ! Declaration of local variables.
  ! -------------------------------

    integer :: k, kB



  ! Returning the inner row.
  ! ------------------------

    do k = 1, rowsize

       kB = k + BOUNDARY

       density   (i_fixed, j_fixed, k) = density_1D   (kB)
       x_velocity(i_fixed, j_fixed, k) = y_velocity_1D(kB)
       y_velocity(i_fixed, j_fixed, k) = z_velocity_1D(kB)
       z_velocity(i_fixed, j_fixed, k) = x_velocity_1D(kB)
       x_momentum(i_fixed, j_fixed, k) = y_momentum_1D(kB)
       y_momentum(i_fixed, j_fixed, k) = z_momentum_1D(kB)
       z_momentum(i_fixed, j_fixed, k) = x_momentum_1D(kB)
       energy    (i_fixed, j_fixed, k) = energy_1D    (kB)
       gamma     (i_fixed, j_fixed, k) = gamma_1D     (kB)
       pressure  (i_fixed, j_fixed, k) = pressure_1D  (kB)

       x_magfield(i_fixed, j_fixed, k) = y_magfield_1D(kB)
       y_magfield(i_fixed, j_fixed, k) = z_magfield_1D(kB)
       z_magfield(i_fixed, j_fixed, k) = x_magfield_1D(kB)
       
          if (DIMENSIONS > 1) then

          delta_bx(i_fixed, j_fixed, k) = delta_bx(i_fixed, j_fixed, k) + (2.0D0 * flux_y2(kB))
          delta_by(i_fixed, j_fixed, k) = delta_by(i_fixed, j_fixed, k) + (2.0D0 * flux_z2(kB))

          if (edgecheck(1, 2)) then
             delta_bz(i_fixed-1, j_fixed, k) = delta_bz(i_fixed-1, j_fixed, k) + flux_y1(kB)
             delta_bx(i_fixed-1, j_fixed, k) = delta_bx(i_fixed-1, j_fixed, k) + flux_y2(kB)
          else
             if (BOUNDTYPE == 2) then
                delta_bz(XSIZE, j_fixed, k) = delta_bz(XSIZE, j_fixed, k) + flux_y1(kB)
                delta_bx(XSIZE, j_fixed, k) = delta_bx(XSIZE, j_fixed, k) + flux_y2(kB)
             end if
          end if
          
          if (edgecheck(2, 2)) then
             delta_bz(i_fixed+1, j_fixed, k) = delta_bz(i_fixed+1, j_fixed, k) - flux_y1(kB)
             delta_bx(i_fixed+1, j_fixed, k) = delta_bx(i_fixed+1, j_fixed, k) + flux_y2(kB)
          else
             if (BOUNDTYPE == 2) then
                delta_bz(1, j_fixed, k) = delta_bz(1, j_fixed, k) - flux_y1(kB)
                delta_bx(1, j_fixed, k) = delta_bx(1, j_fixed, k) + flux_y2(kB)
             end if
          end if
          
          if (edgecheck(1, 3)) then
             delta_bz(i_fixed, j_fixed-1, k) = delta_bz(i_fixed, j_fixed-1, k) + flux_z1(kB)
             delta_by(i_fixed, j_fixed-1, k) = delta_by(i_fixed, j_fixed-1, k) + flux_z2(kB)
          else
             if (BOUNDTYPE == 2) then
                delta_bz(i_fixed, YSIZE, k) = delta_bz(i_fixed, YSIZE, k) + flux_z1(kB)
                delta_by(i_fixed, YSIZE, k) = delta_by(i_fixed, YSIZE, k) + flux_z2(kB)
             end if
          end if
          
          if (edgecheck(2, 3)) then
             delta_bz(i_fixed, j_fixed+1, k) = delta_bz(i_fixed, j_fixed+1, k) - flux_z1(kB)
             delta_by(i_fixed, j_fixed+1, k) = delta_by(i_fixed, j_fixed+1, k) + flux_z2(kB)
          else
             if (BOUNDTYPE == 2) then
                delta_bz(i_fixed, 1, k) = delta_bz(i_fixed, 1, k) - flux_z1(kB)
                delta_by(i_fixed, 1, k) = delta_by(i_fixed, 1, k) + flux_z2(kB)
             end if
          end if

       end if
       
    end do

    return

       

! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Return_row_z
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Deallocate the Grids --------------------------------------------------------------------------------- [GRI11] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Deallocate_grids ()


  ! Deallocating the working 1D arrays.
  ! -----------------------------------

    deallocate (density_1D)
    deallocate (x_velocity_1D)
    deallocate (y_velocity_1D)
    deallocate (z_velocity_1D)
    deallocate (x_momentum_1D)
    deallocate (y_momentum_1D)
    deallocate (z_momentum_1D)
    deallocate (energy_1D)
    deallocate (pressure_1D)
    deallocate (gamma_1D)

    deallocate (x_magfield_1D)
    deallocate (y_magfield_1D)
    deallocate (z_magfield_1D)

    deallocate (flux_y1)
    deallocate (flux_y2)
    deallocate (flux_z1)
    deallocate (flux_z2)



  ! Deallocating the other variables.
  ! ---------------------------------

    deallocate (dx)
    deallocate (dy)
    deallocate (dz)
    deallocate (dt)
    deallocate (dtdx)
    deallocate (dtdy)
    deallocate (dtdz)

    return

  

! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Deallocate_grids
! ----------------------------------------------------------------------------------------------------------------------------------



end module grids
