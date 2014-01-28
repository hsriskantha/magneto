! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
!   
!     C: The Tests Module
!     ===================
!     -- Main actions: define initial conditions.
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



module tests

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! [TES_A] Guide: Boundary Conditions             
  ! --------------------------------------------------------------------------------------------------------------------------------
  !   When setting initial conditions, the BOUNDTYPE flag below determines the boundary conditions.
  !
  !     Set BOUNDTYPE = 1 for EXTENDED boundaries, i.e. ---123 ... 789---  >>  111123 ... 789999
  !     Set BOUNDTYPE = 2 for PERIODIC boundaries, i.e. ---123 ... 789---  >>  789123 ... 789123
  ! --------------------------------------------------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! [TES_B] Guide: Magnetohydrodynamic Units       
  ! --------------------------------------------------------------------------------------------------------------------------------
  !   When setting initial conditions, the MHDUNITS flag below determines what units are used for the magnetic field.
  !
  !     Set MHDUNITS = 1 so that MHDF(1) returns 1.0D0 / dsqrt(4.0D0 * PI). Used for CGS units.
  !     Set MHDUNITS = 2 so that MHDF(1) returns 1.0D0. Used for SI units.
  ! --------------------------------------------------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! [TES_C] GUIDE: More Test Cases...              
  ! --------------------------------------------------------------------------------------------------------------------------------
  !   A comphrensive selection of tests can be found on the ATHENA code website:
  !     Standard (C): http://www.astro.princeton.edu/~jstone/Athena/tests/index.html
  !     Fortran:      http://www.astro.virginia.edu/VITA/ATHENA/athena_testsuite.html
  !   ...from which several of these tests have been taken. 
  ! --------------------------------------------------------------------------------------------------------------------------------


  use setup
  use fluid
  implicit none

contains



! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Set Initial Conditions of Fluid ---------------------------------------------------------------------- [TES01] ---
! ----------------------------------------------------------------------------------------------------------------------------------
! --- Version A: Brio and Wu's 1D Shock Tube Test ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Set_initial_conditions_A ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, x_zero
    real (PREC) :: XZERO
    
    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd



  ! Time and space.
  ! ---------------

    DELTAT   = 0.0001D0
    FULLTIME = 0.08D0

    COURANT    = 0.25D0
    MAX_DELTAT = 0.01D0

    DIMENSIONS = 1
    BOUNDARY   = 4

    XSIZE = 512
    YSIZE = 1
    ZSIZE = 1

    XDOMN = 1.0D0
    YDOMN = 1.0D0 * (real(YSIZE) / real(XSIZE))
    ZDOMN = 1.0D0 * (real(ZSIZE) / real(XSIZE))

    MHDUNITS = 2     ! see <TES_B>

    XZERO = 0.5D0



  ! Initialising fluid.
  ! -------------------

    BASE_FILE_NO = 0
    BASE_NTIME   = 0.00D0

    call Initialise_fluid ()     ! --> <FLU01> 



  ! Setting initial conditions.
  ! ---------------------------

    x_zero = ceiling(XZERO * real(XSIZE))

    do i = 1, x_zero

       ! Note: these values need setting. 

       density(i, 1, 1)  = 1.0D0
       pressure(i, 1, 1) = 1.0D0

       x_velocity(i, 1, 1) = 0.0D0
       y_velocity(i, 1, 1) = 0.0D0
       z_velocity(i, 1, 1) = 0.0D0

       gamma(i, 1, 1) = 2.0D0

       x_magfield(i, 1, 1) = 0.75D0 
       y_magfield(i, 1, 1) = 1.00D0 
       z_magfield(i, 1, 1) = 0.0D0


       ! Note: these values are set automatically.

       x_momentum(i, 1, 1) = density(i, 1, 1) * x_velocity(i, 1, 1)
       y_momentum(i, 1, 1) = density(i, 1, 1) * y_velocity(i, 1, 1)
       z_momentum(i, 1, 1) = density(i, 1, 1) * z_velocity(i, 1, 1)

       magf = Create_vector(x_magfield(i, 1, 1), y_magfield(i, 1, 1), z_magfield(i, 1, 1))
       velc = Create_vector(x_velocity(i, 1, 1), y_velocity(i, 1, 1), z_velocity(i, 1, 1))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure(i, 1, 1) = Calculate_total_pressure (pressure(i, 1, 1), magf_squrd)
       energy(i, 1, 1) = Calculate_energy_EOS (density(i, 1, 1), pressure(i, 1, 1), gamma(i, 1, 1), velc_squrd, magf_squrd)

    end do


    do i = x_zero + 1, XSIZE

       ! Note: these values need setting. 

       density(i, 1, 1)  = 0.125D0
       pressure(i, 1, 1) = 0.1D0

       x_velocity(i, 1, 1) = 0.0D0
       y_velocity(i, 1, 1) = 0.0D0
       z_velocity(i, 1, 1) = 0.0D0

       gamma(i, 1, 1) = 2.0D0

       x_magfield(i, 1, 1) = 0.75D0 
       y_magfield(i, 1, 1) = -1.0D0 
       z_magfield(i, 1, 1) = 0.0D0


       ! Note: these values are set automatically.

       x_momentum(i, 1, 1) = density(i, 1, 1) * x_velocity(i, 1, 1)
       y_momentum(i, 1, 1) = density(i, 1, 1) * y_velocity(i, 1, 1)
       z_momentum(i, 1, 1) = density(i, 1, 1) * z_velocity(i, 1, 1)

       magf = Create_vector(x_magfield(i, 1, 1), y_magfield(i, 1, 1), z_magfield(i, 1, 1))
       velc = Create_vector(x_velocity(i, 1, 1), y_velocity(i, 1, 1), z_velocity(i, 1, 1))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure(i, 1, 1) = Calculate_total_pressure (pressure(i, 1, 1), magf_squrd)
       energy(i, 1, 1) = Calculate_energy_EOS (density(i, 1, 1), pressure(i, 1, 1), gamma(i, 1, 1), velc_squrd, magf_squrd)

    end do



  ! Set the boundary conditions.
  ! ----------------------------

    BOUNDTYPE = 1     ! see <TES_A>

    return
    
    

! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Set_initial_conditions_A
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Set Initial Conditions of Fluid ---------------------------------------------------------------------- [TES02] ---
! ----------------------------------------------------------------------------------------------------------------------------------
! --- Version B: Ryu and Jones' Shock Tube Test (2A) -------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Set_initial_conditions_B ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, x_zero
    real (PREC) :: XZERO
    
    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd



  ! Time and space.
  ! ---------------

    DELTAT   = 0.0001D0
    FULLTIME = 0.2D0

    COURANT    = 0.25D0
    MAX_DELTAT = 0.01D0

    DIMENSIONS = 1
    BOUNDARY   = 4

    XSIZE = 512
    YSIZE = 1
    ZSIZE = 1

    XDOMN = 1.0D0
    YDOMN = 1.0D0 * (real(YSIZE) / real(XSIZE))
    ZDOMN = 1.0D0 * (real(ZSIZE) / real(XSIZE))

    MHDUNITS = 2     ! see <TES_B>

    XZERO = 0.5D0



  ! Initialising fluid.
  ! -------------------

    BASE_FILE_NO = 0
    BASE_NTIME   = 0.00D0

    call Initialise_fluid ()     ! --> <FLU01> 



  ! Setting initial conditions.
  ! ---------------------------

    x_zero = ceiling(XZERO * real(XSIZE))

    do i = 1, x_zero

       ! Note: these values need setting. 

       density(i, 1, 1)  = 1.08D0
       pressure(i, 1, 1) = 0.95D0

       x_velocity(i, 1, 1) = 1.2D0
       y_velocity(i, 1, 1) = 0.01D0
       z_velocity(i, 1, 1) = 0.5D0

       gamma(i, 1, 1) = 5.0D0 / 3.0D0

       x_magfield(i, 1, 1) = 2.0D0 / sqrt(4.0D0 * PI)
       y_magfield(i, 1, 1) = 3.6D0 / sqrt(4.0D0 * PI)
       z_magfield(i, 1, 1) = 2.0D0 / sqrt(4.0D0 * PI)


       ! Note: these values are set automatically.

       x_momentum(i, 1, 1) = density(i, 1, 1) * x_velocity(i, 1, 1)
       y_momentum(i, 1, 1) = density(i, 1, 1) * y_velocity(i, 1, 1)
       z_momentum(i, 1, 1) = density(i, 1, 1) * z_velocity(i, 1, 1)

       magf = Create_vector(x_magfield(i, 1, 1), y_magfield(i, 1, 1), z_magfield(i, 1, 1))
       velc = Create_vector(x_velocity(i, 1, 1), y_velocity(i, 1, 1), z_velocity(i, 1, 1))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure(i, 1, 1) = Calculate_total_pressure (pressure(i, 1, 1), magf_squrd)
       energy(i, 1, 1) = Calculate_energy_EOS (density(i, 1, 1), pressure(i, 1, 1), gamma(i, 1, 1), velc_squrd, magf_squrd)

    end do


    do i = x_zero + 1, XSIZE

       ! Note: these values need setting. 

       density(i, 1, 1)  = 1.0D0
       pressure(i, 1, 1) = 1.0D0

       x_velocity(i, 1, 1) = 0.0D0
       y_velocity(i, 1, 1) = 0.0D0
       z_velocity(i, 1, 1) = 0.0D0

       gamma(i, 1, 1) = 5.0D0 / 3.0D0

       x_magfield(i, 1, 1) = 2.0D0 / sqrt(4.0D0 * PI)
       y_magfield(i, 1, 1) = 4.0D0 / sqrt(4.0D0 * PI)
       z_magfield(i, 1, 1) = 2.0D0 / sqrt(4.0D0 * PI)


       ! Note: these values are set automatically.

       x_momentum(i, 1, 1) = density(i, 1, 1) * x_velocity(i, 1, 1)
       y_momentum(i, 1, 1) = density(i, 1, 1) * y_velocity(i, 1, 1)
       z_momentum(i, 1, 1) = density(i, 1, 1) * z_velocity(i, 1, 1)

       magf = Create_vector(x_magfield(i, 1, 1), y_magfield(i, 1, 1), z_magfield(i, 1, 1))
       velc = Create_vector(x_velocity(i, 1, 1), y_velocity(i, 1, 1), z_velocity(i, 1, 1))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure(i, 1, 1) = Calculate_total_pressure (pressure(i, 1, 1), magf_squrd)
       energy(i, 1, 1) = Calculate_energy_EOS (density(i, 1, 1), pressure(i, 1, 1), gamma(i, 1, 1), velc_squrd, magf_squrd)

    end do



  ! Set the boundary conditions.
  ! ----------------------------

    BOUNDTYPE = 1     ! see <TES_A>

    return
    
    

! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Set_initial_conditions_B
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Set Initial Conditions of Fluid ---------------------------------------------------------------------- [TES03] ---
! ----------------------------------------------------------------------------------------------------------------------------------
! --- Version C: Ryu and Jones' Shock Tube Test (4D) -------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Set_initial_conditions_C ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, x_zero
    real (PREC) :: XZERO
    
    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd



  ! Time and space.
  ! ---------------

    DELTAT   = 0.0001D0
    FULLTIME = 0.16D0

    COURANT    = 0.25D0
    MAX_DELTAT = 0.01D0

    DIMENSIONS = 1
    BOUNDARY   = 4

    XSIZE = 512
    YSIZE = 1
    ZSIZE = 1

    XDOMN = 1.0D0
    YDOMN = 1.0D0 * (real(YSIZE) / real(XSIZE))
    ZDOMN = 1.0D0 * (real(ZSIZE) / real(XSIZE))

    MHDUNITS = 2     ! see <TES_B>

    XZERO = 0.5D0



  ! Initialising fluid.
  ! -------------------

    BASE_FILE_NO = 0
    BASE_NTIME   = 0.00D0

    call Initialise_fluid ()     ! --> <FLU01> 



  ! Setting initial conditions.
  ! ---------------------------

    x_zero = ceiling(XZERO * real(XSIZE))

    do i = 1, x_zero

       ! Note: these values need setting. 

       density(i, 1, 1)  = 1.00D0
       pressure(i, 1, 1) = 1.00D0

       x_velocity(i, 1, 1) = 0.0D0
       y_velocity(i, 1, 1) = 0.0D0
       z_velocity(i, 1, 1) = 0.0D0

       gamma(i, 1, 1) = 5.0D0 / 3.0D0

       x_magfield(i, 1, 1) = 0.7D0
       y_magfield(i, 1, 1) = 0.0D0
       z_magfield(i, 1, 1) = 0.0D0


       ! Note: these values are set automatically.

       x_momentum(i, 1, 1) = density(i, 1, 1) * x_velocity(i, 1, 1)
       y_momentum(i, 1, 1) = density(i, 1, 1) * y_velocity(i, 1, 1)
       z_momentum(i, 1, 1) = density(i, 1, 1) * z_velocity(i, 1, 1)

       magf = Create_vector(x_magfield(i, 1, 1), y_magfield(i, 1, 1), z_magfield(i, 1, 1))
       velc = Create_vector(x_velocity(i, 1, 1), y_velocity(i, 1, 1), z_velocity(i, 1, 1))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure(i, 1, 1) = Calculate_total_pressure (pressure(i, 1, 1), magf_squrd)
       energy(i, 1, 1) = Calculate_energy_EOS (density(i, 1, 1), pressure(i, 1, 1), gamma(i, 1, 1), velc_squrd, magf_squrd)

    end do


    do i = x_zero + 1, XSIZE

       ! Note: these values need setting. 

       density(i, 1, 1)  = 0.3D0
       pressure(i, 1, 1) = 0.2D0

       x_velocity(i, 1, 1) = 0.0D0
       y_velocity(i, 1, 1) = 0.0D0
       z_velocity(i, 1, 1) = 1.0D0

       gamma(i, 1, 1) = 5.0D0 / 3.0D0

       x_magfield(i, 1, 1) = 0.7D0
       y_magfield(i, 1, 1) = 1.0D0
       z_magfield(i, 1, 1) = 0.0D0


       ! Note: these values are set automatically.

       x_momentum(i, 1, 1) = density(i, 1, 1) * x_velocity(i, 1, 1)
       y_momentum(i, 1, 1) = density(i, 1, 1) * y_velocity(i, 1, 1)
       z_momentum(i, 1, 1) = density(i, 1, 1) * z_velocity(i, 1, 1)

       magf = Create_vector(x_magfield(i, 1, 1), y_magfield(i, 1, 1), z_magfield(i, 1, 1))
       velc = Create_vector(x_velocity(i, 1, 1), y_velocity(i, 1, 1), z_velocity(i, 1, 1))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure(i, 1, 1) = Calculate_total_pressure (pressure(i, 1, 1), magf_squrd)
       energy(i, 1, 1) = Calculate_energy_EOS (density(i, 1, 1), pressure(i, 1, 1), gamma(i, 1, 1), velc_squrd, magf_squrd)

    end do



  ! Set the boundary conditions.
  ! ----------------------------

    BOUNDTYPE = 1     ! see <TES_A>

    return
    
    

! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Set_initial_conditions_C
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Set Initial Conditions of Fluid ---------------------------------------------------------------------- [TES05] ---
! ----------------------------------------------------------------------------------------------------------------------------------
! --- Version E: Circularly Polarised Alfven Waves ---------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Set_initial_conditions_E ()


  ! Declaration of local variables.
  ! -------------------------------
    
    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd

    real (PREC), dimension(:,:), allocatable :: v_prll, v_perp
    real (PREC), dimension(:,:), allocatable :: B_prll, B_perp

    real (PREC) :: x, y

    integer :: i, j



  ! Time and space.
  ! ---------------

    DELTAT   = 0.001D0                  
    FULLTIME = 5.0D0        

    COURANT    = 0.25D0
    MAX_DELTAT = 0.01D0    
  
    DIMENSIONS = 2
    BOUNDARY   = 4

    XSIZE = 32     
    YSIZE = 32   
    ZSIZE = 1    

    alpha = PI / 6.0D0   

    XDOMN = 1.0D0 / cos(alpha)
    YDOMN = 1.0D0 / sin(alpha)           
    ZDOMN = 1.0D0 * (real(ZSIZE) / real(XSIZE))

    MHDUNITS = 2     ! see <TES_B>



  ! Initialising fluid.
  ! -------------------

    BASE_FILE_NO = 0
    BASE_NTIME   = 0.00D0

    call Initialise_fluid ()     ! --> <FLU01>


  ! Allocating variables. 
  ! ---------------------

    allocate (v_prll(XSIZE, YSIZE))
    allocate (v_perp(XSIZE, YSIZE))
    allocate (B_prll(XSIZE, YSIZE))
    allocate (B_perp(XSIZE, YSIZE))



  ! Setting the velocity and magnetic field.
  ! ----------------------------------------

    do j = 1, YSIZE
       do i = 1, XSIZE

          x = (real(i) - 0.5D0) * DELTAX
          y = (real(j) - 0.5D0) * DELTAY

          v_prll(i, j) = 0.0D0
          v_perp(i, j) = 0.1 * sin(2.0D0 * PI * (x * cos(alpha) + y * sin(alpha)))
          B_prll(i, j) = 1.0D0
          B_perp(i, j) = v_perp(i, j)

       end do
    end do



  ! Setting initial conditions.
  ! ---------------------------

    do j = 1, YSIZE
       do i = 1, XSIZE

          x = (real(i) - 0.5D0) * DELTAX
          y = (real(j) - 0.5D0) * DELTAY


          ! Note: these values need setting. 

          density(i, j, 1)  = 1.0D0      
          pressure(i, j, 1) = 0.1D0  

          x_velocity(i, j, 1) = (v_prll(i, j) * cos(alpha)) - (v_perp(i, j) * sin(alpha))
          y_velocity(i, j, 1) = (v_prll(i, j) * sin(alpha)) + (v_perp(i, j) * cos(alpha))       
          z_velocity(i, j, 1) = 0.1 * cos(2.0D0 * PI * (x * cos(alpha) + y * sin(alpha)))           

          gamma(i, j, 1) = 5.0D0 / 3.0D0     

          x_magfield(i, j, 1) = ((B_prll(i, j) * cos(alpha)) - (B_perp(i, j) * sin(alpha)))
          y_magfield(i, j, 1) = ((B_prll(i, j) * sin(alpha)) + (B_perp(i, j) * cos(alpha))) 
          z_magfield(i, j, 1) = (0.1 * cos(2.0D0 * PI * (x * cos(alpha) + y * sin(alpha))))   

          ! Note: these values are set automatically.

          x_momentum(i, j, 1) = density(i, j, 1) * x_velocity(i, j, 1)
          y_momentum(i, j, 1) = density(i, j, 1) * y_velocity(i, j, 1)
          z_momentum(i, j, 1) = density(i, j, 1) * z_velocity(i, j, 1)

          magf = Create_vector(x_magfield(i, j, 1), y_magfield(i, j, 1), z_magfield(i, j, 1))
          velc = Create_vector(x_velocity(i, j, 1), y_velocity(i, j, 1), z_velocity(i, j, 1))

          magf_squrd = Dotproduct (magf, magf)
          velc_squrd = Dotproduct (velc, velc)

          pressure(i, j, 1) = Calculate_total_pressure (pressure(i, j, 1), magf_squrd)
          energy(i, j, 1) = Calculate_energy_EOS (density(i, j, 1), pressure(i, j, 1), gamma(i, j, 1), velc_squrd, magf_squrd)

       end do
    end do



  ! Set the boundary conditions.
  ! ----------------------------

    BOUNDTYPE = 2     ! see <TES_A>



  ! Deallocating variables.
  ! -----------------------

    deallocate (v_prll)
    deallocate (v_perp)
    deallocate (B_prll)
    deallocate (B_perp)

    return
    
    
! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Set_initial_conditions_E
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Set Initial Conditions of Fluid ---------------------------------------------------------------------- [TES06] ---
! ----------------------------------------------------------------------------------------------------------------------------------
! --- Version F: The Blast Problem -------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Set_initial_conditions_F ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, j

    real (PREC) :: x0, y0
    real (PREC) :: radius
    
    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd



  ! Time and space.
  ! ---------------

    DELTAT   = 0.00001D0
    FULLTIME = 0.2D0

    COURANT    = 0.25D0
    MAX_DELTAT = 0.01D0

    DIMENSIONS = 2
    BOUNDARY   = 4

    XSIZE = 200
    YSIZE = 300
    ZSIZE = 1

    XDOMN = 1.0D0
    YDOMN = 1.5D0
    ZDOMN = 1.0D0 * (real(ZSIZE) / real(XSIZE))

    MHDUNITS = 1     ! see <TES_B>



  ! Initialising fluid.
  ! -------------------

    BASE_FILE_NO = 0
    BASE_NTIME   = 0.00D0

    call Initialise_fluid ()     ! --> <FLU01>



  ! Setting initial conditions.
  ! ---------------------------

    x0 = 0.5D0 * real(XSIZE)
    y0 = 0.5D0 * real(YSIZE)


    do j = 1, YSIZE
       do i = 1, XSIZE

          ! Note: these values need setting. 

          density(i, j, 1)  = 1.0D0      
          pressure(i, j, 1) = 0.1D0  

          x_velocity(i, j, 1) = 0.00D0   
          y_velocity(i, j, 1) = 0.00D0     
          z_velocity(i, j, 1) = 0.00D0 

          gamma(i, j, 1) = 5.0D0 / 3.0D0     

          x_magfield(i, j, 1) = (1.0D0 / dsqrt(2.0D0)) * dsqrt(4.0D0 * PI)
          y_magfield(i, j, 1) = (1.0D0 / dsqrt(2.0D0)) * dsqrt(4.0D0 * PI)
          z_magfield(i, j, 1) = 0.0D0    

          radius = dsqrt( (((real(i)-0.5D0) - x0) * DELTAX)**2.0D0 + &
                          (((real(j)-0.5D0) - y0) * DELTAY)**2.0D0 )
 
          if (radius <= 0.1) pressure(i, j, 1) = 10.0D0


          ! Note: these values are set automatically.

          x_momentum(i, j, 1) = density(i, j, 1) * x_velocity(i, j, 1)
          y_momentum(i, j, 1) = density(i, j, 1) * y_velocity(i, j, 1)
          z_momentum(i, j, 1) = density(i, j, 1) * z_velocity(i, j, 1)

          magf = Create_vector(x_magfield(i, j, 1), y_magfield(i, j, 1), z_magfield(i, j, 1))
          velc = Create_vector(x_velocity(i, j, 1), y_velocity(i, j, 1), z_velocity(i, j, 1))

          magf_squrd = Dotproduct (magf, magf)
          velc_squrd = Dotproduct (velc, velc)

          pressure(i, j, 1) = Calculate_total_pressure (pressure(i, j, 1), magf_squrd)
          energy(i, j, 1) = Calculate_energy_EOS (density(i, j, 1), pressure(i, j, 1), gamma(i, j, 1), velc_squrd, magf_squrd)

       end do
    end do



  ! Set the boundary conditions.
  ! ----------------------------

    BOUNDTYPE = 2     ! see <TES_A>

    return
    
    

! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Set_initial_conditions_F
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Set Initial Conditions of Fluid ---------------------------------------------------------------------- [TES07] ---
! ----------------------------------------------------------------------------------------------------------------------------------
! --- Version G: The Orszag-Tang Vortex --------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Set_initial_conditions_G ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, j
    real (PREC) :: x, y, B0
    
    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd



  ! Time and space.
  ! ---------------

    DELTAT   = 0.00001D0
    FULLTIME = 0.25D0

    COURANT    = 0.25D0
    MAX_DELTAT = 0.01D0
 
    DIMENSIONS = 2
    BOUNDARY   = 4

    XSIZE = 256
    YSIZE = 256
    ZSIZE = 1

    XDOMN = 1.0D0
    YDOMN = 1.0D0
    ZDOMN = 1.0D0 * (real(ZSIZE) / real(XSIZE))

    MHDUNITS = 1     ! see <TES_B> 



  ! Initialising fluid.
  ! -------------------

    BASE_FILE_NO = 0
    BASE_NTIME   = 0.00D0

    call Initialise_fluid ()     ! --> <FLU01>



  ! Setting initial conditions.
  ! ---------------------------

    do j = 1, YSIZE
       do i = 1, XSIZE

          x = (real(i)-0.5D0) * DELTAX
          y = (real(j)-0.5D0) * DELTAY

          B0 = 1.0D0 / sqrt(4.0D0 * PI)


          ! Note: these values need setting. 

          density(i, j, 1)  = 25.0D0 / (36.0D0 * PI)
          pressure(i, j, 1) = 5.0D0 / (12.0D0 * PI)

          x_velocity(i, j, 1) = -sin(2.0D0 * PI * y)
          y_velocity(i, j, 1) =  sin(2.0D0 * PI * x)   
          z_velocity(i, j, 1) = 0.00D0 

          gamma(i, j, 1) = 5.0D0 / 3.0D0      

          x_magfield(i, j, 1) = -B0 * sin(2.0D0 * PI * y) * sqrt(4.0D0 * PI)  
          y_magfield(i, j, 1) =  B0 * sin(4.0D0 * PI * x) * sqrt(4.0D0 * PI)
          z_magfield(i, j, 1) = 0.0D0   


          ! Note: these values are set automatically.

          x_momentum(i, j, 1) = density(i, j, 1) * x_velocity(i, j, 1)
          y_momentum(i, j, 1) = density(i, j, 1) * y_velocity(i, j, 1)
          z_momentum(i, j, 1) = density(i, j, 1) * z_velocity(i, j, 1)

          magf = Create_vector(x_magfield(i, j, 1), y_magfield(i, j, 1), z_magfield(i, j, 1))
          velc = Create_vector(x_velocity(i, j, 1), y_velocity(i, j, 1), z_velocity(i, j, 1))

          magf_squrd = Dotproduct (magf, magf)
          velc_squrd = Dotproduct (velc, velc)

          pressure(i, j, 1) = Calculate_total_pressure (pressure(i, j, 1), magf_squrd)
          energy(i, j, 1) = Calculate_energy_EOS (density(i, j, 1), pressure(i, j, 1), gamma(i, j, 1), velc_squrd, magf_squrd)

       end do
    end do



  ! Set the boundary conditions.
  ! ----------------------------

    BOUNDTYPE = 2     ! see <TES_A>

    return
    
    

! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Set_initial_conditions_G
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Set Initial Conditions of Fluid ---------------------------------------------------------------------- [TES08] ---
! ----------------------------------------------------------------------------------------------------------------------------------
! --- Version H: The MHD Rotor Problem ---------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Set_initial_conditions_H ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, j

    real (PREC) :: x0, y0
    real (PREC) :: radius
    real (PREC) :: r0, r1, v0
    real (PREC) :: x, y, f
    
    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd



  ! Time and space.
  ! ---------------

    DELTAT   = 0.00001D0
    FULLTIME = 0.15D0  

    COURANT    = 0.25D0
    MAX_DELTAT = 0.01D0    
  
    DIMENSIONS = 2
    BOUNDARY   = 4

    XSIZE = 256
    YSIZE = 256
    ZSIZE = 1

    XDOMN = 1.0D0
    YDOMN = 1.0D0
    ZDOMN = 1.0D0 * (real(ZSIZE) / real(XSIZE))

    MHDUNITS = 1     ! see <TES_B>



  ! Initialising fluid.
  ! -------------------

    BASE_FILE_NO = 0
    BASE_NTIME   = 0.00D0

    call Initialise_fluid ()     ! --> <FLU01>



  ! Setting initial conditions.
  ! ---------------------------

    r0 = 0.1D0
    r1 = 0.115D0

    v0 = 2.0D0

    x0 = 0.5D0 * XDOMN
    y0 = 0.5D0 * YDOMN


    do j = 1, YSIZE
       do i = 1, XSIZE

          x = (real(i)-0.5D0) * DELTAX
          y = (real(j)-0.5D0) * DELTAY


          ! Note: these values need setting. 

          density(i, j, 1)  =  1.0D0          
          pressure(i, j, 1) =  1.0D0        

          x_velocity(i, j, 1) = 0.0D0         
          y_velocity(i, j, 1) = 0.0D0           
          z_velocity(i, j, 1) = 0.0D0         

          gamma(i, j, 1) = 1.4D0                     

          x_magfield(i, j, 1) = 5.0D0                                 
          y_magfield(i, j, 1) = 0.0D0                                          
          z_magfield(i, j, 1) = 0.0D0              

          radius = dsqrt( (x - x0)**2.0D0 + &
                          (y - y0)**2.0D0 )

          if (radius <= r0) then          ! the rotor

             density(i, j, 1) = 10.0D0

             x_velocity(i, j, 1) = -v0 * (y - 0.5D0) / r0
             y_velocity(i, j, 1) =  v0 * (x - 0.5D0) / r0

          else if (radius <= r1) then     ! the taper

             f = (r1 - radius)/(r1 - r0)

             density(i, j, 1) = 1.0D0 + (9.0D0 * f)

             x_velocity(i, j, 1) = -f * v0 * (y - 0.5D0) / radius
             y_velocity(i, j, 1) =  f * v0 * (x - 0.5D0) / radius

          end if


          ! Note: these values are set automatically.

          x_momentum(i, j, 1) = density(i, j, 1) * x_velocity(i, j, 1)
          y_momentum(i, j, 1) = density(i, j, 1) * y_velocity(i, j, 1)
          z_momentum(i, j, 1) = density(i, j, 1) * z_velocity(i, j, 1)

          magf = Create_vector(x_magfield(i, j, 1), y_magfield(i, j, 1), z_magfield(i, j, 1))
          velc = Create_vector(x_velocity(i, j, 1), y_velocity(i, j, 1), z_velocity(i, j, 1))

          magf_squrd = Dotproduct (magf, magf)
          velc_squrd = Dotproduct (velc, velc)

          pressure(i, j, 1) = Calculate_total_pressure (pressure(i, j, 1), magf_squrd)
          energy(i, j, 1) = Calculate_energy_EOS (density(i, j, 1), pressure(i, j, 1), gamma(i, j, 1), velc_squrd, magf_squrd)

       end do
    end do



  ! Set the boundary conditions.
  ! ----------------------------

    BOUNDTYPE = 2     ! see <TES_A>

    return
    
    

! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Set_initial_conditions_H
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Set Initial Conditions of Fluid ---------------------------------------------------------------------- [TES09] ---
! ----------------------------------------------------------------------------------------------------------------------------------
! --- Version I: The Kelvin-Helmholtz Instability ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
 
  subroutine Set_initial_conditions_I


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, j, qtop, qbot

    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd
 


  ! Time and space.
  ! ---------------

    DELTAT   = 0.00001D0
    FULLTIME = 1.0D0

    COURANT    = 0.25D0
    MAX_DELTAT = 0.01D0

    DIMENSIONS = 2
    BOUNDARY   = 4

    XSIZE = 240
    YSIZE = 240
    ZSIZE = 1

    XDOMN = 1.0D0
    YDOMN = 1.0D0
    ZDOMN = 1.0D0 * (real(ZSIZE) / real(XSIZE))

    MHDUNITS = 1     ! see <TES_B>             



  ! Initialising fluid.
  ! -------------------

    BASE_FILE_NO = 0
    BASE_NTIME   = 0.00D0

    call Initialise_fluid ()     ! --> <FLU01>



  ! Setting initial conditions.
  ! ---------------------------

    qtop = int(0.75D0 * real(YSIZE))
    qbot = int(0.25D0 * real(YSIZE))


    do j = 1, YSIZE
       do i = 1, XSIZE

          ! Note: these values need setting. 

          density(i, j, 1)  = 2.0D0          
          pressure(i, j, 1) = 2.5D0 

          x_velocity(i, j, 1) = 0.5D0 + (0.01D0 * rand())
          y_velocity(i, j, 1) = 0.0D0 + (0.01D0 * rand())       
          z_velocity(i, j, 1) = 0.0D0 

          gamma(i, j, 1) = 1.4D0              

          x_magfield(i, j, 1) = 0.5D0  
          y_magfield(i, j, 1) = 0.0D0   
          z_magfield(i, j, 1) = 0.0D0   

          if ((j > qtop) .or. (j < qbot)) then

             x_velocity(i, j, 1) = -0.5D0 + (0.01D0 * rand())              
             density(i, j, 1)    =  1.0D0

          end if


          ! Note: these values are set automatically.

          x_momentum(i, j, 1) = density(i, j, 1) * x_velocity(i, j, 1)
          y_momentum(i, j, 1) = density(i, j, 1) * y_velocity(i, j, 1)
          z_momentum(i, j, 1) = density(i, j, 1) * z_velocity(i, j, 1)

          magf = Create_vector(x_magfield(i, j, 1), y_magfield(i, j, 1), z_magfield(i, j, 1))
          velc = Create_vector(x_velocity(i, j, 1), y_velocity(i, j, 1), z_velocity(i, j, 1))

          magf_squrd = Dotproduct (magf, magf)
          velc_squrd = Dotproduct (velc, velc)

          pressure(i, j, 1) = Calculate_total_pressure (pressure(i, j, 1), magf_squrd)
          energy(i, j, 1) = Calculate_energy_EOS (density(i, j, 1), pressure(i, j, 1), gamma(i, j, 1), velc_squrd, magf_squrd) 

       end do
    end do



  ! Set the boundary conditions.
  ! ----------------------------

    BOUNDTYPE = 2     ! see <TES_A>

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Set_initial_conditions_I
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Load a Restart File  --------------------------------------------------------------------------------- [TES10] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Load_restart_files()


  ! Declaration of local variables.
  ! -------------------------------

    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd

    integer :: i, j, k



  ! Reading the restart files.
  ! --------------------------

    open (50, file = "output/restart1.bin", form = 'unformatted')
    read (50) BASE_FILE_NO, BASE_NTIME, DELTAT, FULLTIME, COURANT, MAX_DELTAT, XSIZE, &
         YSIZE, ZSIZE, DIMENSIONS, BOUNDARY, XDOMN, YDOMN, ZDOMN, MHDUNITS, BOUNDTYPE
    close (50)

    call Initialise_fluid ()     ! --> <FLU01>

    open (60, file = "output/restart2.bin", form = 'unformatted')
    read (60) density, pressure, x_velocity, y_velocity, z_velocity, gamma, &
         x_magfield, y_magfield, z_magfield
    close (60)



  ! Calculating the secondary values.
  ! ---------------------------------

    do k = 1, ZSIZE
       do j = 1, YSIZE
          do i = 1, XSIZE

             x_momentum(i, j, k) = density(i, j, k) * x_velocity(i, j, k)
             y_momentum(i, j, k) = density(i, j, k) * y_velocity(i, j, k)
             z_momentum(i, j, k) = density(i, j, k) * z_velocity(i, j, k)
             
             magf = Create_vector(x_magfield(i, j, k), y_magfield(i, j, k), z_magfield(i, j, k))
             velc = Create_vector(x_velocity(i, j, k), y_velocity(i, j, k), z_velocity(i, j, k))
             
             magf_squrd = Dotproduct (magf, magf)
             velc_squrd = Dotproduct (velc, velc)
             
             energy(i, j, k) = Calculate_energy_EOS (density(i, j, k), pressure(i, j, k), gamma(i, j, k), velc_squrd, magf_squrd)

          end do
       end do
    end do

    return
    


! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Load_restart_files
! ----------------------------------------------------------------------------------------------------------------------------------



end module tests
