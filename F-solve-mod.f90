! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
! 
!     F: The Riemann Solver Module
!     ============================
!     -- Main variables: the inter-cell fluxes for the state variables.
!     -- Main actions: determine these inter-cell fluxes.
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



module solve

  use setup
  use recon
  implicit none


  ! The inter-cell fluxes.
  ! ----------------------

  real (PREC), dimension(:), allocatable :: density_flux
  real (PREC), dimension(:), allocatable :: x_momentum_flux
  real (PREC), dimension(:), allocatable :: y_momentum_flux
  real (PREC), dimension(:), allocatable :: z_momentum_flux
  real (PREC), dimension(:), allocatable :: energy_flux

  real (PREC), dimension(:), allocatable :: x_magfield_flux
  real (PREC), dimension(:), allocatable :: y_magfield_flux
  real (PREC), dimension(:), allocatable :: z_magfield_flux



  ! Working variables.
  ! ------------------

  real (PREC), dimension(:,:), allocatable :: flux_L, flux_R
  
  real (PREC), dimension(:), allocatable :: SL, SR, SM, SLS, SRS
  
  real (PREC), dimension(:), allocatable :: pressure_S1, x_velocity_S1, x_magfield_S1
  
  real (PREC), dimension(:), allocatable :: density_SL1,    energy_SL1
  real (PREC), dimension(:), allocatable :: y_velocity_SL1, y_magfield_SL1
  real (PREC), dimension(:), allocatable :: z_velocity_SL1, z_magfield_SL1
  
  real (PREC), dimension(:), allocatable :: density_SR1,    energy_SR1
  real (PREC), dimension(:), allocatable :: y_velocity_SR1, y_magfield_SR1
  real (PREC), dimension(:), allocatable :: z_velocity_SR1, z_magfield_SR1
  
  real (PREC), dimension(:), allocatable :: pressure_S2
  real (PREC), dimension(:), allocatable :: x_velocity_S2,  x_magfield_S2
  real (PREC), dimension(:), allocatable :: y_velocity_S2,  y_magfield_S2
  real (PREC), dimension(:), allocatable :: z_velocity_S2,  z_magfield_S2
  
  real (PREC), dimension(:), allocatable :: density_SL2,    energy_SL2
  real (PREC), dimension(:), allocatable :: density_SR2,    energy_SR2
  
  integer, dimension(:), allocatable :: flux_flag
  
contains



! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Initialise the Solver -------------------------------------------------------------------------------- [SOL01] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Initialise_solver ()


  ! Allocating the inter-cell fluxes.
  ! ---------------------------------

    allocate (density_flux   (MFULL))
    allocate (x_momentum_flux(MFULL))
    allocate (y_momentum_flux(MFULL))
    allocate (z_momentum_flux(MFULL))
    allocate (energy_flux    (MFULL))

    allocate (x_magfield_flux(MFULL))
    allocate (y_magfield_flux(MFULL))
    allocate (z_magfield_flux(MFULL))



  ! Allocating the working variables.
  ! ---------------------------------

    allocate (flux_L(MFULL, 8)); allocate (flux_R(MFULL, 8))

    allocate (SL(MFULL));  allocate (SR(MFULL));  allocate (SM(MFULL))
    allocate (SLS(MFULL)); allocate (SRS(MFULL))

    allocate (pressure_S1(MFULL)); allocate (x_velocity_S1(MFULL)); allocate (x_magfield_S1(MFULL))

    allocate (density_SL1(MFULL));    allocate (energy_SL1(MFULL))
    allocate (y_velocity_SL1(MFULL)); allocate (y_magfield_SL1(MFULL))
    allocate (z_velocity_SL1(MFULL)); allocate (z_magfield_SL1(MFULL))

    allocate (density_SR1(MFULL));    allocate (energy_SR1(MFULL))
    allocate (y_velocity_SR1(MFULL)); allocate (y_magfield_SR1(MFULL))
    allocate (z_velocity_SR1(MFULL)); allocate (z_magfield_SR1(MFULL))

    allocate (pressure_S2(MFULL))
    allocate (x_velocity_S2(MFULL));  allocate (x_magfield_S2(MFULL))
    allocate (y_velocity_S2(MFULL));  allocate (y_magfield_S2(MFULL))
    allocate (z_velocity_S2(MFULL));  allocate (z_magfield_S2(MFULL))

    allocate (density_SL2(MFULL));    allocate (energy_SL2(MFULL))
    allocate (density_SR2(MFULL));    allocate (energy_SR2(MFULL))

    allocate (flux_flag(MFULL))

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Initialise_solver
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Calculating the Inter-cell Fluxes (HLLD Edition) ----------------------------------------------------- [SOL02] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Calculate_intercell_fluxes_HLLD ()

    ! This subroutine is based on the algorithm of Miyoshi and Kusano, J. Comput. Phys., 208, 315 (2005)


  ! Declaration of local variables.
  ! -------------------------------

    real (PREC), dimension(3) :: velc, magf

    real (PREC) :: gas_pressure, a2, ca2, cax2, s1, s2, cfL, cfR

    real (PREC) :: cL, cR, bL, bR, qL, qR, wL, wR, magf_squrd, scratch, roe_temp

    real (PREC), dimension(3) :: uX, BX, uS1, BS1, uS2, BS2

    real (PREC) ::  x_momentum_SX1, y_momentum_SX1
    real (PREC) ::  z_momentum_SX1, denRS, denLS

    real (PREC), parameter :: epsilon = 1.0D-14

    integer :: m



! --- PART ONE: Preliminary Calculations -------------------------------------------------------------------------------------------


    ! Obtaining the left and right flux functions.
    ! --------------------------------------------

    do m = BOUNDARY, rowsize + BOUNDARY

       ! The left states.

       velc = Create_vector (x_velocity_L(m), y_velocity_L(m), z_velocity_L(m))
       magf = Create_vector (x_magfield_L(m), y_magfield_L(m), z_magfield_L(m))
       magf_squrd = Dotproduct (magf, magf)

       flux_L(m, 1) = density_L(m) * x_velocity_L(m)
       flux_L(m, 2) = flux_L(m, 1) * x_velocity_L(m) + pressure_L(m) - MHDF(2) * (x_magfield_L(m)**2.0D0)
       flux_L(m, 3) = flux_L(m, 1) * y_velocity_L(m) - MHDF(2) * (x_magfield_L(m) * y_magfield_L(m))
       flux_L(m, 4) = flux_L(m, 1) * z_velocity_L(m) - MHDF(2) * (x_magfield_L(m) * z_magfield_L(m))
       flux_L(m, 5) = (energy_L(m) + pressure_L(m)) * x_velocity_L(m) - MHDF(2) * (x_magfield_L(m) * Dotproduct (magf, velc))

       flux_L(m, 6) = 0.00D0
       flux_L(m, 7) = x_velocity_L(m) * y_magfield_L(m) - y_velocity_L(m) * x_magfield_L(m)
       flux_L(m, 8) = x_velocity_L(m) * z_magfield_L(m) - z_velocity_L(m) * x_magfield_L(m) 


       ! The right states.

       velc = Create_vector (x_velocity_R(m), y_velocity_R(m), z_velocity_R(m))
       magf = Create_vector (x_magfield_R(m), y_magfield_R(m), z_magfield_R(m))

       flux_R(m, 1) = density_R(m) * x_velocity_R(m)
       flux_R(m, 2) = flux_R(m, 1) * x_velocity_R(m) + pressure_R(m) - MHDF(2) * (x_magfield_R(m)**2.0D0)
       flux_R(m, 3) = flux_R(m, 1) * y_velocity_R(m) - MHDF(2) * (x_magfield_R(m) * y_magfield_R(m))
       flux_R(m, 4) = flux_R(m, 1) * z_velocity_R(m) - MHDF(2) * (x_magfield_R(m) * z_magfield_R(m))
       flux_R(m, 5) = (energy_R(m) + pressure_R(m)) * x_velocity_R(m) - MHDF(2) * (x_magfield_R(m) * Dotproduct (magf, velc))

       flux_R(m, 6) = 0.00D0
       flux_R(m, 7) = x_velocity_R(m) * y_magfield_R(m) - y_velocity_R(m) * x_magfield_R(m)
       flux_R(m, 8) = x_velocity_R(m) * z_magfield_R(m) - z_velocity_R(m) * x_magfield_R(m)

    end do



  ! Estimating the wavespeeds: SL and SR.
  ! -------------------------------------

    do m =  BOUNDARY, rowsize + BOUNDARY

       ! Wavespeed estimate of Miyoshi and Kusano, J. Comput. Phys., 208, 315 (2005) -- eq. (67)

       if (WAVESPEED_TYPE == 'M') then

          magf = Create_vector (x_magfield_L(m), y_magfield_L(m), z_magfield_L(m))
          magf_squrd = Dotproduct (magf, magf)

          gas_pressure = Calculate_gas_pressure (pressure_L(m), magf_squrd)

          a2 = (gamma_L(m) * gas_pressure) / density_L(m)
          
          ca2  = (MHDF(2) * magf_squrd) / density_L(m)
          cax2 = (MHDF(2) * x_magfield_L(m)**2.0D0) / density_L(m)

          s1 = a2 + ca2
          s2 = 4.0D0 * a2 * cax2

          cfL = dsqrt(0.5D0 * (s1 + dsqrt(s1**2.0D0 - s2)))


          magf = Create_vector (x_magfield_R(m), y_magfield_R(m), z_magfield_R(m))
          magf_squrd = Dotproduct (magf, magf)

          gas_pressure = Calculate_gas_pressure (pressure_R(m), magf_squrd)

          a2 = (gamma_R(m) * gas_pressure) / density_R(m)
          
          ca2  = (MHDF(2) * magf_squrd) / density_R(m)
          cax2 = (MHDF(2) * x_magfield_R(m)**2.0D0) / density_R(m)

          s1 = a2 + ca2
          s2 = 4.0D0 * a2 * cax2

          cfR = dsqrt(0.5D0 * (s1 + dsqrt(s1**2.0D0 - s2)))

          
          SL(m) = min(x_velocity_L(m), x_velocity_R(m)) - max(cfL, cfR)
          SR(m) = max(x_velocity_L(m), x_velocity_R(m)) + max(cfL, cfR)

          
          ! DEBUG MODE: Checking the CFL stability condition

          if (DEBUG_MODE) then

             scratch = (max(abs(SL(m)), abs(SR(m))) * dt(m)) / dx(m)

             if (scratch > COURANT) CFL_violations = CFL_violations + 1

          end if

       end if


       ! Wavespeed estimate of Janhunen, J. Comput. Phys., 160, 649 (2000) -- eq. (27)

       if (WAVESPEED_TYPE == 'J') then

          a2 = gamma_1D(m) * max(pressure_L(m), pressure_R(m)) / min(density_L(m), density_R(m))

          magf = Create_vector (x_magfield_L(m), y_magfield_L(m), z_magfield_L(m))
          s1 = Dotproduct (magf, magf) * MHDF(2)

          magf = Create_vector (x_magfield_R(m), y_magfield_R(m), z_magfield_R(m))
          s2 = Dotproduct (magf, magf) * MHDF(2)

          ca2  = max(s1, s2) / min(density_L(m), density_R(m))
          cax2 = max(x_magfield_L(m)**2.0D0, x_magfield_R(m)**2.0D0) * MHDF(2) / min(density_L(m), density_R(m))


          s1 = (a2 - ca2)**2.0D0 + (4.0D0 * a2 * (ca2 - cax2))
          s2 = 0.5D0 * (a2 + ca2 + sqrt(s1))

          SL(m) = min(x_velocity_L(m), x_velocity_R(m)) - sqrt(s2)
          SR(m) = max(x_velocity_L(m), x_velocity_R(m)) + sqrt(s2)


          ! DEBUG MODE: Checking the CFL stability condition

          if (DEBUG_MODE) then

             scratch = (max(abs(SL(m)), abs(SR(m))) * dt(m)) / dx(m)

             if (scratch > COURANT) CFL_violations = CFL_violations + 1

          end if
       
       end if
          
    end do


    
  ! Estimating the wavespeeds: SM.
  ! ------------------------------

    do m = BOUNDARY, rowsize + BOUNDARY

       qL = (density_L(m) * x_velocity_L(m) * (SL(m) - x_velocity_L(m))) - pressure_L(m)
       qR = (density_R(m) * x_velocity_R(m) * (SR(m) - x_velocity_R(m))) - pressure_R(m)

       wL = density_L(m) * (SL(m) - x_velocity_L(m))
       wR = density_R(m) * (SR(m) - x_velocity_R(m))

       SM(m) = (qR - qL) / (wR - wL)

    end do




! --- PART TWO: Finding the Intermediate States ------------------------------------------------------------------------------------


  ! Finding U*.
  ! -----------

    do m = BOUNDARY, rowsize + BOUNDARY

       ! Calculating the pressure

       scratch = (density_L(m) * density_R(m)) * (SR(m) - x_velocity_R(m)) 
       scratch = scratch * (SL(m) - x_velocity_L(m)) * (x_velocity_R(m) - x_velocity_L(m))

       pressure_S1(m) = scratch        + (SR(m) - x_velocity_R(m)) * density_R(m) * pressure_L(m)
       pressure_S1(m) = pressure_S1(m) - (SL(m) - x_velocity_L(m)) * density_L(m) * pressure_R(m)

       scratch = ((SR(m) - x_velocity_R(m)) * density_R(m)) - ((SL(m) - x_velocity_L(m)) * density_L(m))
       pressure_S1(m) = pressure_S1(m) / scratch


       ! The x-velocity and x-magnetic-field

       x_velocity_S1(m) = SM(m)
       x_magfield_S1(m) = x_magfield_L(m)

    end do



  ! Finding UL*.
  ! ------------

    do m = BOUNDARY, rowsize + BOUNDARY

       density_SL1(m) = density_L(m) * (SL(m) - x_velocity_L(m)) / (SL(m) - SM(m))


       scratch = (density_L(m) * (SL(m) - x_velocity_L(m)) * (SL(m) - SM(m))) - (MHDF(2) * x_magfield_S1(m)**2.0D0)

       if (abs(scratch) < epsilon) then

          y_velocity_SL1(m) = y_velocity_L(m)
          z_velocity_SL1(m) = z_velocity_L(m)
          y_magfield_SR1(m) = 0.0D0
          z_magfield_SR1(m) = 0.0D0

       else

          y_velocity_SL1(m) = y_velocity_L(m) - ((MHDF(2)*x_magfield_S1(m)*y_magfield_L(m)*(SM(m) - x_velocity_L(m))) / scratch)
          z_velocity_SL1(m) = z_velocity_L(m) - ((MHDF(2)*x_magfield_S1(m)*z_magfield_L(m)*(SM(m) - x_velocity_L(m))) / scratch)
          
          y_magfield_SL1(m) = y_magfield_L(m)*(density_L(m) * (SL(m) - x_velocity_L(m))**2.0D0 - MHDF(2)*x_magfield_S1(m)**2.0D0)
          z_magfield_SL1(m) = z_magfield_L(m)*(density_L(m) * (SL(m) - x_velocity_L(m))**2.0D0 - MHDF(2)*x_magfield_S1(m)**2.0D0)
          y_magfield_SL1(m) = y_magfield_SL1(m) / scratch
          z_magfield_SL1(m) = z_magfield_SL1(m) / scratch
          
       end if


       uX  = Create_vector (x_velocity_L(m),  y_velocity_L(m),   z_velocity_L(m))
       BX  = Create_vector (x_magfield_L(m),  y_magfield_L(m),   z_magfield_L(m))
       uS1 = Create_vector (x_velocity_S1(m), y_velocity_SL1(m), z_velocity_SL1(m))
       BS1 = Create_vector (x_magfield_S1(m), y_magfield_SL1(m), z_magfield_SL1(m))

       energy_SL1(m) = ((SL(m) - x_velocity_L(m)) * energy_L(m)) - (pressure_L(m) * x_velocity_L(m)) + (pressure_S1(m) * SM(m))
       energy_SL1(m) = (energy_SL1(m) + MHDF(2) * x_magfield_L(m) * (Dotproduct(uX, BX) - Dotproduct(uS1, BS1))) / (SL(m) - SM(m))

    end do



  ! Finding UR*.
  ! ------------

    do m = BOUNDARY, rowsize + BOUNDARY

       density_SR1(m) = density_R(m) * (SR(m) - x_velocity_R(m)) / (SR(m) - SM(m))


       scratch = (density_R(m) * (SR(m) - x_velocity_R(m)) * (SR(m) - SM(m))) - (MHDF(2) * x_magfield_S1(m)**2.0D0)

       if (abs(scratch) < epsilon) then

          y_velocity_SR1(m) = y_velocity_R(m)
          z_velocity_SR1(m) = z_velocity_R(m)
          y_magfield_SR1(m) = 0.0D0
          z_magfield_SR1(m) = 0.0D0

       else

          y_velocity_SR1(m) = y_velocity_R(m) - ((MHDF(2)*x_magfield_S1(m)*y_magfield_R(m) * (SM(m) - x_velocity_R(m))) / scratch)
          z_velocity_SR1(m) = z_velocity_R(m) - ((MHDF(2)*x_magfield_S1(m)*z_magfield_R(m) * (SM(m) - x_velocity_R(m))) / scratch)
          
          y_magfield_SR1(m) = y_magfield_R(m) * (density_R(m)*(SR(m) - x_velocity_R(m))**2.0D0 - MHDF(2)*x_magfield_S1(m)**2.0D0)
          z_magfield_SR1(m) = z_magfield_R(m) * (density_R(m)*(SR(m) - x_velocity_R(m))**2.0D0 - MHDF(2)*x_magfield_S1(m)**2.0D0)
          y_magfield_SR1(m) = y_magfield_SR1(m) / scratch
          z_magfield_SR1(m) = z_magfield_SR1(m) / scratch
       
       end if


       uX  = Create_vector (x_velocity_R(m),  y_velocity_R(m),   z_velocity_R(m))
       BX  = Create_vector (x_magfield_R(m),  y_magfield_R(m),   z_magfield_R(m))
       uS1 = Create_vector (x_velocity_S1(m), y_velocity_SR1(m), z_velocity_SR1(m))
       BS1 = Create_vector (x_magfield_S1(m), y_magfield_SR1(m), z_magfield_SR1(m))

       energy_SR1(m) = ((SR(m) - x_velocity_R(m)) * energy_R(m)) - (pressure_R(m) * x_velocity_R(m)) + (pressure_S1(m) * SM(m))
       energy_SR1(m) = (energy_SR1(m) + MHDF(2) * x_magfield_R(m) * (Dotproduct(uX, BX) - Dotproduct(uS1, BS1))) / (SR(m) - SM(m))

    end do




  ! Estimating the wavespeeds: SL* and SR*.
  ! ---------------------------------------

    do m = BOUNDARY, rowsize + BOUNDARY

       SLS(m) = SM(m) - abs(x_magfield_S1(m)) / (dsqrt(density_SL1(m)) * MHDF(-1))
       SRS(m) = SM(m) + abs(x_magfield_S1(m)) / (dsqrt(density_SR1(m)) * MHDF(-1))

    end do



  ! Finding U**.
  ! ------------

    pressure_S2   = pressure_S1
    x_velocity_S2 = x_velocity_S1
    x_magfield_S2 = x_magfield_S1

    
    do m = BOUNDARY, rowsize + BOUNDARY

       scratch = MHDF(1) * (dsqrt(density_SL1(m)) + dsqrt(density_SR1(m)))

       roe_temp         = roeavg(MHDF(2) * density_SL1(m), MHDF(2) * density_SR1(m), y_velocity_SL1(m), y_velocity_SR1(m))
       y_velocity_S2(m) = MHDF(2) * (y_magfield_SR1(m) - y_magfield_SL1(m)) * sign(1.0D0, x_magfield_S1(m))
       y_velocity_S2(m) = roe_temp + (y_velocity_S2(m) / scratch)

       roe_temp         = roeavg(MHDF(2) * density_SL1(m), MHDF(2) * density_SR1(m), z_velocity_SL1(m), z_velocity_SR1(m))
       z_velocity_S2(m) = MHDF(2) * (z_magfield_SR1(m) - z_magfield_SL1(m)) * sign(1.0D0, x_magfield_S1(m))
       z_velocity_S2(m) = roe_temp + (z_velocity_S2(m) / scratch)


       scratch = MHDF(-1) * (dsqrt(density_SL1(m)) + dsqrt(density_SR1(m)))

       roe_temp         = roeavg(MHDF(-2) * density_SL1(m), MHDF(-2) * density_SR1(m), y_magfield_SR1(m), y_magfield_SL1(m))
       y_magfield_S2(m) = MHDF(-2) * dsqrt(density_SL1(m) * density_SR1(m)) * (y_velocity_SR1(m) - y_velocity_SL1(m))
       y_magfield_S2(m) = y_magfield_S2(m) * sign(1.0D0, x_magfield_S1(m))
       y_magfield_S2(m) = roe_temp + (y_magfield_S2(m) / scratch)

       roe_temp         = roeavg(MHDF(-2) * density_SL1(m), MHDF(-2) * density_SR1(m), z_magfield_SR1(m), z_magfield_SL1(m))
       z_magfield_S2(m) = MHDF(-2) * dsqrt(density_SL1(m) * density_SR1(m)) * (z_velocity_SR1(m) - z_velocity_SL1(m))
       z_magfield_S2(m) = z_magfield_S2(m) * sign(1.0D0, x_magfield_S1(m))
       z_magfield_S2(m) = roe_temp + (z_magfield_S2(m) / scratch)

    end do



  ! Finding UL**.
  ! -------------

    do m = BOUNDARY, rowsize + BOUNDARY

       density_SL2(m) = density_SL1(m)

       uS1 = Create_vector (x_velocity_S1(m), y_velocity_SL1(m), z_velocity_SL1(m))
       BS1 = Create_vector (x_magfield_S1(m), y_magfield_SL1(m), z_magfield_SL1(m))
       uS2 = Create_vector (x_velocity_S2(m), y_velocity_S2(m),  z_velocity_S2(m))
       BS2 = Create_vector (x_magfield_S2(m), y_magfield_S2(m),  z_magfield_S2(m))

       energy_SL2(m) = dsqrt(density_SL1(m)) * MHDF(1) * (Dotproduct (uS1, BS1) - Dotproduct(uS2, BS2)) 
       energy_SL2(m) = energy_SL2(m) * sign(1.0D0, x_magfield_S1(m))
       energy_SL2(m) = energy_SL1(m) - energy_SL2(m)

    end do



  ! Finding UR**.
  ! -------------

    do m = BOUNDARY, rowsize + BOUNDARY

       density_SR2(m) = density_SR1(m)

       uS1 = Create_vector (x_velocity_S1(m), y_velocity_SR1(m), z_velocity_SR1(m))
       BS1 = Create_vector (x_magfield_S1(m), y_magfield_SR1(m), z_magfield_SR1(m))
       uS2 = Create_vector (x_velocity_S2(m), y_velocity_S2(m),  z_velocity_S2(m))
       BS2 = Create_vector (x_magfield_S2(m), y_magfield_S2(m),  z_magfield_S2(m))

       energy_SR2(m) = dsqrt(density_SR1(m)) * MHDF(1) * (Dotproduct (uS1, BS1) - Dotproduct(uS2, BS2)) 
       energy_SR2(m) = energy_SR2(m) * sign(1.0D0, x_magfield_S1(m))
       energy_SR2(m) = energy_SR1(m) + energy_SR2(m)

    end do



    
! --- PART THREE: Finding the Inter-cell Fluxes ------------------------------------------------------------------------------------


  ! Determining correct flux.
  ! -------------------------

    do m = BOUNDARY, rowsize + BOUNDARY

       if (SL(m) > 0.00D0) flux_flag(m) = 30

       if ((SL(m) <= 0.00D0)  .and. (0.00D0 <= SLS(m))) flux_flag(m) = 20

       if ((SLS(m) <= 0.00D0) .and. (0.00D0 <= SM(m)))  flux_flag(m) = 10


       if ((SM(m) <= 0.00D0)  .and. (0.00D0 <= SRS(m))) flux_flag(m) = 60
       
       if ((SRS(m) <= 0.00D0) .and. (0.00D0 <= SR(m)))  flux_flag(m) = 50

       if (SR(m) < 0.00D0) flux_flag(m) = 40

    end do



  ! Finding the inter-cell fluxes.
  ! ------------------------------

    do m = BOUNDARY, rowsize + BOUNDARY


       if (flux_flag(m) < 35) then

          density_flux(m)    = flux_L(m, 1)
          x_momentum_flux(m) = flux_L(m, 2)
          y_momentum_flux(m) = flux_L(m, 3)
          z_momentum_flux(m) = flux_L(m, 4)
          energy_flux(m)     = flux_L(m, 5)

          x_magfield_flux(m) = flux_L(m, 6)
          y_magfield_flux(m) = flux_L(m, 7)
          z_magfield_flux(m) = flux_L(m, 8)

       end if


       if (flux_flag(m) < 25) then

          x_momentum_SX1 = density_SL1(m) * x_velocity_S1(m)
          y_momentum_SX1 = density_SL1(m) * y_velocity_SL1(m)
          z_momentum_SX1 = density_SL1(m) * z_velocity_SL1(m)

          density_flux(m)    = density_flux(m)    + SL(m) * (density_SL1(m) - density_L(m))
          x_momentum_flux(m) = x_momentum_flux(m) + SL(m) * (x_momentum_SX1 - x_momentum_L(m))
          y_momentum_flux(m) = y_momentum_flux(m) + SL(m) * (y_momentum_SX1 - y_momentum_L(m))
          z_momentum_flux(m) = z_momentum_flux(m) + SL(m) * (z_momentum_SX1 - z_momentum_L(m))
          energy_flux(m)     = energy_flux(m)     + SL(m) * (energy_SL1(m)  - energy_L(m))

          x_magfield_flux(m) = x_magfield_flux(m) + SL(m) * (x_magfield_S1(m)  - x_magfield_L(m))
          y_magfield_flux(m) = y_magfield_flux(m) + SL(m) * (y_magfield_SL1(m) - y_magfield_L(m))
          z_magfield_flux(m) = z_magfield_flux(m) + SL(m) * (z_magfield_SL1(M) - z_magfield_L(m))

       end if


       if (flux_flag(m) < 15) then

          density_flux(m)    = density_flux(m)    + SLS(m) * (density_SL2(m) - density_SL1(m))
          x_momentum_flux(m) = x_momentum_flux(m) + SLS(m) * (density_SL2(m) * x_velocity_S2(m) - x_momentum_SX1)
          y_momentum_flux(m) = y_momentum_flux(m) + SLS(m) * (density_SL2(m) * y_velocity_S2(m) - y_momentum_SX1)
          z_momentum_flux(m) = z_momentum_flux(m) + SLS(m) * (density_SL2(m) * z_velocity_S2(m) - z_momentum_SX1)
          energy_flux(m)     = energy_flux(m)     + SLS(m) * (energy_SL2(m)  - energy_SL1(m))

          x_magfield_flux(m) = x_magfield_flux(m) + SLS(m) * (x_magfield_S2(m) - x_magfield_S1(m))
          y_magfield_flux(m) = y_magfield_flux(m) + SLS(m) * (y_magfield_S2(m) - y_magfield_SL1(m))
          z_magfield_flux(m) = z_magfield_flux(m) + SLS(m) * (z_magfield_S2(m) - z_magfield_SL1(m))

       end if


       if (flux_flag(m) > 35) then

          density_flux(m)    = flux_R(m, 1)
          x_momentum_flux(m) = flux_R(m, 2)
          y_momentum_flux(m) = flux_R(m, 3)
          z_momentum_flux(m) = flux_R(m, 4)
          energy_flux(m)     = flux_R(m, 5)

          x_magfield_flux(m) = flux_R(m, 6)
          y_magfield_flux(m) = flux_R(m, 7)
          z_magfield_flux(m) = flux_R(m, 8)

       end if


       if (flux_flag(m) > 45) then

          x_momentum_SX1 = density_SR1(m) * x_velocity_S1(m)
          y_momentum_SX1 = density_SR1(m) * y_velocity_SR1(m)
          z_momentum_SX1 = density_SR1(m) * z_velocity_SR1(m)

          density_flux(m)    = density_flux(m)    + SR(m) * (density_SR1(m) - density_R(m))
          x_momentum_flux(m) = x_momentum_flux(m) + SR(m) * (x_momentum_SX1 - x_momentum_R(m))
          y_momentum_flux(m) = y_momentum_flux(m) + SR(m) * (y_momentum_SX1 - y_momentum_R(m))
          z_momentum_flux(m) = z_momentum_flux(m) + SR(m) * (z_momentum_SX1 - z_momentum_R(m))
          energy_flux(m)     = energy_flux(m)     + SR(m) * (energy_SR1(m)  - energy_R(m))

          x_magfield_flux(m) = x_magfield_flux(m) + SR(m) * (x_magfield_S1(m)  - x_magfield_R(m))
          y_magfield_flux(m) = y_magfield_flux(m) + SR(m) * (y_magfield_SR1(m) - y_magfield_R(m))
          z_magfield_flux(m) = z_magfield_flux(m) + SR(m) * (z_magfield_SR1(M) - z_magfield_R(m))

       end if



       if (flux_flag(m) > 55) then

          density_flux(m)    = density_flux(m)    + SRS(m) * (density_SR2(m) - density_SR1(m))
          x_momentum_flux(m) = x_momentum_flux(m) + SRS(m) * (density_SR2(m) * x_velocity_S2(m) - x_momentum_SX1)
          y_momentum_flux(m) = y_momentum_flux(m) + SRS(m) * (density_SR2(m) * y_velocity_S2(m) - y_momentum_SX1)
          z_momentum_flux(m) = z_momentum_flux(m) + SRS(m) * (density_SR2(m) * z_velocity_S2(m) - z_momentum_SX1)
          energy_flux(m)     = energy_flux(m)     + SRS(m) * (energy_SR2(m)  - energy_SR1(m))

          x_magfield_flux(m) = x_magfield_flux(m) + SRS(m) * (x_magfield_S2(m) - x_magfield_S1(m))
          y_magfield_flux(m) = y_magfield_flux(m) + SRS(m) * (y_magfield_S2(m) - y_magfield_SR1(m))
          z_magfield_flux(m) = z_magfield_flux(m) + SRS(m) * (z_magfield_S2(m) - z_magfield_SR1(m))

       end if

    end do

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Calculate_intercell_fluxes_HLLD
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Calculating the Roe Averages --------------------------------------------------------------------------- [SOL03] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  function roeavg (rhoL, rhoR, varL, varR)

    real (PREC), intent (in) :: rhoL, rhoR
    real (PREC), intent (in) :: varL, varR
    real (PREC) :: roeavg

    roeavg = (dsqrt(rhoL) * varL) + (dsqrt(rhoR) * varR)
    roeavg = roeavg / (dsqrt(rhoL) + dsqrt(rhoR))

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end function roeavg
! ----------------------------------------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Calculating the HLL State ------------------------------------------------------------------------------ [SOL04] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  function HLL_State (SL, SR, FL, FR, UL, UR)

    real (PREC), intent (in) :: SL, SR, FL, FR, UL, UR
    real (PREC) :: scratch, HLL_State

    scratch = SR - SL

    HLL_State = (SR * UR) - (SL * UL) - (FR - FL)
    HLL_State = HLL_State / scratch

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end function HLL_State
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Deallocate the Solver -------------------------------------------------------------------------------- [SOL05] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Deallocate_solver ()


  ! Deallocating the inter-cell fluxes.
  ! ---------------------------------

    deallocate (density_flux)
    deallocate (x_momentum_flux)
    deallocate (y_momentum_flux)
    deallocate (z_momentum_flux)
    deallocate (energy_flux)

    deallocate (x_magfield_flux)
    deallocate (y_magfield_flux)
    deallocate (z_magfield_flux)



  ! Deallocating the working variables.
  ! -----------------------------------

    deallocate (flux_L); deallocate (flux_R)

    deallocate (SL);  deallocate (SR);  deallocate (SM)
    deallocate (SLS); deallocate (SRS);

    deallocate (pressure_S1); deallocate (x_velocity_S1); deallocate (x_magfield_S1)

    deallocate (density_SL1);    deallocate (energy_SL1)
    deallocate (y_velocity_SL1); deallocate (y_magfield_SL1)
    deallocate (z_velocity_SL1); deallocate (z_magfield_SL1)

    deallocate (density_SR1);    deallocate (energy_SR1)
    deallocate (y_velocity_SR1); deallocate (y_magfield_SR1)
    deallocate (z_velocity_SR1); deallocate (z_magfield_SR1)

    deallocate (pressure_S2)
    deallocate (x_velocity_S2);  deallocate (x_magfield_S2)
    deallocate (y_velocity_S2);  deallocate (y_magfield_S2)
    deallocate (z_velocity_S2);  deallocate (z_magfield_S2)

    deallocate (density_SL2);    deallocate (energy_SL2);
    deallocate (density_SR2);    deallocate (energy_SR2);

    deallocate (flux_flag)

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Deallocate_solver
! ----------------------------------------------------------------------------------------------------------------------------------



end module solve
