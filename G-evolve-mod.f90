! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
!
!     G: The Evolve Module
!     ====================
!     -- Main actions: evolve the entire grid for one timestep; evolve individual rows.
!
! ----------------------------------------------------------------------------------------------------------------------------------
! --- This code is best viewed with a window at least 135 characters wide. ---------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------



module evolve

  use setup
  use fluid
  use tests
  use grids
  use recon
  use solve
  implicit none


  ! The updated state variables.
  ! ----------------------------

  real (PREC), dimension(:), allocatable :: density_new
  real (PREC), dimension(:), allocatable :: x_velocity_new
  real (PREC), dimension(:), allocatable :: y_velocity_new
  real (PREC), dimension(:), allocatable :: z_velocity_new
  real (PREC), dimension(:), allocatable :: energy_new
  
  real (PREC), dimension(:), allocatable :: x_magfield_new
  real (PREC), dimension(:), allocatable :: y_magfield_new
  real (PREC), dimension(:), allocatable :: z_magfield_new

contains



! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Initialise the Simulation ---------------------------------------------------------------------------- [EVO01] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Initialise_simulation ()


  ! Setting initial conditions.
  ! ---------------------------

    if (TEST_PROBLEM == 'A') call Set_initial_conditions_A ()     ! --> <TES01>
    if (TEST_PROBLEM == 'B') call Set_initial_conditions_B ()     ! --> <TES02>
    if (TEST_PROBLEM == 'C') call Set_initial_conditions_C ()     ! --> <TES03>
    if (TEST_PROBLEM == 'D') call Set_initial_conditions_D ()     ! --> <TES04>
    if (TEST_PROBLEM == 'E') call Set_initial_conditions_E ()     ! --> <TES05>
    if (TEST_PROBLEM == 'F') call Set_initial_conditions_F ()     ! --> <TES06>
    if (TEST_PROBLEM == 'X') call Load_restart_files ()           ! --> <TES07>


  ! Allocating arrays.
  ! ------------------

    call Initialise_grids ()            ! --> <GRI01>
    call Initialise_input_states ()     ! --> <REC01>
    call Initialise_solver ()           ! --> <SOL01>



  ! Allocating the updated state variables.
  ! ---------------------------------------

    allocate (density_new   (MFULL))
    allocate (x_velocity_new(MFULL))
    allocate (y_velocity_new(MFULL))
    allocate (z_velocity_new(MFULL))
    allocate (energy_new    (MFULL))

    allocate (x_magfield_new(MFULL))
    allocate (y_magfield_new(MFULL))
    allocate (z_magfield_new(MFULL))

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Initialise_simulation
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Evolve the Grid -------------------------------------------------------------------------------------- [EVO02] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Evolve_grid ()


  ! --------------------------------------------------------------------------------------------------------------------------------
  ! [EVO_A] Guide: Evolving the Grid: 1D vs 3D     
  ! --------------------------------------------------------------------------------------------------------------------------------
  !   The method used for evolving the magnetic field in 3D (while maintaining the divergence-free condition) does not reduce
  !     to the 1D method when evolving 1D systems. Hence, we use "if (DIMENSIONS > 1)" conditionals.
  ! --------------------------------------------------------------------------------------------------------------------------------


  ! First dimensional passes.
  ! -------------------------

    if (DIMENSIONS > 1) call Initialise_magnetic_field_arrays ()     ! --> <GRI03>

    call Evolve_x_direction ()     ! --> <EVO03>

    if (DIMENSIONS > 1) then

       call Evolve_y_direction ()     ! --> <EVO04>
       call Evolve_z_direction ()     ! --> <EVO05>
       
       call Finish_magnetic_field_evolution ()     ! --> <EVO07>

    end if



  ! Second dimensional sweeps.
  ! --------------------------

    if (DIMENSIONS > 1) then

       call Initialise_magnetic_field_arrays ()

       call Evolve_z_direction ()
       call Evolve_y_direction ()

    end if

    call Evolve_x_direction ()

    if (DIMENSIONS > 1) call Initialise_magnetic_field_arrays ()

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Evolve_grid
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Evolve the x-Direction ------------------------------------------------------------------------------- [EVO03] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Evolve_x_direction ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, j, k



  ! Evolve fluid.
  ! -------------

    do k = 1, ZSIZE
       do j = 1, YSIZE

          call Clear_1D_data ()         ! --> <GRI02>
          call Extract_row_x (j, k)     ! --> <GRI04>
          call Evolve_row ()            ! --> <EVO06>
          call Return_row_x (j, k)      ! --> <GRI08>

       end do
    end do

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Evolve_x_direction
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Evolve the y-Direction ------------------------------------------------------------------------------- [EV004] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Evolve_y_direction ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, j, k



  ! Evolve y-Direction. 
  ! -------------------

    do k = 1, ZSIZE
       do i = 1, XSIZE

          call Clear_1D_data ()         ! --> <GRI02>
          call Extract_row_y (i, k)     ! --> <GRI05>
          call Evolve_row ()            ! --> <EVO06>
          call Return_row_y (i, k)      ! --> <GRI09>

       end do
    end do

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Evolve_y_direction
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Evolve the z-Direction ------------------------------------------------------------------------------- [EV005] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Evolve_z_direction ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: i, j, k



  ! Evolve fluid.
  ! -------------

    do j = 1, YSIZE
       do i = 1, XSIZE

          call Clear_1D_data ()         ! --> <GRI02>
          call Extract_row_z (i, j)     ! --> <GRI06>
          call Evolve_row ()            ! --> <EVO06>
          call Return_row_z (i, j)      ! --> <GRI10>

       end do
    end do

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Evolve_z_direction
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Evolve the Row --------------------------------------------------------------------------------------- [EV006] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Evolve_row ()


  ! Declaration of local variables.
  ! -------------------------------

    real (PREC), dimension(3) :: velc, magf
    real (PREC) :: velc_squrd, magf_squrd

    integer :: m



  ! Calculating the inter-cell fluxes.
  ! ----------------------------------

    if (RECONSTRUCT_TYPE == 'T') call Determine_Riemann_input_TVD       ! --> <REC02>
    if (RECONSTRUCT_TYPE == 'P') call Determine_Riemann_input_PPM       ! --> <REC03>
    if (RECONSTRUCT_TYPE == 'C') call Determine_Riemann_input_PPMCS     ! --> <REC04>

    call Calculate_intercell_fluxes_HLLD ()     ! --> <SOL02>



  ! Determining new fluid quantities.
  ! ---------------------------------

    do m = 1 + BOUNDARY, rowsize + BOUNDARY

       density_new(m)    = density_1D(m) - dtdx(m) * (density_flux(m) - density_flux(m-1))

       x_velocity_new(m) = density_1D(m) * x_velocity_1D(m) - dtdx(m) * (x_momentum_flux(m) - x_momentum_flux(m-1))
       x_velocity_new(m) = x_velocity_new(m) / density_new(m)

       y_velocity_new(m) = density_1D(m) * y_velocity_1D(m) - dtdx(m) * (y_momentum_flux(m) - y_momentum_flux(m-1))
       y_velocity_new(m) = y_velocity_new(m) / density_new(m)

       z_velocity_new(m) = density_1D(m) * z_velocity_1D(m) - dtdx(m) * (z_momentum_flux(m) - z_momentum_flux(m-1))
       z_velocity_new(m) = z_velocity_new(m) / density_new(m)

       energy_new(m)     = energy_1D(m) - dtdx(m) * (energy_flux(m) - energy_flux(m-1))

    end do



  ! Determining new field quantities.
  ! ---------------------------------

    do m = 1 + BOUNDARY, rowsize + BOUNDARY
       
       x_magfield_new(m) = x_magfield_1D(m) - dtdx(m) * (x_magfield_flux(m) - x_magfield_flux(m-1))
       y_magfield_new(m) = y_magfield_1D(m) - dtdx(m) * (y_magfield_flux(m) - y_magfield_flux(m-1))
       z_magfield_new(m) = z_magfield_1D(m) - dtdx(m) * (z_magfield_flux(m) - z_magfield_flux(m-1))
       
    end do
    
    if (DIMENSIONS > 1) then
       
       do m = 1 + BOUNDARY, rowsize + BOUNDARY
          
          flux_y1(m) = (0.125D0 * dtdy(m)) * ( y_magfield_flux(m) + y_magfield_flux(m-1))
          flux_y2(m) = (0.125D0 * dtdx(m)) * (-y_magfield_flux(m) + y_magfield_flux(m-1))
          
          flux_z1(m) = (0.125D0 * dtdz(m)) * ( z_magfield_flux(m) + z_magfield_flux(m-1))
          flux_z2(m) = (0.125D0 * dtdx(m)) * (-z_magfield_flux(m) + z_magfield_flux(m-1))
          
       end do

    end if



  ! Saving new quantities.
  ! ----------------------

    do m = 1 + BOUNDARY, rowsize + BOUNDARY

       density_1D(m)    = density_new(m)

       x_velocity_1D(m) = x_velocity_new(m)
       y_velocity_1D(m) = y_velocity_new(m)
       z_velocity_1D(m) = z_velocity_new(m)

       x_momentum_1D(m) = x_velocity_new(m) * density_new(m)
       y_momentum_1D(m) = y_velocity_new(m) * density_new(m)
       z_momentum_1D(m) = z_velocity_new(m) * density_new(m)
  
       energy_1D(m)     = energy_new(m)

       x_magfield_1D(m) = x_magfield_new(m)
       y_magfield_1D(m) = y_magfield_new(m)
       z_magfield_1D(m) = z_magfield_new(m)
       
       velc = Create_vector (x_velocity_1D(m), y_velocity_1D(m), z_velocity_1D(m))
       magf = Create_vector (x_magfield_1D(m), y_magfield_1D(m), z_magfield_1D(m))
       
       velc_squrd = Dotproduct (velc, velc)
       magf_squrd = Dotproduct (magf, magf)
       
       pressure_1D(m) = Calculate_pressure_EOS (density_1D(m), energy_1D(m), gamma_1D(m), velc_squrd, magf_squrd)

    end do

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Evolve_row
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Finishing the Magnetic Field Evolution --------------------------------------------------------------- [EVO07] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Finish_magnetic_field_evolution ()

  ! Declaration of local variables.
  ! -------------------------------

    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd, scratch

    integer :: i, j, k



  ! Subtracting magnetic field from pressure and energy.
  ! ----------------------------------------------------

    if (PRESSURE_FIX) then

       do k = 1, ZSIZE
          do j = 1, YSIZE
             do i = 1, XSIZE
                
                magf = Create_vector (x_magfield(i, j, k), y_magfield(i, j, k), z_magfield(i, j, k))
                magf_squrd = Dotproduct (magf, magf)
                
                scratch = (0.5D0 * magf_squrd * MHDF(2))
                
                pressure(i, j, k) = pressure(i, j, k) - (scratch * (2.0D0 - gamma(i, j, k)))
                energy(i, j, k) = energy(i, j, k) - scratch
                
             end do
          end do
       end do
       
    end if


  ! Adding temporary magnetic field.
  ! --------------------------------

    do k = 1, ZSIZE
       do j = 1, YSIZE
          do i = 1, XSIZE

             x_magfield(i, j, k) = x_magfield_stored(i, j, k) + delta_bx(i, j, k)
             y_magfield(i, j, k) = y_magfield_stored(i, j, k) + delta_by(i, j, k)
             z_magfield(i, j, k) = z_magfield_stored(i, j, k) + delta_bz(i, j, k)

          end do
       end do
    end do



   ! Updating the pressure.
   ! ----------------------

    if (PRESSURE_FIX) then

       do k = 1, ZSIZE
          do j = 1, YSIZE
             do i = 1, XSIZE
                
                magf = Create_vector (x_magfield(i, j, k), y_magfield(i, j, k), z_magfield(i, j, k))
                magf_squrd = Dotproduct (magf, magf)
                
                scratch = (0.5D0 * magf_squrd * MHDF(2))
                
                pressure(i, j, k) = pressure(i, j, k) + (scratch * (2.0D0 - gamma(i, j, k)))
                energy(i, j, k) = energy(i, j, k) + scratch
                
                
             end do
          end do
       end do
       
    end if

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Finish_magnetic_field_evolution
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Shut Down the Simulation ----------------------------------------------------------------------------- [EVO08] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Shut_down_simulation ()


  ! Deallocating arrays.
  ! --------------------

    call Deallocate_fluid ()           ! --> <FLU12>
    call Deallocate_grids ()           ! --> <GRI11>
    call Deallocate_input_states()     ! --> <REC06>
    call Deallocate_solver ()          ! --> <SOL05>



  ! Deallocating the updated state variables.
  ! -----------------------------------------

    deallocate (density_new)
    deallocate (x_velocity_new)
    deallocate (y_velocity_new)
    deallocate (z_velocity_new)
    deallocate (energy_new)
    
    deallocate (x_magfield_new)
    deallocate (y_magfield_new)
    deallocate (z_magfield_new)

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Shut_down_simulation
! ----------------------------------------------------------------------------------------------------------------------------------



end module evolve
