! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
! 
!     B: The Fluid Module
!     ===================
!     -- Main variables: the 3D fluid and magnetic field variables.
!     -- Main actions: EOS calculations; determine fluid properties; output data to file.
!
! ----------------------------------------------------------------------------------------------------------------------------------
! --- This code is best viewed with a window at least 135 characters wide. ---------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------



module fluid

  use setup
  implicit none


  ! The fluid variables.
  ! --------------------

  real (PREC), dimension(:,:,:), allocatable :: density
  real (PREC), dimension(:,:,:), allocatable :: x_velocity
  real (PREC), dimension(:,:,:), allocatable :: y_velocity
  real (PREC), dimension(:,:,:), allocatable :: z_velocity
  real (PREC), dimension(:,:,:), allocatable :: x_momentum
  real (PREC), dimension(:,:,:), allocatable :: y_momentum
  real (PREC), dimension(:,:,:), allocatable :: z_momentum
  real (PREC), dimension(:,:,:), allocatable :: energy
  real (PREC), dimension(:,:,:), allocatable :: pressure
  real (PREC), dimension(:,:,:), allocatable :: gamma



  ! The field variables.
  ! --------------------

  real (PREC), dimension(:,:,:), allocatable :: x_magfield
  real (PREC), dimension(:,:,:), allocatable :: y_magfield
  real (PREC), dimension(:,:,:), allocatable :: z_magfield

  real (PREC), dimension(:,:,:), allocatable :: x_magfield_stored
  real (PREC), dimension(:,:,:), allocatable :: y_magfield_stored
  real (PREC), dimension(:,:,:), allocatable :: z_magfield_stored

  real (PREC), dimension(:,:,:), allocatable :: delta_bx
  real (PREC), dimension(:,:,:), allocatable :: delta_by
  real (PREC), dimension(:,:,:), allocatable :: delta_bz

contains



! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Initialise Fluid ------------------------------------------------------------------------------------- [FLU01] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Initialise_fluid ()


  ! Setting constants.
  ! ------------------

    DELTAX = XDOMN / real(XSIZE)
    DELTAY = YDOMN / real(YSIZE)
    DELTAZ = ZDOMN / real(ZSIZE)

    XFULL = (2 * BOUNDARY) + XSIZE
    YFULL = (2 * BOUNDARY) + YSIZE
    ZFULL = (2 * BOUNDARY) + ZSIZE

    NTIME  = BASE_NTIME
    NSTEPS = ceiling(FULLTIME/DELTAT)

    MSIZE = max(XSIZE, YSIZE, ZSIZE)
    MFULL = max(XFULL, YFULL, ZFULL)



  ! Allocating the fluid variables.
  ! -------------------------------

    allocate (density   (XSIZE, YSIZE, ZSIZE))
    allocate (x_velocity(XSIZE, YSIZE, ZSIZE))
    allocate (y_velocity(XSIZE, YSIZE, ZSIZE))
    allocate (z_velocity(XSIZE, YSIZE, ZSIZE))
    allocate (x_momentum(XSIZE, YSIZE, ZSIZE))
    allocate (y_momentum(XSIZE, YSIZE, ZSIZE))
    allocate (z_momentum(XSIZE, YSIZE, ZSIZE))
    allocate (energy    (XSIZE, YSIZE, ZSIZE))
    allocate (pressure  (XSIZE, YSIZE, ZSIZE))
    allocate (gamma     (XSIZE, YSIZE, ZSIZE))



  ! Allocating the field variables.
  ! -------------------------------
    
    allocate (x_magfield(XSIZE, YSIZE, ZSIZE))
    allocate (y_magfield(XSIZE, YSIZE, ZSIZE))
    allocate (z_magfield(XSIZE, YSIZE, ZSIZE))

    if (DIMENSIONS > 1) then

       allocate (x_magfield_stored(XSIZE, YSIZE, ZSIZE))
       allocate (y_magfield_stored(XSIZE, YSIZE, ZSIZE))
       allocate (z_magfield_stored(XSIZE, YSIZE, ZSIZE))
       
       allocate (delta_bx(XSIZE, YSIZE, ZSIZE))
       allocate (delta_by(XSIZE, YSIZE, ZSIZE))
       allocate (delta_bz(XSIZE, YSIZE, ZSIZE))
       
    end if
    
    return
   


! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Initialise_fluid
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Calculate the Energy using an Ideal Equation of State -------------------------------------------------- [FLU02] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  function Calculate_energy_EOS (dens_val, pres_val, gamma_val, velc_squrd, magf_squrd)

    ! Note: this function expects the total pressure (gas + magnetic). 

    real (PREC), intent (in) :: dens_val, pres_val, gamma_val, velc_squrd, magf_squrd
    real (PREC)              :: Calculate_energy_EOS, scratch1, scratch2, scratch3

    scratch1 = Calculate_gas_pressure (pres_val, magf_squrd) / (gamma_val - 1.0D0)
    scratch2 = 0.5D0 * dens_val * velc_squrd
    scratch3 = 0.5D0 * magf_squrd * MHDF(2)

    Calculate_energy_EOS = scratch1 + scratch2 + scratch3

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end function Calculate_energy_EOS
! ----------------------------------------------------------------------------------------------------------------------------------


! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Calculate the Pressure using an Ideal Equation of State ------------------------------------------------ [FLU03] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  function Calculate_pressure_EOS (dens_val, ergy_val, gamma_val, velc_squrd, magf_squrd)

    ! Note: this function returns the total pressure (gas + magnetic). 

    real (PREC), intent (in) :: dens_val, ergy_val, gamma_val, velc_squrd, magf_squrd
    real (PREC)              :: Calculate_pressure_EOS, scratch1, scratch2, scratch3
 
    scratch1 = gamma_val - 1.0D0
    scratch2 = ergy_val - (0.5D0 * dens_val * velc_squrd) - (0.5D0 * magf_squrd * MHDF(2))
    scratch3 = scratch1 * scratch2

    Calculate_pressure_EOS = Calculate_total_pressure (scratch3, magf_squrd)

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end function Calculate_pressure_EOS
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Calculate Total Pressure ------------------------------------------------------------------------------- [FLU04] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  function Calculate_total_pressure (thermal_pressure, magf_squrd)

    real (PREC), intent (in) :: thermal_pressure, magf_squrd
    real (PREC)              :: Calculate_total_pressure

    Calculate_total_pressure = thermal_pressure + 0.5D0 * magf_squrd * MHDF(2)

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end function Calculate_total_pressure
! ----------------------------------------------------------------------------------------------------------------------------------


! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Calculate Gas Pressure (from total pressure) ----------------------------------------------------------- [FLU05] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  function Calculate_gas_pressure (total_pressure, magf_squrd)

    real (PREC), intent (in) :: total_pressure, magf_squrd
    real (PREC)              :: Calculate_gas_pressure

    Calculate_gas_pressure = total_pressure - 0.5D0 * magf_squrd * MHDF(2)

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end function Calculate_gas_pressure
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Determine MHD Factor (depending on units used) --------------------------------------------------------- [FLU06] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  function MHDF (power)

    integer, intent (in) :: power
    real (PREC) :: MHDF

    ! Note: see <SET_C> for explanation of how this works. The value of
    !   MHDUNITS is set with the initial conditions in the <TESxx> subroutines.

    if (MHDUNITS == 1) then

       MHDF = 1.0D0 / dsqrt(4.0D0 * PI)
       MHDF = MHDF**real(power)

    else

       MHDF = 1.0D0

    end if

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end function MHDF
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Determine the Timestep (using CFL condition) --------------------------------------------------------- [FLU07] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Determine_timestep ()


  ! Declaration of local variables.
  ! -------------------------------

    real (PREC) :: x_max_sigspeed
    real (PREC) :: y_max_sigspeed
    real (PREC) :: z_max_sigspeed

    real (PREC) :: gas_pressure
    real (PREC) :: sigspeed

    real (PREC), dimension(3) :: magf
    real (PREC) :: magf_squrd

    real (PREC) :: a2, Ca2
    real (PREC) :: Cax2, Cay2, Caz2

    real (PREC), parameter :: epsilon = 1.0D-10
    integer :: i, j, k



  ! Finding the largest signal speed.
  ! ---------------------------------

    x_max_sigspeed = -1.0D-10
    y_max_sigspeed = -1.0D-10
    z_max_sigspeed = -1.0D-10

    do k = 1, ZSIZE
       do j = 1, YSIZE
          do i = 1, XSIZE

             magf = Create_vector (x_magfield(i, j, k), y_magfield(i, j, k), z_magfield(i, j, k))
             magf_squrd = Dotproduct (magf, magf)

             gas_pressure = Calculate_gas_pressure (pressure(i, j, k), magf_squrd)


             ! a2 is the speed of sound squared

             a2 = (gamma(i, j, k) * gas_pressure) / density(i, j, k)


             ! Cax2 is the x-component of the Alfven wavespeed squared

             Cax2 = (x_magfield(i, j, k)**2.0D0 / density(i, j, k)) * MHDF(2)
             Cay2 = (y_magfield(i, j, k)**2.0D0 / density(i, j, k)) * MHDF(2)
             Caz2 = (z_magfield(i, j, k)**2.0D0 / density(i, j, k)) * MHDF(2)

             Ca2 = Cax2 + Cay2 + Caz2

             
             ! The x-direction signal speed

             sigspeed = (a2 + Ca2) + dsqrt((a2 + Ca2)**2.0D0 - (4.0D0 * a2 * Cax2))
             sigspeed = dsqrt(0.5D0 * sigspeed)

             sigspeed = sigspeed + abs(x_velocity(i, j, k))

             if (sigspeed == sigspeed) then     ! check for NaNs
                if (sigspeed > x_max_sigspeed) x_max_sigspeed = sigspeed
             end if


             ! The y-direction signal speed

             sigspeed = (a2 + Ca2) + dsqrt((a2 + Ca2)**2.0D0 - (4.0D0 * a2 * Cay2))
             sigspeed = dsqrt(0.5D0 * sigspeed)

             sigspeed = sigspeed + abs(y_velocity(i, j, k))

             if (sigspeed == sigspeed) then     ! check for NaNs
                if (sigspeed > y_max_sigspeed) y_max_sigspeed = sigspeed
             end if


             ! The z-direction signal speed

             sigspeed = (a2 + Ca2) + dsqrt((a2 + Ca2)**2.0D0 - (4.0D0 * a2 * Caz2))
             sigspeed = dsqrt(0.5D0 * sigspeed)

             sigspeed = sigspeed + abs(z_velocity(i, j, k))

             if (sigspeed == sigspeed) then     ! check for NaNs
                if (sigspeed > z_max_sigspeed) z_max_sigspeed = sigspeed
             end if

          end do
       end do
    end do



  ! Calculating the new timestep.
  ! -----------------------------

    x_max_sigspeed = max(x_max_sigspeed, epsilon)
    y_max_sigspeed = max(y_max_sigspeed, epsilon)
    z_max_sigspeed = max(z_max_sigspeed, epsilon)
             
    DELTAT = COURANT * min((DELTAX/x_max_sigspeed), (DELTAY/y_max_sigspeed), (DELTAZ/z_max_sigspeed))

    if (DELTAT > MAX_DELTAT) DELTAT = MAX_DELTAT

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Determine_timestep
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Output Data ------------------------------------------------------------------------------------------ [FLU08] ---
! ----------------------------------------------------------------------------------------------------------------------------------
! --- Version A: 1D Output (density; pressure; y-velocity; y-magnetic field) -------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Output_data_A (i, j, k, num)

    integer, intent (in) :: i, j, k, num


  ! --------------------------------------------------------------------------------------------------------------------------------
  ! [FLU_A] Guide: Output Data - Version A         
  ! --------------------------------------------------------------------------------------------------------------------------------
  !   Quantities written: density; (thermal) pressure; y-velocity and y-magnetic field.                                
  !
  !     Usage: choose slice direction by setting appropriate parameter to -1. 
  !     e.g. "call Output_data_A (-1, j_fixed, k_fixed, num)" for x-direction slice.
  ! --------------------------------------------------------------------------------------------------------------------------------


  ! Declaration of local variables.
  ! -------------------------------

    character :: dens_file*(19)
    character :: pres_file*(19)
    character :: yvel_file*(19)
    character :: ymag_file*(19)
    character :: chari*(4)

    real (PREC), dimension(3) :: magf
    real (PREC) :: magf_squrd
    real (PREC) :: gas_pressure, pos
    integer :: m



  ! Obtaining file names.
  ! ---------------------

    write(chari, '(i4.4)') num
    chari = adjustl(chari)

    dens_file = "output/dens"//chari//".txt"
    pres_file = "output/pres"//chari//".txt"
    yvel_file = "output/yvel"//chari//".txt"
    ymag_file = "output/ymag"//chari//".txt"



  ! Opening files.
  ! --------------

    open(10, file = dens_file, form = 'formatted')
    open(20, file = pres_file, form = 'formatted')
    open(30, file = yvel_file, form = 'formatted')
    open(40, file = ymag_file, form = 'formatted')



  ! The x-direction slice...
  ! ------------------------

    if(i < 0) then

       do m = 1, XSIZE

          magf = Create_vector (x_magfield(m, j, k), y_magfield(m, j, k), z_magfield(m, j, k))
          magf_squrd = Dotproduct (magf, magf)

          gas_pressure = Calculate_gas_pressure (pressure(m, j, k), magf_squrd)

          pos  = ((real(m) - 0.5D0)/real(XSIZE)) * XDOMN

          write(10, '(f12.4, a, f12.4)'), pos, ' ', density(m, j, k)
          write(20, '(f12.4, a, f12.4)'), pos, ' ', gas_pressure
          write(30, '(f12.4, a, f12.4)'), pos, ' ', y_velocity(m, j, k)
          write(40, '(f12.4, a, f12.4)'), pos, ' ', y_magfield(m, j, k)

       end do

    end if



  ! The y-direction slice...
  ! ------------------------

    if(j < 0) then

       do m = 1, YSIZE

          magf = Create_vector (x_magfield(i, m, k), y_magfield(i, m, k), z_magfield(i, m, k))
          magf_squrd = Dotproduct (magf, magf)

          gas_pressure = Calculate_gas_pressure (pressure(i, m, k), magf_squrd)

          pos  = ((real(m) - 0.5D0)/real(YSIZE)) * YDOMN

          write(10, '(f12.4, a, f12.4)'), pos, ' ', density(i, m, k)
          write(20, '(f12.4, a, f12.4)'), pos, ' ', gas_pressure
          write(30, '(f12.4, a, f12.4)'), pos, ' ', y_velocity(i, m, k)
          write(40, '(f12.4, a, f12.4)'), pos, ' ', y_magfield(i, m, k)

       end do

    end if



  ! The z-direction slice...
  ! ------------------------

    if(k < 0) then

       do m = 1, ZSIZE

          magf = Create_vector (x_magfield(i, j, m), y_magfield(i, j, m), z_magfield(i, j, m))
          magf_squrd = Dotproduct (magf, magf)

          gas_pressure = Calculate_gas_pressure (pressure(i, j, m), magf_squrd)

          pos  = ((real(m) - 0.5D0)/real(ZSIZE)) * ZDOMN

          write(10, '(f12.4, a, f12.4)'), pos, ' ', density(i, j, m)
          write(20, '(f12.4, a, f12.4)'), pos, ' ', gas_pressure
          write(30, '(f12.4, a, f12.4)'), pos, ' ', y_velocity(i, j, m)
          write(40, '(f12.4, a, f12.4)'), pos, ' ', y_magfield(i, j, m)

       end do

    end if



  ! Closing files.
  ! --------------

    close (10)
    close (20)
    close (30)
    close (40)

    return
    


! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Output_data_A
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Output Data ------------------------------------------------------------------------------------------ [FLU09] ---
! ----------------------------------------------------------------------------------------------------------------------------------
! --- Version B: 2D Output (density; thermal pressure; magnetic pressure; velocity squared) ----------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Output_data_B (num)

    integer, intent (in) :: num


  ! Declaration of local variables.
  ! -------------------------------

    character :: dens_file*(19)
    character :: pthe_file*(19)
    character :: pmag_file*(19)
    character :: vels_file*(19)
    character :: chari*(4)

    real (PREC), dimension(3) :: velc, magf
    real (PREC) :: velc_squrd, magf_squrd
    real (PREC) :: gas_pressure, mag_pressure
    real (PREC) :: xpos, ypos
    integer :: i, j



  ! Obtaining file names.
  ! ---------------------

    write(chari, '(i4.4)') num
    chari = adjustl(chari)

    dens_file = "output/dens"//chari//".txt"
    pthe_file = "output/pthe"//chari//".txt"
    pmag_file = "output/pmag"//chari//".txt"
    vels_file = "output/vels"//chari//".txt"



  ! Opening files.
  ! --------------

    open(10, file = dens_file, form = 'formatted')
    open(20, file = pthe_file, form = 'formatted')
    open(30, file = pmag_file, form = 'formatted')
    open(40, file = vels_file, form = 'formatted')



  ! The 2D slice...
  ! ---------------

    do j = 1, YSIZE
       do i = 1, XSIZE

          xpos = ((real(i) - 0.5D0)/real(XSIZE)) * XDOMN
          ypos = ((real(j) - 0.5D0)/real(YSIZE)) * YDOMN

          velc = Create_vector (x_velocity(i, j, 1), y_velocity(i, j, 1), z_velocity(i, j, 1))
          magf = Create_vector (x_magfield(i, j, 1), y_magfield(i, j, 1), z_velocity(i, j, 1))

          velc_squrd = Dotproduct (velc, velc)
          magf_squrd = Dotproduct (magf, magf)

          gas_pressure = Calculate_gas_pressure (pressure(i, j, 1), magf_squrd)
          mag_pressure = 0.5D0 * magf_squrd * MHDF(2)

          write(10, '(f12.4, a, f12.4, a, f12.4)'), xpos, ' ', ypos, ' ', density(i, j, 1)
          write(20, '(f12.4, a, f12.4, a, f12.4)'), xpos, ' ', ypos, ' ', gas_pressure
          write(30, '(f12.4, a, f12.4, a, f12.4)'), xpos, ' ', ypos, ' ', mag_pressure
          write(40, '(f12.4, a, f12.4, a, f12.4)'), xpos, ' ', ypos, ' ', velc_squrd

       end do
    end do



  ! Closing files.
  ! --------------

    close (10)
    close (20)
    close (30)
    close (40)

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Output_data_B
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Output Data ------------------------------------------------------------------------------------------ [FLU10] ---
! ----------------------------------------------------------------------------------------------------------------------------------
! --- Version C: Mixed Output (magnetic field components: z- in 2D; transverse in 1D) ----------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Output_data_C (num)

    integer, intent (in) :: num


  ! Declaration of local variables.
  ! -------------------------------

    character :: magt_file*(19)
    character :: magz_file*(19)
    character :: chari*(4)

    real (PREC) :: xpos, ypos
    integer :: m, i, j, k



  ! Obtaining file names.
  ! ---------------------

    write(chari, '(i4.4)') num
    chari = adjustl(chari)

    magt_file = "output/magt"//chari//".txt"
    magz_file = "output/magz"//chari//".txt"



  ! Opening files.
  ! --------------

    open(10, file = magt_file, form = 'formatted')
    open(20, file = magz_file, form = 'formatted')



  ! Storing the 1D magnetic field file.
  ! -----------------------------------

    do m = 1, XSIZE

       xpos  = ((real(m) - 0.5D0)/real(XSIZE)) * XDOMN
       
       write(10, '(f12.4, a, f12.4)'), xpos, ' ', -(sin(alpha) * x_magfield(m, 1, 1)) + (cos(alpha) * y_magfield(m, 1, 1))
       
    end do
    


  ! Storing the 2D magnetic field file.
  ! -----------------------------------

    do j = 1, YSIZE
       do i = 1, XSIZE

          xpos = ((real(i) - 0.5D0)/real(XSIZE)) * XDOMN
          ypos = ((real(j) - 0.5D0)/real(YSIZE)) * YDOMN

          write(20, '(f12.4, a, f12.4, a, f12.4)'), xpos, ' ', ypos, ' ', z_magfield(i, j, 1)

       end do
    end do



  ! Closing files.
  ! --------------

    close (10)
    close (20)

    return
    


! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Output_data_C
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Save Restart Files  ---------------------------------------------------------------------------------- [FLU11] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Save_restart_files ()

    BASE_NTIME = BASE_NTIME + FULLTIME

    open (50, file = "output/restart1.bin", form = 'unformatted')
    write (50) BASE_FILE_NO, BASE_NTIME, DELTAT, FULLTIME, COURANT, MAX_DELTAT, XSIZE, &
         YSIZE, ZSIZE, DIMENSIONS, BOUNDARY, XDOMN, YDOMN, ZDOMN, MHDUNITS, BOUNDTYPE
    close (50)

    open (60, file = "output/restart2.bin", form = 'unformatted')
    write (60) density, pressure, x_velocity, y_velocity, z_velocity, gamma, &
         x_magfield, y_magfield, z_magfield
    close (60)

    return
    


! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Save_restart_files
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Deallocate the Fluid --------------------------------------------------------------------------------- [FLU12] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Deallocate_fluid ()


  ! Deallocating the fluid variables.
  ! ---------------------------------

    deallocate (density)
    deallocate (x_velocity)
    deallocate (y_velocity)
    deallocate (z_velocity)
    deallocate (x_momentum)
    deallocate (y_momentum)
    deallocate (z_momentum)
    deallocate (energy)
    deallocate (pressure)
    deallocate (gamma)



  ! Deallocating the magnetic field variables.
  ! ------------------------------------------

    deallocate (x_magfield)
    deallocate (y_magfield)
    deallocate (z_magfield)

    if (DIMENSIONS > 1) then

       deallocate (x_magfield_stored)
       deallocate (y_magfield_stored)
       deallocate (z_magfield_stored)
       
       deallocate (delta_bx)
       deallocate (delta_by)
       deallocate (delta_bz)

    end if
    
    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Deallocate_fluid
! ----------------------------------------------------------------------------------------------------------------------------------



end module fluid
