! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
! --- The Main Program -------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
! --- This code is best viewed with a window at least 135 characters wide. ---------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------



program magneto

  use start
  use setup
  use fluid
  use tests
  use grids
  use recon
  use solve
  use evolve
  use error
  implicit none


  ! Declaration of local variables.
  ! -------------------------------

  integer :: n, file_no
  real (PREC) :: percentage
  real (PREC) :: scratch

  character :: bksp*8
  do n = 1, 8
     bksp(n:n) = char(8)
  end do



  ! Initialising the simulation.
  ! ----------------------------

  call Check_initial_errors ()      ! --> <ERR01>
  call Initialise_simulation ()     ! --> <EVO01>

  if (OUTPUT_TYPE == 'A') call Output_data_A (-1, 1, 1, BASE_FILE_NO)     ! --> <FLU08>
  if (OUTPUT_TYPE == 'B') call Output_data_B (BASE_FILE_NO)               ! --> <FLU09>
  if (OUTPUT_TYPE == 'C') call Output_data_C (BASE_FILE_NO)               ! --> <FLU10>



! ----------------------------------------------------------------------------------------------------------------------------------
! Running simulation: variable DELTAT.
! ----------------------------------------------------------------------------------------------------------------------------------

  if (VARIABLE_DELTAT) then

     if (PRINT_DELTAT) open(50, file = "output/timestep.txt", form = 'formatted', position = 'append')

     write (*, '(a)', advance = 'no') 'Percentage complete...    0.00%'

     n = 0

     do while (NTIME < (BASE_NTIME + FULLTIME))

        ! Evolving the grid

        call Evolve_grid ()     ! --> <EVO02>
        n = n + 1


        ! Updating timestep

        if (PRINT_DELTAT) write (50, '(f12.4, a, f16.8)'), NTIME, ' ', DELTAT

        NTIME = NTIME + DELTAT
        call Determine_timestep ()     ! --> <FLU07>

        scratch = (BASE_NTIME + FULLTIME) - NTIME
        if (scratch < DELTAT) DELTAT = scratch


        ! Displaying percentage indicator

        percentage = ((NTIME - BASE_NTIME)/FULLTIME) * 100
        write (*, '(a, f7.3, a)', advance = 'no') bksp, percentage, '%'


        ! Writing output files

        if (mod(n, PRINT_FREQ) == 0) then

           file_no = BASE_FILE_NO + (n / PRINT_FREQ)

           if (OUTPUT_TYPE == 'A') call Output_data_A (-1, 1, 1, file_no)     ! --> <FLU08>
           if (OUTPUT_TYPE == 'B') call Output_data_B (file_no)               ! --> <FLU09>
           if (OUTPUT_TYPE == 'C') call Output_data_C (file_no)               ! --> <FLU10>

        end if

     end do

     write (*, '(a)') ' '

     if (PRINT_DELTAT) close (50)



! ----------------------------------------------------------------------------------------------------------------------------------
! Running simulation: fixed DELTAT.
! ----------------------------------------------------------------------------------------------------------------------------------

  else

     write (*, '(a)', advance = 'no') 'Percentage complete...    0.00%'

     do n = 1, NSTEPS

        ! Evolving the grid
        
        DELTAT = MAX_DELTAT
        call Evolve_grid ()     ! --> <EVO02>
        
        
        ! Displaying percentage indicator
        
        percentage = (real(n)/real(NSTEPS)) * 100
        write (*, '(a, f7.3, a)', advance = 'no') bksp, percentage, '%'
          
        
        ! Writing output files

        if (mod(n, PRINT_FREQ) == 0) then

           file_no = BASE_FILE_NO + (n / PRINT_FREQ)

           if (OUTPUT_TYPE == 'A') call Output_data_A (-1, 1, 1, file_no)     ! --> <FLU08>
           if (OUTPUT_TYPE == 'B') call Output_data_B (file_no)               ! --> <FLU09>
           if (OUTPUT_TYPE == 'C') call Output_data_C (file_no)               ! --> <FLU10>

        end if

     end do

     write (*, '(a)') ' '

  end if



! ----------------------------------------------------------------------------------------------------------------------------------

  ! Print final results.
  ! --------------------

  if (OUTPUT_TYPE == 'A') call Output_data_A (-1, 1, 1, 9999)     ! --> <FLU08>
  if (OUTPUT_TYPE == 'B') call Output_data_B (9999)               ! --> <FLU09>
  if (OUTPUT_TYPE == 'C') call Output_data_C (9999)               ! --> <FLU10>



  ! Shutting down simulation.
  ! -------------------------

  BASE_FILE_NO = file_no
  call Save_restart_files ()       ! --> <FLU11>
  call Shut_down_simulation ()     ! --> <EVO08>
          
  return

end program magneto
