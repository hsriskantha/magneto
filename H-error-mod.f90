! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.1 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
!
!     H: The Error Module 
!     ===================
!     -- Main actions: stop program for some errors.
!
! ----------------------------------------------------------------------------------------------------------------------------------
! --- This code is best viewed in a window at least 135 characters wide. -----------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------



module error

  use setup
  implicit none

contains



! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Check for Initial Errors ----------------------------------------------------------------------------- [ERR01] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Check_initial_errors ()


  ! Declaration of local variables.
  ! -------------------------------

    integer :: string_test
    logical :: error_flag, output_flag
    logical :: restart_flag1, restart_flag2

    error_flag = .false.



  ! Testing for output folder.
  ! --------------------------

    inquire(file = './output/.', exist = output_flag)
    if (output_flag .eqv. .false.) call system('mkdir output')
    


  ! Testing input in A_setup_mod
  ! ----------------------------

    ! TEST_PROBLEM

    string_test = 0
    string_test = verify(TEST_PROBLEM, 'ABCDEFX')

    if (string_test > 0) then
       write (*, *) 'ERROR: Invalid value of TEST_PROBLEM in A-setup-mod'
       write (*, *) ' --> See <SET_A> for more information.'
       write (*, *) ' '
       error_flag = .true.
    end if


    ! OUTPUT_TYPE

    string_test = 0
    string_test = verify(OUTPUT_TYPE, 'ABC')

    if (string_test > 0) then
       write (*, *) 'ERROR: Invalid value of OUTPUT_TYPE in A-setup-mod'
       write (*, *) ' --> See <SET_A> for more information.'
       write (*, *) ' '
       error_flag = .true.
    end if


    ! RECONSTRUCT_TYPE

    string_test = 0
    string_test = verify(RECONSTRUCT_TYPE, 'TPL')

    if (string_test > 0) then
       write (*, *) 'ERROR: Invalid value of RECONSTRUCT_TYPE in A-setup-mod'
       write (*, *) ' --> See <SET_A> for more information.'
       write (*, *) ' '
       error_flag = .true.
    end if


    ! WAVESPEED_TYPE

    string_test = 0
    string_test = verify(WAVESPEED_TYPE, 'MJ')

    if (string_test > 0) then
       write (*, *) 'ERROR: Invalid value of WAVESPEED_TYPE in A-setup-mod'
       write (*, *) ' --> See <SET_A> for more information.'
       write (*, *) ' '
       error_flag = .true.
    end if





  ! Testing presence of restart file.
  ! ---------------------------------

    if (TEST_PROBLEM == 'X') then

       inquire(file='./output/restart1.bin', exist = restart_flag1)
       inquire(file='./output/restart2.bin', exist = restart_flag2)

       if ((restart_flag1 .eqv. .false.) .or. (restart_flag2 .eqv. .false.)) then
          write (*, *) 'ERROR: You have chosen to use a restart file in A-setup-mod.'
          write (*, *) 'However, at least one of the restart files cannot be found.'
          write (*, *) 'Make sure restart1.bin and restart2.bin are saved in the'
          write (*, *) 'folder ./output/ or alternatively choose a new test problem.'
          write (*, *) ' '
          error_flag = .true.
       end if

    end if



  ! Terminating program.
  ! --------------------

    if (error_flag) stop

    return
     


! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Check_initial_errors
! ----------------------------------------------------------------------------------------------------------------------------------



end module error
