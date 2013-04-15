! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
!
!     H: The Error Module 
!     ===================
!     -- Main actions: stop program for some errors.
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
    


  ! Testing input in start.txt
  ! --------------------------

    ! TEST_PROBLEM

    string_test = 0
    string_test = verify(TEST_PROBLEM, 'ABCDEFX')

    if (string_test > 0) then
       write (*, *) 'ERROR: Invalid value of TEST_PROBLEM in start.txt'
       write (*, *) ' --> See <STA_A> for more information.'
       write (*, *) ' '
       error_flag = .true.
    end if


    ! OUTPUT_TYPE

    string_test = 0
    string_test = verify(OUTPUT_TYPE, 'ABC')

    if (string_test > 0) then
       write (*, *) 'ERROR: Invalid value of OUTPUT_TYPE in start.txt'
       write (*, *) ' --> See <STA_A> for more information.'
       write (*, *) ' '
       error_flag = .true.
    end if


    ! RECONSTRUCT_TYPE

    string_test = 0
    string_test = verify(RECONSTRUCT_TYPE, 'TPC')

    if (string_test > 0) then
       write (*, *) 'ERROR: Invalid value of RECONSTRUCT_TYPE in start.txt'
       write (*, *) ' --> See <STA_A> for more information.'
       write (*, *) ' '
       error_flag = .true.
    end if


    ! WAVESPEED_TYPE

    string_test = 0
    string_test = verify(WAVESPEED_TYPE, 'MJ')

    if (string_test > 0) then
       write (*, *) 'ERROR: Invalid value of WAVESPEED_TYPE in start.txt'
       write (*, *) ' --> See <STA_A> for more information.'
       write (*, *) ' '
       error_flag = .true.
    end if



  ! Testing presence of restart file.
  ! ---------------------------------

    if (TEST_PROBLEM == 'X') then

       inquire(file='./output/restart1.bin', exist = restart_flag1)
       inquire(file='./output/restart2.bin', exist = restart_flag2)

       if ((restart_flag1 .eqv. .false.) .or. (restart_flag2 .eqv. .false.)) then
          write (*, *) 'ERROR: You have chosen to use a restart file in start.txt.'
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





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Check for Errors in the Test Subroutines ------------------------------------------------------------- [ERR02] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Check_test_errors ()      


  ! Declaration of local variables.
  ! -------------------------------

    logical :: error_flag = .false.



  ! Testing boundary cells.
  ! -----------------------

    if (RECONSTRUCT_TYPE == 'C') then

       if (BOUNDARY < 4) then

          write (*, *) 'ERROR: Need at least four boundary cells if using the PPM'
          write (*, *) '(CS Edition) reconstruction. Adjust the variable BOUNDARY'
          write (*, *) 'in the appropriate <TESxx> subroutine.'
          write (*, *) ' '
          error_flag = .true.

       end if

    end if
 


  ! Terminating program.
  ! --------------------

    if (error_flag) stop

    return
     


! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Check_test_errors
! ----------------------------------------------------------------------------------------------------------------------------------



end module error
