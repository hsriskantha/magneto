! ----------------------------------------------------------------------------------------------------------------------------------
! --- Quarkwood Magneto 8.2 --------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------
! 
!     E: The Reconstruction Module
!     ============================
!     -- Main variables: the left and right Riemann input states for the state variables.
!     -- Main actions: determine these input states.
!
! ----------------------------------------------------------------------------------------------------------------------------------
! --- This code is best viewed with a window at least 135 characters wide ----------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------



module recon

  use setup
  use grids
  implicit none


  ! --------------------------------------------------------------------------------------------------------------------------------
  ! [REC_A] Guide: Cell-Centred vs Face-Centred Values
  ! --------------------------------------------------------------------------------------------------------------------------------
  !   In this module, we introduce face-centred arrays. In the literature, face-centred values are often denoted using '+1/2'. 
  !   For example, the face (i+1/2, j, k) lies between (i, j, k) and (i+1, j, k). 
  !
  !   We used starred notation, such that i* = i + 1/2.
  !
  !   In the code, we still use (i, j k) for both cell-centred and face-centred values; the meaning now depends on the context.
  !   (However, to reference the right element, we just drop the stars. So (i*, j, k) is found at element (i, j, k).)
  !           ___ ___ ___ ___ ___ ___ ___
  !     e.g. |___|___|___|___|___|___|___ ...
  !           
  !            1   2   3   4   5   6      <-- cell-centred values
  !              1   2   3   4   5   6    <-- face-centred values
  !
  !   Note that we ignore the first face (which is not used anyway) to preserve array length.  
  ! --------------------------------------------------------------------------------------------------------------------------------


  ! The Riemann problem input states.
  ! ---------------------------------

  real (PREC), dimension(:), allocatable :: density_L,     density_R
  real (PREC), dimension(:), allocatable :: x_velocity_L,  x_velocity_R
  real (PREC), dimension(:), allocatable :: y_velocity_L,  y_velocity_R
  real (PREC), dimension(:), allocatable :: z_velocity_L,  z_velocity_R
  real (PREC), dimension(:), allocatable :: x_momentum_L,  x_momentum_R
  real (PREC), dimension(:), allocatable :: y_momentum_L,  y_momentum_R
  real (PREC), dimension(:), allocatable :: z_momentum_L,  z_momentum_R
  real (PREC), dimension(:), allocatable :: energy_L,      energy_R
  real (PREC), dimension(:), allocatable :: pressure_L,    pressure_R
  real (PREC), dimension(:), allocatable :: gamma_L,       gamma_R
  real (PREC), dimension(:), allocatable :: sound_speed_L, sound_speed_R

  real (PREC), dimension(:), allocatable :: x_magfield_L,  x_magfield_R
  real (PREC), dimension(:), allocatable :: y_magfield_L,  y_magfield_R
  real (PREC), dimension(:), allocatable :: z_magfield_L,  z_magfield_R



  ! Working variables.
  ! ------------------

  real (PREC), dimension(:), allocatable   :: gas_pressure_1D
  real (PREC), dimension(:), allocatable   :: a, bx, by, bz
  real (PREC), dimension(:), allocatable   :: cf, ca, cs, cax
  
  real (PREC), dimension(:,:), allocatable :: lambda, w
  real (PREC), dimension(:,:,:), allocatable :: L, R
  
  real (PREC), dimension(:,:), allocatable :: w_hat_L, w_hat_R
  real (PREC), dimension(:,:), allocatable :: w_new_L, w_new_R
  real (PREC), dimension(:,:), allocatable :: delta_w
  real (PREC), dimension(:,:), allocatable :: w_L, w_R


  ! Specific for TVD

  real (PREC), dimension(:,:), allocatable :: wp, w0, wm
  real (PREC), dimension(:,:), allocatable :: rp_p, rm_p, rd_p
  real (PREC), dimension(:,:), allocatable :: rp_m, rm_m, rd_m
  real (PREC), dimension(:,:), allocatable :: delta_V
  real (PREC), dimension(:,:), allocatable :: sigma, dw


  ! Specific for PPM(L)

  real (PREC), dimension(:,:), allocatable :: wh
  
  real (PREC), dimension(:,:), allocatable :: dwL, dwR, dwC, dwG
  real (PREC), dimension(:,:), allocatable :: daL, daR, daC, daG
  
  real (PREC), dimension(:,:), allocatable :: PPM_dw, PPM_w6
  real (PREC), dimension(:,:), allocatable :: delta_a

contains



! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Initialise the Riemann Problem Input States ---------------------------------------------------------- [REC01] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Initialise_input_states ()


  ! Allocating the Riemann problem input states.
  ! --------------------------------------------

    allocate (density_L    (MFULL));   allocate (density_R    (MFULL))
    allocate (x_velocity_L (MFULL));   allocate (x_velocity_R (MFULL))
    allocate (y_velocity_L (MFULL));   allocate (y_velocity_R (MFULL))
    allocate (z_velocity_L (MFULL));   allocate (z_velocity_R (MFULL))
    allocate (x_momentum_L (MFULL));   allocate (x_momentum_R (MFULL))
    allocate (y_momentum_L (MFULL));   allocate (y_momentum_R (MFULL))
    allocate (z_momentum_L (MFULL));   allocate (z_momentum_R (MFULL))
    allocate (energy_L     (MFULL));   allocate (energy_R     (MFULL))
    allocate (pressure_L   (MFULL));   allocate (pressure_R   (MFULL))
    allocate (gamma_L      (MFULL));   allocate (gamma_R      (MFULL))
    allocate (sound_speed_L(MFULL));   allocate (sound_speed_R(MFULL))

    allocate (x_magfield_L (MFULL));   allocate (x_magfield_R (MFULL))
    allocate (y_magfield_L (MFULL));   allocate (y_magfield_R (MFULL))
    allocate (z_magfield_L (MFULL));   allocate (z_magfield_R (MFULL))



  ! Allocating the working variables.
  ! ---------------------------------

    allocate (gas_pressure_1D(MFULL))

    allocate (a(MFULL));    allocate (cax(MFULL))
    allocate (bx(MFULL));   allocate (by(MFULL));   allocate (bz(MFULL))
    allocate (cf(MFULL));   allocate (ca(MFULL));   allocate (cs(MFULL))

    allocate (w(7, MFULL));      allocate (lambda(7, MFULL))
    allocate (L(7, 7, MFULL));   allocate (R(7, 7, MFULL))

    allocate (delta_w(7, MFULL))

    allocate (w_hat_L(7, MFULL));   allocate (w_hat_R(7, MFULL))
    allocate (w_new_L(7, MFULL));   allocate (w_new_R(7, MFULL))
    allocate (w_L(7, MFULL));       allocate (w_R(7, MFULL))


    ! Specific for TVD

    if (RECONSTRUCT_TYPE == 'T') then

       allocate (wp(7, MFULL));   allocate (w0(7, MFULL));   allocate (wm(7, MFULL));

       allocate (rp_p(7, MFULL));   allocate (rm_p(7, MFULL));   allocate (rd_p(7, MFULL));
       allocate (rp_m(7, MFULL));   allocate (rm_m(7, MFULL));   allocate (rd_m(7, MFULL))

       allocate (delta_V(7, MFULL))
       allocate (sigma(7, MFULL))
       allocate (dw(7, MFULL))

    end if


    ! Specific for PPM(L)

    if ((RECONSTRUCT_TYPE == 'P') .or. (RECONSTRUCT_TYPE == 'L')) then

       allocate (wh(7, MFULL))

       allocate (dwL(7, MFULL));   allocate (dwR(7, MFULL));   allocate (dwC(7, MFULL));   allocate (dwG(7, MFULL))
       allocate (daL(7, MFULL));   allocate (daR(7, MFULL));   allocate (daC(7, MFULL));   allocate (daG(7, MFULL))
       
       allocate (PPM_dw(7, MFULL))   
       allocate (PPM_w6(7, MFULL))
       allocate (delta_a(7, MFULL))  

    end if

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Initialise_input_states
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Determine the Input States for the Riemann Problem (using Balsara's TVD Scheme) ---------------------- [REC02] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Determine_Riemann_input_TVD

  ! This subroutine uses the algorithm of Balsara, Astrophys. J. Suppl. S., 116, 133 (1998)


  ! Declaration of local variables.
  ! -------------------------------

    real (PREC), dimension(7) :: psi

    real (PREC), dimension(7) :: rpl_p, rmi_p, del_p
    real (PREC), dimension(7) :: rpl_m, rmi_m, del_m

    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd

    real (PREC) :: beta, mu, qa

    real (PREC) :: scratch1, scratch2

    integer :: m     ! m is used to label spatial position (i.e. m in [1, MFULL])
    integer :: p     ! p and q are used to label elements of state-based vectors (i.e. with 7 elements).  
    integer :: q     !   Usually, p labels eigenvalues (e.g. p = 1 corresponds to v_x - c_f).



  ! Step one: calculating eigens-system.
  ! ------------------------------------

    call Calculate_eigensystem ()



  ! Step two: project the state vectors onto the characteristic variables.
  ! ----------------------------------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

       do q = 1, 7

          wp(q, m) = 0.0D0
          w0(q, m) = 0.0D0
          wm(q, m) = 0.0D0

       end do

       do q = 1, 7
          do p = 1, 7

             wp(p, m) = wp(p, m) + (L(p, q, m) * w(q, m+1))
             w0(p, m) = w0(p, m) + (L(p, q, m) * w(q, m  ))
             wm(p, m) = wm(p, m) + (L(p, q, m) * w(q, m-1))

          end do
       end do

    end do



  ! Step three: apply monontonicity constraints.
  ! --------------------------------------------

    psi = [2.0D0, 1.0D0, 2.0D0, 1.0D0, 2.0D0, 1.0D0, 2.0D0]

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          delta_w(q, m) = 0.0D0

          if (((wp(q, m) - w0(q, m)) * (w0(q, m) - wm(q, m))) .gt. 0.0D0) then

             delta_w(q, m) = min (0.50D0 * abs(wp(q, m) - wm(q, m)), &
                                  psi(q) * abs(wp(q, m) - w0(q, m)), &
                                  psi(q) * abs(w0(q, m) - wm(q, m))  )

             delta_w(q, m) = delta_w(q, m) * sign(1.0D0, (wp(q, m) - w0(q, m)))

          end if

       end do
    end do



  ! Step four: project the monotonized differences back onto the primitive variables.
  ! ---------------------------------------------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

       do q = 1, 7

          delta_V(q, m) = 0.0D0

       end do

       do q = 1, 7
          do p = 1, 7

             delta_V(p, m) = delta_V(p, m) + (delta_w(q, m) * R(p, q, m))

          end do
       end do

    end do



  ! Step five: steepening algorithm.
  ! --------------------------------

    ! Calculating sigma

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          beta = 0.0D0
          mu   = 0.0D0

          sigma(q, m) = 0.0D0


          if (psi(q) .lt. 1.5D0) then     ! i.e. only linearly degenerate characteristic fields

             if (((wp(q, m) - w0(q, m)) * (w0(q, m) - wm(q, m))) .gt. 0.0D0) then

                scratch1 = abs((w0(q, m) - wm(q, m)) / (wp(q, m) - w0(q, m)))
                scratch2 = 4.3D0 * (lambda(q, m) * dtdx(m)) - 0.5D0 * sign(1.0D0, lambda(q, m))

                beta = scratch1**scratch2

             end if

             scratch1 = wp(q, m) - 2.0D0 * w0(q, m) + wm(q, m)
             scratch2 = abs(w0(q, m) - wm(q, m)) + abs(wp(q, m) - w0(q, m))

             if (abs(scratch2) .gt. 0.0D0) mu = (scratch1/scratch2)**2.0D0

             sigma(q, m) = 33.0D0 * mu * beta

          end if

       end do
    end do


    ! Projections

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

       do q = 1, 7

          rp_p(q, m) = 0.0D0
          rm_p(q, m) = 0.0D0
          rd_p(q, m) = 0.0D0

          rp_m(q, m) = 0.0D0
          rm_m(q, m) = 0.0D0
          rd_m(q, m) = 0.0D0


          rpl_p(q) = w(q, m+1) - 0.5D0 * delta_V(q, m+1)
          rpl_m(q) = w(q, m  ) - 0.5D0 * delta_V(q, m  )

          rmi_p(q) = w(q, m  ) + 0.5D0 * delta_V(q, m  )
          rmi_m(q) = w(q, m-1) + 0.5D0 * delta_V(q, m-1)

          del_p(q) = rpl_p(q) - rmi_p(q)
          del_m(q) = rpl_m(q) - rmi_m(q)

       end do

       do q = 1, 7
          do p = 1, 7

             rp_p(p, m) = rp_p(p, m) + (L(p, q, m) * rpl_p(q))
             rm_p(p, m) = rm_p(p, m) + (L(p, q, m) * rmi_p(q))
             rd_p(p, m) = rd_p(p, m) + (L(p, q, m) * del_p(q))

             rp_m(p, m) = rp_m(p, m) + (L(p, q, m) * rpl_m(q))
             rm_m(p, m) = rm_m(p, m) + (L(p, q, m) * rmi_m(q))
             rd_m(p, m) = rd_m(p, m) + (L(p, q, m) * del_m(q))

          end do
       end do

    end do


    ! Slope modifier

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          dw(q, m) = 0.0D0

          if (psi(q) .lt. 1.5D0) then

             scratch1 = minmod(rd_m(q, m), rd_p(q, m))

             scratch2 = minmod( (wp(q, m) - rm_p(q, m)),  (rp_m(q, m-1) - wm(q, m)) )

             dw(q, m) = 2.0D0 * minmod(sigma(q, m) * scratch1, scratch2)

          end if

          delta_w(q, m) = delta_w(q, m) + dw(q, m)

       end do
    end do


    ! Applying steepening

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       
       do q = 1, 7
          do p = 1, 7
             
             delta_V(p, m) = delta_V(p, m) + (dw(q, m) * R(p, q, m))
             
          end do
       end do
       
    end do



  ! Step six: compute the left and right interface values.
  ! ------------------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

        scratch1 = 0.5D0 * max(lambda(7, m  ), 0.0D0) * dtdx(m)
        scratch2 = 0.5D0 * min(lambda(1, m+1), 0.0D0) * dtdx(m)

        do q = 1, 7

           w_hat_L(q, m) = w(q, m  ) + (0.5D0 - scratch1) * delta_V(q, m)
           w_hat_R(q, m) = w(q, m+1) - (0.5D0 + scratch2) * delta_V(q, m+1)

        end do

    end do



  ! Step seven: the characteristic tracing.
  ! ---------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          w_new_L(q, m) = w_hat_L(q, m)
          w_new_R(q, m) = w_hat_R(q, m)

       end do
    end do


    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do p = 1, 7

          if (lambda(p, m) .gt. 0.0D0) then

             qa = 0.5D0 * dtdx(m) * (lambda(7, m) - lambda(p, m))

             do q = 1, 7

                w_new_L(q, m) = w_new_L(q, m) + qa * delta_w(p, m) * R(q, p, m)

             end do

          end if


          if (lambda(p, m+1) .lt. 0.0D0) then

             qa = 0.5D0 * dtdx(m+1) * (lambda(1, m+1) - lambda(p, m+1))

             do q = 1, 7

                w_new_R(q, m) = w_new_R(q, m) + qa * delta_w(p, m+1) * R(q, p, m+1)

             end do

          end if

       end do
    end do



  ! Step eight: the final left and right states.
  ! --------------------------------------------
  
    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

       density_L(m)    = w_new_L(1, m)
       x_velocity_L(m) = w_new_L(2, m)
       y_velocity_L(m) = w_new_L(3, m)
       z_velocity_L(m) = w_new_L(4, m)
       x_momentum_L(m) = w_new_L(1, m) * w_new_L(2, m)
       y_momentum_L(m) = w_new_L(1, m) * w_new_L(3, m)
       z_momentum_L(m) = w_new_L(1, m) * w_new_L(4, m)
       gamma_L(m)      = gamma_1D(m)

       x_magfield_L(m) = 0.5D0 * (x_magfield_1D(m) + x_magfield_1D(m+1))
       y_magfield_L(m) = w_new_L(6, m) * MHDF(-1)
       z_magfield_L(m) = w_new_L(7, m) * MHDF(-1)


       magf = Create_vector (x_magfield_L(m), y_magfield_L(m), z_magfield_L(m))
       velc = Create_vector (x_velocity_L(m), y_velocity_L(m), z_velocity_L(m))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure_L(m)   = Calculate_total_pressure (w_new_L(5, m), magf_squrd)
       energy_L(m)     = Calculate_energy_EOS (density_L(m), pressure_L(m), gamma_L(m), velc_squrd, magf_squrd)


       density_R(m)    = w_new_R(1, m)
       x_velocity_R(m) = w_new_R(2, m)
       y_velocity_R(m) = w_new_R(3, m)
       z_velocity_R(m) = w_new_R(4, m)
       x_momentum_R(m) = w_new_R(1, m) * w_new_R(2, m)
       y_momentum_R(m) = w_new_R(1, m) * w_new_R(3, m)
       z_momentum_R(m) = w_new_R(1, m) * w_new_R(4, m)
       gamma_R(m)      = gamma_1D(m+1)

       x_magfield_R(m) = x_magfield_L(m)
       y_magfield_R(m) = w_new_R(6, m) * MHDF(-1)
       z_magfield_R(m) = w_new_R(7, m) * MHDF(-1)


       magf = Create_vector (x_magfield_R(m), y_magfield_R(m), z_magfield_R(m))
       velc = Create_vector (x_velocity_R(m), y_velocity_R(m), z_velocity_R(m))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure_R(m)   = Calculate_total_pressure (w_new_R(5, m), magf_squrd)
       energy_R(m)     = Calculate_energy_EOS (density_R(m), pressure_R(m), gamma_R(m), velc_squrd, magf_squrd)

    end do

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Determine_Riemann_input_TVD
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Determine the Input States for the Riemann Problem (using Piecewise Parabolic Method) ------------------ [REC03] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Determine_Riemann_input_PPM ()

  ! This subroutine uses the PPM algorithm of Collela and Woodward, J. Comput. Phys., 54, 174 (1984)
  !   --> the characteristic tracing is based on S4.2.3 of Stone et al., Astrophys. J. Suppl. S., 178, 137 (2008)


  ! Declaration of local variables.
  ! -------------------------------

    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd

    logical :: cond1, cond2, cond3, cond4

    real (PREC) :: A, B, C
    real (PREC) :: scratch1, scratch2

    integer :: m     ! m is used to label spatial position (i.e. m in [1, MFULL])
    integer :: p     ! p and q are used to label elements of state-based vectors (i.e. with 7 elements).  
    integer :: q     !   Usually, p labels eigenvalues (e.g. p = 1 corresponds to v_x - c_f).



  ! Step one: calculating eigen-system.
  ! -----------------------------------

    call Calculate_eigensystem ()



  ! Step two: compute left, right and centred differences.
  ! ------------------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          dwG(q, m) = 0.0D0

          dwL(q, m) =  w(q, m  ) - w(q, m-1)
          dwR(q, m) =  w(q, m+1) - w(q, m  )
          dwC(q, m) = (w(q, m+1) - w(q, m-1))/2.0D0

          if ((dwL(q, m) * dwR(q, m)) > 0.0D0) then 
             dwG(q, m) = 2.0D0 * dwL(q, m) * dwR(q, m) / (dwL(q, m) + dwR(q, m))
          end if

       end do
    end do



  ! Step three: project differences onto characteristic variables.
  ! --------------------------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
    
       do q = 1, 7

          daL(q, m) = 0.0D0
          daR(q, m) = 0.0D0
          daC(q, m) = 0.0D0
          daG(q, m) = 0.0D0

       end do

       do q = 1, 7
          do p = 1, 7

             daL(p, m) = daL(p, m) + (L(p, q, m) * dwL(q, m))
             daR(p, m) = daR(p, m) + (L(p, q, m) * dwR(q, m))
             daC(p, m) = daC(p, m) + (L(p, q, m) * dwC(q, m))
             daG(p, m) = daG(p, m) + (L(p, q, m) * dwG(q, m))

          end do
       end do

    end do



  ! Step four: apply monotonicity constraints.
  ! ------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          delta_a(q, m) = 0.0D0

          if ((daL(q, m) * daR(q, m)) >= 0.0D0) then

             scratch1 = min(abs(daL(q, m)), abs(daR(q, m)))
             scratch2 = min(abs(daC(q, m)), abs(daG(q, m)))

             scratch2 = abs(daC(q, m))
             
             delta_a(q, m) = sign(1.0D0, daC(q, m)) * min(2.0D0 * scratch1, scratch2)

          end if

       end do
    end do



  ! Step five: project back onto the primitive variables.
  ! -----------------------------------------------------
    

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

       do q = 1, 7

          delta_w(q, m) = 0.0D0

       end do

       do q = 1, 7
          do p = 1, 7

             delta_w(p, m) = delta_w(p, m) + (delta_a(q, m) * R(p, q, m))

          end do
       end do

    end do



  ! Step six: parabolic interpolation.
  ! ----------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          w_L(q, m) = ((w(q, m  ) + w(q, m-1))/2.0D0) - ((delta_w(q, m  ) - delta_w(q, m-1))/6.0D0)
          w_R(q, m) = ((w(q, m+1) + w(q, m  ))/2.0D0) - ((delta_w(q, m+1) - delta_w(q, m  ))/6.0D0)

       end do
    end do



  ! Step seven: further monotonicity constraints.
  ! ---------------------------------------------

    ! Taken from the Athena (Fortran) source code

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          w_L(q, m) = max(min(w(q, m), w(q, m-1)), w_L(q, m))
          w_L(q, m) = min(max(w(q, m), w(q, m-1)), w_L(q, m))
          w_R(q, m) = max(min(w(q, m), w(q, m+1)), w_R(q, m))
          w_R(q, m) = min(max(w(q, m), w(q, m+1)), w_R(q, m))

       end do
    end do


    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          if (((w_R(q, m) - w(q, m)) * (w(q, m) - w_L(q, m))) <= 0.0D0) then

             w_L(q, m) = w(q, m)
             w_R(q, m) = w(q, m)
             
          end if
          
          
          scratch1 = (w_R(q, m) - w_L(q, m))
          scratch2 = (w_R(q, m) + w_L(q, m)) / 2.0D0
          
          if ((6.0D0 * scratch1) * (w(q, m) - scratch2) > scratch1**2.0D0) then
             
             w_L(q, m) = (3.0D0 * w(q, m)) - (2.0D0 * w_R(q, m))
             
          else if ((6.0D0 * scratch1) * (w(q, m) - scratch2) < -(scratch1)**2.0D0) then
             
             w_R(q, m) = (3.0D0 * w(q, m)) - (2.0D0 * w_L(q, m))
             
          end if

       end do
    end do



  ! Step eight: compute PPM coefficients.
  ! -------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          PPM_dw(q, m) = w_R(q, m) - w_L(q, m)
          PPM_w6(q, m) = (6.0D0 * w(q, m)) - 3.0D0 * (w_L(q, m) + w_R(q, m))

       end do
    end do



  ! Step nine: compute the left and right interface values.
  ! -------------------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

        scratch1 =  max(lambda(7, m), 0.0D0) * dtdx(m)
        scratch2 = -min(lambda(1, m), 0.0D0) * dtdx(m)

        do q = 1, 7

           w_hat_L(q, m) = w_R(q, m) - (0.5D0 * scratch1) * (PPM_dw(q, m) - ((1.0D0 - (scratch1*(2.0D0/3.0D0))) * PPM_w6(q, m)))
           w_hat_R(q, m) = w_L(q, m) + (0.5D0 * scratch2) * (PPM_dw(q, m) + ((1.0D0 - (scratch2*(2.0D0/3.0D0))) * PPM_w6(q, m)))

        end do

    end do



  ! Step ten: the characteristic tracing.
  ! -------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          w_new_L(q, m) = w_hat_L(q, m)
          w_new_R(q, m) = w_hat_R(q, m+1)

       end do
    end do


    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do p = 1, 7

          if (lambda(p, m) > 0.0D0) then

             C = 0.0D0

             A = 0.5D0 * dtdx(m) * (lambda(7, m) - lambda(p, m))
             B = (1.0D0/3.0D0) * dtdx(m)**2.0D0 * (lambda(7, m)**2.0D0 - lambda(p, m)**2.0D0)

             do q = 1, 7
                C = C + L(p, q, m) * (A * (PPM_dw(q, m) - PPM_w6(q, m)) + (B * (PPM_w6(q, m))))
             end do

             do q = 1, 7
                w_new_L(q, m) = w_new_L(q, m) + C * R(q, p, m)
             end do

          end if


          if (lambda(p, m) < 0.0D0) then

             C = 0.0D0

             A = 0.5D0 * dtdx(m) * (lambda(1, m) - lambda(p, m))
             B = (1.0D0/3.0D0) * dtdx(m)**2.0D0 * (lambda(1, m)**2.0D0 - lambda(p, m)**2.0D0)

             do q = 1, 7
                 C = C + L(p, q, m) * (A * (PPM_dw(q, m) + PPM_w6(q, m)) + (B * (PPM_w6(q, m))))
             end do

             do q = 1, 7
               w_new_R(q, m) = w_new_R(q, m) + C * R(q, p, m)
             end do

          end if

       end do
    end do     



  ! Step eleven: the final left and right states.
  ! ---------------------------------------------                              

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

       density_L(m)    = w_new_L(1, m)
       x_velocity_L(m) = w_new_L(2, m)
       y_velocity_L(m) = w_new_L(3, m)
       z_velocity_L(m) = w_new_L(4, m)
       x_momentum_L(m) = w_new_L(1, m) * w_new_L(2, m)
       y_momentum_L(m) = w_new_L(1, m) * w_new_L(3, m)
       z_momentum_L(m) = w_new_L(1, m) * w_new_L(4, m)
       gamma_L(m)      = gamma_1D(m)

       x_magfield_L(m) = 0.5D0 * (x_magfield_1D(m) + x_magfield_1D(m+1))
       y_magfield_L(m) = w_new_L(6, m) * MHDF(-1)
       z_magfield_L(m) = w_new_L(7, m) * MHDF(-1)


       magf = Create_vector (x_magfield_L(m), y_magfield_L(m), z_magfield_L(m))
       velc = Create_vector (x_velocity_L(m), y_velocity_L(m), z_velocity_L(m))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure_L(m)   = w_new_L(5, m) + (0.5D0 * magf_squrd * MHDF(2))
       energy_L(m)     = Calculate_energy_EOS (density_L(m), pressure_L(m), gamma_L(m), velc_squrd, magf_squrd)


       density_R(m)    = w_new_R(1, m)
       x_velocity_R(m) = w_new_R(2, m)
       y_velocity_R(m) = w_new_R(3, m)
       z_velocity_R(m) = w_new_R(4, m)
       x_momentum_R(m) = w_new_R(1, m) * w_new_R(2, m)
       y_momentum_R(m) = w_new_R(1, m) * w_new_R(3, m)
       z_momentum_R(m) = w_new_R(1, m) * w_new_R(4, m)
       gamma_R(m)      = gamma_1D(m+1)

       x_magfield_R(m) = x_magfield_L(m)
       y_magfield_R(m) = w_new_R(6, m) * MHDF(-1)
       z_magfield_R(m) = w_new_R(7, m) * MHDF(-1)


       magf = Create_vector (x_magfield_R(m), y_magfield_R(m), z_magfield_R(m))
       velc = Create_vector (x_velocity_R(m), y_velocity_R(m), z_velocity_R(m))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure_R(m)   = w_new_R(5, m) + (0.5D0 * magf_squrd * MHDF(2))
       energy_R(m)     = Calculate_energy_EOS (density_R(m), pressure_R(m), gamma_R(m), velc_squrd, magf_squrd)

    end do

    return
    


! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Determine_Riemann_input_PPM
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: Determine the Input States for the Riemann Problem (using Piecewise Parabolic Method with CS limiter)--- [REC04] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Determine_Riemann_input_PPML ()

  ! This subroutine uses the PPM algorithm of Collela and Woodward, J. Comput. Phys., 54, 174 (1984)
  !   --> ... and the limiter of Colella and Sekora, J. Comput. Phys., 227, 7069 (2008)


  ! Declaration of local variables.
  ! -------------------------------

    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd

    logical :: cond1, cond2, cond3, cond4

    real (PREC) :: A, B, C
    real (PREC) :: scratch1, scratch2, wlim

    integer :: m     ! m is used to label spatial position (i.e. m in [1, MFULL])
    integer :: p     ! p and q are used to label elements of state-based vectors (i.e. with 7 elements).  
    integer :: q     !   Usually, p labels eigenvalues (e.g. p = 1 corresponds to v_x - c_f).

    real (PREC), parameter :: epsilon = 1D-10
    real (PREC), parameter :: Clim = 1.25D0



  ! Step one: calculating eigen-system.
  ! -----------------------------------

    call Calculate_eigensystem ()



  ! Step two: compute left, right and centred differences.
  ! ------------------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          dwG(q, m) = 0.0D0

          dwL(q, m) =  w(q, m  ) - w(q, m-1)
          dwR(q, m) =  w(q, m+1) - w(q, m  )
          dwC(q, m) = (w(q, m+1) - w(q, m-1))/2.0D0

          if ((dwL(q, m) * dwR(q, m)) > 0.0D0) then 
             dwG(q, m) = 2.0D0 * dwL(q, m) * dwR(q, m) / (dwL(q, m) + dwR(q, m))
          end if

       end do
    end do



  ! Step three: project differences onto characteristic variables.
  ! --------------------------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
    
       do q = 1, 7

          daL(q, m) = 0.0D0
          daR(q, m) = 0.0D0
          daC(q, m) = 0.0D0
          daG(q, m) = 0.0D0

       end do

       do q = 1, 7
          do p = 1, 7

             daL(p, m) = daL(p, m) + (L(p, q, m) * dwL(q, m))
             daR(p, m) = daR(p, m) + (L(p, q, m) * dwR(q, m))
             daC(p, m) = daC(p, m) + (L(p, q, m) * dwC(q, m))
             daG(p, m) = daG(p, m) + (L(p, q, m) * dwG(q, m))

          end do
       end do

    end do



  ! Step four: apply monotonicity constraints.
  ! ------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          delta_a(q, m) = 0.0D0

          if ((daL(q, m) * daR(q, m)) >= 0.0D0) then

             scratch1 = min(abs(daL(q, m)), abs(daR(q, m)))
             scratch2 = min(abs(daC(q, m)), abs(daG(q, m)))

             scratch2 = abs(daC(q, m))
             
             delta_a(q, m) = sign(1.0D0, daC(q, m)) * min(2.0D0 * scratch1, scratch2)

          end if

       end do
    end do



  ! Step five: project back onto the primitive variables.
  ! -----------------------------------------------------
    
    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

       do q = 1, 7

          delta_w(q, m) = 0.0D0

       end do

       do q = 1, 7
          do p = 1, 7

             delta_w(p, m) = delta_w(p, m) + (delta_a(q, m) * R(p, q, m))

          end do
       end do

    end do



  ! Step six: parabolic interpolation.
  ! ----------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          w_L(q, m) = ((w(q, m  ) + w(q, m-1))/2.0D0) - ((delta_w(q, m  ) - delta_w(q, m-1))/6.0D0)
          w_R(q, m) = ((w(q, m+1) + w(q, m  ))/2.0D0) - ((delta_w(q, m+1) - delta_w(q, m  ))/6.0D0)

       end do
    end do



  ! Step seven: further monotonicity constraints.
  ! ---------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          cond1 = ((w_R(q, m) - w(q, m)) * (w(q, m) - w_L(q, m)) <= 0.0D0)
          cond2 = ((w(q, m-1) - w(q, m)) * (w(q, m) - w(q, m+1)) <= 0.0D0)

          if (cond1 .or. cond2) then

             PPM_w6(q, m) = (6.0D0 * w(q, m)) - 3.0D0 * (w_L(q, m) + w_R(q, m))

             daG(q, m) = -2.0D0 * PPM_w6(q, m)

             daC(q, m) = w(q, m-1) - 2.0D0 * w(q, m  ) + w(q, m+1)
             daL(q, m) = w(q, m-2) - 2.0D0 * w(q, m-1) + w(q, m  )
             daR(q, m) = w(q, m  ) - 2.0D0 * w(q, m+1) + w(q, m+2)

             cond3 = ((daG(q, m) < 0.0D0) .and. (daC(q, m) < 0.0D0) .and. (daL(q, m) < 0.0D0) .and. (daR(q, m) < 0.0D0))
             cond4 = ((daG(q, m) > 0.0D0) .and. (daC(q, m) > 0.0D0) .and. (daL(q, m) > 0.0D0) .and. (daR(q, m) > 0.0D0))

             wlim = 0.0D0

             if (cond3 .or. cond4) then

                scratch1 = min(abs(daL(q, m)), abs(daR(q, m)), abs(daC(q, m)))

                wlim = sign(1.0D0, daG(q, m)) * min(Clim * scratch1, abs(daG(q, m)))

             end if

             if (abs(daG(q, m)) < epsilon) then

                w_L(q, m) = w(q, m)
                w_R(q, m) = w(q, m)

             else

                w_L(q, m) = w(q, m) + (w_L(q, m) - w(q, m)) * (wlim/daG(q, m))
                w_R(q, m) = w(q, m) + (w_R(q, m) - w(q, m)) * (wlim/daG(q, m))
                
             end if

          else


             scratch1 = (w_R(q, m) - w_L(q, m))
             scratch2 = (w_R(q, m) + w_L(q, m)) / 2.0D0
             
             if ((6.0D0 * scratch1) * (w(q, m) - scratch2) > scratch1**2.0D0) then
                
                w_L(q, m) = (3.0D0 * w(q, m)) - (2.0D0 * w_R(q, m))
                
             else if ((6.0D0 * scratch1) * (w(q, m) - scratch2) < -(scratch1)**2.0D0) then
                
                w_R(q, m) = (3.0D0 * w(q, m)) - (2.0D0 * w_L(q, m))
                
             end if

          end if


       end do
    end do



  ! Step eight: compute PPM coefficients.
  ! -------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          PPM_dw(q, m) = w_R(q, m) - w_L(q, m)
          PPM_w6(q, m) = (6.0D0 * w(q, m)) - 3.0D0 * (w_L(q, m) + w_R(q, m))

       end do
    end do



  ! Step nine: compute the left and right interface values.
  ! -------------------------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

        scratch1 =  max(lambda(7, m), 0.0D0) * dtdx(m)
        scratch2 = -min(lambda(1, m), 0.0D0) * dtdx(m)

        do q = 1, 7

           w_hat_L(q, m) = w_R(q, m) - (0.5D0 * scratch1) * (PPM_dw(q, m) - ((1.0D0 - (scratch1*(2.0D0/3.0D0))) * PPM_w6(q, m)))
           w_hat_R(q, m) = w_L(q, m) + (0.5D0 * scratch2) * (PPM_dw(q, m) + ((1.0D0 - (scratch2*(2.0D0/3.0D0))) * PPM_w6(q, m)))

        end do

    end do



  ! Step ten: the characteristic tracing.
  ! -------------------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do q = 1, 7

          w_new_L(q, m) = w_hat_L(q, m)
          w_new_R(q, m) = w_hat_R(q, m+1)

       end do
    end do


    do m = BOUNDARY, (BOUNDARY + rowsize) + 1
       do p = 1, 7

          if (lambda(p, m) > 0.0D0) then

             C = 0.0D0

             A = 0.5D0 * dtdx(m) * (lambda(7, m) - lambda(p, m))
             B = (1.0D0/3.0D0) * dtdx(m)**2.0D0 * (lambda(7, m)**2.0D0 - lambda(p, m)**2.0D0)

             do q = 1, 7
                C = C + L(p, q, m) * (A * (PPM_dw(q, m) - PPM_w6(q, m)) + (B * (PPM_w6(q, m))))
             end do

             do q = 1, 7
                w_new_L(q, m) = w_new_L(q, m) + C * R(q, p, m)
             end do

          end if


          if (lambda(p, m) < 0.0D0) then

             C = 0.0D0

             A = 0.5D0 * dtdx(m) * (lambda(1, m) - lambda(p, m))
             B = (1.0D0/3.0D0) * dtdx(m)**2.0D0 * (lambda(1, m)**2.0D0 - lambda(p, m)**2.0D0)

             do q = 1, 7
                 C = C + L(p, q, m) * (A * (PPM_dw(q, m) + PPM_w6(q, m)) + (B * (PPM_w6(q, m))))
             end do

             do q = 1, 7
               w_new_R(q, m) = w_new_R(q, m) + C * R(q, p, m)
             end do

          end if

       end do
    end do     



  ! Step eleven: the final left and right states.
  ! ---------------------------------------------                              

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

       density_L(m)    = w_new_L(1, m)
       x_velocity_L(m) = w_new_L(2, m)
       y_velocity_L(m) = w_new_L(3, m)
       z_velocity_L(m) = w_new_L(4, m)
       x_momentum_L(m) = w_new_L(1, m) * w_new_L(2, m)
       y_momentum_L(m) = w_new_L(1, m) * w_new_L(3, m)
       z_momentum_L(m) = w_new_L(1, m) * w_new_L(4, m)
       gamma_L(m)      = gamma_1D(m)

       x_magfield_L(m) = 0.5D0 * (x_magfield_1D(m) + x_magfield_1D(m+1))
       y_magfield_L(m) = w_new_L(6, m) * MHDF(-1)
       z_magfield_L(m) = w_new_L(7, m) * MHDF(-1)


       magf = Create_vector (x_magfield_L(m), y_magfield_L(m), z_magfield_L(m))
       velc = Create_vector (x_velocity_L(m), y_velocity_L(m), z_velocity_L(m))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure_L(m)   = w_new_L(5, m) + (0.5D0 * magf_squrd * MHDF(2))
       energy_L(m)     = Calculate_energy_EOS (density_L(m), pressure_L(m), gamma_L(m), velc_squrd, magf_squrd)


       density_R(m)    = w_new_R(1, m)
       x_velocity_R(m) = w_new_R(2, m)
       y_velocity_R(m) = w_new_R(3, m)
       z_velocity_R(m) = w_new_R(4, m)
       x_momentum_R(m) = w_new_R(1, m) * w_new_R(2, m)
       y_momentum_R(m) = w_new_R(1, m) * w_new_R(3, m)
       z_momentum_R(m) = w_new_R(1, m) * w_new_R(4, m)
       gamma_R(m)      = gamma_1D(m+1)

       x_magfield_R(m) = x_magfield_L(m)
       y_magfield_R(m) = w_new_R(6, m) * MHDF(-1)
       z_magfield_R(m) = w_new_R(7, m) * MHDF(-1)


       magf = Create_vector (x_magfield_R(m), y_magfield_R(m), z_magfield_R(m))
       velc = Create_vector (x_velocity_R(m), y_velocity_R(m), z_velocity_R(m))

       magf_squrd = Dotproduct (magf, magf)
       velc_squrd = Dotproduct (velc, velc)

       pressure_R(m)   = w_new_R(5, m) + (0.5D0 * magf_squrd * MHDF(2))
       energy_R(m)     = Calculate_energy_EOS (density_R(m), pressure_R(m), gamma_R(m), velc_squrd, magf_squrd)

    end do

    return


            
! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Determine_Riemann_input_PPML
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- FUNCTION: General Minmod Function -------------------------------------------------------------------------------- [REC05] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  function minmod (x, y)

    real (PREC), intent (in) :: x, y
    real (PREC) :: minmod

    if (x * y .gt. 0.0D0) then

       minmod = sign(1.0D0, x) * min(abs(x), abs(y))

    else

       minmod = 0.0D0
       
    end if

    return

! ----------------------------------------------------------------------------------------------------------------------------------
  end function minmod
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Calculate the Eigen-system --------------------------------------------------------------------------- [REC06] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Calculate_eigensystem ()

    ! See Roe and Balsara, SIAM J. Appl. Math., 56, 57 (1996)


  ! Declaration of variables.
  ! -------------------------

    real (PREC), dimension(7) :: wt

    real (PREC), dimension(3) :: magf, velc
    real (PREC) :: magf_squrd, velc_squrd

    real (PREC) :: scratch1, scratch2, phi, S, asq
    real (PREC) :: alpha_f, alpha_s, beta_y, beta_z

    integer :: m     ! m is used to label spatial position (i.e. m in [1, MFULL])
    integer :: p     ! p and q are used to label elements of state-based vectors (i.e. with 7 elements).  
    integer :: q     !   Usually, p labels eigenvalues (e.g. p = 1 corresponds to v_x - c_f).



  ! Some initial calculations.
  ! --------------------------

    do m = 1, rowfull

       ! Calculating the gas pressure

       magf = Create_vector (x_magfield_1D(m), y_magfield_1D(m), z_magfield_1D(m))
       magf_squrd = Dotproduct (magf, magf)

       gas_pressure_1D(m) = Calculate_gas_pressure (pressure_1D(m), magf_squrd)

       scratch1 = dsqrt(gas_pressure_1D(m))


       ! Storing the state vector

       wt = Create_vector7 (density_1D(m), x_velocity_1D(m), y_velocity_1D(m), z_velocity_1D(m), gas_pressure_1D(m), &
                            (y_magfield_1D(m)*MHDF(1)), (z_magfield_1D(m)*MHDF(1)))

       do q = 1, 7
          w(q, m) = wt(q)
       end do

    end do


    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

       ! Calculating wavespeeds

       a(m) = dsqrt((gamma_1D(m) * gas_pressure_1D(m)) / density_1D(m))

       bx(m) = x_magfield_1D(m) * MHDF(1)
       by(m) = y_magfield_1D(m) * MHDF(1)
       bz(m) = z_magfield_1D(m) * MHDF(1)

       ca(m)  = dsqrt((bx(m)**2.0D0 + by(m)**2.0D0 + bz(m)**2.0D0) / density_1D(m))
       cax(m) = dsqrt(bx(m)**2.0D0 / density_1D(m))

       scratch1 = a(m)**2.0D0 + ca(m)**2.0D0
       scratch2 = 4.0D0 * a(m)**2.0D0 * cax(m)**2.0D0

       cf(m) = dsqrt(0.5D0 * (scratch1 + dsqrt(scratch1**2.0D0 - scratch2)))
       cs(m) = dsqrt(0.5D0 * (scratch1 - dsqrt(scratch1**2.0D0 - scratch2)))

    end do



  ! Calculating eigens-system.
  ! --------------------------

    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

            !(p, m)
       lambda(1, m) = x_velocity_1D(m) - cf(m)
       lambda(2, m) = x_velocity_1D(m) - cax(m)
       lambda(3, m) = x_velocity_1D(m) - cs(m)
       lambda(4, m) = x_velocity_1D(m)
       lambda(5, m) = x_velocity_1D(m) + cs(m)
       lambda(6, m) = x_velocity_1D(m) + cax(m)
       lambda(7, m) = x_velocity_1D(m) + cf(m)

    end do


    do m = BOUNDARY, (BOUNDARY + rowsize) + 1

       ! First, some constants

       scratch1 = cf(m)**2.0D0 - cs(m)**2.0D0

       if (scratch1 == 0.0D0) then

          phi = 0.5D0 * atan((ca(m) - cax(m)) / (abs(bx(m)) - a(m)))

          alpha_f = sin(phi)
          alpha_s = cos(phi)

       else if ((a(m)**2.0D0 - cs(m)**2.0D0) .le. 0.0D0) then

          alpha_f = 0.0D0
          alpha_s = dsqrt((cf(m)**2.0D0 - a(m)**2.0D0)/scratch1) 

       else if ((cf(m)**2.0D0 - a(m)**2.0D0) .le. 0.0D0) then

          alpha_f = dsqrt((a(m)**2.0D0 - cs(m)**2.0D0)/scratch1)
          alpha_s = 0.0D0

       else

          alpha_f = dsqrt((a(m)**2.0D0 - cs(m)**2.0D0)/scratch1)
          alpha_s = dsqrt((cf(m)**2.0D0 - a(m)**2.0D0)/scratch1)

       end if


       S = sign(1.0D0, x_magfield_1D(m))
       asq = 1.0D0 / (2.0D0 * a(m)**2.0D0)


       scratch2 = dsqrt(by(m)**2.0D0 + bz(m)**2.0D0)

       if (scratch2 .gt. 0.0D0) then

          beta_y = by(m) / scratch2
          beta_z = bz(m) / scratch2

       else

          beta_y = 1.0D0 / dsqrt(2.0D0)
          beta_z = 1.0D0 / dsqrt(2.0D0)
          
       end if
       

       ! The right eigenvectors

       R(1, 1, m) =  alpha_f * density_1D(m)                          
       R(2, 1, m) = -alpha_f * cf(m)
       R(3, 1, m) = +alpha_s * cs(m) * beta_y * S
       R(4, 1, m) = +alpha_s * cs(m) * beta_z * S
       R(5, 1, m) =  alpha_f * density_1D(m) * a(m)**2.0D0
       R(6, 1, m) =  alpha_s * MHDF(-1) * a(m) * beta_y
       R(7, 1, m) =  alpha_s * MHDF(-1) * a(m) * beta_z
       
       R(1, 2, m) =  0.0D0
       R(2, 2, m) =  0.0D0
       R(3, 2, m) = -beta_z     
       R(4, 2, m) = +beta_y     
       R(5, 2, m) =  0.0D0
       R(6, 2, m) = -MHDF(-1) * beta_z * S  
       R(7, 2, m) =  MHDF(-1) * beta_y * S                               
       
       R(1, 3, m) =  alpha_s * density_1D(m)                                    
       R(2, 3, m) = -alpha_s * cs(m)                
       R(3, 3, m) = -alpha_f * cf(m) * beta_y * S                 
       R(4, 3, m) = -alpha_f * cf(m) * beta_z * S               
       R(5, 3, m) =  alpha_s * density_1D(m) * a(m)**2.0D0                                           
       R(6, 3, m) = -alpha_f * MHDF(-1) * a(m) * beta_y                                               
       R(7, 3, m) =  alpha_f * MHDF(-1) * a(m) * beta_z                
       
       R(1, 4, m) =  1.0D0
       R(2, 4, m) =  0.0D0  
       R(3, 4, m) =  0.0D0  
       R(4, 4, m) =  0.0D0  
       R(5, 4, m) =  0.0D0 
       R(6, 4, m) =  0.0D0  
       R(7, 4, m) =  0.0D0  
       
       R(1, 5, m) =  R(1, 3, m)                                
       R(2, 5, m) = -R(2, 3, m)             
       R(3, 5, m) = -R(3, 3, m)                 
       R(4, 5, m) = -R(4, 3, m)                
       R(5, 5, m) =  R(5, 3, m)                                       
       R(6, 5, m) =  R(6, 3, m)                     
       R(7, 5, m) =  R(7, 3, m)                
       
       R(1, 6, m) =  R(1, 2, m)
       R(2, 6, m) =  R(2, 2, m)
       R(3, 6, m) = -R(3, 2, m)      
       R(4, 6, m) = -R(4, 2, m)      
       R(5, 6, m) =  R(5, 2, m)       
       R(6, 6, m) =  R(6, 2, m)                       
       R(7, 6, m) =  R(7, 2, m)                       
       
       R(1, 7, m) =  R(1, 1, m)             
       R(2, 7, m) = -R(2, 1, m)
       R(3, 7, m) = -R(3, 1, m) 
       R(4, 7, m) = -R(4, 1, m)       
       R(5, 7, m) =  R(5, 1, m)                            
       R(6, 7, m) =  R(6, 1, m)    
       R(7, 7, m) =  R(7, 1, m)     


       ! The left eigenvectors 

       L(1, 1, m) =  asq * 0.0D0  
       L(2, 1, m) =  0.5D0 * 0.0D0 
       L(3, 1, m) =  asq * 0.0D0                    
       L(4, 1, m) =  1.0D0    
       L(5, 1, m) =  L(3, 1, m)         
       L(6, 1, m) =  L(2, 1, m)
       L(7, 1, m) =  L(1, 1, m) 

       L(1, 2, m) = -asq * alpha_f * cf(m)
       L(2, 2, m) =  0.5D0 * 0.0D0     
       L(3, 2, m) = -asq * alpha_s * cs(m)              
       L(4, 2, m) =  0.0D0      
       L(5, 2, m) = -L(3, 2, m)           
       L(6, 2, m) =  L(2, 2, m)
       L(7, 2, m) = -L(1, 2, m)              

       L(1, 3, m) = +asq * alpha_s * cs(m) * beta_y * S
       L(2, 3, m) = -0.5D0 * beta_z         
       L(3, 3, m) = -asq * alpha_f * cf(m) * beta_y * S                       
       L(4, 3, m) =  0.0D0     
       L(5, 3, m) = -L(3, 3, m)                    
       L(6, 3, m) = -L(2, 3, m)               
       L(7, 3, m) = -L(1, 3, m)                     

       L(1, 4, m) = +asq * alpha_s * cs(m) * beta_z * S
       L(2, 4, m) = +0.5D0 * beta_y         
       L(3, 4, m) = -asq * alpha_f * cf(m) * beta_z * S                       
       L(4, 4, m) =  0.0D0   
       L(5, 4, m) = -L(3, 4, m)                   
       L(6, 4, m) = -L(2, 4, m)       
       L(7, 4, m) = -L(1, 4, m)                    

       L(1, 5, m) =  asq * alpha_f / density_1D(m)
       L(2, 5, m) =  0.5D0 * 0.0D0
       L(3, 5, m) =  asq * alpha_s / density_1D(m)                                  
       L(4, 5, m) = -1.0D0 / a(m)**2.0D0     
       L(5, 5, m) =  L(3, 5, m)                   
       L(6, 5, m) =  L(2, 5, m)
       L(7, 5, m) =  L(1, 5, m)                                 

       L(1, 6, m) =  asq * alpha_s * a(m) * beta_y * MHDF(1)       
       L(2, 6, m) = -0.5D0 * beta_z * S * MHDF(1)                      
       L(3, 6, m) = -asq * alpha_f * a(m) * beta_y * MHDF(1)                                       
       L(4, 6, m) =  0.0D0              
       L(5, 6, m) =  L(3, 6, m)                         
       L(6, 6, m) =  L(2, 6, m)                                       
       L(7, 6, m) =  L(1, 6, m)                                   

       L(1, 7, m) =  asq * alpha_s * a(m) * beta_z * MHDF(1)                                              
       L(2, 7, m) =  0.5D0 * beta_y * S * MHDF(1)                   
       L(3, 7, m) = -asq * alpha_s * a(m) * beta_z * MHDF(1)                                            
       L(4, 7, m) =  0.0D0     
       L(5, 7, m) =  L(3, 7, m)                        
       L(6, 7, m) =  L(2, 7, m)                                   
       L(7, 7, m) =  L(1, 7, m)                                  
       
    end do

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Calculate_eigensystem
! ----------------------------------------------------------------------------------------------------------------------------------





! ----------------------------------------------------------------------------------------------------------------------------------
! --- SUBROUTINE: Deallocate the Input States -------------------------------------------------------------------------- [REC07] ---
! ----------------------------------------------------------------------------------------------------------------------------------

  subroutine Deallocate_input_states ()


  ! Deallocating the Riemann problem input states.
  ! ----------------------------------------------

    deallocate (density_L    );   deallocate (density_R    )
    deallocate (x_velocity_L );   deallocate (x_velocity_R )
    deallocate (y_velocity_L );   deallocate (y_velocity_R )
    deallocate (z_velocity_L );   deallocate (z_velocity_R )
    deallocate (x_momentum_L );   deallocate (x_momentum_R )
    deallocate (y_momentum_L );   deallocate (y_momentum_R )
    deallocate (z_momentum_L );   deallocate (z_momentum_R )
    deallocate (energy_L     );   deallocate (energy_R     )
    deallocate (pressure_L   );   deallocate (pressure_R   )
    deallocate (gamma_L      );   deallocate (gamma_R      )
    deallocate (sound_speed_L);   deallocate (sound_speed_R)

    deallocate (x_magfield_L );   deallocate (x_magfield_R )
    deallocate (y_magfield_L );   deallocate (y_magfield_R )
    deallocate (z_magfield_L );   deallocate (z_magfield_R )



  ! Deallocating the working variables.
  ! -----------------------------------

    deallocate (gas_pressure_1D)

    deallocate (a);    deallocate (cax)
    deallocate (bx);   deallocate (by);   deallocate (bz)
    deallocate (cf);   deallocate (ca);   deallocate (cs)

    deallocate (w);   deallocate (lambda)
    deallocate (L);   deallocate (R)

    deallocate (delta_w)

    deallocate (w_hat_L);   deallocate (w_hat_R)
    deallocate (w_new_L);   deallocate (w_new_R)
    deallocate (w_L);       deallocate (w_R)


    ! Specific for TVD

    if (RECONSTRUCT_TYPE == 'T') then

       deallocate (wp);   deallocate (w0);   deallocate (wm);

       deallocate (rp_p);   deallocate (rm_p);   deallocate (rd_p);
       deallocate (rp_m);   deallocate (rm_m);   deallocate (rd_m);

       deallocate (delta_V)
       deallocate (sigma)
       deallocate (dw)

    end if


    ! Specific for PPM(L)

    if ((RECONSTRUCT_TYPE == 'P') .or. (RECONSTRUCT_TYPE == 'L')) then

       deallocate (wh)

       deallocate (dwL);   deallocate (dwR);   deallocate (dwC);   deallocate (dwG)
       deallocate (daL);   deallocate (daR);   deallocate (daC);   deallocate (daG)
       
       deallocate (PPM_dw)   
       deallocate (PPM_w6)
       deallocate (delta_a)  

    end if

    return



! ----------------------------------------------------------------------------------------------------------------------------------
  end subroutine Deallocate_input_states
! ----------------------------------------------------------------------------------------------------------------------------------



end module recon
