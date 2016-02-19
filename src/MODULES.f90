!    Copyright (C) 2014 The SPEED FOUNDATION
!    Author: Ilario Mazzieri
!
!    This file is part of SPEED.
!
!    SPEED is free software; you can redistribute it and/or modify it
!    under the terms of the GNU Affero General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    SPEED is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with SPEED.  If not, see <http://www.gnu.org/licenses/>.



!> @brief Set maximal bounds.
!! @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
module max_var
        integer, parameter :: max_dim           = 3      !<dimension of the problem
        integer, parameter :: max_el_conf       = 200    !<max number of neighbouring elements 
        integer, parameter :: nofqp             = 6      !<max number of 1-D quadrature point per element 
        integer, parameter :: nofinr            = 500    !<max number of newton rapson iterations
        integer, parameter :: max_quad_points   = 8000   !<max number of quadrature nodes on a DG surface
        real,    parameter :: tol_nl            = 0.001  !<tolerance on nonlinear calculations
end module max_var


!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains mesh structure (scratch)
module str_mesh_scratch

  use max_var , only: nofqp !<max number of 1-D quadrature point per element

  type scratch_ELEMENT !<contains coordinates of quadrature nodes
     real*8, dimension(nofqp**2) :: x_nq !<x coordinate
     real*8, dimension(nofqp**2) :: y_nq !<y coordinate
     real*8, dimension(nofqp**2) :: z_nq !<z coordinate
  end type scratch_ELEMENT

end module str_mesh_scratch


!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains mesh structure for DG interface elements 
module str_mesh

 use max_var

 type ELEMENT  !< Interface DG Element (quad)
   
   integer*4:: mat          !< block id (material)
   integer*4:: ind_el             !< element index (global numeration)
   integer*4:: face_el      !< face index (from 1 to 6)
   integer*4:: spct_deg     !< polynomila degree inside the element 
   integer*4:: quad_rule    !< number of quadrature point in 1D
   integer*4:: nofne        !< number of neighbouring elements
   integer*4:: proj_yn      !< 1 if the element project quad nodes 0 otherwise 
   real*8   :: nx,ny    !< normal to the element
   

   real*8, dimension(nofqp**2) :: wx_pl  !< weights of the quadrature rule (x)
   real*8, dimension(nofqp**2) :: wy_pl  !< weights of the quadrature rule (y)


 end type ELEMENT

END MODULE str_mesh


!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains mesh structure for DG elements after pre-processing 
module str_mesh_after

 use max_var,  only : max_quad_points, max_el_conf

 type ELEMENT_after
   
   integer*4:: mat          !< block id (material)
   integer*4:: ind_el             !< element index (global numeration)
   integer*4:: face_el      !< face index (from 1 to 6)
   integer*4:: spct_deg     !< polynomial degree 
   integer*4:: quad_rule    !< number of quadrature point in 1D
   integer*4:: nofne        !< number of neighbouring elements
   real*8   :: nx,ny    !< normal to the element

   real*8, dimension(max_quad_points) :: x_pl  !< quadrature points x (+,+)
   real*8, dimension(max_quad_points) :: y_pl  !< quadrature points y (+,+)

   real*8, dimension(max_quad_points) :: x_mn  !< quadrature points x (+,-)
   real*8, dimension(max_quad_points) :: y_mn  !< quadrature points y (+,-)


   real*8, dimension(max_quad_points) :: wx_pl !< quadrature weights x
   real*8, dimension(max_quad_points) :: wy_pl !< quadrature weights y

 
   integer*4, dimension(max_quad_points,0:3) :: omega_minus !< matrix containing neigh el. info (quad points)
   integer*4, dimension(max_el_conf,0:2) :: conf            !< matrix containing neigh el. info (mat,el,ind,face)


 end type ELEMENT_after


END MODULE str_mesh_after



!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains structure for jump matrices
module DGJUMP

use max_var , only: max_el_conf

type matrix 

   real*8, dimension(:,:), pointer :: MJUMP !< jump matrix   
   real*8, dimension(:,:), pointer :: MJUMP_only_uv  !<jump matrix fot testmode (only [u][v])

end type matrix


type el4loop   !< element structure for time loop (RAM saving)

   integer*4:: ind                                   !< element index (gl. num.)
   integer*4:: face                                  !< element face  (1<-->6)
   integer*4:: deg                                   !< pol.degree
   integer*4:: num_of_ne                             !< number of neigh. elements
   integer*4:: nnz_plus                              !< nonzero els. jump matrix (+,+)
   integer*4:: nnz_minus                             !< nonzero els. jump matrix (+,-) 
   integer*4:: nnz_col                               !< length of u_m vector in timeloop
   integer*4:: nnz_plus_only_uv                      !< nonzero els. jump matrix (+,+) [u][v] - testmode
   integer*4:: nnz_minus_only_uv                     !< nonzero els. jump matrix (+,-) [u][v] - testmode
   integer*4:: nnz_col_only_uv                       !< length of u_m vector in timeloop [u][v] - testmode
   
   integer*4, dimension(:), pointer :: IPlus, JPlus  !< RCS format rows
   integer*4, dimension(:), pointer :: IMin, JMin    !< RCS format columns
   integer*4, dimension(max_el_conf,0:2) :: el_conf  !< matrix for neigh. elements
   real*8, dimension(:), pointer :: matPlus, matMin  !< RCS format matrix
   real*8, dimension(:,:), pointer :: matP           !< jump matrix (+,+)
   type(matrix), dimension(:), pointer :: matM       !< jump matrix (+,-)

   integer*4, dimension(:), pointer :: IPlus_only_uv, JPlus_only_uv  !< RCS format rows [u][v] - testmode
   integer*4, dimension(:), pointer :: IMin_only_uv, JMin_only_uv    !< RCS format columns [u][v] - testmode
   real*8, dimension(:), pointer :: matPlus_only_uv, matMin_only_uv  !< RCS format matrix [u][v] - testmode
   real*8, dimension(:,:), pointer :: matP_only_uv                   !< jump matrix (+,+) [u][v] - testmode
   type(matrix), dimension(:), pointer :: matM_only_uv               !< jump matrix (+,-) [u][v] - testmode

end type 

end module





!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief Contains SPEED paramters (used in  MAKE_DG_INTERFACE_CONDITIONS) 
module speed_par_dg

      use max_var
      use str_mesh 
      use str_mesh_scratch                   
      use DGJUMP
!      use speed_par, only: nelem_dg

      type(el4loop), dimension(:), allocatable :: el_new                        
      type(ELEMENT), dimension(:), allocatable :: dg_els
      type(scratch_ELEMENT), dimension(:), allocatable :: scratch_dg_els

end module speed_par_dg

!> @author Ilario Mazzieri
!> @date September, 2013 - Creation
!> @version 1.0
!> @brief  Quick-sort algorithm
module qsort_c_module    

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A)
  integer, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  integer, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  integer :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition
end module qsort_c_module

!> @author Filippo Gatti
!> @date February, 2016
!> @version 1.0
!> @brief  Module for non linear calculation in 2D

module nonlinear2d
    use max_var, only : tol_nl
    contains
        ! COMPUTE ELASTIC STIFFNESS MATRIX
        subroutine stiff_matrix(lambda,mu,DEL_ijhk)
            
            implicit none
            
            real*8, intent(in)                      :: lambda,mu
            real*8, dimension(4,4), intent(inout)   :: DEL_ijhk

            DEL_ijhk(:,:) = 0d0
            DEL_ijhk(1:3,1:3)   = DEL_ijhk(1:3,1:3) + lambda
            DEL_ijhk(1,1)       = DEL_ijhk(1,1) + 2*mu
            DEL_ijhk(2,2)       = DEL_ijhk(2,2) + 2*mu
            DEL_ijhk(3,3)       = DEL_ijhk(3,3) + 2*mu
            DEL_ijhk(4,4)       = DEL_ijhk(4,4) + mu
            
            return
        
        end subroutine stiff_matrix
        
        ! MISES YIELDING LOCUS AND GRADIENT
        subroutine mises_yld_locus(Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)
            
            implicit none
            
            real*8,               intent(in)    :: R,sigma_yld
            real*8, dimension(4), intent(in)    :: Sigma_ij,X_ij
            real*8, dimension(4), intent(out)   :: gradF_mises
            real*8,               intent(out)   :: F_mises
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,2.0/)
            real*8, dimension(4)                :: Sigma_ij_dev
            real*8                              :: tau_eq_mises
            integer                             :: k
            
            call tensor_components (Sigma_ij, Sigma_ij_dev)
            call tau_mises(Sigma_ij_dev-X_ij, tau_eq_mises)
            F_mises             =   tau_eq_mises-sigma_yld-R
            gradF_mises         =   0d0
            do k=1,4
                gradF_mises(k) = 0.5*A(k)*(Sigma_ij_dev(k)-X_ij(k))/tau_eq_mises
            end do
        
            return

        end subroutine mises_yld_locus
        
        ! OCTAHEDRAL SHEAR STRESS
        subroutine tau_mises(tensor,tau_eq_mises)
            
            implicit none
            
            real*8, dimension(4), intent(in)    :: tensor
            real*8,               intent(out)   :: tau_eq_mises
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,2.0/) 
            integer                             :: k
            
            tau_eq_mises = 0.0d0
            do k=1,4
                tau_eq_mises = tau_eq_mises+A(k)*(tensor(k)**2)
            end do
            tau_eq_mises = sqrt(0.5*tau_eq_mises)

            return

        end subroutine
        
        ! TENSOR COMPONENTS (SPHERICAL & DEVIATORIC)
        subroutine tensor_components(Sigma_ij, Sigma_ij_dev)
            
            implicit none
            
            real*8, dimension(3), intent(in)    :: Sigma_ij
            real*8, dimension(3), intent(out)   :: Sigma_ij_dev
            real*8                              :: Sigma_P
            integer                             :: k

            Sigma_ij_dev = Sigma_ij
            Sigma_P=sum(Sigma_ij(1:3))/3
            do k=1,3
                Sigma_ij_dev(k) = Sigma_ij_dev(k)-Sigma_P
            end do

            return

        end subroutine tensor_components
            
        ! CHECK PLASTIC CONSISTENCY (KKT CONDITIONS)
        subroutine check_plasticity (dSigma_ij_trial, Sigma_ij_start, X_ij, R, sigma_yld, &
            st_elp, alpha_elp)
            
            implicit none
            
            integer*4                               :: k
            real*8                                  :: Fstart, Ftrial,checkload
            real*8, dimension(4)                    :: gradFstart, gradFtrial, Sigma_ij_trial
            real*8,                 intent(in)      :: R,sigma_yld
            real*8, dimension(4),   intent(in)      :: X_ij,Sigma_ij_start
            real*8,                 intent(out)     :: alpha_elp
            logical,                intent(out)     :: st_elp
            real*8, dimension(4),   intent(inout)   :: dSigma_ij_trial
            real*8, dimension(4),   parameter       :: A=(/1.0,1.0,1.0,2.0/)
            ! Yield function at Sigma_ij_start
            call mises_yld_locus(Sigma_ij_start, X_ij, R, sigma_yld, &
                Fstart, gradFstart)
            ! Yield function at Sigma_ij_trial
            Sigma_ij_trial=Sigma_ij_start+dSigma_ij_trial
            call mises_yld_locus(Sigma_ij_trial, X_ij, R, sigma_yld, &
                Ftrial, gradFtrial)

            write(*,*) "*********************************"
            write(*,*) "Fstart:",Fstart,"Ftrial:",Ftrial
            ! LOADING CONDITION    
            checkload=0d0
            do k=1,4
                checkload=checkload+10*gradFstart(k)*dSigma_ij_trial(k)
            end do
            ! KKT CONDITION
            if (abs(Fstart) .le. tol_nl) then
                if (checkload .ge. 0) then
                    alpha_elp = 0d0
                    st_elp    = .true.
                    write(*,*) "PURE-PLASTIC STEP"
                else
                    if (Ftrial .gt. tol_nl) then
                        alpha_elp = 2*sigma_yld/(2*sigma_yld+Ftrial)
                        st_elp    = .true.
                        write(*,*) "ELASTO-PLASTIC WITH REVERSAL STEP"
                    else
                        alpha_elp = 1d0
                        st_elp    = .false. 
                        write(*,*) "PURE ELASTIC WITH REVERSAL STEP"
                    endif
                end if
            elseif (Fstart .lt. -tol_nl) then
                if (Ftrial .lt. tol_nl) then
                    alpha_elp = 1d0
                    st_elp    = .false.
                    write(*,*) "PURE ELASTIC STEP"
                else
                    call gotoFtan(Sigma_ij_start,dSigma_ij_trial,Fstart,Ftrial,X_ij,R,sigma_yld,alpha_elp)
                    st_elp    = .true.
                    write(*,*) "ELASTO-PLASTIC STEP"
                end if
            elseif (Fstart .gt. tol_nl) then
                write(*,*) "ERROR!"
                write(*,*) "Fstart: ",Fstart,">",tol_nl,"!!!!"
                write(*,*) "ERROR!"
                stop
            end if
            dSigma_ij_trial=Sigma_ij_start+dSigma_ij_trial*alpha_elp
            call mises_yld_locus(dSigma_ij_trial,X_ij,R,sigma_yld,Fstart,gradFstart)
            write(*,*) "check plasticity:",Fstart,"<=",tol_nl
            write(*,*) "ALPHA:",alpha_elp
            write(*,*) "*********************************"
            write(*,*) ""
        end subroutine check_plasticity
        
        ! NR ALGORITHM AND DRIFT CORRECTION
        subroutine plastic_corrector (dEpsilon_ij_alpha, Sigma_ij, sigma_yld, X_ij, &
            R, b_lmc, Rinf_lmc, C_lmc, kapa_lmc, mu, lambda, dEpsilon_ij_pl)
            
            implicit none
            
            real*8,               intent(in)    :: sigma_yld            ! first yield limit
            real*8,               intent(in)    :: b_lmc, Rinf_lmc      ! Lamaitre and Chaboche parameters (isotropic hardening)
            real*8,               intent(in)    :: C_lmc, kapa_lmc      ! Lamaitre and Chaboche parameters (kinematic hardening)
            real*8,               intent(in)    :: mu, lambda           ! elastic parameters
            real*8,               intent(inout) :: R                    ! starting mises radius
            real*8, dimension(4), intent(inout) :: Sigma_ij             ! starting stress state
            real*8, dimension(4), intent(inout) :: X_ij                 ! starting back stress
            real*8, dimension(4), intent(inout) :: dEpsilon_ij_alpha    ! percentage of elastic-plastic strain
            real*8, dimension(4), intent(inout) :: dEpsilon_ij_pl       ! percentage of plastic strain

            real*8, dimension(4)                :: dX_ij, gradF_0, gradF_mises
            real*8, dimension(4)                :: temp_vec
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)              :: DEL_ijhk
            integer,            parameter       :: N_incr = 10
            integer                             :: i,j,k
            real*8                              :: dR, dPlastMult, F_mises_0, F_mises,Ffinal

            ! COMPUTE ELASTIC STIFFNESS MATRIX
            call stiff_matrix(lambda,mu,DEL_ijhk)

            ! ELASTO-PLASTIC SUB-STEPPING
            dEpsilon_ij_alpha = dEpsilon_ij_alpha/N_incr

            do i = 1,N_incr

                ! PREDICTION
                call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises_0, gradF_0)
                write(*,*) "prediction",F_mises_0
                write(*,*) "grad",gradF_0
                ! COMPUTE PLASTIC MULTIPLIER
                call compute_plastic_modulus(dEpsilon_ij_alpha, Sigma_ij, X_ij, R, mu, lambda, sigma_yld, &
                    b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dPlastMult)
                
                ! HARDENING INCREMENTS
                call hardening_increments(Sigma_ij, R, X_ij, sigma_yld, &
                    b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR, dX_ij)
                dR = dPlastMult*dR
                dX_ij = dPlastMult*dX_ij

                ! VARIABLE UPDATE
                do k=1,4 ! plastic strains
                    dEpsilon_ij_pl(k)=dEpsilon_ij_pl(k)+dPlastMult*gradF_0(k)*A(k)
                end do
                R=R+dR          ! isotropic hardening update
                X_ij=X_ij+dX_ij ! back-stress update
                ! stress update
                do j = 1,4
                    do k = 1,4
                        Sigma_ij(j)=Sigma_ij(j)+DEL_ijhk(k,j)*(dEpsilon_ij_alpha(k)-A(k)*dPlastMult*gradF_0(k))
                    end do
                end do
                ! DRIFT CORRECTION
                call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, Ffinal, gradF_0)
                write(*,*) "Ffinal (BD): ",Ffinal
                if (Ffinal .gt. tol_nl) then
                    call drift_corr(.true.,Sigma_ij, X_ij, R, sigma_yld,&
                        b_lmc, Rinf_lmc, C_lmc, kapa_lmc, lambda, mu, dEpsilon_ij_pl)
                    call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, Ffinal, gradF_0)
                end if
                write(*,*) "Ffinal (AD): ",Ffinal
                write(*,*) ""
            end do
        end subroutine plastic_corrector
        
        ! HARDENING MODULUS AND PLASTIC MULTIPLIER INCREMENT
        subroutine compute_plastic_modulus(dEpsilon_ij, Sigma_ij, X_ij, R, mu, lambda, &
            sigma_yld, b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dPlastMult)
            implicit none
            real*8, dimension(4), intent(in)   :: dEpsilon_ij         ! percentage of elastic-plastic strain
            real*8, dimension(4), intent(in)   :: Sigma_ij            ! starting stress state
            real*8, dimension(4), intent(in)   :: X_ij                ! starting back stress
            real*8,                 intent(in) :: R                   ! starting mises radius
            real*8,                 intent(in) :: b_lmc, Rinf_lmc     ! Lamaitre and Chaboche parameters (isotropic hardening)
            real*8,                 intent(in) :: C_lmc, kapa_lmc     ! Lamaitre and Chaboche parameters (kinematic hardening)
            real*8,                 intent(in) :: mu, lambda          ! elastic parameters
            real*8,                 intent(in) :: sigma_yld           ! first yielding limit
            real*8,                 intent(out):: dPlastMult          ! plastic multiplier increment
            real*8, dimension(4)               :: gradF_mises
            real*8, dimension(4),   parameter  :: A =(/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)             :: DEL_ijhk
            real*8                            :: h_iso, h_kin        ! isotropic and kinematic hardening modula
            real*8                            :: F_mises
            real*8                            :: h_lmc               ! total hardening modulus
            real*8                            :: temp_vec
            integer                         :: j,k
            
            ! COMPUTE ELASTIC STIFFNESS MATRIX
            call stiff_matrix(lambda,mu,DEL_ijhk)

            ! COMPUTE HARDENING MODULUS
            call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)

            h_iso = b_lmc*(Rinf_lmc-R)
            h_kin = 0d0
            do k=1,4
                h_kin=h_kin+kapa_lmc*X_ij(k)*gradF_mises(k)
            end do
            h_kin=C_lmc-h_kin
            h_lmc = h_iso+h_kin
            ! COMPUTE PLASTIC MULTIPLIER
            temp_vec   = 0d0
            dPlastMult = 0d0

            do k=1,4
                do j=1,4
                    temp_vec   = temp_vec+gradF_mises(j)*DEL_ijhk(j,k)*gradF_mises(k)*A(k)
                    dPlastMult = dPlastMult+gradF_mises(j)*DEL_ijhk(j,k)*dEpsilon_ij(k)
                end do
            end do

            write(*,*) "nominatore",dPlastMult
            write(*,*) "deps",dEpsilon_ij
            write(*,*) "grad",gradF_mises
            dPlastMult = dPlastMult/(h_lmc+temp_vec)
            if (dPlastMult .lt. 0) then
                write(*,*) "DPLAST NEGATIVE"
                stop
            end if

            return

        end subroutine compute_plastic_modulus

        ! INCREMENTS OF INTRINSIC STATIC VARIABLES
        subroutine hardening_increments(Sigma_ij, R, X_ij, sigma_yld, &
            b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR, dX_ij)
            
            implicit none
            
            real*8, dimension(4), intent(in)    :: Sigma_ij        ! actual stress state
            real*8, dimension(4), intent(in)    :: X_ij            ! actual back stress state
            real*8,               intent(in)    :: sigma_yld       ! first yielding limit
            real*8,               intent(in)    :: R               ! actual mises radius
            real*8,               intent(in)    :: b_lmc, Rinf_lmc ! Lamaitre and Chaboche parameters (isotropic hardening)
            real*8,               intent(in)    :: C_lmc, kapa_lmc ! Lamaitre and Chaboche parameters (kinematic hardening)
            real*8, dimension(4), intent(out)   :: dX_ij           ! back stress increment
            real*8,               intent(out)   :: dR              ! mises radius increment
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,0.5/)
            real*8, dimension(4)                :: gradF_mises
            real*8                              :: F_mises
            integer                             :: k

            ! INCREMENT IN ISOTROPIC HARDENING VARIABLES (R)
            dR = b_lmc*(Rinf_lmc-R)
            ! INCREMENT IN KINEMATIC HARDENING VARIABLES (Xij)
            call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)
            dX_ij=2*A*gradF_mises*C_lmc/3-X_ij*kapa_lmc

            return

        end subroutine hardening_increments
        
        ! DRIFT CORRECTION (RADIAL RETURN)
        subroutine drift_corr(drift,Sigma_ij, X_ij, R, sigma_yld, &
            b_lmc, Rinf_lmc, C_lmc, kapa_lmc, lambda, mu, dEpsilon_ij_pl)
            
            implicit none

            logical, intent(in)                     :: drift
            real*8, dimension(4), intent(inout)     :: Sigma_ij, X_ij, dEpsilon_ij_pl
            real*8,                 intent(inout)   :: R
            real*8,                 intent(in)      :: lambda, mu, sigma_yld, b_lmc, Rinf_lmc, C_lmc, kapa_lmc
            real*8, dimension(4)                    :: gradF_mises,gradF0,Sigma_temp,Sigma_dev_temp,Sigma_dev_ij
            real*8, dimension(4), parameter         :: A = (/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)                  :: DEL_ijhk
            real*8:: Fmises,Fmises0,beta,h_kin,h_iso,h_lmc,dbeta,err0,err1
            integer:: k,j
            
            ! INITIAL PLASTIC CONDITION
            call mises_yld_locus(Sigma_ij, X_ij, R, sigma_yld, Fmises0, gradF0)
            beta=0d0
            do 
                call mises_yld_locus(Sigma_ij,X_ij,R,sigma_yld,Fmises,gradF_mises)
                ! COMPUTE BETA FOR DRIFT CORRECTION
                dbeta=0d0
                do k=1,4
                    dbeta = dbeta+gradF0(k)*gradF_mises(k)
                end do
                beta=beta-Fmises/dbeta
                ! STRESS CORRECTION (RADIAL RETURN)
                Sigma_temp=Sigma_ij+beta*gradF0*A
                call mises_yld_locus(Sigma_temp,X_ij,R,sigma_yld,Fmises,gradF_mises)
                call tensor_components(Sigma_temp, Sigma_dev_temp)
                call tensor_components(Sigma_ij,Sigma_dev_ij)
                call tau_mises(Sigma_dev_ij-Sigma_dev_temp,err0)
                call tau_mises(Sigma_dev_ij-X_ij,err1)
                err0=err0/err1
                Sigma_ij=Sigma_temp
                !if (abs(Fmises) .le. tol_nl) then
                if (err0 .le. 0.0000001 .or. abs(Fmises) .le. tol_nl) then
                    exit
                end if
            end do

            return

        end subroutine drift_corr

        subroutine gotoFtan(start0,dtrial0,F0,Ftrial,center,radius,s0,alpha)
            
            implicit none

            real*8, dimension(3), intent(in)        :: start0,dtrial0,center
            real*8,               intent(in)        :: radius,s0
            real*8,               intent(inout)     :: F0,Ftrial
            real*8,               intent(out)       :: alpha
            real*8, dimension(3)                    :: start,gradF,dev,dev_temp,temp,dev0
            real*8                                  :: Fstart,err0,err1,beta
            integer                                 :: counter
            
            call tensor_components(start0,dev0)
            alpha=F0/(F0-Ftrial)
            start=start0
            temp=start+alpha*dtrial0
            do counter=1,10
                call tensor_components(temp, dev_temp)
                call tensor_components(start,dev)
                call tau_mises(-dev+dev_temp,err0)
                call tau_mises(dev-dev0,err1)
                write(*,*) "err0",err0
                write(*,*) "err1",err1
                err0=err0/err1
                start=temp
                call mises_yld_locus(start,center,radius,s0,Fstart,gradF)
                if (abs(Fstart).le.tol_nl .or. err0.lt.0.000001) then
                    exit
                end if
                beta=sum(gradF*dtrial0)
                alpha=alpha-Fstart/beta
                temp=start-(Fstart/beta)*dtrial0
            end do
            write(*,*) "F:",Fstart
            write(*,*) "err0",err0
            write(*,*) "gotoF: stress:",start
            write(*,*) ""

            return

        end subroutine gotoFtan
         
        subroutine gotoFsec(start0,dtrial0,F0,Ftrial,center,radius,s0,alpha)
            
            implicit none

            real*8, dimension(4), intent(in)        :: start0,dtrial0,center
            real*8,               intent(in)        :: radius,s0
            real*8,               intent(inout)     :: F0,Ftrial
            real*8,               intent(out)       :: alpha
            real*8, dimension(4)                    :: start,gradF,dev,dev_temp,temp,dev0
            real*8                                  :: alphanew,err0,err1,beta
            integer                                 :: counter
            
            call tensor_components(start0,dev0)
            alpha=F0/(F0-Ftrial)
            beta=0d0
            alphanew=0d0
            start=start0
            temp=start+alpha*dtrial0
            do counter=1,10
                call tensor_components(temp, dev_temp)
                call tensor_components(start,dev)
                call tau_mises(-dev+dev_temp,err0)
                call tau_mises(dev-dev0,err1)
                err0=err0/err1
                call mises_yld_locus(start,center,radius,s0,F0,gradF)
                call mises_yld_locus(temp,center,radius,s0,Ftrial,gradF)
                start=temp
                if (abs(Ftrial).le.tol_nl .or. err0.lt.0.0000001) exit
                alphanew=alpha-(Ftrial/(Ftrial-F0))*(alpha-beta)
                beta=alpha
                alpha=alphanew
                temp=start-(alpha-beta)*dtrial0
                
            enddo

            return

        end subroutine gotoFsec
end module nonlinear2d
