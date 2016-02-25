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
       real*8, parameter :: FTOL = 1.0000000000D0
       real*8, parameter :: LTOL = 0.010000000000D0
       real*8, parameter :: STOL = 0.0010000000000000D0

    contains
        ! COMPUTE ELASTIC STIFFNESS MATRIX
        subroutine stiff_matrix(lambda,mu,DEL)
            
            implicit none
            
            real*8, intent(in)                      :: lambda,mu
            real*8, dimension(4,4), intent(inout)   :: DEL

            DEL(:,:) = 0d0
            DEL(1:3,1:3)   = DEL(1:3,1:3) + lambda
            DEL(1,1)       = DEL(1,1) + 2*mu
            DEL(2,2)       = DEL(2,2) + 2*mu
            DEL(3,3)       = DEL(3,3) + 2*mu
            DEL(4,4)       = DEL(4,4) + mu
            
            return
        
        end subroutine stiff_matrix
        
        ! MISES YIELDING LOCUS AND GRADIENT
        subroutine mises_yld_locus(stress, center, radius, syld, FM, gradF)
            
            implicit none
            
            real*8,               intent(in)    :: radius,syld
            real*8, dimension(4), intent(in)    :: stress,center
            real*8, dimension(4), intent(out)   :: gradF
            real*8,               intent(out)   :: FM
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,2.0/)
            real*8, dimension(4)                :: dev
            real*8                              :: tau_eq
            integer*4                           :: k
            
            call tensor_components (stress, dev)
            call tau_mises(dev-center, tau_eq)
            FM             =   tau_eq-syld-radius
            gradF = 0.5d0*A*(dev-center)/tau_eq
        
            return

        end subroutine mises_yld_locus
        
        ! OCTAHEDRAL SHEAR STRESS
        subroutine tau_mises(stress,tau_eq)
            
            implicit none
            
            real*8, dimension(4), intent(in)    :: stress
            real*8,               intent(out)   :: tau_eq
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,2.0/) 
            integer*4                           :: k
            
            tau_eq = 0.0d0
            do k=1,4
                tau_eq = tau_eq+A(k)*(stress(k)**2)
            end do
            tau_eq = sqrt(0.5*tau_eq)

            return

        end subroutine
        
        ! TENSOR COMPONENTS (SPHERICAL & DEVIATORIC)
        subroutine tensor_components(stress, dev)
            
            implicit none
            
            real*8, dimension(4), intent(in)    :: stress
            real*8, dimension(4), intent(out)   :: dev
            real*8                              :: press 
            integer*4                           :: k

            dev = stress
            press = sum(stress(1:3))/3
            dev(1:3) = dev(1:3)-press
            return

        end subroutine tensor_components
            
        ! CHECK PLASTIC CONSISTENCY (KKT CONDITIONS)
        subroutine check_plasticity (dtrial, stress0, center, radius, syld, &
            st_elp, alpha_elp)
            
            implicit none
            
            integer*4                               :: k
            real*8                                  :: FS, FT,checkload
            real*8, dimension(4)                    :: gradFS, gradFT, stress1
            real*8,                 intent(in)      :: radius,syld
            real*8, dimension(4),   intent(in)      :: center,stress0
            real*8,                 intent(out)     :: alpha_elp
            logical,                intent(out)     :: st_elp
            real*8, dimension(4),   intent(inout)   :: dtrial
            
            stress1=stress0+dtrial
            call mises_yld_locus(stress0,center,radius,syld,FS,gradFS)
            call mises_yld_locus(stress1,center,radius,syld,FT,gradFT)
            checkload=sum(gradFS*dtrial)/sum(gradFS**2)/sum(dtrial**2)

            if (abs(FS).le.FTOL) then
                if (checkload.ge.-LTOL) then
                    alpha_elp = 0d0
                    st_elp    = .true.
                else
                    if (FT.lt.-FTOL) then
                        alpha_elp = 1d0
                        st_elp    = .false.
                    elseif(FT.gt.FTOL) then
                        call gotoFpegasus(stress0,dtrial,center,radius,syld,10,alpha_elp)
                        st_elp    = .true.
                    endif
                end if
            elseif (FS.lt.-FTOL) then
                if (FT.le.FTOL) then
                    alpha_elp = 1d0
                    st_elp    = .false.
                else
                    call gotoFpegasus(stress0,dtrial,center,radius,syld,1,alpha_elp)
                    st_elp    = .true.
                end if
            elseif (FS.gt.FTOL) then
                write(*,*) 'ERROR: FS: ',FS,'>',FTOL,'!!!!'
                alpha_elp  = 0d0
                st_elp = .true.
            end if
            dtrial=stress0+dtrial*alpha_elp
            call mises_yld_locus(dtrial,center,radius,syld,FS,gradFS)
        end subroutine check_plasticity
        
        subroutine plastic_corrector(dEps_alpha,stress,center,syld, &
            radius,biso,Rinf,Ckin,kkin,mu,lambda,dEpl)

            implicit none
            real*8, dimension(4), intent(inout) :: dEps_alpha,stress,center,dEpl 
            real*8,                 intent(inout) :: radius 
            real*8,                 intent(in)    :: syld,biso,Rinf,Ckin,kkin,mu,lambda
            real*8, dimension(4,4)            :: DEL
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,0.5/)
            real*8                              :: Ttot,deltaTk,qq,R1,R2,dR1,dR2,err0,err1
            real*8                              :: FM,hard1,hard2,deltaTmin
            real(8)                             :: Resk
            logical                             :: flag_fail
            real*8, dimension(4)                :: gradFM,S1,S2,X1,X2
            real*8, dimension(4)                :: dS1,dS2,dX1,dX2,dEpl1,dEpl2

            call stiff_matrix(lambda,mu,DEL)
            deltaTk = 1.0d0
            Ttot    = 0.0d0
            deltaTmin = 0.01d0
            do while (Ttot.lt.1d0-FTOL)
                Resk     = 0d0
                dS1(1:4) = 0d0
                dX1(1:4) = 0d0
                dS2(1:4) = 0d0
                dX2(1:4) = 0d0
                dR1      = 0d0
                dR2      = 0d0
                dEpl1(1:4) = 0d0
                dEpl2(1:4) = 0d0

                ! FIRST ORDER COMPUTATION
                call ep_integration(dEps_alpha*deltaTk,stress,center,radius,syld,&
                    mu,lambda,biso,Rinf,Ckin,kkin,dS1,dX1,dR1,dEpl1,hard1)

                S1 = stress + dS1
                X1 = center + dX1 
                R1 = radius + dR1
               
                ! SECOND ORDER COMPUTATION
                call ep_integration(dEps_alpha*deltaTk,S1,X1,R1,syld,mu,lambda,&
                    biso,Rinf,Ckin,kkin,dS2,dX2,dR2,dEpl2,hard2)
                
                ! TEMPORARY VARIABLES
                S1 = stress + 0.5d0*(dS1+dS2)
                X1 = center + 0.5d0*(dX1+dX2)
                R1 = radius + 0.5d0*(dR1+dR2)
                dEpl1 = 0.5d0*(dEpl1+dEpl2)

                ! ERROR
                call tau_mises(dS2-dS1,err0)
                call tau_mises(S1,err1)

                Resk=0.5d0*max(epsilon(Resk),err0/err1)
                
                if (Resk.le.STOL) then
                    
                    stress = S1
                    center = X1
                    radius = R1

                    call mises_yld_locus (stress, center,radius,syld,FM,gradFM)
                    if (FM.gt.FTOL) then
                        call drift_corr(stress,center,radius,syld,&
                                biso,Rinf,Ckin,kkin,lambda,mu,dEpl)
                    endif
                   
                    Ttot=Ttot+deltaTk
                    qq = min(0.9d0*sqrt(STOL/Resk),1.1d0) 
                    if (flag_fail) then
                        qq = min(qq,1.0d0)
                    endif
                    flag_fail=.false.
                    deltaTk=qq*deltaTk
                    deltaTk=max(qq*deltaTk,deltaTmin)
                    deltaTk=min(deltaTk,1d0-Ttot)

                else
                    qq=max(0.9d0*sqrt(STOL/Resk),deltaTmin)
                    deltaTk=qq*deltaTk
                    flag_fail=.true.
                end if
            end do
        end subroutine plastic_corrector

        subroutine ep_integration(dStrain,Stress,center,radius,syld,mu,lambda,biso,Rinf,&
            Ckin,kapakin,dStress,dcenter,dradius,dEplast,hard)
            
            implicit none
            real*8                , intent(in) :: radius,syld,mu,lambda,biso,Rinf,Ckin,kapakin
            real*8, dimension(4), intent(in) :: dStrain,Stress,center
            real*8, dimension(4), intent(inout):: dStress,dcenter,dEplast
            real*8,                 intent(inout):: dradius
            real*8,                 intent(out)  :: hard
            real*8, dimension(4)             :: gradF
            real*8, dimension(4), parameter  :: A = (/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)         :: DEL
            real*8                             :: Fmises,dPlast
            integer*4                        :: j,k
            
            ! PREDICTION
            call mises_yld_locus (Stress,center,radius,syld,Fmises,gradF)
            call stiff_matrix(lambda,mu,DEL)

            ! PLASTIC MULTIPLIER
            call compute_plastic_modulus(dStrain,Stress,center,radius,mu,lambda,syld, &
                biso,Rinf,Ckin,kapakin,dPlast,hard)
            
            ! INCREMENTS
            call hardening_increments(Stress,radius,center,syld, &
                biso,Rinf,Ckin,kapakin,dradius,dcenter)
            
            dradius = dPlast*dradius
            dcenter = dPlast*dcenter
            dEplast = 0d0
            do k=1,4 
                dEplast(k)=dEplast(k)+dPlast*gradF(k)*A(k)
            end do
            do j = 1,4 ! stress increment
                do k = 1,4
                    dstress(j)=DEL(k,j)*(dstrain(j)-A(k)*dPlast*gradF(k))  
                end do
            end do
            
            return

        end subroutine ep_integration


        subroutine compute_plastic_modulus(dEps,stress,center,radius,mu,lambda, &
            syld,biso,Rinf,Ckin,kkin,dPlastMult,hard)

            implicit none
            real*8,                 intent(in) :: mu,lambda,syld   
            real*8,                 intent(in) :: radius,biso,Rinf,Ckin,kkin
            real*8, dimension(4), intent(in) :: dEps,stress,center
            real*8,                 intent(out):: dPlastMult,hard
            real*8                             :: temp_vec,FM
            real*8, dimension(4)             :: gradF
            real*8, dimension(4), parameter  :: A =(/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)         :: DEL
            integer*4                          :: j,k
            
            call stiff_matrix(lambda,mu,DEL)
            call mises_yld_locus(stress,center,radius,syld,FM,gradF)
            
            hard = Ckin+biso*(Rinf-radius)
            hard = hard-kkin*sum(center*gradF)
            
            temp_vec   = 0d0
            dPlastMult = 0d0

            do k = 1,4
                do j = 1,4
                    temp_vec   = temp_vec+gradF(j)*DEL(j,k)*gradF(k)*A(k)
                    dPlastMult = dPlastMult+gradF(j)*DEL(j,k)*dEps(k)
                end do
            end do
            dPlastMult = max(0d0,dPlastMult/(hard+temp_vec))

        end subroutine compute_plastic_modulus


        subroutine hardening_increments(Sigma_ij, R, X_ij, sigma_yld, &
            b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR, dX_ij)

            ! INCREMENTS OF INTRINSIC STATIC VARIABLES

            real*8, dimension(4), intent(in) :: Sigma_ij        ! actual stress state
            real*8, dimension(4), intent(in) :: X_ij            ! actual back stress state
            real*8,                 intent(in) :: sigma_yld       ! first yielding limit
            real*8,                 intent(in) :: R               ! actual mises radius
            real*8,                 intent(in) :: b_lmc, Rinf_lmc ! Lamaitre and Chaboche parameters (isotropic hardening)
            real*8,                 intent(in) :: C_lmc, kapa_lmc ! Lamaitre and Chaboche parameters (kinematic hardening)

            real*8,                 intent(out):: dR              ! mises radius increment
            real*8, dimension(4), intent(out):: dX_ij           ! back stress increment
            
            real*8                             :: F_mises
            real*8, dimension(4)             :: gradF_mises
            real*8, dimension(4), parameter  :: A = (/1.0,1.0,1.0,0.5/)
            integer*4                        :: k

            ! INCREMENT IN ISOTROPIC HARDENING VARIABLES (R)
            dR = b_lmc*(Rinf_lmc-R)

            ! INCREMENT IN KINEMATIC HARDENING VARIABLES (Xij)
            call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)
            dX_ij=2*A*gradF_mises*C_lmc/3-X_ij*kapa_lmc

        end subroutine hardening_increments

        subroutine drift_corr(stress,center,radius,syld, &
            biso,Rinf,Ckin,kkin,lambda,mu,dEplastic)

            ! DRIFT CORRECTION (RADIAL RETURN)
            real*8, dimension(4), intent(inout) :: stress,center,dEplastic
            real*8,                 intent(inout) :: radius
            real*8,                 intent(in)    :: lambda,mu,syld,biso,Rinf,Ckin,kkin
            real*8                                :: F0,F1,beta,hard,radiust
            real*8, dimension(4)                :: gradF0,gradF1,dstress,stresst,centert
            real*8, dimension(4),     parameter :: A = (/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)            :: DEL
            integer*4                        ::counter,k,j
            
            ! INITIAL PLASTIC CONDITION
            call mises_yld_locus(stress,center,radius,syld,F0,gradF0)
            call stiff_matrix(lambda,mu,DEL)
            do counter=1,5 
                ! COMPUTE HARDENING INCREMENTS
                hard = biso*(Rinf-radius)
                hard = hard + Ckin
                hard = hard - kkin*sum(gradF0*center)
                ! COMPUTE BETA FOR DRIFT CORRECTION
                beta = 0d0
                do j=1,4
                    do k=1,4
                        beta=beta+gradF0(k)*DEL(k,j)*A(j)*gradF0(j)
                    end do
                end do
                beta=F0/(hard+beta)
                ! STRESS-STRAIN-HARDENING CORRECTION
                dstress=0d0
                do k=1,4
                    do j=1,4
                        dstress(k)=dstress(k)-beta*DEL(j,k)*A(j)*gradF0(j)
                    end do
                end do
                stresst = stress+dstress
                centert = center+beta*(2*A*gradF0*Ckin/3-center*kkin)
                radiust = radius+beta*(Rinf-radius)*biso
                
                ! CHECK DRIFT
                call mises_yld_locus(stresst,centert,radiust,syld,F1,gradF1)
                if (abs(F1).gt.abs(F0)) then
                    beta   = F0/sum(gradF0*gradF0)
                    stress = stress-beta*gradF0
                else
                    stress = stresst
                    center = centert
                    radius = radiust
                    dEplastic = dEplastic+beta*A*gradF0
                    if (abs(F1).le.FTOL) then
                        exit
                    else
                        F0     = F1
                        gradF0 = gradF1
                    endif
                endif
                
            enddo
            return
        end subroutine drift_corr

        subroutine gotoFpegasus(start0,dtrial,center,radius,s0,nsub,alpha)
            implicit none
            real*8, dimension(4), intent(in)    :: start0,dtrial,center
            real*8,                 intent(in)    :: radius,s0
            integer*4,              intent(in)    :: nsub
            real*8,                 intent(out)   :: alpha
            real*8, dimension(4)                :: stress0,stress1,stress,gradF
            real*8                                :: dalpha,alpha0,alpha1,F0,F1,FM,Fsave
            integer*4                             :: counter0,counter1
            logical                             :: flagxit

            alpha1  = 1d0
            alpha0  = 0d0
            stress0 = start0+alpha0*dtrial
            stress1 = start0+alpha1*dtrial
            call mises_yld_locus(stress0,center,radius,s0,F0,gradF)
            call mises_yld_locus(stress1,center,radius,s0,F1,gradF)
            if (nsub.gt.1)then
                Fsave=F0
                do counter0=1,10
                    dalpha = (alpha1-alpha0)/nsub
                    do counter1=1,nsub
                        alpha=alpha0+dalpha
                        stress=start0+alpha*dtrial
                        call mises_yld_locus(stress,center,radius,s0,F1,gradF)
                        if (FM.gt.FTOL) then
                            alpha1=alpha
                            if (F0.lt.-FTOL) then
                                F1=FM
                                flagxit=.true.
                                exit
                            else
                                alpha0=0d0
                                F0=Fsave
                            endif
                        else
                            alpha0=alpha
                            F0=FM
                        endif
                    end do
                    if (flagxit) then
                        exit
                    endif
                end do
                if (.not.flagxit) then
                    write(*,*) "ERROR IN FINDING F=0 (REVERSAL)"
                    alpha1=1d0
                    alpha0=0d0
                endif
            end if

            do counter0=1,10
                alpha  = alpha1-F1*(alpha1-alpha0)/(F1-F0)
                stress = start0+alpha*dtrial
                call mises_yld_locus(stress,center,radius,s0,FM,gradF)
                if (abs(FM).le.FTOL) then
                    exit
                else
                    alpha0  =   alpha1
                    alpha1  =   alpha
                    F0      =   F1
                    F1      =   FM
                endif

            end do

            if (FM.gt.FTOL) then
                write(*,*) "WARNING: F>TOL"
            endif
        end subroutine gotoFpegasus

        subroutine tau_deepsoil(gamma_dp,gamma_R,gamma_rev,gamma_max,&
            flag_rev,tau_rev,mu,p1_dp,p2_dp,p3_dp,beta_dp,S_dp,tau_dp)
            real*8,       intent(in)  :: gamma_dp,gamma_R,gamma_max,gamma_rev
            real*8,       intent(in)  :: tau_rev,mu,p1_dp,p2_dp,p3_dp
            logical,    intent(in)  :: flag_rev
            real*8,       intent(out) :: tau_dp
            real*8                    :: Fred,temp

            if (.not.flag_rev) then ! LOADING
                tau_dp=(mu*gamma_dp)/(1+beta_dp*(gamma_dp/gamma_R)**S_dp)
            else                    ! UNLOADING
                tau_dp=tau_rev
                tau_dp=tau_dp+(mu*(gamma_dp-gamma_rev))/&
                    (1+beta_dp*(gamma_max-gamma_R)**S_dp)
                call reduction_factor_deepsoil(gamma_max,p1_dp,p2_dp,p3_dp,mu,Fred)
                temp=(mu*(gamma_dp-gamma_rev))/&
                    (1+beta_dp*(0.5*(gamma_dp-gamma_rev)*gamma_R)**S_dp)
                temp=temp-mu*(gamma_dp-gamma_rev)/(1+beta_dp*(gamma_max-gamma_R)**S_dp)
                tau_dp=tau_dp+Fred*temp
            endif
            
        end subroutine tau_deepsoil
        
        subroutine gamma_ref(stress,A_dp,C_dp,gamma_R)
            real*8, dimension(4), intent(in)    :: stress
            real*8,                 intent(in)    :: C_dp
            real*8,                 intent(out)   :: gamma_R
            real*8                                :: press
            press = sum(stress(1:3))/3 
            gamma_R=A_dp*(press/101325.0)**C_dp
        end subroutine gamma_ref
        
        subroutine reduction_factor_deepsoil(gamma_max,p1_dp,p2_dp,p3_dp,mu,Fred)
            real*8, intent(in) :: gamma_max,p1_dp,p2_dp,p3_dp,mu
            real*8, intent(out):: Fred
            Fred=p1_dp-p2_dp*(1-Ggamma/mu)**p3_dp
        end subroutine
end module nonlinear2d
