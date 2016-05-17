!   Copyright (C) 2014 The SPEED FOUNDATION
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
    real*8, parameter :: FTOL = 0.00100000000D0
    real*8, parameter :: LTOL = 0.0010000000000D0
    real*8, parameter :: STOL = 0.0010000000000000D0
    type nl_element   !< element structure for time loop (RAM saving)
        ! nl stored variables
        real*8, dimension(:,:),   allocatable    :: radius
        real*8, dimension(:,:,:), allocatable    :: stress
        real*8, dimension(:,:,:), allocatable    :: strain
        real*8, dimension(:,:,:), allocatable    :: center
        real*8, dimension(:,:,:), allocatable    :: plastic_strain
        ! elastic parameters
        real*8, dimension(:,:), allocatable      :: lambda,mu
        ! nonlinear parameters
        real*8, dimension(:,:),    allocatable   :: syld,rinf,biso
        real*8, dimension(:,:),    allocatable   :: ckin,kkin
    
    end type 
    contains
        
        !****************************************************************************
        ! MAKE NL INTERNAL FORCES
        !****************************************************************************
        
        subroutine MAKE_INTERNAL_FORCES_NL(nnt,ne,nm,cs_nnz,cs,sdeg_mat,snl,&
            alfa1,alfa2,beta1,beta2,gamma1,gamma2,displ,fk,mvec)
            !
            implicit none
            ! intent IN
            integer*4, intent(in)                       :: nnt,ne,nm,cs_nnz
            integer*4, intent(in), dimension(nm)        :: sdeg_mat
            integer*4, intent(in), dimension(0:cs_nnz)  :: cs
            real*8,    intent(in), dimension(ne)        :: alfa1,beta1,gamma1
            real*8,    intent(in), dimension(ne)        :: alfa2,beta2,gamma2
            real*8,    intent(in), dimension(2*nnt)     :: displ,mvec
            ! intent INOUT
            real*8,    intent(inout), dimension(2*nnt)  :: fk
            type(nl_element), intent(inout), dimension(ne) :: snl 
            ! 
            real*8,     dimension(:),  allocatable      :: ct,ww
            real*8,     dimension(:),  allocatable      :: dxdx,dydy,dxdy,dydx
            real*8,     dimension(:,:),allocatable      :: dd,fx,fy
            real*8,     dimension(:,:,:), allocatable   :: dstrain
            real*8,     dimension(4)                    :: dtrial,dstrain_alpha
            real*8,     dimension(4)                    :: dplastic_strain
            real*8                                      :: alpha_epl,t1ux,t1uy,t2ux
            real*8                                      :: t2uy,t1fx,t1fy,t2fx,t2fy
            integer*4                                   :: ie,ip,iq,il,im,nn,is,in
            logical                                     :: st_epl
            !
            do ie = 1,ne
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                !*********************************************************************************
                ! COMPUTE STRAIN 
                !*********************************************************************************
                
                call ALLOINIT_LOC(ie,nnt,cs,cs_nnz,nn,ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain, &
                    fx,fy,displ,alfa1(ie),alfa2(ie),beta1(ie),beta2(ie),gamma1(ie),gamma2(ie))
                
                dstrain(:,:,:) = dstrain(:,:,:) - snl(ie)%strain(:,:,:)
                
                write(*,*) "dstrain :",dstrain
                write(*,*) "-------"
                write(*,*) "strain0 :",snl(ie)%strain
                read(*,*)
                
                dplastic_strain = 0.d0
                
                !*********************************************************************************
                ! COMPUTE STRESS
                !*********************************************************************************
                
                do iq = 1,nn
                    do ip = 1,nn
                        ! COMPUTE TRIAL STRESS INCREMENT
                        dtrial(:) = 0.d0
                        dstrain_alpha(:) = 0.d0
                        dstrain_alpha(1:2) = dstrain(1:2,ip,iq)
                        dstrain_alpha(4)   = dstrain(3,ip,iq)
                        call MAKE_STRESS_LOC(snl(ie)%lambda(ip,iq),snl(ie)%mu(ip,iq),&
                            dstrain_alpha,dtrial)
                        ! CHECK PLASTICITY
                        call check_plasticity(dtrial,snl(ie)%stress(:,ip,iq),snl(ie)%center(:,ip,iq),&
                            snl(ie)%radius(ip,iq),snl(ie)%syld(ip,iq),st_epl,alpha_epl,ip,iq)
                        ! PLASTIC CORRECTION 
                        if (st_epl) then
                            write(*,*) "plastic step"
                            dstrain_alpha = (1.d0-alpha_epl)*dstrain_alpha
                            call plastic_corrector(dstrain_alpha,dtrial,snl(ie)%center(:,ip,iq),  &
                                snl(ie)%radius(ip,iq),snl(ie)%syld(ip,iq),snl(ie)%biso(ip,iq),&
                                snl(ie)%rinf(ip,iq),snl(ie)%ckin(ip,iq),snl(ie)%kkin(ip,iq),  &
                                snl(ie)%mu(ip,iq),snl(ie)%lambda(ip,iq),dplastic_strain,ie)
                        endif
                        ! STRAIN VECTOR
                        snl(ie)%strain(:,ip,iq)         = snl(ie)%strain(:,ip,iq)         + &
                            dstrain(:,ip,iq)
                        ! PLASTIC STRAIN VECTOR
                        snl(ie)%plastic_strain(:,ip,iq) = snl(ie)%plastic_strain(:,ip,iq) + &
                            (/dplastic_strain(1:2),dplastic_strain(4)/)
                        ! STRESS VECTOR
                        snl(ie)%stress(:,ip,iq)         = dtrial(:)
                    enddo
                enddo
                
                !*********************************************************************************
                ! COMPUTE INTERNAL FORCES
                !*********************************************************************************
                
                call MAKE_INTERNAL_FORCE_EL(nn,ww,dd,dxdx,dxdy,dydx,dydy,snl(ie)%stress(1,:,:),&
                    snl(ie)%stress(2,:,:),snl(ie)%stress(4,:,:),fx,fy)
                do iq = 1,nn
                    do ip = 1,nn
                        is = nn*(iq-1)+ip
                        in = cs(cs(ie-1)+is)
                        fk(in)     = fk(in)     + fx(ip,iq)/mvec(in)
                        fk(in+nnt) = fk(in+nnt) + fy(ip,iq)/mvec(in+nnt)
!                        
!                        t1fx = 0.0d0
!                        t2fx = 0.0d0
!                        t1fy = 0.0d0
!                        t2fy = 0.0d0
!       
!                        ! derivatives with respect to eta ( there is a delta(iq,im) )
!                        do il=1,nn
!                            t1fx = t1fx+dd(il,ip)*ww(il)*ww(iq)*&
!                                (snl(ie)%stress(1,il,iq)*dydy(il)-snl(ie)%stress(4,il,iq)*dxdy(il))
!                            t1fy = t1fy+dd(il,ip)*ww(il)*ww(iq)*&
!                                (snl(ie)%stress(4,il,iq)*dydy(il)-snl(ie)%stress(2,il,iq)*dxdy(il))
!                        enddo
!                        ! derivatives with respect to xi ( there is a delta(ip,il) )
!                        do im=1,nn
!                            t2fx = t2fx+dd(im,iq)*ww(ip)*ww(im)*&
!                                (snl(ie)%stress(1,ip,im)*dydx(im)-snl(ie)%stress(4,ip,im)*dxdx(im))
!                            t2fy = t2fy+dd(im,iq)*ww(ip)*ww(im)*&
!                                (snl(ie)%stress(4,ip,im)*dydx(im)-snl(ie)%stress(2,ip,im)*dxdx(im))
!                        enddo
!                        ! ELEMENT-WISE FORCES
!                        fx(ip,iq) = t1fx - t2fx
!                        fy(ip,iq) = t1fy - t2fy
!                        ! GLOBAL NODAL FORCES
!

                    enddo
                enddo
                ! DEALLOCATE ELEMENT-WISE VARIABLES
                call DEALLOCATE_LOC(ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain,fx,fy)
            enddo
            ! 
            return
            !
        end subroutine MAKE_INTERNAL_FORCES_NL
        
        !****************************************************************************
        ! ALLOCATE LOCAL VARIABLES
        !****************************************************************************
        
        subroutine ALLOINIT_LOC(ie,nnt,cs,cs_nnz,nn,ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain, &
            fx,fy,displ,alfa1,alfa2,beta1,beta2,gamma1,gamma2)
            !
            implicit none    
            ! intent IN
            integer*4,                              intent(in)      :: ie,nnt,cs_nnz,nn
            integer*4, dimension(0:cs_nnz),         intent(in)      :: cs
            real*8, dimension(2*nnt),               intent(in)      :: displ
            real*8,    intent(in)                                   :: alfa1,beta1,gamma1
            real*8,    intent(in)                                   :: alfa2,beta2,gamma2 
            ! intent OUT 
            real*8,     dimension(:),  allocatable, intent(out)     :: ct,ww
            real*8,     dimension(:),  allocatable, intent(out)     :: dxdx,dydy,dxdy,dydx
            real*8,     dimension(:,:),allocatable, intent(out)     :: dd,fx,fy
            real*8,     dimension(:,:,:), allocatable, intent(out)  :: dstrain
            real*8,     dimension(nn,nn)                            :: ux,uy
            real*8,     dimension(nn,nn)                            :: duxdx,duxdy,duydx,duydy
            integer*4                                               :: i,j,is,in
            ! ALLOCATION
            allocate(ct(nn))
            allocate(ww(nn))
            allocate(dd(nn,nn))
            allocate(dxdx(nn))
            allocate(dxdy(nn))
            allocate(dydx(nn))
            allocate(dydy(nn))
            allocate(fx(nn,nn))
            allocate(fy(nn,nn))
            allocate(dstrain(3,nn,nn))
            ! INITIALIZATION
            call lgl(nn,ct,ww,dd)
            ! LOOP OVER GLL
            ux = 0.0d0
            uy = 0.0d0
            duxdx = 0.d0
            duxdy = 0.d0
            duydx = 0.d0
            duydy = 0.d0
            fx = 0.d0
            fy = 0.d0
            dxdx = 0.d0
            dxdy = 0.d0
            dydx = 0.d0
            dydy = 0.d0
            dstrain(:,:,:) = 0.d0
            do j = 1,nn
                do i = 1,nn
                    is = nn*(j -1) +i
                    in = cs(cs(ie -1) + is)
                    ux(i,j) = displ(in)
                    uy(i,j) = displ(in+nnt)
                enddo
            enddo
            ! MAKE DERIVATIVES
            call MAKE_DERIVATIVES_LOC(nn,alfa1,alfa2,beta1,beta2,gamma1,gamma2,ct,&
                dxdx,dxdy,dydx,dydy)
             
            ! MAKE STRAIN
            call MAKE_STRAIN(nn,dd,dxdx,dxdy,dydx,dydy,ux,uy,duxdx,duxdy,duydx,duydy)
            do j = 1,nn
                do i = 1,nn
                    dstrain(1,i,j) = duxdx(i,j)
                    dstrain(2,i,j) = duydy(i,j)
                    dstrain(3,i,j) = duxdy(i,j)+duydx(i,j)
                enddo
            enddo

            return
        end subroutine ALLOINIT_LOC
        
        !****************************************************************************
        ! MAKE DERIVATIVES
        !****************************************************************************
        
        subroutine MAKE_DERIVATIVES_LOC(nn,alfa1,alfa2,beta1,beta2,gamma1,gamma2,ct,&
            dxdx,dxdy,dydx,dydy)
            !           
            implicit none
            ! intent IN
            integer*4, intent(in)               :: nn
            real*8, intent(in)                  :: alfa1,alfa2
            real*8, intent(in)                  :: beta1,beta2
            real*8, intent(in)                  :: gamma1,gamma2
            real*8, intent(in)   , dimension(nn):: ct
            ! intent INOUT
            real*8, intent(inout), dimension(nn):: dxdx,dxdy
            real*8, intent(inout), dimension(nn):: dydx,dydy
            ! 
            integer*4 :: i
            ! COMPUTE SHAPE-FUNCTION DERIVATIVES
            do i = 1,nn
                dxdy(i) = beta1 + gamma1 * ct(i)
                dydy(i) = beta2 + gamma2 * ct(i)
                dxdx(i) = alfa1 + gamma1 * ct(i)
                dydx(i) = alfa2 + gamma2 * ct(i)
            enddo
            ! 
            return
        end subroutine MAKE_DERIVATIVES_LOC
        
        !****************************************************************************
        ! COMPUTE ELEMENT-WISE STRAIN TENSOR
        !****************************************************************************
        
!        subroutine MAKE_STRAIN_LOC(nn,dd,ux,uy,dxdx,dxdy,dydx,dydy,dstrain)
!            ! 
!            implicit none
!            ! intent IN
!            integer*4, intent(in)                       :: nn
!            real*8, dimension(nn),    intent(in)        :: dxdx,dxdy,dydx,dydy
!            real*8, dimension(nn,nn), intent(in)        :: dd,ux,uy
!            ! intent INOUT
!            real*8, dimension(3,nn,nn), intent(inout)   :: dstrain
!            real*8                                      :: t1ux,t1uy,t2ux,t2uy
!            real*8                                      :: t1fx,t1fy,t2fx,t2fy,det_j
!            !
!            integer*4                                   :: ip,iq,il,im
!            ! COMPUTE STRAIN
!            do iq = 1,nn
!                do ip = 1,nn
!                    t1ux = 0.d0
!                    t1uy = 0.d0
!                    t2ux = 0.d0
!                    t2uy = 0.d0
!                    
!                    det_j = dxdx(iq)*dydy(ip)-dxdy(ip)*dydx(iq)
!                    
!                    do il = 1,nn
!                       t1ux = t1ux + ux(il,iq) * dd(ip,il)
!                       t1uy = t1uy + uy(il,iq) * dd(ip,il)
!                    enddo
!
!                    do im = 1,nn
!                       t2ux = t2ux + ux(ip,im) * dd(iq,im)
!                       t2uy = t2uy + uy(ip,im) * dd(iq,im)
!                    enddo
!                    dstrain(1,ip,iq) = ( 1.d0 / det_j) * ((dydy(ip) * t1ux) - (dydx(iq) * t2ux))
!                    dstrain(2,ip,iq) = (-1.d0 / det_j) * ((dxdy(ip) * t1uy) - (dxdx(iq) * t2uy))
!                    dstrain(3,ip,iq) = ( 1.d0 / det_j) * ((dydy(ip) * t1uy) - (dydx(iq) * t2uy))
!                    dstrain(3,ip,iq) = dstrain(3,ip,iq)+(-1.d0 / det_j) * &
!                        ((dxdy(ip) * t1ux) - (dxdx(iq) * t2ux))
!                enddo
!            enddo
!            !
!            return
!        end subroutine MAKE_STRAIN_LOC
        
        !****************************************************************************
        ! DEALLOCATE LOCAL VARIABLES
        !****************************************************************************
        
        subroutine DEALLOCATE_LOC(ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain,fx,fy)
            !
            implicit none    
            ! intent INOUT 
            real*8,     dimension(:),  allocatable, intent(inout) :: ct,ww
            real*8,     dimension(:),  allocatable, intent(inout) :: dxdx,dydy,dxdy,dydx
            real*8,     dimension(:,:),allocatable, intent(inout) :: dd,fx,fy
            real*8,     dimension(:,:,:), allocatable, intent(inout) :: dstrain
            ! DEALLOCATION
            deallocate(ct)
            deallocate(ww)
            deallocate(dd)
            deallocate(dxdx)
            deallocate(dxdy)
            deallocate(dydx)
            deallocate(dydy)
            deallocate(fx)
            deallocate(fy)
            deallocate(dstrain)
            !
            return
        end subroutine DEALLOCATE_LOC
        
        subroutine STIFF_MATRIX(lambda,mu,DEL)
            implicit none
            ! intent IN
            real*8, intent(in) :: lambda,mu
            ! intent OUT
            real*8, intent(out), dimension(4,4) :: DEL
            !
            DEL(:,:) = 0.d0
            DEL(1,1) = 2*mu
            DEL(2,2) = 2*mu
            DEL(3,3) = 2*mu
            DEL(4,4) = mu
            DEL(1:3,1:3) = DEL(1:3,1:3) + lambda
            !
            return
        end subroutine STIFF_MATRIX

        !****************************************************************************
        ! COMPUTE STRESS STATE
        !****************************************************************************
        
        subroutine MAKE_STRESS_LOC(lambda,mu,dE,dS)
            !
            implicit none
            ! input IN
            real*8, intent(in)                  :: lambda,mu
            real*8, intent(in), dimension(4)    :: dE
            ! input INOUT
            real*8, intent(inout), dimension(4) :: dS
            !
            real*8, dimension(4,4)              :: DEL
            integer*4                           :: ip 
            ! STIFFNESS MATRIX (PLANE STRAIN VERSION)
            call STIFF_MATRIX(lambda,mu,DEL)
            ! COMPUTE STRESS 
            do ip=1,4
                dS(ip) = DEL(ip,1)*dE(1)+DEL(ip,2)*dE(2)+DEL(ip,3)*dE(3)+DEL(ip,4)*dE(4)
            enddo
            return
            !
        end subroutine MAKE_STRESS_LOC
        
        !****************************************************************************
        ! MISES YIELDING LOCUS AND GRADIENT
        !****************************************************************************

        subroutine mises_yld_locus(stress, center, radius, syld, FM, gradF)
            ! 
            implicit none
            ! intent IN
            real*8,               intent(in)    :: radius,syld
            real*8, dimension(4), intent(in)    :: stress,center
            ! intent OUT
            real*8, dimension(4), intent(out)   :: gradF
            real*8,               intent(out)   :: FM
            !
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,2.0/)
            real*8, dimension(4)                :: dev
            real*8                              :: tau_eq
            integer*4                           :: k
            ! COMPUTE TENSOR COMPONENTS 
            call tensor_components (stress, dev)
            ! COMPUTE MISES FUNCTION
            call tau_mises(dev-center, tau_eq)
            FM = tau_eq-syld-radius
            ! COMPUTE MISES FUNCTION GRADIENT
            gradF = 1.5d0*A*(dev-center)/tau_eq
            !  
            return
            !
        end subroutine mises_yld_locus
        
        !****************************************************************************
        ! OCTAHEDRAL SHEAR STRESS
        !****************************************************************************
        
        subroutine tau_mises(stress,tau_eq)
            !
            implicit none
            ! intent IN
            real*8, dimension(4), intent(in)    :: stress
            ! intent OUT
            real*8,               intent(out)   :: tau_eq
            !
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,2.0/) 
            integer*4                           :: k
            ! COMPUTE OCTAHEDRAL SHEAR STRESS 
            tau_eq = 0.0d0
            do k=1,4
                tau_eq = tau_eq+A(k)*(stress(k)**2)
            end do
            tau_eq = sqrt(1.5*tau_eq)
            !
            return
            ! 
        end subroutine
        
        !****************************************************************************
        ! TENSOR COMPONENTS (SPHERICAL & DEVIATORIC)
        !****************************************************************************
        
        subroutine tensor_components(stress, dev)
            !     
            implicit none
            ! intent IN
            real*8, dimension(4), intent(in)    :: stress
            ! intent OUT
            real*8, dimension(4), intent(out)   :: dev
            !
            real*8                              :: press 
            !
            dev = stress
            press = sum(stress(1:3))/3
            dev(1:3) = dev(1:3)-press
            !
            return
            !
        end subroutine tensor_components
            
        !****************************************************************************
        ! CHECK PLASTIC CONSISTENCY (KKT CONDITIONS)
        !****************************************************************************
        
        subroutine check_plasticity (dtrial, stress0, center, radius, syld, &
            st_elp, alpha_epl, ip,iq)
            !     
            implicit none
            ! intent IN
            integer*4, intent(in)                   :: ip,iq 
            real*8,                 intent(in)      :: radius,syld
            real*8, dimension(4),   intent(in)      :: center,stress0
            ! intent INOUT
            real*8, dimension(4),   intent(inout)   :: dtrial
            ! intent OUT
            real*8,                 intent(out)     :: alpha_epl
            logical,                intent(out)     :: st_elp
            !
            integer*4                               :: k
            real*8                                  :: FS, FT,checkload
            real*8, dimension(4)                    :: gradFS, gradFT, stress1
            !
            stress1= dtrial+stress0
            call mises_yld_locus(stress0,center,radius,syld,FS,gradFS)
            call mises_yld_locus(stress1,center,radius,syld,FT,gradFT)
            checkload=sum(gradFS*dtrial)/sum(gradFS**2)/sum(dtrial**2)

            if (abs(FS).le.FTOL) then
                if (checkload.ge.-LTOL) then
                    alpha_epl = 0d0
                    st_elp    = .true.
                else
                    if (FT.lt.-FTOL) then
                        alpha_epl = 1d0
                        st_elp    = .false.
                    elseif(FT.gt.FTOL) then
                        call gotoFpegasus(stress0,dtrial,center,radius,syld,10,alpha_epl)
                        st_elp    = .true.
                    endif
                end if
            elseif (FS.lt.-FTOL) then
                if (FT.le.FTOL) then
                    alpha_epl = 1d0
                    st_elp    = .false.
                else
                    call gotoFpegasus(stress0,dtrial,center,radius,syld,1,alpha_epl)
                    st_elp    = .true.
                end if
            elseif (FS.gt.FTOL) then
                write(*,*) 'ERROR: FS: ',FS,'>',FTOL,'!!!!'
                write(*,*) 'stress0: ',stress0
                write(*,*) 'dtrial: ',dtrial
                write(*,*) 'center: ',center
                write(*,*) 'syld: ',syld
                write(*,*) 'glls: ',ip,iq
                stop
                alpha_epl  = 0d0
                st_elp = .true.
            end if
            dtrial=stress0+dtrial*alpha_epl
            call mises_yld_locus(dtrial,center,radius,syld,FS,gradFS)
        end subroutine check_plasticity
        
        subroutine plastic_corrector(dEps_alpha,stress,center,syld, &
            radius,biso,Rinf,Ckin,kkin,mu,lambda,dEpl,nel)

            implicit none
            integer*4, intent(in) :: nel
            real*8, dimension(4), intent(inout) :: dEps_alpha,stress,center,dEpl 
            real*8,               intent(inout) :: radius 
            real*8,               intent(in)    :: syld,biso,Rinf,Ckin,kkin,mu,lambda
            real*8, dimension(4,4)              :: DEL
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
                Resk  = 0.0d0
                dS1   = 0.0d0
                dX1   = 0.0d0
                dS2   = 0.0d0
                dX2   = 0.0d0
                dR1   = 0.0d0
                dR2   = 0.0d0
                dEpl1 = 0.0d0
                dEpl2 = 0.0d0

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
                    dEpl   = dEpl+dEpl1
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
            !
            implicit none
            ! intent IN
            real*8              , intent(in)    :: radius,syld,mu,lambda,biso,Rinf,Ckin,kapakin
            real*8, dimension(4), intent(in)    :: dStrain,Stress,center
            ! intent INOUT
            real*8, dimension(4), intent(inout) :: dStress,dcenter,dEplast
            real*8,               intent(inout) :: dradius
            ! intent OUT
            real*8,               intent(out)   :: hard
            !
            real*8, dimension(4)                :: gradF
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)              :: DEL
            real*8                              :: Fmises,dPlast
            integer*4                           :: j,k
            
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
            !
            implicit none
            ! intent IN
            real*8,                 intent(in) :: mu,lambda,syld   
            real*8,                 intent(in) :: radius,biso,Rinf,Ckin,kkin
            real*8, dimension(4),   intent(in) :: dEps,stress,center
            ! intent OUT
            real*8,                 intent(out):: dPlastMult,hard
            !
            real*8                             :: temp_vec,FM
            real*8, dimension(4)               :: gradF
            real*8, dimension(4),   parameter  :: A =(/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)             :: DEL
            integer*4                          :: j,k
            
            call stiff_matrix(lambda,mu,DEL)
            call mises_yld_locus(stress,center,radius,syld,FM,gradF)
            
            hard = Ckin+biso*(Rinf-radius)
            hard = hard-kkin*sum(center*gradF)
            
            temp_vec   = 0.0d0
            dPlastMult = 0.0d0

            do k = 1,4
                do j = 1,4
                    temp_vec   = temp_vec+gradF(j)*DEL(j,k)*gradF(k)*A(k)
                    dPlastMult = dPlastMult+gradF(j)*DEL(j,k)*dEps(k)
                end do
            end do
            dPlastMult = max(0.0d0,dPlastMult/(hard+temp_vec))

        end subroutine compute_plastic_modulus


        subroutine hardening_increments(Sigma_ij, R, X_ij, sigma_yld, &
            b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR, dX_ij)

            ! INCREMENTS OF INTRINSIC STATIC VARIABLES

            real*8,               intent(in) :: sigma_yld       ! first yielding limit
            real*8,               intent(in) :: R               ! actual mises radius
            real*8,               intent(in) :: b_lmc, Rinf_lmc ! Lamaitre and Chaboche parameters (isotropic hardening)
            real*8,               intent(in) :: C_lmc, kapa_lmc ! Lamaitre and Chaboche parameters (kinematic hardening)
            real*8, dimension(4), intent(in) :: Sigma_ij        ! actual stress state
            real*8, dimension(4), intent(in) :: X_ij            ! actual back stress state

            real*8,               intent(out):: dR              ! mises radius increment
            real*8, dimension(4), intent(out):: dX_ij           ! back stress increment
            
            real*8                           :: F_mises
            real*8, dimension(4)             :: gradF_mises
            real*8, dimension(4), parameter  :: A = (/1.0,1.0,1.0,0.5/)
            integer*4                        :: k

            ! INCREMENT IN ISOTROPIC HARDENING VARIABLES (R)
            dR = b_lmc*(Rinf_lmc-R)

            ! INCREMENT IN KINEMATIC HARDENING VARIABLES (Xij)
            call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)
            dX_ij(:)=2*A(:)*gradF_mises(:)*C_lmc/3-X_ij(:)*kapa_lmc

        end subroutine hardening_increments

        subroutine drift_corr(stress,center,radius,syld, &
            biso,Rinf,Ckin,kkin,lambda,mu,dEplastic)

            ! DRIFT CORRECTION (RADIAL RETURN)
            real*8, dimension(4), intent(inout) :: stress,center,dEplastic
            real*8,               intent(inout) :: radius
            real*8,               intent(in)    :: lambda,mu,syld,biso,Rinf,Ckin,kkin
            real*8                              :: F0,F1,beta,hard,radiust
            real*8, dimension(4)                :: gradF0,gradF1,dstress,stresst,centert
            real*8, dimension(4),     parameter :: A = (/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)              :: DEL
            integer*4                           :: counter,k,j
            
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
            ! intent IN
            real*8, dimension(4), intent(in)    :: start0,dtrial,center
            real*8,               intent(in)    :: radius,s0
            integer*4,            intent(in)    :: nsub
            ! intent OUT
            real*8,               intent(out)   :: alpha
            real*8, dimension(4)                :: stress0,stress1,stress,gradF
            real*8                              :: dalpha,alpha0,alpha1,F0,F1,FM,Fsave
            integer*4                           :: counter0,counter1
            logical                             :: flagxit

            alpha1  = 1.0d0
            alpha0  = 0.0d0
            stress0 = start0+alpha0*dtrial
            stress1 = start0+alpha1*dtrial
            call mises_yld_locus(stress0,center,radius,s0,F0,gradF)
            call mises_yld_locus(stress1,center,radius,s0,F1,gradF)
            if (nsub.gt.1)then
                Fsave=F0
                do counter0=1,3
                    dalpha = (alpha1-alpha0)/nsub
                    flagxit=.false.
                    do counter1=1,nsub
                        alpha=alpha0+dalpha
                        stress=start0+alpha*dtrial
                        call mises_yld_locus(stress,center,radius,s0,FM,gradF)
                        if (FM.gt.FTOL) then
                            alpha1=alpha
                            if (F0.lt.-FTOL) then
                                F1=FM
                                flagxit=.true.
                            else
                                alpha0=0.0d0
                                F0=Fsave
                            endif
                            exit
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
                    alpha1=1.0d0
                    alpha0=0.0d0
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
                read(*,*)
            endif
        end subroutine gotoFpegasus

!        subroutine tau_deepsoil(gamma_dp,gamma_R,gamma_rev,gamma_max,&
!            flag_rev,tau_rev,mu,p1_dp,p2_dp,p3_dp,beta_dp,S_dp,tau_dp)
!            real*8,       intent(in)  :: gamma_dp,gamma_R,gamma_max,gamma_rev
!            real*8,       intent(in)  :: tau_rev,mu,p1_dp,p2_dp,p3_dp
!            logical,    intent(in)  :: flag_rev
!            real*8,       intent(out) :: tau_dp
!            real*8                    :: Fred,temp
!
!            if (.not.flag_rev) then ! LOADING
!                tau_dp=(mu*gamma_dp)/(1+beta_dp*(gamma_dp/gamma_R)**S_dp)
!            else                    ! UNLOADING
!                tau_dp=tau_rev
!                tau_dp=tau_dp+(mu*(gamma_dp-gamma_rev))/&
!                    (1+beta_dp*(gamma_max-gamma_R)**S_dp)
!                call reduction_factor_deepsoil(gamma_max,p1_dp,p2_dp,p3_dp,mu,Fred)
!                temp=(mu*(gamma_dp-gamma_rev))/&
!                    (1+beta_dp*(0.5*(gamma_dp-gamma_rev)*gamma_R)**S_dp)
!                temp=temp-mu*(gamma_dp-gamma_rev)/(1+beta_dp*(gamma_max-gamma_R)**S_dp)
!                tau_dp=tau_dp+Fred*temp
!            endif
!            
!        end subroutine tau_deepsoil
!        
!        subroutine gamma_ref(stress,A_dp,C_dp,gamma_R)
!            real*8, dimension(4), intent(in)    :: stress
!            real*8,                 intent(in)    :: C_dp
!            real*8,                 intent(out)   :: gamma_R
!            real*8                                :: press
!            press = sum(stress(1:3))/3 
!            gamma_R=A_dp*(press/101325.0)**C_dp
!        end subroutine gamma_ref
!        
!        subroutine reduction_factor_deepsoil(gamma_max,p1_dp,p2_dp,p3_dp,mu,Fred)
!            real*8, intent(in) :: gamma_max,p1_dp,p2_dp,p3_dp,mu
!            real*8, intent(out):: Fred
!            Fred=p1_dp-p2_dp*(1-Ggamma/mu)**p3_dp
!        end subroutine
end module nonlinear2d

module seismic
    
    use nonlinear2d
    implicit none
    contains
        
        subroutine MAKE_SEISMIC_FORCES(nnt,nm,ne,nf,cs_nnz,cs,sdeg_mat,nfunc_data,nl_sism,length_cns,&
            check_node_sism,check_dist_node_sism,func_data,func_type,tag_func,func_indx,facsmom,tt1,&
            snl,alfa1,alfa2,beta1,beta2,gamma1,gamma2,displ,mvec,sism,fe)
            !
            implicit none
            ! INTENT IN
            integer*4, intent(in)                               :: nnt,nm,ne,nf,cs_nnz
            integer*4, intent(in)                               :: nfunc_data,nl_sism,length_cns
            integer*4, intent(in), dimension(nm)                :: sdeg_mat
            integer*4, intent(in), dimension(nf)                :: func_type,tag_func
            integer*4, intent(in), dimension(nf+1)              :: func_indx
            integer*4, intent(in), dimension(0:cs_nnz)          :: cs
            integer*4, intent(in), dimension(length_cns,5)      :: check_node_sism
            real*8,    intent(in)                               :: tt1
            real*8,    intent(in), dimension(ne)                :: alfa1,beta1,gamma1
            real*8,    intent(in), dimension(ne)                :: alfa2,beta2,gamma2
            real*8,    intent(in), dimension(2*nnt)             :: displ,mvec
            real*8,    intent(in), dimension(nl_sism,3)         :: facsmom       
            real*8,    intent(in), dimension(nfunc_data)        :: func_data         
            real*8,    intent(in), dimension(length_cns,1)      :: check_dist_node_sism
            ! INTENT INOUT
            type(nl_element), intent(inout), dimension(ne)      :: snl 
            real*8,           intent(inout), dimension(2*nnt)   :: sism,fe
            ! 
            real*8,     dimension(:),  allocatable      :: ct,ww
            real*8,     dimension(:),  allocatable      :: dxdx,dydy,dxdy,dydx
            real*8,     dimension(:,:),allocatable      :: dd,det_j,fxs,fys
            real*8,     dimension(:,:),allocatable      :: sxxs,syys,szzs,sxys
            real*8,     dimension(:,:,:), allocatable   :: dstrain
            real*8                                      :: t1ux,t1uy,t2ux
            real*8                                      :: t2uy,t1fx,t1fy,t2fx,t2fy
            integer*4                                   :: ie,ip,iq,il,im,nn,is,in
            !
            fe = 0.d0
            do ie = 1,ne 
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                
                ! ALLOCATION/INITIALIZATION LOCAL VARIABLES 
                call ALLOINIT_LOC(ie,nnt,cs,cs_nnz,nn,ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain, &
                    fxs,fys,displ,alfa1(ie),alfa2(ie),beta1(ie),beta2(ie),gamma1(ie),gamma2(ie))
                
                allocate(sxxs(nn,nn))
                allocate(syys(nn,nn))
                allocate(szzs(nn,nn))
                allocate(sxys(nn,nn))
                sxxs = 0.d0
                syys = 0.d0
                sxys = 0.d0
                szzs = 0.d0
                ! LOOP OVER GLL
                do iq = 1,nn
                    do ip = 1,nn
                        call MAKE_SEISMIC_MOMENT_NEW(nn,sxxs,syys,szzs,sxys,&
                            check_node_sism,check_dist_node_sism,length_cns,ie,facsmom, &
                            nl_sism,func_type,func_indx,func_data,nf,tt1, &
                            nfunc_data,tag_func)
                        call MAKE_INTERNAL_FORCE_EL(nn,ww,dd,dxdx,dxdy,dydx,dydy,&
                            sxxs,syys,sxys,fxs,fys)
                    enddo
                enddo
                ! ASSIGN TO GLOBAL NODAL EXTERNAL FORCES
                do iq = 1,nn
                    do ip = 1,nn
                        is = nn*(iq-1)+ip
                        in = cs(cs(ie-1)+is)
                        sism(in)    = sism(in)      + fxs(ip,iq)
                        sism(in+nnt)= sism(in+nnt)  + fys(ip,iq)
                        fe(in)      = fe(in)        + sism(in)/mvec(in)
                        fe(in+nnt)  = fe(in+nnt)    + sism(in+nnt)/mvec(in+nnt)  
                    enddo
                enddo
            enddo
    end subroutine 
end module seismic

module write_output

    use nonlinear2d
    implicit none
    
    type nodepatched
        ! nl stored variables
        real*8, dimension(:,:), allocatable    :: stress
        real*8, dimension(:,:), allocatable    :: strain
        real*8, dimension(:,:), allocatable    :: plastic_strain
    end type nodepatched 
    contains
        subroutine MAKE_MONITOR_FILES(nnt,nmonit,option_out_var,nnode_TOT,tagstep,&
                unit_disp,unit_vel,unit_acc,unit_stress,unit_strain,unit_omega,unit_uDRM)

            implicit none

            integer*4 ,                   intent(in)    :: nnt,nmonit,nnode_TOT,tagstep
            integer*4, dimension(nmonit), intent(inout) :: unit_disp
            integer*4, dimension(nmonit), intent(inout) :: unit_vel
            integer*4, dimension(nmonit), intent(inout) :: unit_acc
            integer*4, dimension(nmonit), intent(inout) :: unit_stress
            integer*4, dimension(nmonit), intent(inout) :: unit_strain
            integer*4, dimension(nmonit), intent(inout) :: unit_omega
            integer*4, dimension(nnode_TOT), intent(inout) :: unit_uDRM
            integer*4, dimension(6)     , intent(in)    :: option_out_var    
            character*70                                :: file_disp
            character*70                                :: file_vel
            character*70                                :: file_acc
            character*70                                :: file_stress 
            character*70                                :: file_strain 
            character*70                                :: file_omega 
            character*70                                :: file_uDRM
            integer*4                                   :: i 
            if (nmonit.ge.1) then
                if (option_out_var(1).eq.1) then  
                    file_disp = 'monitorXXXXX.d'  
                    do i = 1,nmonit
                        unit_disp(i)=50000+i
                        if (i.lt.10) then
                            write(file_disp(8:11),'(a4)')'0000'
                            write(file_disp(12:12),'(i1)')i
                        else if (i.le.99) then
                            write(file_disp(8:10),'(a3)')'000'
                            write(file_disp(11:12),'(i2)')i
                        else if (i.le.999) then
                            write(file_disp(8:9),'(a2)')'00'
                            write(file_disp(10:12),'(i3)')i
                        else if (i.le.9999) then
                            write(file_disp(8:8),'(a1)')'0'
                            write(file_disp(9:12),'(i4)')i    
                        else if (i.le.99999) then
                            write(file_disp(8:12),'(i5)')'i'     
                        endif
                        open(unit_disp(i),file=file_disp)
                    enddo
                endif

                if (option_out_var(2).eq.1) then  
                    file_vel = 'monitorXXXXX.v' 
                    do i = 1,nmonit
                        unit_vel(i)=100000+i
                        if (i.lt.10) then
                            write(file_vel(8:11),'(a4)')'0000'
                            write(file_vel(12:12),'(i1)')i
                        else if (i.le.99) then
                            write(file_vel(8:10),'(a3)')'000'
                            write(file_vel(11:12),'(i2)')i
                        else if (i.le.999) then
                            write(file_vel(8:9),'(a2)')'00'
                            write(file_vel(10:12),'(i3)')i
                        else if (i.le.9999) then
                            write(file_vel(8:8),'(a1)')'0'
                            write(file_vel(9:12),'(i4)')i    
                        else if (i.le.99999) then
                            write(file_vel(8:12),'(i5)')'i'     
                        endif
                        open(unit_vel(i),file=file_vel)
                    enddo
                endif


                if (option_out_var(3).eq.1) then  

                    file_acc = 'monitorXXXXX.a' 
                    do i = 1,nmonit
                        unit_acc(i)=200000+i
                        if (i.lt.10) then
                            write(file_acc(8:11),'(a4)')'0000'
                            write(file_acc(12:12),'(i1)')i
                        else if (i.le.99) then
                            write(file_acc(8:10),'(a3)')'000'
                            write(file_acc(11:12),'(i2)')i
                        else if (i.le.999) then
                            write(file_acc(8:9),'(a2)')'00'
                            write(file_acc(10:12),'(i3)')i
                        else if (i.le.9999) then
                            write(file_acc(8:8),'(a1)')'0'
                            write(file_acc(9:12),'(i4)')i    
                        else if (i.le.99999) then
                            write(file_acc(8:12),'(i5)')'i'     
                        endif
                        open(unit_acc(i),file=file_acc)
                    enddo 
                endif

                if (option_out_var(4).eq.1) then
                    file_stress = 'monitorXXXXX.s'
                    do i = 1,nmonit
                        unit_stress(i)=300000+i
                        if (i.lt.10) then
                            write(file_stress(8:11),'(a4)')'0000'
                            write(file_stress(12:12),'(i1)')i
                        else if (i.le.99) then
                            write(file_stress(8:10),'(a3)')'000'
                            write(file_stress(11:12),'(i2)')i
                        else if (i.le.999) then
                            write(file_stress(8:9),'(a2)')'00'
                            write(file_stress(10:12),'(i3)')i
                        else if (i.le.9999) then
                            write(file_stress(8:8),'(a1)')'0'
                            write(file_stress(9:12),'(i4)')i    
                        else if (i.le.99999) then
                            write(file_stress(8:12),'(i5)')'i'     
                        endif         
                        open(unit_stress(i),file=file_stress)
                    enddo
                endif

                if (option_out_var(5).eq.1) then  
                    file_strain = 'monitorXXXXX.e'
                    do i = 1,nmonit
                        unit_strain(i)=400000+i
                        if (i.lt.10) then
                            write(file_strain(8:11),'(a4)')'0000'
                            write(file_strain(12:12),'(i1)')i
                        else if (i.le.99) then
                            write(file_strain(8:10),'(a3)')'000'
                            write(file_strain(11:12),'(i2)')i
                        else if (i.le.999) then
                            write(file_strain(8:9),'(a2)')'00'
                            write(file_strain(10:12),'(i3)')i
                        else if (i.le.9999) then
                            write(file_strain(8:8),'(a1)')'0'
                            write(file_strain(9:12),'(i4)')i    
                        else if (i.le.99999) then
                            write(file_strain(8:12),'(i5)')'i'     
                        endif                   
                        open(unit_strain(i),file=file_strain)
                    enddo
                endif


                if (option_out_var(6).eq.1) then  
                    file_omega = 'monitorXXXXX.w'
                    do i = 1,nmonit
                        unit_omega(i)=500000+i
                        if (i.lt.10) then
                            write(file_omega(8:11),'(a4)')'0000'
                            write(file_omega(12:12),'(i1)')i
                        else if (i.le.99) then
                            write(file_omega(8:10),'(a3)')'000'
                            write(file_omega(11:12),'(i2)')i
                        else if (i.le.999) then
                            write(file_omega(8:9),'(a2)')'00'
                            write(file_omega(10:12),'(i3)')i
                        else if (i.le.9999) then
                            write(file_omega(8:8),'(a1)')'0'
                            write(file_omega(9:12),'(i4)')i    
                        else if (i.le.99999) then
                            write(file_omega(8:12),'(i5)')'i'     
                        endif                   
                        open(unit_omega(i),file=file_omega)
                    enddo
                endif

            endif

            !---DRM---------------------------------------------------------------------
            !Open output files for DRM I step                     !DRM Scandella 25.11.2005

            if ((nnode_TOT.ne.0).and.(tagstep.eq.1)) then   !DRM Scandella 25.11.2005 
                file_uDRM = 'monDRMXXXXX.d'                  !DRM Scandella 25.11.2005 
                do i = 1,nnode_TOT                           !DRM Scandella 25.11.2005 
                    unit_uDRM(i)=600000+i                    !DRM Scandella 25.11.2005
                    if (i.lt.10) then                      !DRM Scandella 25.11.2005 
                        write(file_uDRM(7:10),'(a4)')'0000'     !DRM Scandella 25.11.2005
                        write(file_uDRM(11:11),'(i1)')i     !DRM Scandella 25.11.2005
                    else if (i.le.99) then                 !DRM Scandella 25.11.2005
                        write(file_uDRM(7:9),'(a3)')'000'     !DRM Scandella 25.11.2005
                        write(file_uDRM(10:11),'(i2)')i     !DRM Scandella 25.11.2005
                    else if (i.le.999) then                !DRM Scandella 25.11.2005
                        write(file_uDRM(7:8),'(a2)')'00'     !DRM Scandella 25.11.2005
                        write(file_uDRM(9:11),'(i3)')i      !DRM Scandella 25.11.2005
                    else if (i.le.9999) then               !DRM Scandella 25.11.2005
                        write(file_uDRM(7:7),'(a1)')'0'     !DRM Scandella 25.11.2005 
                        write(file_uDRM(8:11),'(i4)')i      !DRM Scandella 25.11.2005
                    else if (i.le.99999) then              !DRM Scandella 25.11.2005
                        write(file_uDRM(7:11),'(i5)')i      !DRM Scandella 25.11.2005  
                    endif                                  !DRM Scandella 25.11.2005              
                    open(unit_uDRM(i),file=file_uDRM)         !DRM Scandella 25.11.2005
                    !close(unit_uDRM)                      !DRM Scandella 13.12.2005 			       
                enddo                                        !DRM Scandella 25.11.2005
            endif                                           !DRM Scandella 25.11.2005
            return
        end  subroutine MAKE_MONITOR_FILES

        subroutine WRITE_MONITOR_EL(cs_nnz,cs,ne,nm,sdeg_mat,prop_mat,nmonit,option_out_var, &
                unit_disp,unit_vel,unit_acc,unit_strain,unit_stress,unit_omega,node_m,nnt,         &
                tt1,dis,vel,acc,nodal_counter,alfa1,alfa2,beta1,beta2,gamma1,gamma2)
            
            implicit none
            real*8,                         intent(in)      :: tt1
            integer*4,                      intent(in)      :: nmonit,nnt,nm,ne,cs_nnz
            integer*4,  dimension(nm),      intent(in)      :: sdeg_mat
            integer*4,  dimension(0:cs_nnz),intent(in)      :: cs
            real*8,     dimension(2*nnt) ,  intent(inout)   :: dis,vel,acc
            real*8,     dimension(nm,8)  ,  intent(in)      :: prop_mat
            real*8,     dimension(ne)    ,  intent(in)      :: alfa1,beta1,gamma1
            real*8,     dimension(ne)    ,  intent(in)      :: alfa2,beta2,gamma2 
            integer*4,  dimension(6),       intent(in)      :: option_out_var
            integer*4,  dimension(nnt),     intent(in)      :: nodal_counter
            integer*4,  dimension(nmonit),  intent(in)      :: node_m
            integer*4,  dimension(nmonit),  intent(inout)   :: unit_disp
            integer*4,  dimension(nmonit),  intent(inout)   :: unit_vel
            integer*4,  dimension(nmonit),  intent(inout)   :: unit_acc
            integer*4,  dimension(nmonit),  intent(inout)   :: unit_stress
            integer*4,  dimension(nmonit),  intent(inout)   :: unit_strain
            integer*4,  dimension(nmonit),  intent(inout)   :: unit_omega
            integer*4                                       :: ii,jj,i,ie,im,is,in,nn

            real*8, dimension(:),   allocatable :: ct,ww,dxdx_el,dxdy_el,dydx_el,dydy_el
            real*8, dimension(:,:), allocatable :: dd,ux_el,uy_el
            real*8, dimension(:,:), allocatable :: duxdx_el,duxdy_el,duydx_el,duydy_el
            real*8, dimension(:,:), allocatable :: sxx_el,syy_el,szz_el,sxy_el
            real*8, dimension(:,:), allocatable :: det_j,mu_el,lambda_el 
            real*8, dimension(nnt)              :: sxx,syy,sxy,szz 
            real*8, dimension(nnt)              :: duxdx,duydy,duxdy,duydx
            
           ! Displacements
            if (option_out_var(1).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(dis(in)).lt.(1.0d-99))     dis(in)= 0.d0
                    if (dabs(dis(in+nnt)).lt.(1.0d-99)) dis(in+nnt)=0.d0
                    write(unit_disp(i),'(3E16.6)') tt1,dis(in),dis(in+nnt)
                enddo
            endif
            ! Velocity
            if (option_out_var(2).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(vel(in)).lt.(1.0d-99))     vel(in)= 0.d0
                    if (dabs(vel(in+nnt)).lt.(1.0d-99)) vel(in+nnt)=0.d0
                    write(unit_vel(i),'(3E16.6)') tt1,vel(in),vel(in+nnt)
                enddo
            endif
            ! Acceleration
            if (option_out_var(3).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(acc(in)).lt.(1.0d-99))     acc(in)= 0.d0
                    if (dabs(acc(in+nnt)).lt.(1.0d-99)) acc(in+nnt)=0.d0
                    write(unit_acc(i),'(3E16.6)') tt1,acc(in),acc(in+nnt)
                enddo
            endif
            ! Compute strain
            if (sum(option_out_var(4:6)).gt.0) then
                duxdx = 0.d0
                duydy = 0.d0
                duxdy = 0.d0
                duydx = 0.d0
                
                if (option_out_var(4).eq.1) then
                    sxx = 0.d0
                    syy = 0.d0 
                    sxy = 0.d0
                    szz = 0.d0  
                endif

                do ie=1,ne
                    im = cs(cs(ie -1) +0)   
                    nn = sdeg_mat(im) +1

                    allocate(ct(nn))
                    allocate(ww(nn))
                    allocate(dd(nn,nn))
                    allocate(dxdx_el(nn))
                    allocate(dxdy_el(nn))
                    allocate(dydx_el(nn))
                    allocate(dydy_el(nn))
                    allocate(ux_el(nn,nn))
                    allocate(uy_el(nn,nn))
                    allocate(duxdx_el(nn,nn))
                    allocate(duxdy_el(nn,nn))
                    allocate(duydx_el(nn,nn))
                    allocate(duydy_el(nn,nn))
                    
                    
                    call LGL(nn,ct,ww,dd)
                    
                    call MAKE_DERIVATIVES_LOC(nn,alfa1(ie),alfa2(ie),beta1(ie),beta2(ie),gamma1(ie),gamma2(ie),ct,&
                        dxdx_el,dxdy_el,dydx_el,dydy_el)
                    
                    do jj=1,nn
                        do ii=1,nn
                            is = nn*(jj -1) + ii
                            in = cs(cs(ie -1) + is)
                            
                            ux_el(ii,jj)=dis(in)
                            uy_el(ii,jj)=dis(in+nnt)
                        enddo
                    enddo

                    call MAKE_STRAIN(nn,dd,dxdx_el,dxdy_el,dydx_el,dydy_el,ux_el,uy_el, &
                        duxdx_el,duxdy_el,duydx_el,duydy_el)
                    
                    do jj=1,nn
                        do ii=1,nn
                            is = nn*(jj -1) + ii
                            in = cs(cs(ie -1) + is)
                            
                            duxdx(in) = duxdx(in) + duxdx_el(ii,jj)
                            duydy(in) = duydy(in) + duydy_el(ii,jj)
                            duxdy(in) = duxdy(in) + duxdy_el(ii,jj)
                            duydx(in) = duydx(in) + duydx_el(ii,jj)
                        enddo
                    enddo
                    

                    if (option_out_var(4).eq.1) then
                        allocate(sxx_el(nn,nn))
                        allocate(syy_el(nn,nn))
                        allocate(szz_el(nn,nn))
                        allocate(sxy_el(nn,nn))
                        allocate(mu_el(nn,nn))
                        allocate(lambda_el(nn,nn))
                        lambda_el = prop_mat(im,2)
                        mu_el     = prop_mat(im,3)
                        
                        call MAKE_STRESS(lambda_el,mu_el,nn,     &
                            duxdx_el,duxdy_el,duydx_el,duydy_el, &
                            sxx_el,syy_el,szz_el,sxy_el)
                        do jj=1,nn
                            do ii=1,nn
                                is = nn*(jj -1) + ii
                                in = cs(cs(ie -1) + is)
                                sxx(in) = sxx(in) + sxx_el(ii,jj)
                                syy(in) = syy(in) + syy_el(ii,jj)
                                szz(in) = szz(in) + szz_el(ii,jj)
                                sxy(in) = sxy(in) + sxy_el(ii,jj)
                            enddo
                        enddo
                        
                    endif
                    if(allocated(ct)) deallocate(ct)
                    if(allocated(ww)) deallocate(ww)
                    if(allocated(dd)) deallocate(dd)
                    if(allocated(dxdx_el)) deallocate(dxdx_el)
                    if(allocated(dxdy_el)) deallocate(dxdy_el)
                    if(allocated(dydx_el)) deallocate(dydx_el)
                    if(allocated(dydy_el)) deallocate(dydy_el)
                    if(allocated(ux_el   )) deallocate(ux_el   )   
                    if(allocated(uy_el   )) deallocate(uy_el   )
                    if(allocated(duxdx_el)) deallocate(duxdx_el)
                    if(allocated(duxdy_el)) deallocate(duxdy_el)
                    if(allocated(duydx_el)) deallocate(duydx_el)
                    if(allocated(duydy_el)) deallocate(duydy_el)
                    if(allocated(sxx_el   )) deallocate(sxx_el   )
                    if(allocated(syy_el   )) deallocate(syy_el   )  
                    if(allocated(szz_el   )) deallocate(szz_el   )   
                    if(allocated(sxy_el   )) deallocate(sxy_el   )
                    if(allocated(mu_el    )) deallocate(mu_el    )
                    if(allocated(lambda_el)) deallocate(lambda_el)

                enddo            
            endif

            if (option_out_var(4) .eq. 1) then

                do in = 1, nnt  
                    sxx(in) = sxx(in) / nodal_counter(in)
                    syy(in) = syy(in) / nodal_counter(in)
                    sxy(in) = sxy(in) / nodal_counter(in)
                    szz(in) = szz(in) / nodal_counter(in)
                enddo
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(sxx(in)).lt.(1.0d-99)) sxx(in)=0.0
                    if (dabs(syy(in)).lt.(1.0d-99)) syy(in)=0.0
                    if (dabs(sxy(in)).lt.(1.0d-99)) sxy(in)=0.0
                    if (dabs(szz(in)).lt.(1.0d-99)) szz(in)=0.0

                    write(unit_stress(i),'(5E16.8)') tt1,sxx(in),syy(in),sxy(in),szz(in)
                enddo
            endif


            if (option_out_var(5) .eq. 1) then  !Out Options Scandella 02.07.2007
                do in = 1, nnt  
                    duxdx(in) = duxdx(in) / nodal_counter(in)
                    duydy(in) = duydy(in) / nodal_counter(in)
                    duxdy(in) = duxdy(in) / nodal_counter(in)
                    duydx(in) = duydx(in) / nodal_counter(in) 
                enddo 

                do i = 1,nmonit 
                    in = node_m(i) 
                    if (dabs(duxdx(in)).lt.(1.0d-99)) duxdx(in)=0.0
                    if (dabs(duydy(in)).lt.(1.0d-99)) duydy(in)=0.0
                    if (dabs(duxdy(in)).lt.(1.0d-99)) duxdy(in)=0.0
                    if (dabs(duydx(in)).lt.(1.0d-99)) duydx(in)=0.0

                    write(unit_strain(i),'(4E16.8)') tt1,duxdx(in),duydy(in),0.5*(duxdy(in)+duydx(in)) 
                enddo
            endif

            if (option_out_var(6) .eq. 1) then  !Out Options Scandella 02.07.2007
                if (option_out_var(5) .ne. 1) then
                    do in = 1, nnt  
                        duxdx(in) = duxdx(in) / nodal_counter(in)
                        duydy(in) = duydy(in) / nodal_counter(in)
                        duxdy(in) = duxdy(in) / nodal_counter(in)
                        duydx(in) = duydx(in) / nodal_counter(in) 
                    enddo 
                endif

                do i = 1,nmonit 
                in = node_m(i) 
                if (dabs(duxdy(in)).lt.(1.0d-99)) duxdy(in)=0.0
                if (dabs(duydx(in)).lt.(1.0d-99)) duydx(in)=0.0

                write(unit_omega(i),'(2E16.8)') tt1, 0.5*(duxdy(in)-duydx(in)) 
                enddo
            endif

            return
        end subroutine WRITE_MONITOR_EL
        
        ! nonlinear output 
        subroutine WRITE_MONITOR_NL(unit_disp,unit_vel,unit_acc,unit_strain,unit_stress,unit_omega,&
            option_out_var,nmonit,ndt_monitor,node_m,nm,ne,nnt,cs,cs_nnz,sdeg_mat,snl,its,tt1,&
            dis,vel,acc,nodal_counter,disout)
            ! 
            implicit none
            ! intent IN
            real*8,                         intent(in)  :: tt1,ndt_monitor
            integer*4,                      intent(in)  :: nm,ne,nnt,cs_nnz,nmonit,its
            integer*4, dimension(6),        intent(in)  :: option_out_var
            integer*4, dimension(nm),       intent(in)  :: sdeg_mat
            integer*4, dimension(nnt),      intent(in)  :: nodal_counter
            integer*4, dimension(nmonit),   intent(in)  :: node_m
            integer*4, dimension(0:cs_nnz), intent(in)  :: cs
            type(nl_element), dimension(ne),intent(in)  :: snl
            ! intent INOUT
            integer*4, dimension(nmonit), intent(inout) :: unit_disp
            integer*4, dimension(nmonit), intent(inout) :: unit_vel
            integer*4, dimension(nmonit), intent(inout) :: unit_acc
            integer*4, dimension(nmonit), intent(inout) :: unit_stress
            integer*4, dimension(nmonit), intent(inout) :: unit_strain
            integer*4, dimension(nmonit), intent(inout) :: unit_omega
            real*8,    dimension(2*nnt) , intent(inout) :: dis,vel,acc
            type(nodepatched),            intent(inout) :: disout
            !
            integer*4                                     ::  in,i
            real*8                                        :: sxx_out,syy_out,szz_out,sxy_out
            real*8                                        :: exx_out,eyy_out,gxy_out
            
            if (option_out_var(1).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(dis(in)).lt.(1.0d-99))     dis(in)= 0.d0
                    if (dabs(dis(in+nnt)).lt.(1.0d-99)) dis(in+nnt)=0.d0
                    write(unit_disp(i),'(3ES16.6)') tt1,dis(in),dis(in+nnt)
                enddo
            endif

            if (option_out_var(2).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(vel(in)).lt.(1.0d-99))     vel(in)= 0.d0
                    if (dabs(vel(in+nnt)).lt.(1.0d-99)) vel(in+nnt)=0.d0
                    write(unit_vel(i),'(3ES16.6)') tt1,vel(in),vel(in+nnt)
                enddo
            endif

            if (option_out_var(3).eq.1) then   
                do i = 1,nmonit
                    in = node_m(i)
                    if (dabs(acc(in)).lt.(1.0d-99))     acc(in)= 0.d0
                    if (dabs(acc(in+nnt)).lt.(1.0d-99)) acc(in+nnt)=0.d0
                    write(unit_acc(i),'(3ES16.6)') tt1,acc(in),acc(in+nnt)
                enddo
            endif
            ! STRESS OUTPUT
            if (option_out_var(4).eq.1) then

                call UPDATE_OUT_STRESS(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
                do i = 1,nmonit
                    in = node_m(i)
                    sxx_out = disout%stress(1,in)/nodal_counter(in)
                    syy_out = disout%stress(2,in)/nodal_counter(in)
                    szz_out = disout%stress(3,in)/nodal_counter(in)
                    sxy_out = disout%stress(4,in)/nodal_counter(in)
                    if (dabs(disout%stress(1,in)).lt.(1.0d-99)) sxx_out=0.d0
                    if (dabs(disout%stress(2,in)).lt.(1.0d-99)) syy_out=0.d0
                    if (dabs(disout%stress(3,in)).lt.(1.0d-99)) szz_out=0.d0
                    if (dabs(disout%stress(4,in)).lt.(1.0d-99)) sxy_out=0.d0
                    write(unit_stress(i),'(5ES16.8)') tt1,sxx_out,syy_out,sxy_out,szz_out
                enddo
            endif
            ! STRAIN OUTPUT
            if (option_out_var(5).eq.1) then

                call UPDATE_OUT_STRAIN(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
                do i = 1,nmonit
                    in = node_m(i) 
                    exx_out = disout%strain(1,in) / nodal_counter(in)
                    eyy_out = disout%strain(2,in) / nodal_counter(in)
                    gxy_out = disout%strain(3,in) / nodal_counter(in)
                    if (dabs(disout%strain(1,in)).lt.(1.0d-99)) exx_out=0.d0
                    if (dabs(disout%strain(2,in)).lt.(1.0d-99)) eyy_out=0.d0
                    if (dabs(disout%strain(3,in)).lt.(1.0d-99)) gxy_out=0.d0
                    write(unit_strain(i),'(4ES16.8)') tt1,exx_out,eyy_out,gxy_out 
                enddo
            endif
            return
        end subroutine WRITE_MONITOR_NL
        !
        subroutine UPDATE_OUT_STRESS(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
            !
            implicit none
            ! intent IN 
            integer*4, intent(in)                       :: ne,nnt,cs_nnz,nm
            integer*4, intent(in), dimension(nm)        :: sdeg_mat
            integer*4, intent(in), dimension(0:cs_nnz)  :: cs
            type(nl_element), intent(in), dimension(ne) :: snl
            ! intent INOUT
            type(nodepatched), intent(inout)            :: disout
            !
            integer*4                                   :: ie,im,nn,in,is,i,j 
            
            do ie = 1,ne
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                do j = 1,nn
                    do i = 1,nn
                        is = nn*(j -1) +i
                        in = cs(cs(ie -1) + is)
                        
                        disout%stress(:,in) = disout%stress(:,in) + &
                            snl(ie)%stress(:,i,j)

                    enddo
                enddo
            enddo
            return
        end subroutine UPDATE_OUT_STRESS
        !
        subroutine UPDATE_OUT_STRAIN(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
            !
            implicit none
            ! intent IN 
            integer*4, intent(in)                       :: ne,nnt,cs_nnz,nm
            integer*4, intent(in), dimension(nm)        :: sdeg_mat
            integer*4, intent(in), dimension(0:cs_nnz)  :: cs
            type(nl_element), intent(in), dimension(ne) :: snl
            ! intent INOUT
            type(nodepatched), intent(inout)            :: disout
            !
            integer*4                                   :: ie,im,nn,in,is,i,j 
            
            do ie = 1,ne
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                do j = 1,nn
                    do i = 1,nn
                        is = nn*(j -1) +i
                        in = cs(cs(ie -1) + is)
                        
                        disout%strain(:,in) = disout%strain(:,in) + &
                            snl(ie)%strain(:,i,j)

                    enddo
                enddo
            enddo
            return
        end subroutine UPDATE_OUT_STRAIN
        !
        subroutine UPDATE_OUT_PLASTIC_STRAIN(ne,nnt,cs_nnz,cs,nm,sdeg_mat,snl,disout)
            !
            implicit none
            ! intent IN 
            integer*4, intent(in)                       :: ne,nnt,cs_nnz,nm
            integer*4, intent(in), dimension(nm)        :: sdeg_mat
            integer*4, intent(in), dimension(0:cs_nnz)  :: cs
            type(nl_element), intent(in), dimension(ne) :: snl
            ! intent INOUT
            type(nodepatched), intent(inout)            :: disout
            !
            integer*4                                   :: ie,im,nn,in,is,i,j 
            
            do ie = 1,ne
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                do j = 1,nn
                    do i = 1,nn
                        is = nn*(j -1) +i
                        in = cs(cs(ie -1) + is)
                        
                        disout%plastic_strain(:,in) = disout%plastic_strain(:,in) + &
                            snl(ie)%plastic_strain(:,i,j)

                    enddo
                enddo
            enddo
            return
    end subroutine UPDATE_OUT_PLASTIC_STRAIN
end module write_output
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

