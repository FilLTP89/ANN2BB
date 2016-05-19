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


module seismic
    use fields
    !    
    implicit none
    !
    contains
        
        subroutine MAKE_SEISMIC_FORCES(nnt,nm,ne,nf,cs_nnz,cs,sdeg_mat,nfunc_data,nl_sism,length_cns,&
            check_node_sism,check_dist_node_sism,func_data,func_type,tag_func,func_indx,facsmom,tt1,&
            alfa1,alfa2,beta1,beta2,gamma1,gamma2,displ,mvec,sism,fe)
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
            real*8,           intent(inout), dimension(2*nnt)   :: sism,fe
            ! 
            real*8,     dimension(:),  allocatable      :: ct,ww
            real*8,     dimension(:),  allocatable      :: dxdx,dydy,dxdy,dydx
            real*8,     dimension(:,:),allocatable      :: dd,det_j,fxs,fys
            real*8,     dimension(:,:),allocatable      :: sxxs,syys,szzs,sxys
            real*8,     dimension(:,:,:), allocatable   :: dstrain,dstrial
            real*8                                      :: t1ux,t1uy,t2ux
            real*8                                      :: t2uy,t1fx,t1fy,t2fx,t2fy
            integer*4                                   :: ie,ip,iq,il,im,nn,is,in
            !
            fe = 0.d0
            do ie = 1,ne 
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                
                ! ALLOCATION/INITIALIZATION LOCAL VARIABLES 
                call ALLOINIT_LOC_EL(ie,nnt,cs,cs_nnz,nn,ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain, &
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
                        call MAKE_INTERNAL_FORCE(nn,ww,dd,dxdx,dxdy,dydx,dydy,&
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
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
