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
!    along with SPEED.  If not, see <http://wqw.gnu.org/licenses/>.

!> @brief Allocate/Initialize all variables necessary for step by step calculations 
!! @author Filippo Gatti
!> @date February,2016
!> @version 1.0
!> @param[inout] sxx,syy,sxy,szz nodal values for the stress tensor (on element)

subroutine ALLOINIT_NL_ALL(ne,sdeg_mat,nm,nnt,cs_nnz,cs,u1,u2,vel,acc,v1,fk,fe,fd,&
    duxdx,duxdy,duydx,duydy,dduxdx,dduxdy,dduydx,dduydy,sxx,syy,szz,sxy,dsxx,dsyy,dszz,dsxy,&
    xkin_all,dxkin_all,riso_all,driso_all,epl_all,depl_all,option_out_var,nl_sism,sism,&
    update_index_el_az,nodal_counter)  

    implicit none
    ! intent IN
    integer*4,                              intent(in)      :: ne,nm,nnt,cs_nnz,nl_sism
    integer*4,  dimension(6),               intent(in)      :: option_out_var
    integer*4,  dimension(nm),              intent(in)      :: sdeg_mat
    integer*4,  dimension(0:cs_nnz),        intent(in)      :: cs
    real*8,     dimension(2*nnt),           intent(in)      :: v1
    ! intent INOUT
    real*8,     dimension(:), allocatable,  intent(inout)   :: u1,u2,vel,acc,fk,fe,fd,sism
    real*8,     dimension(:), allocatable,  intent(inout)   :: xkin_all,epl_all,riso_all
    real*8,     dimension(:), allocatable,  intent(inout)   :: dxkin_all,depl_all,driso_all
    real*8,     dimension(:), allocatable,  intent(inout)   :: duxdx,duydy,duydx,duxdy
    real*8,     dimension(:), allocatable,  intent(inout)   :: dduxdx,dduydy,dduydx,dduxdy
    real*8,     dimension(:), allocatable,  intent(inout)   :: sxx,syy,szz,sxy,dsxx,dsyy,dszz,dsxy
    integer*4,  dimension(:), allocatable,  intent(inout)   :: update_index_el_az,nodal_counter
    ! counters
    integer*4                                               :: nn,im,iaz,ie,in,is,i,j
    ! current stress state (to be saved for nonlinear calculations & outputs)
    allocate(sxx(nnt));   sxx = 0.d0
    allocate(syy(nnt));   syy = 0.d0 
    allocate(sxy(nnt));   sxy = 0.d0
    allocate(szz(nnt));   szz = 0.d0
    ! time-step stress increment (to be stored to update current stress state)
    allocate(dsxx(nnt)); dsxx = 0.d0
    allocate(dsyy(nnt)); dsyy = 0.d0 
    allocate(dsxy(nnt)); dsxy = 0.d0
    allocate(dszz(nnt)); dszz = 0.d0  
    ! time-step strain rate state (to be saved for nonlinear calculations) 
    allocate(duxdx(nnt)); duxdx = 0.d0
    allocate(duydy(nnt)); duydy = 0.d0
    allocate(duxdy(nnt)); duxdy = 0.d0
    allocate(duydx(nnt)); duydx = 0.d0
    ! time-step strain increment (to be stored to update current stress state)
    allocate(dduxdx(nnt)); dduxdx = 0.d0
    allocate(dduydy(nnt)); dduydy = 0.d0 
    allocate(dduxdy(nnt)); dduxdy = 0.d0
    allocate(dduydx(nnt)); dduydx = 0.d0  
    ! current displacement/velocity/acceleration vectors
    allocate(u1(2*nnt));  u1 = 0.d0
    allocate(u2(2*nnt));  u2 = 0.d0
    allocate(vel(2*nnt)); vel = v1 
    allocate(acc(2*nnt)); acc = 0.d0
    ! current non linear hardening variables & increments (to be saved for nonlinear calculations)
    allocate(xkin_all(4*nnt));   xkin_all = 0.0d0
    allocate(dxkin_all(4*nnt)); dxkin_all = 0.0d0
    allocate(epl_all(4*nnt));     epl_all = 0.0d0
    allocate(depl_all(4*nnt));   depl_all = 0.0d0
    allocate(riso_all(nnt));     riso_all = 0.0d0
    allocate(driso_all(nnt));   driso_all = 0.0d0
    ! current internal & external forces 
    allocate(fk(2*nnt)); fk = 0.0d0
    allocate(fd(2*nnt)); fd = 0.0d0
    allocate(fe(2*nnt)); fe = 0.0d0
    
    ! external tensor moment load
    if(nl_sism.gt.0) then
        allocate(sism(2*nnt))
    endif
    ! update index 
    allocate(update_index_el_az(2*nnt))
    do iaz = 1,2*nnt
        update_index_el_az(iaz) = iaz
    enddo
    ! nodal counter (for results' average)
    if(sum(option_out_var(4:6)).ge.1) then 
        allocate(nodal_counter(nnt)) 
        nodal_counter = 0
        do ie = 1,ne
            im = cs(cs(ie-1)+0)
            nn = sdeg_mat(im)+1
            do j = 1,nn
                do i = 1,nn 
                    is = nn*(j -1) +i
                    in = cs(cs(ie -1) + is)
                    nodal_counter(in) = nodal_counter(in)+1
                enddo 
            enddo
        enddo
    endif  
    return
end subroutine ALLOINIT_NL_ALL
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
