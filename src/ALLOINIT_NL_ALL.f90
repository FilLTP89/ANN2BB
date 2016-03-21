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
    duxdx,duxdy,duydx,duydy,sxx,syy,szz,sxy,dsxx,dsyy,dszz,dsxy,&
    xkin_all,riso_all,epl_all,option_out_var,nl_sism,sism,update_index_el_az,nodal_counter)  
    
    implicit none
    integer*4,                              intent(in)      :: ne,nm,nnt,cs_nnz
    integer*4,  dimension(6),               intent(in)      :: option_out_var
    integer*4,  dimension(nm),              intent(in)      :: sdeg_mat
    integer*4,  dimension(0:cs_nnz),        intent(in)      :: cs
    real*8,     dimension(2*nnt),           intent(in)      :: v1
    real*8,     dimension(:), allocatable,  intent(inout)   :: u1,u2,vel,acc,fk,fe,fd
    real*8,     dimension(:), allocatable,  intent(inout)   :: sism
    real*8,     dimension(:), allocatable,  intent(inout)   :: xkin_all,epl_all,riso_all
    real*8,     dimension(:), allocatable,  intent(inout)   :: duxdx,duydy,duydx,duxdy
    real*8,     dimension(:), allocatable,  intent(inout)   :: sxx,syy,szz,sxy
    real*8,     dimension(:), allocatable,  intent(inout)   :: dsxx,dsyy,dszz,dsxy
    integer*4,  dimension(:), allocatable,  intent(inout)   :: update_index_el_az
    integer*4,  dimension(:), allocatable,  intent(inout)   :: nodal_counter
    integer*4, intent(in)                                   :: nl_sism
    integer*4                                               :: nn,im,iaz,ie,in,is,i,j

    ! global stress state to be saved for nonlinear calculations
    allocate(sxx(nnt)); sxx = 0.d0
    allocate(syy(nnt)); syy = 0.d0 
    allocate(sxy(nnt)); sxy = 0.d0
    allocate(szz(nnt)); szz = 0.d0  
    allocate(dsxx(nnt)); dsxx = 0.d0
    allocate(dsyy(nnt)); dsyy = 0.d0 
    allocate(dsxy(nnt)); dsxy = 0.d0
    allocate(dszz(nnt)); dszz = 0.d0  
    ! global strain state to be saved for nonlinear calculations (verify)
    allocate(duxdx(nnt)); duxdx = 0.d0
    allocate(duydy(nnt)); duydy = 0.d0
    allocate(duxdy(nnt)); duxdy = 0.d0
    allocate(duydx(nnt)); duydx = 0.d0
    ! global displacement/velocity/acceleration vector
    allocate(u1(2*nnt)); u1 = 0.d0
    allocate(u2(2*nnt)); u2 = 0.d0
    allocate(vel(2*nnt)); vel = v1 
    allocate(acc(2*nnt)); acc = 0.d0
    ! global non linear variables
    allocate(xkin_all(4*nnt)); xkin_all = 0.0d0
    allocate(epl_all(4*nnt));  epl_all  = 0.0d0
    allocate(riso_all(nnt));   riso_all = 0.0d0
    ! global forces
    allocate(fk(2*nnt)); fk = 0.0d0
    allocate(fd(2*nnt)); fd = 0.0d0
    allocate(fe(2*nnt)); fe = 0.0d0
    
    if(nl_sism.gt.0) then
        allocate(sism(2*nnt))
    endif
    
    allocate(update_index_el_az(2*nnt))
    do iaz = 1,2*nnt
        update_index_el_az(iaz) = iaz
    enddo
    
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
