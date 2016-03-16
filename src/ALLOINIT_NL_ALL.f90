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

subroutine ALLOINIT_NL_ALL(ne,sdeg_mat,nm,nnt,cs_nnz,cs,u1,u2,fk,fe,fd,sism,vel,&
        acc,v1,update_index_el_az,duxdx,duxdy,duydx,duydy,sxx,syy,szz,sxy,&
        xkin_all,riso_all,epl_all,option_out_var,nodal_counter)  
    implicit none
    integer*4,                              intent(in)  :: ne,nm,nnt,cs_nnz
    integer*4,  dimension(6),               intent(in)  :: option_out_var
    integer*4,  dimension(nm),              intent(in)  :: sdeg_mat
    integer*4,  dimension(0:cs_nnz),        intent(in)  :: cs
    real*8,     dimension(2*nnt),           intent(in)  :: v1
    integer*4,  dimension(:),  allocatable, intent(out) :: update_index_el_az
    integer*4,  dimension(:),  allocatable, intent(out) :: nodal_counter
    real*8,     dimension(:),  allocatable, intent(out) :: u1,u2,fk,fe,fd,sism
    real*8,     dimension(:),  allocatable, intent(out) :: xkin_all,epl_all,riso_all
    real*8,     dimension(:),  allocatable, intent(out) :: duxdx,duydy,duydx,duxdy
    real*8,     dimension(:),  allocatable, intent(out) :: sxx,syy,szz,sxy,vel,acc
    integer*4                                           :: nn,im,iaz,ie,in,is,i,j

    allocate(u1(2*nnt))
    allocate(u2(2*nnt))
    allocate(fk(2*nnt))
    allocate(fe(2*nnt))
    allocate(fd(2*nnt))
    allocate(sism(2*nnt))
    allocate(vel(2*nnt))
    allocate(acc(2*nnt))
    allocate(update_index_el_az(2*nnt))
    allocate(riso_all(nnt))
    allocate(xkin_all(4*nnt))
    allocate(epl_all(4*nnt))

    do iaz = 1,2*nnt
        update_index_el_az(iaz) = iaz
    enddo
    
    ! global stress state to be saved for nonlinear calculations
    allocate(sxx(nnt))
    allocate(syy(nnt))
    allocate(sxy(nnt))
    allocate(szz(nnt))           
    sxx = 0.d0
    syy = 0.d0 
    sxy = 0.d0
    szz = 0.d0  
    ! global strain state to be saved for nonlinear calculations (verify)
    allocate(duxdx(nnt))
    allocate(duydy(nnt))
    allocate(duxdy(nnt))
    allocate(duydx(nnt)) 
    duxdx = 0.d0
    duydy = 0.d0
    duxdy = 0.d0
    duydx = 0.d0
    u1          = 0.d0
    u2          = 0.d0
    vel         = v1
    acc         = 0.0d0
    xkin_all    = 0.0d0
    epl_all     = 0.0d0
    riso_all    = 0.0d0

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
