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

subroutine ALLOINIT_NL_ALL(ne,sdeg_mat,nm,nnt,cs_nnz,cs,prop_mat,u1,u2,vel,acc,v1,fk,fe,fd,&
    snl,option_out_var,nl_sism,sism,update_index_el_az,nodal_counter)  
    
    use nonlinear2d

    implicit none
    
    ! intent IN
    integer*4,                              intent(in)      :: ne,nm,nnt,cs_nnz,nl_sism
    real*8,     dimension(nm,9)             intent(in)      :: prop_mat
    integer*4,  dimension(6),               intent(in)      :: option_out_var
    integer*4,  dimension(nm),              intent(in)      :: sdeg_mat
    integer*4,  dimension(0:cs_nnz),        intent(in)      :: cs
    real*8,     dimension(2*nnt),           intent(in)      :: v1
    ! intent INOUT
    real*8,     dimension(:), allocatable,  intent(inout)   :: u1,u2,vel,acc,fk,fe,fd,sism
    integer*4,  dimension(:), allocatable,  intent(inout)   :: update_index_el_az,nodal_counter
    type(nl_element), dimension(:), allocatable, intent(inout)   :: snl
    ! counters
    integer*4                                               :: nn,im,iaz,ie,in,is,i,j
    
    ! allocate element-wise variables to be stored
    allocate(snl(ne))
    do ie = 1,ne
        im = cs(cs(ie-1)+0)
        nn = sdeg_mat(im)+1 
        ! allocation
        allocate(snl(ie)%lambda(nn,nn))
        allocate(snl(ie)%mu(nn,nn))
        allocate(snl(ie)%syld(nn,nn))
        allocate(snl(ie)%ckin(nn,nn))
        allocate(snl(ie)%kkin(nn,nn))
        allocate(snl(ie)%rinf(nn,nn))
        allocate(snl(ie)%biso(nn,nn))
        !
        allocate(snl(ie)%radius(nn,nn))
        allocate(snl(ie)%stress(4,nn,nn))
        allocate(snl(ie)%strain(3,nn,nn))
        allocate(snl(ie)%center(4,nn,nn))
        allocate(snl(ie)%plastic_strain(3,nn,nn))
        ! initialization
        snl(ie)%lambda = prop_mat(im,2)
        snl(ie)%mu     = prop_mat(im,3)
        snl(ie)%syld   = prop_mat(im,5)
        snl(ie)%ckin   = prop_mat(im,6)
        snl(ie)%kkin   = prop_mat(im,7)
        snl(ie)%rinf   = prop_mat(im,8)
        snl(ie)%biso   = prop_mat(im,9)
        !
        snl(ie)%radius(:,:)             = 0.d0
        snl(ie)%stress(:,:,:)           = 0.d0
        snl(ie)%strain(:,:,:)           = 0.d0
        snl(ie)%center(:,:,:)           = 0.d0
        snl(ie)%plastic_strain(:,:,:)   = 0.d0

    enddo
    
    ! current displacement/velocity/acceleration vectors
    allocate(u1(2*nnt));  u1 = 0.d0
    allocate(u2(2*nnt));  u2 = 0.d0
    allocate(vel(2*nnt)); vel = v1 
    allocate(acc(2*nnt)); acc = 0.d0
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
