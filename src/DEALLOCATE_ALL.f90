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

subroutine DEALLOCATE_ALL(ne,u1,u2,vel,acc,fk,fe,fd,snl,disout,update_index_el_az,nodal_counter)  
    ! 
    use fields 
    !
    implicit none
    ! intent IN
    integer*4, intent(in)                                      :: ne
    ! intent INOUT
    real*8,     dimension(:), allocatable,  intent(inout)      :: u1,u2,vel,acc,fk,fe,fd
    integer*4,  dimension(:), allocatable,  intent(inout)      :: update_index_el_az,nodal_counter
    type(nl_element), dimension(:), allocatable, intent(inout) :: snl
    type(nodepatched), intent(inout)                           :: disout
    ! counters
    integer*4                                                  :: ie 
    
    ! deallocate element-wise variables to be stored
    do ie = 1,ne
        ! allocation
        deallocate(snl(ie)%lambda)
        deallocate(snl(ie)%mu)
        deallocate(snl(ie)%syld)
        deallocate(snl(ie)%ckin)
        deallocate(snl(ie)%kkin)
        deallocate(snl(ie)%rinf)
        deallocate(snl(ie)%biso)
        !
        deallocate(snl(ie)%radius)
        deallocate(snl(ie)%stress)
        deallocate(snl(ie)%strain)
        deallocate(snl(ie)%center)
        deallocate(snl(ie)%plastic_strain)

    enddo
    deallocate(snl)
    
    ! current displacement/velocity/acceleration vectors
    deallocate(u1)  
    deallocate(u2)  
    deallocate(vel) 
    deallocate(acc)
    ! current internal & external forces 
    deallocate(fk)
    deallocate(fd)
    deallocate(fe)
    ! update index 
    deallocate(update_index_el_az)
    ! nodal counter (for results' average)
    if (allocated(nodal_counter)) deallocate(nodal_counter) 
    if (allocated(disout%stress)) deallocate(disout%stress)
    if (allocated(disout%strain)) deallocate(disout%strain)
    if (allocated(disout%plastic_strain)) deallocate(disout%plastic_strain)
    return
end subroutine DEALLOCATE_ALL
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

