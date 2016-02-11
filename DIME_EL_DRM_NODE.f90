!    Copyright (C) 2014 The SPEED FOUNDATION
!    Author: 
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

!> @brief Computes number of DRM elements nodes without DRM boundary nodes and duplicates
!! @author 
!> @date August, 2015
!> @version 1.0

!> @param[in] nnode_MDRM number of DRM element nodes without boundary DRM nodes
!> @param[in] node_MDRM DRM element nodes without boundary DRM nodes
!> @param[out] nnode_EL number of DRM elements nodes without DRM boundary nodes and duplicates


    subroutine DIME_EL_DRM_NODE(nnode_MDRM,node_MDRM,nnode_EL)

!   © POLIMI, 2006, All Rights Reserved
!   Author: Laura Scandella      
	
	implicit none

	integer*4 :: nnode_MDRM
	integer*4, dimension(nnode_MDRM) :: node_MDRM

	integer*4 :: nnode_EL

	integer*4 :: i,j,check

    nnode_EL = 0
	do i = 1,nnode_MDRM
	   check = 0
	   do j = i,nnode_MDRM
	      if (node_MDRM(i).eq.node_MDRM(j)) then
		     check = check+1 
		  endif
       enddo
	   if (check.eq.1) then
	   nnode_EL = nnode_EL+1
	   endif
	enddo

	return

    end subroutine 	DIME_EL_DRM_NODE