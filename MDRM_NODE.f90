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

!> @brief Computes DRM elements nodes without DRM boundary nodes
!! @author 
!> @date August, 2015
!> @version 1.0

!> @param[in] ii_node_DRM number of nodes of DRM domain
!> @param[in] node_DRM node's numbers of DRM domain
!> @param[in] nnode_BD size of boundary DRM nodes without duplicates
!> @param[in] node_BD list of boundary DRM nodes without duplicates
!> @param[in] nnode_MDRM number of DRM element nodes without boundary DRM nodes
!> @param[out] node_MDRM DRM element nodes without boundary DRM nodes


    subroutine MDRM_NODE(ii_node_DRM,node_DRM,nnode_BD,& 
	                           node_BD,nnode_MDRM,node_MDRM)

!   © POLIMI, 2006, All Rights Reserved
!   Author: Laura Scandella     
	
	implicit none

	integer*4 :: ii_node_DRM,nnode_BD,nnode_MDRM
	integer*4, dimension(ii_node_DRM) :: node_DRM
    integer*4, dimension(nnode_BD) :: node_BD
    integer*4, dimension(nnode_MDRM) :: node_MDRM

    integer*4 :: nnode_EL,nnod_TOT
	integer*4 :: i,j,pnod,check,k

!   DRM element nodes without boundary DRM nodes (node_MDRM)
    k = 0
    do j = 1,ii_node_DRM
	   pnod = node_DRM(j) 
	   check = 0
	   do i = 1,nnode_BD
	      if (pnod.eq.node_BD(i)) then
             check = 1
			 exit
		  endif 
	   enddo
	   if (check.ne.1) then
	       k = k+1
		   node_MDRM(k) = pnod
	   endif
	enddo

	return

    end subroutine MDRM_NODE