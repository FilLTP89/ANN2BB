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

!> @brief Computes DRM boundary nodes without duplicates
!! @author 
!> @date August, 2015
!> @version 1.0

!> @param[in] ii_node_BDRM number of nodes of DRM boundary domain
!> @param[in] node_BDRM node's numbers of DRM internal boundary
!> @param[in] nnode_BD size of boundary DRM nodes without duplicates
!> @param[out] node_BD list of boundary DRM nodes without duplicates


      subroutine EXTRACT_BDRM_NODE(node_BDRM,ii_node_BDRM, &
	                              nnode_BD,node_BD)
	                      	                						 
!    © POLIMI, 2006, All Rights Reserved
!    Author: Laura Scandella

     implicit none
      
	  integer*4 :: ii_node_BDRM
      integer*4, dimension (ii_node_BDRM) :: node_BDRM 

      integer*4 :: nnode_BD
      integer*4, dimension (nnode_BD) :: node_BD

	  integer*4 :: i,j,k   


!    List of boundary DRM nodes without duplicates
     k = 0
	 do i = 1,ii_node_BDRM-1
		if (node_BDRM(i+1).ne.node_BDRM(i)) then
		k = k+1
		node_BD(k) = node_BDRM(i)
		endif
	 enddo
	 k = k+1
	 node_BD(k) = node_BDRM(ii_node_BDRM)  

     return

     end subroutine EXTRACT_BDRM_NODE    