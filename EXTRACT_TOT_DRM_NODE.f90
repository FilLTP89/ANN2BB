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

!> @brief Computes Node's coordinates of total DRM nodes and writing PDRM in *. mat file
!! @author 
!> @date August, 2015
!> @version 1.0

!> @param[in] nnode number of nodes
!> @param[in] nnode_MDRM number of DRM element nodes without boundary DRM nodes
!> @param[in] node_MDRM DRM element nodes without boundary DRM nodes
!> @param[in] nnode_BD size of boundary DRM nodes without duplicates
!> @param[in] node_BD list of boundary DRM nodes without duplicates
!> @param[in] xs x coordinate
!> @param[in] ys  y coordinate
!> @param[in] nnode_EL number of DRM elements nodes without DRM boundary nodes and duplicates
!> @param[in] nnode_TOT number of DRM total nodes
!> @param[in] node_TOT nodes of total DRM domain
!> @param[out] xx_TOT node's coordinates of total DRM domain
!> @param[out] yy_TOT node's coordinates of total DRM domain


    subroutine EXTRACT_TOT_DRM_NODE(nnode,nnode_MDRM,node_MDRM, &
	           nnode_BD,node_BD,xs,ys,nnode_EL,nnode_TOT,node_TOT,&
			   xx_TOT,yy_TOT)   
			   
			    
!   © POLIMI, 2006, All Rights Reserved
!   Author: Laura Scandella 
	
	implicit none

	integer*4 :: nnode,nnode_MDRM,nnode_BD,nnode_EL,nnode_TOT
    real*8, dimension(nnode) :: xs,ys
	integer*4, dimension(nnode_MDRM) :: node_MDRM
	integer*4, dimension(nnode_BD) :: node_BD
	integer*4, dimension(nnode_TOT) :: node_TOT
	real*8, dimension(nnode_TOT) :: xx_TOT,yy_TOT

	integer*4 :: i,j,pnod,k,check

    node_TOT(1:nnode_BD) = node_BD(1:nnode_BD)

    k = 0
	do i = 1,nnode_MDRM
	   pnod = node_MDRM(i)
	   if (k.ge.1) then
	      check = 0
	         do j = nnode_BD+1,nnode_BD+k
	            if (pnod.eq.node_TOT(j)) then
		           check = 1
				 exit
		        endif
	         enddo
	         if (check.ne.1) then
                k = k+1
	            node_TOT(nnode_BD+k)=pnod	   	              
             endif
		else
           k = k+1
	       node_TOT(nnode_BD+k)=pnod
		endif
	enddo

	      
!    Node's coordinates of DRM internal boundary (xx_BDRM, yy_BDRM)  
	 do i = 1,nnode
	    do j = 1,nnode_TOT
		    if (i.eq.node_TOT(j)) then
			   xx_TOT(j) = xs(i)
			   yy_TOT(j) = ys(i)
			   exit
			endif
	    enddo
     enddo


    return

    end subroutine EXTRACT_TOT_DRM_NODE