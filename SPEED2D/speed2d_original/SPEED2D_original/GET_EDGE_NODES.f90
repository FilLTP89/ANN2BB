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
!    along with SPEED.  If not, see <http://www.gnu.org/licenses/>.

!> @brief Computes total number of boundary nodes 
!! and second derivative on a given point x.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] nnode number of nodes
!> @param[in] nnz length of cs
!> @param[in] cs  spectral connectivity vector
!> @param[in] nl_bc  load number for boundary condition
!> @param[in] tag_bc  label for boundary condition 
!> @param[in] node_index  node_index(i) = 1 if the node belongs to the boundary
!> @param[out] nedge_nodes number of nodes lying on the boundary

      subroutine GET_EDGE_NODES(nnode,nnz,cs,nl_bc,tag_bc,&
                                nedge_nodes,node_index)
      
      
      implicit none
      
      integer*4 :: nnode,nnz,nl_bc
      integer*4, dimension(0:nnz) :: cs
      integer*4, dimension(nl_bc) :: tag_bc
      integer*4, dimension(nnode) :: node_index
      integer*4 :: nedge_nodes
      
      integer*4 :: i,j,ie,ne,nn,check
      
      do i = 1,nnode
         node_index(i) = 0
      enddo
      
      nedge_nodes = 0
      
      if (nnz.gt.0) then
         ne = cs(0) -1
         do ie = 1,ne
            nn = cs(ie) - cs(ie -1) -1
            
            check = 0
            do j = 1,nl_bc
               if (cs(cs(ie -1) +0).eq.tag_bc(j)) check = 1
            enddo
            
            if (check.ne.0) then
               do i = 1,nn
                  node_index(cs(cs(ie -1) +i)) = 1
               enddo
            endif
         enddo
         
         do i = 1,nnode
            if (node_index(i).ne.0) then
               nedge_nodes = nedge_nodes +1
               node_index(i) = nedge_nodes
            endif
         enddo
      endif
      
      return
      
      end subroutine GET_EDGE_NODES
