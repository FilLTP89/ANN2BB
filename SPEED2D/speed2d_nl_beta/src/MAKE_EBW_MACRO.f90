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

!> @brief Computes multeplicity for mesh nodes. 
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
!> @param[in] nnode  grid nodes
!> @param[in] nelem  number of elements (hexes)
!> @param[in] con con(i,1) -> element,  
!!                con(i,2),...,con(i,5) -> grid nodes of element con_matrix(i,1) 
!> @param[out] bw  bw(i) matrix containing the multiplicity of the i-th grid node 
!!                  (e.g. if the node i is shared by 4 elements bw(i) = 4  
!> @param[out] nnz number of grid nodes + number of grid nodes including repetitions + 1

      subroutine MAKE_EBW_MACRO(nnode,nelem,con,bw,nnz)
      
      implicit none
      
      integer*4 :: nnode,nelem,nnz
      integer*4, dimension(nelem,5) :: con
      integer*4, dimension(nnode) :: bw
      
      integer*4 :: i,in,ie,nn
      
      nnz = nnode +1
      
      do in = 1,nnode
         bw(in) = 0
      enddo
      
      do ie = 1,nelem
         do i = 1,4
            in = con(ie,i +1)
            
            bw(in) = bw(in) +1
            nnz = nnz +1
         enddo
      enddo
      
      return
      
      end subroutine MAKE_EBW_MACRO

