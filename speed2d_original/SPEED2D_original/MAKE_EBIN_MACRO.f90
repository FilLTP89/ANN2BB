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

!> @brief Makes pointer for connectivity of mesh nodes. 
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
!> @param[in] nnode  grid nodes  
!> @param[in] nelem  number of elements (quad)
!> @param[in] con   con(i,1) -> element 
!!                  con(i,2),...,con(i,5) -> grid nodes of element con(i,1) 
!> @param[in] bw  bw(i) contains the multiplicity of the grid node i
!!                  (e.g. if the node node i is shared by 4 elements bw(i) = 4  
!> @param[in] nnz  number of grid nodes + number of grid nodes including repetitions + 1
!!
!> @param[out] bin  pointer for mesh nodes connectivity 
!  
!    E.G. two elements grid   
!    1------2------3
!    |      |      |       nnz = 15
!    | el.1 | el.2 |       bw(1) = bw(3) = bw(4) = bw(6) = 1
!    4------5------6       bw(2) = bw(5) = 2
!      
!    BEFORE
!    bin(0) = 7                      i4count(1) = 7     
!    bin(1) = 8                      i4count(2) = 8 
!    bin(2) = 10                     i4count(3) = 10
!    bin(3) = 11                     i4count(4) = 11
!    bin(4) = 12                     i4count(5) = 12
!    bin(5) = 14                     i4count(6) = 14
!    bin(6) = 15                     i4count(7) = 15
!    bin(7) = ... = bin(15) = 0
!
!
!    AFTER
!    bin(0),...,bin(6) unchanged
!    bin(7) = 1                      i4count(1) = 8 
!    node 1 el. 1    
!    bin(8) = 1                      i4count(2) = 10 
!    bin(9) = 2                      i4count(3) = 11
!    node 2 el. 1 & 2
!    bin(10) = 2                     i4count(4) = 12
!    node 3 el. 2
!    bin(11) = 1                     i4count(5) = 14
!    node 4 el. 1    
!    bin(12) = 1                     i4count(6) = 15
!    bin(13) = 2                     i4count(7) = 16
!    node 5 el. 1 & 2
!    bin(14) = 2
!    node 6 el. 2        
!    bin(15) = 2        
!    node 7 el. 2
!
!**************************************************************************************************  

      subroutine MAKE_EBIN_MACRO(nnode,nelem,con,bw,nnz,bin)
      
      implicit none
      
      integer*4 :: nnode,nelem,nnz
      integer*4, dimension(nelem,5) :: con
      integer*4, dimension(nnode) :: bw
      integer*4, dimension(0:nnz) :: bin
      
      integer*4, dimension(:), allocatable :: ic
      integer*4 :: i,j,k,in,ie,nn2
      integer*4 :: bj
      
      
      allocate(ic(nnode))
      
      bin(0) = nnode +1
      do in = 1,nnode
         bin(in) = bin(in-1) + bw(in)
      enddo
      
      do in = 1,nnode
         ic(in) = bin(in-1)
      enddo
      
      do ie = 1,nelem
         do i = 1,4
            in = con(ie,i +1)
            
            bin(ic(in)) = ie
            ic(in) = ic(in) +1
         enddo
      enddo
      
      deallocate(ic)
      
      return
      
      end subroutine MAKE_EBIN_MACRO
