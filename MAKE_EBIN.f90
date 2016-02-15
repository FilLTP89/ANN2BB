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

!> @brief Makes a pointer for spectral connectivity vector.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
!> @param[in] nnode  GLL nodes
!> @param[in] cs(0:cs_nnz) spectral connectivity vector
!> @param[in] bw bw(i) contains the multiplicity of the GLL node i
!!                  (e.g. if the GLL node i is shared by 4 elements bw(i) = 4  
!> @param[in] nnz  number of GLL nodes + number of GLL nodes including repetitions + 1
!> @param[out] bin as explained in MAKE_EBIN_MACRO.f90


      subroutine MAKE_EBIN(nnode,cs_nnz,cs,bw,nnz,bin)
      
      implicit none
      
      integer*4 :: nnode,cs_nnz,nnz
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nnode) :: bw
      integer*4, dimension(0:nnz) :: bin
      
      integer*4, dimension(:), allocatable :: ic
      integer*4 :: i,j,k,in,ie,ne,nn2
      integer*4 :: bj
      
      
      allocate(ic(nnode))
      
      bin(0) = nnode +1
      do in = 1,nnode
         bin(in) = bin(in -1) + bw(in)
      enddo
      
      do in = 1,nnode
         ic(in) = bin(in -1)
      enddo
      
      ne = cs(0) -1
      do ie = 1,ne
         nn2 = cs(ie) - cs(ie -1) -1
         do i = 1,nn2
            in = cs(cs(ie -1) +i)
            
            bin(ic(in)) = ie
            ic(in) = ic(in) +1
         enddo
      enddo
      
      deallocate(ic)
      
      return
      
      end subroutine MAKE_EBIN

