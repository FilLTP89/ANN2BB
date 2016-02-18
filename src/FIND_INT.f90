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


!> @brief Find element in a matrix
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
 
!> @param[in] M matrix
!> @param[in] dim1 number of rows for M
!> @param[in] dim2 number of columns for M
!> @param[in] ind_i row index
!> @param[in] ind_j column index
!> @param[out] tt 1/0 element (ind_i,ind_j) is found/not found

     subroutine FIND_INT(M,dim1,dim2,ind_i,ind_j, tt)

      implicit none

      integer, intent(in) :: ind_i, ind_j
      integer, intent(out) :: tt      
      integer :: ij
      integer*4 :: dim1, dim2
      integer*4, dimension(dim1, dim2) :: M

      tt = 0
      if(M(ind_i,1) .eq. 0 ) then
          tt = 0
          return
      endif 
      
      ij = 1
      do while (ij .le. dim2 )
          if(M(ind_i,ij) .eq. ind_j) then
            tt = 1
            return
          endif
          ij = ij + 1           
      enddo  

     return
     end subroutine FIND_INT
