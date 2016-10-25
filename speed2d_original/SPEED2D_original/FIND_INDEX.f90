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


!> @brief Find next position where storing a column index
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
 
!> @param[in] M matrix
!> @param[in] dim1 number of rows for M
!> @param[in] dim2 number of columns for M
!> @param[in] ind_i row index
!> @param[out] ind_j column index
 
     subroutine FIND_INDEX(M,dim1,dim2,ind_i, ind_j)

      implicit none

      integer, intent(in) :: ind_i
      integer, intent(out) :: ind_j      
      integer*4 :: dim1, dim2
      integer*4, dimension(dim1, dim2) :: M
      
      
      ind_j = 1

      if(M(ind_i,ind_j) .eq. 0 ) then
          return
      endif 
      
      do while (M(ind_i,ind_j) .ne. 0 )
          ind_j = ind_j + 1           
      enddo  

     return
     end subroutine FIND_INDEX

