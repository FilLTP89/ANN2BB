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

!> @brief Find the position for storing the column index for the i-th row
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] M matrix
!> @param[in] dim1 number of rows for M
!> @param[in] dim2 number of columns for M
!> @param[in] ind_i row id
!> @param[out] ind_j position where storing the next column id

     subroutine COUNT_NNZ_EL(M,dim1,dim2,ind_i,ind_j)

      
      implicit none

      integer, intent(out) :: ind_j      
      integer, intent(in) :: ind_i
      integer*4 :: dim1, dim2
      integer*4, dimension(dim1, dim2) :: M

      ind_j = 0

      do while (M(ind_i,ind_j+1) .ne. 0 )
          ind_j = ind_j + 1           
      enddo  
     
      return

     end subroutine COUNT_NNZ_EL
