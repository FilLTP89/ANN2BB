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


!> @brief Find index position of a number in a vector 
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
 
!> @param[in] Jsparse vector containing column indeces
!> @param[in] dim1 dimension of Jsparse
!> @param[in] ind_s starting index for the search
!> @param[in] ind_e ending index for the search
!> @param[in] ind_f index to be found
!> @param[out] ind_j position of the index in the matrix

subroutine FIND_POSITION(Jsparse, dim1, ind_s, ind_e, ind_f, ind_j)

      
      implicit none

      integer*4, intent(in) :: ind_s, ind_e, ind_f
      integer*4, intent(out) :: ind_j      
      integer*4 :: par, dim1
      integer*4, dimension(dim1) :: Jsparse
      
      
      par = ind_s
      do while (par .le. ind_e)
         if(Jsparse(par) .eq. ind_f) then
            ind_j = par
            return
         endif

         par = par + 1
      enddo

     return
end subroutine FIND_POSITION
