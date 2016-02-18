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

!> @brief Find element index from verteces number.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] Ennz number of element
!> @param[in] Ebin vector containing spectral node repetition
!> @param[in] n1 node index of first vertex of the element (line) 
!> @param[in] n2 node index of second vertex of the element (line)
!> @param[out] ie element index

      subroutine GET_EDGE_ELEMENT(Ennz,Ebin,n1,n2,ie)
      
      implicit none
      
      integer*4 :: Ennz,n1,n2,ie
      integer*4, dimension(0:Ennz) :: Ebin
      
      integer*4 :: i,j
      
      ie = 0
      do i = Ebin(n1 -1),Ebin(n1) -1
         do j = Ebin(n2 -1),Ebin(n2) -1
            if (Ebin(i).eq.Ebin(j)) then
               ie = Ebin(i)
            endif
         enddo
      enddo
      
      return
      
      end subroutine GET_EDGE_ELEMENT

