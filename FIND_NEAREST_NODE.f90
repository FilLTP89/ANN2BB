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


!> @brief Find nearest node to monitored point xt,yt
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
 
!> @param[in] n number of nodes
!> @param[in] xs,ys coordinates of spectral nodes
!> @param[in] xt,yt monitored point
!> @param[out] nt id of the nearest node 

      subroutine FIND_NEAREST_NODE(n,xs,ys,xt,yt,nt)
      
      implicit none
      
      integer*4 :: n,nt
      real*8, dimension(n) :: xs,ys
      real*8 :: xt,yt
      
      integer*4 :: i
      real*8 :: dx,dy,d2,d2min
      
      d2min = 1.0d30
      
      do i = 1,n
         dx = xs(i) - xt
         dy = ys(i) - yt
         d2 = dx*dx + dy*dy
         if (d2.lt.d2min) then
            d2min = d2
            nt = i
         endif
      enddo
      
      return
      end subroutine FIND_NEAREST_NODE
      

