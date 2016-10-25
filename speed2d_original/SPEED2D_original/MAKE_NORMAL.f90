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

!> @brief Compute the normal vector
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] ind label for the edge
!> @param[in] xs1,ys1 coordinate for the 1st vertex 
!> @param[in] xs2,ys2 coordinate for the 2nd vertex
!> @param[in] nx,ny component of the normal vector

    subroutine MAKE_NORMAL(ind, xs1, xs2, ys1, ys2, nx, ny)

        implicit none
        
        integer*4 :: ind
        real*8 :: edge_lx, edge_ly, edge_ll 

        real*8, intent(in)  :: xs1, xs2, ys1, ys2                                   
        real*8, intent(out) :: nx,ny

        edge_lx = xs2 - xs1
        edge_ly = ys2 - ys1
        edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)
        nx = edge_ly / edge_ll
        ny = -1.0d0 * edge_lx / edge_ll                  
                        
                                    
end subroutine MAKE_NORMAL
