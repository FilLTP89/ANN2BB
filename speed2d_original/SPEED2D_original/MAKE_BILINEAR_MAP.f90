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

!> @brief Assign values from nodes array. 
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[out] nodes costant values for the bilinear map to be read
!> @param[out] c_alfa1 costant values for the bilinear map
!> @param[out] c_alfa2 costant values for the bilinear map
!> @param[out] c_beta1 costant values for the bilinear map
!> @param[out] c_beta2 costant values for the bilinear map
!> @param[out] c_gamma1 costant values for the bilinear map
!> @param[out] c_gamma2 costant values for the bilinear map
!> @param[out] c_delta1 costant values for the bilinear map
!> @param[out] c_delta2 costant values for the bilinear map

     subroutine MAKE_BILINEAR_MAP(nodes, c_alfa1, c_alfa2, c_beta1, c_beta2, &
                                         c_gamma1, c_gamma2, &
                                         c_delta1, c_delta2)

     implicit none
     
     real*8, intent(out) :: c_alfa1, c_alfa2
     real*8, intent(out) :: c_beta1, c_beta2
     real*8, intent(out) :: c_gamma1, c_gamma2, c_delta1, c_delta2

     real*8, dimension(8) :: nodes
     

     c_alfa1 = nodes(1)
     c_alfa2 = nodes(2)

     c_beta1 = nodes(3)
     c_beta2 = nodes(4)

     c_gamma1 = nodes(5)
     c_gamma2 = nodes(6)

     c_delta1 = nodes(7)
     c_delta2 = nodes(8)

 
    end subroutine MAKE_BILINEAR_MAP
