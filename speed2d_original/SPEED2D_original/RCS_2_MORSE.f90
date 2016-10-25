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

!> @brief Coonvert CRS matrix to morse matrix 
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] NNZ_K nnzero elements for the CRS matrix
!> @param[in] K_TOT,IK_TOT,JK_TOT CRS matrix
!> @param[in] NNZ_K_MORSE nnzero elements for the morse matrix
!> @param[out] IK_MORSE,JK_MORSE,K_MORSE morse format

      subroutine RCS_2_MORSE(nnod, K_TOT, IK_TOT, JK_TOT, NNZ_K, &
			IK_MORSE, JK_MORSE, K_MORSE, NNZ_K_MORSE )
			
      implicit none		
	
      integer*4 :: nnod, NNZ_K, i,j,k, NNZ_K_MORSE
      integer*4, dimension(0:2*nnod) :: IK_TOT
      integer*4, dimension(NNZ_K) :: JK_TOT
      integer*4, dimension(NNZ_K_MORSE) :: IK_MORSE, JK_MORSE				
      real*8, dimension(NNZ_K) :: K_TOT
      real*8, dimension(NNZ_K_MORSE) :: K_MORSE

      
      k = 1
      do i = 1, 2*nnod
         do j = IK_TOT(i-1), IK_TOT(i)-1
            if (K_TOT(j) .ne. 0.d0) then
               IK_MORSE(k) = i;
               JK_MORSE(k) = JK_TOT(j);
               K_MORSE(k) = K_TOT(j);
               k = k + 1
            endif   
         enddo
      enddo   
      
			
       end subroutine RCS_2_MORSE			
