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

!> @brief Coonvert morse matrix to CRS matrix 
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] nnod dimension of the matrix
!> @param[in] nnz_dg nnzero elements of the matrix
!> @param[in] IDG_MORSE,JDG_MORSE,MDG_MORSE matrix in morse format
!> @param[out] IDG,JDG,MDG matrix in CRS format

         subroutine MORSE_2_RCS(nnod, nnz_dg, IDG_MORSE, JDG_MORSE, MDG_MORSE, &
                         IDG, JDG, MDG)
                         
         implicit none
         integer*4 :: nnod, nnz_dg
         integer*4, dimension(0:2*nnod) :: IDG
         integer*4, dimension(nnz_dg) :: IDG_MORSE, JDG_MORSE, JDG
         
         real*8, dimension(nnz_dg)  :: MDG, MDG_MORSE
         
         integer*4 :: i,j,k
         
         
         IDG = 0;
         k = 1;
         do i = 1, 2*nnod
            do j = 1, nnz_dg
               
               if(IDG_MORSE(j) .eq. i) then
                 IDG(i) = IDG(i) + 1
                 JDG(k) = JDG_MORSE(j)
                 MDG(k) = MDG_MORSE(j)

                 k = k + 1
               endif
            enddo
         enddo                 
                 
         do i = 1, 2*nnod 
            IDG(i) = IDG(i) + IDG(i-1)
         enddo    
                         
         end subroutine MORSE_2_RCS
