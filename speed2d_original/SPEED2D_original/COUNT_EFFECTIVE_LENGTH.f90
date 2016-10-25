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

!> @brief Compute effective non-zero elements for matrix MDG (without duplicates)
!! @author Ilario Mazzieri
!> @date April, 2014 
!> @version Beta
!> @param[in] nnz_dg_total length for vectros IDG, JDG, MDG
!> @param[in] IDG,JDG,MDG sparse matrix in morse format
!> @param[out] nnz_dg effective length for matrix MDG

        subroutine COUNT_EFFECTIVE_LENGTH(nnz_dg_total,IDG,JDG,MDG,nnz_dg)
        
        implicit none
        
        integer*4 :: nnz_dg_total
        integer*4, intent(out) :: nnz_dg
        integer*4, dimension(nnz_dg_total) :: IDG, JDG
        real*8, dimension(nnz_dg_total) :: MDG
        integer*4 :: i, j 
        
        i = 1;
        do while (i .lt. nnz_dg_total)
           do  j = i + 1, nnz_dg_total
               if (IDG(j) .eq. IDG(i) .and. JDG(j) .eq. JDG(i) .and. IDG(j) .ne. 0) then 
                  MDG(i) = MDG(i) + MDG(j)
                  IDG(j) = 0
                  JDG(j) = 0
                endif
           enddo
          i = i + 1
        enddo
        
        nnz_dg = 0
        do i = 1, nnz_dg_total        
           if (IDG(i) .ne. 0) nnz_dg = nnz_dg + 1
        enddo    
        
        end subroutine COUNT_EFFECTIVE_LENGTH
