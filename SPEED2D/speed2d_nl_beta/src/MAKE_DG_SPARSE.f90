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

!> @brief Build DG matrix in sparse morse format
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] nnz_dg_total nnzero elements in MDG_TOTAL
!> @param[in] IDG_TOTAL, JDG_TOTAL, MDG_TOTAL  old dg matrix 
!> @param[in] nnz_dg nnzero elements in MDG
!> @param[out] IDG, JDG, MDG  new dg matrix

        subroutine MAKE_DG_SPARSE(nnz_dg_total,IDG_TOTAL,JDG_TOTAL,MDG_TOTAL,nnz_dg,IDG,JDG,MDG)
        
        
        implicit none
        
        integer*4 :: nnz_dg_total
        integer*4 :: nnz_dg
        integer*4, dimension(nnz_dg_total) :: IDG_TOTAL, JDG_TOTAL
        real*8, dimension(nnz_dg_total) :: MDG_TOTAL
       
        integer*4, dimension(nnz_dg) :: IDG, JDG
        real*8, dimension(nnz_dg) :: MDG
       
        integer*4 :: i, j 
        
        
        j = 1
        do i = 1, nnz_dg_total
           if (IDG_TOTAL(i) .ne. 0) then
              IDG(j) = IDG_TOTAL(i);
              JDG(j) = JDG_TOTAL(i);
              MDG(j) = MDG_TOTAL(i);
              j = j + 1;
           endif
        enddo      
        
        end subroutine MAKE_DG_SPARSE
