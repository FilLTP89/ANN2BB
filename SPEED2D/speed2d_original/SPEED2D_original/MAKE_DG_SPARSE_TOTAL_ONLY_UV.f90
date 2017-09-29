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

!> @brief Build DG matrix in sparse morse format from local DG matrices
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
!> @param[in] nelem_dg number of dg elements
!> @param[in] el_new structure containing dg elements
!> @param[in] nnz_dg_total nnzero elements in MDG_TOTAL
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[in] nnod number of nodes
!> @param[in] nmat number of materials
!> @param[in] sdeg_mat polynomial degrees
!> @param[out] IDG_TOTAL, JDG_TOTAL, MDG_TOTAL  dg matrix in sparse morse format

     subroutine MAKE_DG_SPARSE_TOTAL_ONLY_UV(nelem_dg, el_new, nnz_dg_total, cs_nnz, cs, &
                                     nnod, nmat, sdeg_mat, IDG_TOTAL, &
                                     JDG_TOTAL, MDG_TOTAL)


     use max_var
     use DGJUMP

     implicit none
  

     type(el4loop), dimension(nelem_dg), intent(inout):: el_new      

     integer*4 :: nnz_dg_total, nelem_dg, cs_nnz, start, finish,nmat
     integer*4, dimension(0:cs_nnz) :: cs
     integer*4, dimension(nmat) :: sdeg_mat
 
     integer*4, dimension(nnz_dg_total), intent(inout) :: IDG_TOTAL, JDG_TOTAL
     real*8, dimension(nnz_dg_total), intent(inout) :: MDG_TOTAL
                      
     integer*4 :: i,j,k, i_total, len_plus, len_minus, mm, ie, ie_c, nnod
     integer*4, dimension(:), allocatable :: IDG_PLUS, IDG_MINUS, nodes_minus, nodes_plus
     
     i_total = 1
   
     do i = 1, nelem_dg
        len_plus = 0;
        len_minus = 0;
        
        ie = el_new(i)%ind
        
        len_plus = 2*el_new(i)%deg**2
        do j = 1, el_new(i)%num_of_ne
           mm = sdeg_mat(el_new(i)%el_conf(j,0)) + 1
           len_minus = len_minus + 2*mm**2
        enddo
        
        allocate(nodes_plus(len_plus), nodes_minus(len_minus))
        nodes_plus(1:len_plus/2) = cs(cs(ie-1)+1: cs(ie)-1);
        nodes_plus(len_plus/2+1:len_plus) = cs(cs(ie-1)+1: cs(ie)-1) + nnod
                
        start = 1
        finish = 0;
        do j = 1, el_new(i)%num_of_ne
           ie_c = el_new(i)%el_conf(j,1)
           mm = sdeg_mat(el_new(i)%el_conf(j,0)) + 1                      
           finish = finish + mm**2
           
           nodes_minus(start:finish) = cs(cs(ie_c-1)+1: cs(ie_c)-1);

           start = start + mm**2;
           finish = finish + mm**2;
           
           nodes_minus(start:finish) = cs(cs(ie_c-1)+1: cs(ie_c)-1) + nnod
           start = start + mm**2;
              
     
        enddo
          
                  
        allocate(IDG_PLUS(el_new(i)%nnz_plus_only_uv),IDG_MINUS(el_new(i)%nnz_minus_only_uv))
        
        do j = 1, 2*(el_new(i)%deg)**2 
          
           do k = el_new(i)%IPLUS_only_uv(j-1) + 1, el_new(i)%IPLUS_only_uv(j)
               IDG_PLUS(k) = j
           enddo
           do k = el_new(i)%IMIN_only_uv(j-1) + 1, el_new(i)%IMIN_only_uv(j)
               IDG_MINUS(k) = j
           enddo
        enddo
        
        do j = 1,  el_new(i)%nnz_plus_only_uv
           IDG_TOTAL(i_total) = nodes_plus(IDG_PLUS(j))
           JDG_TOTAL(i_total) = nodes_plus(el_new(i)%JPLUS_only_uv(j))
           MDG_TOTAL(i_total) = el_new(i)%matPlus_only_uv(j)
           i_total = i_total + 1
        enddo   
        
        do j = 1,  el_new(i)%nnz_minus_only_uv
           IDG_TOTAL(i_total) = nodes_plus(IDG_MINUS(j))
           JDG_TOTAL(i_total) = nodes_minus(el_new(i)%JMIN_only_uv(j))
           MDG_TOTAL(i_total) = el_new(i)%matMin_only_uv(j)
           i_total = i_total + 1
        enddo   


        deallocate(nodes_plus, nodes_minus, IDG_PLUS, IDG_MINUS)
        
    enddo       
     
     

 
     end subroutine MAKE_DG_SPARSE_TOTAL_ONLY_UV
