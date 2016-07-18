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

!> @brief Counts number of DG elements (local and global)
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] nm = number of materials
!> @param[in] sd(nm) = spectral degree vector
!> @param[in] tag_mat(nm) = tag for materials
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc spectral connectivity vector
!> @param[in] nn_loc number of  nodes 
!> @param[in] xs x-coord local spectral nodes 
!> @param[in] ys y-coord local spectral nodes 
!> @param[in] i4count  vector identifying nodes where DG conditions are applied
!> @param[out] nel_dg  number of dg elements

     subroutine GET_DIME_DG(nm, sd, tag_mat, cs_nnz_loc, cs_loc, &
                           nn_loc,xs,ys,nel_dg,i4count)

     implicit none


     integer*4 :: nm, cs_nnz_loc, nn_loc, nel_dg
     integer*4 :: im, nn, ie, ned
     integer*4 :: ne1, ne2

     integer*4, dimension(nm) :: tag_mat, sd
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(nn_loc) :: i4count
     
     real*8, dimension(nn_loc) :: xs,ys

     
     nel_dg = 0
     ned = cs_loc(0) - 1

      do im = 1,nm
      
         nn = sd(im) +1         

         do ie = 1,ned
            if (cs_loc(cs_loc(ie -1) + 0) .eq. tag_mat(im)) then
               !1-2
               ne1 = cs_loc(cs_loc(ie -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn)
               
                                           
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0)) then
                  nel_dg = nel_dg +1
               endif
               
               !2-3
               ne1 = cs_loc(cs_loc(ie -1) +nn)             
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn)
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0)) then
                  nel_dg = nel_dg +1
               endif
               
               !3-4
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn)
               ne2 = cs_loc(cs_loc(ie -1) +nn*(nn -1) +1)
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) ) then
                  nel_dg = nel_dg +1
               endif
               
               !4-1
               ne1 = cs_loc(cs_loc(ie -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*(nn -1) +1)
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0)) then
                  nel_dg = nel_dg +1
               endif
               
            endif
         enddo
      enddo    
      
      

      end subroutine GET_DIME_DG

