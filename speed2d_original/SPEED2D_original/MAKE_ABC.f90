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

!> @brief Set the number of edges where ABC is applied (nedge_abc).
!! Set the edge-id where ABC is applied (iedge_abc).
!! Set the elements to which the edges belong to (ielem_abc).
 
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
 
!> @param[in] nedge number of boundary edges
!> @param[in] cs_nnz_bc legnth of cs_bc
!> @param[in] cs_bc vector for boundary connectivity
!> @param[in] cs_nnz length of cs 
!> @param[in] cs vector for grid connectivity
!> @param[in] nl_abc number of abc conditions
!> @param[in] tag_abc label for abc conditions
!> @param[inout] i4count work vector for counting the number of ab edges  
!> @param[out] nedge_abc number of ab edges

 
     subroutine MAKE_ABC(nedge_abc, nedge, i4count, cs_nnz_bc,cs_bc,&
                         cs_nnz,cs,nl_abc,tag_abc)
     
     implicit none
     integer*4 :: nedge_abc,cs_nnz_bc,cs_nnz,nl_abc,nelem_abc
     integer*4 :: nedge,i,iedge, ied1,ied2,ie,ne,nn
     integer*4 :: iel1,iel2,iel3,iel4
     integer*4, dimension(nl_abc) :: tag_abc
     integer*4, dimension(0:cs_nnz_bc) :: cs_bc
     integer*4, dimension(0:cs_nnz) :: cs
     
     integer*4, dimension(nedge) :: i4count
     
      
     ne = cs(0)-1;
     
     do i = 1, nl_abc
        do iedge = 1, nedge
  
           if (cs_bc(cs_bc(iedge -1) +0) .eq. tag_abc(i)) then 
               ied1 = cs_bc(cs_bc(iedge -1) +1)
               ied2 = cs_bc(cs_bc(iedge) -1)
                    
               do ie = 1, ne
                   nn = cs_bc(iedge) - cs_bc(iedge -1) -1 
                   iel1 = cs(cs(ie -1) +1)
                   iel2 = cs(cs(ie -1) +nn)
                   iel3 = cs(cs(ie -1) +nn*nn)
                   iel4 = cs(cs(ie -1) +nn*(nn -1) +1)
                   if (((ied1.eq.iel1).or.(ied1.eq.iel2).or. &
                        (ied1.eq.iel3).or.(ied1.eq.iel4)).and. &
                        ((ied2.eq.iel1).or.(ied2.eq.iel2).or. &
                        (ied2.eq.iel3).or.(ied2.eq.iel4))) then
                        i4count(iedge) = iedge
                   endif
               enddo
           endif
         enddo
      enddo
         
   
         
      do iedge = 1,nedge
         if (i4count(iedge).gt.0) then
            nedge_abc = nedge_abc +1
            i4count(iedge) = nedge_abc
         endif
      enddo
      
      end subroutine MAKE_ABC
