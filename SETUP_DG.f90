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

!> @brief Setup for DG faces. Stores data on arrays faces and area_nodes. 
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] nm  number of materials
!> @param[in] sd  polynomial degree vector
!> @param[in] tag_mat labels for materials
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local spectral connectivity vector
!> @param[in] nn_loc number of nodes
!> @param[in] ne_loc number of elements
!> @param[in] xs x-coord of spectral nodes
!> @param[in] ys y-coord of spectral nodes
!> @param[in] nel_dg_loc  number of  dg elements
!> @param[in] i4count vector identifying nodes where DG conditions are applied
!> @param[in] alfa1 costant for bilinear mapping
!> @param[in] alfa2 costant for bilinear mapping
!> @param[in] beta1 costant for bilinear mapping
!> @param[in] beta2 costant for bilinear mapping
!> @param[in] gamma1 costant for bilinear mapping
!> @param[in] gamma2 costant for bilinear mapping
!> @param[in] delta1 costant for bilinear mapping
!> @param[in] delta2 costant for bilinear mapping
!> @param[out] faces  material, element, face for a DG face
!> @param[out] area_nodes info about DG faces.
!!      area_nodes(1,i) = area of the face, 
!!      area_nodes(2,i),...,area_nodes(9,i) -> costants for the bilinear map 

      subroutine SETUP_DG(nm, sd, tag_mat, cs_nnz_loc, cs_loc, &
                           nn_loc, ne_loc, xs,ys,&
                           nel_dg_loc, i4count,&
                           alfa1, alfa2, beta1, beta2, gamma1, gamma2, delta1, delta2, & 
                           faces, area_nodes)


     implicit none

     integer*4 :: nm, cs_nnz_loc, nn_loc, ne_loc, nel_dg_loc
     integer*4 :: im, nn, ie, ned, err_out
     integer*4 :: ne1, ne2, ne3, ne4, ic1, ic2, ic3, ic4
     integer*4 :: ne5, ne6, ne7, ne8, ic5, ic6, ic7, ic8

     integer*4, dimension(nm) :: tag_mat, sd
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(nn_loc) :: i4count

     integer*4, dimension(3,nel_dg_loc), intent(inout) :: faces

     real*8 :: surf

     real*8, dimension(nn_loc) :: xs,ys
     real*8, dimension(ne_loc) :: alfa1,alfa2
     real*8, dimension(ne_loc) :: beta1,beta2
     real*8, dimension(ne_loc) :: gamma1,gamma2, delta1, delta2
     
     real*8, dimension(9,nel_dg_loc), intent(inout) :: area_nodes
     real*8 :: x1,x2,y1,y2

!*****************************************************************************************
!loading of data structure faces containing the description for DG faces
!*****************************************************************************************      

      nel_dg_loc = 0   
      ned = cs_loc(0) - 1
      
      do im = 1,nm
      
         nn = sd(im) +1         

         do ie = 1,ned
            if (cs_loc(cs_loc(ie -1) + 0) .eq. tag_mat(im)) then
               
               !1st edge
               ne1 = cs_loc(cs_loc(ie -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn)
               
               x1 = xs(ne1); x2 = xs(ne2);
               y1 = ys(ne1); y2 = ys(ne2);
                       
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) ) then
                  nel_dg_loc = nel_dg_loc +1
                  
                  surf = dsqrt((x1-x2)**2 + (y1-y2)**2);
                  
                  !1st edge
                                              
                  faces(1,nel_dg_loc) = tag_mat(im)
                  faces(2,nel_dg_loc) = ie
                  faces(3,nel_dg_loc) = 1

                  area_nodes(1,nel_dg_loc) = surf
                  area_nodes(2,nel_dg_loc) = alfa1(ie)
                  area_nodes(3,nel_dg_loc) = alfa2(ie)                
                  area_nodes(4,nel_dg_loc) = beta1(ie)
                  area_nodes(5,nel_dg_loc) = beta2(ie)                 
                  area_nodes(6,nel_dg_loc) = gamma1(ie)
                  area_nodes(7,nel_dg_loc) = gamma2(ie)
                  area_nodes(8,nel_dg_loc) = delta1(ie)
                  area_nodes(9,nel_dg_loc) = delta2(ie)

                     
               endif               

               !2nd edge 
               ne1 = cs_loc(cs_loc(ie -1) +nn)             
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) ) then
                  nel_dg_loc = nel_dg_loc +1
                  
                  surf = dsqrt((x1-x2)**2 + (y1-y2)**2);
                  
                  !2n edge
                                              
                  faces(1,nel_dg_loc) = tag_mat(im)
                  faces(2,nel_dg_loc) = ie
                  faces(3,nel_dg_loc) = 2

                  area_nodes(1,nel_dg_loc) = surf
                  area_nodes(2,nel_dg_loc) = alfa1(ie)
                  area_nodes(3,nel_dg_loc) = alfa2(ie)                
                  area_nodes(4,nel_dg_loc) = beta1(ie)
                  area_nodes(5,nel_dg_loc) = beta2(ie)                 
                  area_nodes(6,nel_dg_loc) = gamma1(ie)
                  area_nodes(7,nel_dg_loc) = gamma2(ie)
                  area_nodes(8,nel_dg_loc) = delta1(ie)
                  area_nodes(9,nel_dg_loc) = delta2(ie)
                     
               endif               

               !3rd edge
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn)
               ne2 = cs_loc(cs_loc(ie -1) +nn*(nn -1) +1)

               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) ) then
                  nel_dg_loc = nel_dg_loc +1
                  
                  surf = dsqrt((x1-x2)**2 + (y1-y2)**2);
                  
                  !3rd edge
                                              
                  faces(1,nel_dg_loc) = tag_mat(im)
                  faces(2,nel_dg_loc) = ie
                  faces(3,nel_dg_loc) = 3

                  area_nodes(1,nel_dg_loc) = surf
                  area_nodes(2,nel_dg_loc) = alfa1(ie)
                  area_nodes(3,nel_dg_loc) = alfa2(ie)                
                  area_nodes(4,nel_dg_loc) = beta1(ie)
                  area_nodes(5,nel_dg_loc) = beta2(ie)                 
                  area_nodes(6,nel_dg_loc) = gamma1(ie)
                  area_nodes(7,nel_dg_loc) = gamma2(ie)
                  area_nodes(8,nel_dg_loc) = delta1(ie)
                  area_nodes(9,nel_dg_loc) = delta2(ie)
                     
               endif               

               !4th edge
               ne1 = cs_loc(cs_loc(ie -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn*(nn -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) ) then
                  nel_dg_loc = nel_dg_loc +1
                  
                  surf = dsqrt((x1-x2)**2 + (y1-y2)**2);
                  
                  !4th edge
                                              
                  faces(1,nel_dg_loc) = tag_mat(im)
                  faces(2,nel_dg_loc) = ie
                  faces(3,nel_dg_loc) = 4

                  area_nodes(1,nel_dg_loc) = surf
                  area_nodes(2,nel_dg_loc) = alfa1(ie)
                  area_nodes(3,nel_dg_loc) = alfa2(ie)                
                  area_nodes(4,nel_dg_loc) = beta1(ie)
                  area_nodes(5,nel_dg_loc) = beta2(ie)                 
                  area_nodes(6,nel_dg_loc) = gamma1(ie)
                  area_nodes(7,nel_dg_loc) = gamma2(ie)
                  area_nodes(8,nel_dg_loc) = delta1(ie)
                  area_nodes(9,nel_dg_loc) = delta2(ie)
                     
               endif               
      
            endif
         enddo
      enddo
      


      return

      end subroutine SETUP_DG
