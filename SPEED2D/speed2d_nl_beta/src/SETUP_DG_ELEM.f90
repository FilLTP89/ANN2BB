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

!> @brief Setup for DG data structure.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] nm  number of materials
!> @param[in] sd spectral degree vector 
!> @param[in] tag_mat label for materials
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local spectral connectivity
!> @param[in] nn_loc number of nodes
!> @param[in] ne_loc number of elements
!> @param[in] xs x-coordinate of GLL nodes
!> @param[in] ys y-coordinate of GLL nodes
!> @param[in] nel_dg_loc number of DG elements
!> @param[in] i4count  vector containing info for DG interface (1 in the i-th position if the i-th 
!!                node lies on a DG interface)
!> @param[in] alfa1 bilinear mapping constant
!> @param[in] alfa2 bilinear mapping constant
!> @param[in] beta1 bilinear mapping constant
!> @param[in] beta2 bilinear mapping constant
!> @param[in] gamma1 bilinear mapping constant
!> @param[in] gamma2 bilinear mapping constant
!> @param[in] delta1 bilinear mapping constant
!> @param[in] delta2 bilinear mapping constant
!> @param[in] tag_dg_el label for DG interfaces
!> @param[in] tag_dg_yn label for projecting nodes from a DG surface
!> @param[in] nload_dg number of DG interfaces
!> @param[in] con_bc connectivity matrix for boundary faces
!> @param[in] nface number of boundary faces
!> @param[out] dg_els  data structure for dg interface elements 
!> @param[out] scratch_dg_els  temporary data structure for dg interface elements 

     subroutine SETUP_DG_ELEM(nm, sd, tag_mat, cs_nnz_loc, cs_loc, &
                           nn_loc, ne_loc,&
                           xs,ys, nel_dg_loc, i4count, &
                           alfa1, alfa2, beta1, beta2, &
                           gamma1, gamma2, &
                           delta1, delta2, &
                           dg_els, scratch_dg_els, &
                           tag_dg_el, tag_dg_yn, nload_dg, &
                           con_bc, nface)
                     

     use max_var
     use str_mesh 
     use str_mesh_scratch                   


     implicit none
                   
     
     type(ELEMENT), dimension(nel_dg_loc), intent(inout) :: dg_els
     type(scratch_ELEMENT), dimension(nel_dg_loc), intent(inout) :: scratch_dg_els
  
     character*70 :: filempi, filename
     
     integer*4 :: nm, cs_nnz_loc, nn_loc, ne_loc, nel_dg_loc, nel_dg_glo
     integer*4 :: nload_dg, tag_ind, nface
     integer*4 :: im, nn, ie, ned
     integer*4 :: ne1, ne2
     integer*4 :: el_conf, face_conf, face_found, imate, iele, iface
     integer*4 :: ip, k, j, i, it, ih, ik, tt, indice, node_not_ass

     integer*4, dimension(nm) :: tag_mat, sd
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(nn_loc) :: i4count

     integer*4, dimension(nload_dg) :: tag_dg_el, tag_dg_yn

     integer*4, dimension(:,:), allocatable :: ielem_dg 
     integer*4, dimension(:,:), allocatable :: mat_el_face
     integer*4, dimension(nface,5) :: con_bc

     real*8 :: normal_x, normal_y, normal_z

     real*8, dimension(:), allocatable :: ctgl,wwgl
     real*8, dimension(:), allocatable :: ct, ww
     real*8, dimension(nn_loc) :: xs,ys
     real*8, dimension(ne_loc) :: alfa1,alfa2
     real*8, dimension(ne_loc) :: beta1,beta2
     real*8, dimension(ne_loc) :: gamma1,gamma2, delta1, delta2
     
     real*8, dimension(:,:), allocatable :: dd
     real*8, dimension(:,:), allocatable :: normalxyz

      nel_dg_loc = 0     
      ned = cs_loc(0) - 1


      do im = 1,nm

         nn = sd(im) +1         

         do ie = 1,ned
            if (cs_loc(cs_loc(ie -1) + 0) .eq. tag_mat(im)) then
               
               !1st edge
               ne1 = cs_loc(cs_loc(ie -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +nn)
                       
                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0)) then

                  call GET_TAG_BC(con_bc, nface, ne1, ne2, tag_dg_el, nload_dg, tag_ind)
                  
                  if (tag_ind == 0) print*, 'ERROR IN GET_TAG_BC'
                   
                  nel_dg_loc = nel_dg_loc + 1                    
                  dg_els(nel_dg_loc)%ind_el = ie
                  dg_els(nel_dg_loc)%face_el = 1
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = nofqp
 
                  dg_els(nel_dg_loc)%proj_yn = tag_dg_yn(tag_ind)
 
                 
                  call MAKE_NORMAL(xs(ne1), xs(ne2), ys(ne1), ys(ne2), normal_x, normal_y)

                  dg_els(nel_dg_loc)%nx = normal_x
                  dg_els(nel_dg_loc)%ny = normal_y
                              
                              
                  allocate(ct(dg_els(nel_dg_loc)%quad_rule), ww(dg_els(nel_dg_loc)%quad_rule), & 
                           dd(dg_els(nel_dg_loc)%quad_rule, dg_els(nel_dg_loc)%quad_rule))
                           
                  call MAKE_LGL_NW(dg_els(nel_dg_loc)%quad_rule, ct, ww, dd)  
                  
                  allocate(ctgl(dg_els(nel_dg_loc)%quad_rule), wwgl(dg_els(nel_dg_loc)%quad_rule))
                  
                  call MAKE_GL_NW(dg_els(nel_dg_loc)%quad_rule, ctgl, wwgl)  
                  
                  ip = 0
                  do j = 1,1
                     do i = 1,dg_els(nel_dg_loc)%quad_rule

                        ip = ip + 1
                        scratch_dg_els(nel_dg_loc)%x_nq(ip) = alfa1(ie)*ctgl(i) + beta1(ie)*ct(j) &
                           + gamma1(ie)*ctgl(i)*ct(j) + delta1(ie)
                        
                        scratch_dg_els(nel_dg_loc)%y_nq(ip) = alfa2(ie)*ctgl(i) + beta2(ie)*ct(j) &
                           + gamma2(ie)*ctgl(i)*ct(j) + delta2(ie)
                        
                      
                        dg_els(nel_dg_loc)%wx_pl(ip) = 1.d0
                        dg_els(nel_dg_loc)%wy_pl(ip) = wwgl(i)

                      enddo
                   enddo
                   deallocate(ct,ww,dd)
                   deallocate(ctgl,wwgl)

               endif               
               
               !2nd edge
               ne1 = cs_loc(cs_loc(ie -1) +nn)             
               ne2 = cs_loc(cs_loc(ie -1) +nn*nn)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0)) then
 
                  call GET_TAG_BC(con_bc, nface, ne1, ne2, tag_dg_el, nload_dg, tag_ind)
                  
                  if (tag_ind == 0) print*, 'ERROR IN GET_TAG_BC'

                  nel_dg_loc = nel_dg_loc +1
                  dg_els(nel_dg_loc)%ind_el = ie
                  dg_els(nel_dg_loc)%face_el = 2
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = nofqp
 
                  dg_els(nel_dg_loc)%proj_yn = tag_dg_yn(tag_ind)
                 
                  call MAKE_NORMAL(xs(ne1), xs(ne2), ys(ne1), ys(ne2), normal_x, normal_y)

                  dg_els(nel_dg_loc)%nx = normal_x
                  dg_els(nel_dg_loc)%ny = normal_y
                              
                  allocate(ct(dg_els(nel_dg_loc)%quad_rule), ww(dg_els(nel_dg_loc)%quad_rule), &
                           dd(dg_els(nel_dg_loc)%quad_rule,dg_els(nel_dg_loc)%quad_rule))
                           
                  call MAKE_LGL_NW(dg_els(nel_dg_loc)%quad_rule,ct,ww,dd)  
                  
                  allocate(ctgl(dg_els(nel_dg_loc)%quad_rule),wwgl(dg_els(nel_dg_loc)%quad_rule))
                  
                  call MAKE_GL_NW(dg_els(nel_dg_loc)%quad_rule,ctgl,wwgl)  
                  
                  ip = 0
                    do j = 1, dg_els(nel_dg_loc)%quad_rule
                      do i = dg_els(nel_dg_loc)%quad_rule, dg_els(nel_dg_loc)%quad_rule 
                        ip = ip + 1

                        scratch_dg_els(nel_dg_loc)%x_nq(ip) = alfa1(ie)*ct(i) + beta1(ie)*ctgl(j) &
                           + gamma1(ie)*ct(i)*ctgl(j) + delta1(ie)
                        
                        scratch_dg_els(nel_dg_loc)%y_nq(ip) = alfa2(ie)*ct(i) + beta2(ie)*ctgl(j) &
                           + gamma2(ie)*ct(i)*ctgl(j) + delta2(ie)                        

                        dg_els(nel_dg_loc)%wx_pl(ip) = wwgl(j)
                        dg_els(nel_dg_loc)%wy_pl(ip) = 1.d0


                     enddo
                   enddo
                   


                   
                   deallocate(ct,ww,dd)
                   deallocate(ctgl,wwgl)
               endif               

               !3rd edge
               ne1 = cs_loc(cs_loc(ie -1) +nn*nn)
               ne2 = cs_loc(cs_loc(ie -1) +nn*(nn -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0) ) then

                  call GET_TAG_BC(con_bc, nface, ne1, ne2, tag_dg_el, nload_dg, tag_ind)
                  
                  if (tag_ind == 0) print*, 'ERROR IN GET_TAG_BC'


                  nel_dg_loc = nel_dg_loc +1
                  dg_els(nel_dg_loc)%ind_el = ie
                  dg_els(nel_dg_loc)%face_el = 3
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = nofqp

                  dg_els(nel_dg_loc)%proj_yn = tag_dg_yn(tag_ind)
                  
                  call MAKE_NORMAL( xs(ne1), xs(ne2), ys(ne1), ys(ne2), normal_x, normal_y)

                  dg_els(nel_dg_loc)%nx = normal_x
                  dg_els(nel_dg_loc)%ny = normal_y
                              
                 
                  allocate(ct(dg_els(nel_dg_loc)%quad_rule), ww(dg_els(nel_dg_loc)%quad_rule), & 
                            dd(dg_els(nel_dg_loc)%quad_rule,dg_els(nel_dg_loc)%quad_rule))
                            
                  call MAKE_LGL_NW(dg_els(nel_dg_loc)%quad_rule, ct, ww, dd) 
                  
                  allocate(ctgl(dg_els(nel_dg_loc)%quad_rule), wwgl(dg_els(nel_dg_loc)%quad_rule))
                  
                  call MAKE_GL_NW(dg_els(nel_dg_loc)%quad_rule, ctgl, wwgl)   
                  
                  ip = 0
                    do j = dg_els(nel_dg_loc)%quad_rule, dg_els(nel_dg_loc)%quad_rule
                      do i = 1, dg_els(nel_dg_loc)%quad_rule
                        ip = ip + 1

                        scratch_dg_els(nel_dg_loc)%x_nq(ip) = alfa1(ie)*ctgl(i) + beta1(ie)*ct(j) &
                           + gamma1(ie)*ctgl(i)*ct(j) + delta1(ie)
                        
                        scratch_dg_els(nel_dg_loc)%y_nq(ip) = alfa2(ie)*ctgl(i) + beta2(ie)*ct(j) &
                           + gamma2(ie)*ctgl(i)*ct(j) + delta2(ie)                        
                                                    

                        dg_els(nel_dg_loc)%wx_pl(ip) = wwgl(i)
                        dg_els(nel_dg_loc)%wy_pl(ip) = 1.d0
                                               
                     enddo
                   enddo
                   deallocate(ct,ww,dd)
                   deallocate(ctgl,wwgl)
             endif               


               !4th edge
               ne1 = cs_loc(cs_loc(ie -1) +nn*(nn -1) +1)
               ne2 = cs_loc(cs_loc(ie -1) +1)

                            
               if ((i4count(ne1).ne.0) .and. (i4count(ne2).ne.0)) then

                  call GET_TAG_BC(con_bc, nface, ne1, ne2, tag_dg_el, nload_dg, tag_ind)
                  
                  if (tag_ind == 0) print*, 'ERROR IN GET_TAG_BC'

                  nel_dg_loc = nel_dg_loc +1
                  dg_els(nel_dg_loc)%ind_el = ie
                  dg_els(nel_dg_loc)%face_el = 4
                  dg_els(nel_dg_loc)%mat = tag_mat(im)
                  dg_els(nel_dg_loc)%spct_deg = nn-1
                  dg_els(nel_dg_loc)%quad_rule = nofqp

                  dg_els(nel_dg_loc)%proj_yn = tag_dg_yn(tag_ind)
 
                  call MAKE_NORMAL( xs(ne1), xs(ne2), ys(ne1), ys(ne2), normal_x, normal_y)

                  dg_els(nel_dg_loc)%nx = normal_x
                  dg_els(nel_dg_loc)%ny = normal_y
                              
                 
                  allocate(ct(dg_els(nel_dg_loc)%quad_rule),ww(dg_els(nel_dg_loc)%quad_rule), & 
                           dd(dg_els(nel_dg_loc)%quad_rule,dg_els(nel_dg_loc)%quad_rule))
                           
                  call MAKE_LGL_NW(dg_els(nel_dg_loc)%quad_rule, ct, ww, dd) 
                  
                  allocate(ctgl(dg_els(nel_dg_loc)%quad_rule),wwgl(dg_els(nel_dg_loc)%quad_rule))
                  
                  call MAKE_GL_NW(dg_els(nel_dg_loc)%quad_rule, ctgl, wwgl)   
                  
                  ip = 0
                    do j = 1, dg_els(nel_dg_loc)%quad_rule
                      do i = 1,1
                        ip = ip + 1

                        scratch_dg_els(nel_dg_loc)%x_nq(ip) = alfa1(ie)*ct(i) + beta1(ie)*ctgl(j) &
                           + gamma1(ie)*ctgl(i)*ct(j) + delta1(ie)
                        
                        scratch_dg_els(nel_dg_loc)%y_nq(ip) = alfa2(ie)*ct(i) + beta2(ie)*ctgl(j) &
                           + gamma2(ie)*ct(i)*ctgl(j) + delta2(ie)                        
                     
                        dg_els(nel_dg_loc)%wx_pl(ip) = 1.d0
                        dg_els(nel_dg_loc)%wy_pl(ip) = wwgl(j)

                     enddo
                   enddo
                   
                   
                   deallocate(ct,ww,dd)
                   deallocate(ctgl,wwgl)
               endif             


            endif
            
            
         enddo
         
      enddo

     


       open(600,file='NORMALL.input')
       write(600,*) nel_dg_loc 
       do i = 1, nel_dg_loc
          write(600,"(1I2,1X,1I12,1X,1I2,1X,3(1X,ES12.4))") &
                       dg_els(i)%mat, dg_els(i)%ind_el, dg_els(i)%face_el, &
                       dg_els(i)%nx, dg_els(i)%ny
       enddo
       close(600)     
     
          
      
      return
      
      end subroutine SETUP_DG_ELEM
