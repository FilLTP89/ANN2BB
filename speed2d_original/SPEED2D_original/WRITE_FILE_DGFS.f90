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

!> @brief Writes file DGFS.input, containing infos for computing integrals on DG interfaces.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] nm  number of materials
!> @param[in] sd polynomial degree vector 
!> @param[in] tag_mat label for materials
!> @param[in] cs_nnz_loc length of cs_loc
!> @param[in] cs_loc local spectral connectivity
!> @param[in] nn_loc number of nodes
!> @param[in] ne_loc number of elements
!> @param[in] xs x-coordinate of GLL nodes
!> @param[in] ys y-coordinate of GLL nodes
!> @param[in] nel_dg_loc  number of DG elements
!> @param[in] alfa1 bilinear mapping from (-1,1)^2 to the current el
!> @param[in] alfa2 bilinear mapping from (-1,1)^2 to the current el
!> @param[in] beta1 bilinear mapping from (-1,1)^2 to the current el
!> @param[in] beta2 bilinear mapping from (-1,1)^2 to the current el
!> @param[in] gamma1 bilinear mapping from (-1,1)^2 to the current el
!> @param[in] gamma2 bilinear mapping from (-1,1)^2 to the current el
!> @param[in] delta1 bilinear mapping from (-1,1)^2 to the current el
!> @param[in] delta2 bilinear mapping from (-1,1)^2 to the current el
!> @param[in] faces  data structure containing info about DG faces (material, element, face)
!> @param[in] area_nodes  data structure containing info about DG faces (area, constants for the 
!!                   bilinear mapping)
!> @param[in] dg_els  data structure for dg interface elements --> see module.f90
!> @param[in] scratch_dg_els  temporary data structure for dg interface elements --> see module.f90!
!> @param[out] filename (DGFS.input)  file containing info for DG interface integrals


     subroutine WRITE_FILE_DGFS(nm, sd, tag_mat, cs_nnz_loc, cs_loc, &
                           nn_loc, ne_loc,  &
                           xs,ys,&
                           nel_dg_loc, &
                           alfa1, alfa2, &
                           beta1, beta2, &
                           gamma1, gamma2, &
                           delta1, delta2, &
                           faces, area_nodes, dg_els, scratch_dg_els, &
                           filename)

     use max_var
     use str_mesh 
     use str_mesh_scratch                   

     implicit none

     
     type(ELEMENT), dimension(nel_dg_loc), intent(inout) :: dg_els
     type(scratch_ELEMENT), dimension(nel_dg_loc), intent(inout) :: scratch_dg_els
 
     character*70 :: filename
     
     integer*4 :: nm, cs_nnz_loc, nn_loc, ne_loc, nel_dg_loc

     integer*4 :: im, nn, ie, ned, yon
     integer*4 :: el_conf, face_conf, face_found
     integer*4 :: ip, k, j, i, it, ih, ik, tt, indice, node_not_ass, ic, dim1, dim2, jstart
     integer*4 :: error, ncol
     integer*4 :: mpierror, nofel

     integer*4, dimension(2) :: dims, dimsfi
     integer*4, dimension(nm) :: tag_mat, sd
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc


     integer*4, dimension(3,nel_dg_loc) :: faces
     integer*4, dimension(nel_dg_loc,3) :: mat_el_fac

     real*8 :: normal_x, normal_y
     real*8 :: c_alfa1, c_alfa2
     real*8 :: c_beta1, c_beta2
     real*8 :: c_gamma1, c_gamma2, c_delta1, c_delta2
     real*8 :: xnod, ynod, csi, eta

     real*8, dimension(nn_loc) :: xs,ys,zs
     
     real*8, dimension(ne_loc) :: alfa1,alfa2
     real*8, dimension(ne_loc) :: beta1,beta2
     real*8, dimension(ne_loc) :: gamma1,gamma2
     real*8, dimension(ne_loc) :: delta1,delta2

     real*8, dimension(9,nel_dg_loc) :: area_nodes

     real*8, dimension(nel_dg_loc,2) :: normalxyz
     integer :: ios

     open(40,file='NORMALL.input')                        
     read(40,*) nofel
     do i = 1, nofel
       read(40,*) mat_el_fac(i,1),mat_el_fac(i,2),mat_el_fac(i,3) ,&
                       normalxyz(i,1), normalxyz(i,2)
     enddo
     
     close(40)
     
                        
     open(500,file='DGFS.input')  
       
     ic = 0
             
     ned = cs_loc(0) - 1
     node_not_ass = 0

     do it = 1, nel_dg_loc   
                  
        if(dg_els(it)%proj_yn .eq. 1) then 


           
           do i = 1, dg_els(it)%quad_rule

              tt = 0
              ih = 1
                             
              do while (tt.eq.0 .and. ih.le. nel_dg_loc)   
              
                 if( faces(1,ih) .ne. dg_els(it)%mat) then
                              
                    !ik = faces(2,ih)
                    !il = faces(3,ih)
                      
                    !CHECK THE NORMAL !!!
                    call CHECK_NORMAL(dg_els(it)%nx, dg_els(it)%ny, &
                                       faces(1,ih), faces(2,ih), faces(3,ih), &
                                       nel_dg_loc, normalxyz, mat_el_fac, yon)
                    
                    
                    if(yon .eq. 1) then 
   
                       xnod =  scratch_dg_els(it)%x_nq(i)
                       ynod =  scratch_dg_els(it)%y_nq(i)

                       call MAKE_BILINEAR_MAP(area_nodes(2:9,ih), &
                                             c_alfa1, c_alfa2, c_beta1, c_beta2, & 
                                             c_gamma1, c_gamma2, c_delta1, c_delta2)
 
                       call NEWTON_RAPSON(xnod, ynod, &
                                          c_alfa1, c_alfa2, c_beta1, c_beta2, & 
                                          c_gamma1, c_gamma2, c_delta1, c_delta2, &
                                          tt, csi, eta, nofinr,&
                                          dg_els(it)%ind_el, faces(2,ih), 1.d-6, 1.01d0, 1)
                                           
                                           
                        if(tt == 1) then
                       
                           ic = ic + 1
                            
                           write(500,"(1I2,1X,1I12,1X,1I2,1X,1I2,1X,1I12,1X,1I2,2(1X,ES25.16),2(1X,ES25.16))") &
                           dg_els(it)%mat, dg_els(it)%ind_el, dg_els(it)%face_el, &
                           faces(1,ih), faces(2,ih), faces(3,ih), &
                           xnod, ynod, &
                           dg_els(it)%wx_pl(i),dg_els(it)%wy_pl(i)
                           
                        endif ! if (tt == 1)                  
                                           
 
 
                    else
                    
                       tt = 0 
                    
                  
                    endif ! (yon .eq. 1)      
                   
                 endif ! (faces(ih,1) .ne. dg_els(it)%mat)         
              
                 ih = ih + 1
                  
              enddo  ! while tt == 0
    
              !  CHECK ON NODE NOT ASSIGNED
              if(tt .eq. 0 .and. ih .gt. nel_dg_loc) then
                  node_not_ass = node_not_ass + 1
                  write(*,*) 'NODE', i,' NOT ASSIGNED', ' el', dg_els(it)%ind_el
              endif
                 

           enddo ! i         
       
        endif !if proj_yn == 1
       
     enddo ! it 

         
     dims(1) = 6
     dims(2) = ic
     
     close(500)

     
     if(ic .ne. 0) then 
        write(*,'(A,I5,A,I6)') 'Not assigned nodes : ', node_not_ass
     endif                                                              
        
    return
    
    
     end subroutine WRITE_FILE_DGFS
