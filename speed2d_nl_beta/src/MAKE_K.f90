!    Copyright (C) 2014 The SPEED FOUNDATION
!    Author: 
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

!> @brief Computes stiffnes matrix only for DRM elements and without assembling
!! @author 
!> @date August, 2015
!> @version 1.0

!> @param[in] nn number of nodes of the element
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[in] ne number of local elements
!> @param[in] nMDRM
!> @param[in] tag_MDRM labels for DRM blocks
!> @param[in] n_el_DRM number of DRM elements
!> @param[in] el_DRM "Macro" nodes matrix of DRM domain
!> @param[in] a1 costant values for the bilinear map
!> @param[in] a2 costant values for the bilinear map
!> @param[in] b1 costant values for the bilinear map
!> @param[in] b2 costant values for the bilinear map
!> @param[in] g1 costant values for the bilinear map
!> @param[in] g2 costant values for the bilinear map
!> @param[in] nm number of materials
!> @param[in] prop_mat material properties (rho,lambda,mu,gamma)
!> @param[in] tag_mat material labels
!> @param[in] sdeg polynomial degree vector
!> @param[out] K_DRM stiffnes matrix of DRM elements


    subroutine  MAKE_K(nn,cs_nnz,cs,ne,nMDRM,tag_MDRM,n_el_DRM,el_DRM, &
                         a1,a2,b1,b2,g1,g2,&
						 nm,prop_mat,tag_mat,sdeg,&
						 K_DRM)
			        

!   © POLIMI, 2006, All Rights Reserved
!   Author: Laura Scandella    
	
	implicit none

    integer*4 :: cs_nnz,ne,nMDRM,nm,nn,n_el_DRM
    integer*4, dimension(0:cs_nnz) :: cs
    integer*4, dimension(nMDRM) :: tag_MDRM  
	real*8, dimension(ne) :: a1,b1,g1,a2,b2,g2
    integer*4, dimension(nm) :: sdeg
    integer*4, dimension(nm) :: tag_mat
	real*8, dimension(nm,4) :: prop_mat
    integer*4, dimension(n_el_DRM,5) :: el_DRM             

	real*8, dimension(:), allocatable :: ct,ww
    real*8, dimension(:,:), allocatable :: dd
	real*8, dimension(:), allocatable :: dxdx,dydx,dxdy,dydy
    real*8 :: det_j
    real*8, dimension(2*nn*nn,2*nn*nn) :: stiff_el
    real*8, dimension(2*nn*nn*n_el_DRM,2*nn*nn) :: K_DRM

	real*8, dimension(2,2) :: stiff_ij
	real*8 :: dfi_i_dcsi_l,dfi_j_deta_l,dfi_i_dcsi_m,dfi_j_deta_m
	real*8 :: dfi_i_dx,dfi_i_dy,dfi_j_dx,dfi_j_dy
    real*8 :: lambda,mu
    integer*4 :: imDRM,im,ie,i,j,lj,li,mj,mi,kj,ki,l,m,cont,el


    allocate(ct(nn),ww(nn),dd(nn,nn))
    allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn))
    
    call lgl(nn,ct,ww,dd)

 !	cont = 0
 !   do  imDRM = 1,nMDRM   !For each MDRM group
 !      do ie = 1,ne       !For each element
 !         if (cs(cs(ie -1) +0).eq.tag_MDRM(imDRM)) then 

     do el = 1,n_el_DRM
	    ie = el_DRM(el,1)
		do  imDRM = 1,nMDRM   
          if (cs(cs(ie -1) +0).eq.tag_MDRM(imDRM)) then 
		    do im = 1,nm  !For each material
			   if (tag_mat(im).eq.tag_MDRM(imDRM)) then				           
                      lambda = prop_mat(im,2)
                      mu = prop_mat(im,3)
			   endif
			enddo
            !Calculation of transform functions
            do i = 1,nn !Calculated with csi; equal to: dxdcj(1,2),dxdcj(2,2)
                  dxdy(i) = b1(ie) + g1(ie) * ct(i)
                  dydy(i) = b2(ie) + g2(ie) * ct(i)
            enddo                       
            do j = 1,nn !Calculated with eta; equal to: dxdcj(1,1),dxdcj(2,1)
                  dxdx(j) = a1(ie) + g1(ie) * ct(j)
                  dydx(j) = a2(ie) + g2(ie) * ct(j)
            enddo

			do lj = 1,nn
			   do li = 1,nn
			      l = (lj-1)*nn+li  !Node number in the reference element
				  do mj =1,nn
				     do mi =1,nn
                        m = (mj-1)*nn+mi  !Node number in the reference element

					    !Partial stiff matrix initialization
                        do i=1,2
						   do j= 1,2
                               stiff_ij(i,j) = 0
						   enddo
						enddo
						!Integration loop on nodes K
					    do kj = 1,nn
						   do ki = 1,nn
						   !Derivation of shape functions respect to local riferimrnt
						   !(equal to spectral derivate if no zero) 
							  if (kj.eq.lj) then
							      dfi_i_dcsi_l = dd(ki,li)
							  else
							      dfi_i_dcsi_l = 0
							  endif
							  if (ki.eq.li) then
							      dfi_j_deta_l = dd(kj,lj)
							  else
							      dfi_j_deta_l = 0
							  endif
							  if (kj.eq.mj) then
							      dfi_i_dcsi_m = dd(ki,mi)
							  else
							      dfi_i_dcsi_m = 0
							  endif
							  if (ki.eq.mi) then
							      dfi_j_deta_m = dd(kj,mj)
							  else
							      dfi_j_deta_m = 0
							  endif
							  !Jacobian determinant
                              det_j = dxdx(kj)*dydy(ki)-dydx(kj)*dxdy(ki)
							  !Derivation of shape functions respect to global riferiment    
							   dfi_i_dx = (dydy(ki)*dfi_i_dcsi_l-dydx(kj)*dfi_j_deta_l)/det_j
                               dfi_i_dy = (dxdx(kj)*dfi_j_deta_l-dxdy(ki)*dfi_i_dcsi_l)/det_j
							   dfi_j_dx = (dydy(ki)*dfi_i_dcsi_m-dydx(kj)*dfi_j_deta_m)/det_j
							   dfi_j_dy = (dxdx(kj)*dfi_j_deta_m-dxdy(ki)*dfi_i_dcsi_m)/det_j

                              !Stiff matrix components
                               stiff_ij(1,1) = stiff_ij(1,1) + det_j*ww(ki)*ww(kj) &
							                   *((lambda+2*mu)*dfi_i_dx*dfi_j_dx &
							                   + mu*dfi_i_dy*dfi_j_dy)

							   stiff_ij(1,2) = stiff_ij(1,2) + det_j*ww(ki)*ww(kj) &
							                   *(lambda*dfi_i_dx*dfi_j_dy &
							                   + mu*dfi_i_dy*dfi_j_dx)

					           stiff_ij(2,1) = stiff_ij(2,1) + det_j*ww(ki)*ww(kj) &
							                   *(lambda*dfi_i_dy*dfi_j_dx &
							                   + mu*dfi_i_dx*dfi_j_dy)

			                   stiff_ij(2,2) = stiff_ij(2,2) + det_j*ww(ki)*ww(kj) &
							                   *((lambda+2*mu)*dfi_i_dy*dfi_j_dy &
							                   + mu*dfi_i_dx*dfi_j_dx)

							 enddo
						  enddo
                          !Element stiff matrix
                          stiff_el(2*l-1,2*m-1) = stiff_ij(1,1)
						  stiff_el(2*l-1,2*m)   = stiff_ij(1,2)
						  stiff_el(2*l,2*m-1)   = stiff_ij(2,1)
						  stiff_el(2*l,2*m)     = stiff_ij(2,2)

						enddo
					 enddo
				  enddo
			   enddo 
               !Stiff matrix of DRM elements without assemble
		!      K_DRM(1+2*nn*nn*cont:2*nn*nn*(1+cont),1:2*nn*nn) = stiff_el
		!	   cont = cont + 1	
			   K_DRM(1+(el-1)*2*nn*nn:el*2*nn*nn,1:2*nn*nn) = stiff_el	 					
          endif
       enddo
   enddo

   return
               
   end subroutine MAKE_K 