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

!> @brief Makes sitffness matrix in CRS sparse format.
!! @author Ilario Mazzieri
!> @date April,2014
!> @version 1.0
!> @param[in] nnod number of nodes
!> @param[in] length nnzero elements for stiff. matrix
!> @param[in] nelem number of elements
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[in] nm number of materials
!> @param[in] tag_mat label for materials
!> @param[in] prop_mat material properties
!> @param[in] sd pol. degrees
!> @param[in] alfa,beta,gamma constatn for bilinear map
!> @param[inout] I_STIFF,J_STIFF,M_STIFF stiffness matrix in CRS format


      subroutine MAKE_STIFF_MATRIX(nnod,length,I_STIFF,J_STIFF,M_STIFF, &
                                   nelem, cs_nnz,cs, nm, tag_mat, prop_mat, sd, &
                                   alfa1, alfa2, beta1, beta2, gamma1, gamma2)
   
      implicit none
      integer*4 :: nelem, cs_nnz,nm,nn, i, s, r, ll, n,m, kk, ind_c
      integer*4 :: in, is, it, j,h, icol1, icol2, ie, im, jm, lenght
      integer*4 :: irow1, irow2, nnod, length
      integer*4, dimension(nm) :: tag_mat, sd
      real*8, dimension(nm,4) :: prop_mat
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(0:2*nnod) :: I_STIFF 
      integer*4, dimension(length) :: J_STIFF
      
      real*8 :: rho, lambda, mu
      real*8 :: delta_mi, delta_ri, delta_nj, delta_sj
      real*8, dimension(nelem) :: alfa1, alfa2, beta1, beta2, gamma1, gamma2
      real*8, dimension(length) ::  M_STIFF
      real*8, dimension(:), allocatable :: ct, ww
      real*8, dimension(:), allocatable :: dxdx_el, dxdy_el, dydx_el, dydy_el
      real*8, dimension(:,:), allocatable :: dd, det_j
      
      real*8, dimension(:,:), allocatable :: mat_A,mat_B,mat_C,mat_D
      real*8, dimension(:), allocatable :: vec_a,vec_b,vec_c,vec_d,vec_U,vec_V,vec_W,vec_Z
      
      real*8 :: c11,c33,c55,c13,c15,c35
      
      do ie = 1,nelem
      
           im = cs(cs(ie -1) +0)
	   rho = prop_mat(im,1)
	   lambda = prop_mat(im,2)
           mu = prop_mat(im,3)

           c11 = lambda + 2*mu; 
           c55 = mu;
           c15 = 0;
           c35 = 0;

           !if (im .eq. 1) then
              c33 = c11  
              c13 = lambda 
           !else
           !   c33 = c11 * 6.20/16.50;
           !   c13 = lambda * 5/16.50;
           !endif
                         
           !write(*,*) c11, c33, c55, c13, c35
           !read(*,*)



           nn = sd(im) + 1  
           allocate(ct(nn),ww(nn),dd(nn,nn))
           allocate(dxdx_el(nn),dxdy_el(nn),dydx_el(nn),dydy_el(nn))
           allocate(det_j(nn,nn))

           dxdx_el = 0.d0;   dxdy_el = 0.d0
           dydx_el = 0.d0;  dydy_el = 0.d0
           det_j = 0.d0

            call LGL(nn,ct,ww,dd)
               
           do i = 1,nn
              dxdy_el(i) = beta1(ie) + gamma1(ie) * ct(i)
              dydy_el(i) = beta2(ie) + gamma2(ie) * ct(i)
           enddo
                      
           do j = 1,nn
              dxdx_el(j) = alfa1(ie) + gamma1(ie) * ct(j)
              dydx_el(j) = alfa2(ie) + gamma2(ie) * ct(j)
           enddo

 	   do j = 1,nn
              do i = 1,nn
                 det_j(i,j) = dxdx_el(j)*dydy_el(i) - dxdy_el(i)*dydx_el(j)
              enddo
           enddo
              
           allocate(mat_A(nn**2,nn**2),mat_B(nn**2,nn**2),mat_C(nn**2,nn**2),mat_D(nn**2,nn**2))
           allocate(vec_a(nn**2),vec_b(nn**2),vec_c(nn**2),vec_d(nn**2))
           allocate(vec_U(nn**2),vec_V(nn**2),vec_W(nn**2),vec_Z(nn**2))
           vec_U = 0.d0;   vec_V = 0.d0;    vec_W = 0.d0;      vec_Z = 0.d0;
           vec_a = 0.d0;   vec_b = 0.d0;    vec_c = 0.d0;      vec_d = 0.d0;
           mat_A = 0.d0;   mat_B = 0.d0;    mat_C = 0.d0;      mat_D = 0.d0;
            
           do s = 1, nn 
              do r = 1, nn
                 ll = (s-1)*nn + r 

                 do n = 1, nn
                    do m = 1, nn
                       kk = (n-1)*nn +m               
                     
                       do j = 1, nn
                          do i = 1, nn 
                             h = (j-1)*nn + i 
                             if(m.eq.i) then 
                               delta_mi = 1.d0
                             else 
                               delta_mi = 0.d0
                             endif
             
                             if(r.eq.i) then 
                               delta_ri = 1.d0
                             else 
                               delta_ri = 0.d0
                             endif
                             if(n.eq.j) then 
                               delta_nj = 1.d0
                             else 
                               delta_nj = 0.d0
                             endif

                             if(s.eq.j) then 
                               delta_sj = 1.d0
                             else 
                               delta_sj = 0.d0
                             endif
                             vec_U(h) = (dydy_el(i)*dd(i,m)*delta_nj - dydx_el(j)*dd(j,n)*delta_mi )  
                             vec_V(h) = (dydy_el(i)*dd(i,r)*delta_sj - dydx_el(j)*dd(j,s)*delta_ri )  
                             vec_W(h) = (dxdx_el(j)*dd(j,s)*delta_ri - dxdy_el(i)*dd(i,r)*delta_sj )
                             vec_Z(h) = (dxdx_el(j)*dd(j,n)*delta_mi - dxdy_el(i)*dd(i,m)*delta_nj )
                             vec_a(h) = ww(i)*ww(j)*(1.d0/det_j(i,j)) * vec_U(h)*vec_V(h)
                             vec_b(h) = ww(i)*ww(j)*(1.d0/det_j(i,j)) * vec_U(h)*vec_W(h)
                             vec_c(h) = ww(i)*ww(j)*(1.d0/det_j(i,j)) * vec_Z(h)*vec_V(h)
                             vec_d(h) = ww(i)*ww(j)*(1.d0/det_j(i,j)) * vec_Z(h)*vec_W(h)
                          enddo
                      enddo
                   mat_A(ll,kk) = sum(vec_a)
                   mat_B(ll,kk) = sum(vec_b)
                   mat_C(ll,kk) = sum(vec_c)   
                   mat_D(ll,kk) = sum(vec_d)
                 enddo
               enddo
                
            enddo
          enddo
                       
          do j = 1,nn
             do i = 1,nn
                is = nn*(j -1) + i
                in = cs(cs(ie -1) + is)
                irow1 = in
                irow2 = in + nnod


                do n = 1,nn
                   do m = 1,nn
                      it = nn*(n -1) + m
                      jm = cs(cs(ie -1) + it)

                      icol1 = jm;    icol2 = jm + nnod 
                                   
                      call FIND_POSITION(J_STIFF, length, I_STIFF(irow1-1) + 1, I_STIFF(irow1), icol1, ind_c)
!                      M_STIFF(ind_c) = M_STIFF(ind_c) + (lambda+2.d0*mu)*mat_A(is,it) + mu*mat_D(is,it)
                      M_STIFF(ind_c) = M_STIFF(ind_c) + c11*mat_A(is,it) + c15*mat_B(is,it) &
                                                      + c15*mat_C(is,it) + c55*mat_D(is,it)


                      call FIND_POSITION(J_STIFF, length, I_STIFF(irow1-1) + 1, I_STIFF(irow1), icol2, ind_c)
!                      M_STIFF(ind_c) = M_STIFF(ind_c) + lambda*mat_C(is,it) + mu*mat_B(is,it)
                      M_STIFF(ind_c) = M_STIFF(ind_c) + c15*mat_A(is,it) + c55*mat_B(is,it) &
                                                      + c13*mat_C(is,it) + c35*mat_D(is,it) 
 
  
                      call FIND_POSITION(J_STIFF, length, I_STIFF(irow2-1) + 1, I_STIFF(irow2), icol1, ind_c)
!                      M_STIFF(ind_c) = M_STIFF(ind_c)  + lambda*mat_B(is,it) + mu*mat_C(is,it)
                      M_STIFF(ind_c) = M_STIFF(ind_c) + c13*mat_B(is,it) + c15*mat_A(is,it) &
                                                      + c35*mat_D(is,it) + c55*mat_C(is,it) 


                      call FIND_POSITION(J_STIFF, length, I_STIFF(irow2-1) + 1, I_STIFF(irow2), icol2, ind_c)
!                      M_STIFF(ind_c) = M_STIFF(ind_c) + mu*mat_A(is,it) + (lambda+2.d0*mu)*mat_D(is,it)
                      M_STIFF(ind_c) = M_STIFF(ind_c) + c35*mat_B(is,it) + c55*mat_A(is,it) &
                                                      + c33*mat_D(is,it) + c35*mat_C(is,it) 


                                      
                  enddo
               enddo
             enddo
          enddo


         deallocate(mat_A,mat_B,mat_C,mat_D,vec_a,vec_b,vec_c,vec_d,vec_U,vec_V,vec_W,vec_Z)
         deallocate(ct,ww,dd,dxdx_el,dxdy_el,dydx_el,dydy_el,det_j)
    enddo              
    
    end subroutine MAKE_STIFF_MATRIX
