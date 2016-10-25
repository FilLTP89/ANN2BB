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

!> @brief Computes norms of error in test mode case.
!! @author Ilario Mazzieri
!> @date April, 2014 
!> @version 1.0
!> @param[in] ... 
!> @param[out] ...

     subroutine COMPUTE_ENERGY_ERROR(nnod, u1, v1, time, ne, cs_loc, cs_nnz_loc,&
                              nm,prop_mat, sdeg_mat,tag_mat, &
                              nelem_dg, &
                              alfa1,alfa2,beta1,beta2, &
                              gamma1,gamma2,delta1,delta2, &
                              xs,ys, &
                              IDG_only_uv, JDG_only_uv, MDG_only_uv, nnz_dg_only_uv)
 

     !use max_var
     !use DGJUMP
     
     implicit none
      
     !type(el4loop), dimension(ne_loc), intent(in) :: el_new 
     
     character*100000 :: input_line
     character*4 :: keyword
     
     integer*4 :: nnod, nfunc, cs_nnz_loc, nnz_dg_only_uv
     integer*4 :: nm, ne, status, ic, id, jb, ib, isb
     integer*4 :: ielem, mm, nelem_dg
     integer*4 :: ie, im, nn, is, in, i, j, ileft,iright,nval
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     integer*4, dimension(0:2*nnod) :: IDG_only_uv     
     integer*4, dimension(nnz_dg_only_uv) :: JDG_only_uv
     real*8, dimension(nnz_dg_only_uv) :: MDG_only_uv


     integer*4, dimension(nm) :: sdeg_mat
     integer*4, dimension(nm) :: tag_mat
               
     real*8 :: time, L2_err, EN_err, term_ener,term_l2 
     real*8 :: Linf_err, EN_err_tot, L2_err_tot, Linf_err_tot, xp, yp, phi

     real*8 :: L2_err_vel, H1_err, term_H1, H1_err_tot, term_l2_vel
     real*8 :: L2_err_vel_tot
     real*8 :: e11_ex, e12_ex, e22_ex
     real*8 :: sxx_ex, syy_ex, sxy_ex


     real*8 :: rho, lambda, mu, pi, x, y, sqrt2
     real*8 :: u1_ctm, u2_ctm, v1_ctm, v2_ctm
     real*8 :: dxdx, dxdy, dydx, dydy, det_j, u1_ex, u2_ex, v1_ex, v2_ex

     real*8, dimension(3) :: Linf_array, Linf_array_vel
     real*8, dimension(nm,4) :: prop_mat
     real*8, dimension(ne) :: alfa1,alfa2
     real*8, dimension(ne) :: beta1,beta2
     real*8, dimension(ne) :: gamma1,gamma2,delta1,delta2
     real*8, dimension(2*nnod) :: u1, v1
     real*8, dimension(nnod) :: xs, ys
     
     real*8, dimension(:), allocatable :: ct, ww, ctm, wwm, u_p, jump
     real*8, dimension(:,:), allocatable :: dd, ddm     
     real*8, dimension(:), allocatable :: val_data, basis_functions
     
     real*8, dimension(:), allocatable :: dxdx_el,dxdy_el,dydx_el,dydy_el
     real*8, dimension(:,:), allocatable :: ux_el,uy_el,uz_el
     real*8, dimension(:,:), allocatable :: duxdx_el,duydx_el 
     real*8, dimension(:,:), allocatable :: duxdy_el,duydy_el 
     real*8, dimension(:,:), allocatable :: sxx_el,syy_el,sxy_el,szz_el
     real*8, dimension(:,:), allocatable :: lambda_el,mu_el 

     pi = 4.d0*datan(1.d0)
     sqrt2 = dsqrt(2.d0)


     L2_err = 0.d0; EN_err = 0.d0; Linf_err = 0.d0;
     L2_err_tot = 0.d0; EN_err_tot = 0.d0; Linf_err_tot = 0.d0;
     
     L2_err_vel = 0.d0;  H1_err = 0.d0;       L2_err_vel_tot = 0.d0;  H1_err_tot = 0.d0; 


     do ie = 1,ne
 
         im = cs_loc(cs_loc(ie -1)) 
         nn = sdeg_mat(im) +1;
         !mm = nn ! 8 ! gll nodes for computing norms
         allocate(ct(nn),ww(nn),dd(nn,nn));  !allocate(ctm(mm),wwm(mm),ddm(mm,mm));
         
         allocate(ux_el(nn,nn),uy_el(nn,nn))
     	 allocate(dxdx_el(nn),dxdy_el(nn),dydx_el(nn),dydy_el(nn))          
     	 allocate(duxdx_el(nn,nn),duydx_el(nn,nn)) 
    	 allocate(duxdy_el(nn,nn),duydy_el(nn,nn)) 
     	 allocate(sxx_el(nn,nn),syy_el(nn,nn),sxy_el(nn,nn),szz_el(nn,nn))
      	 allocate(lambda_el(nn,nn),mu_el(nn,nn)) 

         call LGL(nn,ct,ww,dd)
         !call LGL(mm,ctm,wwm,ddm)
 
                      
         rho = prop_mat(im,1); lambda = prop_mat(im,2); mu = prop_mat(im,3)

         do j = 1,nn
            do i = 1,nn
               is =  nn*(j -1) +i
                  in = cs_loc(cs_loc(ie -1) + is)

                  ux_el(i,j) = u1(in) 
                  uy_el(i,j) = u1(in+nnod)  

            enddo
         enddo

         do i = 1,nn
            dxdy_el(i) = beta1(ie) + gamma1(ie) * ct(i)
            dydy_el(i) = beta2(ie) + gamma2(ie) * ct(i)
            dxdx_el(i) = alfa1(ie) + gamma1(ie) * ct(i)
            dydx_el(i) = alfa2(ie) + gamma2(ie) * ct(i)
         enddo



         call MAKE_STRAIN(nn,ct,ww,dd,&				
                          dxdx_el,dxdy_el,dydx_el,dydy_el, &		
			  ux_el,uy_el, &					
                          duxdx_el,duxdy_el,duydx_el,duydy_el)

         lambda_el = lambda; mu_el = mu;
        
        
         call MAKE_STRESS(lambda_el,mu_el,nn,&
			  duxdx_el,duxdy_el,duydx_el,duydy_el,&
			  sxx_el,syy_el,szz_el,sxy_el) 
			     

         do j = 1,nn
            do i = 1,nn
          
               dxdx = alfa1(ie) + gamma1(ie)*ct(j)
               dydx = alfa2(ie) + gamma2(ie)*ct(j)
                              
               dxdy = beta1(ie) + gamma1(ie)*ct(i)
               dydy = beta2(ie) + gamma2(ie)*ct(i)
                                 
               det_j = (dxdx*dydy - dydx*dxdy)
                                 
               !Local to global node transformation
               !xp = alfa1(ie) * ct(i) + beta1(ie) * ct(j) &
               !         + gamma1(ie) * ct(i) * ct(j) + delta1(ie)
               !yp = alfa2(ie) * ct(i) + beta2(ie) * ct(j) &
               !         + gamma2(ie) * ct(i) * ct(j) + delta2(ie)

               is = nn*(j -1) +i
               in = cs_loc(cs_loc(ie -1) +is)

               xp = xs(in); yp = ys(in);
               
               !PHD THESIS
               !EXACT SOLUTION
               u1_ex = - dsin(sqrt(2.d0)*pi*time) * (dsin(pi*xp)**2) * dsin(2.d0*pi*yp) 
               u2_ex =   dsin(sqrt(2.d0)*pi*time) * dsin(2.d0*pi*xp) * dsin(pi*yp)**2

               v1_ex = - sqrt(2.d0)*pi*dcos(sqrt(2.d0)*pi*time) * (dsin(pi*xp)**2) * dsin(2.d0*pi*yp) 
               v2_ex =   sqrt(2.d0)*pi*dcos(sqrt(2.d0)*pi*time) * dsin(2.d0*pi*xp) * dsin(pi*yp)**2


               !EXACT STRAIN TENSOR
               e11_ex = -2.d0*pi * dcos(pi*xp) * dsin(pi*xp) * dsin(2.d0*pi*yp) * dsin(pi*sqrt2*time)
               e22_ex =  2.d0*pi * dcos(pi*yp) * dsin(2.d0*pi*xp) * dsin(pi*yp) * dsin(pi*sqrt2*time)
               e12_ex =  pi * dcos(2.d0*pi*xp) * dsin(pi*yp)**2 * dsin(pi*sqrt2*time) &
                       - pi * dcos(2.d0*pi*yp) * dsin(pi*xp)**2 * dsin(pi*sqrt2*time)
              
                            
                       
               sxx_ex = - lambda * (2.d0*pi*dcos(pi*xp) * dsin(pi*xp) * dsin(2.d0*pi*yp) * dsin(pi*sqrt2*time) &
                        - 2.d0*pi * dcos(pi*yp) * dsin(2.d0*pi*xp) * dsin(pi*yp) * dsin(pi*sqrt2*time)) &
                        - 4.d0*pi * mu*dcos(pi*xp) * dsin(pi*xp) * dsin(2.d0*pi*yp) * dsin(pi*sqrt2*time)

               syy_ex = 4.d0*pi*mu * dcos(pi*yp) * dsin(2.d0*pi*xp) * dsin(pi*yp) * dsin(pi*sqrt2*time) - &
                        lambda * (2.d0*pi*dcos(pi*xp) * dsin(pi*xp) * dsin(2.d0*pi*yp) * dsin(pi*sqrt2*time) - &
                        2.d0*pi * dcos(pi*yp) * dsin(2.d0*pi*xp) * dsin(pi*yp) * dsin(pi*sqrt2*time))
                        
               sxy_ex = mu*(2.d0*pi * dcos(2.d0*pi*xp) * dsin(pi*yp)**2 * dsin(pi*sqrt2*time) &
                      - 2.d0*pi  *dcos(2.d0*pi*yp) * dsin(pi*xp)**2 * dsin(pi*sqrt2*time))
 

               u1_ctm = 0.d0; u2_ctm = 0.d0; 
               v1_ctm = 0.d0; v2_ctm = 0.d0; 
              
               u1_ctm = u1(in); u2_ctm = u1(in+nnod);
               v1_ctm = v1(in); v2_ctm = v1(in+nnod);  
                      
               !allocate(basis_functions(mm**3)) ! basis function evaluation
               !do jb = 1, nn
               !   do ib =1, nn
               !      isb = nn*(jb -1) +ib
               !      in = cs_loc(cs_loc(ie -1) +isb)

               !      basis_functions(isb) = phi(nn-1, ct(ib), ctm(i)) * &
               !                             phi(nn-1, ct(jb), ctm(j))

               !      u1_ctm =  u1_ctm +  basis_functions(isb)*u1(in);
               !      u2_ctm =  u2_ctm +  basis_functions(isb)*u1(in+nnod);

               !      v1_ctm =  v1_ctm +  basis_functions(isb)*v1(in);
               !      v2_ctm =  v2_ctm +  basis_functions(isb)*v1(in+nnod);
                                                                                    
               !   enddo
               !enddo
               !deallocate(basis_functions)


               !L2 NORM VELOCITY
               !write(*,*) time,xp,yp
               !write(*,*) u1_ex, u1_ctm, u2_ex, u2_ctm
               !write(*,*) v1_ex, v1_ctm, v2_ex, v2_ctm
               
               term_l2_vel =  det_j * ww(i) * ww(j) * ((v1_ex - v1_ctm)**2.d0 + (v2_ex - v2_ctm)**2.d0)   
               L2_err_vel = L2_err_vel + term_l2_vel 
                      

               !L2 NORM DISPLACEMENT
               term_l2 =  det_j * ww(i) * ww(j) * ((u1_ex - u1_ctm)**2.d0 + (u2_ex - u2_ctm)**2.d0)                     
               L2_err = L2_err + term_l2 


               !LINF NORM - DISPL & VEL
               Linf_array = (/ Linf_err,abs(u1_ex - u1_ctm),abs(u2_ex - u2_ctm) /) 
               Linf_err = maxval(Linf_array)


               !ENERGY NORM
               term_ener =  det_j * ww(i) * ww(j) * ( &
                                   rho * (v1_ex - v1_ctm)**2.d0 &
                                 + rho * (v2_ex - v2_ctm)**2.d0 & 
                                 + (sxx_el(i,j) - sxx_ex) * (duxdx_el(i,j) - e11_ex) &
                                 + (syy_el(i,j) - syy_ex) * (duydy_el(i,j) - e22_ex) &
                                 + 2.d0*(sxy_el(i,j) - sxy_ex) * (0.5d0*duxdy_el(i,j) + 0.5d0*duydx_el(i,j) - e12_ex) &
                                )
                                  
                                  
                EN_err = EN_err + term_ener 
                      
                term_H1 =  det_j * ww(i) * ww(j) * ( &
                                 + (sxx_el(i,j) - sxx_ex) * (duxdx_el(i,j) - e11_ex) &
                                 + (syy_el(i,j) - syy_ex) * (duydy_el(i,j) - e22_ex) &
                                 + 2.d0*(sxy_el(i,j) - sxy_ex) * (0.5d0*duxdy_el(i,j) + 0.5d0*duydx_el(i,j) - e12_ex) &
                                )

                 H1_err = H1_err + term_H1
            
            enddo            
         enddo

         deallocate(ct,ww,dd) !(ctm, wwm, ddm)
         deallocate(ux_el,uy_el)
     	 deallocate(dxdx_el,dydx_el,dxdy_el,dydy_el) 
     	 deallocate(duxdx_el,duydx_el,duxdy_el,duydy_el) 
     	 deallocate(sxx_el,syy_el,sxy_el,szz_el)
         deallocate(lambda_el,mu_el) 

     
     enddo     
     
               
!----------------------------------------------------------------------------------------------------------
!    NOW COMPUTE NON CONFORMING TERM S_e int_e [u-uh][u-uh] 
!----------------------------------------------------------------------------------------------------------     
     if(nelem_dg .gt. 0) then
        allocate(u_p(2*nnod)); u_p = 0.d0
        allocate(jump(2*nnod)); jump = 0.d0
              
        do ie = 1, ne

           do j = 1,nn
              do i = 1,nn
                 is = nn*(j -1) +i
                 in = cs_loc(cs_loc(ie -1) + is)
                 xp = xs(in); yp = ys(in);
                  
                 !PHD THESIS        
                 u1_ex = - dsin(sqrt(2.d0)*pi*time) * (dsin(pi*xp)**2) * dsin(2.d0*pi*yp) 
                 u2_ex =   dsin(sqrt(2.d0)*pi*time) * dsin(2.d0*pi*xp) * dsin(pi*yp)**2

                 u_p(in) = abs(u1(in) - u1_ex)
                 u_p(in+nnod) = abs(u1(in+nnod) - u2_ex)

              enddo
           enddo
        enddo
       
        call MATMUL_SPARSE(MDG_only_uv, nnz_dg_only_uv, JDG_only_uv, IDG_only_uv, jump, 2*nnod, u_p, 2*nnod, 1) 

        EN_err = EN_err + dot_product(jump,u_p)
     
     endif      
      
     Linf_err_tot = dabs(Linf_err)

     L2_err_tot = dabs(L2_err)
     L2_err_vel_tot = dabs(L2_err_vel)

     H1_err_tot = dabs(H1_err)
     EN_err_tot = dabs(En_err)

     if (EN_err_tot .ge. 100.d0) stop

     open(40,file='EN.ERR', position='APPEND')
     write(40,*) time, Linf_err_tot, sqrt(L2_err_tot), sqrt(L2_err_vel_tot)  !, sqrt(H1_err_tot), sqrt(EN_err_tot)
     !read(*,*)


     
     
     
     
     
     
     
     end subroutine COMPUTE_ENERGY_ERROR
