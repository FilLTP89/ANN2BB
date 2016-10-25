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

!> @brief Computes DG matrices for jumps.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
!> @param[in] x_pl  x-coord of Gauss Legendre node for the integration in Omega+
!> @param[in] y_pl  y-coord of Gauss Legendre node for the integration in Omega+
!> @param[in] wx_pl Gauss Legendre weight
!> @param[in] wy_pl Gauss Legendre weight
!> @param[in] x_mn  x-coord. of Gauss Legendre nodes for the integration in Omega-
!> @param[in] y_mn  y-coord. of Gauss Legendre nodes for the integration in Omega-
!> @param[in] nq  number of quadrature pints
!> @param[in] nn  polynomial degree in Omega+
!> @param[in] mm  polynomial degree in Omega-
!> @param[in] o_minus  matrix containing info about Omega- o_minus(0,i) = material, 
!!            o_minus(1,i) = el index, o_minus(2,i) = face, o_minus(3,1) = neighbouring element in 
!!            a local numeration
!> @param[in] alfa11 costant values for the bilinear map of Omega+
!> @param[in] alfa12 costant values for the bilinear map of Omega+
!> @param[in] beta11 costant values for the bilinear map of Omega+
!> @param[in] beta12 costant values for the bilinear map of Omega+
!> @param[in] gamma1 costant values for the bilinear map of Omega+
!> @param[in] gamma2 costant values for the bilinear map of Omega+
!> @param[in] delta1 costant values for the bilinear map of Omega+
!> @param[in] delta2 costant values for the bilinear map of Omega+
!> @param[in] ine index for neighbouring element
!> @param[in] alfa11_mn costant values for the bilinear map of Omega-
!> @param[in] alfa12_mn costant values for the bilinear map of Omega-
!> @param[in] beta11_mn costant values for the bilinear map of Omega-
!> @param[in] beta12_mn costant values for the bilinear map of Omega-
!> @param[in] gamma1_mn costant values for the bilinear map of Omega-
!> @param[in] gamma2_mn costant values for the bilinear map of Omega-
!> @param[in] delta1_mn costant values for the bilinear map of Omega-
!> @param[in] delta2_mn costant values for the bilinear map of Omega-
!> @param[in] cp_a constant appearing in the interface integrals
!> @param[in] cp_b constant appearing in the interface integrals
!> @param[in] cp_c constant appearing in the interface integrals
!> @param[in] cp_e constant appearing in the interface integrals
!> @param[in] cp_f constant appearing in the interface integrals
!> @param[in] cp_g constant appearing in the interface integrals
!> @param[in] pen  penalty constant
!> @param[in] dg_cnst  -1 SIPG, 0 IIPG, 1 NIPG
!> @param[in] det_scal  determinant of the transformation used for computing interface integrals
!> @param[out] JP  jump matrix for the part +,+
!> @param[out] JM  jump matrix for the part +,-

      subroutine MAKE_LOC_MATRIX_DG(x_pl, y_pl, &
                       wx_pl,wy_pl,&
                       x_mn, y_mn,  &
                       nq, nn, mm, &
                       o_minus, &
                       alfa1,alfa2, &
                       beta1,beta2,&
                       gamma1,gamma2,&
                       delta1,delta2,&
                       ine,&
                       alfa1_mn,alfa2_mn,&
                       beta1_mn,beta2_mn,&
                       gamma1_mn,gamma2_mn,&
                       delta1_mn,delta2_mn,&
                       cp_a,cp_b,cp_c,cp_e,cp_f,cp_g, &
                       pen, dg_cnst,JP,JM,det_scal,testmode,&
                       JP_uv, JM_uv)


     use max_var
     use str_mesh
     
     implicit none
     
     
     integer*4 :: i, testmode
     integer*4 :: ip1 , i1, j1, h1, m1, n1, p1, L1, K1
     integer*4 :: ip2 , i2, j2, h2, m2, n2, p2, L2, K2
     integer*4 :: ip3 , i3, j3, h3, m3, n3, p3, L3, K3
     !integer*4 :: nthreads,  OMP_GET_NUM_THREADS

     integer*4, intent(in) :: nn, nq, mm, ine  

     integer*4, dimension(1:nq,0:3), intent(in) :: o_minus

     real*8 :: phi, dphi
     real*8 :: dxdx,dxdy,dydx, dydy, det_j, det_s

     real*8 :: dxdx_m,dxdy_m, dydx_m,dydy_m,det_j_m

     real*8 :: dpdc1, dpde1,  dpdx1, dpdy1,  plx1, pkx1
     real*8 :: dpdc2, dpde2,  dpdx2, dpdy2,  plx2, pkx2
     real*8 :: dpdc3, dpde3,  dpdx3, dpdy3,  plx3, pkx3


     real*8, intent(in) :: alfa1, alfa2
     real*8, intent(in) :: beta1, beta2
     real*8, intent(in) :: gamma1, gamma2, delta1, delta2
     real*8, intent(in) :: alfa1_mn, alfa2_mn
     real*8, intent(in) :: beta1_mn, beta2_mn
     real*8, intent(in) :: gamma1_mn, gamma2_mn
     real*8, intent(in) :: delta1_mn, delta2_mn
     real*8, intent(in) :: cp_a, cp_b, cp_c, cp_e, cp_f, cp_g
     real*8, intent(in) :: pen, dg_cnst, det_scal
     real*8, dimension(nq), intent(in) :: x_pl, y_pl, wx_pl, wy_pl, x_mn, y_mn 
     real*8, dimension(nn) :: ctp
     real*8, dimension(mm) :: ctm

     real*8, dimension(:,:), allocatable :: DD,EE,GG,HH,QQ
     real*8, dimension(:,:), allocatable :: AA,BB,PP
     real*8, dimension(2*nn**2,2*nn**2), intent(out)  :: JP, JP_uv  
     real*8, dimension(2*nn**2,2*mm**2), intent(out)  :: JM, JM_uv
    
    call MAKE_LGL_NODES(nn,ctp)
    call MAKE_LGL_NODES(mm,ctm)

!********************************************************************************************!   
!   MATRICES (+,+)                                                                           !
!********************************************************************************************!
     
     
     
   allocate(AA(nn**2,nn**2),BB(nn**2,nn**2),PP(nn**2,nn**2))

    AA = 0.d0; BB = 0.d0; PP = 0.d0;
 

    allocate(DD(nn**2,mm**2),EE(nn**2,mm**2),QQ(nn**2,mm**2))

    DD = 0.d0; EE = 0.d0; QQ = 0.d0



    allocate(GG(nn**2,mm**2),HH(nn**2,mm**2))
    GG = 0.d0; HH = 0.d0;

        do n1 = 1, nn
          do m1 = 1, nn
             L1 = m1 + (n1-1)*nn 
     
     
                do j1 = 1, nn
                   do i1 = 1, nn
                      K1 = i1 + (j1-1)*nn 
                        
                      do  ip1 = 1, nq   
                      
                          if( o_minus(ip1,1) .eq. ine) then
       
                        !ip1
                        dxdy = beta1 + gamma1 * x_pl(ip1)
                        dydy = beta2 + gamma2 * x_pl(ip1)
                        dxdx = alfa1 + gamma1 * y_pl(ip1)
                        dydx = alfa2 + gamma2 * y_pl(ip1)

               
                        det_j = dxdx * dydy - dxdy * dydx

              dpdc1 = dphi(nn-1, ctp(i1), x_pl(ip1)) * phi(nn-1, ctp(j1), y_pl(ip1)) 
              dpde1 = phi(nn-1, ctp(i1), x_pl(ip1)) * dphi(nn-1, ctp(j1), y_pl(ip1)) 

               plx1 = phi(nn-1, ctp(m1), x_pl(ip1)) * phi(nn-1, ctp(n1), y_pl(ip1)) 
               pkx1 = phi(nn-1, ctp(i1), x_pl(ip1)) * phi(nn-1, ctp(j1), y_pl(ip1)) 


              dpdx1 = (dydy*dpdc1 - dydx*dpde1)/det_j
 
              dpdy1 = (dxdx*dpde1 - dxdy*dpdc1)/det_j

               det_s = det_scal/2.d0
              

              AA(L1,K1) = AA(L1,K1)  +  det_s * wx_pl(ip1)*wy_pl(ip1) * dpdx1 * plx1  
              BB(L1,K1) = BB(L1,K1)  +  det_s * wx_pl(ip1)*wy_pl(ip1) * dpdy1 * plx1 
              PP(L1,K1) = PP(L1,K1)  +  det_s * wx_pl(ip1)*wy_pl(ip1) * pkx1 * plx1  
              
              

                           endif
               
  
                      enddo   ! end loop on ip1

                  enddo
               enddo
            
           
        enddo
     enddo
         
!********************************************************************************************!   
!   MATRICES (-,+)                                                                           !
!********************************************************************************************!

        do n2 = 1, nn
          do m2 = 1, nn
             L2 = m2 + (n2-1)*nn 
     
     
                do j2 = 1, mm
                   do i2 = 1, mm
                      K2 = i2 + (j2-1)*mm 
                    
      
     
     
                      do  ip2 = 1, nq  

                          if( o_minus(ip2,1) .eq. ine) then
       
                       dxdy = beta1 + gamma1 * x_pl(ip2)
                       dydy = beta2 + gamma2 * x_pl(ip2)
                       dxdx = alfa1 + gamma1 * y_pl(ip2)
                       dydx = alfa2 + gamma2 * y_pl(ip2)


                       dxdy_m = beta1_mn + gamma1_mn * x_mn(ip2)
                       dydy_m = beta2_mn + gamma2_mn * x_mn(ip2)
                       dxdx_m = alfa1_mn + gamma1_mn * y_mn(ip2)
                       dydx_m = alfa2_mn + gamma2_mn * y_mn(ip2)
               


                        det_j_m = dxdx_m * dydy_m - dxdy_m * dydx_m


              dpdc2 = dphi(mm-1, ctm(i2), x_mn(ip2)) * phi(mm-1, ctm(j2), y_mn(ip2)) 
              dpde2 = phi(mm-1, ctm(i2), x_mn(ip2)) * dphi(mm-1, ctm(j2), y_mn(ip2)) 

                plx2 = phi(nn-1, ctp(m2), x_pl(ip2)) * phi(nn-1, ctp(n2), y_pl(ip2)) 
                pkx2 = phi(mm-1, ctm(i2), x_mn(ip2)) * phi(mm-1, ctm(j2), y_mn(ip2)) 

              dpdx2 = (dydy_m*dpdc2 - dydx_m*dpde2)/det_j_m
 
              dpdy2 = (dxdx_m*dpde2 - dxdy_m*dpdc2)/det_j_m

              det_s = det_scal/2.d0            

              DD(L2,K2) = DD(L2,K2)  +  det_s * wx_pl(ip2)*wy_pl(ip2) * dpdx2 * plx2  
              EE(L2,K2) = EE(L2,K2)  +  det_s * wx_pl(ip2)*wy_pl(ip2) * dpdy2 * plx2  
              QQ(L2,K2) = QQ(L2,K2)  +  det_s * wx_pl(ip2)*wy_pl(ip2) * pkx2 * plx2  
              
              

                           endif

  
                      enddo   ! end loop on ip2


               enddo
            enddo  !end loop on K
            
     enddo
  enddo    !end loop on L

              
!********************************************************************************************!   
!   MATRICES (+,-)                                                                           !
!********************************************************************************************!

!!$OMP SECTION

        do n3 = 1, nn
          do m3 = 1, nn
             L3 = m3 + (n3-1)*nn
     
     
                do j3 = 1, mm
                   do i3 = 1, mm
                      K3 = i3 + (j3-1)*mm 
                    
                      do  ip3 = 1, nq   
                      
                     
                          if( o_minus(ip3,1) .eq. ine) then
       
                        dxdy = beta1 + gamma1 * x_pl(ip3)
                        dydy = beta2 + gamma2 * x_pl(ip3)
                        dxdx = alfa1 + gamma1 * y_pl(ip3)
                        dydx = alfa2 + gamma2 * y_pl(ip3)

               
                        det_j = dxdx * dydy - dxdy * dydx


              dpdc3 = dphi(nn-1, ctp(m3), x_pl(ip3)) * phi(nn-1, ctp(n3), y_pl(ip3)) 
              dpde3 = phi(nn-1, ctp(m3), x_pl(ip3)) * dphi(nn-1, ctp(n3), y_pl(ip3)) 

                pkx3 = phi(mm-1, ctm(i3), x_mn(ip3)) * phi(mm-1, ctm(j3), y_mn(ip3)) 


              dpdx3 = (dydy*dpdc3 - dydx*dpde3)/det_j
              dpdy3 = (dxdx*dpde3 - dxdy*dpdc3)/det_j
                      
              det_s = det_scal/2.d0
 
              GG(L3,K3) = GG(L3,K3)  +  det_s * wx_pl(ip3)*wy_pl(ip3) * dpdx3 * pkx3  
              HH(L3,K3) = HH(L3,K3)  +  det_s * wx_pl(ip3)*wy_pl(ip3) * dpdy3 * pkx3  

                           endif

  
                      enddo   ! end loop on ip3

               enddo
            enddo  !end loop on K
            
     enddo
  enddo    !end loop on L
                  



  JP = 0.d0          
  JP(1 : nn**2, 1 : nn**2) =  pen*PP - cp_a*AA - cp_c*BB &
               + dg_cnst*cp_a*transpose(AA) + dg_cnst*cp_c*transpose(BB) 
  
  JP(1 : nn**2, nn**2+1 : 2*nn**2) = - cp_c*AA - cp_b*BB &
               + dg_cnst*cp_g*transpose(AA) + dg_cnst*cp_e*transpose(BB) 
               
  JP(nn**2+1 : 2*nn**2, 1 : nn**2) = - cp_g*AA - cp_e*BB &
               + dg_cnst*cp_c*transpose(AA) + dg_cnst*cp_b*transpose(BB) 
  
  JP(nn**2+1 : 2*nn**2, nn**2+1 : 2*nn**2) = pen*PP - cp_e*AA - cp_f*BB  &
               + dg_cnst*cp_e*transpose(AA) + dg_cnst*cp_f*transpose(BB) 
  

  JM = 0.d0
  JM(1 : nn**2, 1 : mm**2) = - pen*QQ - cp_a*DD - cp_c*EE  &
                             - dg_cnst*cp_a*GG -dg_cnst*cp_c*HH 
  
  JM(1 : nn**2, mm**2+1 : 2*mm**2) = - cp_c*DD - cp_b*EE - dg_cnst*cp_g*GG - dg_cnst*cp_e*HH  
             
  
  JM(nn**2+1 : 2*nn**2, 1 : mm**2) = - cp_g*DD - cp_e*EE - dg_cnst*cp_c*GG - dg_cnst*cp_b*HH
  
  JM(nn**2+1 : 2*nn**2, mm**2+1 : 2*mm**2) = - pen*QQ - cp_e*DD - cp_f*EE &
                                             - dg_cnst*cp_e*GG - dg_cnst*cp_f*HH 
  
  

   if(testmode .eq. 1) then

     JP_uv = 0.d0          
     JP_uv(1 : nn**2, 1 : nn**2) =  pen*PP 
     JP_uv(nn**2+1 : 2*nn**2, nn**2+1 : 2*nn**2) = pen*PP 
  
     JM_uv = 0.d0
     JM_uv(1 : nn**2, 1 : mm**2) = - pen*QQ 
     JM_uv(nn**2+1 : 2*nn**2, mm**2+1 : 2*mm**2) = - pen*QQ 
     
   endif



   deallocate(AA,BB,PP)            
   deallocate(DD,EE,QQ)                  
   deallocate(GG,HH) 


   end subroutine MAKE_LOC_MATRIX_DG
