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

!> @brief Performs the Newton Rapson algorithm.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] xnod x-coord of the point
!> @param[in] ynod y-coord of the point
!> @param[in] alfa1 costant for the bilinear map
!> @param[in] alfa2 costant for the bilinear map
!> @param[in] beta1 costant for the bilinear map
!> @param[in] beta2 costant for the bilinear map
!> @param[in] gamma1  costant for the bilinear map
!> @param[in] gamma2  costant for the bilinear map
!> @param[in] delta1  costant for the bilinear map
!> @param[in] delta2  costant for the bilinear map
!> @param[in] tt   1 if the node belongs to the element, 0 otherwise
!> @param[in] nofi  number of iterations (dummy)
!> @param[in] ipiu control parameter (dummy)
!> @param[in] imeno control parameter (dummy)
!> @param[in] toll1 fixed tolerance
!> @param[in] toll2  fixed tolerance
!> @param[in] flag  control parameter
!> @param[out] csi_s x-coordinate in the reference element of xnod,ynod
!> @param[out] eta_s y-coordinate in the reference element of xnod,ynod


     subroutine NEWTON_RAPSON(xnod, ynod,  &
                           alfa1,alfa2, &
                           beta1,beta2, & 
                           gamma1,gamma2,&
                           delta1,delta2,& 
                           tt, csi_s, eta_s, nofi, &
                           ipiu, imeno, toll1, toll2, flag)

      implicit none
                      
      integer*4 :: nofiter
      integer*4 :: per_csi, per_eta
      integer*4 ::  iic, il,ih
      integer*4, intent(out) :: tt 
      integer*4, intent(in) :: nofi,ipiu,imeno, flag
      
      real*8 :: ll,lm ,xvt1,yvt1,xvt2,yvt2,xvc1,yvc1,xvc2,yvc2
      real*8 :: a1,a2,b1,b2,c1,c2,d1,d2
      real*8 :: delta_csi, delta_eta 
      real*8 :: A,B,C,D, det
      real*8 :: FF1, FF2
      real*8 :: alfa1,alfa2
      real*8 :: beta1,beta2
      real*8 :: gamma1,gamma2,delta1,delta2
      real*8, intent(in) :: xnod, ynod, toll1, toll2
      real*8, intent(out) ::  csi_s, eta_s


      a1 = 4.d0*gamma1
      a2 = 4.d0*gamma2
      b1 = 4.d0*alfa1
      b2 = 4.d0*alfa2
      c1 = 4.d0*beta1
      c2 = 4.d0*beta2
      d1 = 4.d0*delta1
      d2 = 4.d0*delta2

      csi_s = 0.d0
      eta_s = 0.d0
      delta_csi = 0.0d0
      delta_eta = 0.0d0 
      
      tt = 0
      nofiter = 2000

      per_csi = 0
      per_eta = 0


      do il = 1, nofiter
                
        A = a1*eta_s + b1
        B = a1*csi_s + c1 
        C = a2*eta_s + b2
        D = a2*csi_s + c2

        det = A*D-B*C
        
        if(det .eq. 0.d0) then
          return
        endif

        
        FF1 = a1*csi_s*eta_s + b1*csi_s + c1*eta_s + d1 - 4.d0*xnod 
        
        FF2 = a2*csi_s*eta_s + b2*csi_s + c2*eta_s + d2 - 4.d0*ynod 
        
        delta_csi = (-1.d0/det)*(D*FF1-B*FF2)

        delta_eta = (-1.d0/det)*(-C*FF1+A*FF2)



         
        if (abs(delta_csi).le. toll1 .and. abs(delta_eta).le. toll1 ) then
            tt = 1
           
           if(flag.eq.1) then 
              if (dabs(csi_s) .gt. toll2 .or. dabs(eta_s) .gt. toll2 ) then
              tt = 0
              return
            endif
          endif
         
          return
        endif
        
       
        csi_s = csi_s + delta_csi
        eta_s = eta_s + delta_eta
        
        if (csi_s .gt. 1.d0 .AND. csi_s .lt. 2.d0 .AND. per_csi .eq. 0) then
            csi_s = 1.d0 
            per_csi = 1
        endif
             
        if (csi_s .lt. -1.d0 .AND. csi_s .gt. -2.d0 .AND. per_csi .eq. 0) then
            csi_s = -1.d0
            per_csi = 1
        endif
        
        
        if (eta_s .gt. 1.d0 .AND. eta_s .lt. 2.d0 .AND. per_eta .eq. 0) then
            eta_s = 1.d0 
            per_eta = 1
        endif
            
        if (eta_s .lt. -1.d0 .AND. eta_s .gt. -2.d0 .AND. per_eta .eq. 0) then
            eta_s = -1.d0
            per_eta = 1
        endif    
            

      enddo   
      
      if(il .ge. nofiter) then      
        tt = 0
        return
      endif

     end subroutine NEWTON_RAPSON
