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
!    along with SPEED.  If not, see <http://wqw.gnu.org/licenses/>.

!> @brief Makes internal forces in non-linear case.
!! @author Ilario Mazzieri and Filippo Gatti
!> @date February,2016
!> @version 1.0
!> @param[in] nn number of 1-D GLL nodes
!> @param[in] xq GLL nodes
!> @param[in] wq GLL weights
!> @param[in] dd spectral derivative matrix
!> @param[in] dxdx,dxdy,dydx,dydy jacobian of the bilinear map
!> @param[in] sxx nodal values for the stress tensor
!> @param[in] syy nodal values for the stress tensor
!> @param[in] szz nodal values for the stress tensor
!> @param[in] sxy nodal values for the stress tensor
!> @param[out] fx x-componnent for internal forces
!> @param[out] fy y-componnent for internal forces


subroutine MAKE_INTERNAL_FORCE_NL(nn,ct,ww,dd,&             
    dUxdx,dUxdy,dUydx,dUydy,&
    Xkin_ij_lmc,Riso_lmc,&
    sxx,syy,szz_el,sxy,&
    dSxx,dSyy,dSzz,dSxy,&
    mu,lambda,sigma_yld,&
    Ckin,kapakin,Rinf,biso_el,&
    fx,fy)                               
    implicit none
    use nonlinear2d

    integer*4 :: nn
    real*8, dimension(0:nn-1),            intent(in)    :: ct,ww,dUxdx,dUxdy,dUydx,dUydy  
    real*8, dimension(0:nn-1,0:nn-1),     intent(in)    :: dd
    real*8, dimension(0:nn-1,0:nn-1),     intent(inout) :: Riso_lmc
    real*8, dimension(0:nn-1,0:nn-1),     intent(inout) :: dSxx,dSyy,dSzz,dSxy,sxx,syy,szz,sxy,fx,fy
    real*8, dimension(0:2,0:nn-1,0:nn-1), intent(inout) :: Xkin_ij_lmc_el
    real*8, dimension(0:2)                              :: Sigma_ij_start,dSigma_ij_trial
    real*8, dimension(0:2)                              :: dEpsilon_ij_pl,dEpsilon_ij_alpha
    integer*4                                           :: ip,iq,il,im,i
    real*8 :: det_j,t1ux,t1uy,t2ux,t2uy,t1fx,t1fy,t2fx,t2fy
    
    do iq = 1,nn
        do ip = 1,nn
            Sigma_ij_start      = (/sxx(ip,iq),syy(ip,iq),sxy(ip,iq)/)
            Sigma_ij_trial      = (/dSxx(ip,iq),dSyy(ip,iq),dSxy(ip,iq)/)
            dEpsilon_ij_alpha   = (/dUxdx(ip,iq),dUydy(ip,iq),(dUxdy(ip,iq)+dUydx(ip,iq))/)
            dEpsilon_ij_pl(0:2) = 0d0
            call check_plasticity (Sigma_ij_trial, Sigma_ij_start, Xkin_ij_N, Riso_N, &
                sigma_yld,st_epl,alpha_elp,i,j,k,nelement)
            if (st_epl == 1) then                                                                                                                              write(*,*) "1-alpha",1-alpha_elp 
                call plastic_corrector(dEpsilon_ij_alpha, Sigma_ij_trial, Xkin_ij_N, sigma_yld, &                                                                  Riso_N, b_iso, Rinf_iso, C_kin, kapa_kin, xmu, xla, dEpsilon_ij_pl)
            end if
            sxx = Sigma_ij_trial(0)
            syy = Sigma_ij_trial(1)
            sxy = Sigma_ij_trial(3)
            szz = Sigma_ij_trial(4)
            ! FORCE CALCULATION

            t1fx = 0.0d0; t1fy = 0.0d0
            t2fx = 0.0d0; t2fy = 0.0d0

            do il = 1,nn
               t1fx = t1fx + dd(il,ip) * ww(il)*ww(iq) &
                    * (sxx(il,iq)*dydy(il) - sxy(il,iq)*dxdy(il))
               t1fy = t1fy + dd(il,ip) * ww(il)*ww(iq) &
                    * (sxy(il,iq)*dydy(il) - syy(il,iq)*dxdy(il))
            enddo

            do im = 1,nn
               t2fx = t2fx + dd(im,iq) * ww(ip)*ww(im) &
                    * (sxx(ip,im)*dydx(im) - sxy(ip,im)*dxdx(im))
               t2fy = t2fy + dd(im,iq) * ww(ip)*ww(im) &
                    * (sxy(ip,im)*dydx(im) - syy(ip,im)*dxdx(im))
            enddo

            det_j = dxdx(iq)*dydy(ip) - dxdy(ip)*dydx(iq)

            fx(ip,iq) = t1fx - t2fx
            fy(ip,iq) = t1fy - t2fy
        enddo
    enddo


    return

end subroutine MAKE_INTERNAL_FORCE

