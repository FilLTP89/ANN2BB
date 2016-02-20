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
!> @param[in] ct LGL nodes
!> @param[in] ww LGL weights
!> @param[in] dd matrix of spectral derivatives
!> @param[in] dUxdx,dUxdy,dUydx,dUydy nodal values for strain tensor (on element)
!> @param[inout] sxx,syy,sxy,szz nodal values for the stress tensor (on element)
!> @param[inout] Xkin_el,Riso_el hardening variables on element LGL
!> @param[in] lambda_el,mu_el Lamè parameters on element LGL
!> @param[in] syld_el yield limit on element LGL
!> @param[in] Ckin_el,kkin_el kinematic hardening parameters on element LGL
!> @param[in] Rinf_el,biso_el isotropic hardening parameters on element LGL
!> @param[inout] dEpl_el plastic strain increment on element LGL
!> @param[out] fx x-componnent for internal forces
!> @param[out] fy y-componnent for internal forces


subroutine MAKE_INTERNAL_FORCE_NL(nn,ct,ww,dd,      &             
    dUxdx,dUxdy,dUydx,dUydy,sxx,syy,szz,sxy,        &
    Xkin_el,Riso_el,mu_el,lambda_el,syld_el,        &
    Ckin_el,kkin_el,Rinf_el,biso_el,dEpl_el,        &
    dxdx,dxdy,dydx,dydy,fx,fy)                               
    
    use nonlinear2d
    
    implicit none

    integer*4                                 :: ip,iq,il,im,i
    real*8                                    :: t1fx,t1fy,t2fx,t2fy
    real*8                                    :: det_j,t1ux,t1uy,t2ux,t2uy
    
    integer*4,                  intent(in)    :: nn
    logical                                   :: st_epl
    real*8                                    :: alpha_elp
    real*8, dimension(4)                      :: dEalpha,Sstart,dStrial,Strial
    real*8, dimension(nn),      intent(inout) :: ct,ww,dxdx,dxdy,dydx,dydy
    real*8, dimension(nn,nn),   intent(inout) :: Ckin_el,kkin_el
    real*8, dimension(nn,nn),   intent(inout) :: Rinf_el,biso_el,Riso_el
    real*8, dimension(nn,nn),   intent(inout) :: dUxdx,dUxdy,dUydx,dUydy  
    real*8, dimension(nn,nn),   intent(inout) :: dd,lambda_el,mu_el,syld_el
    real*8, dimension(nn,nn),   intent(inout) :: sxx,syy,sxy,szz,fx,fy
    real*8, dimension(4,nn,nn), intent(inout) :: Xkin_el,dEpl_el

    do iq = 1,nn
        do ip = 1,nn
            
            ! FIRST MECHANISM XY
            Sstart  = (/sxx(ip,iq),syy(ip,iq),szz(ip,iq),sxy(ip,iq)/)
            dEalpha = (/dUxdx(ip,iq),dUydy(ip,iq),0d0,(dUxdy(ip,iq)+dUydx(ip,iq))/)

            ! COMPUTE TRIAL STRESS INCREMENT
            call MAKE_STRESS_LOC(lambda_el(ip,iq),mu_el(ip,iq),dEalpha,Strial)
            dStrial=Strial-Sstart

            call check_plasticity (Strial, Sstart, Xkin_el(:,ip,iq), Riso_el(ip,iq), &
                syld_el(ip,iq),st_epl,alpha_elp)
            
            if (st_epl) then
                write(*,*) "1-alpha",1-alpha_elp 
                call plastic_corrector(dEalpha,Strial,syld_el(ip,iq), &
                    Xkin_el(:,ip,iq),Riso_el(ip,iq),biso_el(ip,iq),   &
                    Rinf_el(ip,iq),Ckin_el(ip,iq),kkin_el(ip,iq),     &
                    mu_el(ip,iq),lambda_el(ip,iq),dEpl_el(:,ip,iq))
            end if
            sxx(ip,iq) = Strial(1)
            syy(ip,iq) = Strial(2)
            szz(ip,iq) = Strial(3)
            sxy(ip,iq) = Strial(4)
        end do
    end do

    do iq = 1,nn
        do ip = 1,nn

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

end subroutine MAKE_INTERNAL_FORCE_NL
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

