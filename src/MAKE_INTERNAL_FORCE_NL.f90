!    Copyright (C) 2014 The SPEED FOuNDATION
!    Author: Ilario Mazzieri
!
!    This file is part of SPEED.
!
!    SPEED is free software; you can redistribute it and/or modify it
!    under the terms of the GNu Affero General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    SPEED is distributed in the hope that it will be useful, but
!    WITHOuT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICuLAR PURPOSE.  See the GNU
!    Affero General Public License for more details.
!
!    You should have received a copy of the GNu Affero General Public License
!    along with SPEED.  If not, see <http://wqw.gnu.org/licenses/>.

!> @brief Makes internal forces in non-linear case.
!! @author Ilario Mazzieri and Filippo Gatti
!> @date February,2016
!> @version 1.0

subroutine MAKE_INTERNAL_FORCES_NL(dt,u1,ne,cs_nnz,cs,ct,nm,ww,dd,nnt,snl,&
    nl_sism,vel,alfa1,alfa2,beta1,beta2,gamma1,gamma2)
    
    use nonlinear2d
    
    implicit none


    do iq = 1,nn
        do ip = 1,nn
            stress0 =   (/sxx(ip,iq),syy(ip,iq),szz(ip,iq),sxy(ip,iq)/)
            dEalpha =   (/duxdx(ip,iq),duydy(ip,iq),0.0d0,(duxdy(ip,iq)+duydx(ip,iq))/)
            syld    =   syld_el(ip,iq)
            radius  =   riso_el(ip,iq)
            center  =   xkin_el(:,ip,iq)
            dtrial(:) = 0.d0
!
            ! COMPUTE TRIAL STRESS INCREMENT
            call MAKE_STRESS_LOC(lambda_el(ip,iq),mu_el(ip,iq),dEalpha,dtrial)
            call check_plasticity(dtrial,stress0,center,radius,syld,st_epl,alpha_elp,nel)
            ! PLASTIC CORRECTION 
            if (st_epl) then
                write(*,*) "PLASTIC"
                call plastic_corrector(dEalpha,dtrial,center,syld, &
                    radius,biso_el(ip,iq),Rinf_el(ip,iq),Ckin_el(ip,iq),kkin_el(ip,iq),     &
                    mu_el(ip,iq),lambda_el(ip,iq),dEpl_el(:,ip,iq),nel)
            end if
            sxx(ip,iq)  = dtrial(1)-stress0(1)
            syy(ip,iq)  = dtrial(2)-stress0(2)
            szz(ip,iq)  = dtrial(3)-stress0(3)
            sxy(ip,iq)  = dtrial(4)-stress0(4)
            
            xkin_el(:,ip,iq)    = center-Xkin_el(:,ip,iq)

        end do
    end do

    ! FORCE CALCULATION
    do iq = 1,nn
        do ip = 1,nn

            t1fx = 0.0d0
            t2fx = 0.0d0
            t1fy = 0.0d0
            t2fy = 0.0d0

            ! derivatives with respect to eta ( there is a delta(iq,im) )
            do il=1,nn
                t1fx = t1fx+dd(il,ip)*ww(il)*ww(iq)*(sxx(il,iq)*dydy(il)-sxy(il,iq)*dxdy(il))
                t1fy = t1fy+dd(il,ip)*ww(il)*ww(iq)*(sxy(il,iq)*dydy(il)-syy(il,iq)*dxdy(il))
            enddo
            ! derivatives with respect to xi ( there is a delta(ip,il) )
            do im=1,nn
                t2fx=t2fx+dd(im,iq)*ww(ip)*ww(im)*(sxx(ip,im)*dydx(im)-sxy(ip,im)*dxdx(im))
                t2fy=t2fy+dd(im,iq)*ww(ip)*ww(im)*(sxy(ip,im)*dydx(im)-syy(ip,im)*dxdx(im))
            enddo
            fx(ip,iq) = t1fx-t2fx
            fy(ip,iq) = t1fy-t2fy
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

