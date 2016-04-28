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

!> @brief Update all variables in leap-frog time stepping.
!! @author Filippo Gatti
!> @date February,2016
!> @version 1.0
!> @param[in] nnt,nn            
!> @param[]             
!> @param[]             
!> @param[]             
!> @param[]             
!> @param[]             
!> @param[]             
!> @param[]             
!> @param[]             
!> @param[]             
!> @param[]             
!> @param[]             
!> @param[]             
!> @param[]             
!> @param[]             
subroutine UPDATE_ALL_INCREMENTS(ie,nnt,nn,cs_nnz,cs,fk,mvec,dsxx,dsyy,dszz,dsxy,&
    dduxdx,dduxdy,dduydx,dduydy,dxkin_all,depl_all,driso_all,fx_el,fy_el,&
    sxx_el,syy_el,szz_el,sxy_el,duxdx_el,duxdy_el,duydx_el,duydy_el,xkin_el,&
    riso_el,depl_el,fxs_el,fys_el,sism)
    
    implicit none
    ! intent IN 
    integer*4,                          intent(in)          :: nnt,nn,cs_nnz,ie
    integer*4,  dimension(0:cs_nnz),    intent(in)          :: cs
    real*8,     dimension(nn,nn),       intent(in)          :: fx_el,fy_el
    real*8,     dimension(nn,nn),       intent(in)          :: duxdx_el,duxdy_el,duydx_el,duydy_el
    real*8,     dimension(nn,nn),       intent(in)          :: sxx_el,syy_el,szz_el,sxy_el,riso_el
    real*8,     dimension(4,nn,nn),     intent(in)          :: xkin_el,depl_el
    real*8,     dimension(2*nnt),       intent(in)          :: mvec
    ! intent INOUT
    real*8,     dimension(nnt),         intent(inout)       :: dsxx,dsyy,dszz,dsxy
    real*8,     dimension(nnt),         intent(inout)       :: dduxdx,dduxdy,dduydx,dduydy
    real*8,     dimension(nnt),         intent(inout)       :: driso_all
    real*8,     dimension(2*nnt),       intent(inout)       :: fk
    real*8,     dimension(4*nnt),       intent(inout)       :: dxkin_all,depl_all
    real*8,     dimension(nn,nn),       intent(in),    optional   :: fxs_el,fys_el
    real*8,     dimension(2*nnt),       intent(inout), optional   :: sism
    integer*4                                               :: in,is,i,j 
    
    do j = 1,nn
        do i = 1,nn
            is = nn*(j -1) +i
            in = cs(cs(ie -1) + is)

            dsxx(in)                = dsxx(in) + sxx_el(i,j)
            dsyy(in)                = dsyy(in) + syy_el(i,j)
            dszz(in)                = dszz(in) + szz_el(i,j)
            dsxy(in)                = dsxy(in) + sxy_el(i,j)
            dduxdx(in)              = dduxdx(in) + duxdx_el(i,j)
            dduydy(in)              = dduydy(in) + duydy_el(i,j)
            dduxdy(in)              = dduxdy(in) + duxdy_el(i,j)
            dduydx(in)              = dduydx(in) + duydx_el(i,j)
            dxkin_all(in)           = dxkin_all(in)      + xkin_el(1,i,j)
            dxkin_all(in+nnt)       = dxkin_all(in+nnt)  + xkin_el(2,i,j)
            dxkin_all(in+2*nnt)     = dxkin_all(in+2*nnt)+ xkin_el(3,i,j)
            dxkin_all(in+3*nnt)     = dxkin_all(in+3*nnt)+ xkin_el(4,i,j)
            driso_all(in)           = driso_all(in)      + riso_el(i,j)
            depl_all(in)            = depl_all(in)       + depl_el(1,i,j)
            depl_all(in+nnt)        = depl_all(in+nnt)   + depl_el(2,i,j)
            depl_all(in+2*nnt)      = depl_all(in+2*nnt) + depl_el(3,i,j)
            depl_all(in+3*nnt)      = depl_all(in+3*nnt) + depl_el(4,i,j)
        enddo
    enddo
    return
end subroutine UPDATE_ALL_INCREMENTS
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

