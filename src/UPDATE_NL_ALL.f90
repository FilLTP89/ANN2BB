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
subroutine UPDATE_NL_ALL(ie,nnt,nn,cs_nnz,cs,fk,sxx,syy,szz,sxy,&
    duxdx,duxdy,duydx,duydy,xkin_all,epl_all,riso_all,fx_el,fy_el,&
    sxx_el,syy_el,szz_el,sxy_el,duxdx_el,duxdy_el,duydx_el,duydy_el,xkin_el,&
    riso_el,depl_el,fxs_el,fys_el,sism)
    
    implicit none
    integer*4, intent(in)                               :: nnt,nn,cs_nnz,ie
    integer*4, dimension(0:cs_nnz), intent(in)          :: cs
    real*8, dimension(nn,nn), intent(in)                :: fx_el,fy_el
    real*8, dimension(nn,nn), intent(in)                :: duxdx_el,duxdy_el,duydx_el,duydy_el
    real*8, dimension(nn,nn), intent(in)                :: sxx_el,syy_el,szz_el,sxy_el,riso_el
    real*8, dimension(4,nn,nn), intent(in)              :: xkin_el,depl_el
    real*8, dimension(nnt),   intent(inout)             :: sxx,syy,szz,sxy
    real*8, dimension(nnt),   intent(inout)             :: duxdx,duxdy,duydx,duydy
    real*8, dimension(nnt),   intent(inout)             :: riso_all
    real*8, dimension(2*nnt), intent(inout)             :: fk
    real*8, dimension(4*nnt), intent(inout)             :: xkin_all,epl_all
    real*8, dimension(nn,nn), intent(in),    optional   :: fxs_el,fys_el
    real*8, dimension(2*nnt), intent(inout), optional   :: sism
    integer*4                                           :: in,is,i,j 
    
    do j = 1,nn
        do i = 1,nn
            is = nn*(j -1) +i
            in = cs(cs(ie -1) + is)
            fk(in)               = fk(in)      + fx_el(i,j)
            fk(in+nnt)           = fk(in+nnt)  + fy_el(i,j)
            if (in==11) then
                write(*,*) "====== DEBUG (TIME_LOOP_NL.F90) ======="
                write(*,*) sxx_el(i,j),syy_el(i,j),szz_el(i,j),sxy_el(i,j)
                read(*,*)
            endif
            sxx(in) = sxx_el(i,j)
            syy(in) = syy_el(i,j)
            szz(in) = szz_el(i,j)
            sxy(in) = sxy_el(i,j)
            duxdx(in) = duxdx_el(i,j)
            duydy(in) = duydy_el(i,j)
            duxdy(in) = duxdy_el(i,j)
            duydx(in) = duydx_el(i,j)
            xkin_all(in)         = xkin_el(1,i,j)
            xkin_all(in+nnt)     = xkin_el(2,i,j)
            xkin_all(in+2*nnt)   = xkin_el(3,i,j)
            xkin_all(in+3*nnt)   = xkin_el(4,i,j)
            riso_all(in)         = riso_el(i,j)
            epl_all(in)          = epl_all(in)       + depl_el(1,i,j)
            epl_all(in+nnt)      = epl_all(in+nnt)   + depl_el(2,i,j)
            epl_all(in+2*nnt)    = epl_all(in+2*nnt) + depl_el(3,i,j)
            epl_all(in+3*nnt)    = epl_all(in+3*nnt) + depl_el(4,i,j)
            if (present(sism)) then
                sism(in) = sism(in)         + fxs_el(i,j)
                sism(in+nnt) = sism(in+nnt) + fys_el(i,j)
            endif
        enddo
    enddo
    return
end subroutine UPDATE_NL_ALL
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
