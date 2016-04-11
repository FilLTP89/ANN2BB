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

!> @brief Initialize variables for nonlinear calculations 
!! @author Filippo Gatti
!> @date February,2016
!> @version 1.0
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[in] nm number of materials
!> @param[in] nn number of 1-D GLL nodes
!> @param[in] nnt number of nodes
!> @param[in] im material number
!> @param[in] ie element identity
!> @param[in] prop_mat material properties
!> @param[inout] sxx_el,syy_el,sxy_el,szz_el nodal values for the stress tensor (on element)
!> @param[inout] Xkin_el,Riso_el hardening variables on element LGL
!> @param[in] lambda_el,mu_el Lamè parameters on element LGL
!> @param[in] syld_el yield limit on element LGL
!> @param[in] Ckin_el,kkin_el kinematic hardening parameters on element LGL
!> @param[in] Rinf_el,biso_el isotropic hardening parameters on element LGL
!> @param[inout] dEpl_el plastic strain increment on element LGL
!> @param[inout] det_j determinant of jacobian
!> @param[inout] fx_el internal forces along x-direction on element LGL
!> @param[inout] fy_el internal forces along y-direction on element LGL

subroutine INITIAL_NL_EL(cs_nnz,cs,nm,ct,ww,dd,nn,nnt,im,ie,prop_mat, &
    displ,sxx,syy,szz,sxy,ux_el,uy_el,sxx_el,syy_el,szz_el,sxy_el,    &
    fx_el,fy_el,nl_sism,fxs_el,fys_el,sxxs_el,syys_el,szzs_el,sxys_el,&
    lambda_el,mu_el,syld_el,Ckin_el,kkin_el,Rinf_el,biso_el,          &
    riso_all,riso_el,xkin_all,xkin_el)
    
    implicit none
    integer*4, intent(in)                       :: nl_sism,im,nnt,nn,cs_nnz,nm,ie
    integer*4, intent(in),  dimension(0:cs_nnz) :: cs
    real*8, intent(in),     dimension(nm,9)     :: prop_mat
    real*8, intent(in),     dimension(nnt)      :: sxx,syy,szz,sxy,riso_all
    real*8, intent(in),     dimension(2*nnt)    :: displ
    real*8, intent(in),     dimension(4*nnt)    :: xkin_all
    real*8, intent(inout),  dimension(nn)       :: ct,ww,fx_el,fy_el,fxs_el,fys_el
    real*8, intent(inout),  dimension(nn,nn)    :: dd,ux_el,uy_el,lambda_el,mu_el
    real*8, intent(inout),  dimension(nn,nn)    :: sxx_el,syy_el,szz_el,sxy_el
    real*8, intent(inout),  dimension(nn,nn)    :: sxxs_el,syys_el,szzs_el,sxys_el
    real*8, intent(inout),  dimension(nn,nn)    :: syld_el,Riso_el
    real*8, intent(inout),  dimension(nn,nn)    :: Ckin_el,kkin_el
    real*8, intent(inout),  dimension(nn,nn)    :: Rinf_el,biso_el
    real*8, intent(inout),  dimension(4,nn,nn)  :: Xkin_el
    integer*4                                   :: i,j,is,in

    call lgl(nn,ct,ww,dd)
    
    lambda_el = prop_mat(im,2)
    mu_el     = prop_mat(im,3)
    syld_el   = prop_mat(im,5)
    Ckin_el   = prop_mat(im,6)
    kkin_el   = prop_mat(im,7)
    Rinf_el   = prop_mat(im,8)
    biso_el   = prop_mat(im,9)
    fx_el      = 0.d0
    fy_el      = 0.d0
    
    do j = 1,nn
        do i = 1,nn
            is = nn*(j -1) +i
            in = cs(cs(ie -1) + is)
            ux_el(i,j)      = displ(in)
            uy_el(i,j)      = displ(in+nnt)

            sxx_el(i,j)     = sxx(in)
            syy_el(i,j)     = syy(in)
            szz_el(i,j)     = szz(in)
            sxy_el(i,j)     = szz(in)
     
            Riso_el(i,j)    = Riso_all(in)
            Xkin_el(1,i,j)  = Xkin_all(in)     
            Xkin_el(2,i,j)  = Xkin_all(in+nnt)     
            Xkin_el(3,i,j)  = Xkin_all(in+2*nnt)     
            Xkin_el(4,i,j)  = Xkin_all(in+3*nnt)
        enddo
    enddo
    
    if (nl_sism.gt.0) then
        fxs_el=0d0
        fys_el=0d0
        sxxs_el=0d0
        syys_el=0d0
        sxys_el=0d0
        szzs_el=0d0
    endif
    
    return
end subroutine INITIAL_NL_EL
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
