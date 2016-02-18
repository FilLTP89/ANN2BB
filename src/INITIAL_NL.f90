subroutine INITIAL_NL(cs_nnz,cs,nm,nn,nnt,im,prop_mat,Sxx_el,Syy_el,Sxy_el,Szz_el  &
    lambda_el,mu_el,Syld_el,Ckin_el,kkin_el,Riso_el,Rinf_el,biso_el,Xkin_el,    &
    Stress_all,Xkin_all,Riso_all,fx_el,fy_el)
    
    implicit none
    integer*4, intent(in)                     :: im,nnt,nn,cs_nnz,nm
    integer*4, intent(in), dimension(0:cs_nnz):: cs
    real*8, intent(in), dimension(nm,8)       :: prop_mat
    real*8, intent(in), dimension(nnt)        :: Riso_all
    real*8, intent(inout), dimension(nn)      :: fx_el,fy_el
    real*8, intent(in), dimension(4*nnt)      :: Stress_all,Xkin_all
    real*8, intent(inout), dimension(nn,nn)   :: Sxx_el,Syy_el,Szz_el,Sxy_el
    real*8, intent(inout), dimension(nn,nn)   :: lamdbda_el,mu_el,syld_el,Riso_el
    real*8, intent(inout), dimension(nn,nn)   :: Ckin_el,kkin_el,Rinf_el,biso_el
    real*8, intent(inout), dimension(4,nn,nn) :: Xkin_el
    integer*4,                                :: i,j,is,in

    lambda_el = prop_mat(im,2)
    mu_el     = prop_mat(im,3)
    syld_el   = prop_mat(im,4)
    Ckin_el   = prop_mat(im,5)
    kkin_el   = prop_mat(im,6)
    Rinf_el   = prop_mat(im,7)
    biso_el   = prop_mat(im,8)

    fk_el      = 0.d0
    fx_el      = 0.d0
    fy_el      = 0.d0
    
    do j = 1,nn
        do i = 1,nn
            is = nn*(j -1) +i
            in = cs(cs(ie -1) + is)
            Riso_el(i,j)    = Riso_all(in)
            sxx_el(i,j)     = Stress_all(in)
            syy_el(i,j)     = Stress_all(in+nnt)
            sxy_el(i,j)     = Stress_all(in+2*nnt)
            szz_el(i,j)     = Stress_all(in+3*nnt)
            Xkin_el(1,i,j)  = Xkin_all(in)     
            Xkin_el(2,i,j)  = Xkin_all(in+nnt)     
            Xkin_el(3,i,j)  = Xkin_all(in+2*nnt)     
            Xkin_el(4,i,j)  = Xkin_all(in+3*nnt)
        enddo
    enddo

    return
end subroutine INITIAL_NL
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

