module fields    
    !
    type nl_element   !< element structure for time loop (RAM saving)
        ! nl stored variables
        real*8, dimension(:,:),   allocatable    :: radius
        real*8, dimension(:,:,:), allocatable    :: stress
        real*8, dimension(:,:,:), allocatable    :: strain
        real*8, dimension(:,:,:), allocatable    :: center
        real*8, dimension(:,:,:), allocatable    :: pstrain
        ! elastic parameters
        real*8, dimension(:,:), allocatable      :: lambda,mu
        ! nonlinear parameters
        real*8, dimension(:,:),    allocatable   :: syld,rinf,biso
        real*8, dimension(:,:),    allocatable   :: ckin,kkin
    
    end type
    ! 
    type nodepatched
        ! nl stored variables
        real*8, dimension(:,:), allocatable    :: stress
        real*8, dimension(:,:), allocatable    :: strain
        real*8, dimension(:,:), allocatable    :: pstrain
    end type nodepatched 
    !
    contains
    
        !****************************************************************************
        ! ALLOCATE LOCAL VARIABLES (LINEAR CASE)
        !****************************************************************************
            
        subroutine ALLOINIT_LOC_EL(ie,nnt,cs,cs_nnz,nn,ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain, &
            fx,fy,displ,alfa1,alfa2,beta1,beta2,gamma1,gamma2)
            !
            implicit none    
            ! intent IN
            integer*4,                              intent(in)      :: ie,nnt,cs_nnz,nn
            integer*4, dimension(0:cs_nnz),         intent(in)      :: cs
            real*8, dimension(2*nnt),               intent(in)      :: displ
            real*8,    intent(in)                                   :: alfa1,beta1,gamma1
            real*8,    intent(in)                                   :: alfa2,beta2,gamma2 
            ! intent OUT 
            real*8, dimension(:),     allocatable, intent(out)      :: ct,ww
            real*8, dimension(:),     allocatable, intent(out)      :: dxdx,dydy,dxdy,dydx
            real*8, dimension(:,:),   allocatable, intent(out)      :: dd,fx,fy
            real*8, dimension(:,:,:), allocatable, intent(out)      :: dstrain
            !
            real*8,     dimension(nn,nn)                            :: duxdx,duxdy,duydx,duydy
            real*8,     dimension(nn,nn)                            :: ux,uy,sxx,syy,szz,sxy
            integer*4                                               :: i,j,is,in
            
            ! ALLOCATION
            allocate(ct(nn),ww(nn),dd(nn,nn))
            allocate(dxdx(nn),dxdy(nn),dydx(nn),dydy(nn))
            allocate(fx(nn,nn),fy(nn,nn))
            allocate(dstrain(3,nn,nn))
            ! INITIALIZATION
            call lgl(nn,ct,ww,dd)
            ! LOOP OVER GLL
            dxdx = 0.d0
            dxdy = 0.d0
            dydx = 0.d0
            dydy = 0.d0
            !
            ux = 0.0d0
            uy = 0.0d0
            !
            duxdx = 0.d0
            duxdy = 0.d0
            duydx = 0.d0
            duydy = 0.d0
            !
            sxx = 0.d0
            syy = 0.d0
            szz = 0.d0
            sxy = 0.d0
            !
            fx = 0.d0
            fy = 0.d0
            !
            dstrain(:,:,:) = 0.d0
            !
            do j = 1,nn
                do i = 1,nn
                    is = nn*(j -1) +i
                    in = cs(cs(ie -1) + is)
                    ux(i,j) = displ(in)
                    uy(i,j) = displ(in+nnt)
                enddo
            enddo
            
            ! MAKE DERIVATIVES
            call MAKE_DERIVATIVES(alfa1,alfa2,beta1,beta2,gamma1,gamma2,nn,ct,&
                dxdx,dxdy,dydx,dydy)
             
            ! MAKE STRAIN
            call MAKE_STRAIN(nn,dd,dxdx,dxdy,dydx,dydy,ux,uy,duxdx,duxdy,duydx,duydy)
            do j = 1,nn
                do i = 1,nn
                    dstrain(1,i,j) = duxdx(i,j)
                    dstrain(2,i,j) = duydy(i,j)
                    dstrain(3,i,j) = duxdy(i,j)+duydx(i,j)
                enddo
            enddo

            return
        end subroutine ALLOINIT_LOC_EL
        
        !****************************************************************************
        ! ALLOCATE LOCAL VARIABLES (NONLINEAR CASE)
        !****************************************************************************
        
        subroutine ALLOINIT_LOC_NL(ie,nnt,cs,cs_nnz,nn,ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain,dstrial, &
            fx,fy,displ,alfa1,alfa2,beta1,beta2,gamma1,gamma2,lambda,mu)
            !
            implicit none    
            ! intent IN
            integer*4,                              intent(in)      :: ie,nnt,cs_nnz,nn
            integer*4, dimension(0:cs_nnz),         intent(in)      :: cs
            real*8, dimension(2*nnt),               intent(in)      :: displ
            real*8, dimension(nn,nn),               intent(in)      :: lambda,mu
            real*8,    intent(in)                                   :: alfa1,beta1,gamma1
            real*8,    intent(in)                                   :: alfa2,beta2,gamma2 
            ! intent OUT 
            real*8, dimension(:),     allocatable, intent(out)      :: ct,ww
            real*8, dimension(:),     allocatable, intent(out)      :: dxdx,dydy,dxdy,dydx
            real*8, dimension(:,:),   allocatable, intent(out)      :: dd,fx,fy
            real*8, dimension(:,:,:), allocatable, intent(out)      :: dstrain
            real*8, dimension(:,:,:), allocatable, intent(out)      :: dstrial
            !
            real*8,     dimension(nn,nn)                            :: duxdx,duxdy,duydx,duydy
            real*8,     dimension(nn,nn)                            :: ux,uy,sxx,syy,szz,sxy
            integer*4                                               :: i,j,is,in
            
            ! ALLOCATION
            allocate(ct(nn),ww(nn),dd(nn,nn))
            allocate(dxdx(nn),dxdy(nn),dydx(nn),dydy(nn))
            allocate(fx(nn,nn),fy(nn,nn))
            allocate(dstrain(3,nn,nn))
            allocate(dstrial(4,nn,nn))
            ! INITIALIZATION
            call LGL(nn,ct,ww,dd)
            ! LOOP OVER GLL
            dxdx = 0.0d0
            dxdy = 0.0d0
            dydx = 0.0d0
            dydy = 0.0d0
            !
            ux = 0.0d0
            uy = 0.0d0
            !
            duxdx = 0.0d0
            duxdy = 0.0d0
            duydx = 0.0d0
            duydy = 0.0d0
            !
            sxx = 0.0d0
            syy = 0.0d0
            szz = 0.0d0
            sxy = 0.0d0
            !
            fx = 0.0d0
            fy = 0.0d0
            !
            dstrain(:,:,:) = 0.0d0
            dstrial(:,:,:) = 0.0d0
            !
            do j = 1,nn
                do i = 1,nn
                    is = nn*(j -1) +i
                    in = cs(cs(ie -1) + is)
                    ux(i,j) = displ(in)
                    uy(i,j) = displ(in+nnt)
                enddo
            enddo
            
            ! MAKE DERIVATIVES
            call MAKE_DERIVATIVES(alfa1,alfa2,beta1,beta2,gamma1,gamma2,nn,ct,&
                dxdx,dxdy,dydx,dydy)
             
            ! MAKE STRAIN
            call MAKE_STRAIN(nn,dd,dxdx,dxdy,dydx,dydy,ux,uy,duxdx,duxdy,duydx,duydy)
            do j = 1,nn
                do i = 1,nn
                    dstrain(1,i,j) = duxdx(i,j)
                    dstrain(2,i,j) = duydy(i,j)
                    dstrain(3,i,j) = duxdy(i,j)+duydx(i,j)
                enddo
            enddo

            ! MAKE STRESS
            call MAKE_STRESS(nn,lambda,mu,duxdx,duxdy,duydx,duydy,sxx,syy,szz,sxy)
            do j = 1,nn
                do i = 1,nn
                    dstrial(1,i,j) = sxx(i,j)
                    dstrial(2,i,j) = syy(i,j)
                    dstrial(3,i,j) = szz(i,j)
                    dstrial(4,i,j) = sxy(i,j)
                enddo
            enddo
            !
            return
            !
        end subroutine ALLOINIT_LOC_NL
        
        !****************************************************************************
        ! DEALLOCATE LOCAL VARIABLES
        !****************************************************************************
        
        subroutine DEALLOCATE_LOC(ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain,dstrial,fx,fy)
            !
            implicit none    
            ! intent INOUT 
            real*8,     dimension(:),  allocatable, intent(inout) :: ct,ww
            real*8,     dimension(:),  allocatable, intent(inout) :: dxdx,dydy,dxdy,dydx
            real*8,     dimension(:,:),allocatable, intent(inout) :: dd,fx,fy
            real*8,     dimension(:,:,:), allocatable, intent(inout) :: dstrain
            real*8,     dimension(:,:,:), allocatable, intent(inout) :: dstrial
            ! DEALLOCATION
            deallocate(ct)
            deallocate(ww)
            deallocate(dd)
            deallocate(dxdx)
            deallocate(dxdy)
            deallocate(dydx)
            deallocate(dydy)
            deallocate(fx)
            deallocate(fy)
            deallocate(dstrain)
            deallocate(dstrial)
            !
            return
        end subroutine DEALLOCATE_LOC

end module fields
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
