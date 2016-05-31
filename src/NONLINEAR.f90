module nonlinear2d
    !
    use fields
    !
    implicit none
    !
    real*8, parameter                   :: zero=0.d0,one=1.0d0
    real*8, parameter                   :: half=0.5d0,two=2.0d0,three=3.0d0
    !
    real*8, parameter :: FTOL = 0.000001D0
    real*8, parameter :: LTOL = 0.000001D0
    real*8, parameter :: STOL = 1.0D0
    !
    real*8, parameter, dimension(4,4)   :: MM = reshape((/ &
        one , zero, zero, zero, &
        zero, one , zero, zero, &
        zero, zero, one , zero, &
        zero, zero, zero, two   &
        /), (/4,4/))
    
    real*8, parameter, dimension(4,4)   :: MM1 = reshape((/ &   
        one , zero, zero, zero, &
        zero, one , zero, zero, &
        zero, zero, one , zero, &
        zero, zero, zero, half   &
        /), (/4,4/))
    
    real*8, parameter, dimension(4)     :: m = (/one,one,one,zero/)

    contains
        
        !****************************************************************************
        ! MAKE NL INTERNAL FORCES
        !****************************************************************************
        
        subroutine MAKE_INTERNAL_FORCES_NL(nnt,ne,nm,cs_nnz,cs,sdeg_mat,snl,&
            alfa1,alfa2,beta1,beta2,gamma1,gamma2,displ,fk,mvec,dt)
            !
            implicit none
            ! intent IN
            integer*4, intent(in)                       :: nnt,ne,nm,cs_nnz
            integer*4, intent(in), dimension(nm)        :: sdeg_mat
            integer*4, intent(in), dimension(0:cs_nnz)  :: cs
            real*8,    intent(in), dimension(ne)        :: alfa1,beta1,gamma1
            real*8,    intent(in), dimension(ne)        :: alfa2,beta2,gamma2
            real*8,    intent(in), dimension(2*nnt)     :: displ,mvec
            real*8,    intent(in)                       :: dt
            ! intent INOUT
            real*8,    intent(inout), dimension(2*nnt)  :: fk
            type(nl_element), intent(inout), dimension(ne) :: snl 
            ! 
            real*8,     dimension(:),  allocatable      :: ct,ww
            real*8,     dimension(:),  allocatable      :: dxdx,dydy,dxdy,dydx
            real*8,     dimension(:,:),allocatable      :: dd,fx,fy
            real*8,     dimension(:,:,:), allocatable   :: dstrain,dstrial
            !
            real*8, dimension(4)                        :: stress_,dstrial_,center_
            real*8, dimension(4)                        :: dstrain_,dpstrain_ 
            real*8                                      :: radius_,syld_ 
            real*8                                      :: lambda_,mu_     
            real*8                                      :: ckin_,kkin_  
            real*8                                      :: biso_,rinf_   
            !
            real*8                                      :: alpha_epl,t1ux,t1uy,t2ux
            real*8                                      :: t2uy,t1fx,t1fy,t2fx,t2fy
            integer*4                                   :: ie,ip,iq,il,im,nn,is,in
            logical                                     :: st_epl
            !
            do ie = 1,ne
                im = cs(cs(ie-1) + 0)
                nn = sdeg_mat(im)+1
                !*********************************************************************************
                ! COMPUTE STRAIN 
                !*********************************************************************************
                
                call ALLOINIT_LOC_NL(ie,nnt,cs,cs_nnz,nn,ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain,dstrial, &
                    fx,fy,displ,alfa1(ie),alfa2(ie),beta1(ie),beta2(ie),gamma1(ie),gamma2(ie),&
                    snl(ie)%lambda,snl(ie)%mu)
!                    dstrain = dstrain*dt
!                    dstrial = dstrial*dt
!                dstrain(:,:,:) = dstrain(:,:,:) - snl(ie)%strain(:,:,:)
!                dstrial(:,:,:) = dstrial(:,:,:) - snl(ie)%stress(:,:,:)!
                snl(ie)%strain(:,:,:) = dstrain(:,:,:)
                snl(ie)%stress(:,:,:) = dstrial(:,:,:) 
                !*********************************************************************************
                ! COMPUTE STRESS
                !*********************************************************************************
                
!                do iq = 1,nn
!                    do ip = 1,nn
!                        ! STARTING POINT
!                        stress_ = snl(ie)%stress(:,ip,iq)
!                        center_ = snl(ie)%center(:,ip,iq)
!                        radius_ = snl(ie)%radius(ip,iq)
!                        lambda_ = snl(ie)%lambda(ip,iq)
!                        mu_     = snl(ie)%mu(ip,iq)
!                        syld_   = snl(ie)%syld(ip,iq)
!                        ckin_   = snl(ie)%ckin(ip,iq)
!                        kkin_   = snl(ie)%kkin(ip,iq)
!                        biso_   = snl(ie)%biso(ip,iq)
!                        rinf_   = snl(ie)%rinf(ip,iq)
!                        
!
!                        ! STRAIN INCREMENT
!                        dstrain_(:)   = zero
!                        dstrain_(1:2) = dstrain(1:2,ip,iq)
!                        dstrain_(4)   = dstrain(3,ip,iq)
!                        dpstrain_(:)  = zero
!                        ! TRIAL STRESS INCREMENT
!                        dstrial_(:)   = zero
!                        dstrial_(:)   = dstrial(:,ip,iq)
!                        
!                        ! CHECK PLASTICITY
!                        call check_plasticity(dstrial_,stress_,center_,radius_,&
!                            syld_,st_epl,alpha_epl)
!                        ! PLASTIC CORRECTION
!                        if (st_epl) then
!                            write(*,*) "alpha_epl",alpha_epl
!                            dstrain_(:) = (one-alpha_epl)*dstrain_(:)
!                            call plastic_corrector(dstrain_,dstrial_,center_,radius_,syld_,&
!                                biso_,rinf_,ckin_,kkin_,mu_,lambda_,dpstrain_)
!                        endif
!                        ! STRESS VECTOR
!                        snl(ie)%stress(:,ip,iq) = dstrial_(:)
!                        ! CENTER
!                        snl(ie)%center(:,ip,iq) = center_(:)
!                        ! RADIUS
!                        snl(ie)%radius(ip,iq)   = radius_
!                        ! PLASTIC STRAIN VECTOR
!                        snl(ie)%plastic_strain(:,ip,iq) = snl(ie)%plastic_strain(:,ip,iq) + &
!                            dpstrain_(:)
!                        
!                    enddo
!                enddo
                
                !*********************************************************************************
                ! COMPUTE LOCAL INTERNAL FORCES
                !*********************************************************************************
                call MAKE_INTERNAL_FORCE(nn,ww,dd,dxdx,dxdy,dydx,dydy,snl(ie)%stress(1,:,:),&
                    snl(ie)%stress(2,:,:),snl(ie)%stress(4,:,:),fx,fy)
                !*********************************************************************************
                ! COMPUTE LOCAL INTERNAL FORCES
                !*********************************************************************************
                do iq = 1,nn
                    do ip = 1,nn
                        is = nn*(iq-1)+ip
                        in = cs(cs(ie-1)+is)
                        fk(in)     = fk(in)     + fx(ip,iq)
                        fk(in+nnt) = fk(in+nnt) + fy(ip,iq)
                    enddo
                enddo
                ! DEALLOCATE ELEMENT-WISE VARIABLES
                call DEALLOCATE_LOC(ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain,dstrial,fx,fy)
            enddo

            return
            !
        end subroutine MAKE_INTERNAL_FORCES_NL
        
        !*********************************************************************************
        ! STIFFNESS MATRIX 
        !*********************************************************************************
        
        subroutine STIFF_MATRIX(lambda,mu,DEL)
            implicit none
            ! intent IN
            real*8, intent(in) :: lambda,mu
            ! intent OUT
            real*8, intent(out), dimension(4,4) :: DEL
            !
            DEL(:,:) = 0.d0
            DEL(1,1) = lambda+2*mu
            DEL(2,2) = lambda+2*mu
            DEL(3,3) = lambda+2*mu
            DEL(4,4) = mu
            DEL(1,2) = lambda
            DEL(1,3) = lambda
            DEL(2,1) = lambda
            DEL(2,3) = lambda
            DEL(3,1) = lambda
            DEL(3,2) = lambda
            !
            return
        end subroutine STIFF_MATRIX
        
        !****************************************************************************
        ! MISES YIELDING LOCUS AND GRADIENT
        !****************************************************************************

        subroutine mises_yld_locus(stress, center, radius, syld, FM, gradFM)
            ! 
            implicit none
            ! intent IN
            real*8,               intent(in)    :: radius,syld
            real*8, dimension(4), intent(in)    :: stress,center
            ! intent OUT
            real*8, dimension(4), intent(out)   :: gradFM
            real*8,               intent(out)   :: FM
            !
            real*8, dimension(4)                :: dev
            real*8                              :: tau_eq
            integer*4                           :: k
            ! COMPUTE TENSOR COMPONENTS 
            call tensor_components (stress, dev)
            ! COMPUTE MISES FUNCTION
            call tau_mises(dev-center,tau_eq)
            FM = tau_eq - syld - radius
            ! COMPUTE MISES FUNCTION GRADIENT
            gradFM(1) = three/two*(dev(1)-center(1))/tau_eq
            gradFM(2) = three/two*(dev(2)-center(2))/tau_eq
            gradFM(3) = 1.5*(dev(3)-center(3))/tau_eq
!            gradFM(3) = zero
            gradFM(4) = three*(dev(4)-center(4))/tau_eq
            !  
            return
            !
        end subroutine mises_yld_locus
        
        !****************************************************************************
        ! OCTAHEDRAL SHEAR STRESS
        !****************************************************************************
        
        subroutine tau_mises(dev,J2M)
            ! correspondence rule:
            !  2-d case     3-d case  
            ! (11) <=> sxx   (11) <=> sxx 
            ! (22) <=> syy   (22) <=> syy 
            ! (33) <=> szz   (33) <=> szz 
            ! (12) <=> sxy   (12) <=> sxy
            !                (23) <=> syz
            !                (31) <=> szx
            ! and continuum mechanics sign conventions
            !
            implicit none
            ! intent IN
            real*8, dimension(4), intent(in)    :: dev
            ! intent OUT
            real*8,               intent(out)   :: J2M
            !
            real*8, dimension(4)                :: temp
            real*8                              :: J2M2
            
            ! COMPUTE OCTAHEDRAL SHEAR STRESS 
            temp = three/two*matmul(MM,dev)
            J2M2 = dot_product(dev,temp)
            J2M  = sqrt(J2M2)
            !
            return
            ! 
        end subroutine
        
        !****************************************************************************
        ! TENSOR COMPONENTS (SPHERICAL & DEVIATORIC)
        !****************************************************************************
        
        subroutine tensor_components(stress, dev)
            !     
            implicit none
            ! intent IN
            real*8, dimension(4), intent(in)    :: stress
            ! intent OUT
            real*8, dimension(4), intent(out)   :: dev
            !
            real*8                              :: I1 
            !
            I1 = dot_product(stress,m)/three
            !
            dev = stress - I1*m
            !
            return
            !
        end subroutine tensor_components
            
        !****************************************************************************
        ! CHECK PLASTIC CONSISTENCY (KKT CONDITIONS)
        !****************************************************************************
        
        subroutine check_plasticity (dtrial, stress0, center, radius, syld, &
            st_elp, alpha_epl)
            !     
            implicit none
            ! intent IN
            real*8,                 intent(in)      :: radius,syld
            real*8, dimension(4),   intent(in)      :: center,stress0
            ! intent INOUT
            real*8, dimension(4),   intent(inout)   :: dtrial
            ! intent OUT
            real*8,                 intent(out)     :: alpha_epl
            logical,                intent(out)     :: st_elp
            !
            integer*4                               :: k
            real*8                                  :: FS,FT,checkload
            real*8, dimension(4)                    :: gradFS, gradFT, stress1
            integer*4                               :: plasticity
            !
            ! PREDICTION STRESS
            stress1= dtrial+stress0
            ! 
            ! CHECK MISES FUNCTION
            call mises_yld_locus(stress0,center,radius,syld,FS,gradFS)
            call mises_yld_locus(stress1,center,radius,syld,FT,gradFT)
            
            ! CHECK LOAD DIRECTION 
            checkload = dot_product(gradFS,dtrial)/sum(gradFS**2)/sum(dtrial**2)
            
            ! CHECK PLASTICITiY 
            if (abs(FS).le.FTOL) then ! FS = 0
                !
                if (checkload.ge.-LTOL) then ! PLASTIC LOADING
                    alpha_epl = zero
                    st_elp    = .true.
                    plasticity = 1
                else ! GENERALIZED UNLOADING
                    
                    if (FT.lt.-FTOL) then ! ELASTIC UNLOADING
                        alpha_epl = one
                        st_elp    = .false.
                        plasticity = 2
                    
                    elseif(FT.gt.FTOL) then ! PLASTIC UNLOADING
                        call gotoFpegasus(stress0,dtrial,center,radius,syld,10,alpha_epl)
                        st_elp    = .true.
                        plasticity=3
                    endif
                
                end if

            elseif (FS.lt.-FTOL) then ! FS<0

                if (FT.le.FTOL) then ! ELASTIC LOADING
                    alpha_epl = one
                    st_elp    = .false.
                    plasticity=4
                else ! ELASTO-PLASTIC LOADING
                    call gotoFpegasus(stress0,dtrial,center,radius,syld,1,alpha_epl)
                    st_elp    = .true.
                    plasticity=5
                end if

            elseif (FS.gt.FTOL) then
                write(*,*) "ERROR FS:",FS,">FTOL"
                alpha_epl  = zero
                st_elp = .true.
                plasticity=6
            end if
            ! ON-LOCUS STRESS STATE 
            dtrial=stress0+dtrial*alpha_epl
            call mises_yld_locus(dtrial,center,radius,syld,FS,gradFS)
            write(*,*) "plasticity",plasticity
        end subroutine check_plasticity
        
        subroutine plastic_corrector(dEps_alpha,stress,center,syld, &
            radius,biso,Rinf,Ckin,kkin,mu,lambda,dEpl)
            !
            implicit none
            !
            real*8, dimension(4), intent(inout) :: dEps_alpha,stress,center,dEpl 
            real*8,               intent(inout) :: radius 
            real*8,               intent(in)    :: syld,biso,Rinf,Ckin,kkin,mu,lambda
            real*8, dimension(4,4)              :: DEL
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,0.5/)
            real*8                              :: Ttot,deltaTk,qq,R1,R2,dR1,dR2,err0,err1,err2,err3
            real*8                              :: hard0,hard1,hard2,deltaTmin
            real(8)                             :: FM,Resk
            logical                             :: flag_fail
            real*8, dimension(4)                :: gradFM,S1,S2,X1,X2
            real*8, dimension(4)                :: dS1,dS2,dX1,dX2,dEpl1,dEpl2
            
            call stiff_matrix(lambda,mu,DEL)
            deltaTk = one
            Ttot    = zero
            deltaTmin = 0.001d0
            flag_fail =.true.
            do while (Ttot.lt.one)
                write(*,*) "BEFORE------"
                write(*,*) "flag_fail",flag_fail
                write(*,*) "deltaTk",deltaTk
                write(*,*) "Ttot",Ttot
                Resk  = zero
                dS1   = zero
                dX1   = zero
                dS2   = zero
                dX2   = zero
                dR1   = zero
                dR2   = zero
                dEpl1 = zero
                dEpl2 = zero
                
                ! FIRST ORDER COMPUTATION
                
                call ep_integration(dEps_alpha*deltaTk,stress,center,radius,syld,&
                    mu,lambda,biso,Rinf,Ckin,kkin,dS1,dX1,dR1,dEpl1,hard1)

                S1 = stress + dS1
                X1 = center + dX1 
                R1 = radius + dR1
                
                ! SECOND ORDER COMPUTATION
                call ep_integration(dEps_alpha*deltaTk,S1,X1,R1,syld,mu,lambda,&
                    biso,Rinf,Ckin,kkin,dS2,dX2,dR2,dEpl2,hard2)
                
                ! TEMPORARY VARIABLES
                S1 = stress + half*(dS1+dS2)
                X1 = center + half*(dX1+dX2)
                R1 = radius + half*(dR1+dR2)
                dEpl1 =       half*(dEpl1+dEpl2)
                
                ! ERROR
                call tau_mises(dS2-dS1,err0)
                call tau_mises(S1,err1)
                
                call tau_mises(dX1-dX2,err2)
                call tau_mises(X1,err3)
                Resk = half*max(epsilon(Resk),err0/err1,err2/err3)
                
                if (Resk.le.STOL) then
                    write(*,*) "RESK",Resk,"< STOL",STOL 
                    stress = S1
                    center = X1
                    radius = R1
                    dEpl   = dEpl+dEpl1

                    call mises_yld_locus (stress, center,radius,syld,FM,gradFM)
                    write(*,*) "before drift - FM:",FM
                    if (FM.gt.FTOL) then
                        call drift_corr(stress,center,radius,syld,&
                                biso,Rinf,Ckin,kkin,lambda,mu,dEpl)
                    endif
                    qq = min(0.9d0*sqrt(STOL/Resk),1.1d0)
                    write(*,*) "qq",qq
                    if (flag_fail) then
                        qq = min(qq,one)
                        write(*,*) "qq",qq
                    endif
                    flag_fail=.false.
                    Ttot=Ttot+deltaTk
                    deltaTk=qq*deltaTk
                    write(*,*) "deltaT",deltaTk
                    deltaTk=max(qq*deltaTk,deltaTmin)
                    write(*,*) "deltaT",deltaTk
                    deltaTk=min(deltaTk,one-Ttot)
                    write(*,*) "deltaT",deltaTk
                else
                    qq=max(0.9d0*sqrt(STOL/Resk),0.1d0)
                    deltaTk=max(qq*deltaTk,deltaTmin)
                    flag_fail=.true.

                end if
                write(*,*) "BEFORE------"
                write(*,*) "flag_fail",flag_fail
                write(*,*) "deltaTk",deltaTk
                write(*,*) "Ttot",Ttot
            end do
        end subroutine plastic_corrector

        subroutine ep_integration(dStrain,Stress,center,radius,syld,mu,lambda,biso,Rinf,&
            Ckin,kapakin,dStress,dcenter,dradius,dEplast,hard)
            !
            implicit none
            ! intent IN
            real*8              , intent(in)    :: radius,syld,mu,lambda,biso,Rinf,Ckin,kapakin
            real*8, dimension(4), intent(in)    :: dStrain,Stress,center
            ! intent INOUT
            real*8, dimension(4), intent(inout) :: dStress,dcenter,dEplast
            real*8,               intent(inout) :: dradius
            ! intent OUT
            real*8,               intent(out)   :: hard
            !
            real*8, dimension(4)                :: gradF
            !real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)              :: DEL
            real*8                              :: Fmises,dPlast
            integer*4                           :: j,k
            
            ! PREDICTION
            call mises_yld_locus (Stress,center,radius,syld,Fmises,gradF)
            call stiff_matrix(lambda,mu,DEL)

            ! PLASTIC MULTIPLIER
            call compute_plastic_modulus(dStrain,Stress,center,radius,mu,lambda,syld, &
                biso,Rinf,Ckin,kapakin,dPlast,hard)
            
            ! INCREMENTS
            call hardening_increments(Stress,radius,center,syld, &
                biso,Rinf,Ckin,kapakin,dradius,dcenter)
            
            dradius = dPlast*dradius
            dcenter = dPlast*dcenter
!            do k=1,4 
!                dEplast(k)=dEplast(k)+dPlast*gradF(k)*A(k)
!            end do
            dEplast = dPlast*matmul(MM1,gradF)
            dstress = matmul(DEL,dstrain-dEplast) 
!            do j = 1,4 ! stress increment
!                do k = 1,4
!                    dstress(j)=DEL(k,j)*(dstrain(j)-A(k)*dPlast*gradF(k))  
!                end do
!            end do
            
            return

        end subroutine ep_integration


        subroutine compute_plastic_modulus(dEps,stress,center,radius,mu,lambda, &
            syld,biso,Rinf,Ckin,kkin,dPlastMult,hard)
            !
            implicit none
            ! intent IN
            real*8,                 intent(in) :: mu,lambda,syld   
            real*8,                 intent(in) :: radius,biso,Rinf,Ckin,kkin
            real*8, dimension(4),   intent(in) :: dEps,stress,center
            ! intent OUT
            real*8,                 intent(out):: dPlastMult,hard
            !
            real*8                             :: Ah,FM
            real*8, dimension(4)               :: gradF,tempv
            real*8, dimension(4,4)             :: DEL
!            integer*4                          :: j,k
            
            call stiff_matrix(lambda,mu,DEL)
            call mises_yld_locus(stress,center,radius,syld,FM,gradF)
            
            hard = Ckin+biso*(Rinf-radius)
            hard = hard-kkin*dot_product(center,gradF)
            
            dPlastMult = zero 
            
            tempv = matmul(MM1,gradF)
            tempv = matmul(DEL,tempv)
            Ah    = dot_product(tempv,gradF)
            tempv = matmul(DEL,dEps)
            dPlastMult = dot_product(gradF,tempv)
!            Ah=zero
!            do k = 1,4
!                do j = 1,4
!                    Ah   = Ah+gradF(j)*DEL(j,k)*gradF(k)*A(k)
!                    dPlastMult = dPlastMult+gradF(j)*DEL(j,k)*dEps(k)
!                end do
!            end do
            dPlastMult = max(0.0d0,dPlastMult/(hard+Ah))
        end subroutine compute_plastic_modulus


        subroutine hardening_increments(stress, radius, center, syld, &
            biso, rinf, ckin, kkin, dradius, dcenter)

            ! INCREMENTS OF INTRINSIC STATIC VARIABLES

            real*8,               intent(in) :: syld,radius
            real*8,               intent(in) :: biso,rinf,ckin,kkin
            real*8, dimension(4), intent(in) :: stress,center
            !
            real*8,               intent(out):: dradius
            real*8, dimension(4), intent(out):: dcenter
            
            real*8                           :: FM
            real*8, dimension(4)             :: gradFM
            !real*8, dimension(4), parameter  :: A = (/1.0,1.0,1.0,0.5/)
            !integer*4                        :: k

            ! INCREMENT IN ISOTROPIC HARDENING VARIABLES (R)
            dradius = biso*(rinf-radius)

            ! INCREMENT IN KINEMATIC HARDENING VARIABLES (Xij)
            call mises_yld_locus (stress,center,radius,syld,FM,gradFM)
            dcenter = (ckin*two/three)*matmul(MM1,gradFM)-center*kkin
            !dX_ij(:)=2*A(:)*gradF_mises(:)*C_lmc/3-X_ij(:)*kapa_lmc

        end subroutine hardening_increments

        subroutine drift_corr(stress,center,radius,syld, &
            biso,Rinf,Ckin,kkin,lambda,mu,dEplastic)

            ! DRIFT CORRECTION (RADIAL RETURN)
            real*8, dimension(4), intent(inout) :: stress,center,dEplastic
            real*8,               intent(inout) :: radius
            real*8,               intent(in)    :: lambda,mu,syld,biso,Rinf,Ckin,kkin
            real*8                              :: F0,F1,beta,hard,radiust
            real*8, dimension(4)                :: tempv,gradF0,gradF1,dstress,stresst,centert
            !real*8, dimension(4),     parameter :: A = (/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)              :: DEL
            integer*4                           :: counter,k,j
            real*8, parameter :: FTOL_DRIFT =   0.0000001D0 
            ! INITIAL PLASTIC CONDITION
            call mises_yld_locus(stress,center,radius,syld,F0,gradF0)
            call stiff_matrix(lambda,mu,DEL)
            do counter=1,5 
                ! COMPUTE HARDENING INCREMENTS
                hard = biso*(Rinf-radius)
                hard = hard + Ckin
                hard = hard - kkin*dot_product(gradF0,center)
                ! COMPUTE BETA FOR DRIFT CORRECTION
                beta  = zero
                tempv = matmul(MM1,gradF0)
                tempv = matmul(DEL,tempv)
                beta  = dot_product(tempv,gradF0)
!                do j=1,4
!                    do k=1,4
!                        beta=beta+gradF0(k)*DEL(k,j)*A(j)*gradF0(j)
!                    end do
!                end do
                beta = F0/(hard+beta)
                ! STRESS-STRAIN-HARDENING CORRECTION
                dstress = zero
                dstress = -beta*tempv
                
!                do k=1,4
!                    do j=1,4
!                        dstress(k)=dstress(k)-beta*DEL(j,k)*A(j)*gradF0(j)
!                    end do
!                end do
                stresst = stress+dstress
                centert = center+beta*((two*ckin/three)*matmul(MM,gradF0)-center*kkin)
                radiust = radius+beta*(Rinf-radius)*biso
                
                ! CHECK DRIFT
                call mises_yld_locus(stresst,centert,radiust,syld,F1,gradF1)
                if (abs(F1).gt.abs(F0)) then
                    beta    = F0/dot_product(gradF0,gradF0)
                    dstress = -beta*gradF0
                    stresst = stress+dstress
                    centert = center
                    radiust = radius
                    call mises_yld_locus(stresst,centert,radiust,syld,F1,gradF1)
                endif
                stress = stresst
                center = centert
                radius = radiust
                dEplastic = dEplastic+beta*matmul(MM1,gradF0)

                if (abs(F1).le.FTOL_DRIFT) then
                    exit
                else
                    F0     = F1
                    gradF0 = gradF1
                endif
            enddo
            if (abs(F1).le.FTOL) then
                write(*,*) "drift corrected!",F1
            else 
                write(*,*) "DRIFT NOT CORRECTED"
                read(*,*)
            endif 
            return
        end subroutine drift_corr

        subroutine gotoFpegasus(start0,dtrial,center,radius,s0,nsub,alpha)
            implicit none
            ! intent IN
            real*8, dimension(4), intent(in)    :: start0,dtrial,center
            real*8,               intent(in)    :: radius,s0
            integer*4,            intent(in)    :: nsub
            ! intent OUT
            real*8,               intent(out)   :: alpha
            real*8, dimension(4)                :: stress0,stress1,stress,gradF
            real*8                              :: dalpha,alpha0,alpha1,F0,F1,FM,Fsave
            integer*4                           :: counter0,counter1
            logical                             :: flagxit

            alpha1  = one
            alpha0  = zero
            stress0 = start0+alpha0*dtrial
            stress1 = start0+alpha1*dtrial
            call mises_yld_locus(stress0,center,radius,s0,F0,gradF)
            call mises_yld_locus(stress1,center,radius,s0,F1,gradF)
            if (nsub.gt.1)then
                Fsave=F0
                do counter0=1,3
                    dalpha = (alpha1-alpha0)/nsub
                    flagxit=.false.
                    do counter1=1,nsub
                        alpha=alpha0+dalpha
                        stress=start0+alpha*dtrial
                        call mises_yld_locus(stress,center,radius,s0,FM,gradF)
                        if (FM.gt.FTOL) then
                            alpha1=alpha
                            if (F0.lt.-FTOL) then
                                F1=FM
                                flagxit=.true.
                            else
                                alpha0=zero
                                F0=Fsave
                            endif
                            exit ! exit loop counter1=1,nsub
                        else
                            alpha0=alpha
                            F0=FM
                        endif
                    end do
                    if (flagxit) then
                        exit ! exit loop counter0=1,3
                    endif
                end do
            end if
            do counter0=1,10
                alpha  = alpha1 - F1*(alpha1 - alpha0)/(F1-F0)
                stress = start0 + alpha*dtrial
                call mises_yld_locus(stress,center,radius,s0,FM,gradF)
                if (abs(FM).le.FTOL) then
                    exit
                else
                    if(FM*F1.lt.zero) then
                        alpha0=alpha1
                        F0=F1
                    else
                        F0=F0*half
                    endif
                    F1=FM
                    alpha1=alpha
                endif
            end do

            if (FM.gt.FTOL) then
                write(*,*) "WARNING: F=",FM,">FTOL!!!!!!"
            else
                write(*,*) "PEGASUS: INTERSECTION FOUND",FM
            endif
        end subroutine gotoFpegasus
        
end module nonlinear2d
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
