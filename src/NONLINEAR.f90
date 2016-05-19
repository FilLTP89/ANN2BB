module nonlinear2d
    !
    use fields
    !
    implicit none
    !
    real*8, parameter :: FTOL = 0.0010000000000d0
    real*8, parameter :: LTOL = 0.0010000000000d0
    real*8, parameter :: STOL = 0.0010000000000d0
    !
    contains
        
        !****************************************************************************
        ! MAKE NL INTERNAL FORCES
        !****************************************************************************
        
        subroutine MAKE_INTERNAL_FORCES_NL(nnt,ne,nm,cs_nnz,cs,sdeg_mat,snl,&
            alfa1,alfa2,beta1,beta2,gamma1,gamma2,displ,fk,mvec)
            !
            implicit none
            ! intent IN
            integer*4, intent(in)                       :: nnt,ne,nm,cs_nnz
            integer*4, intent(in), dimension(nm)        :: sdeg_mat
            integer*4, intent(in), dimension(0:cs_nnz)  :: cs
            real*8,    intent(in), dimension(ne)        :: alfa1,beta1,gamma1
            real*8,    intent(in), dimension(ne)        :: alfa2,beta2,gamma2
            real*8,    intent(in), dimension(2*nnt)     :: displ,mvec
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
            real*8                                      :: FS,alpha_epl,t1ux,t1uy,t2ux
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
                
                dstrain(:,:,:) = dstrain(:,:,:) - snl(ie)%strain(:,:,:)
                dstrial(:,:,:) = dstrial(:,:,:) - snl(ie)%stress(:,:,:)
                snl(ie)%strain(:,:,:) = snl(ie)%strain(:,:,:) + dstrain(:,:,:) 
                !*********************************************************************************
                ! COMPUTE STRESS
                !*********************************************************************************
                
                do iq = 1,nn
                    do ip = 1,nn
                        write(*,*) "=================================="
                        write(*,*) "STRESS N: ",SNL(IE)%STRESS(:,IP,IQ)
                        ! STARTING POINT
                        stress_ = snl(ie)%stress(:,ip,iq)
                        center_ = snl(ie)%center(:,ip,iq)
                        radius_ = snl(ie)%radius(ip,iq)
                        lambda_ = snl(ie)%lambda(ip,iq)
                        mu_     = snl(ie)%mu(ip,iq)
                        syld_   = snl(ie)%syld(ip,iq)
                        ckin_   = snl(ie)%ckin(ip,iq)
                        kkin_   = snl(ie)%kkin(ip,iq)
                        biso_   = snl(ie)%biso(ip,iq)
                        rinf_   = snl(ie)%rinf(ip,iq)
                        
                        ! STRAIN INCREMENT
                        dstrial_(:)   = 0.d0
                        dstrain_(:)   = 0.d0
                        dstrain_(1:2) = dstrain(1:2,ip,iq)
                        dstrain_(4)   = dstrain(3,ip,iq)
                        dpstrain_(:)  = 0.d0
                        
                        ! TRIAL STRESS INCREMENT
                        dstrial_(:)    = dstrial(:,ip,iq)

                        ! CHECK PLASTICITY
                        call check_plasticity(dstrial_,stress_,center_,radius_,&
                            syld_,st_epl,alpha_epl,ip,iq,FS)
                        WRITE(*,*) "ALPHA: ",alpha_epl
                        ! PLASTIC CORRECTION
                        if (st_epl) then
                            dstrain_ = (1.d0-alpha_epl)*dstrain_
                            write(*,*) "plastic"
                            call plastic_corrector(dstrain_,dstrial_,center_,radius_,syld_,&
                                biso_,rinf_,ckin_,kkin_,mu_,lambda_,dpstrain_,ie,FS)
                        endif
                        write(*,*) "STRIAL N: ",DSTRIAL_ 
                        
                        ! STRESS VECTOR
                        snl(ie)%stress(:,ip,iq) = dstrial_(:)
                        ! CENTER
                        snl(ie)%center(:,ip,iq) = center_(:)
                        ! RADIUS
                        snl(ie)%radius(ip,iq)   = radius_
                        ! PLASTIC STRAIN VECTOR
                        snl(ie)%plastic_strain(:,ip,iq) = snl(ie)%plastic_strain(:,ip,iq) + &
                            (/dpstrain_(1:2),dpstrain_(4)/)
                        
                        write(*,*) "STRESS N+1: ",SNL(IE)%STRESS(:,IP,IQ)
                        write(*,*) "=================================="

                    enddo
                enddo
                
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
                        fk(in)     = fk(in)     + fx(ip,iq)/mvec(in)
                        fk(in+nnt) = fk(in+nnt) + fy(ip,iq)/mvec(in+nnt)
                    enddo
                enddo
                ! DEALLOCATE ELEMENT-WISE VARIABLES
                call DEALLOCATE_LOC(ct,ww,dd,dxdx,dxdy,dydx,dydy,dstrain,dstrial,fx,fy)
            enddo
            return
            !
        end subroutine MAKE_INTERNAL_FORCES_NL
        
        
        subroutine STIFF_MATRIX(lambda,mu,DEL)
            implicit none
            ! intent IN
            real*8, intent(in) :: lambda,mu
            ! intent OUT
            real*8, intent(out), dimension(4,4) :: DEL
            !
            DEL(:,:) = 0.d0
            DEL(1,1) = 2*mu
            DEL(2,2) = 2*mu
            DEL(3,3) = 2*mu
            DEL(4,4) = mu
            DEL(1:3,1:3) = DEL(1:3,1:3) + lambda
            !
            return
        end subroutine STIFF_MATRIX

        !****************************************************************************
        ! MISES YIELDING LOCUS AND GRADIENT
        !****************************************************************************

        subroutine mises_yld_locus(stress, center, radius, syld, FM, gradF)
            ! 
            implicit none
            ! intent IN
            real*8,               intent(in)    :: radius,syld
            real*8, dimension(4), intent(in)    :: stress,center
            ! intent OUT
            real*8, dimension(4), intent(out)   :: gradF
            real*8,               intent(out)   :: FM
            !
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,2.0/)
            real*8, dimension(4)                :: dev
            real*8                              :: tau_eq
            integer*4                           :: k
            ! COMPUTE TENSOR COMPONENTS 
            call tensor_components (stress, dev)
            ! COMPUTE MISES FUNCTION
            call tau_mises(dev-center, tau_eq)
            FM = tau_eq-syld-radius
            ! COMPUTE MISES FUNCTION GRADIENT
            gradF = 1.5d0*A*(dev-center)/tau_eq
            !  
            return
            !
        end subroutine mises_yld_locus
        
        !****************************************************************************
        ! OCTAHEDRAL SHEAR STRESS
        !****************************************************************************
        
        subroutine tau_mises(stress,tau_eq)
            !
            implicit none
            ! intent IN
            real*8, dimension(4), intent(in)    :: stress
            ! intent OUT
            real*8,               intent(out)   :: tau_eq
            !
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,2.0/) 
            integer*4                           :: k
            ! COMPUTE OCTAHEDRAL SHEAR STRESS 
            tau_eq = 0.0d0
            do k=1,4
                tau_eq = tau_eq+A(k)*(stress(k)**2)
            end do
            tau_eq = sqrt(1.5*tau_eq)
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
            real*8                              :: press 
            !
            dev = stress
            press = sum(stress(1:3))/3
            dev(1:3) = dev(1:3)-press
            !
            return
            !
        end subroutine tensor_components
            
        !****************************************************************************
        ! CHECK PLASTIC CONSISTENCY (KKT CONDITIONS)
        !****************************************************************************
        
        subroutine check_plasticity (dtrial, stress0, center, radius, syld, &
            st_elp, alpha_epl, ip,iq,FS)
            !     
            implicit none
            ! intent IN
            integer*4, intent(in)                   :: ip,iq 
            real*8,                 intent(in)      :: radius,syld
            real*8, dimension(4),   intent(in)      :: center,stress0
            ! intent INOUT
            real*8, dimension(4),   intent(inout)   :: dtrial
            ! intent OUT
            real*8,                 intent(out)     :: FS,alpha_epl
            logical,                intent(out)     :: st_elp
            !
            integer*4                               :: k
            real*8                                  :: FT,checkload
            real*8, dimension(4)                    :: gradFS, gradFT, stress1
            !
            stress1= dtrial+stress0
            call mises_yld_locus(stress0,center,radius,syld,FS,gradFS)
            call mises_yld_locus(stress1,center,radius,syld,FT,gradFT)
            checkload=sum(gradFS*dtrial)/sum(gradFS**2)/sum(dtrial**2)
           
            if (abs(FS).le.FTOL) then
                if (checkload.ge.-LTOL) then
                    alpha_epl = 0d0
                    st_elp    = .true.
                else
                    if (FT.lt.-FTOL) then
                        alpha_epl = 1d0
                        st_elp    = .false.
                    elseif(FT.gt.FTOL) then
                        call gotoFpegasus(stress0,dtrial,center,radius,syld,10,alpha_epl)
                        st_elp    = .true.
                    endif
                end if
            elseif (FS.lt.-FTOL) then
                if (FT.le.FTOL) then
                    alpha_epl = 1d0
                    st_elp    = .false.
                else
                    call gotoFpegasus(stress0,dtrial,center,radius,syld,1,alpha_epl)
                    st_elp    = .true.
                end if
            elseif (FS.gt.FTOL) then
                write(*,*) "ERROR FS>0"
                alpha_epl  = 0d0
                st_elp = .true.
            end if
            dtrial=stress0+dtrial*alpha_epl
            call mises_yld_locus(dtrial,center,radius,syld,FS,gradFS)
        end subroutine check_plasticity
        
        subroutine plastic_corrector(dEps_alpha,stress,center,syld, &
            radius,biso,Rinf,Ckin,kkin,mu,lambda,dEpl,nel,FM)

            implicit none
            integer*4, intent(in) :: nel
            real*8, dimension(4), intent(inout) :: dEps_alpha,stress,center,dEpl 
            real*8,               intent(inout) :: radius 
            real*8,               intent(in)    :: syld,biso,Rinf,Ckin,kkin,mu,lambda
            real*8, dimension(4,4)              :: DEL
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,0.5/)
            real*8                              :: Ttot,deltaTk,qq,R1,R2,dR1,dR2,err0,err1
            real*8                              :: hard1,hard2,deltaTmin
            real*8, intent(out)                 :: FM
            real(8)                             :: Resk
            logical                             :: flag_fail
            real*8, dimension(4)                :: gradFM,S1,S2,X1,X2
            real*8, dimension(4)                :: dS1,dS2,dX1,dX2,dEpl1,dEpl2
            
            call stiff_matrix(lambda,mu,DEL)
            deltaTk = 1.0d0
            Ttot    = 0.0d0
            deltaTmin = 0.01d0
            do while (Ttot.lt.1d0-FTOL)
                Resk  = 0.0d0
                dS1   = 0.0d0
                dX1   = 0.0d0
                dS2   = 0.0d0
                dX2   = 0.0d0
                dR1   = 0.0d0
                dR2   = 0.0d0
                dEpl1 = 0.0d0
                dEpl2 = 0.0d0

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
                S1 = stress + 0.5d0*(dS1+dS2)
                X1 = center + 0.5d0*(dX1+dX2)
                R1 = radius + 0.5d0*(dR1+dR2)
                dEpl1 = 0.5d0*(dEpl1+dEpl2)
                
                ! ERROR
                call tau_mises(dS2-dS1,err0)
                call tau_mises(S1,err1)

                Resk=0.5d0*max(epsilon(Resk),err0/err1)
                
                if (Resk.le.STOL) then
                    
                    stress = S1
                    center = X1
                    radius = R1
                    dEpl   = dEpl+dEpl1
                    call mises_yld_locus (stress, center,radius,syld,FM,gradFM)
                    if (FM.gt.FTOL) then
                        call drift_corr(stress,center,radius,syld,&
                                biso,Rinf,Ckin,kkin,lambda,mu,dEpl)
                    endif
                   
                    Ttot=Ttot+deltaTk
                    qq = min(0.9d0*sqrt(STOL/Resk),1.1d0) 
                    if (flag_fail) then
                        qq = min(qq,1.0d0)
                    endif
                    flag_fail=.false.
                    deltaTk=qq*deltaTk
                    deltaTk=max(qq*deltaTk,deltaTmin)
                    deltaTk=min(deltaTk,1d0-Ttot)

                else
                    qq=max(0.9d0*sqrt(STOL/Resk),deltaTmin)
                    deltaTk=qq*deltaTk
                    flag_fail=.true.
                end if
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
            real*8, dimension(4), parameter     :: A = (/1.0,1.0,1.0,0.5/)
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
            do k=1,4 
                dEplast(k)=dEplast(k)+dPlast*gradF(k)*A(k)
            end do
            do j = 1,4 ! stress increment
                do k = 1,4
                    dstress(j)=DEL(k,j)*(dstrain(j)-A(k)*dPlast*gradF(k))  
                end do
            end do
            
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
            real*8                             :: temp_vec,FM
            real*8, dimension(4)               :: gradF
            real*8, dimension(4),   parameter  :: A =(/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)             :: DEL
            integer*4                          :: j,k
            
            call stiff_matrix(lambda,mu,DEL)
            call mises_yld_locus(stress,center,radius,syld,FM,gradF)
            
            hard = Ckin+biso*(Rinf-radius)
            hard = hard-kkin*sum(center*gradF)
            
            temp_vec   = 0.0d0
            dPlastMult = 0.0d0

            do k = 1,4
                do j = 1,4
                    temp_vec   = temp_vec+gradF(j)*DEL(j,k)*gradF(k)*A(k)
                    dPlastMult = dPlastMult+gradF(j)*DEL(j,k)*dEps(k)
                end do
            end do
            dPlastMult = max(0.0d0,dPlastMult/(hard+temp_vec))
        end subroutine compute_plastic_modulus


        subroutine hardening_increments(Sigma_ij, R, X_ij, sigma_yld, &
            b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR, dX_ij)

            ! INCREMENTS OF INTRINSIC STATIC VARIABLES

            real*8,               intent(in) :: sigma_yld       ! first yielding limit
            real*8,               intent(in) :: R               ! actual mises radius
            real*8,               intent(in) :: b_lmc, Rinf_lmc ! Lamaitre and Chaboche parameters (isotropic hardening)
            real*8,               intent(in) :: C_lmc, kapa_lmc ! Lamaitre and Chaboche parameters (kinematic hardening)
            real*8, dimension(4), intent(in) :: Sigma_ij        ! actual stress state
            real*8, dimension(4), intent(in) :: X_ij            ! actual back stress state

            real*8,               intent(out):: dR              ! mises radius increment
            real*8, dimension(4), intent(out):: dX_ij           ! back stress increment
            
            real*8                           :: F_mises
            real*8, dimension(4)             :: gradF_mises
            real*8, dimension(4), parameter  :: A = (/1.0,1.0,1.0,0.5/)
            integer*4                        :: k

            ! INCREMENT IN ISOTROPIC HARDENING VARIABLES (R)
            dR = b_lmc*(Rinf_lmc-R)

            ! INCREMENT IN KINEMATIC HARDENING VARIABLES (Xij)
            call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)
            dX_ij(:)=2*A(:)*gradF_mises(:)*C_lmc/3-X_ij(:)*kapa_lmc

        end subroutine hardening_increments

        subroutine drift_corr(stress,center,radius,syld, &
            biso,Rinf,Ckin,kkin,lambda,mu,dEplastic)

            ! DRIFT CORRECTION (RADIAL RETURN)
            real*8, dimension(4), intent(inout) :: stress,center,dEplastic
            real*8,               intent(inout) :: radius
            real*8,               intent(in)    :: lambda,mu,syld,biso,Rinf,Ckin,kkin
            real*8                              :: F0,F1,beta,hard,radiust
            real*8, dimension(4)                :: gradF0,gradF1,dstress,stresst,centert
            real*8, dimension(4),     parameter :: A = (/1.0,1.0,1.0,0.5/)
            real*8, dimension(4,4)              :: DEL
            integer*4                           :: counter,k,j
            
            ! INITIAL PLASTIC CONDITION
            call mises_yld_locus(stress,center,radius,syld,F0,gradF0)
            call stiff_matrix(lambda,mu,DEL)
            do counter=1,5 
                ! COMPUTE HARDENING INCREMENTS
                hard = biso*(Rinf-radius)
                hard = hard + Ckin
                hard = hard - kkin*sum(gradF0*center)
                ! COMPUTE BETA FOR DRIFT CORRECTION
                beta = 0d0
                do j=1,4
                    do k=1,4
                        beta=beta+gradF0(k)*DEL(k,j)*A(j)*gradF0(j)
                    end do
                end do
                beta=F0/(hard+beta)
                ! STRESS-STRAIN-HARDENING CORRECTION
                dstress=0d0
                do k=1,4
                    do j=1,4
                        dstress(k)=dstress(k)-beta*DEL(j,k)*A(j)*gradF0(j)
                    end do
                end do
                stresst = stress+dstress
                centert = center+beta*(2*A*gradF0*Ckin/3-center*kkin)
                radiust = radius+beta*(Rinf-radius)*biso
                
                ! CHECK DRIFT
                call mises_yld_locus(stresst,centert,radiust,syld,F1,gradF1)
                if (abs(F1).gt.abs(F0)) then
                    beta   = F0/sum(gradF0*gradF0)
                    stress = stress-beta*gradF0
                else
                    stress = stresst
                    center = centert
                    radius = radiust
                    dEplastic = dEplastic+beta*A*gradF0
                    if (abs(F1).le.FTOL) then
                        exit
                    else
                        F0     = F1
                        gradF0 = gradF1
                    endif
                endif
                
            enddo
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

            alpha1  = 1.0d0
            alpha0  = 0.0d0
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
                                alpha0=0.0d0
                                F0=Fsave
                            endif
                            exit
                        else
                            alpha0=alpha
                            F0=FM
                        endif
                    end do
                    if (flagxit) then
                        exit
                    endif
                end do
                if (.not.flagxit) then
                    write(*,*) "ERROR IN FINDING F=0 (REVERSAL)"
                    alpha1=1.0d0
                    alpha0=0.0d0
                endif
            end if

            do counter0=1,10
                alpha  = alpha1-F1*(alpha1-alpha0)/(F1-F0)
                stress = start0+alpha*dtrial
                call mises_yld_locus(stress,center,radius,s0,FM,gradF)
                if (abs(FM).le.FTOL) then
                    exit
                else
                    alpha0  =   alpha1
                    alpha1  =   alpha
                    F0      =   F1
                    F1      =   FM
                endif

            end do

            if (FM.gt.FTOL) then
                write(*,*) "WARNING: F>TOL"
            endif
        end subroutine gotoFpegasus
        
end module nonlinear2d
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
