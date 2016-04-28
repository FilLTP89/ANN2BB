        subroutine MAKE_INTERNAL_FORCES_NL(nnt,ne,cs_nnz,cs,sdeg_mat,snl,&
            alfa1,alfa2,beta1,beta2,gamma1,gamma2,dt,displ,fk,mvec,fe)
            real*8,     dimension(:,:),allocatable      :: fxs,fys
            real*8,     dimension(:,:),allocatable      :: sxxs,syys,sxys,szzs
            if (nl_sism.gt.0) then
                allocate(fxs(nn,nn))
                allocate(fys(nn,nn))
                allocate(sxxs(nn,nn))
                allocate(syys(nn,nn))
                allocate(szzs(nn,nn))
                allocate(sxys(nn,nn))
            endif
            if (nl_sism.gt.0) then
                fxs  = 0.d0
                fys  = 0.d0
                sxxs = 0.d0
                syys = 0.d0
                sxys = 0.d0
                szzs = 0.d0
            endif

            implicit none
            ! INTENT IN
            integer*4, intent(in)                       :: nnt,ne,nm,cs_nnz,nl_sism
                
                ! COMPUTE SEISMIC EXTERNAL FORCES 
                if (nl_sism.gt.0) then
                    fe = 0.0d0
                    do iq = 1,nn
                        do ip = 1,nn
                            is = nn*(iq-1)+ip
                            in = cs(cs(ie-1)+is)

                    call MAKE_SEISMIC_MOMENT_NEW(nn,sxxs,syys,szzs,sxys,&
                        check_node_sism,check_dist_node_sism,length_cns,ie,facsmom, &
                        nl_sism,func_type,func_indx,func_data,nf,tt1, &
                        nfunc_data,tag_func)
                    call MAKE_INTERNAL_FORCE_EL(nn,ct,ww,dd,dxdx,dxdy,dydx,dydy,&
                        sxxs,syys,szzs,sxys,fxs,fys)
                    sism(in) = sism(in)         + fxs_el(i,j)
                    sism(in+nnt) = sism(in+nnt) + fys_el(i,j)
                endif

                
                endif

