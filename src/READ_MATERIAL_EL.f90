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
!    along with SPEED.  If not, see <http://www.gnu.org/licenses/>.

!> @brief Reads and stores data from filemate (*.mate)
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0


!> @param[in] file_mat  file name (*.mate)
!> @param[in] nm  number of materials
!> @param[in] propm  material properties (rho,lambda,mu,gamma)
!> @param[in] typem spectral degree use within the material
!> @param[in] tagm  label for materials
!> @param[in] ndX  number of Dirichlet boundary conditions (x-dir)
!> @param[in] ndY  number of Dirichlet boundary conditions (y-dir)
!> @param[in] nnX  number of Neumann boundary conditions (x-dir)
!> @param[in] nnY  number of Neumann boundary conditions (y-dir)
!> @param[in] npX  number of point load volume force (x-dir) 
!> @param[in] npY  number of point load volume force (y-dir) 
!> @param[in] nplX  number of plane wave loads (x-dir)
!> @param[in] nplY  number of plane wave loads (y-dir)
!> @param[in] ntest  test case mode
!> @param[in] na  number of abc boundary conditions
!> @param[in] nMDRM
!> @param[in] nBDRM
!> @param[in] tag_MDRM labels for DRM blocks
!> @param[in] tag_BDRM
!> @param[in] tagstep
!> @param[in] nPDRM
!> @param[in] valPDRM
!> @param[in] fPDRM
!> @param[in] nb_dg  number od DG interface conditions
!> @param[in] nsism  number of sismic loads
!> @param[in] nf  number of functions
!> @param[in] nf_drm  number of DRM functions
!> @param[out] valdX  val_ue of the displacement the Dirichlet boundary (x-dir)
!> @param[out] fdX  number of the time function for the Dir. b.c. (x-dir)
!> @param[out] tagdX  label of the Dir. b.c. (x-dir)
!> @param[out] vdY  val_ue of the displacement the Dirichlet boundary (y-dir)
!> @param[out] fdY  number of the time function for the Dir. b.c. (y-dir)
!> @param[out] tagdY  label of the Dir. b.c. (y-dir)
!> @param[out] valnX amplitude of the Neumann load (x-dir)
!> @param[out] fnX  number of the time function for the Neu. b.c. (x-dir)
!> @param[out] tagnX  label of the Neu. b.c. (x-dir)
!> @param[out] valnY amplitude of the Neumann load (y-dir)
!> @param[out] fnY  number of the time function for the Neu. b.c. (y-dir)
!> @param[out] tagnY  label of the Neu. b.c. (y-dir)
!> @param[out] valpX  amplitude of point load volume force (x-dir)
!> @param[out] fpX  number of the time function for point load volume force (x-dir) 
!> @param[out] valpY  amplitude of point load volume force (y-dir)
!> @param[out] fpY  number of the time function for point load volume force (y-dir) 
!> @param[out] valplX  amplitude for plane wave loads (x-dir)
!> @param[out] fplX  number of the time function for plane wave loads (x-dir)
!> @param[out] tagplX  label of the plane wave load (x-dir)
!> @param[out] valplY  amplitude for plane wave loads (y-dir)
!> @param[out] fplY  number of the time function for plane wave loads (y-dir)
!> @param[out] tagplY  label of the plane wave load (y-dir)
!> @param[out] ftest  function for test mode 
!> @param[out] taga  label for abc boundary conditions
!> @param[out] lab_dg  label for DG interface conditions
!> @param[out] valsism  val_ues for seismic loads
!> @param[out] fsism  time function for seismic loads
!> @param[out] tagsism  label for seismic loads
!> @param[out] func_type  function type 
!> @param[out] func_type_drm  DRM function type 
!> @param[out] func_indx  pointers for data functions
!> @param[out] func_indx_drm  pointers for DRM data functions
!> @param[out] func_data  data functions
!> @param[out] func_data_drm  DRM data functions
!> @param[out] tag_func  label for time functions
!> @param[out] tag_func_drm  label for DRM time functions
!> @param[out] fmax  max frequency of the plane wave load 
!> @param[out] lab_dg_yn label for projection of interface nodes

subroutine READ_MATERIAL_EL(file_mat,nm,propm, typem, tagm, &
    ndX,valdX,fdX,tagdX,ndY,valdY,fdY,tagdY,nnX,valnX,fnX,tagnX,nnY,valnY,fnY,tagnY,    &
    npX,valpX,fpX,npY,valpY,fpY,nplX,valplX,fplX,tagplX,nplY,valplY,fplY,tagplY,        &
    nsism,valsism,fsism,tagsism,na,taga,nMDRM,nBDRM,tag_MDRM,tag_BDRM,tagstep,          & !DRM Scandella 27.09.2005
    nPDRM,valPDRM,fPDRM,nf,func_type,func_indx,func_data,tag_func, & 
    nf_drm,func_type_drm,func_indx_drm,func_data_drm,tag_func_drm,  & !DRM Scandella 11.04.2006
    fmax, ntest, ftest,nb_dg,lab_dg,lab_dg_yn,NLFLAG)

    implicit none

    character*70 :: file_mat      
    integer*4 :: nm,ndT,nnT,ncT,ndX,ndY,nnX,nnY,npX,npY,nplX,nplY,nsism,na,npdT,nf
    integer*4 :: nMDRM,nBDRM,nPDRM       !DRM Scandella 20.10.2005
    integer*4 :: tagstep                 !DRM Scandella 17.10.2005
    integer*4 :: nf_drm                  !DRM Scandella 11.04.2006

    integer*4 :: nb_dg, idg
    integer*4, dimension(nb_dg) :: lab_dg, lab_dg_yn

    integer*4, dimension(nm) :: typem
    real*8, dimension(nm,9) :: propm
    integer*4, dimension(nm) :: tagm
    integer*4, dimension(nf) :: func_type
    integer*4, dimension(nf +1) :: func_indx
    real*8, dimension(*) :: func_data
    integer*4, dimension(nf) :: tag_func

    integer*4, dimension(nf_drm) :: func_type_drm    !DRM Scandella 11.04.2006 
    integer*4, dimension(nf_drm +1) :: func_indx_drm !DRM Scandella 11.04.2006  
    real*8, dimension(*) :: func_data_drm            !DRM Scandella 11.04.2006 
    integer*4, dimension(nf_drm) :: tag_func_drm     !DRM Scandella 11.04.2006


    real*8, dimension(ndX,*) :: valdX
    real*8, dimension(ndY,*) :: valdY
    integer*4, dimension(*) :: fdX,tagdX,fdY,tagdY

    real*8, dimension(nnX,*) :: valnX
    real*8, dimension(nnY,*) :: valnY
    integer*4, dimension(*) :: fnX,tagnX,fnY,tagnY

    real*8, dimension(npX,*) :: valpX
    real*8, dimension(npY,*) :: valpY
    integer*4, dimension(*) :: fpX,fpY

    real*8, dimension(nplX,*) :: valplX
    real*8, dimension(nplY,*) :: valplY
    integer*4, dimension(*) :: fplX,tagplX,fplY,tagplY

    integer*4, dimension(*) :: fPDRM          ! DRM Scandella 20.10.2005
    real*8, dimension(nPDRM,*) :: valPDRM     ! DRM Scandella 20.10.2005

    real*8, dimension(nsism,*) :: valsism  
    integer*4, dimension(*) :: fsism,tagsism

    integer*4, dimension(*) :: taga
    integer*4, dimension(nMDRM) :: tag_MDRM   ! DRM Scandella 27.09.2005
    integer*4, dimension(nBDRM) :: tag_BDRM   ! DRM Scandella 27.09.2005
    integer*4 :: ntest
    integer*4, dimension(ntest) :: ftest

    real*8 :: fmax

    integer*4 :: im,ifunc,idf,func_nd
    integer*4 :: ifunc_drm, func_nd_drm                        ! DRM Scandella 11.04.2006

    integer*4 :: iml                                        
    integer*4 :: idX,idY,inX,inY,ipX,ipY,iplX,iplY,isism,ia
    integer*4 :: iMDRM,iBDRM,iPDRM                             !DRM Scandella 20.10.2005
    integer*4 :: ileft,iright, itest
    integer*4 :: i,j,trash,status

    character*2000 :: input_line
    character*4 :: keyword
    logical,intent(in) :: NLFLAG

    open(23,file=file_mat)

    im = 0

    idX = 0
    idY = 0
    inX = 0
    inY = 0
    ipX = 0
    ipY = 0
    iplX = 0
    iplY = 0
    isism = 0
    ia = 0
    iMDRM = 0  ! DRM Scandella 27.09.2005
    iBDRM = 0  ! DRM Scandella 27.09.2005
    iPDRM = 0  ! DRM Scandella 20.10.2005
    idg = 0
    ifunc = 0
    ifunc_drm = 0 !DRM Scandella 11.04.2006
    itest = 0

    fmax = 0


    if (nf.gt.0) then
        func_indx(1) = 1
    endif

    if (nf_drm.gt.0) then      !DRM Scandella 11.04.2006
        func_indx_drm(1) = 1    !DRM Scandella 11.04.2006
    endif                      !DRM Scandella 11.04.2006

    do 
        read(23,'(A)',IOSTAT = status) input_line

        if (status.ne.0) exit

        keyword = input_line(1:4)

        ileft = 0
        iright = len(input_line)
        do i = 1,iright
            if (input_line(i:i).eq.' ') exit
        enddo
        ileft = i 
        write(*,*) 'material no: ',im
        if (keyword.eq.'MATE') then
            im = im + 1
            if (NLFLAG) then
                read(input_line(ileft:iright),*) tagm(im),typem(im), &
                    propm(im,1),propm(im,2),propm(im,3),propm(im,4),&
                    propm(im,5),propm(im,6),propm(im,7),propm(im,8),&
                    propm(im,9)
            else
                read(input_line(ileft:iright),*) tagm(im),typem(im),&
                    propm(im,1),propm(im,2),propm(im,3),propm(im,4)
            endif
        elseif (keyword.eq.'DIRX') then
            idX = idX + 1
            read(input_line(ileft:iright),*)tagdX(idX),fdX(idX),&
                valdX(idX,1),valdX(idX,2)
        elseif (keyword.eq.'DIRY') then
            idY = idY + 1
            read(input_line(ileft:iright),*)tagdY(idY),fdY(idY),&
                valdY(idY,1),valdY(idY,2)
        elseif (keyword.eq.'NEUX') then
            inX = inX + 1

            read(input_line(ileft:iright),*)tagnX(inX),fnX(inX),&
                valnX(inX,1),valnX(inX,2)
        elseif (keyword.eq.'NEUY') then
            inY = inY + 1
            read(input_line(ileft:iright),*)tagnY(inY),fnY(inY),&
             valnY(inY,1),valnY(inY,2)
        elseif (keyword.eq.'PLOX') then
            ipX = ipX + 1
            read(input_line(ileft:iright),*)fpX(ipX),&
             valpX(ipX,1),valpX(ipX,2),valpX(ipX,3)
        elseif (keyword.eq.'PLOY') then
            ipY = ipY + 1
            read(input_line(ileft:iright),*)fpY(ipY),&
             valpY(ipY,1),valpY(ipY,2),valpY(ipY,3)
        elseif (keyword.eq.'PLAX') then
            iplX = iplX + 1
            read(input_line(ileft:iright),*)fplX(iplX),&
             tagplX(iplX),valplX(iplX,1)

        elseif (keyword.eq.'PLAY') then
            iplY = iplY + 1
            read(input_line(ileft:iright),*)fplY(iplY),&
             tagplY(iplY),valplY(iplY,1)

        elseif (keyword.eq.'SISM') then
            isism = isism + 1
            read(input_line(ileft:iright),*)fsism(isism),&
                tagsism(isism),valsism(isism,1),valsism(isism,2),&
                valsism(isism,3),valsism(isism,4),valsism(isism,5),&
                valsism(isism,6),valsism(isism,7),valsism(isism,8),&
                valsism(isism,9),valsism(isism,10),valsism(isism,11),&
                valsism(isism,12)
        elseif (keyword.eq.'ABSO') then
            ia = ia + 1
            read(input_line(ileft:iright),*)taga(ia)

        elseif (keyword.eq.'SDRM') then                           ! DRM Scandella 20.10.2005  
            read(input_line(ileft:iright),*) tagstep               ! DRM Scandella 20.10.2005 
        elseif (keyword.eq.'MDRM') then                           ! DRM Scandella 27.09.2005   
            iMDRM = iMDRM + 1                                      ! DRM Scandella 27.09.2005
            read(input_line(ileft:iright),*)tag_MDRM(iMDRM)        ! DRM Scandella 27.09.2005  
        !
        elseif (keyword.eq.'BDRM') then                           ! DRM Scandella 27.09.2005 
            iBDRM = iBDRM + 1                                      ! DRM Scandella 27.09.2005    
            read(input_line(ileft:iright),*)tag_BDRM(iBDRM)        ! DRM Scandella 27.09.2005

        elseif (keyword.eq.'PDRM') then                           ! DRM Scandella 20.10.2005 
            iPDRM = iPDRM + 1                                      ! DRM Scandella 20.10.2005    
            read(input_line(ileft:iright),*)fPDRM(iPDRM),&         ! DRM Scandella 20.10.2005 
                valPDRM(iPDRM,1),valPDRM(iPDRM,2),valPDRM(iPDRM,3)! DRM Scandella 20.10.2005

        elseif (keyword.eq.'DGIC') then 
            idg = idg + 1
            read(input_line(ileft:iright),*) lab_dg(idg), lab_dg_yn(idg)            

        elseif (keyword.eq.'TEST') then                                        
            itest = itest + 1                                                        
            read(input_line(ileft:iright),*) ftest(itest)            

        elseif (keyword.eq.'FUNC') then
            ifunc = ifunc + 1
            read(input_line(ileft:iright),*) tag_func(ifunc),&
                func_type(ifunc)
            if (func_type(ifunc).eq.0) then
                func_indx(ifunc +1) = func_indx(ifunc) + 0 
            elseif (func_type(ifunc).eq.1) then
                func_indx(ifunc +1) = func_indx(ifunc) + 2
                read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
            elseif (func_type(ifunc).eq.2) then
                func_indx(ifunc +1) = func_indx(ifunc) + 2
                read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
            elseif (func_type(ifunc).eq.3) then
                read(input_line(ileft:iright),*)trash,trash,func_nd
                func_indx(ifunc +1) = func_indx(ifunc) + 2*func_nd
                read(input_line(ileft:iright),*)trash,trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
            elseif (func_type(ifunc).eq.4) then
                func_indx(ifunc +1) = func_indx(ifunc) + 2
                read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
            elseif (func_type(ifunc).eq.5) then
                func_indx(ifunc +1) = func_indx(ifunc) + 2
                read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
            elseif (func_type(ifunc).eq.12) then
                func_indx(ifunc +1) = func_indx(ifunc) + 3
                read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
            elseif (func_type(ifunc).eq.31) then
                func_indx(ifunc +1) = func_indx(ifunc) + 3
                read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
            elseif (func_type(ifunc).eq.60) then                                   
                read(input_line(ileft:iright),*)trash,trash,func_nd                 
                func_indx(ifunc +1) = func_indx(ifunc) + 2*func_nd                  
                read(input_line(ileft:iright),*)trash,trash,trash,&                 
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)    
            elseif (func_type(ifunc).eq.61) then                                   
                read(input_line(ileft:iright),*)trash,trash,func_nd                 
                func_indx(ifunc +1) = func_indx(ifunc) + 2*func_nd                  
                read(input_line(ileft:iright),*)trash,trash,trash,&                 
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)     
            endif
        elseif (keyword.eq.'FDRM') then                                                           !DRM Scandella 11.04.2006
            ifunc_drm = ifunc_drm + 1                                                              !DRM Scandella 11.04.2006 
            read(input_line(ileft:iright),*)tag_func_drm(ifunc_drm),&                              !DRM Scandella 11.04.2006
                func_type_drm(ifunc_drm)                                                          !DRM Scandella 11.04.2006 
            if (func_type_drm(ifunc_drm).eq.50) then                                               !DRM Scandella 11.04.2006
                read(input_line(ileft:iright),*)trash,trash,func_nd_drm                             !DRM Scandella 11.04.2006
                func_indx_drm(ifunc_drm +1) = func_indx_drm(ifunc_drm) + 3*func_nd_drm              !DRM Scandella 11.04.2006
                read(input_line(ileft:iright),*)trash,trash,trash,&                                 !DRM Scandella 11.04.2006
                    (func_data_drm(j), j = func_indx_drm(ifunc_drm),func_indx_drm(ifunc_drm +1) -1)!DRM Scandella 11.04.2006
            endif

        elseif (keyword.eq.'FMAX') then
            read(input_line(ileft:iright),*)fmax

        endif

    enddo


    close(23)

    return

end subroutine READ_MATERIAL_EL
