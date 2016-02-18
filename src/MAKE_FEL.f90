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

!> @brief Computes external loads.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0
!> @param[in] nnode number of nodes
!> @param[in] xs x-coordinate of GLL nodes
!> @param[in] ys y-coordinate of GLL nodes
!> @param[in] cs_nnz length of cs
!> @param[in] cs spectral connectivity vector
!> @param[in] nm number of materials
!> @param[in] tag_mat material labels
!> @param[in] type_mat material type (dummy)
!> @param[in] sdeg_mat polynomial degree vector
!> @param[in] tref_mat dummy
!> @param[in] prop_mat  material properties (rho, lambda, mu, gamma)
!> @param[in] ne number of elements
!> @param[in] a1 costant values for the bilinear map
!> @param[in] a2 costant values for the bilinear map
!> @param[in] b1 costant values for the bilinear map
!> @param[in] b2 costant values for the bilinear map
!> @param[in] g1 costant values for the bilinear map
!> @param[in] g2 costant values for the bilinear map
!> @param[in] cs_nnz_bc length of cs_bc
!> @param[in] cs_bc local spectral boundary connectivity vector
!> @param[in] nl_***  number of load of type ***
!> @param[in] val_***  value for load type ***
!> @param[in] fun_*** func value for load type ***
!> @param[in] tag_***  tag for load type ***
!> @param[in] node_TOT nodes of total DRM domain
!> @param[in] nMDRM
!> @param[in] tag_MDRM labels for DRM blocks
!> @param[in] ntest  number of functions in test mode
!> @param[in] fun_test function for test mode
!> @param[in] nfunc  number of functions
!> @param[in] tag_func label for functions
!> @param[in] nquad number of hexahedras
!> @param[in] con_quad connectivity matrix
!> @param[in] nline number of hexahedras
!> @param[in] con_line connectivity matrix
!> @param[in] length_cns  lenght check seismic nodes
!> @param[in] max_num_node_sism max number of seismic nodes
!> @param[in] num_node_sism number of seismic loads
!> @param[in] sour_node_sism infos about  seismic nodes  
!> @param[in] facsmom  seismic moment factor
!> @param[in] tagstep label of the step of DRM Analysis
!> @param[in] nf_drm number of DRM functions
!> @param[in] glob_x
!> @param[in] glob_y
!> @param[in] dist_sour_node_sism distance from hypo
!> @param[in] test 1 if test mode is active, 0 otherwise
!> @param[out] fmat vector of external applied loads

      subroutine MAKE_FEL(nnode,xs,ys,cs_nnz,cs,&
                          nm,tag_mat,sdeg_mat,prop_mat,&
                          ne,a1,b1,g1,a2,b2,g2,&
                          cs_nnz_bc,cs_bc,&
                          nl_dirX,val_dirX,fun_dirX,tag_dirX,&
                          nl_dirY,val_dirY,fun_dirY,tag_dirY,&
                          nl_neuX,val_neuX,fun_neuX,tag_neuX,&
                          nl_neuY,val_neuY,fun_neuY,tag_neuY,&
                          nl_poiX,val_poiX,fun_poiX,&
                          nl_poiY,val_poiY,fun_poiY,&
                          nl_plaX,val_plaX,fun_plaX,tag_plaX,&
                          nl_plaY,val_plaY,fun_plaY,tag_plaY,&
                          nl_sism,val_sism,fun_sism,tag_sism,&
						  nl_PDRM,val_TOT,node_TOT,&             !DRM Scandella 12.04.2006
						  nMDRM,tag_MDRM,&                       !DRM Scandella 26.10.2005  
                          nfunc,tag_func,Fmat,&
                          xx,yy,con_quad,nquad,con_line,nline,&
                          num_node_sism,max_num_node_sism,&
                          sour_node_sism,dist_sour_node_sism,&
                          length_cns,facsmom, &
						  tagstep,&                              !DRM Scandella 21.10.2005
						  nf_drm,glob_x,glob_y, &                !DRM Scandella 11.04.2006
						  test, ntest, fun_test)      


      implicit none
      
      integer*4 :: nnode,cs_nnz,nm,ne,cs_nnz_bc,nquad,nline
      integer*4 :: nl_dirX,nl_dirY,nl_neuX,nl_neuY, test
	  integer*4 :: nl_PDRM,nMDRM                                         !DRM Scandella 02.11.2005
      integer*4, dimension(nMDRM) :: tag_MDRM                            !DRM Scandella 21.10.2005 
	  integer*4 :: nfunc, nf_drm

      real*8, dimension(nnode) :: xs,ys

      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tag_mat,sdeg_mat
      real*8, dimension(nm,4) :: prop_mat

      integer*4 :: nl_poiX, nl_poiY, nl_plaX, nl_plaY,nl_sism
      real*8, dimension(ne) :: a1,b1,g1,a2,b2,g2
      integer*4, dimension(0:cs_nnz_bc) :: cs_bc
      real*8, dimension(nl_dirX,2) :: val_dirX
      real*8, dimension(nl_dirY,2) :: val_dirY
      integer*4, dimension(nl_dirX) :: fun_dirX
      integer*4, dimension(nl_dirY) :: fun_dirY
      integer*4, dimension(nl_dirX) :: tag_dirX
      integer*4, dimension(nl_dirY) :: tag_dirY
      real*8, dimension(nl_neuX,2) :: val_neuX
      real*8, dimension(nl_neuY,2) :: val_neuY
      integer*4, dimension(nl_neuX) :: fun_neuX
      integer*4, dimension(nl_neuY) :: fun_neuY
      integer*4, dimension(nl_neuX) :: tag_neuX
      integer*4, dimension(nl_neuY) :: tag_neuY
      real*8, dimension(nl_plaX,1) :: val_plaX
      real*8, dimension(nl_plaY,1) :: val_plaY
      integer*4, dimension(nl_plaX) :: fun_plaX
      integer*4, dimension(nl_plaY) :: fun_plaY
      integer*4, dimension(nl_plaX) :: tag_plaX
      integer*4, dimension(nl_plaY) :: tag_plaY
      real*8, dimension(nl_poiX,3) :: val_poiX
      real*8, dimension(nl_poiY,3) :: val_poiY
      integer*4, dimension(nl_poiX) :: fun_poiX
      integer*4, dimension(nl_poiY) :: fun_poiY
      real*8, dimension(nl_sism,12) :: val_sism
      integer*4, dimension(nl_sism) :: fun_sism
      integer*4, dimension(nl_sism) :: tag_sism
	  real*8, dimension(nl_PDRM) ::     val_TOT   !DRM Scandella 21.10.2005
      integer*4 :: ntest
      integer*4, dimension(ntest) :: fun_test 


      integer*4, dimension(nfunc) :: tag_func
      real*8, dimension(nl_sism,3) :: facsmom
      real*8, dimension(nl_sism) :: term_vet
      real*8, dimension(nfunc,2*nnode) :: Fmat
      
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8, dimension(:), allocatable :: dxdx,dydx,dxdy,dydy
      real*8, dimension(:,:), allocatable :: det_j
      integer*4, dimension(nl_poiX) :: node_poiX
      integer*4, dimension(nl_poiY) :: node_poiY
	  integer*4, dimension(nl_PDRM) :: node_TOT !DRM Scandella 12.04.2006 
	  integer*4, dimension(nl_PDRM,2) ::glob_x  !DRM Scandella 11.04.2006 
	  integer*4, dimension(nl_PDRM,2) ::glob_y  !DRM Scandella 11.04.2006 

      integer*4, dimension(nl_sism) :: num_node_sism
      integer*4 :: max_num_node_sism
      integer*4, dimension(max_num_node_sism,nl_sism) :: sour_node_sism
      real*8, dimension(max_num_node_sism,nl_sism) :: dist_sour_node_sism
      integer*4 :: length_cns

      integer*4 :: im,if,ie,ip,ipl,isism,il,nedge,nn,fn
      integer*4 :: is,in,id1,id2
      integer*4 :: i,j,k,h,l,m,ii
      real*8 :: rho,lambda,mu,x,y
      real*8 :: lx,ly,ll,v1,v2,v,term
      real*8 :: rp2,t1x,t1y,t2x,t2y, ellex
      real*8, dimension(nnode) :: xx,yy
      integer*4, dimension(nquad,5) :: con_quad
      integer*4, dimension(nline,3) :: con_line
      integer*4 :: C,sit
      integer*4, dimension(4) :: num_sit
      integer*4, dimension(8) :: coord_sit
      real*8 :: slip1,slip2,norm1,norm2,amp_sism
      integer*4 :: tagstep                               !DRM Scandella 21.10.2005 
      real*8 :: pi
 
      pi = 4.d0*datan(1.d0)

      nn = 2
      
      allocate(ct(nn),ww(nn),dd(nn,nn))
      call lgl(nn,ct,ww,dd)
      allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))
      
      num_sit = (/1,2,2,1/)
      coord_sit = (/3,2,5,4,4,3,2,5/)
      length_cns = 0

      Fmat = 0.d0
      
      do i = 1,nl_poiX
         call FIND_NEAREST_NODE(nnode,xs,ys,val_poiX(i,1),val_poiX(i,2),node_poiX(i))
      enddo
      do i = 1,nl_poiY
         call FIND_NEAREST_NODE(nnode,xs,ys,val_poiY(i,1),val_poiY(i,2),node_poiY(i))
      enddo

      ne = cs(0) -1
      
      do im = 1,nm
         if ((sdeg_mat(im) +1) .ne. nn) then
            deallocate(ct,ww,dd)
            deallocate(dxdx,dydx,dxdy,dydy,det_j)
            
            nn = sdeg_mat(im) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call lgl(nn,ct,ww,dd)
            allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))
         endif
         
         rho = prop_mat(im,1)
         lambda = prop_mat(im,2)
         mu = prop_mat(im,3)
         
         do ie = 1,ne
            if (cs(cs(ie -1) +0).eq.tag_mat(im)) then
               do i = 1,nn
                  dxdy(i) = b1(ie) + g1(ie) * ct(i)
                  dydy(i) = b2(ie) + g2(ie) * ct(i)
               enddo
               
               do j = 1,nn
                  dxdx(j) = a1(ie) + g1(ie) * ct(j)
                  dydx(j) = a2(ie) + g2(ie) * ct(j)
               enddo
               
               do j = 1,nn
                  do i = 1,nn
                     det_j(i,j) = dxdx(j)*dydy(i) - dxdy(i)*dydx(j)
                  enddo
               enddo
         
               
               if (ntest .gt. 0) then 
                  do ip = 1,ntest
                     fn = 0
                     do if = 1,nfunc
                        if (fun_test(ip).eq.tag_func(if)) fn = if
                     enddo
                     
                    ! write(*,*) 'fn test', fn
                    ! read(*,*)
                     
                     if (fn.gt.0) then

                           do j = 1,nn
                              do i = 1,nn
                                 
                                 is = nn*(j -1) +i
                                 in = cs(cs(ie -1) +is)
                                 x =  xs(in); y = ys(in);           
          
                                !fmat(fn,in) = fmat(fn,in) + det_j(i,j) * ww(i) * ww(j) * ( &
                                !              - mu*(2.d0*x**2 - 2.d0) - (2.d0*y**2 - 2.d0)*(lambda + 2.d0*mu))
 
                                !fmat(fn,in+nnode) = fmat(fn,in+nnode) + det_j(i,j) * ww(i) * ww(j) * (&
                                !                   - 4.d0*x*y*(lambda + mu)) 
                
                                fmat(fn,in) = fmat(fn,in) + det_j(i,j) * ww(i) * ww(j) * ( &   
                                                  -2.d0*pi**2*dcos(pi*y)*dsin(pi*y)*( &
                                                  rho*2.d0*dcos(pi*x)**2 - 2.d0*rho - 8.d0*mu*dcos(pi*x)**2 + 6.d0*mu) &
                                                  )
                                
                                fmat(fn,in+nnode) = fmat(fn,in+nnode) + det_j(i,j) * ww(i) * ww(j) * (& 
                                                  2.d0*pi**2*dcos(pi*x)*dsin(pi*x)*( &
                                                  rho*2.d0*dcos(pi*y)**2 - 2.d0*rho - 8.d0*mu*dcos(pi*y)**2 + 6.d0*mu) &
                                                  )
                                !fmat(fn,in) = fmat(fn,in) + sin(2.d0*pi*ys(in)) * ( 1.d0 - 3.d0*(sin(pi*xs(in)))**2.d0 ) &
                                !        * (det_j(i,j) * ww(i) * ww(j))
		                !fmat(fn,in+nnode) = fmat(fn,in+nnode) -sin(2.d0*pi*xs(in)) &
		                !                   * ( 1.d0 - 3.d0*(sin(pi*ys(in)))**2.d0 ) &
                                !                   * (det_j(i,j) * ww(i) * ww(j))

                                
                                
  
                              enddo
                           enddo               
                      endif
                  enddo
               endif                  
                              
! Point load X
                  
               if (nl_poiX.gt.0) then
                  do ip = 1,nl_poiX
                     fn = 0
                     do if = 1,nfunc
                        if (fun_poiX(ip).eq.tag_func(if)) fn = if
                     enddo
                     
                     if (fn.gt.0) then
                        do j = 1,nn
                           do i = 1,nn
                              is = nn*(j -1) +i
                              in = cs(cs(ie -1) +is)
                              
                                 
                              if (in .eq. node_poiX(ip)) then
                                    term = val_poiX(ip,3) * (det_j(i,j) * ww(i) * ww(j))
                                    fmat(fn,in) = fmat(fn,in) + term
						
                              endif
                           enddo
                        enddo
                     endif
                  enddo
               endif
              
               
! Point load Y
               
               if (nl_poiY.gt.0) then
                  do ip = 1,nl_poiY
                     fn = 0
                     do if = 1,nfunc
                        if (fun_poiY(ip).eq.tag_func(if)) fn = if
                     enddo
                     
                     if (fn.gt.0) then
                        do j = 1,nn
                           do i = 1,nn
                              is = nn*(j -1) +i
                              in = cs(cs(ie -1) +is)

                              if (in .eq. node_poiY(ip)) then
                                    term = val_poiY(ip,3) * (det_j(i,j) * ww(i) * ww(j))
                                    fmat(fn,in + nnode) = fmat(fn,in + nnode) + term
                              endif
                           enddo
                        enddo
                     endif
                  enddo
               endif
!---------------------------------------------------------------------------
! DRM Scandella 21.10.2005  

! Point PDRM
               if (tagstep.eq.2) then                             ! DRM Scandella 21.10.2005  
                 if (nl_PDRM.gt.0) then                           ! DRM Scandella 21.10.2005 
                   do ip = 1,nl_PDRM                              ! DRM Scandella 21.10.2005                     
                        do j = 1,nn                               ! DRM Scandella 21.10.2005 
                           do i = 1,nn                            ! DRM Scandella 21.10.2005 
                              is = nn*(j -1) +i                   ! DRM Scandella 21.10.2005 
                              in = cs(cs(ie -1) +is)              ! DRM Scandella 21.10.2005 
                                 
                                 if (in.eq.node_TOT(ip)) then     ! DRM Scandella 21.10.2005 
							         term = val_TOT(ip)           ! DRM Scandella 21.10.2005 
								     glob_x(ip,1) = in            ! DRM Scandella 21.10.2005 
                                     glob_y(ip,1) = in + nnode    ! DRM Scandella 21.10.2005 
								     glob_x(ip,2) = term          ! DRM Scandella 21.10.2005 
                                     glob_y(ip,2) = term          ! DRM Scandella 21.10.2005                               
                                 endif                            ! DRM Scandella 21.10.2005 
 
                            enddo                                 ! DRM Scandella 21.10.2005 
                         enddo                                    ! DRM Scandella 21.10.2005 
			        enddo                                         ! DRM Scandella 21.10.2005 
                 endif                                            ! DRM Scandella 21.10.2005 
			   endif                                              ! DRM Scandella 21.10.2005 

!End DRM Scandella 21.10.2005
!---------------------------------------------------------------------------

! Plane Wave X
               if (nl_plaX.gt.0) then


                 ! Find out plane wave condition load

                  do ipl = 1,nl_plaX

                     !Check on the Plane Wave Material - start

                     if (tag_mat(im).eq.tag_plaX(ipl)) then
                        C=dsqrt(mu*rho)
                        fn = 0
                        do if = 1,nfunc
                           if (fun_plaX(ipl).eq.tag_func(if)) fn = if
                        enddo
                        if (fn.gt.0) then

                           !Table of the possible situation:
                           !hypothesis:
                           !- the plane wave load source is placed on a rectangular element, aligned 
                           !  along the principal axis (x,y). It does NOT work properly on a deformed 
                           !  element!!!
                           !- the node orientation is ALWAYS counter-clockwise!!!
                           !
                           !    S1        S2        S3        S4
                           ! 2------1  4------3  3------2  1------4
                           ! |      |  |      |  |      |  |      |
                           ! |      |  |      |  |      |  |      |
                           ! |      |  |      |  |      |  |      |
                           ! 3------4  1------2  4------1  2------3

                           !****** Situation S1 or S2 - begin ******

                           if (yy(con_quad(ie,2)).eq.yy(con_quad(ie,3))) then
                            
                              !*** Situation S1 - this means y1 = y2 and y1 > y4 ***
                              !
                              ! e.g.: spectral degree = 3, nn = 4
                              ! [#] = macro nodes
                              !  #  = micro nodes
                              !
                              !          S1
                              !   [2]4--3--2--1[1]
                              !    |            |
                              !  ->|8---7--6---5|<- Applied plane wave load
                              !    |            |
                              !    |12--11-10--9|
                              !    |            |
                              !   [3]16-15-14-13[4]
                              !
                              ! e.g.: spectral degree = 4, nn = 5
                              !
                              !           S1
                              !   [2]5--4--3--2--1[1]
                              !    |               |
                              !    |10--9--8--7---6|
                              !    |               |
                              !  ->|15-14-13-12--11|<- Applied plane wave load
                              !    |               |
                              !    |20-19-18-17--16|
                              !    |               |
                              !   [3]25-24-23-22-21[4]

                              if (yy(con_quad(ie,2)).gt.yy(con_quad(ie,5))) then
                              
                                 sit = 1
                           
                              !*** Situation S2 - this means y1=y2 and y4 > y1 ***

                              else
                              
                                 sit = 2
                           
                              endif

                              ellex=dabs(xx(con_quad(ie,coord_sit((sit-1)*2+1)))&
                                    - xx(con_quad(ie,coord_sit((sit-1)*2+2))))

                              if (mod(nn-1,2).eq.0) then
                                j = ((nn-1)/2)+1
                              else
                                j = (int((nn-1)/2))+num_sit(sit)
                              endif
                                
                              do i = 1,nn
                                 is = nn*(j -1) +i
                                 in = cs(cs(ie -1) +is)
                              
                                 id1 = in
                                 
                                 term = C*val_plaX(ipl,1) * ellex*ww(i)
                                 fmat(fn,id1) = fmat(fn,id1) + term

                              enddo
                              
                           !****** Situation S1 or S2 - end ******

                           !****** Situation S3 or S4 - begin ******

                           else
                            
                              !*** Situation S3 - this means y2=y3 and y2 > y1 ***

                              if (yy(con_quad(ie,3)).gt.yy(con_quad(ie,2))) then
                            
                                 sit = 3

                              !*** Situation S4 - this means y2=y3 and y1 > y2 ***
                           
                              else

                                 sit = 4

                              endif

                              ellex=dabs(xx(con_quad(ie,coord_sit((sit-1)*2+1)))&
                                    - xx(con_quad(ie,coord_sit((sit-1)*2+2))))

                              if (mod(nn-1,2).eq.0) then
                                i = ((nn-1)/2)+1
                              else
                                i = (int((nn-1)/2))+num_sit(sit)
                              endif

                              do j = 1,nn
                                 is = (j-1)*nn+i
                                 in = cs(cs(ie -1) +is)
                              
                                 id1 = in
                                 term = C*val_plaX(ipl,1) *ellex*ww(j)
                                 fmat(fn,id1) = fmat(fn,id1) + term

                              enddo                          

               
                           endif

                           !****** Situation S3 or S4 - end ******

                        endif
                    
                     endif

                     !Check on the Plane Wave Material - end

                  enddo
               endif

! Plane Wave Y

               if (nl_plaY.gt.0) then


                 ! Find out plane wave condition load

                  do ipl = 1,nl_plaY

                     !Check on the Plane Wave Material - start

                     if (tag_mat(im).eq.tag_plaY(ipl)) then
                        C=dsqrt((lambda+2*mu)*rho)
                        fn = 0
                        do if = 1,nfunc
                           if (fun_plaY(ipl).eq.tag_func(if)) fn = if
                        enddo
                        if (fn.gt.0) then

                           !****** Situation S1 or S2 - begin ******

                          if (yy(con_quad(ie,2)).eq.yy(con_quad(ie,3))) then
                            
                              !*** Situation S1 - this means y1 = y2 and y1 > y4 ***

                              if (yy(con_quad(ie,2)).gt.yy(con_quad(ie,5))) then
                              
                                 sit = 1
                           
                              !*** Situation S2 - this means y1=y2 and y4 > y1 ***

                              else
                              
                                 sit = 2
                           
                              endif

                              ellex=dabs(xx(con_quad(ie,coord_sit((sit-1)*2+1)))&
                                    - xx(con_quad(ie,coord_sit((sit-1)*2+2))))

                              if (mod(nn-1,2).eq.0) then
                                j = ((nn-1)/2)+1
                              else
                                j = (int((nn-1)/2))+num_sit(sit)
                              endif
                                
                              do i = 1,nn
                                 is = nn*(j -1) +i
                                 in = cs(cs(ie -1) +is)
                              
                                 id2 = in + nnode
                                 
                                 term = C*val_plaY(ipl,1) * ellex*ww(i)
                                 fmat(fn,id2) = fmat(fn,id2) + term

                              enddo
                              
                           !****** Situation S1 or S2 - end ******

                           !****** Situation S3 or S4 - begin ******

                           else
                            
                              !*** Situation S3 - this means y2=y3 and y2 > y1 ***

                              if (yy(con_quad(ie,3)).gt.yy(con_quad(ie,2))) then
                            
                                 sit = 3

                              !*** Situation S4 - this means y2=y3 and y1 > y2 ***
                           
                              else

                                 sit = 4

                              endif

                              ellex=dabs(xx(con_quad(ie,coord_sit((sit-1)*2+1)))&
                                    - xx(con_quad(ie,coord_sit((sit-1)*2+2))))

                              if (mod(nn-1,2).eq.0) then
                                i = ((nn-1)/2)+1
                              else
                                i = (int((nn-1)/2))+num_sit(sit)
                              endif

                              do j = 1,nn
                                 is = (j-1)*nn+i
                                 in = cs(cs(ie -1) +is)
                              
                                 id2 = in + nnode
                                 term = C*val_plaY(ipl,1) *ellex*ww(j)
                                 fmat(fn,id2) = fmat(fn,id2) + term
                                  
                              enddo                          

               
                           endif

                           !****** Situation S3 or S4 - end ******

                        endif
                    
                     endif

                     !Check on the Plane Wave Material - end

                  enddo
               endif



! Seismic Moment Load

! Seismic moment scale factor - begin

               if (nl_sism.gt.0) then

                  do isism = 1,nl_sism

          	     do k = 1,num_node_sism(isism)        

			do j = 1,nn
                           do i = 1,nn
                              is = nn*(j -1) +i
                              in = cs(cs(ie -1) +is)
				
                                 if (in.eq.sour_node_sism(k,isism)) then
				    facsmom(isism,1) = facsmom(isism,1) + det_j(i,j) * ww(i) * ww(j)
                                    facsmom(isism,2) = facsmom(isism,2) + det_j(i,j) * ww(i) * ww(j)
                                    facsmom(isism,3) = facsmom(isism,3) + det_j(i,j) * ww(i) * ww(j)
                                    length_cns = length_cns + 1
                                 endif
                                 
			    enddo
		        enddo

                      enddo
	          enddo
	       endif
                                    
! Seismic moment scale factor - end 
       
               
            endif
         enddo
      enddo

      if (nl_sism.gt.0) then
                              
          do isism = 1,nl_sism
             slip1 = val_sism(isism,7)
             slip2 = val_sism(isism,8)
             norm1 = val_sism(isism,9)
             norm2 = val_sism(isism,10)
             amp_sism = val_sism(isism,11)

             facsmom(isism,1) = 1/facsmom(isism,1) &
                                       * (slip1*norm1+slip1*norm1) &
                                       * amp_sism
             facsmom(isism,2) = 1/facsmom(isism,2) &
                                       * (slip2*norm2+slip2*norm2) &
                                       * amp_sism
             facsmom(isism,3) = 1/facsmom(isism,3) &
                                       * (slip1*norm2+slip2*norm1) &
                                       * amp_sism
      
          enddo
      endif
      
      if (cs_nnz_bc.gt.0) then
         nedge = cs_bc(0) -1
         
         do ie = 1,nedge
            if ((cs_bc(ie) - cs_bc(ie -1) -1).ne.nn) then
               deallocate(ct,ww,dd)
               nn = cs_bc(ie) - cs_bc(ie -1) -1
               allocate(ct(nn),ww(nn),dd(nn,nn))
               call lgl(nn,ct,ww,dd)
            endif
            
            lx = xs(cs_bc(cs_bc(ie -1) +nn)) - xs(cs_bc(cs_bc(ie -1) +1))
            ly = ys(cs_bc(cs_bc(ie -1) +nn)) - ys(cs_bc(cs_bc(ie -1) +1))
            ll = dsqrt(lx*lx + ly*ly)
            
 
! Neumann load X
            
            if (nl_neuX.gt.0) then
               do il = 1,nl_neuX
                  if (tag_neuX(il).eq.cs_bc(cs_bc(ie -1) +0)) then
                     fn = 0
                     do if = 1,nfunc
                        if (fun_neuX(il).eq.tag_func(if)) fn = if
                     enddo
                     
                     if (fn.gt.0) then
                        v1 = val_neuX(il,1)
                        v2 = val_neuX(il,2)
                        
                        do i = 1,nn
                           in = cs_bc(cs_bc(ie -1) +i)
                           
                           id1 = in
                           v = 0.5d0*((ct(i) + 1.0d0)*v2 - (ct(i) - 1.0d0)*v1)
                           term = 0.5d0 * ll * ww(i) * v
                              
                           fmat(fn,id1) = fmat(fn,id1) + term

                        enddo
                     endif
                  endif
               enddo
            endif
            
            
! Neumann load Y
            
            if (nl_neuY.gt.0) then
               do il = 1,nl_neuY
                  if (tag_neuY(il).eq.cs_bc(cs_bc(ie -1) +0)) then
                     fn = 0
                     do if = 1,nfunc
                        if (fun_neuY(il).eq.tag_func(if)) fn = if
                     enddo
                     
                     if (fn.gt.0) then
                        v1 = val_neuY(il,1)
                        v2 = val_neuY(il,2)
                        
                        do i = 1,nn
                           in = cs_bc(cs_bc(ie -1) +i)
                           
                           id2 = in + nnode
                           v = 0.5d0*((ct(i) + 1.0d0)*v2 - (ct(i) - 1.0d0)*v1)
                           term = 0.5d0 * ll * ww(i) * v
                              
                           fmat(fn,id2) = fmat(fn,id2) + term

                        enddo
                     endif
                  endif
               enddo
            endif
            
            
           ! write(*,*) fmat
           ! read(*,*)
            
! Dirichlet X
            
            if (nl_dirX.gt.0) then
               do il = 1,nl_dirX
                  if (tag_dirX(il).eq.cs_bc(cs_bc(ie -1) +0)) then
                     fn = 0
                     do if = 1,nfunc
                        if (fun_dirX(il).eq.tag_func(if)) fn = if
                     enddo
                     
                     if (fn.gt.0) then
                        v1 = val_dirX(il,1)
                        v2 = val_dirX(il,2)
                        
                        do i = 1,nn
                           in = cs_bc(cs_bc(ie -1) +i)
                           
                           id1 = in
                           !do if = 1,nfunc
                           !   fmat(if,id1) = 0.0d0
                           !enddo
                              
                           v = 0.5d0*((ct(i) + 1.0d0)*v2 - (ct(i) - 1.0d0)*v1)
                           fmat(fn,id1) = v

                        enddo
                     endif
                  endif
               enddo
            endif
            
            
! Dirichlet Y
            
            if (nl_dirY.gt.0) then
               do il = 1,nl_dirY
                  if (tag_dirY(il).eq.cs_bc(cs_bc(ie -1) +0)) then
                     fn = 0
                     do if = 1,nfunc
                        if (fun_dirY(il).eq.tag_func(if)) fn = if
                     enddo
                     
                     if (fn.gt.0) then
                        v1 = val_dirY(il,1)
                        v2 = val_dirY(il,2)
                        
                        do i = 1,nn
                           in = cs_bc(cs_bc(ie -1) +i)
                           
                           id2 = in + nnode
                           !do if = 1,nfunc
                           !      fmat(if,id2) = 0.0d0
                           !enddo
                              
                           v = 0.5d0*((ct(i) + 1.0d0)*v2 - (ct(i) - 1.0d0)*v1)
                           fmat(fn,id2) = v

                        enddo
                     endif
                  endif
               enddo
            endif
            
         enddo
      endif

           ! write(*,*) fmat
           ! read(*,*)
      
      
      deallocate(ct,ww,dd)
      deallocate(dxdx,dydx,dxdy,dydy,det_j)
      
      return
      end subroutine MAKE_FEL

