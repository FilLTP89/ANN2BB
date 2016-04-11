
      program SEelse2D

!***********************************************************************
!     
!     Spectral Element solver for 2D elasto-viscoplastic dynamics
!     
!     F. Maggio, L. Massidda, J. Sabadell, G. Siddi  -  CRS4 MCI/SSM  
!     http://www.crs4.it/Areas/wms.html
!
!     Marco Stupazzini & Clara Zambelli  -  POLIMI - DIS
!     http://stru.polimi.it
!
!     
!     © CRS4, 1999-2005, All Rights Reserved
!     
!     Disclaimer:
!     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
!     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
!     NONINFRINGEMENT. IN NO EVENT SHALL CRS4 OR POLIMI BE LIABLE FOR ANY 
!     CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
!     TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
!     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!     
!     Last revision LucaM -> ??-??-????
!     Last revision M&C   -> 29-10-2004
!     
!***********************************************************************
      
      
      implicit none
      
      integer*4, dimension(3)  :: clock
      integer*4 :: start,finish
      real*4    :: time_in_seconds
      
      integer*4 :: nargs
      integer*4 :: iargc
      
      real*8 :: deltat, xtime
      integer*4 :: opt_out_form,opt_out_data
      
      real*8 :: deltat_cfl,fmax
      character*3 :: deltat_fixed        
      
!     MODEL VARIABLES
      
      real*8, dimension (:), allocatable :: xx_macro,yy_macro
      real*8, dimension (:), allocatable :: xx_spx,yy_spx
      integer*4 :: nnod_macro,nnod
      
      integer*4, dimension (:,:), allocatable :: con
      integer*4, dimension (:), allocatable :: con_spx
      integer*4 :: nelem,con_nnz
      
      integer*4, dimension (:,:), allocatable :: con_bc
      integer*4, dimension (:), allocatable :: con_spx_bc
      integer*4 ::nedge,con_nnz_bc
      
      integer*4, dimension (:), allocatable :: Ebw,Ebin,Nbw,Nbin
      integer*4 :: Ennz,Nnnz
      
      real*8, dimension (:), allocatable :: alfa1,beta1,gamma1,delta1
      real*8, dimension (:), allocatable :: alfa2,beta2,gamma2,delta2
      
      
      
!     MECHANICAL VARIABLES
      
      integer*4, dimension (:), allocatable :: sdeg_mat

      integer*4, dimension (:), allocatable :: type_mat

      real*8, dimension (:), allocatable :: tref_mat

      real*8, dimension (:,:), allocatable :: prop_mat

      integer*4, dimension(:), allocatable :: tag_mat

      integer*4 :: nmat



      integer*4, dimension (:), allocatable :: type_matg ! Mexico Paco grad_lin 29.11.2004

      real*8, dimension (:), allocatable :: tref_matg ! Mexico Paco grad_lin 29.11.2004

      real*8, dimension (:,:), allocatable :: prop_matg ! Mexico Paco grad_lin 29.11.2004

      integer*4, dimension(:), allocatable :: tag_matg ! Mexico Paco grad_lin 29.11.2004

      integer*4 :: nmatg ! Mexico Paco grad_lin 29.11.2004
      
      real*8, dimension (:), allocatable :: u0,u_1
      real*8, dimension (:), allocatable :: Mel,Cel,KCel
      real*8, dimension (:,:), allocatable :: Fel
      
      integer*4, dimension (:),   allocatable :: Kbin_el,Kbw_el
      integer*4 :: Knnz_el
      
      real*8, dimension (:,:), allocatable :: val_dirX_el,val_dirY_el
      real*8, dimension (:,:), allocatable :: val_neuX_el,val_neuY_el

	  real*8, dimension (:,:), allocatable :: val_neuN_el,val_neuT_el ! M Mexico Arroyo
      real*8, dimension (:,:), allocatable :: val_poiX_el,val_poiY_el

	  real*8, dimension (:,:), allocatable :: val_palX_el,val_palY_el ! Marco

	  real*8, dimension (:,:), allocatable :: val_dipX_el,val_dipY_el ! Marco
      real*8, dimension (:,:), allocatable :: val_plaX_el,val_plaY_el
      real*8, dimension (:,:), allocatable :: val_sism_el
      real*8, dimension (:,:), allocatable :: val_vpcl_el !Clara

	  real*8, dimension (:,:), allocatable :: val_vpsa_el !MCNAP

	  real*8, dimension (:,:), allocatable :: val_vpsl_el !MCNAP

	  real*8, dimension (:,:), allocatable :: val_vpsd_el !MCNAP

	  real*8, dimension (:,:), allocatable :: val_maps_el !Marco
      integer*4, dimension (:), allocatable :: fun_dirX_el,fun_neuX_el

	  integer*4, dimension (:), allocatable :: fun_neuN_el ! M Mexico Arroyo
      integer*4, dimension (:), allocatable :: fun_dirY_el,fun_neuY_el

	  integer*4, dimension (:), allocatable :: fun_neuT_el ! M Mexico Arroyo
      integer*4, dimension (:), allocatable :: fun_poiX_el,fun_poiY_el

	  integer*4, dimension (:), allocatable :: fun_palX_el,fun_palY_el ! Marco

	  integer*4, dimension (:), allocatable :: fun_dipX_el,fun_dipY_el ! Marco
      integer*4, dimension (:), allocatable :: fun_plaX_el,fun_plaY_el
      integer*4, dimension (:), allocatable :: fun_sism_el
	  integer*4, dimension (:), allocatable :: fun_vpcl_el !Clara

	  !integer*4, dimension (:), allocatable :: fun_vpsa_el !MCNAP

	  !integer*4, dimension (:), allocatable :: fun_vpsl_el !MCNAP

	  !integer*4, dimension (:), allocatable :: fun_vpsd_el !MCNAP

	  integer*4, dimension (:), allocatable :: fun_maps_el !Marco
      integer*4, dimension(:), allocatable :: tag_dirX_el,tag_neuX_el,tag_plaX_el

	  integer*4, dimension(:), allocatable :: tag_neuN_el ! M Mexico Arroyo
      integer*4, dimension(:), allocatable :: tag_dirY_el,tag_neuY_el,tag_plaY_el

	  integer*4, dimension(:), allocatable :: tag_neuT_el ! M Mexico Arroyo

	  integer*4, dimension(:), allocatable :: tag_palX_el,tag_palY_el ! M Mexico _A_
      integer*4, dimension(:), allocatable :: tag_sism_el
	  integer*4, dimension(:), allocatable :: tag_vpcl_el !Clara

	  integer*4, dimension(:), allocatable :: tag_vpsa_el !MCNAP

	  integer*4, dimension(:), allocatable :: tag_vpsl_el !MCNAP

      integer*4, dimension(:), allocatable :: tag_vpsd_el !MCNAP

	  integer*4, dimension(:), allocatable :: tag_maps_el !Marco
      integer*4, dimension(:), allocatable :: tag_abc_el

	  integer*4, dimension(:), allocatable :: inode_dipX_el,inode_dipY_el ! Marco
      integer*4 :: nload_dirX_el,nload_dirY_el
      integer*4 :: nload_neuX_el,nload_neuY_el

	  integer*4 :: nload_neuN_el,nload_neuT_el ! M Mexico Arroyo
      integer*4 :: nload_poiX_el,nload_poiY_el

	  integer*4 :: nload_palX_el,nload_palY_el ! Marco

	  integer*4 :: nload_dipX_el,nload_dipY_el ! Marco
      integer*4 :: nload_plaX_el,nload_plaY_el
      integer*4 :: nload_sism_el
      integer*4 :: nload_vpcl_el !Clara

	  integer*4 :: nload_vpsa_el !MCNAP

	  integer*4 :: nload_vpsl_el !MCNAP

	  integer*4 :: nload_vpsd_el !MCNAP

	  integer*4 :: nload_maps_el !Marco
      integer*4 :: nload_abc_el
      
      integer*4, dimension (:), allocatable :: tag_func
      integer*4, dimension (:), allocatable :: func_type
      integer*4, dimension (:), allocatable :: func_indx
      real*8, dimension (:), allocatable :: func_data
      integer*4 :: nfunc,nfunc_data
      
      real*8, dimension (:), allocatable :: tsnap

	  real*8, dimension (:), allocatable :: tsnapps ! Marco
      real*8, dimension (:), allocatable :: x_monitor,y_monitor
      
      integer*4, dimension (:),   allocatable :: itersnap

	  integer*4, dimension (:),   allocatable :: itersnapps ! Marco
      integer*4, dimension (:),   allocatable :: n_monitor
      
      integer*4 :: ns,nn,nn2
      integer*4 :: im,ie,i,j,in,ic,id
      integer*4 :: nts,nsnaps,nmonitors

	  integer*4 :: nsnapsps ! marco
      integer*4 :: err_out,trash
      
      character*70 :: head_file,grid_file,mat_file,out_file
      
      real*8 :: eps
      
      integer*4 :: myid,nnode_dom,nprocs
      integer*4, dimension(:), allocatable :: node_index
     
 
!   SEISMIC MOMENT VARIABLE

      real*8, dimension(:,:), allocatable :: facsmom
      integer*4, dimension(:), allocatable :: num_node_sism
      integer*4, dimension(:,:), allocatable :: sour_node_sism  
      real*8, dimension(:,:), allocatable :: dist_sour_node_sism
      integer*4, dimension(:,:), allocatable :: check_node_sism
      real*8, dimension(:,:), allocatable :: check_dist_node_sism
      integer*4 :: max_num_node_sism,length_check_node_sism
      integer*4 :: conta
 

!   DAMPING MATRIX VARIABLES
      integer*4 :: make_damping_yes_or_not

      real*8 :: ndt_mon



!   PRE-STRESS VARIABLES

!	  real*8, dimension(:,:,:), allocatable :: sxx_ps,syy_ps ! Marco 
      
!**********************************************************
      
      call system_clock(COUNT=start,COUNT_RATE=clock(2))
      
      err_out = 0
      myid = 0
      nprocs = 1
      
      head_file = 'else2_input.d'
      
      if (myid.eq.0) then
         write(*,'(A)')''
         write(*,'(A)')''
         write(*,'(A)')'****************************************'
         write(*,'(A)')'*                                      *'
         write(*,'(A)')'*  ****** **      ****  ******  ****   *'
         write(*,'(A)')'*  ****** **     ****** ****** ******  *'
         write(*,'(A)')'*  **     **     **  ** **     **  **  *'
         write(*,'(A)')'*  **     **     **     **         **  *'
         write(*,'(A)')'*  *****  **     *****  ******  *****  *'
         write(*,'(A)')'*  *****  **      ***** ****** *****   *'
         write(*,'(A)')'*  **     **         ** **     **      *'
         write(*,'(A)')'*  **     **     **  ** **     **      *'
         write(*,'(A)')'*  ****** ****** ****** ****** ******  *'
         write(*,'(A)')'*  ****** ******  ****  ****** ******  *'
         write(*,'(A)')'*                                      *'
         write(*,'(A)')'****************************************'
         write(*,'(A)')''
         write(*,'(A)')''
         write(*,'(A)')'****************************************'
         write(*,'(A)')'*                                      *'
         write(*,'(A)')'* SPECTRAL ELEMENT CODE FOR SOLVING 2D *'
         write(*,'(A)')'* ELASTICITY PROBLEMS                  *'
         write(*,'(A)')'*                                      *'
         write(*,'(A)')'*  © CRS4, 2002, All Rights Reserved   *'
         write(*,'(A)')'*                                      *'
         write(*,'(A)')'****************************************'
         write(*,'(A)')''
         write(*,'(A)')''
      endif
      
      if (myid.eq.0) write(*,'(A)')
      if (myid.eq.0) write(*,'(A,A20)')'Reading Header file : ',head_file
      
      call check_file(head_file,err_out)
      if (err_out.ne.0) then
         stop
      endif
      
      opt_out_form = 1
      opt_out_data = 1
      
      call read_dime_header(head_file,nmonitors,nsnaps,nsnapsps) ! Marco
      
      if (nmonitors.gt.0) then
         allocate(x_monitor(nmonitors),y_monitor(nmonitors))
         allocate(n_monitor(nmonitors))
      endif
      
      if (nsnaps.gt.0) then

        allocate(tsnap(nsnaps),itersnap(nsnaps))

      endif



	  if (nsnapsps.gt.0) then

        allocate(tsnapps(nsnapsps),itersnapps(nsnapsps))

      endif
      
      call read_header(head_file,grid_file,mat_file,out_file,&
                       ns,deltat,xtime,opt_out_data,opt_out_form,&
                       nmonitors,x_monitor,y_monitor,&
                       nsnaps,tsnap,ndt_mon,&

					   nsnapsps,tsnapps) ! Marco
      
      if (myid.eq.0) write(*,'(A)')'Read.'
      
      mat_file = mat_file(1:len_trim(mat_file)) // '.mat'
      
      if (myid.eq.0) write(*,'(A)')
      if (myid.eq.0) write(*,'(A,A20)')'Reading Material file : ',mat_file
      
      call check_file(mat_file,err_out)
      if (err_out.ne.0) then
         stop
      endif
      
      call read_dime_mat_el(mat_file,nmat, &

							nmatg, &  ! Mexico Paco grad_lin 29.11.2004
                            nload_dirX_el,nload_dirY_el, &
                            nload_neuX_el,nload_neuY_el, &

							nload_neuN_el,nload_neuT_el, & ! M Mexico Arroyo
                            nload_poiX_el,nload_poiY_el, &

							nload_palX_el,& ! Marco

							nload_palY_el,& ! Marco

							nload_dipX_el,nload_dipY_el, & ! Marco
                            nload_plaX_el,nload_plaY_el, &
                            nload_sism_el, &
							nload_vpcl_el, &  !Clara

							nload_vpsa_el, &  !MCNAP

							nload_vpsl_el, &  !MCNAP

							nload_vpsd_el, &  !MCNAP

							nload_maps_el, &  !Marco
                            nload_abc_el, &
                            nfunc,nfunc_data)
      
      if (myid.eq.0) write(*,'(A,I8)')'Materials ........: ',nmat

	  if (myid.eq.0) write(*,'(A,I8)')'Materials ........: ',nmatg ! Mexico Paco grad_lin 29.11.2004
      if (myid.eq.0) write(*,'(A,I8)')'Dichlet X B.C.....: ',nload_dirX_el
      if (myid.eq.0) write(*,'(A,I8)')'Dichlet Y B.C.....: ',nload_dirY_el

	  if (myid.eq.0) write(*,'(A,I8)')'DichletP X B.C....: ',nload_dipX_el ! Marco

      if (myid.eq.0) write(*,'(A,I8)')'DichletP Y B.C....: ',nload_dipY_el ! Marco
      if (myid.eq.0) write(*,'(A,I8)')'Neumann X B.C.....: ',nload_neuX_el
      if (myid.eq.0) write(*,'(A,I8)')'Neumann Y B.C.....: ',nload_neuY_el

	  if (myid.eq.0) write(*,'(A,I8)')'Neumann N B.C.....: ',nload_neuN_el ! M Mexico Arroyo

      if (myid.eq.0) write(*,'(A,I8)')'Neumann T B.C.....: ',nload_neuT_el ! M Mexico Arroyo
      if (myid.eq.0) write(*,'(A,I8)')'Absorbing B.C.....: ',nload_abc_el
      if (myid.eq.0) write(*,'(A,I8)')'Point Loads X ....: ',nload_poiX_el
      if (myid.eq.0) write(*,'(A,I8)')'Point Loads Y ....: ',nload_poiY_el

	  if (myid.eq.0) write(*,'(A,I8)')'Point Loads All X : ',nload_palX_el ! Marco

	  if (myid.eq.0) write(*,'(A,I8)')'Point Loads All Y : ',nload_palY_el ! Marco
      if (myid.eq.0) write(*,'(A,I8)')'Plane Loads X ....: ',nload_plaX_el
      if (myid.eq.0) write(*,'(A,I8)')'Plane Loads Y ....: ',nload_plaY_el
	  if (myid.eq.0) write(*,'(A,I8)')'Materials VP Clay.: ',nload_vpcl_el !Clara

	  if (myid.eq.0) write(*,'(A,I8)')'Materials VP Sand.: ',nload_vpsa_el !MCNAP

	  if (myid.eq.0) write(*,'(A,I8)')'MaterialsPS ......: ',nload_maps_el !Marco
      if (myid.eq.0) write(*,'(A,I8)')'Moment Loads .....: ',nload_sism_el
      if (myid.eq.0) write(*,'(A,I8)')'Functions ........: ',nfunc
      
      if (nmat.le.0) then
         write(*,*)'Error ! nmat = 0'
         stop
      endif
      
      allocate (type_mat(nmat))
      allocate (tref_mat(nmat))
      allocate (prop_mat(nmat,4))
      allocate (sdeg_mat(nmat))
      allocate (tag_mat(nmat))





        ! MATG -> MATerial with linear Gradient

	    !

        ! Properties are stored as follows:

		! 1 - rho1

		! 2 - lambda1

		! 3 - mu1

		! 4 - gamma1

        ! 5 - x1

	    ! 6 - y1

		! 7 - rho2

		! 8 - lambda2

		! 9 - mu2

		! 10 - gamma2

        ! 11 - x2

	    ! 12 - y2

		!  

		! Right now implementation allows only a linear gradient of velocity along

		! y direction; it seems to be easy to generalize at least for a linear gradient

		! in a generic 2D direction (this is the reason whay it was decided to store

		! both coordinates of point P1 (x1,y1) and point P2 (x2,y2))

		!

		!                  P1 (x1,y1)

		!

        !   +--------------*---------------+.. \  Vs1 Vp1

		!   |                              |    \ 

		!   |                              |     \ 

		!   |                              |      \

		!   |                              |       \

		!   |                              |        \

		!   |                              |         \

		!   |                              |          \

		!   +--------------*---------------+...........\ Vs2 Vp2 

		!   

		!                   P2 (x2,y2)

		!

		!

	  if (nmatg.gt.0) then ! Mexico Paco grad_lin 29.11.2004

		allocate (type_matg(nmatg)) ! Mexico Paco grad_lin 29.11.2004

		!allocate (tref_matg(nmatg)) ! Mexico Paco grad_lin 29.11.2004

		allocate (prop_matg(nmatg,18)) ! Mexico Paco grad_lin 29.11.2004

		!allocate (sdeg_matg(nmatg)) ! Mexico Paco grad_lin 29.11.2004

		allocate (tag_matg(nmatg)) ! Mexico Paco grad_lin 29.11.2004

	  endif ! Mexico Paco grad_lin 29.11.2004
      
      
      if (nload_dirX_el.gt.0) then
         allocate (val_dirX_el(nload_dirX_el,2))
         allocate (fun_dirX_el(nload_dirX_el))
         allocate (tag_dirX_el(nload_dirX_el))
      endif
      
      if (nload_dirY_el.gt.0) then
         allocate (val_dirY_el(nload_dirY_el,2))
         allocate (fun_dirY_el(nload_dirY_el))
         allocate (tag_dirY_el(nload_dirY_el))
      endif
      
      if (nload_neuX_el.gt.0) then
         allocate (val_neuX_el(nload_neuX_el,2))
         allocate (fun_neuX_el(nload_neuX_el))
         allocate (tag_neuX_el(nload_neuX_el))
      endif
      
      if (nload_neuY_el.gt.0) then
         allocate (val_neuY_el(nload_neuY_el,2))
         allocate (fun_neuY_el(nload_neuY_el))
         allocate (tag_neuY_el(nload_neuY_el))
      endif



      if (nload_neuN_el.gt.0) then ! M Mexico Arroyo

         allocate (val_neuN_el(nload_neuN_el,2)) ! M Mexico Arroyo

         allocate (fun_neuN_el(nload_neuN_el)) ! M Mexico Arroyo

         allocate (tag_neuN_el(nload_neuN_el)) ! M Mexico Arroyo

      endif ! M Mexico Arroyo

      

      if (nload_neuT_el.gt.0) then ! M Mexico Arroyo

         allocate (val_neuT_el(nload_neuT_el,2)) ! M Mexico Arroyo

         allocate (fun_neuT_el(nload_neuT_el)) ! M Mexico Arroyo

         allocate (tag_neuT_el(nload_neuT_el)) ! M Mexico Arroyo

      endif ! M Mexico Arroyo
      
      if (nload_poiX_el.gt.0) then

         allocate (val_poiX_el(nload_poiX_el,3))

         allocate (fun_poiX_el(nload_poiX_el))

      endif



      if (nload_poiY_el.gt.0) then

         allocate (val_poiY_el(nload_poiY_el,3))

         allocate (fun_poiY_el(nload_poiY_el))

      endif



	  if (nload_palX_el.gt.0) then ! Marco

         allocate (val_palX_el(nload_palX_el,1)) ! Marco

         allocate (fun_palX_el(nload_palX_el)) ! Marco

		 allocate (tag_palX_el(nload_palX_el)) ! M Mexico _A_

      endif ! Marco


	  if (nload_palY_el.gt.0) then ! Marco

         allocate (val_palY_el(nload_palY_el,1)) ! Marco

         allocate (fun_palY_el(nload_palY_el)) ! Marco

		 allocate (tag_palY_el(nload_palY_el)) ! M Mexico _A_

      endif ! Marco



      if (nload_dipX_el.gt.0) then ! Marco

         allocate (val_dipX_el(nload_dipX_el,3)) ! Marco

         allocate (fun_dipX_el(nload_dipX_el)) ! Marco

		 allocate (inode_dipX_el(nload_dipX_el)) ! Marco

      endif ! Marco

      

      if (nload_dipY_el.gt.0) then ! Marco

         allocate (val_dipY_el(nload_dipY_el,3)) ! Marco

         allocate (fun_dipY_el(nload_dipY_el)) ! Marco

		 allocate (inode_dipY_el(nload_dipX_el)) ! Marco

      endif ! Marco


      if (nload_plaX_el.gt.0) then
         allocate (val_plaX_el(nload_plaX_el,1))
         allocate (fun_plaX_el(nload_plaX_el))
         allocate (tag_plaX_el(nload_plaX_el))
      endif
      
      if (nload_plaY_el.gt.0) then
         allocate (val_plaY_el(nload_plaY_el,1))
         allocate (fun_plaY_el(nload_plaY_el))
         allocate (tag_plaY_el(nload_plaY_el))
      endif

      if (nload_sism_el.gt.0) then
         allocate (val_sism_el(nload_sism_el,12))
         allocate (fun_sism_el(nload_sism_el))
         allocate (tag_sism_el(nload_sism_el))  
      endif

	  if (nload_vpcl_el.gt.0) then !Clara
         allocate (val_vpcl_el(nload_vpcl_el,6)) !Clara
         allocate (fun_vpcl_el(nload_vpcl_el)) !Clara
         allocate (tag_vpcl_el(nload_vpcl_el))  !Clara 
      endif !Clara



	  if (nload_vpsa_el.gt.0) then !MCNAP

         allocate (val_vpsa_el(nload_vpsa_el,5)) !MCNAP

         !allocate (fun_vpsa_el(nload_vpsa_el)) !MCNAP

         allocate (tag_vpsa_el(nload_vpsa_el))  !MCNAP 

      endif !MCNAP



	  if (nload_vpsl_el.gt.0) then !MCNAP

         allocate (val_vpsl_el(nload_vpsl_el,16)) !MCNAP

         !allocate (fun_vpsl_el(nload_vpsl_el)) !MCNAP

         allocate (tag_vpsl_el(nload_vpsl_el))  !MCNAP 

      endif !MCNAP



	  if (nload_vpsd_el.gt.0) then !MCNAP

         allocate (val_vpsd_el(nload_vpsd_el,16)) !MCNAP

         !allocate (fun_vpsd_el(nload_vpsd_el)) !MCNAP

         allocate (tag_vpsd_el(nload_vpsd_el))  !MCNAP 

      endif !MCNAP



	  if (nload_maps_el.gt.0) then !Marco

         allocate (val_maps_el(nload_maps_el,1)) !Marco

         allocate (fun_maps_el(nload_maps_el)) !Marco

         allocate (tag_maps_el(nload_maps_el))  !Marco 

      endif !Marco

      
      if (nload_abc_el.gt.0) then
         allocate (tag_abc_el(nload_abc_el))
      endif
      
      if (nfunc.gt.0) then
         allocate (tag_func(nfunc))
         allocate (func_type(nfunc))
         allocate (func_indx(nfunc +1))
         allocate (func_data(nfunc_data))
      endif
      
      
      call read_material_el(mat_file,nmat,prop_mat,type_mat,tref_mat,tag_mat,&

				nmatg,prop_matg,type_matg,tag_matg, & ! Mexico Paco grad_lin 29.11.2004
                nload_dirX_el,val_dirX_el,fun_dirX_el,tag_dirX_el, &
                nload_dirY_el,val_dirY_el,fun_dirY_el,tag_dirY_el, &
                nload_neuX_el,val_neuX_el,fun_neuX_el,tag_neuX_el, &
                nload_neuY_el,val_neuY_el,fun_neuY_el,tag_neuY_el, &

				nload_neuN_el,val_neuN_el,fun_neuN_el,tag_neuN_el, & ! M Mexico Arroyo

                nload_neuT_el,val_neuT_el,fun_neuT_el,tag_neuT_el, & ! M Mexico Arroyo
                nload_poiX_el,val_poiX_el,fun_poiX_el, &
                nload_poiY_el,val_poiY_el,fun_poiY_el, &

				nload_palX_el,val_palX_el,fun_palX_el,tag_palX_el, & ! Marco ! M Mexico _A_

				nload_palY_el,val_palY_el,fun_palY_el,tag_palY_el, & ! Marco ! M Mexico _A_

				nload_dipX_el,val_dipX_el,fun_dipX_el, & ! Marco

                nload_dipY_el,val_dipY_el,fun_dipY_el, & ! Marco
                nload_plaX_el,val_plaX_el,fun_plaX_el,tag_plaX_el, &
                nload_plaY_el,val_plaY_el,fun_plaY_el,tag_plaY_el, &
		        nload_sism_el,val_sism_el,fun_sism_el,tag_sism_el, &
				nload_vpcl_el,val_vpcl_el,fun_vpcl_el,tag_vpcl_el, & !Clara

				nload_vpsa_el,val_vpsa_el,tag_vpsa_el, & !MCNAP

				nload_vpsl_el,val_vpsl_el,tag_vpsl_el, & !MCNAP

				nload_vpsd_el,val_vpsd_el,tag_vpsd_el, & !MCNAP

				nload_maps_el,val_maps_el,fun_maps_el,tag_maps_el, & !Marco
                nload_abc_el,tag_abc_el, &
                nfunc,func_type,func_indx,func_data,nfunc_data,tag_func, &
                fmax, &
		        err_out)
      
      do i = 1,nmat
         sdeg_mat(i) = ns
      enddo
      
      if (myid.eq.0) write(*,'(A)')'Read.'
      
            
      if (myid.eq.0) then
         write(*,*)
         do im = 1,nmat
            write(*,'(A,I8)') 'MATERIAL : ',tag_mat(im)
            write(*,'(A,I8)') 'TYPE : ',type_mat(im)
            write(*,'(A,E12.4)') 'rho : ',prop_mat(im,1)
            write(*,'(A,E12.4)') 'Vp : ',((prop_mat(im,2)+ &
                                         2*prop_mat(im,3))/prop_mat(im,1))**0.5
            write(*,'(A,E12.4)') 'Vs : ',(prop_mat(im,3)/prop_mat(im,1))**0.5
            write(*,'(A,E12.4)') 'gamma : ',prop_mat(im,4)
            write(*,*)
         enddo
      endif



	  if (myid.eq.0) then ! Mexico Paco grad_lin 29.11.2004

         write(*,*) ! Mexico Paco grad_lin 29.11.2004

         do im = 1,nmatg ! Mexico Paco grad_lin 29.11.2004

            write(*,'(A,I8)') 'MATERIAL : ',tag_matg(im) ! Mexico Paco grad_lin 29.11.2004

            write(*,'(A,I8)') 'TYPE : ',type_matg(im) ! Mexico Paco grad_lin 29.11.2004

            write(*,'(A,E12.4)') 'rho1 : ',prop_matg(im,1) ! Mexico Paco grad_lin 29.11.2004

            write(*,'(A,E12.4)') 'Vp1 : ',((prop_matg(im,2)+ & ! Mexico Paco grad_lin 29.11.2004

                                         2*prop_matg(im,3))/prop_matg(im,1))**0.5 ! Mexico Paco grad_lin 29.11.2004

            write(*,'(A,E12.4)') 'Vs1: ',(prop_matg(im,3)/prop_matg(im,1))**0.5 ! Mexico Paco grad_lin 29.11.2004

            write(*,'(A,E12.4)') 'gamma1 : ',prop_matg(im,4) ! Mexico Paco grad_lin 29.11.2004

			write(*,'(A,E12.4)') 'x1 : ',prop_matg(im,5) ! Mexico Paco grad_lin 29.11.2004

			write(*,'(A,E12.4)') 'y1 : ',prop_matg(im,6) ! Mexico Paco grad_lin 29.11.2004

            write(*,*) ! Mexico Paco grad_lin 29.11.2004

			write(*,'(A,E12.4)') 'rho2 : ',prop_matg(im,7) ! Mexico Paco grad_lin 29.11.2004

            write(*,'(A,E12.4)') 'Vp2 : ',((prop_matg(im,8)+ & ! Mexico Paco grad_lin 29.11.2004

                                         2*prop_matg(im,9))/prop_matg(im,7))**0.5 ! Mexico Paco grad_lin 29.11.2004

            write(*,'(A,E12.4)') 'Vs2: ',(prop_matg(im,9)/prop_matg(im,7))**0.5 ! Mexico Paco grad_lin 29.11.2004

            write(*,'(A,E12.4)') 'gamma2 : ',prop_matg(im,10) ! Mexico Paco grad_lin 29.11.2004

			write(*,'(A,E12.4)') 'x2 : ',prop_matg(im,11) ! Mexico Paco grad_lin 29.11.2004

			write(*,'(A,E12.4)') 'y2 : ',prop_matg(im,12) ! Mexico Paco grad_lin 29.11.2004

			write(*,*) ! Mexico Paco grad_lin 29.11.2004

         enddo ! Mexico Paco grad_lin 29.11.2004

      endif ! Mexico Paco grad_lin 29.11.2004

      !Clara begin
	  if (myid.eq.0) then
         write(*,*)
         do im = 1,nload_vpcl_el
            write(*,'(A,I8)') 'MATERIAL VP Clay........: ',tag_vpcl_el(im)
            write(*,'(A,I8)') 'TYPE.... ...............: ',fun_vpcl_el(im)
            write(*,'(A,E12.4)') 'Initial Yield Stress : ',val_vpcl_el(im,1)
            write(*,'(A,E12.4)') 'Softening modulus ...: ',val_vpcl_el(im,2)
            write(*,'(A,E12.4)') 'gamma (Perzyna) .....: ',val_vpcl_el(im,3)
            write(*,'(A,E12.4)') 'N (Perzyna) .........: ',val_vpcl_el(im,4)
			write(*,'(A,E12.4)') 'friction angle (DEG).: ',val_vpcl_el(im,5)

			write(*,'(A,E12.4)') 'friction angle (DEG).: ',val_vpcl_el(im,6)
            write(*,*)
         enddo
      endif
      !Clara end



	  !MCNAP begin

	  if (myid.eq.0) then

         write(*,*)

         do im = 1,nload_vpsl_el

            write(*,'(A,I8)') 'MATERIAL VP Sand .......: ',tag_vpsa_el(im)

            !write(*,'(A,I8)') 'TYPE.... ..............: ',fun_vpsa_el(im)

            write(*,'(A,E12.4)') 'densità rel..Dr(%)...: ',val_vpsa_el(im,1)

            write(*,'(A,E12.4)') 'XXXXXXXXXXXXXXXXXXXX.: ',val_vpsa_el(im,2)

            write(*,'(A,E12.4)') 'XXXXXXXXXXXXXXXXXXXX.: ',val_vpsa_el(im,3)

            write(*,'(A,E12.4)') 'XXXXXXXXXXXXXXXXXXXX.: ',val_vpsa_el(im,4)

			write(*,'(A,E12.4)') 'XXXXXXXXXXXXXXXXXXXX.: ',val_vpsa_el(im,5)

            write(*,*)

         enddo

      endif

      !MCNAP end



	  !MCNAP begin

	  if (myid.eq.0) then

         write(*,*)

         do im = 1,nload_vpsl_el

            write(*,'(A,I8)') 'MATERIAL VP Sand Loose..: ',tag_vpsl_el(im)

            !write(*,'(A,I8)') 'TYPE.... ..............: ',fun_vpsl_el(im)

            write(*,'(A,E12.4)') 'ga_l............gamma: ',val_vpsl_el(im,1)

            write(*,'(A,E12.4)') 'la_l...............Bp: ',val_vpsl_el(im,2)

            write(*,'(A,E12.4)') 'bef0_l.......beta_f_0: ',val_vpsl_el(im,3)

            write(*,'(A,E12.4)') 'csi_l...dilat_compres: ',val_vpsl_el(im,4)

			write(*,'(A,E12.4)') 'psi_l....dilat_estens: ',val_vpsl_el(im,5)

			write(*,'(A,E12.4)') 'gw_l...............tp: ',val_vpsl_el(im,6)

			write(*,'(A,E12.4)') 'po_l..............rc0: ',val_vpsl_el(im,7)

			write(*,'(A,E12.4)') 'pq_l.... teta_compres: ',val_vpsl_el(im,8)

			write(*,'(A,E12.4)') 'pw_l...............cp: ',val_vpsl_el(im,9)

			write(*,'(A,E12.4)') 'bm_l......teta_estens: ',val_vpsl_el(im,10)

			write(*,'(A,E12.4)') 'befl_l.....beta_f_lim: ',val_vpsl_el(im,11)

			write(*,'(A,E12.4)') 'gam_l..g_visco[kPa/s]: ',val_vpsl_el(im,12)

			write(*,'(A,E12.4)') 'gav_l......alfa_visco: ',val_vpsl_el(im,13)

			write(*,'(A,E12.4)') 'emax.................: ',val_vpsl_el(im,14)

			write(*,'(A,E12.4)') 'XXXXXXXXXXXXXXXXXXXX.: ',val_vpsl_el(im,15)

			write(*,'(A,E12.4)') 'XXXXXXXXXXXXXXXXXXXX.: ',val_vpsl_el(im,16)

            write(*,*)

         enddo

      endif

      !MCNAP end

	 

	  !MCNAP begin

	  if (myid.eq.0) then

         write(*,*)

         do im = 1,nload_vpsd_el

            write(*,'(A,I8)') 'MATERIAL VP Sand Dense..: ',tag_vpsd_el(im)

            !write(*,'(A,I8)') 'TYPE.... ..............: ',fun_vpsd_el(im)

            write(*,'(A,E12.4)') 'ga_d............gamma: ',val_vpsd_el(im,1)

            write(*,'(A,E12.4)') 'la_d...............Bp: ',val_vpsd_el(im,2)

            write(*,'(A,E12.4)') 'bef0_d.......beta_f_0: ',val_vpsd_el(im,3)

            write(*,'(A,E12.4)') 'csi_d...dilat_compres: ',val_vpsd_el(im,4)

			write(*,'(A,E12.4)') 'psi_d....dilat_estens: ',val_vpsd_el(im,5)

			write(*,'(A,E12.4)') 'gw_d...............tp: ',val_vpsd_el(im,6)

			write(*,'(A,E12.4)') 'po_d..............rc0: ',val_vpsd_el(im,7)

			write(*,'(A,E12.4)') 'pq_d.... teta_compres: ',val_vpsd_el(im,8)

			write(*,'(A,E12.4)') 'pw_d...............cp: ',val_vpsd_el(im,9)

			write(*,'(A,E12.4)') 'bm_d......teta_estens: ',val_vpsd_el(im,10)

			write(*,'(A,E12.4)') 'befl_d.....beta_f_lim: ',val_vpsd_el(im,11)

			write(*,'(A,E12.4)') 'gam_d..g_visco[kPa/s]: ',val_vpsd_el(im,12)

			write(*,'(A,E12.4)') 'gav_d......alfa_visco: ',val_vpsd_el(im,13)

			write(*,'(A,E12.4)') 'emin.................: ',val_vpsd_el(im,14)

			write(*,'(A,E12.4)') 'XXXXXXXXXXXXXXXXXXXX.: ',val_vpsd_el(im,15)

			write(*,'(A,E12.4)') 'XXXXXXXXXXXXXXXXXXXX.: ',val_vpsd_el(im,16)

            write(*,*)

         enddo

      endif

      !MCNAP end



	  !Marco begin



!     Set initial pre-stress condition

      

      !allocate (sxx_ps(nn,nn,nelem),syy_ps(nn,nn,nelem))

      

	  !do ie=1,nelem

	  !	do i = 1,nn

	  !		do j = 1,nn

	  !			sxx_ps(i,j,ie) = 0.0d0

	  !			syy_ps(i,j,ie) = 0.0d0

	  !		enddo

	  !	enddo

	  !enddo



	  if (myid.eq.0) then

         write(*,*)

         do im = 1,nload_maps_el

            write(*,'(A,I8)')    'MATERIAL ...............: ',tag_maps_el(im)

            write(*,'(A,I8)')    'TYPE.... ...............: ',fun_maps_el(im)

            write(*,'(A,E12.4)') 'Prestressed with phi :    ',val_maps_el(im,1)

            write(*,*)



            !do ie = 1,nelem

			!	do j = 1,nmat

			!		if (tag_mat(j).eq.tag_maps_el(im)) then

			!			!which_maps = im ! Clara

			!			call make_pre_stress(prop_mat(im,1),prop_mat(im,2),prop_mat(im,3),&

			!								 val_maps_el(im,1),&

			!								 sxx_ps,syy_ps,&

			!								 nn,ie,nelem,&

			!					             nnod,xx_spx,yy_spx,&

			!								 con_nnz,con_spx)

			!	    endif

			!	enddo

			!enddo



         enddo

      endif

      !Marco end
      
      grid_file = grid_file(1:len_trim(grid_file)) // '.inp'
      
      if (myid.eq.0) write(*,'(A)')
      if (myid.eq.0) write(*,'(A,A20)')'Reading Grid file : ',grid_file
      
      call check_file(grid_file,err_out)
      if (err_out.ne.0) then
         stop
      endif
      
      call read_dime_grid_el(grid_file,nmat,tag_mat,&
                     nload_dirX_el,tag_dirX_el,nload_dirY_el,tag_dirY_el,&
                     nload_neuX_el,tag_neuX_el,nload_neuY_el,tag_neuY_el,&

					 nload_neuN_el,tag_neuN_el,nload_neuT_el,tag_neuT_el,& ! M Mexico Arroyo
                     nload_abc_el,tag_abc_el,&
                     nnod_macro,nelem,nedge)
      
      if (myid.eq.0) write(*,'(A,I8)')'Nodes : ',nnod_macro
      if (myid.eq.0) write(*,'(A,I8)')'Elements : ',nelem
      if (myid.eq.0) write(*,'(A,I8)')'Edges : ',nedge
      
      if (nnod_macro.gt.0) then
         allocate (xx_macro(nnod_macro),yy_macro(nnod_macro))
      else
         write(*,*)'Error ! nnod_macro = 0'
         stop
      endif
      
      if (nelem.gt.0) then
         allocate (con(nelem,5))
      else
         write(*,*)'Error ! nelem = 0'
         stop
      endif
      
      if (nedge.gt.0) allocate (con_bc(nedge,3))
      
      call read_grid_el(grid_file,nmat,tag_mat,prop_mat,&
                        nload_dirX_el,tag_dirX_el,nload_dirY_el,tag_dirY_el,&
                        nload_neuX_el,tag_neuX_el,nload_neuY_el,tag_neuY_el,&

						nload_neuN_el,tag_neuN_el,nload_neuT_el,tag_neuT_el,& ! M Mexico Arroyo
                        nload_abc_el,tag_abc_el,&
                        nnod_macro,xx_macro,yy_macro,nelem,con,nedge,con_bc)
      
      if (myid.eq.0) write(*,'(A)')'Read.'
      
      if (myid.eq.0) write(*,'(A)')
      if (myid.eq.0) write(*,'(A)')'Making Spectral connectivities'
      
      
!     Make spectral connectivities for the elements
      
      allocate(Ebw(nnod_macro))
      
      call make_Ebw_macro(nnod_macro,nelem,con,Ebw,Ennz)
      
      allocate(Ebin(0:Ennz))
      
      call make_Ebin_macro(nnod_macro,nelem,con,Ebw,Ennz,Ebin)
      
      deallocate(Ebw)
      
      con_nnz = nelem +1
      do ie = 1,nelem
         do j = 1,nmat
            if (tag_mat(j).eq.con(ie,1)) nn = sdeg_mat(j) +1
         enddo
         con_nnz = con_nnz + nn*nn +1
      enddo
      
      allocate(con_spx(0:con_nnz))
      
      call make_spectral_connectivity(nelem,con,nmat,tag_mat,sdeg_mat,&
                                      Ennz,Ebin,con_nnz,con_spx,nnod)
      
      
      deallocate(Ebin)
      
! Dual connectivity for spectral nodes
      
      allocate(Ebw(nnod))
      
      call make_Ebw(nnod,con_nnz,con_spx,Ebw,Ennz)
      
      allocate(Ebin(0:Ennz))
      
      call make_Ebin(nnod,con_nnz,con_spx,Ebw,Ennz,Ebin)
      
      deallocate(Ebw)
      
      allocate(xx_spx(nnod),yy_spx(nnod))
      allocate(alfa1(nelem),beta1(nelem),gamma1(nelem),delta1(nelem))
      allocate(alfa2(nelem),beta2(nelem),gamma2(nelem),delta2(nelem))
      
      if (myid.eq.0) write(*,'(A,I8)')'Spectral Nodes : ',nnod
      
      call make_spectral_grid(nnod_macro,xx_macro,yy_macro,con_nnz,con_spx,&
                              nmat,tag_mat,sdeg_mat,nelem,&
                              alfa1,beta1,gamma1,delta1,&
                              alfa2,beta2,gamma2,delta2,&
                              nnod,xx_spx,yy_spx)
      
!     Make the spectral connectivities for the boundary
      
      con_nnz_bc = 0
      
      if (nedge.gt.0) then
         con_nnz_bc = nedge +1
         do i = 1,nedge
            call get_edge_element(Ennz,Ebin,&
                          con_bc(i,2),con_bc(i,3),ie)
            do j = 1,nmat
               if (tag_mat(j).eq.con(ie,1)) nn = sdeg_mat(j) +1
            enddo
            con_nnz_bc = con_nnz_bc +nn +1
         enddo
         
         allocate(con_spx_bc(0:con_nnz_bc))
         
         call make_spectral_boundary(con_nnz,con_spx,nedge,con_bc,&
                                     nmat,tag_mat,sdeg_mat,Ennz,Ebin,&
                                     con_nnz_bc,con_spx_bc)
      endif
      
      if (myid.eq.0) write(*,'(A)')'Made.'
      
      
      allocate(node_index(nnod))
      
      do in = 1,nnod
         node_index(in) = in
      enddo
      nnode_dom = nnod
      
!*************************************************************
      
      if (myid.eq.0) write(*,'(A)')
      if (myid.eq.0) write(*,'(A)')
      if (myid.eq.0) write(*,'(A)')'Building the ELASTIC matrices ...'
      
      
!     ELASTIC MATRICES 
      
      Knnz_el = 0
      allocate (Kbin_el(0:Knnz_el))
      
      
      if (myid.eq.0) write(*,'(A)')
      if (myid.eq.0) write(*,'(A)')'Building the inertia matrix...'
      
      allocate (Mel(2*nnode_dom))

      call make_Mel(nnod,node_index,con_nnz,con_spx,&
                    nmat,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                    nelem,alfa1,beta1,gamma1,alfa2,beta2,gamma2,&
                    nnode_dom,Mel)
                    
      if (myid.eq.0) write(*,'(A)')'Inertia matrix built'


      !Check if there is any damping factor between materials characteristics
      
      make_damping_yes_or_not = 0
      
      do im = 1,nmat
	 !write(50,*),im,'| prop_mat =',prop_mat(im,4)
         if (abs(prop_mat(im,4)).gt.10e-5) then
            !write(50,*)'YES!!!'
            make_damping_yes_or_not = 1
         endif
      enddo

      if (make_damping_yes_or_not.eq.1) then
      
         if (myid.eq.0) write(*,'(A)')
         if (myid.eq.0) write(*,'(A)')'Building the damping matrix...'
      
         allocate (Cel(2*nnode_dom))
         allocate (KCel(2*nnode_dom))

         call make_Cel_KCel(nnod,node_index,con_nnz,con_spx,&
                            nmat,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                            nelem,alfa1,beta1,gamma1,alfa2,beta2,gamma2,&
                            nnode_dom,Cel,KCel)
                            
         if (myid.eq.0) write(*,'(A)')'Damping matrix built'
  
      else
      
         if (myid.eq.0) write(*,'(A)')
         if (myid.eq.0) write(*,'(A)')'Damping matrix... NOT BUILT!!!'
         if (myid.eq.0) write(*,'(A)')'There are no materials with damping defined on'
      
      endif
      
      
      
      if (myid.eq.0) write(*,'(A)')
      if (myid.eq.0) write(*,'(A)')'Building the load matrix...'


      ! Dimensioning vector 'num_node_sism'(nodes number generating each single fault)
      
      if (nload_sism_el.gt.0) then

      allocate (num_node_sism(nload_sism_el))

         do i = 1,nload_sism_el
            if (((val_sism_el(i,1).eq.val_sism_el(i,3)).and.(val_sism_el(i,3).eq.val_sism_el(i,5))) &
               .and.((val_sism_el(i,2).eq.val_sism_el(i,4)).and.(val_sism_el(i,4).eq.val_sism_el(i,6))))  then
               
               num_node_sism(i)=1
            
            else  
               call dime_sism_nodes(val_sism_el(i,1),val_sism_el(i,2),val_sism_el(i,3),val_sism_el(i,4),&
                                    val_sism_el(i,5),val_sism_el(i,6),val_sism_el(i,7),val_sism_el(i,8),&
                                    nnod,xx_spx,yy_spx,num_node_sism(i))
            endif
         enddo

         !Checking the maximum number of fault nodes
         max_num_node_sism = num_node_sism(1)
         do i = 1,nload_sism_el
                if (num_node_sism(i).gt.max_num_node_sism) then
                        max_num_node_sism = num_node_sism(i)
                endif
         enddo
      
         allocate (sour_node_sism(max_num_node_sism,nload_sism_el))
         allocate (dist_sour_node_sism(max_num_node_sism,nload_sism_el))
      
         !Searching the node 'id' in the global numeration for each fault.
         !sour_node_sism = node id (global numeration) generating the fault 'i'
         do i = 1,nload_sism_el
                if (num_node_sism(i).eq.1) then
                   call find_nearest_node(nnod,xx_spx,yy_spx,&
                                          val_sism_el(i,1),val_sism_el(i,2),sour_node_sism(1,i))
                      dist_sour_node_sism(1,i) = 0
                
                else 
                   call read_sism_nodes(val_sism_el(i,1),val_sism_el(i,2),val_sism_el(i,3),val_sism_el(i,4),&
                                        val_sism_el(i,5),val_sism_el(i,6),val_sism_el(i,7),val_sism_el(i,8),&
                                        nnod,xx_spx,yy_spx,num_node_sism(i),sour_node_sism,i,&
                                        dist_sour_node_sism,nload_sism_el,&
                                        max_num_node_sism)
                endif
                if (myid.eq.0) then
                   if (i.eq.1) then
                      write(*,'(A)')'Sesmic moment & fault'
                      write(*,'(A)')
                   endif
                   if (num_node_sism(i).eq.1) then
                         write(*,'(A,I6,A)')'Sism ',i,' is located on:'
                         write(*,'(I6,2E14.5,I6)')(j,xx_spx(sour_node_sism(j,i)),&
                               yy_spx(sour_node_sism(j,i)),sour_node_sism(j,i),&
                               j=1,num_node_sism(i))
                   else
                      write(*,'(A,I6,A,I6,A)')'Fault ',i,' is generated by ',num_node_sism(i),' nodes'
                      write(*,'(I6,2E14.5,I6)')(j,xx_spx(sour_node_sism(j,i)),&
                               yy_spx(sour_node_sism(j,i)),sour_node_sism(j,i),&
                               j=1,num_node_sism(i)) 
                   endif
                   write(*,'(A)')
                endif
         enddo


      endif

      if (nfunc.le.0) nfunc = 1
     
         if (nload_sism_el.gt.0) then 
            allocate (facsmom(nload_sism_el,3))
         endif

      allocate (Fel(nfunc,2*nnode_dom))
      
      call make_Fel(nnod,xx_spx,yy_spx,node_index,con_nnz,con_spx,&
                    nmat,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                    nelem,alfa1,beta1,gamma1,alfa2,beta2,gamma2,&
                    con_nnz_bc,con_spx_bc,&
                    nload_dirX_el,val_dirX_el,fun_dirX_el,tag_dirX_el,&
                    nload_dirY_el,val_dirY_el,fun_dirY_el,tag_dirY_el,&
                    nload_neuX_el,val_neuX_el,fun_neuX_el,tag_neuX_el,&
                    nload_neuY_el,val_neuY_el,fun_neuY_el,tag_neuY_el,&

					nload_neuN_el,val_neuN_el,fun_neuN_el,tag_neuN_el,& ! M Mexico Arroyo

                    nload_neuT_el,val_neuT_el,fun_neuT_el,tag_neuT_el,& ! M Mexico Arroyo
                    nload_poiX_el,val_poiX_el,fun_poiX_el,&
                    nload_poiY_el,val_poiY_el,fun_poiY_el,&

					nload_palX_el,val_palX_el,fun_palX_el,tag_palX_el,& ! Marco ! M Mexico _A_

 					nload_palY_el,val_palY_el,fun_palY_el,tag_palY_el,& ! Marco ! M Mexico _A_

					nload_dipX_el,val_dipX_el,fun_dipX_el,inode_dipX_el,& ! Marco

                    nload_dipY_el,val_dipY_el,fun_dipY_el,inode_dipY_el,& ! Marco
                    nload_plaX_el,val_plaX_el,fun_plaX_el,tag_plaX_el,&
                    nload_plaY_el,val_plaY_el,fun_plaY_el,tag_plaY_el,&
                    nload_sism_el,val_sism_el,fun_sism_el,tag_sism_el,&
                    nfunc,tag_func,nnode_dom,Fel,&
                    xx_macro,yy_macro,con,nelem,con_bc,nedge,&
                    num_node_sism,max_num_node_sism,&  
                    sour_node_sism,dist_sour_node_sism,&
                    length_check_node_sism,facsmom)

      if (nload_sism_el.gt.0) then
         allocate (check_node_sism(length_check_node_sism,5)) 
         allocate (check_dist_node_sism(length_check_node_sism,1))
      
         call check_sism(con_nnz,con_spx,&
                         nmat,tag_mat,sdeg_mat,&   
                         nelem,&
                         nload_sism_el,&
                         num_node_sism,max_num_node_sism,&
                         sour_node_sism,dist_sour_node_sism,&
                         check_node_sism,check_dist_node_sism,&
                         length_check_node_sism,&
                         fun_sism_el,nfunc,tag_func,val_sism_el)

         deallocate (sour_node_sism)
         deallocate (dist_sour_node_sism)

         !controllo se l'aver portato fuori la routine funziona
         !MS

                          
         !write(51,*),length_check_node_sism
         !do conta = 1,length_check_node_sism
         !   write(51,*),check_node_sism(conta,1),' | ',&                                   
         !               check_node_sism(conta,2),' | ',&                           
         !               check_node_sism(conta,3),' | ',&                           
         !               check_node_sism(conta,4),' | ',&
         !               check_node_sism(conta,5),' | ',&                           
         !               check_dist_node_sism(conta,1),' | '                           
         !enddo                                     
      
         !MS - fine
      endif


      if (myid.eq.0) write(*,'(A)')'Load matrix built'
      
      
!*************************************************************
!*************************************************************
      
      deallocate (Ebin)
      
      if (myid.eq.0) write(*,'(A)')
      if (myid.eq.0) write(*,'(A)')'Matrices construction terminated.'
      
!******************************************************************
      
      
      if (myid.eq.0) write(*,'(A)')
      if (myid.eq.0) write(*,'(A)')'Setting calculation parameters'
      
      
!     Set initial condition ELASTIC
      
      allocate (u0(2*nnode_dom),u_1(2*nnode_dom))
      
      do in = 1,2*nnode_dom
         u0(in) = 0.0d0
         u_1(in) = 0.0d0
      enddo
      
      
      
      
!     TIME STEPS
      
      if (myid.eq.0) write(*,'(A,E12.4)')'deltat chosen = ',deltat
      call deltat_max(deltat,nnod,nmat,tag_mat,prop_mat,sdeg_mat,&
        xx_macro,yy_macro,nelem,con,deltat_cfl,fmax)
      nts = int(xtime / deltat)
      if (myid.eq.0) write(*,'(A,I8)')'Number of time-steps = ',nts
      xtime = dfloat(nts) * deltat
      if (myid.eq.0) write(*,'(A,E12.4)')'Final time = ',xtime
      
      
!     SNAPSHOTS
      
      if (nsnaps.gt.0) then
        do i = 1,nsnaps
          itersnap(i) = int(tsnap(i)/deltat)
          if (itersnap(i).gt.nts) itersnap(i) = nts
        enddo
      endif





!     SNAPSHOTS PRESTRESS

      

      if (nsnapsps.gt.0) then

        do i = 1,nsnapsps

          itersnapps(i) = int(tsnapps(i)/deltat)

          if (itersnapps(i).gt.nts) itersnapps(i) = nts

        enddo

      endif
      
      
!     NODI MONITORATI
      
      if (nmonitors.ge.1) then
         do i = 1,nmonitors
            call find_nearest_node(nnod,xx_spx,yy_spx,&
                                   x_monitor(i),y_monitor(i),n_monitor(i))
         enddo
      endif
      
      if (myid.eq.0) then
         write(*,'(A)')'Monitored nodes - user coordinates'
         write(*,'(I6,2E14.5,I6)')(i,x_monitor(i),y_monitor(i),n_monitor(i),&
              i=1,nmonitors)
      endif
      
       if (myid.eq.0) then
         write(*,'(A)')'Monitored nodes - spectral coordinates'
         write(*,'(I6,2E14.5,I6)')(i,xx_spx(n_monitor(i)),yy_spx(n_monitor(i)),n_monitor(i),&
              i=1,nmonitors)
      endif

      
      call system_clock(COUNT=finish)
      time_in_seconds = float(finish-start)/float(clock(2))
      
      if (myid.eq.0) then
         write(*,'(A)')
         write(*,'(A)')'**************************************************'
         write(*,'(A,F8.4,A)')'Set-up time = ',time_in_seconds,' s'
         write(*,'(A)')'**************************************************'
         write(*,'(A)')
         
         write(*,'(A)')
         write(*,'(A)')'*************************************************'
         write(*,'(A)')'Beginning of the time-loop...'
         write(*,'(A)')
      endif
      
      
      call time_loop_el(nnod,xx_spx,yy_spx,node_index,con_nnz,con_spx,&
                        nmat,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&

						nmatg,prop_matg,tag_matg,type_matg,& ! Mexico Paco grad_lin 29.11.2004
                        nelem,alfa1,beta1,gamma1,alfa2,beta2,gamma2,&
                        con_nnz_bc,con_spx_bc,&
                        nload_dirX_el,tag_dirX_el,nload_dirY_el,tag_dirY_el,&

						nload_dipX_el,inode_dipX_el,nload_dipY_el,inode_dipY_el,& ! Marco
                        nload_abc_el,tag_abc_el,&
                        nfunc,func_type,func_indx,func_data,nfunc_data,&
                        nnode_dom,Mel,Cel,KCel,Knnz_el,Kbin_el,&
                        Fel,u0,u_1,&
                        nts,deltat,nmonitors,n_monitor,nsnaps,itersnap,&
                        out_file,myid,&
                        check_node_sism,check_dist_node_sism,&
                        length_check_node_sism,facsmom,&
                        nload_sism_el,&
                        make_damping_yes_or_not,ndt_mon,&
						nload_vpcl_el,val_vpcl_el,fun_vpcl_el,tag_vpcl_el,&  !Clara

						nload_vpsa_el,val_vpsa_el,tag_vpsa_el,&  !MCNAP

						nload_vpsl_el,val_vpsl_el,tag_vpsl_el,&  !MCNAP

						nload_vpsd_el,val_vpsd_el,tag_vpsd_el,&  !MCNAP

						nload_maps_el,val_maps_el,fun_maps_el,tag_maps_el,&  ! Marco

						nsnapsps,itersnapps,nprocs)  !Marco
      
      if (myid.eq.0) then
         if (nsnaps.gt.0) then
            write(*,'(A)')
            write(*,'(A)')'Print output'
            
            if (opt_out_form.eq.1) then
               !call write_UCD_output_el(nnod,xx_spx,yy_spx,con_nnz,con_spx,&
               !           nmat,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
               !           nelem,alfa1,alfa2,beta1,beta2,gamma1,gamma2,&
               !           opt_out_data,nsnaps,tsnap,out_file,nprocs)
			   !


               ! DA SISTEMARE!!!!!! 
			   call write_UCD_GID_output_el(nnod,xx_spx,yy_spx,con_nnz,con_spx,& !Clara
                          nmat,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,& !Clara
                          nelem,alfa1,alfa2,beta1,beta2,gamma1,gamma2,& !Clara
                          opt_out_data,nsnaps,tsnap,out_file,nprocs,& !Clara

						  nload_maps_el,& ! Marco

						  nsnapsps,itersnapps,fun_maps_el) ! Marco
               ! DA SISTEMARE!!!!!!

            elseif (opt_out_form.eq.2) then
               do i = 1,nsnaps
                  if ((opt_out_data.eq.1).or.(opt_out_data.eq.2) & 
                       .or.(opt_out_data.eq.3)) then
                     call write_nodal_output_el(nnod,xx_spx,yy_spx,&
                                      con_nnz,con_spx,nmat,tag_mat,sdeg_mat,&
                                      out_file,i,tsnap(i),nprocs)
                  endif
                  if ((opt_out_data.eq.2).or.(opt_out_data.eq.3)) then
                     call write_element_output_el(nnod,xx_spx,yy_spx,&
                          con_nnz,con_spx,&
                          nmat,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                          nelem,alfa1,alfa2,beta1,beta2,gamma1,gamma2,&
                          out_file,i,tsnap(i),nprocs)
                  endif
               enddo
            endif

			

			if (nsnapsps.gt.0) then



			endif


         endif
         
         write(*,'(A)')'**************************************'
         write(*,'(A)')
         write(*,'(A)')'Bye.'
      endif
      
      
      if (nmonitors.gt.0) deallocate(x_monitor,y_monitor,n_monitor)
      if (nsnaps.gt.0) deallocate(tsnap,itersnap)
      
      deallocate (tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat)
      
      if (nload_dirX_el.gt.0) deallocate (val_dirX_el,fun_dirX_el,tag_dirX_el)
      if (nload_dirY_el.gt.0) deallocate (val_dirY_el,fun_dirY_el,tag_dirY_el)
      if (nload_neuX_el.gt.0) deallocate (val_neuX_el,fun_neuX_el,tag_neuX_el)
      if (nload_neuY_el.gt.0) deallocate (val_neuY_el,fun_neuY_el,tag_neuY_el)

	  if (nload_neuN_el.gt.0) deallocate (val_neuN_el,fun_neuN_el,tag_neuN_el) ! M Mexico Arroyo

      if (nload_neuT_el.gt.0) deallocate (val_neuT_el,fun_neuT_el,tag_neuT_el) ! M Mexico Arroyo
      if (nload_poiX_el.gt.0) deallocate (val_poiX_el,fun_poiX_el)
      if (nload_poiY_el.gt.0) deallocate (val_poiY_el,fun_poiY_el)

	  if (nload_palX_el.gt.0) deallocate (val_palX_el,fun_palX_el,tag_palX_el) ! Marco ! M Mexico _A_

	  if (nload_palY_el.gt.0) deallocate (val_palY_el,fun_palY_el,tag_palY_el) ! Marco ! M Mexico _A_

	  if (nload_dipX_el.gt.0) deallocate (val_dipX_el,fun_dipX_el) ! Marco

      if (nload_dipY_el.gt.0) deallocate (val_dipY_el,fun_dipY_el) ! Marco
      if (nfunc.gt.0) deallocate (func_type,func_indx,func_data)
      
      
      deallocate (xx_macro,yy_macro,con)
      
      deallocate (xx_spx,yy_spx,con_spx)
      deallocate (alfa1,beta1,gamma1,delta1,alfa2,beta2,gamma2,delta2)
      
      if (nedge.gt.0) deallocate (con_bc,con_spx_bc)
      
      deallocate (node_index)
      
      deallocate (Kbin_el,Mel)
      deallocate (u0,u_1)
      
      
      end program SEelse2D
      
      
!     *****************************************************************
      
      subroutine legendref (p2,p2der,p1,p1der,n,x)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio
!     
!     COMPUTES THE LEGENDRE POLYNOMIAL (p2) AND ITS DERIVATIVE (p2der)
!     OF DEGREE n AT x AND
!     THE LEGENDRE POLYNOMIAL (p1) AND ITS DERIVATIVE (p1der)
!     OF DEGREE n-1 AT x                                              
                                                                     
      implicit real*8 (a-h,o-z)                            

      p2=1.d0                                                    
      p2der=0.d0                                                    
      if (n.eq.0) return                                          
      p1=p2
      p2=x
      p1der=p2der
      p2der=1.d0
      if (n.eq.1) return                               
      do k=1,n-1
        p0=p1
        p1=p2
        p0der=p1der
        p1der=p2der
        dk=dfloat(k)                                         
        a1=(2.d0*dk+1.d0)/(dk+1.d0)
        a2=-dk/(dk+1.d0)
        p2=a1*p1*x+a2*p0
        p2der=a1*p1+a1*p1der*x+a2*p0der
      enddo
      return                    
      end subroutine legendref


!     *****************************************************************


      subroutine lgl(np,ct,ww,dd)

!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio
!     
      implicit real*8 (a-h,o-z)

      parameter (nstep=1000,acc=1.d-15)

      dimension ct(*),ww(*),dd(np,*)

      n=np-1
      ct(1)=-1.d0
      ct(np)=1.d0
      n2=idint(0.5d0*np)
      dx=2.d0/dfloat(nstep)
      x=-1.d0
      iroot=1
      call legendref (p2,a1,p1,p1der,n,x)
      do while (iroot.lt.n2)
        x=x+dx
        call legendref (p2,a2,p1,p1der,n,x)
        if (dabs(a2).le.acc) then
          iroot=iroot+1
          ct(iroot)=a2
        else
          aa=a1*a2
          if (aa.lt.0.d0) then
            iroot=iroot+1
            ct(iroot)=rtbis(x-dx,x,acc,n)
          endif
        endif
        a1=a2
      enddo

      n_filt=2*n2
      if (n_filt.ne.np) then         ! np odd
        ct(n2+1)=0.d0
        do i=1,n2-1
          ct(np-i)=-ct(i+1)
        enddo
      else                           ! np even
        do i=1,n2-1
          ct(np-i)=-ct(i+1)
        enddo
      endif

      xn=dfloat(n)
      acost=2.d0/(xn*(xn+1.d0))

      do j=1,np
         call legendref (p2,p2der,p1,p1der,n,ct(j))
         den=p2*p2
         ww(j)=acost/den
         do i=1,np
            if (i.ne.j) then
               call legendref (pnum,p2der,p1,p1der,n,ct(i))
               den=p2*(ct(i)-ct(j))
               dd(i,j)=pnum/den
            endif
         enddo
      enddo
      do j=2,np-1
         dd(j,j)=0.d0
      enddo
      dd(1,1)=-0.25d0*xn*(xn+1.d0)
      dd(np,np)=0.25d0*xn*(xn+1.d0)

      return

      end subroutine lgl


!     *****************************************************************


      real*8 function rtbis(x1,x2,xacc,n)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio
!     
      implicit real*8(a-h,o-z)
      
      parameter (jmax=1000)

      call legendref(p2,fmid,p1,p1der,n,x2)
      call legendref(p2,f,p1,p1der,n,x1)
      if (f*fmid.ge.0.d0) then
        write(*,*) 'root non-bracketed !'
        stop
      endif
      if (f.lt.0.d0) then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif

      do j=1,jmax
        dx=0.5d0*dx
        xmid=rtbis+dx
        call legendref(p2,fmid,p1,p1der,n,xmid)
        if (fmid.le.0.d0) rtbis=xmid
        if (dabs(dx).lt.xacc.and.dabs(fmid).le.xacc) then
          rtbis=xmid
          return
        endif
      enddo

      end function rtbis


!     *********************************************************************

      subroutine make_Ebw_macro(nnode,nelem,con,bw,nnz)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Luca Massidda, Javier Sabadell, Giuliana Siddi
!     
      implicit none
      
      integer*4 :: nnode,nelem,nnz
      integer*4, dimension(nelem,5) :: con
      integer*4, dimension(nnode) :: bw
      
      integer*4 :: i,in,ie,nn
      
      nnz = nnode +1
      
      do in = 1,nnode
         bw(in) = 0
      enddo
      
      do ie = 1,nelem
         do i = 1,4
            in = con(ie,i +1)
            
            bw(in) = bw(in) +1
            nnz = nnz +1
         enddo
      enddo
      
      return
      
      end subroutine make_Ebw_macro
      
      
      
!     *********************************************************************
      
      subroutine make_Ebin_macro(nnode,nelem,con,bw,nnz,bin)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Luca Massidda, Javier Sabadell, Giuliana Siddi
!     
      implicit none
      
      integer*4 :: nnode,nelem,nnz
      integer*4, dimension(nelem,5) :: con
      integer*4, dimension(nnode) :: bw
      integer*4, dimension(0:nnz) :: bin
      
      integer*4, dimension(:), allocatable :: ic
      integer*4 :: i,j,k,in,ie,nn2
      integer*4 :: bj
      
      
      allocate(ic(nnode))
      
      bin(0) = nnode +1
      do in = 1,nnode
         bin(in) = bin(in-1) + bw(in)
      enddo
      
      do in = 1,nnode
         ic(in) = bin(in-1)
      enddo
      
      do ie = 1,nelem
         do i = 1,4
            in = con(ie,i +1)
            
            bin(ic(in)) = ie
            ic(in) = ic(in) +1
         enddo
      enddo
      
      deallocate(ic)
      
      return
      
      end subroutine make_Ebin_macro
      
      
      
!     ****************************************************************
      
      subroutine make_spectral_connectivity(nelem,con_mac,nmat,tag_mat,sdeg,&
                                            Ennz,Ebin,con_nnz,con_spx,nnode)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Luca Massidda, Javier Sabadell, Giuliana Siddi
!     
      implicit none
      
      integer*4 :: nelem,nmat,Ennz,con_nnz,nnode
      integer*4, dimension(nelem,5) :: con_mac
      integer*4, dimension(nmat) :: tag_mat
      integer*4, dimension(nmat) :: sdeg
      integer*4, dimension(0:Ennz) :: Ebin
      integer*4, dimension(0:con_nnz) :: con_spx
      
      integer*4 :: nnode_mac,nn,imat,ie,i,j,k,an,bn,check
      
!     Makes spectral connectivity for counter clockwise oriented elements
      
      
      con_spx(0) = nelem +1
      do ie = 1,nelem
         do imat = 1,nmat
            if (con_mac(ie,1).eq.tag_mat(imat)) then
               nn = sdeg(imat) +1
               con_spx(ie) = con_spx(ie -1) + nn*nn +1
            endif
         enddo
      enddo
      
      nnode_mac = Ebin(0) -1
      
      nnode = nnode_mac
      
      do ie = 1,nelem
         do imat = 1,nmat
            if (con_mac(ie,1).eq.tag_mat(imat)) then
               nn = sdeg(imat) +1
               
! First put material index
               con_spx(con_spx(ie -1) +0) = con_mac(ie,1)
               
! Then put vertices
               con_spx(con_spx(ie -1) +1) = con_mac(ie,2)
               con_spx(con_spx(ie -1) + nn) = con_mac(ie,3)
               con_spx(con_spx(ie -1) + nn*(nn-1) +1) = con_mac(ie,5)
               con_spx(con_spx(ie -1) + nn*nn) = con_mac(ie,4)
               
! Then construct edge connectivity
               
! First edge (down)
               
               an = con_spx(con_spx(ie -1) +1)
               bn = con_spx(con_spx(ie -1) +nn)
               
               check = 0
               do i = Ebin(an -1),Ebin(an) -1
                  do j = Ebin(bn -1),Ebin(bn) -1
                     if (Ebin(i).eq.Ebin(j)) then
                        if ((Ebin(i).lt.ie).and.((nn*nn +1)&
                        .eq.(con_spx(Ebin(i)) - con_spx(Ebin(i) -1)))) then
                           
                           if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) + k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(k -1) +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) + k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn -k +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*(nn -1) +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) + k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -1) +k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) + k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -k +1))
                              enddo
                              check = 1
                              
                           endif
                        endif
                     endif
                  enddo
               enddo
               
               if (check.eq.0) then
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     con_spx(con_spx(ie -1) + k) = nnode
                  enddo
               endif
               
               
! Fourth edge (left)
               
               an = con_spx(con_spx(ie -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*(nn -1) +1)
               
               check = 0
               do i = Ebin(an -1),Ebin(an) -1
                  do j = Ebin(bn -1),Ebin(bn) -1
                     if (Ebin(i).eq.Ebin(j)) then
                        if ((Ebin(i).lt.ie).and.((nn*nn +1)&
                        .eq.(con_spx(Ebin(i)) - con_spx(Ebin(i) -1)))) then
                           
                           if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(k -1) +1) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(k -1) +1) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*(nn -1) +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(k -1) +1) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -k) +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(k -1) +1) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*nn -k +1)
                              enddo
                              check = 1
                              
                           endif
                        endif
                     endif
                  enddo
               enddo
               
               if (check.eq.0) then
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     con_spx(con_spx(ie -1) +nn*(k -1) +1) = nnode
                  enddo
               endif
               
               
! Third edge (up)
               
               an = con_spx(con_spx(ie -1) +nn*(nn -1) +1)
               bn = con_spx(con_spx(ie -1) +nn*nn)
               
               check = 0
               do i = Ebin(an -1),Ebin(an) -1
                  do j = Ebin(bn -1),Ebin(bn) -1
                     if (Ebin(i).eq.Ebin(j)) then
                        if ((Ebin(i).lt.ie).and.((nn*nn +1)&
                        .eq.(con_spx(Ebin(i)) - con_spx(Ebin(i) -1)))) then
                           
                           if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(nn -1) +k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(nn -1) +k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*(nn -1) +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(nn -1) +k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -k) +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*(nn -1) +k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*nn -k +1)
                              enddo
                              check = 1
                              
                           endif
                        endif
                     endif
                  enddo
               enddo
               
               if (check.eq.0) then
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     con_spx(con_spx(ie -1) +nn*(nn -1) +k) = nnode
                  enddo
               endif
               
               
! Second edge (rigth)
               
               an = con_spx(con_spx(ie -1) +nn)
               bn = con_spx(con_spx(ie -1) +nn*nn)
               
               check = 0
               do i = Ebin(an -1),Ebin(an) -1
                  do j = Ebin(bn -1),Ebin(bn) -1
                     if (Ebin(i).eq.Ebin(j)) then
                        if ((Ebin(i).lt.ie).and.((nn*nn +1)&
                        .eq.(con_spx(Ebin(i)) - con_spx(Ebin(i) -1)))) then
                           
                           if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(k -1) +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn -k +1)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*(nn -1) +1)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -1) +k)
                              enddo
                              check = 1
                              
                           else if (an.eq.con_spx(con_spx(Ebin(i) -1) &
                                +nn*nn)) then
                              do k = 2,(nn -1)
                                 con_spx(con_spx(ie -1) +nn*k) = &
                                      con_spx(con_spx(Ebin(i) -1) &
                                      +nn*(nn -k +1))
                              enddo
                              check = 1
                              
                           endif
                        endif
                     endif
                  enddo
               enddo
               
               if (check.eq.0) then
                  do k = 2,(nn -1)
                     nnode = nnode +1
                     con_spx(con_spx(ie -1) +nn*k) = nnode
                  enddo
               endif
               
               
! End of edge connectivity
               
! Finally fill inside the elements
               
               do j = 2,(nn -1)
                  do i = 2,(nn -1)
                     nnode = nnode +1
                     con_spx(con_spx(ie -1) +nn*(j -1) +i) = nnode
                  enddo
               enddo
            endif
         enddo
      enddo
      
      return
      
      end subroutine make_spectral_connectivity
      
      
      
!     ******************************************************
      
      subroutine make_spectral_grid(nn_m,xx_m,yy_m,cs_nnz,cs,nm,tm,sd,ne,&
                                    alfa1,beta1,gamma1,delta1,&
                                    alfa2,beta2,gamma2,delta2,&
                                    nn_s,xx_s,yy_s)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Luca Massidda, Javier Sabadell, Giuliana Siddi
!     
      implicit none
      
      integer*4 :: nn_m,cs_nnz,nm,ne,nn_s
      real*8, dimension(nn_m) :: xx_m,yy_m
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tm
      integer*4, dimension(nm) :: sd
      real*8, dimension(ne) :: alfa1,beta1,gamma1,delta1
      real*8, dimension(ne) :: alfa2,beta2,gamma2,delta2
      real*8, dimension(nn_s) :: xx_s,yy_s
      
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8 :: xp,yp,x1,x2,x3,x4,y1,y2,y3,y4,jac
      integer*4 :: im,ie,i,j,nn,ip
      
      
      nn = 2
      allocate(ct(nn),ww(nn),dd(nn,nn))
      call lgl(nn,ct,ww,dd)
      
      do im = 1,nm
         if ((sd(im) +1).ne.nn) then
            deallocate(ct,ww,dd)
            nn = sd(im) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call lgl(nn,ct,ww,dd)
         endif
         
         do ie = 1,ne
            if (cs(cs(ie -1) +0).eq.tm(im)) then
               x1 = xx_m(cs(cs(ie -1) +1))
               x2 = xx_m(cs(cs(ie -1) +nn))
               x3 = xx_m(cs(cs(ie -1) +nn*nn))
               x4 = xx_m(cs(cs(ie -1) +nn*(nn -1) +1))
               
               y1 = yy_m(cs(cs(ie -1) +1))
               y2 = yy_m(cs(cs(ie -1) +nn))
               y3 = yy_m(cs(cs(ie -1) +nn*nn))
               y4 = yy_m(cs(cs(ie -1) +nn*(nn -1) +1))
               
               alfa1(ie) = 0.25d0*(-x1 +x2 +x3 -x4)
               beta1(ie) = 0.25d0*(-x1 -x2 +x3 +x4)
               gamma1(ie) = 0.25d0*(+x1 -x2 +x3 -x4)
               delta1(ie) = 0.25d0*(+x1 +x2 +x3 +x4)
               
               alfa2(ie) = 0.25d0*(-y1 +y2 +y3 -y4)
               beta2(ie) = 0.25d0*(-y1 -y2 +y3 +y4)
               gamma2(ie) = 0.25d0*(+y1 -y2 +y3 -y4)
               delta2(ie) = 0.25d0*(+y1 +y2 +y3 +y4)
               
               jac = alfa1(ie)*beta2(ie) - alfa2(ie)*beta1(ie)
               
               if (jac.le.0.d0) then
                  write(*,*)'Error ! Orientation non-conforming !'
                  write(*,*)'(element ',ie,' is clockwise oriented)'
                  stop
               endif
               
               do j = 1,nn
                  do i = 1,nn
                     xp = alfa1(ie) * ct(i) + beta1(ie) * ct(j) &
                        + gamma1(ie) * ct(i) * ct(j) + delta1(ie)
                     yp = alfa2(ie) * ct(i) + beta2(ie) * ct(j) &
                        + gamma2(ie) * ct(i) * ct(j) + delta2(ie)
                     
                     ip = nn*(j -1) +i
                     xx_s(cs(cs(ie -1) + ip)) = xp
                     yy_s(cs(cs(ie -1) + ip)) = yp
                  enddo
               enddo
            endif
         enddo
      enddo
      
      return
      
      end subroutine make_spectral_grid
      
      
      
!     ****************************************************************
      
      subroutine make_spectral_boundary(cs_nnz,cs,ne_bc,cm_bc,nm,tm,sd,&
                                        Ennz,Ebin,cs_nnz_bc,cs_bc)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      
      implicit none
      
      integer*4 :: cs_nnz,ne_bc,nm,Ennz,cs_nnz_bc
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(ne_bc,*) :: cm_bc
      integer*4, dimension(nm) :: tm 
      integer*4, dimension(nm) :: sd 
      integer*4, dimension(0:Ennz) :: Ebin
      integer*4, dimension(0:cs_nnz_bc) :: cs_bc
      
      integer*4 :: im,ie,ielem,i,j,k,nn
      integer*4 :: an,bn,n1,n2,n3,n4
      
      
      cs_bc(0) = ne_bc +1
      
      do ie = 1,ne_bc
         an = cm_bc(ie,1 +1)
         bn = cm_bc(ie,2 +1)
         
         call get_edge_element(Ennz,Ebin,&
                               cm_bc(ie,2),cm_bc(ie,3),ielem)
         do im = 1,nm
            if (tm(im).eq.cs(cs(ielem -1) +0)) nn = sd(im) +1
         enddo
         cs_bc(ie) = cs_bc(ie -1) + nn +1
      enddo
      
      
      do ie = 1,ne_bc
         an = cm_bc(ie,1 +1)
         bn = cm_bc(ie,2 +1)
         
         cs_bc(cs_bc(ie -1) +0) = cm_bc(ie,1)
         
         call get_edge_element(Ennz,Ebin,&
                               cm_bc(ie,2),cm_bc(ie,3),ielem)
         do im = 1,nm
            if (tm(im).eq.cs(cs(ielem -1) +0)) nn = sd(im) +1
         enddo
         
         n1 = cs(cs(ielem -1) +1)
         n2 = cs(cs(ielem -1) +nn)
         n3 = cs(cs(ielem -1) +nn*nn)
         n4 = cs(cs(ielem -1) +nn*(nn -1) +1)
         
         if (((an.eq.n1).and.(bn.eq.n2)) &
              .or.((bn.eq.n1).and.(an.eq.n2))) then
            do k = 1,nn
               cs_bc(cs_bc(ie -1) + k) = &
                    cs(cs(ielem -1) +k)
            enddo
            
         else if (((an.eq.n2).and.(bn.eq.n3)) &
              .or.((bn.eq.n2).and.(an.eq.n3))) then
            do k = 1,nn
               cs_bc(cs_bc(ie -1) + k) = &
                    cs(cs(ielem -1) +nn*k)
            enddo
            
         else if (((an.eq.n3).and.(bn.eq.n4)) &
              .or.((bn.eq.n3).and.(an.eq.n4))) then
            do k = 1,nn
               cs_bc(cs_bc(ie -1) + k) = &
                    cs(cs(ielem -1) +nn*(nn -1) +k)
            enddo
            
         else if (((an.eq.n4).and.(bn.eq.n1)) &
              .or.((bn.eq.n4).and.(an.eq.n1))) then
            do k = 1,nn
               cs_bc(cs_bc(ie -1) + k) = &
                    cs(cs(ielem -1) +nn*(k -1) +1)
            enddo
            
         endif
         
      enddo
      
      return
      
      end subroutine make_spectral_boundary
      
      
      
!     *********************************************************************

      subroutine make_Ebw(nnode,cs_nnz,cs,bw,nnz)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Luca Massidda, Javier Sabadell, Giuliana Siddi
!     
      implicit none
      
      integer*4 :: nnode,cs_nnz,nnz
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nnode) :: bw
      
      integer*4 :: i,in,ie,ne,nn2
      
      nnz = nnode +1
      
      do in = 1,nnode
         bw(in) = 0
      enddo
      
      ne = cs(0) -1
      do ie = 1,ne
         nn2 = cs(ie) - cs(ie -1) -1
         do i = 1,nn2
            in = cs(cs(ie -1) +i)
            
            bw(in) = bw(in) +1
            nnz = nnz +1
         enddo
      enddo
      
      return
      
      end subroutine make_Ebw
      
      
      
!     *********************************************************************
      
      subroutine make_Ebin(nnode,cs_nnz,cs,bw,nnz,bin)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Luca Massidda, Javier Sabadell, Giuliana Siddi
!     
      implicit none
      
      integer*4 :: nnode,cs_nnz,nnz
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nnode) :: bw
      integer*4, dimension(0:nnz) :: bin
      
      integer*4, dimension(:), allocatable :: ic
      integer*4 :: i,j,k,in,ie,ne,nn2
      integer*4 :: bj
      
      
      allocate(ic(nnode))
      
      bin(0) = nnode +1
      do in = 1,nnode
         bin(in) = bin(in -1) + bw(in)
      enddo
      
      do in = 1,nnode
         ic(in) = bin(in -1)
      enddo
      
      ne = cs(0) -1
      do ie = 1,ne
         nn2 = cs(ie) - cs(ie -1) -1
         do i = 1,nn2
            in = cs(cs(ie -1) +i)
            
            bin(ic(in)) = ie
            ic(in) = ic(in) +1
         enddo
      enddo
      
      deallocate(ic)
      
      return
      
      end subroutine make_Ebin
      
      
      
!     ***********************************************************************
      
      subroutine make_Mel(nnode,node_index,cs_nnz,cs,&
                          nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                          ne,alfa1,beta1,gamma1,alfa2,beta2,gamma2,&
                          nnode_dom,mvec)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Luca Massidda, Javier Sabadell, Giuliana Siddi
!     

      implicit none
      
      integer*4 :: nnode,cs_nnz,nm,ne,nnode_dom
      integer*4, dimension(nnode) :: node_index
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat
      real*8, dimension(nm) :: tref_mat
      real*8, dimension(nm,4) :: prop_mat
      real*8, dimension(ne) :: alfa1,beta1,gamma1
      real*8, dimension(ne) :: alfa2,beta2,gamma2
      real*8, dimension(2*nnode_dom) :: mvec
      
      real*8, dimension(:), allocatable :: dxdy,dydy,dxdx,dydx
      real*8, dimension(:,:), allocatable :: det_j
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8 :: rho,term
      integer*4 :: im,ie,i,j,is,in,id1,id2,nn,nn2
      
      
      nn = 2
      
      allocate(ct(nn),ww(nn),dd(nn,nn))
      call lgl(nn,ct,ww,dd)
      allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))
      
      do i = 1,2*nnode_dom
         mvec(i) = 0.d0
      enddo
      
      ne = cs(0) -1
      
      do im = 1,nm
         rho = prop_mat(im,1)
         
         if ((sdeg_mat(im) +1).ne.nn) then
            deallocate(ct,ww,dd)
            deallocate(dxdx,dydx,dxdy,dydy,det_j)
            nn = sdeg_mat(im) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call lgl(nn,ct,ww,dd)
            allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))
         endif
         
         do ie = 1,ne
            if (cs(cs(ie -1) +0).eq.tag_mat(im)) then
               do i = 1,nn
                  dxdy(i) = beta1(ie) + gamma1(ie) * ct(i)
                  dydy(i) = beta2(ie) + gamma2(ie) * ct(i)
               enddo
               
               do j = 1,nn
                  dxdx(j) = alfa1(ie) + gamma1(ie) * ct(j)
                  dydx(j) = alfa2(ie) + gamma2(ie) * ct(j)
               enddo
               
               do j = 1,nn
                  do i = 1,nn
                     det_j(i,j) = dxdx(j)*dydy(i) - dxdy(i)*dydx(j)
                  enddo
               enddo
               
               do j = 1,nn
                  do i = 1,nn
                     is = nn*(j -1) +i
                     in = cs(cs(ie -1) +is)
                     
                     if (node_index(in).ne.0) then
                        id1 = node_index(in)
                        id2 = node_index(in) +nnode_dom
                        
                        term = rho * det_j(i,j) * ww(i) * ww(j)
                        mvec(id1) = mvec(id1) + term
                        mvec(id2) = mvec(id2) + term
                     endif
                  enddo
               enddo
            endif
         enddo
         
      enddo
      
      
!     Diagonal inertia matrix has been built (vector mvec(:))
      
      return
      end subroutine make_Mel
      
!     ***********************************************************************
      
      subroutine make_Cel_KCel(nnode,node_index,cs_nnz,cs,&
                               nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                               ne,alfa1,beta1,gamma1,alfa2,beta2,gamma2,&
                               nnode_dom,cvec,kcvec)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Marco Stupazzini
!     

      implicit none
      
      integer*4 :: nnode,cs_nnz,nm,ne,nnode_dom
      integer*4, dimension(nnode) :: node_index
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat
      real*8, dimension(nm) :: tref_mat
      real*8, dimension(nm,4) :: prop_mat
      real*8, dimension(ne) :: alfa1,beta1,gamma1
      real*8, dimension(ne) :: alfa2,beta2,gamma2
      real*8, dimension(2*nnode_dom) :: cvec
      real*8, dimension(2*nnode_dom) :: kcvec
      
      real*8, dimension(:), allocatable :: dxdy,dydy,dxdx,dydx
      real*8, dimension(:,:), allocatable :: det_j
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8 :: rho,term_Cel,term_KCel,gamma
      integer*4 :: im,ie,i,j,is,in,id1,id2,nn,nn2
      
      
      nn = 2
      
      allocate(ct(nn),ww(nn),dd(nn,nn))
      call lgl(nn,ct,ww,dd)
      allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))
      
      do i = 1,2*nnode_dom
         cvec(i) = 0.d0
         kcvec(i) = 0.d0
      enddo
      
      ne = cs(0) -1
      
      do im = 1,nm
         rho = prop_mat(im,1)
         gamma = prop_mat(im,4)
         
         if ((sdeg_mat(im) +1).ne.nn) then
            deallocate(ct,ww,dd)
            deallocate(dxdx,dydx,dxdy,dydy,det_j)
            nn = sdeg_mat(im) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call lgl(nn,ct,ww,dd)
            allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))
         endif
         
         do ie = 1,ne
            if (cs(cs(ie -1) +0).eq.tag_mat(im)) then
               do i = 1,nn
                  dxdy(i) = beta1(ie) + gamma1(ie) * ct(i)
                  dydy(i) = beta2(ie) + gamma2(ie) * ct(i)
               enddo
               
               do j = 1,nn
                  dxdx(j) = alfa1(ie) + gamma1(ie) * ct(j)
                  dydx(j) = alfa2(ie) + gamma2(ie) * ct(j)
               enddo
               
               do j = 1,nn
                  do i = 1,nn
                     det_j(i,j) = dxdx(j)*dydy(i) - dxdy(i)*dydx(j)
                  enddo
               enddo
               
               do j = 1,nn
                  do i = 1,nn
                     is = nn*(j -1) +i
                     in = cs(cs(ie -1) +is)
                     
                     if (node_index(in).ne.0) then
                        id1 = node_index(in)
                        id2 = node_index(in) +nnode_dom
                        
                        term_Cel = 2 * gamma * rho &
                                   * det_j(i,j) * ww(i) * ww(j)
                        cvec(id1) = cvec(id1) + term_Cel
                        cvec(id2) = cvec(id2) + term_Cel
                                                
                        term_KCel = (gamma**2) * rho &
                                    * det_j(i,j) * ww(i) * ww(j)
                        kcvec(id1) = kcvec(id1) + term_KCel
                        kcvec(id2) = kcvec(id2) + term_KCel
                        
                     endif
                  enddo
               enddo
            endif
         enddo
         
      enddo
      
      
!     Diagonal damping  matrix has been built (vector cvec(:) and kcvec(:))
      
      return
      end subroutine make_Cel_KCel
            
      
! *********************************************************
      
      subroutine make_Fel(nnode,xs,ys,node_index,cs_nnz,cs,&
                          nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                          ne,a1,b1,g1,a2,b2,g2,&
                          cs_nnz_bc,cs_bc,&
                          nl_dirX,val_dirX,fun_dirX,tag_dirX,&
                          nl_dirY,val_dirY,fun_dirY,tag_dirY,&
                          nl_neuX,val_neuX,fun_neuX,tag_neuX,&
                          nl_neuY,val_neuY,fun_neuY,tag_neuY,&

						  nl_neuN,val_neuN,fun_neuN,tag_neuN,& ! M Mexico Arroyo

                          nl_neuT,val_neuT,fun_neuT,tag_neuT,& ! M Mexico Arroyo
                          nl_poiX,val_poiX,fun_poiX,&
                          nl_poiY,val_poiY,fun_poiY,&

						  nl_palX,val_palX,fun_palX,tag_palX,& ! Marco

						  nl_palY,val_palY,fun_palY,tag_palY,& ! Marco

						  nl_dipX,val_dipX,fun_dipX,node_dipX,& ! Marco

                          nl_dipY,val_dipY,fun_dipY,node_dipY,& ! Marco
                          nl_plaX,val_plaX,fun_plaX,tag_plaX,&
                          nl_plaY,val_plaY,fun_plaY,tag_plaY,&
                          nl_sism,val_sism,fun_sism,tag_sism,&
                          nfunc,tag_func,nnode_dom,fmat,&
                          xx,yy,con_quad,nquad,con_line,nline,&
                          num_node_sism,max_num_node_sism,&
                          sour_node_sism,dist_sour_node_sism,&
                          length_cns,facsmom)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      implicit none
      
      integer*4 :: nnode,cs_nnz,nm,ne,cs_nnz_bc,nquad,nline
      integer*4 :: nl_dirX,nl_dirY,nl_neuX,nl_neuY

	  integer*4 :: nl_neuN,nl_neuT ! M Mexico Arroyo
      integer*4 :: nl_poiX,nl_poiY,nl_plaX,nl_plaY,nl_sism

	  integer*4 :: nl_palX,nl_palY ! Marco

	  integer*4 :: nl_dipX,nl_dipY ! Marco
      integer*4 :: nfunc,nnode_dom
      real*8, dimension(nnode) :: xs,ys
      integer*4, dimension(nnode) :: node_index
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat
      real*8, dimension(nm) :: tref_mat
      real*8, dimension(nm,4) :: prop_mat
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

	  real*8, dimension(nl_neuN,2) :: val_neuN ! M Mexico Arroyo

      real*8, dimension(nl_neuT,2) :: val_neuT ! M Mexico Arroyo
      integer*4, dimension(nl_neuX) :: fun_neuX
      integer*4, dimension(nl_neuY) :: fun_neuY

	  integer*4, dimension(nl_neuN) :: fun_neuN ! M Mexico Arroyo

      integer*4, dimension(nl_neuT) :: fun_neuT ! M Mexico Arroyo
      integer*4, dimension(nl_neuX) :: tag_neuX
      integer*4, dimension(nl_neuY) :: tag_neuY

	  integer*4, dimension(nl_neuN) :: tag_neuN ! M Mexico Arroyo

      integer*4, dimension(nl_neuT) :: tag_neuT ! M Mexico Arroyo
      real*8, dimension(nl_plaX,1) :: val_plaX
      real*8, dimension(nl_plaY,1) :: val_plaY
      integer*4, dimension(nl_plaX) :: fun_plaX
      integer*4, dimension(nl_plaY) :: fun_plaY
      integer*4, dimension(nl_plaX) :: tag_plaX
      integer*4, dimension(nl_plaY) :: tag_plaY
      real*8, dimension(nl_poiX,3) :: val_poiX
      real*8, dimension(nl_poiY,3) :: val_poiY

	  real*8, dimension(nl_palX,1) :: val_palX ! Marco

      real*8, dimension(nl_palY,1) :: val_palY ! Marco
      integer*4, dimension(nl_poiX) :: fun_poiX
      integer*4, dimension(nl_poiY) :: fun_poiY

	  integer*4, dimension(nl_palX) :: fun_palX ! Marco

      integer*4, dimension(nl_palY) :: fun_palY ! Marco

	  integer*4, dimension(nl_palX) :: tag_palX ! M Mexico _A_

	  integer*4, dimension(nl_palY) :: tag_palY ! M Mexico _A_

	  real*8, dimension(nl_dipX,3) :: val_dipX ! Marco

      real*8, dimension(nl_dipY,3) :: val_dipY ! Marco

      integer*4, dimension(nl_dipX) :: fun_dipX ! Marco

      integer*4, dimension(nl_dipY) :: fun_dipY ! Marco
      real*8, dimension(nl_sism,12) :: val_sism
      integer*4, dimension(nl_sism) :: fun_sism
      integer*4, dimension(nl_sism) :: tag_sism
      integer*4, dimension(nfunc) :: tag_func
      real*8, dimension(nl_sism,3) :: facsmom
      real*8, dimension(nl_sism) :: term_vet
      real*8, dimension(nfunc,2*nnode_dom) :: fmat
      
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8, dimension(:), allocatable :: dxdx,dydx,dxdy,dydy
      real*8, dimension(:,:), allocatable :: det_j
      integer*4, dimension(nl_poiX) :: node_poiX
      integer*4, dimension(nl_poiY) :: node_poiY

	  integer*4, dimension(nl_palX) :: node_palX ! Marco

      integer*4, dimension(nl_palY) :: node_palY ! Marco

	  integer*4, dimension(nl_dipX) :: node_dipX ! Marco

      integer*4, dimension(nl_dipY) :: node_dipY ! Marco

      integer*4, dimension(nl_sism) :: num_node_sism
      integer*4 :: max_num_node_sism
      integer*4, dimension(max_num_node_sism,nl_sism) :: sour_node_sism
      real*8, dimension(max_num_node_sism,nl_sism) :: dist_sour_node_sism
      integer*4 :: length_cns

      integer*4 :: im,ifn,ie,ip,idp,ipl,isism,il,nedge,nn,fn ! Marco

	  integer*4 :: ineuN,ineuT ! M Mexico Arroyo
      integer*4 :: is,in,id1,id2
      integer*4 :: i,j,k,h,l,m
      real*8 :: rho,lambda,mu
      real*8 :: lx,ly,ll,v1,v2,v,term
      real*8 :: rp2,t1x,t1y,t2x,t2y, ellex
      real*8, dimension(nnode) :: xx,yy
      integer*4, dimension(nquad,5) :: con_quad
      integer*4, dimension(nline,3) :: con_line
      integer*4 :: C,sit
      integer*4, dimension(4) :: num_sit
      integer*4, dimension(8) :: coord_sit
      real*8 :: slip1,slip2,norm1,norm2,amp_sism



	  integer*4, dimension(:), allocatable :: i4count ! M Mexico Arroyo

      integer*4, dimension(:), allocatable :: ielem_N,iedge_N ! M Mexico Arroyo

	  integer*4, dimension(:), allocatable :: ielem_T,iedge_T ! M Mexico Arroyo

	  integer*4, dimension(:), allocatable :: ielem_ebe ! M Mexico Arroyo

	  integer*4, dimension(:), allocatable :: inode_ebe

	  integer*4, dimension(:), allocatable :: node_ind_ebe ! M Mexico Arroyo

	  real*8, dimension(:), allocatable :: edge_nx_neuN,edge_ny_neuN ! M Mexico Arroyo

	  real*8, dimension(:), allocatable :: edge_nx_neuT,edge_ny_neuT ! M Mexico Arroyo

	  integer*4 :: nelem_N,nedge_N ! M Mexico Arroyo

	  integer*4 :: nelem_T,nedge_T ! M Mexico Arroyo

	  integer*4 :: iedge,ied1,ied2,iel1,iel2,iel3,iel4 ! M Mexico Arroyo

	  integer*4 :: nelem_ebe,nnode_ebe ! M Mexico Arroyo

	  real*8 :: edge_lx,edge_ly,edge_ll,edge_nx,edge_ny ! M Mexico Arroyo

	  !integer*4 :: ie

      !integer*4 :: edge_ia,edge_ja,edge_ib,edge_jb

      !real*8 :: edge_lx,edge_ly,edge_ll,edge_nx,edge_ny



	  !allocate(i4count(nnode))

      nn = 2
      
      allocate(ct(nn),ww(nn),dd(nn,nn))
      call lgl(nn,ct,ww,dd)
      allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))
      
      num_sit = (/1,2,2,1/)
      coord_sit = (/3,2,5,4,4,3,2,5/)
      length_cns = 0







	  do fn = 1,nfunc

         do id1 = 1,nnode_dom

            Fmat(fn,id1) = 0.0d0

         enddo

         do id2 = nnode_dom+1,2*nnode_dom

            Fmat(fn,id2) = 0.0d0

         enddo

      enddo





	 ! *************************************************

	 ! *************************************************



     ! M Mexico Arroyo



	 if (nl_neuN.gt.0) then



		do im = 1,nm

			if ((sdeg_mat(im) +1).ne.nn) then

				deallocate(ct,ww,dd)

				deallocate(dxdx,dydx,dxdy,dydy,det_j)

            

				nn = sdeg_mat(im) +1

				allocate(ct(nn),ww(nn),dd(nn,nn))

				call lgl(nn,ct,ww,dd)

				allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))

			endif

		enddo



		!if (cs_nnz_bc.gt.0) then

		!	nedge = cs_bc(0) -1

        ! 

		!	do ie = 1,nedge

		!		if ((cs_bc(ie) - cs_bc(ie -1) -1).ne.nn) then

		!			deallocate(ct,ww,dd)

		!			nn = cs_bc(ie) - cs_bc(ie -1) -1

		!			allocate(ct(nn),ww(nn),dd(nn,nn))

		!			call lgl(nn,ct,ww,dd)

		!		endif

		!	enddo

        !

		!endif





		allocate(i4count(ne))

      

		nelem_ebe = 0 

		do ie = 1,ne

			i4count(ie) = 0

		enddo

      

		do im = 1,nm

			if (type_mat(im).eq.2) then

				do ie = 1,ne

					if (cs(cs(ie -1) +0).eq.tag_mat(im)) then

						do is = 1,(cs(ie) - cs(ie -1) -1)

							in = cs(cs(ie -1) +is)

							i4count(ie) = i4count(ie) + node_index(in)

						enddo

					endif

				enddo

			endif

		enddo

      

		do ie = 1,ne

			if (i4count(ie).gt.0) then

				nelem_ebe = nelem_ebe +1

				i4count(ie) = nelem_ebe

			endif

		enddo

      

		if (nelem_ebe.gt.0) then

			allocate(ielem_ebe(nelem_ebe))

         

			do ie = 1,ne

				if (i4count(ie).ne.0) then

				ielem_ebe(i4count(ie)) = ie

				endif

			enddo

		endif

      

		deallocate(i4count)









		nedge_N = 0 

      

		if (cs_nnz_bc.gt.0) then

			nedge = cs_bc(0) -1

         

			allocate(i4count(nedge))

				do iedge = 1,nedge

					i4count(iedge) = 0

				enddo

         

			if (nl_neuN.gt.0) then

				do i = 1,nl_neuN

					do iedge = 1,nedge

						if (cs_bc(cs_bc(iedge -1) +0).eq.tag_neuN(i)) then

							ied1 = cs_bc(cs_bc(iedge -1) +1)

							ied2 = cs_bc(cs_bc(iedge) -1)

                     

							do j = 1,nelem_ebe

								ie = ielem_ebe(j)

                        

								nn = cs_bc(iedge) - cs_bc(iedge -1) -1

                        

								iel1 = cs(cs(ie -1) +1)

								iel2 = cs(cs(ie -1) +nn)

								iel3 = cs(cs(ie -1) +nn*nn)

								iel4 = cs(cs(ie -1) +nn*(nn -1) +1)

                        

								if (((ied1.eq.iel1).or.(ied1.eq.iel2).or. &

									(ied1.eq.iel3).or.(ied1.eq.iel4)).and. &

									((ied2.eq.iel1).or.(ied2.eq.iel2).or. &

									(ied2.eq.iel3).or.(ied2.eq.iel4))) then

									i4count(iedge) = iedge

								endif

							enddo

						endif

					enddo

				enddo

			endif

         

			do iedge = 1,nedge

				if (i4count(iedge).gt.0) then

					nedge_N = nedge_N +1

					i4count(iedge) = nedge_N

				endif

			enddo

      

			if (nedge_N.gt.0) then

				allocate(iedge_N(nedge_N))

            

				do iedge = 1,nedge

					if (i4count(iedge).ne.0) then

						iedge_N(i4count(iedge)) = iedge

					endif

				enddo

			endif

         

			deallocate(i4count)

		endif

      

      

		nelem_N = nedge_N

      

		if (nelem_N.gt.0) then

			allocate(ielem_N(nelem_N))

         

			do i = 1,nedge_N

				iedge = iedge_N(i)

            

				ied1 = cs_bc(cs_bc(iedge -1) +1)

				ied2 = cs_bc(cs_bc(iedge) -1)

            

				do j = 1,nelem_ebe

					ie = ielem_ebe(j)

               

					nn = cs_bc(iedge) - cs_bc(iedge -1) -1

               

					iel1 = cs(cs(ie -1) +1)

					iel2 = cs(cs(ie -1) +nn)

					iel3 = cs(cs(ie -1) +nn*nn)

					iel4 = cs(cs(ie -1) +nn*(nn -1) +1)

               

					if (((ied1.eq.iel1).or.(ied1.eq.iel2).or. &

						(ied1.eq.iel3).or.(ied1.eq.iel4)).and. &

						((ied2.eq.iel1).or.(ied2.eq.iel2).or. &

						(ied2.eq.iel3).or.(ied2.eq.iel4))) then

						ielem_N(i) = ie

					endif

				enddo

			enddo

		endif ! M Mexico Arroyo





		allocate(node_ind_ebe(nnode))

      

		nnode_ebe = 0

		do in = 1,nnode

			node_ind_ebe(in) = 0

		enddo

      

		do i = 1,nelem_ebe

			ie = ielem_ebe(i)

			do is = 1,(cs(ie) - cs(ie -1) -1)

				in = cs(cs(ie -1) +is)

				node_ind_ebe(in) = 1

			enddo

		enddo

      

		do im = 1,nm

         if (type_mat(im).ne.2) then

            do ie = 1,ne

               if (cs(cs(ie -1) +0).eq.tag_mat(im)) then

                  do is = 1,(cs(ie) - cs(ie -1) -1)

                     in = cs(cs(ie -1) +is)

                     node_ind_ebe(in) = 0

                  enddo

               endif

            enddo

         endif

		enddo

      

		do in = 1,nnode

         if (node_ind_ebe(in).ne.0) then

            if (node_index(in).ne.0) then

               nnode_ebe = nnode_ebe +1

               node_ind_ebe(in) = nnode_ebe

            else

               node_ind_ebe(in) = 0

            endif

         endif

		enddo

      

		if (nnode_ebe.gt.0) then

         allocate(inode_ebe(nnode_ebe))

         

         do in = 1,nnode

            if (node_ind_ebe(in).ne.0) then

               inode_ebe(node_ind_ebe(in)) = in

            endif

         enddo

		endif





     ! *************************************************



	 if (nelem_N.gt.0) then

	        

			allocate(edge_nx_neuN(nedge))

			allocate(edge_ny_neuN(nedge))



			do i = 1,nedge

				edge_nx_neuN(i)=0.0d0

				edge_ny_neuN(i)=0.0d0

			enddo



	        ne = cs(0) -1



            do im = 1,nm



               !rho = prop_mat(im,1)

               !lambda = prop_mat(im,2)

               !mu = prop_mat(im,3)

               

               if (type_mat(im).eq.2) then

                  do k = 1,nelem_N

                     ie = ielem_N(k)

                     iedge = iedge_N(k)

                     

                     if (cs(cs(ie -1) +0).eq.tag_mat(im)) then

                        

                        ied1 = cs_bc(cs_bc(iedge -1) +1)

                        ied2 = cs_bc(cs_bc(iedge) -1)

                        

                        iel1 = cs(cs(ie -1) +1)

                        iel2 = cs(cs(ie -1) +nn)

                        iel3 = cs(cs(ie -1) +nn*nn)

                        iel4 = cs(cs(ie -1) +nn*(nn -1) +1)

                        

                        

! First edge

                        if (((ied1.eq.iel1).and.(ied2.eq.iel2))&

                             .or.((ied2.eq.iel1).and.(ied1.eq.iel2))) then

                           edge_lx = xs(iel2) - xs(iel1)

                           edge_ly = ys(iel2) - ys(iel1)

                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)

                           edge_nx = edge_ly / edge_ll

                           edge_ny = -1.0d0 * edge_lx / edge_ll

						   edge_nx_neuN(iedge) = edge_nx

						   edge_ny_neuN(iedge) = edge_ny

                           !edge_ia = 1; edge_ib = nn

                           !edge_ja = 1; edge_jb = 1

                        endif

                        

! Second edge

                        if (((ied1.eq.iel2).and.(ied2.eq.iel3))&

                             .or.((ied2.eq.iel2).and.(ied1.eq.iel3))) then

                           edge_lx = xs(iel3) - xs(iel2)

                           edge_ly = ys(iel3) - ys(iel2)

                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)

                           edge_nx = edge_ly / edge_ll

                           edge_ny = -1.0d0 * edge_lx / edge_ll

						   edge_nx_neuN(iedge) = edge_nx

						   edge_ny_neuN(iedge) = edge_ny

                           !edge_ia = nn; edge_ib = nn

                           !edge_ja = 1;  edge_jb = nn

                        endif

                        

! Third edge

                        if (((ied1.eq.iel3).and.(ied2.eq.iel4))&

                             .or.((ied2.eq.iel3).and.(ied1.eq.iel4))) then

                           edge_lx = xs(iel4) - xs(iel3)

                           edge_ly = ys(iel4) - ys(iel3)

                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)

                           edge_nx = edge_ly / edge_ll; 

                           edge_ny = -1.0d0 * edge_lx / edge_ll

						   edge_nx_neuN(iedge) = edge_nx

						   edge_ny_neuN(iedge) = edge_ny

                           !edge_ia = 1;  edge_ib = nn

                           !edge_ja = nn; edge_jb = nn

                        endif

                        

! Fourth edge

                        if (((ied1.eq.iel4).and.(ied2.eq.iel1))&

                             .or.((ied2.eq.iel4).and.(ied1.eq.iel1))) then

                           edge_lx = xs(iel1) - xs(iel4)

                           edge_ly = ys(iel1) - ys(iel4)

                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)

                           edge_nx = edge_ly / edge_ll

                           edge_ny = -1.0d0 * edge_lx / edge_ll

						   edge_nx_neuN(iedge) = edge_nx

						   edge_ny_neuN(iedge) = edge_ny

						   !iedge

                           !edge_ia = 1; edge_ib = 1

                           !edge_ja = 1; edge_jb = nn

                        endif

                        

                        

                        !do i = 1,nn

                        !   dxdy_el(i) = beta1(ie) + gamma1(ie) * ct(i)

                        !   dydy_el(i) = beta2(ie) + gamma2(ie) * ct(i)

                        !enddo

                        

                        !do j = 1,nn

                        !   dxdx_el(j) = alfa1(ie) + gamma1(ie) * ct(j)

                        !   dydx_el(j) = alfa2(ie) + gamma2(ie) * ct(j)

                        !enddo

                        

                        !do j = 1,nn

                        !   do i = 1,nn

                        !      is = nn*(j -1) +i

                        !      in = cs(cs(ie -1) + is)

                        !      

                        !      iaz = global_index_el_az(in)

                        !      ux_el(i,j) = u1_az(iaz)

                        !      vx_el(i,j) = (3.0d0*u1_az(iaz) -4.0d0*u0_az(iaz)&

                        !                   +1.0d0*u_1_az(iaz)) / (2.0d0*dt) 

                        !      

                        !      iaz = global_index_el_az(in +nnt)

                        !      uy_el(i,j) = u1_az(iaz)

                        !      vy_el(i,j) = (3.0d0*u1_az(iaz) -4.0d0*u0_az(iaz)&

                        !                   +1.0d0*u_1_az(iaz)) / (2.0d0*dt) 

                        !   enddo

                        !enddo

                        

                        !call make_neumannN_load(lambda,mu,rho,nn,ct,ww,dd,&

                        !          dxdx_el,dxdy_el,dydx_el,dydy_el,&

                        !          edge_nx,edge_ny,edge_ll,&

                        !          edge_ia,edge_ja,edge_ib,edge_jb,&

                        !          ux_el,uy_el,vx_el,vy_el,fx_el,fy_el)





                        !do ineuN = 1,nl_neuN

						!fn = 0

						!do ifn = 1,nfunc

						!	if (fun_neuN(ineuN).eq.tag_func(ifn)) fn = ifn

						!enddo

                        !

						!if (fn.gt.0) then

                        !

						!	do j = 1,nn

						!		do i = 1,nn

						!			is = nn*(j -1) +i

						!			in = cs(cs(ie -1) + is)

                        !     

						!			if (node_ind_ebe(in).ne.0) then

						!				id1 = node_index(in)

						!				id2 = node_index(in) +nnode_dom

                        !

						!				v1 = val_neuN(ineuN,1)

						!				v2 = val_neuN(ineuN,2)

                        !

						!				!do i = 1,nn

						!				!	in = cs_bc(cs_bc(ie -1) +i)

                        !                

						!				!	if (node_index(in).ne.0) then

						!				!		id1 = node_index(in)

						!				!		id2 = node_index(in) +nnode_dom

                        !

						!						v = 0.5d0*((ct(i) + 1.0d0)*v2 &

						!							- (ct(i) - 1.0d0)*v1)

						!						term = 0.5d0 * edge_ll * ww(i) * v

                        !      

						!				!		fmat(fn,id1) = fmat(fn,id1) + term

						!				!	endif

						!				!enddo

                        !        

						!				!iaz = update_index_el_az(id1 -1)

						!				!fk_az(iaz) = fk_az(iaz) +fx_el(i,j)

                        !         

						!				!iaz = update_index_el_az(id2 -1)

						!				!fk_az(iaz) = fk_az(iaz) +fy_el(i,j)

                        !

						!				!fmat(fn,id1) = fmat(fn,id1) + term * edge_nx !+fx_el(i,j)

						!				!fmat(fn,id2) = fmat(fn,id2) + term * edge_ny !+fy_el(i,j)

						!			endif

						!		enddo

						!	enddo

						!  endif ! if (fn.gt.0) then

						!  enddo !do ineuN = 1,nl_neuN



						

						endif

                  enddo

               endif

               

            enddo

         endif







	  endif ! if (nl_neuN.gt.0) then ! M Mexico Arroyo



	  ! *************************************************

	  ! *************************************************





	 ! *************************************************

	 ! *************************************************



     ! M Mexico Arroyo



	 if (nl_neuT.gt.0) then



		do im = 1,nm

			if ((sdeg_mat(im) +1).ne.nn) then

				deallocate(ct,ww,dd)

				deallocate(dxdx,dydx,dxdy,dydy,det_j)

            

				nn = sdeg_mat(im) +1

				allocate(ct(nn),ww(nn),dd(nn,nn))

				call lgl(nn,ct,ww,dd)

				allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))

			endif

		enddo



		!if (cs_nnz_bc.gt.0) then

		!	nedge = cs_bc(0) -1

        ! 

		!	do ie = 1,nedge

		!		if ((cs_bc(ie) - cs_bc(ie -1) -1).ne.nn) then

		!			deallocate(ct,ww,dd)

		!			nn = cs_bc(ie) - cs_bc(ie -1) -1

		!			allocate(ct(nn),ww(nn),dd(nn,nn))

		!			call lgl(nn,ct,ww,dd)

		!		endif

		!	enddo

        !

		!endif



		allocate(i4count(ne))

      

		nelem_ebe = 0 

		do ie = 1,ne

			i4count(ie) = 0

		enddo

      

		do im = 1,nm

			if (type_mat(im).eq.2) then

				do ie = 1,ne

					if (cs(cs(ie -1) +0).eq.tag_mat(im)) then

						do is = 1,(cs(ie) - cs(ie -1) -1)

							in = cs(cs(ie -1) +is)

							i4count(ie) = i4count(ie) + node_index(in)

						enddo

					endif

				enddo

			endif

		enddo

      

		do ie = 1,ne

			if (i4count(ie).gt.0) then

				nelem_ebe = nelem_ebe +1

				i4count(ie) = nelem_ebe

			endif

		enddo

      

		if (nelem_ebe.gt.0) then

			allocate(ielem_ebe(nelem_ebe))

         

			do ie = 1,ne

				if (i4count(ie).ne.0) then

				ielem_ebe(i4count(ie)) = ie

				endif

			enddo

		endif

      

		deallocate(i4count)









		nedge_T = 0 

      

		if (cs_nnz_bc.gt.0) then

			nedge = cs_bc(0) -1

         

			allocate(i4count(nedge))

				do iedge = 1,nedge

					i4count(iedge) = 0

				enddo

         

			if (nl_neuT.gt.0) then

				do i = 1,nl_neuT

					do iedge = 1,nedge

						if (cs_bc(cs_bc(iedge -1) +0).eq.tag_neuT(i)) then

							ied1 = cs_bc(cs_bc(iedge -1) +1)

							ied2 = cs_bc(cs_bc(iedge) -1)

                     

							do j = 1,nelem_ebe

								ie = ielem_ebe(j)

                        

								nn = cs_bc(iedge) - cs_bc(iedge -1) -1

                        

								iel1 = cs(cs(ie -1) +1)

								iel2 = cs(cs(ie -1) +nn)

								iel3 = cs(cs(ie -1) +nn*nn)

								iel4 = cs(cs(ie -1) +nn*(nn -1) +1)

                        

								if (((ied1.eq.iel1).or.(ied1.eq.iel2).or. &

									(ied1.eq.iel3).or.(ied1.eq.iel4)).and. &

									((ied2.eq.iel1).or.(ied2.eq.iel2).or. &

									(ied2.eq.iel3).or.(ied2.eq.iel4))) then

									i4count(iedge) = iedge

								endif

							enddo

						endif

					enddo

				enddo

			endif

         

			do iedge = 1,nedge

				if (i4count(iedge).gt.0) then

					nedge_T = nedge_T +1

					i4count(iedge) = nedge_T

				endif

			enddo

      

			if (nedge_T.gt.0) then

				allocate(iedge_T(nedge_T))

            

				do iedge = 1,nedge

					if (i4count(iedge).ne.0) then

						iedge_T(i4count(iedge)) = iedge

					endif

				enddo

			endif

         

			deallocate(i4count)

		endif

      

      

		nelem_T = nedge_T

      

		if (nelem_T.gt.0) then

			allocate(ielem_T(nelem_T))

         

			do i = 1,nedge_T

				iedge = iedge_T(i)

            

				ied1 = cs_bc(cs_bc(iedge -1) +1)

				ied2 = cs_bc(cs_bc(iedge) -1)

            

				do j = 1,nelem_ebe

					ie = ielem_ebe(j)

               

					nn = cs_bc(iedge) - cs_bc(iedge -1) -1

               

					iel1 = cs(cs(ie -1) +1)

					iel2 = cs(cs(ie -1) +nn)

					iel3 = cs(cs(ie -1) +nn*nn)

					iel4 = cs(cs(ie -1) +nn*(nn -1) +1)

               

					if (((ied1.eq.iel1).or.(ied1.eq.iel2).or. &

						(ied1.eq.iel3).or.(ied1.eq.iel4)).and. &

						((ied2.eq.iel1).or.(ied2.eq.iel2).or. &

						(ied2.eq.iel3).or.(ied2.eq.iel4))) then

						ielem_T(i) = ie

					endif

				enddo

			enddo

		endif ! M Mexico Arroyo





		allocate(node_ind_ebe(nnode))

      

		nnode_ebe = 0

		do in = 1,nnode

			node_ind_ebe(in) = 0

		enddo

      

		do i = 1,nelem_ebe

			ie = ielem_ebe(i)

			do is = 1,(cs(ie) - cs(ie -1) -1)

				in = cs(cs(ie -1) +is)

				node_ind_ebe(in) = 1

			enddo

		enddo

      

		do im = 1,nm

         if (type_mat(im).ne.2) then

            do ie = 1,ne

               if (cs(cs(ie -1) +0).eq.tag_mat(im)) then

                  do is = 1,(cs(ie) - cs(ie -1) -1)

                     in = cs(cs(ie -1) +is)

                     node_ind_ebe(in) = 0

                  enddo

               endif

            enddo

         endif

		enddo

      

		do in = 1,nnode

         if (node_ind_ebe(in).ne.0) then

            if (node_index(in).ne.0) then

               nnode_ebe = nnode_ebe +1

               node_ind_ebe(in) = nnode_ebe

            else

               node_ind_ebe(in) = 0

            endif

         endif

		enddo

      

		if (nnode_ebe.gt.0) then

         allocate(inode_ebe(nnode_ebe))

         

         do in = 1,nnode

            if (node_ind_ebe(in).ne.0) then

               inode_ebe(node_ind_ebe(in)) = in

            endif

         enddo

		endif





     ! *************************************************



	 if (nelem_T.gt.0) then



		allocate(edge_nx_neuT(nedge))

		allocate(edge_ny_neuT(nedge))



		do i = 1,nelem_T

			edge_nx_neuT(i)=0.0d0

			edge_ny_neuT(i)=0.0d0

		enddo



		ne = cs(0) -1



	 

            do im = 1,nm

               !rho = prop_mat(im,1)

               !lambda = prop_mat(im,2)

               !mu = prop_mat(im,3)

               

               if (type_mat(im).eq.2) then

                  do k = 1,nelem_T

                     ie = ielem_T(k)

                     iedge = iedge_T(k)

                     

                     if (cs(cs(ie -1) +0).eq.tag_mat(im)) then

                        

                        ied1 = cs_bc(cs_bc(iedge -1) +1)

                        ied2 = cs_bc(cs_bc(iedge) -1)

                        

                        iel1 = cs(cs(ie -1) +1)

                        iel2 = cs(cs(ie -1) +nn)

                        iel3 = cs(cs(ie -1) +nn*nn)

                        iel4 = cs(cs(ie -1) +nn*(nn -1) +1)

                        

                        

! First edge

                        if (((ied1.eq.iel1).and.(ied2.eq.iel2))&

                             .or.((ied2.eq.iel1).and.(ied1.eq.iel2))) then

                           edge_lx = xs(iel2) - xs(iel1)

                           edge_ly = ys(iel2) - ys(iel1)

                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)

                           edge_nx = edge_ly / edge_ll

                           edge_ny = -1.0d0 * edge_lx / edge_ll

						   edge_nx_neuT(iedge) = edge_nx

						   edge_ny_neuT(iedge) = edge_ny

                           !edge_ia = 1; edge_ib = nn

                           !edge_ja = 1; edge_jb = 1

                        endif

                        

! Second edge

                        if (((ied1.eq.iel2).and.(ied2.eq.iel3))&

                             .or.((ied2.eq.iel2).and.(ied1.eq.iel3))) then

                           edge_lx = xs(iel3) - xs(iel2)

                           edge_ly = ys(iel3) - ys(iel2)

                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)

                           edge_nx = edge_ly / edge_ll

                           edge_ny = -1.0d0 * edge_lx / edge_ll

						   edge_nx_neuT(iedge) = edge_nx

						   edge_ny_neuT(iedge) = edge_ny

                           !edge_ia = nn; edge_ib = nn

                           !edge_ja = 1;  edge_jb = nn

                        endif

                        

! Third edge

                        if (((ied1.eq.iel3).and.(ied2.eq.iel4))&

                             .or.((ied2.eq.iel3).and.(ied1.eq.iel4))) then

                           edge_lx = xs(iel4) - xs(iel3)

                           edge_ly = ys(iel4) - ys(iel3)

                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)

                           edge_nx = edge_ly / edge_ll; 

                           edge_ny = -1.0d0 * edge_lx / edge_ll

						   edge_nx_neuT(iedge) = edge_nx

						   edge_ny_neuT(iedge) = edge_ny

                           !edge_ia = 1;  edge_ib = nn

                           !edge_ja = nn; edge_jb = nn

                        endif

                        

! Fourth edge

                        if (((ied1.eq.iel4).and.(ied2.eq.iel1))&

                             .or.((ied2.eq.iel4).and.(ied1.eq.iel1))) then

                           edge_lx = xs(iel1) - xs(iel4)

                           edge_ly = ys(iel1) - ys(iel4)

                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)

                           edge_nx = edge_ly / edge_ll

                           edge_ny = -1.0d0 * edge_lx / edge_ll

						   edge_nx_neuT(iedge) = edge_nx

						   edge_ny_neuT(iedge) = edge_ny

                           !edge_ia = 1; edge_ib = 1

                           !edge_ja = 1; edge_jb = nn

                        endif

                        

                        

                        !do i = 1,nn

                        !   dxdy_el(i) = beta1(ie) + gamma1(ie) * ct(i)

                        !   dydy_el(i) = beta2(ie) + gamma2(ie) * ct(i)

                        !enddo

                        

                        !do j = 1,nn

                        !   dxdx_el(j) = alfa1(ie) + gamma1(ie) * ct(j)

                        !   dydx_el(j) = alfa2(ie) + gamma2(ie) * ct(j)

                        !enddo

                        

                        !do j = 1,nn

                        !   do i = 1,nn

                        !      is = nn*(j -1) +i

                        !      in = cs(cs(ie -1) + is)

                        !      

                        !      iaz = global_index_el_az(in)

                        !      ux_el(i,j) = u1_az(iaz)

                        !      vx_el(i,j) = (3.0d0*u1_az(iaz) -4.0d0*u0_az(iaz)&

                        !                   +1.0d0*u_1_az(iaz)) / (2.0d0*dt) 

                        !      

                        !      iaz = global_index_el_az(in +nnt)

                        !      uy_el(i,j) = u1_az(iaz)

                        !      vy_el(i,j) = (3.0d0*u1_az(iaz) -4.0d0*u0_az(iaz)&

                        !                   +1.0d0*u_1_az(iaz)) / (2.0d0*dt) 

                        !   enddo

                        !enddo

                        

                        !call make_neumannN_load(lambda,mu,rho,nn,ct,ww,dd,&

                        !          dxdx_el,dxdy_el,dydx_el,dydy_el,&

                        !          edge_nx,edge_ny,edge_ll,&

                        !          edge_ia,edge_ja,edge_ib,edge_jb,&

                        !          ux_el,uy_el,vx_el,vy_el,fx_el,fy_el)





                        !do ineuT = 1,nl_neuT

						!fn = 0

						!do ifn = 1,nfunc

						!	if (fun_neuT(ineuT).eq.tag_func(ifn)) fn = ifn

						!enddo

                     

						!if (fn.gt.0) then

                        !

						!	!do j = 1,nn

						!	!	do i = 1,nn

						!	!		is = nn*(j -1) +i

						!	!		in = cs(cs(ie -1) + is)

                        !   !  

						!	!		if (node_ind_ebe(in).ne.0) then

						!	!			id1 = node_index(in)

						!	!			id2 = node_index(in) +nnode_dom

                        !

						!				v1 = val_neuT(ineuT,1)

						!				v2 = val_neuT(ineuT,2)

                        !

						!				do i = 1,nn

						!					in = cs_bc(cs_bc(ie -1) +i)

                        !   

						!					if (node_index(in).ne.0) then

						!						id1 = node_index(in)

                        !

						!						v = 0.5d0*((ct(i) + 1.0d0)*v2 &

						!							- (ct(i) - 1.0d0)*v1)

						!						term = 0.5d0 * edge_ll * ww(i) * v

                        !      

						!						!fmat(fn,id1) = fmat(fn,id1) + term

						!					!endif

						!				!enddo

                        !         

						!				!iaz = update_index_el_az(id1 -1)

						!				!fk_az(iaz) = fk_az(iaz) +fx_el(i,j)

                        !         

						!				!iaz = update_index_el_az(id2 -1)

						!				!fk_az(iaz) = fk_az(iaz) +fy_el(i,j)

                        !

						!				fmat(fn,id1) = fmat(fn,id1) + term * edge_nx !+fx_el(i,j)

						!				fmat(fn,id2) = fmat(fn,id2) + term * edge_ny !+fy_el(i,j)

						!			endif

						!		enddo

						!	!enddo

						!  endif ! if (fn.gt.0) then

						!  enddo !do ineuT = 1,nl_neuT



						

						endif

                  enddo

               endif

               

            enddo

         endif







	  endif ! if (nl_neuT.gt.0) then ! M Mexico Arroyo



	  ! *************************************************

	  ! *************************************************






      
      do i = 1,nl_poiX
         call find_nearest_node(nnode,xs,ys,val_poiX(i,1),val_poiX(i,2),&
                                node_poiX(i))
      enddo
      do i = 1,nl_poiY
         call find_nearest_node(nnode,xs,ys,val_poiY(i,1),val_poiY(i,2),&
                                node_poiY(i))
      enddo



	  do i = 1,nl_dipX ! Marco

         call find_nearest_node(nnode,xs,ys,val_dipX(i,1),val_dipX(i,2),& ! Marco

                                node_dipX(i)) ! Marco

		 !inode_dipx(i) = node_dipX(i)

      enddo ! Marco

      do i = 1,nl_dipY ! Marco

         call find_nearest_node(nnode,xs,ys,val_dipY(i,1),val_dipY(i,2),& ! Marco

                                node_dipY(i)) ! Marco

      enddo ! Marco

      ne = cs(0) -1
      
      do im = 1,nm
         if ((sdeg_mat(im) +1).ne.nn) then
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
         
               
               
! Point load X
                  
               if (nl_poiX.gt.0) then


                  do ip = 1,nl_poiX
                     fn = 0
                     do ifn = 1,nfunc
                        if (fun_poiX(ip).eq.tag_func(ifn)) fn = ifn
                     enddo
                     
                     if (fn.gt.0) then
                        do j = 1,nn
                           do i = 1,nn
                              is = nn*(j -1) +i
                              in = cs(cs(ie -1) +is)
                              
                              if (node_index(in).ne.0) then
                                 id1 = node_index(in)
                                 
                                 if (in.eq.node_poiX(ip)) then
                                    term = val_poiX(ip,3) &
                                         * (det_j(i,j) * ww(i) * ww(j))
                                    fmat(fn,id1) = fmat(fn,id1) + term
                                 endif
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
                     do ifn = 1,nfunc
                        if (fun_poiY(ip).eq.tag_func(ifn)) fn = ifn
                     enddo
                     
                     if (fn.gt.0) then
                        do j = 1,nn
                           do i = 1,nn
                              is = nn*(j -1) +i
                              in = cs(cs(ie -1) +is)
                              
                              if (node_index(in).ne.0) then
                                 id2 = node_index(in) +nnode_dom
                                 
                                 if (in.eq.node_poiY(ip)) then
                                    term = val_poiY(ip,3) &
                                         * (det_j(i,j) * ww(i) * ww(j))
                                    fmat(fn,id2) = fmat(fn,id2) + term
                                 endif
                              endif
                           enddo
                        enddo
                     endif
                  enddo
               endif



! Point load All X

                  

               if (nl_palX.gt.0) then



			      

                  do ip = 1,nl_palX



				     if (tag_mat(im).eq.tag_palX(ip)) then



                     fn = 0

                     do ifn = 1,nfunc

                        if (fun_palX(ip).eq.tag_func(ifn)) fn = ifn

                     enddo

                     

                     if (fn.gt.0) then

                        do j = 1,nn

                           do i = 1,nn

                              is = nn*(j -1) +i

                              in = cs(cs(ie -1) +is)

                              

                              if (node_index(in).ne.0) then

                                 id1 = node_index(in)

                                 

                                 !if (in.eq.node_poiX(ip)) then

                                    term = rho * val_palX(ip,1) &

                                         * (det_j(i,j) * ww(i) * ww(j))

                                    fmat(fn,id1) = fmat(fn,id1) + term

                                 !endif

                              endif

                           enddo

                        enddo

                     endif

					 endif

                  enddo

               endif



! Point load All Y

               

               if (nl_palY.gt.0) then



                  do ip = 1,nl_palY



					 if (tag_mat(im).eq.tag_palY(ip)) then



                     fn = 0

                     do ifn = 1,nfunc

                        if (fun_palY(ip).eq.tag_func(ifn)) fn = ifn

                     enddo

                     

                     if (fn.gt.0) then

                        do j = 1,nn

                           do i = 1,nn

                              is = nn*(j -1) +i

                              in = cs(cs(ie -1) +is)

                              

                              if (node_index(in).ne.0) then

                                 id2 = node_index(in) +nnode_dom

                                 

                                 !if (in.eq.node_palY(ip)) then

                                    term = rho * val_palY(ip,1) &

                                         * (det_j(i,j) * ww(i) * ww(j))

                                    fmat(fn,id2) = fmat(fn,id2) + term

                                 !endif

                              endif

                           enddo

                        enddo

                     endif

					 endif

                  enddo

               endif



! Dirichlet Point X

                  

               if (nl_dipX.gt.0) then

                  do idp = 1,nl_dipX

                     fn = 0

                     do ifn = 1,nfunc

                        if (fun_dipX(idp).eq.tag_func(ifn)) fn = ifn

                     enddo

                     

                     if (fn.gt.0) then

                        do j = 1,nn

                           do i = 1,nn

                              is = nn*(j -1) +i

                              in = cs(cs(ie -1) +is)

                              

                              if (node_index(in).ne.0) then

                                 id1 = node_index(in)

                                 

                                 if (in.eq.node_dipX(idp)) then

								    term = val_dipX(idp,3)

                                    !term = val_dipX(idp,3) &

                                    !     * (det_j(i,j) * ww(i) * ww(j))

                                    !fmat(fn,id1) = fmat(fn,id1) + term

									fmat(fn,id1) = term

                                 endif

                              endif

                           enddo

                        enddo

                     endif

                  enddo

               endif

              

               

! Dirichlet Point Y

               

               if (nl_dipY.gt.0) then

                  do idp = 1,nl_dipY

                     fn = 0

                     do ifn = 1,nfunc

                        if (fun_dipY(idp).eq.tag_func(ifn)) fn = ifn

                     enddo

                     

                     if (fn.gt.0) then

                        do j = 1,nn

                           do i = 1,nn

                              is = nn*(j -1) +i

                              in = cs(cs(ie -1) +is)

                              

                              if (node_index(in).ne.0) then

                                 id2 = node_index(in) +nnode_dom

                                 

                                 if (in.eq.node_dipY(idp)) then

									term = val_dipY(idp,3)

                                    !term = val_dipY(idp,3) &

                                    !     * (det_j(i,j) * ww(i) * ww(j))

                                    !fmat(fn,id2) = fmat(fn,id2) + term

									fmat(fn,id2) = term

                                 endif

                              endif

                           enddo

                        enddo

                     endif

                  enddo

               endif

! Plane Wave X
               if (nl_plaX.gt.0) then


                 ! Find out plane wave condition load

                  do ipl = 1,nl_plaX

                     !Check on the Plane Wave Material - start

                     if (tag_mat(im).eq.tag_plaX(ipl)) then
                        C=dsqrt(mu*rho)
                        fn = 0
                        do ifn = 1,nfunc
                           if (fun_plaX(ipl).eq.tag_func(ifn)) fn = ifn
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
                              
                                 if (node_index(in).ne.0) then
                                    id1 = node_index(in)
                                 
                                    term = C*val_plaX(ipl,1) &
                                            * ellex*ww(i)
                                    fmat(fn,id1) = fmat(fn,id1) + term
                                 endif
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
                              
                                 if (node_index(in).ne.0) then
                                    id1 = node_index(in)
                                       term = C*val_plaX(ipl,1) &
                                            *ellex*ww(j)
                                       fmat(fn,id1) = fmat(fn,id1) + term
                                  
                                 endif
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
                        do ifn = 1,nfunc
                           if (fun_plaY(ipl).eq.tag_func(ifn)) fn = ifn
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
                              
                                 if (node_index(in).ne.0) then
                                    id2 = node_index(in) + nnode_dom
                                 
                                    term = C*val_plaY(ipl,1) &
                                            * ellex*ww(i)
                                    fmat(fn,id2) = fmat(fn,id2) + term
                                 endif
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
                              
                                 if (node_index(in).ne.0) then
                                    id2 = node_index(in) + nnode_dom
                                       term = C*val_plaY(ipl,1) &
                                            *ellex*ww(j)
                                       fmat(fn,id2) = fmat(fn,id2) + term
                                  
                                 endif
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
                     !write(49,*)'--------------------------'
                     !write(49,*)'isism =',isism
                     !write(49,*)'--------------------------'
                     !write(49,*),slip1,slip2,norm1,norm2
          	     do k = 1,num_node_sism(isism)        

			do j = 1,nn
                           do i = 1,nn
                              is = nn*(j -1) +i
                              in = cs(cs(ie -1) +is)
				
                                 if (in.eq.sour_node_sism(k,isism)) then
				    !write(49,*)'TROVATO isism=',isism,'/',nl_sism
				    facsmom(isism,1) = facsmom(isism,1) + det_j(i,j) * ww(i) * ww(j)
                                    facsmom(isism,2) = facsmom(isism,2) + det_j(i,j) * ww(i) * ww(j)
                                    facsmom(isism,3) = facsmom(isism,3) + det_j(i,j) * ww(i) * ww(j)
				    !write(49,*)'facsmom=',facsmom(isism,1)
                                    !write(49,*)'facsmom=',facsmom(isism,2)
                                    !write(49,*)'facsmom=',facsmom(isism,3)

				    !write(49,*)'in     =',in,'/',sour_node_sism(k,isism)
                                    length_cns = length_cns + 1
                                 endif
                                 
			    enddo
		        enddo

                      enddo
!		    write(49,*)'isism_facsmom =',facsmom(isism,1),' |isism =',isism
!		    write(49,*)'IE =',ie,' | nodo:',in

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
            
            lx = xs(cs_bc(cs_bc(ie -1) +nn)) &
                 - xs(cs_bc(cs_bc(ie -1) +1))
            ly = ys(cs_bc(cs_bc(ie -1) +nn)) &
                 - ys(cs_bc(cs_bc(ie -1) +1))
            ll = dsqrt(lx*lx + ly*ly)
            


            
! Neumann load X
            
            if (nl_neuX.gt.0) then
               do il = 1,nl_neuX
                  if (tag_neuX(il).eq.cs_bc(cs_bc(ie -1) +0)) then
                     fn = 0
                     do ifn = 1,nfunc
                        if (fun_neuX(il).eq.tag_func(ifn)) fn = ifn
                     enddo
                     
                     if (fn.gt.0) then
                        v1 = val_neuX(il,1)
                        v2 = val_neuX(il,2)
                        
                        do i = 1,nn
                           in = cs_bc(cs_bc(ie -1) +i)
                           
                           if (node_index(in).ne.0) then
                              id1 = node_index(in)
                     
                              v = 0.5d0*((ct(i) + 1.0d0)*v2 &
                                   - (ct(i) - 1.0d0)*v1)
                              term = 0.5d0 * ll * ww(i) * v
                              
                              fmat(fn,id1) = fmat(fn,id1) + term
                           endif
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
                     do ifn = 1,nfunc
                        if (fun_neuY(il).eq.tag_func(ifn)) fn = ifn
                     enddo
                     
                     if (fn.gt.0) then
                        v1 = val_neuY(il,1)
                        v2 = val_neuY(il,2)
                        
                        do i = 1,nn
                           in = cs_bc(cs_bc(ie -1) +i)
                           
                           if (node_index(in).ne.0) then
                              id2 = node_index(in) +nnode_dom
                              
                              v = 0.5d0*((ct(i) + 1.0d0)*v2 &
                                   - (ct(i) - 1.0d0)*v1)
                              term = 0.5d0 * ll * ww(i) * v
                              
                              fmat(fn,id2) = fmat(fn,id2) + term
                           endif
                        enddo
                     endif
                  endif
               enddo
            endif





! Neumann load N ! M Mexico Arroyo



            

            if (nl_neuN.gt.0) then ! M Mexico Arroyo

               do il = 1,nl_neuN ! M Mexico Arroyo

                  if (tag_neuN(il).eq.cs_bc(cs_bc(ie -1) +0)) then ! M Mexico Arroyo

                     fn = 0 ! M Mexico Arroyo

                     do ifn = 1,nfunc ! M Mexico Arroyo

                        if (fun_neuN(il).eq.tag_func(ifn)) fn = ifn ! M Mexico Arroyo

                     enddo ! M Mexico Arroyo

                     

                     if (fn.gt.0) then ! M Mexico Arroyo

                        v1 = val_neuN(il,1) ! M Mexico Arroyo

                        v2 = val_neuN(il,2) ! M Mexico Arroyo

                        

                        do i = 1,nn ! M Mexico Arroyo

                           in = cs_bc(cs_bc(ie -1) +i) ! M Mexico Arroyo

                           

                           if (node_index(in).ne.0) then ! M Mexico Arroyo



                              id1 = node_index(in) ! M Mexico Arroyo

                              v = 0.5d0*((ct(i) + 1.0d0)*v2 & ! M Mexico Arroyo

                                   - (ct(i) - 1.0d0)*v1) ! M Mexico Arroyo

                              term = 0.5d0 * ll * ww(i) * v ! M Mexico Arroyo

                              fmat(fn,id1) = fmat(fn,id1) + term * edge_nx_neuN(ie) ! M Mexico Arroyo



							  id2 = node_index(in) +nnode_dom ! M Mexico Arroyo

                              v = 0.5d0*((ct(i) + 1.0d0)*v2 & ! M Mexico Arroyo

                                   - (ct(i) - 1.0d0)*v1) ! M Mexico Arroyo

                              term = 0.5d0 * ll * ww(i) * v ! M Mexico Arroyo

                              fmat(fn,id2) = fmat(fn,id2) + term * edge_ny_neuN(ie) ! M Mexico Arroyo



                           endif ! M Mexico Arroyo

                        enddo ! M Mexico Arroyo

                     endif ! M Mexico Arroyo

                  endif ! M Mexico Arroyo

               enddo ! M Mexico Arroyo

            endif ! M Mexico Arroyo

            

            

! Neumann load T ! M Mexico Arroyo

            

            if (nl_neuT.gt.0) then ! M Mexico Arroyo

               do il = 1,nl_neuT ! M Mexico Arroyo

                  if (tag_neuT(il).eq.cs_bc(cs_bc(ie -1) +0)) then ! M Mexico Arroyo

                     fn = 0 ! M Mexico Arroyo

                     do ifn = 1,nfunc ! M Mexico Arroyo

                        if (fun_neuT(il).eq.tag_func(ifn)) fn = ifn ! M Mexico Arroyo

                     enddo ! M Mexico Arroyo

                     

                     if (fn.gt.0) then ! M Mexico Arroyo

                        v1 = val_neuT(il,1) ! M Mexico Arroyo

                        v2 = val_neuT(il,2) ! M Mexico Arroyo

                        

                        do i = 1,nn ! M Mexico Arroyo

                           in = cs_bc(cs_bc(ie -1) +i) ! M Mexico Arroyo

                           

                           if (node_index(in).ne.0) then ! M Mexico Arroyo



                              id1 = node_index(in) ! M Mexico Arroyo

                              v = 0.5d0*((ct(i) + 1.0d0)*v2 & ! M Mexico Arroyo

                                   - (ct(i) - 1.0d0)*v1) ! M Mexico Arroyo

                              term = 0.5d0 * ll * ww(i) * v ! M Mexico Arroyo

                              fmat(fn,id1) = fmat(fn,id1) + term * edge_ny_neuT(ie)  ! M Mexico Arroyo!



							  id2 = node_index(in) +nnode_dom ! M Mexico Arroyo

                              v = 0.5d0*((ct(i) + 1.0d0)*v2 & ! M Mexico Arroyo

                                   - (ct(i) - 1.0d0)*v1) ! M Mexico Arroyo

                              term = 0.5d0 * ll * ww(i) * v ! M Mexico Arroyo

                              fmat(fn,id2) = fmat(fn,id2) + term * edge_nx_neuT(ie)  ! M Mexico Arroyo



                           endif ! M Mexico Arroyo

                        enddo ! M Mexico Arroyo

                     endif ! M Mexico Arroyo

                  endif ! M Mexico Arroyo

               enddo ! M Mexico Arroyo

            endif ! M Mexico Arroyo
            
            
! Dirichlet X
            
            if (nl_dirX.gt.0) then
               do il = 1,nl_dirX
                  if (tag_dirX(il).eq.cs_bc(cs_bc(ie -1) +0)) then
                     fn = 0
                     do ifn = 1,nfunc
                        if (fun_dirX(il).eq.tag_func(ifn)) fn = ifn
                     enddo
                     
                     if (fn.gt.0) then
                        v1 = val_dirX(il,1)
                        v2 = val_dirX(il,2)
                        
                        do i = 1,nn
                           in = cs_bc(cs_bc(ie -1) +i)
                           
                           if (node_index(in).ne.0) then
                              id1 = node_index(in)
                              
                              do ifn = 1,nfunc
                                 fmat(ifn,id1) = 0.0d0
                              enddo
                              
                              v = 0.5d0*((ct(i) + 1.0d0)*v2 &
                                   - (ct(i) - 1.0d0)*v1)
                              
                              fmat(fn,id1) = v
                           endif
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
                     do ifn = 1,nfunc
                        if (fun_dirY(il).eq.tag_func(ifn)) fn = ifn
                     enddo
                     
                     if (fn.gt.0) then
                        v1 = val_dirY(il,1)
                        v2 = val_dirY(il,2)
                        
                        do i = 1,nn
                           in = cs_bc(cs_bc(ie -1) +i)
                           
                           if (node_index(in).ne.0) then
                              id2 = node_index(in) +nnode_dom
                              
                              do ifn = 1,nfunc
                                 fmat(ifn,id2) = 0.0d0
                              enddo
                              
                              v = 0.5d0*((ct(i) + 1.0d0)*v2 &
                                   - (ct(i) - 1.0d0)*v1)
                              
                              fmat(fn,id2) = v
                           endif
                        enddo
                     endif
                  endif
               enddo
            endif
            
         enddo
      endif
      
      !do isism = 1,nl_sism
      !	 write(49,*)'facsmom=',facsmom(isism,1),facsmom(isism,2),facsmom(isism,3)
      !enddo	

      
      deallocate(ct,ww,dd)
      deallocate(dxdx,dydx,dxdy,dydy,det_j)
      
      return
      end subroutine make_Fel
                  
      
      
!     ********************************************************
      
      subroutine time_loop_el(nnt,xs,ys,node_index,cs_nnz,cs,&
                              nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&

							  nmatg,valmatg,tagmatg,typematg,& ! Mexico Paco grad_lin 29.11.2004
                              ne,alfa1,beta1,gamma1,alfa2,beta2,gamma2,&
                              cs_nnz_bc,cs_bc,&
                              nl_dirX,tag_dirX,nl_dirY,tag_dirY,&

							  nl_dipX,in_dipX,nl_dipY,in_dipY,& ! Marco
                              nl_abc,tag_abc,&
                              nf,func_type,func_indx,func_data,nfunc_data,&
                              nnd,mvec,cvec,kcvec,Kel_nnz,Kel_bin,&
                              Fmat,u0,u_1,&
                              nts,dt,nmonit,node_m,nsnap,itersnap,&
                              file_out,myid,&

                              check_node_sism,check_dist_node_sism,&  
                              length_cns,facsmom,nl_sism,&
                              make_damping_yes_or_not,ndt_monitor,&
							  
							  nvpcl,valvpcl,fvpcl,tagvpcl,& ! Clara

							  nvpsa,valvpsa,tagvpsa,& ! MCNAP

							  nvpsl,valvpsl,tagvpsl,& ! MCNAP

							  nvpsd,valvpsd,tagvpsd,& ! MCNAP

							  nmaps,valmaps,fmaps,tagmaps,& ! Marco

							  nsnapps,itersnapps,nprocs) ! Marco
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Luca Massidda
!     
      implicit none
      
      integer*4 :: nnt,cs_nnz,nm,ne,cs_nnz_bc,nf
      integer*4 :: nl_dirX,nl_dirY,nl_abc

	  integer*4 :: nl_dipX,nl_dipY ! Marco
      integer*4 :: nnd,Kel_nnz,Kth_nnz,nts,nmonit,nsnap

	  integer*4 :: nsnapps ! Marco
      integer*4 :: myid
      real*8, dimension(nnt) :: xs,ys
      integer*4, dimension(nnt) :: node_index
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat
      real*8, dimension(nm) :: tref_mat
      real*8, dimension(nm,4) :: prop_mat
      real*8, dimension(ne) :: alfa1,beta1,gamma1,alfa2,beta2,gamma2
      integer*4, dimension(0:cs_nnz_bc) :: cs_bc
      integer*4, dimension(nl_dirX) :: tag_dirX
      integer*4, dimension(nl_dirY) :: tag_dirY

	  integer*4, dimension(nl_dipX) :: in_dipX ! Marco

      integer*4, dimension(nl_dipY) :: in_dipY ! Marco
      integer*4, dimension(nl_abc) :: tag_abc
      integer*4, dimension(nf) :: func_type
      integer*4, dimension(nf +1) :: func_indx
      integer*4 :: nfunc_data

	  real*8, dimension(nfunc_data) :: func_data
      real*8, dimension(0:(2*nnd -1)) :: mvec,cvec,kcvec,new
      integer*4, dimension(0:Kel_nnz) :: Kel_bin
      real*8, dimension(nf,2*nnd) :: Fmat
      real*8, dimension(2*nnd) :: u0,u_1
      integer*4, dimension(nmonit) :: node_m
      integer*4, dimension(nsnap) :: itersnap

	  integer*4, dimension(nsnapps) :: itersnapps ! Marco
      real*8 :: dt
      character*70 :: file_out
      
      
      integer*4, dimension(3) :: clock
      integer*4 :: clock_start,clock_finish
      integer*4 :: isnap,its,fn,ie,nn,ip,im

      integer*4 :: isnapps ! Marco
      integer*4 :: i,is,in,id1,id2,idt,j,jn,k,kn,iaz,jaz,kaz
      real*4 :: time_in_seconds,time_total
      real*8 :: tt0,tt1,tt2,dt2
      real*8 :: r8t,eps
      integer*4 :: i4t
      
      real*8 :: get_func_value
      real*8, dimension(:), allocatable :: func_value
      
      character*70 :: file_timing,file_monitor,file_stress
      character*70 :: file_outE
      character*70 :: file_outV ! Clara
	  character*70 :: file_outT1 ! Clara
	  character*70 :: file_outT2 ! Clara
	  character*70 :: file_outT3 ! Clara

	  character*70 :: file_outS1 ! Clara

	  character*70 :: file_outS2 ! Clara

	  character*70 :: file_outS3 ! Clara

	  !character*70 :: file_outDvp1 ! MCNAP

	  !character*70 :: file_outDvp2 ! MCNAP

      !character*70 :: file_outDvp3 ! MCNAP

	  character*70 :: file_out_u1 ! M Mexico _A_ 

	  character*70 :: file_out_u0 ! M Mexico _A_ 

	  character*70 :: file_out_Fe ! M Mexico _A_ 

	  character*70 :: file_out_Fk ! M Mexico _A_ 

	  character*70 :: file_out_vivel ! M Mexico _A_

	  character*70 :: file_out_vistn1_0 ! M Mexico _A_ 

	  character*70 :: file_out_vistn2_0  ! M Mexico _A_

	  character*70 :: file_out_vistn3_0  ! M Mexico _A_  




      integer*4 :: unit_monitor,lname

      
!*************************************************
!
!             BOUNDARY CONDITIONS
!
!*************************************************
      
      integer*4, dimension(:), allocatable :: i4count
      
      integer*4, dimension(:), allocatable :: inode_dirX
      integer*4, dimension(:), allocatable :: inode_dirY
      integer*4, dimension(:), allocatable :: inode_neuX
      integer*4, dimension(:), allocatable :: inode_neuY

	  integer*4, dimension(:), allocatable :: inode_neuN ! M Mexico Arroyo

      integer*4, dimension(:), allocatable :: inode_neuT ! M Mexico Arroyo
      integer*4, dimension(:), allocatable :: ielem_abc,iedge_abc

      integer*4, dimension(:), allocatable :: ielem_ebe
      integer*4, dimension(:), allocatable :: inode_ebe
      
      integer*4, dimension(:), allocatable :: node_ind_ebe
      
      integer*4 :: nnode_dirX,nnode_dirY,nnode_neuX,nnode_neuY

	  integer*4 :: nnode_neuN,nnode_neuT ! M Mexico Arroyo
      integer*4 :: nelem_abc,nedge_abc
      
      integer*4 :: iedge,nedge,ied1,ied2,iel1,iel2,iel3,iel4
      integer*4 :: edge_ia,edge_ja,edge_ib,edge_jb
      real*8 :: edge_lx,edge_ly,edge_ll,edge_nx,edge_ny
      
      
      
!*************************************************
!
!              ELEMENT BY ELEMENT
!
!*************************************************
      
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8, dimension(:), allocatable :: dxdx_el,dxdy_el,dydx_el,dydy_el
      real*8, dimension(:,:), allocatable :: ux_el,uy_el
      real*8, dimension(:,:), allocatable :: vx_el,vy_el
      real*8, dimension(:,:), allocatable :: duxdx_el,duxdy_el
      real*8, dimension(:,:), allocatable :: duydx_el,duydy_el
      real*8, dimension(:,:), allocatable :: sxx_el,syy_el,sxy_el,fx_el,fy_el
      real*8 :: lambda,mu,rho
      integer*4 :: nelem_ebe,nnode_ebe
      
      
      
!*************************************************
!
!                    AZTEC
!
!*************************************************
      
      real*8, dimension(:), allocatable :: u_1_az,u0_az,u1_az,u2_az
      real*8, dimension(:), allocatable :: fk_az,fe_az
      
      integer*4, dimension(:), allocatable :: out_index_el_az
      
      integer*4, dimension(:), allocatable :: proc_config_az
      real*8, dimension(:), allocatable :: status_az,params_az
      integer*4, dimension(:), allocatable :: options_az
     
      integer*4, dimension(:), allocatable :: update_el_az,update_index_el_az
      integer*4, dimension(:), allocatable :: extern_el_az,extern_index_el_az
      integer*4, dimension(:), allocatable :: data_org_el_az
      integer*4, dimension(:), allocatable :: extern_proc_el_az
      integer*4, dimension(:), allocatable :: global_index_el_az
      integer*4 :: N_intern_el_az,N_border_el_az,N_update_el_az,N_extern_el_az



!*************************************************
!
!                  SEISMIC MOMENT
!
!*************************************************
      
      integer*4 :: length_cns
	  integer*4 :: nl_sism 
	  integer*4, dimension(length_cns,5) :: check_node_sism
      real*8, dimension(length_cns,1) :: check_dist_node_sism
      real*8, dimension(nl_sism,3) :: facsmom
      
      
!*************************************************
!
!                      DAMPING
!
!*************************************************      
      
     integer*4 :: make_damping_yes_or_not

	 real*8 :: ndt_monitor



!*************************************************

!

!              LINEAR GRADIENT MATERIAL

!

!*************************************************



      !!! SISTEMARE LE VARIABILI CHE NON SERVONO !!!

      integer*4 :: nmatg,imatg, which_matg ! Mexico Paco grad_lin 29.11.2004

      real*8, dimension(nmatg,*) :: valmatg ! Mexico Paco grad_lin 29.11.2004

      integer*4, dimension(*) :: tagmatg ! Mexico Paco grad_lin 29.11.2004

	  integer*4, dimension(*) :: typematg ! Mexico Paco grad_lin 29.11.2004

!*************************************************
!
!              VISCO-PLASTIC MATERIAL
!
!*************************************************

      !!! SISTEMARE LE VARIABILI CHE NON SERVONO !!!
      integer*4 :: nvpcl,ivpcl, which_vpcl ! Clara
      real*8, dimension(nvpcl,*) :: valvpcl ! Clara  
      integer*4, dimension(*) :: fvpcl,tagvpcl ! Clara



	  integer*4 :: nvpsa,ivpsa, which_vpsa ! MCNAP

      real*8, dimension(nvpsa,*) :: valvpsa ! MCNAP

      integer*4, dimension(*) :: tagvpsa ! MCNAP



	  integer*4 :: nvpsl,ivpsl, which_vpsl ! MCNAP

      real*8, dimension(nvpsl,*) :: valvpsl ! MCNAP

      integer*4, dimension(*) :: tagvpsl ! MCNAP



	  integer*4 :: nvpsd,ivpsd, which_vpsd ! MCNAP

      real*8, dimension(nvpsd,*) :: valvpsd ! MCNAP

      integer*4, dimension(*) :: tagvpsd ! MCNAP



	  integer*4 :: nmaps,imaps, which_maps ! Marco

      real*8, dimension(nmaps,*) :: valmaps ! Marco  

      integer*4, dimension(*) :: fmaps,tagmaps ! Marco 


	  integer*4 :: vivel_length,ipoin ! Clara 
	  real*8, dimension(:,:,:), allocatable :: vivel_el ! Clara 
	  real*8, dimension(:), allocatable :: vivel_az ! Clara 
	  !integer*4, dimension(:), allocatable :: check_vivel_az ! Clara 
      integer*4, dimension(:), allocatable :: nodal_counter ! Clara 
	  real*8, dimension(:,:,:), allocatable :: vistn1_0,vistn2_0,vistn3_0 ! Clara 

        

      real*8, dimension(:,:,:), allocatable :: a1n1,a1n2,a1n3,a1n4,a1n5,a1n6 !MCNAP vpsa

      real*8, dimension(:,:,:), allocatable :: epn1,epn2,epn3,epn4,epn5,epn6 !MCNAP vpsa

	  real*8, dimension(:,:,:), allocatable :: z14n,u7n,densit               !MCNAP vpsa

	  real*8, dimension(:), allocatable :: sxx_az,syy_az,sxy_az ! Clara 


	  real*8, dimension(:), allocatable :: vivel_az_1 ! M Mexico _A_
	  real*8, dimension(:), allocatable :: vistn1_0_az, vistn2_0_az, vistn3_0_az ! M Mexico _A_

	  real*8, dimension(:), allocatable :: vistn1_0_az_1, vistn2_0_az_1, vistn3_0_az_1 ! M Mexico _A_







!*************************************************

!

!                  PRE-STRESS

!

!*************************************************



	  real*8, dimension(:,:,:), allocatable :: sxx_ps,syy_ps,sxy_ps,szz_ps

	  real*8, dimension(:), allocatable :: sxx_ps_read,syy_ps_read,sxy_ps_read,szz_ps_read ! Clara

	  integer*4 :: nprocs 



!*************************************************

!

!               INITIAL-CONDITION

!

!*************************************************



	  real*8, dimension(:), allocatable :: u1_read,u0_read ! M Mexico _A_

	  real*8, dimension(:), allocatable :: Fe_read,Fk_read ! M Mexico _A_

	  real*8, dimension(:), allocatable :: vivel_read ! M Mexico _A_

	  real*8, dimension(:), allocatable :: vistn1_0_read,vistn2_0_read,vistn3_0_read ! M Mexico _A_

	  real*8, dimension(:), allocatable :: u2_star ! M Mexico _A_

	  real*8, dimension(:), allocatable :: fe_star ! M Mexico _A_






      ne = cs(0) -1
      
      eps = 1.0d3 * dabs(epsilon(mvec(0)))
      
      allocate(func_value(nf))
      
      allocate(i4count(nnt))

	  allocate(nodal_counter(nnt))

	  do in = 1,nnt
		nodal_counter(in) = 0
	  enddo
      
      nnode_neuX = 0
      do in = 1,nnt
         i4count(in) = 0
      enddo
      
      do in = 1,nnt
         if (node_index(in).ne.0) then
            id1 = node_index(in)
            
            if (nf.gt.0) then
               r8t = 0.0d0
               do fn = 1,nf
                  r8t = r8t + dabs(Fmat(fn,id1))
               enddo
               if (r8t.gt.eps) then
                  nnode_neuX = nnode_neuX +1
                  i4count(in) = nnode_neuX
               endif
            endif
         endif
      enddo
      
      if (nnode_neuX.gt.0) then
         allocate(inode_neuX(nnode_neuX))
         
         do in = 1,nnt
            if (i4count(in).ne.0) then
               inode_neuX(i4count(in)) = in
            endif
         enddo
      endif
      
      
      nnode_neuY = 0
      do in = 1,nnt
         i4count(in) = 0
      enddo
      
      do in = 1,nnt
         if (node_index(in).ne.0) then
            id2 = node_index(in) +nnd
            
            if (nf.gt.0) then
               r8t = 0.0d0
               do fn = 1,nf
                  r8t = r8t + dabs(Fmat(fn,id2))
               enddo
               if (r8t.gt.eps) then
                  nnode_neuY = nnode_neuY +1
                  i4count(in) = nnode_neuY
               endif
            endif
         endif
      enddo
      
      if (nnode_neuY.gt.0) then
         allocate(inode_neuY(nnode_neuY))
         
         do in = 1,nnt
            if (i4count(in).ne.0) then
               inode_neuY(i4count(in)) = in
            endif
         enddo
      endif



      ! TUTTO DA CONTROLLARE: CREDO CI VOGLIA UNA COSA COME QUELLA PER LE ABC



	  nnode_neuN = 0 ! M Mexico Arroyo

      do in = 1,nnt ! M Mexico Arroyo

         i4count(in) = 0 ! M Mexico Arroyo

      enddo ! M Mexico Arroyo

      

      do in = 1,nnt ! M Mexico Arroyo

         if (node_index(in).ne.0) then ! M Mexico Arroyo

            id1 = node_index(in) ! M Mexico Arroyo

            

            if (nf.gt.0) then ! M Mexico Arroyo

               r8t = 0.0d0 ! M Mexico Arroyo

               do fn = 1,nf ! M Mexico Arroyo

                  r8t = r8t + dabs(Fmat(fn,id1)) ! M Mexico Arroyo

               enddo ! M Mexico Arroyo

               if (r8t.gt.eps) then ! M Mexico Arroyo

                  nnode_neuN = nnode_neuN +1 ! M Mexico Arroyo

                  i4count(in) = nnode_neuN ! M Mexico Arroyo

               endif ! M Mexico Arroyo

            endif ! M Mexico Arroyo

         endif ! M Mexico Arroyo

      enddo ! M Mexico Arroyo

      

      if (nnode_neuN.gt.0) then ! M Mexico Arroyo

         allocate(inode_neuN(nnode_neuN)) ! M Mexico Arroyo

         

         do in = 1,nnt ! M Mexico Arroyo

            if (i4count(in).ne.0) then ! M Mexico Arroyo

               inode_neuN(i4count(in)) = in ! M Mexico Arroyo

            endif ! M Mexico Arroyo

         enddo ! M Mexico Arroyo

      endif ! M Mexico Arroyo

      

      

      nnode_neuT = 0 ! M Mexico Arroyo

      do in = 1,nnt ! M Mexico Arroyo

         i4count(in) = 0 ! M Mexico Arroyo

      enddo ! M Mexico Arroyo

      

      do in = 1,nnt ! M Mexico Arroyo

         if (node_index(in).ne.0) then ! M Mexico Arroyo

            id2 = node_index(in) +nnd ! M Mexico Arroyo

            

            if (nf.gt.0) then ! M Mexico Arroyo

               r8t = 0.0d0 ! M Mexico Arroyo

               do fn = 1,nf ! M Mexico Arroyo

                  r8t = r8t + dabs(Fmat(fn,id2)) ! M Mexico Arroyo

               enddo ! M Mexico Arroyo

               if (r8t.gt.eps) then ! M Mexico Arroyo

                  nnode_neuT = nnode_neuT +1 ! M Mexico Arroyo

                  i4count(in) = nnode_neuT ! M Mexico Arroyo

               endif ! M Mexico Arroyo

            endif ! M Mexico Arroyo

         endif ! M Mexico Arroyo

      enddo ! M Mexico Arroyo

      

      if (nnode_neuT.gt.0) then ! M Mexico Arroyo

         allocate(inode_neuT(nnode_neuT)) ! M Mexico Arroyo

         

         do in = 1,nnt ! M Mexico Arroyo

            if (i4count(in).ne.0) then ! M Mexico Arroyo

               inode_neuT(i4count(in)) = in ! M Mexico Arroyo

            endif ! M Mexico Arroyo

         enddo ! M Mexico Arroyo

      endif ! M Mexico Arroyo
      
      
      
      nnode_dirX = 0
      do in = 1,nnt
         i4count(in) = 0
      enddo
      
      call get_edge_nodes(nnt,cs_nnz_bc,cs_bc,nl_dirX,tag_dirX,&
                          nnode_dirX,i4count)
      
      if (nnode_dirX.gt.0) then
         allocate(inode_dirX(nnode_dirX))
         
         do in = 1,nnt
            if (i4count(in).ne.0) then
               inode_dirX(i4count(in)) = in
            endif
         enddo
      endif



	  !if (nnode_dipX.gt.0) then ! Marco

      !   allocate(inode_dipX(nnode_dipX))  ! Marco

      !   

      !   i = 0

      !   do in = 1,nnt  ! Marco

      !      if (i4count(in).ne.0) then ! Marco "find_nearest"

	  !         i = 1 + 1

      !         inode_dirX(i) = in ! Marco

      !      endif ! Marco

      !   enddo ! Marco

      !endif ! Marco
      
      
      nnode_dirY = 0
      do in = 1,nnt
         i4count(in) = 0
      enddo
      
      call get_edge_nodes(nnt,cs_nnz_bc,cs_bc,nl_dirY,tag_dirY,&
                          nnode_dirY,i4count)
      
      if (nnode_dirY.gt.0) then
         allocate(inode_dirY(nnode_dirY))
         
         do in = 1,nnt
            if (i4count(in).ne.0) then
               inode_dirY(i4count(in)) = in
            endif
         enddo
      endif
      
      
      deallocate(i4count)
      
      
      
      allocate(i4count(ne))
      
      nelem_ebe = 0 
      do ie = 1,ne
         i4count(ie) = 0
      enddo
      
      do im = 1,nm
         if (type_mat(im).eq.2) then
            do ie = 1,ne
               if (cs(cs(ie -1) +0).eq.tag_mat(im)) then
                  do is = 1,(cs(ie) - cs(ie -1) -1)
                     in = cs(cs(ie -1) +is)
                     i4count(ie) = i4count(ie) + node_index(in)
                  enddo
               endif
            enddo
         endif
      enddo
      
      do ie = 1,ne
         if (i4count(ie).gt.0) then
            nelem_ebe = nelem_ebe +1
            i4count(ie) = nelem_ebe
         endif
      enddo
      
      if (nelem_ebe.gt.0) then
         allocate(ielem_ebe(nelem_ebe))
         
         do ie = 1,ne
            if (i4count(ie).ne.0) then
               ielem_ebe(i4count(ie)) = ie
            endif
         enddo
      endif
      
      deallocate(i4count)
      
      
      
      nedge_abc = 0 
      
      if (cs_nnz_bc.gt.0) then
         nedge = cs_bc(0) -1
         
         allocate(i4count(nedge))
         do iedge = 1,nedge
            i4count(iedge) = 0
         enddo
         
         if (nl_abc.gt.0) then
            do i = 1,nl_abc
               do iedge = 1,nedge
                  if (cs_bc(cs_bc(iedge -1) +0).eq.tag_abc(i)) then
                     ied1 = cs_bc(cs_bc(iedge -1) +1)
                     ied2 = cs_bc(cs_bc(iedge) -1)
                     
                     do j = 1,nelem_ebe
                        ie = ielem_ebe(j)
                        
                        nn = cs_bc(iedge) - cs_bc(iedge -1) -1
                        
                        iel1 = cs(cs(ie -1) +1)
                        iel2 = cs(cs(ie -1) +nn)
                        iel3 = cs(cs(ie -1) +nn*nn)
                        iel4 = cs(cs(ie -1) +nn*(nn -1) +1)
                        
                        if (((ied1.eq.iel1).or.(ied1.eq.iel2).or. &
                             (ied1.eq.iel3).or.(ied1.eq.iel4)).and. &
                            ((ied2.eq.iel1).or.(ied2.eq.iel2).or. &
                             (ied2.eq.iel3).or.(ied2.eq.iel4))) then
                           i4count(iedge) = iedge
                        endif
                     enddo
                  endif
               enddo
            enddo
         endif
         
         do iedge = 1,nedge
            if (i4count(iedge).gt.0) then
               nedge_abc = nedge_abc +1
               i4count(iedge) = nedge_abc
            endif
         enddo
      
         if (nedge_abc.gt.0) then
            allocate(iedge_abc(nedge_abc))
            
            do iedge = 1,nedge
               if (i4count(iedge).ne.0) then
                  iedge_abc(i4count(iedge)) = iedge
               endif
            enddo
         endif
         
         deallocate(i4count)
      endif
      
      
      nelem_abc = nedge_abc
      
      if (nelem_abc.gt.0) then
         allocate(ielem_abc(nelem_abc))
         
         do i = 1,nedge_abc
            iedge = iedge_abc(i)
            
            ied1 = cs_bc(cs_bc(iedge -1) +1)
            ied2 = cs_bc(cs_bc(iedge) -1)
            
            do j = 1,nelem_ebe
               ie = ielem_ebe(j)
               
               nn = cs_bc(iedge) - cs_bc(iedge -1) -1
               
               iel1 = cs(cs(ie -1) +1)
               iel2 = cs(cs(ie -1) +nn)
               iel3 = cs(cs(ie -1) +nn*nn)
               iel4 = cs(cs(ie -1) +nn*(nn -1) +1)
               
               if (((ied1.eq.iel1).or.(ied1.eq.iel2).or. &
                    (ied1.eq.iel3).or.(ied1.eq.iel4)).and. &
                   ((ied2.eq.iel1).or.(ied2.eq.iel2).or. &
                    (ied2.eq.iel3).or.(ied2.eq.iel4))) then
                  ielem_abc(i) = ie
               endif
            enddo
         enddo
      endif
      
      
      
      
      allocate(node_ind_ebe(nnt))
      
      nnode_ebe = 0
      do in = 1,nnt
         node_ind_ebe(in) = 0
      enddo
      
      do i = 1,nelem_ebe
         ie = ielem_ebe(i)
         do is = 1,(cs(ie) - cs(ie -1) -1)
            in = cs(cs(ie -1) +is)
            node_ind_ebe(in) = 1
         enddo
      enddo
      
      do im = 1,nm
         if (type_mat(im).ne.2) then
            do ie = 1,ne
               if (cs(cs(ie -1) +0).eq.tag_mat(im)) then
                  do is = 1,(cs(ie) - cs(ie -1) -1)
                     in = cs(cs(ie -1) +is)
                     node_ind_ebe(in) = 0
                  enddo
               endif
            enddo
         endif
      enddo
      
      do in = 1,nnt
         if (node_ind_ebe(in).ne.0) then
            if (node_index(in).ne.0) then
               nnode_ebe = nnode_ebe +1
               node_ind_ebe(in) = nnode_ebe
            else
               node_ind_ebe(in) = 0
            endif
         endif
      enddo
      
      if (nnode_ebe.gt.0) then
         allocate(inode_ebe(nnode_ebe))
         
         do in = 1,nnt
            if (node_ind_ebe(in).ne.0) then
               inode_ebe(node_ind_ebe(in)) = in
            endif
         enddo
      endif
      
      
      lname = len_trim(file_out)
      file_outE = file_out(1:lname) // 'E'
      file_outV = file_out(1:lname) // 'V' ! Clara
	  file_outT1 = file_out(1:lname) // '1' ! Clara
	  file_outT2 = file_out(1:lname) // '2' ! Clara
	  file_outT3 = file_out(1:lname) // '3' ! Clara

	  file_outS1 = file_out(1:lname) // '4' ! Clara

	  file_outS2 = file_out(1:lname) // '5' ! Clara

	  file_outS3 = file_out(1:lname) // '6' ! Clara

	  !file_outDvp1 = file_out(1:lname) // '7' ! MCNAP

	  !file_outDvp2 = file_out(1:lname) // '8' ! MCNAP

	  !file_outDvp3 = file_out(1:lname) // '9' ! MCNAP

	  file_out_u1 = file_out(1:lname) // '_u1' ! M Mexico _A_

	  file_out_u0 = file_out(1:lname) // '_u0' ! M Mexico _A_

	  file_out_Fe = file_out(1:lname) // '_Fe' ! M Mexico _A_

	  file_out_Fk = file_out(1:lname) // '_Fk' ! M Mexico _A_

	  file_out_vivel = file_out(1:lname) // '_vivel' ! M Mexico _A_

	  file_out_vistn1_0 = file_out(1:lname) // '_vistn1_0' ! M Mexico _A_

	  file_out_vistn2_0 = file_out(1:lname) // '_vistn2_0' ! M Mexico _A_

	  file_out_vistn3_0 = file_out(1:lname) // '_vistn3_0' ! M Mexico _A_
     
      file_timing = 'timingXX.d'
      if (myid.lt.10) then
         write(file_timing(7:7),'(a1)')'0'
         write(file_timing(8:8),'(i1)')myid
      else if (myid.le.99) then
         write(file_timing(7:8),'(i2)')myid
      endif
      open(40,file=file_timing)
      
      
      if (nmonit.ge.1) then 
         file_monitor = 'monitorXX.d'
         do i = 1,nmonit
            unit_monitor = 40 + i
            in = node_m(i)
            if (node_index(in).ne.0) then
               if (i.lt.10) then
                  write(file_monitor(8:8),'(a1)')'0'
                  write(file_monitor(9:9),'(i1)')i
               else if (i.le.99) then
                  write(file_monitor(8:9),'(i2)')i
               endif
               
               open(unit_monitor,file=file_monitor)
            endif
         enddo

         ! SISTEMARE CON LUCA GLI STRESS
         ! Clara - begin
		 file_stress = 'stressXX.d'
         do i = 1,nmonit
            unit_monitor = 10000 + i
            in = node_m(i)
            if (node_index(in).ne.0) then
               if (i.lt.10) then
                  write(file_stress(7:7),'(a1)')'0'
                  write(file_stress(8:8),'(i1)')i
               else if (i.le.99) then
                  write(file_stress(7:8),'(i2)')i
               endif
               
               open(unit_monitor,file=file_stress)
            endif
         enddo
		 ! Clara  - end
      endif
      
      
      
!**************************************************
!  
!  AZTEC INITIALIZATION
!
!*************************************************
     
      N_update_el_az = 2*nnt
      N_extern_el_az = 0
      
      allocate(update_el_az(0:(N_update_el_az -1)))
      do in = 1,nnt
         if (node_index(in).ne.0) then
            id1 = node_index(in)
            id2 = node_index(in) +nnd
            
            update_el_az(id1 -1) = in -1
            update_el_az(id2 -1) = in +nnt -1
         endif
      enddo
      
      allocate(update_index_el_az(0:N_update_el_az -1))
      allocate(global_index_el_az(2*nnt))
     
      do iaz = 0,N_update_el_az -1
         update_index_el_az(iaz) = iaz
      enddo
      
      do i = 1,2*nnt
         global_index_el_az(i) = 0
      enddo
      
      do iaz = 0,N_update_el_az -1
         global_index_el_az(update_el_az(iaz) +1) = update_index_el_az(iaz)
      enddo
      
      
      allocate(out_index_el_az(0:(N_update_el_az -1)))
      
      do iaz = 0,N_update_el_az -1
         out_index_el_az(update_index_el_az(iaz)) = update_el_az(iaz) +1
      enddo
      
      
      
!**************************************************
!  
!     UNKNOWNS INITIALIZATION
!
!*************************************************
      
      
      allocate(u_1_az(0:(N_update_el_az + N_extern_el_az -1)))
      allocate(u0_az(0:(N_update_el_az + N_extern_el_az -1)))
      allocate(u1_az(0:(N_update_el_az + N_extern_el_az -1)))
      allocate(u2_az(0:(N_update_el_az + N_extern_el_az -1)))
      allocate(fk_az(0:(N_update_el_az + N_extern_el_az -1)))
      allocate(fe_az(0:(N_update_el_az -1)))
	  allocate(vivel_az(0:(N_update_el_az -1))) ! Clara 
	  allocate(sxx_az(0:(N_update_el_az -1))) ! Clara 
	  allocate(syy_az(0:(N_update_el_az -1))) ! Clara 
	  allocate(sxy_az(0:(N_update_el_az -1))) ! Clara
	  allocate(vistn1_0_az(0:(N_update_el_az -1))) ! Clara
	  allocate(vistn2_0_az(0:(N_update_el_az -1))) ! Clara
	  allocate(vistn3_0_az(0:(N_update_el_az -1))) ! Clara

	  allocate(vivel_az_1(0:(N_update_el_az -1))) ! M Mexico _A_

	  allocate(vistn1_0_az_1(0:(N_update_el_az -1))) ! M Mexico _A_

	  allocate(vistn2_0_az_1(0:(N_update_el_az -1))) ! M Mexico _A_

	  allocate(vistn3_0_az_1(0:(N_update_el_az -1))) ! M Mexico _A_

	  allocate(sxx_ps_read(0:(N_update_el_az -1))) ! Marco 

	  allocate(syy_ps_read(0:(N_update_el_az -1))) ! Marco 

	  allocate(sxy_ps_read(0:(N_update_el_az -1))) ! Marco

      allocate(szz_ps_read(0:(N_update_el_az -1))) ! Clara 12/01/2005



	  allocate(u1_read(0:(N_update_el_az + N_extern_el_az -1))) ! M Mexico _A_

	  allocate(u0_read(0:(N_update_el_az + N_extern_el_az -1))) ! M Mexico _A_

	  allocate(Fe_read(0:(N_update_el_az + N_extern_el_az -1))) ! M Mexico _A_

	  allocate(Fk_read(0:(N_update_el_az + N_extern_el_az -1))) ! M Mexico _A_

	  allocate(vivel_read(0:(N_update_el_az + N_extern_el_az -1))) ! M Mexico _A_

	  allocate(vistn1_0_read(0:(N_update_el_az + N_extern_el_az -1))) ! M Mexico _A_
      allocate(vistn2_0_read(0:(N_update_el_az + N_extern_el_az -1))) ! M Mexico _A_

	  allocate(vistn3_0_read(0:(N_update_el_az + N_extern_el_az -1))) ! M Mexico _A_



	  allocate(u2_star(0:(N_update_el_az + N_extern_el_az -1))) ! M Mexico _A_

	  allocate(fe_star(0:(N_update_el_az + N_extern_el_az -1))) ! M Mexico _A_




      
      do iaz = 0,N_update_el_az -1
         fe_az(iaz) = 0.0d0

		 sxx_ps_read(iaz) = 0.0d0 ! Marco

		 syy_ps_read(iaz) = 0.0d0 ! Marco

		 sxy_ps_read(iaz) = 0.0d0 ! Marco

         szz_ps_read(iaz) = 0.0d0 ! Clara 12/01/2005

		 

		 u1_read(iaz) = 0.0d0 ! M Mexico _A_

		 u0_read(iaz) = 0.0d0 ! M Mexico _A_

		 Fe_read(iaz) = 0.0d0 ! M Mexico _A_

		 Fk_read(iaz) = 0.0d0 ! M Mexico _A_

		 vivel_read(iaz) = 0.0d0 ! M Mexico _A_

		 vistn1_0_read(iaz) = 0.0d0 ! M Mexico _A_

		 vistn2_0_read(iaz) = 0.0d0 ! M Mexico _A_

		 vistn3_0_read(iaz) = 0.0d0 ! M Mexico _A_



		 u2_star(iaz) = 0.0d0 ! M Mexico _A_

		 fe_star(iaz) = 0.0d0 ! M Mexico _A_




      enddo
      
      dt2 = dt*dt
      
      do i = 0,N_update_el_az -1
         mvec(i) = mvec(i) / dt2
         if (make_damping_yes_or_not.gt.0) then
            cvec(i) = cvec(i) / ( 2 * dt )
            new(i) =  mvec(i) + cvec(i)
         endif
         
      enddo
      
      do i = 0,N_update_el_az -1
         u_1_az(i) = u_1(i +1)
         u0_az(i) = u0(i +1)
      enddo
      
      
      nn = sdeg_mat(1) +1
      allocate(ct(nn),ww(nn),dd(nn,nn))
      allocate(dxdx_el(nn),dxdy_el(nn),dydx_el(nn),dydy_el(nn))
      allocate(ux_el(nn,nn),uy_el(nn,nn))
      allocate(vx_el(nn,nn),vy_el(nn,nn))
      allocate(duxdx_el(nn,nn),duxdy_el(nn,nn),duydx_el(nn,nn),duydy_el(nn,nn))
      allocate(sxx_el(nn,nn),syy_el(nn,nn),sxy_el(nn,nn))
      allocate(fx_el(nn,nn),fy_el(nn,nn))
      
      call lgl(nn,ct,ww,dd)

      
!**************************************************
!  
!     VISCO-PLASTIC INITIALIZATION 
!
!*************************************************	  
	  
	  vivel_length = nnt ! Clara 
	  allocate(vivel_el(nn,nn,nelem_ebe)) ! Clara 
	  allocate(vistn1_0(nn,nn,nelem_ebe)) ! Clara 
	  allocate(vistn2_0(nn,nn,nelem_ebe)) ! Clara 
	  allocate(vistn3_0(nn,nn,nelem_ebe)) ! Clara 



      allocate(a1n1(nn,nn,nelem_ebe))     !MCNAP vpsa

      allocate(a1n2(nn,nn,nelem_ebe))     !MCNAP vpsa

	  allocate(a1n3(nn,nn,nelem_ebe))     !MCNAP vpsa

	  allocate(a1n4(nn,nn,nelem_ebe))     !MCNAP vpsa

	  allocate(a1n5(nn,nn,nelem_ebe))     !MCNAP vpsa

	  allocate(a1n6(nn,nn,nelem_ebe))     !MCNAP vpsa



      allocate(epn1(nn,nn,nelem_ebe))     !MCNAP vpsa

      allocate(epn2(nn,nn,nelem_ebe))     !MCNAP vpsa

	  allocate(epn3(nn,nn,nelem_ebe))     !MCNAP vpsa

	  allocate(epn4(nn,nn,nelem_ebe))     !MCNAP vpsa

	  allocate(epn5(nn,nn,nelem_ebe))     !MCNAP vpsa

	  allocate(epn6(nn,nn,nelem_ebe))     !MCNAP vpsa

  

      allocate(z14n(nn,nn,nelem_ebe))     !MCNAP vpsa

	  allocate(u7n(nn,nn,nelem_ebe))      !MCNAP vpsa

      allocate(densit(nn,nn,nelem_ebe))   !MCNAP vpsa







	  

!*************************************************

!

!                      PRE-STRESS

!

!*************************************************





	 allocate(sxx_ps(nn,nn,nelem_ebe)) ! Marco

	 allocate(syy_ps(nn,nn,nelem_ebe)) ! Marco

	 allocate(sxy_ps(nn,nn,nelem_ebe)) ! Marco

     allocate(szz_ps(nn,nn,nelem_ebe)) ! Clara 12/01/2005



      
      
!**************************************************
!  
!     FIRST STEP
!
!*************************************************



      write(*,'(A)')'FIRST STEP'
      
      isnap = 1

	  isnapps = 1 ! Marco
      
      time_total = 0.d0

	  do i = 1,nn ! Clara 
		do j = 1,nn ! Clara 
			do ie = 1,nelem_ebe ! Clara 



				vistn1_0(i,j,ie) = 0.0d0 ! Clara 

				vistn2_0(i,j,ie) = 0.0d0 ! Clara 

				vistn3_0(i,j,ie) = 0.0d0 ! Clara 

				vivel_el(i,j,ie) = 0.0d0 ! Clara

				sxx_ps(i,j,ie) = 0.0d0 ! Marco

				syy_ps(i,j,ie) = 0.0d0 ! Marco

				sxy_ps(i,j,ie) = 0.0d0 ! Marco

                szz_ps(i,j,ie) = 0.0d0 ! Clara 12/01/2005



				if (nelem_ebe.gt.0) then



					do im = 1,nm



					! NOT VERY EFFICIENT CHOICE!!!

					! SHOULD BE CHANGED!!!

               

					which_vpsa = 0 ! Clara 



					do ivpsa = 1,nvpsa ! Clara 

						if (tag_mat(im).eq.tagvpsa(ivpsa)) then ! Clara 

							which_vpsa = ivpsa ! Clara 



							!isotropic initialization of yield function axes

							a1n1(i,j,ie) =  0.00001773             !MCNAP vpsa

							a1n2(i,j,ie) =  0.00001773             !MCNAP vpsa

							a1n3(i,j,ie) =  0.00001773             !MCNAP vpsa

							a1n4(i,j,ie) =  0.0d0                  !MCNAP vpsa

							a1n5(i,j,ie) =  0.0d0                  !MCNAP vpsa

							a1n6(i,j,ie) =  0.0d0                  !MCNAP vpsa

							!initialization of viscoplastic deformation!!

							!controllare se non è la stessa cosa usare le vistn1_0....

							epn1(i,j,ie) =  0.0d0                  !MCNAP vpsa

							epn2(i,j,ie) =  0.0d0                  !MCNAP vpsa

							epn3(i,j,ie) =  0.0d0                  !MCNAP vpsa

							epn4(i,j,ie) =  0.0d0                  !MCNAP vpsa

							epn5(i,j,ie) =  0.0d0                  !MCNAP vpsa

							epn6(i,j,ie) =  0.0d0                  !MCNAP vpsa

							!initialization of plastic deformation accumulated before

							z14n(i,j,ie) = 0.0d0                   !MCNAP vpsa

							!initialization of rc (=sqrt(3.)*pc)

							!!CONTROLLARE PERCHè NON DOVREBBE ESSERE ZERO

							u7n(i,j,ie) = 0.0d0                    !MCNAP vpsa

							!initialization of densità

							densit(i,j,ie)=valvpsa(which_vpsa,1)   !MCNAP vpsa

				

						endif   !MCNAP vpsa

					enddo   !MCNAP vpsa



				enddo !im

			  endif !nelem_ebe

				


			enddo ! Clara 
		enddo ! Clara 
	  enddo ! Clara 





	  write(*,'(A)')'FIRST STEP1'





	  if (nmaps.gt.0) then





		do im = 1,nm



		if (type_mat(im).eq.2) then



			do k = 1,nelem_ebe

				ie = ielem_ebe(k)



				which_maps = 0 ! Marco 

				do imaps = 1,nmaps ! Marco

			   

					if (fmaps(imaps).le.3) then ! M Mexico _A_



						if (tag_mat(im).eq.tagmaps(imaps)) then ! Marco 

							which_maps = imaps ! Marco 



							call make_pre_stress(prop_mat(imaps,1),prop_mat(imaps,2),prop_mat(imaps,3),&

					 					  valmaps(imaps,1),fmaps(imaps),&

					 					  sxx_ps,syy_ps,sxy_ps,szz_ps,&

					 					  nn,ie,nelem_ebe,&

					 				   	  nnt,xs,ys,&

					 					  cs_nnz,cs,node_index,&

										  N_update_el_az,update_index_el_az,&

								          file_out,nsnapps,itersnapps,&

						                  sxx_ps_read,syy_ps_read,sxy_ps_read,szz_ps_read,&

										  !vistnxx_ps_read,vistnyy_ps_read,vistnxy_ps_read,& ! MCNAP

										  nprocs)

						endif ! Marco					 

                

					endif ! M Mexico _A_

 

				enddo !imaps ! Marco



			   enddo ! k

			endif



 

		enddo ! im



		do imaps = 1,nmaps ! Marco

			if (fmaps(imaps).eq.4) then ! M Mexico _A_

					write(*,'(A)')'FIRST STEP1a'



					    call read_initial_condition(nnt,file_out,nsnapps,itersnapps,& ! M Mexico _A_

								 u1_read,u0_read,& ! M Mexico _A_

								 Fe_read,Fk_read,& ! M Mexico _A_

								 vivel_read,& ! M Mexico _A_

								 vistn1_0_read,vistn2_0_read,vistn3_0_read,& ! M Mexico _A_

								 nprocs) ! M Mexico _A_	



				     write(*,'(A)')'FIRST STEP1b'

			endif

		enddo



	  endif



	  write(*,'(A)')'FIRST STEP2'







         if (nmaps.gt.0) then ! M Mexico _A_

			do imaps=1,nmaps ! M Mexico _A_

				if (fmaps(nmaps).eq.4) then ! M Mexico _A_

					

					do iaz = 0,N_update_el_az -1

					! 

					! I don't like this solution!!! Not general: 

					! It should be modified ASAP!!!

					! si dovrebbe poter pre-stressare solo alcuni materiali e non tutti!!!



						u0_az(iaz)  = u1_read(iaz) ! M Mexico _A_

						u_1_az(iaz) = u0_read(iaz) ! M Mexico _A_

						fe_az(iaz)  = Fe_read(iaz) ! M Mexico _A_

						fk_az(iaz)  = Fk_read(iaz) ! M Mexico _A_

						vivel_az(iaz)  = vivel_read(iaz) ! M Mexico _A_

						vistn1_0_az(iaz)  = vistn1_0_read(iaz) ! M Mexico _A_

						vistn2_0_az(iaz)  = vistn2_0_read(iaz) ! M Mexico _A_

						vistn3_0_az(iaz)  = vistn3_0_read(iaz) ! M Mexico _A_



					enddo ! iaz



				endif ! M Mexico _A_

			enddo ! M Mexico _A_

		 endif ! M Mexico _A_





	  write(*,'(A)')'FIRST STEP3'





      if (nmaps.gt.0) then ! M Mexico _A_

		do imaps=1,nmaps ! M Mexico _A_

			if (fmaps(nmaps).eq.4) then ! M Mexico _A_

      

	  		do i = 1,nn ! M Mexico _A_

	  			do j = 1,nn ! M Mexico _A_

	  				do ie = 1,nelem_ebe ! M Mexico _A_

      

	  					is = nn*(j -1) +i ! M Mexico _A_

	  					in = cs(cs(ie -1) + is) ! M Mexico _A_

	  			

	  					id1 = node_index(in) ! M Mexico _A_

	  					iaz = update_index_el_az(id1 -1) ! M Mexico _A_



						vivel_el(i,j,ie) = vivel_az(iaz) ! M Mexico _A_



	  					vistn1_0(i,j,ie) = vistn1_0_az(iaz) ! M Mexico _A_

	  					vistn2_0(i,j,ie) = vistn2_0_az(iaz) ! M Mexico _A_

	  					vistn3_0(i,j,ie) = vistn3_0_az(iaz) ! M Mexico _A_

						 

	  				enddo ! M Mexico _A_

	  			enddo ! M Mexico _A_

	  		enddo ! M Mexico _A_



			endif  ! M Mexico _A_ 

		enddo

	  endif



	  write(*,'(A)')'FIRST STEP4'

      
      do fn = 1,nf
         func_value(fn) = get_func_value(nf,func_type,func_indx,func_data,nfunc_data,&
                                         fn,0.0d0,0.0d0)
      enddo
      
      
      if (nnode_neuX.gt.0) then
         do i = 1,nnode_neuX
            in = inode_neuX(i)
            
            if (node_index(in).ne.0) then
               id1 = node_index(in)
               iaz = update_index_el_az(id1 -1)



			   !fe_az(iaz) = 0.0d0




               !if (nmaps.gt.0) then ! M Mexico _A_

				!	do imaps=1,nmaps ! M Mexico _A_

				!		if (fmaps(nmaps).eq.4) then ! M Mexico _A_
				!			fe_az(iaz) = Fe_read(iaz) ! M Mexico _A_

				!		endif ! M Mexico _A_

				!	enddo ! M Mexico _A_

				!endif ! M Mexico _A_





               fe_az(iaz) = Fe_read(iaz) ! M Mexico _A_
               do fn = 1,nf
                  fe_az(iaz) = fe_az(iaz) + Fmat(fn,id1) * func_value(fn)
               enddo
            endif
         enddo
      endif
      
      if (nnode_neuY.gt.0) then
         do i = 1,nnode_neuY
            in = inode_neuY(i)
            
            if (node_index(in).ne.0) then
               id2 = node_index(in) +nnd
               iaz = update_index_el_az(id2 -1)


               !fe_az(iaz) = 0.0d0


               !if (nmaps.gt.0) then ! M Mexico _A_

				!	do imaps=1,nmaps ! M Mexico _A_

				!		if (fmaps(nmaps).eq.4) then ! M Mexico _A_

				!			fe_az(iaz) = Fe_read(iaz) ! M Mexico _A_

				!		endif ! M Mexico _A_ 

				!	enddo ! M Mexico _A_

				!endif ! M Mexico _A_



			   fe_az(iaz) = Fe_read(iaz) ! M Mexico _A_		   
               do fn = 1,nf
                  fe_az(iaz) = fe_az(iaz) + Fmat(fn,id2) * func_value(fn)
               enddo
            endif
         enddo
      endif





	if (nnode_neuN.gt.0) then ! M Mexico Arroyo

         do i = 1,nnode_neuN ! M Mexico Arroyo

            in = inode_neuN(i) ! M Mexico Arroyo

            

            if (node_index(in).ne.0) then ! M Mexico Arroyo

               

			   id1 = node_index(in) ! M Mexico Arroyo

               iaz = update_index_el_az(id1 -1) ! M Mexico Arroy

               fe_az(iaz) = Fe_read(iaz) ! M Mexico Arroyo

               do fn = 1,nf ! M Mexico Arroyo

                  fe_az(iaz) = fe_az(iaz) + Fmat(fn,id1) * func_value(fn) ! M Mexico Arroyo

               enddo ! M Mexico Arroyo



			   id2 = node_index(in) +nnd ! M Mexico Arroyo

               iaz = update_index_el_az(id2 -1) ! M Mexico Arroyo

			   fe_az(iaz) = Fe_read(iaz) ! M Mexico Arroyo  

               do fn = 1,nf ! M Mexico Arroyo

                  fe_az(iaz) = fe_az(iaz) + Fmat(fn,id2) * func_value(fn) ! M Mexico Arroyo

               enddo ! M Mexico Arroyo



            endif ! M Mexico Arroyo

         enddo ! M Mexico Arroyo

      endif ! M Mexico Arroyo

      

      if (nnode_neuT.gt.0) then ! M Mexico Arroyo

         do i = 1,nnode_neuT ! M Mexico Arroyo

            in = inode_neuT(i) ! M Mexico Arroyo

            

            if (node_index(in).ne.0) then ! M Mexico Arroyo

               

			   id1 = node_index(in) ! M Mexico Arroyo

               iaz = update_index_el_az(id1 -1) ! M Mexico Arroy

               fe_az(iaz) = Fe_read(iaz) ! M Mexico Arroyo

               do fn = 1,nf ! M Mexico Arroyo

                  fe_az(iaz) = fe_az(iaz) + Fmat(fn,id1) * func_value(fn) ! M Mexico Arroyo

               enddo ! M Mexico Arroyo



			   id2 = node_index(in) +nnd ! M Mexico Arroyo

               iaz = update_index_el_az(id2 -1) ! M Mexico Arroyo

			   fe_az(iaz) = Fe_read(iaz) ! M Mexico Arroyo  

               do fn = 1,nf ! M Mexico Arroyo

                  fe_az(iaz) = fe_az(iaz) + Fmat(fn,id2) * func_value(fn) ! M Mexico Arroyo

               enddo ! M Mexico Arroyo



            endif ! M Mexico Arroyo

         enddo ! M Mexico Arroyo

      endif ! M Mexico Arroyo



	  write(*,'(A)')'FIRST STEP4'
      
      do iaz = 0,N_update_el_az -1
         u_1_az(iaz) = u0_az(iaz)
      enddo
      
      
      do iaz = 0,N_update_el_az -1

		!fk_az(iaz) = 0.0d0

		!if (nmaps.gt.0) then ! M Mexico _A_

		!	do imaps=1,nmaps ! M Mexico _A_

		!		if (fmaps(nmaps).eq.4) then ! M Mexico _A_

		!			fk_az(iaz) = Fk_read(iaz) ! M Mexico _A_

		!		endif ! M Mexico _A_

		!	enddo ! M Mexico _A_

		 !endif ! M Mexico _A

		 fk_az(iaz) = Fk_read(iaz) ! M Mexico _A_
		 vivel_az(iaz) = 0.0d0 ! Clara 
		 sxx_az(iaz) = 0.0d0 ! Clara 
		 syy_az(iaz) = 0.0d0 ! Clara 
		 sxy_az(iaz) = 0.0d0 ! Clara

		 vistn1_0_az(iaz)  = 0.0d0 ! M Mexico _A_

		 vistn2_0_az(iaz)  = 0.0d0 ! M Mexico _A_ 

		 vistn3_0_az(iaz)  = 0.0d0 ! M Mexico _A_ 

	  enddo





	  !do iaz = 0,N_update_el_az -1 ! M Mexico _A_

      !   if (nmaps.gt.0) then ! M Mexico _A_

	  !		do imaps=1,nmaps ! M Mexico _A_

	  !			if (fmaps(nmaps).eq.4) then ! M Mexico _A_

      !

	  !				!u0_az(iaz)  = u1_read(iaz) ! M Mexico _A_

	  !				u_1_az(iaz) = u0_read(iaz) ! M Mexico _A_

	  !				!fk_az(iaz)  = Fk_read(iaz) ! M Mexico _A_

	  !				!fe_az(iaz)  = 0.0d0 ! M Mexico _A_

	  !				!fe_az(iaz)  = Fe_read(iaz) ! M Mexico _A_

	  !				vivel_az(iaz)  = vivel_read(iaz) ! M Mexico _A_

	  !				vistn1_0_az(iaz)  = vistn1_0_read(iaz) ! M Mexico _A_

	  !				vistn2_0_az(iaz)  = vistn2_0_read(iaz) ! M Mexico _A_

	  !				vistn3_0_az(iaz)  = vistn3_0_read(iaz) ! M Mexico _A_

      !

	  !		    endif ! M Mexico _A_

	  !	   enddo ! M Mexico _A_

	  !	 endif ! M Mexico _A_

	  !enddo ! iaz ! M Mexico _A_



	  write(*,'(A)')'FIRST STEP5'


      
      if (make_damping_yes_or_not.eq.0) then
         do iaz = 0,N_update_el_az -1
            u1_az(iaz) = 2.0d0 * u0_az(iaz) - u_1_az(iaz) &
                         + (fe_az(iaz) - fk_az(iaz)) / mvec(iaz)
         enddo
      else
         do iaz = 0,N_update_el_az -1
         u1_az(iaz) = ( fe_az(iaz) - fk_az(iaz) - kcvec(iaz) * u0_az(iaz) &
                       + cvec(iaz) * u_1_az(iaz) + mvec(iaz) &
                       * ( 2.0d0 * u0_az(iaz) - u_1_az(iaz) ) ) &
                       / new(iaz)
         enddo
      endif



     



		if (nmaps.gt.0) then ! M Mexico _A_

			do imaps=1,nmaps ! M Mexico _A_

				if (fmaps(nmaps).eq.4) then ! M Mexico _A_



				    do iaz = 0,N_update_el_az -1

						u2_star(iaz) = u0_az(iaz) ! M Mexico _A_

						fe_star(iaz) = fe_az(iaz) ! M Mexico _A_

					enddo



				endif ! M Mexico _A_

			enddo ! M Mexico _A_

		 endif ! M Mexico _A_





	write(*,'(A)')'FIRST STEP6'






      
      do fn = 1,nf
         func_value(fn) = get_func_value(nf,func_type,func_indx,func_data,nfunc_data,&
                                         fn,dt,0.d0)
      enddo
      
!     Dirichlet b.c.
      if (nnode_dirX.gt.0) then
         do i = 1,nnode_dirX
            in = inode_dirX(i)
            
            if (node_index(in).ne.0) then
               id1 = node_index(in)
               iaz = update_index_el_az(id1 -1)
               

			   u1_az(iaz) = u2_star(iaz) ! M Mexico _A_
               !u1_az(iaz) = 0.0d0
               do fn = 1,nf
                  u1_az(iaz) = u1_az(iaz) + Fmat(fn,id1) * func_value(fn)
               enddo
            endif
         enddo
      endif



!     Dirichlet b.c. ! Marco

      if (nl_dipX.gt.0) then ! Marco

         do i = 1,nl_dipX ! Marco

            in = in_dipX(i) ! Marco

            

            if (node_index(in).ne.0) then ! Marco

               id1 = node_index(in) ! Marco

               iaz = update_index_el_az(id1 -1) ! Marco



               u1_az(iaz) = u2_star(iaz) ! M Mexico _A_

               !u1_az(iaz) = 0.0d0 ! Marco

               do fn = 1,nf ! Marco

                  u1_az(iaz) = u1_az(iaz) + Fmat(fn,id1) * func_value(fn) ! Marco

               enddo ! Marco

            endif ! Marco

         enddo ! Marco

      endif ! Marco


      
      if (nnode_dirY.gt.0) then
         do i = 1,nnode_dirY
            in = inode_dirY(i)
            
            if (node_index(in).ne.0) then
               id2 = node_index(in) +nnd
               iaz = update_index_el_az(id2 -1)



			   u1_az(iaz) = u2_star(iaz) ! M Mexico _A_
               !u1_az(iaz) = 0.0d0
               do fn = 1,nf
                  u1_az(iaz) = u1_az(iaz) + Fmat(fn,id2) * func_value(fn)
               enddo
            endif
         enddo
      endif



	  if (nl_dipY.gt.0) then ! Marco

         do i = 1,nl_dipY ! Marco

            in = in_dipY(i) ! Marco

            

            if (node_index(in).ne.0) then ! Marco

               id2 = node_index(in) +nnd ! Marco

               iaz = update_index_el_az(id2 -1) ! Marco

               

			   u1_az(iaz) = u2_star(iaz) ! M Mexico _A_

               !u1_az(iaz) = 0.0d0 ! Marco

               do fn = 1,nf ! Marco

                  u1_az(iaz) = u1_az(iaz) + Fmat(fn,id2) * func_value(fn) ! Marco

               enddo ! Marco

            endif ! Marco

         enddo ! Marco

      endif ! Marco
!     End Dirichlet b.c.



write(*,'(A)')'FIRST STEP7'






      
      
      
! *********************************************************         
      
      if ((nsnap.gt.0).and.(isnap.le.nsnap)) then
         if (its.ge.itersnap(isnap)) then
            call write_bin_file(file_outE,isnap,myid,&
                                N_update_el_az,out_index_el_az,u1_az)

			call write_bin_file(file_outV,isnap,myid,& ! Clara 
                                N_update_el_az,out_index_el_az,vivel_az) ! Clara

			call write_bin_file(file_outT1,isnap,myid,& ! Clara 
                                N_update_el_az,out_index_el_az,vistn1_0_az) ! Clara

			call write_bin_file(file_outT2,isnap,myid,& ! Clara 
                                N_update_el_az,out_index_el_az,vistn2_0_az) ! Clara

			call write_bin_file(file_outT3,isnap,myid,& ! Clara 
                                N_update_el_az,out_index_el_az,vistn3_0_az) ! Clara

            isnap = isnap + 1
            
            if (myid.eq.0) write(*,'(A)')'PRESTRESS'
         endif
      endif



	  if ((nsnapps.gt.0).and.(isnapps.le.nsnapps)) then

		if (its.ge.itersnapps(isnapps)) then



			call write_bin_file(file_outS1,isnapps,myid,& ! Clara 

                                N_update_el_az,out_index_el_az,sxx_az) ! Clara

		    call write_bin_file(file_outS2,isnapps,myid,& ! Clara 

                                N_update_el_az,out_index_el_az,syy_az) ! Clara

            call write_bin_file(file_outS3,isnapps,myid,& ! Clara 

                                N_update_el_az,out_index_el_az,sxy_az) ! Clara



			isnapps = isnapps + 1



		endif

	  endif
      

      
      
!**************************************************
!  
!     ALL STEPS
!
!*************************************************



write(*,'(A)')'ALL STEPS'
      

    ! SISTEMARE - begin

    do ie = 1,ne ! Clara 

		do j = 1,nn ! Clara 

			do i = 1,nn ! Clara 

				is = nn*(j -1) +i ! Clara 

				in = cs(cs(ie -1) + is) ! Clara 

				nodal_counter(in) = nodal_counter(in)+1 ! Clara 

			enddo ! Clara 

		enddo ! Clara 

	enddo ! Clara 


      do its = 1,nts



	     !write(*,*)its

		 ipoin = 0
         

         call system_clock(COUNT=clock_start,COUNT_RATE=clock(2))
         
         tt0 = dfloat(its -1) * dt
         tt1 = dfloat(its) * dt
         tt2 = dfloat(its +1) * dt
         

		do iaz = 0,N_update_el_az -1



			vivel_az_1(iaz) = vivel_az(iaz)! M Mexico _A_

			vistn1_0_az_1(iaz) = vistn1_0_az(iaz)! M Mexico _A_

	  		vistn2_0_az_1(iaz) = vistn2_0_az(iaz)! M Mexico _A_

	  		vistn3_0_az_1(iaz) = vistn3_0_az(iaz)! M Mexico _A_

			

			vivel_az(iaz) = 0.0d0! M Mexico _A_ 

			vistn1_0_az(iaz)  = 0.0d0 ! M Mexico _A_ 

			vistn2_0_az(iaz)  = 0.0d0 ! M Mexico _A_

			vistn3_0_az(iaz)  = 0.0d0 ! M Mexico _A_ 





		enddo ! M Mexico _A_


         
         if (nnode_neuX.gt.0) then
            do i = 1,nnode_neuX
               in = inode_neuX(i)
               
               if (node_index(in).ne.0) then
                  id1 = node_index(in)
                  iaz = update_index_el_az(id1 -1)
                  
                  !fe_az(iaz) = 0.0d0

				  fe_az(iaz) = fe_star(iaz) ! M Mexico _A_
                  do fn = 1,nf
                     fe_az(iaz) = fe_az(iaz) + Fmat(fn,id1) * func_value(fn)
                  enddo
               endif
            enddo
         endif
         
         if (nnode_neuY.gt.0) then
            do i = 1,nnode_neuY
               in = inode_neuY(i)
               
               if (node_index(in).ne.0) then
                  id2 = node_index(in) +nnd
                  iaz = update_index_el_az(id2 -1)
                  
                  

				  !fe_az(iaz) = 0.0d0

				  fe_az(iaz) = fe_star(iaz) ! M Mexico _A_
                  do fn = 1,nf
                     fe_az(iaz) = fe_az(iaz) + Fmat(fn,id2) * func_value(fn)
                  enddo
               endif
            enddo
         endif



		 if (nnode_neuN.gt.0) then ! M Mexico Arroyo

            do i = 1,nnode_neuN ! M Mexico Arroyo

               in = inode_neuN(i) ! M Mexico Arroyo

               

               if (node_index(in).ne.0) then ! M Mexico Arroyo

                  

				  id1 = node_index(in) ! M Mexico Arroyo

                  iaz = update_index_el_az(id1 -1) ! M Mexico Arroyo

                  !fe_az(iaz) = 0.0d0 ! M Mexico Arroyo

				  fe_az(iaz) = fe_star(iaz) ! M Mexico Arroyo

                  do fn = 1,nf ! M Mexico Arroyo

                     fe_az(iaz) = fe_az(iaz) + Fmat(fn,id1) * func_value(fn) ! M Mexico Arroyo

                  enddo ! M Mexico Arroyo



				  id2 = node_index(in) +nnd ! M Mexico Arroyo

                  iaz = update_index_el_az(id2 -1) ! M Mexico Arroyo                 

				  !fe_az(iaz) = 0.0d0 ! M Mexico Arroyo

				  fe_az(iaz) = fe_star(iaz) ! M Mexico Arroyo

                  do fn = 1,nf ! M Mexico Arroyo

                     fe_az(iaz) = fe_az(iaz) + Fmat(fn,id2) * func_value(fn) ! M Mexico Arroyo

                  enddo ! M Mexico Arroyo



               endif ! M Mexico Arroyo

            enddo ! M Mexico Arroyo

         endif ! M Mexico Arroyo

         

         if (nnode_neuT.gt.0) then ! M Mexico Arroyo

            do i = 1,nnode_neuT ! M Mexico Arroyo

               in = inode_neuT(i) ! M Mexico Arroyo

               

               if (node_index(in).ne.0) then ! M Mexico Arroyo

                  

				  id1 = node_index(in) ! M Mexico Arroyo

                  iaz = update_index_el_az(id1 -1) ! M Mexico Arroyo

                  !fe_az(iaz) = 0.0d0 ! M Mexico Arroyo

				  fe_az(iaz) = fe_star(iaz) ! M Mexico Arroyo

                  do fn = 1,nf ! M Mexico Arroyo

                     fe_az(iaz) = fe_az(iaz) + Fmat(fn,id1) * func_value(fn) ! M Mexico Arroyo

                  enddo ! M Mexico Arroyo



				  id2 = node_index(in) +nnd ! M Mexico Arroyo

                  iaz = update_index_el_az(id2 -1) ! M Mexico Arroyo                 

				  !fe_az(iaz) = 0.0d0 ! M Mexico Arroyo

				  fe_az(iaz) = fe_star(iaz) ! M Mexico Arroyo

                  do fn = 1,nf ! M Mexico Arroyo

                     fe_az(iaz) = fe_az(iaz) + Fmat(fn,id2) * func_value(fn) ! M Mexico Arroyo

                  enddo ! M Mexico Arroyo



               endif ! M Mexico Arroyo

            enddo ! M Mexico Arroyo

         endif ! M Mexico Arroyo
         
         
         if (nnode_ebe.gt.0) then
            do i = 1,nnode_ebe
               in = inode_ebe(i)
               
               if (node_index(in).ne.0) then
                  id1 = node_index(in)
                  iaz = update_index_el_az(id1 -1)
                  fk_az(iaz) = 0.0d0
				  vivel_az(iaz) = 0.0d0 ! Clara 
				  sxx_az(iaz) = 0.0d0 ! Clara 
				  syy_az(iaz) = 0.0d0 ! Clara 
				  sxy_az(iaz) = 0.0d0 ! Clara





				  ! INUTILE GIA' FATTO AL PASSO PRECEDENTE!!!

				  !if (nmaps.eq.0) then ! MCNAP

				  !	 vistn1_0_az(iaz) = 0.0d0 ! Clara 

			      !   vistn2_0_az(iaz) = 0.0d0 ! Clara 

			      !   vistn3_0_az(iaz) = 0.0d0 ! Clara 

		          !elseif (nmaps.gt.0) then ! MCNAP

			      !   vistn1_0_az(iaz) = vistnxx_ps_read(iaz) ! MCNAP

			      !   vistn2_0_az(iaz) = vistnyy_ps_read(iaz) ! MCNAP

			      !   vistn3_0_az(iaz) = vistnxy_ps_read(iaz) ! MCNAP

		          !endif 
                  
                  id2 = node_index(in) +nnd
                  iaz = update_index_el_az(id2 -1)
                  fk_az(iaz) = 0.0d0
               endif
            enddo
         endif
         
         if (nelem_ebe.gt.0) then
            do im = 1,nm

			   ! NOT VERY EFFICIENT CHOICE!!!
			   ! SHOULD BE CHANGED!!!
               
			   which_vpcl = 0 ! Clara 
			   do ivpcl = 1,nvpcl ! Clara 
                  if (tag_mat(im).eq.tagvpcl(ivpcl)) then ! Clara 
				     which_vpcl = ivpcl ! Clara 
				  endif ! Clara 
			   enddo ! Clara



			   which_vpsa = 0 ! MCNAP 

			   do ivpsa = 1,nvpsa ! MCNAP 

                  if (tag_mat(im).eq.tagvpsa(ivpsa)) then ! MCNAP 

				     which_vpsa = ivpsa ! MCNAP 

				  endif ! MCNAP 

			   enddo ! MCNAP



			   which_vpsl = 0 ! MCNAP 

			   do ivpsl = 1,nvpsl ! MCNAP 

                  if (tag_mat(im).eq.tagvpsl(ivpsl)) then ! MCNAP 

				     which_vpsl = ivpsl ! MCNAP 

				  endif ! MCNAP 

			   enddo ! MCNAP



			   which_vpsd = 0 ! MCNAP 

			   do ivpsd = 1,nvpsd ! MCNAP 

                  if (tag_mat(im).eq.tagvpsd(ivpsd)) then ! MCNAP 

				     which_vpsd = ivpsd ! MCNAP 

				  endif ! MCNAP 

			   enddo ! MCNAP



			   which_maps = 0 ! Marco 

			   do imaps = 1,nmaps ! Marco 

                  if (tag_mat(im).eq.tagmaps(imaps)) then ! Marco 

				     which_maps = imaps ! Marco 

				  endif ! Marco 

			   enddo ! Marco

			

			   which_matg = 0 ! Mexico Paco grad_lin 29.11.2004

			   do imatg = 1,nmatg ! Mexico Paco grad_lin 29.11.2004

                  if (tag_mat(im).eq.tagmatg(imatg)) then ! Mexico Paco grad_lin 29.11.2004

				     which_matg = imatg ! Mexico Paco grad_lin 29.11.2004

				  endif ! Mexico Paco grad_lin 29.11.2004

			   enddo ! Mexico Paco grad_lin 29.11.2004



			   ! NOT VERY EFFICIENT CHOICE!!!
			   ! SHOULD BE CHANGED!!!
			   
			   
			   rho = prop_mat(im,1)
			   lambda = prop_mat(im,2)
               mu = prop_mat(im,3)
               
               if (type_mat(im).eq.2) then
                  do k = 1,nelem_ebe
                     ie = ielem_ebe(k)
                     
                     if (cs(cs(ie -1) +0).eq.tag_mat(im)) then
                        do i = 1,nn
                           dxdy_el(i) = beta1(ie) + gamma1(ie) * ct(i)
                           dydy_el(i) = beta2(ie) + gamma2(ie) * ct(i)
                        enddo
                        
                        do j = 1,nn
                           dxdx_el(j) = alfa1(ie) + gamma1(ie) * ct(j)
                           dydx_el(j) = alfa2(ie) + gamma2(ie) * ct(j)
                        enddo
                        
                        do j = 1,nn
                           do i = 1,nn
                              is = nn*(j -1) +i
                              in = cs(cs(ie -1) + is)

							  !if (node_ind_ebe(in).ne.0) then
							  !
                              !   id1 = node_index(in)
                              !   iaz = update_index_el_az(id1 -1)
							  !	 vivel_az(iaz) = 0.0d0 ! Clara 
							  !	 sxx_az(iaz) = 0.0d0 ! Clara 
							  !	 syy_az(iaz) = 0.0d0 ! Clara 
							  !	 sxy_az(iaz) = 0.0d0 ! Clara 
                              !
							  !endif
                              
                              iaz = global_index_el_az(in)
                              ux_el(i,j) = u1_az(iaz)

							  !vivel_el(i,j,ie) = vivel_az(iaz) 
							  !nodal_counter(in) = nodal_counter(in)+1
                              
                              iaz = global_index_el_az(in +nnt)
                              uy_el(i,j) = u1_az(iaz)


                           enddo
                        enddo

                        call make_internal_force(lambda,mu,rho,
                            nmatg,valmatg,which_matg,tagmatg,typematg,& ! Mexico Paco grad_lin 29.11.2004
                            nn,ct,ww,dd,&
                            dxdx_el,dxdy_el,dydx_el,dydy_el,&
                            ux_el,uy_el,&
                            duxdx_el,duxdy_el,duydx_el,duydy_el,&
                            sxx_el,syy_el,sxy_el,fx_el,fy_el,&
                            check_node_sism,check_dist_node_sism,&
                            length_cns,ie,facsmom,&
                            nl_sism,&
                            func_type,func_indx,func_data,nfunc_data,nf,tt2,&

                            nvpcl,fvpcl,valvpcl,which_vpcl,& ! Clara
                            nvpsa,valvpsa,which_vpsa,& ! MCNAP
                            nvpsl,valvpsl,which_vpsl,& ! MCNAP
                            nvpsl,valvpsd,which_vpsd,& ! MCNAP
                            nmaps,fmaps,valmaps,which_maps,& ! Marco 
                            cs_nnz,cs,vivel_length,vivel_el,dt,nelem_ebe,& ! Clara 
                            vistn1_0,vistn2_0,vistn3_0,& ! Clara 
                            tagvpcl,tagvpsa,tagvpsl,tagvpsd,tag_mat(im),& ! Clara &  MCNAP
                            tagmaps,& ! Marco
                            xs,ys,nnt,in,&
                            sxx_ps,syy_ps,sxy_ps,szz_ps,its,& ! Clara   
                            a1n1,a1n2,a1n3,a1n4,a1n5,a1n6,&  !MCNAP vpsa
                            epn1,epn2,epn3,epn4,epn5,epn6,&  !MCNAP vpsa
                            z14n,u7n,densit)                 !MCNAP vpsa

                        do j = 1,nn
                           do i = 1,nn
                              is = nn*(j -1) +i
                              in = cs(cs(ie -1) + is)
                              
                              if (node_ind_ebe(in).ne.0) then
                                 id1 = node_index(in)
                                 id2 = node_index(in) +nnd
                                 
                                 iaz = update_index_el_az(id1 -1)
                                 fk_az(iaz) = fk_az(iaz) +fx_el(i,j)
								 sxx_az(iaz) = sxx_az(iaz) + sxx_el(i,j)! Clara 
								 syy_az(iaz) = syy_az(iaz) + syy_el(i,j)! Clara 
								 sxy_az(iaz) = sxy_az(iaz) + sxy_el(i,j)! Clara
								 vivel_az(iaz) = vivel_az(iaz) + vivel_el(i,j,ie)
                                 vistn1_0_az(iaz) = vistn1_0_az(iaz) + vistn1_0(i,j,ie) ! Clara
								 vistn2_0_az(iaz) = vistn2_0_az(iaz) + vistn2_0(i,j,ie) ! Clara
								 vistn3_0_az(iaz) = vistn3_0_az(iaz) + vistn3_0(i,j,ie) ! Clara

                                 iaz = update_index_el_az(id2 -1)
                                 fk_az(iaz) = fk_az(iaz) +fy_el(i,j)

                              endif
                           enddo
                        enddo
                     endif


                  enddo
               endif
                     


               
            enddo
         endif

		 do in = 1,nnt ! Clara 
			id1 = node_index(in) ! Clara 
            iaz = update_index_el_az(id1 -1) ! Clara 
			vivel_az(iaz) = vivel_az(iaz) / nodal_counter(in) ! Clara
			vistn1_0_az(iaz) = vistn1_0_az(iaz) / nodal_counter(in) ! Clara
			vistn2_0_az(iaz) = vistn2_0_az(iaz) / nodal_counter(in) ! Clara
			vistn3_0_az(iaz) = vistn3_0_az(iaz) / nodal_counter(in) ! Clara

			
			!CAPIRE PERCHE' la variabile "vivel" va normalizzata, mentre gli stress no!!! (???)
			!CHIEDERE LUCA
			 
			sxx_az(iaz) = sxx_az(iaz) / nodal_counter(in) ! Clara 
			syy_az(iaz) = syy_az(iaz) / nodal_counter(in) ! Clara 
			sxy_az(iaz) = sxy_az(iaz) / nodal_counter(in) ! Clara 
		enddo ! Clara 
         

         
         if (nelem_abc.gt.0) then
            do im = 1,nm
               rho = prop_mat(im,1)
               lambda = prop_mat(im,2)
               mu = prop_mat(im,3)
               
               if (type_mat(im).eq.2) then
                  do k = 1,nelem_abc
                     ie = ielem_abc(k)
                     iedge = iedge_abc(k)
                     
                     if (cs(cs(ie -1) +0).eq.tag_mat(im)) then
                        
                        ied1 = cs_bc(cs_bc(iedge -1) +1)
                        ied2 = cs_bc(cs_bc(iedge) -1)
                        
                        iel1 = cs(cs(ie -1) +1)
                        iel2 = cs(cs(ie -1) +nn)
                        iel3 = cs(cs(ie -1) +nn*nn)
                        iel4 = cs(cs(ie -1) +nn*(nn -1) +1)
                        
                        
! First edge
                        if (((ied1.eq.iel1).and.(ied2.eq.iel2))&
                             .or.((ied2.eq.iel1).and.(ied1.eq.iel2))) then
                           edge_lx = xs(iel2) - xs(iel1)
                           edge_ly = ys(iel2) - ys(iel1)
                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)
                           edge_nx = edge_ly / edge_ll
                           edge_ny = -1.0d0 * edge_lx / edge_ll
                           edge_ia = 1; edge_ib = nn
                           edge_ja = 1; edge_jb = 1
                        endif
                        
! Second edge
                        if (((ied1.eq.iel2).and.(ied2.eq.iel3))&
                             .or.((ied2.eq.iel2).and.(ied1.eq.iel3))) then
                           edge_lx = xs(iel3) - xs(iel2)
                           edge_ly = ys(iel3) - ys(iel2)
                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)
                           edge_nx = edge_ly / edge_ll
                           edge_ny = -1.0d0 * edge_lx / edge_ll
                           edge_ia = nn; edge_ib = nn
                           edge_ja = 1;  edge_jb = nn
                        endif
                        
! Third edge
                        if (((ied1.eq.iel3).and.(ied2.eq.iel4))&
                             .or.((ied2.eq.iel3).and.(ied1.eq.iel4))) then
                           edge_lx = xs(iel4) - xs(iel3)
                           edge_ly = ys(iel4) - ys(iel3)
                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)
                           edge_nx = edge_ly / edge_ll; 
                           edge_ny = -1.0d0 * edge_lx / edge_ll
                           edge_ia = 1;  edge_ib = nn
                           edge_ja = nn; edge_jb = nn
                        endif
                        
! Fourth edge
                        if (((ied1.eq.iel4).and.(ied2.eq.iel1))&
                             .or.((ied2.eq.iel4).and.(ied1.eq.iel1))) then
                           edge_lx = xs(iel1) - xs(iel4)
                           edge_ly = ys(iel1) - ys(iel4)
                           edge_ll = dsqrt(edge_lx*edge_lx + edge_ly*edge_ly)
                           edge_nx = edge_ly / edge_ll
                           edge_ny = -1.0d0 * edge_lx / edge_ll
                           edge_ia = 1; edge_ib = 1
                           edge_ja = 1; edge_jb = nn
                        endif
                        
                        
                        do i = 1,nn
                           dxdy_el(i) = beta1(ie) + gamma1(ie) * ct(i)
                           dydy_el(i) = beta2(ie) + gamma2(ie) * ct(i)
                        enddo
                        
                        do j = 1,nn
                           dxdx_el(j) = alfa1(ie) + gamma1(ie) * ct(j)
                           dydx_el(j) = alfa2(ie) + gamma2(ie) * ct(j)
                        enddo
                        
                        do j = 1,nn
                           do i = 1,nn
                              is = nn*(j -1) +i
                              in = cs(cs(ie -1) + is)
                              
                              iaz = global_index_el_az(in)
                              ux_el(i,j) = u1_az(iaz)
                              vx_el(i,j) = (3.0d0*u1_az(iaz) -4.0d0*u0_az(iaz)&
                                           +1.0d0*u_1_az(iaz)) / (2.0d0*dt) 
                              
                              iaz = global_index_el_az(in +nnt)
                              uy_el(i,j) = u1_az(iaz)
                              vy_el(i,j) = (3.0d0*u1_az(iaz) -4.0d0*u0_az(iaz)&
                                           +1.0d0*u_1_az(iaz)) / (2.0d0*dt) 
                           enddo
                        enddo
                        
                        call make_absorbing_force(lambda,mu,rho,nn,ct,ww,dd,&
                                  dxdx_el,dxdy_el,dydx_el,dydy_el,&
                                  edge_nx,edge_ny,edge_ll,&
                                  edge_ia,edge_ja,edge_ib,edge_jb,&
                                  ux_el,uy_el,vx_el,vy_el,fx_el,fy_el)
                        
                        do j = 1,nn
                           do i = 1,nn
                              is = nn*(j -1) +i
                              in = cs(cs(ie -1) + is)
                              
                              if (node_ind_ebe(in).ne.0) then
                                 id1 = node_index(in)
                                 id2 = node_index(in) +nnd
                                 
                                 iaz = update_index_el_az(id1 -1)
                                 fk_az(iaz) = fk_az(iaz) +fx_el(i,j)
                                 
                                 iaz = update_index_el_az(id2 -1)
                                 fk_az(iaz) = fk_az(iaz) +fy_el(i,j)
                              endif
                           enddo
                        enddo
                     endif
                     
                  enddo
               endif
               
            enddo
         endif
        
        
        if (make_damping_yes_or_not.eq.0) then 
           do iaz = 0,N_update_el_az -1
              u2_az(iaz) = 2.0d0 * u1_az(iaz) - u0_az(iaz) &
                           + (fe_az(iaz) - fk_az(iaz)) / mvec(iaz)
           enddo
        else 
           do iaz = 0,N_update_el_az -1
              u2_az(iaz) = ( fe_az(iaz) - fk_az(iaz) - kcvec(iaz) * u1_az(iaz) &
                           + cvec(iaz) * u0_az(iaz) + mvec(iaz) &
                           * ( 2.0d0 * u1_az(iaz) - u0_az(iaz) ) ) &
                           / new(iaz)
           enddo
         endif

         
         do fn = 1,nf
            func_value(fn) = get_func_value(nf,func_type,func_indx,func_data,nfunc_data,&
                                            fn,tt2,0.0d0)
         enddo
         
         
!        Dirichlet b.c.
         if (nnode_dirX.gt.0) then
            do i = 1,nnode_dirX
               in = inode_dirX(i)
               
               if (node_index(in).ne.0) then
                  id1 = node_index(in)
                  iaz = update_index_el_az(id1 -1)
                  
                  !u2_az(iaz) = 0.0d0

				  u2_az(iaz) = u2_star(iaz) ! M Mexico _A_
                  do fn = 1,nf
                     u2_az(iaz) = u2_az(iaz) + Fmat(fn,id1) * func_value(fn)
                  enddo
               endif
            enddo
         endif



!        Dirichlet Point b.c. ! Marco

         if (nl_dipX.gt.0) then  ! Marco

            do i = 1,nl_dipX ! Marco

               in = in_dipX(i) ! Marco

               

               if (node_index(in).ne.0) then ! Marco

                  id1 = node_index(in) ! Marco

                  iaz = update_index_el_az(id1 -1) ! Marco

                  

                  !u2_az(iaz) = 0.0d0

				  u2_az(iaz) = u2_star(iaz) ! M Mexico _A_

                  do fn = 1,nf ! Marco

                     u2_az(iaz) = u2_az(iaz) + Fmat(fn,id1) * func_value(fn) ! Marco

                  enddo ! Marco

               endif ! Marco

            enddo ! Marco

         endif ! Marco
         
         if (nnode_dirY.gt.0) then
            do i = 1,nnode_dirY
               in = inode_dirY(i)
               
               if (node_index(in).ne.0) then
                  id2 = node_index(in) +nnd
                  iaz = update_index_el_az(id2 -1)


                  !u2_az(iaz) = 0.0d0

				  u2_az(iaz) = u2_star(iaz) ! M Mexico _A_
                  do fn = 1,nf
                     u2_az(iaz) = u2_az(iaz) + Fmat(fn,id2) * func_value(fn)
                  enddo
               endif
            enddo
         endif



		 if (nl_dipY.gt.0) then ! Marco

            do i = 1,nl_dipY ! Marco

               in = in_dipY(i) ! Marco

               

               if (node_index(in).ne.0) then ! Marco

                  id2 = node_index(in) +nnd ! Marco

                  iaz = update_index_el_az(id2 -1) ! Marco

                  

                  !u2_az(iaz) = 0.0d0

				  u2_az(iaz) = u2_star(iaz) ! M Mexico _A_

                  do fn = 1,nf ! Marco

                     u2_az(iaz) = u2_az(iaz) + Fmat(fn,id2) * func_value(fn) ! Marco

                  enddo ! Marco

               endif ! Marco

            enddo ! Marco

         endif ! Marco


!        End Dirichlet b.c.
         
         
         call system_clock(COUNT=clock_finish)
         
         time_in_seconds = float(clock_finish - clock_start) &
                         / float(clock(2))
         
         time_total = time_total + time_in_seconds
         
         write(40,*)'Time-step ',(its-1),' done; CPU-time = ',         &
              time_in_seconds,'s '
         
         if (nmonit.ge.1) then
		    if (int(real(its)/ndt_monitor).eq.(real(its)/ndt_monitor)) then 
            do i = 1,nmonit
               unit_monitor = 40 + i
               in = node_m(i)
               if (node_index(in).ne.0) then
                  id1 = node_index(in)
                  id2 = node_index(in) +nnd
                  iaz = update_index_el_az(id1 -1)
                  jaz = update_index_el_az(id2 -1)
                  
                  write(unit_monitor,'(3E16.8)') &
                       tt1,u1_az(iaz),u1_az(jaz)
               endif
            enddo


            ! DA SISTEMARE!!!
			do i = 1,nmonit ! Clara 
               unit_monitor = 10000 + i ! Clara 
               in = node_m(i) ! Clara 
               if (node_index(in).ne.0) then ! Clara 
                  id1 = node_index(in) ! Clara 
                  iaz = update_index_el_az(id1 -1) ! Clara 
                  
                  write(unit_monitor,'(4E16.8)') & ! Clara 
                       tt1,sxx_az(iaz),syy_az(iaz),sxy_az(iaz) ! Clara 
               endif ! Clara 
            enddo ! Clara 

			endif
         endif





		 !do iaz = 0,N_update_el_az -1



		!	vivel_az_1(iaz) = vivel_az(iaz)! M Mexico _A_

	  	!	vistn1_0_az_1(iaz) = vistn1_0_az(iaz)! M Mexico _A_

	  	!	vistn2_0_az_1(iaz) = vistn2_0_az(iaz)! M Mexico _A_

	  	!	vistn3_0_az_1(iaz) = vistn3_0_az(iaz)! M Mexico _A_ 



		!enddo ! M Mexico _A_
         
         
! *********************************************************         
         
         if ((nsnap.gt.0).and.(isnap.le.nsnap)) then
            if (its.ge.itersnap(isnap)) then
               call write_bin_file(file_outE,isnap,myid,&
                                   N_update_el_az,out_index_el_az,u1_az)
			   
			   call write_bin_file(file_outV,isnap,myid,& ! Clara 
                                 N_update_el_az,out_index_el_az,vivel_az) ! Clara 
			   

				call write_bin_file(file_outT1,isnap,myid,& ! Clara 
                                N_update_el_az,out_index_el_az,vistn1_0_az) ! Clara

				call write_bin_file(file_outT2,isnap,myid,& ! Clara 
                                N_update_el_az,out_index_el_az,vistn2_0_az) ! Clara

				call write_bin_file(file_outT3,isnap,myid,& ! Clara 
                                N_update_el_az,out_index_el_az,vistn3_0_az) ! Clar
			   
               isnap = isnap + 1
 
               if (myid.eq.0) write(*,'(A,E12.4)')'TIME = ',tt1
            endif
         endif



	  if ((nsnapps.gt.0).and.(isnapps.le.nsnapps)) then

		if (its.ge.itersnapps(isnapps)) then



			call write_bin_file(file_outS1,isnapps,myid,& ! Clara 

                                N_update_el_az,out_index_el_az,sxx_az) ! Clara

		    call write_bin_file(file_outS2,isnapps,myid,& ! Clara 

                                N_update_el_az,out_index_el_az,syy_az) ! Clara

            call write_bin_file(file_outS3,isnapps,myid,& ! Clara 

                                N_update_el_az,out_index_el_az,sxy_az) ! Clara





			call write_bin_file_new(file_out_u1,isnapps,myid,& ! M Mexico _A_ 

                                N_update_el_az,out_index_el_az,u1_az) ! M Mexico _A_

			call write_bin_file_new(file_out_u0,isnapps,myid,& ! M Mexico _A_ 

                                N_update_el_az,out_index_el_az,u0_az) ! M Mexico _A_

			call write_bin_file_new(file_out_Fe,isnapps,myid,& ! M Mexico _A_ 

                                N_update_el_az,out_index_el_az,fe_az) ! M Mexico _A_

			call write_bin_file_new(file_out_Fk,isnapps,myid,& ! M Mexico _A_ 

                                N_update_el_az,out_index_el_az,fk_az) ! M Mexico _A_

			call write_bin_file_new(file_out_vivel,isnapps,myid,& ! M Mexico _A_ 

                                N_update_el_az,out_index_el_az,vivel_az_1) ! M Mexico _A_

			call write_bin_file_new(file_out_vistn1_0,isnapps,myid,& ! M Mexico _A_

                                N_update_el_az,out_index_el_az,vistn1_0_az_1) ! M Mexico _A_

			call write_bin_file_new(file_out_vistn2_0,isnapps,myid,& ! M Mexico _A_ 

                                N_update_el_az,out_index_el_az,vistn2_0_az_1) ! M Mexico _A_

			call write_bin_file_new(file_out_vistn3_0,isnapps,myid,& ! M Mexico _A_ 

                                N_update_el_az,out_index_el_az,vistn3_0_az_1) ! M Mexico _A_





			isnapps = isnapps + 1



		endif

	  endif
         
! *********************************************************
         
         
         
         do i = 0,N_update_el_az -1
            u_1_az(i) = u0_az(i)
            u0_az(i) = u1_az(i)
            u1_az(i) = u2_az(i)
         enddo
         
      enddo
      
      
      close(40)
      
      if (nmonit.ge.1) then 
         do i = 1,nmonit
            unit_monitor = 40 + i
            in = node_m(i)
            if (node_index(in).ne.0) then
               close(unit_monitor)
            endif
         enddo

		 do i = 1,nmonit
            unit_monitor = 10000 + i
            in = node_m(i)
            if (node_index(in).ne.0) then
               close(unit_monitor)
            endif
         enddo
      endif
      
      if (myid.eq.0) write(*,'(A,F8.4,A)')'Mean time-step CPU time= ', &
           time_total / dfloat(nts - 1),' s'
      
      
      deallocate(update_el_az,update_index_el_az)
      deallocate(global_index_el_az)
      
      deallocate(u_1_az,u0_az,u1_az,u2_az,fk_az,fe_az)
      deallocate(out_index_el_az)
      
      deallocate(func_value)
      
      if (nnode_dirX.gt.0) deallocate(inode_dirX)
      if (nnode_dirY.gt.0) deallocate(inode_dirY)
      if (nnode_neuX.gt.0) deallocate(inode_neuX)
      if (nnode_neuY.gt.0) deallocate(inode_neuY)

	  if (nnode_neuN.gt.0) deallocate(inode_neuN) ! M Mexico Arroyo

      if (nnode_neuT.gt.0) deallocate(inode_neuT) ! M Mexico Arroyo
      if (nedge_abc.gt.0) deallocate(iedge_abc)
      if (nelem_abc.gt.0) deallocate(ielem_abc)
      if (nnode_ebe.gt.0) deallocate(inode_ebe)
      if (nelem_ebe.gt.0) deallocate(ielem_ebe)
      
	  ! SISTEMARE I DEALLOCATE SULLE VARIABILI INTRODOTTE DA CLARA
      
      end subroutine time_loop_el
      
      
      
!     ***********************************************************************
      
      subroutine count_edge_nodes(nnode,nnz,cs,nedge_nodes)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      implicit none
      
      integer*4 :: nnode,nnz,nedge_nodes
      integer*4, dimension(0:nnz) :: cs
      
      integer*4, dimension(:),allocatable :: node_index
      integer*4 :: i,ie,ne,nn
      
      allocate (node_index(nnode))
      
      do i = 1,nnode
         node_index(i) = 0
      enddo
      
      nedge_nodes = 0
      
      ne = cs(0) -1
      
      if (nnz.gt.0) then
         do ie = 1,ne
            nn = cs(ie) - cs(ie -1) -1
            do i = 1,nn
               node_index(cs(cs(ie -1) +i)) = 1
            enddo
         enddo
         
         do i = 1,nnode
            if (node_index(i).ne.0) then
               nedge_nodes = nedge_nodes +1
               node_index(i) = nedge_nodes
            endif
         enddo
      endif
      
      deallocate (node_index)
      
      return
      
      end subroutine count_edge_nodes
      
      
      
!     ***********************************************************************
      
      subroutine get_edge_nodes(nnode,nnz,cs,nl_bc,tag_bc,&
                                nedge_nodes,node_index)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      
      implicit none
      
      integer*4 :: nnode,nnz,nl_bc
      integer*4, dimension(0:nnz) :: cs
      integer*4, dimension(nl_bc) :: tag_bc
      integer*4, dimension(nnode) :: node_index
      integer*4 :: nedge_nodes
      
      integer*4 :: i,j,ie,ne,nn,check
      
      do i = 1,nnode
         node_index(i) = 0
      enddo
      
      nedge_nodes = 0
      
      if (nnz.gt.0) then
         ne = cs(0) -1
         do ie = 1,ne
            nn = cs(ie) - cs(ie -1) -1
            
            check = 0
            do j = 1,nl_bc
               if (cs(cs(ie -1) +0).eq.tag_bc(j)) check = 1
            enddo
            
            if (check.ne.0) then
               do i = 1,nn
                  node_index(cs(cs(ie -1) +i)) = 1
               enddo
            endif
         enddo
         
         do i = 1,nnode
            if (node_index(i).ne.0) then
               nedge_nodes = nedge_nodes +1
               node_index(i) = nedge_nodes
            endif
         enddo
      endif
      
      return
      
      end subroutine get_edge_nodes
      
      
      
!     ****************************************************************
      
      subroutine get_edge_element(Ennz,Ebin,n1,n2,ie)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      implicit none
      
      integer*4 :: Ennz,n1,n2,ie
      integer*4, dimension(0:Ennz) :: Ebin
      
      integer*4 :: i,j
      
      ie = 0
      do i = Ebin(n1 -1),Ebin(n1) -1
         do j = Ebin(n2 -1),Ebin(n2) -1
            if (Ebin(i).eq.Ebin(j)) then
               ie = Ebin(i)
            endif
         enddo
      enddo
      
      return
      
      end subroutine get_edge_element
      
      
      
! *********************************************************
      
      real*8 function get_func_value(nf,func_type,func_indx,func_data,nfunc_data,n,t,&
                                    t0_sism)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     Modified: Marco Stupazzini (introduced new functions since 2002)
!     
      implicit none
      
      integer*4 :: nf,n,nfunc_data
      integer*4, dimension(nf) :: func_type
      integer*4, dimension(nf +1) :: func_indx
      real*8, dimension(nfunc_data) :: func_data
      real*8 :: t, beta2, t0_sism
      
      integer*4 :: i
      real*8 :: val,PI,tt_t0,t0,t1,v0,v1
	  real*8 :: mean,sigma,amp
      real*8 :: tp,ts
      
      val = 0.0d0
      
      if (func_type(n).eq.0) then
         val = 1.0d0

      ! 1 - Ricker "beta" type
      elseif (func_type(n).eq.1) then
         tt_t0 = t - func_data(func_indx(n) +1)
         val = (1.0d0 - 2.0d0*func_data(func_indx(n))*tt_t0*tt_t0) &
              * dexp(-1.0d0*func_data(func_indx(n))*tt_t0*tt_t0)
      
      ! 2 - Ricker "cos" type
      elseif (func_type(n).eq.2) then
         PI = 4.0d0 * datan(1.0d0)
         tt_t0 = t - func_data(func_indx(n) +1)
         val = dcos(PI*func_data(func_indx(n))*tt_t0) &
              * dexp(-0.5d0*func_data(func_indx(n)) &
              * func_data(func_indx(n))*tt_t0*tt_t0)

      ! 3 - Force history
      elseif (func_type(n).eq.3) then
         do i = func_indx(n),func_indx(n+1) -3,2
            t0 = func_data(i)
            t1 = func_data(i +2)
            v0 = func_data(i +1)
            v1 = func_data(i +3)
            if ((t.ge.t0).and.(t.le.t1)) then
               val = (v1 - v0) / (t1 - t0) * (t - t0)  + v0
            endif
         enddo
      
      ! 4 - First derivative of the Ricker ( d(Ricker(t))/dt )
      elseif (func_type(n).eq.4) then
         tt_t0 = t - func_data(func_indx(n) +1)
         beta2=func_data(func_indx(n))
         val = 2.0d0*beta2*tt_t0 &
              *(-3.0d0 + 2.0d0*beta2*tt_t0*tt_t0) &
              * dexp(-beta2*tt_t0*tt_t0)

      ! 5 - Ricker "beta" type for seismic moment
      elseif (func_type(n).eq.5) then
         tt_t0 = t - func_data(func_indx(n) +1) - t0_sism
         val = (1.0d0 - 2.0d0*func_data(func_indx(n))*tt_t0*tt_t0) &
              * dexp(-1.0d0*func_data(func_indx(n))*tt_t0*tt_t0)
         !write(55,*),func_type(n),func_data(func_indx(n)),&
         !            func_data(func_indx(n)+1),val
        

      ! - Double impulse for seismic moment.. not yet implemented
      !elseif (func_type(n).eq.6) then
      !tt_t0 = t - func_data(func_indx(n) +1)
      !   beta2=func_data(func_indx(n))
      !   val = tt_t0*dexp(-2.0d0*beta2*tt_t0*tt_t0)

      
	  ! 6 - Ramp function:
	  !
	  !
	  !   /\ F(t)
	  !   |
	  ! S1|...............  ______________
	  !   |               /.
      !   |             /  .
	  !   |           /    .
	  ! S0|________ /      .
	  !   |         .      .
	  !   +-------------------------------> Time
	  !             t0     t1
	  !
	  !
	  !   if  t0 > t then F(t) = S0
	  !
	  !   if  t0 <= t <= t1 then 
	  !                          _               _
	  !             S1 - S0     |      S1 - S0    |
	  !      F(t) = ------- t + | S0 - ------- t0 |
	  !             t1 - t0     |_     t1 - t0   _|
	  !
	  !   if  t1 < t then F(t) = S1
	  !

      elseif (func_type(n).eq.6) then

         tt_t0 = t - func_data(func_indx(n) +2)

		 if (func_data(func_indx(n) + 2).gt.t) then
			val = func_data(func_indx(n))
	     elseif ((func_data(func_indx(n) + 2).le.t).and.(func_data(func_indx(n) + 3).ge.t)) then
			
			val = ( func_data(func_indx(n) + 1) - func_data(func_indx(n)) ) / &
		          ( func_data(func_indx(n) + 3) - func_data(func_indx(n) + 2) )	* &
				  tt_t0 + func_data(func_indx(n))
          
         elseif (func_data(func_indx(n) + 3).lt.t) then
		    val = func_data(func_indx(n) + 1)
	     endif

      ! 7 - Linear function:
	  !
	  !
	  !   /\ F(t)
	  !   |
	  ! S1|.........________
	  !   |        |        |
      !   |        |        | 
	  !   |        |        |
	  ! S0|________|        |______________
	  !   |         .      .
	  !   +-------------------------------> Time
	  !             t0     t1
	  !
	  !
	  !   if  t0 > t then F(t) = S0
	  !
	  !   if  t0 <= t <= t1 then F(t) = S1
	  !
	  !   if  t1 < t then F(t) = S0
	  !

      elseif (func_type(n).eq.7) then

         tt_t0 = t - func_data(func_indx(n) +2)

		 if (func_data(func_indx(n) + 3).lt.t) then

			val = func_data(func_indx(n))

	     elseif ((func_data(func_indx(n) + 2).le.t).and.(func_data(func_indx(n) + 3).ge.t)) then
			
			val = func_data(func_indx(n) + 1) 
          
         elseif (func_data(func_indx(n) + 2).gt.t) then

		    val = func_data(func_indx(n))

	     endif
	  

      ! 8 - sin(w*t): PLANE WAVE (what is implemented is the derivative of sin(wt) )
	  elseif (func_type(n).eq.8) then

		tt_t0 = t - func_data(func_indx(n) +1)

		if ((1/func_data(func_indx(n))).lt.tt_t0) then

			val = 0.0d0

	    elseif (0.le.tt_t0) then

			PI = 4.0d0 * datan(1.0d0)
			val = cos(2*PI*func_data(func_indx(n))*tt_t0)*(2*PI*func_data(func_indx(n)))

		else
			
			val = 0.0d0
		
		endif

      ! 9 - Gaussian distribution (parameters: mean, sigma, amp)
      elseif (func_type(n).eq.9) then

		mean = func_data(func_indx(n))
		sigma = func_data(func_indx(n) +1)
		amp = func_data(func_indx(n) +2)

		PI = 4.0d0 * datan(1.0d0)

        val = amp * 1/((2*PI*sigma**2)**.5) * exp(-(t-mean)**2/(2*sigma**2))



	  ! 10 - sin(w1*t) * sin(w2*t)  -> sinus func. modulated by a sencond sinus func. with
	  !                                characteristic period T2 = 1/(2*n*T1), where T1 = 1/f
	  !      w1 = 2*pi*f
	  !      w2 = 2*pi*(f/(2*n))
      ! 
      ! parameters:
	  ! * t0 = func_data(func_indx(n) +1)
	  ! * f  = func_data(func_indx(n))
	  ! * n  = func_data(func_indx(n) +2)
	  !
      !  /\  sin(w1*t)
      !  |                                                                                                
      !  |         ***                     ***                     ***                     ***                
      !  |       **   **                 **   **                 **   **                 **   **              
      !  |      *       *               *       *               *       *               *       *             
      !  |     *         *             *         *             *         *             *         *            
      !  o----*-----------*-----------*-----------*-----------*-----------*-----------*-----------*-----------*-----------> Time
      !       |            *         *|            *         *             *         *             *         *|
      !       |             *       * |             *       *               *       *               *       * |
      !       |              **   **  |              **   **                 **   **                 **   **  |
      !       |                ***    |                ***                     ***                     ***    |
      !       |_______________________|                                                                       |
      !       |       T1 = 1/f        |                                                                       |
      !       t0                                                                                             n*(1/f)+t0
	  !
	  !

	  elseif (func_type(n).eq.10) then

		tt_t0 = t - func_data(func_indx(n) +1)

		if (tt_t0.gt.(func_data(func_indx(n) +2)*(1/func_data(func_indx(n))))) then

			val = 0.0d0

	    elseif (tt_t0.ge.0) then

			PI = 4.0d0 * datan(1.0d0)
			val = sin(2*PI*func_data(func_indx(n))*tt_t0) &
			      * sin(2*PI*func_data(func_indx(n))/(2*func_data(func_indx(n) +2))*tt_t0)
			      
		else
			
			val = 0.0d0
		
		endif

	  ! 11 - w1*cos(w1*t) * sin(w2*t) +
	  !      sin(w1*t) * w2 * cos(w2*t)    -> sinus func. modulated by a sencond sinus func. with
	  !                                     characteristic period T2 = 1/(2*n*T1), where T1 = 1/f
	  !      w1 = 2*pi*f
	  !      w2 = 2*pi*(f/(2*n))
	  !
	  ! PLANE WAVE
	  ! sin(w1*t) * sin(w2*t)
      ! 
      ! parameters:
	  ! * t0 = func_data(func_indx(n) +1)
	  ! * f  = func_data(func_indx(n))
	  ! * n  = func_data(func_indx(n) +2)
	  !
      !  /\  sin(w1*t)
      !  |                                                                                                
      !  |         ***                     ***                     ***                     ***                
      !  |       **   **                 **   **                 **   **                 **   **              
      !  |      *       *               *       *               *       *               *       *             
      !  |     *         *             *         *             *         *             *         *            
      !  o----*-----------*-----------*-----------*-----------*-----------*-----------*-----------*-----------*-----------> Time
      !       |            *         *|            *         *             *         *             *         *|
      !       |             *       * |             *       *               *       *               *       * |
      !       |              **   **  |              **   **                 **   **                 **   **  |
      !       |                ***    |                ***                     ***                     ***    |
      !       |_______________________|                                                                       |
      !       |       T1 = 1/f        |                                                                       |
      !       t0                                                                                             n*(1/f)+t0
	  !
	  !

	  elseif (func_type(n).eq.11) then

		tt_t0 = t - func_data(func_indx(n) +1)

		if (tt_t0.gt.(func_data(func_indx(n) +2)*(1/func_data(func_indx(n))))) then

			val = 0.0d0

	    elseif (tt_t0.ge.0) then

			PI = 4.0d0 * datan(1.0d0)
			val = 2*PI*func_data(func_indx(n)) * cos(2*PI*func_data(func_indx(n))*tt_t0) &
			      * sin(2*PI*func_data(func_indx(n))/(2*func_data(func_indx(n) +2))*tt_t0) &
				  + sin(2*PI*func_data(func_indx(n))*tt_t0) &
			      * 2*PI*func_data(func_indx(n))/(2*func_data(func_indx(n) +2)) &
				  * cos(2*PI*func_data(func_indx(n))/(2*func_data(func_indx(n) +2))*tt_t0)
			      
		else
			
			val = 0.0d0
		
		endif


      ! 12 - sigmf(t,[a c]) = amp*(1/(1+exp(-a*(t-c))))
	  ! 
	  !
      !  /\  
      !  |                                                                                                
      !  |............*************************......amp      
      !  |          ** : 
	  !  |     a   *   :
	  !  |        *    :
	  !  |      **     :   
      !  o******----------------------------------> Time
      !       |        |
      !       |____c___|                                                                       |
      !       |        |                                                                       |
      !                                                                                                         n*(1/f)+t0
	  !
	  !

	  elseif (func_type(n).eq.12) then

		tt_t0 = t - func_data(func_indx(n) +2)
        !tt_t0 = t - 1

		val = func_data(func_indx(n)) &
		      * (1/(1+exp(-func_data(func_indx(n) +1)*(tt_t0))))



      ! 14 - sin(w*t) INFINITE PERIOD!!!:

	  elseif (func_type(n).eq.14) then



		tt_t0 = t - func_data(func_indx(n) +1)



		!if ((1/func_data(func_indx(n))).lt.tt_t0) then



        !	val = 0.0d0



	    if (tt_t0.gt.0) then



			PI = 4.0d0 * datan(1.0d0)

			val = sin(2*PI*func_data(func_indx(n))*tt_t0)



		else

			

			val = 0.0d0

		

		endif



      ! 6 - Triangular function in order to test with PACO analytical sol (prg. "Gradelan")

	  !

	  ! func_data(func_indx(n))    = tp

	  ! func_data(func_indx(n) +1) = ts

	  ! func_data(func_indx(n) +2) = AMP

	  !

	  !

	  !   /\ F(t)

	  !   |

	  !AMP|...............  

	  !   |               /.\

      !   |             /  .  \ 

	  !   |           /    .    \

	  !   +----------             ---------> Time

	  !   |          :     :     :

	  !   

	  !         tp-ts/2    tp    tp+ts/2 

	  !

	  !

	  !   if  t > tp+ts/2 then F(t) = 0

	  !

	  !   if  tp < t <= tp+ts/2 then 

	  !                         

	  !             -2 * AMP

	  !      F(t) = -------   (t-tp) + AMP 

	  !                ts 

	  !   if  tp-ts/2 < t <= tp then 

	  !                         

	  !             2 * AMP

	  !      F(t) = -------   (t-tp) + AMP 

	  !                ts 

	  !

	  !

	  !   else  F(t) = 0

	  !



      elseif (func_type(n).eq.20) then



         tt_t0 = t - func_data(func_indx(n))



		 if (tt_t0.gt.func_data(func_indx(n) + 1)/2) then



			val = 0.0d0



	     elseif ((tt_t0.gt.0).and.(tt_t0.le.func_data(func_indx(n) + 1)/2)) then

			

			val = -func_data(func_indx(n)+2) * 2 / func_data(func_indx(n)+1) *&

			      tt_t0 + func_data(func_indx(n)+2)

          

         elseif ((tt_t0.le.0).and.(tt_t0.gt.-func_data(func_indx(n) + 1)/2)) then

			

			val = func_data(func_indx(n)+2) * 2 / func_data(func_indx(n)+1) *&

			      tt_t0 + func_data(func_indx(n)+2)

         else

		    

            val = 0.0d0



	     endif



	  ! 6 - SIN CARDINAL

	  !

	  !  disp = sin(pi * f * tt_t0)/ (pi * f * tt_t0) exp(-beta * tt_t0)

	  !

	  ! func_data(func_indx(n))    = f

	  ! func_data(func_indx(n) +1) = t0

	  ! func_data(func_indx(n) +2) = beta

	  ! func_data(func_indx(n) +3) = AMP

	  !

	  !

	  !



      elseif (func_type(n).eq.21) then



         tt_t0 = t - func_data(func_indx(n)+1)



		 if (tt_t0.gt.0) then



			val = sin (pi * func_data(func_indx(n)) * tt_t0) / &

			      (pi * func_data(func_indx(n)) * tt_t0) * &

				  exp(-func_data(func_indx(n) +2) * tt_t0) * &

				  func_data(func_indx(n) +3)



	     elseif (tt_t0.eq.0) then

			

			val = func_data(func_indx(n)+3)

          

         else

		    

            val = sin (pi * func_data(func_indx(n)) * tt_t0) / &

			      (pi * func_data(func_indx(n)) * tt_t0) * &

				  exp(func_data(func_indx(n) +2) * tt_t0) * &

				  func_data(func_indx(n) +3)

	     endif



	  ! 23 - Triangular function in order to test with PACO analytical sol (prg. "Gradelan")

	  !

	  ! func_data(func_indx(n))    = tp

	  ! func_data(func_indx(n) +1) = ts

	  ! func_data(func_indx(n) +2) = AMP

	  !

	  !

	  !   /\ F(t)

	  !   |

	  !AMP|...............  

	  !   |               /.\

      !   |             /  .  \ 

	  !   |           /    .    \

	  !   +----------             ---------> Time

	  !   |          :     :     :

	  !   

	  !         tp-ts/2    tp    tp+ts/2 

	  !

	  !

	  !   if  t > tp+ts/2 then F(t) = 0

	  !

	  !   if  tp < t <= tp+ts/2 then 

	  !                         

	  !             -2 * AMP

	  !      F(t) = -------   (t-tp) + AMP 

	  !                ts 

	  !   if  tp-ts/2 < t <= tp then 

	  !                         

	  !             2 * AMP

	  !      F(t) = -------   (t-tp) + AMP 

	  !                ts 

	  !

	  !

	  !   else  F(t) = 0

	  !



      elseif (func_type(n).eq.22) then



	     tp = 0.3d0

         ts = 0.4d0

		 AMP = 1.0d0



         tt_t0 = t - tp



		 if (tt_t0.gt.ts/2) then



			val = 0.0d0



	     elseif ((tt_t0.gt.0).and.(tt_t0.le.ts/2)) then

			

			val = -AMP * 2 / ts *&

			      tt_t0 + AMP

          

         elseif ((tt_t0.le.0).and.(tt_t0.gt.-ts/2)) then

			

			val = AMP * 2 / ts *&

			      tt_t0 + AMP

         else

		    

            val = 0.0d0



	     endif



	  ! 23 - sin(w*tt_t0) * exp(-beta*tt_t0)

	  ! PLANE WAVE (what is implemented is the derivative of sin(wt) )

	  !

	  elseif (func_type(n).eq.23) then



		tt_t0 = t - func_data(func_indx(n) +1)



	    if (0.le.tt_t0) then



			PI = 4.0d0 * datan(1.0d0)

			!val = sin(2*PI*func_data(func_indx(n))*tt_t0) * &

			!      exp(-func_data(func_indx(n)+2) * tt_t0)



		    val = cos(2*PI*func_data(func_indx(n))*tt_t0) * &

			     (2*PI*func_data(func_indx(n)) ) * &

				 exp(-func_data(func_indx(n)+2) * tt_t0) - &

				 sin(2*PI*func_data(func_indx(n))*tt_t0) * &

			     (func_data(func_indx(n)+2) ) * &

				 exp(-func_data(func_indx(n)+2) * tt_t0)

			      



		else

			

			PI = 4.0d0 * datan(1.0d0)

			val = cos(2*PI*func_data(func_indx(n))*tt_t0) * &

			     (2*PI*func_data(func_indx(n)) ) * &

				 exp(func_data(func_indx(n)+2) * tt_t0) + &

				 sin(2*PI*func_data(func_indx(n))*tt_t0) * &

			     (func_data(func_indx(n)+2) ) * &

				 exp(func_data(func_indx(n)+2) * tt_t0)

		

		endif



	  ! 24 - Morlet wavelt (http://www.eso.org/projects/esomidas/doc/user/98NOV/volb/node312.html) 

	  ! PLANE WAVE 

	  !

	  elseif (func_type(n).eq.24) then



		tt_t0 = t - func_data(func_indx(n) +1)





			PI = 4.0d0 * datan(1.0d0)



		    val = 1/(2*PI)**.5 * cos(2*PI*func_data(func_indx(n))*tt_t0) * &

				 exp((-tt_t0**2)/func_data(func_indx(n) +2))

				

		



	  ! 25 - first derivative of Morlet wavelet 

	  ! (http://www.eso.org/projects/esomidas/doc/user/98NOV/volb/node312.html) 

	  ! PLANE WAVE 

	  !

	  elseif (func_type(n).eq.25) then



		tt_t0 = t - func_data(func_indx(n) +1)





			PI = 4.0d0 * datan(1.0d0)



		    val = 1/(2*PI)**.5  * &

				 exp((-tt_t0**2)/func_data(func_indx(n) +2)) * &

				 ( - 2*tt_t0 / func_data(func_indx(n) +2) * &

				   cos(2*PI*func_data(func_indx(n))*tt_t0) - &

				   2 * PI *func_data(func_indx(n)) * &

				   sin(2*PI*func_data(func_indx(n))*tt_t0)) / &

				   (2* PI *func_data(func_indx(n)))

				

		



	  endif

      
      get_func_value = val
      
      return
      
      end function get_func_value
      
      
      
! *********************************************************
      
      subroutine find_nearest_node(n,xs,ys,xt,yt,nt)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      implicit none
      
      integer*4 :: n,nt
      real*8, dimension(*) :: xs,ys
      real*8 :: xt,yt
      
      integer*4 :: i
      real*8 :: dx,dy,d2,d2min
      
      d2min = 1.0d30
      
      do i = 1,n
         dx = xs(i) - xt
         dy = ys(i) - yt
         d2 = dx*dx + dy*dy
         if (d2.lt.d2min) then
            d2min = d2
            nt = i
         endif
      enddo
      
      return
      end subroutine find_nearest_node
      
! *********************************************************
         
      subroutine dime_sism_nodes(Xipo,Yipo,X1,Y1,&
				  X2,Y2,slip1,slip2,nnod,xs,ys,&
                                  node_sism)
         
!     © CRS4, 2002, All Rights Reserved
!     Authors: Marco Stupazzini
!

      implicit none
      
      integer*4 :: nnod,node_sism
      real*8, dimension(*) :: xs,ys
      integer*4 :: isn
      real*8 :: Xipo,Yipo,X1,Y1,X2,Y2,slip1,slip2,tol
      real*8 :: Xmax,Xmin,Ymax,Ymin,cost

      node_sism = 0
      tol = 1.0d-2
      cost = xipo*slip2-yipo*slip1

	 !Check on the relative position of fault end points
	 !In this preliminary part we check if P1 point is the lower-left point or not
	 !If this check would fail the two points are swapped
  
	 if (X1.gt.X2) then
	    Xmax=X1
	    Xmin=X2
	 else
	    Xmax=X2
            Xmin=X1
	 endif
         if (Y1.gt.Y2) then
            Ymax=Y1
            Ymin=Y2
         else
            Ymax=Y2 
            Ymin=Y1
         endif


  
	 do isn = 1,nnod

	    if (dabs(ys(isn)*slip1-xs(isn)*slip2+cost).le.tol) then
!		write(47,*)'dabs =',dabs(ys(isn)*slip1-xs(isn)*slip2+cost)
		if ((xs(isn).ge.Xmin).and.(xs(isn).le.Xmax)) then
!		   write(47,*)'CHECK sulla X'
		   if ((ys(isn).ge.Ymin).and.(ys(isn).le.Ymax)) then
!			write(47,*)'CHECK sulla Y' 
 			node_sism = node_sism+1
!                   write(47,*)'node_sism =',node_sism
!  		   write(47,*)'Node =',isn,' | x =',xs(isn),' | y =',ys(isn)
	          endif
	        endif
            endif
             
	 enddo
!	write(47,*)'node_sism =',node_sism

      return
      end subroutine dime_sism_nodes
      
! *********************************************************
         
      subroutine read_sism_nodes(Xipo,Yipo,X1,Y1,&
				  X2,Y2,slip1,slip2,nnod,xs,ys,&
                                  num_node_sism,sour_node_sism,i,&
				  dist_sour_node_sism,nl_sism,&
				  max_num_node_sism)
         
!     © CRS4, 2002, All Rights Reserved
!     Authors: Marco Stupazzini
!

      implicit none
      
      integer*4 :: nnod,node_sism,num_node_sism,nl_sism,max_num_node_sism
      real*8, dimension(*) :: xs,ys
      integer*4, dimension(max_num_node_sism,nl_sism) :: sour_node_sism
      real*8, dimension(max_num_node_sism,nl_sism) :: dist_sour_node_sism
      integer*4 :: i,isn
      real*8 :: Xipo,Yipo,X1,Y1,X2,Y2,slip1,slip2,tol
      real*8 :: Xmax,Xmin,Ymax,Ymin,cost

      node_sism = 0
      tol = 1.0d-2
      cost = xipo*slip2-yipo*slip1

      

	 !Check on the relative position of fault end points
	 !In this preliminary part we check if P1 point is the lower-left point or not
	 !If this check would fail the two points are swapped
  
	 if (X1.gt.X2) then
            Xmax=X1
            Xmin=X2
         else
            Xmax=X2 
            Xmin=X1
         endif
         if (Y1.gt.Y2) then
            Ymax=Y1
            Ymin=Y2
         else
            Ymax=Y2
            Ymin=Y1
         endif
 
   
	 !write(47,*)'READ_SISM_NODE'
         !write(47,*)'X1 =',X1,' | X2 =',X2
         !write(47,*)'Y1 =',Y1,' | Y2 =',Y2

	 do isn = 1,nnod
	    if (dabs(ys(isn)*slip1-xs(isn)*slip2+cost).le.tol) then
		if ((xs(isn).ge.Xmin).and.(xs(isn).le.Xmax)) then
		   if ((ys(isn).ge.Ymin).and.(ys(isn).le.Ymax)) then
                      node_sism = node_sism + 1	
                      sour_node_sism(node_sism,i) = isn
		      dist_sour_node_sism(node_sism,i) = sqrt((Xipo - xs(isn))**2 +&
			 (Yipo - ys(isn))**2)
                      
                      !write(47,*)'sour_node_sism =',sour_node_sism(node_sism,i),' |dist',&
		      !		dist_sour_node_sism(node_sism,i)
	           endif
	        endif
            endif
             
	 enddo

	!write(47,*)'node_sism =',node_sism
     
      return
      end subroutine read_sism_nodes
      

! *********************************************************
                
      subroutine check_sism(cs_nnz,cs,&
                          nm,tag_mat,sdeg_mat,&
                          ne,&
                          nl_sism,&
                          num_node_sism,max_num_node_sism,&
                          sour_node_sism,dist_sour_node_sism,&
                          check_node_sism,check_dist_node_sism,&
                          length_cns,&
                          fun_sism,nf,tag_func,valsism)


!     © CRS4, 2002, All Rights Reserved
!     Authors: Marco Stupazzini
!
         
      implicit none

      integer*4 :: cs_nnz,nm,ne,nl_sism,nf
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tag_mat,sdeg_mat
      integer*4, dimension(nl_sism) :: num_node_sism
      integer*4, dimension(nl_sism) :: fun_sism
      real*8, dimension(nl_sism,12) :: valsism
      real*8 :: vel_rup 
      integer*4, dimension(nf) :: tag_func
      integer*4 :: max_num_node_sism
      integer*4, dimension(max_num_node_sism,nl_sism) :: sour_node_sism
      real*8, dimension(max_num_node_sism,nl_sism) :: dist_sour_node_sism
      integer*4 :: length_cns
      integer*4, dimension(length_cns,5) :: check_node_sism
      real*8, dimension(length_cns,1) :: check_dist_node_sism
      integer*4 :: im,ie,isism,nn
      integer*4 :: is,in
      integer*4 :: i,j,k,h
      integer*4 :: fn,fun_sism_k
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8, dimension(:), allocatable :: dxdx,dydx,dxdy,dydy
      real*8, dimension(:,:), allocatable :: det_j

      !write(50,*),length_cns
      nn = 2
      h = 1

      allocate(ct(nn),ww(nn),dd(nn,nn))
      call lgl(nn,ct,ww,dd)
      allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))
   
      ne = cs(0) -1  
         
      do im = 1,nm
         if ((sdeg_mat(im) +1).ne.nn) then
            deallocate(ct,ww,dd)
            deallocate(dxdx,dydx,dxdy,dydy,det_j)
         
            nn = sdeg_mat(im) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call lgl(nn,ct,ww,dd)
            allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn),det_j(nn,nn))   
         endif

         do ie = 1,ne
            if (cs(cs(ie -1) +0).eq.tag_mat(im)) then
            
            !h = 1
                                                
            do isism = 1,nl_sism
               do fn = 1,nf
                  if (tag_func(fn).eq.fun_sism(isism)) then
                     vel_rup = valsism(isism,12)
                     fun_sism_k = fn
                  endif
               enddo
               do k = 1,num_node_sism(isism)
                  do j = 1,nn   
                     do i = 1,nn
                        is = nn*(j -1) +i
                        in = cs(cs(ie -1) +is)
                        if (in.eq.sour_node_sism(k,isism)) then
                            check_node_sism(h,1) = ie
                            check_node_sism(h,2) = i
                            check_node_sism(h,3) = j
                            check_node_sism(h,4) = isism
                            check_node_sism(h,5) = fun_sism_k
                            !check_node_sism(h,5) = 5
                            check_dist_node_sism(h,1) = dist_sour_node_sism(k,isism) &
                                                        / vel_rup

                            !write(50,*),check_node_sism(h,1),' | ',&
                            !            check_node_sism(h,2),' | ',&
                            !            check_node_sism(h,3),' | ',&
                            !            check_node_sism(h,4),' | ',&
                            !            check_node_sism(h,5),' | ',& 
                            !            check_dist_node_sism(h,1),' | ',&
                            !            vel_rup
                            h = h + 1
                        endif
                     enddo
                  enddo

               enddo
            enddo                           

! Seismic moment scale factor - end
                                    
      endif
    enddo
  enddo

      !write(50,*),h
      !do isism =  1,h-1
      !	                             write(52,*),check_node_sism(isism,1),' | ',&
      !                                  check_node_sism(isism,2),' | ',&
      !                                  check_node_sism(isism,3),' | ',&
      !                                  check_node_sism(isism,4),' | ',&
      !                                  check_dist_node_sism(isism,1),' | '
      !enddo
       
      return
      end subroutine check_sism

      
!     *************************************************************
      
      subroutine make_node_index(cs_nnz,cs,nm,sd,&
                                 nnode_macro,node_d,nelem,elem_d,&
                                 ndom,myid,nnode,node_i)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      implicit none
      
      integer*4 :: cs_nnz,nm,nnode_macro,nelem,ndom,myid,nnode
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: sd
      integer*4, dimension(nnode_macro) :: node_d
      integer*4, dimension(nelem) :: elem_d
      integer*4, dimension(nnode) :: node_i
      
      integer*4 :: i,j,in,is,ie,ic,id,nn,ai,bi
      
      
      do in = 1,nnode
         node_i(in) = 0
      enddo
      
      nelem = cs(0) -1
      
!     Vertices
      do in = 1,nnode_macro
         if (node_d(in).eq.myid) then
            node_i(in) = 1
         endif
      enddo
      
      do ie = 1,nelem
         nn = sd(cs(cs(ie -1) +0)) +1
         
         if (nn.gt.2) then
            
!           First edge
            ai = cs(cs(ie -1) +1)
            bi = cs(cs(ie -1) +nn)
            
            if (node_d(ai).eq.node_d(bi)) then
               id = node_d(ai)
            else
               id = elem_d(ie)
            endif
            
            if (id.eq.myid) then
               do i = 2,nn -1
                  is = nn*0 +i
                  in = cs(cs(ie -1) +is)
                  node_i(in) = 1
               enddo
            endif
            
!           Second edge
            ai = cs(cs(ie -1) +nn)
            bi = cs(cs(ie -1) +nn*nn)
            
            if (node_d(ai).eq.node_d(bi)) then
               id = node_d(ai)
            else
               id = elem_d(ie)
            endif
            
            if (id.eq.myid) then
               do i = 2,nn -1
                  is = nn*(i -1) +nn
                  in = cs(cs(ie -1) +is)
                  node_i(in) = 1
               enddo
            endif
            
!           Third edge
            ai = cs(cs(ie -1) +nn*nn)
            bi = cs(cs(ie -1) +nn*(nn -1) +1)
            
            if (node_d(ai).eq.node_d(bi)) then
               id = node_d(ai)
            else
               id = elem_d(ie)
            endif
            
            if (id.eq.myid) then
               do i = 2,nn -1
                  is = nn*(nn -1) +i
                  in = cs(cs(ie -1) +is)
                  node_i(in) = 1
               enddo
            endif
            
!           Fourth edge
            ai = cs(cs(ie -1) +nn*(nn -1) +1)
            bi = cs(cs(ie -1) +1)
            
            if (node_d(ai).eq.node_d(bi)) then
               id = node_d(ai)
            else
               id = elem_d(ie)
            endif
            
            if (id.eq.myid) then
               do i = 2,nn -1
                  is = nn*(i -1) +1
                  in = cs(cs(ie -1) +is)
                  node_i(in) = 1
               enddo
            endif
            
!           Internal nodes
            if (elem_d(ie).eq.myid) then
               do j = 2,nn -1
                  do i = 2,nn -1
                     is = nn*(j -1) +i
                     in = cs(cs(ie -1) +is)
                     node_i(in) = 1
                  enddo
               enddo
            endif
         endif
      enddo
      
      ic = 1
      do in = 1,nnode
         if (node_i(in).ne.0) then
            node_i(in) = ic
            ic = ic +1
         endif
      enddo
      
      return
      
      end subroutine make_node_index
      
      
!     *****************************************************************
      
      subroutine read_dime_header(file_head,nmonitors,nsnapshots,nsnapsprestress) ! Marco
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!
      implicit none
      
      integer*4 :: nmonitors,nsnapshots,nsnapsprestress ! Marco
      character*70 :: file_head
      
      character*80 :: input_line
      character*8 :: keyword
      integer*4 :: status
      
      
      nmonitors = 0
      nsnapshots = 0

	  nsnapsprestress = 0 ! Marco
      
      open(20,file=file_head)
      
      do
         read(20,'(A)',IOSTAT = status) input_line
         
         if (status.ne.0) exit
         
         keyword = input_line(1:8)
         
         if (keyword(1:7).eq.'MONITOR') then
            nmonitors = nmonitors + 1
         elseif (keyword(1:8).eq.'SNAPSHOT') then

            nsnapshots = nsnapshots + 1

	     elseif (keyword(1:8).eq.'SNAPPRES') then ! Marco

            nsnapsprestress = nsnapsprestress + 1 ! Marco
         endif
      enddo
      
      close(20)
      
      return
      end subroutine read_dime_header
      
      
      
!     *****************************************************************
      
      subroutine read_header(file_head,file_grid,file_mat,file_out,&
                             spectral_degree,time_step,stop_time,&
                             option_out_data,option_out_form,&
                             nmonitors,x_monitor,y_monitor,&
                             nsnapshots,t_snapshot,ndt_monitor,&

							 nsnapshotprestress,t_snapshotprestress) ! Marco
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     Modified: Marco Stupazzini (16 Feb 2003)
!     ndt_monitor
!
     
      implicit none
      
      integer*4 :: spectral_degree,nmonitors,nsnapshots

      integer*4 :: nsnapshotprestress ! Marco
      integer*4 :: option_out_data,option_out_form
      real*8 :: time_step,stop_time
      real*8, dimension(*) :: x_monitor,y_monitor,t_snapshot

	  real*8, dimension(*) :: t_snapshotprestress ! Marco
      character*70 :: file_head,file_grid,file_mat,file_out
      
      character*80 :: input_line
      character*8 :: keyword
      integer*4 :: status
      integer*4 :: ileft,iright
      integer*4 :: i,j,im,is

	  integer*4 :: ips ! Marco
      real*8 :: val
      
      character*3 :: deltat_fixed

      real*8 :: ndt_monitor
     
      
      ndt_monitor = 1
      im = 0
      is = 0

	  ips = 0 ! Marco
      
      open(20,file=file_head)
      
      do
         read(20,'(A)',IOSTAT = status) input_line
         
         if (status.ne.0) exit
         
         keyword = input_line(1:8)
         
         ileft = 0
         iright = len(input_line)
         do i = 1,iright
            if (input_line(i:i).eq.' ') exit
         enddo
         ileft = i
         
         if (keyword(1:8).eq.'GRIDFILE') then
            read(input_line(ileft:iright),*) file_grid
         elseif (keyword(1:7).eq.'MATFILE') then
            read(input_line(ileft:iright),*) file_mat
         elseif (keyword(1:7).eq.'OUTFILE') then
            read(input_line(ileft:iright),*) file_out
         elseif (keyword(1:6).eq.'DEGREE') then
            read(input_line(ileft:iright),*) spectral_degree
         elseif (keyword(1:8).eq.'TIMESTEP') then
            read(input_line(ileft:iright),*) time_step
	        deltat_fixed = 'not'
		 elseif (keyword(1:8).eq.'TIMEFIXE') then
            read(input_line(ileft:iright),*) time_step
	        deltat_fixed = 'yes'
         elseif (keyword(1:8).eq.'STOPTIME') then
            read(input_line(ileft:iright),*) stop_time
         elseif (keyword(1:7).eq.'MONITOR') then
            im = im +1
            read(input_line(ileft:iright),*) x_monitor(im),y_monitor(im)
         elseif (keyword(1:8).eq.'SNAPSHOT') then
            is = is +1
            read(input_line(ileft:iright),*) t_snapshot(is)

		 elseif (keyword(1:8).eq.'SNAPPRES') then

            ips = ips +1

            read(input_line(ileft:iright),*) t_snapshotprestress(ips)
         elseif (keyword(1:7).eq.'OPTIOUT') then
            read(input_line(ileft:iright),*) option_out_form,option_out_data
	 elseif (keyword(1:8).eq.'TMONITOR') then
            read(input_line(ileft:iright),*) ndt_monitor
	        !deltat_monitor = 'yes'
         endif
      enddo
      
      
! Riordino gli snapshots
      
      if (nsnapshots.gt.1) then
         do i = 1,nsnapshots-1
            do j = i+1, nsnapshots
               if (t_snapshot(i).gt.t_snapshot(j)) then
                  val = t_snapshot(i)
                  t_snapshot(i) = t_snapshot(j)
                  t_snapshot(j) = val
               endif
            enddo
         enddo
      endif



	  if (nsnapshotprestress.gt.1) then

         do i = 1,nsnapshotprestress-1

            do j = i+1, nsnapshotprestress

               if (t_snapshotprestress(i).gt.t_snapshotprestress(j)) then

                  val = t_snapshotprestress(i)

                  t_snapshotprestress(i) = t_snapshotprestress(j)

                  t_snapshotprestress(j) = val

               endif

            enddo

         enddo

      endif
      
      close(20)
      
      return
      end subroutine read_header
      
      
      
!     *****************************************************************
      
      subroutine read_dime_mat_el(file_mat,nmat, &

								  nmatg, &  ! Mexico Paco grad_lin 29.11.2004
                                  nload_dirX,nload_dirY, &
                                  nload_neuX,nload_neuY, &

								  nload_neuN,nload_neuT, & ! M Mexico Arroyo
                                  nload_poiX,nload_poiY, &

								  nload_palX,nload_palY, & ! Marco

								  nload_dipX,nload_dipY, & ! Marco
                                  nload_plaX,nload_plaY, &
                                  nload_sism, &
                                  nload_vpcl, & ! Clara

								  nload_vpsa, & ! MCNAP

								  nload_vpsl, & ! MCNAP

								  nload_vpsd, & ! MCNAP 

								  nload_maps, & ! Marco 
                                  nload_abc, &
                                  nfunc,nfunc_data)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     Modified by: Marco Stupazzini (02 august 2002)
!     Modification: Plane Wave input
!
      implicit none
      
      integer*4 :: nmat

	  integer*4 :: nmatg  ! Mexico Paco grad_lin 29.11.2004
      integer*4 :: nload_dirX,nload_dirY,nload_neuX,nload_neuY

	  integer*4 :: nload_neuN,nload_neuT ! M Mexico Arroyo
      integer*4 :: nload_poiX,nload_poiY,nload_dipX,nload_dipY,nload_plaX,nload_plaY,nload_sism,nload_vpcl,nload_abc ! Clara

	  integer*4 :: nload_vpsa,nload_vpsl,nload_vpsd ! MCNAP

	  integer*4 :: nload_palX,nload_palY ! Marco

	  integer*4 :: nload_maps ! Marco 
      integer*4 :: nfunc,nfunc_data
      character*70   :: file_mat
      
      character*80 :: input_line
      character*4 :: keyword
      integer*4 :: status
      integer*4 :: tag_func,func_type,func_nd
      
      nmat = 0

	  nmatg = 0  ! Mexico Paco grad_lin 29.11.2004
      nload_dirX = 0
      nload_dirY = 0
      nload_neuX = 0

      nload_neuN = 0 ! M Mexico Arroyo
      nload_neuY = 0

	  nload_neuT = 0 ! M Mexico Arroyo 
      nload_poiX = 0 
      nload_poiY = 0

	  nload_palX = 0 ! Marco 

      nload_palY = 0 ! Marco

	  nload_dipX = 0 ! Marco 

      nload_dipY = 0 ! Marco
      nload_plaX = 0
      nload_plaY = 0 
      nload_sism = 0
	  nload_vpcl = 0 ! Clara

	  nload_vpsa = 0 ! MCNAP

	  nload_vpsl = 0 ! MCNAP

	  nload_vpsd = 0 ! MCNAP

	  nload_maps = 0 ! Marco 
      nload_abc = 0
      nfunc = 0
      nfunc_data = 0
      
      open(23,file=file_mat)
      
      do
         read(23,'(A)',IOSTAT = status) input_line
         
         if (status.ne.0) exit
         
         keyword = input_line(1:4)
         
         if (keyword.eq.'MATE') then
            nmat = nmat + 1

	     elseif (keyword.eq.'MATG') then  ! Mexico Paco grad_lin 29.11.2004

            nmatg = nmatg + 1  ! Mexico Paco grad_lin 29.11.2004
         elseif (keyword.eq.'DIRX') then
            nload_dirX = nload_dirX + 1
         elseif (keyword.eq.'DIRY') then
            nload_dirY = nload_dirY + 1
         elseif (keyword.eq.'NEUX') then
            nload_neuX = nload_neuX + 1
         elseif (keyword.eq.'NEUY') then
            nload_neuY = nload_neuY + 1

		 elseif (keyword.eq.'NEUN') then ! M Mexico Arroyo

            nload_neuN = nload_neuN + 1 ! M Mexico Arroyo

         elseif (keyword.eq.'NEUT') then ! M Mexico Arroyo

            nload_neuT = nload_neuT + 1 ! M Mexico Arroyo
         elseif (keyword.eq.'PLOX') then
            nload_poiX = nload_poiX + 1
         elseif (keyword.eq.'PLOY') then
            nload_poiY = nload_poiY + 1

	     elseif (keyword.eq.'PALX') then ! Marco

            nload_palX = nload_palX + 1 ! Marco

         elseif (keyword.eq.'PALY') then ! Marco

            nload_palY = nload_palY + 1 ! Marco

		 elseif (keyword.eq.'DIPX') then ! Marco

            nload_dipX = nload_dipX + 1  ! Marco

         elseif (keyword.eq.'DIPY') then ! Marco

            nload_dipY = nload_dipY + 1  ! Marco
         elseif (keyword.eq.'PLAX') then
            nload_plaX = nload_plaX + 1
         elseif (keyword.eq.'PLAY') then
            nload_plaY = nload_plaY + 1
         elseif (keyword.eq.'SISM') then
            nload_sism = nload_sism + 1
	     elseif (keyword.eq.'VPCL') then ! Clara 
            nload_vpcl = nload_vpcl + 1 ! Clara

	     elseif (keyword.eq.'VPSA') then ! MCNAP 

            nload_vpsa = nload_vpsa + 1 ! MCNAP

		 elseif (keyword.eq.'VPSL') then ! MCNAP 

            nload_vpsl = nload_vpsl + 1 ! MCNAP

		 elseif (keyword.eq.'VPSD') then ! MCNAP 

            nload_vpsd = nload_vpsd + 1 ! MCNAP

	     elseif (keyword.eq.'MAPS') then ! Marco 

            nload_maps = nload_maps + 1 ! Marco
         elseif (keyword.eq.'ABSO') then
            nload_abc = nload_abc + 1
         elseif (keyword.eq.'FUNC') then
            nfunc = nfunc + 1
            
            read(input_line(5:),*)tag_func,func_type
            
            if (func_type.eq.0) then
               nfunc_data = nfunc_data + 0
            elseif (func_type.eq.1) then
               nfunc_data = nfunc_data + 2
            elseif (func_type.eq.2) then
               nfunc_data = nfunc_data + 2
            elseif (func_type.eq.3) then
               read(input_line(5:),*)tag_func,func_type,func_nd
               nfunc_data = nfunc_data + 2*func_nd
            elseif (func_type.eq.4) then
               nfunc_data = nfunc_data + 2
		    elseif (func_type.eq.5) then
               nfunc_data = nfunc_data + 2
		    elseif (func_type.eq.6) then
               nfunc_data = nfunc_data + 4
		    elseif (func_type.eq.7) then
               nfunc_data = nfunc_data + 4
		    elseif (func_type.eq.8) then
               nfunc_data = nfunc_data + 2
			elseif (func_type.eq.9) then
               nfunc_data = nfunc_data + 3
			elseif (func_type.eq.10) then
               nfunc_data = nfunc_data + 3
			elseif (func_type.eq.11) then
               nfunc_data = nfunc_data + 3
		    elseif (func_type.eq.12) then
               nfunc_data = nfunc_data + 3
            elseif (func_type.eq.14) then

               nfunc_data = nfunc_data + 2

            elseif (func_type.eq.20) then

               nfunc_data = nfunc_data + 3

            elseif (func_type.eq.21) then

               nfunc_data = nfunc_data + 3

            elseif (func_type.eq.22) then

               nfunc_data = nfunc_data + 3

            elseif (func_type.eq.23) then

               nfunc_data = nfunc_data + 3

            elseif (func_type.eq.24) then

               nfunc_data = nfunc_data + 3

            elseif (func_type.eq.25) then

               nfunc_data = nfunc_data + 3

            endif
         endif
      enddo
      
      close(23)
      
      return
      end subroutine read_dime_mat_el
      
      
      
! *********************************************************

      subroutine read_material_el(file_mat,nm,propm,typem,trefm,tagm, &

	                              nmg,valmg,typemg,tagmg, & ! Mexico Paco grad_lin 29.11.2004
                                  ndX,valdX,fdX,tagdX, &
                                  ndY,valdY,fdY,tagdY, &
                                  nnX,valnX,fnX,tagnX, &
                                  nnY,valnY,fnY,tagnY, &

								  nnN,valnN,fnN,tagnN, & ! M Mexico Arroyo

                                  nnT,valnT,fnT,tagnT, & ! M Mexico Arroyo
                                  npX,valpX,fpX, &
                                  npY,valpY,fpY, &

								  npaX,valpaX,fpaX,tagpaX, & ! Marco

                                  npaY,valpaY,fpaY,tagpaY, & ! Marco

								  ndpX,valdpX,fdpX, & ! Marco

                                  ndpY,valdpY,fdpY, & ! Marco
                                  nplX,valplX,fplX,tagplX, &
                                  nplY,valplY,fplY,tagplY, &
                                  nsism,valsism,fsism,tagsism, &
                                  nvpcl,valvpcl,fvpcl,tagvpcl, & ! Clara

								  nvpsa,valvpsa,tagvpsa, & ! MCNAP

								  nvpsl,valvpsl,tagvpsl, & ! MCNAP

								  nvpsd,valvpsd,tagvpsd, & ! MCNAP

								  nmaps,valmaps,fmaps,tagmaps, & ! Marco 
								  na,taga, &
                                  nf,func_type,func_indx,func_data,nfunc_data,tag_func, &
                                  fmax, &
                                  err_out)
      
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     Modified by: Marco Stupazzini (02 august 2002)
!     Modification: Plane Wave input
!     
      implicit none
      
      character*70 :: file_mat      
      integer*4 :: nm,ndT,ncT,ndX,ndY,ndpX,ndpY,nnX,nnY,npX,npY,nplX,nplY,nsism,nvpcl,nmaps,na,npdT,nf ! Clara&Marco
      integer*4 :: nmg  ! Mexico Paco grad_lin 29.11.2004

	  integer*4 :: nnN,nnT ! M Mexico Arroyo

	  integer*4 :: nvpsa,nvpsl,nvpsd ! MCNAP

	  integer*4 :: npaX,npaY ! Marco

	  character*12, dimension(nm) :: namem
      integer*4, dimension(nm) :: typem
      real*8, dimension(nm) :: trefm
      real*8, dimension(nm,4) :: propm
      integer*4, dimension(nm) :: tagm





      real*8, dimension(nmg,*) :: valmg ! Mexico Paco grad_lin 29.11.2004

	  integer*4, dimension(nmg) :: typemg ! Mexico Paco grad_lin 29.11.2004

      integer*4, dimension(nmg) :: tagmg ! Mexico Paco grad_lin 29.11.2004

	  real*8 :: Vs1,Vs2,Vp1,Vp2 ! Mexico Paco grad_lin 29.11.2004

	  real*8 :: alpha_Vs,beta_Vs,alpha_Vp,beta_Vp ! Mexico Paco grad_lin 29.11.2004
      

	  integer*4 :: nfunc_data
      integer*4, dimension(nf) :: func_type
      integer*4, dimension(nf +1) :: func_indx
      real*8, dimension(nfunc_data) :: func_data
      integer*4, dimension(nf) :: tag_func
      
      real*8, dimension(ndX,*) :: valdX
      real*8, dimension(ndY,*) :: valdY
      integer*4, dimension(*) :: fdX,tagdX,fdY,tagdY
      
      real*8, dimension(nnX,*) :: valnX
      real*8, dimension(nnY,*) :: valnY
      integer*4, dimension(*) :: fnX,tagnX,fnY,tagnY



	  real*8, dimension(nnN,*) :: valnN ! M Mexico Arroyo

      real*8, dimension(nnT,*) :: valnT ! M Mexico Arroyo

      integer*4, dimension(*) :: fnN,tagnN,fnT,tagnT ! M Mexico Arroyo
      
      real*8, dimension(npX,*) :: valpX
      real*8, dimension(npY,*) :: valpY
      integer*4, dimension(*) :: fpX,fpY



	  real*8, dimension(npaX,*) :: valpaX ! Marco

      real*8, dimension(npaY,*) :: valpaY ! Marco

      integer*4, dimension(*) :: fpaX,tagpaX,fpaY,tagpaY ! Marco



	  real*8, dimension(ndpX,*) :: valdpX ! Marco

      real*8, dimension(ndpY,*) :: valdpY ! Marco

      integer*4, dimension(*) :: fdpX,fdpY ! Marco

      real*8, dimension(nplX,*) :: valplX
      real*8, dimension(nplY,*) :: valplY
      integer*4, dimension(*) :: fplX,tagplX,fplY,tagplY
      
      real*8, dimension(nsism,*) :: valsism  
      integer*4, dimension(*) :: fsism,tagsism

	  real*8, dimension(nvpcl,*) :: valvpcl   ! Clara 
      integer*4, dimension(*) :: fvpcl,tagvpcl ! Clara 



	  real*8, dimension(nvpsa,*) :: valvpsa   ! MCNAP 

      integer*4, dimension(*) :: tagvpsa ! MCNAP



	  real*8, dimension(nvpsl,*) :: valvpsl   ! MCNAP 

      integer*4, dimension(*) :: tagvpsl ! MCNAP



	  real*8, dimension(nvpsd,*) :: valvpsd   ! MCNAP 

      integer*4, dimension(*) :: tagvpsd ! MCNAP



	  real*8, dimension(nmaps,*) :: valmaps   ! Marco 

      integer*4, dimension(*) :: fmaps,tagmaps ! Marco 

      integer*4, dimension(*) :: taga
      
      real*8 :: fmax
      
      integer*4 :: im,ifunc,idf,func_nd

	  integer*4 :: img  ! Mexico Paco grad_lin 29.11.2004
      integer*4 :: idX,idY,idpX,idpY,inX,inY,ipX,ipY,iplX,iplY,isism,ivpcl,imaps,ia ! Clara&Marco
      integer*4 :: inN,inT ! M Mexico Arroyo

	  integer*4 :: ivpsa,ivpsl,ivpsd ! MCNAP

	  integer*4 :: ipaX,ipaY ! Marco

	  integer*4 :: ileft,iright
      integer*4 :: i,j,err_out,trash,status
      character*100000 :: input_line
      character*4 :: keyword
      
      
      open(23,file=file_mat)
      
      err_out=0
      
      im = 0

	  img = 0  ! Mexico Paco grad_lin 29.11.2004
      idX = 0
      idY = 0
      inX = 0
      inY = 0

	  inN = 0 ! M Mexico Arroyo

      inT = 0 ! M Mexico Arroyo
      ipX = 0
      ipY = 0

	  ipaX = 0 ! Marco

      ipaY = 0 ! Marco

	  idpX = 0 ! Marco

      idpY = 0 ! Marco
      iplX = 0
      iplY = 0
      isism = 0
	  ivpcl = 0 ! Clara

	  ivpsa = 0 ! MCNAP

	  ivpsl = 0 ! MCNAP

	  ivpsd = 0 ! MCNAP

	  imaps = 0 ! Marco 
      ia = 0
      ifunc = 0
      
      if (nf.gt.0) then
         func_indx(1) = 1
      endif
      
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
         
         if (keyword.eq.'MATE') then
            im = im + 1
            read(input_line(ileft:iright),*) tagm(im),typem(im),&
                 propm(im,1),propm(im,2),propm(im,3),propm(im,4)
!

		 elseif (keyword.eq.'MATG') then



		    Vs1 = 0.0d0

			Vs2 = 0.0d0

			Vp1 = 0.0d0

			Vp2 = 0.0d0



			alpha_Vs = 0.0d0

			beta_Vs = 0.0d0

			alpha_Vp = 0.0d0

			beta_Vp = 0.0d0



            img = img + 1

            read(input_line(ileft:iright),*) tagmg(img),typemg(img),&

                 valmg(img,1),valmg(img,2),valmg(img,3),valmg(img,4),&

				 valmg(img,5),valmg(img,6),valmg(img,7),valmg(img,8),&

				 valmg(img,9),valmg(img,10),valmg(img,11),valmg(img,12)



				 valmg(img,13) = 0.0d0

				 valmg(img,14) = 0.0d0

				 valmg(img,15) = 0.0d0

				 valmg(img,16) = 0.0d0

				 valmg(img,17) = 0.0d0

				 valmg(img,18) = 0.0d0





                 !***********************************************************

				 Vs1 = (valmg(img,3)/valmg(img,1))**.5

				 Vp1 = ((valmg(img,2) + 2 * valmg(img,3))/valmg(img,1))**.5

				 !y1 = propmg(img,6)

				 Vs2 = (valmg(img,9)/valmg(img,7))**.5

				 Vp2 = ((valmg(img,8) + 2 * valmg(img,9))/valmg(img,7))**.5

				 !y2 = propmg(img,12)



				 !***********************************************************



				 if (typemg(img).eq.1) then

                 

					alpha_Vs = (Vs2 - Vs1)/(valmg(img,12) - valmg(img,6))

					beta_Vs = Vs1 - valmg(img,6) * alpha_Vs



					!rho_alpha2 = propmg(img,1) * alpha_vs **2

					!2rho_alpha_beta = propmg(img,1) * alpha_vs * beta_vs

					!rho_beta2 = propmg(img,1) * beta_vs **2



					valmg(img,13) = valmg(img,1) * alpha_Vs **2

					valmg(img,14) = 2 * valmg(img,1) * alpha_Vs * beta_Vs

					valmg(img,15) = valmg(img,1) * beta_Vs **2



					!***********************************************************



					alpha_Vp = (Vp2 - Vp1)/(valmg(img,12) - valmg(img,6))

					beta_Vp = Vp1 - valmg(img,6) * alpha_Vp



					valmg(img,16) = valmg(img,1) * alpha_Vp **2

					valmg(img,17) = 2 * valmg(img,1) * alpha_Vp * beta_Vp

					valmg(img,18) = valmg(img,1) * beta_Vp **2



					!***********************************************************

				

				elseif (typemg(img).eq.2) then



				   valmg(img,13) = Vs2 - Vs1

                   valmg(img,14) = Vs1



				   valmg(img,15) = Vp2 - Vp1

                   valmg(img,16) = Vp1



				   !valmg(img,17) = 0.05d0

				   valmg(img,17) = valmg(img,11)



				   valmg(img,18) = - valmg(img,6) - (valmg(img,12) - valmg(img,6))/2



				endif



!
         elseif (keyword.eq.'DIRX') then
            idX = idX + 1
            read(input_line(ileft:iright),*)tagdX(idX),fdX(idX),&
                 valdX(idX,1),valdX(idX,2)
!
         elseif (keyword.eq.'DIRY') then
            idY = idY + 1
            read(input_line(ileft:iright),*)tagdY(idY),fdY(idY),&
                 valdY(idY,1),valdY(idY,2)
!
         elseif (keyword.eq.'NEUX') then
            inX = inX + 1
            read(input_line(ileft:iright),*)tagnX(inX),fnX(inX),&
                 valnX(inX,1),valnX(inX,2)
!
         elseif (keyword.eq.'NEUY') then
            inY = inY + 1
            read(input_line(ileft:iright),*)tagnY(inY),fnY(inY),&
                 valnY(inY,1),valnY(inY,2)

!

         elseif (keyword.eq.'NEUN') then ! M Mexico Arroyo

            inN = inN + 1 ! M Mexico Arroyo

            read(input_line(ileft:iright),*)tagnN(inN),fnN(inN),& ! M Mexico Arroyo

                 valnN(inN,1),valnN(inN,2) ! M Mexico Arroyo

!

         elseif (keyword.eq.'NEUT') then ! M Mexico Arroyo

            inT = inT + 1 ! M Mexico Arroyo

            read(input_line(ileft:iright),*)tagnT(inT),fnT(inT),& ! M Mexico Arroyo

                 valnT(inT,1),valnT(inT,2) ! M Mexico Arroyo
!
         elseif (keyword.eq.'PLOX') then
            ipX = ipX + 1
            read(input_line(ileft:iright),*)fpX(ipX),&
                 valpX(ipX,1),valpX(ipX,2),valpX(ipX,3)
!
         elseif (keyword.eq.'PLOY') then
            ipY = ipY + 1
            read(input_line(ileft:iright),*)fpY(ipY),&
                 valpY(ipY,1),valpY(ipY,2),valpY(ipY,3)

!

         elseif (keyword.eq.'PALX') then

            ipaX = ipaX + 1

            read(input_line(ileft:iright),*)tagpaX(ipaX),fpaX(ipaX),&

                 valpaX(ipaX,1)

!

         elseif (keyword.eq.'PALY') then

            ipaY = ipaY + 1

            read(input_line(ileft:iright),*)tagpaY(ipaY),fpaY(ipaY),&

                 valpaY(ipaY,1)

!

         elseif (keyword.eq.'DIPX') then ! Marco

            idpX = idpX + 1 ! Marco

            read(input_line(ileft:iright),*)fdpX(idpX),& ! Marco

                 valdpX(idpX,1),valdpX(idpX,2),valdpX(idpX,3) ! Marco

!

         elseif (keyword.eq.'DIPY') then ! Marco

            idpY = idpY + 1 ! Marco

            read(input_line(ileft:iright),*)fdpY(idpY),& ! Marco

                 valdpY(idpY,1),valdpY(idpY,2),valdpY(idpY,3) ! Marco

!
         elseif (keyword.eq.'PLAX') then
            iplX = iplX + 1
            read(input_line(ileft:iright),*)fplX(iplX),&
                 tagplX(iplX),valplX(iplX,1)
!

        elseif (keyword.eq.'PLAY') then
            iplY = iplY + 1
            read(input_line(ileft:iright),*)fplY(iplY),&
                 tagplY(iplY),valplY(iplY,1)
!
                 
        elseif (keyword.eq.'SISM') then
            isism = isism + 1
            read(input_line(ileft:iright),*)fsism(isism),&
                 tagsism(isism),valsism(isism,1),valsism(isism,2),&
		 valsism(isism,3),valsism(isism,4),valsism(isism,5),&
		 valsism(isism,6),valsism(isism,7),valsism(isism,8),&
		 valsism(isism,9),valsism(isism,10),valsism(isism,11),&
		 valsism(isism,12)
!
                 
        elseif (keyword.eq.'VPCL') then ! Clara 
            ivpcl = ivpcl + 1 ! Clara 
            read(input_line(ileft:iright),*)tagvpcl(ivpcl),& ! Clara 
                 fvpcl(ivpcl),valvpcl(ivpcl,1),valvpcl(ivpcl,2),& ! Clara 
		 valvpcl(ivpcl,3),valvpcl(ivpcl,4),valvpcl(ivpcl,5), valvpcl(ivpcl,6) ! Clara  !Clara_non_ass
!

                 

        elseif (keyword.eq.'VPSA') then ! MCNAP 

            ivpsa = ivpsa + 1 ! MCNAP 

            read(input_line(ileft:iright),*)tagvpsa(ivpsa),& ! MCNAP 

                 valvpsa(ivpsa,1),valvpsa(ivpsa,2),& ! MCNAP 

		 valvpsa(ivpsa,3),valvpsa(ivpsa,4),valvpsa(ivpsa,5) ! MCNAP

!

                 

        elseif (keyword.eq.'VPSL') then ! MCNAP 

            ivpsl = ivpsl + 1 ! MCNAP 

            read(input_line(ileft:iright),*)tagvpsl(ivpsl),& ! MCNAP 

                 valvpsl(ivpsl,1),valvpsl(ivpsl,2),& ! MCNAP 

		 valvpsl(ivpsl,3),valvpsl(ivpsl,4),valvpsl(ivpsl,5), valvpsl(ivpsl,6),& ! MCNAP

		 valvpsl(ivpsl,7),valvpsl(ivpsl,8),valvpsl(ivpsl,9), valvpsl(ivpsl,10),& ! MCNAP

		 valvpsl(ivpsl,11),valvpsl(ivpsl,12),valvpsl(ivpsl,13), valvpsl(ivpsl,14),& ! MCNAP

		 valvpsl(ivpsl,15),valvpsl(ivpsl,16) ! MCNAP



!

                 

        elseif (keyword.eq.'VPSD') then ! MCNAP 

            ivpsd = ivpsd + 1 ! MCNAP 

            read(input_line(ileft:iright),*)tagvpsd(ivpsd),& ! MCNAP 

                 valvpsd(ivpsd,1),valvpsd(ivpsd,2),& ! MCNAP 

		 valvpsd(ivpsd,3),valvpsd(ivpsd,4),valvpsd(ivpsd,5),valvpsd(ivpsd,6),& ! MCNAP

		 valvpsd(ivpsd,7),valvpsd(ivpsd,8),valvpsd(ivpsd,9),valvpsd(ivpsd,10),& ! MCNAP

		 valvpsd(ivpsd,11),valvpsd(ivpsd,12),valvpsd(ivpsd,13),valvpsd(ivpsd,14),& ! MCNAP

		 valvpsd(ivpsd,15),valvpsd(ivpsd,16) ! MCNAP



!

                 

        elseif (keyword.eq.'MAPS') then ! Marco 

            imaps = imaps + 1 ! Marco 

            read(input_line(ileft:iright),*)tagmaps(imaps),& ! Marco 

                 fmaps(imaps),valmaps(imaps,1) ! Marco 

!

         elseif (keyword.eq.'ABSO') then
            ia = ia + 1
            read(input_line(ileft:iright),*)taga(ia)
!
         elseif (keyword.eq.'FUNC') then
            ifunc = ifunc + 1
            read(input_line(ileft:iright),*)tag_func(ifunc),&
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
		    elseif (func_type(ifunc).eq.6) then
               func_indx(ifunc +1) = func_indx(ifunc) + 4
               read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
		    elseif (func_type(ifunc).eq.7) then
               func_indx(ifunc +1) = func_indx(ifunc) + 4
               read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
		    elseif (func_type(ifunc).eq.8) then
               func_indx(ifunc +1) = func_indx(ifunc) + 2
               read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
		    elseif (func_type(ifunc).eq.9) then
               func_indx(ifunc +1) = func_indx(ifunc) + 3
               read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
		    elseif (func_type(ifunc).eq.10) then
               func_indx(ifunc +1) = func_indx(ifunc) + 3
               read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
		    elseif (func_type(ifunc).eq.11) then
               func_indx(ifunc +1) = func_indx(ifunc) + 3
               read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
		    elseif (func_type(ifunc).eq.12) then
               func_indx(ifunc +1) = func_indx(ifunc) + 3
               read(input_line(ileft:iright),*)trash,trash,&
                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)
            elseif (func_type(ifunc).eq.14) then

               func_indx(ifunc +1) = func_indx(ifunc) + 2

               read(input_line(ileft:iright),*)trash,trash,&

                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)

            elseif (func_type(ifunc).eq.20) then

               func_indx(ifunc +1) = func_indx(ifunc) + 3

               read(input_line(ileft:iright),*)trash,trash,&

                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)

            elseif (func_type(ifunc).eq.21) then

               func_indx(ifunc +1) = func_indx(ifunc) + 3

               read(input_line(ileft:iright),*)trash,trash,&

                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)

            elseif (func_type(ifunc).eq.22) then

               func_indx(ifunc +1) = func_indx(ifunc) + 3

               read(input_line(ileft:iright),*)trash,trash,&

                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)

            elseif (func_type(ifunc).eq.23) then

               func_indx(ifunc +1) = func_indx(ifunc) + 3

               read(input_line(ileft:iright),*)trash,trash,&

                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)

            elseif (func_type(ifunc).eq.24) then

               func_indx(ifunc +1) = func_indx(ifunc) + 3

               read(input_line(ileft:iright),*)trash,trash,&

                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)

            elseif (func_type(ifunc).eq.25) then

               func_indx(ifunc +1) = func_indx(ifunc) + 3

               read(input_line(ileft:iright),*)trash,trash,&

                    (func_data(j), j = func_indx(ifunc),func_indx(ifunc +1) -1)

            endif
            
!
         elseif (keyword.eq.'FMAX') then
            read(input_line(ileft:iright),*)fmax
         endif
!
      enddo
      
      
      close(23)
      
      return
      
      end subroutine read_material_el
      
      
      
      
!     *****************************************************************
      
      subroutine read_dime_grid_el(file_grid,nmat,tag_mat,&
                                   ndirX,tag_dirX,ndirY,tag_dirY,&
                                   nneuX,tag_neuX,nneuY,tag_neuY,&

								   nneuN,tag_neuN,nneuT,tag_neuT,& ! M Mexico Arroyo
                                   nabc,tag_abc,&
                                   nnode,nquad,nline)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      implicit none
      
      character*70 :: file_grid
      integer*4 :: nmat,ndirX,ndirY,nneuX,nneuY,nabc

	  integer*4 :: nneuN,nneuT ! M Mexico Arroyo
      integer*4, dimension(nmat) :: tag_mat
      integer*4, dimension(ndirX) :: tag_dirX
      integer*4, dimension(ndirY) :: tag_dirY
      integer*4, dimension(nneuX) :: tag_neuX
      integer*4, dimension(nneuY) :: tag_neuY

	  integer*4, dimension(nneuN) :: tag_neuN ! M Mexico Arroyo

      integer*4, dimension(nneuT) :: tag_neuT ! M Mexico Arroyo
      integer*4, dimension(nabc) :: tag_abc
      integer*4 :: nnode,nquad,nline
      
      integer*4 :: nelem,ie,i,mcode,ileft,iright,sl,trash,status,check
      character*80 :: input_line
      character*20 :: ecode
      
      
      nnode = 0
      nquad = 0
      nline = 0
      
      open(23,file=file_grid)
      
      do 
         read(23,'(A)') input_line
         if (input_line(1:1) .ne. '#') exit
      enddo
      
      read(input_line,*)nnode,nelem
      
      do i = 1,nnode
         read(23,'(A)')input_line
      enddo
      
      do ie = 1,nelem
         read(23,'(A)')input_line
         
         sl = len(input_line)
         ileft = 0
         iright = 0 
         do i = 1,sl
            if (input_line(i:i).ge.'A') exit
         enddo
         ileft = i
         do i = ileft,sl
            if (input_line(i:i).lt.'A') exit
         enddo
         iright = i
         
         ecode = input_line(ileft:iright)
         
         read(input_line(1:ileft),*)trash,mcode
         
         if ((ecode.eq.'quad').or.(ecode.eq.'QUAD')) then
            check = 0
            do i = 1,nmat
               if (tag_mat(i).eq.mcode) check = 1
            enddo
            !
            if (check.ne.0) nquad = nquad +1
         !
         elseif ((ecode.eq.'line').or.(ecode.eq.'LINE')) then
            check = 0
            do i = 1,ndirX
               if (tag_dirX(i).eq.mcode) check = 1
            enddo
            do i = 1,ndirY
               if (tag_dirY(i).eq.mcode) check = 1
            enddo
            do i = 1,nneuX
               if (tag_neuX(i).eq.mcode) check = 1
            enddo
            do i = 1,nneuY
               if (tag_neuY(i).eq.mcode) check = 1
            enddo

			do i = 1,nneuN ! M Mexico Arroyo

               if (tag_neuN(i).eq.mcode) check = 1 ! M Mexico Arroyo

            enddo ! M Mexico Arroyo

            do i = 1,nneuT ! M Mexico Arroyo

               if (tag_neuT(i).eq.mcode) check = 1 ! M Mexico Arroyo

            enddo ! M Mexico Arroyo
            do i = 1,nabc
               if (tag_abc(i).eq.mcode) check = 1
            enddo
            !
            if (check.ne.0) nline = nline +1
         endif
      enddo
      
      close(23)
      
      return
      end subroutine read_dime_grid_el
      
      
      
!     ********************************************************
      
      subroutine read_grid_el(file_grid,nmat,tag_mat,prop_mat,&
                              ndirX,tag_dirX,ndirY,tag_dirY,&
                              nneuX,tag_neuX,nneuY,tag_neuY,nabc,tag_abc,&

							  nneuN,tag_neuN,nneuT,tag_neuT,& ! M Mexico Arroyo
                              nnode,xx,yy,nquad,con_quad,nline,con_line)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      implicit none
      
      character*70 :: file_grid
      integer*4 :: nmat,ndirX,ndirY,nneuX,nneuY,nabc

	  integer*4 :: nneuN,nneuT ! M Mexico Arroyo
      integer*4, dimension(nmat) :: tag_mat
      real*8, dimension(nmat,4) :: prop_mat
      integer*4, dimension(ndirX) :: tag_dirX
      integer*4, dimension(ndirY) :: tag_dirY
      integer*4, dimension(nneuX) :: tag_neuX
      integer*4, dimension(nneuY) :: tag_neuY

	  integer*4, dimension(nneuN) :: tag_neuN ! M Mexico Arroyo

      integer*4, dimension(nneuT) :: tag_neuT ! M Mexico Arroyo
      integer*4, dimension(nabc) :: tag_abc
      integer*4 :: nnode,nquad,nline
      real*8, dimension(nnode) :: xx,yy
      integer*4, dimension(nquad,5) :: con_quad
      integer*4, dimension(nline,3) :: con_line
      
      integer*4 :: inode,iquad,iline
      
      integer*4 :: nelem,ie,i,j,mcode,ileft,iright,sl,trash,status,check
      character*80 :: input_line
      character*20 :: ecode
      

      inode = 0
      iquad = 0
      iline = 0
      
      status = 0 

      open(23,file=file_grid)
      
      do 
         read(23,'(A)') input_line
         if (input_line(1:1) .ne. '#') exit
      enddo
      
      read(input_line,*)nnode,nelem
      
      do i = 1,nnode
        read(23,*)inode,xx(inode),yy(inode)
        if (inode.ne.i) then
          status = 1
        endif
      enddo
      
      do ie = 1,nelem
         read(23,'(A)')input_line
         
         sl = len(input_line)
         ileft = 0
         iright = 0 
         do i = 1,sl
            if (input_line(i:i).ge.'A') exit
         enddo
         ileft = i
         do i = ileft,sl
            if (input_line(i:i).lt.'A') exit
         enddo
         iright = i
         
         ecode = input_line(ileft:iright)
         
         read(input_line(1:ileft),*)trash,mcode
         
         if ((ecode.eq.'quad').or.(ecode.eq.'QUAD')) then
            check = 0
            do i = 1,nmat
               if (tag_mat(i).eq.mcode) check = 1
            enddo
            !
            if (check.ne.0) then
               iquad = iquad + 1
               con_quad(iquad,1) = mcode
               read(input_line(iright:sl),*)(con_quad(iquad,j),j=2,5)
            endif
        !
         elseif ((ecode.eq.'line').or.(ecode.eq.'LINE')) then
            check = 0
            do i = 1,ndirX
               if (tag_dirX(i).eq.mcode) check = 1
            enddo
            do i = 1,ndirY
               if (tag_dirY(i).eq.mcode) check = 1
            enddo
            do i = 1,nneuX
               if (tag_neuX(i).eq.mcode) check = 1
            enddo
            do i = 1,nneuY
               if (tag_neuY(i).eq.mcode) check = 1
            enddo

			do i = 1,nneuN ! M Mexico Arroyo

               if (tag_neuN(i).eq.mcode) check = 1 ! M Mexico Arroyo

            enddo ! M Mexico Arroyo

            do i = 1,nneuT ! M Mexico Arroyo

               if (tag_neuT(i).eq.mcode) check = 1 ! M Mexico Arroyo

            enddo ! M Mexico Arroyo
            do i = 1,nabc
               if (tag_abc(i).eq.mcode) check = 1
            enddo
            !
            if (check.ne.0) then
               iline = iline +1
               con_line(iline,1) = mcode
               read(input_line(iright:sl),*)(con_line(iline,j),j=2,3)
            endif
         !
         endif
      enddo
        
      
      close(23)
      
      return
      
      end subroutine read_grid_el
      
      
      
      
!     ********************************************************

      subroutine check_file(fichier,err_out)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Javier Sabadell
!     
      character*70 :: fichier
      integer*4 :: err_out
      
      integer*4 :: i,fich_a,fich_b
      logical*4 :: f_ex
      
      do i=1,70
        if (fichier(i:i).ne.' ') exit
      enddo
      
      fich_a=i
      do i=70,1,-1
        if (fichier(i:i).ne.' ') exit
      enddo
      fich_b=i
      
      f_ex=.false. ; err_out=0
      inquire (file=fichier,exist=f_ex)
      if (f_ex) then
        continue
      else
        err_out = 101
        write(*,*)'Input file ',fichier(fich_a:fich_b), & 
                       ' does not exist!'
      endif
      
      return
      
      end subroutine check_file



!     ********************************************************
      
      subroutine write_nodal_output_el(nnode,xs,ys,cs_nnz,cs,nmat,tag_mat,sd,&
                                       file_name,count,time,np)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Luca Massidda, Javier Sabadell
!     
      implicit none
      
      integer*4 :: nnode,cs_nnz,nmat,count,np
      real*8 :: time
      real*8, dimension(nnode) :: xs,ys
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nmat) :: tag_mat
      integer*4, dimension(nmat) :: sd
      character*70 :: file_name
      
      character*70 :: in_fileE,out_file
      character*40 :: subtitle
      integer*4 :: ip,i,in,nd
      integer*4 :: lname
      real*8 :: val
      real*8, dimension(:), allocatable :: ux,uy,ud
      
      allocate (ux(nnode),uy(nnode),ud(nnode))
      
      lname = len_trim(file_name)
      in_fileE = file_name(1:lname) // 'E_xxx_xxx.bin'
      out_file = file_name(1:lname) // '_xx_x.mtv'
      
      if (count.lt.10) then
         write(in_fileE(lname+3:lname+3),'(a1)')'0'
         write(in_fileE(lname+4:lname+4),'(i1)')count
         write(out_file(lname+2:lname+2),'(a1)')'0'
         write(out_file(lname+3:lname+3),'(i1)')count
      else if (count.le.99) then
         write(in_fileE(lname+3:lname+4),'(i2)')count
         write(out_file(lname+2:lname+3),'(i2)')count
      endif
      
      do ip = 0,np-1
         if (ip.lt.10) then
            write(in_fileE(lname+6:lname+6),'(a1)')'0'
            write(in_fileE(lname+7:lname+7),'(i1)')ip
         else if (ip.le.99) then
            write(in_fileE(lname+6:lname+7),'(i2)')ip
         endif
         
         open(20,file=in_fileE)
         
         read(20,*)nd
         
         do i = 1,nd
            read(20,*)in,val
            if (dabs(val).lt.1.d-30) val=0.0d0
            if (in.le.nnode) then
               ux(in) = val
            else
               uy(in -nnode) = val
            endif
         enddo
         
         close(20)
      enddo
      
      do in = 1,nnode
         ud(in) = sqrt(ux(in)**2 + uy(in)**2)
         if (dabs(ud(in)).lt.1.d-30) ud(in) = 0.d0
      enddo
      
      
      write(out_file(lname+5:lname+5),'(a1)')'X'
      subtitle = 'X displacements'
      call write_plotmtv_file(nnode,xs,ys,cs_nnz,cs,nmat,tag_mat,sd,&
                              out_file,subtitle,time,ux)
      
      write(out_file(lname+5:lname+5),'(a1)')'Y'
      subtitle = 'Y displacements'
      call write_plotmtv_file(nnode,xs,ys,cs_nnz,cs,nmat,tag_mat,sd,&
                              out_file,subtitle,time,uy)
      
      write(out_file(lname+5:lname+5),'(a1)')'D'
      subtitle = 'D displacements'
      call write_plotmtv_file(nnode,xs,ys,cs_nnz,cs,nmat,tag_mat,sd,&
                              out_file,subtitle,time,ud)
      
      end subroutine write_nodal_output_el
      
      
      
!     ********************************************************
      
!DA SISTEMARE CON CLARA - 10 Sep 2004 

      subroutine write_element_output_el(nnode,xs,ys,cs_nnz,cs,&
                       nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                       ne,alfa1,alfa2,beta1,beta2,gamma1,gamma2,&
                       file_name,count,time,np)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      implicit none
      
      integer*4 :: nnode,cs_nnz,nm,ne,count,np
      real*8 :: time
      real*8, dimension(nnode) :: xs,ys
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat
      real*8, dimension(nm) :: tref_mat
      real*8, dimension(nm,4) :: prop_mat
      real*8, dimension(ne) :: alfa1,alfa2,beta1,beta2,gamma1,gamma2
      character*70 :: file_name
      
      real*8, dimension(:), allocatable :: ux,uy
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8, dimension(:), allocatable :: dxdx,dxdy,dydx,dydy
      real*8, dimension(:,:), allocatable :: duxdx,duxdy,duydx,duydy
      real*8, dimension(:,:), allocatable :: sxx,sxy,syy,szz,spr,svm
      
      real*8 :: lambda,mu,tref,det_j
      integer*4 :: i,j,nn
      integer*4 :: il,im,ip,iq,ie,imat
      integer*4 :: is,in,is_lq,in_lq,is_pm,in_pm
      integer*4 :: n1,n2,n3,n4
      real*8 :: t1ux,t1uy,t2ux,t2uy
      real*8 :: val,k
      character*70 :: in_fileE
      character*70, dimension(6) :: out_file
      character*40, dimension(6) :: title
      integer,dimension(7) :: iunit
      integer*4 :: nd,lname
      character :: opt
      integer*4 :: tag
      
      
      
      nn = 2
      
      allocate(ct(nn),ww(nn),dd(nn,nn))
      call lgl(nn,ct,ww,dd)
      allocate (dxdx(nn),dydx(nn),dxdy(nn),dydy(nn))
      allocate (duxdx(nn,nn),duxdy(nn,nn),duydx(nn,nn),duydy(nn,nn))
      allocate (sxx(nn,nn),sxy(nn,nn),syy(nn,nn),szz(nn,nn))
      allocate (spr(nn,nn),svm(nn,nn))
      
      allocate (ux(nnode),uy(nnode))
      
      lname = len_trim(file_name)
      in_fileE = file_name(1:lname) // 'E_xxx_xxx.bin'
      
      if (count.lt.10) then
         write(in_fileE(lname+3:lname+3),'(a1)')'0'
         write(in_fileE(lname+4:lname+4),'(i1)')count
      else if (count.le.99) then
         write(in_fileE(lname+3:lname+4),'(i2)')count
      endif
      
      do ip = 0,np-1
         if (ip.lt.10) then
            write(in_fileE(lname+6:lname+6),'(a1)')'0'
            write(in_fileE(lname+7:lname+7),'(i1)')ip
         else if (ip.le.99) then
            write(in_fileE(lname+6:lname+7),'(i2)')ip
         endif
         
         open(20,file=in_fileE)
         
         read(20,*)nd
         
         do i = 1,nd
            read(20,*)in,val
            if (dabs(val).lt.1.d-30) val=0.0d0
            if (in.le.nnode) then
               ux(in) = val
            else
               uy(in -nnode) = val
            endif
         enddo
         
         close(20)
      enddo
      
      out_file(1) = file_name(1:lname) // '_xx_SXX.mtv'
      out_file(2) = file_name(1:lname) // '_xx_SYY.mtv'
      out_file(3) = file_name(1:lname) // '_xx_SXY.mtv'
      out_file(4) = file_name(1:lname) // '_xx_SZZ.mtv'
      out_file(5) = file_name(1:lname) // '_xx_SPR.mtv'
      out_file(6) = file_name(1:lname) // '_xx_SVM.mtv'
      
      title(1) = 'SXX'
      title(2) = 'SYY'
      title(3) = 'SXY'
      title(4) = 'SZZ'
      title(5) = 'SPR'
      title(6) = 'SVM'
      
      if (count.lt.10) then
         do i = 1,6
            write(out_file(i)(lname+2:lname+2),'(a1)')'0'
            write(out_file(i)(lname+3:lname+3),'(i1)')count
         enddo
      else if (count.le.99) then
         do i = 1,6
            write(out_file(i)(lname+2:lname+3),'(i2)')count
         enddo
      endif
      do i = 1,6
         iunit(i) = 30 + i
      enddo
      
      do i = 1,6
         open(iunit(i),file=out_file(i))
      enddo
      
      opt='2'
      tag=0
      
      do i = 1,6
         write(iunit(i),*)'$DATA=CONTCURVE'
         write(iunit(i),*)'%contstyle = '//opt
         write(iunit(i),*)'%toplabel = "',title(i),'"'
         write(iunit(i),299)time
         write(iunit(i),*)'%xlabel= "lateral dimension (m)"'
         write(iunit(i),*)'%ylabel= "depth (m)"'
         write(iunit(i),*)'%nsteps = 50'
         write(iunit(i),*)'%eyepos.x = 0.1d0'
         write(iunit(i),*)'%eyepos.y = 0.5d0'
         write(iunit(i),*)'%eyepos.z = 0.25d0'
         write(iunit(i),*)'%leftworld = FALSE'
         write(iunit(i),*)'%fitpage = FALSE'
         write(iunit(i),*)
      enddo
      
      
      ne = cs(0) -1
      
      do imat = 1,nm
         lambda = prop_mat(imat,2)
         mu = prop_mat(imat,3)
         
         if ((sdeg_mat(imat) +1).ne.nn) then
            deallocate(ct,ww,dd)
            deallocate(dxdx,dydx,dxdy,dydy)
            deallocate(duxdx,duxdy,duydx,duydy)
            deallocate(sxx,sxy,syy,szz)
            deallocate(spr,svm)
            
            nn = sdeg_mat(imat) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call lgl(nn,ct,ww,dd)
            allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn))
            allocate(duxdx(nn,nn),duxdy(nn,nn),duydx(nn,nn),duydy(nn,nn))
            allocate(sxx(nn,nn),sxy(nn,nn),syy(nn,nn),szz(nn,nn))
            allocate(spr(nn,nn),svm(nn,nn))
         endif
         
         do ie = 1,ne
            if (cs(cs(ie -1) +0).eq.tag_mat(imat)) then
               do i = 1,nn
                  dxdy(i) = beta1(ie) + gamma1(ie) * ct(i)
                  dydy(i) = beta2(ie) + gamma2(ie) * ct(i)
               enddo
               
               do j = 1,nn
                  dxdx(j) = alfa1(ie) + gamma1(ie) * ct(j)
                  dydx(j) = alfa2(ie) + gamma2(ie) * ct(j)
               enddo
               
               do iq = 1,nn
                  do ip = 1,nn
                     t1ux = 0.d0; t1uy = 0.d0
                     t2ux = 0.d0; t2uy = 0.d0
                     
                     det_j = dxdx(iq)*dydy(ip) - dxdy(ip)*dydx(iq)
                     
                     do il = 1,nn
                        is_lq = (iq - 1) * nn + il
                        in_lq = cs(cs(ie -1) + is_lq)
                        
                        t1ux = t1ux + ux(in_lq) * dd(ip,il)
                        t1uy = t1uy + uy(in_lq) * dd(ip,il)
                     enddo
                     
                     do im = 1,nn
                        is_pm = (im - 1) * nn + ip
                        in_pm = cs(cs(ie -1) +is_pm)
                        
                        t2ux = t2ux + ux(in_pm) * dd(iq,im)
                        t2uy = t2uy + uy(in_pm) * dd(iq,im)
                     enddo
                     
                     duxdx(ip,iq) = 1.0d0 / det_j &
                          * ((dydy(ip) * t1ux) - (dydx(iq) * t2ux))
                     duydx(ip,iq) = 1.0d0 / det_j &
                          * ((dydy(ip) * t1uy) - (dydx(iq) * t2uy))
                     duxdy(ip,iq) = -1.0d0 / det_j &
                          * ((dxdy(ip) * t1ux) - (dxdx(iq) * t2ux))
                     duydy(ip,iq) = -1.0d0 / det_j &
                          * ((dxdy(ip) * t1uy) - (dxdx(iq) * t2uy))
                  enddo
               enddo
               
               do iq = 1,nn
                  do ip = 1,nn
                     is = (iq - 1) * nn + ip
                     in = cs(cs(ie -1) +is)
                     
                     sxx(ip,iq) = (lambda + 2.0d0 * mu) * duxdx(ip,iq) &
                          + lambda * duydy(ip,iq)
                     syy(ip,iq) = (lambda + 2.0d0 * mu) * duydy(ip,iq) &
                          + lambda * duxdx(ip,iq)
                     sxy(ip,iq) = mu * (duxdy(ip,iq) + duydx(ip,iq))
                     szz(ip,iq) = lambda * (duxdx(ip,iq) + duydy(ip,iq))
                     
                     spr(ip,iq) = (-1.0d0 / 3.0d0) &
                          * (sxx(ip,iq) + syy(ip,iq) + szz(ip,iq))
                     
                     svm(ip,iq) = sxx(ip,iq)*sxx(ip,iq) &
                          + syy(ip,iq)*syy(ip,iq) + szz(ip,iq)*szz(ip,iq) &
                          - syy(ip,iq)*szz(ip,iq) - szz(ip,iq)*sxx(ip,iq) &
                          - sxx(ip,iq)*syy(ip,iq) + 3.0d0*sxy(ip,iq)*sxy(ip,iq)
                     svm(ip,iq) = dsqrt(dabs(svm(ip,iq)))
                  enddo
               enddo
               
               do im = 1,(nn -1)
                  do il =1,(nn -1)
                     n1 = cs(cs(ie -1) +((im - 1) * nn + il))
                     n2 = cs(cs(ie -1) +((im - 1) * nn + il + 1))
                     n3 = cs(cs(ie -1) +(im * nn + il + 1))
                     n4 = cs(cs(ie -1) +(im * nn + il))
                     
                     write(iunit(1),311)xs(n1),ys(n1),sxx(il,im),tag
                     write(iunit(1),311)xs(n2),ys(n2),sxx(il+1,im),tag
                     write(iunit(1),311)xs(n3),ys(n3),sxx(il+1,im+1),tag
                     write(iunit(1),311)xs(n4),ys(n4),sxx(il,im+1),tag
                     write(iunit(1),*)
                     
                     write(iunit(2),311)xs(n1),ys(n1),syy(il,im),tag
                     write(iunit(2),311)xs(n2),ys(n2),syy(il+1,im),tag
                     write(iunit(2),311)xs(n3),ys(n3),syy(il+1,im+1),tag
                     write(iunit(2),311)xs(n4),ys(n4),syy(il,im+1),tag
                     write(iunit(2),*)
                     
                     write(iunit(3),311)xs(n1),ys(n1),sxy(il,im),tag
                     write(iunit(3),311)xs(n2),ys(n2),sxy(il+1,im),tag
                     write(iunit(3),311)xs(n3),ys(n3),sxy(il+1,im+1),tag
                     write(iunit(3),311)xs(n4),ys(n4),sxy(il,im+1),tag
                     write(iunit(3),*)
                     
                     write(iunit(4),311)xs(n1),ys(n1),szz(il,im),tag
                     write(iunit(4),311)xs(n2),ys(n2),szz(il+1,im),tag
                     write(iunit(4),311)xs(n3),ys(n3),szz(il+1,im+1),tag
                     write(iunit(4),311)xs(n4),ys(n4),szz(il,im+1),tag
                     write(iunit(4),*)
                     
                     write(iunit(5),311)xs(n1),ys(n1),spr(il,im),tag
                     write(iunit(5),311)xs(n2),ys(n2),spr(il+1,im),tag
                     write(iunit(5),311)xs(n3),ys(n3),spr(il+1,im+1),tag
                     write(iunit(5),311)xs(n4),ys(n4),spr(il,im+1),tag
                     write(iunit(5),*)
                     
                     write(iunit(6),311)xs(n1),ys(n1),svm(il,im),tag
                     write(iunit(6),311)xs(n2),ys(n2),svm(il+1,im),tag
                     write(iunit(6),311)xs(n3),ys(n3),svm(il+1,im+1),tag
                     write(iunit(6),311)xs(n4),ys(n4),svm(il,im+1),tag
                     write(iunit(6),*)
                  enddo
               enddo
               
            endif
         enddo
      enddo
      
      do i = 1,6
         write(iunit(i),*)'$END'
         close(iunit(i))
      enddo
      
299   format(t2,'%subtitle = "t = ',f6.3,' "')
311   format(5x,e14.6,5x,e14.6,5x,e14.6,5x,i2)
      
      
      deallocate (ct,ww,dd)
      deallocate (ux,uy)
      deallocate (dxdx,dxdy,dydx,dydy)
      deallocate (duxdx,duxdy,duydx,duydy)
      deallocate (sxx,sxy,syy,szz)
      deallocate (spr,svm)
      
      return
      
      end subroutine write_element_output_el
      
      

!     ********************************************************

      subroutine write_plotmtv_file(nnode,xx,yy,cs_nnz,cs,nmat,tag_mat,sd,&
                                    out_file,title,time,vn)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Luca Massidda, Javier Sabadell
!     
      integer*4 :: nnode,cs_nnz,nmat
      character*70 :: out_file
      character*40 :: title
      real*8 :: time
      real*8, dimension(nnode) :: xx,yy,vn
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nmat) :: tag_mat
      integer*4, dimension(nmat) :: sd
      
      character*1 :: opt
      integer*4 :: nn,tag
      integer*4 :: imat,ie,i,j,n1,n2,n3,n4
      
      
      opt='2'

      open(21,file=out_file)
      
      tag=0;
      
      write(21,*)'$DATA=CONTCURVE'
      write(21,*)'%contstyle = '//opt
      write(21,*)'%toplabel = "',title,'"'
      write(21,199)time
      write(21,*)'%xlabel= "lateral dimension (m)"'
      write(21,*)'%ylabel= "depth (m)"'
      write(21,*)'%nsteps = 50'
      write(21,*)'%eyepos.x = 0.1d0'
      write(21,*)'%eyepos.y = 0.5d0'
      write(21,*)'%eyepos.z = 0.25d0'
      write(21,*)'%leftworld = FALSE'
      write(21,*)'%fitpage = FALSE'
      write(21,*)
      
      
      ne = cs(0) -1
      
      do imat = 1,nmat
         nn = sd(imat) +1
         
         do ie = 1,ne
            if (cs(cs(ie -1) +0).eq.tag_mat(imat)) then
               do j = 1,nn -1
                  do i = 1,nn -1
                     n1 = cs(cs(ie -1) +nn*(j -1) +i)
                     n2 = cs(cs(ie -1) +nn*(j -1) +i +1)
                     n3 = cs(cs(ie -1) +nn*j +i +1)
                     n4 = cs(cs(ie -1) +nn*j +i)
                     
                     write(21,211)xx(n1),yy(n1),vn(n1),tag
                     write(21,211)xx(n2),yy(n2),vn(n2),tag
                     write(21,211)xx(n3),yy(n3),vn(n3),tag
                     write(21,211)xx(n4),yy(n4),vn(n4),tag
                     write(21,*)
                  enddo
               enddo
            endif
         enddo
      enddo
      
      
199   format(t2,'%subtitle = "t = ',E12.4,' "')
211   format(5x,e14.6,5x,e14.6,5x,e14.6,5x,i2)
      
      write(21,*)'$END'
      
      close(21)
      
      
      return
      
      end subroutine write_plotmtv_file
      
      
      
!     ********************************************************
      
      subroutine calculate_element_results_el(nnode,xs,ys,cs_nnz,cs,&
                           nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                           ne,alfa1,alfa2,beta1,beta2,gamma1,gamma2,&
                           ux,uy,sxx,syy,szz,sxy,spr,svm,&
						   duxdx_average,duydx_average,duxdy_average,duydy_average,&
						   vistn1_0_read,vistn2_0_read,vistn3_0_read)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      implicit none
      
      integer*4 :: nnode,cs_nnz,nm,ne
      real*8, dimension(nnode) :: xs,ys
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat
      real*8, dimension(nm) :: tref_mat
      real*8, dimension(nm,4) :: prop_mat
      real*8, dimension(ne) :: alfa1,alfa2,beta1,beta2,gamma1,gamma2
      real*8, dimension(nnode) :: ux,uy,ut
      real*8, dimension(nnode) :: sxx,syy,szz,sxy
      real*8, dimension(nnode) :: spr,svm
	  real*8, dimension(nnode) :: vistn1_0_read,vistn2_0_read,vistn3_0_read
      
      real*8 :: lambda,mu,tref
      integer*4 :: nn
      integer*4 :: i,j,il,im,ip,iq,ie,imat
      integer*4 :: is,in,is_lq,in_lq,is_pm,in_pm
      real*8, dimension(:), allocatable :: ct,ww
      real*8, dimension(:,:), allocatable :: dd
      real*8, dimension(:), allocatable :: dxdx,dxdy,dydx,dydy
      real*8, dimension(:,:), allocatable :: duxdx,duxdy,duydx,duydy
      integer*4, dimension(:), allocatable :: nodal_counter
      real*8 :: det_j,t1ux,t1uy,t2ux,t2uy

      real*8, dimension(nnode) :: duxdx_average,duydx_average,duxdy_average,duydy_average
      
      
      nn = 2
      
      allocate(ct(nn),ww(nn),dd(nn,nn))
      call lgl(nn,ct,ww,dd)
      allocate (dxdx(nn),dydx(nn),dxdy(nn),dydy(nn))
      allocate (duxdx(nn,nn),duxdy(nn,nn),duydx(nn,nn),duydy(nn,nn))
      
      allocate(nodal_counter(nnode))
      
      
      do in = 1,nnode
         nodal_counter(in) = 0
         sxx(in) = 0.0d0
         syy(in) = 0.0d0
         szz(in) = 0.0d0
         sxy(in) = 0.0d0
         
         spr(in) = 0.0d0
         svm(in) = 0.0d0


		 duxdx_average(in) = 0.0d0
		 duydx_average(in) = 0.0d0
		 duxdy_average(in) = 0.0d0
		 duydy_average(in) = 0.0d0

      enddo
      
      
      ne = cs(0) -1
      
      do imat = 1,nm
         lambda = prop_mat(imat,2)
         mu = prop_mat(imat,3)
         
         if ((sdeg_mat(imat) +1).ne.nn) then
            deallocate(ct,ww,dd)
            deallocate(dxdx,dydx,dxdy,dydy)
            deallocate(duxdx,duxdy,duydx,duydy)
            
            nn = sdeg_mat(imat) +1
            allocate(ct(nn),ww(nn),dd(nn,nn))
            call lgl(nn,ct,ww,dd)
            allocate(dxdx(nn),dydx(nn),dxdy(nn),dydy(nn))
            allocate(duxdx(nn,nn),duxdy(nn,nn),duydx(nn,nn),duydy(nn,nn))
         endif
         
         do ie = 1,ne
            if (cs(cs(ie -1) +0).eq.tag_mat(imat)) then
               do i = 1,nn
                  dxdy(i) = beta1(ie) + gamma1(ie) * ct(i)
                  dydy(i) = beta2(ie) + gamma2(ie) * ct(i)
               enddo
               
               do j = 1,nn
                  dxdx(j) = alfa1(ie) + gamma1(ie) * ct(j)
                  dydx(j) = alfa2(ie) + gamma2(ie) * ct(j)
               enddo
               
               do iq = 1,nn
                  do ip = 1,nn
                     t1ux = 0.d0; t1uy = 0.d0
                     t2ux = 0.d0; t2uy = 0.d0
                     
                     det_j = dxdx(iq)*dydy(ip) - dxdy(ip)*dydx(iq)
                     
                     do il = 1,nn
                        is_lq = (iq - 1) * nn + il
                        in_lq = cs(cs(ie -1) + is_lq)
                        
                        t1ux = t1ux + ux(in_lq) * dd(ip,il)
                        t1uy = t1uy + uy(in_lq) * dd(ip,il)
                     enddo
                     
                     do im = 1,nn
                        is_pm = (im - 1) * nn + ip
                        in_pm = cs(cs(ie -1) +is_pm)
                        
                        t2ux = t2ux + ux(in_pm) * dd(iq,im)
                        t2uy = t2uy + uy(in_pm) * dd(iq,im)
                     enddo
                     
                     duxdx(ip,iq) = 1.0d0 / det_j &
                          * ((dydy(ip) * t1ux) - (dydx(iq) * t2ux))
                     duydx(ip,iq) = 1.0d0 / det_j &
                          * ((dydy(ip) * t1uy) - (dydx(iq) * t2uy))
                     duxdy(ip,iq) = -1.0d0 / det_j &
                          * ((dxdy(ip) * t1ux) - (dxdx(iq) * t2ux))
                     duydy(ip,iq) = -1.0d0 / det_j &
                          * ((dxdy(ip) * t1uy) - (dxdx(iq) * t2uy))
                  enddo
               enddo
               
               do iq = 1,nn
                  do ip = 1,nn
                     is = (iq - 1) * nn + ip
                     in = cs(cs(ie -1) +is)
                     
                     nodal_counter(in) = nodal_counter(in) +1
                     
                     sxx(in) = sxx(in) + (lambda + 2.0d0*mu) * (duxdx(ip,iq) - vistn1_0_read(in)) &
                          + lambda * (duydy(ip,iq) - vistn2_0_read(in))
                     syy(in) = syy(in) + (lambda + 2.0d0*mu) * (duydy(ip,iq) - vistn2_0_read(in)) &
                          + lambda * (duxdx(ip,iq) - vistn1_0_read(in))
                     szz(in) = szz(in) + lambda &
                          * (duxdx(ip,iq) - vistn1_0_read(in) + duydy(ip,iq) - vistn2_0_read(in))
                     sxy(in) = sxy(in) + mu * (duxdy(ip,iq) + duydx(ip,iq) - 2 * vistn3_0_read(in))

                    
					 duxdx_average(in) = duxdx_average(in) + duxdx(ip,iq)
                     duydx_average(in) = duydx_average(in) + duydx(ip,iq)
					 duxdy_average(in) = duxdy_average(in) + duxdy(ip,iq)
					 duydy_average(in) = duydy_average(in) + duydy(ip,iq)                  
 
                  enddo
               enddo
               
            endif
         enddo
      enddo
      
      
      do in = 1,nnode
         sxx(in) = sxx(in) / nodal_counter(in)
         syy(in) = syy(in) / nodal_counter(in)
         szz(in) = szz(in) / nodal_counter(in)
         sxy(in) = sxy(in) / nodal_counter(in)
         
         spr(in) = (-1.0d0 / 3.0d0) * (sxx(in) + syy(in) + szz(in))
         svm(in) = sxx(in)*sxx(in) + syy(in)*syy(in) + szz(in)*szz(in) &
                 - syy(in)*szz(in) - szz(in)*sxx(in) - sxx(in)*syy(in) &
                 + 3.0d0*sxy(in)*sxy(in)
         svm(in) = dsqrt(dabs(svm(in)))

         duxdx_average(in) = duxdx_average(in) / nodal_counter(in)
		 duydx_average(in) = duydx_average(in) / nodal_counter(in)
		 duxdy_average(in) = duxdy_average(in) / nodal_counter(in)
		 duydy_average(in) = duydy_average(in) / nodal_counter(in)

      enddo
      
      deallocate(ct,ww,dd)
      deallocate(dxdx,dydx,dxdy,dydy)
      deallocate(duxdx,duxdy,duydx,duydy)
      deallocate(nodal_counter)
      
      return
      
      end subroutine calculate_element_results_el
      
      
      
!!     ********************************************************
!      
!      subroutine write_UCD_output_el(nnode,xs,ys,cs_nnz,cs,&
!                           nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
!                           ne,alfa1,alfa2,beta1,beta2,gamma1,gamma2,&
!                           opt,nsnaps,tsnaps,file_name,np)
!      
!!     © CRS4, 2002, All Rights Reserved
!!     Authors: Luca Massidda
!!     
!      implicit none
!      
!      integer*4 :: nnode,cs_nnz,nm,ne,opt,nsnaps,np
!      real*8, dimension(nnode) :: xs,ys
!      integer*4, dimension(0:cs_nnz) :: cs
!      integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat
!      real*8, dimension(nm) :: tref_mat
!      real*8, dimension(nm,4) :: prop_mat
!      real*8, dimension(ne) :: alfa1,alfa2,beta1,beta2,gamma1,gamma2
!      real*8, dimension(nsnaps) :: tsnaps
!      character*70 :: file_name
!      
!      real*8, dimension(:), allocatable :: ux,uy
!      real*8, dimension(:), allocatable :: sxx,syy,szz,sxy
!      real*8, dimension(:), allocatable :: spr,svm
!      integer*4 :: ip,i,j,in,nd,nn,nelem
!      integer*4 :: imat,ie,mn,n1,n2,n3,n4,count,ic
!      real*8 :: val
!      integer*4 :: istep,nsteps      
!      character*70 :: in_fileE,out_file
!      character*40 :: elcode
!      integer*4 :: lname
!
!      real*8, dimension(:), allocatable :: duxdx_average,duydx_average,duxdy_average,duydy_average
!      
!      
!      lname = len_trim(file_name)
!      in_fileE = file_name(1:lname) // 'E_xxx_xxx.bin'
!      out_file = file_name(1:lname) // '.inp'
!      
!      istep = 1
!
!      allocate(ux(nnode),uy(nnode))
!      allocate(sxx(nnode),syy(nnode),szz(nnode),sxy(nnode))
!      allocate(spr(nnode),svm(nnode))
!	  allocate(duxdx_average(nnode),duydx_average(nnode),duxdy_average(nnode),duydy_average(nnode))
!      
!      if (opt.eq.1) then
!         nsteps = nsnaps
!      elseif (opt.eq.2) then
!         nsteps = nsnaps*2
!      elseif (opt.eq.3) then
!         nsteps = nsnaps*3
!      endif
!      
!      
!      open(21,file=out_file)
!      
!      write(21,'(A)')'# Elastic else output'
!      write(21,'(A)')'# CRS4'
!      write(21,'(A)')'# Center for Advanced Studies, '
!      write(21,'(A)')'# Research and Development in Sardinia'
!      write(21,'(A)')'#'
!      write(21,'(I2)')nsteps
!      write(21,'(A)')'data'
!      
!      do count = 1,nsnaps
!         if (istep.lt.10) then
!            write(21,'(A4,I1,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!         else
!            write(21,'(A4,I2,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!         endif
!         
!         if (istep.eq.1) then
!            ne = cs(0) -1
!            
!            nelem = 0
!            do imat = 1,nm
!               nn = sdeg_mat(imat) +1
!               do ie = 1,ne
!                  if (cs(cs(ie -1) +0).eq.tag_mat(imat)) then
!                     nelem = nelem + (nn -1)*(nn -1)
!                  endif
!               enddo
!            enddo
!            
!            write(21,'(I8,2X,I8)')nnode,nelem
!            
!            do i = 1,nnode
!               write(21,'(I8,2X,E14.6,2X,E14.6,2X,E14.6)')i,xs(i),ys(i),0.0d0
!            enddo
!            
!            ic = 0
!            elcode = 'quad'
!            
!            do imat = 1,nm
!               nn = sdeg_mat(imat) +1
!               
!               do ie = 1,ne
!                  if (cs(cs(ie -1) +0).eq.tag_mat(imat)) then
!                     mn = cs(cs(ie -1) +0)
!                     
!                     do j = 1,nn -1
!                        do i = 1,nn -1
!                           n1 = cs(cs(ie -1) +nn*(j -1) +i)
!                           n2 = cs(cs(ie -1) +nn*(j -1) +i +1)
!                           n3 = cs(cs(ie -1) +nn*j +i +1)
!                           n4 = cs(cs(ie -1) +nn*j +i)
!                           
!                           ic = ic +1
!                           write(21,'(I8,2X,I4,2X,A6,2X,I8,2X,I8,2X,I8,2X,I8)')&
!                                ic,mn,elcode,n1,n2,n3,n4
!                        enddo
!                     enddo
!                  endif
!               enddo
!            enddo
!            
!         endif
!         
!         
!         if (count.lt.10) then
!            write(in_fileE(lname+3:lname+4),'(a2)')'00'
!            write(in_fileE(lname+5:lname+5),'(i1)')count
!	     else if (count.lt.100) then
!            write(in_fileE(lname+3:lname+3),'(a1)')'0'
!            write(in_fileE(lname+4:lname+5),'(i2)')count
!         else if (count.le.999) then
!            write(in_fileE(lname+3:lname+5),'(i3)')count
!         endif
!         
!         do ip = 0,np -1
!            if (ip.lt.10) then
!               write(in_fileE(lname+7:lname+8),'(a2)')'00'
!               write(in_fileE(lname+9:lname+9),'(i1)')ip
!			else if (ip.lt.100) then
!               write(in_fileE(lname+7:lname+7),'(a1)')'0'
!               write(in_fileE(lname+8:lname+9),'(i2)')ip
!            else if (ip.le.999) then
!               write(in_fileE(lname+7:lname+9),'(i3)')ip
!            endif
!            
!            open(20,file=in_fileE)
!            
!            read(20,*)nd
!            
!            do i = 1,nd
!               read(20,*)in,val
!               if (dabs(val).lt.1.d-30) val=0.d0
!               if (in.le.nnode) then
!                  ux(in) = val
!               else
!                  uy(in -nnode) = val
!               endif
!            enddo
!            
!            close(20)
!
!
!
!         enddo
!         
!         
!         write(21,*)'3 0 '
!         write(21,*)'1 3 '
!         write(21,*)'Displacement, NO_UNITS'
!         
!         do i = 1,nnode
!            write(21,'(I8,2X,E14.6,2X,E14.6,2X,E14.6)')&
!                 i,ux(i),uy(i),0.0d0
!         enddo
!         
!         istep = istep + 1
!         
!         if ((opt.eq.2).or.(opt.eq.3)) then 
!            call calculate_element_results_el(nnode,xs,ys,cs_nnz,cs,&
!                           nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
!                           ne,alfa1,alfa2,beta1,beta2,gamma1,gamma2,&
!                           ux,uy,sxx,syy,szz,sxy,spr,svm,&
!						   duxdx_average,duydx_average,duxdy_average,duydy_average)
!            
!            if (istep.lt.10) then
!               write(21,'(A4,I1,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!            else
!               write(21,'(A4,I2,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!            endif
!            
!            write(21,*)'4 0 '
!            write(21,*)'4 1 1 1 1 '
!            write(21,*)'SXX, NO_UNITS'
!            write(21,*)'SYY, NO_UNITS'
!            write(21,*)'SZZ, NO_UNITS'
!            write(21,*)'SXY, NO_UNITS'
!            
!            do i = 1,nnode
!               write(21,'(I8,2X,E14.6,2X,E14.6,2X,E14.6,2X,E14.6)')&
!                    i,sxx(i),syy(i),szz(i),sxy(i)
!            enddo
!            
!            istep = istep + 1
!         endif
!         
!         if (opt.eq.3) then
!            if (istep.lt.10) then
!               write(21,'(A4,I1,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!            else
!               write(21,'(A4,I2,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!            endif
!            
!            write(21,*)'2 0 '
!            write(21,*)'2 1 1 '
!            write(21,*)'Pressure, NO_UNITS'
!            write(21,*)'VonMises, NO_UNITS'
!            
!            do i = 1,nnode
!               write(21,'(I8,2X,E14.6,2X,E14.6)')&
!                    i,spr(i),svm(i)
!            enddo
!            
!            istep = istep + 1
!         endif
!         
!      enddo
!      
!      close(21)
!      
!      deallocate(ux,uy)
!      deallocate(sxx,syy,szz,sxy)
!      deallocate(spr,svm)
!      deallocate(duxdx_average,duydx_average,duxdy_average,duydy_average)
!
!      return
!      
!      end subroutine write_UCD_output_el
!      
!      
!      
!!     ********************************************************
      
      subroutine write_bin_file(file_name,count,proc,nu,ui,uv)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     
      integer*4 :: count,proc,nu
      real*8, dimension(*) :: uv
      integer*4, dimension(*) :: ui
      character*70 :: file_name

      character*70 :: out_file
      integer*4 :: i,lname
      
      lname = len_trim(file_name)
      out_file = file_name(1:lname) // '_xxx_xxx.bin'
      
      if (count.lt.10) then
         write(out_file(lname+2:lname+3),'(a2)')'00'
         write(out_file(lname+4:lname+4),'(i1)')count
      else if (count.lt.100) then
		 write(out_file(lname+2:lname+2),'(a1)')'0'
         write(out_file(lname+3:lname+4),'(i2)')count
	  else if (count.le.999) then
         write(out_file(lname+2:lname+4),'(i3)')count
      endif
      
      if (proc.lt.10) then
         write(out_file(lname+6:lname+7),'(a2)')'00'
         write(out_file(lname+8:lname+8),'(i1)')proc
      else if (proc.lt.100) then
		 write(out_file(lname+6:lname+6),'(a1)')'0'
         write(out_file(lname+7:lname+8),'(i2)')proc
      else if (proc.le.999) then
         write(out_file(lname+6:lname+8),'(i3)')proc
      endif
      
      open(20,file=out_file)
      
      write(20,*)nu
      
      do i = 1,nu
         write(20,*)ui(i),uv(i)
      enddo
      
      close(20)
      
      return
      
      end subroutine write_bin_file
      
      
      
      
!     ********************************************************
      
    subroutine make_internal_force(lambda0,mu0,rho0,&

        nmatg,valmatg,which_matg,tagmatg,typematg,& ! Mexico Paco grad_lin 29.11.2004
        nn,ct,ww,dd,&
        dxdx,dxdy,dydx,dydy,ux,uy, &
        duxdx,duxdy,duydx,duydy,sxx,syy,sxy,&
        fx,fy,&
        check_node_sism,check_dist_node_sism,&
        length_cns,ielem,facsmom,nl_sism,& 
        func_type,func_indx,func_data,nfunc_data,nf,tt2,&
        nvpcl,fvpcl,valvpcl,which_vpcl,& ! Clara
        nvpsa,valvpsa,which_vpsa,& ! MCNAP
        nvpsl,valvpsl,which_vpsl,& ! MCNAP
        nvpsd,valvpsd,which_vpsd,& ! MCNAP 
        nmaps,fmaps,valmaps,which_maps,& ! Marco 
        cs_nnz,cs,vivel_length,vivel_el,dt,nelem,& ! Clara 
        vistn1_0,vistn2_0,vistn3_0,& ! Clara 
        tagvpcl,tagvpsa,tagvpsl,tagvpsd,tagmat,& ! Clara & MCNAP
        tagmaps,& ! Marco 
        xs,ys,nnt,tempo,&
        sxx_ps,syy_ps,sxy_ps,szz_ps,its,& ! Clara
        a1n1,a1n2,a1n3,a1n4,a1n5,a1n6,&  !MCNAP vpsa
        epn1,epn2,epn3,epn4,epn5,epn6,&  !MCNAP vpsa
        z14n,u7n,densit)                 !MCNAP vpsa
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     Modified by: Marco Stupazzini (27/09/2002)
!     Modification: introduction of moment tensor
!     
      implicit none
      
      real*8 :: lambda0,mu0,rho0
	  real*8 :: lambda,mu
	  integer*4 :: imatg,nmatg,which_matg ! Mexico Paco grad_lin 29.11.2004
	  integer*4, dimension(*) :: tagmatg ! Mexico Paco grad_lin 29.11.2004
	  integer*4, dimension(*) :: typematg ! Mexico Paco grad_lin 29.11.2004
      real*8, dimension(nmatg,18) :: valmatg ! Mexico Paco grad_lin 29.11.2004
      integer*4 :: nn,nelem
      real*8, dimension(nn) :: ct,ww
      real*8, dimension(nn,nn) :: dd
      real*8, dimension(nn) :: dxdx,dxdy,dydx,dydy
      real*8, dimension(nn,nn) :: ux,uy
      real*8, dimension(nn,nn) :: duxdx,duxdy,duydx,duydy
      real*8, dimension(nn,nn) :: sxx,syy,sxy,fx,fy

	  real*8, dimension(nn,nn,nelem) :: sxx_ps,syy_ps,sxy_ps,szz_ps ! Marco
      
      integer*4 :: ip,iq,il,im,ielem,i
      real*8 :: det_j,t1ux,t1uy,t2ux,t2uy,t1fx,t1fy,t2fx,t2fy
      character*12 :: name_prop
      
      integer*4 :: length_cns,nl_sism,nf
	  integer*4, dimension(length_cns,5) :: check_node_sism
      real*8, dimension(length_cns,1) :: check_dist_node_sism
      real*8, dimension(nl_sism,3) :: facsmom
      real*8 :: tt2
      integer*4, dimension(nf) :: func_type
      integer*4, dimension(nf +1) :: func_indx



	  integer*4 :: nfunc_data
      real*8, dimension(nfunc_data) :: func_data
      real*8 :: get_func_value

	  integer*4 :: nvpcl,which_vpcl ! Clara 
      real*8, dimension(nvpcl,*) :: valvpcl  ! Clara  
      integer*4, dimension(*) :: fvpcl ! Clara



	  integer*4 :: nvpsa,which_vpsa ! MCNAP 

      real*8, dimension(nvpsa,*) :: valvpsa  ! MCNAP  

      !integer*4, dimension(*) :: fvpsa ! MCNAP



	  integer*4 :: nvpsl,which_vpsl ! MCNAP 

      real*8, dimension(nvpsl,*) :: valvpsl  ! MCNAP  

      !integer*4, dimension(*) :: fvpsl ! MCNAP



	  integer*4 :: nvpsd,which_vpsd ! MCNAP :-)

      real*8, dimension(nvpsd,*) :: valvpsd  ! MCNAP  

      !integer*4, dimension(*) :: fvpsd ! MCNAP 



	  integer*4 :: nmaps,which_maps ! Marco 

      real*8, dimension(nmaps,*) :: valmaps  ! Marco  

      integer*4, dimension(*) :: fmaps ! Marco 

	  !real*8 :: fdatm, hards, gamma, ddelta, frict ! Clara 
	  real*8 :: devia1,devia2,devia3 ! Clara 
      real*8 :: avect1,avect2,avect3 ! Clara 
      real*8 :: vivel1,vivel2,vivel3 ! Clara 
	  real*8 :: vistn1,vistn2,vistn3 ! Clara 
	  real*8 :: steff,theta, varj2, yield ! Clara 

      integer*4 :: vivel_length ! Clara 
      integer*4 :: cs_nnz,is,in
	  integer*4, dimension(0:cs_nnz) :: cs

	  real*8, dimension(nn,nn,nelem) :: vivel_el ! Clara 
	  real*8, dimension(nn,nn,nelem) :: vistn1_0 ! Clara 
	  real*8, dimension(nn,nn,nelem) :: vistn2_0 ! Clara 
	  real*8, dimension(nn,nn,nelem) :: vistn3_0 ! Clara 



      real*8, dimension(nn,nn,nelem) :: a1n1,a1n2,a1n3,a1n4,a1n5,a1n6 !MCNAP vpsa

      real*8, dimension(nn,nn,nelem) :: epn1,epn2,epn3,epn4,epn5,epn6 !MCNAP vpsa

	  real*8, dimension(nn,nn,nelem) :: z14n,u7n,densit               !MCNAP vpsa



      real*8, dimension(6)  :: s,ep,dep,a1  !MCNAP vpsa

	  real*8, dimension(14) :: u            !MCNAP vpsa

      real*8, dimension(3,3):: a,ak,D,DV    !MCNAP vpsa

      real*8 :: al,csi,z14,bef,pc,dens  !MCNAP vpsa

	  integer*4 :: j,k                      !MCNAP vpsa

	
	  real*8 :: destn1,destn2,destn3,debar,dt ! Clara 

	  integer*4 :: nnt,iaz ! Clara 

	  integer*4 :: ivpcl,tagmat ! Clara

	  integer*4 :: ivpsa ! MCNAP

	  integer*4 :: ivpsl ! MCNAP

	  integer*4 :: ivpsd ! MCNAP

	  integer*4 :: imaps ! Marco
	  integer*4, dimension(*) :: tagvpcl ! Clara

	  integer*4, dimension(*) :: tagvpsa ! MCNAP

	  integer*4, dimension(*) :: tagvpsl ! MCNAP

	  integer*4, dimension(*) :: tagvpsd ! MCNAP

	  integer*4, dimension(*) :: tagmaps ! Marco 
	  real*8, dimension(nnt) :: xs,ys ! Clara 
	  integer*4 :: go,tempo ! Clara 

	  integer*4 :: its ! Marco



	  !real*8 :: rho_alpha2,rho_alpha_beta,rho_beta2 ! M Mexico Paco

	  !real*8 :: rho_gamma2,rho_gamma_delta,rho_delta2 ! M Mexico Paco

	  !integer*4 :: paco ! M Mexico Paco



	  !paco = 0 ! M Mexico Paco


      
!   DERIVATIVE CALCULATION
      
      do iq = 1,nn
         do ip = 1,nn
            t1ux = 0.d0; t1uy = 0.d0
            t2ux = 0.d0; t2uy = 0.d0
            
            det_j = dxdx(iq)*dydy(ip) - dxdy(ip)*dydx(iq)
            
            do il = 1,nn
               t1ux = t1ux + ux(il,iq) * dd(ip,il)
               t1uy = t1uy + uy(il,iq) * dd(ip,il)
            enddo
            
            do im = 1,nn
               t2ux = t2ux + ux(ip,im) * dd(iq,im)
               t2uy = t2uy + uy(ip,im) * dd(iq,im)
            enddo
            
            duxdx(ip,iq) = 1.0d0 / det_j &
                 * ((dydy(ip) * t1ux) - (dydx(iq) * t2ux))
            duydx(ip,iq) = 1.0d0 / det_j &
                 * ((dydy(ip) * t1uy) - (dydx(iq) * t2uy))
            duxdy(ip,iq) = -1.0d0 / det_j &
                 * ((dxdy(ip) * t1ux) - (dxdx(iq) * t2ux))
            duydy(ip,iq) = -1.0d0 / det_j &
                 * ((dxdy(ip) * t1uy) - (dxdx(iq) * t2uy))
         enddo
      enddo
      
      
!   STRESS CALCULATION
      
                 
      do iq = 1,nn   
         do ip = 1,nn


!! LINEAR VARIATION OF THE INTERNAL SOIL MECHANICS PROPERTIES:

!! VS & VP INCREASE LINEARLY INSIDE THE ELMENT

!! IMPORTANT: see validation tests into laptop Marco:

!! E:\Users\stupa\Crs4\Selse2d\point_load\sesma\

!! LOOK AT "mate.xls" in order to understand parameters:

!! rho_alpha2, rho_alpha_beta, etc. etc.



if (which_matg.ne.0) then ! Mexico Paco grad_lin 29.11.2004

    do imatg = 1,nmatg ! Mexico Paco grad_lin 29.11.2004

        if (tagmat.eq.tagmatg(imatg)) then ! Mexico Paco grad_lin 29.11.2004


					  is = nn*(iq -1) +ip ! Mexico Paco grad_lin 29.11.2004

					  in = cs(cs(ielem -1) + is) ! Mexico Paco grad_lin 29.11.2004



					  ! LINEAR GRADIENT

					  if (typematg(imatg).eq.1) then



						mu = valmatg(imatg,13) * ys(in)**2 + & ! Mexico Paco grad_lin 29.11.2004

							 valmatg(imatg,14) * ys(in) + & ! Mexico Paco grad_lin 29.11.2004

					  		 valmatg(imatg,15) ! Mexico Paco grad_lin 29.11.2004



						lambda = valmatg(imatg,16) * ys(in)**2 + & ! Mexico Paco grad_lin 29.11.2004

								 valmatg(imatg,17) * ys(in) + & ! Mexico Paco grad_lin 29.11.2004

				      			 valmatg(imatg,18) - 2*mu ! Mexico Paco grad_lin 29.11.2004



					  ! SIGMOIDAL GRADIENT

					  elseif (typematg(imatg).eq.2) then



					    mu = ( valmatg(imatg,13) * &

						     (1/(1+exp(-valmatg(imatg,17)*(ys(in) + valmatg(imatg,18)))))  + &

						       valmatg(imatg,14) )**2 * rho0



					    lambda = ( valmatg(imatg,15) * &

						         (1/(1+exp(-valmatg(imatg,17)*(ys(in) + valmatg(imatg,18)))))  + &

						          valmatg(imatg,16) )**2 * rho0 - 2 * mu





					  endif







					endif ! Mexico Paco grad_lin 29.11.2004

				enddo ! Mexico Paco grad_lin 29.11.2004

		   	else ! Mexico Paco grad_lin 29.11.2004

			   

			   mu = mu0 ! Mexico Paco grad_lin 29.11.2004

			   lambda = lambda0 ! Mexico Paco grad_lin 29.11.2004



			endif ! Mexico Paco grad_lin 29.11.2004




            sxx(ip,iq) = (lambda +2.0d0*mu) * (duxdx(ip,iq) - vistn1_0(ip,iq,ielem)) & ! Clara 
                         + lambda * (duydy(ip,iq) - vistn2_0(ip,iq,ielem)) ! Clara 
            syy(ip,iq) = (lambda +2.0d0*mu) * (duydy(ip,iq) - vistn2_0(ip,iq,ielem)) & ! Clara 
                         + lambda * (duxdx(ip,iq) - vistn1_0(ip,iq,ielem)) ! Clara 
            sxy(ip,iq) = mu * (duxdy(ip,iq) + duydx(ip,iq) - 2*vistn3_0(ip,iq,ielem)) ! Clara

			

			!if (its.eq.3) then

				sxx(ip,iq) = sxx(ip,iq) - sxx_ps(ip,iq,ielem) ! Marco	

				syy(ip,iq) = syy(ip,iq) - syy_ps(ip,iq,ielem) ! Marco

				sxy(ip,iq) = sxy(ip,iq) - sxy_ps(ip,iq,ielem) ! Marco

            !    szz(ip,iq) = szz(ip,iq) - szz_ps(ip,iq,ielem) ! Clara 12/01/2005

			!endif

			!sxx(ip,iq) = sxx(ip,iq) - ((lambda +2.0d0*mu) * vistn1_0(ip,iq,ielem) &
            !             + lambda * vistn2_0(ip,iq,ielem)) 
            !syy(ip,iq) = syy(ip,iq) - ((lambda +2.0d0*mu) * vistn2_0(ip,iq,ielem) &
            !             + lambda * vistn1_0(ip,iq,ielem)) 
            !sxy(ip,iq) = sxy(ip,iq) - (mu * (vistn3_0(ip,iq,ielem) + vistn3_0(ip,iq,ielem)))


         enddo
      enddo



	  !write(*,*)sqrt((lambda+2*mu)/rho0)


! **********************

! ***** PRE-STRESS ***** 

!

!	if (which_maps.ne.0) then ! Marco

!

!		do iq = 1,nn   ! Marco

!			do ip = 1,nn ! Marco

!                is = nn*(iq -1) +ip ! Marco

!                in = cs(cs(ielem -1) + is) ! Marco

!				syy(ip,iq) = syy(ip,iq) + 2000 * 9.81 * abs(ys(in)) ! Marco

!				sxx(ip,iq) = sxx(ip,iq) + (1-sin(valmaps(which_maps,1)*0.017453292)) * 2000 * 9.81 * abs(ys(in)) ! Marco

!			enddo ! Marco

!		enddo ! Marco

!

!	endif ! Marco

!

! ***** PRE-STRESS *****

! **********************






!!   STRESS CALCULATION!
!	
!     ! Clara - Geostatico begin
!	 ! NON FUNZIONA!!!
!
!	go = 1
!
!	!if (tempo.eq.1) then
!	if (go.eq.1) then
!		
!		do iq = 1,nn   
!			do ip = 1,nn
!
!               is = nn*(iq -1) +ip
!                in = cs(cs(ielem -1) + is)
!
!				sxx(ip,iq) = sxx(ip,iq) + (1-sin(27*0.017453292)) * 2000 * 9.81 * abs(ys(in))
!				syy(ip,iq) = syy(ip,iq) + 2000 * 9.81 * abs(ys(in))
!				!sxy(ip,iq) = sxy(ip,iq) - 1500 * 9.81 * xs(in)
!				!write(54,*) iq, ip,in, ys(in),sxx(ip,iq),syy(ip,iq),(1-sin(27*0.017453292)) * 2000 * 9.81 * abs(ys(in)),2000 * 9.81 * abs(ys(in))
!
!			enddo
!		enddo
!
!	endif
!
!    ! Clara - Geostatico end
!
!	!endif
      


!   STRESS CALCULATION - SEISMIC MOMENT
!   OCCHIO AGLI STRESS
      
      if (nl_sism.gt.0) then
         if ((ielem.ge.check_node_sism(1,1)).and. &
            (ielem.le.check_node_sism(length_cns,1))) then

            do i = 1,length_cns
               if (ielem.eq.check_node_sism(i,1)) then              
               

                   do iq = 1,nn
                      do ip = 1,nn
                         if ((check_node_sism(i,3).eq.iq).and. &
                            (check_node_sism(i,2).eq.ip)) then
                            sxx(ip,iq) = sxx(ip,iq) &
                                         - get_func_value(nf,func_type,func_indx,func_data,nfunc_data, &
                                         check_node_sism(i,5),tt2,check_dist_node_sism(i,1)) &
                                      * facsmom(check_node_sism(i,4),1)
                            syy(ip,iq) = syy(ip,iq) &
                                         - get_func_value(nf,func_type,func_indx,func_data,nfunc_data, &
                                         check_node_sism(i,5),tt2,check_dist_node_sism(i,1)) &
                                         * facsmom(check_node_sism(i,4),2)
                            sxy(ip,iq) = sxy(ip,iq) &
                                         - get_func_value(nf,func_type,func_indx,func_data,nfunc_data, &
                                         check_node_sism(i,5),tt2,check_dist_node_sism(i,1)) &
                                         * facsmom(check_node_sism(i,4),3)

                         endif
                      enddo
                   enddo

               endif
             enddo
          endif  
       endif !if (nl_sism.gt.0) then


       !   STRESS CALCULATION - MATERIAL VISCO-PLASTIC

       ! Clara - begin

       !               ************
       ! VISCO PLASTIC *** CLAY ***

       !               ************

        if (which_vpcl.ne.0) then
        ! ****************************
            do ivpcl = 1,nvpcl
                if (tagmat.eq.tagvpcl(ivpcl)) then

                !fdatm     = valvpcl(which_vpcl,1)   ! Yield stress         
                !hards     = valvpcl(which_vpcl,2)   ! Softening modulus      
                !gamma     = valvpcl(which_vpcl,3)   ! First parameter of Perzyn
                !ddelta    = valvpcl(which_vpcl,4)   ! Second parameter -------
                !frict     = valvpcl(which_vpcl,5)   ! Friction angle in grades of plastic yield      !!Clara_non_ass
                !frictpot  = valvpcl(which_vpcl,6)   ! Friction angle in grades of plastic potential  !!Clara_non_ass
                    do iq = 1,nn   
                        do ip = 1,nn
                        !vivel1 = 0.0d0
                        !vivel2 = 0.0d0
                        !vivel3 = 0.0d0

                            call invarg (devia1,devia2,devia3, &
                                valvpcl(which_vpcl,1),valvpcl(which_vpcl,2),valvpcl(which_vpcl,3), &
                                valvpcl(which_vpcl,4),valvpcl(which_vpcl,5), &
                                steff, sxx(ip,iq),syy(ip,iq),sxy(ip,iq), &
                                theta, varj2, yield, fvpcl(which_vpcl))

                            call yieldf (avect1,avect2,avect3,devia1,devia2,devia3, &
                                valvpcl(which_vpcl,1),valvpcl(which_vpcl,2),valvpcl(which_vpcl,3), &
                                valvpcl(which_vpcl,4),valvpcl(which_vpcl,6), &  !!Clara_non_ass
                                steff, theta,varj2,fvpcl(which_vpcl))
                            call flowvp (valvpcl(which_vpcl,1),valvpcl(which_vpcl,2),valvpcl(which_vpcl,3), &
                                valvpcl(which_vpcl,4),valvpcl(which_vpcl,5), &
                                devia1,devia2,devia3, avect1,avect2,avect3, vivel1,vivel2,vivel3, &
                                vivel_el(ip,iq,ielem),yield,vivel_length,fvpcl(which_vpcl))

                            ! Incremental viscoplastic strain rate and the total viscopl.strain

                            destn1 = vivel1 * dt
                            destn2 = vivel2 * dt
                            destn3 = vivel3 * dt

            !vistn1 = vistn1 + destn1
            !vistn2 = vistn2 + destn2
            !vistn3 = vistn3 + destn3


            ! Accumulate the absolute value of the viscoplastic strain increment


            !write(54,*)ip, iq, ielem, debar

            debar=dsqrt((2.0*(destn1**2 + destn2**2 + destn3**2)/3.0))

            vivel_el(ip,iq,ielem) =  vivel_el(ip,iq,ielem) + debar


            vistn1_0(ip,iq,ielem) = vistn1_0(ip,iq,ielem) + destn1 
            vistn2_0(ip,iq,ielem) = vistn2_0(ip,iq,ielem) + destn2 
            vistn3_0(ip,iq,ielem) = vistn3_0(ip,iq,ielem) + destn3 

            enddo !ip
            enddo !iq


        endif
        enddo !ivpcl

        ! ****************************

    endif ! which_vpcl

    ! Clara - end







	 ! MCNAPvpsa   - begin

	 !

	 !               ******************************************************

	 ! VISCO PLASTIC *** di Prisco Viscoplastic model for SAND material ***

	 !               ******************************************************



	 if (which_vpsa.ne.0) then



	 ! ****************************



	 do ivpsa = 1,nvpsa

	    if (tagmat.eq.tagvpsa(ivpsa)) then

               

		do ivpsl = 1,nvpsl

			if (tagmat.eq.tagvpsl(ivpsl)) then



			do ivpsd = 1,nvpsd

				if (tagmat.eq.tagvpsd(ivpsd)) then

	

				do iq = 1,nn   

					do ip = 1,nn

                       

!                     conversion of stress or strain vectors according to conventions:

!					   s = stress/strain vector defined according to the following

!                              correspondence rule:

!                         2-d case       

!                         (11) <=> 1    

!                         (22) <=> 2    

!                         (33) <=> 3    

!                         (12) <=> 4   

!                        and soil mechanics sign conventions (compression positive)

!                      sxx = stress/strain vector defined according to GEO-ELSE

!                           correspondence rule:

!                         2-d case       

!                         (11) <=> 1    

!                         (22) <=> 2    

!                         (33) <=> 4    

!                         (12) <=> 3   

!                        and continuum mechanics sign conventions (tractions positive)

!



				    	s(1)=-sxx(ip,iq)

                        s(2)=-syy(ip,iq)

                        s(3)=0.0d0     !controllare se è zero

                        s(4)=-sxy(ip,iq)

	                    s(5)=0.0d0

	                    s(6)=0.0d0



						a1(1)=a1n1(ip,iq,ielem)

                        a1(2)=a1n2(ip,iq,ielem)

						a1(3)=a1n3(ip,iq,ielem)

						a1(4)=a1n4(ip,iq,ielem)

						a1(5)=a1n5(ip,iq,ielem)

						a1(6)=a1n6(ip,iq,ielem)



						ep(1)=epn1(ip,iq,ielem)

						ep(2)=epn2(ip,iq,ielem)

                        ep(3)=epn3(ip,iq,ielem)

                        ep(4)=epn4(ip,iq,ielem)

						ep(5)=epn5(ip,iq,ielem)

                        ep(6)=epn6(ip,iq,ielem)

	

                        u(7)=u7n(ip,iq,ielem)  

						z14=z14n(ip,iq,ielem)

						dens=densit(ip,iq,ielem)





                        call PAR(dens,u,& 

	                             valvpsl(which_vpsl,1),valvpsl(which_vpsl,2),& 

		                         valvpsl(which_vpsl,3),valvpsl(which_vpsl,4),&

					             valvpsl(which_vpsl,5), valvpsl(which_vpsl,6),& 

		                         valvpsl(which_vpsl,7),valvpsl(which_vpsl,8),&

		                         valvpsl(which_vpsl,9), valvpsl(which_vpsl,10),& 

		                         valvpsl(which_vpsl,11),valvpsl(which_vpsl,12),&

		                         valvpsl(which_vpsl,13),& 

		                         valvpsd(which_vpsd,1),valvpsd(which_vpsd,2),& 

	                          	 valvpsd(which_vpsd,3),valvpsd(which_vpsd,4),&

		                         valvpsd(which_vpsd,5),valvpsd(which_vpsd,6),& 

		                         valvpsd(which_vpsd,7),valvpsd(which_vpsd,8),&

	                           	 valvpsd(which_vpsd,9),valvpsd(which_vpsd,10),& 

		                         valvpsd(which_vpsd,11),valvpsd(which_vpsd,12),&

		                         valvpsd(which_vpsd,13))

       

       



                       bef=u(11)+(u(3)-u(11))*exp(-z14*u(5))

			           pc=u(7)



                       call ten(a1,ak)

                       call scal(ak,ak,al)

                       al=sqrt(al)



					   do j = 1,3   

                         do k = 1,3

                            ak(j,k)=ak(j,k)/al

                         enddo

                       enddo

                          

					   do j = 1,3   

                          do k = 1,3

						     d(j,k)=0.0d0

                             a(j,k)=ak(j,k)

                          enddo

                       enddo



					   do j = 1,3   

                          d(j,j)=1.

                       enddo





                     !      MATRICES CLEARING

					   do j = 1,6   

                          dep(j)=0.0d0

                       enddo

     

                       call leg(s,ak,dep,D,DV,dt,csi,bef,pc,u)



	                   call upd(a,ep,ak,dep,D,DV,csi,z14,bef,pc,u)



					   call DENSITA(dens,dep,valvpsl(which_vpsl,14),valvpsd(which_vpsd,14))

                      

					   do j = 1,3   

                          a1(j)=ak(j,j)

                       enddo

                       a1(4)=ak(1,2)

                       a1(5)=ak(1,3)

                       a1(6)=ak(2,3)



                       a1n1(ip,iq,ielem)=a1(1)

                       a1n2(ip,iq,ielem)=a1(2)

					   a1n3(ip,iq,ielem)=a1(3)

					   a1n4(ip,iq,ielem)=a1(4)

					   a1n5(ip,iq,ielem)=a1(5)

					   a1n6(ip,iq,ielem)=a1(6)



					   epn1(ip,iq,ielem)=ep(1)

					   epn2(ip,iq,ielem)=ep(2)

                       epn3(ip,iq,ielem)=ep(3)

                       epn4(ip,iq,ielem)=ep(4)

					   epn5(ip,iq,ielem)=ep(5)

                       epn6(ip,iq,ielem)=ep(6)

	

                        u7n(ip,iq,ielem)=u(7)  

						z14n(ip,iq,ielem)=z14

						densit(ip,iq,ielem)=dens



  

                   ! conversion of stress or strain vectors according to conventions

      

                        destn1=-dep(1)

						destn2=-dep(2)

						destn3=-dep(4)

       

           

    				! Incremental viscoplastic strain rate and the total viscopl.strain

        			! Accumulate the absolute value of the viscoplastic strain increment



                

					!write(54,*)ip, iq, ielem, debar



					debar=dsqrt((2.0*(destn1**2 + destn2**2 + destn3**2)/3.0))

				

					vivel_el(ip,iq,ielem) =  vivel_el(ip,iq,ielem) + debar





					vistn1_0(ip,iq,ielem) = vistn1_0(ip,iq,ielem) + destn1 

					vistn2_0(ip,iq,ielem) = vistn2_0(ip,iq,ielem) + destn2 

					vistn3_0(ip,iq,ielem) = vistn3_0(ip,iq,ielem) + destn3 



					enddo !ip

				enddo !iq



				endif

			  enddo !ivpsd



			endif

		enddo !ivpsl



		endif

	enddo !ivpsa



	! ****************************



	endif ! which_vpsa



	 

	! MCNAPvpsa   - end







	 do iq = 1,nn   

         do ip = 1,nn



		 	sxx(ip,iq) = sxx(ip,iq) + sxx_ps(ip,iq,ielem) ! Marco	

		    syy(ip,iq) = syy(ip,iq) + syy_ps(ip,iq,ielem) ! Marco

			sxy(ip,iq) = sxy(ip,iq) + sxy_ps(ip,iq,ielem) ! Marco

!			szz(ip,iq) = szz(ip,iq) + szz_ps(ip,iq,ielem) ! Clara 12/01/2005



         enddo

      enddo

      
! FORCE CALCULATION
      
      do iq = 1,nn
         do ip = 1,nn
            t1fx = 0.0d0; t1fy = 0.0d0
            t2fx = 0.0d0; t2fy = 0.0d0
            
            do il = 1,nn
               t1fx = t1fx + dd(il,ip) * ww(il)*ww(iq) &
                    * (sxx(il,iq)*dydy(il) - sxy(il,iq)*dxdy(il))
               t1fy = t1fy + dd(il,ip) * ww(il)*ww(iq) &
                    * (sxy(il,iq)*dydy(il) - syy(il,iq)*dxdy(il))
            enddo
            
            do im = 1,nn
               t2fx = t2fx + dd(im,iq) * ww(ip)*ww(im) &
                    * (sxx(ip,im)*dydx(im) - sxy(ip,im)*dxdx(im))
               t2fy = t2fy + dd(im,iq) * ww(ip)*ww(im) &
                    * (sxy(ip,im)*dydx(im) - syy(ip,im)*dxdx(im))
            enddo
            
            det_j = dxdx(iq)*dydy(ip) - dxdy(ip)*dydx(iq)
            
            fx(ip,iq) = t1fx - t2fx
            fy(ip,iq) = t1fy - t2fy
         enddo
      enddo
      
      
      return
      
      end subroutine make_internal_force
      
      
!       ********************************************************
      
        subroutine deltat_max(time_step,nnode,nm,tm,pm,sdeg,&
                            xx,yy,nquad,con_quad,time_step_cfl,fmax)   
      
!       © CRS4, 2002, All Rights Reserved  
!       Authors: Marco Stupazzini
!     
        implicit none
         
        integer*4 :: nm,im
        integer*4, dimension(nm) :: tm,sdeg
        real*8, dimension(nm,4) :: pm
        integer*4 :: nnode,nquad
        real*8, dimension(nnode) :: xx,yy
        integer*4, dimension(nquad,5) :: con_quad
        integer*4, dimension(12):: vet_node
        real*8 :: length_min,length,vs_length_min,vs_length,length_vp_min,length_vp
        real*8 :: time_step,time_step_cfl,percent_deltat
        real*8 :: num_of_points,fmax,vs_npoints,vp_deltat_cfl
        real*8 :: rho,lambda,mu,vs,vp
        
        integer*4 :: iquad,conta
        integer*4 :: mat_deltat_cfl,sdeg_deltat_cfl
        integer*4 :: mat__npoints,sdeg_npoints
        
        integer*4 :: ie,i,j,mcode,smcode,nn

        real*8, dimension(:), allocatable :: ct,ww
        real*8, dimension(:,:), allocatable :: dd
	
	character*3 :: deltat_fixed
	      
        vet_node= (/2,2,2,3,3,4,3,4,5,4,5,5/)
        conta=1
!        fmax=10

!ACHTUNG - era commentato da qui fino a......................
       iquad=1

       do j = 1,nm
               if (tm(j).eq.con_quad(iquad,1)) im=j
        enddo
        
	smcode=sdeg(im)        
        rho = pm(im,1)   
        lambda = pm(im,2)
        mu = pm(im,3)
        vs = dsqrt(mu/rho)
        vp = dsqrt((lambda+2*mu)/rho)


       !The program assumes the first value of length_min the distance between
       !       nodes 1 2 belonging to the first quad

       length=dsqrt((xx(con_quad(iquad,vet_node(conta))) &
               -xx(con_quad(iquad,vet_node(conta+6))))**2 &
                +(yy(con_quad(iquad,vet_node(conta))) &
                -yy(con_quad(iquad,vet_node(conta+6))))**2)
                        
        length_vp_min = length/vp
	vp_deltat_cfl = vp
	sdeg_deltat_cfl = smcode

        vs_length_min = vs/length
	vs_npoints = vs
	sdeg_npoints = smcode


!ACHTUNG  - ....................................fino a qui
                       
      
        do iquad=1,nquad
                mcode=con_quad(iquad,1)

!                write(47,'(A,I5,A,I5,A,I5)')'iquad =',iquad,' | mcode =',mcode,' |smcode =',smcode
!                write(48,'(A,I5,A,I5,A,I5)')'iquad =',iquad,' | mcode =',mcode,' | smcode =',smcode 
 
                do j = 1,nm
                    if (tm(j).eq.con_quad(iquad,1)) im=j
                enddo
                
                smcode = sdeg(im)
                rho = pm(im,1)
                lambda = pm(im,2)
                mu = pm(im,3)
                vs = dsqrt(mu/rho)
                vp = dsqrt((lambda+2*mu)/rho)
 
                !Min Length check

                do conta=1,6
!                       write(47,'(A,I5)')'conta :',conta
!			write(48,'(A,I5)')'conta :',conta

                        length=dsqrt((xx(con_quad(iquad,vet_node(conta))) &
                                -xx(con_quad(iquad,vet_node(conta+6))))**2 &
                                +(yy(con_quad(iquad,vet_node(conta))) &
                                -yy(con_quad(iquad,vet_node(conta+6))))**2)
                        
                        
                        length_vp = length/vp

                        vs_length = vs/length

                        !Check on the minimum length/vp of the element                  

                        if (length_vp_min.gt.length_vp) then
                                length_vp_min = length_vp
                                vp_deltat_cfl = vp
                                sdeg_deltat_cfl = smcode
!				write(47,'(A,I5)')'*** conta :',conta
                        endif

                        !Check on the number of points per element
                        !vs_length_min = Minimun ratio between the maximum element length (length) and
                        !                minimum shear velocity (vs)
                        !mat_length_npoints = label of the material associated with the vs_lenght_min
                        !sdeg_mat_length_npoints = spectral degree of the material associated with the vs_lenght_min
                        
                        if (vs_length_min.gt.vs_length) then
                                vs_length_min = vs_length
                                vs_npoints = vs
                                sdeg_npoints = smcode
!				write(48,'(A,I5)')'*** conta :',conta
                        endif
                                
               enddo
!               write(47,'(A,I5,A,E12.4,A,E12.4,A,E12.4)')'iquad =',iquad,&
!                       ' |vp_deltat_cfl:',vp_deltat_cfl,&
!                       ' |length_vp_min:',length_vp_min,&
!                       ' |lenght:',length_vp_min*vp_deltat_cfl
!                       
!		write(48,'(A,I5,A,E12.4,A,E12.4,A,E12.4)')'iquad =',iquad,&
!                       ' |vs_npoints:',vs_npoints,&
!                       ' |vs_length_min:',vs_length_min,&
!                       ' |lenght:',1/vs_length_min*vs_npoints

        enddo

        
        
        nn=sdeg_deltat_cfl+1

        allocate(ct(nn),ww(nn),dd(nn,nn))
        call lgl(nn,ct,ww,dd)   
        
        time_step_cfl=length_vp_min*0.5d0*(ct(2)-ct(1))
        write(*,'(A)')' '
        write(*,'(A)')'Stability for time advancing algorithm '
        write(*,'(A,E12.4)')'Min. element length :',length_vp_min*vp_deltat_cfl
        write(*,'(A,E12.4)')'Min. pressure wave velocity :',vp_deltat_cfl
        write(*,'(A)') '*******************************************************'
        if (time_step.le.time_step_cfl) then
                percent_deltat=time_step/time_step_cfl*100
                if (percent_deltat.le.1) then
                        write(*,'(A,E12.4)')'WARNING!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                        write(*,'(A)')'This time step is excedingly lower the deltat CFL'
			if (deltat_fixed.eq.'not') then	 
                        	write(*,'(A)')'deltat chosen will be substituted with 1% of deltat CFL'
                        	time_step=time_step_cfl*0.01
                        	write(*,'(A,E12.4)')'deltat chosen :',time_step
			endif
                elseif (percent_deltat.le.25) then
                        write(*,'(A,E12.4)')'OK!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                else
                        write(*,'(A,E12.4)')'WARNING!!! deltat CFL =',&
                        time_step_cfl
                        write(*,'(A,F6.2,A)')'You are using ',percent_deltat,&
                                '% of critical time step.'
                        write(*,'(A)')'This could be not enough for the correct working of ABC'
			if (deltat_fixed.eq.'not') then
                        	write(*,'(A)')'deltat chosen will be substituted with 25% of deltat CFL'
                        	time_step=time_step_cfl*0.25
                        	write(*,'(A,E12.4)')'deltat chosen :',time_step
			endif
                endif
        elseif (time_step.gt.time_step_cfl) then
                write(*,'(2X,A,E12.4)')'ERROR!!! deltat CFL = ',&
                time_step_cfl,' < deltat = ',time_step
                write(*,'(A)')'The time advancing scheme will be unstable!'
		if (deltat_fixed.eq.'not') then
                	write(*,'(A)')'deltat chosen will be substituted with 25% of deltat CFL'
                	time_step=time_step_cfl*0.25
                	write(*,'(A,E12.4)')'deltat chosen :',time_step
		endif
        endif
        write(*,'(A)')'*******************************************************'
        write(*,'(A)')' '

        !Writing of the results for the maximum number of points per wave length
        
        nn=sdeg_npoints+1
        deallocate(ct,ww,dd)
	allocate(ct(nn),ww(nn),dd(nn,nn))
        call lgl(nn,ct,ww,dd)

        num_of_points = vs_length_min*(ct(1)-ct(nn))/(ct(int(nn/2))-ct(int(nn/2)+1))/fmax

        write(*,'(A)')'Number of points per wave length '
        write(*,'(A,E12.4)')'Max. element length :',1/vs_length_min*vs_npoints
        write(*,'(A,E12.4)')'Max. shear wave velocity :',vs_npoints
        write(*,'(A,E12.4)')'Number of points per wave length :',num_of_points
        write(*,'(A)')' '

!	write(48,'(A,E12.4,E12.4,E12.4,E12.4)')'ct =',ct(1),&
!		ct(nn),ct(int(nn/2)),ct(int(nn/2)+1)

        return
               
        end subroutine deltat_max


!     ********************************************************
      
      subroutine make_absorbing_force(lambda,mu,rho,nn,ct,ww,dd,&
                                      dxdx,dxdy,dydx,dydy,&
                                      nnx,nny,ll,ia,ja,ib,jb,&
                                      ux,uy,vx,vy,fx,fy)
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Fabio Maggio, Luca Massidda, Javier Sabadell
!     
      implicit none
      
      real*8 :: lambda,mu,rho
      integer*4 :: nn
      real*8, dimension(nn) :: ct,ww
      real*8, dimension(nn,nn) :: dd
      real*8, dimension(nn) :: dxdx,dxdy,dydx,dydy
      real*8 :: nnx,nny,ll
      integer*4 :: ia,ja,ib,jb
      real*8, dimension(nn,nn) :: ux,uy,vx,vy
      real*8, dimension(nn,nn) :: fx,fy
      
      integer*4 :: ip,iq,il,im,i,j
      real*8 :: det_j,nnx2,nny2,term,term11,term12,term21,term22,term_k
      
      do j = 1,nn
         do i = 1,nn
            fx(i,j) = 0.0d0
            fy(i,j) = 0.0d0
         enddo
      enddo
      
      nnx2 = nnx * nnx
      nny2 = nny * nny
      
      term = 0.5d0*ll * (2.0d0*mu - dsqrt(mu * (lambda+2.0d0*mu)))
      
      term11 = nnx2 * dsqrt(lambda +2.0d0*mu) + nny2 * dsqrt(mu)
      term11 = -0.5d0 * ll * dsqrt(rho) * term11
      
      term12 = (dsqrt(mu) - dsqrt(lambda +2.0d0*mu)) * nnx * nny
      term12 = 0.5d0 * ll * dsqrt(rho) * term12
      
      term21 = nny2 * dsqrt(lambda +2.0d0*mu) + nnx2 * dsqrt(mu)
      term21 = -0.5d0 * ll * dsqrt(rho) * term21
      
      term22 = (dsqrt(mu) - dsqrt(lambda +2.0d0*mu)) * nnx * nny
      term22 = 0.5d0 * ll * dsqrt(rho) * term22
      
      j = 0
      
      do il = ia,ib
         do im = ja,jb
            j = j +1
            
            det_j = dxdx(im)*dydy(il) - dxdy(il)*dydx(im)
            
            term_k = term * (dydy(il)*nny + dxdy(il)*nnx) &
                   * ww(j) / det_j
            
            iq = im
            do ip = 1,nn
               fx(il,im) = fx(il,im) -(term_k * dd(il,ip) * uy(ip,iq))
               fy(il,im) = fy(il,im) +(term_k * dd(il,ip) * ux(ip,iq))
            enddo
            
            term_k = -1.0d0 * term * (dydx(im)*nny + dxdx(im)*nnx) &
                    * ww(j) / det_j
            
            ip = il
            do iq = 1,nn
               fx(il,im) = fx(il,im) -(term_k * dd(im,iq) * uy(ip,iq))
               fy(il,im) = fy(il,im) +(term_k * dd(im,iq) * ux(ip,iq))
            enddo
            
            fx(il,im) = fx(il,im) -term11*ww(j) *vx(il,im)
            fx(il,im) = fx(il,im) -term12*ww(j) *vy(il,im)
            
            fy(il,im) = fy(il,im) -term21*ww(j) *vy(il,im)
            fy(il,im) = fy(il,im) -term22*ww(j) *vx(il,im)
            
         enddo
      enddo
      
      return
      
      end subroutine make_absorbing_force
      
      

!**************************************************************************
    subroutine invarg(devia1,devia2,devia3,fdatm, hards, gamma, ddelta,&
            frict, steff, stemp1,stemp2,stemp3, theta, varj2, yield, ncrit)
      
        ! © Politecnico di Milano  - DIS, 2003, All Rights Reserved
        ! Authors: Manollo Pastor, Clara Zambelli, Marco Stupazzini
        ! Comments:
        ! Obtain: (a)  Deviatoric stresses DEVIA
        !         (b)  Invariants
        !         (c)  Size of F passing through Stress point 
    
        implicit none
      
        real*8 :: devia1,devia2,devia3
        real*8 :: fdatm,hards,gamma,ddelta,frict
        real*8 :: steff
        real*8 :: stemp1,stemp2,stemp3
        real*8 :: theta, varj2, yield
        integer*4 :: ncrit
        real*8 :: root3,smean,varj3,sint3,devia4,phira,snphi
        sint3 = 0.0d0
        root3 = sqrt(3.00)
        ! First stress invariant (plane stress)
        smean  = (stemp1 + stemp2 )/3.0
        ! Deviatoric stress state (plane stress)
        devia1 = stemp1 - smean
        devia2 = stemp2 - smean
        devia3 = stemp3
        ! Second deviatoric stress invariant (plane stress)
        varj2 = devia3**2 + 0.5 * (devia1**2 + devia2**2)
        ! Third deviatoric stress invariant (plane stress)
        varj3 = 0.0d0

!        ! Second deviatoric stress invariant (plane strain)
!        varj2 = 0.5 * (devia1**2 + devia2**2 + devia3**2)
!        ! Third deviatoric stress invariant (plane strain)
!        varj3 = devia1 * devia2 * devia3

        steff=dsqrt(varj2)

        if (steff.ne.(0.0))  then
            sint3=-3.0*root3*varj3/(2.0*varj2*steff)
            if(sint3.gt.1.0) then 
                sint3=1.0
            else
                sint3=0.0
            endif
        endif

        if(sint3.lt.(-1.0)) then
            sint3=-1.0
        endif
        if(sint3.gt.(1.0)) then
            sint3= 1.0
        endif
        
        theta=asin(sint3)/3.0
        ! Tresca     
        if (ncrit.eq.1) then
            yield=2.0*cos(theta)*steff
        ! Von Mises
        elseif (ncrit.eq.2) then
            yield=root3*steff
        ! Mohr-Coulomb
        elseif (ncrit.eq.3) then
            phira=frict*0.017453292
            snphi=sin(phira)
            yield=smean*snphi+steff*(cos(theta)-sin(theta)*snphi/root3)
        ! Drucker-prager
        elseif (ncrit.eq.4) then
            phira=frict*0.017453292
            snphi=sin(phira)
            yield=6.0*smean*snphi/(root3*(3.0-snphi))+steff
        endif
        return
    end subroutine invarg



!**************************************************************************
    subroutine yieldf(avect1,avect2,avect3,devia1,devia2,devia3, &
        fdatm,hards,gamma,ddelta,frictpot,steff,theta,varj2,ncrit)

        ! © Politecnico di Milano  - DIS, 2003, All Rights Reserved
        ! Authors: Manollo Pastor, Clara Zambelli, Marco Stupazzini
        ! Comments:
        ! Obtain: plastic vectors df/ds

        implicit none

        real*8 :: avect1,avect2,avect3
        real*8 :: devia1,devia2,devia3
        real*8 :: fdatm,hards,gamma,ddelta,frictpot  !!Clara_non_ass
        real*8 :: steff
        real*8 :: theta, varj2
        integer*4 :: ncrit

        real*8 :: veca11,veca12,veca13
        real*8 :: veca21,veca22,veca23
        real*8 :: veca31,veca32,veca33

        real*8 :: sinth,costh,tanth,tant3,cost3,root3
        real*8 :: cons1,cons2,cons3,abthe,plumi,snphi



        if(steff.eq.0.0) return

            ! 2D and 3D

            tanth=tan(theta)
            tant3=tan(3.0*theta)
            sinth=sin(theta)
            costh=cos(theta)
            cost3=cos(3.0*theta)
            root3=sqrt(3.0)

            ! 2D situations
            !
            ! Get vector a1=d(Sm)/ds

            veca11=1.0
            veca12=1.0
            veca13=0.0

            ! Get vector a2=d(Steff)/ds

            veca21=devia1/(2.0*steff)
            veca22=devia2/(2.0*steff)
            veca23=devia3/steff

            ! Get vector a3

            veca31=varj2/3.0
            veca32=varj2/3.0
            veca33=0.0


            ! Tresca     
            if (ncrit.eq.1) then
                cons1=0.0
                abthe=abs(theta*57.29577951308)

                if(abthe.lt.29.0) then
                    cons2=2.0*(costh+sinth*tant3)
                    cons3=root3*sinth/(varj2*cost3)
                else
                    cons2=root3
                    cons3=0.0
                endif

            ! Von Mises
            elseif (ncrit.eq.2) then
                cons1=0.0
                cons2=root3
                cons3=0.0

            ! Mohr-Coulomb
            elseif (ncrit.eq.3) then
                cons1=sin(frictpot*0.017453292)/3.0    !!Clara_non_ass
                abthe=abs(theta*57.29577951308)

                if(abthe.lt.29.0) then
                    cons2=costh*((1.0+tanth*tant3)+cons1*(tant3-tanth)*root3)
                    cons3=(root3*sinth+3.0*cons1*costh)/(2.0*varj2*cost3)
                else
                    cons3=0.0
                    plumi=1.0

                    if(theta.gt.0.0) then
                        plumi=-1.0
                    endif

                    cons2=0.5*(root3+plumi*cons1*root3)
                endif


            ! Drucker-prager
            elseif (ncrit.eq.4) then

                snphi=sin(frictpot*0.017453292)      !!Clara_non_ass
                cons1=2.0*snphi/(root3*(3.0-snphi))
                cons2=1.0
                cons3=0.0

            endif

        avect1=cons1 * veca11 + cons2 * veca21 + cons3 * veca31
        avect2=cons1 * veca12 + cons2 * veca22 + cons3 * veca32
        avect3=cons1 * veca13 + cons2 * veca23 + cons3 * veca33

        return

    end subroutine yieldf



!     ********************************************************
      
      subroutine flowvp(fdatm,hards,gamma, &
                        ddelta,frict, &
			            devia1,devia2,devia3, &
			            avect1,avect2,avect3, &
						vivel1,vivel2,vivel3, &
						vivel,yield,npoin,ncrit)
      
!     © Politecnico di Milano  - DIS, 2003, All Rights Reserved
!     Authors: Manollo Pastor, Clara Zambelli, Marco Stupazzini
!

      implicit none
      
	  real*8 :: devia1,devia2,devia3
	  real*8 :: avect1,avect2,avect3
      real*8 :: fdatm,hards,gamma,ddelta,frict
	  integer*4 :: npoin
	  real*8 :: vivel1,vivel2,vivel3
      real*8 :: vivel
      real*8 :: yield,fdatm1
      integer*4 :: ncrit
      
	  real*8 :: allow,root3,fcurr,fnorm,istr1,cmult



		allow=0.01
		root3=sqrt(3.0)
	!	frict=frict*0.017453292
	    fdatm1=fdatm

		if(ncrit.eq.3) then		 
			fdatm1=fdatm1*cos(frict*0.017453292)
		elseif(ncrit.eq.4) then
			fdatm1=6*fdatm1*cos(frict*0.017453292)/(root3*(3.0-sin(frict*0.017453292)))
	    endif

		fdatm1=fdatm1 + vivel * hards

	    fcurr=yield-fdatm1
		fnorm=fcurr/fdatm1

      if (fnorm.lt.0.0) then
	     vivel1 = 0.0
         vivel2 = 0.0
		 vivel3 = 0.0
	  else
	      cmult= gamma*(fnorm**ddelta)
		  avect1 = avect1 * cmult
		  avect2 = avect2 * cmult
		  avect3 = avect3 * cmult

          vivel1 = avect1
		  vivel2 = avect2
		  vivel3 = avect3

      endif

      return

      end subroutine flowvp



!!     ********************************************************
!      
!      subroutine write_GID_output_el(nnode,xs,ys,cs_nnz,cs,&
!                           nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
!                           ne,alfa1,alfa2,beta1,beta2,gamma1,gamma2,&
!                           opt,nsnaps,tsnaps,file_name,np)
!      
!!     © CRS4, 2002, All Rights Reserved
!!     Authors: Clara Zambelli, Marco S.
!!     
!      implicit none
!      
!      integer*4 :: nnode,cs_nnz,nm,ne,opt,nsnaps,np
!      real*8, dimension(nnode) :: xs,ys
!      integer*4, dimension(0:cs_nnz) :: cs
!      integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat
!      real*8, dimension(nm) :: tref_mat
!      real*8, dimension(nm,4) :: prop_mat
!      real*8, dimension(ne) :: alfa1,alfa2,beta1,beta2,gamma1,gamma2
!      real*8, dimension(nsnaps) :: tsnaps
!      character*70 :: file_name
!      
!      real*8, dimension(:), allocatable :: ux,uy,vivel_read
!      real*8, dimension(:), allocatable :: sxx,syy,szz,sxy
!      real*8, dimension(:), allocatable :: spr,svm
!      integer*4 :: ip,i,j,in,nd,nn,nelem
!      integer*4 :: imat,ie,mn,n1,n2,n3,n4,count,ic
!      real*8 :: val
!      integer*4 :: istep,nsteps      
!      character*70 :: in_fileE,in_fileV,out_file
!      character*40 :: elcode
!      integer*4 :: lname
!      
!      real*8, dimension(:), allocatable :: duxdx_average,duydx_average,duxdy_average,duydy_average
!
!      lname = len_trim(file_name)
!      in_fileE = file_name(1:lname) // 'E_xxx_xxx.bin'
!	  in_fileV = file_name(1:lname) // 'V_xxx_xxx.bin'
!      out_file = file_name(1:lname) // '.inp'
!
!      istep = 1
!      allocate(ux(nnode),uy(nnode),vivel_read(nnode))
!      allocate(sxx(nnode),syy(nnode),szz(nnode),sxy(nnode))
!      allocate(spr(nnode),svm(nnode))
!	  allocate(duxdx_average(nnode),duydx_average(nnode),duxdy_average(nnode),duydy_average(nnode))
!      
!      if (opt.eq.1) then
!         nsteps = nsnaps
!      elseif (opt.eq.2) then
!         nsteps = nsnaps*2
!      elseif (opt.eq.3) then
!         nsteps = nsnaps*3
!	  elseif (opt.eq.7) then
!         nsteps = nsnaps*2
!      endif
!
!
!
!      !write(21,'(I2)')nsteps
!      !write(21,'(A)')'data'
!      
!      do count = 1,nsnaps
!
!!         if (istep.lt.10) then
!!           write(21,'(A4,I1,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!!         else
!!            write(21,'(A4,I2,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!!         endif
!!         
!         if (istep.eq.1) then
!		     open(unit=1021,file='sperem.dat')
!
!             write(1021,'(A)') ' line 1: output of bending'
!             write(1021,'(A)') ' line 2: output of bending'
!	         write(1021,'(A)') ' line 3: output of bending'
!	         write(1021,'(A)') ' line 4: output of bending'
!	         write(1021,'(A)') ' line 5: output of bending'
!	         write(1021,'(A)') 'nelem   npoin   nelem_type'
!
!            ne = cs(0) -1
!            
!            nelem = 0
!            do imat = 1,nm
!               nn = sdeg_mat(imat) +1
!               do ie = 1,ne
!                  if (cs(cs(ie -1) +0).eq.tag_mat(imat)) then
!                     nelem = nelem + (nn -1)*(nn -1)
!                  endif
!               enddo
!            enddo
!            
!
!            write(1021,'(I8,2X,I8,2X,I8)') nelem, nnode, 4
!			
!			write(1021,'(A)')' --- coordenadas ---'
!
!            do i = 1,nnode
!               write(1021,'(I8,2X,E14.6,2X,E14.6)') i,xs(i),ys(i)
!            enddo
!
!            write(1021,*)' --- ielem    n1   n2    n3    n4'
!            
!            ic = 0
!            !elcode = ''
!            
!            do imat = 1,nm
!               nn = sdeg_mat(imat) +1
!               
!               do ie = 1,ne
!                  if (cs(cs(ie -1) +0).eq.tag_mat(imat)) then
!                     mn = cs(cs(ie -1) +0)
!                     
!                     do j = 1,nn -1
!                        do i = 1,nn -1
!                           n1 = cs(cs(ie -1) +nn*(j -1) +i)
!                           n2 = cs(cs(ie -1) +nn*(j -1) +i +1)
!                           n3 = cs(cs(ie -1) +nn*j +i +1)
!                           n4 = cs(cs(ie -1) +nn*j +i)
!                           
!                           ic = ic +1
!                           write(1021,'(I8,2X,I8,2X,I8,2X,I8,2X,I8)')&
!                                 ic,n1,n2,n3,n4
!                        enddo
!                     enddo
!                  endif
!               enddo
!            enddo
!
!			close(1021)
!            
!         endif
!         
!
!         open(unit=1022,file='sperem.res')
!
!	!	 write(1022,'(A15,E14.6,A8)')'   displ      1',istep,'2   1  0'      
!		 write(1022,*) '   displ      1 ', istep, '  2  1    0'
!
!
!
!		if (count.lt.10) then
!            write(in_fileE(lname+3:lname+4),'(a2)')'00'
!            write(in_fileE(lname+5:lname+5),'(i1)')count
!
!			write(in_fileV(lname+3:lname+4),'(a2)')'00'
!            write(in_fileV(lname+5:lname+5),'(i1)')count
!	     else if (count.lt.100) then
!            write(in_fileE(lname+3:lname+3),'(a1)')'0'
!            write(in_fileE(lname+4:lname+5),'(i2)')count
!
!			write(in_fileV(lname+3:lname+3),'(a1)')'0'
!            write(in_fileV(lname+4:lname+5),'(i2)')count
!
!         else if (count.le.999) then
!            write(in_fileE(lname+3:lname+5),'(i3)')count
!
!			write(in_fileV(lname+3:lname+5),'(i3)')count
!         endif
!         
!         do ip = 0,np -1
!            if (ip.lt.10) then
!               write(in_fileE(lname+7:lname+8),'(a2)')'00'
!               write(in_fileE(lname+9:lname+9),'(i1)')ip
!
!               write(in_fileV(lname+7:lname+8),'(a2)')'00'
!               write(in_fileV(lname+9:lname+9),'(i1)')ip
!
!			else if (ip.lt.100) then
!               write(in_fileE(lname+7:lname+7),'(a1)')'0'
!               write(in_fileE(lname+8:lname+9),'(i2)')ip
!
!			   write(in_fileV(lname+7:lname+7),'(a1)')'0'
!               write(in_fileV(lname+8:lname+9),'(i2)')ip
!
!            else if (ip.le.999) then
!               write(in_fileE(lname+7:lname+9),'(i3)')ip
!			   write(in_fileV(lname+7:lname+9),'(i3)')ip
!
!            endif
!
!            
!            open(20,file=in_fileE)
!            
!            read(20,*)nd
!            
!            do i = 1,nd
!               read(20,*)in,val
!               if (dabs(val).lt.1.d-30) val=0.d0
!               if (in.le.nnode) then
!                  ux(in) = val
!               else
!                  uy(in -nnode) = val
!               endif
!            enddo
!            
!            close(20)
!
!
!            open(20,file=in_fileV)
!            
!            read(20,*)nd
!            
!            do i = 1,nnode
!               read(20,*)in,val
!               if (dabs(val).lt.1.d-30) val=0.d0
!               vivel_read(i) = val
!
!            enddo
!            
!            close(20)
!
!         enddo
!         
!         
!         !write(21,*)'3 0 '
!         !write(21,*)'1 3 '
!         !write(21,*)'Displacement, NO_UNITS'
!         
!         do i = 1,nnode
!            write(1022,'(I8,2X,E14.6,2X,E14.6)')&
!                 i,ux(i),uy(i)
!         enddo
!         
!         istep = istep + 1
!         
!         if ((opt.eq.2).or.(opt.eq.3)) then
!
!         !   write(1022,'(A13,E14.6,A8)')'   str      1',istep,'3   1  0'
!			write(1022,*) '   str         1 ', istep-1, '  3   1   0'
!
!            call calculate_element_results_el(nnode,xs,ys,cs_nnz,cs,&
!                           nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
!                           ne,alfa1,alfa2,beta1,beta2,gamma1,gamma2,&
!                           ux,uy,sxx,syy,szz,sxy,spr,svm,&
!						   duxdx_average,duydx_average,duxdy_average,duydy_average)
!            
!            if (istep.lt.10) then
!               !write(21,'(A4,I1,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!            else
!               !write(21,'(A4,I2,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!            endif
!            
!            !write(21,*)'4 0 '
!            !write(21,*)'4 1 1 1 1 '
!            !write(21,*)'SXX, NO_UNITS'
!            !write(21,*)'SYY, NO_UNITS'
!            !write(21,*)'SZZ, NO_UNITS'
!            !write(21,*)'SXY, NO_UNITS'
!            
!            do i = 1,nnode
!               write(1022,'(I8,2X,E14.6,2X,E14.6,2X,E14.6)')&
!                    i,sxx(i),syy(i),sxy(i)
!            enddo
!
!            write(1022,*) '   strain       1 ', istep-1, '  3   1   0'
!
!			do i = 1,nnode
!               write(1022,'(I8,2X,E14.6,2X,E14.6,2X,E14.6)')&
!                    i,duxdx_average(i),duydy_average(i),duxdy_average(i)
!            enddo
!
!			write(1022,*) '   vivel       1 ', istep-1, '  1   1   0'
!
!			do i = 1,nnode
!               write(1022,'(I8,2X,E14.6)')&
!                    i,vivel_read(i)
!            enddo
!            
!            istep = istep + 1
!         endif
!         
!         if (opt.eq.3) then
!            if (istep.lt.10) then
!               write(21,'(A4,I1,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!            else
!               write(21,'(A4,I2,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!            endif
!            
!            write(21,*)'2 0 '
!            write(21,*)'2 1 1 '
!            write(21,*)'Pressure, NO_UNITS'
!            write(21,*)'VonMises, NO_UNITS'
!            
!            do i = 1,nnode
!               write(21,'(I8,2X,E14.6,2X,E14.6)')&
!                    i,spr(i),svm(i)
!            enddo
!            
!            istep = istep + 1
!		  endif
!
!		  if (opt.eq.7) then
!            if (istep.lt.10) then
!               write(21,'(A4,I1,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!            else
!               write(21,'(A4,I2,A3,E14.6)')'step',istep,' , ',tsnaps(count)
!            endif
!            
!            write(21,*)'1 0 '
!            write(21,*)'1 1 '
!            write(21,*)'Vivel, NO_UNITS'
!            
!            do i = 1,nnode
!               write(21,'(I8,2X,E14.6)')&
!                    i,vivel_read(i)
!            enddo
!            
!            istep = istep + 1
!         endif
!         
!      enddo
!      
!      close(1022)
!      
!      deallocate(ux,uy)
!      deallocate(sxx,syy,szz,sxy)
!      deallocate(spr,svm)
!	  deallocate(duxdx_average,duydx_average,duxdy_average,duydy_average)
!      
!      return
!      
!      end subroutine write_GID_output_el
!      
!      
!      
!!     ********************************************************
      
      subroutine write_UCD_GID_output_el(nnode,xs,ys,cs_nnz,cs,&
                           nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                           ne,alfa1,alfa2,beta1,beta2,gamma1,gamma2,&
                           opt,nsnaps,tsnaps,file_name,np,&

						   nmaps,nsnapps,tsnapsps,fmaps) ! Marco
      
!     © CRS4, 2002, All Rights Reserved
!     Authors: Luca Massidda
!     © POLIMI, 2004, All Rights Reserved
!     Authors: Clara Zambelli, Marco S.
!     Controllare output variabili geotecniche:p,q,dev,dep --->>CLARA!
     
      implicit none
      
      integer*4 :: nnode,cs_nnz,nm,ne,opt,nsnaps,np
      real*8, dimension(nnode) :: xs,ys
      integer*4, dimension(0:cs_nnz) :: cs
      integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat
      real*8, dimension(nm) :: tref_mat
      real*8, dimension(nm,4) :: prop_mat
      real*8, dimension(ne) :: alfa1,alfa2,beta1,beta2,gamma1,gamma2
      real*8, dimension(nsnaps) :: tsnaps
      character*70 :: file_name
      
      real*8, dimension(:), allocatable :: ux,uy,vivel_read,vistn1_0_read,vistn2_0_read,vistn3_0_read
      real*8, dimension(:), allocatable :: sxx,syy,szz,sxy

	  real*8, dimension(:), allocatable :: sxx_ps,syy_ps,sxy_ps,szz_ps ! Marco

	  !real*8, dimension(:), allocatable :: vistnxx_ps,vistnyy_ps,vistnxy_ps ! MCNAP
      real*8, dimension(:), allocatable :: spr,svm
      integer*4 :: ip,i,j,in,nd,nn,nelem
      integer*4 :: imat,ie,mn,n1,n2,n3,n4,count,ic
      real*8 :: val
      integer*4 :: istep,nsteps      
      character*70 :: in_fileE,in_fileV,in_fileT1,in_fileT2,in_fileT3,out_file
      character*40 :: elcode
      integer*4 :: lname
      

      real*8 ::p,q    !Clara

	  real*8, dimension (:), allocatable::value  !Clara


      real*8, dimension(:), allocatable :: duxdx_average,duydx_average,duxdy_average,duydy_average


	  integer*4 :: nmaps,imaps ! Marco

	  integer*4, dimension(1:nmaps) :: fmaps

	  integer*4:: nsnapps ! Marco

	  integer*4, dimension(nsnapps) :: tsnapsps ! Marco


      lname = len_trim(file_name)
      in_fileE = file_name(1:lname) // 'E_xxx_xxx.bin'
	  in_fileV = file_name(1:lname) // 'V_xxx_xxx.bin'
      in_fileT1 = file_name(1:lname) // '1_xxx_xxx.bin'
	  in_fileT2 = file_name(1:lname) // '2_xxx_xxx.bin'
	  in_fileT3 = file_name(1:lname) // '3_xxx_xxx.bin'
      out_file = file_name(1:lname) // '.inp'

      istep = 1
      allocate(ux(nnode),uy(nnode),vivel_read(nnode),vistn1_0_read(nnode),vistn2_0_read(nnode),vistn3_0_read(nnode))
      allocate(sxx(nnode),syy(nnode),szz(nnode),sxy(nnode))

	  allocate(sxx_ps(nnode),syy_ps(nnode),sxy_ps(nnode),szz_ps(nnode)) ! Marco

	  !allocate(vistnxx_ps(nnode),vistnyy_ps(nnode),vistnxy_ps(nnode)) ! MCNAP
      allocate(spr(nnode),svm(nnode))
	  allocate(duxdx_average(nnode),duydx_average(nnode),duxdy_average(nnode),duydy_average(nnode))
      

      allocate(value(4)) !Clara


      if (opt.eq.1) then
         nsteps = nsnaps
      elseif (opt.eq.2) then
         nsteps = nsnaps*2
      elseif (opt.eq.3) then
         nsteps = nsnaps*3
	  elseif (opt.eq.7) then
         nsteps = nsnaps*2
      endif





	  if (nmaps.gt.0) then

		do imaps = 1,nmaps

			if (fmaps(imaps).eq.3) then

				call read_pre_stress(nnode,file_name,nsnapps,tsnapsps,&

									 sxx_ps,syy_ps,sxy_ps,szz_ps,&

									 !vistnxx_ps,vistnyy_ps,vistnxy_ps,& !MCNAP

									 np)

		    endif

		enddo

	  endif

	     


      open(21,file=out_file)
      
      write(21,'(A)')'# Elastic else output'
      write(21,'(A)')'# CRS4'
      write(21,'(A)')'# Center for Advanced Studies, '
      write(21,'(A)')'# Research and Development in Sardinia'
      write(21,'(A)')'#'
      write(21,'(I2)')nsteps
      write(21,'(A)')'data'


      
      do count = 1,nsnaps

         if (istep.lt.10) then
           write(21,'(A4,I1,A3,E14.6)')'step',istep,' , ',tsnaps(count)
         else
            write(21,'(A4,I2,A3,E14.6)')'step',istep,' , ',tsnaps(count)
         endif

         
         if (istep.eq.1) then
		     open(unit=1021,file='sperem.dat')

             write(1021,'(A)') ' line 1: output of bending'
             write(1021,'(A)') ' line 2: output of bending'
	         write(1021,'(A)') ' line 3: output of bending'
	         write(1021,'(A)') ' line 4: output of bending'
	         write(1021,'(A)') ' line 5: output of bending'
	         write(1021,'(A)') 'nelem   npoin   nelem_type'

            ne = cs(0) -1
            
            nelem = 0
            do imat = 1,nm
               nn = sdeg_mat(imat) +1
               do ie = 1,ne
                  if (cs(cs(ie -1) +0).eq.tag_mat(imat)) then
                     nelem = nelem + (nn -1)*(nn -1)
                  endif
               enddo
            enddo
            
			write(21,'(I8,2X,I8)')nnode,nelem
            write(1021,'(I8,2X,I8,2X,I8)') nelem, nnode, 4
			
			write(1021,'(A)')' --- coordenadas ---'

            do i = 1,nnode
				write(21,'(I8,2X,E14.6,2X,E14.6,2X,E14.6)') i,xs(i),ys(i),0.0d0
				write(1021,'(I8,2X,E14.6,2X,E14.6)') i,xs(i),ys(i)
            enddo

            write(1021,*)' --- ielem    n1   n2    n3    n4'
            
            ic = 0
            elcode = 'quad'
            
            do imat = 1,nm
               nn = sdeg_mat(imat) +1
               
               do ie = 1,ne
                  if (cs(cs(ie -1) +0).eq.tag_mat(imat)) then
                     mn = cs(cs(ie -1) +0)
                     
                     do j = 1,nn -1
                        do i = 1,nn -1
                           n1 = cs(cs(ie -1) +nn*(j -1) +i)
                           n2 = cs(cs(ie -1) +nn*(j -1) +i +1)
                           n3 = cs(cs(ie -1) +nn*j +i +1)
                           n4 = cs(cs(ie -1) +nn*j +i)
                           
                           ic = ic +1
						   write(21,'(I8,2X,I4,2X,A6,2X,I8,2X,I8,2X,I8,2X,I8)')&
                                ic,mn,elcode,n1,n2,n3,n4
                           write(1021,'(I8,2X,I8,2X,I8,2X,I8,2X,I8)')&
                                 ic,n1,n2,n3,n4
                        enddo
                     enddo
                  endif
               enddo
            enddo

			close(1021)
            
         endif
         

         open(unit=1022,file='sperem.res')

    
		 write(1022,*) '   displ      1 ', tsnaps(count), '  2  1    0'

		if (count.lt.10) then
            write(in_fileE(lname+3:lname+4),'(a2)')'00'
            write(in_fileE(lname+5:lname+5),'(i1)')count

			write(in_fileV(lname+3:lname+4),'(a2)')'00'
            write(in_fileV(lname+5:lname+5),'(i1)')count

			write(in_fileT1(lname+3:lname+4),'(a2)')'00'
            write(in_fileT1(lname+5:lname+5),'(i1)')count

			write(in_fileT2(lname+3:lname+4),'(a2)')'00'
            write(in_fileT2(lname+5:lname+5),'(i1)')count

			write(in_fileT3(lname+3:lname+4),'(a2)')'00'
            write(in_fileT3(lname+5:lname+5),'(i1)')count



	     else if (count.lt.100) then
            write(in_fileE(lname+3:lname+3),'(a1)')'0'
            write(in_fileE(lname+4:lname+5),'(i2)')count

			write(in_fileV(lname+3:lname+3),'(a1)')'0'
            write(in_fileV(lname+4:lname+5),'(i2)')count

			write(in_fileT1(lname+3:lname+3),'(a1)')'0'
            write(in_fileT1(lname+4:lname+5),'(i2)')count

			write(in_fileT2(lname+3:lname+3),'(a1)')'0'
            write(in_fileT2(lname+4:lname+5),'(i2)')count

			write(in_fileT3(lname+3:lname+3),'(a1)')'0'
            write(in_fileT3(lname+4:lname+5),'(i2)')count




         else if (count.le.999) then
            write(in_fileE(lname+3:lname+5),'(i3)')count

			write(in_fileV(lname+3:lname+5),'(i3)')count

			write(in_fileT1(lname+3:lname+5),'(i3)')count

			write(in_fileT2(lname+3:lname+5),'(i3)')count

			write(in_fileT3(lname+3:lname+5),'(i3)')count


         endif
         
         do ip = 0,np -1
            if (ip.lt.10) then
               write(in_fileE(lname+7:lname+8),'(a2)')'00'
               write(in_fileE(lname+9:lname+9),'(i1)')ip

               write(in_fileV(lname+7:lname+8),'(a2)')'00'
               write(in_fileV(lname+9:lname+9),'(i1)')ip

			   write(in_fileT1(lname+7:lname+8),'(a2)')'00'
               write(in_fileT1(lname+9:lname+9),'(i1)')ip

			   write(in_fileT2(lname+7:lname+8),'(a2)')'00'
               write(in_fileT2(lname+9:lname+9),'(i1)')ip

			   write(in_fileT3(lname+7:lname+8),'(a2)')'00'
               write(in_fileT3(lname+9:lname+9),'(i1)')ip




			else if (ip.lt.100) then
               write(in_fileE(lname+7:lname+7),'(a1)')'0'
               write(in_fileE(lname+8:lname+9),'(i2)')ip

			   write(in_fileV(lname+7:lname+7),'(a1)')'0'
               write(in_fileV(lname+8:lname+9),'(i2)')ip

			   write(in_fileT1(lname+7:lname+7),'(a1)')'0'
               write(in_fileT1(lname+8:lname+9),'(i2)')ip

			   write(in_fileT2(lname+7:lname+7),'(a1)')'0'
               write(in_fileT2(lname+8:lname+9),'(i2)')ip

			   write(in_fileT3(lname+7:lname+7),'(a1)')'0'
               write(in_fileT3(lname+8:lname+9),'(i2)')ip



            else if (ip.le.999) then
               write(in_fileE(lname+7:lname+9),'(i3)')ip

				write(in_fileV(lname+7:lname+9),'(i3)')ip

				write(in_fileT1(lname+7:lname+9),'(i3)')ip
				write(in_fileT2(lname+7:lname+9),'(i3)')ip
				write(in_fileT3(lname+7:lname+9),'(i3)')ip



            endif

            
            open(20,file=in_fileE)
            
            read(20,*)nd
            
            do i = 1,nd
               read(20,*)in,val
               if (dabs(val).lt.1.d-30) val=0.d0
               if (in.le.nnode) then
                  ux(in) = val
               else
                  uy(in -nnode) = val
               endif
            enddo
            
            close(20)


            open(20,file=in_fileV)
            
            read(20,*)nd
            
            do i = 1,nnode
               read(20,*)in,val
               if (dabs(val).lt.1.d-30) val=0.d0
               vivel_read(i) = val

            enddo
            
            close(20)


			open(20,file=in_fileT1)
            
            read(20,*)nd
            
            do i = 1,nnode
               read(20,*)in,val
               if (dabs(val).lt.1.d-30) val=0.d0
               vistn1_0_read(i) = val

            enddo
            
            close(20)

			open(20,file=in_fileT2)
            
            read(20,*)nd
            
            do i = 1,nnode
               read(20,*)in,val
               if (dabs(val).lt.1.d-30) val=0.d0
               vistn2_0_read(i) = val

            enddo
            
            close(20)

			open(20,file=in_fileT3)

            

            read(20,*)nd

            

            do i = 1,nnode

               read(20,*)in,val

               if (dabs(val).lt.1.d-30) val=0.d0

               vistn3_0_read(i) = val



            enddo

            

            close(20)

         enddo
         
         
         write(21,*)'3 0 '
         write(21,*)'1 3 '
         write(21,*)'Displacement, NO_UNITS'
         
         do i = 1,nnode
		    write(21,'(I8,2X,E14.6,2X,E14.6,2X,E14.6)')&
                 i,ux(i),uy(i),0.0d0
            write(1022,'(I8,2X,E14.6,2X,E14.6)')&
                 i,ux(i),uy(i)
         enddo
         
         istep = istep + 1
         
         if ((opt.eq.2).or.(opt.eq.3).or.(opt.eq.7)) then

			write(1022,*) '   str         1 ', tsnaps(count), '  4   1   0'

            call calculate_element_results_el(nnode,xs,ys,cs_nnz,cs,&
                           nm,tag_mat,type_mat,sdeg_mat,tref_mat,prop_mat,&
                           ne,alfa1,alfa2,beta1,beta2,gamma1,gamma2,&
                           ux,uy,sxx,syy,szz,sxy,spr,svm,&
						   duxdx_average,duydx_average,duxdy_average,duydy_average,&
						   vistn1_0_read,vistn2_0_read,vistn3_0_read)
            
            if (istep.lt.10) then
               write(21,'(A4,I1,A3,E14.6)')'step',istep,' , ',tsnaps(count)
            else
               write(21,'(A4,I2,A3,E14.6)')'step',istep,' , ',tsnaps(count)
            endif
            
            write(21,*)'4 0 '
            write(21,*)'4 1 1 1 1 '
            write(21,*)'SXX, NO_UNITS'
            write(21,*)'SYY, NO_UNITS'
            write(21,*)'SZZ, NO_UNITS'
            write(21,*)'SXY, NO_UNITS'



			if (nmaps.gt.0) then

				do imaps = 1,nmaps

					if (fmaps(imaps).eq.3) then

						do i = 1,nnode

							sxx(i) = sxx(i) + sxx_ps(i)

		    				syy(i) = syy(i) + syy_ps(i)

							sxy(i) = sxy(i) + sxy_ps(i)

							szz(i) = szz(i) + szz_ps(i)    !Clara 12/01/2005 

						enddo

					endif

				enddo

		    endif
            
            do i = 1,nnode
			   write(21,'(I8,2X,E14.6,2X,E14.6,2X,E14.6,2X,E14.6)')&
                    i,sxx(i),syy(i),szz(i),sxy(i)
               write(1022,'(I8,2X,E14.6,2X,E14.6,2X,E14.6,2X,E14.6)')&
                    i,sxx(i),syy(i),sxy(i),szz(i)
            enddo

            write(1022,*) '   strain       1 ', tsnaps(count), '  3   1   0'

			do i = 1,nnode
               write(1022,'(I8,2X,E14.6,2X,E14.6,2X,E14.6)')&
                    i,duxdx_average(i),duydy_average(i),duxdy_average(i)
            enddo

!Clara begin



            write(1022,*) '   mean stress       1 ', tsnaps(count), '  1   1   0'

            do i = 1,nnode

			   value= [-sxx(i),-syy(i),-szz(i),-sxy(i)]   !definire value

               call inv_nc(4,value,p,q,1)  

               write(1022,'(I8,2X,E14.6)')i,p

			enddo

            

			write(1022,*) '   dev stress       1 ', tsnaps(count), '  1   1   0'

            do i = 1,nnode

			   value= [-sxx(i),-syy(i),-szz(i),-sxy(i)]   !definire value

               call inv_nc(4,value,p,q,1)  

               write(1022,'(I8,2X,E14.6)')i,q

			enddo

			 

            write(1022,*) '   Eps_Vol       1 ', tsnaps(count), '  1   1   0'

            do i = 1,nnode

			   value= [-duxdx_average(i),-duydy_average(i),0.00d0,-duxdy_average(i)]   !definire le deformazioni

               call inv_nc(4,value,p,q,2)  

               write(1022,'(I8,2X,E14.6)')i,p

			enddo

            

			write(1022,*) '   Eps_Dev      1 ', tsnaps(count), '  1   1   0'

            do i = 1,nnode

			   value= [-duxdx_average(i),-duydy_average(i),0.00d0,-duxdy_average(i)]   !definire le deformazioni

               call inv_nc(4,value,p,q,2)  

               write(1022,'(I8,2X,E14.6)')i,q

			enddo 



!Clara end
		 endif

         
         if (opt.eq.3) then
            if (istep.lt.10) then
               write(21,'(A4,I1,A3,E14.6)')'step',istep,' , ',tsnaps(count)
            else
               write(21,'(A4,I2,A3,E14.6)')'step',istep,' , ',tsnaps(count)
            endif
            
            write(21,*)'2 0 '
            write(21,*)'2 1 1 '
            write(21,*)'Pressure, NO_UNITS'
            write(21,*)'VonMises, NO_UNITS'
            
            do i = 1,nnode
               write(21,'(I8,2X,E14.6,2X,E14.6)')&
                    i,spr(i),svm(i)
            enddo
            
            istep = istep + 1
         endif

		 if (opt.eq.7) then
			
			write(1022,*) '   vivel       1 ', tsnaps(count), '  1   1   0'

            if (istep.lt.10) then
               write(21,'(A4,I1,A3,E14.6)')'step',istep,' , ',tsnaps(count)
            else
               write(21,'(A4,I2,A3,E14.6)')'step',istep,' , ',tsnaps(count)
            endif
            
            write(21,*)'1 0 '
            write(21,*)'1 1 '
            write(21,*)'Vivel, NO_UNITS'
            
            do i = 1,nnode
               write(21,'(I8,2X,E14.6)')&
                    i,vivel_read(i)
			   write(1022,'(I8,2X,E14.6)')&
                    i,vivel_read(i)
            enddo
            
            istep = istep + 1
         endif
         
      enddo
      
	  close(21)
      close(1022)
      
      deallocate(ux,uy,vivel_read)
      deallocate(sxx,syy,szz,sxy)
      deallocate(spr,svm)
	  deallocate(duxdx_average,duydx_average,duxdy_average,duydy_average)



	  deallocate(value)  !Clara
      
      return
      
      end subroutine write_UCD_GID_output_el
      
      
      
!     ********************************************************


!     ********************************************************

      

      subroutine make_pre_stress(rho,lambda,mu,&

	                             phi,fmaps,&

								 sxx_ps,syy_ps,sxy_ps,szz_ps,&

								 nn,ielem,nelem,&

								 nnt,xs,ys,&

								 cs_nnz,cs,node_index,&

								 N_update_el_az,update_index_el_az,&

								 file_name,nsnapps,tsnapsps,&

						         sxx_ps_read,syy_ps_read,sxy_ps_read,szz_ps_read,&

								 !vistnxx_ps_read,vistnyy_ps_read,vistnxy_ps_read,& ! MCNAP

								 nprocs)



      

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Marco Stupazzini (22/09/2004)

!     

      implicit none



	  integer*4:: nsnapps

	  integer*4, dimension(nsnapps) :: tsnapsps ! Marco

      real*8 :: rho,lambda,mu,phi,k0

      integer*4 :: nn

      integer*4 :: iq,ip,is,in



	  integer*4 :: nnt,ielem,nelem



	  real*8, dimension(nn,nn,nelem) :: sxx_ps 

	  real*8, dimension(nn,nn,nelem) :: syy_ps

	  real*8, dimension(nn,nn,nelem) :: sxy_ps

      real*8, dimension(nn,nn,nelem) :: szz_ps !Clara 12/01/2005



	  real*8, dimension(nnt) :: xs,ys



	  integer*4 :: cs_nnz

	  integer*4, dimension(0:cs_nnz) :: cs



	  integer*4 :: fmaps



      integer*4 :: N_update_el_az

      integer*4, dimension(nnt) :: node_index

	  integer*4, dimension(0:N_update_el_az -1) :: update_index_el_az

	  real*8, dimension(0:N_update_el_az -1) :: sxx_ps_read,syy_ps_read,sxy_ps_read,szz_ps_read

	  !real*8, dimension(0:N_update_el_az -1) :: vistnxx_ps_read,vistnyy_ps_read,vistnxy_ps_read ! MCNAP

	  !real*8, dimension(nnt) :: sxx_ps_read,syy_ps_read,sxy_ps_read

	  character*70 :: file_name



	  integer*4 :: id1,iaz,nprocs



  



	! "Jaky" empyrical relationship (clay and sand)

	! "Lezioni di geotecnica", R.Nova, CittàStudiEdizioni, 1995, pag.74



	if (fmaps.eq.1) then



		k0 = (1-sin(phi*0.017453292))



		do iq = 1,nn   

			do ip = 1,nn 

                is = nn*(iq -1) +ip 

                in = cs(cs(ielem -1) + is) 



				syy_ps(ip,iq,ielem) = rho * 9.81 * abs(ys(in)) 

				sxx_ps(ip,iq,ielem) = k0 * syy_ps(ip,iq,ielem)

				sxy_ps(ip,iq,ielem) = 0.0d0

				szz_ps(ip,iq,ielem) = 0.0d0  !Clara 12/01/2005 !!Controllare



			enddo 

		enddo 





	! Linear elastic approximation

	!"Lezioni di geotecnica", R.Nova, CittàStudiEdizioni, 1995, pag.71)



    elseif (fmaps.eq.2) then



	    k0 = lambda / ( 2 * (lambda + mu))



		do iq = 1,nn   

			do ip = 1,nn 

                is = nn*(iq -1) +ip 

                in = cs(cs(ielem -1) + is) 



				syy_ps(ip,iq,ielem) = rho * 9.81 * abs(ys(in)) 

				sxx_ps(ip,iq,ielem) = k0 * syy_ps(ip,iq,ielem)

				sxy_ps(ip,iq,ielem) = 0.0d0

				szz_ps(ip,iq,ielem) = 0.0d0  !Clara 12/01/2005 !!Controllare



			enddo 

		enddo 





    ! Reading the prestress - ELASTIC ONLY & withouth deformation

    elseif (fmaps.eq.3) then

		

        call read_pre_stress(nnt,file_name,nsnapps,tsnapsps,&

						    sxx_ps_read,syy_ps_read,sxy_ps_read,szz_ps_read,&

							!vistnxx_ps_read,vistnyy_ps_read,vistnxy_ps_read,& ! MCNAP

							nprocs)



	    do iq = 1,nn   

			do ip = 1,nn 

        

                is = nn*(iq -1) +ip 

                in = cs(cs(ielem -1) + is) 

        

				id1 = node_index(in) 

                iaz = update_index_el_az(id1 -1)

        

				sxx_ps(ip,iq,ielem) = sxx_ps_read(iaz) 

				syy_ps(ip,iq,ielem) = syy_ps_read(iaz) 

                sxy_ps(ip,iq,ielem) = sxy_ps_read(iaz) 

				szz_ps(ip,iq,ielem) = szz_ps_read(iaz) !Clara 12/01/2005

        

			enddo  

		enddo 





    endif

      

      

      return

      

      end subroutine make_pre_stress

      

      

!     ********************************************************







!********************************************************************

      subroutine inv_nc(nstre,vect,i1,j2,iopt)

!********************************************************************

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Clara Zambelli (28//09/2004)

!    

!

! ... routine for the evaluation of invariants of stress and strain for output purpose

!

!     nstre    = dimension of stress/strain vectors (4 for 2-d, 6 for 3-d)

!     vect = stress or strain vector defined

!                correspondence rule:

!                 2-d case     3-d case  

!                (11) <=> sxx   (11) <=> sxx 

!                (22) <=> syy   (22) <=> syy 

!                (33) <=> szz   (33) <=> szz 

!                (12) <=> sxy   (12) <=> sxy

!                               (23) <=> syz

!                               (31) <=> szx

!                and continuum mechanics sign conventions

!

!     iopt     = 1: vect => sigma

!                2: vect => epsilon

!

      data zero,one,two,three &

           /0.0d0,1.0d0,2.0d0,3.0d0/

      data half/0.5d0/

	 

	  integer*4 i,j,iopt,nstre

	  real*8, dimension(nstre):: vect

      real*8 q2,epss2,i1,j2

      real*8 temp1,temp2,temp3,temp4   

      real*8, dimension(nstre):: tempv,dev,m

      real*8, dimension (nstre,nstre):: MM,MM1



!

! ... change convention

!

!      call move_vect(nstre,vect,vect_ghm,1)

!

! ... identity tensors

!

      MM=zero

      MM1=zero

      m=zero

!

      m(1)=one

      m(2)=one

      m(3)=one

!

      MM(1,1)=one

      MM(2,2)=one

      MM(3,3)=one

      MM1(1,1)=one

      MM1(2,2)=one

      MM1(3,3)=one

!

      do i=4,nstre

!

       MM(i,i)=two

       MM1(i,i)=half

!

      end do

!

      if(iopt.eq.1) then 

!

! ... stress invariants (i1,j2)

!

       i1=(vect(1)+vect(2)+vect(3))/three

!

       dev=vect-i1*m

!

       temp1=three/two

       tempv=matmul(MM,dev)

       temp2=dot_product(dev,tempv)

       q2=temp1*temp2

!

       j2=dsqrt(q2)

!

      else

!

! ... strain invariants (epsve,epsse)

!

       i1=vect(1)+vect(2)+vect(3)

!

       dev=vect-(i1/three)*m

!

       temp1=two/three

       tempv=matmul(MM1,dev)

       temp2=dot_product(dev,tempv)

       epss2=temp1*temp2

!

       j2=dsqrt(epss2)

!

      end if

!

      end subroutine inv_nc





!!     ********************************************************



      subroutine read_pre_stress(nnode,file_name,nsnapsps,tsnapsps,&

						         sxx_ps_read,syy_ps_read,sxy_ps_read,szz_ps_read,&

								 !vistnxx_ps_read,vistnyy_ps_read,vistnxy_ps_read,& ! MCNAP

								 np)

      

      

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Clara Zambelli, Marco S.

!

     

      implicit none

      

      !integer*4 :: nnode,cs_nnz,nm,ne,opt,nsnaps,np

      !real*8, dimension(nnode) :: xs,ys

      !integer*4, dimension(0:cs_nnz) :: cs

      !integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat

      !real*8, dimension(nm) :: tref_mat

      !real*8, dimension(nm,4) :: prop_mat

      !real*8, dimension(ne) :: alfa1,alfa2,beta1,beta2,gamma1,gamma2

      !real*8, dimension(nsnaps) :: tsnaps

      !character*70 :: file_name

      

      !real*8, dimension(:), allocatable :: ux,uy,vivel_read,vistn1_0_read,vistn2_0_read,vistn3_0_read

      !real*8, dimension(:), allocatable :: sxx,syy,szz,sxy

      !real*8, dimension(:), allocatable :: spr,svm

      !integer*4 :: ip,i,j,in,nd,nn,nelem

      !integer*4 :: imat,ie,mn,n1,n2,n3,n4,count,ic

      !real*8 :: val

      !integer*4 :: istep,nsteps      

      !character*70 :: in_fileE,in_fileV,in_fileT1,in_fileT2,in_fileT3,out_file

	  !character*70 :: in_fileS1,in_fileS2,in_fileS3

      !character*40 :: elcode

      !integer*4 :: lname

      

      !real*8 ::p,q    !Clara

	  !real*8, dimension (:), allocatable::value  !Clara



      integer*4 :: i,ip,np,nd,in

      integer*4 :: nnode,lname,opt

      integer*4 :: nsnapsps

	  integer*4, dimension(nsnapsps) :: tsnapsps ! Marco

	  integer*4 :: count,id1,iaz

      real*8, dimension(nnode) :: sxx_ps_read,syy_ps_read,sxy_ps_read,szz_ps_read

	  !real*8, dimension(nnode) :: vistnxx_ps_read,vistnyy_ps_read,vistnxy_ps_read ! MCNAP

	  character*70 :: file_name,in_fileS4,in_fileS5,in_fileS6,in_fileS7

	  !character*70 :: in_fileDvp4,in_fileDvp5,in_fileDvp6 ! MCNAP

	  real*8 :: val



      lname = len_trim(file_name)

	  in_fileS4 = file_name(1:lname) // '4_xxx_xxx.bin'

	  in_fileS5 = file_name(1:lname) // '5_xxx_xxx.bin'

	  in_fileS6 = file_name(1:lname) // '6_xxx_xxx.bin'

      in_fileS7 = file_name(1:lname) // '7_xxx_xxx.bin'  !Clara 12/01/2005

	  !in_fileDvp4 = file_name(1:lname) // '7_xxx_xxx.bin' ! MCNAP

	  !in_fileDvp5 = file_name(1:lname) // '8_xxx_xxx.bin' ! MCNAP

	  !in_fileDvp6 = file_name(1:lname) // '9_xxx_xxx.bin' ! MCNAP

      !out_file = file_name(1:lname) // '.inp'







      !istep = 1

      !allocate(sxx_ps_read(nnode),syy_ps_read(nnode),sxy_ps_read(nnode))

      



      !allocate(value(4)) !Clara



      !if (opt.eq.1) then

      !   nsteps = nsnapsps

      !elseif (opt.eq.2) then

      !   nsteps = nsnapsps*2

      !elseif (opt.eq.3) then

      !   nsteps = nsnapsps*3

	  !elseif (opt.eq.7) then

      !   nsteps = nsnapsps*2

      !endif





	  count = 1

      

      !do count = 1,nsnapsps



	!	if (count.lt.10) then



			write(in_fileS4(lname+3:lname+4),'(a2)')'00'

            write(in_fileS4(lname+5:lname+5),'(i1)')count



			write(in_fileS5(lname+3:lname+4),'(a2)')'00'

            write(in_fileS5(lname+5:lname+5),'(i1)')count



			write(in_fileS6(lname+3:lname+4),'(a2)')'00'

            write(in_fileS6(lname+5:lname+5),'(i1)')count



			write(in_fileS7(lname+3:lname+4),'(a2)')'00'  !Clara 12/01/2005

            write(in_fileS7(lname+5:lname+5),'(i1)')count  !Clara 12/01/2005



			!write(in_fileDvp4(lname+3:lname+4),'(a2)')'00' ! MCNAP

            !write(in_fileDvp4(lname+5:lname+5),'(i1)')count ! MCNAP



			!write(in_fileDvp5(lname+3:lname+4),'(a2)')'00' ! MCNAP

            !write(in_fileDvp5(lname+5:lname+5),'(i1)')count ! MCNAP



			!write(in_fileDvp6(lname+3:lname+4),'(a2)')'00' ! MCNAP

            !write(in_fileDvp6(lname+5:lname+5),'(i1)')count ! MCNAP



	 !    else if (count.lt.100) then

!

!			write(in_fileS4(lname+3:lname+3),'(a1)')'0'

!            write(in_fileS4(lname+4:lname+5),'(i2)')count

!

!			write(in_fileS5(lname+3:lname+3),'(a1)')'0'

!            write(in_fileS5(lname+4:lname+5),'(i2)')count

!

!			write(in_fileS6(lname+3:lname+3),'(a1)')'0'

!            write(in_fileS6(lname+4:lname+5),'(i2)')count

!

!         else if (count.le.999) then

!

!			write(in_fileS4(lname+3:lname+5),'(i3)')count

!

!

!			write(in_fileS5(lname+3:lname+5),'(i3)')count

!

!			write(in_fileS6(lname+3:lname+5),'(i3)')count

!

!         endif

         

         do ip = 0,np -1

 !           if (ip.lt.10) then



			   write(in_fileS4(lname+7:lname+8),'(a2)')'00'

               write(in_fileS4(lname+9:lname+9),'(i1)')ip



			   write(in_fileS5(lname+7:lname+8),'(a2)')'00'

               write(in_fileS5(lname+9:lname+9),'(i1)')ip



			   write(in_fileS6(lname+7:lname+8),'(a2)')'00'

               write(in_fileS6(lname+9:lname+9),'(i1)')ip



			   write(in_fileS7(lname+7:lname+8),'(a2)')'00'  !Clara 12/01/2005

               write(in_fileS7(lname+9:lname+9),'(i1)')ip    !Clara 12/01/2005



			   !write(in_fileDvp4(lname+7:lname+8),'(a2)')'00' ! MCNAP

               !write(in_fileDvp4(lname+9:lname+9),'(i1)')ip ! MCNAP



			   !write(in_fileDvp5(lname+7:lname+8),'(a2)')'00' ! MCNAP

               !write(in_fileDvp5(lname+9:lname+9),'(i1)')ip ! MCNAP



			   !write(in_fileDvp6(lname+7:lname+8),'(a2)')'00' ! MCNAP

               !write(in_fileDvp6(lname+9:lname+9),'(i1)')ip ! MCNAP







!			else if (ip.lt.100) then

!

!			   write(in_fileS4(lname+7:lname+7),'(a1)')'0'

!               write(in_fileS4(lname+8:lname+9),'(i2)')ip

!

!			   write(in_fileS5(lname+7:lname+7),'(a1)')'0'

!               write(in_fileS5(lname+8:lname+9),'(i2)')ip

!

!			   write(in_fileS6(lname+7:lname+7),'(a1)')'0'

!               write(in_fileS6(lname+8:lname+9),'(i2)')ip

!

!            else if (ip.le.999) then

!

!				write(in_fileS4(lname+7:lname+9),'(i3)')ip

!				write(in_fileS5(lname+7:lname+9),'(i3)')ip

!				write(in_fileS6(lname+7:lname+9),'(i3)')ip

!

!            endif





			open(20,file=in_fileS4)

            

            read(20,*)nd

            

            do i = 1,nnode

               read(20,*)in,val

               if (dabs(val).lt.1.d-30) val=0.d0

               sxx_ps_read(i) = val



            enddo

            

            close(20)



			open(20,file=in_fileS5)

            

            read(20,*)nd

            

            do i = 1,nnode

               read(20,*)in,val

               if (dabs(val).lt.1.d-30) val=0.d0

               syy_ps_read(i) = val



            enddo

            

            close(20)



			open(20,file=in_fileS6)

            

            read(20,*)nd

            

            do i = 1,nnode

               read(20,*)in,val

               if (dabs(val).lt.1.d-30) val=0.d0

               sxy_ps_read(i) = val



            enddo

            

            close(20)



            !begin  Clara 12/01/2005

			open(20,file=in_fileS7)

            

            read(20,*)nd

            

            do i = 1,nnode

               read(20,*)in,val

               if (dabs(val).lt.1.d-30) val=0.d0

               szz_ps_read(i) = val



            enddo

            

            close(20)



             !end  Clara 12/01/2005

            

			!open(20,file=in_fileDvp4) ! MCNAP

            

            !read(20,*)nd ! MCNAP

            

            !do i = 1,nnode ! MCNAP

               !read(20,*)in,val ! MCNAP

               !if (dabs(val).lt.1.d-30) val=0.d0 ! MCNAP

               !vistnxx_ps_read(i) = val ! MCNAP



            !enddo ! MCNAP

            

            !close(20) ! MCNAP





			!open(20,file=in_fileDvp5) ! MCNAP

            

            !read(20,*)nd ! MCNAP

            

            !do i = 1,nnode ! MCNAP

               !read(20,*)in,val ! MCNAP

               !if (dabs(val).lt.1.d-30) val=0.d0 ! MCNAP

               !vistnyy_ps_read(i) = val ! MCNAP



            !enddo ! MCNAP

            

            !close(20) ! MCNAP





			!open(20,file=in_fileDvp6) ! MCNAP

            

            !read(20,*)nd ! MCNAP

            

            !do i = 1,nnode ! MCNAP

               !read(20,*)in,val ! MCNAP

               !if (dabs(val).lt.1.d-30) val=0.d0 ! MCNAP

               !vistnxy_ps_read(i) = val ! MCNAP



            !enddo ! MCNAP

            

            !close(20) ! MCNAP





		enddo

         

         

!      enddo

      

      

      !deallocate(ux,uy,vivel_read)

      !deallocate(sxx,syy,szz,sxy)

      !deallocate(spr,svm)

	  !deallocate(duxdx_average,duydx_average,duxdy_average,duydy_average)



	  !deallocate(value)  !Clara

      

      return

      

      end subroutine read_pre_stress

      



!********************************************************************

      subroutine LEG(s,ak,dep,D,DV,dt,csi,bef,pc,u)

!********************************************************************

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Clara Zambelli (28//10/2004)





      implicit none



	  real*8 :: z6,csi,bef,pc

      real*8 :: J3E,J2E,J7,P1,SQ3

	  real*8, dimension(6) ::s,dep

	  real*8, dimension(3,3) :: ak,D,DV

      real*8, dimension(3,3) :: S1,ET,SET,DPP,SS

	  integer*4 :: i,j,k

      

	  real*8, dimension(14) ::u 

	  real*8 :: gam,gav,ga,z51,ff,dt



      gam=u(13)   

      gav=u(14)   



      call TEN(s,S1)



      do i = 1,3

	    do j = 1,3

	   	   DPP(i,j)=0.0d0

 !          S1(i,j)=S1(i,j)-UA*D(i,j)  !don't use because there isn't water -> ua=0

        enddo

      enddo

	  	 

  	!VERIFICARE modifica perche' il vettore s1 e' tutto nullo pongo s1(1,1)=1   !ClaraZ&begin

 	do i=1,3

      do j=1,3

   		if (s1(i,j).eq.0) then

			s1(1,1) = 1.0d0

		endif

	  enddo

	enddo                                                !ClaraZ&end



      SQ3=SQRT(3.)

      ga=u(1)   !cambiare



      call SCAL(S1,ak,P1)



      do i = 1,3

	    do j = 1,3

	   	  ET(i,j)=(S1(i,j)/P1-AK(i,j))*SQ3

        enddo

      enddo

	  

	  do i = 1,3

	    do j = 1,3

	   	 SET(i,j)=0.0d0

		 do k = 1,3

		   SET(i,j)=SET(i,j)+ET(i,k)*ET(k,j)

         enddo

        enddo

      enddo

	  

	  call SCAL(ET,ET,J2E)

      call SCAL(SET,ET,J3E)

      call SCAL(SET,AK,J7)     

     

!      FORMATION OF OMEGA 1



      Z51=(9.*(ga-3.)+3.*ga*J3E+3.*SQ3*ga*J7-9.*(ga-1.)*J2E/2.)/P1



!     calcolus of plastic potential



      do i = 1,3

	    do j = 1,3

	   	  SS(i,j)=Z51*AK(I,J)-3.*ga*SQ3/P1*SET(I,J)+9./2.*(ga-1.)*SQ3/P1*ET(I,J)

        enddo

      enddo

	 

      call LIM(U,AK,S1,DV,D,csi)

!     DEFINITION OF TENSORIAL COMPLIANCES





      FF=(P1*EXP((9.*(GA-1.)*J2E/4.-GA*J3E)/(3.*BEF*(GA-3.)))-PC)/pc



      FF=1.0d0

!      FF=(3.*BEF*(GA-3.))*log(P1/pc)+ 9.*(GA-1.)*J2E/4.-GA*J3E

!      write(*,*)ff,p1,j2e,j3e

      

	  do i = 1,3

	    do j = 1,3

	   	  DPP(i,j)=SS(i,j)*gam*p1*exp(gav*FF)*DT

        enddo

      enddo

     

      do i = 1,3

	   	  DEP(i)=DPP(i,i)

      enddo

      do i = 1,2

	   	  DEP(i+3)=DPP(1,i+1)*2.

      enddo

      DEP(6)=DPP(2,3)*2.

   

   !   write(*,*)dep



     return



     end subroutine LEG	   







!********************************************************************

      subroutine LIM(U,AK,S1,DV,D,csi)

!********************************************************************

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Clara Zambelli (28//10/2004)



      

      implicit none



	  real*8 :: csi 

      real*8 :: B4,RTT,RN,Z7

	  real*8, dimension(3,3) :: AW,DV,R,AK,st2

      real*8, dimension(3,3) :: SV0,S1,RV,ST,ST1,D

	  integer*4 :: i,j



      real*8, dimension(14) :: u 

	  real*8 :: cs1,cs2,PQ,BM,TOL,SQ3

	  real*8 :: R0,ST0,sq1,sq2,B7,B5,B9,PT4,PT5

      



      cs1=u(4)  

      cs2=u(12)  

      PQ=u(8)    

      BM=u(9)    

      TOL=.0000001

      SQ3=SQRT(3.)

	  do i = 1,3

	    do j = 1,3

	   	   st2(i,j)=0.0d0

           ST1(i,j)=0.0d0 

        enddo

      enddo



      ST1(1,1)=.8164965809

      ST1(2,2)=-.4082482905

      ST1(3,3)=-.4082482905





      call SCAL(S1,S1,Z7)

      Z7=SQRT(Z7)



      do i = 1,3

	    do j = 1,3

	   	   SV0(i,j)=S1(i,j)/Z7

        enddo

      enddo

	



      R0=0.0d0

      ST0=0.0d0

      do i = 1,3

	    do j = 1,3

	   	   R0=R0+AK(i,i)

           ST0=ST0+SV0(i,i)

        enddo

      enddo



      do i = 1,3

	    do j = 1,3

	   	   RV(i,j)=AK(i,j)-R0*D(i,j)/3.

           ST(j,j)=SV0(i,j)-ST0*D(i,j)/3.

        enddo

      enddo

	 

      do i = 1,3

	    do j = 1,3

	   	   R(i,j)=ST(i,j)-RV(i,j) 

        enddo

      enddo

	

      call SCAL(R,R,RN)

      call SCAL(RV,RV,RTT)

      RN=SQRT(RN)



      if(RN.lt.TOL)then

	    RN=TOL

	  endif

	  do i = 1,3

	    do j = 1,3

	   	   R(I,J)=R(I,J)/RN

        enddo

      enddo

	  

      sq1=(r(1,1)+r(2,2))/2.

      sq2=sqrt(((r(1,1)-r(2,2))/2.)**2+r(1,2)**2)

      st2(1,1)=sq1+sq2

      st2(2,2)=sq1-sq2

      st2(3,3)=r(3,3)



      call SCAL(st2,ST1,B4)



      B7=ABS(B4)

      if(B7.gt.1.)B7=1.

      if(B4.gt.1.)B4=1.

      if(B4.lt.0.)B4=-B7

      B5=ACOS(B4)



	  if(B5.ge.1.046666)then

          if(B5.ge.3.14)then

		     if (B5.lt.5.2333332)B5=B5-4.1866664

		  else 

            B5=B5-2.0933332

		  endif

      endif

     

	  B5=ABS(B5)



      B9=(pq-bm)/2.*(2.*b5**2-4.*b5+1.)+(pq+bm)/2.

      csi=(cs1-cs2)/2.*(2.*b5**2-4.*b5+1.)+(cs1+cs2)/2.



      PT4=COS(B9)

      PT5=SIN(B9)



     do i = 1,3

	    do j = 1,3

	   	   AW(i,j)=PT4*D(i,j)/SQ3+PT5*R(i,j)

           DV(i,j)=AW(i,j)-AK(i,j)

        enddo

      enddo

	 



    !  call SCAL(AK,DV,HP)

    ! PW=U(10)  !cambiare

	!  do i = 1,3

	!    do j = 1,3

	!  	   DEL(i,j)=(DV(i,j)-AK(i,j)*HP)*PW/Z6 

    !    enddo

    !  enddo

    

	  return



      end subroutine LIM		   



!********************************************************************

      subroutine UPD(a,ep,ak,dep,d,dv,csi,z14,bef,pc,u)

!********************************************************************

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Clara Zambelli (28//10/2004)

!

      

      implicit none

      

	  real*8 :: z6,csi,z14,bef,pc

	  real*8 :: depv

      real*8 :: z13,ell,vx,epv,pp,vz

	  real*8 :: tol1,pcd

	  real*8, dimension(6) :: ep,dep

	  real*8, dimension(3,3) :: ak,d,dv

      real*8, dimension(3,3) :: dep5,ep7,ep8,dep3,dep4,a

	  integer*4 :: i,j



      real*8, dimension(14) :: u 

    

      tol1=.0000001



 !     gw=U(5)    !!!cambiare e controllare se serve

 !     ww=1./3.    !controllare se serve

 !     vz=tol1   ! a che serve??



!

!     UPDATING OF STRESSES AND STRAINS

!

      do i = 1,3

	   ep(i)=ep(i)+dep(i)

      enddo

    

!

!     UPDATING A

!

      call ten(dep,dep3)

      call ten(ep,ep7)



      do i = 1,3

	    do j = 1,3

		   if (i.ne.j) then

	        dep3(i,j)=dep3(i,j)/2.

            ep7(i,j)=ep7(i,j)/2. 

           endif

        enddo

      enddo

	 



      call SCAL(ep7,d,epv)

	  epv=epv/3.

	  call SCAL(dep3,ak,vx)

	  do i = 1,3

	    do j = 1,3

	   	   ep8(i,j)=ep7(i,j)-epv*d(i,j)

           dep5(i,j)=dep3(i,j)-vx*ak(i,j)  

        enddo

      enddo

     

      call SCAL(ep8,ep8,pp)  ! PIANTA!!!

      pp=sqrt(pp)





      call SCAL(dep5,dep5,vz)

      vz=sqrt(vz)



 

      depv=0.0d0

	  do i = 1,3

	     depv=dep3(i,i)/3.+depv

      enddo



      

      call SCAL(dep3,dep3,ell)

      ell=sqrt(ell)



      pcd=pc*((vx+csi*vz)/U(2))   !!cambiare u(2)

      pc=pc+pcd



      do i = 1,3

	    do j = 1,3

	   	   dep4(i,j)=dep3(i,j)-depv*d(i,j)

           a(i,j)=a(i,j)+dv(i,j)*ell*U(10)   !cambiare u(10)

        enddo

      enddo



      call SCAL(a,a,z6)

	  call SCAL(dep4,ep8,z13)

      if(pp.lt.tol1)then

	    pp=tol1

	  endif



      z13=z13/pp



      z6=sqrt(z6)

      z14=z14+z13

      bef=u(11)+(u(3)-u(11))*exp(-z14*u(5))    !cambiare le varie U



      do i = 1,3

	    do j = 1,3

	   	    ak(i,j)=a(i,j)/z6

        enddo

      enddo

	 



      return



      end subroutine UPD





!********************************************************************

      subroutine SCAL(A,B,C)

!********************************************************************

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Clara Zambelli (28//10/2004)

!

! Scalar product of two arrays A(n,n) , B(n,n) ----> C





      implicit none



      real*8 :: C

	  real*8, dimension(3,3) :: A,B 

	  integer*4 :: i,j

          



      C=0.0d0



	  do i = 1,3

	   do j = 1,3

		  C=C+A(i,j)*B(j,i)  

       enddo

      enddo



	  return



      end subroutine SCAL







!********************************************************************

      subroutine TEN(A,B)

!********************************************************************

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Clara Zambelli (28//10/2004)

!

! from vector A(2n) to array B(n,n)

 

      implicit none

      

	  real*8, dimension(6) :: A

	  real*8, dimension(3,3) :: B 

	  integer*4 :: i,j,ni



      do i = 1,3

	    do j = 1,3

		   ni=i+j+1

		   if (i.eq.j) then

	         B(i,i)=A(i)

	       else

	         B(i,j)=A(ni)

           endif

        enddo

      enddo

      

      return



      end subroutine TEN





!********************************************************************

   	subroutine PAR(dens,u,&

		 ga_l,la_l,bef0_l,csi_l,psi_l,gw_l,po_l,pq_l,pw_l,bm_l,befl_l,gam_l,gav_l,&

		 ga_d,la_d,bef0_d,csi_d,psi_d,gw_d,po_d,pq_d,pw_d,bm_d,befl_d,gam_d,gav_d)

!********************************************************************

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Clara Zambelli (28//10/2004)

!

! 

      implicit none

      

	  real*8 :: dens

      real*8 :: densr1

      real*8 :: ga_l,la_l,bef0_l,csi_l,psi_l,gw_l,po_l,pq_l,pw_l,bm_l,befl_l,gam_l,gav_l

	  real*8 ::	ga_d,la_d,bef0_d,csi_d,psi_d,gw_d,po_d,pq_d,pw_d,bm_d,befl_d,gam_d,gav_d



	  real*8 :: ga,la,bef0,csi,psi,gw,po,pq,pw,bm,befl,gam,gav



	  real*8, dimension(14) :: u

 

          

      densr1=dens-20.



	  if(dens.le.20.)densr1=0.

      if(dens.ge.100.)densr1=80.

      

          

      ga=ga_l+(ga_d-ga_l)*(densr1/80.) 

      la=la_l+(la_d-la_l)*(densr1/80.)   

      bef0=bef0_l+(bef0_d-bef0_l)*(densr1/80.) 

      csi=csi_l+(csi_d-csi_l)*(densr1/80.) 

      psi=psi_l+(psi_d-psi_l)*(densr1/80.) 

      gw=gw_l+(gw_d-gw_l)*(densr1/80.) 

      po=po_l+(po_d-po_l)*(densr1/80.) 

      pq=pq_l+(pq_d-pq_l)*(densr1/80.) 

      pw=pw_l+(pw_d-pw_l)*(densr1/80.) 

      bm=bm_l+(bm_d-bm_l)*(densr1/80.) 

      befl=befl_l+(befl_d-befl_l)*(densr1/80.) 

      gam=gam_l+(gam_d-gam_l)*(densr1/80.)                   

      gav=gav_l+(gav_d-gav_l)*(densr1/80.)   



!    identificazione dei parametri

      u(1)=GA

      u(2)=LA

      u(3)=BEF0

      u(4)=CSI

      u(5)=GW

      u(6)=PO

      u(8)=PQ

      u(9)=BM

      u(10)=PW

      u(11)=befl

      u(12)=psi

      u(13)=GAM

      u(14)=GAV



      return



      end subroutine PAR



!********************************************************************

      subroutine DENSITA(dens,dep,emax,emin)

!********************************************************************

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Clara Zambelli (28//10/2004)

!

      implicit none

      

	  real*8 :: dens,emax,emin,dvol

      real*8 :: rrr,ecor,einc,efin

	  real*8, dimension(6) :: dep

	 



	  rrr=(emax-emin)/100

      ecor=emax-rrr*dens



      dvol=dep(1)+dep(2)+dep(3)



      einc=-(1+ecor)*dvol

      efin=ecor+einc

      dens=(emax-efin)/rrr

      

      return



      end subroutine DENSITA





!!     ********************************************************



      subroutine read_initial_condition(nnode,file_name,nsnapsps,tsnapsps,&

								 u1_read,u0_read,&

								 Fe_read,Fk_read,&

								 vivel_read,&

								 vistn1_0_read,vistn2_0_read,vistn3_0_read,&

								 !vistnxx_ps_read,vistnyy_ps_read,vistnxy_ps_read,& ! MCNAP

								 np)

      

      

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Marco S.

!     Mexico _A_

     

      implicit none

      

      !integer*4 :: nnode,cs_nnz,nm,ne,opt,nsnaps,np

      !real*8, dimension(nnode) :: xs,ys

      !integer*4, dimension(0:cs_nnz) :: cs

      !integer*4, dimension(nm) :: tag_mat,type_mat,sdeg_mat

      !real*8, dimension(nm) :: tref_mat

      !real*8, dimension(nm,4) :: prop_mat

      !real*8, dimension(ne) :: alfa1,alfa2,beta1,beta2,gamma1,gamma2

      !real*8, dimension(nsnaps) :: tsnaps

      !character*70 :: file_name

      

      !real*8, dimension(:), allocatable :: ux,uy,vivel_read,vistn1_0_read,vistn2_0_read,vistn3_0_read

      !real*8, dimension(:), allocatable :: sxx,syy,szz,sxy

      !real*8, dimension(:), allocatable :: spr,svm

      !integer*4 :: ip,i,j,in,nd,nn,nelem

      !integer*4 :: imat,ie,mn,n1,n2,n3,n4,count,ic

      !real*8 :: val

      !integer*4 :: istep,nsteps      

      !character*70 :: in_fileE,in_fileV,in_fileT1,in_fileT2,in_fileT3,out_file

	  !character*70 :: in_fileS1,in_fileS2,in_fileS3

      !character*40 :: elcode

      !integer*4 :: lname

      

      !real*8 ::p,q    !Clara

	  !real*8, dimension (:), allocatable::value  !Clara



      integer*4 :: i,ip,np,nd,in

      integer*4 :: nnode,lname,opt

      integer*4 :: nsnapsps

	  integer*4, dimension(nsnapsps) :: tsnapsps ! Marco

	  integer*4 :: count,id1,iaz

	  real*8, dimension(nnode*2) :: u1_read,u0_read

	  real*8, dimension(nnode*2) :: Fe_read,Fk_read

	  real*8, dimension(nnode*2) :: vivel_read

	  real*8, dimension(nnode*2) :: vistn1_0_read,vistn2_0_read,vistn3_0_read

	  character*70 :: file_name

	  character*70 :: in_file_u1,in_file_u0

	  character*70 :: in_file_Fe,in_file_Fk

	  character*70 :: in_file_vivel

	  character*70 :: in_file_vistn1_0,in_file_vistn2_0,in_file_vistn3_0

	  real*8 :: val



      lname = len_trim(file_name)

	  in_file_u1 = file_name(1:lname) // '_u1.bin'

	  in_file_u0 = file_name(1:lname) // '_u0.bin'

	  in_file_Fe = file_name(1:lname) // '_Fe.bin'

	  in_file_Fk = file_name(1:lname) // '_Fk.bin'

	  in_file_vivel = file_name(1:lname) // '_vivel.bin'

	  in_file_vistn1_0 = file_name(1:lname) // '_vistn1_0.bin'

	  in_file_vistn2_0 = file_name(1:lname) // '_vistn2_0.bin'

	  in_file_vistn3_0 = file_name(1:lname) // '_vistn3_0.bin'





	  open(20,file=in_file_u1)           

      read(20,*)nd           

      do i = 1,nd

		read(20,*)in,val

			if (dabs(val).lt.1.d-30) val=0.d0

		    u1_read(i) = val   

	  enddo ! i        

      close(20)



	  open(20,file=in_file_u0)           

      read(20,*)nd           

      do i = 1,nd

		read(20,*)in,val

			if (dabs(val).lt.1.d-30) val=0.d0

		    u0_read(i) = val   

	  enddo ! i        

      close(20)



	  open(20,file=in_file_Fe)           

      read(20,*)nd           

      do i = 1,nd

		read(20,*)in,val

			if (dabs(val).lt.1.d-30) val=0.d0

		    Fe_read(i) = val   

	  enddo ! i        

      close(20)



	  open(20,file=in_file_Fk)           

      read(20,*)nd           

      do i = 1,nd

		read(20,*)in,val

			if (dabs(val).lt.1.d-30) val=0.d0

		    Fk_read(i) = val   

	  enddo ! i        

      close(20)



	  open(20,file=in_file_vivel)           

      read(20,*)nd           

      do i = 1,nd

		read(20,*)in,val

			if (dabs(val).lt.1.d-30) val=0.d0

		    vivel_read(i) = val   

	  enddo ! i        

      close(20)



	  open(20,file=in_file_vistn1_0)           

      read(20,*)nd           

      do i = 1,nd

		read(20,*)in,val

			if (dabs(val).lt.1.d-30) val=0.d0

		    vistn1_0_read(i) = val   

	  enddo ! i        

      close(20)



	  open(20,file=in_file_vistn2_0)           

      read(20,*)nd           

      do i = 1,nd

		read(20,*)in,val

			if (dabs(val).lt.1.d-30) val=0.d0

		    vistn2_0_read(i) = val   

	  enddo ! i        

      close(20)



	  open(20,file=in_file_vistn3_0)           

      read(20,*)nd           

      do i = 1,nd

		read(20,*)in,val

			if (dabs(val).lt.1.d-30) val=0.d0

		    vistn3_0_read(i) = val   

	  enddo ! i        

      close(20)





      

      return

      

      end subroutine read_initial_condition



!!     ********************************************************

      

      subroutine write_bin_file_new(file_name,count,proc,nu,ui,uv)

      

!     © POLIMI, 2004, All Rights Reserved

!     Authors: Marco Stupazzini

!     

      integer*4 :: count,proc,nu

      real*8, dimension(*) :: uv

      integer*4, dimension(*) :: ui

      character*70 :: file_name



      character*70 :: out_file

      integer*4 :: i,lname

      

      lname = len_trim(file_name)

      out_file = file_name(1:lname) // '.bin'

      

      open(20,file=out_file)

      

      write(20,*)nu

      

      do i = 1,nu

         write(20,*)ui(i),uv(i)

      enddo

      

      close(20)

      

      return

      

      end subroutine write_bin_file_new

      

      

      

      

!     ********************************************************

