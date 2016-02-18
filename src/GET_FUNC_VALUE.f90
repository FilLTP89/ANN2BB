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

!> @brief Computes time evolution function.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] nf number of functions
!> @param[in] func_type function type
!> @param[in] func_indx indices for the data 
!> @param[in] func_data data for the calculation (depending on type_fnc) 
!> @param[in] n  number of the function
!> @param[in] t  instant time
!> @param[in] t0_sism time delay
!> @param[out] GET_FUNC_VALUE value of the time function

      real*8 function GET_FUNC_VALUE(nf,func_type,func_indx,func_data,n,t,t0_sism)
     
      implicit none
      
      integer*4 :: nf,n
      integer*4, dimension(nf) :: func_type
      integer*4, dimension(nf +1) :: func_indx
      real*8, dimension(*) :: func_data
      real*8 :: t, beta2, t0_sism
      
      integer*4 :: i
      real*8 :: val,PI,tt_t0,t0,t1,v0,v1
      real*8 :: ts, tp, AMP
      
      
      val = 0.0d0
      
      if (func_type(n).eq.0) then
         val = 1.0d0

      elseif (func_type(n) .eq. 100) then
         PI = 4.0d0 * datan(1.0d0)
         val = dsin(dsqrt(2.d0) * PI * t)

      ! - Ricker "beta" type
      elseif (func_type(n).eq.1) then
         tt_t0 = t - func_data(func_indx(n) +1)
         val = (1.0d0 - 2.0d0*func_data(func_indx(n))*tt_t0*tt_t0) &
              * dexp(-1.0d0*func_data(func_indx(n))*tt_t0*tt_t0)
      
      ! - Ricker "cos" type
      elseif (func_type(n).eq.2) then
         PI = 4.0d0 * datan(1.0d0)
         tt_t0 = t - func_data(func_indx(n) +1)
         val = dcos(PI*func_data(func_indx(n))*tt_t0) &
              * dexp(-0.5d0*func_data(func_indx(n)) &
              * func_data(func_indx(n))*tt_t0*tt_t0)

      ! - Force history
      elseif (func_type(n).eq.3) then
         do i = func_indx(n),func_indx(n+1) -3,2  !Mette -3 per eliminare l'ultimo step che viene considerato
            t0 = func_data(i)                     !prendendo  (i +2) e (i +3) 
            t1 = func_data(i +2)
            v0 = func_data(i +1)
            v1 = func_data(i +3)
            if ((t.ge.t0).and.(t.le.t1)) then
               val = (v1 - v0) / (t1 - t0) * (t - t0)  + v0
            endif
         enddo
      
      ! - First derivative of the Ricker ( d(Ricker(t))/dt )
      elseif (func_type(n).eq.4) then
         tt_t0 = t - func_data(func_indx(n) +1)
         beta2=func_data(func_indx(n))
         val = 2.0d0*beta2*tt_t0 &
              *(-3.0d0 + 2.0d0*beta2*tt_t0*tt_t0) &
              * dexp(-beta2*tt_t0*tt_t0)

      ! - Ricker "beta" type for seismic moment
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

           ! 31 - INTEGRAL OF FUNCTION 23:
	  !      Triangular function in order to test with PACO analytical sol (prg. "Gradelan")
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

      elseif (func_type(n).eq.31) then

	     tp = func_data(func_indx(n))	 !0.30d0
         ts = func_data(func_indx(n) +1) ! 0.5d0
		 AMP = func_data(func_indx(n) +2)! 1.0d0

         tt_t0 = t - tp - t0_sism
		 !tt_t0 = t - func_data(func_indx(n) +1) - t0_sism


		 if (tt_t0.gt.ts/2) then

			val = -AMP/ts * (ts/2)**2 + AMP * (ts/2)+ 1/4*AMP*ts;

	     elseif ((tt_t0.gt.0).and.(tt_t0.le.ts/2)) then
			
			val = -AMP/ts * (tt_t0)**2 + AMP * (tt_t0)+ 1/4*AMP*ts;
          
         elseif ((tt_t0.le.0).and.(tt_t0.gt.-ts/2)) then
			
			val = AMP/ts * (tt_t0)**2 + AMP * (tt_t0) + 1/4*AMP*ts;
         else
		    
            val = 0.0d0

	     endif

!-----------------------------------------------------
	  ! FUNCTION FOR G/G0
	  elseif (func_type(n).eq.60) then              ! Non-Linear Elasticity 03.10.2006
		do i = func_indx(n),func_indx(n+1) -3,2     ! Non-Linear Elasticity 03.10.2006
			t0 = func_data(i)                       ! Non-Linear Elasticity 03.10.2006
            t1 = func_data(i +2)                    ! Non-Linear Elasticity 03.10.2006
            v0 = func_data(i +1)                    ! Non-Linear Elasticity 03.10.2006
            v1 = func_data(i +3)                    ! Non-Linear Elasticity 03.10.2006
			if (abs(t).le.func_data(func_indx(n))) then           ! Non-Linear Elasticity 03.10.2006
				val = func_data(func_indx(n)+1)                   ! Non-Linear Elasticity 03.10.2006
            elseif ((abs(t).ge.t0).and.(abs(t).le.t1)) then       ! Non-Linear Elasticity 03.10.2006               
				val = (v1 - v0) / (t1 - t0) * (abs(t) - t0)  + v0 ! Non-Linear Elasticity 03.10.2006
			elseif (abs(t).ge.func_data(func_indx(n+1)-2)) then   ! Non-Linear Elasticity 03.10.2006                                              ! Non-Linear Elasticity 03.10.2006
				val = func_data(func_indx(n+1)-1)                 ! Non-Linear Elasticity 03.10.2006
            endif                                                 ! Non-Linear Elasticity 03.10.2006
		enddo
!-----------------------------------------------------
	  ! FUNCTION FOR DAMPING
	  elseif (func_type(n).eq.61) then              ! Non-Linear Elasticity 03.10.2006
		do i = func_indx(n),func_indx(n+1) -3,2     ! Non-Linear Elasticity 03.10.2006
			t0 = func_data(i)                       ! Non-Linear Elasticity 03.10.2006
            t1 = func_data(i +2)                    ! Non-Linear Elasticity 03.10.2006
            v0 = func_data(i +1)                    ! Non-Linear Elasticity 03.10.2006
            v1 = func_data(i +3)                    ! Non-Linear Elasticity 03.10.2006
			if (abs(t).le.func_data(func_indx(n))) then           ! Non-Linear Elasticity 03.10.2006
				val = func_data(func_indx(n)+1)                   ! Non-Linear Elasticity 03.10.2006
            elseif ((abs(t).ge.t0).and.(abs(t).le.t1)) then       ! Non-Linear Elasticity 03.10.2006               
				val = (v1 - v0) / (t1 - t0) * (abs(t) - t0)  + v0 ! Non-Linear Elasticity 03.10.2006
			elseif (abs(t).ge.func_data(func_indx(n+1)-2)) then   ! Non-Linear Elasticity 03.10.2006                                              ! Non-Linear Elasticity 03.10.2006
				val = func_data(func_indx(n+1)-1)                 ! Non-Linear Elasticity 03.10.2006
            endif                                                 ! Non-Linear Elasticity 03.10.2006
		enddo
      endif
      
      get_func_value = val
      
      return
      
      end function GET_FUNC_VALUE

