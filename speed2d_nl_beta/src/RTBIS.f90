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

!> @brief Find a root of a polynomial of degree n.
!! @author Ilario Mazzieri
!> @date April, 2014
!> @version 1.0

!> @param[in] x1 left bound of the interval [x1,x2]
!> @param[in] x2 right bound of the interval [x1,x2]
!> @param[in] xacc  acceleration parameter
!> @param[in] n  polynomial degree
!> @param[out] RTBIS  zero of the polynomial function of order n 

      real*8 function RTBIS(x1,x2,xacc,n)

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

      end function RTBIS

