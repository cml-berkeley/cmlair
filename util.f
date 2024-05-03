c=======================================
      real(8) function fit1(x)
c=======================================
      common/var/s,sc
      
      real(8) :: x, s, sc
      
      fit1=(x-s)*exp(-0.5*x*x)

      return
      end
c=======================================
      real(8) function fitb1(x)
c=======================================
      
      common/var/s,sc
      
      real(8) :: x, s, sc
      
      fitb1=(2.d0*(x-s)-sc)*exp(-0.5*x*x)

      return
      end
c========================================
      real(8) function fit32(x)
c========================================
      
      common/var/s,sc
      
      real(8) :: x, s, sc

c ****sometimes fails when elastic-plastic model is chosen
      fit32=(x-s)*dsqrt(x-s)*exp(-0.5*x*x)
      return
      end
c==========================================
      subroutine QROMB(func,a,b,ss)
c==========================================
      implicit none
      
      real(8) :: a, b, dss, ss
      
      external func
      real(8), parameter :: eps=1.d-6
      
      
      integer, parameter :: jmax=40, jmaxp=jmax+1, k=5, km=k-1
      integer :: j
      real(8) :: s(jmaxp), h(jmaxp)
      
      
      h(1)=1.d0
      do j=1,jmax
	  call trapzd(func,a,b,s(j),j)
	  if(j.ge.k) then
	    call polint(h(j-km),s(j-km),k,0.d0,ss,dss)
c	    write (*,*) 'qromb iteration: ', j, 'error: ', dabs(dss)
	    if(dabs(dss).lt.eps*dabs(ss)) return
	  endif
	  s(j+1)=s(j)
	  h(j+1)=0.25d0*h(j)
      enddo
	write(*,*) 'too many steps in qromb.  Exiting'
	write(*,*) 'this error was caused by the'
	write(*,*) 'Elastic-Plastic contact model'
c	Changed from pause to stop by BC so it won't crash the interface
	stop
c      pause 'too many steps in qromb'

      return
      end
c==========================================
      subroutine trapzd(func,a,b,s,n)
c==========================================
      implicit none
      
      real(8) :: a, b, s, tnm, del, x, sum
      real(8) :: func
      integer :: j, n, it
      
      external func
      
      if(n.eq.1) then
	s=0.5d0*(b-a)*(func(a)+func(b))
	it=1
      else
	tnm=it
	del=(b-a)/tnm
	x=a+0.5d0*del
	sum=0.d0
	do j=1,it
	  sum=sum+func(x)
	  x=x+del
	enddo
	s=0.5d0*(s+(b-a)*sum/tnm)
	it=2*it
      endif

      return
      end
      
c==========================================
      subroutine polint(xa,ya,n,x,y,dy)
c==========================================
      implicit none
      
      real(8) :: xa,ya,x,y,dy
      integer :: i, m, nmax, n, ns
      real(8) :: dif, dift, c, d, ho, hp, w, den
      
      parameter (nmax=10)
      dimension xa(n),ya(n),c(nmax),d(nmax)
      ns=1
      dif=dabs(x-xa(1))
      do i=1,n
	  dift=dabs(x-xa(i))
	  if(dift.lt.dif) then
	    ns=i
	    dif=dift
	  endif
	  c(i)=ya(i)
	  d(i)=ya(i)
      enddo
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
	  do i=1,n-m
	    ho=xa(i)-x
	    hp=xa(i+m)-x
	    w=c(i+1)-d(i)
	    den=ho-hp
	    if(den.eq.0.d0) stop
	    den=w/den
	    d(i)=hp*den
	    c(i)=ho*den
	  enddo
	  if(2*ns.lt.n-m) then
	    dy=c(ns+1)
	  else
	    dy=d(ns)
	    ns=ns-1
	  endif
	  y=y+dy
      enddo

      return
      end

c     two dimensional polynomial interpolation as implemented
c     in numerical recipies in fortran/c
c==========================================================
c      subroutine polint2(x1a, x2a, ya, m, n, x1, x2, y, dy)
c==========================================================
c      INTEGER m, n, NMAX, MMAX
c	REAL dy, x1, x2, y, x1a(m), x2a(n), ya(m,n)
c	PARAMETER (NMAX=20, MMAX=20)
c
c	INTEGER j, k
c	REAL ymtmp(MMAX), yntmp(NMAX)
c
c     do j=1, m
c	  do k=1, n
c	    yntmp(k)=ya(j,k)
c	  enddo
c	  call polint(x2a, yntmp, n, x2, ymtmp(j), dy)
c	enddo
c      call polint(x1a, ymtmp, m, x1, y, dy)
c 	return
c      END

c=================================================
      subroutine finteg(tgs,fgs1,fgs2,n,x,y1,y2)
c=================================================
      implicit none
      
      real(8) :: tgs,fgs1,fgs2,x,y1,y2
      real(8) :: a, b
      integer :: n, j
      
      dimension tgs(n),fgs1(n),fgs2(n)

      call locate(tgs,n,x,j)
      if(tgs(j).eq.tgs(j+1)) write(*,*)'In finteg, j = ',j
      a=(x-tgs(j+1))/(tgs(j)-tgs(j+1))
      b=(x-tgs(j))/(tgs(j+1)-tgs(j))


      y1=a*fgs1(j)+b*fgs1(j+1)
      y2=a*fgs2(j)+b*fgs2(j+1)

      return
      end
      
c=================================================
      subroutine finteg1(tgs,fgs,n,x,y)
c=================================================
c     linear interpolation subroutine
c     tgs: a data array of X
c     fgs: a data array of Y
c     n: size of the array
c     x: input
c     y: output, interpolated value at x.

      implicit none
      
      real(8) :: tgs,fgs,x,y, a, b
      integer :: n, j
      dimension tgs(n),fgs(n)

      if((x.gt.tgs(n).and.x.gt.tgs(1)).
     &	    or.(x.lt.tgs(n).and.x.lt.tgs(1)))then
	  write(*,*)'In finteg1 overflow data'
	  stop
	endif
      call locate(tgs,n,x,j)
      if(tgs(j).eq.tgs(j+1)) write(*,*)'In finteg1 equal value, j = ',j
      a=(x-tgs(j+1))/(tgs(j)-tgs(j+1))
      b=(x-tgs(j))/(tgs(j+1)-tgs(j))

      y=a*fgs(j)+b*fgs(j+1)
      return
      end      
      
      
c===================================
      subroutine locate(xx,n,x,j)
c===================================
      implicit none
      real(8) :: xx, x
      integer :: n, j, jl, ju, jm
      
      dimension xx(n)

      jl=1
      ju=n

      do while ((ju-jl).gt.1)
	  jm=(ju+jl)/2
	  if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
	    jl=jm
	  else
	    ju=jm
	  endif
      enddo

      j=jl

      return
      end


c=========================================================
      subroutine interp(n1,n2,nx,nxn,ny,nyn,x,xn,y,yn,p)
c==========================================================

      use Q4_Sizes
      implicit none
      
      real(8) :: x, xn, y, yn, p
      integer :: n1,n2,nx,nxn,ny,nyn, ixn(nxx), jyn(nyx)
      real(8) :: cox(nxx), coy(nyx), pn(nxx,nyx)
      integer :: ix, j, jy, i

      dimension x(nx),xn(nxn),y(ny),yn(nyn), p(nx,ny) !  p(n1,n2)

!     find the neighboring grid index
      ixn(1)=1
      ixn(nxn)=nx-1
      do i=2,nxn-1
	  ix=ixn(i-1)
	  do while(.not.(x(ix).le.xn(i) .and. x(ix+1).gt.xn(i)))
	    ix=ix+1
	  enddo
	  ixn(i)=ix
!	interpolation coefficient
	  cox(i)=(x(ix+1)-xn(i))/(x(ix+1)-x(ix))
      enddo
      cox(1)=1.d0
      cox(nxn)=0.d0
      
      jyn(1)=1
      jyn(nyn)=ny-1
      do j=2,nyn-1
	  jy=jyn(j-1)
	  do while(.not.(y(jy).le.yn(j) .and. y(jy+1).gt.yn(j)))
	    jy=jy+1
	  enddo
	  jyn(j)=jy
	  coy(j)=(y(jy+1)-yn(j))/(y(jy+1)-y(jy))
      enddo
      coy(1)=1.d0
      coy(nyn)=0.d0

      do i=1,nxn
	  do j=1,nyn
	    pn(i,j)=cox(i)*(p(ixn(i),jyn(j))*coy(j)+
     &	  	    p(ixn(i),jyn(j)+1)*(1.d0-coy(j)))+(1.d0-cox(i))
     &		    *(p(ixn(i)+1,jyn(j))*coy(j)+
     &		    p(ixn(i)+1,jyn(j)+1)*(1.d0-coy(j)))
	  enddo
      enddo

      do i=1,nxn
	  do j=1,nyn
	    p(i,j)=pn(i,j)
	  enddo
      enddo

      return
      end
c=========================================================
      subroutine interp1(n1,nx,nxn,x,xn,p)
c==========================================================
c     BC -- function seems to stretch (interpolate) p out to dimensions of nxn
      implicit none
      real(8) :: x,xn,p
      integer :: n1,nx,nxn
      integer :: i, ix, ixn
      real(8) :: cox, pn
      dimension cox(1001),pn(1001),x(nx),xn(nxn),p(n1),ixn(1001)

c     find the neighbouring grid index
      ixn(1)=1
      ixn(nxn)=nx-1
      do i=2,nxn-1
	  ix=ixn(i-1)
	  do while(.not.(x(ix).le.xn(i) .and. x(ix+1).gt.xn(i)))
	    ix=ix+1
	  enddo
	  ixn(i)=ix
c	  interpolation coefficient
	  cox(i)=(x(ix+1)-xn(i))/(x(ix+1)-x(ix))
      enddo
      cox(1)=1.d0
      cox(nxn)=0.d0
      do i=1,nxn
	  pn(i)=cox(i)*p(ixn(i))+(1.d0-cox(i))*p(ixn(i)+1)
      enddo
      do i=1,nxn
	  p(i)=pn(i)
      enddo

      return
      end

c======================================
      subroutine sort(iorder,eps,x,ix)
c======================================
c iorder = 1: strictly ascending, 
c iorder = 2: strictly descending
c iorder = 3: ascending
c iorder = 4: descending

      implicit none
      
      integer :: iorder
      integer :: k, i, icount
      real(8) :: xtp
      real(8) :: eps,x
      integer :: ix
      
      dimension x(ix)

c inefficient bubble sort for small arrays
      if(mod(iorder,2).eq.1)then
         do k=1,ix-1
            do i=1,ix-k
               if(x(i).gt.x(i+1)) then
                  xtp=x(i)
                  x(i)=x(i+1)
                  x(i+1)=xtp
               endif
            enddo
         enddo
      else
         do k=1,ix-1
            do i=1,ix-k
               if(x(i).lt.x(i+1)) then
                  xtp=x(i)
                  x(i)=x(i+1)
                  x(i+1)=xtp
               endif
            enddo
         enddo
	endif

	if(iorder.ge.3) return

c     get rid of redundancy 
      icount=1

      do i=2,ix
	   if(dabs(x(i)-x(i-1)).gt.eps)then
	     icount=icount+1
	     x(icount)=x(i)
	   endif
      enddo

      ix=icount
      return
      end


c========================================
      subroutine abs_val_sort(iorder,eps,x,y,ix)
c========================================
      implicit none
      
      integer :: iorder
      real(8) :: eps,x
      integer :: k, i, ix
      real(8) :: y, xtp, ytp
      
      dimension x(ix), y(ix)

         do k=1,ix-1
            do i=1,ix-k
c --------------------------------------------
c this line modified to allow non-unique radii
c change by rdg 12/10/98, quick419e
c --------------------------------------------
               if(dabs(x(i)) .gt. dabs(x(i+1)).
     &           or. (dabs(x(i)) .eq. dabs(x(i+1)).
     &           and. y(i) .lt. y(i+1))) then
c ------------- end modification -------------
                  xtp=x(i)
                  x(i)=x(i+1)
                  x(i+1)=xtp
                  ytp=y(i)
                  y(i)=y(i+1)
                  y(i+1)=ytp
               endif
            enddo
         enddo
	return
	end



c========================================
      subroutine sort2(iorder,eps,x,y,ix)
c========================================
c sort ix-element-arrays x, y according to x
c iorder = 1: strictly ascending, 
c iorder = 2: strictly descending
c iorder = 3: ascending (sort by x then y)
c iorder = 4: descending (sort by x then y)
      implicit none
      integer :: iorder, k, i, icount, ix
      real(8) :: eps, x(ix), y(ix), xtp, ytp      

c inefficient bubble sort for small arrays
      if(mod(iorder,2).eq.1)then
         do k=1,ix-1
            do i=1,ix-k
c --------------------------------------------
c this line modified to allow non-unique radii
c change by rdg 12/10/98, quick419e
c --------------------------------------------
               if(x(i).gt.x(i+1).
     &           or.(x(i).eq.x(i+1).
     &           and.y(i).gt.y(i+1))) then
c ------------- end modification -------------
                  xtp=x(i)
                  x(i)=x(i+1)
                  x(i+1)=xtp
                  ytp=y(i)
                  y(i)=y(i+1)
                  y(i+1)=ytp
               endif
            enddo
         enddo
      else
         do k=1,ix-1
            do i=1,ix-k
c --------------------------------------------
c this line modified to allow non-unique radii
c change by rdg 12/10/98, quick419e
c --------------------------------------------
               if(x(i).lt.x(i+1).
     &           or.(x(i).eq.x(i+1).
     &           and.y(i).lt.y(i+1)))then
c -------------- end modification ------------
                  xtp=x(i)
                  x(i)=x(i+1)
                  x(i+1)=xtp
                  ytp=y(i)
                  y(i)=y(i+1)
                  y(i+1)=ytp
               endif
            enddo
         enddo
	endif

	if(iorder.ge.3) return

      icount=1

      do i=2,ix
	   if(dabs(x(i)-x(i-1)).gt.eps)then
	     icount=icount+1
	     x(icount)=x(i)
	     y(icount)=y(i)
	   endif
      enddo

      ix=icount
      return
      end
      
c==============================================
      subroutine matrix_multi(n1,n2,n3,a,b,c)
c=============================================

      implicit none
      integer :: n1, n2, n3, i, j, k
      real(8) :: a, b, c
      dimension a(n1,n2), b(n2,n3), c(n1,n3)

      do i=1,n1
	  do j=1,n3
	      c(i,j)=0.d0
	      do k=1,n2
	        c(i,j)=c(i,j)+ a(i,k)*b(k,j)
	      enddo
	  enddo
      enddo
      return
      end
      
c======================================
      subroutine matrix33_inverse(a)
c======================================
      implicit none
      real(8) :: a(3,3)
      real(8) :: deter, ai11, ai12, ai13, ai21, ai22, ai23, ai31, 
     &          ai32, ai33
      
      deter=a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) +
     &	    a(2,1)*a(3,2)*a(1,3) - a(1,3)*a(2,2)*a(3,1) -
     &	    a(1,2)*a(2,1)*a(3,3) - a(3,2)*a(2,3)*a(1,1)

      ai11=a(2,2)*a(3,3)-a(2,3)*a(3,2)
      ai12=-(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      ai13=a(2,1)*a(3,2)-a(3,1)*a(2,2)
      ai21=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
      ai22=a(1,1)*a(3,3)-a(1,3)*a(3,1)
      ai23=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))
      ai31=a(1,2)*a(2,3)-a(1,3)*a(2,2)
      ai32=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
      ai33=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      a(1,1)=ai11/deter
      a(1,2)=ai21/deter
      a(1,3)=ai31/deter
      a(2,1)=ai12/deter
      a(2,2)=ai22/deter
      a(2,3)=ai32/deter
      a(3,1)=ai13/deter
      a(3,2)=ai23/deter
      a(3,3)=ai33/deter

      return

      end
c==========================================
      subroutine equate(ndx,ndy, nx,ny, xl, xr)
c==========================================
      implicit none
      
      integer :: ndx,ndy,nx,ny
      !real*8 xl(ndx,ndy),xr(ndx,ndy)
      real*8 xl(nx,ny),xr(ndx,ndy)
      integer :: i, j
      
      do i=1,nx
	  do j=1,ny
	    xl(i,j)=xr(i,j)
	  end do
      end do

      return
      end
      
      
c==========================================
      subroutine equateToLarger(ndx,ndy, nx,ny, xl, xr)
c==========================================
      implicit none
      
      integer :: ndx,ndy,nx,ny
      !real*8 xl(ndx,ndy),xr(ndx,ndy)
      real*8 xl(nx,ny),xr(ndx,ndy)
      integer :: i, j
      
      do i=1,ndx
	  do j=1,ndy
	    xl(i,j)=xr(i,j)
	  end do
      end do

      return
      end
      
      
c===============================================
      subroutine tridag(n,a,b,c,d)
c===============================================

      use Q4_globals
      implicit none
      
      integer :: n, i, im1
      real(8) :: a,b,c,d
      real(8) :: beta, gamma
      dimension a(n),b(n),c(n),d(n)
      dimension beta(nmx),gamma(nmx)

c     this subroutine inverts a tri-diagonal matrix

      beta(1)=b(1)
      gamma(1)=d(1)/beta(1)

      do i=2,n
	im1=i-1
	beta(i)=b(i)-a(i)*c(im1)/beta(im1)
	gamma(i)=(d(i)-a(i)*gamma(im1))/beta(i)
      enddo

      d(n)=gamma(n)

      do i=n-1,1,-1
	d(i)=gamma(i)-c(i)*d(i+1)/beta(i)
      enddo

      return
      end

c===================================================================
      subroutine restrict(pfine,pcors,nxf,nyf,ndxf,ndyf,ndxc,ndyc)
c===================================================================
      implicit none
      real(8) :: pfine, pcors
      integer :: nxf, nyf, ndxf, ndyf, ndxc, ndyc
      integer :: nxc, nyc, i, j
      dimension pfine(ndxf,ndyf), pcors(ndxc,ndyc)

      nxc=(nxf+1)/2
      nyc=(nyf+1)/2

      do i=1,nxc
	  do j=1,nyc
	    pcors(i,j)=pfine(2*i-1,2*j-1)
	  enddo
      enddo


!
!     BC - Full weighting was found to have negligable impact on solution speed and convergence
!
!      do i = 1, nxc
!	  do j = 1, nyc
!	    if (i .eq. 1 .or. i .eq. nxc .or. 
!     &        j .eq. 1 .or. j .eq. nyc) then
!              pcors(i,j) = pfine(2*i-1, 2*j-1) 
!          else
!	      pcors(i,j) = 0.25*pfine(2*i-1, 2*j-1) + 
!     &      0.0625 * (pfine(2*i-2, 2*j-2)+pfine(2*i, 2*j-2)) + 
!     &      0.0625 * (pfine(2*i-2, 2*j)+pfine(2*i, 2*j)) + 
!     &      0.125 * (pfine(2*i-1, 2*j-2)+pfine(2*i-1, 2*j)) +  !top, bottom
!     &      0.125 * (pfine(2*i-2, 2*j-1)+pfine(2*i, 2*j-1))    !left, right
!          endif
!	  enddo
!      enddo

      return
      end

c==================================================================
      subroutine prolong(pcors,pfine,nxc,nyc,ndxf,ndyf,ndxc,ndyc)
c==================================================================
      implicit none

      real(8) :: pcors,pfine
      integer :: nxc,nyc,ndxf,ndyf,ndxc,ndyc
      integer :: nxf, nyf, i, j
      dimension pfine(ndxf,ndyf),pcors(ndxc,ndyc)

      !BC - Restriction and prolongation operators from Shyy 
      !and Sun (1992)
      !     Simple and effective...
      nxf=2*nxc-1
      nyf=2*nyc-1
      
      do i=1,nxc
	  do j=1,nyc
	
	    pfine(2*i-1,2*j-1)=pcors(i,j)
	
	    if(j.ne.nyc)then
	      pfine(2*i-1,2*j)=(pcors(i,j)+pcors(i,j+1))/2.d0
	
	      if(i.ne.nxc) then
	        pfine(2*i,2*j)=(pcors(i+1,j+1)+pcors(i,j+1)+
     &				         pcors(i+1,j)+pcors(i,j))/4.d0
            endif
	    endif
	  
	    if(i.ne.nxc)pfine(2*i,2*j-1)=(pcors(i,j)+pcors(i+1,j))/2.d0
	  
	  enddo
      enddo      


!     BC - testing out bilinear interpolation.  It doesn't change convergence much or speed things up much at all
!
!      do j = 1, nyf
!        do i = 1, nxf
!          pfine(i, j) = 0.0
!        enddo
!      enddo
!
!
!      do i=1, nxc
!	  do j=1, nyc
!	
!	    !center
!	    pfine(2*i-1,2*j-1) = pcors(i,j)
!	    
!	    !upperLeft
!	    if (i .ne. 1 .and. j .ne. 1) then
!	      pfine(2*i-2,2*j-2) = pfine(2*i-2,2*j-2) + 0.25*pcors(i,j)
!	    endif
!	
!	    !upperCenter
!	    if(j .ne. 1)then
!	      pfine(2*i-1,2*j-2) = pfine(2*i-1,2*j-2) + 0.50*pcors(i,j)
!	    endif
!	    
!	    !upperRight
!	    if (i .ne. nxc .and. j .ne. 1) then
!	      pfine(2*i, 2*j-2) = pfine(2*i, 2*j-2) + 0.25*pcors(i,j)
!	    endif
!	    
!	    !left
!	    if (i .ne. 1) then
!	      pfine(2*i-2, 2*j-1) = pfine(2*i-2, 2*j-1) + 0.50*pcors(i,j)
!	    endif
!	    
!	    !right
!	    if (i .ne. nxc) then
!	      pfine(2*i,2*j-1) = pfine(2*i, 2*j-1) + 0.50*pcors(i,j)
!	    endif
!	    
!	    !lowerLeft
!	    if (i .ne. 1 .and. j .ne. nyc) then
!	      pfine(2*i-2, 2*j) = pfine(2*i-2, 2*j) + 0.25*pcors(i,j)
!	    endif
!	    !lower Center
!	    if (j .ne. nyc) then
!	      pfine(2*i-1, 2*j) = pfine(2*i-1, 2*j) + 0.50*pcors(i,j)
!	    endif
!	    
!	    !lower Right
!	    if (i .ne. nxc .and. j .ne. nyc) then
!	      pfine(2*i, 2*j) = pfine(2*i, 2*j) + 0.25*pcors(i,j)
!	    endif
!	  
!	  enddo
!      enddo

      return
      end
      
c===================================================
      subroutine matrixsub(dif,p1,p2,nx,ny,ndx,ndy)
c===================================================
      implicit none
      real(8) :: dif,p1,p2
      integer :: nx,ny,ndx,ndy, i, j
      dimension dif(ndx,ndy),p1(ndx,ndy),p2(ndx,ndy)

      do i=1,nx
	  do j=1,ny
	    dif(i,j)=p1(i,j)-p2(i,j)
	  enddo
      enddo

      return
      end

c===================================================
      subroutine matrixinc(p,dp,nx,ny,ndx,ndy,ml,negp)
c===================================================
      implicit none
      
      real(8) :: p,dp, pnew
      integer :: nx,ny,ndx,ndy,ml,negp
      integer :: i, j
      dimension p(ndx,ndy),dp(ndx,ndy)

      ! boudary is not included here
      negp=0
      do i=2,nx-1
	  do j=2,ny-1
            ! no negative increment allowed here
	      pnew=p(i,j)+dp(i,j)
	      if(pnew.lt.0.d0)then
	          if(ml.eq.0) then
	              p(i,j)=0.05d0
	          else
	              p(i,j)=pnew
	              negp=1
	          endif
	      else
	          p(i,j)=pnew
	      endif
        enddo
      enddo
      
      return
      end
      
c==========================================
      subroutine wmatrix(p, nx, ny, fn)
c==========================================

      use Q4_sizes
      implicit none

      real(8) :: p(nxx, nyx)
      integer :: nx, ny
      integer :: i, j
      
      character*12 fn

      call openout(13, fn, 23 * nmx)

      do i = 1, nx
	    write(13, *) (p(i, j), j = 1, ny)
      enddo

      close(13)

      return
      end
      
c===================================================
      subroutine matInit(p, nxx, nyx, nx, ny, value)
c===================================================

      implicit none
      real(8) :: p, value
      integer :: nxx, nyx, nx, ny
      integer :: i, j

      dimension p(nxx, nyx)

      do i = 1, nx
         do j = 1, ny
            p(i, j) = value
         enddo
      enddo

      return
      end
      
c===================================================
      subroutine endTime()
c===================================================
      
      use Q4_globals
      implicit none
      
	
      CHARACTER (LEN = 12) REAL_CLOCK (3)
      CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2),
     &  REAL_CLOCK (3), end_time )
	
	write(*,*)
	write(*,10) 'SIMULATION STARTED AT: ', START_TIME(5), START_TIME(6),
     &  START_TIME(7)
10	format(1X, A, I2.2, ':', I2.2, ':', I2.2) 
	write(*,15) 'SIMULATION ENDED AT:   ', end_time(5), end_time(6),
     &  end_time(7)
15	format(1X, A, I2.2, ':', I2.2, ':', I2.2) 
	write(*,*)
	end