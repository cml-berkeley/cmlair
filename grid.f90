
!===========================================
    subroutine adaptive
!===========================================
    use NestingInfo
    use Q4_globals
    !use AdaptiveGridParams
    implicit none
    integer :: i
    
    !interp p onto coarser meshes
    call restrictPressure(p, p1, xref, yref, xref1, yref1, &
                          nx, ny, nx1, ny1, enclosingFineRectX1, enclosingFineRectY1)
    call restrictPressure(p1, p2, xref1, yref1, xref2, yref2, &
                          nx1, ny1, nx2, ny2, enclosingFineRectX2, enclosingFineRectY2)
    call restrictPressure(p2, p3, xref2, yref2, xref3, yref3, &
                          nx2, ny2, nx3, ny3, enclosingFineRectX3, enclosingFineRectY3)
    call restrictPressure(p3, p4, xref3, yref3, xref4, yref4, &
                          nx3, ny3, nx4, ny4, enclosingFineRectX4, enclosingFineRectY4)

    !call adaptNonNested on each level
    call adaptNonNested(p, himax, hnew, xref, yref, nx, ny, xl, yl, 1)
    call snapGridNonNested(xref, nx, icntlx, xcontrol, 1.d0)
    call snapGridNonNested(yref, ny, icntly, ycontrol, yl)
    
    call adaptNonNested(p1, himax1, hnew1, xref1, yref1, nx1, ny1, xl, yl, 1)
    call snapGridCoarse(xref1, nx1, icntlx, xcontrol, 1.d0)
    call snapGridCoarse(yref1, ny1, icntly, ycontrol, yl)
    
    call adaptNonNested(p2, himax2, hnew2, xref2, yref2, nx2, ny2, xl, yl, 1)
    call snapGridCoarse(xref2, nx2, icntlx, xcontrol, 1.d0)
    call snapGridCoarse(yref2, ny2, icntly, ycontrol, yl)
    
    call adaptNonNested(p3, himax3, hnew3, xref3, yref3, nx3, ny3, xl, yl, 1)
    call snapGridCoarse(xref3, nx3, icntlx, xcontrol, 1.d0)
    call snapGridCoarse(yref3, ny3, icntly, ycontrol, yl)
    
    call adaptNonNested(p4, himax4, hnew4, xref4, yref4, nx4, ny4, xl, yl, 1)
    call snapGridCoarse(xref4, nx4, icntlx, xcontrol, 1.d0)
    call snapGridCoarse(yref4, ny4, icntly, ycontrol, yl)
    
    hflag = ON
    uflag = ON

    call outputGrid
    
    return
    end subroutine adaptive
    

!==========================================================
    subroutine get_grids
!==========================================================
!   Set up the initial grid based on the Geometric series data
!   Mostly this isn't used these days so we just have a regularly
!   spaced grid over the slider

    use Q4_globals
    implicit none
    
    integer :: i, j, k
    real(8) :: const

    real(8) :: dxref(nxx), dyref(nyx)

    ! boundary grid information
    xnt(1)=0.d0
    nxt(1)=1
    xnt(nsx+1)=1.d0
    nxt(nsx+1)=nx

    ynt(1)=0.d0
    nyt(1)=1

    if(isymmetry.eq.1) then
        ! for symmetric y grid, calculate half first
	    ynt(nsy+1)=yl / 2
	    nyt(nsy+1)=int(ny/2)+1
    else
	    ynt(nsy+1)=yl
	    nyt(nsy+1)=ny
    endif


    ! start calculate x grid
    do i=1,nsx
	    xref(nxt(i))=xnt(i)
	    dxref(nxt(i))=1.d0

	    do k=nxt(i)+1,nxt(i+1)
	        dxref(k)=dxref(k-1)*dxr(i)
	        xref(k)=xref(k-1)+dxref(k-1)
	    enddo 

        ! calculate the constant for adjusting the grid step length
	    const=(xnt(i+1)-xnt(i))/(xref(nxt(i+1))-xnt(i))

	    do k=nxt(i)+1,nxt(i+1)
	        dxref(k-1)=dxref(k-1)*const
	        xref(k)=xref(k-1)+dxref(k-1)
	    enddo
    enddo

    ! start calculate y grid
    do j=1,nsy
	    yref(nyt(j))=ynt(j)
	    dyref(nyt(j))=1.d0

	    do k=nyt(j)+1,nyt(j+1)
	        dyref(k)=dyref(k-1)*dyr(j)
	        yref(k)=yref(k-1)+dyref(k-1)
	    enddo

	    if(isymmetry.eq.1.and.j.eq.nsy .and. mod(int(ny), 2).eq.0)then
        ! for symmetric y grid and even number of total y grid point,
        ! half grid at center
	        const=(ynt(j+1)-ynt(j))/(yref(nyt(j+1))-ynt(j)-dyref(nyt(nsy+1)-1)/2.d0)
	    else
	        const=(ynt(j+1)-ynt(j))/(yref(nyt(j+1))-ynt(j))
	    endif

	    do k=nyt(j)+1,nyt(j+1)
	        dyref(k-1)=dyref(k-1)*const
	        yref(k)=yref(k-1)+dyref(k-1)
	    enddo
    enddo

    ! for symmetric y grid, generate the other half
    if(isymmetry.eq.1)then
	    do k=ny,int(ny/2)+1,-1
	        yref(k) = yl - yref(ny-k+1)
    	enddo
    endif

    return
    end subroutine get_grids
    

!============================================================
    subroutine adpt_grid(nx,nxp,xl,x,dnp,xp,decay,xt,axt,dxt)
!============================================================
!nx - size of enlarged grid
!nxp - size of actual grid
!xl - length of grid (1.0 or yl)
!x - evenly spaced grid from 0 to length
!dnp - adaptation curve/pressure gradient weigthing function
!xp - new grid locations
!grid decay
!taper stuff
!???    
    
    implicit none
    
    integer :: nx,nxp, k, j
    real(8) :: a, b, c
    real(8) :: xl, decay, xt, axt, dxt, const
    
    real*8  n(1001),dn(1001),dx(1001),x(nx),dnp(nx),xp(nxp)

    !calculate space between existing grid points
    do k=1,nx-1
	    dx(k) = x(k+1) - x(k)
    enddo

    !incorporate decay into grid weighting function
    !weighting = weighting + sum (weighting at other points * exponential distance from other points)
    do k=1, nx !calculate the gradiant based on neighboring gradients
	    dn(k) = 0.d0
	    do j=1, nx
	        dn(k) = dn(k) + dnp(j) * dexp( -dabs(x(k)-x(j))*decay )
	    enddo
    enddo

    ! calculate scaling constant for dn(k)
    const = 0.d0
    do k = 1, nx-1
	    const = const + dx(k) * (dn(k)+dn(k+1)) / 2.d0
    enddo

    do k=1,nx !does nothing since axt is 0 in second itr of adaptive grid
	    dn(k) = dn(k) + const * axt * dexp(-dabs(x(k)-xt)*dxt)
    enddo

    const=0.d0 !calculate scaling constant for dn(k) again...not really necessary since axt is 0
    do k=1,nx-1
	    const = const + dx(k)*(dn(k)+dn(k+1))/2.d0
    enddo
    const=dfloat(nxp-1)/const

    ! scaling dn(k) with const
    do k=1,nx
	    dn(k) = dn(k) * const
    enddo

    ! calculate fractional grid number at old grid points
    !basically at each point, 1 through 1001, figure out how many grid points
    !were supposed to come before this point.  n(1001) will be equal to the number of grid points
    n(1) = 1.d0
    n(nx) = dfloat(nxp)
    do k = 2, nx
	    n(k) = n(k-1) + dx(k-1)*(dn(k-1)+dn(k))/2.d0
    enddo

    ! calculate new grid location (grid location is xp)
    do k=1,nx-1
	    a = (dn(k+1) - dn(k)) / 2.d0 / dx(k)
	    b = (x(k+1) * dn(k) - x(k) * dn(k+1)) / 2.d0 / dx(k)
	    c = n(k) - x(k) * ((2.d0*x(k+1) - x(k)) * dn(k) - x(k)*dn(k+1)) / 2.d0 / dx(k)
	    do j = int(n(k)+1), int(n(k+1))
	        if (dabs(a) .gt. dabs(b)) then
	            xp(j) = (dsqrt(b**2 - a*(c-dfloat(j))) - b) / a
	        else
	            xp(j) = (dfloat(j) - c) / (dsqrt(b**2 - a*(c-dfloat(j))) + b)
	        endif
	    enddo
    enddo

    xp(1)=0.d0
    xp(nxp)=xl

    return
    end subroutine adpt_grid
    
    
!===========================================
    subroutine adaptNonNested(p, himax, hnew, xref, yref, nx, ny, xl, yl, meshLevel)
!===========================================    
    use Q4_sizes
    use AdaptiveGridParams
    implicit none
    
    integer nx, ny, meshLevel
    real(8) :: xl, yl
    real(8) :: p(nx, ny), himax(nx, ny), hnew(nx, ny), xref(nx), yref(ny)
    
    !local data
    real(8) :: dnx(1001), dny(1001), xrefn(nxx), yrefn(nyx)
    real(8) :: xtemp(1001), ytemp(1001)
    real(8) :: ptemp(nxx,nyx)

	integer :: ngrad, nxn, nyn, i, j, im, jm1, ienough, iround, ind, ninc
	real(8) :: axt, dxt, dnxx, dnxmin, dnyx, dnymin, hip1j, him1j, hij, rcnt
	real(8) :: pmaxind, pmaxi    
    
    ngrad=1001
    
    nxn=nx
    nyn=ny

    ! amplitude for grid density at taper end xt
    axt=2.d0

    !  exponential decay rate for xt grid density
    dxt=10000.d0
    ! dxt=100.d0

    dnxx=0.001d0
    dnxmin=1.d20
    dnyx=0.001d0
    dnymin=1.d20

19  continue

!   calculate average gradient in x direction
    do i=1,nx
        ! Changed on 7/23/03 for geometric adaptive grid, p->h.
        ! avoid zero gradient
        dnx(i)=0.001d-9
        ! dnx(i)=0.001
         	
        do j=2,ny-1
            !geometry adaptation
            if (iadaptMode_x .eq. 0) then
	            if(i.gt.1.and.i.lt.nx) then
	                hip1j=himax(i+1,j)+hnew(i+1,j)
	                him1j=himax(i-1,j)+hnew(i-1,j)
	                hij=himax(i,j)+hnew(i,j)
	            else
	                if(i.eq.1) then
	                    hip1j=himax(i+1,j)+hnew(i+1,j)
	                    him1j=himax(1,j)+hnew(1,j)
	                    hij=himax(i,j)+hnew(i,j)
	                else
	                    hip1j=himax(i,j)+hnew(i,j)
	                    him1j=himax(i-1,j)+hnew(i-1,j)
	                    hij=himax(i,j)+hnew(i,j)
 	                endif
	            endif
	            if(ipmax_x.eq.1)then
	                if(i.eq.1) then
		                dnx(i)=dmax1(dnx(i),+dabs(hip1j-hij))
	                else
		                dnx(i)=dmax1(dnx(i),+dabs(hij-him1j))
	                endif
	            else
	                if(i.eq.1)then
		                dnx(i)=dnx(i)+dabs(hip1j-hij)
	                else
		                dnx(i)=dnx(i)+dabs(hij-him1j)
	                endif
	            endif !ipmax_x
	        else
	            if(ipmax_x.eq.1)then
	                if(i.eq.1) then
	   	                dnx(i)=dmax1(dnx(i),+dabs(p(i+1,j)-p(i,j)))
	                else
		                dnx(i)=dmax1(dnx(i),+dabs(p(i,j)-p(i-1,j)))
	                endif
	            else
	                if(i.eq.1)then
		                dnx(i)=dnx(i)+dabs(p(i+1,j)-p(i,j))
	                else
		                dnx(i)=dnx(i)+dabs(p(i,j)-p(i-1,j))
	                endif 
	            endif !ipmax_x .eq. 1
	        endif
	    enddo

	    if(i.eq.nx)then
	        dnx(i)=dnx(i)/(xref(i)-xref(i-1))
	    else if(i.eq.1) then
	        dnx(i)=dnx(i)/(xref(i+1)-xref(i))
	    else
	        dnx(i)=dnx(i)/(xref(i+1)-xref(i-1))
	    endif

        ! maximum and minimum gradient
	    if(dnxx.lt.dnx(i)) then
	        dnxx=dnx(i)
	    endif
	    if(dnxmin.gt.dnx(i)) then
	        dnxmin=dnx(i)
	    endif
    enddo


    ! calculate average pressure gradient in y direction

    do j=1,ny
        ! avoid zero gradient
	    dny(j)=0.001d0
	    
	    do i=2,nx-1
	        if(ipmax_y.eq.1)then
	            if(j.eq.1)then
		            dny(j)=dmax1(dny(j),dabs(p(i,j+1)-p(i,j)))
	            else
		            dny(j)=dmax1(dny(j),dabs(p(i,j)-p(i,j-1)))
	            endif
	        else
	            if(j.eq.1)then
	                dny(j)=dny(j)+dabs(p(i,j+1)-p(i,j))
	            else
	                dny(j)=dny(j)+dabs(p(i,j)-p(i,j-1))
	            endif
	        endif
	    enddo

	    if(j.eq.ny)then
	        dny(j)=dny(j)/(yref(j)-yref(j-1))
	    else if(j.eq.1) then
	        dny(j)=dny(j)/(yref(j+1)-yref(j))
	    else
	        dny(j)=dny(j)/(yref(j+1)-yref(j-1))
	    endif

        ! maximum and minimum gradient
	    if(dnyx.lt.dny(j))dnyx=dny(j)
	    if(dnymin.gt.dny(j))dnymin=dny(j)
    enddo

    ! limit the difference between maximum and minimum gradient
    ! updated for geometry adaptation
    if (iadaptMode_x .eq. 0) then
        do i=1,nx
	        dnx(i) = dmin1(dnx(i), dnxmin*(difmax_x / meshLevel))
        enddo
    else
        do i=1,nx
	        dnx(i) = dmax1(dnx(i), dnxx/(difmax_x / meshLevel))
        enddo       
    endif

    do j=1,ny
	    dny(j) = dmax1(dny(j), dnyx/(difmax_y / meshLevel))
    enddo

!	8/7/2003
!     -------------ADAPTIVE MESH REFINEMENT---------------------------
!	nxmr: number of regions to be refined in x direction.
!	nymr: number of regions to be refined in y direction.	 
!   rlevelx(i): mesh refinement level for ith region in x direction.
!   rlevely(i): mesh refinement level for ith region in y direction.
!	xginit(i): starting of refinement region in x direction.
!	xgend(i): end of refinement region in x direction.
!	yginit(i): starting of refinement region in y direction.
!	ygend(i): end of refinement region in y direction.

!	Normalize mesh refinement regions by xl.
    do i=1,nxmr
	    xginit(i)=xginit(i)/xl
        xgend(i)=xgend(i)/xl
    enddo
    do i=1,nymr
	    yginit(i)=yginit(i)/xl
        ygend(i)=ygend(i)/xl
    enddo
    ! mesh refinement in x direction
    do i=1,nxmr
	    do im=1,nx
	        rcnt=rlevelx(i)
	        if(xref(im).gt.xginit(i).and.xref(im).lt.xgend(i)) then
	            dnx(im)=dmax1(dnx(im), rcnt*dnxmin*(difmax_x / meshLevel))
	        endif
	    enddo
    enddo
    ! mesh refinement in y direction
    do j=1,nymr
	    do jm1=1,ny
	        rcnt=rlevely(j)
	        if(yref(jm1).gt.yginit(j).and.yref(jm1).lt.ygend(j)) then
	            dny(jm1)=dmax1(dny(jm1), rcnt*dnyx)
	        endif
	    enddo
    enddo
!----End of adaptive mesh refinement----------------------

    !xtemp is evenly spaced values from 0 to 1
    !ytemp is evenly spaced values from 0 to yl
    do i=1,ngrad
	    xtemp(i)=dfloat(i-1)*(1.d0/dfloat(ngrad-1))
	    ytemp(i)=xtemp(i)*yl
    enddo

    ! interpolate dnx and dny to xtemp and ytemp grids
    ! this stretches dnx out to the size of xtemp (see comment above about ngrad)
    call interp1(ngrad,nx,ngrad,xref,xtemp,dnx)
    call interp1(ngrad,ny,ngrad,yref,ytemp,dny)


    ! create adaptive grids in x and y directions
    call adpt_grid(ngrad,nyn,yl,ytemp,dny,yrefn,decay_y,0.d0,0.d0,0.d0)
	call adpt_grid(ngrad,nxn,1.d0,xtemp,dnx,xrefn,decay_x,0.d0,0.d0,dxt)
	
	!interp back to nx and ny
	call interp(nxx,nyx,nx,nxn,ny,nyn,xref,xrefn,yref,yrefn,p)
!    endif

    ! use the old names for the grid coordinates
    do i=1,nxn
	    xref(i)=xrefn(i)
    enddo

    do j=1,nyn
	    yref(j)=yrefn(j)
    enddo

    nx=nxn
    ny=nyn

    !call snapGridNonNested
    
    
    return
    end subroutine adaptNonNested
    
    
!============================================
	subroutine InitialGridGen
!============================================
!	set initial grid

    use Q4_globals
    implicit none
    
    logical, parameter :: importAllLevels = .false.
    integer :: i

    if(ioldgrid .eq. 0)then
        ! create the regular finite difference mesh
        call get_grids
        !if(iadpt .eq. 1) call snapGridInitial
    else
	    open(31,err=999, file='x.dat', status='old')
	    read(31, *, err = 998) (xref(i), i=1,nx)
	    close(31)

	    if (importAllLevels .eq. .true.) then
	        open(31,err=999, file='x2.dat', status='old')
	        read(31, *, err = 998) (xref1(i), i=1,nx1)
	        close(31)
    	    
	        open(31,err=999, file='x3.dat', status='old')
	        read(31, *, err = 998) (xref2(i), i=1,nx2)
	        close(31)
    	    
	        open(31,err=999, file='x4.dat', status='old')
	        read(31, *, err = 998) (xref3(i), i=1,nx3)
	        close(31)
    	    
	        open(31,err=999, file='x5.dat', status='old')
	        read(31, *, err = 998) (xref4(i), i=1,nx4)
	        close(31)
        endif	    
	    
	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	    	    	    

	    open(31,err=999, file='y.dat', status='old')
	    read(31, *, err = 999) (yref(i), i=1,ny)
	    close(31)
	    
	    if (importAllLevels .eq. .true.) then
	        open(31,err=999, file='y2.dat', status='old')
	        read(31, *, err = 999) (yref1(i), i=1,ny1)
	        close(31)
    	    
	        open(31,err=999, file='y3.dat', status='old')
	        read(31, *, err = 999) (yref2(i), i=1,ny2)
	        close(31)
    	    
	        open(31,err=999, file='y4.dat', status='old')
	        read(31, *, err = 999) (yref3(i), i=1,ny3)
	        close(31)
    	    
	        open(31,err=999, file='y5.dat', status='old')
	        read(31, *, err = 999) (yref4(i), i=1,ny4)
	        close(31)	
        endif
        
        call CheckInitialPressureFile   	    	    
    endif

	return

998	write(*,*)'Error reading x.dat'
	write(numRes,*) 'Error reading x.dat'
	write(numRes,*) 'errorcode: xgriderr'
	stop

999	write(*,*)'Error reading y.dat'
	write(numRes,*) 'Error reading y.dat'
	write(numRes,*) 'errorcode: ygriderr'
	stop

	end subroutine InitialGridGen
	
	
!=======================================================================	
	subroutine SnapGridNonNested(grid, n, icntl, control, sliderSize)
!=======================================================================	
    use Q4_sizes
    implicit none
    
    integer :: n, icntl
    real(8) :: grid(n), control(nxx), sliderSize
    
    real(8) :: co
    integer :: index, indexold, i, ip, ik
    real*8 temp(nxx)

    do i=1,n
	    temp(i)=grid(i)
    enddo
    
    if (icntl .le. n) then
	    index=2
	    do i = 1, icntl-2
	        indexold=index
	        if (grid(index) .gt. control(i+1)) then
	            grid(index) = control(i+1)
	            index=index+1
	        else
	            do ip = index, n-1
		            if (control(i+1).ge.grid(ip) .and. control(i+1).lt.grid(ip+1)) then !control(i+1) is straddled by xref(ip) and xref(ip+1)
		                !find the closest grid point to the control point
		                if ( dabs(control(i+1)-grid(ip)) .le. dabs(control(i+1)-grid(ip+1)) ) then
			                index = ip
		                else
			                index = ip+1
		                endif
		                index = min( int(n+1-icntl+i), index)
		                grid(index)=control(i+1)
		                index=index+1
		                goto 10
		            endif
	            enddo
10	            continue
	        endif
	        co=(grid(index-1)-grid(indexold-1)) / (temp(index-1)-temp(indexold-1))
	        do ik = indexold, index-2  !I think this loop serves to regularly space the grid from the previous control point to the current control point
	            grid(ik)=grid(indexold-1)+co*(grid(ik) - temp(indexold-1))
	        enddo
	    enddo
	    co=(sliderSize-grid(index-1))/(sliderSize-temp(index-1))
	    do ik=index,n-1
	        grid(ik)=grid(index-1)+co*(grid(ik)-temp(index-1))
	    enddo
    endif
    
    return
    end subroutine snapGridNonNested	
	
!	
!!=======================================================================	
!	subroutine snapGridFinestCoarse(grid, n, icntl, control, sliderSize)
!!=======================================================================	
!    use Q4_sizes
!    implicit none
!    
!    integer :: n, icntl
!    real(8) :: grid(n), control(nxx)
!    real(8) :: co, sliderSize
!    integer :: index, indexold, i, ip, ik
!    real(8) :: temp(nxx)
!    real(8) :: controlSelected(nxx)
!    real(8) :: meshSpacing, curControlLocation
!    integer :: numControlPtsX
!    
!    !for each control point
!    index = 2
!    controlSelected(1) = control(1) !take the first point
!    numControlPtsX = 1
!    do i = 2, icntl-1
!        curControlLocation = control(i)
!        do while (grid(index) .lt. curControlLocation)
!            index = index+1
!        enddo
!        !at this point, grid(index) holds the location of the grid point one past the control point at control(i)
!        meshSpacing = grid(index) - grid(index-1) !this is why we start with index = 2...
!        !control points can't be TOO close together, keep control points 
!        if ((curControlLocation - controlSelected(numControlPtsX)) .gt. (0.4 * meshSpacing)) then
!            controlSelected(numControlPtsX + 1) = control(i)
!            numControlPtsX = numControlPtsX + 1
!        endif
!    enddo
!    !take last point
!    controlSelected(numControlPtsX + 1) = control(icntl)
!    numControlPtsX = numControlPtsX + 1
!    
!    do i=1,n
!	    temp(i)=grid(i)
!    enddo
!
!    if (numControlPtsX .le. n) then
!	    index=2
!	    do i = 1, numControlPtsX-2
!	        indexold=index
!	        if (grid(index) .gt. controlSelected(i+1)) then
!	            grid(index) = controlSelected(i+1)
!	            index=index+1
!	        else
!	            do ip = index, n-1
!		            if (controlSelected(i+1).ge.grid(ip) .and. controlSelected(i+1).lt.grid(ip+1)) then !controlSelected(i+1) is straddled by grid(ip) and grid(ip+1)
!		                !find the closest grid point to the control point
!		                if ( dabs(controlSelected(i+1)-grid(ip)) .le. dabs(controlSelected(i+1)-grid(ip+1)) ) then
!			                index = ip
!		                else
!			                index = ip+1
!		                endif
!		                index = min( int(n+1-numControlPtsX+i), index)
!		                grid(index) = controlSelected(i+1)
!		                index=index+1
!		                goto 10
!		            endif
!	            enddo
!10	            continue
!	        endif
!	        co=(grid(index-1)-grid(indexold-1)) / (temp(index-1)-temp(indexold-1))
!	        do ik = indexold, index-2  !I think this loop serves to regularly space the grid from the previous control point to the current control point
!	            grid(ik)=grid(indexold-1)+co*(grid(ik) - temp(indexold-1))
!	        enddo
!	    enddo
!	    co=(sliderSize-grid(index-1))/(sliderSize-temp(index-1))
!	    do ik=index,n-1
!	        grid(ik)=grid(index-1)+co*(grid(ik)-temp(index-1))
!	    enddo
!    endif
!    
!    return
!    end subroutine snapGridFinestCoarse	
!	
	
	
!=======================================================================	
	subroutine snapGridCoarse(grid, n, icntl, control, sliderSize)
!=======================================================================	
    use Q4_sizes
    implicit none
    
    integer :: n, icntl
    real(8), parameter :: SNAP_TOL = 0.6
    real(8) :: grid(n), control(nxx)
    real(8) :: co, sliderSize
    integer :: index, indexold, i, ip, ik
    real(8) :: temp(nxx)
    real(8) :: controlSelected(nxx)
    real(8) :: meshSpacing, curControlLocation
    integer :: numControlPtsX
    
    !for each control point
    index = 2
    controlSelected(1) = control(1) !take the first point
    numControlPtsX = 1
    do i = 2, icntl-1
        curControlLocation = control(i)
        do while (grid(index) .lt. curControlLocation)
            index = index+1
        enddo
        !at this point, grid(index) holds the location of the grid point one past the control point at control(i)
        meshSpacing = grid(index) - grid(index-1) !this is why we start with index = 2...
        !control points can't be TOO close together, keep control points at least SNAP_TOL apart
        if ((curControlLocation - controlSelected(numControlPtsX)) .gt. (SNAP_TOL * meshSpacing)) then
            controlSelected(numControlPtsX + 1) = control(i)
            numControlPtsX = numControlPtsX + 1
        endif
    enddo
    !take last point
    controlSelected(numControlPtsX + 1) = control(icntl)
    numControlPtsX = numControlPtsX + 1
    
    do i=1,n
	    temp(i)=grid(i)
    enddo

    if (numControlPtsX .le. n) then
	    index=2
	    do i = 1, numControlPtsX-2
	        indexold=index
	        if (grid(index) .gt. controlSelected(i+1)) then
	            grid(index) = controlSelected(i+1)
	            index=index+1
	        else
	            do ip = index, n-1
		            if (controlSelected(i+1).ge.grid(ip) .and. controlSelected(i+1).lt.grid(ip+1)) then !controlSelected(i+1) is straddled by grid(ip) and grid(ip+1)
		                !find the closest grid point to the control point
		                if ( dabs(controlSelected(i+1)-grid(ip)) .le. dabs(controlSelected(i+1)-grid(ip+1)) ) then
			                index = ip
		                else
			                index = ip+1
		                endif
		                index = min( int(n+1-numControlPtsX+i), index)
		                grid(index) = controlSelected(i+1)
		                index=index+1
		                goto 10
		            endif
	            enddo
10	            continue
	        endif
	        co=(grid(index-1)-grid(indexold-1)) / (temp(index-1)-temp(indexold-1))
	        do ik = indexold, index-2  !I think this loop serves to regularly space the grid from the previous control point to the current control point
	            grid(ik)=grid(indexold-1)+co*(grid(ik) - temp(indexold-1))
	        enddo
	    enddo
	    co=(sliderSize-grid(index-1))/(sliderSize-temp(index-1))
	    do ik=index,n-1
	        grid(ik)=grid(index-1)+co*(grid(ik)-temp(index-1))
	    enddo
    endif
    
    return
    end subroutine snapGridCoarse	

!=======================================================================
    subroutine snapGridInitial
!=======================================================================
!   We snap the grid to the control points for the initial (regularly spaced grid)
!   and for the adapted grid

!notice that we don't snap if icntlx is greater than nx...

    use Q4_globals
    implicit none
    
    real(8) :: co
    integer :: index, indexold, i, ip, ik
    real*8 xtemp(nxx),ytemp(nyx)

    do i=1,nx
	    xtemp(i)=xref(i)
    enddo

    do i=1,ny
	    ytemp(i)=yref(i)
    enddo


    if (icntlx .le. nx) then
	    index=2
	    do i = 1, icntlx-2
	        indexold=index
	        if (xref(index) .gt. xControl(i+1)) then
	            xref(index) = xControl(i+1)
	            index=index+1
	        else
	            do ip = index, nx-1
		            if (xControl(i+1).ge.xref(ip) .and. xControl(i+1).lt.xref(ip+1)) then !xControl(i+1) is straddled by xref(ip) and xref(ip+1)
		                !find the closest grid point to the control point
		                if ( dabs(xControl(i+1)-xref(ip)) .le. dabs(xControl(i+1)-xref(ip+1)) ) then
			                index = ip
		                else
			                index = ip+1
		                endif
		                index = min( int(nx+1-icntlx+i), index)
		                xref(index)=xControl(i+1)
		                index=index+1
		                goto 10
		            endif
	            enddo
10	            continue
	        endif
	        co=(xref(index-1)-xref(indexold-1)) / (xtemp(index-1)-xtemp(indexold-1))
	        do ik = indexold, index-2  !I think this loop serves to shift the grid and make it somewhat more uniform from the previous control point to the current control point
	            xref(ik) = xref(indexold-1) + co*(xref(ik) - xtemp(indexold-1))
	        enddo
	    enddo
	    co=(1.d0-xref(index-1))/(1.d0-xtemp(index-1))  !regularly space the final section of the grid
	    do ik = index, nx-1
	        xref(ik)=xref(index-1)+co*(xref(ik)-xtemp(index-1))
	    enddo
    endif

    if(icntly.le.ny)then
	    index=2
	    do i=1,icntly-2
	        indexold=index
	        if(yref(index).gt.ycontrol(i+1))then
	            yref(index)=ycontrol(i+1)
	            index=index+1
	        else
	            do ip=index,ny-1
		            if(ycontrol(i+1).ge.yref(ip) .and. ycontrol(i+1).lt.yref(ip+1))then
		                if(dabs(ycontrol(i+1)-yref(ip)) .le. dabs(ycontrol(i+1)-yref(ip+1)))then
			                index=ip
		                else
			                index=ip+1
		                endif
		                index=min( int(ny+1-icntly+i),index)
		                yref(index)=ycontrol(i+1)
		                index=index+1
		                goto 20
		            endif
	            enddo
20	            continue
	        endif
	        co=(yref(index-1)-yref(indexold-1)) / (ytemp(index-1)-ytemp(indexold-1))
	        do ik=indexold,index-2
	        yref(ik)=yref(indexold-1) + co*(yref(ik)-ytemp(indexold-1))
	        enddo
	    enddo
	    co=(yl-yref(index-1))/(yl-ytemp(index-1))
	    do ik=index,ny-1
	        yref(ik)=yref(index-1)+co*(yref(ik)-ytemp(index-1))
	    enddo
    endif
    
    call outputGrid 
    
    return
    end subroutine snapGridInitial
    
   
!=========================================
    subroutine getControl
!=========================================
!  set x and y grid control points

    use Q4_globals
    implicit none
    
    integer :: k, i

    icntlx = 1
    icntly = 1

    xcontrol(1) = 0.d0
    ycontrol(1) = 0.d0

    !ok, here we need to make sure that we don't go over 1000 points on the slider
    do k=1,numRails
	    do i=1, npoints(k)
	        if(xrail(i,k).gt.0.d0.and.xrail(i,k).lt.1.d0)then
	            icntlx = icntlx+1
	            xcontrol(icntlx) = xrail(i,k)
	        endif
	        if(yrail(i,k).gt.0.d0.and.yrail(i,k).lt.yl)then
	            icntly = icntly+1
	            ycontrol(icntly) = yrail(i,k)
	        endif
	    enddo
!	    icntlx = icntlx+1
!        xcontrol(icntlx) = 0.9146
!	    icntly = icntly+1
!	    ycontrol(icntly) = .3263
	
!c	  if(nwpoint(k).ge.2)then
!c	     if(wpoint(k,nwpoint(k)).gt.0.d0)then
!c	      do i=1, npoints(k)
!c		  if(xwallo(i,k).gt.0.d0.and.
!c     &		     xwallo(i,k).lt.1.d0)then
!c		    icntlx=icntlx+1
!c		    xcontrol(icntlx)=xwallo(i,k)
!c		  endif
!c		  if(ywallo(i,k).gt.0.d0.and.
!c     &		     ywallo(i,k).lt.yl)then
!c		    icntly=icntly+1
!c		    ycontrol(icntly)=ywallo(i,k)
!c		  endif
!c	       enddo
!c	     endif
!c	     if(wpoint(k,1).lt.0.d0)then
!c	      do i=1, npoints(k)
!c		  if(xwalli(i,k).gt.0.d0.and.
!c     &		     xwalli(i,k).lt.1.d0)then
!c		    icntlx=icntlx+1
!c		    xcontrol(icntlx)=xwalli(i,k)
!c		  endif
!c		  if(ywalli(i,k).gt.0.d0.and.
!c     &		     ywalli(i,k).lt.yl)then
!c		    icntly=icntly+1
!c		    ycontrol(icntly)=ywalli(i,k)
!c		  endif
!c	       enddo
!c	     endif
!c	   endif

    enddo

    if(xt .gt. 0.d0)then
	    icntlx = icntlx+1
	    xcontrol(icntlx) = xt
    endif

    icntlx = icntlx + 1
    xcontrol(icntlx) = 1.d0
    icntly = icntly + 1
    ycontrol(icntly) = yl

    if(icntlx .ge. 4) call sort(1, 1.d-3, xcontrol, icntlx)
    if(icntly .ge. 4) call sort(1, 1.d-3, ycontrol, icntly)

	return
	end subroutine getControl
 
!=======================================
	subroutine xtGrid(xtold)
!=======================================
! adjust x grids according to new xt (taper length)
    
    use Q4_globals
    implicit none
    
    real(8) :: xtold
    integer :: i, j, j1, j2, iold, inew
	real(8) :: xctemp(nxx), co

    ! create new control points with new xt point

    do i = 2, icntlx-1
	    if(dabs(xcontrol(i)-xtold).lt.1d-5)then
	        iold = i
	        goto 10
	    endif
    enddo
10  continue

    do i = 1, icntlx
	    if(xcontrol(i).gt.xt)then
	        inew = i
	        goto 20
	    endif
    enddo
20  continue


    j = 1
    do i = 1, icntlx
	    if(i.eq.inew)then
	        xctemp(j) = xt
	        j = j + 1
	    endif
	    if(i.ne.iold) then
	        xctemp(j) = xcontrol(i)
	        j = j + 1
	    endif
    enddo

    !adjust grid in segments
    j1 = 2
    do i = 1, icntlx-1
	    j2 = nx-1
	    do j = j1,  nx
	        if(xref(j).ge.xcontrol(i+1))then
	            j2 = j - 1
	            goto 30
	        endif
	    enddo
30	    continue
        ! adjust if either end has changed
	    if(xctemp(i).ne.xcontrol(i) .or. xctemp(i+1).ne.xcontrol(i+1)) then
	        co = (xctemp(i+1) - xctemp(i))/(xcontrol(i+1)-xcontrol(i))
	    do	j = j1, j2
	        xref(j) = xctemp(i)+(xref(j)-xcontrol(i))*co
	    enddo
	endif
	j1 = j2 + 1
    enddo


    do i = 1, icntlx
	    xcontrol(i) = xctemp(i)
    enddo

    return
    end subroutine xtGrid
	