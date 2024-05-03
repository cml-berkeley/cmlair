!========================================
      subroutine initws
!========================================
!
    use Q4_globals
    implicit none

!c     input parameters from data file : run.dat
!c
!c     in summary:
!c
!c	variables	usage
!c	---------	-----
!c	nx,ny		# of points in x & y axis
!c	isolv		0=normal solution,1=inverse solution, 2=geometry only
!c	hm		initial flying height of slider at TEC(m)
!c	h0	    pitch of slider
!c	hs		roll angle (rad) of slider
!c	xl		length of slider (m)
!c	xg		norm. position of cg from front of slider
!c	yl	       width of slider
!
!c	ityact		 0=no actuator,1=rotary,2=linear
!c	dact		angular or linear position of the actuator arm
!c	vact		velocity of the actuator arm
!c	xact		rotary actuator dimension in x-direction
!c
!c	yact		rotary actuator dimension in y-direction
!c	ske		skew angle (deg)
!c	ra		radial location of slider (m)
!c	rpm,	      surf. vel. (rpm)
!c	akmax		convergence parameter max for reynolds eq solver
!c			 (suggested value 0.1e-04)
!c	iqpo		 poiseuille flow,0 cont,1 lin,2 quad,3 boltzman
!c			asymptotic,4 boltzman exact
!c	pos		circumferential position of slider
!c	flh		disk runout
!c	f0		load acting on slider (used when isolv=1)
!c	xf0		position of f0 from the cg (x-dir.)
!c	yf0		position of f0 from the cg (y-dir.)
!c	emax		convergence parameter for inverse solution
!c
!c -------------------------------------------------------------
!c Add code for expiration date (added 11/15/99, rdg)
!c -------------------------------------------------------------
    !INTEGER DATE_TIME

    integer :: i
    
    CHARACTER (LEN = 12) REAL_CLOCK (3)
    CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), REAL_CLOCK (3), START_TIME)
    if (START_TIME(1).ge.2018) then
	    if (START_TIME(2).ge.9)then
	        write(*,*) 'This solver has expired.'
	        write(*,*) 'Please download the latest copy of CMLAir'
	        stop
	    endif
    endif
! -------------------------------------------------------------
! End Code for expiration
! -------------------------------------------------------------

!   iterations in multigrid cycles
    mitr(1) = 4
    mitr(2) = 6
    mitr(3) = 8
    mitr(4) = 8
    mitp(1) = 2
    mitp(2) = 4
    mitp(3) = 6
    mitp(4) = 8
    
    hflag = ON  !ON if the recess height of our slider has changed (crown camber etc.)
    uflag = ON  !ON if we need to update the bearing number
    tflag = ON  !ON if this is the first time through the traversal subroutine

    !iFirstRun = 1 !see doCase for a description of this param
    iRunCounter = 0
    write(*,*) '*******************************************'
    write(*,*) 'CML VERSION ',version,' Quick Solver'
    write(*,*) '*******************************************'
	write(*,*)
	write(*,10) 'Simulation started at:  ', START_TIME(5), START_TIME(6), &
	             START_TIME(7), START_TIME(2), START_TIME(3), START_TIME(1)
10	format(1X,A,I2.2,':',I2.2,':',I2.2,' -- ',I2.2,'/',I2.2,'/',I4.4) 
	write(*,*)

	numRes = 12
	open(numRes,err=999,file='result.dat',status='unknown')

    call readRail
    call checkRail
    call READRUN
    call AllocateRun
    call LoadInitialPressure
    call normalize


    if(ialt.ge.1) then
	    call setAlt(alts(1))
    else
	    call setAlt(0.0)
    endif

    ! Write out header in result.dat
    call initResultDat
 
    write(*,'(a22,i3,a3,i3)')' FINEST GRID SIZE   = ',nx,' x ',ny

    ! coarse grid size for multigrid
    nxm(1) = nx1
    nym(1) = ny1
    nxm(2) = nx2
    nym(2) = ny2
    nxm(3) = nx3
    nym(3) = ny3
    nxm(4) = nx4
    nym(4) = ny4
    
    nxm1 = nx-1
    nxm2 = nx-2
    nym1 = ny-1
    nym2 = ny-2

!c   define the system constants
!c
!c   in summary:
!c
!c	variable       usage
!c	--------       -----
!c	vis1		dynamic viscosity (kg/m/s)
!c	p0	        atmospheric pressure (pa)
!c	al	        mean free path of air (m)
!c	u0	        linear velocity at center of slider
!c	ske	        skew angle
!c	d0	        modified inverse knudsen number
!c	icoe	    index for boltzman flow
!c	nter	    number of terms for boltzman flow

    gama = 0.57721566490153286d0
    pir  = 1.77245385090551603d0  !sqrt(pi)
    pit  = 1.023326708d0
    nter = 11
    icoe = 0

    return
999 write(*,*)'TROUBLE OPENING FILES.'
    stop
    end


!=====================================
    subroutine READRUN
!=====================================

    use Q4_globals
    use AdaptiveGridParams
    implicit none
    
    integer :: i
    real(8) :: pitchinitial
    
	character*70 line,dversion
	
!     No longer need to use dversion to identify our data
!     use string manipulation routines to determine
!     number of parameters we expect to read in
    character *250 header !holds parameter header
    character *250 words(10) !holds header split into individual words
    integer nwords !number of words in the parameter header
    
    open(2, err=999,file='run.dat',status='old')

!	grab the version from the first line
	read(2,'(a70)')line
	dversion=line(13:17)
	
	!init optional parameters
	call InitOptionalParams

!      read(2, *)
    read(2, *)
!     ***************Solution Control***************
    read(2, *)
    read(2, *)
    read(2, *) istiff, isolv, ioldgrid, iadpt, isave
    if (iadpt .eq. 2) then
        iUseNewGridEachRun = 1
        iadpt = 1
    endif
    
!     ***************Intial Attitude***************
    read(2, *)
    read(2, *)
    read(2, *) hm, h0, hs
    pitchinitial = h0 !needed for computing syspension stiffnesses

!----------------------------------------------
!          Added in quick418
! ----------------------------------------------
    irailgeom = OFF
    gapExtrapolationRun = .false.

    !isolv:  0 = Fixed Attitude; 1 = Fly Height; 2 = Geometry Only
    if (isolv .eq. 2) then
        irailgeom = ON
        iadpt = OFF
        isave = OFF
	    istiff = OFF
    endif
!    Added by Ling Huang For Gap Extrapolation.
    if (isolv.eq.3) then
        gapExtrapolationRun = .true.
        irailgeom = OFF
	    iadpt = OFF
	    isave = OFF
	    istiff = OFF
! 23 Oct 2001 read the pint heights for gap extrapolation by Richard Blanco
        open(33,err=999, file='gap.dat', status='old')
	    read(33,*)(hintge(i), i = 1,4)
! 23 Oct 2001  end code modification by Richard Blanco
    endif


!     Modified by Xinjiang Shen and Lion Huang. Less fly height can be
!     less than 5 nm.
    if(hm.lt.5.d-9) then
	    hmin = hm / 5.d-9
	    hm = 5.d-9
    else
	    hmin = 1.d0  !normalized
    endif

    hx0 = - h0 * xl / hm  !(h0*xl) is the arc length. Then we normalize to hm
    h0  =  hmin - hx0
    hy  =   hs * xl / hm

!     ***************Runs***************
    read(2, *)
    !user defined sort of RPM and Radii added 6/5/05
    read(2,'(a250)') header
    call parse(header, words, nwords)
    iRadiiSort = 0
    if (nwords .lt. 4) then
        read(2, *) irad, irpm, ialt
    else
        read(2, *) irad, irpm, ialt, iRadiiSort
    endif
       
!     RDG --------  Error Check: Too many radii, rpm, altitudes? ---------
    if (irad .gt. maxrad) then
        write(*,*) 'Number of radii exceeds maximum of ', maxrad
	    write(numRes,*) 'Number of radii exceeds maximum of ', maxrad
	    write(numRes,*) 'errorcode: maxrad ', maxrad, irad
        stop
    endif
    if (irpm .gt. maxrpm) then
        write(*,*) 'Number of rpms exceeds maximum of ', maxrpm
	    write(numRes,*) 'Number of rpms exceeds maximum of ', maxrpm
	    write(numRes,*) 'errorcode: maxrpm ', maxrpm, irpm
        stop
    endif
    if (ialt .gt. maxalt) then
        write(*,*) 'Number of altitudes exceeds maximum of ', maxalt
	    write(numRes,*) 'Number of altitudes exceeds maximum of ', maxalt
	    write(numRes,*) 'errorcode: maxalt ', maxalt, ialt
        stop
    endif
!     ----------------------------------------------------

    if(irad .lt. 1) then
	    write(*,*) 'At least one radial position must be used!'
	    write(numRes,*) 'At least one radial position must be used!'
	    write(numRes,*) 'errorcode: ialtlt1'
	    stop
    endif

    if(irpm .lt. 1) then
	    write(*,*) 'At least one RPM must be used!'
	    write(numRes,*) 'At least one RPM must be used!'
	    write(numRes,*) 'errorcode: irpmlt1'
	    stop
    endif

    read(2, *)
    read(2, *) (radii(i), i = 1, irad)
    read(2, *)
    read(2, *) (skews(i), i = 1, irad)  

    ! Modified to allow sort of radii by absolute value of skew
    if (irad.ge.2 .and. iRadiiSort.eq.1) then
        call sort2(4, 1.d-4, radii, skews, irad)
    else if(irad.ge.2) then !iRadiiSort.eq.0
        call abs_val_sort(4, 1.d-4, skews, radii, irad)
    endif
	
    call setRadius(radii(1), skews(1))

    read(2, *)
    read(2, *) (rpms(i), i = 1, irpm)

    if(irpm.ge.2) call sort(2, 1.d-4, rpms, irpm)

    call setRPM(rpms(1))

    read(2, *)

!     ***************Air Parameters***************
!     Go with the old way of doing the alts
    if (ialt.ge.1) then
	    read(2, *) (alts(i), i = 1, ialt)
!	  sort altitude in ascending order
	    if(ialt.ge.2) call sort(1, 1.d0, alts, ialt)
	    read(2, *)
	    read(2, *)
	    read(2, *) p0, al, vis1
	    !read(2, *)
    else
	    read(2, *)
	    read(2, *)
	    read(2, *)
!           p0, al and vis1 are currently overwritten in setAlt for higher altitudes
	    read(2, *) p0, al, vis1
    endif
      
!     ***************Load Parameters***************
!     load parameters modified 10/21/04 in v4.23 to accomodate stiffnesses
    read(2, *)
    read(2,'(a250)') header
    call parse(header, words, nwords)
    if (nwords .lt. 4) then !no suspension stiffness parameters
        iUseStiffnesses = 0
        read(2, *) f0, xf0, yf0
        read(2, *)
        read(2, *) xfs, yfs, emax
    else
        read(2, *) f0, xf0, yf0, PSA, RSA, iUseStiffnesses
        read(2, *)
        read(2, *) xfs, yfs, emax, Pitch_Stiffness, Roll_Stiffness
    endif
      
!     easiest method is to update xfs/yfs with correct torque every iteration
    if (iUseStiffnesses .eq. 1) then
!       initial update for start of simulation
        xfs = (PSA - pitchinitial) * Pitch_Stiffness * 1000000.0
        yfs = (RSA - hs) * Roll_Stiffness * 1000000.0
    endif

!     ***************Grid Control***************
    read(2, *)
    read(2, *)
    read(2, *) nx, ny

!     RDG --------  Error Check: Grid number too high? ---------
    nx= max0(  min0( int((nx+7)/16)*16+1, int((nxx-2)/16)*16+1)  , 33  )
    ny= max0(  min0( int((ny+7)/16)*16+1, int((nyx-2)/16)*16+1)  , 33  )

    if (nx .gt. nxx) then
        write(*,*) 'X grid size exceeds maximum of ', nxx
	    write(numRes,*) 'X grid size exceeds maximum of ', nxx
	    write(numRes,*) 'errorcode: nxx ', nxx, nx
        stop
    endif
    if (ny .gt. nyx) then
        write(*,*) 'Y grid size exceeds maximum of ', nyx
	    write(numRes,*) 'Y grid size exceeds maximum of ', nyx
	    write(numRes,*) 'errorcode: nyx ', nyx, ny
        stop
    endif
    
    read(2, *)
    read(2, *) nsx, nsy, isymmetry

    if (nsx.gt.1) then
	    read(2, *) !read end coordinate of each x grid section
	    read(2, *) (xnt(i), i = 2, nsx)
	    read(2, *) !read end grid number of each x grid section
	    read(2, *) (nxt(i), i = 2, nsx)
    else
	    read(2, *)
	    read(2, *)
	    read(2, *)
	    read(2, *)
    endif

    read(2, *) !read expansion ratio of each x grid section
    read(2, *) (dxr(i), i = 1, nsx)

    if (nsy.gt.1) then
	    read(2, *) !read end coordinate of each y grid section
	    read(2, *) (ynt(i), i = 2, nsy)
	    read(2, *) !read end grid number of each y grid section
	    read(2, *) (nyt(i), i = 2, nsy)
    else
	    read(2, *)
	    read(2, *)
	    read(2, *)
	    read(2, *)
    endif

!     **************Adaptive Grid****************
    read(2, *) !read expansion ratio of each y grid section
    read(2, *) (dyr(i), i = 1, nsy)
    read(2, *)
!	The following input is added by XJ. Shen on 7/29/03 to adapt the y grid separately.
    read(2, *) !read in adaptive grid params for the x direction
    read(2, *) difmax_x, decay_x, ipmax_x, iadaptMode_x
    read(2, *) !read in adaptive grid params for the y direction
    read(2, *) difmax_y, decay_y, ipmax_y  

!     ********Mesh Refinement in X Direction*******
!	Mesh refinement parameters
    read(2,*)
    read(2,*) !NumXRegions
    read(2,*)  nxmr
    read(2,*) !Lower Upper Amount
    if(nxmr.ge.1) then
	    do i=1,nxmr
	        read(2,*) xginit(i),xgend(i),rlevelx(i)
	    enddo
    endif
!    ********Mesh Refinement in Y Direction*******
    read(2,*)
    read(2,*) !NumYRegions
    read(2,*)  nymr
    read(2,*) !Lower Upper Amount
    if(nymr.ge.1) then
	    do i=1,nymr
	        read(2,*) yginit(i),ygend(i),rlevely(i)
	    enddo
    endif

!     ***************Reynolds Equation***************      
    read(2, *)
    read(2,'(a250)') header
    call parse(header, words, nwords)
    mg_nest = 4
    accom = 1.0
	if (nwords .lt.6) then
        read(2, *) idisc, iqpo, akmax, slip_beta, slip_gamma
	else if (nwords .lt.7) then
        read(2, *) idisc, iqpo, akmax, slip_beta, slip_gamma, accom
    else
        read(2, *) idisc, iqpo, akmax, slip_beta, slip_gamma,accom, max_mg_nest
	endif
	if (accom .lt. 0.99) then
	    write(*,*) 'WARNING: Accomodation coefficient reset to 1.0'
	endif
	
	mg_nest = max_mg_nest
	call SetGSIterations()
!	corCoef is the Surface Correction coefficient = (2-accommodation)/accommodation
	write(*,'(a21, f5.2)') ' Accommodation Coef: ', accom
	corCoef = (2 - accom)/accom

    if (iqpo .eq. 3) iqpo = 5 !FK model
    
!     ***************Partial Contact***************
    read(2, *)
!     temporary contact input
!	hamacker consts ahc and bhc added in CMLAir 6.20
    read(2,'(a250)') header
    call parse(header, words, nwords)
    if (nwords .lt. 4) then
	    read(2, *) icmod, rsik, cta
	else
	  read(2, *) icmod, rsik, cta, ahc, bhc
	endif

    !in the past icmod == 3 was used to turn IMF on, now icmod is only used to specify contact
    !IMF has been moved to additional parameters section
    if (icmod .eq. 3) then
        iUseIMF = 1
        icmod = 0
    endif

	!elecpot (electrostatic force) added 10/04
    read(2,'(a250)') header
    call parse(header, words, nwords)
    if (nwords .lt. 4) then
        elecpot = 0.0
        read(2, *) rasp, eyoung, ydst
    else
        read(2, *) rasp, eyoung, ydst, elecpot
    endif
    read(2, *)
    read(2, *) frcoe , pratio
    ydcoe=1.282d0 + 1.158d0 * pratio
    
    
!	***************Sensitivities***************
    read(2, *)
    read(2, *)
    read(2, *) crninc, cbrinc, twstinc
    read(2, *)
    read(2, *) tlnginc, tanginc, sldinc
    read(2, *)
    read(2, *) ptqinc, rtqinc, rcsinc
    read(2, *)
    read(2, *) iwscale

    call GetOptionalParameters(2)
    call CheckAlts()

    if (iUseIMF .eq. 1) then
        write(*,*)'INTERMOLECULAR FORCES ARE: ON'
    else
        write(*,*)'INTERMOLECULAR FORCES ARE: OFF'
    endif

    close(2)
    return

999 write(*,*)'trouble opening files.'
    write(numRes,*) 'trouble opening files.'
    stop
    end subroutine READRUN
    
    
!===================================================================
    subroutine CheckAlts()
!===================================================================
    use OptionalParams
    use Q4_globals
    implicit none

    !if run.dat has a specified non-zero altitude (we'll assume it's positive)
    !and the user wanted to 
    if (useStandardAir .eq. .false. .and. ialt .gt. 0) then
        !specified 1 altitude and it's not zero or more than 1 altitude
        if ((ialt .eq. 1 .and. abs(alts(1)) .gt. 0.01) .or. ialt .gt. 1) then
            write(*,*) 'Warning: cannot perform altitude simulation using'
            write(*,*) 'user defined air parameters (pressure, mean free path, viscosity)'
            write(*,*) 'continuing simulation as if no altitudes were specified'
            ialt = 0
        endif
    endif
    
    end subroutine CheckAlts
    
    
!===================================================================
    subroutine InitOptionalParams
!===================================================================
    use OptionalParams
    use Q4_globals
    use TemperatureHumidity
    implicit none

   	iOutputShear = 0
	iUseIMF = 0
	
	doHumidity = .false.
	humidity = 50.d0
    temperature = 25.d0
    
    useStandardAir = .true.
	
    useInitialPressure = .false.
    return
    end    


!===================================================================
    subroutine ProcessParameter(strParamName, paramValue)
!===================================================================
    use OptionalParams
    use Q4_globals
    use TemperatureHumidity
    implicit none
    
    character *PARAMSIZE :: strParamName
    real(8) :: paramValue, temp
    
    if (INDEX(strParamName, 'output_shear') .ne. 0) then
        iOutputShear = paramValue
    else if (INDEX(strParamName, 'use_imf') .ne. 0) then
        iUseIMF = paramValue
        if (iUseIMF .eq. 1 .and. icmod .ne. 0) then
            bhc = 0.d0
        endif
    else if (INDEX(strParamName, 'do_humidity') .ne. 0) then
        doHumidity = .false.
        if (paramValue .ne. 0.d0) doHumidity = .true.
    else if (INDEX(strParamName, 'humidity') .ne. 0) then
        humidity = paramValue
    else if (INDEX(strParamName, 'temperature') .ne. 0) then
        temperature = paramValue
    else if (INDEX(strParamName, 'standard_air') .ne. 0) then
        temp = paramValue
        if (paramValue .eq. 0) then
            useStandardAir = .false.
        else
            useStandardAir = .true.
        endif
    
        
!    else if (INDEX(strParamName, 'use_init_press') .ne. 0) then
!        useInitialPressure = .true.
!    else if (INDEX(strParamName, 'stagger_shear') .ne. 0) then
!        staggerShear = .true.
    else
        write(*,*) 'Unknown parameters encoutered in input file:'
        write(*,*) trim(strParamName), ' = ', paramValue
    endif
    
    return
    end subroutine ProcessParameter

!===================================================================
    subroutine CheckInitialPressureFile
!===================================================================
    use Q4_globals
    implicit none
    
    useInitialPressure = .false.
    if (iOldGrid .eq. 0) goto 656

    OPEN(UNIT=73, ERR=656, FILE='press_init.dat', STATUS='OLD')
    close(73)
    useInitialPressure = .true.
    
656 continue 
    
    return
    end subroutine CheckInitialPressureFile


!===================================================================
    subroutine LoadInitialPressure
!===================================================================
    use Q4_globals
    use NestingInfo
    implicit none
    integer iSuccess

    if (useInitialPressure .eq. .true.) then
        !try and use initial pressure
        call LoadMatrix(p, nx, ny, 'press_init.dat', iSuccess)

        if (iSuccess .eq. 0) goto 654

        !restrict (interpolate) pressure to all of the multigrid levels
        call restrictPressure(p, p1, xref, yref, xref1, yref1, &
                              nx, ny, nx1, ny1, enclosingFineRectX1, enclosingFineRectY1)
        call restrictPressure(p1, p2, xref1, yref1, xref2, yref2, &
                              nx1, ny1, nx2, ny2, enclosingFineRectX2, enclosingFineRectY2)
        call restrictPressure(p2, p3, xref2, yref2, xref3, yref3, &
                              nx2, ny2, nx3, ny3, enclosingFineRectX3, enclosingFineRectY3)
        call restrictPressure(p3, p4, xref3, yref3, xref4, yref4, &
                              nx3, ny3, nx4, ny4, enclosingFineRectX4, enclosingFineRectY4)
        return
    endif

654 continue
    call matInit(p, nx, ny, nx, ny, 1.d0)
    call matInit(p1, nx1, ny1, nx1, ny1, 1.d0)
    call matInit(p2, nx2, ny2, nx2, ny2, 1.d0)
    call matInit(p3, nx3, ny3, nx3, ny3, 1.d0)
    call matInit(p4, nx4, ny4, nx4, ny4, 1.d0)
    return
    end subroutine LoadInitialPressure


!=================================== 
    subroutine readRail
!===================================
!   Read in Rail.dat file
    use Q4_globals
    implicit none
    integer :: i, j, k

    !first go through and allocate space to hold our rail data    
    call GetMaxRailPts
    call allocateRails
    
    !now read everything in
    open(1,err=999,file='rail.dat',status='old')

    read(1, *)
    read(1, *)
    read(1, *) xl, yl, zl

    !now for each rail
    read(1,*) numRails, numWalls

    do k=1,numRails
	    read(1,*) npoints(k),istep(k)

	    do j=1,npoints(k)
	        read(1,*) xrail(j,k),yrail(j,k),indexw(k,j)
	    end do

	    if(istep(k).eq.0) then
	        read(1,*) (hramp(i,k),i=1,3)
	    else
	        read(1,*) hramp(1,k)
	    endif
    end do		      ! k=1,numRails...

    if(numWalls .gt. 0) read(1, *)(nwpoint(i),i=1, numwalls)

    do k = 1, numWalls
	    read(1,*) (wpoint(k,i),i=1,nwpoint(k))
	    read(1,*) (wrecess(k,i),i=1,nwpoint(k))
    enddo

    ! read the different step height
    read(1,*) xt,ht,rebase  !taper length, taper angle and base recess
    read(1,*) crown,camber,twist
    read(1,*) (xintNew(i),i=1,numPOI)
    read(1,*) (yintNew(i),i=1,numPOI)

    close(1)
    
    call ReadUDG !user defined geometry
    
    return

999 write(*,*)'trouble opening file: rail.dat'
    stop
    end subroutine readRail


! used to figure out how much space to allocate for the rail and wall data
!======================================
    subroutine GetMaxRailPts()
!======================================
    use Q4_globals
    implicit none
    
    integer :: i, j, k
    integer :: numRailPts
    integer :: iJunk, nwords
    character *4000 header
    !character *4000 words(110)
    character *4000 words(300)
    integer, dimension(:), allocatable :: numWallPts

    open(1,err=999,file='rail.dat',status='old')  

    read(1, *)
    read(1, *)
    read(1, *) !xl, yl, zl

    !read in the number of rails and walls
    read(1,*) numRails, numWalls
    
    do k = 1,numRails
	    read(1,*) numRailPts, iJunk  !numPointsOnRail and rail type (step or ramp)

        if (numRailPts > maxRailPts) then
            maxRailPts = numRailPts
        endif
        
        !now read in all the nodes for this rail and junk em
	    do j=1,numRailPts
	        read(1,*)
	    end do
        
        !read in the recess/step info
        read(1,*)
    end do ! k=1,numRails...

    ! read in the number of points in each wall
    if(numWalls .gt. 0) then
        allocate (numWallPts (numWalls))
        read(1, *)(numWallPts(i),i=1, numWalls)
    endif
    
    do k = 1, numWalls
        if (numWallPts(k) .gt. maxWallPts) then
            maxWallPts = numWallPts(k)
        endif
    enddo
    
    !finally pass by each wall, taper, baserecess and crown/camber info an see how many POIs we have
    do k = 1, numWalls
	    read(1,*)
	    read(1,*)
    enddo

    
    read(1,*)! pass by the taper, base reacess
    read(1,*)! pass by the crown camber and twist info

    !get number of Points of Interest    
    read(1, '(a4000)') header
    call parse(header, words, nwords)
    if (nwords > 250) then
          write(*,*) ' no. of POI''s exceeded the max value of 250'
          write(*,*) ' taking the first 250 points only'
    nwords = 250
    endif
    numPOI = nwords
    if (numPOI < 4) numPOI = 4
    
    
    close(1)
    return

999 write(*,*)'trouble opening file: rail.dat'
    stop
    end subroutine GetMaxRailPts
    
    
!======================================    
    subroutine checkRail()
!======================================
    !subroutine added to perform some error checking on the rail data
    !before allowing it to pass onto the main computations
    use Q4_globals
    implicit none
    
    integer :: k, j
    
    !make sure each wall index for each rail doesn't 
    !point to a wall that doesn't exist
    !if it does, just replace the wall with no profile
    do k = 1, numRails
	    do j=1, npoints(k)
            if (indexw(k,j) .gt. numWalls) then
                indexw(k,j) = 0
            endif
        enddo
    enddo
    
    return
    end subroutine checkRail
    
    
!======================================    
    subroutine allocateRails()
!======================================
    
    use Q4_globals
    implicit none
    
    !numRails, maxRailPts
    
    !allocate rail info
    allocate (xrail (maxRailPts+2, numRails))
    allocate (yrail (maxRailPts+2, numRails))
    allocate (xdiff (maxRailPts+2, numRails))
    allocate (ydiff (maxRailPts+2, numRails))
    allocate (cramp (3,numRails))
    allocate (hramp (3, numRails))
    allocate (xpt1 (numRails))
    allocate (xpt2 (numRails))
    allocate (ypt1 (numRails))
    allocate (ypt2 (numRails))
    allocate (npoints (numRails))
    allocate (istep (numRails))
    
    !allocate wall info
    allocate (wpoint (numWalls, maxWallPts))
    allocate (wrecess (numWalls, maxWallPts))
    allocate (xwalli (maxRailPts+1, numRails))
    allocate (ywalli (maxRailPts+1, numRails))
    allocate (xwallo (maxRailPts+1, numRails))
    allocate (ywallo (maxRailPts+1, numRails))
    allocate (flush_lower_xwallo (maxRailPts, numRails))
    allocate (flush_lower_ywallo (maxRailPts, numRails))
    allocate (flush_upper_xwallo (maxRailPts, numRails))
    allocate (flush_upper_ywallo (maxRailPts, numRails))
    allocate (xw1 (maxRailPts, numRails))
    allocate (yw1 (maxRailPts, numRails))
    allocate (xw2 (maxRailPts, numRails))
    allocate (yw2 (maxRailPts, numRails))
    allocate (indexw (numRails, maxRailPts+1))
    allocate (nwpoint (numWalls))
    
    !allocate space for POIs
    allocate (xintNew (numPOI))
    allocate (yintNew (numPOI))
    allocate (hintNew (numPOI))
    
    return
    end subroutine allocateRails
    

!======================================
    subroutine colinear(k)
!======================================

    use Q4_globals
    implicit none
    
    integer :: n, j, k, nrind1, nrind2, nrind3
    real(8) :: a1, a2, b1, b2, determ0

! check for and remove colinear points
! added 5/25/99 by RDG

    !only remove if at least 3 points in rail
    if(npoints(k).gt.2) then
    ! loop through all points in rail, starting with 3rd
        j = 3
        do while(j.le.npoints(k))
            !make sure we haven't exceeded array bounds
            nrind1 = j
            nrind2 = j+1
            nrind3 = j+2
            if(nrind2.gt.npoints(k)) then
                nrind2 = 1
                nrind3 = 2
            endif
            if(nrind3.gt.npoints(k)) then
                nrind3 = 1
            endif

            !Are the current 3 points colinear?
            a1 = (yrail(nrind2,k) - yrail(nrind1,k))
            a2 = (yrail(nrind3,k) - yrail(nrind2,k))
            b1 = (xrail(nrind1,k) - xrail(nrind2,k))
            b2 = (xrail(nrind2,k) - xrail(nrind3,k))
            determ0 = a1*b2-a2*b1
            if(dabs(determ0).lt.1.e-08) then
                write(*, '(a34, i2)') 'Removing colinear point in rail ',k
                j = 3
                do n=nrind2,npoints(k)-1
                    xrail(n,k) = xrail(n+1,k)
                    yrail(n,k) = yrail(n+1,k)
                enddo
                npoints(k) = npoints(k)-1
            endif
            j = j+1
        enddo
10  endif

    return
    end subroutine colinear
    

!======================================
      subroutine normalize
!======================================

    use Q4_globals
    implicit none
    
    integer :: i, ip1, j, k
    real(8) :: a(3,3), b(3)

    yl = yl / xl
    zl = zl / xl
    
    xf0 = xf0 / xl
    yf0 = yf0 / xl

    xfs = xfs * 1.d-6 / xl / 9.81d0  !convert to length-normalized kg
    yfs = yfs * 1.d-6 / xl / 9.81d0
    ptqinc = ptqinc * 1.d-6 / xl / 9.81d0
    rtqinc = rtqinc * 1.d-6 / xl / 9.81d0

    do k=1,numRails
	    do j=1,npoints(k)
	        xrail(j,k) = xrail(j,k) / xl
	        yrail(j,k) = yrail(j,k) / xl
	    end do

        call colinear(k)

	    xpt1(k)=1.d0
	    xpt2(k)=0.d0
	    ypt1(k)=1.d0
	    ypt2(k)=0.d0
	    do j=1,npoints(k)
	        xpt1(k)=dmin1(xpt1(k),xrail(j,k))
	        xpt2(k)=dmax1(xpt2(k),xrail(j,k))
	        ypt1(k)=dmin1(ypt1(k),yrail(j,k))
	        ypt2(k)=dmax1(ypt2(k),yrail(j,k))
	    end do
    end do		      ! k=1,numRails...

    do k=1,numWalls
	    do i=1, nwpoint(k)
	        wpoint(k,i)  = wpoint(k,i)  / xl
	        wrecess(k,i) = wrecess(k,i) / hm
	    end do
    enddo

    do i = 1, numPOI
	   xintNew(i) = xintNew(i) / xl
	   yintNew(i) = yintNew(i) / xl
	enddo

    ht = xt * ht
    ht = ht / hm
    xt = xt / xl
    tlnginc = tlnginc / xl

    crown  = crown  / hm
    crninc = crninc / hm
    camber = camber / hm
    cbrinc = cbrinc / hm
    twist  = twist  / hm
    twstinc= twstinc/ hm
    rebase = rebase / hm
    rcsinc = rcsinc / hm

    al = al / hm

    do k = 1,numRails
        if(istep(k).eq.0)then
	        a(1,1) = xrail(1,k)
	        a(1,2) = yrail(1,k)
	        a(1,3) = 1.d0
	        a(2,1) = xrail(2,k)
	        a(2,2) = yrail(2,k)
	        a(2,3) = 1.d0
	        a(3,1) = xrail(3,k)
	        a(3,2) = yrail(3,k)
	        a(3,3) = 1.d0
	        do i =1,3
	            hramp(i,k) = hramp(i,k)/hm
	            b(i) = hramp(i,k)
	        enddo
	        call matrix33_inverse(a)
	        call matrix_multi(3,3,1,a,b,cramp(1,k))
        else
	        hramp(1,k)=hramp(1,k)/hm
        endif

        ! calculate the slope of each lines
        do i = 1,npoints(k)
	        ip1 = i+1
	        if (i.eq.npoints(k)) ip1 = 1
	        xdiff(i,k) = xrail(ip1,k)-xrail(i,k)
	        ydiff(i,k) = yrail(ip1,k)-yrail(i,k)
        end do
    enddo   !numRails

    do i = 2, nsx
        xnt(i) = xnt(i) / xl
    enddo

    do i = 2, nsy
        ynt(i) = ynt(i) / xl
    enddo

    !computes the wallprofile for each rail
    call wallprofile

    !get grid control points
    call getControl

    return
    end subroutine normalize


!=========================================
	subroutine SetGSIterations
!=========================================
    !sets the number of Line by Line Gauss-Seidel iterations
    !used to invert the matrix in reyneq.
    !choices are based on speed tests conducted on a few sliders
    !at various values for mg_nest
    
    use Q4_globals
    implicit none
    
    if (mg_nest .eq. 4) then
        GS_Iterations = 3
    elseif (mg_nest .eq. 3) then
        GS_Iterations = 4
    elseif (mg_nest .eq. 2) then
        GS_Iterations = 4
    elseif (mg_nest .eq. 1) then
        GS_Iterations = 6
    else
        GS_Iterations = 10
    endif

    return
    end subroutine SetGSIterations
    
    
!=========================================
	subroutine IncreaseNesting
!=========================================
    use Q4_globals
    implicit none
    
    if (mg_nest .lt. max_mg_nest) then
        write(*,*) 'increasing nesting'
        mg_nest = max_mg_nest
    endif
    
    call SetGSIterations
    
    return
    end subroutine IncreaseNesting
    
    
!=========================================
	subroutine ReduceNesting
!=========================================
    use Q4_globals
    implicit none

    if (mg_nest .gt. 1) then
        
        mg_nest = mg_nest - 1
        write(*,*)'Automatically reducing multigrid nesting to', mg_nest
        
        call SetGSIterations
    endif

    return
    end subroutine ReduceNesting
    
    
!============================================
    subroutine allocateRun
!============================================

    use Q4_globals
    use NestingInfo
    use ShearArrays
    use TemperatureHumidity
    implicit none
    
    integer :: i, j, k

    nx1 = nx/2+1
    ny1 = ny/2+1
    
    nx2 = nx1/2+1
    ny2 = ny1/2+1
    
    nx3 = nx2/2+1
    ny3 = ny2/2+1
    
    nx4 = nx3/2+1
    ny4 = ny3/2+1

    !pressure
    allocate (p (nx, ny))
    allocate (p1 (nx1, ny1))
    allocate (p2 (nx2, ny2))
    allocate (p3 (nx3, ny3))
    allocate (p4 (nx4, ny4))
    
    allocate (savedPressure (nx, ny))
    
    !van der walls forces
    allocate (vdwMolecularForceMap (nx, ny))
    
    !reynolds equation residuals
    allocate (res (nx, ny))
    allocate (res1 (nx1, ny1))
    allocate (res2 (nx2, ny2))
    allocate (res3 (nx3, ny3))
    allocate (res4 (nx4, ny4))
    allocate (su01 (nx1, ny1))
    allocate (su02 (nx2, ny2))
    allocate (su03 (nx3, ny3))
    allocate (su04 (nx4, ny4))
    
    
    do j = 1, ny
        do i = 1, nx
            res(i, j) = 0.0
        enddo
    enddo
    do j = 1, ny1
        do i = 1, nx1
            res1(i, j) = 0.0
        enddo
    enddo
    do j = 1, ny2
        do i = 1, nx2
            res2(i, j) = 0.0
        enddo
    enddo
    do j = 1, ny3
        do i = 1, nx3
            res3(i, j) = 0.0
        enddo
    enddo
    do j = 1, ny4
        do i = 1, nx4
            res4(i, j) = 0.0
        enddo
    enddo
    
    
    !grids
    allocate (xref (nx))
    allocate (xref1 (nx1))
    allocate (xref2 (nx2))
    allocate (xref3 (nx3))
    allocate (xref4 (nx4))
    allocate (yref (ny))
    allocate (yref1 (ny1))
    allocate (yref2 (ny2))
    allocate (yref3 (ny3))
    allocate (yref4 (ny4))
    
    !bearing numbers
    allocate (bearx (nx, ny))
    allocate (bearx1 (nx1, ny1))
    allocate (bearx2 (nx2, ny2))
    allocate (bearx3 (nx3, ny3))
    allocate (bearx4 (nx4, ny4))
    allocate (beary (nx, ny))
    allocate (beary1 (nx1, ny1))
    allocate (beary2 (nx2, ny2))
    allocate (beary3 (nx3, ny3))
    allocate (beary4 (nx4, ny4))
    
    !shear data
    allocate (PoisilleShearX(nx-1, ny))
    allocate (CouetteShearX(nx-1, ny))
    allocate (PoisilleShearY(nx, ny-1))
    allocate (CouetteShearY(nx, ny-1))
    
    !height data???
    ! I can't remember the difference between h and hnew
    ! hnew is slider plane surface
    allocate (h (nx, ny))
    allocate (h1 (nx1, ny1))
    allocate (h2 (nx2, ny2))
    allocate (h3 (nx3, ny3))
    allocate (h4 (nx4, ny4))
    allocate (hnew (nx, ny))
    allocate (hnew1 (nx1, ny1))
    allocate (hnew2 (nx2, ny2))
    allocate (hnew3 (nx3, ny3))
    allocate (hnew4 (nx4, ny4))
    
    !various parameters needed for the reynolds equation
    allocate (cohimx (nx, ny))
    allocate (cohimx1 (nx1, ny1))
    allocate (cohimx2 (nx2, ny2))
    allocate (cohimx3 (nx3, ny3))
    allocate (cohimx4 (nx4, ny4))
    allocate (cohjmx (nx, ny))
    allocate (cohjmx1 (nx1, ny1))
    allocate (cohjmx2 (nx2, ny2))
    allocate (cohjmx3 (nx3, ny3))
    allocate (cohjmx4 (nx4, ny4))
    
    allocate (himax (nx, ny))
    allocate (himax1 (nx1, ny1))
    allocate (himax2 (nx2, ny2))
    allocate (himax3 (nx3, ny3))
    allocate (himax4 (nx4, ny4))
    allocate (himin (nx, ny))
    allocate (himin1 (nx1, ny1))
    allocate (himin2 (nx2, ny2))
    allocate (himin3 (nx3, ny3))
    allocate (himin4 (nx4, ny4))
    
    allocate (hjmax (nx, ny))
    allocate (hjmax1 (nx1, ny1))
    allocate (hjmax2 (nx2, ny2))
    allocate (hjmax3 (nx3, ny3))
    allocate (hjmax4 (nx4, ny4))
    allocate (hjmin (nx, ny))
    allocate (hjmin1 (nx1, ny1))
    allocate (hjmin2 (nx2, ny2))
    allocate (hjmin3 (nx3, ny3))
    allocate (hjmin4 (nx4, ny4))
    
    allocate (recssi (nx, ny))
    allocate (recssi1 (nx1, ny1))
    allocate (recssi2 (nx2, ny2))
    allocate (recssi3 (nx3, ny3))
    allocate (recssi4 (nx4, ny4))
    allocate (recssj (nx, ny))
    allocate (recssj1 (nx1, ny1))
    allocate (recssj2 (nx2, ny2))
    allocate (recssj3 (nx3, ny3))
    allocate (recssj4 (nx4, ny4))
    
    !multigrid nesting info
    allocate (enclosingCoarseRectX (nx))
    allocate (enclosingCoarseRectY (ny))
    allocate (enclosingCoarseRectX1 (nx1))
    allocate (enclosingCoarseRectY1 (ny1))
    allocate (enclosingCoarseRectX2 (nx2))
    allocate (enclosingCoarseRectY2 (ny2))
    allocate (enclosingCoarseRectX3 (nx3))
    allocate (enclosingCoarseRectY3 (ny3))
    !holds the coords of the lower x and y indexes of the rectangle on the finer mesh that
    !encloses the coarser mesh point
    allocate (enclosingFineRectX1 (nx1))
    allocate (enclosingFineRectY1 (ny1))
    allocate (enclosingFineRectX2 (nx2))
    allocate (enclosingFineRectY2 (ny2))
    allocate (enclosingFineRectX3 (nx3))
    allocate (enclosingFineRectY3 (ny3))
    allocate (enclosingFineRectX4 (nx4))
    allocate (enclosingFineRectY4 (ny4))
    
    call getNumRailNodes
    allocate (xControl (maxControlPts))
    allocate (yControl (maxControlPts))
    
    do j = 1, ny
        do i = 1, nx
            vdwMolecularForceMap(i,j) = 0
        enddo
    enddo
                
    do j = 1, ny
        do i = 1, nx-1
            PoisilleShearX(i, j) = 0
            CouetteShearX(i, j) = 0
        enddo
    enddo

    do j = 1, ny-1
        do i = 1, nx
            PoisilleShearY(i, j) = 0
            CouetteShearY(i, j) = 0
        enddo
    enddo

    
    return
    end subroutine allocateRun
    
    
!============================================
    subroutine getNumRailNodes()
!============================================
    use Q4_globals
    implicit none
    
    integer :: numNodes, i, j
    
    numNodes = 0
    
    do i = 1, numRails
        do j = 1, npoints(i)
            numNodes = numNodes+1
        enddo
    enddo
    maxControlPts = int((numNodes + 2) * 1.5)  !a little breathing room, just in case
    
    return
    end subroutine getNumRailNodes


!!====================================================
!SUBROUTINE ReadLineChars(iunit,line,length,icase)
!
!!    USE IO_stats
!
!    IMPLICIT NONE
!
!    ! Arguments: iunit  <-> IO unit
!    !            line   <-> string
!    !            length <-> actual length ( not declared )
!    !            icase  <-> report what happened
!
!    INTEGER, INTENT(IN) :: iunit
!    CHARACTER(LEN=*) :: line
!    INTEGER, INTENT(OUT) :: length, icase
!
!    ! Local variables
!
!    CHARACTER(LEN=1) :: char
!    INTEGER :: size
!
!    ! Get declared size of string
!
!    size= LEN(line)
!
!    ! Initially no characters have been read
!
!    length = 0
!
!    DO
!        READ(iunit, FMT='(A1)', ADVANCE='NO', IOSTAT=icase) char
!        IF (icase .GE. 0) THEN
!            length = length + 1
!            IF (length <= size) THEN
!                line(length:length) = char
!                CYCLE            ! there may be more to read
!            ELSE
!                icase = size+1   ! ERROR: line as too short string
!                EXIT
!            ENDIF
!        ELSE
!            EXIT
!        ENDIF
!            
!            
!!        SELECT CASE (icase)
!!        CASE(ok)  ! All OK
!!            length = length + 1
!!            IF (length <= size) THEN
!!                line(length:length) = char
!!                CYCLE            ! there may be more to read
!!            ELSE
!!                icase = size+1   ! ERROR: line as too short string
!!                EXIT
!!            ENDIF
!!        CASE(eor) ! End of record encountered
!!            EXIT   ! Record has been read
!!        CASE(eof) ! End of file encountered
!!            EXIT   ! No more records to read
!!        CASE DEFAULT ! I/O error has occured
!!            EXIT      ! Nothing more to be done
!!        END SELECT
!    ENDDO
!
!    END SUBROUTINE ReadLineChars       