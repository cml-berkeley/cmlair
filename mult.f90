!===============================================
    subroutine  setAlt(calt)
!===============================================

    use Q4_globals
    use TemperatureHumidity
    implicit none

    !Note: variable al is the Mean Free Path, not altitude
    real(8), parameter :: pi = 3.1415926535897932384626433832795d0
    real(8) :: calt, p00, ttt
    logical :: bChangedAirParameters
    real(8) :: g0, M, Rstar, R, stdDensity, S, MolecularDiameter, Na, B, tempKelvins
    
!    BCox - Much of this data is from :
!      * U.S. standard atmosphere, 1976 
!      * Allen's Astrophysical Quantities - chapter 11

    g0 = 9.80665d0
    M = 0.0289644d0
    Rstar = 8.31432d0
    R = 287.05d0
    tempKelvins = 273.15d0 + temperature
    stdDensity = 1.2928d0
    S = 110.4d0
    molecularDiameter = 3.65d-10
    Na  = 6.02213d23
    B = 1.458d-6
    if (useStandardAir .eq. .true.) then
        ! our pressure calculation should take temperature lapse rate into account
        ! if our drive is not sealed, the pressure in the drive is equal to the pressure outside the drive
        ! the pressure outside the drive is the pressure of the outside world.  Best way to express that
        ! is using the barometric formula that includes the temperature lapse rate
        p00   = 1.01325d+05
        ttt   = 288.16d0 + (-0.0065d0 * calt)
        p0    = p00 * (ttt/288.16d0)**5.256d0
        
        !p0 = 101325.d0 * exp((-g0 * M * calt) / (Rstar * tempKelvins))
        vis1 = (B * tempKelvins**1.5) / (tempKelvins + 110.4)
        al = (Rstar * tempKelvins) / (dsqrt(2.d0) * pi * molecularDiameter*molecularDiameter * Na * p0);
        al = al / hm
    endif
!
!    if(ialt .gt. 0 .and. (calt .gt. 0.001 .or. calt .lt. -0.001)) then
!        write(*,'(a22,f10.5)') ' ALTITUDE           = ', calt
!!        need to obtain air parameters
!
!        !we should warn the user that we're overwriting the standard air parrameters
!        !people disagree on the standard value for some of these so we'll allow a bit of
!        !leeway before warning the user
!        bChangedAirParameters = .false.
!        if ( abs(p0 - 1.0135d+05) .gt. 26 ) bChangedAirParameters = .true.
!        if ( abs(vis1 - 0.1806d-04) .gt. 0.000001 ) bChangedAirParameters = .true.
!        if ( abs(al - 63.5d-9/hm) .gt. (3.0d-9/hm) ) bChangedAirParameters = .true.
!        if (bChangedAirParameters .eq. .true.) then
!            write(*,*) ' WARNING:'
!            write(*,*) ' Default air parameters will be used for altitudes different than 0'
!        endif
!        
!        p00   = 1.0135d+05
!        ttt   = 288.16d0 - 0.0065d0 * calt
!
!        vis1  = 0.1806d-04
!        p0    = p00 * (ttt/288.16d0)**5.256d0
!        al    = 52.6202393d0 * vis1 * dsqrt(ttt/twopi) / p0
!        al    = al / hm
!    else
!        write(*,*) 'NO ALTITUDE GIVEN. AIR PARAMETERS USED.'
!    endif

    !D0 = 13.657*Pa*h0  (Pa is atmospheric pressure and h0 is characteristic length)

    d0 = 0.8862269254527d0/al     !BC - this is what Lin Wu refered to as the "modified inverse knudsen number"
    !beginning coefficient is 1/2 the square root of pi.  It's required for all slip models.  
    ! See 1990 F-K paper in Journal of Tribology.  
    ! Refered to as 'D'.  They simply refer to it as the "Inverse Knudsen number"
    
!   Modified by Lion Huang. Add new slip models in Quick4.
!   if(iqpo.eq.1) then
    if (iqpo.eq.1 .or. iqpo .eq. 11) then !first order slip
	    t1 = 6.d0*al
    else if (iqpo.eq.2 .or. iqpo.eq.15 .or. iqpo.eq.22 .or. iqpo.eq.6 .or. iqpo.eq.8) then !2nd, 1.5 continuum	
	    t1 = 6.d0*al											!pressure grad.
	    t2 = 6.d0*al*al										!and user def.
    else if (iqpo.eq.3 .or. iqpo.eq.5) then			!FK
	    t1 = 6.81971863421d0*al
	    t2 = 8.31137590768d0*al*al
	    t3 = 0.d0
	    t4 = -10.9086332458d0*al*al*al*al
    end if

    f2p=p0*hm*hm/xl/xl/vis1/12.d0
    p0xl=p0*xl*xl/9.81d0

    uflag = ON

    return
    end

!========================================
    subroutine setRPM(crpm)
!========================================
    use Q4_globals
    implicit none

    real(8) :: crpm
    
    rpm = crpm
    u0 = twopi * ra * crpm / 60.d0

    uflag = ON

    write(*,'(a22,g12.6)') ' ROTATION PER MINUTE= ', crpm
    write(*,'(a22,f10.5)') ' DISK VELOCITY(M/S) = ', u0

    return
    end
    
!========================================
    subroutine setRadius(radius, skew)
!========================================

    use Q4_globals
    implicit none
    
    real(8) :: radius, skew

    ra = radius
    ske = skew * twopi / 360.d0

    u0 = twopi * ra * rpm / 60.d0

    uflag = ON


    return
    end
    
!===================================================================
    subroutine traverse
!===================================================================
    !goes through each radial position/skew and calls doCase
    use Q4_globals
    implicit none
    
    integer :: i

    !first traversal of radial positions
    if(tflag .eq. ON) then
        icurrpm = 1
        if(ialt.eq.0) then
            icuralt = 0
        else
            icuralt = 1
        endif
        icursen = 0
    endif

    do i = 1, irad
        icurrad = i
        call setRadius(radii(i), skews(i))
        write(*,*)
        write(*,'(a22,f10.5)') ' RADIAL POSITION(MM)= ',ra * 1000.d0
        write(*,'(a22,f10.5)') ' SKEW ANGLE(DEGREES)= ',ske * 360.d0 / twopi
        write(*,*)

        call doCase(i)
        uflag = ON
        hflag = OFF
    enddo

    icurrad = 1

    !first traversal done, base case
    tflag = OFF
        
    return
    end
    
!======================================================
    subroutine doAlts
!======================================================
    use Q4_globals
    implicit none
    
    integer :: i

    if(ialt.lt.2) return
    write(*,*)
    write(*,*)'RUNNING AT DIFFERENT ALTITUDES...'
    write(*,*)

    do i = 2, ialt
        icuralt = i  
        write(*,*)'CURRENT ALTITUDE ',alts(i),' (M)' 
	    call setAlt(alts(i))
	    call traverse
    enddo
    icuralt = 1

    return
    end
    
!=================================
    subroutine senstv
!=================================

    use Q4_globals
    implicit none

    !No sensitivity for fixed attitude
	if(isolv .eq. 0) return

    call crwnss
    call cmbrss
    call twstss
    call loadss
    call ptquss
    call rtquss
    call tlngss
    call tangss
    call rcssss

    return
    end

!===================================
    subroutine crwnss
!===================================
    use Q4_globals
    implicit none
    
    real(8) :: crnold

    if(crninc .eq. 0) return

    write(*,*)
    write(*,*)'Now doing CROWN tolerence...'

    icursen = -1

    crnold = crown
    crown = crnold - crninc

    !debug
    write(*,*)
    write(*,'(a23, f5.1, a5)') ' Original CROWN      : ', crnold*hm*1.d9, ' (NM)'
    write(*,'(a23, f5.1, a5)') ' After negative. INC.: ', crown*hm*1.d9, ' (NM)'
    write(*,*)


    hflag = 1
    call traverse

    icursen = 1
    crown = crnold + crninc
    hflag = ON

    !debug
    write(*,*)
    write(*,'(a23, f5.1, a5)') ' Original CROWN      : ', crnold*hm*1.d9, ' (NM)'
    write(*,'(a23, f5.1, a5)') ' After positive. INC.: ', crown*hm*1.d9, ' (NM)'
    write(*,*)

    call traverse

    crown = crnold
    hflag = ON

    icursen = 0
    return
    end


!===================================
    subroutine cmbrss
!===================================
    use Q4_globals
    implicit none
    
    real(8) :: cbrold

    if(cbrinc .eq. 0) return

    write(*,*)
    write(*,*)'Now doing CAMBER tolerence...'

    icursen = -2

    cbrold = camber
    camber = cbrold - cbrinc
    hflag = ON

    !debug
    write(*,*)
    write(*,'(a23, f5.1, a5)') ' Original CAMBER     : ', cbrold*hm*1.d9, ' (NM)'
    write(*,'(a23, f5.1, a5)') ' After negative. INC.: ', camber*hm*1.d9, ' (NM)'
    write(*,*)

    call traverse

    icursen = 2
    camber = cbrold + cbrinc

    !debug
    write(*,*)
    write(*,'(a23, f5.1, a5)') ' Original CAMBER     : ', cbrold*hm*1.d9, ' (NM)'
    write(*,'(a23, f5.1, a5)') ' After positive. INC.: ', camber*hm*1.d9, ' (NM)'
    write(*,*)

    hflag = ON
    call traverse

    camber = cbrold
    hflag = ON

    icursen = 0
    return
    end

!===================================
    subroutine twstss
!===================================
    use Q4_globals
    implicit none
    
    real(8) :: twtold

    if(twstinc .eq. 0) return

    write(*,*)
    write(*,*)'Now doing TWIST tolerence...'

    icursen = -3

    twtold = twist
    twist = twtold - twstinc
    hflag = ON

    !debug
    write(*,*)
    write(*,'(a23, f5.1, a5)') ' Original TWIST      : ', twtold*hm*1.d9, ' (NM)'
    write(*,'(a23, f5.1, a5)') ' After negative. INC.: ', twist*hm*1.d9, ' (NM)'
    write(*,*)

    call traverse

    icursen = 3
    twist = twtold + twstinc

    !debug
    write(*,*)
    write(*,'(a23, f5.1, a5)') ' Original TWIST      : ', twtold*hm*1.d9, ' (NM)'
    write(*,'(a23, f5.1, a5)') ' After positive. INC.: ', twist*hm*1.d9, ' (NM)'
    write(*,*)

    hflag = ON

    call traverse

    twist = twtold
    hflag = ON

    icursen = 0
    return
    end

!===================================
    subroutine tlngss
!===================================
    use Q4_globals
    implicit none
    
    real(8) :: xtold, xttemp, htold

    if(tlnginc .eq. 0) return
    if(tlnginc .gt. 0.5d0*xt) then
        write(*,*)
        write(*,*) 'TAPER LENGTH TOLERENCE TOO LARGE.'
        write(*,*)
        return
    endif

    write(*,*)
    write(*,*) 'Now doing TAPER LENGTH tolerence...'

    icursen = -7

    xtold = xt
    htold = ht

    xt = xtold - tlnginc
    ht = htold * xt / xtold

    hflag = ON
    uflag = ON

    !adjust grid using new taper length
    !call xtGrid(xtold)

    !debug
    write(*,*)
    write(*,'(a23, f6.4, a5)') ' Original TAPER LENGTH: ', xtold*xl*1.d3, ' (MM)'
    write(*,'(a23, f6.4, a5)') ' After negative  INC.: ', xt*xl*1.d3, ' (MM)' 
    write(*,*)

    call traverse

    icursen = 7

    xttemp = xt
    xt = xtold + tlnginc
    ht = htold * xt / xtold

    hflag = ON
    uflag = ON
    !adjust grid using new taper length
    !call xtGrid(xttemp)

    !debug
    write(*,*)
    write(*,'(a23, f6.4, a5)') ' Original TAPER LENGTH: ', xtold*xl*1.d3, ' (MM)'
    write(*,'(a23, f6.4, a5)') ' After positive INC.: ', xt*xl*1.d3, ' (MM)' 
    write(*,*)

    call traverse

    xttemp = xt
    xt = xtold
    ht = htold
    !adjust grid using new taper length
    !call xtGrid(xttemp)

    hflag = ON
    uflag = ON

    icursen = 0
    return
    end
    
!=======================
    subroutine tangss
!=======================

    use Q4_globals
    implicit none
    
    real(8) :: htinc, htold

    if (tanginc .eq. 0) return

    htinc = tanginc * xt * xl / hm 

    if(htinc.ge.ht) then
        write(*,*)
        write(*,*)'TAPER ANGLE TOLERENCE TOO LARGE.'
        write(*,*)
    endif

    write(*,*)
    write(*,*)'Now doing TAPER ANGLE tolerence...'

    icursen = -8

    htold = ht

    ht = htold - htinc
    hflag = ON

    !debug
    write(*,*)
    write(*,'(a23, f6.3, a7)') ' Original TAPER ANGLE: ', htold*hm/xt/xl*1000.d0, ' (MRAD)'
    write(*,'(a23, f6.3, a7)') ' After negative  INC.: ', ht*hm/xt/xl*1000.d0, ' (MRAD)' 
    write(*,*)

    call traverse

    icursen = 8
    ht = htold + htinc 
    hflag = ON

    !debug
    write(*,*)
    write(*,'(a23, f6.3, a7)') ' Original TAPER ANGLE: ', htold*hm/xt/xl*1000.d0, ' (MRAD)'
    write(*,'(a23, f6.3, a7)') ' After positive  INC.: ', ht*hm/xt/xl*1000.d0, ' (MRAD)' 
    write(*,*)

    call traverse

    ht = htold
    hflag = ON

    icursen = 0
    return
    end

!======================
    subroutine loadss
!======================
    use Q4_globals
    implicit none
    
    real(8) :: sldold

    if(sldinc .eq. 0) return

    write(*,*)
    write(*,*)'Now doing LOAD tolerence...'

    icursen = -4

    sldold = f0 
    f0 = sldold - sldinc

    !debug
    write(*,*)
    write(*,'(a23, f6.3, a7)') ' Original LOAD       : ', sldold*1000.d0, ' (G)'
    write(*,'(a23, f6.3, a7)') ' After negative  INC.: ', f0*1000.d0, ' (G)' 
    write(*,*)

    call traverse

    icursen = 4
    f0 = sldold + sldinc

    !debug
    write(*,*)
    write(*,'(a23, f6.3, a4)') ' Original LOAD       : ', sldold*1000.d0, ' (G)'
    write(*,'(a23, f6.3, a4)') ' After negative  INC.: ', f0*1000.d0, ' (G)' 
    write(*,*)

    call traverse

    f0 = sldold
    icursen = 0
    return
    end

!======================
    subroutine ptquss
!======================
    use Q4_globals
    implicit none
    
    real(8) :: ptqold

    if(ptqinc .eq. 0) return

    write(*,*)
    write(*,*)'Now doing STATIC PITCH tolerence...'

    icursen = -5

    ptqold = xfs 
    xfs = ptqold - ptqinc

    !debug
    write(*,*)
    write(*,'(a23, f7.3, a7)') ' Original PTORQUE    : ', ptqold*xl*9.81d6, ' (uN-M)'
    write(*,'(a23, f7.3, a7)') ' After negative  INC.: ', xfs*xl*9.81d6, ' (uN-M)' 
    write(*,*)

    call traverse

    icursen = 5
    xfs = ptqold + ptqinc

    !debug
    write(*,*)
    write(*,'(a23, f7.3, a7)') ' Original PTORQUE    : ', ptqold*xl*9.81d6, ' (uN-M)'
    write(*,'(a23, f7.3, a7)') ' After positive  INC.: ', xfs*xl*9.81d6, ' (uN-M)' 
    write(*,*)

    call traverse

    xfs = ptqold
    icursen = 0
    return
    end

!======================
    subroutine rtquss
!======================
    use Q4_globals
    implicit none
    
    real(8) :: rtqold

    if(rtqinc .eq. 0) return

    write(*,*)
    write(*,*)'Now doing STATIC ROLL tolerence...'
    write(*,*)

    icursen = -6

    rtqold = yfs 
    yfs = rtqold - rtqinc

    !debug
    write(*,*)
    write(*,'(a23, f7.3, a7)') ' Original RTORQUE    : ', rtqold*xl*9.81d6, ' (uN-M)'
    write(*,'(a23, f7.3, a7)') ' After negative  INC.: ', yfs*xl*9.81d6, ' (uN-M)' 
    write(*,*)

    call traverse

    icursen = 6
    yfs = rtqold + rtqinc

    !debug
    write(*,*)
    write(*,'(a23, f7.3, a7)') ' Original RTORQUE    : ', rtqold*xl*9.81d6, ' (uN-M)'
    write(*,'(a23, f7.3, a7)') ' After positive  INC.: ', yfs*xl*9.81d6, ' (uN-M)' 
    write(*,*)

    call traverse

    yfs = rtqold
    icursen = 0
    
    return
    end

!======================
    subroutine rcssss
!======================
    use Q4_globals
    implicit none
    
    real(8) :: rcsold, retemp

    if(rcsinc .eq. 0) return

    write(*,*)
    write(*,*)'Now doing RECESS tolerence...'

    icursen = -9

    rcsold = rebase
    rebase = rcsold - rcsinc
    call chgWalls(rcsold, rebase)
    hflag = ON

    !debug
    write(*,*)
    write(*,'(a23, f6.3, a5)') ' Original RECESS     : ', rcsold*hm*1.d6, ' (UM)'
    write(*,'(a23, f6.3, a5)') ' After negative  INC.: ', rebase*hm*1.d6, ' (UM)' 
    write(*,*)

    call traverse

    icursen = 9

    retemp = rebase
    rebase = rcsold + rcsinc

    call chgWalls(retemp, rebase)
    hflag = ON

    !debug
    write(*,*)
    write(*,'(a23, f6.3, a5)') ' Original RECESS     : ', rcsold*hm*1.d6, ' (UM)'
    write(*,'(a23, f6.3, a5)') ' After positive  INC.: ', rebase*hm*1.d6, ' (UM)' 
    write(*,*)

    call traverse

    call chgWalls(rebase, rcsold)
    rebase = rcsold
    hflag = ON
    icursen = 0
    
    return
    end

!=====================================================
    subroutine doRpms
!======================================================
    use Q4_globals
    implicit none
    
    integer :: i

    if(irpm .lt. 2) return

    write(*,*)
    write(*,*)'Running at different RPMs: '
    write(*,*)
    if(ialt.ge.1) then
        write(*,*)'Resetting alititude to base'
        write(*,*)
        call setAlt(alts(1))
    endif

    do i = 2, irpm
        icurrpm = i 
        call setRpm(rpms(i))
        uflag = ON
        call traverse
    enddo

    icurrpm = 1 
    return
    end

!===================================
    subroutine doCase(ir)
!===================================
    use Q4_globals
    implicit none
    
    integer :: ir

    !this sucks
    !2 primary cases: one for flyheight solutions and one for static attitude solutions
    if (isolv .eq. 1) then
        !we need to adapt the grid to the approximate FH at the first radial position and
        !subsequent radial positions if the user wants a new grid at each radial position
        if ((iadpt .eq. 1) .and.    &
            ((iUseNewGridEachRun .eq. 1) .or. (iRunCounter .eq. 0))) then
            call calcReynoldsCoefs
            call approximateFHAdapt
        endif
        call calcReynoldsCoefs !no matter what, we need to call this after initCasePreAdapt
    elseif (isolv .eq. 0) then
        !if the user wants a new grid for each radial position then
        !we've already adapted for the first run, but subsequent runs need an adaptive grid
        if (iUseNewGridEachRun .eq. 1) then
            if (iRunCounter .gt. 0 .and. iadpt .eq. 1) then
                call calcReynoldsCoefs
                call initCasePreAdapt
            endif
        endif
        call calcReynoldsCoefs  !no matter what, we need to call this after initCasePreAdapt
    endif
    
    !if we already have an initial pressure file, don't do a full multigrid, just start on fine level
    if (iadpt .eq. 0 .and. useInitialPressure .eq. .false.) then
        call fullmult(1, mg_nest, 0)
    else
        call fullmult(0, mg_nest, 0)
    endif
    
    if (isolv .eq. ON) call inverse(1)

!   figure out if we should save pressure    
    if (isave .eq. 0) then  !don't save any pressure
        pflag = OFF
    else
        !save for only first run if isave == 1, save for all runs if isave == 2
        if (tflag .eq. ON) then
            pflag = ON
        else if (isave .eq. 2) then
            pflag = ON
        else
            pflag = OFF
        endif
    endif

    !first radial traversal, base case
	if(tflag.eq.ON) then
        sflag = istiff
	else
        sflag = OFF
    endif

	if(sflag .eq. ON) call stiffness


    !run finished, output data
    if (isolv .ne. ON) then
        !not a fly height calculation
        call getPint
	    call getMinH
        err=(dabs(f-f0)+dabs(xf-xfs)+dabs(yf-yfs))/f0  !this got moved here so we can output error for all cases
        call outputResult(numRes)
	endif
	    
	if (pflag .eq. ON) then
	    call mflow
	    call outputPressures(iRunCounter+1) 
    endif

    iRunCounter = iRunCounter + 1

    return
    end


!=====================================
	subroutine chgWalls(rold, rnew)
!=====================================

    use Q4_globals
    implicit none

    real(8) :: rold, rnew, fac
    integer :: i, ip, j

	do i = 1, numWalls
	    if(wrecess(i, nwpoint(i)) .eq. rold) then
            !only change the part outside nominal wall
	        j = 1
	        do while(wpoint(i, j) .le. 0.d0)
		        j = j + 1
	        enddo

	        fac = (rnew - wrecess(i,1)) / (rold - wrecess(i,1))

	        do ip = j, nwpoint(i) - 1
		        wrecess(i,ip) = wrecess(i,ip) * fac
	        enddo
	        wrecess(i, nwpoint(i)) = rnew
	        if(iwscale .eq. 1) then
		        do ip = j, nwpoint(i)
                    wpoint(i,ip) = wpoint(i,ip) * fac
		        enddo
	        endif
	    endif
	enddo

	call wallprofile

	return
	end

!===================================================
 !   subroutine endTime
!===================================================
 !   use timing_module
 !   implicit none
	
  !  CHARACTER (LEN = 12) REAL_CLOCK (3)
   ! CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2),REAL_CLOCK (3), end_time )
	
!	write(*,*)
!	write(*,10) 'SIMULATION STARTED AT: ', START_TIME(5), START_TIME(6), START_TIME(7)
!10	format(1X, A, I2.2, ':', I2.2, ':', I2.2) 
!	write(*,15) 'SIMULATION ENDED AT:   ', end_time(5), end_time(6), end_time(7)
!15	format(1X, A, I2.2, ':', I2.2, ':', I2.2) 
!	write(*,*)
	
!	return
!	end subroutine endTime
	
!! DATE_AND_TIME(DATE, TIME, ZONE, VALUES) gets the corresponding date and time information from the real-time system clock.
!DATE is INTENT(OUT) and has form ccyymmdd. TIME is INTENT(OUT) and has form hhmmss.sss. ZONE is INTENT(OUT) and has form (+-)hhmm, 
!representing the difference with respect to Coordinated Universal Time (UTC). Unavailable time and date parameters return blanks.

!VALUES is INTENT(OUT) and provides the following:

!VALUE(1):	The year 
!VALUE(2):	The month 
!VALUE(3):	The day of the month 
!VALUE(4):	Time difference with UTC in minutes 
!VALUE(5):	The hour of the day 
!VALUE(6):	The minutes of the hour 
!VALUE(7):	The seconds of the minute 
!VALUE(8):	The milliseconds of the second 

!============================
    subroutine tic(t)
!============================
! this subroutine gives you the count of milliseconds from the start of the day
! you run this code
 
CHARACTER (LEN = 12) REAL_CLOCK (3)
integer :: t ! where we save the output in milliseconds
integer :: str_time(8) !start time

CALL DATE_AND_TIME (REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3),str_time)


t = str_time(5)*60*60*1000+str_time(6)*60*1000+str_time(7)*1000+str_time(8)

return
end subroutine tic