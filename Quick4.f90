!  Quick4.f90 
!
!  FUNCTIONS:
!  Quick4      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Quick4
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
!
! file handles
!============================================================================
! 2 = run.dat
! 12 = numRes, result.dat
! 25 = output grid data
! 31 = grid data
! 41 = initial pressure matrix


!=============================================
    program Quick4
!==============================================
    use Q4_globals
    
    implicit none
    
    integer :: tss,tse,ts,tsp,tp
    
    ! for time calc purposes 
    ! added by Dolf 2/9/17
    open(5, err=999,file='renq_tm.dat',status='unknown')
    open(4, err=999,file='subs_tm.dat',status='unknown')
    
    call tic(tss)
    tsp = tss
    call initws
    call tic(tse)
    ts = tse-tss
    write(4,*)'initws time :', ts
    
    !call Create_FK_Tables
    call tic(tss)
    call create_dbase
    call tic(tse)
    ts = tse-tss
    write(4,*)'create_db time :', ts
    
    call tic(tss)
    call initialGridGen
    call tic(tse)
    ts = tse-tss
    write(4,*)'iniGrid time :', ts
    
    call tic(tss)
    call setUpInitialNesting
    call tic(tse)
    ts = tse-tss
    write(4,*)'setUpIniNest time :', ts
    
    !output initial grid for debugging purposes
    call tic(tss)
    call outputGrid
    call tic(tse)
    ts = tse-tss
    write(4,*)'opGrid time :', ts
    
    !simple gap extrapolation and Geometry calculations will call stop
    
    if (gapExtrapolationRun .eq. .true.) then
        call tic(tss)
        call doGapExtrapolationRun
        call tic(tse)
        ts = tse-tss
        write(4,*)'GapExtraRun time :', ts
    endif
    
    if (iRailGeom .eq. ON) then 
        call tic(tss)
        call doRailGeomRun
        call tic(tse)
        ts = tse-tss
        write(4,*)'doRailGeomRun time :', ts
    endif
    
    !flyheight or static attitude calculation
    if (iadpt.eq.ON) then 
        call tic(tss)
        call initCasePreAdapt
        call tic(tse)
        ts = tse-tss
        write(4,*)'iniCasePrAdpt time :', ts
    endif
    
    
    !call traverse for all cases except railgeom.  This is the heart of running simulations
    call tic(tss)
    call traverse
    call tic(tse)
    ts = tse-tss
    write(4,*)'traverse time :', ts
    
    call tic(tss)
    call outputGrid
    call tic(tse)
    ts = tse-tss
    write(4,*)'opGrid time :', ts
    
    !compute sensitivity to manufacturing tolerance
    call tic(tss)
    call senstv 
    call tic(tse)
    ts = tse-tss
    write(4,*)'senstv time :', ts
    
    !run through additional altitudes (2..n)
    call tic(tss)
    call doAlts
    call tic(tse)
    ts = tse-tss
    write(4,*)'doAlts time :', ts

    !run through additional RPMs (2..n)
    call tic(tss)
    call doRPMs
    call tic(tse)
    ts = tse-tss
    write(4,*)'doRPM time :', ts
    
    call tic(tss)
    call output_spacing_2
    call tic(tse)
    ts = tse-tss
    tp = tse-tsp
    write(4,*)'opSpace time :', ts
    
    call endTime()
    
    
    write(4,*)'quick4 time :', tp
    write(*,*)'quick4 time :', tp
    
    write(*,*)'PROGRAM FINISHED.'
    close(5)
    close(4)
    stop
999 write(*,*)'trouble opening files.'

end program Quick4



!============================================
    subroutine initCasePreAdapt
!============================================
    !find the pressure under the slider for the flying attitude specified by the user
    !adapts the grid to that pressure.  After we adapt the grid we go through and (for a standard run)
    !we call approximateFHAdapt.
    
    use Q4_globals
    implicit none

    if(iadpt .eq. ON)then
        hflag = ON
        uflag = ON
        call calcReynoldsCoefs
        
        !solve the problem at given flying atitude

        if (iRunCounter .eq. 0) then  !good ol' cosmetics: don't want to print this out twice
            write(*,'(a22,f10.5)')' RADIAL POSITION(MM)= ',ra*1000.d0
            write(*,'(a22,f10.5)')' SKEW ANGLE(DEGREES)= ',ske*360.d0/twopi
        endif
        write(*,*)
        write(*,*) 'COMPUTING PRESSURE FIELD ON THE CURRENT GRID...'
        write(*,*)
        if (useInitialPressure .eq. .true.) then
            call fullmult(0,mg_nest,0)
        else
            call fullmult(1,mg_nest,0)
        endif
        
        !resolve the problem on adaptive grid
        write(*,*)
        write(*,*) 'ADAPTING GRID TO PRESSURE AT INITAL ATTITUDE...'
        write(*,*)
        call adaptive
        hflag = ON
        uflag = ON
        call calcReynoldsCoefs
        
        call outputGrid
    endif

    return
end subroutine initCasePreAdapt


!=========================================
    subroutine approximateFHAdapt
!=========================================
    !computes the inverse solution with reduced load error tolerance 
    !then adapts the grid to the pressure found at the approximate flyheight

    use Q4_globals
    implicit none

    if (isolv .eq. ON .and. iadpt .eq. ON)then
        write(*,*)
        write(*,*) 'COMPUTING  PRESSURE FIELD ON ADAPTIVE GRID...'
        write(*,*)

        call fullmult(0, mg_nest, 0)
        call inverse(0)

        !if we crashed, just quit
        if(crash) then
            write(numRes,*) 'Slider Crashed'
            call outputCrashedResult(numRes)
            stop
        endif

        write(*,'(/,a60,/)') 'FIRST ROUND QUASI-NEWTON SEARCH FINISHED.'
        !resolve the problem on adaptive grid
        write(*,*)
        write(*,*) 'ADAPTING GRID TO PRESSURE AT APPROXIMATE FLY HEIGHT...'
        write(*,*)
        call adaptive   !don't call getHeight, ave_multi, or calc_bearing_number here
                        !they're called in doCase
        hflag = ON
        uflag = ON
    endif    

	return
    end subroutine approximateFHAdapt

!=========================================
    subroutine mflow
!=========================================
!wrapper routine for massflow subroutine

    use Q4_globals
    implicit none
    
    integer :: incom

    call massflow(bearx,beary,cohimx,cohjmx,f2p,h,himax, &
		          himin,hjmax,hjmin,hm,hnew,p,p0,res, &
     		      recssi,recssj,xl,xref,yref,incom,nx,ny,nx,ny)
	return
    end subroutine mflow


!========================================
    subroutine doRailGeomRun()
!========================================

    use Q4_globals
    implicit none    

    call getHeight
    call outputGrid
    
    if (irailgeom .eq. ON) then
        write(*,*)'GEOMETRY GENERATED.'

        call getPint
        call getMinH

        call outputPressures(1) !not sure why we call this...
        
        icurrad = 1
        icurrpm = 1
        
        call outputResult(numRes)
        
    endif
    
    call endTime()
    write(*,*) 'PROGRAM FINISHED.'
    write(numRes,*) 'From CMLAir, set solution type in the Options menu'
    stop
    
end subroutine doRailGeomRun
 

!========================================
    subroutine doGapExtrapolationRun()
!========================================
!Gap extrap moved to function 6/3/05
!not sure if anyone even uses this func...
!should actually be a function in the interface 
!since it doesn't require a large amount of 
!complicated computation...oh well.  

    use Q4_globals
    implicit none
    
    integer :: i
    real(8) :: c1, c2, c3, y1, y2, y3, a1, a2, a, b, c
    real(8) :: zxgap, zygap
    
    write(*,*)
    write(*,*) 'Gap Extrapolation'
    write(*,*)
    !Calculate recess values for the four 
    !points of interest
    do i=1,4
        call POINTRECESS(rint(i),xintNew(i),yintNew(i))
    enddo
    !Solving linear equation for hy, hx0, and h0
    c1 = hintge(1)/hm-rint(1)
    c2 = hintge(2)/hm-rint(2)
    c3 = hintge(3)/hm-rint(3)
    y1 = yintNew(1)-yl*0.5d0
    y2 = yintNew(2)-yl*0.5d0 
    y3 = yintNew(3)-yl*0.5d0
    a1 = (c3-c1)/(xintNew(3)-xintNew(1))-(c2-c1)/(xintNew(2)-xintNew(1))
    a2 = (y1-y3)/(xintNew(3)-xintNew(1))-(y1-y2)/(xintNew(2)-xintNew(1))
    a = a1/a2
    b = (y1-y2)/(xintNew(2)-xintNew(1))
    c = (c2-c1)/(xintNew(2)-xintNew(1))
    hy = a
	    hx0 = c-b*a
    h0 = c1-xintNew(1)*c+(xintNew(1)*b+y1)*a
    !write(numRes,'(4(a5,g12.4))') ' hm= ',hm,' h0= ',h0,' hx0=',hx0,' hy= ',hy
    zxgap=h0+hx0*xintNew(4)
    zygap=hy*(yintNew(4)-yl*0.5d0)
    !rint = rint(4)
    hgap = (rint(4)+zxgap-zygap)*hm*1d9
    !Transfer data to "result.dat" output file	
    write(numRes,*)
    write(numRes,*) 'Gap Extrapolation Solution (FH, Pitch, and Roll)'
    write(numRes,*)
    write(numres,'(a23,f12.4)') ' FH on 4th Point (nm)= ',hgap(1)
    write(numres,'(a15,f12.4)') ' Pitch (탍ads)=',-hx0*hm/xl*1d6
    write(numres,'(a14,f12.4)') ' Roll (탍ads)=',hy*hm/xl*1d6

    write(*,*)
    write(*,*) 'Gap Extrapolation Solution (FH, Pitch, and Roll)...'
    write(*,*)
    write(*,'(a23,f12.4)') ' FH on 4th Point (nm)= ',hgap(1)
    write(*,'(a15,f12.4)') ' Pitch (탍ads)=',-hx0*hm/xl*1d6
    write(*,'(a14,f12.4)') ' Roll (탍ads)=',hy*hm/xl*1d6
    write(*,*)
    write(*,*) 'Calculation Completed'
    call endTime()
    write(*,*)
    stop
end subroutine doGapExtrapolationRun

