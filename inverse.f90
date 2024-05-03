!===================================
      subroutine calc_deriv (incrementAmount)
!===================================

    use Q4_globals
    implicit none

    real(8), intent(IN) :: incrementAmount
    real(8) :: psave(nx,ny)
    real(8) :: delth, deltp, deltr, hmsave, hx0save, hysave, h0save
    real(8) :: g1, g2, g3, g1dh, g2dh, g3dh, g1dp, g2dp, g3dp, g1dr, g2dr, g3dr
    real(8), parameter :: finite_difference = 1.0d-1
    
    call fullmult(0, mg_nest, 0)
    
    !call getfxp(0)
    g1=f
    g2=xf
    g3=yf

    ! the increment could be adjusted if necessary
    

    if(f.gt.f0) then
        delth = incrementAmount
    else
        delth = -incrementAmount
    endif

    if(xf.gt.xfs) then
        deltp = incrementAmount
    else
        deltp = -incrementAmount
    endif

    if(yf.gt.yfs) then
        deltr = -incrementAmount
    else
        deltr = incrementAmount
    endif

    write(*,'(/)')
    write(6,*) 'Start forming JACOBIAN matrix...'

    hmsave  = hmin
    hx0save = hx0
    hysave  = hy
    h0save  = hmin - hx0  !I don't think this is ever used...

    !save current pressure...
    call EQUATE(nx,ny,nx,ny,psave,p)

    !first, calculate derivatives wrt. him...

    hmin = hmsave + delth
    hx0  = hx0save
    hy   = hysave
    h0   = hmin - hx0save

    write(*,'(/)')
    write(6,*) 'Calculating deriv. wrt. HEIGHT...'

    call gethx
    do while(crash)
        delth = delth / 2.d0
        hmin = hmsave + delth
        h0   = hmin - hx0save
        call gethx
    enddo

    call fullmult(0,mg_nest,0)

    g1dh=f
    g2dh=xf
    g3dh=yf

    !reset pressure and heights...
    call equate(nx,ny,nx,ny,p,psave)

    !next, calculate derivatives wrt. pitch...

    hmin = hmsave
    hx0  = hx0save + deltp 
    hy   = hysave
    h0   = hmsave-hx0

    write(*,'(/)')
    write(6,*) 'Calculating deriv. wrt. PITCH...'
    call gethx
    do while(crash)
        deltp = deltp / 2.d0
        hx0  = hx0save + deltp 
        h0   = hmsave-hx0
        call gethx
    enddo
    call fullmult(0,mg_nest,0)

    g1dp=f
    g2dp=xf
    g3dp=yf


    !reset pressure and heights...
    call equate(nx,ny,nx,ny,p,psave)

    !next, calculate derivatives wrt. roll...

    hmin = hmsave
    hx0  = hx0save
    hy   = hysave + deltr 
    h0   = hmin - hx0save

    write(*,'(/)')
    write(6,*) 'Calculating deriv. wrt. ROLL...'

    call gethx
    do while(crash)
        deltr = deltr / 2.d0
        hy = hysave + deltr 
        call gethx
    enddo
    
    call fullmult(0,mg_nest,0)

    g1dr=f
    g2dr=xf
    g3dr=yf

    !reset pressure and heights...
    call equate(nx, ny, nx, ny, p, psave)

    hmin = hmsave
    hx0  = hx0save
    hy   = hysave
    h0   = hmin - hx0save

    f  = g1
    xf = g2
    yf = g3

    !jacobi matrix

    jac(1, 1) = (g1dh - g1) / delth
    jac(1, 2) = (g1dp - g1) / deltp
    jac(1, 3) = (g1dr - g1) / deltr
    jac(2, 1) = (g2dh - g2) / delth
    jac(2, 2) = (g2dp - g2) / deltp
    jac(2, 3) = (g2dr - g2) / deltr
    jac(3, 1) = (g3dh - g3) / delth
    jac(3, 2) = (g3dp - g3) / deltp
    jac(3, 3) = (g3dr - g3) / deltr

    write (*, '(/)')
    write (6, *) 'Finished forming JACOBIAN matrix...'
    write (*, '(/)')

    return
    end subroutine calc_deriv


!==========================================
      subroutine inverse(iind)
!==========================================

    use Q4_globals
    use PressureWarnings
    implicit none
    
    integer :: iind, icrash, ktot, ilocal, max_line_searches, kreduce, i, j
    integer, parameter :: maxQNIterations = 50 !used to be 20 (8/18/05)
    real(8), parameter :: derivAKReduction = 4.0
    !integer :: numSearchFailures
    logical :: calcDerivInPrevIteration, bHasLastGood
    real(8) :: etemp, errold, fsave, xfsave, yfsave, hmsave, hx0save, hysave
    real(8) :: sn1, sn2, sn3, y1, y2, y3, s1, s2, s3
    !real(8) :: lastGoodIncHeight, lastGoodIncPitch, lastGoodIncRoll
    
    !save all previous good iterations
    real(8) :: fPrev(maxQNIterations), xfPrev(maxQNIterations), yfPrev(maxQNIterations)
    real(8) :: hminPrev(maxQNIterations), hx0Prev(maxQNIterations), hyPrev(maxQNIterations)
    real(8) :: errorPrev(maxQNIterations)
    real(8) lowestErrorsIndex(3), tmpLowestError
    integer :: tmpLowestErrorsIndex
    real(8) :: akMaxOld
    

    !real*8 psave(nx,ny)

!	etemp is the convergence criterion in this routine
!     BC -- iind gets set to 0 the first time we calculate the inverse solution
!		  we reduce the convergence criterion for the first time since we
!		  just want to figure out how much to adapt the grid
!		  (Remember, we call inverse twice.  Once to figure out an approx FH
!		   we then adapt the grid and call it again to get the "exact" FH)

    if(iind.eq.0 .and. iadpt.eq.1)then
        !reduce convergence requirement for the first equilibrium search
        etemp = emax * 4.d0
    else
        etemp=emax
    endif

    errold = (dabs(f-f0) + dabs(xf-xfs) + dabs(yf-yfs)) / f0

    !!!!!!Put in own subroutine
    fPrev(1) = f
    xfPrev(1) = xf
    yfPrev(1) = yf
    hminPrev(1) = hmin
    hx0Prev(1) = hx0
    hyPrev(1) = hy
    errorPrev(1) = errold
    !!!!!!!!!!!!!!!!!!!!!!!!!!!


    !is the convergence criterion met (load error and residual)
    if (errold.lt.etemp .and. ak.lt.akmax) then
        if(iind.gt.0 .or. iadpt.eq.0) call outputResult(numRes)
	    return
    endif


    akmax = akmax / derivAKReduction
    call calc_deriv(5.0d-2)
    akmax = akmax * derivAKReduction
    !call calc_deriv (1.0d-1)
    calcDerivInPrevIteration = .true.

    write(*,'(/,a40,/)')'Start QUASI-NEWTON iteration.'

    call saveState(fsave, xfsave, yfsave, hmsave, hx0save, hysave)

    ! icrash keeps track of CONSECUTIVE crash
	icrash = 0
	!numSearchFailures = 0

    !  -----start iteration	 --------------------------------
    do ktot = 1, maxQNIterations

        write(*,'(/a20,i10)') 'Iteration ',ktot
        
1234    call getInc(sn1, sn2, sn3)
        !write(*,*) 'WARNING: reducing inc by 2 off the bat'
        !call reduceInc(sn1, sn2, sn3)

	    call stateInc(hmsave, hx0save, hysave, sn1, sn2, sn3)

        call gethx

        ilocal=0
        do while(crash) 
	        ilocal=ilocal+1
	        crash=.false.
	        call reduceInc(sn1, sn2, sn3)
	        call stateInc(hmsave, hx0save, hysave, sn1, sn2, sn3)
            call gethx
	    enddo

	    if(ilocal.eq.0) then
            !restart crash counter
            icrash = 0
        else    
            icrash = icrash + 1
        endif

        !slider crash assumed
	    if(icrash.ge.8) goto 99
        max_line_searches = 2
	    do kreduce = 0, max_line_searches
	    
            call fullmult(0, mg_nest, 0)
            
            err = ( dabs(f-f0) + dabs(xf-xfs) + dabs(yf-yfs) ) / f0

            !If we've got a suitable step then update our jacobian and continue
            if (err < errold) goto 89

            !if we did 2 reductions in the step size and we're still not at a lower
            !error then we'll try to re-calculate the derivative.  This seems to do a 
            !good job of improving convergence
            if (kreduce .ge. max_line_searches) then
                !if (err < 5*emax) goto 89
                !we want to make sure that we don't just keep on calling calc_deriv if things are going bad
                if (calcDerivInPrevIteration .eq. .false.) then
                    !if (numSearchFailures .ge. 2) then
                    calcDerivInPrevIteration = .true.
                    akmax = akmax / derivAKReduction
                    call calc_deriv(5.0d-2)
                    akmax = akmax * derivAKReduction
                    call saveState(fsave, xfsave, yfsave, hmsave, hx0save, hysave)	
                    goto 1234
                endif
                goto 89
            endif

            !write(*,*)'**reducing step.  Step reduction number:', kreduce
	        call reduceInc(sn1, sn2, sn3)
	        call stateInc(hmsave, hx0save, hysave, sn1, sn2, sn3)
            call gethx
	    enddo  !kreduce

89      continue

        calcDerivInPrevIteration = .false.

        !points of interest
	    call getPint

        !minimum height
	    call getMinH


        !Modified by Xinjiang Shen and Lion Huang.
        !Fixed the bug in slip model.
		y1=f-fsave 
		y2=xf-xfsave                  
		y3=yf-yfsave

		s1=hmin-hmsave
		s2=hx0-hx0save
		s3=hy-hysave
		
        fPrev(ktot) = f
        xfPrev(ktot) = xf
        yfPrev(ktot) = yf
        hminPrev(ktot) = hmin
        hx0Prev(ktot) = hx0
        hyPrev(ktot) = hy
        errorPrev(ktot) = err
		
		call jac_update(s1, s2, s3, y1, y2, y3, jac)

!       call jac_update(fsave, xfsave, yfsave, hmsave, hx0save, hysave)

	    call display(6)

90      call saveState(fsave, xfsave, yfsave, hmsave, hx0save, hysave)	

        ! check for convergence
        errold=err

        if (err.lt.etemp .and. ak.lt.akmax) then
            ! Success!!!
            if(iind.gt.0 .or. iadpt.eq.0) then
	            !call CoForce  !done, calculate the center of force
                call outputResult(numRes)
                write(*,*)
		        write(*,*)'STEADY STATE FLY HEIGHT FOUND.'
                write(*,*)
            endif
	        return
	    endif
	    
    enddo  !ktot

99  write(*,'(/)')

    if(ktot.ge.maxQNIterations) then
        write(*,*)
        write(*,*)'QUASI-NEWTON iteration terminated before'
        write(*,*)'residual criteria is satified.'
        write(*,*)
    else
        crash=.true.
        write(*,*)
        write(*,*)'Slider Crashed.'
        write(*,*)
    endif

!	write data indicatind error
!     BCox: We should write more values in here to REALLY notify the user that there was a problem

    err = -1.d0

    do i = 1, numPOI
	    hintNew(i) = -1.d0 / hm	
    enddo

    if (iind.gt.0 .or. iadpt.eq.0) call outputResult(numRes)
	
    return
    end
    
!===================================================
	subroutine add_lowest_errors_to_jacobian(fPrev, xfPrev, yfPrev, &
                                             hminPrev, hx0, hy, errorPrev)
!===================================================
    implicit none
    integer, parameter :: maxQNIterations = 50 !used to be 20 (8/18/05)
    real(8) :: fPrev(maxQNIterations), xfPrev(maxQNIterations), yfPrev(maxQNIterations)
    real(8) :: hminPrev(maxQNIterations), hx0(maxQNIterations), hy(maxQNIterations)
    real(8) :: errorPrev(maxQNIterations)    

    

	end subroutine add_lowest_errors_to_jacobian

!=====================================================
	subroutine jac_update(s1, s2, s3, y1, y2, y3, jac)
!=====================================================
!     BC -- For info on Xinjiang's fix, see CML Report 02008.pdf
!           Basically, if I remember right, we were overwriting an important global variable.  Not all that interesting
!           As for info on jac_update, we use Broyden's Method to update our Jacobian
!           so that we don't have to calculate our full Jacobian at each step
!           For more on Broyden's method, see the book "Numerical Methods for Unconstrained 
!           Optimization and Nonlinear Equations", Dennis, J.E and Schnabel, R.B.

!      Modified by Xinjiang shen and Lion Huang. This is the reason why we have
!      a bug in slip model in Quick419
!      include 'common.fi'
	implicit none
	
	real(8) :: s1, s2, s3, y1, y2, y3, jac(3,3), dj(3,3), denom, t1, t2, t3
	integer :: i, j

    denom=s1*s1+s2*s2+s3*s3

    t1=y1-(jac(1,1)*s1+jac(1,2)*s2+jac(1,3)*s3)
    t2=y2-(jac(2,1)*s1+jac(2,2)*s2+jac(2,3)*s3)
    t3=y3-(jac(3,1)*s1+jac(3,2)*s2+jac(3,3)*s3)

    dj(1,1)=t1*s1
    dj(1,2)=t1*s2
    dj(1,3)=t1*s3
    dj(2,1)=t2*s1
    dj(2,2)=t2*s2
    dj(2,3)=t2*s3
    dj(3,1)=t3*s1
    dj(3,2)=t3*s2
    dj(3,3)=t3*s3

    do i=1,3
	    do j=1,3
            jac(i,j) = jac(i,j) + dj(i,j)/denom
	    enddo
    enddo

    return
    end
    
!=====================================
	subroutine getInc(sn1, sn2, sn3)
!=====================================
    use Q4_globals
    implicit none

    real(8) :: sn1, sn2, sn3, jacsave(3,3)
    
    call equate(3, 3, 3, 3, jacsave, jac)
    call matrix33_inverse(jacsave)

    !Newton step
    sn1 = -jacsave(1,1)*(f-f0) - jacsave(1,2)*(xf-xfs) - jacsave(1,3)*(yf-yfs)
    sn2 = -jacsave(2,1)*(f-f0) - jacsave(2,2)*(xf-xfs) - jacsave(2,3)*(yf-yfs)
    sn3 = -jacsave(3,1)*(f-f0) - jacsave(3,2)*(xf-xfs) - jacsave(3,3)*(yf-yfs)

	return
	end
	
!=======================================================================
	subroutine saveState(fsave, xfsave, yfsave, hmsave, hx0save, hysave)
!=======================================================================
    use Q4_globals
    implicit none

    real(8) :: fsave, xfsave, yfsave, hmsave, hx0save, hysave
    
    fsave=f
    xfsave=xf
    yfsave=yf

    hmsave=hmin
    hx0save=hx0
    hysave=hy

	return
	end
	
!====================================
	subroutine getPint
!====================================
!     calculate height at points of interest

    use Q4_globals
    implicit none
    
    real(8) :: rinthl, zxint, zyint
    integer :: i

    do i=1,numPOI
	    call pointrecess(rinthl, xintNew(i), yintNew(i))
	    zxint = h0 + hx0*xintNew(i)
	    zyint = hy*(yintNew(i) - yl*0.5d0)
	    hintNew(i) = rinthl + zxint - zyint
    enddo

	return
	end
	
!==========================================
    subroutine getMinH
!==========================================
    use Q4_globals
    implicit none

    real(8) :: xv, yv, rcs, zxint, zyint, href, RailCenterX, RailCenterY
    integer :: i, j, k


    MinFH=1.d10
    MinFHLocX=.0d0
    MinFHLocY=0.d0

    do i=1,nx
	    do j=1,ny
	        xv=xref(i)
	        yv=yref(j)
	        if(i.eq.1)xv=1.d-8
	        if(i.eq.nx)xv=1.d0-1.d-8
	        if(j.eq.1)yv=1.d-8
	        if(j.eq.ny)yv=yl-1.d-8
	        call pointrecess(rcs,xv,yv)
	        zxint=h0+hx0*xref(i)
	        zyint=hy*(yref(j)-yl*0.5d0)
	        href=rcs+zxint-zyint

	        if(href .lt. MinFH)then
	            MinFH=href
	            MinFHLocX=xref(i)
	            MinFHLocY=yref(j)
	        endif
        enddo
    enddo
    
    !added 4/25/06 to check min rail heights and locations at POI
    !users sometime add a very small rail to model the transducer
    !the rail might be so small that the grid doesn't pick it up
    !we should also check points of interest
    do i = 1, numPOI
        if (hintNew(i) < MinFH) then
            MinFH = hintNew(i)
            MinFHLocX = xintNew(i)
            MinFHLocY = yintNew(i)
        endif
    end do
    
    do k = 1, numRails
        !use the bounding box to test the center point of the rail
        !we make the dangerous assumption that the transducer is a convex polygon
        RailCenterX = xpt2(k) + ((xpt1(k) - xpt2(k)) / 2)
        RailCenterY = ypt2(k) + ((ypt1(k) - ypt2(k)) / 2)
        call pointrecess(rcs, RailCenterX, RailCenterY)
        zxint = h0 + (hx0 * RailCenterX)
        zyint = hy * (RailCenterY - yl*0.5d0)
        href = rcs + zxint - zyint
        if(href .lt. MinFH)then
	        MinFH=href
	        MinFHLocX=RailCenterX
	        MinFHLocY=RailCenterY
	    endif
    enddo

    return
    end
	
!========================================
	subroutine reduceInc(sn1, sn2, sn3)
!========================================
    use Q4_globals
    implicit none
    
    real(8) :: sn1, sn2, sn3

	sn1 = sn1*0.5d0
	sn2 = sn2*0.5d0
	sn3 = sn3*0.5d0

	return
	end
	
!==========================================
	subroutine stateInc(hmsave, hx0save, hysave, sn1, sn2, sn3)
!==========================================

    use Q4_globals
    implicit none
    
    real(8) :: hmsave, hx0save, hysave, sn1, sn2, sn3

    hmin = hmsave + sn1
    hx0 = hx0save + sn2
    hy = hysave + sn3
    h0 = hmin - hx0

	return
	end