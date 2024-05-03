!=======================================================
      subroutine pointrecess(recess,xv,yv)
!=======================================================
!  recess at given point(xv,yv)
    use Q4_globals 
    use UserDefinedGeometry
    implicit none
    
    real(8), intent(INOUT) :: recess
    real(8), intent(IN) :: xv, yv
    
    real(8) :: recessold, xDistToCorner, yDistToCorner, fDistLower
    real(8) :: fDistUpper, fDist, distFromZero, udgheight
    real(8) :: crwn, cmbr, twst
    real(8) :: xw(6),yw(6)
    real(8) :: clambda, ew, ew1, ew2, wr1, wr2
    
    logical :: inside, iPointInFlushBoundary, iPointInWallioBoundary
    logical :: IsPointInQuad
    
    integer :: k, it, ip, ipt, iedge, iw
    
!----- Added 1/8/98 to correct point inclusion test bug -----
    real*8 xtemp(maxRailPts),ytemp(maxRailPts)
    integer tempint

!------------------------------------------------------------
      

    !assume base recess
    recess=rebase

    !point (xv,yv) belong to rail k?
    do k=numRails,1,-1
        !assume outside first
        inside=.false.

        !skip if outside range of this rail
        if(yv.lt.ypt1(k) .or. yv.gt.ypt2(k).or. &
           xv.lt.xpt1(k) .or. xv.gt.xpt2(k))  goto 10

!----- Added 1/8/98 to correct point inclusion test bug -----
!     copy rail points to temporary arrays

        do it=1,npoints(k)
            xtemp(it)=xrail(it,k)
            ytemp(it)=yrail(it,k)
        enddo

!       check whether the pt is inside the boundary.....
        call PNPOLY (xv, yv, xtemp, ytemp, npoints(k), tempint) 
        if (tempint .eq. -1) then
            inside=.false.
            call PTONPOLYBOUNDARY(xv, yv, xtemp, ytemp, npoints(k), tempint)
            if (tempint .eq. 1) then
                inside=.true.
            endif
        else 
            inside=.true.
        endif
!------------------------------------------------------------

	    if (inside) then
	        if(istep(k).ne.0) then
	            recess=hramp(1,k)
	        else
                recess=xv*cramp(1,k)+yv*cramp(2,k)+cramp(3,k)
	        endif
	        goto 20
	    end if   ! if inside
10	    continue
	enddo ! numRails

20  continue

!     check for sloped walls

    recessold = recess
!	point (xv,yv) belong to sloped walls of	 rail k?
	do k=numRails,1,-1
!	    skip if not a raised step
	    if(istep(k).eq.0.or.hramp(1,k).ge.rebase) goto 18
	    do ip=1,npoints(k)
!	        skip if outside range of this rail
	        if(yv.lt.yw1(ip,k).or.yv.gt.yw2(ip,k).or. &
     	        xv.lt.xw1(ip,k).or.xv.gt.xw2(ip,k).or. &
     	        indexw(k,ip).eq.0)  goto 15
	
!...........................................
!     See if we're inside the flush boundary	    
		    xw(1) = xwalli(ip, k)
		    xw(2) = xwalli(ip+1, k)
		    xw(3) = flush_upper_xwallo(ip, k)
		    xw(4) = flush_lower_xwallo(ip, k)
		    xw(5) = xw(1)
		    xw(6) = xw(2)
		    yw(1) = ywalli(ip, k)
		    yw(2) = ywalli(ip+1, k)
		    yw(3) = flush_upper_ywallo(ip, k)
		    yw(4) = flush_lower_ywallo(ip, k)
		    yw(5) = yw(1)
		    yw(6) = yw(2)

		    iPointInFlushBoundary = .false.
		    if(indexw(k,ip) .ne. 0) then
		  	    iPointInFlushBoundary = IsPointInQuad(xv, yv, xw, yw)
		    endif
!...........................................
!           Now see if we're in the total wallo boundary
	        xw(1)=xwalli(ip,k)
	        xw(2)=xwalli(ip+1,k)
	        xw(3)=xwallo(ip+1,k)
	        xw(4)=xwallo(ip,k)
	        xw(5)=xw(1)
	        xw(6)=xw(2)
	        yw(1)=ywalli(ip,k)
	        yw(2)=ywalli(ip+1,k)
	        yw(3)=ywallo(ip+1,k)
	        yw(4)=ywallo(ip,k)
	        yw(5)=yw(1)
	        yw(6)=yw(2)

	        iPointInWallioBoundary = .false.
	        if(indexw(k,ip) .ne. 0) then
	            iPointInWallioBoundary = IsPointInQuad(xv, yv, xw, yw)
	        endif
        
!           Now see what type of corner we're at...
		    if ((iPointInFlushBoundary .eq. .true.) .and. &
                (iPointInWallioBoundary .eq. .false. )) then
!		        convex corner.  We're actually not on a wall
                goto 15 !next point on rail
            else if (iPointInFlushBoundary .eq. .true.) then
!               regular wall boundary calculation
                iedge=ip
                clambda=((yw(2)-yw(1))*(xv-xw(1))-(xw(2)-xw(1))*(yv-yw(1)))/((yw(2)-yw(1))*(xw(3)-xw(1)) &
     		            -(xw(2)-xw(1))*(yw(3)-yw(1)))
		        iw=indexw(k,ip)
	            ew=clambda*wpoint(iw,nwpoint(iw))+(1.d0-clambda) * wpoint(iw,1)
	            do ipt=1,nwpoint(iw)-1
		            ew1=wpoint(iw,ipt)
		            ew2=wpoint(iw,ipt+1)
		            if(ew.ge.ew1.and.ew.lt.ew2)then
		                wr1=wrecess(iw,ipt)
		                wr2=wrecess(iw,ipt+1)
		                recess=wr1+(ew-ew1)/(ew2-ew1)*(wr2-wr1)
		                goto 25
		            endif
	            enddo
25	            continue
                !we can have negative wall distances.  Check to see if our rail point is inside the rail poly or
                !outside.  Also include tolerance for floating point errors
	            if(ew .gt. 1.0d-12)recess=dmin1(recess,recessold)
	            recessold=recess
!               ...............................................
            else if (iPointInWallioBoundary .eq. .true.) then
!               New wall calculation for corner region	   
                recess = rebase

                xDistToCorner = xv - xwalli(ip, k) !corner
                yDistToCorner = yv - ywalli(ip, k)
                fDistLower = sqrt(xDistToCorner*xDistToCorner + yDistToCorner*yDistToCorner)
                xDistToCorner = xv - xwalli(ip+1, k)!corner
                yDistToCorner = yv - ywalli(ip+1, k)
                fDistUpper = sqrt(xDistToCorner*xDistToCorner + yDistToCorner*yDistToCorner)
                fDist = -1.0
                if (fDistUpper < fDistLower) then 
                    fDist = fDistUpper 
                else
                    fDist = fDistLower
                endif
                
!               We'll use the existing trick to determine if we should allow this point
!               to be below the top rail. (related to walls with neg nominal distances 
!               see bug report in PointRecess.h)            	      
	            clambda=((yw(2)-yw(1))*(xv-xw(1))-(xw(2)-xw(1))*(yv-yw(1)))/((yw(2)-yw(1))*(xw(3)-xw(1)) &
                        -(xw(2)-xw(1))*(yw(3)-yw(1)))
!               BC: notice how k changes sides...there is a god, he has a PhD, and he hates me...
		        iw=indexw(k,ip)
	            ew=clambda*wpoint(iw,nwpoint(iw))+(1.d0-clambda) * wpoint(iw,1)
            	
                distFromZero = wpoint(iw,1)
  	            do ipt=1,nwpoint(iw)-1
		            ew1=wpoint(iw,ipt) - distFromZero
		            ew2=wpoint(iw,ipt+1) - distFromZero
		            if(fDist.ge.ew1.and.fDist.lt.ew2)then
		                wr1=wrecess(iw,ipt)
		                wr2=wrecess(iw,ipt+1)
		                recess=wr1+(fDist-ew1)/(ew2-ew1)*(wr2-wr1)
		                goto 26
		            endif
	            enddo
26              continue
	            if(ew .gt. 1.0d-12)recess=dmin1(recess,recessold)
	            recessold=recess
	        endif
15	        continue
	    enddo	  !npoints
18	    continue
    enddo	!numRails


    if(xv.lt.xt)recess=dmax1(ht*(1.d0-xv/xt),recess)

    crwn=4.d0*crown*xv*(1.d0-xv)
    cmbr=4.d0*camber*yv*(yl-yv)/yl**2
    twst=4.d0*twist*(xv-0.5d0)*(yv-.5d0*yl)/yl
    recess=recess-crwn-cmbr-twst
    
    !usergeom changed 2/24/05 (always add to camber crown and twist)
    if (numUDGRegions .gt. 0) then
        !udgamount = udgheight(xv,yv)
        recess = recess - (udgheight(xv,yv) * (0.000000001 / hm))
    endif

    return
    end subroutine pointrecess


!=============================================================
    subroutine wallprofile
!=============================================================
    use Q4_globals
    implicit none
    
    integer :: i, im1, ip1, k, iymax, iymin, im, np, iwm, iw
    real(8) :: clock, x1, x2, x3, y1, y2, y3, ymin, ymax
    real(8) :: di1, do1, di2, do2, x_mag, y_mag, vec_mag, dnom_dist, wall_percent, tempx

    do k=1, numRails
	    xrail(npoints(k)+1,k) = xrail(1,k)
	    yrail(npoints(k)+1,k) = yrail(1,k)
	    xrail(npoints(k)+2,k) = xrail(2,k)
	    yrail(npoints(k)+2,k) = yrail(2,k)

!c       BC - 9/28/04
!c       Bug: We don't put a wall profile on any rail that is lower than the base recess
!c       This is a problem since people use contamination chanels and other etches below the base recess
!c       Keep in mind that most walls below .3 recess don't really matter but it's possible that we
!c       have a wall stretching from 0 recess to below the base recess.
	    if(istep(k).ne.0.and.hramp(1,k).lt.rebase)then
	        ymax=-1.d0
	        ymin=1.d0

            !obtain the highest or lowest point
	        do i=1, npoints(k)
	            if(yrail(i,k).gt.ymax) then
		            ymax=yrail(i,k)
		            iymax=i
	            endif
	            if(yrail(i,k).lt.ymin) then
		            ymin=yrail(i,k)
		            iymin=i
	            endif
	        enddo

	        if(iymax.eq.1)then
	            im=iymin
	        else
	            im=iymax
	        endif
            ! determine if rail polygon is clockwise or counterclockwise
            if((xrail(im,k)-xrail(im-1,k))*(yrail(im+1,k)-  &
     	            yrail(im,k)).gt. (xrail(im+1,k)-xrail(im,k))*  &
     	            (yrail(im,k)-yrail(im-1,k))) then
	            clock=-1.d0
            else
	            clock=1.d0
            endif

	        np=npoints(k)
	        indexw(k,np+1)=indexw(k,1)
            ! since subscript starts from 1, we have to shift the loop
            do i=2,np+1
	            im1=i-1
	            ip1=i+1
	            x1=xrail(im1,k)
	            x2=xrail(i,k)
	            x3=xrail(ip1,k)
	            y1=yrail(im1,k)
	            y2=yrail(i,k)
	            y3=yrail(ip1,k)
	            iwm=indexw(k,im1)
	            iw=indexw(k,i)
                !if no profile, set outer distance to 0
	            if(iwm.eq.0)then
		            di1=0.d0
		            do1=0.d0
	            else
		            di1=wpoint(iwm,1)
		            do1=wpoint(iwm,nwpoint(iwm))
	            endif
                ! same deal for the upper point
	            if(iw.eq.0)then
		            di2=0.d0
		            do2=0.d0
	            else
		            di2=wpoint(iw,1)
		            do2=wpoint(iw,nwpoint(iw))
	            endif

                !obtain intersection of outside wall lines
		        call intersection(clock,x1,y1,x2,y2,x3,y3,di1,do1,di2,do2, &
		                            xwalli(i,k),ywalli(i,k), xwallo(i,k),ywallo(i,k))
	        enddo

            
	        xwalli(1,k)=xwalli(np+1,k)
	        xwallo(1,k)=xwallo(np+1,k)
	        ywalli(1,k)=ywalli(np+1,k)
	        ywallo(1,k)=ywallo(np+1,k)
	        
!........................................................
!     new wall generation algorithm added 11/04
!     see CPointRecess.h in interface for explanation...
!     in short, we round the outer perimeter of the 
!     walls at convex corners
            do i=1, npoints(k)
                ip1=i+1
                iw=indexw(k, i)
                if ((iw .le. 0) .or. ((iw-1) .gt. numWalls)) then
                    flush_lower_xwallo(i,k) = xwalli(i, k)
	                flush_upper_xwallo(i,k) = xwalli(ip1, k)
	                flush_lower_ywallo(i,k) = ywalli(i, k)
	                flush_upper_ywallo(i,k) = ywalli(ip1, k)
                else
                    x_mag = xwalli(ip1, k) - xwalli(i,k)
	                y_mag = ywalli(ip1, k) - ywalli(i,k)
                    ! get length of our line
	                vec_mag = sqrt(x_mag*x_mag + y_mag*y_mag)
                    ! get nominal distance of wall (that's how far out we'll end up going
	                dnom_dist = (wpoint(iw, nwpoint(iw))) - (wpoint(iw, 1))
                    ! Make sure we have a nominal distance and don't divide by 0
	                wall_percent = 0.0
	                if (vec_mag .gt. 0.0) wall_percent = dnom_dist / vec_mag
                    !rotate 90 degrees
                    tempx = x_mag
	                x_mag = -y_mag
	                y_mag = tempx
                    ! shorten/lengthen our vector so it's the nominal wall length
	                x_mag = x_mag * wall_percent * clock
	                y_mag = y_mag * wall_percent * clock
                    ! translate xwalli line out to x_mag
	                flush_lower_xwallo(i,k) = x_mag + xwalli(i,k)
	                flush_upper_xwallo(i,k) = x_mag + xwalli(ip1, k)
	                flush_lower_ywallo(i,k) = y_mag + ywalli(i,k)
	                flush_upper_ywallo(i,k) = y_mag + ywalli(ip1, k)
                endif
            enddo
!.........................................................end 11/04          
        endif
    enddo

    do k=1,numRails
        do i=1,npoints(k)
	        xw1(i,k)=dmin1(xwalli(i,k),xwalli(i+1,k),xwallo(i,k),xwallo(i+1,k))
	        xw2(i,k)=dmax1(xwalli(i,k),xwalli(i+1,k),xwallo(i,k),xwallo(i+1,k))
	        yw1(i,k)=dmin1(ywalli(i,k),ywalli(i+1,k),ywallo(i,k),ywallo(i+1,k))
	        yw2(i,k)=dmax1(ywalli(i,k),ywalli(i+1,k),ywallo(i,k),ywallo(i+1,k))
    	enddo
    enddo

    return
    end subroutine wallprofile
    
          
!=========================================================
    integer function IsPointInQuad(xv, yv, xw, yw)
!=========================================================

    use Q4_globals
    implicit none
    
    REAL(8) :: xv, yv, xw(6),yw(6)
    real(8) :: xcross
    integer :: i, ncross
    logical :: inside
    
    inside=.false.
    ncross=0
    
    do i=1,4
        if(dabs(yw(i+1)-yw(i)).gt.1.d-10)then
	        if(((yw(i).gt.yv).and.(yw(i+1).lt.yv)) .or.               &
     	         ((yw(i).lt.yv).and.(yw(i+1).gt.yv)))then
		        xcross=xw(i)+(xw(i+1)-xw(i))/(yw(i+1)-yw(i)) *(yv-yw(i))
		        if(xv.gt.xcross) ncross=ncross+1
		        if(ncross.gt.1)goto 53
		    else if(yw(i+1).eq.yv)then
		        if(((yw(i+1)-yw(i+2))*(yw(i+1)-yw(i)).lt.0.d0).and.    &
     	            (xv.gt.xw(i+1)))ncross=ncross+1
		        if(ncross.gt.1)goto 53
		    endif
	    else
		    if(((xv.ge.xw(i).and.xv.le.xw(i+1)).or.                 &
     	        (xv.le.xw(i).and.xv.ge.xw(i+1))).and.              &
     	         dabs(yv-yw(i)).lt.1.d-10)then
		        inside=.true.
		        goto 52
		    endif
	    endif
	enddo
52  continue
	if(ncross.eq.1) inside=.true.
53      continue
	IsPointInQuad = inside

	end
      
!=======================================================
    subroutine PNPOLY(PX, PY, XX, YY, N, INOUT) 
!=======================================================

!c  Note from BCox 2/17/05: This subroutine actually doesn't check to see if we're
!c  on the boundary of ALL points.  This is bad since we snap the grid
!c  to our boundary!  From Prof. Franklin himself:
!c  "If you want to know when a point is exactly on the boundary, 
!c   you need another program"
!c   See PTONPOLYBOUNDARY below

!c  Courtesy of Prof. Randolph Franklin, Rensselaer Polytechnic Inst, Troy NY
!C     ..................................................................
!C                                                                       
!C        SUBROUTINE PNPOLY                                              
!C                                                                       
!C        PURPOSE                                                        
!C           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
!C                                                                       
!C        USAGE                                                          
!C           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
!C                                                                       
!C        DESCRIPTION OF THE PARAMETERS                                  
!C           PX      - X-COORDINATE OF POINT IN QUESTION.                
!C           PY      - Y-COORDINATE OF POINT IN QUESTION.                
!C           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
!C                     VERTICES OF POLYGON.                              
!C           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
!C                     VERTICES OF POLYGON.                              
!C           N       - NUMBER OF VERTICES IN THE POLYGON.                
!C           INOUT   - THE SIGNAL RETURNED:                              
!C                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
!C                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
!C                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
!C                                                                       
!C        REMARKS                                                        
!C           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
!C           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
!C           OPTIONALLY BE INCREASED BY 1.                               
!C           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
!C           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
!C           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
!C           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
!C           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
!C           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
!C           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
!C                                                                       
!C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
!C           NONE                                                        
!C                                                                       
!C        METHOD                                                         
!C           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
!C           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
!C           POINT IS INSIDE OF THE POLYGON.                             
!C     ..................................................................
!C                                                                       
    use Q4_globals
    implicit none

    integer :: N
    REAL(8) :: PX, PY, X(maxRailPts),Y(maxRailPts),XX(N),YY(N)
    LOGICAL MX,MY,KX,KY                                              
    INTEGER :: O
    integer i, J, MAXDIM, INOUT
    
!   OUTPUT UNIT FOR PRINTED MESSAGES                                 
    DATA O/6/                                                         

    MAXDIM = maxRailPts+10
    IF (N .LE. MAXDIM) GO TO 6                                            
!    WRITE(O,7)                                                        
!7   FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY. 1 RESULTS INVALID')                                                 
!    RETURN
6   DO 1 I=1,N
        X(I)=XX(I)-PX
1       Y(I)=YY(I)-PY
    INOUT=-1
    DO 2 I=1,N
        J=1+MOD(I,N)
        MX=X(I).GE.0.0
        KX=X(J).GE.0.0
        MY=Y(I).GE.0.0
        KY=Y(J).GE.0.0
        IF(.NOT.((MY.OR.KY).AND.(MX.OR.KX)).OR.(MX.AND.KX)) GO TO 2
        IF(.NOT.(MY.AND.KY.AND.(MX.OR.KX).AND..NOT.(MX.AND.KX))) GO TO 3
        INOUT=-INOUT
        GO TO 2
3       IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5
4       INOUT=0
        RETURN
5       INOUT=-INOUT
2   CONTINUE

    RETURN
    END
      
      
!====================================================================      
    subroutine PTONPOLYBOUNDARY(PX, PY, XX, YY, N, INOUT) 
!====================================================================      
!c     Created 2/17/05 to correct a bug in the (long used) PNPOLY subroutine
!c     PNPOLY will not catch when we're on certain boundaries of our polygon
!c     Usually this wouldn't be a problem but since we snap the grid to rail
!c     edges that are parallel to the X and Y axes, we often end up on the
!c     edge of a polygon.  If we miss it (which happens) we end up with an
!c     entire segment of grids that are at the base recess.  This is bad.
!
!C        DESCRIPTION OF THE PARAMETERS                                  
!C           PX      - X-COORDINATE OF POINT IN QUESTION.                
!C           PY      - Y-COORDINATE OF POINT IN QUESTION.                
!C           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
!C                     VERTICES OF POLYGON.                              
!C           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
!C                     VERTICES OF POLYGON.                              
!C           N       - NUMBER OF VERTICES IN THE POLYGON.                
!C           INOUT   - THE SIGNAL RETURNED:                              
!C                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
!C                      1 IF THE POINT IS ON AN EDGE OR AT A VERTEX,

    implicit none
    
    integer :: N, MAXDIM
    REAL(8) :: PX, PY  
    REAL(8) :: XX(N),YY(N),YMIN,YMAX,XMIN,XMAX
    INTEGER :: INOUT, I, J

    INOUT=-1
    DO I=1,N
        J=1+MOD(I,N)

        !First check if we're on a corner
        if ((PX .eq. XX(i)) .and. (PY .eq. YY(i))) then
            INOUT = 1
            return
        endif

!     now check if we're sitting on a line in the X or Y direction only
        if ((XX(i) - XX(j) .eq. 0.d0) .and. (PX .eq. XX(i))) then
!           we're inline with a vertical line so see if we're between the endpoints
            if (YY(i) .lt. YY(j)) then 
                YMIN = YY(i)
                YMAX = YY(j)
            else
                YMIN = YY(j)
                YMAX = YY(i)
            endif
            if ((PY .gt. YMIN) .and. (PY .lt. YMAX)) then
                INOUT = 1
                return
            endif
        endif
        if ((YY(i) - YY(j) .eq. 0.d0) .and. (PY .eq. YY(i)))then
!           we're inline with a horizontal line, see if we're between the endpoints        
            if (XX(i) .lt. XX(j)) then 
                XMIN = XX(i)
                XMAX = XX(j)
            else
                XMIN = XX(j)
                XMAX = XX(i)
            endif
            if ((PX .gt. XMIN) .and. (PX .lt. XMAX)) then
                INOUT = 1
                return
            endif
        endif
    enddo
!   Not on a vertical or horizontal line so (for now) we'll assume that 
!   we're not on the boundary (DANGEROUS, PLEASE FIX!!!)
    RETURN                                                            
    END      
      

!==========================================================
      subroutine intersection(clock,x1,y1,x2,y2,x3,y3,di1,do1,di2,do2,xi,yi,xo,yo)
!==========================================================

      implicit none
      
      real(8) :: clock,x1,y1,x2,y2,x3,y3,di1,do1,di2,do2,xi,yi,xo,yo
      real(8) :: a1, a2, b1, b2, c1, c2, coe1, coe2, c1i, c2i, c1o, c2o
      real(8) :: determ0, determxi, determyi, determxo, determyo


      a1=y2-y1
      a2=y3-y2
      b1=x1-x2
      b2=x2-x3
      c1=x1*y2-x2*y1
      c2=x2*y3-x3*y2
      coe1=sqrt(a1*a1+b1*b1)
      coe2=sqrt(a2*a2+b2*b2)

      c1i=c1-coe1*di1*clock
      c2i=c2-coe2*di2*clock
      c1o=c1-coe1*do1*clock
      c2o=c2-coe2*do2*clock

      determ0=a1*b2-a2*b1

      if(dabs(determ0) .lt. 1.e-30) then
	    write(*, *)'Divide by zero in Subroutine intersection!'
	    write(*, *)'Colinear points are not allowed in the rail!'
      endif

      determxi=c1i*b2-c2i*b1
      determyi=a1*c2i-a2*c1i
      determxo=c1o*b2-c2o*b1
      determyo=a1*c2o-a2*c1o

      xi=determxi/determ0
      yi=determyi/determ0
      xo=determxo/determ0
      yo=determyo/determ0

      return
      end subroutine intersection