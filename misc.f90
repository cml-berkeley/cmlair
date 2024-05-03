! Restriction and Prolongation
! ----------------------------
! Residual - fine grid residual is linearly distributed to the coarse grid
! node enclosing it.  Use the Finite Element trial functions (weight of barycentric coords)
!
! Prolonging the Correction - linear interpolate the coarse grid at the
! location of the fine grid node

! Pressure restriction- linear interpolation over the fine grid node
! coarse grid node needs to know which fine grid square encloses it

!Thus, each fine grid node needs to know which coarse grid node encloses it
!and each coarse grid nodes need to know which fine grid node encloses it.
!Nesting Info holds this info in enclosingFineRectXY and enclosingCoarseRectXY


!=================================================================
    subroutine restrictPressure(fineMesh, coarseMesh, xrefFine, yrefFine, xrefCoarse, yrefCoarse, &
                                nxFine, nyFine, nxCoarse, nyCoarse, enclosingFineRectX, enclosingFineRectY)
!=================================================================
    implicit none
        
    integer :: nxFine, nxCoarse, nyFine, nyCoarse
    real(8) :: fineMesh(nxFine, nyFine), coarseMesh(nxCoarse, nyCoarse), &
               xrefFine(nxFine), yrefFine(nyFine), xrefCoarse(nxCoarse), yrefCoarse(nyCoarse)
    real(8) :: enclosingFineRectX(nxCoarse), enclosingFineRectY(nyCoarse)
    
    integer :: xIndex, yIndex, i, j
    real(8) :: length, height, b, c, xx, yy, H1, H2, H3, H4
    
    do i = 1, nxCoarse
        do j = 1, nyCoarse
            xIndex = enclosingFineRectX(i);
            yIndex = enclosingFineRectY(j);

            length = xrefFine(xIndex+1) - xrefFine(xIndex);
            height = yrefFine(yIndex+1) - yrefFine(yIndex);
            b = 0.5*length;
            c = 0.5*height;

            !local/barycentric coords of our point on the coarse grid
            xx = xrefCoarse(i) - xrefFine(xIndex) - b;
            yy = yrefCoarse(j) - yrefFine(yIndex) - c;

            H1 = (1/(4*b*c)) * (b-xx) * (c-yy);
            H2 = (1/(4*b*c)) * (b+xx) * (c-yy);
            H3 = (1/(4*b*c)) * (b+xx) * (c+yy);
            H4 = (1/(4*b*c)) * (b-xx) * (c+yy);

            coarseMesh(i, j) =  H1 * fineMesh(xIndex, yIndex) + H2 * fineMesh(xIndex+1, yIndex) + &
                                H3 * fineMesh(xIndex+1, yIndex+1) + H4 * fineMesh(xIndex, yIndex+1);    
        enddo
    enddo
    
    return
    end
    
    
!=================================================================
    subroutine restrictToClosest(fineMesh, coarseMesh, xrefFine, yrefFine, xrefCoarse, yrefCoarse, &
                                nxFine, nyFine, nxCoarse, nyCoarse, enclosingFineRectX, enclosingFineRectY)
!=================================================================
    implicit none
        
    integer :: nxFine, nxCoarse, nyFine, nyCoarse
    real(8) :: fineMesh(nxFine, nyFine), coarseMesh(nxCoarse, nyCoarse), &
               xrefFine(nxFine), yrefFine(nyFine), xrefCoarse(nxCoarse), yrefCoarse(nyCoarse)
    real(8) :: enclosingFineRectX(nxCoarse), enclosingFineRectY(nyCoarse)
    real(8) :: distanceUpperX, distanceLowerX, distanceUpperY, distanceLowerY
    real(8) :: maxTemp1, maxTemp2
    
    integer :: xIndex, yIndex, i, j
    real(8) :: length, height, b, c, xx, yy, H1, H2, H3, H4
    
    do i = 1, nxCoarse
        do j = 1, nyCoarse
            xIndex = enclosingFineRectX(i);
            yIndex = enclosingFineRectY(j);
            
!            distanceUpperX = abs(xrefFine(xIndex+1) - xrefCoarse(i))  !abs is here just in case
!            distanceLowerX = abs(xrefCoarse(i) - xrefFine(xIndex))
!
!            distanceUpperY = abs(yrefFine(yIndex+1) - yrefCoarse(j))
!            distanceLowerY = abs(yrefCoarse(j) - yrefFine(yIndex))
!
!            if (distanceUpperX < distanceLowerX) xIndex = xIndex+1
!            if (distanceUpperY < distanceLowerY) yIndex = yIndex+1
            maxTemp1 = max( abs(fineMesh(xIndex,yIndex)), abs(fineMesh(xIndex+1,yIndex)) )
            maxTemp2 = max( abs(fineMesh(xIndex, yIndex+1)), abs(fineMesh(xIndex+1, yIndex+1)))
            coarseMesh(i,j) = max( maxTemp1, maxTemp2 )
            
        
!            xIndex = enclosingFineRectX(i);
!            yIndex = enclosingFineRectY(j);
!
!            length = xrefFine(xIndex+1) - xrefFine(xIndex);
!            height = yrefFine(yIndex+1) - yrefFine(yIndex);
!            b = 0.5*length;
!            c = 0.5*height;
!
!            !local/barycentric coords of our point on the coarse grid
!            xx = xrefCoarse(i) - xrefFine(xIndex) - b;
!            yy = yrefCoarse(j) - yrefFine(yIndex) - c;
!
!            H1 = (1/(4*b*c)) * (b-xx) * (c-yy);
!            H2 = (1/(4*b*c)) * (b+xx) * (c-yy);
!            H3 = (1/(4*b*c)) * (b+xx) * (c+yy);
!            H4 = (1/(4*b*c)) * (b-xx) * (c+yy);
!
!            coarseMesh(i, j) =  H1 * fineMesh(xIndex, yIndex) + H2 * fineMesh(xIndex+1, yIndex) + &
!                                H3 * fineMesh(xIndex+1, yIndex+1) + H4 * fineMesh(xIndex, yIndex+1);    
        enddo
    enddo
    
    return
    end subroutine restrictToClosest

!=================================================================
    subroutine restrictResidual(fineMesh, coarseMesh, xrefFine, yrefFine, xrefCoarse, yrefCoarse, &
                                nxFine, nyFine, nxCoarse, nyCoarse, enclosingCoarseRectX, enclosingCoarseRectY)
!=================================================================
    implicit none
    
    
    integer :: nxFine, nxCoarse, nyFine, nyCoarse
    real(8) :: fineMesh(nxFine, nyFine), coarseMesh(nxCoarse, nyCoarse), &
               xrefFine(nxFine), yrefFine(nyFine), xrefCoarse(nxCoarse), yrefCoarse(nyCoarse)
    real(8) :: enclosingCoarseRectX(nxFine), enclosingCoarseRectY(nyFine)
    
    integer :: xIndex, yIndex, i, j
    real(8) :: length, height, b, c, xx, yy, H1, H2, H3, H4
    
    do i = 1, nxCoarse
        do j = 1, nyCoarse
            coarseMesh = 0.0
        enddo
    enddo

    do i = 2, nxFine-1
        do j = 2, nyFine-1
            xIndex = enclosingCoarseRectX(i);
            yIndex = enclosingCoarseRectY(j);

            length = xrefCoarse(xIndex+1) - xrefCoarse(xIndex);
            height = yrefCoarse(yIndex+1) - yrefCoarse(yIndex);
            b = 0.5*length;
            c = 0.5*height;

            !local/barycentric coords of our point on the fine grid
            xx = xrefFine(i) - xrefCoarse(xIndex) - b;
            yy = yrefFine(j) - yrefCoarse(yIndex) - c;

            H1 = (1/(4*b*c)) * (b-xx) * (c-yy);
            H2 = (1/(4*b*c)) * (b+xx) * (c-yy);
            H3 = (1/(4*b*c)) * (b+xx) * (c+yy);
            H4 = (1/(4*b*c)) * (b-xx) * (c+yy);

            coarseMesh(xIndex, yIndex) = coarseMesh(xIndex, yIndex) +  H1 * fineMesh(i, j);
            coarseMesh(xIndex+1, yIndex) = coarseMesh(xIndex+1, yIndex) +  H2 * fineMesh(i, j);
            coarseMesh(xIndex+1, yIndex+1) = coarseMesh(xIndex+1, yIndex+1) +  H3 * fineMesh(i, j);
            coarseMesh(xIndex, yIndex+1) = coarseMesh(xIndex, yIndex+1) +  H4 * fineMesh(i, j);
        enddo
    enddo
    
    return
    end subroutine restrictResidual  

!=================================================================
    subroutine prolongCorrection(fineMesh, coarseMesh, xrefFine, yrefFine, xrefCoarse, yrefCoarse, &
                                 nxFine, nyFine, nxCoarse, nyCoarse, enclosingCoarseRectX, enclosingCoarseRectY)
!=================================================================    
    implicit none
    

    integer :: nxFine, nxCoarse, nyFine, nyCoarse
    real(8) :: fineMesh(nxFine, nyFine), coarseMesh(nxCoarse, nyCoarse), &
               xrefFine(nxFine), yrefFine(nyFine), xrefCoarse(nxCoarse), yrefCoarse(nyCoarse)
    real(8) :: enclosingCoarseRectX(nxFine), enclosingCoarseRectY(nyFine)
    
    integer :: xIndex, yIndex, i, j
    real(8) :: length, height, b, c, xx, yy, H1, H2, H3, H4

    do i = 2, nxFine-1
        do j = 2, nyFine-1
            xIndex = enclosingCoarseRectX(i)
            yIndex = enclosingCoarseRectY(j)

            length = xrefCoarse(xIndex+1) - xrefCoarse(xIndex)
            height = yrefCoarse(yIndex+1) - yrefCoarse(yIndex)
            b = 0.5*length
            c = 0.5*height

            !local/barycentric coords of our point on the coarse grid
            xx = xrefFine(i) - xrefCoarse(xIndex) - b
            yy = yrefFine(j) - yrefCoarse(yIndex) - c

            H1 = (1/(4*b*c)) * (b-xx) * (c-yy)
            H2 = (1/(4*b*c)) * (b+xx) * (c-yy)
            H3 = (1/(4*b*c)) * (b+xx) * (c+yy)
            H4 = (1/(4*b*c)) * (b-xx) * (c+yy)

            fineMesh(i, j) = H1 * coarseMesh(xIndex, yIndex) + H2 * coarseMesh(xIndex+1, yIndex) + &
                             H3 * coarseMesh(xIndex+1, yIndex+1) + H4 * coarseMesh(xIndex, yIndex+1)
        enddo
    enddo

    return
    end subroutine prolongCorrection
    
    
!=================================================================    
    subroutine SetUpMGNesting()
!=================================================================
    use Q4_globals
    use NestingInfo
    implicit none
    
    !go through each grid and set up the MG nesting
    call findEnclosingCoarseRects(xref, yref, xref1, yref1, nx, ny, nx1, ny1, enclosingCoarseRectX, enclosingCoarseRectY)
    call findEnclosingCoarseRects(xref1, yref1, xref2, yref2, nx1, ny1, nx2, ny2, enclosingCoarseRectX1, enclosingCoarseRectY1)
    call findEnclosingCoarseRects(xref2, yref2, xref3, yref3, nx2, ny2, nx3, ny3, enclosingCoarseRectX2, enclosingCoarseRectY2)
    call findEnclosingCoarseRects(xref3, yref3, xref4, yref4, nx3, ny3, nx4, ny4, enclosingCoarseRectX3, enclosingCoarseRectY3)
    
    call findEnclosingFineRects(xref, yref, xref1, yref1, nx, ny, nx1, ny1, enclosingFineRectX1, enclosingFineRectY1)
    call findEnclosingFineRects(xref1, yref1, xref2, yref2, nx1, ny1, nx2, ny2, enclosingFineRectX2, enclosingFineRectY2)
    call findEnclosingFineRects(xref2, yref2, xref3, yref3, nx2, ny2, nx3, ny3, enclosingFineRectX3, enclosingFineRectY3)
    call findEnclosingFineRects(xref3, yref3, xref4, yref4, nx3, ny3, nx4, ny4, enclosingFineRectX4, enclosingFineRectY4)
    
    end subroutine SetUpMGNesting
    
    
!=================================================================    
    subroutine findEnclosingCoarseRects(xrefFine, yrefFine, xrefCoarse, yrefCoarse, &
                                        nxFine, nyFine, nxCoarse, nyCoarse, &
                                        enclosingCoarseRectX, enclosingCoarseRectY)
!=================================================================      
    implicit none
    
    integer :: nxFine, nxCoarse, nyFine, nyCoarse
    real(8) :: xrefFine(nxFine), yrefFine(nyFine), xrefCoarse(nxCoarse), yrefCoarse(nyCoarse)
    real(8) :: enclosingCoarseRectX(nxFine), enclosingCoarseRectY(nyFine)
    
    integer xIndex, yIndex, i
    
    
    !for each node in the fine grid, find enclosing coarse grid square 
    do i = 1, nxFine
        xIndex = 1;
        do while (xrefCoarse(xIndex) .le. xrefFine(i) .and. xIndex .lt. nxCoarse)
            xIndex = xIndex+1;
        end do
        !since we used <= we need to subtract 1 from xIndex
        !our index will hold the lower index of the enclosing square
        xIndex = xIndex-1;
        
        enclosingCoarseRectX(i) = xIndex;
    enddo

    !now in the Y direction
    do i = 1, nyFine
        yIndex = 1;
        do while (yrefCoarse(yIndex) .le. yrefFine(i) .and. yIndex .lt. nyCoarse)
            yIndex = yIndex+1;
        end do
        yIndex = yIndex-1;
        
        enclosingCoarseRectY(i) = yIndex;
    enddo
    
        
    return
    end subroutine findEnclosingCoarseRects
    
    
!=================================================================    
    subroutine findEnclosingFineRects(xrefFine, yrefFine, xrefCoarse, yrefCoarse, &
                                        nxFine, nyFine, nxCoarse, nyCoarse, &
                                        enclosingFineRectX, enclosingFineRectY)
!=================================================================      
    implicit none
    
    integer :: nxFine, nxCoarse, nyFine, nyCoarse
    real(8) :: xrefFine(nxFine), yrefFine(nyFine), xrefCoarse(nxCoarse), yrefCoarse(nyCoarse)
    real(8) :: enclosingFineRectX(nxCoarse), enclosingFineRectY(nyCoarse)
    integer xIndex, yIndex, i

    !for each node in the coarse grid, find the rect on the fine grid that encloses it
    do i = 1, nxCoarse
        !find enclosing fine grid square
        xIndex = 1;
        do while (xrefFine(xIndex) .le. xrefCoarse(i) .and. xIndex .lt. nxFine)
            xIndex = xIndex+1;
        end do
        xIndex = xIndex-1;
        enclosingFineRectX(i) = xIndex;
    enddo

    do i = 1, nyCoarse        
        yIndex = 1;
        do while (yrefFine(yIndex) .le. yrefCoarse(i) .and. yIndex .lt. nyFine)
            yIndex = yIndex+1;
        end do
        yIndex = yIndex-1;
        enclosingFineRectY(i) = yIndex;
    enddo
        
    return
    end subroutine findEnclosingFineRects
    
!=================================================================    
    subroutine setUpInitialNesting
!=================================================================    
    use Q4_globals
    implicit none
    
    integer :: i
    real(8) :: stepSize

!    do i=1,nxm(1)
!	    xref1(i)=xref(2*i-1)
!    enddo
!
!    do i=1,nxm(2)
!	    xref2(i)=xref1(2*i-1)
!    enddo
!
!    do i=1,nxm(3)
!	    xref3(i)=xref2(2*i-1)
!    enddo
!
!    do i=1,nxm(4)
!	    xref4(i)=xref3(2*i-1)
!    enddo

    stepSize = 1.0 / (nx1-1)
    do i=1,nxm(1)
	    xref1(i) = stepSize * (i-1)
    enddo

    stepSize = 1.0 / (nx2-1)
    do i=1,nxm(2)
	    xref2(i) = stepSize * (i-1)
    enddo

    stepSize = 1.0 / (nx3-1)
    do i=1,nxm(3)
	    xref3(i) = stepSize * (i-1)
    enddo

    stepSize = 1.0 / (nx4-1)
    do i=1,nxm(4)
	    xref4(i) = stepSize * (i-1)
    enddo



    !semiCoarsening in Y direction
    stepSize = yl / (ny1-1)
    do i=1,ny1
        yref1(i) = stepSize * (i-1)
    enddo

    stepSize = yl / (ny2-1)
    do i=1,ny2
        yref2(i) = stepSize * (i-1)
    enddo

    stepSize = yl / (ny3-1)
    do i=1,ny3
	    yref3(i) = stepSize * (i-1)
    enddo

    stepSize = yl / (ny4-1)
    do i=1,ny4
	    yref4(i) = stepSize * (i-1)
    enddo
!    

!    do i=1,nym(1)
!	    yref1(i)=yref(2*i-1)
!    enddo
!
!    do i=1,nym(2)
!        yref2(i)=yref1(2*i-1)
!    enddo
!
!    do i=1,nym(3)
!	    yref3(i)=yref2(2*i-1)
!    enddo
!
!    do i=1,nym(4)
!	    yref4(i)=yref3(2*i-1)
!    enddo
    
    return
    end subroutine setUpInitialNesting
    
    
!=================================================================
    subroutine ave_multi
!=================================================================

    use Q4_globals
    implicit none
    

    !coordinates for the coarse grids
    call SetUpMGNesting()

    call ave_height(cohimx,cohjmx,himax,himin,hjmax,hjmin, &
     		        recssi,recssj,xref,yref,nx,ny,nx,20,ny,10)

    call ave_height(cohimx1,cohjmx1,himax1,himin1,hjmax1,hjmin1, &
     		        recssi1,recssj1,xref1,yref1,nx1,ny1,nxm(1),20,nym(1),10)

    call ave_height(cohimx2,cohjmx2,himax2,himin2,hjmax2,hjmin2, &
     		        recssi2,recssj2,xref2,yref2,nx2,ny2,nxm(2),20,nym(2),10)

    call ave_height(cohimx3,cohjmx3,himax3,himin3,hjmax3,hjmin3, &
     		        recssi3,recssj3,xref3,yref3,nx3,ny3,nxm(3),20,nym(3),10)

    call ave_height(cohimx4,cohjmx4,himax4,himin4,hjmax4,hjmin4, &
                    recssi4,recssj4,xref4,yref4,nx4,ny4,nxm(4),20,nym(4),10)

    return
    end

!====================================================================
      subroutine ave_height(cohimx0,cohjmx0,himax0,himin0, &
                            hjmax0,hjmin0,recssi0,recssj0, &
                            xref0,yref0,ndx,ndy,nx0,nxd,ny0,nyd)
!====================================================================

    use Q4_globals
    implicit none
    
    real(8) :: cohimx0,cohjmx0,himax0,himin0, &
               hjmax0,hjmin0,recssi0,recssj0, xref0,yref0
    integer :: ndx,ndy,nx0,nxd,ny0,nyd

    integer :: i, j, jd, id
    real(8) :: xv, yv, dxf, dyf, recess, xface, yface
    
    dimension xface(nx+1),yface(ny+1), cohimx0(ndx,ndy),cohjmx0(ndx,ndy), &
     		    himax0(ndx,ndy),himin0(ndx,ndy), hjmax0(ndx,ndy),hjmin0(ndx,ndy), &
     		    recssi0(ndx,ndy), recssj0(ndx,ndy), xref0(ndx),yref0(ndy)

    !xface and y face hold finite volume boundaries/walls
    !xface(i) holds the x location of the face in front of grid point xref0(i)
    !at boundary xface holds boundary
    !thus recess0(i) holds height data for face in front of xref0(i)
    do i = 2,nx0
	    xface(i)=(xref0(i)+xref0(i-1))/2.d0
    enddo
    xface(1) = xref(1)
    xface(nx0+1) = xref0(nx0)

    do i = 2, ny0
	    yface(i)=(yref0(i)+yref0(i-1))/2.d0
    enddo
    yface(1) = yref(1)
    yface(ny0+1) = yref0(ny0)

    do j = 1, ny0
	    dyf=(yface(j+1)-yface(j))/dfloat(nyd)  !nyd = step size  yface goes from 1, ny+1
	    do i=2,nx0
	        recssi0(i,j)=0.d0
	        himin0(i,j)=1.d10
	        himax0(i,j)=-1.d10
	        !microstep along volume boundary, getting height
	        do jd = 1, nyd
	            xv = xface(i)
	            yv = yface(j) + dyf*(dfloat(jd) - 0.5d0)
	            call pointrecess(recess,xv,yv)
	            recssi0(i,j) = recssi0(i,j) + recess   !accumulate recess (for avg recess calculation)
	            if (himax0(i,j) .lt. recess) then
	                himax0(i,j) = recess  !save max and min recess for this volume
	            endif
	            if (himin0(i,j) .gt. recess) then
	                himin0(i,j) = recess
	            endif
	        enddo
	        recssi0(i,j) = recssi0(i,j) / dfloat(nyd)    !average recess for this volume boundary
	        
	        if( (himax0(i,j) - himin0(i,j)) .gt. 1.d-5) then
	            cohimx0(i,j) = (recssi0(i,j) - himin0(i,j)) / (himax0(i,j)-himin0(i,j))
	        else
	            cohimx0(i,j)=0.5d0
	        endif
	    enddo
    enddo

    do i=1,nx0
        dxf=(xface(i+1)-xface(i))/dfloat(nxd)
        do j=2,ny0
	        recssj0(i,j)=0.d0
	        hjmin0(i,j)=1.d10
	        hjmax0(i,j)=-1.d10
	        do id=1,nxd
	            xv=xface(i)+dxf*(dfloat(id)-0.5d0)
	            yv=yface(j)
	            call pointrecess(recess,xv,yv)
	            recssj0(i,j)=recssj0(i,j)+recess
	            if(hjmax0(i,j).lt.recess)hjmax0(i,j)=recess
	            if(hjmin0(i,j).gt.recess)hjmin0(i,j)=recess
	        enddo
	        recssj0(i,j)=recssj0(i,j)/dfloat(nxd)
	        if((hjmax0(i,j)-hjmin0(i,j)).gt.1.d-5)then
	            cohjmx0(i,j)=(recssj0(i,j)-hjmin0(i,j))/(hjmax0(i,j)-hjmin0(i,j))
	        else
	            cohjmx0(i,j)=0.5d0
	        endif
        enddo
    enddo

    return
    end
      

!============================
    subroutine gethx
!============================
    use Q4_globals
    use NestingInfo
    implicit none

    integer :: i, j
    
!   this subroutine calculates the steady state bearing separation
!   used when surface is smooth and slider separation is sought

    do j=1,ny
        do i=1,nx
            hnew(i,j)=h0 + hx0*xref(i) - hy*(yref(j) - yl/2.d0)
        enddo
    enddo

	call checkCrash
    if (crash) return     

!   restrict to coarse grid


    call restrictPressure(hnew, hnew1, xref, yref, xref1, yref1, nx, ny, &
                          nx1, ny1, enclosingFineRectX1, enclosingFineRectY1)
    call restrictPressure(hnew1, hnew2, xref1, yref1, xref2, yref2, nx1, ny1, &
                          nx2, ny2, enclosingFineRectX2, enclosingFineRectY2)
    call restrictPressure(hnew2, hnew3, xref2, yref2, xref3, yref3, nx2, ny2, &
                          nx3, ny3, enclosingFineRectX3, enclosingFineRectY3)
    call restrictPressure(hnew3, hnew4, xref3, yref3, xref4, yref4, nx3, ny3, &
                          nx4, ny4, enclosingFineRectX4, enclosingFineRectY4)

!    call restrict(hnew,hnew1,nx,ny,nx,ny,nx1,ny1)
!    call restrict(hnew1,hnew2,nxm(1),nym(1),nx1,ny1,nx2,ny2)
!    call restrict(hnew2,hnew3,nxm(2),nym(2),nx2,ny2,nx3,ny3)
!    call restrict(hnew3,hnew4,nxm(3),nym(3),nx3,ny3,nx4,ny4)

    return
    end



!=============================
    subroutine getHeight
!=============================
!     height information

    use Q4_globals
    implicit none

    integer :: i, j
    real(8) :: xv, yv

    do i = 1,nx
	    do j = 1,ny
            xv = xref(i)
            yv = yref(j)
            ! make sure inside boundary
            if (i.eq.1)  xv = 1.d-8
            if (i.eq.nx) xv = 1.d0 - 1.d-8
            if (j.eq.1)  yv = 1.d-8
            if (j.eq.ny) yv = yl - 1.d-8
            call pointrecess(h(i,j),xv,yv)
        enddo
    enddo

    return
    end
    
    
!=============================================
    subroutine checkCrash
!=============================================
! make sure no negative value occur.
! otherwise log function in flow calculation
! may generate floating point error

    use Q4_globals
    implicit none

    integer i, j

    do i = 2, nx
        do j = 2, ny - 1
            if (himin(i,j)+(hnew(i,j)+hnew(i-1,j))*0.5d0 .lt. 2.d-10) then
                crash = .true.
                return
            endif
        enddo
    enddo
    
    do j = 2, ny
        do i = 2, nx - 1
            if (hjmin(i,j)+(hnew(i,j)+hnew(i,j-1))*0.5d0 .lt. 2.d-10) then
                crash = .true.
                return
            endif
        enddo
    enddo

    crash = .false.
    
    return
    end    



!====================================
    subroutine calc_bearing_number
!====================================
!
!     Notes from BC:
!     OK, so take this with a grain of salt cause I'm a moronic cretin:
!     bearx and beary are actually NOT the bearing numbers
!     as they are defined in most thesises: ((6*vis*(U or V)*L)/(pa*hm*hm))
!     they seem to be simplified to only include the non-dimensionalized x and y velocity
!     components (U or V)at points on the slider
!

    use Q4_globals
    implicit none
    
    integer :: i, j
    real(8) :: xs0, yc0, yc1, yc2, rnew1, skettl, rad

    rad = ra/xl  !radial position normalized to slider length

    do i = 1,nx
	    xs0 = 0.5-xref(i)
	    do j = 1,ny
	        yc0 = yl/2.d0-yref(j)
	        yc1 = rad - yc0*cos(ske) - xs0*sin(ske)
	        yc2 = xs0*cos(ske) - yc0*sin(ske)
	        rnew1 = dsqrt(yc1*yc1+yc2*yc2)/rad
	        skettl = ske - atan(yc2/yc1)
	        bearx(i,j) = u0*cos(skettl)*rnew1/xl
	        beary(i,j) = - u0*sin(skettl)*rnew1/xl
	    end do
    end do
    
    do i = 1,nx1
	    xs0 = 0.5-xref1(i)
	    do j = 1,ny1
	        yc0 = yl/2.d0-yref1(j)
	        yc1 = rad - yc0*cos(ske) - xs0*sin(ske)
	        yc2 = xs0*cos(ske) - yc0*sin(ske)
	        rnew1 = dsqrt(yc1*yc1+yc2*yc2)/rad
	        skettl = ske - atan(yc2/yc1)
	        bearx1(i,j) = u0*cos(skettl)*rnew1/xl
	        beary1(i,j) = - u0*sin(skettl)*rnew1/xl
	    end do
    end do

    do i = 1,nx2
	    xs0 = 0.5-xref2(i)
	    do j = 1,ny2
	        yc0 = yl/2.d0-yref2(j)
	        yc1 = rad - yc0*cos(ske) - xs0*sin(ske)
	        yc2 = xs0*cos(ske) - yc0*sin(ske)
	        rnew1 = dsqrt(yc1*yc1+yc2*yc2)/rad
	        skettl = ske - atan(yc2/yc1)
	        bearx2(i,j) = u0*cos(skettl)*rnew1/xl
	        beary2(i,j) = - u0*sin(skettl)*rnew1/xl
	    end do
    end do

    do i = 1,nx3
	    xs0 = 0.5-xref3(i)
	    do j = 1,ny3
	        yc0 = yl/2.d0-yref3(j)
	        yc1 = rad - yc0*cos(ske) - xs0*sin(ske)
	        yc2 = xs0*cos(ske) - yc0*sin(ske)
	        rnew1 = dsqrt(yc1*yc1+yc2*yc2)/rad
	        skettl = ske - atan(yc2/yc1)
	        bearx3(i,j) = u0*cos(skettl)*rnew1/xl
	        beary3(i,j) = - u0*sin(skettl)*rnew1/xl
	    end do
    end do

    do i = 1,nx4
	    xs0 = 0.5-xref4(i)
	    do j = 1,ny4
	        yc0 = yl/2.d0-yref4(j)
	        yc1 = rad - yc0*cos(ske) - xs0*sin(ske)
	        yc2 = xs0*cos(ske) - yc0*sin(ske)
	        rnew1 = dsqrt(yc1*yc1+yc2*yc2)/rad
	        skettl = ske - atan(yc2/yc1)
	        bearx4(i,j) = u0*cos(skettl)*rnew1/xl
	        beary4(i,j) = - u0*sin(skettl)*rnew1/xl
	    end do
    end do

    uflag = OFF
    
    
!   restrict to multigrid

!    call restrict(bearx,bearx1,nx,ny,nx,ny,nx1,ny1)
!    call restrict(bearx1,bearx2,nxm(1),nym(1),nx1,ny1,nx2,ny2)
!    call restrict(bearx2,bearx3,nxm(2),nym(2),nx2,ny2,nx3,ny3)
!    call restrict(bearx3,bearx4,nxm(3),nym(3),nx3,ny3,nx4,ny4)
!
!    call restrict(beary,beary1,nx,ny,nx,ny,nx1,ny1)
!    call restrict(beary1,beary2,nxm(1),nym(1),nx1,ny1,nx2,ny2)
!    call restrict(beary2,beary3,nxm(2),nym(2),nx2,ny2,nx3,ny3)
!    call restrict(beary3,beary4,nxm(3),nym(3),nx3,ny3,nx4,ny4)

    return
    end

      
!======================================================
    subroutine calcReynoldsCoefs
!======================================================
    use Q4_globals
    implicit none
    
    if(hflag.eq.ON) then
	    call getHeight
	    call ave_multi
    endif
    
    if(uflag.eq.ON) then
        call calc_bearing_number
    endif  
    
    hflag = OFF
    uflag = OFF
      
    return
    end subroutine calcReynoldsCoefs
      
      
!======================================================
      subroutine stiffness
!======================================================
!c  stiffness matrix is 3x3
!c  The slider motion is decomposed into three degrees
!c  of freedom: 
!c       Z: vertical motion(not just TEC),  
!c       P: pitch about pivot point(xf0, yf0)
!c       R: roll about pivot point   
!c  Units and sign convention for increments:
!c       Z: positive displacement (nm) reduces GAP,         
!c          corresponding to increased LOAD (g) .
!c       P: positive displacement (uRad) increases pitch,
!c          corresponding to increased P-TORQUE (uN-M).
!c       R: positive displacement (uRad) increases roll,
!c          corresponding to increased R-TORQUE (uN-M).
!c  Note the sign of roll:
!c       OD is closer to disk than ID for positive roll.
!c       This has been changed from last version.
!c  Also, the diagonal elements of the stiffness
!c  matrix are now normally positive because of the new
!c  sign convention. 
!c  Let L(1) = LOAD, L(2) = P-TORQUE, L(3) = R-TORQUE
!c  and D(1) = Z, D(2) = P, D(3) = R, then the stiffness
!c  matrix are the following differentials:
!c   
!c           stiff(i,j) = dL(i)/dD(j) 
!c  
!c  For example, stiff(2, 3) = 3 ==> 3 uN-M P-TORQUE required to
!c  increase R (roll) by 1 uRad.

    use Q4_globals
    implicit none

    real(8) :: stiff(3,3)

    if(crash) then
        call matInit(stiff, 3, 3, 3, 3, 0.d0)
    else
        call calc_deriv(5.0d-2)

        stiff(1,1)=-jac(1,1)*1.d-6/hm    
        stiff(1,2)=-(jac(1,2)+(.5d0-xf0)*jac(1,1))*1.d-3*xl/hm 
        stiff(1,3)=(jac(1,3)+yf0*jac(1,1))*1.d-3*xl/hm  
        stiff(2,1)=-jac(2,1)*1.d-3*xl/hm*9.81d0  
        stiff(2,2)=-(jac(2,2)+(.5d0-xf0)*jac(2,1))*xl*xl/hm*9.81d0 
        stiff(2,3)=(jac(2,3)+yf0*jac(2,1))*xl*xl/hm*9.81d0 
        stiff(3,1)=-jac(3,1)*1.d-3*xl/hm*9.81d0 
        stiff(3,2)=-(jac(3,2)+(.5d0-xf0)*jac(3,1))*xl*xl/hm*9.81d0
        stiff(3,3)=(jac(3,3)+yf0*jac(3,1))*xl*xl/hm*9.81d0
    endif

    call wrStiff(numRes, stiff)
    call wrStiff(6, stiff)

    return
    end

!=================================================
    subroutine wrStiff(numfile, stiff)
!=================================================
    use Q4_globals
    implicit none

    integer :: numfile, j
    real(8) :: stiff(3, 3)

    write(numfile,*)
    write(numfile,*)'STIFFNESS MATRIX'
    write(numfile,'(a15,3e15.6)')'LOAD(G)',(stiff(1,j),j=1,3)
    write(numfile,'(a15,3e15.6)')'P-TORQUE(uN-M)',(stiff(2,j),j=1,3)
    write(numfile,'(a15,3e15.6)')'R-TORQUE(uN-M)',(stiff(3,j),j=1,3)
    write(numfile,'(4a15)')' ','HEIGHT(NM)','PITCH(uRAD)','ROLL(uRAD)'

    return
    end


!returns the number of lines in the input file
!input file must not be open but must exist       
!=======================================================
    integer function GetNumLinesInFile(filename)
!=======================================================
    implicit none
    character*(*) :: filename
    character curChar
    integer iFile, numLines
    
    numLines = 0
    iFile = 314
    open(iFile, err=456, file=filename, status='unknown')
    do while (.true.)
        read(iFile, '(a1)', end=456, err=456) curChar
        numLines = numLines + 1
    enddo
    
456 continue
    close(iFile)
    GetNumLinesInFile = numLines
    
    return
    end
    

!===============================================================
    subroutine LoadMatrix(matrix, ix, iy, filename, iSuccess)
!===============================================================
    implicit none
    character*(*) :: filename
    integer :: i, j, ix, iy, numLinesInFile, GetNumLinesInFile
    integer :: iSuccess
    real(8) :: matrix(ix, iy)

    iSuccess = 0
    !first make sure that we have the correct number of lines in the file
    numLinesInFile = GetNumLinesInFile(filename)
    if (numLinesInFile .ne. iy) goto 654
    
    !open file (file must exist)
    open(41, err=654, file=filename, status='unknown')
    do j = 1, iy
        read(41, *, err=654) (matrix(i, j), i = 1, ix)
    enddo       
    
    iSuccess = 1
    return

654 continue  !error handler
    write(*,*) 'Warning: errors encountered when opening or reading:', filename
!    stop 'exiting...'
    iSuccess = 0
    
    end subroutine LoadMatrix