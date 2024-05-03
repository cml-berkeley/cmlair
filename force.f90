!==========================================
    subroutine getfxp(iwrite)
!==========================================
    use Q4_globals
    implicit none

    integer, intent(in) :: iwrite
    
    integer :: i, im1, ip1, j, jm1, jp1, mini, minj, maxi, maxj
    real(8) :: t, r, pmax, pmin, dely, delx, dxdy, a
    real(8) :: fvdw, fvdwdash, tvdw1, tvdw2, rvdw1, rvdw2
    real(8) :: fposvdw, fnegvdw, hamcon1, hamcon2, vdwnormf, tempgap
    real(8) :: vdw1, vdw2, xfvdw, yfvdw

! Added suspension stiffnesses option 10/21/04
! we'll update suspension torques here
! I think this is probably the best place to do this...
    if (iUseStiffnesses .eq. 1) then
        call UpdateSuspensionTorques()
    endif

    f=0.d0
    t=0.d0
    r=0.d0
    fpos =0.d0
    fneg=0.d0
    pmax=-1.d0
    pmin=10.d0

    do j=2,nym1
	    jm1=j-1
	    jp1=j+1
	    dely=(yref(jp1)-yref(jm1))/2.d0
	    do i=2,nxm1
	        im1=i-1
	        ip1=i+1
	        delx=(xref(ip1)-xref(im1))/2.d0
	        dxdy=delx*dely
	        a=(p(i,j)-1.d0)*dxdy
	        if(pmax.lt.p(i,j))then
	            pmax=p(i,j)
	            maxi=i
	            maxj=j
	        endif
	        if(pmin.gt.p(i,j))then
	            pmin=p(i,j)
	            mini=i
	            minj=j
	        endif
	        t=t+a*(xref(i) - 0.5d0 - xf0)
	        r=r+a*(yref(j) - yl/2.d0 - yf0)
	        if (a.le.0.d0) then
	            fneg=fneg+a
	        else
	            fpos=fpos+a
	        endif
	    enddo
    enddo

	f = fpos + fneg

	pmax = pmax - 1.d0
    pmin = pmin - 1.d0

    xf= t*p0xl
    yf= r*p0xl
    f = f*p0xl
    fpos=fpos*p0xl
    fneg=fneg*p0xl

    if (elecpot .ne. 0) call ELECTROSTATIC

      
    !IMF
    !added by vineet to incorporate the effect of vanderwall forces
	if (iUseIMF .eq. 1) then
        
	    fvdw=0.d0
	    fvdwdash=0.d0
        tvdw1=0.d0
	    tvdw2=0.d0
        rvdw1=0.d0
	    rvdw2=0.d0
        fposvdw =0.d0
        fnegvdw=0.d0
	
        ! hamaker constants
	    hamcon1 = ahc/6/3.141592653589793
	    hamcon2 = bhc/45/3.141592653589793

	    vdwnormf = xl*xl/9.81d0 !BC --  this is actually an un-normalization factor
        ! write(*,*)'***AHC is', ahc, 'BHC is ', bhc
        
        do j=2,nym1
	        jm1=j-1
	        jp1=j+1
	        dely=(yref(jp1)-yref(jm1))/2.d0
	        do i=2,nxm1
	            im1=i-1
	            ip1=i+1
	            delx=(xref(ip1)-xref(im1))/2.d0
	            dxdy=delx*dely

		        fvdwdash=0.d0

	            tempgap = hm*(recssi(i,j) + hnew(i,j))

		        vdw1=hamcon1*dxdy/(tempgap**3)
		        vdw2=hamcon2*dxdy/(tempgap**9)

		        fvdwdash=fvdwdash-vdw1+vdw2

                ! pitch torques from van der walls force

	            tvdw1=tvdw1-vdw1*(xref(i) - 0.5d0 - xf0)
	            tvdw2=tvdw2+vdw2*(xref(i) - 0.5d0 - xf0)

                ! roll torques from vdw
	            rvdw1=rvdw1-vdw1*(yref(j) - yl/2.d0 - yf0)
	            rvdw2=rvdw2+vdw2*(yref(j) - yl/2.d0 - yf0)

                ! if (fvdw.le.0.d0) then
                !   fnegvdw=fnegvdw-vdw1
                ! else
                !   fposvdw=fposvdw+vdw2
                ! endif
		        fnegvdw=fnegvdw-vdw1
		        fposvdw=fposvdw+vdw2
		        vdwMolecularForceMap(i,j) = (vdw2 - vdw1) / (dxdy*9.81) !kg/(m*m)
	        enddo
        enddo

        fvdw = fposvdw + fnegvdw
        !pull out only vdw forces before electrostatic is added in
        fvdw_output = fvdw*vdwnormf
      	  
      	  
	    xfvdw = (tvdw1+tvdw2)*vdwnormf
        yfvdw = (rvdw1+rvdw2)*vdwnormf
        fvdw = fvdw*vdwnormf
        fposvdw = fposvdw*vdwnormf
        fnegvdw = fnegvdw*vdwnormf

	    xf = xf + xfvdw
	    yf = yf + yfvdw
	    f = f + fvdw

	    fpos = fpos + fposvdw
  	    fneg = fneg + fnegvdw
	endif !end of IMF

    !we don't need to modify pressure due to humidity and temperature here
    !it was done before calling getfxp
    call shear

    xf=xf+fsp
    yf=yf+fsr
    ! contact force
    call getrcf(1)
    f=f+fcr
    xf=xf+txr
    yf=yf+tyr

    if(iwrite.eq.1) then
	    write(*,*)
	    write(*,*)'POSITIVE FORCE(G) = ', fpos*1000
	    write(*,*)'NEGATIVE FORCE(G) = ', fneg*1000
	    write(*,675) (f-f0)*1.d3
	    write(*,676) (xf-xfs)*xl*9.81d6
	    write(*,677) (yf-yfs)*xl*9.81d6
	    write(*,678) Zmom
	    write(*,615)'PMAX =',pmax,'AT (',xref(maxi)*xl*1.d3,',',yref(maxj)*xl*1.d3,')'
	    write(*,615)'PMIN =',pmin,'AT (',xref(mini)*xl*1.d3,',',yref(minj)*xl*1.d3,')'
	    write(*,'(/)')
    endif

615 format(a7,g12.5,a4,f8.3,a1,f8.3,a1)
655 format(2g12.4/)
675 format(' FORCE ERR IN GRAM    = ',g12.5)
676 format(' P-TORQUE ERR IN uN-M = ',g12.5)
677 format(' R-TORQUE ERR IN uN-M = ',g12.5)
678 format(' Z-MOMENT     IN uN-M = ',g12.5)

    return
    end subroutine getfxp


     
!=============================================
    subroutine UpdateSuspensionTorques
!=============================================
    use Q4_globals
    implicit none
      
    xfs = (PSA - (-hx0*hm/xl)) * Pitch_Stiffness  !  N*M not uN*M
    yfs = (RSA - (hy*hm/xl)) * Roll_Stiffness
    
    xfs = xfs / xl / 9.81d0  !get back to length-normalized kg
    yfs = yfs / xl / 9.81d0
    
    return
    end subroutine UpdateSuspensionTorques



!=============================================
    subroutine shear
!=============================================
    use Q4_globals
    use ShearArrays
    implicit none
    
    real(8) :: fsp1, fsr1, fsp2, fsr2, coeff1, coeff2, dy, dx
    real(8) :: tempfsp1, tempy, tempfsr1, tempx, dxdy, temp, bear, height, blah
    integer :: i, ip1, im1, j, jp1, jm1

    fsp1=0.d0
    fsr1=0.d0
    fsp2=0.d0
    fsr2=0.d0

    coeff1=zl*xl*p0*hm/9.81d0
    coeff2=zl*vis1*xl*xl*xl/hm/9.81d0
!
!     fsp1 is first term in integral ((H/2)*(dp/dx))
!     fsp2 is second term (Bear (actually just the velocity))/(slip stuff)
!     Bearx and beary are actually not the bearing number but the velocity in x and y
!     the other params in the bearing number are included in coeff1 and coeff2
!
    Zmom = 0.d0
    nxm1=nx-1
    nym1=ny-1
    
    ! poiseuille in X Direction
    do i=1,nxm1
	    ip1 = i+1
	    im1 = i-1
	    if(i.eq.1) im1=i
	    
	    dx = xref(ip1) - xref(i)

	    do j=1,ny
	        jm1 = j-1
	        jp1 = j+1
	        if (j .eq. 1) jm1 = j
	        if (j .eq. ny) jp1 = j
	        dy = (yref(jp1) - yref(jm1)) / 2.d0
	        
	        height = recssi(ip1,j) + (hnew(i,j)+hnew(ip1,j)) / 2.d0
            tempfsp1= -(p(ip1,j)-p(i,j)) * dy * height / 2.d0 
	        fsp1=fsp1+tempfsp1
	        
	        !shear stress = (shear_force / area)
	        PoisilleShearX(i,j) = (9.81 * tempfsp1 * coeff1 / zl) / (dx*dy * xl * xl * p0) !Atmospheres
	        
		    tempy = (yref(jp1)+yref(jm1))/2.d0
		    Zmom = Zmom + tempfsp1*(0.5d0*yl+yf0-tempy)*xl*coeff1
	    enddo
    enddo

    ! poiseuille in Y Direction
    do j = 1,nym1
	    jp1=j+1
	    jm1 = j-1
	    if(j.eq.1)jm1=j
	    dy = yref(jp1) - yref(j)
	    
	    do i = 1,nx
	        im1=i-1
	        ip1=i+1
            if (i .eq. 1) im1 = i
            if (i .eq. nx) ip1 = i        
	        dx = (xref(ip1) - xref(im1)) / 2.d0

            height = recssj(i,jp1) + (hnew(i,j)+hnew(i,jp1)) / 2.d0
		    tempfsr1 = -(p(i,jp1)-p(i,j)) * dx * height /2.d0
	        fsr1=fsr1 + tempfsr1
	        
	        !Shear Stress
	        PoisilleShearY(i,j) = (9.81 * tempfsr1 * coeff1 / zl) / (dx*dy * xl * xl * p0) !Atmospheres
		    tempx = (xref(ip1)+xref(im1))/2.d0
		    Zmom = Zmom - tempfsr1*(0.5d0+xf0-tempx)*xl*coeff1
	    enddo
    enddo

    ! couette shear in X Direction
    do i=1,nxm1
	    ip1 = i+1
	    im1 = i-1
	    if(i.eq.1) im1=i
	    
	    dx = xref(ip1) - xref(i)

	    do j=1,ny
	        jm1 = j-1
	        jp1 = j+1
	        if (j .eq. 1) jm1 = j
	        if (j .eq. ny) jp1 = j
	        dy = (yref(jp1) - yref(jm1)) / 2.d0
	        
	        dxdy=dx*dy
	        height = recssi(ip1,j) + (hnew(i,j)+hnew(ip1,j)) / 2.d0
	        temp=dxdy/(height+2*corCoef*al)
	        
	        bear = (bearx(i,j) + bearx(ip1, j)) / 2.d0
	        fsp2=fsp2 + bear * temp
	        
	        CouetteShearX(i,j) = (9.81 * bear * temp * coeff2 / zl) / (dxdy * xl * xl * p0)
	        
	        tempy = (yref(jp1)+yref(jm1))/2.d0
		    Zmom = Zmom + bearx(i,j)*temp*(0.5d0*yl+yf0-tempy)*xl*coeff2
	    enddo
    enddo

    ! couette shear in Y Direction
    do j = 1,nym1
	    jp1=j+1
	    jm1 = j-1
	    if(j.eq.1)jm1=j
	    dy = yref(jp1) - yref(j)
	    
	    do i = 1,nx
	        im1=i-1
	        ip1=i+1
            if (i .eq. 1) im1 = i
            if (i .eq. nx) ip1 = i        
	        dx = (xref(ip1) - xref(im1)) / 2.d0
	        
	        dxdy=dx*dy
            height = recssj(i,jp1) + (hnew(i,j)+hnew(i,jp1)) / 2.d0
            temp=dxdy/(height+2*corCoef*al)
            
            bear = (beary(i,j) + beary(i, jp1)) / 2.d0
            fsr2=fsr2 + bear*temp
            
            CouetteShearY(i,j) = (9.81 * bear * temp * coeff2 / zl) / (dxdy * xl * xl * p0)
            
            tempx = (xref(ip1)+xref(im1))/2.d0
		    Zmom = Zmom - beary(i,j)*temp*(0.5d0+xf0-tempx)*xl*coeff2
	    enddo
    enddo

    fsp1=fsp1*coeff1
    fsr1=fsr1*coeff1
    fsp2=fsp2*coeff2
    fsr2=fsr2*coeff2

    fsp=fsp1+fsp2
    fsr=fsr1+fsr2

    Zmom = Zmom*9.81d6/zl
    return
    end !subroutine shear

!=============================================
    subroutine shearfujitsu
!=============================================
    use Q4_globals
    use ShearArrays
    implicit none
    
    real(8) :: fsp1, fsr1, fsp2, fsr2, coeff1, coeff2
    real(8) :: dx_pressure, dx_area, dy_pressure, dy_area, dpdx, dpdy
    real(8) :: tempfsp1, tempy, tempfsr1, tempx, dxdy, temp, bear, height
    integer :: i, ip1, im1, j, jp1, jm1

    fsp1=0.d0
    fsr1=0.d0
    fsp2=0.d0
    fsr2=0.d0

    coeff1=zl*xl*p0*hm/9.81d0
    coeff2=zl*vis1*xl*xl*xl/hm/9.81d0
!
!     fsp1 is first term in integral ((H/2)*(dp/dx))
!     fsp2 is second term (Bear (actually just the velocity))/(slip stuff)
!     Bearx and beary are actually not the bearing number but the velocity in x and y
!     the other params in the bearing number are included in coeff1 and coeff2
!
    Zmom = 0.d0
    nxm1=nx-1
    nym1=ny-1
    
    ! poiseuille in X Direction
    ! dx_pressure is simply for the dp/dx term in the shear calculation
    ! dx_area is the area over which we're integrating
    ! We have 3 regions: edge of slider at the start, all middle regions, edge of slider at end
    do i=1,nx+1
	    ip1 = i+1
	    im1 = i-1

        if(i.eq.1) then
            dx_pressure = xref(2) - xref(1)
            dx_area = (xref(2) - xref(1)) / 2.d0 
        elseif (i .eq. nx+1) then
            dx_pressure = xref(nx) - xref(nx-1)
            dx_area = (xref(nx) - xref(nx-1)) / 2.d0 
        else
            dx_pressure = xref(i) - xref(im1)
            dx_area = dx_pressure
        endif

	    do j=1,ny
	        jm1 = j-1
	        jp1 = j+1
	        if (j .eq. 1) jm1 = j
	        if (j .eq. ny) jp1 = j
            dy_area = (yref(jp1) - yref(jm1)) / 2.d0
            
            if (i .eq. 1) then
                height = h(1, j) +hnew(1, j)
                dpdx = -(p(2, j)-p(1,j)) / dx_pressure
            elseif (i .eq. nx+1) then
                height = h(nx, j) + hnew(nx, j)
                dpdx = -(p(nx, j)-p(nx-1 ,j)) / dx_pressure
            else
                height = recssi(i,j) + (hnew(i,j)+hnew(im1,j)) / 2.d0
                dpdx = -(p(i,j)-p(im1,j)) / dx_pressure
            endif

            tempfsp1= dpdx * dy_area * dx_area * height / 2.d0
	        if (i .gt. 1 .and. i .lt. nx+1) fsp1=fsp1+tempfsp1
	        
	        !shear stress = (shear_force / area)
	        PoisilleShearX(i, j) = (9.81 * tempfsp1 * coeff1 / zl) / (dy_area * dx_area * xl * xl * p0) !Atmospheres
	        
		    tempy = yref(j)
		    Zmom = Zmom + tempfsp1*(0.5d0*yl+yf0-tempy)*xl*coeff1
	    enddo
    enddo

    ! poiseuille in Y Direction
    do j = 1,ny+1
	    jp1=j+1
	    jm1 = j-1
        
        if(j.eq.1) then
            dy_pressure = yref(2) - yref(1)
            dy_area = (yref(2) - yref(1)) / 2.d0 
        elseif (j .eq. ny+1) then
            dy_pressure = yref(ny) - yref(ny-1)
            dy_area = (yref(ny) - yref(ny-1)) / 2.d0 
        else
            dy_pressure = yref(j) - yref(jm1)
            dy_area = dy_pressure
        endif
	    
	    do i = 1,nx
	        im1=i-1
	        ip1=i+1
            if (i .eq. 1) im1 = i
            if (i .eq. nx) ip1 = i        
	        dx_area = (xref(ip1) - xref(im1)) / 2.d0

            if (j .eq. 1) then
                height = h(i, 1) + hnew(i, 1)
                dpdy = -(p(i, 2)-p(i,1)) / dy_pressure
            elseif (j .eq. ny+1) then
                height = h(i, ny) + hnew(i, ny)
                dpdy = -(p(i, ny)-p(i ,ny-1)) / dy_pressure
            else
                height = recssj(i,j) + (hnew(i,jm1)+hnew(i,j)) / 2.d0
                dpdy = -(p(i,j)-p(i,jm1)) / dy_pressure
            endif

		    tempfsr1 = dpdy * dy_area * dx_area * height /2.d0
		    if (j .gt. 1 .and. j .lt. ny+1) fsr1=fsr1 + tempfsr1
	        
	        !Shear Stress
	        PoisilleShearY(i,j) = (9.81 * tempfsr1 * coeff1 / zl) / (dx_area*dy_area * xl * xl * p0) !Atmospheres
		    tempx = xref(i)
		    Zmom = Zmom - tempfsr1*(0.5d0+xf0-tempx)*xl*coeff1
	    enddo
    enddo

    ! couette shear in X Direction
    do i=1,nx+1
	    ip1 = i+1
	    im1 = i-1

        if(i.eq.1) then
            dx_area = (xref(2) - xref(1)) / 2.d0 
        elseif (i .eq. nx+1) then
            dx_area = (xref(nx) - xref(nx-1)) / 2.d0 
        else
            dx_area = xref(i) - xref(im1)
        endif
        
	    do j=1,ny
	        jm1 = j-1
	        jp1 = j+1
	        if (j .eq. 1) jm1 = j
	        if (j .eq. ny) jp1 = j
	        dy_area = (yref(jp1) - yref(jm1)) / 2.d0
	        
	        dxdy=dx_area*dy_area
            
            if (i .eq. 1) then
                height = h(1, j) + hnew(1, j)
                bear = bearx(1, j)
            elseif (i .eq. nx+1) then
                height = h(nx, j) + hnew(nx, j)
                bear = bearx(nx, j)
            else
                height = recssi(i,j) + (hnew(i,j)+hnew(im1,j)) / 2.d0
                bear = (bearx(im1,j) + bearx(i, j)) / 2.d0
            endif
            
	        temp=dxdy/(height+2*corCoef*al)
	        if (i .gt. 1 .and. i .lt. nx+1) fsp2=fsp2 + bear * temp
	        
	        CouetteShearX(i,j) = (9.81 * bear * temp * coeff2 / zl) / (dxdy * xl * xl * p0)
	        
	        tempy = yref(j)
		    Zmom = Zmom + bearx(i,j)*temp*(0.5d0*yl+yf0-tempy)*xl*coeff2
	    enddo
    enddo

    ! couette shear in Y Direction
    do j = 1,ny+1
	    jp1=j+1
	    jm1 = j-1

        if(j.eq.1) then
            dy_area = (yref(2) - yref(1)) / 2.d0 
        elseif (j .eq. ny+1) then
            dy_area = (yref(ny) - yref(ny-1)) / 2.d0 
        else
            dy_area = yref(j) - yref(jm1)
        endif

	    do i = 1,nx
	        im1=i-1
	        ip1=i+1
            if (i .eq. 1) im1 = i
            if (i .eq. nx) ip1 = i        
	        dx_area = (xref(ip1) - xref(im1)) / 2.d0
	        
	        dxdy=dx_area*dy_area
            if (j .eq. 1) then
                height = h(i, 1) + hnew(i, 1)
                bear = beary(i, 1)
            elseif (j .eq. ny+1) then
                height = h(i, ny) + hnew(i, ny)
                bear = beary(i, ny)
            else
                height = recssj(i,j) + (hnew(i,j)+hnew(i,jm1)) / 2.d0
                bear = (beary(i,jm1) + beary(i, j)) / 2.d0
            endif 

            temp=dxdy/(height+2*corCoef*al)
            if (j .gt. 1 .and. j .lt. ny+1) fsr2=fsr2 + bear*temp
            
            CouetteShearY(i,j) = (9.81 * bear * temp * coeff2 / zl) / (dxdy * xl * xl * p0)
            
            tempx = xref(i)
		    Zmom = Zmom - beary(i,j)*temp*(0.5d0+xf0-tempx)*xl*coeff2
	    enddo
    enddo

    fsp1=fsp1*coeff1
    fsr1=fsr1*coeff1
    fsp2=fsp2*coeff2
    fsr2=fsr2*coeff2

    fsp=fsp1+fsp2
    fsr=fsr1+fsr2

    Zmom = Zmom*9.81d6/zl
    return
    end !subroutine shear
    

!
!!=============================================
!    subroutine shearorig
!!=============================================
!
!    use Q4_globals
!     use ShearArrays
!    implicit none
!    
!    real(8) :: fsp1, fsr1, fsp2, fsr2, coeff1, coeff2, dy, dx
!    real(8) :: tempfsp1, tempy, tempfsr1, tempx, dxdy, temp
!    integer :: i, ip1, im1, j, jp1, jm1
!
!    fsp1=0.d0
!    fsr1=0.d0
!    fsp2=0.d0
!    fsr2=0.d0
!
!    coeff1=zl*xl*p0*hm/9.81d0
!    coeff2=zl*vis1*xl*xl*xl/hm/9.81d0
!!
!!     fsp1 is first term in integral ((H/2)*(dp/dx))
!!     fsp2 is second term (Bear (actually just the velocity))/(slip stuff)
!!     Bearx and beary are actually not the bearing number but the velocity in x and y
!!     the other params in the bearing number are included in coeff1 and coeff2
!!
!    Zmom = 0.d0
!    nxm1=nx-1
!    nym1=ny-1
!    do i=1,nxm1
!	    ip1 = i+1
!	    im1 = i-1
!	    if(i.eq.1)im1=i
!	    dx=(xref(ip1)-xref(im1))/2.d0
!	    do j=2,nym1
!	        jm1 = j-1
!	        jp1 = j+1
!	        dy=(yref(jp1)-yref(jm1))/2.d0
!            ! Modified and added by Ling Huang for Z moment calculation.
!	        tempfsp1= -(p(ip1,j)-p(i,j)) * dy * (recssi(ip1,j)+(hnew(i,j)+hnew(ip1,j))/2.d0)/2.d0 
!	        fsp1=fsp1+tempfsp1
!	        
!	        !shear stress = (shear_force / area)
!	        !PoisilleShearX(i,j) = (tempfsp1 * coeff1 * (1000/zl)) / (dx*dy * xl * xl) !grams
!	        PoisilleShearX(i,j) = (9.81 * tempfsp1 * coeff1 / zl) / (dx*dy * xl * xl * p0) !Atmospheres
!	        
!		    tempy = (yref(jp1)+yref(jm1))/2.d0
!		    Zmom = Zmom + tempfsp1*(0.5d0*yl+yf0-tempy)*xl*coeff1
!	    enddo
!    enddo
!
!    do j=1,nym1
!	    jp1=j+1
!	    jm1 = j-1
!	    if(j.eq.1)jm1=i
!	    dy=(yref(jp1)-yref(jm1))/2.d0 
!	    do i=2,nxm1
!	        im1=i-1
!	        ip1=i+1
!	        dx=(xref(ip1)-xref(im1))/2.d0
!            ! Modified and added by Ling Huang for Z moment calculation.
!		    tempfsr1 = -(p(i,jp1)-p(i,j)) * dx * (recssj(i,jp1)+(hnew(i,j)+hnew(i,jp1))/2.d0)/2.d0
!	        fsr1=fsr1 + tempfsr1
!	        
!	        !PoisilleShearY(i,j) = (tempfsr1 * coeff1 * (1000/zl)) / (dx*dy * xl * xl) !grams
!	        PoisilleShearY(i,j) = (9.81 * tempfsr1 * coeff1 / zl) / (dx*dy * xl * xl * p0) !Atmospheres
!		    tempx = (xref(ip1)+xref(im1))/2.d0
!		    Zmom = Zmom - tempfsr1*(0.5d0+xf0-tempx)*xl*coeff1
!	    enddo
!    enddo
!
!    do i=1,nx
!	    ip1=i+1
!	    im1=i-1
!	    if(i.eq.1)im1=i
!	    if(i.eq.nx)ip1=i
!	    dx=(xref(ip1)-xref(im1))/2.d0
!	    do j=1,ny
!	        jp1=j+1
!	        jm1=j-1
!	        if(j.eq.1) jm1=j
!	        if(j.eq.ny) jp1=j
!	        dy=(yref(jp1)-yref(jm1))/2.d0
!	        dxdy=dx*dy
!	        temp=dxdy/(h(i,j)+hnew(i,j)+2*corCoef*al)
!	        fsp2=fsp2+bearx(i,j)*temp
!	        fsr2=fsr2+beary(i,j)*temp
!	        
!!            CouetteShearX(i,j) = (bearx(i,j) * temp * coeff2 * 1000/zl) / (dxdy * xl * xl)
!!            CouetteShearY(i,j) = (beary(i,j) * temp * coeff2 * 1000/zl) / (dxdy * xl * xl)
!            CouetteShearX(i,j) = (9.81 * bearx(i,j) * temp * coeff2 / zl) / (dxdy * xl * xl * p0)
!            CouetteShearY(i,j) = (9.81 * beary(i,j) * temp * coeff2 / zl) / (dxdy * xl * xl * p0)
!
!            ! Added by Ling Huang for Z moment calculation.
!		    tempy = (yref(jp1)+yref(jm1))/2.d0
!		    Zmom = Zmom + bearx(i,j)*temp*(0.5d0*yl+yf0-tempy)*xl*coeff2
!		    tempx = (xref(ip1)+xref(im1))/2.d0
!		    Zmom = Zmom - beary(i,j)*temp*(0.5d0+xf0-tempx)*xl*coeff2
!	    enddo
!    enddo
!
!    fsp1=fsp1*coeff1
!    fsr1=fsr1*coeff1
!    fsp2=fsp2*coeff2
!    fsr2=fsr2*coeff2
!
!    fsp=fsp1+fsp2
!    fsr=fsr1+fsr2
!
!    Zmom = Zmom*9.81d6/zl
!    return
!    end


!=============================================
    subroutine shearDuCorrect
!=============================================

    use Q4_globals
    implicit none

    real(8) :: fsp1, fsr1, fsp2, fsr2, coeff1, coeff2, dy, dx
    real(8) :: tempfsp1, tempy, tempfsr1, tempx, dxdy, temp
    integer :: i, ip1, im1, j, jp1, jm1

    fsp1=0.d0
    fsr1=0.d0
    fsp2=0.d0
    fsr2=0.d0

    coeff1=zl*xl*p0*hm/9.81d0
    coeff2=zl*vis1*xl*xl*xl/hm/9.81d0
!
!     fsp1 is first term in integral ((H/2)*(dp/dx))
!     fsp2 is second term (Bear (actually just the velocity))/(slip stuff)
!     Bearx and beary are actually not the bearing number but the velocity in x and y
!     the other params in the bearing number are included in coeff1 and coeff2
!
    Zmom = 0.d0
    nxm1=nx-1
    nym1=ny-1
    do i=1,nxm1
	    ip1=i+1
	    do j=2,nym1
	        jm1=j-1
	        jp1=j+1
	        dy=(yref(jp1)-yref(jm1))/2.d0 
    	    
            ! Modified and added by Ling Huang for Z moment calculation.
	        tempfsp1= -(p(ip1,j)-p(i,j))*dy * (recssi(ip1,j)+(hnew(i,j)+hnew(ip1,j))/2.d0)/2.d0 
	        fsp1=fsp1+tempfsp1
		    tempy = (yref(jp1)+yref(jm1))/2.d0
		    Zmom = Zmom + tempfsp1*(0.5d0*yl+yf0-tempy)*xl*coeff1
	    enddo
    enddo

    do j=1,nym1
	    jp1=j+1
	    do i=2,nxm1
	        im1=i-1
	        ip1=i+1
	        dx=(xref(ip1)-xref(im1))/2.d0
            ! Modified and added by Ling Huang for Z moment calculation.
		    tempfsr1 = -(p(i,jp1)-p(i,j))*dx*(recssj(i,jp1)+(hnew(i,j)+hnew(i,jp1))/2.d0)/2.d0
	        fsr1=fsr1 + tempfsr1
		    tempx = (xref(ip1)+xref(im1))/2.d0
		    Zmom = Zmom - tempfsr1*(0.5d0+xf0-tempx)*xl*coeff1
	    enddo
    enddo

    do i=1,nx
	    ip1=i+1
	    im1=i-1
	    if(i.eq.1) im1=i
	    if(i.eq.nx) ip1=i
	    dx=(xref(ip1)-xref(im1))/2.d0
	    do j=1,ny
	        jp1=j+1
	        jm1=j-1
	        if(j.eq.1) jm1=j
	        if(j.eq.ny) jp1=j
	        dy=(yref(jp1)-yref(jm1))/2.d0
	        dxdy=dx*dy
            ! shear1
            ! temp=dxdy/(h(i,j)+hnew(i,j)+2*al)
            ! the above is changed to the right(not local kn),shear 2
            temp=dxdy/(h(i,j)+hnew(i,j)+2*corCoef*al)
            
	        fsp2=fsp2+bearx(i,j)*temp
	        fsr2=fsr2+beary(i,j)*temp

            ! Added by Ling Huang for Z moment calculation.
		    tempy = (yref(jp1)+yref(jm1))/2.d0
		    Zmom = Zmom + bearx(i,j)*temp*(0.5d0*yl+yf0-tempy)*xl*coeff2
		    tempx = (xref(ip1)+xref(im1))/2.d0
		    Zmom = Zmom - beary(i,j)*temp*(0.5d0+xf0-tempx)*xl*coeff2
	    enddo
    enddo

    fsp1=fsp1*coeff1
    fsr1=fsr1*coeff1
    fsp2=fsp2*coeff2
    fsr2=fsr2*coeff2
    !  ****
    fsp=fsp1+fsp2
    fsr=fsr1+fsr2
    ! ****
    Zmom = Zmom*9.81d6/zl
    return
    end      
      
!=============================================
      subroutine shearKang
!=============================================
!
! Uses Soo-Choon Kang's Slip table to calculate the shear force
! See: "A new molecular gas lubrication theory suitable for head-disk interface modeling"
!      and also Kang's thesis from CMU under  Myung S. Jhon
!
    use Q4_globals
    implicit none
    
    integer, parameter :: FK_DATA_SIZE = 49
    integer :: i, j, ip1, im1, jp1, jm1
    real(8) :: fsp1, fsr1, fsp2, fsr2, coeff1, coeff2, zz, coep, temp
    real(8) :: tempfsp1, tempy, tempfsr1, tempx, dy, dx, dxdy, coec
	real(8) :: wc1(FK_DATA_SIZE),wc2(FK_DATA_SIZE),wc3(FK_DATA_SIZE),wc4(FK_DATA_SIZE), &
	           D(FK_DATA_SIZE),wp(FK_DATA_SIZE)
	
    data D/0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.05,0.06,0.07, &
     	    0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7, &
             0.8,0.9,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,9.0, &
             10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,70.0,80.0, &
             90.0,100.0/	         
	data wp/0.9998,0.9999,0.9999,0.9999,0.9999,0.9999,1.0000,1.0000, &
              1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000, &
              1.0000,1.0001,1.0001,1.0001,1.0001,1.0001,1.0001,1.0002, &
              1.0002,1.0001,1.0001,1.0002,1.0002,1.0003,1.0003,1.0005, &
              1.0006,1.0007,1.0009,1.0011,1.0013,1.0007,1.0011,1.0016, &
              1.0022,1.0027,1.0034,1.0048,1.0063,1.0079,1.0097,1.0115, &
              1.0134/
    data wc1/0.0056,0.0084,0.0111,0.0138,0.0165,0.0192,0.0218,0.0271, &
              0.0323,0.0374,0.0424,0.0473,0.0522,0.0758,0.0979,0.1187, &
              0.1385,0.1572,0.1751,0.2083,0.2387,0.2667,0.2926,0.3166, &
              0.3389,0.4313,0.5006,0.5548,0.5984,0.6342,0.6643,0.7117, &
              0.7474,0.7754,0.7978,0.8163,0.8317,0.8811,0.9085,0.9260, &
              0.9381,0.9472,0.9543,0.9648,0.9725,0.9786,0.9837,0.9881, &
              0.9920/

    data wc2/0.0051,0.0076,0.0101,0.0125,0.0150,0.0174,0.0198,0.0246, &
               0.0293,0.0340,0.0386,0.0431,0.0476,0.0692,0.0896,0.1090, &
               0.1274,0.1449,0.1616,0.1930,0.2218,0.2485,0.2733,0.2964, &
               0.3180,0.4083,0.4771,0.5314,0.5756,0.6121,0.6429,0.6920, &
               0.7293,0.7586,0.7822,0.8017,0.8180,0.8709,0.9003,0.9192, &
               0.9324,0.9422,0.9499,0.9613,0.9697,0.9762,0.9817,0.9864, &
               0.9907/

    data wc3/0.0046,0.0069,0.0091,0.0113,0.0136,0.0158,0.0180,0.0223, &
              0.0266,0.0309,0.0351,0.0392,0.0433,0.0631,0.0819,0.0999, &
              0.1170,0.1333,0.1490,0.1784,0.2057,0.2311,0.2548,0.2769, &
              0.2977,0.3856,0.4535,0.5078,0.5523,0.5895,0.6210,0.6715, &
              0.7102,0.7408,0.7657,0.7862,0.8035,0.8598,0.8915,0.9118, &
              0.9261,0.9367,0.9451,0.9575,0.9665,0.9736,0.9795,0.9846, &
              0.9891/
      
    data wc4/0.0037,0.0056,0.0074,0.0093,0.0111,0.0129,0.0147,0.0183, &
              0.0218,0.0253,0.0288,0.0323,0.0357,0.0523,0.0681,0.0833, &
              0.0980,0.1121,0.1256,0.1514,0.1755,0.1982,0.2195,0.2397, &
              0.2587,0.3408,0.4062,0.4597,0.5043,0.5422,0.5747,0.6277, &
              0.6690,0.7021,0.7293,0.7520,0.7713,0.8349,0.8713,0.8949, &
              0.9116,0.9241,0.9338,0.9484,0.9590,0.9672,0.9741,0.9800, &
              0.9852/      


    fsp1=0.d0
    fsr1=0.d0
    fsp2=0.d0
    fsr2=0.d0

    coeff1=zl*xl*p0*hm/9.81d0
    coeff2=zl*vis1*xl*xl*xl/hm/9.81d0

    Zmom=0.d0
    nxm1=nx-1
    nym1=ny-1
    do i=1,nxm1
	    ip1=i+1
	    do j=2,nym1
	        jm1=j-1
	        jp1=j+1
	        dy=(yref(jp1)-yref(jm1))/2.d0
            !  wp    
            zz=d0*(h(i,j)+hnew(i,j))
	        if((zz.gt.100).or.(zz.lt.0.01)) then
	            coep=1.0
	        else
                call finteg1(D,wp,49,zz,coep)
	        endif
	        tempfsp1=-coep*(p(ip1,j)-p(i,j))*dy*(recssi(ip1,j)+(hnew(i,j)+hnew(ip1,j))/2.d0)/2.d0 
	        fsp1=fsp1+tempfsp1
		    tempy = (yref(jp1)+yref(jm1))/2.d0
		    Zmom = Zmom + tempfsp1*(0.5d0*yl+yf0-tempy)*xl*coeff1
	    enddo
    enddo

    do j=1,nym1
	    jp1=j+1
	    do i=2,nxm1
	        im1=i-1
	        ip1=i+1
	        dx=(xref(ip1)-xref(im1))/2.d0
            !  wp    
            zz=d0*(h(i,j)+hnew(i,j))
            if((zz.gt.100).or.(zz.lt.0.01)) then
	            coep=1.0
	        else
                call finteg1(D,wp,49,zz,coep)
	        endif
		    tempfsr1=-coep*(p(i,jp1)-p(i,j))*dx*(recssj(i,jp1)+(hnew(i,j)+hnew(i,jp1))/2.d0)/2.d0
	        fsr1=fsr1+tempfsr1
		    tempx = (xref(ip1)+xref(im1))/2.d0
		    Zmom = Zmom - tempfsr1*(0.5d0+xf0-tempx)*xl*coeff1
	    enddo
    enddo

    do i=1,nx
	    ip1=i+1
	    im1=i-1
	    if(i.eq.1)im1=i
	    if(i.eq.nx)ip1=i
	    dx=(xref(ip1)-xref(im1))/2.d0
	    do j=1,ny
	        jp1=j+1
	        jm1=j-1
	        if(j.eq.1)jm1=j
	        if(j.eq.ny)jp1=j
	        dy=(yref(jp1)-yref(jm1))/2.d0
            !  wc
            zz=d0*(h(i,j)+hnew(i,j))
            if(zz.gt.100)then
	            coec=1.0
		    else if(zz.lt.0.01) then
	            coec=0.0
	        else
            if(accom.eq.1.00d0) call finteg1(D,wc1,49,zz,coec)
	        if(accom.eq.0.95d0) call finteg1(D,wc2,49,zz,coec)
	        if(accom.eq.0.90d0) call finteg1(D,wc3,49,zz,coec)
     	    if(accom.eq.0.80d0) call finteg1(D,wc4,49,zz,coec)
            endif
	        dxdy=dx*dy
	        temp=coec*dxdy/(h(i,j)+hnew(i,j))
	        fsp2=fsp2+bearx(i,j)*temp
	        fsr2=fsr2+beary(i,j)*temp
		    tempy = (yref(jp1)+yref(jm1))/2.d0
		    Zmom = Zmom + bearx(i,j)*temp*(0.5d0*yl+yf0-tempy)*xl*coeff2
		    tempx = (xref(ip1)+xref(im1))/2.d0
		    Zmom = Zmom - beary(i,j)*temp*(0.5d0+xf0-tempx)*xl*coeff2
	    enddo
    enddo

    fsp1=fsp1*coeff1
    fsr1=fsr1*coeff1
    fsp2=fsp2*coeff2
    fsr2=fsr2*coeff2

    fsp=fsp1+fsp2
    fsr=fsr1+fsr2

	Zmom = Zmom*9.81d6/zl
    return
    end subroutine shearKang
    
    
    
!==========================================
      subroutine ELECTROSTATIC
!==========================================      
!    Calculate ELECTROSTATIC FORCES. added 10/8/04
    use Q4_globals
    implicit none

    real(8) :: felecst, relecst, telecst, eps0, elecnormf, tempgap, elecst
    real(8) :: delx, dxdy, xloc, dely
    integer, parameter :: Ke = 1 ! dielectric constant of medium i.e. air
    integer :: i, im1, ip1, j, jm1, jp1

    felecst=0.d0
    relecst=0.d0
    telecst=0.d0     
      	
    ! permittivity constant
    eps0 = 8.85d-12 
    
    ! un-normalization factor
    elecnormf = xl*xl/9.81d0
	
    do j=2,nym1
        jm1=j-1
        jp1=j+1
        dely=(yref(jp1)-yref(jm1))/2.d0
        do i=2,nxm1
            im1=i-1
            ip1=i+1
            delx=(xref(ip1)-xref(im1))/2.d0
            dxdy=delx*dely

            tempgap = hm*(recssi(i,j) + hnew(i,j))
	        elecst=eps0*Ke*elecpot*elecpot*dxdy/2.d0/(tempgap**2) !Electrostatic
            	
            felecst=felecst-elecst
            	
            ! pitch torque from electrostatic potential
            telecst=telecst-elecst*(xref(i) - 0.5d0 - xf0)
            xloc = xloc - elecst * xref(i)
            ! roll torque from electrostatic potential
            relecst=relecst-elecst*(yref(j) - yl/2.d0 - yf0)
        enddo
    enddo

    xloc = xloc / felecst
      
    felecst = felecst * elecnormf  	  
	telecst = telecst * elecnormf
    relecst = relecst * elecnormf
	  
	xf = xf + telecst
	yf = yf + relecst
	f = f + felecst
	
	fneg = fneg + felecst 	
	  	  
    return
    end subroutine ELECTROSTATIC
    

!======================================
    subroutine getrcf(iwrite)
!======================================
    use Q4_globals
    common/var/hnt,dfctn
    
    integer, intent(in) :: iwrite
    
    real(8) :: atot, ac, fmotion, hcm, ssk, csk, hp, hyy, coe1
    real(8) :: dely, delx, dxdy, hnt
    real(8) :: aa, af, x1, y1, coe2, dfctn, hnte, aa1, aa2
    real(8) :: fit1, fc1, fitb1, fcb1, fit32, fc32
    integer :: i, im1, ip1, j, jp1, jm1
    
    external fit1,fitb1,fit32

    integer, parameter :: CONTACT_DATA_SIZE = 110
    real(8) :: tgs(CONTACT_DATA_SIZE), fgs1(CONTACT_DATA_SIZE), fgs2(CONTACT_DATA_SIZE)

    ! create the functions for the greenwood-williamson model

    data tgs/0.00d0,0.05d0,0.10d0,0.15d0,0.20d0,0.25d0,0.30d0, &
     	        0.35d0,0.40d0,0.45d0,0.50d0,0.55d0,0.60d0,0.65d0, &
     	        0.70d0,0.75d0,0.80d0,0.85d0,0.90d0,0.95d0,1.00d0, &
     	        1.05d0,1.10d0,1.15d0,1.20d0,1.25d0,1.30d0,1.35d0, &
     	        1.40d0,1.45d0,1.50d0,1.55d0,1.60d0,1.65d0,1.70d0, &
     	        1.75d0,1.80d0,1.85d0,1.90d0,1.95d0,2.00d0,2.05d0, &
     	        2.10d0,2.15d0,2.20d0,2.25d0,2.30d0,2.35d0,2.40d0, &
     	        2.45d0,2.50d0,2.55d0,2.60d0,2.65d0,2.70d0,2.75d0, &
     	        2.80d0,2.85d0,2.90d0,2.95d0,3.00d0,3.05d0,3.10d0, &
     	        3.15d0,3.20d0,3.25d0,3.30d0,3.35d0,3.40d0,3.45d0, &
     	        3.50d0,3.55d0,3.60d0,3.65d0,3.70d0,3.75d0,3.80d0, &
     	        3.85d0,3.90d0,3.95d0,4.00d0,4.05d0,4.10d0,4.15d0, &
     	        4.20d0,4.25d0,4.30d0,4.35d0,4.40d0,4.45d0,4.50d0, &
     	        4.55d0,4.60d0,4.65d0,4.70d0,4.75d0,4.80d0,4.85d0, &
     	        4.90d0,4.95d0,5.00d0,5.50d0,6.00d0,7.00d0,8.00d0, &
     	        9.00d0,10.0d0,20.0d0,500.d0,1000000000.d0/
      data fgs1/0.3989d0,0.3744d0,0.3509d0,0.3284d0,0.3069d0,0.2863d0, &
     		0.2668d0,0.2481d0,0.2304d0,0.2137d0,0.1978d0,0.1828d0, &
     		0.1687d0,0.1554d0,0.1429d0,0.1312d0,0.1202d0,0.1100d0, &
     		0.1004d0,0.9156d-1,0.8332d-1,0.7568d-1,0.6862d-1, &
     		0.6210d-1,0.5610d-1,0.5059d-1,0.4553d-1,0.4090d-1, &
     		0.3667d-1,0.3281d-1,0.2930d-1,0.2612d-1,0.2324d-1, &
     		0.2064d-1,0.1829d-1,0.1617d-1,0.1428d-1,0.1257d-1, &
     		0.1105d-1,0.9698d-2,0.8490d-2,0.7418d-2,0.6468d-2, &
     		0.5628d-2,0.4887d-2,0.4235d-2,0.3662d-2,0.3159d-2, &
     		0.2720d-2,0.2337d-2,0.2004d-2,0.1715d-2,0.1464d-2, &
     		0.1247d-2,0.1060d-2,0.8992d-3,0.7611d-3,0.6428d-3, &
     		0.5417d-3,0.4555d-3,0.3822d-3,0.3199d-3,0.2673d-3, &
     		0.2228d-3,0.1852d-3,0.1537d-3,0.1273d-3,0.1051d-3, &
     		0.8666d-4,0.7127d-4,0.5848d-4,0.4788d-4,0.3911d-4, &
     		0.3188d-4,0.2592d-4,0.2103d-4,0.1702d-4,0.1375d-4, &
     		0.1108d-4,0.8908d-5,0.7145d-5,0.5718d-5,0.4566d-5, &
     		0.3637d-5,0.2891d-5,0.2292d-5,0.1814d-5,0.1432d-5, &
     		0.1127d-5,0.8857d-6,0.6942d-6,0.5429d-6,0.4236d-6, &
     		0.3297d-6,0.2560d-6,0.1984d-6,0.1533d-6,0.1182d-6, &
     		0.9096d-7,0.6982d-7,0.5346d-7,0.3255d-8,0.1564d-9, &
     		0.1760d-12,0.7550d-16,0.1225d-19,0.7475d-24, &
     		0.1370d-89,0.0d0,0.0d0/
      data fgs2/0.4299d0,0.4000d0,0.3715d0,0.3446d0,0.3191d0,0.2952d0, &
     		0.2725d0,0.2513d0,0.2313d0,0.2127d0,0.1951d0,0.1789d0, &
     		0.1636d0,0.1495d0,0.1363d0,0.1241d0,0.1127d0,0.1023d0, &
     		0.9267d-1,0.8382d-1,0.7567d-1,0.6819d-1,0.6132d-1, &
     		0.5508d-1,0.4935d-1,0.4417d-1,0.3944d-1,0.3517d-1, &
     		0.3129d-1,0.2779d-1,0.2463d-1,0.2180d-1,0.1925d-1, &
     		0.1697d-1,0.1493d-1,0.1311d-1,0.1149d-1,0.1005d-1, &
     		0.8773d-2,0.7646d-2,0.6646d-2,0.5769d-2,0.4995d-2, &
     		0.4319d-2,0.3724d-2,0.3207d-2,0.2754d-2,0.2362d-2, &
     		0.2020d-2,0.1725d-2,0.1469d-2,0.1250d-2,0.1060d-2, &
     		0.8979d-3,0.7587d-3,0.6396d-3,0.5380d-3,0.4518d-3, &
     		0.3784d-3,0.3164d-3,0.2639d-3,0.2197d-3,0.1825d-3, &
     		0.1513d-3,0.1251d-3,0.1032d-3,0.8500d-4,0.6984d-4, &
     		0.5725d-4,0.4684d-4,0.3823d-4,0.3113d-4,0.2529d-4, &
     		0.2051d-4,0.1660d-4,0.1340d-4,0.1079d-4,0.8670d-5, &
     		0.6952d-5,0.5561d-5,0.4438d-5,0.3535d-5,0.2809d-5, &
     		0.2227d-5,0.1762d-5,0.1391d-5,0.1095d-5,0.8603d-6, &
     		0.6743d-6,0.5274d-6,0.4115d-6,0.3204d-6,0.2488d-6, &
     		0.1928d-6,0.1491d-6,0.1150d-6,0.8851d-7,0.6796d-7, &
     		0.5206d-7,0.3979d-7,0.3034d-7,0.1774d-8,0.8204d-10, &
     		0.8621d-13,0.3478d-16,0.5341d-20,0.3101d-24, &
     		0.4059d-90,0.0d0,0.0d0/

    ! calculate the contact force

    atot=0.d0
    ac  =0.d0
    fcr =0.d0
    txr =0.d0
    tyr =0.d0
    if(eyoung.eq.0.d0) return

    if(rsik.eq.0.d0.or.cta.eq.0.d0.or.rasp.eq.0.d0) return

    !if(omega.eq.0.d0) then
        !fmotion=0.d0
    !else
        fmotion=1.d0
    !endif

    if(icmod.eq.1) then
        !**********************************
        ! Greewood-Williamson contact model       
	    hcm=hm/rsik
	    ssk=dsin(ske)
	    csk=dcos(ske)
	    hp=-hx0*hm/xl
	    hyy=hy*hm/xl
	    coe1=4.0d0*dsqrt(rasp*rsik)*cta*eyoung*rsik/3.d0/9.81d0

        do j=1,ny
	        jm1=j-1
	        jp1=j+1
	        if(j.eq.1) then
	            dely=(yref(jp1)-yref(j))/2.d0
	        else if(j.eq.ny) then
	            dely=(yref(j)-yref(jm1))/2.d0
	        else
	            dely=(yref(jp1)-yref(jm1))/2.d0
	        endif
	        do i=1,nx
	            im1=i-1
	            ip1=i+1
	            if(i.eq.1) then
	                delx=(xref(ip1)-xref(i))/2.d0
	            else if(i.eq.nx) then
	                delx=(xref(i)-xref(im1))/2.d0
	            else
	                delx=(xref(ip1)-xref(im1))/2.d0
	            endif

	            dxdy=delx*dely
	            atot=atot+dxdy
	            hnt=hcm*(hnew(i,j)+h(i,j))
	            call finteg(tgs,fgs1,fgs2,110,hnt,aa,af)
	            aa=aa*dxdy
	            cp(i,j)=coe1*af
	            af=cp(i,j)*dxdy
	            ac=ac+aa
	            x1=xref(i)*xl
	            y1=yref(j)*xl

	            fcr=fcr+af
	            txr=txr+af*(xref(i)- 0.5d0 - xf0)*xl + frcoe*af*csk*zl
	            tyr=tyr+af*(yref(j)- 0.5d0*yl - yf0)*xl + frcoe*af*ssk*zl
	        end do
        end do
	    ac=twopi*cta*rasp*rsik*xl*xl*ac/2.d0

    else if(icmod.eq.2) then
        !******************************
        ! Elastic-plastic contact model
	    coe2=cta*eyoung*rsik/9.81d0
	    dfctn=(rasp*(twopi*ydcoe*ydst/eyoung/4.d0)**2)/rsik
	    hcm=hm/rsik
	    ssk=dsin(ske)
	    csk=dcos(ske)
	    hp=-hx0*hm/xl
	    hyy=hy*hm/xl

        do j=1,ny
	        jm1=j-1
	        jp1=j+1

	        if(j.eq.1) then
	            dely=(yref(jp1)-yref(j))/2.d0
	        else if(j.eq.ny) then
	            dely=(yref(j)-yref(jm1))/2.d0
	        else
	            dely=(yref(jp1)-yref(jm1))/2.d0
	        endif

	        do i=1,nx
	            im1=i-1
	            ip1=i+1

	            if(i.eq.1) then
	                delx=(xref(ip1)-xref(i))/2.d0
	            else if(i.eq.nx) then
	                delx=(xref(i)-xref(im1))/2.d0
	            else
	                delx=(xref(ip1)-xref(im1))/2.d0
	            endif

	            dxdy=delx*dely
	            atot=atot+dxdy
	            hnt=hcm*(hnew(i,j)+h(i,j))
	            hnte=hnt+dfctn

	            if(hnt.lt.4.d0) then
	                call qromb(fit1,hnt,hnte,fc1)
	                fc1=fc1/dsqrt(twopi)
	                call qromb(fitb1,hnte,6.d0,fcb1)
	                fcb1=fcb1/dsqrt(twopi)
	                call qromb(fit32,hnt,hnte,fc32)
	                fc32=fc32/dsqrt(twopi)
	                aa1=fc1+fcb1
	                aa2=4.d0*dsqrt(rasp*rsik)*fc32/3.d0+twopi*rasp*ydcoe*ydst*fcb1/eyoung/2.d0
	                ac=ac+twopi*cta*rasp*rsik*xl*xl*aa1*dxdy/2.d0
	                cp(i,j)=coe2*aa2
	                fcr=fcr+cp(i,j)*dxdy
	                x1=xref(i)*xl
	                y1=yref(j)*xl
	                txr=txr+cp(i,j)*((xref(i)- 0.5d0 - xf0)*xl+fmotion*zl*frcoe*csk)*dxdy
	                tyr=tyr+cp(i,j)*((yref(j)-(0.5d0*yl+yf0))*xl+fmotion*zl*frcoe*ssk)*dxdy
	            endif
	        end do
	    end do
    !*****************
	! no contact model
    else 
	    return
    endif
    
    atot=xl*xl*atot
    aratio=ac/atot

    ! convert to KG unit
    fcr=fcr*xl*xl
    txr=txr*xl*xl
    tyr=tyr*xl*xl

    !output the results
    if(iwrite.eq.1) then
	write(*,*)
	write(*,755) aratio
	write(*,756) fcr * 1000.d0
	write(*,759) txr * xl * 9.81d6
	write(*,760) tyr * xl * 9.81d6
    end if

755 format(' CONTACT AREA (NORM. BY APPARENT AREA) = ',g14.8)
756 format(' CONTACT FORCE  IN  GRAMS = ',g14.8)
759 format(' CONTACT P-MOMENT IN uN-M = ',g14.8)
760 format(' CONTACT R-TORQUE IN uN-M = ',g14.8)

    return
    end subroutine getrcf

!==========================================
	subroutine CoForce
!==========================================
! calculates the center positition of positive and negative force on the slider
	
    use Q4_globals
    implicit none
    
    real :: delx(nx-1), dely(ny-1)
    real :: area, force, xPosF, yPosF, xNegF, yNegF, negF, posF
    integer :: i, j

    xPosF = 0
    yPosF = 0
    xNegF = 0
    yNegF = 0
    negF = 0
    posF = 0

    do j = 2, ny-1
	    dely(j) = (yref(j+1) - yref(j-1)) / 2.0;
	    do i = 2, nx-1
	        delx(i) = (xref(i+1) - xref(i-1)) / 2.0;
		    area = delx(i) * dely(j)
		    force = (p(i,j) - 1.0) * area   !p(x,y)==1 is ambient pressure
		    if ( p(i,j) .ge. 1 ) then
		        XPosF = XPosF + (xref(i) * force)
		        yPosF = YPosF + (yref(j) * force)
		        posF = posF + force
		    else
		        xNegF = xNegF + (xref(i) * force)
		        yNegF = yNegF + (yref(j) * force)
   		        negF = negF + force
		    end if
	    enddo
    enddo
	
    xPosLoc = (xPosF/posF) * xl * 1000
    yPosLoc = (yPosF/posF) * xl * 1000
    if (negF .eq. 0.0) then
	    xNegLoc = 0.d0
	    yNegLoc = 0.d0
    else
	    xNegLoc = (xNegF/negF) * xl * 1000
	    yNegLoc = (yNegF/negF) * xl * 1000
    endif

	write(*,*)
	write(*,*) "CENTER POSITIONS OF FORCE"
	write(*, *) "pos X", xPosLoc, " (mm)"
	write(*, *) "neg X", xNegLoc, " (mm)"
	write(*, *) "pos Y", YPosLoc, " (mm)"
	write(*, *) "neg Y", YNegLoc, " (mm)"

	return
	end subroutine CoForce