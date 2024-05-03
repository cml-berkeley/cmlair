!========================================================
    subroutine ReynoldsDirected(xMin, xMax, yMin, yMax)
!========================================================
    use Q4_globals
    use PressureWarnings
    use NestingInfo
    implicit none    
    
    real(8) :: ak_std
    integer :: negp, i, j, ii, jj, im, ip, jm, jp
    integer :: xMin, xMax, yMin, yMax
    real(8) :: resArray(nx, ny)
    real(8) :: pold(nx, ny), su0(nx, ny)
    
    
    integer :: nx_m, ny_m
    integer :: boxSize
    real(8), dimension(:, :), allocatable :: bearx_m, beary_m, cohimx_m,        &
    cohjmx_m, h_m, himax_m, himin_m, hjmax_m, hjmin_m, hnew_m, p_m, pold_m,     &
    recssi_m, recssj_m, res_m, su0_m
    real(8), dimension(:), allocatable :: xref_m, yref_m
    real(8) :: ak0_m, ak_m, maximum
!    bearx_m
!    beary_m
!    cohimx_m
!    cohjmx_m
!    h_m
!    himax_m
!    himin_m
!    hjmax_m
!    hjmin_m
!    hnew_m
!    p_m
!    pold_m
!    recssi_m
!    recssj_m
!    res_m
!    su0_m
!    xref_m
!    yref_m
    
	do i=1,nx
	    do j=1,ny
	        resArray(i, j) = 0.0
	    enddo
	enddo
	
    !find maximum of AK
    !call GetReynoldsResidualAK(ak_std, bearx,beary,cohimx,cohjmx, &
    !             f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,recssi,recssj,resArray, su0,    &
    ! 	         xref,yref, idisc,nx,ny,nx,ny,0,GS_Iterations,negp)
     	         
    
!    call outputArray(p, nx, ny, 'pressOrig') 	         
!    do while (ak0 .gt. akmax)
!        call reyneq_nomax(akmax,1.d0,ak0,bearx,beary,cohimx,cohjmx,           &
!    		    f2p,h,himax,himin,hjmax,hjmin,hm,hnew,              &
!    		    p,pold,p0,recssi,recssj,res,su0,xl,xref,yref,1,     &
!    		    idisc,mitp(1),nx,ny,nx,ny,0,5,negp,1)
!    	write(*,*) ak0
!    enddo
!    call outputArray(p, nx, ny, 'pressNoMax')
!    
    !pause        


!    call outputArray(resArray, nx, ny, 'resArrayDirected')
!    
!    maximum = resArray(2, 2)
!    ii = 2
!    jj = 2
!    do j = 2, ny-1
!        do i = 2, nx-1
!            if (resArray(i, j) .gt. maximum) then
!                maximum = resArray(i, j)
!                ii = i
!                jj = j
!            endif
!        enddo
!    enddo
!    
!    !now expand outwards by 10 on each side or until we hit the edge of our domain
!    boxSize = 60
!    im = max(1, ii-boxSize)
!    ip = min(nx, ii+boxSize)
!    jm = max(1, jj-boxSize)
!    jp = min(ny, jj+boxSize)
!    

    im = xMin
    ip = xMax
    jm = yMin
    jp = yMax

    !now create parameters that will be used in the reynolds equation
    nx_m = ip - im + 1
    ny_m = jp - jm + 1



    allocate(bearx_m (nx_m, ny_m))
    allocate(beary_m (nx_m, ny_m))
    allocate(cohimx_m (nx_m, ny_m))
    allocate(cohjmx_m (nx_m, ny_m))
    allocate(h_m (nx_m, ny_m))
    allocate(himax_m (nx_m, ny_m))
    allocate(himin_m (nx_m, ny_m))
    allocate(hjmax_m (nx_m, ny_m))
    allocate(hjmin_m (nx_m, ny_m))
    allocate(hnew_m (nx_m, ny_m))
    allocate(p_m (nx_m, ny_m))
    allocate(pold_m (nx_m, ny_m))
    allocate(recssi_m (nx_m, ny_m))
    allocate(recssj_m (nx_m, ny_m))
    allocate(res_m (nx_m, ny_m))
    allocate(su0_m (nx_m, ny_m))
    allocate(xref_m (nx_m))
    allocate(yref_m (ny_m))

    bearx_m = bearx(im:ip, jm:jp)
    beary_m = beary(im:ip, jm:jp)
    cohimx_m = cohimx(im:ip, jm:jp)
    cohjmx_m = cohjmx(im:ip, jm:jp)
    h_m = h(im:ip, jm:jp)
    himax_m = himax(im:ip, jm:jp)
    himin_m = himin(im:ip, jm:jp)
    hjmax_m = hjmax(im:ip, jm:jp)
    hjmin_m = hjmin(im:ip, jm:jp)
    hnew_m = hnew(im:ip, jm:jp)
    p_m = p(im:ip, jm:jp)
    pold_m = pold(im:ip, jm:jp)
    recssi_m = recssi(im:ip, jm:jp)
    recssj_m = recssj(im:ip, jm:jp)
    res_m = res(im:ip, jm:jp)
    su0_m = su0(im:ip, jm:jp)
    xref_m = xref(im:ip)
    yref_m = yref(jm:jp)
    

    call outputArray(bearx_m, nx_m, ny_m, 'bearx')
    call outputArray(beary_m, nx_m, ny_m, 'beary')
    call outputArray(cohimx_m, nx_m, ny_m, 'cohimx')
    call outputArray(cohjmx_m, nx_m, ny_m, 'cohjmx')
    call outputArray(h_m, nx_m, ny_m, 'h')
    call outputArray(himax_m, nx_m, ny_m, 'himax')
    call outputArray(himin_m, nx_m, ny_m, 'himin')
    call outputArray(hjmax_m, nx_m, ny_m, 'hjmax')
    call outputArray(hjmin_m, nx_m, ny_m, 'hjmin')
    call outputArray(hnew_m, nx_m, ny_m, 'hnew')
    call outputArray(p_m, nx_m, ny_m, 'p')
    call outputArray(pold_m, nx_m, ny_m, 'pold')
    call outputArray(recssi_m, nx_m, ny_m, 'recssi')
    call outputArray(recssj_m, nx_m, ny_m, 'recssj')
    call outputArray(res_m, nx_m, ny_m, 'res')
    call outputArray(su0_m, nx_m, ny_m, 'su0')
    
    call openout(25, 'y_m.dat       ', 23*nmx)
    write(25,'(10000e23.15)')(yref_m(i),i=1,ny_m)
    close(25)   
    call openout(25, 'x_m.dat       ', 23*nmx)
    write(25,'(10000e23.15)')(xref_m(i),i=1,nx_m)
    close(25)
    
    ak0_m = ak0 * .01
    ak_m = ak
    
   
    call outputArray(p, nx, ny, 'pressIn')

    do i = 1, 100
        call reyneq(akmax, ak0_m, ak_m, bearx_m, beary_m, cohimx_m, cohjmx_m,             &
    		        f2p, h_m, himax_m, himin_m, hjmax_m, hjmin_m, hm, hnew_m,              &
    		        p_m, pold_m, p0, recssi_m, recssj_m, res_m, su0_m, xl, xref_m, yref_m, 1,     &
    		        idisc,mitp(1),nx_m,ny_m,nx_m,ny_m,0,10,negp,1)
        write(*,*) 'Directed AK: ', ak_m
        if (ak_m .lt. 1E-10) goto 111
    enddo

111 continue

    call outputArray(p_m, nx_m, ny_m, 'pressArray')
    
    p(im:ip, jm:jp) = p_m
    res(im:ip, jm:jp) = res_m
    
    call outputArray(p, nx, ny, 'NewPress')
    
    do i = 1, 10
        call reyneq(akmax,1.d0,ak0,bearx,beary,cohimx,cohjmx,           &
    		    f2p,h,himax,himin,hjmax,hjmin,hm,hnew,              &
    		    p,pold,p0,recssi,recssj,res,su0,xl,xref,yref,1,     &
    		    idisc,mitp(1),nx,ny,nx,ny,0,GS_Iterations,negp,1)
        write(*,*) 'Standard AK', ak0
    enddo    
    call outputArray(p, nx, ny, 'PressAgain')
    
    deallocate(bearx_m)
    deallocate(beary_m)
    deallocate(cohimx_m)
    deallocate(cohjmx_m)
    deallocate(h_m)
    deallocate(himax_m)
    deallocate(himin_m)
    deallocate(hjmax_m)
    deallocate(hjmin_m)
    deallocate(hnew_m)
    deallocate(p_m)
    deallocate(pold_m)
    deallocate(recssi_m)
    deallocate(recssj_m)
    deallocate(res_m)
    deallocate(su0_m)
    deallocate(xref_m)
    deallocate(yref_m)
    
    call outputStagnationData
    
    pause
    
    return
    end
    