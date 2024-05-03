!!    Control Codes:
!!   ============================
!!    ml - multigrid level
!!    0 == top level (could be level 5, 4, 3 2, or 1 if we're doing a full multigrid cycle at beginning of simulation) 
!!    negative means we're moving downwards to a coarser mesh
!!    -4 is level 5
!!    -3 is level 4
!!    -2 is level 3
!!    -1 is level 2
!!    positive means we're moving upwards
!!    3 = level 4
!!    2 = level 3
!!    1 = level 2
!!    ===================
!!   nsp = num line by line gauss-siedel sweeps
!

!====================================================================
      subroutine reyneq(akmax,akin,ak,bearx,beary,cohimx,cohjmx,    &
     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
     	   p0,recssi,recssj,res,su0,xl,xref,yref,           &
     	   idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp,curLevel)
!====================================================================
!     this is the reynolds equation solver
      !include 'size.fi'

    use Q4_sizes
    !use KrylovAccel
    use PressureWarnings
    
    implicit none
    
    real(8) :: akmax,akin,ak,bearx,beary,cohimx,cohjmx,     &
     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
     	   p0,recssi,recssj,res,su0,xl,xref,yref
    integer :: idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp, nrey, curLevel
    
    integer :: i, j, im1, ip1, jm1, jp1, nxm1, nxm2, nym1, nym2
    integer :: nsweep
    
    real(8) :: deltax, deltay, delxi, delxim1, delyj, delyjm1, resn
    real(8) :: dpi, dmi, dpj, dmj, fpi, fmi, fpj, fmj, sliderPlaneHeight
    real(8) :: pim0, him, pim, qn1, qn2, qni, pjm0, hjm, pjm, qnj
    real(8) :: ymyp, xmxp, qnim, qnip, qnjm, qnjp
    real(8) :: ppi, pmi, ppj, pmj, api, ami, apj, amj
    real(8) :: apitrm, amitrm, apjtrm, amjtrm
    real(8) :: aw, ae, as, an, ap, su, temp
    real(8) :: ai, bi, ci, di
    real(8) :: ak_std, ak_krylov
    integer :: L,ts,te,t
    
    !date and time info
    !==================
	integer :: start_time(8), end_time(8)
    CHARACTER (LEN = 12) REAL_CLOCK (3)
    
!    real(8) :: q(savedKrylovIterations, savedKrylovIterations), hMatrix(savedKrylovIterations, savedKrylovIterations),  &
!               eta(savedKrylovIterations), beta(savedKrylovIterations), alpha(savedKrylovIterations), alphaSum,   &
!               u_a(nx, ny), r_a(nx, ny), resNorm_a, resNorm
!    real(8) :: InnerProductMatrix
    
    !debugging info
    real(8) :: cur_residual
    
    
    dimension p(nx, ny)
    dimension bearx(ndx,ndy),beary(ndx,ndy),  &
     	pold(ndx,ndy),           &
     	qni(nx,ny),qnj(nx,ny),          &
     	h(ndx,ndy),hnew(ndx,ndy),           &
     	recssi(ndx,ndy), recssj(ndx,ndy),   &
     	xref(ndx), yref(ndy),               &
     	cohimx(ndx,ndy),cohjmx(ndx,ndy),    &
     	himax(ndx,ndy),himin(ndx,ndy),      &
     	hjmax(ndx,ndy),hjmin(ndx,ndy)       
    dimension ae(nx,ny),aw(nx,ny),as(nx,ny),an(nx,ny),    &
     	ap(nx,ny),su(nx,ny),res(ndx,ndy),       &
     	su0(ndx,ndy),ai(nmx),bi(nmx),ci(nmx),di(nmx)
    !dimension resSmoothRes(nx,ny), resSmoothPress(nx,ny)

    nxm1=nx-1
    nxm2=nx-2
    nym1=ny-1
    nym2=ny-2
    ak=1.d0
    nrey=0    

    write(5,*)'nx,ny ', nx,ny
    
    do while((.not.(nrey.gt.itrey .and. ak.lt.akin)) .and. nrey.le.30)
        nrey=nrey+1
        ak=0.d0
        ! tic
        call tic(ts)
        call GetReynoldsStencil(akmax,akin,ak,bearx,beary,cohimx,cohjmx,    &
     	                        f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
     	                        p0,recssi,recssj,res,su0,xl,xref,yref,           &
     	                        idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp,      &
                                aw, ae, as, an, ap, su, nrey)
        ! toc
        call tic(te)
        t = te - ts
        ! outputing time to reneq_tm.dat
        write(5,*)'nrey = ', nrey
        write(5,10)'ReneqStensil time = ',t
        
        !normalized residual
        ak=ak/(nx-2)/(ny-2)

        if (ak.lt.akmax.and.nrey.gt.1) goto 1002 ! return
        if (nrey.ge.(itrey+1) .and. ak.lt.akin) goto 1002 ! return

        nsweep=0
        resn=1.d0
   ! tic 
   call tic(ts)	
	    do while(nsweep.lt.nsp .and. resn.gt.akmax)
	        nsweep=nsweep+1
  
            !sweep in y direction
	        do j=2,nym1
	            jm1=j-1
	            jp1=j+1
	            do i=2,nxm1
	                ip1=i+1
	                im1=i-1
	                ai(im1)=-aw(i,j)
	                bi(im1)=ap(i,j)
	                ci(im1)=-ae(i,j)

	                if(1.eq.nxm2)then
		                di(im1) = su(i,j) + as(i,j)*p(i,jm1) + an(i,j)*p(i,jp1) + aw(i,j)*p(im1,j) + ae(i,j)*p(ip1,j)
	                else if(im1.eq.1)then
		                di(im1) = su(i,j) + as(i,j)*p(i,jm1) + an(i,j)*p(i,jp1) + aw(i,j)*p(im1,j)
	                else if(im1.eq.nxm2) then
		                di(im1) = su(i,j) + as(i,j)*p(i,jm1) + an(i,j)*p(i,jp1) + ae(i,j)*p(ip1,j)
	                else
		                di(im1) = su(i,j) + as(i,j)*p(i,jm1) + an(i,j)*p(i,jp1)
	                endif
	            enddo

	            call tridag(nxm2,ai,bi,ci,di)

	            do i=2,nxm1
	                p(i,j)=di(i-1)
	            enddo
	        enddo
	        
            !sweep in x direction
	        do i=2,nxm1
	            ip1=i+1
	            im1=i-1

	            do j=2,nym1
	                jm1=j-1
	                jp1=j+1
	                ai(jm1)=-as(i,j)
	                bi(jm1)=ap(i,j)
	                ci(jm1)=-an(i,j)
	                if (1.eq.nym2) then
		                di(jm1) = su(i,j) + aw(i,j)*p(im1,j) + ae(i,j)*p(ip1,j) + as(i,j)*p(i,jm1) + an(i,j)*p(i,jp1)
	                else if (jm1.eq.1) then
		                di(jm1) = su(i,j) + aw(i,j)*p(im1,j) + ae(i,j)*p(ip1,j) + as(i,j)*p(i,jm1)
	                else if (jm1.eq.nym2) then
		                di(jm1) = su(i,j) + aw(i,j)*p(im1,j) + ae(i,j)*p(ip1,j) + an(i,j)*p(i,jp1)
	                else
		                di(jm1) = su(i,j) + aw(i,j)*p(im1,j) + ae(i,j)*p(ip1,j)
	                endif
	            enddo

	            call tridag(nym2,ai,bi,ci,di)

	            do j=2,nym1
	                p(i,j)=di(j-1)
	            enddo
	        enddo

	        resn=0.d0

	        negp=0
 

            do i=2,nxm1
	            do j=2,nym1
                    !absolute negative pressure not allowed
	                if(p(i,j).le.0.d0)then
		                if(ml.eq.0)then
		                    p(i,j)=5.d-2
		                else
		                    negp=1.d0
		                endif
	                endif
	                
                    temp = ap(i,j) * p(i,j)
                    !sign error in cur_residual corrected 2/23/07
                    cur_residual = dabs((ae(i,j)*p(i+1,j)+aw(i,j)*p(i-1,j)+          &
            			        an(i,j)*p(i,j+1)+as(i,j)*p(i,j-1) - temp+su(i,j))/temp)

                    resn = resn + cur_residual

	            enddo
	        enddo

	        resn = resn /(nx-2) /(ny-2)
	        if(negp .eq. 1)then
	            write(*,*)'Warning negative pressure obtained.'
	            !numPressureWarnings = numPressureWarnings + 1
	            
	            if (ml .ne. 0) goto 1002 !return2

	        endif
	    enddo !number of sweeps
	    
	    ! toc
        call tic(te)
        t = te - ts
        !output time
        !write(5,*)
        write(5,10)'Sweep time = ',t
10	format(1X, A, I3.3)
    enddo  !fixed point iterations
    
1002 continue


    call GetReynoldsResidual(akmax,akin,ak,bearx,beary,cohimx,cohjmx,    &
     	                        f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
     	                        p0,recssi,recssj,res,su0,xl,xref,yref,           &
     	                        idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp,      &
                                aw, ae, as, an, ap, su, nrey)
    return
    end



!====================================================================
      subroutine GetReynoldsStencil(akmax,akin,ak,bearx,beary,cohimx,cohjmx,    &
     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,p0,recssi,recssj,res,   &
     	   su0,xl,xref,yref, idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp,       &
     	   aw, ae, as, an, ap, su, nrey)
!====================================================================
    use Q4_sizes
    use PressureWarnings
    implicit none
    
    real(8) :: akmax,akin,ak,bearx,beary,cohimx,cohjmx,     &
     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
     	   p0,recssi,recssj,res,su0,xl,xref,yref
    real(8) :: aw, ae, as, an, ap, su
    integer :: nrey
         	   
    integer :: idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp
    integer :: i, j, im1, ip1, jm1, jp1, nxm1, nxm2, nym1, nym2
    integer :: nsweep
    
    real(8) :: deltax, deltay, delxi, delxim1, delyj, delyjm1, resn
    real(8) :: dpi, dmi, dpj, dmj, fpi, fmi, fpj, fmj, sliderPlaneHeight
    real(8) :: pim0, him, pim, qn1, qn2, qni, pjm0, hjm, pjm, qnj
    real(8) :: ymyp, xmxp, qnim, qnip, qnjm, qnjp
    real(8) :: ppi, pmi, ppj, pmj, api, ami, apj, amj
    real(8) :: apitrm, amitrm, apjtrm, amjtrm, temp
    logical, parameter :: bCalculatingResiduals = .false.

    dimension p(nx, ny)
    dimension bearx(ndx,ndy),beary(ndx,ndy),  &
     	pold(ndx,ndy), qni(nx,ny),qnj(nx,ny), h(ndx,ndy), hnew(ndx,ndy),    &
     	recssi(ndx,ndy), recssj(ndx,ndy), xref(ndx), yref(ndy),             &
     	cohimx(ndx,ndy),cohjmx(ndx,ndy), himax(ndx,ndy),himin(ndx,ndy),     &
     	hjmax(ndx,ndy),hjmin(ndx,ndy)
    dimension ae(nx,ny),aw(nx,ny),as(nx,ny),an(nx,ny),    &
     	ap(nx,ny),su(nx,ny),res(ndx,ndy), su0(ndx,ndy)
     	
    
    call GetReynoldsVariables(akmax,akin,ak,bearx,beary,cohimx,cohjmx,    &
     	                        f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
     	                        p0,recssi,recssj,res,su0,xl,xref,yref,           &
     	                        idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp,      &
                                aw, ae, as, an, ap, su, nrey, bCalculatingResiduals)
     	
    return
    end subroutine GetReynoldsStencil



!====================================================================
      subroutine GetReynoldsResidual(akmax,akin,ak,bearx,beary,cohimx,cohjmx,    &
     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,p0,recssi,recssj,res,   &
     	   su0,xl,xref,yref, idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp,       &
     	   aw, ae, as, an, ap, su, nrey)
!====================================================================
    use Q4_sizes
    use PressureWarnings
    implicit none
    
    real(8) :: akmax,akin,ak,bearx,beary,cohimx,cohjmx,     &
     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
     	   p0,recssi,recssj,res,su0,xl,xref,yref
    real(8) :: aw, ae, as, an, ap, su
    integer :: nrey
         	   
    integer :: idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp
    integer :: i, j, im1, ip1, jm1, jp1, nxm1, nxm2, nym1, nym2
    integer :: nsweep
    
    real(8) :: deltax, deltay, delxi, delxim1, delyj, delyjm1, resn
    real(8) :: dpi, dmi, dpj, dmj, fpi, fmi, fpj, fmj, sliderPlaneHeight
    real(8) :: pim0, him, pim, qn1, qn2, qni, pjm0, hjm, pjm, qnj
    real(8) :: ymyp, xmxp, qnim, qnip, qnjm, qnjp
    real(8) :: ppi, pmi, ppj, pmj, api, ami, apj, amj
    real(8) :: apitrm, amitrm, apjtrm, amjtrm, temp
    logical, parameter :: bCalculatingResiduals = .true.
    
    real(8) :: aklocal
    
    dimension p(nx, ny)
    dimension bearx(ndx,ndy),beary(ndx,ndy),  &
     	pold(ndx,ndy), qni(nx,ny),qnj(nx,ny), h(ndx,ndy), hnew(ndx,ndy),    &
     	recssi(ndx,ndy), recssj(ndx,ndy), xref(ndx), yref(ndy),             &
     	cohimx(ndx,ndy),cohjmx(ndx,ndy), himax(ndx,ndy),himin(ndx,ndy),     &
     	hjmax(ndx,ndy),hjmin(ndx,ndy)
    dimension ae(nx,ny),aw(nx,ny),as(nx,ny),an(nx,ny),    &
     	ap(nx,ny),su(nx,ny),res(ndx,ndy), su0(ndx,ndy)
     	
     	
    call GetReynoldsVariables(akmax,akin,aklocal,bearx,beary,cohimx,cohjmx,    &
     	                        f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
     	                        p0,recssi,recssj,res,su0,xl,xref,yref,           &
     	                        idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp,      &
                                aw, ae, as, an, ap, su, nrey, bCalculatingResiduals)
     	
    return
    end subroutine GetReynoldsResidual




!====================================================================
      subroutine GetReynoldsVariables(akmax,akin,ak,bearx,beary,cohimx,cohjmx,    &
     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,p0,recssi,recssj,res,   &
     	   su0,xl,xref,yref, idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp,       &
     	   aw, ae, as, an, ap, su, nrey, bCalculatingResiduals)
!====================================================================
    use Q4_sizes
    use PressureWarnings
    implicit none
    
    real(8) :: akmax,akin,ak,bearx,beary,cohimx,cohjmx,     &
     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
     	   p0,recssi,recssj,res,su0,xl,xref,yref
    real(8) :: aw, ae, as, an, ap, su
    integer :: nrey
    logical :: bCalculatingResiduals
         	   
    integer :: idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp
    integer :: i, j, im1, ip1, jm1, jp1, nxm1, nxm2, nym1, nym2
    integer :: nsweep
    
    real(8) :: deltax, deltay, delxi, delxim1, delyj, delyjm1, resn
    real(8) :: dpi, dmi, dpj, dmj, fpi, fmi, fpj, fmj, sliderPlaneHeight
    real(8) :: pim0, him, pim, qn1, qn2, qni, pjm0, hjm, pjm, qnj
    real(8) :: ymyp, xmxp, qnim, qnip, qnjm, qnjp
    real(8) :: ppi, pmi, ppj, pmj, api, ami, apj, amj
    real(8) :: apitrm, amitrm, apjtrm, amjtrm, temp
    
    !debugging info
    real(8) :: cur_residual
    
    
    dimension p(nx, ny)
    dimension bearx(ndx,ndy),beary(ndx,ndy),  &
     	pold(ndx,ndy), qni(nx,ny),qnj(nx,ny), h(ndx,ndy), hnew(ndx,ndy),    &
     	recssi(ndx,ndy), recssj(ndx,ndy), xref(ndx), yref(ndy),             &
     	cohimx(ndx,ndy),cohjmx(ndx,ndy), himax(ndx,ndy),himin(ndx,ndy),     &
     	hjmax(ndx,ndy),hjmin(ndx,ndy)
    dimension ae(nx,ny),aw(nx,ny),as(nx,ny),an(nx,ny),    &
     	ap(nx,ny),su(nx,ny),res(ndx,ndy), su0(ndx,ndy)

    nxm1=nx-1
    nxm2=nx-2
    nym1=ny-1
    nym2=ny-2
    
	do i=2,nx
	    im1=i-1
	    ip1=i+1
	    do j=2,ny-1
	        jm1=j-1
	        jp1=j+1
	        pim0 = (p(i,j)+p(im1,j))/2.d0         !avg pressure along volume wall
	        sliderPlaneHeight = (hnew(i,j)+hnew(im1,j))/2.d0   !avg slider plane spacing along volume

            !slip at maximum recess along volume boundary
	        him=himax(i,j) + sliderPlaneHeight
	        pim = pim0 * him
	        call flow(pim,1.d0,1.d0,1.d0,qn1)
	        qn1 = qn1 * him**3
            
            !slip at minimum recess along volume boundary
	        him = dmax1(2d-10, himin(i,j) + sliderPlaneHeight)
	        pim = pim0 * him
	        call flow(pim,1.d0,1.d0,1.d0,qn2)
	        qn2 = qn2 * him**3
            
            !weighted average of slip at maximum and minimum recess heights
            !weighting is based on ratio of average to maximum and minimum 
            !height along volume boundary
	        qni(i,j)=qn1*cohimx(i,j)+qn2*(1.d0-cohimx(i,j))
	    enddo
	enddo

    do i=2,nx-1
	    im1=i-1
	    ip1=i+1
	    do j=2,ny
	        jm1=j-1
	        jp1=j+1
	        pjm0=(p(i,j)+p(i,jm1))/2.d0
	        sliderPlaneHeight=(hnew(i,j)+hnew(i,jm1))/2.d0

	        hjm=hjmax(i,j)+sliderPlaneHeight
	        pjm=pjm0*hjm
	        call flow(pjm,1.d0,1.d0,1.d0,qn1)
	        qn1=qn1*hjm**3

	        hjm= dmax1(2d-10, hjmin(i,j)+sliderPlaneHeight)
	        pjm=pjm0*hjm
	        call flow(pjm,1.d0,1.d0,1.d0,qn2)
	        qn2=qn2*hjm**3

	        qnj(i,j)=qn1*cohjmx(i,j)+qn2*(1.d0-cohjmx(i,j))
	    enddo
    enddo

    do j=2,nym1
	    jm1=j-1
	    jp1=j+1
	    delyj=yref(jp1)-yref(j)
	    delyjm1=yref(j)-yref(jm1)
	    ymyp=(yref(jp1)-yref(jm1))/2.d0
	    deltay=ymyp

	    do i=2,nxm1
	        im1=i-1
	        ip1=i+1
	        delxi=xref(ip1)-xref(i)
	        delxim1=xref(i)-xref(im1)
	        xmxp=(xref(ip1)-xref(im1))/2.d0
	        deltax =xmxp

	        qnim=qni(i,j)*f2p
	        qnip=qni(ip1,j)*f2p
	        qnjm=qnj(i,j)*f2p
	        qnjp=qnj(i,jp1)*f2p

	        qnim=qnim*(p(im1,j)+p(i,j))/2.d0
	        qnip=qnip*(p(ip1,j)+p(i,j))/2.d0
	        qnjm=qnjm*(p(i,jm1)+p(i,j))/2.d0
	        qnjp=qnjp*(p(i,jp1)+p(i,j))/2.d0

	        dpi=qnip/delxi*deltay
	        dmi=qnim/delxim1*deltay
	        dpj=qnjp/delyj*deltax
	        dmj=qnjm/delyjm1*deltax

            !recessX vars hold average recess at volume boundary
	        fpi=(bearx(ip1,j)+bearx(i,j))/4.d0*(recssi(ip1,j)+  &
                (hnew(i,j)+hnew(ip1,j))/2.d0)*deltay

	        fmi=(bearx(im1,j)+bearx(i,j))/4.d0*(recssi(i,j) +   &
     	        (hnew(i,j)+hnew(im1,j))/2.d0)*deltay

	        fpj=(beary(i,jp1)+beary(i,j))/4.d0*(recssj(i,jp1) + &
	            (hnew(i,j) + hnew(i,jp1))/2.d0)*deltax

	        fmj=(beary(i,jm1)+beary(i,j))/4.d0*(recssj(i,j) +   &
                (hnew(i,j)+hnew(i,jm1))/2.d0)*deltax

	        ppi=fpi/dpi
	        pmi=fmi/dmi
	        ppj=fpj/dpj
	        pmj=fmj/dmj

            !idisc = 3
	        if(idisc.eq.0)then
                !upwind
		        api=1.d0
		        ami=1.d0
		        apj=1.d0
		        amj=1.d0
	        else if(idisc.eq.1) then
                !hybrid
		        api=dmax1(0.d0,1.d0-0.5d0*dabs(ppi))
		        ami=dmax1(0.d0,1.d0-0.5d0*dabs(pmi))
		        apj=dmax1(0.d0,1.d0-0.5d0*dabs(ppj))
		        amj=dmax1(0.d0,1.d0-0.5d0*dabs(pmj))
		    else if (idisc .eq. 3) then
		        !no upwinding
		        api=1.d0 - 0.5d0*dabs(ppi)
		        ami=1.d0 - 0.5d0*dabs(pmi)
		        apj=1.d0 - 0.5d0*dabs(ppj)
		        amj=1.d0 - 0.5d0*dabs(pmj)
	        else
                !power law, idisc=2 (default)
		        apitrm=(1.d0-0.1d0*dabs(ppi))**5
		        amitrm=(1.d0-0.1d0*dabs(pmi))**5
		        apjtrm=(1.d0-0.1d0*dabs(ppj))**5
		        amjtrm=(1.d0-0.1d0*dabs(pmj))**5

		        api=dmax1(0.d0,apitrm)
		        ami=dmax1(0.d0,amitrm)
		        apj=dmax1(0.d0,apjtrm)
		        amj=dmax1(0.d0,amjtrm)
	        endif

	        aw(i,j)=dmi*ami+dmax1(fmi,0.d0)
	        ae(i,j)=dpi*api+dmax1(-fpi,0.d0)
	        as(i,j)=dmj*amj+dmax1(fmj,0.d0)
	        an(i,j)=dpj*apj+dmax1(-fpj,0.d0)
	        ap(i,j)=ae(i,j)+aw(i,j)+as(i,j)+an(i,j) + dmax1(fpi-fmi+fpj-fmj,0.d0)
	        su(i,j)=dmax1(fmi-fpi+fmj-fpj,0.d0)*p(i,j)


            !multigrid restriction source term
            !coarse grid MG method is to solve A(u) = f + A(v) + r
            !where u is exact pressure and v is restricted pressure from finer grid
            !if we're moving from a finer to coarser coarser grid and it's the first of our fixed point iterations
            !then add in our residual + solution approx from the finer grid pressure approximation
            !remember to multiply residual times area because we divided it by area on the finer level
	        if(ml.lt.0.and.nrey.eq.1 .and. bCalculatingResiduals .eq. .false.) then
	            su0(i,j) = res(i,j)*ymyp*xmxp                           &  
     	                    - aw(i,j)*p(im1,j) - ae(i,j)*p(ip1,j)       &
     	                    - as(i,j)*p(i,jm1) - an(i,j)*p(i,jp1) + ap(i,j)*p(i,j)-su(i,j)
            endif
            
	        if(ml.ne.0) su(i,j) = su(i,j) + su0(i,j)
	        !if(ml.lt.0) su(i,j) = su(i,j) + su0(i,j)

            !calculate residuals here
            !we should probably calculate this at the end with updated coefs

            temp = ap(i,j) * p(i,j)
	        res(i,j)=(aw(i,j)*p(im1,j) + ae(i,j)*p(ip1,j)           &
        		    + as(i,j)*p(i,jm1) + an(i,j)*p(i,jp1) - temp+su(i,j))


	        ak = ak + dabs(res(i,j)/temp)
            !output for the residual restriction
	        res(i,j)=res(i,j)/(ymyp*xmxp)
	    enddo
	enddo

    return
    end
    
    
!===================================================     
    subroutine calcCenterOfResidual (resArray, xref, yref, sizeX, sizeY)
!=================================================== 
    implicit none
    
    real(8) :: resArray, xref, yref
    integer :: sizeX, sizeY, i, j
    dimension :: resArray(sizeX, sizeY)
    dimension :: xref(sizeX), yref(sizeY)
    real(8) :: delx, dely
    real(8) :: area, force, xPosF, yPosF, posF, xposLoc, yposLoc
    real(8), parameter :: xl = 1.24
    
    xPosF = 0.0
    yPosF = 0.0
    posF = 0.0

    do j = 2, sizeY-1
	    dely = (yref(j+1) - yref(j-1)) / 2.0
	    do i = 2, sizeX-1
	        delx = (xref(i+1) - xref(i-1)) / 2.0
		    area = delx * dely
		    force = (resArray(i,j) - 1.0) * area   !p(x,y)==1 is ambient pressure
		    XPosF = XPosF + (xref(i) * force)
		    yPosF = YPosF + (yref(j) * force)
		    posF = posF + force
	    enddo
    enddo
	
    xPosLoc = (xPosF/posF) * xl
    yPosLoc = (yPosF/posF) * xl

	write(*, *) "pos X", xPosLoc, " (mm)"
	write(*, *) "pos Y", YPosLoc, " (mm)"   
    
    return
    end    


!     Multigrid vcycle updated 2/18/05 by BCox
!     Now allows for greater control over nesting
!     Still uses an abundance of goto statements but
!     in the end that's the easiest way to program it
!==========================================================
    subroutine vcycle(mlv,mnestIn)
!==========================================================

    use NestingInfo
    use Q4_globals
    implicit none

    real(8) :: ptemp1(nx1,ny1),ptemp2(nx2,ny2),ptemp3(nx3,ny3),ptemp4(nx4,ny4), &
               dp(nx,ny),dp1(nx1,ny1),dp2(nx2,ny2),dp3(nx3,ny3),dp4(nx4,ny4)
    real(8) :: pold, pold1, pold2, pold3, pold4 !, su0  !su0 isn't used on finest grid so there's no need to store array
    real(8) :: su0(nx, ny)
    integer :: mlv, mnest, ml, negp, mnestIn


    integer LevelsSeen
    integer extraGSIterations
    
    LevelsSeen = 0
    extraGSIterations = 0
    
    ! if we're on the finest level then vary the nesting a bit
    ! every other iteration go down 1 level less
    ! VCycleIteration
    !
    !
    mnest = mnestIn


    if(mlv.eq.2) goto 10
    if(mlv.eq.3) goto 20
    if(mlv.eq.4) goto 30

    !call restrict(p,p1,nx,ny,nx,ny,nx1,ny1)
    call restrictPressure(p, p1, xref, yref, xref1, yref1, &
                                nx, ny, nx1, ny1, enclosingFineRectX1, enclosingFineRectY1)
    call equate(nx1,ny1,nxm(1),nym(1),ptemp1,p1)
    call restrictPressure(res, res1, xref, yref, xref1, yref1, &
                                nx, ny, nx1, ny1, enclosingFineRectX1, enclosingFineRectY1)
!    call restrictResidual(res, res1, xref, yref, xref1, yref1, &
!                                nx, ny, nx1, ny1, enclosingCoarseRectX, enclosingCoarseRectY)
    !call restrict(res,res1,nx,ny,nx,ny,nx1,ny1)
      
    !write(*,*)'********grid level 2***********'
    
    !if (LevelsSeen+1 .ge. mnest) extraGSIterations = 7
    
    call reyneq(akmax,1.d0,ak1,bearx1,beary1,cohimx1,cohjmx1,           &
     		    f2p,h1,himax1,                                          &
     		    himin1,hjmax1,hjmin1,hm,hnew1,p1,pold1,p0,              &
     		    recssi1,recssj1,res1,su01,xl,xref1,yref1,1,             &
     		    idisc,mitr(1),nxm(1),nym(1),nx1,ny1,-1,GS_Iterations + extraGSIterations,negp,2)

    LevelsSeen = LevelsSeen + 1
    if(negp.eq.1)goto 75
    if(ak1.lt.akmax .or. LevelsSeen.ge.mnest) then
        extraGSIterations = 0
        goto 70
    endif
!	if(mnest.eq.1.or.negp.eq.1)goto 70

10  continue

    !call restrict(p1,p2,nxm(1),nym(1),nx1,ny1,nx2,ny2)
    call restrictPressure(p1, p2, xref1, yref1, xref2, yref2, &
                                nx1, ny1, nx2, ny2, enclosingFineRectX2, enclosingFineRectY2)
    call equate(nx2,ny2,nxm(2),nym(2),ptemp2,p2)
    call restrictPressure(res1, res2, xref1, yref1, xref2, yref2, &
                                nx1, ny1, nx2, ny2, enclosingFineRectX2, enclosingFineRectY2)
!    call restrictResidual(res1, res2, xref1, yref1, xref2, yref2, &
!                                nx1, ny1, nx2, ny2, enclosingCoarseRectX1, enclosingCoarseRectY1)
    !call restrict(res1,res2,nxm(1),nym(1),nx1,ny1,nx2,ny2)

    !write(*,*)'********grid level 3***********'
    
    !if (LevelsSeen+1 .ge. mnest) extraGSIterations = 7
    
    call reyneq(akmax,1.d0,ak2,bearx2,beary2,cohimx2,cohjmx2,         &
     		f2p,h2,himax2,                                            &
     		himin2,hjmax2,hjmin2,hm,hnew2,p2,pold2,p0,                &
     		recssi2,recssj2,res2,su02,xl,xref2,yref2,1,               &
     		idisc,mitr(2),nxm(2),nym(2),nx2,ny2,-2,GS_Iterations + extraGSIterations,negp,3)

    LevelsSeen = LevelsSeen + 1
    if(negp.eq.1)goto 65
    if(ak2.lt.akmax .or.  LevelsSeen.ge.mnest) then
        extraGSIterations = 0
        goto 60
    endif
!	if(mnest.eq.2.or.negp.eq.1)goto 60


20  continue

    !call restrict(p2,p3,nxm(2),nym(2),nx2,ny2,nx3,ny3)
    call restrictPressure(p2, p3, xref2, yref2, xref3, yref3, &
                                nx2, ny2, nx3, ny3, enclosingFineRectX3, enclosingFineRectY3)
    call equate(nx3,ny3,nxm(3),nym(3),ptemp3,p3)
    call restrictPressure(res2, res3, xref2, yref2, xref3, yref3, &
                                nx2, ny2, nx3, ny3, enclosingFineRectX3, enclosingFineRectY3)
!    call restrictResidual(res2, res3, xref2, yref2, xref3, yref3, &
!                                nx2, ny2, nx3, ny3, enclosingCoarseRectX2, enclosingCoarseRectY2)
    !call restrict(res2,res3,nxm(2),nym(2),nx2,ny2,nx3,ny3)

    !write(*,*)'********grid level 4***********'
    
    !if (LevelsSeen+1 .ge. mnest) extraGSIterations = 7
    
    call reyneq(akmax,1.d0,ak3,bearx3,beary3,cohimx3,cohjmx3,       &
     		    f2p,h3,himax3,                                      &
     		    himin3,hjmax3,hjmin3,hm,hnew3,p3,pold3,p0,          &
     		    recssi3,recssj3,res3,su03,xl,xref3,yref3,1,         &
     		    idisc,mitr(3),nxm(3),nym(3),nx3,ny3,-3,GS_Iterations + extraGSIterations,negp,4)

    LevelsSeen = LevelsSeen + 1
    if(negp.eq.1)goto 55
    if(ak3.lt.akmax .or. LevelsSeen.ge.mnest) then
        extraGSIterations = 0
        goto 50
    endif
    !if(mnest.eq.3.or.negp.eq.1)goto 50

30  continue

    !call restrict(p3,p4,nxm(3),nym(3),nx3,ny3,nx4,ny4)
    call restrictPressure(p3, p4, xref3, yref3, xref4, yref4, &
                                nx3, ny3, nx4, ny4, enclosingFineRectX4, enclosingFineRectY4)
    call equate(nx4,ny4,nxm(4),nym(4),ptemp4,p4)
    call restrictPressure(res3, res4, xref3, yref3, xref4, yref4, &
                                nx3, ny3, nx4, ny4, enclosingFineRectX4, enclosingFineRectY4)
!    call restrictResidual(res3, res4, xref3, yref3, xref4, yref4, &
!                                nx3, ny3, nx4, ny4, enclosingCoarseRectX3, enclosingCoarseRectY3)
    !call restrict(res3,res4,nxm(3),nym(3),nx3,ny3,nx4,ny4)

    !extraGSIterations = 7
    !write(*,*)'********grid level 5***********'
    LevelsSeen = LevelsSeen + 1
    call reyneq(akmax,1.d0,ak4,bearx4,beary4,cohimx4,cohjmx4,         &
    		    f2p,h4,himax4,                                        &
    		    himin4,hjmax4,hjmin4,hm,hnew4,p4,pold4,p0,            &
    		    recssi4,recssj4,res4,su04,xl,xref4,yref4,1,           &
    		    idisc,mitr(4),nxm(4),nym(4),nx4,ny4,-4,GS_Iterations + extraGSIterations,negp,5)

    extraGSIterations = 0
    if(negp.eq.1)goto 45

    if(mlv.eq.4)then
	    ml=0
    else
	    ml=3
    endif

    call matrixsub(dp4,p4,ptemp4,nxm(4),nym(4),nx4,ny4)
    !call prolong(dp4,dp3,nxm(4),nym(4),nx3,ny3,nx4,ny4)
    call prolongCorrection(dp3, dp4, xref3, yref3, xref4, yref4, &
                                 nx3, ny3, nx4, ny4, enclosingCoarseRectX3, enclosingCoarseRectY3)
    call matrixinc(p3,dp3,nxm(3),nym(3),nx3,ny3,ml,negp)

    if(negp.eq.1) goto 49

45  continue

    if(mlv.eq.4)then
	    ml=0
    else
	    ml=3
    endif

    !write(*,*)'********grid level 4***********'
    call reyneq(akmax,ak3,ak,bearx3,beary3,cohimx3,cohjmx3,         &
     		    f2p,h3,himax3,                                      &
     		    himin3,hjmax3,hjmin3,hm,hnew3,p3,pold3,p0,          &
     		    recssi3,recssj3,res3,su03,xl,xref3,yref3,1,         &
     		    idisc,mitp(4),nxm(3),nym(3),nx3,ny3,ml,GS_Iterations,negp,4)

    ak3=ak

49  continue

    if(mlv.eq.4) return
    if(negp.eq.1)goto 55

50  continue

    if(mlv.eq.3)then
	    ml=0
    else
	    ml=2
    endif

    call matrixsub(dp3,p3,ptemp3,nxm(3),nym(3),nx3,ny3)
    !call prolong(dp3,dp2,nxm(3),nym(3),nx2,ny2,nx3,ny3)
    call prolongCorrection(dp2, dp3, xref2, yref2, xref3, yref3, &
                                 nx2, ny2, nx3, ny3, enclosingCoarseRectX2, enclosingCoarseRectY2)
    call matrixinc(p2,dp2,nxm(2),nym(2),nx2,ny2,ml,negp)
    
    if(negp.eq.1)goto 59
    
55  continue

    if(mlv.eq.3)then
	    ml=0
    else
	    ml=2
    endif

    !write(*,*)'********grid level 3***********'
    call reyneq(akmax,ak2,ak,bearx2,beary2,cohimx2,cohjmx2,     &
     		f2p,h2,himax2,                                      &
     		himin2,hjmax2,hjmin2,hm,hnew2,p2,pold2,p0,          &
     		recssi2,recssj2,res2,su02,xl,xref2,yref2,1,         &
     		idisc,mitp(3),nxm(2),nym(2),nx2,ny2,ml,GS_Iterations,negp,3)

    ak2=ak
    
59  continue

    if(mlv.eq.3) return
    if(negp.eq.1)goto 65
    
60  continue

    if(mlv.eq.2)then
	    ml=0
    else
	    ml=1
    endif

    call matrixsub(dp2,p2,ptemp2,nxm(2),nym(2),nx2,ny2)
    !call prolong(dp2,dp1,nxm(2),nym(2),nx1,ny1,nx2,ny2)
    call prolongCorrection(dp1, dp2, xref1, yref1, xref2, yref2, &
                                 nx1, ny1, nx2, ny2, enclosingCoarseRectX1, enclosingCoarseRectY1)
    call matrixinc(p1,dp1,nxm(1),nym(1),nx1,ny1,ml,negp)

    if(negp.eq.1) goto 69
65  continue

    if(mlv.eq.2)then
	    ml=0
    else
	    ml=1
    endif

    !write(*,*)'********grid level 2***********'
    call reyneq(akmax,ak1,ak,bearx1,beary1,cohimx1,cohjmx1,           &
    		    f2p,h1,himax1,                                        &
    		    himin1,hjmax1,hjmin1,hm,hnew1,p1,pold1,p0,            &
    		    recssi1,recssj1,res1,su01,xl,xref1,yref1,1,           &
    		    idisc,mitp(2),nxm(1),nym(1),nx1,ny1,ml,GS_Iterations,negp,2)

    ak1=ak
69  continue

    if(mlv.eq.2) return
    if(negp.eq.1)goto 75
70  continue

    call matrixsub(dp1,p1,ptemp1,nxm(1),nym(1),nx1,ny1)
    !call prolong(dp1,dp,nxm(1),nym(1),nx,ny, nx1,ny1)
    call prolongCorrection(dp, dp1, xref, yref, xref1, yref1, &
                                 nx, ny, nx1, ny1, enclosingCoarseRectX, enclosingCoarseRectY)
    call matrixinc(p,dp,nx,ny,nx,ny,0,negp)

75  continue

    !write(*,*)'********grid level 1***********'
    call reyneq(akmax,ak0,ak,bearx,beary,cohimx,cohjmx,             &
    		    f2p,h,himax,himin,hjmax,hjmin,hm,hnew,              &
    		    p,pold,p0,recssi,recssj,res,su0,xl,xref,yref,1,     &
    		    idisc,mitp(1),nx,ny,nx,ny,0,GS_Iterations,negp,1)

    ak0=ak

    return
    end
    

!=========================================================
    subroutine fullmult(mstart,mnest,iwrite)
!=========================================================

    use Q4_globals
    use PressureWarnings
    use NestingInfo
    use TemperatureHumidity
    implicit none

    real(8) :: pold, pold1, pold2, pold3, pold4, su0
    integer :: mstart,mnest,iwrite, negp
    integer :: icount
    integer, parameter :: MAX_VCYCLE_ITR = 80
    real(8) :: percentChange, previousAK(MAX_VCYCLE_ITR+20)
    logical :: takenFirst, wantWCycle
    integer :: stagnationCount, numStagnations
    real(8) :: wAK

    
    call gethx
    !check to make sure we're not in contact in the disk
    do while(crash)
	    if(isolv.eq.0) then
	        write(*,*)'CONTACT DETECTED.'
	        stop
	    endif
	    write(*,*)'CONTACT DETECTED, ADJUSTING FLY HEIGHT...'
	    hmin = hmin + 0.2d0
	    h0 = hmin - hx0
	    call gethx
    enddo

    if(mstart.eq.0)goto 65

    write(*,*)'starting grid level 5'
    ak=1.d0

    icount=0

    do while(ak.gt.akmax .and. icount.lt.MAX_VCYCLE_ITR)
	    icount=icount+1

	    call reyneq(akmax,1.d0,ak,bearx4,beary4,cohimx4,cohjmx4,            &
        		    f2p,h4,himax4,                                          &
        		    himin4,hjmax4,hjmin4,hm,hnew4,p4,pold4,p0,              &
        		    recssi4,recssj4,res4,su04,xl,xref4,yref4,1,             &
        		    idisc,mitp(4),nxm(4),nym(4),nx4,ny4,0,GS_Iterations,negp,5)
	    write(*,*)'Normalized residual = ',ak
    enddo

    !call prolong(p4,p3,nxm(4),nym(4),nx3,ny3,nx4,ny4)
    call prolongCorrection(p3, p4, xref3, yref3, xref4, yref4, &
                           nx3, ny3, nx4, ny4, enclosingCoarseRectX3, enclosingCoarseRectY3)

5   continue
    write(*,*)'starting grid level 4'
6   continue
    call reyneq(akmax,1.d0,ak3,bearx3,beary3,cohimx3,cohjmx3,               &
    		  f2p,h3,himax3,                                                &
    		  himin3,hjmax3,hjmin3,hm,hnew3,p3,pold3,p0,                    &
    		  recssi3,recssj3,res3,su03,xl,xref3,yref3,1,                   &
    		  idisc,mitp(4),nxm(3),nym(3),nx3,ny3,0,GS_Iterations,negp,4)

    if(mnest.eq.0) then
	    write(*,*)'Normalized residual = ',ak3
	    if(ak3.gt.akmax) then
	        goto 6
	    else
	        goto 10
	    endif
    endif
    ak=1.d0

    icount=0
    do while(ak.gt.akmax .and. icount.lt.MAX_VCYCLE_ITR)
	    icount=icount+1
	    call vcycle(4,mnest)
	    write(*,*)'Normalized residual = ',ak
    enddo

10  continue

    !call prolong(p3,p2,nxm(3),nym(3),nx2,ny2,nx3,ny3)
    call prolongCorrection(p2, p3, xref2, yref2, xref3, yref3, &
                           nx2, ny2, nx3, ny3, enclosingCoarseRectX2, enclosingCoarseRectY2)
25  continue
    write(*,*)'starting grid level 3'
26  continue
    call reyneq(akmax,1.d0,ak2,bearx2,beary2,cohimx2,cohjmx2,       &
    		  f2p,h2,himax2,                                        &
    		  himin2,hjmax2,hjmin2,hm,hnew2,p2,pold2,p0,            &
    		  recssi2,recssj2,res2,su02,xl,xref2,yref2,1,           &
    		  idisc,mitp(3),nxm(2),nym(2),nx2,ny2,0,GS_Iterations,negp,3)

    if(mnest.eq.0) then
	    write(*,*)'Normalized residual = ',ak2
	    if(ak2.gt.akmax) then
	        goto 26
	    else
	        goto 30
	    endif
    endif

    ak=1.d0
    icount=0
    numPressureWarnings = 0
    stagnationCount = 0
    
    mnest = max_mg_nest
    do while(ak.gt.akmax .and. icount.lt.MAX_VCYCLE_ITR)
	    icount=icount+1
	    call vcycle(3,mnest)
	    write(*,*)'Normalized residual = ',ak
	    previousAK(icount) = ak
        call NestingCheck(previousAK, numPressureWarnings, stagnationCount, numStagnations, iCount, 0)
	    
    enddo

30  continue

    !call prolong(p2,p1,nxm(2),nym(2),nx1,ny1,nx2,ny2)
    call prolongCorrection(p1, p2, xref1, yref1, xref2, yref2, &
                           nx1, ny1, nx2, ny2, enclosingCoarseRectX1, enclosingCoarseRectY1)
45  continue
    write(*,*)'starting grid level 2'
46  continue

    call reyneq(akmax,1.0d0,ak1,bearx1,beary1,cohimx1,cohjmx1,      &
    		    f2p,h1,himax1,                                      &
    		    himin1,hjmax1,hjmin1,hm,hnew1,p1,pold1,p0,          &
    		    recssi1,recssj1,res1,su01,xl,xref1,yref1,1,         &
    		    idisc,mitp(2),nxm(1),nym(1),nx1,ny1,0,GS_Iterations,negp,2)

    if(mnest.eq.0) then
	    write(*,*)'Normalized residual = ',ak1
	    if(ak1.gt.akmax) then
	        goto 46
	    else
	        goto 50
	    endif
    endif

    ak=1.d0
    icount=0
    numPressureWarnings = 0
    stagnationCount = 0
    mnest = max_mg_nest
    
    do while(ak.gt.akmax .and. icount.lt.MAX_VCYCLE_ITR)
	    icount=icount+1
	    call vcycle(2,mnest)
	    write(*,*)'Normalized residual = ',ak
	    previousAK(icount) = ak
        call NestingCheck(previousAK, numPressureWarnings, stagnationCount, numStagnations, iCount, 0)
    enddo

50  continue

    !call prolong(p1,p,nxm(1),nym(1),nx,ny,nx1,ny1)
    call prolongCorrection(p, p1, xref, yref, xref1, yref1, &
                           nx, ny, nx1, ny1, enclosingCoarseRectX, enclosingCoarseRectY)

65  continue
    write(*,*)'starting finest grid, level 1'
    call DisplayAttitude
    call IncreaseNesting
    
54  continue

    
    !call InitKrylovAccel()
    call reyneq(akmax,1.d0,ak0,bearx,beary,cohimx,cohjmx,           &
    		    f2p,h,himax,himin,hjmax,hjmin,hm,hnew,              &
    		    p,pold,p0,recssi,recssj,res,su0,xl,xref,yref,1,     &
    		    idisc,mitp(1),nx,ny,nx,ny,0,GS_Iterations,negp,1)

    if(mnest.eq.0) then
	    write(*,*)'Normalized residual = ',ak0
	    if(ak0.gt.akmax) then
	        goto 54
	    else
	        goto 51
	    endif
    endif


    ak=1.d0
    icount=0
    !akPrevious = 1000000.0
    numPressureWarnings = 0
    stagnationCount = 0
    numStagnations = 0
    wantWCycle = .false.
    do while (ak.gt.akmax .and. icount.lt.MAX_VCYCLE_ITR)
	    icount = icount + 1
	    call vcycle(1, mnest)
	    write(*,*) 'Normalized residual = ', ak
        
        
        
        if (wantWCycle .eq. .true.) then
        
            !Four Levels        
            takenFirst = .false.
            if (ak.gt.akmax .and. mnest .gt. 0) then
                takenFirst = .true.
                mnest = mnest - 1
                call vcycle(1, mnest)
                write(*,*) 'Normalized residual = ', ak, 'Levels: ', mnest
            endif
            
            !Three Levels
            if (ak.gt.akmax .and. mnest .gt. 0) then
                mnest = mnest - 1
                call vcycle(1, mnest)
                write(*,*) 'Normalized residual = ', ak, 'Levels: ', mnest
                mnest = mnest + 1
            endif
            
            !Four Levels
            if (ak.gt.akmax .and. takenFirst .eq. .true.) then
                call vcycle(1, mnest)
                write(*,*) 'Normalized residual = ', ak, 'Levels: ', mnest
            endif

            !Three Levels
            if (ak.gt.akmax .and. mnest .gt. 0) then
                mnest = mnest - 1
                call vcycle(1, mnest)
                write(*,*) 'Normalized residual = ', ak, 'Levels: ', mnest
                mnest = mnest + 1
            endif

            !Four Levels
            if (ak.gt.akmax .and. takenFirst .eq. .true.) then
                call vcycle(1, mnest)
                write(*,*) 'Normalized residual = ', ak, 'Levels: ', mnest
                mnest = mnest + 1
            endif

            wAK = ak
            
	        call vcycle(1, mnest)
	        write(*,*) 'Normalized residual = ', ak
	    endif
	    
	    previousAK(icount) = ak
#ifdef _DEBUG
!	    call outputStagnationData
#endif
	    
!	    if (mod(icount, 4) .eq. 0) then
!	        !do reynolds directed
!            call ReynoldsDirected(89, 181, 5, 65)
!            call ReynoldsDirected(89, 181, 320, 397)
!	    endif
	    
	    if (icount .eq. 16) then
	        write(*,*) 'Turning off w-cycle'
	        wantWCycle = .false.
	        call ReduceNesting()
!	        pause
	        !call outputStagnationData
!	        call ReynoldsDirected
	    endif
            
        if (wantWCycle .eq. .true. .and. wAK .lt. ak .and. mnest .gt. 3) then
            call ReduceNesting()
        else
            call NestingCheck(previousAK, numPressureWarnings, stagnationCount, numStagnations, iCount, 1)
        endif

	    if (icount .eq. 10) then
	        call IncreaseNesting
	        wantWCycle = .true.
	    endif
    enddo

51  continue
	
	if (ak .le. akmax) then
	    write(*,*)'Convergence on the finest grid obtained.'
	    !call outputStagnationData
	else
	    write(*,*)'WARNING: Convergence on the finest grid NOT obtained.'
	endif

    if (doHumidity .eq. .true.) then
	    call SavePressure
	    call ModifyHumidPress
	    call getfxp(1)
	    call RestorePressure
    else
	    call getfxp(1)
    endif

    return
    end
    

!================================================
    subroutine NestingCheck(previousAK, numPressureWarnings, stagnationCount, numStagnations, iCount, wantStagnation)
!================================================
    implicit none
    
    integer, parameter :: MAX_VCYCLE_ITR = 80
    integer, parameter :: MAX_PRESSURE_WARNINGS = 40
    real(8) :: previousAK(MAX_VCYCLE_ITR+20), percentChange
    integer :: stagnationCount, numStagnations, numPressureWarnings
    integer :: iCount, wantStagnation
    
    
	! adjust Multigrid nesting if necessary (8/12/05):
	! if we've gone 13 iterations or 13 iterations since the last stagnation
	stagnationCount = stagnationCount + 1
    if (wantStagnation .eq. 1 .and. stagnationCount .ge. 13) then
	    ! try to figure out if the method is stagnating by % change in the residual over the last X iterations
	    percentChange = abs((previousAK(icount) - previousAK(icount-10)) / previousAK(icount-10))
	    if (percentChange .lt. 0.50) then
	        write(*,*) 'The current iteration seems to be stagnating'
	        !call outputStagnationData()
	        numStagnations = numStagnations + 1
            call ReduceNesting
	        stagnationCount = 0
	        numStagnations = 0
	        numPressureWarnings = 0
	    endif
	!if the residual increases or we get a number of pressure warnings then reduce the nesting
	else if (iCount .ge. 2) then
	    if (iCount .gt. 2 .and. (previousAK(iCount) .gt. previousAK(iCount-1)) .or. (numPressureWarnings .ge. MAX_PRESSURE_WARNINGS)) then
            
            !write(*,*) 'WARNING: not reducing nesting'
	        !return
	        
            call ReduceNesting
            numPressureWarnings = 0
            stagnationCount = 0
        endif
    endif
    
    return
    end
    
    
!========================================================    
    subroutine outputStagnationData
!========================================================
    use Q4_globals
    use PressureWarnings
    use NestingInfo
    implicit none    
    
    real(8) :: ak_std
    integer :: negp, i, j
    real(8) :: resArray(nx, ny), su0
    
	do i=1,nx
	    do j=1,ny
	        resArray(i, j) = 0.0
	    enddo
	enddo
    
!    call reyneq(akmax,1.d0,ak0,bearx,beary,cohimx,cohjmx,           &
!    		f2p,h,himax,himin,hjmax,hjmin,hm,hnew,              &
!    		p,pold,p0,recssi,recssj,res,su0,xl,xref,yref,1,     &
!    		idisc,mitp(1),nx,ny,nx,ny,0,GS_Iterations,negp,1)
    
    call GetReynoldsResidualAK(ak_std, bearx,beary,cohimx,cohjmx, &
                 f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,recssi,recssj,resArray, su0,    &
     	         xref,yref, idisc,nx,ny,nx,ny,0,GS_Iterations,negp)
     	         
     	         
    call outputArray(resArray, nx, ny, 'akArrayStag')
    call outputArray(res, nx, ny, 'resArrayStag')
    call outputArray(p, nx, ny, 'pArrayStag')
    
    return
    end