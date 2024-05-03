!====================================================================
      subroutine GetReynoldsResidualAK(ak, bearx,beary,cohimx,cohjmx,                     &
                 f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,recssi,recssj,res, su0,    &
     	         xref,yref, idisc,nx,ny,ndx,ndy,ml,nsp,negp)
!====================================================================
    use Q4_sizes
    use PressureWarnings
    implicit none
    
    real(8) :: akmax,akin,ak,bearx,beary,cohimx,cohjmx,     &
     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
     	   p0,recssi,recssj,res,su0,xl,xref,yref
    real(8) :: aw, ae, as, an, ap, su
    !integer :: nrey
         	   
    integer :: idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp
    integer :: i, j, im1, ip1, jm1, jp1, nxm1, nxm2, nym1, nym2
    integer :: nsweep
    
    real(8) :: deltax, deltay, delxi, delxim1, delyj, delyjm1, resn
    real(8) :: dpi, dmi, dpj, dmj, fpi, fmi, fpj, fmj, sliderPlaneHeight
    real(8) :: pim0, him, pim, qn1, qn2, qni, pjm0, hjm, pjm, qnj
    real(8) :: ymyp, xmxp, qnim, qnip, qnjm, qnjp
    real(8) :: ppi, pmi, ppj, pmj, api, ami, apj, amj
    real(8) :: apitrm, amitrm, apjtrm, amjtrm, temp
    real(8) :: reynoldsNumber, conv, diff, largestReynoldsNumber
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
     	su0(ndx,ndy)

    nxm1=nx-1
    nxm2=nx-2
    nym1=ny-1
    nym2=ny-2
    largestReynoldsNumber = 0.0
    
    ak=0.d0

    
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

	        if(ml.ne.0) su(i,j) = su(i,j) + su0(i,j)

            !calculate residuals here
            !we should probably calculate this at the end with updated coefs

            temp = ap(i,j) * p(i,j)
	        res(i,j)=(aw(i,j)*p(im1,j) + ae(i,j)*p(ip1,j)           &
        		    + as(i,j)*p(i,jm1) + an(i,j)*p(i,jp1) - temp+su(i,j))


	        res(i, j) = dabs(res(i, j)/temp)
            
            ak = ak + res(i,j)
            !output for the residual restriction
	             
	    enddo
	enddo

    ak=ak

    return
    end    
    
!
!!====================================================================
!      subroutine GetReynoldsStencil(akmax,akin,ak,bearx,beary,cohimx,cohjmx,    &
!     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,p0,recssi,recssj,res,   &
!     	   su0,xl,xref,yref, idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp,       &
!     	   aw, ae, as, an, ap, su, nrey)
!!====================================================================
!    use Q4_sizes
!    use PressureWarnings
!    implicit none
!    
!    real(8) :: akmax,akin,ak,bearx,beary,cohimx,cohjmx,     &
!     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
!     	   p0,recssi,recssj,res,su0,xl,xref,yref
!    real(8) :: aw, ae, as, an, ap, su
!    integer :: nrey
!         	   
!    integer :: idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp
!    integer :: i, j, im1, ip1, jm1, jp1, nxm1, nxm2, nym1, nym2
!    integer :: nsweep
!    
!    real(8) :: deltax, deltay, delxi, delxim1, delyj, delyjm1, resn
!    real(8) :: dpi, dmi, dpj, dmj, fpi, fmi, fpj, fmj, sliderPlaneHeight
!    real(8) :: pim0, him, pim, qn1, qn2, qni, pjm0, hjm, pjm, qnj
!    real(8) :: ymyp, xmxp, qnim, qnip, qnjm, qnjp
!    real(8) :: ppi, pmi, ppj, pmj, api, ami, apj, amj
!    real(8) :: apitrm, amitrm, apjtrm, amjtrm, temp
!    
!    !debugging info
!    real(8) :: cur_residual
!    
!    
!    dimension p(nx, ny)
!    dimension bearx(ndx,ndy),beary(ndx,ndy),  &
!     	pold(ndx,ndy), qni(nx,ny),qnj(nx,ny), h(ndx,ndy), hnew(ndx,ndy),    &
!     	recssi(ndx,ndy), recssj(ndx,ndy), xref(ndx), yref(ndy),             &
!     	cohimx(ndx,ndy),cohjmx(ndx,ndy), himax(ndx,ndy),himin(ndx,ndy),     &
!     	hjmax(ndx,ndy),hjmin(ndx,ndy)
!    dimension ae(nx,ny),aw(nx,ny),as(nx,ny),an(nx,ny),    &
!     	ap(nx,ny),su(nx,ny),res(ndx,ndy), su0(ndx,ndy)
!
!    nxm1=nx-1
!    nxm2=nx-2
!    nym1=ny-1
!    nym2=ny-2
!    
!	do i=2,nx
!	    im1=i-1
!	    ip1=i+1
!	    do j=2,ny-1
!	        jm1=j-1
!	        jp1=j+1
!	        pim0 = (p(i,j)+p(im1,j))/2.d0         !avg pressure along volume wall
!	        sliderPlaneHeight = (hnew(i,j)+hnew(im1,j))/2.d0   !avg slider plane spacing along volume
!
!            !slip at maximum recess along volume boundary
!	        him=himax(i,j) + sliderPlaneHeight
!	        pim = pim0 * him
!	        call flow(pim,1.d0,1.d0,1.d0,qn1)
!	        qn1 = qn1 * him**3
!            
!            !slip at minimum recess along volume boundary
!	        him = dmax1(2d-10, himin(i,j) + sliderPlaneHeight)
!	        pim = pim0 * him
!	        call flow(pim,1.d0,1.d0,1.d0,qn2)
!	        qn2 = qn2 * him**3
!            
!            !weighted average of slip at maximum and minimum recess heights
!            !weighting is based on ratio of average to maximum and minimum 
!            !height along volume boundary
!	        qni(i,j)=qn1*cohimx(i,j)+qn2*(1.d0-cohimx(i,j))
!	    enddo
!	enddo
!
!    do i=2,nx-1
!	    im1=i-1
!	    ip1=i+1
!	    do j=2,ny
!	        jm1=j-1
!	        jp1=j+1
!	        pjm0=(p(i,j)+p(i,jm1))/2.d0
!	        sliderPlaneHeight=(hnew(i,j)+hnew(i,jm1))/2.d0
!
!	        hjm=hjmax(i,j)+sliderPlaneHeight
!	        pjm=pjm0*hjm
!	        call flow(pjm,1.d0,1.d0,1.d0,qn1)
!	        qn1=qn1*hjm**3
!
!	        hjm= dmax1(2d-10, hjmin(i,j)+sliderPlaneHeight)
!	        pjm=pjm0*hjm
!	        call flow(pjm,1.d0,1.d0,1.d0,qn2)
!	        qn2=qn2*hjm**3
!
!	        qnj(i,j)=qn1*cohjmx(i,j)+qn2*(1.d0-cohjmx(i,j))
!	    enddo
!    enddo
!
!    do j=2,nym1
!	    jm1=j-1
!	    jp1=j+1
!	    delyj=yref(jp1)-yref(j)
!	    delyjm1=yref(j)-yref(jm1)
!	    ymyp=(yref(jp1)-yref(jm1))/2.d0
!	    deltay=ymyp
!
!	    do i=2,nxm1
!	        im1=i-1
!	        ip1=i+1
!	        delxi=xref(ip1)-xref(i)
!	        delxim1=xref(i)-xref(im1)
!	        xmxp=(xref(ip1)-xref(im1))/2.d0
!	        deltax =xmxp
!
!	        qnim=qni(i,j)*f2p
!	        qnip=qni(ip1,j)*f2p
!	        qnjm=qnj(i,j)*f2p
!	        qnjp=qnj(i,jp1)*f2p
!
!	        qnim=qnim*(p(im1,j)+p(i,j))/2.d0
!	        qnip=qnip*(p(ip1,j)+p(i,j))/2.d0
!	        qnjm=qnjm*(p(i,jm1)+p(i,j))/2.d0
!	        qnjp=qnjp*(p(i,jp1)+p(i,j))/2.d0
!
!	        dpi=qnip/delxi*deltay
!	        dmi=qnim/delxim1*deltay
!	        dpj=qnjp/delyj*deltax
!	        dmj=qnjm/delyjm1*deltax
!
!            !recessX vars hold average recess at volume boundary
!	        fpi=(bearx(ip1,j)+bearx(i,j))/4.d0*(recssi(ip1,j)+  &
!                (hnew(i,j)+hnew(ip1,j))/2.d0)*deltay
!
!	        fmi=(bearx(im1,j)+bearx(i,j))/4.d0*(recssi(i,j) +   &
!     	        (hnew(i,j)+hnew(im1,j))/2.d0)*deltay
!
!	        fpj=(beary(i,jp1)+beary(i,j))/4.d0*(recssj(i,jp1) + &
!	            (hnew(i,j) + hnew(i,jp1))/2.d0)*deltax
!
!	        fmj=(beary(i,jm1)+beary(i,j))/4.d0*(recssj(i,j) +   &
!                (hnew(i,j)+hnew(i,jm1))/2.d0)*deltax
!
!	        ppi=fpi/dpi
!	        pmi=fmi/dmi
!	        ppj=fpj/dpj
!	        pmj=fmj/dmj
!
!	        if(idisc.eq.0)then
!                !upwind
!		        api=1.d0
!		        ami=1.d0
!		        apj=1.d0
!		        amj=1.d0
!	        else if(idisc.eq.1) then
!                !hybrid
!		        api=dmax1(0.d0,1.d0-0.5d0*dabs(ppi))
!		        ami=dmax1(0.d0,1.d0-0.5d0*dabs(pmi))
!		        apj=dmax1(0.d0,1.d0-0.5d0*dabs(ppj))
!		        amj=dmax1(0.d0,1.d0-0.5d0*dabs(pmj))
!	        else
!                !power law, idisc=2 (default)
!		        apitrm=(1.d0-0.1d0*dabs(ppi))**5
!		        amitrm=(1.d0-0.1d0*dabs(pmi))**5
!		        apjtrm=(1.d0-0.1d0*dabs(ppj))**5
!		        amjtrm=(1.d0-0.1d0*dabs(pmj))**5
!
!		        api=dmax1(0.d0,apitrm)
!		        ami=dmax1(0.d0,amitrm)
!		        apj=dmax1(0.d0,apjtrm)
!		        amj=dmax1(0.d0,amjtrm)
!	        endif
!
!	        aw(i,j)=dmi*ami+dmax1(fmi,0.d0)
!	        ae(i,j)=dpi*api+dmax1(-fpi,0.d0)
!	        as(i,j)=dmj*amj+dmax1(fmj,0.d0)
!	        an(i,j)=dpj*apj+dmax1(-fpj,0.d0)
!	        ap(i,j)=ae(i,j)+aw(i,j)+as(i,j)+an(i,j) + dmax1(fpi-fmi+fpj-fmj,0.d0)
!	        su(i,j)=dmax1(fmi-fpi+fmj-fpj,0.d0)*p(i,j)
!
!
!            !multigrid restriction source term
!            !coarse grid MG method is to solve A(u) = f + A(v) + r
!            !where u is exact pressure and v is restricted pressure from finer grid
!            !if we're moving from a finer to coarser coarser grid and it's the first of our fixed point iterations
!            !then add in our residual + solution approx from the finer grid pressure approximation
!            !remember to multiply residual times area because we divided it by area on the finer level
!	        if(ml.lt.0.and.nrey.eq.1) then
!	            su0(i,j) = res(i,j)*ymyp*xmxp                           &  
!     	                    - aw(i,j)*p(im1,j) - ae(i,j)*p(ip1,j)       &
!     	                    - as(i,j)*p(i,jm1) - an(i,j)*p(i,jp1) + ap(i,j)*p(i,j)-su(i,j)
!            endif
!            
!	        if(ml.ne.0) su(i,j) = su(i,j) + su0(i,j)
!	        !if(ml.lt.0) su(i,j) = su(i,j) + su0(i,j)
!
!            !calculate residuals here
!            !we should probably calculate this at the end with updated coefs
!
!            temp = ap(i,j) * p(i,j)
!	        res(i,j)=(aw(i,j)*p(im1,j) + ae(i,j)*p(ip1,j)           &
!        		    + as(i,j)*p(i,jm1) + an(i,j)*p(i,jp1) - temp+su(i,j))
!
!
!	        ak = ak + dabs(res(i,j)/temp)
!            !output for the residual restriction
!	        res(i,j)=res(i,j)/(ymyp*xmxp)
!	    enddo
!	enddo
!
!    return
!    end
    
!    
!!====================================================================
!      subroutine GetReynoldsResidualSub(bearx,beary,cohimx,cohjmx,                     &
!                 f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,recssi,recssj,res, su0,    &
!     	         xref,yref, idisc,nx,ny,ndx,ndy,ml,nsp,negp)
!!====================================================================
!    use Q4_sizes
!    use PressureWarnings
!    implicit none
!    
!    real(8) :: akmax,akin,ak,bearx,beary,cohimx,cohjmx,     &
!     	   f2p,h,himax,himin,hjmax,hjmin,hm,hnew,p,pold,    &
!     	   p0,recssi,recssj,res,su0,xl,xref,yref
!    real(8) :: aw, ae, as, an, ap, su
!    !integer :: nrey
!         	   
!    integer :: idirec,idisc,itrey,nx,ny,ndx,ndy,ml,nsp,negp
!    integer :: i, j, im1, ip1, jm1, jp1, nxm1, nxm2, nym1, nym2
!    integer :: nsweep
!    
!    real(8) :: deltax, deltay, delxi, delxim1, delyj, delyjm1, resn
!    real(8) :: dpi, dmi, dpj, dmj, fpi, fmi, fpj, fmj, sliderPlaneHeight
!    real(8) :: pim0, him, pim, qn1, qn2, qni, pjm0, hjm, pjm, qnj
!    real(8) :: ymyp, xmxp, qnim, qnip, qnjm, qnjp
!    real(8) :: ppi, pmi, ppj, pmj, api, ami, apj, amj
!    real(8) :: apitrm, amitrm, apjtrm, amjtrm, temp
!    real(8) :: reynoldsNumber, conv, diff, largestReynoldsNumber
!    !debugging info
!    real(8) :: cur_residual
!    
!    dimension p(nx, ny)
!    dimension bearx(ndx,ndy),beary(ndx,ndy),  &
!     	pold(ndx,ndy), qni(nx,ny),qnj(nx,ny), h(ndx,ndy), hnew(ndx,ndy),    &
!     	recssi(ndx,ndy), recssj(ndx,ndy), xref(ndx), yref(ndy),             &
!     	cohimx(ndx,ndy),cohjmx(ndx,ndy), himax(ndx,ndy),himin(ndx,ndy),     &
!     	hjmax(ndx,ndy),hjmin(ndx,ndy)
!    dimension ae(nx,ny),aw(nx,ny),as(nx,ny),an(nx,ny),    &
!     	ap(nx,ny),su(nx,ny),res(ndx,ndy), su0(ndx,ndy)
!
!    nxm1=nx-1
!    nxm2=nx-2
!    nym1=ny-1
!    nym2=ny-2
!    largestReynoldsNumber = 0.0
!	do i=2,nx
!	    im1=i-1
!	    ip1=i+1
!	    do j=2,ny-1
!	        jm1=j-1
!	        jp1=j+1
!	        pim0 = (p(i,j)+p(im1,j))/2.d0         !avg pressure along volume wall
!	        sliderPlaneHeight = (hnew(i,j)+hnew(im1,j))/2.d0   !avg slider plane spacing along volume
!
!            !slip at maximum recess along volume boundary
!	        him=himax(i,j) + sliderPlaneHeight
!	        pim = pim0 * him
!	        call flow(pim,1.d0,1.d0,1.d0,qn1)
!	        qn1 = qn1 * him**3
!            
!            !slip at minimum recess along volume boundary
!	        him = dmax1(2d-10, himin(i,j) + sliderPlaneHeight)
!	        pim = pim0 * him
!	        call flow(pim,1.d0,1.d0,1.d0,qn2)
!	        qn2 = qn2 * him**3
!            
!            !weighted average of slip at maximum and minimum recess heights
!            !weighting is based on ratio of average to maximum and minimum 
!            !height along volume boundary
!	        qni(i,j)=qn1*cohimx(i,j)+qn2*(1.d0-cohimx(i,j))
!	    enddo
!	enddo
!
!    do i=2,nx-1
!	    im1=i-1
!	    ip1=i+1
!	    do j=2,ny
!	        jm1=j-1
!	        jp1=j+1
!	        pjm0=(p(i,j)+p(i,jm1))/2.d0
!	        sliderPlaneHeight=(hnew(i,j)+hnew(i,jm1))/2.d0
!
!	        hjm=hjmax(i,j)+sliderPlaneHeight
!	        pjm=pjm0*hjm
!	        call flow(pjm,1.d0,1.d0,1.d0,qn1)
!	        qn1=qn1*hjm**3
!
!	        hjm= dmax1(2d-10, hjmin(i,j)+sliderPlaneHeight)
!	        pjm=pjm0*hjm
!	        call flow(pjm,1.d0,1.d0,1.d0,qn2)
!	        qn2=qn2*hjm**3
!
!	        qnj(i,j)=qn1*cohjmx(i,j)+qn2*(1.d0-cohjmx(i,j))
!	    enddo
!    enddo
!
!    do j=2,nym1
!	    jm1=j-1
!	    jp1=j+1
!	    delyj=yref(jp1)-yref(j)
!	    delyjm1=yref(j)-yref(jm1)
!	    ymyp=(yref(jp1)-yref(jm1))/2.d0
!	    deltay=ymyp
!
!	    do i=2,nxm1
!	        im1=i-1
!	        ip1=i+1
!	        delxi=xref(ip1)-xref(i)
!	        delxim1=xref(i)-xref(im1)
!	        xmxp=(xref(ip1)-xref(im1))/2.d0
!	        deltax =xmxp
!
!	        qnim=qni(i,j)*f2p
!	        qnip=qni(ip1,j)*f2p
!	        qnjm=qnj(i,j)*f2p
!	        qnjp=qnj(i,jp1)*f2p
!
!	        qnim=qnim*(p(im1,j)+p(i,j))/2.d0
!	        qnip=qnip*(p(ip1,j)+p(i,j))/2.d0
!	        qnjm=qnjm*(p(i,jm1)+p(i,j))/2.d0
!	        qnjp=qnjp*(p(i,jp1)+p(i,j))/2.d0
!
!	        dpi=qnip/delxi*deltay
!	        dmi=qnim/delxim1*deltay
!	        dpj=qnjp/delyj*deltax
!	        dmj=qnjm/delyjm1*deltax
!
!            !recessX vars hold average recess at volume boundary
!	        fpi=(bearx(ip1,j)+bearx(i,j))/4.d0*(recssi(ip1,j)+  &
!                (hnew(i,j)+hnew(ip1,j))/2.d0)*deltay
!
!	        fmi=(bearx(im1,j)+bearx(i,j))/4.d0*(recssi(i,j) +   &
!     	        (hnew(i,j)+hnew(im1,j))/2.d0)*deltay
!
!	        fpj=(beary(i,jp1)+beary(i,j))/4.d0*(recssj(i,jp1) + &
!	            (hnew(i,j) + hnew(i,jp1))/2.d0)*deltax
!
!	        fmj=(beary(i,jm1)+beary(i,j))/4.d0*(recssj(i,j) +   &
!                (hnew(i,j)+hnew(i,jm1))/2.d0)*deltax
!            
!
!	        ppi=fpi/dpi
!	        pmi=fmi/dmi
!	        ppj=fpj/dpj
!	        pmj=fmj/dmj
!
!	        if(idisc.eq.0)then
!                !upwind
!		        api=1.d0
!		        ami=1.d0
!		        apj=1.d0
!		        amj=1.d0
!	        else if(idisc.eq.1) then
!                !hybrid
!		        api=dmax1(0.d0,1.d0-0.5d0*dabs(ppi))
!		        ami=dmax1(0.d0,1.d0-0.5d0*dabs(pmi))
!		        apj=dmax1(0.d0,1.d0-0.5d0*dabs(ppj))
!		        amj=dmax1(0.d0,1.d0-0.5d0*dabs(pmj))
!	        else
!                !power law, idisc=2 (default)
!		        apitrm=(1.d0-0.1d0*dabs(ppi))**5
!		        amitrm=(1.d0-0.1d0*dabs(pmi))**5
!		        apjtrm=(1.d0-0.1d0*dabs(ppj))**5
!		        amjtrm=(1.d0-0.1d0*dabs(pmj))**5
!
!		        api=dmax1(0.d0,apitrm)
!		        ami=dmax1(0.d0,amitrm)
!		        apj=dmax1(0.d0,apjtrm)
!		        amj=dmax1(0.d0,amjtrm)
!	        endif
!
!	        aw(i,j)=dmi*ami+dmax1(fmi,0.d0)
!	        ae(i,j)=dpi*api+dmax1(-fpi,0.d0)
!	        as(i,j)=dmj*amj+dmax1(fmj,0.d0)
!	        an(i,j)=dpj*apj+dmax1(-fpj,0.d0)
!	        ap(i,j)=ae(i,j)+aw(i,j)+as(i,j)+an(i,j) + dmax1(fpi-fmi+fpj-fmj,0.d0)
!	        su(i,j)=dmax1(fmi-fpi+fmj-fpj,0.d0)*p(i,j)
!            
!	        if(ml.ne.0) su(i,j) = su(i,j) + su0(i,j)
!
!            !calculate residuals here
!            !we should probably calculate this at the end with updated coefs
!
!            temp = ap(i,j) * p(i,j)
!	        res(i,j)=(aw(i,j)*p(im1,j) + ae(i,j)*p(ip1,j)           &
!        		    + as(i,j)*p(i,jm1) + an(i,j)*p(i,jp1) - temp+su(i,j))
!
!
!	        ak = ak + dabs(res(i,j)/temp)            
!            
!            !output for the residual restriction
!	        res(i,j)=res(i,j)/(ymyp*xmxp)
!	                
!	    enddo
!	enddo
!
!    return
!    end    
!    
