!c==============================================
      subroutine flow (pn,hbar,h2bar,h3bar,qn)
!c==============================================
!c     this subroutine calculates the poiseuille flow
!c     under the slider based on the boltzman equation
!c     the first and second derivatives are also com-
!c     computed. an asymptotic form is used for large values.

      !include 'common.fi'
      use Q4_globals
      use FKLookupTable
      implicit none
      
      real(8) :: A, AIJ, ALZ, B, C11, C12, C22, DEL, EV, F1, F2, FO, H2BAR, H3BAR, HBAR
      real(8) :: PN, Q2, QN, R1, S1, U1, V, V1, V2, V3, V4, V5, XI, Z, qn1, qn2, z0, z1, z2
      integer :: i, loc
      dimension a(20),b(20)


    !we'll use a large table lookup to compute FK model.
    if (iqpo .eq. 5) then
        z = d0 * pn;
        if (z .ge. zBaseLow .and. z .lt. zBaseMid) then
            loc = floor((z - zBaseLow) / zStepLow) + 1;
            z1=zdatLow(loc);
            qn1=qndatLow(loc);
            z2=zdatLow(loc+1);
            qn2=qndatLow(loc+1);
            qn=(qn2-qn1)*(z-z1)/(z2-z1)+qn1;
            return
        else if (z .ge. zBaseMid .and. z .lt. zBaseHigh) then
            loc = floor((z - zBaseMid) / zStepMid) + 1;
            z1=zdatMid(loc);
            qn1=qndatMid(loc);
            z2=zdatMid(loc+1);
            qn2=qndatMid(loc+1);
            qn=(qn2-qn1)*(z-z1)/(z2-z1)+qn1;
            return
        else if (z .ge. zBaseHigh .and. z .le. zBaseCeil) then
            loc = floor((z - zBaseHigh) / zStepHigh) + 1;
            z1=zdatHigh(loc);
            qn1=qndatHigh(loc);
            z2=zdatHigh(loc+1);
            qn2=qndatHigh(loc+1);
            qn=(qn2-qn1)*(z-z1)/(z2-z1)+qn1;    
            return
        else
            !if we're out of the range of the database then jump to full FK computation
            goto 654
        endif  

!cha   use the continuum model.....
      else if (iqpo.eq.0) then
	qn = h3bar
	return

!cha   use the first order slip model....
      else if (iqpo.eq.1) then
	qn = h3bar + h2bar*t1/pn
	return

!cha   use the second order slip model.....

      else if (iqpo.eq.2) then
	r1 = t2/pn/pn
	qn = h3bar + h2bar*t1/pn + hbar*r1
	return

! Lion  use the new first order slip model....
      else if (iqpo.eq.11) then
	qn = h3bar + (4.0/6.0)*h2bar*t1/pn
	return

! Lion  use the new 1.5 order model
	else if (iqpo.eq.15) then
	r1 = t2/pn/pn
	qn = h3bar + h2bar*t1/pn + (8.0/3.0/6.0)*hbar*r1
	return

! Lion  use the new second order model
	else if (iqpo.eq.22) then
	r1 = t2/pn/pn
	qn = h3bar + (4.0/6.0)*h2bar*t1/pn + 0.5*hbar*r1
	return

! Lion  use the user-defined slip model.....
      else if (iqpo.eq.8) then
	r1 = t2/pn/pn
	qn = h3bar + slip_beta/6.0*h2bar*t1/pn + slip_gamma/6.0*hbar*r1
	return

! Lion  use the pressure gradient model
	else if (iqpo.eq.6) then
	r1 = t2/pn/pn
	qn = h3bar + h2bar*t1/pn + 2.0*hbar*r1
	return

!cha   use the asymtotic fk model.....
      else if (iqpo.eq.3) then
	r1 = t2/pn/pn
	s1 = t3/(pn*pn*pn)
	u1 = t4/(pn*pn*pn*pn)
	qn = 1.d0+t1/pn + r1 - s1 - u1
	return
!c
!cha   use the full fk model.....
    else if (iqpo.eq.4) then
654	    z  = d0*pn
	    if (icoe.eq.1) go to 15
	    a(1) = 0.d0
	    a(2) = -1.d0
	    b(1) = -pir
	    b(2) = 1.5d0*(1.d0-gama)
	    do i = 3,nter
	        xi  = dfloat(i)
	        aij = xi*(xi-1.d0)*(xi-2.d0)
	        a(i) = (-2.d0*a(i-2))/aij
	        b(i) = (-2.d0*b(i-2)-(3.d0*xi*xi-6.d0*xi+2.d0)*a(i))/aij
	    enddo
	    icoe = 1
    15	if (z.ge.1.1d0)	 then
	        v = 1.889881575d0*z**0.66666666666666666d0
	        v1 = 12.d0*v
	        v2 = 24.d0*v*v1
	        v3 = 36.d0*v*v2
	        v4 = 48.d0*v*v3
	        v5 = 60.d0*v*v4
	        ev = pit*exp(-v)
	        fo = ev*(1.d0-1.d0/v1+25.d0/v2-1225.d0/v3+89425.d0/v4-7263025.d0/v5)
	        f1 = ev*dsqrt(v/3.d0)*(1.d0+5.d0/v1-35.d0/v2+665.d0/v3+9625.d0/v4-9284275.d0/v5)
	        f2 = ev*(v/3.d0)*(1.d0+17.d0/v1-35.d0/v2+1925.d0/v3-175175.d0/v4+22247225.d0/v5)
	    else
	        fo   = 0.d0
	        f1   = 0.d0
	        f2   = 0.d0
	        alz  = log(z)
	        do i = nter,1,-1
	            xi	= dfloat(i)
	            aij = 1.d0/(1.d0+xi)
	            if (i.eq.1) then
	            fo = ((alz*xi+1.d0)*a(i) + b(i)*xi + fo)
	            go to 12
	            end if
	            fo = ((alz*xi+1.d0)*a(i) + b(i)*xi + fo)*z
12	            f1 = (a(i)*alz + b(i) + f1)*z
	            f2 = (aij*((alz-aij)*a(i)+b(i)) + f2)*z
	        enddo
	        fo = -0.5d0*fo
	        f1 = 0.5d0*(f1 + 1.d0)
	        f2 = pir/4.d0 - 0.5d0*(f2 + 1.d0)*z
	    endif
	    c11 = 8.d0-z*z*z*(pir/12.d0-z/16.d0)-2.d0*z*(4.d0+z*z)*fo-(16.d0+z*z*(8.d0+z*z/8.d0))*f1-z*(16.d0+z*z)*f2
	    c22 = 1.d0- 2.d0*f1
	    c12 = 2.d0 - (pir/2.d0)*z + z*z/4.d0-2.d0*z*fo-(4.d0 +z*z/2.d0)*f1 - 2.d0*z*f2
	    del = c11*c22 - c12*c12
	    q2 = (pir/del)*(c11 - z*z*(c12/6.d0-z*z*c22/144.d0))
	    qn = (-1.d0/z+q2)*(6.d0/z)
	    return
	    
    endif

!	if inverse knudsen no. is out of range of database
!	use first order slip
55	continue
	r1   = t2/pn/pn
	s1   = t3/(pn*pn*pn)
	u1   = t4/(pn*pn*pn*pn)
	qn   = 1.d0+t1/pn + r1 - s1 - u1
	return
	
56	continue
	qn=-6.d0/pir * log(z0)/z0



    return
    end


!   old subroutine to compute small FK tables.  Has been replaced by large database creation
!===============================================
    subroutine Create_FK_Tables()
!===============================================
    use Q4_globals
    implicit none
      
!     We might want to use the interpolation of FK tables for
!     Accomodation coefs != 1 since the curve fits aren't all that
!     good for lower accomodataion coefs      
      
    integer :: i
    real(8) :: fkFlowLower(19),fkFlowLowerMid(19),fkFlowUpperMid(19),fkFlowUpper(19)
!      data dInvKnud/0.0100,0.0150,0.0200,0.0250,0.0300,0.0350,
!     &  0.0400,0.0500,0.0600,0.0700,0.0800,0.0900,0.1000,0.1500,
!     &  0.2000,0.2500,0.3000,0.3500,0.4000,0.5000,0.6000,0.7000,
!     &  0.8000,0.9000,1.0000,1.5000,2.0000,2.5000,3.0000,3.5000,
!     &  4.0000,5.0000,6.0000,7.0000,8.0000,9.0000,10.0000,15.0000,
!     &  20.0000,25.0000,30.0000,35.0000,40.0000,50.0000,60.0000,
!     &  70.0000,80.0000,90.0000,100.0000/
!      data dPou/3.5640,3.3120,3.1450,3.0210,2.9230,2.8410,
!     &  2.7720,2.6590,2.5700,2.4970,2.4360,2.3840,2.3380,2.1750,2.0720,
!     &  2.0000,1.9470,1.9060,1.8730,1.8260,1.7950,1.7740,1.7600,1.7510,
!     &  1.7460,1.7540,1.7890,1.8340,1.8860,1.9430,2.0070,2.1410,2.2840,
!     &  2.4310,2.5870,2.7370,2.8980,3.6950,4.4940,5.3040,6.1020,6.8880,
!     &  7.6590,9.2290,10.7370,12.2630,13.7250,15.2450,16.7790/

      data fkFlowUpper/2.8980,3.6950,4.4940,5.3040,6.1020,6.8880, &
                       7.6590,8.4440,9.2290,9.9830,10.7370,11.5000,12.2630,12.9940, &
                       13.7250,14.4850,15.2450,16.0120,16.7790/
      data fkFlowUpperMid/1.7460,1.7540,1.7890,1.8340,1.8860,1.9430, &
                          2.0070,2.0740,2.1410,2.2125,2.2840,2.3575,2.4310,2.5090,2.5870, &
                          2.6620,2.7370,2.8175,2.8980/
      data fkFlowLowerMid/2.3380,2.1750,2.0720,2.0000,1.9470,1.9060, &
                          1.8730,1.8495,1.8260,1.8105,1.7950,1.7845,1.7740,1.7670, &
                          1.7600,1.7555,1.7510,1.7485,1.7460/
      data fkFlowLower/3.5640,3.3120,3.1450,3.0210,2.9230,2.8410, &
                       2.7720,2.7155,2.6590,2.6145,2.5700,2.5335,2.4970,2.4665, &
                       2.4360,2.4100,2.3840,2.3610,2.3380/

!      data fkFlowUpper/2.7680,3.5780,4.3980,5.2220,6.0490,
!     &  6.8780,7.7080,8.5390,9.3700,10.2015,11.0330,11.8655,12.6980,
!     &  13.5305,14.3630,15.1955,16.0280,16.8605,17.6930/
!      data fkFlowUpperMid/1.5390,1.5540,1.5950,1.6490,1.7110,
!     &  1.7770,1.8460,1.9185,1.9910,2.0625,2.1340,2.2130,2.2920,2.3705,
!     &  2.4490,2.5285,2.6080,2.6880,2.7680/
!      data fkFlowLowerMid/2.0330,1.8950,1.8080,1.7480,1.7030,
!     &  1.6680,1.6410,1.6215,1.6020,1.5890,1.5760,1.5675,1.5590,
!     &  1.5535,1.5480,1.5450,1.5420,1.5405,1.5390/
!      data fkFlowLower/3.0600,2.8460,2.7070,2.6040,2.5220,
!     &  2.4540,2.3970,2.3495,2.3020,2.2650,2.2280,2.1975,2.1670,2.1410,
!     &  2.1150,2.0930,2.0710,2.0520,2.0330/

    do i = 1,19
        FKUpper(i) = fkFlowUpper(i)
        FKMidUpper(i) = fkFlowUpperMid(i)
        FKMidLower(i) = fkFlowLowerMid(i)
        FKLower(i) = fkFlowLower(i)
    enddo
    
    return 
    end
      
! interpolation of larger FK tables.  No longer used
!===============================================      
    subroutine FKInterp(pn, qn)
!===============================================
    use Q4_globals
    implicit none
    
    real(8) :: pn, hbar, h2bar, h3bar, qn
    real(8) :: Ainterp, Binterp, lowerVal, z, step
    integer iLowerIndex

    
    z = d0 * pn
    if ((z .ge. 10) .and. (z .lt. 100.0)) then
        step = 5.0
        ilowerIndex = idint(z / step) - 1  !get lower index in FK
        lowerVal = 10 + dble(iLowerIndex-1)*step !get the value at that index
        !now linear interp between the two values
        Ainterp = ((lowerVal+step) - z) / step
        Binterp = 1.0 - Ainterp
        qn = Ainterp * FKUpper(iLowerIndex) + Binterp * FKUpper(iLowerIndex+1)
        qn = qn/z * 6.0
    else if ((z .ge. 1.0) .and. (z .lt. 10.0)) then
        step = 0.5
        ilowerIndex = idint(z / step) - 1
        lowerVal = 1.0 + dble(iLowerIndex-1)*step
        Ainterp = ((lowerVal + step) - z) / step
        Binterp = 1.0 - Ainterp
        qn = Ainterp * FKMidUpper(iLowerIndex) + Binterp * FKMidUpper(iLowerIndex+1)
        qn = qn/z * 6.0
    else if ((z .ge. 0.1) .and. (z .lt. 1.0)) then
        step = 0.05
        ilowerIndex = idint(z / step) - 1
        lowerVal = 0.1 + dble(iLowerIndex-1)*step
        Ainterp = ((lowerVal + step) - z) / step
        Binterp = 1.0 - Ainterp
        qn = Ainterp * FKMidLower(iLowerIndex) + Binterp * FKMidLower(iLowerIndex+1)
        qn = qn/z * 6.0
    else if ((z .ge. 0.01) .and. (z .lt. 0.1)) then
        step = 0.005
        ilowerIndex = idint(z / step) - 1
        lowerVal = 0.01 + dble(iLowerIndex-1)*step
        Ainterp = ((lowerVal + step) - z) / step
        Binterp = 1.0 - Ainterp
        qn = Ainterp * FKLower(iLowerIndex) + Binterp * FKLower(iLowerIndex+1)
        qn = qn/z * 6.0      
    else  !Using limit
	    qn=-6.d0/pir*dlog10(z)/z
    endif
	
    return
    end


!   To create the FK datbase we'll use 3 seperate tables (lower, middle and high)
!   this will allow us to achieve fairly good accuracy given the amount of space we'll use
!   The FK tables change somewhat logrythmically so our lower tables will have much smaller 
!   spacing than our upper tables
!==================================================
    subroutine create_dbase
!==================================================

    use Q4_globals
    use FKLookupTable
    implicit none

    integer :: i, zSize
    real(8) :: dz, pn, qn, z

    write(*,*) 'Creating FK database...'
	if(iqpo .ne. 5) return
	!temporarily change iqpo to 4 to create data base
    iqpo=4
    write(*,*)qn
    !ranges for each seperate table
    zBaseLow = 0.01d0;
    zBaseMid = 1.0d0;
    zBaseHigh = 100.0d0;
    zBaseCeil = 10000.0d0;

    !Lower
    dz = zBaseLow * 0.01d0
    zStepLow = dz
    zSize = int((zBaseMid - zBaseLow) / dz) + 3
    allocate(zdatLow (zSize))
    allocate(qndatLow (zSize))
    z = zBaseLow
    i = 1
    do while (z .le. zBaseMid)
        !simply step along by dz and compute the flow at that inverse knudson number
        z = zBaseLow + (dz * (i-1))
        pn=z/d0
        call flow(pn,1.d0,1.d0,1.d0,qn)
        zdatLow(i)=z  !used for interpolation
        qndatLow(i)=qn  !actual flow value
        i = i+1
    enddo

    !Middle
    dz = zBaseMid * 0.01d0
    zStepMid = dz
    zSize = int((zBaseHigh - zBaseMid) / dz) + 3
    allocate(zdatMid (zSize))
    allocate(qndatMid (zSize))
    z = zBaseMid
    i = 1
    do while (z .le. zBaseHigh)
        z = zBaseMid + (dz * (i-1))
        pn=z/d0
        call flow(pn,1.d0,1.d0,1.d0,qn)
        zdatMid(i)=z
        qndatMid(i)=qn
        i = i+1
    enddo
    
    !High    
    dz = zBaseHigh * 0.01d0
    zStepHigh = dz
    zSize = int((zBaseCeil - zBaseHigh) / dz) + 3
    allocate(zdatHigh (zSize))
    allocate(qndatHigh (zSize))
    z = zBaseHigh
    i = 1
    do while (z .le. zBaseCeil)
        z = zBaseHigh + (dz * (i-1))
        pn=z/d0
        call flow(pn,1.d0,1.d0,1.d0,qn)
        zdatHigh(i)=z
        qndatHigh(i)=qn
        i = i+1
    enddo 

!   change iqpo back to 5
    iqpo=5

    write(*,*) 'Finished creating FK database'

    return
    end
    
    
!===============================================================
    subroutine massflow(bearx,beary,cohimx,cohjmx,f2p,  &
     			h,himax,himin,hjmax,hjmin,hm,           &
     			hnew,p, p0,res,recssi,recssj,xl,        &
     			xref,yref,incom,nx,ny,ndx,ndy)      
!===============================================================
    
    use Q4_Sizes
    implicit none
    
    real(8) :: bearx,beary,cohimx,cohjmx,f2p,           &
     			h,himax,himin,hjmax,hjmin,hm,           &
     			hnew,p, p0,res,recssi,recssj,xl, xref,yref
    integer :: incom, nx, ny, ndx, ndy
    
    integer :: i, im1, ip1, j, jm1, jp1, nxm1, nxm2, nym1, nym2
    real(8) :: hfim, him, pim, pim0, delyj, delyjm1, hfjm, hjm, qn1, qn2
    real(8) :: delxi, delxim1, pjm, pjm0, dpi, dmi, dpj, dmj
    real(8) :: fpi, fmi, fpj, fmj, qni, qnj
    real(8) :: qnim, qnip, qnjm, qnjp, xmxp, ymyp

    dimension p(nx,ny)
    dimension bearx(ndx,ndy),beary(ndx,ndy),         &
     	qni(nx,ny),qnj(nx,ny),h(ndx,ndy),hnew(ndx,ndy),     &
     	res(nx,ny),recssi(ndx,ndy),recssj(ndx,ndy),xref(ndx), &
     	yref(ndy),cohimx(ndx,ndy),cohjmx(ndx,ndy),              &
     	himax(ndx,ndy),himin(ndx,ndy),hjmax(ndx,ndy),           &
     	hjmin(ndx,ndy)

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
	        pim0=(p(i,j)+p(im1,j))/2.d0
	        hfim=(hnew(i,j)+hnew(im1,j))/2.d0
	        him=himax(i,j)+hfim
	        pim=pim0*him

	        call flow(pim,1.d0,1.d0,1.d0,qn1)

	        qn1=qn1*him**3
	        him=himin(i,j)+hfim
	        pim=pim0*him

	        call flow(pim,1.d0,1.d0,1.d0,qn2)

	        qn2=qn2*him**3
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
	        hfjm=(hnew(i,j)+hnew(i,jm1))/2.d0
	        hjm=hjmax(i,j)+hfjm
	        pjm=pjm0*hjm	

	        call flow(pjm,1.d0,1.d0,1.d0,qn1)

	        qn1=qn1*hjm**3
	        hjm=hjmin(i,j)+hfjm
	        pjm=pjm0*hjm

	        call flow(pjm,1.d0,1.d0,1.d0,qn2)

	        qn2=qn2*hjm**3
	        qnj(i,j)=qn1*cohjmx(i,j)+qn2*(1.d0-cohjmx(i,j))
	    enddo
    enddo

    res(1,1)=0.d0

    do j=2,nym1
	    jm1=j-1
	    jp1=j+1
	    delyj=yref(jp1)-yref(j)
	    delyjm1=yref(j)-yref(jm1)
	    ymyp=(yref(jp1)-yref(jm1))/2.d0
	    do i=2,nxm1
	        im1=i-1
	        ip1=i+1
	        delxi=xref(ip1)-xref(i)
	        delxim1=xref(i)-xref(im1)
	        xmxp=(xref(ip1)-xref(im1))/2.d0
	        qnim=qni(i,j)*f2p*(p(i,j)-p(im1,j))
	        qnip=qni(ip1,j)*f2p*(p(ip1,j)-p(i,j))
	        qnjm=qnj(i,j)*f2p*(p(i,j)-p(i,jm1))
	        qnjp=qnj(i,jp1)*f2p*(p(i,jp1)-p(i,j))

	        if(incom.eq.0) then
	            qnim=qnim*(p(im1,j)+p(i,j))/2.d0
	            qnip=qnip*(p(ip1,j)+p(i,j))/2.d0
	            qnjm=qnjm*(p(i,jm1)+p(i,j))/2.d0
	            qnjp=qnjp*(p(i,jp1)+p(i,j))/2.d0
	        endif

	        dpi=qnip/delxi
	        dmi=qnim/delxim1
	        dpj=qnjp/delyj
	        dmj=qnjm/delyjm1

	        fpi=(bearx(ip1,j)+bearx(i,j)) / 4.d0*(recssi(ip1,j)+           &
     	        (hnew(i,j)+hnew(ip1,j))/2.d0)*(p(ip1,j)+p(i,j))/2.d0

	        fmi=(bearx(im1,j)+bearx(i,j)) / 4.d0*(recssi(i,j) +            &
     	        (hnew(i,j)+hnew(im1,j))/2.d0)*(p(im1,j)+p(i,j))/2.d0

	        fpj=(beary(i,jp1)+beary(i,j))/4.d0*(recssj(i,jp1)+(hnew(i,j) + &
     	        hnew(i,jp1))/2.d0)*(p(i,jp1)+p(i,j))/2.d0

	        fmj=(beary(i,jm1)+beary(i,j))/4.d0*(recssj(i,j)+               &
     	        (hnew(i,j)+hnew(i,jm1))/2.d0)*(p(i,jm1)+p(i,j))/2.d0

            !qni and qnj form mass flux vector, res is stream function
	        qni(i,j)=((fmi-dmi)+(fpi-dpi))/2.d0
	        qnj(i,j)=((fmj-dmj)+(fpj-dpj))/2.d0

            !define left edge
	        if(j.eq.2) res(i,jm1)=res(im1,jm1)-(fmj-dmj)*xmxp

            !define bottom edge
	        if(i.eq.2) res(im1,j)=res(im1,jm1)+(fmi-dmi)*ymyp

	        res(i,j)=res(im1,j)-(fpj-dpj)*xmxp

	    enddo
    enddo

    do i=2,nxm1
	    qni(i,1) = qni(i,2)
	    qni(i,ny)= qni(i,nym1)
	    qnj(i,1) = qnj(i,2)
	    qnj(i,ny)= qnj(i,nym1)
    enddo

    do j=1,ny
	    qni(1,j) = qni(2,j)
	    qnj(1,j) = qnj(2,j)
	    qni(nx,j)= qni(nxm1,j)
	    qnj(nx,j)= qnj(nxm1,j)
    enddo

    return
    end    