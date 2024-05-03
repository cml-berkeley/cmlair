!=========================================
	subroutine initResultDat
!=========================================
! write some initial info in result.dat

    use Q4_globals
    implicit none
    
    integer :: i

    write(numRes, *)'CML VERSION ', VERSION,' RESULT.DAT'
    write(numRes, *)'REPORT BUGS TO CML: bugs@cml.me.berkeley.edu'
    write(numRes, *)'NUMBER OF DISK RADII: ',irad
    write(numRes, '(a17, 20f8.4)') '   RADII(MM)   :',(radii(i)*1000, i = 1, irad)
    write(numRes, '(a17, 20f8.3)') '   SKEWS(DEG)  :',(skews(i), i = 1, irad)
    write(numRes, *)'NUMBER OF RPMS:       ',irpm
    write(numRes, '(a17, 20f9.2)') '   RPMS	   :',(rpms(i), i = 1, irpm)
    write(numRes, *)'NUMBER OF ALTITUDES:  ',ialt
    write(numRes, '(a17, 20f8.2)') '   ALTITUDES(M):',(alts(i), i = 1, ialt)

    write(numRes,*)
    write(numRes, '(a19, g12.6)') 'LENGTH(MM)     : ',xl*1000 
    write(numRes, '(a19, g12.6)') 'WIDTH(MM)      : ',yl*xl*1000
    write(numRes, '(a19, g12.6)') 'THICKNESS(MM)  : ',zl*xl*1000
    write(numRes, '(a19, g12.6)') 'TAPER LEN(MM)  : ',xt*xl*1000
    if (xt.ne.0)then
        write(numRes, '(a19, g12.6)') 'TAPER ANG(MRAD): ',ht*hm/xt/xl*1000.d0
    else
        write(numRes, *) 'TAPER ANG(MRAD): 0.0'
    endif
    write(numRes, '(a19, g12.6)') 'LOAD(G)        : ',f0*1000


    write(numRes,*)
    write(numRes, *)'SENSITIVITY CASE IDENTIFIERS:'
    write(numRes, *) '   CROWN  : -1,+1	CAMBER: -2,+2	 TWIST: -3,+3'
    write(numRes, *) '   LOAD   : -4,+4  PTORQUE: -5,+5  RTORQUE: -6,+6'
    write(numRes, *) '   TAPER-L: -7,+7  TAPER-A: -8,+8	RECESS: -9,+9'
    write(numRes, '(//)')

	return
	end subroutine initResultDat
	

!============================================
    subroutine outputGrid
!============================================
!write out x.dat, y.dat, geom.dat (added in 4.0.17), xm.dat, ym.dat
    use Q4_globals
    implicit none
    
    integer :: i, j
    logical, parameter :: outputMultilevelGrid = .false.
    
    call openout(25, 'x.dat       ', 23*nmx)
    write(25,'(10000e23.15)')(xref(i),i=1,nx)
    close(25)

    call openout(25, 'y.dat       ', 23*nmx)
    write(25,'(10000e23.15)')(yref(i),i=1,ny)
    close(25)      

    if (outputMultilevelGrid .eq. .true.) then
        call openout(25, 'x2.dat       ', 23*nmx)
        write(25,'(10000e23.15)')(xref1(i),i=1,nx1)
        close(25)
        call openout(25, 'y2.dat       ', 23*nmx)
        write(25,'(10000e23.15)')(yref1(i),i=1,ny1)
        close(25)    
        
        call openout(25, 'x3.dat       ', 23*nmx)
        write(25,'(10000e23.15)')(xref2(i),i=1,nx2)
        close(25)
        call openout(25, 'y3.dat       ', 23*nmx)
        write(25,'(10000e23.15)')(yref2(i),i=1,ny2)
        close(25)    

        call openout(25, 'x4.dat       ', 23*nmx)
        write(25,'(10000e23.15)')(xref3(i),i=1,nx3)
        close(25)
        call openout(25, 'y4.dat       ', 23*nmx)
        write(25,'(10000e23.15)')(yref3(i),i=1,ny3)
        close(25)    
        
        call openout(25, 'x5.dat       ', 23*nmx)
        write(25,'(10000e23.15)')(xref4(i),i=1,nx4)
        close(25)
        call openout(25, 'y5.dat       ', 23*nmx)
        write(25,'(10000e23.15)')(yref4(i),i=1,ny4)
        close(25)                
    endif
        
    if(isave.ne.0 .or. irailgeom.eq.ON) then
        call openout(25,'geom.dat     ', 23*nmx)
        do j=1,ny
            write(25,'(10000e23.15)') ((h(i,j)*(-1.0*hm)), i = 1,nx)
        enddo
        close(25)
    endif

    if(isave.ne.0 .or. isolv.eq.OFF) then
        !we use a shifted grid for mass flow
        call openout(25, 'xm.dat      ', 23*nmx)
        write(25,'(10000e23.15)') ((xref(i) + xref(i+1))*0.5d0, i = 1,nx-1)
        close(25)

        call openout(25, 'ym.dat      ', 23*nmx)
        write(25,'(10000e23.15)') ((yref(i) + yref(i+1))*0.5d0 ,i=1,ny-1)
        close(25)
    endif

	return
    end subroutine outputGrid
	
	
!======================================
    subroutine outputPressures(ifile)
!======================================
! ifile is the case index   

    use Q4_globals
    use ShearArrays
    use TemperatureHumidity
    implicit none

    character*12 fnpr, fncp, fnmf, fvdw, xdat, ydat, xmdat, ymdat
    character*13 cShearX, cShearY, pShearX, pShearY
    character*7 fn 
    integer :: i, j, i1, i2, ifile

    i1 = mod(ifile, 10)	
    i2 = mod(ifile / 10, 10)

    fn = char(48+i2) // char(48+i1) // '.dat'
    fnpr = 'press' // fn 
    fncp = 'cprss' // fn
    fnmf = 'mflow' // fn 
    fvdw = 'imf' // fn
    xdat = 'x' // fn
    ydat = 'y' // fn
    xmdat = 'xm' // fn
    ymdat = 'ym' // fn
    cShearX = 'CShearX' // fn
    cShearY = 'CShearY' // fn
    pShearX = 'PShearX' // fn
    pShearY = 'PShearY' // fn
    

    !adjust output pressure for humidity if necessary
    call openout(23, fnpr, 23*nmx) !pressure    
    if (doHumidity .eq. .true.) then
	    call SavePressure
	    call ModifyHumidPress
        do j=1,ny
            write(23,'(10000e23.15)')(p(i,j)-1.d0,i=1,nx)
        enddo
        call RestorePressure
    else
        do j=1,ny
            write(23,'(10000e23.15)')(p(i,j)-1.d0,i=1,nx)
        enddo
    endif
    close(23)

    
    call openout(23, fncp, 23*nmx) !contact pressure
    do j=1,ny
        write(23,'(10000e23.15)')(cp(i,j),i=1,nx)
    enddo
    close(23)

    call openout(23, fnmf, 23*nmx) !massflow
    do j=1, ny-1
        write(23,'(10000e23.15)')(res(i,j), i = 1, nx-1)
    enddo
    close(23)

    call openout(23, fvdw, 23*nmx) !MolecularForces
    do j=1, ny
        write(23,'(10000e23.15)')(vdwMolecularForceMap(i,j), i = 1, nx)
    enddo
    close(23)
    
    if (iOutputShear .ne. 0) then
        call openout(23, pShearX, 23*nmx)
        do j=1, ny
            write(23,'(10000e23.15)')(PoisilleShearX(i,j), i = 1, nx-1)
        enddo
        close(23)    

        call openout(23, pShearY, 23*nmx)
        do j=1, ny-1
            write(23,'(10000e23.15)')(PoisilleShearY(i,j), i = 1, nx)
        enddo
        close(23) 
        
        call openout(23, cShearX, 23*nmx)
        do j=1, ny
            write(23,'(10000e23.15)')(CouetteShearX(i,j), i = 1, nx-1)
        enddo
        close(23)    

        call openout(23, cShearY, 23*nmx)
        do j=1, ny-1
            write(23,'(10000e23.15)')(CouetteShearY(i,j), i = 1, nx)
        enddo
        close(23) 
    endif   

    !the following was added to accomodate an adaptive grid that adapts at each radial position
    call openout(23, xdat, 23*nmx) !x0x.dat
    write(23,'(10000e23.15)') (xref(i), i = 1,nx)
    close(23)
    
    call openout(23, ydat, 23*nmx) !y0x.dat
    write(23,'(10000e23.15)') (yref(i), i = 1,ny)
    close(23)    

    if(isave.ne.0 .or. isolv.eq.OFF) then
        !we use a shifted grid for mass flow
        call openout(23, xmdat, 23*nmx)
        write(23,'(10000e23.15)') ((xref(i) + xref(i+1))*0.5d0, i = 1,nx-1)
        close(23)

        call openout(23, ymdat, 23*nmx)
        write(23,'(10000e23.15)') ((yref(i) + yref(i+1))*0.5d0 ,i=1,ny-1)
        close(23)
    endif

    return
    end subroutine outputPressures


!=====================================
    subroutine outputResult(nf)
!=====================================
! write the result to file with file number nf.
 
    use Q4_globals
    implicit none
	
	integer :: nf, i
	real(8) :: tempHeight
	
	if (isolv.eq.0 .or. isolv.eq.1) call CoForce()
	
	write(numRes, '(///)')
	write(numRes, '(4(a11,i2,a2))' )    'RADIUS NO.' ,icurrad, '  ', 'RPM    NO.' ,icurrpm, '  ', &
                                        'ALTIT. NO.' ,icuralt, '  ', 'SENSI. NO.' ,icursen
    
    write(numRes,*)
    
    call display(numRes)
	
	write(numRes,'(a20,g12.6)') 'POSITIVE FORCE(G): ', fpos*1000
    write(numRes,'(a20,g12.6)') 'NEGATIVE FORCE(G): ', fneg*1000 
    write(numRes,'(a20,g12.6)') 'CONTACT  FORCE(G): ', fcr*1000
    write(numRes,'(a25,g12.6)') 'VAN DER WALLS FORCE(G): ', fvdw_output*1000
	write(numRes,'(a20,g12.6)') 'X-SHEAR  FORCE(G): ', fsp*1000/zl
    write(numRes,'(a20,g12.6)') 'Y-SHEAR  FORCE(G): ', fsr*1000/zl 
    write(numRes,'(a18,g12.6)') 'Z-MOMENT (uN-M): ', Zmom 
	write(numRes,'(a33,g12.6)') 'X-CENTER OF POSITIVE FORCE(mm): ', xPosLoc  
	write(numRes,'(a33,g12.6)') 'X-CENTER OF NEGATIVE FORCE(mm): ', xNegLoc
	write(numRes,'(a33,g12.6)') 'Y-CENTER OF POSITIVE FORCE(mm): ', yPosLoc
	write(numRes,'(a33,g12.6)') 'Y-CENTER OF NEGATIVE FORCE(mm): ', yNegLoc
	!write out the rest of the POIs
	if (numPOI .gt. 4) then
	    write(numRes,'(a31,g12.6)') 'REMAINING POINTS OF INTEREST: ', numPOI - 4
	    do i = 5, numPOI
            
            tempHeight = hintNew(i)*hm*1d9
            !you know something?  I really HATE fortran format statements
            !Height( X.XXX, X.XXX): XXXX.XXXXXXX
	        write(numRes, '(a13,2(f5.3,a2),a3,g12.6)') 'POI Height( ',xintNew(i)*xl*1000,',',yintNew(i)*xl*1000,' )', ': ', tempHeight
	    enddo
	endif
	    
	return
	end subroutine outputResult



!=====================================
    subroutine display(numfile)
!=====================================
    use Q4_globals
    implicit none
    
    integer :: numfile
    real(8) :: blong
    integer :: i

!      write(numfile,'(4a15)')
!     &                      '   ERROR      ', ' NOMINAL HM(NM)',
!     &                      '  PITCH(URAD) ', '  ROLL(URAD)  '
!	write(numfile,'(e12.4,a3,3(f12.4,a3))')

    write(numfile,'(4a16)') '      ERROR     ', &
                            '  NOMINAL HM(NM)', &
                            '   PITCH(URAD)  ', &
                            '   ROLL(URAD)   '
    write(numfile,'(4G16.8)') err,hmin*hm*1d9,-hx0*hm/xl*1d6,hy*hm/xl*1d6
    write(numfile,*)

    blong = 0
    do i = 1,4
        if (xintNew(i)*xl*1000.ge.10) blong = 1
        if (yintNew(i)*yl*1000.ge.10) blong = 1
    enddo

    if (blong.eq.0) then
        write(numfile,'(4(a3,2(f5.3,a1)))') ('H(',xintNew(i)*xl*1000,',',yintNew(i)*xl*1000,')',i=1,4)
    else
        write(numfile,'(4(a3,2(f5.2,a1)))') ('H(',xintNew(i)*xl*1000,',',yintNew(i)*xl*1000,')',i=1,4)
    endif

 	write(numfile,'(4(g12.6,a3))') (hintNew(i)*hm*1d9, '   ', i=1,4)
    write(numfile,*)

    write(numfile,'(a15,g12.6,a10,f6.3,a2,f6.3,a1)') &
                ' MIN. HEIGHT = ',MinFH*hm*1d9,' (NM) AT (',MinFHLocX*xl*1000,  &
                ', ', MinFHLocY*xl*1000, ')'

    return
    end subroutine display
    

!=====================================
    subroutine output_spacing
    ! Jinglin: this function is added to automatically output the spacing at current grid.
    ! currently not used by the code.
!=====================================
    use Q4_globals
    implicit none
    
    integer :: i, j, i1, i2
    character*24 filename
	character*7 fn

    do icurrad = 1, irad
	    call setRadius(radii(icurrad), skews(icurrad))
        call gethx  
        i1 = mod(icurrad, 10)	
	    i2 = mod(icurrad / 10, 10)
        fn = char(48+i2) // char(48+i1) // '.dat'
	    filename = 'FHSpacingGrid' // fn   
        call openout(12,filename,23*nmx)
        do j=1,ny
            write(12,'(10000e23.15)')( (hnew(i,j)+h(i,j))*hm*1.d09,i=1,nx)
        enddo
        close(12)
        filename = 'xx' // fn
        call openout(12, filename, 23*nmx)
        write(12, '(10000e23.15)')(xref(i),i=1,nx)
        close(12)
        filename = 'yy' // fn
        call openout(12, filename, 23*nmx)
        write(12,'(10000e23.15)')(yref(i),i=1,ny)
        close(12)
    enddo
     
    return
     
    end subroutine output_spacing
    
    
!==============================================================
    subroutine output_spacing_2
!==============================================================
! Dolf: this subroutine added to automatically save 
! the output spacing into a file.

   
   use Q4_globals 
   implicit none 
!  already defined in Q4_Globals module..   
!  real(8), dimension(:,:), allocatable :: sliderdiskspacing
   real(8):: pitch1,roll1,xv1,yv1,resessdif
   real(8):: NomFH1
   real(8):: deltaX,deltaY,deltaZ,deltaZplus,recess
   integer :: i, j, i1, i2
   character*24 filename
   character*7 fn
   
   allocate (sliderdiskspacing(nx+1,ny+1))
	

    do icurrad = 1, irad
	   call setRadius(radii(icurrad), skews(icurrad))
      
        i1 = mod(icurrad, 10)	
	    i2 = mod(icurrad / 10, 10)
        fn = char(48+i2) // char(48+i1) // '.dat'
	    filename = 'FHSpacingGrid' // fn   
        call openout(15,filename,23*nmx)
        
        xv1=0.0
        yv1=0.0
        pitch1 = -hx0*hm/xl*1d6
        roll1 = hy*hm/xl*1d6
        NomFH1 = hmin*hm*1d9
        
     do j=1, ny
     
     yv1=yref(j)
     
     do i=1, nx
     
     xv1=xref(i)
     
	 call pointrecess(recess,xv1,yv1)
     deltaX = 0.0
     deltaY = 0.0
!     sliderdiskspacing(xref(i),yref(j)) = 0.0
     sliderdiskspacing(i,j) = 0.0
     deltaX = h0+hx0*xref(i)
     deltaY = hy*(yref(j)-yl*0.5d0)
     sliderdiskspacing(i,j) = recess+deltaX-deltaY
     
     
     
!     write (15,'(1000G18.12)')((sliderdiskspacing(i,j)*hm*1d9),i=1,nx)
      write (15,'(1000G18.12)') xref(i),yref(j)*xl,sliderdiskspacing(i,j)*hm*1d9
   
   enddo !i=1, nx
        
   enddo   !j=1, ny
   close (15)
        
        filename = 'xx' // fn
        call openout(12, filename, 23*nmx)
        write(12, '(10000e23.15)')(xref(i),i=1,nx)
        close(12)
        filename = 'yy' // fn
        call openout(12, filename, 23*nmx)
        write(12,'(10000e23.15)')(yref(i),i=1,ny)
        close(12)
   
   enddo   !icurrad
   !deallocate (sliderdiskspacing(nx,ny))
     
   return
   end subroutine output_spacing_2
         
   
!==============================================================
    subroutine OutputCrashedResult(nf)
!==============================================================
    use Q4_globals
    implicit none
    integer :: i, nf
    
    err = -1
    hmin = -1
    hx0 = -1
    hy = -1
    fpos = -1
    fneg = -1
    fcr = -1
    fvdw_output = -1
	fsp = -1
    fsr = -1
    Zmom = -1
	xPosLoc = -1
	xNegLoc = -1
	yPosLoc = -1
	yNegLoc = -1
    do i = 1, numPOI
        hintNew(i) = -1.d0 * hm * 1d-9
    enddo
    call outputResult(nf)

    return      
	end subroutine OutputCrashedResult
	
	
!==============================================================
    subroutine outputArray(theArray, sizeX, sizeY, outFileName)
!==============================================================
    use Q4_globals
    implicit none
    
    character*(*) :: outFileName
    real(8) :: theArray
    integer :: sizeX, sizeY, i, j
    dimension :: theArray(sizeX, sizeY)
    
    call openout(25, 'x.dat       ', 23*nmx)
    write(25,'(10000e23.15)')(xref(i),i=1,nx)
    close(25)

    call openout(25, 'y.dat       ', 23*nmx)
    write(25,'(10000e23.15)')(yref(i),i=1,ny)
    close(25)   

    call openout(866, outFileName, 23*1000)
    do j=1,sizeY
        write(866, '(10000e23.15)') (theArray(i,j), i=1, sizeX)
    enddo
    close(866)
    
    return
    end
	


!=====================================
    subroutine DisplayAttitude
!=====================================
    use Q4_globals
    implicit none

    real(8) :: blong
    integer :: i

    write(*,'(3a16)') '  NOMINAL HM(NM)', &
                            '   PITCH(URAD)  ', &
                            '   ROLL(URAD)   '
    write(*,'(3G16.8)') hmin*hm*1d9,-hx0*hm/xl*1d6,hy*hm/xl*1d6

    return
    end subroutine DisplayAttitude
    
    	