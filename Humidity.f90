!==========================================
	subroutine ModifyHumidPress
!==========================================
! Modifies the pressure to account for effects of water condensation due to humidity temperature and 
! high pressures see "Effects of Humid Air on Air-Bearing Flying Height" by Strom, Zhang, Lee, 
! Khurshudov and Tyndall.  This appeared in IEEE Transactions on Magnetics, Vol. 43, No. 7, July 2007.
    use Q4_globals
    use TemperatureHumidity
    implicit none

    integer :: i, j
    real(8) :: satVaporPressure, partialPressWater, pressWater, pressureLoss, pressureLossAtm
    
    satVaporPressure = 611.d0 * dexp( (17.5 * temperature) / (241.0 + temperature) )
    partialPressWater = (humidity/100.d0) * satVaporPressure
    
#ifdef _DEBUG
    call outputArray(savedPressure, nx, ny, 'pOrig')
#endif !#_DEBUG
    
    do j = 2, ny-1
        do i = 2, nx-1
            pressWater = p(i, j) * partialPressWater
            pressureLoss = 0.d0
            if (pressWater .gt. satVaporPressure) then
                pressureLoss = pressWater - satVaporPressure
            endif
            pressureLossAtm = pressureLoss / p0
            p(i, j) = p(i, j) - pressureLossAtm
        enddo
    enddo
    
#ifdef _DEBUG
    call outputArray(p, nx, ny, 'pHumid')
#endif !#_DEBUG
    
	return
	end subroutine ModifyHumidPress
	
!==========================================
	subroutine SavePressure
!==========================================
    !simply copies p to savedPressure
    use Q4_globals
    use TemperatureHumidity
    implicit none
    
    integer :: i, j
    
    do j = 2, ny-1
        do i = 2, nx-1
            savedPressure(i, j) = p(i, j)
        enddo
    enddo
    
	return
	end subroutine SavePressure


!==========================================
	subroutine RestorePressure
!==========================================
    !simply copies savedPressure back to p
    use Q4_globals
    use TemperatureHumidity
    implicit none
   
    integer :: i, j
    
    do j = 2, ny-1
        do i = 2, nx-1
            p(i, j) = savedPressure(i, j)
        enddo
    enddo
	
	return
	end subroutine RestorePressure
