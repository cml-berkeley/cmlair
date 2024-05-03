!   allow for a user defined geometry matrix for the slider
!   instead of using parabolic camber crown and twist
!======================================
    subroutine ReadUDG()
!======================================
    use Q4_globals
    use UserDefinedGeometry
    implicit none

    integer :: i, j, nSizeXTemp, nSizeYTemp, iFile
    real(8) :: dLowerXTemp, dUpperXTemp, dLowerYTemp, dUpperYTemp
    !used for parsing the input file
    character *500 header !holds parameter header
    character *500 words(10) !holds header split into individual words
    integer nwords
    integer :: getNumUDGFiles

    !first determine how many UDG files there are
    !integer :: numUDGAreas
    !numUDGAreas = getNumUDGAreas()
    !allocate the UDG areas array
    !allocate arrays to hold the fXLower, fXUpper, fYLower, fYUpper and range values
    !read in each UDG file and area data and store it in the array


    character*15 UDGFileName
    character*7 fn 
    integer :: i1, i2
    
    numUDGRegions = getNumUDGFiles()
    if (numUDGRegions .eq. 0) return
    
    allocate (UserDefGeom(numUDGRegions))

    do ifile = 1, numUDGRegions
        !open the files one by one
        i1 = mod(ifile, 10)    
        i2 = mod(ifile / 10, 10)
        fn = char(48+i2) // char(48+i1) // '.dat'
        UDGFileName = 'Usergeom' // fn
        
        OPEN(UNIT=66, ERR=555, FILE=UDGFileName, STATUS='OLD')

        
        !read in the header info.  Handle the old style gracefully
        read(66,'(a500)') header
        rewind(66)
        call parse(header, words, nwords)
        if (nwords .eq. 4) then
            read(66, *) dLowerXTemp, dUpperXTemp, dLowerYTemp, dUpperYTemp
            read(66, *) nSizeXTemp, nSizeYTemp
        else if (nwords .eq. 2) then
            read(66, *) nSizeXTemp, nSizeYTemp
            dLowerXTemp = 0.0
            dUpperXTemp = 1.0
            dLowerYTemp = 0.0
            dUpperYTemp = yl
        else
            write(*,*) 'Error encountered when reading in User Defined Geometry Files'
            write(*,*) 'The error occoured on file ', ifile, 'of ', numUDGRegions
            numUDGRegions = ifile-1
            goto 555
        endif
        
        !normalize lengths and store the data
        UserDefGeom(iFile)%dLowerX = dLowerXTemp / (xl * 1000.0)
        UserDefGeom(iFile)%dUpperX = dUpperXTemp / (xl * 1000.0)
        UserDefGeom(iFile)%dLowerY = dLowerYTemp / (xl * 1000.0)
        UserDefGeom(iFile)%dUpperY = dUpperYTemp / (xl * 1000.0)
        UserDefGeom(iFile)%nSizeX = nSizeXTemp
        UserDefGeom(iFile)%nSizeY = nSizeYTemp

        !now that we have the data on the matrix, allocate space for it and read it in


    !you will need to enable the preprocessor in the project settings
    !TO COMPILE UNDER COMPAQ VISUAL FORTRAN UNCOMMENT THE NEXT LINE
!#define COMPAQ_VISUAL_FORTRAN__
#ifndef COMPAQ_VISUAL_FORTRAN__
        allocate (UserDefGeom(iFile)%UDGMatrix (UserDefGeom(iFile)%nSizeX, UserDefGeom(iFile)%nSizeY) )
#endif

        do i = 1, UserDefGeom(iFile)%nSizeX
            READ(66,*)(UserDefGeom(iFile)%UDGMatrix(i,j), j=1, UserDefGeom(iFile)%nSizeY)
        enddo
        CLOSE(66)
    enddo

555 continue

    return
    end subroutine ReadUDG
    
    
    
!======================================
    integer function getNumUDGFiles()
!======================================    
    implicit none
    
    character*15 UDGFileName
    character*7 fn 
    integer :: i1, i2, ifile
    
    getNumUDGFiles = 0
    iFile = 1
    do while (.true.)
        
        i1 = mod(iFile, 10)    
        i2 = mod(iFile / 10, 10)
        fn = char(48+i2) // char(48+i1) // '.dat'
        UDGFileName = 'Usergeom' // fn
        
        OPEN(UNIT=66, ERR=556, FILE=UDGFileName, STATUS='OLD')
        
        getNumUDGFiles = getNumUDGFiles + 1
        ifile = iFile + 1
        
        CLOSE(66)
    enddo
    
556 continue

    return
    end
    

!     determine height at point (xv, yv) on user
!     defined geometry matrix using 2-d interp
!    Added by BC 9/03
!=======================================================
    real(8) function UDGHeight(xv, yv)
!=======================================================
    use Q4_globals
    use UserDefinedGeometry
    implicit none

    INTEGER i, j
    real(8) xv, yv
    real(8) ztmp(3)
    real(8) dx, dy, QuadraticInterp
    integer :: xfloor, yfloor
    real(8) :: dLowerX , dUpperX, dLowerY, dUpperY, dLengthX, dLengthY


!c    QUADRATIC INTERPOLATION OF UDG

    UDGHeight = 0
    do i = 1, numUDGRegions
        dLowerX = UserDefGeom(i)%dLowerX
        dUpperX = UserDefGeom(i)%dUpperX
        dLowerY = UserDefGeom(i)%dLowerY
        dUpperY = UserDefGeom(i)%dUpperY
        if (xv .ge. dLowerX .and. xv .le. dUpperX .and.   &
            yv .ge. dLowerY .and. yv .le. dUpperY) then
            
            dLengthX = UserDefGeom(i)%dUpperX - UserDefGeom(i)%dLowerX
            dLengthY = UserDefGeom(i)%dUpperY - UserDefGeom(i)%dLowerY
            
            !get the width of the space between UDG points
            dx = (dLengthX) / (UserDefGeom(i)%nSizeX-1)
            dy = (dLengthY) / (UserDefGeom(i)%nSizeY-1)
        
            
            if (xv .lt. (dUpperX - dx)) then  !if we're not in the last section of the UDG
                !way to complicated way to find the udg point that is below our point
                xfloor = idint((xv - dLowerX) * ((UserDefGeom(i)%nSizeX - 1.d0) / dLengthX))
            else 
                xfloor = UserDefGeom(i)%nSizeX - 3
            endif

            if (yv .lt. (dUpperY - dy)) then
                yfloor = idint((yv - dLowerY) * ((UserDefGeom(i)%nSizeY - 1.d0) / dLengthY))
            else
                yfloor = UserDefGeom(i)%nSizeY - 3
            endif        
        
            do j=1, 3
                ztmp(j) = QuadraticInterp(yv-dLowerY, yfloor * dy, (yfloor+1) * dy,             &
                                        UserDefGeom(i)%UDGMatrix(xfloor+j, yfloor+1),   &
                                        UserDefGeom(i)%UDGMatrix(xfloor+j, yfloor+2),   &
                                        UserDefGeom(i)%UDGMatrix(xfloor+j, yfloor+3),   &
                                        dy)
            enddo

            UDGheight = UDGheight + QuadraticInterp(xv-dLowerX, xfloor * dx,(xfloor+1) * dx,        &
                                                    ztmp(1),ztmp(2),ztmp(3), dx)
        
        endif
    enddo

    return
    end

!=======================================================
    real(8) function QuadraticInterp(x, x0, x1, y0, y1, y2, h)
!=======================================================
    real(8) :: a, b, c, x, x0, x1, y0, y1, y2, h

    a = y0
    b = (y1 - y0) / h
    c = (y2 - (2.d0*y1) + y0) / (2.d0 * (h**2))
    
    QuadraticInterp = a + b*(x - x0) + c*(x-x0)*(x-x1)
    
    return
    end