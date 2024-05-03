!===================================================================
subroutine GetOptionalParameters(iFile)
!===================================================================
    use OptionalParams
    implicit none

    character *PARAMSIZE :: strParamNames(MAXPARAMS)
    character *PARAMSIZE :: strParamValues(MAXPARAMS)
    character *LINESIZE :: curLine
    integer :: iFile
    real(8) :: paramValues(MAXPARAMS)
    integer :: numNames, numValues, i
    
    integer :: indexVal
    
    !first see if we have an optional parameters section
    read(iFile, '(a1024)', end=10, err=10) curLine
    if (INDEX(curLine, 'End of Input Data') .ne. 0) goto 10
    if (INDEX(curLine, 'omments') .ne. 0) goto 10  !look for 'comments' but disregard capitalization
    if (INDEX(curLine, 'dditional') .eq. 0 .and. &  !to satisfy a few members, make sure we also write "Additional Parameters" (don't worry about caps)
        INDEX(curLine, 'arameters') .eq. 0) goto 10
    !now read in all optional parameters
    do while (1)
        read(iFile, '(a1024)', end=10, err=10) curLine
        
        indexVal = INDEX(curLine, '****') !make sure we're not at the end of the optional parameters section
        if (indexVal .ne. 0) goto 10
        
        call Split(curLine, strParamNames, numNames)
        read(iFile, '(a1024)') curLine              !read in parameter values
        call Split(curLine, strParamValues, numValues)
        
        if (numNames .ne. numValues) then
            write(*,*) 'error encountered when reading parameters:'
            write(*,*) 'parameter header does not match parameter values'
            goto 10
        endif
        
        call StringsToReal(strParamValues, paramValues, numValues)
        
        do i = 1, numValues
            call ProcessParameter(strParamNames(i), paramValues(i))
        enddo
    enddo

10  continue

    return
end subroutine GetOptionalParameters


!===================================================================
subroutine ClearString(str, size)
!===================================================================
    implicit none
    character *(*) :: str
    integer :: i, size
        
    do i = 1, size
        str(i:i) = ' '
    enddo
end subroutine ClearString

!===================================================================
subroutine Split(inputString, outputStrings, numElements)
!===================================================================
    use OptionalParams 
    implicit none
    
    character *PARAMSIZE :: outputStrings(MAXPARAMS), tempString
    character *LINESIZE :: inputString
    integer :: i, j, numElements
    logical :: IsWhiteSpace
    
    numElements = 0
    i = 1
    !first clear output strings
    do j = 1, MAXPARAMS
        call ClearString(outputStrings(i), PARAMSIZE)
    enddo
    
    do while (i .lt. LINESIZE)
    
        !find first non-space/tab character
        do while (IsWhiteSpace(inputString(i:i)) .and. i .lt. LINESIZE)
            i = i+1
        enddo
        if (i .ge. LINESIZE) goto 20
        
        !copy characters from current word to outputStrings until we get to a space or tab
        j = 1
        do while (IsWhiteSpace(inputString(i:i)) .eq. .false. .and. i .lt. LINESIZE .and. j .lt. PARAMSIZE)
            !outputStrings(numElements) = inputString(i)
            tempString(j:j) = inputString(i:i)
            i = i+1
            j = j+1
        enddo
        numElements = numElements+1
        outputStrings(numElements) = tempString(1:j-1)
        
        !gone over, not a necessary statement...
        if (i .ge. LINESIZE) goto 20
        
    enddo !repeat until end of inputString
    
    
20  continue

    return
end subroutine Split

!===================================================================
logical function IsWhiteSpace(targetChar)
!===================================================================
    use OptionalParams
    implicit none
    
    integer :: ichar
    character :: targetChar
        
    if (targetChar .eq. ' ' .or. ichar(targetChar) .eq. 9) then !if space or tab
        IsWhiteSpace = .true.
    else
        IsWhiteSpace = .false.
    endif
    
    return
end

!===================================================================
subroutine StringsToReal(inputStrings, outputReals, numElements)
!===================================================================
    use OptionalParams
    implicit none
    
    character *PARAMSIZE :: inputStrings(MAXPARAMS)
    REAL(8) :: tempReal, outputReals(MAXPARAMS)
    integer :: i, numElements
    
    do i = 1, numElements
        tempReal = dnum(inputStrings(i))
        outputReals(i) = tempReal
    enddo
    
    return
end subroutine StringsToReal