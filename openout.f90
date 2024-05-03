!==============================================
subroutine openout(num, fname, nrec)
!==============================================
! this is wrapper function for opening an file
! for  output.	The record length restriction is
! system dependent.
! 
! Precondition:
!   num:	file number
!	fname:   12 character long file name, if string literal is passed, make sure it's 12 char.
!	nrec:	record length
! Postcondition:
!	file opened for formatted output with nrec record length	

    implicit none
    
    integer :: num, nrec
    character*(*) fname


! for PC(DOS/Windows, Watcom compiler) must use
!     open(num, file = fname, status = 'unknown', 
!    &     form = 'formatted', recl = nrec, err = 999)

! for IBM RS/6000, early version must use
!      open(num, file = fname, status = 'unknown',
!     &     form = 'formatted', err = 999)

! for DEC Alpha and IBM new vesion, both formats are OK.

    open(num, file = fname, status = 'unknown', form = 'formatted', recl = nrec, err = 999)
    
    return

999	write(*,*)'Error opening file ', fname
	stop
	
end subroutine openout


!==============================================
subroutine openin(num, fname, nrec)
!==============================================
! this is wrap up function for opening an file
! for  input.	The record length restriction is
! system dependent.
!     Precondition:
!		num	file number
!		fname   12 character long file name, if string literal is used, make sure it's 12 char.
!		nrec	record length
!     Postcondition:
!		file opened for unformatted input
!		with nrec record length	

    implicit none

    integer :: num, nrec
    character*12 :: fname


! for PC(DOS/Windows, Watcom compiler) must use
!     open(num, file = fname, status = 'old', 
!    &     recl = nrec, err = 999)

! for Unix, must use
!     open(num, file = fname, status = 'old',
!    &     err = 999)

    open(num, file = fname, status = 'old', err = 999)

    return

999 write(*,*)'Error opening file ', fname
	stop
	
end subroutine openin
