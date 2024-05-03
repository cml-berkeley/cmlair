c
c  Free string manipulation library added 12/6/04
c  Hopefully we can clear up some of the problems with parsing the input files
c



!Index
!
!1.	Introduction
!2.	Installation
!3.	Usage
!4.	Concluding Remarks
!
!
!1. Introduction
!===============
!
!This fortran77 library is developed to simplify strings handling
!and i/o parsing operation.
!In fortan the read operation of free format files (like C files) is 
!not very easy, because this language requires fields with formatted 
!position.
!On the other hand the native fortran77 does'nt have directed functions 
!to manage the strings also for simpliest operations.
!In this library you can find functions & routines to read also the 
!free format files and handle the strings.
!
!2. Installation
!===============
!
!In archive file you can find a Makefile to build a libstr.a library
!in order to use routines & functions in your program.
!Check in makefile compilation and link options for your computer.
!Default options are for SGI machines.
!If you can move the libstr.a in your lib directory, make install and 
!link your programs with -lstr option, otherwise use -L -lstr 
!options (PATH is directory who is placed libstr.a).
!
!
!3. Usage
!========
!
!3.1 Functions
!~~~~~~~~~~~~~
!
!3.1.1. Length - Return the string length without the blanks characters
!
!Synopsis:
!
!	integer function length(char)
!	character *N char
!
!Description:
!Instead of the standard len function, this routine returns the 
!real length of a string (without the blanks).
!
!3.1.2 Isnumber - Check if the string argument contain a number 
!
!Synopsis:
!
!	logical function isnumber(char)
!	character *N char
!
!Description:
!This routine returns a true value if char contain a number otherwise
!returns false.
!
!3.1.3 Valnum - Return the real value contained into string argument
!
!Synopsis:
!
!	real function valnum(char)
!	character *N char
!
!Description:
!Obviously this routine can be used also with an integer value making a 
!final casting.
!
!
!3.2 Subroutines
!~~~~~~~~~~~~~~~
!
!3.2.1 Right - return the right string portion.
!
!Synopsis:
!
!	subroutine right(input,nchar,output)
!	character *N input,output
!	integer nchar
!
!Description:
!This routine returns in output argument the right nchar portion 
!of input string.
!
!3.2.2 Union - Join two strings 
!
!Synopsis:
!
!	subroutine union(input1,input2,output)
!	character *N input1,input2,output
!
!Description:
!This routine chains two input strings in output argument.
!
!3.2.3 Uniblk - Join two strings with spacing
!
!Synopsis:
!
!	subroutine uniblk(input1,input2,nblank,output)
!	character *N input1,input2,output
!	integer nblank
!
!Description:
!This routine returns in output argument the chaining of the
!two input strings with nblank blank characters between them
!
!3.2.4 Intstr - Translate a integer value into string
!Synopsis:
!
!	subroutine intstr(int,char,nchar)
!	character *N char
!	integer nchar,int
!
!Description:
!This routine translates the int integer into char string. Nchar is
!the output string length to obtain a fixed size of output (example:
!001,002,....,023,etc). If nchar is zero, the routine generates a 
!string with the exact size of the int number 
!(example:1,2,....,3423,etc).
!
!
!3.2.5 Flostr - Translate a real*8 value into string
!Synopsis:
!
!	subroutine intstr(real,char,ndec)
!	character *N char
!	integer	ndec
!	real*8  real
!
!Description:
!This routine translates the real value into char string. The real value
!must be in double precision. Ndec value is the size of decimal part
!      
!
!3.2.6 Parse - Perform the string parsing
!
!Synopsis:
!
!	subroutine parse(char,word,nword)
!	character *N char,word(*)
!	integer nword
!
!Description:
!This routine parses a string in word defined by blank 
!characters. The routine output is a string array (word(*)) and the 
!nword value is the word number obtained from string parsing.
!
!Example:
!	character *20 char,word(5)
!        integer nword,i
!
!	char='Fortran77 is a beautiful language'
!	call parse(char,word,nword)
!	doi=1,nword 
!	write(6,*)word(i)
!	end do
!	end
!
!This program shows this output:
!	Fortran77	(word(1))
!	is		 ......
!	a		 ......
!	beautiful	 ......
!	language	(word(5))
!
!3.2.7 readf - Parse a string with specified template
!Synopsis:
!
!	subroutine readf(line,format,float,int,char,success)
!	character *N line,format,char(*)
!	integer int(*)
!	real float(*)
!        logical success
!
!Description:
!This routine searchs elements indicated into format string and puts
!them in float, int and char arrays.
!Possibilities for format elements are:
!a	string field
!i	integer field
!f	float field
!x	skip a field
!the logical success is returned false if routine readf finds a format
!error.
!
!Example:
!With readf routine you can read a PDB file with this simple program:
!
!	implicit integer (k-l)
!	parameter (maxatom=10000)
!	character *80 line
!	character *4 str(2)	! atom name and residue name
!	character *4 atmname(maxatom),resname(maxatom)
!	real flo(3)		! x,y,z coordinates
!	real x(maxatom),y(maxatom),z(maxatom)
!	integer in(2)		! atom e residue number
!	integer atmnum(maxatom),resnum(maxatom)
!	logical suc
!
!	open(10,file='file.pdb')
!	i=0
!10	read(10,'(a80)',end=99)line
!	if(line(1:4).eq.'ATOM')then
!	i=i+1
!	call readf(line,'xiaaxifff',flo,in,str,suc)
!	if(.not.suc)stop
!	X(i)=flo(1)
!	Y(i)=flo(2)
!	Z(i)=flo(3)
!	atmnum(i)=in(1)
!	resnum(i)=in(2)
!	atmname(i)=str(1)
!	resname(i)=str(2)
!	end if
!	goto 10
!99	end
!
!Note:
!The utility of this routine is very remarkeable with files written
!in free format such the Quanta/Chemnote mol files or Hyperchem hin
!files, unreadable with native fortran read instruction. 
!  
!4. Concluding Remarks
!=====================
!
!All trademarks and  softwares  directly  or  indirectly  referred  in  this
!document,  are  copyrighted from legal owners. LibString is a program freeware
!and can be spreaded trough Internet, BBS, CD-ROM and each other electronic form.
!Authors of this library accepts no responsabilty for hardware/software
!damages resulting from the use of this package.


c***************************************************************
c                           LibString 1.0			
c             (c) 1998 Giulio Vistoli & Alex Pedretti
c***************************************************************

c==========================================
      integer function length(str)
c     return the string length without the blanks characters

      implicit none
      integer :: i, lmax
      
      character *(*) str    
      
      lmax=len(str)      
       
c     search the last non blank character
      do i=lmax,1,-1
        if(str(i:i).ne.' ')then
          length=i
          return
        end if
      end do

      length=lmax

      return
      end     
      
      
c==========================================
      logical function isnumber(str)
c     check if the string argument contain a number

      implicit none
      integer :: i, length
      character *(*) str   

      isnumber=.true.
	
      do i=1,length(str)
        if((str(i:i).lt.'0'.or.str(i:i).gt.'9').and.
     $      str(i:i).ne.'.'.and.str(i:i).ne.'-'.and.
     $      str(i:i).ne.'+')then
          isnumber=.false.
          return
        end if
      end do

      return  
      end

c==========================================
      subroutine  right(str,nch,res)
c     return the right string portion

      implicit none
      
      integer :: l, nch, length
      character *(*) str,res
      
      l=length(str)
      res=str(l-nch+1:l)

      return
      end
      
c==========================================      
      subroutine union(str,add,res)
c     join two strings

      implicit none
      integer l1, length
      character *(*) str,res,add
      
      l1=length(str)
      
      res=str
      res(l1+1:l1+length(add))=add
      
      return
      end                      
      
      
      
c==========================================      
      subroutine uniblk(str,add,nbl,res)
c     join two strings with spacing

      implicit none
      
      integer :: l1, length, nbl
      character *(*) str,res,add
      
      l1=length(str)
      
      res=str
      res(l1+1+nbl:l1+length(add)+nbl)=add
      
      return
      end      
      
c==========================================
      subroutine readf(line,str,fl,in,ch,success)
c     parse a string with specified template

      implicit none

      character *(*) str,ch(*)
      character *(*) line
      character *30 word(50) 
      real fl(*),flo, valnum
      integer in(*)
      logical success, isnumber
      integer na, ni, nf, nw, nn, length
      
      success=.true.
      na=0
      ni=0
      nf=0

      call parse(line,word,nw)	

      nn=1

c     if nn field is a string
10    if(str(nn:nn).eq.'a')then
        na=na+1
        ch(na)=word(nn)

c     if nn field is a integer
      elseif(str(nn:nn).eq.'i')then
        ni=ni+1     
        success=isnumber(word(nn))
        if(.not.success)return
        flo=valnum(word(nn))
        in(ni)=int(flo)

c     if nn field is a float
      elseif(str(nn:nn).eq.'f')then
        nf=nf+1     
        success=isnumber(word(nn))
        if(.not.success)return
        fl(nf)=valnum(word(nn))
      end if

c     next field
      nn=nn+1
      if(nn.le.nw .and. nn.le.length(str))goto 10

      return
      end
      
c==========================================      
      subroutine parse(str,words,nw)
c     perform the string parsing

      implicit none

      
      integer :: nw, length, l, lf, ltot
 !     character *(*) str,words(*)
      character *(*) str,words(300)
      
      nw=0
      ltot=length(str)
      l=ltot

c     find and skip blank characters
10    if(str(1:1).ne.' ')goto 20     
      str(1:l-1)=str(2:l)
      l=l-1
      goto 10 

20    lf=index(str(1:l),' ')
      nw=nw+1

c     define nw word
      if(lf.eq.0)then
        words(nw)=str(1:l)
      return
      end if

c     define last word
      words(nw)=str(1:lf-1)
      str=str(lf+1:l)
      l=l-lf
      
      goto 10      

      end
      
c==========================================      
      real function valnum(str)
c     return the real value contained into a string

      implicit none
      
      integer :: i, lu, ifr, iin, length, k
      character *(*) str
      logical segno
      

      segno=.false.
      valnum=0.00
      lu=length(str)
      
c     check the number sign
      if(str(1:1).eq.'-')then
        segno=.true.
        str=str(2:lu)
        lu=lu-1
      end if

c     check if number is float or integer
      if(index(str,'.').ne.0)then
        iin=index(str,'.')-1
      else
        iin=lu
      end if

      ifr=lu-(iin+1)

c     translate the integer portion
      do i=1,iin  
        k=ichar(str(i:i))-48
        valnum=valnum+float(k)*10.00**float(iin-i)
      end do         

      if(iin.eq.lu)goto 10
      str=str(iin+2:lu)

c     translate the decimal portion
      do i=1,ifr  
        k=ichar(str(i:i))-48
        valnum=valnum+float(k)/10.00**float(i)
      end do

10    if(segno)valnum=-valnum

      return
      end     


c==========================================
      subroutine intstr(num,str,l)
c     translate a integer value into string

      implicit none

      integer :: lun, num, j, l, n
      character *(*)str
      character *1 cifra(10)
      logical segno

      data cifra /'0','1','2','3','4','5','6','7','8','9'/

      lun=len(str)
      if(lun.gt.30)stop
      segno=.false.

c     check the number sign
      if(num.lt.0)then
        segno=.true.
        num=abs(num)
      end if

c     translate the integer num
      do j=1,lun
        n=num/10**(lun-j)
        num=num-(n*10**(lun-j))
        str(j:j)=cifra(n+1)
      end do

c     if the str length is fixed (l)
      if(l.ne.0)then
        call right(str,l,str)
        str=str(1:l)
        return
      end if

c     else delete zero characters
      l=lun
10    if(str(1:1).ne.'0')goto 20
      str(1:l-1)=str(2:l)
      l=l-1
      goto 10

20    if(segno)then
        str(2:l+1)=str(1:l)
        str(1:1)='-'
        str=str(1:l+1)
      else
        str=str(1:l)
      end if

      return
      end

      
c==========================================
      subroutine flostr(flo,str,ndec)
c     Translate a real*8 value into string

      implicit none

      integer :: l, ndec, length, n1, n
      character *(*)str
      character *20 temp
      real *8 flo
    
c     translate the integer portion
      n=int(flo)
      n1=int(abs(flo))
      call intstr(n,str,0)

c     transform the decimal portion in integer 
      l= length(str)
c
c     Modified by BCox, change FLOAT to REAL...might not work right?
c     Problem was that abs returns a real and float takes an integer
      n=int((REAL(abs(flo)-n1))*10**(ndec))
c     translate the decimal portion
      call intstr(n,temp,ndec)

c     join two portion
      str(l+1:l+1)='.'
      str(l+2:l+ndec+1)=temp
      str=str(1:l+ndec+1)

      return
      end
      
      
      
      