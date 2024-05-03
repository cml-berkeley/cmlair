!***************************************************************
module Q4_sizes
!***************************************************************
    implicit none
    
    integer, parameter :: nxx=1000, nyx=1000, nmx=1000 !GRID SIZES

    integer, parameter :: maxrls=200, maxrpt=300 !rail limits
    integer, parameter :: maxwls=200, maxwpt=50 !wall limits
    integer, parameter :: maxseg=100  !maximum number of grid segments
    integer, parameter :: maxrad=100, maxrpm=100, maxalt=100 !radius, rpm, altitude

end module Q4_sizes


!***************************************************************
module ShearArrays
!***************************************************************
    implicit none
    
    real(8), dimension(:,:), allocatable :: PoisilleShearX, CouetteShearX
    real(8), dimension(:,:), allocatable :: PoisilleShearY, CouetteShearY
end module ShearArrays


!***************************************************************
module Q4_globals
!***************************************************************
    use Q4_sizes
    implicit none
    
    integer, parameter :: ON = 1, OFF = 0
	
    character(7), parameter :: version = '4.032'	
    real(8), parameter :: twopi = 6.283185307179586


    integer :: VCycleIteration

    !Grid
    !==========
    integer :: nx,ny,nxm1,nym1,nxm2,nym2, nsx, nsy
    integer :: nx1, ny1, nx2, ny2, nx3, ny3, nx4, ny4
    
    real(8) :: dxr(maxseg),dyr(maxseg), xnt(maxseg),ynt(maxseg)
    real(8), dimension(:), allocatable :: xcontrol, ycontrol
    integer :: maxControlPts
    
    
    real(8), dimension(:), allocatable :: xref, xref1, xref2, xref3, xref4, &
                                            yref, yref1, yref2, yref3, yref4
    
!   disk spacing..
    real(8), dimension(:,:), allocatable :: sliderdiskspacing
        
    integer :: nxt(maxseg), nyt(maxseg)

    integer :: isymmetry, ioldgrid, iadpt, nxm(4), nym(4), icntlx, icntly
    integer :: iUseNewGridEachRun

!	MinFH is the actual mininum point
!	him is the nominal TEC height
!	hm is reference height (often = initial him)

    !Run Setup and Initial attitude
    real(8) :: ske, ra, rpm, u0, hm, hx0, h0, hs, hy, p0xl, f2p      		    
    
    !slider info
    !xl = x length, yl = y length, zl=**, xg = **,  xt=taper Length, ht = taper angle, rebase = base recess
    real(8) :: xl, yl, zl, xg, xt, ht, rebase, crown, camber, twist


    !air parameters.
    real(8) :: al, vis1, p0
    logical :: useStandardAir
    
    integer :: iOutputShear
    logical :: staggerShear

    !constants
    real(8) :: d0, gama, pir, pit, t1, t2, t3, t4, slip_beta, slip_gamma, accom, corCoef

    !used by the FK database
    integer, parameter :: FK_TABLE_SIZE = 1001
    real(8) :: qndat(FK_TABLE_SIZE), zdat(FK_TABLE_SIZE)
    integer :: icoe, nter, iqpo

    !velocity of disk under the slider in the X and Y direction
    real(8), dimension(:,:), allocatable :: bearx, bearx1, bearx2,      &
    	       bearx3, bearx4, beary, beary1, beary2, beary3, beary4

    !height data (say what hnew and h are)
    real(8), dimension(:,:), allocatable :: h, h1, h2, h3, h4, hnew, hnew1, hnew2, hnew3, hnew4
    

    !pressure
    !========
    real(8), dimension(:,:), allocatable :: p, p1, p2, p3, p4
    logical :: useInitialPressure
    

    !reynolds eqn residuals
    !=====================
    real(8), dimension(:,:), allocatable :: res, res1, res2, res3, res4, &
                                            su01, su02, su03, su04

    !rails
    !=====
    real(8), dimension(:,:), allocatable :: xrail, yrail, xdiff, ydiff, cramp, hramp
    real(8), dimension(:), allocatable ::xpt1, xpt2, ypt1, ypt2
    integer, dimension(:), allocatable :: npoints, istep
    integer :: numRails, maxRailPts, numWalls, maxWallPts
    !walls
    !=====
    real(8), dimension(:,:), allocatable :: wpoint,wrecess, &
                    xwalli, ywalli, xwallo, ywallo,         &
                    flush_lower_xwallo, flush_lower_ywallo, &
                    flush_upper_xwallo, flush_upper_ywallo, &
                    xw1, yw1, xw2, yw2
    integer, dimension(:,:), allocatable :: indexw
    integer, dimension(:), allocatable :: nwpoint

    !reynolds equation
    !=================
    real(8) :: akmax, ak, ak0, ak1, ak2, ak3, ak4
    integer :: itrey, idisc, mg_nest, max_mg_nest, GS_Iterations

    !multigrid nesting
    integer :: mitr(4),mitp(4)


    logical crash

    
    !FORCES
    !================
    real(8) :: emax, err, f, fneg, fpos, fsp, fsr, xf, yf, &
               !xint(4), yint(4), hint(4), &
               hmin, MinFH, MinFHLocX, MinFHLocY, &
               jac(3,3), rint(4), hintge(4),hgap(4), &
               Xmom, Ymom, Zmom, xPosLoc, yPosLoc, xNegLoc, yNegLoc, &
               fvdw_output
    real(8), dimension (:,:), allocatable :: vdwMolecularForceMap
    real(8), dimension (:), allocatable :: xintNew, yintNew, hintNew       
    integer :: numRes, numPOI
    !forces on the slider
    real(8) :: f0, xf0, yf0, xfs, yfs, Pitch_Stiffness, Roll_Stiffness, PSA, RSA
    integer :: iUseStiffnesses


    !contact and surface roughness
    !=============================
    real(8) :: cp(nxx,nyx),rsik,cta,rasp,fcr,txr,tyr, &
    		   aratio,eyoung,ydst,ydcoe,pratio,frcoe,ahc,bhc,elecpot
    integer :: icmod, iUseIMF

    !runtime flags
    !=============
    integer :: hflag !hflag: indicate change of slider geometry or grid
    integer :: uflag !uflag: indicate change of velocity or grid
    integer :: pflag !pflag: indicate save pressure for the case
    integer :: sflag !sflag: indicate calculate stiffness
    integer :: tflag !tflag: indicate first traversal of disk radius


    !various parameters needed by the reynolds equation
    !many of these are set in ave_height()
    real(8), dimension(:,:), allocatable :: cohimx,cohimx1,cohimx2, &
     		    cohimx3,cohimx4,cohjmx,cohjmx1,cohjmx2,cohjmx3,cohjmx4, &
     		    himax,himax1,himax2,himax3,himax4, &
     		    himin,himin1,himin2,himin3,himin4, &
     		    hjmax,hjmax1,hjmax2,hjmax3,hjmax4, &
     		    hjmin,hjmin1,hjmin2,hjmin3,hjmin4, &
     		    recssi,recssi1,recssi2,recssi3,recssi4, &
     		    recssj,recssj1,recssj2,recssj3,recssj4

    !run conditions
    !==============
    real(8) :: radii(maxrad), skews(maxrad), rpms(maxrpm), alts(maxalt)
    integer :: irad, irpm, ialt, iRadiiSort
	integer :: icurrad, icuralt, icurrpm, icursen
	integer :: iRunCounter
	
	!run type
	!========
    integer isolv, isave, istiff, irailgeom
    logical gapExtrapolationRun
    
    !sensitivities
    !=============
    real(8) :: crninc, cbrinc, twstinc, tlnginc, tanginc, sldinc, ptqinc, rtqinc, rcsinc
    integer :: iwscale

    !date and time info
    !==================
	integer :: start_time(8), end_time(8)

    !F-K database arrays
    !=======================
	real(8) :: FKUpper(19), FKMidUpper(19), FKMidLower(19), FKLower(19)
	
end module Q4_globals

!**************************************************************
module AdaptiveGridParams
!**************************************************************
    implicit none
    
    ! Adaptive Grid 
    real(8) :: difmax_x, decay_x
    real(8) :: difmax_y, decay_y
    integer :: iAdaptMode_x, ipmax_y, ipmax_x  !iadaptMode_x: 0 = Pressure Adaptation; 1 = Geometric Adaptation    
    ! mesh refinement
    real(8) :: rlevelx(12), rlevely(12), xginit(12), yginit(12), xgend(12), ygend(12)
	integer :: nxmr, nymr

end module AdaptiveGridParams


!**************************************************************
module PressureWarnings
!***************************************************************
    implicit none
    
    integer numPressureWarnings
    
end module PressureWarnings


!**************************************************************
module UserDefinedGeometry
!**************************************************************
    implicit none
    
    !I really didn't want to use user defined types (b.t.w. who the fuck thought up the % as the character to dereference?)
    !but it just made a lot more sense use one here
    integer :: numUDGRegions
    type UserDefGeomType
        real(8) :: dLowerX, dUpperX, dLowerY, dUpperY
        integer :: nSizeX, nSizeY

  !you will need to enable the preprocessor in the project settings
  !to compile under Compaq visual fortran uncomment the next line
!#define COMPAQ_VISUAL_FORTRAN__

        real(8), dimension(:,:), allocatable :: UDGMatrix

    end type UserDefGeomType
    
    type (UserDefGeomType), dimension(:), allocatable :: UserDefGeom

end module UserDefinedGeometry


!**************************************************************
module NestingInfo
!***************************************************************
    implicit none
    !holds the coords of the lower x and y indexes of the rectangle on the coarser mesh that
    !encloses the finer mesh point
    !numbers after the variable name (1-4) have kept the same convention as the grid names
    !ie Y1 is the second finest grid (#2) and Y is the finest grid (#1)
    real(8), dimension(:), allocatable :: enclosingCoarseRectX, enclosingCoarseRectY
    real(8), dimension(:), allocatable :: enclosingCoarseRectX1, enclosingCoarseRectY1
    real(8), dimension(:), allocatable :: enclosingCoarseRectX2, enclosingCoarseRectY2
    real(8), dimension(:), allocatable :: enclosingCoarseRectX3, enclosingCoarseRectY3
    !holds the coords of the lower x and y indexes of the rectangle on the finer mesh that
    !encloses the coarser mesh point
    real(8), dimension(:), allocatable :: enclosingFineRectX1, enclosingFineRectY1
    real(8), dimension(:), allocatable :: enclosingFineRectX2, enclosingFineRectY2
    real(8), dimension(:), allocatable :: enclosingFineRectX3, enclosingFineRectY3
    real(8), dimension(:), allocatable :: enclosingFineRectX4, enclosingFineRectY4

end module NestingInfo


!**************************************************************
module OptionalParams
!***************************************************************
    implicit none
    
    integer, parameter :: LINESIZE = 1024
    integer, parameter :: PARAMSIZE = 255
    integer, parameter :: MAXPARAMS = 10

end module OptionalParams


!**************************************************************
module FKLookupTable
!**************************************************************
    implicit none

    real(8), dimension (:), allocatable :: zdatLow, qndatLow
    real(8), dimension (:), allocatable :: zdatMid, qndatMid
    real(8), dimension (:), allocatable :: zdatHigh, qndatHigh

    real(8) :: zBaseLow, zBaseMid, zBaseHigh, zBaseCeil
    real(8) :: zStepLow, zStepMid, zStepHigh
    
end module FKLookupTable

!
!!**************************************************************
!module ResidualSmoothingModule
!!**************************************************************
!    implicit none
!    
!    real(8), dimension(:,:), allocatable :: resSmoothRes, resSmoothRes1, resSmoothRes2,         &
!                                            resSmoothRes3, resSmoothRes4
!    real(8), dimension(:,:), allocatable :: resSmoothPress, resSmoothPress1, resSmoothPress2,   &
!                                            resSmoothPress3, resSmoothPress4
!    
!end module ResidualSmoothingModule



!***************************************************************
module TemperatureHumidity
!***************************************************************
    implicit none

    real(8) :: temperature, humidity
    real(8), allocatable, dimension(:,:) :: savedPressure
    logical :: doHumidity
    
    
end module TemperatureHumidity


!***************************************************************
module KrylovAccel
!***************************************************************
    implicit none

    real(8), dimension(:,:,:), allocatable :: r_dot, u_dot
    integer, parameter :: savedKrylovIterations = 7
    integer :: krylovIteration, krylovResidualNorm

end module KrylovAccel