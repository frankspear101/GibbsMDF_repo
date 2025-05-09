!*******/MONIT/*****************************************************
!     THESE ARE DIMENSIONED FOR COMMON BLOCK /MONIT/
!     DIMENSION (16) REFERS TO THE MAXIMUM NUMBER OF MONITOR PARAMETERS
!     DIMENSION (MaxPCMonit) REFERS TO THE MAXIMUM NUMBER OF VARIABLES IN A PROBLEM
!	ncont is the contour routine monitors (initialized at startup)

	integer maxPCMonit; parameter(maxPCMonit=100)		! maximum total phase components in a problem 
	! maxPCMonit must be the same value as maxPC in Assemb.inc -- if you change it here, you must change it there 
      INTEGER*4 nStep,snStep,mon(50),smon(50),iPoint(maxPCMonit),nCont(16),allFractl(maxPCMonit)
      integer*4 molesPhPtr,openPtr,newVarPtr
      REAL*8 deltaX(50),sDel(50),allX(6,maxPCMonit)
      CHARACTER vn1(maxPCMonit)*4,vn2(maxPCMonit)*16
      COMMON /MONIT/ allX,deltaX,sDel,vn1,vn2,nStep,			&
      		iPoint,snStep,mon,smon,nCont,molesPhPtr,openPtr,newVarPtr,allFractl
!     Definition of ALLX variable:
!     1,50   current value
!     2,50     starting value
!     3,50   ref at start of contour
!     4,50   user selected reference
!     5,50   previous finite diff point
!     6,50     (not used)
!********************************************************************
! 	Note: 
!	iPoint contains pointers to the neq non-monitor parameters.
!	iPoint is set in routine that calls Compute (e.g. Steps, Contour)

!	These next 3 are set in Subroutine Names
!	molesPhPtr is a pointer to the first phase moles (Mk) in allX array
!	openPtr is a pointer to the first open system component (mi) in allX array
!	newVarPtr is a pointer to the first new variable in allX array
