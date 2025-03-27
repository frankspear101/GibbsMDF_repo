! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE STEPS
    	use AWE_Interfaces
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "SaveRxn.inc"
	include "TPSave.inc"
! *****************************************
!      local variables
      integer*4 ifdplt,i,isteps,j,L,K,kcur,istep,izero,irefer,jsymb,NeedToPick,iok,Picked(10),myColor,nStop
      real*4 ValuesPicked(10)


	TPCounter = 1
	TPSave(TPCounter,1) = ALLX(1,1)
	TPSave(TPCounter,2) = ALLX(1,2)

      Plotcounter=0
      ifdplt=0
!	Set up initial monitors and deltas
!	Find where the moles of the first MDF phase starts
	j = 2
	do 20001 kcur = 1,numPh
	k = AsmCurrent(kcur)
	j = j + numPhCo(k) - 1
20001	continue
!	The first phase is the grain boundary -- not a monitor parameter
	j = j + 1		! should be M_Grain_boundary
	SMon(3) = j + 1
	SMon(4) = j + 2
!	The deltas must be scaled to conserve oxygen
!	Sdel(3) = -.6/1000.		! this assumes quartz O = 2
!	Sdel(4) = .1/1000.		! this assumes garnet O = 12



1     CONTINUE
      WRITE(*,*)' *********************************'
      write(*,4)TC,PB
4   format (' T = ',f10.1,'  P = ',f10.1)
	call Variance			! write out variance
      WRITE(*,*)' CURRENT MONITOR PARAMETERS ARE'
      if(nvar-neq.gt.0)then
      do 6 i = 1,nvar-neq
      j = smon(i)
      if(VN1(j)(1:2).eq."M_")then
		WRITE(*,5)j,VN1(j),VN2(j),ALLX(1,j)*1000.,SDEL(I)*SNSTEP*1000.
		else
		WRITE(*,5)j,VN1(j),VN2(j),ALLX(1,j),SDEL(I)*SNSTEP
		endif
5     FORMAT(I3,1X,A2,A6,2F13.4)
6     continue
      else
      write(*,*)' No monitor parameters - variance <= 0'
      endif
      WRITE(*,*)'NSTEP =',SNSTEP
      if(idraw.eq.0)write(*,*)'IDRAW = 0 (do not draw the line)'
      if(idraw.eq.1)write(*,*)'IDRAW = 1 (draw the line)'


! ----------------------------------------------------
	if(iternary.eq.0)then			! X-Y plot
		write(*,*)'Current X-Y axes are:'
		write(*,*)'       X             Y'
		do 9 k = 1,NoTieMinerals
		i = PlotIndex(k,1)
		j = PlotIndex(k,2)
9		write(*,7)i,VN1(i),VN2(i),j,VN1(j),VN2(j)
! 9		write(*,7)PlotIndex(k,1),VN1(PlotIndex(k,1)),VN2(PlotIndex(k,1)),PlotIndex(k,2),VN1(PlotIndex(k,2)),VN2(PlotIndex(k,2))
7		FORMAT(4(I4,' ',A2,A8))

	else   					! ternary plot
	      	write(*,*)'Coordinates are:'
	      	write(*,12)(CompPlotName(Apx(i)),i=1,iternary+2)
12		Format(4(A8,8x))
		do 13 k = 1,NoTieMinerals
13		write(*,14)PlotIndex(k,1),PhName(PlotIndex(k,1))
14		FORMAT(4(I4,' ',A12))

	endif
!-----------------------------------------------



      WRITE(*,*)' *********************************'
     	If(inewton.eq.1.and.imass.eq.1)write(*,*)'Newton''s method; Mass balance = ON'
     	If(inewton.eq.1.and.imass.eq.0)write(*,*)'Newton''s method; Mass balance = OFF'
     	If(inewton.eq.0.and.imass.eq.1)write(*,*)'Gibbs'' method; Mass balance = ON'
     	If(inewton.eq.0.and.imass.eq.0)write(*,*)'Gibbs'' method; Mass balance = OFF'
	if(imass.eq.1)then
		if(bulkCompSwitch)then
			write(*,*)'Bulk composition = ',bulkCompTitle
		else
			write(*,*)'Bulk composition = (none - bulk comp may change)'
		endif
	endif
      WRITE(*,*)' *********************************'
      WRITE(*,*)' SINGLE STEPS MENU OPTIONS'
      WRITE(*,*)'    0 = Return'
      !WRITE(*,*)!'    1 = Line or tie-lines (don''t touch unless you want to change plotting)'
      WRITE(*,*)'    2 = Choose monitors/set deltas'
      WRITE(*,*)'    3 = Compute one increment - EQUIL routine'
      WRITE(*,*)'    33 = Compute one increment - EQUIL vacancy routine'
      WRITE(*,*)'    4 = Compute one increment - MDF2 routine'
      WRITE(*,*)'   44 = Compute one increment - MDF3 routine'
      WRITE(*,*)'    5 = Unstep one finite difference'
      WRITE(*,*)'    6 = Reference points'
      WRITE(*,*)'    7 = Print assemblage information to Output window'
      WRITE(*,*)'    8 = Go to global menu'
      write(*,*)'    9 = Go to Plotter options'
      write(*,*)'   10 = Gibbs/Newton switch'
      write(*,*)'   20 = Begin/Save menu'
      write(*,*)'---------------------------------'
      write(*,*)'   -1 = Change mineral assemblage'
      write(*,*)'   -2 = Set output length'
      write(*,*)'   -3 = NewtonStepsOutput switch'
      write(*,*)'   -4 = Pick line color'
      WRITE(*,*)' CHOOSE OPTION'
      read(*,*)isteps
      
      
	Select Case (isteps)

! ---------------------------------------
	case (-1)
!      		if we are changing the assemblage, we must zero out any "new" variables
!      		because there is no way to determine whether they will still be in the 
!      		list after a mineral is removed.
      		NVAR=NVAR-numNew
      		NEQ=NEQ-numNew
      		numNew=0
       		call change
!       	call rexn
		go to 1
! ---------------------------------------
	case (-2)
		isteps = 0
		call AWE_setoutput
		goto 1
! ---------------------------------------
	case (-3)
		write(*,*)'NewtonStepsOutput = ',newtonStepsOutput
		write(*,*)'Input new value (0 or 1)'
		read(*,*)newtonStepsOutput
		if(newtonStepsOutput.eq.0)close(43)
		if(newtonStepsOutput.eq.1)then
			open(43,file='',status='UNKNOWN')
			endif
		goto 1
! ---------------------------------------
	case (-4)
!		myColor = 1		! Keep the same myColor as previously
		call PickAWEColor(myColor)
		CurrentColor = myColor				! CurrentColor is in PlotStuff.inc common block
		goto 1
	
! ------Begin/Save routine-------------
	case(20)
		call GibbsBegin
		go to 1
! ----------------EXIT-----------------
	case(0)
		RETURN

!  ----Reset reaction counter------------------------------------------
	case(11)
!  ----inewton------------------------------------------
	case(10)
		call GibbsorNewton
!  ----Set plot stuff------------------------------------------
	case(1)

!  ----CHOOSE MONITORS/set deltas-------------------------------
	case(2)
		write(*,*)'nvar, neq ',Nvar,neq
		NeedToPick = Nvar-neq
		do i = 1,NeedToPick
			Picked(i) = SMon(i)
			ValuesPicked(i) = SDel(i)*Dfloat(SNstep)
			end do
		call PickMonitors(nvar,NeedToPick,Picked,ValuesPicked,iok)
        	if(iok.eq.0)then
			if(SNStep.le.0)SNStep = 1
			do i = 1,NeedToPick
				SMon(i) = Picked(i)
				SDel(i) = ValuesPicked(i)/DFloat(SNstep)
				end do
			endif

!  -----SET DELTAS----------------------------------
!	case(3)
		! This routine now done in subroutine Pkmon
! ----COMPUTE ONE INCREMENT------------------------
	case(3,33,4,44)
		jsymb=isymb
		if(isymb.eq.0)jsymb=3
!      		set monitor parameters
		NSTEP=SNSTEP
		if(isteps.eq.3.or.isteps.eq.33)then
			nstop = numEqEquns
			else
			nstop = numMDFEquns		! this is the number of MDF equations (numMDFrxn + nc)
			endif
!		do i=1,nvar-neq
		do i=1,nvar-nStop
			deltax(i)=sdel(i)
			mon(i)=smon(i)
			end do
!      		Set up array IPOINT to contain pointers to non-monitor parameters
		J=0
		Do i=1,NVAR
			if(isteps.eq.3.or.isteps.eq.33)then
				nstop = numEqEquns	! This is the number of equations total (equilibrium + mass balance)
				else
				nstop = numMDFEquns		! this is the number of MDF equations (numMDFrxn + nc)
				endif
!			DO L=1,NVAR-NEQ
			DO L=1,NVAR-nstop
				IF(MON(L).eq.i)go to 410
				end do
			J=J+1
			IPOINT(J) = I
410	   		CONTINUE
			end do
		IF(iLong(2).eq.1) write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,numEqEquns)

!      		THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
      		DO I=1,nvar
!	         	ALLX(3,I) IS WHERE calculations STARTED-save
			ALLX(3,I)=ALLX(1,I)
			end do

	      	do istep=1,snstep
			if(inewton.eq.1)then
				!Store new values for each independent variable (monitor parameter)
				!This is where Neton's method increments composition, T or P
				DO J=1,NVAR-NEQ
					i=mon(j)
					ALLX(1,i)=ALLX(1,i) + DELTAX(J)
					end do
				 call SetTPX
				 endif
 
			write(*,461)TC,PB
461   			format (' T = ',f10.1,'  P = ',f10.1)
			write(*,*)' NSTEP = ',snstep,' Counter = ',istep
			IZERO=0
			select case(isteps)
				case(3)
				call Compute4EQUIL(izero)
				case(33)
				call Compute4EQUIL_Vac(izero)
				case(4)
				call Compute4MDF2(izero)
				case(44)
				call Compute4MDF3(izero)
				case default
				end select

			IF(IZERO.EQ.1)GO TO 470
			!Check for mass balance error
			if(imass.eq.1)then
				call CheckMassBal
				endif
			end do

470		continue

		call CalcAffinity(1)
!  ---UNSTEP---------------------------------------
	case(5)
	CALL RESET(5,1)

!  ---REFERENCE POINTS--------------------------
	case(6)
		WRITE(*,*)' REFERENCE POINT OPTIONS'
		WRITE(*,*)'   0 = Return  '
		WRITE(*,*)'   1 = Set reference point'
		WRITE(*,*)'   2 = Return to reference'
		WRITE(*,*)'   3 = Return to starting conditions'
		READ(*,*)IREFER
!		Definition of ALLX variable:
!      		1,50   current value
!      		2,50     starting value
!      		3,50   ref at start of contour
!      		4,50   user selected reference
!      		5,50   previous finite diff point
!      		6,50     (not used)
		IF(IREFER.EQ.1) THEN
			DO I=1,nvar
				ALLX(4,I)=ALLX(1,I)
				end do
			ENDIF
	       IF(IREFER.EQ.2) CALL RESET(4,1)
	       IF(IREFER.EQ.3) CALL RESET(2,1)
!  ----PRINT------------------------------------------
	case(7)
		CALL PRINTT(1)
!  ----GLOBAL------------------------------------------
	case(8)
		CALL GLOBAL(isteps)
! ----------plotter stuff------------------------------
	case(9)
		CALL PLOTIN

! ----------------------------------------
	case default
		go to 1
	end select
! ----------------------------------------
        go to 1

	END


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




	
