! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE AdjustAsm(ichanged)
      implicit none
!      ROUTINE TO adjust the arrays after changing the mineral assemblage
! *****************************************
	include "Assemb.inc"
! *****************************************
!      LOCAL VARIABLES
      	INTEGER*4 k,kCur,iChanged,i,j

!		Calculate total number of phase components
		NP=0
	!	newVarGlobal = 0			! flag that says we have new variables to consider (0 = no, 1 = yes)
      		do 15 kCur = 1,numPh
		k = asmCurrent(kCur)
		NP = numPhCo(K)+NP
		if(iChanged.ne.0)then
		!	if(k.eq.1)then		! set grain boundary volume = 1. initially
		! The model assumes the same reactive volume for each phase. 
		!This is arbitrarily set to 1.0 and includes solid phases and the grain boundary
		! In other words, if the grain boundary is around 1 nm wide, then so is the reactive volume of the soid phases
		! The molar volumes of the GB components is 0.853 and they are all the same
		! The molar volumes of the solid phases are those in the database
		!
		! The fractl code is used in the MDF routines to set the composition to the starting composition for "fractl = 2"
		! It shouldn't change the reactive volume, so these lines of code are commented out
		Vp0(k) = 1.
! 			if(fractl(k).eq.0)then		! phases that fractionate = 1 (i.e. Garnet ± chlorite)
! 				Vp0(k) = 1.
! !				Vp0(k) = 1.0d-9		! cm^3 units -- this is 1nm x 1mm x 1mm volume in cm^3
! 				Vp1(k) = VP0(k)
! 				Vp2(k) = VP0(k)
! 				else		! set the amounts of all phases = 0 initially
! 				Vp0(k) = 0.
! 				Vp1(k) = VP0(k)
! 				Vp2(k) = VP0(k)
! 				endif

			! note that vmol(K) hasn't been reset at this point. 
			!    so the value is from the last time this phase was in a call to GHSV
			! This is probably OK because vmol(K) is not strongly dependent on composition
			! But, to be precise, I should call GHSV before making this calculation
			!mp0(K)=vp0(K)/(vmol(K)*10.0D0)
			mp0(K) = vp0(K)/(vmolStart(K)*10.0D0)	! use the initial value of vmol(k) that is set on file input.
								! this avoids cases where vmol(k) = 0 from the previous set of calculations
			! mp0(k) is the number of moles of a phase that is reacting along the grain boundary
			! It will change depending on whether the phase is consumed or produced
			mp1(k)=mp0(k)
			mp2(k)=mp0(k)
			MPSTART(k)=mp0(k)			! mpStart is the starting moles of each phase
		!	mp0(k) = 0.001		! reset the moles of phases to .1  (100 millimoles) but only if we actually changed the assemblage

			endif
15		CONTINUE
! 	   	Compute the starting number of moles of each system component
		do 310 i=1,nc
			molesStart(i)=0.D0			! molesStart is the starting moles of each system component
			Do 312 kCur = 1,numPh
				k = asmCurrent(kCur)
				do 315 j=1,numPhCo(K)
					molesStart(i) = molesStart(i) + MPStart(K)*xPhCo(k,j)*comp(k,j,i)
					! Note that xPhCo(k,j) is set equal to the nodePhaseComp(iNode,i,j) in the calling routines
			315	   	continue
		312	   	continue
310	   	continue
!     		compute number of independent compositional (dX) variables (before SetAllX)
	!	don't reset the compositions .... it's a pain
	!	xPhCo = xPhCoInitial 	! set all compositions back to initial (input) values
      		NX=NP-numPh		! This is the number of independent compositional variables (first X in each phase is dependent)
		call SetALLX()		! set up ALLX array
		call Names()		! sets up names of variables for asmCurrent
      		call REXN		! calculate linearly independent reactions
!     		compute total number of variables
      		NVAR= TANDP + NX + numPh
		if(EqOrMDFProblem.eq.0)then
			NEQ = numEqRxn + nc
			else
			NEQ = numMDFRxn + nc
			!NEQ = NEQ + NC			! add a mass balance equation for each system component
			! Note that neq is set neq = numMDFrxn in routine REXN. This will give the wrong answer for equilibrium calculations
			endif
      		return
		end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE EQNode(iNodeTemp,numPhT,ph1,ph2,ph3,output,initializeOnly)
	use MatrixArrays
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
!	include "Solute.inc"
	include "Diffuse_GB_Hex.inc"
	include "AutoAffinity.inc"
! *****************************************
	integer*4 iNode,i,j,k,k1,k2,k3,k4,l,izero,iNodeTemp,kk
	!real*8 deltaMoles,deltaMolesdX,dF_dDelta,deltaAffinityInit,deltaAffinityOld,deltaAffinityNew
	!real*8 deltaMolesdX2,deltaMolesdX3,deltaMolesdX4,dF3_dDelta2,dF3_dDelta3,dF4_dDelta2,dF4_dDelta3,		&
	!		deltaAff3_old,deltaAff3_new,deltaAff4_old,deltaAff4_new,F3o,F4o
	!real*8 deta
	integer*4 output,initializeOnly,numPhT,ph1,ph2,ph3
!	if output = 1 then print
!	if output = 0 then don't print


!	do i = 1,numNodesWith3Phases
!		if(nodesWith3PhasesIndex(i).eq.iNodeTemp)then
!			iNode = i
			! iNode is the index for this node in the array nodesWith3PhasesIndex
!			go to 3
!			endif
!		end do
!	call FSS_Alert('ALERT','Did not find the index for this node')
!	return

	iNode = iNodeTemp
3	continue		

	snstep = 1
	numPh = numPhT + 1		! should be 3 or 4
	asmCurrent(1) = 1			! always the grain boundary
	asmCurrent(2) = nodeMIFID(iNode,ph1)
	asmCurrent(3) = nodeMIFID(iNode,ph2)

	k1 = asmCurrent(1)
	k2 = asmCurrent(2)
	k3 = asmCurrent(3)
	if(numPhT.eq.3)then
		asmCurrent(4) = nodeMIFID(iNode,ph3)
		k4 = asmCurrent(4)
		endif

	if(output.eq.1)then
		write(12,*)'  '
		write(12,*)'  '
		write(12,*)'  '
		write(12,*)'********************************************************'
		write(12,*)' Current mineral assemblage'
		do k = 1,numPh
			WRITE(12,1533)MINREC(asmCurrent(K)),PHNAME(asmCurrent(K))
1533  			FORMAT(3(I4,2X,A8,2x))
			end do
		endif	
	do j = 1,numPhCo(1)		! grain boundary
		xPhCo(1,j) = nodeComp(iNode,j)		! numPhCo(1) should be every element
		end do
	do i = 1,numPhT		! do the 2 or 3 solid phases
		k = asmCurrent(i+1)
		kk = nodeReactionPhases(iNode,i)	!Note that the index of phases around the node isn't 1,2,3 necessarily
		select case(fractl(k))
		case(0,1)	! phases are not fractionating or are growing and fractionating
			do j = 1,numPhCo(k)
				xPhCo(k,j) = nodePhaseComp(iNode,kk,j)
				end do
		case(2)		! phase is being consumed - use fixed composition (i.e. chlorite)
			do j = 1,numPhCo(k)
				xPhCo(k,j) = xPhCoInitial(k,j)
				end do
		case default
		end select		
		end do


	CALL AdjustAsm(1)		! we changed assemblage. Adjust volumes, moles etc.

!	The variance should be 2 (Duhem's theorem)
	sdel(1) = 0.0d0			! T
	sdel(2) = 0.0d0			! P
!	sdel(3) = 0.0d0			! first time through, just do calculation to see which phase has larger affinity
!	sdel(4) = 0.0d0
!	sdel(5) = 0.0d0
	smon(1) = 1
	smon(2) = 2
	!        P&T   (   number of dX terms          - 4 dependent dXterms   ) + previous dMterms
!	smon(3) = 2 + numPhCo(k1) + numPhCo(k2) + numPhCo(k3) + numPhCo(k4) - 4 + 2
!	smon(4) = 2 + numPhCo(k1) + numPhCo(k2) + numPhCo(k3) + numPhCo(k4) - 4 + 3
!	smon(5) = 2 + numPhCo(k1) + numPhCo(k2) + numPhCo(k3) + numPhCo(k4) - 4 + 4
!      		set monitor parameters
	NSTEP=SNSTEP
	do i=1,nvar-numEqEquns		
		deltax(i)=sdel(i)	! these are all = 0 at this point
		mon(i)=smon(i)
		end do
!      		Set up array IPOINT to contain pointers to non-monitor parameters
	J=0
	Do i=1,NVAR
		DO L=1,NVAR-numEqEquns	
			IF(MON(L).eq.i)go to 1510
			end do
		J=J+1
		IPOINT(J) = I
1510		CONTINUE
		end do
	if(output.eq.1)then
		call PrinTT(1)
		write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NEQ)
		endif
!      	THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
	DO I=1,nvar
!         	ALLX(3,I) IS WHERE calculations STARTED-save
		ALLX(3,I)=ALLX(1,I)
		end do
	call SetTPX

!	This first call is required because it sets the phase composition to be that of the MDF criteria
	IZERO=0
	if(output.eq.1) then
		write(12,*)' Call to CalcAffinity(1) to get affinities'
		call CalcAffinity(1)	! print this initial calculation
		write(12,*)' '
		write(12,*)' Calling Compute4EQUIL in subroutine EQNode (2201)'
		endif
	! store the initial moles for the phases
	MDFmoles(1,2) = mp0(k2)
	MDFmoles(1,3) = mp0(k3)
	if(numPhT.eq.3)then
		MDFmoles(1,4) = mp0(k4)
		endif
!	pause 'ready to call Compute4EQUIL -- Check the output window -- hit return when ready'
	call Compute4EQUIL(izero)
	if(izero.gt.0)then
		write(*,*)' Something went wrong in Compute4EQUIL '
		endif
		
	! we should now have the equilibrium assemblage

!	Store the final moles
	MDFmoles(2,2) = mp0(k2)
	MDFmoles(2,3) = mp0(k3)
	if(numPhT.eq.3)then
		MDFmoles(2,4) = mp0(k4)
		endif

	!store the new grain boundary composition
	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		GBComp(2,j) = xPhCo(1,j)
		end do
	if(output.eq.1)then
		call Printt(1)
		endif

!	At this point, we need to 
		!(a) save the GB composition and 
		!(b) store the phase compositions to examine zoning
		!(c) store the amount of each phase produced or consumed 
	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		nodeComp(iNode,j) = xPhCo(1,j)
		end do
	do j = 1,numPhCo(k2)
		nodePhaseComp(iNode,1,j) = xPhCo(k2,j)			! nodePhaseComp(maxNodes,3,maxPhCo)
		end do
	do j = 1,numPhCo(k3)
		nodePhaseComp(iNode,2,j) = xPhCo(k3,j)
		end do
	if(numPhT.eq.3)then
		do j = 1,numPhCo(k4)
			nodePhaseComp(iNode,3,j) = xPhCo(k4,j)
			end do
		endif
	! these are the total moles since the model began
	nodePhaseMoles(iNode,1) = nodePhaseMoles(iNode,1) + MDFmoles(2,2)-MDFmoles(1,2)
	nodePhaseMoles(iNode,2) = nodePhaseMoles(iNode,2) + MDFmoles(2,3)-MDFmoles(1,3)


	! these are the change in moles for this iteration
	nodePhaseMolesDelta(iNode,1) = MDFmoles(2,2)-MDFmoles(1,2)
	nodePhaseMolesDelta(iNode,2) = MDFmoles(2,3)-MDFmoles(1,3)
	if(numPhT.eq.3)then
		nodePhaseMoles(iNode,3) = nodePhaseMoles(iNode,3) + MDFmoles(2,4)-MDFmoles(1,4)
		nodePhaseMolesDelta(iNode,3) = MDFmoles(2,4)-MDFmoles(1,4)
		endif
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE EQSegment(iSegTemp,iPointX,output,initializeOnly)
	use MatrixArrays
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
!	include "Solute.inc"
	include "Diffuse_GB_Hex.inc"
	include "AutoAffinity.inc"
! *****************************************
	integer*4 iSeg,i,j,k,k1,k2,k3,l,izero,iSegTemp,iPointX
	integer*4 output,initializeOnly
!	if output = 1 then print
!	if output = 0 then don't print


	iSeg = iSegTemp
3	continue		

	snstep = 1
	numPh = 3
	asmCurrent(1) = 1			! always the grain boundary
	asmCurrent(2) = SegMIFID(iSeg,1)
	asmCurrent(3) = SegMIFID(iSeg,2)

	k1 = asmCurrent(1)
	k2 = asmCurrent(2)
	k3 = asmCurrent(3)

	if(output.eq.1)then
		write(12,*)'  '
		write(12,*)'  '
		write(12,*)'********************************************************'
		write(12,*)' In Subroutine EQSegment '
		write(12,*)'********************************************************'
		write(12,*)' Current mineral assemblage'
		do k = 1,numPh
			WRITE(12,1533)MINREC(asmCurrent(K)),PHNAME(asmCurrent(K))
1533  			FORMAT(3(I4,2X,A8,2x))
			end do
		endif	
	do j = 1,numPhCo(1)		! grain boundary
		xPhCo(1,j) = pointComp(iSeg,iPointX,j)		! numPhCo(1) should be every element
		end do
	do i = 1,2		! do the 2 or 3 solid phases
		k = asmCurrent(i+1)
		select case(fractl(k))
		case(0,1)	! phases are not fractionating or are growing and fractionating
			do j = 1,numPhCo(k)
				xPhCo(k,j) = pointPhaseComp(iSeg,iPointX,i,j)
				end do
		case(2)		! phase is being consumed - use fixed composition (i.e. chlorite)
			do j = 1,numPhCo(k)
				xPhCo(k,j) = xPhCoInitial(k,j)
				end do
		case default
			call FSS_Alert('ALERT','Error in EqSegment')
		end select		
		end do


	CALL AdjustAsm(1)		! we changed assemblage. Adjust volumes, moles etc.

!	The variance should be 2 (Duhem's theorem)
	sdel(1) = 0.0d0			! T
	sdel(2) = 0.0d0			! P
	smon(1) = 1
	smon(2) = 2
	!        P&T   (   number of dX terms          - 4 dependent dXterms   ) + previous dMterms
!      		set monitor parameters
	NSTEP=SNSTEP
	do i=1,nvar-numEqEquns		
		deltax(i)=sdel(i)	! these are all = 0 at this point
		mon(i)=smon(i)
		end do
!      		Set up array IPOINT to contain pointers to non-monitor parameters
	J=0
	Do i=1,NVAR
		DO L=1,NVAR-numEqEquns	
			IF(MON(L).eq.i)go to 1510
			end do
		J=J+1
		IPOINT(J) = I
1510		CONTINUE
		end do
	if(output.eq.1)then
		call PrinTT(1)
		write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NEQ)
		endif
!      	THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
	DO I=1,nvar
!         	ALLX(3,I) IS WHERE calculations STARTED-save
		ALLX(3,I)=ALLX(1,I)
		end do
	call SetTPX

!	This first call is required because it sets the phase composition to be that of the MDF criteria
	IZERO=0
	if(output.eq.1) then
		write(12,*)' Call to CalcAffinity(1) to get affinities'
		call CalcAffinity(1)	! print this initial calculation
		write(12,*)' '
		endif
	! store the initial moles for the phases
	MDFmoles(1,2) = mp0(k2)
	MDFmoles(1,3) = mp0(k3)
!	pause 'ready to call Compute4EQUIL -- Check the output window -- hit return when ready'
	call Compute4EQUIL(izero)
	if(izero.gt.0)then
		write(*,*)' Something went wrong in Compute4EQUIL '
		endif
		
	! we should now have the equilibrium assemblage

!	Store the final moles
	MDFmoles(2,2) = mp0(k2)
	MDFmoles(2,3) = mp0(k3)

	!store the new grain boundary composition
	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		GBComp(2,j) = xPhCo(1,j)
		end do
	if(output.eq.1)then
		call Printt(1)
		endif

!	At this point, we need to 
		!(a) save the GB composition and 
		!(b) store the phase compositions to examine zoning
		!(c) store the amount of each phase produced or consumed 
	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		pointComp(iSeg,iPointX,j) = xPhCo(1,j)
		end do
	do j = 1,numPhCo(k2)
		pointPhaseComp(iSeg,iPointX,1,j) = xPhCo(k2,j)
		end do
	do j = 1,numPhCo(k3)
		pointPhaseComp(iSeg,iPointX,2,j) = xPhCo(k3,j)
		end do

	! these are the total moles since the model began
	pointPhaseMoles(iSeg,iPointX,1) = pointPhaseMoles(iSeg,iPointX,1) + MDFmoles(2,2)-MDFmoles(1,2)
	pointPhaseMoles(iSeg,iPointX,2) = pointPhaseMoles(iSeg,iPointX,2) + MDFmoles(2,3)-MDFmoles(1,3)	

	! these are the change in moles for this iteration
	pointPhaseMolesDelta(iSeg,iPointX,1) = MDFmoles(2,2)-MDFmoles(1,2)
	pointPhaseMolesDelta(iSeg,iPointX,2) = MDFmoles(2,3)-MDFmoles(1,3)	

	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine CalculateGBChange(deltaMolesX,GBChangeX)
! *****************************************
	include "Assemb.inc"
!	include "Monit.inc"
!	include "Diffuse_GB_Hex.inc"
!	include "AutoAffinity.inc"
! *****************************************
	integer*4 i,ii,j,k,kcur,L
	real*8 deltaMolesX(4),GBTempX(20,10),GBChangeX(10)
	! now calculate the change in GB composition based on the above
	! zero out array that will hold changes in GB composition for each reaction
	! index "i" is the reaction
	! index "j" is the component
	! deltaMolesX(k) is the number of moles of solid phase to multiply by the reaction coefficient
	do i = 1,numEqRxn
		do J = 1,numPhCo(1)
			GBTempX(i,j) = 0.0d0
			end do
		end do

	i = 0
	ii = np - numEqRxn			! Starting index of the linearly independent reactions
	do kcur = 2,numPh			! loop through the solid phase components -- one reaction for each
		k = asmCurrent(kcur)
		do j = 1,numPhCo(k)		! phase components for the solid phase
			ii = ii + 1		! index of the independent reaction in array Arx
			i = i + 1		! index of the independent reaction starting at 1
			do L = 1,numPhCo(1)	! phase components for the GB phase
						! NOTE: the reaction coefficients (Arx array) are scaled to 1 oxygen/phase
				GBTempX(i,L) = Arx(ii,L)*deltaMolesX(kcur)*xPhCo(k,j)	
				end do
			end do
		end do


	! add up the changes
	do j = 1,numPhCo(1)
		sum = 0.0d0
		do i = 1,numEqRxn
			sum = sum + GBTempX(i,j)
			end do
		GBChangeX(j) = sum
		xPhCo(1,j) = xPhCo(1,j) + GBChangeX(j)	
		if(xPhCo(1,j).lt.0.0d0)then
			write(12,*)' Phase component is < 0. j = ',j,phCoName(1,j),xPhCo(1,j)
			pause 'Hit return to continue'
			endif
		end do
	! Now we have the composition of the grain boundary for the small increment of moles (deltaMolesdX)
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Save_and_Reset_Comps(iType)
      implicit none
!	Routine to save or reset the compositions of the phases in the current assemblage
!	itype = 1  Save current compositions
!	itype = 2  Reset current compositions to the last save

! *****************************************
	include "Assemb.inc"
! *****************************************
	integer*4 iType,kcur,k,j
	
	select case (iType)
	Case(1)		! save current composition
	do kcur = 1,numPh
		k = AsmCurrent(kcur)
		do j = 1,numPhCo(k)
			xPhCoLast(k,j) = xPhCo(k,j)	
			end do
		end do
	
	Case(2)		! reset current composition to last save
	do kcur = 1,numPh
		k = AsmCurrent(kcur)
		do j = 1,numPhCo(k)
			xPhCo(k,j) = xPhCoLast(k,j)	
			end do
		end do
	
	Case default
	
	end select
	
	return
	end



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE CalcAffinity(iprint)
      implicit none
!	Routine to calculate the affinity between the phases and the grain boundary
!	For each component of each phase
!	A = GatTP(phase component in phase) - GatTP(GB for that phase component)
!	The value of G(GB for that phase component) is calculated from the linearly independent reactions in array Arx:
!	  Linearly independent reactions:
!	     SiGB    AlGB    MgGB    FeGB    abQz    Prp     Alm 
!	   -1.000   0.000   0.000   0.000   1.000   0.000   0.000
!	   -3.000  -2.000  -3.000   0.000   0.000   1.000   0.000
!	   -3.000  -2.000   0.000  -3.000   0.000   0.000   1.000

!	iPrint = 1 then print calculations
!	iPrint = 0 do not print calculations
! *****************************************
	include "Assemb.inc"
	include "Output.inc"
! *****************************************
!      local variables
	integer*4 rxnCount,k,j,kcur,jj,ivol,izero,iPrint
	real*8 mu(MaxPC),muPhCo(phMax,phCoMax)
	real*8 Rjoules,AffTemp
      	DATA Rjoules/8.3144D0/


!	Calculate thermodynamic values at the new compositions	
!	ilong(3) = 1		! list out partial derivatives. It appears to be working fine
	TK=TC+273.15
	izero = 0
	ivol = 0
! 	CALL ALLKduTPX(ivol,izero)
! 	if(izero.gt.0)return
! 	TK = TC + 273.15d0
! 	call CalculateCPTemp(TK)
! 	Do kCur = 1,numPh
! 		k = asmCurrent(kCur)
! 		call duTPX(k,ivol,izero)
! 		if(izero.gt.0)return
! 		call dlnAdX(k,izero)
! 		if(izero.gt.0)return
! 		end do
!	ilong(3) = 0
!	Calculate µ of each system component
	jj = 0
	if(iPrint.eq.1)then
		write(12,*)' '
		write(12,*)'----------------------------------'
		write(12,*)'mu = chemical potential'
		endif
	do kcur = 1,numPh
		k = AsmCurrent(kcur)
		call duTPX(k,ivol,izero)
		if(izero.gt.0)return
		call dlnAdX(k,izero)
		if(izero.gt.0)return
		do j = 1,numPhCo(k)
		!	uZero(k,j) = gattp(k,j)
		!	mu(k,j) = gattp(k,j) + Rjoules*TK*lnAct(k,j)
			jj = jj + 1
			mu(jj) = gattp(k,j) + Rjoules*TK*lnAct(k,j)		! These are stored in the same order as array Arx
			mu(jj) = mu(jj)/numOxygens(k)				! scale the µ to the number of oxygens
			muPhCo(k,j) = mu(jj)					! redundant? Ever used?? 
			if(iPrint.eq.1)then
				write(12,101)minrec(k),PhName(k),PhCoName(k,j),xPhCo(k,j),mu(jj),muPhCo(k,j)
				endif
			end do
		end do
		

	if(iPrint.eq.1)then
! 		write(12,*)' '
! 		write(12,*)' Gsystem = ',gSystem
		write(12,*)' '
		write(12,*)'----------------------------------'
		write(12,*)'Affinities'
		endif
!	calculate the affinity for each phase component from the linearly independent equations (see example above)
!	There is one independent equation for each phase component except the GB phase
!	ii = np-nrx      !note that reaction coefficients are stored in the bottom of array ARX
	!rxnCount = np - nrx
	rxnCount = np - numEqRxn
	do kcur = 2,numPh		! skip the first phase (it's the GB)
		k = AsmCurrent(kcur)
		do j = 1,numPhCo(k)
			AffTemp = 0.0D0
			rxnCount = rxnCount + 1
			do jj = 1,np			! loop through all of the phase components including GB ones
			!	we can loop through all phase components because the all entries for the other phase components are zero
				AffTemp = Afftemp + Arx(rxnCount,jj)*mu(jj)
				end do
! 			Aff(k,j) = AffTemp/numOxygens(k)
			Aff(k,j) = AffTemp
			if(iPrint.eq.1)then
				write(12,101)minrec(k),PhName(k),PhCoName(k,j),xPhCo(k,j),Aff(k,j)
				endif
		101	Format(I5,2x,A8,2x,A4,2x,F12.5,3F12.1)	
			end do
		end do

	if(iPrint.eq.1)then
		write(12,*)'----------------------------------'
		endif
99	continue
	return
	end	


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine SetPointToNodePhaseComps(iunit)
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	integer*4 iNode,iSeg,i,j,MIFID,MIFseg,myPoint,imin,kk,jj,iunit
! --------------------------------------------------------------------

	! set the point phase comp for the point on top of the node to the node phase comp (they are equivalent)
	do iNode = 1,numNodes
! 		do i = 1,3		! cycle through 3 segments that link to this node
		do i = 1,numNodeSegs(iNode)		! cycle through 3 segments that link to this node
			iSeg = nodeSegConnect(iNode,i)
			if(numSegReactionPhases(iSeg).eq.0)cycle		! nothing to set on this segment (only options are 0 and 2)

			if(NodePointOnTop(iNode,i).eq.numPointStart(iSeg))then	! see if the node connects to the start or the end of the seg
				myPoint = numPointStart(iSeg)
				else
				myPoint = numPointEnd(iSeg)
				endif
				
			do kk = 1,2		! cycle through the 2 phases along the segment
				MIFseg = segMIFID(iSeg,kk)

! 		numNodeReactionPhases(maxNodes),	& ! the number of different reacting phases at a node
! 		nodeReactionPhases(maxNodes,3),		& ! index of the reacting phases around the node 
! 							  ! if numNodeReactionPhases = 0 this is irrelevent
! 							  ! if numNodeReactionPhases = 3, then the indicies are 1, 2, 3
! 							  ! if numNodeReactionPhases=2 the indicies are either 1&2, 1&3 or 2&3


				if(numNodeReactionPhases(iNode).eq.0)go to 10
				do jj = 1,numNodeReactionPhases(iNode)	! cycle through the different phases at the node
					imin = nodeReactionPhases(iNode,jj)
					MIFID = nodeMIFID(iNode,imin)
					if(MIFID.eq.MIFseg)then
						do j = 1,numPhCo(MIFID)
						     pointPhaseComp(iSeg,myPoint,kk,j) = nodePhaseComp(iNode,imin,j)
						     end do
						go to 10
						endif
					end do
				write(*,*)' No match found in Sub SetPointToNodePhaseComps()'
				write(*,101)iNode,i,(NodeMIFID(iNode,nodeReactionPhases(iNode,jj)),jj=1,numNodeReactionPhases(iNode))
				write(*,101)iSeg,kk,(segMIFID(iSeg,jj),jj=1,2)
				write(iunit,*)' No match found in Sub SetPointToNodePhaseComps()'
				write(iunit,101)iNode,i,(NodeMIFID(iNode,nodeReactionPhases(iNode,jj)),jj=1,numNodeReactionPhases(iNode))
				write(iunit,101)iSeg,kk,(segMIFID(iSeg,jj),jj=1,2)
101				format(8I8)
				!call fss_alert('Ach','Ach....no match found')
				pause 'Hit return to continue'
				return

10				continue
				end do		! end kk loop over 2 phases along segment
			end do		! end loop over 3 segments
		end do 		! end loop over nodes
			

	return
	end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine SetPointToNodePhaseComps_old()
!	This routine is no longer used-- it can be deleted

	implicit none
! *****************************************
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	integer*4 iNode,iSeg,i,j,k,MIFID,MIFseg,ptStart,ptEnd,imin,kk
! --------------------------------------------------------------------

	! set the point phase comp for the point on top of the node to the node phase comp (they are equivalent)
	do iSeg = 1,numSegs
		if(numSegReactionPhases(iSeg).eq.0)cycle
		do k = 1,2		!cycle through the 2 nodes that this seg connects to
			iNode = segNodes(iSeg,k)
			do kk = 1,2	! cycle through the 2 phases along this seg
				MIFseg = segMIFID(iSeg,kk)
				do i = 1,3 	! cycle through the 3 segs that connect to the node to find the right one
					if(NodeSegConnect(iNode,i).eq.iSeg)then	!This is the correct segment for this node
						if(NodePointOnTop(iNode,i).eq.numPointStart(iSeg))then	! see if the node connects to the start or the end of the seg
							ptStart = numPointStart(iSeg)
							!do imin = 1,3	! cycle through the 3 phases at the node
							do imin = 1,numNodeReactionPhases(iNode)	! cycle through the different phases at the node (2 or 3)
								MIFID = nodeMIFID(iNode,imin)
								if(MIFID.eq.MIFseg)then
									do j = 1,numPhCo(MIFID)
									     pointPhaseComp(iSeg,ptStart,kk,j) = nodePhaseComp(iNode,imin,j)
									     end do
									go to 10
									endif
								end do
							else			! the node connects to the end of the seg
							ptEnd = numPointEnd(iSeg)
							do imin = 1,3		! cycle through the 3 phases at the node
								MIFID = nodeMIFID(iNode,imin)
								if(MIFID.eq.MIFseg)then
									do j = 1,numPhCo(MIFID)
									     pointPhaseComp(iSeg,ptEnd,kk,j) = nodePhaseComp(iNode,imin,j)
									     end do
									go to 10
									endif
								end do
							endif
						endif
					end do	! end i loop over the 3 phases at the node
10				continue
				end do	! end the kk loop over the 2 phases at the segment
			end do	! end k loop over the 2 nodes
		end do		! end iSeg loop


	return
	end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE CheckMassBalance(iunit)
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
! *****************************************
	integer*4 iNode,iSeg,iPointX,i,MIFID,iRx,iunit
	real*8 MPhSum,MPhDeltaSum

	write(iunit,*)' '
	write(iunit,*)'==================================================================='
	write(iunit,*)' Checking mass balance '
	write(iunit,*)'-------Nodes--------'
	write(iunit,*)'                  delta moles sum     moles sum'
	do iNode = 1,numNodes
		!write(*,*)iNode,numNodeReactionPhases(iNode)
		MPhSum = 0.0D0
		MPhDeltaSum = 0.0d0
		do i = 1,numNodeReactionPhases(iNode)
			iRx = nodeReactionPhases(iNode,i)
			MIFID = nodeMIFID(iNode,iRx)
! 			MPhSum      = MphSum      + nodePhaseMoles(iNode,iRx)*numOxygens(MIFID)
! 			MPhDeltaSum = MphDeltaSum + nodePhaseMolesDelta(iNode,iRx)*numOxygens(MIFID)
			MPhSum      = MphSum      + nodePhaseMoles(iNode,iRx)
			MPhDeltaSum = MphDeltaSum + nodePhaseMolesDelta(iNode,iRx)
			end do
		write(iunit,101) iNode,numNodeReactionPhases(iNode),					&
			(phName(nodeMIFID(iNode,nodeReactionPhases(iNode,i))),i=1,numNodeReactionPhases(iNode))	
101		format(2I5,3(2x,A8))
		write(iunit,*)'                          ',MPhDeltasum,MPhsum
		end do
	write(iunit,*)'-------Segs--------'
	write(iunit,*)'        Point        delta moles sum     moles sum'
	do iSeg = 1,numSegs
		if(numSegReactionPhases(iSeg).eq.2)then
			MPhSum = 0.0D0
			MPhDeltaSum = 0.0d0
			write(iunit,102)iSeg,(phName(segMIFID(iSeg,i)),i=1,2)
102			format(I5,3(2x,A8))
			do iPointX = numPointStart(iSeg)+1,numPointEnd(iSeg)-1		! skip the points on the nodes
				do i = 1,2
					MIFID = segMIFID(iSeg,i)
! 					MPhsum = MPhsum + pointPhaseMoles(iSeg,iPointX,i)*numOxygens(MIFID)
! 					MPhDeltaSum = MPhDeltaSum + pointPhaseMolesDelta(iSeg,iPointX,i)*numOxygens(MIFID)
					MPhsum = MPhsum + pointPhaseMoles(iSeg,iPointX,i)
					MPhDeltaSum = MPhDeltaSum + pointPhaseMolesDelta(iSeg,iPointX,i)
					end do
				write(iunit,*)'                   ',iPointX,MPhDeltasum,MPhsum
				end do
			endif
		end do


	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PhaseMassSum()
      implicit none
!	This routine will sum up the change in moles for a group of nodes and segments
!	The goal is to calculate the amount that phases have reacted at various parts of the grid
! *****************************************
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
! *****************************************
	integer*4 iNode,iSeg,iPointX,MIFID,saywhat
	integer*4 i,k,kk,numNodesToSurvey,numSegsToSurvey,NodeToSum(20),segToSum(20)
	real*8 sumMoles(20)

1	continue
	write(*,*)' Sum moles around nodes and segments'
	write(*,*)' 0 = return'
	write(*,*)' 1 = calculate'
	read(*,*)saywhat
	if(saywhat.eq.0)return

	numNodesToSurvey = 0
10	continue
	write(*,*)'Input sequence of nodes to survey. Enter 0 when done'
	read(*,*)i
	if(i.gt.0)then
		numNodesToSurvey= numNodesToSurvey + 1
		nodeToSum(numNodesToSurvey) = i
		go to 10
		endif	

	numSegsToSurvey = 0
11	continue
	write(*,*)'Input sequence of segments to survey. Enter 0 when done'
	read(*,*)i
	if(i.gt.0)then
		numSegsToSurvey= numSegsToSurvey + 1
		segToSum(numSegsToSurvey) = i
		go to 11
		endif	

	do i = 1,numPhMIF
		sumMoles(i) = 0.0d0	! this is where we store the sum of the moles for each phase -- indexed to the MIF
		end do

	do i = 1,numNodesToSurvey
		iNode = nodeToSum(i)
		do kk = 1,numNodeReactionPhases(iNode)
			k = nodeReactionPhases(iNode,kk)
			MIFID = nodeMIFID(iNode,k)
			sumMoles(MIFID) = sumMoles(MIFID) + nodePhaseMoles(iNode,k)
			end do
		end do

	do i = 1,numSegsToSurvey
		iSeg = segToSum(i)
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
			do k = 1,numSegReactionPhases(iSeg)
				MIFID = segMIFID(iSeg,k)
				sumMoles(MIFID) = sumMoles(MIFID) + pointPhaseMoles(iSeg,iPointX,k)
				end do
			end do
		end do


	write(12,*)' '
	write(12,*)'==================================================================='
	write(12,*)' Mass Sums '
	write(12,103)(i,PhName(i),sumMoles(i),i=2,numPhMIF)
103	format(i5,A10,E12.4)
	

	go to 1

	end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE CheckNodeInformation()
	implicit none
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	integer*4 iNode,i,j,k,kk
	integer*4 iSeg,iPoint,MIFID
	

35		continue
		write(*,*)'  '
		write(*,*)'  '
		write(*,*)'Input node to check. 0 to exit'
		read(*,*)iNode
		if(iNode.eq.0)return
		write(*,*)'Node array information'
		write(*,*)' Node XY           ',nodeX(iNode),nodeY(iNode)
		write(*,*)' numNodeSegs       ',numNodeSegs(iNode)
		write(*,*)' nodeNodeConnect   ',(nodeNodeConnect(iNode,j),j=1,numNodeSegs(iNode))
		write(*,*)' nodeSegConnect    ',nodeSegConnect(iNode,1),nodeSegConnect(iNode,2),nodeSegConnect(iNode,3)
		write(*,*)' NodePointOnTop     ',nodePointOnTop(iNode,1),nodePointOnTop(iNode,2),nodePointOnTop(iNode,3)
		write(*,*)' NodePointNextTo '
		do i = 1,3
			iSeg = nodeSegConnect(iNode,i)
			iPoint = nodePointNextTo(iNode,i)
			write(*,*)iSeg,iPoint,PointX(iseg,iPoint),PointY(iSeg,iPoint)
			end do
		write(*,*)'------------------------'
		write(*,*)' Node Phases (should be anticlockwise)'
		!write(*,*)' NodeCrystalIndex    ',nodeCrystalIndex(iNode,1),nodeCrystalIndex(iNode,2),nodeCrystalIndex(iNode,3)
		write(*,*)' i  MIFID    Phase name         '
		do i = 1,3
			write(*,*)i,nodeMIFID(iNode,i),PhName(nodeMIFID(iNode,i))
			end do
		write(*,*)'------------------------'
		write(*,*)' iSeg,segMIFID '
		do i = 1,3
			iSeg = nodeSegConnect(iNode,i)
			write(*,*)iSeg,segMIFID(iSeg,1),segMIFID(iSeg,2)
			end do
		write(*,*)'------------------------'
		write(*,*)' numNodeReactionPhases(iNode) ',numNodeReactionPhases(iNode)
		do i = 1,numNodeReactionPhases(iNode)
			k = nodeReactionPhases(iNode,i)
			j = nodeMIFID(iNode,k)
			write(*,*)k,nodeMIFID(iNode,k),phName(j)	
			end do
		k = nodeSegTheSamePhases(iNode)
		write(*,*)' nodeSegTheSamePhases(iNode)  ',k			
		if(k.gt.0)write(*,*)' Segment number ... ',nodeSegConnect(iNode,k)	! k will be 0 if all 3 phases are the same
		write(*,*)'------------------------'
		write(*,*)' Grain boundary'
		MIFID = 1		! this is the grain boundary
		kk = 0
		write(*,*)MIFID,PhName(MIFID),kk,kk,kk,		&
				numPhCo(MIFID),(phCoName(MIFID,j),nodeComp(iNode,j),j=1,numPhCo(MIFID))


		go to 35

	end	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine OpenOutputFile(status)
	implicit none
	integer*4 status
! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
!	include "MIFarrays.inc"
! --------------------------------------------------------------------
!	open(73,FILE='',status='UNKNOWN',iostat=status)
	open(73,FILE='',status='NEW',iostat=status,action='WRITE')
	if(status.ne.0)return
	inquire(73,NAME=ModelOutputFileBase)
	write(73,*)' This file is intentionally left empty'
	write(73,*)' Do not delete -- it is used as the base file name for creating animations'
	close(73)
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE WriteAllHeader(iunit)
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Diffuse_GB_Hex.inc"
! *****************************************
	integer*4 iunit,i
	character*16 date,time,zone
	integer timevalues(8)

	write(iunit,*)' GibbsMDF output file'				! extra line to preserve compatibility with PseudoForwardModel routine

	!date_and_time (date, time, zone, values)
	!subroutine date_and_time(date,time,zone,values) ! returns date and
							! time information
	!character(*), optional, intent(out)) :: date 	! date in CCYYMODD
							! format
	!character(*), optional, intent(out)) :: time 	! time in
							! HHMMSS.SSS
							! format
	!character(*), optional, intent(out)) :: zone 	! zone in ±HHMM
							! format
	!integer, optional, intent(out)) :: values(:) 	! date, time, and
							! zone in integer
							! form
	!end subroutine
	!year    = timevalues(1)
	!month   = timevalues(2)
	!day     = timevalues(3)
	!zone    = timevalues(4)
	!hours   = timevalues(5)
	!minutes = timevalues(6)
	!seconds = timevalues(7)
	!milliseconds = timevalues(8)
	!write(*,*)'hours, minutes, seconds ',hours,minutes,seconds
	call date_and_time(date,time,zone,timevalues)
	
	
	write(iunit,*)'***************************************'
	write(iunit,*)'Model file      = ',modelFile
	write(iunit,*)'MasterInputFile = ',FileIn
	write(iunit,*)'***************************************'
	write(iunit,101)(timevalues(i),i=1,3),(timevalues(i),i=5,7)
101	format('Year ',I6,'  Month ',I3,'  Day ',I3,' Hours ',I3,' Minutes ',I3,' Seconds ',I3)	
	return
	end	

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE WriteAllNodesAndSegments(iunit)
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
! *****************************************
	integer*4 iunit,iNode,k,MIFID,j,iSeg,iPointX,kk,i,PCN
	real*8 GBmoles
	
	
	write(iunit,*)'-----------------------------'
	write(iunit,94)numDiffIterations
94	format(I10,'    Diffusion iterations')
	write(iunit,95)totalCycles
95	format(I10,'    Cycle number')
!	write(iunit,*)'            Phase      mMoles phase          Phase composition'

	write(iunit,*)'-----------------------------'
	!write(iunit,*)'Nodes with 3 phases'
	write(iunit,*)'Nodes        Phase      ∆mMoles phase   mMoles phase          Phase composition'
	write(iunit,101)numNodes
101	format(I6,'    numNodes ')
	GBmoles = 0.0d0

! 	do iNode = 1,numNodes
! 		write(iunit,102)iNode,nodeX(iNode),nodeY(iNode),numNodeReactionPhases(iNode)
! 102		format(I6,2F12.3,I5,'         Node number,  X, Y    numNodePhases')
! 		MIFID = 1		! this is the grain boundary
! 		write(iunit,103)MIFID,PhName(MIFID),GBmoles,GBmoles,		&
! 				numPhCo(MIFID),(phCoName(MIFID,j),nodeComp(iNode,j),j=1,numPhCo(MIFID))
! 		do kk = 1,numNodeReactionPhases(iNode)
! 			k = nodeReactionPhases(iNode,kk)
! 			MIFID = nodeMIFID(iNode,k)
! 			write(iunit,103)MIFID,PhName(MIFID),nodePhaseMolesDelta(iNode,k),nodePhaseMoles(iNode,k),	&
! 				numPhCo(MIFID),(phCoName(MIFID,j),nodePhaseComp(iNode,k,j),j=1,numPhCo(MIFID))
! 103			format(T8,I4,2x,A8,2E15.5,I8,20(5x,A8,F12.5))
! 			end do

	do iNode = 1,numNodes
		write(iunit,102)iNode,nodeX(iNode),nodeY(iNode)
102		format(I6,2F12.3,'         Node number,  X, Y')
		write(iunit,108)(nodeMIFID(iNode,i),i=1,3)
108		format(T8,3I4,T26,'nodeMIFID')
		select case(numNodeReactionPhases(iNode))
			case(0)
				write(iunit,130)numNodeReactionPhases(iNode)
	130			format(T8,I4,T26,'numNodeReactionPhases, index')
			case(1)
				write(iunit,131)numNodeReactionPhases(iNode),(nodeReactionPhases(iNode,k),k=1,numNodeReactionPhases(iNode))
	131			format(T8,2I4,T26,'numNodeReactionPhases, index')
			case(2)
				write(iunit,132)numNodeReactionPhases(iNode),(nodeReactionPhases(iNode,k),k=1,numNodeReactionPhases(iNode))
	132			format(T8,3I4,T26'numNodeReactionPhases, index')
			case(3)
				write(iunit,133)numNodeReactionPhases(iNode),(nodeReactionPhases(iNode,k),k=1,numNodeReactionPhases(iNode))
	133			format(T8,4I4,T26,'numNodeReactionPhases, index')
			case default
			write(*,*)' Bad numNodeReactionPhases(iNode) iNode = ',iNode
			pause 'Hit return to continue'
			end select
		MIFID = 1		! this is the grain boundary
		select case (outputType)
		case(1)		! no affinity
			write(iunit,103)MIFID,PhName(MIFID),GBmoles,GBmoles,		&
					numPhCo(MIFID),(phCoName(MIFID,j),nodeComp(iNode,j),j=1,numPhCo(MIFID))
			do kk = 1,numNodeReactionPhases(iNode)
				k = nodeReactionPhases(iNode,kk)
				MIFID = nodeMIFID(iNode,k)
				write(iunit,103)MIFID,PhName(MIFID),nodePhaseMolesDelta(iNode,k),nodePhaseMoles(iNode,k),	&
					numPhCo(MIFID),(phCoName(MIFID,j),nodePhaseComp(iNode,k,j),j=1,numPhCo(MIFID))
	103			format(T8,I4,2x,A8,2E15.5,I8,20(5x,A8,F12.5))
				end do
		case(2)		! with affinity
			write(iunit,1032)MIFID,PhName(MIFID),GBmoles,GBmoles,GBmoles,		&
					numPhCo(MIFID),(phCoName(MIFID,j),nodeComp(iNode,j),j=1,numPhCo(MIFID))
			do kk = 1,numNodeReactionPhases(iNode)
				k = nodeReactionPhases(iNode,kk)
				MIFID = nodeMIFID(iNode,k)
				write(iunit,1032)MIFID,PhName(MIFID),nodePhaseMolesDelta(iNode,k),nodePhaseMoles(iNode,k),	&
					nodePhaseAffinity(iNode,k),								&
					numPhCo(MIFID),(phCoName(MIFID,j),nodePhaseComp(iNode,k,j),j=1,numPhCo(MIFID))
	1032			format(T8,I4,2x,A8,3E15.5,I8,20(5x,A8,F12.5))
				end do
		case default
		end select


		write(iunit,140)numNodeSegs(iNode)
	140	format(I11,'          NumNodeSegs')
		write(iunit,141)(nodeNodeConnect(iNode,j),j=1,numNodeSegs(iNode))
	141	format(10x,3I10'    nodeNodeConnect')
		write(iunit,142)(nodeSegConnect(iNode,j),j=1,numNodeSegs(iNode))
	142	format(10x,3I10,'    nodeSegConnect')
		write(iunit,143)(nodePointOnTop(iNode,j),j=1,numNodeSegs(iNode))
	143	format(10x,3I10,'    nodePointOnTop')
		write(iunit,144)(nodePointNextTo(iNode,j),j=1,numNodeSegs(iNode))
	144	format(10x,3I10,'    nodePointNextTo')
	

		end do

	write(iunit,*)'-----------------------------'
	write(iunit,*)'Segments'
	write(iunit,105)numSegs
105	format(I6,'    numSegs')
	do iSeg = 1,numSegs
		!iSeg = segReactionPhases(i)
!		write(iunit,107)iSeg,numPoints(iSeg),segNodes(iSeg,1),segNodes(iSeg,2),numSegReactionPhases(iSeg)
!107		format(5I5,'     iSeg, numPoints(iSeg), NodeStart, NodeEnd   numSegReactionPhases ')
		write(iunit,107)iSeg,numPointStart(iSeg),numPointEnd(iSeg),segNodes(iSeg,1),segNodes(iSeg,2),numSegReactionPhases(iSeg)
107		format(6I6,'     iSeg, numPointStart,  numPointEnd, NodeStart, NodeEnd   numSegReactionPhases ')
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
			write(iunit,106)iSeg,iPointX,pointX(iSeg,iPointX),pointY(iSeg,iPointX)
106			format(2I8,2F12.3,'      Segment, point, X, Y ')
			MIFID = 1		! this is the grain boundary
			write(iunit,124)MIFID,PhName(MIFID),GBmoles,GBmoles,		&
					numPhCo(MIFID),(phCoName(MIFID,j),pointComp(iSeg,iPointX,j),j=1,numPhCo(MIFID))
124			format(T8,I4,2x,A8,2E15.5,I8,20(5x,A8,F12.5))
			!do k = 1,2
			do k = 1,numSegReactionPhases(iSeg)
				MIFID = segMIFID(iSeg,k)
				write(iunit,103)MIFID,PhName(MIFID),pointPhaseMolesDelta(iSeg,iPointX,k),		&
					pointPhaseMoles(iSeg,iPointX,k),						&
					numPhCo(MIFID),(phCoName(MIFID,j),pointPhaseComp(iSeg,iPointX,k,j),j=1,numPhCo(MIFID))
				end do
			end do
		end do

	! Now write out the grain boundary composition for every segment,point
	
		!pointX(maxSegs,maxPoints)
		!pointY(maxSegs,maxPoints)
		!pointPhaseComp(maxSegs,maxPoints,2,maxPhCo),	&! Composition of each of 2 phases at a reaction segment

	write(iunit,*)'-----------------------------'
	write(iunit,*)'Grain Boundary composition'
	write(iunit,111)(elName(j),j=1,numEl)
111	format('   Seg    Point    X       Y           ',20(A3,9x))
	write(iunit,110)numSegs
110	format(I5,'    Number of segments')
	do iSeg = 1,numSegs
		write(iunit,109)iSeg,numPointStart(iSeg),numPointEnd(iSeg)
109		format(4x,3I5,'     iSeg,numPointStart, numPointEnd')
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
!		do iPointX = 1,numPoints(iSeg)
			write(iunit,118)iSeg,iPointX,pointX(iSeg,iPointX),pointY(iSeg,iPointX),		&
				(pointComp(iSeg,iPointX,j),j=1,numEl)
118			format(2I8,2F12.3,20F12.5)
			end do
		end do

	write(iunit,*)'++++++Seg length, shorter, exclusions ++++++++++++++++++++++++++++++++'
	do iSeg = 1,numSegs
		write(iunit,*)iSeg,segLength(iSeg),segGettingShorter(iSeg),segCaptureExclusion(iSeg)
		end do		

!	Write out crystal information
	write(iunit,*)'++++++Crystal Information ++++++++++++++++++++++++++++++++'
	Do i = 1,numPhases
		write(iunit,114)PhaseName(i)
	114	format(A8)
		do j = 1,numPhaseCrystals(i)			! this will loop on a specific phase type
			PCN = PhaseCrystalNumber(i,j)		! Phase Crystal Number
			write(iunit,113)PCN
		113	format(I8,'      Phase crystal number')
			write(iunit,115)numCrystalNodes(PCN)
		115	format(I8,'      Number of crystal nodes (numCrystalNodes)')
			do k = 1,numCrystalNodes(PCN) 
! 				iSeg = CrystalSegs(PCN,k)
! 				iNode = CrystalNodes(PCN,k)
				write(iunit,112)CrystalNodes(PCN,k),CrystalSegs(PCN,k)
		112		format(3I8)
				end do
			end do
		end do
	
	return
	end
	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine OpenModelOutputFile(isFileOpen2,readXY)
!	Open and read the output file for a specific model step
      	implicit none
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Diffuse_GB_Hex.inc"
	integer*4 status,i,j,k,kk,l,iNode,MIFID,iPointX,iSeg,isFileOpen2,node1,node2,readXY,PCN
	integer*4 MIFID1,MIFID2,MIFID3
	character*8 dummy
	character*1 char1,char2
	real*8 GBmoles,x,y

	if(isFileOpen2.eq.0)then		! we need to open the file
		open(40,file='',status = 'OLD',iostat = status)
		if(status.ne.0)then
			call FSS_Alert('Alert','Problem opening model output file')
			return
			endif
		INQUIRE(40, NAME=modelOutputFile)
		write(*,*)'ModelOutputFile name = ',modelOutputFile
		endif
	read(40,*)dummy		!Line 1
	read(40,*)dummy		!Line 2 -- stars
	read(40,*)dummy		!Line 3 -- Model file name
	read(40,*)dummy		!Line 4 -- MIF file name
	read(40,*)dummy		!Line 5 -- stars
	read(40,*)dummy		!Line 6 -- Date and time
	read(40,*)dummy		!Line 7 -- dashes
	read(40,*)numDiffIterations	!Line 8 -- diffusion iterations
	read(40,*)totalCycles	!Line 9 -- Cycle number
	read(40,*)dummy		!Line 10 -- Dashes
	read(40,*)dummy		!Line 11 -- Column header for nodes
!	read(40,*)dummy		!Line 12 -- text 'Nodes'
	

	read(40,*)numNodes
	do i = 1,numNodes
		if(readXY.eq.0)then
			read(40,*)iNode,X,Y	! do not read the node XY (keep the previous positions)
			else
			read(40,*)iNode,nodeX(iNode),nodeY(iNode)
			endif
		if(iNode.ne.i)then
			call fss_Alert('ALERT','Problem reading output file')
			write(*,*)'i , iNode ',i,iNode
			pause 'hit return to continue'
			endif
		read(40,*)(nodeMIFID(iNode,kk),kk=1,3)
		read(40,*)numNodeReactionPhases(iNode),(nodeReactionPhases(iNode,kk),kk=1,numNodeReactionPhases(iNode))

		select case(outputType)
		case(1)		! no affinity
			! this first one is the grain boundary
			read(40,*)MIFID,PhName(MIFID),GBmoles,GBmoles,		&
					numPhCo(MIFID),(phCoName(MIFID,j),nodeComp(iNode,j),j=1,numPhCo(MIFID))
			do kk = 1,numNodeReactionPhases(iNode)
				k = nodeReactionPhases(iNode,kk)
				read(40,*)MIFID,PhName(MIFID),nodePhaseMolesDelta(iNode,k),nodePhaseMoles(iNode,k),		&
					numPhCo(MIFID),(phCoName(MIFID,j),nodePhaseComp(iNode,k,j),j=1,numPhCo(MIFID))
				nodeMIFID(iNode,k) = MIFID
				end do
		case(2)		! with affinity
			read(40,*)MIFID,PhName(MIFID),GBmoles,GBmoles,GBmoles,		&
					numPhCo(MIFID),(phCoName(MIFID,j),nodeComp(iNode,j),j=1,numPhCo(MIFID))
			do kk = 1,numNodeReactionPhases(iNode)
				k = nodeReactionPhases(iNode,kk)
				read(40,*)MIFID,PhName(MIFID),nodePhaseMolesDelta(iNode,k),nodePhaseMoles(iNode,k),		&
					nodePhaseAffinity(iNode,k),								&
					numPhCo(MIFID),(phCoName(MIFID,j),nodePhaseComp(iNode,k,j),j=1,numPhCo(MIFID))
				nodeMIFID(iNode,k) = MIFID
				end do
		case default
		end select

		read(40,*)numNodeSegs(iNode)
		read(40,*)(nodeNodeConnect(iNode,j),j=1,numNodeSegs(iNode))
		read(40,*)(nodeSegConnect(iNode,j),j=1,numNodeSegs(iNode))
		read(40,*)(nodePointOnTop(iNode,j),j=1,numNodeSegs(iNode))
		read(40,*)(nodePointNextTo(iNode,j),j=1,numNodeSegs(iNode))
		
		end do

	do iNode = 1,numNodes

		call AssignReactionPhases(iNode)
		go to 444

		MIFID1 = nodeMIFID(iNode,1)
		MIFID2 = nodeMIFID(iNode,2)
		MIFID3 = nodeMIFID(iNode,3)
		if(MIFID1.eq.MIFID2.and.MIFID1.eq.MIFID3)then	! all phases are the same
			numNodeReactionPhases(iNode) = 0
			go to 20
			endif
		if(MIFID1.ne.MIFID2.and.MIFID1.ne.MIFID3.and.MIFID2.ne.MIFID3)then	! all 3 phases different
			numNodeReactionPhases(iNode) = 3
			nodeReactionPhases(iNode,1) = 1
			nodeReactionPhases(iNode,2) = 2
			nodeReactionPhases(iNode,3) = 3
			nodeSegTheSamePhases(iNode) = 0	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
			go to 20
			endif	
		! if here, then 2 phases must be the same
		numNodeReactionPhases(iNode) = 2
		if(MIFID1.eq.MIFID2)then
			nodeReactionPhases(iNode,1) = 1
			nodeReactionPhases(iNode,2) = 3
			do k = 1,3
				iSeg = nodeSegConnect(iNode,k)
				!write(12,*)iNode,iSeg,segMIFID(iSeg,1),segMIFID(iSeg,2)
				if(segMIFID(iSeg,1).eq.SegMIFID(iSeg,2))then
					nodeSegTheSamePhases(iNode) = k	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
					go to 20
					endif
				end do
! 			go to 20
			endif
		if(MIFID1.eq.MIFID3)then
			nodeReactionPhases(iNode,1) = 2
			nodeReactionPhases(iNode,2) = 3
			do k = 1,3
				iSeg = nodeSegConnect(iNode,k)
				!write(12,*)iNode,iSeg,segMIFID(iSeg,1),segMIFID(iSeg,2)
				if(segMIFID(iSeg,1).eq.SegMIFID(iSeg,2))then
					nodeSegTheSamePhases(iNode) = k	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
					go to 20
					endif
				end do
! 			go to 20
			endif
		if(MIFID2.eq.MIFID3)then
			nodeReactionPhases(iNode,1) = 1
			nodeReactionPhases(iNode,2) = 2
			do k = 1,3
				iSeg = nodeSegConnect(iNode,k)
				!write(12,*)iNode,iSeg,segMIFID(iSeg,1),segMIFID(iSeg,2)
				if(segMIFID(iSeg,1).eq.SegMIFID(iSeg,2))then
					nodeSegTheSamePhases(iNode) = k	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
					go to 20
					endif
				end do
! 			go to 20
			endif
		! if here, we have a problem
		write(12,*)' In routine OpenModelOutputFile'
		write(12,*)' Failure to determine nodeReactionPhases and nodeSegTheSamePhases'
		write(12,*)' iNode = ',iNode
		write(12,227)iNode,(nodeMIFID(iNode,j),j=1,numCrystalsAtNode(iNode)),			&
				(nodeReactionPhases(iNode,j),j=1,numNodeReactionPhases(iNode)),nodeSegTheSamePhases(iNode)
	227	format(20I8)
	
	20	continue	

444		continue	! this is just a dummy jump to see if new subroutine works

		end do





				
	read(40,*)dummy		! dashes
	read(40,*)dummy		! 'Segments'
	read(40,*)numSegs
	do i = 1,numSegs
!		read(40,*)iSeg,numPoints(iSeg),segNodes(iSeg,1),segNodes(iSeg,2),numSegReactionPhases(iSeg)
		read(40,*)iSeg,numPointStart(iSeg),numPointEnd(iSeg),segNodes(iSeg,1),segNodes(iSeg,2),numSegReactionPhases(iSeg)
		!segReactionPhases(i) = iSeg
!		do L = 1,numPoints(iSeg)
		do L = numPointStart(iSeg),numPointEnd(iSeg)
			if(readXY.eq.0)then
				read(40,*)iSeg,iPointX,X,Y
				else
				read(40,*)iSeg,iPointX,pointX(iSeg,iPointX),pointY(iSeg,iPointX)
				endif
			read(40,*)MIFID,PhName(MIFID),GBmoles,GBmoles,		&
					numPhCo(MIFID),(phCoName(MIFID,j),pointComp(iSeg,iPointX,j),j=1,numPhCo(MIFID))
			do k = 1,numSegReactionPhases(iSeg)
				read(40,*)MIFID,PhName(MIFID),pointPhaseMolesDelta(iSeg,iPointX,k),		&
					pointPhaseMoles(iSeg,iPointX,k),					&
					numPhCo(MIFID),(phCoName(MIFID,j),pointPhaseComp(iSeg,iPointX,k,j),j=1,numPhCo(MIFID))
				segMIFID(iSeg,k) = MIFID
				end do
			end do
		segX(iSeg,1) = pointX(iSeg,numPointStart(iSeg))
		segY(iSeg,1) = pointY(iSeg,numPointStart(iSeg))
		segX(iSeg,2) = pointX(iSeg,numPointEnd(iSeg))
		segY(iSeg,2) = pointY(iSeg,numPointEnd(iSeg))

		end do
		
	read(40,*)dummy		! dashes
	read(40,*)dummy		! read header 'Grain boundary composition'
	read(40,*)dummy		! read header lines

	read(40,*)numSegs
	do i = 1,numSegs
!		read(40,*)iSeg,numPoints(iSeg)
		read(40,*)iSeg,numPointStart(iSeg),numPointEnd(iSeg)
!		do k = 1,numPoints(iSeg)
		do k = numPointStart(iSeg),numPointEnd(iSeg)
			if(readXY.eq.0)then
				read(40,*)iSeg,iPointX,X,Y,		&
					(pointComp(iSeg,iPointX,j),j=1,numEl)
				else
				read(40,*)iSeg,iPointX,pointX(iSeg,iPointX),pointY(iSeg,iPointX),		&
					(pointComp(iSeg,iPointX,j),j=1,numEl)
				endif
			end do
		end do

!	Now set the seg points that are next to the nodes
!	This is needed because the number of points in a segment changes with model progress
	do iSeg = 1,numSegs
		node1 = segNodes(iSeg,1)
		!The first point falls on the first node
		do j = 1,numNodeSegs(node1)
			if(nodeSegConnect(node1,j).eq.iSeg)then
				nodePointNextTo(node1,j) = numPointStart(iSeg)+1	! the second point in this segment (the first point is coincident with the node)
				nodePointOnTop(node1,j)  = numPointStart(iSeg)		! the firstpoint in this segment -- coincident with the node)
				endif
			end do
		node2 = segNodes(iSeg,2)
		do j = 1,numNodeSegs(node2)
			if(nodeSegConnect(node2,j).eq.iSeg)then
				nodePointNextTo(node2,j) = numPointEnd(iSeg)-1	! the second to last point in this segment (the last point is coincident with the node)
				nodePointOnTop(node2,j)  = numPointEnd(iSeg)		! the last point in this segment -- coincident with the node)
				endif
			end do
		end do	

!	Read Seg exclusion stuff
	read(40,*)dummy		! dashes
	do i = 1,numSegs
		!read(40,*)iSeg,segLength(iSeg),segGettingShorter(iSeg),segCaptureExclusion(iSeg)
		read(40,*)iSeg,segLength(iSeg),char1,char2
		if(char1.eq.'T')then
			segGettingShorter(iSeg) = .TRUE.
			else
			segGettingShorter(iSeg) = .FALSE.
			endif
		if(char2.eq.'T')then
			segCaptureExclusion(iSeg) = .TRUE.
			else
			segCaptureExclusion(iSeg) = .FALSE.
			endif
		end do

!	read crystal information
	read(40,*)dummy		! dashes
	do i = 1,numPhases
		read(40,*)dummy		! phase name
		do j = 1,numPhaseCrystals(i)
			read(40,*)PCN		! PCN = phase crystal number (DOESN'T CHANGE)
			if(PCN.ne.PhaseCrystalNumber(i,j))then
				call fss_alert('Alert','PCN does not equal PhaseCrystalNumber')				
				endif
			read(40,*)numCrystalNodes(PCN)
			do k = 1,numCrystalNodes(PCN)
				read(40,*)CrystalNodes(PCN,k),CrystalSegs(PCN,k)
				end do
			end do
		end do

	return
	end
			
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine AssignReactionPhases(iNode)
	implicit none
	include "Diffuse_GB_Hex.inc"
	integer*4 iNode,k,j,MIFID1,MIFID2,MIFID3,iSeg
	
	MIFID1 = nodeMIFID(iNode,1)
	MIFID2 = nodeMIFID(iNode,2)
	MIFID3 = nodeMIFID(iNode,3)
	if(MIFID1.eq.MIFID2.and.MIFID1.eq.MIFID3)then	! all phases are the same
		numNodeReactionPhases(iNode) = 0
		go to 20
		endif
	if(MIFID1.eq.0.or.MIFID2.eq.0.or.MIFID3.eq.0)then	! all phases are the same
		numNodeReactionPhases(iNode) = 0
		go to 20
		endif

	if(MIFID1.ne.MIFID2.and.MIFID1.ne.MIFID3.and.MIFID2.ne.MIFID3)then	! all 2 phases different
		numNodeReactionPhases(iNode) = 3
		nodeReactionPhases(iNode,1) = 1
		nodeReactionPhases(iNode,2) = 2
		nodeReactionPhases(iNode,3) = 3
		nodeSegTheSamePhases(iNode) = 0
		go to 20
		endif	
	! if here, then 2 phases must be the same
	numNodeReactionPhases(iNode) = 2
	if(MIFID1.eq.MIFID2)then
		nodeReactionPhases(iNode,1) = 1
		nodeReactionPhases(iNode,2) = 3
		do k = 1,3
			iSeg = nodeSegConnect(iNode,k)
			if(segMIFID(iSeg,1).eq.SegMIFID(iSeg,2))then
				nodeSegTheSamePhases(iNode) = k	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
				go to 20
				endif
			end do
		endif
	if(MIFID1.eq.MIFID3)then
		nodeReactionPhases(iNode,1) = 2
		nodeReactionPhases(iNode,2) = 3
		do k = 1,3
			iSeg = nodeSegConnect(iNode,k)
			if(segMIFID(iSeg,1).eq.SegMIFID(iSeg,2))then
				nodeSegTheSamePhases(iNode) = k	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
				go to 20
				endif
			end do
		endif
	if(MIFID2.eq.MIFID3)then
		nodeReactionPhases(iNode,1) = 1
		nodeReactionPhases(iNode,2) = 2
		do k = 1,3
			iSeg = nodeSegConnect(iNode,k)
			if(segMIFID(iSeg,1).eq.SegMIFID(iSeg,2))then
				nodeSegTheSamePhases(iNode) = k	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
				go to 20
				endif
			end do
		endif
	! if here, we have a problem
	write(12,*)' In Subroutine AssignReactionPhases(iNode)'
	write(12,*)' Failure to determine nodeReactionPhases and nodeSegTheSamePhases'
	write(12,*)' iNode = ',iNode
	write(12,227)iNode,(nodeMIFID(iNode,j),j=1,numCrystalsAtNode(iNode)),			&
			(nodeReactionPhases(iNode,j),j=1,numNodeReactionPhases(iNode)),nodeSegTheSamePhases(iNode)
227	format(20I8)

20	continue	

	return
	end	
	

