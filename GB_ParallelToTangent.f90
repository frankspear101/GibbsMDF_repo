
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine ParallelToTangent(k,output,izero)
!	routine to find composition of the phase where the tangent plane to the phase is parallel to the grain boundary tangent plane
!	Note that the grain boundary phase contains all the system components, so it is essentially the system tangent plane
!	This code is taken from Program Gibbs3 where the system tangent is defined by the equilibrium phase assemblage
!		so everywhere that it says "system tangent" think "grain boundary tangent"
!	That is find F = (µ(j)-µ(1)) - (µOnTan(j)-µOnTan(1)) = 0
!	where µ(j)-µ(1) is the tangent to the Gibbs surface for the phase and
!	µOnTan(j)-µOnTan(1) is the slope of the system tangent in the same composition direction
!
!	We want the Xj where they are equal -- this will be the composition of the phase that we want when the tangents are parallel
!
!	Using Newton's method we set up a system of equations of the form
!
!	Fj(Xo) + d(Fj(atXo))/dXj*∆Xj = 0
!
!	there is one equation for each independent composition variable (Xj).
!	All necessary values and derivatives are computed in subroutine dlnAdX (taken from dudTPX)
!	
!	After a call to dlnAdX, we have for each phase component (even dependent ones)
!	lnAct(jj) = (jj)/R*T		(for ideal solutions)
!	duTPX(jj,2+j) = d/dXj where the derivatives are only for each independent phase component (there are numPhCo terms and numPhCo-1 derivatives)
!		(the first two array elements in the array dudTPX are the T and P derivatives of each  but we don't need those here so we don't calculate them)
!
!	So our equations are calculated as
!	F(Xoj) = (uo(j)-uo(1)) + R*T*(lnAct(j) - lnAct(1)) - (µOnTan(j)-µOnTan(1))
!	d(FXoj)/dXj = dudTPX(j,2+j) - dudTPX(1,2+j)
!
!	so our matrix is
!	d(F2Xo)/dX2*∆X2 + d(F2Xo)/dX3*∆X3 + d(F2Xo)/dX4*∆X4 + ...    = -F2(Xo)
!	d(F3Xo)/dX2*∆X2 + d(F3Xo)/dX3*∆X3 + d(F3Xo)/dX4*∆X4 + ...    = -F3(Xo)
!	d(F4Xo)/dX2*∆X2 + d(F4Xo)/dX3*∆X3 + d(F4Xo)/dX4*∆X4 + ...    = -F4(Xo)
!
!	which we solve for ∆Xj . then compute
!	Xnewj = Xoj + ∆Xj
!	and repeat until we converge
!
!	The compositions of each phase at the minimum are stored in phCoatGMin(noPh,numPhCo)
!	The G of the phase is stored in gOfPhase(k)
!	The G on the tangent at this composition is stored in gOnTan(k)
!	The difference is stored in gDifference(k)
!	
	use MatrixArrays
	implicit none
! *****************************************
	include "Assemb.inc"
	include "GB_Tangent.inc"
	include "Newton.inc"
!	include "MIFarrays.inc"
! --------------------------------------------------------------------
	integer*4 i,k,ivol,j,jjj,izero,nsolve,ier,iflag,notconverged,jdep,j1,output,numeq
	integer*4 loopcount,looping,fixTheStepsize
	real*8 Xnew(phMax,phCoMax),deta,sum,stepsize,xTemp(20),LastYY(phCoMax),LastSlope(phCoMax,phCoMax),slope
	real*8 Rjoules
      	DATA Rjoules/8.3144D0/
! --------------------------------------------------------------------

	TK = TC+273.15
	Pfluid = Pb
	izero=0
	ivol=0            

	newtonStepMax = 0.1

!	If a phase displays a solvus, there are situations where there is no solution to the parallel tangent
!	This occurs because the equation to solve is a cubic and there are either 1 or 3 roots. 
!		When there is only 1 root there is no solution for the other phase
!		Two examples are muscovite-paragonite and alkali feldspars
!	We need to be able to identify these situations and stop trying to find a root
!	We will do this by examining the slope to the equations
!		F(Xoj) = (uo(j)-uo(1)) + R*T*(lnAct(j) - lnAct(1)) - (µOnTan(j)-µOnTan(1))
!		d(FXoj)/dXj = dudTPX(j,2+j) - dudTPX(1,2+j)
!	If the slope changes sign (eg negative to positive) and the function does not change sign, then
!		we don't have a root (remember we're trying to find the composition where the function = 0)
!		if the function is never zero we can't have a root.
!	I suppose there are situations where this algorithm won't work and we'll probably need to look at the matrix of second derivatives
!		to determine if it's positive or negative definite... etc.
!		For example, this might only work on a binary.
!	A more robust solution might be to find the minimum and maximum of the function and see if the function has different signs at these points.
!		That would ensure that there was a root.
!		

	! k is the phase we are working with

	if(output.ge.1)then		
		write(12,*)' '
		write(12,*)'+++++++++++++++'
		write(12,*)' In Sub ParallelToTangent'
		write(12,*)'Phase,PhaseType = ',k,phName(k),minRec(k),PhaseType(k)
		endif

	stepsize = 1.0		! this will be adjusted below so we stay in bounds
!	stepsize = 0.1		! this will be adjusted below so we stay in bounds
	fixTheStepsize = 0	! this will allow the stepsize to be set to 1 each new iteration
		

	call duTPX(1,ivol,izero)	! calculate the value of µzero for each component of the grain boundary
	call dlnAdX(1,izero)	! calculate the value of the activity of each component of the grain boundary
	if(izero.gt.0)then
		write(12,*)' Problem in dlnAdX (line 103 of ParallelToTangent())'
		return
		endif
	do j = 1,numPhCo(1)
		tanPlane(j) = gattp(1,j) + Rjoules*TK*lnAct(1,j)
		end do
				
	! Calculate the value of µzero for each end member in phase k based on the tangent plane
	call CalcMuOnTangent(k)		! this returns values for uOnTan(k,j) and uOnTanDelta(k,j)	

	if(output.eq.1)then
		write(12,*)'TanPlane    '
		do j = 1,numPhCo(1)
			write(12,89)j,phCoName(1,j),xPhCo(1,j),lnAct(1,j),tanPlane(j)
			end do
		write(12,*)' uOnTanAU'
		do j = 1,numPhCo(k)
			write(12,89)j,phCoName(k,j),xPhCo(k,j),lnAct(k,j),uOnTanAU(j)
			end do
		write(12,*)' '
		end if

	call duTPX(k,ivol,izero)	! calculate the value of µzero for each phase component at T and P of interest
	call dlnAdX(k,izero)	! calculate the value of the activity of each component of the phase of interest
	if(izero.gt.0)then
		write(12,*)' Problem in dlnAdX (line 128 of ParallelToTangent())'
		return
		endif


!	For a phase of fixed composition we need only to calculate the G on the tangent
!		at the same composition as the phase(k)
!	This is already done in Subroutine CalcMuOnTangent
!		Result is stored in uOnTan(j)
	if(numPhCo(K).eq.1)then
		!call FSS_Alert('Alert','We should not be here')
		gOnTanAU(K)   = uOnTanAU(1)
		gOfPhaseAU(K) = gattp(k,1)/numOxygens(k)
		gDifferenceAU(K) = gOfPhaseAU(K) - gOnTanAU(K)
		! Note that this should be exactly the same value as we get in Sub CalcAffinity
		! We'll use this as a check for now
		go to 999				! done with this phase
		endif

!	This is code for phases with more than 1 phase component
	jdep = 1
	numEq = numPhCo(K) - 1
	loopcount = 0
	looping = 0
	do j = 1,numEq
		LastYY(j) = 0.0d0
			do i = 1,numEq
			LastSlope(j,i) = 0.0d0
			end do
		end do

	do j = 1,numPhCo(k)
		uzeroPhaseDeltaAU(j) = (gattp(k,j) - gattp(k,1))/numOxygens(k)
		uOnTanDeltaAU(j)   = uOnTanAU(j) - uOnTanAU(1)	
		uPhaseAU(j) = (gattp(k,j) + Rjoules*TK*lnAct(k,j))/numOxygens(k)
		end do
	do j = 1,numPhCo(k)
		uPhaseDeltaAU(j) = uPhaseAU(j) - uPhaseAU(1)
! 		uzeroPh_TanDeltaAU(jjj) = uzeroPhaseDeltaAU(j) - uzeroTanDeltaAU(j)
		end do
	if(output.eq.2)then
		write(12,*)' '
		write(12,*)' Values for this problem '
		write(12,*)'TanPlane    ',(tanPlane(j),j=1,numPhCo(1))
		write(12,*)'  Chemical Potentials (1 oxygen basis)'		
		write(12,*)'  We want uPhaseDelta to equal uOnTanDelta (i.e.difference = 0)'		
		write(12,*)'                uPhase/ox        uPhaseDelta         uOnTanDelta       Diff'		
		do j = 1,numPhCo(K)
			write(12,89)j,phCoName(k,j),uPhaseAU(j),uPhaseDeltaAU(j),uOnTanDeltaAU(j),uPhaseDeltaAU(j)-uOnTanDeltaAU(j)
 	89		format(I4,1x,A10,8F15.5)
			end do

! 		write(12,*)' Deltas '
! 		write(12,*)'                uZeroDeltaAU     uOnTanDeltaAU      Difference'		
! 		do j = 1,numPhCo(k)
! 			write(12,89)j,phCoName(k,j),uzeroDeltaAU(j),uOnTanDeltaAU(j),uPh_TanDeltaAU(j)
! 			end do
		endif
			
!	loop on Newton's method begins here
4000	continue
	loopcount = loopcount+1
	izero = 0
	! note that this routine doesn't work on MIF arrays
	call dLnAdX(K,izero)				! update the activities based on the current phase composition
	if(izero.gt.0.and.output.ge.1)then
      		write(12,*)' error in Subroutine dLnAdX called from subroutine ParallelToTangent in file TangentRoutines.f'
		write(12,*)' Loop count - ',loopcount,TC,PB
		write(12,*)' Phase, comp',K,phName(k),(xPhCo(k,j),j=1,numPhCo(k))
		write(12,*)' Hit return to continue'
      		pause ' Hit return to continue - Line 198 in Sub ParallelToTangent'
      		endif
!	So our equations are calculated as
!	F(Xoj) = (uo(j)-uo(1)) + R*T*(lnAct(j) - lnAct(1)) - (µOnTan(j)-µOnTan(1))
!	d(FXoj)/dXj = dudTPX(j,2+j) - dudTPX(1,2+j)

	! calculate 
	do j = 1,numEq   ! numEq = numPhCo(K) - 1 (calculated above)
		jjj = j + 1	! we skip the first phase component
! 		YY(j) = -(uzeroDeltaAU(k,jjj) + Rjoules*TK*((lnAct(k,jjj)/atomNorm(k,jjj)) - (lnAct(k,jdep)/atomNorm(k,jdep))) &
! 	     	              - uOnTanDeltaAU(k,jjj))
		!YY(j) = -(uzeroDeltaAU(k,jjj) + Rjoules*TK*((lnAct(k,jjj)/numOxygens(k)) - (lnAct(k,jdep)/numOxygens(k)))  - uOnTanDeltaAU(k,jjj))

! 		YY(j) = -(uzeroPh_TanDeltaAU(jjj) + Rjoules*TK*(lnAct(k,jjj) - lnAct(k,jdep))/numOxygens(k))

		YY(j) = -(uzeroPhaseDeltaAU(jjj) - uOnTanDeltaAU(jjj) + Rjoules*TK*(lnAct(k,jjj) - lnAct(k,jdep))/numOxygens(k))
	
		! remember that uzeroDelta(K,j) are differenes - the first j=1 is a dummy
		A(j,numEq+1) = YY(j)
		do i = 1,numEq
			slope = (dudTPX(k,jjj,2+i) - dudTPX(k,jdep,2+i))/numOxygens(k)
			LastSlope(j,i) = slope
			AA(j,i) = slope
			A(j,i) = AA(j,i)	
			end do
		LastYY(j) = YY(j)
		end do

	if(output.eq.2)then		
		write(12,*)' '
		write(12,*)'Y = ',(YY(j),j=1,numEq)
		write(12,*)'X = ',(xPhCo(k,j),j=1,numPhCo(k))
		write(12,*)'A matrix'
		do i = 1,numEq
			write(12,80)(A(i,j),j=1,numEq+1)
	80		format(50F15.5)
			end do
		endif		

!     Find solution
	nsolve = numEq + 1
	DETA=0.D0
	IER=0
	CALL REDUCE (numEq,NSOLVE,DETA,IER)
	IF (IER.EQ.1) then
	      WRITE(12,*)' ************ ERROR **************************'
	      WRITE(12,*)' Matrix failed to invert in SUBROUTINE REDUCE'
	      write(12,*)' We are in Subroutine ParallalToTangent in file line 180 or so'
	      write(12,*)' Phase = ',K,phName(k)
	      izero=1
!	      write(12,*) 'Hit return to continue...'
!	      pause 'Hit return to continue...'
	      return
	      endif

	if(output.eq.2)then		
		write(12,*)'Solution xx'
		write(12,84)(xx(j,1),j=1,numEq)
	84	format(50E15.5)
		endif		

!	Check to see if the stepsize is too large.
!	In Newton's method, the solution finds the roots of the equation. 
!		If the curvature of the function is too great, this can result in the calculated composition being
!			outside of the range 0-1
!		Check this and reduce stepsize until we stay within range.
!		The same code is used in Subroutine C2ompute

!	make initial calculation with original stepSize
5001	continue
	sum = 1.0d0
	do j = 1,numEq
		jjj = jdep+j		! this is component number 2, 3, 4, etc
		j1 = j+1
		xTemp(j1) = xPhCo(k,jjj) + stepsize*xx(j,1)			
		sum = sum - xTemp(j1)
		end do
	xTemp(1) = sum	!dependent composition variable (number 1)
	if(output.eq.2)then
		write(12,*)'stepsize = ',stepsize		
		write(12,*)'   PhaseCo         Xcurrent         Xnew(tempP'		
		do j = 1,numPhCo(K)
			write(12,88)phCoName(k,j),xPhCo(k,j),xTemp(j)
	88		format(A10,5x,2F12.5)
			end do
		endif

	iflag = 0
	call CheckStoichiometry(k,xTemp,iFlag)
	if(iflag.eq.1)then			! failed stoichiometry test
		stepSize = stepSize/2.0d0
		if(stepSize.gt.1.0d-10)go to 5001			! try again with this smaller stepsize
!		stepsize is too small. Abort this phase
		write(12,*)'stepSize = ',stepSize
		write(12,*)'This is too small. Something must be wrong in the code SUB ParallelToTangent - line 224'
		write(12,*)'Phase = ',k
		write(12,*)phName(k)
		do j1 = 1,numPhCo(K)
			write(12,*)phCoName(k,j1),xPhCo(k,j1),xTemp(j1)
			end do
		write(12,*)' '
		write(12,*)'Y = ',(YY(j),j=1,numEq),(xPhCo(k,j),j=1,numPhCo(k))
		write(12,*)'A matrix'
		do i = 1,numEq
			write(12,80)(A(i,j),j=1,numEq+1)
			end do
		write(12,*)'Solution xx'
		write(12,80)(xx(j1,1),j1=1,numEq)
		write(12,*) ' hit return to continue with the next phase in the list'
	!	pause ' hit return to continue with the next phase in the list'
		endif
!	if here, then iflag = 0 and we are OK
! 	calculate new compositions with this stepSize
	sum = 1.0d0
	do j = 1,numEq		!number of independent components
		jjj = jdep+j		! this is component number 2, 3, 4, etc
		xnew(k,jjj) = xPhCo(k,jjj) + stepsize*xx(j,1)
		sum = sum - xnew(k,jjj)
		end do
	xnew(k,jdep) = sum	!dependent composition variable (number 1)
	notConverged = 0
	do j = 1,numPhCo(K)
		jjj = j
		if(Dabs(xnew(k,jjj)-xPhCo(k,jjj)).gt.0.00001)then	! TolNewton is the fraction of X to scale convergence criterion
			notConverged = 1
			endif
		end do


	if(output.ge.1)then
		if(notConverged.eq.0)then			! notConverged = 0 means we have converged
			write(12,*)'Converged     '
			do j = 1,numPhCo(k)
				write(12,83)phCoName(k,j),xPhCo(k,j),xnew(k,j),xnew(k,j)-xPhCo(k,j)
	83			format(A8,T10,20F15.5)
				end do
			else
			write(12,*)'NOT converged '
			do j = 1,numPhCo(k)
				write(12,83)phCoName(k,j),xPhCo(k,j),xnew(k,j),xnew(k,j)-xPhCo(k,j)
				end do
			endif	
		write(*,*) ' Hit 0 to continue, 1 to abort'
		read(*,*)j
		if(j.eq.1)return
		endif

	if(notConverged.eq.1)then 		! we have not converged
		if(loopCount.gt.1000)then	! something is wrong - we should have converged but we havent
						! this can happen when paragonite bounces to muscovite and back (across the solvus)
			write(12,*)'         LoopCount = 1000  Stepsize,Phase = ',stepSize,k
			looping = looping + 1
			if(looping.gt.3)then
				write(12,*)'         Looping = 3 -- abort this calculation ',k
				gDifferenceAU(K) = 1.0d4		! Make = 10000. just to move on
				go to 999
				endif
			do j = 1,numPhCo(K)
	!			xPhCo(k,j) = xPhCoLast(k,j)	! start over with last good compositions (this will probably oscillate)
				xPhCo(k,j) = xPhCoInitial(k,j)	! start over with initial compositions
				end do
			loopCount = 0			! reset loopcount		
			if(looping.eq.1)then
				stepsize = 0.01			! make stepsize small this time
				else
				stepsize = .001
				endif
			else				! Loopcount < 100
			do j = 1,numPhCo(K)
				xPhCo(k,j) = xnew(k,j)
				end do
			endif
		if(fixTheStepSize.eq.0)then		! only make stepsize large if we are not looping
			stepsize = 1.0
!			stepsize = 0.1
			endif
		go to 4000
		endif

	do j = 1,numPhCo(K)
		xPhCo(k,j)    = xnew(k,j)
		end do

!	Converged - now calculate activities one more time to get gPhase
	izero = 0
	call dLnAdX(K,izero)
	! this call should calculate the activities for phase K provided it was loaded as an assemblage
	if(izero.gt.0.and.output.ge.1)then
	      	write(12,*)' error in Subroutine dLnAdX called from subroutine ParallelToTangent in file TangentRoutines.f'
		write(12,*)' Loop count - ',loopcount
      		write(12,*) ' This is second call - Hit return to continue'
      		pause ' This is second call - Hit return to continue'
      		return
      		endif
	if(izero.gt.0)then	! write to the log file
		write(75,*)' Error in ParallelToTangent -> Call dLnAdX after convergence (this should never be a problem) '
		write(75,*)' Phase, comp',K,phName(k),(xPhCo(k,j),j=1,numPhCo(k))
		endif
		
!	Calculate the differences in µ for the phase components and see if they are correct
	do j = 1,numPhCo(k)
	!	uzeroPhaseDeltaAU(j) = (gattp(k,j) - gattp(k,1))/numOxygens(k)  -- already done above and shouldn't change
	!	uOnTanDeltaAU(j)   = uOnTanAU(j) - uOnTanAU(1)	--- shouldn't have changed
		uPhaseAU(j) = (gattp(k,j) + Rjoules*TK*lnAct(k,j))/numOxygens(k)	! should be the new value
		end do
	do j = 1,numPhCo(k)
		uPhaseDeltaAU(j) = uPhaseAU(j) - uPhaseAU(1)		! should be the new value
		end do

	if(output.eq.2)then
		write(12,*)' '
		write(12,*)' Final values for this problem '
		write(12,*)'  Chemical Potentials (1 oxygen basis)'		
		write(12,*)'  We want uPhaseDelta to equal uOnTanDelta (i.e. Difference = 0)'		
		write(12,*)'                uPhase/ox        uPhaseDelta         uOnTanDelta      Diff'		
		do j = 1,numPhCo(K)
			write(12,89)j,phCoName(k,j),uPhaseAU(j),uPhaseDeltaAU(j),uOnTanDeltaAU(j),uPhaseDeltaAU(j)-uOnTanDeltaAU(j)
 	!89		format(I4,1x,A10,3F15.5)
			end do
		endif
	if(output.ge.1)then
		write(12,*)' '
		write(12,*)' Affinities '
		write(12,*)'                uPhase-uOnTangent (should be the same for each phase component'		
		do j = 1,numPhCo(k)
			write(12,89)j,phCoName(k,j),uPhaseAU(j)-uOnTanAU(j)
			end do

		endif




	! Calculate G of the phase at this composition and compare with G on the tangent
	! These were used in program Gibbs3 for MAD calculations
	! I don't use them here -- maybe just omit?
	gOnTanAU(K) = 0.0d0
	gOfPhaseAU(K) = 0.0D0
	do j = 1,numPhCo(K)
		gOnTanAU(K)   = gOnTanAU(K)   + xPhCo(k,j)*uOnTanAU(j)
		gOfPhaseAU(K) = gOfPhaseAU(K) + xPhCo(k,j)*(gattp(k,j) + Rjoules*TK*lnAct(k,j))/numOxygens(k)
		gDifferenceAU(K) = (gOfPhaseAU(K) - gOnTanAU(K))
		end do

!	End routine for this phase
999	continue
	if(output.ge.1)then
		write(12,*)' '
		write(12,*)' Gsystem stuff '
		write(12,*)'                 gOfPhaseAU     gOnTanAu       gDifferenceAU      Phase Comp '
		write(12,81)phName(K),gOfPhaseAU(K),gOnTanAU(K),gDifferenceAU(k),(xPhCo(k,j),j=1,numPhCo(K))
81		format(A10,T15,3F15.2,15F12.5)
		write(12,*)' '
		!write(12,*)' '
		write(12,*)' Leaving Sub ParallelToTangent'
		write(12,*)'+++++++++++++++'
		write(12,*)' '

		endif

!	newtonStepMax = 1.
	newtonStepMax = .4

	return
	end
	



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine CalcMuOnTangent(k)
	implicit none
! *****************************************
	include "Assemb.inc"
	include "GB_Tangent.inc"
!	include "MIFarrays.inc"
	integer*4 j,l,k,i
	real*8 sum
!----------------------------------------------------
! 	Given the current tangent plane, 
!	calculate the value of µ for each phase component on the tangent
!	These calculations are all done in atom units
	do j = 1,numPhCo(k)
		sum = 0.0d0
		do l = 1,nc
			!sum = sum + tanPlane(l)*comp(k,j,l)/atomNorm(k,j)		! comp() is in molar units
			!sum = sum + tanPlane(l)*comp(k,j,l)		! comp() is in molar units
			i = GB_to_Sys_Co(L)
			sum = sum + tanPlane(l)*comp(k,j,i)		! comp() is in molar units
			end do
		uOnTanAU(j) = sum/numOxygens(k)
		! This is the chemical potential of the END MEMBER component on the GB tangent plane
		end do

! 	do j = 1,numPhCo(k)
! 		uOnTanDeltaAU(j) = uOnTanAU(j) - uOnTanAU(1)
! 		! uOnTanDelta = (independent phase component - dependent phase component)
! 		end do
	return
	end




	
	
		
