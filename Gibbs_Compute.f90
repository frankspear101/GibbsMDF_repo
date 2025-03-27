
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Compute4EQUIL(izero)
	use MatrixArrays
      implicit none

! 	Routine to calculate Jacobian for equilibrium model phases
!	This routine is a generalized solution where
!	(a) Reactions are stored in array ARX
!		EQ phases have an entry for every phase component
!		OS phases have an entry for nphco-1 "components" and they are all PhCo(ind)-PhCo(dep)
!		Note that there are npc - 1 such independent equations
!	Basically, the method should provide a matrix solution that allows any number of phases to be EQ or OS
!
!	The variance of this system should be 2 + #OS phases so we can solve for changes in all phases given the growth of the OS phase
!		For example dMph/dMOSphase
!	Original version written June 8, 2018 by F. Spear
!	Generalized version written Oct 31, 2018 by F. Spear
! 	Note that this is the same code as in C2OMPUTE except 
!		generalized for OS phases
!		No contouring option
!		No generic "new" variables
! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "Singular.inc"
!	include "Solute.inc"
!	include "Tangent.inc"
	include "MatrixSingular.inc"
! c*****************************************
!     local variables
      integer*4 i,ii,j,jj,k,L,LL,Lcont,jcol,irow,ier,j1,kCur,			&
     	 izero,imiss,nsolve,newtoncounter,ivol
      REAL*8 DETA,R,SUM,Pi,TEMPT,TEMPP,TEMPX,TEMPH,TEMPS,TEMPK,xTemp(20),allXTemp(100)
      CHARACTER vn51(100)*8
      DATA R,Pi/8.3144D0,3.141592654D0/


! --------------------------------------------

      IF(iLong(2).eq.1)then

         write(12,*)' '
         write(12,*)'------------------------------------------------------'
         write(12,*)'------------------------------------------------------'
         write(12,*)' Variable list '
         do i = 1,nvar
        	write(12,15)I,VN1(I),VN2(I)
15      	FORMAT(I3,2X,A2,A16)
		end do
         write(12,*)'numEqEquns=',numEqEquns,'  NVAR=',NVAR
         write(12,*)
         write(12,*)'VARIANCE =',NVAR-numEqEquns
         write(12,*)
         write(12,*)'MONITOR PARAMETERS    Delta X for monitors'
         if(nvar-numEqEquns.gt.0)then
         do i = 1,nvar-numEqEquns
         	write(12,16)MON(I),VN1(MON(I)),VN2(MON(I)),Deltax(i)*NSTEP
16       	FORMAT(I3,2X,A2,A16,F15.8)
		end do
	!write(12,*)' DELTA Xi FOR MONITOR PHASES'
	!write(12,16)(DELTAX(I)*NSTEP,I=1,NVAR-numEqEquns)
	!16       FORMAT(10F12.3)
         else
         write(12,*)' No monitor parameters (variance <= 0)'
         endif
         write(12,*)'NUMBER OF FINITE DIFF. ITERATIONS=',NSTEP
         write(12,*)
         write(12,*)'INITIAL MINERAL DATA'
         write(12,17)TC,PB
17       format(' T START (DEG C)=',F10.2,'   P START=',F10.2)
	if(PFluidSwitch.gt.0)write(12,18)PFluid
18	format(' PFluid = ',F10.2)
	endif
90    CONTINUE
!     LCONT IS COUNTER IN CONTOUR BLOCK
      LCONT=0
      newtoncounter=0

! -----------------------------------------------------------
!                             
100   CONTINUE


      TK=TC+273.15D0

!     THESE CALLS CALCULATE THE THERMODYNAMIC
!     PROPERTIES OF THE PHASES OF INTEREST

!     SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      	izero=0
      	ivol=0            !compute everything
      	K=0
	! x       write(*,*)' In c2ompute #3-calling duTPX'
	! x       pause
	!	This call will callculate all thermodynamic properties - the TP dependent ones and the composition dependent ones
	!	We need this when we're using T or P as dependent variables. But not when they're independent (because they're already set)
      	CALL ALLKduTPX(ivol,izero)
	! x       write(*,*)' In c2ompute #4-back from duTPX'
	! x        pause
	if(izero.gt.0)go to 9999

	!COMPUTE PHASE VOLUMES
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
            	VP0(K)=MP0(K)*VMOL(K)*10.D0
		end do

	!OUTPUT THERMO DATA TO DISK

      IF(iLong(3).eq.1)then
350   	CONTINUE
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
          	write(12,*)
          	write(12,200)MINREC(K),SITMUL(K),PHNAME(K)
200       	FORMAT(' ',I4,F8.1,4X,A32)
          	write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
          	IF(IMASS.EQ.1) then
          		write(12,507)MP0(K),VP0(K)
          		endif
507       	FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
          	write(12,501)(phCoName(k,j),J=1,numPhCo(K))
501       	FORMAT('          ',12(A4,6X))
          	write(12,504)(xPhCo(k,j),J=1,numPhCo(K))
504       	FORMAT('  COMP ',12(F10.4))
          	write(12,505)(Dexp(lnAct(k,j)),J=1,numPhCo(K))
505       	FORMAT(' Activ  ',12(E10.3))
          	write(12,506)(lnAct(k,j),J=1,numPhCo(K))
506       	FORMAT('  Ln(a) ',12(F10.5))
          	write(12,514)(hPhCoZero(k,j),J=1,numPhCo(K))
514       	Format(' hPhCoZero',12E15.7)
          	write(12,515)(HATTP(k,j),J=1,numPhCo(K))
515       	Format(' HATTP',12E15.7)
          	write(12,509)(sPhCoZero(k,j),J=1,numPhCo(K))
509       	Format(' SZERO',12F15.5)
          	write(12,510)(SATTP(k,j),J=1,numPhCo(K))
510       	Format(' SATTP',12F15.5)
          	write(12,508)(vPhCoZero(k,j),J=1,numPhCo(K))
508       	Format(' VZERO',12F15.5)
          	write(12,513)(VATTP(k,j),J=1,numPhCo(K))
513       	Format(' VATTP',12F15.5)
          	write(12,516)(GATTP(k,j),J=1,numPhCo(K))
516       	Format(' GATTP',12E15.7)
          	write(12,*)'dudTPX ARRAY:'
          	write(12,*)'     dudT   .    dudP   .    dudX2  .    dudX3  .    dudX4...'
           	do j = 1,numPhCo(K)
           		write(12,511)dudTPX(k,j,1),dudTPX(k,j,2),(dudTPX(k,j,L),L=3,numPhCo(K)+1)
511       		FORMAT(15E12.5)
			end do
408       	format(' ',a5,6E10.3)
		end do

	    endif


!     THIS PART OF THE PROGRAM SETS UP THE MASTER MATRIX

!     ZERO MASTER ARRAY
      	DO J=1,NVAR
		DO I=1,numEqEquns
			AA(I,J)=0.0D0
			end do
		end do
!     The partial derivatives in matrix dudTPX are:
!     In storage space 1: dudTPX(jj,1):
!     Temperature:      /T = - S(P,T,X)  -- Partial molar Entropy at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state entropy, plus the entropy of mixing plus any exess entropy
!                                              terms.
!     In storage space 2: dudTPX(jj,2):
!     Pressure   :      /P =   V(P,T,X)  -- Partial molar Volume at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state volume, plus any exess volume
!                                              terms.
!     In storage space 3 to (numPhCo(K)+1): dudTPX(jj,3), dudTPX(jj,3), dudTPX(jj,3) etc:
!           Note that there are only numPhCo(K)-1 partial derivatives and the subscripts range
!           from 3 to numPhCo(K)+1.
!     Composition:      /X = (ideal)/X + (margules)/X + (reciprocal)/X

!       
!     The matrix Arx(numEqrxn,npc) contains the stoichiometric coefficients for each phase component in the EQ reactions
!	Note that these equations define the phase components for the minerals in terms of the phase components of the grain boundary (matrix)
!     in each of the reactions, subscript 1 refers to reaction number, subscript 2 refers to the phase component

!     The matrix AA therefore contains in each cell terms such as
!     
!       (n(i)*u(i))/T
!       (n(i)*u(i))/P
!       (n(i)*u(i))/X2
!       (n(i)*u(i))/X3

!     for each independent variable (where n(j,i) are the stoichiometric coefficients for phase component i in reaction j)
!       Each row of AA contains coefficients for each of the reactions.

!     The only catch to the calculation of AA is that the matrix dudPTX is
!     not stored in true matrix form.  That is, to save space the independent variables are listed for each phase only.
!     This simply requires a bit of extra bookkeeping during the matrix multiplication

	ii = np - numEqRxn
	do i = 1,numEqRxn		! loop for each reaction
		TempT = 0.0d0
		TempP = 0.0d0
		j = 0
		do kcur = 1,numPh		! Loop for every phase
			k = asmCurrent(kCur)
			do jj = 1,numPhCo(k)	! Loop for every phase component in every phase
				j = j+1
				! dT and dP coefficient
				TempT = TempT + Arx(i+ii,j)*(dudTPX(k,jj,1))	! index 1 is partial molar S
				TempP = TempP + Arx(i+ii,j)*(dudTPX(k,jj,2))	! index 2 is partial molar V
				end do
!1015				continue
			end do
!1014			continue
	      	AA(i,1)=TEMPT
	      	AA(i,2)=TEMPP
		end do
!1010		continue

!	Do dX terms
	do i = 1,numEqRxn
		L = 2				! starts in the third column (2 = dT, dP)
		jj = 0
		do kcur = 1,numPh
			k = asmCurrent(kCur)
			if(numPhCo(K).eq.1)then
				jj = jj + 1
				go to 1021	! there are no dX terms for a phase of fixed composition
				endif
			Do LL = 3,numPhCo(K)+1     	! this loops on the number of independent composition derivatives
				tempX = 0.D0
				do j=1,numPhCo(K)          	! this loops on the number of phase components in each K phase
					TEMPX=TEMPX + Arx(i+ii,jj+j)*dudTPX(k,j,LL)
					end do
	!				1023  	continue
				L=L+1
				AA(i,L)=TEMPX
				end do
!				1022  	continue
			jj = jj + numPhCo(k)
1021  			continue
			end do
		end do
!1020  		continue




!     	SET UP MASS BALANCE EQUATIONS
      	if(imass.eq.0)go to 1100
!     	FIRST SET UP dM Columns (M=MOLE FRACTION OF PHASE)
!     	JCOL COUNTS WHERE dM TERMS ARE
      	JCOL = TANDP + NX
!     	loop through each phase (there is a dM term for each phase)
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
	      	JCOL = JCOL + 1
		!LOOP THROUGH ALL mass balance equations for each SYSTEM COMPONENT
      		call bulkcomp				! get total moles of each system component (moles(i) in common block)
      		DO i = 1,nc			! loop through every mass balance equation
      			iRow = i + numEqRxn			! iRow counts the row number for the mass balance equations
      			sum=0.D0
			!LOOP THROUGH ALL PHASE COMPONENTS for this phase
      			DO j=1,numPhCo(K)
				sum = sum + comp(k,J,I)*xPhCo(k,J)
				end do
				!1120  	continue
	
		      	AA(IROW,JCOL)=SUM
		      	YY(IROW) = -(moles(i) - molesStart(i))			! molesStart(i) is the bulk composition we want to fit
      										! moles(i) is the current bulk comp based on phase comp and M(k)
										! when they are equal, we have the correct X and M for the desired bulk comp
										! openFlux are open system components either dependent or independent variables
			end do
			!1105  	CONTINUE
		end do
		!1110  	CONTINUE	!loop for each phase

!     	NOW SET UP COEFFICIENTS FOR dX TERMS in mass balance equations
!     	Note that the dependent dX term is removed as above
      	JCOL=TANDP					!     	JCOL is the column for the dX term
!     	FIRST LOOP THROUGH ALL PHASES TO FIND independent dX TERMS
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
		IF(numPhCo(K).EQ.1) GO TO 1130
		!IF HERE THEN WE HAVE FOUND A PHASE WITH dX TERMS (I.E. MORE THAN 1 COMPONENT)
		!LOOP THROUGH ALL dX TERMS FOR THIS PHASE
			DO J=2,numPhCo(K)
				JCOL=JCOL+1
				!LOOP THROUGH ALL MASS BALANCE EQUATIONS FOR EACH dX TERM
				DO I=1,NC
					!IROW= I + NRX
					IROW= I + numEqRxn
					AA(IROW,JCOL)=MP0(K)*(COMP(k,J,I)-COMP(k,1,I))
					end do
					!1150	continue
				end do
				!1140  	CONTINUE
1130  		CONTINUE
		end do

1100  continue   ! end of mass balance equations


!     If we are in MODE=NEWTON then we need to set up a data vector that
!     contains -(H - TS + R T lnKeq) for each reaction
!     (note there is no V term here because H and S are calculated at T and P
      	if(inewton.eq.1) then
!		ii = np-nrx      ! note that reaction coefficients are stored in the end of array ARX
		DO i=1,numEqRxn
			TEMPH=0.D0
			TEMPS=0.D0
			TEMPK=0.D0
			jj = 0
			Do kCur = 1,numPh
				k = asmCurrent(kCur)
				do J = 1,numPhCo(k)
					jj = jj + 1
					TEMPH =TEMPH + Arx(i+ii,jj)*HatTP(k,j)		! ∆H reaction
					TEMPS =TEMPS + Arx(i+ii,jj)*SatTP(k,j)		! ∆H reaction
					TEMPK =TEMPK + Arx(i+ii,jj)*lnAct(k,j)		! ∆H reaction
					end do
				end do
				!1415    	continue      
 			YY(i)=-(TEMPH - TK * TEMPS + R*TK*TEMPK)
			end do
			!1410     	continue
		endif      

!     	Rearrange MASTER matrix AA
!     	MOVE MONITOR PARAMETERS TO RIGHT SIDE OF EQUNS IN ARRAY A
	nsolve=nvar
	if(inewton.eq.1)nsolve=nvar+1
	DO I=1,NVAR-numEqEquns
		DO J=1,numEqEquns
		  	A(J,numEqEquns+I)=-AA(J,MON(I))
			end do
		end do
!2100	continue
!     	FILL LEFT SIDE OF EQUATIONS IN MATRIX A FROM AA
      	DO I=1,numEqEquns
 	     	DO J=1,numEqEquns
			A(J,I)=AA(J,IPOINT(I))
			end do
		end do
!2110	continue
      	if(inewton.eq.1)then
!           Fill Nvar+1=nsolve array element in A with the Y vector      
            DO J=1,numEqEquns
        	A(J,Nsolve)=YY(J)
		end do	
!2115        A(J,Nsolve)=YY(J)
            endif
            

!     	Output matrix and modified matrix, if desired
	if(iLong(5).eq.1)then
2330  		CONTINUE
      		write(12,*)' '
      		write(12,*)' '
      		write(12,*)' MASTER MATRIX'
	        write(12,2321)(VN1(j),VN2(j),J=1,nvar),'  YY(i)'
2321  		FORMAT('     ',30(A2,A4,6x))
	      	DO I=1,numEqEquns
			!write(12,*)'--------------------'
			write(12,2320)(AA(I,J),J=1,NVAR),YY(I)
			end do
			!2310  		continue
2320  		FORMAT(' ',30E12.4)
		do i = 1,nvar-numEqEquns
			vn51(numEqEquns+i) = trim(vn1(mon(i)))//vn2(mon(i))
			end do
			!2311		continue
		do i = 1,numEqEquns
			vn51(i) = trim(vn1(ipoint(i)))//vn2(ipoint(i))
			end do
			!2312		continue
		write(12,*)' '
		write(12,*)' '
		write(12,*)' MODIFIED MASTER MATRIX'
	        write(12,2313)(VN51(j),J=1,nvar),' YY(i)'
2313  		FORMAT('     ',30(A6,6x))
		DO I=1,numEqEquns
			!write(12,*)'--------------------'
			write(12,2320)(A(I,J),J=1,nsolve)
			end do
			!2315  		continue
		endif

!     Find solution
      	DETA=0.D0
      	IER=0
      	CALL REDUCE (numEqEquns,NSOLVE,DETA,IER)
      	IF (IER.EQ.1) then
		izero=1
!		if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)go to 2509
		WRITE(*,*)' ************ ERROR **************************'
		WRITE(*,*)' Matrix failed to invert in SUBROUTINE REDUCE'
		WRITE(*,*)' The most likely cause is that the variables'
		WRITE(*,*)' specified as monitor parameters are not linearly'
		WRITE(*,*)' independent.  Note that only the number of'
		WRITE(*,*)' intensive variables equal to the phase rule'
		WRITE(*,*)' variance of the assemblage are independent.'
		WRITE(*,*)' Return to previous menu and try again.'
		write(*,*)' '
		pause 'Hit return to continue...'
		go to 9999
		endif

2509	continue


	if(matrixIsSingular.eq.1)go to 9999

	if(iLong(2).eq.1)then
      		write(12,*)'   '
      		write(12,*)' RESULTS OF THIS FINITE DIFFERENCE ITERATION'
      		write(12,*)'DETERMINANT=',DETA
		endif


! -------------------------------
!     	Set ALLX(5,L) equal to current values of ALLX(1,L)
!     	Note that below ALLX(1,L) is incremented by DELX
      	DO L=1,NVAR
		ALLX(5,L)=ALLX(1,L)
		end do
		!3000  	ALLX(5,L)=ALLX(1,L)

! ____________________Newtons method___________________________
      	if(inewton.eq.1)then
!     		XX IS A numEqEqunsx(NVAR-numEqEquns+1) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		The last column (NVAR-numEqEquns+1) is the solution to the x for Newtons method
!     		Calculate the change in the variables
!     		New variables are computed as
!     		Xnew = Xold + X
! 		Note that we are doing these calculations with the values of the "independent" variables constant
! 		i.e. dxi = 0. For example, if it is univariant and we are doing calculations at constant P then dP = 0
! 		This way, we only need to solve for the values of the "dependent" variables (dT, dXi, etc.)
!

		! skip this newton step stuff -- it seems to be causing trouble
!		go to 5010

!       	Composition variables can only range from 0-1. 
!		If the solution vector component (xx) is too large, it may result in the change in mole fraction (xPhCo(i) exceeding this range
!		This block of code checks this and adjusts newtonStep to ensure this doesn't happen
!		newtonStep = 1.0d0

		newtonStep = newtonStepMax
		jCol = NVAR-numEqEquns+1	! index of solution vector in array XX
		do i = 1,nvar
			allXtemp(i) = allX(1,i)	! load up allXTemp
			end do
			!5002		continue
! ----		loop to here for newtonStep iteration --------
5001		continue
		do i = 1,numEqEquns
			ii=ipoint(i)
			!adjust the values of AllXtemp that we have just calculated
			!note that this changes all dependent values, not just the composition ones
			!we could make this a wee bit faster if we just calculated the composition ones but we're not sure which of them are dependent
			!jcol is the column where the solution vector resides in XX
			ALLXtemp(ii)=ALLX(1,ii) + newtonstep*XX(i,jCol)
			end do
			!5004		continue
!		now pick out the composition values and store in xTemp
		l = TandP		! this is where the independent composition derivatives start
		Do kCur = 1,numPh
			k = asmCurrent(kCur)
			if(numPhCo(k).eq.1)go to 5010	! no composition derivatives for 1 component phases
			sum = 1.0d0		! this is for the dependent composition derivative (the first one)
			DO J = 2,numPhCo(K) ! loop from second to last phase component for this mineral
				l = l + 1		! index of independent composition variable in AllX
				!j indexes the phase components of phase K
				xTemp(j) = allXTemp(l)
				sum = sum - xTemp(j)
				end do
				!5020  		continue
			xTemp(1) = sum
			ier = 0
			call CheckStoichiometry(k,xTemp,ier)
			if(ier.eq.1)then			! failed stoichiometry test
				newtonStep = newtonStep/2.0d0
				!if(newtonStep.gt.1.0d-6)go to 5001			! try again with this smaller stepsize
				if(newtonStep.gt.1.0d-15)go to 5001			! try again with this smaller stepsize
				!stepsize is too small. Abort this phase
				!if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)then
				izero = 1
				write(*,*)'stepSize = ',newtonStep
				write(*,*)'This is too small. Something must be wrong in the code. Line 497 of Compute4EQUIL'
				write(*,*)'Phase = ',k
				write(*,*)phName(k)
				do j1 = 1,numPhCo(K)
					write(*,*)phCoName(k,j1),xPhCo(k,j1),xTemp(j1)
					end do
					!5026			continue
				pause ' hit return to continue'
				izero = 1		! set error flag
				return
				endif

5010  			continue
			end do
		!write(*,*)'NewtonStep = ',newtonStep
		!check for convergence
		!for convergence we check to see whether the x value is < tolerance*X value      
      		jj = NVAR-numEqEquns+1		! index of solution vector in array XX
      		imiss=0
      		do i = 1,numEqEquns
	      		ii=ipoint(i)
			!Newtonstep scales the solution vector so we can step towards the solution in smaller steps
	      		ALLX(1,ii)=ALLX(1,ii) + newtonstep*XX(i,jj)
			!AllX(1,ii) is the current value of the dependent variables. xx(i,jj) is the solution vector (variable).
			!When variable is less than a fraction(TolNewton) of the dependent variable, we have converged
	      		if(dabs(xx(i,jj)).gt.dabs(TolNewton*Allx(1,ii)))then
	      			if(allFractl(ii).eq.0)then		! only count imiss if the phase is not fractionating (seeGibbs_FileInput.f line 520 or so)
	      		      		imiss=imiss+1
					endif
				endif
			end do
			!5110  		continue
 		if(iopen.eq.1)then
			do i = 1,iopen
				openFlux(iOpena(i)) = allX(1,openPtr - 1 + i)
				end do
				!5120			continue
 			endif	! end open system check
		!if imiss>0 then one or more values failed tolerance test.  Rearrange
		!things and go do another iteration
      		go to 3200
      		endif	! end iNewton = 1
! ________________________________________________________________
      
! _______________________Gibbs Method__________________________
      	if(inewton.eq.0)then
!     		XX IS A numEqEqunsx(NVAR-numEqEquns) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		COMPUTE new values FOR EACH VARIABLE
!     		note that the old variables DELX and SDELX are now defined as
!           	DELX(i)  = ALLX(1,i) - ALLX(5,i)  :current minus last finite difference
!           	SDELX(i) = ALLX(1,i) - ALLX(2,i)  :current minus starting
!     		Definition of ALLX variable:
!     		1,50   current value
!     		2,50     starting value
!     		3,50   ref at start of contour
!     		4,50   user selected reference
!     		5,50   previous finite diff point
!     		6,50     (not used)
!     		Store new values for each independent variable (monitor parameter)
      		DO J=1,NVAR-numEqEquns
			ii=mon(j)
			ALLX(1,ii)=ALLX(1,ii) + DELTAX(J)
			end do
			!3100  		continue
!     		compute new values for each dependent variable
      		Do i = 1,numEqEquns
	      		ii=ipoint(i)
	      		DO J=1,NVAR-numEqEquns
		      		ALLX(1,ii)=ALLX(1,ii) + XX(i,J) * DELTAX(J)
				end do
			end do
			!3110  		continue    
      		endif       ! gibbs/newton

! ________________________________________________________________
3200  	continue

! ______________BOTH_____________________________________________
!     OUTPUT DERIVATIVES FOR THIS FINITE DIFFERENCE ITERATION
	if(iLong(2).eq.1)then
		call WriteOutJacobian()
	   	endif

! ----------------------------
      	call setTPX
!     	WRITE TO TERMINAL AND DISK
      	IF(iLong(4).eq.1)then
         	write(12,*)
         	write(12,*)
         	write(12,*)'Calculated Deltas'
         	write(12,*)'                      .       Value  .This increment.   Running sum'
         	do i = 1,nvar
         		write(12,4060)I,VN1(I),VN2(I),ALLX(1,i),ALLX(1,i)-ALLX(5,i),ALLX(1,i)-ALLX(2,i)
4060     		FORMAT(I3,1X,A2,A16,F15.4,3E15.5)
!4060     		FORMAT(I3,1X,A2,A16,3F15.4)
			end do
			!4018     	continue
		endif

	if(iLong(2)+iLong(3)+iLong(4)+iLong(5)+iLong(6).ge.1)then
		write(*,*) "Type 0 (enter) to continue with next step"
		write(*,*) "Type 1 (enter) to abort"
		write(*,*) "Type -1 (enter) to continue without extended output"
!		pause "Hit return to continue with next step"
		read(*,*)izero
		if(izero.eq.1)then
			go to 9999
		    	endif
		if(izero.lt.0)then
			iLong(2) = 0
			iLong(3) = 0
			iLong(4) = 0
			iLong(5) = 0
			iLong(6) = 0
			endif
		endif

!     DONE WITH THIS STEP
      	if(inewton.eq.1.and.newtonStepsOutput.eq.1)then
		call NewtonStepSave
		write(*,*)newtonCounter,imiss
		endif

      	if(inewton.eq.1.and.imiss.gt.0)then
 		newtoncounter = newtoncounter+1
		if(newtoncounter.gt.1000)then
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			write(*,*)'Newton Counter > 1000 iterations. Aborting....'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			izero = 1
			go to 9999
			endif		
        	go to 100	! loop for another iteration
		endif
	
9999	continue		! exit routine

      	RETURN
      	END


! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Compute4MDF3(izero)
	use MatrixArrays
      implicit none

! 	Routine to calculate Jacobian for MDF (OS) model phases
!	This routine is a generalized solution where
!	(a) Reactions are stored in array ARX
!		EQ phases have an entry for every phase component
!		OS phases have an entry for nphco-1 "components" and they are all PhCo(ind)-PhCo(dep)
!		Note that there are npc - 1 such independent equations
!	Basically, the method should provide a matrix solution that allows any number of phases to be EQ or OS
!
!	The variance of this system should be 2 + #OS phases so we can solve for changes in all phases given the growth of the OS phase
!		For example dMph/dMOSphase
!	Original version written June 8, 2018 by F. Spear
!	Generalized version written Oct 31, 2018 by F. Spear
! 	Note that this is the same code as in C2OMPUTE except 
!		generalized for OS phases
!		No contouring option
!		No generic "new" variables
! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "Singular.inc"
!	include "Solute.inc"
!	include "Tangent.inc"
	include "MatrixSingular.inc"
! c*****************************************
!     local variables
      integer*4 i,ii,j,jj,k,L,LL,Lcont,jcol,irow,ier,j1,kCur,			&
     	 izero,imiss,nsolve,newtoncounter,ivol
      REAL*8 DETA,R,SUM,Pi,TEMPT,TEMPP,TEMPX,TEMPH,TEMPS,TEMPK,xTemp(20),allXTemp(100)
      CHARACTER vn51(100)*8
      DATA R,Pi/8.3144D0,3.141592654D0/


! --------------------------------------------

      IF(iLong(2).eq.1)then

         write(12,*)' '
         write(12,*)'------------------------------------------------------'
         write(12,*)'------------------------------------------------------'
         write(12,*)' Variable list '
         do i = 1,nvar
        	write(12,15)I,VN1(I),VN2(I)
15      	FORMAT(I3,2X,A2,A16)
		end do
         write(12,*)'NEQ=',NEQ,'  NVAR=',NVAR
         write(12,*)
         write(12,*)'VARIANCE =',NVAR-NEQ
         write(12,*)
         write(12,*)'MONITOR PARAMETERS    Delta X for monitors'
         if(nvar-neq.gt.0)then
         do i = 1,nvar-neq
         	write(12,16)MON(I),VN1(MON(I)),VN2(MON(I)),Deltax(i)*NSTEP
16       	FORMAT(I3,2X,A2,A16,F15.8)
		end do
	!write(12,*)' DELTA Xi FOR MONITOR PHASES'
	!write(12,16)(DELTAX(I)*NSTEP,I=1,NVAR-NEQ)
	!16       FORMAT(10F12.3)
         else
         write(12,*)' No monitor parameters (variance <= 0)'
         endif
         write(12,*)'NUMBER OF FINITE DIFF. ITERATIONS=',NSTEP
         write(12,*)
         write(12,*)'INITIAL MINERAL DATA'
         write(12,17)TC,PB
17       format(' T START (DEG C)=',F10.2,'   P START=',F10.2)
	if(PFluidSwitch.gt.0)write(12,18)PFluid
18	format(' PFluid = ',F10.2)
	endif
90    CONTINUE
!     LCONT IS COUNTER IN CONTOUR BLOCK
      LCONT=0
      newtoncounter=0

! -----------------------------------------------------------
!                             
100   CONTINUE


      TK=TC+273.15D0

!     THESE CALLS CALCULATE THE THERMODYNAMIC
!     PROPERTIES OF THE PHASES OF INTEREST

!     SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      izero=0
      ivol=0            !compute everything
      K=0
!	This call will callculate all thermodynamic properties - the TP dependent ones and the composition dependent ones
!	We need this when we're using T or P as dependent variables. But not when they're independent (because they're already set)
      CALL ALLKduTPX(ivol,izero)

      if(izero.gt.0)go to 9999

!     COMPUTE PHASE VOLUMES
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
            	VP0(K)=MP0(K)*VMOL(K)*10.D0
		end do

!     OUTPUT THERMO DATA TO DISK

      IF(iLong(3).eq.1)then
350   	CONTINUE
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
          	write(12,*)
          	write(12,200)MINREC(K),SITMUL(K),PHNAME(K)
200       	FORMAT(' ',I4,F8.1,4X,A32)
          	write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
          	IF(IMASS.EQ.1) then
          		write(12,507)MP0(K),VP0(K)
          		endif
507       	FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
          	write(12,501)(phCoName(k,j),J=1,numPhCo(K))
501       	FORMAT('          ',12(A4,6X))
          	write(12,504)(xPhCo(k,j),J=1,numPhCo(K))
504       	FORMAT('  COMP ',12(F10.4))
          	write(12,505)(Dexp(lnAct(k,j)),J=1,numPhCo(K))
505       	FORMAT(' Activ  ',12(E10.3))
          	write(12,506)(lnAct(k,j),J=1,numPhCo(K))
506       	FORMAT('  Ln(a) ',12(F10.5))
          	write(12,514)(hPhCoZero(k,j),J=1,numPhCo(K))
514       	Format(' hPhCoZero',12E15.7)
          	write(12,515)(HATTP(k,j),J=1,numPhCo(K))
515       	Format(' HATTP',12E15.7)
          	write(12,509)(sPhCoZero(k,j),J=1,numPhCo(K))
509       	Format(' SZERO',12F15.5)
          	write(12,510)(SATTP(k,j),J=1,numPhCo(K))
510       	Format(' SATTP',12F15.5)
          	write(12,508)(vPhCoZero(k,j),J=1,numPhCo(K))
508       	Format(' VZERO',12F15.5)
          	write(12,513)(VATTP(k,j),J=1,numPhCo(K))
513       	Format(' VATTP',12F15.5)
          	write(12,516)(GATTP(k,j),J=1,numPhCo(K))
516       	Format(' GATTP',12E15.7)
          	write(12,*)'dudTPX ARRAY:'
          	write(12,*)'     dudT   .    dudP   .    dudX2  .    dudX3  .    dudX4...'
           	do j = 1,numPhCo(K)
           		write(12,511)dudTPX(k,j,1),dudTPX(k,j,2),(dudTPX(k,j,L),L=3,numPhCo(K)+1)
511       		FORMAT(15E12.5)
			end do
408       	format(' ',a5,6E10.3)
		end do

	    endif


!     THIS PART OF THE PROGRAM SETS UP THE MASTER MATRIX

!     ZERO MASTER ARRAY
      	DO J=1,NVAR
		DO I=1,NEQ
			AA(I,J)=0.0D0
			end do
		end do
!     The partial derivatives in matrix dudTPX are:
!     In storage space 1: dudTPX(jj,1):
!     Temperature:      /T = - S(P,T,X)  -- Partial molar Entropy at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state entropy, plus the entropy of mixing plus any exess entropy
!                                              terms.
!     In storage space 2: dudTPX(jj,2):
!     Pressure   :      /P =   V(P,T,X)  -- Partial molar Volume at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state volume, plus any exess volume
!                                              terms.
!     In storage space 3 to (numPhCo(K)+1): dudTPX(jj,3), dudTPX(jj,3), dudTPX(jj,3) etc:
!           Note that there are only numPhCo(K)-1 partial derivatives and the subscripts range
!           from 3 to numPhCo(K)+1.
!     Composition:      /X = (ideal)/X + (margules)/X + (reciprocal)/X

!       
!     The matrix ArxMDF(numMDFrxn,npc) contains the stoichiometric coefficients for each phase component in the MDF reactions
!	Note that these equations define the phase components for the minerals in terms of the phase components of the grain boundary (matrix)
!     in each of the reactions, subscript 1 refers to reaction number, subscript 2 refers to the phase component

!	The MDF equations take the form
!	0 = (µ2 - µ1)(in MDF phase) - (µ2 - µ1)(in GB phase)
!	0 = (µ3 - µ1)(in MDF phase) - (µ3 - µ1)(in GB phase)
!	etc.
!	where µ1, µ2, µ3 in the GB phase are defined by the linearly independent reactions in array Arx

!     The matrix AA therefore contains in each cell terms such as
!     
!       (n(i)*u(i))/T
!       (n(i)*u(i))/P
!       (n(i)*u(i))/X2
!       (n(i)*u(i))/X3

!     for each independent variable (where n(j,i) are the stoichiometric coefficients for phase component i in reaction j)
!       Each row of AA contains coefficients for each of the reactions.

!     The only catch to the calculation of AA is that the matrix dudPTX is
!     not stored in true matrix form.  That is, to save space the independent variables are listed for each phase only.
!     This simply requires a bit of extra bookkeeping during the matrix multiplication



!	Here is what the MDF equations look like -- stored in array ArxMDF
!	  Number of MDF equations (numMDFrxn) =   1
!	     SiGB    AlGB    MgGB    FeGB    abQz    Prp     Alm 
!	    0.000   0.000   3.000  -3.000   0.000  -1.000   1.000

	do i = 1,numMDFrxn		! loop for each reaction
		TempT = 0.0d0
		TempP = 0.0d0
		j = 0
		do kcur = 1,numPh		! Loop for every phase
			k = asmCurrent(kCur)
			do jj = 1,numPhCo(k)	! Loop for every phase component in every phase
				j = j+1
				! dT and dP coefficient
				TempT = TempT + ArxMDF(i,j)*(dudTPX(k,jj,1))	! index 1 is partial molar S
				TempP = TempP + ArxMDF(i,j)*(dudTPX(k,jj,2))	! index 2 is partial molar V
				end do
!1015				continue
			end do
!1014			continue
	      	AA(i,1)=TEMPT
	      	AA(i,2)=TEMPP
		end do
!1010		continue

!	Do dX terms
	do i = 1,numMDFrxn
		L = 2				! starts in the third column (2 = dT, dP)
		jj = 0
		do kcur = 1,numPh
			k = asmCurrent(kCur)
			if(numPhCo(K).eq.1)then
				jj = jj + 1
				go to 1021	! there are no dX terms for a phase of fixed composition
				endif
			Do LL = 3,numPhCo(K)+1     	! this loops on the number of independent composition derivatives
				tempX = 0.D0
				do j=1,numPhCo(K)          	! this loops on the number of phase components in each K phase
					TEMPX=TEMPX + ArxMDF(i,jj+j)*dudTPX(k,j,LL)
					end do
	!				1023  	continue
				L=L+1
				AA(i,L)=TEMPX
				end do
!				1022  	continue
			jj = jj + numPhCo(k)
1021  			continue
			end do
		end do
!1020  		continue




!     	SET UP MASS BALANCE EQUATIONS
      	if(imass.eq.0)go to 1100
!     	FIRST SET UP dM Columns (M=MOLE FRACTION OF PHASE)
!     	JCOL COUNTS WHERE dM TERMS ARE
      	JCOL = TANDP + NX
!     	loop through each phase (there is a dM term for each phase)
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
	      	JCOL = JCOL + 1
		!LOOP THROUGH ALL mass balance equations for each SYSTEM COMPONENT
      		call bulkcomp				! get total moles of each system component (moles(i) in common block)
      		DO i = 1,nc			! loop through every mass balance equation
      			iRow = i + numMDFrxn			! iRow counts the row number for the mass balance equations
      			sum=0.D0
			!LOOP THROUGH ALL PHASE COMPONENTS for this phase
      			DO j=1,numPhCo(K)
				sum = sum + comp(k,J,I)*xPhCo(k,J)
				end do
				!1120  	continue
	
		      	AA(IROW,JCOL)=SUM
		      	YY(IROW) = -(moles(i) - molesStart(i))			! molesStart(i) is the bulk composition we want to fit
      										! moles(i) is the current bulk comp based on phase comp and M(k)
										! when they are equal, we have the correct X and M for the desired bulk comp
										! openFlux are open system components either dependent or independent variables
			end do
			!1105  	CONTINUE
		end do
		!1110  	CONTINUE	!loop for each phase

!     	NOW SET UP COEFFICIENTS FOR dX TERMS in mass balance equations
!     	Note that the dependent dX term is removed as above
      	JCOL=TANDP					!     	JCOL is the column for the dX term
!     	FIRST LOOP THROUGH ALL PHASES TO FIND independent dX TERMS
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
		IF(numPhCo(K).EQ.1) GO TO 1130
		!IF HERE THEN WE HAVE FOUND A PHASE WITH dX TERMS (I.E. MORE THAN 1 COMPONENT)
		!LOOP THROUGH ALL dX TERMS FOR THIS PHASE
			DO J=2,numPhCo(K)
				JCOL=JCOL+1
				!LOOP THROUGH ALL MASS BALANCE EQUATIONS FOR EACH dX TERM
				DO I=1,NC
					!IROW= I + NRX
					IROW= I + numMDFrxn
					AA(IROW,JCOL)=MP0(K)*(COMP(k,J,I)-COMP(k,1,I))
					end do
					!1150	continue
				end do
				!1140  	CONTINUE
1130  		CONTINUE
		end do

1100  continue   ! end of mass balance equations


!     If we are in MODE=NEWTON then we need to set up a data vector that
!     contains -(H - TS + R T lnKeq) for each reaction
!     (note there is no V term here because H and S are calculated at T and P
      	if(inewton.eq.1) then
!		ii = np-nrx      ! note that reaction coefficients are stored in the end of array ARX
		DO i=1,numMDFrxn
			TEMPH=0.D0
			TEMPS=0.D0
			TEMPK=0.D0
			jj = 0
			Do kCur = 1,numPh
				k = asmCurrent(kCur)
				do J = 1,numPhCo(k)
					jj = jj + 1
					TEMPH =TEMPH + ArxMDF(i,jj)*HatTP(k,j)		! ∆H reaction
					TEMPS =TEMPS + ArxMDF(i,jj)*SatTP(k,j)		! ∆H reaction
					TEMPK =TEMPK + ArxMDF(i,jj)*lnAct(k,j)		! ∆H reaction
					end do
				end do
				!1415    	continue      
 			YY(i)=-(TEMPH - TK * TEMPS + R*TK*TEMPK)
			end do
			!1410     	continue
		endif      

!     	Rearrange MASTER matrix AA
!     	MOVE MONITOR PARAMETERS TO RIGHT SIDE OF EQUNS IN ARRAY A
	nsolve=nvar
	if(inewton.eq.1)nsolve=nvar+1
	DO I=1,NVAR-NEQ
		DO J=1,NEQ
		  	A(J,NEQ+I)=-AA(J,MON(I))
			end do
		end do
!2100	continue
!     	FILL LEFT SIDE OF EQUATIONS IN MATRIX A FROM AA
      	DO I=1,NEQ
 	     	DO J=1,NEQ
			A(J,I)=AA(J,IPOINT(I))
			end do
		end do
!2110	continue
      	if(inewton.eq.1)then
!           Fill Nvar+1=nsolve array element in A with the Y vector      
            DO J=1,NEQ
        	A(J,Nsolve)=YY(J)
		end do	
!2115        A(J,Nsolve)=YY(J)
            endif
            

!     	Output matrix and modified matrix, if desired
	if(iLong(5).eq.1)then
2330  		CONTINUE
      		write(12,*)' '
      		write(12,*)' '
      		write(12,*)' MASTER MATRIX'
	        write(12,2321)(VN1(j),VN2(j),J=1,nvar),'  YY(i)'
2321  		FORMAT('     ',30(A2,A4,6x))
	      	DO I=1,NEQ
			!write(12,*)'--------------------'
			write(12,2320)(AA(I,J),J=1,NVAR),YY(I)
			end do
			!2310  		continue
2320  		FORMAT(' ',30E12.4)
		do i = 1,nvar-neq
			vn51(neq+i) = trim(vn1(mon(i)))//vn2(mon(i))
			end do
			!2311		continue
		do i = 1,neq
			vn51(i) = trim(vn1(ipoint(i)))//vn2(ipoint(i))
			end do
			!2312		continue
		write(12,*)' '
		write(12,*)' '
		write(12,*)' MODIFIED MASTER MATRIX'
	        write(12,2313)(VN51(j),J=1,nvar),' YY(i)'
2313  		FORMAT('     ',30(A6,6x))
		DO I=1,NEQ
			!write(12,*)'--------------------'
			write(12,2320)(A(I,J),J=1,nsolve)
			end do
			!2315  		continue
		endif

!     Find solution
      	DETA=0.D0
      	IER=0
      	CALL REDUCE (NEQ,NSOLVE,DETA,IER)
      	IF (IER.EQ.1) then
		izero=1
!		if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)go to 2509
		WRITE(*,*)' ************ ERROR **************************'
		WRITE(*,*)' Matrix failed to invert in SUBROUTINE REDUCE'
		WRITE(*,*)' The most likely cause is that the variables'
		WRITE(*,*)' specified as monitor parameters are not linearly'
		WRITE(*,*)' independent.  Note that only the number of'
		WRITE(*,*)' intensive variables equal to the phase rule'
		WRITE(*,*)' variance of the assemblage are independent.'
		WRITE(*,*)' Return to previous menu and try again.'
		write(*,*)' '
		pause 'Hit return to continue...'
		go to 9999
		endif

2509	continue


	if(matrixIsSingular.eq.1)go to 9999

	if(iLong(2).eq.1)then
      		write(12,*)'   '
      		write(12,*)' RESULTS OF THIS FINITE DIFFERENCE ITERATION'
      		write(12,*)'DETERMINANT=',DETA
		endif


! -------------------------------
!     	Set ALLX(5,L) equal to current values of ALLX(1,L)
!     	Note that below ALLX(1,L) is incremented by DELX
      	DO L=1,NVAR
		ALLX(5,L)=ALLX(1,L)
		end do
		!3000  	ALLX(5,L)=ALLX(1,L)

! ____________________Newtons method___________________________
      	if(inewton.eq.1)then
!     		XX IS A NEQx(NVAR-NEQ+1) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		The last column (NVAR-NEQ+1) is the solution to the x for Newtons method
!     		Calculate the change in the variables
!     		New variables are computed as
!     		Xnew = Xold + X
! 		Note that we are doing these calculations with the values of the "independent" variables constant
! 		i.e. dxi = 0. For example, if it is univariant and we are doing calculations at constant P then dP = 0
! 		This way, we only need to solve for the values of the "dependent" variables (dT, dXi, etc.)
!

		! skip this newton step stuff -- it seems to be causing trouble
!		go to 5010

!       	Composition variables can only range from 0-1. 
!		If the solution vector component (xx) is too large, it may result in the change in mole fraction (xPhCo(i) exceeding this range
!		This block of code checks this and adjusts newtonStep to ensure this doesn't happen
!		newtonStep = 1.0d0

		newtonStep = newtonStepMax
		jCol = NVAR-NEQ+1	! index of solution vector in array XX
		do i = 1,nvar
			allXtemp(i) = allX(1,i)	! load up allXTemp
			end do
			!5002		continue
! ----		loop to here for newtonStep iteration --------
5001		continue
		do i = 1,NEQ
			ii=ipoint(i)
			!adjust the values of AllXtemp that we have just calculated
			!note that this changes all dependent values, not just the composition ones
			!we could make this a wee bit faster if we just calculated the composition ones but we're not sure which of them are dependent
			!jcol is the column where the solution vector resides in XX
			ALLXtemp(ii)=ALLX(1,ii) + newtonstep*XX(i,jCol)
			end do
			!5004		continue
!		now pick out the composition values and store in xTemp
		l = TandP		! this is where the independent composition derivatives start
		Do kCur = 1,numPh
			k = asmCurrent(kCur)
			if(numPhCo(k).eq.1)go to 5010	! no composition derivatives for 1 component phases
			sum = 1.0d0		! this is for the dependent composition derivative (the first one)
			DO J = 2,numPhCo(K) ! loop from second to last phase component for this mineral
				l = l + 1		! index of independent composition variable in AllX
				!j indexes the phase components of phase K
				xTemp(j) = allXTemp(l)
				sum = sum - xTemp(j)
				end do
				!5020  		continue
			xTemp(1) = sum
			ier = 0
			call CheckStoichiometry(k,xTemp,ier)
			if(ier.eq.1)then			! failed stoichiometry test
				newtonStep = newtonStep/2.0d0
				!if(newtonStep.gt.1.0d-6)go to 5001			! try again with this smaller stepsize
				if(newtonStep.gt.1.0d-15)go to 5001			! try again with this smaller stepsize
				!stepsize is too small. Abort this phase
				!if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)then
				izero = 1
				write(*,*)'stepSize = ',newtonStep
				write(*,*)'This is too small. Something must be wrong in the code. Line 1144of Compute4MDF3'
				write(*,*)'Phase = ',k
				write(*,*)phName(k)
				do j1 = 1,numPhCo(K)
					write(*,*)phCoName(k,j1),xPhCo(k,j1),xTemp(j1)
					end do
					!5026			continue
				pause ' hit return to continue'
				izero = 1		! set error flag
				return
				endif

5010  			continue
			end do
		!write(*,*)'NewtonStep = ',newtonStep
		!check for convergence
		!for convergence we check to see whether the x value is < tolerance*X value      
      		jj = NVAR-NEQ+1		! index of solution vector in array XX
      		imiss=0
      		do i = 1,NEQ
	      		ii=ipoint(i)
			!Newtonstep scales the solution vector so we can step towards the solution in smaller steps
	      		ALLX(1,ii)=ALLX(1,ii) + newtonstep*XX(i,jj)
			!AllX(1,ii) is the current value of the dependent variables. xx(i,jj) is the solution vector (variable).
			!When variable is less than a fraction(TolNewton) of the dependent variable, we have converged
	      		if(dabs(xx(i,jj)).gt.dabs(TolNewton*Allx(1,ii)))then
	      			if(allFractl(ii).eq.0)then		! only count imiss if the phase is not fractionating (seeGibbs_FileInput.f line 520 or so)
	      		      		imiss=imiss+1
					endif
				endif
			end do
			!5110  		continue
 		if(iopen.eq.1)then
			do i = 1,iopen
				openFlux(iOpena(i)) = allX(1,openPtr - 1 + i)
				end do
				!5120			continue
 			endif	! end open system check
		!if imiss>0 then one or more values failed tolerance test.  Rearrange
		!things and go do another iteration
      		go to 3200
      		endif	! end iNewton = 1
! ________________________________________________________________
      
! _______________________Gibbs Method__________________________
      	if(inewton.eq.0)then
!     		XX IS A NEQx(NVAR-NEQ) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		COMPUTE new values FOR EACH VARIABLE
!     		note that the old variables DELX and SDELX are now defined as
!           	DELX(i)  = ALLX(1,i) - ALLX(5,i)  :current minus last finite difference
!           	SDELX(i) = ALLX(1,i) - ALLX(2,i)  :current minus starting
!     		Definition of ALLX variable:
!     		1,50   current value
!     		2,50     starting value
!     		3,50   ref at start of contour
!     		4,50   user selected reference
!     		5,50   previous finite diff point
!     		6,50     (not used)
!     		Store new values for each independent variable (monitor parameter)
      		DO J=1,NVAR-NEQ
			ii=mon(j)
			ALLX(1,ii)=ALLX(1,ii) + DELTAX(J)
			end do
			!3100  		continue
!     		compute new values for each dependent variable
      		Do i = 1,NEQ
	      		ii=ipoint(i)
	      		DO J=1,NVAR-NEQ
		      		ALLX(1,ii)=ALLX(1,ii) + XX(i,J) * DELTAX(J)
				end do
			end do
			!3110  		continue    
      		endif       ! gibbs/newton

! ________________________________________________________________
3200  	continue

! ______________BOTH_____________________________________________
!     OUTPUT DERIVATIVES FOR THIS FINITE DIFFERENCE ITERATION
	if(iLong(2).eq.1)then
		call WriteOutJacobian()
	   	endif

! ----------------------------
      	call setTPX
!     	WRITE TO TERMINAL AND DISK
      	IF(iLong(4).eq.1)then
         	write(12,*)
         	write(12,*)
         	write(12,*)'Calculated Deltas'
         	write(12,*)'                      .       Value  .This increment.   Running sum'
         	do i = 1,nvar
         		write(12,4060)I,VN1(I),VN2(I),ALLX(1,i),ALLX(1,i)-ALLX(5,i),ALLX(1,i)-ALLX(2,i)
4060     		FORMAT(I3,1X,A2,A16,F15.4,3E15.5)
!4060     		FORMAT(I3,1X,A2,A16,3F15.4)
			end do
			!4018     	continue
		endif

	if(iLong(2)+iLong(3)+iLong(4)+iLong(5)+iLong(6).ge.1)then
		write(*,*) "Type 0 (enter) to continue with next step"
		write(*,*) "Type 1 (enter) to abort"
		write(*,*) "Type -1 (enter) to continue without extended output"
!		pause "Hit return to continue with next step"
		read(*,*)izero
		if(izero.eq.1)then
			go to 9999
		    	endif
		if(izero.lt.0)then
			iLong(2) = 0
			iLong(3) = 0
			iLong(4) = 0
			iLong(5) = 0
			iLong(6) = 0
			endif
		endif

!     DONE WITH THIS STEP
      	if(inewton.eq.1.and.newtonStepsOutput.eq.1)then
		call NewtonStepSave
		write(*,*)newtonCounter,imiss
		endif

      	if(inewton.eq.1.and.imiss.gt.0)then
 		newtoncounter = newtoncounter+1
		if(newtoncounter.gt.1000)then
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			write(*,*)'Newton Counter > 1000 iterations. Aborting....'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			izero = 1
			go to 9999
			endif		
        	go to 100	! loop for another iteration
		endif
	
9999	continue		! exit routine

      	RETURN
      	END


! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Compute4MDF2(izero)
	use MatrixArrays
      implicit none

! 	Routine to calculate Jacobian for MDF (OS) model phases
!	This routine is a generalized solution where
!	(a) Reactions are stored in array ARX
!		EQ phases have an entry for every phase component
!		OS phases have an entry for nphco-1 "components" and they are all PhCo(ind)-PhCo(dep)
!		Note that there are npc - 1 such independent equations
!	Basically, the method should provide a matrix solution that allows any number of phases to be EQ or OS
!
!	The variance of this system should be 2 + #OS phases so we can solve for changes in all phases given the growth of the OS phase
!		For example dMph/dMOSphase
!	Original version written June 8, 2018 by F. Spear
!	Generalized version written Oct 31, 2018 by F. Spear
! 	Note that this is the same code as in C2OMPUTE except 
!		generalized for OS phases
!		No contouring option
!		No generic "new" variables
! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "Singular.inc"
!	include "Solute.inc"
!	include "Tangent.inc"
	include "MatrixSingular.inc"
! c*****************************************
!     local variables
      integer*4 i,ii,j,jj,k,L,LL,Lcont,jcol,irow,ier,j1,kCur,			&
     	 izero,imiss,nsolve,newtoncounter,ivol
      REAL*8 DETA,R,SUM,Pi,TEMPT,TEMPP,TEMPX,TEMPH,TEMPS,TEMPK,xTemp(20),allXTemp(100)
      CHARACTER vn51(100)*8
      DATA R,Pi/8.3144D0,3.141592654D0/


! --------------------------------------------

      IF(iLong(2).eq.1)then

         write(12,*)' '
         write(12,*)'------------------------------------------------------'
         write(12,*)'------------------------------------------------------'
         write(12,*)' Variable list '
         do 11 i = 1,nvar
         write(12,15)I,VN1(I),VN2(I)
15       FORMAT(I3,2X,A2,A16)
11       continue
         write(12,*)'NEQ=',NEQ,'  NVAR=',NVAR
         write(12,*)
         write(12,*)'VARIANCE =',NVAR-NEQ
         write(12,*)
         write(12,*)'MONITOR PARAMETERS    Delta X for monitors'
         if(nvar-neq.gt.0)then
         do 12 i = 1,nvar-neq
         write(12,16)MON(I),VN1(MON(I)),VN2(MON(I)),Deltax(i)*NSTEP
16       FORMAT(I3,2X,A2,A16,F15.8)
12       continue
!         write(12,*)' DELTA Xi FOR MONITOR PHASES'
!         write(12,16)(DELTAX(I)*NSTEP,I=1,NVAR-NEQ)
!16       FORMAT(10F12.3)
         else
         write(12,*)' No monitor parameters (variance <= 0)'
         endif
         write(12,*)'NUMBER OF FINITE DIFF. ITERATIONS=',NSTEP
         write(12,*)
         write(12,*)'INITIAL MINERAL DATA'
         write(12,17)TC,PB
17       format(' T START (DEG C)=',F10.2,'   P START=',F10.2)
	if(PFluidSwitch.gt.0)write(12,18)PFluid
18	format(' PFluid = ',F10.2)
	endif
90    CONTINUE
!     LCONT IS COUNTER IN CONTOUR BLOCK
      LCONT=0
      newtoncounter=0

! -----------------------------------------------------------
!                             
100   CONTINUE


      TK=TC+273.15D0

!     THESE CALLS CALCULATE THE THERMODYNAMIC
!     PROPERTIES OF THE PHASES OF INTEREST

!     SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      izero=0
      ivol=0            !compute everything
      K=0
! x       write(*,*)' In c2ompute #3-calling duTPX'
! x       pause
!	This call will callculate all thermodynamic properties - the TP dependent ones and the composition dependent ones
!	We need this when we're using T or P as dependent variables. But not when they're independent (because they're already set)
      CALL ALLKduTPX(ivol,izero)
! x       write(*,*)' In c2ompute #4-back from duTPX'
! x        pause
      if(izero.gt.0)go to 9999

!     COMPUTE PHASE VOLUMES
	Do 340 kCur = 1,numPh
	k = asmCurrent(kCur)
            VP0(K)=MP0(K)*VMOL(K)*10.D0
340         continue

!     OUTPUT THERMO DATA TO DISK

      IF(iLong(3).eq.1)then
350       CONTINUE
	Do 390 kCur = 1,numPh
	k = asmCurrent(kCur)
          write(12,*)
          write(12,200)MINREC(K),SITMUL(K),PHNAME(K)
200       FORMAT(' ',I4,F8.1,4X,A32)
          write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
          IF(IMASS.EQ.1) then
             write(12,507)MP0(K),VP0(K)
          endif
507       FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
          write(12,501)(phCoName(k,j),J=1,numPhCo(K))
501       FORMAT('          ',12(A4,6X))
          write(12,504)(xPhCo(k,j),J=1,numPhCo(K))
504       FORMAT('  COMP ',12(F10.4))
          write(12,505)(Dexp(lnAct(k,j)),J=1,numPhCo(K))
505       FORMAT(' Activ  ',12(E10.3))
          write(12,506)(lnAct(k,j),J=1,numPhCo(K))
506       FORMAT('  Ln(a) ',12(F10.5))
          write(12,514)(hPhCoZero(k,j),J=1,numPhCo(K))
514       Format(' hPhCoZero',12E15.7)
          write(12,515)(HATTP(k,j),J=1,numPhCo(K))
515       Format(' HATTP',12E15.7)
          write(12,509)(sPhCoZero(k,j),J=1,numPhCo(K))
509       Format(' SZERO',12F15.5)
          write(12,510)(SATTP(k,j),J=1,numPhCo(K))
510       Format(' SATTP',12F15.5)
          write(12,508)(vPhCoZero(k,j),J=1,numPhCo(K))
508       Format(' VZERO',12F15.5)
          write(12,513)(VATTP(k,j),J=1,numPhCo(K))
513       Format(' VATTP',12F15.5)
          write(12,516)(GATTP(k,j),J=1,numPhCo(K))
516       Format(' GATTP',12E15.7)
          write(12,*)'dudTPX ARRAY:'
          write(12,*)'     dudT   .    dudP   .    dudX2  .    dudX3  .    dudX4...'
            do 512 j = 1,numPhCo(K)
           	write(12,511)dudTPX(k,j,1),dudTPX(k,j,2),(dudTPX(k,j,L),L=3,numPhCo(K)+1)
511       	FORMAT(15E12.5)
512   	continue
408       format(' ',a5,6E10.3)
551       CONTINUE
390       CONTINUE

	endif
395       CONTINUE


!     THIS PART OF THE PROGRAM SETS UP THE MASTER MATRIX

!     ZERO MASTER ARRAY
      DO 1000 J=1,NVAR
      DO 1000 I=1,NEQ
1000  AA(I,J)=0.0D0

!     The partial derivatives in matrix dudTPX are:
!     In storage space 1: dudTPX(jj,1):
!     Temperature:      /T = - S(P,T,X)  -- Partial molar Entropy at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state entropy, plus the entropy of mixing plus any exess entropy
!                                              terms.
!     In storage space 2: dudTPX(jj,2):
!     Pressure   :      /P =   V(P,T,X)  -- Partial molar Volume at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state volume, plus any exess volume
!                                              terms.
!     In storage space 3 to (numPhCo(K)+1): dudTPX(jj,3), dudTPX(jj,3), dudTPX(jj,3) etc:
!           Note that there are only numPhCo(K)-1 partial derivatives and the subscripts range
!           from 3 to numPhCo(K)+1.
!     Composition:      /X = (ideal)/X + (margules)/X + (reciprocal)/X

!       
!     The matrix ArxMDF(numMDFrxn,npc) contains the stoichiometric coefficients for each phase component in the MDF reactions
!	Note that these equations define the phase components for the minerals in terms of the phase components of the grain boundary (matrix)
!     in each of the reactions, subscript 1 refers to reaction number, subscript 2 refers to the phase component

!	The MDF equations take the form
!	0 = (µ2 - µ1)(in MDF phase) - (µ2 - µ1)(in GB phase)
!	0 = (µ3 - µ1)(in MDF phase) - (µ3 - µ1)(in GB phase)
!	etc.
!	where µ1, µ2, µ3 in the GB phase are defined by the linearly independent reactions in array Arx

!     The matrix AA therefore contains in each cell terms such as
!     
!       (n(i)*u(i))/T
!       (n(i)*u(i))/P
!       (n(i)*u(i))/X2
!       (n(i)*u(i))/X3

!     for each independent variable (where n(j,i) are the stoichiometric coefficients for phase component i in reaction j)
!       Each row of AA contains coefficients for each of the reactions.

!     The only catch to the calculation of AA is that the matrix dudPTX is
!     not stored in true matrix form.  That is, to save space the independent variables are listed for each phase only.
!     This simply requires a bit of extra bookkeeping during the matrix multiplication



!	Here is what the MDF equations look like -- stored in array ArxMDF
!	  Number of MDF equations (numMDFrxn) =   1
!	     SiGB    AlGB    MgGB    FeGB    abQz    Prp     Alm 
!	    0.000   0.000   3.000  -3.000   0.000  -1.000   1.000

	do 1010 i = 1,numMDFrxn		! loop for each reaction
		TempT = 0.0d0
		TempP = 0.0d0
		j = 0
		do 1014 kcur = 1,numPh		! Loop for every phase
			k = asmCurrent(kCur)
			do 1015 jj = 1,numPhCo(k)	! Loop for every phase component in every phase
				j = j+1
				! dT and dP coefficient
				TempT = TempT + ArxMDF(i,j)*(dudTPX(k,jj,1))	! index 1 is partial molar S
				TempP = TempP + ArxMDF(i,j)*(dudTPX(k,jj,2))	! index 2 is partial molar V
1015				continue
1014			continue
	      	AA(i,1)=TEMPT
	      	AA(i,2)=TEMPP
1010		continue

!	Do dX terms
	do 1020 i = 1,numMDFrxn
	L = 2				! starts in the third column (2 = dT, dP)
	jj = 0
	do 1021 kcur = 1,numPh
	k = asmCurrent(kCur)
	if(numPhCo(K).eq.1)then
		jj = jj + 1
		go to 1021	! there are no dX terms for a phase of fixed composition
		endif
	Do 1022 LL = 3,numPhCo(K)+1     	! this loops on the number of independent composition derivatives
	tempX = 0.D0
	do 1023 j=1,numPhCo(K)          	! this loops on the number of phase components in each K phase
	TEMPX=TEMPX + ArxMDF(i,jj+j)*dudTPX(k,j,LL)
1023  	continue
      	L=L+1
      	AA(i,L)=TEMPX
1022  	continue
	jj = jj + numPhCo(k)
1021  	continue
1020  	continue




!     	SET UP MASS BALANCE EQUATIONS
      	if(imass.eq.0)go to 1100
!     	FIRST SET UP dM Columns (M=MOLE FRACTION OF PHASE)
!     	JCOL COUNTS WHERE dM TERMS ARE
      	JCOL = TANDP + NX
!     	loop through each phase (there is a dM term for each phase)
	Do 1110 kCur = 1,numPh
	k = asmCurrent(kCur)
      	JCOL = JCOL + 1
!     	LOOP THROUGH ALL mass balance equations for each SYSTEM COMPONENT
      	call bulkcomp				! get total moles of each system component (moles(i) in common block)
      	DO 1105 i = 1,nc			! loop through every mass balance equation
      	iRow = i + numMDFrxn			! iRow counts the row number for the mass balance equations
      	sum=0.D0
!     	LOOP THROUGH ALL PHASE COMPONENTS for this phase
      	DO 1120 j=1,numPhCo(K)
1120  	sum = sum + comp(k,J,I)*xPhCo(k,J)
      	AA(IROW,JCOL)=SUM
!      	YY(IROW) = -(moles(i) - openFlux(i) - molesStart(i))	! molesStart(i) is the bulk composition we want to fit
      	YY(IROW) = -(moles(i) - molesStart(i))			! molesStart(i) is the bulk composition we want to fit
      								! moles(i) is the current bulk comp based on phase comp and M(k)
								! when they are equal, we have the correct X and M for the desired bulk comp
								! openFlux are open system components either dependent or independent variables
1105  	CONTINUE
1110  	CONTINUE	!loop for each phase

!     	NOW SET UP COEFFICIENTS FOR dX TERMS in mass balance equations
!     	Note that the dependent dX term is removed as above
      	JCOL=TANDP					!     	JCOL is the column for the dX term
!     	FIRST LOOP THROUGH ALL PHASES TO FIND independent dX TERMS
	Do 1130 kCur = 1,numPh
	k = asmCurrent(kCur)
      	IF(numPhCo(K).EQ.1) GO TO 1130
!     	IF HERE THEN WE HAVE FOUND A PHASE WITH dX TERMS (I.E. MORE THAN 1 COMPONENT)
!     	LOOP THROUGH ALL dX TERMS FOR THIS PHASE
      	DO 1140 J=2,numPhCo(K)
      	JCOL=JCOL+1
!     	LOOP THROUGH ALL MASS BALANCE EQUATIONS FOR EACH dX TERM
      	DO 1150 I=1,NC
!      	IROW= I + NRX
      	IROW= I + numMDFrxn
	AA(IROW,JCOL)=MP0(K)*(COMP(k,J,I)-COMP(k,1,I))		! this is where we remove the dependent variable from the mass balance equations
1150	continue
1140  	CONTINUE
1130  	CONTINUE

1100  continue   ! end of mass balance equations


!     If we are in MODE=NEWTON then we need to set up a data vector that
!     contains -(H - TS + R T lnKeq) for each reaction
!     (note there is no V term here because H and S are calculated at T and P
      	if(inewton.eq.1) then
!		ii = np-nrx      ! note that reaction coefficients are stored in the end of array ARX
		DO 1410 i=1,numMDFrxn
		TEMPH=0.D0
		TEMPS=0.D0
		TEMPK=0.D0
		jj = 0
		Do 1415 kCur = 1,numPh
		k = asmCurrent(kCur)
		do 1415 J = 1,numPhCo(k)
		jj = jj + 1
		TEMPH =TEMPH + ArxMDF(i,jj)*HatTP(k,j)		! ∆H reaction
		TEMPS =TEMPS + ArxMDF(i,jj)*SatTP(k,j)		! ∆H reaction
		TEMPK =TEMPK + ArxMDF(i,jj)*lnAct(k,j)		! ∆H reaction
1415    	continue      
 		YY(i)=-(TEMPH - TK * TEMPS + R*TK*TEMPK)
1410     	continue
         	endif      

!     	Rearrange MASTER matrix AA
!     	MOVE MONITOR PARAMETERS TO RIGHT SIDE OF EQUNS IN ARRAY A
	nsolve=nvar
	if(inewton.eq.1)nsolve=nvar+1
	DO 2100 I=1,NVAR-NEQ
	DO 2100 J=1,NEQ
  	A(J,NEQ+I)=-AA(J,MON(I))
2100	continue
!     	FILL LEFT SIDE OF EQUATIONS IN MATRIX A FROM AA
      	DO 2110 I=1,NEQ
      	DO 2110 J=1,NEQ
	A(J,I)=AA(J,IPOINT(I))
2110	continue
      	if(inewton.eq.1)then
!           Fill Nvar+1=nsolve array element in A with the Y vector      
            DO 2115 J=1,NEQ
2115        A(J,Nsolve)=YY(J)
            endif
            

!     	Output matrix and modified matrix, if desired
	if(iLong(5).eq.1)then
2330  		CONTINUE
      		write(12,*)' '
      		write(12,*)' '
      		write(12,*)' MASTER MATRIX'
	        write(12,2321)(VN1(j),VN2(j),J=1,nvar),'  YY(i)'
2321  		FORMAT('     ',30(A2,A4,6x))
	      	DO 2310 I=1,NEQ
!		write(12,*)'--------------------'
		write(12,2320)(AA(I,J),J=1,NVAR),YY(I)
2310  		continue
2320  		FORMAT(' ',30E12.4)
		do 2311 i = 1,nvar-neq
		vn51(neq+i) = trim(vn1(mon(i)))//vn2(mon(i))
2311		continue
		do 2312 i = 1,neq
		vn51(i) = trim(vn1(ipoint(i)))//vn2(ipoint(i))
2312		continue
		write(12,*)' '
		write(12,*)' '
		write(12,*)' MODIFIED MASTER MATRIX'
	        write(12,2313)(VN51(j),J=1,nvar),' YY(i)'
2313  		FORMAT('     ',30(A6,6x))
		DO 2315 I=1,NEQ
	!         write(12,*)'--------------------'
		write(12,2320)(A(I,J),J=1,nsolve)
2315  		continue
		endif

!     Find solution
      	DETA=0.D0
      	IER=0
      	CALL REDUCE (NEQ,NSOLVE,DETA,IER)
      	IF (IER.EQ.1) then
		izero=1
!		if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)go to 2509
		WRITE(*,*)' ************ ERROR **************************'
		WRITE(*,*)' Matrix failed to invert in SUBROUTINE REDUCE'
		WRITE(*,*)' The most likely cause is that the variables'
		WRITE(*,*)' specified as monitor parameters are not linearly'
		WRITE(*,*)' independent.  Note that only the number of'
		WRITE(*,*)' intensive variables equal to the phase rule'
		WRITE(*,*)' variance of the assemblage are independent.'
		WRITE(*,*)' Return to previous menu and try again.'
		write(*,*)' '
		pause 'Hit return to continue...'
		go to 9999
		endif

2509	continue


	if(matrixIsSingular.eq.1)go to 9999

	if(iLong(2).eq.1)then
      		write(12,*)'   '
      		write(12,*)' RESULTS OF THIS FINITE DIFFERENCE ITERATION'
      		write(12,*)'DETERMINANT=',DETA
		endif


! -------------------------------
!     	Set ALLX(5,L) equal to current values of ALLX(1,L)
!     	Note that below ALLX(1,L) is incremented by DELX
      	DO 3000 L=1,NVAR
3000  	ALLX(5,L)=ALLX(1,L)

! ____________________Newtons method___________________________
      	if(inewton.eq.1)then
!     		XX IS A NEQx(NVAR-NEQ+1) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		The last column (NVAR-NEQ+1) is the solution to the x for Newtons method
!     		Calculate the change in the variables
!     		New variables are computed as
!     		Xnew = Xold + X
! 		Note that we are doing these calculations with the values of the "independent" variables constant
! 		i.e. dxi = 0. For example, if it is univariant and we are doing calculations at constant P then dP = 0
! 		This way, we only need to solve for the values of the "dependent" variables (dT, dXi, etc.)
!

		! skip this newton step stuff -- it seems to be causing trouble
!		go to 5010

!       	Composition variables can only range from 0-1. 
!		If the solution vector component (xx) is too large, it may result in the change in mole fraction (xPhCo(i) exceeding this range
!		This block of code checks this and adjusts newtonStep to ensure this doesn't happen
!		newtonStep = 1.0d0

		newtonStep = newtonStepMax
		jCol = NVAR-NEQ+1	! index of solution vector in array XX
		do 5002 i = 1,nvar
		allXtemp(i) = allX(1,i)	! load up allXTemp
5002		continue
! ----		loop to here for newtonStep iteration --------
5001		continue
		do 5004 i = 1,NEQ
		ii=ipoint(i)
!     		adjust the values of AllXtemp that we have just calculated
!		note that this changes all dependent values, not just the composition ones
!		we could make this a wee bit faster if we just calculated the composition ones but we're not sure which of them are dependent
!		jcol is the column where the solution vector resides in XX
		ALLXtemp(ii)=ALLX(1,ii) + newtonstep*XX(i,jCol)
5004		continue
!		now pick out the composition values and store in xTemp
		l = TandP		! this is where the independent composition derivatives start
!		DO 5010 k = 1,numPh
		Do 5010 kCur = 1,numPh
		k = asmCurrent(kCur)
		if(numPhCo(k).eq.1)go to 5010	! no composition derivatives for 1 component phases
		sum = 1.0d0		! this is for the dependent composition derivative (the first one)
		DO 5020 J = 2,numPhCo(K) ! loop from second to last phase component for this mineral
		l = l + 1		! index of independent composition variable in AllX
!		j indexes the phase components of phase K
		xTemp(j) = allXTemp(l)
		sum = sum - xTemp(j)
5020  		continue
      		xTemp(1) = sum
		ier = 0
		call CheckStoichiometry(k,xTemp,ier)
		if(ier.eq.1)then			! failed stoichiometry test
			newtonStep = newtonStep/2.0d0
			!if(newtonStep.gt.1.0d-6)go to 5001			! try again with this smaller stepsize
			if(newtonStep.gt.1.0d-15)go to 5001			! try again with this smaller stepsize
!			stepsize is too small. Abort this phase
!			if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)then
			izero = 1
!				return
!				endif
			write(*,*)'stepSize = ',newtonStep
			write(*,*)'This is too small. Something must be wrong in the code. Line 1766 of Compute4MDF2'
			write(*,*)'Phase = ',k
			write(*,*)phName(k)
			do 5026 j1 = 1,numPhCo(K)
			write(*,*)phCoName(k,j1),xPhCo(k,j1),xTemp(j1)
5026			continue
			pause ' hit return to continue'
			izero = 1		! set error flag
			return
			endif

5010  		continue
!		write(*,*)'NewtonStep = ',newtonStep
!     		check for convergence
!     		for convergence we check to see whether the x value is < tolerance*X value      
      		jj = NVAR-NEQ+1		! index of solution vector in array XX
      		imiss=0
      		do 5110 i = 1,NEQ
      		ii=ipoint(i)
!     		Newtonstep scales the solution vector so we can step towards the solution in smaller steps
      		ALLX(1,ii)=ALLX(1,ii) + newtonstep*XX(i,jj)
!		AllX(1,ii) is the current value of the dependent variables. xx(i,jj) is the solution vector (variable).
!		When variable is less than a fraction(TolNewton) of the dependent variable, we have converged
      		if(dabs(xx(i,jj)).gt.dabs(TolNewton*Allx(1,ii)))then
      			if(allFractl(ii).eq.0)then		! only count imiss if the phase is not fractionating (seeGibbs_FileInput.f line 520 or so)
      		      		imiss=imiss+1
				endif
			endif
5110  		continue
 		if(iopen.eq.1)then
			do 5120 i = 1,iopen
			openFlux(iOpena(i)) = allX(1,openPtr - 1 + i)
5120			continue
 		endif
!     		if imiss>0 then one or more values failed tolerance test.  Rearrange
!     		things and go do another iteration
      		go to 3200
      		endif
! ________________________________________________________________
      
! _______________________Gibbs Method__________________________
      	if(inewton.eq.0)then
!     		XX IS A NEQx(NVAR-NEQ) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		COMPUTE new values FOR EACH VARIABLE
!     		note that the old variables DELX and SDELX are now defined as
!           	DELX(i)  = ALLX(1,i) - ALLX(5,i)  :current minus last finite difference
!           	SDELX(i) = ALLX(1,i) - ALLX(2,i)  :current minus starting
!     		Definition of ALLX variable:
!     		1,50   current value
!     		2,50     starting value
!     		3,50   ref at start of contour
!     		4,50   user selected reference
!     		5,50   previous finite diff point
!     		6,50     (not used)
!     		Store new values for each independent variable (monitor parameter)
      		DO 3100 J=1,NVAR-NEQ
      		ii=mon(j)
      		ALLX(1,ii)=ALLX(1,ii) + DELTAX(J)
3100  		continue
!     		compute new values for each dependent variable
      		Do 3110 i = 1,NEQ
      		ii=ipoint(i)
      		DO 3110 J=1,NVAR-NEQ
      		ALLX(1,ii)=ALLX(1,ii) + XX(i,J) * DELTAX(J)
3110  		continue    
      		endif       ! gibbs/newton

! ________________________________________________________________
3200  	continue

! ______________BOTH_____________________________________________
!     OUTPUT DERIVATIVES FOR THIS FINITE DIFFERENCE ITERATION
	if(iLong(2).eq.1)then
		call WriteOutJacobian()
	   	endif

! ----------------------------
      	call setTPX
!     	WRITE TO TERMINAL AND DISK
      	IF(iLong(4).eq.1)then
         	write(12,*)
         	write(12,*)
         	write(12,*)'Calculated Deltas'
         	write(12,*)'                      .       Value  .This increment.   Running sum'
         	do 4018 i = 1,nvar
         	write(12,4060)I,VN1(I),VN2(I),ALLX(1,i),ALLX(1,i)-ALLX(5,i),ALLX(1,i)-ALLX(2,i)
4060     		FORMAT(I3,1X,A2,A16,F15.4,3E15.5)
!4060     	FORMAT(I3,1X,A2,A16,3F15.4)
4018     	continue
		endif

	if(iLong(2)+iLong(3)+iLong(4)+iLong(5)+iLong(6).ge.1)then
		write(*,*) "Type 0 (enter) to continue with next step"
		write(*,*) "Type 1 (enter) to abort"
		write(*,*) "Type -1 (enter) to continue without extended output"
!		pause "Hit return to continue with next step"
		read(*,*)izero
		if(izero.eq.1)then
			go to 9999
		    	endif
		if(izero.lt.0)then
			iLong(2) = 0
			iLong(3) = 0
			iLong(4) = 0
			iLong(5) = 0
			iLong(6) = 0
			endif
		endif

!     DONE WITH THIS STEP
      	if(inewton.eq.1.and.newtonStepsOutput.eq.1)then
		call NewtonStepSave
		write(*,*)newtonCounter,imiss
		endif

      	if(inewton.eq.1.and.imiss.gt.0)then
 		newtoncounter = newtoncounter+1
		if(newtoncounter.gt.1000)then
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			write(*,*)'Newton Counter > 1000 iterations. Aborting....'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			izero = 1
			go to 9999
			endif		
        	go to 100	! loop for another iteration
		endif
	
9999	continue		! exit routine

      	RETURN
      	END


! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Compute4MDF4(izero)
	use MatrixArrays
      implicit none

!	This is modified after the old Compute4MDF2/3
!	The difference is that here we calculate the dependent compositional variable in the GB using the linearly independent reactions
!		rather than using the stoichiometric constraint
!	In this way, the GB composition doesn't sum to 1.0

! 	Routine to calculate Jacobian for MDF (OS) model phases
!	This routine is a generalized solution where
!	(a) Reactions are stored in array ARX
!		EQ phases have an entry for every phase component
!		OS phases have an entry for nphco-1 "components" and they are all PhCo(ind)-PhCo(dep)
!		Note that there are npc - 1 such independent equations
!	Basically, the method should provide a matrix solution that allows any number of phases to be EQ or OS
!
!	The variance of this system should be 2 + #OS phases so we can solve for changes in all phases given the growth of the OS phase
!		For example dMph/dMOSphase
!	Original version written June 8, 2018 by F. Spear
!	Generalized version written Oct 31, 2018 by F. Spear
! 	Note that this is the same code as in C2OMPUTE except 
!		generalized for OS phases
!		No contouring option
!		No generic "new" variables
! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "Singular.inc"
!	include "Solute.inc"
!	include "Tangent.inc"
	include "MatrixSingular.inc"
! c*****************************************
!     local variables
      integer*4 i,ii,j,jj,k,L,LL,Lcont,jcol,irow,ier,j1,kCur,			&
     	 izero,imiss,nsolve,newtoncounter,ivol
      REAL*8 DETA,R,SUM,Pi,TEMPT,TEMPP,TEMPX,TEMPH,TEMPS,TEMPK,xTemp(20),allXTemp(100),TempDeltaGB
      CHARACTER vn51(100)*8
      DATA R,Pi/8.3144D0,3.141592654D0/


! --------------------------------------------

      IF(iLong(2).eq.1)then

         write(12,*)' '
         write(12,*)'------------------------------------------------------'
         write(12,*)'------------------------------------------------------'
         write(12,*)' Variable list '
         do i = 1,nvar
        	write(12,15)I,VN1(I),VN2(I)
15      	FORMAT(I3,2X,A2,A16)
		end do
         write(12,*)'NEQ=',NEQ,'  NVAR=',NVAR
         write(12,*)
         write(12,*)'VARIANCE =',NVAR-NEQ
         write(12,*)
         write(12,*)'MONITOR PARAMETERS    Delta X for monitors'
         if(nvar-neq.gt.0)then
         do i = 1,nvar-neq
         	write(12,16)MON(I),VN1(MON(I)),VN2(MON(I)),Deltax(i)*NSTEP
16       	FORMAT(I3,2X,A2,A16,F15.8)
		end do
         else
         write(12,*)' No monitor parameters (variance <= 0)'
         endif
         write(12,*)'NUMBER OF FINITE DIFF. ITERATIONS=',NSTEP
         write(12,*)
         write(12,*)'INITIAL MINERAL DATA'
         write(12,17)TC,PB
17       format(' T START (DEG C)=',F10.2,'   P START=',F10.2)
	if(PFluidSwitch.gt.0)write(12,18)PFluid
18	format(' PFluid = ',F10.2)
	endif
90    CONTINUE
!     LCONT IS COUNTER IN CONTOUR BLOCK
      LCONT=0
      newtoncounter=0

! -----------------------------------------------------------
!                             
100   CONTINUE


      TK=TC+273.15D0

!     THESE CALLS CALCULATE THE THERMODYNAMIC
!     PROPERTIES OF THE PHASES OF INTEREST

!     SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      izero=0
      ivol=0            !compute everything
      K=0
!	This call will callculate all thermodynamic properties - the TP dependent ones and the composition dependent ones
!	We need this when we're using T or P as dependent variables. But not when they're independent (because they're already set)
      CALL ALLKduTPX(ivol,izero)

      if(izero.gt.0)go to 9999

!     COMPUTE PHASE VOLUMES
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
            	VP0(K)=MP0(K)*VMOL(K)*10.D0
		end do

!     OUTPUT THERMO DATA TO DISK

      IF(iLong(3).eq.1)then
350   	CONTINUE
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
          	write(12,*)
          	write(12,200)MINREC(K),SITMUL(K),PHNAME(K)
200       	FORMAT(' ',I4,F8.1,4X,A32)
          	write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
          	IF(IMASS.EQ.1) then
          		write(12,507)MP0(K),VP0(K)
          		endif
507       	FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
          	write(12,501)(phCoName(k,j),J=1,numPhCo(K))
501       	FORMAT('          ',12(A4,6X))
          	write(12,504)(xPhCo(k,j),J=1,numPhCo(K))
504       	FORMAT('  COMP ',12(F10.4))
          	write(12,505)(Dexp(lnAct(k,j)),J=1,numPhCo(K))
505       	FORMAT(' Activ  ',12(E10.3))
          	write(12,506)(lnAct(k,j),J=1,numPhCo(K))
506       	FORMAT('  Ln(a) ',12(F10.5))
          	write(12,514)(hPhCoZero(k,j),J=1,numPhCo(K))
514       	Format(' hPhCoZero',12E15.7)
          	write(12,515)(HATTP(k,j),J=1,numPhCo(K))
515       	Format(' HATTP',12E15.7)
          	write(12,509)(sPhCoZero(k,j),J=1,numPhCo(K))
509       	Format(' SZERO',12F15.5)
          	write(12,510)(SATTP(k,j),J=1,numPhCo(K))
510       	Format(' SATTP',12F15.5)
          	write(12,508)(vPhCoZero(k,j),J=1,numPhCo(K))
508       	Format(' VZERO',12F15.5)
          	write(12,513)(VATTP(k,j),J=1,numPhCo(K))
513       	Format(' VATTP',12F15.5)
          	write(12,516)(GATTP(k,j),J=1,numPhCo(K))
516       	Format(' GATTP',12E15.7)
          	write(12,*)'dudTPX ARRAY:'
          	write(12,*)'     dudT   .    dudP   .    dudX2  .    dudX3  .    dudX4...'
           	do j = 1,numPhCo(K)
           		write(12,511)dudTPX(k,j,1),dudTPX(k,j,2),(dudTPX(k,j,L),L=3,numPhCo(K)+1)
511       		FORMAT(15E12.5)
			end do
408       	format(' ',a5,6E10.3)
		end do

	    endif


!     THIS PART OF THE PROGRAM SETS UP THE MASTER MATRIX

!     ZERO MASTER ARRAY
      	DO J=1,NVAR
		DO I=1,NEQ
			AA(I,J)=0.0D0
			end do
		end do
!     The partial derivatives in matrix dudTPX are:
!     In storage space 1: dudTPX(jj,1):
!     Temperature:      /T = - S(P,T,X)  -- Partial molar Entropy at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state entropy, plus the entropy of mixing plus any exess entropy
!                                              terms.
!     In storage space 2: dudTPX(jj,2):
!     Pressure   :      /P =   V(P,T,X)  -- Partial molar Volume at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state volume, plus any exess volume
!                                              terms.
!     In storage space 3 to (numPhCo(K)+1): dudTPX(jj,3), dudTPX(jj,3), dudTPX(jj,3) etc:
!           Note that there are only numPhCo(K)-1 partial derivatives and the subscripts range
!           from 3 to numPhCo(K)+1.
!     Composition:      /X = (ideal)/X + (margules)/X + (reciprocal)/X

!       
!     The matrix ArxMDF(numMDFrxn,npc) contains the stoichiometric coefficients for each phase component in the MDF reactions
!	Note that these equations define the phase components for the minerals in terms of the phase components of the grain boundary (matrix)
!     in each of the reactions, subscript 1 refers to reaction number, subscript 2 refers to the phase component

!	The MDF equations take the form
!	0 = (µ2 - µ1)(in MDF phase) - (µ2 - µ1)(in GB phase)
!	0 = (µ3 - µ1)(in MDF phase) - (µ3 - µ1)(in GB phase)
!	etc.
!	where µ1, µ2, µ3 in the GB phase are defined by the linearly independent reactions in array Arx

!     The matrix AA therefore contains in each cell terms such as
!     
!       (n(i)*u(i))/T
!       (n(i)*u(i))/P
!       (n(i)*u(i))/X2
!       (n(i)*u(i))/X3

!     for each independent variable (where n(j,i) are the stoichiometric coefficients for phase component i in reaction j)
!       Each row of AA contains coefficients for each of the reactions.

!     The only catch to the calculation of AA is that the matrix dudPTX is
!     not stored in true matrix form.  That is, to save space the independent variables are listed for each phase only.
!     This simply requires a bit of extra bookkeeping during the matrix multiplication



!	Here is what the MDF equations look like -- stored in array ArxMDF
!	  Number of MDF equations (numMDFrxn) =   1
!	     SiGB    AlGB    MgGB    FeGB    abQz    Prp     Alm 
!	    0.000   0.000   3.000  -3.000   0.000  -1.000   1.000

	do i = 1,numMDFrxn		! loop for each reaction
		TempT = 0.0d0
		TempP = 0.0d0
		j = 0
		do kcur = 1,numPh		! Loop for every phase
			k = asmCurrent(kCur)
			do jj = 1,numPhCo(k)	! Loop for every phase component in every phase
				j = j+1
				! dT and dP coefficient
				TempT = TempT + ArxMDF(i,j)*(dudTPX(k,jj,1))	! index 1 is partial molar S
				TempP = TempP + ArxMDF(i,j)*(dudTPX(k,jj,2))	! index 2 is partial molar V
				end do
			end do
	      	AA(i,1)=TEMPT
	      	AA(i,2)=TEMPP
		end do

!	Do dX terms
	do i = 1,numMDFrxn
		L = 2				! starts in the third column (2 = dT, dP)
		jj = 0
		do kcur = 1,numPh
			k = asmCurrent(kCur)
			if(numPhCo(K).eq.1)then
				jj = jj + 1
				go to 1021	! there are no dX terms for a phase of fixed composition
				endif
			Do LL = 3,numPhCo(K)+1     	! this loops on the number of independent composition derivatives
				tempX = 0.D0
				do j=1,numPhCo(K)          	! this loops on the number of phase components in each K phase
					TEMPX=TEMPX + ArxMDF(i,jj+j)*dudTPX(k,j,LL)
					end do
				L=L+1
				AA(i,L)=TEMPX
				end do
			jj = jj + numPhCo(k)
1021  			continue
			end do
		end do




!     	SET UP MASS BALANCE EQUATIONS
	! imass should always be 1
!      	if(imass.eq.0)go to 1100
!     	FIRST SET UP dM Columns (M=MOLE FRACTION OF PHASE)
!     	JCOL COUNTS WHERE dM TERMS ARE
      	JCOL = TANDP + NX
!     	loop through each phase (there is a dM term for each phase)
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
	      	JCOL = JCOL + 1
		!LOOP THROUGH ALL mass balance equations for each SYSTEM COMPONENT
      		call bulkcomp				! get total moles of each system component (moles(i) in common block)
      		DO i = 1,nc			! loop through every mass balance equation
      			iRow = i + numMDFrxn			! iRow counts the row number for the mass balance equations
      			sum=0.D0
			!LOOP THROUGH ALL PHASE COMPONENTS for this phase
      			DO j=1,numPhCo(K)
				sum = sum + comp(k,J,I)*xPhCo(k,J)
				end do
	
		      	AA(IROW,JCOL)=SUM
		      	YY(IROW) = -(moles(i) - molesStart(i))			! molesStart(i) is the bulk composition we want to fit
      										! moles(i) is the current bulk comp based on phase comp and M(k)
										! when they are equal, we have the correct X and M for the desired bulk comp
										! openFlux are open system components either dependent or independent variables
			end do
		end do

!     	NOW SET UP COEFFICIENTS FOR dX TERMS in mass balance equations
!     	Note that the dependent dX term is removed as above
      	JCOL=TANDP					!     	JCOL is the column for the dX term
!     	FIRST LOOP THROUGH ALL PHASES TO FIND independent dX TERMS
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
		IF(numPhCo(K).EQ.1) GO TO 1130
		!IF HERE THEN WE HAVE FOUND A PHASE WITH dX TERMS (I.E. MORE THAN 1 COMPONENT)
		!LOOP THROUGH ALL dX TERMS FOR THIS PHASE
			DO J=2,numPhCo(K)
				JCOL=JCOL+1
				!LOOP THROUGH ALL MASS BALANCE EQUATIONS FOR EACH dX TERM
				DO I=1,NC
					!IROW= I + NRX
					IROW= I + numMDFrxn
					AA(IROW,JCOL)=MP0(K)*(COMP(k,J,I)-COMP(k,1,I))
					end do
				end do
1130  		CONTINUE
		end do

	! end of mass balance equations


!     If we are in MODE=NEWTON then we need to set up a data vector that
!     contains -(H - TS + R T lnKeq) for each reaction
!     (note there is no V term here because H and S are calculated at T and P
!	iNewton is always 1 for MDF routines	
!      	if(inewton.eq.1) then
!		ii = np-nrx      ! note that reaction coefficients are stored in the end of array ARX
	DO i=1,numMDFrxn
		TEMPH=0.D0
		TEMPS=0.D0
		TEMPK=0.D0
		jj = 0
		Do kCur = 1,numPh
			k = asmCurrent(kCur)
			do J = 1,numPhCo(k)
				jj = jj + 1
				TEMPH =TEMPH + ArxMDF(i,jj)*HatTP(k,j)		! ∆H reaction
				TEMPS =TEMPS + ArxMDF(i,jj)*SatTP(k,j)		! ∆H reaction
				TEMPK =TEMPK + ArxMDF(i,jj)*lnAct(k,j)		! ∆H reaction
				end do
			end do
		YY(i)=-(TEMPH - TK * TEMPS + R*TK*TEMPK)
		end do

!     	Rearrange MASTER matrix AA
!     	MOVE MONITOR PARAMETERS TO RIGHT SIDE OF EQUNS IN ARRAY A
	nsolve=nvar
	nsolve=nvar+1
	DO I=1,NVAR-NEQ
		DO J=1,NEQ
		  	A(J,NEQ+I)=-AA(J,MON(I))
			end do
		end do
!     	FILL LEFT SIDE OF EQUATIONS IN MATRIX A FROM AA
      	DO I=1,NEQ
 	     	DO J=1,NEQ
			A(J,I)=AA(J,IPOINT(I))
			end do
		end do
	! Fill Nvar+1=nsolve array element in A with the Y vector      
	DO J=1,NEQ
		A(J,Nsolve)=YY(J)
		end do	
            
!     	Output matrix and modified matrix, if desired
	if(iLong(5).eq.1)then
2330  		CONTINUE
      		write(12,*)' '
      		write(12,*)' '
      		write(12,*)' MASTER MATRIX'
	        write(12,2321)(VN1(j),VN2(j),J=1,nvar),'  YY(i)'
2321  		FORMAT('     ',30(A2,A4,6x))
	      	DO I=1,NEQ
			!write(12,*)'--------------------'
			write(12,2320)(AA(I,J),J=1,NVAR),YY(I)
			end do
			!2310  		continue
2320  		FORMAT(' ',30E12.4)
		do i = 1,nvar-neq
			vn51(neq+i) = trim(vn1(mon(i)))//vn2(mon(i))
			end do
			!2311		continue
		do i = 1,neq
			vn51(i) = trim(vn1(ipoint(i)))//vn2(ipoint(i))
			end do
			!2312		continue
		write(12,*)' '
		write(12,*)' '
		write(12,*)' MODIFIED MASTER MATRIX'
	        write(12,2313)(VN51(j),J=1,nvar),' YY(i)'
2313  		FORMAT('     ',30(A6,6x))
		DO I=1,NEQ
			!write(12,*)'--------------------'
			write(12,2320)(A(I,J),J=1,nsolve)
			end do
			!2315  		continue
		endif

!     Find solution
      	DETA=0.D0
      	IER=0
      	CALL REDUCE (NEQ,NSOLVE,DETA,IER)
      	IF (IER.EQ.1) then
		izero=1
!		if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)go to 2509
		WRITE(*,*)' ************ ERROR **************************'
		WRITE(*,*)' Matrix failed to invert in SUBROUTINE REDUCE'
		WRITE(*,*)' The most likely cause is that the variables'
		WRITE(*,*)' specified as monitor parameters are not linearly'
		WRITE(*,*)' independent.  Note that only the number of'
		WRITE(*,*)' intensive variables equal to the phase rule'
		WRITE(*,*)' variance of the assemblage are independent.'
		WRITE(*,*)' Return to previous menu and try again.'
		write(*,*)' '
		pause 'Hit return to continue...'
		go to 9999
		endif

2509	continue


	if(matrixIsSingular.eq.1)go to 9999

	if(iLong(2).eq.1)then
      		write(12,*)'   '
      		write(12,*)' RESULTS OF THIS FINITE DIFFERENCE ITERATION'
      		write(12,*)'DETERMINANT=',DETA
		endif


! -------------------------------
!     	Set ALLX(5,L) equal to current values of ALLX(1,L)
!     	Note that below ALLX(1,L) is incremented by DELX
      	DO L=1,NVAR
		ALLX(5,L)=ALLX(1,L)
		end do

! ____________________Newtons method___________________________
!      	if(inewton.eq.1)then
!	iNewton is always 1 in MDF routines
!     	XX IS A NEQx(NVAR-NEQ+1) MATRIX WITH THE SOLUTIONS OF THE
!     	SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     	The last column (NVAR-NEQ+1) is the solution to the x for Newtons method
!     	Calculate the change in the variables
!     	New variables are computed as
!     	Xnew = Xold + X
! 	Note that we are doing these calculations with the values of the "independent" variables constant
! 	i.e. dxi = 0. For example, if it is univariant and we are doing calculations at constant P then dP = 0
! 	This way, we only need to solve for the values of the "dependent" variables (dT, dXi, etc.)
!
!       Composition variables can only range from 0-1. 
!	If the solution vector component (xx) is too large, it may result in the change in mole fraction (xPhCo(i) exceeding this range
!	This block of code checks this and adjusts newtonStep to ensure this doesn't happen
!	newtonStep = 1.0d0

	newtonStep = newtonStepMax
	jCol = NVAR-NEQ+1	! index of solution vector in array XX
	do i = 1,nvar
		allXtemp(i) = allX(1,i)	! load up allXTemp
		end do
! ----		loop to here for newtonStep iteration --------
5001	continue
	do i = 1,NEQ
		ii=ipoint(i)
		!adjust the values of AllXtemp that we have just calculated
		!note that this changes all dependent values, not just the composition ones
		!we could make this a wee bit faster if we just calculated the composition ones but we're not sure which of them are dependent
		!jcol is the column where the solution vector resides in XX
		ALLXtemp(ii)=ALLX(1,ii) + newtonstep*XX(i,jCol)
		end do

!	now pick out the composition values and store in xTemp
	l = TandP		! this is where the independent composition derivatives start

	! The definitions of the GB components in relation to the solids are in the linearly independent reaction set
	! We know the change in moles of the solid phases, so we can calculate the change in moles of the GB components
	! We already have all of the independent GB compositions so we only need to calculate the dependent one
	! do the grain boundary			
	tempDeltaGB = 0.0d0
	ii = np - numMDFRxn	! this is one shy of the index for the first linearly independent reaction
	j1 = numPhCo(1)		! number of phase components in the GB phase
	do i = 1,numMDFRxn		! number of linearly independent reactions
	      	J = 2	 	! this is where the ∆M terms start for the solid phases in the deltax array
		do jj = 1,numPh-1		! loop on all of the solid phases
			j = j + 1	! index for deltaX array
		!	xTemp(1) = xTemp(1) + MP0(jj)*Arx(ii+i,j)
		! not quite right. We want the CHANGE in the moles, which is newtonstep*xx(??) for this iteration
			tempDeltaGB = tempDeltaGB + DeltaX(j)*Arx(ii+i,j1+jj)			
			end do
		end do
	xTemp(1) = xPhCo(1,1) + tempDeltaGB
	GBdependentX = xTemp(1)		! temporary storage to use in subroutine SetTPX
      	IF(iLong(4).eq.1)then
		write(12,*)'***** Dependent GB composition = ',xTemp(1)			
		endif
		
	Do kCur = 2,numPh	! do all the phases except the grain boundary
		k = asmCurrent(kCur)
		if(numPhCo(k).eq.1)go to 5010	! no composition derivatives for 1 component phases
		! here is where the dependent variable is calculated as X(dep) = 1 - sum of others
		sum = 1.0d0		! this is for the dependent composition derivative (the first one)
		DO J = 2,numPhCo(K) ! loop from second to last phase component for this mineral
			l = l + 1		! index of independent composition variable in AllX
			!j indexes the phase components of phase K
			xTemp(j) = allXTemp(l)
			sum = sum - xTemp(j)
			end do
			!5020  		continue
		xTemp(1) = sum
		ier = 0
		call CheckStoichiometry(k,xTemp,ier)
		if(ier.eq.1)then			! failed stoichiometry test
			newtonStep = newtonStep/2.0d0
			!if(newtonStep.gt.1.0d-6)go to 5001			! try again with this smaller stepsize
			if(newtonStep.gt.1.0d-15)go to 5001			! try again with this smaller stepsize
			!stepsize is too small. Abort this phase
			izero = 1
			write(*,*)'stepSize = ',newtonStep
			write(*,*)'This is too small. Something must be wrong in the code. Line 1144of Compute4MDF3'
			write(*,*)'Phase = ',k
			write(*,*)phName(k)
			do j1 = 1,numPhCo(K)
				write(*,*)phCoName(k,j1),xPhCo(k,j1),xTemp(j1)
				end do
			pause ' hit return to continue'
			izero = 1		! set error flag
			return
			endif

5010  		continue
		end do
	!write(*,*)'NewtonStep = ',newtonStep
	!check for convergence
	!for convergence we check to see whether the x value is < tolerance*X value      
	jj = NVAR-NEQ+1		! index of solution vector in array XX
	imiss=0
	do i = 1,NEQ
		ii=ipoint(i)
		!Newtonstep scales the solution vector so we can step towards the solution in smaller steps
		ALLX(1,ii)=ALLX(1,ii) + newtonstep*XX(i,jj)
		!AllX(1,ii) is the current value of the dependent variables. xx(i,jj) is the solution vector (variable).
		!When variable is less than a fraction(TolNewton) of the dependent variable, we have converged
		if(dabs(xx(i,jj)).gt.dabs(TolNewton*Allx(1,ii)))then
			if(allFractl(ii).eq.0)then		! only count imiss if the phase is not fractionating (seeGibbs_FileInput.f line 520 or so)
				imiss=imiss+1
				endif
			endif
		end do

	!if imiss>0 then one or more values failed tolerance test.  Rearrange
	!things and go do another iteration


! ______________BOTH_____________________________________________
!     OUTPUT DERIVATIVES FOR THIS FINITE DIFFERENCE ITERATION
	if(iLong(2).eq.1)then
		call WriteOutJacobian()
		endif

! ----------------------------
      	call setTPX		! this routine normalizes compositions - fix so it doesn't normalize the grain boundary
!     	WRITE TO TERMINAL AND DISK
      	IF(iLong(4).eq.1)then
         	write(12,*)
         	write(12,*)
         	write(12,*)' Newton counter = ',newtonCounter
         	write(12,*)'Calculated Deltas'
         	write(12,*)'                      .       Value  .This increment.   Running sum'
         	do i = 1,nvar
         		write(12,4060)I,VN1(I),VN2(I),ALLX(1,i),ALLX(1,i)-ALLX(5,i),ALLX(1,i)-ALLX(2,i)
4060     		FORMAT(I3,1X,A2,A16,F15.4,3E15.5)
			end do
		endif

	if(iLong(2)+iLong(3)+iLong(4)+iLong(5)+iLong(6).ge.1)then
		write(*,*) "Type 0 (enter) to continue with next step"
		write(*,*) "Type 1 (enter) to abort"
		write(*,*) "Type -1 (enter) to continue without extended output"
!		pause "Hit return to continue with next step"
		read(*,*)izero
		if(izero.eq.1)then
			go to 9999
		    	endif
		if(izero.lt.0)then
			iLong(2) = 0
			iLong(3) = 0
			iLong(4) = 0
			iLong(5) = 0
			iLong(6) = 0
			endif
		endif

!     DONE WITH THIS STEP
      	if(inewton.eq.1.and.newtonStepsOutput.eq.1)then
		call NewtonStepSave
		write(*,*)newtonCounter,imiss
		endif

      	if(inewton.eq.1.and.imiss.gt.0)then
 		newtoncounter = newtoncounter+1
		if(newtoncounter.gt.1000)then
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			write(*,*)'Newton Counter > 1000 iterations. Aborting....'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			izero = 1
			go to 9999
			endif		
        	go to 100	! loop for another iteration
		endif
	
9999	continue		! exit routine

      	RETURN
      	END


! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Compute4EQUIL_Vac(izero)
	use MatrixArrays
      implicit none

! 	Routine to calculate Jacobian for equilibrium model phases
!	This routine is a generalized solution where
!	(a) Reactions are stored in array ARX
!		EQ phases have an entry for every phase component

!	This routine requires that the GB phase has a vacancy as the first phase component.
!	The composition of the vacancy is a component that isn't used in any other phase.
!		I chose CO2 because it ensures the CO2 mass balance will be the last equation
!	We replace the CO2 mass balance with the oxygen balance equation:
!		0 = ∆Mph1 - ∆MPh2 		(if 2 solid phases)
!		0 = ∆MPh1 + ∆MPh2 - ∆MPh3	(if 3 solid phases) 
!		This works as Ox balance because all phases are written with only 1 oxygen

!	Since the GB is the first phase and Vac is the first component, Xvac will be calculated by differenct in routine Sub SetTPX
!		This allows the vacancy to vary.

! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "Singular.inc"
!	include "Solute.inc"
!	include "Tangent.inc"
	include "MatrixSingular.inc"
! c*****************************************
!     local variables
      integer*4 i,ii,j,jj,k,L,LL,Lcont,jcol,irow,ier,j1,kCur,			&
     	 izero,imiss,nsolve,newtoncounter,ivol,iRowLast
      REAL*8 DETA,R,SUM,Pi,TEMPT,TEMPP,TEMPX,TEMPH,TEMPS,TEMPK,xTemp(20),allXTemp(100)
      CHARACTER vn51(100)*8
      DATA R,Pi/8.3144D0,3.141592654D0/


! --------------------------------------------

!	Set the variables so they are understandable
	! nVar = number of variables
	! NX = number of independent compositional variables (= np - nph because 1 variable is dependent for each phase)
	! numEqRxn = number of equilibrium reactions
	! numEqEquations = total number of equations = numEqRxn + nc (equilibrium reactions plus mass balance)
	! variance = nVar - nunEqEquations
	! numNonMonits = numEqEquns 


      IF(iLong(2).eq.1)then

         write(12,*)' '
         write(12,*)'------------------------------------------------------'
         write(12,*)'------------------------------------------------------'
         write(12,*)' Variable list '
         do i = 1,nvar
        	write(12,15)I,VN1(I),VN2(I)
15      	FORMAT(I3,2X,A2,A16)
		end do
         write(12,*)'numEqEquns=',numEqEquns,'  NVAR=',NVAR
         write(12,*)
         write(12,*)'VARIANCE =',NVAR-numEqEquns
         write(12,*)
         write(12,*)'MONITOR PARAMETERS    Delta X for monitors'
         if(nvar-numEqEquns.gt.0)then
         do i = 1,nvar-numEqEquns
         	write(12,16)MON(I),VN1(MON(I)),VN2(MON(I)),Deltax(i)*NSTEP
16       	FORMAT(I3,2X,A2,A16,F15.8)
		end do
         write(12,*)
	write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NumEqEquns)

	!write(12,*)' DELTA Xi FOR MONITOR PHASES'
	!write(12,16)(DELTAX(I)*NSTEP,I=1,NVAR-numEqEquns)
	!16       FORMAT(10F12.3)
         else
         write(12,*)' No monitor parameters (variance <= 0)'
         endif
         write(12,*)'NUMBER OF FINITE DIFF. ITERATIONS=',NSTEP
         write(12,*)
         write(12,*)'INITIAL MINERAL DATA'
         write(12,17)TC,PB
17       format(' T START (DEG C)=',F10.2,'   P START=',F10.2)
	if(PFluidSwitch.gt.0)write(12,18)PFluid
18	format(' PFluid = ',F10.2)
	endif
90    CONTINUE
!     LCONT IS COUNTER IN CONTOUR BLOCK
      LCONT=0
      newtoncounter=0

! -----------------------------------------------------------
!                             
100   CONTINUE


      TK=TC+273.15D0

!     THESE CALLS CALCULATE THE THERMODYNAMIC
!     PROPERTIES OF THE PHASES OF INTEREST

!     SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      	izero=0
      	ivol=0            !compute everything
      	K=0
	! x       write(*,*)' In c2ompute #3-calling duTPX'
	! x       pause
	!	This call will callculate all thermodynamic properties - the TP dependent ones and the composition dependent ones
	!	We need this when we're using T or P as dependent variables. But not when they're independent (because they're already set)
      	CALL ALLKduTPX(ivol,izero)
	! x       write(*,*)' In c2ompute #4-back from duTPX'
	! x        pause
	if(izero.gt.0)go to 9999

	!COMPUTE PHASE VOLUMES
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
            	VP0(K)=MP0(K)*VMOL(K)*10.D0
		end do

	!OUTPUT THERMO DATA TO DISK

      IF(iLong(3).eq.1)then
350   	CONTINUE
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
          	write(12,*)
          	write(12,200)MINREC(K),SITMUL(K),PHNAME(K)
200       	FORMAT(' ',I4,F8.1,4X,A32)
          	write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
          	IF(IMASS.EQ.1) then
          		write(12,507)MP0(K),VP0(K)
          		endif
507       	FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
          	write(12,501)(phCoName(k,j),J=1,numPhCo(K))
501       	FORMAT('          ',12(A4,6X))
          	write(12,504)(xPhCo(k,j),J=1,numPhCo(K))
504       	FORMAT('  COMP ',12(F10.4))
          	write(12,505)(Dexp(lnAct(k,j)),J=1,numPhCo(K))
505       	FORMAT(' Activ  ',12(E10.3))
          	write(12,506)(lnAct(k,j),J=1,numPhCo(K))
506       	FORMAT('  Ln(a) ',12(F10.5))
          	write(12,514)(hPhCoZero(k,j),J=1,numPhCo(K))
514       	Format(' hPhCoZero',12E15.7)
          	write(12,515)(HATTP(k,j),J=1,numPhCo(K))
515       	Format(' HATTP',12E15.7)
          	write(12,509)(sPhCoZero(k,j),J=1,numPhCo(K))
509       	Format(' SZERO',12F15.5)
          	write(12,510)(SATTP(k,j),J=1,numPhCo(K))
510       	Format(' SATTP',12F15.5)
          	write(12,508)(vPhCoZero(k,j),J=1,numPhCo(K))
508       	Format(' VZERO',12F15.5)
          	write(12,513)(VATTP(k,j),J=1,numPhCo(K))
513       	Format(' VATTP',12F15.5)
          	write(12,516)(GATTP(k,j),J=1,numPhCo(K))
516       	Format(' GATTP',12E15.7)
          	write(12,*)'dudTPX ARRAY:'
          	write(12,*)'     dudT   .    dudP   .    dudX2  .    dudX3  .    dudX4...'
           	do j = 1,numPhCo(K)
           		write(12,511)dudTPX(k,j,1),dudTPX(k,j,2),(dudTPX(k,j,L),L=3,numPhCo(K)+1)
511       		FORMAT(15E12.5)
			end do
408       	format(' ',a5,6E10.3)
		end do

	    endif


!     THIS PART OF THE PROGRAM SETS UP THE MASTER MATRIX

!     ZERO MASTER ARRAY
      	DO J=1,NVAR
		DO I=1,numEqEquns
			AA(I,J)=0.0D0
			end do
		end do
!     The partial derivatives in matrix dudTPX are:
!     In storage space 1: dudTPX(jj,1):
!     Temperature:      /T = - S(P,T,X)  -- Partial molar Entropy at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state entropy, plus the entropy of mixing plus any exess entropy
!                                              terms.
!     In storage space 2: dudTPX(jj,2):
!     Pressure   :      /P =   V(P,T,X)  -- Partial molar Volume at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state volume, plus any exess volume
!                                              terms.
!     In storage space 3 to (numPhCo(K)+1): dudTPX(jj,3), dudTPX(jj,3), dudTPX(jj,3) etc:
!           Note that there are only numPhCo(K)-1 partial derivatives and the subscripts range
!           from 3 to numPhCo(K)+1.
!     Composition:      /X = (ideal)/X + (margules)/X + (reciprocal)/X

!       
!     The matrix Arx(numEqRxn,npc) contains the stoichiometric coefficients for each phase component in the EQ reactions
!	Note that these equations define the phase components for the minerals in terms of the phase components of the grain boundary (matrix)
!     in each of the reactions, subscript 1 refers to reaction number, subscript 2 refers to the phase component

!     The matrix AA therefore contains in each cell terms such as
!     
!       (n(i)*u(i))/T
!       (n(i)*u(i))/P
!       (n(i)*u(i))/X2
!       (n(i)*u(i))/X3

!     for each independent variable (where n(j,i) are the stoichiometric coefficients for phase component i in reaction j)
!       Each row of AA contains coefficients for each of the reactions.

!     The only catch to the calculation of AA is that the matrix dudPTX is
!     not stored in true matrix form.  That is, to save space the independent variables are listed for each phase only.
!     This simply requires a bit of extra bookkeeping during the matrix multiplication

	ii = np - numEqRxn
	do i = 1,numEqRxn		! loop for each reaction
		TempT = 0.0d0
		TempP = 0.0d0
		j = 0
		do kcur = 1,numPh		! Loop for every phase
			k = asmCurrent(kCur)
			do jj = 1,numPhCo(k)	! Loop for every phase component in every phase
				j = j+1
				! dT and dP coefficient
				TempT = TempT + Arx(i+ii,j)*(dudTPX(k,jj,1))	! index 1 is partial molar S
				TempP = TempP + Arx(i+ii,j)*(dudTPX(k,jj,2))	! index 2 is partial molar V
				end do
!1015				continue
			end do
!1014			continue
	      	AA(i,1)=TEMPT
	      	AA(i,2)=TEMPP
		end do
!1010		continue

!	Do dX terms
	do i = 1,numEqRxn
		L = 2				! starts in the third column (2 = dT, dP)
		jj = 0
		do kcur = 1,numPh
			k = asmCurrent(kCur)
			if(numPhCo(K).eq.1)then
				jj = jj + 1
				go to 1021	! there are no dX terms for a phase of fixed composition
				endif
			Do LL = 3,numPhCo(K)+1     	! this loops on the number of independent composition derivatives
				tempX = 0.D0
				do j=1,numPhCo(K)          	! this loops on the number of phase components in each K phase
					TEMPX=TEMPX + Arx(i+ii,jj+j)*dudTPX(k,j,LL)
					end do
				L=L+1
				AA(i,L)=TEMPX
				end do
			jj = jj + numPhCo(k)
1021  			continue
			end do
		end do




!     	SET UP MASS BALANCE EQUATIONS

!     	FIRST SET UP dM Columns (M=MOLE FRACTION OF PHASE)
!     	JCOL COUNTS WHERE dM TERMS ARE
      	JCOL = TANDP + NX
!     	loop through each phase (there is a dM term for each phase)
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
	      	JCOL = JCOL + 1
		!LOOP THROUGH ALL mass balance equations for each SYSTEM COMPONENT
      		call bulkcomp				! get total moles of each system component (moles(i) in common block)
      		DO i = 1,nc			! loop through every mass balance equation
      			iRow = i + numEqRxn			! iRow counts the row number for the mass balance equations
      			sum=0.D0
			!LOOP THROUGH ALL PHASE COMPONENTS for this phase
      			DO j=1,numPhCo(K)
				sum = sum + comp(k,J,I)*xPhCo(k,J)
				end do
				!1120  	continue
	
		      	!AA(IROW,JCOL)=SUM
		      	AA(IROW,JCOL)=SUM/numOxygens(k)
		      	YY(IROW) = -(moles(i) - molesStart(i))			! molesStart(i) is the bulk composition we want to fit
      										! moles(i) is the current bulk comp based on phase comp and M(k)
										! when they are equal, we have the correct X and M for the desired bulk comp
										! openFlux are open system components either dependent or independent variables
			end do
		end do


!     	NOW SET UP COEFFICIENTS FOR dX TERMS in mass balance equations
!     	Note that the dependent dX term is removed as above
      	JCOL=TANDP					!     	JCOL is the column for the dX term
!     	FIRST LOOP THROUGH ALL PHASES TO FIND independent dX TERMS
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
		IF(numPhCo(K).EQ.1) GO TO 1130
		!IF HERE THEN WE HAVE FOUND A PHASE WITH dX TERMS (I.E. MORE THAN 1 COMPONENT)
		!LOOP THROUGH ALL dX TERMS FOR THIS PHASE
			DO J=2,numPhCo(K)
				JCOL=JCOL+1
				!LOOP THROUGH ALL MASS BALANCE EQUATIONS FOR EACH dX TERM
				DO I=1,NC
					!IROW= I + NRX
					IROW= I + numEqRxn
					!AA(IROW,JCOL)=MP0(K)*(COMP(k,J,I)-COMP(k,1,I))
					AA(IROW,JCOL)=MP0(K)*(COMP(k,J,I)-COMP(k,1,I))/numOxygens(k)
					end do
					!1150	continue
				end do
				!1140  	CONTINUE
1130  		CONTINUE
		end do


	! end of mass balance equations
	iRowLast = numEqRxn + nc

	! now substitute the oxygen balance equation for the last mass balance equation
!		0 = ∆Mph1 - ∆MPh2 		(if 2 solid phases)
!		0 = ∆MPh1 + ∆MPh2 - ∆MPh3	(if 3 solid phases) 
	do i = 1,nVar
		AA(iRowLast,i) = 0.0d0
		end do
	YY(iRowLast) = 0.0d0
	! Find the MPhase variables
	jCol = TandP + nx + 1	! This should be the first ∆M variable = MGB
	select case (numPh)	! should be 3 or 4
	case(3)
		AA(iRowLast,jCol+1) =  1.0d0 
		AA(iRowLast,jCol+2) = -1.0d0 
	case(4)
		AA(iRowLast,jCol+1) =  1.0d0 
		AA(iRowLast,jCol+2) =  1.0d0 
		AA(iRowLast,jCol+3) = -1.0d0 
	case default
		call FSS_Alert('Alert','Number of phases is not 3 or 4')
	end select
	
!     If we are in MODE=NEWTON then we need to set up a data vector that
!     contains -(H - TS + R T lnKeq) for each reaction
!     (note there is no V term here because H and S are calculated at T and P
!		ii = np-nrx      ! note that reaction coefficients are stored in the end of array ARX
	DO i=1,numEqRxn
		TEMPH=0.D0
		TEMPS=0.D0
		TEMPK=0.D0
		jj = 0
		Do kCur = 1,numPh
			k = asmCurrent(kCur)
			do J = 1,numPhCo(k)
				jj = jj + 1
				TEMPH =TEMPH + Arx(i+ii,jj)*HatTP(k,j)		! ∆H reaction
				TEMPS =TEMPS + Arx(i+ii,jj)*SatTP(k,j)		! ∆H reaction
				TEMPK =TEMPK + Arx(i+ii,jj)*lnAct(k,j)		! ∆H reaction
				end do
			end do
			!1415    	continue      
		YY(i)=-(TEMPH - TK * TEMPS + R*TK*TEMPK)
		end do


!     	Rearrange MASTER matrix AA
!     	MOVE MONITOR PARAMETERS TO RIGHT SIDE OF EQUNS IN ARRAY A
	nsolve=nvar
	if(inewton.eq.1)nsolve=nvar+1
	DO I=1,NVAR-numEqEquns
		DO J=1,numEqEquns
		  	A(J,numEqEquns+I)=-AA(J,MON(I))
			end do
		end do

!     	FILL LEFT SIDE OF EQUATIONS IN MATRIX A FROM AA
      	DO I=1,numEqEquns			!This should be the number of non-monitor parameters
 	     	DO J=1,numEqEquns	! number of equations
			A(J,I)=AA(J,IPOINT(I))
			end do
		end do

!       Fill Nvar+1=nsolve array element in A with the Y vector      
	DO J=1,numEqEquns
		A(J,Nsolve)=YY(J)
		end do	

            

!     	Output matrix and modified matrix, if desired
	if(iLong(5).eq.1)then
2330  		CONTINUE
      		write(12,*)' '
      		write(12,*)' '
      		write(12,*)' MASTER MATRIX'
	        write(12,2321)(VN1(j),VN2(j),J=1,nvar),'  YY(i)'
2321  		FORMAT('     ',30(A2,A4,6x))
	      	DO I=1,numEqEquns
			!write(12,*)'--------------------'
			write(12,2320)(AA(I,J),J=1,NVAR),YY(I)
			end do
			!2310  		continue
2320  		FORMAT(' ',30E12.4)
		do i = 1,nvar-numEqEquns
			vn51(numEqEquns+i) = trim(vn1(mon(i)))//vn2(mon(i))
			end do
			!2311		continue
		do i = 1,numEqEquns
			vn51(i) = trim(vn1(ipoint(i)))//vn2(ipoint(i))
			end do
			!2312		continue
		write(12,*)' '
		write(12,*)' '
		write(12,*)' MODIFIED MASTER MATRIX'
	        write(12,2313)(VN51(j),J=1,nvar),' YY(i)'
2313  		FORMAT('     ',30(A6,6x))
		DO I=1,numEqEquns
			!write(12,*)'--------------------'
			write(12,2320)(A(I,J),J=1,nsolve)
			end do
			!2315  		continue
		endif

!     Find solution
      	DETA=0.D0
      	IER=0
      	CALL REDUCE (numEqEquns,NSOLVE,DETA,IER)
      	IF (IER.EQ.1) then
		izero=1
!		if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)go to 2509
		WRITE(*,*)' ************ ERROR **************************'
		WRITE(*,*)' Matrix failed to invert in SUBROUTINE REDUCE'
		WRITE(*,*)' The most likely cause is that the variables'
		WRITE(*,*)' specified as monitor parameters are not linearly'
		WRITE(*,*)' independent.  Note that only the number of'
		WRITE(*,*)' intensive variables equal to the phase rule'
		WRITE(*,*)' variance of the assemblage are independent.'
		WRITE(*,*)' Return to previous menu and try again.'
		write(*,*)' '
		pause 'Hit return to continue...'
		go to 9999
		endif

2509	continue


	if(matrixIsSingular.eq.1)go to 9999

	if(iLong(2).eq.1)then
      		write(12,*)'   '
      		write(12,*)' RESULTS OF THIS FINITE DIFFERENCE ITERATION'
      		write(12,*)'DETERMINANT=',DETA
		endif


! -------------------------------
!     	Set ALLX(5,L) equal to current values of ALLX(1,L)
!     	Note that below ALLX(1,L) is incremented by DELX
      	DO L=1,NVAR
		ALLX(5,L)=ALLX(1,L)
		end do
		!3000  	ALLX(5,L)=ALLX(1,L)

! ____________________Newtons method___________________________
      	if(inewton.eq.1)then
!     		XX IS A numEqEqunsx(NVAR-numEqEquns+1) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		The last column (NVAR-numEqEquns+1) is the solution to the x for Newtons method
!     		Calculate the change in the variables
!     		New variables are computed as
!     		Xnew = Xold + X
! 		Note that we are doing these calculations with the values of the "independent" variables constant
! 		i.e. dxi = 0. For example, if it is univariant and we are doing calculations at constant P then dP = 0
! 		This way, we only need to solve for the values of the "dependent" variables (dT, dXi, etc.)
!

		! skip this newton step stuff -- it seems to be causing trouble
!		go to 5010

!       	Composition variables can only range from 0-1. 
!		If the solution vector component (xx) is too large, it may result in the change in mole fraction (xPhCo(i) exceeding this range
!		This block of code checks this and adjusts newtonStep to ensure this doesn't happen
!		newtonStep = 1.0d0

		newtonStep = newtonStepMax
		jCol = NVAR-numEqEquns+1	! index of solution vector in array XX
		do i = 1,nvar
			allXtemp(i) = allX(1,i)	! load up allXTemp
			end do
			!5002		continue
! ----		loop to here for newtonStep iteration --------
5001		continue
		do i = 1,numEqEquns
			ii=ipoint(i)
			!adjust the values of AllXtemp that we have just calculated
			!note that this changes all dependent values, not just the composition ones
			!we could make this a wee bit faster if we just calculated the composition ones but we're not sure which of them are dependent
			!jcol is the column where the solution vector resides in XX
			ALLXtemp(ii)=ALLX(1,ii) + newtonstep*XX(i,jCol)
			end do
			!5004		continue
!		now pick out the composition values and store in xTemp
		l = TandP		! this is where the independent composition derivatives start
		Do kCur = 1,numPh
			k = asmCurrent(kCur)
			if(numPhCo(k).eq.1)go to 5010	! no composition derivatives for 1 component phases
			sum = 1.0d0		! this is for the dependent composition derivative (the first one)
			DO J = 2,numPhCo(K) ! loop from second to last phase component for this mineral
				l = l + 1		! index of independent composition variable in AllX
				!j indexes the phase components of phase K
				xTemp(j) = allXTemp(l)
				sum = sum - xTemp(j)
				end do
				!5020  		continue
			xTemp(1) = sum
			ier = 0
			call CheckStoichiometry(k,xTemp,ier)
			if(ier.eq.1)then			! failed stoichiometry test
				newtonStep = newtonStep/2.0d0
				!if(newtonStep.gt.1.0d-6)go to 5001			! try again with this smaller stepsize
				if(newtonStep.gt.1.0d-15)go to 5001			! try again with this smaller stepsize
				!stepsize is too small. Abort this phase
				!if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)then
				izero = 1
				write(*,*)'stepSize = ',newtonStep
				write(*,*)'This is too small. Something must be wrong in the code. Line 497 of Compute4EQUIL'
				write(*,*)'Phase = ',k
				write(*,*)phName(k)
				do j1 = 1,numPhCo(K)
					write(*,*)phCoName(k,j1),xPhCo(k,j1),xTemp(j1)
					end do
					!5026			continue
				pause ' hit return to continue'
				izero = 1		! set error flag
				return
				endif

5010  			continue
			end do
		!write(*,*)'NewtonStep = ',newtonStep
		!check for convergence
		!for convergence we check to see whether the x value is < tolerance*X value      
      		jj = NVAR-numEqEquns+1		! index of solution vector in array XX
      		imiss=0
      		do i = 1,numEqEquns
	      		ii=ipoint(i)
			!Newtonstep scales the solution vector so we can step towards the solution in smaller steps
	      		ALLX(1,ii)=ALLX(1,ii) + newtonstep*XX(i,jj)
			!AllX(1,ii) is the current value of the dependent variables. xx(i,jj) is the solution vector (variable).
			!When variable is less than a fraction(TolNewton) of the dependent variable, we have converged
	      		if(dabs(xx(i,jj)).gt.dabs(TolNewton*Allx(1,ii)))then
	      			if(allFractl(ii).eq.0)then		! only count imiss if the phase is not fractionating (seeGibbs_FileInput.f line 520 or so)
	      		      		imiss=imiss+1
					endif
				endif
			end do
			!5110  		continue
 		if(iopen.eq.1)then
			do i = 1,iopen
				openFlux(iOpena(i)) = allX(1,openPtr - 1 + i)
				end do
				!5120			continue
 			endif	! end open system check
		!if imiss>0 then one or more values failed tolerance test.  Rearrange
		!things and go do another iteration
      		go to 3200
      		endif	! end iNewton = 1
! ________________________________________________________________
      
! _______________________Gibbs Method__________________________
      	if(inewton.eq.0)then
!     		XX IS A numEqEqunsx(NVAR-numEqEquns) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		COMPUTE new values FOR EACH VARIABLE
!     		note that the old variables DELX and SDELX are now defined as
!           	DELX(i)  = ALLX(1,i) - ALLX(5,i)  :current minus last finite difference
!           	SDELX(i) = ALLX(1,i) - ALLX(2,i)  :current minus starting
!     		Definition of ALLX variable:
!     		1,50   current value
!     		2,50     starting value
!     		3,50   ref at start of contour
!     		4,50   user selected reference
!     		5,50   previous finite diff point
!     		6,50     (not used)
!     		Store new values for each independent variable (monitor parameter)
      		DO J=1,NVAR-numEqEquns
			ii=mon(j)
			ALLX(1,ii)=ALLX(1,ii) + DELTAX(J)
			end do
			!3100  		continue
!     		compute new values for each dependent variable
      		Do i = 1,numEqEquns
	      		ii=ipoint(i)
	      		DO J=1,NVAR-numEqEquns
		      		ALLX(1,ii)=ALLX(1,ii) + XX(i,J) * DELTAX(J)
				end do
			end do
			!3110  		continue    
      		endif       ! gibbs/newton

! ________________________________________________________________
3200  	continue

! ______________BOTH_____________________________________________
!     OUTPUT DERIVATIVES FOR THIS FINITE DIFFERENCE ITERATION
	if(iLong(2).eq.1)then
		call WriteOutJacobian()
	   	endif

! ----------------------------
      	call setTPX
!     	WRITE TO TERMINAL AND DISK
      	IF(iLong(4).eq.1)then
         	write(12,*)
         	write(12,*)
         	write(12,*)'Calculated Deltas'
         	write(12,*)'                      .       Value  .This increment.   Running sum'
         	do i = 1,nvar
         		write(12,4060)I,VN1(I),VN2(I),ALLX(1,i),ALLX(1,i)-ALLX(5,i),ALLX(1,i)-ALLX(2,i)
4060     		FORMAT(I3,1X,A2,A16,F15.4,3E15.5)
!4060     		FORMAT(I3,1X,A2,A16,3F15.4)
			end do
			!4018     	continue
		endif

	if(iLong(2)+iLong(3)+iLong(4)+iLong(5)+iLong(6).ge.1)then
		write(*,*) "Type 0 (enter) to continue with next step"
		write(*,*) "Type 1 (enter) to abort"
		write(*,*) "Type -1 (enter) to continue without extended output"
!		pause "Hit return to continue with next step"
		read(*,*)izero
		if(izero.eq.1)then
			go to 9999
		    	endif
		if(izero.lt.0)then
			iLong(2) = 0
			iLong(3) = 0
			iLong(4) = 0
			iLong(5) = 0
			iLong(6) = 0
			endif
		endif

!     DONE WITH THIS STEP
      	if(inewton.eq.1.and.newtonStepsOutput.eq.1)then
		call NewtonStepSave
		write(*,*)newtonCounter,imiss
		endif

      	if(inewton.eq.1.and.imiss.gt.0)then
 		newtoncounter = newtoncounter+1
		if(newtoncounter.gt.1000)then
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			write(*,*)'Newton Counter > 1000 iterations. Aborting....'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			izero = 1
			go to 9999
			endif		
        	go to 100	! loop for another iteration
		endif
	
9999	continue		! exit routine

      	RETURN
      	END


