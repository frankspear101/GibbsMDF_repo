! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine Variance
      implicit none
!     Subroutine to calculate the thermodynamic variance
! ****************************************
 	include "Assemb.inc"
! ****************************************

      	gibbsVariance = TANDP + NX - NRX
      	duhemsVariance = nvar - neq
      	MDFVariance = nvar - neq + 2
		      write(*,*)' Gibbs variance  = ',gibbsVariance
      	if(imass.eq.1)Write(*,*)' Duhems variance = ',duhemsVariance
		      write(*,*)' MDF   variance  = ',MDFVariance

	return
	
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine ZeroiLong
      implicit none
!     Subroutine to zero the iLong array (output options
	include "Output.inc"
      	ilong(1)=0
      	ilong(2)=0
      	ilong(3)=0
      	ilong(4)=0
      	ilong(5)=0
      	ilong(6)=0
		iLong(7)=0
		iLong(8)=0
		iLong(9)=0
		iLong(10)=0
		ilong(11)=0
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine FluidPresent
      implicit none
!     Subroutine to determine if one of the minerals is a fluid
! 	returns value of KFlu (in Assemb.inc)
! 	KFlu = 0 if no fluid is present
! 	KFlu = mineral number in list if a  fluid is present
! ****************************************
	include "Assemb.inc"
! ****************************************
	integer*4 K,kCur
!--------------------------------------------------------
! 	check to see if a fluid is present from the fluid numbers list
! 	Fluid numbers in data base are
! 	2	H2O only (WHAAR equation)
! 	3	Mixed H2O-CO2 fluid (Kerrick and Jacobs w/WHAAR for H2O)
! 	4	Mixed H2O-CO2 fluid (not working yet)

!       Check to see if a fluid is really in data file
	Do 10 kCur = 1,numPh
	k = asmCurrent(kCur)
	IF(MINREC(K).EQ.2.or.Minrec(K).eq.3.or.Minrec(K).eq.4.or.Minrec(K).eq.2002)GO TO 20
10	continue
! 	A fluid is not present
 	! call FSS_alert('ALERT!!','A fluid is not present in the data file.')
	KFlu = 0
	return

! 	if here, then a fluid is present.  Go to appropriate subroutine to calculate volume
20	continue
	KFlu = K
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine PhaseComp
      implicit none
!     Subroutine to determine the phase composition 
! ****************************************
	include "Assemb.inc"
! ****************************************

!     local variables
      integer*4 i,j,k,kCur
      real*4 sum
! *******************************************************************
!     Figure out the compositions of the minerals
	Do 431 kCur = 1,numPh
	k = asmCurrent(kCur)
      do 432 i = 1,NC
      PhComp(K,i) = 0.
      do 430 j=1,numPhCo(K)
      PhComp(K,i) = PhComp(K,i) + xPhCo(k,j)*comp(k,j,i)
430   continue
432   continue
431   continue

	Do 150 kCur = 1,numPh
	k = asmCurrent(kCur)
!     convert moles to wt % oxides
      sum=0.D0
      do 120 i=1,nc
      PhWtOx(K,i)=PhComp(K,i)*molwt(i)/NumCatInOxide(i)
120   sum=sum + PhWtOx(k,i)
      if(sum.le.0.)go to 140
       sum=100.D0/sum
      do 130 i=1,nc     
 	PhWtOx(k,i)=PhWtOx(k,i)*sum
	PhWtEl(k,i)=PhWtOx(k,i)*OxToElWt(i)
130	continue
140   continue
150	continue
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine CheckMassBal
      implicit none
!     Subroutine to check that mass balance is still within tolerance
! ****************************************
	include "Assemb.inc"
	include "Output.inc"
! ****************************************

!     local variables
      integer*4 i,j,k,kCur
!      real*8 pcterror(24)
! *******************************************************************
!     Compute number of moles of each system component
	do 100 i=1,nc
	moles(i)=0.D0
	Do 110 kCur = 1,numPh
	k = asmCurrent(kCur)
	do 115 j=1,numPhCo(K)
	moles(i) = moles(i) + MP0(K)*xPhCo(k,j)*comp(k,j,i)
115	continue
110	continue
! 	Check difference from starting value
	go to 100
!	pctError(i) = Dabs(   (moles(i)-molesStart(i))/molesStart(i)  )
!	if(pctError(i).gt.0.1)then
!		write(12,*)' element ',coname(i),moles(i),pctError(i)
!		if(masserroralert.eq.0)then
!	call fss_alert('Mass balance error is greater than 10%.Try using smaller step sizes. (This message will not appear again.)')
!			masserroralert = 1
!			endif
!		endif
100	continue
	return
	END



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine BULKCOMP ()
      implicit none
!     Subroutine to determine the bulk composition of a rock given the 
!     mineral compositions and modes
! ****************************************
	include "Assemb.inc"
! ****************************************

!     local variables
      integer*4 i,j,k,kCur
      real*8 sum
! *******************************************************************
!     Compute number of moles of each system component
      do 100 i=1,nc
      wtpct(i)=0.D0
      moles(i)=0.D0
	Do 100 kCur = 1,numPh
	k = asmCurrent(kCur)
      do 100 j=1,numPhCo(K)
100   moles(i) = moles(i) + MP0(K)*xPhCo(k,j)*comp(k,j,i)
!     convert moles to wt % oxides
      sum=0.D0
      do 120 i=1,nc
      wtpct(i)=moles(i)*molwt(i)/NumCatInOxide(i)
120   sum=sum + wtpct(i)
      if(sum.le.0.)go to 140
      sum=100.D0/sum
! 	sum = 1.D0		! don't normalize
      do 130 i=1,nc     
130   wtpct(i)=wtpct(i)*sum
140   continue
	return
      END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE names
      implicit none
!     THIS SUBROUTINE IS DESIGNED TO SET UP variable names

! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
! ****************************************

!     Local variables
      integer*4 i,j,k,L,kCur
      CHARACTER VNAME(15)*2,BLANK*4
      common /local3/vname,blank

      BLANK='    '
!      VNAME(1)='T'
!      vname(2)='P'
!      vname(3)='u'
!      vname(4)='X'
!      vname(5)='M'
!      vname(6)='m'
!      vname(7)='1'
!      vname(8)='2'
!      vname(9)='3'

      VNAME(1)='T_'
      vname(2)='P_'
      vname(3)='u_'
      vname(4)='X_'
      vname(5)='M_'
      vname(6)='m_'
      vname(7)='1_'
      vname(8)='2_'
      vname(9)='3_'

!     THIS SECTION SETS UP THE VARIABLE NAMES

      VN1(1)=VNAME(1)
      VN2(1)=BLANK
      VN1(2)=VNAME(2)
      VN2(2)=BLANK
      IF(KFLU.NE.0)THEN
         VN1(3)=VNAME(2)
         VN2(3)=PHNAME(KFLU)
      ENDIF
5006  L=TANDP
	Do 5004 kCur = 1,numPh
	k = asmCurrent(kCur)
            DO 5005 J=2,numPhCo(K)
             L=L+1
             VN1(L)=VNAME(4)
5005         VN2(L)=phCoName(k,J)
5004  CONTINUE
      IF (IMASS.EQ.1) then
	 molesPhPtr = L + 1
	Do 5001 kCur = 1,numPh
	k = asmCurrent(kCur)
         L=L+1
         VN1(L)=VNAME(5)
5001     VN2(L)=PHNAME(K)
         ENDIF
!     names of open system components
      if(iopen.ne.0) then
	    openPtr = L + 1
            do 5100 i=1,iopen
            L=L+1
            vn1(L)=vname(6)
            vn2(L)=coname(iOpena(i))
5100        continue
      endif
!     NAMES  of new variables
      newVarPtr = L + 1			! pointer to new variables in AllX array - KEEP!
!	do 20 kCur = 1,numPh
!	k = asmCurrent(kCur)
!	if(includeNew(k).eq.1)then
!		do 21 i = 1,numNewInPhK(k)
	!      do 20 I = 1,numNew
!	        L=L+1
!       	 jj = Jnew(i)
!       	  VN1(L)=' '
!        	VN1(L)='__'
!        	VN2(L)=newVarName(k,i)
!       	 VN2(L)=NewName(Jnew(i))
!21		continue
!		endif
!20    continue
! 	now place dashes in all unused names, so they don't show up by accident
	do 30 i = L+1,50
	VN1(i) = '--'
	vn2(i) = '----'
30	continue
      return
      end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE RESET(IRESET,ISUB)
    	use AWE_Interfaces
	use MyCanvas
      implicit none
!     SUBROUTINE TO RESET THE MINERAL COMPOSITIONS TO THEIR STARTING
!     VALUES
! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
! ****************************************

!     local variables
      integer*4 j,k,L,izero,ireset,isub,kCur
      real*4 XXX,YYY
      character*128 ppp
!
!
!
!     VALUES FOR IRESET
!        2 = STARTING CONDITIONS
!        3 = Reference point for Contour routine
!        4 = REFERENCE POINT-user selected
!        5 = PREVIOUS finite difference POINT
!        6 = IN sub CONTOUR ( go to 1/2 of last finite
!             difference step)

!                  allx(6,L) not used

!     Values of ISUB specify from where this subroutine was called
!     ISUB = 1  called from STEP
!     ISUB = 2  called from CONTOUR
!     ISUB = 3  called from TERNARY
!     ISUB = 4  called from CONTOUR
!
!
! --------------------------------------------------
!
	DO 1490 L=1,nvar
	IF(IRESET.EQ.6)THEN
		ALLX(1,L)=(ALLX(1,L)+ALLX(5,L))/2
		ELSE
		ALLX(1,L)=ALLX(IRESET,L)
		ALLX(5,L)=ALLX(1,L)     !set previous finite difference point equal to current point
		ENDIF
1490	CONTINUE
	L=TANDP
	Do 25 kCur = 1,numPh
	k = asmCurrent(kCur)
	IF(numPhCo(K).EQ.1)GO TO 424
	xPhCo(k,1)=1
	DO 24 J=2,numPhCo(K)
	L=L+1
	xPhCo(k,J)=ALLX(1,L)
	xPhCo(k,1)=xPhCo(k,1)-xPhCo(k,J)
24	CONTINUE
424	CONTINUE
25	CONTINUE
	!     RESET MINERAL ABUNDANCE MOLES
	IF(IMASS.EQ.0) GO TO 30
	L=TANDP+NX
	Do 35 kCur = 1,numPh
	k = asmCurrent(kCur)
	L=L+1
	MP0(K)=ALLX(1,L)
	MP1(K)=MP0(K)
	MP2(K)=MP0(K)
35	CONTINUE
30	CONTINUE
!
	TC=ALLX(1,1)
	PB=ALLX(1,2)
	IF(KFLU.NE.0)PFLUID=ALLX(1,3)

1000	continue
	if(isub.eq.3)go to 1020
	IF(IRESET.NE.6)THEN
	!     calls to plot only to move plotter pen back to starting conditions
		xxx=ALLX(1,NXPLT)
		YYY=ALLX(1,NYPLT)
		IF (NYPLT.EQ.2)yyy=YYY/1000.
		call symb(XYPlot,XXX,YYY,1.,1.)    !draw a symbol on the screen just to show we got there
		ENDIF
1020		continue

	if(imass.eq.1.and.ireset.ne.6)then
	!        PHASE VOLUMES RECOMPUTED HERE
		izero=0
		call AllKduTPX(1,izero)
	!        returns in J/bar; display in cc/mole
	!        moles per cm**3 (note VMOL is in joule/bar =CC/10)
	!        note: if this code is changed, it must also be changed in routine CHANGE
		Do 300 kCur = 1,numPh
		k = asmCurrent(kCur)
		vp0(k)=mp0(K)*vmol(k)*10.D0
		vp1(k)=vp0(k)
		vp2(k)=vp0(k)
		if(FRACTL(k).eq.1.OR.FRACTL(k).eq.2)then
			ppp = 'Warning: The phase '//Trim(phname(K))//' is fractionating. Volumes and moles may not be reset correctly.'
			endif
300		continue
		masserroralert = 0  	!switch for massbalance error alert (reset to start conditions)
		endif

	RETURN
	END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine SavRxn
      
!     Subroutine to save to a file a digitized reaction curve
      implicit none
! ****************************************
	include "TPSave.inc"
! ****************************************


      integer*4 iunit,iopt,i
      character*32 flname
      character*32 fspearhead
      character*64 rxnname
      character*32 prompt

      integer*4 iopen
      data iopen /0/
      data iunit /76/

      data fspearhead /'FSpear digitized data file      '/
! ************main menu options******************************

	save
      if(iopen.eq.0)then      !a data file is not yet opened
1           continue
            write(*,*)'Data file is not yet opened.  Make a choice'
            write(*,*)'0 = return'
            write(*,*)'1 = open new file'
            write(*,*)'2 = open old file (for append)'
	    read(*,*)(iopt)
            if(iopt.eq.0)return

            if(iopt.eq.1)then
                prompt = 'Input new file name'
                flname = 'T-P output.dig'
		open(iunit,File='',status='NEW')
                write(iunit,*)FSpearhead
                iopen=1
 
 	        elseif(iopt.eq.2)then
      	        OPEN(IUNIT,FILE='',STATUS='OLD')
                do
                read(iunit,*,end=9)
                repeat
 9              continue
                Backspace(iunit)
                iopen=1

                else
                go to 1  !only if user typed in bad option number
                endif

          endif

      write(*,*)'Input a title for this reaction'
      read(*,'(A64)')Rxnname
      write(iunit,'(A64)')Rxnname
      do 100 i = 1,TPcounter
      write(iunit,1000)i,TPSave(i,1),TPSave(i,2)
100   continue
1000  format(I5,2f12.2)
      write(iunit,1000)0,0,0

      return
      end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SetTPX
      implicit none
!
! 	This routine was part of old Subroutine Compute and was moved to
! 	here when Newton's method was introduced.  
! 	The purpose is to take the new values of variables as stored in ALLX
! 	and place them into the correct compositional arrays, P and T.
!
! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
 	include "Newton.inc"
! ****************************************
!     local variables
      integer*4 j,k,L,kCur
      real*8 DELX,SUM


	! This routine sets the arrays after calculations in Compute
	! We need to
	! Set new T and P (if not isothermao and isobaric
	! Set compositions in array xPhCo and calculate the dependent variable assuming sumX = 1
	!  ... or, for grain boundary phase use the reaction stoichiometries
	! Set the phase mole arrays (MP0 etc).
	! Set NEW P AND T
	TC=ALLX(1,1)
	TK=TC+273.15D0
	PB=ALLX(1,2)
	PFluid = PB
	!set new values of mineral composition in phase component array X
	L=TANDP                 ! Start counting ALLX array at first independent component
	Do kCur = 1,numPh
	!Do kCur = 2,numPh	! only do the solid phases
		k = asmCurrent(kCur)
		sum=0.0d0
		DO J = 2,numPhCo(K) ! loop from second to last phase component for this mineral
			L=L+1
			xPhCo(k,j)=ALLX(1,L)
			sum=sum+xPhCo(k,j)
			end do
		xPhCo(k,1)=1.0D0-sum		! only set the dependent variable for solid phases so sum = 1
		end do

	!call CalculateGBdependent()

	! Code for the grain boundary -- this is NOT correct
	! GB phase
! 	xPhCo(k,1) = GBdependentX
! 	write(12,*)'GBdependentX = ',GBdependentX


!     compute new moles of phases
	! imass is always = 1
!      IF(IMASS.EQ.0)GO TO 3620
!        COMPUTE NEW Mi FOR EACH PHASE
	L=TANDP+NX
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
		L=L+1
		DELX  =ALLX(1,L) - ALLX(5,L)
		MP0(K)=MP0(K)+DELX
		mp1(k)=mp1(k)+delx
		mp2(k)=mp2(k)+delx
		VP0(K)=MP0(K)*VMOL(K)*10.D0   !note that molar volume of phase (VMOL) is in J/mol
		vp1(k)=mp1(k)*vmol(k)*10.D0 ! and VP0, VP1, and VP2 are in cm^3
		vp2(k)=mp2(k)*vmol(k)*10.D0
		end do    


      RETURN
      END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE CalculateGBdependent()
      implicit none

! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
 	include "Newton.inc"
	include "Output.inc"
! ****************************************
!     local variables
      integer*4 i,j,k,L,kCur
      real*8 DELX,SUM,temp,delXdep


	!calculate the value for the dependent GB phase component (usually H) based
	! on the amount of each phase component (xPhCo) and the amount of the dependent X in each phase
	! and on the amount of moles of each phase that is produced or consumed
	IF(iLong(4).eq.1)then
		write(12,*)'---Calc dependent----'
		endif
	delXdep = 0.0d0
	L = TANDP + NX + 1		! should be index just before the first solid phase
	Do kCur = 2,numPh	! do all phases except the GB phase
		k = asmCurrent(kCur)
		i = nc			! the last system component = the first GB component

		sum=0.D0
		!LOOP THROUGH ALL PHASE COMPONENTS for this phase
		DO j=1,numPhCo(K)
			sum = sum + comp(k,j,i)*xPhCo(k,j)
			end do
		! now "sum" contains the moles of sysCo(i) in phase k
		! Get the ∆moles of phase k
		L=L+1
		DELX  =ALLX(1,L) - ALLX(5,L)
		temp = -delX*sum
		delXdep = delXdep + temp	
		IF(iLong(4).eq.1)then
			write(12,*)kcur,k,sum,delX,temp,delXdep
			endif
		end do
	xPhCo(1,1) = XPhCo(1,1) + delXdep
	IF(iLong(4).eq.1)then
		write(12,*)' Xdep = ',xPhCo(1,1)
		endif
	

	return
	end
		

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Solve_for_H(ThmoFile)
	use MatrixArrays
      implicit none
!
!     written July, 1994 by F. Spear
!     ROUTINE TO calculate enthalpies of phase components from an input dataset
!
! ****************************************
	include "Assemb.inc"
	include "Output.inc"
!	include "Solute.inc"
! ****************************************
!     local variables
      integer*4 i,ii,j,jj,k,ier,izero,ivol,numSolve,kCur,kSolve(20),jSolve(20),kk
      REAL*8 DETA,R,TEMPH,TEMPS,TEMPK
      Character*32 ThmoFile         !thermodynamic data file name
!
      DATA R/8.3144D0/
!

          call printt(1)        !option 1 prints to the file and screen

!         RECOMPUTE REACTION MATRIX

	  iLong(7) = 1
          CALL REXN
          iLong(7) = 0





!
      TK=TC+273.15D0
!
!     THESE CALLS CALCULATE THE THERMODYNAMIC
!     PROPERTIES OF THE PHASES OF INTEREST
!
!     SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      izero=0
      ivol=0            !compute everything
      K=0

      CALL AllKduTPX(ivol,izero)
      if(izero.gt.0)return
!
!     OUTPUT THERMO DATA TO DISK
!


!      iout=9
350       CONTINUE
 	Do 390 kCur = 1,numPh
	k = asmCurrent(kCur)
          write(12,*)
          write(12,200)MINREC(K),SITMUL(K),PHNAME(K)
200       FORMAT(' ',I2,F8.1,4X,A32)
          write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
          IF(IMASS.EQ.1)then
          write(12,507)MP0(K),VP0(K)
          endif
507       FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
          write(12,501)(phCoName(k,J),J=1,numPhCo(K))
501       FORMAT('          ',6(A4,6X))
          write(12,504)(xPhCo(k,J),J=1,numPhCo(K))
504       FORMAT('  COMP ',6(F10.4))
          write(12,505)(Dexp(lnAct(k,J)),J=1,numPhCo(K))
505       FORMAT(' Activ  ',6(E10.3))
          write(12,506)(lnAct(k,J),J=1,numPhCo(K))
506       FORMAT('  Ln(a) ',6(F10.5))
          write(12,514)(hPhCoZero(k,J),J=1,numPhCo(K))
514       Format(' hPhCoZero',6E15.7)
          write(12,515)(HATTP(k,J),J=1,numPhCo(K))
515       Format(' HATTP',6E15.7)
          write(12,509)(sPhCoZero(k,J),J=1,numPhCo(K))
509       Format(' SZERO',6F15.5)
          write(12,510)(SATTP(k,J),J=1,numPhCo(K))
510       Format(' SATTP',6F15.5)
          write(12,508)(vPhCoZero(k,J),J=1,numPhCo(K))
508       Format(' VZERO',6F15.5)
          write(12,513)(VATTP(k,J),J=1,numPhCo(K))
513       Format(' VATTP',6F15.5)
          write(12,516)(GATTP(k,J),J=1,numPhCo(K))
516       Format(' GATTP',6E15.7)


511       FORMAT(12E12.5)
512	continue
!         write(12,408)' ACP',(aPhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' BCP',(bPhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' CCP',(cPhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' DCP',(dPhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' ECP',(ePhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' FCP',(fPhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' GCP',(gPhCoCp(k,J),J=1,numPhCo(K))
408       format(' ',a5,6E10.3)
551       CONTINUE
390       CONTINUE
!
!

!
!     THIS PART OF THE PROGRAM SETS UP THE MATRIX
!
!     ZERO MASTER ARRAY
!
      DO 1000 J=1,NVAR
      DO 1000 I=1,NEQ
1000  AA(I,J)=0.0D0
!
!
!     Finished with extra equations
!
!
!     Set up a data vector (Y) that contains H - TS + R T lnKeq for each reaction
!

      	ii = np-nrx      !note that reaction coefficients are stored in the end of array ARX
      	DO 1010 i=1,NRX
      	TEMPH=0.D0
      	TEMPS=0.D0
      	TEMPK=0.D0
      	jj = 0
	Do 1015 kCur = 1,numPh
	k = asmCurrent(kCur)
      	do 1015 J = 1,numPhCo(k)
!     	S reaction
 	jj = jj + 1
      	TEMPS =TEMPS + ARX(i+ii,jj)*SatTP(k,j)
!     	LnK reaction
      	TEMPK = TEMPK + ARX(i+ii,jj)*lnAct(k,j)
!     	H reaction
      	TEMPH =TEMPH + ARX(i+ii,jj)*HatTP(k,j)
1015  	continue      
      	YY(i)=(TEMPH - TK * TEMPS + R*TK*TEMPK)
1010	continue

!     Now choose the phase components to solve for H

      WRITE(*,*)'Phase components in this problem'
!	ier=0
	Do 105 kCur = 1,numPh
	k = asmCurrent(kCur)
      Do 105 j=1,numPhCo(k)      
!	ier = ier + 1
      WRITE(*,*)k,j,phCoName(k,j)
105   continue

      write(*,*)' Choose phase components to solve for H'
      write(*,*)' You must choose this many components: ',nrx
      do 106 i=1,nrx
!      write(*,*)'Component ',i
!	read(*,*)j
!      Nsolve(i)=j
!      write(*,*)'Component ',i,j
      write(*,*)'Phase-K, Component-J (2 integers)'
      read(*,*)kSolve(i),jSolve(i)
!	read(*,*)j
!      Nsolve(i)=j
      write(*,*)'Phase, Component ',kSolve(i),jSolve(i)
106   continue

!     Now construct a matrix (AA) with the appropriate coefficients for each phase component to solve.

!	pause 'this code will not work.....solve for H in Gibbs_misc.f'
      numSolve = nrx+1
      do 1105 i=1,nrx
      TempH = 0
 !	The problem is that this loop is for the pointer arrays, which no longer exist.
 !	I need to correctly index the right phase component.
	jj = 0
	do 1110 kCur = 1,numPh
	k = asmCurrent(kcur)
	do 1112 j = 1,numPhCo(k)
	jj = jj + 1
	do 1120 kk = 1,nrx
	if(kSolve(kk).eq.k.and.jSolve(kk).eq.j)then
	        A(i,kk) = -ARX(ii+i,jj)       !the A matrix coefficient is the reaction coefficient
	        TempH = TempH + ARX(ii+i,jj)*hPhCoZero(k,j)    !Subtract hPhCoZero from data vector because we need to solve for this
                                                !Note it was added in above already
		endif
!      do 1110 j = 1,np
!      do 1120 k = 1,nrx
!      if(nsolve(k).eq.j)then        !this is one of the components to solve for H
!        A(i,k) = -ARX(ii+i,j)       !the A matrix coefficient is the reaction coefficient
!        TempH = TempH + ARX(ii+i,j)*hPhCoZero(k,j)    !Subtract hPhCoZero from data vector because we need to solve for this
                                                !Note it was added in above already
!        endif
1120  	continue
1112	continue
1110  	continue
      	A(i,numSolve) = YY(i) - TempH      
1105  	continue
            
!
!     Output matrix and modified matrix, if desired
!
1320  FORMAT(' ',10E16.4)
!
      write(12,*)' '
      write(12,*)' '
      write(12,*)' MASTER MATRIX'
      write(12,1321)(phCoName(kSolve(i),jSolve(i)),i=1,nrx)
1321  format(6x,10A16)
      DO 1312 I=1,NRX
      write(12,*)'--------------------'
      write(12,1320)(A(I,J),J=1,numSolve),YY(i)
1312  continue
!
!     Find solution
!
      DETA=0.D0
      IER=0
      CALL REDUCE (NRX,numSolve,DETA,IER)
      IF (IER.EQ.1) then

      WRITE(*,*)' ************ ERROR **************************'
      WRITE(*,*)' Matrix failed to invert in SUBROUTINE REDUCE'
      WRITE(*,*)' The most likely cause is that the variables'
      WRITE(*,*)' specified as monitor parameters are not linearly'
      WRITE(*,*)' independent.  Note that only the number of'
      WRITE(*,*)' intensive variables equal to the phase rule'
      WRITE(*,*)' variance of the assemblage are independent.'
      WRITE(*,*)' Return to previous menu and try again.'
      write(*,*)' '
      izero=1
      pause 'Hit return to continue...'
      return
      endif
!
!
      write(12,*)'   '
      write(12,*)' RESULTS '
      write(12,*)'DETERMINANT=',DETA

!     Write out solution vector
      write(12,*)
      write(12,*)'Phase components        H'
      Do 1605 i=1,nrx
!      write(12,*)phCoName(k,Nsolve(i)),XX(i,1)
      write(12,*)phCoName(kSolve(i),jSolve(i)),XX(i,1)
1605  continue
      write(12,*)

      RETURN
      END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	subroutine timdat(iunit)
!	integer *4 MM,DD,YY,SEC,Hrs,Min
!        character*255 ppp


!	call date(MM,DD,YY)
!	Call TIME(sec)
!	Hrs=INT((sec+1)/3600)
!	sec=sec-3600*hrs
!	Min=int(SEC)/60
!	sec=sec-60*min
!	write(12,10)MM,DD,YY,Hrs,min,sec
!        call printf(ppp)
!10	format('Date: ',I2,'-',I2,'-',I2,'    Time: ',I2,':',I2,':',I2)
!	return
!	end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine XtoSite(K)
      implicit none
!     ROUTINE to compute site occupanicies of multisite minerals from mole fractions of components
! ****************************************
	include "Assemb.inc"
! ****************************************
!     local variables
      integer L,j,K

!     calculate atoms on site information
      if(numSiteAtom(K).eq.0)go to 699
      do 611 L = 1,numSiteAtom(K)
      SiteAtom(k,L) = 0
      do 610 j = 1,numPhCo(K)
      SiteAtom(k,L) = SiteAtom(k,L) + XToSiteAtom(k,L,j)*xPhCo(k,j)
!     Note that the difficult subscripting is due to:
!     (1)  X and phCoUsedinFile are linear arrays of only those phase components that are used
!            phCoUsedinFile contains integers refering to the original list of components from the thermodynamic data file
!     (2)  SiteAtom is a linear array of all the site fractions (used and unused) from the thermo file
!     (3)  XtoSiteAtom is a 2-D array
!           subscript 1 is the same as SiteAtom subscripts
!           subscript 2 refers to the original list of components in the thermo file, and hence needs
!                 the index from phCoUsedinFile
610   continue
611   continue

699   continue
      return
      end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine WriteOutJacobian
	use MatrixArrays
	implicit none
	include "Assemb.inc"
 	include "Monit.inc"
 	include "Newton.inc"
! 	include "Solute.inc"
	integer*4 i,j,jj,ll,nsolve
	character soln1*2,soln2*4

	write(12,*)
	write(12,3507)
3507    FORMAT('    Jacobian matrix:')
	if(inewton.eq.0)then
		write(12,3502)(VN1(MON(I)),VN2(MON(I)),i=1,nvar-neq)
		nsolve = nvar
		else
		soln1 = 'So'
		soln2 = 'ln V'
		write(12,3502)(VN1(MON(I)),VN2(MON(I)),i=1,nvar-neq),soln1,soln2
		nsolve = nvar + 1
		endif		

3502    FORMAT(' Monitors:',5x,15(A2,A4,9x))
	write(12,3506)(deltax(I),i=1,nvar-neq)
3506    format('      ',15(f15.4))
	JJ=0
	DO 3520 I=1,NVAR
	DO 3523 J=1,NVAR-NEQ
3523    IF(I.EQ.MON(J))GO TO 3520
	JJ=JJ+1
	write(12,3501)VN1(I),VN2(I),(XX(JJ,LL),ll=1,nsolve-neq)	! writes out soln in last XX
3501    FORMAT(' ',A2,A4,' = ',15E15.5)
3520    CONTINUE
	return
	end

