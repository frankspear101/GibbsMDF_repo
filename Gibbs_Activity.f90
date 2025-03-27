! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      Subroutine Activity(K,uMin,Xtemp,TKtemp,PBtemp,R,iflag)
!     Activity of a generic multi site mineral, as defined in data file
! 	Subroutine returns R*T*ln(act) of components in uMin
      implicit none
! ****************************************
	include "Assemb.inc"
	include "Output.inc"
! ****************************************
      	integer*4 K,iflag,j,L,n,i1,i2,j1,k1,jj,ii
      	real*8 uMin(PhCoMax),Xtemp(PhCoMax),TKtemp,Pbtemp,R,RTK,lnatemp,uexcess,sum,pkzero1,pkzero2,TR
	real*8 SiteFraction(numSiteAtomMax),XtempASF(PhCoMax),tempsum,ASFtotal(PhCoMax),WGHP11(50)
	real*8 WGrecip,uRecip,uRecipTemp

	real*8 MargulesWG(numMargulesWTermsMax,numMargulesWTermsMax,numMargulesWTermsMax),XMarg(numMargulesSiteCatsMax)
	real*8 uexcessBerman(numSiteAtomMax),uexcessTemp(numMargulesSiteCatsMax),uexcessMargules
	integer*4 WIndex(numMargulesWTermsMax,numMargulesWTermsMax,numMargulesWTermsMax),index,index1
      	data TR/298.15D0/

	if(iLong(6).eq.1)then
		write(12,*)
		write(12,*)'TK = ',TKtemp,'    Pbar = ',Pbtemp
		write(12,*)'Min ',Minrec(K),phname(K)
		endif

	RTK = R*TKtemp


!-------------------------
!	This section sets the appropriate composition terms into the array Xtem2 for use in calculating the excess energy
	select case (DatasetKey(k))
	case (1)		! Spac
!		Ideal activity calculations
	!      	Calculate mole fractions of cations on sites
	!     	in terms of mole fractions of thermodynamic components
	! 	This is the same code as in subroutine XtoSite except here I am
	! 	using Xtemp for compositions instead of X because I need to calculate
	! 	the finite difference derivative w/r/t X and Xtemp is passed through from calling routine
		do 1011 L = 1,numSiteAtom(K)
		SiteFraction(L) = 0
		do 1010 j = 1,numPhCo(K)
		SiteFraction(L) = SiteFraction(L) + XToSiteAtom(k,L,j)*Xtemp(j)/SiteMultiplicity(k,L)
1010   		continue
1011   		continue
! 		Check to see that no atoms in use are less than zero
		do 1015 L = 1,numSiteAtom(K)
        	SiteFraction(L) = SiteFraction(L)+1.d-20	! if zero, make very slightly larger
 		if(SiteFraction(L).le.1.0D-30)go to 999
1015		continue
		if(iLong(6).eq.1)then
			write(12,*)' Component X values'
			write(12,1016)(Xtemp(j),j=1,numPhCo(K))
1016			format(15E15.8)
			write(12,*)' Site Fractions '
			write(12,*)(SiteFraction(L),L=1,numSiteAtom(K))
			endif
! 		If here, then all site atoms are > 0.0

!		Calculate Margules excess values for every site
		if(numMargulesSites(k).gt.0)then
			do 1205 jj = 1,numSiteAtom(k)
			uexcessBerman(jj) = 0.0d0
1205			continue

			do 1210 jj = 1,numMargulesSites(k)
	
			do 1212 i1 = 1,numMargulesSiteCats(k,jj)
			uexcessBerman(i1) = 0.0d0
			do 1212 j1 = 1,numMargulesSiteCats(k,jj)
			do 1212 k1 = 1,numMargulesSiteCats(k,jj)
			MargulesWG(i1,j1,k1) = 0.0D0
			Windex(i1,j1,k1) = 0				! integer
1212			continue		
			do 1215 L = 1,numMargulesSiteCats(k,jj)
			Xmarg(L) = SiteFraction(MargulesSiteCats(k,jj,L))
			uexcessTemp(L) = 0.0d0
1215			continue
			!write(*,*)' Margules calcs. Phase = ',k,phName(k)
			do 1220 L = 1,numMargulesWterms(k,jj)
			index = MargulesWindex(k,jj,L)
			index1 = index
!			write(*,*)L,index
			i1 = index1/100
			index1 = index1 - i1*100
			j1 = index1/10
			index1 = index1 - j1*10
			k1 = index1/1
!			MargulesWG(i1,j1,k1) = MargulesWH(k,jj,L) - MargulesWS(k,jj,L)*(TKtemp - TR) + MargulesWV(k,jj,L)*PbTemp
			MargulesWG(i1,j1,k1) = MargulesWH(k,jj,L) - MargulesWS(k,jj,L)*(TKtemp) + MargulesWV(k,jj,L)*PbTemp
			WIndex(i1,j1,k1) = index
			!write(12,*)index,MargulesWG(i1,j1,k1)
1220			continue
			!pause ' take a look and hit return'
			uexcessTemp = 0.0d0
			i1 = numMargulesSiteCats(k,jj)
			j1 = numMargulesWterms(k,jj)
!			,MargulesWG,WIndex,Xmarg,uexcessBerman(jj)
!			call BermanExcess(numMargulesSiteCats(k,jj),numMargulesWterms(k,jj),MargulesWG,WIndex,Xmarg,uexcessBerman(jj))
			call BermanExcess(i1,j1,MargulesWG,WIndex,Xmarg,uexcessTemp)
			do 1225 i1 = 1,numMargulesSiteCats(k,jj)
			!MargulesSiteCats(phMax,numMargulesSitesMax,numMargulesSiteCatsMax)
			uexcessBerman(MargulesSiteCats(k,jj,i1)) = uexcessTemp(i1)
1225			continue
			!write(12,*)'Phase,site,uexcess = ',k,jj,uexcessBerman(jj)
1210			continue
			if(iLong(6).eq.1)then
				write(12,*)'uexcessBerman '
				write(12,1256)(SiteAtomName(k,jj),jj=1,numSiteAtom(k))
1256				Format(20(7x,A8))
				write(12,1257)(uexcessBerman(jj),jj=1,numSiteAtom(k))
1257				Format(20F15.3)
				endif

			endif		! end Margules calculations


!		****************************** phase component loop
! 		loop for every phase component
		Do 1020 j = 1,numPhCo(K)	! calculate  for each phase component
!		Calculate ideal activities
! 		If ActivityConstant = 0, then we are using a molecular model
		if(ActivityConstant(k,j).eq.0.0)then
			Lnatemp = Dlog(Xtemp(j))
			else
! 			Ideal activity part is made from expression
! 			R*TK*(ln(ActivityConstant) + alpha1*ln(X1) + alpha2*ln(X2) + alpha3*ln(X3) + .....)
			Lnatemp = Dlog(ActivityConstant(k,j))
			do 1030 L = 1,numSiteAtom(K)
			lnatemp = lnatemp + Alpha(k,j,L)*Dlog(SiteFraction(L))
1030			continue
			endif
	

!	Reciprocal free energy calculations
			uRecip = 0.0d0
			do  1110 ii = 1,numReciprocal(k)
			WGrecip = WHrecip(k,ii) - WSrecip(k,ii)*(TKtemp-TR) + WVrecip(k,ii)*Pbtemp
			do 1130 L = 1,numSiteAtom(K)
			if(RecipIndex(k,j,ii,L).eq.0)go to 1130
			if(RecipIndex(k,j,ii,L).eq.1)then		! either +1 or -1
				uRecipTemp = WGrecip * SiteFraction(L)
				else
				uRecipTemp = WGrecip * (1.0d0 - SiteFraction(L))
				endif				
1130			continue
			uRecipTemp = uRecipTemp * float(RecipConstant(k,j,ii))		! either +1 or -1
			uRecip = uRecip + uRecipTemp	
1110			continue


!		Add in Margules contribution
		uexcessMargules = 0.0d0
		do 1040 jj = 1,numMargulesSites(k)
		do 1046 i1 = 1,numMargulesSiteCats(k,jj)
		L = MargulesSiteCats(k,jj,i1)				! this is the cation number
		if(XtoSiteAtom(k,L,j).ne.0)then				! L is the site atom
			uexcessMargules = uexcessMargules + uexcessBerman(L)
			go to 1046					! loop out because we only add uexess once for each site
			endif
1046		continue
1045		continue
1040		continue



1250		continue		! jump here for melt models
		uMin(j) = lnatemp*RTK + uexcessMargules + uRecip

		if(iLong(6).eq.1)then
!			write(12,*)'Ph comp,lnatemp,lnatemp*RTK,uRecip,uexcess,utotal'
			write(12,1055)j,lnatemp,lnatemp*RTK,uRecip,uexcessMargules,uMin(j)
1055			Format('Component,lnatemp,lnatemp*RTK,uRecip,uexcess,utotal',I5,8E15.5)
			endif


1020		continue		! end of loop for this phase component

		iflag = 0
		return


	case (2)	! HP98 activity models

	!      	Calculate mole fractions of cations on sites
	!     	in terms of mole fractions of thermodynamic components
	! 	This is the same code as in subroutine XtoSite except here I am
	! 	using Xtemp for compositions instead of X because I need to calculate
	! 	the finite difference derivative w/r/t X and Xtemp is passed through from calling routine
		do 2011 L = 1,numSiteAtom(K)
		SiteFraction(L) = 0
		do 2010 j = 1,numPhCo(K)
		SiteFraction(L) = SiteFraction(L) + XToSiteAtom(k,L,j)*Xtemp(j)/SiteMultiplicity(k,L)
2010   		continue
2011   		continue
! 		Check to see that no atoms in use are less than zero
		do 2015 L = 1,numSiteAtom(K)
        	SiteFraction(L) = SiteFraction(L)+1.d-20	! if zero, make very slightly larger
 		if(SiteFraction(L).le.1.0D-30)go to 999
2015		continue
! 		If here, then all site atoms are > 0.0


		do 2030 n = 1,numSFterms(k)
!		WGHP11(n) = SFWH(k,n) - SFWS(k,n)*TKtemp + SFWV(k,n)*Pbtemp/1000.0d0	Wv is KJ/Kb or J/bar -- all my calcs are in J/bar
		WGHP11(n) = SFWH(k,n) - SFWS(k,n)*(TKtemp-TR) + SFWV(k,n)*Pbtemp
2030		continue

		if(iLong(6).eq.1)then
			write(12,*)'Min ',Minrec(K),' Nonideal terms = ',numSFterms(k)
			write(12,*)'  WH   WS    WV    ID1   ID2    WG'
			do 2045 n = 1,numSFterms(k)
			write(12,3043)SFWH(k,n),SFWS(k,n),SFWV(k,n),SFID(k,n,1),SFID(k,n,2),WGHP11(n)
2043			format(3F12.3,2I8,F15.3)
2045			continue
			endif

! 		loop for every phase component in the current problem
		Do 2020 j = 1,numPhCo(K)	! calculate  for each phase component
		select case (PhaseType(k))
		case(105)		! HP11 melt model ideal activities
!			This code REQUIRES that the first phase component in the melt is always H2O
!			ASF parameters for HP11 melt are all = 1, so Xtemp(j) = XtempASF(j)
!			For the melt the site fraction is equal to the mole fraction of the component
			if(j.eq.1)then
				lnatemp = Xtemp(1)**2			! this one for H2OL
				else
				lnatemp = Xtemp(j)*(1.0d0-Xtemp(1))	! Xtemp2(1) is H2OL
				endif
				lnatemp = Dlog(lnatemp)
		case(103)		! HP11 ideal ionic activities
!			Calculate ideal activities
! 			If ActivityConstant = 0, then we are using a molecular model
			if(ActivityConstant(k,j).eq.0.0)then
				Lnatemp = Dlog(Xtemp(j))
				else
! 				Ideal activity part is made from expression
!	 			R*TK*(ln(ActivityConstant) + alpha1*ln(X1) + alpha2*ln(X2) + alpha3*ln(X3) + .....)
				Lnatemp = Dlog(ActivityConstant(k,j))
				do 2021 L = 1,numSiteAtom(K)
				lnatemp = lnatemp + Alpha(k,j,L)*Dlog(SiteFraction(L))
2021				continue
				endif
		case default
			write(*,*)' PhaseType error in thermodyanmic data file for phase:'
			write(*,*)k,phName(k)
			write(*,*)' PhaseType must be >200 for HP11 files. PhaseType for this phase = ',PhaseType(k)
			write(*,*)' This must be fixed and the program restarted'
			pause 'Hit return to continue'
		end select
	
	
		uexcess = 0.0d0
		! Equation 20 of Powell and Holland 1993 Am Mineral
		! and Holland and Powell 2003

!		----- We are calculating the uexcess for phase component j
!		Note that for the set of independent components, pkzero is either 0 or 1
!		However, for dependent components it may be a fraction 
!			(pkzero is the mole fraction of the dependent component in the component of interest)
!		But since Gibbs only deals with independent component sets, this coding should be fine
		tempsum = 0.0d0
		do 2050 n = 1,numSFterms(k)
		i1 = SFID(k,n,1)
		pkzero1 = 0.0d0
		if(j.eq.i1)pkzero1 = 1.0d0
		i2 = SFID(k,n,2)
		pkzero2 = 0.0d0
		if(j.eq.i2)pkzero2 = 1.0d0
!		tempsum = tempsum + WGHP11(n)*ASFtotal(j)*(pkzero1 - XtempASF(i1))*(pkzero2 - XtempASF(i2))
		tempsum = tempsum + WGHP11(n)*(pkzero1 - Xtemp(i1))*(pkzero2 - Xtemp(i2))
2050		continue
		uexcess = - tempsum   

2250		continue		! jump here for melt models
		uMin(j) = lnatemp*RTK + uexcess

		if(iLong(6).eq.1)then
			write(12,2055)j,exp(lnatemp),lnatemp,lnatemp*RTK,uexcess,uMin(j)
2055			Format('Component,act,lnatemp,lnatemp*RTK,uexcess,utotal',I5,5E15.5)
			endif

2020	continue	! end of loop for every phase component


		iflag = 0
		return



	case(3)		! for HoPo11  we need the size-paramater adjusted mole fractions of components
! 		This code is for excess models that use symmetric formalism (moles of components rather than cations on sites)

	!      	Calculate mole fractions of cations on sites
	!     	in terms of mole fractions of thermodynamic components
	! 	This is the same code as in subroutine XtoSite except here I am
	! 	using Xtemp for compositions instead of X because I need to calculate
	! 	the finite difference derivative w/r/t X and Xtemp is passed through from calling routine
		do 3011 L = 1,numSiteAtom(K)
		SiteFraction(L) = 0
		do 3010 j = 1,numPhCo(K)
		SiteFraction(L) = SiteFraction(L) + XToSiteAtom(k,L,j)*Xtemp(j)/SiteMultiplicity(k,L)
3010   		continue
3011   		continue
! 		Check to see that no atoms in use are less than zero
		do 3015 L = 1,numSiteAtom(K)
        	SiteFraction(L) = SiteFraction(L)+1.d-20	! if zero, make very slightly larger
 		if(SiteFraction(L).le.1.0D-30)go to 999
3015		continue
! 		If here, then all site atoms are > 0.0

		sum = 0.0d0
		do 3025 j = 1,numPhCo(k)	! Loop through only the phase components being used
!		is this + or - ASF(2)*TK??  It's an entropy, so I think it should be -
!		ASFtotal(j) = (ASF(k,j,1) - ASF(k,j,2)*TKtemp + ASF(k,j,3)*Pbtemp)
		ASFtotal(j) = (ASF(k,j,1) - ASF(k,j,2)*(TKtemp-TR) + ASF(k,j,3)*Pbtemp)
!		ASFtotal(j) = (ASF(k,j,1) + ASF(k,j,2)*TKtemp + ASF(k,j,3)*Pbtemp/1000.0d0)
		sum = sum + Xtemp(j)*ASFtotal(j)		! Scale the composition to the ASF size
3025		continue
!		Now scale according to the size parameter
		do 3027 j = 1,numPhCo(k)
		XtempASF(j) = Xtemp(j)*ASFtotal(j)/sum
3027		continue
		if(iLong(6).eq.1)then
			Write(12,*)' ASF parameters (1,2,3, ASFtotal, Xtemp, XtempASFScaled)'
			do 3023 j = 1,numPhCo(k)
			write(12,2024)ASF(k,j,1),ASF(k,j,2),ASF(k,j,3),ASFtotal(j),xtemp(j),XtempASF(j)
2024			format(3F10.3,3F15.7)
3023			continue
			endif
!		At this point the array XtempASF(j) contains the size-parameter scaled compositions of each phase component


		do 3030 n = 1,numSFterms(k)
!		WGHP11(n) = SFWH(k,n) - SFWS(k,n)*TKtemp + SFWV(k,n)*Pbtemp/1000.0d0	Wv is KJ/Kb or J/bar -- all my calcs are in J/bar
		WGHP11(n) = SFWH(k,n) - SFWS(k,n)*(TKtemp-TR) + SFWV(k,n)*Pbtemp
3030		continue
!		Apply size parameter scaling to W terms Wij* = Wij/(ASF(i)+ASF(j))	= page 493 of HP 2003
! 			note that we will multiply by 2*ASF(phase component of interest) below
		do 3035 n = 1,numSFterms(k)
!		WGHP11(n) = WGHP11(n) /(ASFtotal(SFID(k,n,1)) + ASFtotal(SFID(k,n,2)))
		WGHP11(n) = 2.0d0*WGHP11(n) /(ASFtotal(SFID(k,n,1)) + ASFtotal(SFID(k,n,2)))
3035		continue

		if(iLong(6).eq.1)then
			write(12,*)'Min ',Minrec(K),' Nonideal terms = ',numSFterms(k)
			write(12,*)'  WH   WS    WV    ID1   ID2    WG'
			do 3045 n = 1,numSFterms(k)
			write(12,3043)SFWH(k,n),SFWS(k,n),SFWV(k,n),SFID(k,n,1),SFID(k,n,2),WGHP11(n)
3043			format(3F12.3,2I8,F15.3)
3045			continue
			endif

! 		loop for every phase component in the current problem
		Do 3020 j = 1,numPhCo(K)	! calculate  for each phase component
		! I need to change this toXS model type so new melts can be added in the data file rather than the code here
		select case (PhaseType(k))
		case(205)		! HP11 melt model ideal activities
!			This code REQUIRES that the first phase component in the melt is always H2O
!			ASF parameters for HP11 melt are all = 1, so Xtemp(j) = XtempASF(j)
!			For the melt the site fraction is equal to the mole fraction of the component
			if(j.eq.1)then
				lnatemp = XtempASF(1)**2			! this one for H2OL
				else
				lnatemp = XtempASF(j)*(1.0d0-XtempASF(1))	! Xtemp2(1) is H2OL
				endif
				lnatemp = Dlog(lnatemp)
		case(203)		! HP11 ideal ionic activities
!			Calculate ideal activities
! 			If ActivityConstant = 0, then we are using a molecular model
			if(ActivityConstant(k,j).eq.0.0)then
				Lnatemp = Dlog(Xtemp(j))
				else
! 				Ideal activity part is made from expression
!	 			R*TK*(ln(ActivityConstant) + alpha1*ln(X1) + alpha2*ln(X2) + alpha3*ln(X3) + .....)
				Lnatemp = Dlog(ActivityConstant(k,j))
				do 3021 L = 1,numSiteAtom(K)
				lnatemp = lnatemp + Alpha(k,j,L)*Dlog(SiteFraction(L))
3021				continue
				endif
		case default
			write(*,*)' PhaseType error in thermodyanmic data file for phase:'
			write(*,*)k,phName(k)
			write(*,*)' PhaseType must be >200 for HP11 files. PhaseType for this phase = ',PhaseType(k)
			write(*,*)' This must be fixed and the program restarted'
			pause 'Hit return to continue'
		end select
	
	
		uexcess = 0.0d0
!		This code is for the ASF model of
		! Equation 20 of Powell and Holland 1993 Am Mineral
		! and Holland and Powell 2003

!		----- We are calculating the uexcess for phase component j
!		Note that for the set of independent components, pkzero is either 0 or 1
!		However, for dependent components it may be a fraction 
!			(pkzero is the mole fraction of the dependent component in the component of interest)
!		But since Gibbs only deals with independent component sets, this coding should be fine
		tempsum = 0.0d0
		do 3050 n = 1,numSFterms(k)
		i1 = SFID(k,n,1)
		pkzero1 = 0.0d0
		if(j.eq.i1)pkzero1 = 1.0d0
		i2 = SFID(k,n,2)
		pkzero2 = 0.0d0
		if(j.eq.i2)pkzero2 = 1.0d0
!		tempsum = tempsum + WGHP11(n)*2.0d0*ASFtotal(j)*(pkzero1 - XtempASF(i1))*(pkzero2 - XtempASF(i2))
		tempsum = tempsum + WGHP11(n)*ASFtotal(j)*(pkzero1 - XtempASF(i1))*(pkzero2 - XtempASF(i2))
3050		continue
		uexcess = - tempsum   

3250		continue		! jump here for melt models
		uMin(j) = lnatemp*RTK + uexcess

		if(iLong(6).eq.1)then
			write(12,3055)j,exp(lnatemp),lnatemp,lnatemp*RTK,uexcess,uMin(j)
3055			Format('Component,act,lnatemp,lnatemp*RTK,uexcess,utotal',I5,5E15.5)
			endif

3020	continue	! end loop on every phase component

	iflag = 0
	return

	case default

	write(*,*)' Bad dataset type (DatasetKey). Please fix the thermodynamic data file'
	Pause 'Abort program'
	stop 

	end select


250	continue



999	continue	! a site atom has a zero or negative value
        iflag=1
	if(iLong(1).eq.1)then

	      	write(*,*)' In subroutine Activity'
		write(*,*)' A cation on a site has a zero or negative value'
		Write(*,*)PhName(K)
		write(*,*)'Atom      value'
		do 910 L = 1,numSiteAtom(K)
		write(*,*)SiteAtomName(k,l),SiteFraction(L)
910		continue
		write(*,*)
		write(*,*)'Component    value'
		do 915 j = 1,numPhCo(K)
		write(*,*)phCoName(k,j),Xtemp(j),xPhCo(k,j)
915		continue
	 	write(*,*)'Tc,Pb= ',TC,PB
	      	write(*,*) 'Hit return to continue'
	      	pause
		endif
        return


      end
