	real*8 tanPlane(sysCoMax)
	real*8 uPhaseAu(phCoMax),uzeroPhaseDeltaAU(phCoMax),uPhaseDeltaAU(phCoMax)
	real*8 uOnTanAU(phCoMax),uOnTanDeltaAU(phCoMax)
	real*8 gOnTanAU(phMax),gOfPhaseAU(phMax),gDifferenceAU(phMax)

	common/tangentCom/	  					&
     		tanPlane, 						&
     		uPhaseAu,uZeroPhaseDeltaAU,uPhaseDeltaAU,		&
     		uOnTanAU,uOnTanDeltaAU,					&
     		gOnTanAU,gOfPhaseAU,gDifferenceAU
	! variables with "au" at the end are in atom units. All others in molar units


