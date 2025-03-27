	integer*4 numPairs,Pairs(24,3),numTrios,Trios(24,4)
	! Subscripts for GBComp
	!	24 = Pair number
	!	2  = start, finish (GB), phase 2 or 3 (moles)
	!	10 = number of components in GB
	real*8 GBComp(2,10),MDFmoles(2,4),MDFAff(2,4)
	common /AutoAff/numPairs,Pairs,GBComp,MDFmoles,MDFAff,numTrios,Trios	