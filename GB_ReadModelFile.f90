!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine ReadModelFile()
	implicit none
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				! required for color file information
	include "Assemb.inc"
	integer*4 i,j,k,num,R,G,B
	character*128 dummy
	character*16 dum
	character*256 GridFileName
!	open(16,file='',status = 'OLD',iostat = status,err=99)
!	if(status.ne.0)then
!		call FSS_Alert('Alert','Problem opening model file')
!		return
!		endif
!        INQUIRE (16, NAME=filein)
!        logfile = trim(filein)//'.log'
!	open(66,file=logFile,status = 'NEW',iostat = status,err=99)
!	if(status.ne.0)then
!		call FSS_Alert('Alert','Problem opening log file for output')
!		return
!		endif
!	open(95,file=logFile,status = 'OLD',iostat = status,err=99)
409	format(128A)
	read(16,*)GridFileName
	call ReadGridFile(GridFileName)
	read(16,409)dummy				!dashes
	write(66,*)dummy			
	read(16,*)dXgrid,dGridRemove,dGridAdd
	write(66,*)dXgrid,dGridRemove,dGridAdd,'     dXgrid, dGridRemove,dGridAdd'
	! this doesn't work because entire grid is done with values ca 1-1000. Grid spacing must be similar (i.e. 10)
	! These are nominally microns (10-6 meters)
! 	dXgrid = dXgrid * 1.0d-6	! convert from µm to meters
! 	dgridRemove= dgridRemove * 1.0d-6	! convert from µm to meters
! 	dgridAdd = dgridAdd * 1.0d-6	! convert from µm to meters
! 	write(66,*)dXgrid,dGridRemove,dGridAdd,'     dXgrid, dGridRemove,dGridAdd'
	read(16,409)dummy				!dashes
	write(66,*)dummy			
	read(16,*)convergenceJoules			! convergence criteria in Joules
	write(66,*)convergenceJoules,'  Convergence criteria in Joules'			! convergence criteria in Joules
	read(16,409)dummy				!dashes
	write(66,*)dummy			
	read(16,409)dummy				!Title line
	write(66,*)dummy			
	read(16,*)numEl
	write(66,*)'numEl = ',numEl
	read(16,409)dummy				!Title line
	write(66,*)dummy			
	do i = 1,numEl
		read(16,*)ElName(i),ElCharge(i),Dzero(i),DActE(i),DActVol(i)
		write(66,*)ElName(i),ElCharge(i),Dzero(i),DActE(i),DActVol(i)
		DactE(i) = DactE(i)*1000.0d0		! convert from kJ to Joules
		Dzero(i) = Dzero(i) * 1.D12		! convert from m^2/sec to µm^2/sec
		write(66,*)ElName(i),ElCharge(i),Dzero(i),DActE(i),DActVol(i)
		end do
	read(16,409)dummy
	write(66,*)dummy			!dashes
	read(16,409)dummy
	write(66,*)dummy			!Title line
	read(16,409)dummy
	write(66,*)dummy			!Title line

	! Read enthalpies for each element
	! These need to be the equilibrium values calculated for the assemblage without the phase that is to nucleate
	! These are calculated for the GB composition and assemblage in the MIF using program Gibbs3
	read(16,409)dummy
	write(66,*)dummy			!dashes
	read(16,409)dummy
	write(66,*)dummy			!Title line
	do i = 1,numEl
		read(16,*)dum,hPhCoZero(1,i)
		write(66,*)dum,hPhCoZero(1,i)
		end do


	! Read colors for each element
	read(16,409)dummy
	write(66,*)dummy			!dashes
	read(16,409)dummy
	write(66,*)dummy			!Title line
	read(16,409)dummy
	write(66,*)dummy			!Title line

! 	AweNumColors = numEl + 1		! black is always the first color
	do i = 100,100+numEl			! this will store the element colors in the AWEColor arrays starting at #101
		read(16,*)dum,AWEColorName(i),AWEColorNumber(i),R,G,B,(CMYK(i,j),j=1,4)			! Note: Colors are in RGB Hexadecimal
		write(66,*)dum,AWEColorName(i),AWEColorNumber(i),R,G,B,(CMYK(i,j),j=1,4)			! Note: Colors are in RGB Hexadecimal
		end do
! 	do i = 1,numEl + 1		! black is always the first color
! 		read(16,*)dum,ElColorName(i),ElColorNumber(i),R,G,B,(ElCMYK(i,j),j=1,4)			! Note: Colors are in RGB Hexadecimal
! 		write(66,*)dum,ElColorName(i),ElColorNumber(i),R,G,B,(ElCMYK(i,j),j=1,4)			! Note: Colors are in RGB Hexadecimal
! 		end do

!	Read phases and node identities
	read(16,409)dummy			! dashes
	write(66,*)dummy			!Dashes
	read(16,409)dummy
	write(66,*)dummy			!Title line
	read(16,*)numPhases,dummy		! number of phases
	write(66,*)numPhases,dummy
	do i = 1,numPhases
		read(16,409)dummy			! dashes
		write(66,*)dummy			!Dashes
		read(16,*)PhaseMIF(i),PhaseName(i)
		write(66,*)PhaseMIF(i),PhaseName(i)
		read(16,*) PhaseColorName(i),PhaseColorNumber(i),R,G,B,(PhaseCMYK(i,j),j=1,4)	!,R,G,B,(CMYK(i,j),j=1,4)			! Note: Colors are in RGB Hexadecimal
		write(66,*)PhaseColorName(i),PhaseColorNumber(i),R,G,B,(PhaseCMYK(i,j),j=1,4)	!,R,G,B,(CMYK(i,j),j=1,4)			! Note: Colors are in RGB Hexadecimal
		read(16,*)numPhaseCrystals(i)
		write(66,*)numPhaseCrystals(i)
		if(numPhaseCrystals(i).eq.0)then		! this is probably quartz. Fill in all crystals with quarta
			do j = 1,numPolys
				PhaseCrystalNumber(i,j) = j 
				CrystalMIFID(PhaseCrystalNumber(i,j)) = PhaseMIF(i)
				CrystalPhaseName(PhaseCrystalNumber(i,j)) = PhaseName(i)
				end do
			else
			do j = 1,numPhaseCrystals(i)
				read(16,*)PhaseCrystalNumber(i,j)		! read nodes that define this phase
				write(66,*)PhaseCrystalNumber(i,j)
				CrystalMIFID(PhaseCrystalNumber(i,j)) = PhaseMIF(i)
				CrystalPhaseName(PhaseCrystalNumber(i,j)) = PhaseName(i)
				end do
			endif
		end do			

	! Read nodes to plot composition
	read(16,409)dummy
	write(66,*)dummy			!Dashes
	read(16,409)dummy
	write(66,*)dummy			!Title line
	read(16,*)numNodesToPlot
	write(66,*)numNodesToPlot
	do i = 1,numNodesToPlot
		read(16,*)nodeToPlot(i)
		write(66,*)nodeToPlot(i)
		end do

!	Read color gradients
	read(16,409)dummy
	write(66,*)dummy			!Dashes
	read(16,409)dummy
	write(66,*)dummy			!Title line
	read(16,409)dummy
	write(66,*)dummy			!Title line
	read(16,*)numPositiveColors
	write(66,*)numPositiveColors
	do i = 1,numPositiveColors
		read(16,102)positiveColors(i)
		write(66,103)positiveColors(i),positiveColors(i)
		end do
102	format(z)
103	format(z12,I12)
	read(16,409)dummy
	write(66,*)dummy			!Title line
	read(16,*)numNegativeColors
	write(66,*)numNegativeColors
	do i = 1,numNegativeColors
		read(16,102)negativeColors(i)
		write(66,103)negativeColors(i),negativeColors(i)
		end do		

!	Read grain boundary perturbation
! 	read(16,409)dummy
! 	write(66,*)dummy			!Dashes
! 	read(16,409)dummy
! 	write(66,*)dummy			!Title line
! 	read(16,*)nodeToPerturb
! 	write(66,*)nodeToPerturb
! 	read(16,409)dummy
! 	write(66,*)dummy			!Title line
! 	do i = 1,numEl
! 		read(16,*)dum,PerturbComp(i)
! 		write(66,*)dum,PerturbComp(i)
! 		end do

!	Read node capture exclusion values
!	Actually, these are the segments that connect the nodes that should never be captured
	read(16,409)dummy
	write(66,*)dummy			!Dashes
	read(16,409)dummy
	write(66,*)dummy			!Title line
	! The array segCaptureExclusion is initialized to all .FALSE. in subroutine ReadGridFile Line 220
	! this code reads in the Crystal numbers and then determines which segments to exclude
	! I changed it (May 11, 2023) so that it actually reads in the numbers of the segments that will never be flipped
! 	read(16,*)num
! 	write(66,*)num,'   Number of exclusions to read'
! 	do k = 1,num
! 		read(16,*)i	! crystal number
! 		write(66,*)i,'   Excluded crystal'
! 		do j = 1,numPolyNodes(i)
! 			segCaptureExclusion(PolySegs(i,j)) = .TRUE.	
! 			end do
! 		write(66,*)(PolySegs(i,j),j=1,numPolyNodes(i))
! 		end do

	read(16,*)num			! number of segs to be excluded
	write(66,*)num,'   Number of segment exclusions to read'
	write(66,*)'Segments to exclude'
	do i = 1,num
		read(16,*)k
		segCaptureExclusion(k) = .TRUE.
		write(66,*)i,k,segCaptureExclusion(k)
		end do




	close (16)

	call SetUpCrystals()
	call SetUpSegPoints()

	!call CalculateAtomUnits()
	
	dtime = 40.		! should ensure R < .5
        !R = Dt/x^2
        ! x = .1, D(H) = 1e-4
        ! so t = 50 at R = .5 

	return

99	continue
	call FSS_Alert('Alert','Error opening Input file...abort')	
	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine ReadGridFile(GridFileName)
	! Routine to read the grid file created from the SVG file
	implicit none
	include "Diffuse_GB_Hex.inc"
	!include "SVG_Stuff.inc"

	character*256 GridFileName
	character*24 dummy
	integer*4 status,i,j,ii,k


	open(17,file=GridFileName,status = 'OLD',iostat = status)
	if(status.ne.0)then
		call FSS_Alert('Alert','Problem opening Grid file')
		write(*,*)' GridFileName = ',GridFileName
		stop
		return
		endif

	read(17,*)dummy			! title
	read(17,*)dummy			! dashes
	read(17,*)dummy			! title (segments)
	read(17,*)numSegs
	read(17,*)dummy			! header    x1                y1                x2              y2       length
	do i = 1,numSegs
		read(17,*)ii,x1(i),y1(i),x2(i),y2(i),segLength(i)
		segGettingShorter(i) = .FALSE.				! initialize to .FALSE.
		segCaptureExclusion(i) = .FALSE.			! initialize to .FALSE.
		end do
	read(17,*)minX,maxX		!,' = minX, maxX '
	read(17,*)minY,maxY		!,' = minY, maxY '
	read(17,*)dummy			! dashes
	read(17,*)dummy			! title (Nodes)
	read(17,*)numNodes
	read(17,*)dummy			! header --- NodeX     NodeY        num       nodeSegConnect nodeSegEnd(not used)
	do i = 1,numNodes
! 		read(17,*)ii,nodeX(i),nodeY(i),numNodeSegs(i), &
! 			((nodeNodeConnect(i,j),nodeSegConnect(i,j),nodeSegEnd(i,j)),j=1,numNodeSegs(i))
		read(17,*)ii,nodeX(i),nodeY(i),numNodeSegs(i), &
			((nodeNodeConnect(i,j),nodeSegConnect(i,j),k),j=1,numNodeSegs(i))
		end do
	read(17,*)dummy			! dashes
	read(17,*)dummy			! Title   Nodes at ends of segments and segment xy
	read(17,*)dummy			! header   Seg   SegNode1   segX    segY      SegNode2     segX     segY
	do j = 1,numSegs
		read(17,*)ii,segNodes(j,1),segX(j,1),segY(j,1),segNodes(j,2),segX(j,2),segY(j,2)
		end do
	read(17,*)dummy			! dashes
	read(17,*)dummy			! Title   Polygon list'
	read(17,*)numPolys		!, ' = numPolys '
	read(17,*)dummy			! header   Poly   num    PolyNodes PolySegs
	do i = 1,numPolys
		read(17,*)ii,numPolyNodes(i),(PolyNodes(i,k),PolySegs(i,k),k=1,numPolyNodes(i))
		end do
	read(17,*)dummy			! dashes
	read(17,*)dummy			! Title   Center of Polys
	read(17,*)dummy			! header   Poly    Xcenter    Ycenter
	do i = 1,numPolys
		read(17,*)ii,PolyCenterX(i),PolyCenterY(i)
		end do
	read(17,*)dummy			! dashes

	close(17)

	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine SetUpSegPoints()
	! the segments read in from the model file have this information
! 	integer*4
! 		numSegs,				&! = number of total segments
! 		segNodes(maxSegs,2),			&! = index of the 2 nodes that connect each segment
! 		numPoints(maxSegs)			 ! = number of points along the indexed segment
! 	Real*8						&!
! 		segX(maxSegs,2),			&
! 		segY(maxSegs,2),			&!
! 		segLength(maxSegs), 			&! = length of each segment between two nodes
! 		pointX(maxSegs,maxPoints),		&!
! 		pointY(maxSegs,maxPoints),		& ! = x,y coordinates of each point along a given segment
! 		segdXgrid(maxSegs),			& ! The grid spacing for each segment
! 		dXGrid					 ! = grid spacing between points (constant throughout the entire grid)  

	! this routine creates the point arrays
	
	implicit none
	include "Diffuse_GB_Hex.inc"
	integer*4 j,i,k,iNode,numPointsT,iSeg
	real*8 dX,dY


	! Now create the diffusion grid - this one uses even spacing for all segments
	! At the end, we will scale certain crystals and adjust the number of points and grid spacing for each segment
	!	dXgrid = sideSize/10.d0		
	Do iSeg = 1,numSegs
		numPointStart(iSeg) = 1			! the initial value to start a segment
		numPointsT = int(1.0d0 + 1.0001d0*segLength(iSeg)/dXgrid)			! round up just a bit
		if(numPointsT.eq.1)numPointsT = 2						! 2 points minimum = endpoints at nodes
		numPointEnd(iSeg) = numPointsT
		segdXgrid(iSeg) = dXgrid							! here I'm using a constant grid spacing, but this could be changed
		dX = (segX(iSeg,2) - segX(iSeg,1))/(numPointsT-1)
		dY = (segY(iSeg,2) - segY(iSeg,1))/(numPointsT-1)
		pointX(iSeg,1) = segX(iSeg,1)
		pointY(iSeg,1) = segY(iSeg,1)
		do j = 2,numPointsT - 1
			pointX(iSeg,j) = segX(iSeg,1) + dX*float(j-1)
			pointY(iSeg,j) = segY(iSeg,1) + dY*float(j-1)
			end do
		!	The last point falls on the last node
		pointX(iSeg,numPointsT) = segX(iSeg,2)
		pointY(iSeg,numPointsT) = segY(iSeg,2)
		end do
! 	Do iSeg = 1,numSegs
! 		numPoints(iSeg) = int(1.0d0 + 1.0001d0*segLength(iSeg)/dXgrid)			! round up just a bit
! 		if(numPoints(iSeg).eq.1)numPoints(iSeg) = 2						! 2 points minimum = endpoints at nodes
! 		segdXgrid(iSeg) = dXgrid								! here I'm using a constant grid spacing, but this could be changed
! 		!dX = (nodeX(segNodes(iSeg,2)) - nodeX(segNodes(iSeg,1)))/(numPoints(iSeg)-1)
! 		!dY = (nodeY(segNodes(iSeg,2)) - nodeY(segNodes(iSeg,1)))/(numPoints(iSeg)-1)
! 		dX = (segX(iSeg,2) - segX(iSeg,1))/(numPoints(iSeg)-1)
! 		dY = (segY(iSeg,2) - segY(iSeg,1))/(numPoints(iSeg)-1)
! 		pointX(iSeg,1) = segX(iSeg,1)
! 		pointY(iSeg,1) = segY(iSeg,1)
! 		do j = 2,numPoints(iSeg)-1
! 			pointX(iSeg,j) = segX(iSeg,1) + dX*float(j-1)
! 			pointY(iSeg,j) = segY(iSeg,1) + dY*float(j-1)
! 			end do
! 		!	The last point falls on the last node
! 		pointX(iSeg,numPoints(iSeg)) = segX(iSeg,2)
! 		pointY(iSeg,numPoints(iSeg)) = segY(iSeg,2)
! 		end do

	!	The reason the first and last points fall on the nodes is that these are used in calculation of
	!	diffusion for points 2 to numPoints-1 in Subroutine DiffuseSegment
	!	The first and last point compositions are then set to the Node composition after DiffuseNodes is done
	write(66,*)'======================='
	write(66,*)'Segments'
	write(66,*)'      j      node1    node2    length    numPtsStart   numPtsEnd in_Seg'
		  !       1       1       2       1.000      11
	do j = 1,numSegs	
		write(66,65)j,segNodes(j,1),segNodes(j,2),segLength(j),numPointStart(j),numPointEnd(j)
65		format(3I8,F12.3,2I8)
		end do
	write(66,*)'Number of segments = ',numSegs

	do iNode = 1,numNodes
		do j = 1,numNodeSegs(iNode)
			i = nodeNodeConnect(iNode,j)
			do k = 1,numSegs
				if(segNodes(k,1).eq.i.or.segNodes(k,2).eq.i)then
					if(segNodes(k,1).eq.iNode)then
						! nodeSegConnect already has this info
						!nodePointSeg(iNode,j) = k			! index of segment that attaches to this node
						nodePointNextTo(iNode,j) = 2			! The initial value: the second point in this segment (the first point is coincident with the node)
						nodePointOnTop(iNode,j) = 1			! The initial value: the first point in this segment -- coincident with the node
						go to 72
						endif
					if(segNodes(k,2).eq.iNode)then
						!nodePointSeg(iNode,j) = k			! index of segment that attaches to this node
						nodePointNextTo(iNode,j) = numPointEnd(k)-1	! the second to last point in this segment (the last point is coincident with the node)
						nodePointOnTop(iNode,j)  = numPointEnd(k)	! the last point in this segment -- coincident with the node)
						go to 72
						endif
					endif
				end do
			write(*,*)' Did not find segment for this node'
			write(*,*)' node    numNodeSegs  '
			write(*,*)i,j
			pause ' hit return to continue'
72			continue
			end do
		end do


	write(66,*)' '
	write(66,*)'   node    connect        OnTop    connect        OnTop    connect        OnTop'
	write(66,*)'  index     node    seg   point     node    seg   point     node    seg   point'
	do iNode = 1,numNodes
		write(66,77)iNode,(nodeNodeConnect(iNode,j),nodeSegConnect(iNode,j),nodePointOnTop(iNode,j),j = 1,numNodeSegs(iNode))
77		format(15I8)
		end do

	write(66,*)'End setting up points in Segments'
	write(66,*)'======================='

	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine SetUpCrystals()
!	Set up crystals for the grid

	implicit none
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	integer*4 j,index,i,k,iNode,iCount,iseg,ii,iii,iConnect1,iConnect2,jj,kcur
	integer*4 MIFID1,MIFID2,MIFID3,MIF1temp,MIF2temp,nodeStart,Xltemp1,Xltemp2

!	Index the crystals
!	A crystal is made of a polygon
!	Polygon vertices are read in from the Grid_Maker file 

!	Crystal node indicies go clockwise from lower left corner.

	numCrystals = numPolys

	do i = 1,numCrystals
		numCrystalNodes(i) = numPolyNodes(i)
		do j = 1,numPolyNodes(i)
			crystalNodes(i,j) = PolyNodes(i,j)
			end do
		end do
		

	write(66,*)' ----'
	write(66,*)'numCrystals = ',numCrystals
	write(66,*)' ----'

	write(66,*)' '
	write(66,*)' List all phases for all crystals'
	write(66,*)' k MIDFID   Name'
	do k = 1,numCrystals
!		write(66,*)k,CrystalMIFID(k),PhaseName(CrystalMIFID(k))
		write(66,*)k,CrystalMIFID(k),CrystalPhaseName(k)
		end do

	write(66,*)' '
	write(66,*)' Phases at each node - '
	
	! Determine the phases around each node and store in nodeMIFID (to be used in MDF routines)				
!	Now determine the 3 phases at each reaction node
!	Note that some nodes have only 2 or 1 phases. These will not be reaction nodes
! 	do iNode = 1,numNodes			! loop through every node
! 		icount = 0
! 		do k = 1,numCrystals
! 			do i = 1,numCrystalNodes(k)			! number of nodes for each crystal = numPolyNodes
! 				if(CrystalNodes(k,i).eq.iNode)then	! this node is part of this crystal
! 					icount = icount + 1
! 					nodeMIFID(iNode,icount) = CrystalMIFID(k)
! 					endif
! 				end do
! 			end do
! 		numCrystalsAtNode(iNode) = icount		! This would be the same number as numPolysAtNode (but there is no such variable)
! 		end do
	! this code will assign crystals to nodes.
	! Crystal order should be the same as nodeConnect order and nodeSegConnect order
	do iNode = 1,numNodes
		icount = 0
		do j = 1,3	! check node connect order
			jj = j+1
			iConnect1 = nodeNodeConnect(iNode,j)
			if(j.eq.3)jj = 1
			iConnect2 = nodeNodeConnect(iNode,jj)
			do k = 1,numCrystals
				do i = 1,numCrystalNodes(k)
					if(CrystalNodes(k,i).eq.iNode)then		! this crystal poly has this node
						do ii = 1,numCrystalNodes(k)
							if(CrystalNodes(k,ii).eq.iConnect1)then	! this crystal poly has the connected node as well
								do iii = 1,numCrystalNodes(k)
									if(CrystalNodes(k,iii).eq.iConnect2)then
										icount = icount + 1
										nodeMIFID(iNode,icount) = CrystalMIFID(k)
										nodeCrystalIndex(iNode,icount) = k	! k is the crystal number
										endif
									end do
								end if
							end do
						endif
					end do
				end do
			end do
		numCrystalsAtNode(iNode) = icount		! This would be the same number as numPolysAtNode (but there is no such variable)
		end do


	write(66,*)' --------------------'
	write(66,*)' Phases around each node (all 3)'
	write(66,*)' iNode    NodeMIFID(3)   NumCrystalsAtNode'
	do iNode = 1,numNodes
		write(66,*)iNode,(NodeMIFID(iNode,i),i=1,3),numCrystalsAtNode(iNode)
		end do
	! check to see how many phases around each node are different
	! This will determine which MDF routine is invoked
	write(66,*)' --------------------'
	write(66,*)' Unique phases around each node'
                    !       1          5      11      60       2       7       2
	write(66,*)'  Node       MIF IDs                Reaction phases'
	do iNode = 1,numNodes			! loop through every node
		do k = 1,numCrystalsAtNode(iNode)
			if(numCrystalsAtNode(iNode).eq.1)then
				numNodeReactionPhases(iNode) = 0	! no reaction
				endif
			if(numCrystalsAtNode(iNode).eq.2)then
				if(nodeMIFID(iNode,1).eq.nodeMIFID(iNode,2))then
					numNodeReactionPhases(iNode) = 0	! no reaction
					else
					numNodeReactionPhases(iNode) = 2	! 2-phase reaction
					nodeReactionPhases(iNode,1) = 1
					nodeReactionPhases(iNode,2) = 2
					endif
				endif					
			if(numCrystalsAtNode(iNode).eq.3)then
				! with 3 phases at a node the options are
				! All phases the same (no reaction)
				! All phases different (3-phase reaction)
				! 2 phases the same -- either 1 & 2 or 1 & 3 (2-phase reaction)

			!	call AssignReactionPhases(iNode)
			!	go to 444
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
					nodeSegTheSamePhases(iNode) = 0
					go to 20
					endif	
				! if here, then 2 phases must be the same
				numNodeReactionPhases(iNode) = 2
				if(MIFID1.eq.MIFID2)then
					nodeReactionPhases(iNode,1) = 1
					nodeReactionPhases(iNode,2) = 3
					nodeSegTheSamePhases(iNode) = 2	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
					go to 20
					endif
				if(MIFID1.eq.MIFID3)then
					nodeReactionPhases(iNode,1) = 2
					nodeReactionPhases(iNode,2) = 3
					nodeSegTheSamePhases(iNode) = 1	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
					go to 20
					endif
				if(MIFID2.eq.MIFID3)then
					nodeReactionPhases(iNode,1) = 1
					nodeReactionPhases(iNode,2) = 2
					nodeSegTheSamePhases(iNode) = 3	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
					go to 20
					endif
				! if here, we have a problem
				write(12,*)' In routine SetUpCrystals'
				write(12,*)' Failure to determine nodeReactionPhases and nodeSegTheSamePhases'
				write(12,*)' iNode = ',iNode
				write(12,227)iNode,(nodeMIFID(iNode,j),j=1,numCrystalsAtNode(iNode)),			&
						(nodeReactionPhases(iNode,j),j=1,numNodeReactionPhases(iNode)),nodeSegTheSamePhases(iNode)
			227	format(20I8)
			
			20	continue	


				endif		
		
444			continue	! this is just a dummy jump to see if new subroutine works


			end do   ! end crystal loop		
		select case (numNodeReactionPhases(iNode))
		case(0)
			write(66,*)iNode,(nodeMIFID(iNode,j),j=1,numCrystalsAtNode(iNode)), 'No reaction phases'
		case(2)
			if(nodeSegTheSamePhases(iNode).eq.0)then
				write(66,230)iNode,(nodeMIFID(iNode,j),j=1,numCrystalsAtNode(iNode)),			&
					numNodeReactionPhases(iNode),(nodeReactionPhases(iNode,j),j=1,numNodeReactionPhases(iNode)),&
					nodeSegTheSamePhases(iNode)
230				format(7I8,8x,I8,'  Edge node')
				else
				write(66,229)iNode,(nodeMIFID(iNode,j),j=1,numCrystalsAtNode(iNode)),			&
					numNodeReactionPhases(iNode),(nodeReactionPhases(iNode,j),j=1,numNodeReactionPhases(iNode)),&
					nodeSegTheSamePhases(iNode),nodeSegConnect(iNode,nodeSegTheSamePhases(iNode))
				endif
229				format(7I8,8x,2I8)
		case(3)
			write(66,228)iNode,(nodeMIFID(iNode,j),j=1,numCrystalsAtNode(iNode)),			&
				numNodeReactionPhases(iNode),(nodeReactionPhases(iNode,j),j=1,numNodeReactionPhases(iNode)),	&
				nodeSegTheSamePhases(iNode)
228			format(20I8)
		case default
		write(66,*)iNode,'   Something went wrong'
		end select
		end do	! end this node


!	Now determine the 2 phases along each segment
!	Note that some segments have only 1 phases.
!	Also, some segments have the same 2 phases.
!	 These will not be reaction segments
	write(66,*)' '
	! first set the CrystalSegs to the PolySegs
	do i = 1,numPolys
		do j = 1,numPolyNodes(i)		! this is the same as numPolySegs
			CrystalSegs(i,j) = PolySegs(i,j)
			end do
		end do


	! Now find the 1 or 2 crystals associated with each segment

!	write(66,*)' Phases at each segment'
	!numSegsWith2Phases = 0
	do iSeg = 1,numSegs			! loop through every segment
!	write(66,*)'iseg, nodes', iseg,segNodes(iSeg,1),segNodes(iseg,2)
		icount = 0
		do i = 1,numPolys
			do j = 1,numPolyNodes(i)		! The number of crystalsegs = number of poly nodes
				if(CrystalSegs(i,j).eq.iSeg)then
					icount = icount+1
					segMIFID(iSeg,iCount) = CrystalMIFID(i)		! This is the MIFID for the crystal associated with the seg	
					segCrystalID(iSeg,iCount) = i			! This is one of the 2 crystal numbers associated with this segment -- same as the poly number
					go to 15		
					endif
				end do
15			continue
			end do
		! now determine whether the phases are the same or not
		numSegReactionPhases(iSeg) = 0
		if(icount.eq.2)then
			if(segMIFID(iSeg,1).ne.segMIFID(iSeg,2))then		! they are the same
				numSegReactionPhases(iSeg) = 2
				endif
		! Note: if icount=0 or icount=1 then there is no reaction along this segment	
			endif
		end do

        ! Now go through every segment to be sure that segMIF1 is to the right of the segment (clockwise)
        !   If not, then swap segMIF1 and segMIF2 (this is required in the Subroutine GrowSegs so we know which direction to move the points
	! The phases and segments around every node are anticlockwise
	! and the segments are ordered so that they are "behind" (i.e. clockwise) to the phase
! 	
! 	          Seg1
! 	           |
! 	           |
! 	     Ph1   |
! 	           |     Ph3
! 	          / \
! 	         /   \
! 	        /     \
! 	       /  Ph2  \
!              Seg2      Seg3
! 	
	
	
	do iSeg = 1,numSegs
		nodeStart = segNodes(iSeg,1)
		! we only want to look around the starting node for each segment, not the ending node
		do j = 1,numNodeSegs(nodeStart)		
			if(nodeSegConnect(nodeStart,j).eq.iSeg)then		! this is the one we want
				if(j.eq.1)then
					jj = numNodeSegs(nodeStart)
					else
					jj = j-1
					endif
				if(segMIFID(iSeg,1).eq.nodeMIFID(nodeStart,jj))then	
					go to 55		! this one is set up correctly
					else			!swap the segMIFID and the segCrystalID
					MIF2temp = segMIFID(iSeg,1)
					MIF1temp = segMIFID(iSeg,2)
					segMIFID(iSeg,1) = MIF1temp
					segMIFID(iSeg,2) = MIF2temp
					Xltemp2 = segCrystalID(iSeg,1)
					Xltemp1 = segCrystalID(iSeg,2)
					segCrystalID(iSeg,1) = Xltemp1
					segCrystalID(iSeg,2) = Xltemp2
					endif
				endif
			end do
55		continue
		end do

	write(66,*)' '
	write(66,*)'iSeg,numSegReactionPhases,segMIFID,segCrystalID '
	do iSeg = 1,numSegs
		write(66,*)iSeg,numSegReactionPhases(iSeg),(segMIFID(iSeg,j),j=1,numSegReactionPhases(iSeg)),	&
						       (segCrystalID(iSeg,j),j=1,numSegReactionPhases(iSeg))
		end do
				
	! calculate and store the center of each crystal
	write(66,*)' '
	write(66,*)' Crystal centers, nodes and segments'
	write(66,*)'   i           X            Y             Nodes               Segments'
	Do index = 1,numCrystals
		CrystalCenterX(index) = PolyCenterX(index)
		CrystalCenterY(index) = PolyCenterY(index)
		write(66,1066)index,CrystalCenterX(index),CrystalCenterY(index),			&
			(CrystalNodes(index,j),j=1,numCrystalNodes(index)),(CrystalSegs(index,j),j=1,numCrystalNodes(index))
1066		format(I5,2F12.5,40I5)
		end do


	! Array segChange(maxsegs) is a flag that indicates whether the points along a segment will change after MDF reactions
	! This is necessary because when a node moves, the length of segments and the positions of points along a segment change
	! This array is only used by subroutine ResetSegPoints()
	! Set the array segChange = 0 (no change) = 1 (change)

! 	do iSeg = 1,numSegs
! 		segChange(iSeg) = 0		! initialize to all zero
! 		end do
! 	do iNode = 1,numNodes
! 		if(numNodeReactionPhases(iNode).eq.0)cycle	! there is no reaction at this node, so segs don't change
! 		do i = 1,numNodeSegs(iNode)
! 			j = nodeSegConnect(iNode,i)
! 			segChange(j) = 1
! 			end do
! 		end do 
	
	! Cycle through all of the nodes --
	! For nodes with reaction phases, calculate the number of MDF reactions
	! We do this because in cases where numMDFRxnsNodes = 0, we will use an equilibrium model
	write(66,*)' '
	write(66,*)' Number of MDF reactions in each node'
	do iNode = 1,numNodes
!		if(numNodeReactionPhases(iNode).eq.0)cycle	! there is no reaction at this node, so segs don't change
		select case(numNodeReactionPhases(iNode))
			case(0)
				cycle
			case(3)
				numPh = 4
				asmCurrent(1) = 1			! always the grain boundary
				asmCurrent(2) = nodeMIFID(iNode,nodeReactionPhases(iNode,1))
				asmCurrent(3) = nodeMIFID(iNode,nodeReactionPhases(iNode,2))
				asmCurrent(4) = nodeMIFID(iNode,nodeReactionPhases(iNode,3))

			case(2)	
				numPh = 3
				asmCurrent(1) = 1			! always the grain boundary
				asmCurrent(2) = nodeMIFID(iNode,nodeReactionPhases(iNode,1))
				asmCurrent(3) = nodeMIFID(iNode,nodeReactionPhases(iNode,2))
			case default
			end select	
		NP=0
      		do kCur = 1,numPh
			k = asmCurrent(kCur)
			NP = numPhCo(K)+NP
			end do
      		NX=NP-numPh		! This is the number of independent compositional variables (first X in each phase is dependent)
!		call SetALLX()		! set up ALLX array
		call Names()		! sets up names of variables for asmCurrent
      		call REXN		! calculate linearly independent reactions
!		numMDFRxnsNodes(iNode) = NEQ		!NEQ is either the number of EQ reactions OR number of MSDF reactions
		numMDFRxnsNodes(iNode) = numMDFrxn
		write(66,*)iNode,numMDFRxnsNodes(iNode)
		end do

	write(66,*)' Number of MDF reactions in each segment'
	do iSeg = 1,numSegs
		if(numSegReactionPhases(iSeg).eq.0)cycle	! there is no reaction at this node, so segs don't change
		numPh = 3
		asmCurrent(1) = 1			! always the grain boundary
		asmCurrent(2) = segMIFID(iSeg,1)
		asmCurrent(3) = segMIFID(iSeg,2)
		NP=0
      		do kCur = 1,numPh
			k = asmCurrent(kCur)
			NP = numPhCo(K)+NP
			end do
      		NX=NP-numPh		! This is the number of independent compositional variables (first X in each phase is dependent)
!		call SetALLX()		! set up ALLX array
		call Names()		! sets up names of variables for asmCurrent
      		call REXN		! calculate linearly independent reactions
		!numMDFRxnsSegs(iSeg) = NEQ	
		numMDFRxnsSegs(iSeg) = numMDFrxn	
		write(66,*)iSeg,numMDFRxnsSegs(iSeg)
		end do


	return
	end



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE InitializePhaseComps()
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
! *****************************************
	integer*4 k1,k,j,l,iNode,jj,j1,j2,iSeg

	! Initialize all nodes and segments with phase compositions from input MIF

	write(66,*)' '
	write(66,*)'*********************************************'
	write(66,*)'Initialize nodes with phase compositionss'
	write(66,*)'iNode,numNodeReactionPhases(iNode),nodeMIFID(iNode,1),nodeMIFID(iNode,2),nodeMIFID(iNode,3)'
	write(66,*)'iNode,numNodeReactionPhases(iNode),nodeMIFID(iNode,j1),nodeMIFID(iNode,j2)'
	write(66,*)'iNode,numNodeReactionPhases(iNode)'
	do iNode = 1,numNodes
		select case(numNodeReactionPhases(iNode))
		case(3)
			!if(numNodeReactionPhases(iNode).eq.3)then
			write(66,*)iNode,numNodeReactionPhases(iNode),nodeMIFID(iNode,1),nodeMIFID(iNode,2),nodeMIFID(iNode,3)
			do k1 = 1,3
				k = nodeMIFID(iNode,k1)			! this is the phase number from the MIF
				nodePhaseMoles(iNode,k1) = 0.0d0
				nodePhaseMolesDelta(iNode,k1) = 0.0d0
				do j = 1,numPhCo(k)
					nodePhaseComp(iNode,k1,j) = xPhCo(k,j)
					end do
				end do
		case(2)
			j1 = nodeReactionPhases(iNode,1)	! should be 1
			j2 = nodeReactionPhases(iNode,2)	! should be 2 or 3
			
			write(66,*)iNode,numNodeReactionPhases(iNode),nodeMIFID(iNode,j1),nodeMIFID(iNode,j2)
			do k1 = 1,2
				jj = nodeReactionPhases(iNode,k1)
				k = nodeMIFID(iNode,jj)			! this is the phase number from the MIF
				nodePhaseMoles(iNode,jj) = 0.0d0
				nodePhaseMolesDelta(iNode,k1) = 0.0d0
				do j = 1,numPhCo(k)
					nodePhaseComp(iNode,jj,j) = xPhCo(k,j)
					end do
				end do
		case default
			write(66,*)iNode,numNodeReactionPhases(iNode)
		! no reaction if numNodeReaction phases isn't 2 or 3
		end select
		end do		! end loop on nodes
	


	write(66,*)'Initialize compositions of all segs with 2 phases and points'
	write(66,*)' iSeg    MIFID1   MIFID2'
	do iSeg = 1,numSegs
		select case (numSegReactionPhases(iSeg))
			case(2)
				write(66,*)iSeg,numSegReactionPhases(iSeg),segMIFID(iSeg,1),segMIFID(iSeg,2)
				do k1 = 1,2
					k = segMIFID(iSeg,k1)
!					do l = 1,numPoints(iSeg)
					do l = numPointStart(iSeg),numPointEnd(iSeg)
						pointPhaseMoles(iSeg,l,k1) = 0.0d0
						do j = 1,numPhCo(k)
							pointPhaseComp(iSeg,l,k1,j) = xPhCo(k,j)
							end do
						end do
					end do
			case default
				write(66,*)iSeg,numSegReactionPhases(iSeg)
			end select
		end do

		! 	do i1 = 1,numSegsWith2Phases
		! 		i = segReactionPhases(i1)
		! 		write(66,*)i1,i,segMIFID(i,1),segMIFID(i,2)
		! 		do k1 = 1,2
		! 			k = segMIFID(i,k1)
		! 			do l = 1,numPoints(i)
		! 				pointPhaseMoles(i,l,k1) = 0.0d0
		! 				do j = 1,numPhCo(k)
		! 					pointPhaseComp(i,l,k1,j) = xPhCo(k,j)
		! 					end do
		! 				end do
		! 			end do
		! 		end do

	write(66,*)'Initialize all nodes and points with the GB composition'
	write(66,*)' Node, Mn comp'
	do iNode = 1,numNodes
		do j = 1,numPhCo(1)
			NodeComp(iNode,j) = xPhCo(1,j)
			end do
		write(66,*)iNode,NodeComp(iNode,6)
		end do
	write(66,*)' Seg, point, Mn comp'
	do iSeg = 1,numSegs
!		do k = 1,numPoints(iSeg)
		do k = numPointStart(iSeg),numPointEnd(iSeg)
			do j = 1,numPhCo(1)
				pointComp(iSeg,k,j) = xPhCo(1,j)
				end do
			write(66,*)iSeg,k,pointComp(iSeg,k,6)
			end do
		end do			

	write(66,*)'End Initialization'
	write(66,*)'*********************************************'

	return
	end	


		
		
		
