	Subroutine Grow3Nodes(iNode,moleScaleFactor,debug)
	USE AWE_INTERFACES
	! calling routine for Sub MoveNode to get the relevant information for a specific node
	! debug = 1: generate output
	! debug = 0: no output
	implicit none	
	TYPE(AWE_Canvas) :: NodePlot
	include "Diffuse_GB_Hex.inc"
	include "GibbsFiles.inc"
	include "Assemb.inc"

	integer*4 iRx1,iRx2,iRx3,iNode,iSeg1,iSeg2,iSeg3,iPoint1,iPoint2,iPoint3,MIFID1,MIFID2,MIFID3,MIF1
	integer*4 debug,err
	real*8 Nold(2),P1(2),P2(2),P3(2),MPh1,MPh2,MPh3,Nnew(2),moleScaleFactor,MPh1T,MPh2T,MPh3T

! 	X,Y of the node
	Nold(1) = NodeX(iNode)
	Nold(2) = NodeY(iNode)

	! connecting segments and points that form the triangle around the node
 	iseg1 = nodeSegConnect(iNode,1)
 	iPoint1 = nodePointNextTo(iNode,1)
	P1(1) = PointX(iSeg1,iPoint1)
	P1(2) = PointY(iSeg1,iPoint1)
 	iseg2 = nodeSegConnect(iNode,2)
 	iPoint2 = nodePointNextTo(iNode,2)
	P2(1) = PointX(iSeg2,iPoint2)
	P2(2) = PointY(iSeg2,iPoint2)
 	iseg3 = nodeSegConnect(iNode,3)
 	iPoint3 = nodePointNextTo(iNode,3)
	P3(1) = PointX(iSeg3,iPoint3)
	P3(2) = PointY(iSeg3,iPoint3)


	if(debug.eq.1)then
		write(12,*)'   '
		write(12,*)'   '
		write(12,*)'*&*&*&*&*&*&*&*&*&*&**&*&*&*&*&**'
		write(12,*)' Test Grow3Nodes routine  '
		write(12,*)' Node  ',iNode,NodeX(iNode),NodeY(iNode)
		write(12,*)' Segments and points '
		write(12,104)iSeg1,iPoint1,P1(1),P1(2)
		write(12,104)iSeg2,iPoint2,P2(1),P2(2)
		write(12,104)iSeg3,iPoint3,P3(1),P3(2)
104		format(2I8,2F12.5)
		
		! plot the results if we're debugging
		call PlotNode(NodePlot,iNode,iSeg1,iPoint1,P1,iSeg2,iPoint2,P2,iSeg3,iPoint3,P3)



		endif


	! Get the MIF IDs so we know what the 3 phases are. 
	! These are needed so we can scale the moles to the number of oxygens
! 	MIFID1 = nodeMIFID(iNode,1)
! 	MPh1 = nodePhaseMolesDelta(iNode,1)*numOxygens(MIFID1)
! 	MIFID2 = nodeMIFID(iNode,2)
! 	MPh2 = nodePhaseMolesDelta(iNode,2)*numOxygens(MIFID2)
! 	MIFID3 = nodeMIFID(iNode,3)
! 	MPh3 = nodePhaseMolesDelta(iNode,3)*numOxygens(MIFID3)

	iRx1 = nodeReactionPhases(iNode,1)
	MIFID1 = nodeMIFID(iNode,iRx1)
	!MPh1T = nodePhaseMolesDelta(iNode,iRx1)*numOxygens(MIFID1)
	MPh1T = nodePhaseMolesDelta(iNode,iRx1)
	iRx2 = nodeReactionPhases(iNode,2)
	MIFID2 = nodeMIFID(iNode,iRx2)
	!MPh2T = nodePhaseMolesDelta(iNode,iRx2)*numOxygens(MIFID2)
	MPh2T = nodePhaseMolesDelta(iNode,iRx2)
	iRx3 = nodeReactionPhases(iNode,3)
	MIFID3 = nodeMIFID(iNode,iRx3)
	!MPh3T = nodePhaseMolesDelta(iNode,iRx3)*numOxygens(MIFID3)
	MPh3T = nodePhaseMolesDelta(iNode,iRx3)

	! This just scales to moles -- not sure what the scaling factor should be
	MPh1T = MPh1T*moleScaleFactor
	MPh2T = MPh2T*moleScaleFactor
!	MPh3 = MPh3*moleScaleFactor
	Mph3T = -(MPh1T + MPh2T)		! this ensures that the sum = 0
					! Probably not needed during a model run, but when you read in from a .GBM model there is round off errors

	if(debug.eq.1)then
		write(12,*)'Node reactionPhasese(iNode,i)=iRxn'
		write(12,*)iRx1,iRx2,iRx3
		write(12,*)' '
		write(12,*)'     MIFID    PhName    Del_moles      Del_area (delta_moles*scale)'
		write(12,103)MIFID1,PhName(MIFID1),nodePhaseMolesDelta(iNode,iRx1),MPh1T
		write(12,103)MIFID2,PhName(MIFID2),nodePhaseMolesDelta(iNode,iRx2),MPh2T
		write(12,103)MIFID3,PhName(MIFID3),nodePhaseMolesDelta(iNode,iRx3),MPh3T
		write(12,*)'Sum = ',Mph1+MPh2+MPh3
103		format(T8,I4,2x,A8,2E15.5,5x,I8,2E5.5)
		endif


!	The phases and segments around a node are anti-clockwise
!	However, they phases and segment sequences may start in different places
!	This code will ensure that phase 1 is anti-clockwise from segment 1, 
!		which will ensure that all nodes are passed to subroutine MoveNode in the same order

!	Find the 2 phases that straddle segments 1 and 2
	! Phases that straddle Seg1

	if(segMIFID(iSeg1,1).eq.segMIFID(iSeg2,1).or.segMIFID(iSeg1,1).eq.segMIFID(iSeg2,2))then
		MIF1 = segMIFID(iSeg1,1)
		go to 5
		endif
	if(segMIFID(iSeg1,2).eq.segMIFID(iSeg2,1).or.segMIFID(iSeg1,2).eq.segMIFID(iSeg2,2))MIF1 = segMIFID(iSeg1,2)
	
5	continue
	! now order the phases accordingly
	if(MIF1.eq.MIFID1)then
		MPh1 = MPh1T
		MPh2 = MPh2T
		MPh3 = MPh3T
		go to 6
		endif
	if(MIF1.eq.MIFID2)then
		MPh1 = MPh2T
		MPh2 = MPh3T
		MPh3 = MPh1T
		go to 6
		endif
	if(MIF1.eq.MIFID2)then
		MPh1 = MPh2T
		MPh2 = MPh3T
		MPh3 = MPh1T
		go to 6
		endif
	if(MIF1.eq.MIFID3)then
		MPh1 = MPh3T
		MPh2 = MPh1T
		MPh3 = MPh2T
		go to 6
		endif

6	continue	

	if(debug.eq.1)then
		write(12,*)' '
		write(12,*)' Adjusted moles around the node'
		write(12,*)' Phase 1 falls between seg1 and seg2 and all go anti-clockwise'
		write(12,*)'               segs       Scaled moles of phases'
		write(12,*)' seg1-seg2 ',iSeg1,iSeg2,MPh1
		write(12,*)' seg2-seg3 ',iSeg2,iSeg3,MPh2
		write(12,*)' seg3-seg1 ',iSeg3,iSeg1,MPh3

		write(12,*)' ---------------------------------'
		write(12,*)' Calling MoveNode '
		write(12,*)' ---------------------------------'
		endif

!11	continue
	err = 0
! 	Call MoveNode(iNode,P1,P2,P3,Nold,MPh1,MPh2,MPh3,Nnew,debug,err)
	Call MoveNode2(iNode,P1,P2,P3,Nold,MPh1,MPh2,MPh3,Nnew,debug,err)
	if(err.eq.1)then
		write(75,*)'Failure in Sub MoveNode called from Grow3Nodes. Node = ',iNode
		write(12,*)'Failure in Sub MoveNode called from Grow3Nodes. Node = ',iNode
		write(12,*)' Cycle = ',totalCycles
		!pause 'pausing'
		return		!Bail out
		endif
	NodeX(iNode) = Nnew(1)
	NodeY(iNode) = Nnew(2)
	if(debug.eq.1)then
		write(12,*)' Old node coordinates: ',Nold(1),Nold(2)
		write(12,*)' New node coordinates: ',Nnew(1),Nnew(2)
		write(12,*)' NodeX and nodeY:      ',NodeX(iNode),NodeY(iNode)
		Call PlotNewNode(NodePlot,Nnew)
		endif
		
!	Try to figure out why the node doesn't move
! 	if(iNode.eq.11)then
! 		write(12,*)NodeX(iNode),NodeY(iNode)
! 		endif

	! I need to set segX and segY when I move the nodes
		if(segNodes(iSeg1,1).eq.iNode)then
			segX(iSeg1,1) = NodeX(iNode)
			segY(iSeg1,1) = NodeY(iNode)
			pointX(iSeg1,numPointStart(iSeg1)) = segX(iSeg1,1)
			pointY(iSeg1,numPointStart(iSeg1)) = segY(iSeg1,1)
			else
			segX(iSeg1,2) = NodeX(iNode)
			segY(iSeg1,2) = NodeY(iNode)
			pointX(iSeg1,numPointEnd(iSeg1)) = segX(iSeg1,2)
			pointY(iSeg1,numPointEnd(iSeg1)) = segY(iSeg1,2)
			endif
		if(segNodes(iSeg2,1).eq.iNode)then
			segX(iSeg2,1) = NodeX(iNode)
			segY(iSeg2,1) = NodeY(iNode)
			pointX(iSeg2,numPointStart(iSeg2)) = segX(iSeg2,1)
			pointY(iSeg2,numPointStart(iSeg2)) = segY(iSeg2,1)
			else
			segX(iSeg2,2) = NodeX(iNode)
			segY(iSeg2,2) = NodeY(iNode)
			pointX(iSeg2,numPointEnd(iSeg2)) = segX(iSeg2,2)
			pointY(iSeg2,numPointEnd(iSeg2)) = segY(iSeg2,2)
			endif
		if(segNodes(iSeg3,1).eq.iNode)then
			segX(iSeg3,1) = NodeX(iNode)
			segY(iSeg3,1) = NodeY(iNode)
			pointX(iSeg3,numPointStart(iSeg3)) = segX(iSeg3,1)
			pointY(iSeg3,numPointStart(iSeg3)) = segY(iSeg3,1)
			else
			segX(iSeg3,2) = NodeX(iNode)
			segY(iSeg3,2) = NodeY(iNode)
			pointX(iSeg3,numPointEnd(iSeg3)) = segX(iSeg3,2)
			pointY(iSeg3,numPointEnd(iSeg3)) = segY(iSeg3,2)
			endif

	return
	end
	


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine MoveNode(iNode,P1,P2,P3,Nold,MPh1T,MPh2T,MPh3T,Nnew,debug,err)
	! Routine to calculate the new position of a node located inside of a triangle
	! given the amount the areas of the 3 sub triangles change
	use MatrixArrays
	implicit none	
	External Area
	real*8 Area
	integer*4 icount,i,j,nsolve,numEq,ier,pass,output,debug,err,iNode,PlotsOutside,PlotsNearTheEdge
	real*8 P1(2),P2(2),P3(2)					! 3 points that define the outer triangle
	real*8 Nold(2),Nnew(2)						! node in the center of the triangle
	real*8 MPh1,MPh2,MPh3,MPh1T,MPh2T,MPh3T
	real*8 Area123,Area12Nold,Area13Nold,Area23Nold,Area12Nnew,Area13Nnew,Area23Nnew
	real*8 Area12NdXplus, Area13NdXplus, Area23NdXplus, Area12NdYplus, Area13NdYplus, Area23NdYplus, AreaSum
	real*8 Area12NdXminus,Area13NdXminus,Area23NdXminus,Area12NdYminus,Area13NdYminus,Area23NdYminus
	real*8 F12N,F13N,F23N,dF12NdX,dF12NdY,dF13NdX,dF13NdY,dF23NdX,dF23NdY,Target1,Target2,Target3
	real*8 delta,deta,tol,damp,R2,abit,tempDelta
	real*8 Asign12,Asign23,Asign13
	
	tol = 0.001d0	! tolerance for convergence of areas. The zero functions are Target - area. So units are in areas, which are on the order of 1-10
	damp = 1.0d0
	delta = 0.0001d0

	! Assign moles 1,2,3 to the appropriate triangles
	! this is the correct assignment if....
	! the 3 phases are always anticlockwise around the node in the same order as the segments.
	! There is code in Subroutine Grow3Nodes that is supposed to order these correctly.
	MPh1 = MPh1T	! area12
	MPh2 = MPh2T	! area23
	MPh3 = MPh3T	! area13
	! Calculate the initial (original) areas of the 3 subtriangles
	Area123 = Area(P1,P2,P3)
	! These are initially set = 1 because the first call to Get Areas we need only the area of each triangle
	Asign12 = 1.0d0
	Asign23 = 1.0d0
	Asign13 = 1.0d0
	! on this first pass, all signs = 1 so all areas are returned as positive values
	call GetAreas(P1,P2,P3,Nold,Area12Nold,Area23Nold,Area13Nold,Asign12,Asign23,Asign13,debug)
	AreaSum = Area12Nold+Area23Nold+Area13Nold
	Target1 = MPh1 + Area12Nold
	Target2 = MPh2 + Area23Nold
	Target3 = Mph3 + Area13Nold
	if((AreaSum-Area123).gt.0.0d0)then
		PlotsOutside = 1
		else
		PlotsOutside = 0
		endif
	tempDelta = abs(MPh1)+abs(MPh2)+abs(MPh3)
! 	if(abs(Area123-AreaSum).lt.1e-12)then		! the node is very near the edge of the 123 triangle
	if(abs(Area123-AreaSum).lt.tempDelta)then		! the node is very near the edge of the 123 triangle
		PlotsNearTheEdge = 1
		else
		PlotsNearTheEdge = 0
		endif	
	if(debug.eq.1)then
		write(12,*)' '
		write(12,*)'----Before adjusting areasigns---------------- '		
		write(12,*)' Nold        ',Nold(1),Nold(2)
		write(12,*)'Area123        ',Area123
		write(12,*)'AreaSum        ',AreaSum
		write(12,*)'AreaSum-Area123',AreaSum-Area123
		!if(abs(Area123-AreaSum).gt.1.0d-10)then
		if(PlotsOutside.eq.1)then
			write(12,*)'Node plots outside of triangle'
			else
			write(12,*)'Node plots inside of triangle'
			endif
		if(PlotsNearTheEdge.eq.1)then
			write(12,*)'Node plots near the edge of triangle'
			else
			write(12,*)'Node plots far from the edge of triangle'
			endif
		write(12,*)'              Area 12N             Area 23N              Area 13N'
		write(12,*)'Old_area  ',Area12Nold,Area23Nold,Area13Nold
		write(12,*)'AreaDelta ',MPh1,MPh2,MPh3
		write(12,*)'Targets   ',Target1,Target2,Target3
		write(12,*)' ---------------------------------------------------'
		write(12,*)' '
		endif
	! There are 4 cases to consider
	!	1 Plots inside far from the edge
	!	2 Plots inside near the edge (does it move outside?)
	!	3 Plots outside far from the edge
	!	4 Plots outside near the edge (does it move inside?)

!	if(abs(Area123-AreaSum).lt.1e-12)then		! the node is very near the edge of the 123 triangle
	if(PlotsOutside.eq.0.and.PlotsNearTheEdge.eq.1)then
		if(debug.eq.1)then
			write(12,*)' Plots near the edge inside the triangle'
			write(12,*)'Working on determining if node will move outside of triangle'
			endif
		abit = 0.5			! this is the amount we are going to move along the line
		if(Target1.lt.0.0d0)then
			Call MoveOut(iNode,P1,P2,P3,Nold,abit,Nnew,debug,err)
			Nold(1) = Nnew(1)
			Nold(2) = NNew(2)
			Asign12 = -1.*Asign12
			! This call will get the new "old" areas with the new "Nold"
			call GetAreas(P1,P2,P3,Nold,Area12Nold,Area23Nold,Area13Nold,Asign12,Asign23,Asign13,debug)
			AreaSum = Area12Nold+Area23Nold+Area13Nold
			endif
		if(Target2.lt.0.0d0)then
			Call MoveOut(iNode,P2,P3,P1,Nold,abit,Nnew,debug,err)
			Nold(1) = Nnew(1)
			Nold(2) = NNew(2)
			Asign23 = -1.*Asign23
			! This call will get the new "old" areas with the new "Nold"
			call GetAreas(P1,P2,P3,Nold,Area12Nold,Area23Nold,Area13Nold,Asign12,Asign23,Asign13,debug)
			AreaSum = Area12Nold+Area23Nold+Area13Nold
			endif
		if(Target3.lt.0.0d0)then
			Call MoveOut(iNode,P1,P3,P2,Nold,abit,Nnew,debug,err)
			Nold(1) = Nnew(1)
			Nold(2) = NNew(2)
			Asign13 = -1.*Asign13
			! This call will get the new "old" areas with the new "Nold"
			call GetAreas(P1,P2,P3,Nold,Area12Nold,Area23Nold,Area13Nold,Asign12,Asign23,Asign13,debug)
			AreaSum = Area12Nold+Area23Nold+Area13Nold
			endif
		go to 5
		endif		

	if(PlotsOutside.eq.1.and.PlotsNearTheEdge.eq.1)then
		if(debug.eq.1)then
			write(12,*)' Plots near the edge outside the triangle'
			write(12,*)'Working on determining if node will move inside of triangle'
			endif
		abit = -0.5			! this is the amount we are going to move along the line
		if(Target1.lt.0.0d0)then
			Call MoveOut(iNode,P1,P2,P3,Nold,abit,Nnew,debug,err)
			Nold(1) = Nnew(1)
			Nold(2) = NNew(2)
			Asign12 = -1.*Asign12
			! This call will get the new "old" areas with the new "Nold"
			call GetAreas(P1,P2,P3,Nold,Area12Nold,Area23Nold,Area13Nold,Asign12,Asign23,Asign13,debug)
			AreaSum = Area12Nold+Area23Nold+Area13Nold
			endif
		if(Target2.lt.0.0d0)then
			Call MoveOut(iNode,P2,P3,P1,Nold,abit,Nnew,debug,err)
			Nold(1) = Nnew(1)
			Nold(2) = NNew(2)
			Asign23 = -1.*Asign23
			! This call will get the new "old" areas with the new "Nold"
			call GetAreas(P1,P2,P3,Nold,Area12Nold,Area23Nold,Area13Nold,Asign12,Asign23,Asign13,debug)
			AreaSum = Area12Nold+Area23Nold+Area13Nold
			endif
		if(Target3.lt.0.0d0)then
			Call MoveOut(iNode,P1,P3,P2,Nold,abit,Nnew,debug,err)
			Nold(1) = Nnew(1)
			Nold(2) = NNew(2)
			Asign13 = -1.*Asign13
			! This call will get the new "old" areas with the new "Nold"
			call GetAreas(P1,P2,P3,Nold,Area12Nold,Area23Nold,Area13Nold,Asign12,Asign23,Asign13,debug)
			AreaSum = Area12Nold+Area23Nold+Area13Nold
			endif
		go to 5
		endif		

	if(PlotsOutside.eq.1.and.PlotsNearTheEdge.eq.0)then
		if(debug.eq.1)then
			write(12,*)' Plots outside the triangle but not near the edge'
			write(12,*)'Working on determining which side we are outside of'
			endif
		Call FindIntersection(P1,P2,P3,Nold,Asign12,Asign23,Asign13,debug)   ! returns Asign = -1 depending on which side we are outside of
		Area12Nold = Area12Nold*Asign12
		Area23Nold = Area23Nold*Asign23
		Area13Nold = Area13Nold*Asign13
		AreaSum = Area12Nold+Area23Nold+Area13Nold
		Target1 = MPh1 + Area12Nold
		Target2 = MPh2 + Area23Nold
		Target3 = Mph3 + Area13Nold
		if(debug.eq.1)then
			write(12,*)' The Asign changes depending on which side we are outside of - see that this makes sense'
			write(12,*)' Targets should have changed -- check below to see if this is correct'
			write(12,*)'Asigns    ',Asign12,Asign23,Asign13
			write(12,*)'Old area  ',Area12Nold,Area23Nold,Area13Nold
			write(12,*)'AreaDelta ',MPh1,MPh2,MPh3
			write(12,*)'Targets   ',Target1,Target2,Target3
			endif
		! after this call, one area should be negative
		! Now see if the node plots near the edge and adjust if it will move inside
		if(PlotsNearTheEdge.eq.1)then
			if(debug.eq.1)then
				write(12,*)' '
				write(12,*)' Plots outside the triangle and near the edge'
				write(12,*)'Working on determining if node will move inside of triangle'
				endif
			abit = -0.5			! this is the amount we are going to move along the line
			if(Target1.lt.0.0d0)then
				Call MoveOut(iNode,P1,P2,P3,Nold,abit,Nnew,debug,err)
				Nold(1) = Nnew(1)
				Nold(2) = NNew(2)
				Asign12 = -1.*Asign12		! changes the sign
				!Asign12 = -1.
				! This call will get the new "old" areas with the new "Nold"
				call GetAreas(P1,P2,P3,Nold,Area12Nold,Area23Nold,Area13Nold,Asign12,Asign23,Asign13,debug)
				AreaSum = Area12Nold+Area23Nold+Area13Nold
				endif
			if(Target2.lt.0.0d0)then
				Call MoveOut(iNode,P2,P3,P1,Nold,abit,Nnew,debug,err)
				Nold(1) = Nnew(1)
				Nold(2) = NNew(2)
				Asign23 = -1.*Asign23		! changes the sign
				!Asign23 = -1.
				! This call will get the new "old" areas with the new "Nold"
				call GetAreas(P1,P2,P3,Nold,Area12Nold,Area23Nold,Area13Nold,Asign12,Asign23,Asign13,debug)
				AreaSum = Area12Nold+Area23Nold+Area13Nold
				endif
			if(Target3.lt.0.0d0)then
				Call MoveOut(iNode,P1,P3,P2,Nold,abit,Nnew,debug,err)
				Nold(1) = Nnew(1)
				Nold(2) = NNew(2)
				Asign13 = -1.*Asign13		! changes the sign
				!Asign13 = -1.
				! This call will get the new "old" areas with the new "Nold"
				call GetAreas(P1,P2,P3,Nold,Area12Nold,Area23Nold,Area13Nold,Asign12,Asign23,Asign13,debug)
				AreaSum = Area12Nold+Area23Nold+Area13Nold
				endif
			if(debug.eq.1)then
				write(12,*)' Plots near the edge outside the triangle'
				write(12,*)'Working on determining if node will move inside of triangle'
				write(12,*)' The Asign changes depending on which side we are outside of - see that this makes sense'
				write(12,*)' Targets should have changed -- check below to see if this is correct'
				write(12,*)'Asigns    ',Asign12,Asign23,Asign13
				write(12,*)'Old area  ',Area12Nold,Area23Nold,Area13Nold
				write(12,*)'Moles     ',MPh1,MPh2,MPh3
				write(12,*)'Targets   ',Target1,Target2,Target3
				endif

			endif	
		go to 5
		endif


5	continue

	if(debug.eq.1)then
		write(12,*)' '
		write(12,*)'----After adjusting area signs---------------- '		
		write(12,*)' Nold        ',Nold(1),Nold(2)
		write(12,*)'Area123        ',Area123
		write(12,*)'AreaSum        ',AreaSum
		write(12,*)'AreaSum-Area123',AreaSum-Area123
		write(12,*)'              Area 12N             Area 23N              Area 13N'
		write(12,*)'Old area  ',Area12Nold,Area23Nold,Area13Nold
		write(12,*)'areaDelta ',MPh1,MPh2,MPh3
		write(12,*)'Targets   ',Target1,Target2,Target3
		endif
	Nnew(1) = Nold(1)
	Nnew(2) = Nold(2)

	pass = 1


10	continue
	call GetAreas(P1,P2,P3,Nnew,Area12Nnew,Area23Nnew,Area13Nnew,Asign12,Asign23,Asign13,debug)

	if(debug.eq.1)then
		write(12,*)' '
		write(12,*)' '
		write(12,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
! 		write(12,*)' Pass = ',pass
		select case (pass)
			case(1)
				write(12,*)'Functions = F12 & F13'
			case(2)
				write(12,*)'Functions = F12 & F23'
			case(3)
				write(12,*)'Functions = F13 & F23'
			case default
			end select

		write(12,*)' Nnew     ',Nnew(1),Nnew(2)
		write(12,*)'New area  ',Area12Nnew,Area23Nnew,Area13Nnew
		write(12,*)'New-Old   ',Area12Nnew-Area12Nold,Area23Nnew-Area23Nold,Area13Nnew-Area13Nold
		F12N =  Target1 - Area12Nnew
		F23N =  Target2 - Area23Nnew
		F13N =  Target3 - Area13Nnew
		write(12,*)'F12N =    ',F12N
		write(12,*)'F23N =    ',F23N
		write(12,*)'F13N =    ',F13N

		endif

	icount = 0
20	continue
	icount = icount + 1
	if(icount.gt.20)damp = 0.5d0
	if(icount.gt.50)damp = 0.2d0
	if(icount.gt.75)damp = 0.1d0
	if(icount.gt.100)then
		err = 1
		write(75,*)
		write(75,*)' ****--->>> FAILURE TO CONVERGE in Grow3Nodes (Sub MoveNode - line 514)'
		write(12,*)' ****--->>> FAILURE TO CONVERGE in Grow3Nodes (Sub MoveNode - line 514)'
		if(debug.eq.1)then
			write(12,*)
			write(12,*)' ****--->>> FAILURE TO CONVERGE in Grow3Nodes (Sub MoveNode - line 514)'
			endif
		return		!Bail out
		endif
!	write(12,*)'------------------- icount ',icount
	! Change the position of the node and see how it affects the 3 areas
	Nnew(1) = Nnew(1) + delta
	Nnew(2) = Nnew(2)
	call GetAreas(P1,P2,P3,Nnew,Area12NdXplus,Area23NdXplus,Area13NdXplus,Asign12,Asign23,Asign13,debug)
	Nnew(1) = Nnew(1) - 2.0d0*delta
!	Nnew(2) = Nnew(2)
	call GetAreas(P1,P2,P3,Nnew,Area12NdXminus,Area23NdXminus,Area13NdXminus,Asign12,Asign23,Asign13,debug)
	
	Nnew(1) = Nnew(1) + delta	! set back to the original
	Nnew(2) = Nnew(2) + delta
	call GetAreas(P1,P2,P3,Nnew,Area12NdYplus,Area23NdYplus,Area13NdYplus,Asign12,Asign23,Asign13,debug)
!	Nnew(1) = Nnew(1)
	Nnew(2) = Nnew(2) - 2.0d0*delta
	call GetAreas(P1,P2,P3,Nnew,Area12NdYminus,Area23NdYminus,Area13NdYminus,Asign12,Asign23,Asign13,debug)
	Nnew(2) = Nnew(2) + delta	!set back to original


	! Functions to solve using Newton's method to get the new node position:
	! 0 = ∆MPh1 - ∆A12N = ∆MPh1 - (A12Nnew - A12Nold) = ∆MPh1 + A12Nold - A12Nnew =  ∆MPh1 + A12Nold - A12Nnew + dF12N/dX*∆X + dF12N/dY*∆Y
	! 0 = ∆MPh3 - ∆A13N = ∆MPh3 - (A13Nnew - A13Nold) = ∆MPh3 + A13Nold - A13Nnew =  ∆MPh3 + A13Nold - A13Nnew + dF13N/dX*∆X + dF13N/dY*∆Y

	! Calculate functions
	F12N =  Target1 - Area12Nnew
	F23N =  Target2 - Area23Nnew
	F13N =  Target3 - Area13Nnew
	
	! calculate derivatives
	dF12NdX =  -(Area12NdXplus - Area12NdXminus)/(2.0d0*delta)
	dF12NdY =  -(Area12NdYplus - Area12NdYminus)/(2.0d0*delta)
	dF23NdX =  -(Area23NdXplus - Area23NdXminus)/(2.0d0*delta)
	dF23NdY =  -(Area23NdYplus - Area23NdYminus)/(2.0d0*delta)
	dF13NdX =  -(Area13NdXplus - Area13NdXminus)/(2.0d0*delta)
	dF13NdY =  -(Area13NdYplus - Area13NdYminus)/(2.0d0*delta)

	!solve for ∆X and ∆Y using Newton's method
	select case (pass)
		case(1)
			A(1,1) = dF12NdX
			A(1,2) = dF12NdY
			A(1,3) = -F12N
			A(2,1) = dF13NdX
			A(2,2) = dF13NdY
			A(2,3) = -F13N
		case(2)
			A(1,1) = dF12NdX
			A(1,2) = dF12NdY
			A(1,3) = -F12N
			A(2,1) = dF23NdX
			A(2,2) = dF23NdY
			A(2,3) = -F23N
		case(3)
			A(1,1) = dF13NdX
			A(1,2) = dF13NdY
			A(1,3) = -F13N
			A(2,1) = dF23NdX
			A(2,2) = dF23NdY
			A(2,3) = -F23N
		case default
		end select

	do j=1,3
		do i=1,2
			AA(i,j) = A(i,j)
			end do
		end do
	!Find solution
	nsolve = 3
	DETA=0.D0
	IER=0
	numEq = 2
	CALL REDUCE (numEq,NSOLVE,DETA,IER)
	IF (IER.GT.0) then
		write(12,*)' Matrix singular in Subroutine MoveNode (line 590) in file GB_GrowNodes.f90: iNode= ',iNode
		write(75,*)' Matrix singular in Subroutine MoveNode (line 590) in file GB_GrowNodes.f90: iNode= ',iNode
		err = 1
		return
! 	      write(12,*)' ************ ERROR **************************'
! 	      write(12,*)' Matrix failed to invert in SUBROUTINE REDUCE'
! 	      write(12,*)' We are in Subroutine MDFNodesWith3Phases in file GB_MDF_Routines.f90 '
! 	      izero=1
! !	      write(12,*) 'Hit return to continue...'
! 	      pause 'Hit return to continue...'
!		      return
	      endif
	output = 0
	if(debug.eq.1)then
		write(12,*)' ----'
		write(12,*)'counter =  ',icount
		write(12,*)' '
		write(12,*)'A matrix'
		write(12,*)'    dF_dDeltaX          dF_dDeltaY       -Fo'
		do i = 1,2
			write(12,1580)(AA(i,j),j=1,3)
1580			format(50F15.5)
			end do
		write(12,*)' '
		write(12,*)'Solution xx'
		write(12,*)'    deltaX           deltaY '
		write(12,1581)(xx(j,1),j=1,numEq)
1581		format(2E15.5)
		endif

	! Calculate new values of the functions
	Nnew(1) = Nnew(1) + damp*xx(1,1)
	Nnew(2) = Nnew(2) + damp*xx(2,1)
	call GetAreas(P1,P2,P3,Nnew,Area12Nnew,Area23Nnew,Area13Nnew,Asign12,Asign23,Asign13,debug)
! 	call GetSign(Area123,Area12Nnew,Area23Nnew,Area13Nnew)		! returns the signed area. Area<0 if N is outside of the triangle
	! check to see if the node has moved outside of the original triangle and adjust if necessary
	AreaSum = Area12Nnew+Area23Nnew+Area13Nnew

	F12N =  Target1 - Area12Nnew
	F23N =  Target2 - Area23Nnew
	F13N =  Target3 - Area13Nnew
	if(debug.eq.1)then
		write(12,*)'------------------- icount ',icount
		write(12,*)' Nnew     ',Nnew(1),Nnew(2)
		write(12,*)'Area new  ',Area12Nnew,Area23Nnew,Area13Nnew
		write(12,*)'F12N =    ',F12N
		write(12,*)'F23N =    ',F23N
		write(12,*)'F13N =    ',F13N
		pause 'hit return to continue'
		endif	

	select case (pass)
	case(1)		! Use areas 12 and 13
		if(abs(F12N).lt.tol.and.abs(F13N).lt.tol)then
			R2 = (F12N)**2 + (F23N)**2 + (F13N)**2
			if(debug.eq.1)then	
				write(12,*)' Pass = (Fn12&Fn13)',pass
				write(12,*)'------------------- icount ',icount
				write(12,*)' Nnew     ',Nnew(1),Nnew(2)
				write(12,*)'Area new  ',Area12Nnew,Area23Nnew,Area13Nnew
				write(12,*)'F12N =    ',F12N
				write(12,*)'F23N =    ',F23N
				write(12,*)'F13N =    ',F13N
				write(12,*)'R2   =    ',R2
				endif
			if(R2.LT.1D-10)return		! all 3 triangles are mathed
			!pause 'hit return to continue'
			pass = pass + 1
			!if(pass.gt.3)go to 99		! do another pass
			go to 10
			endif
	case (2)		! Use areas 12 and 23
		if(abs(F12N).lt.tol.and.abs(F23N).lt.tol)then
			R2 = (F12N)**2 + (F23N)**2 + (F13N)**2
			if(debug.eq.1)then	
				write(12,*)' Pass = (Fn12&Fn23)',pass
				write(12,*)'------------------- icount ',icount
				write(12,*)' Nnew     ',Nnew(1),Nnew(2)
				write(12,*)'Area new  ',Area12Nnew,Area23Nnew,Area13Nnew
				write(12,*)'F12N =    ',F12N
				write(12,*)'F23N =    ',F23N
				write(12,*)'F13N =    ',F13N
				write(12,*)'R2   =    ',R2
				endif
			if(R2.LT.1D-10)return		! all 3 triangles are matched
			!pause 'hit return to continue'
			pass = pass + 1
			!if(pass.gt.3)go to 99		! do another pass
			go to 10
			endif
	case (3)		! Use areas 13 and 23
		if(abs(F13N).lt.tol.and.abs(F23N).lt.tol)then
			R2 = (F12N)**2 + (F23N)**2 + (F13N)**2
			if(debug.eq.1)then	
				write(12,*)' Pass = (Fn13&Fn23)',pass
				write(12,*)'------------------- icount ',icount
				write(12,*)' Nnew     ',Nnew(1),Nnew(2)
				write(12,*)'Area new  ',Area12Nnew,Area23Nnew,Area13Nnew
				write(12,*)'F12N =    ',F12N
				write(12,*)'F23N =    ',F23N
				write(12,*)'F13N =    ',F13N
				write(12,*)'R2   =    ',R2
				endif
			!pause 'hit return to continue'
			go to 99
			endif
	case default
	end select
		
	go to 20		! do another iteration

99	continue
!	pause ' hit return to end'

	end
	
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine MoveNode2(iNode,P1,P2,P3,Nold,MPh1T,MPh2T,MPh3T,Nnew,debug,err)
	! Routine to calculate the new position of a node located inside of a triangle
	! given the amount the areas of the 3 sub triangles change
	use MatrixArrays
	implicit none	
	External Area
	real*8 Area
	integer*4 icount,i,j,nsolve,numEq,ier,pass,output,debug,err,iNode,PlotsOutside,outsideTriangle
	real*8 P1(2),P2(2),P3(2)					! 3 points that define the outer triangle
	real*8 Nold(2),Nnew(2)						! node in the center of the triangle
	real*8 MPh1,MPh2,MPh3,MPh1T,MPh2T,MPh3T
	real*8 Area123,Area12Nold,Area13Nold,Area23Nold,Area12Nnew,Area13Nnew,Area23Nnew
	real*8 Area12NdXplus, Area13NdXplus, Area23NdXplus, Area12NdYplus, Area13NdYplus, Area23NdYplus, AreaSum
	real*8 Area12NdXminus,Area13NdXminus,Area23NdXminus,Area12NdYminus,Area13NdYminus,Area23NdYminus
	real*8 F12N,F13N,F23N,dF12NdX,dF12NdY,dF13NdX,dF13NdY,dF23NdX,dF23NdY,Target1,Target2,Target3
	real*8 delta,deta,tol,damp,R2
	real*8 Asign12,Asign23,Asign13
	
	tol = 0.001d0	! tolerance for convergence of areas. The zero functions are Target - area. So units are in areas, which are on the order of 1-10
	damp = 1.0d0
	delta = 0.0001d0

	! Assign moles 1,2,3 to the appropriate triangles
	! this is the correct assignment if....
	! the 3 phases are always anticlockwise around the node in the same order as the segments.
	! There is code in Subroutine Grow3Nodes that is supposed to order these correctly.
	MPh1 = MPh1T	! area12
	MPh2 = MPh2T	! area23
	MPh3 = MPh3T	! area13
	! Calculate the initial (original) areas of the 3 subtriangles
	Area123 = Area(P1,P2,P3)
	! These are initially set = 1 because the first call to Get Areas we need only the area of each triangle
	Asign12 = 1.0d0
	Asign23 = 1.0d0
	Asign13 = 1.0d0
	! on this first pass, all signs = 1 so all areas are returned as positive values
	call GetAreas(P1,P2,P3,Nold,Area12Nold,Area23Nold,Area13Nold,Asign12,Asign23,Asign13,debug)
	AreaSum = Area12Nold+Area23Nold+Area13Nold
	Target1 = MPh1 + Area12Nold
	Target2 = MPh2 + Area23Nold
	Target3 = Mph3 + Area13Nold
	if(AreaSum.gt.Area123)then
		PlotsOutside = 1
		! figure out which side of triangle 123 the node sits outside of
		if(Area12Nold+Area23Nold.ge.AreaSum)then
			outsideTriangle = 13
			go to 16
			endif
		if(Area12Nold+Area13Nold.ge.AreaSum)then
			outsideTriangle = 23
			go to 16
			endif
		if(Area13Nold+Area23Nold.ge.AreaSum)then
			outsideTriangle = 12
			go to 16
			endif
		else
		PlotsOutside = 0
		outsideTriangle = 0
		endif
16 	continue
		
	! use the largest 2 triangles for the solution (only 2 are independent
	! The hope here is that if the node is outside or if it might move outside, then the 2 largest triangles won't change sign
	! Note that if the node moves from inside to outside or outside to inside, the sign of the area for this node+seg must change
	! this creates a lot of extra code to keep track of.
	if(Area23Nold.le.Area12Nold.and.Area23Nold.le.Area13Nold)then
		! Area23N is the smallest. Use the other Area12 and Area13
		pass = 1
		go to 15
		endif
	if(Area13Nold.le.Area23Nold.and.Area13Nold.le.Area12Nold)then
		! Area13N is the smallest. Use the other Area12 and Area23
		pass = 2
		go to 15
		endif
	if(Area12Nold.le.Area23Nold.and.Area12Nold.le.Area13Nold)then
		! Area12N is the smallest. Use the other Area13 and Area23
		pass = 3
		go to 15
		endif

15	continue
	if(debug.eq.1)then
		write(12,*)' '
		write(12,*)' '
		write(12,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
		write(12,*)' Set up before starting iterations'
 		write(12,*)' Pass = ',pass
		select case (pass)
			case(1)
				write(12,*)'Functions = F12 & F13'
			case(2)
				write(12,*)'Functions = F12 & F23'
			case(3)
				write(12,*)'Functions = F13 & F23'
			case default
				write(12,*)' Pass not set properly'
			end select
		write(12,*)' Nold          ',Nold(1),Nold(2)
		write(12,*)'Area123        ',Area123
		write(12,*)'AreaSum        ',AreaSum
		write(12,*)'AreaSum-Area123',AreaSum-Area123

		if(AreaSum.gt.Area123)then
			write(12,*)'Node plots outside of triangle'
			write(12,*)' OutsideTriangle = ',outsideTriangle
			select case(outsideTriangle)
				case(13)
					write(12,*)' Outside triangle = 13N'
				case(23)
					write(12,*)' Outside triangle = 23N'
				case(12)
					write(12,*)' Outside triangle = 12N'
				case default
					Write(12,*)' Outside triangle FAILURE'
				end select				
			else
			write(12,*)'Node plots inside of triangle'
			endif

		write(12,*)'              Area 12N             Area 23N              Area 13N'
		write(12,*)'Old_area  ',Area12Nold,Area23Nold,Area13Nold
		write(12,*)'AreaDelta ',MPh1,MPh2,MPh3
		write(12,*)'Targets   ',Target1,Target2,Target3
		write(12,*)' ---------------------------------------------------'
		write(12,*)' '
		endif

	Nnew(1) = Nold(1)
	Nnew(2) = Nold(2)


! 10	continue
	call GetAreas(P1,P2,P3,Nnew,Area12Nnew,Area23Nnew,Area13Nnew,Asign12,Asign23,Asign13,debug)


	icount = 0
20	continue		! Iteration on Newton's method loops here
	icount = icount + 1
	if(icount.gt.20)damp = 0.5d0
	if(icount.gt.50)damp = 0.2d0
	if(icount.gt.75)damp = 0.1d0
	if(icount.gt.100)then
		err = 1
		write(75,*)
		write(75,*)' ****--->>> FAILURE TO CONVERGE in Grow3Nodes (Sub MoveNode - line 514)'
		write(12,*)' ****--->>> FAILURE TO CONVERGE in Grow3Nodes (Sub MoveNode - line 514)'
		if(debug.eq.1)then
			write(12,*)
			write(12,*)' ****--->>> FAILURE TO CONVERGE in Grow3Nodes (Sub MoveNode - line 514)'
			endif
		return		!Bail out
		endif

	! Calculate functions
	F12N =  Target1 - Area12Nnew
	F23N =  Target2 - Area23Nnew
	F13N =  Target3 - Area13Nnew
	
	! Change the position of the node and see how it affects the 3 areas
	! This is for calculation of derivatives
	Nnew(1) = Nnew(1) + delta
	Nnew(2) = Nnew(2)
	call GetAreas(P1,P2,P3,Nnew,Area12NdXplus,Area23NdXplus,Area13NdXplus,Asign12,Asign23,Asign13,debug)
	Nnew(1) = Nnew(1) - 2.0d0*delta
!	Nnew(2) = Nnew(2)
	call GetAreas(P1,P2,P3,Nnew,Area12NdXminus,Area23NdXminus,Area13NdXminus,Asign12,Asign23,Asign13,debug)
	
	Nnew(1) = Nnew(1) + delta	! set back to the original
	Nnew(2) = Nnew(2) + delta
	call GetAreas(P1,P2,P3,Nnew,Area12NdYplus,Area23NdYplus,Area13NdYplus,Asign12,Asign23,Asign13,debug)
!	Nnew(1) = Nnew(1)
	Nnew(2) = Nnew(2) - 2.0d0*delta
	call GetAreas(P1,P2,P3,Nnew,Area12NdYminus,Area23NdYminus,Area13NdYminus,Asign12,Asign23,Asign13,debug)
	Nnew(2) = Nnew(2) + delta	!set back to original


	! Functions to solve using Newton's method to get the new node position:
	! 0 = ∆MPh1 - ∆A12N = ∆MPh1 - (A12Nnew - A12Nold) = ∆MPh1 + A12Nold - A12Nnew =  ∆MPh1 + A12Nold - A12Nnew + dF12N/dX*∆X + dF12N/dY*∆Y
	! 0 = ∆MPh3 - ∆A13N = ∆MPh3 - (A13Nnew - A13Nold) = ∆MPh3 + A13Nold - A13Nnew =  ∆MPh3 + A13Nold - A13Nnew + dF13N/dX*∆X + dF13N/dY*∆Y

	! calculate derivatives
	dF12NdX =  -(Area12NdXplus - Area12NdXminus)/(2.0d0*delta)
	dF12NdY =  -(Area12NdYplus - Area12NdYminus)/(2.0d0*delta)
	dF23NdX =  -(Area23NdXplus - Area23NdXminus)/(2.0d0*delta)
	dF23NdY =  -(Area23NdYplus - Area23NdYminus)/(2.0d0*delta)
	dF13NdX =  -(Area13NdXplus - Area13NdXminus)/(2.0d0*delta)
	dF13NdY =  -(Area13NdYplus - Area13NdYminus)/(2.0d0*delta)

	!solve for ∆X and ∆Y using Newton's method
	select case (pass)
		case(1)
			A(1,1) = dF12NdX
			A(1,2) = dF12NdY
			A(1,3) = -F12N
			A(2,1) = dF13NdX
			A(2,2) = dF13NdY
			A(2,3) = -F13N
		case(2)
			A(1,1) = dF12NdX
			A(1,2) = dF12NdY
			A(1,3) = -F12N
			A(2,1) = dF23NdX
			A(2,2) = dF23NdY
			A(2,3) = -F23N
		case(3)
			A(1,1) = dF13NdX
			A(1,2) = dF13NdY
			A(1,3) = -F13N
			A(2,1) = dF23NdX
			A(2,2) = dF23NdY
			A(2,3) = -F23N
		case default
		end select

	do j=1,3
		do i=1,2
			AA(i,j) = A(i,j)
			end do
		end do
	!Find solution
	nsolve = 3
	DETA=0.D0
	IER=0
	numEq = 2
	CALL REDUCE (numEq,NSOLVE,DETA,IER)
	IF (IER.GT.0) then
		write(12,*)' Matrix singular in Subroutine MoveNode (line 590) in file GB_GrowNodes.f90: iNode= ',iNode
		write(75,*)' Matrix singular in Subroutine MoveNode (line 590) in file GB_GrowNodes.f90: iNode= ',iNode
		err = 1
		return
! 	      write(12,*)' ************ ERROR **************************'
! 	      write(12,*)' Matrix failed to invert in SUBROUTINE REDUCE'
! 	      write(12,*)' We are in Subroutine MDFNodesWith3Phases in file GB_MDF_Routines.f90 '
! 	      izero=1
! !	      write(12,*) 'Hit return to continue...'
! 	      pause 'Hit return to continue...'
!		      return
	      endif
	output = 0
	if(debug.eq.1)then
		write(12,*)' ----'
		write(12,*)'counter =  ',icount
		write(12,*)' '
		write(12,*)'A matrix'
		write(12,*)'    dF_dDeltaX          dF_dDeltaY       -Fo'
		do i = 1,2
			write(12,1580)(AA(i,j),j=1,3)
1580			format(50F15.5)
			end do
		write(12,*)' '
		write(12,*)'Solution xx'
		write(12,*)'    deltaX           deltaY '
		write(12,1581)(xx(j,1),j=1,numEq)
1581		format(2E15.5)
		write(12,*)' Damp = ',damp
		endif

	! Calculate new values of the functions
	Nnew(1) = Nnew(1) + damp*xx(1,1)
	Nnew(2) = Nnew(2) + damp*xx(2,1)
	call GetAreas(P1,P2,P3,Nnew,Area12Nnew,Area23Nnew,Area13Nnew,Asign12,Asign23,Asign13,debug)
! 	call GetSign(Area123,Area12Nnew,Area23Nnew,Area13Nnew)		! returns the signed area. Area<0 if N is outside of the triangle
	! check to see if the node has moved outside of the original triangle and adjust if necessary
	AreaSum = Area12Nnew+Area23Nnew+Area13Nnew

	F12N =  Target1 - Area12Nnew
	F23N =  Target2 - Area23Nnew
	F13N =  Target3 - Area13Nnew
	if(debug.eq.1)then
		write(12,*)'------------------- icount ',icount
		write(12,*)' Nnew     ',Nnew(1),Nnew(2)
		write(12,*)'Area new  ',Area12Nnew,Area23Nnew,Area13Nnew
		write(12,*)'F12N =    ',F12N
		write(12,*)'F23N =    ',F23N
		write(12,*)'F13N =    ',F13N
		pause 'hit return to continue'
		endif	

	select case (pass)
	case(1)
		if(abs(F12N).lt.tol.and.abs(F13N).lt.tol)then
			R2 = (F12N)**2 + (F23N)**2 + (F13N)**2
			if(debug.eq.1)then	
				write(12,*)' Pass = (Fn12&Fn13)',pass
				write(12,*)'------------------- icount ',icount
				write(12,*)' Nnew     ',Nnew(1),Nnew(2)
				write(12,*)'Area new  ',Area12Nnew,Area23Nnew,Area13Nnew
				write(12,*)'F12N =    ',F12N
				write(12,*)'F23N =    ',F23N
				write(12,*)'F13N =    ',F13N
				write(12,*)'R2   =    ',R2
				endif
			return
! 			if(R2.LT.1D-10)return		! all 3 triangles are matched
			!pause 'hit return to continue'

			endif
	case (2)
		if(abs(F12N).lt.tol.and.abs(F23N).lt.tol)then
			R2 = (F12N)**2 + (F23N)**2 + (F13N)**2
			if(debug.eq.1)then	
				write(12,*)' Pass = (Fn12&Fn23)',pass
				write(12,*)'------------------- icount ',icount
				write(12,*)' Nnew     ',Nnew(1),Nnew(2)
				write(12,*)'Area new  ',Area12Nnew,Area23Nnew,Area13Nnew
				write(12,*)'F12N =    ',F12N
				write(12,*)'F23N =    ',F23N
				write(12,*)'F13N =    ',F13N
				write(12,*)'R2   =    ',R2
				endif
			return
! 			if(R2.LT.1D-10)return		! all 3 triangles are matched
			!pause 'hit return to continue'
			endif
	case (3)
		if(abs(F13N).lt.tol.and.abs(F23N).lt.tol)then
			R2 = (F12N)**2 + (F23N)**2 + (F13N)**2
			if(debug.eq.1)then	
				write(12,*)' Pass = (Fn13&Fn23)',pass
				write(12,*)'------------------- icount ',icount
				write(12,*)' Nnew     ',Nnew(1),Nnew(2)
				write(12,*)'Area new  ',Area12Nnew,Area23Nnew,Area13Nnew
				write(12,*)'F12N =    ',F12N
				write(12,*)'F23N =    ',F23N
				write(12,*)'F13N =    ',F13N
				write(12,*)'R2   =    ',R2
				endif
			return
! 			if(R2.LT.1D-10)return		! all 3 triangles are matched
			!pause 'hit return to continue'
			endif
	case default
	end select
		
	go to 20		! do another iteration

99	continue
!	pause ' hit return to end'

	end
	
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Real*8 Function Area(P1,P2,P3)
! 	Routine to calculate the area of a triangle using Heron's formula
!	The area of each triangle is calculated using Heron's formula:
!
!	A = squareRoot (s*(s-a)*(s-b)*(s-c)) where
!	s = (a+b+c)/2
! 	the sides of the triangle are a, b, c
	implicit none
	real*8 P1(2),P2(2),P3(2),a,b,c,s
	a = sqrt((P1(1)-P2(1))**2 + (P1(2)-P2(2))**2)
	b = sqrt((P1(1)-P3(1))**2 + (P1(2)-P3(2))**2)
	c = sqrt((P2(1)-P3(1))**2 + (P2(2)-P3(2))**2)

	s = (a+b+c)/2.0d0
	Area = sqrt(s*(s-a)*(s-b)*(s-c))
	end function Area
	

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine GetAreas(P1,P2,P3,N,A12,A23,A13,Asign12,Asign23,Asign13,debug)
	! Routine needs P1, P2,P3,N and returns A12,A23,A13 signed	
	implicit none
	External Area
	real*8 Area
	real*8 A12,A13,A23,A12T,A13T,A23T
	real*8 P1(2),P2(2),P3(2),N(2)			
	real*8 Asign12,Asign23,Asign13
	integer*4 debug

	A12T = Area(P1,P2,N)
	A23T = Area(P2,P3,N)
	A13T = Area(P1,P3,N)
	A12 = A12T*Asign12
	A23 = A23T*Asign23
	A13 = A13T*Asign13
			
	if(debug.eq.1)then
		write(12,*)'In routine GetAreas '
		write(12,*)'Areas(unsigned)',A12T,A23T,A13T
		write(12,*)'Areas(signed)  ',A12,A23,A13
		endif
	return
	end			
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine FindIntersection(P1,P2,P3,N,Asign12,Asign23,Asign13,debug)
	! Routine to determine which line (P1-P2, P1-P3 or P2-P3 is intersected by the line P?-N
	Implicit none
	
	real*8 P1(2),P2(2),P3(2),N(2)			
	real*8 P1_P2_Slope,P1_P2_Intcpt,P2_P3_Slope,P2_P3_Intcpt,P1_P3_Slope,P1_P3_Intcpt
	real*8 P1_N_Slope,P1_N_intcpt,P2_N_Slope,P2_N_intcpt,P3_N_Slope,P3_N_intcpt
	real*8 P1_P2_X,P1_P2_Y,P2_P3_X,P2_P3_Y,P1_P3_X,P1_P3_Y
	real*8 P1_P2_Diff_X,P1_P2_Diff_Y,P2_P3_Diff_X,P2_P3_Diff_Y,P1_P3_Diff_X,P1_P3_Diff_Y
	real*8 Asign12,Asign23,Asign13
	integer*4 debug
	
	P1_P2_Slope = (P1(2)-P2(2))/(P1(1)-P2(1))
	P2_P3_Slope = (P2(2)-P3(2))/(P2(1)-P3(1))
	P1_P3_Slope = (P1(2)-P3(2))/(P1(1)-P3(1))
	P1_P2_Intcpt = P1(2) - P1_P2_Slope * P1(1)
	P2_P3_Intcpt = P2(2) - P2_P3_Slope * P2(1)
	P1_P3_Intcpt = P1(2) - P1_P3_Slope * P1(1)

	P1_N_Slope = (P1(2)-N(2))/(P1(1)-N(1))
	P2_N_Slope = (P2(2)-N(2))/(P2(1)-N(1))
	P3_N_Slope = (P3(2)-N(2))/(P3(1)-N(1))
	P1_N_Intcpt = P1(2) - P1_N_Slope * P1(1)
	P2_N_Intcpt = P2(2) - P2_N_Slope * P2(1)
	P3_N_Intcpt = P3(2) - P3_N_Slope * P3(1)

	
	!Solve for intersections
	! only need to solve for P1P2-P3N, P2P3-P1N, and P1P3-P2N
	P1_P2_X = -(P3_N_Intcpt - P1_P2_Intcpt)/(P3_N_Slope - P1_P2_Slope)
	P1_P2_Y = P1_P2_Slope * P1_P2_X + P1_P2_Intcpt
	P2_P3_X = -(P1_N_Intcpt - P2_P3_Intcpt)/(P1_N_Slope - P2_P3_Slope)
	P2_P3_Y = P2_P3_Slope * P2_P3_X + P2_P3_Intcpt
	P1_P3_X = -(P2_N_Intcpt - P1_P3_Intcpt)/(P2_N_Slope - P1_P3_Slope)
	P1_P3_Y = P1_P3_Slope * P1_P3_X + P1_P3_Intcpt

	if(debug.eq.1)then
		write(12,*)' '
		write(12,*)' In routine FindIntersection'
		write(12,*)'            Slope           Intercept'
		write(12,*)'P1_P2 ',P1_P2_Slope,P1_P2_Intcpt
		write(12,*)'P2_P3 ',P2_P3_Slope,P2_P3_Intcpt
		write(12,*)'P1_P3 ',P1_P3_Slope,P1_P3_Intcpt
		write(12,*)'P1_N  ',P1_N_Slope,P1_N_Intcpt
		write(12,*)'P2_N  ',P2_N_Slope,P2_N_Intcpt
		write(12,*)'P3_N  ',P3_N_Slope,P3_N_Intcpt
		write(12,*)' Points            X     Y'
		write(12,*)'P1    ',P1(1),P1(2)
		write(12,*)'P2    ',P2(1),P2(2)
		write(12,*)'P3    ',P3(1),P3(2)
		write(12,*)'N     ',N(1),N(2)
	
		write(12,*)' Intersections     X        Y'
		write(12,*)'P1_P2 ',P1_P2_X,P1_P2_Y
		write(12,*)'P2_P3 ',P2_P3_X,P2_P3_Y
		write(12,*)'P1_P3 ',P1_P3_X,P1_P3_Y
		write(12,*)' '
		endif		

	! Determine whether any of the 3 points lies between the two end points
	! e.g. does P1_P2(X,Y) lie between points P1 and P2
	! subtract the 2 points from the intersection. If...
	! If product > 0 the point lies outside (both differences have the same sign)
	! If product < 0 the point lies inside (the differences have opposite signs)
	P1_P2_Diff_X = (P1_P2_X - P1(1)) * (P1_P2_X - P2(1))
	P1_P2_Diff_Y = (P1_P2_Y - P1(2)) * (P1_P2_Y - P2(2))
	if(P1_P2_Diff_X.le.0.0d0.and.P1_P2_Diff_Y.le.0.0d0)then
		if(debug.eq.1)then
			write(12,*)'---------------------'
			write(12,*)'P3_N intersects P1_P2'
			write(12,*)'---------------------'
			endif
		Asign12 = -1.0d0*Asign12	! change the sign
		endif
	P2_P3_Diff_X = (P2_P3_X - P2(1)) * (P2_P3_X - P3(1))
	P2_P3_Diff_Y = (P2_P3_Y - P2(2)) * (P2_P3_Y - P3(2))
	if(P2_P3_Diff_X.le.0.0d0.and.P2_P3_Diff_Y.le.0.0d0)then
		if(debug.eq.1)then
			write(12,*)'---------------------'
			write(12,*)'P1_N intersects P2_P3'
			write(12,*)'---------------------'
			endif
		Asign23 = -1.0d0*Asign23	! change the sign
		endif
	P1_P3_Diff_X = (P1_P3_X - P1(1)) * (P1_P3_X - P3(1))
	P1_P3_Diff_Y = (P1_P3_Y - P1(2)) * (P1_P3_Y - P3(2))
	if(P1_P3_Diff_X.le.0.0d0.and.P1_P3_Diff_Y.le.0.0d0)then
		if(debug.eq.1)then
			write(12,*)'---------------------'
			write(12,*)'P2_N intersects P1_P3'
			write(12,*)'---------------------'
			endif
		Asign13 = -1.0d0*Asign13	! change the sign
		endif

	return
	end		


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine MoveOut(iNode,P1,P2,P3,Nold,abit,Nnew,debug,err)
	! Routine to move the starting node outside of the triangle
	!   when the AreaOld + Mph value = Target is negative
	! Moving outside the triangle to start will enable the "outside the triangle" algorithm to work
	! In this code
	!	P1 - P2 is the line across which the node will move
	!	The node will move along the line Nold - P3
	!  Note that in reality, the node could move across any side of the triangle
	!  Rather than code each one, that will be handled by the calling routine
		! Call MoveOut(P1,P2,P3,Nold,Nnew)	Move across side P1-P2
		! Call MoveOut(P2,P3,P1,Nold,Nnew)	Move across side P2-P3
		! Call MoveOut(P1,P3,P2,Nold,Nnew)	Move across side P1-P3

	use MatrixArrays		! common block used with matrix inversion routine Sub Reduce
	Implicit none

	integer*4 i,j,counter,nsolve,numEq,ier,debug,iNode,err
	real*8 Deta
	real*8 P1(2),P2(2),P3(2),Nold(2),Nnew(2) ! Points P1,P2,P3 and Nold. Routine returns Nnew			
	real*8 P1P2_S,P1P2_I			! Slope and intersection of line P1-P2
	real*8 P3Nold_S,P3Nold_I		! Slope and intersection of line P3-Nold
	real*8 P12_X,P12_Y			! Intersection of line P1-P22 and line P3-Nold
	real*8 LP3_XY				! length of line from the intersection along P1-P2 to P3
	real*8 abit				! How far outside the triangle we move
	real*8 Target				! Final length of line Nnew-P3 = L_P3XY + abit
	real*8 F1,dF1dX,dF1dY			! First function: 0 = L_P3XY - Target = desired length of line
	real*8 F2,dF2dX,dF2dY			! Second function: 0 = NnewY - P3Nold_S * NnewX - P3Nold_I (equation of the line P3-Nold but here using Nnew)
	! Note that this is a quadradic so there are 2 roots. Hopefully, starting outside the triangle will yield the correct root	
	real*8 tol,damp
	
	tol = 1.0d-7		! tolerance for convergence	
	damp = 1.0d0		! damping the solution, if needed

	! calculate the line P1-P2	
	P1P2_S = (P1(2)-P2(2))/(P1(1)-P2(1))
	P1P2_I = P1(2) - P1P2_S * P1(1)
	! calculate the line P3-Nold
	P3Nold_S = (P3(2)-Nold(2))/(P3(1)-Nold(1))
	P3Nold_I = P3(2) - P3Nold_S * P3(1)

	
	!Solve for intersection P1-P2 and P3-Nold
	P12_X = -(P3Nold_I - P1P2_I)/(P3Nold_S - P1P2_S)
	P12_Y = P1P2_S * P12_X + P1P2_I

	if(debug.eq.1)then
		write(12,*)'            Slope           Intercept'
		write(12,*)'P1_P2 ',P1P2_S,P1P2_I
		write(12,*)'P3_N  ',P3Nold_S,P3Nold_I
		write(12,*)' Points            X     Y'
		write(12,*)'P1    ',P1(1),P1(2)
		write(12,*)'P2    ',P2(1),P2(2)
		write(12,*)'P3    ',P3(1),P3(2)
		write(12,*)'Nold  ',Nold(1),Nold(2)
		write(12,*)' Intersections     X        Y'
		write(12,*)'P1_P2 ',P12_X,P12_Y
		endif	
		
	! Calculate the length L_P3XY
	LP3_XY = sqrt((P12_X - P3(1))**2 + (P12_Y - P3(2))**2)
	Target = LP3_XY + abit
	
	! Initialize Nnew as the value of the intersection
	Nnew(1) = P12_X
	Nnew(2) = P12_Y
	
	! Loop to here
	counter = 0
	F1 = ((Nnew(1) - P3(1))**2 + (Nnew(2) - P3(2))**2) - Target**2		! note we have squared both sides
	F2 = Nnew(2) - P3Nold_S*Nnew(1) - P3Nold_I

10	continue
	counter = counter + 1
	if(counter.gt.50)then
		write(12,*)' Counter > 50'
		write(12,*)' Something is wrong'
		pause 'pausing'
		return
		endif
		
	! Functions
	dF1dX = 2.0d0*(Nnew(1) - P3(1))
	dF1dY = 2.0d0*(Nnew(2) - P3(2))
	dF2dX = -P3Nold_S
	dF2dY = 1.0d0

	A(1,1) = dF1dX
	A(1,2) = dF1dY
	A(1,3) = -F1
	A(2,1) = dF2dX
	A(2,2) = dF2dY
	A(2,3) = -F2

	do j=1,3
		do i=1,2
			AA(i,j) = A(i,j)
			end do
		end do
	!Find solution
	nsolve = 3
	DETA=0.D0
	IER=0
	numEq = 2
	CALL REDUCE (numEq,NSOLVE,DETA,IER)
	IF (IER.GT.0) then
		write(12,*)' Matrix singular in Subroutine MoveOut(967) in file GB_GrowNodes.f90: iNode= ',iNode
		write(75,*)' Matrix singular in Subroutine MoveOut(967) in file GB_GrowNodes.f90: iNode= ',iNode
		err = 1
		return
! 	      write(12,*)' ************ ERROR **************************'
! 	      write(12,*)' Matrix failed to invert in SUBROUTINE REDUCE'
! 	      write(12,*)' We are in Subroutine MDFNodesWith3Phases in file GB_MDF_Routines.f90 '
! 	      izero=1
! !	      write(12,*) 'Hit return to continue...'
! 	      pause 'Hit return to continue...'
! !		      return
	      endif

	if(debug.eq.1)then
		write(12,*)' ----'
		write(12,*)'counter =  ',counter
		write(12,*)' '
		write(12,*)'A matrix'
		write(12,*)'    dF_dDeltaX          dF_dDeltaY       -Fo'
		do i = 1,2
			write(12,1580)(AA(i,j),j=1,3)
1580			format(50F15.5)
			end do
		write(12,*)' '
		write(12,*)'Solution xx'
		write(12,*)'    deltaX           deltaY '
		write(12,1581)(xx(j,1),j=1,numEq)
1581		format(2E15.5)
		endif

	! Calculate new values of the functions
	Nnew(1) = Nnew(1) + damp*xx(1,1)
	Nnew(2) = Nnew(2) + damp*xx(2,1)

	F1 = ((Nnew(1) - P3(1))**2 + (Nnew(2) - P3(2))**2) - Target**2		! note we have squared both sides
	F2 = Nnew(2) - P3Nold_S*Nnew(1) - P3Nold_I

	if(debug.eq.1)then
		write(12,*)'------------------- icount ',counter
		write(12,*)' Nnew     ',Nnew(1),Nnew(2)
		write(12,*)'F1 =    ',F1
		write(12,*)'F2 =    ',F2
		pause 'hit return to continue'
		endif	

	if(abs(F1).lt.tol.and.abs(F2).lt.tol)then
		if(debug.eq.1)then	
			write(12,*)' We have converged'
			endif
		return
		endif
	go to 10		! do another iteration
	

	end		


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine Grow2Nodes(iNode,moleScaleFactor,debug)
	! Routine to move nodes where 2 phases are the same
	! calling routine for Sub MoveNode to get the relevant information for a specific node
	USE AWE_INTERFACES
	
	implicit none	
	TYPE(AWE_Canvas) :: NodePlot
	include "Diffuse_GB_Hex.inc"
	include "GibbsFiles.inc"
	include "Assemb.inc"

	External Area
	real*8 Area
	integer*4 i,iNode,iSeg1,iSeg2,iSeg3,iPoint1,iPoint2,iPoint3,MIFID1,MIFID2,MIFIDuse
	integer*4 iRx1,iRx2,indexSegSame,iSegSame,icount,debug,igo,useDeltaXorY
	real*8 Nold(2),P1(2),P2(2),P3(2),MPh1,MPh2,Nnew(2),moleScaleFactor,MPhuse,Noriginal(2),delta,Ntemp(2)
	real*8 intercept,slope,AreaNNew,AreaNOld,tol,AreaNtemp,dArea,derivative,deltacalc
	real*8 AreaError,Target,Fn,AreaNoriginal

	tol = 0.0001d0		! this is the difference between the target area and the actual area. It should be scaled to the grid size

! 	X,Y of the node
	Nold(1) = NodeX(iNode)
	Nold(2) = NodeY(iNode)
	Noriginal(1) = NodeX(iNode)
	Noriginal(2) = NodeY(iNode)
	indexSegSame = nodeSegTheSamePhases(iNode)	  ! either 1, 2 or 3 -- index of the segment along which both phases are the same
	if(indexSegSame.eq.0)return			! this is the case for nodes on the boundary of the grid -- don't try to move
	iSegSame = nodeSegConnect(iNode,indexsegSame)	! segment number along which both phases are the same
	! this code simply sets the correct points -- see diagram below

	! Here are identities of the point
	! In this version we are going to use the 2 triangles
	!	Nold-P1-P3 + Nold-P2-P3
	! We do this because often Nold is nearly co-linear with P1 and P2 and as a result of moving the node
	!	it moves outside of the triangle N-P1-P2. This causes all sorts of problems
	! By using the above, we should avoid any situations where the sign of the area N-P1-P2 would have to change 
	!	(version _old_v1 tries to handle this)
	!	(version _old_v2 doesn't consider this and it causes issues for a significant number of nodes

	!
	!         P3
	!	  |
	!	  |Nold
	!	/  \
	!     /     \
	!   P1-------P2
	!

	! Code to get points P1, P2 and P3
	select case(indexSegSame)
		case(1)
			iSeg1 = nodeSegConnect(iNode,2)
			iPoint1 = nodePointNextTo(iNode,2)
			P1(1) = PointX(iSeg1,iPoint1)
			P1(2) = PointY(iSeg1,iPoint1)
			iSeg2 = nodeSegConnect(iNode,3)
			iPoint2 = nodePointNextTo(iNode,3)
			P2(1) = PointX(iSeg2,iPoint2)
			P2(2) = PointY(iSeg2,iPoint2)
		 	iseg3 = nodeSegConnect(iNode,1)
			iPoint3 = nodePointNextTo(iNode,1)
			P3(1) = PointX(iSeg3,iPoint3)
			P3(2) = PointY(iSeg3,iPoint3)
			!MIFIDuse = nodeMIFID(iNode,2)	! we use the phase that is opposite the segment where the 2 phases are the same to do the grow

		case(2)
			iSeg1 = nodeSegConnect(iNode,1)
			iPoint1 = nodePointNextTo(iNode,1)
			P1(1) = PointX(iSeg1,iPoint1)
			P1(2) = PointY(iSeg1,iPoint1)
			iSeg2 = nodeSegConnect(iNode,3)
			iPoint2 = nodePointNextTo(iNode,3)
			P2(1) = PointX(iSeg2,iPoint2)
			P2(2) = PointY(iSeg2,iPoint2)
		 	iseg3 = nodeSegConnect(iNode,2)
			iPoint3 = nodePointNextTo(iNode,2)
			P3(1) = PointX(iSeg3,iPoint3)
			P3(2) = PointY(iSeg3,iPoint3)
			!MIFIDuse = nodeMIFID(iNode,3)	! we use the phase that is opposite the segment where the 2 phases are the same to do the grow

		case(3)
			iSeg1 = nodeSegConnect(iNode,1)
			iPoint1 = nodePointNextTo(iNode,1)
			P1(1) = PointX(iSeg1,iPoint1)
			P1(2) = PointY(iSeg1,iPoint1)
			iSeg2 = nodeSegConnect(iNode,2)
			iPoint2 = nodePointNextTo(iNode,2)
			P2(1) = PointX(iSeg2,iPoint2)
			P2(2) = PointY(iSeg2,iPoint2)
		 	iseg3 = nodeSegConnect(iNode,3)
			iPoint3 = nodePointNextTo(iNode,3)
			P3(1) = PointX(iSeg3,iPoint3)
			P3(2) = PointY(iSeg3,iPoint3)
			!MIFIDuse = nodeMIFID(iNode,1)	! we use the phase that is opposite the segment where the 2 phases are the same to do the grow

		case default
			call fss_alert('Alert',' indexSegSame is not 1, 2 or 3')
			pause ' Figure out why'
		end select			


	! iSegSame = number of the seg where both phases are the same
	! segMIFID(iSegSame,1) = segMIFID(iSegSame,2) is the MIFID of the 2 phases that are in common across iSegSame
	! We will use the change in moles of the unique phase (the one that is not the same) as the monitor of how the grid will change	
	do  i = 1,3
		if(segMIFID(iSegSame,1).ne.nodeMIFID(iNode,i))then
			MIFIDuse = nodeMIFID(iNode,i)	! find the "other" MIFID = the one we want to use
			MPhuse = nodePhaseMolesDelta(iNode,i)*moleScaleFactor
			go to 20
			endif
		end do
	! if we get here, we didn't find what we wanted
	write(12,*)' ^$^$^$^$^$^$^$^$^$^$^'
	write(12,*)' In routine Grow2Nodes'
	write(12,*)' Failure to determine the correct MIFID'


20	continue	
	! get the initial areal of the triangle N-P1-P2
	! Note that this routine always returns a positive area
!	AreaNold = Area(P1,P2,Nold)		! this is the old code (v2)
	AreaNold = Area(P1,P3,Nold) + Area(P2,P3,Nold)
	AreaNOriginal = AreaNOld
!	Target = MPhuse + AreaNOriginal		! Target is the new area that we want -- this is the old code (v2)
	Target = -MPhuse + AreaNOriginal		! Target is the new area that we want
	if(debug.eq.1)then
		write(12,*)'   '
		write(12,*)'   '
		write(12,*)'--------------------------------------'
		write(12,*)'*&*&*&*&*&*&*&*&*&*&**&*&*&*&*&**'
		write(12,*)' Test Grow2Nodes routine  '
		write(12,*)' Node  ',iNode,NodeX(iNode),NodeY(iNode)
		write(12,*)'   '
		write(12,*)' Segments and points '
		write(12,104)iSeg1,iPoint1,P1(1),P1(2)
		write(12,105)iSeg2,iPoint2,P2(1),P2(2)
		write(12,106)iSeg3,iPoint3,P3(1),P3(2)
		write(12,*)'indexSegSame = ',indexSegSame
		write(12,*)'iSegSame     = ',iSegSame
104		format('seg1      ',2I8,2f12.3)
105		format('seg2      ',2I8,2f12.3)
106		format('seg3      ',2I8,2f12.3)
		write(12,*)' '
		write(12,*)'     MIFID    PhName    Del_moles      Del_area (del_moles*scale)'
		iRx1 = nodeReactionPhases(iNode,1)
		MIFID1 = nodeMIFID(iNode,iRx1)
		MPh1 = nodePhaseMolesDelta(iNode,iRx1)
		iRx2 = nodeReactionPhases(iNode,2)
		MIFID2 = nodeMIFID(iNode,iRx2)
		MPh2 = nodePhaseMolesDelta(iNode,iRx2)
		write(12,103)MIFID1,PhName(MIFID1),nodePhaseMolesDelta(iNode,iRx1),MPh1*moleScaleFactor
		write(12,103)MIFID2,PhName(MIFID2),nodePhaseMolesDelta(iNode,iRx2),MPh2*moleScaleFactor
		write(12,*)'Sum = ',Mph1+MPh2
103		format(T8,I4,2x,A8,2E15.5,5x,I8,2E5.5)
		write(12,*)' '
		write(12,*)'-------------------- '		
		write(12,*)'MIFIDuse  ',MIFIDuse
		write(12,*)'Moles     ',MPhuse
		write(12,*)'AreaNold  ',AreaNold
		write(12,*)'Target    ',Target
!		write(12,*)'Target    ',Target1

		! plot the results if we're debugging -- I don't actually think I need to plot this one -- it seems to work OK
		 call PlotNode(NodePlot,iNode,iSeg1,iPoint1,P1,iSeg2,iPoint2,P2,iSeg3,iPoint3,P3)
		endif

	! We define Target = AreaNoriginal - MPhuse --- 
		! that is, Target is the area of the triangle P1-P3-Nnew + area P2-P3-Nnew 

	! The node point (N) is constrained to move along the line Nold-P3 (this is the segment iSegSame) 
	! because the phases on either side of this line are the same
	! Calculate equation of the line along the segment with the common phases (line Nold-P3)
	! check to be sure slope isn't vertical
	if(abs(Nold(1) - P3(1)).lt.1.d-10)then
		slope = 1.d10		! nearly vertical -- close enough
		else
		slope = (Nold(2) - P3(2))/(Nold(1) - P3(1))
		endif
	intercept = Nold(2) - slope*Nold(1)

	! note that we must be careful if the slope is nearly infinite (∆x for points P1 and P2 is zero)
	! We solve this by finding the solution by using either ∆x or ∆y depending on whether the slope is 1<slope<1
	! That is, if the slope < 1, we use ∆x
	!          if the slope > 1, we use ∆y
	if(abs(slope).le.1.0d0)then
		useDeltaXorY = 1	! use ∆X
		else
		useDeltaXorY = 2	! use ∆Y
		endif

	! Loop here for Newton's method
	icount = 1
10	continue
	! Start Newton's method
	! Get the derivative using finite difference
	! Move the node along the line Nold-P3 and see if the area gets larger or smaller
	delta = 5.d-3		! a small increment to move the node
	select case (useDeltaXorY)
		case(1)		! we are using ∆X
			Ntemp(1) = Nold(1) + delta
			Ntemp(2) = slope*Ntemp(1) + intercept
		case(2)		! we are using ∆Y
			Ntemp(2) = Nold(2) + delta
			Ntemp(1) = (Ntemp(2) - intercept)/slope
		case default
			call FSS_Alert('Alert_1', 'Something went wrong in routine Grow2Nodes Alert #1')
		end select

	! calculate new area (temporary)
!	AreaNtemp = Area(P1,P2,Ntemp)		(v2) code
	AreaNtemp = Area(P1,P3,Ntemp) + Area(P2,P3,Ntemp)

	! calculate the derivative for Newton's method
	dArea = AreaNtemp - AreaNold
	Fn = AreaNOld - Target 		! we want this to be zero
	derivative = dArea/delta	! note that delta is either ∆X or ∆Y so derivative is either dArea/dX or dArea/dY
	deltaCalc = -Fn/derivative

	select case (useDeltaXorY)
		case(1)		! we are using ∆X
			Nnew(1) = Nold(1) + deltaCalc
			Nnew(2) = slope*Nnew(1) + intercept
		case(2)		! we are using ∆Y
			Nnew(2) = Nold(2) + deltaCalc
			Nnew(1) = (Nnew(2) - intercept)/slope
		case default
			call FSS_Alert('Alert_2', 'Something went wrong in routine Grow2Nodes Alert #2')
		end select

	! Calculate new area
!	AreaNnew = Area(P1,P2,Nnew)	(v2) code
	AreaNnew = Area(P1,P3,Nnew) + Area(P2,P3,Nnew)
	
	! This new area (AreaNnew) should be what we want
	! Check the area change against the ∆moles
	AreaError = Target - AreaNNew
	
	if(debug.eq.1)then
		write(12,*)'-------------------- '		
		write(12,*)'iCount     ',iCount
		write(12,*)' Derivatives etc.'
		write(12,*)'Ntemp      ',Ntemp(1),Ntemp(2)
		write(12,*)'AreaNTemp  ',AreaNTemp
		write(12,*)'slope      ',slope
		if(useDeltaXorY.eq.1)then
			write(12,*)'Using delta X'
			else
			write(12,*)'Using delta Y'
			endif			
		write(12,*)'dA,dA/dXorY',dArea,derivative
		write(12,*)'deltaCalc  ',deltaCalc
		write(12,*)'Noriginal  ',Noriginal(1),Noriginal(2)
		write(12,*)'Nold       ',Nold(1),Nold(2)
		write(12,*)'NNew       ',Nnew(1),Nnew(2)
		write(12,*)'Moles      ',MPhuse
		write(12,*)'AreaNold   ',AreaNold
		write(12,*)'AreaNnew   ',AreaNnew
		write(12,*)'TargetArea ',Target
		write(12,*)'Difference ',AreaError

		! plot the results if we're debugging -- I don't actually think I need to plot this one -- it seems to work OK
		 call PlotNewNode(NodePlot,Nnew)
		endif
	

	! check for convergence

	if(abs(AreaError).lt.tol)then
		NodeX(iNode) = Nnew(1)
		NodeY(iNode) = Nnew(2)
		if(debug.eq.1)then
			write(12,*)'---------------------------------------'
			write(12,*)'We have converged'
			write(12,*)'Moles     ',MPhuse
			write(12,*)'AreaNold  ',AreaNold
			write(12,*)'AreaNnew  ',AreaNnew
			write(12,*)'Target    ',Target
			write(12,*)'AreaError ',AreaError
			write(12,*)'Node X,Y original ',Noriginal(1),Noriginal(2)
			write(12,*)'Node X,Y    new   ',Nnew(1),Nnew(2)
			write(12,*)'iCount ',iCount
			pause 'take a look (hit return when done)'
			endif
		! I need to set segX and segY when I move the nodes
		if(segNodes(iSeg1,1).eq.iNode)then
			segX(iSeg1,1) = NodeX(iNode)
			segY(iSeg1,1) = NodeY(iNode)
			pointX(iSeg1,numPointStart(iSeg1)) = segX(iSeg1,1)
			pointY(iSeg1,numPointStart(iSeg1)) = segY(iSeg1,1)
			else
			segX(iSeg1,2) = NodeX(iNode)
			segY(iSeg1,2) = NodeY(iNode)
			pointX(iSeg1,numPointEnd(iSeg1)) = segX(iSeg1,2)
			pointY(iSeg1,numPointEnd(iSeg1)) = segY(iSeg1,2)
			endif
		if(segNodes(iSeg2,1).eq.iNode)then
			segX(iSeg2,1) = NodeX(iNode)
			segY(iSeg2,1) = NodeY(iNode)
			pointX(iSeg2,numPointStart(iSeg2)) = segX(iSeg2,1)
			pointY(iSeg2,numPointStart(iSeg2)) = segY(iSeg2,1)
			else
			segX(iSeg2,2) = NodeX(iNode)
			segY(iSeg2,2) = NodeY(iNode)
			pointX(iSeg2,numPointEnd(iSeg2)) = segX(iSeg2,2)
			pointY(iSeg2,numPointEnd(iSeg2)) = segY(iSeg2,2)
			endif
		if(segNodes(iSeg3,1).eq.iNode)then
			segX(iSeg3,1) = NodeX(iNode)
			segY(iSeg3,1) = NodeY(iNode)
			pointX(iSeg3,numPointStart(iSeg3)) = segX(iSeg3,1)
			pointY(iSeg3,numPointStart(iSeg3)) = segY(iSeg3,1)
			else
			segX(iSeg3,2) = NodeX(iNode)
			segY(iSeg3,2) = NodeY(iNode)
			pointX(iSeg3,numPointEnd(iSeg3)) = segX(iSeg3,2)
			pointY(iSeg3,numPointEnd(iSeg3)) = segY(iSeg3,2)
			endif
		return		! we have converged

		else		! do another iteration
		icount = icount + 1
		if(icount.gt.100)then
			write(75,*)'Failure in Grow2Nodes. Node = ',iNode
			!pause 'pausing'
			return		!Bail out
			endif

		if(debug.eq.1)then
			write(*,*)' Type 0 to continue with next iteration; 1 = ABORT'
			read(*,*)igo
			if(igo.eq.1)return
			endif
		Nold(1) = Nnew(1)
		Nold(2) = Nnew(2)
		AreaNOld = AreaNnew 
		go to 10
		endif


	end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine Grow2Nodes_old_v2(iNode,moleScaleFactor,debug)
	! Routine to move nodes where 2 phases are the same
	! calling routine for Sub MoveNode to get the relevant information for a specific node
	USE AWE_INTERFACES
	
	implicit none	
	TYPE(AWE_Canvas) :: NodePlot
	include "Diffuse_GB_Hex.inc"
	include "GibbsFiles.inc"
	include "Assemb.inc"

	External Area
	real*8 Area
	integer*4 i,iNode,iSeg1,iSeg2,iSeg3,iPoint1,iPoint2,iPoint3,MIFID1,MIFID2,MIFIDuse
	integer*4 iRx1,iRx2,indexSegSame,iSegSame,icount,debug,igo,useDeltaXorY
	real*8 Nold(2),P1(2),P2(2),P3(2),MPh1,MPh2,Nnew(2),moleScaleFactor,MPhuse,Noriginal(2),delta,Ntemp(2)
	real*8 intercept,slope,AreaNNew,AreaNOld,tol,AreaNtemp,dArea,derivative,deltacalc
	real*8 AreaError,Target,Fn,AreaNoriginal

	tol = 0.0001d0		! this is the difference between the target area and the actual area. It should be scaled to the grid size

! 	X,Y of the node
	Nold(1) = NodeX(iNode)
	Nold(2) = NodeY(iNode)
	Noriginal(1) = NodeX(iNode)
	Noriginal(2) = NodeY(iNode)
	indexSegSame = nodeSegTheSamePhases(iNode)	  ! either 1, 2 or 3 -- index of the segment along which both phases are the same
	if(indexSegSame.eq.0)return			! this is the case for nodes on the boundary of the grid -- don't try to move
	iSegSame = nodeSegConnect(iNode,indexsegSame)	! segment number along which both phases are the same
	! this code simply sets the correct points -- see diagram below

	! Here are identities of the point
	! It is the triangle N-P1-P2 that we monitor to change size according to the change in moles
	! of the phase in this triangle
	!
	!         P3
	!	  |
	!	  |Nold
	!	/  \
	!     /     \
	!   P1-------P2
	!

	! Code to wet points P1, P2 and P3
	select case(indexSegSame)
		case(1)
			iSeg1 = nodeSegConnect(iNode,2)
			iPoint1 = nodePointNextTo(iNode,2)
			P1(1) = PointX(iSeg1,iPoint1)
			P1(2) = PointY(iSeg1,iPoint1)
			iSeg2 = nodeSegConnect(iNode,3)
			iPoint2 = nodePointNextTo(iNode,3)
			P2(1) = PointX(iSeg2,iPoint2)
			P2(2) = PointY(iSeg2,iPoint2)
		 	iseg3 = nodeSegConnect(iNode,1)
			iPoint3 = nodePointNextTo(iNode,1)
			P3(1) = PointX(iSeg3,iPoint3)
			P3(2) = PointY(iSeg3,iPoint3)
			!MIFIDuse = nodeMIFID(iNode,2)	! we use the phase that is opposite the segment where the 2 phases are the same to do the grow

		case(2)
			iSeg1 = nodeSegConnect(iNode,1)
			iPoint1 = nodePointNextTo(iNode,1)
			P1(1) = PointX(iSeg1,iPoint1)
			P1(2) = PointY(iSeg1,iPoint1)
			iSeg2 = nodeSegConnect(iNode,3)
			iPoint2 = nodePointNextTo(iNode,3)
			P2(1) = PointX(iSeg2,iPoint2)
			P2(2) = PointY(iSeg2,iPoint2)
		 	iseg3 = nodeSegConnect(iNode,2)
			iPoint3 = nodePointNextTo(iNode,2)
			P3(1) = PointX(iSeg3,iPoint3)
			P3(2) = PointY(iSeg3,iPoint3)
			!MIFIDuse = nodeMIFID(iNode,3)	! we use the phase that is opposite the segment where the 2 phases are the same to do the grow

		case(3)
			iSeg1 = nodeSegConnect(iNode,1)
			iPoint1 = nodePointNextTo(iNode,1)
			P1(1) = PointX(iSeg1,iPoint1)
			P1(2) = PointY(iSeg1,iPoint1)
			iSeg2 = nodeSegConnect(iNode,2)
			iPoint2 = nodePointNextTo(iNode,2)
			P2(1) = PointX(iSeg2,iPoint2)
			P2(2) = PointY(iSeg2,iPoint2)
		 	iseg3 = nodeSegConnect(iNode,3)
			iPoint3 = nodePointNextTo(iNode,3)
			P3(1) = PointX(iSeg3,iPoint3)
			P3(2) = PointY(iSeg3,iPoint3)
			!MIFIDuse = nodeMIFID(iNode,1)	! we use the phase that is opposite the segment where the 2 phases are the same to do the grow

		case default
			call fss_alert('Alert',' indexSegSame is not 1, 2 or 3')
			pause ' Figure out why'
		end select			


	! iSegSame = number of the seg where both phases are the same
	! segMIFID(iSegSame,1) = segMIFID(iSegSame,2) is the MIFID of the 2 phases that are in common across iSegSame
	! We will use the change in moles of the unique phase (the one that is not the same) as the monitor of how the grid will change	
	do  i = 1,3
		if(segMIFID(iSegSame,1).ne.nodeMIFID(iNode,i))then
			MIFIDuse = nodeMIFID(iNode,i)	! find the "other" MIFID = the one we want to use
			MPhuse = nodePhaseMolesDelta(iNode,i)*moleScaleFactor
			go to 20
			endif
		end do
	! if we get here, we didn't find what we wanted
	write(12,*)' ^$^$^$^$^$^$^$^$^$^$^'
	write(12,*)' In routine Grow2Nodes'
	write(12,*)' Failure to determine the correct MIFID'


20	continue	
	! get the initial areal of the triangle N-P1-P2
	! Note that this routine always returns a positive area
	AreaNold = Area(P1,P2,Nold)
	AreaNOriginal = AreaNOld
	Target = MPhuse + AreaNOriginal		! Target is the new area that we want
	if(debug.eq.1)then
		write(12,*)'   '
		write(12,*)'   '
		write(12,*)'--------------------------------------'
		write(12,*)'*&*&*&*&*&*&*&*&*&*&**&*&*&*&*&**'
		write(12,*)' Test Grow2Nodes routine  '
		write(12,*)' Node  ',iNode,NodeX(iNode),NodeY(iNode)
		write(12,*)'   '
		write(12,*)' Segments and points '
		write(12,104)iSeg1,iPoint1,P1(1),P1(2)
		write(12,105)iSeg2,iPoint2,P2(1),P2(2)
		write(12,106)iSeg3,iPoint3,P3(1),P3(2)
		write(12,*)'indexSegSame = ',indexSegSame
		write(12,*)'iSegSame     = ',iSegSame
104		format('seg1      ',2I8,2f12.3)
105		format('seg2      ',2I8,2f12.3)
106		format('seg3      ',2I8,2f12.3)
		write(12,*)' '
		write(12,*)'     MIFID    PhName    Del_moles      Del_area (del_moles*scale)'
		iRx1 = nodeReactionPhases(iNode,1)
		MIFID1 = nodeMIFID(iNode,iRx1)
		MPh1 = nodePhaseMolesDelta(iNode,iRx1)
		iRx2 = nodeReactionPhases(iNode,2)
		MIFID2 = nodeMIFID(iNode,iRx2)
		MPh2 = nodePhaseMolesDelta(iNode,iRx2)
		write(12,103)MIFID1,PhName(MIFID1),nodePhaseMolesDelta(iNode,iRx1),MPh1*moleScaleFactor
		write(12,103)MIFID2,PhName(MIFID2),nodePhaseMolesDelta(iNode,iRx2),MPh2*moleScaleFactor
		write(12,*)'Sum = ',Mph1+MPh2
103		format(T8,I4,2x,A8,2E15.5,5x,I8,2E5.5)
		write(12,*)' '
		write(12,*)'-------------------- '		
		write(12,*)'MIFIDuse  ',MIFIDuse
		write(12,*)'Moles     ',MPhuse
		write(12,*)'AreaNold  ',AreaNold
		write(12,*)'Target    ',Target
!		write(12,*)'Target    ',Target1

		! plot the results if we're debugging -- I don't actually think I need to plot this one -- it seems to work OK
		 call PlotNode(NodePlot,iNode,iSeg1,iPoint1,P1,iSeg2,iPoint2,P2,iSeg3,iPoint3,P3)
		endif

	! Loop here for Newton's method
	icount = 1
10	continue
	! The difference in area (AreaNnew - AreaNoriginal) must = Mphuse = Target (the change in moles of a phase)
	! We define Target = AreaNoriginal + MPhuse --- 
		! that is, Target is the area of the triangle P1-P2-Nnew 

	! The node point (N) is constrained to move along the line Nold-P3 (this is the segment iSegSame) 
	! because the phases on either side of this line are the same
	! Calculate equation of the line along the segment with the common phases (line Nold-P3)
	! check to be sure slope isn't vertical
	if(abs(Nold(1) - P3(1)).lt.1.d-10)then
		slope = 1.d10		! nearly vertical -- close enough
		else
		slope = (Nold(2) - P3(2))/(Nold(1) - P3(1))
		endif
	intercept = Nold(2) - slope*Nold(1)

	! note that we must be careful if the slope is nearly infinite (∆x for points P1 and P2 is zero)
	! We solve this by finding the solution by using either ∆x or ∆y depending on whether the slope is 1<slope<1
	! That is, if the slope < 1, we use ∆x
	!          if the slope > 1, we use ∆y
	if(abs(slope).le.1.0d0)then
		useDeltaXorY = 1	! use ∆X
		else
		useDeltaXorY = 2	! use ∆Y
		endif

	! Start Newton's method
	! Get the derivative using finite difference
	! Move the node along the line Nold-P3 and see if the area gets larger or smaller
	delta = 5.d-3		! a small increment to move the node
	select case (useDeltaXorY)
		case(1)		! we are using ∆X
			Ntemp(1) = Nold(1) + delta
			Ntemp(2) = slope*Ntemp(1) + intercept
		case(2)		! we are using ∆Y
			Ntemp(2) = Nold(2) + delta
			Ntemp(1) = (Ntemp(2) - intercept)/slope
		case default
			call FSS_Alert('Alert_1', 'Something went wrong in routine Grow2Nodes Alert #1')
		end select

	! calculate new area (temporary)
	AreaNtemp = Area(P1,P2,Ntemp)

	! calculate the derivative for Newton's method
	dArea = AreaNtemp - AreaNold
	Fn = AreaNOld - Target 		! we want this to be zero
	derivative = dArea/delta	! note that delta is either ∆X or ∆Y so derivative is either dArea/dX or dArea/dY
	deltaCalc = -Fn/derivative

	select case (useDeltaXorY)
		case(1)		! we are using ∆X
			Nnew(1) = Nold(1) + deltaCalc
			Nnew(2) = slope*Nnew(1) + intercept
		case(2)		! we are using ∆Y
			Nnew(2) = Nold(2) + deltaCalc
			Nnew(1) = (Nnew(2) - intercept)/slope
		case default
			call FSS_Alert('Alert_2', 'Something went wrong in routine Grow2Nodes Alert #2')
		end select

	! Calculate new area
	AreaNnew = Area(P1,P2,Nnew)
	
	! This new area (AreaNnew) should be what we want
	! Check the area change against the ∆moles
	AreaError = Target - AreaNNew
	
	if(debug.eq.1)then
		write(12,*)'-------------------- '		
		write(12,*)'iCount     ',iCount
		write(12,*)' Derivatives etc.'
		write(12,*)'Ntemp      ',Ntemp(1),Ntemp(2)
		write(12,*)'AreaNTemp  ',AreaNTemp
		write(12,*)'slope      ',slope
		if(useDeltaXorY.eq.1)then
			write(12,*)'Using delta X'
			else
			write(12,*)'Using delta Y'
			endif			
		write(12,*)'dA,dA/dXorY',dArea,derivative
		write(12,*)'deltaCalc  ',deltaCalc
		write(12,*)'Noriginal  ',Noriginal(1),Noriginal(2)
		write(12,*)'Nold       ',Nold(1),Nold(2)
		write(12,*)'NNew       ',Nnew(1),Nnew(2)
		write(12,*)'Moles      ',MPhuse
		write(12,*)'AreaNold   ',AreaNold
		write(12,*)'AreaNnew   ',AreaNnew
		write(12,*)'TargetArea ',Target
		write(12,*)'Difference ',AreaError

		! plot the results if we're debugging -- I don't actually think I need to plot this one -- it seems to work OK
		 call PlotNewNode(NodePlot,Nnew)
		endif
	

	! check for convergence

	if(abs(AreaError).lt.tol)then
		NodeX(iNode) = Nnew(1)
		NodeY(iNode) = Nnew(2)
		if(debug.eq.1)then
			write(12,*)'---------------------------------------'
			write(12,*)'We have converged'
			write(12,*)'Moles     ',MPhuse
			write(12,*)'AreaNold  ',AreaNold
			write(12,*)'AreaNnew  ',AreaNnew
			write(12,*)'Target    ',Target
			write(12,*)'AreaError ',AreaError
			write(12,*)'Node X,Y original ',Noriginal(1),Noriginal(2)
			write(12,*)'Node X,Y    new   ',Nnew(1),Nnew(2)
			write(12,*)'iCount ',iCount
			pause 'take a look (hit return when done)'
			endif
		! I need to set segX and segY when I move the nodes
		if(segNodes(iSeg1,1).eq.iNode)then
			segX(iSeg1,1) = NodeX(iNode)
			segY(iSeg1,1) = NodeY(iNode)
			pointX(iSeg1,numPointStart(iSeg1)) = segX(iSeg1,1)
			pointY(iSeg1,numPointStart(iSeg1)) = segY(iSeg1,1)
			else
			segX(iSeg1,2) = NodeX(iNode)
			segY(iSeg1,2) = NodeY(iNode)
			pointX(iSeg1,numPointEnd(iSeg1)) = segX(iSeg1,2)
			pointY(iSeg1,numPointEnd(iSeg1)) = segY(iSeg1,2)
			endif
		if(segNodes(iSeg2,1).eq.iNode)then
			segX(iSeg2,1) = NodeX(iNode)
			segY(iSeg2,1) = NodeY(iNode)
			pointX(iSeg2,numPointStart(iSeg2)) = segX(iSeg2,1)
			pointY(iSeg2,numPointStart(iSeg2)) = segY(iSeg2,1)
			else
			segX(iSeg2,2) = NodeX(iNode)
			segY(iSeg2,2) = NodeY(iNode)
			pointX(iSeg2,numPointEnd(iSeg2)) = segX(iSeg2,2)
			pointY(iSeg2,numPointEnd(iSeg2)) = segY(iSeg2,2)
			endif
		if(segNodes(iSeg3,1).eq.iNode)then
			segX(iSeg3,1) = NodeX(iNode)
			segY(iSeg3,1) = NodeY(iNode)
			pointX(iSeg3,numPointStart(iSeg3)) = segX(iSeg3,1)
			pointY(iSeg3,numPointStart(iSeg3)) = segY(iSeg3,1)
			else
			segX(iSeg3,2) = NodeX(iNode)
			segY(iSeg3,2) = NodeY(iNode)
			pointX(iSeg3,numPointEnd(iSeg3)) = segX(iSeg3,2)
			pointY(iSeg3,numPointEnd(iSeg3)) = segY(iSeg3,2)
			endif
		return		! we have converged

		else		! do another iteration
		icount = icount + 1
		if(icount.gt.100)then
			write(75,*)'Failure in Grow2Nodes. Node = ',iNode
			!pause 'pausing'
			return		!Bail out
			endif

		if(debug.eq.1)then
			write(*,*)' Type 0 to continue with next iteration; 1 = ABORT'
			read(*,*)igo
			if(igo.eq.1)return
			endif
		Nold(1) = Nnew(1)
		Nold(2) = Nnew(2)
		AreaNOld = AreaNnew 
		go to 10
		endif


	end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine Grow2Nodes_old_v1(iNode,moleScaleFactor,debug)
	! Routine to move nodes where 2 phases are the same
	! calling routine for Sub MoveNode to get the relevant information for a specific node
	
	implicit none	
	include "Diffuse_GB_Hex.inc"
	include "GibbsFiles.inc"
	include "Assemb.inc"

	External Area
	real*8 Area
	integer*4 iNode,iSeg1,iSeg2,iSeg3,iPoint1,iPoint2,iPoint3,MIFID1,MIFID2,MIFIDuse
	integer*4 iRx1,iRx2,indexSegSame,iSegSame,icount,debug,igo
	real*8 Nold(2),P1(2),P2(2),P3(2),MPh1,MPh2,Nnew(2),moleScaleFactor,MPhuse
	real*8 target1,intercept,slope,AreaNdX,delX,AreaNNew,AreaNOld,tol,dX,Fn,FndX,dFn_dX
	real*8 slopeP1P2,InterceptP1P2,Xi,Yi,DistP3_XiYi,DistP3_Nold,abit,Dist,PI(2)

	tol = 0.0001d0		! this is the difference between the target area and the actual area. It should be scaled to the grid size

! 	X,Y of the node
	Nold(1) = NodeX(iNode)
	Nold(2) = NodeY(iNode)

!	Adjacent points
	!write(12,*)'P,  iSeg, iPoint    X    Y'
	indexSegSame = nodeSegTheSamePhases(iNode)	  ! either 1, 2 or 3 -- index of the segment along which both phases are the same
	if(indexSegSame.eq.0)return			! this is the case for nodes on the boundary of the grid -- don't try to move
	iSegSame = nodeSegConnect(iNode,indexsegSame)	! segment number along which both phases are the same
	select case(indexSegSame)
		case(1)
			iSeg1 = nodeSegConnect(iNode,2)
			iPoint1 = nodePointNextTo(iNode,2)
			P1(1) = PointX(iSeg1,iPoint1)
			P1(2) = PointY(iSeg1,iPoint1)
			iSeg2 = nodeSegConnect(iNode,3)
			iPoint2 = nodePointNextTo(iNode,3)
			P2(1) = PointX(iSeg2,iPoint2)
			P2(2) = PointY(iSeg2,iPoint2)
		 	iseg3 = nodeSegConnect(iNode,1)
			iPoint3 = nodePointNextTo(iNode,1)
			P3(1) = PointX(iSeg3,iPoint3)
			P3(2) = PointY(iSeg3,iPoint3)
			MIFIDuse = nodeMIFID(iNode,2)	! we use the phase that is opposite the segment where the 2 phases are the same to do the grow

		case(2)
			iSeg1 = nodeSegConnect(iNode,1)
			iPoint1 = nodePointNextTo(iNode,1)
			P1(1) = PointX(iSeg1,iPoint1)
			P1(2) = PointY(iSeg1,iPoint1)
			iSeg2 = nodeSegConnect(iNode,3)
			iPoint2 = nodePointNextTo(iNode,3)
			P2(1) = PointX(iSeg2,iPoint2)
			P2(2) = PointY(iSeg2,iPoint2)
		 	iseg3 = nodeSegConnect(iNode,2)
			iPoint3 = nodePointNextTo(iNode,2)
			P3(1) = PointX(iSeg3,iPoint3)
			P3(2) = PointY(iSeg3,iPoint3)
			MIFIDuse = nodeMIFID(iNode,3)	! we use the phase that is opposite the segment where the 2 phases are the same to do the grow

		case(3)
			iSeg1 = nodeSegConnect(iNode,1)
			iPoint1 = nodePointNextTo(iNode,1)
			P1(1) = PointX(iSeg1,iPoint1)
			P1(2) = PointY(iSeg1,iPoint1)
			iSeg2 = nodeSegConnect(iNode,2)
			iPoint2 = nodePointNextTo(iNode,2)
			P2(1) = PointX(iSeg2,iPoint2)
			P2(2) = PointY(iSeg2,iPoint2)
		 	iseg3 = nodeSegConnect(iNode,3)
			iPoint3 = nodePointNextTo(iNode,3)
			P3(1) = PointX(iSeg3,iPoint3)
			P3(2) = PointY(iSeg3,iPoint3)
			MIFIDuse = nodeMIFID(iNode,1)	! we use the phase that is opposite the segment where the 2 phases are the same to do the grow

		case default
			call fss_alert('Alert',' indexSegSame is not 1, 2 or 3')
			pause ' Figure out why'
		end select			


	! Get the MIF IDs so we know what the 3 phases are. 
	! These are needed so we can scale the moles to the number of oxygens
	iRx1 = nodeReactionPhases(iNode,1)
	MIFID1 = nodeMIFID(iNode,iRx1)
! 	MPh1 = nodePhaseMolesDelta(iNode,iRx1)*numOxygens(MIFID1)
	MPh1 = nodePhaseMolesDelta(iNode,iRx1)
	iRx2 = nodeReactionPhases(iNode,2)
	MIFID2 = nodeMIFID(iNode,iRx2)
! 	MPh2 = nodePhaseMolesDelta(iNode,iRx2)*numOxygens(MIFID2)
	MPh2 = nodePhaseMolesDelta(iNode,iRx2)

	! This just scales to moles -- not sure what the scaling factor should be
	MPh1 = MPh1*moleScaleFactor
	MPh2 = MPh2*moleScaleFactor
	if(MIFIDuse.eq.MIFID1)then
		MPhuse = MPh1
		else		! MIFIDuse must = MPh2
		MPhuse = MPh2
		endif

!11	continue

	AreaNold = Area(P1,P2,Nold)

	! The "Target" is the area of the new triangle AreaNnew = P1-P2-Nnew
	! Target can only be positive because the Sub Area only returns positive areas
	! There are 6 possible ways to calculate Target depending on
		! (a) The orientation of the original triangle (AreaNold = P1-P2-Nold) relative to the 
			!segment that borders the 2 phases that are the same (iSegSame)
		! (b) The magnitude of MPhuse relative to magnitude of AreaNold

	! Here are the 2 configurations 
	!(a)
	!         P3
	!	  |
	!	  |Nold
	!	/  \
	!     /     \
	!   P1-------P2
	!
	! (b)	
	!         P3
	!	  |
	!    P1---|--P2
	!      \  |  /
	!       \ | /
	!	  |Nold
	
	! We determine which configuration we have by comparing the distances between Nold-P3 to to the 
	!       distance Nold- and the point of intersection between Nold-P3 and the line P1-P2
	! Calculate equation of the line along the segment with the common phases (line Nold-P3)
!	slope = (Nold(2) - PointY(iSegSame,iPoint3))/(Nold(1) - PointX(iSegSame,iPoint3))
	slope = (Nold(2) - P3(2))/(Nold(1) - P3(1))
	intercept = Nold(2) - slope*Nold(1)
	! calculate the equation of the lline P1-P2
	if(abs(P2(1)-P1(1)).lt.1d-10)then
		slopeP1P2 = 1e-10
		interceptP1P2 = P2(2)
		else
		slopeP1P2 = (P2(2)-P1(2))/(P2(1)-P1(1))
		interceptP1P2 = P2(2) - slopeP1P2*P2(1)
		endif
	! Find the intersection of the 2 lines. Set each equation (y = mx+b) equal to each other and solve for x
	Xi = -(interceptP1P2 - intercept)/(slopeP1P2 - slope)
	Yi = slope*Xi + intercept
	! now compute distance Nold to P3 and Nold to XiYi and compare (actually just distance squared -- no need to do square root)
	DistP3_XiYi = (P3(1) - Xi)**2      + (P3(2) - Yi)**2
	DistP3_Nold = (P3(1) - Nold(1))**2 + (P3(2) - Nold(2))**2

	abit = 0.5d0			! distance to nudge the Nold over the P1-P2 line		
	if(DistP3_Nold.le.DistP3_XiYi)then		! configuration (a)
		if(MPhuse.ge.0.0d0)then		! Nold moves towards P3 and triangle gets larger
			Target1 = AreaNold + MPhuse
			go to 15
			endif
		! if here then MPhuse is negative
		if(abs(MPhuse).le.AreaNold)then		! The node is inside the original triangle AreaNold (P1-P2-Nold)
			Target1 = AreaNold + MPhuse	! MPhuse is negative so Target1 is smaller than AreaNold
			go to 15
			else
			Target1 = abs(MPhuse) - AreaNold	! Nold will be below P1-P2 (outside the original triangle)
			! The Sub Area returns a positive area and will solve for Nnew inside the original triangle
			!	unless we give Nold(prime) an initial value that is below P1-P2
			! We solve for this new point by solving the 2 equations 
			!	Distance^2 = (Xnew-Xold)^2 + (Ynew-Yold)^2
			!	Ynew = slope*Xnew + Intcpt
			! This is done in subroutine FindPoint
			Dist = DistP3_XiYi + abit
			PI(1) = Xi
			PI(2) = Yi
			call FindPoint(Dist,slope,Intercept,P3,Nnew,PI,debug)
			if(debug.eq.1)then
				write(12,*)'+++after call to FindPoint+++ Config = (a)'
				write(12,*)'P3     ',P3(1),P3(2)
				write(12,*)'P1     ',P1(1),P1(2)
				write(12,*)'P2     ',P2(1),P2(2)
				write(12,*)'PI     ',PI(1),PI(2)
				write(12,*)'Nold   ',Nold(1),Nold(2)
				write(12,*)'NNew   ',Nnew(1),Nnew(2)
				write(12,*)' DistP3_Nold    DistP3_XiYi     Dist'
				write(12,*)DistP3_Nold,DistP3_XiYi,Dist
				write(12,*)'++++++++++++++++++++++++++++'
				endif
			Nold(1) = Nnew(1)
			Nold(2) = Nnew(2)
			go to 15
			endif
		
		else		! configuration (b)

		if(MPhuse.le.0.0d0)then		! Nold moves away from P3
			Target1 = AreaNold + abs(MPhuse)
			go to 15
			endif
		! if here, then MPhuse is positive and node moves towards P3
		if(MPhuse.le.AreaNold)then		! the node moves to inside the triangle Nold-P1-P2
			Target1 = AreaNOld - MPhuse
			go to 15
			else			! The node modes outside the triangle Nold-P1-P2 (above line P1-P2
			Target1 = MPhuse - AreaNold
			! The Sub Area returns a positive area and will solve for Nnew inside the original triangle
			!	unless we give Nold(prime) an initial value that is below P1-P2
			! We solve for this new point by solving the 2 equations 
			!	Distance^2 = (Xnew-Xold)^2 + (Ynew-Yold)^2
			!	Ynew = slope*Xnew + Intcpt
			! This is done in subroutine FindPoint
			Dist = DistP3_XiYi - abit	! move a bit closer to P3
			PI(1) = Xi
			PI(2) = Yi
			call FindPoint(Dist,slope,Intercept,P3,Nnew,PI,debug)
			if(debug.eq.1)then
				write(12,*)'+++after call to FindPoint+++ Config = (b)'
				write(12,*)'P3     ',P3(1),P3(2)
				write(12,*)'Nold   ',Nold(1),Nold(2)
				write(12,*)'Nold   ',Nnew(1),Nnew(2)
				write(12,*)'PI     ',PI(1),PI(2)
				write(12,*)' DistP3_Nold    DistP3_XiYi     Dist'
				write(12,*)DistP3_Nold,DistP3_XiYi,Dist
				write(12,*)'++++++++++++++++++++++++++++'
				endif
			Nold(1) = Nnew(1)
			Nold(2) = Nnew(2)
			go to 15
			endif
		endif

15	continue
	if(target1.lt.0.0d0)then
		write(12,*)'Target1 < 0 -- this should not happen'
		pause 'Hit return to continue'
		endif

	if(debug.eq.1)then
		write(12,*)'   '
		write(12,*)'   '
		write(12,*)'--------------------------------------'
		write(12,*)'*&*&*&*&*&*&*&*&*&*&**&*&*&*&*&**'
		write(12,*)' Test Grow2Nodes routine  '
		write(12,*)' Node  ',iNode,NodeX(iNode),NodeY(iNode)
		write(12,*)'   '
		write(12,*)' Segments and points '
		write(12,104)iSeg1,iPoint1,P1(1),P1(2)
		write(12,105)iSeg2,iPoint2,P2(1),P2(2)
		write(12,*)'indexSegSame = ',indexSegSame
		write(12,*)'iSegSame     = ',iSegSame
104		format('seg1      ',2I8,2f12.3)
105		format('seg2      ',2I8,2f12.3)
		write(12,*)' '
		write(12,*)'     MIFID    PhName    Del_moles      Del_area (del_moles*scale)'
		write(12,103)MIFID1,PhName(MIFID1),nodePhaseMolesDelta(iNode,iRx1),MPh1
		write(12,103)MIFID2,PhName(MIFID2),nodePhaseMolesDelta(iNode,iRx2),MPh2
		write(12,*)'Sum = ',Mph1+MPh2
103		format(T8,I4,2x,A8,2E15.5,5x,I8,2E5.5)

		write(12,*)' '
		write(12,*)'-------------------- '		
		write(12,*)'AreaNold  ',AreaNold
		write(12,*)'Moles     ',MPh1
		write(12,*)'Target    ',Target1

		! plot the results if we're debugging -- I don't actually think I need to plot this one -- it seems to work OK
!		call PlotNode(iNode,iSeg1,iPoint1,P1,iSeg2,iPoint2,P2,iSeg3,iPoint3,P3)


		endif



	! The new node position must be along the segment iSegSame Y = intercept + slope*X
	! I think there's an analytical solution to this but looking at the formula for the area of a triangle, it would involve a lot of algebra
	! I'll use Newton's method instead
	! The function to find the root of is
	! 0 = (Target1 - Area) = Fn
	! 0 = F + (dF/dx)∆X
	dX = 0.01d0
	icount = 1
10	continue		! loop to here for Newton's method
	Fn = Target1 - AreaNold
	Nnew(1) = Nold(1) + dX
	Nnew(2) = intercept + slope*Nnew(1)
	AreaNdX = Area(P1,P2,Nnew)
	FndX = Target1 - AreaNdX
	dFn_dX = (Fn - FndX)/dX
	delX = (Fn/dFn_dX)
	Nold(1) = Nold(1) + delX
	Nold(2) = intercept + slope*Nold(1)
	AreaNnew = Area(P1,P2,Nold)

	if(abs(Target1 - AreaNnew).lt.tol)then
		NodeX(iNode) = Nold(1)
		NodeY(iNode) = Nold(2)
		if(debug.eq.1)then
			write(12,*)'---------------------------------------'
			write(12,*)'We have converged'
			write(12,*)'Node X,Y original ',nodeX(iNode),nodeY(iNode)
			write(12,*)'Node X,Y    new   ',Nold(1),Nold(2)
			write(12,*)'AreaOld     ',AreaNold
			write(12,*)'Areanew     ',AreaNnew
			write(12,*)'Targ1 - Area',Target1-AreaNnew
			write(12,*)'NewNode ',iNode,NodeX(iNode),NodeY(iNode)
			write(12,*)'iCount ',iCount
			pause 'take a look (hit return when done)'
			endif
		! I need to set segX and segY when I move the nodes
		if(segNodes(iSeg1,1).eq.iNode)then
			segX(iSeg1,1) = NodeX(iNode)
			segY(iSeg1,1) = NodeY(iNode)
			pointX(iSeg1,numPointStart(iSeg1)) = segX(iSeg1,1)
			pointY(iSeg1,numPointStart(iSeg1)) = segY(iSeg1,1)
			else
			segX(iSeg1,2) = NodeX(iNode)
			segY(iSeg1,2) = NodeY(iNode)
			pointX(iSeg1,numPointEnd(iSeg1)) = segX(iSeg1,2)
			pointY(iSeg1,numPointEnd(iSeg1)) = segY(iSeg1,2)
			endif
		if(segNodes(iSeg2,1).eq.iNode)then
			segX(iSeg2,1) = NodeX(iNode)
			segY(iSeg2,1) = NodeY(iNode)
			pointX(iSeg2,numPointStart(iSeg2)) = segX(iSeg2,1)
			pointY(iSeg2,numPointStart(iSeg2)) = segY(iSeg2,1)
			else
			segX(iSeg2,2) = NodeX(iNode)
			segY(iSeg2,2) = NodeY(iNode)
			pointX(iSeg2,numPointEnd(iSeg2)) = segX(iSeg2,2)
			pointY(iSeg2,numPointEnd(iSeg2)) = segY(iSeg2,2)
			endif
		if(segNodes(iSeg3,1).eq.iNode)then
			segX(iSeg3,1) = NodeX(iNode)
			segY(iSeg3,1) = NodeY(iNode)
			pointX(iSeg3,numPointStart(iSeg3)) = segX(iSeg3,1)
			pointY(iSeg3,numPointStart(iSeg3)) = segY(iSeg3,1)
			else
			segX(iSeg3,2) = NodeX(iNode)
			segY(iSeg3,2) = NodeY(iNode)
			pointX(iSeg3,numPointEnd(iSeg3)) = segX(iSeg3,2)
			pointY(iSeg3,numPointEnd(iSeg3)) = segY(iSeg3,2)
			endif
		return		! we have converged

		else		! do another iteration
		icount = icount + 1
		if(icount.gt.100)then
			write(75,*)'Failure in Grow2Nodes. Node = ',iNode
			!pause 'pausing'
			return		!Bail out
			endif
		AreaNold = Area(P1,P2,Nold)
		if(debug.eq.1)then
			write(12,*)'Node X,Y original ',nodeX(iNode),nodeY(iNode)
			write(12,*)'Node X,Y    new   ',Nold(1),Nold(2)
			write(12,*)'AreaOld     ',AreaNold
			write(12,*)'Areanew     ',AreaNnew
			write(12,*)'Targ1 - Area',Target1-AreaNnew
			write(12,*)'iCount ',iCount
			write(12,*)' Type 0 to continue with next iteration; 1 = ABORT'
			read(*,*)igo
			if(igo.eq.1)return
			endif
		go to 10
		endif


	end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine FindPoint(Dist,Slope,Intcpt,P3,Nnew,PI,debug)
	! Routine to calculate the X,Y coordinates of a point along a line given that X,Y are a specified distance from another point (P3)
	implicit none
	Real*8 Slope	! slope of line
	real*8 Intcpt	! intercept of line
	Real*8 P3(2)	! X,Y of "another point"
	Real*8 Nnew(2)	! X,Y of point on line
	real*8 PI(2)	! A point on the line P1-P2 that will be used to check which solution (+ or -) is correct
	real*8 Dist	! Distance along line
	real*8 a,b,c	! Coefs in quadratic formula
	real*8 Xnewplus,Ynewplus,Xnewminus,Ynewminus,Dplus,Dminus
	integer*4 debug
	
	a = 1 + slope**2
	b = -2.0d0*P3(1) + 2.0d0*slope*Intcpt - 2.0d0*slope*P3(2)
	c = P3(1)**2 + Intcpt**2 - 2.0d0*Intcpt*P3(2) + P3(2)**2 - Dist		! note that the distance is already squared with abit added
	! There are 2 roots to the equation -- one on either side of the point P3

	Xnewplus = (-b + sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
	Ynewplus = slope*Xnewplus + Intcpt
	Xnewminus = (-b - sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
	Ynewminus = slope*Xnewminus + Intcpt

	!calculate distance from intersection point to the new point
	Dplus  = (Xnewplus  - PI(1))**2 + (YnewPlus  - PI(2))**2		
	Dminus = (Xnewminus - PI(1))**2 + (YnewMinus - PI(2))**2		

	if(debug.eq.1)then
		write(12,*)'---------------------------'
		write(12,*)'     In Sub FindPoint'
		write(12,*)'Dist   Dist^2',Dist
		write(12,*)'P3           ',P3(1),P3(2)
		write(12,*)'PI           ',PI(1),PI(2)
		write(12,*)'Slope  Intcpt',slope,Intcpt
		write(12,*)' a,b,c       ',a,b,c
		write(12,*)'Plus X,Y,D**2',Xnewplus,Ynewplus,Dplus
		write(12,*)'Minu X,Y,D**2',Xnewminus,Ynewminus,Dminus
		write(12,*)
		pause 'pause'
		endif
	! The correct root is the one that is closest to PI
	if(Dplus.lt.Dminus) then
		Nnew(1) = XnewPlus
		Nnew(2) = YnewPlus
		else
		Nnew(1) = XnewMinus
		Nnew(2) = YnewMinus
		endif		
	return
	end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine ListSegInfo()
	! Routine to list every segment and associated information
	implicit none
	include "Diffuse_GB_Hex.inc"
	integer*4 iSeg,segXl1,segXl2		
	
	write(12,*)' '
	write(12,*)' '
	write(12,*)'&&&&&&&& Listing Seg Information &&&&&&&&&&&&'
	write(12,*)' '
	write(12,*)'   iSeg     N1      Pst     Pen    N2      Xl1     MIF1    SgMIF1    Xl2     MIF2   SgMIF2 '
!                     97     147       1       4     148      93       2      94       5       2       2	
  	do iSeg = 1,numSegs
		segXl1 = segCrystalID(iSeg,1)
		segXl2 = segCrystalID(iSeg,2)
! 		if(CrystalMIFID(segXl1).eq.CrystalMIFID(segXl2))cycle		! skip where both phases are the same (i.e. 2 quqrtz grains
		if(numSegReactionPhases(iSeg).eq.0)cycle			! skip when there is no seg reaction
		write(12,65)iSeg,segNodes(iSeg,1),numPointStart(iSeg),numPointEnd(iSeg),segNodes(iSeg,2),		&
				segXl1,CrystalMIFID(segXl1),segMIFID(iSeg,1),						&
				segXl2,CrystalMIFID(segXl2),segMIFID(iSeg,2)				
				
		end do
		
65	format(15I8)
	return
	end
			
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	Subroutine GrowSegs(GridPlot,iSeg,moleScaleFactor,debug)
	! calling routine for Sub MoveNode to get the relevant information for a specific node
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	
	include "Diffuse_GB_Hex.inc"
	include "GibbsFiles.inc"
	include "Assemb.inc"

	integer*4 iSeg,iPoint,MIFID1,MIFID2,segXl1,segXl2,debug
	integer*4 lastPoint
	real*8 P1X,P1Y,P2X,P2Y,P4X,P4Y,				&
		MPh1,MPh2,moleScaleFactor,a,b,c,area,AdG,AdG2,slope,PS,PS2,b2_4ac,dGridLength
	real*8 CX1,CY1,CX2,CY2,pointXnew(maxPoints),pointYnew(maxPoints)
	real*8 deltaX,deltaY,P3Xplus,P3Yplus,P3Xminus,P3Yminus,Dist1plus,Dist1minus,PI,P1_P2,P1_P4

	segXl1 = segCrystalID(iSeg,1)
	segXl2 = segCrystalID(iSeg,2)
	CX1 = CrystalCenterX(segXl1)
	CY1 = CrystalCenterY(segXl1)
	CX2 = CrystalCenterX(segXl2)
	CY2 = CrystalCenterY(segXl2)
	if(debug.eq.1)then
		write(12,*)' '
		write(12,*)' '
		write(12,*)' '
		write(12,*)' '
		write(12,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
		write(12,*)' Crystal IDs and MIFIDs ',segXl1,CrystalMIFID(segXl1),PhName(CrystalMIFID(segXl1))
		write(12,*)' Crystal IDs and MIFIDs ',segXl2,CrystalMIFID(segXl2),PhName(CrystalMIFID(segXl2))
		write(12,*)' Crystal IDs and Centers - 1:',segXl1,CrystalCenterX(segXl1),CrystalCenterY(segXl1)
		write(12,*)' Crystal IDs and Centers - 2:',segXl2,CrystalCenterX(segXl2),CrystalCenterY(segXl2)
		write(12,*)' '
		write(12,*)' SegNodes: ',segNodes(iSeg,1),numPointStart(iSeg),numPointEnd(iSeg),segNodes(iSeg,2)
		write(12,*)' '
		write(12,*)' '
		endif

	lastPoint = numPointEnd(iSeg) - 1
	do iPoint = numPointStart(iSeg)+1,numPointEnd(iSeg)-1		! we don't move the endpoints because they are coincident with the node

		! Get the MIF IDs so we know what the 2 phases are. 
		! These are needed so we can scale the moles to the number of oxygens
		MIFID1 = segMIFID(iSeg,1)
		MPh1 = pointPhaseMolesDelta(iSeg,iPoint,1)
		MIFID2 = segMIFID(iSeg,2)
		MPh2 = pointPhaseMolesDelta(iSeg,iPoint,2)
		! This just scales to moles -- the scaling factor is set when the model is initiated
		MPh1 = MPh1*moleScaleFactor
		MPh2 = MPh2*moleScaleFactor
		if(debug.eq.1)then
			write(12,*)' -----------------------'
			write(12,*)' -----------------------'
			write(12,*)' Seg, point X,Y ',iSeg,iPoint,pointX(iSeg,iPoint),pointY(iSeg,iPoint)
			write(12,*)' MIFID   Phase   deltaMoles     deltaMoles*Scale'
			write(12,103)MIFID1,PhName(MIFID1),pointPhaseMolesDelta(iSeg,iPoint,1),MPh1
			write(12,103)MIFID2,PhName(MIFID2),pointPhaseMolesDelta(iSeg,iPoint,2),MPh2
103			format(T8,I4,2x,A8,2E15.5,5x,I8,2E5.5)
			write(*,*)'Input value for MPh1 -- note that MPh2 = -MPh1 (oxygen balance)'
			read(*,*)MPh1
			write(12,*)' MPh1 value used:  ',MPh1
			write(12,*)' -----------------------'
			endif
		if(MPH1.eq.0)then
			pointXnew(iPoint) = pointX(iSeg,iPoint)
			pointYnew(iPoint) = pointY(iSeg,iPoint)
			cycle		! skip this point if the ∆moles = 0
					! note that MPH1 = -MPH2 (mass conservation)
			endif

		! use either the forward (or backward) point to get the slope 
		! Don't use the node points (which are actually part of the segment) because these have already been moved in the move node routine
		P1X = pointX(iSeg,iPoint)
		P1Y = pointY(iSeg,iPoint)
		if(iPoint.eq.lastPoint)then
			P2X = pointX(iSeg,iPoint-1)		! use the backward point to get the slope
			P2Y = pointY(iSeg,iPoint-1)
			P4X = pointX(iSeg,iPoint+1)
			P4Y = pointY(iSeg,iPoint+1)
			else
			P2X = pointX(iSeg,iPoint+1)		! use the frontward point to get the slope
			P2Y = pointY(iSeg,iPoint+1)
			P4X = pointX(iSeg,iPoint-1)
			P4Y = pointY(iSeg,iPoint-1)
			endif

!	Geometry we are solving
! 
!                        P3plus
!                .................
!                .       |       .
!                .       |       .
!                .       |       .
!        .---------------.---------------.
!       P4              P1           P2
! 
! 
! 	The area inside the dotted box needs to match the changes in moles
! 	We determine how far to move point P1 along the line perpendicular to line P4-P2 to get the correct area
! 	Then we move the segment point of interest (P1) the same distance along the same line to point P3
! 	Note that segment point (P1) is not necessarily the same position as point Midpoint between P4 and P2 
! 	We don't know whether the MidP needs to be above the line or below the line so we calculate P3plus and P3minus


! 			AdG = area/dXGrid
		! dGrid = distance between P4 and P2
		! dGridLength = (sqrt((P4Y-P2Y)**2 + (P4X-P2X)**2))/2.0d0
		! Grid length is the reaction area between points P1-P2 and P1-P4 (code added July 5, 2024)
! 		P1_P2 = (sqrt((P1Y-P2Y)**2 + (P1X-P2X)**2))
! 		P1_P4 = (sqrt((P1Y-P4Y)**2 + (P1X-P4X)**2))
! 		dGridLength = (P1_P2 + P1_P4)/2.0d0

		dGridLength = (sqrt((P1Y-P2Y)**2 + (P1X-P2X)**2))

!		if(abs(P4Y-P2Y).lt.1.0d-3)then
		if(abs(P1Y-P2Y).lt.1.0d-3)then
			! line is horizontal -- special case
			if(debug.eq.1)	write(12,*)' Horizontal slope case'
			pointXnew(iPoint) = P1X
			P3Xplus  = P1X
			P3Xminus = P1X
			area = abs(MPH1)			! the change in area is scaled to the change in moles of phase 1
			P3Yplus  = P1Y + Area/abs((P2X-P1X))
			P3Yminus = P1Y - Area/abs((P2X-P1X))
! 			P3Yplus  = P1Y + Area/dGridLength
! 			P3Yminus = P1Y - Area/dGridLength

			
			else
			if(debug.eq.1)	write(12,*)' Quadratic formula case'
			! line is not horizontal -- calculate new Y value along the slope			
!			slope = (P4Y-P2Y)/(P4X-P2X)
			slope = (P1Y-P2Y)/(P1X-P2X)
			if (abs(slope).lt.1d-20)slope = 1.D-20		! just to ensure we don't have a zero divide
			PS = -1.0d0/slope				! perpendicular slope
			PI = P1Y - PS*P1X				! Intercept of perpendicular line
			area = abs(MPH1)			! the change in area is scaled to the change in moles of phase 1
			AdG = area/dGridLength
			AdG2 = AdG**2
			! quadratic formula solution
			! 0 = ax^2 + bX + c
			PS2 = 1.0d0 + PS**2	! = a
			a = PS2
			b = 2.0d0*(PS*PI - P1X - P1Y*PS) 
!			b = -2.0d0*P1X*PS2
			c = P1X**2 + P1Y**2 + PI**2 -2.0d0*P1Y*PI- AdG2
!			c = -AdG2 + PS2*P1X**2 
			b2_4ac = b**2 - 4.0d0*a*c
			if(b2_4ac.le.0.0d0)then
				! if the value is lt zero then the new point essentially falls on the old point (this should never happen!)
				P3Xplus  = P1X
				P3Xminus = P1X
				P3Yplus  = P1Y
				P3Yminus = P1Y
				else
! 				P3Xplus =  (-b + sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
! 				P3Xminus = (-b - sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
				P3Xplus =  (-b + sqrt(b2_4ac))/(2.0d0*a)
				P3Xminus = (-b - sqrt(b2_4ac))/(2.0d0*a)
				P3Yplus  = PI + P3Xplus*PS
				P3Yminus = PI + P3Xminus*PS
! 				P3Yplus  = P1Y - PS*(P1X - P3Xplus)
! 				P3Yminus = P1Y - PS*(P1X - P3Xminus)
				endif
			endif
					

		! Calculate the offset of P3 from the midpoint
		deltaX  = abs(P3XPlus - P1X)
		deltaY  = abs(P3YPlus - P1Y)


		! Now we have to figure out which way to move the original point
		! We want to move it according to the value of ∆MPh1 < 0 or ∆MPh1 > 0			
		! Distance from Crystal center 1 to P3
		! The distance is used to compare with the moles produced to figure out which direction Y must go		
		Dist1plus  = (CX1-P3Xplus)**2  + (CY1-P3Yplus)**2		! no need to take the square root
		Dist1minus = (CX1-P3Xminus)**2 + (CY1-P3Yminus)**2
! 		if(MPH1.gt.0.0d0)then  ! this is the ∆moles of the first of the 2 segment phases	
! 			! we need to move away from the center of crystal 1
! 			if(Dist1Plus.gt.Dist1minus)then
! 				pointXnew(iPoint) = pointX(iSeg,iPoint) + deltaX
! 				pointYnew(iPoint) = pointY(iSeg,iPoint) + deltaY
! 				else
! 				pointXnew(iPoint) = pointX(iSeg,iPoint) - deltaX
! 				pointYnew(iPoint) = pointY(iSeg,iPoint) - deltaY
! 				endif
! 			else	!∆moles is negative = move towards crystal 1
! 			if(Dist1Plus.gt.Dist1minus)then
! 				pointXnew(iPoint) = pointX(iSeg,iPoint) - deltaX
! 				pointYnew(iPoint) = pointY(iSeg,iPoint) - deltaY
! 				else
! 				pointXnew(iPoint) = pointX(iSeg,iPoint) + deltaX
! 				pointYnew(iPoint) = pointY(iSeg,iPoint) + deltaY
! 				endif
! 			endif
		if(MPH1.gt.0.0d0)then  ! this is the ∆moles of the first of the 2 segment phases	
			! we need to move away from the center of crystal 1
			if(Dist1Plus.gt.Dist1minus)then
				pointXnew(iPoint) = P3XPlus
				pointYnew(iPoint) = P3YPlus
				else
				pointXnew(iPoint) = P3XMinus
				pointYnew(iPoint) = P3YMinus
				endif
			else	!∆moles is negative = move towards crystal 1
			if(Dist1Plus.gt.Dist1minus)then
				pointXnew(iPoint) = P3XMinus
				pointYnew(iPoint) = P3YMinus
				else
				pointXnew(iPoint) = P3XPlus
				pointYnew(iPoint) = P3YPlus
				endif
			endif


		if(debug.eq.1)then
			write(12,*)'P4             ',P4X,P4Y		
			write(12,*)'P1             ',P1X,P1Y		
			write(12,*)'P2             ',P2X,P2Y		
			if(abs(P1Y-P2Y).lt.1.0d-3)then
				write(12,*)' Horizontal slope case'
				write(12,*)'GridLength     ',dGridLength		
				write(12,*)'Area           ',Area		
				else
				write(12,*)' Quadratic formula case'
				write(12,*)'Slope          ',slope 
				write(12,*)'PerpSlope      ',PS 
				write(12,*)'PerpIntcp      ',PI 
				write(12,*)'GridLength     ',dGridLength		
				write(12,*)'Area           ',Area		
				write(12,*)'AdG            ',AdG		
				write(12,*)'b2_4ac         ',b2_4ac		
				endif
			write(12,*)'---plus-minus results---'
			write(12,*)'P3Plus         ',P3Xplus,P3Yplus		
			write(12,*)'P1             ',P1X,P1Y		
			write(12,*)'P3Minus        ',P3XMinus,P3YMinus		
			write(12,*)'delta_XY       ',deltaX,deltaY		
			write(12,*)'MPh1           ',MPh1
			write(12,*)'Dist+  Dist-   ',Dist1Plus,Dist1Minus
			write(12,*)'Old points     ',pointX(iSeg,iPoint),pointY(iSeg,iPoint)
			write(12,*)'New points     ',pointXnew(iPoint),pointYnew(iPoint)
			pause ' Pausing... do the next point'
			endif

		end do		! loop to the next point in the segment

	! set the new point positions into the proper array
	do iPoint = numPointStart(iSeg)+1,numPointEnd(iSeg)-1	! we don't move the endpoints because they are coincident with the node
		pointX(iSeg,iPoint) = pointXnew(iPoint)
		pointY(iSeg,iPoint) = pointYnew(iPoint)
		end do

	return
	end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	Subroutine GrowSegs_old_v4(GridPlot,iSeg,moleScaleFactor,debug)
	! calling routine for Sub MoveNode to get the relevant information for a specific node
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	
	include "Diffuse_GB_Hex.inc"
	include "GibbsFiles.inc"
	include "Assemb.inc"

	integer*4 iSeg,iPoint,MIFID1,MIFID2,segXl1,segXl2,debug
	integer*4 lastPoint
	real*8 P1X,P1Y,P2X,P2Y,P3X,P3Y,				&
		MPh1,MPh2,moleScaleFactor,a,b,c,area,AdG,AdG2,slope,PS,PS2,b2_4ac
	real*8 CX1,CY1,CX2,CY2,pointXnew(maxPoints),pointYnew(maxPoints)
	real*8 MidPX,MidPY,P1P2Distance,deltaX,deltaY,MidPxP3,P3Xplus,P3Yplus,P3Xminus,P3Yminus,Dist1plus,Dist1minus

	segXl1 = segCrystalID(iSeg,1)
	segXl2 = segCrystalID(iSeg,2)
	CX1 = CrystalCenterX(segXl1)
	CY1 = CrystalCenterY(segXl1)
	CX2 = CrystalCenterX(segXl2)
	CY2 = CrystalCenterY(segXl2)
	if(debug.eq.1)then
		write(12,*)' '
		write(12,*)' '
		write(12,*)' '
		write(12,*)' '
		write(12,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
		write(12,*)' Crystal IDs and MIFIDs ',segXl1,CrystalMIFID(segXl1),PhName(CrystalMIFID(segXl1))
		write(12,*)' Crystal IDs and MIFIDs ',segXl2,CrystalMIFID(segXl2),PhName(CrystalMIFID(segXl2))
		write(12,*)' Crystal IDs and Centers - 1:',segXl1,CrystalCenterX(segXl1),CrystalCenterY(segXl1)
		write(12,*)' Crystal IDs and Centers - 2:',segXl2,CrystalCenterX(segXl2),CrystalCenterY(segXl2)
		write(12,*)' '
		write(12,*)' SegNodes: ',segNodes(iSeg,1),numPointStart(iSeg),numPointEnd(iSeg),segNodes(iSeg,2)
		write(12,*)' '
		write(12,*)' '
		endif

	lastPoint = numPointEnd(iSeg) - 1
	do iPoint = numPointStart(iSeg)+1,numPointEnd(iSeg)-1		! we don't move the endpoints because they are coincident with the node

		! Get the MIF IDs so we know what the 2 phases are. 
		! These are needed so we can scale the moles to the number of oxygens
		MIFID1 = segMIFID(iSeg,1)
		MPh1 = pointPhaseMolesDelta(iSeg,iPoint,1)
		MIFID2 = segMIFID(iSeg,2)
		MPh2 = pointPhaseMolesDelta(iSeg,iPoint,2)
		! This just scales to moles -- the scaling factor is set when the model is initiated
		MPh1 = MPh1*moleScaleFactor
		MPh2 = MPh2*moleScaleFactor
		if(debug.eq.1)then
			write(12,*)' -----------------------'
			write(12,*)' -----------------------'
			write(12,*)' Seg, point X,Y ',iSeg,iPoint,pointX(iSeg,iPoint),pointY(iSeg,iPoint)
			write(12,*)' MIFID   Phase   deltaMoles     deltaMoles*Scale'
			write(12,103)MIFID1,PhName(MIFID1),pointPhaseMolesDelta(iSeg,iPoint,1),MPh1
			write(12,103)MIFID2,PhName(MIFID2),pointPhaseMolesDelta(iSeg,iPoint,2),MPh2
103			format(T8,I4,2x,A8,2E15.5,5x,I8,2E5.5)
			write(*,*)'Input value for MPh1 -- note that MPh2 = -MPh1 (oxygen balance)'
			read(*,*)MPh1
			write(12,*)' MPh1 value used:  ',MPh1
			write(12,*)' -----------------------'
			endif
		if(MPH1.eq.0)then
			pointXnew(iPoint) = pointX(iSeg,iPoint)
			pointYnew(iPoint) = pointY(iSeg,iPoint)
			cycle		! skip this point if the ∆moles = 0
					! note that MPH1 = -MPH2 (mass conservation)
			endif

		! use points on either side of the point in question to get the slope
		P1X = pointX(iSeg,iPoint+1)
		P1Y = pointY(iSeg,iPoint+1)
		P2X = pointX(iSeg,iPoint-1)	
		P2Y = pointY(iSeg,iPoint-1)

!	Geometry we are solving
! 
!                        P3plus
!                .................
!                .       |       .
!                .       |       .
!                .       |       .
!        .---------------.---------------.
!       P2              MidP           P`
! 
! 
! 	The area inside the dotted box needs to match the changes in moles
! 	We determine how far to move point MidP along the line perpendicular to line P1-P2 to get the correct area
! 	Then we move the segment point of interest (P) the same distance along the same line to point P3
! 	Note that segment point (P) is not necessarily the same position as point MidP 
! 	We don't know whether the MidP needs to be above the line or below the line so we calculate P3plus and P3minus

		slope = (P1Y-P2Y)/(P1X-P2X)	! slope of line P1-P2
! 		area = MPH1			! the change in area is scaled to the change in moles of phase 1
		area = abs(MPH1)			! the change in area is scaled to the change in moles of phase 1
		MidPX = (P1X + P2X)/2.0d0
		MidPY = (P1Y + P2Y)/2.0d0
		if(debug.eq.1)then
			write(12,*)'P1_XY       ',P1X,P1Y		
			write(12,*)'midPoint    ',midPX,midPY
			write(12,*)'P2_XY       ',P2X,P2Y		
			endif
		if(abs(P1Y-P2Y).lt.1.0d-3)then
			! line is nearly horizontal -- special case
			pointXnew(iPoint) = midPX
			pointYnew(iPoint) = midPY
			P3X  = MidPX
			P3Y  = MidPY + Area/((P2X-P1X)/2.0d0)
			
			else
			! line is not horizontal -- calculate new Y value along the slope			
			if (abs(slope).lt.1d-20)slope = 1.D-20		! just to ensure we don't have a zero divide
			PS = -1.0d0/slope				! slope of the perpendicular line to points P1-P2
			! distance between P1 and P2
			P1P2distance = sqrt((P2X-P1X)**2 + (P2Y-P1Y)**2) 
			AdG = area/(P1P2distance/2.0d0)
			AdG2 = AdG**2
			! quadratic formula solution

			PS2 = 1.0d0 + PS**2
			a = PS2
			b = -2.0d0*MidPX*PS2
			c = -AdG2 + PS2*MidPX**2
			
			b2_4ac = b**2 - 4.0d0*a*c
			if(b2_4ac.le.0.0d0)then
				! if the value is lt zero then the new point essentially falls on the old point (this should never happen!)
				P3X  = MidPX
				P3Y  = MidPY
				else
				P3X =  (-b + sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
				P3Y  = MidPY - PS*(MidPX - P3X)
				endif
			endif
					
		

		! We need to determine whether vector P2-->P3 is clockwise or anticlockwise from the vector P2-->MidP
		! We do this using the cross product
 		MidPxP3 = (midPX - P2X)*(P3Y - P2Y) - (midPY - P2Y)*(P3X - P2X)
		! if MidPxP3 > 0 then P3 is counter clockwise from P2-->MidP. Otherwise it is clockwise

		! Calculate the offset of P3 from the midpoint
		deltaX  = abs(P3X - midPX)
		deltaY  = abs(P3Y - midPY)
		P3Xplus = midPX + deltaX
		P3Yplus = midPY + deltaY
 		P3Xminus = midPX - deltaX
 		P3Yminus = midPY - deltaY


		! Now we have to figure out which way to move the original point
		! We want to move it according to the value of ∆MPh1 < 0 or ∆MPh1 > 0			
		! Distance from Crystal center 1 to P3
		! The distance is used to compare with the moles produced to figure out which direction Y must go		
		Dist1plus  = (CX1-P3Xplus)**2  + (CY1-P3Yplus)**2		! no need to take the square root
		Dist1minus = (CX1-P3Xminus)**2 + (CY1-P3Yminus)**2
		if(MPH1.gt.0.0d0)then  ! this is the ∆moles of the first of the 2 segment phases	
			! we need to move away from the center of crystal 1
			if(Dist1Plus.gt.Dist1minus)then
				pointXnew(iPoint) = pointX(iSeg,iPoint) + deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) + deltaY
				else
				pointXnew(iPoint) = pointX(iSeg,iPoint) - deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) - deltaY
				endif
			else	!∆moles is negative = move towards crystal 1
			if(Dist1Plus.gt.Dist1minus)then
				pointXnew(iPoint) = pointX(iSeg,iPoint) - deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) - deltaY
				else
				pointXnew(iPoint) = pointX(iSeg,iPoint) + deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) + deltaY
				endif
			endif


		if(debug.eq.1)then
			write(12,*)'P2             ',P2X,P2Y		
			write(12,*)'midPoint       ',midPX,midPY
			write(12,*)'P3             ',P3X,P3Y		
			write(12,*)'delta_XY       ',deltaX,deltaY		
			write(12,*)'MPh1           ',MPh1
			write(12,*)'Dist+  Dist-   ',Dist1Plus,Dist1Minus
			write(12,*)'Old points     ',pointX(iSeg,iPoint),pointY(iSeg,iPoint)
			write(12,*)'New points     ',pointXnew(iPoint),pointYnew(iPoint)
			pause ' Pausing... do the next point'
			endif

		end do		! loop to the next point in the segment

	! set the new point positions into the proper array
	do iPoint = numPointStart(iSeg)+1,numPointEnd(iSeg)-1	! we don't move the endpoints because they are coincident with the node
		pointX(iSeg,iPoint) = pointXnew(iPoint)
		pointY(iSeg,iPoint) = pointYnew(iPoint)
		end do

	return
	end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	Subroutine GrowSegs_old_v3(GridPlot,iSeg,moleScaleFactor,debug)
	! calling routine for Sub MoveNode to get the relevant information for a specific node
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	
	include "Diffuse_GB_Hex.inc"
	include "GibbsFiles.inc"
	include "Assemb.inc"

	integer*4 iSeg,iPoint,MIFID1,MIFID2,segXl1,segXl2,debug
	integer*4 lastPoint
	real*8 P1X,P1Y,P2X,P2Y,P3X,P3Y,				&
		MPh1,MPh2,moleScaleFactor,a,b,c,area,AdG,AdG2,slope,PS,PS2,b2_4ac
	real*8 CX1,CY1,CX2,CY2,pointXnew(maxPoints),pointYnew(maxPoints)
	real*8 MidPX,MidPY,P1P2Distance,deltaX,deltaY,MidPxP3,P3Xplus,P3Yplus

	segXl1 = segCrystalID(iSeg,1)
	segXl2 = segCrystalID(iSeg,2)
	CX1 = CrystalCenterX(segXl1)
	CY1 = CrystalCenterY(segXl1)
	CX2 = CrystalCenterX(segXl2)
	CY2 = CrystalCenterY(segXl2)
	if(debug.eq.1)then
		write(12,*)' '
		write(12,*)' '
		write(12,*)' '
		write(12,*)' '
		write(12,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
		write(12,*)' Crystal IDs and MIFIDs ',segXl1,CrystalMIFID(segXl1),PhName(CrystalMIFID(segXl1))
		write(12,*)' Crystal IDs and MIFIDs ',segXl2,CrystalMIFID(segXl2),PhName(CrystalMIFID(segXl2))
		write(12,*)' Crystal IDs and Centers - 1:',segXl1,CrystalCenterX(segXl1),CrystalCenterY(segXl1)
		write(12,*)' Crystal IDs and Centers - 2:',segXl2,CrystalCenterX(segXl2),CrystalCenterY(segXl2)
		write(12,*)' '
		write(12,*)' SegNodes: ',segNodes(iSeg,1),numPointStart(iSeg),numPointEnd(iSeg),segNodes(iSeg,2)
		write(12,*)' '
		write(12,*)' '
		endif

	lastPoint = numPointEnd(iSeg) - 1
	do iPoint = numPointStart(iSeg)+1,numPointEnd(iSeg)-1		! we don't move the endpoints because they are coincident with the node

		! Get the MIF IDs so we know what the 2 phases are. 
		! These are needed so we can scale the moles to the number of oxygens
		MIFID1 = segMIFID(iSeg,1)
		MPh1 = pointPhaseMolesDelta(iSeg,iPoint,1)
		MIFID2 = segMIFID(iSeg,2)
		MPh2 = pointPhaseMolesDelta(iSeg,iPoint,2)
		! This just scales to moles -- the scaling factor is set when the model is initiated
		MPh1 = MPh1*moleScaleFactor
		MPh2 = MPh2*moleScaleFactor
		if(debug.eq.1)then
			write(12,*)' -----------------------'
			write(12,*)' -----------------------'
			write(12,*)' Seg, point X,Y ',iSeg,iPoint,pointX(iSeg,iPoint),pointY(iSeg,iPoint)
			write(12,*)' MIFID   Phase   deltaMoles     deltaMoles*Scale'
			write(12,103)MIFID1,PhName(MIFID1),pointPhaseMolesDelta(iSeg,iPoint,1),MPh1
			write(12,103)MIFID2,PhName(MIFID2),pointPhaseMolesDelta(iSeg,iPoint,2),MPh2
103			format(T8,I4,2x,A8,2E15.5,5x,I8,2E5.5)
			write(*,*)'Input value for MPh1 -- note that MPh2 = -MPh1 (oxygen balance)'
			read(*,*)MPh1
			write(12,*)' MPh1 value used:  ',MPh1
			write(12,*)' -----------------------'
			endif
		if(MPH1.eq.0)then
			pointXnew(iPoint) = pointX(iSeg,iPoint)
			pointYnew(iPoint) = pointY(iSeg,iPoint)
			cycle		! skip this point if the ∆moles = 0
					! note that MPH1 = -MPH2 (mass conservation)
			endif

		! use points on either side of the point in question to get the slope
		P1X = pointX(iSeg,iPoint+1)
		P1Y = pointY(iSeg,iPoint+1)
		P2X = pointX(iSeg,iPoint-1)	
		P2Y = pointY(iSeg,iPoint-1)

!	Geometry we are solving
! 
!                        P3plus
!                .................
!                .       |       .
!                .       |       .
!                .       |       .
!        .---------------.---------------.
!       P2              MidP           P`
! 
! 
! 	The area inside the dotted box needs to match the changes in moles
! 	We determine how far to move point MidP along the line perpendicular to line P1-P2 to get the correct area
! 	Then we move the segment point of interest (P) the same distance along the same line to point P3
! 	Note that segment point (P) is not necessarily the same position as point MidP 
! 	We don't know whether the MidP needs to be above the line or below the line so we calculate P3plus and P3minus

		slope = (P1Y-P2Y)/(P1X-P2X)	! slope of line P1-P2
! 		area = MPH1			! the change in area is scaled to the change in moles of phase 1
		area = abs(MPH1)			! the change in area is scaled to the change in moles of phase 1
		MidPX = (P1X + P2X)/2.0d0
		MidPY = (P1Y + P2Y)/2.0d0
		if(debug.eq.1)then
			write(12,*)'P1_XY       ',P1X,P1Y		
			write(12,*)'midPoint    ',midPX,midPY
			write(12,*)'P2_XY       ',P2X,P2Y		
			endif
		if(abs(P1Y-P2Y).lt.1.0d-3)then
			! line is nearly horizontal -- special case
			pointXnew(iPoint) = midPX
			pointYnew(iPoint) = midPY
			P3X  = MidPX
			P3Y  = MidPY + Area/((P2X-P1X)/2.0d0)
			
			else
			! line is not horizontal -- calculate new Y value along the slope			
			if (abs(slope).lt.1d-20)slope = 1.D-20		! just to ensure we don't have a zero divide
			PS = -1.0d0/slope				! slope of the perpendicular line to points P1-P2
			! distance between P1 and P2
			P1P2distance = sqrt((P2X-P1X)**2 + (P2Y-P1Y)**2) 
			AdG = area/(P1P2distance/2.0d0)
			AdG2 = AdG**2
			! quadratic formula solution

			PS2 = 1.0d0 + PS**2
			a = PS2
			b = -2.0d0*MidPX*PS2
			c = -AdG2 + PS2*MidPX**2
			
			b2_4ac = b**2 - 4.0d0*a*c
			if(b2_4ac.le.0.0d0)then
				! if the value is lt zero then the new point essentially falls on the old point (this should never happen!)
				P3X  = MidPX
				P3Y  = MidPY
				else
				P3X =  (-b + sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
				P3Y  = MidPY - PS*(MidPX - P3X)
				endif
			endif
					
		

		! We need to determine whether vector P2-->P3 is clockwise or anticlockwise from the vector P2-->MidP
		! We do this using the cross product
 		MidPxP3 = (midPX - P2X)*(P3Y - P2Y) - (midPY - P2Y)*(P3X - P2X)
		! if MidPxP3 > 0 then P3 is counter clockwise from P2-->MidP. Otherwise it is clockwise

		! Calculate the offset of P3 from the midpoint
		deltaX  = abs(P3X - midPX)
		deltaY  = abs(P3Y - midPY)
		P3Xplus = midPX + deltaX
		P3Yplus = midPY + deltaY
! 		P3Xminus = midPX - deltaX
! 		P3Yminus = midPY - deltaY

		! We need to determine whether vector P2-->P3 is clockwise or anticlockwise from the vector P2-->MidP
! 		Calculate the cross product of the vector 
 		MidPxP3 = (midPX - P2X)*(P3Yplus - P2Y) - (midPY - P2Y)*(P3Xplus - P2X)

		! if MidPxP3 > 0 then P3plus is counterclockwise from P2-->MidP. Otherwise it is clockwise

		! Now we have to figure out which way to move the original point
		! We want to move it according to the value of ∆MPh1 < 0 or ∆MPh1 > 0			
		if(MPH1.gt.0.0d0)then  ! this is the ∆moles of the first of the 2 segment phases	
			! We need to increase the area of the first phase (clockwise
			! We do this by moving the point in a counter-clockwise direction
			! ∆moles is positive
			if(MidPxP3.lt.0)then 	! P3plus is clockwise and 
				pointXnew(iPoint) = pointX(iSeg,iPoint) - deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) - deltaY
				else
				pointXnew(iPoint) = pointX(iSeg,iPoint) + deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) + deltaY
				endif
							
			else  ! the ∆moles is negative
			if(MidPxP3.lt.0)then
				! We need to decrease the area of the first phase
				! We do this by moving the point in a clockwise direction
				pointXnew(iPoint) = pointX(iSeg,iPoint) + deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) + deltaY
				else
				pointXnew(iPoint) = pointX(iSeg,iPoint) - deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) - deltaY
				endif

			endif


		if(debug.eq.1)then
			write(12,*)'P2             ',P2X,P2Y		
			write(12,*)'midPoint       ',midPX,midPY
			write(12,*)'P3             ',P3X,P3Y		
! 			write(12,*)'Cross product  ',MidPxP3
			write(12,*)'delta_XY       ',deltaX,deltaY		
			write(12,*)'Old points     ',pointX(iSeg,iPoint),pointY(iSeg,iPoint)
			write(12,*)'New points     ',pointXnew(iPoint),pointYnew(iPoint)
			pause ' Pausing... do the next point'
			endif

		end do		! loop to the next point in the segment

	! set the new point positions into the proper array
	do iPoint = numPointStart(iSeg)+1,numPointEnd(iSeg)-1	! we don't move the endpoints because they are coincident with the node
		pointX(iSeg,iPoint) = pointXnew(iPoint)
		pointY(iSeg,iPoint) = pointYnew(iPoint)
		end do

	return
	end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	Subroutine GrowSegs_old_v2(GridPlot,iSeg,moleScaleFactor,debug)
	! calling routine for Sub MoveNode to get the relevant information for a specific node
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	
	include "Diffuse_GB_Hex.inc"
	include "GibbsFiles.inc"
	include "Assemb.inc"

	integer*4 iSeg,iPoint,MIFID1,MIFID2,segXl1,segXl2,debug
	integer*4 lastPoint
	real*8 P1X,P1Y,P2X,P2Y,P3X,P3Y,				&
		MPh1,MPh2,moleScaleFactor,a,b,c,area,AdG,AdG2,slope,PS,PS2,b2_4ac
	real*8 CX1,CY1,CX2,CY2,dist_to_P3,dist_to_midP,pointXnew(maxPoints),pointYnew(maxPoints)
	real*8 MidPX,MidPY,P1P2Distance,deltaX,deltaY

	segXl1 = segCrystalID(iSeg,1)
	segXl2 = segCrystalID(iSeg,2)
	CX1 = CrystalCenterX(segXl1)
	CY1 = CrystalCenterY(segXl1)
	CX2 = CrystalCenterX(segXl2)
	CY2 = CrystalCenterY(segXl2)
	if(debug.eq.1)then
		write(12,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
		write(12,*)' Crystal IDs and MIFIDs ',segXl1,CrystalMIFID(segXl1),PhName(CrystalMIFID(segXl1))
		write(12,*)' Crystal IDs and MIFIDs ',segXl2,CrystalMIFID(segXl2),PhName(CrystalMIFID(segXl2))
		write(12,*)' Crystal IDs and Centers - 1:',segXl1,CrystalCenterX(segXl1),CrystalCenterY(segXl1)
		write(12,*)' Crystal IDs and Centers - 2:',segXl2,CrystalCenterX(segXl2),CrystalCenterY(segXl2)
		write(12,*)' '
		write(12,*)' SegNodes: ',segNodes(iSeg,1),numPointStart(iSeg),numPointEnd(iSeg),segNodes(iSeg,2)
		write(12,*)' '
		write(12,*)' '
		endif

	lastPoint = numPointEnd(iSeg) - 1
	do iPoint = numPointStart(iSeg)+1,numPointEnd(iSeg)-1		! we don't move the endpoints because they are coincident with the node

		! Get the MIF IDs so we know what the 2 phases are. 
		! These are needed so we can scale the moles to the number of oxygens
		MIFID1 = segMIFID(iSeg,1)
		MPh1 = pointPhaseMolesDelta(iSeg,iPoint,1)
		MIFID2 = segMIFID(iSeg,2)
		MPh2 = pointPhaseMolesDelta(iSeg,iPoint,2)
		! This just scales to moles -- the scaling factor is set when the model is initiated
		MPh1 = MPh1*moleScaleFactor
		MPh2 = MPh2*moleScaleFactor
		if(debug.eq.1)then
			write(12,*)' -----------------------'
			write(12,*)' -----------------------'
			write(12,*)' Seg, point X,Y ',iSeg,iPoint,pointX(iSeg,iPoint),pointY(iSeg,iPoint)
			write(12,*)' MIFID   Phase   deltaMoles     deltaMoles*Scale'
			write(12,103)MIFID1,PhName(MIFID1),pointPhaseMolesDelta(iSeg,iPoint,1),MPh1
			write(12,103)MIFID2,PhName(MIFID2),pointPhaseMolesDelta(iSeg,iPoint,2),MPh2
103			format(T8,I4,2x,A8,2E15.5,5x,I8,2E5.5)
			write(*,*)'Input value for MPh1 -- note that MPh2 = -MPh1 (oxygen balance)'
			read(*,*)MPh1
			write(12,*)' MPh1 value used:  ',MPh1
			write(12,*)' -----------------------'
			endif
		if(MPH1.eq.0)then
			pointXnew(iPoint) = pointX(iSeg,iPoint)
			pointYnew(iPoint) = pointY(iSeg,iPoint)
			cycle		! skip this point if the ∆moles = 0
					! note that MPH1 = -MPH2 (mass conservation)
			endif

		! use points on either side of the point in question to get the slope
		P1X = pointX(iSeg,iPoint+1)
		P1Y = pointY(iSeg,iPoint+1)
! 		if(iPoint.eq.lastPoint)then
		P2X = pointX(iSeg,iPoint-1)	
		P2Y = pointY(iSeg,iPoint-1)
! 			else
! 			P2X = pointX(iSeg,iPoint+1)		! use the frontward point to get the slope
! 			P2Y = pointY(iSeg,iPoint+1)
! 			endif

!	Geometry we are solving
! 
!                        P3plus
!                .................
!                .       |       .
!                .       |       .
!                .       |       .
!        .---------------.---------------.
!       P1              MidP           P2
! 
! 
! 	The area inside the dotted box needs to match the changes in moles
! 	We determine how far to move point MidP along the line perpendicular to line P1-P2 to get the correct area
! 	Then we move the segment point of interest (P) the same distance along the same line to point P3
! 	Note that segment point (P) is not necessarily the same position as point MidP 
! 	We don't know whether the MidP needs to be above the line or below the line so we calculate P3plus and P3minus

		slope = (P1Y-P2Y)/(P1X-P2X)	! slope of line P1-P2
		area = MPH1			! the change in area is scaled to the change in moles of phase 1
		MidPX = (P1X + P2X)/2.0d0
		MidPY = (P1Y + P2Y)/2.0d0
		if(debug.eq.1)then
			write(12,*)'P1_XY       ',P1X,P1Y		
			write(12,*)'midPoint    ',midPX,midPY
			write(12,*)'P2_XY       ',P2X,P2Y		
			endif
		if(abs(P1Y-P2Y).lt.1.0d-3)then
			! line is nearly horizontal -- special case
			pointXnew(iPoint) = midPX
			pointYnew(iPoint) = midPY
			P3X  = MidPX
! 			P3Xminus = MidPX
			P3Y  = MidPY + Area/((P2X-P1X)/2.0d0)
! 			P3Yminus = MidPY - Area/((P2X-P1X)/2.0D0)
			
			else
			! line is not horizontal -- calculate new Y value along the slope			
			!slope = (P1Y-P2Y)/(P1X-P2X)
			if (abs(slope).lt.1d-20)slope = 1.D-20		! just to ensure we don't have a zero divide
			PS = -1.0d0/slope				! slope of the perpendicular line to points P1-P2
			! dGrid = distance between P1 and P2
			! distance between P1 and P2
			P1P2distance = sqrt((P2X-P1X)**2 + (P2Y-P1Y)**2) 
!			AdG = area/dXGrid
			AdG = area/(P1P2distance/2.0d0)
			AdG2 = AdG**2
			! quadratic formula solution
! 			PS2 = 1.0d0 + PS**2
! 			c = -AdG2 + PS2*P1X**2
! 			b = -2.0d0*P1X*PS2
! 			a = PS2

			PS2 = 1.0d0 + PS**2
			a = PS2
			b = -2.0d0*MidPX*PS2
			c = -AdG2 + PS2*MidPX**2
			
			b2_4ac = b**2 - 4.0d0*a*c
			if(b2_4ac.le.0.0d0)then
				! if the value is lt zero then the new point essentially falls on the old point (this should never happen!)
! 				P3Xplus  = P1X
! 				P3Xminus = P1X
! 				P3Yplus  = P1Y
! 				P3Yminus = P1Y
				P3X  = MidPX
! 				P3Xminus = MidPX
				P3Y  = MidPY
! 				P3Yminus = MidPY
				else
! 				P3Xplus =  (-b + sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
! 				P3Xminus = (-b - sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
! 				P3Yplus  = P1Y - PS*(P1X - P3Xplus)
! 				P3Yminus = P1Y - PS*(P1X - P3Xminus)
				P3X =  (-b + sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
! 				P3Xminus = (-b - sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
				P3Y  = MidPY - PS*(MidPX - P3X)
! 				P3Yminus = MidPY - PS*(MidPX - P3Xminus)
				endif
			endif
					
		
		! Calculate the offset of P3 from the midpoint
! 		deltaXminus = midPX - P3Xminus
! 		deltaYminus = midPY - P3Yminus
! 		deltaXplus  = P3Xplus - midPX
! 		deltaYplus  = P3Yplus - midPY
		deltaX  = P3X - midPX
		deltaY  = P3Y - midPY
		
		! Distance from Crystal center 1 to P3
		! The distance is used to compare with the moles produced to figure out which direction Y must go		
! 		Dist1plus  = (CX1-P3X)**2  + (CY1-P3Y)**2		! no need to take the square root
! 		Dist1minus = (CX1-P3X)**2 + (CY1-P3Yminus)**2
		Dist_to_P3  = (CX1-P3X)**2  + (CY1-P3Y)**2		! no need to take the square root
		Dist_to_midP = (CX1-midPX)**2 + (CY1-midPY)**2
		if(Dist_to_P3.gt.Dist_to_midP)then
			if(MPH1.gt.0.0d0)then	
				! Dist1Plus is the correct solution
! 				pointXnew(iPoint) = P3Xplus
! 				pointYnew(iPoint) = P3Yplus
				pointXnew(iPoint) = pointX(iSeg,iPoint) + deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) + deltaY
				if(debug.eq.1)then
					write(12,*)' Dist_to_P3 > Dist_to_midP'
					write(12,*)' MPH1 >= 0 '
! 					write(12,*)'midPoint    ',midPX,midPY
! 					write(12,*)'P3plus      ',P3Xplus,P3Yplus		
! 					write(12,*)'P3minus     ',P3Xminus,P3Yminus		
! 					write(12,*)'delta_XY    ',deltaX,deltaY		
! 					write(12,*)'delta_minus ',deltaXminus,deltaYminus		
! 					write(12,*)' Dist_to_P3   ',Dist_to_P3
! 					write(12,*)' Dist_to_midP ',Dist_to_midP
					write(12,*)'Plus is correct'
! 					write(12,*)'Old points  ',pointX(iSeg,iPoint),pointY(iSeg,iPoint)
! 					write(12,*)'New points  ',pointXnew(iPoint),pointYnew(iPoint)
! 					pause ' Pausing... do the next point'
					endif
				else
				! Dist1minus is the correct solution
! 				pointXnew(iPoint) = P3Xminus
! 				pointYnew(iPoint) = P3Yminus
				pointXnew(iPoint) = pointX(iSeg,iPoint) - deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) - deltaY
				if(debug.eq.1)then
					write(12,*)' Dist_to_P3 > Dist_to_midP'
					write(12,*)' MPH1 < 0 '
! 					write(12,*)'midPoint    ',midPX,midPY
! 					write(12,*)'P3plus      ',P3Xplus,P3Yplus		
! 					write(12,*)'P3minus     ',P3Xminus,P3Yminus		
! 					write(12,*)'delta_plus  ',deltaX,deltaY		
! 					write(12,*)'delta_minus ',deltaXminus,deltaYminus		
! 					write(12,*)' Dist_to_P3   ',Dist_to_P3
! 					write(12,*)' Dist_to_midP ',Dist_to_midP
					write(12,*)'Minus is correct'
! 					write(12,*)'Old points  ',pointX(iSeg,iPoint),pointY(iSeg,iPoint)
! 					write(12,*)'New points  ',pointXnew(iPoint),pointYnew(iPoint)
! 					pause ' Pausing... do the next point'
					endif
				endif
			else
			if(MPH1.gt.0.0d0)then	
				! Dist1minus is the correct solution
! 				pointXnew(iPoint) = P3Xminus
! 				pointYnew(iPoint) = P3Yminus
				pointXnew(iPoint) = pointX(iSeg,iPoint) - deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) - deltaY
				if(debug.eq.1)then
					write(12,*)' Dist_to_P3 < Dist_to_midP'
					write(12,*)' MPH1 >= 0 '
! 					write(12,*)'midPoint    ',midPX,midPY
! 					write(12,*)'P3plus      ',P3Xplus,P3Yplus		
! ! 					write(12,*)'P3minus     ',P3Xminus,P3Yminus		
! 					write(12,*)'delta_plus  ',deltaX,deltaY		
! 					write(12,*)'delta_minus ',deltaXminus,deltaYminus		
! 					write(12,*)' Dist_to_P3   ',Dist_to_P3
! 					write(12,*)' Dist_to_midP ',Dist_to_midP
					write(12,*)'Minus is correct'
! 					write(12,*)'Old points  ',pointX(iSeg,iPoint),pointY(iSeg,iPoint)
! 					write(12,*)'New points  ',pointXnew(iPoint),pointYnew(iPoint)
! 					pause ' Pausing... do the next point'
					endif
				else
				! Dist1plus is the correct solution
! 				pointXnew(iPoint) = P3Xplus
! 				pointYnew(iPoint) = P3Yplus
				pointXnew(iPoint) = pointX(iSeg,iPoint) + deltaX
				pointYnew(iPoint) = pointY(iSeg,iPoint) + deltaY
				if(debug.eq.1)then
					write(12,*)' Dist_to_P3 < Dist_to_midP'
					write(12,*)' MPH1 < 0 '
! 					write(12,*)'midPoint    ',midPX,midPY
! 					write(12,*)'P3plus      ',P3Xplus,P3Yplus		
! ! 					write(12,*)'P3minus     ',P3Xminus,P3Yminus		
! 					write(12,*)'delta_plus  ',deltaX,deltaY		
! 					write(12,*)'delta_minus ',deltaXminus,deltaYminus		
! 					write(12,*)' Dist_to_P3   ',Dist_to_P3
! 					write(12,*)' Dist_to_midP ',Dist_to_midP
					write(12,*)'Plus is correct'
! 					write(12,*)'Old points  ',pointX(iSeg,iPoint),pointY(iSeg,iPoint)
! 					write(12,*)'New points  ',pointXnew(iPoint),pointYnew(iPoint)
! 					pause ' Pausing... do the next point'
					endif
				endif
			endif

		if(debug.eq.1)then
			write(12,*)' Dist_to_P3   ',Dist_to_P3
			write(12,*)' Dist_to_midP ',Dist_to_midP
			write(12,*)'midPoint    ',midPX,midPY
			write(12,*)'P3plus      ',P3X,P3Y		
			write(12,*)'delta_XY    ',deltaX,deltaY		
			write(12,*)'Old points  ',pointX(iSeg,iPoint),pointY(iSeg,iPoint)
			write(12,*)'New points  ',pointXnew(iPoint),pointYnew(iPoint)
			pause ' Pausing... do the next point'
			endif

				

		end do		! loop to the next point in the segment

	! set the new point positions into the proper array
	do iPoint = numPointStart(iSeg)+1,numPointEnd(iSeg)-1	! we don't move the endpoints because they are coincident with the node
		pointX(iSeg,iPoint) = pointXnew(iPoint)
		pointY(iSeg,iPoint) = pointYnew(iPoint)
		end do

! 	write(*,*)'************************'
! 	write(*,*)' Plot the segment? 0 = no, 1 = yes'
! 	read(*,*)iPlot
! 	if(iPlot.eq.1)then
! 		! plot the segment
! 		call PlotGrid(GridPlot,1)
! 		endif
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	Subroutine GrowSegs_old_v1(iSeg,moleScaleFactor,debug)
	! calling routine for Sub MoveNode to get the relevant information for a specific node
	
	implicit none	
	include "Diffuse_GB_Hex.inc"
	include "GibbsFiles.inc"
	include "Assemb.inc"

	integer*4 iSeg,iPoint,MIFID1,MIFID2,segXl1,segXl2,debug
	integer*4 lastPoint
	real*8 P1X,P1Y,P2X,P2Y,P3Xplus,P3Yplus,P3Xminus,P3Yminus,	&
		MPh1,MPh2,moleScaleFactor,a,b,c,area,AdG,AdG2,slope,PS,PS2,b2_4ac
	real*8 CX1,CY1,CX2,CY2,dist1Plus,dist1Minus,pointXnew(maxPoints),pointYnew(maxPoints)


	segXl1 = segCrystalID(iSeg,1)
	segXl2 = segCrystalID(iSeg,2)
	CX1 = CrystalCenterX(segXl1)
	CY1 = CrystalCenterY(segXl1)
	CX2 = CrystalCenterX(segXl2)
	CY2 = CrystalCenterY(segXl2)
	if(debug.eq.1)then
		write(12,*)' Crystal IDs ',segXl1,segXl2
		write(12,*)' Crystal IDs and Centers - 1:',segXl1,CrystalCenterX(segXl1),CrystalCenterY(segXl1)
		write(12,*)' Crystal IDs and Centers - 2:',segXl2,CrystalCenterX(segXl2),CrystalCenterY(segXl2)
		endif

	lastPoint = numPointEnd(iSeg) - 1
	do iPoint = numPointStart(iSeg)+1,numPointEnd(iSeg)-1		! we don't move the endpoints because they are coincident with the node

		! Get the MIF IDs so we know what the 2 phases are. 
		! These are needed so we can scale the moles to the number of oxygens
		MIFID1 = segMIFID(iSeg,1)
		MPh1 = pointPhaseMolesDelta(iSeg,iPoint,1)
		MIFID2 = segMIFID(iSeg,2)
		MPh2 = pointPhaseMolesDelta(iSeg,iPoint,2)
		! This just scales to moles -- the scaling factor is set when the model is initiated
		MPh1 = MPh1*moleScaleFactor
		MPh2 = MPh2*moleScaleFactor
		if(debug.eq.1)then
			write(12,*)' -----------------------'
			write(12,*)' Seg, point X,Y ',iSeg,iPoint,pointX(iSeg,iPoint),pointY(iSeg,iPoint)
			write(12,*)' MIFID   Phase   deltaMoles     deltaMoles*Scale'
			write(12,103)MIFID1,PhName(MIFID1),pointPhaseMolesDelta(iSeg,iPoint,1),MPh1
			write(12,103)MIFID2,PhName(MIFID2),pointPhaseMolesDelta(iSeg,iPoint,2),MPh2
103			format(T8,I4,2x,A8,2E15.5,5x,I8,2E5.5)
			write(*,*)'Input value for MPh1 -- note that MPh2 = -MPh1 (oxygen balance)'
			read(*,*)MPh1
			write(12,*)' -----------------------'
			endif
		if(MPH1.eq.0)then
			pointXnew(iPoint) = pointX(iSeg,iPoint)
			pointYnew(iPoint) = pointY(iSeg,iPoint)
			cycle		! skip this point if the ∆moles = 0
					! note that MPH1 = -MPH2 (mass conservation)
			endif
		P1X = pointX(iSeg,iPoint)
		P1Y = pointY(iSeg,iPoint)
		if(iPoint.eq.lastPoint)then
			P2X = pointX(iSeg,iPoint-1)		! use the backward point to get the slope
			P2Y = pointY(iSeg,iPoint-1)
			else
			P2X = pointX(iSeg,iPoint+1)		! use the frontward point to get the slope
			P2Y = pointY(iSeg,iPoint+1)
			endif

		if(abs(P1Y-P2Y).lt.1.0d-3)then
			! line is horizontal -- special case
			pointXnew(iPoint) = P1X
			P3Xplus  = P1X
			P3Xminus = P1X
			area = MPH1
			P3Yplus  = P1Y + Area/(P2X-P1X)
			P3Yminus = P1Y - Area/(P2X-P1X)
			
			else
			! line is not horizontal -- calculate new Y value along the slope			
			slope = (P1Y-P2Y)/(P1X-P2X)
			if (abs(slope).lt.1d-20)slope = 1.D-20		! just to ensure we don't have a zero divide
			PS = -1.0d0/slope				! perpendicular slope
			! dGrid = distance between P1 and P2
			area = MPH1
			AdG = area/dXGrid
			AdG2 = AdG**2
			! quadratic formula solution
			PS2 = 1.0d0 + PS**2
			c = -AdG2 + PS2*P1X**2
			b = -2.0d0*P1X*PS2
			a = PS2
			b2_4ac = b**2 - 4.0d0*a*c
			if(b2_4ac.le.0.0d0)then
				! if the value is lt zero then the new point essentially falls on the old point (this should never happen!)
				P3Xplus  = P1X
				P3Xminus = P1X
				P3Yplus  = P1Y
				P3Yminus = P1Y
				else
				P3Xplus =  (-b + sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
				P3Xminus = (-b - sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
				P3Yplus  = P1Y - PS*(P1X - P3Xplus)
				P3Yminus = P1Y - PS*(P1X - P3Xminus)
				endif
			endif
					
		
		! Distance from Crystal center 1 to P3
		! The distance is used to compare with the moles produced to figure out which direction Y must go		
		Dist1plus  = (CX1-P3Xplus)**2  + (CY1-P3Yplus)**2		! no need to take the square root
		Dist1minus = (CX1-P3Xminus)**2 + (CY1-P3Yminus)**2
		if(Dist1Plus.gt.Dist1minus)then
			if(MPH1.gt.0.0d0)then	
				! Dist1Plus is the correct solution
				pointXnew(iPoint) = P3Xplus
				pointYnew(iPoint) = P3Yplus
				if(debug.eq.1)then
					write(12,*)' Dist1Plus > Dist1minus'
					write(12,*)' MPH1 >= 0 '
					write(12,*)'P3plus  ',P3Xplus,P3Yplus		
					write(12,*)'P3minus ',P3Xminus,P3Yminus		
					write(12,*)' Dist1plus  ',Dist1plus
					write(12,*)' Dist1minus ',Dist1minus
					write(12,*)'Plus is correct'
					pause ' Pausing... do the next point'
					endif
				else
				! Dist1minus is the correct solution
				pointXnew(iPoint) = P3Xminus
				pointYnew(iPoint) = P3Yminus
				if(debug.eq.1)then
					write(12,*)' Dist1Plus > Dist1minus'
					write(12,*)' MPH1 < 0 '
					write(12,*)'P3plus  ',P3Xplus,P3Yplus		
					write(12,*)'P3minus ',P3Xminus,P3Yminus		
					write(12,*)' Dist1plus  ',Dist1plus
					write(12,*)' Dist1minus ',Dist1minus
					write(12,*)'Minus is correct'
					pause ' Pausing... do the next point'
					endif
				endif
			else
			if(MPH1.gt.0.0d0)then	
				! Dist1minus is the correct solution
				pointXnew(iPoint) = P3Xminus
				pointYnew(iPoint) = P3Yminus
				if(debug.eq.1)then
					write(12,*)' Dist1Plus < Dist1minus'
					write(12,*)' MPH1 >= 0 '
					write(12,*)'P3plus  ',P3Xplus,P3Yplus		
					write(12,*)'P3minus ',P3Xminus,P3Yminus		
					write(12,*)' Dist1plus  ',Dist1plus
					write(12,*)' Dist1minus ',Dist1minus
					write(12,*)'Minus is correct'
					pause ' Pausing... do the next point'
					endif
				else
				! Dist1plus is the correct solution
				pointXnew(iPoint) = P3Xplus
				pointYnew(iPoint) = P3Yplus
				if(debug.eq.1)then
					write(12,*)' Dist1Plus < Dist1minus'
					write(12,*)' MPH1 < 0 '
					write(12,*)'P3plus  ',P3Xplus,P3Yplus		
					write(12,*)'P3minus ',P3Xminus,P3Yminus		
					write(12,*)' Dist1plus  ',Dist1plus
					write(12,*)' Dist1minus ',Dist1minus
					write(12,*)'Plus is correct'
					pause ' Pausing... do the next point'
					endif
				endif
			endif


				

		end do

	! set the new point positions into the proper array
	do iPoint = numPointStart(iSeg)+1,numPointEnd(iSeg)-1	! we don't move the endpoints because they are coincident with the node
		pointX(iSeg,iPoint) = pointXnew(iPoint)
		pointY(iSeg,iPoint) = pointYnew(iPoint)
		end do

	return
	end
	
			
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	Subroutine RemoveSegPoints(iNode,debug)
	! Routine to remove a segment point when a node point gets too close
	! nodeSegConnect(maxNodes,3),		&! = the segment index that each node is adjacent to (either 2 or 3 segments)
	! nodePointOnTop(maxNodes,3),		&! = index of the coincident point to the node in the point arrays (either first or last) 
	! nodePointNextTo(maxNodes,3),		&! = index of the adjacent point in the point arrays  
	! nodeX(maxNodes),nodeY(maxNodes),	& ! = x,y coordinates of the nodes of the hex grid

	! segX(maxSegs,2),			& !The X,Y values of the endpoints of the segment -- same as the node X,Y
	! segY(maxSegs,2),			&!

	! nodePointOnTop always has the same XY as the node XY
	! As the crystal grows, nodePointNextTo gets close to nodePointOnTop
	! Check the distance between the node (or nodePointOnTop) and the nodePointNextTo
	! If the distance is too small (determined by the value of dXgrid) then
	!     Make nodePointNextTo = the next point in the segment (with either a greater or lesser index)
	!     If the node is at the start of the segment, then this just changes numPointStart(iSeg) -- this index can go negative (see include file)
	!     If the node is at the end   of the segment, then this just changes numPointEnd(iSeg)

	
	implicit none	
	include "Diffuse_GB_Hex.inc"
	include "Assemb.inc"

	integer*4 iNode,iSeg,iPtNextTo,iNode2,kNode2
	integer*4 i,j,k,L,debug,istart,iend,inew
	real*8 length,Ltol,LtolMax,xmid,ymid
			
! 	LtolMax = 2.0d0*dXgrid			
! 	Ltol = 0.4*dXgrid				! dXgrid is read in from the input file (5µm at present)
	LtolMax = dGridAdd*dXgrid			! read in from input file
	Ltol = dGridRemove*dXgrid			! dXgrid and dGridRemove are read in from the input file (5µm at present)
	if(debug.eq.1)then
		write(12,*)'Node ',iNode
		write(12,*)' Ltol (minimum) = ',Ltol
		write(12,*)' LtolMax (max)  = ',LtolMax
		endif
	do i = 1,numNodeSegs(iNode)
		iSeg = nodeSegConnect(iNode,i)			! one of 2 or 3 segments around a node
		iPTNextTo = nodePointNextTo(iNode,i)		! point next to the node along iSeg either at beginning or end of the segment
		length = sqrt((nodeX(iNode)-pointX(iSeg,iPtNextTo))**2 + (nodeY(iNode)-pointY(iSeg,iPtNextTo))**2)
		if(debug.eq.1)then
			write(12,*)'iNode,Seg, PtNextTo,Length ',iNode,iSeg,iPtNextTo,length
			endif
		if(length.lt.Ltol)then				! if dXgrid=5 then the Ltol value is 1 and we should reset the points
			if(abs(numPointEnd(iSeg)-numPointStart(iSeg)).eq.1)cycle	! if there are only 2 points on a segment, don't remove either
			if(debug.eq.1)then
				write(12,*)' Removing a point '
				endif

			if(nodePointOnTop(iNode,i).eq.numPointStart(iSeg))then		!the node connects to the beginning of the segment
				nodePointOnTop(iNode,i)  = nodePointOnTop(iNode,i) + 1
				nodePointNextTo(iNode,i) = nodePointNextTo(iNode,i) + 1
				numPointStart(iSeg) = numPointStart(iSeg) + 1
				segX(iSeg,1) = nodeX(iNode)				! X of segment endpoint = node X
				segY(iSeg,1) = nodeY(iNode)				
				pointX(iSeg,numPointStart(iSeg)) = nodeX(iNode)		! first point in segment is the node
				pointY(iSeg,numPointStart(iSeg)) = nodeY(iNode)
				! set the composition of the nodePointOnTop to the node composition
				do j = 1,numPhCo(1)
					pointComp(iSeg,numPointStart(iSeg),j) = nodeComp(iNode,j)
					end do						
				if(debug.eq.1)then
					write(12,*)'Remove Pt start: iNode Seg Start End onTop nextTo',	&
					iNode,iSeg,numPointStart(iSeg),numPointEnd(iSeg),nodePointOnTop(iNode,i),nodePointNextTo(iNode,i)
					endif
				else								! the node connects to the end of the segment
				nodePointOnTop(iNode,i)  = nodePointOnTop(iNode,i) - 1
				nodePointNextTo(iNode,i) = nodePointNextTo(iNode,i) - 1
				numPointEnd(iSeg) = numPointEnd(iSeg) - 1
				segX(iSeg,2) = nodeX(iNode)				! X of segment endpoint = node X
				segY(iSeg,2) = nodeY(iNode)				
				pointX(iSeg,numPointEnd(iSeg)) = nodeX(iNode)		! last point in segment is the node
				pointY(iSeg,numPointEnd(iSeg)) = nodeY(iNode)
				! set the composition of the nodePointOnTop to the node composition
				do j = 1,numPhCo(1)
					pointComp(iSeg,numPointEnd(iSeg),j)   = nodeComp(iNode,j)
					end do						
				if(debug.eq.1)then
					write(12,*)'Remove Pt End: iNode Seg Start End onTop nextTo',	&
					iNode,iSeg,numPointStart(iSeg),numPointEnd(iSeg),nodePointOnTop(iNode,i),nodePointNextTo(iNode,i)
					endif
				go to 10
				endif
			endif
		if(length.gt.LtolMax)then		! we need to add a segment point
			if(nodePointOnTop(iNode,i).eq.numPointStart(iSeg))then		!the node connects to the beginning of the segment
				if(debug.eq.1)then
					write(12,*)' Adding a point to the beginning of the segment '
					endif
				! the node connects to the beginning of the segment
				! we just need to add a point to the start + 1 of the segment
				! and shift everything else over one point
				! start= 3             4     5     6     7     8   ....... end    -- this becomes
				! start= 3      4      5     6     7     8     9   ....... end+1

				! first shift everything over
				! shift end point over
				! Shift the POT and PNT for the node at the other end of the segment -- only if we're adding at the start
				if(segNodes(iSeg,1).eq.iNode)then
					iNode2 = segNodes(iSeg,2)
					else
					iNode2 = segNodes(iSeg,1)
					endif
				do kNode2 = 1,numNodeSegs(iNode2)
					if(iSeg.eq.nodeSegConnect(iNode2,kNode2))then
						nodePointOnTop(iNode2,kNode2)  = nodePointOnTop(iNode2,kNode2) + 1
						nodePointNextTo(iNode2,kNode2) = nodePointNextTo(iNode2,kNode2) + 1
						go to 12
						endif
					end do
12				continue

				numPointEnd(iSeg) = numPointEnd(iSeg) + 1
				! shift intermediate points over
				istart = numPointStart(iSeg)		! this is the new end point
				iend = numPointEnd(iSeg)
				inew = istart + 1		! index of the point we are adding
!				do k = istart+2,iend WRONG!!!
				do k = iend,istart+2,-1		! need to cycle backwards starting with the end point
					pointX(iSeg,k) = pointX(iSeg,k-1)		
					pointy(iSeg,k) = pointY(iSeg,k-1)		
					do j = 1,numPhCo(1)
						pointComp(iSeg,k,j)    = pointComp(iSeg,k-1,j)
						end do
					do L = 1,2
						pointPhaseMoles(iSeg,k,L) = pointPhaseMoles(iSeg,k-1,L)
						pointPhaseMolesDelta(iSeg,k,L) = pointPhaseMolesDelta(iSeg,k-1,L)
						do j = 1,numPhCo(segMIFID(iSeg,L))
							pointPhaseComp(iSeg,k,L,j) = pointPhaseComp(iSeg,k-1,L,j)
							end do
						end do

					end do
				! this is the new point we're adding
				xmid = (PointX(iSeg,istart+2) + nodeX(iNode))/2.0d0	! midpoint between node and former end-1
				ymid = (PointY(iSeg,istart+2) + nodeY(iNode))/2.0d0
				pointX(iSeg,inew) = xmid		! last point in segment is the node
				pointY(iSeg,inew) = ymid
				! set the composition of the nodePointOnTop to the node composition
				do j = 1,numPhCo(1)
					pointComp(iSeg,inew,j) = nodeComp(iNode,j)
					end do						
				do L = 1,2
					pointPhaseMoles(iSeg,inew,L) = pointPhaseMoles(iSeg,inew-1,L)
					pointPhaseMolesDelta(iSeg,inew,L) = pointPhaseMolesDelta(iSeg,inew-1,L)
					do j = 1,numPhCo(segMIFID(iSeg,L))
						pointPhaseComp(iSeg,inew,L,j) = pointPhaseComp(iSeg,inew-1,L,j)
						end do
					end do


				if(debug.eq.1)then
					write(12,*)' Xmid, Ymid = ',xmid,ymid
					write(12,*)'Add Pt at start: iNode Seg Start end onTop nextTo',	&
					iNode,iSeg,numPointStart(iSeg),numPointEnd(iSeg),nodePointOnTop(iNode,i),nodePointNextTo(iNode,i)
					write(12,*)'    Other end of seg: Node2,kNode2,seg,NPOT,NPNT ',	&
					iNode2,kNode2,nodeSegConnect(iNode2,kNode2),nodePointOnTop(iNode2,kNode2),nodePointNextTo(iNode2,kNode2)
					write(12,*)' ----------'
					endif

				else	
				if(debug.eq.1)then
					write(12,*)' Adding a point to the end of the segment '
					endif
				! the node connects to the end of the segment
				! we just need to add a point to the end of the segment
				! start ......  6     7           end=8    becomes
				! start ......  6     7     8     end=9
				nodePointOnTop(iNode,i)  = nodePointOnTop(iNode,i) + 1
				nodePointNextTo(iNode,i) = nodePointNextTo(iNode,i) + 1
				numPointEnd(iSeg) = numPointEnd(iSeg) + 1		! + 1 expands the range
				!segX(iSeg,2) = nodeX(iNode)				! X of segment endpoint = node X
				!segY(iSeg,2) = nodeY(iNode)				! these don't change
				iend = numPointEnd(iSeg)			! 
				pointX(iSeg,iend) = nodeX(iNode)		! last point in segment is the node
				pointY(iSeg,iend) = nodeY(iNode)
				! this is the new point we're adding
				xmid = (PointX(iSeg,iend-2) + nodeX(iNode))/2.0d0	! midpoint between node and former end-1
				ymid = (PointY(iSeg,iend-2) + nodeY(iNode))/2.0d0
				pointX(iSeg,iend-1) = xmid		! last point in segment is the node
				pointY(iSeg,iend-1) = ymid

				! set the composition of the nodePointOnTop to the node composition
				! leave the composition of nodePointNextTo the same (i.e. point 8 was the pointOnTop and now it's the PointNextTo)
				do j = 1,numPhCo(1)
					! set the composition of the nodePointOnTop to the node composition
					pointComp(iSeg,iend,j)    = nodeComp(iNode,j)	! this is the new point
					! set the composition of the nodePointNextTo (the new one) to the node composition
!					pointComp(iSeg,iend-1,j)   = nodeComp(iNode,j)	! this is not needed
					end do						
				do L = 1,2
					pointPhaseMoles(iSeg,iend,L) = pointPhaseMoles(iSeg,iend-1,L)
					pointPhaseMolesDelta(iSeg,iend,L) = pointPhaseMolesDelta(iSeg,iend-1,L)
					do j = 1,numPhCo(segMIFID(iSeg,L))
						pointPhaseComp(iSeg,iend,L,j) = pointPhaseComp(iSeg,iend-1,L,j)
						end do
					end do

				if(debug.eq.1)then
					write(12,*)' Xmid, Ymid = ',xmid,ymid
					write(12,*)'Add Pt at end: iNode Seg start End onTop nextTo',	&
					iNode,iSeg,numPointStart(iSeg),numPointEnd(iSeg),nodePointOnTop(iNode,i),nodePointNextTo(iNode,i)
					write(12,*)' ----------'
					endif
				endif		
			endif
					
10		continue			

		end do

		

	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
			
	Subroutine SegLengthCalculate()
	! routine to calculate the length of every segment
	! Segment length is the distance between 2 nodes
	! It is used to determine whether 2 nodes will collide (i.e. NodeCapture)
      	implicit none

	include "Diffuse_GB_Hex.inc"
	integer*4 iSeg,node1,node2
	real*8 length
	
	do iSeg = 1,numSegs
		node1 = segNodes(iSeg,1)
		node2 = segNodes(iSeg,2)
		length = sqrt((nodeX(node1)-nodeX(node2))**2 + (nodeY(node1)-nodeY(node2))**2)
		if(length.lt.segLength(iSeg))then
			segGettingShorter(iSeg) = .TRUE.
			else
			segGettingShorter(iSeg) = .FALSE.
			endif
		segLength(iSeg) = length
		end do
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
			
	Subroutine NodeCapture(iunit,debug)
	! routine to reset the info in 2 nodes when they get too close together

      	implicit none
! *****************************************
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
! *****************************************

	integer*4 iSeg,node1,node2,debug,i,k,inode,kcur,j,iPointX,kk,kkk,k1(3),k2(3),PoT1(3),PoT2(3),PnT1(3),PnT2(3),pt,iSegT
	integer*4 MIFID1,MIFID2,MIFID3,MIFID4,Ph1,Ph2,Ph3,Ph4,iPh1N1,iPh1N2,iPh2N1,iPh2N2,iPh3N1,iPh4N2,nodeT(2),MIFID,ph_temp
	integer*4 indexNode1,indexNode2,ksave,iunit,kSeg
	real*8 slope,intcpt,midX,midY,dist,nXplus,NXminus,nYplus,nYminus
	real*8 PS,PI,PS2,Temp,a,b,c,b2_4ac
! 	,P1X,P1Y,P3Xplus,P3Xminus
	real*8 GBmoles,N2XPlus
	integer*4 iSeg1N1,iSeg2N1,iSeg3N1,iSeg1N2,iSeg2N2,iSeg3N2,seg1N1,seg2N1,seg3N1,seg1N2,seg2N2,seg3N2
	integer*4 node1N1,node2N1,node3N1,node1N2,node2N2,node3N2,j1,j2,j3,jT
! 	logical clockwise
	
	! the distance between 2 nodes is stored in segLength(iSeg)
	do iSeg = 1,numSegs

		!if(iSeg.eq.33)write(12,*)iSeg,numMDFRxnsSegs(iSeg),segLength(iSeg),segGettingShorter(iSeg)

		if(segCaptureExclusion(iSeg))cycle		! if .TRUE. this will skip any segment on the exclusion list
		if(segLength(iSeg).lt.dXgrid*0.8)then	! dXgrid is the nominal spacing between points along a segment
			if(segGettingShorter(iSeg))then		! only reset the nodes if the 2 nodes are getting closer together
				! the 2 nodes we are resetting
				node1 = segNodes(iSeg,1)
				node2 = segNodes(iSeg,2)
				if(node1.eq.node2)cycle		! This happens when a phase closes on itself. This is a temporary fix
				write(12,*)'      '
				write(12,*)'      '
				write(12,*)'      '
				write(12,*)' ********** NodeCapture routine *************     '
				write(12,*)'Nodes:  ',node1,node2
				write(12,*)'Cycle number: ',TotalCycles
				write(12,*)' ********************************************'
				write(12,*)'Seg number:       ',iSeg
				write(12,*)'segLength:         ',segLength(iSeg)
				write(12,*)'segGettingShorter: ',segGettingShorter(iSeg)
				write(12,*)' ********************************************'

				write(iunit,*)' ********** NodeCapture routine *************     '
				write(iunit,*)'Nodes:             ',node1,node2
				write(iunit,*)'Cycle number:      ',TotalCycles
				write(iunit,*)'Seg number:        ',iSeg
				write(iunit,*)'segLength:         ',segLength(iSeg)
				write(iunit,*)'segGettingShorter: ',segGettingShorter(iSeg)
				write(iunit,*)' ********************************************'

!	Node capture geometry
! 		see AI drawing for labeling of phases around each node

!	Nodes are moved in an anticlockwise direction
!
!	\                       /
!        \        Ph2          /
!         \                   /	
!   Ph3  N1o-----------------o N2     Ph4
!         /                   \
!        /        Ph1          \
!       /                       \

!
!          \  Ph2 /
!           \    /
!            \  /
!              o N2
!              |
!     Ph3      |    Ph4
!              |
!              o N1
!             / \
!            /   \
!           /     \
!          /  Ph1  \



		! We will locate each node along a line that perpendicular to the line that connects the 2 nodes initially
		!
		!                  N2f
		!                   |
		!                   |
		!                   |
		!     N1i------------------------N2i
		!                   |
		!                   |
		!                   |
		!                   N1f
		!
		! we will place N1f---N2f just a bit over dXgrid apart

	! what needs resetting (for each node)?
	! Node changes
	! NodeX,NodeY
	! numNodeReactionPhases(maxNodes),	& ! the number of different reacting phases at a node
	! nodeReactionPhases(maxNodes,3),		& ! index of the reacting phases around the node 
						  ! if numNodeReactionPhases = 0 this is irrelevent
						  ! if numNodeReactionPhases = 3, then the indicies are 1, 2, 3
						  ! if numNodeReactionPhases=2 the indicies are either 1&2, 1&3 or 2&3
	! numMDFRxnsNodes(maxNodes),		& ! number of MDF reactions at a node - calculated on start-up. If numMDFRxnsNodes =0 then an equilibrium model is used
	! nodeSegTheSamePhases(maxNodes),		& ! This is the segment along which the 2 phases are the same for nodes where 2 phases are the same (set in ReadModelFile line 500 or so)
	! nodeMIFID(maxNodes,3),			& ! the 3 phases (MIF ID) at the node
	! nodeCrystalIndex(maxNodes,3)		 ! the crystal number (index) for each of the 3 crystals at the node
	! nodePointOnTop(maxNodes,3),		&! = index of the coincident point to the node in the point arrays (either first or last) 
	! nodePointNextTo(maxNodes,3),		&! = index of the adjacent point in the point arrays  

	! stays the same
	! nodeNodeConnect(maxNodes,3), 		&! = index of the nodes that connect to 2 or 3 other nodes
	! nodeSegConnect(maxNodes,3),		&!   = index of the segments connected to a node. Each node can have up to 3 segments attached -- counterclockwise order
	! nodeSegEnd(maxNodes,3),		&!  1 = start at beginning of segment, 2 = start at end of segment (no longer used)
	
	! seg changes
	! segX(iSeg,2)
	! segY(iSeg,2)
	! pointX(maxSegs,minPoints:maxPoints),	&!
	! pointY(maxSegs,minPoints:maxPoints),	& ! = x,y coordinates of each point along a given segment
	! segLength(iSeg)
	! segCrystalID(iSeg,2) must change
		! After the reset, they are the 2 "other" crystals involved in N1i and N2i
	! segMIFID(iSeg,2) must change with the crystals
	! numMDFRxnsSegs(iSeg) must change
	! numSegReactionPhases(iSeg) must change

	! numCrystalNodes(crystal number)
	! CrystalSegs
	! CrystalNodes

	! nodePhaseComp(maxNodes,3,maxPhCo),			& ! Composition of each of 3 phases at a reaction node
	! nodePhaseMoles(maxNodes,3),				&! Moles of each phase at a reaction node - cumulative
	! nodePhaseMolesDelta(maxNodes,3),			&! Moles of each phase at a reaction node - Only the last set of calculations - temporary to use for crystal growth

	! The tricky part is figuring out which phase is which. Each node has 3 phases around it, but we don't know which is which to start
	! We'll use the crystal numbers, which are unique for each node and segment.

	! The Crystals that will be common across the segment are the 2 crystals that are NOT common across the segment initially
	
	! reset the X,Y for the 2 nodes
	! first get the equation of the line connecting the 2 nodes
	if(abs(nodeX(node2)-nodeX(node1)).lt.1d-10)then		! slope is nearly vertical
		slope = 1.0d10
		else
		slope = (nodeY(node2)-nodeY(node1))/(nodeX(node2)-nodeX(node1))
		if(abs(slope).lt.1.0d-20)slope = 1.0d-20		! in case slope is horizontal
		endif	
	intcpt = nodeY(node2) - slope*nodeX(node2)
	! find the midpoint
	midX = (nodeX(node2)+nodeX(node1))/2.0d0
	midY = (nodeY(node2)+nodeY(node1))/2.0d0
!	Dist = dXgrid/2.0d0 + .1		! set the distance to just a bit longer than the nominal grid spacing
	!Dist = dXgrid/10.0d0 + .1		! set the distance to something small
	Dist = dXgrid*0.4		! This is the distance the new nodes will be apart

	if(Abs(nodeX(node2)-nodeX(node1)).lt.1.0d-3)then		! original segment is nearly vertical -- special case
		if(debug.eq.1)write(12,*)' Nearly vertical'
		nXPlus = midX + dist
		nYPlus = midY
		nXMinus = midX - dist
		nYMinus = midY
		! we need to figure out which goes with which node below
		go to 10
		endif
	if(Abs(nodeY(node2)-nodeY(node1)).lt.1.0d-3)then		! original segment is nearly horizontal -- special case
		if(debug.eq.1)write(12,*)' Nearly horizontal'
		nXPlus = midX
		nYPlus = midY + dist
		nXMinus = midX
		nYMinus = midY - dist
		! we need to figure out which goes with which node below
		go to 10
		endif
	! if here, then the new nodes are along a line that is not horizontal or vertical
	! the code here is identical to that used for growing segments.
	! Get the line that is perpendicular to the line connecting the 2 nodes 
	! and solve for points that are the distance dXgrid/2 above and below the midpoint
	! This requires solving a quadratic in X
	PS = -1.0d0/slope		! slope of the line perpendicular to the original node1-node2 line
	PI = midY - PS*midX		! intercept of the perpendicular line
	! quadratic formula solution
!	0 = a*X**2 + b*X + c	
	PS2 = 1.0d0 + PS**2
!	temp = PI + MidY
	temp = PI - MidY
!	c = temp**2 - dist**2		! constant
	c = midX**2 + temp**2 - dist**2		! constant
!	b = -2.0d0*(midX + PS*temp)	
	b = 2.0d0*(PS*temp - midX)	
	a = PS2
	b2_4ac = b**2 - 4.0d0*a*c
	if(b2_4ac.le.0.0d0)then
		! if the value is lt zero then the new point essentially falls on the old point (this should never happen!)
		write(12,*)' (b^2 - 4ac) <= 0 -- this should not happen'
		write(12,*)'(b^2 - 4ac)  ', b2_4ac
		write(12,*)' a = ',a
		write(12,*)' b = ',b
		write(12,*)' c = ',c
		write(12,*)' Node1 ',node1,nodeX(node1),nodeY(node1)
		write(12,*)' Node2 ',node2,nodeX(node2),nodeY(node2)
		write(12,*)' MidXY ',midX,midY
		write(12,*)' Slope,Int  ',slope,intcpt
		pause 'Hit return to continue'
		else
		nXplus =  (-b + sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
		nXminus = (-b - sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
! 		nYplus  = P1Y - PS*(P1X - P3Xplus)
! 		nYminus = P1Y - PS*(P1X - P3Xminus)
		nYplus  = midY - PS*(midX - nXplus)
		nYminus = midY - PS*(midX - nXminus)
		endif

10	continue
	if(debug.eq.1)then
		write(12,*)
		write(12,*)'----- XY calculations --------'
		write(12,*)' Original XY '
		write(12,*)' Node1 ',node1,nodeX(node1),nodeY(node1)
		write(12,*)' Node2 ',node2,nodeX(node2),nodeY(node2)
		write(12,*)' MidXY ',midX,midY
		write(12,*)' Slope,Int  ',slope,intcpt
		write(12,*)'Xplus,YPlus   ',nXPlus,nYPlus
		write(12,*)'Xminus,Yminus ',nXMinus,nYMinus
		endif

	! we now have X,Y for the new node positions
	! we need to figure out which is which	

	! see AI drawing for labeling of phases around each node

	! The code requires that the sequence of phases P1, P2 and P3 around node1 is counterclockwise
	! Set the sequence of phases around the node.

	do i = 1,3
		! first find the node index (i) that matches the first segCrystal index
		! This code ensures that ph1 --> ph2 is counterclockwise around node1
		! note that sometimes this will be phase 1 in iSeg and sometimes it will be phase 2 in iSeg
		if(nodeCrystalIndex(node1,i).eq.segCrystalID(iSeg,1))then
			k = i + 1
			if(k.gt.3)k = 1
			if(nodeCrystalIndex(node1,k).eq.segCrystalID(iSeg,2))then		
				! we have counterclockwise sequence 1 (ph1-ph2-ph3)
				! note that these are crystal numbers, not indexes
				ph1 = segCrystalID(iSeg,1)
				ph2 = segCrystalID(iSeg,2)
				else
				ph1 = segCrystalID(iSeg,2)
				ph2 = segCrystalID(iSeg,1)
				endif
			go to 25
			endif
		end do
25	continue

! 	MIFID1 = CrystalMIFID(Ph1)
! 	MIFID2 = CrystalMIFID(Ph2)

!	Determine the indexes for phases Ph1 and Ph2 around nodes 1 and 2
	do i = 1,3
		if(nodeCrystalIndex(node1,i).eq.ph1)iPh1N1 = i
		if(nodeCrystalIndex(node2,i).eq.ph1)iPh1N2 = i
		if(nodeCrystalIndex(node1,i).eq.ph2)iPh2N1 = i
		if(nodeCrystalIndex(node2,i).eq.ph2)iPh2N2 = i
		end do

!	Determine the indexes for phases Ph3 and Ph4
!	Do Ph3
	do i = 1,3	!loop on the 3 crystals around Node1
		ph_temp = nodeCrystalIndex(node1,i)
		if(ph_temp.ne.Ph1.and.ph_temp.ne.Ph2)then	! we want the phase that is not common across the segment
			ph3 = ph_temp
			iph3N1 = i		! index for this phase for node1
			MIFID3 = CrystalMIFID(Ph3)
			go to 20
			endif
		end do
	write(12,*)'Could not match the crystals around node1 across iSeg '
	write(12,*)'node1, iSeg ',node1,iSeg
	write(12,*)'Ph1, Ph2 ',Ph1,Ph2
	write(12,*)'Node Ph  ',(nodeCrystalIndex(node1,i),i=1,3)
	pause 'Hit return to continue'
20	continue
!	Do Ph4
	do i = 1,3	!loop on the 3 crystals around node2
		ph_temp = nodeCrystalIndex(node2,i)
		if(ph_temp.ne.Ph1.and.ph_temp.ne.Ph2)then	! we want the phase that is not common across the segment
		!		if(nodeCrystalIndex(node2,i).eq.Ph1.or.nodeCrystalIndex(node2,i).eq.Ph2)cycle
			ph4 = ph_temp
			iph4N2 = i		! index for this phase for node2
			MIFID4 = CrystalMIFID(Ph4)
			go to 22
			endif
		end do
	write(12,*)'Could not match the crystals around node2 across iSeg '
	write(12,*)'node2, iSeg ',node2,iSeg
	write(12,*)'Ph1, Ph2 ',Ph1,Ph2
	write(12,*)'Node Ph  ',(nodeCrystalIndex(node2,i),i=1,3)
	pause 'Hit return to continue'
22	continue

!	Determine the identities and the indexes of the of the segments around each node
!	iSeg will be seg1 in each case (seg1N1 and seg1N2)
!	Note that this code assumes the segments are counterclockwise around each node
	! do node 1
	do i = 1,3
		if(nodeSegConnect(node1,i).eq.iSeg)then
			k = i
			iseg1N1 = i
			seg1N1 = iSeg
			node1N1 = node2		! iSeg connects node1 and node2
			k = k + 1
			if(k.gt.3)k = 1
			seg2N1 = nodeSegConnect(node1,k)
			iseg2N1 = k
			node2N1 = nodeNodeConnect(node1,k)
			k = k + 1
			if(k.gt.3)k = 1
			seg3N1 = nodeSegConnect(node1,k)
			iseg3N1 = k
			node3N1 = nodeNodeConnect(node1,k)
			endif
		end do	
	! do node 2
	do i = 1,3
		if(nodeSegConnect(node2,i).eq.iSeg)then
			k = i
			iseg1N2 = i
			seg1N2 = iSeg
			node1N2 = node1		! iSeg connects node1 and node2
			k = k + 1
			if(k.gt.3)k = 1
			seg2N2 = nodeSegConnect(node2,k)
			iseg2N2 = k
			node2N2 = nodeNodeConnect(node2,k)
			k = k + 1
			if(k.gt.3)k = 1
			seg3N2 = nodeSegConnect(node2,k)
			iseg3N2 = k
			node3N2 = nodeNodeConnect(node2,k)
			endif
		end do	




	if(debug.eq.1)then
		write(12,*)'  -----Before resetting ----' 
		write(12,*)'node1 ',node1
		write(12,*)'      Xl index     ',(nodeCrystalIndex(node1,i),i=1,3)
		write(12,*)'      MIFID        ',(nodeMIFID(node1,i),i=1,3)
		write(12,*)' node seg connect  ',(nodeSegConnect(node1,i),i=1,3)
		write(12,*)' node node connect ',(nodeNodeConnect(node1,i),i=1,3)
		write(12,*)' Reindexed starting with iSeg'
		write(12,*)'   indexes         ',iseg1N1,iseg2N1,iseg3N1
		write(12,*)'   segs            ',seg1N1,seg2N1,seg3N1
		write(12,*)'   nodes           ',node1N1,node2N1,node3N1
		write(12,*)' '
		write(12,*)'node2 ',node2
		write(12,*)'      Xl index     ',(nodeCrystalIndex(node2,i),i=1,3)
		write(12,*)'      MIFID        ',(nodeMIFID(node2,i),i=1,3)
		write(12,*)' node seg connect  ',(nodeSegConnect(node2,i),i=1,3)
		write(12,*)' node node connect ',(nodeNodeConnect(node2,i),i=1,3)
		write(12,*)' Reindexed starting with iSeg'
		write(12,*)'   indexes         ',iseg1N2,iseg2N2,iseg3N2
		write(12,*)'   segs            ',seg1N2,seg2N2,seg3N2
		write(12,*)'   nodes           ',node1N2,node2N2,node3N2
		write(12,*)' '
		write(12,*)'Ph1, Ph2     ',Ph1,Ph2
		write(12,*)'Ph3, Ph4     ',Ph3,Ph4
		write(12,*)' '
		write(12,*)'        iSeg ',iSeg
		write(12,*)' segXlID     ',segCrystalID(iSeg,1),segCrystalID(iSeg,2)
		write(12,*)' segMIDFID   ',segMIFID(iSeg,1),segMIFID(iSeg,2)

		write(12,*)' '
		write(12,*)' '
		write(12,*)' Initial values'
		write(12,*)' '
		write(12,*)' CrystalNodes and CrystalSegs'
		write(12,*)'Ph1   = ',Ph1
		write(12,*)'iPh1N1 = ',iPh1N1
		write(12,*)'iPh1N2 = ',iPh1N2
		write(12,*)'Nodes  Segs'	
		do k = 1,numCrystalNodes(Ph1)
			write(12,*)CrystalNodes(Ph1,k),CrystalSegs(Ph1,k)
			end do
		write(12,*)'Ph2 = ',Ph2
		write(12,*)'iPh2N1 = ',iPh2N1
		write(12,*)'iPh2N2 = ',iPh2N2
		write(12,*)'Nodes  Segs'	
		do k = 1,numCrystalNodes(Ph2)
			write(12,*)CrystalNodes(Ph2,k),CrystalSegs(Ph2,k)
			end do
		write(12,*)'Ph3 =    ',Ph3
		write(12,*)'iPh3N1 = ',iPh3N1
		write(12,*)'Nodes  Segs'	
		do k = 1,numCrystalNodes(Ph3)
			write(12,*)CrystalNodes(Ph3,k),CrystalSegs(Ph3,k)
			end do
		write(12,*)'Ph4 =    ',Ph4
		write(12,*)'iPh4N2 = ',iPh4N2
		write(12,*)'Nodes  Segs'	
		do k = 1,numCrystalNodes(Ph4)
			write(12,*)CrystalNodes(Ph4,k),CrystalSegs(Ph4,k)
			end do
		endif


! 	! we use the distance between the nodePlus and nodeMinus locations and the crystal centers to determine which is which
! 	! this code will put node1 closer to the crystal center for Ph1
! 	CXPhn = CrystalCenterX(Ph1)
! 	CYPhn = CrystalCenterY(Ph1)
! 	distPhnPlus  = (nXplus - CXPhn)**2  + (nYplus - CYPhn)**2		! no need to take the square root
! 	distPhnMinus = (nXMinus - CXPhn)**2 + (nYMinus - CYPhn)**2
! 	if(distPhnPlus.lt.distPhnMinus)then	
! 		! Xplus,Yplus are the new locations of node1
! 		nodeX(node1) = nXplus
! 		nodeY(node1) = nYplus
! 		! node2 is Xminus,Yminus by default
! 		nodeX(node2) = nXminus
! 		nodeY(node2) = nYminus
! 		else
! 		! Xplus,Yplus are the new locations of node2
! 		nodeX(node2) = nXplus
! 		nodeY(node2) = nYplus
! 		! node1 is Xminus,Yminus by default
! 		nodeX(node1) = nXminus
! 		nodeY(node1) = nYminus
! 		endif
!	The above code is TERRIBLE and fails with crystals of odd shape

!	We want the new node to be anticlockwise from the old node.
!	The condition to ensure this is that the cross product of the 2 vectors (center-->Node2 and center-->Node2new) is positive
!	We need to first subtract the center coordinates (midX,midY) so the 2 vectors emanate from the origin
! 	nXplusO = nXplus - midX etc.
! 	nYplusO = nYplus - midY
! 	nXminusO = nXminus - midX
! 	nYminusO = nYminus - midY
!	now do the cross product. Note that the 3rd component of each vector (the "z" component) is always 0 because this is a 2-D plane
!	so we only need to compute the third component of the cross product (the other 2 are both zero)
	N2xPlus = (nodeX(node2)-midX)*(nYplus-midY) - (nodeY(node2)-midY)*(nXplus-midX)
!	if N2xPlus > 0 then Plus mode is clockwise from node2. Otherwise it is counterclockwise
	if(N2xPlus.gt.0.)then
		nodeX(node2) = nXplus
		nodeY(node2) = nYplus
		! node1 is Xminus,Yminus by default
		nodeX(node1) = nXminus
		nodeY(node1) = nYminus
		else		
		! Xplus,Yplus are the new locations of node1
		nodeX(node1) = nXplus
		nodeY(node1) = nYplus
		! node2 is Xminus,Yminus by default
		nodeX(node2) = nXminus
		nodeY(node2) = nYminus
		endif


	! Adjust the variables for the nodes
	nodeT(1) = node1	! temporary storage for node IDs for do loops
	nodeT(2) = node2


	! Noe that if one of the 2 nodes has only 2 Node reaction phases (ie numNodeReactionPhases(iNode) = 2)
	!   Then the compositions are not set for the 3rd phase at this node (only the 2 reacting phases are actually set)
	!   To account for this, we need to set the compositions of the 3rd (non-reacting) phase at this node
	do i = 1,2
		iNode = nodeT(i)
! 		if(i.eq.1)then
! 			iNode = node1
! 			else
! 			iNode = node2
! 			endif
		if(numNodeReactionPhases(iNode).eq.2)then
			! figure out which phase needs to be set
			j1 = nodeReactionPhases(iNode,1)	! should always be = 1
			j2 = nodeReactionPhases(iNode,2)	! may be 2 or 3
			if(j2.eq.2)then
				j3 = 3		! we need to set phase 3
				else
				j3 = 2		! we need to set phase 2
				endif
			! Which phase is the same?
			if(nodeMIFID(iNode,j3).eq.nodeMIFID(iNode,j1))then
				jT = j1
				else
				jT = j2
				endif
			do j = 1,numPhCo(nodeMIFID(iNode,jT))
				nodePhaseComp(iNode,j3,j) = nodePhaseComp(iNode,jT,j)
				end do
			nodePhaseMoles(iNode,j3) = nodePhaseMoles(iNode,jT)
			nodePhaseMolesDelta(iNode,j3) = nodePhaseMolesDelta(iNode,jT)

			endif
		end do

				


	! swap phases at node1 Ph2--> Ph4
	! swap phase compositions Ph2 at node1 is the former Ph4 at node2
	! iPh2N1 is the index for Ph2 at node 1 
	! iPh4N2 is the index for Ph4 at node 2 
	do j = 1,numPhCo(nodeMIFID(node2,iPh4N2))
		nodePhaseComp(node1,iPh2N1,j) = nodePhaseComp(node2,iPh4N2,j)
		end do
	! swap phases at  node2 Ph1--> Ph3
	! swap phase compositions Ph1 at node2 is the former Ph3 at node1
	! iPh1N2 is the index for Ph1 at node 2 
	! iPh3N1 is the index for Ph3 at node 1 
	do j = 1,numPhCo(nodeMIFID(node1,iPh3N1))
		nodePhaseComp(node2,iPh1N2,j) = nodePhaseComp(node1,iPh3N1,j)
		end do


	! Swap phase numbers and MIFIDs
	nodeCrystalIndex(node1,iPh2N1) = Ph4
	nodeMIFID(node1,iPh2N1) = CrystalMIFID(Ph4)
	nodeCrystalIndex(node2,iPh1N2) = Ph3
	nodeMIFID(node2,iPh1N2) = CrystalMIFID(Ph3)
	


	! Adjust CrystalSegs and CrystalNodes
	! Note that CrystalNodes and CrystalSegs go clockwise from start around the entire crystal
	! Ph1 and Ph2 both lose a seg (iSeg) and a node
	! which node is lost depends on whether the sequence is clockwise or counterclockwise such that
	!    node 1 moves towards ph1 or ph2
	! The sequence should be counterclockwise as set up when the grid is created

	! Ph1 loses Node2
	do k = 1,numCrystalNodes(Ph1)
		if(CrystalNodes(Ph1,k).eq.node1)then	! we lose the next node
			if(k.eq.numCrystalNodes(Ph1))then	! this is the last one -- just get rid of the last one
				kk = 1
				kkk = k-1
				else
				kk = k+1
				kkk = k
				endif
			j = CrystalSegs(Ph1,kk)			! temporary storage for the segment next clockwise. This is the segment that goes with the
		write(12,*)'+++++++++++++++'
		write(12,*)Ph1,node1,k,kk,j
		write(12,*)'+++++++++++++++'
			do i = kk,numCrystalNodes(Ph1)-1
				CrystalNodes(Ph1,i)=CrystalNodes(Ph1,i+1)
				CrystalSegs(Ph1,i) =CrystalSegs(Ph1,i+1)
				end do
			CrystalSegs(Ph1,kkk) = j	! this is the seg clockwise from the next (old) node

			go to 34
			endif
		end do		
34	continue
	numCrystalNodes(Ph1) = numCrystalNodes(Ph1) - 1

	! Ph2 loses Node1
	do k = 1,numCrystalNodes(Ph2)
		if(CrystalNodes(Ph2,k).eq.node2)then	! we lose the next node
			if(k.eq.numCrystalNodes(Ph2))then	! this is the last one -- just get rid of the last one
				kk = 1
				kkk = k-1
				else
				kk = k + 1
				kkk = k
				endif
			j = CrystalSegs(Ph2,kk)			! temporary storage for the segment next clockwise. This is the segment that goes with the
		write(12,*)'+++++++++++++++'
		write(12,*)Ph2,node2,k,kk,j
		write(12,*)'+++++++++++++++'
			!CrystalSegs(Ph2,k)=CrystalSegs(Ph2,k+1)	! this is the seg clockwise from the next (old) node
			do i = kk,numCrystalNodes(Ph2)-1
				CrystalNodes(Ph2,i)=CrystalNodes(Ph2,i+1)
				CrystalSegs(Ph2,i)=CrystalSegs(Ph2,i+1)
				end do
			CrystalSegs(Ph2,kkk) = j	! this is the seg clockwise from the next (old) node
			go to 35
			endif
		end do		
35	continue
	numCrystalNodes(Ph2) = numCrystalNodes(Ph2) - 1



	! Ph3 and Ph4 both gain a seg (iSeg) and a node
	! Ph3 gains node2
	do k = 1,numCrystalNodes(Ph3)
		if(CrystalNodes(Ph3,k).eq.node1)then	! we gain inode2 after this node
			ksave = k
			do i = numCrystalNodes(Ph3)+1,k+1,-1
				CrystalNodes(Ph3,i)=CrystalNodes(Ph3,i-1)
				CrystalSegs(Ph3,i) =CrystalSegs(Ph3,i-1)
				end do
			CrystalNodes(Ph3,ksave)=node2
			CrystalSegs(Ph3,ksave)=iSeg
			go to 38
			endif
		end do		
38	continue
	numCrystalNodes(Ph3) = numCrystalNodes(Ph3) + 1
	! Ph4 gains node1
	do k = 1,numCrystalNodes(Ph4)
		if(CrystalNodes(Ph4,k).eq.node2)then	! we gain inode1 after this node
			ksave = k
			do i = numCrystalNodes(Ph4)+1,k+1,-1
				CrystalNodes(Ph4,i)=CrystalNodes(Ph4,i-1)
				CrystalSegs(Ph4,i)=CrystalSegs(Ph4,i-1)
				end do
			CrystalNodes(Ph4,ksave)=node1
			CrystalSegs(Ph4,ksave)=iSeg
			go to 39
			endif
		end do		
39	continue
	numCrystalNodes(Ph4) = numCrystalNodes(Ph4) + 1

	if(debug.eq.1)then
		write(12,*)' '
		write(12,*)' ---------------------'
		write(12,*)' After node capture '
		write(12,*)' CrystalNodes and CrystalSegs'
		write(12,*)'Ph1 = ',Ph1
		write(12,*)'Nodes  Segs'	
		do k = 1,numCrystalNodes(Ph1)
			write(12,*)CrystalNodes(Ph1,k),CrystalSegs(Ph1,k)
			end do
		write(12,*)'Ph2 = ',Ph2
		write(12,*)'Nodes  Segs'	
		do k = 1,numCrystalNodes(Ph2)
			write(12,*)CrystalNodes(Ph2,k),CrystalSegs(Ph2,k)
			end do
		write(12,*)'Ph3 = ',Ph3
		write(12,*)'Nodes  Segs'	
		do k = 1,numCrystalNodes(Ph3)
			write(12,*)CrystalNodes(Ph3,k),CrystalSegs(Ph3,k)
			end do
		write(12,*)'Ph4 = ',Ph4
		write(12,*)'Nodes  Segs'	
		do k = 1,numCrystalNodes(Ph4)
			write(12,*)CrystalNodes(Ph4,k),CrystalSegs(Ph4,k)
			end do
		endif
			


	
	! This code is directly from sub SetUpCrystals in GB_ReadModelFile.f90 around line 515
	! Assume there are 3 phases at a node: the options are
	! All phases the same (no reaction) -- not likely with node capture
	! All phases different (3-phase reaction)
	! 2 phases the same -- either 1 & 2 or 1 & 3 (2-phase reaction)
	do i = 1,2
		iNode = nodeT(i)
		call AssignReactionPhases(iNode)
		go to 444

		MIFID1 = nodeMIFID(iNode,1)
		MIFID2 = nodeMIFID(iNode,2)
		MIFID3 = nodeMIFID(iNode,3)
		if(MIFID1.eq.MIFID2.and.MIFID1.eq.MIFID3)then	! all phases are the same
			numNodeReactionPhases(iNode) = 0
			go to 23
			endif
		if(MIFID1.ne.MIFID2.and.MIFID1.ne.MIFID3.and.MIFID2.ne.MIFID3)then	! all 2 phases different
			numNodeReactionPhases(iNode) = 3
			nodeReactionPhases(iNode,1) = 1
			nodeReactionPhases(iNode,2) = 2
			nodeReactionPhases(iNode,3) = 3
			nodeSegTheSamePhases(iNode) = 0
			go to 23
			endif	
		! if here, then 2 phases must be the same
		numNodeReactionPhases(iNode) = 2
		if(MIFID1.eq.MIFID2)then
			nodeReactionPhases(iNode,1) = 1
			nodeReactionPhases(iNode,2) = 3
			do k = 1,3
				kSeg = nodeSegConnect(iNode,k)
				if(segMIFID(kSeg,1).eq.SegMIFID(kSeg,2))then
					nodeSegTheSamePhases(iNode) = k	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
					go to 23
					endif
				end do
! 			go to 23
			endif
		if(MIFID1.eq.MIFID3)then
			nodeReactionPhases(iNode,1) = 2
			nodeReactionPhases(iNode,2) = 3
			do k = 1,3
				kSeg = nodeSegConnect(iNode,k)
				if(segMIFID(kSeg,1).eq.SegMIFID(kSeg,2))then
					nodeSegTheSamePhases(iNode) = k	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
					go to 23
					endif
				end do
! 			go to 23
			endif
		if(MIFID2.eq.MIFID3)then
			nodeReactionPhases(iNode,1) = 1
			nodeReactionPhases(iNode,2) = 2
			do k = 1,3
				kSeg = nodeSegConnect(iNode,k)
				if(segMIFID(kSeg,1).eq.SegMIFID(kSeg,2))then
					nodeSegTheSamePhases(iNode) = k	! this is the segment along which the 2 phases are identical nodeSegTheSamePhases(iNode)
					go to 23
					endif
				end do
! 			go to 23
			endif
		! if here, we have a problem
		write(12,*)' In routine OpenModelOutputFile'
		write(12,*)' Failure to determine nodeReactionPhases and nodeSegTheSamePhases'
		write(12,*)' iNode = ',iNode
		write(12,227)iNode,(nodeMIFID(iNode,j),j=1,numCrystalsAtNode(iNode)),			&
				(nodeReactionPhases(iNode,j),j=1,numNodeReactionPhases(iNode)),nodeSegTheSamePhases(iNode)
	227	format(20I8)
	
	23	continue	

444		continue	! this is just a dummy jump to see if new subroutine works
		
		end do

	do k = 1,2
!		if(numNodeReactionPhases(iNode).eq.0)cycle	! there is no reaction at this node, so segs don't change
		iNode = nodeT(k)
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
			kk = asmCurrent(kCur)
			NP = numPhCo(kk)+NP
			end do
      		NX=NP-numPh		! This is the number of independent compositional variables (first X in each phase is dependent)
!		call SetALLX()		! set up ALLX array
		call Names()		! sets up names of variables for asmCurrent
      		call REXN		! calculate linearly independent reactions
!		numMDFRxnsNodes(iNode) = NEQ		!NEQ is either the number of EQ reactions OR number of MSDF reactions
		numMDFRxnsNodes(iNode) = numMDFrxn
		end do


	! Now adjust the segment that connects the 2 nodes
	! Set the XY seg endpoints
	! This code sets the endpoints for the segment that links the 2 nodes that are swapped (captured)
	if(segNodes(iSeg,1).eq.node1)then
		segX(iSeg,1) = nodeX(node1)
		segY(iSeg,1) = nodeY(node1)
		pointX(iSeg,numPointStart(iSeg)) = SegX(iSeg,1)		
		pointY(iSeg,numPointStart(iSeg)) = SegY(iSeg,1)		
		segX(iSeg,2) = nodeX(node2)
		segY(iSeg,2) = nodeY(node2)
		pointX(iSeg,numPointEnd(iSeg)) = SegX(iSeg,2)		
		pointY(iSeg,numPointEnd(iSeg)) = SegY(iSeg,2)		
		else
		segX(iSeg,1) = nodeX(node2)
		segY(iSeg,1) = nodeY(node2)
		pointX(iSeg,numPointStart(iSeg)) = SegX(iSeg,2)		
		pointY(iSeg,numPointStart(iSeg)) = SegY(iSeg,2)		
		segX(iSeg,2) = nodeX(node1)
		segY(iSeg,2) = nodeY(node1)
		pointX(iSeg,numPointEnd(iSeg)) = SegX(iSeg,1)		
		pointY(iSeg,numPointEnd(iSeg)) = SegY(iSeg,1)		
		endif

	segLength(iSeg) = sqrt((nodeX(node1)-nodeX(node2))**2 + (nodeY(node1)-nodeY(node2))**2)
	segGettingShorter(iSeg) = .FALSE.

	segCrystalID(iSeg,1) = Ph3		
	segCrystalID(iSeg,2) = Ph4		
	segMIFID(iSeg,1) = CrystalMIFID(Ph3)
	segMIFID(iSeg,2) = CrystalMIFID(Ph4)
	
	if(segMIFID(iSeg,1).ne.segMIFID(iSeg,2))then		! they are not the same phase so there is a reaction
		numSegReactionPhases(iSeg) = 2
		else
		numSegReactionPhases(iSeg) = 0			! the phases are the same so there is no reaction
		numMDFRxnsSegs(iSeg) = 0	
		endif
	if(numSegReactionPhases(iSeg).eq.2)then			! there are 2 phases and hence a reaction. Set the appropriate variables
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
		call Names()		! sets up names of variables for asmCurrent
		call REXN		! calculate linearly independent reactions
		!numMDFRxnsSegs(iSeg) = NEQ	
		numMDFRxnsNodes(iSeg) = numMDFrxn
		! set the phase compositions for the 2 phases at the end points of the segment
		! We get the phase compositions from the values in the nodes
		! We just need to be sure to use the correct phase

		! Loop for the first reaction phase (iSeg,1)
		do i = 1,3								! loop through all 3 phases around node1 (node2 has the same 3 phases
			if(segMIFID(iSeg,1).eq.nodeMIFID(node1,i))then			!find the phase with the matching MIFID
				do iPointX = numPointStart(iSeg),numPointEnd(iSeg)	! do all points along this segment
					do j = 1,numPhCo(segMIFID(iSeg,1))		! do all phase components for this phase
						pointPhaseComp(iSeg,iPointX,1,j) = nodePhaseComp(node1,i,j)	!The phase composition at node1 is the same as at node2			
						end do
					end do
				end if
			end do
		! Loop for the second reaction phase (iSeg,2)
		do i = 1,3	
			if(segMIFID(iSeg,2).eq.nodeMIFID(node1,i))then
				do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
					do j = 1,numPhCo(segMIFID(iSeg,2))
						pointPhaseComp(iSeg,iPointX,2,j) = nodePhaseComp(node1,i,j)				
						end do
					end do
				end if
			end do

		endif

	! I also need to reset the segment endpoints for 2 other segments that attach to these 2 nodes
	!    One gets removed and one gets added to each node

	! Adjust the other 2 segments around each node so the endpoints coincide with the node position
	! Note that the connecting segs switch
	! Note that segs are ordered counter counterclockwise around a node
	! node1
	! 	the new segs are iSeg, iSeg-1 from node2 and iSeg+1 from node1
	do indexNode1 = 1,3		! determine the index of iSeg around node1
		if(nodeSegConnect(Node1,indexNode1).eq.iSeg)then
			go to 40
			endif
		end do
		call fss_alert('ALERT','Did not find the segment for node1')
40	continue
	do indexNode2 = 1,3		! determine the index of iSeg around node1
		if(nodeSegConnect(Node2,indexNode2).eq.iSeg)then
			go to 41
			endif
		end do
		call fss_alert('ALERT','Did not find the segment for node2')
41	continue

	! reset node 1 segments
	! note that the segments are ordered anti-clockwise around the node
	k = indexNode1		! This is the index of the segment iSeg around the node (1,2 or 3)
	k1(k) = iSeg
	PoT1(k) = nodePointOnTop(node1,k)
	PnT1(k) = nodePointNextTo(node1,k)
	k = indexNode1 + 1
	if(k.gt.3)k = 1
	j = indexNode1 - 1
	if(j.le.0)j = 3
	k1(k) = nodeSegConnect(node1,j)
	PoT1(k) = nodePointOnTop(node1,j)
	PnT1(k) = nodePointNextTo(node1,j)
	k = k + 1
	if(k.gt.3)k = 1
	j = indexNode2 + 1		! this is the index next seg counterclockwise from iSeg around node 1
	if(j.gt.3)j = 1			
	k1(k) = nodeSegConnect(node2,j)
	PoT1(k) = nodePointOnTop(node2,j)
	PnT1(k) = nodePointNextTo(node2,j)

	! reset node 2 segments
	k = indexNode2
	k2(k) = iSeg
	PoT2(k) = nodePointOnTop(node2,k)
	PnT2(k) = nodePointNextTo(node2,k)
	k = indexNode2 + 1
	if(k.gt.3)k = 1
	j = indexNode2 - 1
	if(j.le.0)j = 3
	k2(k) = nodeSegConnect(node2,j)
	PoT2(k) = nodePointOnTop(node2,j)
	PnT2(k) = nodePointNextTo(node2,j)
	k = k + 1
	if(k.gt.3)k = 1
	j = indexNode1 + 1		! this is the index next seg counterclockwise from iSeg around node 1
	if(j.gt.3)j = 1			
	k2(k) = nodeSegConnect(node1,j)
	PoT2(k) = nodePointOnTop(node1,j)
	PnT2(k) = nodePointNextTo(node1,j)

	do k = 1,3
		! node 1
		kk = k1(k)		! this is the seg number index around node1
		nodeSegConnect(node1,k) = kk		! this is the new segment connection
		pt = PoT1(k)
		nodePointOnTop(node1,k) = pt
		nodePointNextTo(node1,k) = PnT1(k)
		pointX(kk,pt) = nodeX(node1)
		pointY(kk,pt) = nodeY(node1)
		! node 2
		kk = k2(k)		! this is the seg number
		nodeSegConnect(node2,k) = kk
		pt = PoT2(k)
		nodePointOnTop(node2,k) = pt
		nodePointNextTo(node2,k) = PnT2(k)
		pointX(kk,pt) = nodeX(node2)
		pointY(kk,pt) = nodeY(node2)
		end do

			

	! Set the segNodes for the 2 segments that have swapped end nodes
	! For each of the 2 nodes, the segment we want is the one that is anticlockwise from the segment that joins the 2 nodes (iSeg value)
	! At this point in the code, indexNode1 is the index of iSeg in node1 and indexNode2 is the index of iSeg in node2
	! 
	! do node1
	! The nodeSegConnect are now in the NEW arrangement, so we want to go clockwise to find the correct segment
	k = indexNode1 - 1	! this should be index of the segment we want
	if(k.le.0)k = 3
	iSegT = nodeSegConnect(Node1,k)	! this should be the segment we want
	! we don't know whether the node we want to switch is the first or second in the array SegNodes
	if(segNodes(iSegT,1).eq.node2)segNodes(isegT,1) = node1
	if(segNodes(iSegT,2).eq.node2)segNodes(isegT,2) = node1

	! do node 2
	k = indexNode2 - 1	! this should be index of the segment we want
	if(k.le.0)k = 3
	iSegT = nodeSegConnect(Node2,k)	! this should be the segment we want
	! we don't know whether the node we want to switch is the first or second in the array SegNodes
	if(segNodes(iSegT,1).eq.node1)segNodes(isegT,1) = node2
	if(segNodes(iSegT,2).eq.node1)segNodes(isegT,2) = node2

	! At this point, the segNodes are properly set
	! We will use this to set the correct nodeNodeConnect
	do i = 1,3
		k = nodeSegConnect(node1,i)		! this is a segment number
		if(segNodes(k,1).eq.node1)then
			nodeNodeConnect(node1,i) = segNodes(k,2)
			else
			nodeNodeConnect(node1,i) = segNodes(k,1)
			endif
		k = nodeSegConnect(Node2,i)		! this is a segment number
		if(segNodes(k,1).eq.node2)then
			nodeNodeConnect(node2,i) = segNodes(k,2)
			else
			nodeNodeConnect(node2,i) = segNodes(k,1)
			endif
		end do



	segCaptureExclusion(iSeg) = .TRUE.	! This will ensure that a node pair is only captured onec.

	if(debug.eq.1)then


		write(12,*)'  -----After resetting ----' 
		write(12,*)'node1 ',node1
		write(12,*)'      Xl index     ',(nodeCrystalIndex(node1,i),i=1,3)
		write(12,*)'      MIFID        ',(nodeMIFID(node1,i),i=1,3)
		write(12,*)' node seg connect  ',(nodeSegConnect(node1,i),i=1,3)
		write(12,*)' node node connect ',(nodeNodeConnect(node1,i),i=1,3)
		write(12,*)'node2 ',node2
		write(12,*)'      Xl index     ',(nodeCrystalIndex(node2,i),i=1,3)
		write(12,*)'      MIFID        ',(nodeMIFID(node2,i),i=1,3)
		write(12,*)' node seg connect  ',(nodeSegConnect(node2,i),i=1,3)
		write(12,*)' node node connect ',(nodeNodeConnect(node2,i),i=1,3)
		write(12,*)'        iSeg ',iSeg
		write(12,*)' segXlID     ',segCrystalID(iSeg,1),segCrystalID(iSeg,2)
		write(12,*)' segMIFID    ',segMIFID(iSeg,1),segMIFID(iSeg,2)
		write(12,*)'     '
		write(12,*)' ------ NODE information------'

		do kkk = 1,2
			iNode = nodeT(kkk)
			!iNode = nodesWith3PhasesIndex(i)
			write(12,102)iNode,nodeX(iNode),nodeY(iNode),numNodeReactionPhases(iNode)
102			format(I6,2F12.3,I5,'         Node number,  X, Y    numNodePhases')
			MIFID = 1		! this is the grain boundary
			write(12,103)MIFID,PhName(MIFID),GBmoles,GBmoles,		&
					numPhCo(MIFID),(phCoName(MIFID,j),nodeComp(iNode,j),j=1,numPhCo(MIFID))
			!do k = 1,3
			do kk = 1,numNodeReactionPhases(iNode)
				k = nodeReactionPhases(iNode,kk)
				MIFID = nodeMIFID(iNode,k)
				write(12,103)MIFID,PhName(MIFID),nodePhaseMolesDelta(iNode,k),nodePhaseMoles(iNode,k),	&
					numPhCo(MIFID),(phCoName(MIFID,j),nodePhaseComp(iNode,k,j),j=1,numPhCo(MIFID))
103				format(T8,I4,2x,A8,2E15.5,I8,20(5x,A8,F12.5))
				end do
			write(12,*)'   '
			write(12,228)iNode,(nodeMIFID(iNode,j),j=1,numCrystalsAtNode(iNode)),			&
					(nodeReactionPhases(iNode,j),j=1,numNodeReactionPhases(iNode)),nodeSegTheSamePhases(iNode)
228			format(20I8)
			write(12,*)iNode,numMDFRxnsNodes(iNode)
			write(12,*)' ---------- Segments around this node -----------'
			write(12,*)' nodeSegConnect,segNodes(start),segNodes(end),nodePointOnTop,nodePointNextTo'
			do i = 1,3
				iSegT = nodeSegConnect(iNode,i)
				write(12,227)nodeSegConnect(iNode,i),segNodes(iSegT,1),segNodes(iSegT,2),nodePointOnTop(iNode,i),nodePointNextTo(iNode,i)
				end do
			write(12,*)' -------------------------------'
			write(12,*)

			end do

		write(12,*)' ---------- Segments around the nodes -----------'

		do kkk=1,2
			iNode=NodeT(kkk)
			write(12,*)'----------------------------'
			write(12,*)'Node = ',iNode
			do i = 1,3 	!loop on all segments around a node
				iSegT = nodeSegConnect(iNode,i)
					
				write(12,*)iSegT,numMDFRxnsSegs(iSegT)

				write(12,107)iSegT,numPointStart(iSegT),numPointEnd(iSegT),segNodes(iSegT,1),segNodes(iSegT,2),numSegReactionPhases(iSegT)
107				format(6I6,'     iSeg, numPointStart,  numPointEnd, segNodes(start), segNodes(end)   numSegReactionPhases ')
				do iPointX = numPointStart(iSegT),numPointEnd(iSegT)
					write(12,106)iSegT,iPointX,pointX(iSegT,iPointX),pointY(iSegT,iPointX)
106					format(2I8,2F12.3,'      Segment, point, X, Y ')
					MIFID = 1		! this is the grain boundary
					write(12,104)MIFID,PhName(MIFID),GBmoles,GBmoles,		&
						numPhCo(MIFID),(phCoName(MIFID,j),pointComp(iSegT,iPointX,j),j=1,numPhCo(MIFID))
104					format(T8,I4,2x,A8,2E15.5,I8,20(5x,A8,F12.5))
					do k = 1,numSegReactionPhases(iSegT)
						MIFID = segMIFID(iSegT,k)
						write(12,103)MIFID,PhName(MIFID),pointPhaseMolesDelta(iSegT,iPointX,k),		&
							pointPhaseMoles(iSegT,iPointX,k),						&
							numPhCo(MIFID),(phCoName(MIFID,j),pointPhaseComp(iSegT,iPointX,k,j),j=1,numPhCo(MIFID))
						end do
					end do
				end do
			end do

		write(12,*)'-----------------------------'
		write(12,*)'Grain Boundary composition'
		write(12,111)(elName(j),j=1,numEl)
111		format('   Seg    Point    X       Y           ',20(A3,9x))
		!write(12,110)numSegs
110		format(I5,'    Number of segments')

		write(12,109)iSeg,numPointStart(iSeg),numPointEnd(iSeg)
109		format(4x,3I5,'     iSeg,numPointStart, numPointEnd')
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
			write(12,108)iSeg,iPointX,pointX(iSeg,iPointX),pointY(iSeg,iPointX),		&
				(pointComp(iSeg,iPointX,j),j=1,numEl)
108			format(2I8,2F12.3,20F12.5)
			end do

		write(12,*)'+++++++++SegLengths+++++++++++++++++++++++++++++'
		write(12,*)iSeg,segLength(iSeg),segGettingShorter(iSeg),segCaptureExclusion(iSeg)



		endif		! end debug = 1
	
	endif		! end  if(segGettingShorter(iSeg)) = true

	endif		! end  if(segLength(iSeg).lt.dXgrid)then	! dXgrid is the nominal spacing between points along a segment


	end do		! end loop on all segments
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
			
	Subroutine PlotNode(NodePlot,iNode,iSeg1,iPoint1,P1,iSeg2,iPoint2,P2,iSeg3,iPoint3,P3)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: NodePlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	integer*4 iNode,iSeg1,iSeg2,iSeg3,iPoint1,iPoint2,iPoint3,iup
	real*8 P1(2),P2(2),P3(2),rangeX,rangeY
!	real*8 minX,minY,maxX,maxY
	real*4 x,y,dx,dy
!	common /NodePlotCanvas/NodePlot
	character*128 NodePlotTitle
	Character*5 NodeNumber


	maxX = -1.
	maxY = -1.
	minX = 1.d10
	minY = 1.d10
	if(nodeX(iNode).gt.maxX)maxX = nodeX(iNode)
	if(P1(1).ge.maxX)maxX = P1(1)
	if(P2(1).ge.maxX)maxX = P2(1)
	if(P3(1).ge.maxX)maxX = P3(1)
	if(nodeY(iNode).ge.maxY)maxY = nodeY(iNode)
	if(P1(2).ge.maxY)maxY = P1(2)
	if(P2(2).ge.maxY)maxY = P2(2)
	if(P3(2).ge.maxY)maxY = P3(2)

	if(nodeX(iNode).le.minX)minX = nodeX(iNode)
	if(P1(1).le.minX)minX = P1(1)
	if(P2(1).le.minX)minX = P2(1)
	if(P3(1).le.minX)minX = P3(1)
	if(nodeY(iNode).le.minY)minY = nodeY(iNode)
	if(P1(2).le.minY)minY = P1(2)
	if(P2(2).le.minY)minY = P2(2)
	if(P3(2).le.minY)minY = P3(2)
	
	write(12,*)' ---------------- '
	write(12,*)' Plot coordinates '
	write(12,*)' Min-X,Y  ',minX,minY
	write(12,*)' Max-X,Y  ',maxX,maxY
	rangeX = maxX - minX
	rangeY = maxY - minY
	write(12,*)'Range     ',rangeX,rangeY		
	CurrentColor = 1		! black is the default
	xmin = minX
	ymin = minY
	if(rangeY.gt.rangeX)then
		xmax = minX + rangeY
		ymax = maxY
		else
		xmax = maxX
		ymax = minY + rangeX
		endif
! 	xmin = minX
! 	ymin = minY
! 	xmax = maxX
! 	ymax = maxY

	xlen = 20		! values are in cm
	ylen = 20
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	NodePlot%width  = 1.2*(xmax-xmin)*xconv			! *xScale scales the X dimension
	NodePlot%height = 1.2*(ymax-ymin)*yconv
	xor = 20.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = NodePlot%height - 30

	NodePlotTitle = 'Node number = '
	write(NodeNumber,'(I5)')iNode
	NodePlotTitle = trim(NodePlotTitle)//NodeNumber
	NodePlot%title = NodePlotTitle
	CALL AWE_createCanvas(NodePlot)	

	pen%penStyle = CanvasPenStyle_SolidLine
	iup = 0
	x = P1(1)
	y = P1(2)
	call plot(NodePlot,x,y,iup)
	iup = 1
	x = P2(1)
	y = P2(2)
	call plot(NodePlot,x,y,iup)
	x = P3(1)
	y = P3(2)
	call plot(NodePlot,x,y,iup)
	x = P1(1)
	y = P1(2)
	call plot(NodePlot,x,y,iup)
	
	!plot the node
	x = nodeX(iNode)
	y = nodeY(iNode)
	dx = .005*(xmax-xmin)
	dy = dx
	pen%penColor = AWE_black
	brush%brushColor = AWE_black
	pen%penStyle = CanvasPenStyle_SolidLine
	Call PlotCenteredEllipseOnScreen(NodePlot,X,dX,Y,dY,brush,pen)

	pause ' See how this works before doing calculations'
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine PlotNewNode(NodePlot,Nnew)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: NodePlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
!	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	real*8 Nnew(2)
!	real*8 minX,minY,maxX,maxY
	real*4 x,y,dx,dy
!	common /NodePlotCanvas/NodePlot

	
	!plot the new node
	x = Nnew(1)
	y = Nnew(2)
	dx = .005*(xmax-xmin)
	dy = dx
	pen%penColor = AWE_red
	brush%brushColor = AWE_red
	pen%penStyle = CanvasPenStyle_SolidLine
	Call PlotCenteredEllipseOnScreen(NodePlot,X,dX,Y,dY,brush,pen)

	return
	end

	
			
			
			
			
			
			
			
			