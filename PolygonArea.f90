	Subroutine PolygonArea()

	! Area of an irregular polygon
	! The area is then given by the formula
	! ((x1*y2 − y1*x2) + (x2*y3 − y2*x3).....+(xn*y1 − yn*x1))/2
	!  
	! Where xn is the x coordinate of vertex n, 
	! yn is the y coordinate of the nth vertex etc.


	! For each crystal:
	! First sort the points that connect the nodes so they are sequential going around the crystal
	! i.e. Make an array that contains all of the polygon points in a clockwise loop
	! numNodes = numSegs and they are stored in a clockwise loop around the polygon
	! However, the points in the segments may be numbered high to low or low to high.

	implicit none
	include "Diffuse_GB_Hex.inc"
	integer*4 iXl,numPts,i,iSeg,iNode,iStart,k
	real*8 Area(maxCrystals),polyArrayX(200),polyArrayY(200)


	do iXl = 1,numCrystals
		numPts = 0	! counts the number of points in the polyArray
		do i = 1,numCrystalNodes(iXl)	! this is also the number of segments around the polygon
			iNode = CrystalNodes(iXl,i)
			iSeg = CrystalSegs(iXl,i)
			! Which point starts the segment from the node?
			iStart = nodePointOnTop(iNode,iSeg)
			if(iStart.eq.numPointStart(iSeg))then
				do k = numPointStart(iSeg),numPointEnd(iSeg)
					numPts = numPts + 1
					polyArrayX(numPts) = pointX(iSeg,k)
					polyArrayY(numPts) = pointY(iSeg,k)
					end do
				else	
				do k = numPointEnd(iSeg),numPointStart(iSeg),-1
					numPts = numPts + 1
					polyArrayX(numPts) = pointX(iSeg,k)
					polyArrayY(numPts) = pointY(iSeg,k)
					end do
				endif
			end do
		! Close the loop back to the first point
		iNode = CrystalNodes(iXl,1)
		iSeg = CrystalSegs(iXl,1)
		k = nodePointOnTop(iNode,iSeg)
		polyArrayX(numPts+1) = pointX(iSeg,k)
		polyArrayY(numPts+1) = pointY(iSeg,k)
		
		! now we should have the polyArray with all of the points around crystal iXl			
		! Calculate the areal of the polygon
		Area(iXl) = 0.0d0
		do i = 1,numPts
			Area(iXl) = Area(iXl) + polyArrayX(i)*polyArrayY(i+1) - polyArrayY(i)*polyArrayX(i+1)
			end do		
		Area(iXl) = Area(iXl)/2.0d0

		write(*,*)'iXl,Area ',iXl,CrystalPhaseName(iXl),Area(iXl)
		end do
		
	end