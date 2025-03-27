	Subroutine AffinityGrid(GridPlot)
! 	Routine to calculate the affinity for garnet nucleation at every point in the grid

	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen

	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	include "GB_Tangent.inc"
	INCLUDE "PlotStuff.inc"				

	integer*4 i,iPointX,iSeg,iGarnet,iZero,j
	real*8 minAff,maxAff
	real*4 pointAff(maxSegs,maxPoints)
	
	common /tempAff/pointAff


1	continue
	write(*,*)' Which phase is Garnet?'
	do i = 1,numPhMIF
		write(*,*)i,phName(i)
		end do
	read(*,*)iGarnet		! this is the number of the garnet in the MIF

	if(iGarnet.eq.0)return

	! We work on points in every segment
	! The endpoints of the segments are the nodes, so this will include everything

	! Calculate the Affinities for every point and determine range of affinity values
	maxAff = -1.0
	minAff = 1.0d5		! make very large to start
	write(12,*)'**********************************************************'
	write(12,*)'**********************************************************'
	do iSeg = 1,numSegs
		write(12,*)'**********************************************************'
		write(12,*)'iSeg = ',iSeg
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
			Call ParallelToTangent(iGarnet,0,izero)
			pointAff(iSeg,iPointX) = gDifferenceAU(iGarnet)
			write(12,81)iPointX,gDifferenceAU(iGarnet),(xPhCo(iGarnet,j),j=1,numPhCo(iGarnet))
81			format(5X,I5,T15,15F15.5)

			if(pointAff(iSeg,iPointX).gt.maxAff)then
				maxAff = pointAff(iSeg,iPointX)
				endif
			if(pointAff(iSeg,iPointX).lt.minAff)then
				minAff = pointAff(iSeg,iPointX)
				endif
			end do
		end do
	write(*,*)' Minimum = ',minAff
	write(*,*)' Maximum = ',maxAff

	pause 'Hit ENTER to exit'
	
	return
	end
	