!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine SegPlotRoutines(SegPlot1,SegPlot2,SegPlot3)
!	Routines to plot compositions along segments defined by a set of nodes
!	Nodes must be contiguous or the plots won't make much sense
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: SegPlot1
	TYPE(AWE_Canvas) :: SegPlot2
	TYPE(AWE_Canvas) :: SegPlot3
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"		

!	integer*4 iPlot,iNode,iSegs(20)
	integer*4 iPlot,iNode,reusePlot,idefault,idid
	
!	real*4 yMinMaxSi(2),yMinMaxAl(2),yMinMaxFM(2)
!	common /yMinMax/yMinMaxSi,yMinMaxAl,yMinMaxFM	

1	continue
	write(*,*)' 0 = Return'
	write(*,*)' 1 = Choose nodes to plot'
	write(*,*)' 2 = Setup SegPlot'
	write(*,*)' 3 = Plot segments'
	write(*,*)' 4 = Do all 3 (Choose nodes, setup plot, Plot segments)'
	write(*,*)'20 = Open model input file'
	write(*,*)'21 = Open model output file'
	write(*,*)'22 = Plot compositions and generate AI file'

	read(*,*)iPlot
	select case(iplot)
	case(0)
		return
	case(1)
		call FSS_Alert('ALERT','Open model input file')
		call GibbsBegin()
		call FSS_Alert('ALERT','Open model output file (pick a single cycle -- not the base')
		call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the seg and node xy coordinates
		numNodesToPlot = 0
		write(*,*)' Choose nodes to plot -- they must be adjacent!!'
		write(*,*)' end list with 0'
10		continue
		read(*,*)iNode
		if(iNode.eq.0)go to 1
		numNodesToPlot = numNodesToPlot + 1
		nodeToPlot(numNodesToPlot) = iNode
		go to 10			
		
	case(2)
		call SetupSegPlot(SegPlot1,SegPlot2,SegPlot3)

	case(3)
		Call PlotSegments(SegPlot1,SegPlot2,SegPlot3)

	case(4)
		!call FSS_Alert('ALERT','Open model input file')
		call GibbsBegin()
		call FSS_Alert('ALERT','Open model output file (pick a single cycle -- not the base')
		call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the seg and node xy coordinates
		write(*,*)' Use default nodes = 0 or input new list of nodes = 1'
		read(*,*)idefault
		if(idefault.eq.0)go to 21
		numNodesToPlot = 0
		write(*,*)' Choose nodes to plot -- they must be adjacent!!'
		write(*,*)' end list with 0'
20		continue
		read(*,*)iNode
		if(iNode.eq.0)go to 21
		numNodesToPlot = numNodesToPlot + 1
		nodeToPlot(numNodesToPlot) = iNode
		go to 20			
21		continue

		write(*,*)' Use the same plot = 0; New plot = 1'
		read(*,*)reusePlot
		if(reusePlot.eq.1) then
			call SetupSegPlot(SegPlot1,SegPlot2,SegPlot3)
			endif	
		Call PlotSegments(SegPlot1,SegPlot2,SegPlot3)
	
	case(20)
! 		call FSS_Alert('ALERT','Open model input file')
		call GibbsBegin()
	case(21)
		call FSS_Alert('ALERT','Open model output file')
		call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the seg and node xy coordinates
	case(22)
		write(*,*)' I hope you opened the base file (option 20) and the model file (option 22) already'
		write(*,*)' 0 = no I did not -- please send me back; 1 = yes I did -- get on with it'
		read(*,*)idid
		if(idid.eq.0)go to 1
		call SegPlotandAIoutput()
	case default
		call FSS_Alert('ALERT',' You chose poorly')
	end select
	go to 1
	end


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine SetupSegPlot(SegPlot1,SegPlot2,SegPlot3)
!	Routine to plot compositions along segments defined by a set of nodes
!	Nodes must be contiguous or the plots won't make much sense
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: SegPlot1
	TYPE(AWE_Canvas) :: SegPlot2
	TYPE(AWE_Canvas) :: SegPlot3
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"		
	integer*4 i,j,k,iok,iNode1,iSeg,i1,i2,istart,numPointsT
	real*4 x,y,dist,xplot(maxNodes)
	!yMinMaxSi(2),yMinMaxAl(2),yMinMaxFM(2)
	integer textSize
	character*5 text

	CurrentColor = 1		! black is the default

!	write(*,*)sideSize,numNodesPerRow,numRows

!	First sort through nodeToPlot array and determine the segments to plot
!		and whether they plot forwards or backwards
	write(*,*)' Nodes to plot '
	write(*,*)' Node1   Node2   SegtoPlot ForB  Start   End'
	numSegsToPlot = 0
	do i = 1,numNodesToPlot-1
		do j = 1,numSegs
			if(segNodes(j,1).eq.nodeToPlot(i).and.segNodes(j,2).eq.nodeToPlot(i+1))then
				numSegsToPlot = numSegsToPlot + 1
				segToPlot(numSegsToPlot) = j
				segToPlotBackwards(numSegsToPlot) = 0
				go to 10
				elseif (segNodes(j,2).eq.nodeToPlot(i).and.segNodes(j,1).eq.nodeToPlot(i+1))then
				numSegsToPlot = numSegsToPlot + 1
				segToPlot(numSegsToPlot) = j
				segToPlotBackwards(numSegsToPlot) = 1
				go to 10
				endif
			end do
			write(*,*)' Did not find seg to plot between Nodes = ',nodeToPlot(i),nodeToPlot(i+1)
			pause 'pause'
10		continue
!		write(*,*)' Node1, Node2, segToPlot ',nodeToPlot(i),nodeToPlot(i+1),segToPlot(numSegsToPlot),segToPlotBackwards(numSegsToPlot)
		iSeg = segToPlot(numSegsToPlot)
!		write(*,202)nodeToPlot(i),nodeToPlot(i+1),segToPlot(numSegsToPlot),segToPlotBackwards(numSegsToPlot)
		write(*,202)nodeToPlot(i),nodeToPlot(i+1),iSeg,segToPlotBackwards(numSegsToPlot),numPointStart(iSeg),numPointEnd(iSeg)
	202	format(6I8)
		end do


	! calculate plotting positions for node labels
	xPlot(1) = - .1		! shift a bit to the left
201	format(I5)

	write(*,*)' '
	write(*,*)' Segment plotting points and distances'

	do k = 1,numSegsToPlot
		iSeg = segToPlot(k)
		numPointsT = numPointEnd(iSeg) - numPointStart(iSeg)  ! This is actually 1 less than the number of points -- we don't want to plot the same point twice
		x = 0
		write(*,*)' iSeg, numPoints ',iSeg,numPointsT,segToPlotBackwards(k)
		do i = 1,numPointsT
			if(segToPlotBackwards(k).eq.0)then		
				! Plots frontwards - plot all but the last
				istart = numPointStart(iSeg)-1
				i1 = istart+i
				i2 = istart+i+1
				dist = sqrt((pointX(iseg,i1)-pointX(iseg,i2))**2 	&
					  + (pointY(iseg,i1)-pointY(iseg,i2))**2)
! 				dist = sqrt((pointX(iseg,numPointStart(iSeg))-pointX(iseg,numPointStart(iSeg)+1))**2 	&
! 					  + (pointY(iseg,numPointStart(iSeg))-pointY(iseg,numPointStart(iSeg)+1))**2)
				else
				! Plots backwards to front
				istart = numPointEnd(iSeg)+1
				i1 = istart-i
				i2 = istart-i-1
				dist = sqrt((pointX(iseg,i1)-pointX(iseg,i2))**2 	&
					  + (pointY(iseg,i1)-pointY(iseg,i2))**2)
				endif
			write(*,*)i,i1,i2,dist
			x = x + dist	
			end do
		xPlot(k+1) = xPlot(k) + x
! 		x = dist*float(abs((numPointend(iSeg) - numPointStart(iSeg))))
		end do		
 	xmax = xPlot(numSegsToPlot+1)

	write(*,*)' '
	write(*,*)' Nodes to plot'	
	write(*,*)' i   Node   xPlot'
	do i = 1,numNodesToPlot
		write(*,*)i,nodeToPlot(i),xPlot(i)
		end do
			
	write(*,*)' '
	write(*,*)' Xmax = (set this in plot setup)',xmax
	write(*,*)' '
!	numP = 0
!	do i = 1,numSegsToPlot
		!numP is the number of points to plot
!		numP = numP + numPoints(segToPlot(i)) - 1
!		numP = numP + (numPointEnd(segToPlot(i)) - numPointStart(segToPlot(i)))
!		end do
	PlTitle = ' Seg_Traverse'
	xmin = 0.
!	xmax = numP/10.
!	xmax = 1.
	xlen = 25.
	ylen = 20.
	nXdec = 1.
	nXstep = 20.
	xLab = 'Dist'
	yLab = 'Comp'
!	ylen = 20.
!	GridPlot%width = 1.25*(xlen * 28.34646)			! *xScale scales the X dimension
!	GridPlot%height = 1.25*(ylen * 28.34646)
	! Si plot
	ymin = 0.5
	ymax = 0.7
	nYstep = 10.
	nYdec = 2.
	SegPlot1%width =  1.25*(xLen * 28.34646)			! *xScale scales the X dimension
	SegPlot1%height = 1.25*(yLen * 28.34646)
	xor = 70.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = SegPlot1%height - 70.
	PlTitle = 'Si'
	call Setplot(iOK)
	yMinMaxSi(1) = ymin
	yMinMaxSi(2) = ymax
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	CALL AWE_createCanvas(SegPlot1)	
      	call axis(SegPlot1)			! draws the axis in AWE
	y = ymax - 0.008*(ymax-ymin)	! move a little down from the top of the plot
	textsize = 12
	do i = 1,numNodesToPlot
		iNode1 = nodeToPlot(i)
		write(text,201)iNode1
		text = adjustL(text)
		x = xPlot(i)
!		write(*,*)'x,y,iNode1,text ',x,y,iNode1,text
		CALL TextOnPlot(SegPlot1,x,y,text,textSize)
!		CALL TextOnPlot(SegPlot2,x,y,text,textSize)
!		CALL TextOnPlot(SegPlot3,x,y,text,textSize)
		end do



	! Al-H plot
	ymin = 0.1
	ymax = 0.2
	nYstep = 10.
	nYdec = 2.
	SegPlot2%title = ''
	SegPlot2%maxwindowwidth = -1
	SegPlot2%maxwindowheight = -1
	SegPlot2%width =  SegPlot1%width			! *xScale scales the X dimension
	SegPlot2%height = SegPlot1%height
	!SegPlot2%width =  int(1.25*(xLen * 28.34646))			! *xScale scales the X dimension
	!SegPlot2%height = int(1.25*(yLen * 28.34646))
	PlTitle = 'Al_H'
	call Setplot(iOK)
	yMinMaxAl(1) = ymin
	yMinMaxAl(2) = ymax
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	CALL AWE_createCanvas(SegPlot2)	
      	call axis(SegPlot2)			! draws the axis in AWE
	y = ymax - 0.008*(ymax-ymin)	! move a little down from the top of the plot
	textsize = 12
	do i = 1,numNodesToPlot
		iNode1 = nodeToPlot(i)
		write(text,201)iNode1
		text = adjustL(text)
		x = xPlot(i)
!		write(*,*)'x,y,iNode1,text ',x,y,iNode1,text
!		CALL TextOnPlot(SegPlot1,x,y,text,textSize)
		CALL TextOnPlot(SegPlot2,x,y,text,textSize)
!		CALL TextOnPlot(SegPlot3,x,y,text,textSize)
		end do
	! all others plot
	ymin = 0.
	ymax = 0.03
	nYstep = 12.
	nYdec = 3.
	SegPlot3%width =  1.25*(xLen * 28.34646)			! *xScale scales the X dimension
	SegPlot3%height = 1.25*(yLen * 28.34646)
	PlTitle = 'MgFeMnCaNaK '
	call Setplot(iOK)
	yMinMaxFM(1) = ymin
	yMinMaxFM(2) = ymax
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	CALL AWE_createCanvas(SegPlot3)	
      	call axis(SegPlot3)			! draws the axis in AWE
	y = ymax - 0.008*(ymax-ymin)	! move a little down from the top of the plot
	textsize = 12
	do i = 1,numNodesToPlot
		iNode1 = nodeToPlot(i)
		write(text,201)iNode1
		text = adjustL(text)
		x = xPlot(i)
!		write(*,*)'x,y,iNode1,text ',x,y,iNode1,text
!		CALL TextOnPlot(SegPlot1,x,y,text,textSize)
!		CALL TextOnPlot(SegPlot2,x,y,text,textSize)
		CALL TextOnPlot(SegPlot3,x,y,text,textSize)
		end do

	write(*,*)' Done drawing plots'
	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine PlotSegments(SegPlot1,SegPlot2,SegPlot3)
!	Routine to plot compositions along segments defined by a set of nodes
!	Nodes must be contiguous or the plots won't make much sense
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: SegPlot1
	TYPE(AWE_Canvas) :: SegPlot2
	TYPE(AWE_Canvas) :: SegPlot3
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	real*4 x,y,dX,dY,dist
	!,yMinMaxSi(2),yMinMaxAl(2),yMinMaxFM(2)
	integer*4 iseg,i,j,iup,k,i1,i2,numPointsT,iStart
	CurrentColor = 1		! black is the default



	pen%penStyle = CanvasPenStyle_SolidLine
! 	CurrentColor = 102		! start with the first element color (for H)
	dx = .002*(xmax-xmin)
	dy = .002*(ymax-ymin)
	Do 10 j = 1,numEl
	currentColor = 100 + j		! Note that element colors are stored in the AWEColor arrays starting at 100 = black
	select case(j)
		case(2)
			! Si
			ymin = yMinMaxSi(1)
			ymax = yMinMaxSi(2)
			call USER()
		case(1,3)
			! Al, H
			ymin = yMinMaxAl(1)
			ymax = yMinMaxAl(2)
			call USER()
		case default
			! all others
			ymin = yMinMaxFm(1)
			ymax = yMinMaxFM(2)
			call USER()
		end select
	iup = 0
!	x = 0
!	iseg = segToPlot(1)
!	dist = sqrt((pointX(iseg,2)-pointX(iseg,1))**2 + (pointY(iseg,2)-pointY(iseg,1))**2)
!	x =  -dist		! set initial x value so the plot starts at zero
	x = 0
	do k = 1,numSegsToPlot
		iseg = segToPlot(k)
		numPointsT = numPointEnd(iSeg) - numPointStart(iSeg)  ! This is actually 1 less than the number of points -- we don't want to plot the same point twice
		if(k.eq.numSegsToPlot) numPointsT = numPointsT + 1	! this will force plotting of the very last point on the last node
		do i = 1,numPointsT
			if(segToPlotBackwards(k).eq.0)then		! Plots frontwards
				iStart = numPointStart(iSeg) - 1
				i1 = iStart + i
				i2 = iStart + i + 1
				dist = sqrt((pointX(iseg,i1)-pointX(iseg,i2))**2 	&
					  + (pointY(iseg,i1)-pointY(iseg,i2))**2)
				y = pointComp(iSeg,i1,j)
				else					! plots backwards
				iStart = numPointEnd(iSeg) + 1
				i1 = iStart - i
				i2 = iStart - i -1
				dist = sqrt((pointX(iseg,i1)-pointX(iseg,i2))**2 	&
					  + (pointY(iseg,i1)-pointY(iseg,i2))**2)
				y = pointComp(iSeg,i1,j)
				endif			
			select case(j)
				case(2)
					call plot(SegPlot1,x,y,iup)
				case(1,3)
					call plot(SegPlot2,x,y,iup)
				case default
					call plot(SegPlot3,x,y,iup)
				end select
			!Call PlotCenteredEllipseOnScreen(SegPlot,X,dX,Y,dY,brush,pen)
			x = x + dist
			iup = 1
		!	x = x + dXgrid
			end do
		end do

	
! 	CurrentColor = CurrentColor + 1
10	continue
	return
	end


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine SegPlotandAIoutput()
!	Routines to plot compositions along segments defined by a set of nodes
!	Nodes must be contiguous or the plots won't make much sense
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: SegPlot1
! 	TYPE(AWE_Canvas) :: SegPlot2
! 	TYPE(AWE_Canvas) :: SegPlot3
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"		
	integer*4 idefault,j,i,iNode,iElem,iSeg,k,numPointSt,istart,i1,i2,textsize,inode1,iup,iok
	real*4 dist,x,y,xplot(maxNodes),YmaxT,YminT
	character*5 text

	write(*,*)' Use default nodes = 0 or input new list of nodes = 1'
	read(*,*)idefault
	if(idefault.eq.0)go to 21
	numNodesToPlot = 0
	write(*,*)' Choose nodes to plot -- they must be adjacent!!'
	write(*,*)' end list with 0'
20		continue
	read(*,*)iNode
	if(iNode.eq.0)go to 21
	numNodesToPlot = numNodesToPlot + 1
	nodeToPlot(numNodesToPlot) = iNode
	go to 20			
21		continue



!	First sort through nodeToPlot array and determine the segments to plot
!		and whether they plot forwards or backwards
	write(*,*)' Nodes to plot '
	write(*,*)' Node1   Node2   SegtoPlot ForB  Start   End'
	numSegsToPlot = 0
	do i = 1,numNodesToPlot-1
		do j = 1,numSegs
			if(segNodes(j,1).eq.nodeToPlot(i).and.segNodes(j,2).eq.nodeToPlot(i+1))then
				numSegsToPlot = numSegsToPlot + 1
				segToPlot(numSegsToPlot) = j
				segToPlotBackwards(numSegsToPlot) = 0
				go to 10
				elseif (segNodes(j,2).eq.nodeToPlot(i).and.segNodes(j,1).eq.nodeToPlot(i+1))then
				numSegsToPlot = numSegsToPlot + 1
				segToPlot(numSegsToPlot) = j
				segToPlotBackwards(numSegsToPlot) = 1
				go to 10
				endif
			end do
			write(*,*)' Did not find seg to plot between Nodes = ',nodeToPlot(i),nodeToPlot(i+1)
			pause 'pause'
10		continue
!		write(*,*)' Node1, Node2, segToPlot ',nodeToPlot(i),nodeToPlot(i+1),segToPlot(numSegsToPlot),segToPlotBackwards(numSegsToPlot)
		iSeg = segToPlot(numSegsToPlot)
!		write(*,202)nodeToPlot(i),nodeToPlot(i+1),segToPlot(numSegsToPlot),segToPlotBackwards(numSegsToPlot)
		write(*,202)nodeToPlot(i),nodeToPlot(i+1),iSeg,segToPlotBackwards(numSegsToPlot),numPointStart(iSeg),numPointEnd(iSeg)
	202	format(6I8)
		end do


	! calculate plotting positions for node labels
	xPlot(1) = - .1		! shift a bit to the left
201	format(I5)

	write(*,*)' '
	write(*,*)' Segment plotting points and distances'

	do k = 1,numSegsToPlot
		iSeg = segToPlot(k)
		numPointsT = numPointEnd(iSeg) - numPointStart(iSeg)  ! This is actually 1 less than the number of points -- we don't want to plot the same point twice
		x = 0
		write(*,*)' iSeg, numPoints ',iSeg,numPointsT,segToPlotBackwards(k)
		do i = 1,numPointsT
			if(segToPlotBackwards(k).eq.0)then		
				! Plots frontwards - plot all but the last
				istart = numPointStart(iSeg)-1
				i1 = istart+i
				i2 = istart+i+1
				dist = sqrt((pointX(iseg,i1)-pointX(iseg,i2))**2 	&
					  + (pointY(iseg,i1)-pointY(iseg,i2))**2)
! 				dist = sqrt((pointX(iseg,numPointStart(iSeg))-pointX(iseg,numPointStart(iSeg)+1))**2 	&
! 					  + (pointY(iseg,numPointStart(iSeg))-pointY(iseg,numPointStart(iSeg)+1))**2)
				else
				! Plots backwards to front
				istart = numPointEnd(iSeg)+1
				i1 = istart-i
				i2 = istart-i-1
				dist = sqrt((pointX(iseg,i1)-pointX(iseg,i2))**2 	&
					  + (pointY(iseg,i1)-pointY(iseg,i2))**2)
				endif
			write(*,*)i,i1,i2,dist,(pointComp(iSeg,i1,j),j=1,numEl)
			x = x + dist	
			end do
		xPlot(k+1) = xPlot(k) + x
! 		x = dist*float(abs((numPointend(iSeg) - numPointStart(iSeg))))
		end do		
 	xmax = xPlot(numSegsToPlot+1)

	write(*,*)' '
	write(*,*)' Nodes to plot'	
	write(*,*)' i   Node   xPlot'
	do i = 1,numNodesToPlot
		write(*,*)i,nodeToPlot(i),xPlot(i)
		end do
			
	write(*,*)' '
	write(*,*)' Xmax = (set this in plot setup)',xmax
	write(*,*)' '

!	Default plotting values
	PlTitle = 'SegPlot'
	xmin = 0.
	xlen = 25.
	ylen = 20.
	nXdec = 1.
	nXstep = 20.
	xLab = 'Dist'
	yLab = 'Comp'
!	ylen = 20.
!	GridPlot%width = 1.25*(xlen * 28.34646)			! *xScale scales the X dimension
!	GridPlot%height = 1.25*(ylen * 28.34646)
	! Si plot
	ymin = 0.
	ymax = 1.
	nYstep = 10.
	nYdec = 3.


!	Loop here for doing the next element
100	continue
	do i = 1,numEl
		write(*,*)i,ElName(i)
		end do
	write(*,*)' Pick element to plot - 0- to exit'
	read(*,*)iElem
	if(iElem.eq.0)return

	! find max and min in dataset
	YmaxT = -1.e6
	YminT = 1.e6
	j = iElem
! 	write(*,*)'Ymin,Ymax ',yminT,ymaxT
! 	write(*,*)' '
! 	write(*,*)' '
	do k = 1,numSegsToPlot
		iseg = segToPlot(k)
		do i = numPointStart(iSeg),numPointEnd(iSeg)
			y = pointComp(iSeg,i,j)
			if(y.gt.YmaxT)YmaxT = y
			if(y.lt.YminT)YminT = y
! 			write(*,*)'Y,Ymin,Ymax ',y,yminT,ymaxT
			end do
		end do
	write(*,*)' '
	write(*,*)' Ymin = ',YminT
	write(*,*)' Ymax = ',YmaxT		
	write(*,*)' '
! 	call PickAWEColor(myColor)
	CurrentColor = 1		!black for the axes
	SegPlot1%width =  1.25*(xLen * 28.34646)			! *xScale scales the X dimension
	SegPlot1%height = 1.25*(yLen * 28.34646)
	xor = 70.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = SegPlot1%height - 70.
	PlTitle = 'SegPlot'
	call Setplot(iOK)
!    	Call Setplot(iok)
	close(PSUNIT)         			! close old PostScript scratch file
	call psopcl(4)    			! open a new scratch file for the Illustrator output

      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	CALL AWE_createCanvas(SegPlot1)	
      	call axis(SegPlot1)			! draws the axis in AWE
	call paxis()

	y = ymax - 0.008*(ymax-ymin)	! move a little down from the top of the plot
	textsize = 12
	do i = 1,numNodesToPlot
		iNode1 = nodeToPlot(i)
		write(text,201)iNode1
		text = adjustL(text)
		x = xPlot(i)
!		write(*,*)'x,y,iNode1,text ',x,y,iNode1,text
		CALL TextOnPlot(SegPlot1,x,y,text,textSize)
		call pslab(text,x,y,14,0)
		end do


!	Here is where we plot the data
	pen%penStyle = CanvasPenStyle_SolidLine
! 	CurrentColor = myColor
	CurrentColor = iElem + 1	! colors are in the GB model input file for each element
	j = iElem	! element to plot
	iup = 0
	x = 0
	do k = 1,numSegsToPlot
		iseg = segToPlot(k)
		numPointsT = numPointEnd(iSeg) - numPointStart(iSeg)  ! This is actually 1 less than the number of points -- we don't want to plot the same point twice
		if(k.eq.numSegsToPlot) numPointsT = numPointsT + 1	! this will force plotting of the very last point on the last node
		do i = 1,numPointsT
			if(segToPlotBackwards(k).eq.0)then		! Plots frontwards
				iStart = numPointStart(iSeg) - 1
				i1 = iStart + i
				i2 = iStart + i + 1
				dist = sqrt((pointX(iseg,i1)-pointX(iseg,i2))**2 	&
					  + (pointY(iseg,i1)-pointY(iseg,i2))**2)
				y = pointComp(iSeg,i1,j)
				else					! plots backwards
				iStart = numPointEnd(iSeg) + 1
				i1 = iStart - i
				i2 = iStart - i -1
				dist = sqrt((pointX(iseg,i1)-pointX(iseg,i2))**2 	&
					  + (pointY(iseg,i1)-pointY(iseg,i2))**2)
				y = pointComp(iSeg,i1,j)
				endif			

			call plot(SegPlot1,x,y,iup)
			call pplot(x,y,iup)
			x = x + dist
			iup = 1
			end do
		end do
	call ppenup()	

	write(*,*)' Save plot as an Illustrator file?'
	write(*,*)' 0 = no'
	write(*,*)' 1 = yes'
	read(*,*)i
	if(i.eq.1)then
		call psopcl(5)    				! save old PostScript scratch file
		endif
	CurrentColor = 1		!black for the axes
	go to 100

	end


