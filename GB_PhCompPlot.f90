!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine PhCompPlotRoutines(SegPlot1,SegPlot2,SegPlot3)
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
	include "GibbsFiles.inc"
	include "Assemb.inc"

	integer*4 iNode,phaseToPlot,iok,i,iPlot,plotType

1	continue
	write(*,*)' 0 = Return'
	!write(*,*)' 1 = Choose node to plot'
	!write(*,*)' 2 = Setup PhCompPlots'
	!write(*,*)' 3 = Plot compositions'
	write(*,*)' 4 = Plot compositions - old code'
	write(*,*)'44 = Plot compositions - new code'
	write(*,*)' 5 = Plot Affinities at Garnet node '
	write(*,*)' 6 = Plot Affinities for Garnet nucleation '
	write(*,*)'20 = Open model input file'
	write(*,*)'21 = Open model output file'

	read(*,*)iPlot
	select case(iplot)
	case(0)
		return
	case(1)
		numNodesToPlot = 0
		write(*,*)' Choose node to plot -- the node must contain the phase of interest!!'
		read(*,*)iNode
		if(iNode.eq.0)go to 1
		write(*,*)' Choose phase to plot '
		do i = 1,3
		write(*,*)i,NodeMIFID(iNode,i),phName(nodeMIFID(iNode,i))
		end do
		read(*,*)phaseToPlot
		if(phaseToPlot.eq.0)go to 1
		
	case(2)
		write(*,*)' You need to know the maximum radius or moles to plot on the x-axis before you continue'
		write(*,*)' Do you want to continue? 0 = no; 1 = yes'
		read(*,*)iok
		if(iok.eq.0)go to 1
		call SetupPhCompPlot(SegPlot1,SegPlot2,SegPlot3,plotType)
	
	case(3)

		!Call PlotPhComp(SegPlot1,SegPlot2,SegPlot3,iNode,phaseToPlot,plotType)
	
	case(4)

		Call PlotPhComp2(SegPlot1,SegPlot2,SegPlot3,iNode,phaseToPlot,plotType)

	case(44)

		Call PlotPhComp3(SegPlot1,SegPlot2,SegPlot3,iNode,phaseToPlot,plotType)
	
	case(5)

		Call PlotAffinities(SegPlot1,iNode)
	
	case(6)

		Call PlotAffinity_At_Node(SegPlot1,iNode)
	
	case(20)
		call GibbsBegin()
	case(21)
		call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the node and seg XY coordinates
	case default
		call FSS_Alert('ALERT',' You chose poorly')
	end select
	go to 1
	end


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine SetupPhCompPlot(SegPlot1,SegPlot2,SegPlot3,plotType)
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
	integer*4 iok,plotType

	CurrentColor = 1		! black is the default

1	continue
	xmin = 0.
	xmax = 1.		! just a place holder
	xlen = 25.
	ylen = 20.
	nXdec = 1.
	nXstep = 10.
	if(plotType.eq.1)then
		xLab = 'Radius'
		else
		xLab = 'Moles'
		endif
	yLab = 'Comp'

	! Fe/(Fe+Mg) and almandine plot
	ymin = 0.5
	ymax = 1.0
	nYstep = 10.
	nYdec = 2.
	SegPlot1%width =  1.25*(xLen * 28.34646)			! *xScale scales the X dimension
	SegPlot1%height = 1.25*(yLen * 28.34646)
	xor = 70.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = SegPlot1%height - 70.
	PlTitle = 'Alm_FeMg'
	call Setplot(iOK)
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	CALL AWE_createCanvas(SegPlot1)	
      	call axis(SegPlot1)			! draws the axis in AWE

	! Prp, Sps, Grs plot
	ymin = 0.0
	ymax = 0.3
	nYstep = 15.
	nYdec = 2.
	SegPlot2%title = ''
	SegPlot2%maxwindowwidth = -1
	SegPlot2%maxwindowheight = -1
	SegPlot2%width =  SegPlot1%width			! *xScale scales the X dimension
	SegPlot2%height = SegPlot1%height
	!SegPlot2%width =  int(1.25*(xLen * 28.34646))			! *xScale scales the X dimension
	!SegPlot2%height = int(1.25*(yLen * 28.34646))
	PlTitle = 'PrpSpsGrs'
	call Setplot(iOK)
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	CALL AWE_createCanvas(SegPlot2)	
      	call axis(SegPlot2)			! draws the axis in AWE

	return
	! all others plot
	ymin = 0.
	ymax = 0.03
	nYstep = 12.
	nYdec = 3.
	SegPlot3%width =  1.25*(xLen * 28.34646)			! *xScale scales the X dimension
	SegPlot3%height = 1.25*(yLen * 28.34646)
	PlTitle = 'MgFeMnCaNaK '
	call Setplot(iOK)
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	CALL AWE_createCanvas(SegPlot3)	
      	call axis(SegPlot3)			! draws the axis in AWE

	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
! 	Subroutine PlotPhComp(SegPlot1,SegPlot2,SegPlot3,iNode,phaseToPlot,plotType)
! !	Routine to plot compositions along segments defined by a set of nodes
! !	Nodes must be contiguous or the plots won't make much sense
! 	USE AWE_INTERFACES
!       	implicit none
! 	TYPE(AWE_Canvas) :: SegPlot1
! 	TYPE(AWE_Canvas) :: SegPlot2
! 	TYPE(AWE_Canvas) :: SegPlot3
! 	Type(AWE_CanvasBrush) :: brush
! 	TYPE(AWE_CanvasPen) :: pen
! 	include "Assemb.inc"
! 	include "Diffuse_GB_Hex.inc"
! 	INCLUDE "PlotStuff.inc"				
! 	include "GibbsFiles.inc"
! 	real*4 x,y,dX,dY
! 	integer*4 i,j,iup,iNode,phaseToPlot,plotType,status,stepStart,stepEnd,isteps,numToPlot,MIFID,iOK
! 	real*4 phMolesToPlot(4000),phCompToPlot(4000,5),maxMoles,maxRadius
! 	character*16 stepText
! 
! 	CurrentColor = 1		! black is the default
! 
! 
! 
! 	! Open and read the model information file
! 	write(*,*)' Open the model input file'
! 	call FSS_Alert('ALERT','Open model input file')
! 	call GibbsBegin()
! 
! 
! 	
! 	! We already have a base. Just open this base file name
! 	write(*,*)' Open the model file base name'
! 	call FSS_Alert('ALERT','Open model file base name')
! 	open(73,FILE='',status='OLD',iostat=status,action='WRITE')
! 	if(status.ne.0)then
! 		call FSS_Alert('Alert','Problem opening model base file')
! 		return
! 		endif
! 	inquire(73,NAME=ModelOutputFileBase)
! 	close(73)	! we don't need this file for anything
! 
! 1	continue
! 	numNodesToPlot = 0
! 	write(*,*)' Choose node to plot -- the node must contain the phase of interest!!'
! 	read(*,*)iNode
! 	if(iNode.eq.0)return
! 	write(*,*)' Choose phase to plot '
! 	do i = 1,3
! 		write(*,*)i,NodeMIFID(iNode,i),phName(nodeMIFID(iNode,i))
! 		end do
! 	read(*,*)phaseToPlot
! 	if(phaseToPlot.eq.0)return
! 
! 	write(*,*)' Plot type:'
! 	write(*,*)' 0 = return'
! 	write(*,*)' 1 = radius vs composition'
! 	write(*,*)' 2 = moles vs composition'
! 	read(*,*)plotType
! 	if(plotType.eq.0)return
! 
! 	write(*,*)'Input starting step for making the plot (0 to begin at the beginning)'
! 	read(*,*)stepStart
! 	write(*,*)'Input ending step for making the plot '
! 	read(*,*)stepEnd
! 
! 	if(stepEnd-stepStart.ge.1999)then
! 		call FSS_Alert('ALERT',' Number of steps exceeds array allocation of 4000')
! 		return
! 		endif
! 
! !	read all the models and collect data for the plots
! 	maxMoles = 0
! 	maxRadius = 0
! 	numToPlot = 0
! 	do isteps = stepStart,stepEnd
! 		write(*,*)' Step number ',isteps
! 		write(steptext,55)isteps
! 55		format(I5)
! 		steptext = adjustL(steptext)
! 		ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
! 		write(*,*)'ModelOutputFile name ',ModelOutputFile
! 		open(40,FILE=ModelOutputFile,status = 'UNKNOWN')
! 
! 		! Now read the model information
! 		isFileOpen = 1		! file is already open
! 		call OpenModelOutputFile(isFileOpen)
! 
! 		! find the correct node and phase
! 		MIFID = NodeMIFID(iNode,phaseToPlot)
! 		do i = 1,numNodesWith3Phases
! 			if(iNode.eq.nodesWith3PhasesIndex(i))then
! 			! we have found the correct node. Now get the phase info
! 			numToPlot = numToPlot + 1
! 			phMolesToPlot(numToPlot) = nodePhaseMoles(iNode,phaseToPlot)
! 			do j = 1,numPhCo(MIFID)
! 				phCompToPlot(numToPlot,j) = nodePhaseComp(iNode,phaseToPlot,j)
! 				end do
! 			!Fe/Fe+Mg
! 			phCompToPlot(numToPlot,5) = phCompToPlot(numToPlot,2)/(phCompToPlot(numToPlot,2)+phCompToPlot(numToPlot,1))
! 			go to 10
! 			endif
! 			end do
! 		do i = 1,numNodesWith2Phases
! 			if(iNode.eq.nodesWith2PhasesIndex(i))then
! 			! we have found the correct node. Now get the phase info
! 			numToPlot = numToPlot + 1
! 			phMolesToPlot(numToPlot) = nodePhaseMoles(iNode,phaseToPlot)
! 			do j = 1,numPhCo(MIFID)
! 				phCompToPlot(numToPlot,j) = nodePhaseComp(iNode,phaseToPlot,j)
! 				end do
! 			!Fe/Fe+Mg
! 			phCompToPlot(numToPlot,5) = phCompToPlot(numToPlot,2)/(phCompToPlot(numToPlot,2)+phCompToPlot(numToPlot,1))
! 			go to 10
! 			endif
! 			end do
! 		call FSS_Alert('ALERT','Did not find the phase in the node')
! 
! 10		continue			
! 		if(phMolesToPlot(numToPlot).gt.maxMoles)maxMoles = phMolesToPlot(numToPlot)
! 		end do		! go get the next step				
! 
! 	! All steps are now in the arrays
! 	write(*,*)' Number of points to plot = ',numToPlot
! 	write(*,*)' Max mmoles in dataset    = ',maxMoles*1000.
! 	!x = (phMolesToPlot(i)*112.)		! volume in cubic centimeters
! 	!x = ((3./(4.*3.14159))*x)**0.33333	! radius in cm
! 	!x = x*1.0e4				! radius in µm
! 
! 	maxRadius = 1.0e4*(maxMoles*112.*((3./(4.*3.14159))))**0.3333333	!radius in µm
! 	write(*,*)' Max radius in dataset     = ',maxRadius
! 
! 	xmin = 0.
! 	if(plotType.eq.1)then
! 		xmax = maxRadius
! 		xLab = 'Radius'
! 		else
! 		xmax = maxMoles*1000
! 		xLab = 'milliMoles'
! 		endif
! 	xlen = 25.
! 	ylen = 20.
! 	nXdec = 1.
! 	nXstep = 10.
! 	yLab = 'Comp'
! 
! 	! Fe/(Fe+Mg) and almandine plot
! 	ymin = 0.5
! 	ymax = 1.0
! 	nYstep = 10.
! 	nYdec = 2.
! 	SegPlot1%width =  1.25*(xLen * 28.34646)			! *xScale scales the X dimension
! 	SegPlot1%height = 1.25*(yLen * 28.34646)
! 	xor = 70.				! origin for X-Y plot (xmin,ymin) in pixels
! 	yor = SegPlot1%height - 70.
! 	PlTitle = 'Alm_FeMg'
! 	call Setplot(iOK)
!       	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
! 	CALL AWE_createCanvas(SegPlot1)	
!       	call axis(SegPlot1)			! draws the axis in AWE
! 
! 	! Prp, Sps, Grs plot
! 	ymin = 0.0
! 	ymax = 0.3
! 	nYstep = 15.
! 	nYdec = 2.
! 	SegPlot2%title = ''
! 	SegPlot2%maxwindowwidth = -1
! 	SegPlot2%maxwindowheight = -1
! 	SegPlot2%width =  SegPlot1%width			! *xScale scales the X dimension
! 	SegPlot2%height = SegPlot1%height
! 	!SegPlot2%width =  int(1.25*(xLen * 28.34646))			! *xScale scales the X dimension
! 	!SegPlot2%height = int(1.25*(yLen * 28.34646))
! 	PlTitle = 'PrpSpsGrs'
! 	call Setplot(iOK)
!       	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
! 	CALL AWE_createCanvas(SegPlot2)	
!       	call axis(SegPlot2)			! draws the axis in AWE
! 	
! 
! 	dx = .002*(xmax-xmin)
! 	dy = .002*(ymax-ymin)
! 
! 	pen%penStyle = CanvasPenStyle_SolidLine
! 	do j = 1,numPhCo(MIFID)+1
! 		iup = 0	
! 
! 		! set the user scale for the component to plot
! 		select case(j)
! 		case(1)		!Prp
! 			currentColor = 5	! Cyan
! 			pen%penColor = AWEColorNumber(CurrentColor)
! 			brush%brushColor = AWEColorNumber(CurrentColor)
! !			pen%penStyle = CanvasPenStyle_SolidLine
! 			ymin = 0.0
! 			ymax = 0.3
! 			call USER()
! 		case(3)		!Sps
! 			currentColor = 7	! Magenta
! 			pen%penColor = AWEColorNumber(CurrentColor)
! 			brush%brushColor = AWEColorNumber(CurrentColor)
! 			ymin = 0.0
! 			ymax = 0.3
! 			call USER()
! 		case(4)		!Grs
! 			currentColor = 8	! green
! 			pen%penColor = AWEColorNumber(CurrentColor)
! 			brush%brushColor = AWEColorNumber(CurrentColor)
! 			ymin = 0.0
! 			ymax = 0.3
! 			call USER()
! 		case(2)		! Alm
! 			currentColor = 3	! red
! 			pen%penColor = AWEColorNumber(CurrentColor)
! 			brush%brushColor = AWEColorNumber(CurrentColor)
! 			ymin = 0.5
! 			ymax = 1.0
! 			call USER()
! 		case(5)		! Fe/Fe+Mg
! 			currentColor = 4	! teal
! 			pen%penColor = AWEColorNumber(CurrentColor)
! 			brush%brushColor = AWEColorNumber(CurrentColor)
! 			ymin = 0.5
! 			ymax = 1.0
! 			call USER()
! 		case default
! 		end select								
! 
! 		do i = 1,numToPlot
! 			y = phCompToPlot(i,j)
! 
! 			select case (plotType)
! 			case(1)			! plot radius on X
! 				x = (phMolesToPlot(i)*112.)		! volume in cubic centimeters
! 				x = ((3./(4.*3.14159))*x)**0.33333	! radius in cm
! 				x = x*1.0e4				! radius in µm
! 
! 			case(2)			! plot moles on X
! 				x = 1000.*phMolesToPlot(i)
! 
! 			case default
! 			end select
! 
! 			select case(j)
! 			case(1,3,4)		!Prp, Sps, Alm
! 				call plot(SegPlot2,x,y,iup)
! 				Call PlotCenteredEllipseOnScreen(SegPlot2,X,dX,Y,dY,brush,pen)
! 
! 			case(2,5)
! 				call plot(SegPlot1,x,y,iup)
! 				Call PlotCenteredEllipseOnScreen(SegPlot1,X,dX,Y,dY,brush,pen)
! 			case default
! 			end select								
! 			iup = 1
! 			end do
! 		end do
! 	currentColor = 1
! 	go to 1
! 	end
		
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine PlotPhComp2(SegPlot1,SegPlot2,SegPlot3,iNode,phaseToPlot,plotType)
!	Routine to plot compositions along segments defined by a set of nodes
!	Nodes must be contiguous or the plots won't make much sense
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: SegPlot1
	TYPE(AWE_Canvas) :: SegPlot2
	TYPE(AWE_Canvas) :: SegPlot3
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	include "GibbsFiles.inc"
	real*4 x,y,dX,dY,maxUnits
	integer*4 i,j,iup,k,kk,iNode,phaseToPlot,plotType,status,stepStart,stepEnd,isteps,numToPlot,MIFID,iOK
	real*4 phMolesToPlot(maxNodes,10000),phCompToPlot(maxNodes,10000,5),maxMoles,maxRadius
	real*4 nodeDistance(maxNodes,10000),nodeXstart(maxNodes),nodeYstart(maxNodes),nodeDistanceMax
	character*16 stepText
	integer*4 phaseToPlotMIFID,numNodesWithPhaseToPlot,NodesWithPhaseToPlot(maxNodes),reusePlot
	integer*4 totalCyclesToPlot(10000),maxCycles,minCycles
	common /plotstorage/phMolesToPlot,phCompToPlot,totalCyclesToPlot,NodesWithPhaseToPlot,maxMoles,maxRadius,nodeDistance
	
	CurrentColor = 1		! black is the default



	! Open and read the model information file
	write(*,*)' Open the model input file'
	!call FSS_Alert('ALERT','Open model input file')
	call GibbsBegin()


	
	! We already have a base. Just open this base file name
	write(*,*)' Open the model file base name'
	call FSS_Alert('ALERT','Open model file base name')
	open(73,FILE='',status='OLD',iostat=status,action='WRITE')
	if(status.ne.0)then
		call FSS_Alert('Alert','Problem opening model base file')
		return
		endif
	inquire(73,NAME=ModelOutputFileBase)
	close(73)	! we don't need this file for anything

1	continue
	write(*,*)' Which phase is Garnet?'
	do i = 1,numPhMIF
		write(*,*)i,phName(i)
		end do
	read(*,*)phaseToPlotMIFID		! this is the number of the phase in the MIF
	!phaseToPlotName = phName(phaseToPlotMIFID)
	if(phaseToPlotMIFID.eq.0)return

	numNodesWithPhaseToPlot = 0

	write(*,*)'Input starting step for making the plot (0 to begin at the beginning)'
	read(*,*)stepStart
	write(*,*)'Input ending step for making the plot '
	read(*,*)stepEnd

	if(stepEnd-stepStart.ge.9999)then
		call FSS_Alert('ALERT',' Number of steps exceeds array allocation of 10000')
		return
		endif

!	read all the models and collect data for the plots
	maxMoles = 0
	maxRadius = 0
	numToPlot = 0
	nodeDistanceMax = 0
	minCycles = 100000
	maxCycles = -10
	do isteps = stepStart,stepEnd
		numToPlot = numToPlot + 1
		write(steptext,55)isteps
55		format(I5)
		steptext = adjustL(steptext)
		ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
		if(isteps.eq.stepStart.or.isteps.eq.stepEnd)then
			write(*,*)' Step number ',isteps
			write(*,*)'ModelOutputFile name ',ModelOutputFile		! only write out first and last file names
			endif
		open(40,FILE=ModelOutputFile,status = 'UNKNOWN')

		! Now read the model information
		isFileOpen = 1		! file is already open
!		call OpenModelOutputFile_2(isFileOpen,1)		! This one reads affinities. 1 means read the node and seg XY coordinates
		call OpenModelOutputFile(isFileOpen,1)			! No affinities in this one. 1 means read the node and seg XY coordinates

		totalCyclesToPlot(numToPlot) = TotalCycles
		if(TotalCycles.gt.maxCycles)maxCycles = totalCycles
		if(TotalCycles.lt.minCycles)minCycles = totalCycles
		! find the correct phase and store in arrays
! 		do i = 1,numNodesWith3Phases
! 			iNode = nodesWith3PhasesIndex(i)
		do iNode = 1,numNodes
			if(numToPlot.eq.1)then
				nodeXstart(iNode) = nodeX(iNode)		! This could also be the crystal center
				nodeYstart(iNode) = nodeY(iNode)
				endif
			nodeDistance(iNode,numToPlot) = sqrt((nodeX(iNode) - nodeXstart(iNode))**2 + 	&
				         (nodeY(iNode) - nodeYstart(iNode))**2)
			if(nodeDistance(iNode,numToPlot).gt.nodeDistanceMax)then
				nodeDistanceMax = nodeDistance(iNode,numToPlot)
				endif
			do j = 1,5		! zero out array (in case there is no grossular)
				phCompToPlot(iNode,numToPlot,j) = 0.0			
				end do
				
			select case(numNodeReactionPhases(iNode))
			case(3)
				do k = 1,3		! 4 phases in nodesWith3Phases. #1 is always the GB.
				!nodeMIFID(iNode,k) = MIFID
					MIFID = NodeMIFID(iNode,k)

					if(MIFID.eq.phaseToPlotMIFID)then
						if(isteps.eq.1)then		! only set these the first time through
							numNodesWithPhaseToPlot = numNodesWithPhaseToPlot + 1
							nodesWithPhaseToPlot(numNodesWithPhaseToPlot) = iNode
							endif
						phMolesToPlot(iNode,numToPlot) = nodePhaseMoles(iNode,k)
						do j = 1,numPhCo(MIFID)					! Note that this assumes order is prp Alm Sps Grs
							phCompToPlot(iNode,numToPlot,j) = nodePhaseComp(iNode,k,j)
							end do
						phCompToPlot(iNode,numToPlot,5) = &		! Fe/Fe+Mg assuming this is garnet
							nodePhaseComp(iNode,k,2)/(nodePhaseComp(iNode,k,2)+nodePhaseComp(iNode,k,1))
							!phCompToPlot(numToPlot,2)/(phCompToPlot(numToPlot,2)+phCompToPlot(numToPlot,1))
						if(phMolesToPlot(iNode,numToPlot).gt.maxMoles)maxMoles = phMolesToPlot(iNode,numToPlot)
						endif
					end do
			case(2)
	! 		do i = 1,numNodesWith2Phases
	! 			iNode = nodesWith2PhasesIndex(i)
				do kk = 1,2		! 3 phases in nodesWith2Phases. #1 is always the GB.
				!nodeMIFID(iNode,k) = MIFID
					k = nodeReactionPhases(iNode,kk)
					MIFID = NodeMIFID(iNode,k)

					if(MIFID.eq.phaseToPlotMIFID)then
						if(isteps.eq.1)then		! only set these the first time through
							numNodesWithPhaseToPlot = numNodesWithPhaseToPlot + 1
							nodesWithPhaseToPlot(numNodesWithPhaseToPlot) = iNode
							endif
						phMolesToPlot(iNode,numToPlot) = nodePhaseMoles(iNode,k)
						do j = 1,numPhCo(MIFID)
							phCompToPlot(iNode,numToPlot,j) = nodePhaseComp(iNode,k,j)
							end do
						phCompToPlot(iNode,numToPlot,5) = &		! Fe/Fe+Mg
							nodePhaseComp(iNode,k,2)/(nodePhaseComp(iNode,k,2)+nodePhaseComp(iNode,k,1))
							!phCompToPlot(numToPlot,2)/(phCompToPlot(numToPlot,2)+phCompToPlot(numToPlot,1))
						if(phMolesToPlot(iNode,numToPlot).gt.maxMoles)maxMoles = phMolesToPlot(iNode,numToPlot)
						endif
					end do
			case default
			end select
			end do

		end do		! go get the next step				

	! All steps are now in the arrays
	write(*,*)' Number of points to plot = ',numToPlot
	write(*,*)' Max mmoles in dataset    = ',maxMoles*1000.
	!x = (phMolesToPlot(i)*112.)		! volume in cubic centimeters
	!x = ((3./(4.*3.14159))*x)**0.33333	! radius in cm
	!x = x*1.0e4				! radius in µm

	maxRadius = 1.0e4*(maxMoles*112.*((3./(4.*3.14159))))**0.3333333	!radius in µm
	write(*,*)' Max radius in dataset     = ',maxRadius
	write(*,*)' Max Distance in dataset     = ',nodeDistanceMax
	write(*,*)' Min Cycles in dataset     = ',minCycles
	write(*,*)' Max Cycles in dataset     = ',maxCycles

11	continue
	write(*,*)' Plot type:'
	write(*,*)' 0 = return'
	write(*,*)' 1 = radius vs composition'
	write(*,*)' 2 = moles vs composition'
	write(*,*)' 3 = Model step vs composition'
	write(*,*)' 4 = Node position vs composition'
	write(*,*)' 5 = Cycle number vs composition'
	write(*,*)' 6 = Write garnet composition to output window.'
	read(*,*)plotType
	if(plotType.eq.0)go to 1


10	continue			
	do i = 1,numNodesWithPhaseToPlot
		if(nodesWithPhaseToPlot(i).gt.0)then
			write(*,*)nodesWithPhaseToPlot(i)
			endif
		end do
	write(*,*)' Choose the the node to plot --- enter 0 to exit'
	read(*,*)iNode
	if(iNode.eq.0)go to 11
	!iNode = nodesWithPhaseToPlot(iNodeIndex)
	
	write(*,*)' Use the same plot = 0; New plot = 1'
	read(*,*)reusePlot
	if(reusePlot.eq.0) go to 20

	xmin = 0.
	select case(plotType)
		case(1)
			xmax = maxRadius
			xLab = 'Radius'
		case(2)
			xmax = maxMoles*1000
			xLab = 'milliMoles'
		case(3)
			xmax = numToPlot
			xLab = 'Model_Step'
		case(4)
			xmax = 40.
			xLab = 'Node_Position'
		case(5)
			xmax = float(maxCycles)
			xLab = 'Cycles'
		case(6)		! write out garnet compositions -- actually, only Mn is coded here
			! calculate the maximum arbitrary units for the plot
			x = (phMolesToPlot(iNode,numToPlot)*112.)		! volume in cubic centimeters
			x = ((3./(4.*3.14159))*x)**0.33333	! radius in cm
			x = x*1.0e4				! radius in µm -- these are actually arbitrary units
			maxUnits = x
			write(*,*)' Input actual maximum radius in µm (calculated from NODE positions)'
			read(*,*)maxRadius
			write(12,*)'  '
			write(12,*)'  '
			write(12,*)' Spessartine composition vs radius'
			write(12,*)' Radius   Xsps'
			do i = 1,numToPlot
				y = phCompToPlot(iNode,i,3)		! This will only list spessartine composition (for diffusion modeling)
	
				x = (phMolesToPlot(iNode,i)*112.)		! volume in cubic centimeters
				x = ((3./(4.*3.14159))*x)**0.33333	! radius in cm
				x = x*1.0e4				! radius in µm -- these are actually arbitrary units
				write(12,81)x,x*maxRadius/maxUnits,y
	81			format(3F15.5)
				end do
			write(12,*)'  '
			write(12,*)'  '
			go to 11


		case default
			call FSS_Alert('ALERT','You chose poorly')
			go to 11
		end select
	xlen = 25.
	ylen = 20.
	nXdec = 1.
	nXstep = 10.
	yLab = 'Comp'

	! Fe/(Fe+Mg) and almandine plot
	ymin = 0.5
	ymax = 1.0
	nYstep = 10.
	nYdec = 2.
	SegPlot1%width =  1.25*(xLen * 28.34646)			! *xScale scales the X dimension
	SegPlot1%height = 1.25*(yLen * 28.34646)
	xor = 70.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = SegPlot1%height - 70.
	PlTitle = 'Alm_FeMg'
	call Setplot(iOK)
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	CALL AWE_createCanvas(SegPlot1)	
      	call axis(SegPlot1)			! draws the axis in AWE

	! Prp, Sps, Grs plot
	ymin = 0.0
	ymax = 0.3
	nYstep = 15.
	nYdec = 2.
	SegPlot2%title = ''
	SegPlot2%maxwindowwidth = -1
	SegPlot2%maxwindowheight = -1
	SegPlot2%width =  SegPlot1%width			! *xScale scales the X dimension
	SegPlot2%height = SegPlot1%height
	!SegPlot2%width =  int(1.25*(xLen * 28.34646))			! *xScale scales the X dimension
	!SegPlot2%height = int(1.25*(yLen * 28.34646))
	PlTitle = 'PrpSpsGrs'
	call Setplot(iOK)
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	CALL AWE_createCanvas(SegPlot2)	
      	call axis(SegPlot2)			! draws the axis in AWE
	
20	continue
	dx = .002*(xmax-xmin)
	dy = .002*(ymax-ymin)

	pen%penStyle = CanvasPenStyle_SolidLine
!	do j = 1,numPhCo(phaseToPlotMIFID)+1
	do j = 1,5			! 4 garnet + Fe/Fe+Mg
		iup = 0	

		! set the user scale for the component to plot
		select case(j)
		case(1)		!Prp
			currentColor = 5	! Cyan
			pen%penColor = AWEColorNumber(CurrentColor)
			brush%brushColor = AWEColorNumber(CurrentColor)
!			pen%penStyle = CanvasPenStyle_SolidLine
			ymin = 0.0
			ymax = 0.3
			call USER()
		case(3)		!Sps
			currentColor = 7	! Magenta
			pen%penColor = AWEColorNumber(CurrentColor)
			brush%brushColor = AWEColorNumber(CurrentColor)
			ymin = 0.0
			ymax = 0.3
			call USER()
		case(4)		!Grs
			currentColor = 4	! green
			pen%penColor = AWEColorNumber(CurrentColor)
			brush%brushColor = AWEColorNumber(CurrentColor)
			ymin = 0.0
			ymax = 0.3
			call USER()
		case(2)		! Alm
			currentColor = 2	! red
			pen%penColor = AWEColorNumber(CurrentColor)
			brush%brushColor = AWEColorNumber(CurrentColor)
			ymin = 0.5
			ymax = 1.0
			call USER()
		case(5)		! Fe/Fe+Mg
			currentColor = 1	! Black
			pen%penColor = AWEColorNumber(CurrentColor)
			brush%brushColor = AWEColorNumber(CurrentColor)
			ymin = 0.5
			ymax = 1.0
			call USER()
		case default
		end select								

		do i = 1,numToPlot
			y = phCompToPlot(iNode,i,j)

			select case (plotType)
			case(1)			! plot radius on X
				x = (phMolesToPlot(iNode,i)*112.)		! volume in cubic centimeters
				x = ((3./(4.*3.14159))*x)**0.33333	! radius in cm
				x = x*1.0e4				! radius in µm

			case(2)			! plot moles on X
				x = 1000.*phMolesToPlot(iNode,i)

			case(3)			! plot cycle number (i.e. time) on X
				x = i
			case(4)			! node position
				x = nodeDistance(iNode,i)
			case(5)			! Cycle number
				x = totalCyclesToPlot(i)

			case default
			end select

			select case(j)
			case(1,3,4)		!Prp, Sps, Alm
				call plot(SegPlot2,x,y,iup)
				Call PlotCenteredEllipseOnScreen(SegPlot2,X,dX,Y,dY,brush,pen)

			case(2,5)
				call plot(SegPlot1,x,y,iup)
				Call PlotCenteredEllipseOnScreen(SegPlot1,X,dX,Y,dY,brush,pen)
			case default
			end select								
			iup = 1
			end do
		end do
	currentColor = 1
	go to 10
	end
		
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine PlotPhComp3(SegPlot1,SegPlot2,SegPlot3,iNode,phaseToPlot,plotType)
!	Routine to plot compositions along segments defined by a set of nodes
!	Nodes must be contiguous or the plots won't make much sense
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: SegPlot1
	TYPE(AWE_Canvas) :: SegPlot2
	TYPE(AWE_Canvas) :: SegPlot3
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	include "GibbsFiles.inc"
	real*4 x,y,dX,dY,maxUnits
	integer*4 i,j,iup,k,kk,iNode,phaseToPlot,plotType,status,stepStart,stepEnd,isteps,numToPlot,MIFID,iOK,whatToPlot
	real*4 phMolesToPlot(maxNodes,10000),phCompToPlot(maxNodes,10000,5),maxMoles,maxRadius
	real*4 nodeDistance(maxNodes,10000),nodeXstart(maxNodes),nodeYstart(maxNodes),nodeDistanceMax,xToPlot(10000)
	character*16 stepText
	integer*4 phaseToPlotMIFID,numNodesWithPhaseToPlot,NodesWithPhaseToPlot(maxNodes),reusePlot
	integer*4 totalCyclesToPlot(10000),maxCycles,minCycles
	common /plotstorage/phMolesToPlot,phCompToPlot,totalCyclesToPlot,NodesWithPhaseToPlot,maxMoles,maxRadius,nodeDistance
	
	CurrentColor = 1		! black is the default



	! Open and read the model information file
	write(*,*)' Open the model input file'
	!call FSS_Alert('ALERT','Open model input file')
	call GibbsBegin()


	
	! We already have a base. Just open this base file name
	write(*,*)' Open the model file base name'
	call FSS_Alert('ALERT','Open model file base name')
	open(73,FILE='',status='OLD',iostat=status,action='WRITE')
	if(status.ne.0)then
		call FSS_Alert('Alert','Problem opening model base file')
		return
		endif
	inquire(73,NAME=ModelOutputFileBase)
	close(73)	! we don't need this file for anything

1	continue
	write(*,*)' Which phase is Garnet?'
	do i = 1,numPhMIF
		write(*,*)i,phName(i)
		end do
	read(*,*)phaseToPlotMIFID		! this is the number of the phase in the MIF
	!phaseToPlotName = phName(phaseToPlotMIFID)
	if(phaseToPlotMIFID.eq.0)return

	numNodesWithPhaseToPlot = 0

	write(*,*)'Input starting step for making the plot (0 to begin at the beginning)'
	read(*,*)stepStart
	write(*,*)'Input ending step for making the plot '
	read(*,*)stepEnd

	if(stepEnd-stepStart.ge.9999)then
		call FSS_Alert('ALERT',' Number of steps exceeds array allocation of 10000')
		return
		endif

!	read all the models and collect data for the plots
	maxMoles = 0
	maxRadius = 0
	numToPlot = 0
	nodeDistanceMax = 0
	minCycles = 100000
	maxCycles = -10
	do isteps = stepStart,stepEnd
		numToPlot = numToPlot + 1
		write(steptext,55)isteps
55		format(I5)
		steptext = adjustL(steptext)
		ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
		if(isteps.eq.stepStart.or.isteps.eq.stepEnd)then
			write(*,*)' Step number ',isteps
			write(*,*)'ModelOutputFile name ',ModelOutputFile		! only write out first and last file names
			endif
		open(40,FILE=ModelOutputFile,status = 'UNKNOWN')

		! Now read the model information
		isFileOpen = 1		! file is already open
!		call OpenModelOutputFile_2(isFileOpen,1)		! This one reads affinities. 1 means read the node and seg XY coordinates
		call OpenModelOutputFile(isFileOpen,1)			! No affinities in this one. 1 means read the node and seg XY coordinates

		totalCyclesToPlot(numToPlot) = TotalCycles
		if(TotalCycles.gt.maxCycles)maxCycles = totalCycles
		if(TotalCycles.lt.minCycles)minCycles = totalCycles
		! find the correct phase and store in arrays
! 		do i = 1,numNodesWith3Phases
! 			iNode = nodesWith3PhasesIndex(i)
		do iNode = 1,numNodes
			if(numToPlot.eq.1)then
				nodeXstart(iNode) = nodeX(iNode)		! This could also be the crystal center
				nodeYstart(iNode) = nodeY(iNode)
				endif
			nodeDistance(iNode,numToPlot) = sqrt((nodeX(iNode) - nodeXstart(iNode))**2 + 	&
				         (nodeY(iNode) - nodeYstart(iNode))**2)
			if(nodeDistance(iNode,numToPlot).gt.nodeDistanceMax)then
				nodeDistanceMax = nodeDistance(iNode,numToPlot)
				endif
			do j = 1,5		! zero out array (in case there is no grossular)
				phCompToPlot(iNode,numToPlot,j) = 0.0			
				end do
				
			select case(numNodeReactionPhases(iNode))
			case(3)
				do k = 1,3		! 4 phases in nodesWith3Phases. #1 is always the GB.
				!nodeMIFID(iNode,k) = MIFID
					MIFID = NodeMIFID(iNode,k)

					if(MIFID.eq.phaseToPlotMIFID)then
						if(isteps.eq.1)then		! only set these the first time through
							numNodesWithPhaseToPlot = numNodesWithPhaseToPlot + 1
							nodesWithPhaseToPlot(numNodesWithPhaseToPlot) = iNode
							endif
						phMolesToPlot(iNode,numToPlot) = nodePhaseMoles(iNode,k)
						do j = 1,numPhCo(MIFID)					! Note that this assumes order is prp Alm Sps Grs
							phCompToPlot(iNode,numToPlot,j) = nodePhaseComp(iNode,k,j)
							end do
						phCompToPlot(iNode,numToPlot,5) = &		! Fe/Fe+Mg assuming this is garnet
							nodePhaseComp(iNode,k,2)/(nodePhaseComp(iNode,k,2)+nodePhaseComp(iNode,k,1))
							!phCompToPlot(numToPlot,2)/(phCompToPlot(numToPlot,2)+phCompToPlot(numToPlot,1))
						if(phMolesToPlot(iNode,numToPlot).gt.maxMoles)maxMoles = phMolesToPlot(iNode,numToPlot)
						endif
					end do
			case(2)
	! 		do i = 1,numNodesWith2Phases
	! 			iNode = nodesWith2PhasesIndex(i)
				do kk = 1,2		! 3 phases in nodesWith2Phases. #1 is always the GB.
				!nodeMIFID(iNode,k) = MIFID
					k = nodeReactionPhases(iNode,kk)
					MIFID = NodeMIFID(iNode,k)

					if(MIFID.eq.phaseToPlotMIFID)then
						if(isteps.eq.1)then		! only set these the first time through
							numNodesWithPhaseToPlot = numNodesWithPhaseToPlot + 1
							nodesWithPhaseToPlot(numNodesWithPhaseToPlot) = iNode
							endif
						phMolesToPlot(iNode,numToPlot) = nodePhaseMoles(iNode,k)
						do j = 1,numPhCo(MIFID)
							phCompToPlot(iNode,numToPlot,j) = nodePhaseComp(iNode,k,j)
							end do
						phCompToPlot(iNode,numToPlot,5) = &		! Fe/Fe+Mg
							nodePhaseComp(iNode,k,2)/(nodePhaseComp(iNode,k,2)+nodePhaseComp(iNode,k,1))
							!phCompToPlot(numToPlot,2)/(phCompToPlot(numToPlot,2)+phCompToPlot(numToPlot,1))
						if(phMolesToPlot(iNode,numToPlot).gt.maxMoles)maxMoles = phMolesToPlot(iNode,numToPlot)
						endif
					end do
			case default
			end select
			end do

		end do		! go get the next step				

	! All steps are now in the arrays
	write(*,*)' Number of points to plot = ',numToPlot
	write(*,*)' Max mmoles in dataset    = ',maxMoles*1000.
	!x = (phMolesToPlot(i)*112.)		! volume in cubic centimeters
	!x = ((3./(4.*3.14159))*x)**0.33333	! radius in cm
	!x = x*1.0e4				! radius in µm

	maxRadius = 1.0e4*(maxMoles*112.*((3./(4.*3.14159))))**0.3333333	!radius in µm
	write(*,*)' Max radius in dataset     = ',maxRadius
	write(*,*)' Max Distance in dataset     = ',nodeDistanceMax
	write(*,*)' Min Cycles in dataset     = ',minCycles
	write(*,*)' Max Cycles in dataset     = ',maxCycles

11	continue
	write(*,*)' Plot type:'
	write(*,*)' 0 = return'
	write(*,*)' 1 = radius vs composition'
	write(*,*)' 2 = moles vs composition'
	write(*,*)' 3 = Model step vs composition'
	write(*,*)' 4 = Node position vs composition'
	write(*,*)' 5 = Cycle number vs composition'
	write(*,*)' 6 = Write garnet composition to output window.'
	read(*,*)plotType
	if(plotType.eq.0)go to 1

	xmin = 0.
	select case(plotType)
		case(1)
			xmax = maxRadius
			xLab = 'Radius'
		case(2)
			xmax = maxMoles*1000
			xLab = 'milliMoles'
		case(3)
			xmax = numToPlot
			xLab = 'Model_Step'
		case(4)
			xmax = 40.
			xLab = 'Node_Position'
		case(5)
			xmax = float(maxCycles)
			xLab = 'Cycles'
		case(6)		! write out garnet compositions -- actually, only Mn is coded here
			! calculate the maximum arbitrary units for the plot
			x = (phMolesToPlot(iNode,numToPlot)*112.)		! volume in cubic centimeters
			x = ((3./(4.*3.14159))*x)**0.33333	! radius in cm
			x = x*1.0e4				! radius in µm -- these are actually arbitrary units
			maxUnits = x
			write(*,*)' Input actual maximum radius in µm (calculated from NODE positions)'
			read(*,*)maxRadius
			write(12,*)'  '
			write(12,*)'  '
			write(12,*)' Spessartine composition vs radius'
			write(12,*)' Radius   Xsps'
			do i = 1,numToPlot
				y = phCompToPlot(iNode,i,3)		! This will only list spessartine composition (for diffusion modeling)
	
				x = (phMolesToPlot(iNode,i)*112.)		! volume in cubic centimeters
				x = ((3./(4.*3.14159))*x)**0.33333	! radius in cm
				x = x*1.0e4				! radius in µm -- these are actually arbitrary units
				write(12,81)x,x*maxRadius/maxUnits,y
	81			format(3F15.5)
				end do
			write(12,*)'  '
			write(12,*)'  '
			go to 11


		case default
			call FSS_Alert('ALERT','You chose poorly')
			go to 11
		end select
	xlen = 25.
	ylen = 20.
	nXdec = 1.
	nXstep = 10.


12	continue
	write(*,*)' Select components to plot'
	write(*,*)' 0 = exit (done)'
	write(*,*)' 1,3, or 4 = Prp + Sps + Grs'
	write(*,*)' 2 or 5 = Alm + Fe/Fe+Mg'
	write(*,*)' 200 = save Illustrator file'
	read(*,*)whatToPlot

	if(whatToPlot.eq.0)go to 11
	
	if(whatToPlot.eq.200)then
		call SaveIllustratorFile()
		go to 12
		endif

10	continue			
	do i = 1,numNodesWithPhaseToPlot
		if(nodesWithPhaseToPlot(i).gt.0)then
			write(*,*)nodesWithPhaseToPlot(i)
			endif
		end do
	write(*,*)' Choose the the node to plot --- enter 0 to exit'
	read(*,*)iNode
	if(iNode.eq.0)go to 12
	!iNode = nodesWithPhaseToPlot(iNodeIndex)

!	Set the X plotting variable -- I need to specify the node first
	do i = 1,numToPlot
		select case (plotType)
			case(1)			! plot radius on X
				x = (phMolesToPlot(iNode,i)*112.)		! volume in cubic centimeters
				x = ((3./(4.*3.14159))*x)**0.33333	! radius in cm
				xToPlot(i) = x*1.0e4				! radius in µm
			case(2)			! plot moles on X
				xToPlot(i) = 1000.*phMolesToPlot(iNode,i)
			case(3)			! plot cycle number (i.e. time) on X
				xToPlot(i) = i
			case(4)			! node position
				xToPlot(i) = nodeDistance(iNode,i)
			case(5)			! Cycle number
				xToPlot(i) = totalCyclesToPlot(i)
			case default
			end select
		write(12,*)i,xToPlot(i)
		end do

	
	write(*,*)' Use the same plot = 0; New plot = 1'
	read(*,*)reusePlot

	pen%penStyle = CanvasPenStyle_SolidLine

	select case (whatToPlot)

	case(2,5)	! Alm + Fe/Fe+Mg
		! Fe/(Fe+Mg) and almandine plot
		if(reusePlot.eq.1)then
			yLab = 'Alm_FeMg'
			ymin = 0.5
			ymax = 1.0
			nYstep = 10.
			nYdec = 2.
			SegPlot1%width =  1.25*(xLen * 28.34646)			! *xScale scales the X dimension
			SegPlot1%height = 1.25*(yLen * 28.34646)
			xor = 70.				! origin for X-Y plot (xmin,ymin) in pixels
			yor = SegPlot1%height - 70.
			PlTitle = 'Alm_FeMg'
			call Setplot(iOK)
			close(PSUNIT)         			! close old PostScript scratch file
			call psopcl(4)    			! open a new scratch file
			CALL USER()			! sets the user coordinates that were input in Sub SetPlot
			CALL AWE_createCanvas(SegPlot1)	
			call axis(SegPlot1)			! draws the axis in AWE
			call paxis()
			dx = .002*(xmax-xmin)
			dy = .002*(ymax-ymin)
			endif
		
		! Alm
		currentColor = 2	! red
		pen%penColor = AWEColorNumber(CurrentColor)
		brush%brushColor = AWEColorNumber(CurrentColor)
		j = 2
		iup = 0
		do i = 1,numToPlot
			y = phCompToPlot(iNode,i,j)
			x = xToPlot(i)
			call plot(SegPlot1,x,y,iup)
			Call PlotCenteredEllipseOnScreen(SegPlot1,X,dX,Y,dY,brush,pen)
			call pplot(x,y,iup)
			iup = 1
			end do
		call ppenup()

		! Fe/Fe+Mg
		currentColor = 1	! Black
		pen%penColor = AWEColorNumber(CurrentColor)
		brush%brushColor = AWEColorNumber(CurrentColor)
		j = 5
		iup = 0
		do i = 1,numToPlot
			y = phCompToPlot(iNode,i,j)
			x = xToPlot(i)
			call plot(SegPlot1,x,y,iup)
			Call PlotCenteredEllipseOnScreen(SegPlot1,X,dX,Y,dY,brush,pen)
			call pplot(x,y,iup)
			iup = 1
			end do
		call ppenup()


	case(1,3,4)	! Prp + Sps + Grs	
		! Prp, Sps, Grs plot
		if(reusePlot.eq.1)then
			yLab = 'Prp Sps Grs'
			ymin = 0.0
			ymax = 0.3
			nYstep = 15.
			nYdec = 2.
! 			SegPlot2%title = ''
! 			SegPlot2%maxwindowwidth = -1
! 			SegPlot2%maxwindowheight = -1
			SegPlot2%width =  int(1.25*(xLen * 28.34646))			! *xScale scales the X dimension
			SegPlot2%height = int(1.25*(yLen * 28.34646))
			xor = 70.				! origin for X-Y plot (xmin,ymin) in pixels
			yor = SegPlot2%height - 70.
			PlTitle = 'PrpSpsGrs'
			call Setplot(iOK)
			close(PSUNIT)         			! close old PostScript scratch file
			call psopcl(4)    			! open a new scratch file
			CALL USER()			! sets the user coordinates that were input in Sub SetPlot
			CALL AWE_createCanvas(SegPlot2)	
			call axis(SegPlot2)			! draws the axis in AWE
			call paxis()
			dx = .002*(xmax-xmin)
			dy = .002*(ymax-ymin)
			endif	

		! Prp
		currentColor = 5	! Cyan
		pen%penColor = AWEColorNumber(CurrentColor)
		brush%brushColor = AWEColorNumber(CurrentColor)
		j = 1
		iup = 0
		do i = 1,numToPlot
			y = phCompToPlot(iNode,i,j)
			x = xToPlot(i)
			call plot(SegPlot2,x,y,iup)
			Call PlotCenteredEllipseOnScreen(SegPlot2,X,dX,Y,dY,brush,pen)
			call pplot(x,y,iup)
			iup = 1
			end do
		call ppenup()

		!Sps
		currentColor = 7	! Magenta
		pen%penColor = AWEColorNumber(CurrentColor)
		brush%brushColor = AWEColorNumber(CurrentColor)
		j = 3
		iup = 0
		do i = 1,numToPlot
			y = phCompToPlot(iNode,i,j)
			x = xToPlot(i)
			call plot(SegPlot2,x,y,iup)
			Call PlotCenteredEllipseOnScreen(SegPlot2,X,dX,Y,dY,brush,pen)
			call pplot(x,y,iup)
			iup = 1
			end do
		call ppenup()

		!Grs
		currentColor = 4	! green
		pen%penColor = AWEColorNumber(CurrentColor)
		brush%brushColor = AWEColorNumber(CurrentColor)
		j = 4
		iup = 0
		do i = 1,numToPlot
			y = phCompToPlot(iNode,i,j)
			x = xToPlot(i)
			call plot(SegPlot2,x,y,iup)
			Call PlotCenteredEllipseOnScreen(SegPlot2,X,dX,Y,dY,brush,pen)
			call pplot(x,y,iup)
			iup = 1
			end do
		call ppenup()


	case default
		call FSS_Alert('ALERT','You chose poorly')
		go to 12
	end select
		
	currentColor = 1
	go to 10

	end
		
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine PlotAffinities(SegPlot1,iNode)
!	Routine to plot compositions along segments defined by a set of nodes
!	Nodes must be contiguous or the plots won't make much sense
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: SegPlot1
! 	TYPE(AWE_Canvas) :: SegPlot2
! 	TYPE(AWE_Canvas) :: SegPlot3
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	include "GibbsFiles.inc"
	real*4 x,y,dX,dY
	integer*4 i,iup,k,kk,iNode,status,stepStart,stepEnd,isteps,numToPlot,MIFID,iOK
	real*4 phAffToPlot(maxNodes,4000,4),maxAffinity,minAffinity,Au_to_Ox
	character*16 stepText
	integer*4 reusePlot,plottype
	integer*4 totalCyclesToPlot(4000),maxCycles,minCycles
	common /plotstorage/phAffToPlot,totalCyclesToPlot
	
	CurrentColor = 1		! black is the default



	! Open and read the model information file
	write(*,*)' Open the model input file'
	!call FSS_Alert('ALERT','Open model input file')
	call GibbsBegin()


	
	! We already have a base. Just open this base file name
	write(*,*)' Open the model file base name'
	call FSS_Alert('ALERT','Open model file base name')
	open(73,FILE='',status='OLD',iostat=status,action='WRITE')
	if(status.ne.0)then
		call FSS_Alert('Alert','Problem opening model base file')
		return
		endif
	inquire(73,NAME=ModelOutputFileBase)
	close(73)	! we don't need this file for anything

1	continue


	write(*,*)'Input starting step for making the plot (0 to begin at the beginning)'
	read(*,*)stepStart
	write(*,*)'Input ending step for making the plot '
	read(*,*)stepEnd

	if(stepEnd-stepStart.ge.3999)then
		call FSS_Alert('ALERT',' Number of steps exceeds array allocation of 4000')
		return
		endif

!	read all the models and collect data for the plots
	numToPlot = 0
	do isteps = stepStart,stepEnd
		numToPlot = numToPlot + 1
		write(steptext,55)isteps
55		format(I5)
		steptext = adjustL(steptext)
		ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
		if(isteps.eq.stepStart.or.isteps.eq.stepEnd)then
			write(*,*)' Step number ',isteps
			write(*,*)'ModelOutputFile name ',ModelOutputFile		! only write out first and last file names
			endif
		open(40,FILE=ModelOutputFile,status = 'UNKNOWN')

		! Now read the model information
		isFileOpen = 1		! file is already open
!		call OpenModelOutputFile_2(isFileOpen,1)		! This one reads affinities. 1 means read the node and seg XY coordinates
		call OpenModelOutputFile(isFileOpen,1)			! No affinities in this one. 1 means read the node and seg XY coordinates

		totalCyclesToPlot(numToPlot) = TotalCycles
		! find the correct phase and store in arrays
! 		do i = 1,numNodesWith3Phases
! 			iNode = nodesWith3PhasesIndex(i)
		do iNode = 1,numNodes
! 			select case(numNodeReactionPhases(iNode))
! 			case(3)
! 				do kk = 1,3		! 4 phases in nodesWith3Phases. #1 is always the GB.
! 					k = kk+1
! 					phAffToPlot(iNode,numToPlot,kk) = nodePhaseAffinity(iNode,k) 
! 					end do
! 			case(2)
! 				do kk = 1,2		! 3 phases in nodesWith2Phases. #1 is always the GB.
! 					k = kk+1
! 					phAffToPlot(iNode,numToPlot,kk) = nodePhaseAffinity(iNode,k) 
! 					end do
! 			case default
! 			end select
! 			end do

! 			do k = 1,numNodeReactionPhases(iNode)
! 				kk = nodeReactionPhases(iNode,k)
! 				phAffToPlot(iNode,numToPlot,k) = nodePhaseAffinity(iNode,kk)
! 				end do 
! 			end do

			do kk = 1,numNodeReactionPhases(iNode)
				k = nodeReactionPhases(iNode,kk)
				MIFID = nodeMIFID(iNode,k)
				if(MINREC(MIFID).eq.32)then	! garnet
					AU_to_Ox = 0.66666	! 8/12
					go to 15
					endif
				if(MINREC(MIFID).eq.123)then	! chlorite
					AU_to_Ox = 0.77777	! 14/18
					go to 15
					endif
				if(MINREC(MIFID).eq.1)then	! quartz
					AU_to_Ox = 0.5		! 1/2
					go to 15
					endif
				if(MINREC(MIFID).eq.17)then	! muscovite
					AU_to_Ox = 0.66666	! 16/24
					go to 15
					endif
				if(MINREC(MIFID).eq.122)then	! biotite
					AU_to_Ox = 0.66666	! 16/24
					go to 15
					endif
				if(MINREC(MIFID).eq.93)then	! plagioclase
					AU_to_Ox = 0.625	! 5/8
					go to 15
					endif


				call FSS_Alert('ALERT','Did not find MINREC(MIFID) in routine PlotAffinities line 1372 or so')
				write(*,*) 'Error dump'
				write(*,*)'iNode,RxnPh=k,MIFID,MINREC(MIFID)'
				write(*,*)iNode,k,MIFID,MINREC(MIFID)
				pause 'pausing --- error abort'

15				continue					
				phAffToPlot(iNode,numToPlot,kk) = -nodePhaseAffinity(iNode,k)*AU_to_Ox
				end do 
			end do




		end do		! go get the next step				

	! All steps are now in the arrays
	write(*,*)' Number of points to plot = ',numToPlot
	minCycles = totalCyclesToPlot(1)
	maxCycles = totalCyclesToPlot(numToPlot)
	if(TotalCycles.lt.minCycles)minCycles = totalCycles
	write(*,*)' Min Cycles in dataset     = ',minCycles
	write(*,*)' Max Cycles in dataset     = ',maxCycles

11	continue
	write(*,*)' Plot type:'
	write(*,*)' 0 = return'
	write(*,*)' 1 = Model step vs Affinity'
	write(*,*)' 2 = Cycle number vs Affinity'
	write(*,*)' 200 = save Illustrator file'
! 	case(200)
! 		call SaveIllustratorFile()

	read(*,*)plotType
	if(plotType.eq.0)return

	if(plotType.eq.200)then
		call SaveIllustratorFile()
		go to 11
		endif

10	continue

	write(*,*)' Pick node to plot -- enter 0 to exit'
	write(*,*)' Number of nodes total = ',numNodes
	read(*,*)iNode
	if(iNode.eq.0)go to 11

	write(*,*)' '
	write(*,*)'iNode = ',iNode
	write(*,*)' MIFID        PhName'
	do k = 1,numNodeReactionPhases(iNode)
		kk = nodeReactionPhases(iNode,k)
		MIFID = nodeMIFID(iNode,kk)
		write(*,*)MIFID,phName(MIFID)
		end do
	write(*,*)' Look OK? 0 = NO, 1 = YES (continue)'
	read(*,*)iok
	if(iok.eq.0)go to 10
	
	maxAffinity = -1.0e-10
	minAffinity =  1.e10
	write(12,*)' '
	write(12,*)' '
	write(12,*)' Node = ',iNode
	write(12,*)(phName(nodeMIFID(iNode,nodeReactionPhases(iNode,k))),k=1,numNodeReactionPhases(iNode))	
	do i = 1,numToPlot
		do k = 1,numNodeReactionPhases(iNode)
			if(phAffToPlot(iNode,i,k).gt.maxAffinity)maxAffinity = phAffToPlot(iNode,i,k)
			if(phAffToPlot(iNode,i,k).lt.minAffinity)minAffinity = phAffToPlot(iNode,i,k)
			end do
		write(12,*)(phAffToPlot(iNode,i,k),k=1,numNodeReactionPhases(iNode))
		end do

	write(*,*)' '
	write(*,*)' minAffinity = ',minAffinity
	write(*,*)' maxAffinity = ',maxAffinity
	write(*,*)' '
	
	write(*,*)' Use the same plot = 0; New plot = 1'
	read(*,*)reusePlot
	if(reusePlot.eq.0) go to 20

!	Initialize Illustrator File
	
	xmin = 0.
	select case(plotType)
		case(1)
			xmax = numToPlot - 1
			xLab = 'Model_Step'
		case(2)
			xmax = float(maxCycles)
			xLab = 'Cycles'
		case default
			call FSS_Alert('ALERT','You chose poorly')
			go to 11
		end select

	xlen = 25.
	ylen = 20.
	nXdec = 0
	nXstep = 10.
	yLab = 'Affinity'

	ymin = minAffinity
	ymax = maxAffinity
	nYstep = 10.
	nYdec = 0
	SegPlot1%width =  1.25*(xLen * 28.34646)			! *xScale scales the X dimension
	SegPlot1%height = 1.25*(yLen * 28.34646)
	xor = 70.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = SegPlot1%height - 70.
	PlTitle = 'Affinity'
	call Setplot(iOK)

	close(PSUNIT)         			! close old PostScript scratch file
	call psopcl(4)    			! open a new scratch file

      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	CALL AWE_createCanvas(SegPlot1)	
      	call axis(SegPlot1)			! draws the axis in AWE
	call paxis()				! Draw axis for illustrator file

	!pause 'Hit return to continue '
	
20	continue
	dx = .002*(xmax-xmin)
	dy = .002*(ymax-ymin)

	currentColor = 1	! Black
	call USER()

	! Note that these colors correspond to the colors in the "Element color list" in the model input file -- not such a great way to do things.
	do k = 1,numNodeReactionPhases(iNode)
		kk = nodeReactionPhases(iNode,k)
		MIFID = nodeMIFID(iNode,kk)
		write(*,*)MIFID,phName(MIFID)
		if(trim(phName(MIFID)).eq.'Garnet')then
			currentColor = 7		!magenta
			write(*,*)MIFID,phName(MIFID),currentColor
			endif
		if(trim(phName(MIFID)).eq.'Quartz')then
			currentColor = 5		! Cyan	
			write(*,*)MIFID,phName(MIFID),currentColor
			endif
		if(trim(phName(MIFID)).eq.'Chlorite')then
			currentColor = 4		! Green
			write(*,*)MIFID,phName(MIFID),currentColor
			endif
		if(trim(phName(MIFID)).eq.'Biotite')then
			currentColor = 10		! Brown
			write(*,*)MIFID,phName(MIFID),currentColor
			endif
		if(trim(phName(MIFID)).eq.'Muscovite')then
			currentColor = 24		! Pink
			write(*,*)MIFID,phName(MIFID),currentColor
			endif
		if(trim(phName(MIFID)).eq.'Plagioclase')then
			currentColor = 3		! Yellow
			write(*,*)MIFID,phName(MIFID),currentColor
			endif
		
		pen%penStyle = CanvasPenStyle_SolidLine
		pen%penColor = AWEColorNumber(CurrentColor)
		brush%brushColor = AWEColorNumber(CurrentColor)

		iup = 0	
		do i = 1,numToPlot
			y = phAffToPlot(iNode,i,k)

			select case (plotType)
			case(1)			! plot radius on X
				x = i
			case(2)			! plot moles on X
				x = totalCyclesToPlot(i)
			case default
			end select

			call plot(SegPlot1,x,y,iup)
			Call PlotCenteredEllipseOnScreen(SegPlot1,X,dX,Y,dY,brush,pen)	
			call pplot(x,y,iup)
			iup = 1
			end do
		call ppenup()		! puts an "S" at the end of the line in illustrator file
		end do
	currentColor = 1
	go to 10
	end
		
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine PlotAffinity_At_Node(SegPlot1,iNode)
!	Routine to plot the affinity for garnet nucleation at a node that does not include garnet
!	This does the same calculation of affinity as in routine PlotAffinityGrid except here
!		the sequence of steps is read in for the node so the plot can be against the cycle number
!	The purpose is to examine how the affinity decreases away from the garnet reaction 
!		where the affinity is entirely controlled by GB diffusion
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: SegPlot1
! 	TYPE(AWE_Canvas) :: SegPlot2
! 	TYPE(AWE_Canvas) :: SegPlot3
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	include "GibbsFiles.inc"
	include "GB_Tangent.inc"
	real*4 x,y,dX,dY
	integer*4 i,j,iup,k,kk,iNode,status,stepStart,stepEnd,isteps,numToPlot,MIFID,iOK,iGarnet,izero
	real*4 phAffToPlot(4000),maxAff,minAff,range
	character*16 stepText
	integer*4 reusePlot,plottype
	integer*4 totalCyclesToPlot(4000),maxCycles,minCycles
	common /plotstorage/phAffToPlot,totalCyclesToPlot
	
	CurrentColor = 1		! black is the default



	! Open and read the model information file
	write(*,*)' Open the model input file'
	!call FSS_Alert('ALERT','Open model input file')
	call GibbsBegin()


	
	! We already have a base. Just open this base file name
	write(*,*)' Open the model file base name'
	call FSS_Alert('ALERT','Open model file base name')
	open(73,FILE='',status='OLD',iostat=status,action='WRITE')
	if(status.ne.0)then
		call FSS_Alert('Alert','Problem opening model base file')
		return
		endif
	inquire(73,NAME=ModelOutputFileBase)
	close(73)	! we don't need this file for anything

	write(*,*)' Which phase is Garnet?'
	do i = 1,numPhMIF
		write(*,*)i,phName(i)
		end do
	read(*,*)iGarnet		! this is the number of the phase in the MIF
	!phaseToPlotName = phName(phaseToPlotMIFID)
	if(iGarnet.eq.0)return

1	continue

	write(*,*)'Input starting step for making the plot (0 to begin at the beginning)'
	read(*,*)stepStart
	write(*,*)'Input ending step for making the plot '
	read(*,*)stepEnd

	if(stepEnd-stepStart.ge.3999)then
		call FSS_Alert('ALERT',' Number of steps exceeds array allocation of 4000')
		return
		endif

10	continue

	write(*,*)' Pick node to plot -- enter 0 to exit'
	write(*,*)' Number of nodes total = ',numNodes
	read(*,*)iNode
	if(iNode.eq.0)return


!	read all the models and collect data for the plots
	! Calculate the Affinities for every point and determine range of affinity values
	maxAff = -1.0e5
	minAff = 1.0d5		! make very large to start
	write(12,*)'**********************************************************'
	write(12,*)'**********************************************************'
	numToPlot = 0
	write(*,*)' Reading model files...please be patient'
	do isteps = stepStart,stepEnd
		numToPlot = numToPlot + 1
		write(steptext,55)isteps
55		format(I5)
		steptext = adjustL(steptext)
		ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
		if(isteps.eq.stepStart.or.isteps.eq.stepEnd)then
			write(*,*)' Step number ',isteps
			write(*,*)'ModelOutputFile name ',ModelOutputFile		! only write out first and last file names
			endif
		open(40,FILE=ModelOutputFile,status = 'UNKNOWN')

		! Now read the model information
		isFileOpen = 1		! file is already open
!		call OpenModelOutputFile_2(isFileOpen,1)		! This one reads affinities. 1 means read the node and seg XY coordinates
		call OpenModelOutputFile(isFileOpen,1)			! No affinities in this one. 1 means read the node and seg XY coordinates

		totalCyclesToPlot(numToPlot) = TotalCycles

		do j = 1,numPhCo(1)		! grain boundary
			xPhCo(1,j) = nodeComp(iNode,j)		! numPhCo(1) should be every element
			end do

		Call ParallelToTangent(iGarnet,0,izero)
		phAffToPlot(numToPlot) = -gDifferenceAU(iGarnet)*0.666666	! This converts from AU to J/mol-O	
		write(12,81)isteps,-gDifferenceAU(iGarnet)*0.66666,(xPhCo(iGarnet,j),j=1,numPhCo(iGarnet))
81		format(I8,T15,15F15.5)

		if(phAffToPlot(numToPlot).gt.maxAff)then
			maxAff = phAffToPlot(numToPlot)
			endif
		if(phAffToPlot(numToPlot).lt.minAff)then
			minAff = phAffToPlot(numToPlot)
			endif
		end do		! go get the next step				
		
		write(*,*)' Minimum = ',minAff
		write(*,*)' Maximum = ',maxAff
		range = (maxAff-minAff)
		write(*,*)' Range = ',range



	! All steps are now in the arrays
	write(*,*)' Number of points to plot = ',numToPlot
	minCycles = totalCyclesToPlot(1)
	maxCycles = totalCyclesToPlot(numToPlot)
	if(TotalCycles.lt.minCycles)minCycles = totalCycles
	write(*,*)' Min Cycles in dataset     = ',minCycles
	write(*,*)' Max Cycles in dataset     = ',maxCycles

11	continue
	write(*,*)' Plot type:'
	write(*,*)' 0 = return'
	write(*,*)' 1 = Model step vs Affinity'
	write(*,*)' 2 = Cycle number vs Affinity'
	read(*,*)plotType
	if(plotType.eq.0)go to 10

	write(*,*)' '
	write(*,*)'iNode = ',iNode
	write(*,*)' MIFID        PhName'
	do k = 1,numNodeReactionPhases(iNode)
		kk = nodeReactionPhases(iNode,k)
		MIFID = nodeMIFID(iNode,kk)
		write(*,*)MIFID,phName(MIFID)
		end do
	write(*,*)' Look OK? 0 = NO, 1 = YES (continue)'
	read(*,*)iok
	if(iok.eq.0)go to 10
	
	
	write(*,*)' Use the same plot = 0; New plot = 1'
	read(*,*)reusePlot
	if(reusePlot.eq.0) go to 20

	xmin = 0.
	select case(plotType)
		case(1)
			xmax = numToPlot - 1
			xLab = 'Model_Step'
		case(2)
			xmax = float(maxCycles)
			xLab = 'Cycles'
		case default
			call FSS_Alert('ALERT','You chose poorly')
			go to 11
		end select

	xlen = 25.
	ylen = 20.
	nXdec = 0
	nXstep = 10.
	yLab = 'Affinity'

	ymin = minAff
	ymax = maxAff
	nYstep = 10.
	nYdec = 0
	SegPlot1%width =  1.25*(xLen * 28.34646)			! *xScale scales the X dimension
	SegPlot1%height = 1.25*(yLen * 28.34646)
	xor = 70.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = SegPlot1%height - 70.
	PlTitle = 'Affinity'
	call Setplot(iOK)
	close(PSUNIT)         			! close old PostScript scratch file
	call psopcl(4)    			! open a new scratch file
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	CALL AWE_createCanvas(SegPlot1)	
      	call axis(SegPlot1)			! draws the axis in AWE
	call paxis()

	!pause 'Hit return to continue '
	
20	continue
	dx = .002*(xmax-xmin)
	dy = .002*(ymax-ymin)

	currentColor = 1	! Black
!	call USER()

	! Note that these colors correspond to the colors in the "Element color list" in the model input file -- not such a great way to do things.
	currentColor = 7		!magenta
		
	pen%penStyle = CanvasPenStyle_SolidLine
	pen%penColor = AWEColorNumber(CurrentColor)
	brush%brushColor = AWEColorNumber(CurrentColor)

	iup = 0	
	do i = 1,numToPlot
		y = phAffToPlot(i)

		select case (plotType)
		case(1)			! plot radius on X
			x = i
		case(2)			! plot moles on X
			x = totalCyclesToPlot(i)
		case default
		end select

		call plot(SegPlot1,x,y,iup)
		Call PlotCenteredEllipseOnScreen(SegPlot1,X,dX,Y,dY,brush,pen)	
		call pplot(x,y,iup)
		iup = 1
		end do
	call ppenup()

	currentColor = 1

	call SaveIllustratorFile()

	go to 10
	end
