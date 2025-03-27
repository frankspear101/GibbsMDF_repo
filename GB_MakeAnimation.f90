	Subroutine MakeAnimation()
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen

	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	include "GibbsFiles.inc"
	integer*4 whatToPlot,i,elToPlot,status,stepStart,stepEnd,isteps,iMin,iNode,debug,iSeg,NodeToCheck
	real*4 minEl,maxEl,minMoles,maxMoles
	real*8 moleScaleFactor
	character*16 stepText

	! Routine to open a series of output files, generate GridPlots and save as png files.
	! This will allow animations


1	continue
	write(*,*)'Choose what to plot:'
	write(*,*)' 0 = exit'
	write(*,*)' 1 = moles of phases'
	write(*,*)' 2 = grain boundary composition'
	write(*,*)' 3 = Plot grid evolution'
	write(*,*)' 4 = Test grow all nodes and segs'
	write(*,*)' 5 = Test grow a single node'
	read(*,*)whatToPlot

	select case (whatToPlot)
	case(0)
		return
	case(1)
		write(*,*)'Specify minimum and maximum for scaling (better look at the last model file first)'
		write(*,*)'Minimum = '
		read(*,*)minMoles
		write(*,*)'Maximum = '
		read(*,*)maxMoles
	case(2)
		write(*,*)' Which element to plot?'
		do i = 1,numEl
			write(*,*)i,elName(i)
			end do
		read(*,*)elToPlot
		write(*,*)'Specify minimum and maximum for scaling (better look at the last model file first)'
		write(*,*)'Minimum = '
		read(*,*)minEl
		write(*,*)'Maximum = '
		read(*,*)maxEl
	case(3)
! 		write(*,*)' Do you want to plot points or just the grid?'
! 		write(*,*)' 0 = grid only; 1 = grid and points'
! 		read(*,*)pointsOK
		
	case(4,5)
		call FSS_Alert('Alert',' Open initial model output file')			
		call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the node and seg XY coordinates

		write(*,*)'Input moleScaleFactor (e.g. 10-100 -- but test first to see if it will work)'
		write(*,*)' You can test in the MDF Run Model option 16, 17 or 18'
		write(*,*)' Input 0 to abort'
		read(*,*)moleScaleFactor
		if(moleScaleFactor.eq.0)go to 1
		if(whatToPlot.eq.5)then
			write(*,*)'Input node to check. '
			read(*,*)nodeToCheck
			endif
		write(*,*)'Debug? 0 = no, 1 = yes'
		read(*,*)debug
	case default
		call FSS_Alert('ALERT','You chose poorly')
		go to 1
	end select

	! Open and read the model information file
	write(*,*)' Open the model input file'
	!call FSS_Alert('ALERT','Open model input file')
	call GibbsBegin()
	
	! We already have a base. Just open this base file name
	write(*,*)' Open the model file base name'
	call FSS_Alert('ALERT','Open model file base name')
	open(73,FILE='',status='OLD',iostat=status,action='WRITE')
	if(status.ne.0)then
		call FSS_Alert('Alert','Problem opening model file')
		return
		endif
	inquire(73,NAME=ModelOutputFileBase)
	close(73)	! we don't need this file for anything
	write(*,*)'Input starting step for animation (0 to begin at the beginning)'
	read(*,*)stepStart
	write(*,*)'Input ending step for animation '
	read(*,*)stepEnd
	if(whatToPlot.eq.4)then
		stepStart = stepStart + 1	! we don't want to read the first one again
		endif

	! try clearing canvas rather than closing and reopening
	call SetUpGridPlot(GridPlot,maxX,maxY)
	pause 'Should be a blank canvas'
	do isteps = stepStart,stepEnd
		write(*,*)' Step number ',isteps
		write(steptext,55)isteps
55		format(I5)
		steptext = adjustL(steptext)
		ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
		write(*,*)'ModelOutputFile name ',ModelOutputFile
		open(40,FILE=ModelOutputFile,status = 'UNKNOWN')

		! Now read the model information
		isFileOpen = 1		! file is already open
		if(whatToPlot.eq.4)then
			call OpenModelOutputFile(isFileOpen,0)		! 0 = do not read the node or seg XY coordinates (keep the last)
			else
			call OpenModelOutputFile(isFileOpen,1)		! 1 = do read the node or seg XY coordinates (overwrite the last)
			endif				
		select case (whatToPlot)
		case(1)
			! now draw the Grid and plot the moles
			call DrawCrystals(GridPlot)
			call LabelCrystals(GridPlot)
			call PlotGrid(GridPlot,0)
			iMin = 2
			call Plot_Mineral_Moles_Scaled(GridPlot,1,iMin,minMoles,maxMoles)
			!pause 'Take a look -- then hit return'
			iMin = 4
			call Plot_Mineral_Moles_Scaled(GridPlot,1,iMin,minMoles,maxMoles)
			!pause 'Take a look -- then hit return'
			iMin = 5
			call Plot_Mineral_Moles_Scaled(GridPlot,1,iMin,minMoles,maxMoles)
			!pause 'Take a look -- then hit return'
			iMin = 6
			call Plot_Mineral_Moles_Scaled(GridPlot,1,iMin,minMoles,maxMoles)
			!pause 'Take a look -- then hit return'
			iMin = 7
			call Plot_Mineral_Moles_Scaled(GridPlot,1,iMin,minMoles,maxMoles)
			!pause 'Take a look -- then hit return'
			iMin = 8	
			call Plot_Mineral_Moles_Scaled(GridPlot,1,iMin,minMoles,maxMoles)
			!pause 'Take a look -- then hit return'
		case(2)
			! now draw the Grid and plot the moles
			call DrawCrystals(GridPlot)
			call LabelCrystals(GridPlot)
			call PlotGrid(GridPlot,0)
			call Plot_GB_Comp(GridPlot,1,ElToPlot,minEl,maxEl)
		case(3)
			! We need to go through every node (and segs after code is written) and plot new position
			write(*,*)'iStep ',iSteps
			!call PlotGrid(GridPlot,0)
			!if(pointsOK.eq.1)then
			call PlotPoints(GridPlot,0,1,1)		! 0 = do not pick color, 1 =do connect the dots, 1 = do draw dots
			!endif
		case(4)
			! We need to go through every node (and segs after code is written) and calculate new position
			write(*,*)'iStep ',iSteps
			do iNode = 1,numNodes
				select case (numNodeReactionPhases(iNode))
					case(3)
						write(*,*)''
						write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
						write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
						write(*,*)'3Node ',iNode
						call Grow3Nodes(iNode,moleScaleFactor,debug)
					case(2)
						write(*,*)''
						write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
						write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
						write(*,*)'2Node ',iNode
						call Grow2Nodes(iNode,moleScaleFactor,debug)
					case default
					end select
				end do
			do iSeg = 1,numSegs
				if(numSegReactionPhases(iSeg).eq.2)then
					write(*,*)'iSeg  ',iSeg
					call GrowSegs(iSeg,moleScaleFactor,debug)
					endif
				end do			
			do iNode = 1,numNodes
				call RemoveSegPoints(iNode,debug)
				end do			
			call PlotGrid(GridPlot,0)
		case(5)
			! We need to go through every step for this node
			write(*,*)'iStep ',iSteps
			iNode = NodeToCheck
				select case (numNodeReactionPhases(iNode))
					case(3)
						write(*,*)''
						write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
						write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
						write(*,*)'3Node ',iNode
						call Grow3Nodes(iNode,moleScaleFactor,debug)
					case(2)
						write(*,*)''
						write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
						write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
						write(*,*)'2Node ',iNode
						call Grow2Nodes(iNode,moleScaleFactor,debug)
					case default
					end select

! 			do iSeg = 1,numSegs
! 				if(numSegReactionPhases(iSeg).eq.2)then
! 					write(*,*)'iSeg  ',iSeg
! 					call GrowSegs(iSeg,moleScaleFactor,debug)
! 					endif
! 				end do			

			call RemoveSegPoints(iNode,debug)

!			call PlotGrid(GridPlot,0)
			
		case default
		end select
			
		! now save the canvas
		!pause 'Take a look -- then hit return'
		CanvasFileName = trim(ModelOutputFileBase)//'__'//trim(steptext)
		!CALL AWE_saveCanvas(GridPlot,”canvas.png”)
		write(*,*)'Canvas File name ',CanvasFileName
		CALL AWE_saveCanvas(GridPlot,trim(CanvasFileName))
		if(isteps.ne.stepend)then			! don't erase the last one
			CALL AWE_clearCanvas(GridPlot)
			endif
		!CALL AWE_CloseCanvas(GridPlot)

		end do

	go to 1
	end
		
	
	
	
