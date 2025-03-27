	Subroutine GridPlotRoutines(GridPlot)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	include "Diffuse_GB_Hex.inc"
	include "Assemb.inc"
	include "GibbsFiles.inc"
	integer*4 option,iunit,i,k,iNode,iSeg,another,PCN,debug
	integer*4 status,stepstart,stepend,isteps,stepInterval
	character*5 steptext
	real*4 plotWidth
	real*8 moleScaleFactor




1	continue
	write(*,*)' 0 = return'
	write(*,*)' 1 = SetUpGridPlot'
	write(*,*)' 10 = Auto GridPlot'
	write(*,*)' 101 = Auto GridPlot - super automated'
	write(*,*)' 11 = SetUpGridPlot2 (allows for variable scaling on X and Y)'
	write(*,*)' 12 = SetUpGridPlot3 around a specific node (auto scaling on X and Y)'
	write(*,*)' 2 = plot nodes'
	write(*,*)' 3 = plot points'
	write(*,*)' 4 = label nodes'
	write(*,*)' 5 = label segments'
	write(*,*)' 6 = Draw Crystals'
	write(*,*)' 7 = Label Crystals'
	write(*,*)' 8 = Plot grid compositions'
	write(*,*)'80 = Plot grid compositions -- version 2 (Open model output file(21) then 1 and 6 first)'
	write(*,*)' 9 = Plot grid moles'
	write(*,*)'13 = Calculate and Plot GARNET affinities on GRID (Open model output file(21) then 1 and 6 first)'
	write(*,*)'14 = Make Garnet Zoning Map'
	write(*,*)'17 = Move Segments test'
	write(*,*)'18 = List segment information'
	write(*,*)'20 = Open model input file'
	write(*,*)'21 = Open model output file'
	write(*,*)'22 = Animate grid plots'
	write(*,*)'23 = DumpGrid'
	write(*,*)'30 = Scale crystals'
	write(*,*)'31 = Check crystal segs'
	write(*,*)'32 = Calculate bulk composition from crystal areas'
	write(*,*)'33 = Calculate crystal moles'
	write(*,*)'34 = Test moles to area routine'
	write(*,*)'35 = Check node information'
	write(*,*)'36 = Calculate crystal areas'
	write(*,*)'37 = Check crystal segs and nodes'
	write(*,*)'50 = Muscovite/Biotite ratios at nodes'
	write(*,*)'200 = Save Illustrator file for grid'
	
	read(*,*)option
	select case (option)
	case(0)
		return
	case(1)
		call SetUpGridPlot(GridPlot,maxX,maxY)
	case(10)
		write(*,*)' Input value for plot width in cm (30-50)'
		read(*,*)plotWidth
!		write(*,*)'Open output file to plot '
1010		continue
		call FSS_Alert('Alert','Select output file to plot')
		call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the node and seg XY coordinates
		call SetUpGridPlot10(GridPlot,maxX,maxY,plotWidth)
		call DrawCrystals2(GridPlot)
		call PlotPoints(GridPlot,0,1,1)		! one means pick color; 0 = ask to connect the points
		write(*,*)'Plot another? 0 = no; 1 = yes'
		read(*,*,err=1)another
		if(another.eq.1)then
			go to 1010
			else
			go to 1
			endif
	case(101)
		call GibbsBegin()
		write(*,*)' Input value for plot width in cm (30-50)'
		read(*,*)plotWidth
!		write(*,*)'Open output file to plot '
		call FSS_Alert('Alert','Select base output file ')
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
		write(*,*)'Input step interval'
		read(*,*)stepInterval
		call SetUpGridPlot10(GridPlot,maxX,maxY,plotWidth)

		do isteps = stepStart,stepEnd,stepInterval
			write(*,*)' Step number ',isteps
			write(steptext,55)isteps
	55		format(I5)
			steptext = adjustL(steptext)
			ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
			write(*,*)'ModelOutputFile name ',ModelOutputFile
			open(40,FILE=ModelOutputFile,status = 'UNKNOWN')	! open the model output file
			! Now read the model information
			isFileOpen = 1		! file is already open
			call OpenModelOutputFile(isFileOpen,1)		! 1 = do read the node or seg XY coordinates (overwrite the last)
			!call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the node and seg XY coordinates
			!call SetUpGridPlot10(GridPlot,maxX,maxY,plotWidth)
			call DrawCrystals2(GridPlot)
			call PlotPoints(GridPlot,0,1,1)		! one means pick color; 0 = ask to connect the points
			! now save the canvas
			CanvasFileName = trim(ModelOutputFileBase)//'___'//trim(steptext)
			!CALL AWE_saveCanvas(GridPlot,”canvas.png”)
			!CanvasFileName = trim(ModelOutputFile)
			!write(*,*)'Canvas File name ',ModelOutputFile
			write(*,*)'Canvas File name ',CanvasFileName
			!CALL AWE_saveCanvas(GridPlot,trim(ModelOutputFile))
			CALL AWE_saveCanvas(GridPlot,trim(CanvasFileName))
			if(isteps.ne.stepend)then			! don't erase the last one
				CALL AWE_clearCanvas(GridPlot)
				endif
			!CALL AWE_CloseCanvas(GridPlot)
			end do

	case(11)
		call SetUpGridPlot2(GridPlot,minX,maxX,minY,maxY)
	case(12)
!		call SetUpGridPlot3(GridPlot,minX,maxX,minY,maxY)
		call SetUpGridPlot3(GridPlot)
	case(2)
		call PlotGrid(GridPlot,1)		! this plots nodes only
	case(3)
		call PlotPoints(GridPlot,1,0,0)		! one means pick color; 0 = ask to connect the points
	case(4)
		call LabelNodes(GridPlot)
	case(5)
		call LabelSegs(GridPlot)
	case(6)
		call DrawCrystals2(GridPlot)
	case(7)
		call LabelCrystals(GridPlot)
	case(8)
		call Plot_GB_Comp(GridPlot,0,0,0.,0.)
	case(80)
		call Plot_GB_Comp_2(GridPlot,0,0,0.,0.)
	case(9)
		call Plot_Mineral_Moles_Scaled(GridPlot,0,0,0,0)
	case(13)
		call PlotAffinityGrid(GridPlot)
	case(14)
		call GarnetZoningMap(GridPlot)

	case(17)
170		continue
		write(*,*)'Input Seg to move'
		read(*,*)iSeg
		if(iSeg.eq.0)go to 1
		write(*,*)' Input mole scale factor (around 10-100-1000)'
		read(*,*)moleScaleFactor
		debug = 1
		call GrowSegs(GridPlot,iSeg,moleScaleFactor,debug)
		go to 170
	case(18)	! List segment information
		call ListSegInfo()


	case(20)
		call GibbsBegin()
	case(21)
		call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the node and seg XY coordinates
		close(40)
	case(22)
		call MakeAnimation()
	case(23)
		write(*,*)'Input logical unit number'
		read(*,*)iunit
		call DumpGrid(iunit)
		
	case(30)
! 		write(*,*)' Input crystal number to scale'
! 		read(*,*)iXl
! 		if(iXl.eq.0)go to 1
! 		write(*,*)' Input scale factor'
! 		read(*,*)scaleFact
! 		call ScaleCrystalNodes(iXl,scaleFact)

!		call ScaleCrystals(0,0,0.)
	case(31)
! 		write(*,*)' Input crystal number to check segs'
! 		read(*,*)iXl
! 		if(iXl.eq.0)go to 1
! 		call CheckCrystalSegs(iXL)
	case(32)
		write(*,*)' Do this routine only after opening the initial model file'
		write(*,*)' Note that areas are accurate for any model step, but the phase compositions do not account for zoning'
		write(*,*)' do you want to continue? 0 = no, 1 = yes '
		read(*,*)i
		if(i.eq.0)go to 1
		call GridBulkComposition()
		write(*,*)' Note that areas are accurate for any model step, but the phase compositions do not account for zoning'
		write(*,*)'     so the bulk composition is only valid for the initial model step'
		pause ' Output is in output window ... hit return to continue'
	case(33)
! 		call CrystalMolesTotal()
! 		write(*,*)'Input logical unit number'
! 		read(*,*)iunit
! 		call WriteCrystalMoles(iunit)
! 		close(iunit)
	case(34)
! 34		continue
! 		do i = 1,numCrystalsToScale
! 			write(*,*)i,CrystalsToScale(i),CrystalScaleFact(CrystalsToScale(i))
! 			end do
! 		write(*,*)' Input the one to scale and ∆moles'
! 		read(*,*)iXl,deltaMoles
! 		if(i.eq.0)go to 1
! 		write(*,*)' Crystal number and area ',iXL,CrystalArea(iXl)
! 		Call MolesToAreaRoutine(iXL,deltaMoles)
! 		write(*,*)' Crystal number and area ',iXL,CrystalArea(iXl)
! 		go to 34
	case(35)

		call CheckNodeInformation()


	case(36)
		write(*,*)' Routine to calculate Crystal (polygon) areas for any grid configuration'
		call PolygonArea()
		pause ' hit return to continue'

	case(37)
37		continue
		write(*,*)'------------------'
		write(*,*)'Input crystal you wish to check (0 to exit)'
		read(*,*)PCN
		if(PCN.eq.0)go to 1
		write(*,*)'CrystalNode  iNode  iSeg'
		do k = 1,numCrystalNodes(PCN) 
			iSeg = CrystalSegs(PCN,k)
			iNode = CrystalNodes(PCN,k)
			write(*,*)k,iNode,iSeg
			end do
		go to 37
		
	case(50)
		Call Muscovite_Biotite()
	case(200)
		call SaveIllustratorFile()
	case default
		call FSS_Alert('Alert', 'You chose poorly')
	end select
	go to 1
	end
	

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine SetUpGridPlot(GridPlot,maxX,maxY)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
!	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	INCLUDE "GibbsFiles.inc"				
!	real*4 asprat
	real*8 maxX,maxY

	close(PSUNIT)         			! close old PostScript scratch file
	call psopcl(4)    			! open a new scratch file for the Illustrator output


	CurrentColor = 1		! black is the default
	if(maxY.gt.maxX)then
		xmax = maxY		! these are set to the same values so the x and y axes have equal values and the grid isn't distorted
		ymax = maxY
		else
		xmax = maxX
		ymax = maxX
		endif
	xmin = 0.
	ymin = 0.
! 	xlen = 30		! values are in cm
! 	ylen = 30
! 	write(*,*)'Input value for length of plot in cm (e.g. 30-50 or greater)'
! 	read(*,*)xlen
	xlen=50
	ylen = xlen		! the x and y dimensions must be the same
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
! 	GridPlot%width  = 1.2*(xmax)*xconv			! *xScale scales the X dimension
! 	GridPlot%height = 1.2*(ymax)*yconv
	GridPlot%width  = 1.2*maxX*xconv			! *xScale scales the X dimension
	GridPlot%height = 1.2*maxY*yconv
	xor = 20.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = GridPlot%height - 50
!	yor =  1.2*(ymax)*yconv	- 20		! see if this works better
	write(*,*)GridPlot%width,GridPlot%height
	write(*,*)xor,yor

	GridPlot%title = ModelOutputFile
	CALL AWE_createCanvas(GridPlot)	


	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine SetUpGridPlot10(GridPlot,maxX,maxY,plotWidth)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
!	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	INCLUDE "GibbsFiles.inc"				
!	real*4 asprat
	real*8 maxX,maxY
	real*4 plotWidth


	CurrentColor = 1		! black is the default
	if(maxY.gt.maxX)then
		xmax = maxY		! these are set to the same values so the x and y axes have equal values and the grid isn't distorted
		ymax = maxY
		else
		xmax = maxX
		ymax = maxX
		endif
	xmin = 0.
	ymin = 0.
! 	xlen = 30		! values are in cm
! 	ylen = 30
	xlen = plotWidth
!	write(*,*)'Input value for length of plot in cm (e.g. 30 or greater)'
!	read(*,*)xlen
	ylen = xlen		! the x and y dimensions must be the same
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
! 	GridPlot%width  = 1.2*(xmax)*xconv			! *xScale scales the X dimension
! 	GridPlot%height = 1.2*(ymax)*yconv
	GridPlot%width  = 1.2*maxX*xconv			! *xScale scales the X dimension
	GridPlot%height = 1.2*maxY*yconv
	xor = 20.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = GridPlot%height - 50
!	yor =  1.2*(ymax)*yconv	- 20		! see if this works better
	write(*,*)GridPlot%width,GridPlot%height
	write(*,*)xor,yor

	GridPlot%title = ModelOutputFile
	CALL AWE_createCanvas(GridPlot)	

	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine SetUpGridPlot2(GridPlot,minX,maxX,minY,maxY)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
!	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
!	real*4 asprat
	real*8 minX,maxX,minY,maxY
	integer*8 iok
	CurrentColor = 1		! black is the default


	write(*,*)'minX, maxX',minX,maxX
	write(*,*)'minY, maxY',minY,maxY
! 	write(*,*)' Node to check XY?'
! 	read(*,*)iNode
! 	write(*,*)nodeX(iNode),nodeY(iNode)
	pause 'Ready? Hit return to set up plot'

	close(PSUNIT)         			! close old PostScript scratch file
	call psopcl(4)    			! open a new scratch file for the Illustrator output


    	Call Setplot(iok)

! 	if(maxY.gt.maxX)then
! 		xmax = maxY
! 		ymax = maxY
! 		else
! 		xmax = maxX
! 		ymax = maxX
! 		endif
! 	xmin = 0.
! 	ymin = 0.
	xlen = 30		! values are in cm
	ylen = 30
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	GridPlot%width  = 1.2*(xmax-xmin)*xconv			! *xScale scales the X dimension
	GridPlot%height = 1.2*(ymax-ymin)*yconv
! 	GridPlot%width  = 1.2*maxX*xconv			! *xScale scales the X dimension
! 	GridPlot%height = 1.1*maxY*yconv
	xor = 20.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = GridPlot%height - 20
	yor =  1.2*(ymax-ymin)*yconv	- 20		! see if this works better

	CALL AWE_createCanvas(GridPlot)	

	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine SetUpGridPlot3(GridPlot)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	INCLUDE "GibbsFiles.inc"				
	real*4 offset
	integer*4 iNode
	CurrentColor = 1		! black is the default


 	write(*,*)' Input Node to plot?'
1 	read(*,*)iNode
	if(iNode.eq.0)return
	write(*,*)'minX, maxX',minX,maxX
	write(*,*)'minY, maxY',minY,maxY
 	write(*,*)'node X,Y  ',nodeX(iNode),nodeY(iNode)
	write(*,*)'Input value to offset from node (e.g. 50-100) -- input 0 to abort'
	read(*,*)offset
	if(offset.eq.0)return
!	pause 'Ready? Hit return to set up plot'

!    	Call Setplot(iok)
	close(PSUNIT)         			! close old PostScript scratch file
	call psopcl(4)    			! open a new scratch file for the Illustrator output

	

	xmin = nodeX(iNode) - offset
	xmax = nodeX(iNode) + offset
	ymin = nodeY(iNode) - offset
	ymax = nodeY(iNode) + offset

	xlen = 30		! values are in cm
	ylen = 30
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
	GridPlot%width  = 1.2*(xmax-xmin)*xconv			! *xScale scales the X dimension
	GridPlot%height = 1.2*(ymax-ymin)*yconv
! 	GridPlot%width  = 1.2*maxX*xconv			! *xScale scales the X dimension
! 	GridPlot%height = 1.1*maxY*yconv
	xor = 20.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = GridPlot%height - 20
	yor =  1.2*(ymax-ymin)*yconv	- 20		! see if this works better

	GridPlot%title = ModelOutputFile
	CALL AWE_createCanvas(GridPlot)	

	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine PlotGrid(GridPlot,pickColor)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	real*4 xnode,ynode,x,y,dX,dY
	integer*4 i,j,index,iup,pickColor,connectThePoints
	CurrentColor = 1		! black is the default
1	continue
	if(pickColor.eq.1)then
		Write(*,*)'Select color for grid'
		write(*,*)' 1 = black'
		write(*,*)' 2 = red'
		write(*,*)' 3 = green'
		write(*,*)' 4 = blue'
		read(*,*)icolor
		if(icolor.eq.0)return
		select case(icolor)
			case(1)
				CurrentColor = 1		! Colors from Input problem file = colors for elements
			case(2)
				CurrentColor = 3
			case(3)
				CurrentColor = 8
			case(4)
				CurrentColor = 5
			case default
				call fss_alert('ALERT',' You chose poorly')
				go to 1
			end select

		write(*,*)' Connect the points? 0 = no, 1 = yes'
		read(*,*)connectThePoints

		endif

	pen%penStyle = CanvasPenStyle_SolidLine
	Do index = 1,numNodes
		x = nodeX(index)
		y = nodeY(index)
		dx = .005*(xmax-xmin)
		dy = dx
		!Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
	!	write(*,*)' '
		if(connectThePoints.eq.1)then
			do j = 1, numNodeSegs(index)
			!	write(*,*)index,j
				iup = 0
				xnode = nodeX(index) 		
				ynode = nodeY(index)
				call plot(GridPlot,xnode,ynode,iup)
				iup = 1	
				i = nodeNodeConnect(index,j)
				if(i.lt.index)cycle		! skip if connection is to a node with a lower index
				x = nodeX(i)
				y = nodeY(i)
				call plot(GridPlot,x,y,iup)
				end do
			else
			do j = 1, numNodeSegs(index)
				xnode = nodeX(index) 		
				ynode = nodeY(index)
				Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
				end do
			
			endif	
		end do
	CurrentColor = 1	! back to black
	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine PlotPoints(GridPlot,pickColor,connectPointsSwitch,drawPointsSwitch)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	real*4 x,y,dX,dY,radius
	integer*4 j,iSeg,pickColor,connectThePoints,iup,connectPointsSwitch,drawPointsSwitch,drawThePoints

	CurrentColor = 1		! Black is the default
1	continue
	if(pickColor.eq.1)then
		Write(*,*)'Select color for point ellipses'
		write(*,*)' 1 = black'
		write(*,*)' 2 = red'
		write(*,*)' 3 = green'
		write(*,*)' 4 = blue'
		read(*,*)icolor
		select case(icolor)
			case(1)
				pen%penColor = AWE_black
				brush%brushColor = AWE_black
				CurrentColor = 1		! Colors from Input problem file = colors for elements
			case(2)
				pen%penColor = AWE_red
				brush%brushColor = AWE_red
				CurrentColor = 3
			case(3)
				pen%penColor = AWE_green
				brush%brushColor = AWE_green
				CurrentColor = 8
			case(4)
				pen%penColor = AWE_blue
				brush%brushColor = AWE_blue
				CurrentColor = 5
			case default
				call fss_alert('ALERT',' You chose poorly')
				go to 1
			end select
		endif
	if(connectPointsSwitch.eq.0)then
		write(*,*)' Connect the points? 0 = no, 1 = yes'
		read(*,*)connectThePoints
		else
		connectThePoints = 1
		endif
	if(drawPointsSwitch.eq.0)then
		write(*,*)' Draw the points? 0 = no, 1 = yes'
		read(*,*)drawThePoints
		else
		drawThePoints = 1
		endif
	pen%penStyle = CanvasPenStyle_SolidLine

! 	dx = .002*(xmax-xmin)	! sets the size of the ellipse 
	dx = 0.5	! keep the points the same size
	dy = dx
	radius = dx
	if(drawThePoints.eq.1)then
		Do iSeg = 1,numSegs
			!do j = 1,numPoints(i)
			do j = numPointStart(iSeg),numPointEnd(iSeg)
				x = pointX(iSeg,j)
				y = pointY(iSeg,j)
				Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
				!call psCircle(x,y,radius,StrokeFillBoth)	!SFB = 1 stroke only; =2 Fill only; =3 both
				call psCircle(x,y,radius,3)	!SFB = 1 stroke only; =2 Fill only; =3 both
				end do
			end do
		endif
	if(connectThePoints.eq.1)then
		Do iSeg = 1,numSegs
			iup = 0
			do j = numPointStart(iSeg),numPointEnd(iSeg)
				x = pointX(iSeg,j)
				y = pointY(iSeg,j)
				call plot(GridPlot,x,y,iup)
				call pPlot(x,y,iup)
				iup = 1
				end do
			call pPenup
			end do
		endif
	CurrentColor = 1		! Black is the default
	
	
	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine LabelNodes(GridPlot)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	character*5 text
	real*4 x,y
	integer*4 index,textsize,nodeToPlot2

1	continue
	write(*,*)'Number of nodes = ',numNodes
	write(*,*)' Input node(s) to label'
	write(*,*)'  0 = exit'
	write(*,*)'  n = node number'
	write(*,*)' -1 = all nodes'
	read(*,*)nodeToPlot2
	if(nodeToPlot2.gt.numNodes)then
		call FSS_Alert('Alert','Node number is out of range')
		go to 1
		endif
	if(nodeToPlot2.eq.0)return

	CurrentColor = 1		! black is the default
	textsize = 12
	pen%penStyle = CanvasPenStyle_SolidLine
101	format(I4)
	if(nodeToPlot2.eq.-1)then
		Do index = 1,numNodes
		x = nodeX(index) + .005*(xmax-xmin)
		y = nodeY(index) + .005*(ymax-ymin)
		write(text,101)index
		text = adjustL(text)
		CALL TextOnPlot(GridPlot,x,y,text,textSize)
	!	SUBROUTINE TextOnPlot(GridPlot,xplot,yplot,text,textSize)
		end do
		else
		index = nodeToPlot2
		x = nodeX(index) + .01*(xmax-xmin)
		y = nodeY(index) + .01*(ymax-ymin)
		write(text,101)index
		text = adjustL(text)
		CALL TextOnPlot(GridPlot,x,y,text,textSize)
		endif
	go to 1	

	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine LabelSegs(GridPlot)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	character*5 text
	real*4 x,y
	integer*4 j,textsize,segToPlot2

1	continue
	write(*,*)'Number of segs = ',numSegs
	write(*,*)' Input seg(s) to label'
	write(*,*)'  0 = exit'
	write(*,*)'  n = seg number'
	write(*,*)' -1 = all segs'
	read(*,*)segToPlot2
	if(segToPlot2.gt.numSegs)then
		call FSS_Alert('Alert','Seg number is out of range')
		go to 1
		endif
	if(segToPlot2.eq.0)return

	CurrentColor = 1		! black is the default

	textsize = 12
	pen%penStyle = CanvasPenStyle_SolidLine
	101	format(I4)
	if(segToPlot2.eq.-1)then
		Do j = 1,numSegs
			!	ipoint = numPoints(j)/2 + 1
			!	x = pointX(j,ipoint) + .01*(xmax-xmin)
			x = (segX(j,1) + segX(j,2))/2.
			y = (segY(j,1) + segY(j,2))/2.
			!	x = pointX(j,ipoint) 
			!	y = pointY(j,ipoint) + .01*(ymax-ymin)
			write(text,101)j
			text = 's'//adjustL(text)
			CALL TextOnPlot(GridPlot,x,y,text,textSize)
			!	SUBROUTINE TextOnPlot(GridPlot,xplot,yplot,text,textSize)
			end do
		else
		j = segToPlot2
		x = (segX(j,1) + segX(j,2))/2.
		y = (segY(j,1) + segY(j,2))/2.
		write(text,101)j
		text = 's'//adjustL(text)
		CALL TextOnPlot(GridPlot,x,y,text,textSize)
		endif
	go to 1	

	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine LabelCrystals(GridPlot)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	character*15 text
	real*4 x,y
	integer*4 index,textsize,whatToPlot
	CurrentColor = 1		! black is the default

	write(*,*)' Plot phase names or crystal index number?'
	write(*,*)' 1 = phase name, 2 = crystal index number, 3 = both'
	read(*,*)whatToPlot
	if(whatToPlot.eq.0)return
	textsize = 12
	pen%penStyle = CanvasPenStyle_SolidLine
101	format(I5)
	Do index = 1,numCrystals
	!	if(CrystalMIFID(index).eq.2)cycle		!don't label quartz	
		x = CrystalCenterX(index)
		y = CrystalCenterY(index)
		select case (whatToPlot)
			case(1)
				text = adjustL(CrystalPhaseName(index))
			case(2)
				write(text,101)index
				!write(*,*)index,text
			case(3)
				write(text,101)index
				text = trim(adjustL(text))//' '//adjustL(CrystalPhaseName(index))
			case default
			end select
	!	write(*,*)index,CrystalPhaseName(index),text
		text = adjustL(text)
		CALL TextOnPlot(GridPlot,x,y,text,textSize)
	!	SUBROUTINE TextOnPlot(GridPlot,xplot,yplot,text,textSize)
		end do

	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
! 	Subroutine MakeSegs()
!       	implicit none
! 	include "Diffuse_GB_Hex.inc"
! 	integer*4 i,j,index
! 	numSegs = 0
! 	Do 10 index = 1,numNodes
! !	write(*,*)' '
! 	do 12 j = 1, numNodeSegs(index)
! !	write(*,*)index,j
! 	i = iConnect(index,j)
! 	if(i.lt.index)go to 12		! skip if connection is to a node with a lower index
! 	numSegs = numSegs+1
! 	segNodes(numSegs,1) = index
! 	segNodes(numSegs,2) = i
! 	segLength(numSegs) = sqrt((nodeX(i)-nodeX(index))**2 + (nodeY(i)-nodeY(index))**2)
! 12	continue
! 10	continue
! 	return
! 	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine DrawCrystals(canvas)
!	routine to draw all of the model crystals on the grid
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	real*4 pointArray(400,2)
	integer*4 i,j,k,numPts,PCN

!123	Chlorite
!green                           32768    0  128    0   1.000   0.000   1.000   0.498
!1
!94,95,107,119,118,106

	Do i = 1,numPhases
		do j = 1,numPhaseCrystals(i)
			PCN = PhaseCrystalNumber(i,j)		! Phase Crystal Number
			do k = 1,numCrystalNodes(PCN) 
				pointArray(k,1) = nodeX(CrystalNodes(PCN,k))
				pointArray(k,2) = nodeY(CrystalNodes(PCN,k))
				end do
			numPts = numCrystalNodes(PCN)
			pen%penColor = PhaseColorNumber(i)
			pen%penStyle = CanvasPenStyle_SolidLine
			brush%brushColor = PhaseColorNumber(i)
			call DrawPolygonOnScreen(canvas,numPts,pointArray,brush,pen)
			end do
		end do
	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine DrawCrystals2(canvas)
!	routine to draw all of the model crystals on the grid
!	I need to define the polygon using all of the points along the segments
!	The other code only used the nodes, which only works for the initial grid
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	real*4 pointArray(400,2),CMYKTemp(4)
	integer*4 i,j,k,numPts,PCN,L,m,iSeg,iNode

!123	Chlorite
!green                           32768    0  128    0   1.000   0.000   1.000   0.498
!1
!94,95,107,119,118,106

	Do i = 1,numPhases
		do j = 1,numPhaseCrystals(i)			! this will loop on a specific phase type
			PCN = PhaseCrystalNumber(i,j)		! Phase Crystal Number
			numPts = 0
			do k = 1,numCrystalNodes(PCN) 
				iSeg = CrystalSegs(PCN,k)
				iNode = CrystalNodes(PCN,k)
				do L = 1,3
					if(iSeg.eq.nodeSegConnect(iNode,L))then		! this is the segment we want
!						if(nodeSegEnd(iNode,L).eq.1)then	! we start at the beginning of the segment
						if(segNodes(iSeg,1).eq.iNode)then
							do m = numPointStart(iSeg),numPointEnd(iSeg)
								numPts = numPts + 1
								pointArray(numPts,1) = pointX(iSeg,m)
								pointArray(numPts,2) = pointY(iSeg,m)
								end do
							else				! we start at the end of the segment
							do m = numPointEnd(iSeg),numPointStart(iSeg),-1
								numPts = numPts + 1
								pointArray(numPts,1) = pointX(iSeg,m)
								pointArray(numPts,2) = pointY(iSeg,m)
								end do
							endif															
						endif
					end do				
				end do
			pen%penColor = PhaseColorNumber(i)
			pen%penStyle = CanvasPenStyle_SolidLine
			brush%brushColor = PhaseColorNumber(i)
			call DrawPolygonOnScreen(canvas,numPts,pointArray,brush,pen)
			CMYKTemp(1) = PhaseCMYK(i,1)
			CMYKTemp(2) = PhaseCMYK(i,2)
			CMYKTemp(3) = PhaseCMYK(i,3)
			CMYKTemp(4) = PhaseCMYK(i,4)
			call DrawPSPolygon(numPts,pointArray,CMYKTemp)
			end do
		end do
	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine DrawPolygonOnScreen(canvas,numPts,pointArray,brush,pen)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
!	TYPE(AWE_Point) :: poly_points(:)
	TYPE(AWE_Point) :: poly_points(numPts)

	real xpix,ypix			! functions to convert X,Y to pixels
	real*4 pointArray(400,2)
	integer*4 numPts,i
	!	write(*,*)T,dT,P,dP

	do i = 1,numPts
		poly_points(i)%x = xpix(pointArray(i,1))
		poly_points(i)%y = ypix(pointArray(i,2))
		end do

	call AWE_canvasDrawPolygon(canvas,poly_points,pen,brush)	
	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine Plot_GB_Comp(GridPlot,autoPlot,iEl,minEl,maxEl)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	real*4 x,y,dX,dY,increment,range,comp,minEl,maxEl,autoPlot
	integer*4 i,iEl,iSeg,iPointX

	if(autoPlot.eq.1)go to 10
1	continue
	write(*,*)'Pick element to plot (0 to exit)'
	do i = 1,numEl
		write(*,*)i,elName(i)
		end do
	write(*,*)' Input element to plot'
	read(*,*)iEl
	if(iEl.eq.0)return
	!determine range of composition values
	maxEl = -1.0
	minEl = 1.0d5		! make very large to start
	do iSeg = 1,numSegs
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
			if(pointComp(iseg,iPointX,iEl).gt.maxEl)then
				maxEl = pointComp(iseg,iPointX,iEl)
				endif
			if(pointComp(iseg,iPointX,iEl).lt.minEl)then
				minEl = pointComp(iseg,iPointX,iEl)
				endif
			end do
		end do
	write(*,*)' Minimum = ',minEl
	write(*,*)' Maximum = ',maxEl

10	continue
	range = (maxEl-minEl)
	increment = numPositiveColors - 1
	dx = .002*(xmax-xmin)
	dy = dx
	do iSeg = 1,numSegs
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
			x = pointX(iSeg,iPointX)
			y = pointY(iSeg,iPointX)
			comp = pointComp(iseg,iPointX,iEl)	
			if(range.eq.0.0)then
				i = 1
				else
				i = 1 + int(increment*(comp - minEl)/range)			! fractional amount above minimum
				endif
			write(*,*)iSeg,iPointX,comp,i
			pen%penColor = positiveColors(i)
			pen%penStyle = CanvasPenStyle_SolidLine
			brush%brushColor = positiveColors(i)

			Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
			end do
		end do
	go to 1
	
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine Plot_Mineral_Moles(GridPlot)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	real*4 x,y,dX,dY,increment,rangePlus,rangeMinus,minMoles,maxMoles,plotMoles,shrink
	integer*4 i,iSeg,iPointX,iNode,iMin,k,iC

	write(*,*)'Input shrink amount (in fraction)'
	read(*,*)shrink
1	continue
	write(*,*)'Pick mineral to plot (0 to exit)'
	do i = 2,numPhMIF		! we 
		write(*,*)i,PhName(i)
		end do
	write(*,*)' Input mineral to plot'
	read(*,*)iMin
	if(iMin.eq.0)return
	!determine range of composition values
	maxMoles = -1.0e5
	minMoles = 1.0e5			

	do iNode = 1,numNodes
!	do j = 1,numNodesWith3Phases
!		iNode = nodesWith3PhasesIndex(j)
		select case(numNodeReactionPhases(iNode))
			case(3)
				do i = 1,3				! 3 phases at each node
					if(iMin.eq.nodeMIFID(iNode,i))then
						if(nodePhaseMoles(iNode,i).gt.maxMoles)then
							maxMoles = nodePhaseMoles(iNode,i)
							endif
						if(nodePhaseMoles(iNode,i).lt.minMoles)then
							minMoles = nodePhaseMoles(iNode,i)
							endif
						endif
					end do
			case(2)	
!				do j = 1,numNodesWith2Phases
!					iNode = nodesWith2PhasesIndex(j)
				do i = 1,2				! 3 phases at each node
					if(iMin.eq.nodeMIFID(iNode,i))then
						if(nodePhaseMoles(iNode,i).gt.maxMoles)then
							maxMoles = nodePhaseMoles(iNode,i)
							endif
						if(nodePhaseMoles(iNode,i).lt.minMoles)then
							minMoles = nodePhaseMoles(iNode,i)
							endif
						endif
					end do
			case default
			end select	
		end do
		
		
!	do j = 1,numSegsWith2Phases
	do iSeg = 1,numSegs
!		iSeg = segReactionPhases(j)
		do i = 1,2				! 2 phases at each point
			if(iMin.eq.segMIFID(iSeg,i))then
				do iPointX = numPointStart(iSeg),numPointEnd(iSeg)			
					if(pointPhaseMoles(iSeg,iPointX,i).gt.maxMoles)then
						maxMoles = pointPhaseMoles(iSeg,iPointX,i)
						endif
					if(pointPhaseMoles(iSeg,iPointX,i).lt.minMoles)then
						minMoles = pointPhaseMoles(iseg,iPointX,i)
						endif
					end do
				endif
			end do
		end do

	write(*,*)' Minimum = ',minMoles
	write(*,*)' Maximum = ',maxMoles
	increment = numPositiveColors - 1
	
	! we are using 0 moles as the min and max, respectively
	rangePlus = maxMoles
	rangeMinus = minMoles
	! prevent divide by zero
	if(rangePlus.eq.0.0)rangePlus = 1.e-20
	if(rangeMinus.eq.0.0)rangeMinus = 1.e-20
	dx = .002*(xmax-xmin)
	dy = dx

	do iNode = 1,numNodes

! 	do j = 1,numNodesWith3Phases
! 		iNode = nodesWith3PhasesIndex(j)
		select case(numNodeReactionPhases(iNode))
			case(3)
				do i = 1,3				! 3 phases at each node
					if(iMin.eq.nodeMIFID(iNode,i))then
						x = nodeX(iNode)
						y = nodeY(iNode)
						iC = nodeCrystalIndex(iNode,i)
						x = x + shrink*(CrystalCenterX(iC) - x)
						y = y + shrink*(CrystalCenterY(iC) - y)
						plotMoles = nodePhaseMoles(iNode,i)
						if(plotMoles.gt.0)then
							k = 1 + int(increment*(plotMoles - 0.)/rangePlus)			! fractional amount above minimum
							pen%penColor = positiveColors(k)
							pen%penStyle = CanvasPenStyle_SolidLine
							brush%brushColor = positiveColors(k)
							else
							k = 1 + int(increment*(plotMoles - 0.)/rangeMinus)			! fractional amount above minimum
							pen%penColor = negativeColors(k)
							pen%penStyle = CanvasPenStyle_SolidLine
							brush%brushColor = negativeColors(k)
							endif
						Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
						endif
					end do
			case(2)	

! 	do j = 1,numNodesWith2Phases
! 		iNode = nodesWith2PhasesIndex(j)
				do i = 1,2			
					if(iMin.eq.nodeMIFID(iNode,i))then
						x = nodeX(iNode)
						y = nodeY(iNode)
						iC = nodeCrystalIndex(iNode,i)
						x = x + shrink*(CrystalCenterX(iC) - x)
						y = y + shrink*(CrystalCenterY(iC) - y)
						plotMoles = nodePhaseMoles(iNode,i)
						if(plotMoles.gt.0)then
							k = 1 + int(increment*(plotMoles - 0.)/rangePlus)			! fractional amount above minimum
							pen%penColor = positiveColors(k)
							pen%penStyle = CanvasPenStyle_SolidLine
							brush%brushColor = positiveColors(k)
							else
							k = 1 + int(increment*(plotMoles - 0.)/rangeMinus)			! fractional amount above minimum
							pen%penColor = negativeColors(k)
							pen%penStyle = CanvasPenStyle_SolidLine
							brush%brushColor = negativeColors(k)
							endif
						Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
						endif
					end do
			case default
			end select	
		end do
		
!	do j = 1,numSegsWith2Phases
	do iSeg = 1,numSegs
!		iSeg = segReactionPhases(j)
		do i = 1,2				! 2 phases at each point
			if(iMin.eq.segMIFID(iSeg,i))then
				iC = segCrystalID(iSeg,i)		! this should be the same crystal as with the node???
				do iPointX = numPointStart(iSeg),numPointEnd(iSeg)			
					x = pointX(iSeg,iPointX)
					y = pointY(iSeg,iPointX)
					x = x + shrink*(CrystalCenterX(iC) - x)
					y = y + shrink*(CrystalCenterY(iC) - y)
					plotMoles = pointPhaseMoles(iseg,iPointX,i)
					if(plotMoles.gt.0)then
						k = 1 + int(increment*(plotMoles - 0.)/rangePlus)			! fractional amount above minimum
						pen%penColor = positiveColors(k)
						pen%penStyle = CanvasPenStyle_SolidLine
						brush%brushColor = positiveColors(k)
						else
						k = 1 + int(increment*(plotMoles - 0.)/rangeMinus)			! fractional amount above minimum
						pen%penColor = negativeColors(k)
						pen%penStyle = CanvasPenStyle_SolidLine
						brush%brushColor = negativeColors(k)
						endif
					Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
					end do
				endif
			end do
		end do


	go to 1
	
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine Plot_Mineral_Moles_Scaled(GridPlot,autoPlot,iMin,minMoles,maxMoles)
	! This routine plots circles on the grid with colors scaled to the total amount of moles of a phase that are produced (positive colors)
	! or consumed (negative colors)
	! The plot is on the grid
	! It uses nodePhaseMoles and segPhaseMoles so that the total number of moles are plotted at a given model cycle.
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	real*4 x,y,dX,dY,increment,rangePlus,rangeMinus,minMoles,maxMoles,plotMoles,shrink,PM
	integer*4 i,iSeg,iPointX,iNode,iMin,k,iC,maxNode,minNode,maxSeg,minSeg,autoPlot

	! if iMin = 0 then we need to choose the mineral to plot
	! if iMin > 0 then the mineral to plot is specified in the subroutine arguments

	if(autoPlot.eq.1)go to 1		! skip the min,max stuff

!	determine max and min of dataset


	!determine range of composition values
	maxMoles = -1.0e5
	minMoles = 1.0e5			
	do iNode = 1,numNodes

! 	do j = 1,numNodesWith3Phases
! 		iNode = nodesWith3PhasesIndex(j)
		select case(numNodeReactionPhases(iNode))
			case(3)
				do i = 1,3				! 3 phases at each node
		!		if(iMin.eq.nodeMIFID(iNode,i))then
					!PM = nodePhaseMoles(iNode,i)*numOxygens(nodeMIFID(iNode,i))
					PM = nodePhaseMoles(iNode,i)
					if(PM.gt.maxMoles)then
						maxMoles = PM
						maxNode = iNode
						endif
					if(PM.lt.minMoles)then
						minMoles = PM
						minNode = iNode
						endif
		!			endif
					end do
			case(2)	
!				do j = 1,numNodesWith2Phases
!				iNode = nodesWith2PhasesIndex(j)
				do i = 1,2				! 3 phases at each node
		!		if(iMin.eq.nodeMIFID(iNode,i))then
					!PM = nodePhaseMoles(iNode,i)*numOxygens(nodeMIFID(iNode,i))
					PM = nodePhaseMoles(iNode,i)
					if(PM.gt.maxMoles)then
						maxMoles = PM
						maxNode = iNode
						endif
					if(PM.lt.minMoles)then
						minMoles = PM
						minNode = iNode
						endif
		!			endif
					end do
			case default
			end select	
		end do
		
	write(*,*)' Minimum node = ',minMoles,minNode
	write(*,*)' Maximum node = ',maxMoles,maxNode

	maxSeg = 0
	minSeg = 0		
!	do j = 1,numSegsWith2Phases
	do iSeg = 1,numSegs
!		iSeg = segReactionPhases(j)
		do i = 1,2				! 2 phases at each point
!		if(iMin.eq.segMIFID(iSeg,i))then
			do iPointX = numPointStart(iSeg)+1,numPointEnd(iSeg)-1		
				!PM = pointPhaseMoles(iSeg,iPointX,i)*numOxygens(segMIFID(iSeg,i))
				PM = pointPhaseMoles(iSeg,iPointX,i)
				if(PM.gt.maxMoles)then
					maxMoles = PM
					maxSeg = iSeg
					endif
				if(PM.lt.minMoles)then
					minMoles = PM
					minSeg = iSeg
					endif
				end do
!			endif
			end do
		end do
	write(*,*)' Minimum seg = ',minMoles,minSeg
	write(*,*)' Maximum seg = ',maxMoles,maxSeg


	write(*,*)' Minimum = ',minMoles
	write(*,*)' Maximum = ',maxMoles
	write(*,*)' Input values for minimum (negative) and maximum (positive) moles'
	write(*,*)'Minimum = '
	read(*,*)minMoles
	write(*,*)'Maximum = '
	read(*,*)maxMoles
		
	
1	continue
	! we are using 0 moles as the min and max, respectively
	rangePlus = maxMoles
	rangeMinus = minMoles
	! prevent divide by zero
	if(rangePlus.eq.0.0)rangePlus = 1.e-20
	if(rangeMinus.eq.0.0)rangeMinus = 1.e-20

!	write(*,*)'Input shrink amount (in fraction)'
!	read(*,*)shrink
	shrink = 0.05
	increment = numPositiveColors - 1
	dx = .002*(xmax-xmin)
	dy = dx

	if(autoPlot.eq.0)then		! if not, then we already specified mineral in sub. arguments
		write(*,*)'Pick mineral to plot (0 to exit)'
		do i = 2,numPhMIF		! we 
			write(*,*)i,PhName(i)
			end do
		write(*,*)' Input mineral to plot'
		read(*,*)iMin
		if(iMin.eq.0)return
		endif
		
	do iNode = 1,numNodes

! 	do j = 1,numNodesWith3Phases
! 		iNode = nodesWith3PhasesIndex(j)
		select case(numNodeReactionPhases(iNode))
		case(3)
			do i = 1,3				! 3 phases at each node
				if(iMin.eq.nodeMIFID(iNode,i))then
					x = nodeX(iNode)
					y = nodeY(iNode)
					iC = nodeCrystalIndex(iNode,i)
					x = x + shrink*(CrystalCenterX(iC) - x)
					y = y + shrink*(CrystalCenterY(iC) - y)
					!plotMoles = nodePhaseMoles(iNode,i)*numOxygens(nodeMIFID(iNode,i))
					plotMoles = nodePhaseMoles(iNode,i)
					if(plotMoles.gt.0)then
						k = 2 + int(increment*(plotMoles - 0.)/rangePlus)			! fractional amount above minimum
						pen%penColor = positiveColors(k)
						pen%penStyle = CanvasPenStyle_SolidLine
						brush%brushColor = positiveColors(k)
						else
						k = 2 + int(increment*(plotMoles - 0.)/rangeMinus)			! fractional amount above minimum
						pen%penColor = negativeColors(k)
						pen%penStyle = CanvasPenStyle_SolidLine
						brush%brushColor = negativeColors(k)
						endif
					Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
					endif
				end do
		case(2)
! 	do j = 1,numNodesWith2Phases
! 		iNode = nodesWith2PhasesIndex(j)
			do i = 1,2				! 2 phases at each node
				if(iMin.eq.nodeMIFID(iNode,i))then
					x = nodeX(iNode)
					y = nodeY(iNode)
					iC = nodeCrystalIndex(iNode,i)
					x = x + shrink*(CrystalCenterX(iC) - x)
					y = y + shrink*(CrystalCenterY(iC) - y)
					!plotMoles = nodePhaseMoles(iNode,i)*numOxygens(nodeMIFID(iNode,i))
					plotMoles = nodePhaseMoles(iNode,i)
					if(plotMoles.gt.0)then
						k = 1 + int(increment*(plotMoles - 0.)/rangePlus)			! fractional amount above minimum
						pen%penColor = positiveColors(k)
						pen%penStyle = CanvasPenStyle_SolidLine
						brush%brushColor = positiveColors(k)
						else
						k = 1 + int(increment*(plotMoles - 0.)/rangeMinus)			! fractional amount above minimum
						pen%penColor = negativeColors(k)
						pen%penStyle = CanvasPenStyle_SolidLine
						brush%brushColor = negativeColors(k)
						endif
					Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
					endif
				end do
		case default
		end select
		end do
		
!	do j = 1,numSegsWith2Phases
	do iSeg = 1,numSegs
!		iSeg = segReactionPhases(j)
		do i = 1,2				! 2 phases at each point
			if(iMin.eq.segMIFID(iSeg,i))then
				iC = segCrystalID(iSeg,i)		! this should be the same crystal as with the node???
				do iPointX = numPointStart(iSeg)+1,numPointEnd(iSeg)-1			
					x = pointX(iSeg,iPointX)
					y = pointY(iSeg,iPointX)
					x = x + shrink*(CrystalCenterX(iC) - x)
					y = y + shrink*(CrystalCenterY(iC) - y)
					!plotMoles = pointPhaseMoles(iSeg,iPointX,i)*numOxygens(segMIFID(iSeg,i))
					plotMoles = pointPhaseMoles(iSeg,iPointX,i)
					if(plotMoles.gt.0)then
						k = 1 + int(increment*(plotMoles - 0.)/rangePlus)			! fractional amount above minimum
						pen%penColor = positiveColors(k)
						pen%penStyle = CanvasPenStyle_SolidLine
						brush%brushColor = positiveColors(k)
						else
						k = 1 + int(increment*(plotMoles - 0.)/rangeMinus)			! fractional amount above minimum
						pen%penColor = negativeColors(k)
						pen%penStyle = CanvasPenStyle_SolidLine
						brush%brushColor = negativeColors(k)
						endif
					Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
					end do
				endif
			end do
		end do

	if(autoPlot.eq.1)then
		return
		else
		go to 1
		endif	
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine SaveIllustratorFile()
	implicit none
	INCLUDE "PlotStuff.inc"				
	integer*4 i
	write(*,*)' Save plot as an Illustrator file?'
	write(*,*)' 0 = no'
	write(*,*)' 1 = yes'
	read(*,*)i
	if(i.eq.1)then
!		PSFileName = Trim(filein)//'.ai'        	! default ps file name
!	write(*,*)filein
!	write(*,*)psfilename
!	pause ' hit return'
		call psopcl(5)    				! save old PostScript scratch file
		endif
	return
	end
	

		
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine Muscovite_Biotite()
!	USE AWE_INTERFACES
	implicit none
!	TYPE(AWE_Canvas) :: GridPlot
	include "Diffuse_GB_Hex.inc"
	include "Assemb.inc"
	integer*4 status,i,j,kk,k,iNode,doit,MsMIFID,BtMIFID,MIFID,iMIFID,iSeg,iPointX
	integer*4 nodeToCheck(100),numNodesToCheck
	character*5 title
	real*8 MuscoviteMoles,BiotiteMoles,MuscovitePercent,MIFmoles(10),SEGmoles(10)
		
	write(*,*)'----------------------------'		
	write(*,*)' Routine to calculate muscovite and biotite area percents at nodes'
	write(*,*)' Open text file with nodes'
      	open(23,file='',status='old')
	read(23,*)title
	i = 1
	do
		read(23,*,end=5)nodeToCheck(i)
		i = i + 1
		end do
5	continue
	close(23)
	NumNodesToCheck = i - 1

10	continue

	write(*,*)'  0 = exit'
	write(*,*)'  1 = NODES: Calculate from the list of nodes '
	write(*,*)'  2 = NODES: Calculate from the list of nodes but only for garnet + biotite + muscovite '
	write(*,*)' 22 = NODES: Calculate from the list of nodes but only for nodes with garnet (±biotite ±muscovite) '
	write(*,*)'  3 = NODES: Calculate all phases from the entire grid'
	write(*,*)'  4 = SEGS:  Calculate all phases from the entire grid'
	read(*,*)doit

	select case (doit)
	case (0)
		return
	case(1,2,22)	 ! calculate from the list of nodes

	write(*,*)'Open model file'
	call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the node and seg XY coordinates
	close(40)
	write(*,*)'-------------------------------'
	write(*,*)' TotalCycles =  ',totalCycles
	write(*,*)'-------------------------------'

	! calculate biotite MIFID = 6
	BtMIFID = 6
	BiotiteMoles = 0.0d0
	write(*,*)' ---do biotite---------------------'
	
	do i = 1,numNodesToCheck
		iNode = nodeToCheck(i)
		if(doit.eq.2)then
			if(numNodeReactionPhases(iNode).lt.3)cycle	! only do nodes with all 3 phases
			endif
		if(doit.eq.22)then
			do j = 1,numNodeReactionPhases(iNode)
				if(nodeMIFID(iNode,j).eq.7)go to 221	! MIFID=7 is garnet this one is OK
				end do
			! if we get here then this node doesn't have a garnet
			cycle
			endif
221		continue	! this node is OK to do
		!write(*,*)iNode,nodeReactionPhases
		do kk = 1,numNodeReactionPhases(iNode)
			k = nodeReactionPhases(iNode,kk)
			do j = 1,3
				MIFID = nodeMIFID(iNode,j)
				if(MIFID.eq.BtMIFID)then
					! this is the one we want
					BiotiteMoles = BiotiteMoles + nodePhaseMoles(iNode,j)
					write(*,*)iNode,kk,k,j,MIFID,nodePhaseMoles(iNode,j)
					!go to the next node (don't count anything twice
					go to 20
					endif
				end do
			end do
		! if we get here then biotite is not one of the reaction phases
		write(*,*)iNode,kk,k,'  No biotite'
20		continue
		end do		! end checking every node for biotite
	write(*,*)'BiotiteMoles = ',BiotiteMoles

	! calculate Muscovite MIFID = 4
	MsMIFID = 4
	MuscoviteMoles = 0.0d0
	write(*,*)' ---do muscovite---------------------'
	do i = 1,numNodesToCheck
		iNode = nodeToCheck(i)
		if(doit.eq.2)then
			if(numNodeReactionPhases(iNode).lt.3)cycle	! only do nodes with all 3 phases
			endif
		if(doit.eq.22)then
			do j = 1,numNodeReactionPhases(iNode)
				if(nodeMIFID(iNode,j).eq.7)go to 222	! MIFID=7 is garnet this one is OK
				end do
			! if we get here then this node doesn't have a garnet
			cycle
			endif
222		continue	! this node is OK to do
		!write(*,*)iNode,nodeReactionPhases
		do kk = 1,numNodeReactionPhases(iNode)
			k = nodeReactionPhases(iNode,kk)
			do j = 1,3
				MIFID = nodeMIFID(iNode,j)
				if(MIFID.eq.MsMIFID)then
					! this is the one we want
					MuscoviteMoles = MuscoviteMoles + nodePhaseMoles(iNode,j)
					write(*,*)iNode,kk,k,j,MIFID,nodePhaseMoles(iNode,j)
					!go to the next node (don't count anything twice
					go to 25
					endif
				end do
			end do
		! if we get here then biotite is not one of the reaction phases
		write(*,*)iNode,kk,k,'  No muscovite'
25		continue
		end do		! end checking every node for muscovite

	write(*,*)'MuscoviteMoles = ',MuscoviteMoles
	write(*,*)' '
	MuscovitePercent = 100.0*MuscoviteMoles/(MuscoviteMoles + BiotiteMoles)
	write(*,*)'Muscovite %    = ',MuscovitePercent	
	write(*,*)' '
	write(*,*)' '
	go to 10

	case(3)	! calculate moles for the entire grid

	write(*,*)'Open model file'
	call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the node and seg XY coordinates
	close(40)
	write(*,*)'-------------------------------'
	write(*,*)' TotalCycles =  ',totalCycles
	write(*,*)'-------------------------------'

	Do iMIFID = 2,numPhMIF		! Phase 1 is always the grain boundary
		MIFmoles(iMIFID) = 0.0d0
		! calculate biotite MIFID = 6
	! 	BtMIFID = 6
	! 	BiotiteMoles = 0.0d0
	! 	write(*,*)' ---do biotite---------------------'
	
		do iNode = 1,numNodes
			! iNode = nodeToCheck(i)
			!write(*,*)iNode,nodeReactionPhases
			do kk = 1,numNodeReactionPhases(iNode)
				k = nodeReactionPhases(iNode,kk)
				do j = 1,3
					MIFID = nodeMIFID(iNode,j)
					if(MIFID.eq.iMIFID)then
						! this is the one we want
						MIFmoles(iMIFID) = MIFmoles(iMIFID) + nodePhaseMoles(iNode,j)
						!write(*,*)iNode,kk,k,j,MIFID,nodePhaseMoles(iNode,j)
						!go to the next node (don't count anything twice
						go to 30
						endif
					end do
				end do
			! if we get here then biotite is not one of the reaction phases
	! 		write(*,*)iNode,kk,k,'  No biotite'
	30		continue
			end do		! end checking every node for biotite

		write(*,*)iMIFID,phName(iMIFID),MIFmoles(iMIFID)
		end do

	go to 10

	case(4)		! calculate all phases in ALL SEGMENTS
	write(*,*)'Open model file'
	call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the node and seg XY coordinates
	close(40)
	write(*,*)'-------------------------------'
	write(*,*)' TotalCycles =  ',totalCycles
	write(*,*)'-------------------------------'

	Do iMIFID = 2,numPhMIF		! Phase 1 is always the grain boundary
		SEGmoles(iMIFID) = 0.0d0
	
		do iSeg = 1,numSegs
			do k = 1,numSegReactionPhases(iSeg)
				MIFID = SegMIFID(iSeg,k)
				if(MIFID.eq.iMIFID)then
					! this is the one we want
					do iPointX = numPointStart(iSeg)+1,numPointEnd(iSeg)-1  ! Note that pointStart and pointEnd are on the nodes
												! the values should be 0.0 but no point in adding zero to the sum
						SEGmoles(iMIFID) = SEGmoles(iMIFID) + pointPhaseMoles(iSeg,iPointX,k)
						end do
!						write(*,*)iSeg,kk,k,j,MIFID,nodePhaseMoles(iNode,j)
					!go to the next Seg (don't count anything twice)
					go to 40
					endif
				end do	! end loop on numReactionPhases

	40		continue
			end do		! end loop on numSegs


		write(*,*)iMIFID,phName(iMIFID),SEGmoles(iMIFID)
		end do		! end loop on numPhMIF

	go to 10




	case default
	end select
	go to 10

	end





