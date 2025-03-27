!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine GarnetZoningMap(GridPlot)
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
	INCLUDE "LUTColorStuff.inc"
	INCLUDE "Gibbsfiles.inc"
	integer*4 i,j,k,iPointX,iSeg,iGarnet,iup,izone,status,stepstart,stepend,whatToPlot,isteps
	real*8 minComp,maxComp
	real*4 x,y,dX,dY,increment,radius,penWidth,range,xCenter,yCenter,XYrange,GrtComp
	character*16 stepText

1	continue
	write(*,*)' Make your selection (0 to exit)(200 to save illustrator file)'
	write(*,*)' 0   = exit'
	write(*,*)' 1   = Plot zoning'
	write(*,*)' 2   = Adjust range'
	write(*,*)' 10  = ListLUT'
	write(*,*)' 100 = change pen width'
	write(*,*)' 200 = save Illustrator file'
	write(*,*)' 5   = Step through model steps'
	read(*,*)iZone
	
	select case (iZone)
	case(0)
		return
	case(2)
		write(*,*)' Specify range of element to map'
		write(*,*)' Minimum = '
		read(*,*)minComp
		write(*,*)' Maximum = '
		read(*,*)maxComp
		range = (maxComp-minComp)
		write(*,*)' Range = ',range

	case(10)
		write(12,*)''
		write(12,*)' i   RGB        CMYK'
		do i = 1,256
			write(12,*)i,LUT_RGB(i),(LUT_CMYK(i,j),j=1,4)
			end do
		write(12,*)''
		write(12,*)''
		go to 1

	case(100)
		write(*,*)'Input pen width (e.g., 0-10)'
		read(*,*)penWidth
		go to 1

	case(200)
		call SaveIllustratorFile()
		go to 1

	case(1,5)
	
	write(*,*)'Input coordinates for garnet center: X,Y'
	read(*,*)xCenter,yCenter
	write(*,*)'Input range of X and Y (limits to plot) (width of area to plot)'
	read(*,*)xyRange
	xmin = xcenter - xyrange/2.0
	xmax = xcenter + xyrange/2.0
	ymin = ycenter - xyrange/2.0
	ymax = ycenter + xyrange/2.0

	! set up grid plot canvas
	close(PSUNIT)         			! close old PostScript scratch file
	call psopcl(4)    			! open a new scratch file for the Illustrator output
	CurrentColor = 1		! black is the default
	xlen=30			! this sets the length (in cm) of the canvas
	ylen = xlen		! the x and y dimensions must be the same
      	CALL USER()			! sets the user coordinates that were input in Sub SetPlot
! 	GridPlot%width  = 1.2*(xmax)*xconv			! *xScale scales the X dimension
! 	GridPlot%height = 1.2*(ymax)*yconv
 	GridPlot%width  = 1.2*xyrange*xconv			! *xScale scales the X dimension
 	GridPlot%height = 1.2*xyrange*yconv
	xor = 20.				! origin for X-Y plot (xmin,ymin) in pixels
	yor = GridPlot%height - 50
!	yor =  1.2*(ymax)*yconv	- 20		! see if this works better
	write(*,*)'xCenter    ',xcenter
	write(*,*)'xMin,xMax  ',xmin,xmax
	write(*,*)'yCenter    ',ycenter
	write(*,*)'yMin,yMax  ',ymin,ymax
	write(*,*)GridPlot%width,GridPlot%height
	write(*,*)xor,yor
	GridPlot%title = 'Garnet zoning map'
	CALL AWE_createCanvas(GridPlot)	


!	Get the info for the output files to use for this plot
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
	read(*,*)iGarnet		! this is the number of the garnet in the MIF
	if(iGarnet.eq.0)return

	write(*,*)'Input starting step for making the plot (0 to begin at the beginning)'
	read(*,*)stepStart
	write(*,*)'Input ending step for making the plot '
	read(*,*)stepEnd


	write(*,*)' Select components to plot'
	write(*,*)' 0 = exit (done)'
	write(*,*)' 1 = Prp'
	write(*,*)' 2 = Alm'
	write(*,*)' 3 = Sps'
	write(*,*)' 4 = Grs'
	write(*,*)' 5 = Fe/Fe+Mg'
	read(*,*)whatToPlot
	if(whatToPlot.eq.0)go to 1

	! specify the limits of the compositional variable
	write(*,*)' specify the minimum and maximum range of the composition to map'
	read(*,*)minComp,maxComp
	range = maxComp-minComp

	pause 'Ready?? Hit return to execute '

	! These variables adjust the size of the symbol to plot
	increment = LUTnumColors - 1
	dx = .002*(xmax-xmin)
	dy = dx
	radius = dx
	! Write out the lines on the canvas
	pen%penStyle = CanvasPenStyle_SolidLine
	pen%penWidth = penWidth		! make the pen thicker
	
	
10	continue

	
	if(izone.eq.1)then

		! Read a model output file and process
		write(12,*)'**********************************************************'
		do isteps = stepStart,stepEnd
	
			write(12,*)'------------------------------'
			write(12,*)'Step number  ',isteps
			write(steptext,55)isteps		! this writes the step number into a character string for concatenating onto the file name
	55		format(I5)
			steptext = adjustL(steptext)
			ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
	!		if(isteps.eq.stepStart.or.isteps.eq.stepEnd)then
				!write(12,*)' Step number ',isteps
				write(12,*)'ModelOutputFile name ',ModelOutputFile		! only write out first and last file names
	!			endif
			open(40,FILE=ModelOutputFile,status = 'UNKNOWN')
	
			! Now read the model information
			isFileOpen = 1		! file is already open
	!		call OpenModelOutputFile_2(isFileOpen,1)		! This one reads affinities. 1 means read the node and seg XY coordinates
			call OpenModelOutputFile(isFileOpen,1)			! No affinities in this one. 1 means read the node and seg XY coordinates
	
			! scan the model information for the composition of garnet along each segment
			! Note that the segment endpoints align with nodes, so we don't need to tabulate the nodes
	
	
			Do iSeg = 1,numSegs
				iup = 0
				do k = 1,numSegReactionPhases(iSeg)
					if(segMIFID(iSeg,k).eq.iGarnet)then
	
						do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
							x = pointX(iSeg,iPointX)
							y = pointY(iSeg,iPointX)
							! check to see that these points are within bounds
							if(x.ge.xmin.and.x.le.xmax.and.y.ge.ymin.and.y.le.ymax)then		! this one is inside the search area
								! see if one of the phases is garnet
								! get the desired garnet composition
								GrtComp = pointPhaseComp(iSeg,iPointX,k,whatToPlot)
								write(12,*)iSeg,iPointX,x,y,k,GrtComp
								go to 14
								endif
							cycle	! This one is not inside the search area -- get the next point				
	
		14					continue
							if(GrtComp.gt.maxComp)then	! set the color to black
								i = 256
								else
								if(range.eq.0.0)then
									i = 1
									else
									i = 1 + int(increment*(GrtComp - minComp)/range)			! fractional amount above minimum
									endif
								endif
							pen%penColor = LUT_RGB(i)
							brush%brushColor = LUT_RGB(i)
							pen%penWidth = penWidth		! make the pen thicker
							call plot_2(GridPlot,x,y,iup,i)
							iup = 1
							end do		! end loop through points
						endif		! No garnet in this segment
					end do		! end loop on numReactionPhases
				end do  	! end loop on numSegs	
	
	
	! 		go to 22		! skip the next stuff
	
	
		
			! write out the lines in the Illustrator file
			! Each small segment between points needs to be a separate "line" so that each can be given a distinct color
			write(psunit,*)'u'	! this should "group" all of the line segments for later sizing
			Do iSeg = 1,numSegs
				iup = 0
				do k = 1,numSegReactionPhases(iSeg)
					if(segMIFID(iSeg,k).eq.iGarnet)then
						do iPointX = numPointStart(iSeg),numPointEnd(iSeg) - 1
							x = pointX(iSeg,iPointX)
							y = pointY(iSeg,iPointX)
							! check to see that these points are within bounds
							if(x.ge.xmin.and.x.le.xmax.and.y.ge.ymin.and.y.le.ymax)then		! this one is OK
								! get the desired garnet composition
								GrtComp = pointPhaseComp(iSeg,iPointX,k,whatToPlot)
								go to 16	! jump out and plot this one
								endif
							cycle
								
		16					continue
							if(GrtComp.gt.maxComp)then	! set the color to black
								i = 256
								else
								if(range.eq.0.0)then
									i = 1
									else
									i = 1 + int(increment*(GrtComp - minComp)/range)			! fractional amount above minimum
									endif
								endif		
							call pPlot_2(x,y,0,i)		! move to the start the line -- the color is set inside the routine
							x = pointX(iSeg,iPointX+1)
							y = pointY(iSeg,iPointX+1)
							call pPlot_2(x,y,1,i)		! draw the line to the next point
							call pPenup			! end the line
							iup = 1
							end do		! end loop through points
						endif			! no garnet in this segment
					end do			! end loop through numReactionPhases
				end do				! end loop on numSegs
	
			write(psunit,*)'U'		! Turn off "group" for this assemblage
		
		
			! this next loop puts circles on the plot -- it isn't really needed so skip
			go to 22		
			do iSeg = 1,numSegs
				do k = 1,numSegReactionPhases(iSeg)
					if(segMIFID(iSeg,k).eq.iGarnet)then
						do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
							x = pointX(iSeg,iPointX)
							y = pointY(iSeg,iPointX)
							! check to see that these points are within bounds
							if(x.ge.xmin.and.x.le.xmax.and.y.ge.ymin.and.y.le.ymax)then		! this one is OK
								! see if one of the phases is garnet
								! get the desired garnet composition
								GrtComp = pointPhaseComp(iSeg,iPointX,k,whatToPlot)
								go to 17	! jump out and plot this one
								endif
							cycle
			17				continue
		
							if(GrtComp.gt.maxComp)then	! set the color to black
								i = 256
								else
								if(range.eq.0.0)then
									i = 1
									else
									i = 1 + int(increment*(GrtComp - minComp)/range)			! fractional amount above minimum
									endif
								endif
							pen%penColor = LUT_RGB(i)
							brush%brushColor = LUT_RGB(i)
							pen%penStyle = CanvasPenStyle_SolidLine
							Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
				!			I need to set up the colors for this
							call psCircle_2(x,y,radius,3,i)	!SFB = 1 stroke only; =2 Fill only; =3 both
							end do		! end loop through points
						endif		! No garnet in this segment
					end do		! end loop on numReactionPhases
				end do  	! end loop on numSegs	
		
	
	22		continue
	
			end do		! end reading and processing each output file
	
		go to 1
	
! ----------------------------------------------------------------------	
		else		! step through model steps
	
		! Read a model output file and process
44		continue
		write(12,*)'**********************************************************'
		write(*,*)' Input step to examine (-1 to exit)'
		read(*,*)isteps
		if(isteps.lt.0)go to 1
!		do isteps = stepStart,stepEnd
	
			write(12,*)'------------------------------'
			write(12,*)'Step number  ',isteps
			write(steptext,55)isteps		! this writes the step number into a character string for concatenating onto the file name
!	55		format(I5)
			steptext = adjustL(steptext)
			ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
	!		if(isteps.eq.stepStart.or.isteps.eq.stepEnd)then
				!write(12,*)' Step number ',isteps
				write(12,*)'ModelOutputFile name ',ModelOutputFile		! only write out first and last file names
	!			endif
			open(40,FILE=ModelOutputFile,status = 'UNKNOWN')
	
			! Now read the model information
			isFileOpen = 1		! file is already open
	!		call OpenModelOutputFile_2(isFileOpen,1)		! This one reads affinities. 1 means read the node and seg XY coordinates
			call OpenModelOutputFile(isFileOpen,1)			! No affinities in this one. 1 means read the node and seg XY coordinates
	
			! scan the model information for the composition of garnet along each segment
			! Note that the segment endpoints align with nodes, so we don't need to tabulate the nodes
	
	
			Do iSeg = 1,numSegs
				iup = 0
				do k = 1,numSegReactionPhases(iSeg)
					if(segMIFID(iSeg,k).eq.iGarnet)then
	
						do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
							x = pointX(iSeg,iPointX)
							y = pointY(iSeg,iPointX)
							! check to see that these points are within bounds
							if(x.ge.xmin.and.x.le.xmax.and.y.ge.ymin.and.y.le.ymax)then		! this one is inside the search area
								! see if one of the phases is garnet
								! get the desired garnet composition
								GrtComp = pointPhaseComp(iSeg,iPointX,k,whatToPlot)
								write(12,*)iSeg,iPointX,x,y,k,GrtComp
								go to 45
								endif
							cycle	! This one is not inside the search area -- get the next point				
	
		45					continue
							if(GrtComp.gt.maxComp)then	! set the color to black
								i = 256
								else
								if(range.eq.0.0)then
									i = 1
									else
									i = 1 + int(increment*(GrtComp - minComp)/range)			! fractional amount above minimum
									endif
								endif
							pen%penColor = LUT_RGB(i)
							brush%brushColor = LUT_RGB(i)
							pen%penWidth = penWidth		! make the pen thicker
							call plot_2(GridPlot,x,y,iup,i)
							iup = 1
							end do		! end loop through points
						endif		! No garnet in this segment
					end do		! end loop on numReactionPhases
				end do  	! end loop on numSegs	
	
	 		go to 44		! skip the next stuff
	
	
		endif

	case default
		go to 1
	end select

	go to 1

	end




	


