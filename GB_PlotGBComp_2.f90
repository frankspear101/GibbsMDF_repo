!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine Plot_GB_Comp_2(GridPlot,autoPlot,iEl,minEl,maxEl)
	USE AWE_INTERFACES
      	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				
	INCLUDE "LUTColorStuff.inc"
! !	**********************************
! 	integer*4 LUTnumColors,LUT_RGB(300)
! 	real*4 LUT_CMYK(300,4)
! 	common /LUTcommon/LUTnumColors,LUT_RGB,LUT_CMYK
! !	**********************************

	real*4 x,y,dX,dY,increment,range,comp,minEl,maxEl,autoPlot,radius,penWidth
	integer*4 i,iEl,iSeg,iPointX,status,R,G,B,j,iup
	character*16 dummy

!	Open LUT file
! 	call FSS_Alert('Alert','Open LUT file')
! 	open(41,file='',iostat=status)
! 	if(status.ne.0)return
! !	Read the file
! 	read(41,*)LUTnumColors
! 	read(41,*)dummy
! 	do i = 1,LUTnumColors
! 		read(41,*)j,LUT_RGB(i),R,G,B,(LUT_CMYK(i,j),j=1,4)
! 		end do


	if(autoPlot.eq.1)go to 10
1	continue
	write(*,*)'Pick element to plot (0 to exit)(200 to save illustrator file)'
	write(*,*)' 0 = exit'
	do i = 1,numEl
		write(*,*)i,elName(i)
		end do
	write(*,*)' 100 = change pen width'
	write(*,*)' 200 = save Illustrator file'

	read(*,*)iEl
	if(iEl.eq.0)return
	if(iEl.eq.100)then
		write(*,*)'Input pen width (e.g., 0-10)'
		read(*,*)penWidth
		go to 1
		endif

	if(iEl.eq.200)then
		call SaveIllustratorFile()
		go to 1
		endif
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
	increment = LUTnumColors - 1
	dx = .002*(xmax-xmin)
	dy = dx
	radius = dx

	! Write out the lines on the canvas
	pen%penStyle = CanvasPenStyle_SolidLine
!	pen%penWidth = 2.0		! make the pen thicker
	pen%penWidth = penWidth		! make the pen thicker
	Do iSeg = 1,numSegs
		iup = 0
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
			x = pointX(iSeg,iPointX)
			y = pointY(iSeg,iPointX)
			comp = pointComp(iseg,iPointX,iEl)	
			if(range.eq.0.0)then
				i = 1
				else
				i = 1 + int(increment*(comp - minEl)/range)			! fractional amount above minimum
				endif
! 			write(*,*)iSeg,iPointX,comp,i
! 			pen%penColor = positiveColors(i)
! 			brush%brushColor = positiveColors(i)
			pen%penColor = LUT_RGB(i)
			pen%penWidth = penWidth		! make the pen thicker
			brush%brushColor = LUT_RGB(i)
			call plot_2(GridPlot,x,y,iup,i)
!			call pPlot_2(x,y,iup,i)
			iup = 1
			end do
!		call pPenup
		end do
	pen%penWidth = 0.0		! reset the pen width (I think 0 is the default, but I'm not sure

	! write out the lines in the Illustrator file
	! Each small segment between points needs to be a separate "line" so that each can be given a distinct color
	write(psunit,*)'u'	! this should "group" all of the line segments for later sizing
	Do iSeg = 1,numSegs
		iup = 0
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg) - 1
			x = pointX(iSeg,iPointX)
			y = pointY(iSeg,iPointX)
			comp = pointComp(iseg,iPointX,iEl)	
			if(range.eq.0.0)then
				i = 1
				else
				i = 1 + int(increment*(comp - minEl)/range)			! fractional amount above minimum
				endif
! 			write(*,*)iSeg,iPointX,comp,i
! 			pen%penColor = positiveColors(i)
! 			brush%brushColor = positiveColors(i)
!			pen%penColor = LUT_RGB(i)
!			brush%brushColor = LUT_RGB(i)
!			pen%penStyle = CanvasPenStyle_SolidLine
!			call plot_2(GridPlot,x,y,iup,i)

			call pPlot_2(x,y,0,i)		! move to the start the line
			x = pointX(iSeg,iPointX+1)
			y = pointY(iSeg,iPointX+1)
			call pPlot_2(x,y,1,i)		! draw the line to the next point
			call pPenup			! end the line
			iup = 1
			end do
		end do
	write(psunit,*)'U'		! Turn off "group" for this assemblage



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
! 			write(*,*)iSeg,iPointX,comp,i
! 			pen%penColor = positiveColors(i)
! 			brush%brushColor = positiveColors(i)
			pen%penColor = LUT_RGB(i)
			brush%brushColor = LUT_RGB(i)
			pen%penStyle = CanvasPenStyle_SolidLine

			Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
!			I need to set up the colors for this
			call psCircle_2(x,y,radius,3,i)	!SFB = 1 stroke only; =2 Fill only; =3 both
			end do
		end do



	go to 1
	
	end

!c*******************************************************
!c*******************************************************
	SUBROUTINE plot_2 (canvas,x,y,ipen,RGBcolor)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: canvas
	TYPE(AWE_Line) :: line(1)
	TYPE(AWE_CanvasPen) :: pen

	integer*4 ix,iy,ipen,ixold,iyold,RGBcolor
	real*4 x,y
	INCLUDE "PlotStuff.inc"	
	INCLUDE "LUTColorStuff.inc"
! !	**********************************
! 	integer*4 LUT_RGB(300)
! 	real*4 LUT_CMYK(300,4)
! 	common /LUTcommon/LUTnumColors,LUT_RGB,LUT_CMYK
! !	**********************************
	save
	ix=int(xor+(x-xmin)*xconv + .5)
	iy=int(yor-(y-ymin)*yconv + .5)

! 	pen%penColor = AWEcolorNumber(CurrentColor)		! Color for AWE line - in PlotStuff.inc common block
	pen%penColor = LUT_RGB(RGBcolor)		! Color for AWE line - in PlotStuff.inc common block

	IF (ipen.EQ.0) then
!		If this is a move, then the start and end of the line are the same
		line%start%x = ix
		line%start%y = iy
		line%end%x = ix
		line%end%y = iy
		CALL AWE_canvasDrawLines(canvas, line, pen)		
		else
!		If this is a line, then the start is the last point called
		line%start%x = ixold
		line%start%y = iyold
		line%end%x = ix
		line%end%y = iy
		CALL AWE_canvasDrawLines(canvas, line, pen)		
		ENDIF
	ixold=ix
	iyold=iy
	return
	END
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE pplot_2 (x,y,ipen,CMYKcolor)
	implicit none
!c	routine to move pen to x,y in user units
!c	ipen = 0 move with pen up
!c	ipen = 1 move with pen down
!c------------------------------------------------
	INCLUDE "PlotStuff.inc"	
	INCLUDE "LUTColorStuff.inc"
!c------------------------------------------------
! !	**********************************
! 	integer*4 LUT_RGB(300)
! 	real*4 LUT_CMYK(300,4)
! 	common /LUTcommon/LUTnumColors,LUT_RGB,LUT_CMYK
! !	**********************************
	integer*4 ipen,i,CMYKcolor
	real*4 x,y
	save

	IF (ipen.EQ.0)then
! 		write(psunit,5)(CMYK(CurrentColor,i),i=1,4)
		write(psunit,5)(LUT_CMYK(CMYKcolor,i),i=1,4)
5		format(4f10.5,'    K')			! Capital "K" sets the line color (lower case "k" sets fill color)
		write(psunit,10)pxor+(x-xmin)*pxconv,-pyor+(y-ymin)*pyconv
10		format(2f10.1,' m')	
		else
		write(psunit,20)pxor+(x-xmin)*pxconv,-pyor+(y-ymin)*pyconv
20		format(2f10.1,' L')
		endif
	return
	END

!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	
 	SUBROUTINE psCircle_2(xp,yp,radp,SFB,CMYKcolor)
!	Routine to draw a circle centered at centerx,centery with radius=radius
!		Units are user units
!		ColorS is color stroke CMYK
!		ColorF is color fill   CMYK
!		SFB 1=stroke = s; 2=fill = f; 3=both = b
!		radius is the radius in the X user direction
!			e.g. if cx = T = 500 and rad = 3 then the radius is 3 degrees
	implicit none
	INCLUDE "PlotStuff.inc"	
	INCLUDE "LUTColorStuff.inc"
! !	**********************************
! 	integer*4 LUT_RGB(300)
! 	real*4 LUT_CMYK(300,4)
! 	common /LUTcommon/LUTnumColors,LUT_RGB,LUT_CMYK
! !	**********************************
	integer*4 SFB,i,CMYKcolor
	real*4 xp,yp,radp,cx,cy,rad,x,y,t1x,t1y,t2x,t2y,tangk
	tangk = 0.55222966

!	write(psunit,*)'     0  0  0  1   K'		!black stroke
!	write(psunit,*)'     .611765 0 1 0 k'		! green fill
!	write(psunit,5)(colorS(i),i=1,4)
!	write(psunit,6)(colorF(i),i=1,4)

! 	write(psunit,5)(CMYK(CurrentColor,i),i=1,4)
! 	write(psunit,6)(CMYK(CurrentColor,i),i=1,4)
	write(psunit,5)(LUT_CMYK(CMYKcolor,i),i=1,4)
	write(psunit,6)(LUT_CMYK(CMYKcolor,i),i=1,4)
5	format(4f10.5,'    K')			! Capital "K" sets the line color (lower case "k" sets fill color)
6	format(4f10.5,'    k')

!												=centerx+radius	=centery	m
!=centerx+radius	=centery+radius*tangk	=centerx+radius*tangk	=centery+radius		=centerx	=centery+radius	c
!=centerx-radius*tangk	=centery+radius		=centerx-radius		=centery+radius*tangk	=centerx-radius	=centery	c
!=centerx-radius	=centery-radius*tangk	=centerx-radius*tangk	=centery-radius		=centerx	=centery-radius	c
!=centerx+radius*tangk	=centery-radius		=centerx+radius		=centery-radius*tangk	=centerx+radius	=centery	c

10	format(2f15.4,' m')	
20	format(6f15.4,' c')
! For things to be scaled correctly, cx, cy and radius need to first be converted to psuser units
!	write(*,*)pxconv,pyconv
!	write(*,*)xp,yp,radp
	cx =  pxor+(xp-xmin)*pxconv
	cy = -pyor+(yp-ymin)*pyconv
	rad = radp*pxconv
!	write(*,*)cx,cy,rad
! =centerx+radius	=centery	m
	x = cx + rad
	y = cy
!	x =  pxor+(x-xmin)*pxconv
!	y = -pyor+(y-ymin)*pyconv
	write(psunit,10)x,y
!=centerx+radius	=centery+radius*tangk	=centerx+radius*tangk	=centery+radius		=centerx	=centery+radius	c
	t1x = cx + rad
	t1y = cy + rad*tangk
	t2x = cx + rad*tangk
	t2y = cy + rad
	x = cx
	y = cy + rad
!	t1x =  pxor+(t1x-xmin)*pxconv
!	t1y = -pyor+(t1y-ymin)*pyconv
!	t2x =  pxor+(t2x-xmin)*pxconv
!	t2y = -pyor+(t2y-ymin)*pyconv
!	x =  pxor+(x-xmin)*pxconv
!	y = -pyor+(y-ymin)*pyconv
	write(psunit,20)t1x,t1y,t2x,t2y,x,y
!=centerx-radius*tangk	=centery+radius		=centerx-radius		=centery+radius*tangk	=centerx-radius	=centery	c
	t1x = cx - rad*tangk
	t1y = cy + rad
	t2x = cx - rad
	t2y = cy + rad*tangk
	x = cx - rad
	y = cy
!	t1x =  pxor+(t1x-xmin)*pxconv
!	t1y = -pyor+(t1y-ymin)*pyconv
!	t2x =  pxor+(t2x-xmin)*pxconv
!	t2y = -pyor+(t2y-ymin)*pyconv
!	x = pxor+(x-xmin)*pxconv
!	y = -pyor+(y-ymin)*pyconv
	write(psunit,20)t1x,t1y,t2x,t2y,x,y
!=centerx-radius	=centery-radius*tangk	=centerx-radius*tangk	=centery-radius		=centerx	=centery-radius	c
	t1x = cx - rad
	t1y = cy - rad*tangk
	t2x = cx - rad*tangk
	t2y = cy - rad
	x = cx
	y = cy - rad
!	t1x =  pxor+(t1x-xmin)*pxconv
!	t1y = -pyor+(t1y-ymin)*pyconv
!	t2x =  pxor+(t2x-xmin)*pxconv
!	t2y = -pyor+(t2y-ymin)*pyconv
!	x = pxor+(x-xmin)*pxconv
!	y = -pyor+(y-ymin)*pyconv
	write(psunit,20)t1x,t1y,t2x,t2y,x,y
!=centerx+radius*tangk	=centery-radius		=centerx+radius		=centery-radius*tangk	=centerx+radius	=centery	c
	t1x = cx + rad*tangk
	t1y = cy - rad
	t2x = cx + rad
	t2y = cy - rad*tangk
	x = cx + rad
	y = cy
!	t1x =  pxor+(t1x-xmin)*pxconv
!	t1y = -pyor+(t1y-ymin)*pyconv
!	t2x =  pxor+(t2x-xmin)*pxconv
!	t2y = -pyor+(t2y-ymin)*pyconv
!	x = pxor+(x-xmin)*pxconv
!	y = -pyor+(y-ymin)*pyconv
	write(psunit,20)t1x,t1y,t2x,t2y,x,y

	select case(SFB)
	case(1)		! stroke only
	write(psunit,*)'s'
	case(2)		! fill only
	write(psunit,*)'f'
	case(3)		! stroke and fill
	write(psunit,*)'b'
	case default
	end select

	return
	end	


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine PlotAffinityGrid(GridPlot)
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
	integer*4 i,j,iPointX,iSeg,iGarnet,iAff,iup,izero
	real*8 minAff,maxAff
	real*4 x,y,dX,dY,increment,Affinity,radius,penWidth,range
	real*4 pointAff(maxSegs,maxPoints)
	common /tempAff/pointAff
! !	**********************************
! 	integer*4 LUTnumColors,LUT_RGB(300)
! 	real*4 LUT_CMYK(300,4)
! 	common /LUTcommon/LUTnumColors,LUT_RGB,LUT_CMYK
! !	**********************************	
 

!	Open LUT file
! 	call FSS_Alert('Alert','Open LUT file')
! 	open(41,file='',iostat=status)
! 	if(status.ne.0)return
! !	Read the file
! 	read(41,*)LUTnumColors
! 	read(41,*)dummy
! 	do i = 1,LUTnumColors
! 		read(41,*)j,LUT_RGB(i),R,G,B,(LUT_CMYK(i,j),j=1,4)
! 		end do

	write(*,*)' Which phase is Garnet (0 to exit)?'
	do i = 1,numPhMIF
		write(*,*)i,phName(i)
		end do
	read(*,*)iGarnet		! this is the number of the garnet in the MIF

	if(iGarnet.eq.0)return

	! We work on points in every segment
	! The endpoints of the segments are the nodes, so this will include everything

	! Calculate the Affinities for every point and determine range of affinity values
	maxAff = -1.0e5
	minAff = 1.0d5		! make very large to start
	write(12,*)'**********************************************************'
	write(12,*)'**********************************************************'
	do iSeg = 1,numSegs
		write(12,*)'**********************************************************'
		write(12,*)'iSeg = ',iSeg
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
			do j = 1,numPhCo(1)		! grain boundary
				xPhCo(1,j) = pointComp(iSeg,iPointX,j)		! numPhCo(1) should be every element
				end do
! 			do j = 1,numPhCo(iGarnet)
! 				xPhCo(iGarnet,j) = pointPhaseComp(iSeg,iPointX,i,j)
! 				end do

			Call ParallelToTangent(iGarnet,0,izero)
			pointAff(iSeg,iPointX) = -gDifferenceAU(iGarnet)*0.666666		! This scales to J/mol-O for garnet	
			! NOTE: this will make affinities positive
			! This is being done to change the color spectrum so blue is low and red is high.
			! This seems more aesthetically pleasing and intuitive
			write(12,81)iPointX,-gDifferenceAU(iGarnet)*0.66666,(xPhCo(iGarnet,j),j=1,numPhCo(iGarnet))
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
	range = (maxAff-minAff)
	write(*,*)' Range = ',range


1	continue
	write(*,*)' Make your selection (0 to exit)(200 to save illustrator file)'
	write(*,*)' 0   = exit'
	write(*,*)' 1   = Plot affinities'
	write(*,*)' 2   = Adjust range'
	write(*,*)' 10  = ListLUT'
	write(*,*)' 100 = change pen width'
	write(*,*)' 200 = save Illustrator file'

	read(*,*)iAff
	if(iAff.eq.0)return
	if(iAff.eq.2)then
		write(*,*)'Input value for minAff'
		read(*,*)minAff
		write(*,*)'Input value for maxAff'
		read(*,*)maxAff
		write(*,*)' Minimum = ',minAff
		write(*,*)' Maximum = ',maxAff
		range = (maxAff-minAff)
		write(*,*)' Range = ',range
		go to 1
		endif		
	if(iAff.eq.10)then
		write(12,*)''
		write(12,*)' i   RGB        CMYK'
		do i = 1,256
			write(12,*)i,LUT_RGB(i),(LUT_CMYK(i,j),j=1,4)
			end do
		write(12,*)''
		write(12,*)''
		go to 1
		endif

	if(iAff.eq.100)then
		write(*,*)'Input pen width (e.g., 0-10)'
		read(*,*)penWidth
		go to 1
		endif

	if(iAff.eq.200)then
		call SaveIllustratorFile()
		go to 1
		endif


10	continue
! 	range = (maxAff-minAff)
! 	write(*,*)' Range = ',range
	increment = LUTnumColors - 1
	dx = .002*(xmax-xmin)
	dy = dx
	radius = dx

	! Write out the lines on the canvas
	pen%penStyle = CanvasPenStyle_SolidLine
	pen%penWidth = penWidth		! make the pen thicker
	Do iSeg = 1,numSegs
		iup = 0
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
			x = pointX(iSeg,iPointX)
			y = pointY(iSeg,iPointX)
			Affinity = pointAff(iseg,iPointX)	
			if(Affinity.gt.maxAff)then	! set the color to black
				i = 256
				else
				if(range.eq.0.0)then
					i = 1
					else
					i = 1 + int(increment*(Affinity - minAff)/range)			! fractional amount above minimum
					endif
				endif
			pen%penColor = LUT_RGB(i)
			brush%brushColor = LUT_RGB(i)
			pen%penWidth = penWidth		! make the pen thicker
			call plot_2(GridPlot,x,y,iup,i)
			iup = 1
			end do
		end do
	pen%penWidth = 0.0		! reset the pen width (I think 0 is the default, but I'm not sure

	! write out the lines in the Illustrator file
	! Each small segment between points needs to be a separate "line" so that each can be given a distinct color
	write(psunit,*)'u'	! this should "group" all of the line segments for later sizing
	Do iSeg = 1,numSegs
		iup = 0
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg) - 1
			x = pointX(iSeg,iPointX)
			y = pointY(iSeg,iPointX)
			Affinity = pointAff(iseg,iPointX)	
			if(Affinity.gt.maxAff)then	! set the color to black
				i = 256
				else
				if(range.eq.0.0)then
					i = 1
					else
					i = 1 + int(increment*(Affinity - minAff)/range)			! fractional amount above minimum
					endif
				endif

			call pPlot_2(x,y,0,i)		! move to the start the line
			x = pointX(iSeg,iPointX+1)
			y = pointY(iSeg,iPointX+1)
			call pPlot_2(x,y,1,i)		! draw the line to the next point
			call pPenup			! end the line
			iup = 1
			end do
		end do
	write(psunit,*)'U'		! Turn off "group" for this assemblage



	do iSeg = 1,numSegs
		do iPointX = numPointStart(iSeg),numPointEnd(iSeg)
			x = pointX(iSeg,iPointX)
			y = pointY(iSeg,iPointX)
			Affinity = pointAff(iseg,iPointX)	
			if(Affinity.gt.maxAff)then	! set the color to black
				i = 256
				else
				if(range.eq.0.0)then
					i = 1
					else
					i = 1 + int(increment*(Affinity - minAff)/range)			! fractional amount above minimum
					endif
				endif
			pen%penColor = LUT_RGB(i)
			brush%brushColor = LUT_RGB(i)
			pen%penStyle = CanvasPenStyle_SolidLine

			Call PlotCenteredEllipseOnScreen(GridPlot,X,dX,Y,dY,brush,pen)
!			I need to set up the colors for this
			call psCircle_2(x,y,radius,3,i)	!SFB = 1 stroke only; =2 Fill only; =3 both
			end do
		end do



	go to 1
	
	end


