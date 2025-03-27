!c*******************************************************
!c*******************************************************
	Subroutine OpenCanvasWindow(MyWindow)
    	use AWE_Interfaces
	implicit none
	TYPE(AWE_Canvas) :: MyWindow		! The plotting canvas


	MyWindow%title = 'Counter'
	MyWindow%width = 200
	MyWindow%height = 100
	MyWindow%backgroundColor = AWE_teal
	CALL AWE_createCanvas(MyWindow)	
	return
	end
!c*******************************************************
!c*******************************************************
	SUBROUTINE TextOnMyWindow(MyWindow,xpix,ypix,text,textSize)
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: MyWindow
	TYPE(AWE_CanvasPen) :: pen
	TYPE(AWE_Rect) :: rect
!	INTEGER, optional :: flags
	INTEGER :: flags
	TYPE(AWE_Font) :: font
	INTEGER :: textColor = AWE_black
	integer textSize
	real width
	integer xpix,ypix                ! locations in pixels
	CHARACTER(LEN=*) text
	
	font%pointSize = textSize
	rect%origin%x = xpix
	rect%origin%y = ypix
	text = adjustL(text)
	width = 1.4*float(textSize*(LEN(TRIM(text))))       !12 pixels per character

	rect%size%width = width
	rect%size%height = textSize
	flags = AWE_TextFlag_AlignLeft
	CALL AWE_canvasDrawText(MyWindow, rect, text, flags, font, textColor)
	return
	end
