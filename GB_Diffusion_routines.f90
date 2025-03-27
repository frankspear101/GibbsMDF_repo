	Subroutine GB_Diffuse_Hex()
	USE AWE_INTERFACES
	implicit none
	TYPE(AWE_Canvas) :: GridPlot
	TYPE(AWE_Canvas) :: SegPlot1
	TYPE(AWE_Canvas) :: SegPlot2
	TYPE(AWE_Canvas) :: SegPlot3
!	TYPE(AWE_Canvas) :: SegPlot(3)
	Type(AWE_CanvasBrush) :: brush
	TYPE(AWE_CanvasPen) :: pen
	include "Diffuse_GB_Hex.inc"
	INCLUDE "PlotStuff.inc"				! required for color file information
	integer*4 isegs(20),ichoose,i,numIter,error
!	real*8 rval
!	real*4 yMinMaxSi(2),yMinMaxAl(2),yMinMaxFM(2)	

!	call OpenModelFile()
	
!	ColorFileName = 'GB_Subs/GBColors.txt'
!	call ReadColorFile()
!	Call ReadDiffCoefs()
!	call GridMaker()
!	call SetInitialComp()
		
!	These are for setting up and plotting compositions along a series of segments
!	numSegsToPlot = 2
!	isegs(1) = 1
!	isegs(2) = 3
!	call SetupSegPlot(SegPlot,numSegsToPlot,isegs)
!	call PerturbNode()
!	Call PlotSegments(SegPlot,numSegsToPlot,isegs)

!	call CalcDatTP()
!	pause ' Look at diffusion coeffs'
!	dtime = 1.		! should ensure R < .1
!	dtime = 0.1		! should ensure R < .1
!	dtime = 0.01		! should ensure R < .1
	dtime = 40.		! should ensure R < .5

! 	dtime = 1.d5
! 	dtime = 0.1d0*dXgrid**2/DatTP(1)		! set rval = 0.1 based on the largest D (H) diffusivity
! 	rval = dtime*DatTP(1)/dXGrid**2
! 	write(*,*)' dtime, rval'
! 	write(*,*)dtime,rval

1	continue
	write(*,*)'Main menu'
	write(*,*)' 1 = list node'
	write(*,*)' 2 = list segment'
	write(*,*)' 3 = Setup plot'
	write(*,*)' 4 = Plot segment'
	write(*,*)' 5 = diffuse node'
	write(*,*)' 6 = diffuse segment'
	write(*,*)' 7 = run n diffusion iterations'
	write(*,*)' 8 = Plot and label grid'
	write(*,*)' 9 = Perturb system'

	read(*,*)ichoose

	select case(ichoose)
	case(0)
		return
	case(1)
		call ListNode()
	case(2)	
		call ListSegment()
	case(3)
		numSegsToPlot = 8  
		iSegs(1) = 44
		iSegs(2) = 62
		iSegs(3) = 79
		iSegs(4) = 97
		iSegs(5) = 115
		iSegs(6) = 116
		iSegs(7) = 133
		iSegs(8) = 151
		go to 31		
		i = 1
		write(*,*)' input segment numbers -- end with a 0'
30		continue
		read(*,*)isegs(i)
		if(isegs(i).eq.0)then
			numSegsToPlot = i - 1
			go to 31
			endif
		go to 30
31		continue
		call SetupSegPlot(SegPlot1,SegPlot2,SegPlot3)
!		call SetupSegPlot(SegPlot,numSegsToPlot,isegs)
	case(4)
		Call PlotSegments(SegPlot1,SegPlot2,SegPlot3)
!		Call PlotSegments(SegPlot,numSegsToPlot,isegs)
	case(5)
		write(*,*)' Input number of iterations to run (0 to abort)'
		read(*,*)numIter
		if(numIter.eq.0)go to 1
		do 50 i = 1,numIter
		call DiffuseNodes(error)
50		continue
	case(6)	
		write(*,*)' Input number of iterations to run (0 to abort)'
		read(*,*)numIter
		if(numIter.eq.0)go to 1
		do 60 i = 1,numIter
		call DiffuseSegments()
60		continue
	case(7)	
		write(*,*)' Input number of iterations to run (0 to abort)'
		read(*,*)numIter
		if(numIter.eq.0)go to 1
		do i = 1,numIter
			call DiffuseNodes(error)
			call DiffuseSegments()
			end do
	case(8)	
80		continue
		write(*,*)' 0 = return'
		write(*,*)' 1 = SetUpGridPlot'
		write(*,*)' 2 = plot grid (nodes)'
		write(*,*)' 3 = plot points'
		write(*,*)' 4 = label nodes'
		write(*,*)' 5 = label segments'
		write(*,*)' 6 = Draw Crystals'
		write(*,*)' 7 = Label Crystals'
		read(*,*)numIter
		if(numIter.eq.0)go to 1
		if(numIter.eq.1)call SetUpGridPlot(GridPlot)
		if(numIter.eq.2)call PlotGrid(GridPlot,1)
		if(numIter.eq.3)call PlotPoints(GridPlot,1,0,0)		! one means pick color; 0 = ask to connect the points
		if(numIter.eq.4)call LabelNodes(GridPlot)
		if(numIter.eq.5)call LabelSegs(GridPlot)
		if(numIter.eq.6)call DrawCrystals(GridPlot)
		if(numIter.eq.7)call LabelCrystals(GridPlot)
		go to 80
	case(9)
		call PerturbNode()
	case default
		call FSS_Alert('Alert!!!','You chose poorly.....')
	end select	
	go to 1
	end


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine PerturbNode()
	implicit none
	integer*4 i,iseg,ipoint,j
	include "Diffuse_GB_Hex.inc"

	write(12,*)' Perturbing a node'
	write(12,*)' Node to perturb ',nodeToPerturb
	do j = 1,numEl
		nodeComp(nodeToPerturb,j) = PerturbComp(j)
		write(12,*)j,nodeComp(nodeToPerturb,j)
		end do

	do i = 1,numNodeSegs(nodeToPerturb)
		iseg = nodeSegConnect(nodeToPerturb,i)
		ipoint = nodePointOnTop(nodeToPerturb,i)
		! iseg and ipoint are now the points that fall on the same place as the node
		write(12,*)i,iseg,ipoint
		do j = 1,numEl
			pointComp(iseg,ipoint,j) = nodeComp(nodeToPerturb,j)
			end do
		end do

	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine ListNode()
	implicit none
	integer*4 iNode,j
	include "Diffuse_GB_Hex.inc"
1	continue
	write(*,*)'Input node number to list. 0 to return'
	read(*,*)iNode
	if(iNode.eq.0)return
	write(*,*)' '
	write(*,*)' Node position and composition'
	write(*,11)(ElName(j),j=1,numEl)
11	format(15x,20A12)
	write(*,12)iNode,nodeX(iNode),nodeY(iNode),(nodeComp(iNode,j),j=1,numEl)
12	format(I5,2F8.3,20F12.4)
10	continue
	go to 1
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine ListSegment()
	implicit none
	integer*4 i,iseg,j
	include "Diffuse_GB_Hex.inc"
1	continue
	write(*,*)'Input segment number to list. 0 to return'
	read(*,*)iseg
	if(iseg.eq.0)return
	write(*,*)'Connecting nodes and positions'
	write(*,*)segNodes(iseg,1),nodeX(segNodes(iseg,1)),nodeY(segNodes(iseg,1))
	write(*,*)segNodes(iseg,2),nodeX(segNodes(iseg,2)),nodeY(segNodes(iseg,1))
	write(*,*)' '
	write(*,*)' Segment point positions'
	write(*,11)(ElName(j),j=1,numEl)
11	format(15x,20A12)
	do i = numPointStart(iSeg),numPointEnd(iSeg)
		!do i = 1,numPoints(iseg)
		write(*,12)i,pointX(iseg,i),pointY(iseg,i),(pointComp(iseg,i,j),j=1,numEl)
12		format(I5,2F8.3,20F12.4)
		end do
	go to 1
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$	
	Subroutine SetInitialComp()
	implicit none
	integer*4 i,j,k,iSeg
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
!	Initialize the composition of the grainboundary.
!	Note that the initial GB composition is read from the MIF

	if(numEl.ne.numPhCo(1))then
		call FSS_Alert('ALERT',' GB phase components does not match number of elements...abort')
		return
		endif
		
	do j = 1,numPhCo(1)		! this is the grain boundary phase
		InitialComp(j) = xPhCo(1,j)
		end do

	do i = 1,numNodes
		do j = 1,numEl
			nodeComp(i,j) = InitialComp(j)
			if(i.eq.93)then
				write(*,*)j,nodeComp(i,j)
				endif
			end do
		end do
	do iSeg = 1,numSegs
		do k = numPointStart(iSeg),numPointEnd(iSeg)	! number of points in segment(i)
		!do k = 1,numPoints(iSeg)	! number of points in segment(i)
			do j = 1,numEl
				pointComp(iSeg,k,j) = InitialComp(j)
				end do
			end do
		end do

	end



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine CalcDatTP(TC,PB)
      implicit none
!*****************************************
	include "Diffuse_GB_Hex.inc"
!*****************************************
	integer*4 i
	real*8 Rjoules,TC,PB,TK,rval
	parameter (Rjoules = 8.3144D0)

!	TC = 500.
!	PB = 1.
	TK = TC + 273.15
	write(*,*)' TC, TK, Pb'
	write(*,*)TC,TK,PB
	write(*,*)' El   DatTP'
	! Note: DactE is read in as kJ and converted to J after it is read in
	do 10 i = 1,numEl
	DatTP(i) = Dzero(i)*exp((-DactE(i)+DactVol(i)*Pb)/(Rjoules*TK))
	write(*,*)ElName(i),DatTP(i)
10	continue

	dtime = 1.d5
	dtime = 0.1d0*dXgrid**2/DatTP(1)		! set rval = 0.1 based on the largest D (H) diffusivity
	rval = dtime*DatTP(1)/dXGrid**2
	write(*,*)' dtime, rval'
	write(*,*)dtime,rval

	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine IdealDIJPoint(iseg,ipoint)
      implicit none
!     Program to calculate the values of multicomponent
!     diffusion coefficient matrix elements at designated P-T
!     conditions employing the multicomponent diffusion model
!     proposed by Lasaga (1979), assuming ideal solution behavior.

!	This routine is called for every point in the grain boundary grid when needed
!*****************************************
	include "Diffuse_GB_Hex.inc"
!*****************************************
!     Local variables
      real*8 D_Ddep(maxElements),DENOM,TEMP
      integer*4 j,k,iseg,ipoint

!     calculate difference between independent and dependent D's
      	do 10 j=2,numEl
      	D_Ddep(j) = DatTP(j) - DatTP(1)
10    	continue
!  Begin calculation of diffusion coefficient matrix elements,
!  (Lasaga (1979), eqn. 14b).  Mg dependent case.
	DENOM = 0.0D0
	do 110 j=1,numEl
	DENOM = DENOM + (ElCharge(j)**2)*pointComp(iseg,ipoint,j)*DatTP(j)
110   	continue
!  Matrix elements:
      	do 120 j=2,numEl
	temp  = DatTP(j)*pointComp(iseg,ipoint,j)/DENOM
      	do 125 k = 2,numEl
      	DMatrix(j,k) =  -ElCharge(j)*ElCharge(k)*temp*D_Ddep(k)     
125   	continue
      	DMatrix(j,j) =  DatTP(j) - DMatrix(j,j)
      	DMatrix(j,j) =  DatTP(j)			! this uses just the tracer diffusivity
120   	continue
	return
	do 150 j = 2,numEl
	write(*,155)(Dmatrix(j,k),k=2,numEl)
150	continue
155	format(14E12.4)
	pause 'look at Dmatrix'
      	RETURN
      	END

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine IdealDIJNode(inode)
      implicit none
!     Program to calculate the values of multicomponent
!     diffusion coefficient matrix elements at designated P-T
!     conditions employing the multicomponent diffusion model
!     proposed by Lasaga (1979), assuming ideal solution behavior.

!	This routine is called for every node in the grain boundary grid when needed
!*****************************************
	include "Diffuse_GB_Hex.inc"
!*****************************************
!     Local variables
      real*8 D_Ddep(maxElements),DENOM,TEMP
      integer*4 j,k,inode

!     calculate difference between independent and dependent D's
      	do 10 j=2,numEl
      	D_Ddep(j) = DatTP(j) - DatTP(1)
10    	continue
!  Begin calculation of diffusion coefficient matrix elements,
!  (Lasaga (1979), eqn. 14b).  Mg dependent case.
	DENOM=0.0D0
	do 110 j=1,numEl
	DENOM = DENOM +  (ElCharge(j)**2)*nodeComp(inode,j)*DatTP(j)
110   	continue
!  Matrix elements:
      	do 120 j=2,numEl
	temp  = DatTP(j)* nodeComp(inode,j)/DENOM
      	do 125 k = 2,numEl
      	DMatrix(j,k) =  -ElCharge(j)*ElCharge(k)*temp*D_Ddep(k)     
125   	continue
      	DMatrix(j,j) =  DatTP(j) - DMatrix(j,j)
      	DMatrix(j,j) =  DatTP(j)			! this uses just the tracer diffusivity
120   	continue

	return
	do 150 j = 2,numEl
	write(*,155)(Dmatrix(j,k),k=2,numEl)
150	continue
155	format(14E12.4)
	pause 'look at Dmatrix'

      	RETURN
      	END

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine DiffuseNodes(error)
      implicit none
!	Routine to calculate the change in composition of all nodes

!*****************************************
	include "Diffuse_GB_Hex.inc"
!*****************************************
!     Local variables
	real*8 CompChange(maxElements),sum,dXgrid2,			&
		curvature12(maxElements),curvature13(maxElements),curvature23(maxElements)
	integer*4 j,jj,iNode,error
	integer*4 iseg1,ipoint1,iseg2,ipoint2,iseg3,ipoint3

	dXgrid2 = dXgrid*dXgrid
!	Note that the grid spacing is fixed at a value of dXgrid

	do 10 iNode = 1,numNodes

!	if(iNode.eq.2)then
!		write(*,*)' '
!		write(*,*)' Node = ',iNode
!		write(*,*)'Before ' ,(NodeComp(iNode,j),j=1,numEl)
!		endif

!	call IdealDIJNode(iNode)		! calculate diffusion coeffs for this node

	select case (numNodeSegs(iNode))

	case(2)			! node connects to only 2 others
!				Code is the same as for a segment diffusion
	iseg1 = nodeSegConnect(iNode,1)
	ipoint1 = nodePointNextTo(iNode,1)
	iseg2 = nodeSegConnect(iNode,2)
	ipoint2 = nodePointNextTo(iNode,2)
!	curvature	
	do 20 j = 2,numEl
!	The first element (typically Si) is the dependent one
!      curvature(i,j) = (U(i+1,j) - 2.0D0*U(i,j) + U(i-1,j))/DELX2
	curvature12(j) = (pointComp(iseg1,ipoint1,j) - 2.0D0*nodeComp(iNode,j) + pointComp(iseg2,ipoint2,j))/dXgrid2
20	continue

	do 25 j = 2,numEl
!	sum = 0.0d0
!	do 26 jj = 2,numEl
!	jj = j			! this will only calculate with diagonal terms
!	sum = sum + curvature12(jj) * Dmatrix(j,jj)
!26	continue
!	CompChange(j) = sum * dTime					! uses full D matrix
!	CompChange(j) = curvature12(j) * Dmatrix(j,j) * dTime		! uses diagonal elements only
	CompChange(j) = curvature12(j) * DatTP(j) * dTime		! uses diagonal elements only
25	continue
	sum = 0.0d0
	do 28 j = 2,numEl
	nodeComp(iNode,j) = nodeComp(iNode,j) + CompChange(j)
	sum = sum + nodeComp(iNode,j)
28	continue
	nodeComp(iNode,1) = 1.0d0 - sum	! the dependent composition variable

	go to 29
	if(iNode.eq.1.or.iNode.eq.2.or.iNode.eq.3)then
		write(*,*)' '
		write(*,*)' Node = ',iNode
		write(*,*)' Seg    point'
		write(*,*)iseg1,ipoint1,(pointComp(iseg1,ipoint1,j),j = 2,numEl)
		write(*,*)iseg2,ipoint2,(pointComp(iseg2,ipoint2,j),j = 2,numEl)
		write(*,*)'C12 ',(curvature12(j),j = 2,numEl)
		write(*,*)'D   ',(Dmatrix(2,j), j = 2,numEl)
		write(*,*)'Ch  ',(CompChange(j),j=2,numEl)
		write(*,*)sum,nodeComp(iNode,1)
		endif
29	continue

	case(3)		! node connects to three others

	iseg1 = nodeSegConnect(iNode,1)
	ipoint1 = nodePointNextTo(iNode,1)
	iseg2 = nodeSegConnect(iNode,2)
	ipoint2 = nodePointNextTo(iNode,2)
	iseg3 = nodeSegConnect(iNode,3)
	ipoint3 = nodePointNextTo(iNode,3)
!	curvature	
	do j = 2,numEl
	!	The first element (typically Si) is the dependent one
	!      curvature(i,j) = (U(i+1,j) - 2.0D0*U(i,j) + U(i-1,j))/DELX2
		curvature12(j) = (pointComp(iseg1,ipoint1,j) - 2.0D0*nodeComp(iNode,j) + pointComp(iseg2,ipoint2,j))/dXgrid2
		curvature13(j) = (pointComp(iseg1,ipoint1,j) - 2.0D0*nodeComp(iNode,j) + pointComp(iseg3,ipoint3,j))/dXgrid2
		curvature23(j) = (pointComp(iseg2,ipoint2,j) - 2.0D0*nodeComp(iNode,j) + pointComp(iseg3,ipoint3,j))/dXgrid2
		end do


	do j = 2,numEl
		CompChange(j) = (curvature12(j) + curvature13(j) + curvature23(j)) * DatTP(j) * dTime		! diagonal elements only
		end do

	sum = 0.0d0
	do j = 2,numEl
		nodeComp(iNode,j) = nodeComp(iNode,j) + CompChange(j)
		sum = sum + nodeComp(iNode,j)
		end do
	nodeComp(iNode,1) = 1.0d0 - sum	! the dependent composition variable

	 go to 39
!	if(iNode.eq.1.or.iNode.eq.2.or.iNode.eq.3)then
	if(iNode.eq.2)then
		write(*,*)' '
		write(*,*)' Node = ',iNode
		write(*,*)' Seg    point'
		write(*,*)iseg1,ipoint1,(pointComp(iseg1,ipoint1,j),j = 2,numEl)
		write(*,*)iseg2,ipoint2,(pointComp(iseg2,ipoint2,j),j = 2,numEl)
		write(*,*)iseg3,ipoint3,(pointComp(iseg3,ipoint3,j),j = 2,numEl)
		write(*,*)'C12  ',(curvature12(j),j = 2,numEl)
		write(*,*)'C13  ',(curvature13(j),j = 2,numEl)
		write(*,*)'C23  ',(curvature23(j),j = 2,numEl)
		write(*,*)'Dmat '
		do 33 jj = 2,numEl
		write(*,*)(Dmatrix(jj,j), j = 2,numEl)
33		continue
		write(*,*)'Change ',(CompChange(j),j=2,numEl)
		write(*,*)'Comp   ',(NodeComp(iNode,j),j=2,numEl)
		write(*,*)sum,nodeComp(iNode,1)
		write(*,*)'Node Comp (aft)   ',(NodeComp(iNode,j),j=1,numEl)
		endif
39	continue
	case default
		call FSS_Alert('Alert',' Error in Node Diffusion routine')
		write(*,*)'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
		write(*,*)' The number of NodeSegs is not 2 or 3 '
		write(*,*)'iNode , numNodeSegs(iNode),nodeSegConnect(iNode,1)-3'
		write(*,*)iNode , numNodeSegs(iNode),nodeSegConnect(iNode,1),nodeSegConnect(iNode,2),nodeSegConnect(iNode,3)
		error = 1
		return
	end select
	
10	continue


	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine DiffuseNodesVardX()
      implicit none
!	Routine to calculate the change in composition of all nodes
!	This version uses a variable grid spacing, which is necessary because
!	every segment has a different grid spacing

!*****************************************
	include "Diffuse_GB_Hex.inc"
!*****************************************
!     Local variables
	real*8 CompChange(maxElements),sum,dXgrid2,			&
		curvature12(maxElements),curvature13(maxElements),curvature23(maxElements),	&
		dXgrid1Squared,dXgrid2Squared,dXgridTotal
	integer*4 j,jj,iNode
	integer*4 iseg1,ipoint1,iseg2,ipoint2,iseg3,ipoint3

	dXgrid2 = dXgrid*dXgrid
!	Note that the grid spacing is fixed at a value of dXgrid

	do 10 iNode = 1,numNodes

!	if(iNode.eq.2)then
!		write(*,*)' '
!		write(*,*)' Node = ',iNode
!		write(*,*)'Before ' ,(NodeComp(iNode,j),j=1,numEl)
!		endif

!	call IdealDIJNode(iNode)		! calculate diffusion coeffs for this node

	select case (numNodeSegs(iNode))

	case(2)			! node connects to only 2 others
!				Code is the same as for a segment diffusion
	dXgrid2 = dXgrid*dXgrid

	iseg1 = nodeSegConnect(iNode,1)
	ipoint1 = nodePointNextTo(iNode,1)
	dXgrid1Squared = segdXgrid(iSeg1)*segdXgrid(iSeg1)

	iseg2 = nodeSegConnect(iNode,2)
	ipoint2 = nodePointNextTo(iNode,2)
	dXgrid2Squared = segdXgrid(iSeg2)*segdXgrid(iSeg2)

	dXgridTotal = segdXgrid(iSeg1) + segdXgrid(iSeg2)
	
!	curvature	
	do 20 j = 2,numEl
!	The first element (typically H) is the dependent one
!      curvature(i,j) = (U(i+1,j) - 2.0D0*U(i,j) + U(i-1,j))/DELX2
!	This code is for fixed grid spacing -- I need to put in the variable grid spacing code (at home in the notes)
	curvature12(j) = (pointComp(iseg1,ipoint1,j) - 2.0D0*nodeComp(iNode,j) + pointComp(iseg2,ipoint2,j))/dXgrid2
	

20	continue

	do 25 j = 2,numEl
!	sum = 0.0d0
!	do 26 jj = 2,numEl
!	jj = j			! this will only calculate with diagonal terms
!	sum = sum + curvature12(jj) * Dmatrix(j,jj)
!26	continue
!	CompChange(j) = sum * dTime					! uses full D matrix
!	CompChange(j) = curvature12(j) * Dmatrix(j,j) * dTime		! uses diagonal elements only
	CompChange(j) = curvature12(j) * DatTP(j) * dTime		! uses diagonal elements only
25	continue
	sum = 0.0d0
	do 28 j = 2,numEl
	nodeComp(iNode,j) = nodeComp(iNode,j) + CompChange(j)
	sum = sum + nodeComp(iNode,j)
28	continue
	nodeComp(iNode,1) = 1.0d0 - sum	! the dependent composition variable

	go to 29
	if(iNode.eq.1.or.iNode.eq.2.or.iNode.eq.3)then
		write(*,*)' '
		write(*,*)' Node = ',iNode
		write(*,*)' Seg    point'
		write(*,*)iseg1,ipoint1,(pointComp(iseg1,ipoint1,j),j = 2,numEl)
		write(*,*)iseg2,ipoint2,(pointComp(iseg2,ipoint2,j),j = 2,numEl)
		write(*,*)'C12 ',(curvature12(j),j = 2,numEl)
		write(*,*)'D   ',(Dmatrix(2,j), j = 2,numEl)
		write(*,*)'Ch  ',(CompChange(j),j=2,numEl)
		write(*,*)sum,nodeComp(iNode,1)
		endif
29	continue

	case(3)		! node connects to three others

	iseg1 = nodeSegConnect(iNode,1)
	ipoint1 = nodePointNextTo(iNode,1)
	iseg2 = nodeSegConnect(iNode,2)
	ipoint2 = nodePointNextTo(iNode,2)
	iseg3 = nodeSegConnect(iNode,3)
	ipoint3 = nodePointNextTo(iNode,3)
!	curvature	
	do 30 j = 2,numEl
!	The first element (typically Si) is the dependent one
!      curvature(i,j) = (U(i+1,j) - 2.0D0*U(i,j) + U(i-1,j))/DELX2
	curvature12(j) = (pointComp(iseg1,ipoint1,j) - 2.0D0*nodeComp(iNode,j) + pointComp(iseg2,ipoint2,j))/dXgrid2
	curvature13(j) = (pointComp(iseg1,ipoint1,j) - 2.0D0*nodeComp(iNode,j) + pointComp(iseg3,ipoint3,j))/dXgrid2
	curvature23(j) = (pointComp(iseg2,ipoint2,j) - 2.0D0*nodeComp(iNode,j) + pointComp(iseg3,ipoint3,j))/dXgrid2
30	continue

	do 35 j = 2,numEl
!	sum = 0.0d0
!	do 36 jj = 2,numEl
!	jj = j			! this will only calculate with diagonal terms
!	sum = sum + (curvature12(jj) + curvature13(jj) + curvature23(jj)) * Dmatrix(j,jj)
!36	continue
!	CompChange(j) = sum * dTime
!	CompChange(j) = (curvature12(j) + curvature13(j) + curvature23(j)) * Dmatrix(j,j) * dTime		! diagonal elements only
	CompChange(j) = (curvature12(j) + curvature13(j) + curvature23(j)) * DatTP(j) * dTime		! diagonal elements only
35	continue
	sum = 0.0d0
	do 38 j = 2,numEl
	nodeComp(iNode,j) = nodeComp(iNode,j) + CompChange(j)
	sum = sum + nodeComp(iNode,j)
38	continue
	nodeComp(iNode,1) = 1.0d0 - sum	! the dependent composition variable

	 go to 39
!	if(iNode.eq.1.or.iNode.eq.2.or.iNode.eq.3)then
	if(iNode.eq.2)then
		write(*,*)' '
		write(*,*)' Node = ',iNode
		write(*,*)' Seg    point'
		write(*,*)iseg1,ipoint1,(pointComp(iseg1,ipoint1,j),j = 2,numEl)
		write(*,*)iseg2,ipoint2,(pointComp(iseg2,ipoint2,j),j = 2,numEl)
		write(*,*)iseg3,ipoint3,(pointComp(iseg3,ipoint3,j),j = 2,numEl)
		write(*,*)'C12  ',(curvature12(j),j = 2,numEl)
		write(*,*)'C13  ',(curvature13(j),j = 2,numEl)
		write(*,*)'C23  ',(curvature23(j),j = 2,numEl)
		write(*,*)'Dmat '
		do 33 jj = 2,numEl
		write(*,*)(Dmatrix(jj,j), j = 2,numEl)
33		continue
		write(*,*)'Change ',(CompChange(j),j=2,numEl)
		write(*,*)'Comp   ',(NodeComp(iNode,j),j=2,numEl)
		write(*,*)sum,nodeComp(iNode,1)
		write(*,*)'Node Comp (aft)   ',(NodeComp(iNode,j),j=1,numEl)
		endif
39	continue
	case default
		call FSS_Alert('Alert',' Error in Node Diffusion routine')
		return
	end select
	
10	continue


	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine DiffuseSegments()
      implicit none
!	Routine to calculate the change in composition of points along a segment

!*****************************************
	include "Diffuse_GB_Hex.inc"
!*****************************************
!     Local variables
	real*8 CompChange(maxPoints,maxElements),curvature(maxElements),dXgrid2,sum
	integer*4 i,j,jj
	integer*4 iseg

!	dXgrid2 = dXgrid*dXgrid
!	Note that the grid spacing is fixed at a value of dXgrid
!      DELX2=DELX**2
	
	do iseg = 1,numSegs
		dXgrid2 = segdXgrid(iSeg)*segdXgrid(iSeg)

		!do i = 2,numPoints(iseg)-1
		do i = numPointStart(iSeg)+1,numPointEnd(iSeg)-1
			!Note that the first and last points on a segment fall on the nodes
			!We don't actually calculate changes in the composition of the first and last points here
			!The change in composition of nodes is in another subroutine and we set the new compositions
			!	of these first and last points there

			!call IdealDIJPoint(iseg,i)
			!curvature	
			do j = 2,numEl
				! The first element (typically Si) is the dependent one
				!      curvature(i,j) = (U(i+1,j) - 2.0D0*U(i,j) + U(i-1,j))/DELX2
				curvature(j) = (pointComp(iseg,i+1,j) - 2.0D0*pointComp(iseg,i,j) + pointComp(iseg,i-1,j))/dXgrid2
				end do

			do j = 2,numEl
				CompChange(i,j) = curvature(j) * DatTP(j) * dTime			! diagonal elements only
				!	store the composition change for every point in this segment
				end do

			go to 10
			if(iseg.eq.3.and.i.eq.2)then
				do jj = 2,numEl
					write(*,*)'D ',(Dmatrix(jj,j), j = 2,numEl)
					end do
				write(*,*)'C ',(curvature(j),j = 2,numEl)
				write(*,*)'Change',(compChange(i,j),j=2,numEl)
				endif
10			continue
			end do

		!do i = 2,numPoints(iseg) - 1
		do i = numPointStart(iSeg)+1,numPointEnd(iSeg)-1
			sum = 0.0d0
			do j = 2,numEl
				pointComp(iseg,i,j) = pointComp(iseg,i,j) + CompChange(i,j)
				sum = sum + pointComp(iseg,i,j)
				end do
			pointComp(iseg,i,1) = 1.0d0 - sum	! the dependent composition variable
			end do
		end do

!5	continue


!	Actually, this is done after the MDFPairs routine so it doesn't need to be done here
!	---After the MDF routines are done, the GB composition of the node is set to the point on top
! Set the first and last point compositions to that of the new node composition
! 	do iNode = 1,numNodes
! 		do i = 1,numNodeSegs(iNode)
! 			iseg   = nodeSegConnect(iNode,i)
! 			ipoint = nodePointOnTop(iNode,i)
! 				do j = 1,numEl
! 					pointComp(iseg,ipoint,j) = nodeComp(iNode,j)
! 					end do
! 			end do
! 		end do
	


	return
	end
	

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine AverageGBComp()
      implicit none
!	Routine to calculate the change in composition of points along a segment

!*****************************************
	include "Diffuse_GB_Hex.inc"
!*****************************************
!     Local variables
	real*8 GBAvg(maxElements)
	integer*4 i,j,iSeg,iNode
	integer*4 numInAvg

	! set initial values to zero
	do j = 1,numEl
		GBAvg(j) = 0.0d0
		end do

	numInAvg = 0

	do iseg = 1,numSegs
		do i = numPointStart(iSeg),numPointEnd(iSeg)
			do j = 1,numEl
				GBAvg(j) = GBAvg(j) + pointComp(iseg,i,j)
				end do
			numInAvg = numInAvg + 1
			end do
		end do

	! calculate averages
	do j = 1,numEl
		GBAvg(j) = GBAvg(j)/float(numInAvg)
		end do

	! Set segment points to average
	do iseg = 1,numSegs
		do i = numPointStart(iSeg),numPointEnd(iSeg)
			do j = 1,numEl
				pointComp(iseg,i,j) = GBAvg(j)
				end do
			end do
		end do

!	set node compositions to average
	do iNode = 1,numNodes
		do j = 1,numEl
			nodeComp(iNode,j) = GBAvg(j) 
			end do
		end do


	return
	end