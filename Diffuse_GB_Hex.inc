!	Include file for Diffuse_GB-Hex
	integer*4 maxNodes,maxSegs,minPoints,maxPoints,maxElements,maxCrystals,maxPhases,maxPhCo,maxPolys,maxPolyNodes
	parameter (maxNodes = 1000, maxSegs = 1000,minPoints=-30,maxPoints = 200,maxElements = 10)
	parameter (maxCrystals = 500, maxPhases = 10,maxPhCo = 10)
	parameter (maxPolys = 500)
	parameter (maxPolyNodes = 40)
	
!	Variables explanation
	
!	node variables
	Integer*4					&!
		numNodeSegs(maxNodes),			&! = number of segments that a node connects to ( 2 or 3)
		nodeNodeConnect(maxNodes,3), 		&! = index of the nodes that connect to 2 or 3 other nodes
		nodeSegConnect(maxNodes,3),		&!   = index of the segments connected to a node. Each node can have up to 3 segments attached -- 
									! counterclockwise order
			! nodeSegEnd(maxNodes,3),			&!  1 = start at beginning of segment, 2 = start at end of segment -- no longer used
		nodePointOnTop(maxNodes,3),		&! = index of the coincident point to the node in the point arrays (either first or last) 
		nodePointNextTo(maxNodes,3),		&! = index of the adjacent point in the point arrays  
		numCrystalsAtNode(maxNodes),		& ! number of crystals (or Polys) around a node. Should be 1, 2 or 3
		numNodes		 		 ! = self evident

	real*8						&
		nodeX(maxNodes),nodeY(maxNodes),	& ! = x,y coordinates of the nodes of the hex grid
		nodeXOrig(maxNodes),nodeYOrig(maxNodes) ! = Original x,y coordinates of the nodes of the hex grid (before scaling) - never change


!	Segment variables
	integer*4					&!
		numSegs,				&! = number of total segments
		segNodes(maxSegs,2),			&! = index of the 2 nodes that connect each segment
		numPointStart(maxSegs),			&! = point number where segment starts
		numPointEnd(maxSegs)			 ! = point number where segment ends (numPoints = end-start + 1)

	Real*8						&!
		segX(maxSegs,2),			& !The X,Y values of the endpoints of the segment -- same as the node X,Y
		segY(maxSegs,2),			&!
		segLength(maxSegs), 			&! = length of each segment between two nodes. Calculated in Sub SegLengthCalculate() at start and after each time nodes are moved (Sub MoveNode in GrowNodes)
		pointX(maxSegs,minPoints:maxPoints),	&!
		pointY(maxSegs,minPoints:maxPoints),	& ! = x,y coordinates of each point along a given segment
		segdXgrid(maxSegs),			& ! The grid spacing for each segment
		dXGrid,					& ! = grid spacing between points (constant throughout the entire grid)  
		dGridRemove,				& ! The fraction of dXGrid that will trigger the removal of a grid point. Only used in Sub RemoveSegPoints
		dGridAdd				 ! The fraction of dXGrid that will trigger the addition of a grid point. Only used in Sub RemoveSegPoints

	logical*4 					&
		segCaptureExclusion(maxSegs),		& ! nodecapture should never happen along these segments. These are read in from the model input file
		segGettingShorter(maxSegs)		! Logical variable to flag whether a seg (the distance between 2 nodes) is getting shorter
							! This is used to flag when node capture will occur (when 2 nodes get too close together)
							! See subroutine NodeCapture
							! This flag is also set when a node is captured so it only will happen one time (sub NodeCapture)

!	Polygon variables
	integer*4					&!
		numPolys,				&
		numPolyNodes(maxPolys),			&
		PolyNodes(maxPolys,maxPolyNodes),	& ! The nodes that define each polygon -- in clockwise direction
		PolySegs(maxPolys,maxPolyNodes)		  ! The segments that define each polygon -- clockwise from every node 
							  !      Note that the number of segments = number of nodes (numPolySegs = numPolyNodes)
	real*8						&
		PolyCenterX(maxPolys),			&
		PolyCenterY(maxPolys)

!	Crystal variables
	integer*4					&!
		numCrystals,				&! Number of crystals in grid = numPolys
		numCrystalNodes(maxCrystals),		&! Number of nodes in each crystal = numPolyNodes
		CrystalNodes(maxCrystals,maxPolyNodes),	&! Nodes for every crystal = PolyNodes -- clockwise from start
		CrystalNodeSegs(maxCrystals,6,2),	&! The segments on either side of a node for a specific crystal -- set in GridMaker routine
		CrystalMIFID(0:maxCrystals),		&! phase associated with the crystal. The integer index is the phaseMIF
		CrystalSegs(maxCrystals,maxPolyNodes)	 !Segments for each crystal = these are the same as the PolySegs around a Poly
	Character*8					&
		CrystalPhaseName(MaxCrystals)
	real*8						&
		CrystalCenterX(maxCrystals),		&! same as PolyCenterX
		CrystalCenterY(maxCrystals),		&! same as PolyCenterY
		CrystalArea(maxCrystals)		 ! Area of each crystal (poly) calculated in sub CrystalAreaRoutine
!	Phase variable (these get associated with every crystal)
	integer*4					&!
		numPhases,				&! number of different phases in the rock
		PhaseMIF(maxPhases),			&! mineral record number to use with Gibbs
		PhaseColorNumber(maxPhases),		&! hex color for AWE color for the phase
		numPhaseCrystals(maxPhases),		&! Number of crystals assigned to this phase
		PhaseCrystalNumber(maxPhases,maxCrystals) ! number of each crystal (polygon) to which this phase is assigned

	real*4						&
		PhaseCMYK(maxPhases,4)			 !CMYK color for Illustrator color for the phase

	Character*8					&
		PhaseName(maxPhases),			&! name of phase
		PhaseColorName(maxPhases)		 ! Name of color for phase

!	Node phase variables - for each node in the interior, there are 3 phases associated
!	We need to know these for the MDF reactions
!	The phase ID will be the MIF code so that Gibbs will know which phase it is
!	Note that ONLY the nodes with 3 or 2 different phases will be indexed here because only those can have reactions
	integer*4 					&
		numNodeReactionPhases(maxNodes),	& ! the number of different reacting phases at a node
		nodeReactionPhases(maxNodes,3),		& ! index of the reacting phases around the node 
							  ! if numNodeReactionPhases = 0 this is irrelevent
							  ! if numNodeReactionPhases = 3, then the indicies are 1, 2, 3
							  ! if numNodeReactionPhases=2 the indicies are either 1&2, 1&3 or 2&3
		numMDFRxnsNodes(maxNodes),		& ! number of MDF reactions at a node - calculated on start-up. If numMDFRxnsNodes =0 then an equilibrium model is used
		nodeSegTheSamePhases(maxNodes),		& ! This is the segment along which the 2 phases are the same for nodes where 2 phases are the same (set in ReadModelFile line 500 or so)
		nodeMIFID(maxNodes,3),			& ! the 3 phases (MIF ID) at the node
		nodeCrystalIndex(maxNodes,3)		 ! the crystal number (index) for each of the 3 crystals at the node
!	Segment phase variables. For each segment that has 2 adjoining phases that are different, a reaction can occur
!	Note that only segments that have 2 different phases will be indexed
	integer*4					&
		numSegReactionPhases(maxSegs),		& ! Number of phases that can react on each seg. Values are 0 or 2 
		numMDFRxnsSegs(maxSegs),		& ! number of MDF reactions at a seg - calculated on start-up. If numMDFRxnsSegss =0 then an equilibrium model is used
		segMIFID(maxSegs,2),			&! the 2 phases (MIF ID) along the segment
		segCrystalID(maxSegs,2)			 ! The 2 crystal IDs associated with a segment. Used in Sub GrowSeg routine to figure out which direction to grow the crystal along a seg

!	Phase compositions and amounts
	real*8								&
		nodePhaseComp(maxNodes,3,maxPhCo),			& ! Composition of each of 3 phases at a reaction node
		nodePhaseMoles(maxNodes,3),				&! Moles of each phase at a reaction node - cumulative
		nodePhaseMolesDelta(maxNodes,3),			&! Moles of each phase at a reaction node - Only the last set of calculations - temporary to use for crystal growth
		nodePhaseAffinity(maxNodes,3),				&! Affinity each phase at a reaction node from the last iteration
		pointPhaseComp(maxSegs,minPoints:maxPoints,2,maxPhCo),	&! Composition of each of 2 phases at a reaction segment
		nodeAreaDelta(maxNodes,3),				& ! the ∆area that corresponds to the ∆moles. Calculated in MolesToAreaRoutine
		pointPhaseMoles(maxSegs,minPoints:maxPoints,2),		&! Moles of each phase at a reaction segment point - cumulative
		pointPhaseMolesDelta(maxSegs,minPoints:maxPoints,2), 	&! Moles of each phase at a reaction segment point - Only the last set of calculations
		pointAreaDelta(maxSegs,minPoints:maxPoints,2) 		 ! the ∆area that corresponds to the ∆moles. Calculated in MolesToAreaRoutine
		

		

!	Composition variables
	character*4 ElName(maxElements)
	Integer*4						&
		numEl						 ! = number of elements
	Real*8							&
		dTime,						&! = time step
		InitialComp(maxElements),			&! = initial composition (read from file GB_Initial_Comp.txt)		
		nodeComp(maxNodes,maxElements),			&! = composition of GB at node
		pointComp(maxSegs,minPoints:maxPoints,maxElements) ! = composition of GB at each point

!	Plot composition variables
	integer*4 						&
		numNodesToPlot,					&! number of nodes to plot
		nodeToPlot(maxNodes),				&! index of nodes to plot
		numSegsToPlot,					&! number of segments to plot
		segToPlot(maxSegs),				&! index of segments to plot
		segToPlotBackwards(maxSegs)			 ! =0 to plot seg forward, =1 to plot seg backwards
	
!	Diffusion variables
	real*8 Dzero(maxElements),DactE(maxElements),DactVol(maxElements),DMatrix(maxElements,maxElements),	&
		DatTP(maxElements),ElCharge(maxElements)

!	Gradient color variables
	integer*4						&
		numPositiveColors,numNegativeColors,		&! number of colors
		positiveColors(20),negativeColors(20)		 ! Arrays for color values (in Hex)

!	Perturb Composition variables
	integer*4 nodeToPerturb				 ! this is the number of the node we will perturb
	Real*8						&
		PerturbComp(maxElements)		 ! = perturbation compositions 		

!	Model values
	integer*4					&
		totalCycles,				&! total number of model cycles (written to output file)
		numDiffIterations			 ! number of diffusion iterations in each model cycle

!	Misc
	integer*4 ComputeRoutine   	!2,3 or 4
	integer*4 OutputType   	! 1(no affinity) or 2(affinity=yes) 
	real*8	  newtonDamp,convergenceJoules,timeStep
	common /Diffuse_GB/numNodes,											&
		nodeX,nodeY,numNodeSegs,nodeNodeConnect,nodeSegConnect, 						&
		numSegs,segNodes,segX,segY,segLength,segCaptureExclusion,segGettingShorter,				&
		pointX,pointY,numPointStart,numPointEnd,								&
		dXgrid,dGridRemove,dGridAdd,										&
		nodeXOrig,nodeYOrig,											&
		segdXGrid,nodePointOnTop,nodePointNextTo,numCrystalsAtNode,						&
		numPolys,numPolyNodes,PolyNodes,PolySegs,PolyCenterX,PolyCenterY,					&
		numCrystals,numCrystalNodes,CrystalNodes,CrystalNodeSegs,CrystalMIFID,CrystalSegs,			&
		CrystalPhaseName,CrystalCenterX,CrystalCenterY,CrystalArea,						&
		numPhases,PhaseMIF,PhaseColorNumber,PhaseCMYK,numPhaseCrystals,PhaseCrystalNumber,PhaseName,		&
		PhaseColorName,												&
		nodeMIFID,nodeCrystalIndex,										&
		numNodeReactionPhases,nodeReactionPhases,numMDFRxnsNodes,numMDFRxnsSegs,nodeSegTheSamePhases,		&
		numSegReactionPhases,											&
		segMIFID,segCrystalID,											&
		nodePhaseComp,nodePhaseMoles,pointPhaseComp,pointPhaseMoles,nodePhaseAffinity,				&
		nodePhaseMolesDelta,nodeAreaDelta,									&
		pointPhaseMolesDelta,pointAreaDelta,									&
		numEl,ElName,												&
		dTime,InitialComp,nodeComp,pointComp,									&
		numNodesToPlot,nodeToPlot,numSegsToPlot,segToPlot,segToPlotBackwards,					&
		Dzero,DactE,DactVol,DMatrix,DatTP,ElCharge,nodeToPerturb,PerturbComp,					&
		numPositiveColors,numNegativeColors,positiveColors,negativeColors,					&
		totalCycles,numDiffIterations,ComputeRoutine,newtonDamp,convergenceJoules,outputType,timeStep

		real*4 yMinMaxSi(2),yMinMaxAl(2),yMinMaxFM(2)
		common /yMinMax/yMinMaxSi,yMinMaxAl,yMinMaxFM	

	integer*4 maxXY
	parameter (maxXY = 1000)
	real*8 x1(maxXY),y1(maxXY),x2(maxXY),y2(maxXY)	! = x and y values of segments read off of the svg file. Used only to make the grid
	real*8 minX,maxX,minY,maxY
	common /SVG_Stuff/ x1,y1,x2,y2,minX,maxX,minY,maxY



! Old common block

!	common /Diffuse_GB/numNodes,											&
!		nodeX,nodeY,numNodeSegs,nodeNodeConnect,nodeSegConnect,nodeSegEnd,					&
!		numSegs,segNodes,segX,segY,segLength,pointX,pointY,numPoints,dXgrid,					&
!		nodeXOrig,nodeYOrig,											&
!		segdXGrid,nodePointSeg,nodePointOnTop,nodePointNextTo,							&
!		numPolys,numPolyNodes,PolyNodes,PolySegs,PolyCenterX,PolyCenterY,					&
!		numCrystals,numCrystalNodes,CrystalNodes,CrystalNodeSegs,CrystalMIFID,CrystalSegs,			&
!		CrystalPhaseName,CrystalCenterX,CrystalCenterY,								&
!		numCrystalsToScale,CrystalsToScale,CrystalScaleFact,							&
!		CrystalArea,CrystalMoles,CrystalMolesNodes,CrystalMolesSegs,						&
!		numPhases,PhaseMIF,PhaseColorNumber,numPhaseCrystals,PhaseCrystalNumber,PhaseName,PhaseColorName,	&
!		numNodesWith3Phases,nodesWith3PhasesIndex,nodeMIFID,nodeCrystalIndex,					&
!		numNodesWith2Phases,nodesWith2PhasesIndex,								&
!		numSegsWith2Phases,segReactionPhases,segMIFID,segCrystalID,						&
!		nodePhaseComp,nodePhaseMoles,pointPhaseComp,pointPhaseMoles,						&
!		nodePhaseMolesDelta,pointPhaseMolesDelta,									&
!		numEl,ElName,dTime,InitialComp,nodeComp,pointComp,							&
!		numNodesToPlot,nodeToPlot,numSegsToPlot,segToPlot,segToPlotBackwards,					&
!		Dzero,DactE,DactVol,DMatrix,DatTP,ElCharge,nodeToPerturb,PerturbComp,					&
!		numPositiveColors,numNegativeColors,positiveColors,negativeColors,					&
!		totalCycles,numDiffIterations


