! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GB_RunModel_3(SegPlot1,SegPlot2,SegPlot3)
	USE AWE_INTERFACES
	implicit none
!	TYPE(AWE_Canvas) :: SegPlot(3)
	TYPE(AWE_Canvas) :: GridPlot
	TYPE(AWE_Canvas) :: SegPlot1
	TYPE(AWE_Canvas) :: SegPlot2
	TYPE(AWE_Canvas) :: SegPlot3
	TYPE(AWE_Canvas) :: MyWindow		! The plotting canvas

! *****************************************
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	include "AutoAffinity.inc"
	include "Gibbsfiles.inc"
! *****************************************
	integer*4 iTest,iNode,iSeg,iPointX,i,j,jj,numIter,numCycles,iunit,isteps,numSteps,iok,startStep,status,numCyclesInput
	integer*4 iNodeStart,iNodeEnd,ph1,ph2,ph3,debug,endstep,lastCycle,calcSegments,igo,L,M
	integer*4 k,iXl,NodeCaptureOn,MoveNodesandSegs,error,abort,segStart,segEnd,output,garnetMIFID
	character*16 date,time,zone,steptext
	integer timevalues(8)
	real hoursInit,minutesInit,secondsInit
	real*8 moleScaleFactor,lengthScaleFactor
	character*32 TPtext

	computeRoutine = 2
	totalCycles = 0
	timeStep = 0.1d0	! this is to calculate ∆Aff(target) = ∆Aff(old)*exp(-timeStep)
10	continue
	write(*,*)' ************************'
	write(*,*)' MDF Run Model routines - Incremental method'
	write(*,*)' 0 = return'
	write(*,*)' 1 = Run a model '
 	WRITE(*,*)'101= Poly Model Routine'
	write(*,*)' 3 = Call WriteAllNodesAndSegments (Save)'
	write(*,*)' 5 = SegmentPlotRoutines'
	WRITE(*,*)' 6 = Grid Plot routines'
	WRITE(*,*)' 7 = Phase composition Plot routines'
	WRITE(*,*)'10 = Test MDF calculations on nodes and segments'
	write(*,*)'11 = Run MDF routine on all nodes'
	write(*,*)'110 = Run MDF routine on all segments'
	write(*,*)'12 = Check mass balance on nodes and segments'
	write(*,*)'13 = Run Diffusion Routine'
	write(*,*)'14 = Plot segments'
	write(*,*)'15 = Test GB_Diffusion'
	write(*,*)'16 = Test GB_GrowNode (2 or 3)'
	write(*,*)'17 = Test GB_GrowSegs'
	write(*,*)'18 = Test RemoveSegPoints(iNode)'
	write(*,*)'19 = Test Node Capture'
	write(*,*)'20 = Open model input file'
	write(*,*)'21 = Open model output file'
	write(*,*)'22 = Calculate sum of moles around nodes and segments'
	write(*,*)'35 = Check Node information'
	write(*,*)'36 = Check Seg information'
	write(*,*)'37 = Check Crystal information'
	WRITE(*,*)'77 = Print assemblage'
	WRITE(*,*)'78 = Echo model output file to a new file (for checking)'
	
	read(*,*) itest
	select case(itest)
	case(0)
		return
	case(1)
		write(*,*)' Input number of model steps (each step will take several minutes)'
		write(*,*)' Note that  -- save and write out every step'
		read(*,*)numSteps
		if(numSteps.eq.0)go to 10
		write(*,*)' Input the starting step. = 1 if this is a new model, or = ## if it is a continuation'
		read(*,*)startStep
		if(startStep.ne.1)then
			call FSS_Alert('ALERT','Be sure you have opened model file and last output file before continuing')
			write(*,*)'OK to continue? 0 = NO -- ; 1 = OK'
			read(*,*)iOK
			if(iOK.eq.0)go to 10
			endif

		if(startStep.eq.1)then
			! select the base file name and write out the initial model info
			call FSS_Alert('ALERT',' Specify name of base output file')
			call OpenOutputFile(iok)
			if(iok.ne.0)go to 10
			isteps = 0
			write(steptext,55)isteps
55			format(I5)
			steptext = adjustL(steptext)
			ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
			open(74,FILE=ModelOutputFile,status = 'UNKNOWN')
			iunit = 74
				call WriteAllHeader(iunit)
				call WriteAllNodesAndSegments(iunit)
			close(iunit)

			else
			! We already have a base. Just open this base file name
			call FSS_Alert('ALERT',' Select name of original base output file')
			open(73,FILE='',status='OLD',iostat=status,action='WRITE')
			if(status.ne.0)then
				call FSS_Alert('Alert','Problem opening model file')
				go to 10
				endif
			inquire(73,NAME=ModelOutputFileBase)
			close(73)
			endif

		! Open log file
		if(startstep.eq.1)then
			write(steptext,55)startstep-1
			else
			write(steptext,55)startstep
			endif			
		steptext = adjustL(steptext)
		ModelLogFile = trim(ModelOutputFileBase)//'__'//trim(steptext)//'.log'
		open(75,FILE=ModelLogFile,status = 'UNKNOWN')

		ModelStuffFile = trim(ModelOutputFileBase)//'__'//trim(steptext)//'.stuff.txt'
		open(76,FILE=ModelStuffFile,status = 'UNKNOWN')



		write(*,*)'Input number of MDF + Diffusion cycles to run for each step'
		write(*,*)'Note that each cycle might take around 1 minute, so plan accordingly'
		read(*,*)numCyclesInput
		if(numCyclesInput.eq.0)go to 10
		write(*,*)' Input number of diffusion iterations (e.g. 50-5000)'
		write(*,*)' If you specify a number =100,000 then the GB composition will be averaged after reaction step '
		write(*,*)'      (e.g. infinite diffusion) '
		read(*,*)numDiffIterations

		lastCycle = 0
		if(startStep.eq.1)then
			lastCycle = lastCycle + 10
			endif
		if(numSteps.gt.10.and.startStep.eq.1)then
			lastCycle = lastCycle + numCyclesInput*(numSteps - 10)
			endif
		if(numSteps.gt.10.and.startStep.ne.1)then
			lastCycle = totalCycles + numCyclesInput*(numSteps)	! Total cycles is from the last run or read from the output file
			endif
		!lastCycle = numSteps*numCyclesInput

		moleScaleFactor = 100
		write(*,*)' Input mole scale factor (around 1-10-100-1000?)'
		read(*,*)moleScaleFactor

! 		write(*,*)'Input value for newtonDamp (0.01-1.0)'
! 		read(*,*)newtonDamp
		write(*,*)' Current time step = ',timeStep
		write(*,*)'Input value for timeStep '
		read(*,*)timeStep
					
		write(*,*)'Do you want to calculate segments? 0 = no; 1 = yes'
		read(*,*)calcSegments
		
		write(*,*)' Do you want to move nodes and segments? 0 = no; 1 = yes'
		read(*,*)MoveNodesandSegs
		
		write(*,*)'Do you want NodeCapture to be on? 0 = no; 1 = yes'
		read(*,*)NodeCaptureOn

		do i = 1,numPhMIF
			if(trim(phName(i)).eq.'Garnet')then
				garnetMIFID = i
				write(*,*)' GarnetMIFID = ',garnetMIFID
				endif
			end do

		debug = 0
		! write model parameters to the log file
! 			Model steps = 100
! 			10 iterations/step 
! 			Diffusion iterations = 1000
! 			moleScaleFactor = 100

		! log file
		write(75,*)'Model name = ',ModelOutputFileBase     
		write(75,*)'Model steps             = ',numSteps
		write(75,*)'Iterations/step         = ',numCyclesInput
		write(75,*)'Diffusion iterations    = ',numDiffIterations
		write(75,*)'Mole scale factor       = ',moleScaleFactor
		write(75,*)'Calculate segments(0-1) = ',calcSegments
		write(75,*)'MoveNodesandSegs (0-1)  = ',MoveNodesandSegs
		write(75,*)'NodeCaptureOn (0-1)     = ',NodeCaptureOn
		write(75,*)'ConvergenceJoules       = ',convergenceJoules
		write(75,*)'timeStep                = ',timeStep
		write(75,*)'dTime (Diffusion)       = ',dTime
		write(75,*)'Time/cycle (sec)        = ',dTime*numDiffIterations
		write(75,*)'Time/cycle (min)        = ',dTime*numDiffIterations/60
		write(75,*)'Time/cycle (hrs)        = ',dTime*numDiffIterations/(60*60)
		write(75,*)'Time/cycle (days)       = ',dTime*numDiffIterations/(60*60*24)
		write(75,*)'Time/cycle (yrs)        = ',dTime*numDiffIterations/(60*60*24*364)
		write(75,*)''
		write(75,*)''
		write(75,*)''
		
		! STUFF file
		write(76,*)'Model name = ',ModelOutputFileBase     
		write(76,*)'Model steps             = ',numSteps
		write(76,*)'Iterations/step         = ',numCyclesInput
		write(76,*)'Diffusion iterations    = ',numDiffIterations
		write(76,*)'Mole scale factor       = ',moleScaleFactor
		write(76,*)'Calculate segments(0-1) = ',calcSegments
		write(76,*)'MoveNodesandSegs (0-1)  = ',MoveNodesandSegs
		write(76,*)'NodeCaptureOn (0-1)     = ',NodeCaptureOn
		write(76,*)'ConvergenceJoules       = ',convergenceJoules
		write(76,*)'timeStep                = ',timeStep
		write(76,*)'dTime (Diffusion)       = ',dTime
		write(76,*)'Time/cycle (sec)        = ',dTime*numDiffIterations
		write(76,*)'Time/cycle (min)        = ',dTime*numDiffIterations/60
		write(76,*)'Time/cycle (hrs)        = ',dTime*numDiffIterations/(60*60)
		write(76,*)'Time/cycle (days)       = ',dTime*numDiffIterations/(60*60*24)
		write(76,*)'Time/cycle (yrs)        = ',dTime*numDiffIterations/(60*60*24*364)
		write(76,*)''
		write(76,*)''
		close(76)	! stuff file
		
		call OpenCanvasWindow(MyWindow)
		endstep = numSteps+startStep-1	
!		do isteps = startStep,numSteps+startStep-1
		do isteps = startStep,endstep

			write(75,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
			write(75,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
			write(75,*)'Step number = ',isteps
			write(12,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
			write(12,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
			write(12,*)'Step number = ',isteps

			call date_and_time(date,time,zone,timevalues)
			write(75,*)date,time,zone,timevalues
			!read(time,*)timeInit
			hoursInit   = timevalues(5)
			minutesInit = timevalues(6)
			secondsInit = timevalues(7)
			!		write(*,*)'hours, minutes, seconds ',hours,minutes,seconds
			if(isteps.le.100)then	
				numCycles = 1		! this should save every MDF+Diffusion iteration for the first step
				else
				numCycles = numCyclesInput	! after the first 10 steps, it saves every numCyclesInput
				endif
			do j = 1,numCycles
				totalCycles = totalCycles + 1
				call date_and_time(date,time,zone,timevalues)
				write(75,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
				write(75,*)' Step, Cycle number, totalCycles = ',iSteps,j,totalCycles
				TPtext = 'Steps       Cycles'
				call TextOnMyWindow(MyWindow,20,18,TPtext,16)
				write(TPtext,390)isteps,endstep,j,numCycles
!				write(TPtext,390)isteps,endstep,totalCycles,lastCycle
390				format(2I5,2x,2I5)
				call TextOnMyWindow(MyWindow,20,40,TPtext,16)
! 				read(time,*)timeNow
! 				hoursNow   = timevalues(5)
! 				minutesNow = timevalues(6)
! 				secondsNow = timevalues(7)
! 				timeNow = 60.*(hoursNow-hoursInit) + minutesNow-minutesInit + (secondsNow-secondsInit)/60.
! 				write(75,*)' dTime(min) ',timeNow
				write(75,*)'---------------'
				write(75,*)'Nodes'
				do iNode = 1,numNodes
					!write(75,*)iNode,numNodeReactionPhases(iNode)
					select case(numNodeReactionPhases(iNode))
						case(3)
							ph1 = nodeReactionPhases(iNode,1)
							ph2 = nodeReactionPhases(iNode,2)
							ph3 = nodeReactionPhases(iNode,3)
! 							if(numMDFRxnsNodes(iNode).gt.0)then
! 								call MDFNodesWith3Phases(iNode,3,ph1,ph2,ph3,0,0)
! 								else
! 								call EQNode(iNode,3,ph1,ph2,ph3,0,0)
! 								endif
							igo = 0
							call MDFNodesWith3Phases_3(iNode,3,ph1,ph2,ph3,0,0,igo)
							if(igo.eq.1)go to 10
						case(2)	
							ph1 = nodeReactionPhases(iNode,1)
							ph2 = nodeReactionPhases(iNode,2)
! 							if(numMDFRxnsNodes(iNode).gt.0)then
! 								call MDFNodesWith2Phases(iNode,2,ph1,ph2,0,0)
! 								else
! 								call EQNode(iNode,2,ph1,ph2,0,0,0)		! Ph3=0 because only 2 phases are used
! 								endif
							igo = 0
							call MDFNodesWith2Phases_3(iNode,2,ph1,ph2,0,0,igo)
							if(igo.eq.1)go to 10
						case default
						end select	
						
					end do	
				write(75,*)'---------------'
				write(75,*)'Segments'
					do iSeg = 1,numSegs
						if(calcSegments.eq.1)then
							if(numSegReactionPhases(iSeg).eq.2)then
								! the first and last points on each segment are on the nodes. 
								! get the GB node composition from the node calculations above
								! we don't bother setting the phase compositions at the first and last point because
								! they are calculated and saved with the node.
								! Only the GB comp is needed for the diffusion routine
								!do iPointX = 2,numPoints(iSeg)-1
								do iPointX = numPointStart(iSeg)+1,numPointEnd(iSeg)-1
									!write(75,*)iSeg,iPointX
									igo = 0
									call MDFPairs_3(iSeg,iPointX,0,0,igo)		! seg, point, output, affinity
									if(igo.eq.1)go to 10
									end do
! 								if(numMDFRxnsSegs(iSeg).gt.0)then
! 									do iPointX = numPointStart(iSeg)+1,numPointEnd(iSeg)-1
! 										!write(75,*)iSeg,iPointX
! 										call MDFPairs(iSeg,iPointX,0,0)		! seg, point, output, affinity
! 										end do
! 									else
! 									do iPointX = numPointStart(iSeg)+1,numPointEnd(iSeg)-1
! 										!write(75,*)iSeg,iPointX
! 										call EQSegment(iSeg,iPointX,0,0)		! seg, point, output, affinity
! 										end do
! 									endif
								endif
							endif
						! Set the grain boundary composition for the start and end of the segment
						iNodeStart = segNodes(iSeg,1)
						iNodeEnd   = segNodes(iSeg,2)
						do jj = 1,numPhCo(1)
							pointComp(iSeg,numPointStart(iSeg),jj) = nodeComp(iNodeStart,jj)
							pointComp(iSeg,numPointEnd(iSeg),jj)   = nodeComp(iNodeEnd,jj)
							end do						
						end do

				Call SetPointToNodePhaseComps(75)		! 75 is the log file

				write(75,*)'---------------'
				write(75,*)' Done with MDF'

				if(MoveNodesandSegs.eq.1)then
					write(75,*)'---------------'
					write(75,*)' Doing moving nodes and segments'
					do iNode = 1,numNodes
						select case(numNodeReactionPhases(iNode))
							case(3)
								!if(numMDFRxnsNodes(iNode).gt.0)then	! The only nodes with 0 MDFRxns are Qtz+H2O
													! Do not grow the nodes at these places -- this will simulate water leaving the system 
								call Grow3Nodes(iNode,moleScaleFactor,debug)
								!endif
							case(2)	
								!if(numMDFRxnsNodes(iNode).gt.0)then	! The only nodes with 0 MDFRxns are Qtz+H2O
													! Do not grow the nodes at these places -- this will simulate water leaving the system 
								call Grow2Nodes(iNode,moleScaleFactor,debug)
								!endif
							case default
							end select	
						end do	

					do iSeg = 1,numSegs
						if(numSegReactionPhases(iSeg).eq.2)then
							!if(segMIFID(iSeg,1).eq.GarnetMIFID.or.segMIFID(iSeg,2).eq.GarnetMIFID)cycle	! Don't move garnet segs or add segment points.
							call GrowSegs(GridPlot,iSeg,moleScaleFactor,debug)
							endif
						end do

					do iNode = 1,numNodes
						! the next line is designed to keep garnet with only nodes (no segment points)
						! Remove/add as desired
! 						if(nodeMIFID(iNode,1).eq.GarnetMIFID.or.nodeMIFID(iNode,2).eq.GarnetMIFID.or.   &
! 							nodeMIFID(iNode,3).eq.GarnetMIFID)cycle	! Don't move garnet segs or add segment points.
						call RemoveSegPoints(iNode,debug)
						end do			
				
					write(75,*)'---------------'
					write(75,*)' Done moving nodes and segments'
					endif

				if(numDiffIterations.lt.60000)then
					write(75,*)'---------------'
					write(75,*)' doing diffusion'
					error = 0
					do i = 1,numDiffIterations
						call DiffuseNodes(error)
						if(error.eq.1)then
							write(*,*)' Error returned from Subroutine DiffuseNodes'
							write(*,*)'Abort? 0 = no, 1 = yes'
							read(*,*)abort
							if(abort.eq.1)then
								close(75)
								go to 10
								endif
							endif
						call DiffuseSegments()
						do iSeg = 1,numSegs
						! Set the grain boundary composition for the start and end of the segment equal to the nodeComp
						! note that nodeComp has changed based on Sub DiffuseNodes
							iNodeStart = segNodes(iSeg,1)
							iNodeEnd   = segNodes(iSeg,2)
							do jj = 1,numPhCo(1)
								pointComp(iSeg,numPointStart(iSeg),jj) = nodeComp(iNodeStart,jj)
								pointComp(iSeg,numPointEnd(iSeg),jj)   = nodeComp(iNodeEnd,jj)
								end do						
							end do
	
						end do	! end diffusion iterations
	
					else
					write(75,*)'---------------'
					write(75,*)' doing Average Grain Boundary'
					call AverageGBComp()
					endif

				if(NodeCaptureOn.eq.1)then
					! See if 2 nodes are close together and flip, if needed
					call SegLengthCalculate()
					if (totalcycles.gt.10)call NodeCapture(75,1)		! 75 is the log file,0 means turn debug off
					endif
				end do		! end this cycle

			! Write this model step out to disk
			write(steptext,55)isteps
			steptext = adjustL(steptext)
			ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
			open(74,FILE=ModelOutputFile,status = 'UNKNOWN')
			iunit = 74
			call WriteAllHeader(iunit)
			call WriteAllNodesAndSegments(iunit)
			close(iunit)


			end do		! end this model step (istep)

		! Write out last entry to the log file and close
		call date_and_time(date,time,zone,timevalues)
		write(*,*)date,time,zone,timevalues
		close(75)
		close(76)
		go to 10

 	case(101)
 		call PolyModelRoutine()

	case(3)
		write(*,*)'Input unit number'
		read(*,*)iunit
     		call WriteAllHeader(iunit)
		call WriteAllNodesAndSegments(iunit)
		close(iunit)
	case(5)
		call SegPlotRoutines(SegPlot1,SegPlot2,SegPlot3)

	! ------------------ 6 --------------------
	case(6)
		call GridPlotRoutines(GridPlot)

	case(7)
		call PhCompPlotRoutines(SegPlot1,SegPlot2,SegPlot3)

!  ----PRINT------------------------------------------
	case(77)
		CALL PRINTT(1)
	case(78)
		! write out the problem to the output window -- just to check
		write(*,*)' '
		write(*,*)'*******************************************************'
		write(*,*)' Writing out problem to output window'
		write(*,*)'*******************************************************'
		write(*,*)' Input logical unit number to use'
		read(*,*)iUnit
!		iunit = 12
		call WriteAllHeader(iunit)
		call WriteAllNodesAndSegments(iunit)
	case(10)
		call TestMDFRoutines_3()

	case(11) 
		!computeRoutine = 2
! 		write(*,*)'Input value for newtonDamp (0.1 - 1.0)'
! 		read(*,*)newtonDamp
		write(*,*)' Current time step = ',timeStep
		write(*,*)'Input value for timeStep '
		read(*,*)timeStep
		write(*,*)'  '
		write(12,*)'$$$$$$$$$$$$$$$$$$$$$$$$$'
		write(12,*)'Nodes (iNode, numNodeReactionPhases)'
		do iNode = 1,numNodes
			write(12,*)iNode,numNodeReactionPhases(iNode)
			select case(numNodeReactionPhases(iNode))
				case(3)
					ph1 = nodeReactionPhases(iNode,1)
					ph2 = nodeReactionPhases(iNode,2)
					ph3 = nodeReactionPhases(iNode,3)
! 					if(numMDFRxnsNodes(iNode).gt.0)then
! !						call MDFNodesWith3Phases(iNode,3,ph1,ph2,ph3,1,0)
! 						call MDFNodesWith3Phases(iNode,3,ph1,ph2,ph3,0,0)
! 						else
! !						call EQNode(iNode,3,ph1,ph2,ph3,1,0)
! 						call EQNode(iNode,3,ph1,ph2,ph3,0,0)
! 						endif
					igo = 0
					call MDFNodesWith3Phases_3(iNode,3,ph1,ph2,ph3,0,0,igo)
					if(igo.eq.1)go to 10
				case(2)	
					ph1 = nodeReactionPhases(iNode,1)
					ph2 = nodeReactionPhases(iNode,2)
! 					if(numMDFRxnsNodes(iNode).gt.0)then
! !						call MDFNodesWith2Phases(iNode,2,ph1,ph2,1,0)
! 						call MDFNodesWith2Phases(iNode,2,ph1,ph2,0,0)
! 						else
! !						call EQNode(iNode,2,ph1,ph2,0,1,0)		! Ph3 = 0 because only 2 phases are used
! 						call EQNode(iNode,2,ph1,ph2,0,0,0)		! Ph3 = 0 because only 2 phases are used
! 						endif
					igo = 0
					call MDFNodesWith2Phases_3(iNode,2,ph1,ph2,0,0,igo)
					if(igo.eq.1)go to 10
				case default
				end select	
			end do	
		write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$'
		write(*,*)' Done'

	case(110) 

! 		write(*,*)'Input value for newtonDamp (0.1 - 1.0)'
! 		read(*,*)newtonDamp
		write(*,*)' Current time step = ',timeStep
		write(*,*)'Input value for timeStep '
		read(*,*)timeStep
		
		write(*,*)'  '
		write(*,*)' Input segment number to start'
		read(*,*)segStart
		write(*,*)' Input segment number to end. Maximum = ',numSegs
		read(*,*)segEnd
		write(*,*)' Do you want debugging output? 0 = no, 1 = yes '
		read(*,*)output
		write(12,*)'$$$$$$$$$$$$$$$$$$$$$$$$$'
		write(12,*)'Segments'
!		do iSeg = 1,numSegs
		do iSeg = segStart,segEnd
			if(numSegReactionPhases(iSeg).eq.2)then
				! the first and last points on each segment are on the nodes. 
				! get the GB node composition from the node calculations above
				! we don't bother setting the phase compositions at the first and last point because
				! they are calculated and saved with the node.
				! Only the GB comp is needed for the diffusion routine
				!do iPointX = 2,numPoints(iSeg)-1
				do iPointX = numPointStart(iSeg)+1,numPointEnd(iSeg)-1
					write(12,*)iSeg,iPointX
					igo = 0
					call MDFPairs_3(iSeg,iPointX,output,0,igo)		! seg, point, output, affinity
					if(igo.eq.1)go to 10
					end do
				else
				write(12,*)iSeg,'no rxn'
				endif
			! this assignment needs to be done even if there are no reaction phases along a segment
			! The end node might have changed
			iNodeStart = segNodes(iSeg,1)
			iNodeEnd   = segNodes(iSeg,2)
			do j = 1,numPhCo(1)
				pointComp(iSeg,numPointStart(iSeg),j) = nodeComp(iNodeStart,j)
				pointComp(iSeg,numPointEnd(iSeg),j)   = nodeComp(iNodeEnd,j)
! 					pointComp(iSeg,1,j)               = nodeComp(iNodeStart,j)
! 					pointComp(iSeg,numPoints(iSeg),j) = nodeComp(iNodeEnd,j)
				end do						

			end do
			
		write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$'
		write(*,*)' Done'



	case(12)
		write(*,*)'---------------'
		write(*,*)' Check mass balance on nodes and segments'
		write(*,*)' The sum of the moles/oxygen should = 0 at every reaction point'
		call CheckMassBalance(12)
		write(*,*)'---------------'
		write(*,*)' Done with checking mass balance'

		
	case(13)
		write(*,*)' Input number of iterations to run (0 to abort)'
		read(*,*)numIter
		if(numIter.eq.0)go to 10
		do i = 1,numIter
			call DiffuseNodes(error)
			call DiffuseSegments()
			do iSeg = 1,numSegs
			! Set the grain boundary composition for the start and end of the segment equal to the nodeComp
			! note that nodeComp has changed based on Sub DiffuseNodes
				iNodeStart = segNodes(iSeg,1)
				iNodeEnd   = segNodes(iSeg,2)
				do jj = 1,numPhCo(1)
					pointComp(iSeg,numPointStart(iSeg),jj) = nodeComp(iNodeStart,jj)
					pointComp(iSeg,numPointEnd(iSeg),jj)   = nodeComp(iNodeEnd,jj)
					end do						
				end do

			end do
	case(14)
		Call PlotSegments(SegPlot1,SegPlot2,SegPlot3)
	case(15)
		call GB_Diffuse_Hex()
	case(16)
160		continue
		write(*,*)'   '
		write(*,*)'   '
		write(*,*)'********************************'
		write(*,*)'Input node to move'
		read(*,*)iNode
		if(iNode.eq.0)go to 10
		write(*,*)' Input mole scale factor (around 100)'
		read(*,*)moleScaleFactor
		debug = 1		
		select case(numNodeReactionPhases(iNode))
			case(3)
				call Grow3Nodes(iNode,moleScaleFactor,debug)
			case(2)		
				call Grow2Nodes(iNode,moleScaleFactor,debug)
			case default
				call fss_alert('Alert',' The number of reaction phases is not 2 or 3. Try again')
				go to 160
			end select
		go to 160
	case(17)
170		continue
		write(*,*)'Input Seg to move'
		read(*,*)iSeg
		if(iSeg.eq.0)go to 10
		write(*,*)' Input mole scale factor (around 10-100-1000)'
		read(*,*)moleScaleFactor
		debug = 1
		call GrowSegs(GridPlot,iSeg,moleScaleFactor,debug)
		go to 170
	case(18)
180		continue
		write(*,*)'Input node to check'
		read(*,*)iNode
		if(iNode.eq.0)then
			debug = 0
			go to 10
			endif
! 		do iNode = 1,numNodes
		debug = 1
		call RemoveSegPoints(iNode,debug)
! 			end do			
		go to 180
	case(19)
		call SegLengthCalculate()
! 		do iSeg = 1,numSegs
! 			write(12,*)iSeg,segLength(iSeg),segGettingShorter(iSeg)
! 			end do		
		write(*,*)' Input segment for capture'
		read(*,*)iSeg
		if(iSeg.eq.0)go to 10
		write(*,*)'SegLength = ',segLength(iSeg),segGettingShorter(iSeg)
		write(*,*)'Input scale factor '
		read(*,*)lengthScaleFactor
		segGettingShorter(iSeg) = .TRUE.
		segLength(iSeg) = lengthScaleFactor*dXgrid
		call NodeCapture(75,1)			! 1 = debug is ON
	case(20)
		call GibbsBegin()
	case(21)
		call OpenModelOutputFile(0,1)	!0 means file is not already open; 1 means read the node and seg XY coordinates
	case(22)
		call PhaseMassSum()	
	case(35)

		call CheckNodeInformation()

	case(36)
36		continue
		call SegLengthCalculate()
		write(*,*)' '
		write(*,*)'dXgrid = ',dXgrid
		write(*,*)'  '
		write(*,*)'  '
		write(*,*)'Input Segment to check. 0 to exit'
		read(*,*)iSeg
		if(iSeg.eq.0)go to 10
! 		write(*,*)'Input point to check'
! 		read(*,*)iPoint
		write(*,*)'Seg nodes            ',segNodes(iSeg,1),segNodes(iSeg,2)
		write(*,*)'Start,end            ',numPointStart(iSeg),numPointEnd(iSeg)
		write(*,*)'Start XY             ',segX(iSeg,1),segY(iSeg,1)
		write(*,*)'End   XY             ',segX(iSeg,2),segY(iSeg,2)
		write(*,*)'Seg length           ',segLength(iSeg)
		write(*,*)'Seg getting shorter  ',segGettingShorter(iSeg)
		write(*,*)'numSegReactionPhases ',numSegReactionPhases(iSeg)
		write(*,*)'numMDFRxnsSegs       ',numMDFRxnsSegs(iSeg)
		write(*,*)'segMIFIDs            ',segMIFID(iSeg,1),segMIFID(iSeg,2)
		write(*,*)'segCrystalID         ',segCrystalID(iSeg,1),segCrystalID(iSeg,2)
		write(*,*)'segCrystalID         ',segCrystalID(iSeg,1),segCrystalID(iSeg,2)

		go to 36
	case(37)
37		continue
		write(*,*)'  '
		write(*,*)'  '
		write(*,*)'Input Crystal to check. 0 to exit'
		read(*,*)iXl
		if(iXl.eq.0)go to 10
			write(*,*)'  Node     seg (clockwise)'
			do k = 1,numCrystalNodes(iXl) 
				iSeg = CrystalSegs(iXl,k)
				iNode = CrystalNodes(iXl,k)
				write(*,*)iNode,iSeg
 				do L = 1,3
 					if(iSeg.eq.nodeSegConnect(iNode,L))then		! this is the segment we want
 				!		if(nodeSegEnd(iNode,L).eq.1)then	! we start at the beginning of the segment
						if(segNodes(iSeg,1).eq.iNode)then
 							do m = numPointStart(iSeg),numPointEnd(iSeg)
								write(*,370)m,pointX(iSeg,m),pointY(iSeg,m)
								end do
370								format(10x,I8,2F10.5)
							else				! we start at the end of the segment
							do m = numPointEnd(iSeg),numPointStart(iSeg),-1
								write(*,370)m,pointX(iSeg,m),pointY(iSeg,m)
								end do
							endif															
						endif
					end do				
				end do
		go to 37

	case default
		call FSS_Alert('ALERT',' You chose poorly')
		
	end select
	go to 10

	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE TestMDFRoutines_3()
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	include "AutoAffinity.inc"
! *****************************************
	integer*4 iTest,iNode,iSeg,iPointX,ph1,ph2,ph3,igo


! 	write(*,*)'Input value for newtonDamp (0.1 - 1.0)'
! 	read(*,*)newtonDamp
		write(*,*)' Current time step = ',timeStep
		write(*,*)'Input value for timeStep '
		read(*,*)timeStep

10	continue

	write(*,*)' 0 = return'
	write(*,*)' 1 = select node to test'
	write(*,*)' 2 = select segment, point to test'
	read(*,*) itest

! 	write(*,*)'Select Compute routine to use: 2, 3 or 4'
! 	read(*,*)ComputeRoutine

	select case(itest)
	case(0)
		return
	case(1) 
		write(*,*)'Input node number to test'
		read(*,*)iNode
		if(iNode.eq.0)go to 10
		if(iNode.lt.0.or.iNode.gt.numNodes)then
			call FSS_Alert('ALERT',' Not a valid node number')
			endif			
		select case (numNodeReactionPhases(iNode))
		case(0,1)
			call FSS_Alert('ALERT',' The node has no reaction phases')
		case(2) 
			ph1 = nodeReactionPhases(iNode,1)
			ph2 = nodeReactionPhases(iNode,2)
			!Call AWE_SetOutput
! 			if(numMDFRxnsNodes(iNode).gt.0)then
! 				call MDFNodesWith2Phases(iNode,2,ph1,ph2,1,0)
! 				else
! 				call EQNode(iNode,2,ph1,ph2,0,1,0)		! Ph3 = 0 because only 2 phases are used
! 				endif
			write(12,*)' '
			write(12,*)' '
			write(12,*)' '
			write(12,*)' '
			write(12,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ '
			write(12,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ '
			write(12,*)' '
			write(12,*)' Node = ',iNode
			write(12,*)' '
			write(12,*)' Calling MDFNodesWith2Phases '
			igo = 0
			call MDFNodesWith2Phases_3(iNode,2,ph1,ph2,1,0,igo)
			if(igo.eq.1)go to 10
			call CalcAffinity(1)
			go to 10

		case(3)
			ph1 = nodeReactionPhases(iNode,1)
			ph2 = nodeReactionPhases(iNode,2)
			ph3 = nodeReactionPhases(iNode,3)
			!Call AWE_SetOutput
! 			if(numMDFRxnsNodes(iNode).gt.0)then
! 				call MDFNodesWith3Phases(iNode,3,ph1,ph2,ph3,1,0)
! 				else
! 				call EQNode(iNode,3,ph1,ph2,ph3,1,0)
! 				endif
			write(12,*)' '
			write(12,*)' '
			write(12,*)' '
			write(12,*)' '
			write(12,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ '
			write(12,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ '
			write(12,*)' '
			write(12,*)' Node = ',iNode
			write(12,*)' '
			write(12,*)' Calling MDFNodesWith3Phases '
			igo = 0
			call MDFNodesWith3Phases_3(iNode,3,ph1,ph2,ph3,1,0,igo)
			if(igo.eq.1)go to 10
			call CalcAffinity(1)
			go to 10
						
		case default
			call FSS_Alert('ALERT',' No such node')
		end select
	case(2)
		write(*,*)'Input segment to test'
		read(*,*)iSeg
!		write(*,*)' Number of points = ',numPoints(iSeg)
		write(*,*)' Start point = ',numPointStart(iSeg)
		write(*,*)' End   point = ',numPointEnd(iSeg)
		write(*,*)'Input point number to test'
		read(*,*)iPointX
		if(iPointX.gt.numPointStart(iSeg).and.iPointX.lt.numPointEnd(iSeg))then
			! this point should be OK
			if(numSegReactionPhases(iSeg).eq.2)then
				!Call AWE_SetOutput
! 				if(numMDFRxnsSegs(iSeg).gt.0)then
! 					call MDFPairs(iSeg,iPointX,1,0)		! seg, point, output, affinity
! 					else
! 					call EQSegment(iSeg,iPointX,1,0)		! seg, point, output, affinity
! 					endif
			write(12,*)' '
			write(12,*)' '
			write(12,*)' '
			write(12,*)' '
			write(12,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ '
			write(12,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ '
			write(12,*)'Segment, point  = ',iSeg,iPointX
			write(12,*)' '
			write(12,*)' Calling MDFPairs '
				igo = 0
				call MDFPairs_3(iSeg,iPointX,1,0,igo)		! seg, point, output, affinity
				if(igo.eq.1)go to 10
				call CalcAffinity(1)
				else
				call FSS_Alert('ALERT',' Not a valid seg with 2 phases')
				endif
			else
			call FSS_Alert('ALERT',' Not a valid point on this seg')
			endif			

		go to 10
	
	case default
		call FSS_Alert('ALERT','You chose poorly....')
	end select
	

	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE MDFNodesWith3Phases_3(iNodeTemp,numPhT,ph1,ph2,ph3,output,initializeOnly,igo)
	use MatrixArrays
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
!	include "Solute.inc"
	include "Diffuse_GB_Hex.inc"
	include "AutoAffinity.inc"
	include "GB_Tangent.inc"
! *****************************************
	integer*4 iNode,i,j,k,k1,k2,k3,k4,izero,icount,iNodeTemp,kk,MIFID,igo
	real*8 deltaMolesdX
	real*8 dF3_dDelta2,dF3_dDelta3,dF4_dDelta2,dF4_dDelta3,		&
			deltaAff3_old,deltaAff3_new,deltaAff4_old,deltaAff4_new,F3o,F4o
	real*8 deta,GBmoles,deltaMolesX(4),GBChangeX(10),deltaAff3_Target,deltaAff4_Target
	integer*4 nsolve,numEq,ier,output,initializeOnly,numPhT,ph1,ph2,ph3
!	if output = 1 then print
!	if output = 0 then don't print


	iNode = iNodeTemp
3	continue		

	numPh = numPhT + 1		! should be 4

	asmCurrent(1) = 1			! always the grain boundary
	asmCurrent(2) = nodeMIFID(iNode,ph1)
	asmCurrent(3) = nodeMIFID(iNode,ph2)
	asmCurrent(4) = nodeMIFID(iNode,ph3)

	k1 = asmCurrent(1)
	k2 = asmCurrent(2)
	k3 = asmCurrent(3)
	k4 = asmCurrent(4)
	if(output.eq.1)then
		write(12,*)'  '
		write(12,*)'  '
		write(12,*)'  '
		write(12,*)'********************************************************'
		write(12,*)' Current mineral assemblage'
		do k = 1,numPh
			WRITE(12,1533)MINREC(asmCurrent(K)),PHNAME(asmCurrent(K))
1533  			FORMAT(3(I4,2X,A8,2x))
			end do
		endif	
!	Set the composition of the GB into the XPhCo array for use in Sub CalculateGBChange and ParallelToTangent
	do j = 1,numPhCo(1)		! grain boundary
		xPhCo(1,j) = nodeComp(iNode,j)		! numPhCo(1) should be every element
		end do
	do i = 1,3		! do the 3 solid phases
		k = asmCurrent(i+1)
		kk = nodeReactionPhases(iNode,i)	!Note that the index of phases around the node isn't 1,2,3 necessarily
		select case(fractl(k))
		case(0,1)	! phases are not fractionating or are growing and fractionating
			do j = 1,numPhCo(k)
				xPhCo(k,j) = nodePhaseComp(iNode,kk,j)
				end do
		case(2)		! phase is being consumed - use fixed composition (i.e. chlorite)
			do j = 1,numPhCo(k)
				xPhCo(k,j) = xPhCoInitial(k,j)
				end do
		case default
		end select		
		end do

	CALL AdjustAsm(1)		! we changed assemblage. Adjust volumes, moles etc.


	!store the initial grain boundary composition
	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		GBComp(1,j) = xPhCo(1,j)
		end do


	!Initial calculation of affinity
	if(output.eq.1)then
		write(12,*)'-------------'
		write(12,*)' Initial call to CalcAffinity to get chemical potentials and affinities'
		call CalcAffinity(1)	! print this initial calculation
		write(12,*)' minrec  Ph            Aff       mmol '
		write(12,592)minrec(k2),PhName(k2),Aff(k2,1),1000.0d0*mp0(k2)
		write(12,592)minrec(k3),PhName(k3),Aff(k3,1),1000.0d0*mp0(k3)
		write(12,592)minrec(k4),PhName(k4),Aff(k4,1),1000.0d0*mp0(k4)
592		Format(3(I5,2x,A8,F12.1,E15.5))	
		else
		call CalcAffinity(0)	! Print
		!call CalcAffinity(0)	! do not print this initial calculation
		endif

	call Save_and_Reset_Comps(1)	! save the current phase compositions

	do i = 1,3
		k = asmCurrent(i+1)
		call paralleltotangent(k,0,izero)
		if(izero.gt.0)then
			write(12,*)' In routine MDFNodesWith3Phases_3. Line 1013'
			write(12,*)' Node = ',iNode
			write(12,*)' i =    ',i
			write(12,*)' k =    ',k
			write(12,*)' asmCurrent(j) = ',(asmCurrent(j),j=2,4)
			endif
		end do
	if(output.eq.1)then
		write(12,*)'                 gDifferenceAU      Phase Comp '
		write(12,81)phName(k2),gDifferenceAU(k2),(xPhCo(k2,j),j=1,numPhCo(K2))
		write(12,81)phName(k3),gDifferenceAU(k3),(xPhCo(k3,j),j=1,numPhCo(K3))
		write(12,81)phName(k4),gDifferenceAU(k4),(xPhCo(k4,j),j=1,numPhCo(K4))
		endif
	Aff(k2,1) = gDifferenceAU(k2)
	Aff(k3,1) = gDifferenceAU(k3)
	Aff(k4,1) = gDifferenceAU(k4)

	if(initializeOnly.eq.1)go to 99		! this call was only to establish starting compositions for all nodes

	
	MDFmoles(1,2) = mp0(k2)
	MDFmoles(1,3) = mp0(k3)
	MDFmoles(1,4) = mp0(k4)
	MDFAff(1,2) = Aff(k2,1)
	MDFAff(1,3) = Aff(k3,1)
	MDFAff(1,4) = Aff(k4,1)

! 	icount = 0

	!deltaMolesdX = 1.0d-4		! increment for numerical derivatives = .001 mmoles
	deltaMolesdX = 1.0d-5		! increment for numerical derivatives = .001 mmoles
	deltaAff3_old = (Aff(k3,1) - Aff(k2,1))	!This is the Aff with the original set of GB compositions
	deltaAff4_old = (Aff(k4,1) - Aff(k2,1))	!This is the Aff with the original set of GB compositions


!  100	continue

! 	icount = icount + 1
! 	if(output.eq.1)then
! 		write(12,*)' ----------------------------------'
! 		write(12,*)' ----------------------------------'
! 		write(12,*)' icount = ',icount
! 		endif


! 	if(icount.gt.100)then
! 		write(*,*)' No convergence in MDFNodesWith3Phases. iNode = ',iNode
! 		write(75,*)' No convergence in MDFNodesWith3Phases. iNode = ',iNode
! 		write(*,*)' 0 = continue'
! 		write(*,*)' 1 = abort'
! 		read(*,*)igo
! 		if(igo.eq.1)return
! 		go to 599	! exit routine
! 		endif

	

!	Calculate the derivative

!	X2 derivative
	deltaMolesX(2) = deltaMolesdX		! scaling for ∆moles for determining derivative of deltaAffinity
	deltaMolesX(3) = 0.0d0			! scaling for ∆moles for determining derivative of deltaAffinity
!	deltaMolesdX4 = -(deltaMolesdX2*numOxygens(k2) + deltaMolesdX3*numOxygens(k3))/numOxygens(k4)	! this is oxygen balance constraint
	deltaMolesX(4) = -(deltaMolesX(2) + deltaMolesX(3))	! this is oxygen balance constraint

	! now calculate the change in GB composition based on the above

	call Save_and_Reset_Comps(2)	! reset the current phase compositions after call to ParallelToTangent
	Call CalculateGBChange(deltaMolesX,GBChangeX)

	if(output.eq.1)then
		write(12,*)' '
		write(12,*)' Calculating derivatives '
		write(12,*)' dX2 derivative '
		write(12,*)' deltaMolesX ',deltaMolesX(2),deltaMolesX(3),deltaMolesX(4)
		write(12,*)' GB CHANGE   ',(GBChangeX(j),j=1,numPhCo(1))	
		endif

	! Now we have the composition of the grain boundary for the small increment of moles (deltaMolesdX)
	! Calculate the affinities
	do i = 1,3
		k = asmCurrent(i+1)
		call paralleltotangent(k,0,izero)
		if(izero.gt.0)then
			write(12,*)' In routine MDFNodesWith3Phases_3. Line 1097'
			write(12,*)' Node = ',iNode
			write(12,*)' i =    ',i
			write(12,*)' k =    ',k
			write(12,*)' asmCurrent(j) = ',(asmCurrent(j),j=2,4)
			endif

		end do
	if(output.eq.1)then
		write(12,*)'                 gDifferenceAU      Phase Comp '
		write(12,81)phName(k2),gDifferenceAU(k2),(xPhCo(k2,j),j=1,numPhCo(K2))
		write(12,81)phName(k3),gDifferenceAU(k3),(xPhCo(k3,j),j=1,numPhCo(K3))
		write(12,81)phName(k4),gDifferenceAU(k4),(xPhCo(k4,j),j=1,numPhCo(K4))
	81	format(A10,T15,15F15.5)
		endif

	Aff(k2,1) = gDifferenceAU(k2)
	Aff(k3,1) = gDifferenceAU(k3)
	Aff(k4,1) = gDifferenceAU(k4)

	deltaAff3_new = (Aff(k3,1) - Aff(k2,1))	!This is the change in ∆Aff with a small change in moles -- to get slope
	deltaAff4_new = (Aff(k4,1) - Aff(k2,1))	!This is the change in ∆Aff with a small change in moles -- to get slope
	! calculate the derivatives with respect to dX2
	dF3_dDelta2 = (deltaAff3_new - deltaAff3_old)/deltaMolesX(2)
	dF4_dDelta2 = (deltaAff4_new - deltaAff4_old)/deltaMolesX(2)


!	X3 derivative
	deltaMolesX(2) = 0.0d0		
	deltaMolesX(3) =  deltaMolesdX		! scaling for ∆moles for determining derivative of deltaAffinity
! 	deltaMolesdX4 = -(deltaMolesdX2*numOxygens(k2) + deltaMolesdX3*numOxygens(k3))/numOxygens(k4)	! this is oxygen balance constraint
	deltaMolesX(4) = -(deltaMolesX(2) + deltaMolesX(3))	! this is oxygen balance constraint

	! now calculate the change in GB composition based on the above
	call Save_and_Reset_Comps(2)	! reset the current phase compositions after call to ParallelToTangent
	Call CalculateGBChange(deltaMolesX,GBChangeX)

	if(output.eq.1)then
		write(12,*)' dX3 derivative '
		write(12,*)' deltaMolesX ',deltaMolesX(2),deltaMolesX(3),deltaMolesX(4)
		write(12,*)' GB CHANGE   ',(GBChangeX(j),j=1,numPhCo(1))	
		endif

	! Now we have the composition of the grain boundary for the small increment of moles (deltaMolesdX)
	! Calculate the affinities
	do i = 1,3
		k = asmCurrent(i+1)
		call paralleltotangent(k,0,izero)
		if(izero.gt.0)then
			write(12,*)' In routine MDFNodesWith3Phases_3. Line 1147'
			write(12,*)' Node = ',iNode
			write(12,*)' i = ',i
			write(12,*)' k =   ',k
			write(12,*)' asmCurrent(j) = ',(asmCurrent(j),j=2,4)
			endif

		end do
	if(output.eq.1)then
		write(12,*)'                 gDifferenceAU      Phase Comp '
		write(12,81)phName(k2),gDifferenceAU(k2),(xPhCo(k2,j),j=1,numPhCo(K2))
		write(12,81)phName(k3),gDifferenceAU(k3),(xPhCo(k3,j),j=1,numPhCo(K3))
		write(12,81)phName(k4),gDifferenceAU(k4),(xPhCo(k4,j),j=1,numPhCo(K4))
		endif
	Aff(k2,1) = gDifferenceAU(k2)
	Aff(k3,1) = gDifferenceAU(k3)
	Aff(k4,1) = gDifferenceAU(k4)

	deltaAff3_new = (Aff(k3,1) - Aff(k2,1))	!This is the change in ∆Aff with a small change in moles -- to get slope
	deltaAff4_new = (Aff(k4,1) - Aff(k2,1))	!This is the change in ∆Aff with a small change in moles -- to get slope
	! calculate the derivatives with respect to dX3
	dF3_dDelta3 = (deltaAff3_new - deltaAff3_old)/deltaMolesX(3)
	dF4_dDelta3 = (deltaAff4_new - deltaAff4_old)/deltaMolesX(3)


	! Calculate the functions to match
	! These are the change in affinities calculated from the time step (∆timeStep)
	deltaAff3_Target = deltaAff3_old*exp(-timestep)
	deltaAff4_Target = deltaAff4_old*exp(-timestep)

! 	F3o = deltaAff3_old			! F0 functions
! 	F4o = deltaAff4_old

	F3o = deltaAff3_Target - deltaAff3_old			! F0 functions
	F4o = deltaAff4_Target - deltaAff4_old

	!solve for ∆moles2 and ∆moles3 using Newton's method
! 	deltaMoles = -deltaAffinityOld/dF_dDelta

	A(1,1) = dF3_dDelta2
	A(1,2) = dF3_dDelta3
! 	A(1,3) = -F3o
	A(1,3) = F3o
	A(2,1) = dF4_dDelta2
	A(2,2) = dF4_dDelta3
! 	A(2,3) = -F4o
	A(2,3) = F4o

	do j=1,3
		do i=1,2
			AA(i,j) = A(i,j)
			end do
		end do
	!Find solution
	nsolve = 3
	DETA=0.D0
	IER=0
	numEq = 2
	CALL REDUCE (numEq,NSOLVE,DETA,IER)
	IF (IER.GT.0) then
		write(*,*)' Matrix singular in Subroutine MDFNodesWith3Phases in file GB_MDF_Routines.f90: iNode= ',iNode
		write(75,*)' Matrix singular in Subroutine MDFNodesWith3Phases in file GB_MDF_Routines.f90: iNode= ',iNode
! 		WRITE(75,*)' ************ ERROR **************************'
! 		WRITE(75,*)' Matrix failed to invert in SUBROUTINE REDUCE'
! 		write(75,*)' We are in Subroutine MDFNodesWith3Phases in file GB_MDF_Routines.f90 '
! 		write(75,*)'iNode= ',iNode		
		write(*,*)' 0 = continue'
		write(*,*)' 1 = abort'
		write(*,*)' 2 = Try to debug'
		read(*,*)igo
		select case (igo)
			case(0)
				izero = 0
			case(1)
				izero = 1
				return
			case(2)
				output = 1
				write(12,*)' '
				write(12,*)' '
				write(12,*)' '
				write(12,*)' ^^&^&^&^&^&^&^&^&^&^&&^&^&^&^&^&^&^^&'
				write(12,*)' ^^&^&^&^&^&^&^&^&^&^&&^&^&^&^&^&^&^^&'
				write(12,*)' Starting debug routine'
				write(*,*)'iNode = ',iNode
				go to 3		
			case default
				call fss_alert('Alert','You chose poorly')
			end select
		endif
	if(output.eq.1)then
		write(12,*)' ----'
		write(12,*)'counter =  ',icount
		write(12,*)' '
		write(12,*)'A matrix'
		write(12,*)'    dF_dDelta2     dF_dDelta3     -Fo'
		do i = 1,2
			write(12,1580)(AA(i,j),j=1,3)	! note: This is the original matrix. Matrix "A" gets modified in Sub Reduce
1580			format(50F15.5)
			end do
		write(12,*)' '
		write(12,*)'Solution xx'
		write(12,*)'    deltaMolesX(2)       deltaMolesX(3) '
		write(12,1581)(xx(j,1),j=1,numEq)
1581		format(2E15.5)
		write(12,*)' ----'
		endif

! 	deltaMolesX(2) = xx(1,1)*newtonDamp 		
! 	deltaMolesX(3) = xx(2,1)*newtonDamp
	deltaMolesX(2) = xx(1,1)
	deltaMolesX(3) = xx(2,1)
! 	deltaMolesdX4 = -(deltaMolesdX2*numOxygens(k2) + deltaMolesdX3*numOxygens(k3))/numOxygens(k4)	! this is oxygen balance constraint
	deltaMolesX(4) = -(deltaMolesX(2) + deltaMolesX(3))	! this is oxygen balance constraint

	mp0(k2) = mp0(k2) + deltaMolesX(2)
	mp0(k3) = mp0(k3) + deltaMolesX(3)
	mp0(k4) = mp0(k4) + deltaMolesX(4)
	
	! now calculate the change in GB composition based on the above root
	call Save_and_Reset_Comps(2)	!reset the compositions
	Call CalculateGBChange(deltaMolesX,GBChangeX)
	if(output.eq.1)then
		write(12,*)' '
		write(12,*)' After Newton''s method solution '
		write(12,*)' deltaMolesX ',deltaMolesX(2),deltaMolesX(3),deltaMolesX(4)
		write(12,*)' GB CHANGE   ',(GBChangeX(j),j=1,numPhCo(1))	
		endif

	! Now we have the composition of the grain boundary for the small increment of moles (deltaMolesdX)
	! Calculate the affinities
	do i = 1,3
		k = asmCurrent(i+1)
		call paralleltotangent(k,0,izero)
		if(izero.gt.0)then
			write(12,*)' In routine MDFNodesWith3Phases_3. Line 1267'
			write(12,*)' Node = ',iNode
			write(12,*)' i = ',i
			write(12,*)' k =   ',k
			write(12,*)' asmCurrent(j) = ',(asmCurrent(j),j=2,4)
			endif
		end do
	if(output.eq.1)then
		write(12,*)'                 gDifferenceAU      Phase Comp '
		write(12,81)phName(k2),gDifferenceAU(k2),(xPhCo(k2,j),j=1,numPhCo(K2))
		write(12,81)phName(k3),gDifferenceAU(k3),(xPhCo(k3,j),j=1,numPhCo(K3))
		write(12,81)phName(k4),gDifferenceAU(k4),(xPhCo(k4,j),j=1,numPhCo(K4))
		endif
	Aff(k2,1) = gDifferenceAU(k2)
	Aff(k3,1) = gDifferenceAU(k3)
	Aff(k4,1) = gDifferenceAU(k4)
! 	call CalcAffinity(1)	! Print
! 	call CalcAffinity(0)
	deltaAff3_new = (Aff(k3,1) - Aff(k2,1))	!This is the change in ∆Aff with a small change in moles -- to get slope
	deltaAff4_new = (Aff(k4,1) - Aff(k2,1))	!This is the change in ∆Aff with a small change in moles -- to get slope

	if(output.eq.1)then
		write(12,*)' DeltaAffNew       ',deltaAff3_new,deltaAff4_new	
		write(12,*)'                      Ph           Delta(mmol)          mmol          Aff '
		Write(12,562)k2,PhName(k2),1000.0d0*deltaMolesX(2),1000.0d0*mp0(k2),Aff(k2,1)
		write(12,562)k3,PhName(k3),1000.0d0*deltaMolesX(3),1000.0d0*mp0(k3),Aff(k3,1)
		write(12,562)k4,PhName(k4),1000.0d0*deltaMolesX(4),1000.0d0*mp0(k4),Aff(k4,1)
562		format(' DeltaM(mmol) ',I5,2x,A8,2x,3(E15.5))
		write(12,*)'  Aff(3-2)          Aff(4-2)  (delta Affinities)'
		write(12,*) Aff(k3,1) - Aff(k2,1),Aff(k4,1) - Aff(k2,1)
		write(12,*)' '
		Write(12,*)'---GB composition changes --'
		write(12,501)(PhCoName(1,j),j=1,numPhCo(1))
501		format(10A15)
		write(12,502)(GBChangeX(j),j=1,numPhCo(1))
		write(12,502)(xPhCo(1,j),j=1,numPhCo(1))
502		format(10E15.4)
		endif



			
	!if we get here, then we need to stop growing
599		continue

!	Store the final moles
	MDFmoles(2,2) = mp0(k2)
	MDFmoles(2,3) = mp0(k3)
	MDFmoles(2,4) = mp0(k4)
	MDFAff(2,2) = Aff(k2,1)
	MDFAff(2,3) = Aff(k3,1)
	MDFAff(2,4) = Aff(k4,1)

!	store the final grain boundary composition
	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		GBComp(2,j) = xPhCo(1,j)
		end do

1500		continue	! loop on every pair		
	
!	write out initial and final GB compositions
	if(output.eq.1)then
		write(12,*)'  '
		write(12,*)'  '
		write(12,*)'xxxxxxxxxxxxxx-Done-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
		write(12,*)' Summary for each phase pair '
		write(12,*)'     phase       Aff_Start   Aff_End     mmol_start  mmol_end    mmol_delta '

465		format(2(A8,4x))
		write(12,1564)PhName(k2),MDFAff(1,2),MDFAff(2,2),		&
			1000.*MDFmoles(1,2),1000.*MDFmoles(2,2),1000.*(MDFmoles(2,2)-MDFmoles(1,2))
		write(12,1564)PhName(k3),MDFAff(1,3),MDFAff(2,3),		&
			1000.*MDFmoles(1,3),1000.*MDFmoles(2,3),1000.*(MDFmoles(2,3)-MDFmoles(1,3))
		write(12,1564)PhName(k4),MDFAff(1,4),MDFAff(2,4),		&
			1000.*MDFmoles(1,4),1000.*MDFmoles(2,4),1000.*(MDFmoles(2,4)-MDFmoles(1,4))
1564		format(5x,A8,5F12.5)
		write(12,*)'     GBcomp        Start       End         Change '
		do j = 1,numPhCo(1)
			write(12,464)PhCoName(1,j),GBComp(1,j),GBComp(2,j),GBComp(2,j)-GBComp(1,j)
464			format(5x,A8,3F12.5)
			end do
		write(12,*)PhName(k2)
		do j = 1,numPhCo(k2)
			write(12,464)PhCoName(k2,j),XPhCo(k2,j)
			end do
		write(12,*)PhName(k3)
		do j = 1,numPhCo(k3)
			write(12,464)PhCoName(k3,j),XPhCo(k3,j)
			end do
		write(12,*)PhName(k4)
		do j = 1,numPhCo(k4)
			write(12,464)PhCoName(k4,j),XPhCo(k4,j)
			end do
		endif

99	continue
!	At this point, we need to 
		!(a) save the GB composition and 
		!(b) store the phase compositions to examine zoning
		!(c) store the amount of each phase produced or consumed 
	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		nodeComp(iNode,j) = xPhCo(1,j)
		end do
	do j = 1,numPhCo(k2)
		nodePhaseComp(iNode,1,j) = xPhCo(k2,j)			! nodePhaseComp(maxNodes,3,maxPhCo)
		end do
	do j = 1,numPhCo(k3)
		nodePhaseComp(iNode,2,j) = xPhCo(k3,j)
		end do
	do j = 1,numPhCo(k4)
		nodePhaseComp(iNode,3,j) = xPhCo(k4,j)
		end do

	nodePhaseMoles(iNode,1) = nodePhaseMoles(iNode,1) + MDFmoles(2,2)-MDFmoles(1,2)
	nodePhaseMoles(iNode,2) = nodePhaseMoles(iNode,2) + MDFmoles(2,3)-MDFmoles(1,3)
	nodePhaseMoles(iNode,3) = nodePhaseMoles(iNode,3) + MDFmoles(2,4)-MDFmoles(1,4)

	nodePhaseMolesDelta(iNode,1) = MDFmoles(2,2)-MDFmoles(1,2)
	nodePhaseMolesDelta(iNode,2) = MDFmoles(2,3)-MDFmoles(1,3)
	nodePhaseMolesDelta(iNode,3) = MDFmoles(2,4)-MDFmoles(1,4)

	nodePhaseAffinity(iNode,1) = MDFAff(2,2)
	nodePhaseAffinity(iNode,2) = MDFAff(2,3)
	nodePhaseAffinity(iNode,3) = MDFAff(2,4)
	

	if(output.eq.1)then
		write(12,*)' '
		write(12,*)'*********************************************'
		write(12,*)' Record to write to output file'
		write(12,102)iNode,nodeX(iNode),nodeY(iNode)
102		format(I6,2F12.3,I5,'         Node number,  X, Y')
		write(12,108)(nodeMIFID(iNode,i),i=1,3)
108		format(3I8,'  nodeMIFID')
		write(12,104)numNodeReactionPhases(iNode),(nodeReactionPhases(iNode,k),k=1,numNodeReactionPhases(iNode))
104		format(4I5,'  numNodeReactionPhases, phaseID')
		MIFID = 1		! this is the grain boundary
		GBmoles = 0.0d0
		write(12,103)MIFID,PhName(MIFID),GBmoles,GBmoles,GBmoles,						&
				numPhCo(MIFID),(phCoName(MIFID,j),nodeComp(iNode,j),j=1,numPhCo(MIFID))
		do kk = 1,numNodeReactionPhases(iNode)
			k = nodeReactionPhases(iNode,kk)
			MIFID = nodeMIFID(iNode,k)
			write(12,103)MIFID,PhName(MIFID),nodePhaseMolesDelta(iNode,k),nodePhaseMoles(iNode,k),		&
				nodePhaseAffinity(iNode,k),								&
				numPhCo(MIFID),(phCoName(MIFID,j),nodePhaseComp(iNode,k,j),j=1,numPhCo(MIFID))
103			format(T8,I4,2x,A8,3E15.5,I8,20(5x,A8,F12.5))
			end do
		write(12,*)'*********************************************'
		endif

	
	return
	end	

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE MDFNodesWith2Phases_3(iNode,numPhT,ph1,ph2,output,initializeOnly,igo)
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Diffuse_GB_Hex.inc"
	include "AutoAffinity.inc"
	include "GB_Tangent.inc"
! *****************************************
	integer*4 iNode,k1,k2,k3,k,i,j,icount,output,initializeOnly,numPhT,ph1,ph2,kk,kk1,kk2,MIFID,igo,izero
	real*8 deltaMolesCalc,deltaMolesdX,dF_dDelta,deltaAffinityInit,deltaAffinityOld,deltaAffinityNew,GBmoles
	real*8 deltaMolesX(4),GBChangeX(10),deltaAffinityTarget

!	if output = 1 then print
!	if output = 0 then don't print


3	continue		
	numPh = numPhT + 1		! should be 3

	asmCurrent(1) = 1			! always the grain boundary
	asmCurrent(2) = nodeMIFID(iNode,ph1)
	asmCurrent(3) = nodeMIFID(iNode,ph2)
	k1 = asmCurrent(1)
	k2 = asmCurrent(2)
	k3 = asmCurrent(3)
	if(output.eq.1)then
		write(12,*)'  '
		write(12,*)'  '
		write(12,*)'  '
		write(12,*)'********************************************************'
		write(12,*)' In Subroutine MDFNodesWith2Phases '
		write(12,*)' Node =  ',iNode		
		write(12,*)'********************************************************'
		write(12,*)' Current mineral assemblage'
		do k = 1,numPh
			WRITE(12,1533)MINREC(asmCurrent(K)),PHNAME(asmCurrent(K))
1533	  		FORMAT(3(I4,2X,A8,2x))
			end do
		endif

!	set the composition of the GB into the XPhCo array for use in Sub CalculateGBChange and ParallelToTangent
	do j = 1,numPhCo(1)		! grain boundary
		xPhCo(1,j) = nodeComp(iNode,j)		! numPhCo(1) should be every element
		end do
	do i = 1,2
		k = asmCurrent(i+1)
		kk = nodeReactionPhases(iNode,i)	!Note that the index of phases around the node isn't 1,2,3 necessarily
		select case(fractl(k))
		case(0,1)	! phases are not fractionating or are growing and fractionating
			do j = 1,numPhCo(k)
				xPhCo(k,j) = nodePhaseComp(iNode,kk,j)
				end do
		case(2)		! phase is being consumed - use fixed composition (i.e. chlorite)
			do j = 1,numPhCo(k)
				xPhCo(k,j) = xPhCoInitial(k,j)
				end do
		case default
		end select		

		end do


	CALL AdjustAsm(1)		! we changed assemblage. Adjust volumes, moles etc.


	! store the initial grain boundary composition
	! this will be used later to determine the change
	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		GBComp(1,j) = xPhCo(1,j)
		end do

	icount = 0


!	Initial calculation of affinity
	if(output.eq.1)then
		write(12,*)'-------------'
		write(12,*)' Initial call to CalcAffinity to get chemical potentials and affinities'
		call CalcAffinity(1)	! print this initial calculation
		write(12,*)' minrec  Ph            Aff       mmol    '
			   !    1  Quartz           0.9    0.43944E+02   32  Garnet        -386.7    0.85951E+01
		write(12,592)minrec(k2),PhName(k2),Aff(k2,1),1000.0d0*mp0(k2)
		write(12,592)minrec(k3),PhName(k3),Aff(k3,1),1000.0d0*mp0(k3)
592		Format(2(I5,2x,A8,F12.1,E15.5))	
		else
		call CalcAffinity(0)	! do not print this initial calculation
		endif

	call Save_and_Reset_Comps(1)	! save the current phase compositions

	do i = 1,2
		k = asmCurrent(i+1)
		call paralleltotangent(k,0,izero)
		if(izero.gt.0)then
			write(12,*)' In routine MDFNodesWith2Phases_3. Line 1540'
			write(12,*)' Node = ',iNode
			write(12,*)'    i = ',i
			write(12,*)'    k = ',k
			write(12,*)' asmCurrent(j) = ',(asmCurrent(j),j=2,3)
			endif
		end do
	if(output.eq.1)then
		write(12,*)'                 gDifferenceAU      Phase Comp '
		write(12,81)phName(k2),gDifferenceAU(k2),(xPhCo(k2,j),j=1,numPhCo(K2))
		write(12,81)phName(k3),gDifferenceAU(k3),(xPhCo(k3,j),j=1,numPhCo(K3))
	81	format(A10,T15,15F15.5)
		endif
	Aff(k2,1) = gDifferenceAU(k2)
	Aff(k3,1) = gDifferenceAU(k3)
	if(initializeOnly.eq.1)go to 99		! this call was only to establish starting compositions for all nodes

	MDFmoles(1,2) = mp0(k2)
	MDFmoles(1,3) = mp0(k3)
	MDFAff(1,2) = Aff(k2,1)
	MDFAff(1,3) = Aff(k3,1)
		
	deltaAffinityInit = Aff(k3,1) - Aff(k2,1)
	deltaAffinityOld = deltaAffinityInit		! This is the function we are trying to zero out
!	deltaMolesdX = 1.0d-4		! initial scaling for ∆moles for determining derivative of deltaAffinity
	deltaMolesdX = 1.0d-5		! initial scaling for ∆moles for determining derivative of deltaAffinity

		
! 100	continue		! loop here for Newton's method on reaction progress
	icount = icount + 1
	if(output.eq.1)then
		write(12,*)' ----------------------------------'
		write(12,*)' ----------------------------------'
		write(12,*)' icount = ',icount
		endif
	if(icount.gt.100)then
		write(*,*)' No convergence in MDFNodesWith2Phases. iNode = ',iNode
		write(75,*)' No convergence in MDFNodesWith2Phases. iNode = ',iNode
		write(*,*)' 0 = continue'
		write(*,*)' 1 = abort'
		read(*,*)igo
		if(igo.eq.1)return
		go to 599	! exit routine
		endif

	
!	Calculate the derivative
	deltaMolesX(2) = deltaMolesdX
	deltaMolesX(3) = -deltaMolesdX	! this must be negative so that sum oxygen = 0


	! now calculate the change in GB composition based on the above
	call Save_and_Reset_Comps(2)	! reset the current phase compositions after call to ParallelToTangent
	Call CalculateGBChange(deltaMolesX,GBChangeX)
	! Note that the array XPhCo(1,j) is now changed

	if(output.eq.1)then
		write(12,*)' '
		write(12,*)' deltaMolesX ',deltaMolesX(2),deltaMolesX(3)
		write(12,*)' GB CHANGE   ',(GBChangeX(j),j=1,numPhCo(1))	
		endif
	! Now we have the composition of the grain boundary for the small increment of moles (deltaMolesdX)
	! Calculate the affinities

	do i = 1,2
		k = asmCurrent(i+1)
		call paralleltotangent(k,0,izero)	! 1 = output summary info only; 2 = output everything
		if(izero.gt.0)then
			write(12,*)' In routine MDFNodesWith2Phases_3. Line 1609'
			write(12,*)' Node = ',iNode
			write(12,*)'    i = ',i
			write(12,*)'    k = ',k
			write(12,*)' asmCurrent(j) = ',(asmCurrent(j),j=2,3)
			endif
		end do
	if(output.eq.1)then
		write(12,*)'                 gDifferenceAU      Phase Comp '
		write(12,81)phName(k2),gDifferenceAU(k2),(xPhCo(k2,j),j=1,numPhCo(K2))
		write(12,81)phName(k3),gDifferenceAU(k3),(xPhCo(k3,j),j=1,numPhCo(K3))
		endif
	Aff(k2,1) = gDifferenceAU(k2)
	Aff(k3,1) = gDifferenceAU(k3)

! 	call CalcAffinity(0)

	deltaAffinityNew = (Aff(k3,1) - Aff(k2,1))	!This is the change in ∆Aff with a small change in moles -- to get slope
! 	if(abs(deltaAffinityNew).lt.1.0d0)then		! 1 Joule convergence
	if(abs(deltaAffinityNew).lt.convergenceJoules)then		! Convergence criteria set in model input file
		! we have converged
		if(output.eq.1)write(12,*)' Jumping to exit '
		go to 599
		endif

	! calculate the derivative
	dF_dDelta = (deltaAffinityNew - deltaAffinityOld)/deltaMolesdX
	if(output.eq.1)then
		write(12,*)' delAffNew     deAffOld     deltamolesdX ',deltaAffinityNew,deltaAffinityOld,deltaMolesdX
		write(12,*)' Slope = ',dF_dDelta
		endif

!	Solve for the desired new deltaAffinity that we want based on the ∆time incerment
	deltaAffinityTarget = deltaAffinityOld*exp(-timestep)
	! solve for the change in moles to find the root
	deltaMolesCalc = (deltaAffinityTarget-deltaAffinityOld)/dF_dDelta	! Newton's method


	! solve for the change in moles to find the root
! 	deltaMolesCalc = -deltaAffinityOld/dF_dDelta	! Newton's method
	! deltaMoles is the change in moles of phase k2 to get the derivative to zero.
	! Now use this deltaMoles to get a new value for the moles 
	! Calculate the new values of the moles of the 2 phases
! 	deltaMolesCalc = deltaMolesCalc*newtonDamp
	deltaMolesX(2) =  deltaMolesCalc
	deltaMolesX(3) = -deltaMolesCalc	! this must be negative so that sum oxygen = 0

	mp0(k2) = mp0(k2) + deltaMolesX(2)
	mp0(k3) = mp0(k3) + deltaMolesX(3)
	
	! now calculate the change in GB composition based on the above root
	call Save_and_Reset_Comps(2)	! reset the current phase compositions after call to ParallelToTangent
	Call CalculateGBChange(deltaMolesX,GBChangeX)

	if(output.eq.1)then
		write(12,*)' '
		write(12,*)' deltaMolesX ',deltaMolesX(2),deltaMolesX(3)
		write(12,*)' GB CHANGE   ',(GBChangeX(j),j=1,numPhCo(1))	
		endif
	! Now we have the composition of the grain boundary for the small increment of moles (deltaMolesdX)
	! Calculate the affinities
	do i = 1,2
		k = asmCurrent(i+1)
		call paralleltotangent(k,0,izero)
		if(izero.gt.0)then
			write(12,*)' In routine MDFNodesWith2Phases_3. Line 1665'
			write(12,*)' Node = ',iNode
			write(12,*)'    i = ',i
			write(12,*)'    k = ',k
			write(12,*)' asmCurrent(j) = ',(asmCurrent(j),j=2,3)
			endif
		end do
	if(output.eq.1)then
		write(12,*)'                 gDifferenceAU      Phase Comp '
		write(12,81)phName(k2),gDifferenceAU(k2),(xPhCo(k2,j),j=1,numPhCo(K2))
		write(12,81)phName(k3),gDifferenceAU(k3),(xPhCo(k3,j),j=1,numPhCo(K3))
		endif
	Aff(k2,1) = gDifferenceAU(k2)
	Aff(k3,1) = gDifferenceAU(k3)
! 	call CalcAffinity(0)
	deltaAffinityNew = (Aff(k3,1) - Aff(k2,1))	!This is the change in ∆Aff with the calculated ∆moles (deltaMolesCalc) 

	if(output.eq.1)then
		write(12,*)' DeltaAffNew       ',deltaAffinityNew	
		write(12,*)'                   Ph                               Delta(mmol)          mmol                Aff'
		Write(12,*)' Deltax(2) (mmoles)',PhName(k2),1000.0d0*deltaMolesX(2),1000.0d0*mp0(k2),Aff(k2,1)
		write(12,*)' Deltax(3) (mmoles)',PhName(k3),1000.0d0*deltaMolesX(3),1000.0d0*mp0(k3),Aff(k3,1)
		write(12,*)' '
		Write(12,*)'---GB composition changes --'
		write(12,501)(PhCoName(1,j),j=1,numPhCo(1))
501		format(10A15)
		write(12,502)(GBChangeX(j),j=1,numPhCo(1))
		write(12,502)(xPhCo(1,j),j=1,numPhCo(1))
502		format(10E15.4)

		endif


			
!	if we get here, then we need to stop growing
599	continue

	! we will need to put code here to account for the change in phase composition if numMDFRxn > 0

!	Store the final moles
	MDFmoles(2,2) = mp0(k2)
	MDFmoles(2,3) = mp0(k3)
	MDFAff(2,2) = Aff(k2,1)
	MDFAff(2,3) = Aff(k3,1)

!	store the final grain boundary composition
	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		GBComp(2,j) = xPhCo(1,j)
		end do

!500	continue	! loop on every pair		
	
	if(output.eq.1)then
		!write out initial and final GB compositions
		write(12,*)'  '
		write(12,*)'  '
		write(12,*)'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
		write(12,*)' Summary for this phase pair '
		write(12,565)PhName(k2),PhName(k3)
565		format(2(A8,4x))
		write(12,*)'     phase       Aff_Start   Aff_End     mmol_start  mmol_end    mmol_delta '
		write(12,564)PhName(k2),MDFAff(1,2),MDFAff(2,2),		&
			1000.*MDFmoles(1,2),1000.*MDFmoles(2,2),1000.*(MDFmoles(2,2)-MDFmoles(1,2))
		write(12,564)PhName(k3),MDFAff(1,3),MDFAff(2,3),		&
			1000.*MDFmoles(1,3),1000.*MDFmoles(2,3),1000.*(MDFmoles(2,3)-MDFmoles(1,3))
564		format(5x,A8,5F12.5)
		write(12,*)'     GBcomp        Start       End         Change '
		do j = 1,numPhCo(1)
			write(12,564)PhCoName(1,j),GBComp(1,j),GBComp(2,j),GBComp(2,j)-GBComp(1,j)
			end do
		write(12,*)PhName(k2)
		do j = 1,numPhCo(k2)
			write(12,464)PhCoName(k2,j),XPhCo(k2,j)
			end do
		write(12,*)PhName(k3)
		do j = 1,numPhCo(k3)
			write(12,464)PhCoName(k3,j),XPhCo(k3,j)
			end do
464		format(5x,A8,3F12.5)
		endif

99	continue
!	At this point, we need to 
		!(a) save the GB composition and 
		!(b) store the phase compositions to examine zoning
		!(c) store the amount of each phase produced or consumed (I don't have an array for this yet)

	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		nodeComp(iNode,j) = xPhCo(1,j)
		end do
	kk1 = nodeReactionPhases(iNode,1)	!Note that the index of 2 unique phases around the node could be 1&2, 1&3 or 2&3
	kk2 = nodeReactionPhases(iNode,2)
	do j = 1,numPhCo(k2)
		nodePhaseComp(iNode,kk1,j) = xPhCo(k2,j)
		end do
	do j = 1,numPhCo(k3)
		nodePhaseComp(iNode,kk2,j) = xPhCo(k3,j)
		end do

	nodePhaseMoles(iNode,kk1) = nodePhaseMoles(iNode,kk1) + MDFmoles(2,2)-MDFmoles(1,2)
	nodePhaseMoles(iNode,kk2) = nodePhaseMoles(iNode,kk2) + MDFmoles(2,3)-MDFmoles(1,3)

	nodePhaseMolesDelta(iNode,kk1) = MDFmoles(2,2)-MDFmoles(1,2)
	nodePhaseMolesDelta(iNode,kk2) = MDFmoles(2,3)-MDFmoles(1,3)

	nodePhaseAffinity(iNode,kk1) = MDFAff(2,2)
	nodePhaseAffinity(iNode,kk2) = MDFAff(2,3)

	if(output.eq.1)then
		write(12,*)' '
		write(12,*)'*********************************************'
		write(12,*)' Record to write to output file'
		write(12,102)iNode,nodeX(iNode),nodeY(iNode)
102		format(I6,2F12.3,I5,'         Node number,  X, Y')
		write(12,108)(nodeMIFID(iNode,i),i=1,3)
108		format(3I8,'  nodeMIFID')
		write(12,104)numNodeReactionPhases(iNode),(nodeReactionPhases(iNode,k),k=1,numNodeReactionPhases(iNode))
104		format(I5,'  numNodeReactionPhases, phaseID')
		MIFID = 1		! this is the grain boundary
		GBmoles = 0.0d0
		write(12,103)MIFID,PhName(MIFID),GBmoles,GBmoles,GBmoles,						&
				numPhCo(MIFID),(phCoName(MIFID,j),nodeComp(iNode,j),j=1,numPhCo(MIFID))
		do kk = 1,numNodeReactionPhases(iNode)
			k = nodeReactionPhases(iNode,kk)
			MIFID = nodeMIFID(iNode,k)
			write(12,103)MIFID,PhName(MIFID),nodePhaseMolesDelta(iNode,k),nodePhaseMoles(iNode,k),		&
				nodePhaseAffinity(iNode,k),								&
				numPhCo(MIFID),(phCoName(MIFID,j),nodePhaseComp(iNode,k,j),j=1,numPhCo(MIFID))
103			format(T8,I4,2x,A8,3E15.5,I8,20(5x,A8,F12.5))
			end do
		write(12,*)'*********************************************'
		endif
	
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE MDFPairs_3(iSegTemp,iPointX,output,initializeOnly,igo)
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Diffuse_GB_Hex.inc"
	include "AutoAffinity.inc"
	include "GB_Tangent.inc"
! *****************************************
	integer*4 iSeg,iPointX,k1,k2,k3,k,i,j,icount,iSegTemp,output,initializeOnly,MIFID,igo,izero
	real*8 deltaMolesCalc,deltaMolesdX,dF_dDelta,deltaAffinityInit,deltaAffinityOld,deltaAffinityNew,GBmoles
	real*8 deltaMolesX(4),GBChangeX(10),deltaAffinityTarget

!	if output = 1 then print
!	if output = 0 then don't print

	iSeg = iSegTemp
3	continue		
	numPh = 3

	asmCurrent(1) = 1			! always the grain boundary
	asmCurrent(2) = segMIFID(iSeg,1)
	asmCurrent(3) = segMIFID(iSeg,2)
	k1 = asmCurrent(1)
	k2 = asmCurrent(2)
	k3 = asmCurrent(3)
	if(output.eq.1)then
		write(12,*)'  '
		write(12,*)'  '
		write(12,*)'********************************************************'
		write(12,*)' In Subroutine MDFPairs '
		write(12,*)' Segment = ',iSeg
		write(12,*)' iPointX = ',iPointX
		write(12,*)'********************************************************'
		write(12,*)' Current mineral assemblage'
		WRITE(12,433)(MINREC(asmCurrent(K)),PHNAME(asmCurrent(K)),k = 1,numPh)
433   		FORMAT(3(I4,2X,A8,2x))
		write(12,*)' '
		write(12,*)'*********************************************'
		write(12,*)' Initial record from input file  '
		MIFID = 1		! this is the grain boundary
		GBmoles = 0.0d0
		write(12,103)MIFID,PhName(MIFID),GBmoles,GBmoles,		&
				numPhCo(MIFID),(phCoName(MIFID,j),pointComp(iSeg,iPointX,j),j=1,numPhCo(MIFID))
		do k = 1,numSegReactionPhases(iSeg)
			MIFID = segMIFID(iSeg,k)
			write(12,103)MIFID,PhName(MIFID),pointPhaseMolesDelta(iSeg,iPointX,k),		&
				pointPhaseMoles(iSeg,iPointX,k),						&
				numPhCo(MIFID),(phCoName(MIFID,j),						&
				pointPhaseComp(iSeg,iPointX,k,j),j=1,numPhCo(MIFID))
			end do
		write(12,*)'*********************************************'
		endif

!	Set the composition of the GB into the XPhCo array for use in Sub CalculateGBChange and ParallelToTangent
	do j = 1,numPhCo(1)		! grain boundary
		xPhCo(1,j) = pointComp(iSeg,iPointX,j)		! numPhCo(1) should be every element
		end do
	do i = 1,2
		k = asmCurrent(i+1)
		!kk = segReactionPhases(iSeg,i)	!Note that the index of phases around the segment isn't 1,2 necessarily
		select case(fractl(k))
		case(0,1)	! phases are not fractionating or are growing and fractionating
			do j = 1,numPhCo(k)
				xPhCo(k,j) = pointPhaseComp(iSeg,iPointX,i,j)
				end do
		case(2)		! phase is being consumed - use fixed composition (i.e. chlorite)
			do j = 1,numPhCo(k)
				xPhCo(k,j) = xPhCoInitial(k,j)
				end do
		case default
		end select		
		end do


	CALL AdjustAsm(1)		! we changed assemblage. Adjust volumes, moles etc.


!		store the initial grain boundary composition
	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		GBComp(1,j) = xPhCo(1,j)
		end do

!	Initial calculation of affinity
	if(output.eq.1)then
		write(12,*)'-------------'
		write(12,*)' Initial call to CalcAffinity to get chemical potentials and affinities'
		call CalcAffinity(1)	! print this initial calculation
		write(12,*)' minrec  Ph            Aff       mmol   '
		write(12,592)minrec(k2),PhName(k2),Aff(k2,1),1000.0d0*mp0(k2)
		write(12,592)minrec(k3),PhName(k3),Aff(k3,1),1000.0d0*mp0(k3)
592		Format(2(I5,2x,A8,F12.1,E15.5))	
		else
		call CalcAffinity(0)	! do not print this initial calculation
		endif

	call Save_and_Reset_Comps(1)	! save the current phase compositions

	do i = 1,2
		k = asmCurrent(i+1)
		call paralleltotangent(k,0,izero)
		if(izero.gt.0)then
			write(12,*)' In routine MDFPairs_3. Line 1929'
			write(12,*)' iSeg = ',iSeg
			write(12,*)'    i = ',i
			write(12,*)'    k = ',k
			write(12,*)' asmCurrent(j) = ',(asmCurrent(j),j=2,3)
			endif

		end do
	if(output.eq.1)then
		write(12,*)'                 gDifferenceAU      Phase Comp '
		write(12,81)phName(k2),gDifferenceAU(k2),(xPhCo(k2,j),j=1,numPhCo(K2))
		write(12,81)phName(k3),gDifferenceAU(k3),(xPhCo(k3,j),j=1,numPhCo(K3))
		endif
	Aff(k2,1) = gDifferenceAU(k2)
	Aff(k3,1) = gDifferenceAU(k3)

	if(initializeOnly.eq.1)go to 99		! this call was only to establish starting compositions for all nodes

	MDFmoles(1,2) = mp0(k2)
	MDFmoles(1,3) = mp0(k3)
 	MDFAff(1,2) = Aff(k2,1)
 	MDFAff(1,3) = Aff(k3,1)

!	determine which phase has greater affinity- that is the phase that will grow
	deltaAffinityInit = Aff(k3,1) - Aff(k2,1)
	deltaAffinityOld = deltaAffinityInit	
! 	deltaMolesdX = 1.0d-6		! initial scaling for ∆moles for determining derivative of deltaAffinity
 	deltaMolesdX = 1.0d-5		! initial scaling for ∆moles for determining derivative of deltaAffinity
! 	deltaMolesdX = 1.0d-3		! initial scaling for ∆moles for determining derivative of deltaAffinity
!  	deltaMolesdX = 1.0d-1		! initial scaling for ∆moles for determining derivative of deltaAffinity
! 	deltaMolesdX = 1.0d0		! initial scaling for ∆moles for determining derivative of deltaAffinity

!  	newtonDamp =0.5d0
! 	newtonDamp =1.0d0
	icount = 0
100	continue
! 550	continue		! loop until affinities are equal
	icount = icount + 1
	if(output.eq.1)then
		write(12,*)' ----------------------------------'
		write(12,*)' ----------------------------------'
		write(12,*)' icount = ',icount
		endif
	if(icount.gt.100)then
		write(*,*)' No convergence in MDFPairs. iSeg,iPoint = ',iSeg,iPointX
		write(75,*)' No convergence in MDFPairs. iSeg,iPoint = ',iSeg,iPointX
		write(12,*)' ************* --- No convergence in MDFPairs. iSeg,iPoint = ',iSeg,iPointX
		write(*,*)' 0 = continue'
		write(*,*)' 1 = abort'
		read(*,*)igo
		if(igo.eq.1)return
		go to 599	! exit routine
		endif

!	Calculate the derivative

!	Set initial deltas to deltaMoles (scaled to number of oxygens) so we can determine dF/delta
	deltaMolesX(2) = deltaMolesdX
	deltaMolesX(3) = -deltaMolesdX	! this must be negative so that sum oxygen = 0

	! now calculate the change in GB composition based on the above
	call Save_and_Reset_Comps(2)	! reset the current phase compositions after call to ParallelToTangent
	Call CalculateGBChange(deltaMolesX,GBChangeX)
	! Note that we have now changed the XPhCo array for the grain boundary
	if(output.eq.1)then
		write(12,*)' '
		write(12,*)' deltaMolesX ',deltaMolesX(2),deltaMolesX(3)
		write(12,*)' GB CHANGE   ',(GBChangeX(j),j=1,numPhCo(1))	
		endif

	! Now we have the composition of the grain boundary for the small increment of moles (deltaMolesdX)
	! Calculate the affinities
	do i = 1,2
		k = asmCurrent(i+1)
		call paralleltotangent(k,0,izero)
		if(izero.gt.0)then
			write(12,*)' In routine MDFPairs_3. Line 2003'
			write(12,*)' iSeg = ',iSeg
			write(12,*)'    i = ',i
			write(12,*)'    k = ',k
			write(12,*)' asmCurrent(j) = ',(asmCurrent(j),j=2,3)
			endif
		end do
	if(output.eq.1)then
		write(12,*)'                 gDifferenceAU      Phase Comp '
		write(12,81)phName(k2),gDifferenceAU(k2),(xPhCo(k2,j),j=1,numPhCo(K2))
		write(12,81)phName(k3),gDifferenceAU(k3),(xPhCo(k3,j),j=1,numPhCo(K3))
		endif		
! 	call CalcAffinity(0)

	Aff(k2,1) = gDifferenceAU(k2)
	Aff(k3,1) = gDifferenceAU(k3)


	deltaAffinityNew = (Aff(k3,1) - Aff(k2,1))	!This is the change in ∆Aff with a small change in moles -- to get slope
	! check to be sure we need a new calculation. Otherwise, if slope--> 0 we get a NAN
! 	if(abs(deltaAffinityNew).lt.1.0d0)then		! 1 Joule convergence
 	if(abs(deltaAffinityNew).lt.convergenceJoules)then	! convergence in Joules is set in model input file

		! we have converged
		if(output.eq.1)write(12,*)' Jumping to exit '
		go to 599
		endif

	! calculate the derivative
	dF_dDelta = (deltaAffinityNew - deltaAffinityOld)/deltaMolesdX		! this is the slope
	if(output.eq.1)then
		write(12,*)' delAffNew     deAffOld     deltamolesdX ',deltaAffinityNew,deltaAffinityOld,deltaMolesdX
		write(12,*)' Slope = ',dF_dDelta
		endif

!	Solve for the desired new deltaAffinity that we want based on the ∆time incerment
	deltaAffinityTarget = deltaAffinityOld*exp(-timestep)
	! solve for the change in moles to find the root
	deltaMolesCalc = (deltaAffinityTarget-deltaAffinityOld)/dF_dDelta	! Newton's method


	! solve for the change in moles to find the root
! 	deltaMolesCalc = -deltaAffinityOld/dF_dDelta	! Newton's method
	! deltaMoles is the change in moles of phase k2 to get the derivative to zero.
	! Now use this deltaMoles to get a new value for the moles 
	! Calculate the new values of the moles of the 2 phases
! 	deltaMolesCalc = deltaMolesCalc*newtonDamp
	deltaMolesX(2) =  deltaMolesCalc
	deltaMolesX(3) = -deltaMolesCalc	! this must be negative so that sum oxygen = 0

	mp0(k2) = mp0(k2) + deltaMolesX(2)
	mp0(k3) = mp0(k3) + deltaMolesX(3)
	
	! now calculate the change in GB composition based on the above root
	call Save_and_Reset_Comps(2)	! reset the current phase compositions after call to ParallelToTangent
	Call CalculateGBChange(deltaMolesX,GBChangeX)
	! Note that we have now changed the XPhCo array for the grain boundary
	if(output.eq.1)then
		write(12,*)' '
		write(12,*)' deltaMolesX ',deltaMolesX(2),deltaMolesX(3)
		write(12,*)' GB CHANGE   ',(GBChangeX(j),j=1,numPhCo(1))	
		endif

	! Now we have the composition of the grain boundary for the small increment of moles (deltaMolesdX)
	! Calculate the affinities
	do i = 1,2
		k = asmCurrent(i+1)
		call paralleltotangent(k,0,izero)
		if(izero.gt.0)then
			write(12,*)' In routine MDFPairs_3. Line 2063'
			write(12,*)' iSeg = ',iSeg
			write(12,*)'    i = ',i
			write(12,*)'    k = ',k
			write(12,*)' asmCurrent(j) = ',(asmCurrent(j),j=2,3)
			endif
		end do
	if(output.eq.1)then
		write(12,*)'                 gDifferenceAU      Phase Comp '
		write(12,81)phName(k2),gDifferenceAU(k2),(xPhCo(k2,j),j=1,numPhCo(K2))
		write(12,81)phName(k3),gDifferenceAU(k3),(xPhCo(k3,j),j=1,numPhCo(K3))
		endif
81	format(A10,T15,15F15.5)

	Aff(k2,1) = gDifferenceAU(k2)
	Aff(k3,1) = gDifferenceAU(k3)
 	deltaAffinityNew = (Aff(k3,1) - Aff(k2,1))	!This is the change in ∆Aff with the calculated ∆moles (deltaMolesCalc) 

598	continue
	if(output.eq.1)then
		write(12,*)' '
		write(12,*)' DeltaAffNew       ',deltaAffinityNew	
		write(12,*)'                   Ph                               Delta(mmol)          mmol                Aff'
		Write(12,*)' Deltax(2) (mmoles)',PhName(k2),1000.0d0*deltaMolesX(2),1000.0d0*mp0(k2),Aff(k2,1)
		write(12,*)' Deltax(3) (mmoles)',PhName(k3),1000.0d0*deltaMolesX(3),1000.0d0*mp0(k3),Aff(k3,1)
		write(12,*)' '
		Write(12,*)'---GB composition changes --'
		write(12,501)(PhCoName(1,j),j=1,numPhCo(1))
501		format(10A15)
		write(12,502)(GBChangeX(j),j=1,numPhCo(1))
		write(12,502)(xPhCo(1,j),j=1,numPhCo(1))
502		format(10E15.4)

		endif


			
!	if we get here, then we need to stop growing
599	continue

!	Store the final moles
	MDFmoles(2,2) = mp0(k2)
	MDFmoles(2,3) = mp0(k3)
	MDFAff(2,2) = Aff(k2,1)
	MDFAff(2,3) = Aff(k3,1)

!	store the final grain boundary composition
	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		GBComp(2,j) = xPhCo(1,j)
		end do

500	continue	! loop on every pair		
	
!		write out initial and final GB compositions
	if(output.eq.1)then
		write(12,492)minrec(k2),PhName(k2),Aff(k2,1),1000.*mp0(k2),minrec(k3),PhName(k3),Aff(k3,1),1000.*mp0(k3)
492		Format(2(I5,2x,A8,F12.1,E15.5))	
		write(12,*)'  '
		write(12,*)'  '
		write(12,*)'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
		write(12,*)' Summary for this phase pair '
		write(12,565)PhName(k2),PhName(k3)
565		format(2(A8,4x))
		write(12,*)'     phase       Aff_Start   Aff_End     mmol_start  mmol_end    mmol_delta '
		write(12,564)PhName(k2),MDFAff(1,2),MDFAff(2,2),		&
			1000.*MDFmoles(1,2),1000.*MDFmoles(2,2),1000.*(MDFmoles(2,2)-MDFmoles(1,2))
		write(12,564)PhName(k3),MDFAff(1,3),MDFAff(2,3),		&
			1000.*MDFmoles(1,3),1000.*MDFmoles(2,3),1000.*(MDFmoles(2,3)-MDFmoles(1,3))
564		format(5x,A8,5F12.5)
		write(12,*)'     GBcomp        Start       End         Change '
		do j = 1,numPhCo(1)
			write(12,564)PhCoName(1,j),GBComp(1,j),GBComp(2,j),GBComp(2,j)-GBComp(1,j)
			end do
		write(12,*)PhName(k2)
		do j = 1,numPhCo(k2)
			write(12,464)PhCoName(k2,j),XPhCo(k2,j)
			pointPhaseComp(iSeg,iPointX,1,j) = xPhCo(k2,j)
			end do
		write(12,*)PhName(k3)
		do j = 1,numPhCo(k3)
			write(12,464)PhCoName(k3,j),XPhCo(k3,j)
			pointPhaseComp(iSeg,iPointX,2,j) = xPhCo(k3,j)
			end do
464		format(5x,A8,3F12.5)
		endif


99	continue
!	At this point, we need to 
		!(a) save the GB composition and 
		!(b) store the phase compositions to examine zoning
		!(c) store the amount of each phase produced or consumed (I don't have an array for this yet)

	do j = 1,numPhCo(1)		! grain boundary is always phase 1
		pointComp(iSeg,iPointX,j) = xPhCo(1,j)
		end do
	do j = 1,numPhCo(k2)
		pointPhaseComp(iSeg,iPointX,1,j) = xPhCo(k2,j)
		end do
	do j = 1,numPhCo(k3)
		pointPhaseComp(iSeg,iPointX,2,j) = xPhCo(k3,j)
		end do

	pointPhaseMoles(iSeg,iPointX,1) = pointPhaseMoles(iSeg,iPointX,1) + MDFmoles(2,2)-MDFmoles(1,2)
	pointPhaseMoles(iSeg,iPointX,2) = pointPhaseMoles(iSeg,iPointX,2) + MDFmoles(2,3)-MDFmoles(1,3)	

	pointPhaseMolesDelta(iSeg,iPointX,1) = MDFmoles(2,2)-MDFmoles(1,2)
	pointPhaseMolesDelta(iSeg,iPointX,2) = MDFmoles(2,3)-MDFmoles(1,3)	
	
	! Do not save or write out affinities for segments -- only nodes.

	if(output.eq.1)then
		write(12,*) MDFmoles(2,2)-MDFmoles(1,2),1000.*(MDFmoles(2,2)-MDFmoles(1,2)),pointPhaseMolesDelta(iSeg,iPointX,1)
		write(12,*) MDFmoles(2,3)-MDFmoles(1,3),1000.*(MDFmoles(2,3)-MDFmoles(1,3)),pointPhaseMolesDelta(iSeg,iPointX,2)
		write(12,*)' '
		write(12,*)'*********************************************'
		write(12,*)' Record to write to output file'
		MIFID = 1		! this is the grain boundary
		GBmoles = 0.0d0
		write(12,103)MIFID,PhName(MIFID),GBmoles,GBmoles,GBmoles,						&
				numPhCo(MIFID),(phCoName(MIFID,j),pointComp(iSeg,iPointX,j),j=1,numPhCo(MIFID))
		do k = 1,numSegReactionPhases(iSeg)
			MIFID = segMIFID(iSeg,k)
			write(12,103)MIFID,PhName(MIFID),pointPhaseMolesDelta(iSeg,iPointX,k),		&
				pointPhaseMoles(iSeg,iPointX,k),						&
				numPhCo(MIFID),(phCoName(MIFID,j),						&
				pointPhaseComp(iSeg,iPointX,k,j),j=1,numPhCo(MIFID))
			end do
103			format(T8,I4,2x,A8,2E15.5,I8,20(5x,A8,F12.5))
		write(12,*)'*********************************************'
		endif
		
	return
	end

	
