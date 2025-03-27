! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PolyModelRoutine()
	USE AWE_INTERFACES
	implicit none
!	TYPE(AWE_Canvas) :: SegPlot(3)
 	TYPE(AWE_Canvas) :: GridPlot
! 	TYPE(AWE_Canvas) :: SegPlot1
! 	TYPE(AWE_Canvas) :: SegPlot2
! 	TYPE(AWE_Canvas) :: SegPlot3
	TYPE(AWE_Canvas) :: MyWindow		! The plotting canvas

! *****************************************
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"
	include "AutoAffinity.inc"
	include "Gibbsfiles.inc"
! *****************************************
! 	integer*4 iTest,iNode,iSeg,iPointX,i,j,jj,numIter,numCycles,iunit,isteps,numSteps,iok,startStep,status,numCyclesInput
! 	integer*4 iNodeStart,iNodeEnd,ph1,ph2,ph3,debug,endstep,lastCycle,calcSegments,igo,L,M
! 	integer*4 k,iXl,NodeCaptureOn,MoveNodesandSegs,error,abort,segStart,segEnd,output,garnetMIFID
! 	character*16 date,time,zone,steptext
! 	integer timevalues(8)
! 	real hoursInit,minutesInit,secondsInit
! 	real*8 moleScaleFactor,lengthScaleFactor
! 	character*32 TPtext

	Integer*4 status,numPolySteps,PolyStartStep,iSteps,iUnit,GarnetMIFID,numSteps,startStep,numCyclesInput
	integer*4 calcsegments,movenodesandsegs,nodecaptureon,debug,endstep,timevalues(8),numCycles
	integer*4 i,j,jj,iPoly,iNode,ph1,ph2,ph3,iSeg,iPointX,iNodeStart,iNodeEnd,error,abort,igo,iok
 	real hoursInit,minutesInit,secondsInit
 	real*8 moleScaleFactor
 	character*16 date,time,zone,steptext,dummy
 	character*32 TPtext

	! this routine is designed to model oscillatory zoning in garnet
	! The approach is to change the number of diffusion iterations (D/R) in the middle of the run
	! the idea being that Mn will decrease with low D/R and increase if the D/R is increases

	computeRoutine = 2
	totalCycles = 0
	timeStep = 0.1d0	! this is to calculate ∆Aff(target) = ∆Aff(old)*exp(-timeStep)

	write(*,*)'Open the file containing the PolyModel information'
	call FSS_Alert('Alert','Open PolyModel info file')
	open(41,file='',status = 'OLD',iostat = status)
	if(status.ne.0)then
		call FSS_Alert('Alert','Problem opening PolyModel info file')
		return
		endif
	INQUIRE(41, NAME=PolyModelInfoFile)
	write(*,*)'PolyModelInfoFile name = ',PolyModelInfoFile
	
	read(41,*)dummy		! just a title line
	read(41,*)numPolySteps
	read(41,*)ModelOutputFileBase
	read(41,*)PolyStartStep
	if(PolyStartStep.ne.1)then
		call FSS_Alert('ALERT','Be sure you have opened model file and last output file before continuing')
		write(*,*)'OK to continue? 0 = NO -- ; 1 = OK'
		read(*,*)iOK
		if(iOK.eq.0)go to 10
		endif

	! Open and close base file for output
	if(PolyStartStep.eq.1)then
		open(73,FILE=ModelOutputFileBase,status='NEW',iostat=status,action='WRITE')
		if(status.ne.0)then
			call FSS_Alert('ALERT','Problem opening base file. Be sure the name is correct and doesnot already exist')
			return
			endif
		inquire(73,NAME=ModelOutputFileBase)
		write(*,*)' ModelOutputFileBase =  ',ModelOutputFileBase
		write(73,*)' This file is intentionally left empty'
		write(73,*)' Do not delete -- it is used as the base file name for creating animations'
		close(73)
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
	if(Polystartstep.eq.1)then
		write(steptext,55)Polystartstep-1
		else
		write(steptext,55)Polystartstep
		endif			
	steptext = adjustL(steptext)
	ModelLogFile = trim(ModelOutputFileBase)//'__'//trim(steptext)//'.log'
	open(75,FILE=ModelLogFile,status = 'UNKNOWN')

	ModelStuffFile = trim(ModelOutputFileBase)//'__'//trim(steptext)//'.stuff.txt'
	open(76,FILE=ModelStuffFile,status = 'UNKNOWN')
	
	if(PolyStartStep.eq.1)then	! write the initial grid step if this is the very first
		isteps = 0
		write(steptext,55)isteps
55		format(I5)
		steptext = adjustL(steptext)
		ModelOutputFile = trim(ModelOutputFileBase)//'_'//trim(steptext)//'.GBM'
		open(74,FILE=ModelOutputFile,status = 'UNKNOWN')
		iunit = 74
		call WriteAllHeader(iunit)
		call WriteAllNodesAndSegments(iunit)
		close(iunit)
		endif

	do i = 1,numPhMIF
		if(trim(phName(i)).eq.'Garnet')then
			garnetMIFID = i
			write(*,*)' GarnetMIFID = ',garnetMIFID
			endif
		end do


	do iPoly = 1,numPolySteps

		read(41,*)dummy		! row of dashes
		read(41,*)numSteps	! number of model steps 
		read(41,*)startStep	! the step number for the first model. This should be one greater than
					!    last step number
! 		write(*,*)' Input number of model steps (each step will take several minutes)'
! 		write(*,*)' Note that  -- save and write out every step'
! 		read(*,*)numSteps
! 		if(numSteps.eq.0)go to 10
! 		write(*,*)' Input the starting step. = 1 if this is a new model, or = ## if it is a continuation'
! 		read(*,*)startStep
! 		if(startStep.ne.1)then
! 			call FSS_Alert('ALERT','Be sure you have opened model file and last output file before continuing')
! 			write(*,*)'OK to continue? 0 = NO -- ; 1 = OK'
! 			read(*,*)iOK
! 			if(iOK.eq.0)go to 10
! 			endif
		
	
		read(41,*)numCyclesInput
!		write(*,*)'Input number of MDF + Diffusion cycles to run for each step'
!		write(*,*)'Note that each cycle might take around 1 minute, so plan accordingly'
!		read(*,*)numCyclesInput
!		if(numCyclesInput.eq.0)go to 10
		read(41,*)numDiffIterations
! 		write(*,*)' Input number of diffusion iterations (e.g. 50-5000)'
! 		write(*,*)' If you specify a number =100,000 then the GB composition will be averaged after reaction step '
! 		write(*,*)'      (e.g. infinite diffusion) '
! 		read(*,*)numDiffIterations
	
!		This code ensures the first 10 steps are 1 cycle each
!		I don't think we need this with PolyModel
! 		lastCycle = 0
! 		if(startStep.eq.1)then
! 			lastCycle = lastCycle + 10
! 			endif
! 		if(numSteps.gt.10.and.startStep.eq.1)then
! 			lastCycle = lastCycle + numCyclesInput*(numSteps - 10)
! 			endif
! 		if(numSteps.gt.10.and.startStep.ne.1)then
! 			lastCycle = totalCycles + numCyclesInput*(numSteps)	! Total cycles is from the last run or read from the output file
! 			endif
		!lastCycle = numSteps*numCyclesInput
	
		read(41,*)moleScaleFactor
! 		moleScaleFactor = 100
! 		write(*,*)' Input mole scale factor (around 1-10-100-1000?)'
! 		read(*,*)moleScaleFactor

	
		read(41,*)timeStep
! 		write(*,*)' Current time step = ',timeStep
! 		write(*,*)'Input value for timeStep '
! 		read(*,*)timeStep
					
! 		write(*,*)'Do you want to calculate segments? 0 = no; 1 = yes'
! 		read(*,*)calcSegments
		calcSegments = 1		
! 		write(*,*)' Do you want to move nodes and segments? 0 = no; 1 = yes'
! 		read(*,*)MoveNodesandSegs
		moveNodesandSegs = 1		
! 		write(*,*)'Do you want NodeCapture to be on? 0 = no; 1 = yes'
! 		read(*,*)NodeCaptureOn
		NodeCaptureOn = 1	
	
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
		write(76,*)'---------------------------------------------'
		!close(76)	! stuff file
		
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
!		go to 10

		end do	! end do numPolySteps

		close(75)
		close(76)

10	return		! either we are done or we are bailing out for some reason

	end
	