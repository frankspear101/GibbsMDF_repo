!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine DumpGrid(iunit)
!	Routine to write out entire grid to logical unit iunit
	implicit none
	include "Diffuse_GB_Hex.inc"
	integer*4 iunit,i,j,N1,N2,index,iXl,k,iNode


	write(iunit,*)'----------------------------'
	write(iunit,*)'segCaptureExclusion for each segment'
	do i = 1,numSegs
		write(iunit,*)i,segCaptureExclusion(i)
		end do

	write(iunit,*)' ----'
	write(iunit,*)'numCrystals = ',numCrystals
	write(iunit,*)' ----'

	write(iunit,*)' '
	write(iunit,*)' List all phases for all crystals'
	write(iunit,*)' k MIDFID   Name'
	do k = 1,numCrystals
!		write(iunit,*)k,CrystalMIFID(k),PhaseName(CrystalMIFID(k))
		write(iunit,*)k,CrystalMIFID(k),CrystalPhaseName(k)
		end do

!----------------------------

	write(iunit,*)' '
	write(iunit,*)' Phases at each node - '
	write(iunit,*)' --------------------'
	write(iunit,*)' Phases around each node (all 3)'
	write(iunit,*)' iNode    NodeMIFID'
	do iNode = 1,numNodes
		write(iunit,*)iNode,'    ',(NodeMIFID(iNode,i),i=1,3)
		end do
	write(iunit,*)' --------------------'
	write(iunit,*)' Unique phases around each node'
                    !       1          5      11      60       2       7       2
	write(iunit,*)'  Node       MIF IDs                Reaction phases'
	write(iunit,227)iNode,(nodeMIFID(iNode,j),j=1,numCrystalsAtNode(iNode)),			&
				(nodeReactionPhases(iNode,j),j=1,numNodeReactionPhases(iNode)),nodeSegTheSamePhases(iNode)
227	format(I6,'   ',20I8)



	write(iunit,*)'=================='
	write(iunit,*)'Number of nodes = ',numNodes
	do i = 1,numNodes
		write(iunit,101)i,nodeX(i),nodeY(i)
		end do
101	format(I5,2F8.3)

	write(iunit,*)' '
	write(iunit,*)'   node    connect        OnTop    connect        OnTop    connect'
	write(iunit,*)'  index     node    seg   point     node    seg   point     node    seg   point'
	do index = 1,numNodes
		write(iunit,77)index,(nodeNodeConnect(index,j),nodeSegConnect(index,j),nodePointOnTop(index,j),j = 1,numNodeSegs(index))
77		format(15I8)
		end do

	write(iunit,*)'======================='
	write(iunit,*)'Segments'
	write(iunit,*)'Number of segments = ',numSegs
	write(iunit,*)'      j      node1    node2    length  dXgrid  numPts_in_Seg   pointX/Y'
		  !       1       1       2       1.000      11
	do j = 1,numSegs
		N1 = segNodes(j,1)
		N2 = segNodes(j,2)	
		write(iunit,65)j,segNodes(j,1),segNodes(j,2),segLength(j),segdXGrid(j),numPointStart(j),numPointEnd(j),  &
				NodeX(N1),(pointX(j,i),i=numPointStart(j),numPointEnd(j)),nodeX(N2)
		write(iunit,66)                                                                     &
				NodeY(N1),(pointY(j,i),i=numPointStart(j),numPointEnd(j)),nodeY(N2)
		write(iunit,*)' '
		end do
65	format(3I8,2F12.3,I8,50F8.3)
66	format(T57,50F8.3)

	write(iunit,*)'    '
	write(iunit,*)'    '
	write(iunit,*)'======================='
	write(iunit,*)'Crystals'
	write(iunit,*)'Number of crystals = ',numCrystals
	do iXl = 1,numCrystals
!		write(iunit,81)iXl,CrystalPhaseName(iXl),CrystalArea(iXl)
		write(iunit,81)iXl,CrystalPhaseName(iXl)
		end do
81	Format(I5,A16,4x,2F12.5)

	write(iunit,*)' '
	write(iunit,*)' Crystal centers, nodes and segments'
	write(iunit,*)'   i           X            Y             Nodes               Segments'
	Do index = 1,numCrystals
		CrystalCenterX(index) = PolyCenterX(index)
		CrystalCenterY(index) = PolyCenterY(index)
		write(iunit,1066)index,CrystalCenterX(index),CrystalCenterY(index),			&
			(CrystalNodes(index,j),j=1,numCrystalNodes(index)),(CrystalSegs(index,j),j=1,numCrystalNodes(index))
1066		format(I5,2F12.5,40I5)
		end do
! 	write(iunit,*)'Compositions of grain boundary at nodes and segments'
! 	write(iunit,*)' Nodes '
! 	do iNode = 1,numNodes
! 		write(iunit,*)iNode,(NodeComp(iNode,j),j=1,numPhCo(1))
! 		end do
! 	write(iunit,*)' Segments and points '
! 	do iSeg = 1,numSegs
! 		do k = 1,numPoints(iSeg)
! 			write(iunit,*)iSeg,k,(pointComp(iNode,j),j=1,numPhCo(1))
! 			end do
! 		end do


	
	close(iunit)
	return
	end




