	Subroutine GridBulkComposition()
	implicit none
	include "Diffuse_GB_Hex.inc"
	include "Assemb.inc"
	integer*4 i,j,k,iXl,ID
	real*8 MIFPhaseArea(maxPhases),MIFPhaseMoles(maxPhases),sum,sumArea

	write(12,*)' Calculate bulk composition from grid file'
	call CrystalAreaRoutine()
	write(12,*)' Crystal     Area'
	do iXl = 1,numCrystals
		write(12,*)iXl,CrystalArea(iXl)
		end do
	do i = 1,numPhMIF
		MIFPhaseArea(i) = 0.0d0
		end do
	do iXl = 1,numCrystals
		ID = CrystalMIFID(iXl)			
		MIFphaseArea(ID) = MIFPhaseArea(ID) + CrystalArea(iXl)		
		end do
	! we should now have the total area for each phase in the grid, indexed to MIFID
	! Convert to moles and calculate mode
	sumArea = 0.0d0
	do i = 1,numPhMIF
		MIFphaseMoles(i) = MIFPhaseArea(i)/Vmol(i)
		sumArea = sumArea + MIFPhaseArea(i)
		end do
	do i = 1,numPhMIF
		mode(i) = 100.0d0*MIFPhaseArea(i)/sumArea
		end do
	write(12,*)' '
	write(12,*)' MIFID     Phase           PhaseArea      PhaseMode       PhaseMoles'
	do i = 1,numPhMIF
		write(12,320)i,phName(i),MIFPhaseArea(i),mode(i),MIFPhaseMoles(i)
320		format(I8,2x,A12,3F15.5)
		end do
	! Now calculate the bulk composition
	!     Compute number of moles of each system component
	do i=1,nc
		wtpct(i)=0.D0
		moles(i)=0.D0
		Do k = 1,numPhMIF
			do j=1,numPhCo(K)
				moles(i) = moles(i) + MIFPhaseMoles(K)*xPhCo(k,j)*comp(k,j,i)
				end do
			end do
		end do
	!     convert moles to wt % oxides
	sum=0.D0
	do i=1,nc
		wtpct(i)=moles(i)*molwt(i)/NumCatInOxide(i)
		sum=sum + wtpct(i)
		end do
	!if(sum.le.0.)go to 140
	sum=100.D0/sum
	do i=1,nc     
		wtpct(i)=wtpct(i)*sum
		end do
	write(12,*)'  '
	write(12,*)'  '
	write(12,*)' Bulk composition - calculated from mineral composition, crystal area, and moles'
	write(12,321)(coname(i),i=1,nc)
321	format(14x,24a8)
!      	compute and print delta moles of system components
!      	write(12,322)(1000.*moles(i),i=1,nc)
!322   	format(' mMoles   ',24F8.3)
      	write(12,323)(wtpct(i),i=1,nc)
323   	format(' Wt%      ',24F8.3)

	return
	end
	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 	Subroutine MolesToAreaRoutine()
	! Routine to convert a ∆moles from the MDF routine into a ∆area for the phase
	! It uses the molar volume (J/bar*10 = cm^3/mole)
	implicit none
	include "Assemb.inc"
	include "Diffuse_GB_Hex.inc"

	integer*4 ph_i,MIFID_i,iSeg,iNode,i,j
	real*8 molarV_i
	
	do iNode = 1,numNodes
		do i = 1,numNodeReactionPhases(iNode)
			ph_i = nodeReactionPhases(iNode,i)
			MIFID_i = nodeMIFID(iNode,ph_i)
			! here I should calculate the molar volume based on the composition
			! But it's not a large change from a generic value
			molarV_i = Vmol(MIFID_i)
			nodeAreaDelta(iNode,ph_i) = molarV_i*nodePhaseMolesDelta(iNode,ph_i)
		!	write(18,*)iNode,i,nodePhaseMolesDelta(iNode,ph_i),nodeAreaDelta(iNode,ph_i)
			end do
		end do

	do iSeg = 1,numSegs
		do j = 1,numSegReactionPhases(iSeg)		! value is 0 (no reactions) or 2 (2 different phases)
			!ph_i = segReactionPhases(iSeg,j)	! this is either 1 or 2 (I think)
			MIFID_i = segMIFID(iSeg,j)
			molarV_i = Vmol(MIFID_i)
!			do i = 1,numPoints(iSeg)
			do i = numPointStart(iSeg),numPointEnd(iSeg)
				pointAreaDelta(iSeg,i,j) = molarV_i*pointPhaseMolesDelta(iSeg,i,j)	
		!		write(18,*)iSeg,j,i,pointPhaseMolesDelta(iSeg,i,j),pointAreaDelta(iSeg,i,j)
				end do
			end do
		end do
	!close(18)
 	return
 	end	



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 	Subroutine CrystalAreaRoutine()
!	routine to determine the area of each crystal
!	The method is to sum up the areas of each triangle from the center to the side of the crystal
!	The crystals are not all regular hexagons, so we can't use the formula for a hexagon
!	The area of each triangle is calculated using Heron's formula:
!
!	A = squareRoot (s*(s-a)*(s-b)*(s-c)) where
!	s = (a+b+c)/2
! 	the sides of the triangle are a, b, c
!
!	the sides a, b, c are determined as the distances 
!	a = center to node 1
!	b = center to node 2
!	c = node 1 to node 2
!
!	The array CrystalNodes(maxCrystals,6) has node numbers clockwise from lower left
!	numCrystalNodes(maxCrystals) is 4 or 6 (the number of nodes for each crystal

!	NOTE -- this routine uses only NODEs to calculate area and thus assumes
!		(1) The crystals are all convex outward
!		(2) The segments are all straight lines
!		It does NOT calculate areas using segments.
 	implicit none
 	include "Diffuse_GB_Hex.inc"
 	integer*4 iXl,i,node1,node2
 	real*8 Atemp,Btemp,Ctemp,a,b,c,s,areaTemp,area,Xcenter,Ycenter
	
	do iXl = 1,numCrystals
		area = 0.0d0
		Xcenter = CrystalCenterX(iXl)
		Ycenter = CrystalCenterY(iXl)
		do i = 1,numCrystalNodes(iXL)-1
			node1 = CrystalNodes(iXl,i)
			node2 = CrystalNodes(iXl,i+1)
			Atemp = (Xcenter - nodeX(node1))**2 + (Ycenter - nodeY(node1))**2
			Btemp = (Xcenter - nodeX(node2))**2 + (Ycenter - nodeY(node2))**2
			Ctemp = (nodeX(node1) - nodeX(node2))**2 + (nodeY(node1) - nodeY(node2))**2
			a = sqrt(Atemp)
			b = sqrt(Btemp)
			c = sqrt(Ctemp)
			s = (a+b+c)/2.0d0
			areaTemp = sqrt(s*(s-a)*(s-b)*(s-c))
			area = area + areaTemp
			end do
	!	Now close the loop with the last triangle
		node1 = node2			! this will be the last one we did
		node2 = CrystalNodes(iXl,1)	! this is the first node
		Atemp = (Xcenter - nodeX(node1))**2 + (Ycenter - nodeY(node1))**2
		Btemp = (Xcenter - nodeX(node2))**2 + (Ycenter - nodeY(node2))**2
		Ctemp = (nodeX(node1) - nodeX(node2))**2 + (nodeY(node1) - nodeY(node2))**2
		a = sqrt(Atemp)
		b = sqrt(Btemp)
		c = sqrt(Ctemp)
		s = (a+b+c)/2.0d0
		areaTemp = sqrt(s*(s-a)*(s-b)*(s-c))
		area = area + areaTemp
		CrystalArea(iXl) = area
		end do
	return
	end	
 	
