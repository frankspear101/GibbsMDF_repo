	Module MyCanvas
	use AWE_Interfaces
	TYPE(AWE_Canvas) :: XYPlot		! The plotting canvas
	end Module
      PROGRAM GIBBS
	USE AWE_INTERFACES
	use MyCanvas
      implicit none
	TYPE(AWE_Canvas) :: GridPlot		! The plotting canvas
	TYPE(AWE_Canvas) :: SegPlot1
	TYPE(AWE_Canvas) :: SegPlot2
	TYPE(AWE_Canvas) :: SegPlot3
!	TYPE(AWE_Canvas) :: XYPlot

!           Coded by F.S. Spear 2/85
!           Macintosh version written by F. Spear 2/89
!           modified for contour routine by T. Menard 6-9/88
!           modified for open system behavior by F. Spear 3/89
!           modified for lnKeq equations by F. Spear 2/90
!           modified for linear arrays by F. Spear 9/90
!           modified for V=f(P,T) by M.J. Kohn 9/90
!           modified for Margules phases by M.J. Kohn 9/90
!               Kerrick and Jacobs mixed H2O-CO2 fluids
!           modified for Adobe Illustrator output by M.J. Kohn 1/91

!           Version 3.0 involves the major modifications necessary
!               to incorporate multisite mixing F. Spear (January, 1993)
!               v 3.0 also calculates Gibbs matrix by using partial
!                 molar quantities (, S and V) multiplied by
!                 coefficients matrix
!           Version 4.0 incorporates a Newton-Rapheson routine to
!               solve for P-T conditions using enthalpy data
!               plus code to solve for partial derivatives of complex
!               solutions using numerical partial derivatives
!	Version 5 is a complete rewriting of the storage of data.
!		Rather than using linear arrays and pointers (efficient storage but a coding nightmare)
!		Version 5 uses indicies of (phase,phaseComponent) - more storage but much easier to keep track of
!		This version may also eventually use pointers to the "currentAssemblage" rather than
!			loading and reloading every time a new assemblage is tried. This will be much more efficient
!			for things like pseudosections
!	November, 2016. Converted code toe FORTRAN 95 and implemented AWE interfaces
!		This allowed retiring the old Carbon interface routines
!     	Logical unit numbers used are:
!     	1 = thermodynamic data file 
!     	3 = configuration file
! 	23 = Gibbs.PlotDefinitions file (contains X-Y coordinate plot definitions)
!     	4 = RESTORE/SAVE file (*.sav)
!     	5 = input unit number in subroutine BEGIN
!     	26 = PT path input file
!     	9 = system console (on Macintosh)
!     	11 = Graphics window unit number
!     	12 = Output window unit number
!     	14 = New Input File window unit number
!     	15 = DiffGibbs output window unit number
!     	19 = DiffGibbs_D_Coef.dat file unit number
!     	22 = Output of *.in type of save file (SUBROUTINE CHANGE)
!     	51 = Postscript scratch file unit number
!     	21 = Postscript unit number (for saveas)
! c*******************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "SaveRxn.inc"
	include "Singular.inc"
!	include "MakeMyGrid2.inc"
!	include "Tangent.inc"
	include "MatrixSingular.inc"	
	include "Diffuse_GB_Hex.inc"
	INCLUDE "LUTColorStuff.inc"
! c*******************************************
!     LOCAL VARIABLES
	logical exist
	integer*4 ioption,i,j,numBC,R,G,B,status
	character*32 ConfigTitle,rockfile(16)
	character*16 dummy
	character*64 thisThmoFile
	integer*4 skipKFlu
	common /skipswitch/skipKFlu
!  !	**********************************
! 	integer*4 LUTnumColors,LUT_RGB(300)
! 	real*4 LUT_CMYK(300,4)
! 	common /LUTcommon/LUTnumColors,LUT_RGB,LUT_CMYK
! !	**********************************	
      	data iternary,NoTieMinerals/0,1/	! initialize tie line plotting to zero

!     	A LISTING OF THE VARIABLES USED IN THIS PROGRAM IS IN THE
!       GIBBS.TXT FILE INCLUDED WITH THE PROGRAM
!	 open configuration (preferences) file
      	CONFIG='./GibbsMDF_Essentials/GibbsMDF.FIG'
! 	Check to see that the configuration file is present.  If not, display a message
	inquire(file = config, EXIST = exist)
	if(.not.exist)then
		write(*,*)'The configuration file is not present.  Abort the mission'
		pause 'Hit return to continue'
		stop
		endif
      	open(3,file=Config,status='old')
      	read(3,'(A)') ConfigTitle
      	READ(3,'(A)')dummy          	! read the row of dashes
      	call ReadThomoFileList(thisThmoFile)
	thmoFile = thisThmoFile
	write(*,*)' ----------------'
      	READ(3,'(A)')dummy          	! read the row of stars
	read(3,*)numBC	! number of bulkrockcomposition files
	do 10 i = 1,numBC
	read(3,*)rockFile(i)			! temporary storage of file names
	write(*,*)i,rockFile(i)
10	continue
	write(*,*)'Pick the bulk composition file you wish to use'
!	read(*,*)numBC
	numBC = 1
	bulkRockCompositionFile = rockFile(numBC)		! this is the file name used in Sub PickBulkComposition
	bulkRockCompositionFile = './GibbsMDF_Essentials/'//BulkRockCompositionFile		! this is the file name used in Sub PickBulkComposition

      	READ(3,'(A)')dummy          	! read the row of stars
      	READ(3,*)inewton,tolnewton,newtonstep
      	newtonStepMax = newtonStep
! 	inewton=1 uses Newton's method to solve integrated equations; 0 is Gibbs method                                  
! 	Tolnewton is default tolerance for convergence using Newton's method 
      	READ(3,'(A)')dummy          	! read the row of stars
!     	read plotting default values
      	Read(3,*)iBig,iTrip,iAlphaBeta     ! number for default plot in config file
!      	READ(3,'(A)')dummy          	! read the row of stars
!      	READ(3,'(A)')dummy          	! read A HEADER
!	read(3,*)TcStart,TcEnd,Tinc
!	read(3,*)PbStart,PbEnd,Pinc
! 	done reading Gibbs.fig file
      	open(23,file='./GibbsMDF_Essentials/Gibbs.PlotDefinitions',status='old')	! open file with plot definitions (unit 23)
      	READ(23,'(A)')dummy          	! read the header
      	do 215 i = 1,(iBig-1)*7     	! skip the ones we don't need
      	read(23,'(A)',end=218)dummy
215   	continue
      	go to 219         		! No end of file (yet)
218   	continue    			! we hit an end of file - problem!
  !    	call FSS_alert('The default plot number exceeds the number of plot def&
 !    	&initions in the configuration file.  Plot number 1 is being used')
      	rewind(23)			! start over and choose iBig = 1 as default plot
      	READ(23,'(A)')dummy          	! read the header
      	iBig = 1
219   	continue
!     	Now read the plot definition
      	READ(23,'(A)')dummy          	! read the row of stars
      	READ(23,'(a)',end=218)PlTitle
      	READ(23,*,end=218)NXPLT,NYPLT
      	READ(23,'(a)',end=218)XLAB
      	read(23,*,end=218)xor,xmin,xmax,xlen,nxstep,nxdec
      	READ(23,'(a)',end=218)YLAB
      	read(23,*,end=218)yor,ymin,ymax,ylen,nystep,nydec
!     	Open thermodynamic data file
	call OpenThermoFile
! 	Check to see that the file is present.  If not, display a message
! 	HydroPFile = './GibbsMDF_Essentials/HydrostaticP.table'
! 	inquire(file = HydroPFile, EXIST = exist)
! 	if(.not.exist)then
! 		call fss_alert('ALERT!!','The file (HydrostaticP.table) is not present in the Gibbs_Essentials folder.')
! 		stop
! 		else
! ! 		load hydrostatic fluid pressure table
! 		call hpload
! 		endif
! 	Check to see that the file is present.  If not, display a message
! 	DiffGibbsFile = './GibbsMDF_Essentials/DiffGibbs_D_Coeffs.dat'
! 	inquire(file = DiffGibbsFile , EXIST = exist)
! 	if(.not.exist)then
! !		call fss_alert('The file (DiffGibbs_D_Coeffs.dat) is not present in the Gibbs_Essentials folder.&
!  !    & If you try to run a DiffGibbs problem, the program will crash')
! 		endif

!     	initialize variables
      	icont=0
      	numNew=0
      	IOPEN=0
	iCareifMIS = 1		! i care if matrix is singular = yes (print error message)
	call ZeroiLong
      	ilong(1)=1		! error output
      	ilong(11)=1		! echo input file reading new assemblage
	newtonStepsOutput = 0			! flag to PRINTT every Newton Step (only for overstepping calculation)
						! Only used in Subroutine dlnAdX in file Gibbs_Thermocalc.f95
! 	Kfluid variables
	skipKFlu = 0			! allows skipping KFluid questions from AFM routine only
	PFluidSwitch = 0
	KFlu = 0
	RockDensity = 2.7		! gm/cm^3
	bulkCompSwitch = .false.	! No bulk composition is open at start up 
!      	initialize postscript output
!      	SET DEFAULTS FOR MONITOR AND DELTAX
	ndep = 1
      	NSTEP=1
      	SNSTEP=1
      	DO 3001 I=1,16
         SMON(I)=I
         MON(I)=I
	 ncont(i) = i		! contour variable
         DELTAX(I)=0.d0
         SDEL(I)=0.d0
3001  	CONTINUE
!     	set defaults for plotting
      	idraw=1
      	ispace=' '
      	isymb=0
      	nopen=0
      	IPEN=3
! 	Open output window
	OPEN(12,FILE = 'Gibbs6_MDF OutPut',ACCESS = 'window, 800, 900')
	!call AWE_tileWindows()
	call AWE_MoveWindow(12,750,0)
!      	draw initial plot with defaults
!      		CALL USER(xor,XMIN,XMAX,xlen,yor,YMIN,YMAX,ylen) -- now in Sub PlotAxes
!     		call axis(nxstep,NXDEC,XLAB,nystep,NYDEC,YLAB,PLTITL)-- now in Sub PlotAxes
	ColorFileName = './GibbsMDF_Essentials/GibbsColors.txt'
	call ReadColorFile		! set up AWE and PS colors
	CurrentColor = 1		! default = Black
	CurrentColor = 1		! Black for the axes

!	Open LUT file
! 	call FSS_Alert('Alert','Open LUT file')
 	open(41,file='./GibbsMDF_Essentials/LUT_Jet_Formatted.txt',iostat=status)
!	open(41,file='./GibbsMDF_Essentials/LUT 16 colors.txt',iostat=status)
	if(status.ne.0)then
		call FSS_Alert('Alert','Failure opening LUT file in Gibbs Essentials')
		endif
!	Read the file
	read(41,*)LUTnumColors
	read(41,*)dummy
	do i = 1,LUTnumColors
		read(41,*)j,LUT_RGB(i),R,G,B,(LUT_CMYK(i,j),j=1,4)
		end do

	close(41)
	
! 	set default bounds for contour operations as plotting limits
	xMinBound = xmin
	xMaxBound = xmax
      	if(NYPlt.eq.2.and.ymin.eq.0.0)then
		yMinBound = 0.1			! If Y axis is P=0, set lower P bound to 100 bars because at P less than 100 bars, H2O model is no good
		else
		yMinBound = ymin
		endif
	yMaxBound = ymax
! 	set defaults for ternary plot
	Llab = 'Left'
	Rlab = 'Right'
	Tlab = 'Top'
	!Pltitle = 'Ternary plot'
	scale = 1.0
	itic = 1
	PlotIndex(1,1) = NXPLT
	PlotIndex(1,2) = NYPLT
	PlotIndex(1,3) = 3
	PlotIndex(1,4) = 4
      	isfileopen = 0
	! these are initial values for grain boundary energies.
	! They can be changed in subroutine ThermoData
	HseGB =  10000.  	! I'm not sure what these numbers should be (or even what the units are)
	HdeGB =1000000.

	EqOrMDFProblem = 0	! default is equilibrium problem
	!EqOrMDFProblem = 1	! default is MDF problem
	outputType = 2		! affinity = yes in output file
! ----- MAIN CONTROL MENU ----------------------------

100   CONTINUE
	ioption=1         					! default
	write(*,*)' ***********************************'
	write(*,*)'Thermo file: ',THMOFILE         		! thermodynamic data file
     	If(inewton.eq.1)write(*,*)'Newton''s method'
     	If(inewton.eq.0)write(*,*)'Gibbs'' method'
	WRITE(*,*)' ***********************************'
	WRITE (*,*)' MAIN MENU OPTIONS:'
	WRITE (*,*)'   1 = Begin/save problem'
	WRITE (*,*)'  -----------------------------'
	WRITE (*,*)'   2 = Single steps'
	write (*,*)'   3 = Call Subroutine GB_Diffuse_Hex()'
	write (*,*)'   4 = Run Model'
	write (*,*)'   42= Run Model 2'
	write (*,*)'   43= Run Model 3'
	write (*,*)'   5 = Segment Plot routines'
	WRITE (*,*)'   6 = Grid Plot routines'
	WRITE (*,*)'   7 = Phase comp Plot routines'
	WRITE (*,*)'  -----------------------------'
	WRITE (*,*)'   8 = Go to global menu'
	WRITE (*,*)'   9 = Plotting menu'
	WRITE (*,*)'  10 = Set EqOrMDFProblem variable'
 	WRITE (*,*)'  11 = Thermodynamic data menu'
	WRITE (*,*)'  12 = Change output type (affinity in output file)'
	WRITE (*,*)'  -----------------------------'
	WRITE (*,*)'  -----------------------------'
	WRITE (*,*)' CHOOSE OPTION'
	read(*,*)ioption

	select case(ioption)
	! ----------------- 1 ---------------------
	case(1)
		call GibbsBegin
	! ------------------ 2 --------------------
	case(2)
		if(isfileopen.eq.0)Then
			call fss_alert('ALERT!!','You must open an input file first')
			go to 100
			endif
		CALL STEPS
	! ------------------ 3 --------------------
	case(3)
		call GB_Diffuse_Hex()

	! ------------------ 4 --------------------
	case(4)
		call FSS_Alert('Achtung!','I''m sorry, Dave. I''m afraid I can''t do that')
		!call GB_RunModel(SegPlot1,SegPlot2,SegPlot3)

	case(42)
		call FSS_Alert('Achtung!','I''m sorry, Dave. I''m afraid I can''t do that')
! 		call GB_RunModel_2(SegPlot1,SegPlot2,SegPlot3)
	case(43)
		call GB_RunModel_3(SegPlot1,SegPlot2,SegPlot3)

	! ------------------ 5 --------------------
	case(5)
		call SegPlotRoutines(SegPlot1,SegPlot2,SegPlot3)

	! ------------------ 6 --------------------
	case(6)
		call GridPlotRoutines(GridPlot)

	case(7)
		call PhCompPlotRoutines(SegPlot1,SegPlot2,SegPlot3)

	! ------------------ 8 --------------------
      	case(8)
      		CALL GLOBAL(0)
	!-------------- 10 ------------------
	case(10)
		write(*,*)' Equilibrium Or MDF Problem: 0 = Eq, 1 = MDF'
		write(*,*)' Current value = ',EqOrMDFProblem
		write(*,*)' Input new value'
		read(*,*)EqOrMDFProblem
		if(EqOrMDFProblem.ne.0.and.EqOrMDFProblem.ne.1)then
			call FSS_Alert('Alert','You chose poorly...')
			go to 100
			endif


	!-------------- 11 ------------------
	case(11)
		call thermodata
	case(12)
		write(*,*)'OutputType = ',outputType
		write(*,*)'Input new value for outputType (=1 no affinity or =2 yes affinity)'
		read(*,*)outputType

	case default
		call FSS_Alert('Alert','You chose poorly.....')
		go to 100
	end select

	go to 100
	END

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine GibbsBegin
    	use AWE_Interfaces
	implicit none
!*******************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "SaveRxn.inc"
	include "Singular.inc"
!*******************************************
	integer*4 ibegin,inew,i,status
	character*24 dummy

	ibegin = 3
	go to 10
100	continue
        ibegin=0       					! default 
        WRITE(*,*)' ********************************'
        WRITE(*,*)'   0 = Return'
        WRITE(*,*)'   1 = Select bulk composition from disk file BulkRockAnalyses'
        WRITE(*,*)'   2 = READ reference assemblage (MIF) from disk file'
        WRITE(*,*)'   3 = Read GB model file'
	WRITE(*,*)'  -----------------------------'
        WRITE(*,*)'   5 = SAVE current problem as an input file'
	WRITE(*,*)'  -----------------------------'
        WRITE(*,*)'   6 = SAVE Graphics as PostScript file'
	WRITE(*,*)'  -----------------------------'
        WRITE(*,*)'   9 = SAVE Assemblage file and open a new one'
        WRITE(*,*)' CHOOSE OPTION'
      	read(*,*)ibegin

10	continue
	select case(ibegin)
	! ---------------------------
	case(0)
		return
	! ---------------------------
	case(1)
		call PickBulkComp
	! ---------------------------
	case(2)
            	INEW=0
            	numNew=0
            	i = 1			! this instructs subroutine to open a new file
            	CALL BEGIN(i,INEW)
            	IF(INEW.EQ.1)then       			! if true, we actually started a new problem
              		isfileopen = 1
			do i = 2,numPhMIF
				OSPhase(i) = 1			! set all phases = MDF initially
				end do
			OSPhase(1) = 0			! this is the first phase in the MIF and must be the grain boundary
							! The grain boundary is the only phase that is EQ all others are MDF
			call PrintMinAssemblage

!			call change				! to pick the assemblage
			call AdjustASM(1)			! (1) = iChanged (we changed assemblage)
			!call Printt(1)
              		endif
	! ---------------------------
	case(3)
		call FSS_Alert('Alert','Open GibbsMDF model file')
		open(16,file='',status = 'OLD',iostat = status)
		if(status.ne.0)then
			call FSS_Alert('Alert','Problem opening model file')
			return
			endif
		INQUIRE (16, NAME=modelFile)

		read(16,*)dummy			! title of model file
		read(16,*)filein		! input file name
		write(*,*)'Input MIF file  ',filein
	      	OPEN(5,FILE=filein,STATUS='OLD',action='read',iostat=status)
	      	if(status.ne.0)then
	      		call FSS_Alert('ALERT','Error opening input MIF file')
	      		return
			endif
		
            	INEW=0
            	numNew=0
            	i = 2			! this instructs subroutine that the file is already open
            	CALL BEGIN(i,INEW)
            	IF(INEW.EQ.1)then       			! if true, we actually started a new problem
              		isfileopen = 1
			do i = 2,numPhMIF
				OSPhase(i) = 1			! set all phases = MDF initially
				end do
			OSPhase(1) = 0			! this is the first phase in the MIF and must be the grain boundary
							! The grain boundary is the only phase that is EQ all others are MDF
			call PrintMinAssemblage
			call AdjustASM(1)			! (1) = iChanged (we changed assemblage)
			!call change				! to pick the assemblage
			!call Printt(1)
              		endif

		write(*,*)' ------------------------'
		write(*,*)' Echo of model data file is in unit fort.66'
		write(*,*)' ------------------------'
		call ReadModelFile()
!		call GridMaker()
		call SetInitialComp()
		call InitializePhaseComps()
		call CalcDatTP(TC,PB)
		close(66)
		write(*,*)'  '
		write(*,*)'  Done with GridMaker'
		write(*,*)'  '

		return	
	! ---------------------------
	case(5)
           	if(isfileopen.eq.0)then
       			write(*,*)'You must open a file before it can be saved.'
             		go to 100
             		endif  
           	CALL SAVEIN(0,0)
	! ---------------------------
	case(6)
           	call psopcl(5)    				! save old PostScript scratch file
	! ---------------------------
	case(9)
		close(16)
		open(unit=16,file="",status="UNKNOWN")
		write(16,*)numSysCoInFile+1
      		WRITE(16,1010)(CoNameInFile(i),i=1,numSysCoInFile)
1010  		format('Min                     ',40(4x,A4,2x))
	! ---------------------------
	case default
		go to 100
	end select

	go to 100
	end
	
	
	
