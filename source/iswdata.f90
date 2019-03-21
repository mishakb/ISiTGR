! module for storing data from ISW
! real_Cl_file(galaxy distribution index,medianl,Cl,Cl_err) &
! inv_covariance_matrix(in same indexing)
! bdndz_file(to be output by bdndz code);
! NOTE: order of the sample: 2MASS0 2MASS1 2MASS2 2MASS3 LRG0 LRG1 QSO0 QSO1 NVSS'
!
! Module originally written by Shirley Ho
!
! This version of the code has been greatly modified by Jason Dossett for
! the public code ISiTGR
!
! Last Modified 02/2015

	module iswdata
    use ObjectLists
    use settings
    use Interpolation
    use CosmologyTypes
    use CosmoTheory
    use likelihood
    use Calculator_Cosmology
    use Likelihood_Cosmology
    use lrg_2dcl
    implicit none
    private

    type, extends(TCosmologyISWLikelihood) :: ISWLikelihood
    contains
    procedure :: ReadIni => ISWReadIni
    procedure :: LogLike => ISWLnLike
    end type ISWLikelihood

    logical :: use_ISW
    character(Len=:), allocatable :: baserootISW, zfn, cosmoparamfile, command_system

    public baserootISW, ISWLikelihood_Add

    contains

	subroutine ISWLikelihood_Add(LikeList, Ini)
	class(TLikelihoodList) :: LikeList
	Type(TSettingIni) :: ini
	Type(ISWLikelihood), pointer :: this

	use_ISW = Ini%Read_Logical('use_ISW',.false.)
	if (.not. use_ISW) return	
	
	CosmoSettings%use_ISW = .true.

	allocate(this)
	this%LikelihoodType = 'ISW'
	this%needs_powerspectra = .true.
	call this%ReadDatasetFile(Ini%ReadFileName('ISW_dataset'))
	call LikeList%Add(this)
	
	end subroutine ISWLikelihood_Add
	
	subroutine ISWReadIni(this,Ini)
	class(ISWLikelihood) this
    class(TSettingIni) :: Ini
    Type(TTextFile) :: F
	character(Len=:), allocatable :: realcl_file,invcov_file
    character(Len=:), allocatable :: templine
	character(Len=256) cwd, tmpdir
	integer:: i,j,io2
	real(mcp) error
	integer, allocatable :: ltmp(:)

	if (feedback>0) write(*,*) 'Initializing ISW-galaxy cross correlation data: '//trim(this%name)
	
	this%max_z = Ini%Read_Double('ISW_max_z',5.d0)
    this%num_z = Ini%Read_Int('ISW_num_z',50)
	this%ndat = Ini%Read_Int('num_datapoints',42)
	this%NSample = Ini%Read_Int('num_galaxy_samples',9)
	
	!allocate zd data type arrays
	allocate(this%zd%zdist(this%Nsample,maxnz))
	allocate(this%zd%zmin(maxnz),this%zd%zmax(maxnz))
	
	allocate(this%ISW_cls_obs(this%ndat),ltmp(this%ndat),this%clinfo(this%ndat,2))		
	!Reading ISW datafiles
	!---reading in the real Cl file
	!realcl_file = 'data/ISWHo_Cl_k_0p05.dat'
	realcl_file = Ini%ReadFileName('realcl_file')
	call F%open(realcl_file)
	io2 = 0
	i = 1
	do while(i<=this%ndat)
		read (F%unit,*,iostat=io2) this%clinfo(i,1),ltmp(i), this%ISW_cls_obs(i), error
		if(IS_IOSTAT_END(io2)) exit
		if (io2 .ne. 0) cycle
		call this%lsample%Add(ltmp(i))
		i=i+1
	end do
	call F%close()
	if(maxval(this%clinfo(:,1))>this%Nsample) call MpiStop("ISWReadIni: Galaxy dist number too high")
	
	!sort l values & remove duplicates
	call this%lsample%sort()
	call this%lsample%RemoveDuplicates()
	!index lvalues into this%clinfo(:,2)
	do i=1, this%ndat
		j=1
		do while((ltmp(i) - this%lsample%Item(j))/=0)
			j=j+1
			if(j>this%lsample%count) call MpiStop("ISWReadIni: could not find lsample index")
		end do
		this%clinfo(i,2) = j
	end do
	allocate(this%Cl_theory(this%Nsample,this%lsample%count))	
	
	!---- reading in the covariance matrix
	allocate(this%ISW_invcov(this%ndat,this%ndat))
	!invcov_file = 'data/ISWHo_invcov_k_0p05.dat'
	invcov_file = Ini%ReadFileName('invcov_file')
	call F%open(invcov_file)
	do i = 1, this%ndat
		read (F%unit,*) this%ISW_invcov(i,:)
	end do
	call F%close()
	
	!Set up commands for later
	call getcwd(cwd)
	!By default store files in user specified environmental variable TMPDIR
	call get_environment_variable("TMPDIR",tmpdir)
	!If $TMPDIR doesn't exist store in bdndz_code/IO
	if (trim(tmpdir)=='')then
		tmpdir=trim(cwd)//'/bdndz_code/IO'
	end if

	!prepare names of files and commands for bdndz code
	cosmoparamfile = trim(tmpdir)//'/'//trim(baserootISW)//'.dat'
	zfn = trim(tmpdir)//'/bdndz.'//trim(baserootISW)//'.dat'
	
	if(CosmoSettings%ISiTGR_BIN) then
		templine = 'cd '//trim(cwd)//'/bdndz_code/; ./dndz_bin.x'
	else if(CosmoSettings%ISiTGR_scale_dep) then
		templine = 'cd '//trim(cwd)//'/bdndz_code/; ./dndz_func_scale.x'
	else
		templine = 'cd '//trim(cwd)//'/bdndz_code/; ./dndz_func.x'
	end if	

	command_system =trim(templine)//' '//trim(zfn)//' '//trim(cosmoparamfile)//'; cd ..'

	end subroutine ISWReadIni

	real(mcp)  function IswLnLike(this,CMB,Theory,DataParams)
	Class(ISWLikelihood) :: this
	Class(CMBParams) CMB
	Class(TCosmoTheoryPredictions), target :: Theory
	Type(TTextFile) :: F
	real(mcp) DataParams(:)
	real(mcp), allocatable :: ISW_cls_theory(:)
	integer :: cosmotype,i,j
	real(mcp) :: omegam, TCMB
	real(mcp) :: dr_0, dr_inf
	integer rfunc
	!JD timing Failsafe variables ----------------
	integer, Dimension(13) ::buff
	integer :: status, cnt, rtime, zmtime, sz
	integer time, stat	
	
	rtime = time()
	call F%CreateFile(trim(cosmoparamfile))
	omegam = CMB%omb + CMB%omc+ CMB%omnu
	TCMB = 2.726 ! This is set and used only in the ISW likelihood code, and it is
	! consistent with what cosmomc uses in other parts of the code

	!writing to cosmoparam file that is read by the bdndz code
	cosmotype = 1
	write (F%unit,'(I9)') cosmotype
	write (F%unit,'(E20.10)') TCMB
	write (F%unit,'(E20.10)') CMB%omb
	write (F%unit,'(E20.10)') CMB%omnu
	write (F%unit,'(E20.10)') omegam
	write (F%unit,'(E20.10)') CMB%omv
	write (F%unit,'(E20.10)') CMB%w
	!currently the isw likelihood code doesn't support w_a (running of dark energy
	!equation of state), however the bdndz code does.
	write (F%unit,'(E20.10)') CMB%wa  !w_a
	write (F%unit,'(E20.10)') CMB%h
	write (F%unit,'(E20.10)') Theory%sigma_8
	write (F%unit,'(E20.10)') CMB%InitPower(2)  ! ns, scalar spectral index
	write (F%unit,'(E20.10)') CMB%InitPower(3)  ! running of scalar spectral index
	if(CosmoSettings%ISiTGR_BIN)then
		!JD ISiTGR_BIN Parameters
		write (F%unit,'(E20.10)') CMB%TGR_Q1
		write (F%unit,'(E20.10)') CMB%TGR_D1
		write (F%unit,'(E20.10)') CMB%TGR_Q2
		write (F%unit,'(E20.10)') CMB%TGR_D2
		write (F%unit,'(E20.10)') CMB%TGR_Q3
		write (F%unit,'(E20.10)') CMB%TGR_D3
		write (F%unit,'(E20.10)') CMB%TGR_Q4
		write (F%unit,'(E20.10)') CMB%TGR_D4
		write (F%unit,'(E20.10)') CosmoSettings%ISiTGR_zdiv
		write (F%unit,'(E20.10)') CMB%TGR_kc
		if(CosmoSettings%ISiTGR_true_bin)then
			write (F%unit,'(I9)') 1
		else
			write (F%unit,'(I9)') 0
		end if
	else
	    !JD ISiTGR Parameters
	    if(CosmoSettings%ISiTGR_Rfunc)then
			dr_0 = CMB%TGR_R0
        	dr_inf = CMB%TGR_Rinf
        	rfunc = 1
   	 	else
    		dr_0 = CMB%TGR_D0
        	dr_inf = CMB%TGR_Dinf
        	rfunc = 0
    	end if
        write (F%unit,'(E20.10)') CMB%TGR_Q0
        write (F%unit,'(E20.10)') dr_0
        if(CosmoSettings%ISiTGR_scale_dep)then
            write (F%unit,'(E20.10)') CMB%TGR_Qinf
            write (F%unit,'(E20.10)') dr_inf
            write (F%unit,'(E20.10)') CMB%TGR_kc
        end if
        write (F%unit,'(E20.10)') CMB%TGR_s
        write (F%unit,'(I9)') CMB%TGR_tdep
        write (F%unit,'(I9)') rfunc
    end if
	call F%close()

	status = stat(trim(zfn),buff)
	zmtime = buff(10)
	call system(trim(command_system))


	!JD Timing Failsafe make sure zdistribution has been written completely
	cnt = 0
	do While (cnt<3)
		status = stat(trim(zfn),buff)
		if (status == 0 .and. zmtime<buff(10) .and. buff(8)==sz) then
			cnt = cnt + 1
		end if
		sz = buff(8)
		if(time()>(rtime+40))then
			IswLnLike = logZero
			return
		end if
	end do
	!End JD timing failsafe

	!-------- Read in zdistribution
	call this%lrg_rdzdist(trim(zfn))
	!Get the Theoretical Cl's
	call this%lrg_compute2dCl_Limber(CMB,Theory,TCMB)

	!Arrange Cl's in to ISW_CL_theory (to match observation vector shape)
	allocate(ISW_Cls_theory(this%ndat))
	do i=1, this%ndat
		ISW_Cls_theory(i) = This%Cl_theory(this%clinfo(i,1),this%clinfo(i,2))
	end do

	!Compute the chi^2
	ISWLnLike = 0.
	do i=1,this%ndat
		do j=1, this%ndat
			ISWLnLike = ISWLnLike &
			+(ISW_Cls_theory(i)-this%ISW_Cls_obs(i))*this%ISW_invcov(i,j)&
			*(ISW_Cls_theory(j)-this%ISW_Cls_obs(j))
		end do
	end do

	ISWLnLike = ISWLnLike/2.

	if (Feedback > 1) write(*,*) trim(this%name)//': ISWLnlike = ', ISWLnLike

	end function IswLnLike

    end module iswdata
	