!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Original module written by Julien Lesgourgues for using weak lensing data in CosmoMC        !
! It was written for the COSMOS 3D data by Massey et al. [astro-ph/0701480],                  !
! was used in [arXiv:0705.0533] by J. Lesgourgues, M. Viel, M.G. Haehnelt, R. Massey          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This modified version was written by Jason Dossett.                                         !
! It allows the use to the CFHTLens weak lensing tomography data of [arXiv:1303.1808]         !
! and was first used in [arXiv:1501.03119] by J. Dossett, M. Ishak, D. Parkinson, and T. Davis!
! to probe deviations from General Relativity.                                                !
! Last Modified 02/2015                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module CFHTLens
    use settings
    use CosmologyTypes
    use CosmoTheory
    use likelihood
    use WeakLen_Common ! Common module includes function converting harmonic power spectrum into angular correlation functions
    use Calculator_Cosmology
    use Likelihood_Cosmology
    implicit none
    private
    
    type, extends(TCosmologyWLLikelihood) :: CFHTLensLikelihood
    contains
    procedure :: ReadIni => CFHTLensReadIni
    procedure :: LogLike => CFHTLensLnLike
    procedure :: ArrangeForLike =>CFHTLens_ArrangeForLike
    end type CFHTLensLikelihood
  
    logical :: use_CFHTLens=.false.
    
    public CFHTLensLikelihood_Add
  
    contains
  
    subroutine CFHTLensLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    Type(TSettingIni) :: ini
    Type(CFHTLensLikelihood), pointer :: this

    use_CFHTLens = Ini%Read_Logical('use_CFHTLens',.false.)
    if (.not. use_CFHTLens) return
    
    CosmoSettings%use_WeakLensing = .true.

    allocate(this)
    this%LikelihoodType = 'WeakLen'
    call this%ReadDatasetFile(Ini%ReadFileName('CFHTLens_dataset'))
    call LikeList%Add(this)

    end subroutine CFHTLensLikelihood_Add

    subroutine CFHTLensReadIni(this,Ini)
    use MatrixUtils
    class(CFHTLensLikelihood) this
    class(TSettingIni) :: Ini
    real(mcp), allocatable :: d_read(:,:)  !temp for correlation vector
    real(mcp) theta_dum,z_dum
    integer zix,k,ios
    character(LEN=:), allocatable :: prename,name
    character(LEN=3) :: number
    Type(TTextFile) :: F

    if(feedback>0)write(*,*)'Initializing CFHTLens WL data set: '//trim(this%name)
    
    this%nlmax     = 130
    this%xstop     = 200._mcp
    this%dlnl      = 0.1
    this%num_z     = 20
    this%kmax      = 100    
    this%nbinmax   = Ini%Read_Int('num_tomography_bins')
    this%nzmax     = Ini%Read_Int('size_galaxy_zdist')
    this%nthetamax = Ini%Read_int('num_theta')
    this%ncl       = this%nbinmax*(this%nbinmax+1)/2
    
    call this%Initialize()
    
    this%use_IA =  Ini%Read_Logical('use_IA',.false.)    
    CosmoSettings%use_IA = this%use_IA .or. CosmoSettings%use_IA
    if(this%use_IA) call this%loadParamNames(Ini%ReadFileName('CFHTLens_paramnames'))    
    
    !Read galaxy distribution files
    prename = Ini%ReadFileName('redshift_dist_prefix_foreground')
    do k=1,this%nbinmax
        write(number,'(I0)')k
        name=trim(prename)//'_'//trim(number)//'.hist'
        call F%open(trim(name))
        zix=1
        do while (zix<=this%nzmax)
            read(F%unit,*,IOSTAT=ios) this%z(zix), this%etaz_f(zix,k)
            if(IS_IOSTAT_END(ios)) exit
            if(ios.ne.0) cycle
            zix=zix+1
        end do
        call F%close()
    end do    
    this%max_z = this%z(zix-1)
    
    !Read galaxy distribution files
    prename = Ini%ReadFileName('redshift_dist_prefix_background',NotFoundFail = .false.)
    if(prename =='') then 
        this%etaz_b = this%etaz_f
    else
		do k=1,this%nbinmax
			write(number,'(I0)')k
			name=trim(prename)//'_'//trim(number)//'.hist'
			call F%open(trim(name))
			zix=1
			do while (zix<=this%nzmax)
				read(F%unit,*,IOSTAT=ios) z_dum, this%etaz_b(zix,k)
				if(IS_IOSTAT_END(ios)) Exit
				if(ios.ne.0) cycle
				zix=zix+1
			end do
			call F%close()
		end do
	end if
    
    
    !read data vector
    allocate(d_read(2*this%nthetamax,this%ncl))
    name=Ini%ReadFileName('measurements_file')
    call F%open(trim(name))
    do k=1,2*this%nthetamax
        if(k .le. this%nthetamax) then
            read(F%unit,*)this%theta(k),d_read(k,:)
            this%theta(k)=this%theta(k)/60.*pi/180.
        else
            read(F%unit,*)theta_dum,d_read(k,:)
        end if
    end do
    call F%close()

    !rearrange datavector to match covariance matrix
    do k=1,2*this%nthetamax
        do zix=1,this%ncl
            this%d_obs(k+2*this%nthetamax*(zix-1))=d_read(k,zix)
        end do
    end do

    ! Read covariance matrix
    name=Ini%ReadFileName('covmat_file')
    call F%open(trim(name))
    do k=1,2*this%nthetamax*this%ncl
        read(F%unit,*)this%invcov(k,:)
    end do
    call F%close()
    !get inverse covariance matrix
    call Matrix_Inverse(this%invcov)
    
    this%ah_rescale = Ini%Read_Logical('ah_rescale',.false.)
    if(this%ah_rescale) then
    	this%ah_realizations = Ini%Read_Double('ah_realizations')
    	!AH Correction factor as seen in Heymans et al arXiv:1303.1808
    	!Section 3.3.1.    	
    	this%invcov = (this%ah_realizations-2*this%nthetamax*this%ncl-2._mcp)/&
    	(this%ah_realizations-1.)*this%invcov
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! integrate eta(z) over z (in view of normalizing it to one) !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    this%eta_norm_f(:)=0.
    this%eta_norm_b(:)=0.
    do zix=2,this%nzmax
        this%eta_norm_f(:)=this%eta_norm_f(:)+0.5*(this%etaz_f(zix-1,:)+this%etaz_f(zix,:))*(this%z(zix)-this%z(zix-1))
        this%eta_norm_b(:)=this%eta_norm_b(:)+0.5*(this%etaz_b(zix-1,:)+this%etaz_b(zix,:))*(this%z(zix)-this%z(zix-1))
    end do
    
    end subroutine CFHTLensReadIni
  
    function CFHTLensLnLike(this,CMB,Theory,DataParams)
    Class(CFHTLensLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) ::  DataParams(:)
    real(mcp) :: CFHTLensLnLike
    !loop indices
    integer zix,zix2
    real(mcp) CFHTLens_A

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Various steps for Likelihood computation start here !!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!! JD: A_CFHTLens for intrinsic alignment correction  !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CFHTLens_A = 0.
    if(this%use_IA) CFHTLens_A=DataParams(1)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Compute lensing Cross-power spectra C_l^shear !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call this%GetLensingCls(CMB,Theory,CFHTLens_A)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Compute correlation function C_plus,minus(theta) !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    call this%Cl_to_Ctheta()
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Arrange all C(theta) values in a vector, following the same order as in the data !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call this%ArrangeForLike()

    if (Feedback>2) write(*,'("CFHTLens nuisnace parameter A = ",F10.5)')CFHTLens_A
    ! Compute the chi^2
    CFHTLensLnLike = 0.
    do zix=1,2*this%nthetamax*this%ncl
        do zix2=1,2*this%nthetamax*this%ncl 
            CFHTLensLnLike=CFHTLensLnLike+&
            (this%d_theo(zix)-this%d_obs(zix))*this%invcov(zix,zix2)&
            *(this%d_theo(zix2)-this%d_obs(zix2))
        end do
    end do

    ! compute -lnL = chi2/2
    CFHTLensLnLike=CFHTLensLnLike/2._mcp

    if (Feedback >1) write(*,'("CFHTLens -LnLike=",F10.5)')CFHTLensLnLike

    end function CFHTLensLnLike
    
    subroutine CFHTLens_ArrangeForLike(this)
    class(CFHTLensLikelihood) this
    integer zix, nl
    
    do zix=1,this%nthetamax
        do nl=1,this%ncl
            this%d_theo(zix+2*this%nthetamax*(nl-1))=this%cplus(zix,nl)
            this%d_theo(zix+this%nthetamax+2*this%nthetamax*(nl-1))=this%cminus(zix,nl)
        end do    
    end do

    end subroutine CFHTLens_ArrangeForLike

    end module CFHTLens

