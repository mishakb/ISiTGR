    module lrg_2dCl
    use ObjectLists
    use settings
    use Interpolation
    use CosmologyTypes
    use CosmoTheory
    use likelihood
    use Calculator_Cosmology
    use Likelihood_Cosmology
    implicit none
    private

    integer, parameter :: maxnz = 1000

    type ZStruct
        integer :: nz
        real, allocatable :: zmin(:), zmax(:)
        real, allocatable :: zdist (:,:)
    end type ZStruct

    type, extends(TCosmoCalcLikelihood) :: TCosmologyISWLikelihood
        integer :: ndat !number of cl datapoints
        real(mcp), allocatable :: ISW_cls_obs(:),ISW_invcov(:,:)
        !clinfo(:,1)-> galaxy distribution; clinfo(:,2)-> lsample index
        integer, allocatable :: clinfo(:,:)
        !CL_theory 1st index -> galaxy distribution, 2nd index -> l-index
        real(mcp), allocatable :: Cl_theory(:,:)
        Type(TIntegerList) :: lsample
        Type(ZStruct) :: zd
        integer :: Nsample  !number of galaxy samples
    contains
    procedure :: lrg_rdzdist
    procedure :: lrg_compute2dCl_Limber
    end type  TCosmologyISWLikelihood

    public TCosmologyISWLikelihood, maxnz

    contains

    ! Read in a redshift distribution
    subroutine lrg_rdzdist(this, fn)
    Class(TCosmologyISWLikelihood) :: this
    character(*), intent(in) :: fn
    Type(TTextFile) :: F
    integer :: i, iot
    real(mcp) :: zm, z(9), dz

    this%zd%zmin = 0.0
    this%zd%zmax = 0.0
    this%zd%zdist = 0.0
    dz = 0.01          !the redshift bin size is hardcoded.
    
    call F%open(trim(fn))
    iot = 0
    i = 0
    do while (iot .eq. 0)
        read (f%unit,*,iostat=iot) zm, z(:)
        if(IS_IOSTAT_END(iot)) exit
        if (iot .ne. 0) cycle
        i = i+1
        this%zd%zmin(i) = zm- (dz/2.0)
        this%zd%zmax(i) = this%zd%zmin(i) + dz
        this%zd%zdist(:,i) = z
        if (i .gt. maxnz) stop 'ERROR - lrg_rdzdist: maxnz exceeded'
    end do
    call F%close()
    if (i .eq. 0) stop 'ERROR : reading file in lrg_rdzdist'
    this%zd%nz = i

    end subroutine lrg_rdzdist

    !--------
    ! Code to compute 2D power spectra given a 3D power spectrum
    !--------
    subroutine lrg_compute2dCl_Limber(this,CMB,Theory,TCMB)
    Class(TCosmologyISWLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp), intent(in) :: TCMB
    real(mcp) :: clval(this%Nsample), w(this%Nsample)
    real(mcp) :: rc, kh, z0, zold, dz, pkz0
    integer :: iz, il, larr
    
    
    This%Cl_theory = 0.
    do il = 1, this%lsample%count
        larr = this%lsample%Item(il)
        clval= 0.0
        w = 0.0
        zold = (this%zd%zmin(1)+this%zd%zmax(1))/4.0
        do iz = 1, this%zd%nz
            z0 = (this%zd%zmin(iz)+this%zd%zmax(iz))/2.0
            dz = z0-zold
            rc = this%Calculator%f_k(Theory%R%value(z0))
            kh = (larr+0.5)/rc/CMB%h   ! k=(l+0.5)/r in units h/Mpc
            pkz0 = Theory%P_ISW%PowerAt(kh,z0)
            w(:) = this%zd%zdist(:,iz)
            w = w*pkz0
            clval = clval +  w*dz            
            zold = z0
        end do
        this%Cl_theory(:,il) = clval(:)/((larr+0.5)**2)*TCMB
    end do
    !want Cls in microK
    this%Cl_theory = this%Cl_theory*1.d6

    end subroutine lrg_compute2dCl_Limber

    end module lrg_2dCl
