!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module containing common routines for WeakLensing Likelihood functions      !
! Written by Jason Dossett                                                    !
! This module contains updated versions of the routines originally included   !
! in the Angular module written by Julien Lesgourgues which provided simple   !
! routines for inferring angular correlation functions from C_l's             !
! Last Modified 03/20/2014                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module WeakLen_Common
    use settings
    use Interpolation
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    implicit none
    private

    real(mcp), parameter :: dbessel = 0.05_mcp

    type, extends(TCosmoCalcLikelihood) :: TCosmologyWLLikelihood
        !Some Bessels and integration settings
        Type(TCubicSpline) :: Bess0, Bess4
        !Precision Parameters
        integer :: nlmax        ! fixes the maximum required l value in C_l(lensing)
        real(mcp) :: dlnl       ! fixes logarithmic l step for C_l(lensing)
        real(mcp) :: xstop      ! fixes maximum value of x =l*theta in angular correlation function
        ! note: integration of C_l J_0,4(l*theta) over l will stop
        ! at l = min( ll(nlmax),  xstop/theta )
        real(mcp), allocatable :: ll(:)  !array of l values

        !Parameters describing WeakLensing data set:
        integer :: nbinmax    ! Number of redshift bins
        integer :: nzmax      ! Number of z values in galaxy distribution file
        integer :: nthetamax  ! Number of angles theta
        integer :: ncl        ! Number of Cross-power spectra to compute

        !Data Arrays
        real(mcp), allocatable :: z(:),etaz_f(:,:),eta_norm_f(:) ! data's selection functions
        real(mcp), allocatable :: etaz_b(:,:),eta_norm_b(:)
        real(mcp), allocatable :: theta(:)                  ! data array of theta values
        !Arrays for calculating -ln(Like)
        real(mcp), allocatable ::  d_theo(:)    ! Theoretical shear cross-correlations, data arranged like a vector
        real(mcp), allocatable ::  d_obs(:)     ! Observed shear cross-correlations, data arranged like a vector
        real(mcp), allocatable ::  invcov(:,:)  ! Inverse Covariance Matrix
        logical :: ah_rescale = .false.         ! Rescale Inverse Covmat eg Anderson-Hartlap
        real(mcp) :: ah_realizations            ! Number of realizations for AH factor
        
        
        logical :: use_IA = .false.     !Use Intrinsic alignment calibration

        Type(TCubicSpline), allocatable :: Cl(:) !Theoretical Lensing cross power spectra

        ! Theoretical 2-point correlation functions
        ! First index is theta bin;
        ! Second index is the redshift cross-power bin, for example with 3 redshift bins
        ! 1-> 1,1; 2->1,2; 3->1,3; 4->2,2; 5 ->2,3; 6->3,3
        real(mcp), allocatable :: cplus(:,:), cminus(:,:)
    contains
    procedure :: Initialize
    procedure :: Cl_to_Ctheta
    procedure :: ArrangeForLike
    procedure :: GetBessels
    procedure :: GetLensingCls
    end type TCosmologyWLLikelihood
    
    public TCosmologyWLLikelihood

    contains

    subroutine Initialize(this)
    class(TCosmologyWLLikelihood) this
    integer dim, il
    
    this%needs_powerspectra = .true.
    this%needs_Weylpower = .true.
    this%needs_nonlinear_pk = .true.
    
    dim=2*this%nthetamax*this%ncl
    !Allocate static Arrays
    allocate(this%ll(this%nlmax))
    allocate(this%z(this%nzmax),this%etaz_f(this%nzmax,this%nbinmax),this%eta_norm_f(this%nbinmax))
    allocate(this%etaz_b(this%nzmax,this%nbinmax),this%eta_norm_b(this%nbinmax))
    allocate(this%theta(this%nthetamax))
    allocate(this%d_theo(dim),this%d_obs(dim),this%invcov(dim,dim))
    !Allocate Correlation function arrays (we zero them out for each likelihood evaluation
    allocate(this%cplus(this%nthetamax,this%ncl),this%cminus(this%nthetamax,this%ncl))
    allocate(this%Cl(this%ncl))

    !Initialize l values
    do il=1, this%nlmax
        this%ll(il)=exp(this%dlnl*(il-1))
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Initialize Bessel interpolation data !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call this%GetBessels()

    end subroutine Initialize

    subroutine ArrangeForLike(this)
    class(TCosmologyWLLikelihood) this

    call MpiStop('ERROR: You have not set up ArrangeForLike for '//trim(this%name))

    end subroutine ArrangeForLike

    !Subroutines below adapted from Angular Module of Julien Lesgourgues
    subroutine GetLensingCls(this,CMB,Theory,IAParam)
    Class(TCosmologyWLLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp), optional ::  IAParam
    !for etar(z)
    real(mcp), allocatable :: r(:),etar_f(:,:),g_f(:,:)
    real(mcp), allocatable :: etar_b(:,:),g_b(:,:)
    !Temporary holders of power spectra
    real(mcp) P_GG, P_GI, P_II
    real(mcp) kvalue, zvalue, rvalue
    real(mcp), allocatable :: Cl_integrand(:,:),Cl_temp(:,:)
    !loop indices
    integer nbin,nl,nc,zix,zix2,clc1,clc2
    real WL_IA

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! Allocate local arrays !!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(r(this%nzmax),etar_f(this%nzmax,this%nbinmax),g_f(this%nzmax,this%nbinmax))
    allocate(etar_b(this%nzmax,this%nbinmax),g_b(this%nzmax,this%nbinmax))
    allocate(Cl_integrand(this%nzmax,this%ncl),Cl_temp(this%ncl,this%nlmax))


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Various steps for C_l computation start here !!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!! JD: WL_IA for intrinsic alignment correction !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WL_IA = 0.
    if(this%use_IA .and. present(IAParam)) WL_IA=IAParam

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!! Convert eta(z) to eta(r) !!!!!!!!!!!!!!!!!!!!
    !!! (first unnormalized, second normalized to one) !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! for each value of z in the data's selection function file, compute r and dz/dr with splint
    ! compute also the selection function eta(r) = eta(z) dz/dr normalized to one
    do zix=1,this%nzmax
        r(zix)=Theory%R%value(this%z(zix))
        etar_f(zix,:)=this%etaz_f(zix,:)*Theory%dzdr%value(this%z(zix))/this%eta_norm_f(:)
        etar_b(zix,:)=this%etaz_b(zix,:)*Theory%dzdr%value(this%z(zix))/this%eta_norm_b(:)
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! For each bin, compute window function !!!!!!!!!!!
    !!!!! g(r) = sum_0^r_s_max dr_s eta(r_s) f_k(r_s-r)/f_k(r_s) !!!!!
    !!!!! !!!!!!!!!(simple trapezoidal integration is enough) !!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    g_f=0.
    g_b=0.
    do zix=2,this%nzmax-1
        do zix2=zix+1,this%nzmax
            g_f(zix,:)=g_f(zix,:)+0.5*(etar_f(zix2,:)*this%Calculator%f_k(r(zix2)-r(zix))/this%Calculator%f_k(r(zix2)) &
                +etar_f(zix2-1,:)*this%Calculator%f_k(r(zix2-1)-r(zix))/this%Calculator%f_k(r(zix2-1)))*(r(zix2)-r(zix2-1))

            g_b(zix,:)=g_b(zix,:)+0.5*(etar_b(zix2,:)*this%Calculator%f_k(r(zix2)-r(zix))/this%Calculator%f_k(r(zix2)) &
                +etar_b(zix2-1,:)*this%Calculator%f_k(r(zix2-1)-r(zix))/this%Calculator%f_k(r(zix2-1)))*(r(zix2)-r(zix2-1))
        end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! start loop over l  for computation of C_l^shear !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Cl_temp = 0.
    do nl=1,this%nlmax
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!! Set up integrand !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! Cl_integrand = P^{k l}_\kappa_integrand = (g_{k}(r)*g_{l}(r))*P(l/r,z(r)) !!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Cl_integrand=0.
        nc=0
        do zix=2,this%nzmax
            !Get k value and see if it is in bounds, if not integrand is zero
            rvalue = this%Calculator%f_k(r(zix))
            kvalue=this%ll(nl)/rvalue/CMB%h! k=l/r in units h/Mpc
            if(.not.(kvalue<exp(Theory%P_GG%x(1)) .or. kvalue>exp(Theory%P_GG%x(Theory%P_GG%nx)))) then
                zvalue = this%z(zix)
                P_GG = Theory%P_GG%PowerAt(kvalue,zvalue)/(CMB%h)**3 ! in units Mpc**3
                if(this%use_IA) then
                    P_GI = Theory%P_GI%PowerAt(kvalue,zvalue)/(CMB%h)**3 ! in units Mpc**3
                    P_II = Theory%P_II%PowerAt(kvalue,zvalue)/(CMB%h)**3 ! in units Mpc**3
                else
                    P_GI = 0
                    P_II = 0
                end if
                nc=0
                do clc1=1,this%nbinmax
                    do clc2=clc1, this%nbinmax
                        nc = nc+1
                        Cl_integrand(zix,nc)=g_f(zix,clc1)*g_b(zix,clc2)*P_GG &
                        -WL_IA*(g_f(zix,clc1)*etar_b(zix,clc2)+g_b(zix,clc2)&
                        *etar_f(zix,clc1))/rvalue*P_GI &
                        +WL_IA**2*etar_f(zix,clc1)*etar_b(zix,clc2)*P_II/rvalue**2
                    end do
                end do
                if(nc .ne. this%ncl) then
                    write(*,'("nc count is ",I0," after loop instead of ",I0)')nc, this%ncl
                    call MPIStop()
                end if
            end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! Integrate over r to get C_l^shear !!!!!!!!!
            !!! C_l^shear = \sum_0^rmax dr Cl_integrand !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Cl_temp(:,nl)=Cl_temp(:,nl)+0.5*(Cl_integrand(zix,:)&
            +Cl_integrand(zix-1,:))*(r(zix)-r(zix-1))
        end do
    end do

    if (Feedback >2) write(*,'("lmax = ",E15.8)')this%ll(this%nlmax)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Compute correlation function C_plus,minus(theta) !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do nc=1, this%ncl
        call this%Cl(nc)%Init(this%ll,CL_temp(nc,:),this%nlmax)
    end do

    end subroutine GetLensingCls

    subroutine Cl_to_Ctheta(this)
    class(TCosmologyWLLikelihood) this
    ! Step of integration such that J(x) cannot
    ! change to more than J(x+dxmax)
    ! 0.2 sufficient, no need to decrease
    real(mcp), parameter :: dxmax=0.2
    integer ntheta
    ! for integration
    integer nl, nc
    real(mcp) lll,l_previous,lmax,x
    real(mcp), allocatable :: integrand1(:),integrand1_previous(:)
    real(mcp), allocatable :: integrand2(:),integrand2_previous(:)
    integer l_increment ! 0 for logarithmic, 1 for linear

    allocate(integrand1(this%ncl))
    allocate(integrand2(this%ncl))
    allocate(integrand1_previous(this%ncl))
    allocate(integrand2_previous(this%ncl))

    ! loop over theta
    do ntheta=1, this%nthetamax

        ! maximum l value in the integral is the smallest of:
        ! * the max l value at which C_l have been computed
        ! * or the value such that l*theta=xstop
        lmax=min(this%ll(this%nlmax),this%xstop/this%theta(ntheta))

        ! initialize the sum
        this%cplus(ntheta,:)=0.
        this%cminus(ntheta,:)=0.

        ! Limit l-->0 : the integrand vanishes (due to factor l)
        l_previous=0.
        integrand1_previous(:)=0.
        integrand2_previous(:)=0.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! First non-zero l value in the integral !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nl=1
        lll=this%ll(nl)
        l_increment=0 ! start trying logarithmic increment

        do while (lll <= lmax)

            x=lll*this%theta(ntheta)

            ! find integrand
            do nc=1,this%ncl
                integrand1(nc)=this%Cl(nc)%value(lll)*this%Bess0%value(x)*lll/twopi
                integrand2(nc)=this%Cl(nc)%value(lll)*this%Bess4%value(x)*lll/twopi
            end do
            ! increment the sum
            this%cplus(ntheta,:)=this%cplus(ntheta,:)&
            +0.5*(integrand1_previous(:)+integrand1(:))*(lll-l_previous)
            this%cminus(ntheta,:)=this%cminus(ntheta,:)&
            +0.5*(integrand2_previous(:)+integrand2(:))*(lll-l_previous)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Next l value in the integral !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            l_previous=lll
            if (l_increment == 0) then
                if (nl+1 <= this%nlmax) then ! test whether integral not already finished
                    if ((this%ll(nl+1)-this%ll(nl))*this%theta(ntheta) <= dxmax) then
                        nl=nl+1
                        lll=this%ll(nl)   ! continue with logarithmic increment
                    else
                        l_increment=0 !switch to linear increment
                        lll=l_previous+dxmax/this%theta(ntheta) ! linear increment
                    end if
                else
                    lll=lmax+1. ! integral already finished. This will stop the loops.
                end if
            else
                lll=l_previous+dxmax/this%theta(ntheta) ! linear increment
            end if
            integrand1_previous(:)=integrand1(:)
            integrand2_previous(:)=integrand2(:)
        end do
    end do

    end subroutine Cl_to_Ctheta

    subroutine GetBessels(this)
    class(TCosmologyWLLikelihood) this
    integer i
    real(mcp), allocatable :: x(:),Bess0(:),Bess4(:)
    integer max_bes_ix

    max_bes_ix = nint(this%xstop/dbessel)+3
    allocate(x(max_bes_ix),Bess0(max_bes_ix),Bess4(max_bes_ix))
    Bess0(1)=1
    Bess4(1)=0
    x(1)=0
    do i=2, max_bes_ix
        x(i) = (i-1)*dbessel
        Bess0(i) = Bessel_J0(x(i))
        Bess4(i) = Bessel_JN(4,x(i))
    end do

    call this%Bess0%Init(x,Bess0,max_bes_ix)
    call this%Bess4%Init(x,Bess4,max_bes_ix)

    end subroutine GetBessels

    end module WeakLen_Common
