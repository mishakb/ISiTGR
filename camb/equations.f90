    ! Equations module for dark energy with constant equation of state parameter w
    ! allowing for perturbations based on a quintessence model
    ! by Antony Lewis (http://cosmologist.info/)

    ! Dec 2003, fixed (fatal) bug in tensor neutrino setup
    ! Changes to tight coupling approximation
    ! June 2004, fixed problem with large scale polarized tensors; support for vector modes
    ! Generate vector modes on their own. The power spectrum is taken from the scalar parameters.
    ! August 2004, fixed reionization term in lensing potential
    ! Nov 2004, change massive neutrino l_max to be consistent with massless if light
    ! Apr 2005, added DoLateRadTruncation option
    ! June 2006, added support for arbitary neutrino mass splittings
    ! Nov 2006, tweak to high_precision transfer function accuracy at lowish k
    ! June 2011, improved radiation approximations from arXiv: 1104.2933; Some 2nd order tight coupling terms
    !            merged fderivs and derivs so flat and non-flat use same equations; more precomputed arrays
    !            optimized neutrino sampling, and reorganised neutrino integration functions
    ! Feb 2013: fixed various issues with accuracy at larger neutrino masses
    ! Mar 2014: fixes for tensors with massive neutrinos

    module LambdaGeneral
    use precision
    implicit none

    real(dl)  :: w_lam = -1_dl !p/rho for the dark energy (assumed constant)
    real(dl) :: cs2_lam = 1_dl
    !comoving sound speed. Always exactly 1 for quintessence
    !(otherwise assumed constant, though this is almost certainly unrealistic)

!    real(dl), parameter :: wa_ppf = 0._dl !Not used here, just for compatibility with e.g. halofit
    real(dl):: wa_ppf = 0._dl

    logical :: w_perturb = .true.
    !If you are tempted to set this = .false. read
    ! http://cosmocoffee.info/viewtopic.php?t=811
    ! http://cosmocoffee.info/viewtopic.php?t=512

    contains

    subroutine DarkEnergy_ReadParams(Ini)
    use IniFile
    Type(TIniFile) :: Ini

    w_lam = Ini_Read_Double_File(Ini,'w', -1.d0)
    cs2_lam = Ini_Read_Double_File(Ini,'cs2_lam',1.d0)

    end subroutine DarkEnergy_ReadParams

    end module LambdaGeneral

	!>ISiTGR MOD START: Module for ISiTGR, different functions based on MG parameters -----------------------------
    module ISiTGR
    use precision
    implicit none
		
    Type ISiTGR_CP
        real(dl) :: Q0,Qinf	!models 3-6
        real(dl) :: DR0,DRinf !models 3-6
        real(dl) :: s !models 3, 5
        real(dl) :: k_c ! The k value where we change k-bins for true binning !models 1-6
        real(dl) :: Q1,Q2,Q3,Q4 !binning models 1 and 2
        real(dl) :: D1,D2,D3,D4 !binning models 1 and 2
        real(dl) :: z_div !This is the size of each z-bin !models 1 and 2
        real(dl) :: z_TGR !The redshift below which we test GR usueally 2*z_div !models 1 and 2
        real(dl) :: z_tw  !transition width between redshift bins !models 1 and 2
        real(dl) :: k_tw  !transition width between k-bins. !models 1 and 2
        logical :: t_dep 
        logical :: R_func  
        logical  :: true_bin  
        logical  :: ISiTGR_BIN
        !CGQ ------------------------------------------------
		!Functional Form Parameters
		real(dl) :: E11, E22 !For mu, eta Planck's parameterization !functional form
		real(dl) :: mu0, Sigma0 !For mu, Sigma DES parameterization !functional form
		!For scale-dependence
		real(dl) :: c1, c2, lambda
		!Binning Form Parameters
		real(dl) :: mu1, mu2, mu3, mu4 !binning models 
        real(dl) :: eta1, eta2, eta3, eta4 !binning models 
        real(dl) :: Sigma1, Sigma2, Sigma3, Sigma4 !binning models 
		!Dark Energy parameterizations
		real(dl) :: w0, wa, wp, a_p!CGQ for 3 new models of DE
		logical :: ISiTGR_mueta !for functional form
		logical :: ISiTGR_muSigma !for functional form
		logical :: ISiTGR_QDR !for functional form
		logical :: ISiTGR_BIN_mueta !for binning method
		logical :: ISiTGR_BIN_muSigma !for binning method
		!logical :: ISiTGR_k_dep !for scale dependence
		integer :: GR = 0!to use default GR without MG
		integer :: DE_eqstate
        !CGQ ------------------------------------------------
    end type ISiTGR_CP
    Type(ISiTGR_CP) TGR

    contains

    !################### Modified Gravity Function of ISiTGR Parameters #####################
	
	!CGQ added extra functions to use mu-eta and mu-Sigma parameterization
	
    function ISiTGR_mu(k,a,adotoa)
	use ModelParams
	use constants
    real(dl), intent(in) :: k, a, adotoa
    real(dl) ISiTGR_mu
	real(dl) :: s1_k, s2_k
	real(dl) :: mu_MG

	!binning method expression
	if((TGR%ISiTGR_BIN_mueta) .or. (TGR%ISiTGR_BIN_muSigma)) then
        ISiTGR_mu = (1+ISiTGR_mu_Z1(k) +(ISiTGR_mu_Z2(k)-ISiTGR_mu_Z1(k))*tanh((1.d0/a-1.d0-TGR%z_div)/TGR%z_tw) &
        +(1-ISiTGR_mu_Z2(k))*tanh((1.d0/a-1.d0-TGR%z_TGR)/TGR%z_tw))/2.d0

	!functional form expression
	else if (TGR%ISiTGR_mueta) then
	
		!adding functions for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0 !s2_k needed to get s1_k
		s1_k = (1.d0+TGR%c1*s2_k)/(1.d0+s2_k) !function s1_k puts scale dependence into the mu parameter
		!computing E11 as in Planck 2015 parameterization
		mu_MG = 1.d0 + TGR%E11 * OmegaDE(a,adotoa) * s1_k
		
		ISiTGR_mu = mu_MG
	
	!functional form expression
	else if (TGR%ISiTGR_muSigma) then
	
		!adding an extra factor for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0
		s1_k = (1.d0+TGR%c1*s2_k)/(1.d0+s2_k)
		!computed mu as in DES 2018 paper adding scale dependence 
		mu_MG = 1.d0 + TGR%mu0 * (OmegaDE(a,adotoa)/CP%omegav) * s1_k !when TGR%Lambda=0 original DES parameterization is recovered
					
		ISiTGR_mu = mu_MG
		
	end if
	
	end function ISiTGR_mu

	function ISiTGR_mu_dot(k,a,adotoa,Hdot)
	use ModelParams
	use constants
    real(dl), intent(in) :: k, a, adotoa, Hdot
    real(dl) ISiTGR_mu_dot
	real(dl) :: mudot_MG
	real(dl) :: s1_k, s2_k, s1_k_dot, s2_k_dot 
	
	!binning method expression
   	if((TGR%ISiTGR_BIN_mueta) .or. (TGR%ISiTGR_BIN_muSigma)) then
        ISiTGR_mu_dot = adotoa/2.d0/TGR%z_tw/a*((ISiTGR_mu_Z1(k)-ISiTGR_mu_Z2(k))/cosh((1.d0/a-1.d0-TGR%z_div)/TGR%z_tw)**2.d0 &
        +(ISiTGR_mu_Z2(k)-1)/cosh((1.d0/a-1.d0-TGR%z_TGR)/TGR%z_tw)**2.d0)

	!functional form expression
	else if (TGR%ISiTGR_mueta) then
	
		!adding an extra factor for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0
		s1_k = (1.d0+TGR%c1*s2_k)/(1.d0+s2_k)
		s2_k_dot = 2.d0*s2_k*(Hdot-adotoa**2.d0)/adotoa
		s1_k_dot = s2_k_dot*(TGR%c1-1.d0)/(1.d0+s2_k)**2.d0
		mudot_MG = TGR%E11*(OmegaDEdot(a,adotoa,Hdot)*s1_k + OmegaDE(a,adotoa)*s1_k_dot)
		
		ISiTGR_mu_dot = mudot_MG
						
	!functional form expression
	else if (TGR%ISiTGR_muSigma) then
		
		!adding an extra factor for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0
		s1_k = (1.d0+TGR%c1*s2_k)/(1.d0+s2_k)
		s2_k_dot = 2.d0*s2_k*(Hdot-adotoa**2.d0)/adotoa
		s1_k_dot = s2_k_dot*(TGR%c1-1.d0)/(1+s2_k)**2.d0
		mudot_MG = TGR%mu0*(OmegaDEdot(a,adotoa,Hdot)/CP%omegav*s1_k + OmegaDE(a,adotoa)/CP%omegav*s1_k_dot)

		ISiTGR_mu_dot = mudot_MG

	end if

	end function ISiTGR_mu_dot
	!---------------------------------------------------------------------
	
	function ISiTGR_eta(k,a,adotoa)
	use ModelParams
	use constants
   	real(dl), intent(in) :: k, a,adotoa
   	real(dl) ISiTGR_eta
	real(dl) :: s1_k, s2_k
	real(dl) :: eta_MG
	
	!binning method expression
	if(TGR%ISiTGR_BIN_mueta) then
        ISiTGR_eta = (1+ISiTGR_eta_Z1(k) +(ISiTGR_eta_Z2(k)-ISiTGR_eta_Z1(k))*tanh((1.d0/a-1.d0-TGR%z_div)/TGR%z_tw) &
        +(1-ISiTGR_eta_Z2(k))*tanh((1.d0/a-1.d0-TGR%z_TGR)/TGR%z_tw))/2.d0
	
	!functional form expression
	else
		
		!adding an extra factor for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0
		s1_k = (1.d0+TGR%c2*s2_k)/(1.d0+s2_k)
		!scale dependent
		eta_MG = 1.d0 + TGR%E22 * OmegaDE(a,adotoa) * s1_k
			
		ISiTGR_eta = eta_MG			
			
	end if
	
	end function ISiTGR_eta
	
	function ISiTGR_eta_dot(k,a,adotoa,Hdot)
	use ModelParams
	use constants
    real(dl), intent(in) :: k, a,adotoa, Hdot
    real(dl) ISiTGR_eta_dot
	real(dl) :: etadot_MG
	real(dl) :: s1_k, s2_k, s1_k_dot, s2_k_dot
	
	!binning method expression
    if(TGR%ISiTGR_BIN_mueta) then
        ISiTGR_eta_dot = adotoa/2.d0/TGR%z_tw/a*((ISiTGR_eta_Z1(k)-ISiTGR_eta_Z2(k))/cosh((1.d0/a-1.d0-TGR%z_div)/TGR%z_tw)**2.d0 &
        +(ISiTGR_eta_Z2(k)-1)/cosh((1.d0/a-1.d0-TGR%z_TGR)/TGR%z_tw)**2.d0)
	
	else

		!adding an extra factor for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0
		s1_k = (1.d0+TGR%c2*s2_k)/(1.d0+s2_k)
		s2_k_dot = 2.d0*s2_k*(Hdot-adotoa**2.d0)/adotoa
		s1_k_dot = s2_k_dot*(TGR%c2-1.d0)/(1+s2_k)**2.d0
		!scale dependent
		etadot_MG = TGR%E22*(OmegaDEdot(a,adotoa,Hdot)*s1_k + OmegaDE(a,adotoa)*s1_k_dot)

		ISiTGR_eta_dot = etadot_MG	
		
	end if
	
	end function ISiTGR_eta_dot
	!---------------------------------------------------------------------
	
	function ISiTGR_Sigma(k,a,adotoa)
	use ModelParams
   	use constants
    real(dl), intent(in) :: k, a, adotoa
    real(dl) ISiTGR_Sigma
	real(dl) :: Sigma_MG
	real(dl) :: s1_k, s2_k
	
	!binning method expression
	if(TGR%ISiTGR_BIN_muSigma) then
        ISiTGR_Sigma = (1+ISiTGR_Sigma_Z1(k) +(ISiTGR_Sigma_Z2(k)-ISiTGR_Sigma_Z1(k))*tanh((1.d0/a-1.d0-TGR%z_div)/TGR%z_tw) &
        +(1-ISiTGR_Sigma_Z2(k))*tanh((1.d0/a-1.d0-TGR%z_TGR)/TGR%z_tw))/2.d0
	
	!functional form expression
	else
		
		!adding an extra factor for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0
		s1_k = (1.d0+TGR%c2*s2_k)/(1.d0+s2_k)
		!scale dependent
		Sigma_MG = 1.d0 + TGR%Sigma0 * OmegaDE(a,adotoa)/CP%omegav * s1_k
	
		ISiTGR_Sigma = Sigma_MG
	
	end if
	
	end function ISiTGR_Sigma
	
	function ISiTGR_Sigma_dot(k,a,adotoa,Hdot)
	use ModelParams
    use constants
    real(dl), intent(in) :: k, a, adotoa, Hdot
    real(dl) ISiTGR_Sigma_dot
	real(dl) :: Sigmadot_MG
	real(dl) :: s1_k, s2_k, s1_k_dot, s2_k_dot
	
	!binning method expression
    if(TGR%ISiTGR_BIN_muSigma) then
        ISiTGR_Sigma_dot = adotoa/2.d0/TGR%z_tw/a*((ISiTGR_Sigma_Z1(k)-ISiTGR_Sigma_Z2(k))/ & 
		cosh((1.d0/a-1.d0-TGR%z_div)/TGR%z_tw)**2.d0+(ISiTGR_Sigma_Z2(k)-1)/cosh((1.d0/a-1.d0-TGR%z_TGR)/TGR%z_tw)**2.d0)
	
	!functional form expression
	else
		
		!adding an extra factor for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0
		s1_k = (1.d0+TGR%c2*s2_k)/(1.d0+s2_k)
		s2_k_dot = 2.d0*s2_k*(Hdot-adotoa**2.d0)/adotoa
		s1_k_dot = s2_k_dot*(TGR%c2-1.d0)/(1+s2_k)**2.d0
		!scale dependent
		Sigmadot_MG = TGR%Sigma0*(OmegaDEdot(a,adotoa,Hdot)/CP%omegav*s1_k + OmegaDE(a,adotoa)/CP%omegav*s1_k_dot)
			
		ISiTGR_Sigma_dot= Sigmadot_MG	

	end if
	
	end function ISiTGR_Sigma_dot
	!---------------------------------------------------------------------
    !End of functions written by CGQ

    function ISiTGR_Q(k,a,adotoa)
	use ModelParams
	use constants
    real(dl), intent(in) :: k, a, adotoa
    real(dl) ISiTGR_Q
	real(dl) :: s1_k, s2_k

    if(TGR%ISiTGR_BIN) then
        ISiTGR_Q = (1+ISiTGR_Q_Z1(k) +(ISiTGR_Q_Z2(k)-ISiTGR_Q_Z1(k))*tanh((1.d0/a-1.d0-TGR%z_div)/TGR%z_tw) &
        +(1-ISiTGR_Q_Z2(k))*tanh((1.d0/a-1.d0-TGR%z_TGR)/TGR%z_tw))/2.d0
		
	else if (TGR%ISiTGR_QDR .eqv. .true.) then

		!adding an extra factor for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0
		s1_k = (1.d0+TGR%c1*s2_k)/(1.d0+s2_k)
		!scale dependent
		ISiTGR_Q = 1.d0 + TGR%Q0 * OmegaDE(a,adotoa)/CP%omegav * s1_k
		
    else if(TGR%t_dep) then
        ISiTGR_Q=(TGR%Q0*dexp(-k/TGR%k_c) + TGR%Qinf*(1.d0-dexp(-k/TGR%k_c))-1.d0)*a**TGR%s+1.d0
		
	else 
	    ISiTGR_Q=TGR%Q0*dexp(-k/TGR%k_c) + TGR%Qinf*(1.d0-dexp(-k/TGR%k_c))
		
	end if

    end function ISiTGR_Q

    function ISiTGR_Q_dot(k,a,adotoa,Hdot)
	use ModelParams
    use constants
    real(dl), intent(in) :: k, a, adotoa, Hdot
    real(dl) ISiTGR_Q_dot
	real(dl) :: s1_k, s2_k, s1_k_dot, s2_k_dot
	
    if(TGR%ISiTGR_BIN) then
        ISiTGR_Q_dot = adotoa/2.d0/TGR%z_tw/a*((ISiTGR_Q_Z1(k)-ISiTGR_Q_Z2(k))/cosh((1.d0/a-1.d0-TGR%z_div)/TGR%z_tw)**2.d0 &
        +(ISiTGR_Q_Z2(k)-1)/cosh((1.d0/a-1.d0-TGR%z_TGR)/TGR%z_tw)**2.d0)
		
	!functional form expression
	else if (TGR%ISiTGR_QDR .eqv. .true.) then
	
		!adding an extra factor for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0
		s1_k = (1.d0+TGR%c1*s2_k)/(1.d0+s2_k)
		s2_k_dot = 2.d0*s2_k*(Hdot-adotoa**2.d0)/adotoa
		s1_k_dot = s2_k_dot*(TGR%c1-1.d0)/(1+s2_k)**2.d0
		!scale dependent
		ISiTGR_Q_dot = TGR%Q0*(OmegaDEdot(a,adotoa,Hdot)/CP%omegav*s1_k + OmegaDE(a,adotoa)/CP%omegav*s1_k_dot)
    
	else if(TGR%t_dep) then
        ISiTGR_Q_dot=(ISiTGR_Q(k,a,adotoa)-1.d0)*TGR%s*adotoa
	
	else
		ISiTGR_Q_dot=0.
		
	end if

    end function ISiTGR_Q_dot

    function ISiTGR_D(k,a,adotoa)
	use ModelParams
	use constants
    real(dl), intent(in) :: k, a, adotoa
    real(dl) ISiTGR_D
	real(dl) :: s1_k, s2_k

    if(TGR%ISiTGR_BIN) then
        ISiTGR_D = (1+ISiTGR_D_Z1(k) +(ISiTGR_D_Z2(k)-ISiTGR_D_Z1(k))*tanh((1.d0/a-1.d0-TGR%z_div)/TGR%z_tw) &
        +(1-ISiTGR_D_Z2(k))*tanh((1.d0/a-1.d0-TGR%z_TGR)/TGR%z_tw))/2.d0
	
	else
		
		if(TGR%R_func)then
    	    ISiTGR_D = ISiTGR_Q(k,a,adotoa)*(1.d0+ISiTGR_DR(k,a,adotoa))/2.d0
		
		else
			ISiTGR_D = ISiTGR_DR(k,a,adotoa)
		
		end if	

    end if
	
    end function ISiTGR_D

    function ISiTGR_D_dot(k,a,adotoa,Hdot)
	use ModelParams
	use constants
    real(dl), intent(in) :: k, a, adotoa, Hdot
    real(dl) ISiTGR_D_dot
	real(dl) :: s1_k, s2_k, s1_k_dot, s2_k_dot

    if(TGR%ISiTGR_BIN) then
        ISiTGR_D_dot = adotoa/2.d0/TGR%z_tw/a*((ISiTGR_D_Z1(k)-ISiTGR_D_Z2(k))/cosh((1.d0/a-1.d0-TGR%z_div)/TGR%z_tw)**2.d0 &
        +(ISiTGR_D_Z2(k)-1)/cosh((1.d0/a-1.d0-TGR%z_TGR)/TGR%z_tw)**2.d0)
    
	else 
	
		if(TGR%R_func)then
        	ISiTGR_D_dot = ISiTGR_Q_dot(k,a,adotoa,Hdot)*(1.d0+ISiTGR_DR(k,a,adotoa))/2.d0 &
        	+ISiTGR_Q(k,a,adotoa)*ISiTGR_DR_dot(k,a,adotoa,Hdot)/2.d0

		else
			ISiTGR_D_dot = ISiTGR_DR_dot(k,a,adotoa,Hdot)
			
		end if
    
	end if
	
    end function ISiTGR_D_dot

    !################### Auxilary functions for Binning methods of ISiTGR #####################

    function ISiTGR_DR(k,a,adotoa)
	use ModelParams
	use constants
    real(dl), intent(in) :: k, a, adotoa
    real(dl) ISiTGR_DR
	real(dl) :: s1_k, s2_k

	if (TGR%ISiTGR_QDR .eqv. .true.) then
	
		!adding an extra factor for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0
		s1_k = (1.d0+TGR%c2*s2_k)/(1.d0+s2_k)
		!scale dependent
		ISiTGR_DR = 1.d0 + TGR%DR0 * OmegaDE(a,adotoa)/CP%omegav * s1_k

	else

	    if(TGR%t_dep) then
    	    ISiTGR_DR=(TGR%dr0*dexp(-k/TGR%k_c) + TGR%drinf*(1.d0-dexp(-k/TGR%k_c))-1.d0)*a**TGR%s+1.d0
	    else
	        ISiTGR_DR=TGR%dr0*dexp(-k/TGR%k_c) + TGR%drinf*(1.d0-dexp(-k/TGR%k_c))
	    end if
	
	end if

    end function ISiTGR_DR

    function ISiTGR_DR_dot(k,a,adotoa,Hdot)
	use ModelParams
	use constants
    real(dl), intent(in) :: k, a, adotoa, Hdot
    real(dl) ISiTGR_DR_dot
	real(dl) :: s1_k, s2_k, s1_k_dot, s2_k_dot
	
	!functional form expression
	if (TGR%ISiTGR_QDR .eqv. .true.) then
			
		!adding an extra factor for scale dependence
		s2_k = (TGR%lambda*(adotoa/a)/k)**2.d0
		s1_k = (1.d0+TGR%c2*s2_k)/(1.d0+s2_k)
		s2_k_dot = 2.d0*s2_k*(Hdot-adotoa**2.d0)/adotoa
		s1_k_dot = s2_k_dot*(TGR%c2-1.d0)/(1+s2_k)**2.d0
		!scale dependent
		ISiTGR_DR_dot = TGR%DR0*(OmegaDEdot(a,adotoa,Hdot)/CP%omegav*s1_k + OmegaDE(a,adotoa)/CP%omegav*s1_k_dot)
		
	else

	    if(TGR%t_dep) then
    	    ISiTGR_DR_dot=(ISiTGR_DR(k,a,adotoa)-1.d0)*TGR%s*adotoa
	    else
	        ISiTGR_DR_dot=0.
	    end if

	end if

    end function ISiTGR_DR_dot

    !Auxilary functions for ISiTGR - Binning forms
    function ISiTGR_Q_Z1(k)
    real(dl), intent(in) :: k
    real(dl) ISiTGR_Q_Z1

    if(TGR%true_bin) then
        ISiTGR_Q_Z1=(TGR%Q2+TGR%Q1)/2.d0 + (TGR%Q2-TGR%Q1)/2.d0*tanh((k-TGR%k_c)/TGR%k_tw)
    else
        ISiTGR_Q_Z1=TGR%Q1*dexp(-k/TGR%k_c) + TGR%Q2*(1.d0-dexp(-k/TGR%k_c))
    end if

    end function ISiTGR_Q_Z1

    function ISiTGR_Q_Z2(k)
    real(dl), intent(in) :: k
    real(dl) ISiTGR_Q_Z2

    if(TGR%true_bin) then
        ISiTGR_Q_Z2=(TGR%Q4+TGR%Q3)/2.d0 + (TGR%Q4-TGR%Q3)/2.d0*tanh((k-TGR%k_c)/TGR%k_tw)
    else
        ISiTGR_Q_Z2=TGR%Q3*dexp(-k/TGR%k_c) + TGR%Q4*(1.d0-dexp(-k/TGR%k_c))
    end if

    end function ISiTGR_Q_Z2

    function ISiTGR_D_Z1(k)
    real(dl), intent(in) :: k
    real(dl) ISiTGR_D_Z1

    if(TGR%true_bin) then
        ISiTGR_D_Z1=(TGR%D2+TGR%D1)/2.d0 + (TGR%D2-TGR%D1)/2.d0*tanh((k-TGR%k_c)/TGR%k_tw)
    else
        ISiTGR_D_Z1=TGR%D1*dexp(-k/TGR%k_c) + TGR%D2*(1.d0-dexp(-k/TGR%k_c))
    end if

    end function ISiTGR_D_Z1

    function ISiTGR_D_Z2(k)
    real(dl), intent(in) :: k
    real(dl) ISiTGR_D_Z2

    if(TGR%true_bin) then
        ISiTGR_D_Z2=(TGR%D4+TGR%D3)/2.d0 + (TGR%D4-TGR%D3)/2.d0*tanh((k-TGR%k_c)/TGR%k_tw)
    else
        ISiTGR_D_Z2=TGR%D3*dexp(-k/TGR%k_c) + TGR%D4*(1.d0-dexp(-k/TGR%k_c))
    end if

    end function ISiTGR_D_Z2

    !Auxilary functions for ISiTGR - Binning forms !CGQ analogous expressions for mueta and muSigma binning forms
	!---------------------------------------------------------------------
    function ISiTGR_mu_Z1(k)
    real(dl), intent(in) :: k
    real(dl) ISiTGR_mu_Z1

    if(TGR%true_bin) then
        ISiTGR_mu_Z1=(TGR%mu2+TGR%mu1)/2.d0 + (TGR%mu2-TGR%mu1)/2.d0*tanh((k-TGR%k_c)/TGR%k_tw)
    else
        ISiTGR_mu_Z1=TGR%mu1*dexp(-k/TGR%k_c) + TGR%mu2*(1.d0-dexp(-k/TGR%k_c))
    end if

    end function ISiTGR_mu_Z1

    function ISiTGR_mu_Z2(k)
    real(dl), intent(in) :: k
    real(dl) ISiTGR_mu_Z2

    if(TGR%true_bin) then
        ISiTGR_mu_Z2=(TGR%mu4+TGR%mu3)/2.d0 + (TGR%mu4-TGR%mu3)/2.d0*tanh((k-TGR%k_c)/TGR%k_tw)
    else
        ISiTGR_mu_Z2=TGR%mu3*dexp(-k/TGR%k_c) + TGR%mu4*(1.d0-dexp(-k/TGR%k_c))
    end if

    end function ISiTGR_mu_Z2
	!---------------------------------------------------------------------
    function ISiTGR_eta_Z1(k)
    real(dl), intent(in) :: k
    real(dl) ISiTGR_eta_Z1

    if(TGR%true_bin) then
        ISiTGR_eta_Z1=(TGR%eta2+TGR%eta1)/2.d0 + (TGR%eta2-TGR%eta1)/2.d0*tanh((k-TGR%k_c)/TGR%k_tw)
    else
        ISiTGR_eta_Z1=TGR%eta1*dexp(-k/TGR%k_c) + TGR%eta2*(1.d0-dexp(-k/TGR%k_c))
    end if

    end function ISiTGR_eta_Z1

    function ISiTGR_eta_Z2(k)
    real(dl), intent(in) :: k
    real(dl) ISiTGR_eta_Z2

    if(TGR%true_bin) then
        ISiTGR_eta_Z2=(TGR%eta4+TGR%eta3)/2.d0 + (TGR%eta4-TGR%eta3)/2.d0*tanh((k-TGR%k_c)/TGR%k_tw)
    else
        ISiTGR_eta_Z2=TGR%eta3*dexp(-k/TGR%k_c) + TGR%eta4*(1.d0-dexp(-k/TGR%k_c))
    end if

    end function ISiTGR_eta_Z2
	!---------------------------------------------------------------------
    function ISiTGR_Sigma_Z1(k)
    real(dl), intent(in) :: k
    real(dl) ISiTGR_Sigma_Z1

    if(TGR%true_bin) then
        ISiTGR_Sigma_Z1=(TGR%Sigma2+TGR%Sigma1)/2.d0 + (TGR%Sigma2-TGR%Sigma1)/2.d0*tanh((k-TGR%k_c)/TGR%k_tw)
    else
        ISiTGR_Sigma_Z1=TGR%Sigma1*dexp(-k/TGR%k_c) + TGR%Sigma2*(1.d0-dexp(-k/TGR%k_c))
    end if

    end function ISiTGR_Sigma_Z1

    function ISiTGR_Sigma_Z2(k)
    real(dl), intent(in) :: k
    real(dl) ISiTGR_Sigma_Z2

    if(TGR%true_bin) then
        ISiTGR_Sigma_Z2=(TGR%Sigma4+TGR%Sigma3)/2.d0 + (TGR%Sigma4-TGR%Sigma3)/2.d0*tanh((k-TGR%k_c)/TGR%k_tw)
    else
        ISiTGR_Sigma_Z2=TGR%Sigma3*dexp(-k/TGR%k_c) + TGR%Sigma4*(1.d0-dexp(-k/TGR%k_c))
    end if

    end function ISiTGR_Sigma_Z2
	
	!################### Dark Energy parameterization of ISiTGR - modeling Omega_Dark_Energy #####################
	
	! CGQ patch for Dark Energy
	function OmegaDE(a, adotoa)
	use ModelParams
	use constants
    real(dl), intent(in) :: a, adotoa
    real(dl) :: OmegaDE

	if ( TGR%DE_eqstate == 0 ) then
		OmegaDE = CP%omegav * ((CP%H0*1000.d0/c)/(adotoa/a))**2.d0
		
	else if ( TGR%DE_eqstate == 1 ) then
		OmegaDE = CP%omegav * a**(-3.d0*(1.d0+TGR%w0)) * ((CP%H0*1000.d0/c)/(adotoa/a))**2.d0
	
	else if ( TGR%DE_eqstate == 2 ) then
		OmegaDE = CP%omegav * ((CP%H0*1000.d0/c)/(adotoa/a))**2.d0 * a**(-3.d0*(1.d0+TGR%w0+TGR%wa)) * exp(3.d0*TGR%wa*(a-1.d0)) 

	else if ( TGR%DE_eqstate == 3 ) then
		OmegaDE = CP%omegav * ((CP%H0*1000.d0/c)/(adotoa/a))**2.d0 * a**(-3.d0*(1.d0+TGR%wp+TGR%a_p*TGR%wa)) * exp(3.d0*TGR%wa*(a-1.d0)) 
		
	end if
	
    end function OmegaDE
	
	!derivative of Omega_Dark_Energy
	function OmegaDEdot(a, adotoa, Hdot)
	use ModelParams
    real(dl), intent(in) :: a, adotoa, Hdot
    real(dl) :: OmegaDEdot

	if ( TGR%DE_eqstate == 0 ) then
		OmegaDEdot = - OmegaDE(a,adotoa) * 2.d0 * (Hdot - adotoa**2.d0)/adotoa 
	
	else if ( TGR%DE_eqstate == 1 ) then
		OmegaDEdot = - OmegaDE(a,adotoa) * ( 2.d0 * (Hdot - adotoa**2.d0)/adotoa + &
		 			3.d0 * (1.d0+TGR%w0) * adotoa)

	else if ( TGR%DE_eqstate == 2 ) then
		OmegaDEdot = -OmegaDE(a,adotoa) * ( 2.d0 * (Hdot - adotoa**2.d0)/adotoa + &
					3.d0 * (1.d0+TGR%w0) * adotoa + 3.d0 * TGR%wa * adotoa * (1.d0-a) )
		
	else if ( TGR%DE_eqstate == 3 ) then
		OmegaDEdot = -OmegaDE(a,adotoa) * ( 2.d0 * (Hdot - adotoa**2.d0)/adotoa + &
					3.d0 * (1.d0+TGR%wp) * adotoa + 3.d0 * TGR%wa * adotoa * (TGR%a_p-a) )
		
	end if
	
    end function OmegaDEdot
	
!	function rhoDE(a) !rhov_t * a**2
!	use ModelParams
!	use constants
!    real(dl), intent(in) :: a
!    real(dl) rhoDE
!	
!	if ( TGR%DE_eqstate == 0 ) then
!		rhoDE = 3.d0*((CP%H0*1000.d0/c)**2)*CP%omegav*a**2.d0
!
!	else if ( TGR%DE_eqstate == 1 ) then
!		rhoDE = 3.d0*((CP%H0*1000.d0/c)**2)*CP%omegav*a**(-1.d0-3.d0*TGR%w0)
!
!	else if ( TGR%DE_eqstate == 2 ) then
!		rhoDE = 3.d0*((CP%H0*1000.d0/c)**2)*CP%omegav*a**(-1.d0-3.d0*(TGR%w0+TGR%wa))*exp(3.d0*TGR%wa*(a-1.d0))
!
!	else if ( TGR%DE_eqstate == 3 ) then
!		rhoDE = 3.d0*((CP%H0*1000.d0/c)**2)*CP%omegav*a**(-1.d0-3.d0*(TGR%wp+TGR%a_p*TGR%wa))*exp(3.d0*TGR%wa*(a-1.d0))
!
!	end if	
!	
!	end function rhoDE
!	
!	function presDE(a,DErho) 
!    real(dl), intent(in) :: a,DErho
!    real(dl) :: presDE
!	
!	if ( TGR%DE_eqstate == 0 ) then
!		presDE = -DErho
!
!	else if ( TGR%DE_eqstate == 1 ) then
!		presDE = TGR%w0*DErho
!
!	else if ( TGR%DE_eqstate == 2 ) then
!		presDE = (TGR%w0+(1.d0-a)*TGR%wa)*DErho
!		
!	else if ( TGR%DE_eqstate == 3 ) then
!		presDE = (TGR%wp+(TGR%a_p-a)*TGR%wa)*DErho
!
!	end if	
!	
!	end function presDE

	! CGQ End of patch for Dark Energy

    end module ISiTGR
	!<ISiTGR MOD END -----------------------------------------------------------------------------------------------------------------

    !Return OmegaK - modify this if you add extra fluid components
    function GetOmegak()
    use precision
    use ModelParams
    real(dl)  GetOmegak
    GetOmegak = 1 - (CP%omegab+CP%omegac+CP%omegav+CP%omegan)

    end function GetOmegak


    subroutine init_background
    !This is only called once per model, and is a good point to do any extra initialization.
    !It is called before first call to dtauda, but after
    !massive neutrinos are initialized and after GetOmegak
    end  subroutine init_background


    !Background evolution
    function dtauda(a)
    !get d tau / d a
	!>ISiTGR MOD START
	use ISiTGR
	!<ISiTGR MOD END
    use precision
    use ModelParams
    use MassiveNu
    use LambdaGeneral
    implicit none
    real(dl) dtauda
    real(dl), intent(IN) :: a
    real(dl) rhonu,grhoa2, a2
    integer nu_i

    a2=a**2

    !  8*pi*G*rho*a**4.
    grhoa2=grhok*a2+(grhoc+grhob)*a+grhog+grhornomass
	
	!>ISiTGR MOD START
	if (TGR%GR==0) then
	    if (w_lam == -1._dl) then
	        grhoa2=grhoa2+grhov*a2**2
	    else
	        grhoa2=grhoa2+grhov*a**(1-3*w_lam)
	    end if
	else
	
		if ( TGR%DE_eqstate == 0 ) then
		grhoa2 = grhoa2 + grhov*a**4.d0

		else if ( TGR%DE_eqstate == 1 ) then
		grhoa2 = grhoa2 + grhov*a**(1.d0-3.d0*TGR%w0)

		else if ( TGR%DE_eqstate == 2 ) then
		grhoa2 = grhoa2 + grhov*a**(1.d0-3.d0*(TGR%w0+TGR%wa))*exp(3.d0*TGR%wa*(a-1.d0))

		else if ( TGR%DE_eqstate == 3 ) then
		grhoa2 = grhoa2 + grhov*a**(1.d0-3.d0*(TGR%wp+TGR%a_p*TGR%wa))*exp(3.d0*TGR%wa*(a-1.d0))

		end if	

	end if	
	!<ISiTGR MOD END
	
    if (CP%Num_Nu_massive /= 0) then
        !Get massive neutrino density relative to massless
        do nu_i = 1, CP%nu_mass_eigenstates
            call Nu_rho(a*nu_masses(nu_i),rhonu)
            grhoa2=grhoa2+rhonu*grhormass(nu_i)
        end do
    end if

    dtauda=sqrt(3/grhoa2)

    end function dtauda

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !Gauge-dependent perturbation equations

    module GaugeInterface
    use precision
    use ModelParams
    use MassiveNu
    use LambdaGeneral
	!>ISiTGR MOD START
    use constants
    use ISiTGR
	!<ISiTGR MOD END
    use Errors
    use Transfer
    implicit none
    public

    !Description of this file. Change if you make modifications.
    character(LEN=*), parameter :: Eqns_name = 'gauge_inv'

    integer, parameter :: basic_num_eqns = 5

    logical :: DoTensorNeutrinos = .true.

    logical :: DoLateRadTruncation = .true.
    !if true, use smooth approx to radition perturbations after decoupling on
    !small scales, saving evolution of irrelevant osciallatory multipole equations

    logical, parameter :: second_order_tightcoupling = .true.

    real(dl) :: Magnetic = 0._dl
    !Vector mode anisotropic stress in units of rho_gamma
    real(dl) :: vec_sig0 = 1._dl
    !Vector mode shear
    integer, parameter :: max_l_evolve = 256 !Maximum l we are ever likely to propagate
    !Note higher values increase size of Evolution vars, hence memoryWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

    !Supported scalar initial condition flags
    integer, parameter :: initial_adiabatic=1, initial_iso_CDM=2, &
        initial_iso_baryon=3,  initial_iso_neutrino=4, initial_iso_neutrino_vel=5, initial_vector = 0
    integer, parameter :: initial_nummodes =  initial_iso_neutrino_vel

    type EvolutionVars
        real(dl) q, q2
        real(dl) k_buf,k2_buf ! set in initial

        integer w_ix !Index of two quintessence equations
        integer r_ix !Index of the massless neutrino hierarchy
        integer g_ix !Index of the photon neutrino hierarchy

        integer q_ix !index into q_evolve array that gives the value q
        logical TransferOnly

        !       nvar  - number of scalar (tensor) equations for this k
        integer nvar,nvart, nvarv

        !Max_l for the various hierarchies
        integer lmaxg,lmaxnr,lmaxnu,lmaxgpol,MaxlNeeded
        integer lmaxnrt, lmaxnut, lmaxt, lmaxpolt, MaxlNeededt
        logical EvolveTensorMassiveNu(max_nu)
        integer lmaxnrv, lmaxv, lmaxpolv

        integer polind  !index into scalar array of polarization hierarchy

        !array indices for massive neutrino equations
        integer nu_ix(max_nu), nu_pert_ix
        integer nq(max_nu), lmaxnu_pert
        logical has_nu_relativistic

        !Initial values for massive neutrino v*3 variables calculated when switching
        !to non-relativistic approx
        real(dl) G11(max_nu),G30(max_nu)
        !True when using non-relativistic approximation
        logical MassiveNuApprox(max_nu)
        real(dl) MassiveNuApproxTime(max_nu)

        !True when truncating at l=2,3 when k*tau>>1 (see arXiv:1104.2933)
        logical high_ktau_neutrino_approx

        !Massive neutrino scheme being used at the moment
        integer NuMethod

        !True when using tight-coupling approximation (required for stability at early times)
        logical TightCoupling, TensTightCoupling
        real(dl) TightSwitchoffTime

        !Numer of scalar equations we are propagating
        integer ScalEqsToPropagate
        integer TensEqsToPropagate
        !beta > l for closed models
        integer FirstZerolForBeta
        !Tensor vars
        real(dl) aux_buf

        real(dl) pig, pigdot
        real(dl) poltruncfac

        logical no_nu_multpoles, no_phot_multpoles
        integer lmaxnu_tau(max_nu)  !lmax for massive neutinos at time being integrated
        logical nu_nonrelativistic(max_nu)

        real(dl) denlk(max_l_evolve),denlk2(max_l_evolve), polfack(max_l_evolve)
        real(dl) Kf(max_l_evolve)

        integer E_ix, B_ix !tensor polarizatisdon indices
        real(dl) denlkt(4,max_l_evolve),Kft(max_l_evolve)
        real, pointer :: OutputTransfer(:) => null()
        real(dl), pointer :: OutputSources(:) => null()
        real(dl), pointer :: CustomSources(:) => null()

    end type EvolutionVars

    ABSTRACT INTERFACE
    SUBROUTINE TSource_func(sources, tau, a, adotoa, grho, gpres,w_lam, cs2_lam,  &
        grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t, &
        k,etak, etakdot, phi, phidot, sigma, sigmadot, &
        dgrho, clxg,clxb,clxc,clxr, clxnu, clxde, delta_p_b, &
        dgq, qg, qr, qde, vb, qgdot, qrdot, vbdot, &
        dgpi, pig, pir, pigdot, pirdot, diff_rhopi, &
        polter, polterdot, polterddot, octg, octgdot, E, Edot, &
        opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau, &
        tau0, tau_maxvis, Kf, f_K)
    real*8, intent(out) :: sources(:)
    real*8, intent(in) :: tau, a, adotoa, grho, gpres,w_lam, cs2_lam,  &
        grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t, &
        k,etak, etakdot, phi, phidot, sigma, sigmadot, &
        dgrho, clxg,clxb,clxc, clxr, clxnu, clxde, delta_p_b, &
        dgq, qg, qr, qde, vb, qgdot, qrdot, vbdot, &
        dgpi, pig, pir, pigdot, pirdot, diff_rhopi, &
        polter, polterdot, polterddot, octg, octgdot, E(2:3), Edot(2:3), &
        opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau, &
        tau0, tau_maxvis
    REAL*8, intent(in) :: Kf(*)
    real*8, external :: f_K
    END SUBROUTINE TSource_func
    END INTERFACE

    procedure(TSource_func), pointer :: custom_sources_func => null()

    !precalculated arrays
    real(dl) polfac(max_l_evolve),denl(max_l_evolve),vecfac(max_l_evolve),vecfacpol(max_l_evolve)

    real(dl), parameter :: ep0=1.0d-2
    integer, parameter :: lmaxnu_high_ktau=4 !Jan2015, increased from 3 to fix mpk for mnu~6eV

    real(dl) epsw
    real(dl) nu_tau_notmassless(nqmax0+1,max_nu), nu_tau_nonrelativistic(max_nu),nu_tau_massive(max_nu)
    contains


    subroutine GaugeInterface_ScalEv(EV,y,tau,tauend,tol1,ind,c,w)
    type(EvolutionVars) EV
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar), tol1, tau, tauend
    integer ind

    call dverk(EV,EV%ScalEqsToPropagate,derivs,tau,y,tauend,tol1,ind,c,EV%nvar,w)
    if (ind==-3) then
        call GlobalError('Dverk error -3: the subroutine was unable  to  satisfy  the  error ' &
            //'requirement  with a particular step-size that is less than or * ' &
            //'equal to hmin, which may mean that tol is too small' &
            //'--- but most likely you''ve messed up the y array indexing; ' &
            //'compiling with bounds checking may (or may not) help find the problem.',error_evolution)
    end if
    end subroutine GaugeInterface_ScalEv

    function next_nu_nq(nq) result (next_nq)
    integer, intent(in) :: nq
    integer q, next_nq

    if (nq==0) then
        next_nq=1
    else
        q = nu_q(nq)
        if (q>=10) then
            next_nq = nqmax
        else
            next_nq = nq+1
        end if
    end if

    end function next_nu_nq

    recursive subroutine GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
    use ThermoData
    type(EvolutionVars) EV, EVout
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar), yout(EV%nvar), tol1, tau, tauend
    integer ind, nu_i
    real(dl) cs2, opacity, dopacity
    real(dl) tau_switch_ktau, tau_switch_nu_massless, tau_switch_nu_massive, next_switch
    real(dl) tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles,tau_switch_nu_nonrel
    real(dl) noSwitch, smallTime

    noSwitch= CP%tau0+1
    smallTime =  min(tau, 1/EV%k_buf)/100

    tau_switch_ktau = noSwitch
    tau_switch_no_nu_multpoles= noSwitch
    tau_switch_no_phot_multpoles= noSwitch

    !Massive neutrino switches
    tau_switch_nu_massless = noSwitch
    tau_switch_nu_nonrel = noSwitch
    tau_switch_nu_massive= noSwitch

    !Evolve equations from tau to tauend, performing switches in equations if necessary.

    if (.not. EV%high_ktau_neutrino_approx .and. .not. EV%no_nu_multpoles ) then
        tau_switch_ktau=  max(20, EV%lmaxnr-4)/EV%k_buf
    end if

    if (CP%Num_Nu_massive /= 0) then
        do nu_i = 1, CP%Nu_mass_eigenstates
            if (EV%nq(nu_i) /= nqmax) then
                tau_switch_nu_massless = min(tau_switch_nu_massless,nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i))
            else if (.not. EV%nu_nonrelativistic(nu_i)) then
                tau_switch_nu_nonrel = min(nu_tau_nonrelativistic(nu_i),tau_switch_nu_nonrel)
            else if (EV%NuMethod==Nu_trunc .and..not. EV%MassiveNuApprox(nu_i)) then
                tau_switch_nu_massive = min(tau_switch_nu_massive,EV%MassiveNuApproxTime(nu_i))
            end if
        end do
    end if
    !>ISiTGR MOD START
    if (DoLateRadTruncation) then
        if (.not. EV%no_nu_multpoles) & !!.and. .not. EV%has_nu_relativistic .and. tau_switch_nu_massless ==noSwitch)  &
            tau_switch_no_nu_multpoles=max(15/EV%k_buf*AccuracyBoost,min(taurend,matter_verydom_tau))
    !>ISiTGR MOD END
        if (.not. EV%no_phot_multpoles .and. (.not.CP%WantCls .or. EV%k_buf>0.03*AccuracyBoost)) &
            tau_switch_no_phot_multpoles =max(15/EV%k_buf,taurend)*AccuracyBoost
    end if

    next_switch = min(tau_switch_ktau, tau_switch_nu_massless,EV%TightSwitchoffTime, tau_switch_nu_massive, &
        tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles, tau_switch_nu_nonrel,noSwitch)

    if (next_switch < tauend) then
        if (next_switch > tau+smallTime) then
            call GaugeInterface_ScalEv(EV, y, tau,next_switch,tol1,ind,c,w)
            if (global_error_flag/=0) return
        end if

        EVout=EV

        if (next_switch == EV%TightSwitchoffTime) then
            !TightCoupling
            EVout%TightCoupling=.false.
            EVout%TightSwitchoffTime = noSwitch
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            EV=EVout
            y=yout
            ind=1
            !Set up variables with their tight coupling values
            y(EV%g_ix+2) = EV%pig
            call thermo(tau,cs2,opacity,dopacity)

            if (second_order_tightcoupling) then
                ! Francis-Yan Cyr-Racine November 2010

                y(EV%g_ix+3) = (3._dl/7._dl)*y(EV%g_ix+2)*(EV%k_buf/opacity)*(1._dl+dopacity/opacity**2) + &
                    (3._dl/7._dl)*EV%pigdot*(EV%k_buf/opacity**2)*(-1._dl)

                y(EV%polind+2) = EV%pig/4 + EV%pigdot*(1._dl/opacity)*(-5._dl/8._dl- &
                    (25._dl/16._dl)*dopacity/opacity**2) + &
                    EV%pig*(EV%k_buf/opacity)**2*(-5._dl/56._dl)
                y(EV%polind+3) = (3._dl/7._dl)*(EV%k_buf/opacity)*y(EV%polind+2)*(1._dl + &
                    dopacity/opacity**2) + (3._dl/7._dl)*(EV%k_buf/opacity**2)*((EV%pigdot/4._dl)* &
                    (1._dl+(5._dl/2._dl)*dopacity/opacity**2))*(-1._dl)
            else
                y(EV%g_ix+3) = 3./7*y(EV%g_ix+2)*EV%k_buf/opacity
                y(EV%polind+2) = EV%pig/4
                y(EV%polind+3) =y(EV%g_ix+3)/4
            end if
        else if (next_switch==tau_switch_ktau) then
            !k tau >> 1, evolve massless neutrino effective fluid up to l=2
            EVout%high_ktau_neutrino_approx=.true.
            EVout%nq(1:CP%Nu_mass_eigenstates) = nqmax
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch == tau_switch_nu_massless) then
            !Mass starts to become important, start evolving next momentum mode
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (EV%nq(nu_i) /= nqmax .and. &
                    next_switch == nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i)) then
                EVOut%nq(nu_i) = next_nu_nq(EV%nq(nu_i))
                call SetupScalarArrayIndices(EVout)
                call CopyScalarVariableArray(y,yout, EV, EVout)
                EV=EVout
                y=yout
                exit
                end if
            end do
        else if (next_switch == tau_switch_nu_nonrel) then
            !Neutrino becomes non-relativistic, don't need high L
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (.not. EV%nu_nonrelativistic(nu_i) .and.  next_switch==nu_tau_nonrelativistic(nu_i) ) then
                    EVout%nu_nonrelativistic(nu_i)=.true.
                    call SetupScalarArrayIndices(EVout)
                    call CopyScalarVariableArray(y,yout, EV, EVout)
                    EV=EVout
                    y=yout
                    exit
                end if
            end do
        else if (next_switch == tau_switch_nu_massive) then
            !Very non-relativistic neutrinos, switch to truncated velocity-weight hierarchy
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (.not. EV%MassiveNuApprox(nu_i) .and.  next_switch== EV%MassiveNuApproxTime(nu_i) ) then
                    call SwitchToMassiveNuApprox(EV,y, nu_i)
                    exit
                end if
            end do
        else if (next_switch==tau_switch_no_nu_multpoles) then
            !Turn off neutrino hierarchies at late time where slow and not needed.
            ind=1
            EVout%no_nu_multpoles=.true.
            EVOut%nq(1:CP%Nu_mass_eigenstates ) = nqmax
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch==tau_switch_no_phot_multpoles) then
            !Turn off photon hierarchies at late time where slow and not needed.
            ind=1
            EVout%no_phot_multpoles=.true.
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        end if

        call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
        return
    end if

    call GaugeInterface_ScalEv(EV,y,tau,tauend,tol1,ind,c,w)

    end subroutine GaugeInterface_EvolveScal

    subroutine GaugeInterface_EvolveTens(EV,tau,y,tauend,tol1,ind,c,w)
    use ThermoData
    type(EvolutionVars) EV, EVOut
    real(dl) c(24),w(EV%nvart,9), y(EV%nvart),yout(EV%nvart), tol1, tau, tauend
    integer ind
    real(dl) opacity, cs2

    if (EV%TensTightCoupling .and. tauend > EV%TightSwitchoffTime) then
        if (EV%TightSwitchoffTime > tau) then
            call dverk(EV,EV%TensEqsToPropagate, derivst,tau,y,EV%TightSwitchoffTime,tol1,ind,c,EV%nvart,w)
        end if
        EVOut=EV
        EVOut%TensTightCoupling = .false.
        call SetupTensorArrayIndices(EVout)
        call CopyTensorVariableArray(y,yout,Ev, Evout)
        Ev = EvOut
        y=yout
        call thermo(tau,cs2,opacity)
        y(EV%g_ix+2)= 32._dl/45._dl*EV%k_buf/opacity*y(3)
        y(EV%E_ix+2) = y(EV%g_ix+2)/4
    end if

    call dverk(EV,EV%TensEqsToPropagate, derivst,tau,y,tauend,tol1,ind,c,EV%nvart,w)


    end subroutine GaugeInterface_EvolveTens

    function DeltaTimeMaxed(a1,a2, tol) result(t)
    real(dl) a1,a2,t
    real(dl), optional :: tol
    if (a1>1._dl) then
        t=0
    elseif (a2 > 1._dl) then
        t = DeltaTime(a1,1.01_dl, tol)
    else
        t= DeltaTime(a1,a2, tol)
    end if
    end function DeltaTimeMaxed

    subroutine GaugeInterface_Init
    !Precompute various arrays and other things independent of wavenumber
    integer j, nu_i
    real(dl) a_nonrel, a_mass,a_massive, time, nu_mass

    epsw = 100/CP%tau0

    if (CP%WantScalars) then
        do j=2,max_l_evolve
            polfac(j)=real((j+3)*(j-1),dl)/(j+1)
        end do
    end if

    if (CP%WantVectors) then
        do j=2,max_l_evolve
            vecfac(j)=real((j+2),dl)/(j+1)
            vecfacpol(j)=real((j+3)*j,dl)*(j-1)*vecfac(j)/(j+1)**2
        end do
    end if

    do j=1,max_l_evolve
        denl(j)=1._dl/(2*j+1)
    end do

    do nu_i=1, CP%Nu_Mass_eigenstates
        nu_mass = max(0.1_dl,nu_masses(nu_i))
        a_mass =  1.e-1_dl/nu_mass/lAccuracyBoost
        !if (HighAccuracyDefault) a_mass=a_mass/4
        time=DeltaTime(0._dl,nu_q(1)*a_mass)
        nu_tau_notmassless(1, nu_i) = time
        do j=2,nqmax
            !times when each momentum mode becomes signficantly nonrelativistic
            time= time + DeltaTimeMaxed(nu_q(j-1)*a_mass,nu_q(j)*a_mass, 0.01_dl)
            nu_tau_notmassless(j, nu_i) = time
        end do

        a_nonrel =  2.5d0/nu_mass*AccuracyBoost !!!Feb13tweak
        nu_tau_nonrelativistic(nu_i) =DeltaTimeMaxed(0._dl,a_nonrel)
        a_massive =  17.d0/nu_mass*AccuracyBoost
        nu_tau_massive(nu_i) =nu_tau_nonrelativistic(nu_i) + DeltaTimeMaxed(a_nonrel,a_massive)
    end do

    end subroutine GaugeInterface_Init


    subroutine SetupScalarArrayIndices(EV, max_num_eqns)
    !Set up array indices after the lmax have been decided
    use MassiveNu
	!>ISiTGR MOD START
	use ISiTGR
	!<ISiTGR MOD END
    !Set the numer of equations in each hierarchy, and get total number of equations for this k
    type(EvolutionVars) EV
    integer, intent(out), optional :: max_num_eqns
    integer neq, maxeq, nu_i

    neq=basic_num_eqns
    maxeq=neq
    if (.not. EV%no_phot_multpoles) then
        !Photon multipoles
        EV%g_ix=basic_num_eqns+1
        if (EV%TightCoupling) then
            neq=neq+2
        else
            neq = neq+ (EV%lmaxg+1)
            !Polarization multipoles
            EV%polind = neq -1 !polind+2 is L=2, for polarizationthe first calculated
            neq=neq + EV%lmaxgpol-1
        end if
    end if
    if (.not. EV%no_nu_multpoles) then
        !Massless neutrino multipoles
        EV%r_ix= neq+1
        if (EV%high_ktau_neutrino_approx) then
            neq=neq + 3
        else
            neq=neq + (EV%lmaxnr+1)
        end if
    end if
    maxeq = maxeq +  (EV%lmaxg+1)+(EV%lmaxnr+1)+EV%lmaxgpol-1

    !Dark energy
	!>ISiTGR MOD START
    if (w_lam /= -1 .and. w_Perturb .and. TGR%GR==0) then
        EV%w_ix = neq+1
        neq=neq+2
        maxeq=maxeq+2
    else
       	EV%w_ix=0
    end if
	!<ISiTGR MOD END

	!Massive neutrinos
    if (CP%Num_Nu_massive /= 0) then
        EV%has_nu_relativistic = any(EV%nq(1:CP%Nu_Mass_eigenstates)/=nqmax)
        if (EV%has_nu_relativistic) then
            EV%lmaxnu_pert=EV%lmaxnu
            EV%nu_pert_ix=neq+1
            neq = neq+ EV%lmaxnu_pert+1
            maxeq=maxeq+ EV%lmaxnu_pert+1
        else
            EV%lmaxnu_pert=0
        end if

        do nu_i=1, CP%Nu_Mass_eigenstates
            if (EV%high_ktau_neutrino_approx) then
                EV%lmaxnu_tau(nu_i) = lmaxnu_high_ktau *lAccuracyBoost
                if (CP%Transfer%accurate_massive_neutrinos) EV%lmaxnu_tau(nu_i) = EV%lmaxnu_tau(nu_i) *3
            else
                EV%lmaxnu_tau(nu_i) =max(min(nint(0.8_dl*EV%q*nu_tau_nonrelativistic(nu_i)*lAccuracyBoost),EV%lmaxnu),3)
                !!!Feb13tweak
                if (EV%nu_nonrelativistic(nu_i)) EV%lmaxnu_tau(nu_i)=min(EV%lmaxnu_tau(nu_i),nint(4*lAccuracyBoost))
            end if
            if (nu_masses(nu_i) > 5000 .and. CP%Transfer%high_precision) EV%lmaxnu_tau(nu_i) = EV%lmaxnu_tau(nu_i)*2 !megadamping
            EV%lmaxnu_tau(nu_i)=min(EV%lmaxnu,EV%lmaxnu_tau(nu_i))

            EV%nu_ix(nu_i)=neq+1
            if (EV%MassiveNuApprox(nu_i)) then
                neq = neq+4
            else
                neq = neq+ EV%nq(nu_i)*(EV%lmaxnu_tau(nu_i)+1)
            endif
            maxeq = maxeq + nqmax*(EV%lmaxnu+1)
        end do
    else
        EV%has_nu_relativistic = .false.
    end if

    EV%ScalEqsToPropagate = neq
    if (present(max_num_eqns)) then
        max_num_eqns=maxeq
    end if

    end subroutine SetupScalarArrayIndices

    subroutine CopyScalarVariableArray(y,yout, EV, EVout)
    type(EvolutionVars) EV, EVOut
    real(dl), intent(in) :: y(EV%nvar)
    real(dl), intent(out) :: yout(EVout%nvar)
    integer lmax,i, nq
    integer nnueq,nu_i, ix_off, ix_off2, ind, ind2
    real(dl) q, pert_scale

    yout=0
    yout(1:basic_num_eqns) = y(1:basic_num_eqns)
    if (w_lam /= -1 .and. w_Perturb) then
        yout(EVout%w_ix)=y(EV%w_ix)
        yout(EVout%w_ix+1)=y(EV%w_ix+1)
    end if

    if (.not. EV%no_phot_multpoles .and. .not. EVout%no_phot_multpoles) then
        if (EV%TightCoupling .or. EVOut%TightCoupling) then
            lmax=1
        else
            lmax = min(EV%lmaxg,EVout%lmaxg)
        end if
        yout(EVout%g_ix:EVout%g_ix+lmax)=y(EV%g_ix:EV%g_ix+lmax)
        if (.not. EV%TightCoupling .and. .not. EVOut%TightCoupling) then
            lmax = min(EV%lmaxgpol,EVout%lmaxgpol)
            yout(EVout%polind+2:EVout%polind+lmax)=y(EV%polind+2:EV%polind+lmax)
        end if
    end if

    if (.not. EV%no_nu_multpoles .and. .not. EVout%no_nu_multpoles) then
        if (EV%high_ktau_neutrino_approx .or. EVout%high_ktau_neutrino_approx) then
            lmax=2
        else
            lmax = min(EV%lmaxnr,EVout%lmaxnr)
        end if
        yout(EVout%r_ix:EVout%r_ix+lmax)=y(EV%r_ix:EV%r_ix+lmax)
    end if

    if (CP%Num_Nu_massive /= 0) then
        do nu_i=1,CP%Nu_mass_eigenstates
            ix_off=EV%nu_ix(nu_i)
            ix_off2=EVOut%nu_ix(nu_i)
            if (EV%MassiveNuApprox(nu_i) .and. EVout%MassiveNuApprox(nu_i)) then
                nnueq=4
                yout(ix_off2:ix_off2+nnueq-1)=y(ix_off:ix_off+nnueq-1)
            else if (.not. EV%MassiveNuApprox(nu_i) .and. .not. EVout%MassiveNuApprox(nu_i)) then
                lmax=min(EV%lmaxnu_tau(nu_i),EVOut%lmaxnu_tau(nu_i))
                nq = min(EV%nq(nu_i), EVOut%nq(nu_i))
                do i=1,nq
                    ind= ix_off + (i-1)*(EV%lmaxnu_tau(nu_i)+1)
                    ind2=ix_off2+ (i-1)*(EVOut%lmaxnu_tau(nu_i)+1)
                    yout(ind2:ind2+lmax) = y(ind:ind+lmax)
                end do
                do i=nq+1, EVOut%nq(nu_i)
                    lmax = min(EVOut%lmaxnu_tau(nu_i), EV%lmaxnr)
                    ind2=ix_off2+ (i-1)*(EVOut%lmaxnu_tau(nu_i)+1)
                    yout(ind2:ind2+lmax) = y(EV%r_ix:EV%r_ix+lmax)

                    !Add leading correction for the mass
                    q=nu_q(i)
                    pert_scale=(nu_masses(nu_i)/q)**2/2
                    lmax = min(lmax,EV%lmaxnu_pert)
                    yout(ind2:ind2+lmax) = yout(ind2:ind2+lmax) &
                        + y(EV%nu_pert_ix:EV%nu_pert_ix+lmax)*pert_scale
                end do
            end if
        end do

        if (EVOut%has_nu_relativistic .and. EV%has_nu_relativistic) then
            lmax = min(EVOut%lmaxnu_pert, EV%lmaxnu_pert)
            yout(EVout%nu_pert_ix:EVout%nu_pert_ix+lmax)=  y(EV%nu_pert_ix:EV%nu_pert_ix+lmax)
        end if
    end if

    end subroutine CopyScalarVariableArray


    subroutine SetupTensorArrayIndices(EV, maxeq)
    type(EvolutionVars) EV
    integer nu_i, neq
    integer, optional, intent(out) :: maxeq
    neq=3
    EV%g_ix = neq-1 !EV%g_ix+2 is quadrupole
    if (.not. EV%TensTightCoupling) then
        EV%E_ix = EV%g_ix + (EV%lmaxt-1)
        EV%B_ix = EV%E_ix + (EV%lmaxpolt-1)
        neq = neq+ (EV%lmaxt-1)+(EV%lmaxpolt-1)*2
    end if
    if (present(maxeq)) then
        maxeq =3 + (EV%lmaxt-1)+(EV%lmaxpolt-1)*2
    end if
    EV%r_ix = neq -1
    if (DoTensorNeutrinos) then
        neq = neq + EV%lmaxnrt-1
        if (present(maxeq)) maxeq = maxeq+EV%lmaxnrt-1
        if (CP%Num_Nu_massive /= 0 ) then
            do nu_i=1, CP%nu_mass_eigenstates
                EV%EvolveTensorMassiveNu(nu_i) = nu_tau_nonrelativistic(nu_i) < 0.8*tau_maxvis*AccuracyBoost
                if (EV%EvolveTensorMassiveNu(nu_i)) then
                    EV%nu_ix(nu_i)=neq-1
                    neq = neq+ nqmax*(EV%lmaxnut-1)
                    if (present(maxeq)) maxeq = maxeq + nqmax*(EV%lmaxnut-1)
                end if
            end do
        end if
    end if

    EV%TensEqsToPropagate = neq

    end  subroutine SetupTensorArrayIndices

    subroutine CopyTensorVariableArray(y,yout, EV, EVout)
    type(EvolutionVars) EV, EVOut
    real(dl), intent(in) :: y(EV%nvart)
    real(dl), intent(out) :: yout(EVout%nvart)
    integer lmaxpolt, lmaxt, nu_i, ind, ind2, i

    yout=0
    yout(1:3) = y(1:3)
    if (.not. EVOut%TensTightCoupling .and. .not.EV%TensTightCoupling) then
        lmaxt = min(EVOut%lmaxt,EV%lmaxt)
        yout(EVout%g_ix+2:EVout%g_ix+lmaxt)=y(EV%g_ix+2:EV%g_ix+lmaxt)
        lmaxpolt = min(EV%lmaxpolt, EVOut%lmaxpolt)
        yout(EVout%E_ix+2:EVout%E_ix+lmaxpolt)=y(EV%E_ix+2:EV%E_ix+lmaxpolt)
        yout(EVout%B_ix+2:EVout%B_ix+lmaxpolt)=y(EV%B_ix+2:EV%B_ix+lmaxpolt)
    end if
    if (DoTensorNeutrinos) then
        lmaxt=min(EV%lmaxnrt,EVOut%lmaxnrt)
        yout(EVout%r_ix+2:EVout%r_ix+lmaxt)=y(EV%r_ix+2:EV%r_ix+lmaxt)
        do nu_i =1, CP%nu_mass_eigenstates
            if (EV%EvolveTensorMassiveNu(nu_i)) then
                lmaxt=min(EV%lmaxnut,EVOut%lmaxnut)
                do i=1,nqmax
                    ind= EV%nu_ix(nu_i) + (i-1)*(EV%lmaxnut-1)
                    ind2=EVOut%nu_ix(nu_i)+ (i-1)*(EVOut%lmaxnut-1)
                    yout(ind2+2:ind2+lmaxt) = y(ind+2:ind+lmaxt)
                end do
            end if
        end do
    end if

    end subroutine CopyTensorVariableArray

    subroutine GetNumEqns(EV)
    use MassiveNu
    !Set the numer of equations in each hierarchy, and get total number of equations for this k
    type(EvolutionVars) EV
    real(dl) scal, max_nu_mass
    integer nu_i,q_rel,j

    if (CP%Num_Nu_massive == 0) then
        EV%lmaxnu=0
        max_nu_mass=0
    else
        max_nu_mass = maxval(nu_masses(1:CP%Nu_mass_eigenstates))
        do nu_i = 1, CP%Nu_mass_eigenstates
            !Start with momentum modes for which t_k ~ time at which mode becomes non-relativistic
            q_rel=0
            do j=1, nqmax
                !two different q's here EV%q ~k
                if (nu_q(j) > nu_masses(nu_i)*adotrad/EV%q) exit
                q_rel = q_rel + 1
            end do

            if (q_rel>= nqmax-2 .or. CP%WantTensors) then
                EV%nq(nu_i)=nqmax
            else
                EV%nq(nu_i)=q_rel
            end if
            !q_rel = nint(nu_masses(nu_i)*adotrad/EV%q) !two dffierent q's here EV%q ~k
            !EV%nq(nu_i)=max(0,min(nqmax0,q_rel)) !number of momentum modes to evolve intitially
            EV%nu_nonrelativistic(nu_i) = .false.
        end do

        EV%NuMethod = CP%MassiveNuMethod
        if (EV%NuMethod == Nu_Best) EV%NuMethod = Nu_Trunc
        !l_max for massive neutrinos
        if (CP%Transfer%high_precision) then
            EV%lmaxnu=nint(25*lAccuracyBoost)
        else
            EV%lmaxnu=max(3,nint(10*lAccuracyBoost))
            if (max_nu_mass>700) EV%lmaxnu=max(3,nint(15*lAccuracyBoost)) !Feb13 tweak
        endif
    end if

    if (CP%closed) then
        EV%FirstZerolForBeta = nint(EV%q*CP%r)
    else
        EV%FirstZerolForBeta=l0max !a large number
    end if

    EV%high_ktau_neutrino_approx = .false.
    if (CP%WantScalars) then
        EV%TightCoupling=.true.
        EV%no_phot_multpoles =.false.
        EV%no_nu_multpoles =.false.
        EV%MassiveNuApprox=.false.

        if (HighAccuracyDefault .and. CP%AccuratePolarization) then
            EV%lmaxg  = max(nint(11*lAccuracyBoost),3)
        else
            EV%lmaxg  = max(nint(8*lAccuracyBoost),3)
        end if
        EV%lmaxnr = max(nint(14*lAccuracyBoost),3)
        if (max_nu_mass>700 .and. HighAccuracyDefault) EV%lmaxnr = max(nint(32*lAccuracyBoost),3) !Feb13 tweak

        EV%lmaxgpol = EV%lmaxg
        if (.not.CP%AccuratePolarization) EV%lmaxgpol=max(nint(4*lAccuracyBoost),3)

        if (EV%q < 0.05) then
            !Large scales need fewer equations
            scal  = 1
            if (CP%AccuratePolarization) scal = 4  !But need more to get polarization right
            EV%lmaxgpol=max(3,nint(min(8,nint(scal* 150* EV%q))*lAccuracyBoost))
            EV%lmaxnr=max(3,nint(min(7,nint(sqrt(scal)* 150 * EV%q))*lAccuracyBoost))
            EV%lmaxg=max(3,nint(min(8,nint(sqrt(scal) *300 * EV%q))*lAccuracyBoost))
            if (CP%AccurateReionization) then
                EV%lmaxg=EV%lmaxg*4
                EV%lmaxgpol=EV%lmaxgpol*2
            end if
        end if

        if (EV%TransferOnly) then
            EV%lmaxgpol = min(EV%lmaxgpol,nint(5*lAccuracyBoost))
            EV%lmaxg = min(EV%lmaxg,nint(6*lAccuracyBoost))
        end if
        if (CP%Transfer%high_precision) then
            if (HighAccuracyDefault) then
                EV%lmaxnr=max(nint(45*lAccuracyBoost),3)
            else
                EV%lmaxnr=max(nint(30*lAccuracyBoost),3)
            endif
            if (EV%q > 0.04 .and. EV%q < 0.5) then !baryon oscillation scales
                EV%lmaxg=max(EV%lmaxg,10)
            end if
        end if

        if (CP%closed) then
            EV%lmaxnu=min(EV%lmaxnu, EV%FirstZerolForBeta-1)
            EV%lmaxnr=min(EV%lmaxnr, EV%FirstZerolForBeta-1)
            EV%lmaxg=min(EV%lmaxg, EV%FirstZerolForBeta-1)
            EV%lmaxgpol=min(EV%lmaxgpol, EV%FirstZerolForBeta-1)
        end if

        EV%poltruncfac=real(EV%lmaxgpol,dl)/max(1,(EV%lmaxgpol-2))
        EV%MaxlNeeded=max(EV%lmaxg,EV%lmaxnr,EV%lmaxgpol,EV%lmaxnu)
        if (EV%MaxlNeeded > max_l_evolve) call MpiStop('Need to increase max_l_evolve')
        call SetupScalarArrayIndices(EV,EV%nvar)
        if (CP%closed) EV%nvar=EV%nvar+1 !so can reference lmax+1 with zero coefficient
        EV%lmaxt=0
    else
        EV%nvar=0
    end if

    if (CP%WantTensors) then
        EV%TensTightCoupling = .true.
        EV%lmaxt=max(3,nint(8*lAccuracyBoost))
        EV%lmaxpolt = max(3,nint(4*lAccuracyBoost))
        ! if (EV%q < 1e-3) EV%lmaxpolt=EV%lmaxpolt+1
        if (DoTensorNeutrinos) then
            EV%lmaxnrt=nint(6*lAccuracyBoost)
            EV%lmaxnut=EV%lmaxnrt
        else
            EV%lmaxnut=0
            EV%lmaxnrt=0
        end if
        if (CP%closed) then
            EV%lmaxt=min(EV%FirstZerolForBeta-1,EV%lmaxt)
            EV%lmaxpolt=min(EV%FirstZerolForBeta-1,EV%lmaxpolt)
            EV%lmaxnrt=min(EV%FirstZerolForBeta-1,EV%lmaxnrt)
            EV%lmaxnut=min(EV%FirstZerolForBeta-1,EV%lmaxnut)
        end if
        EV%MaxlNeededt=max(EV%lmaxpolt,EV%lmaxt, EV%lmaxnrt, EV%lmaxnut)
        if (EV%MaxlNeededt > max_l_evolve) call MpiStop('Need to increase max_l_evolve')
        call SetupTensorArrayIndices(EV, EV%nvart)
    else
        EV%nvart=0
    end if


    if (CP%WantVectors) then
        EV%lmaxv=max(10,nint(8*lAccuracyBoost))
        EV%lmaxpolv = max(5,nint(5*lAccuracyBoost))

        EV%nvarv=(EV%lmaxv)+(EV%lmaxpolv-1)*2+3

        EV%lmaxnrv=nint(30*lAccuracyBoost)

        EV%nvarv=EV%nvarv+EV%lmaxnrv
        if (CP%Num_Nu_massive /= 0 ) then
            call MpiStop('massive neutrinos not supported for vector modes')
        end if
    else
        EV%nvarv=0
    end if

    end subroutine GetNumEqns

    !cccccccccccccccccccccccccccccccccc
    subroutine SwitchToMassiveNuApprox(EV,y, nu_i)
    !When the neutrinos are no longer highly relativistic we use a truncated
    !energy-integrated hierarchy going up to third order in velocity dispersion
    type(EvolutionVars) EV, EVout
    integer, intent(in) :: nu_i

    real(dl) a,a2,pnu,clxnu,dpnu,pinu,rhonu
    real(dl) qnu
    real(dl) y(EV%nvar), yout(EV%nvar)

    a=y(1)
    a2=a*a
    EVout=EV
    EVout%MassiveNuApprox(nu_i)=.true.
    call SetupScalarArrayIndices(EVout)
    call CopyScalarVariableArray(y,yout, EV, EVout)

    !Get density and pressure as ratio to massles by interpolation from table
    call Nu_background(a*nu_masses(nu_i),rhonu,pnu)

    !Integrate over q
    call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu,dpnu,pinu)
    !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
    dpnu=dpnu/rhonu
    qnu=qnu/rhonu
    clxnu = clxnu/rhonu
    pinu=pinu/rhonu

    yout(EVout%nu_ix(nu_i))=clxnu
    yout(EVout%nu_ix(nu_i)+1)=dpnu
    yout(EVout%nu_ix(nu_i)+2)=qnu
    yout(EVout%nu_ix(nu_i)+3)=pinu

    call Nu_Intvsq(EV,y, a, nu_i, EVout%G11(nu_i),EVout%G30(nu_i))
    !Analytic solution for higher moments, proportional to a^{-3}
    EVout%G11(nu_i)=EVout%G11(nu_i)*a2*a/rhonu
    EVout%G30(nu_i)=EVout%G30(nu_i)*a2*a/rhonu

    EV=EVout
    y=yout

    end subroutine SwitchToMassiveNuApprox
	
	!>ISiTGR MOD START: Modifying subroutine to compute extra term
subroutine MassiveNuVarsOut(EV,y,yprime,a,grho,gpres,dgrho,dgq,dgpi,dgpi_diff,pidot_sum,clxnu_all,dgpi_3wplus1,dgpi_3wplus2, &
dgpi_3wplus1plusbetak) !CGQ
	!<ISiTGR MOD END
    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar), yprime(EV%nvar),a
    real(dl), optional :: grho,gpres,dgrho,dgq,dgpi, dgpi_diff,pidot_sum,clxnu_all
	!>ISiTGR MOD START: adding new terms that contributes to MG
	real(dl), optional :: dgpi_3wplus1, dgpi_3wplus2, dgpi_3wplus1plusbetak!CGQ
	!<ISiTGR MOD END
    !grho = a^2 kappa rho
    !gpres = a^2 kappa p
    !dgrho = a^2 kappa \delta\rho
    !dgp =  a^2 kappa \delta p
    !dgq = a^2 kappa q (heat flux)
    !dgpi = a^2 kappa pi (anisotropic stress)
    !dgpi_diff = a^2 kappa (3*p -rho)*pi
	!dgpi_3wplus1 = a^2 kappa (3w+1)*rho*pi !CGQ
	!dgpi_3wplus2 = a^2 kappa (3w+2)*rho*pi !CGQ
	!dgpi_3wplus1plusbetak = a^2 kappa (3w+1+betak)*rho*pi !CGQ

    integer nu_i
    real(dl) pinudot,grhormass_t, rhonu, pnu,  rhonudot
    real(dl) adotoa, grhonu_t,gpnu_t
    real(dl) clxnu, qnu, pinu, dpnu, grhonu, dgrhonu
    real(dl) dtauda

    grhonu=0
    dgrhonu=0
    do nu_i = 1, CP%Nu_mass_eigenstates
        grhormass_t=grhormass(nu_i)/a**2

        !Get density and pressure as ratio to massless by interpolation from table
        call Nu_background(a*nu_masses(nu_i),rhonu,pnu)

        if (EV%MassiveNuApprox(nu_i)) then
            clxnu=y(EV%nu_ix(nu_i))
            !dpnu = y(EV%iq0+1+off_ix)
            qnu=y(EV%nu_ix(nu_i)+2)
            pinu=y(EV%nu_ix(nu_i)+3)
            pinudot=yprime(EV%nu_ix(nu_i)+3)
        else
            !Integrate over q
            call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu,dpnu,pinu)
            !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
            !dpnu=dpnu/rhonu
            qnu=qnu/rhonu
            clxnu = clxnu/rhonu
            pinu=pinu/rhonu
            adotoa = 1/(a*dtauda(a))
            rhonudot = Nu_drho(a*nu_masses(nu_i),adotoa,rhonu)

            call Nu_pinudot(EV,y, yprime, a,adotoa, nu_i,pinudot)
            pinudot=pinudot/rhonu - rhonudot/rhonu*pinu
        endif

        grhonu_t=grhormass_t*rhonu
        gpnu_t=grhormass_t*pnu

        grhonu = grhonu  + grhonu_t
        if (present(gpres)) gpres= gpres + gpnu_t

        dgrhonu= dgrhonu + grhonu_t*clxnu
        if (present(dgq)) dgq  = dgq   + grhonu_t*qnu
        if (present(dgpi)) dgpi = dgpi  + grhonu_t*pinu
        if (present(dgpi_diff)) dgpi_diff = dgpi_diff + pinu*(3*gpnu_t-grhonu_t)
        if (present(pidot_sum)) pidot_sum = pidot_sum + grhonu_t*pinudot
		!>ISiTGR MOD START: computing (3*w+2)*rho*Pi !CGQ
		if (present(dgpi_3wplus1)) dgpi_3wplus1 = dgpi_3wplus1 + grhonu_t*pinu * (3.d0*(pnu/rhonu) + 1.d0) 
        if (present(dgpi_3wplus2)) dgpi_3wplus2 = dgpi_3wplus2 + grhonu_t*pinu * (3.d0*(pnu/rhonu) + 2.d0) 
        if (present(dgpi_3wplus1plusbetak)) dgpi_3wplus1plusbetak = dgpi_3wplus1plusbetak + grhonu_t*pinu * &
			(3.d0*(pnu/rhonu) + 1.d0 + 1.d0/EV%Kf(1)) 
		!<ISiTGR MOD END
	end do
    if (present(grho)) grho = grho  + grhonu
    if (present(dgrho)) dgrho= dgrho + dgrhonu
    if (present(clxnu_all)) clxnu_all = dgrhonu/grhonu

end subroutine MassiveNuVarsOut

    subroutine Nu_Integrate_L012(EV,y,a,nu_i,drhonu,fnu,dpnu,pinu)
    type(EvolutionVars) EV
    !  Compute the perturbations of density and energy flux
    !  of one eigenstate of massive neutrinos, in units of the mean
    !  density of one eigenstate of massless neutrinos, by integrating over
    !  momentum.
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvar)
    real(dl), intent(OUT) ::  drhonu,fnu
    real(dl), optional, intent(OUT) :: dpnu,pinu
    real(dl) tmp, am, aq,v, pert_scale
    integer iq, ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.

    drhonu=0
    fnu=0
    if (present(dpnu)) then
        dpnu=0
        pinu=0
    end if
    am=a*nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    do iq=1,EV%nq(nu_i)
        aq=am/nu_q(iq)
        v=1._dl/sqrt(1._dl+aq*aq)
        drhonu=drhonu+ nu_int_kernel(iq)* y(ind)/v
        fnu=fnu+nu_int_kernel(iq)* y(ind+1)
        if (present(dpnu)) then
            dpnu=dpnu+  nu_int_kernel(iq)* y(ind)*v
            pinu=pinu+ nu_int_kernel(iq)*y(ind+2)*v
        end if
        ind=ind+EV%lmaxnu_tau(nu_i)+1
    end do
    ind = EV%nu_pert_ix
    do iq=EV%nq(nu_i)+1,nqmax
        !Get the rest from perturbatively relativistic expansion
        aq=am/nu_q(iq)
        v=1._dl/sqrt(1._dl+aq*aq)
        pert_scale=(nu_masses(nu_i)/nu_q(iq))**2/2
        tmp = nu_int_kernel(iq)*(y(EV%r_ix)  + pert_scale*y(ind)  )
        drhonu=drhonu+ tmp/v
        fnu=fnu+nu_int_kernel(iq)*(y(EV%r_ix+1)+ pert_scale*y(ind+1))
        if (present(dpnu)) then
            dpnu=dpnu+ tmp*v
            pinu = pinu+ nu_int_kernel(iq)*(y(EV%r_ix+2)+ pert_scale*y(ind+2))*v
        end if
    end do

    if (present(dpnu)) then
        dpnu = dpnu/3
    end if

    end subroutine Nu_Integrate_L012

    subroutine Nu_pinudot(EV,y, ydot, a,adotoa, nu_i,pinudot)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a,adotoa, y(EV%nvar), ydot(EV%nvar)

    !  Compute the time derivative of the mean density in massive neutrinos
    !  and the shear perturbation.
    real(dl) pinudot
    real(dl) aq,q,v,aqdot,vdot
    real(dl) psi2,psi2dot
    real(dl) am, pert_scale
    integer iq,ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    pinudot=0._dl
    ind=EV%nu_ix(nu_i)+2
    am=a*nu_masses(nu_i)
    do iq=1,EV%nq(nu_i)
        q=nu_q(iq)
        aq=am/q
        aqdot=aq*adotoa
        v=1._dl/sqrt(1._dl+aq*aq)
        vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
        pinudot=pinudot+nu_int_kernel(iq)*(ydot(ind)*v+y(ind)*vdot)
        ind=ind+EV%lmaxnu_tau(nu_i)+1
    end do
    ind = EV%nu_pert_ix+2
    do iq=EV%nq(nu_i)+1,nqmax
        q=nu_q(iq)
        aq=am/q
        aqdot=aq*adotoa
        pert_scale=(nu_masses(nu_i)/q)**2/2
        v=1._dl/sqrt(1._dl+aq*aq)
        vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
        psi2dot=ydot(EV%r_ix+2)  + pert_scale*ydot(ind)
        psi2=y(EV%r_ix+2)  + pert_scale*y(ind)
        pinudot=pinudot+nu_int_kernel(iq)*(psi2dot*v+psi2*vdot)
    end do

    end subroutine Nu_pinudot

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function Nu_pi(EV, y, a, nu_i) result(pinu)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvart)
    real(dl) :: am
    real(dl) pinu,q,aq,v
    integer iq, ind

    if (EV%nq(nu_i)/=nqmax) call MpiStop('Nu_pi: nq/=nqmax')
    pinu=0
    ind=EV%nu_ix(nu_i)+2
    am=a*nu_masses(nu_i)
    do iq=1, EV%nq(nu_i)
        q=nu_q(iq)
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        pinu=pinu+nu_int_kernel(iq)*y(ind)*v
        ind =ind+EV%lmaxnut+1
    end do

    end function Nu_pi

    !cccccccccccccccccccccccccccccccccccccccccccccc
    subroutine Nu_Intvsq(EV,y, a, nu_i, G11,G30)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvar)
    real(dl), intent(OUT) ::  G11,G30

    !  Compute the third order variables (in velocity dispersion)
    !by integrating over momentum.
    real(dl) aq,q,v, am
    integer iq, ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    am=a*nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    G11=0._dl
    G30=0._dl
    if (EV%nq(nu_i)/=nqmax) call MpiStop('Nu_Intvsq nq/=nqmax0')
    do iq=1, EV%nq(nu_i)
        q=nu_q(iq)
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        G11=G11+nu_int_kernel(iq)*y(ind+1)*v**2
        if (EV%lmaxnu_tau(nu_i)>2) then
            G30=G30+nu_int_kernel(iq)*y(ind+3)*v**2
        end if
        ind = ind+EV%lmaxnu_tau(nu_i)+1
    end do

    end subroutine Nu_Intvsq


    subroutine MassiveNuVars(EV,y,a,grho,gpres,dgrho,dgq, wnu_arr)
    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar), a, grho,gpres,dgrho,dgq
    real(dl), intent(out), optional :: wnu_arr(max_nu)
    !grho = a^2 kappa rho
    !gpres = a^2 kappa p
    !dgrho = a^2 kappa \delta\rho
    !dgp =  a^2 kappa \delta p
    !dgq = a^2 kappa q (heat flux)
    integer nu_i
    real(dl) grhormass_t, rhonu, qnu, clxnu, grhonu_t, gpnu_t, pnu

    do nu_i = 1, CP%Nu_mass_eigenstates
        grhormass_t=grhormass(nu_i)/a**2

        !Get density and pressure as ratio to massless by interpolation from table
        call Nu_background(a*nu_masses(nu_i),rhonu,pnu)

        if (EV%MassiveNuApprox(nu_i)) then
            clxnu=y(EV%nu_ix(nu_i))
            qnu=y(EV%nu_ix(nu_i)+2)
        else
            !Integrate over q
            call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu)
            !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
            qnu=qnu/rhonu
            clxnu = clxnu/rhonu
        endif

        grhonu_t=grhormass_t*rhonu
        gpnu_t=grhormass_t*pnu

        grho = grho  + grhonu_t
        gpres= gpres + gpnu_t
        dgrho= dgrho + grhonu_t*clxnu
        dgq  = dgq   + grhonu_t*qnu

        if (present(wnu_arr)) then
            wnu_arr(nu_i) =pnu/rhonu
        end if
    end do

    end subroutine MassiveNuVars

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine output(EV,y, tau,sources, num_custom_sources)
    use ThermoData
    type(EvolutionVars) EV
    real(dl), target :: y(EV%nvar), yprime(EV%nvar)
    real(dl) tau
    real(dl), target :: sources(CTransScal%NumSources)
    integer, intent(in) :: num_custom_sources

    yprime = 0
    EV%OutputSources => Sources
    if (num_custom_sources>0) &
        EV%CustomSources => sources(CTransScal%NumSources - num_custom_sources+1:)
    call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)
    nullify(EV%OutputSources, EV%CustomSources)

    end subroutine output



    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine outputt(EV,yt,n,tau,dt,dte,dtb)
    !calculate the tensor sources for open and closed case
    use ThermoData

    implicit none
    integer n
    type(EvolutionVars) :: EV
    real(dl), target :: yt(n), ytprime(n)
    real(dl) tau,dt,dte,dtb,x,polterdot,polterddot,prefac
    real(dl) pig, pigdot, octg, aux, polter, shear, adotoa,a
    real(dl) sinhxr,cothxor
    real(dl) k,k2
    real(dl), dimension(:),pointer :: E,Bprime,Eprime
    real(dl), target :: pol(3),polEprime(3), polBprime(3)
    real(dl) dtauda
    real(dl) opacity, dopacity, ddopacity, &
        visibility, dvisibility, ddvisibility, exptau, lenswindow

    call derivst(EV,EV%nvart,tau,yt,ytprime)

    k2=EV%k2_buf
    k=EV%k_buf
    aux=EV%aux_buf
    shear = yt(3)

    x=(CP%tau0-tau)/CP%r
    call IonizationFunctionsAtTime(tau, opacity, dopacity, ddopacity, &
        visibility, dvisibility, ddvisibility, exptau, lenswindow)

    !  And the electric part of the Weyl.
    if (.not. EV%TensTightCoupling) then
        !  Use the full expression for pigdt
        pig=yt(EV%g_ix+2)
        pigdot=ytprime(EV%g_ix+2)
        E => yt(EV%E_ix+1:)
        Eprime=> ytprime(EV%E_ix+1:)
        Bprime => ytprime(EV%B_ix+1:)
        octg=ytprime(EV%g_ix+3)
    else
        !  Use the tight-coupling approximation
        a =yt(1)
        adotoa = 1/(a*dtauda(a))
        pigdot=32._dl/45._dl*k/opacity*(2._dl*adotoa*shear+ytprime(3))
        pig = 32._dl/45._dl*k/opacity*shear
        pol=0
        polEprime=0
        polBprime=0
        E=>pol
        EPrime=>polEPrime
        BPrime=>polBPrime
        E(2)=pig/4._dl
        EPrime(2)=pigdot/4
        octg=0
    endif

    sinhxr=rofChi(x)*CP%r

    if (EV%q*sinhxr > 1.e-8_dl) then
        prefac=sqrt(EV%q2*CP%r*CP%r-CP%Ksign)
        cothxor=cosfunc(x)/sinhxr

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
        polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*pigdot
        polterddot = 9._dl/15._dl*(-dopacity*(E(2)-polter)-opacity*(  &
            Eprime(2)-polterdot) + k*(2._dl/3._dl*Bprime(2)*aux - 5._dl/27._dl*Eprime(3)*EV%Kft(2))) &
            +0.1_dl*(k*(-octg*EV%Kft(2)/3._dl + 8._dl/15._dl*ytprime(3)) - &
            dopacity*(pig - polter) - opacity*(pigdot-polterdot))

        dt=(shear*exptau + (15._dl/8._dl)*polter*visibility/k)*CP%r/sinhxr**2/prefac

        dte=CP%r*15._dl/8._dl/k/prefac* &
            ((ddvisibility*polter + 2._dl*dvisibility*polterdot + visibility*polterddot)  &
            + 4._dl*cothxor*(dvisibility*polter + visibility*polterdot) - &
            visibility*polter*(k2 -6*cothxor**2))

        dtb=15._dl/4._dl*EV%q*CP%r/k/prefac*(visibility*(2._dl*cothxor*polter + polterdot) + dvisibility*polter)
    else
        dt=0._dl
        dte=0._dl
        dtb=0._dl
    end if

    end subroutine outputt

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine outputv(EV,yv,n,tau,dt,dte,dtb)
    !calculate the vector sources
    use ThermoData

    implicit none
    integer n
    type(EvolutionVars) :: EV
    real(dl), target :: yv(n), yvprime(n)
    real(dl) tau,dt,dte,dtb,x,polterdot
    real(dl) vb,qg, pig, polter, sigma
    real(dl) k,k2
    real(dl), dimension(:),pointer :: E,Eprime
    real(dl) opacity, dopacity, ddopacity, &
        visibility, dvisibility, ddvisibility, exptau, lenswindow


    call derivsv(EV,EV%nvarv,tau,yv,yvprime)

    k2=EV%k2_buf
    k=EV%k_buf
    sigma = yv(2)
    vb  = yv(3)
    qg  = yv(4)
    pig = yv(5)


    x=(CP%tau0-tau)*k

    if (x > 1.e-8_dl) then
        E => yv(EV%lmaxv+3:)
        Eprime=> yvprime(EV%lmaxv+3:)

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)
        polterdot=9._dl/15._dl*Eprime(2) + 0.1_dl*yvprime(5)

        call IonizationFunctionsAtTime(tau, opacity, dopacity, ddopacity, &
            visibility, dvisibility, ddvisibility, exptau, lenswindow)

        if (yv(1) < 1e-3) then
            dt = 1
        else
            dt =0
        end if
        dt= (4*(vb+sigma)*visibility + 15._dl/2/k*( visibility*polterdot + dvisibility*polter) &
            + 4*(exptau*yvprime(2)) )/x

        dte= 15._dl/2*2*polter/x**2*visibility + 15._dl/2/k*(dvisibility*polter + visibility*polterdot)/x

        dtb= -15._dl/2*polter/x*visibility
    else
        dt=0
        dte=0
        dtb=0
    end if

    end subroutine outputv


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine initial(EV,y, tau)
    !  Initial conditions.
    use ThermoData
    implicit none

    type(EvolutionVars) EV
    real(dl) y(EV%nvar)
    real(dl) Rp15,tau,x,x2,x3,om,omtau, &
        Rc,Rb,Rv,Rg,grhonu,chi
    real(dl) k,k2
    real(dl) a,a2, iqg, rhomass,a_massive, ep
    integer l,i, nu_i, j, ind
    integer, parameter :: i_clxg=1,i_clxr=2,i_clxc=3, i_clxb=4, &
        i_qg=5,i_qr=6,i_vb=7,i_pir=8, i_eta=9, i_aj3r=10,i_clxde=11,i_vde=12
    integer, parameter :: i_max = i_vde
    real(dl) initv(6,1:i_max), initvec(1:i_max)

    nullify(EV%OutputTransfer) !Should not be needed, but avoids issues in ifort 14
    nullify(EV%OutputSources)
    nullify(EV%CustomSources)

    if (CP%flat) then
        EV%k_buf=EV%q
        EV%k2_buf=EV%q2
        EV%Kf(1:EV%MaxlNeeded)=1._dl
    else
        EV%k2_buf=EV%q2-CP%curv
        EV%k_buf=sqrt(EV%k2_buf)

        do l=1,EV%MaxlNeeded
            EV%Kf(l)=1._dl-CP%curv*(l*(l+2))/EV%k2_buf
        end do
    end if

    k=EV%k_buf
    k2=EV%k2_buf

    do j=1,EV%MaxlNeeded
        EV%denlk(j)=denl(j)*k*j
        EV%denlk2(j)=denl(j)*k*EV%Kf(j)*(j+1)
        EV%polfack(j)=polfac(j)*k*EV%Kf(j)*denl(j)
    end do

    !Get time to switch off tight coupling
    !The numbers here are a bit of guesswork
    !The high k increase saves time for very small loss of accuracy
    !The lower k ones are more delicate. Nead to avoid instabilities at same time
    !as ensuring tight coupling is accurate enough
    if (EV%k_buf > epsw) then
        if (EV%k_buf > epsw*5) then
            ep=ep0*5/AccuracyBoost
            if (HighAccuracyDefault) ep = ep*0.65
        else
            ep=ep0
        end if
    else
        ep=ep0
    end if
    if (second_order_tightcoupling) ep=ep*2
    EV%TightSwitchoffTime = min(tight_tau,Thermo_OpacityToTime(EV%k_buf/ep))


    y=0

    !  k*tau, (k*tau)**2, (k*tau)**3
    x=k*tau
    x2=x*x
    x3=x2*x
    rhomass =  sum(grhormass(1:CP%Nu_mass_eigenstates))
    grhonu=rhomass+grhornomass

    om = (grhob+grhoc)/sqrt(3*(grhog+grhonu))
    omtau=om*tau
    Rv=grhonu/(grhonu+grhog)

    Rg = 1-Rv
    Rc=CP%omegac/(CP%omegac+CP%omegab)
    Rb=1-Rc
    Rp15=4*Rv+15

    if (CP%Scalar_initial_condition > initial_nummodes) &
        call MpiStop('Invalid initial condition for scalar modes')

    a=tau*adotrad*(1+omtau/4)
    a2=a*a

    initv=0

    !  Set adiabatic initial conditions

    chi=1  !Get transfer function for chi
    initv(1,i_clxg)=-chi*EV%Kf(1)/3*x2*(1-omtau/5)
    initv(1,i_clxr)= initv(1,i_clxg)
    initv(1,i_clxb)=0.75_dl*initv(1,i_clxg)
    initv(1,i_clxc)=initv(1,i_clxb)
    initv(1,i_qg)=initv(1,i_clxg)*x/9._dl
    initv(1,i_qr)=-chi*EV%Kf(1)*(4*Rv+23)/Rp15*x3/27
    initv(1,i_vb)=0.75_dl*initv(1,i_qg)
    initv(1,i_pir)=chi*4._dl/3*x2/Rp15*(1+omtau/4*(4*Rv-5)/(2*Rv+15))
    initv(1,i_aj3r)=chi*4/21._dl/Rp15*x3
    initv(1,i_eta)=-chi*2*EV%Kf(1)*(1 - x2/12*(-10._dl/Rp15 + EV%Kf(1)))

    if (CP%Scalar_initial_condition/= initial_adiabatic) then
        !CDM isocurvature

        initv(2,i_clxg)= Rc*omtau*(-2._dl/3 + omtau/4)
        initv(2,i_clxr)=initv(2,i_clxg)
        initv(2,i_clxb)=initv(2,i_clxg)*0.75_dl
        initv(2,i_clxc)=1+initv(2,i_clxb)
        initv(2,i_qg)=-Rc/9*omtau*x
        initv(2,i_qr)=initv(2,i_qg)
        initv(2,i_vb)=0.75_dl*initv(2,i_qg)
        initv(2,i_pir)=-Rc*omtau*x2/3/(2*Rv+15._dl)
        initv(2,i_eta)= Rc*omtau*(1._dl/3 - omtau/8)*EV%Kf(1)
        initv(2,i_aj3r)=0
        !Baryon isocurvature
        if (Rc==0) call MpiStop('Isocurvature initial conditions assume non-zero dark matter')

        initv(3,:) = initv(2,:)*(Rb/Rc)
        initv(3,i_clxc) = initv(3,i_clxb)
        initv(3,i_clxb)= initv(3,i_clxb)+1

        !neutrino isocurvature density mode

        initv(4,i_clxg)=Rv/Rg*(-1 + x2/6)
        initv(4,i_clxr)=1-x2/6
        initv(4,i_clxc)=-omtau*x2/80*Rv*Rb/Rg
        initv(4,i_clxb)= Rv/Rg/8*x2
        iqg = - Rv/Rg*(x/3 - Rb/4/Rg*omtau*x)
        initv(4,i_qg) =iqg
        initv(4,i_qr) = x/3
        initv(4,i_vb)=0.75_dl*iqg
        initv(4,i_pir)=x2/Rp15
        initv(4,i_eta)=EV%Kf(1)*Rv/Rp15/3*x2

        !neutrino isocurvature velocity mode

        initv(5,i_clxg)=Rv/Rg*x - 2*x*omtau/16*Rb*(2+Rg)/Rg**2
        initv(5,i_clxr)=-x -3*x*omtau*Rb/16/Rg
        initv(5,i_clxc)=-9*omtau*x/64*Rv*Rb/Rg
        initv(5,i_clxb)= 3*Rv/4/Rg*x - 9*omtau*x/64*Rb*(2+Rg)/Rg**2
        iqg = Rv/Rg*(-1 + 3*Rb/4/Rg*omtau+x2/6 +3*omtau**2/16*Rb/Rg**2*(Rg-3*Rb))
        initv(5,i_qg) =iqg
        initv(5,i_qr) = 1 - x2/6*(1+4*EV%Kf(1)/(4*Rv+5))
        initv(5,i_vb)=0.75_dl*iqg
        initv(5,i_pir)=2*x/(4*Rv+5)+omtau*x*6/Rp15/(4*Rv+5)
        initv(5,i_eta)=2*EV%Kf(1)*x*Rv/(4*Rv+5) + omtau*x*3*EV%Kf(1)*Rv/32*(Rb/Rg - 80/Rp15/(4*Rv+5))
        initv(5,i_aj3r) = 3._dl/7*x2/(4*Rv+5)

        !quintessence isocurvature mode
    end if

    if (CP%Scalar_initial_condition==initial_vector) then
        InitVec = 0
        do i=1,initial_nummodes
            InitVec = InitVec+ initv(i,:)*CP%InitialConditionVector(i)
        end do
    else
        InitVec = initv(CP%Scalar_initial_condition,:)
        if (CP%Scalar_initial_condition==initial_adiabatic) InitVec = -InitVec
        !So we start with chi=-1 as before
    end if

    y(1)=a
    y(2)= -InitVec(i_eta)*k/2
    !get eta_s*k, where eta_s is synchronous gauge variable

    !  CDM
    y(3)=InitVec(i_clxc)

    !  Baryons
    y(4)=InitVec(i_clxb)
    y(5)=InitVec(i_vb)

    !  Photons
    y(EV%g_ix)=InitVec(i_clxg)
    y(EV%g_ix+1)=InitVec(i_qg)

    if (w_lam /= -1 .and. w_Perturb) then
        y(EV%w_ix) = InitVec(i_clxde)
        y(EV%w_ix+1) = InitVec(i_vde)
    end if

    !  Neutrinos
    y(EV%r_ix)=InitVec(i_clxr)
    y(EV%r_ix+1)=InitVec(i_qr)
    y(EV%r_ix+2)=InitVec(i_pir)

    if (EV%lmaxnr>2) then
        y(EV%r_ix+3)=InitVec(i_aj3r)
    endif

    if (CP%Num_Nu_massive == 0) return

    do nu_i = 1, CP%Nu_mass_eigenstates
        EV%MassiveNuApproxTime(nu_i) = Nu_tau_massive(nu_i)
        a_massive =  20000*k/nu_masses(nu_i)*AccuracyBoost*lAccuracyBoost
        if (a_massive >=0.99) then
            EV%MassiveNuApproxTime(nu_i)=CP%tau0+1
        else if (a_massive > 17.d0/nu_masses(nu_i)*AccuracyBoost) then
            EV%MassiveNuApproxTime(nu_i)=max(EV%MassiveNuApproxTime(nu_i),DeltaTime(0._dl,a_massive, 0.01_dl))
        end if
        ind = EV%nu_ix(nu_i)
        do  i=1,EV%nq(nu_i)
            y(ind:ind+2)=y(EV%r_ix:EV%r_ix+2)
            if (EV%lmaxnu_tau(nu_i)>2) y(ind+3)=InitVec(i_aj3r)
            ind = ind + EV%lmaxnu_tau(nu_i)+1
        end do
    end do

    end subroutine initial


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine initialt(EV,yt,tau)
    !  Initial conditions for tensors
    use ThermoData
    implicit none
    real(dl) bigR,tau,x,aj3r,elec, pir, rhomass
    integer l
    type(EvolutionVars) EV
    real(dl) k,k2 ,a, omtau
    real(dl) yt(EV%nvart)
    real(dl) tens0, ep, tensfac

    if (CP%flat) then
        EV%aux_buf=1._dl
        EV%k2_buf=EV%q2
        EV%k_buf=EV%q
        EV%Kft(1:EV%MaxlNeededt)=1._dl !initialize for flat case
    else
        EV%k2_buf=EV%q2-3*CP%curv
        EV%k_buf=sqrt(EV%k2_buf)
        EV%aux_buf=sqrt(1._dl+3*CP%curv/EV%k2_buf)
    endif

    k=EV%k_buf
    k2=EV%k2_buf

    do l=1,EV%MaxlNeededt
        if (.not. CP%flat) EV%Kft(l)=1._dl-CP%curv*((l+1)**2-3)/k2
        EV%denlkt(1,l)=k*denl(l)*l !term for L-1
        tensfac=real((l+3)*(l-1),dl)/(l+1)
        EV%denlkt(2,l)=k*denl(l)*tensfac*EV%Kft(l) !term for L+1
        EV%denlkt(3,l)=k*denl(l)*tensfac**2/(l+1)*EV%Kft(l) !term for polarization
        EV%denlkt(4,l)=k*4._dl/(l*(l+1))*EV%aux_buf !other for polarization
    end do

    if (k > 0.06_dl*epsw) then
        ep=ep0
    else
        ep=0.2_dl*ep0
    end if

    !    finished_tightcoupling = ((k/opacity > ep).or.(1._dl/(opacity*tau) > ep))
    EV%TightSwitchoffTime = min(tight_tau,Thermo_OpacityToTime(EV%k_buf/ep))

    a=tau*adotrad
    rhomass =  sum(grhormass(1:CP%Nu_mass_eigenstates))
    omtau = tau*(grhob+grhoc)/sqrt(3*(grhog+rhomass+grhornomass))

    if (DoTensorNeutrinos) then
        bigR = (rhomass+grhornomass)/(rhomass+grhornomass+grhog)
    else
        bigR = 0._dl
    end if

    x=k*tau

    yt(1)=a
    tens0 = 1

    yt(2)= tens0
    !commented things are for the compensated mode with magnetic fields; can be neglected
    !-15/28._dl*x**2*(bigR-1)/(15+4*bigR)*Magnetic*(1-5./2*omtau/(2*bigR+15))

    elec=-tens0*(1+2*CP%curv/k2)*(2*bigR+10)/(4*bigR+15) !elec, with H=1

    !shear
    yt(3)=-5._dl/2/(bigR+5)*x*elec
    !          + 15._dl/14*x*(bigR-1)/(4*bigR+15)*Magnetic*(1 - 15./2*omtau/(2*bigR+15))

    yt(4:EV%nvart)=0._dl

    !  Neutrinos
    if (DoTensorNeutrinos) then
        pir=-2._dl/3._dl/(bigR+5)*x**2*elec
        !           + (bigR-1)/bigR*Magnetic*(1-15./14*x**2/(15+4*bigR))
        aj3r=  -2._dl/21._dl/(bigR+5)*x**3*elec !&
            !           + 3._dl/7*x*(bigR-1)/bigR*Magnetic
            yt(EV%r_ix+2)=pir
        yt(EV%r_ix+3)=aj3r
        !Should set up massive too, but small anyway..
    end if

    end subroutine initialt

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine initialv(EV,yv,tau)
    !  Initial conditions for vectors

    implicit none
    real(dl) bigR,Rc,tau,x,pir
    type(EvolutionVars) EV
    real(dl) k,k2 ,a, omtau
    real(dl) yv(EV%nvarv)

    if (CP%flat) then
        EV%k2_buf=EV%q2
        EV%k_buf=EV%q
    else
        call MpiStop('Vectors not supported in non-flat models')
    endif

    k=EV%k_buf
    k2=EV%k2_buf

    omtau = tau*(grhob+grhoc)/sqrt(3*(grhog+grhornomass))

    a=tau*adotrad*(1+omtau/4)

    x=k*tau

    bigR = (grhornomass)/(grhornomass+grhog)
    Rc=CP%omegac/(CP%omegac+CP%omegab)

    yv(1)=a


    yv(2)= vec_sig0*(1- 15._dl/2*omtau/(4*bigR+15)) + 45._dl/14*x*Magnetic*(BigR-1)/(4*BigR+15)
    !qg
    yv(4)= vec_sig0/3* (4*bigR + 5)/(1-BigR)*(1  -0.75_dl*omtau*(Rc-1)/(bigR-1)* &
        (1 - 0.25_dl*omtau*(3*Rc-2-bigR)/(BigR-1))) &
        -x/2*Magnetic
    yv(3)= 3._dl/4*yv(4)

    yv(5:EV%nvarv) = 0

    !        if (.false.) then
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1) = vec_sig0/6/bigR*x**2*(1+2*bigR*omtau/(4*bigR+15))
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+2) = -2/3._dl*vec_sig0/bigR*x*(1 +3*omtau*bigR/(4*bigR+15))
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+3) = 1/4._dl*vec_sig0/bigR*(5+4*BigR)
    !         yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+4) =1/9.*x*vec_sig0*(5+4*bigR)/bigR
    !         yv(4) = 0
    !         yv(3)= 3._dl/4*yv(4)
    !          return
    !        end if

    !  Neutrinos
    !q_r
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1) = -1._dl/3*vec_sig0*(4*BigR+5)/bigR &
        + x**2*vec_sig0/6/BigR +0.5_dl*x*(1/bigR-1)*Magnetic
    !pi_r
    pir=-2._dl/3._dl*x*vec_sig0/BigR - (1/bigR-1)*Magnetic
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1 +1)=pir
    yv((EV%lmaxv-1+1)+(EV%lmaxpolv-1)*2+3+1 +2)=3._dl/7*x*Magnetic*(1-1/BigR)

    end subroutine initialv


    subroutine outtransf(EV, y,tau, Arr)
    !write out clxc, clxb, clxg, clxn
    implicit none
    type(EvolutionVars) EV
    real(dl), intent(in) :: tau
    real, target :: Arr(:)
    real(dl) y(EV%nvar),yprime(EV%nvar)

    yprime = 0
    EV%OutputTransfer =>  Arr
    call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)
    nullify(EV%OutputTransfer)

    Arr(Transfer_kh+1:Transfer_max) = Arr(Transfer_kh+1:Transfer_max)/EV%k2_buf

    end subroutine outtransf

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine derivs(EV,n,tau,ay,ayprime)
    !  Evaluate the time derivatives of the perturbations
    !  ayprime is not necessarily GaugeInterface.yprime, so keep them distinct
    use ThermoData
    use MassiveNu
    implicit none
    type(EvolutionVars) EV

    integer n,nu_i
    real(dl) ay(n),ayprime(n)
    real(dl) tau,w
    real(dl) k,k2

    !  Internal variables.

    real(dl) opacity
    real(dl) photbar,cs2,pb43,grho,slip,clxgdot, &
        clxcdot,clxbdot,adotdota,gpres,clxrdot,etak
    real(dl) q,aq,v
    real(dl) G11_t,G30_t, wnu_arr(max_nu)

    real(dl) dgq,grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t,sigma,polter
    real(dl) qgdot,qrdot,pigdot,pirdot,vbdot,dgrho,adotoa
    real(dl) a,a2,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir
    real(dl) clxde, qde,  E2, dopacity
    integer l,i,ind, ind2, off_ix, ix
    real(dl) dgs,sigmadot,dz !, ddz
    real(dl) dgpi,dgrho_matter,grho_matter, clxnu, gpres_nu
    !non-flat vars
    real(dl) cothxor !1/tau in flat case
    !Variables for source calculation
    real(dl) diff_rhopi, pidot_sum, dgpi_diff, phi
    real(dl) E(2:3), Edot(2:3)
    real(dl) phidot, polterdot, polterddot, octg, octgdot
    real(dl) ddopacity, visibility, dvisibility, ddvisibility, exptau, lenswindow
    real(dl) ISW, quadrupole_source, doppler, monopole_source, tau0
	!>ISiTGR MOD START: adding new variables
    !JD ISiTGR vars
    real(dl) TGR_f_q, TGR_f_1, TGR_D_T, TGR_D_T_dot,TGR_Q_T,TGR_Q_T_dot
    real(dl) etakdot, Hdot, TGR_rhoDeltadot !pidot_sum
	real(dl) phipluspsi
	!CGQ new variables ---------------------------------
	real(dl) TGR_mu, TGR_mudot, TGR_eta, TGR_etadot, TGR_Sigma, TGR_Sigmadot !for MG functions
	real(dl) TGR_f_mueta, TGR_f_muSigma, TGR_rhoDelta, dgpi_3wplus1, dgpi_3wplus2, dgpi_3wplus1plusbetak !for terms to compute etakdot
	real(dl) TGR_Phi, TGR_Phidot, TGR_Psi, TGR_Psidot !MG source functions in Newtonian gauge
	real(dl) betak, gpresv_t !CGQ for spatial curvature and Dark Energy parametrizations
    !CGQ -----------------------------------------------
    !<ISiTGR MOD END
	
	k=EV%k_buf
    k2=EV%k2_buf
   
   	betak = 1.d0/EV%Kf(1) !betak=1 for flat universe
   
	a=ay(1)
    a2=a*a
	
    etak=ay(2)
	
    !  CDM variables
    clxc=ay(3)

    !  Baryon variables
    clxb=ay(4)
    vb=ay(5)

    !  Compute expansion rate from: grho 8*pi*rho*a**2

    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2
	
	!>ISiTGR MOD START
	if (TGR%GR==0) then
		if (w_lam==-1._dl) then
    	    grhov_t=grhov*a2
    	else
    	    grhov_t=grhov*a**(-1-3*w_lam)
    	end if
	else
	
		if ( TGR%DE_eqstate == 0 ) then
			grhov_t = grhov*a**2.d0

		else if ( TGR%DE_eqstate == 1 ) then
			grhov_t = grhov*a**(-1.d0-3.d0*TGR%w0)

		else if ( TGR%DE_eqstate == 2 ) then
			grhov_t = grhov*a**(-1.d0-3.d0*(TGR%w0+TGR%wa))*exp(3.d0*TGR%wa*(a-1.d0))

		else if ( TGR%DE_eqstate == 3 ) then
			grhov_t = grhov*a**(-1.d0-3.d0*(TGR%wp+TGR%a_p*TGR%wa))*exp(3.d0*TGR%wa*(a-1.d0))
		
		end if
	end if
	!<ISiTGR MOD END

    !  Get sound speed and ionisation fraction.
    if (EV%TightCoupling) then
        call thermo(tau,cs2,opacity,dopacity)
    else
        call thermo(tau,cs2,opacity)
    end if

    gpres_nu=0
    grhonu_t=0

    !total perturbations: matter terms first, then add massive nu, de and radiation
    !  8*pi*a*a*SUM[rho_i*clx_i]
    dgrho_matter=grhob_t*clxb+grhoc_t*clxc
    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=grhob_t*vb

    if (CP%Num_Nu_Massive > 0) then
        call MassiveNuVars(EV,ay,a,grhonu_t,gpres_nu,dgrho_matter,dgq, wnu_arr)
    end if

    grho_matter=grhonu_t+grhob_t+grhoc_t
    grho = grho_matter+grhor_t+grhog_t+grhov_t

	if (CP%flat) then
        adotoa=sqrt(grho/3)
        cothxor=1._dl/tau
    else
        adotoa=sqrt((grho+grhok)/3._dl)
        cothxor=1._dl/tanfunc(tau/CP%r)/CP%r
    end if	
	
    dgrho = dgrho_matter

	if (TGR%GR==0) then
	    if (w_lam /= -1 .and. w_Perturb) then
	        clxde=ay(EV%w_ix)
	        qde=ay(EV%w_ix+1)*(1+w_lam)
	        dgrho=dgrho + clxde*grhov_t
	        dgq = dgq + qde*grhov_t
	    end if
	end if
		!>ISiTGR MOD START 
    if (EV%no_nu_multpoles) then
		if (TGR%GR==0) then !CGQ to work with default GR or with MG models
        !RSA approximation of arXiv:1104.2933, dropping opactity terms in the velocity
        !Approximate total density variables with just matter terms
        	z=(0.5_dl*dgrho/k + etak)/adotoa
        	dz= -adotoa*z - 0.5_dl*dgrho/k
        	clxr=-4*dz/k
        	qr=-4._dl/3*z
        	pir=0
		else
	    	clxr=2*(grhoc_t*clxc+grhob_t*clxb)/3/k**2
        	qr= clxr*k/sqrt((grhoc_t+grhob_t)/3)*(2/3._dl)
        	pir=0
		end if
		!<ISiTGR MOD END
    else
        !  Massless neutrinos
        clxr=ay(EV%r_ix)
        qr  =ay(EV%r_ix+1)
        pir =ay(EV%r_ix+2)
    endif
	!>ISiTGR MOD START
    if (EV%no_phot_multpoles) then
		if (TGR%GR==0) then !CGQ to work with default GR or with MG models
    	    if (.not. EV%no_nu_multpoles) then
	            z=(0.5_dl*dgrho/k + etak)/adotoa
            	dz= -adotoa*z - 0.5_dl*dgrho/k
        	    clxg=-4*dz/k-4/k*opacity*(vb+z)
    	        qg=-4._dl/3*z
	        else
            	clxg=clxr-4/k*opacity*(vb+z)
        	    qg=qr
    	    end if
	        pig=0
		else
			clxg=2*(grhoc_t*clxc+grhob_t*clxb)/3/k**2
    	    qg= clxg*k/sqrt((grhoc_t+grhob_t)/3)*(2/3._dl)
        	pig=0
		end if
	else
        !  Photons
        clxg=ay(EV%g_ix)
        qg=ay(EV%g_ix+1)
        if (.not. EV%TightCoupling) pig=ay(EV%g_ix+2)
    end if
	!<ISiTGR MOD END

    !  8*pi*a*a*SUM[rho_i*clx_i] - radiation terms
    dgrho=dgrho + grhog_t*clxg+grhor_t*clxr

    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    dgq=dgq + grhog_t*qg+grhor_t*qr

    !  Photon mass density over baryon mass density
    photbar=grhog_t/grhob_t
    pb43=4._dl/3*photbar

    ayprime(1)=adotoa*a
	
	!>ISiTGR MOD START-----------------------------------------------------------------
    !all this module was written by CGQ, adding mueta and muSigma parameterizations.
    !Also, CGQ modified Q,D parameterization originally written by JD in order to work with massive neutrinos.
    
!CGQ expanded ISiTGR block for new parameterizations and massive neutrinos. 

! 1) Get sigma
    if (TGR%GR/=0) then
		!CGQ for Dark Energy pressure ----------------------
		if ( TGR%DE_eqstate == 0 ) then
			gpresv_t = -grhov_t

		else if ( TGR%DE_eqstate == 1 ) then
			gpresv_t = TGR%w0*grhov_t

		else if ( TGR%DE_eqstate == 2 ) then
			gpresv_t = (TGR%w0+(1.d0-a)*TGR%wa)*grhov_t
		
		else if ( TGR%DE_eqstate == 3 ) then
			gpresv_t = (TGR%wp+(TGR%a_p-a)*TGR%wa)*grhov_t

		end if	
		!---------------------------------------------------
    gpres = (grhog_t+grhor_t)/3.d0 + gpresv_t + gpres_nu !CGQ
    TGR_rhoDelta = dgrho+3.d0*adotoa*dgq/k !CGQ
   	Hdot = -(grho+3.d0*gpres)/6.d0   !JD
	dgpi = 0.d0
	dgpi_3wplus1 = 0.d0
	dgpi_3wplus2 = 0.d0
	dgpi_3wplus1plusbetak = 0.d0
    
		if ((TGR%ISiTGR_mueta .eqv. .true.) .or. (TGR%ISiTGR_BIN_mueta .eqv. .true.)) then

			!mu-eta parameterization
		    TGR_mu = ISiTGR_mu(k,a,adotoa)
    		TGR_mudot = ISiTGR_mu_dot(k,a,adotoa,Hdot)
    		TGR_eta = ISiTGR_eta(k,a,adotoa)
	    	TGR_etadot = ISiTGR_eta_dot(k,a,adotoa,Hdot)

			!here we get the contributions for massive neutrinos to dgpi and dgpi_3wplus2
			if (CP%Num_Nu_Massive /= 0) then                                                            
				call MassiveNuVarsOut(EV, ay, ayprime, a, dgpi=dgpi, dgpi_3wplus1=dgpi_3wplus1, dgpi_3wplus2=dgpi_3wplus2, &
				dgpi_3wplus1plusbetak=dgpi_3wplus1plusbetak) 
    		end if             
            
            !adding contributions from photons and massless neutrinos
			dgpi = dgpi + grhor_t*pir + grhog_t*pig 
			dgpi_3wplus1 = dgpi_3wplus1 + 2.d0*(grhor_t*pir+grhog_t*pig) !!Note that (3w_rad+1)=(1+1)=2
			dgpi_3wplus2 = dgpi_3wplus2 + 3.d0*(grhor_t*pir+grhog_t*pig) !!Note that (3w_rad+2)=(1+2)=3
			dgpi_3wplus1plusbetak = dgpi_3wplus1plusbetak + (2.d0+betak)*(grhor_t*pir+grhog_t*pig) !!Note that (3w_rad+1+betak)=2+betak
            
			!!getting the newtonian potentials for MG
			TGR_Psi = -0.5d0*TGR_mu*(betak*TGR_rhoDelta+2.d0*dgpi)/k2
			TGR_Phi = (TGR_mu/k2)*dgpi+TGR_eta*TGR_Psi
	
			!Computing sigma_camb, just called sigma here
			sigma = etak/adotoa - k*TGR_Phi/adotoa

		else if ((TGR%ISiTGR_muSigma .eqv. .true.) .or. (TGR%ISiTGR_BIN_muSigma .eqv. .true.)) then
			!mu-Sigma parameterization
   	 		TGR_mu = ISiTGR_mu(k,a,adotoa)
		    TGR_mudot = ISiTGR_mu_dot(k,a,adotoa,Hdot)
		    TGR_Sigma = ISiTGR_Sigma(k,a,adotoa)
		    TGR_Sigmadot = ISiTGR_Sigma_dot(k,a,adotoa,Hdot)
	
			!here we get the contributions for massive neutrinos to dgpi and dgpi_w
			if (CP%Num_Nu_Massive /= 0) then
				call MassiveNuVarsOut(EV, ay, ayprime, a, dgpi=dgpi, dgpi_3wplus1=dgpi_3wplus1)
		    end if
                        
            !adding contributions from photons and massless neutrinos
			dgpi = dgpi + grhor_t*pir + grhog_t*pig 
			dgpi_3wplus1 = dgpi_3wplus1 + 2.d0*(grhor_t*pir+grhog_t*pig) !!Note that (3w_rad+1)=(1+1)=2
	
			!!getting the newtonian potentials for MG
			TGR_Psi = -0.5d0*TGR_mu*(betak*TGR_rhoDelta+2.d0*dgpi)/k2
			TGR_Phi = -TGR_Sigma*(betak*TGR_rhoDelta+dgpi)/k2-TGR_Psi
	
			!Computing sigma_camb
			sigma = etak/adotoa - k*TGR_Phi/adotoa
	
		else
	    	!JD Begin ISiTGR block
    		TGR_Q_T=ISiTGR_Q(k,a,adotoa)
		    TGR_Q_T_dot=ISiTGR_Q_dot(k,a,adotoa,Hdot)
		    TGR_D_T= ISiTGR_D(k,a,adotoa)
		    TGR_D_T_dot=ISiTGR_D_dot(k,a,adotoa,Hdot)		

			!here we get the contributions for massive neutrinos to dgpi and dgpi_w
			if (CP%Num_Nu_Massive /= 0) then
				call MassiveNuVarsOut(EV, ay, ayprime, a, dgpi=dgpi, dgpi_3wplus1=dgpi_3wplus1)
 		    end if
    
			!adding contributions from photons and massless neutrinos
			dgpi = dgpi + grhor_t*pir + grhog_t*pig 
			dgpi_3wplus1 = dgpi_3wplus1 + 2.d0*(grhor_t*pir+grhog_t*pig) !!Note that (3w_rad+1)=(1+1)=2
			
			!!getting the newtonian potentials for MG
			TGR_Phi = -TGR_Q_T/(2.d0*k2)*betak*TGR_rhoDelta
			
    		!kalpha
		    sigma=etak/adotoa - k*TGR_Phi/adotoa	
	
		end if
        
	! 2) get contributions of massless neutrinos, photons, and massive neutrinos to pidot_sum
	!Computing pidot_sum by calling early some parts of the derivs subroutine, to get pidot in terms of sigma
    	pidot_sum = 0.d0
	!contribution by massless neutrinos --------------------------------------------------------
		if (EV%no_nu_multpoles) then
            pirdot = 0.d0
    	else
	        if (EV%lmaxnr>2) then
		        ix=EV%r_ix+2
        	    pirdot=EV%denlk(2)*qr- EV%denlk2(2)*ay(ix+1)+8._dl/15._dl*k*sigma
			else
				pirdot=EV%denlk(2)*qr +8._dl/15._dl*k*sigma
        	end if
        end if
		!here we add contribution by massless neutrinos
		pidot_sum = pidot_sum + grhor_t*pirdot

	!contribution by photons -------------------------------------------------------------------
    	if (EV%no_phot_multpoles) then
        	pigdot=0.d0
		else 
			if (EV%tightcoupling) then
				pigdot=0.d0
				pigdot = EV%pigdot
				
				!Use explicit equations for photon moments if appropriate
	        else
    	        E2=ay(EV%polind+2)
        	    polter = pig/10+9._dl/15*E2 !2/15*(3/4 pig + 9/2 E2)
				if (EV%lmaxg>2) then
                    ix= EV%g_ix+2
                    pigdot=EV%denlk(2)*qg-EV%denlk2(2)*ay(ix+1)-opacity*(pig - polter) &
                            +8._dl/15._dl*k*sigma
				else !closed case
					pigdot=EV%denlk(2)*qg-opacity*(pig - polter) +8._dl/15._dl*k*sigma
				end if
        	end if
		end if
		!here we add contribution by photons
		pidot_sum = pidot_sum + grhog_t*pigdot

		!contribution by massive neutrinos --------------------------------------------------------
        if (CP%Num_Nu_massive >0) then
        
            !DIR$ LOOP COUNT MIN(1), AVG(1)
            
            do nu_i = 1, CP%Nu_mass_eigenstates
            
                if (EV%MassiveNuApprox(nu_i)) then
                    !Now EV%iq0 = clx, EV%iq0+1 = clxp, EV%iq0+2 = G_1, EV%iq0+3=G_2=pinu
                    !see astro-ph/0203507
                    G11_t=EV%G11(nu_i)/a/a2
                    G30_t=EV%G30(nu_i)/a/a2
                    off_ix = EV%nu_ix(nu_i)
                    w=wnu_arr(nu_i)
                    !ayprime(off_ix)=-k*z*(w+1) + 3*adotoa*(w*ay(off_ix) - ay(off_ix+1))-k*ay(off_ix+2) !contains z
                    !ayprime(off_ix+1)=(3*w-2)*adotoa*ay(off_ix+1) - 5._dl/3*k*z*w - k/3*G11_t !contains z
                    ayprime(off_ix+2)=(3*w-1)*adotoa*ay(off_ix+2) - k*(2._dl/3*EV%Kf(1)*ay(off_ix+3)-ay(off_ix+1))
                    ayprime(off_ix+3)=(3*w-2)*adotoa*ay(off_ix+3) + 2*w*k*sigma - k/5*(3*EV%Kf(2)*G30_t-2*G11_t)
                    
                else
                    ind=EV%nu_ix(nu_i)
                    !DIR$ LOOP COUNT MIN(3), AVG(3)
                    do i=1,EV%nq(nu_i)
                        q=nu_q(i)
                        aq=a*nu_masses(nu_i)/q
                        v=1._dl/sqrt(1._dl+aq*aq)
                        !ayprime(ind)=-k*(4._dl/3._dl*z + v*ay(ind+1)) !contains z
                        ind=ind+1
                        ayprime(ind)=v*(EV%denlk(1)*ay(ind-1)-EV%denlk2(1)*ay(ind+1))
                        ind=ind+1
                        if (EV%lmaxnu_tau(nu_i)==2) then
                            ayprime(ind)=-ayprime(ind-2) -3*cothxor*ay(ind)
                        else
                            ayprime(ind)=v*(EV%denlk(2)*ay(ind-1)-EV%denlk2(2)*ay(ind+1)) &
                            +k*8._dl/15._dl*sigma
                            do l=3,EV%lmaxnu_tau(nu_i)-1
                                ind=ind+1
                                ayprime(ind)=v*(EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
                            end do
                            !  Truncate moment expansion.
                            ind = ind+1
                            ayprime(ind)=k*v*ay(ind-1)-(EV%lmaxnu_tau(nu_i)+1)*cothxor*ay(ind)
                        end if
                        ind = ind+1
                    end do
                end if
            end do
            
            if (EV%has_nu_relativistic) then
                ind=EV%nu_pert_ix
                ayprime(ind)=+k*a2*qr -k*ay(ind+1)
                ind2= EV%r_ix
                do l=1,EV%lmaxnu_pert-1
                    ind=ind+1
                    ind2=ind2+1
                    ayprime(ind)= -a2*(EV%denlk(l)*ay(ind2-1)-EV%denlk2(l)*ay(ind2+1)) &
                  		  		+ (EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
                end do
                ind=ind+1
                ind2=ind2+1
                ayprime(ind)= k*(ay(ind-1) -a2*ay(ind2-1)) -(EV%lmaxnu_pert+1)*cothxor*ay(ind)
            end if
            
			!here we add contribution by massive neutrinos to pidot_sum
            call MassiveNuVarsOut(EV, ay, ayprime, a, pidot_sum=pidot_sum)
        end if

	! 3) Get sigmadot
		if ((TGR%ISiTGR_mueta .eqv. .true.) .or. (TGR%ISiTGR_muSigma .eqv. .true.) &
			.or. (TGR%ISiTGR_BIN_mueta .eqv. .true.) .or. (TGR%ISiTGR_BIN_muSigma .eqv. .true.)) then

			!Computing derivative of sigma_camb, just called sigmadot. Here, we use alpha=(sigma_camb/k)
			sigmadot = k*TGR_Psi - adotoa*sigma
		else
		
			TGR_Psi = -1.d0/k2*(betak*TGR_D_T*TGR_rhoDelta + TGR_Q_T*dgpi) - TGR_Phi

			!Computing derivative of sigma_camb, just called sigmadot. Here, we use alpha=(sigma_camb/k)
			sigmadot = k*TGR_Psi - adotoa*sigma
		
		end if

	! 4) compute etakdot and z
		if ((TGR%ISiTGR_mueta .eqv. .true.) .or. (TGR%ISiTGR_BIN_mueta .eqv. .true.)) then
			!Computing etakdot
			!w_masslessnu=1/3, w_massivenu=0
			TGR_f_mueta = k2 + 1.5d0*betak*TGR_mu*TGR_eta*((grhoc_t+grhob_t)+(grhor_t+grhog_t)*4.d0/3.d0+(grhonu_t+gpres_nu))!& + (grhov_t+gpresv_t)) 
			
			TGR_f_1 =  1.d0+3.d0*(adotoa**2.d0-Hdot)/k2
			
			etakdot = k/(2.d0*TGR_f_mueta)*(k*TGR_mu*TGR_eta*betak*TGR_f_1*dgq + betak*TGR_rhoDelta*(adotoa*TGR_mu*(TGR_eta-1.d0) - &
				TGR_mudot*TGR_eta-TGR_mu*TGR_etadot) + 2.d0*TGR_mu*(1.d0-TGR_eta)*pidot_sum + k*sigma*(-2.d0*(adotoa**2.d0-Hdot) &
				+TGR_mu*TGR_eta*betak*((grhoc_t+grhob_t)+(grhor_t+grhog_t)*4.d0/3.d0+(grhonu_t+gpres_nu))) - 2.d0*dgpi &
				*(TGR_mu*TGR_etadot+TGR_mudot*(TGR_eta-1.d0)) - 2.d0*adotoa*TGR_mu*dgpi_3wplus2 + 2.d0 * adotoa * TGR_mu * TGR_eta * &
				dgpi_3wplus1plusbetak)

			ayprime(2) = etakdot

			!Computing z_camb
			z = sigma - 3.d0*etakdot/k2

			TGR_Phidot = (etakdot - adotoa*sigmadot - Hdot*sigma)/k

		else if ((TGR%ISiTGR_muSigma .eqv. .true.) .or. (TGR%ISiTGR_BIN_muSigma .eqv. .true.)) then
			!Computing etakdot
			TGR_f_muSigma = k2 + 1.5d0*betak*(2.d0*TGR_Sigma-TGR_mu)*((grhoc_t+grhob_t)+(grhor_t+grhog_t)*4.d0/3.d0+(grhonu_t+gpres_nu))
            
			TGR_f_1 =  1.d0+3.d0*(adotoa**2.d0-Hdot)/k2
            
			etakdot = k/(2.d0*TGR_f_muSigma)*(k*betak*(2.d0*TGR_Sigma-TGR_mu)*TGR_f_1*dgq + 2.d0*(TGR_mu-TGR_Sigma)*pidot_sum &
					+ TGR_rhoDelta*betak*((TGR_mudot-2.d0*TGR_Sigmadot)+2.d0*adotoa*(TGR_Sigma-TGR_mu)) + 2.d0*dgpi*(TGR_mudot & 
					-TGR_Sigmadot + adotoa*betak*(2.d0*TGR_Sigma-TGR_mu) - adotoa*TGR_mu) + k*sigma*(-2.d0*(adotoa**2.d0-Hdot) + betak* &
					(2.d0*TGR_Sigma-TGR_mu)*((grhoc_t+grhob_t)+(grhor_t+grhog_t)*4.d0/3.d0+(grhonu_t+gpres_nu))) &
                    + 2.d0*adotoa*(TGR_Sigma-TGR_mu)*dgpi_3wplus1)
                    
			ayprime(2) = etakdot
			!Computing z
			z = sigma - 3.d0*etakdot/k2
	
			TGR_Phidot = (etakdot - adotoa*sigmadot - Hdot*sigma)/k
			
		else
    
		    TGR_f_q = k2 + 1.5d0*betak*TGR_Q_T*((grhoc_t+grhob_t)+(grhor_t+grhog_t)*4.d0/3.d0+(grhonu_t+gpres_nu)) 
			TGR_f_1 = 1+3.d0*(adotoa**2.d0-Hdot)/k2
			
		    etakdot = -k/(2.d0*TGR_f_q)*(2.d0*(adotoa**2.d0-Hdot)*k*sigma+(2.d0*adotoa*(TGR_D_T-TGR_Q_T)+TGR_Q_T_dot)*betak*TGR_rhoDelta &
			-k*sigma*betak*TGR_Q_T*(grhoc_t+grhob_t+(grhor_t+grhog_t)*4.d0/3.d0+(grhonu_t+gpres_nu)) &
			- k*betak*TGR_Q_T*TGR_f_1*dgq - 2.d0*adotoa*TGR_Q_T*(betak-1.d0)*dgpi)

		    ayprime(2) = etakdot

		    !JD get z from etakdot and sigma
		    z= sigma-3.d0*etakdot/k2
			
			TGR_Phidot = (etakdot - adotoa*sigmadot - Hdot*sigma)/k

		end if
        
    else !CGQ to work with default GR or with MG models
    	!  Get sigma (shear) and z from the constraints
    	! have to get z from eta for numerical stability
    	z=(0.5_dl*dgrho/k + etak)/adotoa
    	if (CP%flat) then
        	!eta*k equation
        	sigma=(z+1.5_dl*dgq/k2)
        	ayprime(2)=0.5_dl*dgq
    	else
    	    sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)
    	    ayprime(2)=0.5_dl*dgq + CP%curv*z
    	end if
     
	end if	!here it ends for MG patch
	!<ISiTGR MOD END -----------------------------------------------------------------------

	if (TGR%GR==0) then
	    if (w_lam /= -1 .and. w_Perturb) then
	        ayprime(EV%w_ix)= -3*adotoa*(cs2_lam-w_lam)*(clxde+3*adotoa*qde/k) &
	            - k*qde -(1+w_lam)*k*z

	        ayprime(EV%w_ix+1) = (-adotoa*(1-3*cs2_lam)*qde + k*cs2_lam*clxde)/(1+w_lam)
	    end if
	end if
	
    !  CDM equation of motion
    clxcdot=-k*z
    ayprime(3)=clxcdot

    !  Baryon equation of motion.
    clxbdot=-k*(z+vb)
    ayprime(4)=clxbdot
    !  Photon equation of motion
    clxgdot=-k*(4._dl/3._dl*z+qg)

    ! old comment:Small k: potential problem with stability, using full equations earlier is NOT more accurate in general
    ! Easy to see instability in k \sim 1e-3 by tracking evolution of vb

    !  Use explicit equation for vb if appropriate

    if (EV%TightCoupling) then
        !  ddota/a
		!>ISiTGR MOD START
		if (TGR%GR==0) then
        	gpres=gpres_nu + (grhog_t+grhor_t)/3 +grhov_t*w_lam
		else
			gpres=gpres_nu + (grhog_t+grhor_t)/3 + gpresv_t
		end if
		!<ISiTGR MOD END
			adotdota=(adotoa*adotoa-gpres)/2
			pig = 32._dl/45/opacity*k*(sigma+vb)
        !  First-order approximation to baryon-photon splip
        slip = - (2*adotoa/(1+pb43) + dopacity/opacity)* (vb-3._dl/4*qg) &
            +(-adotdota*vb-k/2*adotoa*clxg +k*(cs2*clxbdot-clxgdot/4))/(opacity*(1+pb43))

        if (second_order_tightcoupling) then
            ! by Francis-Yan Cyr-Racine simplified (inconsistently) by AL assuming flat
            !AL: First order slip seems to be fine here to 2e-4
			if (TGR%GR==0) then
            !  8*pi*G*a*a*SUM[rho_i*sigma_i]
            	dgs = grhog_t*pig+grhor_t*pir

            ! Define shear derivative to first order
            	sigmadot = -2*adotoa*sigma-dgs/k+etak
			else
				sigmadot = k*TGR_Psi - adotoa*sigma
			end if
            !Once know slip, recompute qgdot, pig, pigdot
            qgdot = k*(clxg/4._dl-pig/2._dl) +opacity*slip

            pig = 32._dl/45/opacity*k*(sigma+3._dl*qg/4._dl)*(1+(dopacity*11._dl/6._dl/opacity**2)) &
                + (32._dl/45._dl/opacity**2)*k*(sigmadot+3._dl*qgdot/4._dl)*(-11._dl/6._dl)

            pigdot = -(32._dl/45._dl)*(dopacity/opacity**2)*k*(sigma+3._dl*qg/4._dl)*(1 + &
                dopacity*11._dl/6._dl/opacity**2 ) &
                + (32._dl/45._dl/opacity)*k*(sigmadot+3._dl*qgdot/4._dl)*(1+(11._dl/6._dl) &
                *(dopacity/opacity**2))

            EV%pigdot = pigdot

        end if

        !  Use tight-coupling approximation for vb
        !  zeroth order approximation to vbdot + the pig term
        vbdot=(-adotoa*vb+cs2*k*clxb  &
            +k/4*pb43*(clxg-2*EV%Kf(1)*pig))/(1+pb43)

        vbdot=vbdot+pb43/(1+pb43)*slip
        EV%pig = pig
    else
        vbdot=-adotoa*vb+cs2*k*clxb-photbar*opacity*(4._dl/3*vb-qg)
    end if

    ayprime(5)=vbdot

    if (.not. EV%no_phot_multpoles) then
        !  Photon equations of motion
        ayprime(EV%g_ix)=clxgdot
        qgdot=4._dl/3*(-vbdot-adotoa*vb+cs2*k*clxb)/pb43 &
            +EV%denlk(1)*clxg-EV%denlk2(1)*pig
        ayprime(EV%g_ix+1)=qgdot

        !  Use explicit equations for photon moments if appropriate
        if (.not. EV%tightcoupling) then
            E2=ay(EV%polind+2)
            polter = pig/10+9._dl/15*E2 !2/15*(3/4 pig + 9/2 E2)
            ix= EV%g_ix+2
            if (EV%lmaxg>2) then
                pigdot=EV%denlk(2)*qg-EV%denlk2(2)*ay(ix+1)-opacity*(pig - polter) &
                    +8._dl/15._dl*k*sigma
                ayprime(ix)=pigdot
                do  l=3,EV%lmaxg-1
                    ix=ix+1
                    ayprime(ix)=(EV%denlk(l)*ay(ix-1)-EV%denlk2(l)*ay(ix+1))-opacity*ay(ix)
                end do
                ix=ix+1
                !  Truncate the photon moment expansion
                ayprime(ix)=k*ay(ix-1)-(EV%lmaxg+1)*cothxor*ay(ix) -opacity*ay(ix)
            else !closed case
                pigdot=EV%denlk(2)*qg-opacity*(pig - polter) +8._dl/15._dl*k*sigma
                ayprime(ix)=pigdot
            endif
            !  Polarization
            !l=2
            ix=EV%polind+2
            if (EV%lmaxgpol>2) then
                ayprime(ix) = -opacity*(ay(ix) - polter) - k/3._dl*ay(ix+1)
                do l=3,EV%lmaxgpol-1
                    ix=ix+1
                    ayprime(ix)=-opacity*ay(ix) + (EV%denlk(l)*ay(ix-1)-EV%polfack(l)*ay(ix+1))
                end do
                ix=ix+1
                !truncate
                ayprime(ix)=-opacity*ay(ix) + &
                    k*EV%poltruncfac*ay(ix-1)-(EV%lmaxgpol+3)*cothxor*ay(ix)
            else !closed case
                ayprime(ix) = -opacity*(ay(ix) - polter)
            endif
        end if
    end if

    if (.not. EV%no_nu_multpoles) then
        !  Massless neutrino equations of motion.
        clxrdot=-k*(4._dl/3._dl*z+qr)
        ayprime(EV%r_ix)=clxrdot
        qrdot=EV%denlk(1)*clxr-EV%denlk2(1)*pir
        ayprime(EV%r_ix+1)=qrdot
        if (EV%high_ktau_neutrino_approx) then
            !ufa approximation for k*tau>>1, more accurate when there are reflections from lmax
            !Method from arXiv:1104.2933
            !                if (.not. EV%TightCoupling) then
            !                 gpres=gpres+ (grhog_t+grhor_t)/3 +grhov_t*w_lam
            !                 adotdota=(adotoa*adotoa-gpres)/2
            !                end if
            !                ddz=(2*adotoa**2 - adotdota)*z  &
            !                  + adotoa/(2*k)*( 6*(grhog_t*clxg+grhor_t*clxr) + 2*(grhoc_t*clxc+grhob_t*clxb) ) &
            !                   - 1._dl/(2*k)*( 2*(grhog_t*clxgdot+grhor_t*clxrdot) + grhoc_t*clxcdot + grhob_t*clxbdot )
            !                dz= -adotoa*z - 0.5_dl*dgrho/k
            !                pirdot= -3*pir*cothxor + k*(qr+4._dl/3*z)
            pirdot= -3*pir*cothxor - clxrdot
            ayprime(EV%r_ix+2)=pirdot

            !                pirdot=k*(0.4_dl*qr-0.6_dl*ay(EV%lmaxg+10)+8._dl/15._dl*sigma)
            !                ayprime(EV%lmaxg+9)=pirdot
            !                ayprime(3+EV%lmaxg+7)=k*ay(3+EV%lmaxg+6)- &
            !                                      (3+1)*cothxor*ay(3+EV%lmaxg+7)
            !               ayprime(3+EV%lmaxg+7+1:EV%lmaxnr+EV%lmaxg+7)=0
        else
            ix=EV%r_ix+2
            if (EV%lmaxnr>2) then
                pirdot=EV%denlk(2)*qr- EV%denlk2(2)*ay(ix+1)+8._dl/15._dl*k*sigma
                ayprime(ix)=pirdot
                do l=3,EV%lmaxnr-1
                    ix=ix+1
                    ayprime(ix)=(EV%denlk(l)*ay(ix-1) - EV%denlk2(l)*ay(ix+1))
                end do
                !  Truncate the neutrino expansion
                ix=ix+1
                ayprime(ix)=k*ay(ix-1)- (EV%lmaxnr+1)*cothxor*ay(ix)
            else
                pirdot=EV%denlk(2)*qr +8._dl/15._dl*k*sigma
                ayprime(ix)=pirdot
            end if
        end if
    end if ! no_nu_multpoles

    !  Massive neutrino equations of motion.
    if (CP%Num_Nu_massive >0) then
        !DIR$ LOOP COUNT MIN(1), AVG(1)
        do nu_i = 1, CP%Nu_mass_eigenstates
            if (EV%MassiveNuApprox(nu_i)) then
                !Now EV%iq0 = clx, EV%iq0+1 = clxp, EV%iq0+2 = G_1, EV%iq0+3=G_2=pinu
                !see astro-ph/0203507
                G11_t=EV%G11(nu_i)/a/a2
                G30_t=EV%G30(nu_i)/a/a2
                off_ix = EV%nu_ix(nu_i)
                w=wnu_arr(nu_i)
                ayprime(off_ix)=-k*z*(w+1) + 3*adotoa*(w*ay(off_ix) - ay(off_ix+1))-k*ay(off_ix+2)
                ayprime(off_ix+1)=(3*w-2)*adotoa*ay(off_ix+1) - 5._dl/3*k*z*w - k/3*G11_t
                ayprime(off_ix+2)=(3*w-1)*adotoa*ay(off_ix+2) - k*(2._dl/3*EV%Kf(1)*ay(off_ix+3)-ay(off_ix+1))
                ayprime(off_ix+3)=(3*w-2)*adotoa*ay(off_ix+3) + 2*w*k*sigma - k/5*(3*EV%Kf(2)*G30_t-2*G11_t)
            else
                ind=EV%nu_ix(nu_i)
                !DIR$ LOOP COUNT MIN(3), AVG(3)
                do i=1,EV%nq(nu_i)
                    q=nu_q(i)
                    aq=a*nu_masses(nu_i)/q
                    v=1._dl/sqrt(1._dl+aq*aq)

                    ayprime(ind)=-k*(4._dl/3._dl*z + v*ay(ind+1))
                    ind=ind+1
                    ayprime(ind)=v*(EV%denlk(1)*ay(ind-1)-EV%denlk2(1)*ay(ind+1))
                    ind=ind+1
                    if (EV%lmaxnu_tau(nu_i)==2) then
                        ayprime(ind)=-ayprime(ind-2) -3*cothxor*ay(ind)
                    else
                        ayprime(ind)=v*(EV%denlk(2)*ay(ind-1)-EV%denlk2(2)*ay(ind+1)) &
                            +k*8._dl/15._dl*sigma
                        do l=3,EV%lmaxnu_tau(nu_i)-1
                            ind=ind+1
                            ayprime(ind)=v*(EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
                        end do
                        !  Truncate moment expansion.
                        ind = ind+1
                        ayprime(ind)=k*v*ay(ind-1)-(EV%lmaxnu_tau(nu_i)+1)*cothxor*ay(ind)
                    end if
                    ind = ind+1
                end do
            end if
        end do

        if (EV%has_nu_relativistic) then
            ind=EV%nu_pert_ix
            ayprime(ind)=+k*a2*qr -k*ay(ind+1)
            ind2= EV%r_ix
            do l=1,EV%lmaxnu_pert-1
                ind=ind+1
                ind2=ind2+1
                ayprime(ind)= -a2*(EV%denlk(l)*ay(ind2-1)-EV%denlk2(l)*ay(ind2+1)) &
                    +   (EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
            end do
            ind=ind+1
            ind2=ind2+1
            ayprime(ind)= k*(ay(ind-1) -a2*ay(ind2-1)) -(EV%lmaxnu_pert+1)*cothxor*ay(ind)
        end if
    end if

    if (associated(EV%OutputTransfer) .or. associated(EV%OutputSources)) then
        if (EV%TightCoupling .or. EV%no_phot_multpoles) then
            E=0
            Edot=0
        else
            E = ay(EV%polind+2:EV%polind+3)
            Edot = ayprime(EV%polind+2:EV%polind+3)
        end if
        if (EV%no_nu_multpoles) then
            pirdot=0
            qrdot = -4*dz/3
        end if
        if (EV%no_phot_multpoles) then
            pigdot=0
            octg=0
            octgdot=0
           	qgdot = -4*dz/3
        else
            if (EV%TightCoupling) then
                if (second_order_tightcoupling) then
                    octg = (3._dl/7._dl)*pig*(EV%k_buf/opacity)
                    E(2) = pig/4 + pigdot*(1._dl/opacity)*(-5._dl/8._dl)
                    E(3) = (3._dl/7._dl)*(EV%k_buf/opacity)*E(2)
                    Edot(2)= (pigdot/4._dl)*(1+(5._dl/2._dl)*(dopacity/opacity**2))
                else
                    pigdot = -dopacity/opacity*pig + 32._dl/45*k/opacity*(-2*adotoa*sigma  &
                        +etak/EV%Kf(1)-  dgpi/k +vbdot )
                    Edot(2) = pigdot/4
                    E(2) = pig/4
                    octg=0
                end if
                octgdot=0
            else
                octg=ay(EV%g_ix+3)
                octgdot=ayprime(EV%g_ix+3)
            end if
        end if

        dgpi  = grhor_t*pir + grhog_t*pig
        dgpi_diff = 0  !sum (3*p_nu -rho_nu)*pi_nu
        pidot_sum = grhog_t*pigdot + grhor_t*pirdot
        clxnu =0
        if (CP%Num_Nu_Massive /= 0) then
            call MassiveNuVarsOut(EV,ay,ayprime,a, dgpi=dgpi, clxnu_all=clxnu, &
                dgpi_diff=dgpi_diff, pidot_sum=pidot_sum)
        end if
        diff_rhopi = pidot_sum - (4*dgpi+ dgpi_diff)*adotoa
		
		!>ISiTGR MOD START
		if (TGR%GR==0) then
			gpres=gpres_nu+ (grhog_t+grhor_t)/3 +grhov_t*w_lam
		else
			gpres=gpres_nu+ (grhog_t+grhor_t)/3 + gpresv_t
		end	if
		
		!<ISiTGR MOD END
        !CGQ ---------------------------
		if (TGR%GR==0) then
        phi = -((dgrho +3*dgq*adotoa/k)/EV%Kf(1) + dgpi)/(2*k2)
		else
		phi = (TGR_Psi+TGR_Phi)/2._dl
        end if
        !CGQ ---------------------------
		
		if (associated(EV%OutputTransfer)) then
            EV%OutputTransfer(Transfer_kh) = k/(CP%h0/100._dl)
            EV%OutputTransfer(Transfer_cdm) = clxc
            EV%OutputTransfer(Transfer_b) = clxb
            EV%OutputTransfer(Transfer_g) = clxg
            EV%OutputTransfer(Transfer_r) = clxr
            EV%OutputTransfer(Transfer_nu) = clxnu
            EV%OutputTransfer(Transfer_tot) =  dgrho_matter/grho_matter !includes neutrinos
            EV%OutputTransfer(Transfer_nonu) = (grhob_t*clxb+grhoc_t*clxc)/(grhob_t + grhoc_t)
            EV%OutputTransfer(Transfer_tot_de) =  dgrho/grho_matter
            !Transfer_Weyl is k^2Phi, where Phi is the Weyl potential
			!>ISiTGR MOD START: Weyl function for GR or MG
			if (TGR%GR==0) then !CGQ
          		EV%OutputTransfer(Transfer_Weyl) = k2*phi
			else
	        	EV%OutputTransfer(Transfer_Weyl) = k2*(TGR_Phi+TGR_Psi)/2.d0 !CGQ
			end if
			!<ISiTGR MOD END
            EV%OutputTransfer(Transfer_Newt_vel_cdm)=  -k*sigma/adotoa
            EV%OutputTransfer(Transfer_Newt_vel_baryon) = -k*(vb + sigma)/adotoa
            EV%OutputTransfer(Transfer_vel_baryon_cdm) = vb
			!>ISiTGR MOD START: Different outputs for ISW and WL likelihood for MG
	        !(\phi' + \psi' in Newtonian gauge)
			if(TGR%GR==0) then !CGQ
			    phidot = (1.0d0/2.0d0)*(adotoa*(-dgpi - 2*k2*phi) + dgq*k - &
                diff_rhopi+ k*sigma*(gpres + grho))/k2
				EV%OutputTransfer(Transfer_ISW) = 2*phidot
			else
    		    !(\phi' + \psi' in Newtonian gauge)
    		    EV%OutputTransfer(Transfer_ISW) =TGR_Phidot + TGR_Psidot !CGQ
	        end if
			!Logarithmic growth rate d \ln \delta/d \ln a  (k2 here cause it is divided out later)
        	EV%OutputTransfer(Transfer_f) = k2*(clxcdot*grhoc_t + clxbdot*grhob_t)/(clxc*grhoc_t + clxb*grhob_t)/adotoa
    	    !Transfer function for velocity power spectrum (in Conformal Newtonian Guage)
	        !In units consistent with MPT Breeze
        	EV%OutputTransfer(Transfer_vtot) = c/1000*(sigma + grhob_t*vb/(grhob_t+grhoc_t))
		end if
			!<ISiTGR MOD END
           		
        if (associated(EV%OutputSources)) then

            call IonizationFunctionsAtTime(tau, opacity, dopacity, ddopacity, &
                visibility, dvisibility, ddvisibility, exptau, lenswindow)

            tau0 = CP%tau0
            
        !>ISiTGR MOD START: CGQ to work with GR and MG ------------------------------------------
        if (TGR%GR==0) then
        
            phidot = (1.0d0/2.0d0)*(adotoa*(-dgpi - 2*k2*phi) + dgq*k - &
                    diff_rhopi+ k*sigma*(gpres + grho))/k2
                    
           !time derivative of shear
            sigmadot = -adotoa*sigma - 1.0d0/2.0d0*dgpi/k + k*phi
            
           !quadrupole source derivatives; polter = pi_g/10 + 3/5 E_2
            polter = pig/10+9._dl/15*E(2)
            polterdot = (1.0d0/10.0d0)*pigdot + (3.0d0/5.0d0)*Edot(2)
            polterddot = -2.0d0/25.0d0*adotoa*dgq/(k*EV%Kf(1)) - 4.0d0/75.0d0*adotoa* &
                k*sigma - 4.0d0/75.0d0*dgpi - 2.0d0/75.0d0*dgrho/EV%Kf(1) - 3.0d0/ &
                50.0d0*k*octgdot*EV%Kf(2) + (1.0d0/25.0d0)*k*qgdot - 1.0d0/5.0d0 &
                *k*EV%Kf(2)*Edot(3) + (-1.0d0/10.0d0*pig + (7.0d0/10.0d0)* &
                polter - 3.0d0/5.0d0*E(2))*dopacity + (-1.0d0/10.0d0*pigdot &
                + (7.0d0/10.0d0)*polterdot - 3.0d0/5.0d0*Edot(2))*opacity
            !Temperature source terms, after integrating by parts in conformal time

            !2phi' term (\phi' + \psi' in Newtonian gauge), phi is the Weyl potential
            ISW = 2*phidot*exptau
            monopole_source =  (-etak/(k*EV%Kf(1)) + 2*phi + clxg/4)*visibility
            doppler = ((sigma + vb)*dvisibility + (sigmadot + vbdot)*visibility)/k
            quadrupole_source = (5.0d0/8.0d0)*(3*polter*ddvisibility + 6*polterdot*dvisibility &
                + (k**2*polter + 3*polterddot)*visibility)/k**2

            EV%OutputSources(1) = ISW + doppler + monopole_source + quadrupole_source
		
		else

	        if ((TGR%ISiTGR_mueta .eqv. .true.) .or. (TGR%ISiTGR_muSigma .eqv. .true.) &
				.or. (TGR%ISiTGR_BIN_mueta .eqv. .true.) .or. (TGR%ISiTGR_BIN_muSigma .eqv. .true.)) then !CGQ for ISW effect for mueta, muSigma parameterizations
    	        !Computing Psidot for ISW
    	        TGR_Psidot = -TGR_mudot/(2.d0*k2)*(betak*TGR_rhoDelta+2.d0*dgpi) &
						+ TGR_mu/(2.d0*k2)*(adotoa*betak*TGR_rhoDelta - 2.d0*pidot_sum + &
        	    		k*betak*TGR_f_1*dgq + 2.d0*adotoa*dgpi_3wplus1 + 2.d0*betak*adotoa*dgpi &
						+ (TGR_f_1*sigma*k-3.d0*(TGR_Phidot+adotoa*TGR_Psi))* &
						betak*(grhoc_t+grhob_t+(grhor_t+grhog_t)*4.d0/3.d0+grhonu_t+gpres_nu))
					
	        else
			
        	    TGR_Psidot = -TGR_Phidot - 1.d0/k2 *(betak*TGR_D_T_dot*TGR_rhoDelta &
						+ TGR_Q_T_dot*dgpi - betak*TGR_D_T*adotoa*TGR_rhoDelta &
						-TGR_Q_T*adotoa*dgpi_3wplus1+TGR_Q_T*pidot_sum + betak*(3.d0*TGR_D_T*(TGR_Phidot + adotoa*TGR_Psi) &
						-k*sigma*TGR_f_1*TGR_D_T)*(grhoc_t+grhob_t+(grhor_t+grhog_t)*4.d0/3.d0+grhonu_t+gpres_nu) &
						- 2.d0*betak*TGR_D_T*adotoa*dgpi - k*betak*TGR_D_T*TGR_f_1*dgq)

			end if
			
            !phidot plus psidot        
			ISW = (TGR_Psidot+TGR_Phidot)*exptau

 		    EV%OutputSources(1) = ISW + visibility*pig/16.d0+(3.D0/8.D0*E(2)+clxg/4.d0)*visibility+(11.D0/10.D0*dvisibility &
							    *sigma+(-3.D0/8.D0*EV%Kf(2)*E(3)-9.D0/80.D0*EV%Kf(2)*octg+3.D0/40.D0*qg+vb) &
 			  				    *dvisibility+(3.D0/40.D0*qgdot+21.D0/10.D0*sigmadot+vbdot-9.D0/80.D0*EV%Kf(2) &
							    *octgdot-3.D0/8.D0*EV%Kf(2)*Edot(3))*visibility)/k+((3.D0/16.D0*ddvisibility &
 							    -9.D0/160.D0*visibility*dopacity-9.D0/160.D0*dvisibility*opacity)*pig+(9.D0/8.D0 &
    							*Edot(2)+3.D0/16.D0*pigdot-27.D0/80.D0*opacity*E(2))*dvisibility &
    							+((-9.D0/160.D0*pigdot-27.D0/80.D0*Edot(2))*opacity-27.D0/80.D0 &
    							*dopacity*E(2))*visibility+9.D0/8.D0*ddvisibility*E(2))/k**2

        end if
        !<ISiTGR MOD END -------------------------------------------------------------------
        
			!>ISiTGR MOD START
            if (tau < tau0) then
                !E polarization source
                EV%OutputSources(2)=visibility*polter*(15._dl/8._dl)/(f_K(tau0-tau)**2*k2)
                !factor of four because no 1/16 later
            else
                EV%OutputSources(2)=0
            end if
			
            if (size(EV%OutputSources) > 2) then
                !Get lensing sources
                !Can modify this here if you want to get power spectra for other tracer
                if (tau>tau_maxvis .and. tau0-tau > 0.1_dl) then
					if (TGR%GR==0) then !CGQ
       					phi = -((dgrho +3*dgq*adotoa/k)/EV%Kf(1) + dgpi)/(2*k2)
                    	EV%OutputSources(3) = -2*phi*f_K(tau-tau_maxvis)/(f_K(tau0-tau_maxvis)*f_K(tau0-tau))
					else 
						EV%OutputSources(3) = -(TGR_Phi+TGR_Psi)*f_K(tau-tau_maxvis)/(f_K(tau0-tau_maxvis)*f_K(tau0-tau)) !CGQ
					end if
                else
                    EV%OutputSources(3) = 0
                end if

            !<ISiTGR MOD END
			end if
            if (associated(EV%CustomSources)) then
                call custom_sources_func(EV%CustomSources, tau, a, adotoa, grho, gpres,w_lam, cs2_lam, &
                    grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t, &
                    k, etak, ayprime(2), phi, phidot, sigma, sigmadot, &
                    dgrho, clxg,clxb,clxc,clxr, clxnu, clxde, cs2*clxb, &
                    dgq, qg, qr, qde, vb, qgdot, qrdot, vbdot, &
                    dgpi, pig, pir, pigdot, pirdot, diff_rhopi, &
                    polter, polterdot, polterddot, octg, octgdot, E, Edot, &
                    opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau, &
                    tau0, tau_maxvis, EV%Kf,f_K)
            end if
        end if
    end if

    end subroutine derivs



    subroutine derivsv(EV,n,tau,yv,yvprime)
    !  Evaluate the time derivatives of the vector perturbations, flat case
    use ThermoData
    use MassiveNu
    implicit none
    type(EvolutionVars) EV
    integer n,l
    real(dl), target ::  yv(n),yvprime(n)
    real(dl) ep,tau,grho,rhopi,cs2,opacity,gpres
    logical finished_tightcoupling
    real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
    real(dl)  grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,polter
    real(dl) sigma, qg,pig, qr, vb, rhoq, vbdot, photbar, pb43
    real(dl) k,k2,a,a2, adotdota
    real(dl) pir,adotoa

    k2=EV%k2_buf
    k=EV%k_buf

    !E and B start at l=2. Set up pointers accordingly to fill in y arrays
    E => yv(EV%lmaxv+3:)
    Eprime=> yvprime(EV%lmaxv+3:)
    B => E(EV%lmaxpolv:)
    Bprime => Eprime(EV%lmaxpolv:)
    neutprime => Bprime(EV%lmaxpolv+1:)
    neut => B(EV%lmaxpolv+1:)

    a=yv(1)

    sigma=yv(2)

    a2=a*a

    !  Get sound speed and opacity, and see if should use tight-coupling

    call thermo(tau,cs2,opacity)
    if (k > 0.06_dl*epsw) then
        ep=ep0
    else
        ep=0.2_dl*ep0
    end if

    finished_tightcoupling = &
        ((k/opacity > ep).or.(1._dl/(opacity*tau) > ep .and. k/opacity > 1d-4))


    ! Compute expansion rate from: grho=8*pi*rho*a**2
    ! Also calculate gpres: 8*pi*p*a**2
    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2
    grhov_t=grhov*a**(-1-3*w_lam)

    grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
    gpres=(grhog_t+grhor_t)/3._dl+grhov_t*w_lam

    adotoa=sqrt(grho/3._dl)
    adotdota=(adotoa*adotoa-gpres)/2

    photbar=grhog_t/grhob_t
    pb43=4._dl/3*photbar

    yvprime(1)=adotoa*a

    vb = yv(3)
    qg = yv(4)
    qr = neut(1)

    !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
    rhoq=grhob_t*vb+grhog_t*qg+grhor_t*qr
    !  sigma = 2*rhoq/k**2
    !for non-large k this expression for sigma is unstable at early times
    !so propagate sigma equation separately (near total cancellation in rhoq)
    ! print *,yv(2),2*rhoq/k**2

    if (finished_tightcoupling) then
        !  Use explicit equations:

        pig = yv(5)

        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)

        vbdot = -adotoa*vb-photbar*opacity*(4._dl/3*vb-qg) - 0.5_dl*k*photbar*Magnetic

        !  Equation for the photon heat flux stress

        yvprime(4)=-0.5_dl*k*pig + opacity*(4._dl/3*vb-qg)

        !  Equation for the photon anisotropic stress
        yvprime(5)=k*(2._dl/5*qg -8/15._dl*yv(6))+8._dl/15._dl*k*sigma  &
            -opacity*(pig - polter)
        ! And for the moments
        do  l=3,EV%lmaxv-1
            yvprime(l+3)=k*denl(l)*l*(yv(l+2)-   &
                vecfac(l)*yv(l+4))-opacity*yv(l+3)
        end do
        !  Truncate the hierarchy
        yvprime(EV%lmaxv+3)=k*EV%lmaxv/(EV%lmaxv-1._dl)*yv(EV%lmaxv+2)- &
            (EV%lmaxv+2._dl)*yv(EV%lmaxv+3)/tau-opacity*yv(EV%lmaxv+3)

        !E equations

        Eprime(2) = - opacity*(E(2) - polter) + k*(1/3._dl*B(2) - 8._dl/27._dl*E(3))
        do l=3,EV%lmaxpolv-1
            Eprime(l) =-opacity*E(l) + k*(denl(l)*(l*E(l-1) - &
                vecfacpol(l)*E(l+1)) + 2._dl/(l*(l+1))*B(l))
        end do
        !truncate
        Eprime(EV%lmaxpolv)=0._dl

        !B-bar equations

        do l=2,EV%lmaxpolv-1
            Bprime(l) =-opacity*B(l) + k*(denl(l)*(l*B(l-1) - &
                vecfacpol(l)*B(l+1)) - 2._dl/(l*(l+1))*E(l))
        end do
        !truncate
        Bprime(EV%lmaxpolv)=0._dl
    else
        !Tight coupling expansion results

        pig = 32._dl/45._dl*k/opacity*(vb + sigma)

        EV%pig = pig

        vbdot=(-adotoa*vb  -3._dl/8*pb43*k*Magnetic  -3._dl/8*k*pb43*pig &
            - pb43/(1+pb43)/opacity*(0.75_dl*k*adotoa*pb43**2/(pb43+1)*Magnetic + vb*&
            ( 2*pb43*adotoa**2/(1+pb43) + adotdota)) &
            )/(1+pb43)

        !  Equation for the photon heat flux
        ! Get drag from vbdot expression
        yvprime(4)=-0.5_dl*k*pig - &
            (vbdot+adotoa*vb)/photbar - 0.5_dl*k*Magnetic

        !  Set the derivatives to zero
        yvprime(5:n)=0._dl
        yv(5)=pig
        E(2)=  pig/4
    endif

    yvprime(3) = vbdot

    !  Neutrino equations:

    !  Massless neutrino anisotropic stress
    pir=neut(2)
    neutprime(1)= -0.5_dl*k*pir
    neutprime(2)=2._dl/5*k*qr -8._dl/15._dl*k*neut(3)+ 8._dl/15._dl*k*sigma
    !  And for the moments
    do  l=3,EV%lmaxnrv-1
        neutprime(l)=k*denl(l)*l*(neut(l-1)- vecfac(l)*neut(l+1))
    end do

    !  Truncate the hierarchy
    neutprime(EV%lmaxnrv)=k*EV%lmaxnrv/(EV%lmaxnrv-1._dl)*neut(EV%lmaxnrv-1)-  &
        (EV%lmaxnrv+2._dl)*neut(EV%lmaxnrv)/tau


    !  Get the propagation equation for the shear

    rhopi=grhog_t*pig+grhor_t*pir+ grhog_t*Magnetic

    yvprime(2)=-2*adotoa*sigma -rhopi/k

    end subroutine derivsv



    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine derivst(EV,n,tau,ayt,aytprime)
    !  Evaluate the time derivatives of the tensor perturbations.
    use ThermoData
    use MassiveNu
    implicit none
    type(EvolutionVars) EV
    integer n,l,i,ind, nu_i
    real(dl), target ::  ayt(n),aytprime(n)
    real(dl) tau,grho,rhopi,cs2,opacity,pirdt
    real(dl), dimension(:),pointer :: neut,neutprime,E,B,Eprime,Bprime
    real(dl) q,aq,v
    real(dl)  grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,polter
    real(dl) Hchi,pinu, pig
    real(dl) k,k2,a,a2
    real(dl) pir, adotoa, rhonu, shear

    real(dl) cothxor

    k2=EV%k2_buf
    k= EV%k_buf

    a=ayt(1)

    Hchi=ayt(2)

    shear=ayt(3)

    a2=a*a

    ! Compute expansion rate from: grho=8*pi*rho*a**2
    ! Also calculate gpres: 8*pi*p*a**2
    grhob_t=grhob/a
    grhoc_t=grhoc/a
    grhor_t=grhornomass/a2
    grhog_t=grhog/a2
    if (w_lam==-1._dl) then
        grhov_t=grhov*a2
    else
        grhov_t=grhov*a**(-1-3*w_lam)
    end if

    grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t

    !Do massive neutrinos
    if (CP%Num_Nu_Massive >0) then
        do nu_i=1,CP%Nu_mass_eigenstates
            call Nu_rho(a*nu_masses(nu_i),rhonu)
            grho=grho+grhormass(nu_i)*rhonu/a2
        end do
    end if

    if (CP%flat) then
        cothxor=1._dl/tau
        adotoa=sqrt(grho/3._dl)
    else
        cothxor=1._dl/tanfunc(tau/CP%r)/CP%r
        adotoa=sqrt((grho+grhok)/3._dl)
    end if

    aytprime(1)=adotoa*a

    call thermo(tau,cs2,opacity)

    if (.not. EV%TensTightCoupling) then
        !  Don't use tight coupling approx - use explicit equations:
        !  Equation for the photon anisotropic stress


        !E and B start at l=2. Set up pointers accordingly to fill in ayt arrays
        E => ayt(EV%E_ix+1:)
        B => ayt(EV%B_ix+1:)
        Eprime=> aytprime(EV%E_ix+1:)
        Bprime => aytprime(EV%B_ix+1:)

        ind = EV%g_ix+2

        !  Photon anisotropic stress
        pig=ayt(ind)
        polter = 0.1_dl*pig + 9._dl/15._dl*E(2)

        if (EV%lmaxt > 2) then
            aytprime(ind)=-EV%denlkt(2,2)*ayt(ind+1)+k*8._dl/15._dl*shear  &
                -opacity*(pig - polter)

            do l=3, EV%lmaxt -1
                ind = ind+1
                aytprime(ind)=EV%denlkt(1,L)*ayt(ind-1)-EV%denlkt(2,L)*ayt(ind+1)-opacity*ayt(ind)
            end do

            !Truncate the hierarchy
            ind=ind+1
            aytprime(ind)=k*EV%lmaxt/(EV%lmaxt-2._dl)*ayt(ind-1)- &
                (EV%lmaxt+3._dl)*cothxor*ayt(ind)-opacity*ayt(ind)

            !E and B-bar equations

            Eprime(2) = - opacity*(E(2) - polter) + EV%denlkt(4,2)*B(2) - &
                EV%denlkt(3,2)*E(3)

            do l=3, EV%lmaxpolt-1
                Eprime(l) =(EV%denlkt(1,L)*E(l-1)-EV%denlkt(3,L)*E(l+1) + EV%denlkt(4,L)*B(l)) &
                    -opacity*E(l)
            end do
            l= EV%lmaxpolt
            !truncate: difficult, but setting l+1 to zero seems to work OK
            Eprime(l) = (EV%denlkt(1,L)*E(l-1) + EV%denlkt(4,L)*B(l)) -opacity*E(l)

            Bprime(2) =-EV%denlkt(3,2)*B(3) - EV%denlkt(4,2)*E(2)  -opacity*B(2)
            do l=3, EV%lmaxpolt-1
                Bprime(l) =(EV%denlkt(1,L)*B(l-1) -EV%denlkt(3,L)*B(l+1) - EV%denlkt(4,L)*E(l)) &
                    -opacity*B(l)
            end do
            l=EV%lmaxpolt
            !truncate
            Bprime(l) =(EV%denlkt(1,L)*B(l-1) - EV%denlkt(4,L)*E(l))  -opacity*B(l)

        else !lmax=2

            aytprime(ind)=k*8._dl/15._dl*shear-opacity*(pig - polter)
            Eprime(2) = - opacity*(E(2) - polter) + EV%denlkt(4,2)*B(2)
            Bprime(2) = - EV%denlkt(4,2)*E(2)  -opacity*B(2)
        end if

    else  !Tight coupling
        pig = 32._dl/45._dl*k/opacity*shear
    endif

    rhopi=grhog_t*pig


    !  Neutrino equations:
    !  Anisotropic stress
    if (DoTensorNeutrinos) then
        neutprime => aytprime(EV%r_ix+1:)
        neut => ayt(EV%r_ix+1:)

        !  Massless neutrino anisotropic stress
        pir=neut(2)

        rhopi=rhopi+grhor_t*pir

        if (EV%lmaxnrt>2) then
            pirdt=-EV%denlkt(2,2)*neut(3) + 8._dl/15._dl*k*shear
            neutprime(2)=pirdt
            !  And for the moments
            do  l=3, EV%lmaxnrt-1
                neutprime(l)= EV%denlkt(1,L)*neut(l-1) -EV%denlkt(2,L)*neut(l+1)
            end do

            !  Truncate the hierarchy
            neutprime(EV%lmaxnrt)=k*EV%lmaxnrt/(EV%lmaxnrt-2._dl)*neut(EV%lmaxnrt-1)-  &
                (EV%lmaxnrt+3._dl)*cothxor*neut(EV%lmaxnrt)
        else
            pirdt= 8._dl/15._dl*k*shear
            neutprime(2)=pirdt
        end if

        !  Massive neutrino equations of motion and contributions to anisotropic stress.
        if (CP%Num_Nu_massive > 0) then
            do nu_i=1,CP%Nu_mass_eigenstates
                if (.not. EV%EvolveTensorMassiveNu(nu_i)) then
                    rhopi=rhopi+ grhormass(nu_i)/a2*pir !- good approx, note no rhonu weighting
                else
                    ind=EV%nu_ix(nu_i)+2

                    pinu= Nu_pi(EV, ayt, a, nu_i)
                    rhopi=rhopi+ grhormass(nu_i)/a2*pinu

                    do i=1,nqmax
                        q=nu_q(i)
                        aq=a*nu_masses(nu_i)/q
                        v=1._dl/sqrt(1._dl+aq*aq)
                        if (EV%lmaxnut>2) then
                            aytprime(ind)=-v*EV%denlkt(2,2)*ayt(ind+1)+8._dl/15._dl*k*shear
                            do l=3,EV%lmaxnut-1
                                ind=ind+1
                                aytprime(ind)=v*(EV%denlkt(1,L)*ayt(ind-1)-EV%denlkt(2,L)*ayt(ind+1))
                            end do
                            ind = ind+1
                            !Truncate moment expansion.
                            aytprime(ind)=k*v*EV%lmaxnut/(EV%lmaxnut-2._dl)*ayt(ind-1)-(EV%lmaxnut+3)*cothxor*ayt(ind)
                        else
                            aytprime(ind)=8._dl/15._dl*k*shear
                        end if
                        ind=ind+1
                    end do
                end if
            end do
        end if
    end if

    !  Get the propagation equation for the shear

    if (CP%flat) then
        aytprime(3)=-2*adotoa*shear+k*Hchi-rhopi/k
    else
        aytprime(3)=-2*adotoa*shear+k*Hchi*(1+2*CP%curv/k2)-rhopi/k
    endif

    aytprime(2)=-k*shear

    end subroutine derivst



    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    end module GaugeInterface
