module tjlf_max_dimensions

export nb,nbm,nsm,nt0,nkym,maxmodes,max_ELITE,max_fourier,ms,max_plot

#
# global dimensions for shared arrays
#

# not sure which formatting I like better
const global nb = 32 :: Int
const global nvm::Int = 2*32-1
const global nsm::Int = 12
const global nt0::Int = 40
const global nkym::Int = 512
const global maxmodes::Int = 16
const global max_ELITE::Int = 700
const global max_fourier::Int = 24
const global ms::Int = 128  # ms needs to be divisible by 8
const global max_plot::Int = 18*ms/8+1

end




# module tglf_dimensions

# using tglf_max_dimensions

# #
# # dimensions determined by inputs
# #
# !
# INTEGER nx,nbasis,nbasis_max,ns0,ns,nky,iur,nroot

# end




module tjlf_global
using ..tjlf_max_dimensions
#
#  global controls for the tglf driver routine
#

# global constants
xi :: ComplexF64=0.0+0.0im

# internal flow control switches
new_eikonal_in = true
new_start = true
new_matrix = true
new_geometry = true
new_width = true
new_kyspectrum = true
gauher_uncalled = true
gauss_hermite_uncalled = true
eikonal_unsaved= true
igeo::Int = 1
use_default_species = true
trace_path = zeros(Int,8)
# LOGICAL,EXTERNAL :: tglf_isnan
# LOGICAL,EXTERNAL :: tglf_isinf

# ! input units switch
# CHARACTER (len=8) :: units_in = 'GYRO'

# ! Input Gaussian width
# REAL :: width_in=1.65
# REAL :: width_min_in=0.3
# INTEGER :: nwidth_in=21
# LOGICAL :: find_width_in=.TRUE.

# ! Input kys
# REAL :: ky_in=0.3

# ! Input species 
# INTEGER :: ns_in=2, nstotal_in = 2
# REAL,DIMENSION(nsm) :: mass_in
# REAL,DIMENSION(nsm) :: zs_in

# ! input switches
# LOGICAL :: iflux_in=.TRUE.
# LOGICAL :: use_bper_in=.FALSE.
# LOGICAL :: use_bpar_in=.FALSE.
# LOGICAL :: use_mhd_rule_in=.TRUE.
# LOGICAL :: use_bisection_in=.TRUE.
# LOGICAL :: use_inboard_detrapped_in=.FALSE.
# LOGICAL :: use_ave_ion_grid_in=.FALSE.
# INTEGER :: ibranch_in=-1
# INTEGER :: nmodes_in=2
# INTEGER :: nbasis_max_in=4
# INTEGER :: nbasis_min_in=2
# INTEGER :: nxgrid_in=16
# INTEGER :: nky_in=12
# INTEGER :: mainion=2

# ! input rare switches
# REAL :: theta_trapped_in=0.7
# REAL :: wdia_trapped_in=0.0
# REAL :: park_in=1.0
# REAL :: ghat_in=1.0
# REAL :: gchat_in=1.0
# REAL :: wd_zero_in=0.1
# REAL :: Linsker_factor_in=0.0
# REAL :: gradB_factor_in=0.0
# REAL :: filter_in=2.0
# REAL :: damp_psi_in = 0.0
# REAL :: damp_sig_in = 0.0

# Input model paramaters
adiabatic_elec_in=false
alpha_e_in::Real = 1.0
alpha_p_in::Real = 1.0
alpha_mach_in::Real = 0.0
alpha_quench_in::Real = 0.0
alpha_zf_in ::Real = 1.0
xnu_factor_in::Real = 1.0
debye_factor_in::Real = 1.0
etg_factor_in::Real = 1.25
rlnp_cutoff_in::Real = 18.0
sat_rule_in::Int =  0
kygrid_model_in::Int = 1
xnu_model_in::Int = 2
vpar_model_in::Int = 0
vpar_shear_model_in::Int = 1    
# REAL :: alpha_n_in =0.0  !not used
# REAL :: alpha_t_in =0.0  !not used

# ! Input signs
# REAL :: sign_Bt_in = 1.0
# REAL :: sign_It_in = 1.0
# ! Input field gradients
# REAL,DIMENSION(nsm) :: rlns_in=0.0
# REAL,DIMENSION(nsm) :: rlts_in=0.0
# REAL,DIMENSION(nsm) :: vpar_shear_in=0.0
# REAL :: vexb_shear_in=0.0
# ! Input profile shear
# REAL,DIMENSION(nsm) :: vns_shear_in=0.0
# REAL,DIMENSION(nsm) :: vts_shear_in=0.0
# ! Input field averages
# REAL,DIMENSION(nsm) :: as_in=0.0
# REAL,DIMENSION(nsm) :: taus_in=0.0
# REAL,DIMENSION(nsm) :: vpar_in=0.0
# REAL :: vexb_in = 0.0
# REAL :: betae_in=0.0
# REAL :: xnue_in=0.0
# REAL :: zeff_in=1.0
# REAL :: debye_in=0.0
# ! Shifted circle (s-alpha) inputs
# REAL :: rmin_sa=0.5
# REAL :: rmaj_sa=3.0
# REAL :: q_sa=2.0
# REAL :: shat_sa=1.0
# REAL :: alpha_sa=0.0
# REAL :: xwell_sa=0.0
# REAL :: theta0_sa=0.0
# ! Shifted circle flags
# INTEGER :: b_model_sa=1
# INTEGER :: ft_model_sa=1
# ! Miller inputs
# REAL :: rmin_loc=0.5
# REAL :: rmaj_loc=3.0
# REAL :: zmaj_loc=0.0
# REAL :: q_loc=2.0
# !      REAL :: shat_loc=1.0
# REAL :: dlnpdr_loc=0.0
# REAL :: drmindx_loc=1.0
# REAL :: drmajdx_loc=0.0
# REAL :: dzmajdx_loc=0.0
# REAL :: kappa_loc=1.0
# REAL :: s_kappa_loc=0.0
# REAL :: delta_loc=0.0
# REAL :: s_delta_loc=0.0
# REAL :: zeta_loc=0.0
# REAL :: s_zeta_loc=0.0
# REAL :: p_prime_loc=0.0
# REAL :: q_prime_loc=16.0
# REAL :: beta_loc = 0.0
# REAL :: kx0_loc = 0.0
# ! Fourier geometry inputs
# INTEGER :: nfourier_in        = 16
# REAL :: q_fourier_in          = 2.0
# REAL :: q_prime_fourier_in    = 16.0
# REAL :: p_prime_fourier_in    = 0.0
# REAL :: fourier_in(8,0:max_fourier) = 0.0
# ! ELITE geometry inputs
# INTEGER :: n_ELITE
# REAL :: R_ELITE(0:max_ELITE)
# REAL :: Z_ELITE(0:max_ELITE)
# REAL :: Bp_ELITE(0:max_ELITE)
# REAL :: q_ELITE
# REAL :: q_prime_ELITE
# REAL :: p_prime_ELITE
# ! global variables
# REAL :: ft=0.5
# REAL :: ft_min=0.01
# REAL :: ft_test
# REAL :: modB_test
# REAL :: modB_min
# REAL :: ky=0.3
# REAL :: R_unit=3.0
# REAL :: q_unit=2.0
# REAL :: B_unit=1.0
# REAL :: Rmaj_input = 3.0
# REAL :: q_in = 2.0
# REAL :: rmin_input=0.5
# REAL,DIMENSION(maxmodes) :: gamma_reference_kx0 =0.0
# REAL,DIMENSION(maxmodes) :: freq_reference_kx0 =0.0
# REAL,DIMENSION(2,nkym,maxmodes) :: eigenvalue_first_pass =0.0
# REAL :: pol=1.0
# REAL :: U0=0.0
# REAL :: kx0=0.0
# REAL :: kx0_e=0.0
# REAL :: kx0_p = 0.0
# REAL :: midplane_shear=1.0
# REAL :: kx0_factor=1.0
# REAL :: rho_ion=1.0
# REAL :: rho_e=1.0
# ! output
# COMPLEX,DIMENSION(3,nb) :: field_weight_QL_out=0.0
# COMPLEX,DIMENSION(maxmodes,3,nb) :: field_weight_out=0.0
# COMPLEX,DIMENSION(maxmodes,3,max_plot) :: plot_field_out=0.0
# REAL,DIMENSION(max_plot) :: plot_angle_out=0.0
# REAL,DIMENSION(maxmodes,nsm,3) :: particle_QL_out=0.0
# REAL,DIMENSION(maxmodes,nsm,3) :: energy_QL_out=0.0
# REAL,DIMENSION(maxmodes,nsm,3) :: stress_par_QL_out=0.0
# REAL,DIMENSION(maxmodes,nsm,3) :: stress_tor_QL_out=0.0
# REAL,DIMENSION(maxmodes,nsm,3) :: exchange_QL_out=0.0
# REAL,DIMENSION(maxmodes,nsm) :: N_QL_out=0.0,T_QL_out=0.0
# REAL,DIMENSION(maxmodes,nsm) :: U_QL_out=0.0,Q_QL_out=0.0
# REAL,DIMENSION(maxmodes,nsm) :: n_bar_out=0.0,t_bar_out=0.0
# REAL,DIMENSION(maxmodes,nsm) :: u_bar_out=0.0,q_bar_out=0.0
# REAL,DIMENSION(maxmodes,nsm) :: Ns_Ts_phase_out=0.0
# REAL,DIMENSION(nsm,3) :: particle_flux_out=0.0,energy_flux_out=0.0
# REAL,DIMENSION(nsm,3) :: exchange_out=0.0
# REAL,DIMENSION(nsm,3) :: stress_par_out=0.0,stress_tor_out=0.0
# REAL,DIMENSION(maxmodes) :: gamma_out=0.0,freq_out=0.0
# REAL,DIMENSION(maxmodes) :: v_QL_out=0.0,a_par_QL_out=0.0,b_par_QL_out=0.0
# REAL,DIMENSION(maxmodes) :: phi_bar_out=0.0,v_bar_out=0.0
# REAL,DIMENSION(maxmodes) :: a_par_bar_out=0.0,b_par_bar_out=0.0
# REAL,DIMENSION(maxmodes) :: wd_bar_out=0.0,b0_bar_out=0.0
# REAL,DIMENSION(maxmodes) :: ne_te_phase_out=0.0
# REAL,DIMENSION(maxmodes) :: kx_bar_out=0.0,kpar_bar_out=0.0
# REAL,DIMENSION(maxmodes) :: modB_bar_out=0.0
# REAL,DIMENSION(nsm) :: n_bar_sum_out=0.0,t_bar_sum_out=0.0
# REAL,DIMENSION(nsm) :: q_low_out=0.0
# REAL,DIMENSION(4,nkym,maxmodes) :: field_spectrum_out=0.0
# REAL,DIMENSION(4,nkym,maxmodes) :: QL_field_spectrum_out=0.0
# REAL,DIMENSION(4,nsm,nkym,maxmodes) :: intensity_spectrum_out=0.0
# REAL,DIMENSION(4,nsm,nkym,maxmodes) :: QL_intensity_spectrum_out=0.0
# REAL,DIMENSION(5,nsm,3,nkym,maxmodes) :: flux_spectrum_out=0.0
# REAL,DIMENSION(5,nsm,3,nkym,maxmodes) :: QL_flux_spectrum_out=0.0
# REAL,DIMENSION(2,nkym,maxmodes) :: eigenvalue_spectrum_out=0.0
# REAl,DIMENSION(nkym,maxmodes) :: ne_te_phase_spectrum_out=0.0
# REAl,DIMENSION(nsm,nkym,maxmodes) :: nsts_phase_spectrum_out=0.0
# REAL,DIMENSION(nkym) :: spectral_shift_out=0.0
# REAL,DIMENSION(nkym) :: ave_p0_spectrum_out=0.0
# REAL,DIMENSION(nkym) :: width_out=0.0
# REAL :: Vzf_out = 0.0
# REAL :: kymax_out = 0.0
# REAL :: phi_bar_sum_out=0.0
# REAL :: v_bar_sum_out=0.0
# REAL :: gamma_nb_min_out=0.0
# REAL :: B2_ave_out=1.0
# REAL :: R2_ave_out=1.0
# REAL :: B_ave_out=1.0
# REAL :: Bt_ave_out=1.0
# REAL :: Bp0_out = 1.0
# REAL :: RBt_ave_out=1.0
# REAL :: Grad_r_ave_out=1.0
# REAL :: grad_r0_out=1.0
# REAL :: B_geo0_out = 1.0
# REAL :: Bt0_out = 1.0
# REAL :: SAT_geo0_out=1.0
# REAL :: SAT_geo1_out=1.0
# REAL :: SAT_geo2_out=1.0
# REAL :: kx_geo0_out=1.0
# REAL :: DM_out = 0.25
# REAL :: DR_out = 0.0
# REAL :: Bref_out = 1.0
# REAL :: ave_p0_out = 1.0
# INTEGER :: nmodes_out = 2
# INTEGER :: nfields_out = 1
# INTEGER :: jmax_out = 0
# character (len=80) :: error_msg='null' 
# ! NN activation parameters (thresholds)  
# REAL :: nn_max_error_in = -1.0
# LOGICAL :: valid_nn = .FALSE.


# !
end