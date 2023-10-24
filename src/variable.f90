module variable
  use omp_lib
  use datatype
  use HDF5        ! This module contains all necessary modules
  
  implicit none

  REAL(kind=r2),parameter :: &
        pi          = 3.1415926535897932384626433832_r2, &
        kB          = 1.3806485279E-23_r2, &       ! in J/K
        m_proton    = 1.67262189821E-27, &         ! in kg
        m_amu       = 1.66053904020E-27, &         ! in kg, atomic mass unit
        m_electron  = 9.10938356E-31, &         ! in kg
        q_elementar = 1.602176620898E-19, &        ! in C
        ev_unit     = 1.602176620898E-19, &        ! 1 ev, in J
        a_Bohr      = 5.291772106712E-11, &        ! in m
        m_in_pc     = 3.2407557E-17,&              ! in pc/m
        epsilon_0   = 8.854187817E-12              ! in A**2 s**4/(kg m**3), vacuum permittivity (elektr. Feldkonstante)
        
  CHARACTER(len=1) :: i_dusttype_word
  CHARACTER(len=3) :: time_word, i_asize_word        
  CHARACTER(LEN=5) :: i_timepoints_word
  CHARACTER(LEN=7) :: theme_colour1, theme_colour2, theme_colour3  
  CHARACTER(len=9) :: t_timepoints_word

  CHARACTER(len=6), DIMENSION(:), ALLOCATABLE :: a_word  
  CHARACTER(len=8), DIMENSION(:), ALLOCATABLE :: Material_name
        
  INTEGER :: n_asize, n_dusttype, n_timepoints, n_xgrid, n_ygrid, n_zgrid, n_file_entries 
  INTEGER :: i_timepoints, i_run, Mat_A, Mat_B, ios, n_out_sizes, i_out_size
  
  INTEGER, DIMENSION(:), ALLOCATABLE :: i_min_ini, i_max_ini, DB

  REAL(kind=r2), DIMENSION(:,:), ALLOCATABLE :: count_colli  
  Real(kind=r2), DIMENSION(:,:,:,:,:), ALLOCATABLE ::  data_gas      ! Global data of the gas at x,y,z,t for parameter para
  Real(kind=r2), DIMENSION(:,:,:,:), ALLOCATABLE ::  data_grid         ! Position x,y,z in the grid i,j,k
  Real(kind=r2), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE ::  data_dust_main ! Global data of the dust at x,y,z,t, dusttype, asize for parameter dustpara
  Real(kind=r2), DIMENSION(:), ALLOCATABLE :: time                   ! Simulated i_timepoints, from 1 to n_timepoints
  Real(kind=r2), DIMENSION(:), ALLOCATABLE :: m_grainatom, Z_grainatom, U0, k_value ! Dust properties, required for the sputtering
  Real(kind=r2), DIMENSION(:), ALLOCATABLE :: v_Thres_vapo, v_Thres_shat, E_Possion, gamma_surf, Stress_Thres ! Dust properties, required for the grain-grain collisions
  Real(kind=r2), DIMENSION(:), ALLOCATABLE :: Abun_gas, Z_gas, m_gas  ! Gas properties: Abundances compared to hydrogene, atomic of number of the gas and mass of one gas atom
  Real(kind=r2), DIMENSION(:,:,:,:), ALLOCATABLE :: Ion_grade 
  Real(kind=r2), DIMENSION(:,:,:), ALLOCATABLE :: Ion_grade_average 
  
  REAL (kind=r2), DIMENSION(:,:), ALLOCATABLE :: densi_min_tot, densi_max_tot  
  Real(kind=r2), DIMENSION(:), ALLOCATABLE ::   ymin, ymax !Min and Max of the yrange grain size distributiion
  
  Real(kind=r2), DIMENSION(:), ALLOCATABLE :: a_min, a_max, gamma_asize, sigma_log_n, peak_logn_a, rho_bulk, m_proportion ! Dust distribution parameter    
  Real(kind=r2), DIMENSION(:,:), ALLOCATABLE :: a, sput_size                   ! Dust grain size distribution
  Real(kind=r2), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: v_dust, v_dust_MHD_tra, v_dust_MHD_act      ! velocity_vector of the dust grain, (vx,vy,vz), x,y,z, i_timeinter, i_dusttype, i_asize
  Real(kind=r2), DIMENSION(:), ALLOCATABLE :: M_larger_grains 
  Real(kind=r2), DIMENSION(:), ALLOCATABLE :: c_speed, Pl_Hir, Pv_Hir, s_Hir
  Real(kind=r2), DIMENSION(:), ALLOCATABLE :: rp, rp_slope
  Real(kind=r2), DIMENSION(:,:,:,:,:,:), ALLOCATABLE ::  data_dust_next    ! data of the dust at time t, entries: x,y,z, dusttype, asize, parameter dustpara  
  
  REAL(kind=r2) :: unit_length, unit_time, unit_density, unit_velocity, unit_temperature, unit_energy, unit_mag, unit_press
  REAL(kind=r2) :: length_x, length_y, length_z, clump_x, clump_y, clump_z, clump_r, origin_x, origin_y, origin_z, region_r
  REAL(kind=r2) :: Density_dustcell, V_cell, a_min_total, a_max_total, Delta_gd_0
  REAL(kind=r2) :: var, time_unit_scaler, mu_chem, gamma, delta_rad, k_gas_stick
  REAL(kind=r2) :: Delta_t, Delta_t_inter, phi_stat_help, gas_min_tot, gas_max_tot, T_min_tot,T_max_tot,B_max_tot, M_regular_gas   
  REAL(kind=r2) :: Mass_sum_reg_0, Mass_sum_reg
  
  INTEGER, parameter ::          &
         n_paragas        = 11,  &  ! number of entries in the hydro-simulations (8),       
         n_paradust       = 6,   &
         n_alpha_function = 2      ! 1 = Tielens et al. (1994), 2 = Nozawa et al. (2006)
         
  INTEGER :: n_cores, n_plasma, n_gg, n_gg_coulomb, n_sp, n_sp_coulomb, n_sput_size,  n_gas_accretion, l_hydrocode, l_MHD,&
             l_scenario, l_boundary, n_ce

         
  !--------------------------------------------------------------------------!
  
   CHARACTER(LEN=6), PARAMETER   :: AB_file_folder     = "../../"                         ! folder in which astrobear meta files are stored
   CHARACTER(LEN=18), PARAMETER  :: Pen_file_folder     = "../../Pencil_data/"            ! folder in which Pencil meta files are stored   
   CHARACTER(LEN=17), PARAMETER  :: Swi_file_folder     = "../../Swift_data/"             ! folder in which Swift meta files are stored  
   CHARACTER(LEN=17), PARAMETER  :: Are_file_folder     = "../../Arepo_data/"             ! folder in which Arepo meta files are stored     
   CHARACTER(LEN=16), PARAMETER  :: Amun_file_folder     = "../../Amun_data/"              ! folder in which Amun meta files are stored     
   CHARACTER(LEN=16), PARAMETER  :: Dens_map_folder    = "../Results/Maps/"       ! folder where to store and plot the density maps (of gas and dust)
   CHARACTER(LEN=19), PARAMETER  :: Size_dist_folder   = "../Results/GrainSD/"  ! folder where to store and plot the dust results
   CHARACTER(LEN=23), PARAMETER  :: dsetname = "level_0/data:datatype=0"                  ! Dataset name
   CHARACTER(LEN=13), PARAMETER  :: boxesname = "level_0/boxes"                           ! Boxes name   

   INTEGER(HID_T) :: file_id       ! File identifier
   INTEGER(HID_T) :: dset_id       ! Dataset identifier
 
   INTEGER     ::   error, n_para_hdf5
   INTEGER(HSIZE_T), DIMENSION(1) :: data_dims    

  !--------------------------------------------------------------------------!
  ! Set different parameter
  !--------------------------------------------------------------------------!    
  
  REAL(kind=r2), parameter ::          &
         data_accuracy = 0.0
         
  INTEGER, parameter ::       &
         n_timeinter = 10,    &
         n_v_step    = 15,    &
         n_gas_elements = 1
   
end module variable 
