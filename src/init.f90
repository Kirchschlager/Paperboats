module init
  use omp_lib
  use datatype
  use variable
  use functs
  use tools
  use init_astrobear
  use init_pencil
  use init_swift 
  use init_arepo  
  use init_amun  
  use output_data
  use output_plot   
  
  implicit none
    !--------------------------------------------------------------------------!
!     PRIVATE
!         CHARACTER(len=4), PARAMETER :: file_a  = "huhu"
    !--------------------------------------------------------------------------!
    PUBLIC :: init_paperboats, meta_data, hydro_data, dust_init, stat_calc,    &
              gas_properties
    !--------------------------------------------------------------------------!
CONTAINS 

!################################################################################################################

subroutine init_paperboats ! #s1
  use omp_lib
  use datatype
  use variable  
  use output_data  
  use output_plot
  
  implicit none
  
  print*, ''  
  print*, ' 1. Initialisation.'  
  print*, ''

  call meta_data
  call hydro_data
  call dust_init
  i_timepoints = 1
  call output_paperboats
  call output_time_def
  
end subroutine init_paperboats ! #s1

!################################################################################################################  
  
subroutine meta_data ! #s1.1
  use omp_lib
  use datatype
  use variable
  use tools
  use init_astrobear  
  use init_pencil  
  use init_swift
  use init_arepo  
  use init_amun  
  
  implicit none
  
  REAL (kind=r2):: dummy
  
  !===================
  !1. Read Process.tab
  !===================  
  ! Read which dust processes and effects are considered
  
   call parse_real("n_cores          ", dummy, 1, '../Process.tab')     ! number of cores to use in parallel mode
   n_cores         = nint(dummy)
   call parse_real("l_hydrocode      ", dummy, 1, '../Process.tab') 
   l_hydrocode     = nint(dummy)   
   call parse_real("l_MHD            ", dummy, 1, '../Process.tab')
   l_MHD           = nint(dummy)   
   call parse_real("n_plasma         ", dummy, 1, '../Process.tab')
   n_plasma        = nint(dummy)
   call parse_real("n_ce             ", dummy, 1, '../Process.tab')
   n_ce            = nint(dummy)
   call parse_real("n_gg             ", dummy, 1, '../Process.tab') 
   n_gg            = nint(dummy)   
   call parse_real("n_gg_coulomb     ", dummy, 1, '../Process.tab')  
   n_gg_coulomb    = nint(dummy)   
   call parse_real("n_sp             ", dummy, 1, '../Process.tab')  
   n_sp            = nint(dummy)     
   call parse_real("n_sp_coulomb     ", dummy, 1, '../Process.tab')
   n_sp_coulomb    = nint(dummy)      
   call parse_real("n_sput_size      ", dummy, 1, '../Process.tab') 
   n_sput_size     = nint(dummy)      
   call parse_real("n_gas_accretion  ", dummy, 1, '../Process.tab')  
   n_gas_accretion = nint(dummy)   
   call parse_real("l_scenario       ", dummy, 1, '../Process.tab')  
   l_scenario      = nint(dummy)
   call parse_real("l_boundary       ", dummy, 1, '../Process.tab')  
   l_boundary      = nint(dummy)   
   call parse_char("theme_c1         ", theme_colour1, 1, '../Process.tab')  
   call parse_char("theme_c2         ", theme_colour2, 1, '../Process.tab')  
   call parse_char("theme_c3         ", theme_colour3, 1, '../Process.tab')     
   
!============================================================================= 
   IF     (((l_scenario == 1) .and. (l_hydrocode == 1)) .or. &    ! Cloud-Crushing setup for AstroBEAR: Dust in a specified clump only
        & (((l_scenario == 2) .or. (l_scenario == 3)) .and. (l_hydrocode == 2)) .or. &        ! Dust in the full domain for Pencil            
        & (((l_scenario == 2) .or. (l_scenario == 3)) .and. (l_hydrocode == 3)) .or. &        ! Dust in the full domain for Swift           
        & (((l_scenario == 2) .or. (l_scenario == 3)) .and. (l_hydrocode == 4))        .or. & ! Dust in the full domain for Arepo          
        & (((l_scenario == 2) .or. (l_scenario == 3)) .and. (l_hydrocode == 5))        ) THEN ! Dust in the full domain for Amun           
        print*,""
   ELSE     
        print*, 'Paperboats ERROR 22: Scenario not defined for this (M)HD code.'    
        STOP
   END IF      
!============================================================================= 

  ! Differentiate 4 cases: Hydro-input from AstroBEAR (l_hydrocode = 1), Pencil (l_hydrocode = 2), Swift (l_hydrocode = 3), Arepo (l_hydrocode = 4), or Amun (l_hydrocode = 5)
  
  IF (l_hydrocode == 1) THEN

    call meta_data_astrobear
    
  ELSE IF (l_hydrocode == 2) THEN ! Second case: Hydro-input from Pencil
  
    call meta_data_pencil
    
  ELSE IF (l_hydrocode == 3) THEN ! Third case: Hydro-input from Swift
  
    call meta_data_swift
    
  ELSE IF (l_hydrocode == 4) THEN ! Fourth case: Hydro-input from Arepo
  
    call meta_data_arepo     

  ELSE IF (l_hydrocode == 5) THEN ! Fifth case: Hydro-input from Amun
  
    call meta_data_amun         
   
  ELSE 

    print*, 'Paperboats ERROR 17: Hydro code not defined.'

  END if                          ! End of differentiation of hydro-codes

end subroutine meta_data  ! #s1.1

!################################################################################################################

subroutine hydro_data  ! #s1.2
  use omp_lib
  use datatype
  use variable
  use HDF5        ! This module contains all necessary modules 
  use init_astrobear
  use init_pencil
  use init_swift  
  use init_arepo 
  use init_amun   
  
  implicit none
  
  INTEGER :: i_quadrant, n_quadrant,  n_cell_total, dx, dy, x_in_box, y_in_box, z_in_box
  INTEGER, DIMENSION(:), ALLOCATABLE :: n_cell_box
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: data_box, d_box
  INTEGER :: i_xgrid, i_ygrid, i_zgrid
  
  !====================   
  ! Define position in the grid:
  ALLOCATE(data_grid(1:n_xgrid, 1:n_ygrid, 1:n_zgrid, 1:3))  ! Global importance, contains only the position of the grid  
  Do i_xgrid = 1, n_xgrid
    Do i_ygrid = 1, n_ygrid
        Do i_zgrid = 1, n_zgrid
            data_grid(i_xgrid, i_ygrid, i_zgrid,1) = unit_length * length_x * (i_xgrid-0.5)/(n_xgrid*1.0)  ! x-position in the domain, in m, center of the cell
            data_grid(i_xgrid, i_ygrid, i_zgrid,2) = unit_length * length_y * (i_ygrid-0.5)/(n_ygrid*1.0)  ! y-position in the domain, in m, center of the cell
            data_grid(i_xgrid, i_ygrid, i_zgrid,3) = unit_length * length_z * (i_zgrid-0.5)/(n_zgrid*1.0)  ! z-position in the domain, in m, center of the cell
        END DO
    END DO
  END DO  
  !====================   

  ! Differentiate 5 cases: Hydro-input from AstroBEAR (l_hydrocode = 1), Pencil (l_hydrocode = 2), Swift (l_hydrocode = 3), Arepo (l_hydrocode = 4), or Amun (l_hydrocode = 5)

  IF (l_hydrocode == 1) THEN
  
    call hydro_data_astrobear 
  
  ELSE IF (l_hydrocode == 2) THEN ! Second case: Hydro-input from Pencil
  
    call hydro_data_pencil
    
   ELSE IF (l_hydrocode == 3) THEN ! Third case: Hydro-input from Swift
  
    call hydro_data_swift
    
   ELSE IF (l_hydrocode == 4) THEN ! Fourth case: Hydro-input from Arepo
  
    call hydro_data_arepo     
    
  ELSE IF (l_hydrocode == 5) THEN ! Fifth case: Hydro-input from Amun

    call hydro_data_amun         

  ELSE 

    print*, 'Paperboats ERROR 17: Hydro code not defined.'   

  END if                          ! End of differentiation of hydro-codes    
  
  !==========================================================================================
  ! 2.: Define the gas properties (used for sputtering)
  
  ALLOCATE(Abun_gas(1:n_gas_elements), Z_gas(1:n_gas_elements), m_gas(1:n_gas_elements), &
         Ion_grade(1:n_xgrid, 1:n_ygrid, 1:n_zgrid, 1:n_gas_elements), Ion_grade_average(1:n_xgrid, 1:n_ygrid, 1:n_zgrid))
  
  call gas_properties

end subroutine hydro_data  ! #s1.2

!################################################################################################################

subroutine dust_init  ! #s1.3
  use omp_lib
  use datatype
  use variable
  use tools
  
  implicit none
  
  REAL (kind=r2):: alpha_sum, beta_sum, N_0, dummy, mu_asize 

  INTEGER :: i, Mat_x
  INTEGER :: i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize  
 
  call parse_real("n_dusttype  ", dummy, 1, '../Dust.tab')
  n_dusttype  = nint(dummy)
  call parse_real("n_asize     ", dummy, 1, '../Dust.tab')
  n_asize     = nint(dummy)
  call parse_real("n_out_sizes ", dummy, 1, '../Dust.tab')
  n_out_sizes = nint(dummy)
  if (n_out_sizes==1) THEN
    call parse_real("i_out_size  ", dummy, 1, '../Dust.tab')
    i_out_size = nint(dummy)
    if ((i_out_size<1) .or. (i_out_size>n_asize)) THEN
        i_out_size = 1
    end if    
  end if  
  
  call parse_real("a_min_total ", a_min_total, 1, '../Dust.tab')
  call parse_real("a_max_total ", a_max_total, 1, '../Dust.tab')
  call parse_real("Delta_gd ", Delta_gd_0, 1, '../Dust.tab')  

  if (n_dusttype > 2) THEN
    print*, 'Paperboats ERROR 2 in number of grain materials: Include further dust materials.'
    STOP
  END IF   
 
  !a_min_total =     0.6   ! in nanometer, 0.6 nm (0.65 nm) ~ 100 molecules of C (Si), 1.3 nm (1.4 nm) ~ 1000 molecules of C (Si)
  !a_max_total = 10000.0   ! in nanometer
  
  !=================================================  
  ALLOCATE(Material_name(1:n_dusttype), a_min(1:n_dusttype), a_max(1:n_dusttype), gamma_asize(1:n_dusttype),&
    rho_bulk(1:n_dusttype), m_proportion(1:n_dusttype), v_Thres_vapo(1:n_dusttype), v_Thres_shat(1:n_dusttype),&
    m_grainatom(1:n_dusttype), Z_grainatom(1:n_dusttype), U0(1:n_dusttype), Stress_Thres(1:n_dusttype), &
    k_value(1:n_dusttype), sput_size(1:n_dusttype, 1:6), M_larger_grains(1:n_dusttype), &
    sigma_log_n(1:n_dusttype), peak_logn_a(1:n_dusttype), E_Possion(1:n_dusttype), gamma_surf(1:n_dusttype),&
    c_speed(1:n_dusttype), Pl_Hir(1:n_dusttype), Pv_Hir(1:n_dusttype), s_Hir(1:n_dusttype), rp(1:n_dusttype),&
    rp_slope(1:n_dusttype), DB(1:n_dusttype), densi_min_tot(1:n_dusttype, 1: n_asize),& 
    densi_max_tot(1:n_dusttype, 1: n_asize), ymin(1:n_dusttype), ymax(1:n_dusttype))
     
  ALLOCATE(a(1:n_asize, 1:n_dusttype), i_min_ini(1:n_dusttype), i_max_ini(1:n_dusttype), count_colli(1:n_timepoints, 1:4))
 
  call parse_real("Mat_A       ", dummy, 1, '../Dust.tab') 
  Mat_x = nint(dummy)
  
  call parse_real("a_min_A     ", dummy, 1, '../Dust.tab')
  a_min(1)              = dummy                            ! in nanometer
  call parse_real("a_max_A     ", dummy, 1, '../Dust.tab') 
  a_max(1)              = dummy                            ! in nanometer  ! a_max is the "center" of the maximum grain size bin
  call parse_real("DB_A        ", dummy, 1, '../Dust.tab') 
  DB(1)                 = nint(dummy)                      ! Which kind of distribution: (1) powerlaw, (2) Log-Normal  
  call parse_real("gamma_a_A   ", dummy, 1, '../Dust.tab') 
  gamma_asize(1)        = dummy                            ! if powerlaw-distribution: grain-size-exponent
  call parse_real("peak_logn_A ", dummy, 1, '../Dust.tab')
  peak_logn_a(1)        = dummy                            ! if Log-Normal-distribution: a_peak, in nm
  call parse_real("sigma_log_A ", dummy, 1, '../Dust.tab')
  sigma_log_n(1)        = dummy                            ! if Log-Normal-distribution: "sigma". The larger, the flatter.
  call parse_real("m_prop_A    ", dummy, 1, '../Dust.tab') 
  m_proportion(1)       = dummy  
  
  if (n_dusttype == 2) THEN
    call parse_real("a_min_B     ", dummy, 1, '../Dust.tab')
    a_min(2)              = dummy                            ! in nanometer
    call parse_real("a_max_B     ", dummy, 1, '../Dust.tab') 
    a_max(2)              = dummy                            ! in nanometer  ! a_max is the "center" of the maximum grain size bin
    call parse_real("DB_B        ", dummy, 1, '../Dust.tab') 
    DB(2)                 = nint(dummy)                      ! Which kind of distribution: (1) powerlaw, (2) Log-Normal  
    call parse_real("gamma_a_B   ", dummy, 1, '../Dust.tab') 
    gamma_asize(2)        = dummy                            ! if powerlaw-distribution: grain-size-exponent
    call parse_real("peak_logn_B ", dummy, 1, '../Dust.tab')
    peak_logn_a(2)        = dummy                            ! if Log-Normal-distribution: a_peak, in nm
    call parse_real("sigma_log_B ", dummy, 1, '../Dust.tab')
    sigma_log_n(2)        = dummy                            ! if Log-Normal-distribution: "sigma". The larger, the flatter.
    call parse_real("m_prop_B    ", dummy, 1, '../Dust.tab') 
    m_proportion(2)       = dummy    
  end if
  
  If ((Mat_x==1) .or. (n_dusttype ==2)) THEN
    Mat_A = 1  
    IF ((n_dusttype ==2) .and. (Mat_x==2)) THEN
        Mat_A = 2
    END IF    

    !==========================================================
    Material_name(Mat_A)      = "Silicate"

    rho_bulk(Mat_A)           =  3.3       ! in g cm**(-3),  Tielens et al. (1994), Jones et al. (1994)    
    v_Thres_vapo(Mat_A)       = 19595.3    ! in m/s   tbd
    v_Thres_shat(Mat_A)       = 2700.0     ! in m/s, Ref: Hirashita&Yan2009, Jones+1996/ 
    E_Possion(Mat_A)          =  3.919E+10 ! in J/m**3.0, Ref: Hirashita&Yan2009, Chokshi+1993    
    gamma_surf(Mat_A)         =  0.025     ! in J/m**2.0, Ref: Hirashita&Yan2009, Chokshi+1993
    c_speed(Mat_A)            = 5000.0     ! sound speed of the grain material,  Tielens et al. (1994), Jones et al. (1994)
    Pl_Hir(Mat_A)             =  3.0E+10   ! critical pressure of the grain material
    Pv_Hir(Mat_A)             =  5.4E+11   ! critical vaporisation pressure of the grain material, Jones et al. (1994)
    s_Hir(Mat_A)              =  1.23      ! Material constant,  Tielens et al. (1994), Jones et al. (1994)     
    Stress_Thres(Mat_A)       =  10.0E+8   ! in J/m**3.0, Tazaki et al. (2020)
    
    m_grainatom(Mat_A)        = 20.0       ! in atomic mass unit,  Tielens et al. (1994), Jones et al. (1994)
    Z_grainatom(Mat_A)        = 11.0       ! Tielens et al. (1994), Jones et al. (1994)
    U0(Mat_A)                 =  5.7       ! in eV, Tielens et al. (1994), Jones et al. (1994)
    k_value(Mat_A)            =  0.1       ! Nozawa et al. 2006,  Tielens et al. (1994), Jones et al. (1994)
  
    rp(Mat_A)                 =  0.7       ! Factor multiplied with r_p to obtain r_d (Bocchio+ 2016, for MgSiO_3)
    rp_slope(Mat_A)           = -3.34      ! Value to scale the penetration depth r_p, obtained from Bethe-Bloch-Equation, for MgSiO_3)
    sput_size(Mat_A,1)        =  1.0       ! MgSiO_3, Bocchio et al. (2016)
    sput_size(Mat_A,2)        =  0.5       ! MgSiO_3, Bocchio et al. (2016)
    sput_size(Mat_A,3)        =  1.0       ! MgSiO_3, Bocchio et al. (2016)
    sput_size(Mat_A,4)        =  1.8       ! MgSiO_3, Bocchio et al. (2016)
    sput_size(Mat_A,5)        =  2.1       ! MgSiO_3, Bocchio et al. (2016)
    sput_size(Mat_A,6)        =  0.76      ! MgSiO_3, Bocchio et al. (2016)
  !==========================================================  
  end if 
  
  If ((Mat_x==2) .or. (n_dusttype==2))  THEN
      Mat_B = n_dusttype    
      IF ((n_dusttype ==2) .and. (Mat_x==2)) THEN
        Mat_B = 1
      END IF    
    !==========================================================  
    Material_name(Mat_B)    = "Graphite"
    rho_bulk(Mat_B)         =  2.2       ! in g cm**(-3),  Tielens et al. (1994), Jones et al. (1994)    
    v_Thres_vapo(Mat_B)     = 22627.4    ! in m/s tbd
    v_Thres_shat(Mat_B)     = 1200.0     ! in m/s, Ref: Hirashita&Yan2009, Jones+1996/ 
    E_Possion(Mat_B)        =  6.8E+10   ! in J/m**3.0, Ref: Hirashita&Yan2009, Chokshi+1993  
    gamma_surf(Mat_B)       =  0.012     ! in J/m**2.0, Ref: Hirashita&Yan2009, Chokshi+1993
    c_speed(Mat_B)          = 1800.0     ! sound speed of the grain material,  Tielens et al. (1994), Jones et al. (1994)
    Pl_Hir(Mat_B)           =  4.0E+9    ! critical pressure of the grain material
    Pv_Hir(Mat_B)           =  5.8E+11   ! critical vaporisation pressure of the grain material, Tielens et al. (1994), Jones et al. (1994) 
    s_Hir(Mat_B)            =  1.9       ! Material constant,  Tielens et al. (1994), Jones et al. (1994)    
    Stress_Thres(Mat_B)     =  1.0E+8    ! in J/m**3.0, Tazaki et al. (2020)    
    
    m_grainatom(Mat_B)      = 12.0       ! in atomic mass unit
    Z_grainatom(Mat_B)      =  6.0       ! Bocchio et al. (2014), (Nozawa et al. 2006), Jones et al. (1994)
    U0(Mat_B)               =  4.0       ! in eV, Tielens et al. (1994), Bocchio et al. (2014), (Nozawa et al. 2006), Jones et al. (1994)
    k_value(Mat_B)          =  0.65      ! Bocchio et al. (2014), (Nozawa et al. 2006), Jones et al. (1994)
    
    rp(Mat_B)               =  0.5       ! Factor multiplied with r_p to obtain r_d (Bocchio+ 2012, for graphite)
    rp_slope(Mat_B)         = -4.73      ! Value to scale the penetration depth r_p, obtained from Bethe-Bloch-Equation, for graphite
    sput_size(Mat_B,1)      =  4.9       ! AC, Bocchio et al. (2016)
    sput_size(Mat_B,2)      =  0.55      ! AC, Bocchio et al. (2016)
    sput_size(Mat_B,3)      =  0.77      ! AC, Bocchio et al. (2016)
    sput_size(Mat_B,4)      =  4.7       ! AC, Bocchio et al. (2016)
    sput_size(Mat_B,5)      =  3.0       ! AC, Bocchio et al. (2016)
    sput_size(Mat_B,6)      =  1.2       ! AC, Bocchio et al. (2016)
  !==========================================================
  end if
   
  !===============================================================
  ! Convert, calculate Delta_rad, scale m_proportion (if necessary), check values of a_min and a_max
  a_min_total = a_min_total * 1.0E-9      ! in m
  a_max_total = a_max_total * 1.0E-9      ! in m
  a_min       = a_min * 1.0E-9            ! in m
  a_max       = a_max * 1.0E-9            ! in m
  rho_bulk    = rho_bulk * 1000.0         ! in kg * m**(-3)
  U0          = U0 * ev_unit              ! in J
  m_grainatom = m_grainatom * m_amu       ! in kg  

  dummy = real(n_asize)
  IF (n_asize > 1) THEN
    delta_rad = (a_max_total/a_min_total)**(1.0/(dummy-1.0)) ! This defines the bin-size of the dust grains
  ELSE
    delta_rad = 1.1 ! For n_asize = 1, delta_rad has no meaning. But delta_rad is required for the ghost bins in the sputtering and grain-grain collision routines. An arbitrary value of 10% is assumed in this case for the "width" of the bin. The only condition is, that it has a width, so delta_rad =/ 1.0. Otherwise, the distance between ghost bin and grainsize bin would be 0, and all the material would move instantaneously into the ghostcell.
  END IF  

  ! Calculate alpha_sum = Sum(m_proportion), and divide all 'm_proportion' by alpha_sum, so that Sum(alpha*m_proportion) = 1.0
  alpha_sum = Sum(m_proportion(:))
  m_proportion = m_proportion / alpha_sum

  ! Check if a_min_total <= a_min(i_dusttype)< a_max(i_dusttype) <= a_max_total
  Do i_dusttype = 1, n_dusttype
    IF ((a_min_total > a_min(i_dusttype) * 1.0001) .or. (a_min(i_dusttype) > a_max(i_dusttype)* 1.0001)&
        .or. (a_max(i_dusttype) > a_max_total* 1.0001)) THEN  ! Why 1.0001: The accuracy could be too low resulting in an error here.
        print*, 'Paperboats ERROR 3 in the order of dust grain sizes.'
        STOP
    END IF
  END DO
  
  !===============================================================
  ! Define the dust grain sizes:
  
  Do i_dusttype = 1, n_dusttype
        Do i_asize = 1, n_asize
            a(i_asize, i_dusttype) = a_min_total * delta_rad**(1.0*i_asize-1.0) ! This defines the grain sizes, in m.
        END DO
  END DO
  
  !===============================================================   
  
  !=====================
  ! However, initial dust grains are not from 1 to n_asize, but from  1<=i_min_ini to  i_max_ini<=n_asize
  ! Find i_min_ini with a(i_min_ini, i_dusttype) <= a_min(i_dusttype) < a(i_min_ini+1, i_dusttype)                                              
  i_min_ini = 0
  i_max_ini = 0
  Do i_dusttype = 1, n_dusttype
    Do i = 1, n_asize
        if (a_min(i_dusttype) > a(i, i_dusttype)*0.9999) then ! Why 0.9999: Sometimes, the accuracy is not enough, and two values (which should be equal) are not equal.
            i_min_ini(i_dusttype) = i
        else
            exit
        End if
    END DO
    ! a_min(i_dusttype) is now between a(i_min_ini, i_dusttype) and a(i_min_ini+1, i_dusttype)

    ! Find i_max_ini with a(i_max_ini-1, i_dusttype) <= a_max(i_dusttype) < a(i_max_ini, i_dusttype)   
    Do i = 1, n_asize
        if (a_max(i_dusttype) > a(i, i_dusttype)*0.9999) then ! Why 0.9999 Sometimes, the accuracy is not enough, and two values (which should be equal) are not equal.
            i_max_ini(i_dusttype) = i
        else
            exit
        end if
    END DO
  END DO 
  ! a_max(i_dusttype) is now between a(i_max_ini-1, i_dusttype) and a(i_max_ini, i_dusttype)
  !===================== 

  !===============================================================  
  
  V_cell = length_x * length_y * length_z * unit_length**(3.0)/(n_xgrid * n_ygrid * n_zgrid * 1.0)

  ALLOCATE(data_dust_main(1:n_xgrid, 1:n_ygrid, 1:n_zgrid, 1:3, 1:n_dusttype, 0:n_asize, 1:n_paradust)) ! Global importance !  
  
  ! Set all values to zero, for each timepoint and each cell 
  data_dust_main = 0.0_r2

  ! Set the dust in the clump (and no dust outside the clump), calculate the number of dust particles in the cells.
  Do i_xgrid = 1, n_xgrid
      Do i_ygrid = 1, n_ygrid
          Do i_zgrid = 1, n_zgrid
              ! Check if the cell_center of i_xgrid, i_ygrid, i_zgrid is within the clump
              IF (   ((l_scenario == 1) .and. (l_hydrocode == 1) .and. & ! Cloud-Crushing problem              
                    (( (data_grid(i_xgrid, i_ygrid, i_zgrid, 1) - clump_x)**2.0  +  &
                    (data_grid(i_xgrid, i_ygrid, i_zgrid, 2) - clump_y)**2.0  +  &
                    (data_grid(i_xgrid, i_ygrid, i_zgrid, 3) - clump_z)**2.0) <= &
                     clump_r**2.0 ))   .or. &
                     (((l_hydrocode == 2) .or. (l_hydrocode == 3) .or. (l_hydrocode == 4).or. (l_hydrocode == 5)) &
                     & .and. ((l_scenario == 2) .or. (l_scenario == 3)))  ) THEN
                     
                    ! Put dust into this cell
                    Density_dustcell = data_gas(i_xgrid, i_ygrid, i_zgrid, 1, 7) / Delta_gd_0 ! in kg/ m**3     This is the mass in this cell divided by V_cell, so actually the dust density.
 
                    Do i_dusttype = 1, n_dusttype
                    
                        ! Calculate the scaling factor of grain size distribution N0, N(a) = N0 * a**(-gamma_asize), for each material i_dusttype
                        beta_sum = 0.0
                        Do i_asize = i_min_ini(i_dusttype), i_max_ini(i_dusttype)
                            IF (DB(i_dusttype) == 1) THEN !Powerlaw
                                beta_sum = beta_sum + ((a(i_asize, i_dusttype))**(3.0-gamma_asize(i_dusttype)))&
                                *delta_rad**(real(i_asize)-1.0)           ! This defines the grain size distribution
                            ELSE !Log-normal distribution
                                mu_asize = sigma_log_n(i_dusttype)**2.0 + log(peak_logn_a(i_dusttype)*1.0E-9)
                                beta_sum = beta_sum + ((a(i_asize, i_dusttype))**2.0) *&
                                    exp(-((log(a(i_asize, i_dusttype))-mu_asize)**2.0)/&
                                    (2.0*sigma_log_n(i_dusttype)**2.0)) * delta_rad**(real(i_asize)-1.0)
                            END IF        
                        END DO
                        
                        N_0 = m_proportion(i_dusttype) * Density_dustcell *& 
                               (4.0/3.0 * pi * rho_bulk(i_dusttype) * beta_sum)**(-1.0)                       
                        
                        Do i_asize = i_min_ini(i_dusttype), i_max_ini(i_dusttype)

                            !==================================================================
                            ! 1., 2., 3. Initial velocity vx, vy, vz of the dust. Assumption: v_gas = v_dust at the beginning.                       
                            data_dust_main(i_xgrid, i_ygrid, i_zgrid, 1, i_asize, i_dusttype, 1) = &
                                  data_gas(i_xgrid, i_ygrid, i_zgrid, 1, 4) 
                            data_dust_main(i_xgrid, i_ygrid, i_zgrid, 1, i_asize, i_dusttype, 2) = &
                                  data_gas(i_xgrid, i_ygrid, i_zgrid, 1, 5)
                            data_dust_main(i_xgrid, i_ygrid, i_zgrid, 1, i_asize, i_dusttype, 3) = &
                                  data_gas(i_xgrid, i_ygrid, i_zgrid, 1, 6)                          
                        
                            !==================================================================
                            ! 4. Number of dust particles per Volume (m**3) in cell x,y,z, at t=1, with size i_asize and dust material i_dusttype
                            IF (DB(i_dusttype) == 1) THEN !Powerlaw
                                data_dust_main(i_xgrid, i_ygrid, i_zgrid, 1, i_dusttype, i_asize, 4) = &
                                    N_0 * (a(i_asize, i_dusttype))**(-gamma_asize(i_dusttype))*&
                                    delta_rad**(real(i_asize)-1.0)
                            ELSE !Log-normal distribution
                                mu_asize = sigma_log_n(i_dusttype)**2.0 + log(peak_logn_a(i_dusttype)*1.0E-9)
                                data_dust_main(i_xgrid, i_ygrid, i_zgrid, 1, i_dusttype, i_asize, 4) = &
                                    N_0 * (a(i_asize, i_dusttype))**(-1.0)*&
                                    exp(-((log(a(i_asize, i_dusttype))-mu_asize)**2.0)/&
                                    (2.0*sigma_log_n(i_dusttype)**2.0)) * delta_rad**(real(i_asize)-1.0)   
                            END IF                               
                            
                            !==================================================================
                            ! 5. Mass density of the dust particles (in this cell, with this dust parameter)    
                            data_dust_main(i_xgrid, i_ygrid, i_zgrid, 1, i_dusttype, i_asize, 5) =   &
                                data_dust_main(i_xgrid, i_ygrid, i_zgrid, 1, i_dusttype, i_asize, 4) * &
                                    4.0/3.0 * pi * (a(i_asize,  i_dusttype))**3.0 * rho_bulk(i_dusttype)
                                    
                            !==================================================================
                            ! 6. Charge (in e) of the dust grains in cell x,y,z, at t=1, with size i_asize and dust material i_dusttype 
                            data_dust_main(i_xgrid, i_ygrid, i_zgrid, 1, i_dusttype, i_asize, 6) = 0.0                           
                            
                        END DO
                    END DO
              END IF
              ! All other cells contain no dust at the beginning. The velocity of all dust particles is zero at t = 0
              
          END DO
      END DO      
  END DO

  !==============================================================
  ! Initialise the array of the dust grains with a > a_max_total. 
  ! M_larger_grains(i_dusttype) is at the beginning always zero.      
  
  M_larger_grains(:) = 0.0_r2

  !==============================================================
  ! Initialise the regular gas mass counter (for trapping and sticking). 
  M_regular_gas = 0.0_r2
  
  !==============================================================
  ! Initialise the counter for the collisions       
  count_colli = 0.0_r2
  
  !==============================================================
  ! Start the lift. Timepoint 1 is the initial value and is always stored, timestep 2 is 'now', timestep '3' is next

  data_dust_main(:, :, :, 2, :, :, :) = data_dust_main(:, :, :, 1, :, :, :)

  !==============================================================  
  ! Calculate one single time phi_stat_help (required for charge calculation)
  if (n_plasma == 2) THEN
    call stat_calc
  END IF
  
  DEALLOCATE (a_min, a_max, m_proportion)

end subroutine dust_init ! #s1.3  

!################################################################################################################

subroutine stat_calc ! #s1.3.1 
  use omp_lib
  use datatype
  use variable
  
  implicit none
  
  REAL (kind=r2):: res, y_l, y_m, lx, mx, rx
  
    res = fun(0.0_r2)
    if (res < 0.0) THEN
            print*, "Paperboats ERROR 13, calculation of phi_stat not possible, mu_chem too small."
            stop
    end if    

    Iteration1 : Do
                    lx = -10.0
                    rx = 0.0
                    
                    res = fun(lx)
                    IF (res > 0.0) THEN 
                        lx = lx * 10.0
                        if (lx < -10.0**6.0) THEN
                            print*, "Paperboats ERROR 14, calculation of phi_stat takes too long!"
                            stop
                        end if
                    else
                        exit
                    END IF
                 END DO Iteration1
                 
    Iteration2 : Do                 
                    mx = (lx + rx)*0.5
                    
                    y_l  = fun(lx)
                    y_m  = fun(mx)
     
                    IF (y_l * y_m < 0.0) THEN
                        rx = mx  
                        mx = (lx + rx)*0.5
                    ELSE 
                        lx = mx  
                        mx = (lx + rx)*0.5             
                    END IF
    
                    IF (abs(lx-rx) < 0.001) THEN
                        phi_stat_help = mx
                        EXIT
                    END IF
   
                 END DO Iteration2       

end subroutine stat_calc ! #s1.3.1

!################################################################################################################  
  
subroutine gas_properties
    use omp_lib
    use datatype
    use variable

    implicit none

    REAL (kind=r2) :: Abun_sum

    !==================================================
    !Gas material:
    
    If (nint(mu_chem)==1) THEN
        Abun_gas(1) = 100.0   ! #1 Hydrogen    
        m_gas(1)    =   1.0   ! #1 Hydrogen        
        Z_gas(1)    =   1.0   ! #1 Hydrogen        
    ELSE if (nint(mu_chem)==16) THEN   
        Abun_gas(1) = 100.0   ! #5 Oxygen    
        m_gas(1)    =  16.0   ! #5 Oxygen        
        Z_gas(1)    =   8.0   ! #5 Oxygen        
    ELSE
        print*, "Paperboats ERROR 23, gas is not hydrogen or oxygen."
        stop
    END IF
 
!   m_gas(2) = 4.002602 ! #2 Helium
!   m_gas(3) = 12.0107  ! #3 Carbon
!   m_gas(4) = 14.0067  ! #4 Nitrogen
!   m_gas(5) = 1.00794  ! #1 Hydrogene     
!   m_gas(6) = 20.1797  ! #6 Neon
!   m_gas(7) = 24.305   ! #7 Magnesium
!   m_gas(8) = 28.0855  ! #8 Silicon
!   m_gas(9) = 32.065   ! #9 Sulfur
!   m_gas(10) = 55.845  ! #10 Iron

!   Z_gas(2) =  2       ! #2 Helium
!   Z_gas(3) =  6       ! #3 Carbon
!   Z_gas(4) =  7       ! #4 Nitrogen
!   Z_gas(5) =  1       ! #1 Hydrogene
!   Z_gas(6) = 10       ! #6 Neon
!   Z_gas(7) = 12       ! #7 Magnesium
!   Z_gas(8) = 14       ! #8 Silicon
!   Z_gas(9) = 16       ! #9 Sulfur
!   Z_gas(10) = 26      ! #10 Iron

    ! Calculate Abun_sum = Sum(Abun), and divide all 'Abun' by Abun_sum, so that Sum(Abun) = 1.0
    Abun_sum = Sum(Abun_gas(:))
    Abun_gas = Abun_gas / Abun_sum ! 100% =~ 1.0
    
    m_gas = m_gas * m_amu     ! in kg
    
    k_gas_stick = 1.0
    
end subroutine gas_properties  

!################################################################################################################

end module init
