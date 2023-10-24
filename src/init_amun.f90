module init_amun
  use omp_lib
  use datatype
  use variable
  use tools
  
  implicit none
    !--------------------------------------------------------------------------!
!     PRIVATE
!         CHARACTER(len=4), PARAMETER :: file_a  = "huhu"
    !--------------------------------------------------------------------------!
    PUBLIC :: meta_data_amun, hydro_data_amun
    !--------------------------------------------------------------------------!
CONTAINS 

!################################################################################################################  
  
subroutine meta_data_amun ! #s1.1
    use omp_lib
    use datatype
    use variable
    use tools
  
    implicit none

    REAL (kind=r2):: dummy, t_snap 
    CHARACTER(LEN=21), PARAMETER  :: hydro_resultsfolder = "../../Amun_data/data/"             ! folder with stored hdf5-files    
    
    !======================
    !0./4. Read mu_chem, gamma
    !======================
 
    call parse_real("Xmu               ", mu_chem, 1, Amun_file_folder//'meta.amun') 
    call parse_real("gamma             ", gamma,   1, Amun_file_folder//'meta.amun') 
     
    !======================
    !1. Read units
    !======================

    call parse_real("RSCALE       ", unit_density,     1, Amun_file_folder//'meta.amun') ! in cm**(-3)   
    call parse_real("LSCALE       ", unit_length,      1, Amun_file_folder//'meta.amun') ! in cm                                       
    call parse_real("TIMESCALE    ", unit_time,        1, Amun_file_folder//'meta.amun') ! in seconds
 
    !Convert from cgs into SI-units:
    unit_density  = mu_chem * m_amu * unit_density* 1.0E+6                       ! in kg*m**(-3) 
    unit_length   = unit_length  *   0.01                                        ! in m
 
    unit_energy   = unit_density  *unit_length**(5.0) * unit_time**(-2.0)        ! kg m**2/s**2    
    unit_press    = unit_energy/unit_length**(3.0)                               ! kg/(m *s**2 )
    
    !======================
    !2. Read number of timepoints and calculate time points.
    !======================

    open(1, file = hydro_resultsfolder//'totals.dat') 
      DO 
          read(1,*,iostat=ios)
          IF (ios/=0) EXIT
          n_timepoints = n_timepoints + 1 
      END DO
    close(1) 
    
    call parse_real("n_timepoints     ", dummy, 1, '../Process.tab')     ! number of cores to use in parallel mode
    if (nint(dummy) > 0) THEN
        n_timepoints = min(n_timepoints, nint(dummy))
    end if    

    ALLOCATE(time(1:n_timepoints))
  
    open(unit=1, file = hydro_resultsfolder//'totals.dat', action="read", status="unknown", form="formatted")    
        DO i_timepoints = 1, n_timepoints
            read(unit=1,fmt=*) time(i_timepoints)
            time(i_timepoints) = time(i_timepoints) * unit_time
        END DO
    close(unit=1)  

 !======================
 !3. Read grid extension
 !======================  
 ! Read n_xgrid, n_ygrid, n_zgrid, length_x, length_y, le   CHARACTER(LEN=15) :: filename_hydro                   ! Chombo name (Chombo00xyz.hdf)ngth_z
 
  call parse_real("GmX               ", dummy, 1, Amun_file_folder//'meta.amun') 
  n_xgrid = nint(dummy)
  call parse_real("GmX               ", dummy, 2, Amun_file_folder//'meta.amun') 
  n_ygrid = nint(dummy)
  call parse_real("GmX               ", dummy, 3, Amun_file_folder//'meta.amun') 
  n_zgrid = nint(dummy)
  
  IF (n_xgrid == 0) THEN 
    n_xgrid = 1
  END IF
  IF (n_ygrid == 0) THEN 
    n_ygrid = 1
  END IF
  IF (n_zgrid == 0) THEN 
    n_zgrid = 1
  END IF  
 
  call parse_real("GxBounds          ", length_x, 4, Amun_file_folder//'meta.amun')    ! Length domain in x-direction, in code-units
  call parse_real("GxBounds          ", length_y, 5, Amun_file_folder//'meta.amun')    ! Length domain in y-direction, in code-units 
  call parse_real("GxBounds          ", length_z, 6, Amun_file_folder//'meta.amun')    ! Length domain in z-direction, in code-units
  
  IF (n_zgrid == 1) THEN 
    n_para_hdf5 = 7
    length_z = 1
  ELSE
    n_para_hdf5 = 8
  END IF
  n_file_entries = n_xgrid * n_ygrid * n_zgrid * n_para_hdf5
  
 !======================
 !4. Region of interest
 !======================
 
    call parse_real("origin", origin_x, 1, Amun_file_folder//'meta.amun')    ! in pixel
    call parse_real("origin", origin_y, 2, Amun_file_folder//'meta.amun')    ! in pixel
    call parse_real("origin", origin_z, 3, Amun_file_folder//'meta.amun')    ! in pixel
    call parse_real("radius_r",   region_r, 1, Amun_file_folder//'meta.amun')! in pc
  
    ! Convert from pc into SI-units:
    region_r = region_r / m_in_pc ! in m

end subroutine meta_data_amun  ! #s1.1

!################################################################################################################

subroutine hydro_data_amun  ! #s1.2
  use omp_lib
  use datatype
  use variable
  
  implicit none

  INTEGER :: i_xgrid, i_ygrid, i_zgrid, n_cell_total
  INTEGER :: n_quadrant  
  
  Real, DIMENSION(:,:,:), ALLOCATABLE ::  data_hydro_code      ! Data of all h5_files (trial_cube_00xyz.hdf) Make it NOT to kind=r2.    
  Real, DIMENSION(:,:,:), ALLOCATABLE ::  speed_of_sound       ! Speed of sound in each single cell. Used for the calulation of the gas velocity

  CHARACTER(LEN=21), PARAMETER  :: hydro_resultsfolder = Amun_file_folder//"data/"             ! folder with stored hydro-files
  CHARACTER(LEN=16) :: filename_hydro                                                         ! file name of the hydro-files  
  
   CHARACTER(LEN=14), PARAMETER  :: dsetname_amun1 = "variables/dens"                  ! Dataset name  
      
  ! Differentiate 2 cases: Hydro-input from AstroBEAR (l_hydrocode = 1) or from Pencil (l_hydrocode = 2)

  ALLOCATE(data_hydro_code(1:n_xgrid,1:n_ygrid,1:n_zgrid),speed_of_sound(1:n_xgrid,1:n_ygrid,1:n_zgrid))                            ! Only for this subrotine important
  ALLOCATE(data_gas(1:n_xgrid, 1:n_ygrid, 1:n_zgrid, 1:n_timepoints, 4:9)) ! Global importance, contains all hydro_data
  ! data_gas(:,:,:,:,{1,2,3})  not used!
  ! data_gas(:,:,:,:,4) = gas vx
  ! data_gas(:,:,:,:,5) = gas vy
  ! data_gas(:,:,:,:,6) = gas vz
  ! data_gas(:,:,:,:,7) = gas density
  ! data_gas(:,:,:,:,8) = gas temperature
  ! data_gas(:,:,:,:,9) =  Bx
  ! data_gas(:,:,:,:,10) = By
  ! data_gas(:,:,:,:,11) = Bz

  !==================================================

  CALL h5open_f(error)                                                             ! Initialise FORTRAN interface.
  n_cell_total = n_xgrid * n_ygrid * n_zgrid

  Do i_timepoints = 1, n_timepoints

      write(i_timepoints_word, '(i5.5)' ) i_timepoints-1
      
      filename_hydro = "f0"//i_timepoints_word//"_00000.h5" ! File name
      
      !####### Density
      data_hydro_code = 0.0            
      !=================================
      call h5fopen_f (hydro_resultsfolder//filename_hydro, H5F_ACC_RDWR_F, file_id, error)   ! Open an existing file.  
 
      call h5dopen_f(file_id, "variables/dens", dset_id, error)                         ! Open an existing dataset.

      call h5dread_f(dset_id, H5T_NATIVE_REAL, data_hydro_code, data_dims, error)

      data_gas(:,:,:, i_timepoints,7) = data_hydro_code(:,:,:) * unit_density ! in kg*m**(-3)
!       print*, data_hydro_code(4,16,32), data_gas(4,16,32, i_timepoints,7), " in kg*m**(-3)"

      call h5dclose_f(dset_id, error)       ! Close the dataset, file, Fortran interface      
      call h5fclose_f(file_id, error)       ! Close the dataset, file, Fortran interface                 
      !=================================

      !####### Temperature      
      data_hydro_code = 0.0                          
      !=================================
      call h5fopen_f (hydro_resultsfolder//filename_hydro, H5F_ACC_RDWR_F, file_id, error)   ! Open an existing file.  
 
      call h5dopen_f(file_id, "variables/pres", dset_id, error)                         ! Open an existing dataset.

      call h5dread_f(dset_id, H5T_NATIVE_REAL, data_hydro_code, data_dims, error)

     ! 8.: gas T: 
      data_gas(:,:,:, i_timepoints,8) = ( data_hydro_code(:,:,:) * unit_press)/data_gas(:,:,:, i_timepoints,7) &
        * (gamma * mu_chem * m_amu)/kB  ! K    (Temp = Pressure/Density * (kB/(gamma * mu * amu)))  
!       print*, data_hydro_code(4,16,32),  ( data_hydro_code(4,16,32) * unit_press), data_gas(4,16,32, i_timepoints,8), " in K"
        
      call h5dclose_f(dset_id, error)       ! Close the dataset, file, Fortran interface      
      call h5fclose_f(file_id, error)       ! Close the dataset, file, Fortran interface                 
      !=================================       
      
      !####### Calculate the unit_velocity. This is the speed of sound cs= (KB* T/(mu * amu))**0.5      
      
      speed_of_sound(:,:,:) = data_gas(:,:,:, i_timepoints,8)**0.5 *  (kB/(mu_chem * m_amu))**0.5  ! in m/s   
!       print*, speed_of_sound(4,16,32), " in m/s"
      
      data_hydro_code = 0.0                    
      !=================================
      call h5fopen_f (hydro_resultsfolder//filename_hydro, H5F_ACC_RDWR_F, file_id, error)   ! Open an existing file.  
 
      call h5dopen_f(file_id, "variables/velx", dset_id, error)                         ! Open an existing dataset.

      call h5dread_f(dset_id, H5T_NATIVE_REAL, data_hydro_code, data_dims, error)

     ! 4., 5., 6.: gas vx, vy, vz: 
      data_gas(:,:,:, i_timepoints,4) =  data_hydro_code(:,:,:) * speed_of_sound(:,:,:)   ! in m/s             
!       print*, data_hydro_code(4,16,32), data_gas(4,16,32, i_timepoints,4), " in m/s"
      
      call h5dclose_f(dset_id, error)       ! Close the dataset, file, Fortran interface      
      call h5fclose_f(file_id, error)       ! Close the dataset, file, Fortran interface                 
      !=================================

      data_hydro_code = 0.0                          
      !=================================
      call h5fopen_f (hydro_resultsfolder//filename_hydro, H5F_ACC_RDWR_F, file_id, error)   ! Open an existing file.  
 
      call h5dopen_f(file_id, "variables/vely", dset_id, error)                         ! Open an existing dataset.

      call h5dread_f(dset_id, H5T_NATIVE_REAL, data_hydro_code, data_dims, error)

     ! 4., 5., 6.: gas vx, vy, vz: 
      data_gas(:,:,:, i_timepoints,5) =  data_hydro_code(:,:,:) * speed_of_sound(:,:,:)   ! in m/s    
!       print*, data_hydro_code(4,16,32), data_gas(4,16,32, i_timepoints,5), " in m/s"
      
      call h5dclose_f(dset_id, error)       ! Close the dataset, file, Fortran interface      
      call h5fclose_f(file_id, error)       ! Close the dataset, file, Fortran interface                 
      !=================================

      data_hydro_code = 0.0                          
      !=================================
      call h5fopen_f (hydro_resultsfolder//filename_hydro, H5F_ACC_RDWR_F, file_id, error)   ! Open an existing file.  
 
      call h5dopen_f(file_id, "variables/velz", dset_id, error)                         ! Open an existing dataset.

      call h5dread_f(dset_id, H5T_NATIVE_REAL, data_hydro_code, data_dims, error)

     ! 4., 5., 6.: gas vx, vy, vz: 
      data_gas(:,:,:, i_timepoints,6) =  data_hydro_code(:,:,:) * speed_of_sound(:,:,:)   ! in m/s           
!       print*, data_hydro_code(4,16,32), data_gas(4,16,32, i_timepoints,6), " in m/s"
!       print*, "########################"      
      call h5dclose_f(dset_id, error)       ! Close the dataset, file, Fortran interface      
      call h5fclose_f(file_id, error)       ! Close the dataset, file, Fortran interface                 
      !=================================   
     
     IF  (l_MHD == 1) THEN    ! magnetic field on
        ! 9., 10., 11.: gas Bx, By, Bz: 
        data_gas(:,:,:, i_timepoints,9:11) =  0.0 * unit_mag ! in Tesla
        print*,"The magnetic field strength is set to zero."
     END IF      
    !========================================================================================== 

      
      if (n_timepoints <10) THEN
        write (*,'(A,I1,A,I1,A)') '   Hydro-files loaded: ', i_timepoints,'/',n_timepoints,"." //char(27)//'[A'
      else if (n_timepoints <100) THEN
        write (*,'(A,I2,A,I2,A)') '   Hydro-files loaded: ', i_timepoints,'/',n_timepoints,"." //char(27)//'[A'
      else if (n_timepoints <1000) THEN
        write (*,'(A,I3,A,I3,A)') '   Hydro-files loaded: ', i_timepoints,'/',n_timepoints,"." //char(27)//'[A'
      else if (n_timepoints <10000) THEN
        write (*,'(A,I4,A,I4,A)') '   Hydro-files loaded: ', i_timepoints,'/',n_timepoints,"." //char(27)//'[A'
      else if (n_timepoints <100000) THEN
        write (*,'(A,I5,A,I5,A)') '   Hydro-files loaded: ', i_timepoints,'/',n_timepoints,"." //char(27)//'[A'
      else               
        print*, 'Paperboats ERROR 26 in reading in hydro data for hydro code AstroBEAR. Too many timesteps (>=100000).'
        STOP
      end if

  END DO

  CALL h5close_f(error)
  
  DEALLOCATE(data_hydro_code,speed_of_sound)

end subroutine hydro_data_amun  ! #s1.2

!################################################################################################################

end module init_amun
