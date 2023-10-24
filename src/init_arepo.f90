module init_arepo
  use omp_lib
  use datatype
  use variable
  use tools
  
  implicit none
    !--------------------------------------------------------------------------!
!     PRIVATE
!         CHARACTER(len=4), PARAMETER :: file_a  = "huhu"
    !--------------------------------------------------------------------------!
    PUBLIC :: meta_data_arepo, hydro_data_arepo
    !--------------------------------------------------------------------------!
CONTAINS 

!################################################################################################################  
  
subroutine meta_data_arepo ! #s1.1
    use omp_lib
    use datatype
    use variable
    use tools
  
    implicit none

    REAL (kind=r2):: dummy, t_snap 
    CHARACTER(LEN=22), PARAMETER  :: hydro_resultsfolder = "../../Arepo_data/data/"             ! folder with stored hdf5-files    
    
    !======================
    !0./4. Read mu_chem, gamma
    !======================
 
    call parse_real("Xmu               ", mu_chem, 1, Are_file_folder//'meta.arepo') 
    call parse_real("gamma             ", gamma,   1, Are_file_folder//'meta.arepo') 
     
    !======================
    !1. Read units
    !======================

    call parse_real("MSCALE       ", unit_density,     1, Are_file_folder//'meta.arepo') ! in g, is a mass not a density so far
    call parse_real("LSCALE       ", unit_length,      1, Are_file_folder//'meta.arepo') ! in cm                                       
    call parse_real("TIMESCALE    ", unit_time,        1, Are_file_folder//'meta.arepo') ! in seconds
    call parse_real("VSCALE       ", unit_velocity,    1, Are_file_folder//'meta.arepo') ! in cm/s
    call parse_real("TEMPSCALE    ", unit_temperature, 1, Are_file_folder//'meta.arepo') ! in Kelvin  
 
    !Convert from cgs into SI-units:
    unit_density  = 1000.0 * unit_density * unit_length**(-3.0)                  ! in kg*m**(-3) 
    unit_length   = unit_length   * 0.01                                         ! in m
    unit_velocity = unit_velocity * 0.01                                         ! in m/s
  
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
 
  call parse_real("GmX               ", dummy, 1, Are_file_folder//'meta.arepo') 
  n_xgrid = nint(dummy)
  call parse_real("GmX               ", dummy, 2, Are_file_folder//'meta.arepo') 
  n_ygrid = nint(dummy)
  call parse_real("GmX               ", dummy, 3, Are_file_folder//'meta.arepo') 
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
 
 length_x = n_xgrid    ! Length domain in x-direction, in code-units
 length_y = n_ygrid    ! Length domain in y-direction, in code-units
 length_z = n_zgrid    ! Length domain in z-direction, in code-units
  
  IF (n_zgrid == 1) THEN 
    n_para_hdf5 = 7
  ELSE
    n_para_hdf5 = 8
  END IF
  n_file_entries = n_xgrid * n_ygrid * n_zgrid * n_para_hdf5
 
 !======================
 !4. Region of interest
 !======================
 
    call parse_real("origin", origin_x, 1, Are_file_folder//'meta.arepo')    ! in pixel
    call parse_real("origin", origin_y, 2, Are_file_folder//'meta.arepo')    ! in pixel
    call parse_real("origin", origin_z, 3, Are_file_folder//'meta.arepo')    ! in pixel
    call parse_real("radius_r",   region_r, 1, Are_file_folder//'meta.arepo')! in pc
  
    ! Convert from pc into SI-units:
    region_r = region_r / m_in_pc ! in m

end subroutine meta_data_arepo  ! #s1.1

!################################################################################################################

subroutine hydro_data_arepo  ! #s1.2
  use omp_lib
  use datatype
  use variable
  
  implicit none

  INTEGER :: i_xgrid, i_ygrid, i_zgrid, n_cell_total
  INTEGER :: n_quadrant  
  
  Real, DIMENSION(:,:,:), ALLOCATABLE ::  data_hydro_code      ! Data of all h5_files (trial_cube_00xyz.hdf) Make it NOT to kind=r2.    

  CHARACTER(LEN=22), PARAMETER  :: hydro_resultsfolder = Are_file_folder//"data/"             ! folder with stored hydro-files
  CHARACTER(LEN=19) :: filename_hydro                                                         ! file name of the hydro-files  
      
  ! Differentiate 2 cases: Hydro-input from AstroBEAR (l_hydrocode = 1) or from Pencil (l_hydrocode = 2)
  
  ALLOCATE(data_hydro_code(1:n_xgrid,1:n_ygrid,1:n_zgrid))                            ! Only for this subrotine important
  ALLOCATE(data_gas(1:n_xgrid, 1:n_ygrid, 1:n_zgrid, 1:n_timepoints, 4:n_paragas)) ! Global importance, contains all hydro_data
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
      
      filename_hydro = "trial_cube_"//i_timepoints_word//".h5" ! File name
      
      !####### Density
      data_hydro_code = 0.0    

      !=================================
      call h5fopen_f (hydro_resultsfolder//filename_hydro, H5F_ACC_RDWR_F, file_id, error)   ! Open an existing file.  

      call h5dopen_f(file_id, "rho", dset_id, error)                         ! Open an existing dataset.

      call h5dread_f(dset_id, H5T_NATIVE_REAL, data_hydro_code, data_dims, error)
      
      data_gas(:,:,:, i_timepoints,7) = data_hydro_code(:,:,:) * unit_density ! in kg*m**(-3)
!       print*, data_hydro_code(4,16,32), data_gas(4,16,32, i_timepoints,7), " in kg*m**(-3)"

      call h5dclose_f(dset_id, error)       ! Close the dataset, file, Fortran interface      
      call h5fclose_f(file_id, error)       ! Close the dataset, file, Fortran interface                 
      !=================================   
      
      data_hydro_code = 0.0    
      !=================================
      call h5fopen_f (hydro_resultsfolder//filename_hydro, H5F_ACC_RDWR_F, file_id, error)   ! Open an existing file.  

      call h5dopen_f(file_id, "T", dset_id, error)                         ! Open an existing dataset.

      call h5dread_f(dset_id, H5T_NATIVE_REAL, data_hydro_code, data_dims, error)

      data_gas(:,:,:, i_timepoints,8) = data_hydro_code(:,:,:) * unit_temperature ! in K
!       print*, data_hydro_code(4,16,32), data_gas(4,16,32, i_timepoints,8), " in K

      call h5dclose_f(dset_id, error)       ! Close the dataset, file, Fortran interface      
      call h5fclose_f(file_id, error)       ! Close the dataset, file, Fortran interface                 
      !=================================         
      
      data_hydro_code = 0.0                    
      !=================================
      call h5fopen_f (hydro_resultsfolder//filename_hydro, H5F_ACC_RDWR_F, file_id, error)   ! Open an existing file.  

      call h5dopen_f(file_id, "v1", dset_id, error)                         ! Open an existing dataset.

      call h5dread_f(dset_id, H5T_NATIVE_REAL, data_hydro_code, data_dims, error)

     ! 4., 5., 6.: gas vx, vy, vz: 
      data_gas(:,:,:, i_timepoints,4) =  data_hydro_code(:,:,:) * unit_velocity   ! in m/s             
!       print*, data_hydro_code(4,16,32), data_gas(4,16,32, i_timepoints,4), " in m/s"
      
      call h5dclose_f(dset_id, error)       ! Close the dataset, file, Fortran interface      
      call h5fclose_f(file_id, error)       ! Close the dataset, file, Fortran interface                 
      !=================================

      data_hydro_code = 0.0                          
      !=================================
      call h5fopen_f (hydro_resultsfolder//filename_hydro, H5F_ACC_RDWR_F, file_id, error)   ! Open an existing file.  

      call h5dopen_f(file_id, "v2", dset_id, error)                         ! Open an existing dataset.

      call h5dread_f(dset_id, H5T_NATIVE_REAL, data_hydro_code, data_dims, error)

     ! 4., 5., 6.: gas vx, vy, vz: 
      data_gas(:,:,:, i_timepoints,5) =  data_hydro_code(:,:,:) * unit_velocity   ! in m/s    
!       print*, data_hydro_code(4,16,32), data_gas(4,16,32, i_timepoints,5), " in m/s"
      
      call h5dclose_f(dset_id, error)       ! Close the dataset, file, Fortran interface      
      call h5fclose_f(file_id, error)       ! Close the dataset, file, Fortran interface                 
      !=================================

      data_hydro_code = 0.0                          
      !=================================
      call h5fopen_f (hydro_resultsfolder//filename_hydro, H5F_ACC_RDWR_F, file_id, error)   ! Open an existing file.  

      call h5dopen_f(file_id, "v3", dset_id, error)                         ! Open an existing dataset.

      call h5dread_f(dset_id, H5T_NATIVE_REAL, data_hydro_code, data_dims, error)

     ! 4., 5., 6.: gas vx, vy, vz: 
      data_gas(:,:,:, i_timepoints,6) =  data_hydro_code(:,:,:) * unit_velocity   ! in m/s           
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
  
  DEALLOCATE(data_hydro_code)


end subroutine hydro_data_arepo  ! #s1.2

!################################################################################################################

end module init_arepo
