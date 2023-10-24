module init_astrobear
  use omp_lib
  use datatype
  use variable
  use tools
  
  implicit none
    !--------------------------------------------------------------------------!
!     PRIVATE
!         CHARACTER(len=4), PARAMETER :: file_a  = "huhu"
    !--------------------------------------------------------------------------!
    PUBLIC :: meta_data_astrobear, hydro_data_astrobear
    !--------------------------------------------------------------------------!
CONTAINS 

!################################################################################################################  
  
subroutine meta_data_astrobear ! #s1.1
    use omp_lib
    use datatype
    use variable
    use tools
  
    implicit none
  
    REAL (kind=r2):: dummy
    CHARACTER(LEN=10), PARAMETER  :: hydro_resultsfolder = "../../out/"             ! folder with stored hdf5-files
  
    !===================
    !1. Read scales.data
    !===================  
    ! Read units
  
    call parse_real(" RSCALE  ", unit_density, 1, AB_file_folder//'scales.data')             ! in g*cm**(-3)   
    call parse_real(" LSCALE  ", unit_length, 1, AB_file_folder//'scales.data')              ! in cm                                       
    call parse_real(" TIMESCALE       ", unit_time, 1, AB_file_folder//'scales.data')        ! in seconds
    call parse_real(" TEMPSCALE       ", unit_temperature, 1, AB_file_folder//'scales.data') ! in Kelvin  
    call parse_real(" BSCALE  ", unit_mag, 1, AB_file_folder//'scales.data')                 ! in Gauss
  
    ! Convert from cgs into SI-units:
    unit_density = unit_density * 1000.0                                        ! in kg*m**(-3) 
    unit_length  = unit_length  *   0.01                                        ! in m
    unit_mag     = unit_mag * 1.0E-4                                            ! in Tesla
  
    unit_velocity = unit_length/unit_time
    unit_energy = unit_density  *unit_length**(5.0) * unit_time**(-2.0)         ! kg m**2/s**2    
  
    !===================
    !2. Read totals.data
    !===================
    ! Read number of timepoints, and time points.
  
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
      read(unit=1,fmt=*)
      time(1) = 0.0
      DO i_timepoints = 2, n_timepoints
          read(unit=1,fmt=*) time(i_timepoints)
          time(i_timepoints) = time(i_timepoints) * unit_time
      END DO
    close(unit=1)
  
    !===================
    !3. Read global.data
    !===================  
    ! Read n_xgrid, n_ygrid, n_zgrid, length_x, length_y, length_z
  
    call parse_real("GmX ", dummy, 1, AB_file_folder//'global.data')
    n_xgrid = nint(dummy)
    call parse_real("GmX ", dummy, 2, AB_file_folder//'global.data')
    n_ygrid = nint(dummy)
    call parse_real("GmX ", dummy, 3, AB_file_folder//'global.data')
    n_zgrid = nint(dummy)
  
    call parse_real("GxBounds ", length_x, 4, AB_file_folder//'global.data')
    call parse_real("GxBounds ", length_y, 5, AB_file_folder//'global.data')
    call parse_real("GxBounds ", length_z, 6, AB_file_folder//'global.data')
    
  IF (n_xgrid == 0) THEN 
    n_xgrid = 1
  END IF
  IF (n_ygrid == 0) THEN 
    n_ygrid = 1
  END IF
  IF (n_zgrid == 0) THEN 
    n_zgrid = 1
  END IF    
    
    IF (n_zgrid == 1) THEN 
        length_z = 1
    END IF
  
    ! Magnetic field: off (0) or on (1)
    IF (l_MHD == 0) THEN  
        IF (n_zgrid == 1) THEN 
            n_para_hdf5 = 7
        ELSE
            n_para_hdf5 = 8
        END IF
    ELSE
        n_para_hdf5 = 11
    END IF
    n_file_entries = n_xgrid * n_ygrid * n_zgrid * n_para_hdf5  
  
    !====================
    !4. Read problem.data
    !==================== 
    ! Read clump_x, clump_y, clump_z, clump_r
  
    call parse_real("position", clump_x, 1, AB_file_folder//'problem.data')
    call parse_real("position", clump_y, 2, AB_file_folder//'problem.data')
    call parse_real("position", clump_z, 3, AB_file_folder//'problem.data')
    call parse_real("radius",   clump_r, 1, AB_file_folder//'problem.data')
  
    ! Convert from code-units into SI-units:
    clump_x = clump_x * unit_length
    clump_y = clump_y * unit_length
    clump_z = clump_z * unit_length
    clump_r = clump_r * unit_length
  
    !====================
    !5. Read physics.data
    !==================== 
    ! Read mu_chem, gamma
  
    call parse_real("Xmu", mu_chem, 1, AB_file_folder//'physics.data')  
    call parse_real("gamma", gamma, 1, AB_file_folder//'physics.data')
  
end subroutine meta_data_astrobear  ! #s1.1

!################################################################################################################

subroutine hydro_data_astrobear  ! #s1.2
  use omp_lib
  use datatype
  use variable
  use HDF5        ! This module contains all necessary modules
  
  implicit none
  
  INTEGER :: i_quadrant, n_quadrant,  n_cell_total, dx, dy, x_in_box, y_in_box, z_in_box
  INTEGER, DIMENSION(:), ALLOCATABLE :: n_cell_box
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: data_box, d_box
  INTEGER :: i_xgrid, i_ygrid, i_zgrid
  
  Real, DIMENSION(:,:), ALLOCATABLE ::  data_hydro_code      ! Data of all hdf5_files (Chombo00xyz.hdf) Make it NOT to kind=r2.  
  
  CHARACTER(LEN=10), PARAMETER  :: hydro_resultsfolder = "../../out/"             ! folder with stored hdf5-files
  CHARACTER(LEN=15) :: filename_hydro                                             ! Chombo name (Chombo00xyz.hdf)  
    
  ! Differentiate 2 cases: Hydro-input from AstroBEAR (l_hydrocode = 1) or from Pencil (l_hydrocode = 2)
  
  ALLOCATE(data_hydro_code(1:n_file_entries, 1:n_timepoints))                            ! Only for this subrotine important
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
  ! Read the quadrants

  CALL h5open_f(error)                                                                         ! Initialise FORTRAN interface.
  call h5fopen_f (hydro_resultsfolder//"chombo00000.hdf", H5F_ACC_RDONLY_F, file_id, error)     ! Open an existing file.
  call h5dopen_f(file_id, "domain_level/childbox_count", dset_id, error)                       ! Open an existing dataset. 
     
  ! Read the number of quadrants.
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, n_quadrant, data_dims, error) 
 
  call h5dclose_f(dset_id, error)  ! Close the dataset
  call h5fclose_f(file_id, error)  ! Close the file
  call h5close_f(error)            ! Close the Fortran interface
   
  ALLOCATE(data_box(1:n_quadrant, 1:6), d_box(1:n_quadrant,1:3), n_cell_box(1:n_quadrant))
  
  open(unit=1, file = hydro_resultsfolder//'child_boxes.txt', action="read", status="unknown", form="formatted")    
      DO i_quadrant = 1, n_quadrant
          read(unit=1,fmt=*) data_box(i_quadrant, 1), data_box(i_quadrant, 2),&
              data_box(i_quadrant, 3), data_box(i_quadrant, 4),&
              data_box(i_quadrant, 5), data_box(i_quadrant, 6)
          
              d_box(i_quadrant,1) = data_box(i_quadrant, 4) - data_box(i_quadrant, 1)+1
              d_box(i_quadrant,2) = data_box(i_quadrant, 5) - data_box(i_quadrant, 2)+1
              d_box(i_quadrant,3) = data_box(i_quadrant, 6) - data_box(i_quadrant, 3)+1  
              n_cell_box(i_quadrant) = d_box(i_quadrant,1) * d_box(i_quadrant,2) * d_box(i_quadrant,3)
   END DO
  close(unit=1)
  data_box = data_box + 1
  !==================================================  
  
  CALL h5open_f(error)                                                             ! Initialise FORTRAN interface.

  Do i_timepoints = 1, n_timepoints
      
      write(i_timepoints_word, '(i5.5)' ) i_timepoints-1
      
      filename_hydro = "chombo"//i_timepoints_word//".hdf" ! File name
      
      call h5fopen_f (hydro_resultsfolder//filename_hydro, H5F_ACC_RDWR_F, file_id, error)   ! Open an existing file.
      call h5dopen_f(file_id, dsetname, dset_id, error)                         ! Open an existing dataset.

      !=================
      ! Read the dataset.
      call h5dread_f(dset_id, H5T_NATIVE_REAL, data_hydro_code(:,i_timepoints), data_dims, error)
      !=================

      !==========================================================================================
      !Convert data_hydro_code in new structure, data_hydro_code: 0 density, 1 px, 2 py, 3 E, 4 0?, 5 ?, 6 ?
      Do i_zgrid = 1, n_zgrid
           Do i_ygrid = 1, n_ygrid
               Do i_xgrid = 1, n_xgrid
               
                  !Find the correct quadrant
                  Do i_quadrant = 1, n_quadrant
                      if ((data_box(i_quadrant, 1) <= i_xgrid) .and. (i_xgrid <= data_box(i_quadrant, 4)) .and. &
                      (data_box(i_quadrant, 2) <= i_ygrid) .and. (i_ygrid <= data_box(i_quadrant, 5)) .and. &
                      (data_box(i_quadrant, 3) <= i_zgrid) .and. (i_zgrid <= data_box(i_quadrant, 6)) ) then

                          n_cell_total = sum(n_cell_box(1:i_quadrant-1)) * n_para_hdf5
                          dx = d_box(i_quadrant,1)
                          dy = d_box(i_quadrant,2)
                          x_in_box = i_xgrid - data_box(i_quadrant, 1) + 1
                          y_in_box = i_ygrid - data_box(i_quadrant, 2) + 1
                          z_in_box = i_zgrid - data_box(i_quadrant, 3) + 1
                          
                          i_run = n_cell_total + (z_in_box - 1) * dx * dy + (y_in_box - 1) * dx  + x_in_box 
                          !print*, i_run,"In Echt:", i_xgrid, i_ygrid, i_zgrid,"In BOx:", x_in_box, y_in_box, z_in_box
                          
                          
                          !============================== Magnetic field: off (0) or on (1)
                          IF (l_MHD == 0) THEN
                          
                            ! 4., 5., 6.: gas vx, vy, vz: 
                            data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4) =  &
                                data_hydro_code(i_run + n_cell_box(i_quadrant), i_timepoints) * unit_velocity &
                                /(data_hydro_code(i_run,i_timepoints))                                 ! in m/s
                            data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5) = &
                                data_hydro_code(i_run + 2*n_cell_box(i_quadrant), i_timepoints) * unit_velocity &
                                /(data_hydro_code(i_run,i_timepoints))                                 ! in m/s 
                            if (n_zgrid == 1) then 
                            ! 2D Simulation
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6) = 0
                            else
                            ! 3D Simulation
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6) = &
                                data_hydro_code(i_run+3*n_cell_box(i_quadrant),i_timepoints) * unit_velocity &
                                /(data_hydro_code(i_run,i_timepoints))                                ! in m/s
                            end if    
                          
                            ! 7.: gas density
                            data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,7) = data_hydro_code(i_run,i_timepoints)*unit_density ! in kg*m**(-3) 
                          
                            ! 8.: gas temperature
                            data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,8) = &
                                (data_hydro_code(i_run + (n_para_hdf5 - 4) * n_cell_box(i_quadrant),i_timepoints) * &  
                                unit_energy * (unit_length)**(-3.0) - &
                                0.5 * data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,7) * &
                                (data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4)**2.0 +   &
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5)**2.0 +   &
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6)**2.0)) * &
                                (gamma - 1.0)*(m_amu * mu_chem)/(data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,7) * kB) ! in K 
                                !Why (n_para_hdf5 - 4)? Because for the 2D case, it should be the 3rd block, and for the 3D case the 4th block. And n_para_hdf5 is increasing by 1 from the 2D to 3D case.
                                
                           ELSE IF  (l_MHD == 1) THEN    ! magnetic field on
                                ! 4., 5., 6.: gas vx, vy, vz: 
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4) =  &
                                    data_hydro_code(i_run + n_cell_box(i_quadrant), i_timepoints) * unit_velocity &
                                    /(data_hydro_code(i_run,i_timepoints))                                 ! in m/s
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5) = &
                                    data_hydro_code(i_run + 2*n_cell_box(i_quadrant), i_timepoints) * unit_velocity &
                                    /(data_hydro_code(i_run,i_timepoints))                               
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6) = &
                                    data_hydro_code(i_run + 3*n_cell_box(i_quadrant), i_timepoints) * unit_velocity &
                                    /(data_hydro_code(i_run,i_timepoints)) 
                          
                                ! 7.: gas density
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,7) = data_hydro_code(i_run,i_timepoints)&
                                    *unit_density ! in kg*m**(-3)                                     
 
                                ! 9., 10., 11.: gas Bx, By, Bz: 
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,9) =  &
                                    data_hydro_code(i_run + 5*n_cell_box(i_quadrant), i_timepoints) * unit_mag ! in Tesla
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,10) = &
                                    data_hydro_code(i_run + 6*n_cell_box(i_quadrant), i_timepoints) * unit_mag ! in Tesla
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,11) = &
                                    data_hydro_code(i_run + 7*n_cell_box(i_quadrant), i_timepoints) * unit_mag ! in Tesla
 
                                ! 8.: gas temperature
                                data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,8) = &
                                    (data_hydro_code(i_run + 4 * n_cell_box(i_quadrant),i_timepoints) * &  
                                    unit_energy * (unit_length)**(-3.0) - &
                                    0.5 * data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,7) * &
                                    (data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4)**2.0 +   &
                                    data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5)**2.0 +   &
                                    data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6)**2.0) -  &
                                    0.5 * (data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,9)**2.0 +   &
                                    data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,10)**2.0 +   &
                                    data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,11)**2.0)) * &
                                    (gamma - 1.0)*(m_amu * mu_chem)/(data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,7) * kB) ! in K                               
                              
                              
                          END IF
                          !============================== Magnetic field: off (0) or on (1)                          
                              
                          EXIT
                      END IF
                  END DO
                  
               END DO 
           END DO        
      END DO        
      !==========================================================================================
      
      ! Close the dataset, file, Fortran interface
      call h5dclose_f(dset_id, error)
      call h5fclose_f(file_id, error)
      
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

  END DO ! i_timepoints

  CALL h5close_f(error)
  
  DEALLOCATE(data_hydro_code, data_box, d_box, n_cell_box)
      
end subroutine hydro_data_astrobear  ! #s1.2

!################################################################################################################

end module init_astrobear
