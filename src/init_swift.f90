module init_swift
  use omp_lib
  use datatype
  use variable
  use tools
  
  implicit none
    !--------------------------------------------------------------------------!
!     PRIVATE
!         CHARACTER(len=4), PARAMETER :: file_a  = "huhu"
    !--------------------------------------------------------------------------!
    PUBLIC :: meta_data_swift, hydro_data_swift
    !--------------------------------------------------------------------------!
CONTAINS 

!################################################################################################################  
  
subroutine meta_data_swift ! #s1.1
    use omp_lib
    use datatype
    use variable
    use tools
  
    implicit none

    REAL (kind=r2):: dummy, t_snap   
    
    !======================
    !0./4. Read mu_chem, gamma
    !======================
 
    call parse_real("Xmu               ", mu_chem, 1, Swi_file_folder//'meta.swift') 
    call parse_real("gamma             ", gamma,   1, Swi_file_folder//'meta.swift') 
     
    !======================
    !1. Read units
    !======================

    call parse_real("RSCALE       ", unit_density,     1, Swi_file_folder//'meta.swift') ! in cm**(-3)   
    call parse_real("LSCALE       ", unit_length,      1, Swi_file_folder//'meta.swift') ! in cm                                       
    call parse_real("TIMESCALE    ", unit_time,        1, Swi_file_folder//'meta.swift') ! in seconds
    call parse_real("TEMPSCALE    ", unit_temperature, 1, Swi_file_folder//'meta.swift') ! in Kelvin  
 
    !Convert from cgs into SI-units:
    unit_density  = mu_chem * m_amu * unit_density* 1.0E+6                       ! in kg*m**(-3) 
    unit_length   = unit_length  *   0.01                                        ! in m
 
    unit_velocity = unit_length/unit_time * 1.0E+5                               ! in m/s
    unit_energy   = unit_density  *unit_length**(5.0) * unit_time**(-2.0)        ! kg m**2/s**2    
  
    !======================
    !2. Read number of timepoints and calculate time points.
    !======================
  
     open(1, file = Swi_file_folder//'Ascii_data/totals.dat')
         n_timepoints = 0
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
  
    open(unit=1, file = Swi_file_folder//'Ascii_data/totals.dat', action="read", status="unknown", form="formatted")    
        DO i_timepoints = 1, n_timepoints
            read(unit=1,fmt=*) time(i_timepoints)
            time(i_timepoints) = time(i_timepoints) * unit_time
        END DO
    close(unit=1)  
 
 !======================
 !3. Read grid extension
 !======================  
 ! Read n_xgrid, n_ygrid, n_zgrid, length_x, length_y, le   CHARACTER(LEN=15) :: filename_hydro                   ! Chombo name (Chombo00xyz.hdf)ngth_z
 
  call parse_real("GmX               ", dummy, 1, Swi_file_folder//'meta.swift') 
  n_xgrid = nint(dummy)
  call parse_real("GmX               ", dummy, 2, Swi_file_folder//'meta.swift') 
  n_ygrid = nint(dummy)
  call parse_real("GmX               ", dummy, 3, Swi_file_folder//'meta.swift') 
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
 
  call parse_real("GxBounds          ", length_x, 4, Swi_file_folder//'meta.swift')    ! Length domain in x-direction, in code-units
  call parse_real("GxBounds          ", length_y, 5, Swi_file_folder//'meta.swift')    ! Length domain in y-direction, in code-units 
  call parse_real("GxBounds          ", length_z, 6, Swi_file_folder//'meta.swift')    ! Length domain in z-direction, in code-units
  
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
 
    call parse_real("origin", origin_x, 1, Swi_file_folder//'meta.swift')    ! in pixel
    call parse_real("origin", origin_y, 2, Swi_file_folder//'meta.swift')    ! in pixel
    call parse_real("origin", origin_z, 3, Swi_file_folder//'meta.swift')    ! in pixel
    call parse_real("radius_r",   region_r, 1, Swi_file_folder//'meta.swift')! in pc
  
    ! Convert from pc into SI-units:
    region_r = region_r / m_in_pc ! in m

end subroutine meta_data_swift  ! #s1.1

!################################################################################################################

subroutine hydro_data_swift  ! #s1.2
  use omp_lib
  use datatype
  use variable
  
  implicit none

  INTEGER :: i_xgrid, i_ygrid, i_zgrid

  CHARACTER(LEN=28), PARAMETER  :: hydro_resultsfolder = Swi_file_folder//"Ascii_data/"             ! folder with stored hydro-files
  CHARACTER(LEN=12) :: filename_hydro                                                         ! file name of the hydro-files  
  
  ALLOCATE(data_gas(1:n_xgrid, 1:n_ygrid, 1:n_zgrid, 1:n_timepoints, 4:n_paragas))            ! Global importance, contains all hydro_data
  ! data_gas(:,:,:,:,{1,2,3})  not used!
  ! data_gas(:,:,:,:,4) = gas vx
  ! data_gas(:,:,:,:,5) = gas vy
  ! data_gas(:,:,:,:,6) = gas vz
  ! data_gas(:,:,:,:,7) = gas density
  ! data_gas(:,:,:,:,8) = gas temperature

  Do i_timepoints = 1, n_timepoints
    
    if (i_timepoints < 11) THEN
        write(i_timepoints_word, '(i1.1)' ) i_timepoints-1
    elseif (i_timepoints < 101) THEN    
        write(i_timepoints_word, '(i2.2)' ) i_timepoints-1
    elseif (i_timepoints < 1001) THEN    
        write(i_timepoints_word, '(i3.3)' ) i_timepoints-1
    elseif (i_timepoints < 10001) THEN    
        write(i_timepoints_word, '(i4.4)' ) i_timepoints-1
    elseif (i_timepoints < 100001) THEN    
        write(i_timepoints_word, '(i5.5)' ) i_timepoints-1        
    else 
        print*, 'Paperboats ERROR 21 in reading in hydro data for hydro code Swift. Too many timesteps (>=100000).'
        STOP  
    end if    
        
    filename_hydro = "sedov_slice_" ! File name

    !=================    
    ! 4.: gas vx
    open(unit=1, file = hydro_resultsfolder//filename_hydro//TRIM(i_timepoints_word)//'_2.dat',&
        action="read", status="unknown", form="formatted")    
        DO i_zgrid = 1, n_zgrid 
            DO i_xgrid = 1, n_xgrid      
                read(unit=1,fmt=*) (data_gas(i_xgrid,i_ygrid,  i_zgrid, i_timepoints,4),&
                i_ygrid = 1, n_ygrid)
            END DO
        END DO
        data_gas(:,:,:, i_timepoints,4) = data_gas(:,:,:, i_timepoints,4) * unit_velocity ! in m/s 
    close(unit=1)
 
    !=================
    ! 5.: gas vy      
    open(unit=1, file = hydro_resultsfolder//filename_hydro//TRIM(i_timepoints_word)//'_3.dat',&
        action="read", status="unknown", form="formatted")    
        DO i_zgrid = 1, n_zgrid 
            DO i_xgrid = 1, n_xgrid      
                read(unit=1,fmt=*) (data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5),&
                i_ygrid = 1, n_ygrid) 
            END DO
        END DO
        data_gas(:,:,:, i_timepoints,5) = data_gas(:,:,:, i_timepoints,5) * unit_velocity ! in m/s         
    close(unit=1)
    !=================
    ! 6.: gas vz     
    open(unit=1, file = hydro_resultsfolder//filename_hydro//TRIM(i_timepoints_word)//'_4.dat',&
        action="read", status="unknown", form="formatted")    
        DO i_zgrid = 1, n_zgrid 
            DO i_xgrid = 1, n_xgrid      
                read(unit=1,fmt=*) (data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6),&
                i_ygrid = 1, n_ygrid)
            END DO
        END DO
!         data_gas(:,:,:, i_timepoints,6) = data_gas(:,:,:, i_timepoints,6) * unit_velocity ! in m/s 
        data_gas(:,:,:, i_timepoints,6) = 0.0            
    close(unit=1)
    !=================
    ! 7.: gas density      
    open(unit=1, file = hydro_resultsfolder//filename_hydro//TRIM(i_timepoints_word)//'_0.dat',&
        action="read", status="unknown", form="formatted")    
        DO i_zgrid = 1, n_zgrid 
            DO i_xgrid = 1, n_xgrid      
                read(unit=1,fmt=*) (data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,7),&
                i_ygrid = 1, n_ygrid)  
            END DO
        END DO            
        data_gas(:,:,:, i_timepoints,7) = data_gas(:,:,:, i_timepoints,7) * unit_density ! in kg/m**3          
    close(unit=1)
    !=================    
    ! 8.: gas temperature  
    open(unit=1, file = hydro_resultsfolder//filename_hydro//TRIM(i_timepoints_word)//'_1.dat',&
        action="read", status="unknown", form="formatted")    
        DO i_zgrid = 1, n_zgrid 
            DO i_xgrid = 1, n_xgrid      
                read(unit=1,fmt=*) (data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,8),&
                i_ygrid = 1, n_ygrid)  
            END DO
        END DO            
        data_gas(:,:,:, i_timepoints,8) = data_gas(:,:,:, i_timepoints,8) * unit_temperature ! in K       
    close(unit=1)    
    !=================   
    
!     DO i_zgrid = 1, n_zgrid 
!         DO i_ygrid = 1, n_ygrid    
!             DO i_xgrid = 1, n_xgrid          
! !                 data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4) = &
! !                    10000.0* real(i_xgrid)/(real(i_xgrid**2.0 + i_ygrid**2.0)**0.5)* unit_velocity 
! !                 data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5) = &
! !                    10000.0* real(i_ygrid)/(real(i_xgrid**2.0 + i_ygrid**2.0)**0.5)* unit_velocity
!                 data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4) = data_gas(i_xgrid, 1, i_zgrid, i_timepoints,4)                   
!                 data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5) = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5)   
! !                 data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6) = 0.0*unit_velocity
!                  data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,7) =   data_gas(i_xgrid, 1, i_zgrid, i_timepoints,7) 
! !                 data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,8) = 2183.0 * unit_temperature
!             END DO
!         END DO
!     END DO
            
    
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
        print*, 'Paperboats ERROR 25 in reading in hydro data for hydro code Swift. Too many timesteps (>=100000).'
        STOP
    end if
    
  END DO

end subroutine hydro_data_swift  ! #s1.2

!################################################################################################################

end module init_swift
