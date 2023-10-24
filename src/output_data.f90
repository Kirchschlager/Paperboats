module output_data
    use omp_lib
    use datatype
    use variable
  
    implicit none
    !--------------------------------------------------------------------------!
    PUBLIC :: output_paperboats, output_gasdata, output_dustdata,              &
    output_dustdistributiondata 
    !--------------------------------------------------------------------------!
CONTAINS 

!##################################################################################################################################
!##################################################################################################################################

subroutine output_paperboats        ! #s3
  use omp_lib
  use datatype
  use variable

  implicit none  

   call output_gasdata
   call output_dustdata 
   call output_dustdistributiondata

end subroutine output_paperboats    ! #s3

!##################################################################################################################################
!##################################################################################################################################

subroutine output_gasdata           ! #s3.1
    use omp_lib
    use datatype
    use variable
  
    implicit none  
  
    INTEGER :: i_xgrid, i_ygrid, i_xgrid_h, i_ygrid_h, n_pol_x, n_pol_y
    REAL (kind=r2) ::  B_max, B_here,  B_scale
 
    write(i_timepoints_word, '(i5.5)' ) i_timepoints
                        
    !======================================================== 
    ! 1. Store the gas mass distribution

    1 format (*(G0,1X))    
    open(unit=1, file="./"//Dens_map_folder//"data_gas/Density_"//i_timepoints_word//".dat", &
        &action="write", status="unknown", form="formatted")  
        write(unit=1,fmt=1) "# Matrix 1. number density [cm^{-3}] /// Matrix 2. T [K]" 
        Do i_xgrid = 1, n_xgrid
                  write(unit=1,fmt=1) (data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 7)/(mu_chem * m_amu) * 1.0E-6&
                   &, i_ygrid = 1, n_ygrid) ! in cm**(-3)   
 
!                 write(unit=1,fmt=1) (data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 4),&
!                 i_ygrid = 1, n_ygrid)                                                                    ! vx, in m/s 
                
!                  write(unit=1,fmt=1) (data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 5),&
!                  i_ygrid = 1, n_ygrid)                                                                   ! vy, in m/s      

!                  write(unit=1,fmt=1) (((data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 4))**2.0+&
!                   &(data_gas(i_xgrid,i_ygrid,(n_zgrid+1)/2,i_timepoints,5))**2.0)**0.5,i_ygrid=1,n_ygrid) ! v, in m/s    
        END DO
        write(unit=1,fmt=1) ""
        write(unit=1,fmt=1) ""        
        Do i_xgrid = 1, n_xgrid
                write(unit=1,fmt=1) (data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 8), i_ygrid = 1, n_ygrid)                            ! in K                                    
        END DO
        IF  (l_MHD == 1) THEN 
            write(unit=1,fmt=1) ""
            write(unit=1,fmt=1) ""        
            Do i_xgrid = 1, n_xgrid
                write(unit=1,fmt=1) ((data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 9)**2.0 + &
                                     &data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 10)**2.0  + &
                                     &data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 11)**2.0)**0.5 * 10.0**10.0,&
                                     &i_ygrid = 1, n_ygrid)   ! in Mikro-Gauss                                    
            END DO
        END IF    
    close(unit=1) 
    
    !======================================================== 
    ! 2. Find the maximum of the magnetic field
    
    IF  (l_MHD == 1) THEN  
        Do i_xgrid = 1, n_xgrid
            Do i_ygrid = 1, n_ygrid 
                B_here =  (data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 9)**2.0 + &
                    data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 10)**2.0 + &
                    data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 11)**2.0 )**0.5 * 10.0**10.0
                IF ((i_xgrid == 1) .and. (i_ygrid == 1)) THEN
                    B_max  = B_here       
                ELSE     
                    If (B_here > B_max) THEN
                        B_max = B_here
                    END IF                      
                END IF  
            END DO
        END DO 
        IF (i_timepoints > 1)  THEN
            If (B_max > B_max_tot) THEN
                B_max_tot = B_max
            END IF         
        ELSE        
            B_max_tot   = B_max        
        END IF
    END IF    

    !======================================================== 
    ! 3. Magnetic field vector field 
    
    IF  (l_MHD == 1) THEN
     open(unit=1, file="./"//Dens_map_folder//"data_gas/Mag_field_"//i_timepoints_word//".dat", &
        &action="write", status="unknown", form="formatted")  
        write(unit=1,fmt=1) "#   x                  y  " 
        
        n_pol_x = nint(real(n_xgrid)/20.0)
        n_pol_y = nint(real(n_ygrid)/20.0)
        Do i_xgrid_h = 1, n_pol_x      
            Do i_ygrid_h = 1, n_pol_y         
                i_xgrid = nint(real(n_xgrid)/real(n_pol_x)) * (i_xgrid_h -1) + nint(0.5 * real(n_xgrid)/real(n_pol_x))
                i_ygrid = nint(real(n_ygrid)/real(n_pol_y)) * (i_ygrid_h -1) + nint(0.5 * real(n_xgrid)/real(n_pol_x))
                
                If ((data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 9)**2.0 + &
                        &data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 10)**2.0) <1.0E-40) THEN
                    B_scale = 0.0
                ELSE
                    B_scale = nint(real(n_xgrid)/real(n_pol_x)) * 0.4 * &
                        & (data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 9)**2.0 + &
                        &data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 10)**2.0)**(-0.5)
                END IF
 
                write(unit=1,fmt=1) i_xgrid + B_scale * data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 9),&
                    & i_ygrid + B_scale * data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 10) 
                write(unit=1,fmt=1) i_xgrid - B_scale * data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 9),&
                    & i_ygrid - B_scale * data_gas(i_xgrid, i_ygrid, (n_zgrid+1)/2, i_timepoints, 10)            
                write(unit=1,fmt=1) ""
            END DO
        END DO            
    END IF    

end subroutine output_gasdata      ! #s3.1

!##################################################################################################################################
!##################################################################################################################################

subroutine output_dustdata         ! #s3.2
    use omp_lib
    use datatype
    use variable
  
    implicit none
    
    INTEGER :: i_xgrid, i_ygrid, i_dusttype, i_asize  
    INTEGER :: i_asize_reduced, n_asize_reduced, i_dum
  
    REAL (kind=r2), DIMENSION(:,:), ALLOCATABLE :: densi_min, densi_max
  
    ALLOCATE(a_word(1:n_asize), densi_min(1:n_dusttype, 1: n_asize), densi_max(1:n_dusttype, 1: n_asize))

    Do i_dusttype = 1, n_dusttype
  
        write(i_timepoints_word, '(i5.5)' ) i_timepoints
        write(i_dusttype_word,   '(i1.1)' ) i_dusttype
        
        !=================================================================================================
        ! The data sets can be very large (GBs). To reduce the amount of data, store only max. 4 grain sizes.
        n_asize_reduced = min(n_out_sizes, i_max_ini(i_dusttype))
        
        Do i_asize_reduced = 1, n_asize_reduced
            if (n_asize_reduced == 1) THEN
                i_asize = i_out_size
            else 
                if (i_max_ini(i_dusttype) < n_out_sizes) Then
                    i_asize = i_asize_reduced
                else 
                    i_asize = nint((real((i_max_ini(i_dusttype) - 1) * i_asize_reduced)&
                    + real(n_asize_reduced - i_max_ini(i_dusttype))) / real(n_asize_reduced - 1))
                end if
            end if    
            
            write(a_word(i_asize), '(f6.1)' ) a(i_asize, i_dusttype) * 1.0E+9
            
            if (a(i_asize, i_dusttype) > 0.999E-5) THEN
                write(a_word(i_asize), '(f6.0)' ) a(i_asize, i_dusttype) * 1.0E+9
            end if            
        END DO
            
        !======================================================== 
        ! 1. Store the dust mass distribution

        1 format (*(G0,1X))        
        open(unit=1, file="./"//Dens_map_folder//"data_dust/Density_"//i_dusttype_word//"_"//i_timepoints_word//".dat", &
                & action="write", status="unknown", form="formatted")  
            write(unit=1,fmt=1) "# Number density [cm^{-3}]. Number of rows = number of cells in x-direction. &
                & Each single grain size is in a block, from small to big radius, starting from the top of the file." 

            if (n_asize_reduced == 1) THEN  
                Do i_xgrid = 1, n_xgrid
                      write(unit=1,fmt=1) (data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_out_size, 4)* 1.0E-6,&
                                               i_ygrid = 1, n_ygrid)                                                                           ! number density in cm**(-3)
                                               
!                    write(unit=1,fmt=1) ((data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_out_size, 1)**2.0+&   
!                                         data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_out_size, 2)**2.0)**0.5,&
!                                                i_ygrid = 1, n_ygrid)                                                                         ! velocity (vdustx**2 + vdusty**2)**0.5, in m/s  

!                    write(unit=1,fmt=1) ((data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_out_size, 1)),&
!                                                i_ygrid = 1, n_ygrid)                                                                           ! velocity (vdustx), in m/s  

!                    write(unit=1,fmt=1) ((data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_out_size, 2)),&
!                                                i_ygrid = 1, n_ygrid)                                                                           ! velocity (vdusty), in m/s                                                 
                END DO
            else                        
                if (i_max_ini(i_dusttype) < n_out_sizes) THEN
                    Do i_asize = 1, i_max_ini(i_dusttype)
                        Do i_xgrid = 1, n_xgrid
                            write(unit=1,fmt=1) (data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) &
                                            * 1.0E-6, i_ygrid = 1, n_ygrid)                                                                 ! number density in cm**(-3)   
                        END DO
                        write(unit=1,fmt=1) ""
                        write(unit=1,fmt=1) ""
                    END DO    
                ELSE !  -> i_max_ini(i_dusttype) >= n_out_sizes  
                    Do  i_asize_reduced = 1, n_asize_reduced
                        i_dum = nint((real((i_max_ini(i_dusttype) - 1) * i_asize_reduced) + real(n_asize_reduced &
                                                - i_max_ini(i_dusttype)) )/ real(n_asize_reduced - 1)) 
                        Do i_xgrid = 1, n_xgrid                    
                            write(unit=1,fmt=1)  (data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_dum, 4) * 1.0E-6,&       ! number density in cm**(-3) 
                                                i_ygrid = 1, n_ygrid)                 
                        END DO
                        write(unit=1,fmt=1) ""
                        write(unit=1,fmt=1) ""
                    END DO
                end if                            
            END IF                          
        close(unit=1)
 
!         1 format (*(G0,1X)) 
!         open(unit=1, file="./"//Dens_map_folder//"data_dust/Density_"//i_dusttype_word//"_"//i_timepoints_word//".dat", &
!                 & action="write", status="unknown", form="formatted")  
!             write(unit=1,fmt=1) "# x [pc]                    y [pc]                    &
!                 &number density [cm^{-3}] for each single grain size" 
!             Do i_xgrid = 1, n_xgrid
!                 Do i_ygrid = 1, n_ygrid
!                     if (n_asize_reduced == 1) THEN              
!                         write(unit=1,fmt=1) data_grid(i_xgrid, i_ygrid, (n_zgrid+1)/2, 1) * m_in_pc,&                               ! x in pc
!                                             data_grid(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2) * m_in_pc, &                              ! y in pc
!                                             data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_out_size, 4)* 1.0E-6    ! number density in cm**(-3)
!                                            !(data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_out_size, 1)**2.0+&    ! velocity (vdustx**2 + vdusty**2)**0.5, in m/s
!                                            ! data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_out_size, 2)**2.0)**0.5                                             
!                     else                        
!                                                                
!                     if (i_max_ini(i_dusttype) < n_out_sizes) THEN
!                         write(unit=1,fmt=1) data_grid(i_xgrid, i_ygrid, (n_zgrid+1)/2, 1) * m_in_pc,&                               ! x in pc
!                                             data_grid(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2) * m_in_pc, &                              ! y in pc
!                                            (data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4)* 1.0E-6, &    ! number density in cm**(-3)
!                                                 i_asize = 1, i_max_ini(i_dusttype)) 
!                     ELSE !  -> i_max_ini(i_dusttype) >= n_out_sizes                    
!                         write(unit=1,fmt=1) data_grid(i_xgrid, i_ygrid, (n_zgrid+1)/2, 1) * m_in_pc,&                               ! x in pc
!                                             data_grid(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2) * m_in_pc, &                              ! y in pc
!                                             (data_dust_main(i_xgrid, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype, &                       ! number density in cm**(-3)
!                                                 nint((real((i_max_ini(i_dusttype) - 1) * i_asize_reduced) + real(n_asize_reduced &
!                                                 - i_max_ini(i_dusttype)) )/ real(n_asize_reduced - 1)), 4) * 1.0E-6,&
!                                                 i_asize_reduced = 1, n_asize_reduced)      ! in cm**(-3)
!                     end if                            
!                     END IF
!                 END Do
!                 write(unit=1,fmt=1)  ""                             
!             END DO
!         close(unit=1)        
    End do 
    !======================================================== 
    ! 2. Find the maximum of the dust mass distribution
  
    Do i_dusttype = 1, n_dusttype
        DO i_asize = 1, n_asize
            densi_max(i_dusttype, i_asize) =  max(1.0E-13,maxval(data_dust_main(:,:,(n_zgrid+1)/2,2,i_dusttype,i_asize,4)))*1.0E-6 ! 1.0E-10, otherwise, cbrange empty in gnuplot files.
            densi_min(i_dusttype, i_asize) =  max(1.0E-16,minval(data_dust_main(:,:,(n_zgrid+1)/2,2,i_dusttype,i_asize,4)))*1.0E-6 ! 1.0E-10, otherwise, cbrange empty in gnuplot files.
            
            IF (i_timepoints > 1)  THEN
                If (densi_max(i_dusttype, i_asize) > densi_max_tot(i_dusttype, i_asize)) THEN
                    densi_max_tot(i_dusttype, i_asize) = densi_max(i_dusttype, i_asize)
                END IF      
                If (densi_min(i_dusttype, i_asize) < densi_min_tot(i_dusttype, i_asize)) THEN
                    densi_min_tot(i_dusttype, i_asize) = densi_min(i_dusttype, i_asize)
                END IF                 
            ELSE
                densi_max_tot(i_dusttype, i_asize) = densi_max(i_dusttype, i_asize)
                densi_min_tot(i_dusttype, i_asize) = densi_min(i_dusttype, i_asize)                
            END IF  
        END DO
    END DO  
  
    DEALLOCATE(a_word, densi_min, densi_max)
  
end subroutine output_dustdata          ! #s3.2

!##################################################################################################################################
!##################################################################################################################################

subroutine output_dustdistributiondata ! #s3.3
    use omp_lib
    use datatype
    use variable

    implicit none  
    INTEGER :: i_dusttype, i_asize  
  
    Real (kind=r2), DIMENSION(:,:), ALLOCATABLE :: distrib_minmax
    Real (kind=r2), DIMENSION(:), ALLOCATABLE :: alpha_sum
    
    REAL (kind=r2) :: help
  
    INTEGER :: i, i_xgrid, i_ygrid, i_zgrid
  
    logical :: exist  

    ALLOCATE(alpha_sum(1:n_dusttype), distrib_minmax(1:n_asize, 1:n_dusttype))

    write(i_timepoints_word, '(i5.5)' ) i_timepoints  
  
    !======================================================== 
    ! 1. Store the dust grain size distribution data

    1 format (*(G0,1X))
    open(unit=1, file="./"//Size_dist_folder//"data/Particlenumbers_"//i_timepoints_word//".dat", action="write",&
        status="unknown", form="formatted")    
        write(unit=1,fmt=1)  ("#                 ", Material_name(i_dusttype), i_dusttype = 1, n_dusttype), "          "
        write(unit=1,fmt=1)  ("# a [m]                      N [m**(-3)]               ", i_dusttype = 1, n_dusttype) 
        Do i_asize = 1, n_asize
            Do i_dusttype = 1, n_dusttype
                alpha_sum(i_dusttype) = sum(data_dust_main(:, :, :, 2, i_dusttype, i_asize, 4))/&
                    (n_xgrid*n_ygrid*n_zgrid* 1.0E+9*a(1, i_dusttype)*Delta_rad**(real(i_asize)-1.0))
            END DO
            write(unit=1,fmt=1)  (a(i_asize, i_dusttype), alpha_sum(i_dusttype), i_dusttype = 1, n_dusttype)
        END DO
    close(unit=1)
  
    DEALLOCATE(alpha_sum)
  
    !======================================================== 
    ! 2. Store the dust mass as a function of time
  
    inquire(file="./"//Size_dist_folder//"data/Dustmass_evolution_total.dat", exist=exist)
    if (exist) then
    open(unit=1, file="./"//Size_dist_folder//"data/Dustmass_evolution_total.dat", status="old", position="append", action="write")
    else
    open(unit=1, file="./"//Size_dist_folder//"data/Dustmass_evolution_total.dat", status="new", action="write")
        write(unit=1,fmt=1) "# t ", ("                     M/M(t=0)                  M(without large grains)/M(t=0)   &
            &(Dust + Dusty gas)  ", i_dusttype = 1, n_dusttype), "Number of collisions (1, 2, 3, 4)                       &
            &M_regular_gas/M(t=0)"    
    end if
        write(unit=1,fmt=1) time_unit_scaler * time(i_timepoints), &                                                              ! t
                            ((sum(data_dust_main(:,:,:,2,i_dusttype,1:n_asize, 5)) + M_larger_grains(i_dusttype))&  ! (Mass(t)+Mass_of_larger_grains)/ Mass(t=0)
                                / sum(data_dust_main(:, :, :, 1, i_dusttype, 1:n_asize, 5)), &
                            (sum(data_dust_main(:, :, :, 2, i_dusttype, 1:n_asize, 5)) &                                            ! (Mass(t))/ Mass(t=0)
                                / sum(data_dust_main(:, :, :, 1, i_dusttype, 1:n_asize, 5))), &
                            (sum(data_dust_main(:, :, :, 2,i_dusttype,1:n_asize, 5)) + M_larger_grains(i_dusttype)&  ! (Mass(t)+Mass_of_larger_grains+ Dusty gas)/ Mass(t=0) = 1.0000
                                + sum(data_dust_main(:, :, :, 2, i_dusttype, 0, 4)) *m_grainatom(i_dusttype)) &
                                / sum(data_dust_main(:, :, :, 1, i_dusttype, 1:n_asize, 5)), i_dusttype = 1, n_dusttype),&          ! for i_dusttype = 1, n_dusttype    
                                (count_colli(i_timepoints, i), i = 1, 4), &
                                M_regular_gas/ sum(data_dust_main(:, :, :, 1, i_dusttype, 1:n_asize, 5))        
    close(unit=1)
    
    ! In scenario 2, consider only the region within region_r around the origin origin_x, origin_y. Assuming that nothing is happening outside.
    IF (((l_scenario == 2) .or. (l_scenario==3)) .and. ((l_hydrocode == 2) .or. (l_hydrocode == 3) &
        .or. (l_hydrocode == 4).or. (l_hydrocode == 5))) THEN      ! FK, 15/04/2021, 03/01/2022
        inquire(file="./"//Size_dist_folder//"data/Dustmass_evolution_total_reg.dat", exist=exist)
        if (exist) then
        open(unit=1, file="./"//Size_dist_folder//"data/Dustmass_evolution_total_reg.dat", status="old", position="append",&
            &action="write")
        else
        open(unit=1, file="./"//Size_dist_folder//"data/Dustmass_evolution_total_reg.dat", status="new", action="write")
            write(unit=1,fmt=1) "# t                     M/M(t=0) "  
        end if
    

        Mass_sum_reg = 0.0_r2
        IF (i_timepoints == 1) THEN
            Mass_sum_reg_0 = 0.0_r2
        END IF          
         
         Do i_xgrid = 1, n_xgrid
              Do i_ygrid = 1, n_ygrid     
                    Do i_zgrid = 1, n_zgrid
                        IF ( region_r**2.0 > &
                        ((data_grid(i_xgrid, i_ygrid, i_zgrid, 1) - data_grid(nint(origin_x),nint(origin_y),&
                        nint(origin_z), 1))**2.0+&
                         (data_grid(i_xgrid, i_ygrid, i_zgrid, 2) - data_grid(nint(origin_x),nint(origin_y),&
                         nint(origin_z), 2))**2.0+&
                         (data_grid(i_xgrid, i_ygrid, i_zgrid, 3) - data_grid(nint(origin_x),nint(origin_y),&
                         nint(origin_z), 3))**2.0))THEN
                            IF (i_timepoints == 1) THEN
                                Mass_sum_reg_0 = &
                                    & Mass_sum_reg_0 + sum(data_dust_main(i_xgrid, i_ygrid, i_zgrid, 1, 1:n_dusttype, 1:n_asize, 5))
                            END IF               
                            Mass_sum_reg = &
                                & Mass_sum_reg + sum(data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, 1:n_dusttype, 1:n_asize, 5))                               
                        END IF
                    END DO
              END DO
        END DO
        
        write(unit=1,fmt=1) time_unit_scaler * time(i_timepoints), Mass_sum_reg/Mass_sum_reg_0                                
    close(unit=1)    
    END IF
  
    !======================================================== 
    ! 3. Find the maximum and minimum of the grain size distribution  
  
    Do i_dusttype = 1, n_dusttype
        Do i_asize = 1, n_asize
            distrib_minmax(i_asize, i_dusttype) = sum(data_dust_main(:, :, :, 2, i_dusttype, i_asize, 4))&
                / (n_xgrid*n_ygrid*n_zgrid* 1.0E+9*a(1, i_dusttype)*Delta_rad**(real(i_asize)-1.0))     
        END DO
        help = maxval(distrib_minmax(:, i_dusttype))
        
        IF (i_timepoints == 1 ) THEN
            ymin(i_dusttype) = help * 1.0E-6
            ymax(i_dusttype) = help * 10.0           
        ELSE     
            If (help * 10.0 > ymax(i_dusttype)) THEN
                ymax(i_dusttype) = help * 10.0             
            END IF
        END IF            
    END DO
  
    DEALLOCATE(distrib_minmax)
  
end subroutine output_dustdistributiondata ! #s3.3

!##################################################################################################################################

end module output_data

