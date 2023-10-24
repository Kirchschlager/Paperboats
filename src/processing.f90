module processing
    use omp_lib
    use datatype
    use variable
    use tools
    use functs    
    use output_data
  
    implicit none
    !------------------------------------------------------------------------------------------------!
    PUBLIC :: run_paperboats, cell_values, gas_ionisation, charge_calculation, velocity_calculation, &
    gas_drag, coll_drag, plasma_drag, Lorentz_force, coulomb_explosions, grain_grain_collisions,     &
    sputtering, dadt_sputtering, dust_spreading, lift
    !------------------------------------------------------------------------------------------------!
CONTAINS    

!##################################################################################################################################
!##################################################################################################################################

subroutine run_paperboats ! #s2
    use omp_lib
    use datatype
    use variable
    use output_data  
  
    implicit none  
    
    REAL :: t_cpu_r, t_save, t_remain
    REAL (kind=r2) :: rho_gas, Temp_gas

    INTEGER :: i_xgrid, i_ygrid, i_zgrid    

    LOGICAL :: Condition_v, Condition_T, Condition_den1, Condition_den2

    ALLOCATE(data_dust_next(1:n_xgrid, 1:n_ygrid, 1:n_zgrid, 1:n_dusttype, 0:n_asize, 4:5)) !  only for one timepoint!  
  
    
    ALLOCATE(v_dust(1:3,1:n_xgrid,1:n_ygrid,1:n_zgrid, 0:n_timeinter, 1: n_dusttype, 1: n_asize))
    
    IF (l_MHD == 1) THEN
    ALLOCATE(  v_dust_MHD_tra(1:3,1:n_xgrid,1:n_ygrid,1:n_zgrid, 0:n_timeinter, 1: n_dusttype, 1: n_asize), &
             & v_dust_MHD_act(1:3,1:n_xgrid,1:n_ygrid,1:n_zgrid, 0:n_timeinter, 1: n_dusttype, 1: n_asize))
                  ! v_dust         : dust velocity (translation), when magnetic field is off
                  ! v_dust_MHD_tra : dust velocity (translation), when Lorentz force (due to magnetic field) is acting on charged grains
                  ! v_dust_MHD_act : actual dust velocity, when Lorentz force (due to magnetic field) is acting on charged grains
    END IF              

    t_save = secnds(0.0)
 
    Do i_timepoints = 2, n_timepoints
        Delta_t = time(i_timepoints) - time(i_timepoints-1)
        Delta_t_inter = Delta_t * (1.0*n_timeinter)**(-1.0)  
        
         !==================================================
         ! Start Main loop of the program
      
         !$omp parallel num_threads(n_cores) &
         !$omp private(i_xgrid, i_ygrid, i_zgrid, Temp_gas, rho_gas) &
         !$omp private(Condition_v, Condition_T, Condition_den1, Condition_den2) &       
         !$omp shared(data_gas, data_dust_main, data_dust_next, v_dust, v_dust_MHD_tra, v_dust_MHD_act)
         !$omp do schedule(static,1)
         Do i_xgrid = 1, n_xgrid
              Do i_ygrid = 1, n_ygrid
                    Do i_zgrid = 1,n_zgrid
                    !====================================== 
                    ! Define Conditions
                    IF (data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints, 4)**2.0 &
                       & + data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints, 5)**2.0 &
                       & + data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints, 6)**2.0 > 100.0) THEN
                       ! Gasvelocity is larger than 10 m/s
                       Condition_v = .true.
                    else
                       Condition_v = .false.
                    end if

                    IF (l_scenario == 1) THEN
                        IF (data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints, 8) > &
                             & 100.0) THEN
                            ! Gastemperature is larger than 100 K
                            Condition_T = .true.
                        else
                            Condition_T = .false.                        
                        end if   
                    ELSE IF  ((l_scenario == 2) .or. (l_scenario == 3)) THEN
                        IF (data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints, 8) > &
                            ! & 10.0) THEN
                            & data_gas(n_xgrid, n_ygrid, n_zgrid, i_timepoints, 8)) THEN
                            Condition_T = .true.
                        else
                            Condition_T = .false.                        
                        end if                     
                    ELSE
                        print*, "Paperboats ERROR 14, scenario not defined."
                        STOP
                    END IF

                    IF (sum(data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, :,0:n_asize, 4)) > 0.0) THEN
                       ! Dust is in this cell (grains or dusty gas)
                       Condition_den1 = .true.
                    else
                       Condition_den1 = .false.                        
                    end if  
                    IF (sum(data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, :,1:n_asize, 4)) > 0.0) THEN
                       ! Dust is in this cell (grains)
                       Condition_den2 = .true.
                    else
                       Condition_den2 = .false.                           
                    end if                                         
                    !======================================                    

                    !======================================                    
                    IF (Condition_den1) THEN                    
                          !==== 
                          ! 0. Initialise cell values (Temp_gas, rho_gas, Ionisation_grade
                          call cell_values(i_xgrid, i_ygrid, i_zgrid, Temp_gas, rho_gas)

                          !==== 
                          ! 1. Initialise temporal dust parameters, including dust velocity.
                          call velocity_calculation(i_xgrid, i_ygrid, i_zgrid, Temp_gas, rho_gas)
                    end if
                    !======================================                      
                                                        
                    !======================================                   
                    IF ((Condition_den2) .and. (Condition_v .or. Condition_T)) THEN   
                         !==== 
                         ! 2. Calculate coulomb explosions
                         if (n_ce == 1) then
                            call coulomb_explosions(i_xgrid, i_ygrid, i_zgrid) 
                         end if     
                         !==== 
                         ! 3. Calculate grain-grain collisions                     
                         if (n_gg > 0) then
                            call grain_grain_collisions(i_xgrid, i_ygrid, i_zgrid) 
                         end if                             
                         !==== 
                         ! 4. Calculate sputtering
                        if (n_sp > 0) then
                            call sputtering(i_xgrid, i_ygrid, i_zgrid, Temp_gas, rho_gas) 
                         end if                 
                    end if
                    !======================================                      
               
                    !======================================
                    ! 5. Calculate the dust motion in new cells                                                           
                    IF ((Condition_den1)) THEN
                    ! This routine is important, when this cell contains dust of any size, because it might create 
                    ! dust grains of the other sizes, which have to be included here into the data_dust_main-array
                        call dust_spreading(i_xgrid, i_ygrid, i_zgrid) 
                    end if
                    !======================================                      

                    !======================================  
                    ! 6. Calculate the dust grain charge                                                      
                    IF (sum(data_dust_main(i_xgrid, i_ygrid, i_zgrid, 3, :,1:n_asize, 4)) > 0.0) THEN
                    ! Timepoint i_timepoints+1!
                        if (n_plasma > 0) then
                            call charge_calculation(i_xgrid, i_ygrid, i_zgrid, Temp_gas)
                        END if    
                    end if                    
                   !======================================                        
                  
                END DO
            END DO
        END DO
        !$omp end do
        !$omp end parallel
        !==================================================

        call lift
      
        call output_paperboats
      
        !==================================================
        ! Display the time
        t_cpu_r = secnds(t_save)     
        t_remain = t_cpu_r * real(n_timepoints-i_timepoints)
        t_save = secnds(0.0)

        if (t_remain < 60.0) THEN
             write (*,'(A,I3,A, F5.1,A,F6.2,A)') '  2. Processing: ', int((i_timepoints-1)/real(n_timepoints-1)*100.0),&
                 & ' % done. Est. remain. time: ', t_remain, ' s.   Survival: ',&
                 (100.0*(sum(data_dust_main(:, :, :, 2, 1:n_dusttype, 1:n_asize, 5)) &
                 + M_larger_grains(1:n_dusttype))/ sum(data_dust_main(:, :, :, 1, 1:n_dusttype, 1:n_asize, 5))),&
                 " %." //char(27)//'[A'
        ELSEIF (t_remain < 3600.0) THEN
             write (*,'(A,I3,A, F5.1,A,F6.2,A)') '  2. Processing: ', int((i_timepoints-1)/real(n_timepoints-1)*100.0),&
                 & ' % done. Est. remain. time: ', t_remain/60.0, ' min. Survival: ',&
                 (100.0*(sum(data_dust_main(:, :, :, 2, 1:n_dusttype, 1:n_asize, 5)) &
                 + M_larger_grains(1:n_dusttype))/ sum(data_dust_main(:, :, :, 1, 1:n_dusttype, 1:n_asize, 5))),&
                 " %." //char(27)//'[A'
        ELSE        
             write (*,'(A,I3,A, F5.1,A,F6.2,A)') '  2. Processing: ', int((i_timepoints-1)/real(n_timepoints-1)*100.0),&
                 & ' % done. Est. remain. time: ', t_remain/3600.0, ' h.   Survival: ',&
                 (100.0*(sum(data_dust_main(:, :, :, 2, 1:n_dusttype, 1:n_asize, 5)) &
                 + M_larger_grains(1:n_dusttype))/ sum(data_dust_main(:, :, :, 1, 1:n_dusttype, 1:n_asize, 5))),&
                 " %." //char(27)//'[A'
        END IF  

        !==================================================  
    END DO ! end i_timepoints = 2, n_timepoints 
  
    !==================================================   
    ! Display final time
    write (*,'(A,I3,A,F6.2,A)') '  2. Processing: ', nint(100.0), ' % done.                    Final dust survival: ',   &
        (100.0*(sum(data_dust_main(:,:,:,2, 1:n_dusttype, 1:n_asize, 5)) + M_larger_grains(1:n_dusttype))  &
            /sum(data_dust_main(:, :, :, 1, 1:n_dusttype, 1:n_asize, 5))), ' %.'  //char(27)//'[A'
  !==================================================                 
  
end subroutine run_paperboats ! #s2 

!##################################################################################################################################
!##################################################################################################################################

subroutine cell_values(i_xgrid, i_ygrid, i_zgrid, Temp_gas, rho_gas) ! #s2.0
    use omp_lib
    use datatype
    use variable

    implicit none  
    
    INTEGER, intent(in) :: i_xgrid, i_ygrid, i_zgrid
    REAL (kind=r2), intent(out) :: Temp_gas, rho_gas

    rho_gas  = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,7)    ! in kg/m**(-3) Density
    Temp_gas = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,8)    ! in K

    ! Call the routine gas_ionisation. From now on, the values for the ionisation are stored.   
    call gas_ionisation(i_xgrid, i_ygrid, i_zgrid, Temp_gas)    

end subroutine cell_values ! #s2.0

!##################################################################################################################################
!##################################################################################################################################

subroutine gas_ionisation(i_xgrid, i_ygrid, i_zgrid, Temp_gas) ! #s2.0.1
    use omp_lib
    use datatype
    use variable

    implicit none   
    
    INTEGER,        intent(in) :: i_xgrid, i_ygrid, i_zgrid
    REAL (kind=r2), intent(in) :: Temp_gas
    
    INTEGER :: i, var_h    
    REAL  (kind=r2):: c1fit, c2fit, c3fit, c4fit, c5fit, c6fit, c7fit, c8fit, c9fit, c10fit, c11fit, c12fit
    REAL  (kind=r2):: f1fit, f2fit, f3fit, f4fit
    
    !For partially ionised gas-elements, we should treat the gas-elements as different species (the one with ionised, and the one with unionised gas, even for the same chemical element). For example, H, H+, He and He+ at the same time (!) are all individual gas elements.        
    Do i = 1, n_gas_elements
        var_h = nint(Z_gas(i))
        If (var_h==8) THEN
            c1fit =    3.3448266983032227     
            c2fit =    2.9100341042063804E-029
            c3fit =    24.530826568603516     
            c4fit =    7.2129998207092285     
            c5fit =    1.5999999674579037E-030
            c6fit =    27.353000640869141     
            c7fit =    3.2454819679260254     
            c8fit =    1.0843657264490041E-021
            c9fit =    19.567268371582031     
            c10fit =    198827.59068833734     
            c11fit =    1.2886632265990414E-003
            c12fit =    4.0653022647928454   
            
            f1fit = max(0.0, 2.0 - c1fit*exp(-c2fit* (log(Temp_gas))**c3fit))
            f2fit = max(0.0, 2.0 - c4fit*exp(-c5fit* (log(Temp_gas))**c6fit))
            f3fit = max(0.0, 3.0 - c7fit*exp(-c8fit* (log(Temp_gas))**c9fit)) 
            !f4fit = max(0.0, 1.0 - c10fit*exp(-c11fit* (log(Temp_gas))**c12fit))  
            
            f4fit = 1.0 ! photoemission is causing free electrons low temperatures (MJB)
            
            Ion_grade(i_xgrid, i_ygrid, i_zgrid,i)    = f1fit+f2fit+f3fit+ f4fit    
        ELSE
            Ion_grade(i_xgrid, i_ygrid, i_zgrid,i)    = 1.0 ! Each gas element is 100% single-ionized.
        END IF 
    END DO
    Ion_grade_average(i_xgrid, i_ygrid, i_zgrid) = sum(Ion_grade(i_xgrid, i_ygrid, i_zgrid,:))/(1.0 * n_gas_elements) 
 
end subroutine gas_ionisation ! #s2.0.1

!################################################################################################################
!################################################################################################################

subroutine velocity_calculation(i_xgrid,i_ygrid,i_zgrid, Temp_gas, rho_gas) ! #s2.1
  use omp_lib
  use datatype
  use variable

  implicit none

  INTEGER,       intent(in) :: i_xgrid, i_ygrid, i_zgrid          ! input  
  REAL(kind=r2), intent(in) :: Temp_gas, rho_gas                  ! input
  Real(kind=r2), DIMENSION(:), ALLOCATABLE ::  v_gas(:), v_drag(:)
  Real(kind=r2), DIMENSION(3) :: v_tra, v_act
  Real(kind=r2) :: Phase_Phi, Angular
  
  INTEGER :: i_timeinter,  i_dusttype, i_asize
  
  ALLOCATE(v_gas(1:3), v_drag(1:3))

  v_gas(1) = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4)    ! in m/s
  v_gas(2) = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5)    ! in m/s   
  v_gas(3) = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6)    ! in m/s 
  
  Do i_dusttype = 1, n_dusttype
    Do i_asize = 0, n_asize 
   
        IF (i_asize > 0) THEN
        
            v_drag = 0.0
            
            !$omp CRITICAL
            v_dust(1:3,i_xgrid,i_ygrid,i_zgrid,:,i_dusttype, i_asize) = 0.0
            
            v_dust(1:3,i_xgrid,i_ygrid,i_zgrid,0, i_dusttype, i_asize) = &
                data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i_dusttype, i_asize,1:3) ! vx, vy, vz in m/s, Initialise
                
            IF (l_MHD == 1) THEN 
                v_dust_MHD_act(1:3,i_xgrid,i_ygrid,i_zgrid,:,i_dusttype, i_asize) = 0.0
                v_dust_MHD_tra(1:3,i_xgrid,i_ygrid,i_zgrid,:,i_dusttype, i_asize) = 0.0
                Phase_Phi = 0.0 !          
                ! Calculate v_dust_MHD_act and v_dust_MHD_tra for each inter-timestep i_timeinter = 1, n_timeinter on the basis of v_dust(i_timeinter), 
                ! The value at i_timeinter=0 is not important.
            END IF
            !$omp end CRITICAL
            
           !====================================================================================================================================== 
           ! Here starts the the subdivion of the time_interval (length of one i_timesteps) into n_timeinter time_steps (lengths: Delta_t_inter)
           !i_timeinter = 0, n_timeinter. i_timeinter = 0 is the initial dust velocity, calculated at the timepoint before.        
            DO i_timeinter = 1, n_timeinter

                !____________________________            
                ! 1. Gas drag (collisions and plasma drag)
                call gas_drag(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize,&
                    v_dust(1,i_xgrid,i_ygrid,i_zgrid,i_timeinter-1, i_dusttype, i_asize), &
                    v_dust(2,i_xgrid,i_ygrid,i_zgrid,i_timeinter-1, i_dusttype, i_asize), &
                    v_dust(3,i_xgrid,i_ygrid,i_zgrid,i_timeinter-1, i_dusttype, i_asize), &
                    v_gas(1), v_gas(2), v_gas(3), Temp_gas, rho_gas, v_drag(1), v_drag(2), v_drag(3))             

                !___________________________                    
                ! 2. Lorentz force                   
                IF (l_MHD ==  1) THEN                   
                    call Lorentz_force(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize,-v_drag, Phase_Phi, Angular, v_tra, v_act)
                    Phase_Phi = Phase_Phi + Angular
                    IF (Phase_Phi < 0) THEN
                        Phase_Phi= modulo(Phase_Phi, -2.0_r2*pi)
                    ELSE
                        Phase_Phi= modulo(Phase_Phi, 2.0_r2*pi)                    
                    END IF
                    
                END IF
                !___________________________
                ! 3. Result: dust velocities
                !$omp FLUSH (v_dust, v_dust_MHD_tra, v_dust_MHD_act)
                v_dust(1:3,i_xgrid,i_ygrid,i_zgrid,i_timeinter, i_dusttype, i_asize) = v_gas(1:3) - v_drag(1:3)     ! in m/s, vector, v_drag_a is getting smaller and smaller 
                IF (l_MHD ==  1) THEN                  
                    v_dust_MHD_tra(:,i_xgrid,i_ygrid,i_zgrid,i_timeinter, i_dusttype, i_asize) = v_tra + v_gas               
                    v_dust_MHD_act(:,i_xgrid,i_ygrid,i_zgrid,i_timeinter, i_dusttype, i_asize) = v_act + v_gas
                END IF                
                !$omp FLUSH (v_dust, v_dust_MHD_tra, v_dust_MHD_act)     
            END DO            
            !======================================================================================================================================                     
           
            ! Save it in the data_dust_next array, data only for the NEXT timepoint.
            !$omp CRITICAL
!             data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, 1:3) = &
!                 v_dust(1:3,i_xgrid,i_ygrid,i_zgrid,n_timeinter, i_dusttype, i_asize)  
                
            data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, 4) = &
                data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i_dusttype, i_asize, 4) 
            !$omp end CRITICAL
        
        
        ELSE ! i_asize = 0.0
            !$omp CRITICAL   
!             data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, 1:3) = v_gas  
            data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, 4:5) = &    
                data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i_dusttype, i_asize, 4:5)            
            !$omp end CRITICAL        
        END IF  
    
     END DO
   END DO

   DEALLOCATE(v_gas, v_drag)     
   
end subroutine velocity_calculation ! #s2.1 

!################################################################################################################
!################################################################################################################

subroutine gas_drag(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, &
        vx_dust, vy_dust, vz_dust, vx_gas,vy_gas,vz_gas, Temp_gas, rho_gas, vx_drag, vy_drag, vz_drag) ! #s2.1.1
    use omp_lib
    use datatype
    use variable
    implicit none
    
    INTEGER,       intent(in) :: i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize                      ! input   
    Real(kind=r2), intent(in)  :: vx_dust, vy_dust,vz_dust, vx_gas,vy_gas,vz_gas, Temp_gas, rho_gas  ! input
    Real(kind=r2), intent(out) :: vx_drag, vy_drag, vz_drag                                          ! output
 
    Real(kind=r2) :: n_gas, beta, beta_elek, v_rel, v_rel0, F_s, a_acc, t_run, t_step_step
    Real(kind=r2) :: S, S_elek, collision_drag, coulomb_drag
    
    n_gas = rho_gas/(mu_chem * m_amu)       ! in 1/m**(-3) Density

    beta = (mu_chem * m_amu - Ion_grade_average(i_xgrid, i_ygrid, i_zgrid) * m_electron)&
        * (2.0 * kB * Temp_gas)**(-1.0)   ! beta = (1/(thermal velocity))**2.0  
    beta_elek =  m_electron * (2.0 * kB * Temp_gas)**(-1.0)  
    
    v_rel = ((vx_dust -vx_gas)**2.0 + (vy_dust-vy_gas)**2.0 + (vz_dust -vz_gas)**2.0)**0.5  ! in m/s
    v_rel0 = v_rel ! Inititial v_rel0, in m/s  

    ! The aim is to reduce now v_rel step by step to zero
     IF (v_rel < 0.1) THEN     ! Dont remove. Otherwise: (v_rel -> 0) -> (S -> 0) -> (a_acc -> infinity)
          vx_drag = 0.0
          vy_drag = 0.0
          vz_drag = 0.0
     Else   
         ! Calculate t_step_step, using a first step.
         S = (beta)**0.5 * (v_rel)
         S_elek = (beta_elek)**0.5 * (v_rel)            
         
         !###################################################################         
         !### There are two components that might contribute to the gas drag: 
         !### The drag due to collisions with the ions, and plasma drag
         ! 1. Drag due to collisions
         call coll_drag(S, S_elek, Ion_grade_average(i_xgrid, i_ygrid, i_zgrid), collision_drag)
         ! 2. Drag due to charged gas and dust            
         IF (n_plasma == 0) THEN
            coulomb_drag = 0.0
         ELSE
            call plasma_drag(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, S, S_elek, Temp_gas, rho_gas, coulomb_drag)
          END IF
         !################################################################### 
 
         ! Calculate the force F_s and the acceleration a_acc acting on the dust grains
         F_s = 2.0 * pi**0.5 * n_gas * a(i_asize,i_dusttype)**2.0 * kB * Temp_gas * &
                 (collision_drag + coulomb_drag)
         
         a_acc = F_s/(4.0/3.0*pi * a(i_asize,i_dusttype)**3.0 * rho_bulk(i_dusttype))
 
         t_run = 0.0
     
         Do while ( (t_run <= Delta_t_inter-0.5) .and. (v_rel >= 0.1) .and. (a_acc>0))
             t_step_step = min(0.2*abs(v_rel/a_acc), 0.2*Delta_t_inter)
             
             t_run = t_run + t_step_step
 
             S = (beta)**0.5 * (v_rel)
             S_elek = (beta_elek)**0.5 * (v_rel)            
             
             ! 1. Drag due to collisions
             call coll_drag(S, S_elek, Ion_grade_average(i_xgrid, i_ygrid, i_zgrid), collision_drag)
 
             ! 2. Drag due to charged gas and dust            
             coulomb_drag = 0.0
              IF (n_plasma > 0) THEN
                  call plasma_drag(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, S, S_elek, Temp_gas, rho_gas, coulomb_drag)
              END IF
 
             ! 3. Calculate the force and acceleration:                                         
             F_s = 2.0 * pi**0.5 * n_gas *a(i_asize,i_dusttype)**2.0 * kB * Temp_gas * &
                 (collision_drag + coulomb_drag)
             a_acc = F_s/(4.0/3.0*pi * a(i_asize,i_dusttype)**3.0 * rho_bulk(i_dusttype))
    
             if (a_acc > 0 ) THEN
                 v_rel = v_rel - a_acc * t_step_step
             END IF
         END DO     

        if (v_rel < 0.1) THEN
            v_rel = v_rel*t_run/Delta_t_inter  ! Assuminginear decreasing slope
        END IF      

        vx_drag = (vx_gas - vx_dust) * v_rel/v_rel0 ! in m/s
        vy_drag = (vy_gas - vy_dust) * v_rel/v_rel0 ! in m/s
        vz_drag = (vz_gas - vz_dust) * v_rel/v_rel0 ! in m/s
     END IF

end subroutine gas_drag  ! #s2.1.1

!################################################################################################################
!################################################################################################################

subroutine coll_drag(S, S_elek, Io_g_a, collision_drag)  ! #s2.1.1.1
    use omp_lib
    use datatype
    use variable
    implicit none
    
    Real(kind=r2), intent(in)  :: S, S_elek, Io_g_a      ! input
    Real(kind=r2), intent(out) :: collision_drag ! output
    
    Real(kind=r2) :: coll_drag_ion, coll_drag_elek

    ! Contribution from the ions 
    coll_drag_ion = (S + 0.5*S**(-1.0)) *exp(-1.0*S**2.0) + &
        pi**0.5 *(S**2.0 + 1.0 - 0.25*(S**(-2.0)))* erf(S)                      ! Baines et al. (1965)
        
    ! Contribution from the electrons 
    coll_drag_elek = (S_elek + 0.5*S_elek**(-1.0)) *exp(-1.0*S_elek**2.0) + &
        pi**0.5 *(S_elek**2.0 + 1.0 - 0.25*(S_elek**(-2.0)))* erf(S_elek)       ! Baines et al. (1965)

    collision_drag = coll_drag_ion + (Io_g_a * coll_drag_elek)

end subroutine coll_drag  ! #s2.1.1.1

!################################################################################################################
!################################################################################################################

subroutine plasma_drag(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, S, S_elek, Temp_gas, rho_gas, coulomb_drag) ! #s2.1.1.2
    use omp_lib
    use datatype
    use variable
    implicit none

    INTEGER,      intent(in)   :: i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize   ! input      
    Real(kind=r2), intent(in)  :: S, S_elek, Temp_gas, rho_gas    ! input
    Real(kind=r2), intent(out) :: coulomb_drag                    ! output
    
    REAL(kind=r2) :: Z_grain, phi_pl, lambda_pl, n_gas, n_electron
    Real(kind=r2) :: colou_drag_ion, colou_drag_elek    
    
    n_gas = rho_gas/(mu_chem * m_amu)       ! in 1/m**(-3) Density    
    
    !====================================================================
    ! 1. Define the charge number, phi_pl, n_electron and lambda_pl
    Z_grain = data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i_dusttype, i_asize, 6) ! Charge number, dimensionsless
    
    ! Make a case discrimination. The charge value of abs(Z_grain) = 1.0E-6 is arbitrary. For smaller values, ln(lambda_pl) is huge, but 
    ! coulomb_drag is negligible, because coulomb_drag ~ phi_pl**2.0 * ln(lambda_pl) ~ Z_grain**2.0 * ln(lambda_pl).
    ! For Q = 0, ln(lambda_pl) = infinty, because phi_pl = 0.
    
    If (abs(Z_grain) > 1.0E-6) THEN
        phi_pl = (Z_grain * q_elementar**2.0)/ &
            ((4.0 * pi * epsilon_0) * a(i_asize,i_dusttype) * kB * Temp_gas)
            
        IF (Ion_grade_average(i_xgrid, i_ygrid, i_zgrid) > 0.0) THEN
            n_electron = n_gas * Ion_grade_average(i_xgrid, i_ygrid, i_zgrid)   ! electrons per m**3
    
            ! ln(lambda_pl) should be in the range 30 - 40 (Dwek&Arendt1992, -> Spitzer 1962)
            ! Assume: Each ionised gas atom has a charge of +1. No double or higher ionisation states. 
            lambda_pl = (3.0 * (epsilon_0 * kB * Temp_gas)**0.5 )/(a(i_asize,i_dusttype) * &
                abs(phi_pl) * q_elementar * n_electron**0.5)
        
            !====================================================================
            ! 2. Calculate the term coulomb_drag. Electrons and ions
            ! see e.g. Dwek & Arendt (1992), Fry et al. (2018)            
            colou_drag_ion = Ion_grade_average(i_xgrid, i_ygrid, i_zgrid)**2.0 * phi_pl**2.0 * &
                log(abs(lambda_pl/Ion_grade_average(i_xgrid, i_ygrid, i_zgrid))) * &
                (pi**0.5 *erf(S)*S**(-2.0) - 2.0 * exp(-1.0*S**2.0)*S**(-1.0))
              
            colou_drag_elek = phi_pl**2.0 * log(abs(lambda_pl)) * &
                (pi**0.5 *erf(S_elek)*S_elek**(-2.0) - 2.0 * exp(-1.0*S_elek**2.0)*S_elek**(-1.0))              
              
            coulomb_drag = colou_drag_ion + (Ion_grade_average(i_xgrid, i_ygrid, i_zgrid) * colou_drag_elek)
            !print*, colou_drag_ion, colou_drag_elek, (Ion_grade_average(i_xgrid, i_ygrid, i_zgrid) ), Temp_gas
        ELSE
            coulomb_drag = 0.0        
        END IF
    ELSE
        coulomb_drag = 0.0
    END IF   

end subroutine plasma_drag  ! #s2.1.1.2

!################################################################################################################
!################################################################################################################

subroutine Lorentz_force(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, v_dust_here, Phase_phi, Angular, v_tra, v_act)  ! #s2.1.1.1
    use omp_lib
    use datatype
    use variable
    implicit none
    
    INTEGER,       intent(in) :: i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize                      ! input   
    Real(kind=r2), intent(in) ::  Phase_phi
    Real(kind=r2), intent(in), dimension(3) :: v_dust_here                                           ! input  This is the dust velocity relative to the gas velocity, so v_dust - v_gas. 
                                                                                                     ! Why "-v_gas"? Because the magnetic field moves with the gas
    Real(kind=r2), intent(out), dimension(3) :: v_tra, v_act                                                                                                    
    Real(kind=r2), intent(out) ::  Angular
    
    REAL(kind=r2), dimension(3) :: B_here, v_para, v_perp, u_perp, Delta_r, v_perp_unit, u_perp_unit
    REAL(kind=r2) :: Z_grain, R_gyro, Omega, Spin 

    LOGICAL :: Condition_Lorentz    
      
    !====================================================================
    ! 1. Define the charge number and the magnetic field
    Z_grain = data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i_dusttype, i_asize, 6) ! Charge number, dimensionsless
    
    B_here(1)  = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,9)      ! in Tesla
    B_here(2)  = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,10)     ! in Tesla
    B_here(3)  = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,11)     ! in Tesla
    
    !====================================================================    
    ! 2a. Is there a Lorentz force?    
 
    Condition_Lorentz = .false.
    IF (abs(Z_grain) * f_norm(B_here) * f_norm(v_dust_here) > 0.0) THEN
        IF (abs(dot_product(v_dust_here, B_here))  < (f_norm(B_here) * f_norm(v_dust_here))) THEN
                       ! Otherwise, abs(dot_product(v_dust_here, B_here)) is equal to (f_norm(B_here) * f_norm(v_dust_here)),
                       ! which means that v_dust_here and B_here are parallel or antiparallel, and no Lorentz force is acting on the dust grains.       
            Condition_Lorentz = .true.
        END IF
    END IF   
    
    IF (Condition_Lorentz) THEN
        
        ! Component of v_dust_here that is parallel to B_here
        v_para = (dot_product(v_dust_here, B_here)/(f_norm(B_here))**2.0)   * B_here
        
        ! Component of v_dust_here that is perpendicular to B_here    
        v_perp = v_dust_here - v_para          

        ! Vector that is perpendicular to v_para AND v_perp. At the same time, u_perp is parallel to the Lorentz acceleration.          
        u_perp = cross_product(B_here, v_perp)
        
        ! v_perp and u_perp are the basis vectors that define the Gyration-circle.  
        ! v_perp_unit and u_perp_unit are the normalised basis vectors that define the Gyration-circle with radius 1.  
        v_perp_unit = v_perp/(f_norm(v_perp))  
        u_perp_unit = u_perp/(f_norm(u_perp))  
        
        
        ! gyration frequency, gyration radius, and rotation direction (left or right)   
        Omega  = abs(Z_grain) * q_elementar * f_norm(B_here) / (4.0/3.0 * pi* a(i_asize, i_dusttype)**3.0 * rho_bulk(i_dusttype))     
        R_gyro = f_norm(v_perp)/Omega    
        Spin   = sign(1.0_r2, Z_grain)
        Angular = Omega * Delta_t_inter * Spin
        
        Delta_r =  Spin  * (sin(Phase_phi + Angular) - sin(Phase_phi)) * R_gyro  * v_perp_unit + &
                   Spin  * (cos(Phase_phi + Angular) - cos(Phase_phi)) * R_gyro  * u_perp_unit + &
                            v_para *  Delta_t_inter
                      
        !====================================================================         
        v_tra =      Delta_r / Delta_t_inter
        
        v_act =          f_norm(v_perp) * cos(Phase_phi + Angular) * v_perp_unit + &
                 (-1.0)* f_norm(v_perp) * sin(Phase_phi + Angular) * u_perp_unit + &
                v_para
        !====================================================================          

    ELSE
    ! 2b. If there is no Lorentz force:          
        !====================================================================     
         v_tra = v_dust_here
         v_act = v_dust_here
        !====================================================================            
    END IF
 
end subroutine Lorentz_force

!################################################################################################################
!################################################################################################################

subroutine coulomb_explosions(i_xgrid, i_ygrid, i_zgrid) ! #s2.2
    use omp_lib
    use datatype
    use variable
    use tools

    implicit none
  
    INTEGER, intent(in) :: i_xgrid, i_ygrid, i_zgrid   ! input      
    INTEGER ::             i_dusttype, i_asize 
    REAL (kind =r2) ::     Stress_ce, Z_grain, n_exploded, n_produced

    Do i_dusttype = 1, n_dusttype
        Do i_asize = 1, n_asize
            !====================================================================================================
            ! 1. Calculate the stress tensile using the grain charge number and the grain radius.
            ! The stress is a force per area. Therefore, the Coulomb-stress on a sphere is the Coulomb-force divided by the sphere-surface. In SI-units:

            Z_grain = data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i_dusttype, i_asize, 6) ! Charge number, dimensionsless
            Stress_ce = ((Z_grain * q_elementar)**2.0_r2)/&
                        ((4.0_r2 * pi * epsilon_0) * 4.0_r2 * pi* (a(i_asize,i_dusttype)**4.0_r2)) ! Draine et al (1977,ApJ,231,77, Sects. 3d and 4g) Tazaki et al. (2020), equation (5)
                        
            !====================================================================================================
            ! 2. If the tensile stress is larger than the tensile strength, then the grains are destroyed.
            if (Stress_ce > Stress_Thres(i_dusttype)) Then  

                !============================================== 
                ! Number of exploding dust grains and number of produced particles in the "dusty gas" 
                n_exploded = data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, 4)              
                n_produced = (4.0_r2/3.0_r2 * pi * rho_bulk(i_dusttype)) * &
                    (n_exploded * a(i_asize,i_dusttype)**3.0_r2)/m_grainatom(i_dusttype)
                !============================================== 
                
                !================================================  
                ! Update data_dust_next 
                !$omp FLUSH (data_dust_main, data_dust_next) 
                data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, 4) = &
                    data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, 4) - n_exploded
                !$omp FLUSH (data_dust_main, data_dust_next) 
                data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, 0, 4) = &
                    data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, 0, 4)       + n_produced   
                !================================================  
                
                !================================================
                !Calculate the mass of the "dusty gas"
                !$omp FLUSH (data_dust_main, data_dust_next) 
                data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, 0, 5) = &
                    data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, 0, 4) * m_grainatom(i_dusttype)          
                !======================================== 

            end if
            !====================================================================================================
 
        END DO
    END DO   

end subroutine coulomb_explosions ! #s2.2    

!################################################################################################################
!################################################################################################################

subroutine grain_grain_collisions(i_xgrid, i_ygrid, i_zgrid) ! #s2.2
  use omp_lib
  use datatype
  use variable
  use tools

  implicit none
  
  INTEGER, intent(in)   :: i_xgrid, i_ygrid, i_zgrid   ! input      
  
  INTEGER :: j_asize, j_dusttype, i_max, i_min, i_rem, i_floor, i, k1, k2, i_timeinter, i_dusttype, i_asize
  
  REAL (kind=r2) :: v_col, alpha_1, alpha_0, P_ab, amax_shat, a_new, a_search, amin_shat, a_remnant
  REAL (kind=r2) :: sigma, sum_alpha, n_A, tau, P_left_bin, P_left_bin1, P_left_bin2, P_left_bin3
  REAL (kind=r2) :: n_collide, E_kin_gg_0, E_pot_gg, Rij, v_Thres_boun, m_save
  REAL (kind=r2) :: R_Hir, phi_Hir, M1_Hir, Mr_Hir, sigma1, sigma1i, M_part, m_frag, sti_grains
  REAL (kind=r2) :: a_pro, rho_pro, c_Pro, s_Pro, Pl_Pro, Pv_Pro, a_Tar, rho_Tar, c_Tar, s_Tar, Pl_Tar, Pv_Tar
  REAL (kind=r2) :: Ecol, Ebin_i, Ebin_j, Ngrain_i, Ngrain_j, N_surv_i, N_surv_j
  Real (kind=r2), DIMENSION(:), ALLOCATABLE :: c_h, M_lg_here_change, M_lg_here_init
  Real (kind=r2), DIMENSION(:,:), ALLOCATABLE ::   data_dust_here_change, data_dust_here_init
  
  LOGICAL :: Condition_vapo, Condition_frag, Condition_bounce, Condition_stick
  
  Do i_dusttype = 1, n_dusttype
    Do i_asize = 1, n_asize
  
  ALLOCATE(data_dust_here_change(1:n_dusttype, 0:n_asize), data_dust_here_init(1:n_dusttype, 0:n_asize),&
    M_lg_here_change(1:n_dusttype), M_lg_here_init(1:n_dusttype))
  
  DO  j_dusttype = 1, n_dusttype
        DO j_asize = 0, n_asize
            data_dust_here_init(j_dusttype,j_asize) = data_dust_next(i_xgrid, i_ygrid, i_zgrid, j_dusttype, j_asize, 4)
            data_dust_here_change(j_dusttype, j_asize)      = data_dust_here_init(j_dusttype,j_asize)
        END DO
        M_lg_here_init(j_dusttype) = M_larger_grains(j_dusttype)
        M_lg_here_change(j_dusttype)       = M_lg_here_init(j_dusttype)
  END DO  

  !========================================
  ! Calculate the probability of collisions
        
  DO i_timeinter = 1, n_timeinter ! Divide Delta_t into n_timeinter time steps (lengths: Delta_t_inter)  
    Do j_dusttype = 1, n_dusttype ! Projectile-material
        DO j_asize = 1, n_asize   ! Projectile-size

            !================================================================================================         
            IF (l_MHD ==  0) THEN        
            v_col = ( (v_dust(1,i_xgrid,i_ygrid,i_zgrid,i_timeinter, i_dusttype, i_asize) - &
                v_dust(1,i_xgrid,i_ygrid,i_zgrid,i_timeinter, j_dusttype, j_asize))**2.0_r2 + &
                (v_dust(2,i_xgrid,i_ygrid,i_zgrid,i_timeinter, i_dusttype, i_asize) - &
                v_dust(2,i_xgrid,i_ygrid,i_zgrid,i_timeinter, j_dusttype, j_asize))**2.0_r2 + &
                (v_dust(3,i_xgrid,i_ygrid,i_zgrid,i_timeinter, i_dusttype, i_asize) - &
                v_dust(3,i_xgrid,i_ygrid,i_zgrid,i_timeinter, j_dusttype, j_asize))**2.0_r2)**0.5_r2 ! in m/s 
               
           ELSE    
            v_col = ( (v_dust_MHD_act(1,i_xgrid,i_ygrid,i_zgrid,i_timeinter, i_dusttype, i_asize) - &
                v_dust_MHD_act(1,i_xgrid,i_ygrid,i_zgrid,i_timeinter, j_dusttype, j_asize))**2.0_r2 + &
                (v_dust_MHD_act(2,i_xgrid,i_ygrid,i_zgrid,i_timeinter, i_dusttype, i_asize) - &
                v_dust_MHD_act(2,i_xgrid,i_ygrid,i_zgrid,i_timeinter, j_dusttype, j_asize))**2.0_r2 + &
                (v_dust_MHD_act(3,i_xgrid,i_ygrid,i_zgrid,i_timeinter, i_dusttype, i_asize) - &
                v_dust_MHD_act(3,i_xgrid,i_ygrid,i_zgrid,i_timeinter, j_dusttype, j_asize))**2.0_r2)**0.5_r2 ! in m/s                 
            END If
            !================================================================================================                 

            sigma = pi * ((a(i_asize,i_dusttype) + a(j_asize, j_dusttype)))**2.0_r2

            !==========================================================================
            ! Consider repulsion of charges, approach of Rutherford scattering
            if ((n_gg_coulomb == 1) .and. (v_col > 0.0)) then 
                ! 1. Calculate the energy of the collision (without grain charges, simple E = 0.5 * m1*m2/(m1+m2) * v**2)
                E_kin_gg_0 = 0.5_r2 * v_col**2.0_r2 * a(i_asize, i_dusttype)**3.0_r2 * a(j_asize, j_dusttype)**3.0_r2 *&
                    rho_bulk(i_dusttype) * rho_bulk(j_dusttype) * &
                    (rho_bulk(i_dusttype) * a(i_asize, i_dusttype)**3.0_r2 + &
                    rho_bulk(j_dusttype) * a(j_asize, j_dusttype)**3.0_r2)**(-1.0_r2)
                
                ! 2. Calculate the pot. energy E_pot_gg = Z1 * Z2 * e**2 /(4 * pi * epsilon_0) * 1.0/ (a(i_asize,i_dusttype) + a(j_asize, j_dusttype)) at the distance (a1 + a2) (minimum distance) 
                E_pot_gg = data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i_dusttype, i_asize, 6) * &    ! in J
                    data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, j_dusttype, j_asize, 6) * &
                    q_elementar**2.0/(4.0 * pi * epsilon_0 * (a(i_asize,i_dusttype) + a(j_asize, j_dusttype)))
                
                ! 3. Using 1. and 2., calculate the ratio between potential and kinetic energy
                alpha_1 = E_pot_gg/(E_kin_gg_0)
                
                if (alpha_1 <= 1.0) then                ! Otherwise, the Coulomb-force prevents a collision   
                    sigma = sigma * (1.0 - alpha_1) ! Rutherford scattering experiment
                    v_col = v_col * (1.0 - alpha_1)**0.5   
                end if
            end if
            !==========================================================================

            n_A   = data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, j_dusttype, j_asize, 4) * v_col * Delta_t_inter ! number of particles per area          
            tau = n_A * sigma
            
            ! Calculation of the collision probability:
            ! In general, this would be a summation of the cross sections of all grains, sigma * N                     
            ! P_ab_old = pi * ((a(i_asize,i_dusttype) + a(j_asize, j_dusttype)))**2.0 * v_col * Delta_t_inter * data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, j_dusttype, j_asize, 4)
            ! However, this formula might result in values P_ab >~ 1. In reality, self-shielding of the grains occurs, which limits P_ab to values <1.
            ! This can be realised by multiplication of each single dust grain cross section in an area A, with n_A as the number density per area.
            ! P_ab = 1.0 - (1- sigma/A)**(N)
            ! We can use the general binomical formula to convert the expression in the brackets. This results in:
            
            P_ab = 1.0_r2 - exp(-1.0_r2 * tau)

            ! tau represents the number of particles, that are (in average) at each position of the area. 
            ! The equation is analog to the optical depth in radiative transfer simulations. High optical depths (tau >> 1) represent a high probability for collision (scattering), and low optical depths (tau << 1) a low probability.
            ! For tau << 1, exp(-1.0_r2 * tau) =~ tau, and it follows P_ab = 1 - tau = 1 - n_A * sigma = P_ab_old.
            !=========================
            IF ((P_ab < 0.0) .or. (P_ab > 1.0)) THEN
                print*, "Paperboats ERROR 4, (P_ab < 0) or (1 < P_ab):", P_ab, a(i_asize, i_dusttype), a(j_asize,j_dusttype),tau
                STOP
            END IF
            !=========================
            
            n_collide = data_dust_here_change(i_dusttype, i_asize) * P_ab ! number of collding particles.

            ! IF collisions occur, then evaluate them:
            If ((n_collide > 0.0_r2)) THEN
            
            
                !======================================  
                ! Define conditions for vaporization, fragmentation, bouncing, and sticking
                IF ((v_Thres_vapo(i_dusttype) <= v_col) .and. ((n_gg == 1) .or. (n_gg == 2))) THEN
                       Condition_vapo = .true.
                else
                	Condition_vapo = .false.                        
                end if    
                
                IF ((v_Thres_shat(i_dusttype) <= v_col) .and. (v_col < v_Thres_vapo(i_dusttype)) .and. &
                    & ((n_gg == 1) .or. (n_gg == 3))) THEN
                       Condition_frag = .true.
                else
                	Condition_frag = .false.                        
                end if 
                
                
                ! Caclulate v_Thres_boun (also required for the sticking case)
                
                                
                 Rij = (a(i_asize, i_dusttype) * a(j_asize, j_dusttype)) / (a(i_asize, i_dusttype) + a(j_asize, j_dusttype)) 
                 
                 v_Thres_boun = 21.4_r2 * &
                     (((a(i_asize, i_dusttype)**3.0_r2 + a(j_asize, j_dusttype)**3.0_r2)**0.5_r2) /&
                     ((a(i_asize, i_dusttype) + a(j_asize, j_dusttype))**1.5_r2)) * (gamma_surf(i_dusttype)**(5.0_r2/6.0_r2))/&
                      ((E_Possion(i_dusttype))**(1.0_r2/3.0_r2) * Rij**(5.0_r2/6.0_r2) * (rho_bulk(i_dusttype))**(0.5_r2))                 
                IF ((v_col < v_Thres_shat(i_dusttype)) .and. &
                     ((v_Thres_boun <= v_col) .or. ((i_dusttype .ne. j_dusttype)))) THEN
                       Condition_bounce = .true.
                else
                	Condition_bounce = .false.                        
                end if                                            
            
                IF ((0.0 < v_col) .and. (v_col < v_Thres_boun) .and. (i_dusttype == j_dusttype) .and. (n_gg == 1)) THEN
                       Condition_stick = .true.
                else
                	Condition_stick = .false.                        
                end if 
                !======================================                                                               

                !====================================================================================================
                ! 1. Vaporization

                if (Condition_vapo) then              
!                   count_colli(i_timepoints,1) = count_colli(i_timepoints,1) + n_collide ! counter for vaporization                
!                   data_dust_here_change(i_dusttype, i_asize) =  data_dust_here_change(i_dusttype, i_asize) - n_collide            ! remove particles from the bin i_asize. 
                   
                   Ecol = 2.0/3.0*pi*v_col**2.0 *&
                   	& (rho_bulk(i_dusttype) * a(i_asize, i_dusttype)**3.0 * rho_bulk(j_dusttype) * a(j_asize, j_dusttype)**3.0)/&
                   	& (rho_bulk(i_dusttype) * a(i_asize, i_dusttype)**3.0 + rho_bulk(j_dusttype)* a(j_asize, j_dusttype)**3.0)   ! Collision energy, in SI
                   
                   Ebin_i = 0.74 * m_grainatom(i_dusttype) * ev_unit        ! Binding energy of 1 dust grain atom of grain i
                   Ebin_j = 0.74 * m_grainatom(j_dusttype) * ev_unit        ! Binding energy of 1 dust grain atom of grain j
                   Ngrain_i = (4.0/3.0* pi * rho_bulk(i_dusttype)* a(i_asize, i_dusttype)**3.0)/(m_grainatom(i_dusttype)*m_amu) ! Number of dust atoms in grain i
                   Ngrain_j = (4.0/3.0* pi * rho_bulk(j_dusttype)* a(j_asize, j_dusttype)**3.0)/(m_grainatom(j_dusttype)*m_amu) ! Number of dust atoms in grain j
                                     
                   ! Assumption: The energy Ecol is split into two halfs for the grains i and j                   
                   N_surv_j = Ngrain_j - 0.5* Ecol/Ebin_j ! Number of not vaporized atoms of grain j
                   
                   IF (N_surv_j > 0.0) THEN
                   	N_surv_i = Ngrain_i - 0.5* Ecol/Ebin_i  !number of not vaporized atoms of grain i
                   ELSE	
                   	N_surv_i = Ngrain_i - (Ecol - Ebin_j*Ngrain_j)/Ebin_i  !number of not vaporized atoms of grain i if grain j is fully vapoized.                	
                   END IF
                   
                   IF (N_surv_i < 0.0) THEN
                   	count_colli(i_timepoints,1) = count_colli(i_timepoints,1) + n_collide ! counter for vaporization                
                   	data_dust_here_change(i_dusttype, i_asize) =  data_dust_here_change(i_dusttype, i_asize) - n_collide            ! remove particles from the bin i_asize. 
                   ELSE
                   	count_colli(i_timepoints,1) = count_colli(i_timepoints,1) + n_collide ! counter for vaporization                                    
                   
                   	a_new =  (N_surv_i * (m_grainatom(i_dusttype)*m_amu) * 3.0/(4.0*pi* rho_bulk(i_dusttype)))**(1.0_r2/3.0_r2)
 
                     	call Bin_search(i_dusttype, a_new, i_floor, P_left_bin) 

                     	!================================================
                     	! Reduce number of particles that are vaporised 
                     
                     	data_dust_here_change(i_dusttype, i_asize) = data_dust_here_change(i_dusttype, i_asize) - n_collide                   
                  
                        IF (a_new*delta_rad > a(1,i_dusttype)) THEN ! Otherwise, the remaining bit is smaller than the smallest grain size bin.
                            !================================================
                            ! Add particles to the lower grain size boundary: 
                            IF ((P_left_bin * n_collide) > data_accuracy) THEN
                                data_dust_here_change(i_dusttype, i_floor) = data_dust_here_change(i_dusttype, i_floor) +&
                                    P_left_bin * n_collide
                            END IF
                            !================================================
                            ! Add particles to the upper grain size boundary:
                            if ((i_floor < n_asize) .and. (((1.0_r2 - P_left_bin) * n_collide) > data_accuracy)) THEN
                                data_dust_here_change(i_dusttype, i_floor+1) = &
                                data_dust_here_change(i_dusttype, i_floor+1) + (1.0_r2 - P_left_bin) * n_collide
                            END IF  
                            !================================================                            
                     	END IF

                   END IF
                                    
                end if

                !====================================================================================================
                ! 2. Shattering            
   
                if (Condition_frag) then    
                
                    data_dust_here_change(i_dusttype, i_asize) = data_dust_here_change(i_dusttype, i_asize)   - n_collide
                    count_colli(i_timepoints,2)         = count_colli(i_timepoints,2)           + n_collide
                    
                    a_remnant = 0.0                    
                    if (a(i_asize,i_dusttype) >= a(j_asize,j_dusttype)) THEN ! i is >= j
                        ! Projectile parameter (smaller particle)
                        a_Pro = a(j_asize,j_dusttype)
                        rho_Pro = rho_bulk(j_dusttype)
                        
                        ! Target parameter (larger particle)                     
                        a_Tar   = a(i_asize,i_dusttype)
                        rho_Tar = rho_bulk(i_dusttype)
                        c_Tar   = c_speed(i_dusttype)
                        s_Tar   = s_Hir(i_dusttype)                        
                        Pl_Tar  = Pl_Hir(i_dusttype)
                        Pv_Tar  = Pv_Hir(i_dusttype)

                        R_Hir = ((a_Tar* rho_Tar)/(a_pro* rho_pro))**0.5_r2                      ! Tielens+1994, Eq. (2.13)
                    
                        phi_Hir = Pl_Tar/(rho_Tar*c_Tar**2.0)                                    ! Tielens+1994, Eq. (2.33)
                    
                        M1_Hir = 2.0 * phi_Hir/(1.0 + (1.0 + 4.0 * s_Tar * phi_Hir)**0.5_r2)     ! Tielens+1994, Eq. (2.32)
                                            
                        Mr_Hir = v_col/c_Tar
                    
                        sigma1 = (0.3 * (s_Tar +1.0/M1_Hir - 0.11)**1.3)/(s_Tar+1.0/M1_Hir - 1.0) ! Tielens+1994, Eq. (2.29)
                            
                        sigma1i = (0.3 * (s_Tar+(1.0+R_Hir)/Mr_Hir - 0.11)**1.3)/&
                            (s_Tar+(1.0+R_Hir)/Mr_Hir - 1.0) 
                    
                        ! Mass of particle i, that is shocked.
                        M_part = (4.0/3.0*pi*a_Pro**3.0 * rho_Pro)* (1.0+2.0*R_Hir)& ! Tielens+1994, Eq. (2.25)
                            /(2.0 * sigma1i**(1.0/9.0) * (1.0+R_Hir)**(9.0/16.0)) * &
                            (Mr_Hir/(M1_Hir * sigma1**(0.5)))**(16.0/9.0)
                        
                        if (M_part > (2.0/3.0 * pi * a_Tar**3.0 * rho_Tar)) THEN
                            m_frag =  4.0/3.0 * pi * a_Tar**3.0 * rho_Tar
                            amax_shat = 0.22 * a_tar * (c_Tar/v_col)  * (1.0 + R_Hir) * M1_Hir &
                            * ((a_Tar/a_Pro)**3.0 * (rho_Tar/rho_Pro) * (1.0+2.0*R_Hir)**(-1.0))**(9.0/16.0) &
                            * sigma1**(0.5_r2) * sigma1i**(1.0/16.0)
                            amin_shat = 0.03 * amax_shat 
                        else ! M_part <= (2.0/3.0 * pi * a_Tar**3.0 * rho_Tar) This is the cratering case.
                            m_frag    = 0.4 * M_part
                            a_remnant =  (a(i_asize,i_dusttype)**3.0 - 0.3/pi * M_part/rho_Tar)**(1.0/3.0)
                            
                            IF (a_remnant > a(n_asize,i_dusttype)) THEN
                                print*, 'Paperboats ERROR 8: a_remnant > a_max_total.'
                                stop
                            END IF
                            
                            amax_shat = (m_frag/(64.43839286 * pi * rho_Tar))**(1.0/3.0)
                            amin_shat = amax_shat * (Pl_Tar/Pv_Tar)**1.47
                            
                            call Bin_search(i_dusttype, a_remnant, i_rem, P_left_bin3)                    
                        end if                       
                    else ! a(i_asize,i_dusttype) < a(j_asize,j_dusttype), i is < j.
                        ! Projectile parameter (smaller particle)
                        a_Pro = a(i_asize,i_dusttype)
                        rho_Pro = rho_bulk(i_dusttype)
                        c_Pro   = c_speed(i_dusttype)
                        s_Pro   = s_Hir(i_dusttype)                        
                        Pl_Pro  = Pl_Hir(i_dusttype)
                        Pv_Pro  = Pv_Hir(i_dusttype)                          
                        
                        ! Target parameter (larger particle)                     
                        a_Tar   = a(j_asize,j_dusttype)
                        rho_Tar = rho_bulk(j_dusttype)
                        
                        R_Hir = ((a_Tar* rho_Tar)/(a_pro* rho_pro))**0.5_r2                      ! Tielens+1994, Eq. (2.13)
                        phi_Hir = Pl_Pro/(rho_Pro*c_Pro**2.0) 
                        M1_Hir = 2.0 * phi_Hir/(1.0 + (1.0 + 4.0 * s_Pro * phi_Hir)**0.5_r2)
                        Mr_Hir = v_col/c_Pro
                        sigma1 = (0.3 * (s_Pro +1.0/M1_Hir - 0.11)**1.3)/(s_Pro+1.0/M1_Hir - 1.0)
                        sigma1i = (0.3 * (s_Pro+(1.0+R_Hir)/Mr_Hir - 0.11)**1.3)/&
                            (s_Pro+(1.0+R_Hir)/Mr_Hir - 1.0)
                        
                        m_frag = 4.0/3.0 * pi * a_Pro**3.0 * rho_Pro 
                        
                        amax_shat = 0.22 * a_Pro * (c_Pro/v_col)  * (1.0 + R_Hir) * M1_Hir &
                            * (((a_Pro/a_Tar)**3.0 * (rho_Pro/rho_Tar)) * (1.0+2.0*R_Hir)**(-1.0))**(9.0/16.0) &
                            * sigma1**(0.5_r2) * sigma1i**(1.0/16.0)
                        amin_shat = 0.03 * amax_shat      
                    end if
                    
                    !=======================================================
                    ! Test, if amax_shat smaller than a(i_asize,i_dusttype)
                    If (amax_shat > a(i_asize,i_dusttype)) THEN
                        print*, 'Paperboats ERROR 2211. a_max_shat larger than the initial grain.'                                            
                        print*, amax_shat,">", a(i_asize, i_dusttype)
                        stop
                    END IF                    
                    
                    !==============================================================================================
                    
                    call Bin_search(i_dusttype, amin_shat, i_min, P_left_bin1)
                    call Bin_search(i_dusttype, amax_shat, i_max, P_left_bin2)                                      
                    ! a(i_max) is the largest grain size below amax_shat (or equal to).
                    
                    ! The grains are shattered to a distribution of smaller particles. We assume a powerlaw, n(a(k_asize)) ~ a(k_asize)**(-(3*eta-2))
                    ! Plausible values of η from impact experiments are 1.5, ... 2.0; a “classical” value is 11/6 = 1.83 (corresponds to
                    ! a differential size distribution with the index -(3*(11/6) -2) = 3.5), Ref: e.g. Krivov et al. (2005)
                                        
                    IF (amax_shat > a(1,i_dusttype)/delta_rad) THEN    ! Otherwise, all fragments are smaller than a(1).
                        
                        If (i_min > 0) THEN
                            k1 = i_min
                        ELSE
                            k1 = 1 - ceiling(log(a(1,i_dusttype)/amin_shat)/log(delta_rad))                      
                        END IF
                        
                        !k2 is an integer between 1 and n.
                        k2 = i_max+1

                        ALLOCATE(c_h(k1:k2))
                        sum_alpha = 0.0_r2

                        c_h = 1.0_r2
                        If (k1 == k2-1) THEN ! all the dust is between 2 grain sizes
                            a_search = (amin_shat**3.0_r2 + amax_shat**3.0_r2)**(1.0_r2/3.0_r2)
                            call Bin_search(i_dusttype, a_search, i_floor, P_left_bin1)
                            c_h(k1) = P_left_bin1
                            c_h(k2) = 1.0_r2 - P_left_bin1
                        ELSE
                            c_h(k1) = P_left_bin1
                            c_h(k2) = 1.0_r2 - P_left_bin2
                        END IF

                        sum_alpha = 0.0_r2                    
                        Do i = k1, k2 
                            sum_alpha = sum_alpha + c_h(i) * &
                                (delta_rad**(real(i,kind=r2)-1.0_r2))**(4.0_r2-gamma_asize(i_dusttype))                                 
                        END DO
                        sum_alpha = sum_alpha * (a(1,i_dusttype))**(3.0_r2 - gamma_asize(i_dusttype)) 
                        
                        alpha_0 = n_collide * (a(i_asize, i_dusttype)**3.0_r2-a_remnant**3.0_r2) / sum_alpha   ! n(a(k_asize)) = alpha_0 * a((k_asize))**(-gamma)
                
                        !=====================================================
                        ! Add particles to the grain size bins i, with max(1, k1) =< i <= k2: 
                
                        Do i = max(1, k1), k2
                            data_dust_here_change(i_dusttype, i) = data_dust_here_change(i_dusttype, i) + &
                                alpha_0 * c_h(i) * (a(1,i_dusttype) * delta_rad**(real(i,kind=r2)-1.0_r2))&
                                **(-gamma_asize(i_dusttype)) * delta_rad**(real(i,kind=r2)-1.0_r2)   
                        End Do

                        DEALLOCATE(c_h)

                        !================================================
                        ! Consider the dust remnant (cratering case)
                        IF (a_remnant > a(1,i_dusttype)/delta_rad) THEN
                            ! Add particles to the lower grain size boundary: 
                            data_dust_here_change(i_dusttype, i_rem) = data_dust_here_change(i_dusttype, i_rem) + &
                                P_left_bin3 * n_collide   
                            ! Add particles to the upper grain size boundary:
                            data_dust_here_change(i_dusttype, i_rem+1) = data_dust_here_change(i_dusttype, i_rem+1) + &
                                (1.0_r2 - P_left_bin3) * n_collide
                        END IF
                        !================================================                        
                    end if    
                end if

		!====================================================================================================                  
		! 3. Bouncing
                
               ! No change in the number distribution.
               ! Assume, that the velocity distribution doesn't change. After bouncing of two grains, the velocities of these two might be different, but it is 
               ! rearranged by the gas flow. Otherwise, it would depend on the collision parameters, in particular the collision angle (not only ahead
               ! collisions), which would make it much more complicated.   
 
               ! Counter for bouncing
		If (Condition_bounce) THEN
                     ! The last case only exists because we don't wanna stick grains of different materials.
                     count_colli(i_timepoints,3) = count_colli(i_timepoints,3) + n_collide
		END IF    
        
                 !====================================================================================================            
                 ! 4. Sticking (only for collisions of grains of the same material (i_dusttype == j_dusttype))
                 
                 m_save = 0.0_r2
                 if (Condition_stick) then
                     count_colli(i_timepoints,4) = count_colli(i_timepoints,4) + n_collide ! counter for bouncing

                     a_new = ((a(i_asize, i_dusttype))**3.0_r2 + (a(j_asize, j_dusttype))**3.0_r2)**(1.0_r2/3.0_r2)
 
                     call Bin_search(i_dusttype, a_new, i_floor, P_left_bin3)

                     sti_grains = 0.5_r2 * n_collide
 
                     !================================================
                     ! Reduce number of particles that stick together (with j_dusttype,j_asize) to larger ones
                     
                     data_dust_here_change(i_dusttype, i_asize) =  data_dust_here_change(i_dusttype, i_asize) - sti_grains
                     data_dust_here_change(j_dusttype, j_asize) =  data_dust_here_change(j_dusttype, j_asize) - sti_grains                      
                  
                     !================================================
                     ! Add particles to the lower grain size boundary: 
                     IF ((P_left_bin3 * sti_grains) > data_accuracy) THEN
                        data_dust_here_change(i_dusttype, i_floor) = data_dust_here_change(i_dusttype, i_floor) +&
                            P_left_bin3 * sti_grains
                     END IF
                     
                     !================================================
                     ! Add particles to the upper grain size boundary:

                     if ((i_floor < n_asize) .and. (((1.0_r2 - P_left_bin3) * sti_grains) > data_accuracy)) THEN
                        data_dust_here_change(i_dusttype, i_floor+1) = &
                        data_dust_here_change(i_dusttype, i_floor+1) + (1.0_r2 - P_left_bin3) * sti_grains
                     elseif (i_floor == n_asize ) THEN                               
                        !================================================
                        ! Save the dust grains, which are too large for the size distribution, in M_larger_grains(:)
                        ! These dust is not sputtered or gg collided any more, but also not transported, so it can't move outside the domain.
                        IF (a_new < a(n_asize, i_dusttype) * Delta_rad ) THEN 
                            IF (((1.0_r2 - P_left_bin3) * sti_grains) > data_accuracy) THEN
                                m_save =  (1.0_r2 - P_left_bin3) * sti_grains * 4.0_r2/3.0_r2 * &
                                pi * (a(n_asize, i_dusttype) * delta_rad)**3.0_r2 * rho_bulk(i_dusttype)
                            end if 
                        ELSE ! The new particles are too large
                            m_save =  sti_grains* 4.0_r2/3.0_r2 * pi * a_new**3.0_r2 * rho_bulk(i_dusttype)
                        END IF
                        !$omp FLUSH  
                        M_lg_here_change(i_dusttype) = M_lg_here_change(i_dusttype) + m_save
                        !$omp FLUSH  
                     END if 
                 end if ! End 4. Sticking
                !====================================================================================================                    
                ! End of the evaluation of the cases vaporization, shattering, bouncing and sticking.                
                !====================================================================================================                           
             END IF  ! (n_collide > 0.0_r2)
        END DO       ! END projectile-size of the material j_dusttype
    END DO           ! END projectile-material j_dusttype
  END DO             ! END i_timeinter

  !====================================================================================================  
  ! Calculation of the number of particles in the "dusty gas"
  !$omp CRITICAL    
  data_dust_here_change(i_dusttype, 0) =  data_dust_here_init(i_dusttype, 0)  + &
     ((M_lg_here_init(i_dusttype) - M_lg_here_change(i_dusttype)) + &
     ((4.0_r2/3.0_r2 * pi * rho_bulk(i_dusttype)) * &
      sum((data_dust_here_init(i_dusttype, 1:n_asize)-data_dust_here_change(i_dusttype, 1:n_asize)) *&
      a(1:n_asize,i_dusttype)**3.0_r2))) * 1.0_r2/m_grainatom(i_dusttype)   
   
  DO  j_dusttype = 1, n_dusttype
    DO j_asize = 0, n_asize
      data_dust_next(i_xgrid, i_ygrid, i_zgrid, j_dusttype, j_asize, 4) = &
            data_dust_next(i_xgrid, i_ygrid, i_zgrid, j_dusttype, j_asize, 4) + &
            (data_dust_here_change(j_dusttype, j_asize) - data_dust_here_init(j_dusttype, j_asize))
    END DO
    
    M_larger_grains(i_dusttype) = M_larger_grains(i_dusttype) + &
        (M_lg_here_change(i_dusttype)- M_lg_here_init(i_dusttype))
  END DO
  !$omp end CRITICAL     
 
  DEALLOCATE(data_dust_here_change, data_dust_here_init, M_lg_here_change, M_lg_here_init)

  !================================================
  !Calculate the mass of the "dusty gas"
  !$omp FLUSH  
  data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, 0, 5) = &
        data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, 0, 4) * m_grainatom(i_dusttype)
  !$omp FLUSH          
  !======================================== 
  
    END DO
  END DO

end subroutine grain_grain_collisions ! #s2.2

!################################################################################################################
!################################################################################################################

subroutine sputtering(i_xgrid, i_ygrid, i_zgrid, Temp_gas, rho_gas) ! #s2.3
  use omp_lib
  use datatype
  use variable
  use tools

  implicit none 
  INTEGER, intent(in) :: i_xgrid, i_ygrid, i_zgrid                  !input  
  REAL(kind=r2), intent(in) :: Temp_gas, rho_gas                    !input  
  Real(kind=r2) ::  v_rel, dadt, dadt_rg, a_new, P_sput_left_bin, n_sput, a_without_rg, dummy
  
  INTEGER :: i_lower, i_dusttype, i_asize, j_dusttype, j_asize
  
  Real (kind=r2) , DIMENSION(:), ALLOCATABLE ::     M_lg_here_change, M_lg_here_init, v_gas(:) 
  Real (kind=r2) , DIMENSION(:,:), ALLOCATABLE ::   data_dust_here_change, data_dust_here_init
  Real(kind=r2) :: M_regular_gas_here 
  
  ALLOCATE(v_gas(1:3))

  v_gas(1) = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4)    ! in m/s
  v_gas(2) = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5)    ! in m/s  
  v_gas(3) = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6)    ! in m/s 

  Do i_dusttype = 1, n_dusttype
    Do i_asize = 1, n_asize
  
        ALLOCATE(data_dust_here_change(1:n_dusttype, 0:n_asize), data_dust_here_init(1:n_dusttype, 0:n_asize),&
            M_lg_here_change(1:n_dusttype), M_lg_here_init(1:n_dusttype))  
  
        DO  j_dusttype = 1, n_dusttype
            DO j_asize = 0, n_asize
                data_dust_here_init(j_dusttype,j_asize)     = data_dust_next(i_xgrid, i_ygrid, i_zgrid, j_dusttype, j_asize, 4)
                data_dust_here_change(j_dusttype, j_asize)  = data_dust_here_init(j_dusttype,j_asize)
            END DO
            M_lg_here_init(j_dusttype)   = M_larger_grains(j_dusttype)
            M_lg_here_change(j_dusttype) = M_lg_here_init(j_dusttype)
        END DO 
        
        !=============================================================================================================
            
        IF (l_MHD == 0) THEN ! Calculate the sputtering with the values at i_timeinter = 0
            v_rel =  ((data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i_dusttype, i_asize, 1) - v_gas(1))**2.0 + &
                          (data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i_dusttype, i_asize, 2) - v_gas(2))**2.0 + &
                          (data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i_dusttype, i_asize, 3) - v_gas(3))**2.0)**0.5 ! in m/s         
        ELSE !(l_MHD == 1) THEN ! Calculate the sputtering with the values at i_timeinter = 1
            v_rel = ((v_dust_MHD_act(1,i_xgrid,i_ygrid,i_zgrid,1,i_dusttype, i_asize) - v_gas(1))**2.0 + &
                         (v_dust_MHD_act(2,i_xgrid,i_ygrid,i_zgrid,1,i_dusttype, i_asize) - v_gas(2))**2.0 + &
                         (v_dust_MHD_act(3,i_xgrid,i_ygrid,i_zgrid,1,i_dusttype, i_asize) - v_gas(3))**2.0)**0.5 ! in m/s   
        END IF
 
        !========================
        ! Number conservation by n_sput
        n_sput = data_dust_here_change(i_dusttype, i_asize)            

        call dadt_sputtering(i_xgrid, i_ygrid, i_zgrid,i_dusttype, i_asize, Temp_gas, rho_gas, v_rel, n_sput, dadt, dadt_rg)          
  
        a_new        = a(i_asize, i_dusttype) - dadt * Delta_t ! in m, including regular and dusty gas, accretion and trapping
        a_without_rg = a(i_asize, i_dusttype) - (dadt - dadt_rg) * Delta_t ! in m, without trapping and accretion of regular gas. But the regular gas can still sputter. This quantity is required for reasons of mass conversation.
                
        call Bin_search(i_dusttype, a_new, i_lower, P_sput_left_bin)

        !==============================================
        ! Reduce number of particles that are sputtered   
        data_dust_here_change(i_dusttype, i_asize) = 0.0_r2
  
        !===============================================
        ! Add particles to the lower grain size boundary: 
        IF (P_sput_left_bin * n_sput > data_accuracy) THEN
            if (i_lower > 0) THEN
                data_dust_here_change(i_dusttype, i_lower) = data_dust_here_change(i_dusttype, i_lower) + P_sput_left_bin * n_sput       
            end if
        end if
  
        !================================================
        ! Add particles to the upper grain size boundary:
        IF ((1.0_r2 - P_sput_left_bin * n_sput) > data_accuracy) THEN        
        if (i_lower < n_asize) THEN  
            data_dust_here_change(i_dusttype, i_lower+1) = data_dust_here_change(i_dusttype, i_lower+1) +&
                (1.0_r2 - P_sput_left_bin) * n_sput
        else !i_lower == n_asize
            M_lg_here_change(i_dusttype) = M_lg_here_change(i_dusttype) + &
                (1.0_r2 - P_sput_left_bin) * n_sput * &
                (4.0_r2/3.0_r2*pi*(a(n_asize,i_dusttype)*Delta_rad)**3.0_r2 * rho_bulk(i_dusttype))
        end if
        end if
        !================================================   
        ! Calculation of the number of particles in the "dusty gas"            
        data_dust_here_change(i_dusttype, 0) =  data_dust_here_init(i_dusttype, 0)  + &
            ((M_lg_here_init(i_dusttype) - M_lg_here_change(i_dusttype)) + &
            ((4.0_r2/3.0_r2 * pi * rho_bulk(i_dusttype)) * &
            sum((data_dust_here_init(i_dusttype, 1:n_asize)-data_dust_here_change(i_dusttype, 1:n_asize)) *&
            a(1:n_asize,i_dusttype)**3.0_r2))) * 1.0_r2/m_grainatom(i_dusttype)
            
        ! But a part of the dust is produced from the regular gas! Correct this mass (the number of dusty gas particles):    
        dummy = (4.0_r2/3.0_r2 * pi * rho_bulk(i_dusttype)) * (a_new**3.0_r2-a_without_rg**3.0_r2)  ! Mass [kg] of trapped gas in one grain
        M_regular_gas_here = dummy  * (n_sput)                                                      ! Total mass [kg] of trapped gas per V_cell
        IF (M_regular_gas_here  > rho_gas) THEN
           print*, 'Paperboats ERROR 9. It is more gas per cell trapped/accreted than it contains.'
           print*, "M_regular_gas_here/rho_gas, rho_gas/(mu_chem * m_amu), n_sput, # ,i_asize,&
           &   a(i_asize,i_dusttype), a_new, a_without_rg, #, dadt, dadt_rg, Delta_t  "             
           print*, M_regular_gas_here/rho_gas, rho_gas/(mu_chem * m_amu), n_sput, "#",i_asize,&
           &  a(i_asize,i_dusttype), a_new, a_without_rg, "#", dadt, dadt_rg, Delta_t
           print*, ""           
!            STOP  
        END IF
        data_dust_here_change(i_dusttype, 0) = data_dust_here_change(i_dusttype, 0) + M_regular_gas_here/m_grainatom(i_dusttype) 
    
        !================================================  
        ! Update data_dust_next and M_larger_grains
        
        !$omp CRITICAL         
        DO  j_dusttype = 1, n_dusttype
            DO j_asize = 0, n_asize
            
            !$omp FLUSH (data_dust_main, data_dust_next) 
            data_dust_next(i_xgrid, i_ygrid, i_zgrid, j_dusttype, j_asize, 4) = &
                data_dust_next(i_xgrid, i_ygrid, i_zgrid, j_dusttype, j_asize, 4) + &
                (data_dust_here_change(j_dusttype, j_asize) - data_dust_here_init(j_dusttype, j_asize))
            END DO
    
            M_larger_grains(i_dusttype) = M_larger_grains(i_dusttype) + &
                (M_lg_here_change(i_dusttype)- M_lg_here_init(i_dusttype))
        END DO
        M_regular_gas =  M_regular_gas + M_regular_gas_here         
        !$omp end CRITICAL           

        DEALLOCATE(data_dust_here_change,data_dust_here_init, M_lg_here_change, M_lg_here_init)      
  
        !================================================
        !Calculate the mass of the "dusty gas"
        
        !$omp FLUSH (data_dust_main, data_dust_next) 
        data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, 0, 5) = &
            data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, 0, 4) * m_grainatom(i_dusttype)          
        !======================================== 
        
    END DO
  END DO
  
  DEALLOCATE(v_gas)
  
end subroutine sputtering ! #s2.3

!################################################################################################################
!################################################################################################################

subroutine dadt_sputtering (i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, Temp_gas, rho_gas, v_rel,n_sput, dadt, dadt_rg) ! #s2.3.1
    use omp_lib
    use datatype
    use variable

    implicit none
    
    INTEGER, intent(in)   :: i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize    
    Real(kind=r2), intent(in)  :: Temp_gas, rho_gas, v_rel, n_sput  ! input
    Real(kind=r2), intent(out) :: dadt, dadt_rg                         ! output

    Real(kind=r2), DIMENSION(:), ALLOCATABLE :: n(:), Yield_avg(:), Yield_avg_rg(:), Depl

    Real(kind=r2) :: Yield, Yield_rg, fskm, v, v_before, Delta_v, v_hp, ard, r_d, alpha_1, Z_gas_i, n_gas, E_abs 
    Real(kind=r2) :: alpha_Y, sn_Y, e12_Y, E_threshold_Y, E_Y, g_Y, Zratio_Y, f_size, m_gas_i, Ion_grade_i, tau_Depl, t_ratio 

    INTEGER :: i, j, i_dummy, n_dummy, num_lay, n_v_step_here
    
    ALLOCATE(n(1:n_gas_elements+ n_dusttype), Yield_avg(1:n_gas_elements+ n_dusttype), Yield_avg_rg(1:n_gas_elements+ n_dusttype)) ! + n_dusttype, because each grain species produces dusty gas which can sputter.
    ALLOCATE(Depl(1:n_gas_elements+ n_dusttype))
    
    n_gas = rho_gas/(mu_chem * m_amu)          ! in 1/m**(-3) Density 

    IF ((n_gas_accretion == 0) .or. (n_gas_accretion == 1) .or. (n_gas_accretion == 2)) THEN ! Dusty gas sputters only when regular gas trapping or accretion is off.
        n_dummy = n_gas_elements + n_dusttype ! Do Loop gas_element i  
    ELSE 
        n_dummy = n_gas_elements
    END IF
    
    num_lay = 5 ! number of atomic layers that are required to trap the impinging gas particle
    
    Yield_avg    = 0.0_r2
    Yield_avg_rg = 0.0_r2
    n            = 0.0_r2
    Depl         = 0.0_r2
    
    Do i = 1, n_dummy
   
        IF (i < n_gas_elements+1) THEN ! Regular gas
            n(i)        = n_gas * Abun_gas(i)
            m_gas_i     = m_gas(i)
            Ion_grade_i = Ion_grade(i_xgrid, i_ygrid, i_zgrid, i)
            Z_gas_i     = Z_gas(i)
        ELSE ! Dusty gas
            n(i)        = data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i - n_gas_elements, 0, 4)
            m_gas_i     = m_grainatom(i - n_gas_elements)
            Ion_grade_i = Ion_grade(i_xgrid, i_ygrid, i_zgrid, 1)
            Z_gas_i     = Z_grainatom(i - n_gas_elements)           
        END IF
 
        v_before = 0.0_r2
        
        IF  (n_sp == 1) THEN ! kinetic and thermal sputtering
            v_hp = sqrt(2.0_r2*kB*Temp_gas/m_gas_i + v_rel**2.0) ! velocity with the highest probability in a skewed Maxwellian Distribution. Added the summand v_rel**2.0 on the 11/1/2022 to consider the shift of the maximum due to high v_rel
            n_v_step_here = n_v_step
        else if (n_sp == 2) THEN ! kinetic sputtering
            v_hp = v_rel
            n_v_step_here = 1            
        else if (n_sp == 3) THEN ! thermal sputtering
            v_hp = sqrt(2.0_r2*kB*Temp_gas/m_gas_i)
            n_v_step_here = n_v_step
        END IF
        
        Do j = 1, n_v_step_here ! Do Loop j; n_v_step_here is either 1 or n_v_step (which is fixed in variable.f90)
            
            IF ((n_sp == 1) .or. (n_sp == 3)) THEN
                v = 2.5 * v_hp * (real(j)-0.6)/real(n_v_step-0.6) * 1.1878**((real(j)-0.6)/real(n_v_step-0.6))
                ! for n_v_step = 15, v ranges from 0.07 * v_hp up to 2.97 * v_hp, and for j= 6  the distribution has its maximum.
                Delta_v = v - v_before
            ELSE If (n_sp == 2) THEN ! This is the case: solely kinetic sputtering             
                v = v_hp
                v_before = v ! Needed because we integrate over "0.5_r2 * (v + v_before) * Delta_v" when averaging the Yield                 
                Delta_v  = 1.0 ! Because the distribution is a constant, fskm = 1.0
            END IF
            
            IF (n_sp == 1) THEN           ! kinetic and thermal sputtering 
                IF (v_rel == 0.0_r2) THEN ! Maxwellian distribution
                    fskm =  (m_gas_i/(2.0_r2 * pi * kB * Temp_gas))**1.5_r2 * (4.0_r2 * pi * v**2.0_r2) * &
                        exp(-m_gas_i/(2.0_r2 * kB * Temp_gas) * v**2.0_r2)
                ELSE                      ! skewed Maxwellian distribution       
                    fskm =  (m_gas_i/(2.0_r2 * pi * kB * Temp_gas))**0.5_r2 * (v/v_rel) *  &
                        & (exp(-m_gas_i/(2.0_r2 * kB * Temp_gas) * (v - v_rel)**2.0_r2) - exp(-m_gas_i/ &
                        (2.0_r2 * kB * Temp_gas) * (v + v_rel)**2.0_r2))
                END IF        
            ELSE IF (n_sp == 2) THEN      ! kinetic sputtering 
                fskm = 1.0 ! no distribution, single value.
            ELSE IF (n_sp == 3) THEN      ! thermal sputtering 
                fskm =  (m_gas_i/(2.0_r2 * pi * kB * Temp_gas))**1.5_r2 * (4.0_r2 * pi * v**2.0_r2) * &
                    exp(-m_gas_i/(2.0_r2 * kB * Temp_gas) * v**2.0_r2)
            END IF
            
            !________________________________________________________________________ 
            ! Beginning here: Calculation of the sputtering yield

            !================================================            
            ! 1. Calculate the energy of the gas particle
            E_Y = 0.5_r2 * m_gas_i* v**2.0_r2 ! As m_grain >> m_gas, we assume that the center of the grain is the center of total mass. Energy in Joule.           

             if (n_sp_coulomb == 1) then ! Consider the Coulomb force: 0 = no, 1 = yes
 
                alpha_1 = data_dust_main(i_xgrid, i_ygrid, i_zgrid, 2, i_dusttype, i_asize, 6) * &
                     Ion_grade_i * q_elementar**2.0_r2/(4.0_r2 * pi * epsilon_0 * a(i_asize,i_dusttype)) ! Repulsive: > 0, Attractive: < 0   
                ! In general here should be used not Ion_grade(i), but something like Z_gas_chargr(i). The gas must not be 100% ionised, and also not exact 1.0 ionised per atom (double ionisation etc.).
                ! For partially ionised gas, we should treat the gas as different species (the one with ionised, and the one with unionised gas, even for the same chemical element).

                 if (alpha_1  < E_Y) then  ! Otherwise, the Coulomb-force prevents a collision
                     E_Y = E_Y - alpha_1
                 else     
                     E_Y = 0.0
                 end if
             end if    
            !================================================       
            ! 2. Calculate the threshold energy E_threshold_Y in J
            if ((m_gas_i/ m_grainatom(i_dusttype) ) .le. 0.3_r2) THEN
                g_Y = 4.0_r2 * m_grainatom(i_dusttype) * m_gas_i/(m_grainatom(i_dusttype) + m_gas_i)**2.0_r2  ! no unit 
                E_threshold_Y =  U0(i_dusttype)/(g_Y *(1.0_r2 - g_Y))                                         ! in J
            else    
                E_threshold_Y = 8.0_r2 * U0(i_dusttype) * (m_grainatom(i_dusttype)/m_gas_i)**(-1.0_r2/3.0_r2) ! in J
            end if
            
            !================================================       
            ! 3. Calculate the penetration depth, in nm, Bethe-Bloch-result              
            r_d  = rp(i_dusttype)* (E_Y/ev_unit)*10.0**(2.8*((max(0.1,Ion_grade_i))**(-0.21))+ rp_slope(i_dusttype))                

            !=================================================================================================================================   
            ! 4) We have to distinguish two processes: 4a) Sputtering of the target, 4b) Accretion or trapping of the projectile
            
            ! 4a) What happens with the target? It is sputtered, if the energy is high enough.
            if (E_Y .ge. E_threshold_Y) THEN
        
                ! 4a) Calculate the alpha-function of the sputtering
                ! 1 = Tielens et al. (1994), 2 = Nozawa et al. (2006)
                i_dummy = n_alpha_function
                if (i_dummy == 1) THEN
                    if (m_grainatom(i_dusttype)/m_gas_i > 0.5_r2) THEN
                        alpha_Y =  0.3_r2 * (m_grainatom(i_dusttype)/m_gas_i)**(2.0_r2/3.0_r2)
                    else    
                        alpha_Y = 0.2_r2
                    end if
                else if (i_dummy == 2) THEN
                    if (m_grainatom(i_dusttype)/m_gas_i > 1) THEN
                        alpha_Y =  0.3_r2 * (m_grainatom(i_dusttype)/m_gas_i - 0.6_r2)**(2.0_r2/3.0_r2)
                    else if ((1 >= m_grainatom(i_dusttype)/m_gas_i) .and. (m_grainatom(i_dusttype)/m_gas_i > 0.5_r2) ) THEN
                        alpha_Y =  0.1_r2 * m_gas_i/m_grainatom(i_dusttype) + &
                            0.25_r2 * (m_grainatom(i_dusttype)/m_gas_i - 0.5_r2)**2.0_r2
                    else    
                        alpha_Y = 0.2_r2
                    end if
                else
                    print*, 'Paperboats ERROR 8 in the choice of alpha-function in the sputtering routine.'
                    STOP    
                end if
 
                ! 4b) Calculate the function e12_Y ((4.0 * pi * epsilon_0) is important in SI-units)
                e12_Y = m_grainatom(i_dusttype)/(m_gas_i + m_grainatom(i_dusttype)) * (4.0_r2 * pi * epsilon_0) *&
                    (0.885_r2 * a_Bohr * (Z_gas_i**(2.0_r2/3.0_r2) + Z_grainatom(i_dusttype)**(2.0_r2/3.0_r2))**(-0.5_r2))/ &
                    (Z_gas_i * Z_grainatom(i_dusttype) * q_elementar**2.0_r2) * E_Y ! in J*m/C**2     
            
                ! 4c) Calculate the function sn_Y
                sn_Y = 3.441_r2 * (e12_Y)**0.5 * log(e12_Y + 2.718_r2)/ (1.0_r2 + 6.35_r2 * (e12_Y)**0.5_r2 + e12_Y * &
                    (-1.708_r2 + 6.882_r2 * (e12_Y)**0.5_r2)) 
            
                ! 4d) Calculate the function of the atomic numbers of gas and dust
                Zratio_Y = (Z_gas_i * Z_grainatom(i_dusttype))/ &
                    (Z_gas_i**(2.0_r2/3.0_r2) + Z_grainatom(i_dusttype)**(2.0_r2/3.0_r2))**0.5_r2
                
                !=======================================================================
                ! 5.) Calculate the size-dependent correction factior f_size = Y_a/Y_inf 
                ! (Jurac et al. 1998, Serra Diaz-Cano & Jones 2008, Bocchio et al. 2012, 2014, 2016)
                f_size = 1.0_r2
                
                if (n_sput_size == 1) then
                    ard = a(i_asize,i_dusttype)/(r_d * 1.0E-9) ! Bocchio et al (2012, 2016)
                    
                    f_size = 1.0_r2 + sput_size(i_dusttype, 1) * exp(-((log(ard/sput_size(i_dusttype, 2)))**2.0_r2)/&
                        (2.0_r2*sput_size(i_dusttype, 3)**2.0_r2) ) - sput_size(i_dusttype, 4) &
                        * exp(-(sput_size(i_dusttype, 5) * ard - sput_size(i_dusttype, 6))**2.0_r2)  
                    if (f_size < 0.0_r2) THEN
                        f_size = 0.0_r2 ! The parameter in Bocchio et al (2012) give f_size < 0.0 for a small range of ard (=x), which is not physically
                    end if
                end if

                !================================================                
                ! 6) Finally, calculate the sputtering yield Yield
                Yield = 3.56_r2 * ev_unit /U0(i_dusttype) * (m_gas_i/(m_gas_i + m_grainatom(i_dusttype))) * Zratio_Y  * alpha_Y  * & 
                   (k_value(i_dusttype) * m_grainatom(i_dusttype)/m_gas_i + 1.0_r2)**(-1.0_r2) * sn_Y * (1.0_r2 -    &
                    (E_threshold_Y/E_Y) **(2.0_r2/3.0_r2)) * (1.0_r2 - (E_threshold_Y/E_Y))**2.0_r2 * f_size
                Yield = 2.0 * Yield ! This factor takes care of nonnormal incedent of the gas molecules    
                Yield_rg = 0.0_r2      
                !================================================
    
            else
                !=====================================================
                ! If the enrgy is not high enough, the dust is not sputtered 
                Yield = 0.0_r2 
                Yield_rg = 0.0_r2
            end if 
            
            !=================================================================================================================================               
            ! 4b) What happens with the projectile? It is trapped if the penetration depth is larger than num_lay layers, and smaller than 4/3 of the grain radius.
            If ((1.18416* real(num_lay)* 1.0E+9 * (m_grainatom(i_dusttype)/rho_bulk(i_dusttype))**0.33333 < &  ! real(num_lay) layers ! Trapping!
                (r_d/rp(i_dusttype))) .and. (r_d/rp(i_dusttype) < 4.0/3.0*a(i_asize,i_dusttype)* 1.0E+9)) THEN ! The penetration depth has to be below 4/3 of the grain radius, which is the mean

                IF (i - n_gas_elements == i_dusttype) THEN ! Dusty gas only                    
                    IF (n_gas_accretion == 2) THEN
                        Yield = Yield - 1.0_r2
                        Yield_rg = 0.0_r2
                    END IF                
                ELSE                                       ! Regular gas only
                    IF (n_gas_accretion == 4) THEN
                        Yield    =  Yield -(m_gas_i /m_grainatom(i_dusttype))
                        Yield_rg =  -(m_gas_i /m_grainatom(i_dusttype))     ! Required to account for the trapped gas mass  
                    END IF 
                END IF                
            ELSE ! No trapping  but accretion
!                !===================================================== Ref.: Kirchschlager et al. (2019)
!                 IF  (E_Y < E_threshold_Y)) THEN   ! No trapping and no sputtering, but accretion
!                     IF (i - n_gas_elements == i_dusttype) THEN ! Dusty gas only     
!                         IF ((n_gas_accretion == 1) .or. (n_gas_accretion == 2)) THEN
!                             Yield    = - k_gas_stick * (1.0_r2 - E_Y/E_threshold_Y)         
!                             Yield_rg = 0.0_r2
!                         END IF                
!                     ELSE                                       ! Regular gas only
!                         IF ((n_gas_accretion == 3) .or. (n_gas_accretion == 4)) THEN
!                             Yield    =  - k_gas_stick * (1.0_r2 - E_Y/E_threshold_Y) * (m_gas_i /m_grainatom(i_dusttype)) ! Gas accretion  of the regular gas, Ref.: Kirchschlager et al. (2019)
!                             Yield_rg =  - k_gas_stick * (1.0_r2 - E_Y/E_threshold_Y) * (m_gas_i /m_grainatom(i_dusttype))      ! Required to account for the accreted gas mass, Ref.: Kirchschlager et al. (2019)  
!                         END IF                 
!                     END IF
!                 !ELSE no accretion and no trapping    
!                 END IF
!                !=====================================================   
                
                !=====================================================  Ref.: Kirchschlager et al. (2022)
                E_abs = 1.45 * ev_unit     ! in J, Molpeceres et al. (2019, MNRAS, 482, 5389)
                IF  (E_Y < E_abs) THEN     
                    k_gas_stick  = 1.0
                ELSE
                    k_gas_stick  = (E_abs/E_Y)**2.0
                END IF    
                
                IF (i - n_gas_elements == i_dusttype) THEN ! Dusty gas only     
                    IF ((n_gas_accretion == 1) .or. (n_gas_accretion == 2)) THEN
                        Yield    = Yield - k_gas_stick * 1.0   
                        Yield_rg = 0.0_r2
                    END IF
                ELSE                                       ! Regular gas only
                    IF ((n_gas_accretion == 3) .or. (n_gas_accretion == 4)) THEN
                        Yield    = Yield - k_gas_stick * (m_gas_i /m_grainatom(i_dusttype))  ! Gas accretion  of the regular gas
                        Yield_rg =       - k_gas_stick * (m_gas_i /m_grainatom(i_dusttype))  ! Required to account for the accreted gas mass 
                    END IF 
                END IF               
                !=====================================================     
            END IF
            !=================================================================================================================================                     
                                     
            !===================================================             
            ! 7. Average over the skewed Maxwellian distribution
            Yield_avg(i)    = Yield_avg(i)    + Yield    * fskm * 0.5_r2 * (v + v_before) * Delta_v  
            Yield_avg_rg(i) = Yield_avg_rg(i) + Yield_rg * fskm * 0.5_r2 * (v + v_before) * Delta_v ! Gas accretion  of the regular gas
            !=========================n_gg=========================          
            ! 8. Prepare velocity for the next Calculation step
            v_before = v
            
            ! Ending here: Calculation of the sputtering yield
            !________________________________________________________________________ 
        END Do ! Do Loop j
        
        
        IF (i - n_gas_elements == i_dusttype) THEN ! Dusty gas only
            IF ((n_gas_accretion == 2) .and. (Yield_avg(i) * n_sput < -1.0E-35)) THEN            
                tau_Depl = pi *  (a(i_asize, i_dusttype))**2.0 *  n_sput *  (-Yield_avg(i))**(-1.0) ! Timescale for a single dusty gas atom to be trapped by the dust.
                t_ratio = Delta_t/ tau_Depl
                
                if (t_ratio<0.001) THEN
                    Depl(i) = 1.0_r2 - t_ratio
                elseif (t_ratio>20.0) THEN
                    Depl(i) = 1.0_r2/t_ratio
                else                  
                    Depl(i) = 1.0_r2/t_ratio * (1.0_r2 - exp(-t_ratio)) 
                end if
              
            else 
                Depl(i) = 1.0_r2  
            end if  
        ELSE                                       ! Regular gas only
            IF ((n_gas_accretion == 4) .and. (Yield_avg_rg(i) * n_sput < -1.0E-35)) THEN
                !The dust absorbs gas. To take into account the depletion of the gas in the cell, calculate a Depletion factor:

                tau_Depl = (pi *  (a(i_asize, i_dusttype))**2.0 *  n_sput *  (-Yield_avg_rg(i)))**(-1.0) ! Timescale for a single gas atom to be trapped by the dust.
                t_ratio = Delta_t/ tau_Depl
                
                if (t_ratio<0.001) THEN
                    Depl(i) = 1.0_r2 - t_ratio
                elseif (t_ratio>20.0) THEN
                    Depl(i) = 1.0_r2/t_ratio
                else                  
                    Depl(i) = 1.0_r2/t_ratio * (1.0_r2 - exp(-t_ratio)) 
                end if
                
!               OLD:  tau_Depl = (pi *  (a(i_asize, i_dusttype))**2.0 *  n_sput *  (-Yield_avg_rg(i)))**(-1.0) ! Timescale for a single gas atom to be trapped by the dust.
!               OLD:  Depl(i) = 0.5_r2 * (1.0_r2 - exp(-Delta_t/tau_Depl)) * mu_chem * m_amu/m_grainatom(i_dusttype)  * tau_Depl/Delta_t                    
            ELSE
                Depl(i) = 1.0_r2
            END IF 
        END IF    
        
        IF (Depl(i)> 1.0) THEN
            print*, "Paperboats ERROR 9: Depletion factor is larger than 1!", &
            Depl(i), i_asize, n_sput, Yield_avg_rg(i), tau_Depl/Delta_t
        END IF
    
    END DO     ! Do Loop gas_element i
        
    dadt    = m_grainatom(i_dusttype)/(4.0_r2*rho_bulk(i_dusttype)) * sum(n(:)* Depl(:)*Yield_avg(:))    ! Eq. in Tielens et al. (1994), e.g. Eq. 4.18
    dadt_rg = m_grainatom(i_dusttype)/(4.0_r2*rho_bulk(i_dusttype)) * sum(n(:)* Depl(:)*Yield_avg_rg(:)) ! Same as dadt, but the regular gas that accretes.
    ! Note: The factor of 2 (nonnormal incident of the gas) is removed, so the factor is 1/4, not 1/2
    
    
Deallocate(n, Yield_avg, Yield_avg_rg, Depl)

end subroutine dadt_sputtering ! #s2.3.1

!################################################################################################################
!################################################################################################################

subroutine dust_spreading(i_xgrid, i_ygrid, i_zgrid)  ! #s2.4
  use omp_lib
  use datatype
  use variable

  implicit none
  
  INTEGER, intent(in)  :: i_xgrid, i_ygrid, i_zgrid                    ! input  
  REAL(kind=r2), DIMENSION(:), allocatable :: r_old, r_new, SubVolume  
  REAL(kind=r2) :: Dx, Dy, Dz
  REAL(kind=r2) :: N_old, vx_old, vy_old, vz_old
  REAL(kind=r2) :: N_move, vx_move, vy_move, vz_move, v_press, v_press2, length
  REAL(kind=r2) :: r_out, r_in, wider , factor, Extra_volume, aaa
  REAL(kind=r2) :: V_lo, V_lu, V_ro, V_ru  
  INTEGER :: Dx_int, Dy_int, Dz_int, i, i_new, j_new, k_new, i_timeinter, i_dusttype, i_asize

  Do i_dusttype = 1, n_dusttype
    Do i_asize = 0, n_asize  
    
        ALLOCATE(r_old(1:3), r_new(1:3), SubVolume(1:8))
        
        !=====================================================
        !Choose the center of the current cube
        r_old(1) = data_grid(i_xgrid, i_ygrid, i_zgrid,1) !-unit_length * length_x * (n_xgrid*2.0)**(-1.0)
        r_old(2) = data_grid(i_xgrid, i_ygrid, i_zgrid,2) !-unit_length * length_y * (n_ygrid*2.0)**(-1.0)
        r_old(3) = data_grid(i_xgrid, i_ygrid, i_zgrid,3) !-unit_length * length_z * (n_zgrid*2.0)**(-1.0)

        IF (i_asize > 0) THEN
             IF (l_MHD == 0) THEN
                r_new = r_old + 0.5_r2*Delta_t_inter*(v_dust(:,i_xgrid,i_ygrid,i_zgrid,0, i_dusttype, i_asize) + &
                    v_dust(:,i_xgrid,i_ygrid,i_zgrid,n_timeinter, i_dusttype, i_asize))
                Do i_timeinter = 1,n_timeinter-1
                    r_new = r_new + v_dust(:,i_xgrid,i_ygrid,i_zgrid, i_timeinter, i_dusttype, i_asize) * Delta_t_inter
                END DO
             ELSE
                r_new = r_old            
                Do i_timeinter = 1,n_timeinter
                    r_new = r_new + v_dust_MHD_tra(:,i_xgrid,i_ygrid,i_zgrid, i_timeinter, i_dusttype, i_asize) * Delta_t_inter
                ! v_dust_MHD_act and v_dust_MHD_tra are caclulated for each inter-timestep i_timeinter = 1, n_timeinter on the basis of v_dust(i_timeinter). 
                ! Therfore, the value at i_timeinter=0 is not important.                    
                END DO
             END IF
        ELSE  ! i_asize = 0
            r_new(1) = r_old(1) + data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4) * Delta_t_inter
            r_new(2) = r_old(2) + data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5) * Delta_t_inter
            r_new(3) = r_old(3) + data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6) * Delta_t_inter    
        END IF

        !=====================================================  
        ! Passiv advection
        ! r_new(1) = r_old(1) + data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4) * Delta_t_inter
        ! r_new(2) = r_old(2) + data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5) * Delta_t_inter
        ! r_new(3) = r_old(3) + data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6) * Delta_t_inter 
        !=====================================================    
        
        IF (n_zgrid == 1) THEN 
            r_new(3) = r_old(3)
        END IF         
        
        ! Find i_new, j_new, k_new, use (r_new-r_old)
        Dx = (r_new(1)- r_old(1)) * (n_xgrid*1.0_r2) /( unit_length * length_x)    
        Dy = (r_new(2)- r_old(2)) * (n_ygrid*1.0_r2) /( unit_length * length_y)    
        Dz = (r_new(3)- r_old(3)) * (n_zgrid*1.0_r2) /( unit_length * length_z)   

        Dx_int = int(Dx)   ! int(0.9) = 0, int(1.0)=1, int(-0.1) = 0, int(-0.9)=0, int(-1.0) = -1
        Dy_int = int(Dy)   
        Dz_int = int(Dz)
      
        SubVolume(1) =  (1.0_r2-abs(Dx-Real(Dx_int) ))  * (1.0_r2-abs(Dy-Real(Dy_int) ))  * (1.0_r2-abs(Dz-Real(Dz_int) )) ! l * u
        SubVolume(2) =          abs(Dx-Real(Dx_int) )   * (1.0_r2-abs(Dy-Real(Dy_int) ))  * (1.0_r2-abs(Dz-Real(Dz_int) )) ! r * u 
        SubVolume(3) =          abs(Dx-Real(Dx_int) )   *         abs(Dy-Real(Dy_int) )   * (1.0_r2-abs(Dz-Real(Dz_int) )) ! r * o 
        SubVolume(4) =  (1.0_r2-abs(Dx-Real(Dx_int) ))  *         abs(Dy-Real(Dy_int) )   * (1.0_r2-abs(Dz-Real(Dz_int) )) ! l * o
        SubVolume(5) =  (1.0_r2-abs(Dx-Real(Dx_int) ))  * (1.0_r2-abs(Dy-Real(Dy_int) ))  * abs((Dz-Real(Dz_int) ))
        SubVolume(6) =          abs(Dx-Real(Dx_int) )   * (1.0_r2-abs(Dy-Real(Dy_int) ))  * abs((Dz-Real(Dz_int) ))
        SubVolume(7) =          abs(Dx-Real(Dx_int) )   *         abs(Dy-Real(Dy_int) )   * abs((Dz-Real(Dz_int) ))
        SubVolume(8) =  (1.0_r2-abs(Dx-Real(Dx_int) ))  *         abs(Dy-Real(Dy_int) )   * abs((Dz-Real(Dz_int) ))        
                
         IF (l_scenario == 2) THEN  
            IF ( (i_xgrid .ne. origin_x) .and. (i_ygrid .ne. origin_y)) THEN 
              r_in = ((r_old(1)* (n_xgrid*1.0_r2) /( unit_length * length_x)  -origin_x)**2.0 + &
                    & (r_old(2)* (n_ygrid*1.0_r2) /( unit_length * length_y)  -origin_y)**2.0 + &
                    & (r_old(3)* (n_zgrid*1.0_r2) /( unit_length * length_z)  -origin_z)**2.0 )**0.5
              r_out = ((r_new(1)* (n_xgrid*1.0_r2) /( unit_length * length_x)  -origin_x)**2.0 + &
                     & (r_new(2)* (n_ygrid*1.0_r2) /( unit_length * length_y)  -origin_y)**2.0 + &
                     & (r_new(3)* (n_zgrid*1.0_r2) /( unit_length * length_z)  -origin_z)**2.0 )**0.5  ! in pixel 
                  
              ! In 2D: wider = (r_out/r_in) is the strechting of the transported cell due to the expansion of the material in the SN explosion.
              ! In 3D: wider**2 = (r_out/r_in)**2
              ! Eleminate the cases where divided by "zero" (in the origin)  or where  the material is compressed (r_out<r_in). wider must be >= 1 to be realistic.
              IF (r_in>0) THEN    
                  wider = max(1.0, real(r_out)/real(r_in)) 
              ELSE 
                  wider = max(1.0, real(r_out))
              END IF
  
          
              SubVolume(:) = SubVolume(:)/wider        ! dilution due to the stretching.

              !==============================================================        
              ! Version from 3.12.2020            
              ! All 4 corner (in 2D)
                     If ((1.0 - SubVolume(1))>0) THEN
                         factor = (1.0 - SubVolume(1))/(SubVolume(2) + SubVolume(3) + SubVolume(4))
                         If (factor > 0) THEN
                             SubVolume(2) = SubVolume(2) * factor  !SubVolume(2) + (wider -1.0) * (1.0_r2-abs(Dy-Real(Dy_int) )) 
                             SubVolume(3) = SubVolume(3) * factor  !SubVolume(3) + (wider -1.0) * abs(Dy-Real(Dy_int) )  
                             SubVolume(4) = SubVolume(4) * factor !SubVolume(1) + (wider -1.0) * (1.0_r2-abs(Dy-Real(Dy_int) ))                          
                         END IF    
                     END IF
                 END IF
              !==============================================================            
          END IF   

        IF (l_scenario == 2) THEN
            IF ( (i_xgrid == origin_x) .or. (i_ygrid == origin_y)) THEN
                r_in = ((r_old(1)* (n_xgrid*1.0_r2) /( unit_length * length_x)  -origin_x)**2.0 + &
                   & (r_old(2)* (n_ygrid*1.0_r2) /( unit_length * length_y)  -origin_y)**2.0 + &
                   & (r_old(3)* (n_zgrid*1.0_r2) /( unit_length * length_z)  -origin_z)**2.0 )**0.5
                r_out = ((r_new(1)* (n_xgrid*1.0_r2) /( unit_length * length_x)  -origin_x)**2.0 + &
                    & (r_new(2)* (n_ygrid*1.0_r2) /( unit_length * length_y)  -origin_y)**2.0 + &
                    & (r_new(3)* (n_zgrid*1.0_r2) /( unit_length * length_z)  -origin_z)**2.0 )**0.5  ! in pixel 
                 
                ! In 2D: wider = (r_out/r_in) is the strechting of the transported cell due to the expansion of the material in the SN explosion.
                ! In 3D: wider**2 = (r_out/r_in)**2
                ! Eleminate the cases where divided by "zero" (in the origin)  or where  the material is compressed (r_out<r_in). wider must be >= 1 to be realistic.
                IF (r_in>0) THEN    
                    factor = min(1.0, real(r_in)/real(r_out)) 
                ELSE 
                    factor = 1.0
                END IF
                
                SubVolume(:) = SubVolume(:)* factor       ! dilution due to the stretching and expansion.
                
                ! (1.0 - factor) is the missing volume which has to be expanded to the neighbouring cells.
                  IF ((i_xgrid == origin_x) .and. (i_ygrid > origin_y)) THEN
                     V_lu = (1.0 - factor)*0.5 * (1.0_r2-abs(Dy-Real(Dy_int) ))
                     V_ru = (1.0 - factor)*0.5 * (1.0_r2-abs(Dy-Real(Dy_int) ))
                     V_ro = (1.0 - factor)*0.5 * (       abs(Dy-Real(Dy_int) ))
                     V_lo = (1.0 - factor)*0.5 * (       abs(Dy-Real(Dy_int) ))                    
                     
                     SubVolume(2) =  SubVolume(2) + V_ru
                     SubVolume(3) =  SubVolume(3) + V_ro
                     
                  ELSE IF ((i_xgrid == origin_x) .and. (i_ygrid < origin_y)) THEN               
                     V_lo = (1.0 - factor)*0.5 * (1.0_r2-abs(Dy-Real(Dy_int) ))   
                     V_ro = (1.0 - factor)*0.5 * (1.0_r2-abs(Dy-Real(Dy_int) )) 
                     V_lu = (1.0 - factor)*0.5 * (       abs(Dy-Real(Dy_int) ))                    
                     V_ru = (1.0 - factor)*0.5 * (       abs(Dy-Real(Dy_int) ))                     
                     
                     SubVolume(2) =  SubVolume(2) + V_ro
                     SubVolume(3) =  SubVolume(3) + V_ru
                    
                 ELSE IF ((i_ygrid == origin_y) .and. (i_xgrid > origin_x)) THEN
                    V_lu = (1.0 - factor)*0.5 * (1.0_r2-abs(Dx-Real(Dx_int) ))                 
                    V_ru = (1.0 - factor)*0.5 * (       abs(Dx-Real(Dx_int) ))                 
                    V_lo = (1.0 - factor)*0.5 * (1.0_r2-abs(Dx-Real(Dx_int) ))                 
                    V_ro = (1.0 - factor)*0.5 * (       abs(Dx-Real(Dx_int) )) 
                    
                    SubVolume(3) =  SubVolume(3) + V_ro
                    SubVolume(4) =  SubVolume(4) + V_lo                    
                    
                 ELSE IF ((i_ygrid == origin_y) .and. (i_xgrid < origin_x)) THEN
                    V_lu = (1.0 - factor)*0.5 * (       abs(Dx-Real(Dx_int) ))                 
                    V_ru = (1.0 - factor)*0.5 * (1.0_r2-abs(Dx-Real(Dx_int) ))                 
                    V_lo = (1.0 - factor)*0.5 * (       abs(Dx-Real(Dx_int) ))                 
                    V_ro = (1.0 - factor)*0.5 * (1.0_r2-abs(Dx-Real(Dx_int) ))  
                    
                    SubVolume(3) =  SubVolume(3) + V_lo
                    SubVolume(4) =  SubVolume(4) + V_ro                        
                  END IF
            
            END IF
        END IF
        !========================================================================
        Do i = 1, 8 
            
            i_new = i_xgrid + Dx_int 
            j_new = i_ygrid + Dy_int
            k_new = i_zgrid + Dz_int
            
            IF ((i==2) .or. (i==3) .or. (i==6) .or. (i==7)) THEN
                i_new = i_new + nint(sign(1.0_r2, Dx)) 
            END IF
            
            IF ((i==3) .or. (i==4) .or. (i==7) .or. (i==8)) THEN           
                j_new = j_new + nint(sign(1.0_r2, Dy)) 
            END IF
            
            IF (( i > 4)) THEN 
                 k_new = k_new + nint(sign(1.0_r2, Dz)) 
            END IF
            
            !========================================================================   
            !========================================================================               
            ! Periodic boundaries

            IF (l_boundary == 1) THEN 
                    i_new = i_new -  floor(real(i_new-1)/real(n_xgrid))*  n_xgrid
                    j_new = j_new -  floor(real(j_new-1)/real(n_ygrid))*  n_ygrid
                    k_new = k_new -  floor(real(k_new-1)/real(n_zgrid))*  n_zgrid                     
            END IF           

            If ( (0 < i_new) .and. (i_new < n_xgrid+1) .and.  (0 < j_new) .and. &
                (j_new < n_ygrid+1) .and. (0 < k_new) .and. (k_new < n_zgrid+1)) THEN

                !============================================================================== 
                !$omp CRITICAL
                N_old =  data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize, 4)      
                vx_old = data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize, 1)     
                vy_old = data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize, 2)    
                vz_old = data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize, 3)        
            
                N_move  = data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, 4) * SubVolume(i)
!                 vx_move = data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, 1)      
!                 vy_move = data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, 2)    
!                 vz_move = data_dust_next(i_xgrid, i_ygrid, i_zgrid, i_dusttype, i_asize, 3)                    
            
                IF (i_asize == 0) THEN
                    vx_move = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4)
                    vy_move = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5)
                    vz_move = data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6)
                ELSE
                    IF (l_MHD == 0) THEN 
                        vx_move =  v_dust(1,i_xgrid,i_ygrid,i_zgrid,n_timeinter, i_dusttype, i_asize)                   
                        vy_move =  v_dust(2,i_xgrid,i_ygrid,i_zgrid,n_timeinter, i_dusttype, i_asize)                   
                        vz_move =  v_dust(3,i_xgrid,i_ygrid,i_zgrid,n_timeinter, i_dusttype, i_asize)                    
                    ELSE
                        vx_move = v_dust_MHD_tra(1,i_xgrid,i_ygrid,i_zgrid,n_timeinter, i_dusttype, i_asize)
                        vy_move = v_dust_MHD_tra(2,i_xgrid,i_ygrid,i_zgrid,n_timeinter, i_dusttype, i_asize)
                        vz_move = v_dust_MHD_tra(3,i_xgrid,i_ygrid,i_zgrid,n_timeinter, i_dusttype, i_asize)
                    END IF
                END IF
            
                !=======>
                If ((N_move > 0.0_r2)) THEN      
                    data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize,1) =  &
                    (N_old*vx_old + N_move* vx_move ) /(N_old + N_move) 
                    data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize,2) =  &
                    (N_old*vy_old + N_move* vy_move ) /(N_old + N_move)      
                    data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize,3) =  &
                    (N_old*vz_old + N_move* vz_move ) /(N_old + N_move) 
                    data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize,4) = (N_old + N_move)   
                ELSE
                    data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize,1) = vx_old 
                    data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize,2) = vy_old
                    data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize,3) = vz_old          
                    data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize,4) = N_old
                END IF
            
                IF (i_asize > 0) THEN
                    data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize, 5) =   &
                        data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize, 4) * &
                            4.0_r2/3.0_r2 * pi * (a(i_asize,  i_dusttype))**3.0_r2 * rho_bulk(i_dusttype)
                ! Take here the data of timepoint '3', because it is the mass of the number
                ! of dust particles at  timepoint '3'.
                ELSE ! i_asize = 0
                    data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize, 5) = &
                        data_dust_main(i_new, j_new, k_new, 3, i_dusttype, i_asize, 4) * m_grainatom(i_dusttype)
                END IF    
                !$omp end CRITICAL
            END IF
        END DO  ! end i = 1, 8

          !========================================================================         

        DEALLOCATE(r_old, r_new, SubVolume)

    END DO ! END i_dusttype = 1, n_dustty
  END DO   ! END i_asize = 0, n_asize  
  
end subroutine dust_spreading ! #s2.4

!################################################################################################################
!################################################################################################################

subroutine charge_calculation(i_xgrid, i_ygrid, i_zgrid, Temp_gas) ! #s2.5
    use omp_lib
    use datatype
    use variable

    implicit none
    
    INTEGER, intent(in)  :: i_xgrid, i_ygrid, i_zgrid  ! input      
    REAL (kind=r2), intent(in)  :: Temp_gas 
    
    INTEGER :: i_dusttype, i_asize 
    
    REAL(kind=r2) :: Z_grain, v_cgs, a_cgs
    REAL(kind=r2) :: phi_imp, phi_stat, phi_see1, phi_see2, phi_tran, phi_therm, phi_tot
    REAL(kind=r2) :: psi_A, psi_B, psi_C, psi_D, psi_E, psi_F
    REAL(kind=r2) :: lambda_esc, T_rel, T_imp, T_tran, T_cr, T_see1, T_see2, T_therm
    REAL(kind=r2) :: Theta_A, Theta_B, Theta_C
    REAL(kind=r2) :: phi_min, phi_max
    
    
    Do i_dusttype = 1, n_dusttype
       Do i_asize = 1, n_asize      
       
       Z_grain = 0.0_r2 ! if n)plasma = 0, Z_grain is always zero.
    !======================================================================================           
    ! Charge of a grain in rest, simple equation (Spitzer 1941) 
   
    If ((n_plasma == 1)) THEN
        Z_grain = -2.504_r2 * a(i_asize,i_dusttype) * (kB * Temp_gas) *&
            (4.0_r2 * pi * epsilon_0) / q_elementar**2.0_r2     !
    END IF
    
    !======================================================================================           
    ! Charge of a moving grain (relative to gas), plus secondary electron emission,  (see Fry et al. 2018) 
    If ((n_plasma == 2)) THEN
        ! velocity of the dust relative to the gas, in cgs. We want to consider the dust 
        ! is (by the gas at timestep i_timepoints) sputtered, gg-collided, rearranged. Thatswhy
        ! we choose the dust velocity at i_timepoints+1, but the gas-velocity is still at i_timepoint
        v_cgs = (((data_dust_main(i_xgrid, i_ygrid, i_zgrid, 3,i_dusttype, i_asize, 1) -&
            data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,4))**2.0_r2 &
            + (data_dust_main(i_xgrid, i_ygrid, i_zgrid, 3,i_dusttype, i_asize, 2) - &
            data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,5))**2.0_r2 &
            + (data_dust_main(i_xgrid, i_ygrid, i_zgrid, 3,i_dusttype, i_asize, 3) - &
            data_gas(i_xgrid, i_ygrid, i_zgrid, i_timepoints,6))**2.0_r2)**0.5_r2) *100.0_r2 ! in cm/s
            
        a_cgs = a(i_asize,i_dusttype)* 100.0  
    
        lambda_esc = 7.1477_r2 * (0.4_r2/1.1611_r2)**1.5935_r2 ! for Fe2O3, Fry et al. (2018), Eq. E10. (low impact on the grain charge)
        T_rel   = 0.2506_r2 * (v_cgs * 1.0D-7)**2.0_r2
        T_imp   = 0.3433_r2 * (v_cgs * 1.0E-7)**2.0_r2
        T_tran  = 10.57_r2/(1.0_r2 - exp(-(lambda_esc/a_cgs)**0.75_r2))
        T_see1  = 3.404_r2 * (v_cgs * 1.0E-7)**2.0_r2 
        T_see2  = 34.82_r2 * (v_cgs * 1.0E-7)**1.223_r2
        T_cr    = max(T_tran, T_imp)        
        T_therm = max(703.8_r2, 9.964_r2 * (v_cgs * 1.0E-7)**2.0_r2)
        
        Theta_A = 36.99_r2
        Theta_B = 38.48_r2
        Theta_C = 0.3545_r2 * log(v_cgs * 1.0E-7) + 1.563_r2
        ! if v_cgs = 0, Theta_C is infinty. However, psi_C is then 0.0, and this doesn't matter.
        
        phi_imp   = -0.084_r2 + 1.112E-3 * (v_cgs * 1.0E-7)**2.0_r2 + (T_rel/(Temp_gas * 1.0E-5))**0.75_r2   
        phi_stat  = phi_stat_help ! Calculated in #s1.3.1
        phi_see1  =  2.745_r2  - 1.740_r2 * exp(-0.1037_r2 * (v_cgs*1.0E-7)**2.0_r2)
        phi_see2  = max(0.0_r2, -0.2267_r2 * (v_cgs * 1.0E-7)**2.0_r2 + 1.430_r2)
        phi_tran  =  0.1953_r2 * (Temp_gas * 1.0E-5)**(-0.162_r2)
        phi_therm =  0.1862_r2 * log(Temp_gas * 1.0E-5) - 1.756_r2

        psi_A = ((1.0_r2 + (T_imp /(Temp_gas * 1.0E-5))**Theta_A)**(-1.0_r2)) * max(0.0_r2,(T_tran - T_imp)/abs(T_tran - T_imp))        
        psi_B = ((1.0_r2 + (T_see1/(Temp_gas * 1.0E-5))**Theta_B)**(-1.0_r2))
        psi_C = ((1.0_r2 + (T_see2/(Temp_gas * 1.0E-5))**Theta_C)**(-1.0_r2))
        psi_D = exp(-(T_cr/(Temp_gas * 1.0E-5))**4.0_r2)
        psi_E = exp(-1.0E-4 * (a_cgs/lambda_esc)**4.0_r2)
        psi_F = exp(-((Temp_gas * 1.0E-5)/T_therm)**4.0_r2)
        
        ! Only as a note: phi_tot and phi_pl in the plasma-drag routine are the same quantity. 
        phi_tot =  phi_imp   *                 (1.0_r2 - psi_B) * psi_A * (1.0_r2 - psi_D) * psi_F + &
                   phi_stat  *                                 psi_A * (1.0_r2 - psi_D) * psi_F + &
                  (phi_see1 + phi_see2 * psi_C)     * psi_B  * psi_A * (1.0_r2 - psi_D) * psi_F + &
                   phi_tran  *           psi_E      *                         psi_D          + &   
                   phi_therm *                        psi_B  *         (1.0_r2 - psi_D) * (1.0_r2 - psi_F)
        
        ! Lastly, we also establish potential limits to account for field emission (McKee et al. 1987):
        phi_min = -11.6_r2 * (a_cgs * 1.0E+5) * (Temp_gas * 1.0E-5)**(-1.0_r2)
        phi_max = 348.0_r2 * (a_cgs * 1.0E+5) * (Temp_gas * 1.0E-5)**(-1.0_r2)
        
         if (phi_tot < phi_min) then
             phi_tot = phi_min
         elseif (phi_tot > phi_max) then
             phi_tot = phi_max             
         end if
     
        Z_grain = phi_tot * a(i_asize,i_dusttype) * (kB * Temp_gas) *&
            (4.0_r2 * pi * epsilon_0) / q_elementar**2.0_r2
    end if
        
    data_dust_main(i_xgrid, i_ygrid, i_zgrid, 3, i_dusttype, i_asize, 6) = Z_grain
    
        END DO
    END DO
    
end subroutine charge_calculation ! #s2.5

!################################################################################################################
!################################################################################################################

subroutine lift ! #s2.6
    use omp_lib
    use datatype
    use variable
    
    INTEGER :: i_xgrid, i_ygrid, i_dusttype, i_asize
    REAL(kind=r2), DIMENSION(:), allocatable :: m_dummy(:), m_dummx(:)
    
 
      data_dust_main(:, :, :, 2, :, :, :)  = data_dust_main(:, :, :, 3, :, :, :) 
      data_dust_main(:, :, :, 3, :, :, :)  = 0.0_r2
      
        !========================================================         
        ! Correct the symmetry effect:
        
          IF ((l_scenario == 2) .and. (l_hydrocode == 2)) THEN      ! FK, 01/12/2020
            ALLOCATE(m_dummx(1:n_xgrid), m_dummy(1:n_ygrid))
            Do i_dusttype = 1, n_dusttype
              Do  i_asize = 1, n_asize
                  !===============================================================================================               
                  Do i_ygrid = 1, n_ygrid
                        m_dummx(i_ygrid) = (&
                        data_dust_main(nint(origin_x)-2, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) + &
                        data_dust_main(nint(origin_x)-1, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) + &
                        data_dust_main(nint(origin_x)  , i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) + &
                        data_dust_main(nint(origin_x)+1, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) +&
                        data_dust_main(nint(origin_x)+2, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) &
                        )/5.0
  
                       
                       data_dust_main(nint(origin_x)-2, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) = m_dummx(i_ygrid)
                       data_dust_main(nint(origin_x)-1, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) = m_dummx(i_ygrid)
                       data_dust_main(nint(origin_x)  , i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) = m_dummx(i_ygrid)
                       data_dust_main(nint(origin_x)+1, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) = m_dummx(i_ygrid)
                       data_dust_main(nint(origin_x)+2, i_ygrid, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) = m_dummx(i_ygrid)                  
                  END DO
                  !===============================================================================================
  
                  !===============================================================================================               
                  Do i_xgrid = 1, n_xgrid
                     m_dummy(i_xgrid) = (&
                          data_dust_main(i_xgrid, nint(origin_y)-2, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) + &
                          data_dust_main(i_xgrid, nint(origin_y)-1, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) + &
                          data_dust_main(i_xgrid, nint(origin_y)  , (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) + &
                          data_dust_main(i_xgrid, nint(origin_y)+1, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) +&
                          data_dust_main(i_xgrid, nint(origin_y)+2, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) &
                        )/5.0
  
                       
                       data_dust_main(i_xgrid,nint(origin_y)-2, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) = m_dummy(i_xgrid)
                       data_dust_main(i_xgrid,nint(origin_y)-1, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) = m_dummy(i_xgrid)
                       data_dust_main(i_xgrid,nint(origin_y)  , (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) = m_dummy(i_xgrid)
                       data_dust_main(i_xgrid,nint(origin_y)+1, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) = m_dummy(i_xgrid)
                       data_dust_main(i_xgrid,nint(origin_y)+2, (n_zgrid+1)/2, 2, i_dusttype,i_asize, 4) = m_dummy(i_xgrid)   
                  END DO                      
                  !===============================================================================================     

               END DO
             END DO
            DEALLOCATE(m_dummx, m_dummy)    
          END IF   

end subroutine lift    

end module processing
