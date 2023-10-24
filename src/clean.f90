module clean
  use omp_lib
  use datatype
  use variable
  
  implicit none
  !--------------------------------------------------------------------------!
  PUBLIC :: cleanup
  !--------------------------------------------------------------------------!
CONTAINS 

!################################################################################################################

subroutine cleanup ! #s4
  use omp_lib
  use datatype
  use variable
  
  implicit none
    
   DEALLOCATE(v_dust, densi_min_tot, densi_max_tot, ymin, ymax, data_dust_next) 
   DEALLOCATE(a, M_larger_grains, count_colli, c_speed, Pl_Hir, Pv_Hir, s_Hir)
   DEALLOCATE(time, data_grid, data_gas, data_dust_main, k_value, sput_size, Ion_grade, Ion_grade_average)
   DEALLOCATE(Abun_gas, m_gas, Z_gas, Material_name, gamma_asize, sigma_log_n, peak_logn_a, m_grainatom, Z_grainatom, U0)
   DEALLOCATE(v_Thres_vapo, v_Thres_shat, Stress_Thres, E_Possion, gamma_surf, i_min_ini, i_max_ini, rp, rp_slope, DB)
   
   IF (l_MHD == 1) THEN
     DEALLOCATE(v_dust_MHD_tra, v_dust_MHD_act)
   END IF
   
end subroutine cleanup ! #s4

!################################################################################################################

end module clean