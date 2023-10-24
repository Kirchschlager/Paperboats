PROGRAM Paperboats
use omp_lib
use variable
use datatype
use init
use processing
use output_plot
use clean

implicit none

REAL    :: t_cpu_0,t_cpu_1, t_cpu_second, t_coreh1, t_coreh2 
INTEGER :: t_cpu_hour, t_cpu_minute

print*, " "
print*, " #######################################################"
print*, "   Welcome to Paperboats!                               "
print*, "                                      ___/\___          "
print*, "                                      \______/          "
print*, "                              ___/\___                  "
print*, "                              \______/   ___/\___       " 
print*, "                                         \______/       "
print*, " "
print*, "   Version 3.0                                          "
print*, " #######################################################"
print*, " "

    t_cpu_0 = secnds(0.0)
    call cpu_time(t_coreh1)

    !--------------------------------------------------------------------------!
    ! 1. initiate code 
    !--------------------------------------------------------------------------!
    call init_paperboats

    !--------------------------------------------------------------------------!
    ! 2. time evolution 
    !--------------------------------------------------------------------------!
 
    call run_paperboats

    !--------------------------------------------------------------------------!
    ! 3. store and plot results 
    !--------------------------------------------------------------------------!
  
    call output_paperboats_plot
    
    !--------------------------------------------------------------------------!
    ! 4. clean up
    !--------------------------------------------------------------------------!

    call cleanup

    !--------------------------------------------------------------------------!
    ! 5. time measurement
    !--------------------------------------------------------------------------!
    
    !==============================================
    ! 5.1. Real time
    t_cpu_1 = secnds(t_cpu_0)
    t_cpu_hour   = int((t_cpu_1)/real(3600))
    t_cpu_minute = int((t_cpu_1)/real(60)) - t_cpu_hour*60
    t_cpu_second = t_cpu_1 - t_cpu_hour * 3600 - t_cpu_minute*60

    print*, ' ______________________________________________________'         
    print*, ''
    print '(a, i4, a,i2, a,f4.1, a, f10.1,a)', "  Done. Simulation took",t_cpu_hour, ' h ', &
    t_cpu_minute, ' min ', t_cpu_second, ' s ( =', t_cpu_1 ,' s).'
    !==============================================
    
    !==============================================
    ! 5.2. Core time    
    call cpu_time(t_coreh2)    
    t_cpu_hour   = int((t_coreh2-t_coreh1)/real(3600))
    t_cpu_minute = int((t_coreh2-t_coreh1)/real(60)) - t_cpu_hour*60
    t_cpu_second = (t_coreh2-t_coreh1) - t_cpu_hour * 3600 - t_cpu_minute*60
      
    print*, ''
    print '(a, i6, a,i2, a,f4.1, a, f10.1,a)', "             Core time:",t_cpu_hour, ' core h ', &
    t_cpu_minute, ' core min ', t_cpu_second, ' core s.'
    !==============================================    
    
    print*, " "
    print*, " Bye bye."
    print*          
          
END PROGRAM
