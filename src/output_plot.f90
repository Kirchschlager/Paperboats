
      module output_plot
  use omp_lib
  use datatype
  use variable

  implicit none
    !--------------------------------------------------------------------------!
!     PRIVATE
!         CHARACTER(len=4), PARAMETER :: file_a  = "huhu"
    !--------------------------------------------------------------------------!
    PUBLIC :: output_paperboats_plot, output_time_def, output_gasdata_plot,    &
    output_dustdata_plot, output_dustdistributiondata_plot, output_movies
    !--------------------------------------------------------------------------!
CONTAINS

!################################################################################################################

subroutine output_paperboats_plot ! #s3
  use omp_lib
  use datatype
  use variable

  implicit none

  print*, ''
  print*, ''
  print*, ' 3. Save data.'

  call output_gasdata_plot
  call output_dustdata_plot
  call output_dustdistributiondata_plot
  call output_movies

end subroutine output_paperboats_plot ! #s3

!################################################################################################################

subroutine output_time_def ! #s3.1
  use omp_lib
  use datatype
  use variable

  implicit none

  !=================================================
  !1. What is the appropriate unit for the time-output? yr in most cases/

    IF (time(n_timepoints) > 3600.1*24.0*365.25*10000.0) THEN                                                    ! Kiloyears
        time_unit_scaler  = 1.0/(3600.0*24.0*365.25*1000.0)
        time_word = "kyr"
    else if ((time(n_timepoints) >= 3600.1*24.0*365.25) .and. (time(n_timepoints) < 3600.1*24.0*3652500.0)) THEN ! years
            time_unit_scaler  = 1.0/(3600.0*24.0*365.25)
            time_word = "yr"
    else if ((time(n_timepoints) >= 3600.1*24.0) .and. (time(n_timepoints) < 3600.1*24.0*365.25)) THEN           ! days
            time_unit_scaler    = 1.0/(3600.0*24.0)
            time_word = "d"
    else if ((time(n_timepoints) >= 3600.1) .and. (time(n_timepoints) < 3600.1*24.0)) THEN                       ! hours
            time_unit_scaler    = 1.0/3600.0
            time_word = "h"
    else if ((time(n_timepoints) >= 60.01) .and. (time(n_timepoints) < 3600.1)) THEN                             ! minutes
            time_unit_scaler    = 1.0/60.0
            time_word = "min"
    else if (time(n_timepoints) < 60.01) THEN                                                                    ! seconds
            time_unit_scaler   = 1.0
            time_word = "s"
    END IF

end subroutine output_time_def ! #s3.1

!################################################################################################################

subroutine output_gasdata_plot ! #s3.2
  use omp_lib
  use datatype
  use variable
  use tools

  implicit none

  REAL(kind=r2) :: lab_step

   T_min_tot = minval(data_gas(:, :, (n_zgrid+1)/2, :, 8))
   T_max_tot = maxval(data_gas(:, :, (n_zgrid+1)/2, :, 8))

   gas_min_tot = minval(data_gas(:, :, (n_zgrid+1)/2, :, 7))/(mu_chem * m_amu * 1.0E+6)
   gas_max_tot = maxval(data_gas(:, :, (n_zgrid+1)/2, :, 7))/(mu_chem * m_amu * 1.0E+6)

  !=============================================================
  ! Files with parameter to plot the gas density and temperature
  !
  !=============gnu1
  1 format (*(G0,1X))
  open(unit=2, file="./"//Dens_map_folder//"data_gas/gnu.dat", &
        & action="write", status="unknown", form="formatted")
        write(unit=2,fmt=1)  " # Parameter to plot gas density"
        write(unit=2,fmt=1)  " set terminal postscript enhanced color"
        write(unit=2,fmt=1)  " set encoding iso_8859_1"
        write(unit=2,fmt=1)  " ################"
        write(unit=2,fmt=1)  " set size ratio", real(n_ygrid)/real(n_xgrid)
        write(unit=2,fmt=1)  " unset key"
        write(unit=2,fmt=1)  " set origin -0.04,0.01"
        write(unit=2,fmt=1)  " ################"
        write(unit=2,fmt=1)  " set grid front "
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " set xlabel '{/Oblique x} (pc)' font 'Helvetica,22'"
        write(unit=2,fmt=1)  " set ylabel '{/Oblique y} (pc)' offset -1.0,0.0 font 'Helvetica,22'"
        write(unit=2,fmt=1)  " set cblabel '{/Oblique n}_{/=14 gas} (cm^{-3})' offset 3.0,0.0 font 'Helvetica,22'"
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " set xrange [0:",(length_x*unit_length*m_in_pc),"]"
        write(unit=2,fmt=1)  " set yrange [0:",(length_y*unit_length*m_in_pc),"]"
        write(unit=2,fmt=1)  " "

        call label_generator(1,lab_step)

        write(unit=2,fmt=1)  " set xtics ",lab_step," font 'Helvetica,18'"
        write(unit=2,fmt=1)  " set ytics ",lab_step," font 'Helvetica,18'"
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " set mxtics 4"
        write(unit=2,fmt=1)  " set mytics 4"
        write(unit=2,fmt=1)  " "

        IF (gas_max_tot < 10.0* gas_min_tot) THEN
            write(unit=2,fmt=1)  " unset logscale cb "
            write(unit=2,fmt=1)  " set cbrange [",gas_min_tot,":", gas_max_tot,"]"
            write(unit=2,fmt=1)  " set format cb  '%5.3t{/Symbol \327}10^{%L}'"
            write(unit=2,fmt=1)  " set cbtics (",gas_min_tot,",",gas_min_tot/3.0 + gas_max_tot/1.5,",",&
              &  gas_min_tot/1.5 + gas_max_tot/3.0,",", gas_max_tot,") font 'Helvetica,18'"
        ELSE
            !gas_min_tot = max(gas_min_tot,gas_max_tot*0.25E-6)

            write(unit=2,fmt=1)  " set logscale cb "
            write(unit=2,fmt=1)  " set cbrange [",gas_min_tot,":", 0.25*gas_max_tot,"]"
            write(unit=2,fmt=1)  " set format cb '10^{%L}'  "
            write(unit=2,fmt=1)  " set cbtics font 'Helvetica,18'"
            write(unit=2,fmt=1)  " set mcbtics 10"
        END IF
!         write(unit=2,fmt=1)  " set cbrange [*:*]"

        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " unset grid"
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " load '../../../misc/"//theme_colour1//".col'"
  close(unit=2)

  !=============gnu2
  open(unit=2, file="./"//Dens_map_folder//"data_gas/gnu2.dat", &
        & action="write", status="unknown", form="formatted")
        write(unit=2,fmt=1)  " # Parameter to plot gas temperature"
        write(unit=2,fmt=1)  " set cblabel '{/Oblique T} (K)' offset 3.0,0.0 font 'Helvetica,22'"
        write(unit=2,fmt=1)  " set logscale cb "
!         write(unit=2,fmt=1)  " set cbrange [10**3:10**9]"
        write(unit=2,fmt=1)  " set cbrange [",T_min_tot,":",T_max_tot,"]"
        write(unit=2,fmt=1)  " set format cb '10^{%L}'"
        write(unit=2,fmt=1)  " unset cbtics"
        write(unit=2,fmt=1)  " set cbtics font 'Helvetica,18'"
        write(unit=2,fmt=1)  " set mcbtics 10"
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " load '../../../misc/"//theme_colour2//".col'"
  close(unit=2)

  !=============gnu3
  IF  (l_MHD == 1) THEN
    open(unit=2, file="./"//Dens_map_folder//"data_gas/gnu3.dat", &
        & action="write", status="unknown", form="formatted")
        write(unit=2,fmt=1)  " # Parameter to plot magnetic field"
        write(unit=2,fmt=1)  " set cblabel '{/Oblique B} ({/Symbol m}G)' offset 3.0,0.0 font 'Helvetica,22'"
        write(unit=2,fmt=1)  " unset logscale cb"
!         write(unit=2,fmt=1)  " set cbrange [0:450]"
        write(unit=2,fmt=1)  " set cbrange [0:", max(0.5*B_max_tot,0.001),"]"
        write(unit=2,fmt=1)  " unset format cb"
        write(unit=2,fmt=1)  " unset cbtics"
        write(unit=2,fmt=1)  " set cbtics font 'Helvetica,18'"
        write(unit=2,fmt=1)  " set mcbtics 4"
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " load '../../../misc/"//theme_colour3//".col'"
    close(unit=2)
  END IF
  !          write(unit=2,fmt=1)  " splot 'Density_"//i_timepoints_word//".dat' u 1:2:5 w pm3d"
  !
  !=============================================================

  Do i_timepoints = 1,n_timepoints

    write (*,'(A,I3,A)') '      A) Gas flow data:  ', int(i_timepoints/real(n_timepoints)*100.0),&
        & ' % done.' //char(27)//'[A'

    write(i_timepoints_word, '(i5.5)' ) i_timepoints
    write(t_timepoints_word, '(f9.2)' ) (time(i_timepoints)) * (time_unit_scaler)

    !================================================================
    ! Plot the gas number density, temperature and magnetic field as a function of time

    open(unit=2, file="./"//Dens_map_folder//"data_gas/plot_Density_"//i_timepoints_word//".dat", &
        & action="write", status="unknown", form="formatted")
        write(unit=2,fmt=1)  " reset"
        write(unit=2,fmt=1)  " set output '|ps2pdf - ../pics_gas/Density_"//i_timepoints_word//".pdf' "
        write(unit=2,fmt=1)  " set title '{/Oblique t} = "//t_timepoints_word//" "//time_word//"' offset 0,-0.8 font 'Helvetica,18'"
        write(unit=2,fmt=1)  " "

        write(unit=2,fmt=1)  " do for [i=0:",1+l_MHD,"] { "
        write(unit=2,fmt=1)  " if (i == 0)   { load 'gnu.dat'; plot 'Density_"//i_timepoints_word//".dat' u (",&
            &  (length_x*unit_length*m_in_pc)/real(n_xgrid)," * (($2)+0.5)):(",&
            &  (length_y*unit_length*m_in_pc)/real(n_ygrid)," * (($1)+0.5)):3 index i matrix w image } else {"
        write(unit=2,fmt=1)  " if (i == 1)   { load 'gnu2.dat'; plot 'Density_"//i_timepoints_word//".dat' u (",&
            &  (length_x*unit_length*m_in_pc)/real(n_xgrid)," * (($2)+0.5)):(",&
            &  (length_y*unit_length*m_in_pc)/real(n_ygrid)," * (($1)+0.5)):3 index i matrix w image } else {"
        write(unit=2,fmt=1)  " if (i == 2)   { load 'gnu3.dat'; plot 'Density_"//i_timepoints_word//".dat' u (",&
            &  (length_x*unit_length*m_in_pc)/real(n_xgrid)," * (($2)+0.5)):(",&
            &  (length_y*unit_length*m_in_pc)/real(n_ygrid)," * (($1)+0.5)):3 index i matrix w image,\"
        write(unit=2,fmt=1)  "                                 'Mag_field_"//i_timepoints_word//".dat' u (",&
            &  (length_y*unit_length*m_in_pc)/real(n_ygrid)," * (($1)+0.5)):(",&
            &  (length_x*unit_length*m_in_pc)/real(n_xgrid)," * (($2)+0.5)) w l lc rgb 'white' lw 3 }"
        write(unit=2,fmt=1)  " }}}"
    close(unit=2)
  END DO

  print*, ''

end subroutine output_gasdata_plot ! #s3.2

!################################################################################################################

subroutine output_dustdata_plot ! #s3.3
  use omp_lib
  use datatype
  use variable
  use tools

  implicit none
  INTEGER :: i_dusttype, i_asize
  INTEGER :: i_asize_reduced, n_asize_reduced
  REAL(kind=r2) :: lab_step, cb_min, cb_max

  ALLOCATE(a_word(1:n_asize))

  !=============================================================
  ! File with parameters to plot the dust density
  !
  1 format (*(G0,1X))
  open(unit=2, file="./"//Dens_map_folder//"data_dust/gnu.dat", action="write", status="unknown", form="formatted")
        write(unit=2,fmt=1)  " # Parameter to plot dust density"
        write(unit=2,fmt=1)  " set terminal postscript enhanced color"
        write(unit=2,fmt=1)  " set encoding iso_8859_1"
        write(unit=2,fmt=1)  " ################"
        write(unit=2,fmt=1)  " set size ratio", real(n_ygrid)/real(n_xgrid)
        write(unit=2,fmt=1)  " unset key"
        write(unit=2,fmt=1)  " set origin -0.04,0.01"
        write(unit=2,fmt=1)  " ################"
        write(unit=2,fmt=1)  " set grid front "
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " set xlabel '{/Oblique x} (pc)' font 'Helvetica,22'"
        write(unit=2,fmt=1)  " set ylabel '{/Oblique y} (pc)' offset -1.0,0.0 font 'Helvetica,22'"
        write(unit=2,fmt=1)  " set cblabel '{/Oblique n}_{/=14 dust} (cm^{-3})' offset 3.0,0.0 font 'Helvetica,22'"
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " set xrange [0:",(length_x*unit_length*m_in_pc),"]"
        write(unit=2,fmt=1)  " set yrange [0:",(length_y*unit_length*m_in_pc),"]"
        write(unit=2,fmt=1)  " "

        call label_generator(1,lab_step)

        write(unit=2,fmt=1)  " set xtics ",lab_step," font 'Helvetica,18'"
        write(unit=2,fmt=1)  " set ytics ",lab_step," font 'Helvetica,18'"
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " set mxtics 4"
        write(unit=2,fmt=1)  " set mytics 4"
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " unset grid "
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " load '../../../misc/"//theme_colour1//".col'"
  close(unit=2)

  !=============================================================
  ! File with parameters to plot the dust density colourbar
  !

  Do i_dusttype = 1, n_dusttype
    write(i_dusttype_word, '(i1.1)' ) i_dusttype
    open(unit=2, file="./"//Dens_map_folder//"data_dust/gnu_cb_"//i_dusttype_word//".dat", &
            action="write", status="unknown", form="formatted")
        write(unit=2,fmt=1)  " # Parameter to plot dust density colourbar"
        write(unit=2,fmt=1)  " "

        n_asize_reduced = min(n_out_sizes, i_max_ini(i_dusttype))
        Do i_asize_reduced = 1, n_asize_reduced
            write(unit=2,fmt=1)  " if (i == ",i_asize_reduced-1,")   {"

            !=================================================================================================
            ! The data sets can be very large (GBs). To reduce the data amount, store only max. 4 grain sizes.

            if (n_asize_reduced == 1) THEN
                i_asize = i_out_size
            else
                if (i_max_ini(i_dusttype) < n_out_sizes) Then
                    i_asize = i_asize_reduced
                else
                    i_asize = nint((real((i_max_ini(i_dusttype) - 1) * i_asize_reduced)&
                    + real(n_asize_reduced - i_max_ini(i_dusttype)) )&
                    / real(n_asize_reduced - 1))
                end if
            end if
            !=================================================================================================

            cb_min =  max(densi_min_tot(i_dusttype, i_asize), densi_max_tot(i_dusttype, i_asize) * 0.25E-3)
            cb_max =  densi_max_tot(i_dusttype, i_asize)

            IF (cb_max< 10.0*cb_min) THEN
                write(unit=2,fmt=1)  " set cbrange [",cb_min,":",cb_max,"]"
                write(unit=2,fmt=1)  " set cbtics (",cb_min,",",cb_min/1.5 + cb_max/3.0,",",&
                    &  cb_min/3.0 + cb_max/1.5,",", cb_max,") font 'Helvetica,18'"
                write(unit=2,fmt=1)  " unset logscale cb "
                write(unit=2,fmt=1)  " set format cb  '%4.2t{/Symbol \327}10^{%L}'"
            ELSE
                IF ((l_scenario == 1) .or. (l_scenario == 3)) THEN
                    write(unit=2,fmt=1)  " set cbrange [",cb_min,":",cb_max,"]"
                    write(unit=2,fmt=1)  " set cbtics font 'Helvetica,18'"
                    write(unit=2,fmt=1)  " set logscale cb "
                    write(unit=2,fmt=1)  " set format cb '10^{%L}'  "
                    write(unit=2,fmt=1)  " set mcbtics 10"
                ELSE IF (l_scenario == 2) THEN
                    write(unit=2,fmt=1)  " set cbrange [0:",cb_max,"]"
                    write(unit=2,fmt=1)  " set cbtics (0,", cb_max/3.0,",", cb_max/1.5,",", cb_max,") font 'Helvetica,18'"
                    write(unit=2,fmt=1)  " unset logscale cb "
                    write(unit=2,fmt=1)  " set format cb  '%4.2t{/Symbol \327}10^{%L}'"
                ELSE
                    print*, "Error. Scenario not defined."
                END IF
            END IF
            write(unit=2,fmt=1)  "}"
            write(unit=2,fmt=1)  " "
        END DO
    close(unit=2)
  END DO

  !
  !=============================================================

  Do i_timepoints = 1,n_timepoints
    Do i_dusttype = 1, n_dusttype

        write (*,'(A,I3,A)') '      B) Dust flow data: ', int((real(i_timepoints-1)*n_dusttype+i_dusttype)/ &
            real((n_timepoints)*n_dusttype)*100.0), ' % done.' //char(27)//'[A'

        write(i_timepoints_word, '(i5.5)' ) i_timepoints
        write(i_dusttype_word, '(i1.1)' ) i_dusttype
        write(t_timepoints_word, '(f9.2)' ) (time(i_timepoints)) * (time_unit_scaler)

        !=================================================================================================
        ! The data sets can be very large (GBs). To reduce the data amount, store only max. 4 grain sizes.
        n_asize_reduced = min(n_out_sizes, i_max_ini(i_dusttype))

        Do i_asize_reduced = 1, n_asize_reduced
            if (n_asize_reduced == 1) THEN
                i_asize = i_out_size
            else
                if (i_max_ini(i_dusttype) < n_out_sizes) Then
                    i_asize = i_asize_reduced
                else
                    i_asize = nint((real((i_max_ini(i_dusttype) - 1) * i_asize_reduced)&
                    + real(n_asize_reduced - i_max_ini(i_dusttype)) )&
                    / real(n_asize_reduced - 1))
                end if
            end if

            write(a_word(i_asize), '(f6.1)' ) a(i_asize, i_dusttype)*1.0E+9
            if (a(i_asize, i_dusttype)>0.999* 10.0**(-5)) THEN
              write(a_word(i_asize), '(f6.0)' ) a(i_asize, i_dusttype)*1.0E+9
            end if
        END DO

        !========================================================
        ! Plot the dust number density as a function of time

       open(unit=2, file="./"//Dens_map_folder//"data_dust/plot_Density_"//i_dusttype_word//"_"&
                //i_timepoints_word//".dat", action="write", status="unknown", form="formatted")
            write(unit=2,fmt=1)  " reset"
            write(unit=2,fmt=1)  " set output '|ps2pdf - ../pics_dust/Density_"&
                &//i_dusttype_word//"_"//i_timepoints_word//".pdf' "
            write(unit=2,fmt=1)  " load 'gnu.dat' # load parameter to plot dust density"
            write(unit=2,fmt=1)  " "
            ! For each dust grain size one single plot
            Do i_asize_reduced = 1, n_asize_reduced
                if (n_asize_reduced == 1) THEN
                    i_asize = i_out_size
                else
                    if (i_max_ini(i_dusttype) < n_asize_reduced + 1) Then
                        i_asize = i_asize_reduced
                    else
                        i_asize = nint((real((i_max_ini(i_dusttype) - 1) * i_asize_reduced) + &
                        real(n_asize_reduced - i_max_ini(i_dusttype)) )/ real(n_asize_reduced - 1))
                    end if
                end if

                write(unit=2,fmt=1)  " i = ", i_asize_reduced - 1
                write(unit=2,fmt=1)  " load 'gnu_cb_"//i_dusttype_word//".dat' # load parameter to plot dust density colourbar"
                write(unit=2,fmt=1)  " set  title '{/Oblique t} = "//t_timepoints_word//" "//time_word// &
                    & ", ", Material_name(i_dusttype),", {/Oblique a} = "//a_word(i_asize)//" nm' offset 0,-0.8 font 'Helvetica,18'"

                write(unit=2,fmt=1)  " plot  \"
                write(unit=2,fmt=1)  " 'Density_"//i_dusttype_word//"_"&
                    //i_timepoints_word//".dat'  u (",&
            &  (length_x*unit_length*m_in_pc)/real(n_xgrid)," * (($2)+0.5)):(",&
            &  (length_y*unit_length*m_in_pc)/real(n_ygrid)," * (($1)+0.5)):3 index", i_asize_reduced-1, " matrix w image  "
                write(unit=2,fmt=1) ""

            END DO
        close(unit=2)

    End do
  END DO

  print*, ''

  DEALLOCATE(a_word)

end subroutine output_dustdata_plot ! #s3.3


!################################################################################################################
!################################################################################################################

subroutine output_dustdistributiondata_plot ! #s3.4
  use omp_lib
  use datatype
  use variable

  implicit none
  INTEGER :: i_dusttype

  Real (kind=r2), DIMENSION(:), ALLOCATABLE :: alpha_rad_min, alpha_rad_max

  ALLOCATE(alpha_rad_min(1:n_dusttype), alpha_rad_max(1:n_dusttype))

  !========================================================
  ! 1. Plot the dust mass as a function of time for each material ( 2a), 2b), ..., 2c) )

  1 format (*(G0,1X))
  open(unit=2, file="./"//Size_dist_folder//"data/plot_Evolution", action="write", status="unknown", form="formatted")
    write(unit=2,fmt=1)  " reset"
    write(unit=2,fmt=1)  " set terminal postscript enhanced color"
    write(unit=2,fmt=1)  " set output '|ps2pdf - ../pics/Evolution_total.pdf' "
    write(unit=2,fmt=1)  " set encoding iso_8859_1"
    write(unit=2,fmt=1)  " ################"
    write(unit=2,fmt=1)  " set size 0.95,0.98"
    write(unit=2,fmt=1)  " unset key"
    write(unit=2,fmt=1)  " set origin 0.03,0.01"
    write(unit=2,fmt=1)  " ################"
    write(unit=2,fmt=1)  " set xlabel 't ("//time_word//")' font 'Helvetica,22'"
    write(unit=2,fmt=1)  " set ylabel '{/Symbol-Oblique h} = {/Oblique M} / &
                    &{/Oblique M}_{0} (per cent)' offset -2,0 font 'Helvetica,22' # {/Symbol \045}"
    write(unit=2,fmt=1)  " "
    write(unit=2,fmt=1)  " set xrange [*:*]"
    write(unit=2,fmt=1)  " set yrange [0:103]"
    write(unit=2,fmt=1)  " "
    write(unit=2,fmt=1)  " set xtics font 'Helvetica,18'"
    write(unit=2,fmt=1)  " unset ytics"
    write(unit=2,fmt=1)  " set ytics font 'Helvetica,18'"
    write(unit=2,fmt=1)  " "
    write(unit=2,fmt=1)  " set mxtics 10"
    write(unit=2,fmt=1)  " set mytics 10"
    write(unit=2,fmt=1)  " "

    !###############
    !2a), 2b) For each dust material one single plot
    Do i_dusttype = 1, n_dusttype
        write(unit=2,fmt=1)  " set title '", Material_name(i_dusttype), "' offset 0,-0.8 font 'Helvetica,18'"
        write(unit=2,fmt=1)  " plot 'Dustmass_evolution_total.dat' u ($1):(100.0*$",&
            & i_dusttype*2,")  w l  lc ",8-i_dusttype," lw 6  dt 1"
        write(unit=2,fmt=1)  " "
    END Do

    If (n_dusttype > 1) THEN
        ! 2c) Multiplot with all dust materials
        write(unit=2,fmt=1)  " set title 'All materials'  font 'Helvetica,18'"
        write(unit=2,fmt=1)  " plot  \"
        Do i_dusttype = 1, n_dusttype
            write(unit=2,fmt=1)  " 'Dustmass_evolution_total.dat' u ($1):(100.0*$",&
                & i_dusttype*2,") w l lc ",8-i_dusttype," lw 6  dt 1,\"
        END Do
            write(unit=2,fmt=1) "x+200"
            write(unit=2,fmt=1)  " "
    END if
    !###############

    !###############
    ! In scenario 2 and 3, consider only the region within region_r around the origin origin_x, origin_y.
    IF ((l_scenario == 2) .or. (l_scenario == 3)) THEN      ! FK, 15/04/2021, 19/10/2021

        !2a), 2b) For all dust materials one single plot
        Do i_dusttype = 1, n_dusttype
            write(unit=2,fmt=1)  " set  title ' inner", real(region_r*m_in_pc)," pc '&
                & offset 0,-0.8  font 'Helvetica,18'"
            write(unit=2,fmt=1)  " plot 'Dustmass_evolution_total_reg.dat' u ($1):(100.0*$",&
                & 1*2,") w l lc ",8-1," lw 6 dt 1"
        END Do
    END IF
    !###############


    !###############
    if ((n_gg > 0)) then
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " ################"
        write(unit=2,fmt=1)  " set key box inside top right width -7 spacing 1.2  font ',20'"
        write(unit=2,fmt=1)  " ################"
        write(unit=2,fmt=1)  " set yrange [*:*]"
        write(unit=2,fmt=1)  " set ylabel '{/Oblique N}_{collide}/{/Oblique t}   (yr^{-1})' &
            & offset -2,0 font 'Helvetica,22'"
        !write(unit=2,fmt=1)  " set logscale y"
        !write(unit=2,fmt=1)  " set format y '10^{%L}'"
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " V_cell =", V_cell
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " plot 'Dustmass_evolution_total.dat' u ($1):(V_cell *($", &
            1 + n_dusttype*2 +2, ")) w l lc rgb 'yellow4' lw 8 dt 1 t 'Vaporization   ' ,\"
        write(unit=2,fmt=1)  "      'Dustmass_evolution_total.dat' u ($1):(V_cell *($", &
            1 + n_dusttype*2 +3, ")) w l lc rgb 'dark-yellow' lw 8 dt 1 t 'Fragmentation' ,\"
        write(unit=2,fmt=1)  "      'Dustmass_evolution_total.dat' u ($1):(V_cell *($", &
            1 + n_dusttype*2 +4, ")) w l lc rgb 'goldenrod' lw 8 dt 1 t 'Bouncing        ' ,\"
        write(unit=2,fmt=1)  "      'Dustmass_evolution_total.dat' u ($1):(V_cell *($", &
            1 + n_dusttype*2 +5, ")) w l lc rgb 'sandybrown' lw 8 dt 1 t 'Sticking          '"
        write(unit=2,fmt=1)  " "

        IF (n_timepoints > 3) THEN !  approximation splines, need at least 4 points
            write(unit=2,fmt=1)  " plot 'Dustmass_evolution_total.dat' u ($1):(V_cell *($", &
                1 + n_dusttype*2 +2, ")) smooth acsplines lc rgb 'yellow4' lw 8 dt 1 t 'Vaporization   ' ,\"
            write(unit=2,fmt=1)  "      'Dustmass_evolution_total.dat' u ($1):(V_cell *($", &
                1 + n_dusttype*2 +3, ")) smooth acsplines lc rgb 'dark-yellow' lw 8 dt 1 t 'Fragmentation' ,\"
            write(unit=2,fmt=1)  "      'Dustmass_evolution_total.dat' u ($1):(V_cell *($", &
                1 + n_dusttype*2 +4, ")) smooth acsplines lc rgb 'goldenrod' lw 8 dt 1 t 'Bouncing        ' ,\"
            write(unit=2,fmt=1)  "      'Dustmass_evolution_total.dat' u ($1):(V_cell *($", &
                1 + n_dusttype*2 +5, ")) smooth acsplines lc rgb 'sandybrown' lw 8 dt 1 t 'Sticking          '"
        END IF
    end if
    !###############
    close(unit=2)

  !=============================================================
  ! 2. File with parameters to plot the dust density
  !
    open(unit=2, file="./"//Size_dist_folder//"data/gnu.dat", action="write", status="unknown", form="formatted")
    write(unit=2,fmt=1)  " # Parameter to plot grain size distribution "
    write(unit=2,fmt=1)  " set terminal postscript enhanced color "
    write(unit=2,fmt=1)  " set encoding iso_8859_1 "
    write(unit=2,fmt=1)  " ################"
    write(unit=2,fmt=1)  " set size 0.95,0.98"
    write(unit=2,fmt=1)  " unset key"
    write(unit=2,fmt=1)  " set origin 0.03,0.01"
    write(unit=2,fmt=1)  " ################"
    write(unit=2,fmt=1)  " set logscale xy"
    write(unit=2,fmt=1)  " "
    write(unit=2,fmt=1)  " set xlabel '{/Oblique a} (nm)' font 'Helvetica,22'"
    write(unit=2,fmt=1)  " set ylabel 'Number of dust grains (m^{-3} nm^{-1})' offset -2,0 font 'Helvetica,22'"
    write(unit=2,fmt=1)  " "
    write(unit=2,fmt=1)  " set format y '10^{%L}'  "
    write(unit=2,fmt=1)  " set xtics font 'Helvetica,18'"
    write(unit=2,fmt=1)  " set ytics font 'Helvetica,18'"
    write(unit=2,fmt=1)  " "
    write(unit=2,fmt=1)  " set mxtics 10"
    write(unit=2,fmt=1)  " set mytics 10"
    write(unit=2,fmt=1)  " "
    write(unit=2,fmt=1)  " delta_rad = ", delta_rad
    write(unit=2,fmt=1)  " Fak = 1.0E+9"
  close(unit=2)
  !
  !=============================================================

  Do i_timepoints = 1,n_timepoints

    write (*,'(A,I3,A)') '      C) Grain size distribution data:  ',&
        int(i_timepoints/real(n_timepoints)*100.0), ' % done.' //char(27)//'[A'

    write(i_timepoints_word, '(i5.5)' ) i_timepoints
    write(t_timepoints_word, '(f9.2)' ) (time(i_timepoints)) * (time_unit_scaler)

    !========================================================
    ! 4. Plot the dust grain size distribution for each single time and material ( 4a), 4b), ..., 4c) )

    open(unit=2, file="./"//Size_dist_folder//"data/plot_Particlenumbers_"//i_timepoints_word//".dat",&
        action="write", status="unknown", form="formatted")
    write(unit=2,fmt=1)  " set output '|ps2pdf - ../pics/Particlenumbers_"//i_timepoints_word//".pdf' "

    !4a), 4b) For each dust material one single plot
    Do i_dusttype = 1, n_dusttype
       write(unit=2,fmt=1)  ""
       write(unit=2,fmt=1)  " load 'gnu.dat'"
       write(unit=2,fmt=1)  " ################"
       write(unit=2,fmt=1)  " set title '{/Oblique t} = "//t_timepoints_word//" "//time_word// &
            & ", ", Material_name(i_dusttype),"' offset 0,-0.8 font 'Helvetica,18'"
       write(unit=2,fmt=1)  " set xrange [",0.7*a(1,i_dusttype)*1.0E+9*delta_rad**(-0.5),":",&
            & 1.25*a(n_asize,1)*1.0E+9*delta_rad**(0.5),"]"
       write(unit=2,fmt=1)  " set yrange [",ymin(i_dusttype),":",ymax(i_dusttype),"]"
       write(unit=2,fmt=1)  " plot  'Particlenumbers_00001.dat' u (Fak*$",i_dusttype*2-1,"):(($",i_dusttype*2,"&
            & )+1.0E-300== 0 ? NaN : $",i_dusttype*2,") w boxes  lc rgb 'grey' lw 3  dt 1 , &
            &'Particlenumbers_"//i_timepoints_word//".dat' u (Fak*$", i_dusttype*2-1,"):(($",i_dusttype*2,"&
            &)+1.0E-300 == 0 ? NaN : $",i_dusttype*2,") w boxes  lc ",8-i_dusttype,"  lw 6  dt 1  "
       write(unit=2,fmt=1)  " ################"
       write(unit=2,fmt=1)  ""
       write(unit=2,fmt=1)  " set key inside right  width 0 spacing 1.3  font ',22'"
       write(unit=2,fmt=1)  " set yrange [*:*]"
       write(unit=2,fmt=1)  " set ylabel 'd{/Symbol-Oblique r}/d{/Oblique a} [10^{-27} kg m^{-3} nm^{-1}]' &
                             & offset -3,0 font 'Helvetica,22'"
       write(unit=2,fmt=1)  " Fak2 = ", 4.0/3.0* pi * rho_bulk(i_dusttype)
       write(unit=2,fmt=1)  " plot  'Particlenumbers_00001.dat' u (Fak*$",i_dusttype*2-1,"):(($",i_dusttype*2,"&
            & )+1.0E-300== 0 ? NaN : (Fak2*(Fak*$1)**3.0)*$",i_dusttype*2,") w l  lc rgb 'grey' lw 10  dt 1 t 'initial',&
            &'Particlenumbers_"//i_timepoints_word//".dat' u (Fak*$", i_dusttype*2-1,"):(($",i_dusttype*2,"&
            &)+1.0E-300 == 0 ? NaN : (Fak2*(Fak*$1)**3.0)*$",i_dusttype*2,") w l lc ",8-i_dusttype,"  lw 10  dt 1 t 'final'"
    END DO

    ! 4c) Multiplot with all dust materials
    If (n_dusttype > 1) THEN
        Do i_dusttype = 1, n_dusttype
            alpha_rad_min(i_dusttype) = a(1,i_dusttype)*delta_rad**(-0.5)
            alpha_rad_max(i_dusttype) = a(n_asize,i_dusttype)*delta_rad**(0.5)
        END DO

        write(unit=2,fmt=1)  ""
        write(unit=2,fmt=1)  " load 'gnu.dat'"
        write(unit=2,fmt=1)  " ################"
        write(unit=2,fmt=1)  " set  title '{/Oblique t} = "//t_timepoints_word//" "//time_word//",&
            & all materials' offset 0,-0.8 font 'Helvetica,18'"
        write(unit=2,fmt=1)  " set xrange [",0.7*minval(alpha_rad_min(:))*1.0E+9,":",&
            & 1.25*maxval(alpha_rad_max(:))*1.0E+9,"]"
        write(unit=2,fmt=1)  " set yrange [", minval(ymin(:)),":", maxval(ymax(:)),"]"

        write(unit=2,fmt=1)  " plot  \"
        Do i_dusttype = 1, n_dusttype
            write(unit=2,fmt=1)  "'Particlenumbers_00001.dat' u (Fak*delta_rad**(0.2*",(i_dusttype-1),")*$"&
            &,i_dusttype*2-1,"):($",i_dusttype*2," == 0 ? NaN : $",i_dusttype*2,") w boxes  lc rgb 'grey'&
            & lw 3  dt ",i_dusttype,", 'Particlenumbers_"//i_timepoints_word//".dat' u &
            & (Fak*delta_rad**(0.2*",(i_dusttype-1),")*$", i_dusttype*2-1,"):($",i_dusttype*2,"&
            & == 0 ? NaN : $",i_dusttype*2,") w boxes  lc ",8-i_dusttype,"  lw 6  dt 1  ,\"
        END Do
        write(unit=2,fmt=1)  "''"
        write(unit=2,fmt=1)  " ################"
        write(unit=2,fmt=1)  ""
        write(unit=2,fmt=1)  " set key inside right  width 0 spacing 1.3  font ',22'"
        write(unit=2,fmt=1)  " set yrange [*:*]"
        write(unit=2,fmt=1)  " set ylabel 'd{/Symbol-Oblique r}/d{/Oblique a} (10^{-27} kg m^{-3} nm^{-1})' &
                             & offset -3,0 font 'Helvetica,22'"
        write(unit=2,fmt=1)  " plot  \"
        Do i_dusttype = 1, n_dusttype
            write(i_dusttype_word, '(i1.1)' ) i_dusttype
            write(unit=2,fmt=1)  "'Particlenumbers_00001.dat' u (Fak*$"&
            &,i_dusttype*2-1,"):($",i_dusttype*2," == 0 ? NaN : (Fak2*(Fak*$1)**3.0)*$",i_dusttype*2,") w l&
            & lc rgb 'grey' lw 3  dt ",i_dusttype,"  t 'initial, mat. "//i_dusttype_word//"', &
            & 'Particlenumbers_"//i_timepoints_word//".dat' u &
            & (Fak*$", i_dusttype*2-1,"):($",i_dusttype*2,"&
            & == 0 ? NaN : (Fak2*(Fak*$1)**3.0)*$",i_dusttype*2,") w l  lc ",8-i_dusttype,"&
            & lw 6  dt 1 t 'final, mat. "//i_dusttype_word//"',\"
        END Do
            write(unit=2,fmt=1) "'' t ''"
            write(unit=2,fmt=1)  " ################"
    End if
    close(unit=2)

    !========================================================

    END DO

    print*, ''

    DEALLOCATE(alpha_rad_min, alpha_rad_max)

end subroutine output_dustdistributiondata_plot ! #s3.4

!################################################################################################################

subroutine output_movies ! #s3.5
  use omp_lib
  use datatype
  use variable

  implicit none
  INTEGER :: i_dusttype, i_asize

  REAL :: Movielength

  INTEGER :: i_asize_reduced, n_asize_reduced, dummy_int

!   Movielength = 200.0/real(n_timepoints) !min((n_timepoints-1.0)/30.0, 20.0) ! in 0.1 Sekunden
  Movielength =  max(3.0, real(n_timepoints)/60.0) !min((n_timepoints-1.0)/30.0, 20.0)    ! in frames per second

  !========================================================
  ! 1. Generate movie of gas flow

  1 format (*(G0,1X))
  open(unit=2, file="./"//Dens_map_folder//"data_gas/pro_movie.sh", action="write", status="unknown", form="formatted")
    write(unit=2,fmt=1)  " mogrify -density 500 -flatten -resize 800 -format png ./../pics_gas/Density_*.pdf[000]"
    write(unit=2,fmt=1)  " ffmpeg -framerate ", Movielength, " -pattern_type glob -i '../pics_gas/*.png'&
                            & -c:v libx264 ../pics_gas/Density_gas.mp4"
    write(unit=2,fmt=1)  " rm ../pics_gas/Density_*.png"
    write(unit=2,fmt=1)  " "
    write(unit=2,fmt=1)  " mogrify -density 500 -flatten -resize 800 -format png ./../pics_gas/Density_*.pdf[001]"
    write(unit=2,fmt=1)  " ffmpeg -framerate ", Movielength, " -pattern_type glob -i '../pics_gas/*.png'&
                            & -c:v libx264 ../pics_gas/Temperature_gas.mp4"
    write(unit=2,fmt=1)  " rm ../pics_gas/Density_*.png"

  IF  (l_MHD == 1) THEN
        write(unit=2,fmt=1)  " "
        write(unit=2,fmt=1)  " mogrify -density 500 -flatten -resize 800 -format png ./../pics_gas/Density_*.pdf[002]"
        write(unit=2,fmt=1)  " ffmpeg -framerate ", Movielength, " -pattern_type glob -i '../pics_gas/*.png'&
                            & -c:v libx264 ../pics_gas/Magnetic_field.mp4"
        write(unit=2,fmt=1)  " rm ../pics_gas/Density_*.png"
    END IF

  close(unit=2)

  !========================================================
  ! 2. Generate movie of dust flow

  open(unit=2, file="./"//Dens_map_folder//"data_dust/pro_movie.sh", action="write", status="unknown", form="formatted")
    Do  i_dusttype = 1, n_dusttype
        write(i_dusttype_word, '(i1.1)' ) i_dusttype
        !=================================================================================================
        ! The data sets can be very large (GBs). To reduce the data amount, store only max. 4 grain sizes.
        n_asize_reduced = min(n_out_sizes, i_max_ini(i_dusttype))
        Do i_asize_reduced = 1, n_asize_reduced
            i_asize = i_asize_reduced
            write(i_asize_word, '(i1.1)' ) (i_asize-1)

            if (i_asize-1 < 10) THEN
              i_asize_word = "00"//i_asize_word
            else if ((9 < i_asize-1) .and. (i_asize-1 < 100)) THEN
              i_asize_word = "0"//i_asize_word
            else if (999 < i_asize-1) THEN
                print*, 'Paperboats ERROR 9, number of grain sizes too large.'
                STOP
            END IF

            write(unit=2,fmt=1)  " mogrify -density 500 -flatten -resize 800 -format png &
                &./../pics_dust/Density_"//i_dusttype_word//"*.pdf["//i_asize_word//"]"

            write(i_asize_word, '(i1.1)' ) i_asize
            if (i_asize < 10) THEN
              i_asize_word = "00"//i_asize_word
            else if ((9 < i_asize) .and. (i_asize < 100)) THEN
              i_asize_word = "0"//i_asize_word
            else if (999 < i_asize) THEN
                print*, 'Paperboats ERROR 10, number of grain sizes too large.'
                STOP
            END IF

!             write(unit=2,fmt=1)  " convert -delay ", Movielength, " ../pics_dust/*.png &
!                 &../pics_dust/Density_dust_"//i_dusttype_word//"_"//i_asize_word//".mp4"
            write(unit=2,fmt=1)  " ffmpeg -framerate ", Movielength, " -pattern_type glob -i '../pics_dust/*.png'&
                            & -c:v libx264 ../pics_dust/Density_dust_"//i_dusttype_word//"_"//i_asize_word//".mp4"
            write(unit=2,fmt=1)  " rm ../pics_dust/Density_*.png"
            write(unit=2,fmt=1)  " "
        END DO
    END DO
  close(unit=2)

  !========================================================
  ! 3. Generate movie of dust grain distribution

  open(unit=2, file="./"//Size_dist_folder//"data/pro_movie.sh", action="write", status="unknown", form="formatted")

    IF (n_dusttype > 1) THEN
       dummy_int = 2*n_dusttype + 1
    ELSE
       dummy_int = 1
    END IF
        Do i_dusttype = 0, dummy_int
        write(i_dusttype_word, '(i1.1)' ) i_dusttype

        write(unit=2,fmt=1)  " mogrify -density 500 -flatten -resize 800 -format png &
            & ./../pics/Particlenumbers_*.pdf["//i_dusttype_word//"]"

        write(i_dusttype_word, '(i1.1)' ) i_dusttype+1
!         write(unit=2,fmt=1)  " convert -delay ", Movielength, " ../pics/*.png ../pics/Particlenumbers__"//i_dusttype_word//".mp4"
        write(unit=2,fmt=1)  " ffmpeg -framerate ", Movielength, " -pattern_type glob -i '../pics/*.png'&
                            & -c:v libx264  ../pics/Particlenumbers__"//i_dusttype_word//".mp4"
        write(unit=2,fmt=1)  " rm ../pics/Particlenumbers_*.png"
        write(unit=2,fmt=1)  " "
    END DO
  close(unit=2)

end subroutine output_movies ! #s3.5
!################################################################################################################

end module output_plot

