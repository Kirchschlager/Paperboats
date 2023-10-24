module tools
  use omp_lib
  use datatype
  use variable
  
  implicit none
    !--------------------------------------------------------------------------!
!     PRIVATE
!         CHARACTER(len=4), PARAMETER :: file_a  = "huhu"
    !--------------------------------------------------------------------------!
    PUBLIC :: parse_real, parse_char, Bin_search, label_generator
    !--------------------------------------------------------------------------!
CONTAINS 

!################################################################################################################

subroutine parse_real(keyword, var, word_number, data_file)
    use omp_lib
    use datatype
    implicit none
    ! This routine seaches each line for the keyword and returns the value as a string

    !========================================================    
    CHARACTER(len=*),INTENT(IN)      :: keyword, data_file
    REAL(kind=r2),INTENT(INOUT)      :: var     
    
    CHARACTER(len=256)               :: line, val_st
    INTEGER                          :: io, key_found, word_number 
    INTEGER                          :: i, line_end, border_l, border_r, counter
    !========================================================

    open(unit=1, file=data_file, action="read", status="unknown", form="formatted")
        key_found = 0
        line = ''
        DO WHILE (key_found == 0)
           
            read(unit=1,fmt='(A)',iostat=io) line  ! read the next line of the input file
        
            ! find end of file
            IF (io < 0) THEN 
                EXIT 
            END IF
            line_end = len_trim(line)
        
            ! search for the '=' symbol in this line
            DO border_l=1,line_end
                IF (line(border_l:border_l) == '=') THEN
                    EXIT
                END IF
            END DO

            ! test if the keyword belongs to this entry
            IF (line(1:border_l-1) == keyword) THEN
                key_found = 1

                counter = 0
                DO i = border_l+1,line_end
                   IF (((line(i:i) == ',') .or. (line(i:i) == '!'))) THEN
                      counter = counter + 1
                      IF (counter == word_number) THEN
                          EXIT
                      ELSE
                         border_l = i
                      END IF
                   END IF
                END DO
                border_r = i
                
                IF (border_l==border_r) THEN
                    key_found = 0
                    EXIT
                ELSE                     
                    !now set the file value (between border_l and border_r) to the given var
                    val_st = line(border_l+1:border_r)
                    read(val_st, fmt=*) var         ! Gives an error, when val_st is empty space
                END IF
            END IF
        END DO
        
        IF (key_found == 0) THEN
                print*, 'Paperboats ERROR 11 in input file ('//TRIM(data_file)//') with keyword: ', keyword
                STOP
        END IF
        
    close(unit=1)

end subroutine parse_real


!################################################################################################################

subroutine parse_char(keyword, var, word_number, data_file)
    use omp_lib
    use datatype
    implicit none
    ! This routine seaches each line for the keyword and returns the value as a string

    !========================================================    
    CHARACTER(len=*),INTENT(IN)      :: keyword, data_file
    CHARACTER(len=7),INTENT(INOUT)   :: var       
    
    CHARACTER(len=256)               :: line, val_st
    INTEGER                          :: io, key_found, word_number 
    INTEGER                          :: i, line_end, border_l, border_r, counter
    !========================================================

    open(unit=1, file=data_file, action="read", status="unknown", form="formatted")
        key_found = 0
        line = ''
        DO WHILE (key_found == 0)
           
            read(unit=1,fmt='(A)',iostat=io) line  ! read the next line of the input file
        
            ! find end of file
            IF (io < 0) THEN 
                EXIT 
            END IF
            line_end = len_trim(line)
        
            ! search for the '=' symbol in this line
            DO border_l=1,line_end
                IF (line(border_l:border_l) == '=') THEN
                    EXIT
                END IF
            END DO

            ! test if the keyword belongs to this entry
            IF (line(1:border_l-1) == keyword) THEN
                key_found = 1

                counter = 0
                DO i = border_l+1,line_end
                   IF (((line(i:i) == ',') .or. (line(i:i) == '!'))) THEN
                      counter = counter + 1
                      IF (counter == word_number) THEN
                          EXIT
                      ELSE
                         border_l = i
                      END IF
                   END IF
                END DO
                border_r = i
                
                IF (border_l==border_r) THEN
                    key_found = 0
                    EXIT
                ELSE                     
                    !now set the file value (between border_l and border_r) to the given var
                    val_st = line(border_l+1:border_r)
                    read(val_st, fmt=*) var         ! Gives an error, when val_st is empty space
                END IF
            END IF
        END DO
        
        IF (key_found == 0) THEN
                print*, 'Paperboats ERROR 11 in input file ('//TRIM(data_file)//') with keyword: ', keyword
                STOP
        END IF
        
    close(unit=1)

end subroutine parse_char

!################################################################################################################
!################################################################################################################

subroutine Bin_search(i_dusttype, a_search, i_floor, P_left_bin)
    use omp_lib
    use datatype
    use variable
                          
    implicit none

    INTEGER, intent(in)        :: i_dusttype   ! input    
    Real(kind=r2), intent(in)  :: a_search     ! input
    
    INTEGER, intent(out)       :: i_floor      ! output
    Real(kind=r2), intent(out) :: P_left_bin   ! output
     
    INTEGER :: i
    
    i_floor = 0
    Do i = 1, n_asize
        if (a_search >= a(n_asize + 1 - i, i_dusttype)) then
            i_floor = n_asize + 1 - i
            exit        
        End if
    END DO
    
    !====================  
    ! Mass conservation by P_left_bin: P_left_bin in the left bin, (1.0 - P_left_bin) in the right bin.
    !===================== 
    if ((i_floor > 0) .and. (i_floor < n_asize)) THEN 
        P_left_bin = ((a(i_floor+1, i_dusttype))**3.0_r2 - a_search**3.0_r2)/&
            ((a(i_floor+1, i_dusttype))**3.0_r2 - (a(i_floor, i_dusttype))**3.0_r2)
    elseif (i_floor == 0) THEN
        if (a_search > a(1, i_dusttype)/delta_rad) THEN
            P_left_bin = (1.0_r2 - (a_search/a(1, i_dusttype))**3.0_r2)/&
            (1.0_r2 - delta_rad**(-3.0_r2))
        else
            P_left_bin = 1.0_r2
        end if 
    else ! i_floor = n_asize 
        if (a_search < (a(n_asize, i_dusttype) * delta_rad)) THEN
            P_left_bin = (delta_rad**3.0_r2 - (a_search/a(n_asize, i_dusttype))**3.0_r2)/&
            (delta_rad**3.0_r2 - 1.0_r2)
        else
            P_left_bin = 0.0_r2
        end if     
    end if

    !=====================     
    ! Test P_left_bin
    !===================== 
    if ((P_left_bin < 0.0) .or. (P_left_bin > 1.0)) THEN
        print*, 'Paperboats ERROR 6, P_left_bin not well-defined:', P_left_bin, i_floor, &
            a(i_floor, i_dusttype), a_search, a(i_floor+1, i_dusttype)  
        STOP
    end if    

end subroutine Bin_search

!################################################################################################################
!################################################################################################################

subroutine label_generator(i_coor, lab_step) 
    use omp_lib
    use datatype
    use variable
                          
    implicit none

    INTEGER, intent(in)        :: i_coor         ! input        
    Real(kind=r2), intent(out) :: lab_step       ! output
   
    Integer       :: i_dummy
    Real(kind=r2) :: var1, var2, factor
     
    var1 =  data_grid(n_xgrid, n_ygrid, n_zgrid,i_coor) * m_in_pc

    IF (var1 > 1.0) THEN
        Do i_dummy = 0,15 
            IF ((var1   > 10.0**(real(i_dummy))) .and. (var1 <= 10.0**(real(i_dummy+1)))) THEN
                factor = 10.0**(real(i_dummy))
            EXIT    
            END IF
        END DO
    ELSE
        Do i_dummy = 0,15 
            IF ((var1   <= 0.1**(real(i_dummy))) .and. (var1 > 0.1**(real(i_dummy+1)))) THEN
                factor = 0.1**(real(i_dummy+1))
            EXIT    
            END IF
        END DO
    END IF
    var2 = var1/factor ! value between 1.0000 < var2 <= 10.0
    IF (var2 < 1.5)      THEN
        lab_step = 0.2
    ELSE IF (var2 < 2.2) THEN
        lab_step = 0.3    
    ELSE IF (var2 < 3.0) THEN
        lab_step = 0.4
    ELSE IF (var2 < 3.9) THEN
        lab_step = 0.5        
    ELSE IF (var2 < 4.5) THEN
        lab_step = 0.6    
    ELSE IF (var2 < 6.0) THEN   
        lab_step = 0.8
    ELSE IF (var2 < 7.5) THEN 
        lab_step = 1.0    
    ELSE
        lab_step = 1.5    
    END IF
     
    lab_step = lab_step * factor 
     
end subroutine label_generator

!################################################################################################################

end module tools
