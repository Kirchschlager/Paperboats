module datatype
 implicit none
 integer, parameter, public :: r1=selected_real_kind(1)
 integer, parameter, public :: r2=selected_real_kind(9)
end module datatype

PROGRAM Umrechnen
use datatype
IMPLICIT NONE


Real(kind=r2), DIMENSION(:), ALLOCATABLE :: Gasdensity(:,:), Tem(:,:), Mag(:,:)
Real(kind=r2), DIMENSION(:), ALLOCATABLE :: MinDen(:), MaxDen(:),MinTem(:), MaxTem(:)&
    & ,MinMag(:), MaxMag(:), Time(:) 

INTEGER :: nx,ny,j,i,it,nt

CHARACTER(LEN=5) :: it_word

ny= 600
nx= 1500

nt = 125

ALLOCATE(Gasdensity(1:nx,1:ny), Tem(1:nx,1:ny), Mag(1:nx,1:ny))
ALLOCATE(MinDen(1:nt), MaxDen(1:nt),MinTem(1:nt), MaxTem(1:nt),MinMag(1:nt), MaxMag(1:nt), Time(1:nt))
 


PRINT*, "Diese Programm rechnet eine Wertetabelle in eine andere um. Die Tabellen liegen in Dateiformat vor."

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Einlesen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=================== read Time

open(unit=1, file="./Results/GrainSD/data/Dustmass_evolution_total.dat",&
    & action="read",status="unknown", form="formatted")
read(unit=1,fmt=*) 
Do it = 1, nt
    read(unit=1,fmt=*) Time(it) 
END DO    
close(unit=1)

!=================== read and find in data
Do it = 1, nt
    write(it_word, '(i5.5)' ) it
    open(unit=1, file="./Results/Maps/data_gas/Density_"//it_word//".dat",&
        & action="read",status="unknown",form="formatted")
    read(unit=1,fmt=*)     
    Do i=1,nx
        read(unit=1,fmt=*) (Gasdensity(i,j), j=1,ny) 
    END DO    
    read(unit=1,fmt=*)
    read(unit=1,fmt=*)     
    Do i=1,nx
        read(unit=1,fmt=*) (Tem(i,j), j=1,ny) 
    END DO     
    read(unit=1,fmt=*) 
    read(unit=1,fmt=*)     
    Do i=1,nx
        read(unit=1,fmt=*) (Mag(i,j), j=1,ny) 
    END DO   
    
    close(unit=1)
 
    MaxDen(it) = maxval(Gasdensity)
    MinDen(it) = minval(Gasdensity)
    
    MaxTem(it) = maxval(Tem)*1.0E-6
    MinTem(it) = minval(Tem)*1.0E-6
    
    MaxMag(it) = maxval(Mag)
    MinMag(it) = minval(Mag)    
    

    print*, it,Time(it), MinDen(it), MaxDen(it), MinTem(it), MaxTem(it), MinMag(it), MaxMag(it)
END DO
!=================== read and find in data


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Auslesen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(unit=1, file="./neu.dat", action="write", status="unknown", form="formatted")
   Do it = 1, nt
             write(unit=1,fmt=*) Time(it), MinDen(it), MaxDen(it), MaxTem(it), MinMag(it), MaxMag(it)
   End do
close(unit=1)

DEALLOCATE(Gasdensity, Tem, Mag, MinDen, MaxDen,MinTem, MaxTem,MinMag, MaxMag, Time)
END PROGRAM Umrechnen
