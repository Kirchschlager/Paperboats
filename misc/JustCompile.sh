#!/bin/bash

lin_or_para=2 # 1=linear, 2=parallel

echo 
echo "  Working directory: "$(pwd)
echo

#---------------------------------------------------------------------------------
# Compile

if [ "$lin_or_para" = "1" ]; then
    #h5fc -o paperboats Paperboats.f90
    echo "  Compiling: sequential."    
    echo   
    cd tmp/
    h5fc -c ../src/datatype.f90     -o datatype.o
    h5fc -c ../src/variable.f90     -o variable.o
    h5fc -c ../src/tools.f90        -o tools.o    
    h5fc -c ../src/functs.f90       -o functs.o        
    h5fc -c ../src/output_data.f90  -o output_data.o
    h5fc -c ../src/output_plot.f90  -o output_plot.o  
    h5fc -c ../src/init_astrobear.f90  -o init_astrobear.o         
    h5fc -c ../src/init_pencil.f90  -o init_pencil.o   
    h5fc -c ../src/init_swift.f90   -o init_swift.o 
    h5fc -c ../src/init_arepo.f90   -o init_arepo.o 
    h5fc -c ../src/init_amun.f90   -o init_amun.o       
    h5fc -c ../src/init.f90         -o init.o 
    h5fc -c ../src/processing.f90   -o processing.o  
    h5fc -c ../src/clean.f90        -o clean.o        
    h5fc -c ../src/Paperboats.f90   -o paperboats.o  
   
    h5fc datatype.o variable.o tools.o functs.o output_data.o output_plot.o init_astrobear.o init_pencil.o init_swift.o init_arepo.o  init_amun.o  init.o processing.o clean.o paperboats.o -o paperboats
    cd ../
else
    # h5fc -o paperboats -fopenmp  Paperboats.f90
    echo "  Compiling: parallel."
    echo
    cd tmp/
    h5fc -c -fopenmp ../src/datatype.f90     -o datatype.o
    h5fc -c -fopenmp ../src/variable.f90     -o variable.o
    h5fc -c -fopenmp ../src/tools.f90        -o tools.o    
    h5fc -c -fopenmp ../src/functs.f90       -o functs.o     
    h5fc -c -fopenmp ../src/output_data.f90  -o output_data.o
    h5fc -c -fopenmp ../src/output_plot.f90  -o output_plot.o
    h5fc -c -fopenmp ../src/init_astrobear.f90  -o init_astrobear.o      
    h5fc -c -fopenmp ../src/init_pencil.f90  -o init_pencil.o   
    h5fc -c -fopenmp ../src/init_swift.f90   -o init_swift.o    
    h5fc -c -fopenmp ../src/init_arepo.f90   -o init_arepo.o    
    h5fc -c -fopenmp ../src/init_amun.f90    -o init_amun.o     
    h5fc -c -fopenmp ../src/init.f90         -o init.o        
    h5fc -c -fopenmp ../src/processing.f90   -o processing.o   
    h5fc -c -fopenmp ../src/clean.f90        -o clean.o        
    h5fc -c -fopenmp ../src/Paperboats.f90   -o paperboats.o  
   
    h5fc datatype.o variable.o tools.o functs.o output_data.o output_plot.o init_astrobear.o init_pencil.o init_swift.o init_arepo.o init_amun.o init.o processing.o clean.o paperboats.o -fopenmp -o paperboats
    cd ../
fi

echo " "
echo "  Done."
echo " ______________________________________________________"
# rm  tmp/*
