#!/bin/bash

lin_or_para=2 # 1=linear, 2=parallel

echo 
echo "  Working directory: "$(pwd)
echo


#---------------------------------------------------------------------------------
# Prepare compilation

if [ "$lin_or_para" = "1" ]; then
    echo "  Compiling: sequential."    
    parallelisation_command=" " 
else
    echo "  Compiling: parallel."
    parallelisation_command=" -fopenmp "     
fi

echo
echo "  Use the compiler:" 
if  command -v h5fc ; then
    compiler_command="h5fc"
else if  command -v h5pfc ; then
        compiler_command="h5pfc"
     else
        echo  '  Error. Neither h5fc nor h5pfc compiler can be found.'   
     fi
fi
echo

#---------------------------------------------------------------------------------
# Compile

cd tmp/
$compiler_command -c $parallelisation_command ../src/datatype.f90     -o datatype.o
$compiler_command -c $parallelisation_command ../src/variable.f90     -o variable.o
$compiler_command -c $parallelisation_command ../src/tools.f90        -o tools.o    
$compiler_command -c $parallelisation_command ../src/functs.f90       -o functs.o     
$compiler_command -c $parallelisation_command ../src/output_data.f90  -o output_data.o
$compiler_command -c $parallelisation_command ../src/output_plot.f90  -o output_plot.o
$compiler_command -c $parallelisation_command ../src/init_astrobear.f90  -o init_astrobear.o      
$compiler_command -c $parallelisation_command ../src/init_pencil.f90  -o init_pencil.o   
$compiler_command -c $parallelisation_command ../src/init_swift.f90   -o init_swift.o    
$compiler_command -c $parallelisation_command ../src/init_arepo.f90   -o init_arepo.o    
$compiler_command -c $parallelisation_command ../src/init_amun.f90    -o init_amun.o     
$compiler_command -c $parallelisation_command ../src/init.f90         -o init.o        
$compiler_command -c $parallelisation_command ../src/processing.f90   -o processing.o   
$compiler_command -c $parallelisation_command ../src/clean.f90        -o clean.o        
$compiler_command -c $parallelisation_command ../src/Paperboats.f90   -o paperboats.o  

$compiler_command datatype.o variable.o tools.o functs.o output_data.o output_plot.o init_astrobear.o init_pencil.o init_swift.o init_arepo.o init_amun.o init.o processing.o clean.o paperboats.o $parallelisation_command -o paperboats
cd ../

#---------------------------------------------------------------------------------
# Check whether the directory "Results" exists or not, otherwise create "Results" and subdirectories
wait

RES="/Results"
FOLDER=$(pwd)$RES
 
if   test -d "$FOLDER"; then
    echo "  Write results to the directory '$RES'."
    echo ""    
else     
    mkdir Results/
    mkdir Results/Maps/
    mkdir Results/GrainSD/
    mkdir Results/Maps/data_dust
    mkdir Results/Maps/data_gas
    mkdir Results/Maps/pics_dust
    mkdir Results/Maps/pics_gas
    mkdir Results/GrainSD/data
    mkdir Results/GrainSD/pics
    echo "  Write results to the directory '$RES'."
    echo ""        
fi

#---------------------------------------------------------------------------------
# Check whether results from previous runs exist and delete them

FOLDER=$(pwd)
FILE="/GrainSD/data/Dustmass_evolution_total.dat"
FILE=$(pwd)$RES$FILE

if test -f "$FILE"; then
    rm Results/*/*/* 
    echo "  Old data and figures in '$RES' deleted."
    echo ""     
fi
#---------------------------------------------------------------------------------

wait
cd tmp/
echo "  Start Paperboats."
./paperboats
cd ../
wait
echo " ______________________________________________________"
echo " "
echo "  Gnuplot: Plot data."
cd Results/GrainSD/data/
wait
  echo "    A) Grain size distribution."
    gnuplot plot_*
  wait
  cd ../../Maps/data_gas/
  wait
  echo "    B) Gas flow."
      gnuplot plot_*
  wait
 cd ../data_dust/
 wait
 echo "    C) Dust flow."
 gnuplot plot_*
 wait
cd ../../../
echo " "
echo "  Plotting finished."
echo " ______________________________________________________"
 echo " "
 echo "  Generate movies."
 cd Results/GrainSD/data/
 wait
 chmod u+x pro_movie.sh
 wait
 echo "    A) Grain size distribution."
#    ./pro_movie.sh
 wait
 cd ../../Maps/data_gas/
 wait
 chmod u+x pro_movie.sh
 wait
 echo "    B) Gas flow."
#      ./pro_movie.sh
 wait
 cd ../data_dust/
 wait
 chmod u+x pro_movie.sh
 wait
 echo "    C) Dust flow."
#     ./pro_movie.sh
 wait
 cd ../../../
 echo " "
 echo "  Done."
 echo " ______________________________________________________"
  rm  tmp/*
