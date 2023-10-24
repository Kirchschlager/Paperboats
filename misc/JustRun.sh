#!/bin/bash

lin_or_para=2 # 1=linear, 2=parallel

echo 
echo "  Working directory: "$(pwd)
echo

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
echo "  Done."
echo " ______________________________________________________"
rm  tmp/*