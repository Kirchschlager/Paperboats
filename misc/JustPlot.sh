#!/bin/bash

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
cd ../../
echo " "
echo "  Done."
echo " ______________________________________________________"
  echo " "
  echo "  Generate movies."
  cd GrainSD/data/
  wait
  chmod u+x pro_movie.sh
  wait
  echo "    A) Grain size distribution."
  ./pro_movie.sh
  wait
  cd ../../Maps/data_gas/
  wait
  chmod u+x pro_movie.sh
  wait
  echo "    B) Gas flow."
  ./pro_movie.sh
  wait
  cd ../data_dust/
  wait
  chmod u+x pro_movie.sh
  wait
  echo "    C) Dust flow."
  ./pro_movie.sh
  wait
  cd ../../../
  echo " "
  echo "  Done."
  echo " ______________________________________________________"
  rm  tmp/*