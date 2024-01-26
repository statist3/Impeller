#!/bin/bash

rm joined.stl

rm -rf joined.stl

for filename in *.stl
do
  echo $filename | tee log.sed
  printf '1,$ s/solid.*$/solid '"${filename%\.*}"'/g\nw! '"./${filename%\.*}.stl"'\nq!' | ex - "./${filename%\.*}.stl"
  
  if [[ "$filename" != refine* ]]; 
  then
     cat "./${filename%\.*}.stl" >> "./joined.stl"
     rm "./${filename%\.*}.stl"  # comment line to keep the individual stl files.
  fi
  
done;

