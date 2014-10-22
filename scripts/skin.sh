#!/bin/bash
scene=$1

files=`ls -1 ./output/scene$1*.sim`

for file in $files
do
  tempfile="${file##*/}"
  ./ParticleSkinner/particleskinner 0.05 $file ./render/mesh/${tempfile%.*}.obj
done
