#!/bin/bash

# a small script that converts a .obj list of vertices to a .obj point clound
# including a list of 1 vertex faces, which is what the Assimp library expects.

file=$1

lines=`wc -l $file | awk '{print $1;}'`

temp="tempfile.$$.txt"
for (( i=1; i <= $lines; i++ ))
do
  echo f $i >> $temp
done

cat $temp >> $file
rm -rf $temp
