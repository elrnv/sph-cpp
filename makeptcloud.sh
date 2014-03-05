#!/bin/bash
file=$1

lines=`wc -l $file | awk '{print $1;}'`

temp="tempfile.$$.txt"
for (( i=1; i <= $lines; i++ ))
do
  echo f $i >> $temp
done

cat $temp >> $file
rm -rf $temp
