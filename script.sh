#!/bin/bash

#primeiro script - ola.sh

function myfunc1()
{
return $./mandelbrot
}


i=1
#a=0
while [ $i -le 10 ];do
  #./mandelbrot
  juca = $(myfunc1 10 20)
  i=$[$i+1]
done

#a = $[$a/10]

#echo ola $USER

#pwd

#date
