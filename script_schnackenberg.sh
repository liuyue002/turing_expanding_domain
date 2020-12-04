#!/bin/bash
for X in 0.75 0.6 0.5 0.4 0.3 0.2 0.1
do
	#sed -E -i "s/^(growthrate = )0[.][0-9]+(;.*)$/\1${X}\2/g" schnackenburg_fde_2d.m
	#screen -d -m matlab -nodisplay -nosplash -r "run('./schnackenburg_fde_2d.m');exit"
	#screen -d -m matlab -nodisplay -nosplash -r "schnackenburg_2dfunc(${X});exit"
	screen -d -m matlab -nodisplay -nosplash -r "schnackenberg_2dfunc(${X});exit"
done
