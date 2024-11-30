#!/bin/bash


for ((j = 20; j <= 20; j += 5));
do
	for ((h = 4; h <= 4; h *= 2));
	do
		for i in {0..9}
		do
			./gedUb singlePair/AIDS${j}/$i $j $h >> ResultsUb/AIDS20_${i}_out
		done
	done
done
