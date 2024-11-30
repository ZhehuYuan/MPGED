#!/bin/bash


for ((j = 20; j <= 20; j += 5));
do
	for ((h = 128; h <= 128; h *= 2));
	do
		for ((k = 1; k <= 8; k*=2));
		do
			for i in {0..29}
			do
				./gedTune singlePair/AIDS${j}/$i $j $h $k>> ResultsTune/N${k}_out
			done
		done
	done
done
