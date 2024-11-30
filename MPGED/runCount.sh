#!/bin/bash


for ((j = 20; j <= 20; j += 5));
do
	for ((h = 1; h <= 128; h *= 2));
	do
		for i in {0..29}
		do
			for k in {1..3}
			do
				./gedCount singlePair/AIDS${j}/$i $j $h $k>> ResultsCount/AIDS_len${j}_thread${h}_method${k}_out
			done
		done
	done
done
