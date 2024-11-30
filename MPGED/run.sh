#!/bin/bash


for ((j = 15; j <= 30; j += 5));
do
	for k in {0..3}
	do
		for i in {0..29}
		do
			./ged singlePair/AIDS${j}/$i $j 1 $k>> Results/AIDS_len${j}_thread1_method${k}_out
		done
	done

	for ((h = 4; h <= 4; h *= 2));
	do
		for k in {1..3}
		do
			for i in {0..29}
			do
				./ged singlePair/AIDS${j}/$i $j $h $k>> Results/AIDS_len${j}_thread${h}_method${k}_out
			done
		done
	done
done

for ((j = 15; j <= 30; j += 5));
do
	for k in {0..3}
	do
		for i in {0..29}
		do
			./ged singlePair/PubChem${j}/$i $j 1 $k>> Results/PubChem_len${j}_thread1_method${k}_out
		done
	done

	for ((h = 4; h <= 4; h *= 2));
	do
		for k in {1..3}
		do
			for i in {0..29}
			do
				./ged singlePair/PubChem${j}/$i $j $h $k>> Results/PubChem_len${j}_thread${h}_method${k}_out
			done
		done
	done
done
