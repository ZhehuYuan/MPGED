# MP-GED

This is the code repository for MP-GED

---

## Download the code repository for X-TED

Firstly, use the command line below to download this repository.

```
$ git clone https://github.com/ZhehuYuan/MPGED.git
$ cd MPGED/MPGED/
```

---

## Dataset

We use two real-world application datasets, namely AIDS (https://dtp.cancer.gov/docs/aids/aidsdata.html) and PubChem (https://pubchem.ncbi.nlm.nih.gov/). For each dataset, we randomly select four groups of graph pairs, with each group containing 30 pairs of graphs having the same number of vertices. The graphs in these groups contain 15, 20, 25, and 30 vertices, respectively. The graph pairs are stored in "singlePair/".  We briefly introduce the datasets as follows:

  • AIDS is an antivirus screening dataset containing 42,687 graphs. Each graph in this dataset has an average of 25.6 vertices and 27.6 edges.

  • PubChem is a chemical compound dataset. Due to its ex-tremely large size, we selected the first 21,160 graphs from the dataset, with each graph having an average of 26.8 vertices and 26.6 edges.

---

## Test the MP-GED

Running theses tests, we use a machine with two AMD EPYC 7713 64-Core thread, 250GB of main memory, and the 64-bit Rocky Linux operating system.

1) Complie code

```
$ make
```

2) Tune the heap partitioning parameter (𝑁). The output files will be in "ResultsTune/". The first column is the number of atomicCAS operations, and the second column is the number of heaps searched.

```
$ ./runTune
```

3) Test the GEDs using one and four threads on four groups of graph pairs in each dataset, while the graph in each group contains 15, 20, 25, 30 vertices, respectively. The output files of running time will be in "Results/". 

```
$ ./run
```

4) Test the GEDs using 1 to 128 threads on 30 pairs of graphs in AIDS dataset with 20 vertices. The output files of running time will be in "ResultsMassive/".

```
$ ./runMassive
```

5) Test the GEDs and analyze where does time go using 1 to 128 threads on 30 pairs of graphs in AIDS dataset with 20 vertices. The output files will be in "ResultsTimeGo/". The first and second column of HGED and PGED output files represent the pre-processing and execution time, respectively. The first to fifth column of MP-GED output files represent the read, read search, write, write search, and execution time, respectively.

```
$ ./runTimeGo
```

6) Test the GEDs and count the number of nodes read and nodes explored by each method on 30 pairs of graphs in AIDS dataset with 20 vertices. The output files will be in "ResultsCount/". The first column is the number of nodes read, and the second column is the number of nodes explored.

```
$ ./runCount
```

7) Test the MP-GED and observe the changes of upper bound against time on 10 pairs of graphs in AIDS dataset with 20 vertices. The output files will be in "ResultsUb/". There will be 10 files, each represents the upper bound (the second column) and time (the first column) of one pair of graphs.

```
$./runUb
```


Excepting step 7, all other output files are name as \<dataset name\>_len\<number of vertices in each graph\>_thread\<number of threads\>_method\<method Id, 0 for A-star+, 1 for PGED, 2 for HGED, 3 for MPGED, 4 for CPQ variance of MP-GED\>_out.

## Run the MP-GED with your data and settings

```
./ged <file name> <number of nodes in each graph> <number of threads> <method Id, 0 for A-star+, 1 for PGED, 2 for HGED, 3 for MPGED, 4 for CPQ variance of MP-GED>
```
