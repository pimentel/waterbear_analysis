#!/bin/bash

module load python/3.6.1
snakemake -j 999 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -n {cluster.n} -t {cluster.time} --mem {cluster.memory}"
