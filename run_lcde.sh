#! /usr/bin/env bash
# exit on fail
set -e

mkdir -p workloada
mkdir -p workloadb
mkdir -p workloadc
mkdir -p workloadd
mkdir -p workloade
mkdir -p workloadf
mkdir -p databases
mkdir -p writes

curtime=`date +%m_%d_%H_%M`

nops=10000000

# LCDE
index=LCDE
index_option=INDEX_${index}

rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON -DINTERNAL_TIMER_SWITCH=ON
make -j 8
# books_200M_uint64_shuffle fb_200M_uint64_shuffle wiki_ts_200M_uint64_shuffle osm_cellids_200M_uint64_shuffle
for dataset in fb_200M_uint64; do
  echo "${index} write on ${dataset}"
  # ./ycsb --print -w --dataset ${dataset} > ../writes/write_${index}_${dataset}_${curtime}.txt


  for wl in workloadc; do #workloadc workloada workloadb workloadd workloade workloadf ; do
    sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
    echo "${index} ${wl} on ${dataset}"
    until ./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../${wl}/txn_${index}_${dataset}_${curtime}.txt
    do
      echo "Retrying transaction"
      sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
    done
  done
done

cd ..