#! /usr/bin/env bash
# exit on fail
# set -e

mkdir -p workloada
mkdir -p workloadb
mkdir -p workloadc
mkdir -p workloadd
mkdir -p workloade
mkdir -p workloadf

curtime=`date +%m_%d_%H_%M`

dataset=osm_cellids_10M_uint64

nops=100

for index in  WISCKEY BOURBON LCDE; do
  rm -rf build
  mkdir build
  cd build

  index_option=INDEX_${index}
  cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON
  make -j 8 ycsb

  for wl in workloadc workloadf ; do
    sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
    echo "${index} write"
    until ./ycsb --print -w --dataset ${dataset} > ../${wl}/write_${index}_${wl}_${curtime}.txt
    do
      echo "Retrying write"
    done

    sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
    echo "${index} transaction : ${wl}"
    until ./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../${wl}/txn_${index}_${wl}_${curtime}.txt
    do
      echo "Retrying transaction"
      sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
    done
  done

  cd ..
done

