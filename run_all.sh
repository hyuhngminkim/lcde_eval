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

nops=10000000

for index in  WISCKEY BOURBON LCDE; do #WISCKEY RS PGM CHT LCDE BOURBON
  rm -rf build
  mkdir build
  cd build

  index_option=INDEX_${index}
  cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON
  make -j 8 ycsb

  for dataset in osm_cellids_200M_uint64 wiki_ts_200M_uint64 fb_200M_uint64 books_200M_uint64; do
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
  done
  cd ..
done