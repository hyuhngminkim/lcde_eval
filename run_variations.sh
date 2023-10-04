#! /usr/bin/env bash
# exit on fail
# set -e

mkdir -p results
mkdir -p write_results

curtime=`date +%m_%d_%H_%M`

dataset=osm_cellids_10M_uint64
# select from 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
index_variant=0
# select from 4, 8, 16, 32, 64, 256, 1024, 2048, 4096, 8192
pgm_index_variant=4

wl=workloadc
nops=100000

# WISCKEY
rm -rf build
mkdir build
cd build
index=WISCKEY
index_option=INDEX_${index}
cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON -DINDEX_VARIANT=${index_variant}
make -j 8 ycsb
sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
echo "${index} transaction : ${wl}"
./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../results/txn_${index}_${index_variant}_${wl}_${curtime}.txt
cd ..

sync; echo 3 | sudo tee /proc/sys/vm/drop_caches

# All others
for index in BOURBON RS CHT LINEAR; do
  rm -rf build
  mkdir build
  cd build

  index_option=INDEX_${index}
  cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON -DINDEX_VARIANT=${index_variant}
  make -j 8 ycsb

  sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
  echo "${index} transaction : ${wl}"
  ./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../results/txn_${index}_${index_variant}_${wl}_${curtime}.txt
  cd ..
done

sync; echo 3 | sudo tee /proc/sys/vm/drop_caches

# PGM
rm -rf build
mkdir build
cd build
index=PGM
index_option=INDEX_${index}
cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON -DINDEX_VARIANT=${pgm_index_variant}
make -j 8 ycsb
echo "${index} transaction : ${wl}"
./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../results/txn_${index}_${index_variant}_${wl}_${curtime}.txt

sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
