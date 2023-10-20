#! /usr/bin/env bash
# exit on fail
set -e

mkdir -p results
mkdir -p write_results

rm -rf build
mkdir build
cd build

curtime=`date +%m_%d_%H_%M`

dataset=osm_cellids_10M_uint64

key_size=16
val_size=64

nops=10000

index=CHT

# Main execution code
index_option=INDEX_${index}
cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON
make -j 8 ycsb

for wl in workloada workloadb workloadc workloadd workloade workloadf; do
  sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
  echo "${index} write"
  until ./ycsb --print -w --dataset ${dataset} > ../write_results/write_${index}_${wl}_${curtime}.txt
  do
    echo "Retrying write"
  done

  sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
  echo "${index} transaction : ${wl}"
  until ./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../results/txn_${index}_${wl}_${curtime}.txt
  do
    echo "Retrying transaction"
  done
done