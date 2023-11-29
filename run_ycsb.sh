#! /usr/bin/env bash
# exit on fail
# set -e

mkdir -p workloada
mkdir -p workloadb
mkdir -p workloadc
mkdir -p workloadd
mkdir -p workloade
mkdir -p workloadf
mkdir -p databases
mkdir -p writes

curtime=`date +%m_%d_%H_%M`

nops=500000

# WISCKEY
# Uses no index, simple binary search
rm -rf build
mkdir build
cd build

index=WISCKEY
index_option=INDEX_${index}
cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON
make -j 8
# books_200M_uint64_shuffle fb_200M_uint64_shuffle wiki_ts_200M_uint64_shuffle osm_cellids_200M_uint64_shuffle
for dataset in  books_200M_uint64_shuffle fb_200M_uint64_shuffle wiki_ts_200M_uint64_shuffle osm_cellids_200M_uint64_shuffle; do
  echo "${index} write on ${dataset}"
  # ./ycsb --print -w --dataset ${dataset}
  until ./ycsb --print -w --dataset ${dataset} > ../writes/write_${dataset}_${index}_${curtime}.txt
  do
    echo "Retrying write"
  done

  for wl in workloade; do
    echo "${index} ${wl} on ${dataset}"
    ./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../${wl}/txn_${dataset}_${index}_${curtime}.txt
  done
done

cd ..

# RS
index=RS
index_option=INDEX_${index}
# books_200M_uint64_shuffle,7 fb_200M_uint64_shuffle,8 wiki_ts_200M_uint64_shuffle,5 osm_cellids_200M_uint64_shuffle,4
for i in books_200M_uint64_shuffle,7 fb_200M_uint64_shuffle,8 wiki_ts_200M_uint64_shuffle,5 osm_cellids_200M_uint64_shuffle,4; do
  IFS=","
  set -- $i
  dataset=$1
  variant=$2

  rm -rf build
  mkdir build
  cd build
  cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON -DINDEX_VARIANT=${variant}
  make -j 8
  echo "${index} write on ${dataset}"
  until ./ycsb --print -w --dataset ${dataset} > ../writes/write_${dataset}_${index}_${curtime}.txt
  do
    echo "Retrying write"
  done

  for wl in workloada workloadb workloadc workloadd workloade workloadf; do
    echo "${index} ${wl} on ${dataset}"
    ./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../${wl}/txn_${dataset}_${index}_${curtime}.txt
  done
  cd ..
done

# PGM
index=PGM
index_option=INDEX_${index}
# books_200M_uint64_shuffle,64 fb_200M_uint64_shuffle,32 osm_cellids_200M_uint64_shuffle,8
for i in  books_200M_uint64_shuffle,64 fb_200M_uint64_shuffle,32 osm_cellids_200M_uint64_shuffle,8; do
  IFS=","
  set -- $i
  dataset=$1
  variant=$2

  rm -rf build
  mkdir build
  cd build
  cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON -DINDEX_VARIANT=${variant}
  make -j 8
  echo "${index} write on ${dataset}"
  until ./ycsb --print -w --dataset ${dataset} > ../writes/write_${dataset}_${index}_${curtime}.txt
  do
    echo "Retrying write"
  done

  for wl in workloada workloadb workloadc workloadd workloade workloadf ; do
    echo "${index} ${wl} on ${dataset}"
    ./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../${wl}/txn_${dataset}_${index}_${curtime}.txt
  done
  cd ..
done

# CHT
index=CHT
index_option=INDEX_${index}
# books_200M_uint64_shuffle,1 fb_200M_uint64_shuffle,6 wiki_ts_200M_uint64_shuffle,2 osm_cellids_200M_uint64_shuffle,4
for i in books_200M_uint64_shuffle,1 fb_200M_uint64_shuffle,6 wiki_ts_200M_uint64_shuffle,2 osm_cellids_200M_uint64_shuffle,4; do
  IFS=","
  set -- $i
  dataset=$1
  variant=$2

  rm -rf build
  mkdir build
  cd build
  cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON -DINDEX_VARIANT=${variant}
  make -j 8
  echo "${index} write on ${dataset}"
  until ./ycsb --print -w --dataset ${dataset} > ../writes/write_${dataset}_${index}_${curtime}.txt
  do
    echo "Retrying write"
  done

  for wl in workloada workloadb workloadc workloadd workloade workloadf; do
    echo "${index} ${wl} on ${dataset}"
    ./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../${wl}/txn_${dataset}_${index}_${curtime}.txt
  done
  cd ..
done

# BOURBON
index=BOURBON
index_option=INDEX_${index}

rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON
make -j 8
# books_200M_uint64_shuffle fb_200M_uint64_shuffle wiki_ts_200M_uint64_shuffle 
for dataset in books_200M_uint64_shuffle fb_200M_uint64_shuffle wiki_ts_200M_uint64_shuffle osm_cellids_200M_uint64_shuffle; do
  echo "${index} write on ${dataset}"
  until ./ycsb --print -w --dataset ${dataset} > ../writes/write_${dataset}_${index}_${curtime}.txt
  do
    echo "Retrying write"
  done

  for wl in workloada workloadb workloadc workloadd workloade workloadf; do
    echo "${index} ${wl} on ${dataset}"
    ./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../${wl}/txn_${dataset}_${index}_${curtime}.txt
  done
done

cd ..

# LCDE
index=LCDE
index_option=INDEX_${index}

rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -D${index_option}=ON
make -j 8
# books_200M_uint64_shuffle fb_200M_uint64_shuffle wiki_ts_200M_uint64_shuffle osm_cellids_200M_uint64_shuffle
for dataset in books_200M_uint64_shuffle fb_200M_uint64_shuffle wiki_ts_200M_uint64_shuffle osm_cellids_200M_uint64_shuffle; do
  echo "${index} write on ${dataset}"
  until ./ycsb --print -w --dataset ${dataset} > ../writes/write_${dataset}_${index}_${curtime}.txt
  do
    echo "Retrying write"
  done

  for wl in workloada workloadb workloadc workloadd workloade workloadf; do
    echo "${index} ${wl} on ${dataset}"
    ./ycsb --print --dataset ${dataset} --workload ${wl} -n ${nops} > ../${wl}/txn_${dataset}_${index}_${curtime}.txt
  done
done

cd ..