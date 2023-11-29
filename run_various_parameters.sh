#! /usr/bin/env bash
# exit on fail
# set -e

mkdir -p workloadc

mkdir -p databases

curtime=`date +%m_%d_%H_%M`

nops=1000000
wl=workloadc

for index in WISCKEY ; do #  

  index_option=INDEX_${index}
  index_variant=0

  rm -rf build
  mkdir build
  cd build

  cmake .. -DCMAKE_BUILD_TYPE=Release -DNDEBUG_SWITCH=ON -D${index_option}=ON
  make -j 8 ycsb

  # osm_cellids books wiki_ts fb
  for data in wiki_ts; do
    echo "${index} write ${data}"
    until ./ycsb --print -w --dataset ${data}_10M_uint64 > ../${wl}/write_${index}_${index_variant}_${data}_${curtime}.txt
    do
      echo "Retrying write"
    done

    sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
    echo "${index} transaction / variant ${index_variant} / dataset ${data}"
    until ./ycsb --print --dataset ${data}_10M_uint64 --workload ${wl} -n ${nops} > ../${wl}/read_${index}_${index_variant}_${data}_${curtime}.txt
    do
      sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
      echo "Retrying write"
    done

  done

  cd ..

done


for index in RS CHT ; do #  

  index_option=INDEX_${index}

  for index_variant in 0 1 2 3 4 5 6 7 8 9; do
    rm -rf build
    mkdir build
    cd build

    cmake .. -DCMAKE_BUILD_TYPE=Release -DNDEBUG_SWITCH=ON -D${index_option}=ON -DINDEX_VARIANT=${index_variant}
    make -j 8 ycsb

    for data in wiki_ts; do
      echo "${index} write ${data}"
      until ./ycsb --print -w --dataset ${data}_10M_uint64 > ../${wl}/write_${index}_${index_variant}_${data}_${curtime}.txt
      do
        echo "Retrying write"
      done

      sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
      echo "${index} transaction / variant ${index_variant} / dataset ${data}"
      until ./ycsb --print --dataset ${data}_10M_uint64 --workload ${wl} -n ${nops} > ../${wl}/read_${index}_${index_variant}_${data}_${curtime}.txt
      do
        sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
        echo "Retrying write"
      done

    done

    cd ..

  done

done


for index in PGM ; do #  

  index_option=INDEX_${index}

  for index_variant in 4 8 16 32 64 256 1024 2048 4096 8192; do
    rm -rf build
    mkdir build
    cd build

    cmake .. -DCMAKE_BUILD_TYPE=Release -DNDEBUG_SWITCH=ON -D${index_option}=ON -DINDEX_VARIANT=${index_variant}
    make -j 8 ycsb

    for data in wiki_ts; do
      echo "${index} write ${data}"
      until ./ycsb --print -w --dataset ${data}_10M_uint64 > ../${wl}/write_${index}_${index_variant}_${data}_${curtime}.txt
      do
        echo "Retrying write"
      done

      sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
      echo "${index} transaction / variant ${index_variant} / dataset ${data}"
      until ./ycsb --print --dataset ${data}_10M_uint64 --workload ${wl} -n ${nops} > ../${wl}/read_${index}_${index_variant}_${data}_${curtime}.txt
      do
        sync; echo 3 | sudo tee /proc/sys/vm/drop_caches
        echo "Retrying write"
      done

    done

    cd ..

  done

done

