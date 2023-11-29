#! /usr/bin/env bash

rm -rf results
rm -rf write_results
rm -rf read_results

mkdir -p ycsb_benchmark_results
for workload_dir in workloada workloadb workloadc workloadd workloade workloadf; do
  if [[ -d "${workload_dir}" && -d "ycsb_benchmark_results/${workload_dir}" ]]; then
    mv "${workload_dir}"/* "ycsb_benchmark_results/${workload_dir}/"
  else
    mv ${workload_dir} ycsb_benchmark_results
  fi
  rm -rf ${workload_dir}
done
# mv workloada workloadb workloadc workloadd workloade workloadf ycsb_benchmark_results