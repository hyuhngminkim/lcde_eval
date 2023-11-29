#! /usr/bin/env bash

# https://stackoverflow.com/questions/9713104/loop-over-tuples-in-bash
for index in WISCKEY BOURBON RS PGM CHT LCDE; do
  for i in osm_cellids_200M_uint64,3 wiki_ts_200M_uint64,9 fb_200M_uint64,2 books_200M_uint64,11; do
    IFS=","
    set -- $i
    # echo $1 and $2
    dataset=$1
    variant=$2
    echo $dataset
    echo $variant
  done
done