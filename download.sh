#! /usr/bin/env bash

echo "downloading metadata ..."
JSON=$(curl -s 'https://dataverse.harvard.edu/api/datasets/export?exporter=dataverse_json&persistentId=doi%3A10.7910/DVN/JGVF9A')
echo "done"

function get_checksum() {
  FILE=$1

  if [ -x "$(command -v md5sum)" ]; then
    MD5_RESULT=`md5sum ${FILE} | awk '{print $1}'`
  fi
}

function main() {
  echo "downloading data ..."
  mkdir -p data
  cd data
  suffix="_200M_uint64"
  for data in wiki_ts books osm_cellids fb; do
    DATASET="${data}${suffix}"
    FILE="${DATASET}.zst"
    RES=$(echo "${JSON}|${DATASET}" | python3 ../parse.py)
    ARR=(${RES})
    URL=${ARR[0]}
    CHECKSUM=${ARR[1]}

    if [ -f ${FILE} ]; then
      echo "file already exists"
    else
      wget ${URL} -O "${FILE}"
      get_checksum "${FILE}"
      if [ "${MD5_RESULT}" != "${CHECKSUM}" ]; then
        echo "error checksum does not match: run download again"
        exit -1
      else 
        echo ${FILE} "checksum ok"
        zstd -d "${FILE}"
        rm -f ${FILE}
      fi
      
    fi
  done
  cd ..
  echo "done"
}

main