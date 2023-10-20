import sys
import argparse
import requests
import json

prefix = "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId="
suffix = ".zst"

input_data = sys.stdin.read().strip()

# metadata: the JSON metadata of the entire real-world SOSD dataset
# dataset: the name of the dataset being queried
metadata, dataset = input_data.split('|')
# filename: dataset.zst
filename = dataset + suffix

dataset_url = ""
persistentId = ""
md5 = ""

json_data = json.loads(metadata)

for entry in json_data["datasetVersion"]["files"]:
    if entry["label"] == filename:
        persistentId = entry["dataFile"]["persistentId"]
        md5 = entry["dataFile"]["md5"]

if persistentId == "" or md5 == "":
    sys.exit("Dataset is not valid or metadata is corrupted")

dataset_url = prefix + persistentId

# the dataset url and md5 will be printed to stdout to be captured by the 
# caller script
print(dataset_url)
print(md5)

