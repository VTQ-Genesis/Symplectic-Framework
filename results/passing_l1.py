import json

with open("scan_results_nonempty_l1.json") as file:
    results = json.load(file)


iceberg = results["IcebergCode(6)"]

print(iceberg)


