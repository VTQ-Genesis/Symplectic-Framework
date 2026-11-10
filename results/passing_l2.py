import json
from quantum_anticodes import choi



with open("scan_results_nonempty_l2.json") as file:
    results = json.load(file)




iceberg2 = results["IcebergCode(2)"]

print(iceberg2)