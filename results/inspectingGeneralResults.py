import json
import pprint
from qldpc import codes
from quantum_anticodes import choi
import pprint


# code = codes.SurfaceCode(1, 8, rotated=False)

# print(code.matrix_x)
# print()
# print(code.matrix_z)
# print(code.dimension)
# print(code.get_distance())

# Lx,Lz = choi(code.matrix_z.null_space(), code.matrix_x.null_space())
# print(Lx)
# print(Lz)



with open("scan_results.json") as file:
    results = json.load(file)

results_l1 = results["l1"]
results_l2 = results["l2"]

# # 2-> t, 1- > s
#pprint.pprint(f"s(l=1) result:{results_l1["SurfaceCode(1, 8, rotated=False)"]}")
#print()
#pprint.pprint(f"t(l=2) result:{results_l2["SurfaceCode(1, 8, rotated=False)"]}")

# flag = 0

# # for code_name, data in results_l2.items():
# #     for subset in data["passing_subsets"]:
# #         if len(subset) % 2 != 0:
# #             flag += 1
# #             print(code_name, subset)

# # # 0 if all even
# # print(flag)


# for code_name, data in results_l2.items():
#     if len(data["passing_subsets"]) > 10:
#         print(code_name, len(data["passing_subsets"]))

print(print(results_l2["SurfaceCode(2, 7, rotated=True)"]["passing_subsets"]))

