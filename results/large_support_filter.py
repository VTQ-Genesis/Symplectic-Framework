import json


with open("scan_results_nonempty_l1.json") as file:
    results_l1 = json.load(file)
with open("scan_results_nonempty_l2.json") as file:
    results_l2 = json.load(file)


def has_large_support(result):
    length = result["n"] + result["k"]  # Choi code length, matching J's coordinates.
    min_size = max(4, (length + 1) // 2)
    return any(len(J) >= min_size for J in result["passing_subsets"])


# Filter existing results: at least one passing J must satisfy both size limits.
filtered_l1 = {
    name: result
    for name, result in results_l1.items()
    if has_large_support(result)
}

filtered_l2 = {
    name: result
    for name, result in results_l2.items()
    if has_large_support(result)
}

with open("l1_large_supports.json", "w") as file:
    json.dump(filtered_l1, file, indent=2)

with open("l2_large_supports.json", "w") as file:
    json.dump(filtered_l2, file, indent=2)X

print("l=1:", list(filtered_l1))
print("l=2:", list(filtered_l2))