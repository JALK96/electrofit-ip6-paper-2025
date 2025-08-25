import json

# Load symmetry groups from a JSON file
def load_symmetry_groups(json_path):
    with open(json_path, "r") as file:
        symmetry_groups = json.load(file)
    return symmetry_groups