import json

# === File paths ===
reaction_list_path = "Data/QtFiles/Qt_1/qt_reactions_list.txt"
thresholds_config_path = "thresholds_config.json"
output_path = "reaction_groups.json"

# === Load full Qt reaction list ===
with open(reaction_list_path, 'r') as f:
    reaction_list = [
        line.strip().split(maxsplit=1)[1].replace(" ", "").lower()
        for line in f if line.strip()
    ]

# === Load group names ===
with open(thresholds_config_path, 'r') as f:
    group_names = list(json.load(f).keys())

# === === MANUAL RANGES TO FILL IN === ===
# Each entry is: (group_index, start_index, end_index) — inclusive range
manual_group_ranges = [
    (0, 0, 26),    
    (1, 27, 90),
    (2,91,112),
    (3,113,114),
    (4,115,122),
    (5,123,140),
    (6,141,152),
    (7,153,171),
    (8,172,190),
    (9,191,198),
    (10,199,233),
    (11,234,245),
    (12,246,262),
    (13,264,286),
    (14,287,317),
    (15,318,355),
    (16,356,385),
    (17,386,446),
    (18,447,479),
    (19,480,620),
    (20,621,655)

]

# === Generate reaction_groups dictionary ===
reaction_groups = {}

for group_index, start, end in manual_group_ranges:
    group_name = group_names[group_index]
    reactions = reaction_list[start:end + 1]
    reaction_groups[group_name] = reactions

# === Save output ===
with open(output_path, 'w') as f_out:
    json.dump(reaction_groups, f_out, indent=2)

print("reaction_groups.json has been created")

