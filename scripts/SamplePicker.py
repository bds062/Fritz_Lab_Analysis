import random

# Create a dictionary to store values from the fourth column and their corresponding SRR numbers
column4_dict = {}

# Open and read the SRRINFO.txt file
with open("SRRINFO.txt", "r") as file:
    for line in file:
        # Split the line by tabs and extract the needed columns
        parts = line.strip().split("\t")
        
        if len(parts) < 4:
            continue  # Skip malformed lines
        
        # Assign columns: SRR number (3rd column), Fourth column value (4th column)
        srr_number = parts[2]
        column4_value = parts[3]
        
        # Group the SRR numbers by the fourth column value
        if column4_value in column4_dict:
            column4_dict[column4_value].append(srr_number)
        else:
            column4_dict[column4_value] = [srr_number]

# Now, randomly pick 3 SRR numbers from each fourth column value group
selected_srrs = {}

for column4_value, srr_numbers in column4_dict.items():
    if len(srr_numbers) >= 3:
        # If there are at least 3 SRR numbers, randomly choose 3
        selected_srrs[column4_value] = random.sample(srr_numbers, 3)
    else:
        # If fewer than 3 SRR numbers exist, take all of them
        selected_srrs[column4_value] = srr_numbers

# Print or process the selected SRR numbers from each fourth column value group
for column4_value, srr_numbers in selected_srrs.items():
    print(f"Location: {column4_value} -> Selected SRR Numbers: {' '.join(srr_numbers)}")