import os
import re
import csv

# List of folder names and their abbreviations
folders = {
    "LIFTOVER_3":           "3",
    "LIFTOVER_3_GBWT": "3_GBWT",
    "LIFTOVER_6":           "6",
    "LIFTOVER_6_GBWT": "6_GBWT",
    "LIFTOVER_ALL":       "ALL",
    "LIFTOVER_ALL_GBWT": "ALL_GBWT"
}

metadata_file = "metadata2.csv"
output_file = "get_liftover.csv"

# Initialize sample list
samples = []
with open(metadata_file, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    next(reader)  # Skip header if there is one
    for row in reader:
        samples.append((row[0], row[1]))

# Regex pattern to extract percentage not lifted over
pattern = re.compile(
    r"INFO\s+\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2}\s+LiftoverVcf\s+([\d.]+)% of variants were not successfully lifted over and written to the output"
)

# Function to extract the percentage value from a file
def extract_percentage(file_path):
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            for line in file:
                match = pattern.search(line)
                if match:
                    return float(match.group(1))
    return 0.0  # Return 0 if no match is found

# Header for the CSV output
header = ["Sample", "Library_Type"] + [f"{abbr}_{source}" for abbr in folders.values() for source in ["BCF", "VG", "FB"]]

# Write the header and results to the output CSV file
data_rows = []
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)

    for sample, library_type in samples:
        row = [sample, library_type]
        for folder, abbr in folders.items():
            file_path_bcf = os.path.join(folder, f"{sample}.liftover.bcf.txt")
            file_path_vg = os.path.join(folder, f"{sample}.liftover.vg.txt")
            file_path_fb = os.path.join(folder, f"{sample}.liftover.fb.txt")

            # Append extracted percentages to the row
            row.append(extract_percentage(file_path_bcf))
            row.append(extract_percentage(file_path_vg))
            row.append(extract_percentage(file_path_fb))
            print(f"Processed {folder} for sample {sample}")

        writer.writerow(row)
        data_rows.append(row)  # Collect rows of data

print(f"\nOutput successfully written to {output_file}")
