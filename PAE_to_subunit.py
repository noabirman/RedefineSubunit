import numpy as np
import json


# Load the JSON file
with open("pae.json", "r") as file:
    pae_data = json.load(file)

# Extract the matrix from the JSON
pae_matrix = np.array(pae_data[0]['predicted_aligned_error'], dtype=np.float64)

# Print or use the reconstructed matrix
print(pae_matrix)