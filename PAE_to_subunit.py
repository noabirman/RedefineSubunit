import numpy as np
from sklearn.cluster import SpectralClustering
import json
import matplotlib.pyplot as plt



# Load the JSON file

with open("fold_mll4_1100_end_rbbp5_wdr5_p53x2/fold_mll4_1100_end_rbbp5_wdr5_p53x2_full_data_0.json", "r") as file:
    pae_data = json.load(file)
pae_matrix = np.array(pae_data['pae'])
#use pae as the adjacency matrix
symetrical_matrix = np.minimum(pae_matrix, pae_matrix.T)
# Perform spectral clustering
clustering = SpectralClustering(n_clusters=3, affinity='precomputed').fit(symetrical_matrix)
#show the clustering result as image where each pixel is a node and the color is the cluster
