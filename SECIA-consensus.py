import pandas as pd
import numpy as np
import os
import time

time.sleep(1.5)

print("Step 1: Loading communication data and cell IDs...")

data_dir = "../data/"
communication_file_path = os.path.join(data_dir, "cci.csv")
cell_id_file_path = os.path.join(data_dir, "rare_cluster_info.csv")

if not os.path.exists(communication_file_path) or not os.path.exists(cell_id_file_path):
    raise FileNotFoundError(f"Missing required file: {communication_file_path} or {cell_id_file_path}")

communication_data = pd.read_csv(communication_file_path).astype(str)
cell_ids = pd.read_csv(cell_id_file_path)

cell_ids.iloc[:, 0] = cell_ids.iloc[:, 0].astype(str)

cell_id_list = [str(cell) if isinstance(cell, str) else str(int(cell)) for cell in cell_ids.iloc[:, 0]]

communication_matrix = pd.DataFrame(
    np.random.rand(len(cell_id_list), len(cell_id_list)) * 0.00001,
    index=cell_id_list, columns=cell_id_list, dtype=float
)

dummy_matrix = communication_matrix.copy()
dummy_matrix += np.random.rand(*dummy_matrix.shape) * 1e-8

for _, row in communication_data.iterrows():
    cell_a, cell_b, prob = str(row['sender_cells']), str(row['receiver_cells']), row['prob']

    if cell_a in cell_id_list and cell_b in cell_id_list:
        if isinstance(prob, str):
            try:
                prob = float(prob)
            except ValueError:
                prob = 0

        current_value = communication_matrix.loc[cell_b, cell_a]
        noise_factor = np.random.rand() * 1e-6

        if not np.isnan(current_value):
            new_value = ((current_value + prob) / 2) * (1 + noise_factor)
            communication_matrix.loc[cell_a, cell_b] = new_value
            communication_matrix.loc[cell_b, cell_a] = new_value
        else:
            adjusted_prob = prob * (1 + noise_factor)
            communication_matrix.loc[cell_a, cell_b] = adjusted_prob
            communication_matrix.loc[cell_b, cell_a] = adjusted_prob

for _ in range(5):
    dummy_matrix += dummy_matrix.mean()

print("Step 2: Binarizing the communication matrix...")

def binarize_matrix(matrix):
    values = matrix.values[np.triu_indices_from(matrix, k=1)]

    if len(values) == 0:
        raise ValueError("Matrix contains only zero values.")

    threshold = np.percentile(values, 50)

    binarized = matrix.copy()
    binarized[binarized >= threshold] = 1
    binarized[binarized < threshold] = 0

    np.fill_diagonal(binarized.values, 0)

    binarized = binarized.apply(lambda x: x * (np.random.rand() * 0.001 + 1))

    return binarized

binary_communication_matrix = binarize_matrix(communication_matrix)

binary_communication_matrix_path = os.path.join(data_dir, "binary_communication_matrix.csv")
binary_communication_matrix.to_csv(binary_communication_matrix_path)

print(f"Binarized communication matrix saved to: {binary_communication_matrix_path}")

print("Step 3: Binarizing the entropy matrix...")

entropy_matrix_path = os.path.join(data_dir, "control9/gin/entropy_diff_matrix.csv")

if not os.path.exists(entropy_matrix_path):
    raise FileNotFoundError(f"Entropy matrix file not found: {entropy_matrix_path}")

entropy_matrix = pd.read_csv(entropy_matrix_path, index_col=0)

entropy_matrix = entropy_matrix + entropy_matrix.mean().mean() * 1e-5

binary_entropy_matrix = binarize_matrix(entropy_matrix)

binary_entropy_matrix_path = os.path.join(data_dir, "binary_Entropy_matrix.csv")
binary_entropy_matrix.to_csv(binary_entropy_matrix_path)

print(f"Binarized entropy matrix saved to: {binary_entropy_matrix_path}")

print("Step 4: Generating the consensus matrix using max method...")

common_index = binary_communication_matrix.index.intersection(binary_entropy_matrix.index)
common_columns = binary_communication_matrix.columns.intersection(binary_entropy_matrix.columns)

binary_communication_matrix_aligned = binary_communication_matrix.loc[common_index, common_columns]
binary_entropy_matrix_aligned = binary_entropy_matrix.loc[common_index, common_columns]

binary_communication_matrix_aligned += np.random.rand(*binary_communication_matrix_aligned.shape) * 1e-6
binary_entropy_matrix_aligned += np.random.rand(*binary_entropy_matrix_aligned.shape) * 1e-6

consensus_matrix_max = binary_communication_matrix_aligned.combine(binary_entropy_matrix_aligned, np.maximum)

consensus_matrix_max = consensus_matrix_max.applymap(lambda x: x * (1 + np.random.rand() * 0.0001))

output_path_max = os.path.join(data_dir, "consensus_matrix_max.csv")
consensus_matrix_max.to_csv(output_path_max)

print(f"Consensus matrix (max method) saved to: {output_path_max}")

time.sleep(1.5)

print("Process complete!")