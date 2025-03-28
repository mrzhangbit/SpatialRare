import pandas as pd
import numpy as np
import os
import warnings
from tqdm import tqdm  # 用于进度条显示

# Ignore warnings
warnings.filterwarnings("ignore")

print("Step 1: Loading expression matrix...")
# Load gene expression matrix
data_path = "../data/control9/preprocessed_expression_matrix.csv"
expression_matrix = pd.read_csv(data_path, index_col=0)

# Check for missing values
if expression_matrix.isnull().values.any():
    print("Warning: Missing values detected. Imputing with column means.")
    expression_matrix = expression_matrix.apply(lambda x: x.fillna(x.mean()), axis=0)

print(f"Original expression matrix shape: {expression_matrix.shape}")

# Transpose matrix: rows as genes, columns as cells
expression_matrix = expression_matrix.T
print(f"Transposed expression matrix shape (gene x cell): {expression_matrix.shape}")

print("Step 2: Cleaning data...")
# Remove rows and columns that contain NaN values
expression_matrix = expression_matrix.apply(pd.to_numeric, errors='coerce')
expression_matrix = expression_matrix.dropna(axis=0, how='any').dropna(axis=1, how='any')

# Remove columns that are completely zero
expression_matrix = expression_matrix.loc[:, (expression_matrix != 0).any(axis=0)]
print(f"After cleaning, expression matrix shape: {expression_matrix.shape}")

print("Step 3: Calculating correlation matrix...")
expression_matrix_np = expression_matrix.values
print(f"Input to np.corrcoef dimensions (gene x cell): {expression_matrix_np.shape}")

# Compute gene-to-gene correlation
correlation_matrix = np.corrcoef(expression_matrix_np)
print("Correlation matrix calculation complete.")
print(f"Original GINE matrix dimensions (gene x gene): {correlation_matrix.shape}")

# Set correlation threshold
cor_threshold = 0.9
print(f"Correlation threshold set to: {cor_threshold}")

print("Step 4: Building neighbors dictionary...")

# Find top 100 nearest neighbors for each gene
top_k = 100
neighbors_dict = {
    i: np.argsort(-np.abs(correlation_matrix[i, :]))[:top_k]
    for i in range(correlation_matrix.shape[0])
}
print("Neighbors dictionary built. (Limited to top 100 neighbors per gene)")

print("Step 5: Calculating GINE matrix (gene x gene to cell x gene)...")

def calculate_cell_gene_gine(expression_matrix, correlation_matrix, neighbors_dict):
    expression_array = expression_matrix.to_numpy()
    n_genes, n_cells = expression_array.shape
    gine_matrix = np.zeros((n_cells, n_genes))

    print(f"Expression array shape: {expression_array.shape}")
    print(f"Correlation matrix shape: {correlation_matrix.shape}")

    for i in tqdm(range(n_cells), desc="Processing Cells"):  # Progress bar for better tracking
        for j in range(n_genes):
            neighbors = neighbors_dict[j]
            if len(neighbors) > 1:
                p_ik = correlation_matrix[j, neighbors] * expression_array[neighbors, i]
                p_ik_sum = p_ik.sum()

                if p_ik_sum != 0:
                    p_ik = p_ik / p_ik_sum
                    p_ik = np.clip(p_ik, 1e-10, None)
                    s_ij = -np.sum(p_ik * np.log(p_ik))
                    gine_matrix[i, j] = s_ij / np.log(len(neighbors))
                else:
                    gine_matrix[i, j] = 0
            else:
                gine_matrix[i, j] = 0

    gine_df = pd.DataFrame(gine_matrix, index=expression_matrix.columns, columns=expression_matrix.index)
    return gine_df

# Convert gene-to-gene GINE matrix into cell-to-gene matrix
gine_matrix = calculate_cell_gene_gine(expression_matrix, correlation_matrix, neighbors_dict)
print("Cell x Gene GINE matrix calculation complete.")
print(f"New GINE matrix dimensions (cell x gene): {gine_matrix.shape}")

print("Step 6: Calculating cell entropy values...")
cell_entropy = gine_matrix.mean(axis=1)
cell_entropy_df = pd.DataFrame({
    'Cell Name': gine_matrix.index,
    'Entropy': cell_entropy
})
print("Cell entropy values calculated.")

# Save entropy values
entropy_output_path = "../data/control9/gin/gine.csv"
os.makedirs(os.path.dirname(entropy_output_path), exist_ok=True)
cell_entropy_df.to_csv(entropy_output_path, index=False)
print(f"Cell entropy values saved to: {entropy_output_path}")

print("Step 7: Calculating entropy difference matrix...")
entropy_values = cell_entropy_df['Entropy'].values
cell_ids = cell_entropy_df['Cell Name'].values

# Compute entropy difference matrix
entropy_diff_matrix = np.abs(np.subtract.outer(entropy_values, entropy_values))
entropy_diff_df = pd.DataFrame(entropy_diff_matrix, index=cell_ids, columns=cell_ids)

# Save entropy difference matrix
entropy_diff_output_path = "../data/control9/gin/entropy_diff_matrix.csv"
entropy_diff_df.to_csv(entropy_diff_output_path)
print(f"Entropy difference matrix saved to: {entropy_diff_output_path}")

print("Process complete!")
