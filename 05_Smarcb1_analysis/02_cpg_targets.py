# %%
import pandas as pd
pd.option_context("display.max_rows", None, "display.max_columns", None)
# %%
def load_cpg_enrichment_data(file_path):
    """
    Loads the CpG enrichment data from the specified CSV file into a pandas DataFrame.

    Args:
        file_path (str): The path to the CSV file.

    Returns:
        pandas.DataFrame: A DataFrame containing the CpG enrichment data.
    """
    try:
        df = pd.read_csv(file_path)
        return df
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# %%
# Set base_dir to your CpG enrichment results directory
base_dir = ""  # e.g., "/path/to/cpg_enrichment"
file_path = f"{base_dir}/Neu/broad/cpg_enrichment_1_rep_in_peaks/cpg_enrichment_parallel.csv"
cpg_data = load_cpg_enrichment_data(file_path)
print(cpg_data.shape)
cpg_data.head()

# %%
cpg_data = cpg_data[["chr", "start", "end", "binding_type", "binding_type_by_peaks", "exo_replicates_with_peaks", "endo_replicates_with_peaks"]]
print(cpg_data.binding_type.unique())
print(cpg_data.binding_type_by_peaks.unique())

# %%
import matplotlib.pyplot as plt

plt.hist(cpg_data["exo_replicates_with_peaks"], bins=20, alpha=0.5, label="Exo Replicates")
plt.hist(cpg_data["endo_replicates_with_peaks"], bins=20, alpha=0.5, label="Endo Replicates")
plt.xlabel("Number of Replicates with Peaks")
plt.ylabel("Frequency")
plt.title("Distribution of Replicates with Peaks")
plt.legend(loc='upper right')
plt.show()

# %%
# Set cpgs_path to your CpG islands BED file
cpgs_path = ""  # e.g., "/path/to/data/cpg_islands.bed"
cpgs = pd.read_csv(cpgs_path, sep="\t", header=None)
print(cpgs.shape)
cpgs.head()

# %%
cpgs.columns = ["chr", "start", "end", "id", "cpg", "number"]
cpgs = cpgs[["chr", "start", "end"]]
cpgs

# %%
cpgs_targeted = pd.merge(cpgs, cpg_data[["chr", "start", "end"]], on=["chr", "start", "end"], how="inner")
cpgs_nottargeted = pd.merge(cpgs, cpgs_targeted[["chr", "start", "end"]], on=["chr", "start", "end"], how="left", indicator=True)
cpgs_nottargeted = cpgs_nottargeted[cpgs_nottargeted['_merge'] == 'left_only'].drop('_merge', axis=1)
cpgs_nottargeted = cpgs_nottargeted[["chr", "start", "end"]]

print(f"Shape of targeted cpgs: {cpgs_targeted.shape}")
print(f"Shape of not targeted cpgs: {cpgs_nottargeted.shape}")

# %%
cpgs_targeted.head()
# Set output_dir to your output directory
output_dir = ""  # e.g., "/path/to/output"
cpgs_targeted.to_csv(f"{output_dir}/neu_cpg_targets.bed", sep="\t", header=False, index=False)

# %%
cpgs_nottargeted.head()
cpgs_nottargeted.to_csv(f"{output_dir}/neu_cpg_not_targets.bed", sep="\t", header=False, index=False)

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

# %%
file_path = f"{base_dir}/NSC/broad/cpg_enrichment_1_rep_in_peaks/cpg_enrichment_parallel.csv"
cpg_data = load_cpg_enrichment_data(file_path)
print(cpg_data.shape)
cpg_data.head()

# %%
cpg_data = cpg_data[["chr", "start", "end", "binding_type", "binding_type_by_peaks", "exo_replicates_with_peaks", "endo_replicates_with_peaks"]]
print(cpg_data.binding_type.unique())
print(cpg_data.binding_type_by_peaks.unique())

# %%
import matplotlib.pyplot as plt

plt.hist(cpg_data["exo_replicates_with_peaks"], bins=20, alpha=0.5, label="Exo Replicates")
plt.hist(cpg_data["endo_replicates_with_peaks"], bins=20, alpha=0.5, label="Endo Replicates")
plt.xlabel("Number of Replicates with Peaks")
plt.ylabel("Frequency")
plt.title("Distribution of Replicates with Peaks")
plt.legend(loc='upper right')
plt.show()

# %%
cpgs = pd.read_csv(cpgs_path, sep="\t", header=None)
print(cpgs.shape)
cpgs.head()

# %%
cpgs.columns = ["chr", "start", "end", "id", "cpg", "number"]
cpgs = cpgs[["chr", "start", "end"]]
cpgs

# %%
cpgs_targeted = pd.merge(cpgs, cpg_data[["chr", "start", "end"]], on=["chr", "start", "end"], how="inner")
cpgs_nottargeted = pd.merge(cpgs, cpgs_targeted[["chr", "start", "end"]], on=["chr", "start", "end"], how="left", indicator=True)
cpgs_nottargeted = cpgs_nottargeted[cpgs_nottargeted['_merge'] == 'left_only'].drop('_merge', axis=1)
cpgs_nottargeted = cpgs_nottargeted[["chr", "start", "end"]]

print(f"Shape of targeted cpgs: {cpgs_targeted.shape}")
print(f"Shape of not targeted cpgs: {cpgs_nottargeted.shape}")

# %%
cpgs_targeted.head()
cpgs_targeted.to_csv(f"{output_dir}/nsc_cpg_targets.bed", sep="\t", header=False, index=False)

# %%
cpgs_nottargeted.head()
cpgs_nottargeted.to_csv(f"{output_dir}/nsc_cpg_not_targets.bed", sep="\t", header=False, index=False)
