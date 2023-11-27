import pandas as pd
import seaborn as sns
import numpy as np
import sys

sys.path.insert(0, r"C:\Users\ICN_guest\Desktop\ui")

import gaussianize

def gaussianize_data(data, column_name):
    out = gaussianize.Gaussianize(strategy="brute", max_iter=200, verbose=True)
    out.fit(data[column_name])
    data[f"{column_name}_gaussian"] = out.transform(data[column_name])
    return data

# select files from OneDrive
fpath_fmap = r'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\results_clean\Neurotransmitter_analysis\functional_connectivity\beta_spmT_0001_compound_atlas_HCPex_SUIT_ABGT.csv'
fpath_fdopa = r'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\Analyses\Parcellated_Correlations\Compound_Atlas\FDOPA_fluorodopa_hc12_gomez_compound_atlas_HCPex_SUIT_ABGT.csv'
fpath_fdopa = r'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\Analyses\Parcellated_Correlations\Compound_Atlas\D1_SCH23390_hc13_kaller_compound_atlas_HCPex_SUIT_ABGT.csv'
fpath_fdopa = r'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\PET Data\aggregated\Dopamine_aggregate_map_compound_atlas_HCPex_SUIT_ABGT.csv'

# read table
fmap = pd.read_csv(fpath_fmap, index_col="Name").rename(columns={"Value": "fmap"}).drop(columns="Index")
fdopa = pd.read_csv(fpath_fdopa, index_col="Name").rename(columns={"Value": "fdopa"}).drop(columns=["Index", "structure"])

# preprocess data
data = pd.concat((fmap, fdopa), axis="columns", join="outer")
data = data[data["fmap"] != 0]
data = data[data["fdopa"] != 0]
data = gaussianize_data(data=data, column_name="fmap")
data = gaussianize_data(data=data, column_name="fdopa")


#sns.set_context("poster")
#sns.set_palette("colorblind")
#sns.scatterplot(data=data, x="fmap_gaussian", y="fdopa_gaussian", hue="structure", palette='BuPu')

#sns.jointplot(data=data, x="fmap_gaussian", y="fdopa_gaussian", hue='structure', palette='PuBu')
#sns.jointplot(data=data, x="fmap_gaussian", y="fdopa_gaussian", hue="structure", palette='magma')
# Visualise with seaborn
sns.lmplot(data=data, x="fmap_gaussian", y="fdopa_gaussian", hue="structure", palette='magma')