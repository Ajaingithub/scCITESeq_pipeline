# Performed CellRanger on UCSF server 
# /diazlab/data3/.abhinav/.immune/Mayo/Tcells/analysis/preprocessing
# conda activate rapids_singlecell_gpu
# https://muon-tutorials.readthedocs.io/en/latest/cite-seq/1-CITE-seq-PBMC-5k.html
# https://docs.scvi-tools.org/en/1.0.0/tutorials/notebooks/cite_scrna_integration_w_totalVI.html
import os
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scvi
import rapids_singlecell as rsc
import scanpy as sc

datadir="/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/preprocessing/"
savedir="/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/"
Run01_basedir = datadir+"Run01_CR9/outs/per_sample_outs"
Run04_basedir = datadir+"Run04_CR9/outs/per_sample_outs"
os.chdir(savedir)

# List sample directories
Run01_samples = [os.path.join(Run01_basedir, d) for d in os.listdir(Run01_basedir) if os.path.isdir(os.path.join(Run01_basedir, d))]
Run04_samples = [os.path.join(Run04_basedir, d) for d in os.listdir(Run04_basedir) if os.path.isdir(os.path.join(Run04_basedir, d))]
Run01_Run04_sample = Run01_samples + Run04_samples
rnas = []
proteins = []
run_name= ["Run01"] * 10 + ["Run04"] * 10

for i in range(len(Run01_Run04_sample)):
    # extracting out the RNA
    h5_file = os.path.join(Run01_Run04_sample[i], "count", "sample_filtered_feature_bc_matrix.h5")
    adata = sc.read_10x_h5(h5_file, gex_only=False)
    adata.var_names_make_unique()
    rna = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
    sc.pp.filter_cells(rna, min_genes=100)
    sc.pp.filter_genes(rna, min_counts=3)
    rna.layers["counts"] = rna.X.copy()
    rna.obs_names = [f"{os.path.basename(Run01_Run04_sample[i])}_{run_name[i]}_{bc}" for bc in rna.obs_names]
    rna.obs["sample"] = os.path.basename(Run01_Run04_sample[i]) + "_" + run_name[i]
    rnas.append(rna)

### Trying only only first run
rnas1_adata =rnas[0:9]
keys = [rna.obs["sample"].unique()[0] for rna in rnas1_adata]
rna1_adata = sc.concat(rnas1_adata, label="sample", keys=keys)
rna1_adata.obs["run"] = rna1_adata.obs["sample"].str.rsplit("_", n=1).str[-1]

sc.pp.normalize_total(rna1_adata)
sc.pp.log1p(rna1_adata)
sc.pp.highly_variable_genes(rna1_adata)
sc.pp.scale(rna1_adata)
rsc.tl.pca(rna1_adata)
rsc.pp.neighbors(rna1_adata)
rsc.tl.umap(rna1_adata)

sc.pp.normalize_total(rna1_adata, target_sum=1e4)
sc.pp.log1p(rna1_adata)
sc.pp.highly_variable_genes(rna1_adata, flavor="seurat_v3", n_top_genes=2000)
rna1_adata = rna1_adata[:, rna1_adata.var.highly_variable].copy()
sc.pp.scale(rna1_adata)
rsc.tl.pca(rna1_adata)
rsc.pp.neighbors(rna1_adata, n_neighbors=15, n_pcs=50)
rsc.tl.umap(rna1_adata)

sc.pl.umap(
    rna1_adata,
    color=["run",'sample'],
    save = "_run_and_sample_unintegrated_Run01.pdf",
    legend_loc = "on data"
)

sc.pl.umap(
    rna1_adata,
    color=["run",'sample'],
    save = "_run_and_sample_unintegrated_Run01_2.pdf"
)


# Concatenate into a single AnnData object
keys = [rna.obs["sample"].unique()[0] for rna in rnas]
rna_adata = sc.concat(rnas, label="sample", keys=keys)
rna_adata.obs["run"] = rna_adata.obs["sample"].str.rsplit("_", n=1).str[-1]

sc.pp.normalize_total(rna_adata)
sc.pp.log1p(rna_adata)
sc.pp.highly_variable_genes(rna_adata)
rsc.tl.pca(rna_adata)
rsc.pp.neighbors(rna_adata)
rsc.tl.umap(rna_adata)

sc.pl.umap(
    rna_adata,
    color=["run",'sample'],
    save = "_run_and_sample_unintegrated_2.pdf",
    legend_loc = "on data"
)

scvi.model.SCVI.setup_anndata(
    rna_adata,
    layer="counts",
    batch_key = "run" # covariate
)

vae = scvi.model.SCVI(
    rna_adata,
    n_layers=3, 
    n_latent=30, ## Number of latent space
    encode_covariates=True,
    deeply_inject_covariates=True,  # inject tissue/treatment into all decoder layers
    use_layer_norm="both",
    use_batch_norm="none"
)

vae.train()
vae.save("rna_adata_vae_run_batches")
os.makedirs("saveh5ad", exist_ok = True)
rna_adata.write("saveh5ad/rna_adata_run01_04.h5ad")

SCVI_LATENT_KEY = "X_scVI"
latent = vae.get_latent_representation()
rna_adata.obsm[SCVI_LATENT_KEY] = latent

### GPU is super fast
rsc.tl.pca(rna_adata)
rsc.pp.neighbors(rna_adata,use_rep=SCVI_LATENT_KEY,n_neighbors=10)
rsc.tl.umap(rna_adata, min_dist=0.3)
sc.pl.umap(rna_adata,color=["run",'sample'],frameon=False, legend_loc='on data' ,save = "_run_sample_integrated.pdf")
sc.pl.umap(rna_adata,color=["run",'sample'],frameon=False ,save = "_run_sample_integrated_2.pdf")

scvi.model.SCVI.setup_anndata(
    rna_adata,
    layer="counts",
    batch_key = "sample" # covariate
)

vae = scvi.model.SCVI(
    rna_adata,
    n_layers=3, 
    n_latent=30, ## Number of latent space
    encode_covariates=True,
    deeply_inject_covariates=True,  # inject tissue/treatment into all decoder layers
    use_layer_norm="both",
    use_batch_norm="none"
)

vae.train()
vae.save("rna_adata_vae_run_samples")
os.makedirs("saveh5ad", exist_ok = True)
# rna_adata.write("saveh5ad/rna_adata_run01_04.h5ad")

SCVI_LATENT_KEY = "X_scVI_sample"
latent = vae.get_latent_representation()
rna_adata.obsm[SCVI_LATENT_KEY] = latent

### GPU is super fast
rsc.tl.pca(rna_adata)
rsc.pp.neighbors(rna_adata,use_rep=SCVI_LATENT_KEY,n_neighbors=10)
rsc.tl.umap(rna_adata, min_dist=0.3)
sc.pl.umap(rna_adata,color=["run",'sample'],frameon=False, legend_loc='on data' ,save = "_run_sampleintegrated.pdf")
sc.pl.umap(rna_adata,color=["run",'sample'],frameon=False, save = "_run_sampleintegrated_2.pdf")



# SCVI_CLUSTERS_KEY = "leiden_scVI_res1.4"
# rsc.tl.leiden(immune_adata, key_added=SCVI_CLUSTERS_KEY, resolution=1.4)
# sc.pl.umap(immune_adata,color=[SCVI_CLUSTERS_KEY],frameon=False, legend_loc='on data', save = "_leiden_clustering_rsc.pdf")
# immune_adata.write_h5ad("saveh5ad/immune_adata_8samples_clustering.h5ad")
