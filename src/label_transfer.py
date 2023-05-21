import scanpy as sc
import symphonypy as sp
import matplotlib.pyplot as plt


def build_reference(
    adata_ref, 
    obs_key='leiden', 
    target_group='Cycling', 
    layer='counts',
):
    if layer == "raw":
        adata_ref = adata_ref.raw.to_adata()
    else:
        adata_ref.X = adata_ref.layers[layer].copy()
    adata_ref = adata_ref.copy()
    adata_ref = adata_ref[adata_ref.obs[obs_key] != target_group]
    sc.pp.normalize_total(adata_ref, target_sum=1e5)
    sc.pp.log1p(adata_ref)
    sc.pp.highly_variable_genes(
        adata_ref,
        n_top_genes=3000,
    )
    adata_ref.raw = adata_ref
    adata_ref = adata_ref[:, adata_ref.var.highly_variable]
    sc.pp.scale(adata_ref, max_value=10)
    sc.pp.pca(adata_ref, n_comps=30, zero_center=False)
    return adata_ref


def get_query(
    adata, 
    obs_key='leiden', 
    target_group='Cycling',
    layer='counts',
):
    if layer == "raw":
        adata = adata.raw.to_adata()
    else:
        adata.X = adata.layers[layer].copy()
    adata_query = adata.copy()
    adata_query = adata_query[adata_query.obs[obs_key] == target_group]
    sc.pp.normalize_total(adata_query, target_sum=1e5)
    sc.pp.log1p(adata_query)
    return adata_query


def label_transfer(
    adata_query,
    adata_ref, 
    obs_key='leiden', 
    target_group='Cycling',
):
    sp.tl.map_embedding(
        adata_query=adata_query,
        adata_ref=adata_ref,
    )

    sp.tl.ingest(
        adata_query=adata_query,
        adata_ref=adata_ref,
        use_rep="X_pca_harmony",
    )

    sp.tl.transfer_labels_kNN(
        adata_query=adata_query,
        adata_ref=adata_ref,
        ref_labels="cell_type",
        ref_basis="X_pca",
        query_basis="X_pca_harmony",
    )
    fig, axes = plt.subplots(figsize=(8, 4), ncols=2)
    sc.pl.umap(
        adata_query,
        color="cell_type",
        frameon=False,
        title="Query dataset",
        ax=axes[0],
        show=False,
        legend_loc=None,
        size=120000 / adata_ref.shape[0],
    )

    sc.pl.umap(
        adata_ref,
        color="cell_type",
        frameon=False,
        title="Reference dataset",
        ax=axes[1],
        show=False,
    )
    return adata_query
