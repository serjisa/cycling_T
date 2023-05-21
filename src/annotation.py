import scanpy as sc
from sklearn.neighbors import KNeighborsClassifier


def classify_cells(
    adata: sc.AnnData,
    obs_key: str,
    target_group: str,
    n_neighbors: int = 20,
    use_rep: str = "X_pca",
    key_added: str = "labels_pred",
    use_target: bool = True,
    **kwargs,
):
    knn = KNeighborsClassifier(n_neighbors=n_neighbors, **kwargs)
    if use_target:
        knn.fit(adata.obsm[use_rep], adata.obs[obs_key])
    else:
        knn.fit(
            adata[adata.obs[obs_key] != target_group].obsm[use_rep],
            adata[adata.obs[obs_key] != target_group].obs[obs_key],
        )
    
    labels_pred = knn.predict(adata[adata.obs[obs_key] == target_group].obsm[use_rep])
    adata.obs[key_added] = adata.obs[obs_key].copy()
    adata.obs[key_added][adata[adata.obs[obs_key] == target_group].obs_names] = labels_pred
    
    
def log1p(
    adata: sc.AnnData,
    n_top_genes: int = 3000,
    target_sum: float = 1e4,
    do_scaling: bool = True,
    batch_key: str | None = None,
    n_pcs: int = 30,
) -> sc.AnnData:
    adata = adata.copy()
    
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    if batch_key is None:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    else:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
    adata.raw = adata
    
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=n_pcs)
    
    if batch_key is None:
        use_rep = "X_pca"
    else:
        use_rep = "X_pca_harmony"
        sp.pp.harmony_integrate(adata, key=batch_key, max_iter_harmony=50)
    sc.pp.neighbors(adata, use_rep=use_rep)
    sc.tl.umap(adata)
    
    return adata


def annotation_with_removed_genes(
    adata: sc.AnnData,
    genes: list,
    obs_key: str,
    target_group: str,
    batch_key: str | None = None,
    layer: str = "X",
    n_pcs=30,
    n_neighbors=20,
) -> sc.AnnData:
    adata = adata.copy()

    if layer == "raw":
        adata = adata.raw.to_adata()
    elif layer != "X":
        adata.X = adata.layers[layer].copy()
        
    adata = adata[:, ~adata.var_names.isin(genes)]
    adata = log1p(adata, batch_key=batch_key, n_pcs=n_pcs)
    if batch_key is None:
        use_rep = "X_pca"
    else:
        use_rep = "X_pca_harmony"
    classify_cells(
        adata=adata,
        obs_key=obs_key,
        target_group=target_group,
        n_neighbors=n_neighbors,
    )
    return adata


def annotate_with_cc_scoring(
    adata, 
    obs_key, 
    target_group, 
    layer='counts', 
    batch_key=None, 
    n_pcs=15, 
    n_neighbors=10,
):
    adata = adata.copy()
    if layer == "raw":
        adata = adata.raw.to_adata()
    elif layer != "X":
        adata.X = adata.layers[layer].copy() 
        
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    adata.raw = adata
    sc.pp.scale(adata)
    sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, use_rep="X_pca", n_pcs=n_pcs, n_neighbors=n_neighbors)
    sc.tl.umap(adata)
    classify_cells(
        adata=adata,
        obs_key=obs_key,
        target_group=target_group,
        n_neighbors=n_neighbors,
    )
    return adata


def vanilla_annotation(
    adata: sc.AnnData,
    obs_key: str,
    target_group: str,
    batch_key: str | None = None,
    layer: str = "X",
    n_pcs=30,
    n_neighbors=20,
) -> sc.AnnData:
    adata = adata.copy()

    if layer == "raw":
        adata = adata.raw.to_adata()
    elif layer != "X":
        adata.X = adata.layers[layer].copy()
        
    adata = log1p(adata, batch_key=batch_key, n_pcs=n_pcs)
    if batch_key is None:
        use_rep = "X_pca"
    else:
        use_rep = "X_pca_harmony"
    classify_cells(
        adata=adata,
        obs_key=obs_key,
        target_group=target_group,
        n_neighbors=n_neighbors,
    )
    return adata


