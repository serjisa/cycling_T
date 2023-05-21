import pickle
import scanpy as sc
import symphonypy as sp
import matplotlib.pyplot as plt

from annotation import annotation_with_removed_genes, annotate_with_cc_scoring, vanilla_annotation
from label_transfer import build_reference, get_query, label_transfer


def get_cell_cycle_genes(
    adata, 
    g2m_phase_file, 
    s_phase_file, 
    need_preprocessing=True,
    return_both=False,
):
    with open(g2m_phase_file) as f:
        g2m_genes = f.read().splitlines()

    with open(s_phase_file) as f:
        s_genes = f.read().splitlines()
    
    if need_preprocessing:
        all_genes = list(map(str.title, (g2m_genes + s_genes)))
    else:
        all_genes = s_genes + g2m_genes
    cell_cycle_genes = [x for x in all_genes if x in adata.var_names]
    if return_both:
        return s_genes, g2m_genes
    return cell_cycle_genes


def base_preprocess(adata, batch_key):
    adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=10e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key=batch_key)
    adata.raw = adata

    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=30)
    sp.pp.harmony_integrate(adata, key=batch_key, max_iter_harmony=50)
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)
    
    sc.pl.umap(adata, color='cell_type', title='Initial UMAP', size=9)
    return adata

    
def run_annotation(
    adata,
    batch_name,
    obs_key, 
    target_group,
    layer,
    cell_cycle_genes=None,
    annotation_type='removed_genes',
    n_pcs=30,
    n_neighbors=20,
):
    adatas_corrected = {}

    for batch_key in adata.obs[batch_name].unique():
        if len(adata[adata.obs[batch_name] == batch_key]) != 1:
            batch = adata[adata.obs[batch_name] == batch_key]
            if len(batch[batch.obs[obs_key] == target_group]) != 0:
                if annotation_type == 'removed_genes':
                    adatas_corrected[batch_key] = annotation_with_removed_genes(
                        adata[adata.obs[batch_name] == batch_key],
                        genes=cell_cycle_genes,
                        obs_key=obs_key,
                        target_group=target_group,
                        layer="raw",
                        n_pcs=n_pcs,
                        n_neighbors=n_neighbors,
                    )
                elif annotation_type == 'vanilla':
                    adatas_corrected[batch_key] = vanilla_annotation(
                        adata[adata.obs[batch_name] == batch_key],
                        obs_key=obs_key,
                        target_group=target_group,
                        layer="raw",
                        n_pcs=n_pcs,
                        n_neighbors=n_neighbors,
                    )
                elif annotation_type == 'cc_scoring':
                    adatas_corrected[batch_key] = annotate_with_cc_scoring(
                        adata[adata.obs[batch_name] == batch_key],
                        obs_key=obs_key,
                        target_group=target_group,
                        layer="raw",
                        n_pcs=n_pcs,
                        n_neighbors=n_neighbors,
                    )
                else:
                    try:
                        adata_ref = build_reference(
                            adata[adata.obs[batch_name] == batch_key],
                            obs_key=obs_key,
                            target_group=target_group,
                            layer='raw',
                        )

                        sc.pp.neighbors(adata_ref)
                        sc.tl.umap(adata_ref)
                        sc.tl.leiden(adata_ref)

                        sc.pl.umap(
                            adata_ref,
                            color=[obs_key],
                            frameon=False,
                            title=["Provided cell type", "Clusters"],
                        )

                        adata_query = get_query(
                            adata[adata.obs[batch_name] == batch_key],
                            obs_key=obs_key,
                            target_group=target_group,
                            layer='raw',
                        )
                        adata_query = label_transfer(
                            adata_query,
                            adata_ref,
                            obs_key=obs_key,
                            target_group=target_group,
                        )
                        adatas_corrected[batch_key] = adata_query
                    except:
                        print(f'error for batch {batch_key}')
    return adatas_corrected

                
def plot_umaps(adatas_corrected, obs_key, target):
    for batch in adatas_corrected:
        original_cmap = dict(zip(
            adatas_corrected[batch].obs.cell_type.cat.categories,
            adatas_corrected[batch].uns["cell_type_colors"],
        ))
        adatas_corrected[batch].obs["labels_pred"] = adatas_corrected[
            batch
        ].obs["labels_pred"].astype("category")
        adatas_corrected[batch].uns["labels_pred_colors"] = [
            original_cmap[ct] for ct in adatas_corrected[batch].obs.labels_pred.cat.categories
        ]

        fig, axes = plt.subplots(ncols=2, figsize=(8, 4))
        sc.pl.umap(
            adatas_corrected[batch],
            frameon=False,
            ax=axes[0],
            color=obs_key,
            title=f"{batch}\n(original labels)",
            legend_loc=None,
            show=False,
        )
        sc.pl.umap(
            adatas_corrected[batch],
            frameon=False,
            ax=axes[1],
            show=False,
        )
        sc.pl.umap(
            adatas_corrected[batch][
                adatas_corrected[batch].obs[obs_key] == target
            ],
            frameon=False,
            ax=axes[1],
            color="labels_pred",
            legend_loc=None,
            show=False,
            title=f"{batch}\n(new labels)",
            size=120000 / adatas_corrected[batch].shape[0],
        )
        
def save_pkl(obj, filepath):
    with open(filepath, "wb") as f:
        pickle.dump(obj, f)
        

def load_pkl(filepath):
    with open(filepath, "rb") as f:
        return pickle.load(f)
