# Cell types annotation of cycling T lymphocytes in the tumor microenvironment
This is a repository for scientific project in the Institute of Bioinformatics, 2023

**Students:**

Anastasiia Mikhailichenko ([github](https://github.com/Mikhailichenko), [telegram](https://t.me/A12nastMi))

Kristina Zheltova ([github](https://github.com/pacifikus), [telegram](https://t.me/masterkristall))

**Superviser:**

Sergey Isaev, Medical University of Vienna

## Introduction

Single-cell transcriptomics is a powerful tool for analyzing the tissues of living organisms. One of its applications is the analysis of the population of T-lymphocytes of people suffering from various diseases, in particular cancer.
Based on the analysis of single-cell transcriptomes of T-cells infiltrating tumors, it is possible to make assumptions about progress of the disease and the appropriate therapy for the patient, as well as to study the dynamics of the development of the disease in the course of fundamental research. From this point of view, it is especially important to know which groups of lymphocytes are cycling. However, the signal of cycling stage markers in such cells is stronger than the signal of cell type markers. When defining clusters, сycling cells usually form a separate cluster, which makes it impossible to tell which cell type the cycling cell belongs to.

The main goal of this work is to develop a method for determining the types of cycling cells and its application for the analysis of T cells from the tumor microenvironment.

We have prepared 3 classification pipelines and to validate their performance we used data from single-cell sequencing of T cell receptors (TCR). 

At a certain stage of maturation of a T-lymphocyte, its T cell receptors mature; during this maturation, the receptor gene sequence undergoes recombination and mutational changes, due to which it becomes unique in the hypervariable CDR3 locus [1]. Subsequently, all descendants of a particular lymphocyte can be determined from this unique sequence. The descendants of a single lymphocyte are called a clonotype. It was noted that most of the cells of the same clonotype belong to the same cell type, but some of them may belong to the group of dividing cells. Determining clonotypes in a group of dividing cells and finding representatives of the same clonotypes in clusters of certain cell types, we can reliably say which cell type the dividing cells belong to. We used these properties to validate the piplines we developed.

## Methods

One of the objectives of the project was to develop and test various methods for classifying cycling T cells and identifying the most suitable cell type. To do this, we have prepared 3 classification pipelines:
- The first approach was based on the hypothesis that if we remove the genes associated with the cell cycle, this will help to determine cell type. Then  a classification was made according to the k nearest neighbors method.
- The second pipeline was to regress out cell cycle effect, so we got the genes whose expression is most associated with the cell cycle and then we again classified and annotated the cells using the knn method.
- The third pipeline is label transfer. This method assumes that if we already have an annotated dataset, for which we know the cell types, we can integrate the data of unknown cells with this existing reference dataset and, based on the environment of unknown cells, understand what type they belong to.

We also compared results with the baseline method - knn classification on the dataset without any preprocessing of the cell cycle signal.

### Preparation of datasets for training and validation of pipelines.

The scanpy package was used as the main tool for working with single-sell transcriptomes. Scirpy was used for TCR analysis.

| ![image](https://github.com/serjisa/cycling_T/blob/main/images/scheme.png) | 
|:--:| 
| *Process of preparing datasets for training and validation* |

The preparation of datasets took place in several stages: the first one was downloading data from the GEO NCBI database and uploading it to Jupiter notebook for further analysis. We analyzed two datasets **GSE162500** (11 patients, 17 sample) and **GSE154826** (35 patients, 77 sample). Next, data quality control was performed. Cells with too low and high counts, and with too high percentage of mitochondrial genes were filtered out. For each batch, doublets were filtered out using the scrublet tool. At the end, counts were normalized and log-transformed. 

After that, datasets filtered from unnecessary data went through a multi-stage annotation of cell types. At each stage, 3000 highly variable genes were selected, batch correction using the harmony algorithm, and clustering using the leiden algorithm. First of all, we removed cells that did not express the CD3E T lymphocyte marker. We found **54971** CD3+ in the GSE162500 dataset and **126910** in the GSE154826 dataset After filtering out CD3 cells from the original datasets, we again selected highly variable genes, performed batch correction and clustering. After that, we obtained clusters of three types CD8+, CD4+, CB 4/8 cycling (expressing the division marker MKI67) and CD4/8 IFN (expressing the interferon signal marker gene ISG15). CD8+ and CD4+ cells were divided into 2 datasets, and they were separately selected for highly variable genes, batch correction, and clustering. After that, according to the expression of a number of marker genes taken from the article Gueguen *et al*., 2021 [2]. The cell type was determined for each of the clusters.

| ![image](https://github.com/serjisa/cycling_T/blob/main/images/UMAPS1.png) | 
|:--:| 
| *UMAPs obtained during the annotation of dataset **GSE162500*** |

| ![image](https://github.com/serjisa/cycling_T/blob/main/images/hitmaps1.png) | 
|:--:| 
| *Hitmaps used for annotation of cell types of **GSE162500** dataset* |


| ![image](https://github.com/serjisa/cycling_T/blob/main/images/UMAPS2.png) | 
|:--:| 
| *UMAPs obtained during the annotation of dataset **GSE154826*** |

| ![image](https://github.com/serjisa/cycling_T/blob/main/images/hitmaps2.png) | 
|:--:| 
| *Hitmaps used for annotation of cell types of **GSE154826** dataset* |

The last stage of dataset processing was working with T cell receptors. The scirpy package was used for this. First, the types of T cell receptors for each cell were determined and only cells carrying at least one alpha and beta chain were selected (**27493** in GSE162500 and **11838** in GSE154826 dataset). Clonotypes have also been identified. 

| ![image](https://github.com/serjisa/cycling_T/blob/main/images/clonotypes.png) | 
|:--:| 
| *UMAPs showing cells in largest clonotypes and cell types of all cells for GSE162500 and GSE154826 datasets* |

Datasets for validation contained information on counts for all CD3+ lymphocytes, as well as information on cell types, TCR, and clonotypes.
Moreover, we selected clonotypes that satisfy two conditions: the number of cells in them is more than 10 and at least one of the cells of the clonotype belongs to the cluster of cycling cells. A table was also compiled for validation with the frequencies of occurrence of representatives of the selected clonotypes in all cell types.

## Results


| ![image](https://github.com/serjisa/cycling_T/assets/22592039/1f00f1ba-1b4a-4673-964a-72660e8161ef) | 
|:--:| 
| *Examples of an original UMAP and new labels for **removed genes** pipeline on **GSE162500** dataset* |

| ![image](https://github.com/serjisa/cycling_T/assets/22592039/4e2356ce-cf65-4dba-ac77-3969dac6c01f) | 
|:--:| 
| *Examples of an original UMAP and new labels for **cell cycle scoring** pipeline on **GSE162500** dataset* |


| ![image](https://github.com/serjisa/cycling_T/assets/22592039/6d3b7211-dca3-4109-8551-16008ddb9087) | 
|:--:| 
| *Examples of an original UMAP and new labels for **label transfer** pipeline on **GSE162500** dataset* |

Overall the label transfer method shows a better mean quality and a smaller variance of metrics, however, we would like to test our algorithms on more datasets to obtain more stable results and generalize them.

| ![image](https://github.com/serjisa/cycling_T/assets/22592039/c593ce43-019b-4a41-8c47-1f19b0dc8429) | 
|:--:| 
| *F1-score boxplots for all pipelines* |

| ![image](https://github.com/serjisa/cycling_T/assets/22592039/9dd9ed75-ae7e-40f4-a9a9-24c3446c9a43) | 
|:--:| 
| *F1-score for all pipelines* |



## Repository structure

    ├── datasets_preparation 
    │   ├── data                                 <- Instructions for downloading the source data
    │   ├── GSE154826_analysis_wo_outputs.ipynb  <- Preparation for the GSE154826 dataset
    │   └── GSE162500_analysis_wo_outputs.ipynb  <- Preparation for the GSE162500 dataset
    ├── images                                   <- Images for main README
    ├── src                                      <- Source code of pipelines, validation and visuzation
    └── annotation_pipelines.ipynb               <- Run annotation pipelines on prepared data

## Minimum requirements

- Hardware
    - CPU: 4 CPU Cores
    - RAM: 50 GB
    - System disk space: 50 GB
 - Software
    - Python version: 3.10.5
    - OS: Ubuntu 20.04.4 LTS
    - All used packages are listed in `requirements.txt`

## References
1) Attaf, M., Huseby, E., Sewell, A. αβ T cell receptors as predictors of health and disease. Cell Mol Immunol 12, 391–399 (2015). https://doi.org/10.1038/cmi.2014.134
2) Gueguen P, Metoikidou C, Dupic T, Lawand M, Goudot C, Baulande S, Lameiras S, Lantz O, Girard N, Seguin-Givelet A, Lefevre M, Mora T, Walczak AM, Waterfall JJ, Amigorena S. Contribution of resident and circulating precursors to tumor-infiltrating CD8+ T cell populations in lung cancer. Sci Immunol 6(55), (2021). https://doi.org/10.1126/sciimmunol.abd5778
