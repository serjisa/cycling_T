Download datasets and unzip

```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162500/suppl/GSE162500_RAW.tar
tar -xf GSE162500_RAW.tar
```

Unzip .tar.gz
```
tar -xf GSM4952953_P34_Tumor_raw_feature_bc_matrix.tar.gz
tar -xf GSM4952954_P35_Tumor_raw_feature_bc_matrix.tar.gz
tar -xf GSM4952955_P42_Tumor_raw_feature_bc_matrix.tar.gz
tar -xf GSM4952956_P43_Tumor_raw_feature_bc_matrix.tar.gz
tar -xf GSM4952957_P46_Tumor_raw_feature_bc_matrix.tar.gz
tar -xf GSM4952958_P47_Tumor_raw_feature_bc_matrix.tar.gz
tar -xf GSM4952959_P55_Tumor_raw_feature_bc_matrix.tar.gz
tar -xf GSM4952960_P57_Tumor_raw_feature_bc_matrix.tar.gz
tar -xf GSM4952961_P57_Blood_raw_feature_bc_matrix.tar.gz
tar -xf GSM4952962_P58_Tumor_filtered_feature_bc_matrix.tar.gz
tar -xf GSM4952963_P58_Blood_filtered_feature_bc_matrix.tar.gz
tar -xf GSM4952964_P60_Tumor_filtered_feature_bc_matrix.tar.gz
tar -xf GSM4952965_P60_Juxta_filtered_feature_bc_matrix.tar.gz
tar -xf GSM4952966_P60_Blood_filtered_feature_bc_matrix.tar.gz
tar -xf GSM4952968_P61_Juxta_filtered_feature_bc_matrix.tar.gz
tar -xf GSM4952969_P61_Blood_filtered_feature_bc_matrix.tar.gz
tar -xf GSM4952970_P47_Tumor_filtered_contig_annotations.csv.gz
rm *.tar.gz
```
Unzip .gz
```
zcat GSM4952970_P47_Tumor_filtered_contig_annotations.csv.gz > P47_Tumor_filtered_contig_annotations.csv
zcat GSM4952971_P55_Tumor_filtered_contig_annotations.csv.gz > P55_Tumor_filtered_contig_annotations.csv
zcat GSM4952972_P57_Tumor_filtered_contig_annotations.csv.gz > P57_Tumor_filtered_contig_annotations.csv
zcat GSM4952973_P57_Blood_filtered_contig_annotations.csv.gz > P57_Blood_filtered_contig_annotations.csv
zcat GSM4952974_P58_Tumor_filtered_contig_annotations.csv.gz > P58_Tumor_filtered_contig_annotations.csv
zcat GSM4952975_P58_Blood_filtered_contig_annotations.csv.gz > P58_Blood_filtered_contig_annotations.csv
zcat GSM4952976_P60_Tumor_filtered_contig_annotations.csv.gz > P60_Tumor_filtered_contig_annotations.csv
zcat GSM4952977_P60_Juxta_filtered_contig_annotations.csv.gz > P60_Juxta_filtered_contig_annotations.csv
zcat GSM4952978_P60_Blood_filtered_contig_annotations.csv.gz > P60_Blood_filtered_contig_annotations.csv
zcat GSM4952979_P61_Tumor_filtered_contig_annotations.csv.gz > P61_Tumor_filtered_contig_annotations.csv
zcat GSM4952980_P61_Juxta_filtered_contig_annotations.csv.gz > P61_Juxta_filtered_contig_annotations.csv
zcat GSM4952981_P61_Blood_filtered_contig_annotations.csv.gz > P61_Blood_filtered_contig_annotations.csv
```