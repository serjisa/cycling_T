Load index.html file

```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/

```

Create files_names.txt with names of files which are going to be loaded

```
cat index.html | cut -f 2 -d ">"| cut -f 1 -d "<" | grep GSE154826_ > files_names.txt
```

Create commands_download.txt with all commands needed to load datasets

```
awk ' { print "wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/" $0 } ' files_names.txt > commands_download.txt
```

unzip files from .gz

```
gunzip *gz
```


Create commands_tar.txt with commands to create directories 
```
cat files_names.txt | grep .tar | cut -f 1 -d "." | awk ' { print "mkdir " $0 } ' > commands_tar.txt
```
Create commands_tar.txt with commands to unzip from .tar
```
cat files_names.txt | grep .tar | cut -f 1 -d "." | awk ' { print "tar -xf " $0 ".tar --directory " $0 } '  > commands_tar.txt
```
Create commands_gz.txt with commands to unzip files from .gz in all directories
```
cat files_names.txt | grep .tar | cut -f 1 -d "." | awk ' { print "gzip " $0 "/*" } '  > commands_gz.txt
```
