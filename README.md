# Description

This method is described in the publication from Biorxiv, 2018 available at https://www.biorxiv.org/content/10.1101/426593v2

ICTD web application demo is available at : https://shiny.ph.iu.edu/ICTD/

ICTD web application tutorial is available at : https://github.com/changwn/ICTD/blob/master/vignettes/ICTD_server_tutorial.md

![image](https://github.com/zy26/ICTD/blob/master/img/web_app.png)

Backup link: https://ictd.ccbb.iupui.edu


# Installation

```
#install dependent pkg
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("impute", version = "3.8")
BiocManager::install("GO.db", version = "3.8")
BiocManager::install("sva", version = "3.8")
BiocManager::install("preprocessCore", version = "3.8")

rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)

#install ICTD
install.packages("devtools")
devtools::install_github("changwn/ICTD")
```



Note : For old R version which cannot install 'BiocManager', please use below command to install the dependency.
```
source("https://bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("GO.db")
biocLite("sva")
biocLite("preprocessCore")

rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
```

# Example

```
library(ICTD)

data_bulk = GSE72056_diri_example[[1]]
ictd_result <- ICTD(data_bulk)

#Return value is a list, which the first element is the predicted proportion and 
#the second element is the predicted markers of ICTD
```

# Questions & Problems

If you have any questions or problems when using ICTD, please feel free to open a new issue [here](https://github.com/zy26/ICTD/issues). We will fix the new issue ASAP.  You can also email the maintainers and authors below.

- [Wennan Chang](https://zcslab.github.io/people/wennan/)
(wnchang@iu.edu)

PhD candidate at BDR group, Indiana University School of Medicine

- [Chi Zhang](https://medicine.iu.edu/departments/genetics/faculty/27057/zhang-chi/)
(czhang87@iu.edu)

Assistant Professor

Department of Medical & Molecular Genetics, Indiana University School of Medicine



# Dependencies

We also provide a Docker image to recreate the compute environment. See the Dockerfile for more details.

https://hub.docker.com/r/wnchang/ictd

Using the Docker image could void the conflict issue that R version and several R packages version confict. 

For more details about the Docker, please see Docker documentation page https://docs.docker.com/.
