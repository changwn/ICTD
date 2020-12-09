# Description

This method is described in the publication from Biorxiv, 2018 available at [[https://www.biorxiv.org/content/10.1101/426593v2]](https://www.biorxiv.org/content/10.1101/426593v3)

To cite our software, please use the independent DOI from Zenodo.
[![DOI](https://zenodo.org/badge/224000547.svg)](https://zenodo.org/badge/latestdoi/224000547)


ICTD web application demo is available at : [[https://shiny.ph.iu.edu/ICTD/]](https://shiny.ph.iu.edu/ICTD/)

ICTD web application tutorial is available at : [[https://github.com/changwn/ICTD/blob/master/vignettes/ICTD_server_tutorial.md]](https://github.com/changwn/ICTD/blob/master/vignettes/ICTD_server_tutorial.md)

![[image]](img/web_app.png)

Backup link: [[https://ictd.ccbb.iupui.edu]](https://ictd.ccbb.iupui.edu)

# ICTD Framework
![[fig1]](img/fig1.png)

# Installation

```
#install ICTD
install.packages("devtools")
devtools::install_github("changwn/ICTD")
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

- [Wennan Chang](https://changwn.github.io/)
(wnchang@iu.edu)

PhD candidate at BDR group, Indiana University School of Medicine

- [Chi Zhang](https://medicine.iu.edu/departments/genetics/faculty/27057/zhang-chi/)
(czhang87@iu.edu)

Assistant Professor

Department of Medical & Molecular Genetics, Indiana University School of Medicine



# Dependencies

We also provide a Docker image to recreate the compute environment. See the Dockerfile for more details.

[[https://hub.docker.com/r/wnchang/ictd]](https://hub.docker.com/r/wnchang/ictd)

Using the Docker image could void the conflict issue that R version and several R packages version confict. 

For more details about the Docker, please see Docker documentation page [[https://docs.docker.com/]](https://docs.docker.com/).
