# Human variation in population-wide gene expression data predicts gene function and phenotype
</br>
 
### *Lorenzo Bonaguro, Jonas Schulte-Schrepping, Caterina Carraro, Laura Sun, Benedikt Reiz, Ioanna Gemünd, Adem Saglam, Souad Rahmouni, Michel Georges, Peer Arts, Alexander Hoischen, Leo A. B. Joosten, Frank L. van de Veerdonk, Mihai G. Netea, Kristian Händler, Sach Mukherjee, Thomas Ulas, Joachim L. Schultze, Anna C. Aschenbrenner*

</br></br>


![image](./images/abstract.png)

## Abstract
Population-scale multi-layered datasets assemble extensive experimental data of different types on single healthy individuals in large cohorts, capturing genetic variation and environmental factors influencing gene expression with no clinical evidence of pathology. Variance of gene expression can be exploited to set up a conditional quasi loss- and gain-of-function “in population” experiment if expression values for the gene of interest (GOI) are available. We describe here a novel approach, called huva (human variation), that takes advantage of population-scale multi-layered data to infer gene function and relationships between phenotypes and gene expression. Within a reference dataset, huva derives two experimental groups, i.e. individuals with LOW or HIGH expression of the GOI, enabling the subsequent comparison of their transcriptional profile and functional parameters. We demonstrate that this approach robustly and efficiently identifies the phenotypic relevance of a GOI, allows the stratification of genes according to shared biological functions, and we further generalized this concept to almost 16,000 genes in the human blood transcriptome. Additionally, we describe how huva predicts the phenotype of naturally occurring activating mutations in humans. Here, huva predicts monocytes rather than lymphocytes to be the major cell type in the pathophysiology of STAT1 activating mutations, evidence which was validated in a cohort of clinically characterized patients.

## Why this repository?
This repository provides the original code used to prepared each panel of the manuscript, to ensure the reproducibility of the analysis we provide a docker image for and Rstudio session including all required dependencies and the *huva* package (v. 0.1.4) and the *huva.db* (v. 0.1.4-2).

## How to use this repository
We try to provide access to the analysis in the easiest possible way, the user can follow this few instructions and should be able to be up and running quickly. Withing each folder you can find all the scripts and data required the reproduce the analysis. Just open the 'Figure_X.Rmd' and run it.  

### Requirements
Some of the calculations, especially to reproduce Figure 3 can be quite memory demanding, we suggest a minimum of 32 Gb of system memory.
- [Git](https://git-scm.com/)
- [Docker](https://www.docker.com/)

### Cloning the repository
```sh
# Make a local copy of the repository
git clone https://github.com/lorenzobonaguro/huva_reproducibility

# Navigate to the folder
cd huva_reproducibility
```
### Download the Docker image and start the container
```
docker pull lorenzobonaguro/huva_rep:v4
```

alternatively you can build the docker image, before starting make sure the source code of the *huva* and *huva.db* packages is in the same folder as the dockerfile
```sh
# Build the Docker image
docker build -t huva_rep:v4 . # Skip this step if you want to use the Docker image from DockerHub
```

```sh
# Run a container
docker run -dp 8787:8787 -e USER=mariorossi -e PASSWORD=mariorossi --name rep_huva -v 'your_directory':/home/mariorossi/data/ rep:v01
```

### Open the RStudio session
Enter in your browser `localhost:8787`, this should start a Rstudio session you can use to explore the code and reproduce the analysis.

## How to cite *huva*
If you use *huva* in your research project consider citing us [linktojournal](https://www.cell.com/iscience/fulltext/S2589-0042(22)01600-5#%20).

## Contact or follow us
For any problem of question regrding the *huva* package or this repositoy or you just want to be up to date on what is coming next, send us an [email](mailto:lorenzobonaguro@uni-bonn.de) or follow us:  

<img src="./images/twitter.png" width="12%" style="float: left;">  

[@LorenzoBonaguro](https://twitter.com/LorenzoBonaguro)  
[@AschenbrennerAC](https://twitter.com/AschenbrennerAC)  
[@LabSchultze](https://twitter.com/LabSchultze)