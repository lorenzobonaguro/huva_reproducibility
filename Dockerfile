# Docker container for an Rstudios session for huva analysis
FROM bioconductor/bioconductor_docker:RELEASE_3_12

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    libv8-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages (requitement for the huva package)
RUN Rscript -e 'BiocManager::install(update = T, ask = F)' && \
	Rscript -e 'BiocManager::install(c("ggplot2", "Rmisc", "ggpubr", "reshape2", "ggsci", "plotly", "knitr", "pheatmap", "useful", "rmarkdown"), version = "3.12")' && \
	Rscript -e 'BiocManager::install(c("limma", "GSVA"), version = "3.12")'

COPY fgsea_1.12.0.tar.gz /tmp/fgsea_1.12.0.tar.gz
RUN Rscript -e 'install.packages("/tmp/fgsea_1.12.0.tar.gz", repos = NULL, type = "source")'

# Install R packages (required for other analyses)  
RUN Rscript -e 'BiocManager::install(c("tictoc", "tidyverse", "viridis", "combinat", "ComplexHeatmap", "pcaGoPromoter.Hs.hg19", "ggnetwork", "intergraph", "MCDA", "reactome.db", "ReactomePA", "Hmisc", "gtools", "biomaRt", "foreach", "gplots", "ggbeeswarm", "factoextra", "VennDiagram", "rhdf5", "tximport", "DESeq2", "vsn", "genefilter", "IHW"), version = "3.12")'

COPY org.Hs.eg.db_3.8.2.tar.gz /tmp/org.Hs.eg.db_3.8.2.tar.gz
COPY clusterProfiler_3.12.0.tar.gz /tmp/clusterProfiler_3.12.0.tar.gz
RUN Rscript -e 'install.packages("/tmp/org.Hs.eg.db_3.8.2.tar.gz", repos = NULL, type = "source")' && \
	Rscript -e 'install.packages("/tmp/clusterProfiler_3.12.0.tar.gz", repos = NULL, type = "source")'

# Install huva
COPY huva_0.1.4.tar.gz /tmp/huva_0.1.4.tar.gz
COPY huva.db_0.1.4.tar.gz /tmp/huva.db_0.1.4.tar.gz
RUN Rscript -e 'install.packages("/tmp/huva.db_0.1.4.tar.gz", repos = NULL, type = "source")' && \
	Rscript -e 'install.packages("/tmp/huva_0.1.4.tar.gz", repos = NULL, type = "source")'