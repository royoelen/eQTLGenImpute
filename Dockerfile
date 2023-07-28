FROM nfcore/base:2.1
LABEL authors="urmo.vosa@ut.ee" \
      description="Docker image containing all the software requirements for the eQTLGen imputation pipeline. Based on eQTL Catalogue imputation pipeline."

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/eqtlgen-imputation/bin:$PATH

# Dump the details of the installed packages to a file for posterity
# RUN conda env export --name eqtlgen-imputation > eqtlgen-imputation.yml

#Copy other binary dependencies
RUN apt-get clean && apt-get update && apt-get install -y libgomp1 gawk
COPY bin/eagle /usr/bin/
COPY bin/minimac4 /usr/bin/
COPY bin/GenotypeHarmonizer-1.4.23/ /usr/bin/
