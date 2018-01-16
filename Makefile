project_root := /projects/ezhao_prj/analyses/mutation-signature-pipeline
analysis_root := $(project_root)/analysis

###########################
### Directory Structure ###
###########################

all: \
	scripts/SignIT \
	scripts/pipelines \
	scripts/SignatureEstimation

meta:
	mkdir -p meta

paths:
	mkdir -p paths

scripts:
	mkdir -p scripts

#################
### Load Code ###
#################

scripts/SignIT: scripts
	if [ -d $@ ]; \
	then(cd $@ && git pull); \
	else git clone git@github.com:eyzhao/SignIT.git $@; \
	fi

scripts/pipelines: scripts
	if [ -d $@ ]; \
	then (cd $@ && git pull); \
	else git clone git@github.com:eyzhao/bio-pipelines.git $@; \
	fi

scripts/SignatureEstimation: scripts
	rm -rf $@ \
		&& wget https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.tar.gz \
		&& mv SignatureEstimation.tar.gz scripts \
		&& cd scripts \
		&& tar -xvf SignatureEstimation.tar.gz \
		&& rm -rf SignatureEstimation.tar.gz

