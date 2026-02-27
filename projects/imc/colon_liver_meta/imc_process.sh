#!/bin/bash
SAMPLES=data/*.mcd
for SAMPLE_MCD in ${SAMPLES}
do
	SAMPLE="$(basename "$SAMPLE_MCD" .mcd)"
	echo PROCESSING: ${SAMPLE}
	## output description of acquired data
	imc inspect data/${SAMPLE}.mcd
	## convert MCD to TIFFs and auxiliary files
	imc prepare \
	  --ilastik \
	  --n-crops 0 \
	  --ilastik-compartment nuclear \
	  data/${SAMPLE}.mcd
done

