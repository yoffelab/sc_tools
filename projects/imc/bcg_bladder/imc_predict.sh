#!/bin/bash
SAMPLES=data/*.mcd
TIFFS=processed/*/tiffs/*_full.tiff
for TIFF in ${TIFFS}
do
	BASE_TIFF="$(basename "$TIFF" .tiff)"
	echo PROCESSING: ${BASE_TIFF}
	
	## Output pixel probabilities of nucleus, membrane and background using ilastik
	imc predict $TIFF
	## Segment cell instances with DeepCell
	imc segment \
	  --from-probabilities \
	  --model deepcell \
	  --compartment both $TIFF
	# Quantify channel intensity and morphology for each single cell in every image
	imc quantify $TIFF
	mv processed/quantification.h5ad processed/quantification/quantification_${BASE_TIFF}.h5ad
done
