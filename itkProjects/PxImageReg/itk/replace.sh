#!/bin/bash

for x in itkCorr*
do
	export f=itkEntropy${x#itkCorrelation}
	sed 's/Correlation/Entropy/g' ${x} > ${f}
done;
