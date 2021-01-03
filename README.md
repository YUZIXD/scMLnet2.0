# scMLnet2.0

## Introduction

scMLnet2.0 is a development release of scMLnet R package developed to construct inter-/intracellular multilayer singaling network based on single-cell RNA-seq expression data. scMLnet constructs the multilayer network by integrating intercellular pathways (ligand-receptor interactions) and intracellular subnetworks (receptor-TF pathways and TF-target gene interactions) based on cell-type specific gene expression, prior network information and statistical inference.

## Update

* Expand prior knowledge databases, collecting information from molecular interaction databases, pathway databases, regulatory databases, and tools aimed for cell-cell communication.
* Construct a directed weighted network based on prior information for visualization of multi-layer networks in cell-cell communication.
* Optimized and adjusted the code for screening single-layer network (LigRecTab, RecTFTab, TFTGTab) to adapt to the subsequent analysis steps and improve operating efficiency.
