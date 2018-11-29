==========
ENCODE ChIP-seq Pipeline
==========

ENCODE Uniform processing pipeline for ChIP-seq
===============================================

Current implementation is deployed to the DNAnexus platform.

Mapping
-------

1. Map reads with BWA, mark duplicates Picard, and remove duplicates.
2. Estimate library complexity and calculate calculate NRF (non-redundant fraction), PBC1, PBC2 (PCR bottleneck coefficient).
3. Calculate cross-correlation analysis with spp/phantompeakqualtools.
4. Generate p-value and fold-over-control signal tracks for each replicate and replicates pooled with MACS2.

Peak calling (histone marks)
----------------------------

1. Call peaks with MACS2.
2. Calculate and report overlapping peaks from both replicates.


Peak calling (transcription factors)
------------------------------------

1. Call peaks with SPP.
2. Threshold peaks with IDR.
3. Report IDR-thresholded peak sets, self-consistency ratio, rescue ratio, reproducibility test.

NOTE: Pipeline V2 can be found: https://github.com/ENCODE-DCC/chip-seq-pipeline2
