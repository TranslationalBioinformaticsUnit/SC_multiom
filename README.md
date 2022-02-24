# SC_multiom
Description of the work done and to be done

# DESCRIPTION

SAMPLES: one normal SCRNAseq, and 5 multiom; 1 fresh and 4 frozen. One of the frozens is working good(MM_HR_297353). MM_HR_907052 used the RE sample

ANALYSIS DONE: 
  01. preprocessing of each individual and filter the data. The data use filtered first based on the fresh sample and then with the standar filters (check sc_multiom.R script)
  2. check cluster markers and most variable 2000 genes of each of the samples (fresh and frozen)
  3. Integrate the samples using Seurat
  4. Integrate the samples using harmony
  5. Check ambient RNA with dropletQC, but not very good results
  6. Start with soupX, also a tool to check ambient RNA

NEXT STEPS: 
  1. Check soupX
  2. Run findoublets pipeline to check doublets
  3. Check data of bruno in the folder /mnt/resultados/oncohematologia/resultados/mlagasag/SC_MIREN.tar.gz
  
  # SCRIPTS
  
  01_demultiplexing.sh: script to demultiplex and create fastq files
  
  02_counts.sh: script to create counts matrix
  
  dropletQC: dropletQC  library to eliminate ambient RNA
  
  integration.R : integration pipeline to integrate all the frozen and the fresh sample with seurat pipeline
  
  integration.sh: to run the integration.R in the server
  
  QC_merge: different plots and approches done in the quality control steps
  
  RefGenome.sh: script to create the reference genome for the multiom case
  
  SC_multiom.R: the preprocessing step of each of the sample individually. Here the filters are defined
  
  SC_multiom.sh: script to run the sc_multiom in the server
  
  soupX: soupX pipeline to correct ambient RNA

  Scripts are located in: /mnt/resultados/oncohematologia/resultados/mlagasag/Riney_SC_analysis/SCRIPTS
 
 # DATA
 fastq_files: /mnt/resultados/oncohematologia/resultados/mlagasag/Riney_SC_fastq/
 
 count_matrix: /mnt/resultados/oncohematologia/resultados/mlagasag/Riney_SC_counts/
 
 analysis carried out: /mnt/resultados/oncohematologia/resultados/mlagasag/Riney_SC_analysis/
