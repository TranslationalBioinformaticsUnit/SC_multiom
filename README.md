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
  
