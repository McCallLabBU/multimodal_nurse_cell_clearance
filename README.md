# Multi-modal comparison of molecular programs driving nurse cell death and clearance in Drosophila melanogaster oogenesis 

Scripts for the analysis and visualization of - 

- Translatome (TRAP-seq)
- Secretome (LC-MS/MS)
- transcriptome (integrated analysis of Jevitt 2020 and Rust 2020 single-cell RNAseq)

Paper : 



## Figures

- [x] **Figure 1: Establishing the “translatome”, “secretome”, and “transcriptome” of MBFCs and SFCs**
  - [x] grahical abstract
  - [x] experimental design
  - [x] scRNAseq integration 
- [x] **Figure 2: Disparate genes contribute to overlaps in GO-Terms enriched in SFCs and starved MBFCs**
  - [x] Translatome fold-change quadrant dot plot 
  - [ ] Translatome # of genes congruent/incongruent
  - [x] Translatome GO-term comparison 
  - [x] Secretome fold-change quadrant dot plot 
  - [ ] Secretome # of genes congruent/incongruent 
  - [x] Secretome GO-term comparison
  - [ ] Translatome vs. ScRNAseq fold-change quadrant dot plot 
  - [ ] Translatome vs. ScRNAseq  # of genes congruent/incongruent
  - [ ] Translatome vs. ScRNAseq GO-term comparison
- [ ] **Figure 3: Cytoskeletal dynamics in SFCs and starved MBFCs**
  - [ ] Heatmap of counts of cytoskeletal genes from all 3 datasets
  - [ ] RNAi quantification and images - 
    - [ ] Tm1, Up, Act79B (Shruthi)
    - [ ] TpnC41C, Scp1, cue (Diane Thesis figure 5.13 page 246) 
- [ ] **Figure 4: Hsp23 translation is upregulated in both SFCs and starved MBFCs**
  - [ ] Hsp23 RNAi and OE in conditioned + discs-large staining
  - [ ] Hsp23 OE and RNAi in starved + discs-large staining
- [ ] **Figure 5:** **Toll pathway is activated in SFCs by** **translation** **upregulation of** **Spätzle-Processing Enzyme (SPE)** 
  - [ ] Heatmap
  - [ ] RNAi - SPE, PGRP-LC, Drs, Drsl4  (+discs-large staining)
- [ ] **Figure 6: Starvation elicits translation upregulation of antimicrobial peptide *AttD* in MBFCs**
  - [ ] RNAi - AttD  (+discs-large staining)
  - [ ] Quantification of mid-stage death in AttD + other starvation associated genes
- [ ] **Figure 7: RNAi uncovers both clearance and morphological defects** 
  - [ ] shape analysis 
- [ ] **Figure 8: CG gene validation** 
  - [ ] CG17192, CG6277, CG17633, CG10163, CG31926, CG6508, CG12398, CG31997, CG14259, CG31999, CG14645, CG13618, CG31659, CG13998, CG6403, CG3397, CG11600

## Supplemental figures

- Figure 1
  - Rpl10a-GFP in GR1 and PG150 (Diane thesis Figure 5.2, page 186)
  
  - HRP-KDEL in GR1 and PG150 (Albert thesis Figure 6.3, page 150)
  - HRP-KDEL - streptavidin, a-V5 (Albert thesis Figure 6.4, page 151)
  - Biotinylation western and ponceau (Albert thesis Figure 6.5, page 152)
  - bead enrichment western (Albert thesis Figure 6.6, page 154)



## Supplemental tables

- Translatome Starvation vs phagoptosis - each gene categorized congruency 
- Translatome Starvation vs phagoptosis - go terms in each category
- Compiled list of DE genes/proteins in all 3 datasets in phagoptosis comparison
- Compiled list of DE genes/proteins in all 3 datasets in starvation comparison
- Validation summary
  - clearance defects quantified
  - morphological defects summarized

## Scripts

### Translatome

- **Preprocessing**

  - 0_concatenate_samples.sh
  - 1_qc_raw_fastqs.sh
  - 2_trim_adapters.sh
  - 3_salmon_quantify.sh

- **Differential expression**

  - utils.R 

  - deseq2.R 


### Secretome

- qc.R 
- secretome_differentialabundance.R

### Single-cell RNA-seq

- Integration_utils.R 
- Integrate_rust2020_jevitt2020.R
- mbfc_vs_sfc.R

### Compilation of results

- multimodal_comparison.R 



## To Do List

- [ ] Write GO-terms results section
- [ ] Write Toll pathway activation result
- [ ] Write AttD/starvation genes result 
- [ ] Clean up code base and push to GitHub
  - [x] Clean up translatome
  - [x] secretome
  - [x] ScRNAseq
  - [ ] compilation
- [ ] Shape analysis of morphological defects
- [ ] Write significance statement
- [ ] Write discussion 
  - [ ] summary of findings
  - [ ] future directions
- [ ] Prepare figures 
- [ ] Summarize fly stocks and reagents 
- [ ] Find confocal images of previously validated genes in Diane and Albert's thesis 
- [ ] Find accessory experiment data and images in Diane and Albert's thesis
- [ ] Check if all citations are correct and more are required 
