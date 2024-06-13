# Epidemic dynamics modeling

### **Input**
- 230822_lineage_notes.txt - A lineage description list retrieved from https://github.com/cov-lineages/pango-designation on August 22, 2023
- hap.mut.info.txt - A matrix showing the presence/absence of mutation clusters in each SARS-CoV-2 haplotype; row: haplotype, column: mutation cluster
- lineage_count_matrix.txt - A count matrix for modeling the relationship between amino acid substitutions in SARS-CoV-2 proteins and viral fitness of SARS-CoV-2; row: day, column: haplotype

### **Script**
- multinomial_mut_regression.stan - A stan code for performing epidemic dynamics modeling using a Bayesian hierarchical multinomial logistic model
- summarize_mut_info.ver2_nextclade.py - A Python script to retrieve amino acid mutations in each SARS-CoV-2 genomic sequence compared to the original strain
