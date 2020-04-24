# VAT (Variant Annotation Tool)

This is a variant annotation tool developed using R studio Version 1.2.1335. Each variant is annotated with the following pieces of information:
1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple
possibilities, it was annotated with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation
3. Number of reads supporting the variant
4. Percentage of reads supporting the variant versus those supporting reference reads
5. Allele frequency of variant is collected from Broad Institute ExAC Project API
