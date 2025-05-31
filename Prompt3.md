## Class Test 1 - 2023

### Question 1 [2 marks]
**During Illumina sequencing library preparation, adapters serve two key roles in cluster amplification:**

1. **Forward and reverse primers for bridge amplification**: The adapters contain complementary sequences to oligonucleotides immobilized on the flow cell surface, enabling the DNA fragments to bind and initiate bridge amplification.

2. **Template for cluster generation**: The adapters provide the necessary sequences for repeated cycles of amplification, creating clonal clusters of identical DNA fragments that can be sequenced simultaneously.

### Question 2 [5 marks]
**i) Difference between global and local pairwise alignment (3 marks):**
- **Global alignment**: Attempts to align sequences from end-to-end, considering the entire length of both sequences. Best for comparing sequences of similar length and overall similarity (e.g., Needleman-Wunsch algorithm).
- **Local alignment**: Identifies regions of high similarity within sequences, ignoring poorly matching regions. Optimal for finding conserved domains or motifs within larger, divergent sequences (e.g., Smith-Waterman algorithm used in BLAST).

**ii) E-value explanation (2 marks):**
The E-value represents the **expected number of alignments with scores equal to or better than the observed score that would occur by chance alone** in a database of the given size. Lower E-values indicate more significant matches (e.g., E-value of 1e-10 means this alignment quality would occur by chance only once in 10 billion random searches).

### Question 3 [8 marks]
**i) Importance of TA sites (1 mark):**
TA sites are the **specific insertion sites for Himar-1 transposon**, which only inserts at TA dinucleotide sequences in the genome.

**ii) bioF gene essentiality in lung infection (2 marks):**
The figure suggests **bioF is essential for survival in the lung infection model** because transposon insertions in this gene are depleted in the output pool compared to input, indicating that mutants with disrupted bioF cannot survive or grow in the lung environment.

**iii) bioB gene in tissue culture (2 marks):**
Lower insertion abundance in bioB suggests these mutants have **reduced fitness/growth rate in tissue culture media**. During competitive growth, faster-growing wild-type bacteria outcompete the bioB mutants, leading to their relative depletion in the population.

**iv) Experimental verification for bioB (3 marks):**
1. **Isolate individual bioB mutant clones** from the transposon library
2. **Perform growth curve analysis** comparing bioB mutants vs. wild-type in tissue culture media
3. **Measure doubling times and final cell densities** to quantify the growth defect observed in the pooled experiment

### Question 4 [4 marks]
**i) Why infrequent cutters are used (2 marks):**
Frequent cutters would generate **too many small fragments** that would be difficult to resolve and analyze by PFGE. Infrequent cutters (like NotI, SfiI) create **large fragments (50-1000 kb)** that can be separated by PFGE and provide meaningful fingerprint patterns for strain comparison.

**ii) Different results between fingerprinting and WGS (2 marks):**
- **DNA fingerprinting** has limited resolution and may miss small genetic differences between closely related strains
- **Whole genome sequencing** can detect single nucleotide differences, small indels, and precise genetic relationships, potentially revealing transmission links not detected by fingerprinting or distinguishing strains that appear identical by PFGE.

### Question 5 [1 mark]
The negative control tests for **background ampicillin resistance** in the untransformed bacterial population, ensuring that any ampicillin-resistant colonies observed in transformation experiments are due to successful uptake of the resistance plasmid rather than pre-existing resistance.

---

## Class Test 1 - 2018

### Dr Donaldson's Section [10 marks]

**a) How microarray experiments work (4 marks):**
1. **RNA extraction** from samples and reverse transcription to labeled cDNA (often with fluorescent dyes)
2. **Hybridization** of labeled cDNA to complementary oligonucleotide probes spotted on a solid surface (chip)
3. **Detection** of hybridization signals using fluorescence scanning
4. **Quantification** of signal intensity to determine relative gene expression levels between samples

**b) P-value adjustment importance (1 mark):**
Multiple testing correction is essential to **control false discovery rate** when testing thousands of genes simultaneously, preventing inflation of Type I errors.

**c) Heat map representation (1 mark):**
**Columns**: Different time points (5w, 7w, 9w, 23w)
**Rows**: Individual genes within each cluster

**d) Cluster identification (1 mark):**
**Cluster 4** - shows the pattern of initial decrease followed by increase at 23w.

**e) Correlation analysis power (2 marks):**
Correlation analysis across **diverse experimental conditions** captures functional relationships that persist across different biological contexts, providing stronger evidence for **co-regulation and functional similarity** than clustering within a single experiment.

**f) Functional analysis type (1 mark):**
**Gene Ontology (GO) enrichment analysis** or pathway analysis to identify overrepresented biological processes.

### Dr O'Ryan's Section [20 marks]

**1a) De novo CNV definition (2 marks):**
A **de novo CNV** is a copy number variant present in an individual but **absent in both parents**, indicating it arose as a new mutation. Authors distinguished these by **comparing CNV profiles of affected individuals with their unaffected parents** using family-based analysis.

**1b) Control group features (2 marks):**
Three unique features:
1. **Large sample size** from population-based studies
2. **Matched ethnicity/ancestry** to reduce population stratification
3. **Screened for psychiatric disorders** to ensure true controls

**1c) CNV to postsynaptic signaling flow diagram (3 marks):**
```
CNV identification → Gene mapping → Pathway analysis → Postsynaptic signaling hypothesis
```

**2) ASD genetic variants (4 marks):**
1. **Rare de novo mutations**: High-penetrance variants (CNVs, single gene mutations) found in ~10-15% of cases, revealing critical neurodevelopmental pathways
2. **Common variants**: Low-effect-size polymorphisms identified through GWAS, contributing to polygenic risk and highlighting synaptic function and neuronal development

**3a) GWAS conclusions (2 marks):**
- **Novel locus at 10q24.32** shows genome-wide significant association with ASD
- **Genetic overlap with schizophrenia** suggests shared biological pathways between these psychiatric disorders

**3b) Data shortcoming (2 marks):**
**Limited ancestral diversity** - study focused on European ancestry, limiting generalizability to other populations and potentially missing population-specific risk variants.

**4) Multiple Choice Answers:**
4.1: **c)** Non-random association of alleles in a population or linkage between alleles in families
4.2: **a)** hybridized to short oligonucleotides that correspond to unique DNA flanking specific SNPs
4.3: **a)** rare mutations with high penetrance
4.4: **d)** epigenetic
4.5: **b)** SNPs that identify haplotypes

---

## Class Test 1 - 2020

### Section A - Dr Williams

**Question 1:**
**1) Yellow phenotype explanation (2 marks):**
The M. smegmatis colonies gained the yellow phenotype because they now **contain M. marinum chromosomal DNA** (via the library plasmids) that **encodes the genes responsible for light-induced yellow pigment production**.

**2a) Fragments conferring yellow phenotype (1 mark):**
Fragments **B, C, and D** conferred the yellow phenotype.

**2b) Responsible ORFs (2 marks):**
The data suggests **ORF2 is essential** for the phenotype since:
- Fragment A (containing only ORF1) = white
- Fragments B, C, D (all containing ORF2) = yellow
- ORF3 appears dispensable as fragment D (lacking ORF3) still produces yellow

**2c) Sanger sequencing and BLAST approach (3 marks):**
1. **Sequence ORF2** using Sanger sequencing with gene-specific primers
2. **BLAST search** the ORF2 sequence against Synechocystis genome database
3. **Analyze homologous genes** in Synechocystis for similar function and pathway involvement in pigment biosynthesis

**2d) E-value explanation (2 marks):**
E-value represents the **probability of finding an alignment with equal or better score by chance**. A **good hit would have E-value < 1e-5**, indicating the similarity is highly unlikely to occur randomly.

**Question 2:**
**1a) geneX essentiality (3 marks):**
The absence of transposon insertions in geneX suggests it is **essential for M. tuberculosis survival** under the tested conditions. Since Himar1 can insert at any TA site randomly, the lack of viable mutants indicates that **disruption of geneX is lethal**, preventing recovery of these mutants.

**2) Shuttle vector requirements (2 marks):**
- **Origins of replication** functional in both E. coli and M. tuberculosis
- **Selectable markers** (antibiotic resistance genes) that work in both organisms

### Section B - Dr O'Ryan

**Question 1 (6 marks):**
**a) De novo CNV definition (1 mark):**
A **copy number variant present in an individual but absent in both parents**, representing a new mutation event.

**b) Gene set enrichment meaning (2 marks):**
**Overrepresentation of specific biological pathways or functional gene categories** in schizophrenia-associated CNVs compared to what would be expected by chance, indicating these pathways are preferentially disrupted in the disorder.

**c) Workflow diagram (3 marks):**
```
CNV identification → Gene mapping → Pathway analysis → Biological interpretation
```

**Question 2 (9 marks):**
**a) Type 1 vs Type 2 diabetes differences (2 marks):**
- **Type 1 diabetes**: Shows **stronger, more discrete signals** with higher significance peaks
- **Type 2 diabetes**: Shows **more distributed signals** across the genome with generally lower significance levels

**b) Population structure control (2 marks):**
Authors used **principal component analysis (PCA)** and included **population structure covariates** in their statistical models to account for genetic ancestry differences across the geographic area.

**c) Replication study purpose (2 marks):**
A **replication study validates initial findings** in an independent population to **confirm associations and reduce false positive results**, increasing confidence in the genetic associations.

**d) Correlation between figures (3 marks):**
The data shows **poor correlation** between the British and African studies. This suggests **population-specific genetic architecture** for Type 1 diabetes, with different risk alleles or effect sizes in African vs European populations, highlighting the importance of diverse population studies.

---

## Final Exam - 2015

### Dr O'Ryan Section [17 marks]

**Question 1 (6 marks):**
**a) Incomplete penetrance and screening (2 marks):**
**Incomplete penetrance**: Not all individuals carrying a disease-causing mutation develop the disease (e.g., BRCA1 mutations - ~70% lifetime breast cancer risk). **Screening**: Systematic testing of at-risk individuals enables early detection and intervention.

**b) SSRs and disease (2 marks):**
**Short Tandem Repeats (SSRs)**: Repetitive DNA sequences that can expand, causing disease when exceeding normal ranges (e.g., Huntington's disease - CAG repeat expansion in HTT gene causes neurodegeneration).

**c) CNVs and mental disease (2 marks):**
**Copy Number Variants**: Deletions/duplications associated with psychiatric disorders (e.g., 22q11.2 deletion syndrome causing schizophrenia and developmental delays through dosage-sensitive gene disruption).

**Question 2 (11 marks):**
**a) GWAS explanation (3 marks):**
**Genome-Wide Association Study** compares allele frequencies of hundreds of thousands to millions of SNPs across the genome between cases and controls to identify genetic variants associated with disease susceptibility.

**b) GWAS disadvantages (4 marks):**
1. **Population stratification**: Ancestry differences can create false associations
2. **Missing heritability**: Common variants explain only small fraction of genetic risk
3. **Linkage disequilibrium**: Associated SNPs may not be causal variants
4. **Limited functional insight**: Statistical associations don't reveal biological mechanisms

**c) Figure representation (2 marks):**
**Manhattan plots** showing **-log10(p-values)** for SNP associations across chromosomes, with horizontal line indicating genome-wide significance threshold.

**d) Type 1 vs Type 2 diabetes conclusions (2 marks):**
**Type 1 diabetes** shows **stronger, more localized signals** (higher peaks), while **Type 2 diabetes** shows **more distributed, weaker associations** across multiple loci, reflecting different genetic architectures.

### Dr S. Murray Section [30 marks]

**Question 3 (6 marks):**
**a) Database definitions (3 marks):**
- **CCDS**: Consensus Coding Sequence - high-confidence protein-coding sequences
- **UniProt**: Universal Protein Resource - comprehensive protein sequence and functional information
- **RefSeq**: Reference Sequence - curated, non-redundant sequence database

**b) High-standard transcripts (3 marks):**
Transcripts with **CCDS identifiers** and **complete CDS annotations** represent the highest standard, as they've undergone rigorous curation and validation across multiple databases.

**Question 4 (4 marks):**
**a) GO terms definition (2 marks):**
**Gene Ontology terms** provide standardized vocabulary describing gene function across three domains: **molecular function, biological process, and cellular component**. They enable systematic functional annotation and comparative analysis.

**b) Most confident term (2 marks):**
**IDA (Inferred by Direct Assay)** is more reliable than **ISS (Inferred by Sequence Similarity)** because it's based on **experimental evidence** rather than computational prediction.

**Question 5 (10 marks):**
**a) Log2 fold change utility (3 marks):**
Log2 transformation makes fold changes **symmetric around zero** (upregulation = positive, downregulation = negative), enables **easier statistical analysis**, and **compresses dynamic range** for better visualization.

**b) P-value adjustment necessity (2 marks):**
**Benjamini-Hochberg FDR correction** controls false discovery rate when testing multiple genes simultaneously, **preventing inflation of Type I errors** in large-scale expression studies.

**c) Non-differentially expressed gene (1 mark):**
**Zm.13239.3.S1_at** (adjusted p-value = 0.08 > 0.05) is not significantly differentially expressed.

**d) RNA-Seq additional information (4 marks):**
- **Novel transcript discovery** and splice variant identification
- **Allele-specific expression** analysis
- **Single nucleotide variant detection** in transcribed regions
- **Quantification without prior sequence knowledge**

**Question 6 (10 marks):**
**a) FPKM definition and reasons (3 marks):**
**Fragments Per Kilobase of transcript per Million mapped reads**. Calculated to:
1. **Normalize for transcript length** (longer genes produce more reads)
2. **Normalize for sequencing depth** (enable comparison between samples)

**b) Tophat vs Cufflinks difference (2 marks):**
- **Tophat**: **Read alignment** tool that maps RNA-seq reads to reference genome
- **Cufflinks**: **Transcript assembly and quantification** tool that builds transcript models and estimates expression

**c) Novel transcript analysis (5 marks):**
**i) Structure conclusion (2 marks):** The novel transcript **TCONS_00067999** appears to have **different exon structure** (likely additional exons or alternative splicing) compared to the reference annotation. This could result from **tissue-specific splicing** or **incomplete reference annotation**.

**ii) Validation approach (3 marks):**
1. **RT-PCR with transcript-specific primers** spanning predicted exon junctions
2. **Northern blot analysis** to confirm transcript size
3. **Sanger sequencing** of RT-PCR products to verify exact sequence

### Assoc Prof Reid Section [18 marks]

**Question 7 (7 marks):**
**Quorum sensing in Pseudomonas virulence:**
1. **Density-dependent regulation**: At high cell density, accumulation of autoinducers (3-oxo-C12-HSL, C4-HSL) activates LasR and RhlR regulators
2. **Virulence factor production**: Activated regulators induce expression of elastase, alkaline protease, pyocyanin, and rhamnolipids
3. **Biofilm formation**: Coordinated production of extracellular matrix components
4. **Experimental evidence**: lasR and rhlR mutants show reduced virulence in animal models; complementation restores pathogenicity

**Question 8 (5 marks):**
**i) MS2 A protein regulation (2 marks):**
A protein acts as **translational repressor** of its own mRNA by binding to ribosome binding site, creating **negative feedback loop** that maintains optimal protein levels.

**ii) B. subtilis sigma factor replacement (3 marks):**
1. **Proteolytic degradation** of previous sigma factor
2. **Competitive binding** - higher affinity/concentration of new sigma factor
3. **Anti-sigma factor sequestration** - specific inhibitors remove previous sigma factors

**Question 9 (6 marks):**
**B. subtilis tyrosine operon effector:**
**Effector molecule**: **Tyrosine** itself acts as the inducer
**Experimental proof**: 
1. **In vitro binding studies** showed tyrosine binding to TyrR regulator
2. **Deletion analysis** of tyrR gene eliminated tyrosine-responsive regulation
3. **Complementation experiments** restored tyrosine induction with wild-type tyrR

### Dr M.S. Rafudeen Section [15 marks]

**Question 10 (3 marks):**
**Lambda antisense RNA control:**
**cII antisense RNA** (transcribed from opposite strand) **base-pairs with cII mRNA**, preventing its translation. This **blocks CI repressor synthesis**, ensuring the lytic pathway cannot be established and **maintaining lysogenic commitment**.

**Question 11 (2 marks):**
**Lambda would NOT complete lytic cycle successfully**. Overexpressed integrase would **promote integration of lambda DNA** into the host chromosome, **forcing the lysogenic pathway** even under conditions normally favoring lysis.

**Question 13 (5 marks):**
**Lambda UV induction mechanism:**
```
UV stress → RecA activation → RecA* formation → CI cleavage stimulation → 
CI degradation → Operator derepression → Lytic gene expression
```

**Question 14 (5 marks):**
**CRISPR/Cas9 interference mechanism:**
1. **crRNA:tracrRNA complex** formation with target complementarity
2. **Cas9 recruitment** and DNA binding at PAM-adjacent sites
3. **DNA cleavage** by Cas9 nuclease domains (RuvC and HNH)
4. **Target gene disruption** through double-strand break formation
5. **Cellular repair** often introduces indels, disrupting gene function

### Assoc Prof Vernon Coyne Section [20 marks]

**Question 15 (10 marks):**
**Metagenomic binning:**
**Definition**: Computational process of **grouping contigs/reads** from the same organism based on sequence composition and abundance patterns.

**Usage**: 
- **Genome reconstruction** from mixed communities
- **Taxonomic assignment** of environmental sequences
- **Functional analysis** of uncultured organisms

**Challenges**:
- **Strain-level variation** complicates binning
- **Horizontal gene transfer** creates compositional anomalies
- **Low-abundance organisms** may be missed
- **Contamination** between bins reduces quality

**Question 16 (10 marks):**
**a) INSeq technique basis (6 marks):**
1. **Random transposon mutagenesis** creates insertion library
2. **Selection pressure** applied (gut colonization vs. laboratory growth)
3. **Transposon junction sequencing** identifies insertion sites
4. **Quantitative comparison** of mutant abundance between input and output
5. **Statistical analysis** identifies fitness-affecting genes
6. **Validation** of candidate genes through targeted studies

**b) BT0618 conclusion (4 marks):**
**BT0618 is required for gut colonization** but not laboratory growth. The **significant reduction in output:input ratios** in monoassociated mice (but not chemostat) indicates this gene provides **fitness advantage specifically in the gut environment**, possibly through **sodium/energy metabolism** adaptation to intestinal conditions.

---

## Final Exam - 2019

### Dr L. Donaldson Section [12 marks]

**Question 1 (4 marks):**
**a) Aligned nucleotides (1 mark):** **385 nucleotides** (98% of 393 = 385.14)

**b) Identical nucleotides (1 mark):** **331 nucleotides** (86% of 385 = 331.1)

**c) E-value interpretation (1 mark):** E-value represents **probability of finding alignment by chance**. **4e-75 indicates extremely significant match** - this alignment quality would occur by chance only 4 times in 10^75 random searches.

**d) Max vs Total score difference (1 mark):** **Multiple high-scoring segment pairs (HSPs)** contribute to total score, while max score represents the **single best-scoring alignment segment**.

**Question 2 (8 marks):**
**a) R-value definition (1 mark):** **Pearson correlation coefficient** measuring linear relationship strength between gene expression profiles (-1 to +1).

**b) Gene Ontology definition (2 marks):** **Standardized vocabulary** describing gene function across three domains: **molecular function** (biochemical activity), **biological process** (cellular/organismal functions), and **cellular component** (subcellular localization).

**c) Enrichment explanation (3 marks):** **Enrichment** means the ECGG contains **more auxin-responsive genes than expected by chance** (24/391 vs 282/27906 background). The **highly significant p-value (5.1e-09)** indicates this overrepresentation is **statistically meaningful**, suggesting coordinated regulation.

**d) Hypothesis and test (2 marks):** 
**Hypothesis**: Unknown gene is involved in **auxin signaling/response pathway**
**Experimental test**: **Auxin treatment experiment** - measure unknown gene expression in response to auxin application using qRT-PCR or RNA-seq

### Dr O'Ryan Section [26 marks]

**Question 3 (17 marks):**
**a) De novo CNV identification (2 marks):** Authors ensured de novo status by **comparing CNV profiles between affected individuals and both parents**, including only variants **absent in parental genomes**. They controlled for artifacts by using **matched control populations** and **technical replication**.

**b) Figure interpretation (5 marks):** The figure shows **GO term enrichment** among schizophrenia risk genes. **-log10P(adjusted)** represents **statistical significance** of enrichment (higher values = more significant). The data shows **synaptic transmission and signaling pathways** are most significantly enriched, **consistent with CNV findings** from Kirov et al. that also implicated **postsynaptic signaling**, demonstrating **convergent evidence** across different variant types.

**c) Table interpretation (5 marks):** The table shows **KEGG pathway enrichment** among risk genes. **Seemingly unrelated pathways** (Tuberculosis, Cancer) likely contain **genes with pleiotropic functions** or **shared molecular mechanisms** (e.g., immune signaling, cell cycle regulation) that also affect **neurodevelopment**. This reflects the **interconnected nature** of biological pathways and **shared genetic architecture** across diseases.

**d) Flow diagram (5 marks):**
```
Sample collection → Genotyping → Quality control → Association testing → 
Gene mapping → Pathway analysis → Biological interpretation
```

**Question 4 (9 marks):**
**a) DNA methylation first step (2 marks):** **Bisulfite conversion** - treatment with sodium bisulfite **converts unmethylated cytosines to uracil** while **leaving methylated cytosines unchanged**, enabling methylation detection through sequencing.

**b) ASD pathway classes and hypothesis (5 marks):** 
**Two pathway classes**:
1. **Neurodevelopmental pathways** (synaptic function, neuronal differentiation)
2. **Immune/inflammatory pathways** (cytokine signaling, immune response)

**Hypothesis**: **Immune dysregulation during critical neurodevelopmental periods** contributes to ASD pathogenesis through **neuroinflammation affecting synaptic development**.

**c) Supporting epigenetic evidence (2 marks):** **Histone modification studies** showing altered H3K4me3 and H3K27me3 marks at neurodevelopmental genes, and **microRNA expression changes** affecting immune and synaptic genes.

### Dr P. Meyers Section [22 marks]

**Question 5 (3 marks):**
The reaction depicts **phosphopantetheinyl transferase (PPTase) activity** - transfer of **4'-phosphopantetheine group** from Coenzyme A to **apo-carrier protein**, converting it to **holo-carrier protein** essential for NRPS/PKS function.

**Question 6 (12 marks):**
**a) Amino acid substrates (2 marks):** **Four amino acid substrates** are used, as evidenced by **four A domains** (one per module) in the NRPS system.

**b) Tryptophan isomer (3 marks):** **D-tryptophan** is incorporated by the second module. The **E domain (epimerization)** converts L-tryptophan to D-tryptophan, as indicated by the domain architecture.

**c) TE domain (2 marks):** **Thioesterase domain** - catalyzes **peptide release** from the NRPS through **hydrolysis or cyclization** of the thioester bond.

**d) KisN and KisO proteins (2 marks):** **KisN**: **Flavin reductase** - provides reduced flavin cofactor. **KisO**: **Tryptophan halogenase** - works with KisN to **chlorinate tryptophan residues**.

**e) Tyrosine incorporation (3 marks):** **Three tyrosine molecules** are incorporated into kistamicin structure. **Two are L-tyrosine** (modules 1 and 4 lack E domains), **one is D-tyrosine** (module 3 has E domain for epimerization).

**Question 7 (7 marks):**
**a) EpoB substrate (1 mark):** **Cysteine** is the substrate for the A domain in EpoB.

**b) Cy domain (4 marks):** **Cyclization domain** with functions:
1. **Intramolecular cyclization** of the growing peptide chain
2. **Thioether bond formation** between cysteine and preceding residue
3. **Structural constraint** introduction for bioactivity
4. **Release preparation** for subsequent processing

**c) Ox domain function (2 marks):** **Oxidation domain** - catalyzes **oxidation of thioether to disulfide bond** or **other oxidative modifications** essential for epothilone bioactivity.

### Prof E. Rybicki Section [20 marks]

**Question 8 (5 marks):**
**Virus separation methods:**
1. **Size filtration** (0.22 μm filters) - removes bacteria and larger microbes
2. **Cesium chloride density gradient** - separates viruses based on density
3. **DNase/RNase treatment** - removes free nucleic acids
4. **Chloroform treatment** - removes enveloped viruses if needed
**Rationale**: Each step **selectively enriches viral particles** while removing contaminating cellular material and other microbes.

**Question 9 (10 marks):**
**Marine virus diversity discovery:**
**Advance**: Discovery of **massive viral diversity in ocean ecosystems** through metagenomic sequencing revealed:
- **10^31 viral particles** in oceans globally
- **Novel viral families** with no cultured representatives
- **Auxiliary metabolic genes** in viral genomes affecting host metabolism
- **Viral-host interaction networks** structuring microbial communities

**Impact**: Viruses act as **major drivers of biogeochemical cycles** through host lysis, **horizontal gene transfer vectors**, and **population control agents** maintaining microbial diversity.

**Question 10 (5 marks):**
**Ocean virus implications:**
1. **Biogeochemical cycling** - viral lysis releases nutrients and organic matter
2. **Carbon pump regulation** - affects CO2 sequestration through microbial mortality
3. **Evolutionary pressure** - drives microbial adaptation and diversity
4. **Ecosystem stability** - prevents any single microbial species from dominating
5. **Climate regulation** - influences global carbon and nutrient cycles

### Dr M.S Rafudeen Section [40 marks]

**Question 11 (2 marks):**
**Negatively affect plasmid segregation**. ParM requires **ATP hydrolysis for dynamic instability** - the mutant protein would form **stable, non-dynamic filaments** that cannot undergo the **search-and-capture mechanism** essential for finding and segregating plasmids.

**Question 12 (12 marks):**
**a) 398 vs 441 bp difference (2 marks):** The **398 bp isoform lacks 43 nucleotides** from the 5' end, likely due to **RNase processing** that removes the **mok open reading frame** region.

**b) No 361 bp in lane 1 vs lane 3 (2 marks):** Lane 1 (no treatment) shows **sok RNA inhibition** preventing hok mRNA processing to mature form. Lane 3 (anti-sok PNA) shows **sok RNA sequestration** by PNA, **relieving inhibition** and allowing **hok mRNA maturation** to 361 bp active form.

**c) mok translation confirmation (5 marks):** **Yes, the Northern blot confirms mok won't be translated**. The **398 and 361 bp isoforms lack the 5' region** containing the **mok start codon and ribosome binding site**. Only the **441 bp initial transcript** contains intact mok sequences, but this form is **rapidly processed** and **doesn't accumulate**, preventing mok translation.

**d) Rifampicin and persistence (3 marks):** **Yes, rifampicin may promote persistence**. Lane 2 shows **accumulation of mature 361 bp hok mRNA** similar to anti-sok PNA treatment, suggesting rifampicin **disrupts sok RNA function**, potentially through **transcriptional stress responses** that **activate toxin-antitoxin systems**.

**Question 13 (4 marks):**
**No, the strain will not become persistent**. **RecA is required** for TisB toxin activation through **LexA cleavage**. Without RecA, **LexA remains active** and **represses tisAB expression** even under the constitutive promoter, preventing **TisB-mediated persistence**.

**Question 14 (8 marks):**
**Negative Homeostatic Control Loop:**
```
Environmental stress → Toxin activation → Growth inhibition → 
Reduced antitoxin synthesis → Increased toxin activity → Persistence → 
Stress relief → Antitoxin restoration → Growth resumption
```
**Gene deletion**: **Antitoxin gene deletion** would prevent persistence by eliminating the **regulatory balance** needed for **controlled toxin activation**.

**Question 15 (3 marks):**
**RhlR establishes virulence** by **quorum sensing-dependent activation** of virulence factors including **rhamnolipids, pyocyanin, and elastase** at high cell density, enabling **coordinated pathogenesis** and **biofilm formation** during infection.

**Question 16 (6 marks):**
**a) Xis degradation (1 mark):** **Yes**, Xis is degraded by **both Lon and FtsH proteases** - protein persists longer in mutant strains.

**b) Missing experimental data (1 mark):** **Single lon mutant data** (without ftsH mutation) to determine **individual protease contributions**.

**c) Xis control importance (1 mark):** **Prevents inappropriate excision** - excess Xis could cause **premature prophage excision** during lysogeny.

**d) Int/xis regulation mechanism (3 marks):** **Delayed early transcription** uses **N-mediated antitermination** allowing **read-through transcription** of int and xis genes. **Temporal control** ensures **integration occurs before excision functions** are available, **stabilizing lysogeny**.

**Question 17 (5 marks):**
**a) RNA-OUT function (2 marks):** **Antisense RNA** that **base-pairs with transposase mRNA**, **inhibiting translation** and **controlling transposition frequency**.

**b) Mutation consequences (3 marks):** **Reduced RNA-OUT half-life** leads to **decreased inhibition** of transposase translation, resulting in **increased transposition frequency** and potentially **harmful levels** of transposon activity that could **damage the host genome**.
