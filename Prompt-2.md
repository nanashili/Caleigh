## Final Exam 2021 - MCB3025F

### Section 1: Regulation of Transcription Initiation (30 marks)

**Question 1 [3 marks]**
**General model for nucleosome-free region generation by transcription activators:**

When a sequence-specific transcription activator binds to chromatin, it initiates a cascade of events leading to nucleosome displacement:

1. **Initial binding**: Activator binds to its recognition sequence, potentially competing with nucleosome positioning
2. **Recruitment of chromatin remodeling complexes**: Activators recruit ATP-dependent chromatin remodeling complexes (SWI/SNF family)
3. **Nucleosome sliding/ejection**: Remodeling complexes use ATP hydrolysis to slide or eject nucleosomes
4. **Histone modification**: Recruitment of histone acetyltransferases (HATs) and other modifying enzymes
5. **Stabilization**: Pioneer transcription factors and other proteins stabilize the open chromatin state

**Enzymatic activities involved:**
- **ATP-dependent chromatin remodeling** (SWI/SNF, ISWI, CHD families)
- **Histone acetyltransferase activity** (HATs like CBP/p300)
- **Histone deacetylase activity** (context-dependent)

**Question 2 [4 marks]**
**RACE analysis for transcription start site determination:**

**Experimental steps:**
1. **RNA isolation**: Extract total RNA from cells/tissues of interest
2. **Reverse transcription**: Use gene-specific primer (GSP1) complementary to known sequence downstream of expected TSS
3. **Terminal transferase treatment**: Add homopolymer tail (usually dC) to 3' end of cDNA
4. **PCR amplification**: Use anchor primer complementary to added tail and nested gene-specific primer (GSP2)
5. **Cloning and sequencing**: Clone PCR products and sequence multiple clones

**TSS identification:**
The transcription start site nucleotide is identified by aligning the 5' ends of sequenced cDNA clones with the genomic sequence. The most 5' nucleotide common to multiple independent clones represents the primary TSS.

**Question 3 [6 marks]**

**a) Relative luciferase activity [2 marks]:**
Relative luciferase activity represents the ratio of firefly luciferase activity (from experimental promoter) to Renilla luciferase activity (internal control). This normalization controls for:
- Transfection efficiency variations
- Cell number differences
- General transcriptional/translational variations
- Allows comparison between different promoter constructs

**b) Effect of SRE mutation [2 marks]:**
Mutating SRE transcription factor binding sites dramatically reduces FASN promoter activity (approximately 80-90% reduction compared to wildtype). This indicates SRE sites are critical positive regulatory elements for FASN transcription.

**c) Necessity and sufficiency of SRE sites [2 marks]:**
SRE sites appear **necessary** but **not sufficient** for optimal FASN promoter activity. Evidence:
- **Necessary**: SRE mutation severely impairs activity
- **Not sufficient**: NF-Y mutation also reduces activity, indicating multiple elements are required for optimal function
- Optimal activity requires both SRE and NF-Y sites working synergistically

**Question 4 [3 marks]**

**a) Minor groove base pair discrimination [1 mark]:**
**A-T and T-A base pairs** can be discriminated through direct readout via the minor groove. The minor groove of A-T pairs has distinct hydrogen bonding patterns and electrostatic properties compared to G-C pairs.

**b) Protein dimer type [2 marks]:**
This protein would be a **homodimer**. Reasoning:
- The sequence AGGTCANNNTGACCT shows **palindromic symmetry**
- AGGTCA and TGACCT are reverse complements
- Palindromic recognition sequences typically indicate homodimeric proteins
- Each subunit recognizes one half-site in a symmetrical manner

**Question 5 [6 marks]**

**a) Determining apparent Kd [4 marks]:**
**Approach**: Electrophoretic Mobility Shift Assay (EMSA) with varying protein concentrations

**Experimental conditions:**
1. **Fixed DNA concentration**: Use constant, low concentration of labeled DNA probe
2. **Protein titration**: Vary protein concentration across wide range (typically 10-fold dilutions)
3. **Binding conditions**: Physiological salt, pH 7.4, room temperature, 15-30 min incubation
4. **Analysis**: Quantify bound vs. free DNA, plot fraction bound vs. protein concentration
5. **Kd determination**: Protein concentration giving 50% binding equals apparent Kd

**Specific conditions:**
- Low ionic strength buffer (avoid high salt that disrupts binding)
- Include competitor DNA (poly dI-dC) to reduce non-specific binding
- Use protein concentrations spanning 0.1-10x expected Kd

**b) Apparent vs. true Kd [2 marks]:**
- **Apparent Kd**: Measured under specific experimental conditions, may include effects of competing factors, cofactors, or non-specific binding
- **True Kd**: Intrinsic thermodynamic dissociation constant under ideal conditions
- True Kd is typically lower (higher affinity) than apparent Kd due to experimental artifacts and competing interactions

**Question 6 [5 marks]**

**a) FLAG-epitope tagging approach [3 marks]:**
1. **Transfection**: Introduce FLAG-tagged protein expression vector into cultured cells
2. **Cell lysis**: Harvest cells and prepare nuclear/whole cell extracts
3. **Immunoprecipitation**: Use anti-FLAG antibody to capture tagged protein and associated factors
4. **Analysis**: Analyze co-precipitated proteins by SDS-PAGE and mass spectrometry

**Control**: Cells transfected with empty vector or untagged protein to identify non-specific interactions with FLAG antibody/resin.

**b) Multi-protein complex investigation [2 marks]:**
**Method**: Size exclusion chromatography (gel filtration)
**Procedure**: Apply FLAG-purified complexes to calibrated gel filtration column
**Outcomes**: 
- Large complexes elute early (high molecular weight fractions)
- Individual proteins elute later (lower molecular weight)
- Compare apparent molecular weight to sum of individual protein masses

**Question 7 [3 marks]**

**a) Fold purification calculation [2 marks]:**
Initial specific activity = 8 mg His-protein / 3500 mg total = 0.00229 mg/mg
Final specific activity = 1.8 µg/µL / 3 µg/µL = 0.6 mg/mg
**Fold purification = 0.6 / 0.00229 = 262-fold**

**b) Percent yield calculation [1 mark]:**
Initial His-protein = 8 mg
Final His-protein = 2.5 mL × 2.4 µg/µL = 6 mg
**% Yield = (6/8) × 100 = 75%**

### Section 2: Steroid Receptors (20 marks)

**Question 8 [9 marks]**

**a) Ki calculation for NET-A [2 marks]:**
From graph: IC50 for NET-A ≈ 10 nM
Using Ki = IC50 / [1 + ([ligand]/Kd)]
Ki = 10 nM / [1 + (0.2 nM/1 nM)] = 10 nM / 1.2 = **8.33 nM**

**b) Biocharacter determination [3 marks]:**
**Cannot definitively determine biocharacter** from competitive binding alone. Competitive binding only measures:
- Binding affinity (Ki values)
- Ability to compete for receptor binding

**Examples:**
- MIB, DHT: Known agonists, compete effectively
- RU486: Known antagonist, also competes effectively
- Both agonists and antagonists can have high binding affinity

**Additional assays needed**: Functional assays (transcriptional activation, conformational changes) required to determine agonist vs. antagonist activity.

**c) Questions answered by competitive binding [2 marks]:**
1. **Binding affinity**: Relative Ki values for different ligands
2. **Receptor selectivity**: Cross-reactivity with related receptors
3. **Structure-activity relationships**: How chemical modifications affect binding
4. **Ligand specificity**: Whether compounds bind to target receptor

**d) Importance for structure-function [2 marks]:**
- **Drug development**: Identifies high-affinity compounds for further testing
- **Receptor characterization**: Maps ligand-binding domain requirements
- **Evolutionary relationships**: Compares binding across receptor subtypes
- **Therapeutic targeting**: Guides design of selective modulators

**Question 9 [5 marks]**
**Experimental design for serine residue requirement:**

1. **Construct generation**: Create PR mutant with serine→alanine substitution
2. **Transfection**: Co-transfect cells with:
   - Wild-type or mutant PR expression vectors
   - PRE-luciferase reporter construct
   - Renilla luciferase control
3. **Treatment**: Incubate with/without progesterone
4. **Analysis**: Measure dual luciferase activity

**Controls needed:**
- Empty vector (no PR)
- Reporter alone (no PR)
- Vehicle treatment (no hormone)
- Positive control (known functional PR)

**Expected outcome**: If serine is required, mutant PR will show reduced transactivation compared to wild-type.

**Question 10 [6 marks]**
**Protein-protein interactions illustrating general cell signaling concepts:**

1. **Signal amplification**: Single receptor recruits multiple coactivators
2. **Signal integration**: Multiple signaling pathways converge on shared cofactors
3. **Specificity**: Different receptors recruit distinct cofactor combinations
4. **Regulation**: Post-translational modifications control interaction strength
5. **Cooperativity**: Cofactor binding enhances additional recruitment
6. **Context dependence**: Cell-type specific cofactor expression determines response

### Section 3: Marine Proteomics (30 marks)

**Question 11 [6 marks]**
**Why trypsin is preferred in bottom-up proteomics:**

1. **Predictable cleavage**: Cleaves specifically after basic residues (Arg, Lys)
2. **Optimal peptide size**: Generates peptides ideal for MS analysis (8-25 amino acids)
3. **High activity**: Robust enzyme, works under denaturing conditions
4. **Commercial availability**: Highly pure, modified forms available
5. **Database compatibility**: Protein databases optimized for tryptic peptides
6. **Ionization efficiency**: Tryptic peptides ionize well due to basic C-terminus

**Question 12 [4 marks]**
**Mass calculations:**

Peptide A++: m/z = 376, charge = +2
**Mass of A = (376 × 2) - 2 = 750 Da**

Peptide B+: m/z = 522, charge = +1  
**Mass of B = 522 - 1 = 521 Da**

**Question 13 [10 marks]**
**Shotgun proteomics workflow for heat stress in oysters:**

1. **Sample preparation**: 
   - Collect oysters at control and elevated temperatures
   - Extract proteins using lysis buffer with protease inhibitors
   - **Reason**: Preserve protein integrity and prevent degradation

2. **Protein digestion**:
   - Reduce disulfide bonds, alkylate cysteines
   - Digest with trypsin overnight
   - **Reason**: Generate peptides suitable for MS analysis

3. **Peptide separation**:
   - LC separation using reverse-phase chromatography
   - **Reason**: Reduce sample complexity, improve detection

4. **Mass spectrometry**:
   - ESI-MS/MS analysis with data-dependent acquisition
   - **Reason**: Identify and quantify peptides

5. **Data analysis**:
   - Database searching against oyster proteome
   - Statistical analysis for differential expression
   - **Reason**: Identify heat stress-responsive proteins

6. **Validation**:
   - Western blotting or targeted MS for key proteins
   - **Reason**: Confirm shotgun results

**Question 14 [10 marks]**
**Peptide sequence determination from MS/MS spectrum:**

Working backwards from y-ions:
- y1 = 156.1 (R)
- y2 = 253.15 (P + R) 
- y3 = 350.2 (T + P + R)
- y4 = 481.24 (M + T + P + R)
- y5 = 578.29 (P + M + T + P + R)

**Possible sequence: KPMTPR** (considering two missed cleavages at K and R positions)

### Section 4: Protein Processing and Trafficking (40 marks)

**Question 15 [2 marks]**
**Susceptibility to reducing agents:**

**Vasopressin (CYFQNCPRG)** would be susceptible to beta-mercaptoethanol because it contains **two cysteine residues** (positions 1 and 6) that likely form an intramolecular disulfide bond. Reducing agents break disulfide bonds.

Oxytocin contains only one cysteine, so cannot form disulfide bonds within the peptide.

**Question 16 [12 marks]**

**a) N-glycosylation site prediction [1 mark]:**
The sites are predicted based on the **N-X-S/T consensus sequence** where X is any amino acid except proline. All four sites (153, 172, 223, 354) follow this pattern and are located in the extracellular/luminal domain.

**b) Experimental design for high mannose glycans [9 marks]:**

**Experiment**: Endoglycosidase H (EndoH) sensitivity assay

**Procedure**:
1. Transfect HEK 293 cells with human BACE expression vector
2. Harvest cells and prepare lysates
3. Treat lysates with EndoH or PNGase F
4. Analyze by Western blotting with anti-BACE antibody

**Controls**:
- Untreated lysate
- Mock-transfected cells
- Known high mannose protein (positive control)
- Known complex glycan protein (negative control for EndoH)

**Expected results**:
- **High mannose only**: EndoH treatment reduces molecular weight
- **Complex glycans present**: EndoH has no effect, PNGase F reduces molecular weight
- **Prediction confirmed**: EndoH sensitivity indicates high mannose glycans

**c) Transmembrane domain sequences [2 marks]:**
- **Signal sequence**: N-terminal hydrophobic sequence for ER targeting
- **Stop-transfer sequence**: Hydrophobic transmembrane domain that halts translocation
- **Function**: Signal sequence initiates translocation, stop-transfer anchors protein in membrane

**Question 17 [4 marks]**

**a) EPO length discrepancy [2 marks]:**
- **Signal peptide**: Removed during ER processing (typically 20-30 amino acids)
- **Propeptide sequences**: Additional processing sequences removed
- **5' and 3' UTRs**: Non-coding regions in mRNA not translated
- **Alternative splicing**: Some sequences may be spliced out

**b) SRP degradation effect [2 marks]:**
EPO translocation would be **severely impaired**. SRP is essential for:
- Recognizing signal sequences during translation
- Targeting ribosomes to ER membrane
- **Result**: EPO would be synthesized in cytosol and likely degraded rather than properly processed and secreted

**Question 18 [2 marks]**
**(i) N-terminal signal sequence deleted**: **Cytoplasmic localization** - protein cannot enter ER
**(ii) KDEL sequence deleted**: **Secreted** - protein enters ER but cannot be retained, follows default secretory pathway

**Question 19 [5 marks]**

**a) Type of transport [1 mark]:**
**Retrograde transport** (endosome → Golgi → ER → cytoplasm)

**b) COP I requirements [4 marks]:**
**(i) Uptake**: **COP I not required** - endocytosis uses clathrin
**(ii) Endosome to Golgi**: **COP I not required** - uses different machinery
**(iii) Golgi to ER**: **COP I required** - COP I mediates retrograde Golgi-ER transport
**(iv) ER to cytoplasm**: **Unknown** - experiment doesn't distinguish this step

**Question 20 [8 marks]**

**a) Immunoprecipitation principle [1 mark]:**
Antibody specifically binds target protein, antibody-protein complex captured on protein A/G beads, associated proteins co-precipitate.

**b) Rationale for both antibodies [2 marks]:**
- **Anti-PRC17**: Confirms PRC17 is present in immunoprecipitate
- **Anti-Rab5**: Demonstrates Rab5 co-precipitates with PRC17, indicating interaction

**c) Interaction conclusions [2 marks]:**
**Wild-type PRC17** interacts with Rab5 (both proteins detected). **Mutant PRC17** shows reduced interaction, suggesting GAP activity is important for stable complex formation.

**d) Endogenous PRC17 [1 mark]:**
**No** - no PRC17 detected in control lanes, indicating HEK 293 cells don't express detectable endogenous PRC17.

**e) Rab5 inactivation [2 marks]:**
Target **GEF (Guanine nucleotide Exchange Factor)** for inhibition. GEFs promote GDP→GTP exchange; inhibiting GEF would lock Rab5 in inactive GDP-bound state.

**Question 21 [7 marks]**

**a) Cysteine to tyrosine mutation effect [2 marks]:**
**Reduced lipid-linked proteins** because cysteine residues are required for:
- **Palmitoylation**: Cysteine forms thioester bonds with fatty acids
- **Prenylation**: Some prenylation sites require nearby cysteines
- Loss of cysteine eliminates these lipid modification sites

**b) Ubiquitin fate determination [1 mark]:**
**(ii) The particular lysine amino acid used to link ubiquitin molecules** - K48 linkages target for degradation, K63 linkages for signaling

**c) M-6-P receptor release [1 mark]:**
**Low pH in late endosomes** causes conformational change in M-6-P receptor, reducing affinity for M-6-P and releasing lysosomal enzymes.

**d) Amyloidogenic pathway enhancement [3 marks]:**
**No, ACE 2 does not enhance amyloidogenic pathway**. Evidence:
- **C99 levels decreased** in ACE 2 vs. WT (amyloidogenic pathway reduced)
- **P3 levels increased** in ACE 2 vs. WT (non-amyloidogenic pathway enhanced)
- ACE 2 mutation shifts processing toward non-amyloidogenic pathway

---

## Final Exam 2022 - MCB3025F

### Section 1: Regulation of Transcription Initiation (30 marks)

**Question 1 [3 marks]**

**a) Regulatory promoter features [2 marks]:**
- **Transcription factor binding sites**: Specific DNA sequences recognized by activators
- **Modular organization**: Multiple binding sites work cooperatively
- **Position independence**: Can function at various distances from core promoter
- **Orientation independence**: Function regardless of orientation

**b) Transcription activation mechanism [1 mark]:**
Activators **recruit cofactors and RNA polymerase II machinery** to promoters, facilitating pre-initiation complex formation and transcription initiation.

**Question 2 [3 marks]**
**Primer extension analysis:**

1. **RNA isolation**: Extract total RNA from cells
2. **Primer annealing**: Anneal 32P-labeled gene-specific primer downstream of expected TSS
3. **Reverse transcription**: Extend primer with reverse transcriptase
4. **Gel electrophoresis**: Run products on sequencing gel alongside sequencing ladder
5. **Autoradiography**: Visualize extended products

**TSS identification**: The 5' end of the extended product (determined by comparison with sequencing ladder) indicates the transcription start site nucleotide.

**Question 3 [5 marks]**

**a) Dual luciferase principle [2 marks]:**
**Principle**: Two different luciferases (firefly and Renilla) with distinct substrates allow internal normalization.
- **Firefly luciferase**: Reports promoter activity
- **Renilla luciferase**: Internal control for transfection efficiency
- **Relative activity**: Firefly/Renilla ratio normalizes for experimental variation

**b) PR1 and HRE importance [3 marks]:**
Both PR1 and HRE are **necessary but not sufficient**:
- **PR1 necessity**: pGL3-4 (PR1 mutated) shows ~70% reduction in activity
- **HRE necessity**: pGL3-5 (HRE mutated) shows ~60% reduction in activity  
- **Not sufficient**: Neither element alone supports full activity; both required for optimal promoter function
- **Evidence**: pGL3-3 (lacking PR1-PR5) shows dramatic reduction, indicating multiple elements needed

**Question 4 [3 marks]**

**a) DNA bend vs. kink [1 mark]:**
- **Bend**: Gradual curvature over several base pairs
- **Kink**: Sharp, localized distortion at specific base pair

**b) Homodimer vs. heterodimer [2 marks]:**
**Homodimer** - the sequence AGGTCANNNACTGGA contains **imperfect palindromic symmetry**. AGGTCA and TCCAGT are related but not perfect reverse complements, suggesting a homodimer with some flexibility in recognition.

**Question 5 [8 marks]**

**a) Binding affinity comparison [2 marks]:**
Protein has **higher affinity for competitor A** than competitor B. Evidence:
- **Competitor A**: Effectively competes at low concentrations (5-10 fold excess)
- **Competitor B**: Requires much higher concentrations (50-100 fold) for competition
- Higher affinity = lower concentration needed for competition

**b) Binding affinity and Kd relationship [2 marks]:**
**Inverse relationship**: Higher binding affinity corresponds to lower Kd value. Kd represents the concentration needed for 50% binding; lower Kd = higher affinity.

**c) Further competition with competitor B [2 marks]:**
**Yes**, competition would likely occur at 1000-fold excess. The trend shows increasing competition with higher concentrations; sufficient competitor B should eventually displace the protein-DNA complex.

**d) Identifying important nucleotides [2 marks]:**
**Systematic mutagenesis approach**:
1. Create series of competitors with single nucleotide substitutions
2. Test each mutant in EMSA competition
3. **Loss of competition** indicates critical nucleotide for binding
4. **Maintained competition** suggests nucleotide not essential

**Question 6 [5 marks]**

**a) TFIID affinity for activation domains [2 marks]:**
**AD2 > AD1 > AD3** based on protein band intensities:
- **AD2**: Strongest bands for multiple TFIID subunits
- **AD1**: Moderate interaction
- **AD3**: Weakest interaction
- Right panel confirms specific TFIID interaction

**b) Significance of lane 2 [1 mark]:**
**Negative control** - shows background binding to GST alone, demonstrating specificity of activation domain interactions.

**c) Identifying responsible subunits [2 marks]:**
**Individual subunit pulldown**: Express individual TFIID subunits as GST fusions, test binding to each activation domain separately. Compare results to identify which specific subunit(s) mediate the interactions observed.

**Question 7 [3 marks]**

**a) Fold purification [2 marks]:**
Initial specific activity = 10 mg / 5000 mg = 0.002
Final specific activity = 1.8 µg/µL / 2 µg/µL = 0.9
**Fold purification = 0.9 / 0.002 = 450-fold**

**b) Yield calculation [1 mark]:**
Initial protein = 10 mg
Final protein = 4 mL × 1.8 µg/µL = 7.2 mg
**Yield = (7.2/10) × 100 = 72%**

### Section 2: Steroid Receptors (20 marks)

**Question 8 [10 marks]**

**a) Experimental design for ligand affinity [8 marks]:**

**Parameters to determine experimentally**:
1. **IC50 of ligand X**: Concentration causing 50% inhibition of reference binding
2. **Kd of reference agonist**: Known value
3. **Concentration of reference agonist**: Used in assay

**Experimental procedure**:
1. **Cell preparation**: Use receptor-expressing cells
2. **Binding assay**: Incubate cells with fixed concentration of radiolabeled reference agonist
3. **Competition**: Add increasing concentrations of unlabeled ligand X
4. **Measurement**: Determine bound radioactivity at each concentration
5. **Analysis**: Plot % inhibition vs. ligand X concentration, determine IC50

**Equation**: Ki = IC50 / [1 + ([ligand]/Kd)]

**b) Ki calculation [2 marks]:**
Ki = 1 nM / [1 + (1000 pM / 10 nM)] = 1 nM / [1 + 0.1] = **0.91 nM**

**Question 9 [10 marks]**

**a) Co-immunoprecipitation principles [4 marks]:**
1. **Antibody specificity**: Antibody recognizes target protein specifically
2. **Complex capture**: Antibody-protein complexes captured on protein A/G beads
3. **Co-precipitation**: Associated proteins remain bound during washing
4. **Detection**: Western blotting reveals both target and associated proteins

**b) No GRIP-1 in lane 4 [2 marks]:**
**Possible explanations**:
- Mutant GR (3A) cannot interact with GRIP-1 due to loss of serine residues
- Serine residues required for ligand-dependent conformational change
- Phosphorylation of serines may be necessary for GRIP-1 recruitment

**c) GR-GRIP-1 interaction mechanism [2 marks]:**
Results suggest **ligand-dependent and serine-dependent** interaction:
- Interaction requires dexamethasone treatment
- Serine residues (S203, S211, S226) are essential
- Likely involves ligand-induced conformational change and serine phosphorylation

**d) GRIP-1 function and mechanism [2 marks]:**
**Function**: Transcriptional coactivator
**Mechanism**: 
- Bridges steroid receptors to transcriptional machinery
- Possesses histone acetyltransferase activity
- Recruits additional cofactors to enhance transcription

### Section 3: Marine Proteomics (30 marks)

**Question 10 [8 marks]**
**Mass determination from LC-MS spectrum:**

From the spectrum showing peaks at m/z 211, 363, 584, and 725:

**Protein A**: Multiple charge states visible
- If 725 is [M+2H]2+, then mass = (725 × 2) - 2 = 1448 Da
- If 363 is [M+4H]4+, then mass = (363 × 4) - 4 = 1448 Da

**Protein B**: 
- If 584 is [M+2H]2+, then mass = (584 × 2) - 2 = 1166 Da

**Protein C**:
- If 211 is [M+H]+, then mass = 211 - 1 = 210 Da

**Question 11 [15 marks]**
**Peptide sequence from MS/MS spectrum:**

Working from y-ion series (C-terminus to N-terminus):
- y1 = 156.1 (R)
- y2 = 269.18 (L + R) 
- y3 = 382.26 (L + L + R)
- y4 = 495.34 (I + L + L + R)
- y5 = 608.42 (L + I + L + L + R)

**Sequence (N to C terminus): LLILLR** (no missed cleavages, trypsin cleaves after R)

**Question 12 [3 marks]**
**Three reasons for bottom-up approach:**

1. **Technical feasibility**: Smaller peptides easier to analyze by MS than intact proteins
2. **Database searching**: Extensive peptide databases available for identification
3. **Sensitivity**: Better ionization and detection of peptides vs. whole proteins

**Question 13 [4 marks]**
**Trypsin digestion products:**

MALSTRVATSKLICDVTRPASDTASTVEEKYGDAAS

Trypsin cleaves after K and R:
1. **MALSTR** (cleaved after R)
2. **VATSK** (cleaved after K) 
3. **LICDVTR** (cleaved after R)
4. **PASDTASTVEEK** (cleaved after K)
5. **YGDAAS** (C-terminal fragment)

### Section 4: Protein Processing and Trafficking (40 marks)

**Question 14 [11 marks]**

**a) Mutant gp160 glycan analysis [3 marks]:**
**Both high mannose and complex N-glycans** are present. Evidence:
- **EndoH sensitivity**: Partial molecular weight reduction indicates some high mannose glycans
- **PNGase F sensitivity**: Complete deglycosylation indicates N-linked glycans present
- **Incomplete EndoH digestion**: Suggests mixture of high mannose and complex glycans

**b) EndoH vs. PNGase F molecular weight difference [2 marks]:**
EndoH leaves **GlcNAc residue** attached to asparagine, while PNGase F removes entire glycan. This accounts for the slightly higher molecular weight after EndoH treatment.

**c) Precursor oligosaccharide [2 marks]:**
**Dolichol-linked Glc3Man9GlcNAc2** synthesized in cytoplasm. Enters ER lumen via **flippase activity** that translocates the lipid-linked oligosaccharide across ER membrane.

**d) Translocation mode [4 marks]:**
**Co-translational translocation**:
- Signal sequence recognized by SRP during translation
- Ribosome targeted to ER membrane
- Nascent polypeptide threaded through Sec61 translocon
- Glycosylation occurs during translocation in ER lumen

**Question 15 [7 marks]**

**a) Protein X targeting [2 marks]:**
**Mitochondrial targeting** would predominate. Mitochondrial signal sequences are typically recognized first during translation, before nuclear export signals can function. The protein would be imported into mitochondria.

**b) Transmembrane topology [5 marks]:**
**Required sequences**:
- **NTS (N-terminal signal)**: Required for ER entry, places N-terminus in lumen
- **STA (Stop-transfer anchor)**: Required to halt translocation, creates transmembrane domain
- **SA not required**: Would create different topology

**Reasoning**: NTS initiates translocation (N-terminus in lumen), STA stops translocation leaving C-terminus in cytoplasm, creating the observed topology.

**Question 16 [5 marks]**

**a) pH 7.2 effect [2 marks]:**
**Negative effect** - lysosomes require acidic pH (~4.5-5.0) for optimal function:
- Lysosomal enzymes have acidic pH optima
- Neutral pH reduces enzymatic activity
- Impaired protein degradation and cellular clearance

**b) Pathways to late endosomes [3 marks]:**
1. **Endocytosis**: Plasma membrane → early endosomes → late endosomes
2. **Autophagy**: Autophagosomes fuse with late endosomes
3. **Biosynthetic pathway**: Golgi → late endosomes (newly synthesized lysosomal enzymes)

**Question 17 [6 marks]**

**a) COP-I vs. COP-II roles [2 marks]:**
- **COP-II**: ER → Golgi anterograde transport
- **COP-I**: Golgi → ER retrograde transport and intra-Golgi transport

**b) GAP degradation effects [3 marks]:**
**Impaired vesicle transport** because:
- GAPs inactivate Rab GTPases by promoting GTP hydrolysis
- Without GAPs, Rabs remain constitutively active
- **Result**: Vesicles cannot complete transport cycles, leading to transport defects and vesicle accumulation

**c) Dynamin function [1 mark]:**
**Vesicle scission** - dynamin forms collar around vesicle neck and uses GTP hydrolysis to pinch off vesicles from donor membranes.

**Question 18 [6 marks]**

**a) Phosphorylation reversibility [1 mark]:**
**Yes** - phosphorylation is reversible through kinase (addition) and phosphatase (removal) activities.

**b) No band in wt lane [1 mark]:**
Wild-type mice don't express **human huntingtin transgene**; antibody is human-specific.

**c) Re-probing necessity [1 mark]:**
**Loading control** - confirms equal amounts of huntingtin protein loaded, validates phosphorylation differences are not due to protein level differences.

**d) Akt contribution conclusion [3 marks]:**
**YAC 128 shows increased Akt activation** and **increased huntingtin phosphorylation** compared to YAC 18, suggesting **Akt contributes to pathological phosphorylation** of expanded huntingtin. The correlation between polyglutamine expansion, Akt activation, and huntingtin phosphorylation supports Akt's role in disease pathogenesis.

**Question 19 [5 marks]**
**Experimental design for mutant GCA localization:**

1. **Transfection**: Transfect HEK293 cells with mutant GCA expression vector
2. **Immunofluorescence**: Fix cells, permeabilize, incubate with anti-GCA antibody
3. **Co-localization**: Use organelle-specific markers (ER, Golgi, lysosome markers)
4. **Imaging**: Confocal microscopy to determine subcellular localization
5. **Controls**: Wild-type GCA, empty vector, secondary antibody alone

---

## Class Test 1 - 2019 MCB3025F

**Question 1 [3 marks]**
**Properties and functions of promoter elements:**

**Core promoter regions:**
- **Location**: Immediately surrounding TSS (-40 to +40)
- **Function**: Direct RNA polymerase II binding and positioning
- **Elements**: TATA box, Initiator, DPE, CAAT box

**Upstream regulatory regions:**
- **Location**: -100 to -1000 bp from TSS
- **Function**: Bind gene-specific transcription factors
- **Properties**: Modular, orientation-dependent

**Enhancers:**
- **Location**: Variable distance, can be upstream, downstream, or intronic
- **Function**: Increase transcription rate dramatically
- **Properties**: Position and orientation independent, work over long distances

**Question 2 [5 marks]**

**A. RACE methodology [4 marks]:**
1. **RNA extraction**: Isolate total RNA from cells/tissues
2. **Reverse transcription**: Use gene-specific primer (GSP1) to synthesize cDNA
3. **Tailing**: Add homopolymer tail (dC) to 3' end of cDNA using terminal transferase
4. **PCR amplification**: Use anchor primer (complementary to tail) and nested GSP2
5. **Cloning and sequencing**: Clone products and sequence to identify 5' ends

**B. Importance of TSS mapping [1 mark]:**
TSS position is essential for **defining promoter boundaries** and **identifying functional regulatory elements** relative to the transcription start site.

**Question 3 [3 marks]**
**Experimental strategy for cortisone-responsive elements:**

1. **Promoter cloning**: Clone promoter region into luciferase reporter vector
2. **Deletion analysis**: Create series of 5' deletions, test cortisone responsiveness
3. **Mutagenesis**: Introduce point mutations in candidate regulatory sequences
4. **Transfection assays**: Test constructs in cortisone-responsive cell line

**Controls**: 
- Vehicle treatment (no cortisone)
- Empty vector
- Known cortisone-responsive promoter (positive control)

**Question 4 [5 marks]**

**A. DNA binding sequence determination [3 marks]:**
Based on competition data, the protein binds to **GTGGTA** sequence. 

**Reasoning**:
- **Mut1**: Changes AGA AAT → gtg ggc, **competes poorly** (lane 3)
- **Mut2**: Changes GTG GTA → tct tcg, **competes poorly** (lane 4)  
- **Mut3-5**: Mutations outside this region **compete effectively**
- Core sequence GTGGTA is critical for binding

**B. Testing TEST1 involvement [2 marks]:**
**Supershift assay**: Add anti-TEST1 antibody to EMSA reaction.
**Outcomes**:
- **Supershift**: Confirms TEST1 in complex (higher molecular weight band)
- **No effect**: TEST1 not involved
- **Competition**: Pre-incubate with TEST1 peptide as control

**Question 5 [2 marks]**
**Base readout vs. shape readout:**

**Base readout**: Direct recognition of specific nucleotides through hydrogen bonding and van der Waals contacts in major groove.

**Shape readout**: Recognition of DNA structural features (bending, widening) rather than specific base sequences.

**Question 6 [3 marks]**

**A. Nucleosome architecture [2 marks]:**
**Components**: Histone octamer (2 each of H2A, H2B, H3, H4) + 147 bp DNA wrapped 1.65 times
**Architecture**: DNA wrapped around histone core, connected by linker DNA with H1 histone

**B. Histone acetylation correlation [1 mark]:**
**Positive correlation**: Higher acetylation associated with **increased gene activity** due to chromatin relaxation and enhanced transcription factor access.

**Question 7 [4 marks]**

**A. Spt3-interacting proteins [2 marks]:**
**Proteins B and C** appear to interact with Spt3. Evidence:
- Present in FLAG-HA-Spt3 lane (lane 2) but absent/reduced in control (lane 1)
- Protein A present in both lanes (non-specific)
- Protein D unclear from data shown

**B. Direct interaction testing [2 marks]:**
**GST pulldown assay**: Express individual proteins as GST fusions, test direct binding to purified Spt3 protein in vitro. This eliminates indirect interactions through bridging proteins.

**Question 8 [3 marks]**

**A. Fold purification [1 mark]:**
Initial specific activity = 1 mg / 1000 mg = 0.001
Final specific activity = 0.45 µg/µL / 0.5 µg/µL = 0.9
**Fold purification = 0.9 / 0.001 = 900-fold**

**B. Purity calculation [1 mark]:**
**Purity = (0.45/0.5) × 100 = 90%**

**C. Yield calculation [1 mark]:**
Initial: 1 mg, Final: 1 mL × 0.45 µg/µL = 0.45 mg
**Yield = (0.45/1) × 100 = 45%**

**Question 9 [2 marks]**
**Cryo-EM advantages over X-ray crystallography:**

1. **No crystallization required**: Can study proteins in native-like states
2. **Large complex analysis**: Better suited for analyzing large, flexible complexes
3. **Dynamic information**: Can capture different conformational states
4. **Membrane protein analysis**: Easier to study membrane-embedded proteins

---

## Final Exam 2017 - Sections 1 & 2

### Section 1: Regulation of Transcription Initiation (30 marks)

**Question 1 [3 marks]**
**Proximal regulatory regions vs. enhancers:**

**Similarities**:
- Both bind sequence-specific transcription factors
- Both regulate gene expression
- Both can be targets for regulatory signals

**Differences**:
**Proximal regions**: 
- Located near promoter (-200 to -1000 bp)
- Position-dependent function
- Often orientation-dependent

**Enhancers**:
- Distance-independent (can be >100 kb away)
- Position-independent (upstream, downstream, intronic)
- Orientation-independent
- Stronger activation potential

**Question 2 [2 marks]**
**Nucleosome core particle organization:**

**Structure**: 147 bp DNA wrapped 1.65 times around histone octamer core containing two copies each of histones H2A, H2B, H3, and H4. Histone tails extend outward and are subject to post-translational modifications.

**Question 3 [4 marks]**

**a) RACE experimental steps [3 marks]:**
1. **RNA isolation and reverse transcription**: Use gene-specific primer
2. **Tailing reaction**: Add poly(dC) tail to cDNA 3' end
3. **PCR amplification**: Use anchor primer and nested gene-specific primer
4. **Cloning and sequencing**: Determine 5' end sequences

**b) TSS identification [1 mark]:**
The **5' nucleotide** common to multiple independent cDNA clones represents the transcription start site when aligned with genomic sequence.

**Question 4 [7 marks]**

**a) Identifying GH-responsive elements [5 marks]:**

**Experimental approach**:
1. **Promoter cloning**: Clone upstream region into reporter vector
2. **Deletion analysis**: Create 5' and internal deletions
3. **Transfection assays**: Test GH responsiveness in appropriate cell line
4. **DNase I footprinting**: Identify protein-protected regions
5. **EMSA**: Test specific sequence binding

**Potential outcomes**:
- **Deletion removes responsiveness**: Contains essential element
- **Footprint protection**: Indicates protein binding site
- **EMSA complex formation**: Confirms sequence-specific binding

**b) Testing sufficiency [2 marks]:**
**Heterologous promoter assay**: Clone identified element upstream of minimal promoter (e.g., TATA box only) and test GH responsiveness. Sufficient elements will confer GH responsiveness to minimal promoter.

**Question 5 [6 marks]**

**a) Protein binding sequence [2 marks]:**
Based on competition data, protein binds to **AAAGAAGGGGTA** sequence.

**Reasoning**: Mutations 1, 2, 7, 8, 9 show reduced competition, indicating these regions are critical for binding. The core sequence spans the region where these mutations cluster.

**b) Consensus sequence determination [1 mark]:**
**Systematic mutagenesis**: Create all possible single nucleotide substitutions within the binding site and test each by EMSA competition to determine optimal sequence at each position.

**c) Testing Mark3 involvement [3 marks]:**

**Experiment**: Supershift assay with anti-Mark3 antibody
**Outcomes**:
- **Supershift**: Confirms Mark3 in complex
- **Loss of complex**: Mark3 required for binding
- **No change**: Mark3 not involved

**Controls**:
- Pre-immune serum
- Irrelevant antibody
- Mark3 peptide competition

**Question 6 [2 marks]**
**X-ray crystallography definitions:**

**Unit cell**: Smallest repeating unit of crystal lattice containing complete asymmetric unit(s)
**Asymmetric unit**: Smallest portion of crystal structure from which entire crystal can be generated by symmetry operations

**Question 7 [6 marks]**

**a) Evidence for multi-protein complexes [3 marks]:**
- **Panel A**: Multiple proteins co-purify specifically with FLAG-hSRB10
- **Panel B**: Proteins co-elute in high molecular weight fractions (>669 kDa)
- **Size exclusion**: Complex size much larger than individual protein masses
- **Stoichiometry**: Multiple subunits present in similar amounts

**b) Significance of lane 1 [1 mark]:**
**Negative control** demonstrating specificity - shows background binding to anti-FLAG resin without FLAG-tagged protein.

**c) Direct interaction testing [2 marks]:**
**Binary interaction assay**: Express individual proteins recombinantly, test direct binding between purified hSRB10 and each candidate protein using techniques like surface plasmon resonance or isothermal titration calorimetry.

### Section 2: Steroid Receptors (40 marks)

**Question 8 [4 marks]**

**a) Importance of fractional occupancy [2 marks]:**
Fractional occupancy determines **magnitude of biological response**. It represents the proportion of receptors bound by ligand, which directly correlates with downstream signaling intensity and physiological effects.

**b) Determinants of fractional occupancy [2 marks]:**
- **Ligand concentration [L]**
- **Receptor affinity (Kd)**
- **Relationship**: Fractional occupancy = [L]/(Kd + [L])

**Question 9 [2 marks]**
**Ki calculation:**
Using Cheng-Prusoff equation:
Ki = IC50 / [1 + ([L]/Kd)] = 1 nM / [1 + (1000 pM/10 nM)] = 1 nM / 1.1 = **0.91 nM**

**Question 10 [4 marks]**
**Volume calculation:**
Required: 2 pmoles
Stock: 3 × 10⁸ dpm/L, specific activity: 4 × 10⁸ dpm/nmole

2 pmoles = 0.002 nmoles
Required dpm = 0.002 × 4 × 10⁸ = 8 × 10⁵ dpm
Volume = 8 × 10⁵ dpm / 3 × 10⁸ dpm/L = **2.67 × 10⁻³ L = 2.67 mL**

**Question 11 [7 marks]**
**Questions answered by dose-response analysis:**

1. **Potency (EC50)**: Concentration producing half-maximal response
2. **Efficacy**: Maximum response achievable
3. **Receptor reserve**: Relationship between occupancy and response
4. **Cooperativity**: Hill coefficient indicates binding cooperativity
5. **Partial agonism**: Submaximal efficacy compared to full agonists
6. **Antagonist characterization**: Competitive vs. non-competitive inhibition
7. **Therapeutic window**: Separation between desired and toxic effects

**Question 12 [5 marks]**
**GR conformational changes and transcriptional activity:**

**Agonist binding**:
- Induces conformational change exposing AF-2 domain
- Creates coactivator binding surface
- Allows recruitment of transcriptional machinery
- Results in transcriptional activation

**Antagonist binding**:
- Different conformational change
- AF-2 domain remains buried/altered
- Cannot recruit coactivators effectively
- May recruit corepressors instead
- Results in transcriptional repression/no activation

**Question 13 [6 marks]**

**a) Four protein classes [2 marks]:**
1. **Corepressors** (NCOR, SMRT)
2. **Chromatin remodeling complexes** (SWI/SNF)
3. **General transcription factors** (TFIID, TFIIB)
4. **Kinases/phosphatases** (regulatory enzymes)

**b) Co-IP vs. two-hybrid comparison [4 marks]:**

**Co-IP advantages**:
- Studies interactions in native cellular context
- Detects indirect interactions through bridging proteins
- Physiologically relevant conditions

**Co-IP disadvantages**:
- Cannot distinguish direct vs. indirect interactions
- May miss weak/transient interactions
- Requires specific antibodies

**Two-hybrid advantages**:
- Detects direct protein-protein interactions
- High sensitivity for weak interactions
- No antibody requirement

**Two-hybrid disadvantages**:
- Artificial nuclear environment
- May detect non-physiological interactions
- Cannot study membrane proteins easily

**Question 14 [12 marks]**

**a) ChIP experimental steps [5 marks]:**
1. **Cross-linking**: Formaldehyde treatment to fix protein-DNA interactions
2. **Cell lysis and sonication**: Fragment chromatin to 200-500 bp
3. **Immunoprecipitation**: Use specific antibody to capture protein-DNA complexes
4. **Reverse cross-linking**: Remove formaldehyde, purify DNA
5. **PCR analysis**: Amplify specific promoter regions, quantify enrichment

**b) ChIP controls [3 marks]:**
- **IgG control**: Non-specific antibody
- **Input DNA**: Total chromatin before IP
- **No antibody control**: Background binding
- **Irrelevant promoter**: Specificity control

**c) Repression mechanism [4 marks]:**
Results suggest **active repression mechanism**:
- **GR recruitment**: Rapid binding to promoter
- **Corepressor recruitment**: NCOR1, SMRT recruited with GR
- **HDAC recruitment**: HDAC2/3 recruited for chromatin modification
- **SRC-1 exclusion**: Coactivator not recruited
- **Temporal pattern**: Sequential recruitment suggests coordinated repression complex assembly

---

## Final Exam 2018 - Sections 1 & 2

### Section 1: Regulation of Transcription Initiation (30 marks)

**Question 1 [3 marks]**
**Roles in gene activation:**

**General transcription factors (GTFs)**:
- Form pre-initiation complex at core promoter
- Position RNA polymerase II correctly
- Enable basal transcription

**Gene-specific activators**:
- Bind to regulatory sequences (enhancers/promoters)
- Recruit cofactors and GTFs
- Increase transcription rate above basal level

**Mediator complex**:
- Bridges activators and GTFs
- Transmits regulatory signals to transcription machinery
- Integrates multiple regulatory inputs

**Question 2 [3 marks]**

**a) Chromatin types [2 marks]:**
**Euchromatin**: 
- Loosely packed, transcriptionally active
- Enriched in active histone marks

**Heterochromatin**:
- Tightly packed, transcriptionally silent
- Enriched in repressive histone marks

**b) Histone acetylation effects [1 mark]:**
**Chromatin relaxation** - acetylation neutralizes positive charges on histone tails, reducing DNA-histone interactions and making chromatin more accessible to transcription factors.

**Question 3 [4 marks]**
**RACE methodology:**

1. **RNA extraction**: Isolate total RNA from cells
2. **Reverse transcription**: Use gene-specific primer to synthesize cDNA
3. **Tailing**: Add poly(dC) tail using terminal transferase
4. **PCR**: Amplify using anchor primer and nested gene-specific primer
5. **Sequencing**: Determine 5' end sequences to identify TSS

**TSS identification**: The 5' nucleotide where multiple cDNA sequences begin represents the transcription start site.

**Question 4 [4 marks]**

**a) Renilla luciferase purpose [1 mark]:**
**Internal control** for transfection efficiency and general cellular effects, allowing normalization of firefly luciferase data.

**b) Most important binding sites [2 marks]:**
**SRE sites** appear most important. Evidence:
- Construct (b) lacking SRE shows dramatic reduction (~90%) in activity
- Construct (c) lacking NF-Y shows moderate reduction (~50%)
- SRE removal has greater impact than NF-Y removal

**c) Necessity and sufficiency [1 mark]:**
SRE sites are **necessary but not sufficient**. While essential for high activity, NF-Y sites also contribute significantly, indicating multiple elements required for optimal function.

**Question 5 [2 marks]**

**a) Common DNA-binding motifs [1 mark]:**
- **Helix-turn-helix**
- **Zinc finger**
- **Basic helix-loop-helix**
- **Leucine zipper**

**b) Common structural principle [1 mark]:**
**Alpha helices** fit into major groove of DNA, allowing specific amino acid-base contacts for sequence recognition.

**Question 6 [4 marks]**

**a) Binding affinity comparison [2 marks]:**
Protein has **higher affinity for competitor A** than competitor B. Competitor A effectively competes at lower concentrations (5-10 fold excess), while competitor B requires much higher concentrations (50-100 fold) for equivalent competition.

**b) Identifying important nucleotides [2 marks]:**
**Systematic mutagenesis**: Create competitors with single nucleotide substitutions throughout the binding site. Test each mutant - **loss of competition** indicates critical nucleotides for protein binding.

**Question 7 [3 marks]**
**FLAG-epitope tagging approach:**

1. **Transfection**: Express FLAG-tagged protein in cultured cells
2. **Cell lysis**: Prepare nuclear or whole cell extracts
3. **Immunoprecipitation**: Use anti-FLAG antibody to capture tagged protein
4. **Analysis**: Identify co-precipitated proteins by mass spectrometry

**Control**: Cells transfected with empty vector or untagged protein to identify non-specific interactions.

**Question 8 [4 marks]**

**a) TFIID affinity ranking [2 marks]:**
**AD2 > AD1 > AD3** based on band intensities in both panels. AD2 shows strongest interaction with multiple TFIID subunits, AD1 shows moderate interaction, AD3 shows weakest interaction.

**b) Identifying responsible subunits [2 marks]:**
**Individual subunit analysis**: Express each TFIID subunit separately as GST fusion, test binding to each activation domain. Compare results to identify which specific subunit(s) mediate the observed interactions.

**Question 9 [3 marks]**
**X-ray crystallography workflow:**

1. **Protein purification**: Obtain pure, homogeneous protein
2. **Crystallization**: Screen conditions to grow diffraction-quality crystals
3. **Data collection**: Expose crystals to X-rays, collect diffraction patterns
4. **Phase determination**: Solve phase problem (molecular replacement, heavy atoms)
5. **Model building**: Build atomic model into electron density
6. **Refinement**: Optimize model against experimental data

### Section 2: Protein Processing and Trafficking (40 marks)

**Question 10 [12 marks]**

**a) EndoH vs. PNGase F activity [2 marks]:**
**EndoH**: Cleaves only high mannose and hybrid N-glycans (ER/early Golgi forms)
**PNGase F**: Cleaves all N-linked glycans regardless of processing state

**b) N260Q notation [1 mark]:**
**Asparagine at position 260 mutated to glutamine** - eliminates N-glycosylation site by removing asparagine from N-X-S/T consensus sequence.

**c) Glycan analysis [3 marks]:**
**Wild-type**: Has both high mannose and complex glycans (partial EndoH sensitivity, complete PNGase F sensitivity)
**N260Q mutant**: Has fewer glycosylation sites due to mutation, remaining sites show similar processing pattern

**d) Molecular weight difference [2 marks]:**
EndoH leaves **single GlcNAc residue** attached to asparagine, while PNGase F removes entire glycan including this residue, resulting in lower molecular weight.

**e) Translocation mode [4 marks]:**
**Co-translational translocation**:
- Signal sequence recognized by SRP during translation
- Ribosome targeted to ER membrane
- Nascent chain threaded through Sec61 translocon
- N-glycosylation occurs during translocation in ER lumen

**Question 11 [2 marks]**
**Protein X targeting:**
**Mitochondrial targeting** would predominate. Mitochondrial signal sequences are typically recognized during translation before nuclear export signals can function, directing the protein to mitochondria.

**Question 12 [5 marks]**

**Ran-GTP/GDP cycle [2 marks]:**
- **Nucleus**: RanGEF (RCC1) promotes GDP→GTP exchange
- **Cytoplasm**: RanGAP promotes GTP→GDP hydrolysis
- **Gradient**: High Ran-GTP in nucleus, high Ran-GDP in cytoplasm

**GDP-locked Ran effects [3 marks]:**
**Impaired nuclear transport**:
- **Import**: Reduced cargo release in nucleus
- **Export**: Reduced cargo binding in nucleus
- **Overall**: Bidirectional transport severely compromised

**Question 13 [4 marks]**
**Cysteine to serine mutation effects:**

**Affected modifications**:
1. **Disulfide bond formation**: Loss of cysteine eliminates disulfide bonds
2. **Palmitoylation**: Cannot form thioester bonds with fatty acids
3. **Metal coordination**: Reduced ability to coordinate metal ions
4. **Protein stability**: Altered folding and stability

**Question 14 [3 marks]**
**Anfinsen's second experiment:**

**Result**: **Inactive ribonuclease** with incorrect disulfide bonds
**Reason**: Removing mercaptoethanol in presence of urea allows disulfide formation while protein is still denatured
**Conclusion**: **Correct folding requires proper sequence** - disulfide bonds must form after correct tertiary structure is established

**Question 15 [5 marks]**

**a) Sandwich ELISA principles [2 marks]:**
1. **Capture antibody**: Immobilized on plate, binds target protein
2. **Sample addition**: Target protein binds to capture antibody
3. **Detection antibody**: Binds different epitope on target protein
4. **Signal generation**: Enzyme-linked detection system produces measurable signal

**b) BACE 2 pathway analysis [3 marks]:**
**BACE 2 does NOT enhance amyloidogenic pathway**:
- **C99 (amyloidogenic)**: Decreased in BACE 2 vs. WT BACE
- **C83 (non-amyloidogenic)**: Increased in BACE 2 vs. WT BACE
- **Conclusion**: BACE 2 mutation shifts processing toward non-amyloidogenic pathway

**Question 16 [4 marks]**

**a) pH 7.2 effects [2 marks]:**
**Negative effect** - lysosomes require acidic pH (~4.5) for optimal enzyme function. Neutral pH impairs:
- Lysosomal enzyme activity
- Protein degradation
- Cellular clearance functions

**b) Endocytic pathways [2 marks]:**
**Late endosomes** receive material via:
1. **Endocytosis**: Plasma membrane → early endosomes → late endosomes
2. **Autophagy**: Autophagosomes → late endosomes
3. **Biosynthetic**: Trans-Golgi network → late endosomes

**Question 17 [5 marks]**

**a) Phosphorylation reversibility [1 mark]:**
**Yes** - reversible through kinase (phosphorylation) and phosphatase (dephosphorylation) activities.

**b) No wt band explanation [1 mark]:**
**Human-specific antibody** - wild-type mice don't express human huntingtin transgene.

**c) Re-probing necessity [1 mark]:**
**Loading control** - ensures equal protein amounts, validates that phosphorylation differences reflect true changes rather than loading variations.

**d) Akt contribution [2 marks]:**
**Q111 shows increased Akt activation and huntingtin phosphorylation** compared to Q7, suggesting **Akt contributes to pathological phosphorylation** in polyglutamine-expanded huntingtin, potentially contributing to neurodegeneration.
