## Class Test 3 - [MC3025F] 2017

### Question 1 [3 marks]
Protein X has both high mannose and complex N-glycans based on the gel electrophoresis results.

Reasons:
- In the untreated (UD) lane, Protein X appears at 200 kDa, representing the fully glycosylated form.
- After EndoH treatment (E), Protein X shows a partial shift to 130 kDa, indicating that only some glycans were cleaved. EndoH specifically cleaves high mannose and some hybrid N-glycans but not complex N-glycans.
- After PNGase F treatment (P), Protein X shifts completely to 80 kDa, representing the fully deglycosylated form. PNGase F cleaves all types of N-glycans.
- The partial resistance to EndoH but complete sensitivity to PNGase F confirms that Protein X contains both high mannose glycans (EndoH-sensitive) and complex N-glycans (EndoH-resistant).

### Question 2 [2 marks]
Histone proteins are targeted to the nucleus through the nuclear import pathway:

1. Histones contain nuclear localization signals (NLSs) - short amino acid sequences rich in basic residues (lysine/arginine).
2. In the cytosol, importin-α recognizes and binds to the NLS on the histone proteins.
3. Importin-β then binds to the importin-α/histone complex.
4. This trimeric complex docks at the nuclear pore complex (NPC).
5. The complex translocates through the NPC into the nucleus via interactions with nucleoporins.
6. In the nucleus, RanGTP binds to importin-β, causing the complex to dissociate.
7. The histones are released into the nucleoplasm where they can associate with newly synthesized DNA.

### Question 3 [4 marks]
Experiment to determine if mutant GCA accumulates in the ER:

**Approach**: Fluorescence microscopy with co-localization analysis

**Rationale**: This approach allows direct visualization of the subcellular localization of the mutant GCA protein and confirmation of its retention in the ER through co-localization with known ER markers.

**Main steps**:
1. Generate expression constructs for:
   - Wild-type GCA tagged with GFP (GCA-WT-GFP)
   - Mutant GCA tagged with GFP (GCA-Mut-GFP)
2. Transfect HeLa cells with either construct
3. After 24-48 hours, fix cells and immunostain with antibodies against ER markers (e.g., calnexin or PDI)
4. Perform confocal microscopy to visualize:
   - GFP signal (GCA localization)
   - ER marker signal
   - DAPI staining (nuclei)
5. Analyze co-localization between GCA-GFP and ER markers
6. Compare the degree of ER co-localization between wild-type and mutant GCA

### Question 4 [8 marks]

#### a) Activation states of Rab GTPases and their regulation [2 marks]
Rab GTPases cycle between two states:
- **Active state**: GTP-bound form, which can interact with effector proteins
- **Inactive state**: GDP-bound form, which cannot interact with effectors

This cycling is regulated by:
- **GEFs** (Guanine nucleotide Exchange Factors): Promote exchange of GDP for GTP, activating Rab
- **GAPs** (GTPase-Activating Proteins): Enhance the intrinsic GTPase activity of Rab, accelerating GTP hydrolysis to GDP, inactivating Rab
- **GDIs** (GDP Dissociation Inhibitors): Maintain Rab in the inactive GDP-bound state and regulate membrane association/dissociation

#### b) Rationale for using immunoprecipitation and western blot [2 marks]
The combination of immunoprecipitation (IP) and western blot analysis provides:

1. **Specificity**: IP with anti-Rab5 antibody selectively isolates Rab5 and any proteins physically associated with it from the complex cell lysate.
2. **Sensitivity**: Western blot with anti-PRC17 antibody provides sensitive detection of PRC17 specifically among the immunoprecipitated proteins.
3. **Verification of interaction**: This approach confirms a physical interaction between Rab5 and PRC17 in cellular context.
4. **Quantification**: The intensity of the PRC17 band in the western blot indicates the relative amount of PRC17 associated with Rab5.

#### c) Explanation of experimental results [4 marks]
The western blot shows:

1. **3T3-pLPC lane**: Control cells transfected with empty vector show no PRC17 band after Rab5 immunoprecipitation, confirming specificity of the antibody.

2. **3T3-PRC17wt lane**: Cells expressing wild-type PRC17 show a strong PRC17 band after Rab5 immunoprecipitation, demonstrating that wild-type PRC17 physically interacts with Rab5 in cells.

3. **3T3-PRC17mut lane**: Cells expressing mutant PRC17 (with substitutions that abolish GAP activity) show a PRC17 band of similar intensity to wild-type, indicating that:
   - The mutations that abolish GAP activity do not affect the physical binding of PRC17 to Rab5
   - The interaction between PRC17 and Rab5 is independent of PRC17's GAP activity
   - PRC17 likely has a binding domain for Rab5 that is separate from its catalytic GAP domain

4. **Rab5 bands**: Similar intensity across all lanes confirms equal immunoprecipitation efficiency, validating the comparison of PRC17 binding.

### Question 5 [5 marks]

#### a) Function of pH regulation in the endocytic pathway [1 mark]
pH regulation in the endocytic pathway is crucial for:
- Proper sorting and trafficking of cargo proteins
- Ligand-receptor dissociation in early endosomes
- Activation of hydrolytic enzymes in late endosomes and lysosomes
- Protein conformational changes required for membrane fusion events
- Regulation of protein-protein interactions specific to different compartments

#### b) Reasoning behind using transferrin linked to two different fluorophores [3 marks]
The dual-fluorophore approach enables precise pH measurement because:

1. **pH-dependent fluorescence**: Fluorescein fluorescence is highly pH-sensitive (decreases at acidic pH), while Alexa Fluor 546 fluorescence is pH-stable.

2. **Ratiometric measurement**: By calculating the ratio of fluorescein/Alexa Fluor 546 intensity, researchers can determine the exact pH value independent of:
   - Variations in transferrin uptake between cells
   - Differences in endosome size or number
   - Photobleaching effects
   - Variations in microscope illumination intensity

3. **Calibration curve**: The ratio values can be converted to absolute pH values using a calibration curve generated by exposing cells to buffers of known pH in the presence of ionophores.

#### c) Conclusion from the results [1 mark]
The authors concluded that NHE6 (Na⁺/H⁺ exchanger 6) is essential for maintaining proper endosomal pH. The significantly more acidic pH in both soma and processes of NHE6-null (MUT) neurons compared to wild-type (WT) demonstrates that NHE6 normally functions to alkalinize endosomal compartments by exchanging luminal protons for cytosolic sodium ions.

### Question 6 [2 marks]
The fundamental difference between "top-down" and "bottom-up" proteomic approaches:

**Top-down proteomics**:
- Analyzes intact proteins without prior digestion
- Proteins are separated and then directly introduced into the mass spectrometer
- Preserves post-translational modifications and protein isoforms
- Provides information about the complete protein structure

**Bottom-up proteomics**:
- Analyzes peptides derived from enzymatically digested proteins
- Proteins are first digested (typically with trypsin) into peptides
- Peptides are separated and analyzed by mass spectrometry
- Protein identification occurs by matching peptide mass spectra to database entries

### Question 7 [2 marks]
Two reasons why the "bottom-up" approach is used most often in proteomic studies:

1. **Higher sensitivity and accuracy**: Peptides are more easily ionized, separated, and detected by mass spectrometers than intact proteins, resulting in better sensitivity, mass accuracy, and resolution.

2. **Better chromatographic separation**: Peptides exhibit better chromatographic behavior than intact proteins, allowing for more efficient separation prior to mass spectrometric analysis.

Additional reasons (any two would be sufficient):
- Simpler data interpretation and analysis
- More robust database search algorithms for peptide identification
- Better compatibility with high-throughput workflows
- Lower sample quantity requirements
- Better detection of low-abundance proteins

### Question 8 [4 marks]
Trypsin cleaves peptide bonds at the C-terminal side of lysine (K) and arginine (R) residues, except when followed by proline (P).

Peptides resulting from trypsin digestion of MALSTRVATSKLICDVTRPASDTASTVEEKYGDAAS:

1. MALSTR
2. VATSK
3. LICDVTR
4. PASDTASTVEEK
5. YGDAAS

## Final Exam - [MCB3025F] 2017 - Section 3

### Question 15 [10 marks]

#### a) Schematic diagram of APP processing [4 marks]
```
                   α-secretase        γ-secretase
                        ↓                 ↓
                        |                 |
N-terminus -------|-------|-------------|------- C-terminus
              β-secretase               γ-secretase
                   ↓                     ↓
                   
Non-amyloidogenic pathway (α → γ):
- sAPPα (N-terminal fragment)
- C83 (C-terminal fragment)
- p3 peptide (small fragment)
- AICD (APP Intracellular Domain)

Amyloidogenic pathway (β → γ):
- sAPPβ (N-terminal fragment)
- C99 (C-terminal fragment)
- Aβ peptide (amyloid beta)
- AICD (APP Intracellular Domain)
```

#### b) Formation and function of disulfide bonds [2 marks]
**Formation**:
- Disulfide bonds form between the thiol (-SH) groups of two cysteine residues
- Formation occurs in the oxidizing environment of the ER lumen
- Catalyzed by protein disulfide isomerase (PDI)
- Requires molecular oxygen as the ultimate electron acceptor

**Function**:
- Stabilize tertiary protein structure by covalently linking distant parts of the polypeptide
- Increase protein stability and resistance to denaturation
- Critical for proper protein folding and maintaining native conformation
- Often essential for protein function, especially in secreted proteins

#### c) Experiment to investigate intracellular localization of BACE [4 marks]
**Approach**: Immunofluorescence microscopy with co-localization analysis

**Materials**:
- HEK 293 cells
- Expression vector for BACE with an epitope tag (e.g., FLAG, HA, or GFP)
- Antibodies against the epitope tag
- Antibodies against organelle markers (e.g., calnexin for ER, GM130 for Golgi, LAMP1 for lysosomes, EEA1 for early endosomes)
- Fluorescently-labeled secondary antibodies
- Confocal microscope

**Procedure**:
1. Transfect HEK 293 cells with the BACE expression construct
2. After 24-48 hours, fix cells with paraformaldehyde
3. Permeabilize cells with detergent (e.g., Triton X-100)
4. Block with BSA or serum
5. Incubate with primary antibodies against:
   - The epitope tag on BACE
   - Various organelle markers
6. Wash and incubate with fluorescently-labeled secondary antibodies
7. Counterstain nuclei with DAPI
8. Mount slides and image using confocal microscopy
9. Analyze co-localization between BACE and organelle markers
10. Quantify the degree of co-localization using appropriate software

**Analysis**:
- Calculate Pearson's correlation coefficient between BACE and each organelle marker
- Determine the predominant subcellular localization based on highest co-localization
- Compare results with known localization patterns of BACE from literature

### Question 16 [5 marks]
**Receptor-mediated endocytosis of LDL particles**:

1. **Recognition and binding**:
   - LDL particles in the bloodstream contain apolipoprotein B-100
   - LDL receptors on the cell surface specifically recognize and bind to apolipoprotein B-100
   - Binding occurs in clathrin-coated pits, which are specialized regions of the plasma membrane

2. **Internalization**:
   - Clathrin triskelions assemble on the cytoplasmic side of the membrane, forming a polyhedral lattice
   - Adaptor proteins (AP-2) link clathrin to the membrane and the receptor's cytoplasmic tail
   - The pit invaginates and pinches off to form a clathrin-coated vesicle containing LDL-receptor complexes
   - Dynamin GTPase mediates the scission of the vesicle from the plasma membrane

3. **Uncoating**:
   - The clathrin coat is rapidly removed by the action of Hsc70 and auxilin
   - The uncoated vesicle fuses with early endosomes

4. **Sorting in early endosomes**:
   - The acidic pH (~6.0) of early endosomes causes LDL to dissociate from its receptor
   - The LDL receptor is sorted into tubular extensions of the early endosome
   - These tubules pinch off to form recycling endosomes

5. **Receptor recycling**:
   - Recycling endosomes return the LDL receptor to the plasma membrane
   - The receptor can undergo multiple rounds of endocytosis

6. **LDL processing**:
   - LDL remains in the vesicular portion of the endosome
   - Early endosomes mature into late endosomes with progressively lower pH
   - Late endosomes fuse with lysosomes
   - Lysosomal hydrolases degrade LDL particles:
     - Cholesteryl esters are hydrolyzed by acid lipase to free cholesterol
     - Apolipoprotein B-100 is degraded to amino acids
   - Free cholesterol exits the lysosome and is used by the cell or stored as cholesteryl esters

### Question 17 [2 marks]
To redirect Protein X from the mitochondrial matrix to the ER lumen:

1. **Remove the mitochondrial targeting sequence** from the N-terminus of Protein X
2. **Add an ER signal sequence** to the N-terminus of Protein X
3. **Add an ER retention signal** (KDEL) to the C-terminus to ensure retention in the ER lumen

Annotated diagram:
```
Normal Protein X:
N-terminus [Mitochondrial targeting sequence]---[Protein X core sequence]---C-terminus

Reconstructed Protein X:
N-terminus [ER signal sequence]---[Protein X core sequence]---[KDEL]---C-terminus
```

The ER signal sequence will direct the protein to the ER translocon during translation, resulting in co-translational import into the ER lumen. The KDEL sequence will ensure the protein is retained in the ER lumen by retrieval from the Golgi apparatus.

### Question 18 [6 marks]

#### a) Difference between EndoH and PNGaseF [2 marks]
**Endoglycosidase H (EndoH)**:
- Cleaves high-mannose and some hybrid N-linked glycans
- Cuts between the two N-acetylglucosamine (GlcNAc) residues in the chitobiose core
- Cannot cleave complex N-glycans that have been processed in the Golgi apparatus
- Leaves one GlcNAc residue attached to the asparagine

**Peptide-N-glycosidase F (PNGaseF)**:
- Cleaves all types of N-linked glycans (high-mannose, hybrid, and complex)
- Cuts between the innermost GlcNAc and the asparagine residue
- Removes the entire glycan structure from the protein
- Converts the asparagine to aspartic acid in the process

#### b) Explanation of results in Panel B [4 marks]
Panel B shows the results of EndoH digestion on glucocerebrosidase (GC) from normal cells (Patient 1) and various Gaucher's disease patients.

1. **Normal cells (Patient 1)**:
   - Show two bands after EndoH treatment: an upper EndoH-resistant band and a lower EndoH-sensitive band
   - This indicates that normal GC exists in two forms: one with complex glycans (processed through the Golgi) and one with high-mannose glycans (not fully processed)
   - The presence of EndoH-resistant GC indicates successful trafficking through the Golgi

2. **Gaucher's disease patients**:
   - Different patterns are observed among patients, reflecting the heterogeneity of the disease
   - Some patients (e.g., 13, 6, 5) show predominantly EndoH-sensitive GC, indicating retention in the ER
   - Other patients (e.g., 4, 10) show a mixture of EndoH-sensitive and resistant forms, suggesting partial trafficking through the Golgi
   - The severity of ER retention correlates with disease phenotype

3. **Control proteins**:
   - Erk shows consistent expression across all samples, confirming equal loading
   - β-HexA shows consistent expression and glycosylation pattern, indicating that the glycosylation defect is specific to GC and not a general glycosylation problem

These results demonstrate that many Gaucher's disease mutations cause ER retention of glucocerebrosidase due to misfolding, preventing proper trafficking to lysosomes where it normally functions.

### Question 19 [7 marks]

#### a) Creation of PRC17 mutant [2 marks]
The authors created the PRC17 mutant by:
1. Identifying conserved amino acids in the catalytic domain that are essential for GAP activity
2. Using site-directed mutagenesis to substitute these amino acids with alanine residues
3. Specifically, they substituted two conserved amino acids that have been shown in previous studies to abolish GAP activity
4. The mutations were designed to disrupt the catalytic function while minimally affecting the overall protein structure

#### b) Assay used for Panel A [3 marks]
The assay used to measure GAP activity in Panel A is a **GTPase activity assay**:

1. **Preparation**:
   - Purified Rab5 protein was loaded with radioactive GTP (γ-³²P-GTP)
   - PRC17 (wild-type or mutant) was purified from transfected cells

2. **Reaction**:
   - Rab5-GTP was incubated with or without PRC17 proteins
   - If PRC17 has GAP activity, it will stimulate Rab5 to hydrolyze GTP to GDP + Pi
   - The released radioactive phosphate (³²Pi) was measured

3. **Quantification**:
   - The amount of released ³²Pi was quantified by scintillation counting
   - Results were expressed as percentage of GTP hydrolyzed (% GAP activity)
   - Control reactions without PRC17 established the baseline intrinsic GTPase activity of Rab5

#### c) Conclusion from Panel B [2 marks]
Panel B shows that PRC17 has selective GAP activity toward specific Rab GTPases:

1. PRC17 significantly reduces the GTPase activity of Rab5 (shown by lower % GAP activity when PRC17 is present)
2. PRC17 has no significant effect on the GTPase activity of Rab4 and Rab11
3. This demonstrates that PRC17 is a specific GAP for Rab5 but not for Rab4 or Rab11
4. The specificity suggests that PRC17 regulates specific vesicular trafficking pathways controlled by Rab5 (early endosome function) but not those controlled by Rab4 or Rab11 (recycling endosomes)

## Final Exam - [MCB3025F] 2017 - Section 4

### Question 20 [10 marks]

#### a) Cumulative oyster mortality [1 mark]
The cumulative oyster mortality by the end of the recovery stage (96h) was 64.18%.

This is calculated from the survival ratio at 96h, which is 35.82%. Therefore, mortality = 100% - 35.82% = 64.18%.

#### b) Reasons for fewer mortalities during second stress phase [9 marks]

1. **Selective survival of heat-resistant individuals**:
   - The first heat stress phase (0-24h) caused high mortality (27.36%)
   - The oysters that survived were likely those with inherently greater heat tolerance
   - This natural selection resulted in a more heat-resistant population facing the second stress phase

2. **Acquired thermotolerance through heat shock protein induction**:
   - The heat map shows significant upregulation of multiple HSP genes (HSP70, HSP90, HSP60, HSP40, HSP20) during the first stress phase
   - These heat shock proteins act as molecular chaperones that:
     - Prevent protein denaturation
     - Assist in refolding damaged proteins
     - Target irreparably damaged proteins for degradation
   - The elevated HSP levels persisted through the recovery period, providing protection during the second stress phase

3. **Cellular memory and epigenetic adaptation**:
   - The first heat exposure triggered epigenetic changes that maintained the expression of protective genes
   - This "cellular memory" of the first stress event enabled faster and more robust responses to the second exposure

4. **Metabolic adjustments**:
   - Surviving oysters likely adjusted their metabolic pathways to cope with heat stress
   - The heat map shows upregulation of ATPase genes, suggesting adaptation of energy metabolism
   - These metabolic adjustments would have been maintained or quickly reactivated during the second stress phase

5. **Activation of anti-apoptotic pathways**:
   - The heat map shows upregulation of anti-apoptotic genes (e.g., Bcl2) and modulation of apoptotic pathways
   - This would protect cells from heat-induced programmed cell death during the second stress phase

6. **Antioxidant defense enhancement**:
   - Heat stress typically increases reactive oxygen species (ROS) production
   - The heat map shows upregulation of GST (glutathione S-transferase), an important antioxidant enzyme
   - Enhanced antioxidant defenses would mitigate oxidative damage during the second stress phase

7. **Recovery period benefits**:
   - The 48h recovery period allowed for:
     - Repair of sub-lethal cellular damage
     - Synthesis of protective proteins
     - Restoration of cellular homeostasis
   - These recovery processes strengthened the oysters' resilience to the second heat stress

### Question 21 [10 marks]
The question asks whether the proteome data (Fig. 3) reflect the transcriptome data (Fig. 2) regarding proteasome response to heat stress in two Mytilus species.

**Comparison of transcriptome and proteome data**:

1. **Overall patterns**:
   - Transcriptome data (Fig. 2) shows increased expression of proteasome genes with increasing temperature in both species, with M. galloprovincialis showing higher expression than M. trossulus at all temperatures.
   - Proteome data (Fig. 3) shows more complex patterns with different proteasome isoforms showing distinct responses to temperature in each species.

2. **Species differences**:
   - Transcriptome: M. galloprovincialis consistently shows higher proteasome gene expression than M. trossulus across all temperatures.
   - Proteome: The pattern varies by isoform:
     - Some isoforms (α5, α6) show higher levels in M. galloprovincialis
     - Others (α3, α7, β3) show higher levels in M. trossulus
     - Some (β4) show similar levels in both species

3. **Temperature response**:
   - Transcriptome: Both species show increased proteasome gene expression with increasing temperature, with the highest expression at 32°C.
   - Proteome: The response varies by isoform:
     - Some isoforms increase with temperature (α3, α7 in M. trossulus)
     - Others decrease (α5, α6 in M. galloprovincialis)
     - Some show non-linear responses (β3, β4 in both species)

4. **Conclusion**:
   - The proteome data do not directly reflect the transcriptome data
   - This discrepancy can be explained by:
     - Post-transcriptional regulation (mRNA stability, translation efficiency)
     - Post-translational modifications affecting protein stability
     - Protein turnover rates differing from mRNA turnover
     - Different isoforms being regulated differently despite similar gene expression
     - The proteome reflecting the accumulated history of protein synthesis and degradation, while the transcriptome shows current gene expression

5. **Biological significance**:
   - The complex proteasome response at the protein level suggests sophisticated regulation of protein degradation during heat stress
   - M. galloprovincialis may maintain higher baseline levels of certain proteasome subunits (α5, α6) that contribute to its greater heat tolerance
   - M. trossulus shows more dramatic increases in some subunits (α3, α7) during heat stress, possibly reflecting a more reactive rather than preparative strategy

The data highlight the importance of studying both transcriptome and proteome to understand cellular responses to stress, as mRNA levels alone do not predict protein abundance or activity.

## Final Exam - [MC3025F] 2018 - Section 3

### Question 10 [12 marks]

#### a) Difference between EndoH and PNGaseF [2 marks]
**Endoglycosidase H (EndoH)**:
- Cleaves high-mannose and some hybrid N-linked glycans
- Cuts between the two N-acetylglucosamine (GlcNAc) residues in the chitobiose core
- Cannot cleave complex N-glycans that have been processed in the Golgi apparatus
- Leaves one GlcNAc residue attached to the asparagine

**Peptide-N-glycosidase F (PNGaseF)**:
- Cleaves all types of N-linked glycans (high-mannose, hybrid, and complex)
- Cuts between the innermost GlcNAc and the asparagine residue
- Removes the entire glycan structure from the protein
- Converts the asparagine to aspartic acid in the process

#### b) N260Q notation [1 mark]
The N260Q notation indicates a site-directed mutagenesis where:
- The asparagine (N) at position 260 in the amino acid sequence
- Has been replaced with glutamine (Q)
- This mutation eliminates a potential N-glycosylation site (Asn-X-Ser/Thr) by replacing the asparagine that would normally be glycosylated

#### c) Glycosylation status of WT and mutant gp160 [3 marks]
**Wild-type (WT) gp160**:
- Shows partial sensitivity to EndoH (partial shift in molecular weight)
- Shows complete sensitivity to PNGaseF (complete shift to lower molecular weight)
- This pattern indicates WT gp160 has both high mannose glycans (EndoH-sensitive) and complex N-glycans (EndoH-resistant, PNGaseF-sensitive)
- The presence of complex glycans indicates that WT gp160 has been processed through the Golgi apparatus

**Mutant N260Q gp160**:
- Shows complete sensitivity to both EndoH and PNGaseF (similar shift in molecular weight with both enzymes)
- This pattern indicates N260Q gp160 has only high mannose glycans (fully EndoH-sensitive)
- The absence of EndoH-resistant glycans suggests that the N260Q mutation prevents proper trafficking through the Golgi apparatus where complex glycan formation occurs

#### d) Molecular weight difference after deglycosylation [2 marks]
EndoH-deglycosylated WT gp160 has a slightly higher molecular weight than PNGaseF-deglycosylated WT gp160 because:

1. EndoH cleaves between the two GlcNAc residues in the core of high mannose glycans, leaving one GlcNAc residue still attached to each asparagine
2. PNGaseF cleaves between the innermost GlcNAc and the asparagine, removing the entire glycan structure
3. Additionally, PNGaseF converts each glycosylated asparagine to aspartic acid, which can slightly alter the protein's migration in SDS-PAGE
4. The cumulative effect of these differences across multiple glycosylation sites results in the observed molecular weight difference

#### e) Mode of translocation for gp160 [4 marks]
gp160 utilizes **co-translational translocation** into the ER:

1. **Signal sequence recognition**:
   - gp160 contains an N-terminal signal sequence
   - As translation begins, the signal sequence emerges from the ribosome
   - Signal Recognition Particle (SRP) binds to the signal sequence

2. **Targeting to the ER membrane**:
   - The SRP-ribosome-nascent chain complex docks with the SRP receptor on the ER membrane
   - The ribosome is transferred to the Sec61 translocon complex

3. **Co-translational translocation**:
   - Translation continues as the growing polypeptide is threaded through the Sec61 channel
   - The signal sequence is cleaved by signal peptidase in the ER lumen
   - Transmembrane domains of gp160 are laterally released into the ER membrane through the Sec61 lateral gate

4. **N-linked glycosylation**:
   - As the polypeptide emerges into the ER lumen, oligosaccharyltransferase (OST) adds pre-assembled N-linked glycans to asparagine residues in the consensus sequence Asn-X-Ser/Thr
   - These initial high mannose glycans are later processed in the ER and Golgi

### Question 11 [2 marks]
A protein with both a mitochondrial signal sequence and an internal Nuclear Export Signal (NES) would be targeted to the **mitochondria**.

Reasons:
1. The mitochondrial signal sequence is located at the N-terminus of the protein and is recognized first during translation
2. Mitochondrial targeting is co-translational or occurs immediately after translation, before the protein can fold and expose the internal NES
3. Once the protein is imported into mitochondria, the mitochondrial signal sequence is cleaved by mitochondrial processing peptidase
4. The NES would be irrelevant once the protein is inside the mitochondria, as it would be physically separated from the nuclear export machinery in the cytosol
5. The hierarchical nature of targeting signals typically gives priority to N-terminal signals over internal signals

### Question 12 [5 marks]

#### Mechanism of Ran GTP/GDP cycling [2 marks]
```
                   RanGAP
                      ↓
Cytoplasm:  Ran-GTP → Ran-GDP
               ↑         ↓
               |         |
               |         |
               |         |
Nucleus:    Ran-GDP → Ran-GTP
                      ↑
                   RanGEF
                  (RCC1)
```

- **RanGEF (RCC1)**: Located in the nucleus, bound to chromatin
  - Catalyzes exchange of GDP for GTP on Ran
  - Creates high concentration of Ran-GTP in the nucleus

- **RanGAP**: Located in the cytoplasm, associated with nuclear pore complex
  - Stimulates GTP hydrolysis by Ran
  - Creates high concentration of Ran-GDP in the cytoplasm

- **Nuclear transport factors**:
  - Importins bind cargo in cytoplasm, release it upon binding Ran-GTP in nucleus
  - Exportins bind cargo in nucleus when complexed with Ran-GTP, release it upon GTP hydrolysis in cytoplasm

#### Effect of locking Ran in GDP-bound form [3 marks]
If Ran is locked in its GDP-bound form:

1. **Inhibition of nuclear export**:
   - Nuclear export requires Ran-GTP to form stable exportin-cargo-RanGTP complexes
   - Without Ran-GTP, exportin cannot bind cargo efficiently in the nucleus
   - Proteins and RNAs normally exported from the nucleus would accumulate there

2. **Disruption of nuclear import termination**:
   - Nuclear import cargo release requires Ran-GTP to bind importin-β and cause dissociation of the import complex
   - Without Ran-GTP, imported proteins may not be efficiently released from importins in the nucleus
   - This would eventually deplete free importins and halt further import

3. **Disruption of importin recycling**:
   - Return of importin-α to the cytoplasm requires formation of an exportin (CAS)-importin-α-RanGTP complex
   - Without Ran-GTP, importin-α would accumulate in the nucleus
   - This would further inhibit nuclear import due to lack of available importin-α in the cytoplasm

4. **Overall effect**: Bidirectional nuclear transport would be severely compromised, with more immediate and severe effects on export than import.

### Question 13 [4 marks]
Mutation of cysteine to serine would alter the following aspects of post-translational modification:

1. **Disulfide bond formation**:
   - Cysteine residues form disulfide bonds via their thiol (-SH) groups
   - Serine has a hydroxyl (-OH) group instead of a thiol group and cannot form disulfide bonds
   - This would prevent proper protein folding and stabilization in proteins that rely on disulfide bonds
   - Particularly affects secreted and membrane proteins that typically contain disulfide bonds

2. **Palmitoylation**:
   - Cysteine residues can be palmitoylated (addition of palmitic acid via thioester linkage)
   - Serine cannot undergo palmitoylation
   - This would disrupt membrane association of proteins that rely on palmitoylation
   - Would affect protein localization to lipid rafts and other membrane microdomains

3. **Metal coordination**:
   - Cysteine residues often coordinate metal ions (e.g., zinc in zinc finger proteins)
   - Serine cannot effectively coordinate metals in the same manner
   - This would disrupt the function of metalloproteins and zinc finger transcription factors

4. **Redox sensing and regulation**:
   - Cysteine residues serve as redox sensors through reversible oxidation
   - Serine is not redox-sensitive
   - This would impair redox-regulated protein functions and signaling pathways

5. **Prenylation**:
   - C-terminal cysteine residues in CaaX motifs can be prenylated
   - Serine cannot be prenylated
   - This would prevent membrane association of proteins that rely on prenylation (e.g., small GTPases)

### Question 14 [3 marks]
In Anfinsen's second experiment:

**Result**:
When mercaptoethanol (ME) was removed but 8M urea was maintained, the denatured ribonuclease A formed incorrect disulfide bonds but remained inactive.

**Explanation**:
- Removing ME allowed the formation of disulfide bonds between cysteine residues
- However, the presence of 8M urea kept the protein in an unfolded state
- In this unfolded state, disulfide bonds formed randomly between cysteine residues that were close in the linear sequence or happened to be near each other
- These random disulfide bonds did not correspond to the native disulfide pattern
- The resulting protein had incorrect tertiary structure and lacked enzymatic activity

**Conclusion regarding tertiary structure**:
Anfinsen concluded that the amino acid sequence alone determines the tertiary structure of proteins. Specifically:
1. The native tertiary structure represents the thermodynamically most stable conformation
2. This stable conformation is determined solely by the amino acid sequence
3. Proper folding requires appropriate conditions (absence of denaturants) to allow the protein to explore conformational space and find its lowest energy state
4. The correct disulfide bonds form only when the protein is allowed to fold properly, as they stabilize the native structure rather than direct the folding process

### Question 15 [5 marks]

#### a) Basic principles of sandwich ELISA [2 marks]
1. **Capture antibody immobilization**: A specific antibody (capture antibody) is coated onto a solid surface (typically a microplate well)

2. **Sample addition**: The sample containing the target antigen is added and binds to the capture antibody

3. **Detection antibody binding**: A second antibody (detection antibody) that recognizes a different epitope on the antigen is added, forming a "sandwich" with the antigen between two antibodies

4. **Signal generation**: The detection antibody is linked to an enzyme (typically horseradish peroxidase or alkaline phosphatase)

5. **Substrate addition**: A substrate for the enzyme is added, which is converted to a colored product

6. **Quantification**: The intensity of the color is measured spectrophotometrically and is proportional to the amount of antigen in the sample

#### b) Effect of BACE 2 mutation on amyloidogenic pathway [3 marks]
Based on the ELISA results:

1. **C83 levels** (non-amyloidogenic pathway):
   - Vector only: ~50 ng/ml
   - WT BACE: ~120 ng/ml
   - BACE 2 (mutant): ~320 ng/ml
   - The mutant BACE 2 shows dramatically increased C83 production (2.7-fold higher than WT BACE)

2. **C99 levels** (amyloidogenic pathway):
   - Vector only: ~20 ng/ml
   - WT BACE: ~250 ng/ml
   - BACE 2 (mutant): ~80 ng/ml
   - The mutant BACE 2 shows significantly decreased C99 production (3.1-fold lower than WT BACE)

3. **Conclusion**:
   - The BACE 2 mutation does NOT enhance the amyloidogenic pathway; it actually inhibits it
   - The mutation shifts APP processing strongly toward the non-amyloidogenic pathway (increased C83)
   - This suggests that BACE 2 has reduced β-secretase activity but may enhance α-secretase activity
   - Such a mutation would be expected to reduce Aβ production and could potentially be protective against Alzheimer's disease

### Question 16 [4 marks]

#### a) Effect of pH 7.2 on lysosomal function [2 marks]
A pH of 7.2 in lysosomes would have a **negative effect** on their function:

1. **Enzyme inactivation**: Lysosomal hydrolases (proteases, lipases, nucleases, glycosidases) have acidic pH optima (typically pH 4.5-5.0) and are largely inactive at neutral pH

2. **Substrate processing**: Many substrates require protonation for proper enzyme recognition and processing, which doesn't occur efficiently at neutral pH

3. **Protein denaturation**: The acidic environment helps maintain the proper conformation of lysosomal enzymes; neutral pH may alter their structure

4. **Receptor recycling disruption**: The pH gradient is essential for dissociation of ligands from receptors in the endocytic pathway

5. **Autophagy impairment**: Autophagic degradation requires acidic pH; neutralization would impair this critical cellular process

Normal lysosomal pH is maintained at ~4.5-5.0 by the H⁺ ATPase pump, which actively transports protons into the lysosomal lumen. The mutation in this pump has resulted in failure to acidify the lysosome, severely compromising its degradative functions.

#### b) Endocytic vesicle and pathways [2 marks]
The endocytic vesicle that receives internalized material en route to lysosomes is the **early endosome**.

Three pathways by which the early endosome receives material:

1. **Clathrin-mediated endocytosis**: Receptor-ligand complexes are internalized via clathrin-coated vesicles

2. **Caveolae-mediated endocytosis**: Internalization via flask-shaped invaginations rich in caveolin protein

3. **Macropinocytosis**: Non-selective uptake of extracellular fluid and solutes via large endocytic vesicles

Other pathways include:
- Phagocytosis (cell eating): Uptake of large particles like bacteria
- Clathrin-independent endocytosis: Various mechanisms not dependent on clathrin or caveolin
- Autophagy: Delivery of cytoplasmic components via autophagosomes

### Question 17 [5 marks]

#### a) Basic principles of reversible phosphorylation [2 marks]
Reversible phosphorylation is a post-translational modification that regulates protein function through:

1. **Phosphorylation reaction**:
   - Catalyzed by protein kinases
   - Transfers a phosphate group from ATP to specific amino acid residues (typically serine, threonine, or tyrosine)
   - Results in the addition of a negatively charged phosphate group to the protein
   - Requires ATP as a phosphate donor

2. **Dephosphorylation reaction**:
   - Catalyzed by protein phosphatases
   - Removes the phosphate group from the protein
   - Releases inorganic phosphate
   - Restores the protein to its unphosphorylated state

3. **Functional consequences**:
   - Induces conformational changes in the protein
   - Alters protein activity (activation or inhibition)
   - Modifies protein-protein interactions
   - Changes subcellular localization
   - Affects protein stability and turnover

#### b) Necessity of reprobing with antibody 2166 [1 mark]
Reprobing the immunoblot with antibody 2166 (which detects total Htt) was necessary to:

1. Confirm equal immunoprecipitation efficiency across all samples
2. Verify that differences in phospho-T3 signal are due to differences in phosphorylation status rather than differences in total Htt protein levels
3. Ensure that the 3A mutant protein is expressed and immunoprecipitated at levels comparable to the Q7 and Q111 proteins
4. Provide a denominator for calculating the relative phosphorylation level (pT3/total Htt ratio)

#### c) Conclusion from the data [2 marks]
From the data shown:

1. **Phosphorylation status**:
   - Q7 (normal Htt with 7 glutamines) shows robust phosphorylation at T3
   - Q111 (expanded Htt with 111 glutamines) shows significantly reduced phosphorylation at T3
   - 3A mutant (T3 replaced by alanine) shows no phosphorylation, confirming antibody specificity

2. **Protein expression**:
   - All three constructs (Q7, Q111, and 3A) are expressed and immunoprecipitated at similar levels (2166 blot)
   - GAPDH levels confirm equal loading of cell lysates

3. **Conclusion**:
   - Polyglutamine expansion in Huntingtin (Q111) dramatically reduces phosphorylation at threonine 3
   - This suggests that the disease-causing mutation (expanded polyQ) alters the protein's ability to be phosphorylated at this site
   - Since phosphorylation often regulates protein function, this reduction in T3 phosphorylation may contribute to the pathogenic mechanism of Huntington's disease
   - The specific kinase responsible for T3 phosphorylation may have reduced access to this site in the expanded polyQ form

## Final Exam - [MC3025F] 2018 - Section 4

### Question 18 [8 marks]
To determine the mass of the three proteins from the LC-MS spectrum:

The spectrum shows multiple peaks for each protein, representing different charge states of the same molecule. For a protein with mass M and charge z, the m/z value is:

m/z = (M + z)/z

Rearranging: M = z(m/z) - z

For each protein, we need to identify at least two peaks corresponding to different charge states:

**Protein A**:
- Peak at m/z = 681: If this represents charge state z = 15, then M = 15(681) - 15 = 10,200
- Peak at m/z = 637: If this represents charge state z = 16, then M = 16(637) - 16 = 10,176
- Peak at m/z = 598: If this represents charge state z = 17, then M = 17(598) - 17 = 10,149
- Peak at m/z = 563: If this represents charge state z = 18, then M = 18(563) - 18 = 10,116

The average mass of Protein A is approximately 10,160 Da.

**Protein B**:
- Peak at m/z = 501: If this represents charge state z = 12, then M = 12(501) - 12 = 6,000
- Peak at m/z = 459: If this represents charge state z = 13, then M = 13(459) - 13 = 5,954
- Peak at m/z = 426: If this represents charge state z = 14, then M = 14(426) - 14 = 5,950
- Peak at m/z = 397: If this represents charge state z = 15, then M = 15(397) - 15 = 5,940

The average mass of Protein B is approximately 5,960 Da.

**Protein C**:
- Peak at m/z = 334: If this represents charge state z = 9, then M = 9(334) - 9 = 2,997
- Peak at m/z = 301: If this represents charge state z = 10, then M = 10(301) - 10 = 3,000
- Peak at m/z = 274: If this represents charge state z = 11, then M = 11(274) - 11 = 3,003
- Peak at m/z = 251: If this represents charge state z = 12, then M = 12(251) - 12 = 3,000

The average mass of Protein C is approximately 3,000 Da.

To verify these assignments, we can check that adjacent peaks for the same protein differ by a charge of 1:

For adjacent peaks with m/z values m₁ and m₂ corresponding to charges z and z+1:
m₁ = (M + z)/z
m₂ = (M + z+1)/(z+1)

This relationship confirms our mass assignments for all three proteins.

### Question 19 [10 marks]
**1-D Gel Electrophoresis**:

Principles:
- Separates proteins based on a single physical property, typically molecular weight
- Uses sodium dodecyl sulfate (SDS) to denature proteins and give them a uniform negative charge
- Proteins migrate through a polyacrylamide gel matrix when an electric field is applied
- Smaller proteins migrate faster than larger ones
- The distance traveled is inversely proportional to the logarithm of the molecular weight

Limitations:
- Limited resolution - can only separate a few dozen proteins
- Cannot distinguish proteins of similar molecular weight
- Provides no information about protein isoforms or post-translational modifications
- Not suitable for complex protein mixtures

**2-D Gel Electrophoresis**:

Principles:
- Separates proteins based on two independent physical properties in two sequential steps
- First dimension: Isoelectric focusing (IEF)
  - Separates proteins based on their isoelectric point (pI)
  - Uses a pH gradient in which proteins migrate until they reach the pH where their net charge is zero
  - At this point (pI), proteins stop migrating
- Second dimension: SDS-PAGE
  - Separates proteins based on molecular weight
  - Perpendicular to the first dimension
  - Uses SDS to give proteins uniform negative charge

Advantages of 2-D over 1-D for proteomics:
1. **Superior resolution**: Can resolve thousands of proteins in a single gel
2. **Separation of protein isoforms**: Can distinguish proteins with the same molecular weight but different pI
3. **Detection of post-translational modifications**: Modifications often change pI, creating horizontal "trains" of spots
4. **Visualization of protein expression patterns**: Allows comparison of protein abundance across different conditions
5. **Identification of protein complexes**: Components of complexes can be identified by their migration pattern
6. **Quantitative analysis**: Spot intensity correlates with protein abundance
7. **Compatibility with downstream analysis**: Spots can be excised for mass spectrometry identification

Limitations of 2-D gel electrophoresis:
1. Labor-intensive and time-consuming
2. Poor reproducibility between gels
3. Limited dynamic range
4. Bias against hydrophobic, very basic, very acidic, very large, or very small proteins
5. Difficulty detecting low-abundance proteins

Despite these limitations, 2-D gel electrophoresis remains valuable in proteomics for its ability to provide a visual "map" of the proteome and resolve protein isoforms that might be missed by other techniques.

### Question 20 [12 marks]

#### a) Peptide abundance in heat-stressed abalone [4 marks]
The mass spectrum shows the relative abundance of the peptide in control (labels 114 and 115) versus heat-stressed (labels 116 and 117) abalone haemocytes.

The iTRAQ reporter ion intensities indicate:
- Labels 114 and 115 (controls at 15°C): Low intensity peaks
- Labels 116 and 117 (heat stress at 23°C): High intensity peaks, approximately 3-4 times higher than controls

This indicates that the peptide is significantly upregulated in abalone haemocytes in response to elevated temperature (23°C). The consistent increase across both biological replicates (116 and 117) compared to both control replicates (114 and 115) suggests this is a reliable heat-stress response rather than biological variation.

The upregulation of this peptide suggests it belongs to a protein involved in the cellular response to heat stress, possibly a heat shock protein, chaperone, or another stress-response protein that helps maintain cellular homeostasis under elevated temperature conditions.

#### b) Possible amino acid sequence(s) [5 marks]
To determine the amino acid sequence, we need to analyze the y-ion series in the mass spectrum. Y-ions represent fragments from the C-terminus of the peptide.

Starting from the lowest mass y-ion and moving up:
- 377.24 Da: This likely represents the first 3 amino acids from the C-terminus
- 434.26 Da: Adding 57.02 Da (Glycine, G)
- 581.32 Da: Adding 147.06 Da (Phenylalanine, F)
- 709.39 Da: Adding 128.07 Da (Glutamine, Q or Lysine, K)
- 822.47 Da: Adding 113.08 Da (Leucine, L or Isoleucine, I)
- 919.52 Da: Adding 97.05 Da (Proline, P)
- 1056.57 Da: Adding 137.05 Da (Histidine, H)
- 1153.62 Da: Adding 97.05 Da (Proline, P)
- 1309.72 Da: Adding 156.10 Da (Arginine, R)
- 1410.77 Da: Adding 101.05 Da (Threonine, T)

Working backwards from the C-terminus, the sequence would be:
T-R-P-H-P-L/I-K/Q-F-G-[???]

For the first three amino acids (377.24 Da), possible combinations include:
- K-P-T (128.09 + 97.05 + 101.05 = 326.19 + 51 (water + H) = 377.19)
- L-P-T (113.08 + 97.05 + 101.05 = 311.18 + 51 = 362.18)

The most likely sequence is:
T-R-P-H-P-L-K-F-G-K-P-T or T-R-P-H-P-I-Q-F-G-K-P-T

Since trypsin cleaves after K and R, and we're assuming two missed cleavages, the K and R positions in this sequence are consistent with tryptic digestion.

#### c) Underlying principle of iTRAQ [3 marks]
The iTRAQ (Isobaric Tag for Relative and Absolute Quantification) method is based on the following principles:

1. **Chemical labeling**: Peptides from different samples are labeled with chemically identical tags that differ only in their isotopic composition

2. **Isobaric nature**: Each tag has the same total mass but different distribution of isotopes between the reporter and balance groups
   - Reporter ions: Small fragment released during MS/MS fragmentation
   - Balance group: Maintains equal total mass across all tags

3. **Multiplexing**: Multiple samples can be analyzed in a single MS run (4-plex or 8-plex)
   - Each sample is labeled with a different tag
   - Labeled samples are mixed and analyzed together

4. **Quantification in MS/MS**: 
   - In MS1, peptides from different samples appear as a single peak (same mass)
   - During fragmentation (MS/MS), reporter ions are released
   - Reporter ions appear at different m/z values (114, 115, 116, 117 in 4-plex)
   - Intensity of each reporter ion reflects the relative abundance of the peptide in each sample

5. **Advantages**:
   - Reduces run-to-run variation
   - Increases throughput
   - Improves quantification accuracy
   - Allows comparison of multiple conditions simultaneously
