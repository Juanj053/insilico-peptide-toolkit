# Research Collaboration & Grant Alignment Framework

**Toolkit:** `insilico-peptide-toolkit` (@Juanj053)  
**Lead Investigator & Software Architect:** Juan José Álvarez (`juanjo.alvarez@elevate-biotech.com`)  
**Institutional Hub:** Elevate Biotech — Computational Bioengineering & Translational Peptide Discovery Unit  
**Repository:** [github.com/Juanj053/insilico-peptide-toolkit](https://github.com/Juanj053/insilico-peptide-toolkit)

---

## 1. Executive Summary

`insilico-peptide-toolkit` is an open-source, mathematically rigorous computational platform built to democratize pre-clinical peptide modeling, mutational scanning, and biophysical characterization. This document defines the institutional strategy for international research grant applications, consortia building, experimental wet-lab validation protocols, and open-science governance.

We actively invite academic research groups, biopharmaceutical laboratories, and clinical research institutions to partner with us in translational research initiatives.

---

## 2. Strategic Alignment with International Research Grants

### 2.1. CYTED (Programa Iberoamericano de Ciencia y Tecnología para el Desarrollo)
- **Thematic Area:** *Área 2: Salud y Biomedicina*
- **Action Line:** *Desarrollo de Nuevos Fármacos y Terapias Biotecnológicas para Enfermedades Crónicas no Transmisibles*
- **Consortium Project Title:** *Red Iberoamericana de Diseño In-Silico y Validación Biofísica de Péptidos Terapéuticos y Reguladores Celulares (RED-INPEPT)*
- **Target Regional Partners:** Academic research centers in Ecuador, Colombia, Mexico, Argentina, Chile, Brazil, and Spain.
- **Role of Toolkit:** Provides the common computational engine for virtual screening, alanine scanning, and biophysical benchmarking of regional biodiversity-derived bioactive peptides.

### 2.2. NIH / Fogarty International Center
- **Funding Opportunity:** *Global Non-Communicable Diseases (NCD) and Mobile Health / Computational Biology*
- **Objective:** Developing cost-effective, open-source computational pipelines that substantially reduce pre-clinical drug discovery costs in low- and middle-income countries (LMICs).
- **Translational Focus:** Rational design of long-acting incretin mimetics (GLP-1/GIP receptor co-agonists) and regenerative peptides (gastric cytoprotective analogues) with optimized stability against gastrointestinal and serum proteases.

### 2.3. Chevening Research Fellowships & UK Academic Alliances
- **Academic Exchange Focus:** Structural biology and peptide engineering partnerships with leading UK universities:
  - **University of Oxford (Department of Pharmacology):** Receptor kinetic modeling and peptide-GPCR cleft interaction analysis.
  - **Imperial College London (Institute of Chemical Biology):** Amphiphilic peptide partition thermodynamics and membrane interaction modeling.
  - **University of Dundee (MRC Protein Phosphorylation and Ubiquitylation Unit):** Proteolytic degradation kinetics and ubiquitination signaling pathways.

---

## 3. Wet-Lab Experimental Validation Protocols

To bridge computational predictions and empirical biological evidence, partner laboratories follow standard operating procedures (SOPs):

```
┌────────────────────────────────────────────────────────┐
│ 1. In-Silico Screening (insilico-peptide-toolkit)       │
│    - Compute MW, pI, GRAVY, mu_H, II, AI, Cleavage     │
│    - Alanine Scan & Pareto Multi-Objective Ranking     │
└───────────────────────────┬────────────────────────────┘
                            │ Top Candidates
                            ▼
┌────────────────────────────────────────────────────────┐
│ 2. Solid-Phase Synthesis & Mass Confirmation           │
│    - Fmoc/tBu chemistry (SPPS)                         │
│    - RP-HPLC purity verification (>= 98.0%)            │
│    - ESI-MS / MALDI-TOF exact mass confirmation        │
└───────────────────────────┬────────────────────────────┘
                            │
        ┌───────────────────┴───────────────────┐
        ▼                                       ▼
┌───────────────────────────────┐ ┌───────────────────────────────┐
│ 3. Conformational Characterization│ 4. Proteolytic Stability Assay│
│    - Far-UV CD Spectroscopy   │    - Human plasma incubation  │
│      (190–260 nm, buffer/TFE) │      (37°C; 0, 1, 4, 8, 24 h) │
│    - Helical content vs mu_H  │    - LC-MS/MS half-life (t1/2)│
└───────────────┬───────────────┘ └───────────────┬───────────────┘
                │                                 │
                └───────────────┬─────────────────┘
                                ▼
┌────────────────────────────────────────────────────────┐
│ 5. Receptor Binding Kinetics (SPR / Biacore)          │
│    - Target immobilization (CM5 sensor chip)          │
│    - Kinetic constants: k_on, k_off, K_D              │
└────────────────────────────────────────────────────────┘
```

### Protocol 1: Solid-Phase Peptide Synthesis (SPPS) & RP-HPLC Purity
1. **Synthesis:** Automated microwave-assisted solid-phase synthesis using orthogonal Fmoc/tBu chemistry on Rink Amide AM resin (for amidated C-termini) or Wang resin (for free acids).
2. **Deprotection & Cleavage:** Cleavage cocktail: 92.5% Trifluoroacetic acid (TFA), 2.5% Triisopropylsilane (TIS), 2.5% $\text{H}_2\text{O}$, 2.5% 1,2-Ethanedithiol (EDT) for 3 hours at room temperature.
3. **Purification:** Preparative reversed-phase HPLC on a C18 column ($250 \times 21.2\text{ mm}$, $5\text{ }\mu\text{m}$, 300 Å) with a linear gradient of $0.1\%$ TFA in acetonitrile against $0.1\%$ aqueous TFA.
4. **Purity Standard:** Analytical RP-HPLC must confirm $\ge 98.0\%$ chromatographic purity at 214 nm and 280 nm.

### Protocol 2: High-Resolution Mass Spectrometry (HRMS)
- **Methodology:** Electrospray Ionization Quadrupole Time-of-Flight (ESI-Q-TOF) or MALDI-TOF MS.
- **Verification Criterion:** Measured monoisotopic and average $[M+H]^+$ ions must match `insilico-peptide-toolkit` calculated values within $\pm 0.05\text{ Da}$ (monoisotopic) and $\pm 0.15\text{ Da}$ (average).

### Protocol 3: Circular Dichroism (CD) Secondary Structure Determination
1. **Instrument:** Jasco J-815 or Chirascan spectrometer with $1\text{ mm}$ quartz cuvette.
2. **Buffer Conditions:** $20\text{ }\mu\text{M}$ peptide in $10\text{ mM}$ potassium phosphate buffer ($\text{pH } 7.4$), with and without $30\%\text{ (v/v)}$ 2,2,2-Trifluoroethanol (TFE).
3. **Spectral Range:** $190\text{–}260\text{ nm}$ at $25^\circ\text{C}$ (step $0.5\text{ nm}$, bandwidth $1.0\text{ nm}$).
4. **Analysis:** Deconvolution of mean residue ellipticity $[\theta]$ at $208\text{ nm}$ and $222\text{ nm}$ to determine experimental percentage $\alpha$-helix content, validating the predicted amphiphilic hydrophobic moment ($\mu_H$).

### Protocol 4: Surface Plasmon Resonance (SPR / Biacore) Binding Affinity
1. **Sensor:** Biacore T200 / 8K using a CM5 carboxymethylated dextran chip.
2. **Immobilization:** Recombinant target extracellular domain (e.g. GLP-1R or Insulin-R) immobilized via standard EDC/NHS amine coupling ($~1000\text{ RU}$).
3. **Kinetic Titration:** Multi-cycle or single-cycle kinetics at $25^\circ\text{C}$ across concentrations ($0.1\text{ nM}$ to $1000\text{ nM}$) in HBS-EP+ buffer.
4. **Kinetic Evaluation:** Global fitting to a 1:1 Langmuir binding model to determine association rate constant ($k_{\text{on}}$), dissociation rate constant ($k_{\text{off}}$), and equilibrium dissociation constant ($K_D = k_{\text{off}} / k_{\text{on}}$).

### Protocol 5: In-Vitro Human Plasma Stability Assay
1. **Incubation:** Synthetic peptide ($10\text{ }\mu\text{M}$) spiked into fresh pooled human plasma ($37^\circ\text{C}$).
2. **Sampling:** Aliquots removed at $t = 0, 15\text{ min}, 30\text{ min}, 1\text{ h}, 2\text{ h}, 4\text{ h}, 8\text{ h}, 24\text{ h}$.
3. **Quenching:** Plasma proteins precipitated by addition of ice-cold acetonitrile with internal standard, followed by centrifugation ($14,000 \times g$, $15\text{ min}$, $4^\circ\text{C}$).
4. **Quantification:** Targeted LC-MS/MS multiple reaction monitoring (MRM) to determine residual intact peptide concentration.
5. **Correlation:** Experimental in-vitro half-life ($t_{1/2}$) benchmarked against predicted Guruprasad Instability Index ($II$) and enzymatic vulnerability scores.

---

## 4. Open Science, Intellectual Property & Data Governance

### 4.1. FAIR Principles Compliance
- **Findable:** All release versions assigned a persistent citable Digital Object Identifier (DOI) via Zenodo.
- **Accessible:** Free, unrestricted open-source distribution under the OSI-approved MIT License.
- **Interoperable:** Standard data exchange formats (FASTA, JSON, Markdown, mmCIF, CSV) ensuring integration with existing bioinformatics workflows.
- **Reusable:** Fully documented code, explicit type annotations, reproducible unit test fixtures, and comprehensive API documentation.

### 4.2. Intellectual Property on Novel Designed Analogues
- **Platform vs. Molecules:** The `insilico-peptide-toolkit` codebase is permanently open-source under the MIT license.
- **Novel Chemical Matter:** Any de novo peptide sequence or chemical analogue engineered using the toolkit by collaborative partners remains the exclusive intellectual property of the discovering research team/institution.
- **Co-Development Agreements:** Formal Material Transfer Agreements (MTA) and Joint Intellectual Property Agreements (JIPA) are executed prior to wet-lab synthesis of patentable therapeutic leads.

### 4.3. Authorship & Academic Recognition
All collaborative publications resulting from research utilizing the toolkit adhere strictly to the International Committee of Medical Journal Editors (ICMJE) Vancouver recommendations for authorship and acknowledgment.

---

## 5. Contact & Collaboration Inquiries

To propose a joint grant submission, request access to customized screening pipelines, or initiate pre-clinical validation studies:

- **Juan José Álvarez**  
  Chief Bioengineer & Lead Architect  
  Elevate Biotech & Agromanper Bioengineering  
  Email: `juanjo.alvarez@elevate-biotech.com`  
  GitHub: [@Juanj053](https://github.com/Juanj053)  
  Ambato / Quito, Ecuador
