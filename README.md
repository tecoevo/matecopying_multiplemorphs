# Phenotypic Polymorphism via Mate Copying

**Authors:**  
Srishti Patil¹˒², Sabine Nöbel³, and Chaitanya S. Gokhale⁴*  
¹ Division of Theoretical Systems Biology, German Cancer Research Center (DKFZ), Heidelberg, Germany  
² Indian Institute of Science Education and Research (IISER), Pune, India  
³ Department of Zoology, Animal Ecology, Martin-Luther-University Halle-Wittenberg, Halle (Saale), Germany  
⁴ Center for Computational and Theoretical Biology, Julius-Maximilians University, Würzburg, Germany  

*Corresponding author: [chaitanya.gokhale@uni-wuerzburg.de]*

---

## Overview

This repository contains all simulation data, analysis scripts, and figure-generation files associated with the paper:

> **Patil, S., Nöbel, S., & Gokhale, C. S. (2025).  
> Phenotypic polymorphism via mate copying.**  
>  
> **Abstract:**  
> Classical mate choice theories assume independent decision-making, yet mounting evidence shows that individuals often use social information and copy conspecifics' mate choices — a behaviour termed *mate copying*. While this has been documented across vertebrates and notably in *Drosophila melanogaster*, theoretical and experimental work has largely been restricted to binary choice scenarios.  
>  
> Here, we develop a theoretical model of mate copying in populations with multiple morphs, incorporating both private (inherent) and public (cultural) information in mate choice decisions. The population dynamics of the male morphs depend on varying intensities of conformist and anticonformist mate copying. We demonstrate that mate copying can lead to fixation of low-quality morphs and identify conditions that maintain polymorphism among morphs.  
>  
> This framework reveals how social learning in mate choice can shape evolutionary trajectories and sustain phenotypic diversity, with implications for sexual selection and speciation.

---

## Repository Structure

<pre>
├── LICENSE
├── New_Revised_Figures
│   ├── Fig_1_cartoon
│   ├── Fig_2_confanticonf
│   │   ├── Fig_2_new.nb
│   │   ├── Fig_2confanticonf.graffle
│   │   ├── Fig_2confanticonf.pdf
│   │   ├── type1conf.csv
│   │   ├── typeII_anticonformist_f1.2
│   │   │   ├── image.pdf
│   │   │   ├── out.csv
│   │   │   └── params.csv
│   │   ├── typeII_conformist_f1.2
│   │   │   ├── image.pdf
│   │   │   ├── out.csv
│   │   │   └── params.csv
│   │   ├── typeII_mixed_f1.2
│   │   │   ├── image.pdf
│   │   │   ├── out.csv
│   │   │   └── params.csv
│   │   └── typeI_anticonformist_b-2
│   │       ├── image.pdf
│   │       ├── out.csv
│   │       └── params.csv
│   ├── Fig_3_overlays
│   │   ├── Fig3_overlays.graffle
│   │   ├── Fig3_overlays.pdf
│   │   ├── Fig_3.nb
│   │   ├── Fig_3.pdf
│   │   ├── copying_functions.py
│   │   ├── fixedpoint3morphs.nb
│   │   ├── typeII_conformist_f1.2
│   │   │   ├── out_c0.0.csv
│   │   │   ├── out_c0.05.csv
│   │   │   ├── out_c0.1.csv
│   │   │   ├── out_c0.15.csv
│   │   │   ├── out_c0.2.csv
│   │   │   ├── out_c0.25.csv
│   │   │   ├── out_c0.3.csv
│   │   │   ├── out_c0.35.csv
│   │   │   ├── out_c0.4.csv
│   │   │   ├── out_c0.45.csv
│   │   │   ├── out_c0.5.csv
│   │   │   ├── out_c0.55.csv
│   │   │   ├── out_c0.6.csv
│   │   │   ├── out_c0.65.csv
│   │   │   ├── out_c0.7.csv
│   │   │   ├── out_c0.75.csv
│   │   │   ├── out_c0.8.csv
│   │   │   ├── out_c0.85.csv
│   │   │   ├── out_c0.9.csv
│   │   │   ├── out_c0.95.csv
│   │   │   ├── out_c1.0.csv
│   │   │   └── params.csv
│   │   └── typeI_conformist_b2
│   ├── Fig_4_dynamicGamma
│   │   ├── Fig4_dynamicgamma.pdf
│   │   ├── Fig_4dynamicgamma.pdf
│   │   ├── Fig_dynamicgamma.graffle
│   │   ├── fig_dynamicGamma.nb
│   │   └── fig_dynamicGamma2morphs.nb
│   ├── Fig_SI_anticonf
│   │   ├── Fig_SI_anticonf.nb
│   │   ├── typeII_anticonformist_f1.2
│   │   └── typeI_anticonformist_b-2
│   └── Fig_SI_fixtime
│       ├── Fig_2.nb
│       ├── Fig_fixationtime_23morphs.pdf
│       ├── Untitled.pdf
│       ├── fig2a_data.csv
│       ├── fig2b_data.csv
│       ├── fig2c_data.csv
│       ├── fixationtime3morphs.pdf
│       ├── fixationtime_23morphs.graffle
│       ├── threemorphdynamics_fixtime.mp4
│       └── threemorphdynamics_fixtime.nb
├── OlderFolders
│   ├── Fig_dynamicgamma
│   │   ├── Fig_dynamicgamma.graffle
│   │   ├── Fig_dynamicgamma.pdf
│   │   └── fig_dynamicGamma.nb
│   ├── SI_copyingfunction
│   │   ├── copyingfunction.nb
│   │   └── copyingfunction.pdf
│   ├── SI_fig_fixtime
│   │   ├── Untitled.pdf
│   │   ├── fixationtime3morphs.graffle
│   │   ├── fixationtime3morphs.pdf
│   │   ├── threemorphdynamics_fixtime.mp4
│   │   └── threemorphdynamics_fixtime.nb
│   ├── SI_threemorphs_anticonformism
│   │   ├── Fig_threemorph_SI.nb
│   │   ├── typeII_anticonformist_f1.2
│   │   ├── typeII_mixed_f1.2
│   │   └── typeI_anticonformist_b-2
│   ├── fig2
│   │   └── scripts
│   │       └── simulation-data
│   ├── fig2_new
│   │   ├── Fig_2.graffle
│   │   ├── Fig_2.nb
│   │   ├── Fig_2.pdf
│   │   ├── Fig_2_new.nb
│   │   ├── Fig_2_new.pdf
│   │   ├── fig2a_data.csv
│   │   ├── fig2b_data.csv
│   │   └── fig2c_data.csv
│   ├── fig3_new
│   │   ├── Fig3_new.graffle
│   │   ├── Fig3_new.pdf
│   │   ├── Fig_3.nb
│   │   ├── cdynamics.mov
│   │   ├── plotting_2morphs.ipynb
│   │   ├── simulation_data
│   │   ├── typeII_conformist_f1.2
│   │   └── typeI_conformist_b2
│   ├── fig4_new
│   │   ├── Fig4_numerical.graffle
│   │   ├── Fig4_numerical.pdf
│   │   └── threemorphdynamics.nb
│   ├── plotting_2morphs.ipynb
│   └── supplementary_figuredata
│       ├── 2m_r1.5.csv
│       ├── 2m_r1.5_6_4.csv
│       ├── 2m_r1.8.csv
│       ├── 2m_r1.csv
│       ├── 2m_r4.2.csv
│       ├── Fig_SI.graffle
│       ├── Fig_SI.nb
│       └── Fig_SI.pdf
├── README.md
└── simulations
    ├── data
    │   ├── 2morphs
    │   │   ├── typeII_anticonformist_f1.2
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   ├── typeII_anticonformist_f2
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   ├── typeII_anticonformist_f2_longer
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   ├── typeII_conformist_f1.2
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   ├── typeII_conformist_f2
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   ├── typeII_mixed_f1.2
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   ├── typeII_mixed_f1.2_lowthreshold
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   ├── typeII_mixed_f2_longer
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   ├── typeII_mixed_f2_lowthreshold
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   ├── typeI_anticonformist_b-0.1
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   ├── typeI_anticonformist_b-1
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   ├── typeI_anticonformist_b-2
    │   │   │   ├── image.pdf
    │   │   │   ├── out.csv
    │   │   │   └── params.csv
    │   │   └── typeI_conformist_r2
    │   │       ├── out.csv
    │   │       └── params.csv
    │   └── 3morphs
    │       ├── typeII_anticonformist_f1.2
    │       ├── typeII_anticonformist_f2
    │       ├── typeII_conformist_f1.2
    │       ├── typeII_conformist_f2
    │       ├── typeII_mixed_f1.2
    │       ├── typeII_mixed_f1.2_lowthreshold
    │       ├── typeII_mixed_f2
    │       ├── typeII_mixed_f2_lowthreshold
    │       ├── typeI_anticonformist_b-2
    │       ├── typeI_conformist_b1
    │       ├── typeI_conformist_b2
    ├── gamma_evolution
    │   ├── functions.py
    │   ├── main.py
    │   └── params.py
    ├── matecopying_functions.py
    ├── matecopying_main.py
    └── matecopying_params.py
</pre>
    
---

### 📁 `simulations/`

Contains all numerical simulation data used to generate the figures and results in the paper.  
Simulations model the evolutionary dynamics of male morph frequencies under different mate-copying regimes (conformist, anticonformist, and mixed).

#### Structure
- **`2morphs/`** — Data for two-morph systems.
  - Subdirectories correspond to specific parameter sets (e.g., `typeI_conformist`, `typeII_anticonformist`).
  - Files:
    - `out.csv` — Frequency trajectories or equilibria of morphs.  
    - `params.csv` — Parameters used for that simulation (fitness values, copying intensity, etc.).  
    - `image.pdf` — Visualization of the simulation outcome.

- **`3morphs/`** — Data for three-morph systems.
  - Subdirectories correspond to different conformity intensities (e.g., `c0.0`, `c0.5`, `c1.0`).  
  - Each includes `.csv` outputs for morph frequency dynamics and equilibrium analysis.

These data form the quantitative backbone for all figures in the manuscript.

---

### 📁 `New_Revised_Figures/`

Contains **all scripts, data, and outputs for the revised manuscript figures** (post peer-review).  
Each folder represents one figure, including the Mathematica notebooks and raw data required for reproduction.

| Folder | Description |
|---------|--------------|
| `Fig_1_cartoon/` | Conceptual schematic illustrating the model of mate copying, private vs. public information, and morph interactions. |
| `Fig_2_confanticonf/` | Analytical and numerical exploration of conformist vs. anticonformist dynamics.<br>Includes `.nb` notebooks, `.csv` data files, and figure PDFs. |
| `Fig_3_overlays/` | Overlays showing equilibrium landscapes and stability regions across copying intensities.<br>Includes `copying_functions.py` and Mathematica notebooks for figure generation. |
| `Fig_4_dynamicGamma/` | Figures showing dynamic variation in γ (gamma) — the parameter controlling the relative strength of copying influence. |
| `Fig_SI_anticonf/` | Supplementary analyses of extended anticonformist conditions and their stability outcomes. |
| `Fig_SI_fixtime/` | Supplementary results showing fixation times of morphs and temporal dynamics.<br>Contains `.csv` data and an animation file `threemorphdynamics_fixtime.mp4`. |

All final figure panels are generated as `.pdf` files directly from the notebooks.

---

### 📁 `OlderFolders/`

Archived scripts, data, and figure-generation files from the **original manuscript submission**.  
Preserved for transparency, showing the progression of analysis from submission to revision.

| Folder | Description |
|---------|--------------|
| `Fig_dynamicgamma/`, `fig2_new/`, `fig3_new/`, `fig4_new/` | Earlier versions of the main figures. |
| `SI_threemorphs_anticonformism/`, `SI_fig_fixtime/` | Original supplementary figure data and visuals. |
| `supplementary_figuredata/` | Legacy data for supplementary information. |
| `plotting_2morphs.ipynb` | Python notebook for preliminary 2-morph plotting. |

---

### 📄 `LICENSE`

Specifies reuse conditions for all code and data in this repository.  
Please refer to this file before redistributing or adapting any content.

---

## Reproducing Figures

### Requirements

- **Mathematica** ≥ 13 — For symbolic and numerical analyses (`.nb` notebooks).  
- **Python 3.x** — For selected visualization scripts (e.g. `copying_functions.py`, requiring `numpy`, `matplotlib`, and `pandas`).  
- **OmniGraffle** *(optional)* — For graphical layout of figure panels (`.graffle` files).  

### Instructions

1. Navigate to the relevant directory under `New_Revised_Figures/`.  
2. Open the corresponding Mathematica notebook (e.g., `Fig_3_overlays.nb`).  
3. Evaluate all cells sequentially.  
4. The notebook will automatically load required `.csv` data from the included folders.  
5. The output `.pdf` figure will match the version shown in the manuscript.

---

## Data Summary

Each `.csv` file contains:
- **Time-series or equilibrium data** for morph frequencies.  
- **Model parameters**, including fitness values, conformity intensity, and γ.  
- **Simulation metadata** (replicate count, timestep length, etc.).  

Where relevant, visual outputs (`.pdf`) and videos (`.mp4`) illustrate the population dynamics.

---

## Versioning and Provenance

- **Initial submission:** Data and figure scripts archived in `OlderFolders/`.  
- **Revised manuscript:** Updated analysis, expanded parameter sweeps, and final figures in `New_Revised_Figures/`.  
- **Core simulations:** Shared datasets under `simulations/`, used by both versions.

This structure ensures full transparency and reproducibility across manuscript stages.

---

## Citation

If you use this code or data, please cite:

> Patil, S., Nöbel, S., & Gokhale, C. S. (2025).  
> *Phenotypic polymorphism via mate copying.*  
> Journal will be updated later.
---

## Contact

For questions, collaborations, or reproducibility inquiries, please contact:  
**Chaitanya S. Gokhale**  
Center for Computational and Theoretical Biology, University of Würzburg  
📧 chaitanya.gokhale@uni-wuerzburg.de  

---

© 2025 Patil, Nöbel, & Gokhale.  
All rights reserved under the terms specified in the LICENSE file.