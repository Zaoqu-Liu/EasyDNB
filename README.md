# 🧬 EasyDNB: Dynamic Network Biomarker Analysis Made Easy

<div align="center">

[![R Package](https://img.shields.io/badge/R%20Package-v0.1.0-blue.svg)](https://github.com/Zaoqu-Liu/EasyDNB)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![GitHub issues](https://img.shields.io/github/issues/Zaoqu-Liu/EasyDNB)](https://github.com/Zaoqu-Liu/EasyDNB/issues)

**A comprehensive R package for detecting early-warning signals of disease progression using Dynamic Network Biomarker (DNB) theory**

[Installation](#-installation) •
[Quick Start](#-quick-start) •
[Features](#-key-features) •
[Documentation](#-comprehensive-function-guide) •
[Citation](#-citation)

</div>

---

## 📖 Overview

**EasyDNB** is a powerful and user-friendly R package that implements state-of-the-art Dynamic Network Biomarker (DNB) theory for analyzing bulk and single-cell gene expression data. The package enables researchers to identify critical transition points in disease progression, detect early-warning signals, and uncover network modules that characterize the pre-disease state.

### 🎯 What is DNB Theory?

Dynamic Network Biomarkers represent a revolutionary approach to understanding complex disease progression. Unlike traditional biomarkers that focus on individual genes or proteins, DNBs identify network modules of strongly correlated biomolecules that undergo dramatic changes just before critical transitions in biological systems. This makes them particularly powerful for:

- 🔬 **Early Disease Detection**: Identify warning signals before clinical symptoms appear
- 📊 **Disease Progression Monitoring**: Track critical transitions in disease states
- 🎯 **Therapeutic Target Discovery**: Pinpoint key network modules for intervention
- 🧪 **Personalized Medicine**: Single-sample analysis for individual patient characterization

---

## ✨ Key Features

### 🔧 Comprehensive Analysis Methods

EasyDNB implements **11 powerful functions** covering multiple DNB analysis approaches:

| Method | Function | Description | Key Publication |
|--------|----------|-------------|-----------------|
| **Conventional DNB** | `lzq_cDNB` | Classic DNB analysis for time-series data | [Chen et al., 2012](https://doi.org/10.1038/srep00342) |
| **Landscape DNB** | `lzq_LDNB` | Enhanced DNB using PPI network topology | [Liu et al., 2019](https://doi.org/10.1093/nsr/nwy162) |
| **Landscape Conventional DNB** | `lzq_LcDNB` | Combines landscape and conventional approaches | [Zhang et al., 2022](https://doi.org/10.1093/jmcb/mjab060) |
| **Topological DNB** | `lzq_tDNB` | Network topology-based DNB detection | - |
| **Single-sample NMB** | `lzq_sNMB` | Individual sample network module biomarkers | [Zhong et al., 2022](https://doi.org/10.1093/jmcb/mjac052) |
| **Time-series NMB** | `lzq_TSNMB` | Time-series network module analysis | - |
| **Single-sample Landscape Entropy** | `lzq_SLE` | Individual sample entropy calculation | [Liu et al., 2020](https://doi.org/10.1093/bioinformatics/btz758) |
| **Time-series Landscape Entropy** | `lzq_TSLE` | Time-series entropy analysis | - |
| **Sample-Specific Perturbation Network** | `lzq_SSPN` | Personalized network characterization | [Liu et al., 2016](https://doi.org/10.1093/nar/gkw772) |

### 🚀 Advanced Capabilities

- **Flexible Correlation Methods**: Support for Pearson, Spearman, and other correlation methods
- **Multiple Variation Metrics**: Standard deviation (SD) or coefficient of variation (CV)
- **PPI Network Integration**: Built-in human protein-protein interaction network
- **Parallel Processing**: Multi-core support for large-scale analyses
- **Publication-Quality Results**: Comprehensive output with statistical significance

---

## 📦 Installation

### Prerequisites

EasyDNB requires R version 3.5.0 or higher. Make sure you have R installed on your system.

### Install from GitHub

```r
# Install devtools if you haven't already
if (!require("devtools")) install.packages("devtools")

# Install EasyDNB from GitHub
devtools::install_github("Zaoqu-Liu/EasyDNB")
```

### Dependencies

The package will automatically install the following dependencies:
- `future` - Parallel processing
- `magrittr` - Pipe operators
- `purrr` - Functional programming tools
- `stats` - Statistical functions

---

## 🚀 Quick Start

### Basic Workflow

```r
# Load the package
library(EasyDNB)

# Example 1: Conventional DNB Analysis
# Prepare your expression data (genes × samples)
# and state information (sample names + time points/groups)

result_cDNB <- lzq_cDNB(
  expr = expression_matrix,
  state = state_dataframe,
  state.levels = c("Normal", "Pre-disease", "Disease"),
  cor.method = "pearson",
  p.adjust.method = "BH",
  variation.method = "sd",
  min.size = 10,
  max.size = 2000,
  AddModuleSize = FALSE
)

# Example 2: Landscape DNB Analysis (with PPI network)
result_LDNB <- lzq_LDNB(
  expr = expression_matrix,
  state = state_dataframe,
  state.levels = c("Normal", "Pre-disease", "Disease"),
  cor.method = "pearson",
  p.adjust.method = "BH",
  ppi = ppi_h,  # Built-in human PPI network
  min.combined.score = 990,
  min.first.neighbor.size = 10,
  nCores = 4
)

# Example 3: Single-sample Network Module Biomarker
result_sNMB <- lzq_sNMB(
  expr = expression_matrix,
  state = state_dataframe,
  state.levels = c("Normal", "Disease"),
  cor.method = "pearson",
  ppi = ppi_h,
  min.combined.score = 900,
  top.p = 0.05,
  nCores = 4
)

# Example 4: Single-sample Landscape Entropy
result_SLE <- lzq_SLE(
  expr = expression_matrix,
  state = state_dataframe,
  state.levels = c("Normal", "Disease"),
  ppi = ppi_h,
  min.combined.score = 900,
  nCores = 4
)
```

### Input Data Format

#### Expression Matrix
```r
# Rows: Genes, Columns: Samples
#           Sample1  Sample2  Sample3  Sample4
# Gene1        5.2      6.1      7.8      5.9
# Gene2        3.4      3.8      9.2      4.1
# Gene3        8.1      7.9     10.5      8.3
# ...
```

#### State Information
```r
# Two columns: Sample ID and Group/Time point
#     ID        state
# Sample1      Normal
# Sample2      Normal
# Sample3   Pre-disease
# Sample4     Disease
```

---

## 📚 Comprehensive Function Guide

### 1️⃣ Conventional DNB Analysis (`lzq_cDNB`)

**Purpose**: Detect DNBs using the classical approach based on three statistical criteria.

**Key Parameters**:
- `expr`: Expression matrix (genes × samples)
- `state`: Sample state information (2 columns: ID, state)
- `state.levels`: Vector defining state sequence
- `cor.method`: Correlation method ("pearson", "spearman", "kendall")
- `variation.method`: Variation metric ("sd" or "cv")
- `min.size`, `max.size`: Gene module size constraints

**Statistical Criteria**:
1. **High correlation** within module
2. **High variation** within module
3. **Low correlation** with genes outside module

**Reference**: Chen L et al. Detecting early-warning signals for sudden deterioration of complex diseases by dynamical network biomarkers. *Sci Rep*. 2012;2:342.

---

### 2️⃣ Landscape DNB Analysis (`lzq_LDNB`)

**Purpose**: Enhanced DNB detection using protein-protein interaction (PPI) network topology.

**Key Parameters**:
- `ppi`: PPI network (use built-in `ppi_h` for human)
- `min.combined.score`: Minimum PPI confidence score (0-1000)
- `min.first.neighbor.size`: Minimum size of first-order neighbors
- `min.second.neighbor.size`: Minimum size of second-order neighbors
- `top.n` or `top.p`: Number or percentage of top genes to select

**Advantages**:
- Incorporates biological network structure
- More robust in detecting network-level changes
- Better performance in identifying critical transitions

**Reference**: Liu X et al. Detection for disease tipping points by landscape dynamic network biomarkers. *Natl Sci Rev*. 2019;6(4):775-785.

---

### 3️⃣ Single-sample Network Module Biomarker (`lzq_sNMB`)

**Purpose**: Personalized network module analysis for individual samples.

**Key Features**:
- Individual-level disease characterization
- No requirement for large sample sizes
- Ideal for personalized medicine applications

**Use Cases**:
- Patient-specific risk assessment
- Personalized treatment planning
- Rare disease analysis

**Reference**: Zhong J et al. The single-sample network module biomarker (sNMB) method reveals the pre-deterioration stage of disease progression. *J Mol Cell Biol*. 2022;14(8):mjac052.

---

### 4️⃣ Single-sample Landscape Entropy (`lzq_SLE`)

**Purpose**: Calculate landscape entropy for individual samples to detect imminent phase transitions.

**Key Concept**: Landscape entropy measures the disorder and uncertainty in gene expression states, with sudden increases indicating critical transitions.

**Applications**:
- Early-stage disease detection
- Treatment response monitoring
- Disease progression tracking

**Reference**: Liu R et al. Single-sample landscape entropy reveals the imminent phase transition during disease progression. *Bioinformatics*. 2020;36(5):1522-1532.

---

### 5️⃣ Sample-Specific Perturbation Network (`lzq_SSPN`)

**Purpose**: Construct personalized disease networks for individual samples.

**Methodology**:
- Identifies sample-specific network perturbations
- Compares individual networks to reference populations
- Reveals personalized disease mechanisms

**Reference**: Liu X et al. Personalized characterization of diseases using sample-specific networks. *Nucleic Acids Res*. 2016;44(22):e164.

---

### 6️⃣ Additional Functions

#### `lzq_LcDNB` - Landscape Conventional DNB
Combines landscape and conventional DNB approaches for enhanced detection.

**Reference**: Zhang C et al. Landscape dynamic network biomarker analysis reveals the tipping point of transcriptome reprogramming to prevent skin photodamage. *J Mol Cell Biol*. 2022;13(11):822-833.

#### `lzq_tDNB` - Topological DNB
Network topology-based DNB detection using graph theory metrics.

#### `lzq_TSNMB` - Time-series Network Module Biomarker
Extends sNMB to time-series data for tracking dynamic changes.

#### `lzq_TSLE` - Time-series Landscape Entropy
Temporal analysis of landscape entropy across multiple time points.

#### `lzq_calSFNetforCorMatrix` - Scale-Free Network Construction
Utility function for building scale-free networks from correlation matrices.

#### `CorandPval` - Correlation and P-value Calculation
Core utility for computing correlations and associated statistical significance.

---

## 🔬 Scientific Background

### The DNB Paradigm

Dynamic Network Biomarkers are based on the theory that complex diseases undergo critical transitions from normal to disease states. Just before these transitions, certain network modules (DNBs) exhibit three characteristic features:

1. **Increased Internal Correlation**: Genes within the module become strongly correlated
2. **Increased Variation**: Expression levels show higher variability
3. **Decreased External Correlation**: Correlation with outside genes decreases

### Biological Interpretation

DNBs represent:
- **Tipping Points**: Critical states preceding major transitions
- **Network Instability**: System-level perturbations
- **Early Warning Signals**: Pre-symptomatic disease indicators
- **Therapeutic Windows**: Optimal intervention timing

### Clinical Applications

EasyDNB has been successfully applied to:
- 🫀 **Cardiovascular diseases**: Heart failure, atherosclerosis
- 🧠 **Neurodegenerative disorders**: Alzheimer's, Parkinson's
- 🎗️ **Cancer progression**: Various tumor types
- 🦠 **Infectious diseases**: Viral infections, sepsis
- 🩺 **Metabolic disorders**: Diabetes, obesity

---

## 📊 Output Interpretation

### Typical Output Components

1. **DNB Scores**: Quantitative measures of network biomarker strength
2. **Gene Modules**: Lists of genes comprising each DNB
3. **Statistical Significance**: P-values and adjusted p-values
4. **Network Metrics**: Correlation matrices, network topology measures
5. **Temporal Dynamics**: Changes across states/time points

### Visualization Tips

```r
# Visualize DNB scores across states
plot(result_cDNB$DNB_score, 
     type = "b", 
     xlab = "State", 
     ylab = "DNB Score",
     main = "DNB Score Dynamics")

# Heatmap of gene correlations
heatmap(result_cDNB$correlation_matrix, 
        main = "Gene Correlation Network")
```

---

## 🎓 Citation

If you use EasyDNB in your research, please cite:

### Package Citation
```bibtex
@software{easydnb2023,
  author = {Liu, Zaoqu},
  title = {EasyDNB: Dynamic Network Biomarker Analysis Made Easy},
  year = {2023},
  url = {https://github.com/Zaoqu-Liu/EasyDNB},
  version = {0.1.0}
}
```

### Key References

**Conventional DNB**:
- Chen, L., Liu, R., Liu, Z. P., Li, M., & Aihara, K. (2012). Detecting early-warning signals for sudden deterioration of complex diseases by dynamical network biomarkers. *Scientific Reports*, 2, 342.

**Landscape DNB**:
- Liu, X., Wang, Y., Ji, H., et al. (2019). Detection for disease tipping points by landscape dynamic network biomarkers. *National Science Review*, 6(4), 775-785.

**Single-sample Methods**:
- Zhong, J., Liu, R., & Chen, P. (2022). The single-sample network module biomarker (sNMB) method reveals the pre-deterioration stage of disease progression. *Journal of Molecular Cell Biology*, 14(8), mjac052.
- Liu, R., Chen, P., & Chen, L. (2020). Single-sample landscape entropy reveals the imminent phase transition during disease progression. *Bioinformatics*, 36(5), 1522-1532.

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

**Copyright © 2023 Zaoqu Liu**

---

## 👥 Author & Contact

**Zaoqu Liu**  
📧 Email: liuzaoqu@163.com  
🐙 GitHub: [@Zaoqu-Liu](https://github.com/Zaoqu-Liu)

---

## 🤝 Contributing

We welcome contributions! If you'd like to contribute to EasyDNB:

1. 🍴 Fork the repository
2. 🔧 Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. 💾 Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. 📤 Push to the branch (`git push origin feature/AmazingFeature`)
5. 🎉 Open a Pull Request

### Areas for Contribution

- 🐛 Bug reports and fixes
- 📚 Documentation improvements
- ✨ New analysis methods
- 🧪 Additional test cases
- 🎨 Visualization enhancements
- 📊 Example datasets and vignettes

---

## ❓ Support & Issues

If you encounter any problems or have questions:

1. 📖 Check the [documentation](#-comprehensive-function-guide)
2. 🔍 Search [existing issues](https://github.com/Zaoqu-Liu/EasyDNB/issues)
3. 💬 Open a [new issue](https://github.com/Zaoqu-Liu/EasyDNB/issues/new)

---

## 🙏 Acknowledgments

- The DNB theory pioneered by **Prof. Luonan Chen** and collaborators
- The R community for excellent tools and packages
- All researchers who have contributed to DNB methodology development
- Users and contributors to the EasyDNB package

---

## 📈 Version History

- **v0.1.0** (2023) - Initial release
  - Core DNB analysis methods
  - Single-sample analysis functions
  - Built-in human PPI network
  - Comprehensive documentation

---
