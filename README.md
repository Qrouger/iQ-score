# iQ-score

Code to compute **iQ-score** for protein–protein interaction (PPI) models generated with **AlphaFold2 (AF2)** or **AlphaFold3 (AF3)**.

This repository is designed for large-scale structural PPI analysis workflows and integrates easily with AlphaFold and AlphaPulldown pipelines.

---

## Overview

**iQ-score** is a structural confidence metric developed to evaluate predicted protein–protein interactions.

This repository provides scripts to:

- Compute iQ-score from AF2/AF3 PPI predictions
- Parse structural prediction outputs

---

# Installation

## 1. Clone the repository

```bash
git clone https://github.com/Qrouger/iQ-score
```

---

# Requirements

## Input structure

The scripts assume:

- `.a3m` file of each protein is available
- Each interaction has its own directory
- Directory names follow the format:

```text
prot1_and_prot2
```
---

# CCP4 Dependency

This project requires the **CCP4** suite.

## Download CCP4

https://www.ccp4.ac.uk/download/#os=linux

## Linux installation example

```bash
tar xvzf ccp4-9-setup.tar.gz
./ccp4-9-setup
```
---

# Usage

Basic example:

```bash
python iq_score.py \
    --msa_dir path/to/msa \
    --model_dir path/to/interactions \
    --N_CPU Number of CPU \
    --Path_ccp4 path/to/ccp4
```

---


---

# Citation

```text
Quentin Rouger, Emmanuel Giudice, Damien F Meyer, Kévin Macé,
PPIFold: a tool for analysis of protein–protein interaction from AlphaPullDown,
Bioinformatics Advances, Volume 5, Issue 1, 2025, vbaf090,
https://doi.org/10.1093/bioadv/vbaf090
```
