# Radiative effect of meltwater ponding on the Greenland Ice Sheet

This repository contains the code and data references for the article:

> Ryan, J. C., Cooper, M. G., Cooley, S. W., Rennermalm, Å. K., & Smith, L. C. (2025). **Meltwater ponding has an underestimated radiative effect on the surface of the Greenland Ice Sheet**. Nature Communications, XX, XX–XX. https://doi.org/XX.XXXXX  

---

## 🧊 Summary

Meltwater ponding reduces Greenland Ice Sheet albedo but this process is not included in models. This study uses drone imagery to show that small streams and ponds, often missed by satellites, considerably increases the energy available for melt.

---

![Figure 1](04-figures/fig-1-transect.png)

## 🗂 Simplified repository structure

```bash
meltwater/
├── 01-pre-processing		# Format datasets
├── 02-methods				# Classify drone orthomosaics
├── 03-analysis				# Quantify key metrics
│   ├── 01-variability
│   ├── 02-radiative-effect
│   ├── 03-scale
│   └── 04-extra
├── 04-figures				# Some figures from article
├── LICENSE					# License
└── README.md

