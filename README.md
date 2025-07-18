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
├── 01-pre-processing
├── 02-methods
├── 03-analysis	
│   ├── 01-variability
│   ├── 02-radiative-effect
│   ├── 03-scale
│   └── 04-extra
├── 04-figures
├── LICENSE
└── README.md

## Data availability

Data required to reproduce the findings of this study are available at the [Duke University Libraries Digital Repository](https://doi.org/10.7924/r4ff41j34). 

These datasets include 1) surface water and surface albedo data from satellite remote sensing, 2) near-surface air temperature and downward shortwave radiation data from atmospheric reanalysis, and 3) high-resolution imagery from done surveys over the Greenland Ice Sheet.

The surface water data were derived from Sentinel-2 by [Zhang et al. (2023)](https://www.sciencedirect.com/science/article/abs/pii/S0034425723003322). The albedo data are from the MCD43A3 version 6.1 product derived from the MODerate resolution Imaging Spectroradiometer (MODIS) onboard NASA’s Aqua and Terra satellites. The data represent two distinct time periods: Jul 25–30, 2018 and Jul 29–Aug 5, 2019 when most of the imagery was acquired. We resampled both surface water and surface albedo data onto the ISMIP6 grid by averaging values within each target grid cell. The ISMIP6 grid has an NSIDC Sea Ice Polar Stereographic North (EPSG:3413) projection and a grid cell resolution of 1 × 1 km.

The near-surface air temperature and downward all-sky shortwave radiation data are derived from the MERRA-2 radiation diagnostics product (M2T1NXRAD). The data represent summer means for the 2002-2023 period. We statistically downscaled both variables from 0.625 x 0.5° resolution to the ISMIP6 1 × 1 km grid following an elevation-based approach [(Ryan et al., 2024)](https://www.nature.com/articles/s43247-024-01714-y).

The high-resolution imagery was acquired using a fixed-wing drone similar to that described by [Ryan et al. (2015)](https://tc.copernicus.org/articles/9/1/2015/tc-9-1-2015.html). The drone collected imagery at two separate field sites. The first site includes both Russell Glacier and Isunguata Sermia (180 km2) with an elevational range of 150–660 m a.s.l. During six surveys, the drone collected 9,732 overlapping images of glacier surface on Jul 11–12, 2015. The second site (110 km2) is situated in the dark zone with an elevational range of 1,170–1,290 m a.s.l. Here the drone collected a total of 3,795 images during three surveys between Jul 20–22, 2015. We used Agisoft Metashape Pro v2.2.0 to generate orthomosaics from the overlapping aerial imagery. The final products have a spatial resolution of 0.30 m.

We mapped meltwater ponding in the orthomosaic produced over the dark zone (1,170–1,290 m a.s.l.) using a semi-automated classification approach described in the article. We mapped meltwater over Russell Glacier and Isunguata Sermia by manually digitized ponded areas using the Geo-SAM QGIS plugin, an interactive segmentation tool based on the Segment Anything Model (SAM) foundation AI model. The final water maps are provided as shapefiles.

Ryan, J. C. et al. (2015), UAV photogrammetry and structure from motion to assess calving dynamics at Store Glacier, a large outlet draining the Greenland ice sheet. The Cryosphere 9, 1–11.

Ryan, J. C. (2024), Contribution of surface and cloud radiative feedbacks to Greenland Ice Sheet meltwater production during 2002–2023. Commun. Earth Environ. 5, 1–9.

Zhang, W. et al. (2023), Pan-Greenland mapping of supraglacial rivers, lakes, and water-filled crevasses in a cool summer (2018) and a warm summer (2019). Remote Sens. Environ. 297, 113781.
