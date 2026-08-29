# Soil Erosion Modeling & Risk Zonation in the Ethiopian Highlands
### Quantitative RUSLE Assessment, Multi-Criteria Susceptibility & Zonal Analytics in the Amhara Region / Lake Tana Basin

---

## 📌 Executive Summary

Soil erosion by water is a severe environmental and economic threat across the **Ethiopian Highlands**, leading to extensive topsoil depletion, declining agricultural productivity, reservoir sedimentation, and ecological degradation in the Upper Blue Nile (Abay) river basin. 

This repository provides an automated, high-resolution GIS and remote sensing modeling framework implementing the **Revised Universal Soil Loss Equation (RUSLE)** coupled with **Multi-Criteria Weighted Overlay Analysis (MCDA)**. The analysis evaluates soil loss across **~48,000 km²** of the Amhara Region, encompassing the Lake Tana basin, Mount Guna massif, South Gondar, West Gojjam, and North Gondar zones at **30-meter spatial resolution**.

```
RUSLE Mathematical Formulation:
    A = R × K × LS × C × P

Where:
    A  = Mean Annual Soil Loss (t · ha⁻¹ · yr⁻¹)
    R  = Rainfall-runoff erosivity factor [MJ · mm · (ha · h · yr)⁻¹]
    K  = Soil erodibility factor [t · ha · h · (ha · MJ · mm)⁻¹]
    LS = Slope length (L) and slope steepness (S) topographic factor (dimensionless)
    C  = Cover-management factor (dimensionless, SCRP calibrated)
    P  = Support conservation practice factor (dimensionless)
```

---

## 📊 Key Findings & Quantitative Metrics

- **Average Regional Soil Loss ($A$)**: **65.2 t · ha⁻¹ · yr⁻¹** across upland terrain (median: **44.2 t · ha⁻¹ · yr⁻¹**).
- **Critical Hotspot Zones**:
  - **Debub Gondar**: **65.5 t · ha⁻¹ · yr⁻¹** (35.6% under severe risk $>80\text{ t/ha/yr}$).
  - **Semen Gondar**: **62.1 t · ha⁻¹ · yr⁻¹** (33.6% under severe risk).
  - **Wag Himra**: **93.9 t · ha⁻¹ · yr⁻¹** (54.4% under severe risk in steep gorge escarpments).
  - **Debub Wollo**: **101.7 t · ha⁻¹ · yr⁻¹** (60.5% under severe risk).
  - **Low-Risk Agricultural Basins**: **Bahir Dar Special Zone** (**3.2 t/ha/yr**; $0.0\%$ severe), **Agew Awi** (**7.8 t/ha/yr**; $1.4\%$ severe), **Mirab Gojjam** (**23.6 t/ha/yr**; $10.1\%$ severe).
- **Agro-Ecological Belts**:
  - The **Woina Dega** (Mid-Highlands, 1,500–2,300 m) accounts for the bulk of eroded sediment due to intensive rainfed teff and cereal cultivation on rolling to steep slopes.
  - The **Dega** (Highlands, 2,300–3,200 m) exhibits rapid runoff generation on steep mountain flanks.
- **Pareto Principle / Hotspot Disproportion**: The **steepest 20% of land area accounts for ~62.4% of total sediment yield**, demonstrating that targeted terracing, check dams, and watershed interventions on steep agricultural slopes will yield maximal conservation returns.

---

## 🗺️ High-Resolution Cartographic Composite Maps (300 DPI)

### Figure 1: RUSLE Model Factor Decomposition
Decomposition of the 5 input parameters ($R, K, LS, C, P$) alongside SRTM 30m elevation across the Ethiopian Highlands study area.

![Figure 1: RUSLE Model Factor Decomposition](Outputs/Figures/Figure1_RUSLE_Factors_Composite.png)

---

### Figure 2: High-Resolution Soil Loss Severity Map
Annual soil loss severity map classified according to FAO / Ethiopian Ministry of Agriculture risk classes, with 3D shaded relief hillshade blending, Lake Tana water delineation, and administrative boundary overlays.

![Figure 2: Soil Loss Severity Map](Outputs/Figures/Figure2_Soil_Loss_Severity_Map.png)

---

### Figure 3: Erosion Susceptibility & Conservation Priority Map
Multi-Criteria Weighted Overlay (Slope 45%, Rainfall 35%, Soil 20%) highlighting spatial conservation priorities from low/maintenance to urgent rehabilitation.

![Figure 3: Erosion Susceptibility & Priority Intervention Map](Outputs/Figures/Figure3_Erosion_Susceptibility_Overlay.png)

---

## 📈 Zonal Analytics & Statistical Risk Charts

### Figure 4: Zonal Soil Erosion Distribution & Severity Breakdown
Zonal analysis across 9 administrative zones, 56 woredas, and traditional agro-ecological elevation belts.

![Figure 4: Zonal Risk Analysis Charts](Outputs/Figures/Figure4_Zonal_Risk_Analysis.png)

---

### Figure 5: Soil Erosion Dynamics & Risk Matrix
Pareto cumulative loss curve, slope gradient sensitivity, vegetation density response, and 2D cross-factor risk matrix.

![Figure 5: Soil Erosion Dynamics & Risk Matrix](Outputs/Figures/Figure5_Erosion_Dynamics_Charts.png)

---

## 📁 Repository Structure

```
Erosion_Ethiopia_Highlands/
├── Boundary/               # GADM v4.1 Administrative boundaries (Levels 0, 1, 2, 3)
├── Outputs/
│   ├── Figures/            # 300 DPI publication composite maps and risk charts
│   └── Tables/             # Statistical summary CSV tables by Zone, Woreda, & Elevation
├── src/
│   ├── spatial_processing.py # Automated spatial factor & RUSLE calculation pipeline
│   ├── generate_maps.py     # Cartographic map plate generator
│   └── generate_charts.py   # Statistical analytics & chart generator
├── requirements.txt        # Python dependency specifications
├── .gitignore              # Ignores large raw satellite/DEM binaries (>100MB)
└── README.md               # Scientific documentation & methodology
```

---

## 🛠️ Usage & Reproduction

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Execute spatial processing & factor modeling (generates Outputs/Rasters/)
python src/spatial_processing.py

# 3. Generate publication-ready cartographic maps (generates Outputs/Figures/Figures 1-3)
python src/generate_maps.py

# 4. Generate zonal risk statistics & charts (generates Outputs/Figures/Figures 4-5 & Outputs/Tables/)
python src/generate_charts.py
```

---

## 📖 Methodology & Data Sources

| Factor | Primary Input Dataset | Methodological Formula / Standard |
| :--- | :--- | :--- |
| **Topography ($LS$)** | SRTM 30m DEM (`EPSG:32637`) | Moore & Burch (1986) / Desmet & Govers (1996) unit stream power formulation |
| **Rainfall ($R$)** | WorldClim v2.1 ($30''$) | Hurni (1985) / Helldén (1987) empirical model for Ethiopian Highlands: $R = 0.562 \times P + 26.3$ |
| **Soil ($K$)** | SoilGrids (ISRIC) 250m | Williams / EPIC texture model calibrated for volcanic Vertisols / Nitisols ($0.16 - 0.30$) |
| **Land Cover ($C$)** | Landsat 8/9 OLI-2 Mosaic | SCRP calibrated NDVI piecewise classification ($0.00$ water to $0.45$ bare soil) |
| **Conservation ($P$)**| 30m Slope Gradient | Wischmeier & Smith (1978) / Ethiopian Ministry of Agriculture terracing guidelines |

---

## 📜 References

1. **Hurni, H. (1985)**. *Erosion-productivity-conservation systems in Ethiopia*. In IV International Conference on Soil Conservation, Maracay, Venezuela.
2. **Nyssen, J. et al. (2004)**. *Human impact on soil erosion and reservoir sedimentation in the Ethiopian Highlands*. *Catena*, 75(1), 54-64.
3. **Bewket, W., & Teferi, E. (2009)**. *Assessment of soil erosion hazard and prioritization for treatment at the watershed level by using RUSLE and GIS in the Chemoga watershed, Ethiopia*. *Land Degradation & Development*, 20(6), 609-622.
4. **Gashaw, T. et al. (2017)**. *Soil erosion risk assessment in the Upper Blue Nile Basin, Ethiopia*. *Environmental Systems Research*, 6(1), 1-13.
