"""
Zonal Risk Analytics and Statistical Chart Generator for Ethiopian Highlands Soil Erosion Modeling
Produces:
- Figure 4: Zonal Soil Erosion Distribution & Severity Breakdown (Original 4-panel Layout, 300 DPI)
  - (a) Mean Soil Erosion Rate by Administrative Zone (Horizontal Bar Chart with Value Annotations)
  - (b) Total Annual Sediment Yield Contribution by Zone (Mt/yr)
  - (c) Soil Loss by Agro-Ecological Elevation Belt
  - (d) Total Study Area Severity Class Distribution (Donut Chart)
- Figure 5: Erosion Dynamics, Pareto Cumulative Loss Curve, Slope & Land Cover Risk Matrix (300 DPI)
- CSV Tables: Zone_Soil_Loss_Statistics.csv, Woreda_Soil_Loss_Statistics.csv, Elevation_Soil_Loss_Statistics.csv
"""

import os
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
from shapely.geometry import box

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

def ensure_dirs():
    os.makedirs(r'Outputs\Figures', exist_ok=True)
    os.makedirs(r'Outputs\Tables', exist_ok=True)

def compute_zonal_statistics():
    """Computes zonal statistics of soil loss by Administrative Zone and Woreda using in-memory geometry masking."""
    print("Extracting zonal statistics across administrative boundaries...")
    
    with rasterio.open(r'Outputs\Rasters\soil_loss_rusle_t_ha_yr.tif') as src:
        soil_loss_arr = src.read(1)
        transform = src.transform
        bounds = src.bounds
        crs = src.crs
        res_x, res_y = src.res
        pixel_area_ha = (res_x * res_y) / 10000.0 # 30m x 30m = 0.09 ha
        out_shape = soil_loss_arr.shape
        
    with rasterio.open(r'Outputs\Rasters\water_mask.tif') as w_src:
        water_arr = w_src.read(1).astype(bool)
        
    with rasterio.open(r'Outputs\Rasters\dem_utm37n_30m.tif') as d_src:
        dem_arr = d_src.read(1)
        
    with rasterio.open(r'Outputs\Rasters\slope_deg.tif') as slp_src:
        slope_arr = slp_src.read(1)
        
    with rasterio.open(r'Outputs\Rasters\ndvi.tif') as ndvi_src:
        ndvi_arr = ndvi_src.read(1)

    # 1. ZONAL BY ADMINISTRATIVE ZONE (GADM Level 2)
    gadm2 = gpd.read_file(r'Boundary\gadm41_ETH_2.shp').to_crs(crs)
    bbox_geom = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
    gdf_zones = gadm2[gadm2.intersects(bbox_geom)].copy()
    
    zone_stats = []
    for idx, row in gdf_zones.iterrows():
        try:
            geom = [row.geometry]
            g_mask = geometry_mask(geom, transform=transform, invert=True, out_shape=out_shape)
            
            valid = g_mask & (soil_loss_arr > 0) & (soil_loss_arr < 9000) & (~water_arr)
            vals = soil_loss_arr[valid]
            
            if len(vals) > 1000:
                area_ha = len(vals) * pixel_area_ha
                total_loss_tonnes = float(np.sum(vals) * pixel_area_ha)
                total_loss_mt = total_loss_tonnes / 1e6
                mean_loss = float(np.mean(vals))
                median_loss = float(np.median(vals))
                p90_loss = float(np.percentile(vals, 90))
                
                pct_slight = float(np.mean(vals < 5.0) * 100.0)
                pct_low = float(np.mean((vals >= 5.0) & (vals < 11.0)) * 100.0)
                pct_mod = float(np.mean((vals >= 11.0) & (vals < 25.0)) * 100.0)
                pct_high = float(np.mean((vals >= 25.0) & (vals < 50.0)) * 100.0)
                pct_very_high = float(np.mean((vals >= 50.0) & (vals < 80.0)) * 100.0)
                pct_severe = float(np.mean(vals >= 80.0) * 100.0)
                
                zone_stats.append({
                    'Zone': row['NAME_2'],
                    'Region': row['NAME_1'],
                    'Area_ha': area_ha,
                    'Area_sqkm': area_ha / 100.0,
                    'Mean_Soil_Loss_t_ha_yr': mean_loss,
                    'Median_Soil_Loss_t_ha_yr': median_loss,
                    'P90_Soil_Loss_t_ha_yr': p90_loss,
                    'Total_Displacement_Mt_yr': total_loss_mt,
                    'Slight_pct': pct_slight,
                    'Low_pct': pct_low,
                    'Moderate_pct': pct_mod,
                    'High_pct': pct_high,
                    'Very_High_pct': pct_very_high,
                    'Severe_pct': pct_severe
                })
        except Exception as e:
            continue
            
    df_zones = pd.DataFrame(zone_stats).sort_values('Mean_Soil_Loss_t_ha_yr', ascending=True)
    df_zones.to_csv(r'Outputs\Tables\Zone_Soil_Loss_Statistics.csv', index=False)
    print(f"Saved Zone statistics: {len(df_zones)} zones analyzed.")

    # 2. ZONAL BY WOREDA (GADM Level 3)
    gadm3 = gpd.read_file(r'Boundary\gadm41_ETH_3.shp').to_crs(crs)
    gdf_woredas = gadm3[gadm3.intersects(bbox_geom)].copy()
    
    woreda_stats = []
    for idx, row in gdf_woredas.iterrows():
        try:
            geom = [row.geometry]
            g_mask = geometry_mask(geom, transform=transform, invert=True, out_shape=out_shape)
            valid = g_mask & (soil_loss_arr > 0) & (soil_loss_arr < 9000) & (~water_arr)
            vals = soil_loss_arr[valid]
            
            if len(vals) > 500:
                area_ha = len(vals) * pixel_area_ha
                total_loss_mt = float(np.sum(vals) * pixel_area_ha) / 1e6
                mean_loss = float(np.mean(vals))
                
                if mean_loss > 80.0:
                    priority = 'Critical / Urgent'
                elif mean_loss > 50.0:
                    priority = 'High Priority'
                elif mean_loss > 25.0:
                    priority = 'Moderate Priority'
                else:
                    priority = 'Maintenance / Low'
                    
                woreda_stats.append({
                    'Woreda': row['NAME_3'],
                    'Zone': row['NAME_2'],
                    'Region': row['NAME_1'],
                    'Area_ha': area_ha,
                    'Mean_Soil_Loss_t_ha_yr': mean_loss,
                    'Total_Displacement_Mt_yr': total_loss_mt,
                    'Intervention_Priority': priority
                })
        except Exception as e:
            continue
            
    df_woredas = pd.DataFrame(woreda_stats).sort_values('Mean_Soil_Loss_t_ha_yr', ascending=False)
    df_woredas.to_csv(r'Outputs\Tables\Woreda_Soil_Loss_Statistics.csv', index=False)
    print(f"Saved Woreda statistics: {len(df_woredas)} woredas analyzed.")

    # 3. ZONAL BY TRADITIONAL AGRO-ECOLOGICAL ELEVATION BELTS
    valid_land = (dem_arr > 0) & (soil_loss_arr > 0) & (soil_loss_arr < 9000) & (~water_arr)
    elev_vals = dem_arr[valid_land]
    sl_vals = soil_loss_arr[valid_land]
    
    elev_belts = [
        ('Kolla (Lowlands, <1,500 m)', (elev_vals < 1500)),
        ('Woina Dega (Mid-Highlands, 1,500–2,300 m)', (elev_vals >= 1500) & (elev_vals < 2300)),
        ('Dega (Highlands, 2,300–3,200 m)', (elev_vals >= 2300) & (elev_vals < 3200)),
        ('Wurch (Afro-Alpine, >3,200 m)', (elev_vals >= 3200))
    ]
    
    elev_stats = []
    tot_study_loss_mt = float(np.sum(sl_vals) * pixel_area_ha) / 1e6
    for name, mask_belt in elev_belts:
        b_sl = sl_vals[mask_belt]
        if len(b_sl) > 0:
            area_ha = len(b_sl) * pixel_area_ha
            tot_loss_mt = float(np.sum(b_sl) * pixel_area_ha) / 1e6
            elev_stats.append({
                'Elevation_Belt': name,
                'Area_ha': area_ha,
                'Area_Pct': (len(b_sl) / len(sl_vals)) * 100.0,
                'Mean_Soil_Loss_t_ha_yr': float(np.mean(b_sl)),
                'Median_Soil_Loss_t_ha_yr': float(np.median(b_sl)),
                'Total_Displacement_Mt_yr': tot_loss_mt,
                'Displacement_Pct': (tot_loss_mt / tot_study_loss_mt) * 100.0
            })
            
    df_elev = pd.DataFrame(elev_stats)
    df_elev.to_csv(r'Outputs\Tables\Elevation_Soil_Loss_Statistics.csv', index=False)
    print(f"Saved Elevation belt statistics.")
    
    return df_zones, df_woredas, df_elev, sl_vals, slope_arr[valid_land], ndvi_arr[valid_land], pixel_area_ha

def generate_figure4_zonal_charts(df_zones, df_woredas, df_elev, sl_vals, pixel_area_ha):
    """
    Generates Figure 4: Multi-Panel Zonal Risk Analysis in original clean layout.
    - (a) Mean Soil Erosion Rate by Administrative Zone (Horizontal Bar Chart)
    - (b) Total Annual Sediment Yield Contribution by Zone (Mt/yr)
    - (c) Soil Loss by Agro-Ecological Elevation Belt
    - (d) Total Study Area Severity Class Distribution (Donut Chart)
    """
    print("Generating Figure 4: Zonal Risk Analysis Charts (Original Layout)...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.patch.set_facecolor('#F8F9FA')
    
    # Subplot (a): Mean Soil Loss Rate by Administrative Zone
    ax1 = axes[0, 0]
    top_zones = df_zones.sort_values('Mean_Soil_Loss_t_ha_yr', ascending=True).reset_index(drop=True)
    y_pos = np.arange(len(top_zones))
    
    # Colors mapped directly to severity classes:
    colors_zones = [
        '#7B1113' if v >= 80 else
        '#D7191C' if v >= 50 else
        '#FDAE61' if v >= 25 else
        '#FEE08B' if v >= 11 else
        '#99D8C9' if v >= 5 else
        '#2CA25F'
        for v in top_zones['Mean_Soil_Loss_t_ha_yr']
    ]
    
    bars1 = ax1.barh(y_pos, top_zones['Mean_Soil_Loss_t_ha_yr'], color=colors_zones, edgecolor='black', height=0.62, zorder=3)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(top_zones['Zone'], fontsize=10, fontweight='bold')
    ax1.set_xlabel('Mean Annual Soil Loss (t · ha⁻¹ · yr⁻¹)', fontsize=11, fontweight='bold')
    ax1.set_title('(a) Mean Soil Erosion Rate by Administrative Zone', fontsize=12, fontweight='bold', pad=8)
    ax1.grid(True, linestyle='--', alpha=0.5, axis='x', zorder=1)
    
    for bar in bars1:
        w = bar.get_width()
        ax1.text(w + 1.5, bar.get_y() + bar.get_height()/2, f"{w:.1f} t/ha", va='center', fontsize=9, fontweight='bold')
    ax1.set_xlim(0, max(top_zones['Mean_Soil_Loss_t_ha_yr']) * 1.15)
    
    # Subplot (b): Total Annual Sediment Yield by Zone (Megatonnes / Year)
    ax2 = axes[0, 1]
    top_yield_zones = df_zones.sort_values('Total_Displacement_Mt_yr', ascending=True).reset_index(drop=True)
    y_pos2 = np.arange(len(top_yield_zones))
    
    bars2 = ax2.barh(y_pos2, top_yield_zones['Total_Displacement_Mt_yr'], color='#2B83BA', edgecolor='#1A476F', height=0.62, zorder=3)
    ax2.set_yticks(y_pos2)
    ax2.set_yticklabels(top_yield_zones['Zone'], fontsize=10, fontweight='bold')
    ax2.set_xlabel('Total Annual Soil Displacement (Megatonnes / Year)', fontsize=11, fontweight='bold')
    ax2.set_title('(b) Total Sediment Yield Contribution by Zone (Mt/yr)', fontsize=12, fontweight='bold', pad=8)
    ax2.grid(True, linestyle='--', alpha=0.5, axis='x', zorder=1)
    
    for bar in bars2:
        w = bar.get_width()
        ax2.text(w + 1.0, bar.get_y() + bar.get_height()/2, f"{w:.1f} Mt", va='center', fontsize=9, fontweight='bold')
    ax2.set_xlim(0, max(top_yield_zones['Total_Displacement_Mt_yr']) * 1.18)
    
    # Subplot (c): Soil Loss Across Traditional Elevation Belts
    ax3 = axes[1, 0]
    x_pos3 = np.arange(len(df_elev))
    belts_short = ['Kolla\n(<1,500m)', 'Woina Dega\n(1,500–2,300m)', 'Dega\n(2,300–3,200m)', 'Wurch\n(>3,200m)']
    
    w_bar = 0.35
    b3_mean = ax3.bar(x_pos3 - w_bar/2, df_elev['Mean_Soil_Loss_t_ha_yr'], width=w_bar, color='#E66101', label='Mean Loss (t/ha/yr)', edgecolor='black', zorder=3)
    b3_med = ax3.bar(x_pos3 + w_bar/2, df_elev['Median_Soil_Loss_t_ha_yr'], width=w_bar, color='#FDB863', label='Median Loss (t/ha/yr)', edgecolor='black', zorder=3)
    
    ax3.set_xticks(x_pos3)
    ax3.set_xticklabels(belts_short, fontsize=10, fontweight='bold')
    ax3.set_ylabel('Soil Loss Rate (t · ha⁻¹ · yr⁻¹)', fontsize=11, fontweight='bold')
    ax3.set_title('(c) Soil Loss by Agro-Ecological Elevation Belt', fontsize=12, fontweight='bold', pad=8)
    ax3.legend(framealpha=0.9, fontsize=9.5)
    ax3.grid(True, linestyle='--', alpha=0.5, axis='y', zorder=1)
    
    # Subplot (d): Study Area Severity Class Distribution (Donut Chart)
    ax4 = axes[1, 1]
    
    c_slight = np.sum(sl_vals < 5.0)
    c_low = np.sum((sl_vals >= 5.0) & (sl_vals < 11.0))
    c_mod = np.sum((sl_vals >= 11.0) & (sl_vals < 25.0))
    c_high = np.sum((sl_vals >= 25.0) & (sl_vals < 50.0))
    c_vh = np.sum((sl_vals >= 50.0) & (sl_vals < 80.0))
    c_sev = np.sum(sl_vals >= 80.0)
    
    sev_counts = [c_slight, c_low, c_mod, c_high, c_vh, c_sev]
    sev_labels = [
        f'Slight (<5)\n{c_slight/len(sl_vals)*100:.1f}%',
        f'Low (5–11)\n{c_low/len(sl_vals)*100:.1f}%',
        f'Moderate (11–25)\n{c_mod/len(sl_vals)*100:.1f}%',
        f'High (25–50)\n{c_high/len(sl_vals)*100:.1f}%',
        f'Very High (50–80)\n{c_vh/len(sl_vals)*100:.1f}%',
        f'Severe (>80)\n{c_sev/len(sl_vals)*100:.1f}%'
    ]
    sev_palette = ['#2CA25F', '#99D8C9', '#FEE08B', '#FDAE61', '#D7191C', '#7B1113']
    
    wedges, texts, autotexts = ax4.pie(
        sev_counts, labels=sev_labels, colors=sev_palette, autopct='%1.1f%%',
        startangle=140, pctdistance=0.78, textprops={'fontsize': 8.5, 'fontweight': 'bold'},
        wedgeprops=dict(width=0.42, edgecolor='white', linewidth=1.5)
    )
    for at in autotexts:
        at.set_fontsize(8)
        at.set_color('#111111')
        
    ax4.set_title('(d) Total Study Area Severity Class Distribution (FAO/MoA)', fontsize=12, fontweight='bold', pad=8)
    
    fig.suptitle('Zonal Soil Erosion Distribution & Severity Breakdown — Amhara Region / Lake Tana Basin',
                 fontsize=15, fontweight='bold', y=0.98, color='#1A252C')
    
    out_fig4 = r'Outputs\Figures\Figure4_Zonal_Risk_Analysis.png'
    plt.savefig(out_fig4, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"Figure 4 saved to {out_fig4}")

def generate_figure5_dynamics_charts(sl_vals, slope_vals, ndvi_vals, pixel_area_ha):
    """Generates Figure 5: Soil Erosion Dynamics, Pareto Cumulative Loss Curve, Slope & Land Cover Matrix."""
    print("Generating Figure 5: Soil Erosion Dynamics & Risk Matrix...")
    
    if len(sl_vals) > 500000:
        idx_sub = np.random.choice(len(sl_vals), size=500000, replace=False)
        sub_sl = sl_vals[idx_sub]
        sub_slope = slope_vals[idx_sub]
        sub_ndvi = ndvi_vals[idx_sub]
    else:
        sub_sl = sl_vals
        sub_slope = slope_vals
        sub_ndvi = ndvi_vals
        
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.patch.set_facecolor('#F8F9FA')
    
    # Subplot (a): Pareto Cumulative Soil Loss Distribution (Lorenz Curve)
    ax1 = axes[0, 0]
    sorted_sl = np.sort(sub_sl)[::-1]
    cum_loss = np.cumsum(sorted_sl) / np.sum(sorted_sl) * 100.0
    cum_area = np.linspace(0, 100, len(cum_loss))
    
    ax1.plot(cum_area, cum_loss, color='#D7191C', linewidth=2.8, label='Cumulative Soil Loss Curve', zorder=4)
    ax1.plot([0, 100], [0, 100], color='#666666', linestyle='--', linewidth=1.5, label='1:1 Line of Equality', zorder=2)
    
    idx_20 = int(len(cum_area) * 0.20)
    loss_at_20 = cum_loss[idx_20]
    ax1.axvline(20, color='#2B83BA', linestyle=':', linewidth=1.5, zorder=3)
    ax1.axhline(loss_at_20, color='#2B83BA', linestyle=':', linewidth=1.5, zorder=3)
    ax1.scatter([20], [loss_at_20], color='#2B83BA', s=80, zorder=5)
    ax1.annotate(
        f'Top 20% of Land Area\nAccounts for {loss_at_20:.1f}% of\nTotal Soil Loss',
        xy=(20, loss_at_20), xytext=(35, loss_at_20 - 15),
        arrowprops=dict(facecolor='#2B83BA', shrink=0.08, width=1.5, headwidth=6),
        fontsize=10, fontweight='bold', color='#1A476F',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#2B83BA', alpha=0.9)
    )
    
    ax1.set_xlabel('% of Total Study Land Area (Ranked from Highest to Lowest Loss)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Cumulative % of Total Annual Soil Loss', fontsize=11, fontweight='bold')
    ax1.set_title('(a) Disproportionate Erosion Hotspot Distribution (Pareto Curve)', fontsize=12, fontweight='bold', pad=8)
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 100)
    ax1.grid(True, linestyle='--', alpha=0.5, zorder=1)
    ax1.legend(loc='lower right', framealpha=0.9, fontsize=10)
    
    # Subplot (b): Soil Loss by Slope Gradient Classes
    ax2 = axes[0, 1]
    slope_classes = [
        ('< 5°\n(Gentle Plains)', sub_slope < 5.0),
        ('5° – 10°\n(Undulating)', (sub_slope >= 5.0) & (sub_slope < 10.0)),
        ('10° – 20°\n(Hilly)', (sub_slope >= 10.0) & (sub_slope < 20.0)),
        ('20° – 30°\n(Steep Hills)', (sub_slope >= 20.0) & (sub_slope < 30.0)),
        ('> 30°\n(Very Steep)', sub_slope >= 30.0)
    ]
    
    sc_labels = [sc[0] for sc in slope_classes]
    sc_means = [np.mean(sub_sl[sc[1]]) if np.sum(sc[1]) > 0 else 0.0 for sc in slope_classes]
    sc_medians = [np.median(sub_sl[sc[1]]) if np.sum(sc[1]) > 0 else 0.0 for sc in slope_classes]
    
    x_sc = np.arange(len(sc_labels))
    w_sc = 0.35
    ax2.bar(x_sc - w_sc/2, sc_means, width=w_sc, color='#D7191C', edgecolor='black', label='Mean Loss (t/ha/yr)', zorder=3)
    ax2.bar(x_sc + w_sc/2, sc_medians, width=w_sc, color='#FDAE61', edgecolor='black', label='Median Loss (t/ha/yr)', zorder=3)
    
    ax2.set_xticks(x_sc)
    ax2.set_xticklabels(sc_labels, fontsize=9.5, fontweight='bold')
    ax2.set_ylabel('Soil Loss Rate (t · ha⁻¹ · yr⁻¹)', fontsize=11, fontweight='bold')
    ax2.set_title('(b) Erosion Rate vs. Slope Steepness Classes', fontsize=12, fontweight='bold', pad=8)
    ax2.grid(True, linestyle='--', alpha=0.5, axis='y', zorder=1)
    ax2.legend(framealpha=0.9, fontsize=9.5)
    
    # Subplot (c): Soil Loss vs. Land Cover / NDVI Vegetation Density
    ax3 = axes[1, 0]
    ndvi_classes = [
        ('Dense Forest\n(NDVI ≥ 0.45)', sub_ndvi >= 0.45),
        ('Woodland/Shrub\n(0.30–0.45)', (sub_ndvi >= 0.30) & (sub_ndvi < 0.45)),
        ('Cropland/Mixed\n(0.15–0.30)', (sub_ndvi >= 0.15) & (sub_ndvi < 0.30)),
        ('Degraded Pasture\n(0.05–0.15)', (sub_ndvi >= 0.05) & (sub_ndvi < 0.15)),
        ('Bare Ground\n(< 0.05)', sub_ndvi < 0.05)
    ]
    
    nc_labels = [nc[0] for nc in ndvi_classes]
    nc_means = [np.mean(sub_sl[nc[1]]) if np.sum(nc[1]) > 0 else 0.0 for nc in ndvi_classes]
    nc_colors = ['#2CA25F', '#A1D99B', '#FEE08B', '#FD8D3C', '#7B1113']
    
    bars3 = ax3.bar(np.arange(len(nc_labels)), nc_means, color=nc_colors, edgecolor='black', width=0.55, zorder=3)
    ax3.set_xticks(np.arange(len(nc_labels)))
    ax3.set_xticklabels(nc_labels, fontsize=9, fontweight='bold')
    ax3.set_ylabel('Mean Soil Loss (t · ha⁻¹ · yr⁻¹)', fontsize=11, fontweight='bold')
    ax3.set_title('(c) Soil Loss by Vegetation Cover & Land Use', fontsize=12, fontweight='bold', pad=8)
    ax3.grid(True, linestyle='--', alpha=0.5, axis='y', zorder=1)
    
    for bar in bars3:
        h = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2, h + 2.0, f"{h:.1f}", ha='center', fontsize=9, fontweight='bold')
        
    # Subplot (d): 2D Risk Matrix Heatmap (Slope Class vs. Vegetation Density)
    ax4 = axes[1, 1]
    
    matrix_data = np.zeros((len(slope_classes), len(ndvi_classes)))
    for i, (_, s_mask) in enumerate(slope_classes):
        for j, (_, n_mask) in enumerate(ndvi_classes):
            cell_vals = sub_sl[s_mask & n_mask]
            matrix_data[i, j] = np.mean(cell_vals) if len(cell_vals) > 0 else 0.0
            
    im4 = ax4.imshow(matrix_data, cmap='YlOrRd', aspect='auto', vmin=0, vmax=200)
    ax4.set_xticks(np.arange(len(nc_labels)))
    ax4.set_yticks(np.arange(len(sc_labels)))
    ax4.set_xticklabels(['Forest', 'Woodland', 'Cropland', 'Pasture', 'Bare'], fontsize=10, fontweight='bold')
    ax4.set_yticklabels(['<5°', '5–10°', '10–20°', '20–30°', '>30°'], fontsize=10, fontweight='bold')
    ax4.set_xlabel('Land Cover / Vegetation Class', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Slope Gradient Class', fontsize=11, fontweight='bold')
    ax4.set_title('(d) Soil Erosion Risk Cross-Tabulation Matrix (t/ha/yr)', fontsize=12, fontweight='bold', pad=8)
    
    for i in range(len(slope_classes)):
        for j in range(len(ndvi_classes)):
            val = matrix_data[i, j]
            text_color = 'white' if val > 100 else 'black'
            ax4.text(j, i, f"{val:.0f}", ha='center', va='center', fontsize=10, fontweight='bold', color=text_color)
            
    cb4 = plt.colorbar(im4, ax=ax4, orientation='vertical', fraction=0.046, pad=0.04)
    cb4.set_label('Mean Soil Loss (t · ha⁻¹ · yr⁻¹)', fontsize=10, fontweight='bold')
    
    fig.suptitle('Soil Erosion Sensitivity Dynamics & Cross-Factor Risk Matrix',
                 fontsize=15, fontweight='bold', y=0.98, color='#1A252C')
    
    out_fig5 = r'Outputs\Figures\Figure5_Erosion_Dynamics_Charts.png'
    plt.savefig(out_fig5, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"Figure 5 saved to {out_fig5}")

def main():
    ensure_dirs()
    df_zones, df_woredas, df_elev, sl_vals, slope_vals, ndvi_vals, pixel_area_ha = compute_zonal_statistics()
    generate_figure4_zonal_charts(df_zones, df_woredas, df_elev, sl_vals, pixel_area_ha)
    generate_figure5_dynamics_charts(sl_vals, slope_vals, ndvi_vals, pixel_area_ha)
    print("All zonal charts and statistical tables generated successfully!")

if __name__ == '__main__':
    main()
