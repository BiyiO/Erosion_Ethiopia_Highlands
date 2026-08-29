"""
Publication-Quality Cartographic Map & Analytics Generator for Sediment Delivery Ratio (SDR)
Produces:
- Figure 6: High-Resolution 4-Panel Cartographic Composite (Stream Network, IC, SDR, Net Sediment Export) (300 DPI)
- Figure 7: Lake Tana Tributary Siltation & Sub-Basin Analytics Charts (300 DPI)
"""

import os
import numpy as np
import pandas as pd
import rasterio
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
import matplotlib.patheffects as pe
from shapely.geometry import box

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

def ensure_dirs():
    os.makedirs(r'Outputs\Figures', exist_ok=True)

def read_and_downsample(filepath, target_shape=(1850, 1820)):
    """Reads and smoothly downsamples rasters for fast 300 DPI cartographic rendering."""
    with rasterio.open(filepath) as src:
        data = src.read(1, out_shape=target_shape, resampling=rasterio.enums.Resampling.bilinear)
        bounds = src.bounds
        crs = src.crs
    return data, bounds, crs

def add_north_arrow(ax, x=0.92, y=0.92, size=0.07):
    """Draws a crisp professional cartographic north arrow."""
    ax.annotate(
        'N', xy=(x, y), xytext=(x, y - size),
        xycoords='axes fraction', textcoords='axes fraction',
        arrowprops=dict(facecolor='black', edgecolor='white', lw=1.2, width=3, headwidth=10, headlength=12),
        ha='center', va='bottom', fontsize=12, fontweight='bold',
        color='black', path_effects=[pe.withStroke(linewidth=2, foreground='white')]
    )

def add_scale_bar(ax, extent_km, length_km=40, loc=(0.06, 0.06)):
    """Adds a metric scale bar (km)."""
    xmin, xmax, ymin, ymax = extent_km[0], extent_km[1], extent_km[2], extent_km[3]
    width_km = xmax - xmin
    height_km = ymax - ymin
    
    sb_x0 = xmin + width_km * loc[0]
    sb_y0 = ymin + height_km * loc[1]
    sb_h = height_km * 0.015
    
    ax.add_patch(Rectangle((sb_x0, sb_y0), length_km / 2, sb_h, facecolor='black', edgecolor='black', zorder=10))
    ax.add_patch(Rectangle((sb_x0 + length_km / 2, sb_y0), length_km / 2, sb_h, facecolor='white', edgecolor='black', zorder=10))
    
    ax.text(sb_x0, sb_y0 + sb_h * 1.5, '0', fontsize=8, fontweight='bold', ha='center', va='bottom',
            path_effects=[pe.withStroke(linewidth=2, foreground='white')], zorder=11)
    ax.text(sb_x0 + length_km / 2, sb_y0 + sb_h * 1.5, f'{length_km//2}', fontsize=8, fontweight='bold', ha='center', va='bottom',
            path_effects=[pe.withStroke(linewidth=2, foreground='white')], zorder=11)
    ax.text(sb_x0 + length_km, sb_y0 + sb_h * 1.5, f'{length_km} km', fontsize=8, fontweight='bold', ha='center', va='bottom',
            path_effects=[pe.withStroke(linewidth=2, foreground='white')], zorder=11)

def generate_figure6_sdr_composite():
    """Generates Figure 6: 4-Panel Cartographic Composite for Stream Network, IC, SDR, and Net Export."""
    print("Generating Figure 6: High-Resolution SDR Composite Map...")
    
    # Load rasters
    hs_data, bounds, crs = read_and_downsample(r'Outputs\Rasters\hillshade.tif')
    water_data, _, _ = read_and_downsample(r'Outputs\Rasters\water_mask.tif')
    stream_data, _, _ = read_and_downsample(r'Outputs\Rasters\stream_network_30m.tif')
    ic_data, _, _ = read_and_downsample(r'Outputs\Rasters\sediment_connectivity_index_ic.tif')
    sdr_data, _, _ = read_and_downsample(r'Outputs\Rasters\sediment_delivery_ratio_sdr.tif')
    exp_data, _, _ = read_and_downsample(r'Outputs\Rasters\sediment_export_t_ha_yr.tif')
    
    # Load administrative boundaries
    gadm2 = gpd.read_file(r'Boundary\gadm41_ETH_2.shp').to_crs(crs)
    bbox_geom = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
    gdf_study = gadm2[gadm2.intersects(bbox_geom)].copy()
    
    extent = [bounds.left/1000, bounds.right/1000, bounds.bottom/1000, bounds.top/1000] # in km
    
    fig, axes = plt.subplots(2, 2, figsize=(18, 16))
    fig.patch.set_facecolor('#F8F9FA')
    
    water_masked = np.ma.masked_where(water_data == 0, water_data)
    
    # -------------------------------------------------------------
    # Panel (a): Topographic Stream Network & Lake Tana Basin
    # -------------------------------------------------------------
    ax1 = axes[0, 0]
    ax1.imshow(hs_data, cmap='gray', extent=extent, alpha=0.9, origin='upper')
    ax1.imshow(water_masked, cmap=mcolors.ListedColormap(['#1F78B4']), extent=extent, alpha=0.95, origin='upper')
    
    stream_masked = np.ma.masked_where(stream_data == 0, stream_data)
    ax1.imshow(stream_masked, cmap=mcolors.ListedColormap(['#08519C']), extent=extent, alpha=0.85, origin='upper')
    
    for _, row in gdf_study.iterrows():
        g_proj = row.geometry
        if g_proj.geom_type == 'Polygon':
            coords = np.array(g_proj.exterior.coords) / 1000
            ax1.plot(coords[:, 0], coords[:, 1], color='#333333', linewidth=0.9, linestyle='--')
        elif g_proj.geom_type == 'MultiPolygon':
            for poly in g_proj.geoms:
                coords = np.array(poly.exterior.coords) / 1000
                ax1.plot(coords[:, 0], coords[:, 1], color='#333333', linewidth=0.9, linestyle='--')
                
    ax1.text(350, 1335, 'LAKE\nTANA', fontsize=12, fontweight='bold', color='#FFFFFF', ha='center', va='center',
             bbox=dict(boxstyle='circle,pad=0.4', facecolor='#1F78B4', edgecolor='white', alpha=0.85))
    ax1.text(320, 1260, 'Gilgel Abay\nBasin', fontsize=9, fontweight='bold', color='#1A252C', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.75))
    ax1.text(395, 1300, 'Gumara\nBasin', fontsize=9, fontweight='bold', color='#1A252C', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.75))
    ax1.text(395, 1360, 'Ribb\nBasin', fontsize=9, fontweight='bold', color='#1A252C', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.75))
    ax1.text(345, 1400, 'Megech\nBasin', fontsize=9, fontweight='bold', color='#1A252C', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.75))
    
    ax1.set_title('(a) Topographic Drainage Channels & Lake Tana Tributaries', fontsize=12, fontweight='bold', pad=8)
    ax1.set_xlabel('UTM Easting (km)', fontsize=10, fontweight='bold')
    ax1.set_ylabel('UTM Northing (km)', fontsize=10, fontweight='bold')
    add_north_arrow(ax1)
    add_scale_bar(ax1, extent)
    
    # -------------------------------------------------------------
    # Panel (b): Index of Sediment Connectivity (IC)
    # -------------------------------------------------------------
    ax2 = axes[0, 1]
    ax2.imshow(hs_data, cmap='gray', extent=extent, alpha=0.65, origin='upper')
    
    ic_plot = np.where((ic_data == -9999) | (water_data == 1), np.nan, ic_data)
    im2 = ax2.imshow(ic_plot, cmap='coolwarm', vmin=-4.5, vmax=1.5, extent=extent, alpha=0.75, origin='upper')
    ax2.imshow(water_masked, cmap=mcolors.ListedColormap(['#1F78B4']), extent=extent, alpha=0.95, origin='upper')
    
    for _, row in gdf_study.iterrows():
        g_proj = row.geometry
        if g_proj.geom_type == 'Polygon':
            coords = np.array(g_proj.exterior.coords) / 1000
            ax2.plot(coords[:, 0], coords[:, 1], color='#333333', linewidth=0.7, linestyle=':')
        elif g_proj.geom_type == 'MultiPolygon':
            for poly in g_proj.geoms:
                coords = np.array(poly.exterior.coords) / 1000
                ax2.plot(coords[:, 0], coords[:, 1], color='#333333', linewidth=0.7, linestyle=':')
                
    cb2 = plt.colorbar(im2, ax=ax2, orientation='horizontal', fraction=0.046, pad=0.08, aspect=30)
    cb2.set_label('Index of Sediment Connectivity (IC) [Log₁₀(D_up / D_dn)]', fontsize=9.5, fontweight='bold')
    ax2.set_title('(b) Borselli Index of Sediment Connectivity (IC)', fontsize=12, fontweight='bold', pad=8)
    ax2.set_xlabel('UTM Easting (km)', fontsize=10, fontweight='bold')
    ax2.set_ylabel('UTM Northing (km)', fontsize=10, fontweight='bold')
    add_north_arrow(ax2)
    add_scale_bar(ax2, extent)
    
    # -------------------------------------------------------------
    # Panel (c): Sediment Delivery Ratio (SDR)
    # -------------------------------------------------------------
    ax3 = axes[1, 0]
    ax3.imshow(hs_data, cmap='gray', extent=extent, alpha=0.65, origin='upper')
    
    sdr_plot = np.where((sdr_data <= 0) | (water_data == 1), np.nan, sdr_data)
    im3 = ax3.imshow(sdr_plot, cmap='YlOrRd', vmin=0.0, vmax=0.40, extent=extent, alpha=0.75, origin='upper')
    ax3.imshow(water_masked, cmap=mcolors.ListedColormap(['#1F78B4']), extent=extent, alpha=0.95, origin='upper')
    
    for _, row in gdf_study.iterrows():
        g_proj = row.geometry
        if g_proj.geom_type == 'Polygon':
            coords = np.array(g_proj.exterior.coords) / 1000
            ax3.plot(coords[:, 0], coords[:, 1], color='#333333', linewidth=0.7, linestyle=':')
        elif g_proj.geom_type == 'MultiPolygon':
            for poly in g_proj.geoms:
                coords = np.array(poly.exterior.coords) / 1000
                ax3.plot(coords[:, 0], coords[:, 1], color='#333333', linewidth=0.7, linestyle=':')
                
    cb3 = plt.colorbar(im3, ax=ax3, orientation='horizontal', fraction=0.046, pad=0.08, aspect=30)
    cb3.set_label('Sediment Delivery Ratio (SDR) [Fraction: 0.0 – 0.80]', fontsize=9.5, fontweight='bold')
    ax3.set_title('(c) Calibrated Sediment Delivery Ratio (SDR)', fontsize=12, fontweight='bold', pad=8)
    ax3.set_xlabel('UTM Easting (km)', fontsize=10, fontweight='bold')
    ax3.set_ylabel('UTM Northing (km)', fontsize=10, fontweight='bold')
    add_north_arrow(ax3)
    add_scale_bar(ax3, extent)
    
    # -------------------------------------------------------------
    # Panel (d): Net Annual Sediment Export (E = A * SDR)
    # -------------------------------------------------------------
    ax4 = axes[1, 1]
    ax4.imshow(hs_data, cmap='gray', extent=extent, alpha=0.65, origin='upper')
    
    exp_plot = np.where((exp_data <= 0) | (water_data == 1), np.nan, exp_data)
    
    cmap_exp = mcolors.ListedColormap(['#2CA25F', '#99D8C9', '#FEE08B', '#FDAE61', '#D7191C', '#7B1113'])
    bounds_exp = [0, 2, 5, 10, 20, 40, 150]
    norm_exp = mcolors.BoundaryNorm(bounds_exp, cmap_exp.N)
    
    im4 = ax4.imshow(exp_plot, cmap=cmap_exp, norm=norm_exp, extent=extent, alpha=0.78, origin='upper')
    ax4.imshow(water_masked, cmap=mcolors.ListedColormap(['#1F78B4']), extent=extent, alpha=0.95, origin='upper')
    
    for _, row in gdf_study.iterrows():
        g_proj = row.geometry
        if g_proj.geom_type == 'Polygon':
            coords = np.array(g_proj.exterior.coords) / 1000
            ax4.plot(coords[:, 0], coords[:, 1], color='#333333', linewidth=0.7, linestyle=':')
        elif g_proj.geom_type == 'MultiPolygon':
            for poly in g_proj.geoms:
                coords = np.array(poly.exterior.coords) / 1000
                ax4.plot(coords[:, 0], coords[:, 1], color='#333333', linewidth=0.7, linestyle=':')
                
    cb4 = plt.colorbar(im4, ax=ax4, orientation='horizontal', fraction=0.046, pad=0.08, aspect=30, ticks=bounds_exp)
    cb4.set_label('Net Sediment Export (E) [t · ha⁻¹ · yr⁻¹]', fontsize=9.5, fontweight='bold')
    ax4.set_title('(d) Net Annual Sediment Export into Streams & Lake Tana', fontsize=12, fontweight='bold', pad=8)
    ax4.set_xlabel('UTM Easting (km)', fontsize=10, fontweight='bold')
    ax4.set_ylabel('UTM Northing (km)', fontsize=10, fontweight='bold')
    add_north_arrow(ax4)
    add_scale_bar(ax4, extent)
    
    fig.suptitle('Sediment Connectivity & Delivery Ratio (SDR) Modeling — Lake Tana Basin / Ethiopian Highlands',
                 fontsize=15, fontweight='bold', y=0.98, color='#1A252C')
    
    out_fig6 = r'Outputs\Figures\Figure6_Sediment_Delivery_Ratio_Composite.png'
    plt.savefig(out_fig6, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"Figure 6 saved to {out_fig6}")

def generate_figure7_siltation_charts():
    """Generates Figure 7: Lake Tana Siltation & Tributary Sediment Analytics."""
    print("Generating Figure 7: Lake Tana Siltation Analytics Charts...")
    
    df_sub = pd.read_csv(r'Outputs\Tables\Lake_Tana_Subcatchment_Sediment_Inflow.csv')
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.patch.set_facecolor('#F8F9FA')
    
    # Subplot (a): Annual Sediment Inflow by Major Tributary (Mt/yr)
    ax1 = axes[0, 0]
    df_sort = df_sub.sort_values('Net_Sediment_Export_Mt_yr', ascending=True).reset_index(drop=True)
    y_pos = np.arange(len(df_sort))
    
    bars1 = ax1.barh(y_pos, df_sort['Net_Sediment_Export_Mt_yr'], color='#D7191C', edgecolor='#7B1113', height=0.58, zorder=3)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(df_sort['Tributary_Subcatchment'], fontsize=9.5, fontweight='bold')
    ax1.set_xlabel('Annual Sediment Delivery into Lake Tana (Megatonnes / Year)', fontsize=10.5, fontweight='bold')
    ax1.set_title('(a) Tributary Sediment Inflow into Lake Tana (Mt/yr)', fontsize=12, fontweight='bold', pad=8)
    ax1.grid(True, linestyle='--', alpha=0.5, axis='x', zorder=1)
    
    for bar in bars1:
        w = bar.get_width()
        ax1.text(w + 0.3, bar.get_y() + bar.get_height()/2, f"{w:.2f} Mt/yr", va='center', fontsize=9, fontweight='bold')
    ax1.set_xlim(0, max(df_sort['Net_Sediment_Export_Mt_yr']) * 1.22)
    
    # Subplot (b): Sub-Basin Mean Sediment Export Rate (t/ha/yr)
    ax2 = axes[0, 1]
    df_rate = df_sub.sort_values('Mean_Sediment_Export_t_ha_yr', ascending=True).reset_index(drop=True)
    y_pos2 = np.arange(len(df_rate))
    
    bars2 = ax2.barh(y_pos2, df_rate['Mean_Sediment_Export_t_ha_yr'], color='#FDAE61', edgecolor='#E66101', height=0.58, zorder=3)
    ax2.set_yticks(y_pos2)
    ax2.set_yticklabels(df_rate['Tributary_Subcatchment'], fontsize=9.5, fontweight='bold')
    ax2.set_xlabel('Mean Specific Sediment Export (t · ha⁻¹ · yr⁻¹)', fontsize=10.5, fontweight='bold')
    ax2.set_title('(b) Watershed Sediment Yield Intensity (t/ha/yr)', fontsize=12, fontweight='bold', pad=8)
    ax2.grid(True, linestyle='--', alpha=0.5, axis='x', zorder=1)
    
    for bar in bars2:
        w = bar.get_width()
        ax2.text(w + 0.4, bar.get_y() + bar.get_height()/2, f"{w:.1f} t/ha", va='center', fontsize=9, fontweight='bold')
    ax2.set_xlim(0, max(df_rate['Mean_Sediment_Export_t_ha_yr']) * 1.20)
    
    # Subplot (c): Lake Tana Siltation Contribution Breakdown (Donut)
    ax3 = axes[1, 0]
    colors_pie = ['#2B83BA', '#D7191C', '#FDAE61', '#ABD9E9', '#2CA25F']
    labels_pie = [f"{r['Tributary_Subcatchment'].split('(')[0].strip()}\n({r['Lake_Inflow_Contribution_Pct']:.1f}%)" for _, r in df_sub.iterrows()]
    
    wedges, texts, autotexts = ax3.pie(
        df_sub['Net_Sediment_Export_Mt_yr'], labels=labels_pie, colors=colors_pie[:len(df_sub)],
        autopct='%1.1f%%', startangle=140, pctdistance=0.78,
        textprops={'fontsize': 8.5, 'fontweight': 'bold'},
        wedgeprops=dict(width=0.42, edgecolor='white', linewidth=1.5)
    )
    for at in autotexts:
        at.set_fontsize(8)
        at.set_color('#111111')
        
    ax3.set_title('(c) Proportional Tributary Inflow Contribution to Lake Tana', fontsize=12, fontweight='bold', pad=8)
    
    # Subplot (d): Theoretical Sediment Delivery vs. Distance Decay Curve
    ax4 = axes[1, 1]
    dist_km = np.linspace(0.05, 25.0, 200)
    sdr_decay_steep = 0.80 * np.exp(-dist_km / 3.5)
    sdr_decay_mod = 0.55 * np.exp(-dist_km / 2.0)
    sdr_decay_gentle = 0.25 * np.exp(-dist_km / 1.0)
    
    ax4.plot(dist_km, sdr_decay_steep * 100, color='#D7191C', linewidth=2.5, label='Steep Cultivated Slopes (>20°)', zorder=4)
    ax4.plot(dist_km, sdr_decay_mod * 100, color='#FDAE61', linewidth=2.5, label='Rolling Cropland (10–20°)', zorder=3)
    ax4.plot(dist_km, sdr_decay_gentle * 100, color='#2CA25F', linewidth=2.5, label='Forested / Riparian Buffer (<10°)', zorder=3)
    
    ax4.axvspan(0, 0.5, color='#99D8C9', alpha=0.35, label='Riparian Filter Strip (<500 m)')
    ax4.set_xlabel('Overland Flow Distance to Nearest Stream Channel (km)', fontsize=10.5, fontweight='bold')
    ax4.set_ylabel('Sediment Delivery Ratio (SDR %)', fontsize=10.5, fontweight='bold')
    ax4.set_title('(d) Sediment Delivery Attenuation vs. Stream Distance & Buffers', fontsize=12, fontweight='bold', pad=8)
    ax4.set_xlim(0, 25)
    ax4.set_ylim(0, 85)
    ax4.grid(True, linestyle='--', alpha=0.5, zorder=1)
    ax4.legend(loc='upper right', framealpha=0.92, fontsize=9)
    
    fig.suptitle('Lake Tana Basin Siltation Dynamics & Sub-Catchment Sediment Yield',
                 fontsize=15, fontweight='bold', y=0.98, color='#1A252C')
    
    out_fig7 = r'Outputs\Figures\Figure7_Lake_Tana_Siltation_Analytics.png'
    plt.savefig(out_fig7, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"Figure 7 saved to {out_fig7}")

def main():
    ensure_dirs()
    generate_figure6_sdr_composite()
    generate_figure7_siltation_charts()
    print("All SDR maps and siltation figures generated successfully!")

if __name__ == '__main__':
    main()
