"""
High-Resolution Cartographic Composite Map Generator for Ethiopian Highlands Soil Erosion Modeling
Produces:
- Figure 1: 6-Panel RUSLE Factor Multi-Plot Composite (300 DPI)
- Figure 2: Publication-Quality Soil Loss Severity Map with Hillshade Relief, Zonal Boundaries & Inset
- Figure 3: Multi-Criteria Erosion Susceptibility & Conservation Priority Map
"""

import os
import numpy as np
import rasterio
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
from matplotlib.patches import Rectangle, Polygon, FancyBboxPatch
import matplotlib.patheffects as pe
from shapely.geometry import box

# Configure Matplotlib for high-quality publication output
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

def load_downsampled_raster(path, max_dim=2500):
    """Loads and downsamples raster if needed for crisp memory-efficient plotting."""
    with rasterio.open(path) as src:
        scale = max(src.width // max_dim, src.height // max_dim, 1)
        out_shape = (1, src.height // scale, src.width // scale)
        data = src.read(1, out_shape=(out_shape[1], out_shape[2]), resampling=rasterio.enums.Resampling.bilinear)
        bounds = src.bounds
        crs = src.crs
        nodata = src.nodata
    return data, bounds, crs, nodata

def add_north_arrow(ax, x=0.92, y=0.92, size=0.07):
    """Adds a professional north arrow to map axes."""
    ax.annotate(
        'N',
        xy=(x, y), xytext=(x, y - size),
        xycoords='axes fraction', textcoords='axes fraction',
        arrowprops=dict(facecolor='black', edgecolor='white', lw=1.2, width=3, headwidth=10, headlength=12),
        ha='center', va='bottom', fontsize=12, fontweight='bold',
        color='black', path_effects=[pe.withStroke(linewidth=2, foreground='white')]
    )

def add_scale_bar(ax, bounds, length_km=50, loc=(0.06, 0.06)):
    """Adds a metric scale bar (km) in UTM coordinates."""
    xmin, ymin, xmax, ymax = bounds
    width_m = xmax - xmin
    height_m = ymax - ymin
    
    length_m = length_km * 1000
    sb_x0 = xmin + width_m * loc[0]
    sb_y0 = ymin + height_m * loc[1]
    sb_h = height_m * 0.015
    
    # Scale bar boxes
    ax.add_patch(Rectangle((sb_x0, sb_y0), length_m / 2, sb_h, facecolor='black', edgecolor='black', zorder=10))
    ax.add_patch(Rectangle((sb_x0 + length_m / 2, sb_y0), length_m / 2, sb_h, facecolor='white', edgecolor='black', zorder=10))
    
    # Scale labels
    ax.text(sb_x0, sb_y0 + sb_h * 1.5, '0', fontsize=8, fontweight='bold', ha='center', va='bottom',
            path_effects=[pe.withStroke(linewidth=2, foreground='white')], zorder=11)
    ax.text(sb_x0 + length_m / 2, sb_y0 + sb_h * 1.5, f'{length_km//2}', fontsize=8, fontweight='bold', ha='center', va='bottom',
            path_effects=[pe.withStroke(linewidth=2, foreground='white')], zorder=11)
    ax.text(sb_x0 + length_m, sb_y0 + sb_h * 1.5, f'{length_km} km', fontsize=8, fontweight='bold', ha='center', va='bottom',
            path_effects=[pe.withStroke(linewidth=2, foreground='white')], zorder=11)

def generate_figure1_factors_composite():
    """Generates Figure 1: 6-Panel RUSLE Factors Composite."""
    print("Generating Figure 1: RUSLE Factors Composite...")
    
    dem, bounds, crs, _ = load_downsampled_raster(r'Outputs\Rasters\dem_utm37n_30m.tif')
    slope, _, _, _ = load_downsampled_raster(r'Outputs\Rasters\slope_deg.tif')
    r_factor, _, _, _ = load_downsampled_raster(r'Outputs\Rasters\r_factor.tif')
    k_factor, _, _, _ = load_downsampled_raster(r'Outputs\Rasters\k_factor.tif')
    ls_factor, _, _, _ = load_downsampled_raster(r'Outputs\Rasters\ls_factor.tif')
    c_factor, _, _, _ = load_downsampled_raster(r'Outputs\Rasters\c_factor.tif')
    water, _, _, _ = load_downsampled_raster(r'Outputs\Rasters\water_mask.tif')
    
    extent = [bounds.left / 1000, bounds.right / 1000, bounds.bottom / 1000, bounds.top / 1000] # km
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharex=True, sharey=True)
    fig.patch.set_facecolor('#F8F9FA')
    
    # Subplot 1: DEM Elevation
    ax1 = axes[0, 0]
    im1 = ax1.imshow(np.where(dem > 0, dem, np.nan), extent=extent, cmap='terrain', vmin=1000, vmax=4100)
    ax1.set_title('(a) Elevation (SRTM DEM)', fontsize=13, fontweight='bold', pad=8)
    cb1 = plt.colorbar(im1, ax=ax1, orientation='horizontal', pad=0.06, fraction=0.046)
    cb1.set_label('Elevation (m a.s.l.)', fontsize=10, fontweight='bold')
    
    # Subplot 2: Slope Gradient
    ax2 = axes[0, 1]
    im2 = ax2.imshow(np.where(slope >= 0, slope, np.nan), extent=extent, cmap='YlOrRd', vmin=0, vmax=45)
    ax2.set_title('(b) Slope Gradient', fontsize=13, fontweight='bold', pad=8)
    cb2 = plt.colorbar(im2, ax=ax2, orientation='horizontal', pad=0.06, fraction=0.046)
    cb2.set_label('Slope (degrees)', fontsize=10, fontweight='bold')
    
    # Subplot 3: Rainfall Erosivity R-Factor
    ax3 = axes[0, 2]
    im3 = ax3.imshow(np.where(r_factor > 0, r_factor, np.nan), extent=extent, cmap='Blues', vmin=400, vmax=1100)
    ax3.set_title('(c) Rainfall Erosivity (R-Factor)', fontsize=13, fontweight='bold', pad=8)
    cb3 = plt.colorbar(im3, ax=ax3, orientation='horizontal', pad=0.06, fraction=0.046)
    cb3.set_label('R-Factor [MJ·mm/(ha·h·yr)]', fontsize=10, fontweight='bold')
    
    # Subplot 4: Soil Erodibility K-Factor
    ax4 = axes[1, 0]
    im4 = ax4.imshow(np.where(k_factor > 0, k_factor, np.nan), extent=extent, cmap='YlOrBr', vmin=0.15, vmax=0.32)
    ax4.set_title('(d) Soil Erodibility (K-Factor)', fontsize=13, fontweight='bold', pad=8)
    cb4 = plt.colorbar(im4, ax=ax4, orientation='horizontal', pad=0.06, fraction=0.046)
    cb4.set_label('K-Factor [t·ha·h/(ha·MJ·mm)]', fontsize=10, fontweight='bold')
    
    # Subplot 5: Topographic LS-Factor
    ax5 = axes[1, 1]
    im5 = ax5.imshow(np.where(ls_factor > 0, ls_factor, np.nan), extent=extent, cmap='magma', vmin=0, vmax=25)
    ax5.set_title('(e) Topographic LS-Factor', fontsize=13, fontweight='bold', pad=8)
    cb5 = plt.colorbar(im5, ax=ax5, orientation='horizontal', pad=0.06, fraction=0.046)
    cb5.set_label('LS-Factor (dimensionless)', fontsize=10, fontweight='bold')
    
    # Subplot 6: Cover Management C-Factor
    ax6 = axes[1, 2]
    c_masked = np.where(water == 1, np.nan, c_factor)
    im6 = ax6.imshow(c_masked, extent=extent, cmap='RdYlGn_r', vmin=0.01, vmax=0.60)
    ax6.set_title('(f) Cover Management (C-Factor)', fontsize=13, fontweight='bold', pad=8)
    cb6 = plt.colorbar(im6, ax=ax6, orientation='horizontal', pad=0.06, fraction=0.046)
    cb6.set_label('C-Factor (dimensionless)', fontsize=10, fontweight='bold')
    
    for ax in axes.flat:
        ax.set_facecolor('#E9ECEF')
        ax.grid(True, linestyle='--', alpha=0.4, color='#6C757D')
        ax.tick_params(labelsize=9)
        
    axes[1, 0].set_xlabel('UTM Zone 37N Easting (km)', fontsize=11, fontweight='bold')
    axes[1, 1].set_xlabel('UTM Zone 37N Easting (km)', fontsize=11, fontweight='bold')
    axes[1, 2].set_xlabel('UTM Zone 37N Easting (km)', fontsize=11, fontweight='bold')
    axes[0, 0].set_ylabel('UTM Zone 37N Northing (km)', fontsize=11, fontweight='bold')
    axes[1, 0].set_ylabel('UTM Zone 37N Northing (km)', fontsize=11, fontweight='bold')
    
    fig.suptitle('RUSLE Model Factor Decomposition — Ethiopian Highlands (Amhara Region / Lake Tana Basin)',
                 fontsize=16, fontweight='bold', y=0.98, color='#1A252C')
    
    out_fig1 = r'Outputs\Figures\Figure1_RUSLE_Factors_Composite.png'
    plt.savefig(out_fig1, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"Figure 1 saved to {out_fig1}")

def generate_figure2_soil_loss_map():
    """Generates Figure 2: Soil Loss Severity Map with Hillshade, Administrative Boundaries, Annotations & Inset."""
    print("Generating Figure 2: High-Resolution Soil Loss Severity Map...")
    
    soil_loss, bounds, crs, _ = load_downsampled_raster(r'Outputs\Rasters\soil_loss_rusle_t_ha_yr.tif')
    hillshade, _, _, _ = load_downsampled_raster(r'Outputs\Rasters\hillshade.tif')
    water, _, _, _ = load_downsampled_raster(r'Outputs\Rasters\water_mask.tif')
    
    extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
    
    # Load administrative vector boundaries
    gadm2_path = r'Boundary\gadm41_ETH_2.shp'
    gdf_zones = gpd.read_file(gadm2_path)
    gdf_zones_utm = gdf_zones.to_crs('EPSG:32637')
    
    # Clip vector boundaries to study area bbox
    bbox_geom = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
    gdf_study = gdf_zones_utm[gdf_zones_utm.intersects(bbox_geom)].copy()
    
    # Color mapping for Soil Loss Severity Classes (FAO / Ethiopian MoA):
    # Classes: Slight (<5), Low (5-11), Moderate (11-25), High (25-50), Very High (50-80), Severe (>80)
    sev_colors = [
        '#2B83BA', # Water / None
        '#2CA25F', # 1: Slight (<5)
        '#99D8C9', # 2: Low (5-11)
        '#FEE08B', # 3: Moderate (11-25)
        '#FDAE61', # 4: High (25-50)
        '#D7191C', # 5: Very High (50-80)
        '#7B1113'  # 6: Severe (>80)
    ]
    sev_bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
    sev_cmap = ListedColormap(sev_colors)
    sev_norm = BoundaryNorm(sev_bounds, sev_cmap.N)
    
    # Reclassify soil loss for display
    soil_loss_clean = np.where(soil_loss > 0, soil_loss, 0.0)
    disp_class = np.select(
        [
            water == 1,
            soil_loss_clean < 5.0,
            (soil_loss_clean >= 5.0) & (soil_loss_clean < 11.0),
            (soil_loss_clean >= 11.0) & (soil_loss_clean < 25.0),
            (soil_loss_clean >= 25.0) & (soil_loss_clean < 50.0),
            (soil_loss_clean >= 50.0) & (soil_loss_clean < 80.0),
            soil_loss_clean >= 80.0
        ],
        [0, 1, 2, 3, 4, 5, 6],
        default=0
    )
    
    fig, ax = plt.subplots(figsize=(14, 14), dpi=300)
    fig.patch.set_facecolor('white')
    
    # Base Hillshade Layer
    ax.imshow(hillshade, extent=extent, cmap='gray', vmin=0, vmax=255, alpha=0.65, zorder=1)
    
    # Soil Loss Overlay (Semi-transparent with blend)
    masked_disp = np.ma.masked_where(disp_class == 0, disp_class)
    im = ax.imshow(masked_disp, extent=extent, cmap=sev_cmap, norm=sev_norm, alpha=0.72, zorder=2)
    
    # Water Body Overlay (Lake Tana)
    water_masked = np.ma.masked_where(water == 0, water)
    ax.imshow(water_masked, extent=extent, cmap=ListedColormap(['#3182BD']), alpha=0.9, zorder=3)
    
    # Overlay Administrative Boundaries
    gdf_study.boundary.plot(ax=ax, color='#1A1A1A', linewidth=1.4, linestyle='-', zorder=4)
    
    # Annotate Zones and Major Landmarks
    landmarks = [
        ('Lake Tana', 320000, 1325000, '#084594', 12, 'bold', 'normal'),
        ('Bahir Dar', 325000, 1282000, '#000000', 10, 'bold', 'normal'),
        ('Gondar', 338000, 1390000, '#000000', 10, 'bold', 'normal'),
        ('Debre Tabor', 400000, 1312000, '#000000', 10, 'bold', 'normal'),
        ('Mount Guna\n(4,127 m)', 415000, 1290000, '#67000D', 9, 'bold', 'normal'),
        ('Blue Nile Outlet\n(Abay River)', 342000, 1270000, '#08306B', 9, 'normal', 'italic'),
    ]
    for name, lx, ly, col, fs, fw, fst in landmarks:
        if bounds.left <= lx <= bounds.right and bounds.bottom <= ly <= bounds.top:
            ax.plot(lx, ly, marker='o', markersize=5, color=col, markeredgecolor='white', markeredgewidth=1.2, zorder=6)
            ax.text(lx + 4000, ly + 2000, name, fontsize=fs, fontweight=fw, fontstyle=fst, color=col, zorder=7,
                    path_effects=[pe.withStroke(linewidth=3, foreground='white')])
            
    # Zone Labels from GADM
    for idx, row in gdf_study.iterrows():
        cent = row.geometry.centroid
        zname = row['NAME_2']
        if bounds.left + 15000 <= cent.x <= bounds.right - 15000 and bounds.bottom + 15000 <= cent.y <= bounds.top - 15000:
            ax.text(cent.x, cent.y, zname.upper(), fontsize=10, fontweight='bold', color='#2B2B2B',
                    ha='center', va='center', zorder=5, alpha=0.85,
                    path_effects=[pe.withStroke(linewidth=3, foreground='#F0F0F0')])
            
    # Set Map Extents and Graticule
    ax.set_xlim(bounds.left, bounds.right)
    ax.set_ylim(bounds.bottom, bounds.top)
    
    # Metric tick format
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{int(x/1000):,} km"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, p: f"{int(y/1000):,} km"))
    ax.tick_params(labelsize=10, colors='#333333')
    ax.grid(True, linestyle=':', alpha=0.6, color='#555555', zorder=1)
    
    ax.set_xlabel('Easting — UTM Zone 37N (m)', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_ylabel('Northing — UTM Zone 37N (m)', fontsize=12, fontweight='bold', labelpad=10)
    
    # Title Banner Box
    title_box = dict(boxstyle='round,pad=0.6', facecolor='white', edgecolor='#CCCCCC', alpha=0.92, zorder=8)
    ax.text(0.03, 0.96,
            "ESTIMATED MEAN ANNUAL SOIL LOSS (RUSLE MODEL)\n"
            "Amhara Region & Lake Tana Basin, Ethiopian Highlands",
            transform=ax.transAxes, fontsize=13, fontweight='bold', va='top', ha='left',
            color='#111111', bbox=title_box)
    
    # Cartographic Elements (Scale Bar & North Arrow)
    add_scale_bar(ax, (bounds.left, bounds.bottom, bounds.right, bounds.top), length_km=50, loc=(0.04, 0.04))
    add_north_arrow(ax, x=0.94, y=0.94, size=0.06)
    
    # Legend Elements
    labels = [
        'Water Body (Lake Tana / Reservoirs)',
        'Slight / Tolerable (< 5.0 t/ha/yr)',
        'Low (5.0 – 11.0 t/ha/yr)',
        'Moderate (11.0 – 25.0 t/ha/yr)',
        'High (25.0 – 50.0 t/ha/yr)',
        'Very High (50.0 – 80.0 t/ha/yr)',
        'Severe / Catastrophic (> 80.0 t/ha/yr)'
    ]
    patches = [Rectangle((0, 0), 1, 1, facecolor=sev_colors[i], edgecolor='black', lw=0.6) for i in range(7)]
    leg = ax.legend(patches, labels, loc='lower right', title='Soil Loss Severity Classes (FAO / MoA)',
                    title_fontsize=10, fontsize=9, framealpha=0.92, edgecolor='#CCCCCC')
    leg.get_title().set_fontweight('bold')
    leg.set_zorder(10)
    
    # Inset Map of Ethiopia
    inset_ax = fig.add_axes([0.76, 0.73, 0.20, 0.20], zorder=12)
    inset_ax.set_facecolor('#EBF4FA')
    gadm0_path = r'Boundary\gadm41_ETH_0.shp'
    gdf_eth = gpd.read_file(gadm0_path)
    gdf_eth.plot(ax=inset_ax, facecolor='#D9D9D9', edgecolor='#666666', linewidth=0.8)
    
    # Highlight study area on inset
    bbox_4326 = gpd.GeoDataFrame([1], geometry=[bbox_geom], crs='EPSG:32637').to_crs('EPSG:4326')
    bbox_4326.plot(ax=inset_ax, facecolor='red', edgecolor='darkred', alpha=0.6, linewidth=1.2)
    inset_ax.set_title('Ethiopia Context', fontsize=8, fontweight='bold', pad=3)
    inset_ax.set_xticks([])
    inset_ax.set_yticks([])
    for spine in inset_ax.spines.values():
        spine.set_edgecolor('#888888')
        
    out_fig2 = r'Outputs\Figures\Figure2_Soil_Loss_Severity_Map.png'
    plt.savefig(out_fig2, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"Figure 2 saved to {out_fig2}")

def generate_figure3_susceptibility_map():
    """Generates Figure 3: Multi-Criteria Erosion Susceptibility & Conservation Priority Map."""
    print("Generating Figure 3: Erosion Susceptibility & Priority Intervention Map...")
    
    susceptibility, bounds, crs, _ = load_downsampled_raster(r'Outputs\Rasters\erosion_susceptibility_index.tif')
    hillshade, _, _, _ = load_downsampled_raster(r'Outputs\Rasters\hillshade.tif')
    water, _, _, _ = load_downsampled_raster(r'Outputs\Rasters\water_mask.tif')
    
    extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
    
    gadm2_path = r'Boundary\gadm41_ETH_2.shp'
    gdf_zones = gpd.read_file(gadm2_path)
    gdf_zones_utm = gdf_zones.to_crs('EPSG:32637')
    bbox_geom = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
    gdf_study = gdf_zones_utm[gdf_zones_utm.intersects(bbox_geom)].copy()
    
    # Priority classification:
    # 0: Water
    # 1: Low Risk / Low Priority (0.00 - 0.25)
    # 2: Moderate Risk / Medium Priority (0.25 - 0.45)
    # 3: High Risk / High Priority (0.45 - 0.65)
    # 4: Critical Risk / Urgent Intervention (> 0.65)
    susc_clean = np.where(susceptibility > 0, susceptibility, 0.0)
    prio_class = np.select(
        [
            water == 1,
            susc_clean < 0.25,
            (susc_clean >= 0.25) & (susc_clean < 0.45),
            (susc_clean >= 0.45) & (susc_clean < 0.65),
            susc_clean >= 0.65
        ],
        [0, 1, 2, 3, 4],
        default=0
    )
    
    prio_colors = [
        '#3182BD', # Water
        '#A1D99B', # 1: Low
        '#FED976', # 2: Moderate
        '#FD8D3C', # 3: High
        '#BD0026'  # 4: Critical
    ]
    prio_cmap = ListedColormap(prio_colors)
    prio_norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], prio_cmap.N)
    
    fig, ax = plt.subplots(figsize=(14, 14), dpi=300)
    fig.patch.set_facecolor('white')
    
    ax.imshow(hillshade, extent=extent, cmap='gray', vmin=0, vmax=255, alpha=0.6, zorder=1)
    
    masked_prio = np.ma.masked_where(prio_class == 0, prio_class)
    ax.imshow(masked_prio, extent=extent, cmap=prio_cmap, norm=prio_norm, alpha=0.75, zorder=2)
    
    water_masked = np.ma.masked_where(water == 0, water)
    ax.imshow(water_masked, extent=extent, cmap=ListedColormap(['#3182BD']), alpha=0.95, zorder=3)
    
    gdf_study.boundary.plot(ax=ax, color='#1A1A1A', linewidth=1.4, zorder=4)
    
    # Annotate Zones
    for idx, row in gdf_study.iterrows():
        cent = row.geometry.centroid
        zname = row['NAME_2']
        if bounds.left + 15000 <= cent.x <= bounds.right - 15000 and bounds.bottom + 15000 <= cent.y <= bounds.top - 15000:
            ax.text(cent.x, cent.y, zname.upper(), fontsize=10, fontweight='bold', color='#222222',
                    ha='center', va='center', zorder=5, alpha=0.9,
                    path_effects=[pe.withStroke(linewidth=3, foreground='#FFFFFF')])
            
    ax.set_xlim(bounds.left, bounds.right)
    ax.set_ylim(bounds.bottom, bounds.top)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{int(x/1000):,} km"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, p: f"{int(y/1000):,} km"))
    ax.tick_params(labelsize=10)
    ax.grid(True, linestyle=':', alpha=0.5, color='#666666')
    
    ax.set_xlabel('Easting — UTM Zone 37N (m)', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_ylabel('Northing — UTM Zone 37N (m)', fontsize=12, fontweight='bold', labelpad=10)
    
    title_box = dict(boxstyle='round,pad=0.6', facecolor='white', edgecolor='#CCCCCC', alpha=0.92, zorder=8)
    ax.text(0.03, 0.96,
            "EROSION SUSCEPTIBILITY & CONSERVATION PRIORITY ZONES\n"
            "Multi-Criteria Weighted Overlay Analysis (Slope 45%, Rainfall 35%, Soil 20%)",
            transform=ax.transAxes, fontsize=13, fontweight='bold', va='top', ha='left',
            color='#111111', bbox=title_box)
    
    add_scale_bar(ax, (bounds.left, bounds.bottom, bounds.right, bounds.top), length_km=50, loc=(0.04, 0.04))
    add_north_arrow(ax, x=0.94, y=0.94, size=0.06)
    
    labels = [
        'Water Body (Lake Tana)',
        'Low Susceptibility / Maintenance Priority (0.00 – 0.25)',
        'Moderate Susceptibility / Standard Conservation (0.25 – 0.45)',
        'High Susceptibility / High Conservation Priority (0.45 – 0.65)',
        'Critical Susceptibility / Urgent Rehabilitation (> 0.65)'
    ]
    patches = [Rectangle((0, 0), 1, 1, facecolor=prio_colors[i], edgecolor='black', lw=0.6) for i in range(5)]
    leg = ax.legend(patches, labels, loc='lower right', title='Erosion Susceptibility & Intervention Tiers',
                    title_fontsize=10, fontsize=9, framealpha=0.92, edgecolor='#CCCCCC')
    leg.get_title().set_fontweight('bold')
    leg.set_zorder(10)
    
    out_fig3 = r'Outputs\Figures\Figure3_Erosion_Susceptibility_Overlay.png'
    plt.savefig(out_fig3, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"Figure 3 saved to {out_fig3}")

def main():
    os.makedirs(r'Outputs\Figures', exist_ok=True)
    generate_figure1_factors_composite()
    generate_figure2_soil_loss_map()
    generate_figure3_susceptibility_map()
    print("All high-resolution composite maps generated successfully!")

if __name__ == '__main__':
    main()
