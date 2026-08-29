"""
Sediment Delivery Ratio (SDR) & Stream Sediment Connectivity Pipeline for the Ethiopian Highlands
(Lake Tana Basin / Upper Blue Nile / Amhara Region)

Implements:
1. Topographic Stream Network Extraction & Lake Tana Sink Identification
2. Upslope Sediment Generation Potential (D_up) & Flow Path Impedance (D_dn)
3. Borselli / InVEST Index of Sediment Connectivity (IC)
4. Sediment Delivery Ratio (SDR) via Calibrated Sigmoidal Function
5. Net Annual Sediment Export (E = A * SDR) into Streams and Lake Tana (t/ha/yr)
6. Sub-Catchment Sediment Yield Extraction for Major Lake Tana Tributaries:
   - Gilgel Abay (Little Blue Nile / South)
   - Gumara River (South-East / Mt. Guna)
   - Ribb River (East / Debub Gondar)
   - Megech River (North / Gondar)
   - Direct Lakeshore Littoral Catchments
"""

import os
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from scipy.ndimage import distance_transform_edt, uniform_filter, gaussian_filter
import geopandas as gpd
from shapely.geometry import box, Point

def ensure_dirs():
    os.makedirs(r'Outputs\Rasters', exist_ok=True)
    os.makedirs(r'Outputs\Figures', exist_ok=True)
    os.makedirs(r'Outputs\Tables', exist_ok=True)

def compute_sdr_model():
    """
    Computes Stream Network, Index of Sediment Connectivity (IC), 
    Sediment Delivery Ratio (SDR), and Net Sediment Export (E).
    """
    print("Executing Sediment Delivery Ratio (SDR) & Connectivity Pipeline...")
    
    # 1. Load Base Rasters
    with rasterio.open(r'Outputs\Rasters\dem_utm37n_30m.tif') as dem_src:
        dem = dem_src.read(1).astype(np.float32)
        meta = dem_src.meta.copy()
        bounds = dem_src.bounds
        crs = dem_src.crs
        res_x, res_y = dem_src.res
        nodata = dem_src.nodata
        
    with rasterio.open(r'Outputs\Rasters\slope_deg.tif') as slp_src:
        slope_deg = slp_src.read(1).astype(np.float32)
        
    with rasterio.open(r'Outputs\Rasters\c_factor.tif') as c_src:
        c_factor = c_src.read(1).astype(np.float32)
        
    with rasterio.open(r'Outputs\Rasters\water_mask.tif') as w_src:
        water_mask = w_src.read(1).astype(bool)
        
    with rasterio.open(r'Outputs\Rasters\soil_loss_rusle_t_ha_yr.tif') as sl_src:
        soil_loss = sl_src.read(1).astype(np.float32)
        
    valid_land = (dem > 0) & (dem != nodata) & (slope_deg >= 0) & (~water_mask)
    pixel_area_ha = (res_x * res_y) / 10000.0 # 0.09 ha
    
    slope_rad = np.radians(np.maximum(0.1, slope_deg))
    slope_m_m = np.tan(slope_rad)
    slope_m_m = np.clip(slope_m_m, 0.005, 2.0)
    
    # 2. Extract Topographic Drainage / Stream Network
    # Valley convergence / hollow identification via multi-scale topographic curvature & slope
    print("Delineating topographic stream channels and sediment sinks...")
    gy, gx = np.gradient(dem, res_y, res_x)
    curv_y, _ = np.gradient(gy, res_y, res_x)
    _, curv_x = np.gradient(gx, res_y, res_x)
    prof_curv = curv_x + curv_y
    
    # Stream network: concentrated drainage paths in valleys + Lake Tana boundary
    valley_drainage = (prof_curv > 0.003) & (slope_deg < 15.0)
    smooth_drainage = uniform_filter(valley_drainage.astype(np.float32), size=5) > 0.35
    
    # Sinks = Lake Tana open water + main tributary channels
    sink_mask = water_mask | smooth_drainage
    
    # 3. Downslope Flow Path Impedance (D_dn)
    # D_dn = sum( d_i / (W_i * S_i) )
    # Land cover weighting factor W_i (C-factor based impedance: bare ground transmits sediment faster; dense veg traps it)
    print("Computing downslope flow impedance (D_dn)...")
    W_factor = np.clip(c_factor, 0.005, 0.50)
    travel_impedance = 1.0 / (W_factor * slope_m_m + 1e-4)
    travel_impedance = np.where(valid_land, travel_impedance, 1e4)
    
    # Distance to nearest stream sink
    dist_to_sink_px = distance_transform_edt(~sink_mask)
    dist_to_sink_m = dist_to_sink_px * 30.0
    
    # Downslope cost path proxy (distance-weighted cumulative impedance)
    smooth_impedance = uniform_filter(travel_impedance, size=7)
    D_dn = dist_to_sink_m * smooth_impedance
    D_dn = np.maximum(D_dn, 1.0)
    
    # 4. Upslope Sediment Generation Potential (D_up)
    # D_up = W_mean * S_mean * sqrt(A_contrib)
    print("Computing upslope sediment generation potential (D_up)...")
    # Contributing area proxy from topographic convergence & local upstream slope
    upslope_area_m2 = np.clip((dist_to_sink_m * 25.0) + (res_x * res_y), 900.0, 5e7)
    W_mean = uniform_filter(W_factor, size=9)
    S_mean = uniform_filter(slope_m_m, size=9)
    
    D_up = W_mean * S_mean * np.sqrt(upslope_area_m2)
    D_up = np.maximum(D_up, 1e-3)
    
    # 5. Index of Sediment Connectivity (IC)
    print("Calculating Index of Sediment Connectivity (IC)...")
    IC = np.log10(D_up / D_dn)
    IC = np.where(valid_land, IC, np.nan)
    IC = np.clip(IC, -4.5, 3.5) # Standard Borselli IC range
    
    # 6. Sediment Delivery Ratio (SDR)
    # SDR = SDR_max / (1 + exp( (IC_0 - IC) / k ))
    print("Calculating Sediment Delivery Ratio (SDR)...")
    SDR_max = 0.80
    IC_0 = 0.50
    k_param = 2.0
    
    SDR = SDR_max / (1.0 + np.exp((IC_0 - IC) / k_param))
    SDR = np.where(valid_land, SDR, 0.0)
    SDR = np.where(water_mask, 0.0, np.clip(SDR, 0.0, SDR_max))
    
    # 7. Net Annual Sediment Export (E = A * SDR in t/ha/yr)
    print("Calculating Net Annual Sediment Export (E = A * SDR)...")
    sediment_export = np.where(valid_land, soil_loss * SDR, 0.0)
    sediment_export = np.nan_to_num(sediment_export, nan=0.0)
    sediment_export = np.clip(sediment_export, 0.0, 150.0)
    
    # 8. Save Generated SDR Rasters
    meta.update(dtype='float32', nodata=-9999.0)
    
    ic_path = r'Outputs\Rasters\sediment_connectivity_index_ic.tif'
    with rasterio.open(ic_path, 'w', **meta) as dst:
        dst.write(np.nan_to_num(IC, nan=-9999.0).astype(np.float32), 1)
        
    sdr_path = r'Outputs\Rasters\sediment_delivery_ratio_sdr.tif'
    with rasterio.open(sdr_path, 'w', **meta) as dst:
        dst.write(np.where(valid_land, SDR, -9999.0).astype(np.float32), 1)
        
    exp_path = r'Outputs\Rasters\sediment_export_t_ha_yr.tif'
    with rasterio.open(exp_path, 'w', **meta) as dst:
        dst.write(np.where(valid_land, sediment_export, -9999.0).astype(np.float32), 1)
        
    meta_u8 = meta.copy()
    meta_u8.update(dtype='uint8', nodata=0)
    stream_path = r'Outputs\Rasters\stream_network_30m.tif'
    with rasterio.open(stream_path, 'w', **meta_u8) as dst:
        dst.write(smooth_drainage.astype(np.uint8), 1)
        
    print(f"SDR rasters saved successfully:")
    print(f"  Mean IC: {np.nanmean(IC):.2f}")
    print(f"  Mean SDR: {np.mean(SDR[valid_land]):.3f} (Delivery fraction)")
    print(f"  Mean Net Sediment Export: {np.mean(sediment_export[valid_land]):.2f} t/ha/yr")
    
    # 9. Extract Lake Tana Major Tributary Sub-Catchments & Siltation Metrics
    print("Delineating Lake Tana tributary sub-basins...")
    # Bounding definitions for the 4 primary tributary catchments draining into Lake Tana (UTM 37N coordinates):
    # - Gilgel Abay: South / SW catchment (Easting 280-360 km, Northing 1215-1320 km)
    # - Gumara: South-East catchment (Easting 350-420 km, Northing 1280-1340 km)
    # - Ribb: East catchment (Easting 350-430 km, Northing 1330-1380 km)
    # - Megech: North catchment (Easting 320-380 km, Northing 1360-1425 km)
    # - Direct Littoral: Shoreline perimeter
    
    ny, nx = dem.shape
    x_coords = np.linspace(bounds.left, bounds.right, nx)
    y_coords = np.linspace(bounds.top, bounds.bottom, ny)
    xx, yy = np.meshgrid(x_coords, y_coords)
    
    subcatchments = [
        ('Gilgel Abay (Little Blue Nile)', (xx >= 281455) & (xx <= 360000) & (yy >= 1215943) & (yy <= 1315000) & valid_land),
        ('Gumara River (Mt. Guna / SE)', (xx > 350000) & (xx <= 415000) & (yy >= 1280000) & (yy <= 1335000) & valid_land),
        ('Ribb River (Debub Gondar / East)', (xx > 350000) & (xx <= 425000) & (yy > 1335000) & (yy <= 1380000) & valid_land),
        ('Megech River (North / Gondar)', (xx >= 320000) & (xx <= 380000) & (yy > 1375000) & (yy <= 1425000) & valid_land),
        ('Direct Lake Tana Littoral Catchments', (dist_to_sink_m <= 15000) & valid_land & (~((xx > 350000) & (yy > 1335000))))
    ]
    
    sub_stats = []
    tot_lake_silt_mt = 0.0
    for name, mask_sub in subcatchments:
        sub_sl = soil_loss[mask_sub]
        sub_sdr = SDR[mask_sub]
        sub_exp = sediment_export[mask_sub]
        
        if len(sub_sl) > 1000:
            area_ha = len(sub_sl) * pixel_area_ha
            gross_erosion_mt = float(np.sum(sub_sl) * pixel_area_ha) / 1e6
            net_export_mt = float(np.sum(sub_exp) * pixel_area_ha) / 1e6
            mean_sdr_val = float(np.mean(sub_sdr))
            mean_export_rate = float(np.mean(sub_exp))
            
            # Siltation volume (assuming dry bulk density 1.25 t/m3)
            silt_vol_m3 = (net_export_mt * 1e6) / 1.25
            tot_lake_silt_mt += net_export_mt
            
            sub_stats.append({
                'Tributary_Subcatchment': name,
                'Drainage_Area_ha': area_ha,
                'Drainage_Area_sqkm': area_ha / 100.0,
                'Gross_Soil_Loss_Mt_yr': gross_erosion_mt,
                'Mean_SDR': mean_sdr_val,
                'Net_Sediment_Export_Mt_yr': net_export_mt,
                'Mean_Sediment_Export_t_ha_yr': mean_export_rate,
                'Annual_Siltation_Volume_M_m3_yr': silt_vol_m3 / 1e6
            })
            
    df_sub = pd.DataFrame(sub_stats)
    df_sub['Lake_Inflow_Contribution_Pct'] = (df_sub['Net_Sediment_Export_Mt_yr'] / tot_lake_silt_mt) * 100.0
    df_sub.to_csv(r'Outputs\Tables\Lake_Tana_Subcatchment_Sediment_Inflow.csv', index=False)
    print(f"Saved Lake Tana Sub-Basin Sediment Inflow Table: Total Inflow = {tot_lake_silt_mt:.2f} Mt/yr")
    
    return ic_path, sdr_path, exp_path, stream_path, df_sub

if __name__ == '__main__':
    ensure_dirs()
    compute_sdr_model()
