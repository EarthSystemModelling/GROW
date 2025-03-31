"""
# TODO: docstring überarbeiten
In this script all configs needed for the processing of GROW are summarized.
After, the scripts are initialized one after another.
"""

import pickle
import subprocess


# Processing of IGRACs Groundwater Time series data (02_processing_gw_time_series)
pass

# Processing of IGRACs Groundwater attributes (04_processing_gw_attributes)
pass

# Merge earth system time series and attributes to groundwater data
pass


all_configs = {
    "config_02": {
        "basepath": "/mnt/storage/grow/Groundwater/",
        "wells": "Well_And_Monitoring_Data",
        "country_name_pos": 55,  # position of first country letter in file path to extract country information
        "output": {"name": "_Ricarda",
                   "data": "Wells_timeseries/wells_timeseries",
                   "ts_attributes": "Wells_timeseries/wells_timeseries_attributes",
                   "max_dist": "Statistics/max_distance",
                   "lost_<2_records": "Wells_timeseries/lost_<2_records",
                   "lost_mul_par": "Wells_timeseries/lost_mul_par",
                   "lost_unrealistic_value": "Wells_timeseries/lost_unrealistic_value",
                   "lost_after_aggregation": "Wells_timeseries/lost_after_aggregation",
                   "lost_gap_length": "Wells_timeseries/lost_gap_length",
                   "lost_gap_amount": "Wells_timeseries/lost_gap_amount",
                   "plateau": "Wells_timeseries/wells_plateau",
                   "jumps": "Wells_timeseries/wells_jumps",
                   "depth_all": "Wells_timeseries/wells_depth_all_negative",
                   "depth_any": "Wells_timeseries/wells_depth_some_negative",
                   "par": "Wells_timeseries/wells_mul_par_all",
                   "lost_per": "Statistics/wells_timeseries_drops",
                   "duration": "GGMN_preprocessing_duration.txt"},
        "small": True
        },

    "config_03": {
        "basepath": "/mnt/storage/grow/Groundwater/",
        "wells": "Well_And_Monitoring_Data",
        "timeseries_att": "Wells_timeseries/wells_timeseries_attributes_Ricarda.txt",
        "timeseries": "Wells_timeseries/wells_timeseries_Ricarda.txt",
        "output": {"name": "_Ricarda",
                   "all": "Wells_attributes/wells_attributes_all.csv",
                   "all_dups": "Wells_attributes/wells_attributes_all_dups.txt",
                   "filtered": "Wells_attributes/wells_attributes",
                   "points": "Other_shapes/wells_points",
                   "drops": "Statistics/wells_attributes_drops",
                   "new_ts": "Wells_timeseries/wells_timeseries_final",
                   "coor": "Wells_attributes/wells_wrong_coordinates",
                   "id": "Wells_attributes/wells_dup_id",
                   "merge": "Wells_attributes/wells_pulling_empty",
                   "duploc": "Wells_attributes/wells_duplicate_loc_time"
                   }
        },

    "config_04": {
        "basepath": "/mnt/storage/grow/",
        "wells": "Groundwater/Wells_attributes/wells_attributes_Ricarda.txt",
        "timeseries": "Groundwater/Wells_timeseries/wells_timeseries_final_Ricarda.txt",
        "factors": {"dem": "Topography/MERIT_DEM/MERIT/MERIT_DEM.tif",
                    "slope": "Topography/Slope_MERIT_DEM/dtm_slope_merit.dem_m_250m_s0..0cm_2018_v1.0.tif",
                    "glim": {"data": 'Soils_Geology/GLiM/glim_wgs84_0point5deg.txt.asc',
                             "codes": "Soils_Geology/GLiM/Classnames.txt"},
                    "glymps": "Soils_Geology/GLYHMPS/GLHYMPS.shp",
                    "whymap": "Soils_Geology/WHYMAP/WHYMAP_WOKAM/shp/whymap_karst__v1_poly.shp",
                    "ggde": {"data": "Vegetation/GGDE/Huggins/gde-map.tif",
                             "codes": "Vegetation/GGDE/Huggins/GGDE_names.txt"},
                    "basins": "HydroAtlas/BasinATLAS_v10_shp/BasinATLAS_v10_lev09.shp",
                    "rivers": "HydroAtlas/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp",
                    "drain_den": "HydroAtlas/Drainage_density.shp",
                    "mswep": "Climate/Precipitation/MSWEP/Selection",
                    "gpcc": "Climate/Precipitation/GPCC/",
                    "gleam": "Climate/GLEAM/daily_4_1a",
                    "hydrobelts": {
                        "data": "Climate/Hydroregions/meybeck_et_al_2013_hydrobelts_shp/meybeck_et_al_2013_hydrobelts.shp",
                        "codes": "Climate/Hydroregions/hydrobelts_codes_reduced.txt"},
                    "gk": {"data": "Climate/Climate_zones/CHELSA_kg0_1981-2010_V.2.1.tif",
                           "codes": "Climate/Climate_zones/kg0_names.txt"},
                    "abstract_ind": "Humankind/Abstraction/ISIMIP_industrial",
                    "abstract_dom": "Humankind/Abstraction/ISIMIP_domestic",
                    "ndvi": "Vegetation/NDVI/access",
                    "temperature": "Climate/ERA5-Land_daily/Temperature",
                    "snow_depth": "Climate/ERA5-Land_daily/Snow",
                    "LAI_low": "Climate/ERA5-Land_daily/LAI_low",
                    "LAI_high": "Climate/ERA5-Land_daily/LAI_high",
                    "dist_streams": "Surface_Waters/Distance_perennial_streams/L01_m.tiff",
                    "gw_scapes": "Humankind/Groundwaterscapes/groundwaterscapes.tif",
                    "lu_totals": "Humankind/Land_use/totals",
                    "lu_urban": "Humankind/Land_use/urban",
                    "soil_texture": {"data": "Soils_Geology/HiHydroKlass/STC",
                                     "codes": "Soils_Geology/HiHydroKlass/STC/stc_codes.txt"},
                    "soil_kat": "Soils_Geology/HiHydroKlass/Ksat"},
        "output": "Groundwater/GROW_merge_Ricarda/",
        "output_tables": {"att": "final_grow/grow_attributes_Ricarda.txt", "ts": "final_grow/grow_timeseries_Ricarda.txt"},
        "modules": {"dem": True, "slope": True, "glim": True, "soil_text": True, "soil_kat": True,
                    "dist_streams": True, "gk": True, "glymps": True, "aquifer": True, "soil": False,
                    "lu": True, "mswep": True, "gpcc": True, "ggde": True, "gw_scapes": True,
                    "gleam": True, "dd": True, "dd_calc_map": False,
                    "abstract": True, "ndvi": True, "era5": True, "hydrobelts": True,
                    "static": False, "timeseries": True, "joinall": False},
        "cores_num": 10
        # TODO: test if it still works when only one core is chosen (use on a normal computer and not a server)
        }
}

with open('config.pkl', 'wb') as f:
    pickle.dump(all_configs, f)

# Scripts nacheinander ausführen
scripts = ["02_processing_gw_time_series.py", "03_processing_gw_attributes.py", "04_merge_earth_system_variables.py"]

for script in scripts:
    result = subprocess.run(["python3", script])
    if result.returncode != 0:
        print(f"{script} failed.")
        break
