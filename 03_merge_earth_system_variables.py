## 03_merge_earth_system_variables

'''35 Earth system variables (attributes and time series) are added to the attributes or time series table.
First, 17 attributes are consecutively added to the groundwater data. Afterwards, 18 time series variables
are merged to the groundwater time series within a parallelized process (server with multiple cores needed).
The time series table is split into 100 parts for which the variables are added in parallel. In the end,
all parts are put together again.'''

# Configuration: Path names, output names and other settings are defined here.
config = {
    "basepath" : "/mnt/storage/grow/", # GROW project directory
    "wells" : "Groundwater/02_Attributes/wells_attributes_V07.txt", # path to groundwater attributes table
    "timeseries": "Groundwater/02_Timeseries/wells_timeseries_final_V07.txt", # path to groundwater time series table
    # paths to original data of Earth system variables; path to file for static variables; path to folder containing files for time series variables
    "factors": {"dem": "Topography/MERIT_DEM/MERIT/MERIT_DEM.tif",
                "slope": "Topography/Slope_MERIT_DEM/dtm_slope_merit.dem_m_250m_s0..0cm_2018_v1.0.tif",
                "glim":{"data":'Soils_Geology/GLiM/glim_wgs84_0point5deg.txt.asc',"codes":"Soils_Geology/GLiM/Classnames.txt"},
                "permeability":"Soils_Geology/GLYHMPS2.0/GLHYMPS.shp",
                "porosity":"Soils_Geology/GLHYMPS/GLHYMPS/GLHYMPS.gdb/a00000009.gdbtable",
                "glacier_permafrost":"HydroAtlas/BasinATLAS_v10_shp/BasinATLAS_v10_lev12.shp",
                "whymap": "Soils_Geology/WHYMAP/WHYMAP_WOKAM/shp/whymap_karst__v1_poly.shp",
                "ggde": {"data":"Vegetation/GGDE/Huggins/gde-map.tif","codes":"Vegetation/GGDE/Huggins/GGDE_names.txt"},
                "basins": "HydroAtlas/BasinATLAS_v10_shp/BasinATLAS_v10_lev09.shp",
                "rivers": "HydroAtlas/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp",
                "drain_den": "HydroAtlas/Drainage_density.shp",
                "mswep": "Climate/Precipitation/MSWEP/Selection",
                "gpcc": "Climate/Precipitation/GPCC/",
                "gleam": "Climate/GLEAM/daily_4_1a",
                "hydrobelts": {"data":"Climate/Hydroregions/meybeck_et_al_2013_hydrobelts_shp/meybeck_et_al_2013_hydrobelts.shp", "codes":"Climate/Hydroregions/hydrobelts_codes_reduced.txt"},
                "gk": {"data": "Climate/Climate_zones/CHELSA_kg0_1981-2010_V.2.1.tif", "codes": "Climate/Climate_zones/kg0_names.txt"},
                "abstract_ind": "Humankind/Abstraction/ISIMIP_industrial",
                "abstract_dom": "Humankind/Abstraction/ISIMIP_domestic",
                "ndvi": "Vegetation/NDVI/access",
                "pet": "Climate/ERA5-Land_daily/PET",
                "temperature": "Climate/ERA5-Land_daily/Temperature",
                "snow_depth": "Climate/ERA5-Land_daily/Snow",
                "LAI_low" :"Climate/ERA5-Land_daily/LAI_low",
                "LAI_high" :"Climate/ERA5-Land_daily/LAI_high",
                "dist_streams": "Surface_Waters/Distance_perennial_streams/L01_m.tiff",
                "gw_scapes": "Humankind/Groundwaterscapes/groundwaterscapes.tif",
                "lu_totals": "Humankind/Land_use/totals",
                "lu_urban": "Humankind/Land_use/urban",
                "soil_texture": {"data":"Soils_Geology/HiHydroKlass/STC","codes":"Soils_Geology/HiHydroKlass/STC/stc_codes.txt"},
                "soil_kat": "Soils_Geology/HiHydroKlass/Ksat"},
    "output": "Groundwater/GROW_merge_V07/", # folder in which the interim output files are exported
    # names of final GROW tables
    "output_tables": {"att":"final_grow/V07/grow_attributes_without_lu.csv","att_full":"final_grow/V07/grow_attributes.csv",
                      "att_geo":"final_grow/V07/grow_attributes.json", "ts": "final_grow/V07/grow_timeseries.csv"},
    # With modules, the processing of all static variables, all time series variables or single variables can be enabled (True) or disabled (False)
    "modules": {"gk":True,"hydrobelts":True, "dem":True,"slope":True,"glim":True,"aquifer":True,
                "permeability":True,"porosity":True,"soil_text":True,"soil_kat":True,
                "dist_streams":True,"dd":True,"dd_calc_map":False,"cyro":True,
                "ggde":True,"gw_scapes":True,"mswep":False,"ndvi":True,"gleam":False,
                "era5":False,"gpcc":True,"abstract":True,"lu":True,
                "static":False, "timeseries":True, "joinall":True},
    "cores_num":100 # Number of cores that work in parallel in the time series module
}

# Import packages
'''The package pyarrow needs to be installed as well'''
import os # built-in package
import pandas as pd # imported version: 2.2.3
import geopandas as gpd # imported version: 1.0.1
import multiprocessing # built-in package
import warnings # built-in package
from func_merge_earth_system_variables import merge_vector_point # merge static vector data to groundwater attributes table
from func_merge_earth_system_variables import merge_raster_static # merge static raster data to groundwater attributes table
from func_merge_earth_system_variables import merge_raster_transient # extracts raster values at well locations
from func_merge_earth_system_variables import agg_n_mer # aggregation of time series variables and merge to groundwater time series
from func_merge_earth_system_variables import get_paths # derived file paths
from func_merge_earth_system_variables import calc_dd # calculates global drainage density map

warnings.filterwarnings("ignore")

## Earth system attributes: static

wells = pd.read_csv(config["basepath"] + config["wells"], sep=";") # import groundwater attributes table

if config["modules"]["static"]:

    '''Consecutively, the 17 static variables are merged to the groundwater table.'''

    # Koeppen-Geiger-classification [classes]
    if config["modules"]["gk"]:
        # merge_raster_static: variable (raster data) is merged to groundwater attributes table
        merge_raster_static(wells, config["basepath"] + config["factors"]["gk"]["data"], col_name="koeppen_geiger")
        # convert numeric value to class name
        codes = pd.read_csv(config["basepath"] + config["factors"]["gk"]["codes"], sep=";")
        wells = pd.merge(wells, codes, how="left", left_on="koeppen_geiger", right_on="value").drop(columns={"koeppen_geiger","value","group","class name"})

    # Hydrobelts [classes]
    if config["modules"]["hydrobelts"]:
        hydrobelt = gpd.read_file(config["basepath"] + config["factors"]["hydrobelts"]["data"])
        # convert numeric value to class name
        codes = pd.read_csv(config["basepath"] + config["factors"]["hydrobelts"]["codes"], sep="\t")
        hydrobelt = pd.merge(hydrobelt,codes, how="left", left_on="hydrobelt", right_on="code")
        # merge_vector_point: variable (vector data) is merged to groundwater attributes table
        wells = merge_vector_point(hydrobelt, wells, ["hydrobel"])
        wells.rename(columns={"hydrobel": "hydrobelt_class"}, inplace=True)

    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False) # save interim result

    # MERIT DEM [m]
    if config["modules"]["dem"]:
        merge_raster_static(wells,config["basepath"]+ config["factors"]["dem"], col_name="ground_elevation_m_asl")

    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False) # save interim result

    # Topographic Slope [°]
    if config["modules"]["slope"]:
        merge_raster_static(wells, config["basepath"] + config["factors"]["slope"], col_name="topographic_slope_degree")
        wells["topographic_slope_degree"] = wells["topographic_slope_degree"]/100 # rescale (scaling factor in original data)

    # GLiM - rock type class and aquifer type 1
    if config["modules"]["glim"]:
        merge_raster_static(wells,config["basepath"]+ config["factors"]["glim"]["data"], col_name="rock_type")
        # convert numeric value to class name and add aquifer type based on GLiM rock type (see names table)
        names = pd.read_csv(config["basepath"] + config["factors"]["glim"]["codes"], sep=";")[["Value_","rock_type_class", "aquifer_type_class"]]
        wells = pd.merge(wells, names, how="left", left_on="rock_type", right_on="Value_").drop(columns={"Value_","rock_type"},axis=1)

    # Aquifer type 2: Overwrite karst regions with Whymap [porous, fractured, karst or water_body]
    if config["modules"]["aquifer"]:
        # add whymap
        whymap = gpd.read_file(config["basepath"] + config["factors"]["whymap"]) # no interpolation of na because it means no karst
        wells = merge_vector_point(whymap, wells, ["rock_type"])
        # overwrite aquifer type when whymap indicates karst aquifer
        wells.aquifer_type_class[wells.rock_type.notna()] = "karst"
        wells.drop(["rock_type"],axis=1, inplace=True)

    # GLHYMPS2.0 - Permeability [m²]
    if config["modules"]["permeability"]:
        glyhmps2 = gpd.read_file(config["basepath"]+config["factors"]["permeability"])
        wells = merge_vector_point(glyhmps2,wells,["logK_Ferr_"])
        wells = wells[~wells.duplicated(subset="ID")] # some duplicates were created in merge_vector_point which are removed here
        wells.rename(columns={"logK_Ferr_":"permeability_0-100_m_m-2"}, inplace=True)
        wells["permeability_0-100_m_m-2"] = 10**(wells["permeability_0-100_m_m-2"]/100) # rescale permeability (see GLHYMPS README for more information)
        wells["permeability_0-100_m_m-2"][wells["permeability_0-100_m_m-2"] == 1] = None # outlier that needs to be removed

    # GLHYMPS - Porosity [0-1]
    if config["modules"]["porosity"]:
        glyhmps = gpd.read_file(config["basepath"]+config["factors"]["porosity"])
        wells = merge_vector_point(glyhmps,wells,["Porosity"])
        wells = wells[~wells.duplicated(subset="ID")] # some duplicates were created in merge_vector_point which are removed here
        wells.rename(columns={"Porosity":"total_porosity_0-100_m_fraction"}, inplace=True)

    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False) # save interim results

    # HiHydroSoil - Soil texture [classes]
    if config["modules"]["soil_text"]:
        soilfiles = get_paths(config["basepath"] + config["factors"]["soil_texture"]["data"], "tif") # gets path to topsoil and subsoil file
        # loop for soil texture class for topsoil and subsoil
        for file in soilfiles:
            merge_raster_static(wells, file, col_name=file[-22:-4]) # column name is derived by file name
            # Convert numeric values to class names
            codes = pd.read_csv(config["basepath"] + config["factors"]["soil_texture"]["codes"], sep=";")
            wells = pd.merge(wells, codes, how="left", left_on=file[-22:-4], right_on="value").drop(columns={file[-22:-4], "value"})
            wells.rename(columns={"soil_texture_class": file[-22:-4]}, inplace=True) # rename so that the column is not overwriten in the second run of the loop
        wells.rename(columns={"STC_M_250m_SUBSOIL": "soil_texture_30-200_cm_class",
                              "STC_M_250m_TOPSOIL": "soil_texture_0-30_cm_class"}, inplace=True)

    # HiHydroSoil - Saturated hydraulic conductivity of soil [cm/d]
    if config["modules"]["soil_kat"]:
        soilfiles = get_paths(config["basepath"] + config["factors"]["soil_kat"], "tif") # gets path to topsoil and subsoil file
        for file in soilfiles:
            merge_raster_static(wells, file, col_name=file[-22:-4])
            wells[file[-22:-4]] = wells[file[-22:-4]] * 0.0001  # rescale (scaling factor in original data)
        wells.rename(columns={"sat_M_250m_SUBSOIL": "soil_saturated_conductivity_30-200_cm_cm_d-1",
                              "sat_M_250m_TOPSOIL": "soil_saturated_conductivity_0-30_cm_cm_d-1"}, inplace=True)

    # Distance between perennial streams [m]
    if config["modules"]["dist_streams"]:
        # transform = True means that the coordinate system is reprojected to WGS 84
        merge_raster_static(wells, config["basepath"] + config["factors"]["dist_streams"], col_name="distance_perennial_streams_m", transform = True)

    # Drainage Density [1/m]
    if config["modules"]["dd"]:
        if config["modules"]["dd_calc_map"]:
            # The global drainage density map is calculated
            calc_dd(config)

        # the drainage density map is read and the information is merged to the attributes table
        drain_den = gpd.read_file(config["basepath"] + config["factors"]["drain_den"])
        wells = merge_vector_point(drain_den, wells, ["Drainage_d"])
        wells = wells[~wells.duplicated(subset="ID")] # some duplicates were created in merge_vector_point which are removed here
        wells.rename(columns={"Drainage_d": "drainage_density_m-1"},inplace=True)

    # Glacier extent [fraction] and Permafrost extent [fraction] from BasinAtlas level 12
    if config["modules"]["cyro"]:
        basins = gpd.read_file(config["basepath"] + config["factors"]["glacier_permafrost"])
        wells = merge_vector_point(basins, wells, ["gla_pc_use", "prm_pc_use"]) #
        wells = wells[~wells.duplicated(subset="ID")]
        wells["gla_pc_use"] = wells["gla_pc_use"] / 100 # from percentage to fraction
        wells["prm_pc_use"] = wells["prm_pc_use"] / 100  # from percentage to fraction
        wells.rename(columns={"gla_pc_use": "glacier_cover_fraction", "prm_pc_use": "permafrost_cover_fraction"}, inplace=True)

    # Groundwater-Dependent Ecosystems (GDE) [classes]
    if config["modules"]["ggde"]:
        merge_raster_static(wells,config["basepath"]+ config["factors"]["ggde"]["data"], col_name="groundwater_dependent_ecosystems_class")
        wells["groundwater_dependent_ecosystems_class"] = wells["groundwater_dependent_ecosystems_class"].astype("Int64") # convert to Int64 format (normal "int" gets the error: cannot convert float NaN to integer)
        #convert to numeric values to class names
        names = pd.read_csv(config["basepath"] + config["factors"]["ggde"]["codes"], sep=";")
        wells = pd.merge(wells, names, how="left", left_on="groundwater_dependent_ecosystems_class", right_on="code").drop(columns={"groundwater_dependent_ecosystems_class","code"},axis=1)
        wells.rename(columns={"name":"groundwater_dependent_ecosystems_class"},inplace=True)
        wells.groundwater_dependent_ecosystems_class[wells["groundwater_dependent_ecosystems_class"].isna()] = "No GDE" # No value means that there is no GDE

    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False) # save interim results

    # Groundwaterscapes [classes]
    if config["modules"]["gw_scapes"]:
        merge_raster_static(wells, config["basepath"] + config["factors"]["gw_scapes"],col_name="groundwaterscapes_ID_class")
        wells["groundwaterscapes_ID_class"] = wells["groundwaterscapes_ID_class"].fillna(0).astype("int") # convert to integer to reduce memory
        wells["groundwaterscapes_ID_class"][wells["groundwaterscapes_ID_class"] == 0] = None # remove 0 that are created because of the last step

    # Final export of attributes table
    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False)

## Earth system variables: Time series

if config["modules"]["timeseries"]:

    def total_merge(df,ts,i):

        '''In this function, the time series variables are:
        1) extracted at the well location,
        2) temporally aggregated to the individual groundwater time series resolutions and
        3) merged to the groundwater time series table

        df: part of attributes table (contains well location coordinates)
        ts: groundwater time series table
        i: index of well table part'''

        # trim timeseries by well IDs (one of the 100 parts)
        ts = pd.merge(ts, df["ID"], how="inner", on="ID")

        # MSWEP - 3-hourly precipitation [mm/3 hours]
        if config["modules"]["mswep"]:
            # merge_raster_transient: the raster values are extracted at the well location and exported to a file
            merge_raster_transient(df, config["basepath"] + config["factors"]["mswep"], config,bandname="precipitation", col_name="mswep_" + str(i))
            # agg_n_mer: the extracted values are temporally aggregated and merged with the trimmed time series table
            agg_n_mer(config, ts,"mswep",i, "mswep_join_" + str(i), daily=True)

        # NDVI - daily
        if config["modules"]["ndvi"]:
            merge_raster_transient(df, config["basepath"] + config["factors"]["ndvi"], config, bandname="NDVI",latname="latitude", lonname="longitude", col_name="ndvi_" + str(i), ndvi=True)
            agg_n_mer(config, ts, "ndvi",i, "ndvi_join_" + str(i), daily=True)

        # GLEAM - daily
        if config["modules"]["gleam"]:

            # Potential Evapotranspiration [mm/day], Actual Evapotranspiration [mm/day], Interception loss [mm/day]
            gleam = [["Ep","E_","Ei"],["Ep","E","Ei"]]

            for r in range(len(gleam[0])):
                merge_raster_transient(df, config["basepath"] + config["factors"]["gleam"], config, patt2=gleam[0][r], ful=True,bandname=gleam[1][r], col_name= gleam[1][r]+"_" + str(i))
                agg_n_mer(config, ts, gleam[1][r],i ,gleam[1][r]+"_join_" + str(i), daily=True)

        # ERA5 - daily
        if config["modules"]["era5"]:

            # Potential evapotranspiration [m], Air temperature in 2m [K], Snow depth [m], Leaf area index of low vegetation, Leaf area index of high vegetation
            era5_names = [["pev","t2m","sde","lai_lv","lai_hv"],["pet","temperature","snow_depth","LAI_low","LAI_high"]]

            for r in range(len(era5_names[0])):
                merge_raster_transient(df, config["basepath"] + config["factors"][era5_names[1][r]], config, bandname=era5_names[0][r],latname="latitude", lonname="longitude", col_name=era5_names[1][r] + "_" + str(i), era5 = True)
                agg_n_mer(config, ts, era5_names[1][r], i, era5_names[1][r]+"_join_"+ str(i), daily=True)


        # GPCC - monthly precipitation [mm/month]
        if config["modules"]["gpcc"]:
            merge_raster_transient(df, config["basepath"] + config["factors"]["gpcc"], config, bandname="precip",col_name="gpcc_"+ str(i))
            agg_n_mer(config, ts, "gpcc",i , "gpcc_join_"+ str(i))

        # Water withdrawal [m³/year] - yearly
        if config["modules"]["abstract"]:
            # industrial use
            merge_raster_transient(df, config["basepath"] + config["factors"]["abstract_ind"], config,bandname="indww", col_name="withdrawal_industrial_"+ str(i),isimip=True)
            agg_n_mer(config, ts, "withdrawal_industrial", i,"withdrawal_industrial_join_"+ str(i), yearly=True)

            # domestic use
            merge_raster_transient(df, config["basepath"] + config["factors"]["abstract_dom"], config, bandname="domww",col_name="withdrawal_domestic_" + str(i), isimip=True)
            agg_n_mer(config, ts, "withdrawal_domestic", i, "withdrawal_domestic_join_" + str(i), yearly=True)

        # Land use fraction - yearly
        if config["modules"]["lu"]:

            lu_names = ["cropland_irrigated","cropland_rainfed","pastures","forests_and_natural_vegetation"]

            for name in lu_names:
                merge_raster_transient(df, config["basepath"] + config["factors"]["lu_totals"], config, bandname=name,col_name=name + "_" + str(i), isimip=True)
                agg_n_mer(config, ts, name, i, name + "_join_" + str(i), yearly=True)

            # urban area fraction (is in a different nc than the rest)
            merge_raster_transient(df, config["basepath"] + config["factors"]["lu_urban"], config,bandname="urbanareas", col_name="urbanareas_" + str(i), isimip=True)
            agg_n_mer(config, ts, "urbanareas", i, "urbanareas_join_" + str(i), yearly=True)

        # Join all variables
        # if there already exits a file where all variables are merged in one table ("join_all..."), it is removed
        # Otherwise, this file would be selected in the next step as well
        if os.path.isfile(config["basepath"] + config["output"] + "join_all_" + str(i) +".txt"):
            os.remove(config["basepath"] + config["output"] + "join_all_" + str(i) +".txt")

        # get file paths of all aggregated and merged variable-groundwater tables
        join_files = get_paths(config["basepath"] + config["output"], "join", str(i), True)# to save a merge here

        # Merge each variable as column to the groundwater table
        for ele in join_files:
            par = pd.read_csv(ele, sep=";")
            # only the columns of "ID", "date" and the variable shall be left, the rest can be removed
            par.drop(par.columns[5:15], axis=1, inplace=True)
            par.drop(par.columns[1:4], axis=1, inplace=True)
            par.date = pd.to_datetime(par.date, format='ISO8601') # convert date so that merge works
            ts = pd.merge(ts,par, how="left", on=["ID","date"]) # merge via ID and date of time series

        # export time series table
        ts.to_csv(config["basepath"] + config["output"] + "join_all_" + str(i) +".txt", sep=";", index=False)

    # Load full groundwater time series table
    ts = pd.read_csv(config["basepath"] + config["timeseries"], sep=";")
    ts.date = pd.to_datetime(ts.date, format='ISO8601') # convert date column
    ts.date2 = pd.to_datetime(ts.date2, format='ISO8601') # convert second date column

    # Start multiprocessing (Parallel workflow)
    if __name__ == '__main__':
        # Determine the number of available cores
        num_cores = multiprocessing.cpu_count()
        def run_function(func_arg):
            func, arg = func_arg
            return func(*arg)

        # Create a list of functions and their arguments
        functions = []

        # split attributes table in 100 parts
        n = len(wells)
        wells_parts = [wells.iloc[i * n // config["cores_num"]:(i + 1) * n // config["cores_num"]] for i in range(config["cores_num"])]

        # create a roadmap where for each of the 100 processes the function that shall run and the inserted parameters are given
        for i,ele in enumerate(wells_parts):
            r = str(i) # to name the expoted files with the index
            # for correct order of the exported files, the index is manually put to 01 to 09
            if i < 10:
                r = "0"+ r
            # Add function and parameter to list
            functions.append((total_merge, (ele,ts,r)))

        # Start the 100 processes
        with multiprocessing.Pool(processes=num_cores) as pool:
            results = pool.map(run_function, functions)

## Join all 100 parts to one final time series table
## Main land use is added to attributes table

if config["modules"]["joinall"]:

    # Get file paths of all 100 groundwater-variables time series parts
    files = get_paths(config["basepath"]+config["output"],"join_all")

    # Append them to one table
    parts = []
    for file in files:
        parts.append(pd.read_csv(file, sep=";"))
    ts_fin = pd.concat(parts).reset_index(drop=True)

    # Convert Units
    ts_fin["temperature"] = ts_fin["temperature"] - 273 # from K to C°
    # to mm/year
    ts_fin.gpcc = ts_fin.gpcc * 12 # from mm/month to mm/year
    ts_fin.mswep = ts_fin.mswep * 2920 # from mm/3hours to mm/year
    ts_fin.Ep = ts_fin.Ep * 365 # from mm/day to mm/year
    ts_fin.E = ts_fin.E * 365 # from mm/day to mm/year
    ts_fin.Ei = ts_fin.Ei * 365 # from mm/day to mm/year
    ts_fin.pet = ts_fin.pet * 365 * 1000 * -1 # from m/day to mm/year; in ERA5-Land PET is negative as upward flux --> change the sign

    # Reorder and rename columns
    ts_fin.drop(columns=["date2","country"], inplace = True) # delete second date column as it is not needed anymore
    # reorder columns
    ts_fin = pd.concat([ts_fin.iloc[:,:13],ts_fin.mswep,ts_fin.gpcc,ts_fin.pet,ts_fin.Ep,ts_fin.E,ts_fin.Ei,
                        ts_fin.temperature,ts_fin.snow_depth,ts_fin.ndvi,ts_fin.LAI_low,ts_fin.LAI_high,
                        ts_fin.withdrawal_industrial,ts_fin.withdrawal_domestic,ts_fin.urbanareas,ts_fin.pastures,
                        ts_fin.cropland_rainfed,ts_fin.cropland_irrigated,ts_fin.forests_and_natural_vegetation], axis=1)
    # rename columns
    ts_fin.rename(columns={"ID":"GROW_ID",
                           "mswep":"precipitation_mswep_mm_year-1",
                           "gpcc":"precipitation_gpcc_mm_year-1",
                           "pet":"potential_evapotranspiration_era5_mm_year-1",
                           "Ep":"potential_evapotranspiration_gleam_mm_year-1",
                           "E":"actual_evapotranspiration_mm_year-1",
                           "temperature":"air_temperature_C°",
                           "snow_depth":"snow_depth_m",
                           "Ei":"interception_mm_year-1",
                           "ndvi":"ndvi_ratio",
                           "LAI_low":"lai_low_vegetation_ratio",
                           "LAI_high":"lai_high_vegetation_ratio",
                           "withdrawal_industrial":"withdrawal_industrial_m3_year-1",
                           "withdrawal_domestic":"withdrawal_domestic_m3_year-1",
                           "urbanareas":"urban_area_fraction",
                           "pastures":"pastures_fraction",
                           "cropland_rainfed":"cropland_rainfed_fraction",
                           "cropland_irrigated":"cropland_irrigated_fraction",
                           "forests_and_natural_vegetation":"forests_natural_vegetation_fraction"},inplace=True)

    # final export of full time series table
    ts_fin.to_csv(config["basepath"]+config["output_tables"]["ts"], sep=";", index=False)
    ts_fin.to_parquet(config["basepath"] + config["output_tables"]["ts"], index=False)

    # add main land use to attributes table
    wells = pd.read_csv(config["basepath"] + config["output_tables"]["att"], sep=";")  # import attributes table
    wells.rename(columns={"ID":"GROW_ID","name ":"name"},inplace=True) # rename ID column and name column

    ts_fin_reduced = ts_fin[["GROW_ID","forests_natural_vegetation_fraction","urban_area_fraction","cropland_rainfed_fraction","cropland_irrigated_fraction","pastures_fraction"]] # reduce to land use columns
    ts_agg = ts_fin_reduced.groupby(by="GROW_ID", as_index=False).mean() # aggregate to average land use fraction

    wells_fin = pd.merge(wells, ts_agg, on="GROW_ID", how="outer") # merge mean land use fraction to well attributes
    # find land use with highest fraction per row and set it as main land use
    wells_fin['main_landuse_class'] = wells_fin[["forests_natural_vegetation_fraction","urban_area_fraction","cropland_rainfed_fraction","cropland_irrigated_fraction","pastures_fraction"]].idxmax(axis=1)
    # if every land use fraction is 0, set main land use to NA
    wells_fin.main_landuse_class[(wells_fin[["forests_natural_vegetation_fraction","urban_area_fraction","cropland_rainfed_fraction","cropland_irrigated_fraction","pastures_fraction"]] == 0).all(axis=1)] = None
    # drop single land use fractions
    wells_fin.drop(columns=["forests_natural_vegetation_fraction","urban_area_fraction","cropland_rainfed_fraction","cropland_irrigated_fraction","pastures_fraction"],inplace=True)

    ## export attributes as csv, parquet and json
    # Convert data type to only string so that export as parquet works
    wells_fin['original_ID_groundwater'] = wells_fin['original_ID_groundwater'].astype("str")
    wells_fin['name'] = wells_fin['name'].astype("str")
    wells_fin.to_csv(config["basepath"]+config["output_tables"]["att_full"], sep=";", index=False)
    wells_fin.to_parquet(config["basepath"] + config["output_tables"]["att_full"], index=False, engine='pyarrow')

    # export geo-referenced attributes
    wells_gdf = gpd.GeoDataFrame(wells_fin, geometry=gpd.points_from_xy(wells_fin.longitude, wells_fin.latitude),crs="EPSG:4326")
    wells_gdf.to_file(config["basepath"] + config["output_tables"]["att_geo"], index=False)
