### Merge all static and transient variables to groundwater data

"""In this script..."""


import pickle
import os # does not need to be installed
import pandas as pd
import geopandas as gpd
import multiprocessing # does not need to be installed
import warnings # does not need to be installed
from func_merge_earth_system_variables import merge_vector_point
from func_merge_earth_system_variables import merge_raster_static
from func_merge_earth_system_variables import merge_raster_transient
from func_merge_earth_system_variables import agg_n_mer
from func_merge_earth_system_variables import get_paths
from func_merge_earth_system_variables import calc_dd

warnings.filterwarnings("ignore")

# Configuration
with open('config.pkl', 'rb') as f:
    all_configs = pickle.load(f)

config = all_configs["config_04"]

## Earth system attributes: static

wells = pd.read_csv(config["basepath"] + config["wells"], sep=";")

if config["modules"]["static"]:

    # Koeppen-Geiger-classification [classes]
    if config["modules"]["gk"]:
        merge_raster_static(wells, config["basepath"] + config["factors"]["gk"]["data"], col_name="koeppen_geiger", nodat=0)
        # convert numeric value to class code
        codes = pd.read_csv(config["basepath"] + config["factors"]["gk"]["codes"], sep=";")
        wells = pd.merge(wells, codes, how="left", left_on="koeppen_geiger", right_on="value").drop(columns={"koeppen_geiger","value","class name"})

    # Hydroregions [classes]
    if config["modules"]["hydrobelts"]:
        hydrobelt = gpd.read_file(config["basepath"] + config["factors"]["hydrobelts"]["data"])
        codes = pd.read_csv(config["basepath"] + config["factors"]["hydrobelts"]["codes"], sep="\t")
        hydrobelt = pd.merge(hydrobelt,codes, how="left", left_on="hydrobelt", right_on="code")
        wells = merge_vector_point(hydrobelt, wells, ["hydrobel"])
        wells.rename(columns={"hydrobel": "hydrobelt_class"}, inplace=True)

    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False)

    # MERIT DEM [m]
    if config["modules"]["dem"]:
        merge_raster_static(wells,config["basepath"]+ config["factors"]["dem"], col_name="ground_elevation_m_asl")

    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False)

    # Topographic Slope [°]
    if config["modules"]["slope"]:
        merge_raster_static(wells, config["basepath"] + config["factors"]["slope"], col_name="topographic_slope_degree")
        wells["topographic_slope_degree"] = wells["topographic_slope_degree"]/100 # scaling factor of 100

    # GLiM - rock type class
    if config["modules"]["glim"]:
        merge_raster_static(wells,config["basepath"]+ config["factors"]["glim"]["data"], col_name="rock_type", glim=True)
        # convert numeric value to class name
        names = pd.read_csv(config["basepath"] + config["factors"]["glim"]["codes"], sep=";")[["Value_","rock_type_class", "aquifer_type_class"]]
        wells = pd.merge(wells, names, how="left", left_on="rock_type", right_on="Value_").drop(columns={"Value_","rock_type"},axis=1)

    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False)

    # GLYHMPS - Permeability [m²], Porosity [0-1] and Permafrost
    if config["modules"]["glymps"]:
        glyhmps = gpd.read_file(config["basepath"]+config["factors"]["glymps"])
        wells = merge_vector_point(glyhmps,wells,["logK_Ferr_","Porosity_x"])
        wells = wells[~wells.duplicated(subset="ID")]
        wells.rename(columns={"logK_Ferr_":"permeability_m2","Porosity_x":"porosity_fraction"}, inplace=True)
        wells.permeability_m2 = 10**(wells.permeability_m2/100) # rescale permeability (GLYHMPS readme)
        wells.porosity_fraction = wells.porosity_fraction/100 # rescale porosity

    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False)

    # Estimate aquifer type with GLim and Whymap [porous, fractured, karst or water_body]
    if config["modules"]["aquifer"]:
        # add whymap
        whymap = gpd.read_file(config["basepath"] + config["factors"]["whymap"]) # no interpolation of na because it means no karst
        wells = merge_vector_point(whymap, wells, ["rock_type"])
        # overwrite aquifer type when whymap indicates karst aquifer
        wells.aquifer_type_class[wells.rock_type.notna()] = "karst"
        wells.drop(["rock_type"],axis=1, inplace=True)

    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False)

    # HiHydroSoil - Soil texture [classes]
    if config["modules"]["soil_text"]:
        soilfiles = get_paths(config["basepath"] + config["factors"]["soil_texture"]["data"], "tif")
        for file in soilfiles:
            merge_raster_static(wells, file, col_name=file[-22:-4])
            codes = pd.read_csv(config["basepath"] + config["factors"]["soil_texture"]["codes"], sep=";")
            wells = pd.merge(wells, codes, how="left", left_on=file[-22:-4], right_on="value").drop(
                columns={file[-22:-4], "value"})
            wells.rename(columns={"soil_texture_class": file[-22:-4]}, inplace=True)
        wells.rename(columns={"STC_M_250m_SUBSOIL": "soil_texture_30-200_cm_class",
                              "STC_M_250m_TOPSOIL": "soil_texture_0-30_cm_class"}, inplace=True)

    # HiHydroSoil - Saturated hydraulic conductivity of soil [cm/d]
    if config["modules"]["soil_kat"]:
        soilfiles = get_paths(config["basepath"] + config["factors"]["soil_kat"], "tif")
        for file in soilfiles:
            merge_raster_static(wells, file, col_name=file[-22:-4])
            wells[file[-22:-4]] = wells[file[-22:-4]] * 0.0001  # scaling factor
        wells.rename(columns={"sat_M_250m_SUBSOIL": "soil_saturated_conductivity_30-200_cm_cm_d-1",
                              "sat_M_250m_TOPSOIL": "soil_saturated_conductivity_0-30_cm_cm_d-1"}, inplace=True)

    # Distance between perennial streams [m]
    if config["modules"]["dist_streams"]:
        merge_raster_static(wells, config["basepath"] + config["factors"]["dist_streams"], col_name="distance_perennial_streams_m", nodat=-3.4028230607370965e+38, transform = True)

    # Drainage Density [1/m]
    if config["modules"]["dd"]:
        if config["modules"]["dd_calc_map"]:
            calc_dd(config)

        drain_den = gpd.read_file(config["basepath"] + config["factors"]["drain_den"])
        wells = merge_vector_point(drain_den, wells, ["Drainage_d"])
        wells = wells[~wells.duplicated(subset="ID")]
        wells.rename(columns={"Drainage_d": "drainage_density_m-1"},inplace=True)

    # Global Groundwater Dependent Ecosystems [classes]
    if config["modules"]["ggde"]:
        merge_raster_static(wells,config["basepath"]+ config["factors"]["ggde"]["data"], col_name="groundwater_dependent_ecosystems_class") # no interpolation because NA means no GGDE
        wells["groundwater_dependent_ecosystems_class"] = wells["groundwater_dependent_ecosystems_class"].astype("Int64") # normal "int" gets the error: cannot convert float NaN to integer
        #convert to class names
        names = pd.read_csv(config["basepath"] + config["factors"]["ggde"]["codes"], sep=";")
        wells = pd.merge(wells, names, how="left", left_on="groundwater_dependent_ecosystems_class", right_on="code").drop(columns={"groundwater_dependent_ecosystems_class","code"},axis=1)
        wells.rename(columns={"name":"groundwater_dependent_ecosystems_class"},inplace=True)
        wells.groundwater_dependent_ecosystems_class[wells["groundwater_dependent_ecosystems_class"].isna()] = "No GDE"

    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False)

    # Groundwaterscapes [classes]
    if config["modules"]["gw_scapes"]:
        merge_raster_static(wells, config["basepath"] + config["factors"]["gw_scapes"],col_name="groundwaterscapes_ID_class", nodat=0)
        wells["groundwaterscapes_ID_class"] = wells["groundwaterscapes_ID_class"].astype("int")

    # Final export
    wells.to_csv(config["basepath"]+config["output_tables"]["att"], sep=";", index=False)

## Earth system variables: Time series

if config["modules"]["timeseries"]:

    def total_merge(df,ts,i):

        # trim timeseries by well IDs
        ts = pd.merge(ts, df["ID"], how="inner", on="ID")

        # MSWEP - 3-hourly precipitation [mm/3 hours]
        if config["modules"]["mswep"]:
            merge_raster_transient(df, config["basepath"] + config["factors"]["mswep"], config,bandname="precipitation", col_name="mswep_" + str(i), nodat=-9999)
            agg_n_mer(config, ts,"mswep",i, "mswep_join_" + str(i), daily=True)

        # NDVI - daily
        if config["modules"]["ndvi"]:
            merge_raster_transient(df, config["basepath"] + config["factors"]["ndvi"], config, bandname="NDVI",latname="latitude", lonname="longitude", col_name="ndvi_" + str(i), ndvi=True, nodat=-9999)
            agg_n_mer(config, ts, "ndvi",i, "ndvi_join_" + str(i), daily=True)

        # GLEAM - daily
        if config["modules"]["gleam"]:

            # Potential Evapotranspiration [mm/day], Actual Evapotranspiration [mm/day], Interception loss [mm/day]
            gleam = [["Ep","E_","Ei"],["Ep","E","Ei"]]

            for r in range(len(gleam[0])):
                merge_raster_transient(df, config["basepath"] + config["factors"]["gleam"], config, patt2=gleam[0][r], ful=True,bandname=gleam[1][r], col_name= gleam[1][r]+"_" + str(i),nodat=-999)
                agg_n_mer(config, ts, gleam[1][r],i ,gleam[1][r]+"_join_" + str(i), daily=True)

        # ERA5 - monthly
        if config["modules"]["era5"]:

            # Air temperature in 2m [K], Snow depth [m], Leaf area index of low vegetation, Leaf area index of high vegetation
            era5_names = [["t2m","sde","lai_lv","lai_hv"],["temperature","snow_depth","LAI_low","LAI_high"]]

            for r in range(len(era5_names[0])):
                merge_raster_transient(df, config["basepath"] + config["factors"][era5_names[1][r]], config, bandname=era5_names[0][r],latname="latitude", lonname="longitude", col_name=era5_names[1][r] + "_" + str(i), era5 = True)
                agg_n_mer(config, ts, era5_names[1][r], i, era5_names[1][r]+"_join_"+ str(i), daily=True)


        # GPCC - monthly precipitation [mm/month]
        if config["modules"]["gpcc"]:
            merge_raster_transient(df, config["basepath"] + config["factors"]["gpcc"], config, bandname="precip",col_name="gpcc_"+ str(i),nodat=-99999.9921875)
            agg_n_mer(config, ts, "gpcc",i , "gpcc_join_"+ str(i))

        # Water withdrawal [m³/year]
        if config["modules"]["abstract"]:
            # industrial use
            merge_raster_transient(df, config["basepath"] + config["factors"]["abstract_ind"], config,bandname="indww", col_name="withdrawal_industrial_"+ str(i),isimip=True, nodat=1.0000000200408773e+20)
            agg_n_mer(config, ts, "withdrawal_industrial", i,"withdrawal_industrial_join_"+ str(i), yearly=True)

            # domestic use
            merge_raster_transient(df, config["basepath"] + config["factors"]["abstract_dom"], config, bandname="domww",col_name="withdrawal_domestic_" + str(i), isimip=True, nodat=1.0000000200408773e+20)
            agg_n_mer(config, ts, "withdrawal_domestic", i, "withdrawal_domestic_join_" + str(i), yearly=True)

        # Land use fraction
        if config["modules"]["lu"]:

            lu_names = ["cropland_irrigated","cropland_rainfed","pastures","forests_and_natural_vegetation"]

            for name in lu_names:
                merge_raster_transient(df, config["basepath"] + config["factors"]["lu_totals"], config, bandname=name,col_name=name + "_" + str(i), isimip=True, nodat=1.0000000200408773e+20)
                agg_n_mer(config, ts, name, i, name + "_join_" + str(i), yearly=True)

            # urban area fraction (is in a different nc than the rest)
            merge_raster_transient(df, config["basepath"] + config["factors"]["lu_urban"], config,bandname="urbanareas", col_name="urbanareas_" + str(i), isimip=True, nodat=1.0000000200408773e+20)
            agg_n_mer(config, ts, "urbanareas", i, "urbanareas_join_" + str(i), yearly=True)

        # Join all variables
        if os.path.isfile(config["basepath"] + config["output"] + "join_all_" + str(i) +".txt"):
            os.remove(config["basepath"] + config["output"] + "join_all_" + str(i) +".txt")

        join_files = get_paths(config["basepath"] + config["output"], "join", str(i), True)# to save a merge here

        for ele in join_files:
            par = pd.read_csv(ele, sep=";")
            par.drop(par.columns[5:13], axis=1, inplace=True)
            par.drop(par.columns[1:4], axis=1, inplace=True)
            par.date = pd.to_datetime(par.date, format='ISO8601')
            ts = pd.merge(ts,par, how="left", on=["ID","date"])

        # export
        ts.to_csv(config["basepath"] + config["output"] + "join_all_" + str(i) +".txt", sep=";", index=False)

    # well timeseries
    ts = pd.read_csv(config["basepath"] + config["timeseries"], sep=";")
    ts.date = pd.to_datetime(ts.date, format='ISO8601')
    ts.date2 = pd.to_datetime(ts.date2, format='ISO8601')

    #total_merge(wells, ts, "01")

    if __name__ == '__main__':
        # Determine the number of available cores
        num_cores = multiprocessing.cpu_count()
        def run_function(func_arg):
            func, arg = func_arg
            return func(*arg)

        # Create a list of functions and their arguments
        functions = []

        # split wells in 100 parts because they are 100 cores
        n = len(wells)
        wells_parts = [wells.iloc[i * n // config["cores_num"]:(i + 1) * n // config["cores_num"]] for i in range(config["cores_num"])]

        total_merge(wells_parts[5], ts, "05")

        for i,ele in enumerate(wells_parts):
            r = str(i)
            if i < 10:
                r = "0"+ r
            functions.append((total_merge, (ele,ts,r)))

        # Create a pool of worker processes
        with multiprocessing.Pool(processes=num_cores) as pool:
            # Map the functions to the pool
            results = pool.map(run_function, functions)

if config["modules"]["joinall"]:
    # Join all
    files = get_paths(config["basepath"]+config["output"],"join_all")

    parts = []
    for file in files:
        parts.append(pd.read_csv(file, sep=";"))

    ts_fin = pd.concat(parts).reset_index(drop=True)

    # Units
    ts_fin["temperature"] = ts_fin["temperature"] - 273
    # to mm/year
    ts_fin.gpcc = ts_fin.gpcc * 12
    ts_fin.mswep = ts_fin.mswep * 2920
    ts_fin.Ep = ts_fin.Ep * 365
    ts_fin.E = ts_fin.E * 365
    ts_fin.Ei = ts_fin.Ei * 365

    # Rename columns

    # Export
    ts_fin.to_csv(config["basepath"]+config["output_tables"]["ts"], sep=";", index=False)