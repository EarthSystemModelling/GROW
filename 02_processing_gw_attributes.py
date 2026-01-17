"""The attributes table and the final time series table are created.
Firstly, all files containing static (no temporal dimension) well attributes per country are read and merged together.
A) duplicates by ID and country and b) duplicates by coordinates, starting date, ending date and mean groundwater
table are removed. Further, wells with coordinates outside realistic ranges are removed. Attributes generated in the
time series processing are merged to the attributes from the original data. The time series table is trimmed to
the ID's that are still left in the preprocessed attributes table. An unique GROW ID is created.
"""

import os  # built-in package
import time  # built-in package
import pandas as pd  # imported version: 2.2.3

# Configuration: Path names, output names and other settings are defined here.
config = {
    "basepath" : "/mnt/storage/grow/01_Groundwater/", # GROW project directory for groundwater data
    "wells": "01_IGRAC_data_2025_08_18", # folder in which IGRAC's groundwater data is located
    "timeseries_att": "02_Timeseries/sensitivity_analysis/gaps_10/wells_timeseries_attributes_gaps_10.txt", # time series attributes derived in "01_processing_gw_time_series"
    "timeseries" : "02_Timeseries/sensitivity_analysis/gaps_10/wells_timeseries_gaps_10.txt", # time series table derived in "01_processing_gw_time_series"
    # paths of exported output files
    "output": {"name": "_V09", # name of version
               "all":"02_Attributes/wells_attributes_all.csv",
               "all_dups": "02_Attributes/wells_attributes_all_dups.txt",
               "all_unprocessed": "02_Attributes/wells_attributes_all_unprocessed.csv",
               "filtered": "02_Attributes/wells_attributes",
               "drops": "02_Statistics/wells_attributes_drops",
               "new_ts": "02_Timeseries/wells_timeseries_final",
               "coor": "02_Attributes/wells_wrong_coordinates",
               "id": "02_Attributes/wells_dup_id",
               "duploc": "02_Attributes/wells_duplicate_loc_time",
               },
}

## Merge all attributes excel files and sheets together to one large attributes table and discard "ID" and "Country" duplicates
"""The groundwater time series can only be merged to the attributes via ID and country. If there are duplicates, they cannot be
correctly assigned. That is why "ID" and "country" duplicates must be discarded."""

folders = os.listdir(
    config["basepath"] + config["wells"]
)  # country folders; each contains attribute tables

#folders_filt = [item for item in folders if 'GGMN' not in item]

wells = []
wells_drop = []  # to store discarded well duplicates
all_wells = [] # to merge all unprocessed attribute tables
total_number = 1  # to count all wells
all_kept = 0  # to count all kept wells
all_dropped = 0

# loop over every country folder

for fold in folders:
    # Import well.xlsx table
    file = config["basepath"] + config["wells"] + "/" + fold + "/" + "wells.ods"
    # "General Information" sheet in wells.xlsx
    well = pd.read_excel(file, engine='odf', sheet_name=0, skiprows= [1,1], dtype={"ID":object})
    all_wells.append(well)
    total_number = total_number + len(well)  # to count all existing wells in the dataset
    well.rename(columns={"Ground surface elevation": "provider_ground_elevation_m_asl","DEM elevation based on the GLO_90m dataset":"GLO_90m_elevation_m",
                         "Top of well elevation": "top_of_well_elevation_m_asl","Unnamed: 17":"License_restriction",
                         "Measurement Type": "groundwater_level", "Unnamed: 19": "groundwater_quality",
                         "Measurement Data":"first_date_of_measurement", "Unnamed: 21":"last_date_of_measurement"},inplace=True)  # rename some columns
    well.loc[well["Unnamed: 9"] == "ft", "provider_ground_elevation_m_asl"] = well.loc[well["Unnamed: 9"] == "ft", "provider_ground_elevation_m_asl"] * 0.3048  # convert all feet values to meter
    well.loc[well["Unnamed: 11"] == "ft", "GLO_90m_elevation_m"] = well.loc[well["Unnamed: 11"] == "ft", "GLO_90m_elevation_m"] * 0.3048  # convert all feet values to meter
    well.loc[well["Unnamed: 13"] == "ft", "top_of_well_elevation_m_asl"] = well.loc[well["Unnamed: 13"] == "ft", "top_of_well_elevation_m_asl"] * 0.3048  # convert all feet values to meter
    well.drop(["Unnamed: 9","Unnamed: 11","Unnamed: 13"], axis=1, inplace=True)  # remove unit column because everything is in meter now
    # "Hydrogeology" sheet in wells.xlsx
    hydrogeo = pd.read_excel(file, engine='odf', sheet_name=1, skiprows=[1, 1], dtype={"ID":object})
    hydrogeo.rename(columns={"Unnamed: 10": "Unit_HC", "Unnamed: 12": "Unit_T", "Unnamed: 14": "Unit_Sp",
                                 "Unnamed: 16": "Unit_SC", "Unnamed: 18":"Unit_St"},inplace=True)
    # "Management" sheet in wells.xlsx
    man = pd.read_excel(file, engine='odf', sheet_name=2, skiprows=[1, 1], dtype={"ID":object})
    man.rename(columns={"Unnamed: 3":"manager","Unnamed: 4":"org_description","License":"lic_Name",
                            "Unnamed: 8": "lic_validfrom","Unnamed: 9": "lic_validtil", "Unnamed: 10": "lic_description"}, inplace=True)
    # "Drilling and Construction" sheet in drilling_and_construction.ods; Other sheets in this file contain almost no information (1-2 wells with entries which are test entries)
    file = config["basepath"] + config["wells"] + "/" + fold + "/" + "drilling_and_construction.ods"
    construction = pd.read_excel(file, engine='odf', sheet_name=0, skiprows=[1, 1], dtype={"ID":object})
    construction.rename(columns={"Pump":"pump_installer","Unnamed: 10":"pump_description"," Total depth": "total_drilling_depth_m"},inplace=True)
    construction.loc[construction["Unit"] == "ft", "total_drilling_depth_m"] = construction.loc[construction["Unit"] == "ft", "total_drilling_depth_m"] * 0.3048  # convert all feet values to meter
    construction.drop(columns={"Unit","Year of drilling"},inplace=True) # unit is meters anyway and year of drilling is only provided for 12 wells, both are dropped to keep the table as small as possible
    # merge all attributes tables together
    well_mer1 = pd.merge(well, hydrogeo, on=["ID", "Name "], how="left")
    well_mer2 = pd.merge(well_mer1, man, on=["ID", "Name "], how="left")
    well_full = pd.merge(well_mer2, construction, on=["ID", "Name "], how="left")
    # drop full duplicates and duplicates by ID and country
    well_full = well_full.drop_duplicates()
    subset = well_full[well_full.duplicated(subset=["ID"], keep=False)]
    if not subset.empty:
        # we know that some duplicates are created because records from the Jasechko et al. (2024) study and G3P campaign are duplicates
        # to save some records, we delete the records from Jasechko and G3P from the duplicates subset
        orgs_to_exclude = ["Jasechko", "Wells for G3P evaluation"]
        pattern = "|".join(orgs_to_exclude)
        exclude_mask = subset["Organisation"].str.contains(
            pattern, na=False, regex=True
        )
        excluded_indices = subset[exclude_mask].index
        # drop rows based on organisation
        well_full = well_full.drop(index=excluded_indices, errors="ignore")
        # drop rows based on description
        if subset['Description'].notna().all():
            well_full = well_full.drop(index=subset[subset["Description"].str.contains("G3P", na=False)].index)
        # Drop all remaining duplicates as we can't determine which one to keep
        well_full = well_full.drop_duplicates(subset=["ID"], keep=False)

    all_kept = all_kept + len(well_full) # to count all kept wells in the dataset
    wells.append(well_full)
    wells_drop.append(subset)

# number of lost records in th attributes table due to duplicates
all_dropped = total_number - all_kept
# export partitioned lists of wells
# all wells
pd.concat(all_wells).to_csv(config["basepath"] + config["output"]["all_unprocessed"], sep=";", index=False)
# wells without ID-Country duplicates
# remove final duplicates that appear due to the merge of GGMN and Groundwater Observation Repository
wells = pd.concat(wells)
wells = wells.drop_duplicates(subset=["ID","Country"], keep=False)
wells.to_csv(config["basepath"] + config["output"]["all"], sep=";", index=False)
# dropped duplicate wells
pd.concat(wells_drop).to_csv(config["basepath"] + config["output"]["all_dups"], sep=";", index=False)


## Preprocess attributes table
# Import full attributes table
all_wells = pd.read_csv(
    config["basepath"] +
    config["output"]["all"],
    sep=";",
    dtype={"ID": object},
)
# ID must be in string format so that the merge works
all_wells.columns = all_wells.columns.str.lower()  # lowercase all column names
all_wells.rename(columns={"id": "ID"}, inplace=True)  # rename ID column again

# Merge attributes table with time series attributes by ID and Country
wells_timeseries_attributes = pd.read_csv(
    config["basepath"] +
    config["timeseries_att"],
    sep=";",
    dtype={"ID": object},
)
wells_timeseries_attributes.drop(columns=["length_years_raw","autocorrelation_raw",
       "autocorrelation_agg", "diff_x_percentile ","var_raw",
       "var_agg", "trim_gaplength", "trim_gapamount",
       "sign_cat_raw", "sign_cat_agg", "outliers_raw",
       "outliers_agg","plateaus_raw", "plateaus_agg",
       "trend_raw", "trend_agg"], inplace=True)

# get length of all time series after script 01
all_ts = len(wells_timeseries_attributes)
# ID must be in string format so that the merge works
wells_timeseries_attributes.ID = wells_timeseries_attributes.ID.astype("str")
# remove all duplicates by ID and Country in time series so that a clear assignment is possible
wells_timeseries_attributes = wells_timeseries_attributes.drop_duplicates(subset=["ID", "country"], keep=False)
# actual merge with attributes table
wells_merge = pd.merge(all_wells, wells_timeseries_attributes, on=["ID", "country"], how="inner")

# number of lost time series due to the merge (= lost because of "ID"&"Country" duplicates)
lost_merge = all_ts - len(wells_merge)

# Overwrite -9999.-999 and 0 for provider elevation and top of well elevation
wells_merge.provider_ground_elevation_m_asl[(wells_merge.provider_ground_elevation_m_asl==-999) |
                                            (wells_merge.provider_ground_elevation_m_asl==-9999) ] = None

wells_merge.provider_ground_elevation_m_asl[(wells_merge.top_of_well_elevation_m_asl==-999) |
                                            (wells_merge.top_of_well_elevation_m_asl==-9999) ] = None

# Drop duplicates by coordinates, starting_date, ending_date and mean groundwater table
"""They are duplicate time series in the original data but with different ID's. 
Here, they are found over the exact same coordinates, starting and ending date and mean groundwater table."""

subset = wells_merge[
    wells_merge.duplicated(
        subset=[
            "latitude",
            "longitude",
            "starting_date",
            "ending_date",
            "groundwater_mean_m",
        ]
    )
]
wells_drop_dup = wells_merge.drop(index=subset.index)  # attributes without duplicates

# number of lost time series due to duplicates by coordinates, starting_date, ending_date and mean groundwater table
lost_duplicate_location = (len(wells_merge) - len(wells_drop_dup))

# Drop duplicates from Jasechko study
dummy = wells_drop_dup.copy(deep=True)
dummy["lon_rounded"] = dummy["longitude"].round(1)
dummy["lat_rounded"] = dummy["latitude"].round(1)
dummy["gw_rounded"] = dummy["groundwater_mean_m"].round(0)

subset = dummy[dummy.duplicated(subset=["lon_rounded","lat_rounded","gw_rounded"], keep=False)]
filtered_df = subset.groupby(["lon_rounded","lat_rounded","gw_rounded"]).filter(lambda x: x["organisation"].nunique() > 1)
jasechko_dup = filtered_df[filtered_df.organisation=="Jasechko et al. (2024) Rapid groundwater decline and some cases of recovery in aquifers globally"]

wells_drop_jas = wells_drop_dup.drop(index=jasechko_dup.index)

lost_jasechko_duplicates = len(wells_drop_dup) - len(wells_drop_jas)

# Drop wells in wrong coordinate systems
"""In WGS 84, the latitude coordinates should be between -90 and 90. Longitude coordinates should be between -180 and 180.
Wells with coordinates outside this range are discarded here."""

wells_final = wells_drop_jas.drop(wells_drop_jas[(wells_drop_jas.latitude > 90) | (wells_drop_jas.latitude < -90)].index)
wells_final = wells_final.drop(wells_final[(wells_final.longitude > 180) | (wells_final.longitude < -180)].index)
wells_final.reset_index(inplace=True, drop=True)
wells_final.ID = wells_final.ID.astype("str")

# number of lost time series due to wrong coordinates
lost_wrong_coordinates = len(wells_drop_jas) - len(wells_final)

# make everything clean and tidy
# rename columns
wells_final.rename(columns={
    "ID": "original_ID_groundwater",
    "name ":"name",
    "aquifer name": "aquifer_name",
    "feature type": "feature_type",
    "outliers": "outliers_change_points",
    "sign_cat":"negative_signs_wtd",
}, inplace=True)

# test entries by IGRAC are removed
lost_test_wells = len(wells_final[wells_final.original_ID_groundwater=="whostest123"])
wells_final = wells_final[wells_final.original_ID_groundwater!="whostest123"]
# drop all empty columns
wells_final = wells_final.dropna(axis=1, how='all')
# remove columns that are not wanted in GROW
wells_final.drop(columns=[
    "glo_90m_elevation_m",
    "groundwater_level",
    "groundwater_quality",
    "first_date_of_measurement",
    "last_date_of_measurement",
    "aquifer type"],inplace=True)
# keep column license_restriction only if there is more information than "-"
if len(set(wells_final["license_restriction"])) < 2:
    wells_final.drop(columns=["license_restriction"], inplace=True)
# reorder columns
wells_final = pd.concat([wells_final.iloc[:,:10],
                        wells_final["total_drilling_depth_m"],
                        wells_final.iloc[:,10:12],
                        wells_final.iloc[:,13:17],
                        wells_final["license"],
                        wells_final[["interval","aggregated_from_n_values_median"]],
                        wells_final.iloc[:,19:22],
                        wells_final["autocorrelation"],
                        wells_final.iloc[:,25:31],
                        wells_final["reference_point"],
                        wells_final.iloc[:,31:]],axis=1)

# add unique GROW ID
wells_final.reset_index(inplace=True, drop=True)
wells_final["GROW_ID"] = None
for i in range(len(wells_final)):
    # "GROW-" + a timestamp is used as ID
    wells_final.loc[i,"GROW_ID"] = f"GROW-{str(int(time.process_time_ns()))}"

# change position of ID column
growid = wells_final.pop("GROW_ID")
wells_final.insert(0, "GROW_ID", growid)

# export groundwater attributes table
wells_final.to_csv(
    config["basepath"]
    + config["output"]["filtered"]
    + config["output"]["name"]
    + ".txt",
    sep=";",
    index=False,
)

# export number and percentage of lost time series during attributes processing
drops = pd.DataFrame(
    {
        "name": [
            "n_all_wells",
            "lost_duplicate_ID_country",
            "lost_merge",
            "lost_duplicate_location",
            "lost_jasechko_duplicates",
            "lost_wrong_coordinates",
            "lost_test_wells",
        ],
        "percentage": [
            1,
            ((all_dropped) / total_number),
            (lost_merge / all_ts),
            (lost_duplicate_location / all_ts),
            (lost_jasechko_duplicates / all_ts),
            (lost_wrong_coordinates / all_ts),
            (lost_test_wells / all_ts),
        ],
        "number": [
            total_number,
            all_dropped,
            lost_merge,
            lost_duplicate_location,
            lost_jasechko_duplicates,
            lost_wrong_coordinates,
            lost_test_wells,
        ],
    }
)

drops.to_csv(
    config["basepath"] + config["output"]["drops"] + config["output"]["name"] + ".txt",
    sep=";",
    index=False,
)

## final preprocessing of groundwater time series and clean-up

# trim time series data based on kept attributes
# import groundwater time series from "01_processing_gw_time_series.py"
timeseries = pd.read_csv(
    config["basepath"]
    + config["timeseries"],
    sep=";",
    dtype={"ID":"str"}, # column must be string for the merge to work
)

# merge GROW-ID to time series by original ID and Country
timeseries_trimmed = pd.merge(
    timeseries,
    wells_final[["GROW_ID", "original_ID_groundwater", "country"]],
    left_on=["ID", "country"],
    right_on=["original_ID_groundwater", "country"],
    how="inner",
)
timeseries_trimmed.drop(["original_ID_groundwater", "ID"], axis=1, inplace=True)

# make new ID column first column
growid = timeseries_trimmed.pop('GROW_ID')
timeseries_trimmed.insert(0, 'GROW_ID', growid)

# make three separate columns for the three different reference points
# groundwater depth [from ground /from top of the well] and groundwater level
par_sep = timeseries_trimmed.pivot(columns="parameter", values=["groundwater", "groundwater_filled"])
par_sep.columns = [
    "groundwater_depth_from_ground_elevation_m",
    "groundwater_depth_from_well_top_elevation_m",
    "groundwater_water_level_elevation_m_asl",
    "groundwater_filled_depth_from_ground_elevation_m",
    "groundwater_filled_depth_from_well_top_elevation_m",
    "groundwater_filled_water_level_elevation_m_asl",
]
# change reference_point term to column name
timeseries_trimmed = pd.concat([timeseries_trimmed, par_sep], axis=1).drop(
    columns=["parameter", "groundwater", "groundwater_filled"])

# create datetime columns that are later used for the merge of
# the Earth system variables and groundwater time series
# add second time column so that later Earth system variables with
# monthly resolution can be merged with daily groundwater time series
date2 = timeseries_trimmed["date"]
timeseries_trimmed.insert(1, "date2", date2)
# add month as column
timeseries_trimmed["month"] = (
    pd.to_datetime(timeseries_trimmed["date"], format="%Y-%m-%d")
    .dt.to_period("M")
    .astype("str")
)
timeseries_trimmed.date2[timeseries_trimmed.interval == "d"] = pd.to_datetime(
    timeseries_trimmed.month[timeseries_trimmed.interval == "d"], format="%Y-%m"
)
timeseries_trimmed["date2"] = pd.to_datetime(timeseries_trimmed["date2"])
# add year as column
timeseries_trimmed["year"] = timeseries_trimmed["date2"].dt.year.astype("int")

# Export final time series table
timeseries_trimmed.to_csv(
    config["basepath"] + config["output"]["new_ts"] + config["output"]["name"] + ".txt",
    sep=";",
    index=False,
)
