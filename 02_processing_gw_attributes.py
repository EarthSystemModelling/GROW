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
    "timeseries_att": "02_Timeseries/wells_timeseries_attributes_V08_small.txt", # time series attributes derived in "01_processing_gw_time_series"
    "timeseries" : "02_Timeseries/wells_timeseries_V08_small.txt", # time series table derived in "01_processing_gw_time_series"
    # paths of exported output files
    "output": {"name": "_V08_small", # name of version
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

wells = []
wells_drop = []  # to store discarded well duplicates
all_wells = [] # to merge all unprocessed attribute tables
total_number = 0  # to count all wells
all_kept = 0  # to count all kept wells

# loop over every country folder
for fold in folders:
    # Import well.xlsx table
    file = config["basepath"] + config["wells"] + "/" + fold + "/" + "wells.ods"
    # "General Information" sheet in wells.xlsx
    well = pd.read_excel(file, engine='odf', sheet_name=0, skiprows= [1,1], dtype={"ID":object})
    all_wells.append(well)
    total_number = total_number + len(well)  # to count all existing wells in the dataset
    well.rename(columns={"Ground surface elevation": "provider_elevation_m","DEM elevation based on the GLO_90m dataset":"GLO_90m_elevation_m",
                         "Top of well elevation": "top_of_well_elevation_m","Unnamed: 17":"License_restriction",
                         "Measurement Type": "groundwater_level", "Unnamed: 19": "groundwater_quality",
                         "Measurement Data":"first_date_of_measurement", "Unnamed: 21":"last_date_of_measurement"},inplace=True)  # rename some columns
    well.loc[well["Unnamed: 9"] == "ft", "provider_elevation_m"] = well.loc[well["Unnamed: 9"] == "ft", "provider_elevation_m"] * 0.3048  # convert all feet values to meter
    well.loc[well["Unnamed: 11"] == "ft", "GLO_90m_elevation_m"] = well.loc[well["Unnamed: 11"] == "ft", "GLO_90m_elevation_m"] * 0.3048  # convert all feet values to meter
    well.loc[well["Unnamed: 13"] == "ft", "top_of_well_elevation_m"] = well.loc[well["Unnamed: 13"] == "ft", "top_of_well_elevation_m"] * 0.3048  # convert all feet values to meter
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
    construction.rename(columns={"Pump":"pump_installer","Unnamed: 10":"pump_description"," Total depth": "drilling_total_depth_m"},inplace=True)
    construction.loc[construction["Unit"] == "ft", "drilling_total_depth_m"] = construction.loc[construction["Unit"] == "ft", "drilling_total_depth_m"] * 0.3048  # convert all feet values to meter
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
        # DN: TODO
        if subset['Description'].notna().all(): # Das gibt sonst später eine Fehlermeldung, wenn die ganze Spalte NA ist
            well_full = well_full.drop(index=subset[subset["Description"].str.contains("G3P", na=False)].index)
        # Drop all remaining duplicates as we can't determine which one to keep
        well_full = well_full.drop_duplicates(subset=["ID"], keep=False)

    all_kept = all_kept + len(well_full) # to count all kept wells in the dataset
    wells.append(well_full)
    wells_drop.append(subset)

# number of lost records in th attributes table due to duplicates
all_dropped = total_number - all_kept
# export partitioned lists of wells
pd.concat(wells).to_csv(
    config["basepath"] + config["output"]["all"], sep=";", index=False
)
pd.concat(wells_drop).to_csv(
    config["basepath"] + config["output"]["all_dups"], sep=";", index=False
)

## Preprocess attributes table

all_wells = pd.read_csv(
    config["basepath"] + config["output"]["all"], sep=";", dtype={"ID": object}
)  # import attributes table
all_wells.columns = all_wells.columns.str.lower()  # lowercase all column names
all_wells.rename(columns={"id": "ID"}, inplace=True)  # rename ID column again
len_all = 251398  # total amount of time series after preprocessing in "01_processing_gw_time_series.py"; to calculate percentage loss per processing step

# Merge attributes table with time series attributes by ID and Country
wells_timeseries_attributes = pd.read_csv(
    config["basepath"] +
    config["timeseries_att"]
    , sep=";"
)
# ID must be in string format so that the merge works
wells_timeseries_attributes.ID = wells_timeseries_attributes.ID.astype("str")
# remove all duplicates by ID and Country in time series so that a clear assignment is possible
wells_timeseries_attributes = wells_timeseries_attributes.drop_duplicates(subset=["ID", "country"], keep=False)
wells_merge = pd.merge(all_wells, wells_timeseries_attributes, on=["ID", "country"], how="inner")  # actual merge

# number of lost time series due to the merge (= lost because of "ID"&"Country" duplicates)
lost_merge = (len(wells_timeseries_attributes) - len(wells_merge))

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

# Drop duplicates from Jasechko study that have rounded coordinates
# I will probably add that in the near future

# Drop wells in wrong coordinate systems
"""In WGS 84, the latitude coordinates should be between -90 and 90. Longitude coordinates should be between -180 and 180.
Wells with coordinates outside this range are discarded here."""

wells_filt = wells_drop_dup.drop(wells_drop_dup[
        (wells_drop_dup.latitude > 90) | (wells_drop_dup.latitude < -90)
    ].index
)
wells_filt = wells_filt.drop(wells_filt[(wells_filt.longitude > 180) | (wells_filt.longitude < -180)].index)
wells_filt.reset_index(inplace=True, drop=True)
wells_filt.ID = wells_filt.ID.astype("str")

# number of lost time series due to wrong coordinates
lost_wrong_coordinates = len(wells_drop_dup) - len(wells_filt)

# make everything clean and tidy
# rename columns
wells_filt.rename(columns={
    "ID": "original_ID_groundwater",
    "name ":"name",
    "aquifer name": "aquifer_name",
    "feature type": "feature_type"
}, inplace=True)
# test entries by IGRAC are removed
wells_filt = wells_filt[wells_filt.original_ID_groundwater!="whostest123"]
# drop all empty columns
wells_filt = wells_filt.dropna(axis=1, how='all')
# remove columns that are not wanted in GROW
wells_filt.drop(columns=[
    "glo_90m_elevation_m",
    'groundwater_level',
    'groundwater_quality',
    'first_date_of_measurement',
    'last_date_of_measurement',
    'aquifer type'],inplace=True)
# keep column license_restriction only if there is more information than "-"
if len(set(wells_filt["license_restriction"])) < 2:
    wells_filt.drop(columns=["license_restriction"], inplace=True)
wells_filt.insert(29, "reference_point", wells_filt.pop('reference_point')) # reorder this column

# add unique GROW ID
# TODO@DN: change to uuid?
wells_filt.reset_index(inplace=True, drop=True)
wells_filt["GROW_ID"] = None
for i in range(len(wells_filt)):
    # "GROW-" + a timestamp is used as ID
    wells_filt.loc[i,"GROW_ID"] = f"GROW-{str(int(time.process_time_ns()))}"

# change position of ID column
growid = wells_filt.pop('GROW_ID')
wells_filt.insert(0, 'GROW_ID', growid)

# export groundwater attributes table
wells_filt.to_csv(
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
            "lost_wrong_coordinates",
        ],
        "percentage": [
            1,
            ((all_dropped) / total_number),
            (lost_merge / len_all),
            (lost_duplicate_location / len_all),
            (lost_wrong_coordinates / len_all),
        ],
        "number": [
            total_number,
            all_dropped,
            lost_merge,
            lost_duplicate_location,
            lost_wrong_coordinates,
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
    + config["timeseries"]
    , sep=";"
)
# column must be string for the merge to work
timeseries.ID = timeseries.ID.astype("str")

# merge with GROW-ID by old_ID and Country
timeseries_trimmed = pd.merge(
    timeseries,
    wells_filt[["GROW_ID", "original_ID_groundwater", "country"]],
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
    "groundwater_depth_from_top_elevation_m",
    "groundwater_water_level_m_asl",
    "groundwater_filled_depth_from_ground_elevation_m",
    "groundwater_filled_depth_from_top_elevation_m",
    "groundwater_filled_water_level_m_asl",
]
# change reference_point term to column name
timeseries_trimmed = pd.concat([timeseries_trimmed, par_sep], axis=1).drop(
    columns=["parameter", "groundwater", "groundwater_filled"]
)

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
