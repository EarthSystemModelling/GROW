### Takes all well attributes tables from every country and merges it in one + attributes preprocessing

## Configuration
config = {
    "basepath" : "/mnt/storage/grow/Groundwater/",
    "wells": "Well_And_Monitoring_Data",
    "timeseries_att": "Wells_timeseries/wells_timeseries_attributes_V05.txt",
    "timeseries" : "Wells_timeseries/wells_timeseries_V05.txt",
    "output": {"name": "_V06",
               "all":"Wells_attributes/wells_attributes_all.csv",
               "all_dups": "Wells_attributes/wells_attributes_all_dups.txt",
               "filtered": "Wells_attributes/wells_attributes",
               "points": "Other_shapes/wells_points",
               "drops": "Statistics/wells_attributes_drops",
               "new_ts": "Wells_timeseries/wells_timeseries_final",
               "coor": "Wells_attributes/wells_wrong_coordinates",
               "id": "Wells_attributes/wells_dup_id",
               "merge": "Wells_attributes/wells_pulling_empty",
               "duploc": "Wells_attributes/wells_duplicate_loc_time"}
}

import os
import time
import pandas as pd
import geopandas as gpd

## Merge all attributes together

folders = os.listdir(config["basepath"] + config["wells"])

wells = []
wells_drop = []
total_number = 0
all_kept = 0
for fold in folders:
    file = config["basepath"] + config["wells"] + "/" + fold + "/" + "wells.xlsx"
    # "General Information" sheet in wells.xlsx
    well = pd.read_excel(file, engine='openpyxl', sheet_name=0, skiprows= [1,1], dtype={"ID":object})
    well.rename(columns={"Ground surface elevation\n": "surface_elevation_m", "Top of well elevation": "top_of_well_elevation_m","Unnamed: 15":"License_restriction"},inplace=True)  # rename some columns
    well.loc[well["Unnamed: 9"] == "ft", "surface_elevation_m"] = well.loc[well["Unnamed: 9"] == "ft", "surface_elevation_m"] * 0.3048  # convert all feet values to meter
    well.loc[well["Unnamed: 11"] == "ft", "top_of_well_elevation_m"] = well.loc[well["Unnamed: 11"] == "ft", "top_of_well_elevation_m"] * 0.3048  # convert all feet values to meter
    well.drop(["Unnamed: 9","Unnamed: 11"], axis=1, inplace=True)  # remove unit column because everything is in meter now
    # "Hydrogeology" sheet in wells.xlsx
    hydrogeo = pd.read_excel(file, engine='openpyxl', sheet_name=1, skiprows=[1, 1], dtype={"ID":object})
    hydrogeo.rename(columns={"Unnamed: 9": "Unit_HC", "Unnamed: 11": "Unit_T", "Unnamed: 13": "Unit_Sp",
                                 "Unnamed: 15": "Unit_SC", "Unnamed: 17":"Unit_St"},inplace=True)
    # "Management" sheet in wells.xlsx
    man = pd.read_excel(file, engine='openpyxl', sheet_name=2, skiprows=[1, 1], dtype={"ID":object})
    man.rename(columns={"Unnamed: 2":"Manager","Unnamed: 3":"Org_description","License":"Lic_Name",
                            "Unnamed: 7": "Lic_validfrom","Unnamed: 8": "Lic_validtil", "Unnamed: 9": "Lic_description"}, inplace=True)
    # "Drilling and Construction" sheet in drilling_and_construction.xlsx; Other sheets in this file contain almost no information (1-2 wells with entries which are test entries)
    file = config["basepath"] + config["wells"] + "/" + fold + "/" + "drilling_and_construction.xlsx"
    construction = pd.read_excel(file, engine='openpyxl', sheet_name=0, skiprows=[1, 1], dtype={"ID":object})
    construction.rename(columns={"Pump":"Pump_installer","Unnamed: 9":"Pump_description"," Total depth": "drilling_total_depth_m"},inplace=True)
    construction.drop(columns={"Unit","Year of drilling"},inplace=True) # unit is meters anyway and year of drilling is only provided for 12 wells, both are dropped to keep the table as small as possible
    # merge all together
    well_mer1 = pd.merge(well,hydrogeo,on="ID",how="left")
    well_mer2 = pd.merge(well_mer1, man, on="ID",how="left")
    well_full = pd.merge(well_mer2, construction, on="ID",how="left")
    # drop full duplicates and duplicates by ID and country
    total_number = total_number + len(well_full) # to count all existing wells in the dataset
    well_full = well_full[~well_full.duplicated()] # drop full (identical in all columns) duplictates
    subset = well_full[well_full.duplicated(subset=["ID"], keep=False)]  # subset all duplicates by ID
    if subset.empty == False:
        well_full = well_full.loc[~well_full.index.isin(subset[subset['Organisation'].str.contains("Jasechko")].index),:] # Deleting Jasechko records in duplicates
        well_full = well_full.loc[~well_full.index.isin(subset[subset['Organisation'].str.contains("Wells for G3P evaluation")].index),:] # Deleting G3P records in duplicates
        subset.dropna(subset="Description",inplace= True) # Delete NA records in Description column, necessary for next step
        well_full = well_full.loc[~well_full.index.isin(subset[subset['Description'].str.contains("G3P")].index),:] # Deleting G3P records in duplicates
        subset = well_full[well_full.duplicated(subset=["ID"], keep=False)] # Look for duplicates again
        well_full.drop(subset.index, axis=0, inplace=True)  # Delete all duplicate ID's from well attributes table
    all_kept = all_kept + len(well_full)
    wells.append(well_full)
    wells_drop.append(subset)

all_dropped = total_number - all_kept
pd.concat(wells).to_csv(config["basepath"] + config["output"]["all"],sep=";",index=False)
pd.concat(wells_drop).to_csv(config["basepath"] + config["output"]["all_dups"],sep=";",index=False)

## Preprocess: Drop incorrect/duplicate wells

all_wells = pd.read_csv(config["basepath"] + config["output"]["all"],sep=";", dtype={"ID":object})
all_wells.columns = all_wells.columns.str.lower() # lowercase all column names
all_wells.rename(columns={"id":"ID"},inplace=True)
len_all = 221123 # this should be the total amount of time series

# Merge with time series data by ID and Country
wells_timeseries_attributes = pd.read_csv(config["basepath"] + config["timeseries_att"], sep=";")
wells_timeseries_attributes.ID = wells_timeseries_attributes.ID.astype("str")
wells_timeseries_attributes =  wells_timeseries_attributes[~wells_timeseries_attributes.duplicated(subset=["ID","country"], keep=False)] # remove all duplicates by ID and Country in time series (8 are lost)
wells_merge = pd.merge(all_wells, wells_timeseries_attributes, on = ["ID","country"], how = "inner")

lost_merge = len(wells_timeseries_attributes)-len(wells_merge) # should be the number as the sum of all time series drops

# Drop duplicate by coordinates, time range and mean value
subset = wells_merge[wells_merge.duplicated(subset=["latitude","longitude","firstdate","lastdate","groundwater_mean_[m]"])]
wells_drop_dup = wells_merge.drop(subset.index, axis=0)

lost_duplicate_location = len(wells_merge) - len(wells_drop_dup)

# Drop wells in wrong coordinate systems
wells_filt = wells_drop_dup.drop(wells_drop_dup[(wells_drop_dup.latitude > 90) | (wells_drop_dup.latitude < -90)].index)
wells_filt = wells_filt.drop(wells_filt[(wells_filt.longitude > 180) | (wells_filt.longitude < -180)].index)
wells_filt.reset_index(inplace=True, drop=True)
wells_filt.ID = wells_filt.ID.astype("str")

lost_wrong_coordinates = len(wells_drop_dup) - len(wells_filt)

# make everything clean and tidy
wells_filt.rename(columns={"ID": "ID_old"}, inplace=True)  # rename ID column
wells_filt.dropna(axis=1, how='all', inplace=True) # drop completely empty columns
wells_filt = wells_filt[wells_filt.ID_old!="whostest123"] # this is just a test entry
wells_filt = wells_filt.dropna(axis=1, how='all') # drop all empty columns
wells_filt.loc[(wells_filt["surface_elevation_m"]==-999),"surface_elevation_m"] = None
wells_filt.loc[(wells_filt["top_of_well_elevation_m"]==-999),"top_of_well_elevation_m"] = None

# add GROW ID
wells_filt.reset_index(inplace=True, drop=True)
wells_filt["ID"] = None
for i in range(len(wells_filt)):
    wells_filt.loc[i,"ID"] = f"GROW-{str(int(time.process_time_ns()))}"

# change position of ID column
growid = wells_filt.pop('ID')
wells_filt.insert(0, 'ID', growid)

# export table
wells_filt.to_csv(config["basepath"] + config["output"]["filtered"] + config["output"]["name"] + ".txt", sep=";", index=False)

# export point coordinates
wells_gdf = gpd.GeoDataFrame(wells_filt, geometry= gpd.points_from_xy(wells_filt.longitude, wells_filt.latitude), crs="EPSG:4326")
wells_gdf.to_file(config["basepath"] + config["output"]["points"] + config["output"]["name"] + ".json")

# export drop statistics
drops = pd.DataFrame({"name":["n_all_wells","lost_duplicate_ID_country","lost_merge","lost_duplicate_location","lost_wrong_coordinates"],
                      "percentage":[1,((all_dropped)/total_number),(lost_merge/len_all),(lost_duplicate_location/len_all),(lost_wrong_coordinates/len_all)],
                      "number": [total_number,all_dropped,lost_merge,lost_duplicate_location,lost_wrong_coordinates]})
drops.to_csv(config["basepath"] + config["output"]["drops"] + config["output"]["name"] + ".txt",sep=";", index=False)

## trim time series data based on static preprocessing

timeseries = pd.read_csv(config["basepath"] + config["timeseries"], sep=";")
timeseries.ID = timeseries.ID.astype("str") # so that the merge works

# merge with GROW-ID by old_ID and Country
timeseries_trimmed = pd.merge(timeseries, wells_filt[["ID","ID_old","country"]], left_on=["ID","country"],right_on=["ID_old","country"], how="inner")
timeseries_trimmed.drop(["ID_old","ID_x"],axis=1,inplace=True)

growid = timeseries_trimmed.pop('ID_y')
timeseries_trimmed.insert(0, 'ID', growid)

# make three separate columns for groundwater depth [from ground /from top of the well] and groundwater level
par_sep = timeseries_trimmed.pivot(columns='parameter', values=['groundwater',"groundwater_filled"])
par_sep.columns = ['groundwater_depth_from_ground_elevation_m', 'groundwater_depth_from_top_elevation_m', 'groundwater_water_level_m_asl','groundwater_filled_depth_from_ground_elevation_m', 'groundwater_filled_depth_from_top_elevation_m', 'groundwater_filled_water_level_m_asl']
timeseries_trimmed = pd.concat([timeseries_trimmed, par_sep], axis=1).drop(columns=["parameter","groundwater","groundwater_filled"])

# add second time column
date2 = timeseries_trimmed["date"]
timeseries_trimmed.insert(1, "date2", date2)
timeseries_trimmed["month"] = pd.to_datetime(timeseries_trimmed['date'], format='%Y-%m-%d').dt.to_period('M').astype("str")
timeseries_trimmed.date2[timeseries_trimmed.interval == "d"] = pd.to_datetime(timeseries_trimmed.month[timeseries_trimmed.interval == "d"],format="%Y-%m")
timeseries_trimmed["date2"] = pd.to_datetime(timeseries_trimmed["date2"])
timeseries_trimmed["year"] = timeseries_trimmed["date2"].dt.year.astype("int")

# Export final time serie stable
timeseries_trimmed.to_csv(config["basepath"] + config["output"]["new_ts"] + config["output"]["name"] + ".txt",sep=";",index=False)
