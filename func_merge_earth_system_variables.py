"""This file contains all functions that are used in "03_merge_earth_system_variables.py"."""

# import packages
from datetime import datetime # built-in package
from dateutil.relativedelta import relativedelta # built-in package
import pandas as pd # imported version: 2.2.3
import os # built-in package
import xarray as xr # imported version: 2025.3.0
from netCDF4 import Dataset # imported version: 1.7.2
import geopandas as gpd # imported version: 1.0.1

## get_paths
'''Derives file paths with certain string pattern in folders and subfolders.

path: folder in which to search
pattern: string pattern
pattern2: second string pattern can be given if full=True
full: pattern and pattern2 need to occur in file name'''

def get_paths(path,pattern,pattern2=None,full=False):
    txt_files = []
    # search in every subdirectory for files with pattern
    for root, dirs, files in os.walk(path):
        for file in files:
            if full:
                # if the two pattern occurs in the file name, the file path is appended to txt_files
                if (pattern in file) & (pattern2 in file):
                    txt_files.append(os.path.join(root, file))
            else:
                # if the pattern occurs in the file name, the file path is appended to txt_files
                if pattern in file:
                    txt_files.append(os.path.join(root, file))
    return(txt_files)

## calc_dd
'''Calculates global drainage density (sum of river lengths/ area of basin) with HydroBASINS and HydroRIVERS.

config: configuration dictionary'''

def calc_dd(config):
    # import HydroBASINS shapefile
    basins = gpd.read_file(config["basepath"] + config["factors"]["basins"])
    basins.to_crs(epsg="ESRI:54012", inplace=True) # reproject to metric crs World-Eckert-IV
    basins["area"] = basins.area # recalculate area with the metric crs in m²
    # calculate absolute difference between original area and calculated area in World-Eckert-IV
    basins["area_diff"] = abs(basins["SUB_AREA"] - (basins["area"] / 1000000))

    # import HydroRivers shapefile
    riv = gpd.read_file(config["basepath"] + config["factors"]["rivers"])
    riv.to_crs(epsg="ESRI:54012", inplace=True) # reproject to metric crs World-Eckert-IV
    # intersect rivers and basins (cut river polylines at basin borders)
    riv_cut = gpd.overlay(basins, riv, how="intersection", keep_geom_type=False)
    riv_cut["riv_len"] = riv_cut.length  # calculate length of every river segment in m

    # derive sum of river lengths per basin
    riv_dis = riv_cut.dissolve(by='HYBAS_ID', aggfunc='sum')
    # merge with basin gdf
    dd = basins.merge(riv_dis[["riv_len"]], on="HYBAS_ID", how="left")
    # calculate "drainage density" per basin
    dd["Drainage_den"] = dd["riv_len"] / dd["area"]
    # set drainage density to None, where area difference is 10 km-2 or larger
    wells.loc[wells.area_diff >= 10, "Drainage_den"] = None

    # export as shapefile
    dd.drop(columns={"area_diff","area"}, axis=1, inplace=True)
    dd.to_file(config["basepath"] + config["factors"]["drain_den"])

# merge_raster_static
'''In this function, raster values of a static target variable are extracted at the well locations and 
merged to the well attributes table.

df: dataframe with well coordinates
rasterfile: raster file path of target variable
col_name: column name of extracted target variable
transform: set to True if coordinate system of raster data is not WGS 84'''

def merge_raster_static(df,rasterfile,col_name = None, transform = False):

    # raster of target variable is imported
    raster = xr.open_dataset(rasterfile, engine="rasterio")

    # transform coordinate system if needed
    if transform:
        raster = raster.rio.reproject("EPSG:4326")

    def extract_well_static(row,values):
        extract = raster.sel({"y": row["latitude"], "x": row["longitude"]}, method="nearest") # dataset is trimmed to the well location
        values.append(pd.Series(extract["band_data"].values)[0]) # value is extracted from the xarray dataset and appended to a list

    values=[]
    # for every well location, the raster value is extracted and appended to a list
    df.apply(lambda row: extract_well_static(row,values), axis=1)
    df[col_name]= pd.Series(values) # the list is transformed into a pd.Series and added as column to the dataframe

## merge_vector_point
'''In this function, polygon values of a static target variable are extracted at the well locations and 
merged to the well attribute table via spatial join.

vector: geopandas dataframe of vector product
df: dataframe with well coordinates
columns: attributes that shall be extracted from the vector product'''

def merge_vector_point(vector,df,columns):
    # convert dataframe with well attributes to geopandas dataframe
    gdf = gpd.GeoDataFrame(df, geometry= gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
    # in case the coordinate system of the target variable product is not WGS 84, it is transformed to WGS 84
    if vector.crs != gdf.crs:
       vector.to_crs(gdf.crs,inplace=True)
    # merge vector product to well attributes table
    joined_gdf = gpd.sjoin(gdf, vector, how="left")
    # trim dataframe to specified columns
    column = list(df.columns) + columns
    joined_gdf = joined_gdf[column]

    return(pd.DataFrame(joined_gdf))


## extract_val

'''In this function, raster values of a time-varying target variable are extracted at a specific well location.
The extracted values are appended to an output file.

config: configuration dictionary
row: row of well dataframe with well coordinates (just one well)
raster: xarray dataset with target variable
bandname: name of band for target variable in raster data
latname: name of latitute dimension in xarray dataset with target variable
lonname: name of longitude dimension in xarray dataset with target variable
col_name: column name of extracted target variable
ts: the time dimension of the target variable in case the time could not be decoded in the xarray dataset
'''

def extract_val(config,row,raster,bandname,latname,lonname,col_name,ts):

    # xarray dataset is trimmed to the well location
    extract = raster.sel({latname: row["latitude"], lonname: row["longitude"]}, method="nearest")

    # create dataframe with ID, time column and extracted value/s
    if ts == None:
        # time is taken from the xarray dataset
        ts_df = pd.DataFrame({"ID": row["ID"], "Time": extract.time.values, "Interval": row["interval"],"Value": extract[bandname].values})
    else:
        # for NDVI and ISIMIP, the time column was created manually and is handed over, here
        ts_df = pd.DataFrame({"ID": row["ID"], "Time": ts, "Interval": row["interval"], "Value": extract[bandname].values})

    # cut time - to save storage, only time stamps that actually occur in the groundwater time series are kept in the dataframe
    # substract one year from "starting_date" so that monthly or yearly products which start at the first day of the year/month are still included
    starting_date = datetime.strptime(row["starting_date"],"%Y-%m-%d").date() - relativedelta(years=1)
    # add one year to "ending_date" so that monthly or yearly products which end with the last day of the year/month are still included
    ending_date = datetime.strptime(row["ending_date"],"%Y-%m-%d").date() + relativedelta(years=1)
    ts_df = ts_df[(ts_df.Time >= str(starting_date)) & (ts_df.Time <= str(ending_date))] # trim time to relevant period

    if ts_df.empty:
        # continue if time stamps is not relevant for well
        pass
    else:
        # append dataframe to file
        ts_df.to_csv(config["basepath"] + config["output"] + col_name + "_Extraction.txt", sep=";", mode="a", index=False,header=False)

## merge_raster_transient

'''In this function, raster values of a time-varying target variable are extracted at well locations.
The extracted values per well location are appended to an output file and not stored in memory.

df: well attributes table
rasterpath: path to folder with raster file/files
config: configuration dictionary
patt: string pattern that is searched for when deriving single raster file paths
patt2: second string pattern that is searched for when deriving single raster file paths if ful= True
full: patt and patt2 need to occur in single raster file name
bandname: name of band for target variable in raster data
latname: name of latitute dimension in raster data
lonname: name of longitude dimension in raster data
col_name: column name of extracted target variable
ndvi: set to True if NDVI data is extracted; time can not be automatically decoded, an alternative is used
transform: set to True if coordinate system of raster data is not WGS 84 
era5: set to True if ERA5-Land data is extracted; longitude dimension needs to be adjusted + time dimension name must be renamed 
isimip: set to True if ISIMIP data is extracted; time can not be automatically decoded, an alternative is used
'''

def merge_raster_transient(df,rasterpath,config, patt=".nc",patt2= None,ful = False,bandname= "band_data",latname="lat",lonname="lon",col_name = None, ndvi = False, transform = False, era5 = False, isimip = False):

    # Delete output file because results are appended and not overwritten
    if os.path.isfile(config["basepath"] + config["output"] + col_name + "_Extraction.txt"):
        os.remove(config["basepath"] + config["output"] + col_name + "_Extraction.txt")

    # get paths of every file under that direction that contains specified string pattern
    rasterfiles = get_paths(rasterpath, patt, patt2, ful)

    # loop over every file and extract results stepwise
    for file in rasterfiles:
        # import raster of target variable (and derive time dimension)
        if ndvi == True:
            raster_nc = xr.open_dataset(file, decode_times=False) # time cannot be decoded automatically
            # the newer ndvi raster contain a mask hiding invalid values
            # xarray is not properly reading this mask, consequently invalid values are not filtered
            # To filter them manually, the mask is read via the package netCDF4 and the invalid values are
            # set to NA via indexing the Dataset with the mask from the netCDF4-read raster
            raster_mask = Dataset(file)
            mask = raster_mask.variables[bandname][:].mask
            raster_nc[bandname].values[mask] = None
            ts = datetime(int(file[-27:-23]), int(file[-23:-21]), int(file[-21:-19]))  # deriving time from file name
        elif isimip == True:
            raster_nc = xr.open_dataset(file, decode_times=False) # time cannot be decoded automatically
            # time dimension is always the same for ISIMIP data, 1901 - 2021 and has yearly resolution
            # generating time column, manually
            start_date = datetime(1901, 1, 1)
            num_years = 121
            ts = [start_date + relativedelta(years=x) for x in range(num_years)]
        else:
            raster_nc = xr.open_dataset(file)
            ts=None # time correctly decoded in xarray dataset, no need for extra variable
            if era5 == True:
                # for ERA5-Land data the longitude dimension must be corrected from 0-360 to -180-180
                raster_nc['longitude'] = ((raster_nc['longitude'] + 180) % 360) - 180
                raster_nc = raster_nc.sortby("longitude")
                # time dimension must be renamed because it is addressed with the name "time" for the rest of the function
                raster_nc = raster_nc.rename({'valid_time': 'time'})

        # Transform coordinate system if needed
        if transform:
            raster_nc = raster_nc.rio.reproject("EPSG:4326")

        # Actual extraction is performed for every row in an apply function
        df.apply(lambda row: extract_val(config,row,raster_nc,bandname,latname,lonname,col_name,ts), axis=1)

## Aggregate and merge
'''In this function, the dataframes of the extracted time-varying target variables are imported, aggregated to the temporal resolution of the groundwater data
and merged to the groundwater time series table.

config: configuration dictionary
df: groundwater time series table
imp_name: column name of target variable
i: index of well split
exp_name: name of output file
daily: If true, the resolution of the target variable is daily
yearly: If true, the resolution of the target variable is yearly'''

def agg_n_mer(config,df,imp_name,i,exp_name, daily=False, yearly = False):
    # import dataframe with extracted raster values
    raster = pd.read_csv(config["basepath"] + config["output"] + imp_name + "_" + str(i) + "_Extraction.txt", sep = ";", header=None)
    raster.columns = ["ID","time","intervall","Value"] # set column names
    raster.time = pd.to_datetime(raster.time, format='ISO8601') # generate datetime column

    # snow_depth cannot be negative
    if imp_name == "snow_depth":
        raster.Value[raster.Value<0] = 0

    # create datetime column which shall be used to aggregate and merge the data; different datetime format dependent on temporal resolution of groundwater data
    raster["merge"] = None

    if yearly:
        # target variable has yearly resolution; all data is aggregated to a yearly resolution
        raster["merge"] = raster["time"].dt.year.astype("int")
    elif daily == False:
        # target variable has monthly resolution
        # data assigned to daily and monthly groundwater time series are aggregated to a monthly resolution
        raster["merge"][(raster.intervall == "MS") | (raster.intervall == "d")] = pd.to_datetime(raster["time"][(raster.intervall == "MS") | (raster.intervall == "d")].dt.to_period('M').astype(str), format="%Y-%m")
        # data assigned to yearly groundwater time series are aggregated to a yearly resolution
        raster["merge"][raster.intervall == "YS"] = pd.to_datetime(raster["time"][raster.intervall == "YS"].dt.year.astype("int").astype(str), format="%Y")
    else:
        # target variable has daily resolution
        # data is aggregated to the resolution of the assigned groundwater time series (daily, monthly and yearly)
        raster["merge"][raster.intervall == "d"] = pd.to_datetime(raster["time"][raster.intervall == "d"].dt.to_period('D').astype(str), format="%Y-%m-%d")
        raster["merge"][raster.intervall == "MS"] = pd.to_datetime(raster["time"][raster.intervall == "MS"].dt.to_period('M').astype(str),format="%Y-%m")
        raster["merge"][raster.intervall == "YS"] = pd.to_datetime(raster["time"][raster.intervall == "YS"].dt.year.astype("int").astype(str),format="%Y")

    # aggregation to means
    raster_agg = raster.groupby(["merge","intervall","ID"]).mean()
    raster_agg.reset_index(inplace=True)
    raster_agg.drop(["intervall","time"],axis=1,inplace=True)

    # merge target variable to groundwater time series
    # dependent on the temporal resolution of the target variable, another column in the groundwater data is used for the merge
    if yearly:
        # year: just the year as integer
        df_merged = pd.merge(df, raster_agg, how='outer', left_on=['year', 'ID'], right_on=['merge', 'ID'])
    elif daily:
        # date: normal time stamp column
        df_merged = pd.merge(df, raster_agg,  how='outer', left_on=['date','ID'], right_on = ['merge','ID'])
    else:
        # date2: here the daily data was trimmed to a monthly resolution (first day of a month)
        df_merged = pd.merge(df, raster_agg, how='outer', left_on=['date2', 'ID'], right_on=['merge', 'ID'])

    # Clean-up and export
    df_merged.drop(['merge'],axis=1, inplace=True)
    df_merged.rename(columns={"Value":imp_name}, inplace=True)
    df_merged.to_csv(config["basepath"] + config["output"]  + exp_name +".txt", sep=";", index=False)
    return (trend, slope)