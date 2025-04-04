## func_merge_earth_system_variables

'''This file contains all functions that are used in "03_merge_earth_system_variables.py".'''

# import packages
from datetime import datetime # internal package
from dateutil.relativedelta import relativedelta # internal package
import pandas as pd # imported version: 2.2.3
import os # internal package
import xarray as xr # imported version: 2025.3.0
import geopandas as gpd # imported version: 1.0.1

## get_paths
'''Derives file paths with certain string pattern in folders and subfolders.

path: folder in which is searched for
pattern: string pattern
pattern2: second string pattern can be given if full=True
full: pattern and pattern2 need to occur in file name'''

def get_paths(path,pattern,pattern2=None,full=False):
    txt_files = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if full:
                if (pattern in file) & (pattern2 in file):
                    txt_files.append(os.path.join(root, file))
            else:
                if pattern in file:
                    txt_files.append(os.path.join(root, file))
    return(txt_files)

## Calculate Drainage Density Map
'''Calculate global drainage density map with HydroBASINS and HydroRIVERS

config: configuration dictionary'''

def calc_dd(config):
    # Calculate drainage density (length of rivers/ area of basin) with Hydrobasins and Hydrorivers
    basins = gpd.read_file(config["basepath"] + config["factors"]["basins"])
    basins.to_crs(epsg=4087, inplace=True)
    basins["area"] = basins.area
    #bias berechnen Relation berechnete und angegebene Fläche

    riv = gpd.read_file(config["basepath"] + config["factors"]["rivers"])
    riv.to_crs(epsg=4087, inplace=True)
    riv_cut = gpd.overlay(basins, riv, how="intersection", keep_geom_type=False)
    riv_cut["riv_len"] = riv_cut.length  # length in km

    # derive sum of river lengths per basin
    riv_dis = riv_cut.dissolve(by='HYBAS_ID', aggfunc='sum')
    # Merge with basin gdf
    dd = basins.merge(riv_dis[["riv_len"]], on="HYBAS_ID", how="left")
    # calculate "drainage density" per basin
    dd["Drainage_den"] = dd["riv_len"] / dd["area"]

    # export
    dd.drop("area", axis=1, inplace=True)
    dd.to_file(config["basepath"] + config["factors"]["drain_den"])

# merge_raster_static
''''''

def merge_raster_static(df,rasterfile,col_name = None,glim=False, transform = False):

    raster = xr.open_dataset(rasterfile, engine="rasterio")

    # just for glim (15 means NAN) so that it is interpolated later
    if glim:
        raster = raster.where(raster['band_data'] != 15)

    # Transform coordinate system if needed
    if transform:
        raster = raster.rio.reproject("EPSG:4326")

    def extract_well_static(row,values):
        extract = raster.sel({"y": row["latitude"], "x": row["longitude"]}, method="nearest")
        values.append(pd.Series(extract["band_data"].values)[0])

    values=[]
    df.apply(lambda row: extract_well_static(row,values), axis=1)
    df[col_name]= pd.Series(values)

## merge_vector_point

def merge_vector_point(vector,df,columns,int=False):
    gdf = gpd.GeoDataFrame(df, geometry= gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
    if vector.crs != gdf.crs:
       vector.to_crs(gdf.crs,inplace=True)
    if int:
        joined_gdf = gpd.sjoin_nearest(gdf, vector, how="left") # nearest to avoid Nan
    else:
        joined_gdf = gpd.sjoin(gdf, vector, how="left")
    column = list(df.columns) + columns
    joined_gdf = joined_gdf[column]
    return(pd.DataFrame(joined_gdf))


## extract_val

'''This function takes a raster(-stack) and a dataframe with coordinates and returns a series or
list of series with the extracted values per coordinate pair.'''

def extract_val(config,row,raster,bandname,latname,lonname,col_name,ts):

    extract = raster.sel({latname: row["latitude"], lonname: row["longitude"]}, method="nearest")

    # make df with id, time and extracted value
    if ts == None:
        ts_df = pd.DataFrame({"ID": row["ID"], "Time": extract.time.values, "Interval": row["interval"],"Value": extract[bandname].values})
    else:
        ts_df = pd.DataFrame({"ID": row["ID"], "Time": ts, "Interval": row["interval"], "Value": extract[bandname].values})

    # cut time
    # substract and add one year so that for monthly or yearly products which start at the first day of the year/month, they are still included
    firstdate = datetime.strptime(row["firstdate"],"%Y-%m-%d").date() - relativedelta(years=1)
    lastdate = datetime.strptime(row["lastdate"],"%Y-%m-%d").date() + relativedelta(years=1)
    ts_df = ts_df[(ts_df.Time >= str(firstdate)) & (ts_df.Time <= str(lastdate))]

    # Continue if time step is not relevant for well
    if ts_df.empty:
        pass
    else:
        # export row to file
        ts_df.to_csv(config["basepath"] + config["output"] + col_name + "_Extraction.txt", sep=";", mode="a", index=False,header=False)

## merge_raster_transient

'''This function takes a raster(-stack) and a dataframe with coordinates and returns a series or
list of series with the extracted values per coordinate pair.'''

def merge_raster_transient(df,rasterpath,config, patt=".nc",patt2= None,ful = False,bandname= "band_data",latname="lat",lonname="lon",col_name = None, ndvi = False, transform = False, era5 = False, isimip = False):

    # Delete output file because results are appended and not overwritten
    if os.path.isfile(config["basepath"] + config["output"] + col_name + "_Extraction.txt"):
        os.remove(config["basepath"] + config["output"] + col_name + "_Extraction.txt")

    # get paths of every netcdf file under that direction
    rasterfiles = get_paths(rasterpath, patt, patt2, ful)

    # loop over every nc-file and extract results stepwise
    for file in rasterfiles:
        # load raster (and time)
        if ndvi == True:
            raster_nc = xr.open_dataset(file, decode_times=False)
            ts = datetime(int(file[-27:-23]), int(file[-23:-21]), int(file[-21:-19]))  # taking time from file name
        elif isimip == True:
            # I can also manually
            raster_nc = xr.open_dataset(file, decode_times=False)
            start_date = datetime(1901, 1, 1)
            num_years = 121
            ts = [start_date + relativedelta(years=x) for x in range(num_years)]
        else:
            raster_nc = xr.open_dataset(file)
            ts=None
            if era5 == True:
                raster_nc['longitude'] = ((raster_nc['longitude'] + 180) % 360) - 180
                raster_nc = raster_nc.sortby("longitude")
                raster_nc = raster_nc.rename({'valid_time': 'time'})

        # Transform coordinate system if needed
        if transform:
            raster_nc = raster_nc.rio.reproject("EPSG:4326")

        df.apply(lambda row: extract_val(config,row,raster_nc,bandname,latname,lonname,col_name,ts), axis=1)

## Aggregate and merge

def agg_n_mer(config,df,imp_name,i,exp_name, daily=False, yearly = False):
    # raster extraction
    raster = pd.read_csv(config["basepath"] + config["output"] + imp_name + "_" + str(i) + "_Extraction.txt", sep = ";", header=None)
    raster.columns = ["ID","time","intervall","Value"]
    raster.time = pd.to_datetime(raster.time, format='ISO8601')

    # Create merge table
    raster["merge"] = None

    if yearly:
        raster["merge"] = raster["time"].dt.year.astype("int")
    elif daily == False:
        raster["merge"][(raster.intervall == "MS") | (raster.intervall == "d")] = pd.to_datetime(raster["time"][(raster.intervall == "MS") | (raster.intervall == "d")].dt.to_period('M').astype(str), format="%Y-%m")
        raster["merge"][raster.intervall == "YS"] = pd.to_datetime(raster["time"][raster.intervall == "YS"].dt.year.astype("int").astype(str), format="%Y")
    else:
        raster["merge"][raster.intervall == "d"] = pd.to_datetime(raster["time"][raster.intervall == "d"].dt.to_period('D').astype(str), format="%Y-%m-%d")
        raster["merge"][raster.intervall == "MS"] = pd.to_datetime(raster["time"][raster.intervall == "MS"].dt.to_period('M').astype(str),format="%Y-%m")
        raster["merge"][raster.intervall == "YS"] = pd.to_datetime(raster["time"][raster.intervall == "YS"].dt.year.astype("int").astype(str),format="%Y")

    raster_agg = raster.groupby(["merge","intervall","ID"]).mean()
    raster_agg.reset_index(inplace=True)
    raster_agg.drop(["intervall","time"],axis=1,inplace=True)

    if yearly:
        df_merged = pd.merge(df, raster_agg, how='outer', left_on=['year', 'ID'], right_on=['merge', 'ID'])
    elif daily:
        df_merged = pd.merge(df, raster_agg,  how='outer', left_on=['date','ID'], right_on = ['merge','ID'])
    else:
        df_merged = pd.merge(df, raster_agg, how='outer', left_on=['date2', 'ID'], right_on=['merge', 'ID'])

    df_merged.drop(['merge'],axis=1, inplace=True)
    df_merged.rename(columns={"Value":imp_name}, inplace=True)
    df_merged.to_csv(config["basepath"] + config["output"]  + exp_name +".txt", sep=";", index=False)
    return df_merged
