### usage_example
#
"""This is an example showing how GROW can be easily subset for specific needs. In this example, we might want to investigate
the seasonal pattern of groundwater observations in Brazil. For that, we need to 1.) extract monthly groundwater time series from Brazil,
2.) aggregate the time series and 3.) adjust the units."""

import pandas as pd # imported version: 2.2.3
import pyarrow.dataset as ds # imported version: 19.0.1

path_attributes = "/mnt/storage/grow/09_final_grow/grow_attributes.parquet" # table in parquet-format
path_timeseries = "/mnt/storage/grow/09_final_grow/grow_timeseries.parquet" # table in parquet-format

# load GROW tables
attributes = pd.read_parquet(path_attributes) # load full attribute table

# 1.) subset attribute and time series table
# only time series in Brazil ("BRA")
# only time series which are monthly or daily (not yearly; YS)
attr_bra_monthly = attributes[(attributes.country=="BRA") & (attributes.interval!="YS")] # subset attribute table by data flag columns

# subset time series table BY GROW_ID and load only needed data fields
timeseries_dataset_raw = ds.dataset(path_timeseries, format="parquet") # load table without full memory usage
ts_bra_raw= timeseries_dataset_raw.to_table(filter=ds.field("GROW_ID").isin(attr_bra_monthly.GROW_ID)) # create table based on attribute subset
ts_bra = ts_bra_raw.to_pandas() # convert table object to pandas DataFrame

# 2.) aggregate time series to a monthly resolution
ts_bra.drop(columns=["date","interval","year","aggregated_from_n_values","plateaus"], inplace=True) # drop filter columns that are not needed
ts_bra_monthly = ts_bra.groupby(["GROW_ID","month"],as_index=False).mean() # aggregate all numeric variables per ID and month

# 3.) adjust units
'''All quantity units are given in mm/year. The unit can be adjusted to the daily or monthly scale by the division of 365 or 12. 
All variables with quantity units are arranged in a row (here column 8 to 13)'''
ts_bra_monthly.iloc[:,8:14] = ts_bra_monthly.iloc[:,8:14].apply(lambda x: x/12, axis=1) # adjust quantity units from mm/year to mm/month by dividing with 12
ts_bra_monthly["month"] = pd.to_datetime(ts_bra_monthly["month"],format="%Y-%m") # convert month-column to datetime

'''Now, you are ready to start your analysis :)'''