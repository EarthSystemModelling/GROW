"""This example shows how GROW can be easily subset for specific needs
even on machines that do not have enough RAM to load the whole dataset at once.

Let's assume we want to investigate the seasonal pattern of groundwater observations in Brazil.
For that, we need to
  1.) extract monthly groundwater time series from Brazil,
  2.) aggregate the time series and
  3.) adjust the units.
"""

import pandas as pd  # imported version: 2.2.3
import pyarrow.dataset as ds # imported version: 19.0.1

path_attributes = "/mnt/storage/grow/09_final_grow/grow_attributes.csv"
path_timeseries = "/mnt/storage/grow/09_final_grow/grow_timeseries.parquet"

# load GROW attributes tables
print("loading attributes...")
attributes = pd.read_csv(path_attributes, sep=",")

# 1.) select GROW_IDs based on filters
# only time series in Brazil ("BRA")
# only time series which are monthly or daily (not yearly; YS)
print("filtering attributes...")
attr_bra_monthly = attributes[(attributes.country == "BRA") & (attributes.interval != "YS")]

# load only needed time series based on the filtered attributes GROW-IDs
# Option 1: Using pandas
print("loading timeseries.parquet with filter...")
ts_bra = pd.read_parquet(
    path_timeseries,
    filters=[("GROW_ID", "in", attr_bra_monthly["GROW_ID"])],
)

# Option 2: Using pyarrow
# # load only needed timeseries without full memory load
print("loading timeseries.parquet with filter...")
timeseries_dataset_raw = ds.dataset(path_timeseries, format="parquet")
## create table based on attributes subset
ts_bra_raw= timeseries_dataset_raw.to_table(filter=ds.field("GROW_ID").isin(attr_bra_monthly.GROW_ID))
ts_bra = ts_bra_raw.to_pandas() # convert table object into pandas DataFrame

# 2.) aggregate time series to a monthly resolution
# drop filter columns that are not needed anymore
print("aggregating time series...")
ts_bra.drop(
    columns=["date", "interval", "year", "aggregated_from_n_values", "outliers_change_points", "plateaus"],
    inplace=True,
)
# aggregate all numeric variables per ID and month
ts_bra_monthly = ts_bra.groupby(["GROW_ID", "month"], as_index=False).mean()

print("adjust units in time series to monthly resolution...")
# 3.) adjust units
"""All quantity units are given in mm/year. The unit can be adjusted to the daily or monthly scale by the division of 365 or 12. 
All variables with quantity units are arranged in a row (here column 8 to 14)"""
# adjust quantity units from mm/year to mm/month by dividing with 12
ts_bra_monthly.iloc[:, 8:14] = ts_bra_monthly.iloc[:, 8:14].apply(lambda x: x / 12, axis=1)
# convert month-column to datetime
ts_bra_monthly["month"] = pd.to_datetime(ts_bra_monthly["month"], format="%Y-%m")

"""Now, you are ready to start your analysis :)"""
print(f"filtered {len(attr_bra_monthly)} timeseries with a total of {len(ts_bra_monthly)} data points")