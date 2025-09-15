"""This example shows how to apply filters to query GROW even on machines that do not have enough RAM to load the whole dataset at once.
Let's assume we want to investigate the seasonal pattern of groundwater observations in Brazil.
For that, we need to
  1.) extract monthly groundwater time series from Brazil,
  2.) aggregate the time series and
  3.) adjust the units.
"""

import pandas as pd  # imported version: 2.2.3

# load GROW tables
print("loading csv...")
attributes = pd.read_csv("./data/grow_attributes.csv", sep=";")

print("filtering...")
# 1.) Select GROW_IDs based on filters
# only time series in Brazil ("BRA")
# only time series which are monthly or daily (not yearly; YS)
attr_bra_monthly = attributes[
    (attributes.country == "BRA") & (attributes.interval != "YS")
]


print("loading parquet with filter...")
timeseries = pd.read_parquet(
    "./data/grow_timeseries.parquet",
    filters=[("GROW_ID", "in", attr_bra_monthly["GROW_ID"])],
)

# 2.) aggregate time series to a monthly resolution
# drop filter columns that are not needed
timeseries.drop(
    columns=["date", "interval", "year", "aggregated_from_n_values", "plateaus"],
    inplace=True,
)
# aggregate all numeric variables per ID and month
ts_bra_monthly = timeseries.groupby(["GROW_ID", "month"], as_index=False).mean()

print("processing...")
# 3.) adjust units
"""All quantity units are given in mm/year. The unit can be adjusted to the daily or monthly scale by the division of 365 or 12. 
All variables with quantity units are arranged in a row (here column 8 to 13)"""
ts_bra_monthly.iloc[:, 8:14] = ts_bra_monthly.iloc[:, 8:14].apply(
    lambda x: x / 12, axis=1
)  # adjust quantity units from mm/year to mm/month by dividing with 12
ts_bra_monthly["month"] = pd.to_datetime(
    ts_bra_monthly["month"], format="%Y-%m"
)  # convert month-column to datetime

"""Now, you are ready to start your analysis :)"""
print(
    f"filtered {ts_bra_monthly['GROW_ID'].size} timeseries with a total of {ts_bra_monthly.size} data points"
)
print(ts_bra_monthly.head())