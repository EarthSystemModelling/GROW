### usage_example

'''This is an example showing how GROW can be easily subset for specific needs.
In this example, we might want to investigate decreasing groundwater tables.'''

# import packages
import pandas as pd # imported version: 2.2.3

# load GROW tables
attributes = pd.read_csv("/mnt/storage/grow/final_grow/grow_attributes.csv", sep=";")
timeseries = pd.read_csv("/mnt/storage/grow/final_grow/grow_timeseries_V05.txt", sep=";")

# Subset attribute table
# trend_direction is "decreasing"
# at least 20 years long
# monthly resolution (daily or monthly possible, daily is later aggregated)
attr_decreasing = attributes[(attributes.trend_direction=="decreasing") & (attributes.length_years>=20) & (attributes.interval!="YS")]

# Subset time series table by GROW_ID
ts_decreasing = timeseries[timeseries.GROW_ID.isin(attr_decreasing.GROW_ID)]

# Drop filter columns that are not needed
ts_decreasing.drop(columns=["date","year","aggregated_from_n_values","plateaus"], inplace=True)

# Aggregate time series to a monthly resolution
ts_decreasing_monthly = ts_decreasing.groupby(["GROW_ID","month"],as_index=False).mean() # Aggregate all numeric variables per ID and month
'''All quantity units are given in mm/year. The unit can be adjusted to the daily or monthly scale by the division of 365 or 12.'''
ts_decreasing_monthly.iloc[:,5:10] = ts_decreasing_monthly.iloc[:,5:10].apply(lambda x: x/12, axis=1) # adjust quantity units from mm/year to mm/month by dividing with 12

'''Now, you are ready to start your analysis :)'''

print(len(attr_decreasing[attr_decreasing.main_landuse=="cropland_irrigated"])/len(attr_decreasing))
print(len(attributes[attributes.main_landuse=="cropland_irrigated"])/len(attributes))


