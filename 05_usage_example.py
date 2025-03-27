### Usage example 1

'''This is an example showing how GROW can be easily subset for specific needs. In this example, we might want to investigate
the seasonal pattern of groundwater level observations in Brazil. For that, we need to 1. extract monthly groundwater time series from Brazil,
2. aggregate the time series and 3. adjust the units.'''

import pandas as pd
import scipy
import matplotlib.pyplot as plt

# load GROW tables
attributes = pd.read_csv("/mnt/storage/grow/final_grow/grow_attributes_to_show.txt", sep=";")
timeseries = pd.read_csv("/mnt/storage/grow/final_grow/grow_timeseries_to_show.txt", sep=";")

# Subset attribute table
attr_bra_monthly = attributes[(attributes.country=="BRA") & (attributes.interval!="YS")]

# Subset time series table
ts_bra = timeseries[timeseries.ID.isin(attr_bra_monthly.ID)]

# Drop filter columns that are not needed
ts_bra.drop(columns=["date","methodology","aggregated_from_n_values","plateaus"], inplace=True)

# Aggregate time series to a monthly resolution
ts_bra_monthly = ts_bra.groupby(["ID","month"],as_index=False).mean() # Aggregate all numeric variables per ID and month
'''All quantity units are given in mm/year. The unit can be adjusted to the daily or monthly scale by the division of 365 or 12.'''
ts_bra_monthly.iloc[:,5:10] = ts_bra_monthly.iloc[:,5:10].apply(lambda x: x/12, axis=1) # adjust quantity units from mm/year to mm/month by dividing with 12

'''Now, you are ready to start your analysis :)'''

# Why wait, let's see if there is a correlation between groundwater and precipitation
no_nan = ts_bra_monthly.dropna(subset=["groundwater_depth_m","precipitation_gpcc_mm_year-1","precipitation_mswep_mm_year-1"],how="any")

print(scipy.stats.spearmanr(no_nan["groundwater_depth_m"],no_nan["precipitation_mswep_mm_year-1"]))
print(scipy.stats.spearmanr(no_nan["groundwater_depth_m"],no_nan["precipitation_gpcc_mm_year-1"]))

fig,ax = plt.subplots()
plt.scatter(no_nan["groundwater_depth_m"],no_nan["precipitation_gpcc_mm_year-1"])
ax.set_xlabel("groundwater depth in m")
ax.set_ylabel("Monthly precipitation in mm")
plt.show()