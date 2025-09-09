## func_processing_gw_time_series

'''This file contains all functions that are used in "01_processing_gw_time_series.py".'''

# import packages
import pandas as pd # imported version: 2.2.3
import numpy as np # imported version: 2.2.4
from scipy.spatial.distance import pdist, squareform # imported version: 1.15.2
from scipy.stats import spearmanr # imported version: 1.15.2
from datetime import datetime # internal package
from scipy import interpolate # imported version: 1.15.2
import operator # internal package
import pymannkendall as mk # imported version: 1.4.3


## extract longest sequence
'''Extracts longest continuous sequence of time series in which the time step intervals are below or
equal to a certain threshold.

df: dataframe with groundwater information
series: series with interval length between every time step in days
datatype: data type to which the series is transformed
relate: relative operator to be used for the comparison with the threshold
threshold: threshold value
append: if "True", the next time step of the longest continuous sequence is added to the sequence'''

def extract_seq(df,series,datatype,relate,threshold,append=False):
    rel_ops = {
        '>': operator.gt,
        '<': operator.lt,
        '>=': operator.ge,
        '<=': operator.le,
        '==': operator.eq,
        '!=': operator.ne}
    sequence = np.where(rel_ops[relate](series.astype(datatype), threshold).tolist())[0]  ## get every index without na at value
    longest_seq = max(np.split(sequence, np.where(np.diff(sequence) != 1)[0] + 1), key=len)  # find longest continous sequence in indexes

    if len(longest_seq) == 0:
        return(pd.DataFrame([]))
    else:
        if append:
            longest_seq = np.append(longest_seq, longest_seq[-1] + 1)  # to add last time step
        res = df.iloc[(longest_seq.min()):(longest_seq.max() + 1), :].reset_index(drop=True) # select longest continous sequence and reset index
        return(res)

## get_max_dist: Method and script adapted from Gunnar Lischeid
'''This function calculates the maximum allowed gap length between time steps 
based on autocorrelation. The allowed gap length is as long as the time series 
autocorrelates with a defined correlation coefficient threshold (Spearman r).
The maximum allowed gap length and a flag if the time series autocorrelates is returned.

df: dataframe with groundwater information
threshold: minimum correlation coefficient to label time series as autocorrelated
pairs: minimum number of data pairs (overlap) to perform correlation
pvalue: if True, significance of test is considered
'''

def get_max_dist(df,threshold, pairs, pvalue):
# check if aggregation is allowed
    # date as int
    data = df.copy()
    data["date_int"] = data["date"].apply(lambda x: ((datetime(1800, 1, 1) - x).total_seconds()))/86400
    # Calculate differences and lag classes
    dt = int(data.date_int.diff(periods=-1).median()) # dt is a day, month or year
    if dt == 0:
        count = False
        return(0,count)
    classes = np.arange(dt, 30*dt, dt) # maximum allowed gap is 30 dt steps (probably day, month or year)
    # Calculate correlation per lag class
    dist_matrix = squareform(pdist(data["date_int"].values.reshape(-1, 1), metric='euclidean')) # distance matrix
    dist_matrix[np.triu_indices_from(dist_matrix)] = None # Set upper triangle including diagonal to NA (None in Python)

    acf_table = pd.DataFrame(columns=['From', 'To', 'Autocorrelation',"p-value", 'Number_of_instances'], index=range(len(classes)-1))
    acf_table['From'] = classes[:-1]
    acf_table['To'] = classes[1:]

    # Loop through classes
    for i in range(len(classes) - 1):
    # Find indices where distance is within the current class range
        class_range = (acf_table.loc[i, 'From'] + acf_table.loc[i, 'To'])/2 # to sort the distances to the nearest class
        sel = np.where((dist_matrix >= (class_range-dt)) & (dist_matrix < class_range))

        # Calculate correlation
        if len(sel[0]) > pairs: # at least 30 pairs to calculate correlation
            correlation, p_val = spearmanr(data.loc[sel[0], "Value"], data.loc[sel[1], "Value"])
            acf_table.loc[i, 'Autocorrelation'] = correlation
            if pvalue:
               if p_val > 0.05:
                   acf_table.loc[i, 'Autocorrelation'] = None
            acf_table.loc[i, 'p-value'] = p_val
        else:
            acf_table.loc[i, 'Autocorrelation'] = None
            acf_table.loc[i, 'p-value'] = None

        # Count number of instances
        acf_table.loc[i, 'Number_of_instances'] = len(sel[0])

    # return maximum allowed distance
    subset = acf_table[(acf_table['Autocorrelation']>0).cummin()] # subset all positive correlation coefficients before negative ones
    if any(subset.Autocorrelation>=threshold):
        dists = ((subset.From[subset.Autocorrelation>=threshold]+subset.To[subset.Autocorrelation>=threshold])/2).to_list()[-1] # extract distance where the correlation coefficient is nearest above treshhold
        count = True
    else:
        dists = (acf_table.From[0]+acf_table.To[0])/2 # to include 31-day-months and leap years
        count = False
    return(dists, count)


## trim_max_dist
'''In case the allowed maximum gap length after the autocorrelation test is larger than the 
general gap length threshold, it is set to the threshold.

max_dist: maximum allowed gap length
threshold: general threshold gap length that bases on resolution and/or length of the time series
'''

def trim_max_dist(max_dist,threshold):
    if max_dist > threshold:
        max_dist = threshold
    return(max_dist)

## fill_gaps
'''Gaps in the groundwater time series are linearly interpolated. 
The dataframe with an additional column (filled time series) is returned.

df: dataframe with groundwater information
'''

def fill_gaps(df):
    dummy = df.copy()
    # construct time series with no gaps
    ts = pd.date_range(dummy.loc[0, 'date'], dummy.loc[len(dummy) - 1, 'date'],freq=dummy.interval[0])
    ts_int = (ts - datetime(1800, 1, 1)).total_seconds()
    dummy["dateint"] = dummy["date"].apply(lambda x: ((x - datetime(1800, 1, 1)).total_seconds()))
    gw_func = interpolate.interp1d(dummy["dateint"], dummy["Value"])

    ts_df = pd.DataFrame({"time_int": ts_int, "date": ts})
    raw3 = pd.merge(dummy, ts_df, "outer", left_on="dateint", right_on="time_int")
    raw3.reset_index(inplace=True,drop=True)
    raw3["groundwater_filled"] = gw_func(raw3["time_int"])
    raw3["date_x"] = raw3["date_y"]
    raw3.drop(axis=0, columns=["dateint","time_int","date_y"], inplace=True)
    raw3.rename(columns={"date_x":"date"},inplace=True)
    raw3["ID"] = dummy.ID[0]
    raw3["interval"] = dummy.interval[0]
    raw3["Parameter"] = dummy.Parameter[0]
    return(raw3)

## calc_trend
'''For the groundwater time series, the trend direction and Sen's slope are calculated and returned. 
The Mann-Kendall test is utilized. 

df: dataframe with groundwater information
pre: if True, groundwater time series autocorrelates, and the modified Mann Kendall test
by Hamed and Rao (1998) is used for the trend analysis
'''

def calc_trend(df,pre=False):
    if pre:
        if len(df)>2:
            res = mk.hamed_rao_modification_test(df.groundwater)
        else:
            return("no trend",None)
    else:
        res = mk.original_test(df.groundwater)

    trend = res.trend
    slope = res.slope

    # when there is no significant trend, the calculated slope is ignored; could be quite high but insignificant because of too less data points
    if trend == "no trend":
        return ("no trend", None)

    if df.parameter[0] != "Water level elevation a.m.s.l.":
        slope = slope * (-1)  # switch sign of slope
        if trend == "decreasing":
            trend = "increasing"
        elif trend == "increasing":
            trend = "decreasing"

    if np.unique(df.interval) == "d":
        slope = slope * 365
    elif np.unique(df.interval) == "MS":
        slope = slope * 12

    return (trend, slope)