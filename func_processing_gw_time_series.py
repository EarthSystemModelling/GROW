"""This file contains all functions that are used in "01_processing_gw_time_series.py"."""

# import packages
import pandas as pd  # imported version: 2.2.3
import numpy as np  # imported version: 2.2.4
from scipy.spatial.distance import pdist, squareform  # imported version: 1.16.1
from scipy.stats import spearmanr  # imported version: 1.16.1
from datetime import datetime  # internal package
from scipy import interpolate  # imported version: 1.16.1
import operator  # internal package
import pymannkendall as mk  # imported version: 1.4.3
from sklearn.cluster import DBSCAN # scikit-learn 1.7.2 need to be installed
from sklearn.utils._param_validation import InvalidParameterError # scikit-learn 1.7.2 need to be installed

def extract_seq(
        df: pd.DataFrame,
        series: pd.Series,
        datatype: str,
        relate: str,
        threshold: int,
        append=False,):
    """Extract longest continuous sequence in time series below interval threshold.

    Args:
        df: dataframe with groundwater information
        series: series with interval length between time step in days
        datatype: data type for series conversion
        relate: comparison operator as string ('<=', '>=', etc.)
        threshold: threshold value for comparison
        append: whether to append next time step to sequence

    Returns:
        DataFrame containing the longest valid sequence
    """
    rel_ops = {
        ">": operator.gt,
        "<": operator.lt,
        ">=": operator.ge,
        "<=": operator.le,
        "==": operator.eq,
        "!=": operator.ne,
    }
    ## get index where the interval to the next time step is <= allowed distance
    sequence = np.where(rel_ops[relate](series.astype(datatype), threshold).tolist())[0]
    # find longest continuous sequence in indexes
    longest_seq = max(np.split(sequence, np.where(np.diff(sequence) != 1)[0] + 1), key=len)

    if len(longest_seq) == 0:
        return pd.DataFrame([])
    else:
        if append:
            longest_seq = np.append(longest_seq, longest_seq[-1] + 1)  # to add last record
        # select longest continuous sequence from time series and reset index
        res = df.iloc[(longest_seq.min()) : (longest_seq.max() + 1), :].reset_index(drop=True)
        return res


def get_max_dist(df, threshold, pairs, pvalue):
    """This function calculates the maximum allowed distance between two time steps between time steps
    based on autocorrelation. The allowed distance between two time steps is as long as the time series
    autocorrelates with a defined correlation coefficient threshold (Spearman r).
    The maximum allowed distance between two time steps and a flag if the time series autocorrelates is returned.

    Args:
        df: dataframe with groundwater information
        threshold: minimum correlation coefficient to label time series as autocorrelated
        pairs: minimum number of data pairs (overlap) to perform correlation
        pvalue: if True, significance of test is considered

    Returns:
        Maximum allowed distance between two time steps
        Flag if time series autocorrelates [True/False]
    """
    # date as int
    data = df.copy().reset_index(drop=True)
    data["day"] = pd.to_datetime(data["date"].dt.to_period("D").astype(str))
    data["date_int"] = (
        data["day"].apply(lambda x: ((datetime(1800, 1, 1) - x).total_seconds()))/ 86400
    )
    # set lag distance as median of the time step intervals in the time series
    dt = int(data.date_int.diff(periods=-1).median())  # dt is a day, month or year
    if dt == 0:
        count = False
        return (0, count)
    # generate lag classes
    # a maximum of 30 lags is tested
    classes = np.arange(dt, 30 * dt, dt)
    # calculate correlation per lag class
    # create distance matrix
    dist_matrix = squareform(pdist(data["date_int"].values.reshape(-1, 1), metric="euclidean"))
    # set upper triangle including diagonal to NA (None in Python)
    dist_matrix[np.triu_indices_from(dist_matrix)] = None

    # create table that contains the correlation results per class
    acf_table = pd.DataFrame(
        columns=["From", "To", "Autocorrelation", "p-value", "Number_of_instances"],
        index=range(len(classes) - 1),
    )
    acf_table["From"] = classes[:-1]
    acf_table["To"] = classes[1:]

    # Loop through classes
    for i in range(len(classes) - 1):
        # Find instances where time step distances are within the current class range
        # sort the distances to the nearest class
        class_range = (acf_table.loc[i, "From"] + acf_table.loc[i, "To"]) / 2
        # find index pairs per class
        sel = np.where((dist_matrix >= (class_range - dt)) & (dist_matrix < class_range))

        # Calculate correlation
        if len(sel[0]) >= pairs:  # at least 30 pairs to calculate correlation
            correlation, p_val = spearmanr(data.loc[sel[0], "groundwater"], data.loc[sel[1], "groundwater"])
            acf_table.loc[i, "Autocorrelation"] = correlation
            if pvalue:
                if p_val > 0.05: # not significant
                    acf_table.loc[i, "Autocorrelation"] = None
            acf_table.loc[i, "p-value"] = p_val
        else:
            acf_table.loc[i, "Autocorrelation"] = None
            acf_table.loc[i, "p-value"] = None

        # Count number of instances
        acf_table.loc[i, "Number_of_instances"] = len(sel[0])

    # return maximum allowed distance
    # subset all positive correlation coefficients before negative ones
    subset = acf_table[(acf_table["Autocorrelation"] > 0).cummin()]
    if any(subset.Autocorrelation >= threshold):
        # extract maximum distance where the correlation coefficient is still >= threshold
        dists = ((subset.From[subset.Autocorrelation >= threshold]+subset.To[subset.Autocorrelation >= threshold])/2).to_list()[-1]
        count = True # time series autocorrelates
    else:
        dists = (acf_table.From[0] + acf_table.To[0]) / 2  # to include 31-day-months and leap years
        count = False # time series autocorrelates
    return (dists, count)

def trim_max_dist(max_dist, threshold):
    """In case the allowed maximum gap length after the autocorrelation test is larger than the
    general threshold, it is set to the threshold.

    Args:
        max_dist: maximum allowed distance between two time steps before trimming
        threshold: general threshold for maximum allowed distance that bases on resolution and
        /or length of the time series

    Returns:
        maybe trimmed maximum allowed distance between two time steps
    """
    if max_dist > threshold:
        max_dist = threshold
    return max_dist

def fill_gaps(df):
    """Gaps in the groundwater time series are linearly interpolated.
    The dataframe with an additional column (filled time series) is returned.

    Args:
        df: dataframe with groundwater time series

    Returns:
        Dataframe with additional gap-filled groundwater column
    """
    dummy = df.copy()
    # construct time series with no gaps
    # complete time series
    ts = pd.date_range(dummy.loc[0, "date"], dummy.loc[len(dummy) - 1, "date"], freq=dummy.interval[0])
    # to interpolate the gaps we need a function between the groundwater values and the date
    # the date needs to be numeric for that purpose
    ts_int = (ts - datetime(1800, 1, 1)).total_seconds()
    ts_df = pd.DataFrame({"time_int": ts_int, "date": ts})

    # get interpolation function
    dummy["dateint"] = dummy["date"].apply(lambda x: ((x - datetime(1800, 1, 1)).total_seconds()))
    gw_func = interpolate.interp1d(dummy["dateint"], dummy["groundwater"])

    # merge complete time series with groundwater time series to get a complete date column
    raw3 = pd.merge(dummy, ts_df, "outer", left_on="dateint", right_on="time_int")
    raw3.reset_index(inplace=True, drop=True)

    # use interpolation function to construct a complete groundwater time series
    raw3["groundwater_filled"] = gw_func(raw3["time_int"])
    # I actually want to have date_y as new complete datetime column but to avoid changing column order,
    # I just overwrite the old datetime column
    raw3["date_x"] = raw3["date_y"]
    raw3.drop(axis=0, columns=["dateint", "time_int", "date_y"], inplace=True)
    raw3.rename(columns={"date_x": "date"}, inplace=True)
    # refill unique-value columns to avoid NA in them
    raw3["ID"] = dummy.ID[0]
    raw3["interval"] = dummy.interval[0]
    raw3["parameter"] = dummy.parameter[0]
    return raw3

def flag_negative_signs(df,list_all,list_any):
    sign_cat = None
    if df.loc[0, "parameter"] != "Water level elevation a.m.s.l.":
        dummy = df.groundwater.dropna() # when series contains NA, all() will be skipped
        if (dummy < 0).all():  # all records are negative
            list_all = list_all + [df.copy(deep=True)]
            sign_cat = "All"
        elif (dummy < 0).any():  # only some records are negative
            list_any = list_any + [df.copy(deep=True)]
            sign_cat = "Some"
        else:
            sign_cat = "No"
    return sign_cat, list_all, list_any

def flag_jumps(df,jump_threshold,list_jumps):
    jum = pd.Series(
        np.diff(df.groundwater)
    )  # calculate water level changes between time steps
    jum_num = len(
        jum[abs(jum) >= jump_threshold]
    )  # counts number of water level changes between time steps that are higher than or equal to 50 m
    if (jum_num > 0):  # if there are jumps >50 m list table
        list_jumps = list_jumps + [df.copy(deep=True)]
    return jum_num, list_jumps

def flag_outliers(df,list_outliers):
    ## Outliers detected with DBSCAN
    # Normalize time series
    copy = df.copy(deep=True).reset_index(drop=True)
    # we need to use abs because there are negative and positive records in the dataset
    copy_mean = copy.groundwater.mean()
    copy_std = copy.groundwater.std()
    copy.groundwater = (copy.groundwater - copy_mean) / copy_std
    copy.dropna(subset="groundwater", inplace=True)
    # Document original ts length
    length_start = len(copy)
    # Generate a model
    try:
        model = DBSCAN(min_samples=round(len(copy) * 0.5), eps=2).fit(copy[['groundwater']])
        colors = model.labels_
        colors = abs(colors)
        c = np.bincount(colors).argmax()
        copy['cluster'] = colors
        # To mark outliers in the time series, the cluster is merged to the time series
        df = pd.merge(df,copy[["ID","date","cluster"]],on=["ID","date"] ,how="left")
        # Records which are not part of the main cluster and not NA (gap) are marked as outlier
        df["outliers"] = False
        df.outliers[(df.cluster != c) & (df.cluster.notna())] = True
        df.drop(columns=["cluster"], inplace=True)
        copy = copy.loc[copy['cluster'] == c, :]
        copy['groundwater'] = copy['groundwater'] * copy_std + copy_mean  # de-normalize
        length_end = len(copy)
        if length_end < length_start:
            outlier = "True"
            list_outliers = list_outliers + [df.copy(deep=True)]
        else:
            outlier = "False"
    except InvalidParameterError:
        outlier = "False"
        df["outliers"] = False

    return df, outlier, list_outliers

def flag_plateaus(df,threshold,list_plat):
    # Create groups for consecutive values
    groups = (df["groundwater"] != df["groundwater"].shift()).cumsum()
    # add length per group as column
    df["plateaus"] = df["groundwater"].groupby(groups).transform("size")
    # turn every group under threshold into 0
    df.plateaus[df.plateaus < threshold] = 0
    if (df.plateaus > 0).any():
        plat=True
        list_plat = list_plat + [df.copy(deep=True)]  # store all time series with plateaus
    else:
        plat=False
    return plat, list_plat

def calc_trend(df, slop=True, pre=False):
    """For the groundwater time series, the trend direction and Sen's slope are calculated and returned.
    The Mann-Kendall test is utilized.

    df: dataframe with groundwater information
    pre: if True, groundwater time series autocorrelates, and the modified Mann Kendall test
    by Hamed and Rao (1998) is used for the trend analysis
    """
    if pre:
        if len(df) > 2:
            res = mk.hamed_rao_modification_test(df.groundwater)
        else:
            return ("no trend", None)
    else:
        res = mk.original_test(df.groundwater)

    trend = res.trend

    # when there is no significant trend, the calculated slope is ignored;
    # could be quite high but insignificant because of too less data points
    if trend == "no trend":
        return ("no trend", None)

    slope = res.slope
    if df.parameter[0] != "Water level elevation a.m.s.l.":
        slope = slope * (-1)  # switch sign of slope
        if trend == "decreasing":
            trend = "increasing"
        elif trend == "increasing":
            trend = "decreasing"

    if slop:
        if np.unique(df.interval) == "d":
            slope = slope * 365
        elif np.unique(df.interval) == "MS":
            slope = slope * 12

    return (trend, slope)
