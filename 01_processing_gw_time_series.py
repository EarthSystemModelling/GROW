"""Process groundwater time series data into a unified table.

This module reads individual Excel files containing groundwater time series,
performs quality control, and consolidates them into a single dataset.
Generates data flags and statistics for further analysis.
Discarded data is exported separately for manual quality control.
"""

import os  # built-in package
import shutil  # built-in package
import warnings  # built-in package
from datetime import datetime  # internal package

import numpy as np  # imported version: 2.2.4
import pandas as pd  # imported version: 2.2.3
from dateutil.relativedelta import relativedelta  # built-in package

from func_processing_gw_time_series import (
    calc_trend,  # calculates trend direction and slope
    extract_seq,  # extracts longest sequence in which the gaps do not exceed the individual gap length threshold
    fill_gaps,  # creates column where groundwater data is linearly gap-filled
    get_max_dist,  # derives maximum allowed gap length
    trim_max_dist,  # limits maximum allowed gap length to threshold
)

# Configuration: Path names, output names and other settings are defined here.
config = {
    "basepath" : "/mnt/storage/grow/01_Groundwater/", # GROW project directory for groundwater data
    "wells": "01_IGRAC_data_2025_08_18", # folder in which IGRAC's groundwater data is located
    "country_name_pos": 58, # position of first country letter in file path to extract country information
    # paths of exported output files
    "output": {"name":"_V08", # name of version
               "data":"02_Timeseries/wells_timeseries",
               "ts_attributes": "02_Timeseries/wells_timeseries_attributes",
               "max_dist": "02_Statistics/max_distance",
               "lost_<2_records": "02_Timeseries/lost_<2_records",
               "lost_mul_par": "02_Timeseries/lost_mul_par",
               "lost_unrealistic_value": "02_Timeseries/lost_unrealistic_value",
               "lost_after_aggregation": "02_Timeseries/lost_after_aggregation",
               "lost_gap_length": "02_Timeseries/lost_gap_length",
               "lost_gap_amount": "02_Timeseries/lost_gap_amount",
               "plateau": "02_Timeseries/wells_plateau",
               "jumps": "02_Timeseries/wells_jumps",
               "depth_all":"02_Timeseries/wells_depth_all_negative",
               "depth_any":"02_Timeseries/wells_depth_some_negative",
               "par": "02_Timeseries/wells_mul_par_all",
               "lost_per": "02_Statistics/wells_timeseries_drops",
               "duration": "GGMN_preprocessing_duration.txt"},
    "small": False # If True, only 1/50 of the data is processed to create a smaller test dataset
}

# Constants for temporal classification
DAILY_THRESHOLD_SECONDS = 172800  # 2 days
MONTHLY_THRESHOLD_SECONDS = 3542400  # 41 days
JUMP_THRESHOLD_METERS = 50
GAP_FRACTION_THRESHOLD = 0.2

warnings.filterwarnings("ignore")

# Store start time
startTime = datetime.now()

# Delete every folder in which there is no monitoring folder
folders = os.listdir(config["basepath"] + config["wells"])

full_path = []
for fold in folders:
    folder = config["basepath"] + config["wells"] + "/" + fold + "/" + "monitoring"
    full_path.append(folder)
    if not os.path.exists(
        folder
    ):  # if no monitoring folder exists in the path, delete path
        shutil.rmtree(config["basepath"] + config["wells"] + "/" + fold)

# Initiate variables
counter = 0  # to get the total amount of well time series

data_lost_a = []  # dropped because empty or only one record
data_lost_b = []  # dropped because too short after multiple parameter extraction
data_lost_c = []  # dropped because too short after unrealistic value removal
data_lost_d = []  # dropped because too short after aggregation
data_lost_e = []  # dropped because too short after max gap length extraction
data_lost_f = []  # dropped because too short after max gap fraction extraction

outliers_parameter = []  # all wells with multiple parameters
outliers_depth_any = []  # wells with some depth <= 0
outliers_depth_all = []  # wells with all depth <= 0
outliers_level = []  # wells with level == -999 and -9999
outliers_jumps = []  # wells with 1-2 50 m jumps
outliers_plateaus = []  # wells with plateaus

# Determine functions for aggregation
params = {
    "ID": "unique",
    "Parameter": "unique",
    "Value": "mean",
    "interval": "unique",
    "aggregated_from_n_values": "mean",
}

# loop over every country-folder
for path in full_path:
    # get every xlsx, even in subdirectories
    txt_files = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".ods"):
                txt_files.append(os.path.join(root, file))

    if config["small"]:
        txt_files = [txt_files[i] for i in range(0, len(txt_files), 50)]

    # loop over every well
    for file in txt_files:
        # Save counter
        counter = counter + 1
        pd.Series([str(datetime.now() - startTime), counter, file]).to_csv(
            config["basepath"] + "counter" + config["output"]["name"] + ".txt"
        )

        # Read groundwater time series
        raw = pd.read_excel(file, engine="odf", sheet_name=0, skiprows=[1, 1], dtype={"ID": object})
        raw.dropna(subset=["Value"],inplace = True,  ignore_index= True)

        # Some time series have duplicate records, lets get rid of them
        raw = raw.drop_duplicates().reset_index(drop=True)

        # For large data sets performance might be improved even more with
        # raw = raw.drop_duplicates(inplace=False).reset_index(drop=True)

        # Drop time series that only contain one record or are empty
        if raw.empty:
            data_lost_a = data_lost_a + [raw.copy(deep=True)]
            continue
        elif (
            len(raw) < 2
        ):  # does not work if empty, so there has to be a pre-check if it is empty
            data_lost_a = data_lost_a + [raw.copy(deep=True)]
            continue

        # Convert datetime column
        raw["Date and Time"] = pd.to_datetime(
            raw["Date and Time"], format="%Y-%m-%d %H:%M:%S UTC"
        )
        # flag year,month and day
        raw["year"] = raw["Date and Time"].dt.year.astype("int")
        raw["month"] = raw["Date and Time"].dt.to_period("M")
        raw["day"] = raw["Date and Time"].dt.to_period("D")

        # First clean-up/flagging

        ## Harmonize unit to meter
        raw.loc[raw["Unit"] == "ft", "Value"] = (
            raw.loc[raw["Unit"] == "ft", "Value"] * 0.3048
        )

        ## Check if there is more than one parameter in time series and extract only one (table depth)
        if raw.Parameter.nunique() > 1:
            outliers_parameter = outliers_parameter + [raw.copy(deep=True)]
            # we select water table level because this parameter contains more information than water depth
            raw = raw[raw.Parameter == "Water level elevation a.m.s.l."].reset_index(drop=True)
            # When time series is too short, discard current time series and continue with the next time series
            if len(raw)<2:
                data_lost_b = data_lost_b + [raw.copy(deep=True)]
                continue  ## discard this time series when there is only one record left now

        ## Unrealistic values in water table depth and level
        ### Remove negative groundwater depth or switch sign when all are negative (negative means above ground)
        if raw.loc[0, "Parameter"] != "Water level elevation a.m.s.l.":
            if (raw.Value < 0).all():  # all records are negative
                outliers_depth_all = outliers_depth_all + [raw.copy(deep=True)]
                raw.Value = raw.Value * (-1)  # when all values are negative, we just need a sign switch
            elif (raw.Value < 0).any():  # only some records are negative
                outliers_depth_any = outliers_depth_any + [raw]
                raw = raw.drop(index=raw[raw.Value < 0].index).reset_index(
                    drop=True
                )  # they are dropped
        ### Remove groundwater level records that are -999 or -9999
        else:
            if ((raw.Value == -999).any()) | ((raw.Value == -9999).any()):
                outliers_level = outliers_level + [raw.copy(deep=True)]
                raw = raw.drop(raw[(raw.Value == -999) | (raw.Value == -9999)].index).reset_index(drop=True)

        # When time series is too short, discard current time series and continue with the next time series
        if len(raw)<2:
            data_lost_c = data_lost_c + [raw.copy(deep=True)]
            continue

        # Aggregation

        ## Categorize in yearly, monthly or daily data
        raw["aggregated_from_n_values"] = None
        diff = abs(np.diff(raw["Date and Time"]).astype("timedelta64[s]"))

        if np.percentile(diff, 75) < DAILY_THRESHOLD_SECONDS:  # 2 days
            raw["interval"] = "d"
            raw["merge"] = pd.to_datetime(raw["day"].astype(str), format="%Y-%m-%d")
        elif np.percentile(diff, 75) < MONTHLY_THRESHOLD_SECONDS:  # 41 days
            raw["interval"] = "MS"
            raw["merge"] = pd.to_datetime(raw["month"].astype(str), format="%Y-%m")
        else:
            raw["interval"] = "YS"
            raw["merge"] = pd.to_datetime(raw["year"].astype(str), format="%Y")

        ## Add aggregation flag: Numbers of records per aggregated day, month or year
        ts = set(raw["merge"])
        if len(ts) == len(raw) & (raw["merge"] == raw["Date and Time"]).all():
            raw["aggregated_from_n_values"] = None
        else:
            for t in ts:
                raw["aggregated_from_n_values"][raw["merge"] == t] = len(
                    raw[raw["merge"] == t]
                )  # TODO

        ## Aggregate time series
        raw2 = (
            raw[["ID", "Parameter", "Value", "interval", "aggregated_from_n_values"]]
            .groupby(raw["merge"])
            .aggregate(params)
        )
        raw2.reset_index(inplace=True)
        raw2.rename(
            columns={raw2.columns[0]: "date"}, inplace=True
        )  # rename time column
        ## Aggregated unique parameters were saved as np.array, unpack them
        for col in raw2.columns:
            raw2[col] = raw2.apply(lambda x: pd.Series(x[col]), axis=1)

        ## When time series is too short, discard current time series and continue with the next time series
        if len(raw2) < 2:
            data_lost_d = data_lost_d + [raw2.copy(deep=True)]
            continue

        # Gaps

        ## Get max distance based on autocorrelation
        max_dist, bol = get_max_dist(
            raw2, 0.6, 30, False
        )  # when they are time steps with exact same time, 0 is returned

        ## Check time series for maximum allowed gap
        yn = abs(relativedelta(raw2["date"].iloc[-1], raw2["date"].iloc[0]).years) + 1

        if raw2.interval[0] == "d":
            max_dist = trim_max_dist(
                max_dist, 7
            )  # maximum allowed gap between daily time steps is 6 days --> allowed distance between existing time steps is 7
            threshold = 7  # save for plateaus later
        elif raw2.interval[0] == "MS":
            max_dist = trim_max_dist(
                max_dist, 123
            )  # a distance of 4 months (with possibly 3x 31-days-months) is allowed = 123 days
            threshold = 4
        else:
            if yn < 3:
                max_dist = trim_max_dist(
                    max_dist, 366
                )  # no gaps are allowed, 366 to include leap years
                threshold = 4  # no plateaus can be flagged
            elif (yn >= 3) & (yn < 10):
                max_dist = trim_max_dist(
                    max_dist, 731
                )  # a distance of 2 years (possibly including a leap year) is allowed
                threshold = 2
            elif (yn >= 10) & (yn < 15):
                max_dist = trim_max_dist(
                    max_dist, 1096
                )  # a distance of 3 years (possibly including a leap year) is allowed
                threshold = 3
            elif (yn >= 15) & (yn < 20):
                max_dist = trim_max_dist(
                    max_dist, 1461
                )  # a distance of 4 years (including a leap year) is allowed
                threshold = 4
            else:
                max_dist = trim_max_dist(
                    max_dist, 1827
                )  # a distance of 5 years (including two leap years) is allowed
                threshold = 5

        ## Extract longest sequence with intervals below max_dist
        diff = (
            abs(np.diff(raw2["date"]).astype("timedelta64[s]")).astype("int") / 86400
        )  # interval between every time step in days
        raw3 = extract_seq(raw2, diff, "int", "<=", max_dist, append=True)

        ## When time series is too short, discard current time series and continue with the next time series
        if len(raw3) < 2:
            data_lost_e = data_lost_e + [raw3.copy(deep=True)]
            continue

        ## Fill gaps
        raw4 = fill_gaps(raw3)

        ## Maximum gap length
        isna = raw4.Value.isna()  # get index where groundwater record is missing (NA)
        blocks = (
            ~isna
        ).cumsum()  # creates groups for consecutive values (here only True or False - groups)
        gaps = isna.groupby(blocks).sum()  # calculate length of continuous NA-sequences
        max_gap = gaps.max()
        mean_gap = gaps[gaps != 0].median()

        ## Drop timeseries with more than 20 % gaps
        gap_amount = len(raw4[raw4.Value.isna()]) / len(raw4)
        if gap_amount > GAP_FRACTION_THRESHOLD:
            raw4["year"] = raw4["date"].dt.year.astype("int")
            raw4_group = raw4.groupby("year")
            raw4_list = [raw4_group.get_group(x) for x in raw4_group.groups]
            gaps_year = [len(x[x.Value.isna()]) / len(x) for x in raw4_list]
            sequence = np.where((abs(np.array(gaps_year)) <= 0.2).tolist())[0]
            longest_seq = max(
                np.split(sequence, np.where(np.diff(sequence) != 1)[0] + 1), key=len
            )
            if len(longest_seq) == 0:
                data_lost_f = data_lost_f + [raw4.copy(deep=True)]
                continue
            years_list = list(np.unique(raw4.year))
            years_list = years_list[longest_seq.min() : longest_seq.max() + 1]
            raw5 = raw4[raw4.year.isin(years_list)].reset_index(drop=True)
            # When time series is too short, discard current time series and continue with the next time series
            if len(raw5) < 2:
                data_lost_f = data_lost_f + [raw5.copy(deep=True)]
                continue # data is discarded
            raw5.drop("year", axis=1, inplace=True)
            gap_amount = len(raw5[raw5.Value.isna()]) / len(raw5)
        else:
            raw5 = raw4

        # Flag jumps and spikes
        jum = pd.Series(
            np.diff(raw5.Value)
        )  # calculate water level changes between time steps
        jum_num = len(
            jum[abs(jum) >= JUMP_THRESHOLD_METERS]
        )  # counts number of water level changes between time steps that are higher than or equal to 50 m
        jumps = False
        if (
            jum_num > 0 & jum_num < 3
        ):  # if they are only one or two jumps, it is rather a problem with the device than the natural variance
            outliers_jumps = outliers_jumps + [raw5.copy(deep=True)]
            jumps = True  # Whole time series is flagged to contain suspicious jumps

        # Flag plateaus
        groups = (
            raw5["Value"] != raw5["Value"].shift()
        ).cumsum()  # Create groups for consecutive values
        raw5["plateaus"] = (
            raw5["Value"].groupby(groups).transform("size")
        )  # add length per group as column
        raw5.plateaus[raw5.plateaus < threshold] = (
            0  # turn every group under threshold into 0
        )
        if (raw5.plateaus > 0).any():
            outliers_plateaus = (outliers_plateaus +
                                 [raw5.copy(deep=True)])  # store all time series with plateaus

        # Add year, month and country again
        raw5["year"] = raw5["date"].dt.year.astype("int")
        raw5["month"] = raw5["date"].dt.to_period("M")
        raw5["country"] = file[
            config["country_name_pos"] : config["country_name_pos"] + 3
        ]

        # Clean-up columns
        raw5.rename(
            columns={"Parameter": "parameter", "Value": "groundwater"}, inplace=True
        )
        raw5 = raw5[
            [
                "ID",
                "country",
                "interval",
                "date",
                "year",
                "month",
                "aggregated_from_n_values",
                "plateaus",
                "parameter",
                "groundwater",
                "groundwater_filled",
            ]
        ]  # reorder

        # Derive time series attributes

        ## add trend; trend direction and slope sign are switched when the reference point is water table depth
        trend, slope = calc_trend(
            raw5, bol
        )  # If autocorrelation (bol=True) was detected, Ramed and Hao method is used

        ## number of occuring years
        yn = abs(relativedelta(raw5["date"].iloc[-1], raw5["date"].iloc[0]).years) + 1

        ## Flag for attributes, if time series contains any plateaus
        if raw5.plateaus.any() != 0:
            plat = True
        else:
            plat = False

        ## create attributes table
        attributes = pd.DataFrame(
            {
                "ID": [raw5.loc[0, "ID"]],
                "country": [raw5.loc[0, "country"]],
                "interval": [raw5.loc[0, "interval"]],
                "starting_date": [raw5.loc[0, "date"]],
                "ending_date": [raw5.loc[len(raw5) - 1, "date"]],
                "length_years": [yn],
                "reference_point": [raw5.loc[0, "parameter"]],
                "autocorrelation": bol,
                "aggregated_from_n_values_median": [
                    raw5.aggregated_from_n_values.median()
                ],
                "gap_fraction": [gap_amount],
                "jumps": jumps,
                "plateaus": [plat],
                "trend_direction": [trend],
                "trend_slope_m_year-1": [slope],
                "groundwater_mean_m": [raw5.groundwater.mean()],
                "groundwater_median_m": [raw5.groundwater.median()],
            }
        )
        ## create table with descriptive statistics about the gaps in the time series
        max_distance = pd.DataFrame(
            {
                "ID": [raw5.loc[0, "ID"]],
                "country": [raw5.loc[0, "country"]],
                "interval": [raw5.loc[0, "interval"]],
                "max_dist": [max_dist],
                "max_filled_gap": [max_gap],
                "mean_filled_gap": [mean_gap],
                "gap_count": [gap_amount],
            }
        )

        # append new data to files
        if counter == 1:
            # when this is the first time series
            raw5.to_csv(
                config["basepath"]
                + config["output"]["data"]
                + config["output"]["name"]
                + ".txt",
                index=False,
                sep=";",
            )
            attributes.to_csv(
                config["basepath"]
                + config["output"]["ts_attributes"]
                + config["output"]["name"]
                + ".txt",
                index=False,
                sep=";",
            )
            max_distance.to_csv(
                config["basepath"]
                + config["output"]["max_dist"]
                + config["output"]["name"]
                + ".txt",
                index=False,
                sep=";",
            )
        else:
            raw5.to_csv(
                config["basepath"]
                + config["output"]["data"]
                + config["output"]["name"]
                + ".txt",
                index=False,
                sep=";",
                mode="a",
                header=False,
            )
            attributes.to_csv(
                config["basepath"]
                + config["output"]["ts_attributes"]
                + config["output"]["name"]
                + ".txt",
                index=False,
                sep=";",
                mode="a",
                header=False,
            )
            max_distance.to_csv(
                config["basepath"]
                + config["output"]["max_dist"]
                + config["output"]["name"]
                + ".txt",
                index=False,
                sep=";",
                mode="a",
                header=False,
            )

# export fractions of time series that were dropped or that contain anomalies (outliers)
drops = pd.DataFrame(
    {
        "name": [
            "no_or_one_record",
            "Lost_Multiple_parameters",
            "Unrealistic_Value",
            "After_aggregation",
            "max_gap_length",
            "max_gap_fraction",
            "outliers_multiple_parameters",
            "outliers_depth_any",
            "outliers_depth_all",
            "outliers_level",
            "outliers_jumps",
            "outliers_plateaus",
            "all",
        ],
        "percentage": [
            len(data_lost_a) / counter,
            len(data_lost_b) / counter,
            len(data_lost_c) / counter,
            len(data_lost_d) / counter,
            len(data_lost_e) / counter,
            len(data_lost_f) / counter,
            len(outliers_parameter) / counter,
            len(outliers_depth_any) / counter,
            len(outliers_depth_all) / counter,
            len(outliers_level) / counter,
            len(outliers_jumps) / counter,
            len(outliers_plateaus) / counter,
            counter / counter,
        ],
        "number": [
            len(data_lost_a),
            len(data_lost_b),
            len(data_lost_c),
            len(data_lost_d),
            len(data_lost_e),
            len(data_lost_f),
            len(outliers_parameter),
            len(outliers_depth_any),
            len(outliers_depth_all),
            len(outliers_level),
            len(outliers_jumps),
            len(outliers_plateaus),
            counter,
        ],
    }
)
drops.to_csv(
    config["basepath"]
    + config["output"]["lost_per"]
    + config["output"]["name"]
    + ".txt",
    index=False,
    sep=";",
)

# export time series that were dropped and time series with anomalies (outliers)
lists = [data_lost_a,data_lost_b,data_lost_c,data_lost_d,data_lost_e,
         data_lost_f,outliers_plateaus,outliers_jumps,outliers_depth_all,
         outliers_depth_any,outliers_parameter]
paths = ["lost_<2_records","lost_mul_par","lost_unrealistic_value","lost_after_aggregation","lost_gap_length",
         "lost_gap_amount","plateau","jumps","depth_all",
         "depth_any","par"]

for li,pa in zip(lists,paths):
    try:
        pd.concat(li).to_csv(
            config["basepath"]
            + config["output"][pa]
            + config["output"]["name"]
            + ".txt",
            index=False,
            sep=";",
        )
    except ValueError: # in case the list is empty
        pass

# Print duration of script
np.savetxt(
    config["basepath"]
    + config["output"]["duration"]
    + config["output"]["name"]
    + ".txt",
    np.array([str(datetime.now() - startTime)]),
    fmt="%s",
)  # print time it takes for the whole script to run through
