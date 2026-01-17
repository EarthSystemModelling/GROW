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
import matplotlib.pyplot as plt


from func_processing_gw_time_series import (
    calc_trend,  # calculates trend direction and slope
    extract_seq,  # extracts longest sequence in which the gaps do not exceed the individual gap length threshold
    fill_gaps,  # creates column where groundwater data is linearly gap-filled
    get_max_dist,  # derives maximum allowed gap length
    trim_max_dist, # limits maximum allowed gap length to threshold
    flag_negative_signs, # flags if a time series conatins negative signs
    flag_outliers, # flags if a time series contains outliers or change points
    flag_plateaus # flags if time series contains plateaus
)

# Configuration: Path names, output names and other settings are defined here.
config = {
    "basepath" : "/mnt/storage/grow/01_Groundwater/", # GROW project directory for groundwater data
    "wells": "01_IGRAC_data_2025_08_18", # folder in which IGRAC's groundwater data is located
    "country_name_pos": 58, # position of first country letter in file path to extract country information
    # paths of exported output files
    "output": {"name":"_test", # name of version
               "data":"02_Timeseries/sensitivity_analysis/gaps_10/wells_timeseries",
               "ts_attributes": "02_Timeseries/sensitivity_analysis/gaps_10/wells_timeseries_attributes",
               "max_dist": "02_Statistics/max_distance",
               "lost_<2_records": "02_Timeseries/sensitivity_analysis/gaps_10/lost_<2_records",
               "lost_mul_par": "02_Timeseries/sensitivity_analysis/gaps_10/lost_mul_par",
               "lost_unrealistic_value": "02_Timeseries/sensitivity_analysis/gaps_10/lost_unrealistic_value",
               "lost_after_aggregation": "02_Timeseries/sensitivity_analysis/gaps_10/lost_after_aggregation",
               "lost_gap_length": "02_Timeseries/sensitivity_analysis/gaps_10/lost_gap_length",
               "lost_gap_amount": "02_Timeseries/sensitivity_analysis/gaps_10/lost_gap_amount",
               "plateau": "02_Timeseries/sensitivity_analysis/gaps_10/wells_plateau",
               "outliers": "02_Timeseries/sensitivity_analysis/gaps_10/wells_outliers",
               "depth_all":"02_Timeseries/sensitivity_analysis/gaps_10/wells_depth_all_negative",
               "depth_any":"02_Timeseries/sensitivity_analysis/gaps_10/wells_depth_some_negative",
               "par": "02_Timeseries/sensitivity_analysis/gaps_10/wells_mul_par_all",
               "more_id": "02_Timeseries/sensitivity_analysis/gaps_10/more_than_one_id",
               "lost_per": "02_Statistics/wells_timeseries_drops",
               "duration": "GGMN_preprocessing_duration.txt"},
    "small": False # If True, only 1/50 of the data is processed to create a smaller test dataset
}

# Constants for temporal classification
PERCENTILE_THRESHOLD = 90 # percentile for aggregation
DAILY_THRESHOLD_SECONDS = 172800  # 2 days
MONTHLY_THRESHOLD_SECONDS = 3542400  # 41 days
MAX_DISTANCE_DAYS = 4 # 3 day gap allowed (= 4 days distance)
MAX_DISTANCE_MONTHS = 62 # 1 month gap allowed (= 2 months distance)
MAX_DISTANCE_1_YEARS = 366 # no gaps allowed (= 1 year distance)
MAX_DISTANCE_5_YEARS = 366 # no gaps allowed (= 1 year distance)
MAX_DISTANCE_10_YEARS = 731 # 1 year gap allowed (= 2 year distance)
MAX_DISTANCE_15_YEARS = 731 # 1 year gap allowed (= 2 year distance)
MAX_DISTANCE_20_YEARS = 1096 # 2 year gap allowed (= 3 year distance)
PLAT_THRESHOLD_DAILY_DAYS = 4 # plateaus are flagged starting from 4 days
PLAT_THRESHOLD_MONTHLY_MONTHS = 2 # plateaus are flagged starting from 2 months
PLAT_THRESHOLD_YEARLY_1_YEARS = 2 # plateaus are flagged starting from 2 years
PLAT_THRESHOLD_YEARLY_5_YEARS = 2 # plateaus are flagged starting from 2 years
PLAT_THRESHOLD_YEARLY_10_YEARS = 2 # plateaus are flagged starting from 2 years
PLAT_THRESHOLD_YEARLY_15_YEARS = 2 # plateaus are flagged starting from 2 years
PLAT_THRESHOLD_YEARLY_20_YEARS = 3 # plateaus are flagged starting from 3 years
GAP_FRACTION_THRESHOLD = 0.1 # Total allowed gap fraction

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
outliers_novalue = []  # wells with -999/-9999 records
outliers_dbscan = []  # wells with 1-2 50 m outliers
outliers_plateaus = []  # wells with plateaus
more_than_one_id = [] # wells with more than one ID

# Determine functions for aggregation
params = {
    "ID": "unique",
    "parameter": "unique",
    "groundwater": "mean",
    "interval": "unique",
    "aggregated_from_n_values": "mean",
}

# loop over every country-folder
for path in full_path:
    # get every ods, even in subdirectories
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
        ## rename columns
        raw.rename(columns={"Date and Time": "date",
                            "Parameter": "parameter",
                            "Value": "groundwater"}, inplace=True)
        ## Drop empty records
        raw.dropna(subset=["groundwater"],inplace = True,  ignore_index= True)

        # Some time series have duplicate records, lets get rid of them
        raw = raw.drop_duplicates(inplace=False).reset_index(drop=True)

        # Drop time series that only contain one record or are empty
        if raw.empty:
            data_lost_a = data_lost_a + [raw.copy(deep=True)]
            continue
        elif len(raw)<2: # does not work if empty, so there has to be a pre-check if it is empty
            data_lost_a = data_lost_a + [raw.copy(deep=True)]
            continue

        if len(set(raw.ID))>1:
            more_than_one_id = more_than_one_id + [raw.copy(deep=True)]

        # Convert and rename datetime column
        raw["date"] = pd.to_datetime(raw["date"], format='mixed')
        # flag year,month and day
        raw["year"] = raw["date"].dt.year.astype("int")
        raw["month"] = raw["date"].dt.to_period("M")
        raw["day"] = raw["date"].dt.to_period("D")

        # First clean-up/flagging

        ## Harmonize unit to meter
        raw.loc[raw["Unit"] == "ft", "groundwater"] = raw.loc[raw["Unit"] == "ft", "groundwater"] * 0.3048

        ## Check if there is more than one parameter in time series and extract only one (table depth)
        if raw.parameter.nunique() > 1:
            outliers_parameter = outliers_parameter + [raw.copy(deep=True)]
            # we select water table level because this parameter contains more information than water depth
            raw = raw[raw.parameter == "Water level elevation a.m.s.l."].reset_index(drop=True)
            # When time series is too short, discard current time series and continue with the next time series
            if len(raw)<2:
                data_lost_b = data_lost_b + [raw.copy(deep=True)]
                continue  ## discard this time series when there is only one record left now

        ## Unrealistic values
        ### Overwrite groundwater records that are -999 or -9999 with NA
        if ((raw.groundwater == -999).any()) | ((raw.groundwater == -9999).any()):
                outliers_novalue = outliers_novalue + [raw.copy(deep=True)]
                raw = raw.drop(raw[(raw.groundwater == -999) | (raw.groundwater == -9999)].index).reset_index(drop=True)

        # When time series is too short, discard current time series and continue with the next time series
        if len(raw)<2:
            data_lost_c = data_lost_c + [raw.copy(deep=True)]
            continue

        ## Flag negative water table depth (indicating table is either above ground, artesian or measurement error)
        sign_cat_raw,_,_ = flag_negative_signs(raw.sort_values(by='date'), outliers_depth_all, outliers_depth_any)

        ## Flag autocorrelation before aggregation
        _,autoc_raw = get_max_dist(raw.sort_values(by='date'), 0.6, 30, False)

        ## Flag outliers and spikes before aggregation
        _,outliers_raw,_ = flag_outliers(raw.sort_values(by='date'), outliers_dbscan)

        ## Flag trend and slope before aggregation
        trend_raw, _ = calc_trend(raw.sort_values(by='date')  ,False, autoc_raw)

        ## Variance before aggregation
        var_raw = np.var(raw.groundwater)

        ## number of occuring years
        length_timedelta = raw["date"].iloc[0] - raw["date"].iloc[-1]
        total_days = length_timedelta.days + length_timedelta.seconds / 86400
        yn_raw = round(total_days / 365.2425)
        if yn_raw < 1:
            yn_raw = 1

        # Aggregation

        ## Categorize in yearly, monthly or daily data
        raw["aggregated_from_n_values"] = None
        diff = abs(np.diff(raw["date"]).astype("timedelta64[s]"))
        diff_x_percentile = (np.percentile(diff, PERCENTILE_THRESHOLD) / 86400).astype("int") # in Zukunft lieber nicht runden

        if np.percentile(diff, PERCENTILE_THRESHOLD) < DAILY_THRESHOLD_SECONDS:  # 2 days, da fallen trotzdem noch Lücken rein (1.5 Tage),
            # aber hat keine großen Auswirkungen, nur 9 tägliche zeitreihen werden gekürzt, wegen > 10% Lücken und nur 2 werden gedropped
            raw["interval"] = "d"
            raw["merge"] = pd.to_datetime(raw["day"].astype(str), format="%Y-%m-%d")
            threshold = PLAT_THRESHOLD_DAILY_DAYS
        elif np.percentile(diff, PERCENTILE_THRESHOLD) < MONTHLY_THRESHOLD_SECONDS:  # 41 days
            raw["interval"] = "MS"
            raw["merge"] = pd.to_datetime(raw["month"].astype(str), format="%Y-%m")
            threshold = PLAT_THRESHOLD_MONTHLY_MONTHS
        else:
            raw["interval"] = "YS"
            raw["merge"] = pd.to_datetime(raw["year"].astype(str), format="%Y")
            if yn_raw < 5:
                threshold = PLAT_THRESHOLD_YEARLY_1_YEARS  # no plateaus can be flagged
            elif (yn_raw >= 5) & (yn_raw < 10):
                threshold = PLAT_THRESHOLD_YEARLY_5_YEARS
            elif (yn_raw >= 10) & (yn_raw < 15):
                threshold = PLAT_THRESHOLD_YEARLY_10_YEARS
            elif (yn_raw >= 15) & (yn_raw < 20):
                threshold = PLAT_THRESHOLD_YEARLY_15_YEARS
            else:
                threshold = PLAT_THRESHOLD_YEARLY_20_YEARS

        ## Flag plateaus
        plat_raw,_ = flag_plateaus(raw, threshold, outliers_plateaus)

        ## Add aggregation flag: Numbers of records per aggregated day, month or year
        ts = set(raw["merge"])
        if (len(ts) == len(raw)) & (raw["merge"] == raw["date"]).all():
            raw["aggregated_from_n_values"] = None
        else:
            for t in ts:
                raw["aggregated_from_n_values"][raw["merge"] == t] = (
                    len(raw[raw["merge"] == t]))

        ## Aggregate time series
        raw2 = (
            raw[["ID", "parameter", "groundwater", "interval", "aggregated_from_n_values"]]
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
        max_dist, autoc_agg = get_max_dist(
            raw2, 0.6, 30, False
        )  # when they are time steps with exact same time, 0 is returned

        ## Check time series for maximum allowed gap

        # Calculate non-rounded time series length
        length_timedelta = raw2["date"].iloc[-1] - raw2["date"].iloc[0]
        total_days = length_timedelta.days + length_timedelta.seconds / 86400
        yn = round(total_days / 365.2425)

        if raw2.interval[0] == "d":
            max_dist = trim_max_dist(
                max_dist, MAX_DISTANCE_DAYS
            )  # maximum allowed gap between daily time steps is 3 days --> allowed distance between existing time steps is 4
            threshold = PLAT_THRESHOLD_DAILY_DAYS  # save for plateaus later
        elif raw2.interval[0] == "MS":
            max_dist = trim_max_dist(
                max_dist, MAX_DISTANCE_MONTHS
            )  # a distance of 2 months is allowed
            threshold = PLAT_THRESHOLD_MONTHLY_MONTHS
        else:
            if yn < 5:
                max_dist = trim_max_dist(
                    max_dist, MAX_DISTANCE_1_YEARS
                )  # no gaps are allowed, 366 to include leap years
                threshold = PLAT_THRESHOLD_YEARLY_1_YEARS  # no plateaus can be flagged
            elif (yn >= 5) & (yn < 10):
                max_dist = trim_max_dist(
                    max_dist, MAX_DISTANCE_5_YEARS
                )  # no gaps are allowed, 366 to include leap years
                threshold = PLAT_THRESHOLD_YEARLY_5_YEARS
            elif (yn >= 10) & (yn < 15):
                max_dist = trim_max_dist(
                    max_dist, MAX_DISTANCE_10_YEARS
                )  # a distance of 2 years (possibly including a leap year) is allowed
                threshold = PLAT_THRESHOLD_YEARLY_10_YEARS
            elif (yn >= 15) & (yn < 20):
                max_dist = trim_max_dist(
                    max_dist, MAX_DISTANCE_15_YEARS
                )  # a distance of 2 years (possibly including a leap year) is allowed
                threshold = PLAT_THRESHOLD_YEARLY_15_YEARS
            else:
                max_dist = trim_max_dist(
                    max_dist, MAX_DISTANCE_20_YEARS
                )  # a distance of 3 years (possibly including a leap year) is allowed
                threshold = PLAT_THRESHOLD_YEARLY_20_YEARS

        ## Flag negative water table depth (indicating table is either above ground, artesian or measurement error)
        sign_cat_agg,_,_ = flag_negative_signs(raw2, outliers_depth_all, outliers_depth_any)

        ## Flag outliers and spikes
        _,outliers_agg,_ = flag_outliers(raw2, outliers_dbscan)

        ## Flag plateaus
        plat_agg,_ = flag_plateaus(raw2, threshold, outliers_plateaus)

        ## Flag trend and slope
        trend_agg, _ = calc_trend(raw2, True, autoc_agg)

        ## Variance after aggregation
        var_agg = np.var(raw2.groundwater)

        ## Extract longest sequence with intervals below max_dist
        diff = (
            abs(np.diff(raw2["date"]).astype("timedelta64[s]")).astype("int") / 86400
        )  # interval between every time step in days
        raw3 = extract_seq(raw2, diff, "int", "<=", max_dist, append=True)
        if len(raw3)!= len(raw2):
            trim_gaplength = 1 - ((len(raw2) - len(raw3)) / len(raw2))
        else:
            trim_gaplength = 1

        ## When time series is too short, discard current time series and continue with the next time series
        if len(raw3) < 2:
            data_lost_e = data_lost_e + [raw2.copy(deep=True)]
            continue

        ## Fill gaps
        raw4 = fill_gaps(raw3)

        ## Maximum gap length
        isna = raw4.groundwater.isna()  # get index where groundwater record is missing (NA)
        blocks = (
            ~isna
        ).cumsum()  # creates groups for consecutive values (here only True or False - groups)
        gaps = isna.groupby(blocks).sum()  # calculate length of continuous NA-sequences
        max_gap = gaps.max()
        mean_gap = gaps[gaps != 0].median()

        ## Drop timeseries with more than 10 % gaps
        gap_amount = len(raw4[raw4.groundwater.isna()]) / len(raw4)
        if gap_amount > GAP_FRACTION_THRESHOLD:
            raw4["year"] = raw4["date"].dt.year.astype("int")
            raw4_group = raw4.groupby("year")
            raw4_list = [raw4_group.get_group(x) for x in raw4_group.groups]
            gaps_year = [len(x[x.groundwater.isna()]) / len(x) for x in raw4_list]
            sequence = np.where((abs(np.array(gaps_year)) <= GAP_FRACTION_THRESHOLD).tolist())[0]
            longest_seq = max(
                np.split(sequence, np.where(np.diff(sequence) != 1)[0] + 1), key=len)
            if len(longest_seq) == 0:
                data_lost_f = data_lost_f + [raw3.copy(deep=True)]
                continue
            years_list = list(np.unique(raw4.year))
            years_list = years_list[longest_seq.min() : longest_seq.max() + 1]
            raw5 = raw4[raw4.year.isin(years_list)].reset_index(drop=True)
            # When time series is too short, discard current time series and continue with the next time series
            if len(raw5) < 2:
                data_lost_f = data_lost_f + [raw3.copy(deep=True)]
                continue # data is discarded
            raw5.drop("year", axis=1, inplace=True)
            gap_amount = len(raw5[raw5.groundwater.isna()]) / len(raw5)
        else:
            raw5 = raw4

        if len(raw5)!= len(raw4):
            trim_amount = 1 - ((len(raw4) - len(raw5)) / len(raw4))
        else:
            trim_amount = 1

        # Flags

        ## Flag negative water table depth (indicating table is either above ground, artesian or measurement error)
        sign_cat_final, outliers_depth_all, outliers_depth_any = flag_negative_signs(raw5, outliers_depth_all, outliers_depth_any)

        ## Flag autocorrelation before aggregation
        _, autoc_final = get_max_dist(raw5.dropna(subset="groundwater").reset_index(), 0.6, 30, False)

        ## Flag outliers detected with DBSCAN
        raw5,outliers_final, outliers_dbscan = flag_outliers(raw5, outliers_dbscan)

        ## Flag plateaus
        plat_final, outliers_plateaus = flag_plateaus(raw5,threshold,outliers_plateaus)
        if raw5.plateaus.any() != 0:
            plat = True
        else:
            plat = False

        ## add trend; trend direction and slope sign are switched when the reference point is water table depth
        ## If autocorrelation (autoc_agg=True) was detected, Ramed and Hao method is used
        trend_final, slope_final = calc_trend(raw5,True,autoc_final)

        # Add year, month and country again
        raw5["year"] = raw5["date"].dt.year.astype("int")
        raw5["month"] = raw5["date"].dt.to_period("M")
        raw5["country"] = file[
            config["country_name_pos"] : config["country_name_pos"] + 3
        ]

        ## number of occuring years
        length_timedelta = raw5["date"].iloc[-1] - raw5["date"].iloc[0]
        total_days = length_timedelta.days + length_timedelta.seconds / 86400
        yn_final = round(total_days / 365.2425)
        if yn_final < 1:
            yn_final = 1

        # Reorder columns
        raw5 = raw5[
            [
                "ID",
                "country",
                "interval",
                "date",
                "year",
                "month",
                "aggregated_from_n_values",
                "outliers",
                "plateaus",
                "parameter",
                "groundwater",
                "groundwater_filled",
            ]
        ]

        # Create attributes table
        attributes = pd.DataFrame(
            {
                "ID": [raw5.loc[0, "ID"]],
                "country": [raw5.loc[0, "country"]],
                "interval": [raw5.loc[0, "interval"]],
                "starting_date": [raw5.loc[0, "date"]],
                "ending_date": [raw5.loc[len(raw5) - 1, "date"]],
                "length_years_raw": [yn_raw],
                "length_years": [yn_final],
                "reference_point": [raw5.loc[0, "parameter"]],
                "autocorrelation_raw": autoc_raw,
                "autocorrelation_agg": autoc_agg,
                "autocorrelation": autoc_final,
                "diff_x_percentile ": diff_x_percentile,
                "aggregated_from_n_values_median": [raw5.aggregated_from_n_values.median()],
                "var_raw": var_raw,
                "var_agg": var_agg,
                "gap_fraction": [gap_amount],
                "trim_gaplength": trim_gaplength,
                "trim_gapamount": trim_amount,
                "sign_cat_raw": sign_cat_raw,
                "sign_cat_agg": sign_cat_agg,
                "sign_cat": sign_cat_final,
                "outliers_raw": [outliers_raw],
                "outliers_agg": [outliers_agg],
                "outliers": [outliers_final],
                "plateaus_raw": plat_raw,
                "plateaus_agg": plat_agg,
                "plateaus": plat_final,
                "trend_raw": [trend_raw],
                "trend_agg": [trend_agg],
                "trend_direction": [trend_final],
                "trend_slope_m_year-1": [slope_final],
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
            "outliers_novalue",
            "outliers_dbscan",
            "outliers_plateaus",
            "more_than_one_id",
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
            len(outliers_novalue) / counter,
            len(outliers_dbscan) / counter,
            len(outliers_plateaus) / counter,
            len(more_than_one_id) / counter,
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
            len(outliers_novalue),
            len(outliers_dbscan),
            len(outliers_plateaus),
            len(more_than_one_id),
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
         data_lost_f,outliers_plateaus,outliers_dbscan,outliers_depth_all,
         outliers_depth_any,outliers_parameter,more_than_one_id]
paths = ["lost_<2_records","lost_mul_par","lost_unrealistic_value","lost_after_aggregation","lost_gap_length",
         "lost_gap_amount","plateau","outliers","depth_all",
         "depth_any","par","more_id"]

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
