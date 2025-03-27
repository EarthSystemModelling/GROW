## Takes all well attributes excel files ("wells.xlsx" and "drilling_and_construction.xlsx") from every country and merges the individual sheets it in one containing every well

''' @ me: Add description of whole script'''

# Configuration
config = {
    "basepath" : "/mnt/storage/grow/Groundwater/", # path to directory in which the folder with IGRACs groundwater data is located
    "wells": "Well_And_Monitoring_Data", # folder with IGRACs groundwater data
}

# Import packages
import os
import pandas as pd


# Derive well attributes

folders = os.listdir(config["basepath"] + config["wells"])

hydrogeo = []
management = []
construction = []
water_strike = []
log = []
structure = []

for fold in folders:
    file = config["basepath"] + config["wells"] + "/" + fold + "/" + "wells.xlsx"

    well = pd.read_excel(file, engine = 'openpyxl', sheet_name=1, skiprows= [1,1])
    if (well.empty == False):
        hydrogeo.append(well)

    well = pd.read_excel(file, engine='openpyxl', sheet_name=2, skiprows=[1, 1])
    if (well.empty==False):
        management.append(well)

    file = config["basepath"] + config["wells"] + "/" + fold + "/" + "drilling_and_construction.xlsx"

    well = pd.read_excel(file, engine='openpyxl', sheet_name=0, skiprows=[1, 1])
    if (well.empty == False):
        construction.append(well)

    well = pd.read_excel(file, engine='openpyxl', sheet_name=1, skiprows=[1, 1])
    if (well.empty == False):
        water_strike.append(well)

    well = pd.read_excel(file, engine='openpyxl', sheet_name=2, skiprows=[1, 1])
    if (well.empty == False):
        log.append(well)

    well = pd.read_excel(file, engine='openpyxl', sheet_name=3, skiprows=[1, 1])
    if (well.empty == False):
        structure.append(well)

# Dictionary with dataframe per excel sheet
all_as_df = {"hydrogeo_df":pd.concat(hydrogeo), "management_df": pd.concat(management),
             "construction_df":pd.concat(construction), "water_strike_df":pd.concat(water_strike),
             "log_df":pd.concat(log), "structure_df":pd.concat(structure)}

# Number of wells with information for every excel sheet
for sheet in all_as_df:
    print(sheet)
    print(len(all_as_df[sheet].iloc[:,1:][all_as_df[sheet].iloc[:,1:].isnull().all(1)==False]))

