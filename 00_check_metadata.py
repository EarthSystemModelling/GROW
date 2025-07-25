## 00_check_metadata

'''In the original groundwater data, the attributes per time series are stored in 2 excel files with a total of 7 excel sheets.
In this script, the number of wells that contain information (not NA) per sheet is printed. Based on this, we decided whether
to add the sheet to GROW's attribute table or not.'''

# Configuration: Path names, output names and other settings are defined here.
config = {
    "basepath" : "/mnt/storage/grow/Groundwater/", # directory in which groundwater data and all groundwater related outputs of GROW are located
    "wells": "Well_And_Monitoring_Data", # folder with IGRACs groundwater data
}

# Import packages
import os # built-in package
import pandas as pd # imported version: 2.2.3

# Derive well attributes per sheet

'''In the following, all excel files containing well attributes are sourced and merged into a table per excel sheet.
That makes one large tabe with each well per sheet. Afterwards the lengths of the rows which contain other information than the
wells ID are counted and printed.'''

folders = os.listdir(config["basepath"] + config["wells"]) # directory for each country

# Creating empty lists per excel sheet
hydrogeo = []
management = []
construction = []
water_strike = []
log = []
structure = []

# loop over every country folder
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

