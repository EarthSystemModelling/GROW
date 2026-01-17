<h1 align="center">GROW 1.0 - Data Processing</h1>

<h3 align="center">The python scripts in this repository were used to prepare the GROW 1.0 dataset.</h3>

<p align="center">
<img src="2025_GROW_Blau_Web_150DPI.png" width="200" /> 
</p>
<p align="center">
<em>GROW logo: Bäthge et al. 2026 / CC BY 4.0</em>

<p align="center">GROW (global-scale integrated GROundWater package) is a global-scale, analysis-ready, quality-controlled dataset that combines grounwater depth and level time series from around the world with associated Earth system
variables. The dataset contains > 200,000 time series from 55 countries in a daily, monthly, or yearly temporal 
resolution, accompanied by 36 time series or attributes of meteorological, hydrological, geophysical.
vegetation, and anthropogenic variables (e.g., precipitation, drainage density, aquifer type, NDVI, land use).
34 data flags regarding well features (e.g., location coordinates and country), as well as time series characteristics
(e.g., gap fraction or length), facilitate quick data filtering.</p>

-------------------------------------------------------

<br/>
The data set can be found at: Insert Zendodo link

References and descriptions to the original data can be found in:

(Bäthge et al.2026) ~ here I will cite the paper

**The workflow to create GROW is implemented in the following scripts which should be run sequentially in the order:**

1. 01_processing_gw_time_series.py
2. 02_processing_gw_attributes.py
3. 03_merge_earth_system_variables

A description of each script's function is given below:

### 00_check_metadata.py

In the original groundwater data, the attributes per time series are stored in 2 ods files with a total of 7 sheets.
In this script, the number of wells that contain information (not NA) per sheet is printed. Based on this, we decided whether to add the sheet to GROW's attribute table or not.

### 01_processing_gw_time_series.py

A single table containing all time series is created. 
Every time series of the original data is stored in an individual ods file. In this script, each time series file is consecutively read, processed, quality-checked, eventually discarded or appended to a single table. 
Data flags are generated and either added as column to the time series data or exported to a time series attributes table.
Other than that, discarded time series are collected in extra tables and exported. Time series statistics are derived and exported.

### func_processing_gw_time_series.py

This file contains all functions that are used in "01_processing_gw_time_series.py".

### 02_processing_gw_attributes.py

The attributes table and the final time series table are created.
Firstly, all files containing static (no temporal dimension) well attributes per country are read and merged together.
A) duplicates by ID and country, b) duplicates by coordinates, starting date, ending date and mean groundwater
table and c) duplicates originating from Jasechko et al. 2024 are removed. 
Further, wells with coordinates outside realistic ranges are removed. Attributes generated in the
time series processing are merged to the attributes from the original data. The time series table is trimmed to 
the ID's that are still left in the preprocessed attributes table. 
A unique GROW ID is created.

### 03_merge_earth_system_variables.py

37 Earth system variables (attributes and time series) are added to the attributes or time series table.
First, 20 attributes are consecutively added to the groundwater data. Afterwards, 19 time series variables
are aggregated and merged to the groundwater time series within a parallelized process (server with multiple cores recommended). 
The time series table is split into 100 parts for which the variables are added in parallel. In the end,
all parts are put together again.

### func_merge_earth_system_variables.py

This file contains all functions that are used in "03_merge_earth_system_variables.py".

### usage_example.py

A quick example demonstrates how GROW can be subset and prepared before starting an analysis.
