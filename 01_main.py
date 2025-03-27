''' @Ricarda: Im main sollen alle Skripte, die benötigt wurden, um GROW zu prozessieren, hintereinander zum laufen gebracht werden.
Kannst du dafür bitte hier eine config erstellen, die dann als pickle exportieren (siehe Beispiel unten), auf die dann alle anderen Skripte zugreifen,
indem sie das dictionary wieder einladen.

import pickle

with open('saved_dictionary.pkl', 'wb') as f:
    pickle.dump(dictionary, f)

with open('saved_dictionary.pkl', 'rb') as f:
l
    loaded_dict = pickle.load(f)

Anschließend sollen die einzelnen Skripte nacheinander durchlaufen. Dafür kannst du

# importing subprocess module
import subprocess

# running other file using run()
subprocess.run(["python", "file_1.py"])

nutzen.
'''

# Processing of IGRACs Groundwater Time series data (02_processing_gw_time_series)
pass

# Processing of IGRACs Groundwater attributes (04_processing_gw_attributes)
pass

# Merge earth system time series and attributes to groundwater data
pass

'''In this script there is the possibility to use several computer cores (can be adjusted in config). Can you please test if it still 
works when only one core is chosen (use on a normal computer and not a server)'''