import netCDF4 as nc
import os
import plotnine
import numpy as np
import pandas as pd
from datetime import datetime
os.environ["PROJ_LIB"] = "/vortexfs1/home/akrinos/.conda/envs/scplotenv/share/proj/"
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
np.float = float 

OUTPUTDIR="scenario_B"
os.system("mkdir -p "+OUTPUTDIR)
staticwidthsdir="/proj/omics/alexander/akrinos/2023-Krinos-Ghux-Darwin/data/Darwin_Simulations/StaticWidths_ConRangeTo0C"
morewidthsdir="/proj/omics/alexander/akrinos/2023-Krinos-Ghux-Darwin/data/Darwin_Simulations/ModWidths_ConRangeTo0C_MoreWidths"


### CREATE DATAFRAME TO RELATE PARAMETERS ###
def norberg(a,b,width,opt,Temps):
    return [a*np.exp(b*Temp)*(1-((Temp-opt)/(width/2))**2) for Temp in Temps]
x=list(range(-3,32))
norberg_scenarios = pd.DataFrame()
def return_frame(norberg_scenarios,scenario,\
                 list_of_temp_opt,list_of_a,list_of_b,list_of_w):
    ctr=1
    for opt,a,b,w in zip(list_of_temp_opt,list_of_a,list_of_b,list_of_w):
        curr_frame = pd.DataFrame({"x":x,"y":norberg(a,b,w,opt,x)})
        curr_frame["a"]=a
        curr_frame["b"]=b
        curr_frame["w"]=w
        curr_frame["opt"]=opt
        curr_frame["TRAC"]="TRAC"+str(ctr)
        curr_frame["Tracer"]=ctr
        curr_frame["Scenario"]=scenario
        ctr=ctr+1
        norberg_scenarios=pd.concat([norberg_scenarios,curr_frame])
        
    return norberg_scenarios

### "Simplified" set of scenarios

# scenario A
list_of_temp_opt=[17]*1
list_of_a=[0.26721485]*1
list_of_b=[0.066171077]*1
list_of_w=[24.98224]*1

norberg_scenarios=return_frame(norberg_scenarios,"A",\
                              list_of_temp_opt,list_of_a,list_of_b,list_of_w)

norberg_scenarios
# scenario B
list_of_temp_opt=[0.0000,5.00000,13.00000,14.00000,15.63879,17.00000,18.63879,19.50000,21.63879,23.00000,26.81940,32.00000]
list_of_a=[0.26721485]*12
list_of_b=[0.066171077]*12
list_of_w=[24.98224]*12
norberg_scenarios=return_frame(norberg_scenarios,"B",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

# scenario C
list_of_temp_opt=[17]*6
list_of_a=[0.5858586,0.3838384,0.2828283,0.2222222,0.1818182,0.1515152]
list_of_a.reverse()
list_of_b=[0.066171077]*6
list_of_w=[10.00000,15.00000,20.00000,24.98224,30.00000,35.00000]
norberg_scenarios=return_frame(norberg_scenarios,"C",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)


# Scenario with the full range incl. down to 0C
list_of_temp_opt=[0.00000,5.00000,13.00000,14.00000,15.63879,17.00000,18.63879,19.50000,21.63879,23.00000,26.81940,32.00000]*6
list_of_a=sorted([0.46,0.34,0.27,0.22,0.15]*12)
list_of_a.reverse()
list_of_b=[0.066171077]*60
list_of_w=sorted([15.00000,20.00000,24.98224,30.00000,40.00000]*12)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"ModWidths_FullToptRange",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

# Scenario with two different integrals and hence width sets
list_of_temp_opt=[0.00000,3.555556,
                  7.111111,10.666667,
                  14.222222,17.777778,
                  21.333333,24.888889,
                  28.444444,32.000000]*6
list_of_a=sorted([0.4,0.25,0.1,0.4,0.25,0.1]*10)
list_of_a.reverse()
list_of_b=[0.066171077]*60
list_of_w=sorted([6.711712,10.675676,25.210210,
                  12.207207,19.054054,35.000000]*10)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"ModWidthsTwoIntegral_FullToptRange",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)


# scenario D
list_of_temp_opt=[10,12,14,16,18,20,22,24,26,28,30,32]*5
list_of_a=sorted([0.46,0.34,0.27,0.22,0.15]*12)
list_of_a.reverse()
list_of_b=[0.066171077]*60
list_of_w=sorted([15.00000,20.00000,24.98224,30.00000,40.00000]*12)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"ModWidths_NoLowTopts",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

# scenario D
list_of_temp_opt=[10,11.4545454545455,12.9090909090909,14.3636363636364,15.8181818181818,
                  17.2727272727273,18.7272727272727,20.1818181818182,21.6363636363636,
                  23.0909090909091,24.5454545454545,26]*5
list_of_a=sorted([0.46,0.34,0.27,0.22,0.15]*12)
list_of_a.reverse()
list_of_b=[0.066171077]*60
list_of_w=sorted([15.00000,20.00000,24.98224,30.00000,40.00000]*12)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"ModWidths_ConstrainedRange",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

# scenario D
list_of_temp_opt=[8,10,12,14,16,18,20,22,24,26]*6
list_of_a=sorted([0.585585585585586,0.435435435435435,0.34034034034034,
                  0.28028028028028,0.235235235235235,0.2002002002002]*10)
list_of_a.reverse()
list_of_b=[0.066171077]*60
list_of_w=sorted([12,16,20,24,28,32]*10)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"ModWidths_ConRange_MoreWidths",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

# scenario D
list_of_temp_opt=[4,6.5,9,11.5,14,16.5,19,21.5,24,26.5]*6
list_of_a=sorted([0.435435435435435,0.385,0.34034034034034,
                  0.28028028028028,0.235235235235235,0.2002002002002]*10)
list_of_a.reverse()
list_of_b=[0.066171077]*60
list_of_w=sorted([16,18,20,24,28,32]*10)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"ModWidths_ConRangeTo4C_MoreWidths",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

# scenario D
list_of_temp_opt=[0.0,3.5,7,10.5,14,17.5,21,24.5,28,31.5]*6
list_of_a=sorted([0.435435435435435,0.385,0.34034034034034,
                  0.28028028028028,0.235235235235235,0.2002002002002]*10)
list_of_a.reverse()
list_of_b=[0.066171077]*60
list_of_w=sorted([16,18,20,24,28,32]*10)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"ModWidths_ConRangeTo0C_MoreWidths",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

# scenario D
list_of_temp_opt=[0.0,3.5,7,10.5,14,17.5,21,24.5,28,31.5]
list_of_a=sorted([0.28028028028028]*10)
list_of_a.reverse()
list_of_b=[0.066171077]*10
list_of_w=sorted([24]*10)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"StaticWidths_ConRangeTo0C",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

all_all_depths_ModWidths_ConRange_MoreWidths=pd.read_csv(os.path.join(morewidthsdir,
                                                         "ScenarioB_depths_2Oct2023_temp_tave_monthly_all_biomass.csv"))

def listset(input_list):
    return list(set(input_list))

# filter types that constitute at least 1% of total biomass
total_biomass=all_all_depths_ModWidths_ConRange_MoreWidths.groupby(["Latitude","Longitude","Tracer","Time"]).\
    biomass.sum().reset_index().rename({"biomass":"total_biomass"},axis="columns")
merged_weighted=all_all_depths_ModWidths_ConRange_MoreWidths.merge(total_biomass)
number_coexisting_topts=merged_weighted.\
    loc[merged_weighted.biomass/\
        merged_weighted.total_biomass > 0.01].groupby(["Latitude",
                                                       "Longitude"]).\
    Tracer.agg(listset).reset_index().explode("Tracer")

os.system("mkdir -p num_coexisting_widths")
number_coexisting_topts.to_csv(os.path.join("num_coexisting_widths","coexisting_widths.csv"))