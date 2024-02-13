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

OUTPUTDIR="compare_plots_generalist_v_specialist"
os.system("mkdir -p "+OUTPUTDIR)
staticwidthsdir_generalist="/proj/omics/alexander/akrinos/2023-Krinos-Ghux-Darwin/data/Darwin_Simulations/StaticWidths_FixedHet_MostConservative_Generalist"
staticwidthsdir_specialist="/proj/omics/alexander/akrinos/2023-Krinos-Ghux-Darwin/data/Darwin_Simulations/StaticWidths_FixedHet_MostConservative"


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

# scenario D
list_of_temp_opt=[0.0,3.5,7,10.5,14,17.5,21,24.5,28,31.5]
list_of_a=sorted([0.30]*10)
list_of_a.reverse()
list_of_b=[0.066171077]*10
list_of_w=sorted([16]*10)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"StaticWidths_FixedHet_MostConservative_Specialist",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

# scenario D
list_of_temp_opt=[0.0,3.5,7,10.5,14,17.5,21,24.5,28,31.5]
list_of_a=sorted([0.25]*10)
list_of_a.reverse()
list_of_b=[0.066171077]*10
list_of_w=sorted([26]*10)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"StaticWidths_FixedHet_MostConservative_Generalist",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

p=(plotnine.ggplot(norberg_scenarios.loc[norberg_scenarios.Scenario.isin(["StaticWidths_FixedHet_MostConservative_Specialist",
                                                                       "StaticWidths_FixedHet_MostConservative_Generalist"])])+
    plotnine.geom_point(plotnine.aes(x="w",y="opt",color="Tracer"))+
    plotnine.scale_color_gradient(low="blue",high="red")+
    plotnine.theme_bw(base_size=12)+plotnine.ylab("Thermal optimum (oC)")+plotnine.xlab("Thermal width (oC)")+
    plotnine.facet_wrap("Scenario",nrow=2))

plotnine.ggsave(filename=os.path.join(OUTPUTDIR,"norberg_widths.png"),plot=p)

all_all_depths_StaticRange_Generalist=pd.read_csv(os.path.join(staticwidthsdir_generalist,
                                                    "ScenarioB_depths_2Oct2023_temp_tave_monthly_all_biomass.csv"))

all_all_depths_StaticRange_Specialist=pd.read_csv(os.path.join(staticwidthsdir_specialist,
                                                         "ScenarioB_depths_2Oct2023_temp_tave_monthly_all_biomass.csv"))

small_static=all_all_depths_StaticRange_Generalist.\
    loc[:,["Time","Longitude","Latitude","biomass"]].\
    groupby(["Time","Latitude","Longitude"]).biomass.sum().reset_index().\
    rename({"biomass":"biomass_Generalist"},axis="columns")
print(small_static.head(),flush=True)
merged_two=all_all_depths_StaticRange_Specialist.\
    loc[:,["Time","Longitude","Latitude","biomass"]].\
    groupby(["Time","Latitude","Longitude"]).biomass.sum().reset_index().\
    rename({"biomass":"biomass_Specialist"},axis="columns").\
    merge(small_static,left_on=["Time","Latitude","Longitude"],
          right_on=["Time","Latitude","Longitude"])

p=(plotnine.ggplot(merged_two.loc[(merged_two.biomass_Generalist>0)|(merged_two.biomass_Specialist>0)])+\
    plotnine.geom_point(plotnine.aes(x="biomass_Generalist",y="biomass_Specialist",color="Latitude"))+
    plotnine.geom_abline(plotnine.aes(slope=1,intercept=0)))

plotnine.ggsave(filename=os.path.join(OUTPUTDIR,"merged_comparison.png"),plot=p)

merged_two.to_csv(os.path.join(OUTPUTDIR,"merged_scenarios.csv"))

small_static_all=all_all_depths_StaticRange_Generalist.\
    loc[:,["Time","Longitude","Latitude","biomass"]].\
    groupby(["Latitude","Longitude"]).biomass.sum().reset_index().\
    rename({"biomass":"biomass_Generalist"},axis="columns")
print(small_static.head(),flush=True)
merged_two_all=all_all_depths_StaticRange_Specialist.\
    loc[:,["Time","Longitude","Latitude","biomass"]].\
    groupby(["Latitude","Longitude"]).biomass.sum().reset_index().\
    rename({"biomass":"biomass_Specialist"},axis="columns").\
    merge(small_static_all,left_on=["Latitude","Longitude"],
          right_on=["Latitude","Longitude"])
merged_two_all.to_csv(os.path.join(OUTPUTDIR,"merged_scenarios_all.csv"))

all_all_depths_StaticRange_Generalist["Scenario"] = "StaticWidths_FixedHet_MostConservative_Generalist"

all_all_depths_StaticRange_Generalist = all_all_depths_StaticRange_Generalist.\
    merge(norberg_scenarios.loc[norberg_scenarios.Scenario.isin(["StaticWidths_FixedHet_MostConservative_Specialist",
                                           "StaticWidths_FixedHet_MostConservative_Generalist"]),["w","opt","a",
                                                                                                  "Scenario","Tracer"]]\
          .drop_duplicates(),
          left_on=["Scenario","Tracer"],right_on=["Scenario","Tracer"])


all_all_depths_StaticRange_Specialist["Scenario"] = "StaticWidths_FixedHet_MostConservative_Specialist"

all_all_depths_StaticRange_Specialist = all_all_depths_StaticRange_Specialist.\
    merge(norberg_scenarios.loc[norberg_scenarios.Scenario.isin(["StaticWidths_FixedHet_MostConservative_Specialist",
                                           "StaticWidths_FixedHet_MostConservative_Generalist"]),["w","opt","a",
                                                                                                  "Scenario","Tracer"]]\
          .drop_duplicates(),
          left_on=["Scenario","Tracer"],right_on=["Scenario","Tracer"])

small_static_all=all_all_depths_StaticRange_Generalist.\
    loc[:,["Time","Longitude","Latitude","Tracer",
           "biomass"]].\
    groupby(["Latitude","Longitude","Tracer"]).biomass.sum().reset_index().\
    rename({"biomass":"biomass_Generalist"},axis="columns")
print(small_static.head(),flush=True)
merged_two_all=all_all_depths_StaticRange_Specialist.\
    loc[:,["Time","Longitude","Latitude","Tracer","biomass"]].\
    groupby(["Latitude","Longitude","Tracer"]).biomass.sum().reset_index().\
    rename({"biomass":"biomass_Specialist"},axis="columns").\
    merge(small_static_all,left_on=["Latitude","Longitude","Tracer"],
          right_on=["Latitude","Longitude","Tracer"])
merged_two_all.to_csv(os.path.join(OUTPUTDIR,"merged_scenarios_w_tracer.csv"))

p=(plotnine.ggplot(merged_two_all.loc[(merged_two_all.biomass_Generalist>0)|(merged_two_all.biomass_Specialist>0)])+\
    plotnine.geom_point(plotnine.aes(x="biomass_Generalist",y="biomass_Specialist",color="Latitude"))+
    plotnine.geom_abline(plotnine.aes(slope=1,intercept=0)))

plotnine.ggsave(filename=os.path.join(OUTPUTDIR,"merged_comparison_all.png"),plot=p)

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages(os.path.join(OUTPUTDIR,"merged_map.pdf"))

fig = plt.figure(figsize=(8, 8))
m = Basemap(projection='robin',lon_0=0,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='gray',lake_color='white')
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary(fill_color='white')

for_plot=merged_two_all.loc[(merged_two_all.biomass_Generalist>0)|\
                            (merged_two_all.biomass_Specialist>0)]
for_plot["biomass"] = (for_plot["biomass_Generalist"] - for_plot["biomass_Specialist"])/for_plot["biomass_Specialist"]

# Map (long, lat) to (x, y) for plotting
x, y = m(-122.3, 47.6)
norm = matplotlib.colors.Normalize(vmin=np.min(for_plot.biomass),
                                   vmax=np.max(for_plot.biomass))

cmap=matplotlib.cm.viridis
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

cbar = m.colorbar(sm)
cbar.set_label("Percent difference biomass, generalist vs. specialist")

x, y = m(for_plot.Longitude, for_plot.Latitude)
colormesh = plt.scatter(x,y,
         color=cmap(norm(for_plot.biomass.values)),s=0.5)

pp.savefig(fig)
pp.close()
plt.show()
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

cbar = m.colorbar(sm,
                  location='bottom',pad="5%")
cbar.set_label('mm')
