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

OUTPUTDIR="compare_plots_fixed_het_generalist"
os.system("mkdir -p "+OUTPUTDIR)
staticwidthsdir="/proj/omics/alexander/akrinos/2023-Krinos-Ghux-Darwin/data/Darwin_Simulations/StaticWidths_FixedHet_MostConservative_Generalist"
morewidthsdir="/proj/omics/alexander/akrinos/2023-Krinos-Ghux-Darwin/data/Darwin_Simulations/ModWidths_FixedHet_MostConservative"


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
list_of_temp_opt=[0.0,3.5,7,10.5,14,17.5,21,24.5,28,31.5]*6
list_of_a=sorted([0.30,0.29,0.28,0.27,0.26,0.25]*10)
list_of_a.reverse()
list_of_b=[0.066171077]*60
list_of_w=sorted([16,18,20,22,24,26]*10)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"ModWidths_FixedHet_MostConservative",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

# scenario D
list_of_temp_opt=[0.0,3.5,7,10.5,14,17.5,21,24.5,28,31.5]
list_of_a=sorted([0.30]*10)
list_of_a.reverse()
list_of_b=[0.066171077]*10
list_of_w=sorted([16]*10)
list_of_w=list_of_w
norberg_scenarios=return_frame(norberg_scenarios,"StaticWidths_FixedHet_MostConservative",\
                               list_of_temp_opt,list_of_a,list_of_b,list_of_w)

p=(plotnine.ggplot(norberg_scenarios.loc[norberg_scenarios.Scenario.isin(["StaticWidths_FixedHet_MostConservative",
                                                                       "ModWidths_FixedHet_MostConservative"])])+
    plotnine.geom_point(plotnine.aes(x="w",y="opt",color="Tracer"))+
    plotnine.scale_color_gradient(low="blue",high="red")+
    plotnine.theme_bw(base_size=12)+plotnine.ylab("Thermal optimum (oC)")+plotnine.xlab("Thermal width (oC)")+
    plotnine.facet_wrap("Scenario",nrow=2))

plotnine.ggsave(filename=os.path.join(OUTPUTDIR,"norberg_widths.png"),plot=p)

all_all_depths_StaticRange=pd.read_csv(os.path.join(staticwidthsdir,
                                                    "ScenarioB_depths_2Oct2023_temp_tave_monthly_all_biomass.csv"))

all_all_depths_ModWidths_ConRange_MoreWidths=pd.read_csv(os.path.join(morewidthsdir,
                                                         "ScenarioB_depths_2Oct2023_temp_tave_monthly_all_biomass.csv"))

small_static=all_all_depths_StaticRange.\
    loc[:,["Time","Longitude","Latitude","biomass"]].\
    groupby(["Time","Latitude","Longitude"]).biomass.sum().reset_index().\
    rename({"biomass":"biomass_A"},axis="columns")
print(small_static.head(),flush=True)
merged_two=all_all_depths_ModWidths_ConRange_MoreWidths.\
    loc[:,["Time","Longitude","Latitude","biomass"]].\
    groupby(["Time","Latitude","Longitude"]).biomass.sum().reset_index().\
    rename({"biomass":"biomass_B"},axis="columns").\
    merge(small_static,left_on=["Time","Latitude","Longitude"],
          right_on=["Time","Latitude","Longitude"])

p=(plotnine.ggplot(merged_two.loc[(merged_two.biomass_A>0)|(merged_two.biomass_B>0)])+\
    plotnine.geom_point(plotnine.aes(x="biomass_A",y="biomass_B",color="Latitude"))+
    plotnine.geom_abline(plotnine.aes(slope=1,intercept=0)))

plotnine.ggsave(filename=os.path.join(OUTPUTDIR,"merged_comparison.png"),plot=p)

merged_two.to_csv(os.path.join(OUTPUTDIR,"merged_scenarios.csv"))

small_static_all=all_all_depths_StaticRange.\
    loc[:,["Time","Longitude","Latitude","biomass"]].\
    groupby(["Latitude","Longitude"]).biomass.sum().reset_index().\
    rename({"biomass":"biomass_A"},axis="columns")
print(small_static.head(),flush=True)
merged_two_all=all_all_depths_ModWidths_ConRange_MoreWidths.\
    loc[:,["Time","Longitude","Latitude","biomass"]].\
    groupby(["Latitude","Longitude"]).biomass.sum().reset_index().\
    rename({"biomass":"biomass_B"},axis="columns").\
    merge(small_static_all,left_on=["Latitude","Longitude"],
          right_on=["Latitude","Longitude"])
merged_two.to_csv(os.path.join(OUTPUTDIR,"merged_scenarios_all.csv"))

p=(plotnine.ggplot(merged_two_all.loc[(merged_two_all.biomass_A>0)|(merged_two_all.biomass_B>0)])+\
    plotnine.geom_point(plotnine.aes(x="biomass_A",y="biomass_B",color="Latitude"))+
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

for_plot=merged_two_all.loc[(merged_two_all.biomass_A>0)|\
                            (merged_two_all.biomass_B>0)]
for_plot["biomass"] = (for_plot["biomass_A"] - for_plot["biomass_B"])/for_plot["biomass_B"]

# Map (long, lat) to (x, y) for plotting
x, y = m(-122.3, 47.6)
norm = matplotlib.colors.Normalize(vmin=np.min(for_plot.biomass),
                                   vmax=np.max(for_plot.biomass))

cmap=matplotlib.cm.viridis
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

cbar = m.colorbar(sm)
cbar.set_label("Percent difference biomass, one width vs. multi-width")

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
