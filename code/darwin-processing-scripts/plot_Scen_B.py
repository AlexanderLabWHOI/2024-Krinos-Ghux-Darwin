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

p=(plotnine.ggplot(norberg_scenarios.loc[norberg_scenarios.Scenario.isin(["StaticWidths_ConRangeTo0C",
                                                                       "ModWidths_ConRangeTo0C_MoreWidths"])])+
    plotnine.geom_point(plotnine.aes(x="w",y="opt",color="Tracer"))+
    plotnine.scale_color_gradient(low="blue",high="red")+
    plotnine.theme_bw(base_size=12)+plotnine.ylab("Thermal optimum (oC)")+plotnine.xlab("Thermal width (oC)")+
    plotnine.facet_wrap("Scenario",nrow=2))

plotnine.ggsave(filename=os.path.join(OUTPUTDIR,"norberg_widths.png"),plot=p)

all_all_depths_ModWidths_ConRange_MoreWidths=pd.read_csv(os.path.join(morewidthsdir,
                                                         "ScenarioB_depths_2Oct2023_temp_tave_monthly_all_biomass.csv"))


temps_ModWidths_ConRange_MoreWidths=pd.read_csv(os.path.join(morewidthsdir,
                                    "tave.ScenarioB_no_low_opt.revise.concatenated_temp_phys_tave_monthly_all_biomass.csv"))

temps_ModWidths_ConRange_MoreWidths["DateTime"]=pd.to_datetime(temps_ModWidths_ConRange_MoreWidths.Time)

print("STEP 1 COMPLETE", flush=True)
def calcthermrange(input_series):
    return np.max(input_series)-np.min(input_series)

def calcthermvar(input_series):
    return np.std(input_series)

temps_ModWidths_ConRange_MoreWidths["basin"] = ["Southern" if (lat < -50) else "Arctic" if \
                                             (lat>50) else "Atlantic" if ((long<15)|(long>300))\
                                             else "Pacific" if (long>150)&(long<300) else "Indian"\
                                             for long,lat in zip(temps_ModWidths_ConRange_MoreWidths["Longitude"],
                                                                 temps_ModWidths_ConRange_MoreWidths["Latitude"])]


combine_w_temp=temps_ModWidths_ConRange_MoreWidths.groupby(["Latitude","Longitude","basin"])["Tpot_tave"]\
                 .agg(calcthermvar).reset_index().merge(all_all_depths_ModWidths_ConRange_MoreWidths,
                                                        left_on=["Latitude","Longitude"],
                                       right_on=["Latitude","Longitude"])

combine_all_w_temp=temps_ModWidths_ConRange_MoreWidths.reset_index().loc[:,["Latitude","Longitude",
                                                                            "Time","Tpot_tave"]].\
                    merge(all_all_depths_ModWidths_ConRange_MoreWidths,left_on=["Latitude","Longitude","Time"],
                          right_on=["Latitude","Longitude","Time"])

merged_all_w_temp_scen=combine_all_w_temp.merge(norberg_scenarios[["a","b","w","opt","Tracer","Scenario"]].\
              loc[norberg_scenarios["Scenario"]=="ModWidths_ConRangeTo0C_MoreWidths"].\
              drop_duplicates(),how="inner",left_on=["Tracer"],
                                       right_on=["Tracer"])
merged_all_w_temp_scen["ToptDist"] =  merged_all_w_temp_scen["opt"]-merged_all_w_temp_scen["Tpot_tave"]
merged_all_w_temp_scen.to_csv(os.path.join(OUTPUTDIR,"merged_weighted_widths_step2.csv"))

curr_df = merged_all_w_temp_scen.copy() 
print("STEP 2 COMPLETE", flush=True)

total_biomass=curr_df.groupby(["Latitude","Longitude","Time"]).biomass.sum().\
        reset_index().rename({"biomass":"total_biomass"},axis="columns")

merged_weighted=curr_df.merge(total_biomass)


print("STEP 2.5 COMPLETE", flush=True)

merged_weighted["calculate_offset"]=merged_weighted["w"]*(merged_weighted["biomass"]/\
                                                            merged_weighted["total_biomass"])
merged_weighted["opt_offset"]=merged_weighted["opt"]*(merged_weighted["biomass"]/\
                                                            merged_weighted["total_biomass"])

merged_weighted["calculate_offset"]=[curr if (curr>0)&(curr==curr) else 0 for curr,curr_w in \
                                      zip(merged_weighted["calculate_offset"],merged_weighted["w"])]

merged_weighted["shannon_index"]=[(biom/tot_biom)*np.log(biom/tot_biom) if biom>0 else \
                                  0 for \
                                  biom,tot_biom in zip(merged_weighted["biomass"],
                                                       merged_weighted["total_biomass"])]

merged_weighted["distance_offset"]=merged_weighted["ToptDist"]*(merged_weighted["biomass"]/\
                                                            merged_weighted["total_biomass"])
merged_weighted.to_csv(os.path.join(OUTPUTDIR,"merged_weighted_widths_topt_dist.csv"))

comparison_temp_info=merged_weighted.loc[merged_weighted.calculate_offset!=0,].\
    groupby(["Latitude","Longitude","Time"])\
    [["calculate_offset","opt_offset","distance_offset","shannon_index"]].sum().reset_index()

print("STEP 3 COMPLETE", flush=True)


def lenset(input_list):
    return len(list(set(input_list)))
# filter types that constitute at least 1% of total biomass
number_coexisting_widths=merged_weighted.loc[merged_weighted.biomass/\
                                             merged_weighted.total_biomass > 0.01].groupby(["Latitude",
                                                                                           "Longitude"]).\
    w.agg(lenset).reset_index()#.groupby("w").count()
number_coexisting_widths.to_csv(os.path.join(OUTPUTDIR,"num_coexisting_widths.csv"))

def lenset(input_list):
    return len(list(set(input_list)))
# filter types that constitute at least 1% of total biomass
number_coexisting_topts=merged_weighted.loc[merged_weighted.biomass/\
                                             merged_weighted.total_biomass > 0.01].groupby(["Latitude",
                                                                                           "Longitude"]).\
    opt.agg(lenset).reset_index()#.groupby("w").count()

number_coexisting_topts.to_csv(os.path.join(OUTPUTDIR,"num_coexisting_topts.csv"))

def lenset(input_list):
    return len(list(set(input_list)))

# filter types that constitute at least 1% of total biomass
topt_dist_ranges=merged_weighted.copy(deep=True)
topt_dist_ranges["DistSign"] = ["BelowOpt" if ToptDist<0 else "AboveOpt" for ToptDist in \
                                topt_dist_ranges["ToptDist"]]

topt_dist_ranges["ToptDistAbs"] = np.abs(topt_dist_ranges["ToptDist"])
topt_dist_ranges=topt_dist_ranges.loc[merged_weighted.biomass/\
                                      merged_weighted.total_biomass > 0.05].groupby(["Latitude",
                                                                                           "Longitude",
                                                                                           "DistSign"]).\
    ToptDistAbs.agg([np.min,np.max,np.mean]).reset_index()#.groupby("w").count()



from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages(os.path.join(OUTPUTDIR,'num_coexisting_widths_ScenB.pdf'))

fig = plt.figure(figsize=(8, 8))
m = Basemap(projection='robin',lon_0=0,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='gray',lake_color='white')
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary(fill_color='white')

for_plot=number_coexisting_widths

# Map (long, lat) to (x, y) for plotting
x, y = m(-122.3, 47.6)
norm = matplotlib.colors.Normalize(vmin=np.min(for_plot.w),
                                   vmax=np.max(for_plot.w))

cmap=matplotlib.cm.winter
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

cbar = m.colorbar(sm)
cbar.set_label("Number of coexisting thermal widths")

x, y = m(for_plot.Longitude, for_plot.Latitude)
colormesh = plt.scatter(x,y,
         color=cmap(norm(for_plot.w.values)),s=0.5)

pp.savefig()
pp.close()
plt.show()
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

cbar = m.colorbar(sm,
                  location='bottom',pad="5%")
cbar.set_label('mm')

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages(os.path.join(OUTPUTDIR,"num_coexisting_optima_ScenB.pdf"))

fig = plt.figure(figsize=(8, 8))
m = Basemap(projection='robin',lon_0=0,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='gray',lake_color='white')
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary(fill_color='white')

for_plot=number_coexisting_topts

# Map (long, lat) to (x, y) for plotting
x, y = m(-122.3, 47.6)
norm = matplotlib.colors.Normalize(vmin=np.min(for_plot.opt),
                                   vmax=np.max(for_plot.opt))

cmap=matplotlib.cm.spring
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

cbar = m.colorbar(sm)
cbar.set_label("Number of coexisting thermal optima")

x, y = m(for_plot.Longitude, for_plot.Latitude)
colormesh = plt.scatter(x,y,
         color=cmap(norm(for_plot.opt.values)),s=0.5)

pp.savefig()
pp.close()
plt.show()
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

cbar = m.colorbar(sm,
                  location='bottom',pad="5%")
cbar.set_label('mm')

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages(os.path.join(OUTPUTDIR,'mean_Shannon_index_ScenB.pdf'))

fig = plt.figure(figsize=(8, 8))
m = Basemap(projection='robin',lon_0=0,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='gray',lake_color='white')
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary(fill_color='white')

for_plot=comparison_temp_info.groupby(["Latitude","Longitude"]).shannon_index.mean().reset_index()
for_plot["shannon_index"]=-for_plot["shannon_index"]

# Map (long, lat) to (x, y) for plotting
x, y = m(-122.3, 47.6)
norm = matplotlib.colors.Normalize(vmin=np.min(for_plot.shannon_index),
                                   vmax=np.max(for_plot.shannon_index))

cmap=matplotlib.cm.winter
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

cbar = m.colorbar(sm)
cbar.set_label("Mean Shannon index")

x, y = m(for_plot.Longitude, for_plot.Latitude)
colormesh = plt.scatter(x,y,
         color=cmap(norm(for_plot.shannon_index.values)),s=0.5)

pp.savefig()
pp.close()
plt.show()
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

cbar = m.colorbar(sm,
                  location='bottom',pad="5%")
cbar.set_label('mm')

p=(plotnine.ggplot(merged_weighted.groupby(["Latitude","Longitude","Tracer"]).\
                 ToptDist.mean().reset_index()) + 
     plotnine.geom_point(plotnine.aes(y = "Latitude", x = "Longitude", color="ToptDist"))  + 
     plotnine.theme_bw(base_size=12) + plotnine.scale_color_gradient2(low="blue",
                                                                     high="red",
                                                                     mid="white",
                                                                     midpoint=0,
                                                                     name="Distance between Topt and Tenv") + 
     plotnine.facet_wrap("Tracer"))
plotnine.ggsave(plot=p, filename=os.path.join(OUTPUTDIR,
                                              'tracers_dist_topt_noweight.png'), dpi=100, width=10,height=10,units="in")

p=(plotnine.ggplot(merged_weighted.groupby(["Latitude","Longitude","Tracer"]).\
                 biomass.mean().reset_index()) + 
     plotnine.geom_point(plotnine.aes(y = "Latitude", x = "Longitude", color="biomass"))  + 
     #plotnine.geom_line(plotnine.aes(x = "Latitude", y = "Tpot_tave", color="Longitude",group="Latitude"))+
     plotnine.theme_bw(base_size=12) + plotnine.scale_color_cmap(name="Mean biomass") + 
     plotnine.facet_wrap("Tracer"))
plotnine.ggsave(plot=p, filename=os.path.join(OUTPUTDIR,
                                              'tracers_biomass_lowest.png'), dpi=100, width=10, height=10, units="in")

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages(os.path.join(OUTPUTDIR,'width_weighted_biomass_ScenB.pdf'))

fig = plt.figure(figsize=(8, 8))
m = Basemap(projection='robin',lon_0=0,resolution='c')
m.drawcoastlines()
m.fillcontinents(color='gray',lake_color='white')
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary(fill_color='white')

for_plot=comparison_temp_info

# Map (long, lat) to (x, y) for plotting
x, y = m(-122.3, 47.6)
norm = matplotlib.colors.Normalize(vmin=np.min(for_plot.calculate_offset),
                                   vmax=np.max(for_plot.calculate_offset))

cmap=matplotlib.cm.hot
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

cbar = m.colorbar(sm)
cbar.set_label("Weighted mean thermal width")

x, y = m(for_plot.Longitude, for_plot.Latitude)
colormesh = plt.scatter(x,y,
         color=cmap(norm(for_plot.calculate_offset.values)),s=0.5)

pp.savefig()
pp.close()
plt.show()
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

cbar = m.colorbar(sm,
                  location='bottom',pad="5%")
cbar.set_label('mm')
