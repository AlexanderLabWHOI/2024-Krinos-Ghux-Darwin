# Intraspecific diversity in thermal performance determines phytoplankton ecological niche
### Code and data repository

Code for Krinos et al. _G. huxleyi_ vital rates paper with Darwin model analysis.

## Directory: `code`
Jupyter notebooks used to generate final figures for the paper.

1. `01-map-generation.ipynb` - contains code for creating the map of the isolation locations of the strains of _G huxleyi_ included in the study
2. `02-ehux-save.ipynb` - contains code for creating side scatter plots for supplemental figures and for performing thermal performance curve fitting
3. `02-ehux-vital-rates.ipynb` - contains code for thermal performance curve fitting and for all parts of Figure 2
4. `02-anderson_compare.ipynb` - code for comparing generated data to Anderson et al. (2021) and van Dassow et al. (2021) and generating supplementary figures
5. `03-new-darwin_maps-main-fig.ipynb` - code for creating main text figure 3 & supplementary figures that show generalist-specialist vs. specialist-only scenarios
6. `04-strain-habitat-prediction.ipynb` - code for predicting strain habitat based on Darwin model output
7. `XX-new-darwin_maps-mostconserv-static-generalist.ipynb` - code for assorted supplementary figures related to Darwin model simulations, especially generalist-only vs. specialist-only
8. `XX-darwin-no-v-large-penalty.ipynb` - code to compare the 3 penalty levels for the generalist-specialist simulation as shown in supplementary figure 4

### `darwin-processing-scripts`
Code for processing Darwin model simulation results outside of Jupyter notebooks. Referenced in Darwin simulation Jupyter notebooks above.

## Directory: `data`
Data tables useful for running code, including supplementary data tables included with submitted paper.

1. `00-ehux-all-strains.csv` - metadata associated with _G huxleyi_ strains
2. `01-concentration-table.csv` - daily cell (cells/mL) concentration of _G huxleyi_ cells