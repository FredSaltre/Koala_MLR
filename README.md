# Koala species distribution and fertility-reduction model
<img align="right" src="www/koala.png" alt="koala & joey" width="300" style="margin-top: 20px">
Combining a point-process model, an ensemble species distribution model, and a demographic model to project koala populations in the <a href="https://en.wikipedia.org/wiki/Mount_Lofty_Ranges">Mount Lofty Ranges</a> over the next 25 years to assess the efficiency and cost-effectiveness of fertility-control interventions.<br>
<br>
Accompanies paper:<br>
<br>
<a href="https://www.flinders.edu.au/people/frederik.saltre">Saltré, F</a>, <a href="https://scholars.uow.edu.au/katharina-peters">KJ Peters</a>, <a href="https://au.linkedin.com/in/dan-rogers-a11b4392">DJ Rogers</a>, <a href="https://www.researchgate.net/profile/Joel-Chadoeuf">J Chadœuf</a>, <a href="https://www.flinders.edu.au/people/vera.weisbecker">V Weisbecker</a>, <a href="https://www.flinders.edu.au/people/corey.bradshaw">CJA Bradshaw</a>. Balancing overpopulation and conservation targets to optimize koala management strategies. In review<br>
<br>

### Scripts & Prerequisites
- <a href="https://github.com/FredSaltre/Koala_MLR/tree/main/Demographic_model">Demographic model</a>
- <a href="https://github.com/FredSaltre/Koala_MLR/tree/main/Species_Distribution_Model">Species distribution model</a>
- <a href="https://github.com/FredSaltre/Koala_MLR/tree/main/Point_Process_Model">Point-process model</a>

This project uses **R version 4.0.4**. Ensure you have this version installed for compatibility. If you are using a later version of R, you will need to manually install the package **Imap 1.32**, which is no longer available on CRAN for direct installation. Follow these steps to install it:

  1. **Download the Package:**  
     Visit the [CRAN Archive for    Imap](https://cran.r-project.org/src/contrib/Archive/Imap/) and download the package archive file (`Imap_1.32.tar.gz`).

  2. **Install the Package Manually:**  
     Open R and use the following command to install the downloaded file:
     install.packages("~/Imap_1.32.tar.gz", repos = NULL, type = "source")

Ensure the following R libraries are installed: `geoR`, `sf`, `spatstat`, `MASS`, `proj4`, `sp`, `terra`, `dplyr`, `plotly`,`deSolve`,`biomod2`, `terra`, `ggplot2`, `rgdal`, `sf`, `raster`, `maptools`


---

#### Koala Population Projection Model

##### Overview
This script models the koala population in the Mount Lofty Ranges using a matrix population model. It incorporates stochastic factors, sterilization scenarios, and catastrophic events to project population dynamics over the years 2016–2040. The model also assesses the efficiency and cost-effectiveness of sterilization interventions.

##### Key Features
- **Population Projection**: Simulates future koala population trends under unmanaged and sterilization scenarios.
- **Stochastic Dynamics**: Introduces randomness in fertility and survival rates to account for real-world variability.
- **Catastrophe Scenarios**: Models catastrophic events affecting koala survival.
- **Sterilization Analysis**: Evaluates the impact of sterilization rates (0-90%) on population size and associated costs.
- **Cost Modelling**: Estimates costs for sterilization programs, considering population density and search efforts.

##### Outputs
- **Population Projections**: Median population size, confidence intervals, and extinction probabilities for each scenario.
- **Sterilization Costs**: Total and annual costs of fertility-control programs.
- **Population Metrics**: CSV files summarizing population size, sterilization counts, and intervention costs.

---

#### Species Distribution Model (SDM)

##### Overview
The **Species Distribution Model (SDM)** integrates environmental and species presence/absence data to estimate habitat suitability for koalas across the Mount Lofty Ranges.

##### Key Features
- **Data Preprocessing**:
  - Cleans and filters species presence data (e.g., removes records from Kangaroo Island).
  - Projects species presence data to a consistent geographic coordinate system.
- **Environmental Data Integration**:
  - Loads and processes environmental predictors (e.g., distance to roads, rainfall, vegetation cover).
  - Combines all predictors into a raster stack for analysis.
  - Checks for multicollinearity using the Variance Inflation Factor (VIF).
- **Modeling**:
  - Utilizes the `BIOMOD2` package to create individual species distribution models.
  - Implements cross-validation and evaluates variable importance.
- **Ensemble Modeling**:
  - Combines individual models into an ensemble for robust predictions.
  - Produces evaluation metrics and variable importance plots for ensemble outputs.
- **Population Density Estimation**:
  - Converts habitat suitability scores into population density estimates using regional calibration factors.

##### Outputs
- Habitat suitability raster maps.
- Evaluation metrics, including True Skill Statistic (TSS) and ROC curves.
- Ranked variable importance plots.
- Projected population density maps, including confidence intervals.

##### Script Details
- **Input Data**:
  - Species data: `KoalaData(GKC1&2).csv`
  - Environmental layers: e.g., `Dist2sealRoads.asc`, `MinTemp.asc`.

---

#### Point Process Model for Koala Spatial Analysis

##### Overview

This repository contains the implementation of a Point Process Model designed for analyzing the spatial distribution of koalas in the Mount Lofty Ranges. It integrates multiple spatial datasets and advanced statistical methods to identify factors affecting koala intensity and distribution. The model explores the influence of proximity to parks, roads, and other geographic features. This Point Process Model:
- Uses spatial point pattern analysis to assess koala distributions.
- Employs kernel density estimation (KDE) to identify koala hotspots.
- Fits Poisson models to evaluate the effects of distance from hotspots, parks, and roads.
- Generates intensity maps for visualizing koala density and habitat preferences.

##### Key Features
- **Kernel Density Estimation (KDE):** Highlights regions of high koala intensity.
- **Distance-Based Analysis:** Explores the impact of proximity to hotspots, roads, and parks.
- **Poisson Models:**
  - Model 1: Effect of distance to hotspots.
  - Model 2: Adds the influence of park boundaries.
  - Model 3: Integrates distance to roads.
- **Visualization:**
  - Intensity maps with and without park protection.
  - Graphical representation of the influence of road proximity.

##### Outputs
- **Koala Intensity Maps:**
  - Maps showing koala density across the study area.
  - Separate visualizations for the impact of park protection and road proximity.
- **Statistical Summaries:**
  - Estimated koala density per square kilometer.
  - Comparative intensity within and outside parks.
- **Distance Effects:**
  - Plots illustrating how distance from hotspots and roads influences koala intensity.

##### Script Details
- **Main Script:** `PointProcessModel.R`
- **External Dependencies:**
  - `LatLongToXY.r` script for converting latitude and longitude to Cartesian coordinates.


---

<br>
<br>
<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Stacked_RGB_Master.jpg" alt="Flinders University" width="60" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL" width="120" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://EpicAustralia.org.au"><img align="bottom-left" src="www/CabahFCL.jpg" alt="CABAH logo" height="45" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.environment.sa.gov.au/"><img align="bottom-left" src="www/SADEWlogo.png" alt="DEW" width="100" style="margin-top: 20px"></a></p>


