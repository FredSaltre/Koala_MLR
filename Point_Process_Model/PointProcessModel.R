#############################################################
#### INHOMOGENEOUS POISSON POINT PROCESS MODEL ############
#############################################################
# author: Frederik Saltre & Joel Chadoeuf
# date: 7/09/2024

#This R script performs spatial analysis on koala distribution using a combination of point pattern analysis, kernel density estimation, and Poisson modeling. It integrates various spatial datasets, including koala location data, park boundaries, and distance to roads, to explore how these factors influence koala intensity and distribution.

#Main Steps:
# 1) Loading Spatial Data:
#        The script loads koala location data and park boundary shapefiles.
#        Koala locations are filtered to exclude islands, focusing on a specific geographic area.
# 2) Kernel Density Estimation (KDE):
#        Kernel density estimation is applied to the koala locations to identify areas of high koala intensity.
#        The script identifies the highest intensity point (a proxy for the koala "hotspot") and calculates the distance from each koala to this point.
# 3) Distance-Based Analysis:
#        The script examines how the distance to the hotspot affects koala intensity.
#        It also calculates the distance from each koala to the nearest road, integrating road proximity into the analysis.
# 4) Poisson Models:
#        The script fits several Poisson models to investigate the factors influencing koala intensity:
#        Model 1: Only considers the distance to the hotspot.
#        Model 2: Considers both the distance to the hotspot and whether the koala is inside a park (park effect).
#        Model 3: Incorporates distance to the hotspot, park boundaries, and distance to the nearest road.
# 5) Koala Intensity Estimation:
#        Based on the Poisson models, the script estimates the koala density per square kilometer.
#        The results are visualized through maps showing the estimated koala intensity across the study area, with additional focus on the role of park protection and road proximity.
# 6) Visualization:
#        The script generates plots displaying the spatial distribution of koalas, including intensity maps and distance effect plots.
#        It contrasts the intensity of koalas inside and outside parks, and shows how road proximity impacts koala detection.

#Key Takeaways:
#The script explores how geographic factors like distance to roads, proximity to parks, and distance from the central "hotspot" influence koala density.
#It uses spatial modeling techniques to quantify these effects and produce maps that visualize koala distribution.
#The results provide insights into koala habitat preferences and the potential impact of human infrastructure (roads) on koala populations.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# R libraries
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load necessary libraries
library(geoR)      # For geostatistical analysis
library(sf)        # For handling spatial vector data
library(spatstat)  # For spatial point pattern analysis
library(MASS)      # For statistical functions like kde2d (kernel density estimation)
library(proj4)      # For distance calculations (geodesic distances)
library(sp) 
library(Imap)
library(terra)

# Set the grid size for the Poisson model calculations
ng <- 160  # Grid size for stability in intensity estimation

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
source("LatLongToXY.r")# Load the script to convert latitude/longitude to Cartesian coordinates

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read the shapefile for park boundaries (spatial limits of parks)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
shp <- st_read("./CONSERVATION_NpwsaReserves_SHP-2/CONSERVATION_NpwsaReserves_SHP-2/CONSERVATION_NpwsaReserves_GDA2020.shp", quiet=FALSE)
shp <- st_geometry(shp)  # Extract geometry from the shapefile

# Convert park polygons into observation windows (owin objects) for spatial analysis
lshp <- vector("list", length(shp))
for(i in 1:length(shp)) { lshp[[i]] <- as.owin(shp[[i]]) }

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load and preprocess koala location data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Koala <- read.csv("GKC1.2.long.lat.csv")
Koala <- as.ppp(Koala[,3:2],c(129,140,-38,-26))# Convert koala lat/lon to a point pattern object (ppp)
# Coordinates represent the window of the spatial analysis

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Filter koalas on islands (remove unwanted points based on coordinates)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
selx <- (Koala$x > 138.2591)  & (Koala$x < 138.9909)  # Select x-coordinates within a certain range
sely <- (Koala$y > -35.23021) & (Koala$y < -34.6995)  # Select y-coordinates within a certain range
sel <- selx & sely  # Combine x and y selections
Koala2 <- as.ppp(matrix(c(Koala$x[sel], Koala$y[sel]), ncol=2, byrow=FALSE), c(range(Koala$x[sel]), range(Koala$y[sel])))
# Create a new point pattern object with only selected koala locations

# Create a convex hull around the koala locations (spatial boundary for further analysis)
convKoala2 <- convexhull(Koala2)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Kernel Density Estimation for koala intensity
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

kde.res <- kde2d(Koala2$x, Koala2$y, n=240, lims=c(Koala2$window$x, Koala2$window$y))  # KDE with 240 grid points
kde.res$z <- sqrt(kde.res$z)  # Adjust intensity scale for visualization

# Plot koala intensity (observed density)
plot(Koala2, type="n", main="Observed Koala Intensity")  # Base plot
image(kde.res, add=TRUE)  # Add kernel density estimate to the plot

# Identify the maximum density point (center of highest koala intensity)
mh <- max(kde.res$z)
xh <- rep(kde.res$x, 240)[kde.res$z == mh]  # X-coordinate of the highest intensity point
yh <- sort(rep(kde.res$y, 240))[kde.res$z == mh]  # Y-coordinate of the highest intensity point

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Analyze the effect of distance from the highest density point (xh, yh)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dik <- gdist(Koala2$x, Koala2$y, xh, yh, unit="km")  # Calculate the geodesic distance to (xh, yh)
densdik <- density(c(-dik, dik))  # Create a density estimate of the distances
densdik$y <- log(densdik$y / abs(densdik$x))  # Log-transform the density estimate

# Plot the log-transformed distance effect
plot(densdik$x[densdik$x > 0], densdik$y[densdik$x > 0], xlab="Distance to (xh, yh)", ylab="log-density", main="Distance from Center")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load the distance to roads (optional section for handling spatial rasters)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+# Process distance to roads data
#dist2sroads.asc <- read.asciigrid("Dist2sealRoads.asc", as.image = FALSE, plot.image=T, 
#                                  proj4string = CRS(as.character("+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137.0 +b=6378137.0 +units=m +no_defs")), colname="D2SROADS")
dist2sroads <- rast("Dist2sealRoads.asc")  # Using `terra` for raster reading
crs(dist2sroads) <- "+proj=eqc +lat_ts=60 +lon_0=0 +datum=WGS84"

# Convert the raster to spatial points (SpatVector in terra)
disroute_spatvector <- as.points(dist2sroads)

# Reproject the SpatVector to geographic coordinates (WGS 84)
disroute_spatvector_transformed <- project(disroute_spatvector, "EPSG:4326")

# Extract the transformed coordinates to check
disroute_transformed_coords <- crds(disroute_spatvector_transformed)
# Combine the transformed coordinates with the original data into a data frame
disroute.dat <- data.frame(D2SROADS = values(disroute_spatvector_transformed), 
                           Lon = disroute_transformed_coords[, 1],
                           Lat = disroute_transformed_coords[, 2])

# Print the first few rows to check the transformed data
print(head(disroute.dat))

# Plot convex hull, roads, and koalas
plot(convKoala2, lwd=2, col=5, main="Complete Data Set")  # Plot the koala convex hull
points(disroute.dat$Lon, disroute.dat$Lat, pch=19, col=gray(disroute.dat$Dist2sealRoads**0.3 / max(disroute.dat$Dist2sealRoads**0.3)))  # Plot distance to roads
points(Koala2, pch=19, col=3, cex=0.6)  # Plot koala locations

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create a grid within the convex hull (spatial grid for analysis)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
xc <- range(Koala2$x)  # X-coordinate range
yc <- range(Koala2$y)  # Y-coordinate range

# Generate a grid within the range of koala locations
xx <- matrix(rep(1:ng, ng)-0.5, ncol=ng, byrow=TRUE) / ng
xx <- min(xc) + diff(xc) * as.vector(xx)
yy <- matrix(rep(1:ng, ng)-0.5, ncol=ng, byrow=FALSE) / ng
yy <- min(yc) + diff(yc) * as.vector(yy)

# Filter the grid points inside the convex hull
ok <- inside.owin(xx, yy, convKoala2)  
xx <- xx[ok]
yy <- yy[ok]
zxy <- as.data.frame(list(x=xx, y=yy))  # Store valid grid points

zz_sf <- st_as_sf(zxy, coords = c("x", "y"), crs = 4326)  # 4326 = WGS84
zz_proj_sf <- st_transform(zz_sf, crs = "+proj=eqc +lat_ts=60 +lon_0=0 +datum=WGS84")
print(zz_proj_sf)
zz_proj_coords <- st_coordinates(zz_proj_sf)
print(head(zz_proj_coords))

# Calculate the grid cell sizes (deltax and deltay)
deltax <- max(diff(zz_proj_coords[, 1]))
deltay <- max(diff(zz_proj_coords[order(yy), 2]))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate variables of interest (distance to center and roads)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dhaut <- gdist(xh, yh, xx, yy, unit="km")  # Distance to highest koala density point (xh, yh)


# Identify which grid points are inside the parks (based on shapefile)
ok <- rep(FALSE, length(xx))
for (i in 1:length(shp)) {
  ok[inside.owin(xx, yy, shp[[i]]) == TRUE] <- TRUE
}

# Find nearest road for each grid point and extract road distance
idr <- NULL
for( i in 1:length(xx))
{ idr <-c(idr,order(gdist(xx[i],yy[i], disroute.dat$Lon,disroute.dat$Lat,unit="km"))[1]) }
valdr <- disroute.dat$D2SROADS[idr]  # Distance to the nearest road


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Perform the same calculations for Koala locations
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
khaut <- gdist(xh, yh, Koala2$x, Koala2$y, unit="km")  # Distance of each koala to (xh, yh)
kok <- rep(FALSE, length(Koala2$x))
for (i in 1:length(shp)) {
  kok[inside.owin(Koala2$x, Koala2$y, shp[[i]]) == TRUE] <- TRUE
}

# Find nearest road for each koala location
kidr <- NULL
for( i in 1:length(Koala2$x))
{ kidr <-c(kidr,order(gdist(Koala2$x[i],Koala2$y[i], disroute.dat$Lon,disroute.dat$Lat,unit="km"))[1]) }
kvaldr <- disroute.dat$D2SROADS[kidr]  # Distance to the nearest road for koalas


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Poisson models to evaluate distance effects
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Model 1: Distance to (xh, yh) only
lv1.fct <- function(param, kaut=khaut, dhaut=dhaut, deltax=deltax, deltay=deltay) {
  l1 <- sum(-exp(param[2]) * kaut) + length(kaut) * param[1]
  l2 <- (-1) * sum(exp(-exp(param[2]) * dhaut) * exp(param[1])) * deltax * deltay
  return((-1) * (l1 + l2))
}

Z1 <- nlm(lv1.fct, c(0, 0), kaut=khaut, dhaut=dhaut, deltax=deltax, deltay=deltay)

# Model 2: Distance to (xh, yh) and park effects
lv2.fct <- function(param, kaut=khaut, dhaut=dhaut, kok=kok, ok=ok, deltax=deltax, deltay=deltay) {
  pk <- rep(1, length(kaut))
  pk[kok != TRUE] <- exp(-param[3])  # Adjust for park effect
  po <- rep(1, length(dhaut))
  po[ok != TRUE] <- exp(-param[3])  # Adjust for park effect
  l1 <- sum(-exp(param[2]) * kaut + log(pk)) + length(kaut) * param[1]
  l2 <- (-1) * sum(exp(-exp(param[2]) * dhaut + log(po))) * exp(param[1]) * deltax * deltay
  return((-1) * (l1 + l2))
}

Z2 <- nlm(lv2.fct, c(0, 0, 0), kaut=khaut, dhaut=dhaut, kok=kok, ok=ok, deltax=deltax, deltay=deltay)

# Model 3: Distance to (xh, yh), park, and road effects
lv3.fct <- function(param, kaut=khaut, dhaut=dhaut, kok=kok, ok=ok, kvaldr=kvaldr, valdr=valdr, deltax=deltax, deltay=deltay) {
  pk <- rep(1, length(kaut))
  pk[kok != TRUE] <- exp(-param[3])  # Park effect
  pk <- pk * exp(-exp(param[4]) * kvaldr)  # Road effect
  po <- rep(1, length(dhaut))
  po[ok != TRUE] <- exp(-param[3])  # Park effect
  po <- po * exp(-exp(param[4]) * valdr)  # Road effect
  l1 <- sum(-exp(param[2]) * kaut + log(pk)) + length(kaut) * param[1]
  l2 <- (-1) * sum(exp(-exp(param[2]) * dhaut + log(po))) * exp(param[1]) * deltax * deltay
  return((-1) * (l1 + l2))
}

Z3 <- nlm(lv3.fct, c(0, 0, 0, -10), kaut=khaut, dhaut=dhaut, kok=kok, ok=ok, kvaldr=kvaldr, valdr=valdr, deltax=deltax, deltay=deltay)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Intensity estimates and visualization
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Estimate the koala intensity per square kilometer based on the model
intkm2 <- exp(Z3$est[1]) * 10^6
print(c("Koala intensity per km2", intkm2))

# Plot estimated koala intensity
plot(c(0,1), c(0,1), type="n")
text(0.5, 0.8, "If all effects (parks, roads) affect detection")
text(0.2, 0.5, "Estimated intensity per km2")
text(0.2, 0.4, round(intkm2))

# Map the estimated koala intensity
param <- Z3$est
po <- rep(1, length(dhaut))
po[ok != TRUE] <- exp(-param[3])  # Adjust for park effect
po <- po * exp(-exp(param[4]) * valdr)  # Adjust for road effect
lambobs <- exp(-exp(param[2]) * dhaut + log(po)) * exp(param[1]) * 10^6

plot(convKoala2, lwd=2, col=5, main="Estimated Koala Detection Intensity")
points(xx, yy, col=gray(lambobs**0.3 / max(lambobs**0.3)), pch=19)

# Map koala intensity with park protection effect
lambpres <- po * exp(param[1]) * 10^6
plot(convKoala2, lwd=2, col=5, main="Effect of Park Protection on Koala Intensity")
points(xx, yy, col=gray(lambpres**0.3 / max(lambpres**0.3)), pch=19)

# Display intensity inside and outside parks
plot(c(0,1), c(0,1), type="n")
text(0.5, 0.8, "Intensity within parks")
text(0.5, 0.7, round(unique(lambpres)[2]))
text(0.5, 0.4, "Intensity outside parks")
text(0.5, 0.3, round(unique(lambpres)[1]))