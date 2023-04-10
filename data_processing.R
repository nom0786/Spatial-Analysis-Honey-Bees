#install.packages(c("rgbif"))
library(rgbif)
library(spatstat)
library(maptools)
library(sp)
library(rgdal)

scientificName <- c("Apis mellifera Linnaeus")

# check how many records exist for BC
record_count <- occ_count(scientificName = scientificName,
                          hasCoordinate = TRUE,
                          country = "CA",
                          stateProvince = "British Columbia")

# download occurrence data for the species in BC
gbif_data_bc <- occ_data(scientificName = scientificName,
                           hasCoordinate = TRUE,
                           country = "CA",
                           stateProvince = "British Columbia",
                           limit = record_count) 

# extract the data element from the list
apis_data <- gbif_data_bc$data

# some parameters are list... need to adjust
apis_data <- apply(apis_data,2,as.character)

# save the data as CSV
write.csv(apis_data, "data/apis_data.csv", row.names = FALSE)

# read csv
apis_data <- read.csv("data/apis_data.csv")
head(apis_data)

# get the columns that matter for mapping and cleaning the occurrence data
#used this as reference 'https://www.r-bloggers.com/2021/03/downloading-and-cleaning-gbif-data-with-r/'

apis_data_clean <- apis_data[ , c("decimalLongitude", "decimalLatitude", "country", "stateProvince", "occurrenceStatus", "coordinateUncertaintyInMeters",
                                  "taxonID", "catalogNumber", "institutionCode", "eventTime", "verbatimEventDate", "collectionCode", "gbifID", "verbatimLocality",
                                  "class", "isInCluster", "stateProvince", "year", "month", "day", "eventDate", "modified", "lastInterpreted")]

# remove records of absence or zero-abundance (if any)
names(apis_data_clean)
sort(unique(apis_data_clean$individualCount))  # notice if some points correspond to zero abundance
sort(unique(apis_data_clean$occurrenceStatus))  # check for different indications of "absent", which could be in different languages
absence_rows <- which(apis_data_clean$individualCount == 0 | apis_data_clean$occurrenceStatus %in% c("absent", "Absent", "ABSENT", "ausente", "Ausente", "AUSENTE"))
length(absence_rows)
if (length(absence_rows) > 0) {
  apis_data_clean <- apis_data_clean[-absence_rows, ]
}

#in order to plot the observations the data had to be adjusted so that the 
#decimalLongitude and decimalLatitudevalues in the downloaded dataset matched the 
#format of the covariate data.

# create a SpatialPointsDataFrame object using the apis_data_clean dataset
apis_data_sp <- SpatialPointsDataFrame(coords = apis_data_clean[, c("decimalLongitude", "decimalLatitude")],
                                       data = apis_data_clean)

# assign the WGS84 CRS (EPSG:4326) to the apis_data_sp object
proj4string(apis_data_sp) <- CRS("+proj=longlat +datum=WGS84")

# transform the coordinates to match the CRS of the BC_win object
bc_crs <- "+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs"
apis_data_transformed <- spTransform(apis_data_sp, CRS(bc_crs))

# extract the transformed coordinates and assign them to the original cleaned df
apis_data_clean$decimalLongitude <- apis_data_transformed@coords[, 1]
apis_data_clean$decimalLatitude <- apis_data_transformed@coords[, 2]

head(apis_data_clean)

# save the clean/tranformed data as CSV
write.csv(apis_data_clean, "data/apis_data_clean.csv", row.names = FALSE)

# read cleaned csv
apis_data_clean <- read.csv("data/apis_data_clean.csv")
head(apis_data_clean)

#Visualise the data
plot(decimalLatitude ~ decimalLongitude,
     pch = 16,
     col = "#046C9A",
     data = apis_data_clean)

# BC covariate data
load("data/BC_Covariates.Rda")
BC_win <- DATA$Window
BC_win <- as.owin(DATA$Window)

#Visualise the window
plot(BC_win,
     main = "Observation Window")


apis_data_ppp <- ppp(x = apis_data_clean$decimalLongitude, # X coordinates
                     y = apis_data_clean$decimalLatitude, # Y coordinates
                     window = BC_win) # Observation window


plot(apis_data_ppp)

#Refine the figure
plot(apis_data_ppp, # The dataset to visualise
     col = "darkgreen", #The colour of the window
     cols = 'yellow', #The colours of the points
     pch = 18, # The plotting symbol
     main = "Bees in BC", # The title
     par(bg="grey", cex.main = 2), # Flexible modification of the graphical parameters, 
     cex = 0.6) # Turn of the legend depending on needs`
