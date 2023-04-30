library(spatstat)
library(maptools)
library(sp)

bdata <- read.csv('data/apis_data_clean.csv')
load('data/BC_Covariates.Rda') 

elev <- DATA$Elevation
forest <- DATA$Forest
hfi <- DATA$HFI
dist_water <- DATA$Dist_Water

#setup window and ppp
BC_win <- DATA$Window
BC_win <- as.owin(DATA$Window)

bees_ppp <- ppp(x = bdata$decimalLongitude, # X coordinates
                y = bdata$decimalLatitude, # Y coordinates
                window = BC_win)

#check for collinearity and save into a df for table printing later
cor_df <- as.data.frame(cor.im(elev,forest,hfi,dist_water, use="pairwise.complete.obs"))
colnames(cor_df) <- c('elev','forest','HFI','dist_water')
row.names(cor_df) <- c('elev','forest','HFI','dist_water')

#create and store a collinarity matrix plot
#pairs.im(elev,forest,hfi,dist_water)
#coll_plot = recordPlot()
#replayPlot(coll_plot)

#scale elevation and dist_water
mu <- mean(DATA$Elevation)
stdev <- sd(DATA$Elevation)
DATA$Elevation_scaled <- eval.im((Elevation - mu)/stdev, DATA)
mu <- mean(DATA$Dist_Water)
stdev <- sd(DATA$Dist_Water)
DATA$Dist_Water_scaled <- eval.im((Dist_Water - mu)/stdev, DATA)

elev_scale <- DATA$Elevation_scaled
dist_water_scale <- DATA$Dist_Water_scaled

##################################
#fit a linear and quadratic model
##################################

lin_mod <- ppm(bees_ppp ~ elev_scale + forest + hfi + dist_water_scale
                  , data = DATA)
quad_mod1 <- ppm(bees_ppp ~ elev_scale + I(elev_scale^2) + forest 
                 + I(forest^2) + hfi + I(hfi^2) + dist_water_scale 
                 + I(dist_water_scale^2), data = DATA)
quad_mod1

#visualize quadratic model predictions
plot(quad_mod1,
     se = FALSE,
     superimpose = FALSE,
     log=TRUE,
     n=285,
     main = "Estimated Honey Bee Intensity:\n Quadratic model with all covariates")
#Overlay the bee locations
plot(bees_ppp,
     pch = 16,
     cex = 0.6,
     cols = "white",
     add = TRUE)
plot(bees_ppp,
     pch = 16,
     cex = 0.5,
     cols = "black",
     add = TRUE)

#compute and compare AIC of lin and quad models
lin_AIC <- AIC(lin_mod)
quad1_AIC <- AIC(quad_mod1)
deltaAIC_linquad1 <- AIC(lin_mod) - AIC(quad_mod1)

#Likelihood Ratio Test
anova(lin_mod,quad_mod1,test = 'LRT')

#quadrat test to evaluate deviation between predicted and observed 
quadrat.test(quad_mod1, nx = 2, ny = 3, method="MonteCarlo")

########################################
#new quadratic model without dist_water
#######################################

quad_mod2 <- ppm(bees_ppp ~ elev_scale + I(elev_scale^2) + forest 
                 + I(forest^2) + hfi + I(hfi^2) + I(dist_water_scale^2), data = DATA)
quad_mod2

#visualize quadratic model predictions
plot(quad_mod2,
     se = FALSE,
     superimpose = FALSE,
     log=TRUE,
     n=285,
     main = "Estimated Honey Bee Intensity:\n Quadratic model without dist_water")
#Overlay the bee locations
plot(bees_ppp,
     pch = 16,
     cex = 0.6,
     cols = "white",
     add = TRUE)
plot(bees_ppp,
     pch = 16,
     cex = 0.5,
     cols = "black",
     add = TRUE)

#compute and compare AIC of quad1 and quad2 models
quad1_AIC <- AIC(quad_mod1)
quad2_AIC <- AIC(quad_mod2)
deltaAIC_quad12 <- AIC(quad_mod1) - AIC(quad_mod2)

#Likelihood Ratio Test
anova(quad_mod1,quad_mod2,test = 'LRT')

#quadrat test to evaluate deviation between predicted and observed 
quadrat.test(quad_mod2, nx = 2, ny = 3, method="MonteCarlo")

#diagnostics of the model
diagnose.ppm(quad_mod2)

#############################
#gam model
############################

library(splines)
### fit ppp model

gam_smooth <- ppm(bees_ppp ~ bs(elev_scale,6) + bs(forest, 12) 
                  + bs(dist_water_scale,5) + bs(hfi, 6)
                  , data = DATA, use.gam = TRUE)

#Calculate the partial residuals as a function of elevation
par_res_elev <- parres(gam_smooth, "elev_scale")
#Calculate the relative intensity as a function of forest
par_res_forest <- parres(gam_smooth, "forest")
#Calculate the partial residuals as a function of dist_water
par_res_water <- parres(gam_smooth, "dist_water_scale")
#Calculate the relative intensity as a function of hfi
par_res_hfi <- parres(gam_smooth, "hfi")

#Side by side plotting
par(mfrow = c(2,2))
plot(par_res_elev,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Elevation (m)")
plot(par_res_forest,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Forest Cover")
plot(par_res_water,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Distance to Water (m)")
plot(par_res_hfi,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Human Footprint Index")

#compute and compare AIC of quad2 and gam_smooth models
quad2_AIC <- AIC(quad_mod2)
gam_smooth_AIC <- AIC(gam_smooth)
deltaAIC_quad2gamsmooth <- AIC(quad_mod2) - AIC(gam_smooth)

#Likelihood ratio test
anova(quad_mod2, gam_smooth, test = "LRT")

#Plot the gam predictions
par(mfrow = c(1,1))
plot(gam_smooth,
     se = FALSE,
     superimpose = FALSE,
     log=TRUE,
     n=285,
     main = "Estimated Honey Bee Intensity:\n GAM with all covariates")

#Overlay the bee locations
plot(bees_ppp,
     pch = 16,
     cex = 0.6,
     cols = "white",
     add = TRUE)
plot(bees_ppp,
     pch = 16,
     cex = 0.5,
     cols = "black",
     add = TRUE)

#quadrat test to evaluate deviation between predicted and observed 
quadrat.test(gam_smooth, nx = 2, ny = 3, method="MonteCarlo")

#diagnostics of the model
diagnose.ppm(gam_smooth)

##################################
#gam model with x-y coordinates
#################################

gam_xy_smooth <- ppm(bees_ppp ~ bs(elev_scale,6) + bs(forest, 12) 
                  + bs(dist_water_scale,5) + bs(hfi, 6) + bs(x,6) + bs(y,7)
                  , data = DATA, use.gam = TRUE)

#compute and compare AIC of gam_smooth and gam_smooth_xy models
gam_smooth_AIC <- AIC(gam_smooth)
gam_xy_smooth_AIC <- AIC(gam_xy_smooth)
deltaAIC_gams <- AIC(gam_smooth) - AIC(gam_xy_smooth)

#Likelihood ratio test
anova(gam_smooth, gam_xy_smooth, test = "LRT")

#Calculate the partial residuals as a function of elevation
par_res_elev <- parres(gam_xy_smooth, "elev_scale")
#Calculate the relative intensity as a function of forest
par_res_forest <- parres(gam_xy_smooth, "forest")
#Calculate the partial residuals as a function of dist_water
par_res_water <- parres(gam_xy_smooth, "dist_water_scale")
#Calculate the relative intensity as a function of hfi
par_res_hfi <- parres(gam_xy_smooth, "hfi")
#x
par_res_x <- parres(gam_xy_smooth, "x")
#y
par_res_y <- parres(gam_xy_smooth, "y")

#Side by side plotting
par(mfrow = c(3,2))
plot(par_res_elev,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Elevation (m)")
plot(par_res_forest,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Forest Cover")
plot(par_res_water,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Distance to Water (m)")
plot(par_res_hfi,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Human Footprint Index")
plot(par_res_x,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "X")
plot(par_res_y,
     legend = FALSE,
     lwd = 2,
     main = "",
     xlab = "Y")

#Plot the gam with xy predictions
par(mfrow = c(1,1))
plot(gam_xy_smooth,
     se = FALSE,
     superimpose = FALSE,
     log=TRUE,
     n=285,
     main = "Estimated Honey Bee Intensity:\n GAM with all covariates & x,y coordinates")
#Overlay the bee locations
plot(bees_ppp,
     pch = 16,
     cex = 0.6,
     cols = "white",
     add = TRUE)
plot(bees_ppp,
     pch = 16,
     cex = 0.5,
     cols = "black",
     add = TRUE)

#quadrat test to evaluate deviation between predicted and observed 
quadrat.test(gam_xy_smooth, nx = 2, ny = 3, method="MonteCarlo")

#diagnostics of the model
diagnose.ppm(gam_xy_smooth)

##############################
#Save objects to rdata file
#############################
#add AIC levels to a data-frame
aic_df <- as.data.frame(c(lin_AIC,quad1_AIC,quad2_AIC,gam_smooth_AIC,gam_xy_smooth_AIC))
colnames(aic_df) <- c('AIC')
row.names(aic_df) <- c('linear','quadratic - full','quadratic - no dist_water','gam','gam with xy')
#save objects to 
save(cor_df, gam_smooth, gam_xy_smooth, lin_mod, quad_mod1, quad_mod2, aic_df, file = 'bees.Rdata')
