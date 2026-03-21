##############################################################
##  Project: Chagas Network
##  Title: 4_PredictionModel
##  author: Katie Tseng (katie.tseng@wsu.edu)
##############################################################

################################################################################
### OVERVIEW OF STEPS FOR EACH MODEL: Fit once & predict many scenarios
################################################################################

# Resource: https://punama.github.io/BDI_INLA/

# (0) Setup

# (1) Build a prediction grid over our study area

# (2) Build A_pred, a prediction matrix for our mesh (prediction grid)

# (3) Fit the infestation/infection model once on observed data to:
      # Obtain posterior of fixed effects and spatial field
      # Compute ROC
      # And determine values for fixed effects (e.g., svi and hai) for prediction coordinates

# (4) Model validation
      # 5-fold cross-validation

# (5) Create a prediction function to explore various scenarios, where for each, we:
      # (a) build a new prediction stack with diff covariate values,
      # (b) combine with estimation stack 
      # (c) use the fitted model to compute predictions

# (6) Run scenarios and compare predictions (violinplots)

# (7) Risk maps
   
################################################################################
### DOMICILIARY T. INFESTANS INFESTATION RISK MODEL
################################################################################

# (0) Setup --------------------------------------------------------------------

## clean environment
rm(list=ls()) 
graphics.off()

## Set directory
setwd("/Users/~/Desktop/ChagasNetworkPublic")

## load data
load('data/INLAModel_ObjectsForMapping.RData')

## load formulas as objects in environment
list2env(formulas, envir = .GlobalEnv)

## set seed
set.seed(123)

# (1) Build a prediction grid over our study area ------------------------------
    
# load packages
library(sf)
library(terra)

# Convert boundary coordinates into an sf polygon
    
    ## get coordinates of the inla.nonconvex.hull
    boundary_coords <- boundary$loc[, 1:2]
    
    ## close the polygon ring
    boundary_coords_closed <- rbind(boundary_coords, boundary_coords[1, ])

    ## create polygon from closed coordinates
    boundary_poly <- st_polygon(list(boundary_coords_closed))    

    ## add coordinate reference system
    boundary_sf <- st_sfc(boundary_poly, crs = st_crs(32721))
    
# Make grid of points inside the boundary
    
    ## choose cell size in meters
    grid_pts <- st_make_grid(boundary_sf, cellsize = 10, what = "centers")
    
    ## keep only points inside boundary
    grid_pts <- st_intersection(grid_pts, boundary_sf)
    
    ## get coordinates in matrix form
    grid_coords <- st_coordinates(grid_pts)
    
    ## rename variables
    colnames(grid_coords) <- c("x","y")
    
    ## get number of grids
    n_grid <- nrow(grid_coords)
    n_grid

    
# (2) Build A_pred (remake the A matrix for prediction) ----------------------------
    
    # load packages
    library(INLA)
    
    ## prediction locations as sf in UTM 21S
    Locations_pred <- st_as_sf(
      data.frame(x = grid_coords[,"x"], y = grid_coords[,"y"]),
      coords = c("x", "y"),
      crs = st_crs(32721)
    )
    
    ## projector matrix for prediction locations
    A_pred <- inla.spde.make.A(Mesh, loc = Locations_pred)
    dim(A_pred)
    
# (3) Fit the infestation model once on observed data only ---------------------

# create generic object names for infestation specific objects
formula <- formula_inf ## y ~ svi + hai + mobility + connected + mobility * connected + f(spatial.field, model = spde)
stk_est <- stk_inf
df <- df_inf
covar <- covar_inf
resp <- resp_inf

# Fit estimation model    
model_est <- inla(
  formula,
  data = inla.stack.data(stk_est, spde = spde),
  family = "binomial",
  control.family = list(link = "logit"),
  control.predictor = list(
    link = 1,
    A = inla.stack.A(stk_est),
    compute = TRUE  
    ),
  control.compute = list(dic = TRUE, waic = TRUE)
)

# Compute ROC

    ## get indices for estimation data
    idx_obs <- inla.stack.index(stk_est, "est")$data
    obs  <- inla.stack.data(stk_est)$y[idx_obs]
    pred <- model_est$summary.fitted.values[idx_obs, "mean"]
    
    ## compute ROC
    roc_obj <- pROC::roc(obs, pred, quiet = TRUE)
    opt <- pROC::coords(roc_obj, "best", best.method = "youden", transpose = FALSE)
    
    ## summarize with auc
    roc_summary <- data.frame(
      model = paste(deparse(formula), collapse = ""),        
      auc = as.numeric(roc_obj$auc),
      threshold = opt["threshold"],
      sensitivity = opt["sensitivity"],
      specificity = opt["specificity"]
    )


# (4) Model Validation ---------------------------------------------------------
  
# Load packages
library(INLA)
library(pROC)
library(PresenceAbsence)  # for Kappa

# Set up

    ## k-fold indices
    k <- 5
    n_total <- nrow(df)
    
    ## assign obs to folds stratified by class (positive vs. negative outcome/event vs. non-event)
        
        # outcome vector (0/1)
        y <- df[, resp]
        
        # indices by class
        pos_idx <- which(y == 1)
        neg_idx <- which(y == 0)
        
        # initialize fold vector
        folds <- integer(n_total)
        
        # assign folds separately within each class (balanced within class, then shuffled)
        folds[pos_idx] <- sample(rep(1:k, length.out = length(pos_idx)))
        folds[neg_idx] <- sample(rep(1:k, length.out = length(neg_idx)))    
        
        # check
        print(table(folds, y))
        prop.table(table(folds, y), margin = 1)
        
    ## create df to store cv metrics
    cv_results <- data.frame(
      Fold = 1:k,
      Threshold = NA_real_,
      AUC = NA_real_,
      PCC = NA_real_,
      Kappa = NA_real_,
      Sensitivity = NA_real_,
      Specificity = NA_real_    
      )
    
# Loop over folds
for (i in 1:k) {
  
    ## split train/test based on folds
    train <- df[folds != i, ]
    test <- df[folds == i, ]
    
    ## coordinates
    train_coords <- cbind(train$x, train$y)
    test_coords <- cbind(test$x, test$y)
    
    ## build A matrices
    A_train <- inla.spde.make.A(mesh = Mesh, loc = train_coords)
    A_test <- inla.spde.make.A(mesh = Mesh, loc = test_coords)
    
    ## build training stack
    stk_train <- inla.stack(
      data = list(y = train[, resp]),
      A = list(1, A_train),
      effects = list(
        fixed = train[, covar],
        spatial.field = w.spde$spatial.field
      ),
      tag = "train"
    )
    
    ## build testing stack
    stk_test <- inla.stack(
      data = list(y = NA),
      A = list(1, A_test),
      effects = list(
        fixed = test[, covar],
        spatial.field = w.spde$spatial.field
      ),
      tag = "test"
    )
    
    ## combine stacks
    stk_all <- inla.stack(stk_train, stk_test)
    
    ## fit model
    model_cv <- inla(
      formula,
      data = inla.stack.data(stk_all, spde = spde),
      family = "binomial",
      control.family = list(link = "logit"),
      control.predictor = list(
        link = 1,
        A = inla.stack.A(stk_all),
        compute = TRUE
      ), 
      control.compute = list(dic = TRUE, waic = TRUE),
      verbose = FALSE
    )
    
    ## extract predicted probabilities
    index_test <- inla.stack.index(stk_all, tag = "test")$data
    pred_prob <- model_cv$summary.fitted.values[index_test, "mean"]
    
    ## observed outcomes
    obs <- test[, resp]
    
    ## compute roc to get auc
    roc_obj <- pROC::roc(obs, pred_prob)
    
    ## determine optimal threshold w/ Youden's index
    opt <- pROC::coords(roc_obj, "best", best.method = "youden", transpose = FALSE)
    threshold <- as.numeric(opt["threshold"])  # numeric scalar
    
    ## create df of observed & predicted values where each row represents one plot
    df_accuracy <- data.frame( #confusion matrix object returned by cmx() function
      ID = 1:length(obs), # plot ID
      Obs = obs,          # observed values
      Pred = pred_prob    # predicted probablities
    )
    
    ## get accuracy metrics
    metrics <- presence.absence.accuracy(df_accuracy, threshold = threshold)
    
    ## save metrics
    cv_results[i, "Threshold"]     <- metrics$threshold
    cv_results[i, "AUC"]           <- as.numeric(roc_obj$auc)
    cv_results[i, "PCC"]           <- metrics$PCC
    cv_results[i, "Kappa"]         <- metrics$Kappa
    cv_results[i, "Sensitivity"]   <- metrics$sensitivity
    cv_results[i, "Specificity"]   <- metrics$specificity

}
    
# Save table
write.csv(
  cv_results,
  file = "tables/prediction/infestation_cv_results.csv",
  row.names = FALSE
)

# Boxplot of CV performance metrics

    ## load packages
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    
    ## transform to long
    cv_long <- cv_results %>%
      select(Fold, AUC, PCC, Kappa, Sensitivity, Specificity) %>%
      pivot_longer(
        cols = -Fold,
        names_to = "Metric",
        values_to = "Value"
      )
    
    ## plot and save
    png("figures/prediction/infestation/infestation_cv_boxplot.png", width=9,height=9,units="in",res=600)
    ggplot(cv_long, aes(x = Metric, y = Value)) +
      geom_boxplot(outlier.shape = NA, fill = "grey85") +
      geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
      ylim(0, 1) +
      theme_minimal(base_size = 13) +
      ylab("Cross-validated performance") +
      xlab(NULL) +
      ggtitle("5-fold cross-validation performance metrics")
    dev.off()
    
    
# (5) Create a prediction function for scenario modeling w/ inla -----------------------

# For prediction coordinates, determine values for covariates (e.g., svi and hai)
    
    ## review covariates for infestation model
    print(formula)
    
    ## for svi
    summary(df$svi)
    var(df$svi)
    svi_med <- median(df$svi, na.rm = TRUE)
    
    ## for hai
    summary(df$hai)
    var(df$hai)
    hai_med <- median(df$hai, na.rm = TRUE)

# Create function
predict_scenario <- function(mobility, connected) {
      
    ## make df of covariate data replicating rows n_grid times to match A_pred
    fixed_pred_grid <- data.frame(
      svi = rep(svi_med, n_grid),
      hai = rep(hai_med, n_grid),
      mobility = rep(mobility, n_grid),
      connected = rep(connected, n_grid)
    )
    
    ## stack data
    stk_pred <- inla.stack(
      data = list(y = NA),
      A = list(1, A_pred),
      effects = list(
        fixed = fixed_pred_grid,
        spatial.field = w.spde$spatial.field
      ),
      tag = "pred"
    )
    
    ## combine inf and pred data in stack
    stk_full <- inla.stack(stk_est, stk_pred)
    
    ## predict w/ inla
    model_temp <- inla(
      formula,
      data = inla.stack.data(stk_full, spde = spde),  # extract data part of stk_full
      family = "binomial",                            # likelihood
      control.predictor = list(                       # controls how linear predictor is handled
        link = 1,                                     # uses default link for the binomial family
        A = inla.stack.A(stk_full),                   # supplies projector matrix A, which links spatial field to observed data locations
        compute = TRUE                                # tells INLA to compute posterior summaries for the linear predictor
      ),
      control.mode = list(                            # Tells INLA to reuse the posterior mode of the hyperparameters from model_inf as fixed initial values when fitting a new model, instead of re-optimizing them from scratch.
        theta = model_est$mode$theta,                 # IOW, assumes spatial structure is already known (INLA does not have to refit spatial structure for each prediction run)!
        restart = FALSE                               # This speeds up computation and keeps hyperparameters fixed to ensure comparatibility across scenarios
      ),
      quantiles = c(0.025, 0.5, 0.975)
    )
    
    ## get predictions
    idx <- inla.stack.index(stk_full, "pred")$data
    preds <- model_temp$summary.fitted.values[idx, c("mean", "0.025quant", "0.5quant", "0.975quant")]
    
    ## save predictions as large objects
    saveRDS(preds, file = sprintf("data/prediction/infestation/preds_mobility%s_connected%s.rds", mobility, connected))
    
}    
    

# (6) Run scenarios and compare distribution of predictions ----------------------------

# Run predict_scenario: all dfs should be identical since ROC is based on estimation data
predict_scenario(mobility = 0, connected = 0)
predict_scenario(mobility = 1, connected = 0) 
predict_scenario(mobility = 0, connected = 1) 
predict_scenario(mobility = 1, connected = 1) 

# Load prediction results     
predictions_scenario1 <- readRDS("data/prediction/infestation/preds_mobility0_connected0.rds")
predictions_scenario2 <- readRDS("data/prediction/infestation/preds_mobility1_connected0.rds")
predictions_scenario3 <- readRDS("data/prediction/infestation/preds_mobility0_connected1.rds")
predictions_scenario4 <- readRDS("data/prediction/infestation/preds_mobility1_connected1.rds")

# Create df combining mean and median predictions 
predictions <- data.frame(
    mean = c( 
      predictions_scenario1$mean,
      predictions_scenario2$mean,
      predictions_scenario3$mean,
      predictions_scenario4$mean
    ),
    median = c(
      predictions_scenario1$`0.5quant`,
      predictions_scenario2$`0.5quant`,
      predictions_scenario3$`0.5quant`,
      predictions_scenario4$`0.5quant`
    ),
    scenario = factor(
      rep(c("m0c0", "m1c0", "m0c1", "m1c1"),
          times = c(
            length(predictions_scenario1$mean),
            length(predictions_scenario2$mean),
            length(predictions_scenario3$mean),
            length(predictions_scenario4$mean)
          ))
  )
)

# Transform predictions to long
predictions_long <- reshape(
  predictions,
  varying = c("mean", "median"),  # columns to reshape
  v.names = "value",              # name of new column for values
  timevar = "statistic",          # name of new column for the variable type
  times = c("mean", "median"),    # values to put in 'statistic'
  direction = "long"
)

# Cleanup
str(predictions_long)
predictions_long$id <- NULL  # remove id column
predictions_long$statistic <- factor(predictions_long$statistic)

## reorder factor levels
predictions_long$scenario <- factor(
  predictions_long$scenario,
  levels = c("m0c0", "m0c1", "m1c0", "m1c1")
)

# Boxplot of predicted infestation risk
    
    ## save
    png("figures/prediction/infestation/infestation_predictions_boxplot.png", width=12,height=9,units="in",res=600)
    
    ggplot(predictions_long, aes(x = scenario, y = value, fill = statistic)) +
      geom_boxplot(position = position_dodge(width = 0.75)) +
      scale_x_discrete(labels = c(
        "m0c0" = "Non-mobile / \nSocially unconnected (isolated)",
        "m0c1" = "Non-mobile / \nSocially connected",
        "m1c0" = "Mobile/ \nSocially unconnected (isolated)",
        "m1c1" = "Mobile & \nSocially connected"
      )) +
      labs(x = NULL, y = "Predicted domiciliary T. infestans infestation risk") +
      theme_minimal()
    
    dev.off() 
    
# # Violinplot w/ boxplot overlay of predicted infestation risk
# 
#     ## save
#     png("figures/prediction/infestation/infestation_predictions_violinplot.png",
#         width = 12, height = 9, units = "in", res = 600)
#     
#     ggplot(predictions_long, aes(x = scenario, y = value, fill = statistic)) +
#         geom_violin() +
#         scale_x_discrete(labels = c(
#           "m0c0" = "Non-mobile / Unconnected",
#           "m1c0" = "Mobile only",
#           "m0c1" = "Connected only",
#           "m1c1" = "Mobile + Connected"
#         )) +
#         # coord_cartesian(ylim = c(0, 1)) +
#         labs(x = NULL, y = "Predicted domiciliary T. infestans infestation risk") +
#         theme_minimal()
#     
#     dev.off()

    
# # Violinplot of predicted infestation risk
# 
#     ## Create predictions at observed locations to be added as dots on the violinplot
#     
#         # subset to median only
#         pred_violin <- subset(predictions_long, statistic == "median")
#     
#         # build an A matrix for observed household coords
#         obs_coords <- as.matrix(df[, c("x", "y")])
#         A_obs <- inla.spde.make.A(mesh = Mesh, loc = obs_coords)
#         n_obs <- nrow(df)
#         
#         # function for predicting at observed locations
#         predict_scenario_obs <- function(mobility, connected, tag) {
#           
#             ## set seed for reproducibility
#             set.seed(1234)
#             
#             ## make df of covariate data starting w/ one row of data
#             fixed_obs <- data.frame(
#               svi = rep(svi_med, n_obs),
#               hai = rep(hai_med, n_obs),
#               mobility = rep(mobility, n_obs),
#               connected = rep(connected, n_obs)
#             )
#             
#             ## stack data
#             stk_obs_pred <- inla.stack(
#               data = list(y = NA),
#               A = list(1, A_obs),
#               effects = list(
#                 fixed = fixed_obs,
#                 spatial.field = w.spde$spatial.field
#               ),
#               tag = tag
#             )
#             
#             ## combine
#             stk_full <- inla.stack(stk_est, stk_obs_pred)
#             
#             ## predict w/ inla
#             model_temp <- inla(
#               formula,
#               data = inla.stack.data(stk_full, spde = spde),  # extract data part of stk_full
#               family = "binomial",                            # likelihood
#               control.predictor = list(                       # controls how linear predictor is handled
#                 link = 1,                                     # uses default link for the binomial family
#                 A = inla.stack.A(stk_full),                   # supplies projector matrix A, which links spatial field to observed data locations
#                 compute = TRUE                                # tells INLA to compute posterior summaries for the linear predictor
#               ),
#               control.mode = list(                            # Tells INLA to reuse the posterior mode of the hyperparameters from model_inf as fixed initial values when fitting a new model, instead of re-optimizing them from scratch.
#                 theta = model_est$mode$theta,                 # IOW, assumes spatial structure is already known (INLA does not have to refit spatial structure for each prediction run)!
#                 restart = FALSE                               # This speeds up computation and keeps hyperparameters fixed to ensure comparatibility across scenarios
#               ),
#               quantiles = c(0.025, 0.5, 0.975)
#             )
#             
#             ## get predictions
#             idx <- inla.stack.index(stk_full, tag)$data
#             data.frame(
#               scenario = tag,
#               value = model_temp$summary.fitted.values[idx, "0.5quant"]# median for dots
#             )
#           }
#     
#         # run scenarios and observed-location dots
#         obs_dots <- rbind(
#           predict_scenario_obs(mobility = 0, connected = 0, tag = "m0c0"),
#           predict_scenario_obs(mobility = 1, connected = 0, tag = "m1c0"),
#           predict_scenario_obs(mobility = 0, connected = 1, tag = "m0c1"),
#           predict_scenario_obs(mobility = 1, connected = 1, tag = "m1c1")
#         )
#         
#         # ensure factor
#         obs_dots$scenario <- factor(
#           obs_dots$scenario,
#           levels = levels(pred_violin$scenario)
#         )
#    
#     ## Plot violin + dots (observed locations)
#         
#         # save
#         png("figures/prediction/infestation/infestation_predictions_violinplot(gray).png", width=16,height=9,units="in",res=600)
#       
#         ggplot() +
#           geom_violin(
#             data = pred_violin,
#             aes(x = scenario, y = value),
#             trim = FALSE,
#             fill = "grey80",
#             color = "grey40",
#             alpha = 0.9
#           ) +
#           geom_jitter(
#             data = obs_dots,
#             aes(x = scenario, y = value),
#             position = position_jitter(seed = 123, width = 0.12, height = 0),
#             alpha = 0.3,
#             size = 0.6
#           ) +
#           scale_x_discrete(labels = c(
#             "m0c0" = "Non-mobile / \nSocially unconnected (isolated)",
#             "m0c1" = "Non-mobile / \nSocially connected",
#             "m1c0" = "Mobile / \nSocially unconnected (isolated)",
#             "m1c1" = "Mobile & \nSocially connected"
#           )) +
#           labs(x = NULL, y = "Predicted domiciliary T. infestans infestation risk") +
#           theme_minimal()
#       
#         dev.off()   
         
        
    ## Plot violin by mobility + dots (observed locations)
        
        # color violin plot by infected vector abundance
        pred_violin$mobility <- ifelse(grepl("m1", pred_violin$scenario),
                                    "Mobile",
                                    "Non-mobile")   
    
        # save    
        png("figures/prediction/infestation/infestation_predictions_violinplot(color).png", width=16,height=9,units="in",res=600)
        
        ggplot() +
          geom_violin(
            data = pred_violin,
            aes(x = scenario, y = value, fill = mobility),
            trim = FALSE,
            alpha = 0.7,
            color = "grey30",
          ) +
          geom_jitter(
            data = obs_dots,
            aes(x = scenario, y = value),
            position = position_jitter(seed = 123, width = 0.12, height = 0),
            alpha = 0.3,
            size = 0.6
          ) +
          scale_fill_manual(
            name = "Mobility",
            values = c(
              "Non-mobile" = "#9ecae1",
              "Mobile" = "#08519c"        
              )
          ) +
          scale_color_manual(
            name = "Mobility",
            values = c(
              "Non-mobile" = "#9ecae1",
              "Mobile" = "#08519c"
            )
          ) +
          scale_x_discrete(labels = c(
            "m0c0" = "Non-mobile / \nSocially unconnected (isolated)",
            "m0c1" = "Non-mobile / \nSocially connected",
            "m1c0" = "Mobile / \nSocially unconnected (isolated)",
            "m1c1" = "Mobile & \nSocially connected"
          )) +
          labs(x = NULL, y = "Predicted domiciliary T. infestans infestation risk") +
          theme_minimal()
        
        dev.off()     
     
# (7) Risk maps ----------------------------------------------------------------
    
# Bind predictions to grid_coords

    ## review matrix of coordinates from prediction grid
    head(grid_coords)

    ## function to attach coords
    attach_coords <- function(preds, grid_coords) {
      cbind(data.frame(x = grid_coords[, "x"],
                       y = grid_coords[, "y"]),
            preds)
    }

    ## combine for each scenario
    df_scenario1 <- attach_coords(predictions_scenario1, grid_coords)
    df_scenario2 <- attach_coords(predictions_scenario2, grid_coords)
    df_scenario3 <- attach_coords(predictions_scenario3, grid_coords)
    df_scenario4 <- attach_coords(predictions_scenario4, grid_coords)

    ## check
    head(df_scenario1)

# Convert to raster objects w/ terra
    
    ## load package
    library(terra)
    
    ## create template raster
    r_template <- rast(
      xmin = min(grid_coords[, "x"]),
      xmax = max(grid_coords[, "x"]),
      ymin = min(grid_coords[, "y"]),
      ymax = max(grid_coords[, "y"]),
      resolution = 10, # same as grid cell size
      crs = "EPSG:32721"
    )
    
    ## function to rasterize predictions
    to_raster <- function(df, value_col) {
      
      # convert df to SpatVector with attributes
      v <- vect(
        df,
        geom = c("x", "y"),
        crs = crs(r_template)
      )
      
      # rasterize using attribute column
      r <- rasterize(
        v,
        r_template,
        field = value_col, 
        fun = "mean"
      )
      return(r)
    }
    
    # compute 95% CI width
    df_scenario1$CrI95_width <- df_scenario1$`0.975quant` - df_scenario1$`0.025quant`
    
    ## rasterize
    r_scenario1 <- to_raster(df_scenario1, "0.5quant")
    r_scenario2 <- to_raster(df_scenario2, "0.5quant")
    r_scenario3 <- to_raster(df_scenario3, "0.5quant")
    r_scenario4 <- to_raster(df_scenario4, "0.5quant")
    r_scenario1_CrI95_width <- to_raster(df_scenario1, "CrI95_width")
    
    ## save raster as GeoTIFF (tagged image file format) for import to QGIS
    writeRaster(r_scenario1, "figures/maps/forQGIS/mobility0_connected0.tif", overwrite=TRUE) # overwrite replaces existing file w/ same name
    writeRaster(r_scenario2, "figures/maps/forQGIS/mobility1_connected0.tif", overwrite=TRUE)
    writeRaster(r_scenario3, "figures/maps/forQGIS/mobility0_connected1.tif", overwrite=TRUE)
    writeRaster(r_scenario4, "figures/maps/forQGIS/mobility1_connected1.tif", overwrite=TRUE)
    writeRaster(r_scenario1_CrI95_width, "figures/maps/forQGIS/mobility0_connected0_CrI95width.tif", overwrite = TRUE)
        ## INTERPRETATION ##
        # width of the 95% posterior CrI of predicted infestation risk at each grid cell
        # dark/small IQR: model is confident about the predicted risk (posterior risk is more tightly constrained/predictions are more robust to model uncertainty)
        # light/large IQR: model is uncertain about the predicted risk there

################################################################################
### CHILD T. CRUZI INFECTION RISK MODEL (MOVER*DEGREE)
################################################################################

# (0) Setup --------------------------------------------------------------------

# Clean environment
rm(list = setdiff(ls(), c("boundary_sf", "grid_coords", "n_grid", "A_pred", "roc_summary")))
graphics.off()

## Set directory
setwd("/Users/~/Desktop/ChagasNetworkPublic")

# Load data
load('data/INLAModel_ObjectsForMapping.RData')

# Load formulas as objects in environment
list2env(formulas, envir = .GlobalEnv)

# set.seed
set.seed(123)

# (1) & (2) Skip! We prediction grid and A_pred from infestation prediction model code --------------------------------------------------------------

# (3) Fit the infection model once on observed data only -----------------------

# create generic object names for infestation specific objects
formula <- formula_infect1 ## y ~ age + mom.serology + svi + hai + iva + mover + degree + mover*degree + f(spatial.field, model = spde)
stk_est <- stk_child
df <- df_child
covar <- covar_child
resp <- resp_child

# Fit estimation model    
model_est <- inla(
  formula,
  data = inla.stack.data(stk_est, spde = spde),
  family = "binomial",
  control.family = list(link = "logit"),
  control.predictor = list(
    link = 1,
    A = inla.stack.A(stk_est),
    compute = TRUE  
  ),
  control.compute = list(dic = TRUE, waic = TRUE)
)

# Compute ROC

    ## get indices for estimation data
    idx_obs <- inla.stack.index(stk_est, "est")$data
    obs  <- inla.stack.data(stk_est)$y[idx_obs]
    pred <- model_est$summary.fitted.values[idx_obs, "mean"]
    
    ## compute ROC
    roc_obj <- pROC::roc(obs, pred, quiet = TRUE)
    opt <- pROC::coords(roc_obj, "best", best.method = "youden", transpose = FALSE)
    
    ## summarize with auc
    roc_summary <- rbind(
      roc_summary,
      data.frame(
        model = paste(deparse(formula), collapse = ""),        
        auc = as.numeric(roc_obj$auc),
        threshold = opt["threshold"],
        sensitivity = opt["sensitivity"],
        specificity = opt["specificity"]
      )
    )

# (4) Model Validation ---------------------------------------------------------

# Load packages
library(INLA)
library(pROC)
library(PresenceAbsence)  # for Kappa

# Set up

    ## k-fold indices
    k <- 5
    n_total <- nrow(df)
    
    ## assign obs to folds stratified by class (positive vs. negative outcome/event vs. non-event)
        
        # outcome vector (0/1)
        y <- df[, resp]
        
        # indices by class
        pos_idx <- which(y == 1)
        neg_idx <- which(y == 0)
        
        # initialize fold vector
        folds <- integer(n_total)
        
        # assign folds separately within each class (balanced within class, then shuffled)
        folds[pos_idx] <- sample(rep(1:k, length.out = length(pos_idx)))
        folds[neg_idx] <- sample(rep(1:k, length.out = length(neg_idx)))    
        
        # check
        print(table(folds, y))
        prop.table(table(folds, y), margin = 1)
    
    ## create df to store cv metrics
    cv_results <- data.frame(
      Fold = 1:k,
      Threshold = NA_real_,
      AUC = NA_real_,
      PCC = NA_real_,
      Kappa = NA_real_,
      Sensitivity = NA_real_,
      Specificity = NA_real_
    )

# Loop over folds
for (i in 1:k) {
  
    ## split train/test based on folds
    train <- df[folds != i, ]
    test <- df[folds == i, ]
    
    ## coordinates
    train_coords <- cbind(train$x, train$y)
    test_coords <- cbind(test$x, test$y)
    
    ## build A matrices
    A_train <- inla.spde.make.A(mesh = Mesh, loc = train_coords)
    A_test <- inla.spde.make.A(mesh = Mesh, loc = test_coords)
    
    ## build training stack
    stk_train <- inla.stack(
      data = list(y = train[, resp]),
      A = list(1, A_train),
      effects = list(
        fixed = train[, covar],
        spatial.field = w.spde$spatial.field
      ),
      tag = "train"
    )
    
    ## build testing stack
    stk_test <- inla.stack(
      data = list(y = NA),
      A = list(1, A_test),
      effects = list(
        fixed = test[, covar],
        spatial.field = w.spde$spatial.field
      ),
      tag = "test"
    )
    
    ## combine stacks
    stk_all <- inla.stack(stk_train, stk_test)
    
    ## fit model
    model_cv <- inla(
      formula,
      data = inla.stack.data(stk_all, spde = spde),
      family = "binomial",
      control.family = list(link = "logit"),
      control.predictor = list(
        link = 1,
        A = inla.stack.A(stk_all),
        compute = TRUE
      ), 
      control.compute = list(dic = TRUE, waic = TRUE),
      verbose = FALSE
    )
    
    ## extract predicted probabilities
    index_test <- inla.stack.index(stk_all, tag = "test")$data
    pred_prob <- model_cv$summary.fitted.values[index_test, "mean"]
    
    ## observed outcomes
    obs <- test[, resp]
    
    ## compute roc to get auc
    roc_obj <- pROC::roc(obs, pred_prob)
    
    ## determine optimal threshold w/ Youden's index
    opt <- pROC::coords(roc_obj, "best", best.method = "youden", transpose = FALSE)
    threshold <- as.numeric(opt["threshold"])  # numeric scalar
    
    ## create df of observed & predicted values where each row represents one plot
    df_accuracy <- data.frame( #confusion matrix object returned by cmx() function
      ID = 1:length(obs), # plot ID
      Obs = obs,          # observed values
      Pred = pred_prob    # predicted probablities
    )
    
    ## get accuracy metrics
    metrics <- presence.absence.accuracy(df_accuracy, threshold = threshold)
    
    ## save metrics
    cv_results[i, "Threshold"]     <- metrics$threshold
    cv_results[i, "AUC"]           <- as.numeric(roc_obj$auc)
    cv_results[i, "PCC"]           <- metrics$PCC
    cv_results[i, "Kappa"]         <- metrics$Kappa
    cv_results[i, "Sensitivity"]   <- metrics$sensitivity
    cv_results[i, "Specificity"]   <- metrics$specificity
    
}

# Save table
write.csv(
  cv_results,
  file = "tables/prediction/infection1_cv_results.csv",
  row.names = FALSE
)

# Boxplot of CV performance metrics

    ## load packages
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    
    ## transform to long
    cv_long <- cv_results %>%
      select(Fold, AUC, PCC, Kappa, Sensitivity, Specificity) %>%
      pivot_longer(
        cols = -Fold,
        names_to = "Metric",
        values_to = "Value"
      )
    
    ## plot and save
    png("figures/prediction/infection1/infection1_cv_boxplot.png", width=9,height=9,units="in",res=600)
    ggplot(cv_long, aes(x = Metric, y = Value)) +
      geom_boxplot(outlier.shape = NA, fill = "grey85") +
      geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
      ylim(0, 1) +
      theme_minimal(base_size = 13) +
      ylab("Cross-validated performance") +
      xlab(NULL) +
      ggtitle("5-fold cross-validation performance metrics")
    dev.off()


# (5) Create a prediction function for scenario modeling w/ inla -----------------------

# For prediction coordinates, determine values for covariates (e.g., svi and hai)

    ## review covariates for infestation model
    print(formula)
    
    ## for age
    summary(df$age)
    var(df$age)
    age_med <- median(df$age, na.rm = TRUE)
    
    ## for mom.serology: draw mom.serology as Bernoulli(p)
    p_mom1 <- mean(df$mom.serology == 1, na.rm = TRUE)   # P(mom.serology=1)
        
    ## for svi
    summary(df$svi)
    var(df$svi)
    svi_med <- median(df$svi, na.rm = TRUE)
    
    ## for hai
    summary(df$hai)
    var(df$hai)
    hai_med <- median(df$hai, na.rm = TRUE)
    
    ## for iva
    summary(df$iva)
    var(df$iva)
    iva_med <- median(df$iva, na.rm = TRUE) #median=0
    
    ## for degree
    summary(df$degree)
    var(df$degree)
    degree_med <- median(df$degree, na.rm = TRUE) #median=1

# # Prediction formula without intercept: baseline risk is absorbed by covariates + spatial field
# formula_pred <- y ~ -1 + age + mom.serology + svi + hai + iva + mover + degree + mover * 
#   degree + f(spatial.field, model = spde)

# Create function

    ## set.seed to make random draw of mom.serology reproducible
    set.seed(1234)
    
predict_scenario <- function(mover, degree) {
    
    ## make df of covariate data starting w/ one row of data
    fixed_pred_grid <- data.frame(
      age = rep(age_med, n_grid),
      mom.serology = rbinom(n_grid, size = 1, prob = p_mom1),
      svi = rep(svi_med, n_grid),
      hai = rep(hai_med, n_grid),
      iva = rep(iva_med, n_grid),
      mover = rep(mover, n_grid),
      degree = rep(degree, n_grid)
    )
    
    ## stack data
    stk_pred <- inla.stack(
      data = list(y = NA),
      A = list(1, A_pred),
      effects = list(
        fixed = fixed_pred_grid,
        spatial.field = w.spde$spatial.field
      ),
      tag = "pred"
    )
    
    ## combine inf and pred data in stack
    stk_full <- inla.stack(stk_est, stk_pred)
    
    ## predict w/ inla
    model_temp <- inla(
      formula,
      data = inla.stack.data(stk_full, spde = spde),  # extract data part of stk_full
      family = "binomial",                            # likelihood
      control.predictor = list(                       # controls how linear predictor is handled
        link = 1,                                     # uses default link for the binomial family
        A = inla.stack.A(stk_full),                   # supplies projector matrix A, which links spatial field to observed data locations
        compute = TRUE                                # tells INLA to compute posterior summaries for the linear predictor
      ),
      control.mode = list(                            # Tells INLA to reuse the posterior mode of the hyperparameters from model_inf as fixed initial values when fitting a new model, instead of re-optimizing them from scratch.
        theta = model_est$mode$theta,                 # IOW, assumes spatial structure is already known (INLA does not have to refit spatial structure for each prediction run)!
        restart = FALSE                               # This speeds up computation and keeps hyperparameters fixed to ensure comparatibility across scenarios
      ),
      quantiles = c(0.025, 0.5, 0.975)
    )
    
    ## get predictions...
    idx <- inla.stack.index(stk_full, "pred")$data
    pred_summ <- model_temp$summary.fitted.values[idx,
                                                  c("mean","0.025quant","0.5quant","0.975quant")]
    
    ## ...along with mom.serology
    preds <- data.frame(
      mom.serology = fixed_pred_grid$mom.serology,
      mean  = pred_summ[, "mean"],
      `0.025quant`  = pred_summ[, "0.025quant"],
      `0.5quant`    = pred_summ[, "0.5quant"],
      `0.975quant`  = pred_summ[, "0.975quant"],
      check.names = FALSE   # <--- prevents R from renaming columns
    )
    
    ## save predictions as large objects
    saveRDS(preds, file = sprintf("data/prediction/infection1/preds_mover%s_degree%s.rds", mover, degree))

}    

# (6) Run scenarios and compare distribution of predictions ----------------------------

# Run predict_scenario
predict_scenario(mover = 0, degree = 0)
predict_scenario(mover = 1, degree = 0) 
predict_scenario(mover = 0, degree = degree_med) 
predict_scenario(mover = 1, degree = degree_med) 

# Load prediction results     
predictions_scenario1 <- readRDS("data/prediction/infection1/preds_mover0_degree0.rds")
predictions_scenario2 <- readRDS("data/prediction/infection1/preds_mover1_degree0.rds")
predictions_scenario3 <- readRDS("data/prediction/infection1/preds_mover0_degree1.rds")
predictions_scenario4 <- readRDS("data/prediction/infection1/preds_mover1_degree1.rds")

# Create df combining mean and median predictions 
predictions <- data.frame(
  mean = c(
    predictions_scenario1$mean,
    predictions_scenario2$mean,
    predictions_scenario3$mean,
    predictions_scenario4$mean
  ),
  median = c(
    predictions_scenario1$`0.5quant`,
    predictions_scenario2$`0.5quant`,
    predictions_scenario3$`0.5quant`,
    predictions_scenario4$`0.5quant`
  ),
  mom.serology = c(
    predictions_scenario1$mom.serology,
    predictions_scenario2$mom.serology,
    predictions_scenario3$mom.serology,
    predictions_scenario4$mom.serology
  ),
  scenario = factor(
    rep(c("m0d0", "m1d0", "m0d1", "m1d1"),
        times = c(
          length(predictions_scenario1$mean),
          length(predictions_scenario2$mean),
          length(predictions_scenario3$mean),
          length(predictions_scenario4$mean)
        ))
  )
)

# Transform predictions to long
predictions_long <- predictions %>%
  tidyr::pivot_longer(
    cols = c("mean", "median"),
    names_to = "statistic",
    values_to = "value"
  )

# Cleanup
str(predictions_long)
   
     ## to factor
    predictions_long$statistic <- factor(predictions_long$statistic)

    ## reorder factor levels
    predictions_long$scenario <- factor(
      predictions_long$scenario,
      levels = c("m0d0", "m0d1", "m1d0", "m1d1")
    )
    
    ## to factor
    predictions_long$mom.serology <- factor(
      predictions_long$mom.serology,
      levels = c(0, 1),
      labels = c("Seronegative", "Seropositive")
    )
    

# Boxplot of predicted infection risk
    
    ## save
    png("figures/prediction/infection1/infection1_predictions_boxplot.png", width=12,height=9,units="in",res=600)
    
    ggplot(predictions_long, aes(x = scenario, y = value, fill = statistic)) +
      geom_boxplot(position = position_dodge(width = 0.75)) +
      scale_x_discrete(labels = c(
        "m0d0" = "Non-mover / \nSocially isolated (degree = 0)",
        "m0d1" = "Non-mover / \nSocially connected (median degree = 1)",
        "m1d0" = "Mover / \nSocially isolated (degree = 0)",
        "m1d1" = "Mover & \nSocially connected (median degree = 1)"
      )) +
      labs(x = NULL, y = "Predicted child T. cruzi infection risk") +
      theme_minimal()
    
    dev.off() 
    
# Violinplot of predicted infection risk w/ dots of observed locations

    ## Create predictions at observed locations to be added as dots on the violinplot
    
        # subset to median only
        pred_violin <- subset(predictions_long, statistic == "median")
    
        # build an A matrix for observed household coords
        obs_coords <- as.matrix(df[, c("x", "y")])
        A_obs <- inla.spde.make.A(mesh = Mesh, loc = obs_coords)
        n_obs <- nrow(df)
        
        # function for predicting at observed locations
        predict_scenario_obs <- function(mover, degree, tag) {
          
            ## set seed for reproducibility
            set.seed(1234)
            
            ## make df of covariate data starting w/ one row of data
            fixed_obs <- data.frame(
              age = rep(age_med, n_obs),
              mom.serology = rbinom(n_obs, size = 1, prob = p_mom1),
              svi = rep(svi_med, n_obs),
              hai = rep(hai_med, n_obs),
              iva = rep(iva_med, n_obs),
              mover = rep(mover, n_obs),
              degree = rep(degree, n_obs)
            )
  
            ## stack data
            stk_obs_pred <- inla.stack(
              data = list(y = NA),
              A = list(1, A_obs),
              effects = list(
                fixed = fixed_obs,
                spatial.field = w.spde$spatial.field
              ),
              tag = tag
            )
            
            ## combine
            stk_full <- inla.stack(stk_est, stk_obs_pred)
            
            ## predict w/ inla
            model_temp <- inla(
              formula,
              data = inla.stack.data(stk_full, spde = spde),  # extract data part of stk_full
              family = "binomial",                            # likelihood
              control.predictor = list(                       # controls how linear predictor is handled
                link = 1,                                     # uses default link for the binomial family
                A = inla.stack.A(stk_full),                   # supplies projector matrix A, which links spatial field to observed data locations
                compute = TRUE                                # tells INLA to compute posterior summaries for the linear predictor
              ),
              control.mode = list(                            # Tells INLA to reuse the posterior mode of the hyperparameters from model_inf as fixed initial values when fitting a new model, instead of re-optimizing them from scratch.
                theta = model_est$mode$theta,                 # IOW, assumes spatial structure is already known (INLA does not have to refit spatial structure for each prediction run)!
                restart = FALSE                               # This speeds up computation and keeps hyperparameters fixed to ensure comparatibility across scenarios
              ),
              quantiles = c(0.025, 0.5, 0.975)
            )
            
            ## get predictions
            idx <- inla.stack.index(stk_full, tag)$data
            data.frame(
              scenario = tag,
              value = model_temp$summary.fitted.values[idx, "0.5quant"], # median for dots
              mom.serology = fixed_obs$mom.serology # return the realized mom.serology draw
            )
          }
    
        # run scenarios and observed-location dots
        obs_dots <- rbind(
          predict_scenario_obs(mover = 0, degree = 0, tag = "m0d0"),
          predict_scenario_obs(mover = 1, degree = 0, tag = "m1d0"),
          predict_scenario_obs(mover = 0, degree = 1, tag = "m0d1"),
          predict_scenario_obs(mover = 1, degree = 1, tag = "m1d1")
        )
        
        # ensure scenario is factor
        obs_dots$scenario <- factor(
          obs_dots$scenario,
          levels = levels(pred_violin$scenario)
        )
        
        # ensure mom.serology is factor
        obs_dots$mom.serology <- factor(
          obs_dots$mom.serology,
          levels = c(0, 1),
          labels = c("Seronegative", "Seropositive")
        )

        
    ## Plot violin by mover + dots (observed locations)
        
        # save    
        png("figures/prediction/infection1/infection1_predictions_violinplot(color).png", width=16,height=9,units="in",res=600)
        
        ggplot() +
          geom_violin(
            data = pred_violin,
            aes(x = scenario, y = value, fill = mover),
            trim = FALSE,
            alpha = 0.7,
            color = "grey30",
          ) +
          geom_jitter(
            data = obs_dots,
            aes(x = scenario, y = value, color = mom.serology),
            shape = 21,
            fill = NA,
            stroke = 0.5,
            position = position_jitter(seed = 123, width = 0.12, height = 0),
            alpha = 0.6,
            size = 0.7
          ) +
          scale_fill_manual(
            name = "Mobility",
            values = c(
              "Non-mover" = "#9ecae1",
              "Mover" = "#08519c"        
              )
          ) +
          scale_color_manual(
            name = "Maternal seropositivity",
            values = c(
              "Seronegative" = "grey60",
              "Seropositive" = "#de2d26"
            )
          ) +
          scale_x_discrete(labels = c(
            "m0d0" = "Socially isolated (degree = 0)",
            "m0d1" = "Socially connected (median degree = 1)",
            "m1d0" = "Socially isolated (degree = 0)",
            "m1d1" = "Socially connected (median degree = 1)"
          )) +
          labs(x = NULL, y = "Predicted child T. cruzi infection risk") +
          theme_minimal() +
          guides(
            color = guide_legend(
              override.aes = list(
                size = 3,
                stroke = 1,
                alpha = 1
              )
            )
          )
        
        dev.off()     

    # ## Plot violin by mom.serology (+ mover and dots (observed locations))
    #     
    #     # color violin plot by mobility
    #     pred_violin$mover <- ifelse(grepl("m1", pred_violin$scenario),
    #                                 "Mover",
    #                                 "Non-mover")   
    # 
    #     # create separate dfs for pred_violin by mom.serology
    #     pred_violin0 <- pred_violin[pred_violin$mom.serology == "Seronegative",]
    #     pred_violin1 <- pred_violin[pred_violin$mom.serology == "Seropositive",]
    #     
    #     # create separate dfs for obs_dots by mom.serology
    #     obs_dots0 <- obs_dots[obs_dots$mom.serology == "Seronegative",]
    #     obs_dots1 <- obs_dots[obs_dots$mom.serology == "Seropositive",]
    # 
    #     # save    
    #     png("figures/prediction/infection1/infection1_predictions_violinplot_momseronegative.png", width=16,height=9,units="in",res=600)
    #     
    #     ggplot() +
    #       geom_violin(
    #         data = pred_violin0,
    #         aes(x = scenario, y = value, fill = mover),
    #         trim = FALSE,
    #         alpha = 0.8,
    #         color = "grey30",
    #       ) +
    #       geom_jitter(
    #         data = obs_dots0,
    #         aes(x = scenario, y = value),
    #         position = position_jitter(seed = 123, width = 0.12, height = 0),
    #         alpha = 0.2,
    #         size = 0.6
    #       ) +
    #       scale_fill_manual(
    #         name = "Mobility",
    #         values = c(
    #           "Non-mover" = "#9ecae1",
    #           "Mover" = "#08519c"        
    #           )
    #       ) +
    #       scale_color_manual(
    #         name = "Mobility",
    #         values = c(
    #           "Non-mover" = "#9ecae1",
    #           "Mover" = "#08519c"
    #         )
    #       ) +
    #       scale_x_discrete(labels = c(
    #         "m0d0" = "Socially isolated (degree = 0)",
    #         "m0d1" = "Socially connected (median degree = 1)",
    #         "m1d0" = "Socially isolated (degree = 0)",
    #         "m1d1" = "Socially connected (median degree = 1)"
    #       )) +
    #       labs(x = NULL, y = "Predicted child T. cruzi infection risk") +
    #       theme_minimal()
    #     
    #     dev.off()     
    # 
    #     
    #     # save    
    #     png("figures/prediction/infection1/infection1_predictions_violinplot_momseropositive.png", width=16,height=9,units="in",res=600)
    #     
    #     ggplot() +
    #       geom_violin(
    #         data = pred_violin1,
    #         aes(x = scenario, y = value, fill = mover),
    #         trim = FALSE,
    #         alpha = 0.8,
    #         color = "grey30",
    #       ) +
    #       geom_jitter(
    #         data = obs_dots1,
    #         aes(x = scenario, y = value),
    #         position = position_jitter(seed = 123, width = 0.12, height = 0),
    #         alpha = 0.2,
    #         size = 0.6
    #       ) +
    #       scale_fill_manual(
    #         name = "Mobility",
    #         values = c(
    #           "Non-mover" = "#9ecae1",
    #           "Mover" = "#08519c"        
    #         )
    #       ) +
    #       scale_color_manual(
    #         name = "Mobility",
    #         values = c(
    #           "Non-mover" = "#9ecae1",
    #           "Mover" = "#08519c"
    #         )
    #       ) +
    #       scale_x_discrete(labels = c(
    #         "m0d0" = "Socially isolated (degree = 0)",
    #         "m0d1" = "Socially connected (median degree = 1)",
    #         "m1d0" = "Socially isolated (degree = 0)",
    #         "m1d1" = "Socially connected (median degree = 1)"
    #       )) +
    #       labs(x = NULL, y = "Predicted child T. cruzi infection risk") +
    #       theme_minimal()
    #     
    #     dev.off()     
    #     
    # # Plot split violin
    #     
    #     # Function for split violin
    #         # Source - https://stackoverflow.com/a/45614547
    #         # Posted by jan-glx, modified by community. See post 'Timeline' for change history
    #         # Retrieved 2026-03-02, License - CC BY-SA 4.0
    #     GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
    #                                draw_group = function(self, data, ..., draw_quantiles = NULL) {
    #                                  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    #                                  grp <- data[1, "group"]
    #                                  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
    #                                  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    #                                  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    #                                  
    #                                  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    #                                    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
    #                                                                              1))
    #                                    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    #                                    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    #                                    aesthetics$alpha <- rep(1, nrow(quantiles))
    #                                    both <- cbind(quantiles, aesthetics)
    #                                    quantile_grob <- GeomPath$draw_panel(both, ...)
    #                                    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    #                                  }
    #                                  else {
    #                                    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    #                                  }
    #                                })
    #     
    #     geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
    #                                   draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
    #                                   show.legend = NA, inherit.aes = TRUE) {
    #       layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
    #             position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
    #             params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
    #     }
    #     
    #     # save
    #     png("figures/prediction/infection1/infection1_predictions_violinplot(split).png", width=16,height=9,units="in",res=600)
    #     
    #     ggplot(pred_violin, aes(x = scenario, y = value, fill = mom.serology)) +
    #       geom_split_violin() + 
    #       geom_jitter(data = obs_dots,
    #                   aes(x = scenario, y = value),
    #                   inherit.aes = FALSE,
    #                   width = 0.08,
    #                   height = 0,
    #                   size = 0.8,
    #                   alpha = 0.2,
    #                   color = "black") +
    #       scale_fill_manual(name = "Maternal serology",
    #                         values = c("Seronegative" = "#fc9272",
    #                                    "Seropositive" = "#de2d26")) +
    #       scale_x_discrete(labels = c(
    #         "m0d0" = "Non-mover / \nSocially isolated (degree = 0)",
    #         "m0d1" = "Non-mover / \nSocially connected (median degree = 1)",
    #         "m1d0" = "Mover / \nSocially isolated (degree = 0)",
    #         "m1d1" = "Mover / \nSocially connected (median degree = 1)"
    #       )) +
    #       labs(x = "",
    #            y = expression(paste("Predicted child ", italic("Trypanosoma cruzi"), " infection risk"))) +
    #       theme_minimal() +
    #       theme(legend.position = "bottom",
    #             legend.direction = "horizontal")
    #     
    #     dev.off()
     
# (7) Risk maps ----------------------------------------------------------------

# Bind predictions to grid_coords

    ## review matrix of coordinates from prediction grid
    head(grid_coords)
    
    ## function to attach coords
    attach_coords <- function(preds, grid_coords) {
      cbind(data.frame(x = grid_coords[, "x"],
                       y = grid_coords[, "y"]),
            preds)
    }
    
    ## combine for each scenario
    df_scenario1 <- attach_coords(predictions_scenario1, grid_coords)
    df_scenario2 <- attach_coords(predictions_scenario2, grid_coords)
    df_scenario3 <- attach_coords(predictions_scenario3, grid_coords)
    df_scenario4 <- attach_coords(predictions_scenario4, grid_coords)
    
    ## check
    head(df_scenario1)

# Convert to raster objects w/ terra
    
    ## load package
    library(terra)
    
    ## create template raster
    r_template <- rast(
      xmin = min(grid_coords[, "x"]),
      xmax = max(grid_coords[, "x"]),
      ymin = min(grid_coords[, "y"]),
      ymax = max(grid_coords[, "y"]),
      resolution = 10, # same as grid cell size
      crs = "EPSG:32721"
    )
    
    ## function to rasterize predictions
    to_raster <- function(df, value_col) {
      
      # convert df to SpatVector with attributes
      v <- vect(
        df,
        geom = c("x", "y"),
        crs = crs(r_template)
      )
      
      # rasterize using attribute column
      r <- rasterize(
        v,
        r_template,
        field = value_col, 
        fun = "mean"
      )
      return(r)
    }
    
    ## compute 95% CI width (baseline risk uncertainty)
    df_scenario1$CrI95_width <- df_scenario1$`0.975quant` - df_scenario1$`0.025quant`
    
    ## rasterize
    r_scenario1 <- to_raster(df_scenario1, "0.5quant")
    r_scenario2 <- to_raster(df_scenario2, "0.5quant")
    r_scenario3 <- to_raster(df_scenario3, "0.5quant")
    r_scenario4 <- to_raster(df_scenario4, "0.5quant")
    r_scenario1_CrI95_width <- to_raster(df_scenario1, "CrI95_width")
    
    ## save raster as GeoTIFF (tagged image file format) for import to QGIS
    writeRaster(r_scenario1, "figures/maps/forQGIS/mover0_degree0.tif", overwrite=TRUE) # overwrite replaces existing file w/ same name
    writeRaster(r_scenario2, "figures/maps/forQGIS/mover1_degree0.tif", overwrite=TRUE)
    writeRaster(r_scenario3, "figures/maps/forQGIS/mover0_degree1.tif", overwrite=TRUE)
    writeRaster(r_scenario4, "figures/maps/forQGIS/mover1_degree1.tif", overwrite=TRUE)
    writeRaster(r_scenario1_CrI95_width, "figures/maps/forQGIS/mover0_degree0_CrI95width.tif", overwrite=TRUE)
        ### INTERPRETATION: baseline risk uncertainty map 
        # width of the 95% posterior CrI of predicted infestation risk at each grid cell
        # dark/small range: model is confident about the predicted risk (posterior risk is more tightly constrained/predictions are more robust to model uncertainty)
        # light/large range: model is uncertain about the predicted risk there


################################################################################
### CHILD T. CRUZI INFECTION RISK MODEL (IVA*HARMONIC)
################################################################################

# (0) Setup --------------------------------------------------------------------

# Clean environment
rm(list = setdiff(ls(), c("boundary_sf", "grid_coords", "n_grid", "A_pred", "roc_summary")))
graphics.off()

# Set directory
setwd("/Users/~/Desktop/ChagasNetworkPublic")

# Load data
load('data/INLAModel_ObjectsForMapping.RData')

# Load formulas as objects in environment
list2env(formulas, envir = .GlobalEnv)

# Set seed
set.seed(123)

# (1) & (2) Skip! --------------------------------------------------------------
## We already did this for infestation model and saved those results in setup  
# (3) Fit the infection model once on observed data only -----------------------

# create generic object names for infestation specific objects
formula <- formula_infect2 ## y ~ age + mom.serology + svi + hai + iva + mover + harmonic + iva * harmonic + f(spatial.field, model = spde)
stk_est <- stk_child
df <- df_child
covar <- covar_child
resp <- resp_child

# Fit estimation model    
model_est <- inla(
  formula,
  data = inla.stack.data(stk_est, spde = spde),
  family = "binomial",
  control.family = list(link = "logit"),
  control.predictor = list(
    link = 1,
    A = inla.stack.A(stk_est),
    compute = TRUE  
  ),
  control.compute = list(dic = TRUE, waic = TRUE)
)

# Compute ROC

    ## get indices for estimation data
    idx_obs <- inla.stack.index(stk_est, "est")$data
    obs  <- inla.stack.data(stk_est)$y[idx_obs]
    pred <- model_est$summary.fitted.values[idx_obs, "mean"]
    
    ## compute ROC
    roc_obj <- pROC::roc(obs, pred, quiet = TRUE)
    opt <- pROC::coords(roc_obj, "best", best.method = "youden", transpose = FALSE)
    
    ## summarize with auc
    roc_summary <- rbind(
      roc_summary,
      data.frame(
        model = paste(deparse(formula), collapse = ""),        
        auc = as.numeric(roc_obj$auc),
        threshold = opt["threshold"],
        sensitivity = opt["sensitivity"],
        specificity = opt["specificity"]
      )
    )


# (4) Model Validation ---------------------------------------------------------

# Load packages
library(INLA)
library(pROC)
library(PresenceAbsence)  # for Kappa

# Set up

    ## k-fold indices
    k <- 5
    n_total <- nrow(df)
    
    ## assign obs to folds stratified by class (positive vs. negative outcome/event vs. non-event)
        
        # outcome vector (0/1)
        y <- df[, resp]
        
        # indices by class
        pos_idx <- which(y == 1)
        neg_idx <- which(y == 0)
        
        # initialize fold vector
        folds <- integer(n_total)
        
        # assign folds separately within each class (balanced within class, then shuffled)
        folds[pos_idx] <- sample(rep(1:k, length.out = length(pos_idx)))
        folds[neg_idx] <- sample(rep(1:k, length.out = length(neg_idx)))    
        
        # check
        print(table(folds, y))
        prop.table(table(folds, y), margin = 1)
    
    ## create df to store cv metrics
    cv_results <- data.frame(
      Fold = 1:k,
      Threshold = NA_real_,
      AUC = NA_real_,
      PCC = NA_real_,
      Kappa = NA_real_,
      Sensitivity = NA_real_,
      Specificity = NA_real_
    )

# Loop over folds
for (i in 1:k) {
  
    ## split train/test based on folds
    train <- df[folds != i, ]
    test <- df[folds == i, ]
    
    ## coordinates
    train_coords <- cbind(train$x, train$y)
    test_coords <- cbind(test$x, test$y)
    
    ## build A matrices
    A_train <- inla.spde.make.A(mesh = Mesh, loc = train_coords)
    A_test <- inla.spde.make.A(mesh = Mesh, loc = test_coords)
    
    ## build training stack
    stk_train <- inla.stack(
      data = list(y = train[, resp]),
      A = list(1, A_train),
      effects = list(
        fixed = train[, covar],
        spatial.field = w.spde$spatial.field
      ),
      tag = "train"
    )
    
    ## build testing stack
    stk_test <- inla.stack(
      data = list(y = NA),
      A = list(1, A_test),
      effects = list(
        fixed = test[, covar],
        spatial.field = w.spde$spatial.field
      ),
      tag = "test"
    )
    
    ## combine stacks
    stk_all <- inla.stack(stk_train, stk_test)
    
    ## fit model
    model_cv <- inla(
      formula,
      data = inla.stack.data(stk_all, spde = spde),
      family = "binomial",
      control.family = list(link = "logit"),
      control.predictor = list(
        link = 1,
        A = inla.stack.A(stk_all),
        compute = TRUE
      ), 
      control.compute = list(dic = TRUE, waic = TRUE),
      verbose = FALSE
    )
    
    ## extract predicted probabilities
    index_test <- inla.stack.index(stk_all, tag = "test")$data
    pred_prob <- model_cv$summary.fitted.values[index_test, "mean"]
    
    ## observed outcomes
    obs <- test[, resp]
    
    ## compute roc to get auc
    roc_obj <- pROC::roc(obs, pred_prob)
    
    ## determine optimal threshold w/ Youden's index
    opt <- pROC::coords(roc_obj, "best", best.method = "youden", transpose = FALSE)
    threshold <- as.numeric(opt["threshold"])  # numeric scalar
    
    ## create df of observed & predicted values where each row represents one plot
    df_accuracy <- data.frame( #confusion matrix object returned by cmx() function
      ID = 1:length(obs), # plot ID
      Obs = obs,          # observed values
      Pred = pred_prob    # predicted probablities
    )
    
    ## get accuracy metrics
    metrics <- presence.absence.accuracy(df_accuracy, threshold = threshold)
    
    ## save metrics
    cv_results[i, "Threshold"]     <- metrics$threshold
    cv_results[i, "AUC"]           <- as.numeric(roc_obj$auc)
    cv_results[i, "PCC"]           <- metrics$PCC
    cv_results[i, "Kappa"]         <- metrics$Kappa
    cv_results[i, "Sensitivity"]   <- metrics$sensitivity
    cv_results[i, "Specificity"]   <- metrics$specificity
    
}

# Save table
write.csv(
  cv_results,
  file = "tables/prediction/infection2_cv_results.csv",
  row.names = FALSE
)

# Boxplot of CV performance metrics

    ## load packages
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    
    ## transform to long
    cv_long <- cv_results %>%
      select(Fold, AUC, PCC, Kappa, Sensitivity, Specificity) %>%
      pivot_longer(
        cols = -Fold,
        names_to = "Metric",
        values_to = "Value"
      )
    
    ## plot and save
    png("figures/prediction/infection2/infection2_cv_boxplot.png", width=9,height=9,units="in",res=600)
    ggplot(cv_long, aes(x = Metric, y = Value)) +
      geom_boxplot(outlier.shape = NA, fill = "grey85") +
      geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
      ylim(0, 1) +
      theme_minimal(base_size = 13) +
      ylab("Cross-validated performance") +
      xlab(NULL) +
      ggtitle("5-fold cross-validation performance metrics")
    dev.off()


# (5) Create a prediction function for scenario modeling w/ inla -----------------------

# For prediction coordinates, determine values for covariates (e.g., svi and hai)

    ## review covariates for infestation model
    print(formula)
    
    ## for age
    summary(df$age)
    var(df$age)
    age_med <- median(df$age, na.rm = TRUE)
    
    ## for mom.serology: draw mom.serology as Bernoulli(p)
    p_mom1 <- mean(df$mom.serology == 1, na.rm = TRUE)   # P(mom.serology=1)
        
    ## for svi
    summary(df$svi)
    var(df$svi)
    svi_med <- median(df$svi, na.rm = TRUE)
    
    ## for hai
    summary(df$hai)
    var(df$hai)
    hai_med <- median(df$hai, na.rm = TRUE)
    
    ## for iva: highly zero-inflated (75th% = 0), scenario-based predictions should fix IVA at 0 for baseline and at 1 (discrete nonzero value) rather than using the mean
    summary(df$iva)
    var(df$iva)
    iva1 <- 1
    
    ## for mover: draw mom as Bernoulli(p)
    p_mover1 <- mean(df$mover == 1, na.rm = TRUE)   # P(mom.serology=1)
    
    ## for harmonic, let's explore both 50th (median) and 75th percentiles
        
        # median (50th%)
        summary(df$harmonic)
        var(df$harmonic)
        harmonic_q2 <- median(df$harmonic, na.rm = TRUE) #median=1
        
        # 3rd quartile (75th%)
        harmonic_q3 <- as.numeric(quantile(df$harmonic, 0.75, na.rm = TRUE))
        
# Create function
predict_scenario <- function(iva, harmonic) {
  
    ## set.seed to make random draw of mom.serology reproducible
    set.seed(1234)
    
    ## make df of covariate data starting w/ one row of data
    fixed_pred_grid <- data.frame(
      age = rep(age_med, n_grid),  # replicate rows n_grid times to match A_pred
      mom.serology = rbinom(n_grid, size = 1, prob = p_mom1), # draw a 0/1 for each grid cell
      svi = rep(svi_med, n_grid),
      hai = rep(hai_med, n_grid),
      iva = rep(iva, n_grid),
      mover = rbinom(n_grid, size = 1, prob = p_mover1),
      harmonic = rep(harmonic, n_grid)
    )

    ## stack data
    stk_pred <- inla.stack(
      data = list(y = NA),
      A = list(1, A_pred),
      effects = list(
        fixed = fixed_pred_grid,
        spatial.field = w.spde$spatial.field
      ),
      tag = "pred"
    )
    
    ## combine inf and pred data in stack
    stk_full <- inla.stack(stk_est, stk_pred)
    
    ## predict w/ inla
    model_temp <- inla(
      formula,
      data = inla.stack.data(stk_full, spde = spde),  # extract data part of stk_full
      family = "binomial",                            # likelihood
      control.predictor = list(                       # controls how linear predictor is handled
        link = 1,                                     # uses default link for the binomial family
        A = inla.stack.A(stk_full),                   # supplies projector matrix A, which links spatial field to observed data locations
        compute = TRUE                                # tells INLA to compute posterior summaries for the linear predictor
      ),
      control.mode = list(                            # Tells INLA to reuse the posterior mode of the hyperparameters from model_inf as fixed initial values when fitting a new model, instead of re-optimizing them from scratch.
        theta = model_est$mode$theta,                 # IOW, assumes spatial structure is already known (INLA does not have to refit spatial structure for each prediction run)!
        restart = FALSE                               # This speeds up computation and keeps hyperparameters fixed to ensure comparatibility across scenarios
      ),
      quantiles = c(0.025, 0.5, 0.975)
    )
    
    ## get predictions...
    idx <- inla.stack.index(stk_full, "pred")$data
    pred_summ <- model_temp$summary.fitted.values[idx,
                                                  c("mean","0.025quant","0.5quant","0.975quant")]
    
    ## ...along with mom.serology
    preds <- data.frame(
      mom.serology = fixed_pred_grid$mom.serology,
      mean  = pred_summ[, "mean"],
      `0.025quant`  = pred_summ[, "0.025quant"],
      `0.5quant`    = pred_summ[, "0.5quant"],
      `0.975quant`  = pred_summ[, "0.975quant"],
      check.names = FALSE   # <--- prevents R from renaming columns
    )
    
    ## save predictions as large objects
    saveRDS(preds, file = sprintf("data/prediction/infection2/preds_iva%s_harmonic%.2f.rds", iva, harmonic)) #rounds harmonic value

}    


# (6) Run scenarios and compare distribution of predictions ----------------------------

# Run predict_scenario
predict_scenario(iva = 0, harmonic = 0)
predict_scenario(iva = 0, harmonic = harmonic_q2) 
predict_scenario(iva = 0, harmonic = harmonic_q3) 
predict_scenario(iva = iva1, harmonic = 0) 
predict_scenario(iva = iva1, harmonic = harmonic_q2)
predict_scenario(iva = iva1, harmonic = harmonic_q3)

# Load prediction results: you may have to go to file and manually update the names     
predictions_scenario1 <- readRDS("data/prediction/infection2/preds_iva0_harmonic0.00.rds")
predictions_scenario2 <- readRDS("data/prediction/infection2/preds_iva0_harmonic1.67.rds")
predictions_scenario3 <- readRDS("data/prediction/infection2/preds_iva0_harmonic3.75.rds")
predictions_scenario4 <- readRDS("data/prediction/infection2/preds_iva1_harmonic0.00.rds")
predictions_scenario5 <- readRDS("data/prediction/infection2/preds_iva1_harmonic1.67.rds")
predictions_scenario6 <- readRDS("data/prediction/infection2/preds_iva1_harmonic3.75.rds")

saveRDS(predictions_scenario1, file = sprintf("data/prediction/infection2/preds_iva0_harmonic0.00.rds"))
saveRDS(predictions_scenario2, file = sprintf("data/prediction/infection2/preds_iva0_harmonic1.67.rds"))
saveRDS(predictions_scenario3, file = sprintf("data/prediction/infection2/preds_iva0_harmonic3.75.rds"))
saveRDS(predictions_scenario4, file = sprintf("data/prediction/infection2/preds_iva1_harmonic0.00.rds"))
saveRDS(predictions_scenario5, file = sprintf("data/prediction/infection2/preds_iva1_harmonic1.67.rds"))
saveRDS(predictions_scenario6, file = sprintf("data/prediction/infection2/preds_iva1_harmonic3.75.rds"))

# Create df combining mean and median predictions 
library(dplyr)

predictions <- bind_rows(
  predictions_scenario1 %>% mutate(scenario = "iva0_harm0"),
  predictions_scenario2 %>% mutate(scenario = "iva0_harmq2"),
  predictions_scenario3 %>% mutate(scenario = "iva0_harmq3"),
  predictions_scenario4 %>% mutate(scenario = "iva1_harm0"),
  predictions_scenario5 %>% mutate(scenario = "iva1_harmq2"),
  predictions_scenario6 %>% mutate(scenario = "iva1_harmq3")
)

# Transform predictions to long
predictions_long <- predictions %>%
  tidyr::pivot_longer(
    cols = c("mean", "0.5quant"),
    names_to = "statistic",
    values_to = "value"
  )

# Cleanup
str(predictions_long)

    ## to factor
    predictions_long$statistic <- factor(predictions_long$statistic)
    
    ## reorder factor levels
    predictions_long$scenario <- factor(
      predictions_long$scenario,
      levels = c("iva0_harm0", "iva0_harmq2", "iva0_harmq3", "iva1_harm0", "iva1_harmq2", "iva1_harmq3")
    )
    
    ## to factor
    predictions_long$mom.serology <- factor(
      predictions_long$mom.serology,
      levels = c(0, 1),
      labels = c("Seronegative", "Seropositive")
    )

# Boxplot of predicted infection risk

    ## load package
    library(ggplot2)
    
    ## save
    png("figures/prediction/infection2/infection2_predictions_boxplot.png", width=16,height=9,units="in",res=600)
    
    ggplot(predictions_long, aes(x = scenario, y = value, fill = statistic)) +
      geom_boxplot(position = position_dodge(width = 0.75)) +
      scale_x_discrete(labels = c(
        "iva0_harm0" = "No infected vector / \nSocially isolated (harmonic centrality = 0)",
        "iva0_harmq2" = "No infected vector / \nMedian harmonic centrality",
        "iva0_harmq3" = "No infected vector / \nHigh harmonic centrality (75th perc.)",
        "iva1_harm0" = "≥1 infected vector / \nSocially isolated (harmonic centrality = 0)",
        "iva1_harmq2" = "≥1 infected vector  / \nMedian harmonic centrality",
        "iva1_harmq3" = "≥1 infected vector / \nHigh harmonic centrality (75th perc.)"
      )) +
      labs(x = NULL, y = "Predicted child T. cruzi infection risk") +
      theme_minimal()
    
    dev.off() 

# Violinplot of predicted infection risk

    ## Create predictions at observed locations to be added as dots on the violinplot
    
        # subset to median only
        pred_violin <- subset(predictions_long, statistic == "0.5quant")
    
        # build an A matrix for observed household coords
        obs_coords <- as.matrix(df[, c("x", "y")])
        A_obs <- inla.spde.make.A(mesh = Mesh, loc = obs_coords)
        n_obs <- nrow(df)
        
        # function for predicting at observed locations
        predict_scenario_obs <- function(iva, harmonic, tag) {
          
            ## set seed for reproducibility
            set.seed(1234)
            
            ## make df of covariate data starting w/ one row of data
            fixed_obs <- data.frame(
              age = rep(age_med, n_obs),
              mom.serology = rbinom(n_obs, size = 1, prob = p_mom1),
              svi = rep(svi_med, n_obs),
              hai = rep(hai_med, n_obs),
              iva = rep(iva, n_obs),
              mover = rbinom(n_obs, size = 1, prob = p_mover1),
              harmonic = rep(harmonic, n_obs)
            )
          
            ## stack data
            stk_obs_pred <- inla.stack(
              data = list(y = NA),
              A = list(1, A_obs),
              effects = list(
                fixed = fixed_obs,
                spatial.field = w.spde$spatial.field
              ),
              tag = tag
            )
            
            ## combine
            stk_full <- inla.stack(stk_est, stk_obs_pred)
            
            ## predict w/ inla
            model_temp <- inla(
              formula,
              data = inla.stack.data(stk_full, spde = spde),  # extract data part of stk_full
              family = "binomial",                            # likelihood
              control.predictor = list(                       # controls how linear predictor is handled
                link = 1,                                     # uses default link for the binomial family
                A = inla.stack.A(stk_full),                   # supplies projector matrix A, which links spatial field to observed data locations
                compute = TRUE                                # tells INLA to compute posterior summaries for the linear predictor
              ),
              control.mode = list(                            # Tells INLA to reuse the posterior mode of the hyperparameters from model_inf as fixed initial values when fitting a new model, instead of re-optimizing them from scratch.
                theta = model_est$mode$theta,                 # IOW, assumes spatial structure is already known (INLA does not have to refit spatial structure for each prediction run)!
                restart = FALSE                               # This speeds up computation and keeps hyperparameters fixed to ensure comparatibility across scenarios
              ),
              quantiles = c(0.025, 0.5, 0.975)
            )

            
            ## get predictions
            idx <- inla.stack.index(stk_full, tag)$data
            data.frame(
              scenario = tag,
              value = model_temp$summary.fitted.values[idx, "0.5quant"], # median for dots
              mom.serology = fixed_obs$mom.serology # return the realized mom.serology draw
            )
        }
        
        # run scenarios and observed-location dots
        obs_dots <- rbind(
          predict_scenario_obs(iva = 0, harmonic = 0,           tag = "iva0_harm0"),
          predict_scenario_obs(iva = 0, harmonic = harmonic_q2, tag = "iva0_harmq2"),
          predict_scenario_obs(iva = 0, harmonic = harmonic_q3, tag = "iva0_harmq3"),
          predict_scenario_obs(iva = 1, harmonic = 0,           tag = "iva1_harm0"),
          predict_scenario_obs(iva = 1, harmonic = harmonic_q2, tag = "iva1_harmq2"),
          predict_scenario_obs(iva = 1, harmonic = harmonic_q3, tag = "iva1_harmq3")
        )
        
        # ensure scenario is factor
        obs_dots$scenario <- factor(
          obs_dots$scenario,
          levels = levels(pred_violin$scenario)
        )
        
        # ensure mom.serology is factor
        obs_dots$mom.serology <- factor(
          obs_dots$mom.serology,
          levels = c(0, 1),
          labels = c("Seronegative", "Seropositive")
        )
       
    ## Plot violin by mover + dots (observed locations)
        
        # save
        png("figures/prediction/infection2/infection2_predictions_violinplot(color).png", width=16,height=9,units="in",res=600)
        
        ggplot() +
          geom_violin(
            data = pred_violin,
            aes(x = scenario, y = value, fill = iva),
            trim = FALSE,
            alpha = 0.7,
            color = "grey30",
          ) +
          geom_jitter(
            data = obs_dots,
            aes(x = scenario, y = value, color = mom.serology),
            shape = 21,
            fill = NA,
            stroke = 0.5,
            position = position_jitter(seed = 123, width = 0.12, height = 0),
            alpha = 0.6,
            size = 0.7
          ) +
          scale_fill_manual(
            name = "Infected vector abundance",
            values = c(
              "No infected vector" = "#9ecae1",
              "≥1 infected vector" = "#08519c"        
              )
          ) +

          scale_color_manual(
            name = "Maternal seropositivity",
            values = c(
              "Seronegative" = "grey60",
              "Seropositive" = "#de2d26"
            )
          ) +
          scale_x_discrete(labels = c(
            "iva0_harm0"  = "Socially isolated \n(harmonic centrality = 0)",
            "iva0_harmq2" = "Median harmonic centrality",
            "iva0_harmq3" = "High harmonic centrality \n(75th perc.)",
            "iva1_harm0"  = "Socially isolated \n(harmonic centrality = 0)",
            "iva1_harmq2" = "Median harmonic centrality",
            "iva1_harmq3" = "High harmonic centrality \n(75th perc.)"
          )) +
          labs(x = NULL, y = "Predicted child T. cruzi infection risk") +
          theme_minimal() +
          guides(
            color = guide_legend(
              override.aes = list(
                size = 3,
                stroke = 1,
                alpha = 1
              )
            )
          )
        
        dev.off()     

    # ## Plot violin: one for mom seronegative observations and one for mom seropositive
    # 
    #     # color violin plot by infected vector abundance
    #     pred_violin$iva <- ifelse(grepl("iva1", pred_violin$scenario),
    #                               "≥1 infected vector",
    #                               "No infected vector")                          
    #     
    #     # create separate dfs for pred_violin
    #     pred_violin0 <- pred_violin[pred_violin$mom.serology == "Seronegative",]
    #     pred_violin1 <- pred_violin[pred_violin$mom.serology == "Seropositive",]
    #     
    #     # create separate dfs for obs_dots
    #     obs_dots0 <- obs_dots[obs_dots$mom.serology == "Seronegative",]
    #     obs_dots1 <- obs_dots[obs_dots$mom.serology == "Seropositive",]
    #     
    #     
    #     # save
    #     png("figures/prediction/infection2/infection2_predictions_violinplot(momseronegative).png", width=16,height=9,units="in",res=600)
    #     
    #     ggplot() +
    #       geom_violin(
    #         data = pred_violin0,
    #         aes(x = scenario, y = value, fill = iva),
    #         trim = FALSE,
    #         alpha = 0.8,
    #         color = "grey30",
    #       ) +
    #       geom_jitter(
    #         data = obs_dots0,
    #         aes(x = scenario, y = value),
    #         position = position_jitter(seed = 123, width = 0.12, height = 0),
    #         alpha = 0.2,
    #         size = 0.6
    #       ) +
    #       scale_fill_manual(
    #         name = "Infected vector abundance",
    #         values = c(
    #           "No infected vector" = "#9ecae1",
    #           "≥1 infected vector" = "#08519c"        
    #         )
    #       ) +
    #       scale_color_manual(
    #         name = "Infected vector abundance",
    #         values = c(
    #           "No infected vector" = "#9ecae1",
    #           "≥1 infected vector" = "#08519c"
    #         )
    #       ) +
    #       scale_x_discrete(labels = c(
    #         "iva0_harm0"  = "Socially isolated \n(harmonic centrality = 0)",
    #         "iva0_harmq2" = "Median harmonic centrality",
    #         "iva0_harmq3" = "High harmonic centrality \n(75th perc.)",
    #         "iva1_harm0"  = "Socially isolated \n(harmonic centrality = 0)",
    #         "iva1_harmq2" = "Median harmonic centrality",
    #         "iva1_harmq3" = "High harmonic centrality \n(75th perc.)"
    #       )) +
    #       labs(x = NULL, y = "Predicted child T. cruzi infection risk") +
    #       theme_minimal()
    #     
    #     dev.off()    
    #     
    #     # save
    #     png("figures/prediction/infection2/infection2_predictions_violinplot(momseropositive).png", width=16,height=9,units="in",res=600)
    #     
    #     ggplot() +
    #       geom_violin(
    #         data = pred_violin1,
    #         aes(x = scenario, y = value, fill = iva),
    #         trim = FALSE,
    #         alpha = 0.8,
    #         color = "grey30",
    #       ) +
    #       geom_jitter(
    #         data = obs_dots1,
    #         aes(x = scenario, y = value),
    #         position = position_jitter(seed = 123, width = 0.12, height = 0),
    #         alpha = 0.2,
    #         size = 0.6
    #       ) +
    #       scale_fill_manual(
    #         name = "Infected vector abundance",
    #         values = c(
    #           "No infected vector" = "#9ecae1",
    #           "≥1 infected vector" = "#08519c"        
    #         )
    #       ) +
    #       scale_color_manual(
    #         name = "Infected vector abundance",
    #         values = c(
    #           "No infected vector" = "#9ecae1",
    #           "≥1 infected vector" = "#08519c"
    #         )
    #       ) +
    #       scale_x_discrete(labels = c(
    #         "iva0_harm0"  = "Socially isolated \n(harmonic centrality = 0)",
    #         "iva0_harmq2" = "Median harmonic centrality",
    #         "iva0_harmq3" = "High harmonic centrality \n(75th perc.)",
    #         "iva1_harm0"  = "Socially isolated \n(harmonic centrality = 0)",
    #         "iva1_harmq2" = "Median harmonic centrality",
    #         "iva1_harmq3" = "High harmonic centrality \n(75th perc.)"
    #       )) +
    #       labs(x = NULL, y = "Predicted child T. cruzi infection risk") +
    #       theme_minimal()
    #     
    #     dev.off()   
    #     
    # ## Plot split violin
    #     
    #     # Function for split violin
    #     # Source - https://stackoverflow.com/a/45614547
    #     # Posted by jan-glx, modified by community. See post 'Timeline' for change history
    #     # Retrieved 2026-03-02, License - CC BY-SA 4.0
    #     GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
    #                                draw_group = function(self, data, ..., draw_quantiles = NULL) {
    #                                  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    #                                  grp <- data[1, "group"]
    #                                  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
    #                                  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    #                                  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    #                                  
    #                                  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    #                                    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
    #                                                                              1))
    #                                    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    #                                    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    #                                    aesthetics$alpha <- rep(1, nrow(quantiles))
    #                                    both <- cbind(quantiles, aesthetics)
    #                                    quantile_grob <- GeomPath$draw_panel(both, ...)
    #                                    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    #                                  }
    #                                  else {
    #                                    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    #                                  }
    #                                })
    #     
    #     geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
    #                                   draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
    #                                   show.legend = NA, inherit.aes = TRUE) {
    #       layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
    #             position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
    #             params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
    #     }
    #     
    #     # save
    #     png("figures/prediction/infection2/infection2_predictions_violinplot(split).png", width=16,height=9,units="in",res=600)
    #     
    #     ggplot(pred_violin, aes(x = scenario, y = value, fill = mom.serology)) +
    #       geom_split_violin() + 
    #       geom_jitter(data = obs_dots,
    #                   aes(x = scenario, y = value),
    #                   inherit.aes = FALSE,
    #                   width = 0.08,
    #                   height = 0,
    #                   size = 0.8,
    #                   alpha = 0.2,
    #                   color = "black") +
    #       scale_fill_manual(name = "Maternal serology",
    #                         values = c("Seronegative" = "#fc9272",
    #                                    "Seropositive" = "#de2d26")) +
    #       scale_x_discrete(labels = c(
    #         "iva0_harm0"  = "No infected vector / \nSocially isolated (harmonic centrality = 0)",
    #         "iva0_harmq2" = "No infected vector / \nMedian harmonic centrality",
    #         "iva0_harmq3" = "No infected vector / \nHigh harmonic centrality (75th perc.)",
    #         "iva1_harm0"  = "≥1 infected vector / \nSocially isolated (harmonic centrality = 0)",
    #         "iva1_harmq2" = "≥1 infected vector / \nMedian harmonic centrality",
    #         "iva1_harmq3" = "≥1 infected vector / \nHigh harmonic centrality (75th perc.)"
    #       )) +
    #       labs(x = "",
    #            y = expression(paste("Predicted child ", italic("Trypanosoma cruzi"), " infection risk"))) +
    #       theme_minimal() +
    #       theme(legend.position = "bottom",
    #             legend.direction = "horizontal")
    #     
    #     dev.off()
        
                
# (7) Risk maps ----------------------------------------------------------------

# Bind predictions to grid_coords

    ## review matrix of coordinates from prediction grid
    head(grid_coords)
    
    ## function to attach coords
    attach_coords <- function(preds, grid_coords) {
      cbind(data.frame(x = grid_coords[, "x"],
                       y = grid_coords[, "y"]),
            preds)
    }
    
    ## combine for each scenario
    df_scenario1 <- attach_coords(predictions_scenario1, grid_coords) # iva0_harm0
    df_scenario2 <- attach_coords(predictions_scenario2, grid_coords) # iva0_harmq2
    df_scenario3 <- attach_coords(predictions_scenario3, grid_coords) # iva0_harmq3
    df_scenario4 <- attach_coords(predictions_scenario4, grid_coords) # iva1_harm0
    df_scenario5 <- attach_coords(predictions_scenario5, grid_coords) # iva1_harmq2
    df_scenario6 <- attach_coords(predictions_scenario6, grid_coords) # iva1_harmq3
    
    ## check
    head(df_scenario1)

# Convert to raster objects w/ terra
    
    ## load package
    library(terra)
    
    ## create template raster
    r_template <- rast(
      xmin = min(grid_coords[, "x"]),
      xmax = max(grid_coords[, "x"]),
      ymin = min(grid_coords[, "y"]),
      ymax = max(grid_coords[, "y"]),
      resolution = 10, # same as grid cell size
      crs = "EPSG:32721"
    )
    
    ## function to rasterize predictions
    to_raster <- function(df, value_col) {
      
      # convert df to SpatVector with attributes
      v <- vect(
        df,
        geom = c("x", "y"),
        crs = crs(r_template)
      )
      
      # rasterize using attribute column
      r <- rasterize(
        v,
        r_template,
        field = value_col, 
        fun = "mean"
      )
      return(r)
    }
    
    ## compute 95% CI width (baseline risk uncertainty)
    df_scenario1$CrI95_width <- df_scenario1$`0.975quant` - df_scenario1$`0.025quant`
    
    ## rasterize: this will contain grid cell values, extent, resolution, and CRS
    r_scenario1 <- to_raster(df_scenario1, "0.5quant") # iva0_harm0
    r_scenario2 <- to_raster(df_scenario2, "0.5quant") # iva0_harmq2
    r_scenario3 <- to_raster(df_scenario3, "0.5quant") # iva0_harmq3
    r_scenario4 <- to_raster(df_scenario4, "0.5quant") # iva1_harm0
    r_scenario5 <- to_raster(df_scenario5, "0.5quant") # iva1_harmq2
    r_scenario6 <- to_raster(df_scenario6, "0.5quant") # iva1_harmq3
    r_scenario1_CrI95_width <- to_raster(df_scenario1, "CrI95_width")
    
    ## save raster as GeoTIFF (tagged image file format) for import to QGIS
    writeRaster(r_scenario1, "figures/maps/forQGIS/iva0_harm0.tif", overwrite=TRUE) # overwrite replaces existing file w/ same name
    writeRaster(r_scenario2, "figures/maps/forQGIS/iva0_harmq2.tif", overwrite=TRUE)
    writeRaster(r_scenario3, "figures/maps/forQGIS/iva0_harmq3.tif", overwrite=TRUE)
    writeRaster(r_scenario4, "figures/maps/forQGIS/iva1_harm0.tif", overwrite=TRUE)
    writeRaster(r_scenario5, "figures/maps/forQGIS/iva1_harmq2.tif", overwrite=TRUE)
    writeRaster(r_scenario6, "figures/maps/forQGIS/iva1_harmq3.tif", overwrite=TRUE)   
    writeRaster(r_scenario1_CrI95_width, "figures/maps/forQGIS/iva0_harm0_CrI95width.tif"
     
################################################################################
### ADDITIONAL FIGURES
################################################################################
        
# Combine cross-validation results and save table
        
    ## read in cv
    cv_infestation <- read.csv("tables/prediction/infestation_cv_results.csv")
    cv_infection1 <- read.csv("tables/prediction/infection1_cv_results.csv")
    cv_infection2 <- read.csv("tables/prediction/infection2_cv_results.csv")

    ## add model var
    cv_infestation$model <- "infestation"
    cv_infection1$model <- "infection1"        
    cv_infection2$model <- "infection2"
    
    ## combine
    cv_results <- rbind(cv_infestation, cv_infection1)
    cv_results <- rbind(cv_results, cv_infection2)
    
    ## to factor
    cv_results$model <- factor(cv_results$model, 
                                  levels = c("infestation", "infection1", "infection2"),
                                  labels = c("Domiciliary infestation risk model (household mobility x social connectivity)",
                                             "Child infection risk model (mover x household degree)",
                                             "Child infection risk model (domiciliary infected vector abundance x household harmonic centrality)")
    )
  
    ## save table
    write.csv(
      cv_results,
      file = "tables/prediction/allmodels_cv_results.csv",
      row.names = FALSE
    )

# Get mean ± SD of cv metrics
    
    ## set directory
    setwd("/Users/~/Desktop/ChagasNetworkPublic")
    
    ## read in cv results
    cv_results <- read.csv("tables/prediction/allmodels_cv_results.csv")

    ## transform to long
    cv_long <- cv_results %>%
      pivot_longer(
        cols = -c(Fold, Threshold, model),
        names_to = "Metric",
        values_to = "Value"
      )
    
    # calculate mean ± SD for each model and metric
    cv_summary <- cv_long %>%
      group_by(model, Metric) %>%
      summarise(
        Mean = mean(Value, na.rm = TRUE),
        SD = sd(Value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(`Mean ± SD` = paste0(round(Mean, 2), " ± ", round(SD, 3))) %>%
      select(model, Metric, `Mean ± SD`)
    
    # View the summarized table
    View(cv_summary)

# Boxplot CV results
    
    ## plot and save
    png("figures/prediction/allmodels_cv_boxplot.png", width=9,height=9,units="in",res=600)
    ggplot(cv_long, aes(x = Metric, y = Value, fill = model)) +
      geom_boxplot(
        alpha = 0.7,
        outlier.shape = NA, 
        position = position_dodge(width = 0.9)
        ) +
      geom_jitter(
        aes(fill = model),      # fill matches boxplot
        shape = 21,             # circle with outline
        color = "black",        # outline color
        stroke = 0.5,           # thickness of outline
        size = 2,
        alpha = 0.8,
        position = position_jitterdodge(
          jitter.width = 0.1,
          dodge.width = 0.9
        )
      ) +
      scale_fill_viridis_d(option = "D", end = 0.9) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom",
            legend.direction = "vertical") +
      labs(
        y = "Cross-validated performance",
        x = NULL,
        fill = "Model",
        color = "Model",
        title = "5-fold cross-validation performance metrics"
      )
    dev.off()

# Export Area III layers to QGIS:
    
    ## We needs the following layers from top to bottom for final maps in QGIS
        # Community outlines (lines)
        # Area III outer outline (line)
        # Prediction raster (masked to Area III - polygon)
    
    ## (1) Get polygon of Area III
        
        # load package
        library(sf)
        library(dplyr)
        
        # load shape file    
        sf_comms <- read_sf("/Users/katietseng/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/Chagas Network/GitHub/ChagasNetwork/data/Shapefiles PPI/communities.shp")
        
        # subset shape file to polygons of Area III
        areaIII <- sf_comms %>% filter(tipo == 3)
        
        # fix invalid geometries
        areaIII_valid <- st_make_valid(areaIII) %>%
          st_as_sf()
        
        # project to UTM zone 21S (EPSG:32721) - match your raster
        areaIII_utm <- st_transform(areaIII_valid, 32721)
        
        
        # dissolve all polygons into one geometry
        areaIII_polygon <- areaIII_utm %>%
          st_combine() %>%  # combine geometries but keep class sf
          st_union()        # dissolve boundaries
        
        # check
        plot(areaIII_polygon)
        
        # export for QGIS
        write_sf(
          areaIII_polygon,
          "figures/maps/forQGIS/areaIII_polygon.gpkg",
          layer = "areaIII_polygon",
          delete_layer = TRUE
        )
        dev.off()
    
    ## (2) Using Area III polygon, clip folder of rasters
        
        # load package
        library(terra)
        
        # function
        clip_folder_to_polygon <- function(
          folder,
          polygon_sf,               # your areaIII_polygon (sf) or SpatVector
          pattern = "\\.tif$",
          recursive = FALSE
          ) {
            # list files
            tif_files <- list.files(folder, pattern = pattern, full.names = TRUE, recursive = recursive) # recursive means R will go into directories/subfolders within the folder
            
            # safety check that there are .tif files in the folder
            if (length(tif_files) == 0) stop("No .tif files found in: ", folder)
            
            # convert polygon once (inherit ensures polygon is a terra SpatVector)
            poly_v <- if (inherits(polygon_sf, "SpatVector")) polygon_sf else vect(polygon_sf)
            
            for (f in tif_files) { 
              r <- rast(f)
              
              # make polygon CRS match raster CRS (important):
              poly_use <- poly_v
              if (!is.na(crs(r)) && !is.na(crs(poly_use)) && crs(r) != crs(poly_use)) {  # if (1) raster has a defined CRS, (2) polygon has a CRS, and (3) that CRS are the same, then...
                poly_use <- project(poly_use, crs(r))
              }
              
              # crop raster to bounding box
              r_crop <- crop(r, poly_use)
              
              # mask raster to exact boundary
              r_mask <- mask(r_crop, poly_use)
              
              # overwrite same filename (saves back to the same location)
              writeRaster(r_mask, f, overwrite = TRUE)
              
              message("Clipped: ", basename(f))
            }
            
            invisible(tif_files)
            }  
        
        # run function on folder
        clip_folder_to_polygon(
          folder = "figures/maps/forQGIS",
          polygon_sf = areaIII_polygon
          ) 
        
    ## (3) Get outline of Area III & communities
    
        # assume pdi_outline_projected is in a projected CRS (e.g., UTM)
        areaIII_outline <- st_transform(areaIII_polygon, crs = 32721)
        
        # check
        plot(areaIII_outline)
        
        # export for QGIS
        write_sf(
          areaIII_outline,
          "figures/maps/forQGIS/areaIII_outline.gpkg",
          layer = "areaIII_outline",
          delete_layer = TRUE
        )
        dev.off()
 