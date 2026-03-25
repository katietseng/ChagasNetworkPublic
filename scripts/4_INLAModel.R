##############################################################
##  Project: Chagas Network
##  Title: 3_INLAModel
##  author: Katie Tseng (katie.tseng@wsu.edu)
##############################################################

################################################################################
### PREP DATA
################################################################################

# Setup

    ## clean environment
    rm(list=ls()) 
    graphics.off()
    
    ## set directory
    setwd("/Users/~/Desktop/ChagasNetworkPublic")
    
    ## load data
    load('data/CleanData.RData')

# Create empty dfs for saving GLM and GLMM results

    ## glm results
    glm_results <- data.frame(depvar = character(),
                              covar = character(),
                              level = character(),
                              estimate = numeric(),
                              POR = numeric(),
                              LL = numeric(),
                              UL = numeric(),
                              P = numeric(),
                              AIC = numeric(),
                              nobs = numeric(),
                              ci_method = character(),
                              stringsAsFactors = FALSE
    )

    ## glmm results
    glmm_results <- glm_results

# Create dfs for Infestation and IVA model (inf. vector abundance)
    
    ## gen kin network infestation variable
    df_kin$kin.occ.inf <- ifelse(is.nan(df_kin$kin.inf.prop), NA,
                                 ifelse(df_kin$kin.inf.prop > 0, 1, 0))
    
        # transform to factor
        df_kin$kin.occ.inf <- as.factor(df_kin$kin.occ.inf)
        
        # drop those with NA infestation data
        df <- df_kin[!is.na(df_kin$v.inf),]
    
    ## subset variables for modeling and check type
    df <- subset(df, select = c(hhx.id, x, y, v.inf, reinf, iva, degree, degree.cat, connected, harmonic, svi, hai, mobility, ethnicity, kin.svi.med, kin.svi.relmed, kin.occ.inf, kin.iva.med))
    
        # inspect data (e.g., is response variable numeric binary?)
        str(df)
        
        # save binary factor variables with reference group reversed
        df$connected.bin <- factor(df$connected, levels = c("Unconnected", "Connected"), labels = c("Unconnected", "Connected"))
        df$connected.binrev <- factor(df$connected, levels = c("Connected", "Unconnected"), labels = c("Connected", "Unconnected"))
        df$mobility.bin <- factor(df$mobility, levels = c("non-mover", "mover"), labels = c("non-mover", "mover"))
        df$mobility.binrev <- factor(df$mobility, levels = c("mover", "non-mover"), labels = c("mover", "non-mover"))
        
        # transform binary predictors to numeric
        df$v.inf <- as.numeric(df$v.inf)
        df$v.inf <- ifelse(is.na(df$v.inf), NA, 
                           ifelse(df$v.inf == 2, 1, 0))
        df$reinf <- as.numeric(df$reinf)
        df$reinf <- ifelse(is.na(df$reinf), NA, 
                           ifelse(df$reinf == 2, 1, 0))
        df$connected <- as.numeric(df$connected)
        df$connected <- ifelse(is.na(df$connected), NA, 
                               ifelse(df$connected == 2, 1, 0))
        df$mobility <- as.numeric(df$mobility)
        df$mobility <- ifelse(is.na(df$mobility), NA,
                              ifelse(df$mobility == 2, 1, 0))
        df$ethnicity <- as.numeric(df$ethnicity)
        df$ethnicity <- ifelse(is.na(df$ethnicity), NA, 
                               ifelse(df$ethnicity == 2, 1, 0))
        df$kin.occ.inf <- as.numeric(df$kin.occ.inf)
        df$kin.occ.inf <- ifelse(is.na(df$kin.occ.inf), NA, 
                                 ifelse(df$kin.occ.inf == 2, 1, 0))
      
    ## select model covariates
    spat_inf <- c("hhx.id","x","y") # spatial variables
    resp_inf <- c("v.inf")
    resp_reinf <- c("reinf")
    resp_iva <- c("iva")
    covar_inf <- c("degree",
               "degree.cat",
               "connected",
               "connected.bin",
               "connected.binrev",
               "harmonic",
               "ethnicity",
               "svi",
               "hai",
               "kin.svi.med",
               "kin.svi.relmed",
               "kin.occ.inf",
               "kin.iva.med",
               "mobility",
               "mobility.bin",
               "mobility.binrev")
    
    # get final dfs and drop NAs
    # df <- na.omit(df[,c(spat, resp1, resp2, resp3, covar)]) #336
    df_inf <- na.omit(df[,c(spat_inf, resp_inf, resp_reinf, covar_inf)]) #336
    df_iva <- na.omit(df[,c(spat_inf, resp_iva, covar_inf)]) #324
    
    # how many infestations and reinfestations?
    table(df_inf$v.inf) # n = 88
    table(df_inf$reinf) # n = 6
    
## Create df for Child Infection model
    
    ## get dem data
    df_child <- df_dem
    names(df_child)

    ## get additional hhx-level variables
    df_add <- subset(df_kin, select = c(hhx.id, 
                                        v.inf, 
                                        connected, degree, degree.cat, harmonic, betweenness,
                                        kin.svi.med, kin.svi.relmed,
                                        kin.iva.med, kin.inf.prop,
                                        kin.occ.inf))
    
        ## merge
    df_child <- merge(df_child, df_add, by.x = "hhx.id", all.x = TRUE)
        rm(df_add)
    
        ## checks
            
            # how many households have at least one seropositive person?
            temp <- df_child[!is.na(df_child$serology),] #1368
            temp <- aggregate(serology ~ hhx.id, data = temp, sum) #303
            temp$occ.serology <- ifelse(temp$serology>0, 1, 0)
            mean(temp$occ.serology) #71%
            
            # how many households have at least one seropositive person?
            temp <- df_child[!is.na(df_child$serology.child),]
            temp <- aggregate(serology.child ~ hhx.id, data = temp, sum) #226
            temp$occ.serology.child <- ifelse(temp$serology.child>0, 1, 0)
            mean(temp$occ.serology.child) #22%
            rm(temp)
        
    ## drop those w/ NA serology.child data
    df <- df_child[!is.na(df_child$serology.child),] #705
    
    ## subset variables for modeling and check type
    df <- subset(df, select = c(hhx.id, serology.child, degree, degree.cat, connected, harmonic, age, iva, v.inf, svi, hai, mom.serology, mover, kin.svi.med, kin.svi.relmed, kin.occ.inf, kin.iva.med))
   
        ## inspect data (is response variable numeric binary)
        str(df)
        
        ## merge spatial data
        coords <- subset(df_kin, select = c(hhx.id, x, y))
        df <- merge(df, coords, by.x = "hhx.id", all.x =TRUE)
        
        # save binary factor variables with reference group reversed
        df$connected.bin <- factor(df$connected, levels = c("Unconnected", "Connected"), labels = c("Unconnected", "Connected"))
        df$connected.binrev <- factor(df$connected, levels = c("Connected", "Unconnected"), labels = c("Connected", "Unconnected"))
        df$mover.bin <- factor(df$mover, levels = c("Non-mover", "Mover"), labels = c("Non-mover", "Mover"))
        df$mover.binrev <- factor(df$mover, levels = c("Mover", "Non-mover"), labels = c("Mover", "Non-mover"))
        
        # transform binary predictors to numeric
        df$connected <- as.numeric(df$connected)
        df$connected <- ifelse(df$connected == 2, 1, 0)
        df$mom.serology <- as.numeric(df$mom.serology)
        df$mom.serology <- ifelse(df$mom.serology == 2, 1, 0)
        df$mover <- as.numeric(df$mover)
        df$mover <- ifelse(df$mover == 2, 1, 0)
        df$kin.occ.inf <- as.numeric(df$kin.occ.inf)
        df$kin.occ.inf <- ifelse(df$kin.occ.inf == 2, 1, 0)
        
    # select model variables
    spat_child <- c("hhx.id","x","y") # spatial variables
    resp_child <- c("serology.child")
    covar_child <- c("degree",
               "degree.cat",
               "connected",
               "connected.bin",
               "connected.binrev",
               "harmonic",
               "age",
               "iva",
               "v.inf",
               "svi",
               "hai",
               "mom.serology",
               "mover",
               "mover.bin",
               "mover.binrev",
               "kin.svi.med",
               "kin.svi.relmed",
               "kin.occ.inf",
               "kin.iva.med"
    )
    
    # get final df and drop NAs
    df_child <- na.omit(df[,c(spat_child, resp_child, covar_child)]) #n=547
    
################################################################################
### BUILDING THE SPATIAL FIELD: Mesh + SPDE
################################################################################

# Mesh Construction ------------------------------------------------------------

# load library
library(INLA)

# Get sampling locations (house coordinates)
    
    ## save as dataframe
    Locations = as.data.frame(cbind(df_kin$x, df_kin$y))
    
    # assign UTM Zone 21S CRS (EPSG:32721)
    library(sf)
    Locations <- st_as_sf(Locations, coords = c("V1", "V2"),
                          crs = st_crs(32721)) # UTM Zone 21S: EPSG:32721
    
    # verify the CRS
    st_crs(Locations)

# to determine the maximum edge size of mesh triangles, let's see the distribution of distances between houses
D <- st_distance(Locations, Locations)
mean(as.matrix(D))
temp <- array(D)
hist(temp)

# change dimension names to match our matrix of houses
row.names(D) <- df_kin$hhx.id
colnames(D) <- df_kin$hhx.id

# get boundary based on coordinates of our observations/houses
boundary <- inla.nonconvex.hull(Locations, convex = -0.3) # convex parameter controls how much concavity is allowed (negative val -> tighter fit, pos val -> looser fit)
lines(boundary, add = FALSE)
points(Locations, col = "red")

# save plot w/ transparent background
png("figures/inla/mesh/mesh_boundary.png", width = 6, height = 8, units = "in", res = 600, bg = "transparent")
plot(boundary)
points(Locations, col = "red", cex = 0.7)
dev.off()

# mesh construction w/ boundary: creates meshes w/ diff triangle size
mesh <- inla.mesh.2d(boundary = boundary, max.edge = c(500, 1500), # max.edge: values denote max. allowed triangle edge lengths in region & buffer zone, respectively
                     cutoff = 10) # min. allowed distance b/w points; used to avoid building too many small triangles around clustered data locations
## set inner edge to 1/5-1/10 of the smallest spatial pattern you want to detect:
## occ of seropositive person aggregated at a scale b/w 2-6km (2000-6000m); 
## hotspots of infected vector at 0.2-1.8km (200-1800m); 
## flight range of triatomine varies from 0.1-1.5km (100-1500)
## given the above, let's set max.edge for inner edge at around 500m (min. average flight), the outer edge to 1500m, and the cutoff to 10m (if we imagine triatomine walking)

# plot mesh and coordinate points
png("figures/inla/mesh/mesh_infestation.png",width=6,height=8,units="in",res=600)
plot(mesh, main = "")
points(Locations, col = "red", cex = 0.7)
dev.off()

# select mesh
Mesh <- mesh

# Spatial/SPDE Construction ----------------------------------------------------

# make SPDE for spatial structure # this function will be equivalent to using a GMRF in a spatial model
spde <- inla.spde2.pcmatern(Mesh, 
                            alpha=2, #alpha = nu + delta/2 = 1 + 2/2 = 2
                            # smoothness parameter controls how quickly the spatial effect changes over distance (e.g., neighboring locations are likely to share similar random effects)
                            prior.range = c(1500, 0.95), #SPDE model is operating in meters; so priors should also be in meters
                            # prior.range = c(0.5, 0.5), #assume distribution is skewed to the left (e.g., average flight radius is 500m)
                            ## specifies spatial range, which defines how far apart two points can be while still maintaining a meaningful correlation.
                            ## (range0, Prange) whereby P(range < range0) = Prange
                            ## let's set this to the maximum flight range of a triatomine, 1.5km, where we set an a priori probability of 95% that the range will be greater than 1.5km.
                            ## source: T. infestans may easily fly >550 m (SchoÞeld et al. 1992) or reach 1,500 m (Schweigmann et al. 1988) in an open field
                            prior.sigma = c(.5, .5)
                            #P(sigma > 1.0) = 0.5
                            ## specifies standard deviation (overall variability) of the field: sigma0 and Psigma should reflect how much variation you're willing to allow 
                            ## if the data suggests that high variability in the spatial field is unlikely, you would want a smaller sigma0 and a lower Psigma
                            ## large sigma: If sigma is large, the spatial process is highly variable, meaning there's a lot of fluctuation or noise in the values at different locations. 
                            ## small sigma: If sigma is small, the spatial process is less variable, and neighboring locations are more similar to each other.
                            ## let's set sigma to a larger value b/c flight range of triatomine can vary widely (200-1500m) and spatial aggregation of cases at a higher spatial scale might occur b/c of dispersal of infected vectors
)

    # ## Sensitivity analysis: spatial range prior
    # 
    #       # weaker prior (less informative)
    #       spde_weak_rho <- inla.spde2.pcmatern(
    #         mesh = Mesh,
    #         alpha = 2,
    #         prior.range = c(1500, 0.5),
    #         prior.sigma = c(0.5, 0.5)
    #       )
    #       
    #       # shorter range
    #       spde_short_rho <- inla.spde2.pcmatern(
    #         mesh = Mesh,
    #         alpha = 2,
    #         prior.range = c(500, 0.95),
    #         prior.sigma = c(0.5, 0.5)
    #       )
    #       
    #       # longer range
    #       spde_long_rho <- inla.spde2.pcmatern(
    #         mesh = Mesh,
    #         alpha = 2,
    #         prior.range = c(3000, 0.95),
    #         prior.sigma = c(0.5, 0.5)
    #       )
    #       
    #   ## Sensitivity analysis: sigma (variability) prior
    #       
    #       # more restrictive (low variability)
    #       spde_low_var <- inla.spde2.pcmatern(
    #         mesh = Mesh,
    #         alpha = 2,
    #         prior.range = c(1500, 0.95), 
    #         prior.sigma = c(0.2, 0.5) # 50% chance sigma > 0.2 --> favors smoother field
    #       )
    #       
    #       # more flexible (high variability)
    #       spde_high_var <- inla.spde2.pcmatern(
    #         mesh = Mesh,
    #         alpha = 2,
    #         prior.range = c(1500, 0.95),
    #         prior.sigma = c(1.0, 0.5) # 50% chance sigma > 1.0 --> allows high variability
    #       )
    #       
    #       # more conservative (smaller scale)
    #       spde_cons_var <- inla.spde2.pcmatern(
    #         mesh = Mesh,
    #         alpha = 2,
    #         prior.range = c(1500, 0.95),
    #         prior.sigma = c(0.3, 0.9) # 90% chance sigma < 0.3 --> penalizes high variability
    #       )
    # 
    # ## If running, sensitivity analysis, choose set of hyperparameter priors from the above
    #       
    #       # spde <- spde_weak_rho
    #       # spde <- spde_short_rho
    #       # spde <- spde_long_rho
    #       # spde <- spde_low_var
    #       # spde <- spde_high_var
    #       # spde <- spde_cons_var
          
# create all the required indexes for the SPDE model and name the effect "spatial.field" (for use in creating our formula)
w.spde <- inla.spde.make.index('spatial.field', n.spde = spde$n.spde) # making the w (spatial field)
   

# INLA Stack Construction ------------------------------------------------------

# For Infestation Model
    
    ## get df
    df <- df_inf
    
    ## A matrix: translates spatial locations on the mesh into vectors in the model (maps the GMRF from mesh nodes to the n observation location)
    
        ## get sampling locations - house coordinates
        Locations = as.data.frame(cbind(df$x, df$y))
        
        ## assign UTM Zone 21S CRS (EPSG:32721)
        library(sf)
        Locations <- st_as_sf(Locations, coords = c("V1", "V2"),
                              crs = st_crs(32721)) # UTM Zone 21S: EPSG:32721
        
        ## make the A matrix
        A <- inla.spde.make.A(Mesh, loc = Locations) # Construct observation weight matrices for models based on inla.mesh() objects where the 1st argument is an inla.mesh() object, and the 2nd are the observation coordinates
    
    ## data stacking 
        
        # estimation stack: create an INLA stack with the input data as a list of the response, 
        stk_inf <- inla.stack(
          data = list(y = df[,resp_inf]), #specify resp1 as response variable
          A = list(
            1,         # fixed effects
            A     # spatial field projector
          ),
          effects = list(
            fixed = df[,covar_inf],
            spatial.field = w.spde$spatial.field # spde random effect
          ),
          tag = "est"
        )
     
        
# For Reinfestation Model
        
    ## data stacking
        
        # create an INLA stack with the input data as a list of the response, 
        stk_reinf <- inla.stack(
          data = list(y = df[,resp_reinf]), #specify resp2 as response variable
          A = list(
            1,         # fixed effects
            A     # spatial field projector
          ),
          effects = list(
            fixed = df[,covar_inf],
            spatial.field = w.spde$spatial.field # spde random effect
          ),
          tag = "est"
        )

# For IVA Model
        
    ## get df
    df <- df_iva
        
    ## A matrix: translates spatial locations on the mesh into vectors in the model (maps the GMRF from mesh nodes to the n observation location)
        
        # get sampling locations - house coordinates
        Locations = as.data.frame(cbind(df$x, df$y))
        
        # assign UTM Zone 21S CRS (EPSG:32721)
        library(sf)
        Locations <- st_as_sf(Locations, coords = c("V1", "V2"),
                              crs = st_crs(32721)) # UTM Zone 21S: EPSG:32721
        
        # make the A matrix
        A <- inla.spde.make.A(Mesh, loc = Locations) # Construct observation weight matrices for models based on inla.mesh() objects where the 1st argument is an inla.mesh() object, and the 2nd are the observation coordinates
        
    ## data stacking 
        
        # create an INLA stack with the input data as a list of the response, 
        stk_iva <- inla.stack(
          data = list(y = df[,resp_iva]), #specify resp3 as response variable
          A = list(
            1,         # fixed effects
            A     # spatial field projector
          ),
          effects = list(
            fixed = df[,covar_inf], #same covariates as infestation model
            spatial.field = w.spde$spatial.field # spde random effect
          ),
          tag = "est"
        )
        
# For Child Infection Model
        
    ## get df
    df <- df_child

    # A matrix: translates spatial locations on the mesh into vectors in the model (maps the GMRF from mesh nodes to the n observation location)
    
        ## get sampling locations - house coordinates
        Locations = as.data.frame(cbind(df$x, df$y))
        
        ## assign UTM Zone 21S CRS (EPSG:32721)
        library(sf)
        Locations <- st_as_sf(Locations, coords = c("V1", "V2"),
                              crs = st_crs(32721)) # UTM Zone 21S: EPSG:32721
    
        ## make the A matrix
        A <- inla.spde.make.A(Mesh, loc = Locations) # Construct observation weight matrices for models based on inla.mesh() objects where the 1st argument is an inla.mesh() object, and the 2nd are the observation coordinates
    
    # data stacking 
        
        ## create an INLA stack with the input data as a list of the response, 
        stk_child <- inla.stack(
          data = list(y = df[,resp_child]),
          A = list(
            1,    # fixed effects
            1,    # household random effect
            A     # spatial field projector
          ),
          effects = list(
            fixed = df[,covar_child],    # fixed effects
            hhx.id = df$hhx.id,    # household random effect
            spatial.field = w.spde$spatial.field # spde random effect
          ),
          tag = "est"
        )
        
################################################################################
### INLA: INFESTATION
################################################################################

# Get INLA stack
stk <- stk_inf
        
# Explore Models ---------------------------------------------

# baseline
formula <- y ~ svi + hai + mobility # waic = 366.33

# explore network-constructed vars
formula <- y ~ svi + hai + mobility + kin.svi.med # ns: kin.svi.med+; waic = 364.60
formula <- y ~ svi + hai + mobility + kin.svi.relmed # ns: kin.svi.relmed+, waic = 364.60
formula <- y ~ svi + hai + mobility + kin.occ.inf # ns: kin.occ.inf+, waic = 365.89
formula <- y ~ svi + hai + mobility + kin.iva.med # ns: kin.iva.med+, waic = 368.77

# explore network centrality metrics
formula <- y ~ svi + hai + mobility + degree # ns: mobility-, degree-; waic = 367.59
formula <- y ~ svi + hai + mobility + connected # ns: mobility-, connected-; waic = 367.30
formula <- y ~ svi + hai + mobility + harmonic # ns: mobility-, harmonic-; waic = 268.33
  
# explore interaction w/ mobility
formula <- y ~ svi + hai + mobility + degree + mobility*degree # ns: mobility-, degree-, mobility*degree+; waic = 368.51
formula <- y ~ svi + hai + mobility + connected + mobility*connected # sig: mobility-, connected-, mobility*connected+; waic = 364.02
formula <- y ~ svi + hai + mobility + harmonic + mobility*harmonic # ns: mobility-, harmonic-, mobility*harmonic+; waic = 368.51
    
    # explore spatial random effect
    formula <- y ~ svi + hai + mobility + connected + mobility*connected + f(spatial.field, model=spde) # sig: mobility-, mobility*connected+; ns: connected-; waic = 329.88
    
## fit the model (the INLA function)
model <- inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
              control.family = list(link = "logit"),
              control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  
              control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
              verbose = FALSE)
summary(model)  

# explore alternative binary coding of mobility
formula <- y ~ svi + hai + mobility + connected + mobility*connected + f(spatial.field, model=spde)
formula <- y ~ svi + hai + mobility.binrev + connected + mobility.binrev*connected + f(spatial.field, model=spde)

# Final Models of Infestation --------------------------------------------------

# load libraries
library(INLA)
library(bayesplot)
library(ggplot2)
library(cowplot)
library(lattice)

# list all formulas
formulas <- list(
  infestation0 = y ~ f(spatial.field, model=spde),
  infestation1.0 = y ~ svi + hai + mobility + connected + mobility*connected,
  infestation1.1 = y ~ svi + hai + mobility + connected + mobility*connected + f(spatial.field, model=spde), #BEST FIT#
  # infestation1.1a = y ~ svi + hai + mobility.bin + connected.bin + mobility.bin*connected.bin + f(spatial.field, model=spde), #BEST FIT# sig: svi+, hai+, mobility-, connected*mobility+; ns: connected-; waic=329.20
  # infestation1.1b = y ~ svi + hai + mobility.bin + connected.binrev + mobility.bin*connected.binrev + f(spatial.field, model=spde),
  # infestation1.1c = y ~ svi + hai + mobility.binrev + connected.bin + mobility.binrev*connected.bin + f(spatial.field, model=spde),
  # infestation1.1d = y ~ svi + hai + mobility.binrev + connected.binrev + mobility.binrev*connected.binrev + f(spatial.field, model=spde),
  infestation1.2 = y ~ svi + hai + mobility + connected + f(spatial.field, model=spde) # no interaction term
)

# prepare list to store results
results <- list()

# loop to fit and plot
for (i in seq_along(formulas)) {
  
    name <- names(formulas)[i]
    formula <- formulas[[i]]
    
    # fit model
    model <-inla(formula, 
                 data=inla.stack.data(stk, spde=spde), #extracts dataset from INLA stack object and tells INLA to use a stoch. partial differential equation object to represent the spatial random field
                 family= 'binomial', #specifies likelihood for the outcome variable
                 control.family = list(link = "logit"), #configures details of likelihood family: link="logit" means we are fitting a logistic regression
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE), #controls how the predictor is built: A=...supplies the projection matrix that maps latent effects (only applied if adjusting for spatial effects)
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
  
    # fixed effects posterior plot
        
        # get fixed effects names (minus intercept)
        names_fe <- names(model$marginals.fixed)[-1]
        
          # # transform marginals to odds scale (exp of log-odds)
          # marginals_odds <- lapply(model$marginals.fixed, function(marg) {
          #   inla.tmarginal(function(x) exp(x), marg)
          # })
        
        # sample from marginals (keep as log-odds or use odds scale from marginals_odds above)
        set.seed(4321)
        posterior_samples <- sapply(model$marginals.fixed, function(marg) {
        # posterior_samples <- sapply(marginals_odds, function(marg) {
          inla.rmarginal(10000, marg)
        })
        
        # plot
        g_posterior <- mcmc_areas(posterior_samples,
                                  pars = names_fe,
                                  prob = 0.95) +
          xlab("Posterior Log-Odds") +
          # xlab("Posterior Odds") +
          ggtitle(paste0("WAIC = ", round(model$waic$waic,2))) +
          theme_minimal(base_family = "Arial") +
          theme(plot.title = element_text(size = 16, face = "bold", family = "sans"), 
                axis.title = element_text(size = 16), 
                axis.text = element_text(size = 14))
        
        # save
        fe_count <- length(names(model$marginals.fixed)[-1])
        ggsave(
          filename = paste0("figures/inla/infestation/model_", name, "_fe.png"),
          # filename = paste0("figures/inla/infestation/model_", name, "_fe_odds.png"),
          plot = g_posterior, width = 6, height = fe_count, dpi = 600
        )
        
        # collect model fit metrics
        metrics <- data.frame(
          regression = paste(deparse(formula), collapse = ""),  # store formula as a single string
          WAIC = model$waic$waic,
          pWAIC = model$waic$p.eff,
          DIC = model$dic$dic,
          pDIC = model$dic$p.eff,
          converged = ifelse(model$mode$mode.status == 0, TRUE, FALSE),
          n_params = nrow(model$summary.fixed)
        )
        
        # save results into list
        results[[name]] <- list(
          model = model,
          plot = g_posterior,
          metrics = metrics
        )
    
    # spatial effects plots (if present)
    if ("spatial.field" %in% names(model$summary.random)) {      

        # posterior of range
      
            # get marginals of range
            spde_result <- inla.spde2.result(model, "spatial.field", spde, do.transf = TRUE)
            post_range <- spde_result$marginals.range$range
            df_range <- data.frame(post_range)
            
            # plot
            g_range <- ggplot(df_range, aes(x = x, y = y)) +
              geom_line() +
              geom_vline(xintercept = 1500, color = "red", linetype = "dashed") +
              geom_vline(xintercept = inla.qmarginal(0.5, post_range), color = "blue", linetype = "dashed") +
              labs(x = "Spatial range (meters)",
                   y = "Posterior density",
                   title = "Posterior distribution of spatial range") +
              theme_minimal(base_size = 14) +
              theme(axis.line = element_blank())
            
            # save
            ggsave(
              filename = paste0("figures/inla/infestation/model_", name, "_sfrange.png"),
              plot = g_range, width = 4, height = 4, dpi = 600
            )
        
        # project spatial field mean/sd

            # project mesh, mean, and sd
            gproj <- inla.mesh.projector(mesh, dims = c(300, 300))
            g.mean <- inla.mesh.project(gproj, model$summary.random$spatial.field$mean)
            g.sd   <- inla.mesh.project(gproj, model$summary.random$spatial.field$sd)
            
            # diverging palette for mean
            val_min <- min(model$summary.random$spatial.field$mean, na.rm = TRUE)
            val_max <- max(model$summary.random$spatial.field$mean, na.rm = TRUE)
            prop_neg <- abs(val_max) / (abs(val_min) + val_max)
            n = 200
            n_neg <- max(1, round(n * prop_neg))
            n_pos <- max(1, n - n_neg)
            at <- c(seq(val_min, 0, length.out = n_neg + 1),
                    seq(0, val_max, length.out = n_pos + 1)[-1])
            cols <- c(colorRampPalette(c("blue", "white"))(n_neg),
                      colorRampPalette(c("white", "red"))(n_pos))
            div_mean <- list(at = at, cols = cols)
            
            # red scale for sd
            at_sd <- seq(min(model$summary.random$spatial.field$sd, na.rm = TRUE),
                         max(model$summary.random$spatial.field$sd, na.rm = TRUE), length.out = 200)
            cols_sd <- colorRampPalette(c("red", "white"))(200)
            div_sd <- list(at = at_sd, cols = cols_sd)
            
            my_par <- list(par.xlab.text = list(cex = 1.2, font = 2),
                           par.ylab.text = list(cex = 1.2, font = 2))
            
            # plot
            g_sf <- plot_grid(
              levelplot(g.mean, scales = list(draw = FALSE), xlab = "", ylab = "", main = "mean",
                        region = TRUE, at = div_mean$at, col.regions = div_mean$cols,
                        par.settings = my_par),
              levelplot(g.sd, scales = list(draw = FALSE), xlab = "", ylab = "", main = "sd",
                        region = TRUE, at = div_sd$at, col.regions = div_sd$cols,
                        par.settings = my_par),
              nrow = 1
            )
            
            # save
            ggsave(
              filename = paste0("figures/inla/infestation/model_", name, "_sfmap.png"),
              plot = g_sf, width = 8, height = 4, dpi = 600
            )
    }
}        

# export spatial field map as raster for processing in QGIS
    
    ## formula
    formula = y ~ svi + hai + mobility + connected + mobility*connected + f(spatial.field, model=spde) #BEST FIT#
    
    ## fit model
    model <-inla(formula, 
                 data=inla.stack.data(stk, spde=spde), #extracts dataset from INLA stack object and tells INLA to use a stoch. partial differential equation object to represent the spatial random field
                 family= 'binomial', #specifies likelihood for the outcome variable
                 control.family = list(link = "logit"), #configures details of likelihood family: link="logit" means we are fitting a logistic regression
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE), #controls how the predictor is built: A=...supplies the projection matrix that maps latent effects (only applied if adjusting for spatial effects)
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
    
    ## project spatial field mean/sd on regular grid
    gproj <- inla.mesh.projector(mesh, dims = c(300, 300))
    g.mean <- inla.mesh.project(gproj, model$summary.random$spatial.field$mean)
    g.sd   <- inla.mesh.project(gproj, model$summary.random$spatial.field$sd)

    ## create rasters from projected mean
    # note: terra expects rows from top to bottom, so we flip vertically
    r_mean <- rast(
      nrows = length(gproj$y),
      ncols = length(gproj$x),
      xmin = min(gproj$x),
      xmax = max(gproj$x),
      ymin = min(gproj$y),
      ymax = max(gproj$y),
      vals = as.vector(g.mean[nrow(g.mean):1, ])
    )
    
    r_sd <- rast(
      nrows = length(gproj$y),
      ncols = length(gproj$x),
      xmin = min(gproj$x),
      xmax = max(gproj$x),
      ymin = min(gproj$y),
      ymax = max(gproj$y),
      vals = as.vector(g.sd[nrow(g.sd):1, ])
    )
    
    ## assign CRS
    crs(r_mean) <- "EPSG:32721"
    crs(r_sd) <- "EPSG:32721"
    
    ## export as GeoTIFF
    writeRaster(r_mean, "figures/maps/forQGIS/spatial_field_mean.tif", overwrite = TRUE)
    writeRaster(r_sd, "figures/maps/forQGIS/spatial_field_sd.tif", overwrite = TRUE)

# export house locations used in the model as a point layer to add to QGIS map
    
    ## load library
    library(sf)
    library(terra)
    
    ## print formula and response variable
    print(formula)
    print(resp_inf)
    names(df_inf)
    
    ## identify houses actually used in the infestation model based on variables above
    houses_used <- df_inf[complete.cases(
      df_inf[ ,c("v.inf", "svi", "hai", "mobility", "connected", "x", "y" )]
      ), ]
    
        ## check
        stk_dat <- inla.stack.data(stk, spde = spde)
        sum(!is.na(stk_dat$y))
        nrow(houses_used)
    
    ## turn coordinates to sf points (point object)
    houses_sf <- st_as_sf(
      houses_used,
      coords = c("x", "y"),  
      crs = 32721                  # same CRS as your raster
    )
    
    ## export as geopackage (for QGIS)
    st_write(houses_sf,
      "figures/maps/forQGIS/houses_used_infestation_model.gpkg",
      delete_dsn = TRUE
    )
    
# convert list of summary statistics for posteriors to table
    
    ## fixed effects
    fixed_table <- do.call(rbind, lapply(results, function(x) x$model$summary.fixed))
    fixed_table$kld <- NULL
    fixed_table$type <- "fixedeffects"
    fixed_table$order <- seq.int(nrow(fixed_table))
    
    ## hyperparameters
    hyperpar_table <- do.call(rbind, lapply(results, function(x) x$model$summary.hyperpar))
    hyperpar_table$type <- "hyperparameters"
    hyperpar_table$order <- seq.int(nrow(hyperpar_table))
    
    ## combine
    posteriors_table <- rbind(fixed_table, hyperpar_table)
    
        # save rownames as 'model' and 'parameter'
        library(tibble)
        library(stringr)
        posteriors_table <- posteriors_table %>%
            rownames_to_column("rowname") %>% 
            mutate(
              model     = str_extract(rowname, "^[^.]+\\.[^.]+"),   # everything up to 2nd dot
              parameter = str_remove(rowname, "^[^.]+\\.[^.]+\\.") # everything after 2nd dot
            ) %>%
            select(model, parameter, everything(), -rowname)

    ## convert list of model performance to table
    waic_table <- do.call(rbind, lapply(results, function(x) x$metrics))
    
        # save rowname (model name) as colname
        waic_table$model <- rownames(waic_table)
        rownames(waic_table) <- NULL 
        
    ## merge
    summary_table <- merge(posteriors_table, waic_table, by = "model")

    ## rename '0.5quant' to 'median'
    colnames(summary_table)[colnames(summary_table) == "0.5quant"] <- "median"
    
    ## sort rows
    summary_table <- summary_table[order(summary_table$model, summary_table$type, summary_table$order), ]
    
    ## reorder and remove some vars
    summary_table <- subset(summary_table, select = c("model", "regression", "parameter", "mean", "0.025quant", "median", "0.975quant", "mode", "sd", "WAIC"))
    
# save tables
table_summary_infestation <- summary_table    
# write.csv(summary_table, "tables/inla/table_summary_infestation.csv") 


# Post-hoc calculations --------------------------------------------------------

# Get 95% CrI for the joint posterior distribution of beta_mover_connected (the average effect of connected when mobility = 1)

    ## Why do we need the joint posterior?
        # Interactions are combinations of coefficients.
        # To calculate their posterior distributions, you need the joint posterior (i.e., correlation b/w coefficients).
        # We can't use marginals.fixed() b/c this gives us the marginal posterior distribution for each coefficient, which assumes independence between them.
        # So we use inla.posterior.sample(), which gives us the posterior samples (each draw represents a joint realization of all parameters simultaneously).
        # So if we want the joint posterior of beta_mover_connected (= beta_mobility + beta_connected + beta_mobility*connected),...
        # We can compute beta_mover_connected for each posterior draw.
        # This preserves the joint covariance structure!

    ## get results of final model
    formula <- y ~ svi + hai + mobility + connected + mobility*connected + f(spatial.field, model=spde) #BEST FIT#
    model <-inla(formula, 
                 data=inla.stack.data(stk, spde=spde), #extracts dataset from INLA stack object and tells INLA to use a stoch. partial differential equation object to represent the spatial random field
                 family= 'binomial', #specifies likelihood for the outcome variable
                 control.family = list(link = "logit"), #configures details of likelihood family: link="logit" means we are fitting a logistic regression
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE), #controls how the predictor is built: A=...supplies the projection matrix that maps latent effects (only applied if adjusting for spatial effects)
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))

    ## extract posterior samples for the relevant coefficients: beta_connected and beta_mobility*connected
    set.seed(4321)
    samples <- inla.posterior.sample(n = 2000, model)   # draw 2k posterior samples
    
    ## function that computes beta_mover_connected
    fun_mover_connected <- function() {
      mobility + connected + `mobility:connected`
    }
    
        ## evaluate the function over the samples
        eval_mc <- inla.posterior.sample.eval(fun_mover_connected, samples)
        
        ## summarize
        quants <- quantile(eval_mc, probs = c(0.025, 0.50, 0.975))
        print(quants)
        
        ## summarize in odds scale
        or_draws <- exp(eval_mc)
        or_quants <- quantile(or_draws, probs = c(0.025, 0.50, 0.975))
        print(or_quants)
        
    ## function that computes beta_connected + beta_mobility*connected
    fun_effect_among_movers <- function() {
      connected + `mobility:connected`
    }
        
        ## evaluate the function over the samples
        eval_mc <- inla.posterior.sample.eval(fun_effect_among_movers, samples)
        
        ## summarize
        quants <- quantile(eval_mc, probs = c(0.025, 0.50, 0.975))
        print(quants)
        
        ## summarize in odds scale
        or_draws <- exp(eval_mc)
        or_quants <- quantile(or_draws, probs = c(0.025, 0.50, 0.975))
        print(or_quants)
    
 
################################################################################
### INLA: REINFESTATION 
################################################################################
        
# Get INLA stack
stk <- stk_reinf

# Explore Models -------------------------------------------    

# formula: 
            
    ## centrality metric only
    formula <- y ~ degree # sig: degree-1.437; waic=61.28
    formula <- y ~ connected # sig: connected-2.604; waic=58.20
    formula <- y ~ harmonic # ns: harmonic-; waic=62.23
    formula <- y ~ betweenness # convergence issue
    formula <- y ~ kin.occ.inf # ns: kin.occ.inf-; waic=2275.64
    formula <- y ~ f(degree, model = "rw1") # model as nonlinear w/ rw1: waic=62.20
    
    ## additional covariates
    formula <- y ~ svi + hai + mobility + connected # sig: connected-2.698; ns: svi-, hai+, mobility+; waic=63.12
    formula <- y ~ mobility + connected # sig: -connected; ns: mobility+; waic = 60.09
    
    ## interaction
    formula <- y ~ mobility + connected + mobility*connected #CONVERGENCE WARNING
    formula <- y ~ mobility + degree + mobility*degree #CONVERGENCE WARNING
    
    ## spatial field
    formula <- y ~ mobility + connected + f(spatial.field, model=spde) #CONVERGENCE WARNING
    formula <- y ~ connected + f(spatial.field, model=spde) #CONVERGENCE WARNING (waic not improved)
      
## fit the model (the INLA function)
model <-inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
             control.family = list(link = "logit"),
             control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  
             control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
             verbose = FALSE) 
summary(model)

# Final Models of Reinfestation ------------------------------------------------    
    
# load
library(INLA)
library(bayesplot)
library(ggplot2)

# list all formulas
formulas <- list(
  reinfestation1 = y ~ mobility + connected 
)

# prepare list to store results
results <- list()

# loop to fit and plot
for (i in seq_along(formulas)) {
      
    name <- names(formulas)[i]
    formula <- formulas[[i]]
    
    # fit model
    model <-inla(formula, 
                 data=inla.stack.data(stk, spde=spde), #extracts dataset from INLA stack object and tells INLA to use a stoch. partial differential equation object to represent the spatial random field
                 family= 'binomial', #specifies likelihood for the outcome variable
                 control.family = list(link = "logit"), #configures details of likelihood family: link="logit" means we are fitting a logistic regression
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE), #controls how the predictor is built: A=...supplies the projection matrix that maps latent effects (only applied if adjusting for spatial effects)
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
    
    # fixed effects posterior plot
    
        # get fixed effects names (minus intercept)
        names_fe <- names(model$marginals.fixed)[-1]
        
        # get samples for each fixed effect
        set.seed(4321)
        posterior_samples <- sapply(model$marginals.fixed, function(marg) {
          inla.rmarginal(10000, marg)
        })
        
        # plot
        g_posterior <- mcmc_areas(posterior_samples,
                                  pars = names_fe,
                                  prob = 0.95) +
          ggtitle(paste0("WAIC = ", round(model$waic$waic,2))) +
          theme_minimal(base_family = "Arial") +
          theme(plot.title = element_text(size = 16, face = "bold", family = "sans"), 
                axis.title = element_text(size = 14), 
                axis.text = element_text(size = 12))
        
        # save
        fe_count <- length(names(formulas)[i])
        # fe_count <- length(names(model$marginals.fixed)[-1])
        
        ggsave(
          filename = paste0("figures/inla/reinfestation/model_", name, "_fe.png"),
          plot = g_posterior, width = 6, height = fe_count + 2, dpi = 600
        )
        
    # collect model fit metrics
    metrics <- data.frame(
      regression = paste(deparse(formula), collapse = ""),  # store formula as a single string
      WAIC = model$waic$waic,
      pWAIC = model$waic$p.eff,
      DIC = model$dic$dic,
      pDIC = model$dic$p.eff,
      converged = ifelse(model$mode$mode.status == 0, TRUE, FALSE),
      n_params = nrow(model$summary.fixed)
    )
    
    # save results into list
    results[[name]] <- list(
      model = model,
      plot = g_posterior,
      metrics = metrics
    )

}        
    
# convert list of summary statistics for posteriors to table
    
    ## fixed effects
    fixed_table <- do.call(rbind, lapply(results, function(x) x$model$summary.fixed))
    fixed_table$kld <- NULL

        # save rownames as 'model' and 'parameter'
        library(tibble)
        library(stringr)
        fixed_table <- fixed_table %>%
          rownames_to_column("rowname") %>% 
          mutate(
            model     = str_extract(rowname, "^[^.]+"),   # everything up to 1st dot
            parameter = str_remove(rowname, "^[^.]+\\.") # everything after 1st dot
          ) %>%
          select(model, parameter, everything(), -rowname)
    
    ## convert list of model performance to table
    waic_table <- do.call(rbind, lapply(results, function(x) x$metrics))
        
        # save rowname (model name) as colname
        waic_table$model <- rownames(waic_table)
        rownames(waic_table) <- NULL 
    
    ## merge
    summary_table <- merge(fixed_table, waic_table, by = "model")
    
    ## rename '0.5quant' to 'median'
    colnames(summary_table)[colnames(summary_table) == "0.5quant"] <- "median"
    
    ## reorder and remove some vars
    summary_table <- subset(summary_table, select = c("model", "regression", "parameter", "mean", "0.025quant", "median", "0.975quant", "mode", "sd", "WAIC"))

# save tables
table_summary_reinfestation <- summary_table
# write.csv(summary_table, "tables/inla/table_summary_reinfestation.csv") 



################################################################################
### INLA: INFECTED VECTOR ABUNDANCE 
################################################################################

# Get INLA stack
stk <- stk_iva

# # Select likelihood family -------------------------------------------------------------------------
#   
#     ## get options
#     names(inla.models()$likelihood)
#     
#     ## loop through different likelihood options
#         
#         # define formula
#         formula <- y ~ f(spatial.field, model = spde)
#         
#         # create list of likelihood families
#         families <- c("poisson", "nbinomial", "zeroinflatedpoisson1", "zeroinflatednbinomial1", "zeroinflatednbinomial2")
#         
#         # create list to store model 
#         models <- list()
#       
#         # run model in loop
#         for (i in seq_along(families)) {
#           
#           fam <- families[[i]]
#     
#           models[[fam]] <- inla(formula,
#                          family = fam, 
#                          data = inla.stack.data(stk, spde = spde), 
#                          control.predictor=list(A=inla.stack.A(stk),compute=TRUE),
#                          control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))      
#         }
#     
#     # compare WAIC values
#     models$poisson$waic$waic #5th
#     models$nbinomial$waic$waic #4th
#     models$zeroinflatedpoisson1$waic$waic # Winner winner, chicken dinner!
#     models$zeroinflatednbinomial1$waic$waic #2nd
#     models$zeroinflatednbinomial2$waic$waic #3rd

# Explore Models -------------------------------

# formula
    
    ## baseline
    formula <- y ~ degree # sig: degree-; waic=339.01
    formula <- y ~ connected # ns: connected-; waic=383.70
    formula <- y ~ harmonic # ns: harmonic-; waic=469.33
    formula <- y ~ kin.occ.inf # sig: kin.occ.inf-; waic=312.70
    formula <- y ~ kin.iva.med # sig: kin.iva.med-; waic=313.10
    
    ## combine
    formula <- y ~ degree + kin.occ.inf # sig: kin.occ.inf-; ns: degree+; waic = 337.15
    formula <- y ~ degree + kin.iva.med # ns: degree-, kin.iva.med-; waic = 336.61
    formula <- y ~ kin.occ.inf + kin.iva.med # sig: kin.occ.inf-; ns: kin.iva.med; waic = 312.70
    formula <- y ~ degree + kin.occ.inf + kin.iva.med # ns: degree+, kin.occ.inf-, kin.iva.med-; waic = 337.18
    
    ## interaction?
    formula <- y ~ degree + kin.occ.inf + degree*kin.occ.inf # ns: all; waic = 341.19
    formula <- y ~ mobility + kin.occ.inf + mobility*kin.occ.inf # sig: kin.occ.inf-; ns: mobility+, mobility*kin.occ.inf; waic = 339.19
    formula <- y ~ mobility + degree + mobility*degree # ns: all; waic = 404.95
    formula <- y ~ mobility + connected + mobility*connected # sig: mobility+, mobility*connected-; connected-; waic = 600.86
    
    ## explore additional covariates
    formula <- y ~ svi + hai + degree  # sig: svi+, hai+, degree-; waic=310.77
    formula <- y ~ svi + hai + kin.occ.inf # sig: kin.occ.inf-; waic=305.65
    formula <- y ~ svi + hai + kin.iva.med # sig: kin.iva.med-; waic=307.83
    
    ## explore spatial field
    formula <- y ~ degree + f(spatial.field, model=spde)  # ns: degree-; waic = 233.19
    formula <- y ~ kin.occ.inf + f(spatial.field, model=spde) # ns: kin.occ.inf-; waic=222.22
    formula <- y ~ kin.iva.med + f(spatial.field, model=spde) # ns: kin.iva.med-; waic=234.74
    formula <- y ~ svi + hai + degree + f(spatial.field, model=spde)  # sig: svi+, hai+; ns: degree-; waic=206.69
    
## fit the model (the INLA function)
model <-inla(formula, data=inla.stack.data(stk,spde=spde),
             family = "zeroinflatedpoisson1",
             control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  
             control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
             verbose = FALSE) 
summary(model)

# Final Models of IVA ------------------------------------------------    

# load
library(INLA)
library(bayesplot)
library(ggplot2)

# list all formulas
formulas <- list(
  iva0 = y ~ f(spatial.field, model=spde),
  iva1 = y ~ degree,
  iva2 = y ~ kin.occ.inf,
  iva3 = y ~ kin.iva.med
)

# prepare list to store results
results <- list()

# loop to fit and plot
for (i in seq_along(formulas)) {
    
    name <- names(formulas)[i]
    formula <- formulas[[i]]
    
    # fit model
    model <-inla(formula, 
                 data=inla.stack.data(stk,spde=spde),
                 family = "zeroinflatedpoisson1",
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
    
    # fixed effects posterior plot
    
        # get fixed effects names (minus intercept)
        names_fe <- names(model$marginals.fixed)[-1]
        
        # # transform marginals to odds scale (exp of log-odds)
        # marginals_odds <- lapply(model$marginals.fixed, function(marg) {
        #   inla.tmarginal(function(x) exp(x), marg)
        # })
        
        # sample from marginals (keep as log-odds or use odds scale from marginals_odds above)
        set.seed(4321)
        posterior_samples <- sapply(model$marginals.fixed, function(marg) {
          # posterior_samples <- sapply(marginals_odds, function(marg) {
          inla.rmarginal(10000, marg)
        })
        
        # plot
        g_posterior <- mcmc_areas(posterior_samples,
                                  pars = names_fe,
                                  prob = 0.95) +
          xlab("Posterior Log-Odds") +
          # xlab("Posterior Odds") +
          ggtitle(paste0("WAIC = ", round(model$waic$waic,2))) +
          theme_minimal(base_family = "Arial") +
          theme(plot.title = element_text(size = 16, face = "bold", family = "sans"), 
                axis.title = element_text(size = 16), 
                axis.text = element_text(size = 14))
        
        # save
        fe_count <- length(names(model$marginals.fixed)[-1])
        ggsave(
          filename = paste0("figures/inla/iva/model_", name, "_fe.png"),
          # filename = paste0("figures/inla/iva/model_", name, "_fe_odds.png"),
          plot = g_posterior, width = 6, height = fe_count, dpi = 600
        )
          
    # collect model fit metrics
    metrics <- data.frame(
      regression = paste(deparse(formula), collapse = ""),  # store formula as a single string
      WAIC = model$waic$waic,
      pWAIC = model$waic$p.eff,
      DIC = model$dic$dic,
      pDIC = model$dic$p.eff,
      converged = ifelse(model$mode$mode.status == 0, TRUE, FALSE),
      n_params = nrow(model$summary.fixed)
    )
    
    # save results into list
    results[[name]] <- list(
      model = model,
      plot = g_posterior,
      metrics = metrics
    )
    
    # spatial effects plots (if present)
    
    if ("spatial.field" %in% names(model$summary.random)) {      
      
        # posterior of range
        
        # get marginals of range
        spde_result <- inla.spde2.result(model, "spatial.field", spde, do.transf = TRUE)
        post_range <- spde_result$marginals.range$range
        df_range <- data.frame(post_range)
        
        # plot
        g_range <- ggplot(df_range, aes(x = x, y = y)) +
          geom_line() +
          geom_vline(xintercept = 1500, color = "red", linetype = "dashed") +
          geom_vline(xintercept = inla.qmarginal(0.5, post_range), color = "blue", linetype = "dashed") +
          labs(x = "Spatial range (meters)",
               y = "Posterior density",
               title = "Posterior distribution of spatial range") +
          theme_minimal(base_size = 14) +
          theme(axis.line = element_blank())
        
        # save
        ggsave(
          filename = paste0("figures/inla/iva/model_", name, "_sfrange.png"),
          plot = g_range, width = 4, height = 4, dpi = 600
        )
      
    # project spatial field mean/sd
    
        # project mesh, mean, and sd
        gproj <- inla.mesh.projector(mesh, dims = c(300, 300))
        g.mean <- inla.mesh.project(gproj, model$summary.random$spatial.field$mean)
        g.sd   <- inla.mesh.project(gproj, model$summary.random$spatial.field$sd)
        
        # diverging palette for mean
        val_min <- min(model$summary.random$spatial.field$mean, na.rm = TRUE)
        val_max <- max(model$summary.random$spatial.field$mean, na.rm = TRUE)
        prop_neg <- abs(val_max) / (abs(val_min) + val_max)
        n = 200
        n_neg <- max(1, round(n * prop_neg))
        n_pos <- max(1, n - n_neg)
        at <- c(seq(val_min, 0, length.out = n_neg + 1),
                seq(0, val_max, length.out = n_pos + 1)[-1])
        cols <- c(colorRampPalette(c("blue", "white"))(n_neg),
                  colorRampPalette(c("white", "red"))(n_pos))
        div_mean <- list(at = at, cols = cols)
        
        # red scale for sd
        at_sd <- seq(min(model$summary.random$spatial.field$sd, na.rm = TRUE),
                     max(model$summary.random$spatial.field$sd, na.rm = TRUE), length.out = 200)
        cols_sd <- colorRampPalette(c("red", "white"))(200)
        div_sd <- list(at = at_sd, cols = cols_sd)
        
        my_par <- list(par.xlab.text = list(cex = 1.2, font = 2),
                       par.ylab.text = list(cex = 1.2, font = 2))
        
        # plot
        g_sf <- plot_grid(
          levelplot(g.mean, scales = list(draw = FALSE), xlab = "", ylab = "", main = "mean",
                    region = TRUE, at = div_mean$at, col.regions = div_mean$cols,
                    par.settings = my_par),
          levelplot(g.sd, scales = list(draw = FALSE), xlab = "", ylab = "", main = "sd",
                    region = TRUE, at = div_sd$at, col.regions = div_sd$cols,
                    par.settings = my_par),
          nrow = 1
        )
        
        # save
        ggsave(
          filename = paste0("figures/inla/iva/model_", name, "_sfmap.png"),
          plot = g_sf, width = 8, height = 4, dpi = 600
        )
  }
}        

# convert list of summary statistics for posteriors to table
    
    ## fixed effects
    fixed_table <- do.call(rbind, lapply(results, function(x) x$model$summary.fixed))
    fixed_table$kld <- NULL
    fixed_table$type <- "fixedeffects"
    fixed_table$order <- seq.int(nrow(fixed_table))
    
    ## hyperparameters
    hyperpar_table <- do.call(rbind, lapply(results, function(x) x$model$summary.hyperpar))
    hyperpar_table$type <- "hyperparameters"
    hyperpar_table$order <- seq.int(nrow(hyperpar_table))
    
    ## combine
    posteriors_table <- rbind(fixed_table, hyperpar_table)
        
        # save rownames as 'model' and 'parameter'
        library(tibble)
        library(stringr)
        posteriors_table <- posteriors_table %>%
          rownames_to_column("rowname") %>% 
          mutate(
            model     = str_extract(rowname, "^[^.]+"),   # everything up to 1st dot
            parameter = str_remove(rowname, "^[^.]+\\.") # everything after 1st dot
          ) %>%
          select(model, parameter, everything(), -rowname)
    
    ## convert list of model performance to table
    waic_table <- do.call(rbind, lapply(results, function(x) x$metrics))
        
        # save rowname (model name) as colname
        waic_table$model <- rownames(waic_table)
        rownames(waic_table) <- NULL 
    
    ## merge
    summary_table <- merge(posteriors_table, waic_table, by = "model")
    
    ## rename '0.5quant' to 'median'
    colnames(summary_table)[colnames(summary_table) == "0.5quant"] <- "median"
    
    ## sort rows
    summary_table <- summary_table[order(summary_table$model, summary_table$type, summary_table$order), ]
    
    ## reorder and remove some vars
    summary_table <- subset(summary_table, select = c("model", "regression", "parameter", "mean","0.025quant", "median",  "0.975quant", "mode", "sd", "WAIC"))
    
    # save tables
    table_summary_iva <- summary_table
    # write.csv(summary_table, "tables/inla/table_summary_iva.csv") 
    


################################################################################
### INLA: INFECTION/CHILD SEROPOSITIVITY
################################################################################

# Get A matrix and INLA stack
stk <- stk_child

# Model probability of child infection -----------------------------------------

# baseline
formula <- y ~ age + mom.serology + svi + hai + iva + mover
formula <- y ~ age + mom.serology + svi + hai + iva + mover
formula <- y ~ age + mom.serology + svi + hai + iva + mover
formula <- y ~ age + mom.serology + svi + hai + iva + move

# explore network-constructed vars
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.svi.med # ns: kin.svi.relmed+; waic = 280.23
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.svi.relmed # ns: kin.svi.relmed+; waic = 280.23
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf # sig: kin.occ.inf+; waic = 273.25 !!! WINNER !!!
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.iva.med # ns: kin.iva.med+; waic = 279.50

# explore network centrality metrics
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + degree # ns: degree+; waic = 275.29
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + connected # ns: connected-; waic = 271.89
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + harmonic # ns: harmonic-; waic = 274.31 !!! WINNER !!!

# explore interaction w/ mover
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + degree + mover*degree # sig: mover-, mover*degree+; ns: kin.occ.inf+, degree-; waic = 269.20
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + connected + mover*connected # sig: connected-; ns: mover-, mover*connected+; waic = 272.07
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + harmonic + mover*harmonic # sig: mover-, harmonic-, mover*harmonic+; waic = 267.82
    
    # explore spatial random effect (and removing kin.occ.inf)
    formula <- y ~ age + mom.serology + svi + hai + iva + mover + degree + mover*degree + f(spatial.field, model=spde) # sig: mover-, mover*degree+; ns:, degree-; waic = 262.48 !!! WINNER !!!
    formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + harmonic + mover*harmonic + f(spatial.field, model=spde) # sig: mover-, harmonic-, mover*harmonic+; waic = 262.56
    formula <- y ~ age + mom.serology + svi + hai + iva + mover + harmonic + mover*harmonic + f(spatial.field, model=spde) # sig: mover-, harmonic-, mover*harmonic+; waic = 262.58 !!! NEGLIGIBLE DIFFERENCE !!!
    
    # explore household random effect, but remove mom.serology - INCREASED WAIC
    formula <- y ~ age + svi + hai + iva + mover + degree + mover*degree + f(spatial.field, model=spde) + f(hhx.id, model = "iid") # sig: mover-, mover*degree+; ns: degree-; waic = 263.04
    formula <- y ~ age + svi + hai + iva + mover + kin.occ.inf + harmonic + mover*harmonic + f(spatial.field, model=spde) + f(hhx.id, model = "iid") # sig: mover-, harmonic-, mover*harmonic+; waic = 264.72
    
# explore interaction w/ iva
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + degree + iva*degree # sig: iva*degree-; ns: mover-, degree+; waic = 271.45
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + connected + iva*connected # sig: iva*connected-; ns: mover-, connected-; waic = 273.79
formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + harmonic + iva*harmonic # sig: iva*harmonic-; ns: mover-, harmonic+; waic = 259.36
    
    # explore spatial random effect (and removing kin.occ.inf)
    formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + degree + iva*degree + f(spatial.field, model=spde) # sig: iva+, iva*degree-, kin.occ.inf+; ns: degree+; waic = 265.21
    formula <- y ~ age + mom.serology + svi + hai + iva + mover + kin.occ.inf + harmonic + iva*harmonic + f(spatial.field, model=spde) # sig: iva+, iva*harmonic-, kin.occ.inf+; ns: harmonic+; waic = 256.34 !!! WINNER !!!
    formula <- y ~ age + mom.serology + svi + hai + iva + mover + harmonic + iva*harmonic + f(spatial.field, model=spde) # sig: iva+, iva*harmonic-, kin.occ.inf+; ns: harmonic+; waic = 257.18 !!! NEGLIGIBLE DIFFERENCE !!!
    
    # explore household random effect, but remove mom.serology - INCREASED WAIC
    formula <- y ~ age + svi + hai + iva + mover + kin.occ.inf + degree + iva*degree + f(spatial.field, model=spde) + f(hhx.id, model = "iid") # sig: iva+, iva*degree-; ns: degree+, kin.occ.inf+; waic = 265.93
    formula <- y ~ age + svi + hai + iva + mover + kin.occ.inf + harmonic + iva*harmonic + f(spatial.field, model=spde) + f(hhx.id, model = "iid") # sig: iva+, iva*harmonic-; ns: harmonic+, kin.occ.inf+; waic = 258.32
        
## fit the model
model <-inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
             control.family = list(link = "logit"),
             control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  
             control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), 
             verbose = FALSE)
summary(model)
   
# Final models of infection --------------------------------------------------
    
# load
library(INLA)
library(bayesplot)
library(ggplot2)

# list all formulas
formulas <- list(
  # infection0 = y ~ f(spatial.field, model=spde),
  infection1.0 = y ~ age + mom.serology + svi + hai + iva + mover + degree + mover*degree,
  infection1.1 = y ~ age + mom.serology + svi + hai + iva + mover + degree + mover*degree + f(spatial.field, model=spde), # sig: mover-, mover*degree+; ns:, degree-; waic = 262.48 !!! WINNER !!!
  infection1.2 = y ~ age + mom.serology + svi + hai + iva + mover + degree + f(spatial.field, model=spde),
  infection2.0 = y ~ age + mom.serology + svi + hai + iva + mover + harmonic + mover*harmonic, 
  infection2.1 = y ~ age + mom.serology + svi + hai + iva + mover + harmonic + mover*harmonic + f(spatial.field, model=spde), # sig: mover-, harmonic-, mover*harmonic+; waic = 262.56 !!! WINNER !!!
  infection2.2 = y ~ age + mom.serology + svi + hai + iva + mover + harmonic + f(spatial.field, model=spde),
  infection3.0 = y ~ age + mom.serology + svi + hai + iva + mover + harmonic + iva*harmonic, 
  infection3.1 = y ~ age + mom.serology + svi + hai + iva + mover + harmonic + iva*harmonic + f(spatial.field, model=spde), # sig: iva+, iva*harmonic-, kin.occ.inf+; ns: harmonic+; waic = 256.34 !!! WINNER !!!
  infection3.2 = y ~ age + mom.serology + svi + hai + iva + mover + harmonic + f(spatial.field, model=spde)
)

# prepare list to store results
results <- list()

# loop to fit and plot
for (i in seq_along(formulas)) {
      
    name <- names(formulas)[i]
    formula <- formulas[[i]]
    
    # fit model
    model <-inla(formula, 
                 data=inla.stack.data(stk, spde=spde), #extracts dataset from INLA stack object and tells INLA to use a stoch. partial differential equation object to represent the spatial random field
                 family= 'binomial', #specifies likelihood for the outcome variable
                 control.family = list(link = "logit"), #configures details of likelihood family: link="logit" means we are fitting a logistic regression
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE), #controls how the predictor is built: A=...supplies the projection matrix that maps latent effects (only applied if adjusting for spatial effects)
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
    
    # fixed effects posterior plot
      
        # get fixed effects names (minus intercept)
        names_fe <- names(model$marginals.fixed)[-1]
        
        # get samples for each fixed effect
        set.seed(4321)
        posterior_samples <- sapply(model$marginals.fixed, function(marg) {
          inla.rmarginal(10000, marg)
        })
        
        # plot
        g_posterior <- mcmc_areas(posterior_samples,
                                  pars = names_fe,
                                  prob = 0.95) +
          ggtitle(paste0("WAIC = ", round(model$waic$waic,2))) +
          theme_minimal(base_family = "Arial") +
          theme(plot.title = element_text(size = 16, face = "bold", family = "sans"), 
                axis.title = element_text(size = 14), 
                axis.text = element_text(size = 12))
        
        # save
        fe_count <- length(names(model$marginals.fixed)[-1])
        ggsave(
          filename = paste0("figures/inla/infection/model_", name, "_fe.png"),
          plot = g_posterior, width = 6, height = fe_count, dpi = 600
        )
      
        # collect model fit metrics
        metrics <- data.frame(
          regression = paste(deparse(formula), collapse = ""),  # store formula as a single string
          WAIC = model$waic$waic,
          pWAIC = model$waic$p.eff,
          DIC = model$dic$dic,
          pDIC = model$dic$p.eff,
          converged = ifelse(model$mode$mode.status == 0, TRUE, FALSE),
          n_params = nrow(model$summary.fixed)
        )
        
        # save results into list
        results[[name]] <- list(
          model = model,
          plot = g_posterior,
          metrics = metrics
        )
      
    # spatial effects plots (if present)
    if ("spatial.field" %in% names(model$summary.random)) {      
        
        # posterior of range
        
            # get marginals of range
            spde_result <- inla.spde2.result(model, "spatial.field", spde, do.transf = TRUE)
            post_range <- spde_result$marginals.range$range
            df_range <- data.frame(post_range)
            
            # plot
            g_range <- ggplot(df_range, aes(x = x, y = y)) +
              geom_line() +
              geom_vline(xintercept = 1500, color = "red", linetype = "dashed") +
              geom_vline(xintercept = inla.qmarginal(0.5, post_range), color = "blue", linetype = "dashed") +
              labs(x = "Spatial range (meters)",
                   y = "Posterior density",
                   title = "Posterior distribution of spatial range") +
              theme_minimal(base_size = 14) +
              theme(axis.line = element_blank())
            
            # save
            ggsave(
              filename = paste0("figures/inla/infection/model_", name, "_sfrange.png"),
              plot = g_range, width = 4, height = 4, dpi = 600
            )
        
        # project spatial field mean/sd
        
            # project mesh, mean, and sd
            gproj <- inla.mesh.projector(mesh, dims = c(300, 300))
            g.mean <- inla.mesh.project(gproj, model$summary.random$spatial.field$mean)
            g.sd   <- inla.mesh.project(gproj, model$summary.random$spatial.field$sd)
            
            # diverging palette for mean
            val_min <- min(model$summary.random$spatial.field$mean, na.rm = TRUE)
            val_max <- max(model$summary.random$spatial.field$mean, na.rm = TRUE)
            prop_neg <- abs(val_max) / (abs(val_min) + val_max)
            n = 200
            n_neg <- max(1, round(n * prop_neg))
            n_pos <- max(1, n - n_neg)
            at <- c(seq(val_min, 0, length.out = n_neg + 1),
                    seq(0, val_max, length.out = n_pos + 1)[-1])
            cols <- c(colorRampPalette(c("blue", "white"))(n_neg),
                      colorRampPalette(c("white", "red"))(n_pos))
            div_mean <- list(at = at, cols = cols)
            
            # red scale for sd
            at_sd <- seq(min(model$summary.random$spatial.field$sd, na.rm = TRUE),
                         max(model$summary.random$spatial.field$sd, na.rm = TRUE), length.out = 200)
            cols_sd <- colorRampPalette(c("red", "white"))(200)
            div_sd <- list(at = at_sd, cols = cols_sd)
            
            my_par <- list(par.xlab.text = list(cex = 1.2, font = 2),
                           par.ylab.text = list(cex = 1.2, font = 2))
            
            # plot
            g_sf <- plot_grid(
              levelplot(g.mean, scales = list(draw = FALSE), xlab = "", ylab = "", main = "mean",
                        region = TRUE, at = div_mean$at, col.regions = div_mean$cols,
                        par.settings = my_par),
              levelplot(g.sd, scales = list(draw = FALSE), xlab = "", ylab = "", main = "sd",
                        region = TRUE, at = div_sd$at, col.regions = div_sd$cols,
                        par.settings = my_par),
              nrow = 1
            )
        
            # save
            ggsave(
              filename = paste0("figures/inla/infection/model_", name, "_sfmap.png"),
              plot = g_sf, width = 8, height = 4, dpi = 600
            )
    }
}        

# export spatial field map as raster for processing in QGIS
    
    ## formula
    formula = y ~ age + mom.serology + svi + hai + iva + mover + degree + mover*degree + f(spatial.field, model=spde)
    formula = y ~ age + mom.serology + svi + hai + iva + mover + harmonic + iva*harmonic + f(spatial.field, model=spde)
    
    ## fit model
    model <-inla(formula, 
                 data=inla.stack.data(stk, spde=spde), #extracts dataset from INLA stack object and tells INLA to use a stoch. partial differential equation object to represent the spatial random field
                 family= 'binomial', #specifies likelihood for the outcome variable
                 control.family = list(link = "logit"), #configures details of likelihood family: link="logit" means we are fitting a logistic regression
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE), #controls how the predictor is built: A=...supplies the projection matrix that maps latent effects (only applied if adjusting for spatial effects)
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
    
    ## project spatial field mean/sd on regular grid
    gproj <- inla.mesh.projector(mesh, dims = c(300, 300))
    g.mean <- inla.mesh.project(gproj, model$summary.random$spatial.field$mean)
    g.sd   <- inla.mesh.project(gproj, model$summary.random$spatial.field$sd)

    ## create rasters from projected mean
    # note: terra expects rows from top to bottom, so we flip vertically
    r_mean <- rast(
      nrows = length(gproj$y),
      ncols = length(gproj$x),
      xmin = min(gproj$x),
      xmax = max(gproj$x),
      ymin = min(gproj$y),
      ymax = max(gproj$y),
      vals = as.vector(g.mean[nrow(g.mean):1, ])
    )
    
    r_sd <- rast(
      nrows = length(gproj$y),
      ncols = length(gproj$x),
      xmin = min(gproj$x),
      xmax = max(gproj$x),
      ymin = min(gproj$y),
      ymax = max(gproj$y),
      vals = as.vector(g.sd[nrow(g.sd):1, ])
    )
    
    ## assign CRS
    crs(r_mean) <- "EPSG:32721"
    crs(r_sd) <- "EPSG:32721"
    
    ## export as GeoTIFF
    # writeRaster(r_mean, "figures/maps/forQGIS/infection1_sf_mean.tif", overwrite = TRUE)
    # writeRaster(r_sd, "figures/maps/forQGIS/infection1_sf_sd.tif", overwrite = TRUE)
    writeRaster(r_mean, "figures/maps/forQGIS/infection2_sf_mean.tif", overwrite = TRUE)
    writeRaster(r_sd, "figures/maps/forQGIS/infection2_sf_sd.tif", overwrite = TRUE)

# export house locations used in the model as a point layer to add to QGIS map
    
    ## load library
    library(sf)
    library(terra)
    
    ## print formula and response variable
    print(formula) # formula = y ~ age + mom.serology + svi + hai + iva + mover + degree + mover*degree + f(spatial.field, model=spde)
                   # formula = y ~ age + mom.serology + svi + hai + iva + mover + harmonic + iva*harmonic + f(spatial.field, model=spde)
    print(resp_child)
    names(df_child)
    
    
    ## identify houses actually used in the infestation model based on variables above
    houses_used <- df_child[complete.cases(
      df_child[ ,c("serology.child", "age", "mom.serology", "svi", "hai", "iva", "mover", "degree", "harmonic", "x", "y" )]
      ), ]
    
    ## turn coordinates to sf points (point object)
    houses_sf <- st_as_sf(
      houses_used,
      coords = c("x", "y"),  
      crs = 32721                  # same CRS as your raster
    )
    
    ## export as geopackage (for QGIS)
    st_write(houses_sf,
      "figures/maps/forQGIS/houses_used_infection_models.gpkg",
      delete_dsn = TRUE
    )
              
# convert list of summary statistics for posteriors to table
    
    ## fixed effects
    fixed_table <- do.call(rbind, lapply(results, function(x) x$model$summary.fixed))
    fixed_table$kld <- NULL
    fixed_table$type <- "fixedeffects"
    fixed_table$order <- seq.int(nrow(fixed_table))
    
    ## hyperparameters
    hyperpar_table <- do.call(rbind, lapply(results, function(x) x$model$summary.hyperpar))
    hyperpar_table$type <- "hyperparameters"
    hyperpar_table$order <- seq.int(nrow(hyperpar_table))
    
    ## combine
    posteriors_table <- rbind(fixed_table, hyperpar_table)
        
        # save rownames as 'model' and 'parameter'
        library(tibble)
        library(stringr)
        posteriors_table <- posteriors_table %>%
          rownames_to_column("rowname") %>% 
          mutate(
            model     = str_extract(rowname, "^[^.]+\\.[^.]+"),   # everything up to 2nd dot
            parameter = str_remove(rowname, "^[^.]+\\.[^.]+\\.") # everything after 2nd dot
          ) %>%
          select(model, parameter, everything(), -rowname)
    
          # correct infection0 model and parameter name
          posteriors_table[1,"model"] <- "infection0"
          posteriors_table[1,"parameter"] <- "(Intercept)"
          
    ## convert list of model performance to table
    waic_table <- do.call(rbind, lapply(results, function(x) x$metrics))
        
        # save rowname (model name) as colname
        waic_table$model <- rownames(waic_table)
        rownames(waic_table) <- NULL 
    
    ## merge
    summary_table <- merge(posteriors_table, waic_table, by = "model")
    
    ## rename '0.5quant' to 'median'
    colnames(summary_table)[colnames(summary_table) == "0.5quant"] <- "median"
    
    ## sort rows
    summary_table <- summary_table[order(summary_table$model, summary_table$type, summary_table$order), ]
    
    ## reorder and remove some vars
    summary_table <- subset(summary_table, select = c("model", "regression", "parameter", "mean", "0.025quant", "median", "0.975quant", "mode", "sd", "WAIC"))
    
    ## save tables
    table_summary_infection <- summary_table
    # write.csv(table_summary_infection, "tables/table_summary_infection.csv")

# Post-hoc calculations --------------------------------------------------------

# Get 95% CrI for the joint posterior distribution of beta_mover_connected (the average effect of connected when mobility = 1)

    ## Why do we need the joint posterior?
    # Interactions are combinations of coefficients.
    # To calculate their posterior distributions, you need the joint posterior (i.e., the correlation b/w coefficients).
    # marginals.fixed() gives us the marginal posterior distribution for each coefficient, which assumes independence between them.
    # inla.posterior.sample() gives us posterior samples, whereby each draw represents a joint realization of all parameters simultaneously.
    # So if we want the joint posterior of beta_mover_connected (= beta_mobility + beta_connected + beta_mobility*connected),...
    # ...we can compute beta_mover_connected for each posterior draw to 
    # This preserves the joint covariance structure

    
# Infection Model 1) Interaction Mover*Degree

    ## get model results
    formula <- y ~ age + mom.serology + svi + hai + iva + mover + degree + mover*degree + f(spatial.field, model=spde)
    model <-inla(formula, 
             data=inla.stack.data(stk, spde=spde), #extracts dataset from INLA stack object and tells INLA to use a stoch. partial differential equation object to represent the spatial random field
             family= 'binomial', #specifies likelihood for the outcome variable
             control.family = list(link = "logit"), #configures details of likelihood family: link="logit" means we are fitting a logistic regression
             control.predictor=list(A=inla.stack.A(stk),compute=TRUE), #controls how the predictor is built: A=...supplies the projection matrix that maps latent effects (only applied if adjusting for spatial effects)
             control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
    
    ## extract posterior samples for the relevant coefficients: beta_connected and beta_mobility*connected
    samples <- inla.posterior.sample(n = 2000, model, seed=4321)   # draw 2k posterior samples
    
    ## function that computes beta for one-unit increase in degree among non-movers
    fun_mover_degree <- function() {
      degree
    }
        
        ## evaluate the function over the samples
        eval_mc <- inla.posterior.sample.eval(fun_mover_degree, samples)
        
        ## summarize
        quants <- quantile(eval_mc, probs = c(0.025, 0.50, 0.975))
        print(quants)
        
        ## summarize in odds scale
        or_draws <- exp(eval_mc)
        or_quants <- quantile(or_draws, probs = c(0.025, 0.50, 0.975))
        print(or_quants)
    
    ## function that computes beta_connected + beta_mobility*connected
    fun_effect_among_movers <- function() {
      degree + `mover:degree`
    }
        
        ## evaluate the function over the samples
        eval_mc <- inla.posterior.sample.eval(fun_effect_among_movers, samples)
        
        ## summarize
        quants <- quantile(eval_mc, probs = c(0.025, 0.50, 0.975))
        print(quants)
        
        ## summarize in odds scale
        or_draws <- exp(eval_mc)
        or_quants <- quantile(or_draws, probs = c(0.025, 0.50, 0.975))
        print(or_quants)
        
     ## Check for multicollinearity by fitting a standard glm
        
        # glm for diagnostics only (excluding interaction terms)
        glm_model <- glm(serology.child ~ age + mom.serology + svi + hai + iva + mover + degree,
                         family = binomial, data = df)
        
        # vifs
        library(car)
        vif(glm_model) #all vifs between 1-1.13: suggests no inflation of variance due to correlation.
    
     ## Check for nonlinear and constrained effect
        
        formula <- y ~ age + mom.serology + svi + hai + iva + mover + degree + mover*degree + f(spatial.field, model=spde) # waic = 262.5
        formula <- y ~ f(age, model = "rw1") + mom.serology + svi + hai + iva + mover + degree + mover*degree + f(spatial.field, model=spde) # waic = 267.4
        formula <- y ~ age + mom.serology + svi + hai + f(iva, model = "rw1") + mover + degree + mover*degree + f(spatial.field, model=spde) # waic = 255.2
        formula <- y ~ f(age, model = "clinear", range = c(0, Inf)) + mom.serology + svi + hai + iva + mover + degree + mover*degree + f(spatial.field, model=spde) # 269.1
        formula <- y ~ age + mom.serology + svi + hai + f(iva, model = "clinear", range = c(0, Inf)) + mover + degree + mover*degree + f(spatial.field, model=spde) #273.4
        
        model <-inla(formula, 
                     data=inla.stack.data(stk, spde=spde), #extracts dataset from INLA stack object and tells INLA to use a stoch. partial differential equation object to represent the spatial random field
                     family= 'binomial', #specifies likelihood for the outcome variable
                     control.family = list(link = "logit"), #configures details of likelihood family: link="logit" means we are fitting a logistic regression
                     control.predictor=list(A=inla.stack.A(stk),compute=TRUE), #controls how the predictor is built: A=...supplies the projection matrix that maps latent effects (only applied if adjusting for spatial effects)
                     control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
        summary(model$waic$waic)
            
            # visualize results for age as rw1
            age_effect <- model$summary.random$age
            
                # save plot
                png("figures/inla/infection/model_infection1.1_agerw1.png", width = 8, height = 6, units = "in", res = 600, bg = "transparent")
                ggplot(age_effect, aes(x = ID, y = mean)) +
                  geom_line(color = "blue") +
                  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.3) +
                  labs(title = "Estimated Age Effect on T. cruzi infection (rw1)",
                       x = "Age", y = "Effect Size") +
                  theme_minimal()
                dev.off()
            
            # visualize results for iva as rw1
            iva_effect <- model$summary.random$iva
            
                # save plot
                png("figures/inla/infection/model_infection1.1_ivarw1.png", width = 8, height = 6, units = "in", res = 600, bg = "transparent")
                ggplot(iva_effect, aes(x = ID, y = mean)) +
                  geom_line(color = "blue") +
                  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.3) +
                  labs(title = "Estimated Infected Vector Abundance Effect on T. cruzi infection (rw1)",
                       x = "Abundance of infected vector", y = "Effect Size") +
                  theme_minimal()
                dev.off()
            
# Infection Model 2) Interaction Mover*Harmonic
        
        ## get model results
        formula <- y ~ age + mom.serology + svi + hai + iva + mover + harmonic + mover*harmonic + f(spatial.field, model=spde)
        model <-inla(formula, 
                     data=inla.stack.data(stk, spde=spde), #extracts dataset from INLA stack object and tells INLA to use a stoch. partial differential equation object to represent the spatial random field
                     family= 'binomial', #specifies likelihood for the outcome variable
                     control.family = list(link = "logit"), #configures details of likelihood family: link="logit" means we are fitting a logistic regression
                     control.predictor=list(A=inla.stack.A(stk),compute=TRUE), #controls how the predictor is built: A=...supplies the projection matrix that maps latent effects (only applied if adjusting for spatial effects)
                     control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
        
        ## extract posterior samples for the relevant coefficients: beta_connected and beta_mobility*connected
        set.seed(4321)
        samples <- inla.posterior.sample(n = 2000, model)   # draw 2k posterior samples
        
        ## function that computes beta_mover_harmonic
        fun_mover_harmonic <- function() {
          mover + harmonic + `mover:harmonic`
        }
        
            ## evaluate the function over the samples
            eval_mc <- inla.posterior.sample.eval(fun_mover_harmonic, samples)
            
            ## summarize
            quants <- quantile(eval_mc, probs = c(0.025, 0.50, 0.975))
            print(quants)
            
            ## summarize in odds scale
            or_draws <- exp(eval_mc)
            or_quants <- quantile(or_draws, probs = c(0.025, 0.50, 0.975))
            print(or_quants)
        
        ## function that computes beta_connected + beta_mobility*harmonic
        fun_effect_among_movers <- function() {
          harmonic + `mover:harmonic`
        }
            
            ## evaluate the function over the samples
            eval_mc <- inla.posterior.sample.eval(fun_effect_among_movers, samples)
            
            ## summarize
            quants <- quantile(eval_mc, probs = c(0.025, 0.50, 0.975))
            print(quants)
            
            ## summarize in odds scale
            or_draws <- exp(eval_mc)
            or_quants <- quantile(or_draws, probs = c(0.025, 0.50, 0.975))
            print(or_quants)

        
# Infection Model 3) Interaction IVA*Harmonic
        
        ## get model results
        formula <- y ~ age + mom.serology + svi + hai + iva + mover + degree + iva*degree + f(spatial.field, model=spde)
            
        formula <- y ~ age + mom.serology + svi + hai + iva + mover + harmonic + iva*harmonic + f(spatial.field, model=spde)
        model <-inla(formula, 
                     data=inla.stack.data(stk, spde=spde), #extracts dataset from INLA stack object and tells INLA to use a stoch. partial differential equation object to represent the spatial random field
                     family= 'binomial', #specifies likelihood for the outcome variable
                     control.family = list(link = "logit"), #configures details of likelihood family: link="logit" means we are fitting a logistic regression
                     control.predictor=list(A=inla.stack.A(stk),compute=TRUE), #controls how the predictor is built: A=...supplies the projection matrix that maps latent effects (only applied if adjusting for spatial effects)
                     control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
        
        ## extract posterior samples for the relevant coefficients: beta_connected and beta_mobility*connected
        set.seed(4321)
        samples <- inla.posterior.sample(n = 2000, model)   # draw 2k posterior samples
        
        ## function that computes beta_iva_connected
        fun_iva_harmonic <- function() {
          iva + harmonic + `iva:harmonic`
        }
            
            ## evaluate the function over the samples
            eval_mc <- inla.posterior.sample.eval(fun_iva_harmonic, samples)
            
            ## summarize
            quants <- quantile(eval_mc, probs = c(0.025, 0.50, 0.975))
            print(quants)
            
            ## summarize in odds scale
            or_draws <- exp(eval_mc)
            or_quants <- quantile(or_draws, probs = c(0.025, 0.50, 0.975))
            print(or_quants)
        
        ## function that computes beta_connected + beta_mobility*connected
        fun_effect_among_iva <- function() {
          harmonic + `iva:harmonic`
        }
            
            ## evaluate the function over the samples
            eval_mc <- inla.posterior.sample.eval(fun_effect_among_iva, samples)
            
            ## summarize
            quants <- quantile(eval_mc, probs = c(0.025, 0.50, 0.975))
            print(quants)
            
            ## summarize in odds scale
            or_draws <- exp(eval_mc)
            or_quants <- quantile(or_draws, probs = c(0.025, 0.50, 0.975))
            print(or_quants)
            
# Combine summary tables -------------------------------------------------------

# combine tables
table_summary_allmodels <- rbind(table_summary_infestation, 
                                 table_summary_reinfestation, 
                                 table_summary_iva,
                                 table_summary_infection)
    
# drop intercept observations
table_summary_allmodels <- subset(table_summary_allmodels, parameter != "(Intercept)")
    
# rename observations
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "svi"] <- "Social vulnerability"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "hai"] <- "Host availability"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "mobility"] <- "Mobility"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "iva"] <- "Infected vector abundance"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "kin.occ.inf"] <- "Presence of T. infestans-infested kin household"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "connected"] <- "Connected"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "degree"] <- "Degree"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "harmonic"] <- "Harmonic"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "age"] <- "Age"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "mom.serology"] <- "Mother's serostatus"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "mover"] <- "Mobility"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "Range for spatial.field"] <- "Range for spatial field"
table_summary_allmodels$parameter[table_summary_allmodels$parameter == "Stdev for spatial.field"] <- "SD for spatial field"

# rename columns
colnames(table_summary_allmodels) <- c("Model", "Regression", "Parameter", "Mean", "2.5% percentile", "Median", "97.5% percentile", "Mode", "SD", "WAIC")

# save table
write.csv(table_summary_allmodels, "tables/table_summary_allmodels.csv")

# # save table from sensitivity analysis
#     
#     # spde <- spde_weak_rho
#     write.csv(table_summary_allmodels, "figures/inla/sensitivity analysis/table_inla_spde_weak_rho.csv")
#     
#     # spde <- spde_short_rho
#     write.csv(table_summary_allmodels, "figures/inla/sensitivity analysis/table_inla_spde_short_rho.csv")
#     
#     # spde <- spde_long_rho
#     write.csv(table_summary_allmodels, "figures/inla/sensitivity analysis/table_inla_spde_long_rho.csv")
#     
#     # spde <- spde_low_var
#     write.csv(table_summary_allmodels, "figures/inla/sensitivity analysis/table_inla_spde_low_var.csv")
#     
#     # spde <- spde_high_var
#     write.csv(table_summary_allmodels, "figures/inla/sensitivity analysis/table_inla_spde_high_var.csv")
#     
#     # spde <- spde_cons_var
#     write.csv(table_summary_allmodels, "figures/inla/sensitivity analysis/table_inla_spde_cons_var.csv")


################################################################################
### Investigate Constrained and Nonlinear Effects of Predictors
################################################################################

## see "4.0_PopModel(Amratia).R")

## create empty df for saving results
results <- data.frame("model" = character(),
                      "dic" = numeric())

# Age --------------------------------------------------------------------------

    ## let's explore (1) setting prior estimate, (2) setting gamma prior, (3) nonlinear w/ rw1, (4) constrained linear effect ## clinear - What if the relationship should be non-zero? Let's constrain covariate coefficient by setting "range = c(a, Inf)"

    ## (1) set prior estimate
    formula <- y ~ age 
    model <- inla(formula,
                  data = inla.stack.data(stk, spde = spde),
                  family = 'binomial',
                  control.family = list(link = "logit"),
                  control.fixed = list(
                    mean = list(Intercept = 0, age = 0.0677),
                    prec = list(Intercept = 1e-4, age = 1 / (0.0047684012812^2))
                  ),
                  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
    
    ## save results
    i = 1
    results[i,c("model")] <- "age.prior"
    results[i,c("mean")]  <- model$summary.fixed$mean[[2]]
    results[i,c("mode")]  <- model$summary.fixed$mode[[2]]
    results[i,c("ll")]    <- model$summary.fixed$`0.025quant`[[2]]
    results[i,c("ul")]    <- model$summary.fixed$`0.975quant`[[2]]
    results[i,c("dic")]   <- model$dic$dic
    
    ## (2) set gamma prior
    alpha <- 5
    lambda <- 1
    mu <- alpha / lambda  # mean
    sigma2 <- alpha / (lambda^2)  # variance
    
    formula <- y ~ age 
    inla  <- inla(formula,
                  data = inla.stack.data(stk, spde = spde),
                  family = "binomial", control.family = list(link = "logit"),
                  control.fixed = list(mean=list(age=mu),
                                       prec = list(age = 1/sigma2)), # precision is 1/variance
                  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
    )
    
    ## (3) nonlinear w/ rw1
    formula <- y ~ f(age, model = "rw1") + f(spatial.field, model=spde)
    
    model <-inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
                 # Ntrials = n,
                 control.family = list(link = "logit"),
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log
    
    summary(model)
    
    # visualize results
    age_effect <- model$summary.random$age
    
    ggplot(age_effect, aes(x = ID, y = mean)) +
      geom_line(color = "blue") +
      geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.3) +
      labs(title = "Estimated Age Effect on T. cruzi infection (clinear)",
           x = "Age", y = "Effect Size") +
      theme_minimal()
    
    ## (4) constrained linear effect
    formula <- y ~ f(age, model = "clinear", range = c(0, Inf)) + f(spatial.field, model=spde)
    
    model <-inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
                 # Ntrials = n,
                 control.family = list(link = "logit"),
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log
    
    summary(model)
    
        ## save results
        i = 1
        results[i,c("model")] <- "age.prior"
        results[i,c("mean")]  <- model$summary.fixed$mean[[2]]
        results[i,c("mode")]  <- model$summary.fixed$mode[[2]]
        results[i,c("ll")]    <- model$summary.fixed$`0.025quant`[[2]]
        results[i,c("ul")]    <- model$summary.fixed$`0.975quant`[[2]]
        results[i,c("dic")]   <- model$dic$dic
    
            
    # visualize results
    age_effect <- model$summary.random$age
    
    ggplot(age_effect, aes(x = ID, y = mean)) +
      geom_line(color = "blue") +
      geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.3) +
      labs(title = "Estimated Age Effect on T. cruzi infection (clinear)",
           x = "Age", y = "Effect Size") +
      theme_minimal()

# Infected vector abundance --------------------------------------------------------------------------
    
    ## create empty df for saving results
    results <- data.frame("model" = character(),
                          "mean" = numeric(), 
                          "mode" = numeric(), 
                          "ll" = numeric(), 
                          "ul" = numeric(), 
                          "dic" = numeric())
    
    ## let's explore (1) nonlinear w/ rw1, (2) constrained linear effect ## clinear - What if the relationship should be non-zero? Let's constrain covariate coefficient by setting "range = c(a, Inf)"
    
    ## (1) nonlinear w/ rw1
    formula <- y ~ f(iva, model = "rw1") + f(spatial.field, model=spde)
    model <-inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
                 # Ntrials = n,
                 control.family = list(link = "logit"),
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log
    
    # visualize results
    iva_effect <- model$summary.random$iva
    
    ggplot(iva_effect, aes(x = ID, y = mean)) +
      geom_line(color = "blue") +
      geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.3) +
      labs(title = "IVA Effect on T. cruzi infection (clinear)",
           x = "Age", y = "Effect Size") +
      theme_minimal()
    
    ## (4) constrained linear effect
    formula <- y ~ f(iva, model = "clinear", range = c(0, Inf)) + f(spatial.field, model=spde)
    model <-inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
                 # Ntrials = n,
                 control.family = list(link = "logit"),
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log
    
        # visualize results
        iva_effect <- model$summary.random$iva
        
        ggplot(iva_effect, aes(x = ID, y = mean)) +
          geom_line(color = "blue") +
          geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.3) +
          labs(title = "IVA Effect on T. cruzi infection (clinear)",
               x = "Age", y = "Effect Size") +
          theme_minimal()
    
## Incorporate age w/ clinear constraint

    formula <- y ~ f(age, model = "clinear", range = c(0, Inf)) + iva + svi + hai + mom.serology + degree + mover + degree*mover + f(spatial.field, model=spde) #no divergence warning

    model <-inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
                 # Ntrials = n,
                 control.family = list(link = "logit"),
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log
    
    summary(model) # DIC increased!
    
## Incorporate iva w/ rw1 constraint
      
    formula <- y ~ f(age, model = "clinear", range = c(0, Inf)) + f(iva, model = "rw1") + svi + hai + mom.serology + degree + mover + degree*mover + f(spatial.field, model=spde) #no divergence warning
    
    model <-inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
                 # Ntrials = n,
                 control.family = list(link = "logit"),
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log
    
    summary(model)
    
    # visualize results
    age_effect <- model$summary.random$age
    
    ggplot(age_effect, aes(x = ID, y = mean)) +
      geom_line(color = "blue") +
      geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.3) +
      labs(title = "Estimated Age Effect on T. cruzi infection (clinear)",
           x = "Age", y = "Effect Size") +
      theme_minimal()
    
    # visualize results
    iva_effect <- model$summary.random$iva
    
    ggplot(iva_effect, aes(x = ID, y = mean)) +
      geom_line(color = "blue") +
      geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.3) +
      labs(title = "Estimated IVA Effect on T. cruzi infection (clinear)",
           x = "Infected vector abundance", y = "Effect Size") +
      theme_minimal()

### SUMMARY ###
# incorporating age and iva as potentially non-linear random effects showed that both variables were mostly linear in their effect
# clinear led to increase in DIC
# keep age and iva as is

    
#############################
# ADDITIONAL POST-HOC TESTS? #
#############################
    
# Are mobile households younger? And are younger households more socially vulnerable and thus possibly be more reliant on their kin relations?
    temp <- df_kin[!is.na(df_kin$mobility),]
    
    ggplot(temp, aes(x = mobility, y = age)) +
           geom_boxplot() +
           labs(x = "Mobility", y = "Age") +
           theme_minimal()
    wilcox.test(age ~ mobility, data = temp)
  
# Are mobile households more socially vulnerable? Are younger households more socially vulnerable?
    ggplot(temp, aes(x = mobility, y = svi)) +
      geom_boxplot() +
      labs(x = "Mobility", y = "Social vulnerability") +
      theme_minimal()
    wilcox.test(svi ~ mobility, data = temp)    
    
    cor.test(temp$age, temp$svi, 
             method = "spearman", use = "complete.obs")
    
    
###############################
# SAVE FILES FOR RISK MAPPING #
###############################
    
formulas <- list(
  formula_inf = y ~ svi + hai + mobility + connected + mobility*connected +
    f(spatial.field, model = spde),   # BEST FIT
  
  formula_reinf = y ~ mobility + connected,
  
  formula_infect1 = y ~ age + mom.serology + svi + hai + iva + mover +
    degree + mover*degree +
    f(spatial.field, model = spde),
  
  formula_infect2 = y ~ age + mom.serology + svi + hai + iva + mover +
    harmonic + iva*harmonic +
    f(spatial.field, model = spde)
)
    
save(
  spat_inf, # spatial variables
  resp_inf,
  resp_reinf,
  resp_iva,
  resp_child,
  covar_inf,
  covar_child,
  df_inf,
  df_child,
  boundary,
  Mesh,
  spde,
  w.spde,
  stk_inf,
  stk_reinf,
  stk_child,
  formulas,
  file = "data/INLAModel_ObjectsForMapping.RData"
)
