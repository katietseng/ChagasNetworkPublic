##############################################################
##  Project: Chagas Network
##  Title: 2_GLMModels
##  author: Katie Tseng (katie.tseng@wsu.edu)
##############################################################

################################################################################
### PREP DATA
################################################################################

## clean environment
rm(list=ls()) 
graphics.off()

## set directory
setwd("/Users/~/Desktop/ChagasNetworkPublic")

## load data
load('data/CleanData.RData')

## create df for saving GLM results
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

## create df for saving GLMM results (random effects)
glmm_results <- glm_results

## library
library(glmm)

## gen kin network infestation variable
df_kin$kin.occ.inf <- ifelse(is.nan(df_kin$kin.inf.prop), NA,
                             ifelse(df_kin$kin.inf.prop > 0, 1, 0))

## transform to factor
df_kin$kin.occ.inf <- as.factor(df_kin$kin.occ.inf)

## df for infestation models
df_model_v.inf <- df_kin[!is.na(df_kin$v.inf),]


## df for child seropositivity models
df_model_child <- df_dem
names(df_model_child)

# get additional hhx-level variables
df_add <- subset(df_kin, select = c(hhx.id, 
                                    v.inf, 
                                    connected, degree, degree.cat, harmonic, betweenness,
                                    kin.svi.med, kin.svi.relmed,
                                    kin.iva.med, kin.inf.prop,
                                    kin.occ.inf))

# merge
df_model_child <- merge(df_model_child, df_add, by.x = "hhx.id", all.x = TRUE)
rm(df_add)

# checks
    
    ## how many households have at least one seropositive person?
    temp <- df_model_child[!is.na(df_model_child$serology),] #1368
    temp <- aggregate(serology ~ hhx.id, data = temp, sum) #303
    temp$occ.serology <- ifelse(temp$serology>0, 1, 0)
    mean(temp$occ.serology) #71%
    
    ## how many households have at least one seropositive person?
    temp <- df_model_child[!is.na(df_model_child$serology.child),]
    temp <- aggregate(serology.child ~ hhx.id, data = temp, sum) #226
    temp$occ.serology.child <- ifelse(temp$serology.child>0, 1, 0)
    mean(temp$occ.serology.child) #22%
    rm(temp)
    
# drop those w/ NA serology.child data
df_model_child <- df_model_child[!is.na(df_model_child$serology.child),] #705

################################################################################
### PART 1: BIVARIATE ANALYSIS
################################################################################
        
### GLM INFESTATION (household level) ###
    
# Get df for model at household level 
df_model <- df_model_v.inf
# df_model <- df_model_v.inf[df_model_v.inf$mobility == "non-mover",]
# df_model <- df_model_v.inf[df_model_v.inf$mobility == "mover",]

# Covariates to loop through
covars <- c("degree", 
             "degree.cat",
             "connected", "harmonic", "betweenness",
             "kin.svi.med", "kin.svi.relmed",
             "kin.iva.med", "kin.occ.inf")

# Initialize an empty list to store models
glm_list <- list()

# GLM loop
    
    # glm (save results to list)
    for (covar in covars) {
      
        ## create formula dynamically
        formula <- as.formula(paste("v.inf ~", covar))
        
        ## glm
        glm_list[[covar]] <- glm(formula,                   
                                 family = binomial(link = "logit"),
                                 data = df_model)
    }
    
    ## transform results to table
    for (covar in covars) {
      
      results <- glm_list[[covar]]
      
      ## If covariate is factor...
      if (is.factor(df_model[[covar]])) {
        
        levels_covar <- levels(df_model[[covar]])
        
        # then loop through levels (excluding reference level)
        for (lvl in levels_covar[-1]) {
          coef_name <- paste0(covar, lvl)
          
          if (coef_name %in% names(coef(results))) {
            ci <- confint(results, level = 0.95)
            new_row <- data.frame(
              depvar = "v.inf",
              covar = covar,
              level = lvl,
              estimate = coef(results)[coef_name],
              POR = sprintf("%.2f", exp(coef(results)[coef_name])),
              LL = sprintf("%.2f", exp(ci[coef_name, "2.5 %"])),
              UL = sprintf("%.2f", exp(ci[coef_name, "97.5 %"])),
              P = summary(results)$coefficients[coef_name, "Pr(>|z|)"],
              AIC = AIC(results),
              nobs = nobs(results),
              ci_method = "profile"
            )
            glm_results <- rbind(glm_results, new_row)
          }
        }
        
      } else {
        
        # For numeric (non-factor) covariates
        coef_name <- covar
        if (coef_name %in% names(coef(results))) {
          ci <- confint(results, level = 0.95)
          new_row <- data.frame(
            depvar = "v.inf",
            covar = covar,
            level = NA,
            estimate = coef(results)[coef_name],
            POR = sprintf("%.2f", exp(coef(results)[coef_name])),
            LL = sprintf("%.2f", exp(ci[coef_name, "2.5 %"])),
            UL = sprintf("%.2f", exp(ci[coef_name, "97.5 %"])),
            P = summary(results)$coefficients[coef_name, "Pr(>|z|)"],
            AIC = AIC(results),
            nobs = nobs(results),
            ci_method = "profile"
          )
          glm_results <- rbind(glm_results, new_row)
        }
      }
    }
    
    ## drop row names
    row.names(glm_results) <- NULL

    
### GLM CHILD SEROPOSITIVITY (individual-level) ###

# get df for model at individual-level
df_model <- df_model_child
# df_model <- df_model_child[df_model_child$mover=="Non-mover",]
# df_model <- df_model_child[df_model_child$mover=="Mover",]

# Covariates to loop through
covars <- c("degree", 
            "degree.cat",
            "connected", "harmonic", "betweenness",
            "kin.svi.med", "kin.svi.relmed",
            "kin.iva.med", "kin.occ.inf")
    
# Initialize an empty list to store models
glm_list <- list()
    
# GLM loop
    
    # glm (save results to list)
    for (covar in covars) {
      
      ## create formula dynamically
      formula <- as.formula(paste("serology.child ~", covar))
      
      ## glm
      glm_list[[covar]] <- glm(formula,                   
                               family = binomial(link = "logit"),
                               data = df_model)
    }
    
    ## transform results to table
    for (covar in covars) {
      
      results <- glm_list[[covar]]
      
      ## If covariate is factor...
      if (is.factor(df_model[[covar]])) {
        
        levels_covar <- levels(df_model[[covar]])
        
        # then loop through levels (excluding reference level)
        for (lvl in levels_covar[-1]) {
          coef_name <- paste0(covar, lvl)
          
          if (coef_name %in% names(coef(results))) {
            ci <- confint(results, level = 0.95)
            new_row <- data.frame(
              depvar = "serology.child",
              covar = covar,
              level = lvl,
              estimate = coef(results)[coef_name],
              POR = sprintf("%.2f", exp(coef(results)[coef_name])),
              LL = sprintf("%.2f", exp(ci[coef_name, "2.5 %"])),
              UL = sprintf("%.2f", exp(ci[coef_name, "97.5 %"])),
              P = summary(results)$coefficients[coef_name, "Pr(>|z|)"],
              AIC = AIC(results),
              nobs = nobs(results),
              ci_method = "profile"
            )
            glm_results <- rbind(glm_results, new_row)
          }
        }
        
      } else {
        
        # For numeric (non-factor) covariates
        coef_name <- covar
        if (coef_name %in% names(coef(results))) {
          ci <- confint(results, level = 0.95)
          new_row <- data.frame(
            depvar = "serology.child",
            covar = covar,
            level = NA,
            estimate = coef(results)[coef_name],
            POR = sprintf("%.2f", exp(coef(results)[coef_name])),
            LL = sprintf("%.2f", exp(ci[coef_name, "2.5 %"])),
            UL = sprintf("%.2f", exp(ci[coef_name, "97.5 %"])),
            P = summary(results)$coefficients[coef_name, "Pr(>|z|)"],
            AIC = AIC(results),
            nobs = nobs(results),
            ci_method = "profile"
          )
          glm_results <- rbind(glm_results, new_row)
        }
      }
    }
    
    ## add significance
    glm_results$significance = ifelse(glm_results$P < 0.001, "p<0.001",
                                         ifelse(glm_results$P < 0.01, "p<0.01",
                                                ifelse(glm_results$P < 0.05, "p<0.05", "ns")))
    
    
    ## drop row names
    row.names(glm_results) <- NULL
    
    # Save tables
    write.csv(glm_results, "tables/glmm/bivariate_glm_results.csv") 
    # write.csv(glm_results, "tables/glmm/bivariate_glm_nonmovers_results.csv") 
    # write.csv(glm_results, "tables/glmm/bivariate_glm_movers_results.csv") 
    
### GLMM CHILD SEROPOSITIVITY (w/ Random Effects) ###
    
# Library
library(lme4)    

# get data
df_model <- df_model_child
    
# Initialize an empty list to store models
glmm_list <- list()

## create df for saving GLMM results
glmm_results <- data.frame(depvar = character(),
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

# GLMM loop 

    # Covariates to loop through
    covars <- c("degree", 
                "degree.cat",
                "connected", "harmonic", "betweenness",
                "kin.svi.med", "kin.svi.relmed",
                "kin.iva.med", "kin.occ.inf")

    # glmm (save results to list)
    for (covar in covars) {
      
      ## create formula dynamically
      formula <- as.formula(paste("serology.child ~", covar, " + (1 | hhx.id)"))
      message("Running model for: ", "serology.child ~", covar)
      
      ## glm
      glmm_list[[covar]] <- glmer(formula,                   
                                  family = binomial(link = "logit"),
                                  data = df_model)
    }

    ### NOTE: error in calculating CIs require we calculate CIs one by one and using Wald's method if needed ###

    # Covariates to loop through
    covars <- c(
                # "degree"
                # "degree.cat"
                # "connected"
                # "harmonic"
                # "betweenness"
                # "kin.svi.med"
                # "kin.svi.relmed"
                # "kin.iva.med"
                "kin.occ.inf"
                )

    ## transform results to table
    for (covar in covars) {
          
      results <- glmm_list[[covar]]
          
      ## If covariate is factor...
      if (is.factor(df_model[[covar]])) {
        
          levels_covar <- levels(df_model[[covar]])
          
          # then loop through levels (excluding reference level)
          for (lvl in levels_covar[-1]) {
          fixef_name <- paste0(covar, lvl)
          
          if (fixef_name %in% names(fixef(results))) {
            # ci <- confint(results, level = 0.95) #degree.cat|
            ci <- confint(results, parm = fixef_name, method = "Wald") #connected|kin.occ.inf
            new_row <- data.frame(
              depvar = "serology.child",
              covar = covar,
              level = lvl,
              estimate = fixef(results)[fixef_name],
              POR = sprintf("%.2f", exp(fixef(results)[fixef_name])),
              # LL = sprintf("%.2f", exp(ci[fixef_name, "2.5 %"])),
              # UL = sprintf("%.2f", exp(ci[fixef_name, "97.5 %"])),
              LL = sprintf("%.2f", exp(ci[1, 1])),
              UL = sprintf("%.2f", exp(ci[1, 2])),
              P = summary(results)$coefficients[fixef_name, "Pr(>|z|)"],
              AIC = AIC(results),
              nobs = nobs(results),
              # ci_method = "profile"
              ci_method = "wald"
            )
            glmm_results <- rbind(glmm_results, new_row)
          }
        }
         
      } else {
        # For numeric (non-factor) covariates
        fixef_name <- covar
        if (fixef_name %in% names(fixef(results))) {
          # ci <- confint(results, level = 0.95)
          ci <- confint(results, parm = fixef_name, method = "Wald") #degree|harmonic|betweenness|kin.svi.med|kin.svi.relmed|kin.iva.med
          new_row <- data.frame(
            depvar = "serology.child",
            covar = covar,
            level = NA,
            estimate = fixef(results)[fixef_name],
            POR = sprintf("%.2f", exp(fixef(results)[fixef_name])),
            # LL = sprintf("%.2f", exp(ci[fixef_name, "2.5 %"])),
            # UL = sprintf("%.2f", exp(ci[fixef_name, "97.5 %"])),
            LL = sprintf("%.2f", exp(ci[1, 1])),
            UL = sprintf("%.2f", exp(ci[1, 2])),
            P = summary(results)$coefficients[fixef_name, "Pr(>|z|)"],
            AIC = AIC(results),
            nobs = nobs(results),
            # ci_method = "profile"
            ci_method = "wald"
          )
          glmm_results <- rbind(glmm_results, new_row)
        }
      }
    }
   
    ## add significance
    glmm_results$significance = ifelse(glmm_results$P < 0.001, "p<0.001",
                                      ifelse(glmm_results$P < 0.01, "p<0.01",
                                             ifelse(glmm_results$P < 0.05, "p<0.05", "ns")))
    
    # drop row names
    row.names(glmm_results) <- NULL
    
    # Save tables
    write.csv(glmm_results, "tables/glmm/bivariate_glmm_results.csv")
    
### SUMMARY ###
# v.inf ~ kin.svi.med (p<0.001)
# v.inf ~ kin.occ.inf (p<0.01)
# serology.child ~ degree.cat"1-2" (p<0.05)
# serology.child ~ connected (p<0.05)
# serology.child ~ kin.occ.inf (p<0.05)
# All random effects models were non-significant

    
### GLM MOBILITY (household level) ###
    
    # Get df for model at household level 
    df_model <- df_model_v.inf
    # df_model <- df_model_v.inf[df_model_v.inf$mobility == "non-mover",]
    # df_model <- df_model_v.inf[df_model_v.inf$mobility == "mover",]
    
    # Covariates to loop through
    covars <- c("degree", 
                "degree.cat",
                "connected", "harmonic", "betweenness",
                "kin.svi.med", "kin.svi.relmed",
                "kin.iva.med", "kin.occ.inf")

# glm: mobility ~ degree
    
    results <- glm(mobility ~ degree, 
                   family = binomial(link = "logit"),
                   data = df_model)
    summary(results) #degree -> ns
    
# glm: mobility ~ degree.cat
    
    results <- glm(mobility ~ degree.cat,
                   family = binomial(link = "logit"),
                   data = df_model)
    summary(results) #degree.cat1-2 -> ns
                     #degree.cat3+  -> ns
    
# glm: mobility ~ connected
    results <- glm(mobility ~ connected,
                   family = binomial(link = "logit"),
                   data = df_model)
    summary(results) #connected -> ns
    
# glm: mobility ~ harmonic
    results <- glm(mobility ~ harmonic,
                   family = binomial(link = "logit"),
                   data = df_model)
    summary(results) #harmonic -> ns
    
# glm: mobility ~ betweenness
    results <- glm(mobility ~ betweenness,
                   family = binomial(link = "logit"),
                   data = df_model)
    summary(results) #betweenness -> ns
    
################################################################################
### PART 2: MULTIVARIATE ANALYSIS
################################################################################

# Create df for saving GLM results
glm_results <- data.frame(sample = character(),
                          formula = character(),
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
    

    
# GLM: infestation
    
    # Get df for model at household level 
    df_model <- df_model_v.inf
    
    # Get dfs for movers vs. non-movers
    df_movers <- df_model[df_model$mobility=="mover" & !is.na(df_model$mobility),]
    df_nonmovers <- df_model[df_model$mobility=="non-mover" & !is.na(df_model$mobility),]
    
    # v.inf ~ degree...
        
        ## formula 
        formula <- as.formula(paste("v.inf ~ degree + svi + hai"))
        
        ## glm
        glm1 <- glm(formula, family = binomial(link = "logit"), data = df_model)
        glm2 <- glm(formula, family = binomial(link = "logit"), data = df_movers)
        glm3 <- glm(formula, family = binomial(link = "logit"), data = df_nonmovers)
        
        ## summarize "v.inf ~ degree + svi + hai"
        summary(glm1) # all:    -degree  |+svi***|+hai***; aic = 423.22
        summary(glm2) # movers: +degree  |+svi.  |+hai   ; aic = 94.67 
        summary(glm3) # nonmov: -degree  |+svi***|+hai** ; aic = 277.72
        
    # v.inf ~ connected...
        
        ## formula 
        formula <- as.formula(paste("v.inf ~ connected + svi + hai"))
        
        ## glm
        glm1 <- glm(formula, family = binomial(link = "logit"), data = df_model)
        glm2 <- glm(formula, family = binomial(link = "logit"), data = df_movers)
        glm3 <- glm(formula, family = binomial(link = "logit"), data = df_nonmovers)
        
        ## summarize "v.inf ~ connected + svi + hai"
        summary(glm1) # all:    -connected  |+svi***|+hai***; aic = 422.88
        summary(glm2) # movers: +connected  |+svi.  |+hai   ; aic = 92.65 
        summary(glm3) # nonmov: -connected* |+svi***|+hai** ; aic = 275.04    
    
    # v.inf ~ harmonic...
        
        ## formula 
        formula <- as.formula(paste("v.inf ~ harmonic + svi + hai"))
        
        ## glm
        glm1 <- glm(formula, family = binomial(link = "logit"), data = df_model)
        glm2 <- glm(formula, family = binomial(link = "logit"), data = df_movers)
        glm3 <- glm(formula, family = binomial(link = "logit"), data = df_nonmovers)
        
        ## summarize "v.inf ~ harmonic + svi + hai"
        summary(glm1) # all:    -harmonic  |+svi***|+hai** ; aic = 424.03
        summary(glm2) # movers: +harmonic  |+svi.  |+hai   ; aic = 94.84 
        summary(glm3) # nonmov: -harmonic* |+svi***|+hai** ; aic = 279.2    

# GLM: serology.child ~   
    
    # Get df for model at individual level 
    df_model <- df_model_child
    
    # Get dfs for movers vs. non-movers
    df_movers <- df_model[df_model$mover=="Mover" & !is.na(df_model$mover),]
    df_nonmovers <- df_model[df_model$mover=="Non-mover" & !is.na(df_model$mover),]

    # FULL MODEL from Fernandez et al. 2019
    
        ## glm
        formula <- as.formula(paste("serology.child ~ age + iva + svi + iva*svi + hai + mom.serology + sex + sero.co + sprayed + ethnicity"))
        glm1 <- glm(formula, family = binomial(link = "logit"), data = df_model)
        
        ## summarize 
        summary(glm1) 
        
        ## glmm
        formula <- as.formula(paste("serology.child ~ age + iva + svi + iva*svi + hai + mom.serology + sex + sero.co + sprayed + ethnicity + (1 | hhx.id)"))
        glmm1 <- glmer(formula, family = binomial(link = "logit"), data = df_model)
        
        ## summarize (random effect)
        summary(glmm1) # +age** |+iva*  |+svi   |-hai** |+mom.serology*  |+sex |+sero.co |-sprayed |+ethnicity |-iva*svi. ; aic = 273.3
        
    # ~ degree + age + iva + svi + hai + mom.serology
        
        ## glm
        formula <- as.formula(paste("serology.child ~ degree + age + iva + svi + hai + mom.serology"))
        glm1 <- glm(formula, family = binomial(link = "logit"), data = df_model)
        glm2 <- glm(formula, family = binomial(link = "logit"), data = df_movers)
        glm3 <- glm(formula, family = binomial(link = "logit"), data = df_nonmovers)
        
        ## summarize
        summary(glm1) # all:    +degree. |+age***|+iva** |+svi   |-hai*  |+mom.serology***; aic = 279.02 !!!adding mobility as a predictor lowered aic!!!
        summary(glm2) # movers: +degree**|+age** |+iva   |+svi*  |+hai   |+mom.serology   ; aic = 62.357
        summary(glm3) # nonmov: -degree  |+age*  |+iva** |-svi   |-hai*  |+mom.serology***; aic = 207.46

        ## glmm
        formula <- as.formula(paste("serology.child ~ degree + age + iva + svi + hai + mom.serology + (1 | hhx.id)"))
        glmm1 <- glmer(formula, family = binomial(link = "logit"), data = df_model)
        glmm2 <- glmer(formula, family = binomial(link = "logit"), data = df_movers)
        glmm3 <- glmer(formula, family = binomial(link = "logit"), data = df_nonmovers)
        
        ## summarize (random effect)
        summary(glmm1) # all:    +degree. |+age** |+iva** |+svi   |-hai.  |+mom.serology***; aic = 275.0
        summary(glmm2) # movers: +degree**|+age** |+iva   |+svi*  |+hai   |+mom.serology   ; aic = 64.4
        summary(glmm3) # nonmov: -degree  |+age*  |+iva*  |-svi   |-hai   |+mom.serology** ; aic = 200.5
        
    # ~ connected + age + iva + svi + hai + mom.serology

        ## glm
        formula <- as.formula(paste("serology.child ~ connected + age + iva + svi + hai + mom.serology"))
        glm1 <- glm(formula, family = binomial(link = "logit"), data = df_model)
        glm2 <- glm(formula, family = binomial(link = "logit"), data = df_movers)
        glm3 <- glm(formula, family = binomial(link = "logit"), data = df_nonmovers)
        
        ## summarize
        summary(glm1) # all:    -connected  |+age** |+iva** |+svi   |-hai.  |+mom.serology***; aic = 281.58
        summary(glm2) # movers: +connected  |+age*  |+iva   |+svi*  |+hai   |+mom.serology   ; aic = 74.515
        summary(glm3) # nonmov: -connected  |+age*  |+iva** |-svi   |-hai*  |+mom.serology***; aic = 205.7
        
        ## glmm
        formula <- as.formula(paste("serology.child ~ connected + age + iva + svi + hai + mom.serology + (1 | hhx.id)"))
        glmm1 <- glmer(formula, family = binomial(link = "logit"), data = df_model)
        glmm2 <- glmer(formula, family = binomial(link = "logit"), data = df_movers)
        glmm3 <- glmer(formula, family = binomial(link = "logit"), data = df_nonmovers)
        
        ## summarize (random effect)
        summary(glmm1) # all:    -connected   |+age** |+iva** |+svi   |-hai.  |+mom.serology** ; aic = 276.4
        summary(glmm2) # movers: +connected***|+age***|+iva   |+svi***|+hai***|+mom.serology   ; aic = 68.4
        summary(glmm3) # nonmov: -connected   |+age*  |+iva*  |-svi   |-hai.  |+mom.serology** ; aic = 199.0
        
    # ~ harmonic + age + iva + svi + hai + mom.serology
        
        ## glm
        formula <- as.formula(paste("serology.child ~ harmonic + age + iva + svi + hai + mom.serology"))
        glm1 <- glm(formula, family = binomial(link = "logit"), data = df_model)
        glm2 <- glm(formula, family = binomial(link = "logit"), data = df_movers)
        glm3 <- glm(formula, family = binomial(link = "logit"), data = df_nonmovers)
        
        ## summarize
        summary(glm1) # all:    +harmonic  |+age** |+iva** |+svi   |-hai.  |+mom.serology***; aic = 281.79
        summary(glm2) # movers: +harmonic* |+age** |+iva   |+svi*  |+hai   |+mom.serology   ; aic = 69.791
        summary(glm3) # nonmov: -harmonic  |+age*  |+iva***|-svi   |-hai.  |+mom.serology** ; aic = 205.53

        ## glmm
        formula <- as.formula(paste("serology.child ~ harmonic + age + iva + svi + hai + mom.serology + (1 | hhx.id)"))
        glmm1 <- glmer(formula, family = binomial(link = "logit"), data = df_model)
        glmm2 <- glmer(formula, family = binomial(link = "logit"), data = df_movers)
        glmm3 <- glmer(formula, family = binomial(link = "logit"), data = df_nonmovers)
        
        ## summarize (random effect)
        summary(glmm1) # all:    +harmonic  |+age** |+iva** |+svi   |-hai   |+mom.serology** ; aic = 276.6
        summary(glmm2) # movers: +harmonic  |+age.  |+iva   |+svi   |+hai   |+mom.serology   ; aic = 67.3
        summary(glmm3) # nonmov: -harmonic  |+age*  |+iva*  |+svi   |-hai   |+mom.serology** ; aic = 199.6
        
        
################################################################################
### PART 3A: BAYESIAN SPATIAL ANALYSIS W/ INLA
################################################################################

# Prep -------------------------------------------------------
        
    ## get df for model at individual level 
    df <- df_model_child
        
    ## inspect data (is response variable numeric binary)
    str(df)
    
    ## subset variables for modeling
    df <- subset(df, select = c(hhx.id, serology.child, degree, connected, harmonic, age, iva, svi, hai, mom.serology, mover))
    
    ## merge spatial data
    coords <- subset(df_kin, select = c(hhx.id, x, y, long, lat))
    df <- merge(df, coords, by.x = "hhx.id", all.x =TRUE)
    
    # keep only complete cases
    df <- df[complete.cases(df),]
    
    # transform binary predictors to numeric
    df$connected <- as.numeric(df$connected)
    df$connected <- ifelse(df$connected == 2, 1, 0)
    df$mom.serology <- as.numeric(df$mom.serology)
    df$mom.serology <- ifelse(df$mom.serology == 2, 1, 0)
    df$mover <- as.numeric(df$mover)
    df$mover <- ifelse(df$mover == 2, 1, 0)

# Mesh Construction -------------------------------------------------------
    
    # load library
    library(INLA)
    
    # Get coordinates of houses
    coords = cbind(df$long, df$lat) # location coords used as initial mesh vertices
    
    # Get boundary based on coordinates of our observations/houses
    boundary <- inla.nonconvex.hull(coords, convex = -0.3) # convex parameter controls how much concavity is allowed (negative val -> tighter fit, pos val -> looser fit)
    lines(boundary, add = FALSE)
    points(coords, col = "red")
    
    # Mesh construction w/ boundary
    mesh <- inla.mesh.2d(boundary = boundary, max.edge = c(0.0152, 0.03), # max.edge: values denote max. allowed triangle edge lengths in region & buffer zone, respectively
                         cutoff = 0.001) # min. allowed distance b/w points; used to avoid building too many small triangles around clustered data locations
    ## set inner edge to 1/5-1/10 of the smallest spatial pattern you want to detect:
    ## occ of seropositive person aggregated at a scale b/w 2-6km; 
    ## hotspots of infected vector at 0.2-1.8km; 
    ## flight range of triatomine varies from 0.1-1.5km
    ## given the above, let's set max.edge for inner edge at 0.01-0.02
    
    # Plot mesh and coordinate points
    plot(mesh, main = "")
    points(coords, col = "red")
    
    
# Spatial Construction -----------------------------------------------
    
# Build the A matrix, which will map the Gaussian Markov Random Field (GMRF) from the mesh nodes to the n observation location
A<-inla.spde.make.A(mesh=mesh,loc=as.matrix(coords));dim(A)
    ## you end up with a matrix on number of observations by number of nodes

# Create SPDE for spatial structure # this function will be equivalent to using a GMRF in a spatial model
spde <- inla.spde2.pcmatern(mesh, 
                            alpha=2, #alpha = nu + delta/2 = 1 + 2/2 = 2
                            ## smoothness parameter controls how quickly the spatial effect changes over distance (e.g., neighboring locations are likely to share similar random effects)
                            prior.range = c(1.5, 0.95),
                            # prior.range = c(0.5, 0.5), #assume distribution is skewed to the left (e.g., average flight radius is 500m)
                            ## specifies spatial range, which defines how far apart two points can be while still maintaining a meaningful correlation.
                            ## (range0, Prange) whereby P(range < range0) = Prange
                            ## let's set this to the maximum flight range of a triatomine, 1.5km, where we set an a priori probability of 95% that the range will be less than 1.5km.
                            ## source: T. infestans may easily fly >550 m (SchoÞeld et al. 1992) or reach 1,500 m (Schweigmann et al. 1988) in an open field
                            prior.sigma = c(.5, .5)
                            #P(sigma > 1.0) = 0.5
                            ## specifies standard deviation (overall variability) of the field: sigma0 and Psigma should reflect how much variation you're willing to allow 
                            ## if the data suggests that high variability in the spatial field is unlikely, you would want a smaller sigma0 and a lower Psigma
                            ## large sigma: If sigma is large, the spatial process is highly variable, meaning there's a lot of fluctuation or noise in the values at different locations. 
                            ## small sigma: If sigma is small, the spatial process is less variable, and neighboring locations are more similar to each other.
                            ## let's set sigma to a larger value b/c flight range of triatomine can vary widely (200-1500m) and spatial aggregation of cases at a higher spatial scale might occur b/c of dispersal of infected vectors
)

# Create all the required indexes for the SPDE model and name the effect "spatial.field" (for use in creating our formula)
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)

# Data Stacking -----------------------------------------------------------
## create an INLA stack with the input data as a list of the response, 
## the A matrix (with a 1 added for the list of covariates), and the effects (e.g., spatial effects and fixed effects)

stk <- inla.stack(data=list(y=df$serology.child), #the response
                  A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                  #these are your covariates
                  effects=list(c(list(Intercept=1), #the intercept
                                 iset),  #the spatial index
                               #the covariates
                               list(degree = df$degree,
                                    connected = df$connected,
                                    harmonic = df$harmonic,
                                    age = df$age,
                                    iva = df$iva,
                                    svi = df$svi,
                                    hai = df$hai,
                                    mom.serology = df$mom.serology,
                                    mover = df$mover)
                  ), 
                  #this is a quick name so you can call upon easily
                  tag='dat')

# Model ----------------------------------------------------------

    ## degree|connected|harmonic ONLY
        
        # create formula w/ spatial field (-1 excludes default intercept, which we code ourselves)
        formula <- y ~ -1 + Intercept + degree + f(spatial.field, model=spde)
        formula <- y ~ -1 + Intercept + connected + f(spatial.field, model=spde)
        formula <- y ~ -1 + Intercept + harmonic + f(spatial.field, model=spde)

        # fit the model (the INLA function)
        model <-inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
                     # Ntrials = n,
                     control.family = list(link = "logit"),
                     control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                     control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                     verbose = FALSE) #can include verbose=TRUE to see the log
        
        # model checks
        summary(model) #ns for all

    ## degree
        
        # create formula w/ spatial field (-1 excludes default intercept, which we code ourselves)
        formula <- y ~ -1 + Intercept + age + iva + svi + hai + mom.serology + degree + f(spatial.field, model=spde) 
        formula <- y ~ -1 + Intercept + age + iva + svi + hai + mom.serology + degree + mover + f(spatial.field, model=spde) 
        formula <- y ~ -1 + Intercept + age + iva + svi + hai + mom.serology + degree + mover + degree*mover + f(spatial.field, model=spde) 
        
        # fit the model (the INLA function)
        model <-inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
                     # Ntrials = n,
                     control.family = list(link = "logit"),
                     control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                     control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                     verbose = FALSE) #can include verbose=TRUE to see the log
        
        # model checks
        summary(model)

    ## connected
                
        # create formula w/ spatial field (-1 excludes default intercept, which we code ourselves)
        formula <- y ~ -1 + Intercept + age + iva + svi + hai + mom.serology + connected + f(spatial.field, model=spde)
        formula <- y ~ -1 + Intercept + age + iva + svi + hai + mom.serology + connected + mover + connected*mover + f(spatial.field, model=spde)
        
        # fit the model (the INLA function)
        model <-inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
                     # Ntrials = n,
                     control.family = list(link = "logit"),
                     control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                     control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                     verbose = FALSE) #can include verbose=TRUE to see the log
        
        # model checks
        summary(model) #ns: svi & hai

    # harmonic

        # create formula w/ spatial field (-1 excludes default intercept, which we code ourselves)
        formula <- y ~ -1 + Intercept + age + iva + svi + hai + mom.serology + harmonic + f(spatial.field, model=spde) 
        formula <- y ~ -1 + Intercept + age + iva + svi + hai + mom.serology + harmonic + mover + harmonic*mover + f(spatial.field, model=spde) 
        
        # fit the model (the INLA function)
        model <- inla(formula, data=inla.stack.data(stk,spde=spde),family= 'binomial', 
                     # Ntrials = n,
                     control.family = list(link = "logit"),
                     control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                     control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                     verbose = FALSE) #can include verbose=TRUE to see the log
        
        # model checks
        summary(model)

### SUMMARY ###
# Model stratifying INLA by movers and non-movers failed to converge, so instead, we investigated the interaction between degree*mover
# degree: mean(CI) = -0.140 (-0.636, 0.358) --> Effect of kin connectivity on seropositivity among children of non-movers was not statistically significant
# mover: mean(CI) = -1.626 (-2.958, -0.294) --> Movers have much lower odds of seropositivity when degree = 0
# degree*mover: mean(CI) = +1.016 (0.349, 1.683) --> Among movers, the effect of degree is strongly positive and significant
# When mover = 0: Among children in stable households, being more connected was weakly protective or neutral in terms of seropositivity
# When mover = 1, the effect of degree = -0.140 + 1.016 = +0.876: Among children of mover households, kin connectivity reverses in effect and becomes risk factor such that higher degree is associated with significantly higher odds of seropositivity.
# Household mobility moderates the relationship between social connectivity and disease risk
        ## while kin networks may provide stability and protection in stationary households,
        ## in frequently moving households, those same kin ties may increase the child's exposure to environments where vectors are present, amplifying seropositivity risk among children

