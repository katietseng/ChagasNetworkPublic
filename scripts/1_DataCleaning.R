################################################################################
##  Project: Chagas Network
##  Title: 1_DataCleaning
##  Author: Katie Tseng (katie.tseng@wsu.edu)
################################################################################

# NOTE: The original data contained identifiable information (names and exact household locations) and were de-identified prior to upload; only reproducible, non-identifiable datasets are shared here.

################################################################################
### SETUP ###
################################################################################

# Clean environment
rm(list=ls()) 
graphics.off()

# Set directory
setwd("/Users/~/Desktop/ChagasNetworkPublic")

# Load data
load('data/SourceData_public.RData')

################################################################################
### CREATE ADDITIONAL VARIABLES FOR INDIVIDUALS ###
################################################################################

# Infected vector abundance

    ## binary
    df_dem$iva.bin <- df_dem$iva
    df_dem$iva.bin <- ifelse(df_dem$iva.bin == 0, 0, 1)
    
    ## discrete for only those infested (where iva.bin == 1)
    df_dem$iva.bin1 <- df_dem$iva
    df_dem$iva.bin1 <- ifelse(df_dem$iva.bin1==0, NA, df_dem$iva.bin1)
    
    ## categorical (transform to factor)
    df_dem$iva.cat <- df_dem$iva
    df_dem$iva.cat <- ifelse(is.na(df_dem$iva.cat), NA, 
                             ifelse(df_dem$iva.cat == 0, NA,
                                    ifelse(df_dem$iva.cat >= 1 & df_dem$iva.cat <= 4, 1,
                                           ifelse(df_dem$iva.cat >= 5 & !is.na(df_dem$iva.cat), 2, df_dem$iva.cat)
                                    )
                             )
    )
    table(df_dem$iva.cat, useNA = "ifany")
    df_dem$iva.cat <- factor(df_dem$iva.cat,
                             levels = c(1, 2),
                             labels = c("1-4 low", "5+ high")
    )

# Mover (gen var from mobility)    
    
    ## binary (transform to factor)
    df_dem$mover <- NA
    df_dem$mover <- ifelse(df_dem$mobility == "Non-mover", 0, df_dem$mover)
    df_dem$mover <- ifelse(df_dem$mobility == "Mover" | df_dem$mobility == "Out-migrant" | df_dem$mobility == "In-migrant", 1, df_dem$mover)
    df_dem$mover <- factor(df_dem$mover,
                           levels = c(0, 1),
                           labels = c("Non-mover", "Mover"))

# Seropositive coinhabitants
    
    ## get seropositive data at household level
    temp <- df_dem[, c("hhx.id","serology")]
    
    ## transform serology to numeric
    temp$serology <- as.numeric(temp$serology)
    temp$serology <- ifelse(temp$serology == 1, 0, temp$serology)
    temp$serology <- ifelse(temp$serology == 2, 1, temp$serology)
    
    ## collapse
    temp <- aggregate(serology ~ hhx.id, data = temp, FUN = sum)
    
    ## subtract one
    names(temp)[names(temp) == "serology"] <- "sero.co"
    temp$sero.co <- ifelse(temp$sero.co==0, 0, temp$sero.co - 1)
    
    ## merge
    df_dem <- merge(x = df_dem, y = temp, by.x = "hhx.id", all.x = TRUE)
    rm(temp)
    
    ## create as binary var
    table(df_dem$sero.co, useNA = "ifany")
    df_dem$sero.co.bin <- df_dem$sero.co
    df_dem$sero.co.bin <- ifelse(is.na(df_dem$sero.co.bin), NA,
                                 ifelse(df_dem$sero.co.bin == 0, 0, 1))
    
    ## create discrete w/ NAs for 0s
    df_dem$sero.co.na <- df_dem$sero.co
    df_dem$sero.co.na <- ifelse(df_dem$sero.co.na == 0, NA, df_dem$sero.co.na)
    table(df_dem$sero.co.na, useNA = "ifany")

################################################################################
### CREATE ADDITIONAL VARIABLES FOR KINSHIP ###
################################################################################

### NETWORK CENTRALITY MEASURES ###  
    
    ## degree: the number of edges adjacent/linked to a node
        # A node w/ a high degree centrality is a highly connected one w/ the potential of being influential in a network
    ## betweenness: how often a given node falls along the shortest path b/w two other nodes
        # A node w/ a high betweenness has a large potential for controlling flows thr/ network and is in a position to threaten the network with flow disruptions
        # When network is non-connected, it identifies "bridges" or "brokers" within those smaller connected groups
    ## harmonic (alt to closeness centrality): mean inverse distance to all other vertices
        # Normalized: sum of inverse path lengths to other vertices divided by the number of vertices minus one
        # The inverse distance to an unreachable vertex is considered to be zero, meaning it works for disconnected networks by not returning NA or infinite
        # Identifies nodes that are important for 'information flow' or 'accessibility' even if some paths are unreachable (in a fragmented structure)
        # Nodes that are reachable from many others and with short paths will have high harmonic centrality.
        # Nodes in isolated components or on the fringe of a network will have lower scores, because they are "distant" from most of the network—even if they are central in their tiny component.

    ## build graph from matrix P (how connected are these terms; primary and secondary relations)
    library(igraph)
    g <- graph.adjacency(matrix_kin$P, weighted=T, mode='undirected')
    
    ## summarize
    summary(g)
    
    ## create empty df (here we replace, previous dataframe df_kin)
    df_kin <- as.data.frame(hhx.id)
    colnames(df_kin) <- c("hhx.id")
    
    ## get degree centrality/node degree: number of houses connected to the focal house
    df_kin$degree <- igraph::degree(g)                                              # degree centrality

    ## make degree categorical
    table(df_kin$degree, useNA= "ifany")
    df_kin$degree.cat <- df_kin$degree
    df_kin$degree.cat <- ifelse(df_kin$degree.cat >= 1 & df_kin$degree.cat <= 2, 1,
                                ifelse(df_kin$degree.cat >= 3 & df_kin$degree.cat <= 9, 2, df_kin$degree.cat))
    df_kin$degree.cat <- factor(df_kin$degree.cat,
                                levels = c(0, 1, 2),
                                labels = c("0", "1-2", "3+")
    )

    ## get betweeness
    df_kin$betweenness <- igraph::betweenness(g)                                              # degree centrality
    df_kin$betweenness.norm <- igraph::betweenness(g, normalized=TRUE)                         # normalize
    
    ## get harmonic centrality
    df_kin$harmonic <- igraph::harmonic_centrality(g)
    df_kin$harmonic.norm <- igraph::harmonic_centrality(g, normalized = TRUE) # normalized b/w 0-1
    
    ## get component
    components <- igraph::components(g)
    
    ## get membership id
    membership <- as.data.frame(components$membership)
    membership <- cbind(membership, rownames(membership))
    
    ## rename columns
    colnames(membership) <- c("component.id", "hhx.id")
    
    ## summarize
    table(membership$component.id)
    
    ## plot
    V(g)$name <- components$membership
    
    ## merge with df_kin
    df_kin <- merge(df_kin, membership, by.x = "hhx.id", all.x = TRUE)
    
    ## get whether the focal house is or isn't connected to any other house in the network
    df_kin$connected <- ifelse(df_kin$degree==0,"unconnected", "connected")
    df_kin$connected <- factor(df_kin$connected, 
                               levels = c("unconnected", "connected"),
                               labels = c("Unconnected", "Connected"))

      
### OTHER KINSHIP VARIABLES ###  

    ## kin.inf: Mean occurrence (proportion) of infestation among connected houses
    ## kin.iva.med: Median IVA among connected houses w/ infected vector; *relative to infestation, IVA is more likely to reflect focal conditions (e.g., HAI)
    ## kin.iva.sum: Cumulative IVA of connected houses w/ any infected vector
    ## kin.serology.prop: Mean occurrence (proportion) of any seropositive person among connected houses
    ## kin.serology.med: Median number of seropositive inhabitants among connected houses w/ at least one seropositive person; *speaks to past exposure
    ## kin.serology.sum: Cumulative number of seropositive inhabitants among connected houses w/ at least one seropositive person
    ## kin.svi: Median social vulnerability of connected houses
    ## kin.svi.rel: Median social vulnerability of connected houses relative to the focal house's vulnerability
    ## kin.svi.rel.sum: Cumulative social vulnerability of connected houses relative to the focal house's vulnerability
    ## kin.network.id
    
# Get adjacency list of kin connections for creating additional variables describing characteristics of the houses connected to each focal house
    
    ## For the next set of variables, we will add them one by one to our adjacency list, before aggregating the data at the house-level and taking summary measures
    
    ## convert matrix K into an adjacency list
    mat <- matrix_kin$P
    diag(mat) <- NA
    df_adj <- na.omit(data.frame(as.table(mat)))
    
    ## rename colnames
    head(df_adj)
    colnames(df_adj) <- c("hhx.id", "hhx.id.kin", "weight")
    
    ## drop rows with no kin relation
    df_adj <- df_adj[df_adj$weight!=0,]
    
    ## check: sum of weights in df_adj should be equal to sum of matrix K
    sum(df_adj$weight)
    sum(matrix_kin$P)

# kin.inf.prop: Mean occurrence (proportion) of infestation among connected houses (NOT specific to infected vector)
    
    ## get infestation values from df_base
    temp <- subset(df_base, select = c(hhx.id, v.inf))
    
    ## rename hhx.id to hhx.id.kin b/c these are the values we will use for the houses connected to each focal house
    colnames(temp) <- c("hhx.id.kin", "kin.inf")
    
    ## merge on hhx.id.kin
    df_adj <- merge(df_adj, temp, by.x = "hhx.id.kin", all.x = TRUE)
    rm(temp)
    
    ## transform to binary
    df_adj$kin.inf.prop <- NA
    df_adj$kin.inf.prop <- ifelse(df_adj$kin.inf == "infested", 1, 0)

# kin.iva.sum: Cumulative IVA among connected houses
    
    ## get iva values from df_dem
    df_iva <- aggregate(iva ~ hhx.id, data = df_dem, FUN = median)
    
    ## rename hhx.id to hhx.id.kin b/c these are the values we will use for the houses connected to each focal house
    colnames(df_iva) <- c("hhx.id.kin", "kin.iva")
    
    ## merge with adjacency list of houses w/ kin connection data
    df_iva <- merge(df_adj, df_iva, by.x = "hhx.id.kin", all.x = TRUE)
    
    ## check also iva values from PLOS supplemental file
    temp <- readxl::read_excel("data/pntd.0007430.s006.xlsx", sheet = "household-level data")
    temp_iva <- subset(temp, select = c(`House ID`,`Domiciliary T. cruzi- infected T. infestans abundance`))
    colnames(temp_iva) <- c("hhx.id.kin", "kin.iva.pntd")
    
    ## merge
    df_iva <- merge(df_iva, temp_iva, by.x = "hhx.id.kin", all.x = TRUE)
    rm(temp)
    
    ## replace where kin.iva == NA with value from kin.iva.pntd, if available
    df_iva$kin.iva <- ifelse(is.na(df_iva$kin.iva) & !is.na(df_iva$kin.iva.pntd), df_iva$kin.iva.pntd, df_iva$kin.iva)
    df_iva$kin.iva.pntd <- NULL
    
    ## create cumulative IVA variable
    df_iva_sum <- aggregate(kin.iva ~ hhx.id, data = df_iva, FUN = sum, na.rm = TRUE)
    colnames(df_iva_sum) <- c("hhx.id", "kin.iva.sum")
    
    ## merge cumulative IVA values and missing IVA values with adjacency list by hhx.id
    df_adj <- merge(df_adj, df_iva_sum, by.x = "hhx.id", all.x = TRUE)

# kin.iva.med: Median IVA among connected houses w/ infected vector; *relative to infestation, IVA is more likely to reflect focal conditions (e.g., HAI)

    ## values: NA if none of the connected houses have IVA data; 0 if not connected to any house or if all connected houses have an IVA of 0; >0 if median of those houses w/ IVA is greater than 0
    
    ## create median IVA variable, replacing values of 0 with NA
    df_iva$kin.iva.na <- df_iva$kin.iva
    df_iva$kin.iva.na <- ifelse(df_iva$kin.iva.na==0, NA, df_iva$kin.iva.na)
    
    ## collapse kin.iva.med by hhx.id, taking the median of only those houses with infected vector abundance (>0)
    df_iva_med <- aggregate(kin.iva.na ~ hhx.id, data = df_iva, FUN = median, na.rm = TRUE)
    colnames(df_iva_med) <- c("hhx.id", "kin.iva.med")
    
    ## merge median IVA values with adjacency list by hhx.id
    df_adj <- merge(df_adj, df_iva_med, by.x = "hhx.id", all.x = TRUE)
    
    ## replace NA values that were originally 0's (i.e., all connected houses had an IVA of 0)
    df_adj$kin.iva.med <- ifelse(df_adj$kin.iva.sum == 0, 0, df_adj$kin.iva.med)
    
    ## clean environment
    rm(df_iva, df_iva_med, df_iva_sum)

# kin.sero.prop: Mean occurrence (proportion) of any seropositive person among connected houses
    
    ## get serology values from df_dem
    df_sero <- aggregate(serology ~ hhx.id, data = df_dem, FUN = sum)
    
    ## rename hhx.id to hhx.id.kin b/c these are the values we will use for the houses connected to each focal house
    colnames(df_sero) <- c("hhx.id.kin", "kin.sero")
    
    ## merge on hhx.id.kin
    df_adj <- merge(df_adj, df_sero, by.x = "hhx.id.kin", all.x = TRUE)
    
    ## check for anomalies (e.g., NAs)
    unique(df_adj$kin.sero)
    
    ## add occurrence of seropositive person
    df_adj$kin.sero.prop <- df_adj$kin.sero
    df_adj$kin.sero.prop <- ifelse(df_adj$kin.sero.prop > 0, 1, df_adj$kin.sero.prop)

# kin.sero.sum: Cumulative number of seropositive inhabitants among connected houses w/ at least one seropositive person
    
    ## create cumulative kin serology variable
    df_sero_sum <- aggregate(kin.sero ~ hhx.id, data = df_adj, FUN = sum, na.rm = TRUE)
    colnames(df_sero_sum) <- c("hhx.id", "kin.sero.sum")
    
    ## merge median IVA values and missing IVA values with adjacency list by hhx.id
    df_adj <- merge(df_adj, df_sero_sum, by.x = "hhx.id", all.x = TRUE)
    
    # kin.sero.med: Median number of seropositive inhabitants among connected houses w/ at least one seropositive person; *speaks to past exposure
    
    ## subset serology data from adjacency list of connected houses by hhx.id.kin
    df_sero <- subset(df_adj, select = c(hhx.id, hhx.id.kin, kin.sero))
    
    ## create median kin serology variable, replacing values of 0 with NA
    df_sero$kin.sero.na <- df_sero$kin.sero
    df_sero$kin.sero.na <- ifelse(df_sero$kin.sero.na==0, NA, df_sero$kin.sero.na)
    
    ## collapse kin.sero.na by hhx.id, taking the median of only those houses with infected vector abundance (>0)
    df_sero_med <- aggregate(kin.sero.na ~ hhx.id, data = df_sero, FUN = median, na.rm = TRUE)
    colnames(df_sero_med) <- c("hhx.id", "kin.sero.med")
    
    ## merge median serology values with adjacency list by hhx.id
    df_adj <- merge(df_adj, df_sero_med, by.x = "hhx.id", all.x = TRUE)
    
    ## replace NA values that were originally 0's (i.e., all connected houses had an IVA of 0)
    df_adj$kin.sero.med <- ifelse(df_adj$kin.sero.sum == 0, 0, df_adj$kin.sero.med)
    
    ## clean environment
    rm(df_sero, df_sero_sum, df_sero_med)

# kin.svi.med: Median social vulnerability of connected houses to each focal house
    
    ## get svi values from df_base
    df_svi <- df_base[,c("hhx.id", "svi")]
    
    ## create second df of svi values for connected houses
    df_ksvi <- df_svi
    colnames(df_ksvi) <- c("hhx.id.kin", "kin.svi.med")
    
    ## merge svi data by hhx.id
    df_adj <- merge(df_adj, df_svi, by.x = "hhx.id", all.x = TRUE)
    
    ## merge kin.svi.med data by hhx.id.kin
    df_adj <- merge(df_adj, df_ksvi, by.x = "hhx.id.kin", all.x = TRUE)
    
    ## check for NAs and anomalies 
    unique(df_adj$kin.svi.med)
    df_adj[which(is.na(df_adj$kin.svi.med)),] #PG119/PG68
    df_adj[which(is.na(df_adj$svi)),]
    
    ## check svi in demographic data
    df_dem$svi[which(df_dem$hhx.id=="PG119")] #NA
    df_dem$svi[which(df_dem$hhx.id=="PG68")]  #NA
    df_dem$svi[which(df_dem$hhx.id=="PC109")] #0
    
    ## correct svi value for PC109
    df_adj$kin.svi.med <- ifelse(df_adj$hhx.id.kin=="PC109", 0, df_adj$kin.svi.med)
    df_adj$svi <- ifelse(df_adj$hhx.id == "PC109", 0, df_adj$svi)
    
    ## if kin.svi.med is missing (NA), assign focal house's own SVI value
    df_adj$kin.svi.med <- ifelse(is.na(df_adj$kin.svi.med), df_adj$svi, df_adj$kin.svi.med)

# kin.svi.relmed: Median social vulnerability of connected houses relative to the focal house's vulnerability
    
    ## median(kin.svi – svi)
    df_adj$kin.svi.relmed <- df_adj$kin.svi.med - df_adj$svi
    
    ## check
    which(df_adj$kin.svi.med == df_adj$svi)
    which(df_adj$kin.svi.relmed == 0)

# kin.svi.relsum: Cumulative social vulnerability of connected houses relative to the focal house's vulnerability
    
    ## collapse kin.svi by hhx.id, taking the sum of the svi of connected houses
    df_svi_sum <- aggregate(kin.svi.relmed ~ hhx.id, data = df_adj, FUN = sum, na.rm = TRUE)
    colnames(df_svi_sum) <- c("hhx.id", "kin.svi.relsum")
    
    ## merge 
    df_adj <- merge(df_adj, df_svi_sum, by.x = "hhx.id", all.x = TRUE)
    
    ## clean environment
    rm(df_ksvi, df_svi_sum) #keep df_svi for now, we'll need it when we merge with df_kin data

# Collapse at house-level and combine with df_kin
    
    ## Which summary measure to take?
    # View(df_adj[which(df_adj$hhx.id == "CC40"),])
    # kin.inf.prop         -> take the mean to get the proportion infested(unconnected: NA --> 0)
    # kin.iva.sum          -> take the mean (unconnected: NA --> 0)
    # kin.iva.med       -> take the mean (unconnected: NA --> 0)
    # kin.sero.prop        -> take the mean to get the proportion with any seropositive person (unconnected: NA --> 0)
    # kin.sero.sum         -> take the mean (unconnected: NA --> 0)
    # kin.sero.med      -> take the mean (unconnected: NA --> 0)
    # kin.svi.med       -> take the median (unconnected: NA --> svi of focal house)
    # kin.svi.relmed       -> take the median (unconnected: NA --> 0)
    # kin.svi.relsum       -> take the mean (unconnected: NA --> 0)
    
    ## subset and collapse, taking the mean: 
    df_adj_mean <- subset(df_adj, select = c(hhx.id, kin.inf.prop, kin.iva.sum, kin.iva.med, kin.sero.prop, kin.sero.sum, kin.sero.med, kin.svi.relsum))
    df_adj_mean <- aggregate(. ~ hhx.id, data = df_adj_mean, FUN = mean, na.rm=TRUE, na.action=NULL) # na.action = NULL preserves observations containing NA
    
    ## subset and collapse, taking the median
    df_adj_med <- subset(df_adj, select = c(hhx.id, kin.svi.med, kin.svi.relmed))
    df_adj_med <- aggregate(. ~ hhx.id, data = df_adj_med, FUN = median, na.rm=TRUE, na.action=NULL)
    
    ## combine with df_kin
    df_kin <- merge(df_kin, df_adj_mean, by.x = "hhx.id", all.x = TRUE)
    df_kin <- merge(df_kin, df_adj_med, by.x = "hhx.id", all.x = TRUE)
    
    ## clean environment
    rm(df_adj, df_adj_mean, df_adj_med)
    
    ## replace NA of unconnected houses
    df_kin$kin.inf.prop <- ifelse(df_kin$connected == "Unconnected" & is.na(df_kin$kin.inf.prop), 
                                  0, 
                                  df_kin$kin.inf.prop)
    df_kin$kin.iva.sum <- ifelse(df_kin$connected == "Unconnected" & is.na(df_kin$kin.iva.sum), 0, df_kin$kin.iva.sum)
    df_kin$kin.iva.med <- ifelse(df_kin$connected == "Unconnected" & is.na(df_kin$kin.iva.med), 0, df_kin$kin.iva.med)
    df_kin$kin.sero.prop <- ifelse(df_kin$connected == "Unconnected" & is.na(df_kin$kin.sero.prop), 0, df_kin$kin.sero.prop)
    df_kin$kin.sero.sum <- ifelse(df_kin$connected == "Unconnected" & is.na(df_kin$kin.sero.sum), 0, df_kin$kin.sero.sum)
    df_kin$kin.sero.med <- ifelse(df_kin$connected == "Unconnected" & is.na(df_kin$kin.sero.med), 0, df_kin$kin.sero.med)
    df_kin$kin.svi.relmed <- ifelse(df_kin$connected == "Unconnected" & is.na(df_kin$kin.svi.relmed), 0, df_kin$kin.svi.relmed)
    df_kin$kin.svi.relsum <- ifelse(df_kin$connected == "Unconnected" & is.na(df_kin$kin.svi.relsum), 0, df_kin$kin.svi.relsum)

    ## get svi of unconnected houses
    df_kin <- merge(df_kin, df_svi, by.x = "hhx.id", all.x = TRUE)
    
    ## replace NA of kin.svi.med with svi of focal house
    df_kin$kin.svi.med <- ifelse(is.na(df_kin$kin.svi.med), df_kin$svi, df_kin$kin.svi.med) #if NA, return the svi of focal house
    df_kin$svi = NULL
    
    ## replace NaN with NA
    df_kin$kin.sero.prop <- ifelse(is.nan(df_kin$kin.sero.prop), NA, df_kin$kin.sero.prop)
    df_kin$kin.sero.sum <- ifelse(is.nan(df_kin$kin.sero.sum), NA, df_kin$kin.sero.sum)
    df_kin$kin.sero.med <- ifelse(is.nan(df_kin$kin.sero.med), NA, df_kin$kin.sero.med)

################################################################################
### MERGE ADDITIONAL HOUSEHOLD-LEVEL VARIABLES ###
################################################################################

# merge df_base
df_kin <- merge(df_kin, df_base, by.x = "hhx.id", all.x = TRUE)

names(df_kin)

# Explore whether serology data is missing or partially missing at the household level
    
    ## subset serology data
    df_sero <- df_dem[, c("hhx.id","serology")]
    
    ## identify which hhx have any missing/NA data (including those with partial data)
    df_sero_NA <- df_sero[is.na(df_sero$serology),]
    df_sero_NA$serology.missing <- 1
    df_sero_NA$serology = NULL
    
    ## collapse
    df_sero_NA <- aggregate(serology.missing ~ hhx.id, data = df_sero_NA, FUN = mean)
    df_sero <- aggregate(serology ~ hhx.id, data = df_sero, FUN = sum)
    
    ## merge
    df_kin <- merge(x = df_kin, y = df_sero, by.x = "hhx.id", all.x = TRUE)
    df_kin <- merge(x = df_kin, y = df_sero_NA, by.x = "hhx.id", all.x = TRUE)
    
    ## transform serology.missing to binary
    df_kin$serology.missing <- ifelse(is.na(df_kin$serology.missing), 0, df_kin$serology.missing)
    
    # ## drop rows with no serology data
    # df_kin <- df_kin[which(!is.na(df_kin$serology)),]

# Occurrence of seropositive person (factor): No occurrence (including partial data), occurrence, and no information
    
    # create var
    df_kin$serology.occ <- NA
    
    ## 0 if no occurrence (including partial data), 1 if occurrence
    df_kin$serology.occ <- ifelse(df_kin$serology >=1, 1, 0)
    
    ## to factor
    df_kin$serology.occ <- factor(df_kin$serology.occ,
                                  levels = c(0, 1),
                                  labels = c("No occurrence", "Occurrence"))
    
# Categorical number of seropositive persons (factor): 0, 1-2, 3+, missing
    
    ## create var
    df_kin$serology.cat <- df_kin$serology
    
    ## create cats
    df_kin$serology.cat <- ifelse(df_kin$serology.cat == 1 | df_kin$serology.cat == 2, 1, 
                                  ifelse(df_kin$serology.cat >= 3, 2, df_kin$serology.cat))
    
    ## to factor
    df_kin$serology.cat <- factor(df_kin$serology.cat,
                                  levels = c(0, 1, 2),
                                  labels = c("0", "1-2", "3+"))
    
# Explore whether children's serology data is missing or partially missing at the household level
    
    ## subset serology data for children <15yo
    df_sero <- df_dem[df_dem$age<15 & !is.na(df_dem$age), c("hhx.id","age", "serology.child")]
    
    ## identify which hhx have missing/NA data (including those with partial data)
    df_sero_NA <- df_sero[is.na(df_sero$serology.child),]
    df_sero_NA$serology.child.missing <- 1
    df_sero_NA$serology.child = NULL
    
    ## collapse
    df_sero_NA <- aggregate(serology.child.missing ~ hhx.id, data = df_sero_NA, FUN = mean)
    df_sero <- aggregate(serology.child ~ hhx.id, data = df_sero, FUN = sum)
    
    ## merge
    df_kin <- merge(x = df_kin, y = df_sero, by.x = "hhx.id", all.x = TRUE)
    df_kin <- merge(x = df_kin, y = df_sero_NA, by.x = "hhx.id", all.x = TRUE)
    
    ## transform serology.child.missing to binary
    df_kin$serology.child.missing <- ifelse(is.na(df_kin$serology.child.missing), 0 , df_kin$serology.child.missing)
    
    ## clean environment
    rm(df_sero, df_sero_NA)
    
# Occurrence of seropositive child (factor): 
    
    ## create var (0 if no child, 1 if child but no occ, 2 if occ, NA if no data on child or seropositivity)
    df_kin$serology.child.occ <- NA
    
    ## create category: 0 if no child
    df_kin$serology.child.occ <- ifelse(df_kin$hab_0a5 == 0 & df_kin$hab_6a14 == 0, 0, df_kin$serology.child.occ)
    
    ## create category: 1 if child but no occurrence
    df_kin$serology.child.occ <- ifelse(df_kin$serology.child == 0 & !is.na(df_kin$serology.child), 1, df_kin$serology.child.occ)
    
    ## create category: 2 if child exists and occurrence of seropositive child
    df_kin$serology.child.occ <- ifelse(df_kin$serology.child >= 1 & !is.na(df_kin$serology.child), 2, df_kin$serology.child.occ)
    
    ## to factor
    df_kin$serology.child.occ <- factor(df_kin$serology.child.occ,
                                        levels = c(0, 1, 2),
                                        labels = c("No child", "No occurrence", "Occurrence"))
    
# Add mean age
    
    ## subset age and get mean
    df_age <- df_dem[!is.na(df_dem$age), c("hhx.id", "age")]
    df_age <- aggregate(cbind(age) ~ hhx.id, data = df_age, FUN = mean, na.action = na.pass)
    
    ## merge
    df_kin <- merge(x = df_kin, y = df_age, by = "hhx.id", all.x = TRUE)
    
# Add infected vector abundance data at household level
    
    ## subset inf vector abundance data
    df_iva <- df_dem[!is.na(df_dem$iva), c("hhx.id", "iva", "iva.inf", "iva.bin", "iva.bin1", "iva.cat")]
    df_iva <- aggregate(cbind(iva, iva.inf, iva.bin, iva.bin1, iva.cat) ~ hhx.id, data = df_iva, FUN = mean, na.action = na.pass)

    ## merge
    df_kin <- merge(x = df_kin, y = df_iva, by.x = "hhx.id", all.x = TRUE)

# Add commmunity 
    
    df_kin$community <- NA
    df_kin$community <- ifelse(grepl("BC", df_kin$hhx.id), "BC", df_kin$community)
    df_kin$community <- ifelse(grepl("CC", df_kin$hhx.id), "CC", df_kin$community)
    df_kin$community <- ifelse(grepl("CD", df_kin$hhx.id), "CD", df_kin$community)
    df_kin$community <- ifelse(grepl("PC", df_kin$hhx.id), "PC", df_kin$community)
    df_kin$community <- ifelse(grepl("PG", df_kin$hhx.id), "PG", df_kin$community)
    df_kin$community <- ifelse(grepl("PV", df_kin$hhx.id), "PV", df_kin$community)
    df_kin$community <- ifelse(grepl("RN", df_kin$hhx.id), "RN", df_kin$community)
    df_kin$community <- as.factor(df_kin$community)
    
# Sort by house id
library(dplyr)
my_order <- gtools::mixedsort(df_kin$hhx.id)
df_kin <- df_kin %>% arrange(match(df_kin$hhx.id, my_order))
rm(my_order)

# Reorder columns
names(df_kin)
col_order <- c("hhx.id", "x", "y", "long", "lat","community", # note to public: remove "long" and "lat" if reproducing study
               "hab_tot", "hab_0a5", "hab_6a14", "hab_0a14",
               "serology", "serology.occ", "serology.cat", "serology.missing", 
               "serology.child", "serology.child.occ", "serology.child.missing", 
               "age", "ethnicity", "mobility", "svi", "hai",
               "v.inf", "iva.inf", "iva", "iva.bin", "iva.bin1", "iva.cat", 
               "component.id", "degree", "degree.cat", "betweenness", "betweenness.norm", "harmonic", "harmonic.norm", "connected",
               "kin.sero.prop", "kin.sero.med", "kin.sero.sum", 
               "kin.svi.med", "kin.svi.relmed", "kin.svi.relsum", 
               "kin.inf.prop", "kin.iva.med", "kin.iva.sum")

df_kin <- df_kin[, col_order]
names(df_kin)

################################################################################
### SAVE CLEAN DATA ###
################################################################################

# Save dfs as Rdata file for PRIVATE REPOSITORY
save(df_dem, df_kin, edgelists, matrix_kin, hhx.id, file='data/CleanData.RData')

