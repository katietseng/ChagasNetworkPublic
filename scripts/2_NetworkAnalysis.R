##############################################################
##  Project: Chagas Network
##  Title: 1_NetworkAnalysis
##  Author: Katie Tseng (katie.tseng@wsu.edu)
##############################################################

################################################################################
### PART 1: KINSHIP NETWORK STRUCTURE
################################################################################

# Setup

    ## Clean environment
    rm(list=ls()) 
    graphics.off()
    
    ## Set directory
    setwd("/Users/~/Desktop/ChagasNetworkPublic")
    
    ## Load data
    load('data/CleanData.RData')
    
# Overview of Kinship Network

    ## How many households have any type of kinship tie? 
    df_k <- edgelists$df_k[edgelists$df_k$weights>0,]
    hhx.id_k <- c(df_k$v1, df_k$v2)
    length(unique(hhx.id_k)) #n=281
    
    ## How many households have no ties?
    row_sums <- rowSums(matrix_kin$K)
    sum(row_sums == 0) #n=122
    
    ## How many parent-child kinship ties were identified? 
    sum(edgelists$df_p$weights) #n=270
    
    ## How many sibling kinship ties were identified? 
    sum(edgelists$df_g$weights) #n=566
    
    ## How many total kinship ties were identified?
    sum(edgelists$df_k$weights) #n=836
    
    ## Check
    sum(matrix_kin$K)/2 #n=836
    
    ## Some household pairs are connected by more than one kinship tie; how many unique hhx-to-hhx connections were identified?
    length(which(matrix_kin$K>0))/2 #n=435
    
# Parent-Child Network: Locality
    
    ## How many child-parent ties occur within communities vs. across communities?
  
    ## extract community prefixes (first two) from each hhx.id
    df_p <- edgelists$df_p[edgelists$df_p$weights>0,]
    df_p$comm_from <- substr(df_p$v1, 1, 2)
    df_p$comm_to <- substr(df_p$v2, 1, 2)
    
    ## compare community prefixes to determine if within or between
    df_p$within_community <- df_p$comm_from == df_p$comm_to
    
    ## count how many unique ties are within vs. between
    table(df_p$within_community) #n=150 within; n=34 across
    prop.table(table(df_p$within_community))
    
    ## count how many total ties are within vs. between
    sum(df_p[which(df_p$within_community),]$weights) #n=217
    
# Parent-Child Network: Ethnicity
    
    ## How many parent-child ties are between households of the same ethnicity?
      
    ## get ethnicity of all houses
    df_eth <- subset(df_kin, select = c("hhx.id", "ethnicity"))
    
    ## change column name for merging
    colnames(df_eth) <- c("v1", "eth_from")
    
    ## match ethnicity of 'from' house
    df_p <- merge(df_p, df_eth, by.x = "v1", all.x = TRUE)
    
    ## change column name for merging
    colnames(df_eth) <- c("v2", "eth_to")
    
    ## match ethnicity of 'to' house
    df_p <- merge(df_p, df_eth, by.x = "v2", all.x = TRUE)
    
    ## compare household ethnicity to determine if same or different
    df_p$ethnicity_pair <- with(df_p, ifelse(is.na(eth_from) | is.na(eth_to), NA,
                                            ifelse(eth_from == eth_to, paste(eth_from, eth_to, sep = "-"),
                                                  "Qom-Creole")
                                            )
                                )
    
    ## count the number of ties that are of the same or different ethnicity
    table(df_p$ethnicity_pair) #n=7 Creole-Creole; n=5 Qom-Creole; n=172 Qom-Qom
    prop.table(table(df_p$ethnicity_pair))

    ## Clean environment
    rm(df_k, df_p)
    
# Parent-Child Network: Component Structure
    
    ## get matrix P
    P <- as.matrix(unlist(matrix_kin$P))
    isSymmetric(P)
    
    ## build graph from matrix P (how connected are these terms; primary and secondary relations)
    library(igraph)
    g <- graph.adjacency(P, weighted=T, mode='undirected')
    
    ## summarize
    summary(g)
    
    # Identify connected components
    comps <- components(g)
    
    ## How many connected components are there?
    print(comps$no) #n=223
    
    ## How many of those represent houses with no ties
    length(which(comps$csize==1))
    
    ## How many connect more than one house
    length(which(comps$csize>1))

    ## What is the average and median size of a component excluding isolated houses?
    conn_comps <- comps$csize[which(comps$csize>1)]
    summary(conn_comps) #median=3
    table(comps$csize)    
    hist(comps$csize)
    
    # Find IDs of five largest components
    largest_ids <- order(comps$csize, decreasing = TRUE)[1:5]
    
    ## Extract and plot each component
        
        ## layout for 5 plots in one window 
        par(mfrow = c(2, 3))  
        
        ## plot all 5
        for (i in largest_ids) {
          subg <- induced_subgraph(g, which(comps$membership == i))
          
          ## get matching coordinates
          coords <- as.matrix(df_kin[match(V(subg)$name, df_kin$hhx.id), c("long", "lat")])
          
          plot(subg,
               main = paste("Component", i, "- Size", comps$csize[i]),
               vertex.label = V(subg)$name,
               vertex.size = 25,
               vertex.color = "lightblue",
               layout = coords)
        }
        
    ## Plot just the largest component
    
        ## Reset plotting layout to default
        par(mfrow = c(1, 1))  
        
        ## get largest component and create the subgraph
        largest_comp_id <- which.max(comps$csize)
        g_largest <- induced_subgraph(g, which(comps$membership == largest_comp_id))
        
        ## get matching coordinates
        coords <- as.matrix(df_kin[match(V(g_largest)$name, df_kin$hhx.id), c("long", "lat")])
        
        ## plot largest component only
        plot(g_largest, 
             # layout = coords,
             layout = layout_with_fr(g_largest),
             vertex.size = 4, vertex.label.cex = 0.5, vertex.color = "skyblue"
        )
        
# Parent-Child Network: Other Properties
    
    ## What is the network density? The proportion of possible ties that are realized?
    edge_density(g) #density = (2*number of edges)/(n(n-1)) = 2(184)/(403*402)
    
    ## What is the network diameter? The longest shortest path b/w any two nodes (ignoring edge weights)?
    diameter(g, weights=NA) #11

        ## get diameters of each component
        diameters <- sapply(1:comps$no, function(i) {
          subg <- induced_subgraph(g, which(comps$membership == i))
          diameter(subg)
        })
        
        ## get the index number of the component with the largest diameter
        which_component <- which.max(diameters)
        
        ## get the nodes that are in the component with the largest diameter
        nodes_in_component <- which(comps$membership == which_component)
        
        ## get the subgraph of the component with the largest diameter
        g_diam <- induced_subgraph(g, nodes_in_component)
        
        ## get the diameter path (nodes)
        diam_path <- get_diameter(g_diam)
    
            ## Color nodes and edges in the diameter path
            V(g_diam)$color <- "lightgray"
            V(g_diam)[diam_path]$color <- "tomato"
            
            E(g_diam)$color <- "gray80"
            E(g_diam)$width <- 1
            
            ## get matching coordinates
            coords <- as.matrix(df_kin[match(V(g_diam)$name, df_kin$hhx.id), c("long", "lat")])
            
            ## plot
            set.seed(4321)
            plot(g_diam,
                 vertex.label = NA, vertex.label.cex = 0.7, vertex.label.color = "black",
                 vertex.size = 6,
                 vertex.label.dist = 1.5, vertex.label.degree = 3.5, vertex.label.family = "Arial",
                 layout = layout_with_fr,
                 main = "Diameter Path in the Largest Component (n=22)")

    ## What is the average path length? 
            
        ## For the whole network, where only the lengths of the existing paths are considered and averaged
        average.path.length(g, weights = NA, directed = FALSE, unconnected = TRUE) 
            
        ## For the largest component
        average.path.length(g_largest, weights = NA)

# Parent-Child Network: Visualize
        
    ## Skip down to section on ASSORTATIVE MIXING for network plots
        
##############################################################
### PART 2: SPATIAL DISTRIBUTION OF NETWORK TIES
##############################################################

# Construct the spatial network (distance matrix) for connected houses only (houses w/ at least one parent-child tie)

    ## load library
    library(sf)
        
    ## keep only houses with parent-child relations
    df_hhx <- df_kin[which(df_kin$hhx.id %in% colnames(matrix_kin$P_lim)),]
    
    ## sort by house id
    my_order <- gtools::mixedsort(df_hhx$hhx.id)
    df_hhx <- df_hhx %>% dplyr::arrange(match(df_hhx$hhx.id, my_order))
    
    ## transform df into sf object, for spatial plotting
    sf_hhx <- st_as_sf(df_hhx, coords = c("x", "y"),
                       crs = st_crs(32721)) # UTM Zone 21S: EPSG:32721
    
    ## verify the CRS
    st_crs(sf_hhx$geometry)
    
    ## create distance matrix in kilometers
    D <- st_distance(sf_hhx, sf_hhx)
    D <- units::set_units(D, "km")
    D <- units::drop_units(D)
    
    ## change dimension names to match our matrix of houses
    row.names(D) <- sf_hhx$hhx.id
    colnames(D) <- sf_hhx$hhx.id
    
    ## plot hhx coordinates by community
    library(ggplot2)
    map_hhx_Plim <- ggplot() +
      geom_sf(data = sf_hhx, aes(color=df_hhx$community), size = 2) +
      scale_color_manual(values=rainbow(7))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Get the cumulative and incremental count of neighboring houses across distance intervals

    ## create vector of distances to explore
    max(D)
    distance <- seq(0.2, 6, 0.2) # Fernandez 2019: min cell size/distance of analysis was 200m assuming that each household had at least 3 neighbors; the maximum distance was set at 6km (i.e., half of the dimension of the area)
    
    ## create empty df for cumulative house counts where columns will represent various ranges
    spa_cum <- data.frame(matrix(nrow = ncol(D), ncol = length(distance)))
    
    ## create empty df for incremental house counts where columns will represent various ranges
    spa_inc <- data.frame(matrix(nrow = ncol(D), ncol = length(distance)))
    
    ## set counter to 0 (for referencing the nth object in distance)
    counter = 0
    
    ## for loop over distances of increasing range
    for (j in distance) {
      
      ## set counter equal to the nth iteration (i.e., the nth object in 'distance')
      counter = counter + 1
      
      ## for loop over each house to obtain the cumulative count of houses
      for(i in 1:ncol(D)) {
        
        ## get index of houses whose distance is less than `j' of focal house  
        neighbors = which(D[,i] < j, arr.ind=T)
        
        ## replace value of row `i` and column `counter` with the number of houses (whose proximity was less than `j') subtract one b/c distance b/w house and itself (diagonals in D) is 0
        spa_cum[i,counter] <- ifelse(length(as.vector(neighbors))==0, 0, length(as.vector(neighbors))-1) 
        
        ## change column name to range distance (e.g., 200m)
        colnames(spa_cum)[counter] <- paste0(j)
      }
      
      ## to obtain the incremental count (e.g., the count of neighboring houses for a specific radial band)...
      if (counter>1) {
        
        ## set incremental count equal to that of the cumulative count minus the cumulative count from the previous iteration
        spa_inc[,counter] = spa_cum[,counter] - spa_cum[,counter-1]
        
      } else {
        
        ## set incremental count equal to the cumulative count for the first iteration (where counter = 1)
        spa_inc[,counter] = spa_cum[,counter]
      }
      
      ## change column name to range distance (same as df of cumulative counts)
      colnames(spa_inc) = colnames(spa_cum)
    }  
    
    ## reshape df of cumulative counts 
    spa_cum_long <- data.table::melt(data.table::setDT(spa_cum), variable.name = "radius")
    spa_cum_long$radius <- as.factor(spa_cum_long$radius)
    
    ## reshape df of incremental counts
    spa_inc_long <- data.table::melt(data.table::setDT(spa_inc), variable.name = "radius")
    spa_inc_long$radius <- as.factor(spa_inc_long$radius)
    
# Visualize the spatial distribution of neighboring households
    
    ## plot scatterplot of cumulative count
    scatterplot_spatialhhx_cum_Plim <- ggplot(spa_cum_long, aes(as.numeric(radius), value)) +
      geom_point() +
      geom_smooth(method = "loess") + #local regression fitting (fits a polynomial surface determined by one or more numerical predictors, using local fitting by weighted least squares)
      scale_x_continuous("Distance interval (km)", breaks = seq(0,30,1), label = c("0", as.character(distance))) +
      scale_y_continuous("Number of houses") +
      ggtitle("Scatterplot of cumulative counts of neighboring houses across increasing range") 
    
    ## plot violins+boxplots of incremental count
    violinplot_spatialhhx_inc_Plim <- ggplot(spa_inc_long, aes(x=radius, y=value, fill=radius)) +
      geom_violin(width=3, color="gray50") +
      geom_boxplot(width=0.1, color="black", alpha=0.8) +
      viridis::scale_fill_viridis(discrete = TRUE, alpha=0.5) +
      theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      xlab("Distance interval (km)") +
      ylab("Number of houses") +
      ggtitle("Distribution of incremental counts of neighboring houses across bands of increasing range") 
    
# What is the median number of neighboring houses every 0.2km? Do we notice any patterns in distribution?
    
    ## median or mean?
    hist(spa_inc_long$value) # not normally distributed so use median
    medians_spa_inc <- sapply(spa_inc, median)
    
    ## get summary stats
    summary(medians_spa_inc)
    
    ## is the distribution of the number of neighboring houses non-unimodal?
    library(diptest)
    dip.test(spa_inc_long$value) #Hartigan's dip test for unimodality, p-val sig so non-unimodal/at-least bimodal
    
    ## is the distribution of the number of neighboring houses bimodal?
    library(LaplacesDemon)
    is.unimodal(spa_inc_long$value) #false
    is.bimodal(spa_inc_long$value) #true
    is.trimodal(spa_inc_long$value) #false
    
# (2b) What proportion of those (connected houses) are connected by kinship at various distances?
    
    ## create empty df where columns represent various distances
    soc_cum <- data.frame(matrix(nrow = ncol(D), ncol = length(distance)))
    soc_inc <- data.frame(matrix(nrow = ncol(D), ncol = length(distance)))
    
    ## count the cumulative and incremental number of houses that fall w/in each range for each house 
    counter = 0
    for (j in distance) {
      counter = counter + 1
      for(i in 1:ncol(D)) {
        neighbors = which(D[,i] < j & matrix_kin$P_lim[,i] >= 1, arr.ind=T)
        soc_cum[i,counter] <- length(as.vector(neighbors))  
        colnames(soc_cum)[counter] <- paste0(j)
        #subtract 1 b/c distance between house and itself is 0
      }
      
      if (counter>1) {
        soc_inc[,counter] = soc_cum[,counter] - soc_cum[,counter-1]
      } else {
        soc_inc[,counter] = soc_cum[,counter]
      }
      colnames(soc_inc) = colnames(soc_cum)
    }  
    
    ## convert dfs of cumulative counts to matrix 
    M_social <- as.matrix(soc_cum)
    M_spatial <- as.matrix(spa_cum)
    
    ## divide to get proportions
    M_prop_cum <- ifelse(M_spatial != 0, 
                         M_social/M_spatial,
                         0)
    
    ## convert dfs of incremental counts to matrix 
    M_social <- as.matrix(soc_inc)
    M_spatial <- as.matrix(spa_inc)
    
    ## divide to get proportions
    M_prop_inc <- ifelse(M_spatial != 0, 
                         M_social/M_spatial,
                         0)
    
    ## convert back to df
    prop_cum <- as.data.frame(M_prop_cum)
    prop_inc <- as.data.frame(M_prop_inc)
    
    ## load library for calculating mean and CI
    library(Rmisc)
    
    ## create empty df for storing mean and CI
    prop_cum_mean <- data.frame(matrix(nrow = ncol(prop_cum+prop_inc), ncol = 4))
    colnames(prop_cum_mean) <- c("distance", "mean", "ll", "ul")
    prop_inc_mean <- prop_cum_mean
    
    ## calculate mean and CI at various distances for cumulative values 
    for (i in 1:ncol(prop_cum)) {
      prop_cum_mean[i,"distance"] <- distance[[i]]
      prop_cum_mean[i,"mean"] <- CI(prop_cum[[i]], ci=0.95)["mean"]
      prop_cum_mean[i,"ll"] <- CI(prop_cum[[i]], ci=0.95)["lower"]
      prop_cum_mean[i,"ul"] <- CI(prop_cum[[i]], ci=0.95)["upper"]
    }
    
    ## calculate mean and CI at various distances for incremental values 
    for (i in 1:ncol(prop_inc)) {
      prop_inc_mean[i,"distance"] <- distance[[i]]
      prop_inc_mean[i,"mean"] <- CI(prop_inc[[i]], ci=0.95)["mean"]
      prop_inc_mean[i,"ll"] <- CI(prop_inc[[i]], ci=0.95)["lower"]
      prop_inc_mean[i,"ul"] <- CI(prop_inc[[i]], ci=0.95)["upper"]
    }
    
    ## dot-and-whiskers plot 
    dotplot_prop_cum_Plim <- ggplot(data=prop_cum_mean,
                                   aes(x = distance,y = mean, ymin = ll, ymax = ul ))+
      geom_pointrange(color="black", size=0.3)+
      xlab('Distance interval (km)')+ ylab("Mean proportion of kin-connected houses") +
      geom_errorbar(aes(ymin=ll, ymax=ul),linewidth=0.3,cex=1) +
      geom_vline(xintercept=1.5, linetype="dashed", lwd=0.5, colour="red") +
      geom_hline(yintercept=0,  linetype="dashed", lwd=0.5, colour="black")
    
    
    dotplot_prop_inc_Plim <- ggplot(data=prop_inc_mean,
                                   aes(x = distance,y = mean, ymin = ll, ymax = ul ))+
      geom_pointrange(color="black", size=0.3)+
      xlab('Distance interval (km)')+ ylab("Mean proportion of kin-connected houses") +
      geom_errorbar(aes(ymin=ll, ymax=ul),linewidth=0.3,cex=1) +
      geom_vline(xintercept=1.5, linetype="dashed", lwd=0.5, colour="red") +
      geom_hline(yintercept=0,  linetype="dashed", lwd=0.5, colour="black")

# Save Map & Plots
    
    ## create theme for plotting
    min_theme <- theme_bw() + theme(
      axis.text.x=element_text(size=16),
      axis.text.y=element_text(size=16),
      axis.title=element_text(size=20),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="none",
      plot.title = element_text(size=20)
    )
    
    png("figures/network/spatial/map_hhx_Plim.png",width=8,height=6,units="in",res=600)
    map_hhx_Plim
    dev.off()
    
    png("figures/network/spatial/scatterplot_spatialhhx_cum_Plim.png", width=16,height=6,units="in",res=600)
    scatterplot_spatialhhx_cum_Plim + min_theme
    dev.off()
    
    png("figures/network/spatial/violinplot_spatialhhx_inc_Plim.png",width=16,height=6,units="in",res=600)
    violinplot_spatialhhx_inc_Plim + min_theme
    dev.off()
    
    png("figures/network/spatial/dotplot_prop_cum_Plim.png", width=16,height=6,units="in",res=600)
    dotplot_prop_cum_Plim + min_theme
    dev.off()
    
    png("figures/network/spatial/dotplot_prop_inc_Plim.png",width=16,height=6,units="in",res=600)
    dotplot_prop_inc_Plim + min_theme
    dev.off()

    
# What proportion of household-level kinship ties span a distance greater than 1.5km?

    ## create vector of distances to explore
    max(D)
    distance <- seq(0.1, 17, 0.1) # Fernandez 2019: min cell size/distance of analysis was 200m assuming that each household had at least 3 neighbors; the maximum distance was set at 6km (i.e., half of the dimension of the area)
    
    # (2b) What proportion of those (connected houses) are connected by kinship at various distances?
    
    ## create empty df where columns represent various distances
    soc_cum <- data.frame(matrix(nrow = ncol(D), ncol = length(distance)))
    soc_inc <- data.frame(matrix(nrow = ncol(D), ncol = length(distance)))
    
    ## count the cumulative and incremental number of houses that fall w/in each range for each house 
    counter = 0
    for (j in distance) {
      counter = counter + 1
      for(i in 1:ncol(D)) {
        neighbors = which(D[,i] < j & matrix_kin$P_lim[,i] >= 1, arr.ind=T)
        soc_cum[i,counter] <- length(as.vector(neighbors))  
        colnames(soc_cum)[counter] <- paste0(j)
        #subtract 1 b/c distance between house and itself is 0
      }
      
      if (counter>1) {
        soc_inc[,counter] = soc_cum[,counter] - soc_cum[,counter-1]
      } else {
        soc_inc[,counter] = soc_cum[,counter]
      }
      colnames(soc_inc) = colnames(soc_cum)
    }  
    
    ## get number of ties at each distance
    soc_inc_total <- as.data.frame(colSums(soc_inc))
    soc_inc_total$distance <- rownames(soc_inc_total)
    colnames(soc_inc_total) <- c("kintie", "distance")
    
    
# What proportion of kinship ties are long range (>1.5km ~max flight range of triatomine)
    
    ## create short/long range variable
    soc_inc_total$range <- ifelse(soc_inc_total$distance <= 1.5, "short", "long")
    
    ## get number of long range kinship ties
    long_kintie <- sum(soc_inc_total$kintie[soc_inc_total$distance>1.5])
    
    ## get total number of kinship ties
    total_kintie <- sum(soc_inc_total$kintie)
    
    ## get proportion of kinship ties that are long range
    print(long_kintie/total_kintie) #0.2445652
    

################################################################################
### PART 3. ASSORTATIVE MIXING IN THE KINSHIP NETWORK
################################################################################

# Setup
    
    ## clean environment
    rm(list=ls()) 
    graphics.off()
    
    ## set directory
    setwd("/Users/~/Desktop/ChagasNetworkPublic")
    
    ## load data
    load('data/CleanData.RData')
    
# Get matrix K
P <- as.matrix(unlist(matrix_kin$P))

# Build graph from matrix K (how connected are these terms; primary and secondary relations)
library(igraph)
g <- graph_from_adjacency_matrix(P, weighted=T, mode='undirected')

# Create layouts for plotting network
layout_xy <- cbind(df_kin$long, df_kin$lat) #geocoordinates
layout_fr <- layout_with_fr(g) #fruchterman reingold

# Get outline of PdI communities for plotting network
    
    ## load shape file
    library(sf)
    sf_comms <- read_sf("/Users/katietseng/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Fernandez Lab/Projects (Active)/Chagas Network/GitHub/ChagasNetwork/data/Shapefiles PPI/communities.shp")
    
    ## subset shape file to polygons of Pampa del Indio
    pdi <- sf_comms %>% dplyr::filter(tipo == 3) #subset shape file to polygons of PdI
    plot(st_geometry(pdi), lwd=2, col = c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2")) #test plot
    
    ## fix invalid geometries (e.g., polygon edges  that cross over each other)
    pdi_valid <- st_make_valid(pdi)
    pdi_valid <- st_as_sf(pdi_valid)  # important to keep sf class
    
    ## project to UTM zone 21S (EPSG:32721) - suitable for your area
    pdi_utm <- st_transform(pdi_valid, crs = 32721)
    
    ## dissolve all polygons into one geometry
    pdi_outline <- pdi_utm %>%
      st_combine() %>%  # combine geometries but keep class sf
      st_union()        # dissolve boundaries
    
    # transform from sfc geometry to sf
    plot(pdi_outline)
    
    # Assume pdi_outline_projected is in a projected CRS (e.g., UTM)
    pdi_outline_geographic <- st_transform(pdi_outline, crs = 4326)
    # png("figures/network/pdi_outline.png", width = 6.3, height = 10, units = "in", res = 600, bg = "transparent")
    plot(pdi_outline_geographic)
    # dev.off()
    
# Assign sociodemographic, ecological, and serological (epi) attributes
  
    # check that all vertex names match a row in our data
    sum(V(g)$name %in% df_kin$hhx.id)
    
    # set rownames in df_kin
    rownames(df_kin) <- df_kin$hhx.id
    
    # assign attributes including centrality measures
    names(df_kin)
    V(g)$long <- df_kin[V(g)$name, "long"]
    V(g)$lat <- df_kin[V(g)$name, "lat"]
    V(g)$community <- df_kin[V(g)$name, "community"]
    V(g)$serology <- df_kin[V(g)$name, "serology"]
    V(g)$serology.occ <- df_kin[V(g)$name, "serology.occ"]
    V(g)$serology.cat <- df_kin[V(g)$name, "serology.cat"]
    V(g)$serology.child <- df_kin[V(g)$name, "serology.child"]
    V(g)$serology.child.occ <- df_kin[V(g)$name, "serology.child.occ"]
    V(g)$ethnicity <- df_kin[V(g)$name, "ethnicity"]
    V(g)$mobility <- df_kin[V(g)$name, "mobility"]
    V(g)$svi <- df_kin[V(g)$name, "svi"]
    V(g)$hai <- df_kin[V(g)$name, "hai"]
    V(g)$v.inf <- df_kin[V(g)$name, "v.inf"]
    V(g)$iva.inf <- df_kin[V(g)$name, "iva.inf"]
    V(g)$iva <- df_kin[V(g)$name, "iva"]
    V(g)$degree <- df_kin[V(g)$name, "degree"]
    V(g)$degree.cat <- df_kin[V(g)$name, "degree.cat"]
    V(g)$betweenness <- df_kin[V(g)$name, "betweenness"]
    V(g)$betweenness.norm <- df_kin[V(g)$name, "betweenness.norm"]
    V(g)$harmonic <- df_kin[V(g)$name, "harmonic"]
    V(g)$harmonic.norm <- df_kin[V(g)$name, "harmonic.norm"]

    # alt degree
    degree(g)
           
# Create layout based on location data
layout_xy <- cbind(V(g)$long, V(g)$lat)

# Visualize
    
    ## Serology.cat - ALL
        
        ## check levels of serology.cat
        levels(V(g)$serology.cat)

        ## assign NA as its own category
        V(g)$serology.cat <- as.character(V(g)$serology.cat)
        V(g)$serology.cat[is.na(V(g)$serology.cat)] <- "missing"
        V(g)$serology.cat <- factor(V(g)$serology.cat, levels = c("0", "1-2", "3+", "missing"))
        
        ## set size based on svi_cat
        size_map <- c("0" = 4, "1-2" = 4, "3+" = 4, "missing" = 4)
        V(g)$size <- size_map[as.character(V(g)$serology.cat)]
        
        ## set node fill color and frame color based on data availability
        color_map <- c("0" = "black", "1-2" = "#FF7F7F", "3+" = "#9C0005", "missing" = "gray")
        V(g)$color <- color_map[as.character(V(g)$serology.cat)]
        frame_map <- c("0" = "black", "1-2" = "#FF7F7F", "3+" = "#9C0005", "missing" = "gray")
        V(g)$frame <- frame_map[as.character(V(g)$serology.cat)]
        
        ## plot
        png("figures/network/assortative/g_serology.png",width=6.3,height=10,units="in",res=600)
        set.seed(4321)
        par(mar = c(20, 0, 0, 1), xpd = NA)
        plot(g,
             layout = layout_with_fr(g),
             edge.width = 2,
             edge.color = "black",
             vertex.color = V(g)$color,
             vertex.frame.color = V(g)$frame,
             vertex.label = NA,
             vertex.size = V(g)$size)
        legend("bottom", 
               legend = c("0", "1-2", "3+", "No information"),
               pch = 21, 
               pt.bg = c("black", "#FF7F7F", "#9C0005", "gray"),
               col = c("black", "#FF7F7F", "#9C0005", "gray"),
               pt.cex = c(2, 2, 2, 2), bty = "n",
               x.intersp = 1.5, y.intersp = 1.5,
               cex = 1.5,
               inset = c(0, -0.6),  # Push legend downward outside plot area
               title = "Number of seropositive \npersons per household")
        dev.off()
        
        ## plot with PdI outline
        # png("figures/network/assortative/g_serology_xy.png",width=6.3,height=10,units="in",res=600)
        
            plot(g,
                 layout = layout_xy,
                 # layout = layout_with_fr(g),
                 vertex.color = V(g)$color,
                 vertex.frame.color = V(g)$frame,
                 vertex.label = NA,
                 vertex.size = V(g)$size,
                 edge.color = "gray80")

            # Increase right margin to make room for the legend
            par(mar=c(5, 4, 4, 12))

            # Plot polygon first (outline only)
            plot(pdi_outline_geographic, col = NA, border = "black", lwd = 2)

            # Plot edges manually:
            edges <- as_edgelist(g)
                name_to_index <- setNames(seq_along(V(g)), V(g)$name) # Map vertex names to numeric indices:
                for (i in 1:nrow(edges)) {
                  v1 <- name_to_index[edges[i,1]]
                  v2 <- name_to_index[edges[i,2]]
                  segments(x0 = layout_xy[v1, 1], y0 = layout_xy[v1, 2],
                           x1 = layout_xy[v2, 1], y1 = layout_xy[v2, 2],
                           col = "gray80")
                }

            # Plot vertices on top
            points(layout_xy, col = V(g)$color, pch = 21, bg = V(g)$color, cex = V(g)$size/6) 

            # Add legend
            legend("right", legend = c("0", "1-2", "3+", "No information"),
                   pch = 21,
                   pt.bg = c("black", "red", "darkred", "gray"),
                   col = c("black", "red", "darkred", "gray"),
                   pt.cex = c(2, 4, 6, 2), bty = "n",
                   x.intersp = 2, y.intersp = 2,
                   inset = c(-0.3, 0), xpd = TRUE, #shift legend to the right
                   title = "Number of seropositive \npersons per household")

        # dev.off()
            
        ## reset serology attribute
        V(g)$serology.cat <- df_kin[V(g)$name, "serology.cat"]

## Serology.occ - ALL
        
        ## check values
        table(V(g)$serology.occ, useNA="ifany")
        
        ## assign NA as its own category
        V(g)$serology.occ <- as.character(V(g)$serology.occ)
        V(g)$serology.occ[is.na(V(g)$serology.occ)] <- "missing"
        V(g)$serology.occ <- factor(V(g)$serology.occ, levels = c("No occurrence", "Occurrence", "missing"))
        
        ## set size based on svi_cat
        size_map <- c("No occurrence" = 4, "Occurrence" = 4, "missing" = 4)
        V(g)$size <- size_map[as.character(V(g)$serology.occ)]
        
        ## set node fill color and frame color based on data availability
        color_map <- c("No occurrence" = "black", "Occurrence" = "tomato", "missing" = "gray")
        V(g)$color <- color_map[as.character(V(g)$serology.occ)]
        frame_map <- c("No occurrence" = "black", "Occurrence" = "tomato", "missing" = "gray")
        V(g)$frame <- frame_map[as.character(V(g)$serology.occ)]
        
        ## plot
        png("figures/network/assortative/g_serologyocc.png",width=6.3,height=10,units="in",res=600)
        set.seed(4321)
        par(mar = c(20, 0, 0, 1), xpd = NA)
        plot(g,
             layout = layout_with_fr(g),
             edge.width = 2,
             edge.color = "black",
             vertex.color = V(g)$color,
             vertex.frame.color = V(g)$frame,
             vertex.label = NA,
             vertex.size = V(g)$size)
        legend("bottom", 
               legend = c("No occurrence", "Occurrence", "No information"),
               pch = 21, 
               pt.bg = c("black", "tomato", "gray"),
               col = c("black", "tomato", "gray"),
               pt.cex = c(2, 2, 2), bty = "n",
               x.intersp = 1.5, y.intersp = 1.5,
               cex = 1.5,
               inset = c(0, -0.6),  # Push legend downward outside plot area
               title = "Occurrence of at least \none seropositive inhabitant")
        dev.off()
            
            ## plot w/ geocoordinate layout 
                
                png("figures/network/assortative/g_xy.png",width=6.3,height=10,units="in",res=600)
        
                    # Plot polygon first (outline only)
                    plot(pdi_outline_geographic, col = NA, border = "black", lwd = 2)
                    
                    # Plot edges manually:
                    edges <- as_edgelist(g)
                    name_to_index <- setNames(seq_along(V(g)), V(g)$name) # Map vertex names to numeric indices:
                    for (i in 1:nrow(edges)) {
                      v1 <- name_to_index[edges[i,1]]
                      v2 <- name_to_index[edges[i,2]]
                      segments(x0 = layout_xy[v1, 1], y0 = layout_xy[v1, 2],
                               x1 = layout_xy[v2, 1], y1 = layout_xy[v2, 2],
                               col = "gray80")
                    }
                    
                    # Plot vertices on top
                    points(layout_xy, col = "black", pch = 21, bg = "black", cex = V(g)$size/5) 

                dev.off()
            
            ## plot w/ geocoordinate layout & serologyocc
                
                png("figures/network/assortative/g_serologyocc_xy.png",width=6.3,height=10,units="in",res=600)
                
                    # Plot polygon first (outline only)
                    plot(pdi_outline_geographic, col = NA, border = "black", lwd = 2)
                    
                    # Plot edges manually:
                    edges <- as_edgelist(g)
                    name_to_index <- setNames(seq_along(V(g)), V(g)$name) # Map vertex names to numeric indices:
                    for (i in 1:nrow(edges)) {
                      v1 <- name_to_index[edges[i,1]]
                      v2 <- name_to_index[edges[i,2]]
                      segments(x0 = layout_xy[v1, 1], y0 = layout_xy[v1, 2],
                               x1 = layout_xy[v2, 1], y1 = layout_xy[v2, 2],
                               col = "gray80")
                    }
                    
                    # Plot vertices on top
                    points(layout_xy, col = V(g)$color, pch = 21, bg = V(g)$color, cex = V(g)$size/5)
                    
                dev.off()
                
        ## reset serology attribute
        V(g)$serology.occ <- df_kin[V(g)$name, "serology.occ"]
        
    ## Serology.child.occ - CHILDREN
        
        ## check values
        table(V(g)$serology.child.occ, useNA="ifany")
        
        ## assign NA as its own category
        V(g)$serology.child.occ <- as.character(V(g)$serology.child.occ)
        V(g)$serology.child.occ[is.na(V(g)$serology.child.occ)] <- "missing"
        V(g)$serology.child.occ[V(g)$serology.child.occ=="No child"] <- "missing"
        V(g)$serology.child.occ <- factor(V(g)$serology.child.occ, levels = c("No occurrence", "Occurrence", "missing"))
        
        ## set size based on svi_cat
        size_map <- c("No occurrence" = 4, "Occurrence" = 4, "missing" = 4)
        V(g)$size <- size_map[as.character(V(g)$serology.child.occ)]
        
        ## set node fill color and frame color based on data availability
        color_map <- c("No occurrence" = "black", "Occurrence" = "tomato", "missing" = "gray")
        V(g)$color <- color_map[as.character(V(g)$serology.child.occ)]
        frame_map <- c("No occurrence" = "black", "Occurrence" = "tomato", "missing" = "gray")
        V(g)$frame <- frame_map[as.character(V(g)$serology.child.occ)]
        
        ## plot
        png("figures/network/assortative/g_serologychild.png",width=6.3,height=10,units="in",res=600)
        set.seed(4321)
        par(mar = c(20, 0, 0, 1), xpd = NA)
        plot(g,
             layout = layout_with_fr(g),
             edge.width = 2,
             edge.color = "black",
             vertex.color = V(g)$color,
             vertex.frame.color = V(g)$frame,
             vertex.label = NA,
             vertex.size = V(g)$size)
        legend("bottom", 
               legend = c("No occurrence", "Occurrence", "No child or \nno information"),
               pch = 21, 
               pt.bg = c("black", "tomato", "gray"),
               col = c("black", "tomato", "gray"),
               pt.cex = c(2, 2, 2), bty = "n",
               x.intersp = 1.5, y.intersp = 1.5,
               cex = 1.5,
               inset = c(0, -0.6),  # Push legend downward outside plot area
               title = "Occurrence of at least \none seropositive child")
        dev.off()

        ## plot w/ geocoordinate layout 
        
            png("figures/network/assortative/g_serologychild_xy.png",width=6.3,height=10,units="in",res=600)
                
                # Plot polygon first (outline only)
                plot(pdi_outline_geographic, col = NA, border = "black", lwd = 2)
                
                # Plot edges manually:
                edges <- as_edgelist(g)
                name_to_index <- setNames(seq_along(V(g)), V(g)$name) # Map vertex names to numeric indices:
                for (i in 1:nrow(edges)) {
                  v1 <- name_to_index[edges[i,1]]
                  v2 <- name_to_index[edges[i,2]]
                  segments(x0 = layout_xy[v1, 1], y0 = layout_xy[v1, 2],
                           x1 = layout_xy[v2, 1], y1 = layout_xy[v2, 2],
                           col = "gray80")
                }
                
                # Plot vertices on top
                points(layout_xy, col = V(g)$color, pch = 21, bg = V(g)$color, cex = V(g)$size/5)
                
            dev.off()        
              
        ## reset serology attribute
        V(g)$serology.child.occ <- df_kin[V(g)$name, "serology.child.occ"]
        
    ## Ethnicity
        
        ## set node fill color and frame color
        color_map <- c("Creole" = "maroon1", "Qom" = "maroon4")
        V(g)$color <- color_map[as.character(V(g)$ethnicity)]
        frame_map <- c("Creole" = "maroon1", "Qom" = "maroon4")
        V(g)$frame <- frame_map[as.character(V(g)$ethnicity)]
        
        ## plot
        png("figures/network/assortative/g_ethnicity.png",width=6.3,height=10,units="in",res=600)
        set.seed(4321)
        par(mar = c(20, 0, 0, 1), xpd = NA)
        plot(g, 
             layout = layout_with_fr(g),
             edge.width = 2,
             edge.color = "black",
             vertex.label=NA, vertex.shape="circle", vertex.size=4,
             vertex.color=V(g)$color,
             vertex.frame.color = V(g)$frame)
        legend("bottom", 
               legend = c("Creole", "Qom"),
               pch = 21, pt.bg = c("maroon1", "maroon4"),
               col = c("maroon1", "maroon4"),
               pt.cex = 2, bty = "n",
               x.intersp = 1.5, y.intersp = 1.5,
               cex = 1.5,
               inset = c(0, -0.4),  # Push legend downward outside plot area
               title = "Ethnicity")
        dev.off()

    ## Mobility
        
        ## check levels of mobility
        levels(V(g)$mobility)
        which(is.na(V(g)$mobility))
        
        ## assign NA as its own category
        V(g)$mobility <- as.character(V(g)$mobility)
        V(g)$mobility[is.na(V(g)$mobility)] <- "missing"
        V(g)$mobility <- factor(V(g)$mobility, levels = c("non-mover", "mover", "missing"))
        
        ## set node fill color and frame color
        color_map <- c("non-mover" = "sienna1", "mover" = "sienna4", "missing" = "gray")
        V(g)$color <- color_map[as.character(V(g)$mobility)]
        frame_map <- c("non-mover" = "sienna1", "mover" = "sienna4", "missing" = "gray")
        V(g)$frame <- frame_map[as.character(V(g)$mobility)]
        
        ## plot
        png("figures/network/assortative/g_mobility.png",width=6.3,height=10,units="in",res=600)
        set.seed(4321)
        par(mar = c(20, 0, 0, 1), xpd = NA)
        plot(g, 
             # layout = layout_xy,
             layout = layout_with_fr(g),
             edge.width = 2,
             edge.color = "black",
             vertex.label=NA, vertex.shape="circle", vertex.size=4,
             vertex.color=V(g)$color,
             vertex.frame.color = V(g)$frame)
        legend("bottom", 
               legend = c("Non-mover", "Mover", "No information"),
               pch = 21, pt.bg = c("sienna1", "sienna4", "gray"),
               col = c("sienna1", "sienna4", "gray"), 
               pt.cex = 2, bty = "n",
               x.intersp = 1.5, y.intersp = 1.5,
               cex = 1.5,
               inset = c(0, -0.45),  # Push legend downward outside plot area
               title = "Mobility")
        dev.off()
        
        ## reset mobility attribute
        V(g)$mobility <- df_kin[V(g)$name, "mobility"]

    ## Social vulnerability
        
        ## compute quantiles
        svi_quants <- quantile(V(g)$svi, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
        
        ## make categorical using quantiles
        V(g)$svi_cat <- cut(V(g)$svi, breaks = svi_quants, include.lowest = TRUE, #create cats
                            labels = c("low", "mid", "high"))
        
        ## assign NA as its own category
        V(g)$svi_cat <- as.character(V(g)$svi_cat)
        V(g)$svi_cat[is.na(V(g)$svi_cat)] <- "missing"
        V(g)$svi_cat <- factor(V(g)$svi_cat, levels = c("low", "mid", "high", "missing"))
        
        ## set size based on svi_cat
        size_map <- c("low" = 2, "mid" = 4, "high" = 6, "missing" = 2)
        V(g)$size <- size_map[as.character(V(g)$svi_cat)]
        
        ## set node fill color and frame color based on data availability
        color_map <- c("low" = "lightblue", "mid" = "steelblue3", "high" = "royalblue4", "missing" = "gray")
        V(g)$color <- color_map[as.character(V(g)$svi_cat)]
        frame_map <- c("low" = "lightblue", "mid" = "steelblue3", "high" = "royalblue4", "missing" = "gray")
        V(g)$frame <- frame_map[as.character(V(g)$svi_cat)]
        
        ## plot
        png("figures/network/assortative/g_svicat.png",width=6.3,height=10,units="in",res=600)
        set.seed(4321)
        par(mar = c(20, 0, 0, 1), xpd = NA)
        plot(g,
             # layout = layout_xy,
             layout = layout_with_fr(g),
             edge.width = 2,
             edge.color = "black",
             vertex.color = V(g)$color,
             vertex.frame.color = V(g)$frame,
             vertex.label = NA,
             vertex.size = V(g)$size)
        legend("bottom", 
               legend = c("Low (0-33%)", "Mid (34-66%)", "High (67-100)%", "No information"),
               pch = 21, 
               pt.bg = c("lightblue", "steelblue3", "royalblue4", "gray"),
               col = c("lightblue", "steelblue3", "royalblue4", "gray"), 
               pt.cex = c(1, 2, 3, 1), bty = "n",
               x.intersp = 1.5, y.intersp = 1.5,
               cex = 1.5,
               inset = c(0, -0.6),  # Push legend downward outside plot area
               title = "Social vulnerability index \n(percentile ranges)")
        dev.off()
        
    ## Host availability
        
        ## compute quantiles
        hai_quants <- quantile(V(g)$hai, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
        
        ## make categorical using quantiles
        V(g)$hai_cat <- cut(V(g)$hai, breaks = hai_quants, include.lowest = TRUE, #create cats
                            labels = c("low", "mid", "high"))
        
        ## assign NA as its own category
        V(g)$hai_cat <- as.character(V(g)$hai_cat)
        V(g)$hai_cat[is.na(V(g)$hai_cat)] <- "missing"
        V(g)$hai_cat <- factor(V(g)$hai_cat, levels = c("low", "mid", "high", "missing"))
        
        ## set size based on hai_cat
        size_map <- c("low" = 2, "mid" = 4, "high" = 6, "missing" = 2)
        V(g)$size <- size_map[as.character(V(g)$hai_cat)]
        
        ## set node fill color and frame color based on data availability
        color_map <- c("low" = "plum2", "mid" = "darkorchid1", "high" = "darkmagenta", "missing" = "gray")
        V(g)$color <- color_map[as.character(V(g)$hai_cat)]
        frame_map <- c("low" = "plum2", "mid" = "darkorchid1", "high" = "darkmagenta", "missing" = "gray")
        V(g)$frame <- frame_map[as.character(V(g)$hai_cat)]
        
        ## plot
        png("figures/network/assortative/g_haicat.png",width=6.3,height=10,units="in",res=600)
        set.seed(4321)
        par(mar = c(18, 0, 0, 1), xpd = NA)
        plot(g,
             # layout = layout_xy,
             layout = layout_with_fr(g),
             edge.width = 2,
             edge.color = "black",
             vertex.color = V(g)$color,
             vertex.frame.color = V(g)$frame,
             vertex.label = NA,
             vertex.size = V(g)$size)
        legend("bottom", 
               legend = c("Low (0-33%)", "Mid (34-66%)", "High (67-100)%", "No information"),
               pch = 21, 
               pt.bg = c("plum2", "darkorchid1", "darkmagenta", "gray"),
               col = c("plum2", "darkorchid1", "darkmagenta", "gray"), 
               pt.cex = c(1, 2, 3, 1), bty = "n",
               x.intersp = 1.5, y.intersp = 1.5,
               cex = 1.5,
               inset = c(0, -0.5),  # Push legend downward outside plot area
               title = "Host availability index \n(percentile ranges)")
        dev.off()
        
    ## Infestation with vector
        
        ## assign NA as its own category
        V(g)$v.inf <- as.character(V(g)$v.inf)
        V(g)$v.inf[is.na(V(g)$v.inf)] <- "missing"
        V(g)$v.inf <- factor(V(g)$v.inf, levels = c("non-infested", "infested", "missing"))
        
        ## set node fill color and frame color
        color_map <- c("non-infested" = "lightgreen", "infested" = "darkgreen", "missing" = "gray")
        V(g)$color <- color_map[as.character(V(g)$v.inf)]
        frame_map <- c("non-infested" = "lightgreen", "infested" = "darkgreen", "missing" = "gray")
        V(g)$frame <- frame_map[as.character(V(g)$v.inf)]
        
        ## plot
        png("figures/network/assortative/g_v.inf.png",width=6.3,height=10,units="in",res=600)
        set.seed(4321)
        par(mar = c(20, 0, 0, 1), xpd = NA)
        plot(g,
             layout = layout_with_fr(g),
             edge.width = 2,
             edge.color = "black",
             vertex.color = V(g)$color,
             vertex.frame.color = V(g)$frame,
             vertex.label = NA,
             vertex.size = 4)
        legend("bottom", 
               legend = c("Non-infested", "Infested", "No information"),
               pch = 21, 
               pt.bg = c("lightgreen", "darkgreen", "gray"),
               col = c("lightgreen", "darkgreen", "gray"),
               pt.cex = 2, bty = "n",
               x.intersp = 1.5, y.intersp = 1.5,
               cex = 1.5,
               inset = c(0, -0.45),  # Push legend downward outside plot area
               title = "Occurrence of infestation")
        dev.off()
        
        ## reset v.inf attribute
        V(g)$v.inf <- df_kin[V(g)$name, "v.inf"]

    ## Infestation with INFECTED vector
        
        ## assign NA as its own category
        V(g)$iva.inf <- as.character(V(g)$iva.inf)
        V(g)$iva.inf[is.na(V(g)$iva.inf)] <- "missing"
        V(g)$iva.inf <- factor(V(g)$iva.inf, levels = c("0", "1", "missing"))
        
        ## set node fill color and frame color
        color_map <- c("0" = "#ACFFFC", "1" = "cyan3", "missing" = "gray")
        V(g)$color <- color_map[as.character(V(g)$iva.inf)]
        frame_map <- c("0" = "#ACFFFC", "1" = "cyan3", "missing" = "gray")
        V(g)$frame <- frame_map[as.character(V(g)$iva.inf)]
        
        ## plot
        png("figures/network/assortative/g_iva.inf.png",width=6.3,height=10,units="in",res=600)
        set.seed(4321)
        par(mar = c(20, 0, 0, 1), xpd = NA)
        plot(g,
             layout = layout_with_fr(g),
             edge.width = 2,
             edge.color = "black",
             vertex.color = V(g)$color,
             vertex.frame.color = V(g)$frame,
             vertex.label = NA,
             vertex.size = 4)
        legend("bottom", 
               legend = c("Non-infested", "Infested", "No information"),
               pch = 21, 
               pt.bg = c("#ACFFFC", "cyan3", "gray"),
               col = c("#ACFFFC", "cyan3", "gray"),
               pt.cex = 2, bty = "n",
               x.intersp = 1.5, y.intersp = 1.5,
               cex = 1.5,
               inset = c(0, -0.55),  # Push legend downward outside plot area
               title = "Occurrence of infestation \nwith T. cruzi-infected \nT. infestans")
        dev.off()
        
        ## reset iva.inf attribute
        V(g)$iva.inf <- df_kin[V(g)$name, "iva.inf"]

    ## Connected vs. not connected
        
        ## assign attribute
        V(g)$connected <- df_kin[V(g)$name, "connected"]

        ## check values
        table(V(g)$connected, useNA="ifany")
        
        ## set node fill color and frame color based on data availability
        color_map <- c("Unconnected" = "gray", "Connected" = "#C66B39")
        V(g)$color <- color_map[as.character(V(g)$connected)]
        frame_map <- c("Unconnected" = "gray", "Connected" = "#C66B39")
        V(g)$frame <- frame_map[as.character(V(g)$connected)]
        
        ## plot
        png("figures/network/assortative/g_connected.png",width=6.3,height=10,units="in",res=600)
        set.seed(4321)
        par(mar = c(20, 0, 0, 1), xpd = NA)
        plot(g,
             layout = layout_with_fr(g),
             edge.width = 2,
             edge.color = "black",
             vertex.color = V(g)$color,
             vertex.frame.color = V(g)$frame,
             vertex.label = NA,
             vertex.size = 4)
        legend("bottom", 
               legend = c("No parent-child ties", "At least one parent-child tie"),
               pch = 21, 
               pt.bg = c("gray", "#C66B39"),
               col = c("gray", "#C66B39"),
               pt.cex = 2, bty = "n",
               x.intersp = 1.5, y.intersp = 1.5,
               cex = 1.5,
               inset = c(0, -0.45)  # Push legend downward outside plot area
               )
        dev.off()
        
# Assortativity tests
        
    ## Save results
    df_assortativity <- data.frame(
      variable = character(),
      r = numeric(),
      stringsAsFactors = FALSE
    )
    
    ## review attributes of households
    vertex_attr_names(g)
    
    ## Serology: Do households with higher numbers of seropositive individuals tend to be connected to other similarly infected households within the kinship network?
    g_clean <- delete_vertices(g, which(is.na(V(g)$serology)))
    assort_val <- round(assortativity(g_clean, V(g_clean)$serology, directed = FALSE), 2)
    df_assortativity <- rbind(df_assortativity, data.frame(variable = "No. of seropositive", r = assort_val))
    
    ## Serology.occ: Do households with at least one seropositive inhabitant tend to be connected to other similarly infected households?
    g_clean <- delete_vertices(g, which(is.na(V(g)$serology.occ)))
    assort_val <- round(assortativity_nominal(g, types = V(g)$serology.occ, directed = FALSE), 2)
    df_assortativity <- rbind(df_assortativity, data.frame(variable = "Occ. of seropositive", r = assort_val))    

    ## Serology.child: Do households with higher numbers of seropositive children tend to be connected to other similarly infected households within the kinship network?
    g_clean <- delete_vertices(g, which(is.na(V(g)$serology.child)))
    assort_val <- round(assortativity(g_clean, V(g_clean)$serology.child, directed = FALSE), 2)
    df_assortativity <- rbind(df_assortativity, data.frame(variable = "No. of seropositive children", r = assort_val))
    
    ## Serology.child.occ: Do households with at least one seropositive child tend to be connected to other similarly infected households?
    g_clean <- delete_vertices(g, which(is.na(V(g)$serology.child.occ)))
    assort_val <- round(assortativity_nominal(g, types = V(g)$serology.child.occ, directed = FALSE), 2)
    df_assortativity <- rbind(df_assortativity, data.frame(variable = "Occ. of seropositive child", r = assort_val))    
    
    ## Ethnicitiy: Are Qom households more connected to other Qom households than expected by chance?
    assort_val <- round(assortativity_nominal(g, types = V(g)$ethnicity, directed = FALSE), 2)
    df_assortativity <- rbind(df_assortativity, data.frame(variable = "Ethnicity", r = assort_val))    
    
    ## Mobility: Are mover households more connected to other mover households than expected by chance?
    g_clean <- delete_vertices(g, which(is.na(V(g)$mobility)))
    assort_val <- round(assortativity_nominal(g, types = V(g)$mobility, directed = FALSE), 2)
    df_assortativity <- rbind(df_assortativity, data.frame(variable = "Mobility", r = assort_val))    
    
    ## SVI: Do high-SVI households tend to be connected to other high-SVI households?
    g_clean <- delete_vertices(g, which(is.na(V(g)$svi)))    
    assort_val <- round(assortativity(g_clean, V(g_clean)$svi, directed = FALSE), 2) #r=0.2302431
    df_assortativity <- rbind(df_assortativity, data.frame(variable = "Social vulnerability", r = assort_val))
    
    ## HAI: Do high-HAI households tend to be connected to other high-HAI households?
    g_clean <- delete_vertices(g, which(is.na(V(g)$hai)))
    assort_val <- round(assortativity(g_clean, V(g_clean)$hai, directed = FALSE), 2)
    df_assortativity <- rbind(df_assortativity, data.frame(variable = "Host availability", r = assort_val))
    
    ## Infestation: Are infested households more connected to other infested households than expected by chance?
    g_clean <- delete_vertices(g, which(is.na(V(g)$v.inf)))
    assort_val <- round(assortativity_nominal(g, types = V(g)$v.inf, directed = FALSE), 2)
    df_assortativity <- rbind(df_assortativity, data.frame(variable = "Infestation w/ vector", r = assort_val))    
    
    ## Infestation w/ infected vector: Are households infested with infected vector more connected to other similarly infested households than expected by chance?
    g_clean <- delete_vertices(g, which(is.na(as.factor(V(g)$iva.inf))))
    assort_val <- round(assortativity_nominal(g, types = as.factor(V(g)$iva.inf), directed = FALSE), 2)
    df_assortativity <- rbind(df_assortativity, data.frame(variable = "Infestation w/ infected vector", r = assort_val))    

################################################################################
### PART 4: NETWORK CENTRALITY PATTERNS AND DEGREE ASSORTATIVITY
################################################################################
#(1) Summary statistics: median (IQR), mean (SD), IQR - helps characterize the distribution and skewness (are a few households highly central while most are peripheral?)
#(2) Histograms to show how centrality values are distributed
#(3) Scatterplots to show relationship between centrality measures and correlations - are households that are locally well-connected also globally influential?
#(4) Network-level insights: characterize top-n-most central nodes (ethnicity, mobility, svi)
#(5) Spatial mapping of central nodes: do we see spatial clustering of inluential households?

# Summary statistics
    
    ## degree
    table(V(g)$degree)
    summary(V(g)$degree)
    summary(V(g)$degree[V(g)$degree!=0])
    hist(V(g)$degree)
    
    ## betweenness
    summary(V(g)$betweenness)
    summary(V(g)$betweenness[V(g)$betweenness!=0])
    hist(V(g)$betweenness)
    length(which(V(g)$betweenness==0))

    ## harmonic centrality
    summary(V(g)$harmonic)
    summary(V(g)$harmonic[V(g)$harmonic!=0])
    length(which(V(g)$harmonic==0))

# Summary statistics (for only households that have at least one kinship tie)
    
    ## get df
    df <- df_kin
    
    ## change degree to integer
    df$degree <- as.integer(df$degree)
    
    ## keep only those households with at least one kinship tie
    df <- df[df$degree>0,]
    
    ## create df for results
    df_centrality <- data.frame(
      variable = character(),
      mean_sd = character(),
      median_iqr = character(),
      stringsAsFactors = FALSE
    )
    
    for (var in c("degree", "betweenness", "harmonic")) {
      
      mean_val <- mean(df[[var]])
      sd_val <- sd(df[[var]])
      mean_sd <- paste0(round(mean_val, 1), " (", round(sd_val, 2), ")")
      
      median_val <- median(df[[var]])
      iqr_1st <- round(summary(df[[var]])["1st Qu."], 1)
      iqr_3rd <- round(summary(df[[var]])["3rd Qu."], 1)
      median_iqr = paste0(round(median_val, 1), " (",
                          round(iqr_1st, 2), "-", 
                          round(iqr_3rd, 2), ")")
      
      df_centrality <- rbind(df_centrality, data.frame(variable = var, mean_sd = mean_sd, median_iqr = median_iqr))
    }
    
    print(df_centrality)

# Histograms
    
    ## create my theme
    mytheme <- theme(axis.text=element_text(family = "Arial", size=12),
                     axis.title=element_text(family = "Arial", size=12),
                     plot.title = element_text(family = "Arial", size=14),
                     axis.title.y = element_text(margin = margin(r = 10)),
                     axis.title.x = element_text(margin = margin(t = 10)))

    ## degree
    hist_degree <- ggplot(df, aes(x = degree)) +
                          geom_histogram(binwidth = 0.5, color = "white", alpha = 0.7) +
                          labs(x = "Degree centrality", y = "Frequency") +
                          scale_x_continuous(breaks = seq(floor(min(df$degree)), ceiling(max(df$degree)), by = 1)) +
                          theme_minimal()
    
    ## betweenness
    hist_between <- ggplot(df, aes(x = betweenness)) +
                          geom_histogram(binwidth = 5, color = "white", alpha = 0.7) +
                          labs(x = "Betweenness centrality", y = "Frequency") +
                          # scale_x_continuous(breaks = seq(floor(min(df$betweenness)), ceiling(max(df$betweenness)), by = 1)) +
                          theme_minimal()
    
    ## harmonic closeness
    hist_harmonic <- ggplot(df, aes(x = harmonic)) +
                            geom_histogram(binwidth = 0.5, color = "white", alpha = 0.7) +
                            labs(x = "Harmonic centrality", y = "Frequency") +
                            # scale_x_continuous(breaks = seq(floor(min(df$degree)), ceiling(max(df$degree)), by = 1)) +
                            theme_minimal()

    ## plot together
    library(cowplot)
    png("figures/network/centrality/hist_centrality.png",width=12,height=4,units="in",res=600)
    plot_grid(hist_degree + mytheme,
              hist_between + mytheme,
              hist_harmonic + mytheme,
              nrow = 1, ncol = 3
    )
    dev.off()
    
# Are households that are locally well-connected also globally influential?    
    
    ## scatterplot degree vs. betweenness
    scatter_degree_between <- ggplot(df, aes(x = degree, y = betweenness)) +
                                    geom_point(alpha = 0.6) +
                                    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
                                    theme_minimal(base_size = 14) +
                                    labs(
                                      x = "Degree centrality",
                                      y = "Betweenness centrality",
                                    )
    
    ## scatterplot degree vs. harmonic closeness
    scatter_degree_harmonic <- ggplot(df, aes(x = degree, y = harmonic)) +
                                      geom_point(alpha = 0.6) +
                                      geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
                                      theme_minimal(base_size = 14) +
                                      labs(
                                        x = "Degree centrality",
                                        y = "Harmonic centrality",
                                      )
                                    
    ## scatterplot betweenness and harmonic closeness
    scatter_between_harmonic <- ggplot(df, aes(x = betweenness, y = harmonic)) +
                                       geom_point(alpha = 0.6) +
                                       geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
                                       theme_minimal(base_size = 14) +
                                       labs(
                                         x = "Betweenness centrality",
                                         y = "Harmonic centrality",
                                       )
    
    ## plot together
    library(cowplot)
    png("figures/network/centrality/scatter_centrality.png",width=12,height=4,units="in",res=600)
    plot_grid(scatter_degree_between + mytheme,
              scatter_degree_harmonic + mytheme,
              scatter_between_harmonic + mytheme,
              nrow = 1, ncol = 3
    )
    dev.off()
    
    ## are they correlated?
          
        cor(df$degree, df$betweenness)
        cor(df$degree, df$harmonic)
        cor(df$betweenness, df$harmonic)

# Degree Assortativity
    
    ## assign network centrality measures
    names(df_kin)
    V(g)$degree <- df_kin[V(g)$name, "degree"]
    V(g)$degree.cat <- df_kin[V(g)$name, "degree.cat"]
    V(g)$betweenness <- df_kin[V(g)$name, "betweenness"]
    V(g)$betweenness.norm <- df_kin[V(g)$name, "betweenness.norm"]
    V(g)$harmonic <- df_kin[V(g)$name, "harmonic"]
    V(g)$harmonic.norm <- df_kin[V(g)$name, "harmonic.norm"]
        
    ## test on the whole network (isolated nodes are ignored)    
    
        ## degree
        assort_val <- round(assortativity(g, V(g)$degree, direct = FALSE), 2)
        df_assortativity <- rbind(df_assortativity, data.frame(variable = "Degree centrality", r = assort_val))
        
        ## betweenness
        assort_val <- round(assortativity(g, V(g)$betweenness, direct = FALSE), 2)
        df_assortativity <- rbind(df_assortativity, data.frame(variable = "Betweenness centrality", r = assort_val))
        
        ## harmonic
        assort_val <- round(assortativity(g, V(g)$harmonic, direct = FALSE), 2)
        df_assortativity <- rbind(df_assortativity, data.frame(variable = "Harmonic centrality", r = assort_val))
        
    # test on the largest connected component
    
        ## add network variable to assortativity table
        df_assortativity$network <- "Full network"
        
        ## get the largest connected component
        comp <- components(g)
        lcc <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
    
        ## degree
        assort_val <- round(assortativity(lcc, V(lcc)$degree, direct = FALSE), 2)
        df_assortativity <- rbind(df_assortativity, data.frame(variable = "Degree centrality", r = assort_val, network = "Largest connected component"))
        
        ## betweenness
        assort_val <- round(assortativity(lcc, V(lcc)$betweenness, direct = FALSE), 2)
        df_assortativity <- rbind(df_assortativity, data.frame(variable = "Betweenness centrality", r = assort_val, network = "Largest connected component"))
        
        ## harmonic
        assort_val <- round(assortativity(lcc, V(lcc)$harmonic, direct = FALSE), 2)
        df_assortativity <- rbind(df_assortativity, data.frame(variable = "Harmonic centrality", r = assort_val, network = "Largest connected component"))
        
        ## save
        write.csv(df_assortativity, "tables/network/assortativity_results.csv")
        
# Assortativity Lollipop Plot
        
        # make variable factor
        df_assortativity$variable <- factor(df_assortativity$variable, 
                                            levels = c("Harmonic centrality", "Betweenness centrality", "Degree centrality",
                                                       "Infestation w/ infected vector", "Infestation w/ vector",
                                                       "Host availability", "Social vulnerability",
                                                       "Mobility", "Ethnicity",
                                                       "Occ. of seropositive child", "No. of seropositive children", "Occ. of seropositive", "No. of seropositive"))
        
        # remove seropositivity variables from plot
            
            ## filter
            df_assortativity <- df_assortativity %>%
              dplyr::filter(!variable %in% c(
                "Occ. of seropositive child",
                "No. of seropositive children",
                "Occ. of seropositive",
                "No. of seropositive"
              ))
            
            ## redefine factor
            df_assortativity$variable <- factor(df_assortativity$variable, 
                                                levels = c(
                                                  "Harmonic centrality", "Betweenness centrality", "Degree centrality",
                                                  "Infestation w/ infected vector", "Infestation w/ vector",
                                                  "Host availability", "Social vulnerability",
                                                  "Mobility", "Ethnicity"
                                                )
            )
            
        # plot
        png("figures/network/assortative/lolliplot_assortativity.png",width=8,height=8,units="in",res=600)
        ggplot(df_assortativity, aes(x = variable, y = r, color = network)) +
        geom_segment(aes(xend = variable, xend = variable, y = 0, yend = r),
                     position = position_dodge(width = 0.5)) +
        geom_point(size = 3, 
                   position = position_dodge(width = 0.5)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("gray40", "darkred")) +
        scale_y_continuous(breaks = seq(-0.8, 1.2, 0.2), labels = seq(-0.8, 1.2, 0.2)) +
        theme_minimal(base_size = 14) +
        labs(y = "Assortativity coefficient (r)",
             x = "Node attribute") +
        coord_flip()
        dev.off()
    
# Characterize top-n-most central nodes (ethnicity, mobility, svi)
    
    ## get top 10 indices for each metric
    df$highlight <- "None"
    top_degree <- df[order(df$degree, decreasing = TRUE)[1:10],]
    top_betweenness <- df[order(df$betweenness, decreasing = TRUE)[1:10],]
    top_harmonic <- df[order(df$harmonic, decreasing = TRUE)[1:10],]
    
    ## plot network
    ## occurrence of seropositive inhabitant - assign NA as its own category
    V(g)$occ.case <- as.character(V(g)$serology.occ)
    V(g)$occ.case[is.na(V(g)$occ.case)] <- "missing"
    V(g)$occ.case <- factor(V(g)$occ.case, levels = c("No occurrence", "Occurrence", "missing"))
 
# Plot centrality measures for the largest component
        
    ## get components
    comps <- components(g)
        
    ## get largest component and create the subgraph
    largest_comp_id <- which.max(comps$csize)
    g_largest <- induced_subgraph(g, which(comps$membership == largest_comp_id))
        
    ## get matching coordinates
    coords <- as.matrix(df_kin[match(V(g_largest)$name, df_kin$hhx.id), c("long", "lat")])
        
    ## prep centrality measures for plotting
        
    ## degree
        
        ## quantile bins for non-zero values (3 bins as example)
        summary(V(g_largest)$degree)
        
        ## assign categories
        V(g_largest)$degree.cat <- ifelse(V(g_largest)$degree == 1, 1,
                                          ifelse(V(g_largest)$degree == 2, 2, 3))
        
        
    ## betweenness
        
        ## quantile bins for non-zero values (3 bins as example)
        betweenness_nonzero <- ifelse(V(g_largest)$betweenness > 0, V(g_largest)$betweenness, NA)
        betweenness_breaks <- quantile(betweenness_nonzero, probs = c(0, 0.25, 0.75, 1), na.rm = TRUE)
        
        ## assign categories
        V(g_largest)$betweenness.cat <- ifelse(V(g_largest)$betweenness == 0, 0, 
                                               cut(betweenness_nonzero, breaks = betweenness_breaks, include.lowest = TRUE, 
                                                   labels = c("Low", "Medium", "High")))
        
        
    ## harmonic
        
        ## quantile bins for non-zero values (3 bins as example)
        harmonic_nonzero <- ifelse(V(g_largest)$harmonic > 0, V(g_largest)$harmonic, NA)
        harmonic_breaks <- quantile(harmonic_nonzero, probs = c(0, 0.25, 0.75, 1), na.rm = TRUE)
        
        # Assign other categories
        V(g_largest)$harmonic.cat <- cut(harmonic_nonzero, breaks = harmonic_breaks, include.lowest = TRUE,
                                         labels = c("Low", "Medium", "High")
        )
        
        # Transform to numeric
        V(g_largest)$harmonic.cat <- as.numeric(V(g_largest)$harmonic.cat)
        
    ## prep color variables (serological and infestation vars) for plotting
        
        ## serology: occurrence of seropositive inhabitant - assign NA as its own category
        table(V(g_largest)$serology.occ, useNA="ifany")
        V(g_largest)$serology.occ <- as.character(V(g_largest)$serology.occ)
        V(g_largest)$serology.occ[is.na(V(g_largest)$serology.occ)] <- "missing"
        V(g_largest)$serology.occ <- factor(V(g_largest)$serology.occ, levels = c("No occurrence", "Occurrence", "missing"))
        
        ## serology child
        table(V(g_largest)$serology.child.occ, useNA="ifany")
        V(g_largest)$serology.child.occ <- as.character(V(g_largest)$serology.child.occ)
        V(g_largest)$serology.child.occ[is.na(V(g_largest)$serology.child.occ)] <- "missing"
        V(g_largest)$serology.child.occ[V(g_largest)$serology.child.occ=="No child"] <- "missing"
        V(g_largest)$serology.child.occ <- factor(V(g_largest)$serology.child.occ, levels = c("No occurrence", "Occurrence", "missing"))
        
        ## infestation
        table(V(g_largest)$v.inf, useNA="ifany")
        V(g_largest)$v.inf <- as.character(V(g_largest)$v.inf)
        V(g_largest)$v.inf[is.na(V(g_largest)$v.inf)] <- "missing"
        V(g_largest)$v.inf <- factor(V(g_largest)$v.inf, levels = c("non-infested", "infested", "missing"))
        
    ## define color mappings and variable names
        
        ## color list
        color_map_list <- list(
          serology.occ = c("black", "red", "gray"),
          serology.child.occ = c("black", "pink", "gray"),
          v.inf = c("black", "springgreen", "gray")
        )
        
        ## variable names to use in loop
        color_vars <- names(color_map_list)
        
    ## loop through color variables and plot each centrality
     
          for (color_var in color_vars) {
          
          ## generate base color vector
          color_vals <- color_map_list[[color_var]]
          vertex_colors <- color_vals[as.numeric(as.factor(vertex_attr(g_largest, color_var)))]
          
          ## -------- Degree Centrality Plot --------
          png(paste0("figures/network/centrality/g_lcc_", color_var, "_degree.png"), width = 9, height = 9, units = "in", res = 600)
          par(mar = c(5, 5, 2, 5))  # set  margins to allow space for  legend
          set.seed(4321)
          plot(g_largest, layout = layout_with_fr(g_largest), edge.color = "gray50",
               vertex.label = NA,
               vertex.shape = "circle",
               vertex.size = (V(g_largest)$degree.cat) * 2.5,
               vertex.color = vertex_colors)
          legend("bottomleft",
                 legend = c("1", "2", "3+"),  # categories of vertex size
                 title = "Degree centrality \n(number of kinship ties)",
                 box.lwd = 0, 
                 pch = 21,  # shape of legend symbol
                 pt.cex = c(1, 2, 3),  # adjust size
                 pt.bg = "white",  # fill color
                 cex = 1, 
                 y.intersp = 1.5, x.intersp = 1.5,
                 xpd = TRUE)  
          dev.off()
          
          ## -------- Betweenness Centrality Plot --------
          png(paste0("figures/network/centrality/g_lcc_", color_var, "_betweenness.png"), width = 9, height = 9, units = "in", res = 600)
          par(mar = c(5, 5, 2, 5))
          set.seed(4321)
          plot(g_largest, layout = layout_with_fr(g_largest), edge.color = "gray50",
               vertex.label = NA,
               vertex.shape = "circle",
               vertex.size = (V(g_largest)$betweenness.cat + 1) * 2.5,
               vertex.color = vertex_colors)
          legend("bottomleft",
                 legend = c("0", "Low (0–25th]", "Medium (25–75th]", "High (75–100th]"),
                 title = "Betweenness centrality\n(percentile for nonzero values)",
                 box.lwd = 0,
                 pch = 21,
                 pt.cex = c(1, 2, 3, 4),
                 pt.bg = "white",
                 cex = 1,
                 y.intersp = 1.5, x.intersp = 1.5,
                 xpd = TRUE)
          dev.off()
          
          ## -------- Harmonic Centrality Plot --------
          png(paste0("figures/network/centrality/g_lcc_", color_var, "_harmonic.png"), width = 9, height = 9, units = "in", res = 600)
          par(mar = c(5, 5, 2, 5))
          set.seed(4321)
          plot(g_largest, layout = layout_with_fr(g_largest), edge.color = "gray50",
               vertex.label = NA,
               vertex.shape = "circle",
               vertex.size = V(g_largest)$harmonic.cat * 2.5,
               vertex.color = vertex_colors)
          legend("bottomleft",
                 legend = c("Low (0–25th]", "Medium (25–75th]", "High (75–100th]"),
                 title = "Harmonic centrality\n(percentile)",
                 box.lwd = 0,
                 pch = 21,
                 pt.cex = c(1, 2, 3),
                 pt.bg = "white",
                 cex = 1,
                 y.intersp = 1.5, x.intersp = 1.5,
                 xpd = TRUE)
          dev.off()
        } # end loop
        
        ## save legends for vertex color
        
        ## serology.occ
        png("figures/network/centrality/g_lcc_legend_serology.occ.png", width=6,height=2,units="in",res=600)
        par(mar = c(1, 1, 1, 1))
        plot.new()
        legend("bottomleft", 
               legend = levels(gl(3, 1, labels = c("No occurrence", 
                                                   "Occurrence", 
                                                   "No information"))), 
               title = "Occurrence of at least \none seropositive inhabitant",
               box.lwd = 0,
               fill = c("black", "red", "gray"))
        dev.off()
        
        ## serology.child.occ
        png("figures/network/centrality/g_lcc_legend_serology.child.occ.png", width=6,height=2,units="in",res=600)
        par(mar = c(1, 1, 1, 1))
        plot.new()
        legend("bottomleft", 
               legend = levels(gl(3, 1, labels = c("No occurrence", 
                                                   "Occurrence", 
                                                   "No information"))),
               title = "Occurrence of at least \none seropositive child",
               box.lwd = 0, 
               fill = c("black", "pink", "gray"))
        dev.off()
        
        ## v.inf
        png("figures/network/centrality/g_lcc_legend_v.inf.png", width=6,height=2,units="in",res=600)
        par(mar = c(1, 1, 1, 1))
        plot.new()
        legend("bottomleft", 
               legend = levels(gl(3, 1, labels = c("Non-infested", 
                                                   "Infested", 
                                                   "No information"))), 
               title = "Infestation with \nTriatoma infestans",
               box.lwd = 0, 
               fill = c("black", "springgreen", "gray"))
        dev.off()
        
        
################################################################################
### PART 5: NETWORK STRUCTURE AND SOCIO-DEMOGRAPHIC CHARACTERISTICS
################################################################################

### SETUP ###
    
# Clean environment
rm(list=ls()) 
graphics.off()

# Set directory
setwd("/Users/~/Desktop/ChagasNetworkPublic")

## Load data
load('data/CleanData.RData')

# get df 
df <- df_kin
    
### BARPLOT: Degree centrality vs. binary predictors ###
    
# Prep dataframe and variables
    
    ## subset binary variables
    df_bin <- subset(df, select = c(degree.cat, betweenness, harmonic, serology.occ, serology.child.occ, ethnicity, mobility, v.inf, iva.inf))
    
    ## add degree as binary variable
    df_bin$degree.bin <- factor(ifelse(df_bin$degree.cat == "0", "0", ">0"),
                                levels = c("0", ">0"))
    
    # for serology.child.occ, set observations with 'No child' to NA, and transform to binary
    str(df_bin$serology.child.occ)
    df_bin$serology.child.occ <- ifelse(df_bin$serology.child.occ == "No child", NA, 
                                        ifelse(df_bin$serology.child.occ == "Occurrence", 1, 0))
    
    ## transform remaining factor variables to binary 0/1 (numeric)
    library(dplyr)
    str(df_bin)
    df_bin <- df_bin %>% mutate(serology.occ = case_when(serology.occ == "Occurrence" ~ 1,
                                                         serology.occ == "No occurrence" ~ 0))
    df_bin <- df_bin %>% mutate(ethnicity = case_when(ethnicity == "Qom" ~ 1,
                                                      ethnicity == "Creole" ~ 0))
    df_bin <- df_bin %>% mutate(mobility = case_when(mobility == "mover" ~ 1,
                                                     mobility == "non-mover" ~ 0))
    df_bin <- df_bin %>% mutate(v.inf = case_when(v.inf == "infested" ~ 1,
                                                  v.inf == "non-infested" ~ 0))

# Create lists for looping over and saving results        
    
    ## list of variables to analyze
    variables <- list(
      serology.occ = "Occurrence of \nat least one \nseropositive person",
      serology.child.occ = "Occurrence of \nat least one \nseropositive child",
      ethnicity     = "Ethnicity",
      mobility      = "Mobility",
      v.inf         = "Infestation with \nT. infestans",
      iva.inf       = "Infestation w/ \nT.cruzi-infected \nT. infestans"
    )
    
    ## label replacements for binary categories
    category_labels <- list(
      serology.occ = c("No", "Yes"),
      serology.child.occ = c("No", "Yes"),
      ethnicity    = c("Creole", "Qom"),
      mobility     = c("Non-mover", "Mover"),
      v.inf        = c("Not infested", "Infested"),
      iva.inf      = c("Not infested", "Infested")
    )
    
# Barplot for degree.cat

    # create list for saving results
    results_list <- list()
    
    # Loop over each variable
    library(Hmisc)
    for (var in names(variables)) {
      
        # Frequency table
        tbl <- as.data.frame(table(df_bin$degree.cat, df_bin[[var]]))
        colnames(tbl) <- c("degree.cat", "category", "count")
        
        # Convert factors to characters to allow label assignment later
        tbl$category <- as.character(tbl$category)
        
        # Compute binomial confidence intervals for each (degree.cat, category)
        bin_list <- list()
        for (cat in unique(tbl$category)) {
          for (deg in levels(df_bin$degree.cat)) {
            x <- table(df_bin$degree.cat, df_bin[[var]])[deg, cat]
            n <- sum(table(df_bin$degree.cat, df_bin[[var]])[, cat])
            bin_list[[paste(deg, cat)]] <- binconf(x = x, n = n)
          }
        }
        
        # Combine binconf results
        bin_df <- do.call(rbind, bin_list)
        rownames(bin_df) <- NULL
        bin_df <- as.data.frame(bin_df)
        
        # Combine with count data
        df_out <- cbind(tbl, bin_df)
        colnames(df_out)[4:6] <- c("proportion", "lowerCI", "upperCI")
        
        # Replace category labels
        df_out$category <- recode(df_out$category,
                                  !!!setNames(category_labels[[var]], c("0", "1")))
        
        # Chi-squared test
        df_out$chisq <- NA
        df_out$chisq <- round(chisq.test(table(df_bin$degree.cat, df_bin[[var]]))$statistic, 3)
        df_out$df <- NA
        df_out$df <- chisq.test(table(df_bin$degree.cat, df_bin[[var]]))$parameter
        pval <- round(chisq.test(table(df_bin$degree.cat, df_bin[[var]]))$p.value, 3)
        df_out$p.value <- NA
        df_out$p.value <- pval
        df_out$pval <- NA
        df_out$pval[1] <- ifelse(pval < 0.001, "p<0.001",
                                 ifelse(pval < 0.01, "p<0.01",
                                        ifelse(pval < 0.05, "p<0.05", "ns")))
        
        # Add variable label
        df_out$variable <- variables[[var]]
        
        # Save to list
        results_list[[var]] <- df_out
    }
        
    # Combine all results and check structure
    df_bin_combined <- do.call(rbind, results_list)
    str(df_bin_combined)
        
    # Create factor levels for category (this will ensure that barplots appear in the correct order)
    df_bin_combined$category <- factor(df_bin_combined$category, levels = unique(df_bin_combined$category))
    df_bin_combined$category <- with(df_bin_combined, factor(category, levels = c(
      "No", "Yes",               # serology.occ & serology.child.occ
      "Creole", "Qom",           # ethnicity
      "Non-mover", "Mover",      # mobility
      "Not infested", "Infested" # v.inf and iva.inf
    )))       
        
    # Create barplot
        
        ## plot
        library(ggplot2)
        p <- ggplot(df_bin_combined) +
          geom_bar(aes(x = category, y = proportion, fill = forcats::fct_rev(degree.cat)),
                   position = "stack", stat = "identity",
                   width = 0.6, alpha = 0.7, color = "black", size = 0.1) +
          geom_text(
            aes(x = category, y = proportion, 
                label = scales::percent(proportion, accuracy = 1),
                group = forcats::fct_rev(degree.cat)),
            position = position_stack(vjust = 0.75),
            color = "black", size = 3) +
          facet_grid(.~variable, scales = "free_x") +
          scale_fill_manual(values = c("#E18F5C", "#F6CD73", "#E6E796")) +
          labs(x = "", y = "Proportion",
               fill = "Degrees \nof kinship") +
          theme_minimal() +
          theme(strip.text.x = element_text(size = 10, face = "bold"),
                # strip.text.y = element_text(size = 10, face = "bold"))
                plot.title = element_text(size = 12, vjust=6),
                plot.margin = margin(t = 25, r = 10, b = 10, l = 10),
                panel.grid.major.x = element_blank(),
                axis.title = element_text(size=10),
                axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1.3),
                axis.text.y = element_text(size = 10),
                legend.title=element_text(size=10))
        
        ## create labels for p-value
        df_labels <- data.frame(
          label = c(" p<0.001", "      ns", "      ns", "      ns", "      ns", "  p<0.05"), #have to pad labels w/ space so that they are equal length
          variable   = c("Ethnicity", "Infestation w/ \nT.cruzi-infected \nT. infestans", "Infestation with \nT. infestans", "Mobility", "Occurrence of \nat least one \nseropositive person", "Occurrence of \nat least one \nseropositive child"),
          x     = c("Creole", "Not infested", "Not infested", "Non-mover", "No", "No"),
          y     = c(1, 1, 1, 1, 1, 1)
        )
        
        ## replot w/ additional p-val labels overlaid
        barplot_binvars_degree <- p +
          geom_text(data = df_labels, size = 3.5,           
                    mapping = aes(x = x, y = y, label = label,
                                  hjust = 0, vjust = -0.7))
        
        ## save plot
        png("figures/network/bivariate_descriptive/barplot_binvars_degree.cat.png",width=10,height=6,units="in",res=600)
        barplot_binvars_degree
        dev.off()

# Barplot for degree.bin
        
    # Loop over each variable
    library(Hmisc)
    for (var in names(variables)) {
      
        # frequency table
        tbl <- as.data.frame(table(df_bin$degree.bin, df_bin[[var]]))
        colnames(tbl) <- c("degree.bin", "category", "count")
        
        # convert factors to characters to allow label assignment later
        tbl$category <- as.character(tbl$category)
        
        # compute binomial confidence intervals for each (degree.bin, category)
        bin_list <- list()
        for (cat in unique(tbl$category)) {
          for (deg in levels(df_bin$degree.bin)) {
            x <- table(df_bin$degree.bin, df_bin[[var]])[deg, cat]
            n <- sum(table(df_bin$degree.bin, df_bin[[var]])[, cat])
            bin_list[[paste(deg, cat)]] <- binconf(x = x, n = n)
          }
        }
        
        # combine binconf results
        bin_df <- do.call(rbind, bin_list)
        rownames(bin_df) <- NULL
        bin_df <- as.data.frame(bin_df)
        
        # combine with count data
        df_out <- cbind(tbl, bin_df)
        colnames(df_out)[4:6] <- c("proportion", "lowerCI", "upperCI")
        
        # replace category labels
        df_out$category <- recode(df_out$category,
                                  !!!setNames(category_labels[[var]], c("0", "1")))
        
        # Chi-squared test: save statistic, p-value, and significance level
        df_out$chisq <- NA
        df_out$chisq <- round(chisq.test(table(df_bin$degree.bin, df_bin[[var]]))$statistic, 3)
        df_out$df <- NA
        df_out$df <- chisq.test(table(df_bin$degree.cat, df_bin[[var]]))$parameter
        pval <- round(chisq.test(table(df_bin$degree.bin, df_bin[[var]]))$p.value, 3)
        df_out$p.value <- NA
        df_out$p.value <- pval
        df_out$pval <- NA
        df_out$pval[1] <- ifelse(pval < 0.001, "p<0.001",
                                 ifelse(pval < 0.01, "p<0.01",
                                        ifelse(pval < 0.05, "p<0.05", "ns")))

        # Add variable label
        df_out$variable <- variables[[var]]
        
        # Save to list
        results_list[[var]] <- df_out
    }
    
    # Combine all results and check structure
    df_bin_combined <- do.call(rbind, results_list)
    str(df_bin_combined)
    
    # Create factor levels for category (this will ensure that barplots appear in the correct order)
    df_bin_combined$category <- factor(df_bin_combined$category, levels = unique(df_bin_combined$category))
    df_bin_combined$category <- with(df_bin_combined, factor(category, levels = c(
      "No", "Yes",               # serology.occ & serology.child.occ
      "Creole", "Qom",           # ethnicity
      "Non-mover", "Mover",      # mobility
      "Not infested", "Infested" # v.inf and iva.inf
    )))       
    
    # Create barplot
    
        ## plot
        library(ggplot2)
        p <- ggplot(df_bin_combined) +
          geom_bar(aes(x = category, y = proportion, fill = forcats::fct_rev(degree.bin)),
                   position = "stack", stat = "identity",
                   width = 0.6, alpha = 0.7, color = "black", size = 0.1) +
          geom_text(
            aes(x = category, y = proportion, 
                label = scales::percent(proportion, accuracy = 1),
                group = forcats::fct_rev(degree.bin)),
            position = position_stack(vjust = 0.75),
            color = "black", size = 3) +
          facet_grid(.~variable, scales = "free_x") +
          scale_fill_manual(values = c("#E18F5C", "#E6E796")) +
          labs(x = "", y = "Proportion",
               fill = "Degrees \nof kinship") +
          theme_minimal() +
          theme(strip.text.x = element_text(size = 10, face = "bold"),
                # strip.text.y = element_text(size = 10, face = "bold"))
                plot.title = element_text(size = 12, vjust=6),
                plot.margin = margin(t = 25, r = 10, b = 10, l = 10),
                panel.grid.major.x = element_blank(),
                axis.title = element_text(size=10),
                axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1.3),
                axis.text.y = element_text(size = 10),
                legend.title=element_text(size=10))
        
        ## create labels for p-value
        df_labels <- data.frame(
          label = c(" p<0.001", "      ns", "      ns", "      ns", "      ns", "      ns"), #have to pad labels w/ space so that they are equal length
          variable   = c("Ethnicity", "Infestation w/ \nT.cruzi-infected \nT. infestans", "Infestation with \nT. infestans", "Mobility", "Occurrence of \nat least one \nseropositive person", "Occurrence of \nat least one \nseropositive child"),
          x     = c("Creole", "Not infested", "Not infested", "Non-mover", "No", "No"),
          y     = c(1, 1, 1, 1, 1, 1)
        )
        
        ## replot w/ additional p-val labels overlaid
        barplot_binvars_degree <- p +
          geom_text(data = df_labels, size = 3.5,           
                    mapping = aes(x = x, y = y, label = label,
                                  hjust = 0, vjust = -0.7))
        
        ## save plot
        png("figures/network/bivariate_descriptive/barplot_binvars_degree.bin.png",width=10,height=6,units="in",res=600)
        barplot_binvars_degree
        dev.off()
        
### BOXPLOT: BETWEENNESS & HARMONIC CENTRALITY ###
    
# Subset binary variables
df_bin <- subset(df, select = c(degree.cat, betweenness, harmonic, serology.occ, serology.child.occ, ethnicity, mobility, v.inf, iva.inf))
str(df_bin) 

# Check that all categorical variables are binary

    ## reclassify serology.child.occ: set observations with 'No child' to NA
    df_bin$serology.child.occ <- ifelse(df_bin$serology.child.occ == "No child", NA, 
                                        ifelse(df_bin$serology.child.occ == "Occurrence", 1, 0))
    
    ## convert serology.child.occ to factor
    df_bin$serology.child.occ <- factor(df_bin$serology.child.occ,
                                        levels = c(0, 1),
                                        labels = c("No occurrence", "Occurrence"))
    
    ## reclassify iva.inf: set observations that are non-infested households to NA
    df_bin$iva.inf <- ifelse(df_bin$v.inf == "non-infested", NA,
                             ifelse(df_bin$iva.inf == 1, "Occurrence", "No occurrence"))
    
    ## convert iva.inf to factor
    df_bin$iva.inf <- factor(df_bin$serology.child.occ,
                             levels = c("No occurrence", "Occurrence"),
                             labels = c("No occurrence", "Occurrence"))
    
    ## Rename factor labels
    df_bin$mobility <- factor(df_bin$mobility,
                              levels = c("non-mover", "mover"),
                              labels = c("Non-mover", "Mover"))
    df_bin$v.inf <- factor(df_bin$v.inf,
                           levels = c("non-infested", "infested"),
                           labels = c("No occurrence", "Occurrence"))
     
    ## check
    str(df_bin)
    
# Wilcoxon Rank-Sum test (aka Mann-Whitney U test)

    ## Note: use when (a) comparing exactly two independent groups and (b) if your dependent variable is not normally distributed (non-parametric alternative to a two-sample test)

    ## Dependent variables to test
    dep_vars <- c("betweenness", "harmonic")
    
    ## Binary variables to test
    bin_vars <- c("serology.occ", "serology.child.occ", "ethnicity", "mobility", "v.inf", "iva.inf")
      
    ## create empty df for saving results 
    df_wilcox <- data.frame(dep_var = character(),
                            bin_var = character(),
                            W = numeric(),
                            p.val = numeric(),
                            significance = character(),
                            stringsAsFactors = FALSE  # optional, to avoid factor conversion
                 )
        
    ## loop through variables
    for(dep_var in dep_vars) {
      
      for (bin_var in bin_vars) {
      
      ## create formula dynamically
      formula <- as.formula(paste(dep_var, "~", bin_var))
      
      ## wilcox.test
      test <- wilcox.test(formula, data = df_bin)
      
      ## create new row of data
      new_row <- data.frame(
        dep_var = dep_var,
        bin_var = bin_var,
        W = test$statistic,
        p.val = test$p.value,
        significance = ifelse(test$p.value < 0.001, "p<0.001",
                             ifelse(test$p.value  < 0.01, "p<0.01",
                                    ifelse(test$p.value  < 0.05, "p<0.05", "ns")))
      )
      
      ## add to df_wilcox
      df_wilcox <- rbind(df_wilcox, new_row)
      }
    }
    

# Betweenness boxplots
    
    ## create df w/ betweenness keeping only non-zero betweenness values
    df_between <- subset(df_bin, select = -c(degree.cat, harmonic))
    df_between <- df_between[df_between$betweenness != 0,]

    ## transform df to long format w/ pivot_longer
        
        # add id column if not already present
        df_between <- df_between %>%
          dplyr::mutate(id = row_number())
        
        # pivot to long format
        df_long <- df_between %>%
                    tidyr::pivot_longer(cols = c(serology.occ, serology.child.occ, ethnicity, 
                                          mobility, v.inf, iva.inf),
                                 names_to = "variable",
                                 values_to = "group") %>%
                    dplyr::filter(!is.na(group))  # Optional: remove rows with missing group values
    
    ## create factor levels for category (this will ensure that barplots appear in the correct order)
    df_long$variable <- factor(df_long$variable, levels = unique(df_long$variable))
    df_long$variable <- factor(df_long$variable, 
                               levels = c("serology.child.occ", "serology.occ", "iva.inf", "v.inf", "mobility", "ethnicity"), 
                               labels = c("At least one \nseropositive child", 
                                          "At least one \nseropositive inhabitant", 
                                          "Infestation with \nT. cruzi-infested T. infestans \n(among infested households)", 
                                          "Infestation with \nT. infestans", 
                                          "Mobility: Mover \n(vs. non-mover)", 
                                          "Ethnicity: Qom \n(vs. Creole)")
    )
    
    ## grouped boxplot
    png("figures/network/bivariate_descriptive/boxplot_binvars_betweenness.png", width=9,height=6,units="in",res=600)
    ggplot(df_long, aes(x=variable, y=betweenness, fill=factor(group))) +
           geom_boxplot(width = 0.5) +
           labs(fill = "Group") +
           labs(x = "", y = "Betwenness centrality") +
           scale_y_continuous(breaks = seq(0, 140, by = 20)) +
           scale_fill_manual(values = c(
             "Creole" = alpha("#B0DEA5", 0.6),  
             "Qom" = alpha("#556B2F", 0.8), 
             "Non-mover" = alpha("#B0DEA5", 0.6),  
             "Mover" = alpha("#556B2F", 0.8), 
             "No occurrence" = alpha("#B0DEA5", 0.6),  
             "Occurrence" = alpha("#556B2F", 0.8)     
           )) +
           coord_flip() +
           theme_minimal()
    dev.off()
    
## Harmonic boxplots
    
    ## create df w/ betweenness keeping only non-zero betweenness values
    df_harmonic <- subset(df_bin, select = -c(degree.cat, betweenness))
    df_harmonic <- df_harmonic[df_harmonic$harmonic != 0,]
    
    ## transform df to long format w/ pivot_longer
    
        ## add id column if not already present
        df_harmonic <- df_harmonic %>%
        dplyr::mutate(id = row_number())
    
        ## pivot to long format
        df_long <- df_harmonic %>%
        tidyr::pivot_longer(cols = c(serology.occ, serology.child.occ, ethnicity, 
                                     mobility, v.inf, iva.inf),
                            names_to = "variable",
                            values_to = "group") %>%
        dplyr::filter(!is.na(group))  # Optional: remove rows with missing group values
    
    ## create factor levels for category (this will ensure that barplots appear in the correct order)
    df_long$variable <- factor(df_long$variable, levels = unique(df_long$variable))
    df_long$variable <- factor(df_long$variable, 
                               levels = c("serology.child.occ", "serology.occ", "iva.inf", "v.inf", "mobility", "ethnicity"), 
                               labels = c("At least one \nseropositive child", 
                                          "At least one \nseropositive inhabitant", 
                                          "Infestation with \nT. cruzi-infested T. infestans \n(among infested households)", 
                                          "Infestation with \nT. infestans", 
                                          "Mobility: Mover \n(vs. non-mover)", 
                                          "Ethnicity: Qom \n(vs. Creole)")
                              )
    
    ## grouped boxplot
    png("figures/network/bivariate_descriptive/boxplot_binvars_harmonic.png", width=9,height=6,units="in",res=600)
    ggplot(df_long, aes(x=variable, y=harmonic, fill=factor(group))) +
      geom_boxplot(width = 0.5) +
      labs(fill = "Group") +
      labs(x = "", y = "Harmonic centrality") +
      # scale_y_continuous(breaks = seq(0, 10, by = 2)) +
      scale_fill_manual(values = c(
        "Creole" = alpha("#BCE4D8", 0.6),  
        "Qom" = alpha("#4682B4", 0.8), 
        "Non-mover" = alpha("#BCE4D8", 0.6),  
        "Mover" = alpha("#4682B4", 0.8), 
        "No occurrence" = alpha("#BCE4D8", 0.6),  
        "Occurrence" = alpha("#4682B4", 0.8)     
      )) +
      coord_flip() +
      theme_minimal()
    dev.off()
        

### VIOLIN/BOXPLOTS for continuous and discrete variables ###

# get df 
df <- df_kin

# Significance tests
    
    # Kruskal-Wallis (non-parametric version of one-way ANOVA)
        
        ## Note: (a) comparing three or more independent groups, (b) your dependent variable is not normally distributed
        
        ## List of variable names to test
        vars_to_test <- c("svi", "kin.svi.med", "kin.svi.relmed", "kin.svi.relsum",
                          "hab_tot", "hab_0a14", "age", 
                          "iva", "kin.inf.prop", "kin.iva.med", "kin.iva.sum",
                          "serology", "serology.child", "kin.sero.prop", "kin.sero.med", "kin.sero.sum")
        
        ## Create an empty list to store results
        kruskal_results <- list()
        
        ## Loop over each variable and run Kruskal-Wallis test
        for (var in vars_to_test) {
          
          # Build formula dynamically
          formula <- as.formula(paste(var, "~ degree.cat"))
          
          # Run Kruskal-Wallis test
          test <- kruskal.test(formula, data = df)
          
          # Store results in a list
          kruskal_results[[var]] <- list(
            variable = var,
            statistic = test$statistic,
            p_value = test$p.value,
            df = test$parameter)
          
        }
        
        ## Convert list to data frame
        kruskal_df <- do.call(rbind, lapply(kruskal_results, as.data.frame))
        rownames(kruskal_df) <- NULL
        
        ## note significance
        kruskal_df$significance <- ifelse(kruskal_df$p_value < 0.001, "p<0.001",
                                          ifelse(kruskal_df$p_value < 0.01, "p<0.01",
                                                 ifelse(kruskal_df$p_value < 0.05, "p<0.05", "ns")))
        
        ## View results
        print(kruskal_df)

    # Wilcoxon rank-sum (Mann-Whitney U; non-parametric version of t-test)
        
        ## Note: (a) comparing two independent groups, (b) your dependent variable is not normally distributed
        
        ## List of variable names to test
        vars_to_test <- c("kin.svi.relmed", "kin.svi.relsum",
                          "kin.iva.med", "kin.iva.sum",
                          "kin.sero.med", "kin.sero.sum")
        
        ## Create df excluding houses w/ no kin connections
        df_nonzero <- df[df$degree>0,]
        
        ## Create result list
        wilcox_results <- list()
        
        # Loop over each variable and run Mann-Whitney U test
        for (var in vars_to_test) {
          
          # Extract the two groups
          group1 <- df[df$degree.cat == levels(df$degree.cat)[2], var] #"1-2"
          group2 <- df[df$degree.cat == levels(df$degree.cat)[3], var] #"3+"
          
          # Run Wilcoxon (Mann-Whitney U) test
          test <- wilcox.test(group1, group2, exact = FALSE)
          
          # Store results in a list
          wilcox_results[[var]] <- list(
            variable = var,
            statistic = test$statistic,
            p_value = test$p.value)
        }
        
        ## Convert list to data frame
        wilcox_df <- do.call(rbind, lapply(wilcox_results, as.data.frame))
        rownames(wilcox_df) <- NULL
        
        ## note significance
        wilcox_df$significance <- ifelse(wilcox_df$p_value < 0.001, "p<0.001",
                                          ifelse(wilcox_df$p_value < 0.01, "p<0.01",
                                                 ifelse(wilcox_df$p_value < 0.05, "p<0.05", "ns")))
        
        ## View results
        print(wilcox_df)
    
# Prep for violin-boxplots
    
    ## load libraries
    library(ggplot2)
    library(viridis)
    
    ## set color palette
    mycolors <- viridis(3)
    
    ## create plot theme
    mytheme <- theme(axis.text=element_text(size=8),
                     axis.title=element_text(size=10),
                     plot.title = element_text(size=10))
    
# Age and number of inhabitants by age group
    
    ## hab_tot vs. degree.cat
    length(which(!is.na(df$hab_tot)))
    vp.hab_tot <- ggplot(df, aes(x = degree.cat, y = hab_tot)) +
      geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
      geom_boxplot(width = 0.1, fill = "white") + 
      scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
      labs(title = "Number of inhabitants \nper household (n = 387)", 
           x = "Kinship degree", y = "Number of inhabitants per household",
           fill = "Degrees \nof kinship") +
      theme_minimal()
    
    ## hab_0a14 vs. degree.cat
    length(which(!is.na(df$hab_0a14)))
    vp.hab_0a14 <- ggplot(df, aes(x = degree.cat, y = hab_0a14)) +
      geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
      geom_boxplot(width = 0.1, fill = "white") + 
      scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
      labs(title = "Number of children (<15yo) \nper household (n = 386)", 
           x = "Kinship degree", y = "Number of children per household",
           fill = "Degrees \nof kinship") +
      theme_minimal()
    
    ## age vs. degree.cat
    length(which(!is.na(df$age)))
    vp.age <- ggplot(df, aes(x = degree.cat, y = age)) +
      geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
      geom_boxplot(width = 0.1, fill = "white") + 
      scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
      labs(title = "Mean age in years \nper household (n = 350)", 
           x = "Kinship degree", y = "Mean age in years",
           fill = "Degrees \nof kinship") +
      theme_minimal()
    
    ## get legend
    library(cowplot)
    legend <- get_legend(
      # create some space to the left of the legend
      vp.hab_tot + theme(legend.box.margin = margin(0, 0, 0, 12))
    )
    
    ## combine plots and save
    png("figures/network/bivariate_descriptive/violin_age.png", width=9,height=9,units="in",res=600)
    plot_grid(vp.hab_tot + theme(legend.position="none") + mytheme,
              vp.hab_0a14 + theme(legend.position="none") + mytheme,
              vp.age + theme(legend.position="none") + mytheme,
              legend,
              nrow = 2, ncol = 2
    )
    dev.off()
        
    
# Social vulnerability variables
    
    ## svi vs. degree.cat
    length(which(!is.na(df$svi))) # n=386
    summary(df$svi)
    vp_svi <- ggplot(df, aes(x = degree.cat, y = svi)) +
      geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
      geom_boxplot(width = 0.1, fill = "white") + 
      scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
      ylim(-4,4) +
      labs(title = "Social vulnerability index (n = 386)", 
           x = "Kinship degree", y = "Social vulnerability index",
           fill = "Degrees \nof kinship") +
      theme_minimal()
    
    ## kin.svi.med vs. degree.cat
    length(which(!is.na(df$kin.svi.med))) # n=388
    summary(df$kin.svi.med)
    vp_kin.svi.med <- ggplot(df, aes(x = degree.cat, y = kin.svi.med)) +
      geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
      geom_boxplot(width = 0.1, fill = "white") + 
      scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
      ylim(-4,4) +
      labs(title = "Median social vulnerability index \namong connected houses (n = 388)", 
           x = "Kinship degree", y = "Median social vulnerability index",
           fill = "Degrees \nof kinship") +
      theme_minimal()
    
    ## kin.svi.relmed vs. degree.cat
    length(which(!is.na(df$kin.svi.relmed))) # n=402
    summary(df$kin.svi.relmed)
    vp_kin.svi.relmed <- ggplot(df, aes(x = degree.cat, y = kin.svi.relmed)) +
      geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
      geom_boxplot(width = 0.1, fill = "white") + 
      scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
      ylim(-4,4) +
      labs(title = "Median relative social vulnerability index \namong connected houses (n = 402)", 
           x = "Kinship degree", y = "Median relative social vulnerability index",
           fill = "Degrees \nof kinship") +
      theme_minimal()
    
    ## kin.svi.relsum vs. degree.cat
    length(which(!is.na(df$kin.svi.relsum))) #n=402
    summary(df$kin.svi.relsum[df$degree==0])
    vp_kin.svi.relsum <- ggplot(df, aes(x = degree.cat, y = kin.svi.relsum)) +
      geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
      geom_boxplot(width = 0.1, fill = "white") + 
      scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
      ylim(-6,8) +
      labs(title = "Cumulative relative social vulnerability index \namong connected houses (n = 402)", 
           x = "Kinship degree", y = "Cumulative relative social vulnerability index",
           fill = "Degrees \nof kinship") +
      theme_minimal()
    
    ## get legend
    library(cowplot)
    legend <- get_legend(
      # create some space to the left of the legend
      vp_svi + theme(legend.box.margin = margin(0, 0, 0, 12))
    )
    
    ## combine plots and save
    png("figures/network/bivariate_descriptive/violin_svi.png", width=9,height=9,units="in",res=600)
    plot_grid(vp_svi + theme(legend.position = "none") + mytheme,
              vp_kin.svi.med + theme(legend.position="none") + mytheme,
              vp_kin.svi.relmed + theme(legend.position="none") + mytheme,
              vp_kin.svi.relsum + theme(legend.position="none") + mytheme,
              legend,
              nrow = 3, ncol = 2
    )
    dev.off()
    
    ## Histograms & density plots of SVI of connected vs. unconnected houses
    
        # # svi: overlaid histograms & density plots
        # table(df$connected)
        # hist_svi <- ggplot(df, aes(x=svi, fill=connected)) +
        #   geom_density(aes(group=df$connected, alpha=.1, fill=df$connected)) +
        #   geom_histogram(position="identity", alpha=0.4, aes(y=..density..)) +
        #   ylim(0,1) +
        #   labs(title = "Social vulnerability among kin connected (n=246) \nvs. unconnected houses (n=57)", 
        #        x = "Social vulnerablity index", y = "Density",
        #        fill = "Connectivity") +
        #   theme_minimal()
        # 
        # # kin.svi.med: overlaid histograms & density plots
        # hist_kin.svi.med <- ggplot(df, aes(x=kin.svi.med, fill=connected)) +
        #   geom_density(aes(group=df$connected, alpha=.1, fill=df$connected)) +
        #   geom_histogram(position="identity", alpha=0.4, aes(y=..density..)) +
        #   ylim(0,1) +
        #   labs(title = "Median social vulnerability among houses \nconnected to each kin-connected focal house (n=246) \nvs. the social vulnerability of unconnected houses (n=57)", 
        #        x = "(Median) Social vulnerability index", y = "Density",
        #        fill = "Connectivity") +
        #   theme_minimal()
        # 
        # 
        # ## get legend
        # legend <- get_legend(
        #   # create some space to the left of the legend
        #   hist_svi + theme(legend.box.margin = margin(0, 0, 0, 12))
        # )
        # 
        # # combine plots
        # hist_svi <- plot_grid(hist_svi + theme(legend.position="none") + mytheme,
        #                       hist_kin.svi.med + theme(legend.position="none") + mytheme,
        #                       legend,
        #                       nrow = 1, ncol = 3)
        # 
        # ## save plots
        # png("figures/network/bivariate_descriptive/histogram_svi.png",width=12,height=6,units="in",res=600)
        # hist_svi
        # dev.off()

    # Scatterplots (Loess method) of SVI-related vars vs. network centrality measures among connected houses
    
        ## subset unconnected houses and connected houses
        df_conn <- df[df$connected == "Connected",]
        
        ## degree
            
            # svi
            scatter_degree_svi <- ggplot(df_conn, aes(x = degree, y = svi)) +
              geom_point() +
              geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
              theme_minimal() +
              labs(
                   # title = "Degree vs. Social vulnerability", 
                   x = "Degree centrality", y = "Social vulnerability index")
            
            # kin.svi.med
            scatter_degree_kinsvi <- ggplot(df_conn, aes(x=degree, y=kin.svi.med)) +
              geom_point() +
              geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
              theme_minimal() +
              labs(
                   # title = "Degree vs. Social vulnerability \nof kin network", 
                   x = "Degree centrality", y = "Median social vulnerability index")
            
            # kin.svi.relmed
            scatter_degree_kinsvirel <- ggplot(df_conn, aes(x=degree, y=kin.svi.relmed)) +
              geom_point() +
              geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
              theme_minimal() + 
              labs(
                   # title = "Degree vs. Relative social vulnerability of kin network",
                   x = "Degree centrality", y = "Median social vulnerability index")
              

        ## betweenness

            # svi
            scatter_between_svi <- ggplot(df_conn, aes(x = betweenness, y = svi)) +
              geom_point() +
              geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
              theme_minimal() +
              labs(
                # title = "Betweenness vs. Social vulnerability", 
                x = "Betweenness centrality", y = "Social vulnerability index")
        
            # kin.svi.med
            scatter_between_kinsvi <- ggplot(df_conn, aes(x=betweenness, y=kin.svi.med)) +
              geom_point() +
              geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
              theme_minimal() +
              labs(
                # title = "Betweenness vs. Social vulnerability \nof kin network", 
                x = "Betweenness centrality", y = "Median social vulnerability index")
            
            # kin.svi.relmed
            scatter_between_kinsvirel <- ggplot(df_conn, aes(x=betweenness, y=kin.svi.relmed)) +
              geom_point() +
              geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
              theme_minimal() + 
              labs(
                # title = "Betweenness vs. Relative social vulnerability of kin network",
                x = "Betweenness centrality", y = "Median social vulnerability index")
            
        ## harmonic centrality
        
            # svi
            scatter_harmonic_svi <- ggplot(df_conn, aes(x = harmonic, y = svi)) +
              geom_point() +
              geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
              theme_minimal() +
              labs(
                # title = "Harmonic centrality vs. Social vulnerability", 
                x = "Harmonic centrality", y = "Social vulnerability index")
            
            # kin.svi.med
            scatter_harmonic_kinsvi <- ggplot(df_conn, aes(x=harmonic, y=kin.svi.med)) +
              geom_point() +
              geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
              theme_minimal() +
              labs(
                # title = "Harmonic centrality vs. Social vulnerability \nof kin network", 
                x = "Harmonic centrality", y = "Median social vulnerability index")
            
            # kin.svi.relmed
            scatter_harmonic_kinsvirel <- ggplot(df_conn, aes(x=harmonic, y=kin.svi.relmed)) +
              geom_point() +
              geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
              theme_minimal() + 
              labs(
                # title = "Harmonic centrality vs. Relative social vulnerability of kin network",
                x = "Harmonic centrality", y = "Median social vulnerability index")
        
        # combine plots (adding space in between each plot)
        library(cowplot)
        scatter_networkSVI <- plot_grid(scatter_degree_svi + mytheme,
                                        scatter_degree_kinsvi + mytheme,
                                        scatter_degree_kinsvirel + mytheme,
                                        scatter_between_svi + mytheme,
                                        scatter_between_kinsvi + mytheme,
                                        scatter_between_kinsvirel + mytheme,
                                        scatter_harmonic_svi + mytheme,
                                        scatter_harmonic_kinsvi + mytheme,
                                        scatter_harmonic_kinsvirel + mytheme,
                                        nrow = 3, ncol = 3)
        
        ## save plots
        png("figures/network/bivariate_descriptive/scatter_networkSVI.png",width=12,height=12,units="in",res=600)
        scatter_networkSVI
        dev.off()

        ## significance testing w/ Spearman's Rank Correlation Coefficient (for non-normally distributed continuous variables)
            
            # # create empty df
            # spearman_df <- data.frame(variable1 = character(),
            #                       variable2 = character(),
            #                       correlation = numeric(),
            #                       p_value = numeric(),
            #                       stringsAsFactors = FALSE
            #             )
            
            # List of variable names to test
            vars_to_test <- c("svi", "kin.svi.med", "kin.svi.relmed")
            
            ## Create result list
            spearman_results <- list()
            
            # DEGREE: Loop over each variable and run Spearman's rank correlation test
            for (var in vars_to_test) {
                
                ## create df w/ no NAs
                df_clean <- df_conn[!is.na(df_conn$degree) & !is.na(df_conn[[var]]), ]
                
                ## run test
                test <- cor.test(df_clean$degree, df_clean[[var]], method = "spearman")
                
                # store results in a list
                spearman_results[[var]] <- list(
                  var1 = "degree",
                  var2 = var,
                  correlation = test$estimate,
                  p_value = test$p.value)
            }
            
            ## convert list to data frame
            spearman_df <- do.call(rbind, lapply(spearman_results, as.data.frame))

            # BETWEENNESS: Loop over each variable and run Spearman's rank correlation test
            for (var in vars_to_test) {
              
              ## create df w/ no NAs
              df_clean <- df_conn[!is.na(df_conn$betweenness) & !is.na(df_conn[[var]]), ]
              
              ## run test
              test <- cor.test(df_clean$betweenness, df_clean[[var]], method = "spearman")
              
              # store results in a list
              spearman_results[[var]] <- list(
                var1 = "betweenness",
                var2 = var,
                correlation = test$estimate,
                p_value = test$p.value)
            }
            
            ## convert list to data frame
            spearman_df <- rbind(spearman_df, do.call(rbind, lapply(spearman_results, as.data.frame)))

            # HARMONIC: Loop over each variable and run Spearman's rank correlation test
            for (var in vars_to_test) {
              
              ## create df w/ no NAs
              df_clean <- df_conn[!is.na(df_conn$harmonic) & !is.na(df_conn[[var]]), ]
              
              ## run test
              test <- cor.test(df_clean$harmonic, df_clean[[var]], method = "spearman")
              
              # store results in a list
              spearman_results[[var]] <- list(
                var1 = "harmonic",
                var2 = var,
                correlation = test$estimate,
                p_value = test$p.value)
            }
            
            ## convert list to data frame
            spearman_df <- rbind(spearman_df, do.call(rbind, lapply(spearman_results, as.data.frame)))
            rownames(spearman_df) <- NULL
            
            ## note significance
            spearman_df$significance <- ifelse(spearman_df$p_value < 0.001, "p<0.001",
                                               ifelse(spearman_df$p_value < 0.01, "p<0.01",
                                                      ifelse(spearman_df$p_value < 0.05, "p<0.05", "ns")))
            
            ## View results
            print(spearman_df)

# Host Availability Index
            
    ## hai vs. degree.cat
    length(which(!is.na(df$hai))) # n=388
    summary(df$hai)
    vp_hai <- ggplot(df, aes(x = degree.cat, y = hai)) +
      geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
      geom_boxplot(width = 0.1, fill = "white") + 
      scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
      ylim(-4,4) +
      labs(title = "Host availability index (n = 388)", 
           x = "Kinship degree", y = "Host availability index",
           fill = "Degrees \nof kinship") +
      theme_minimal()
        
        # save plot
        png("figures/network/bivariate_descriptive/violin_hai.png",width=8,height=4,units="in",res=600)
        vp_hai
        dev.off()
        
        # kruskal-wallis test
        kruskal.test(df$hai, df$degree.cat)
            
    # Scatterplots (Loess method) of hai and centrality measures
            
        ## subset unconnected houses and connected houses
        df_conn <- df[df$connected == "Connected",]
        
        ## degree
        scatter_degree_hai <- ggplot(df_conn, aes(x = degree, y = hai)) +
          geom_point() +
          geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
          theme_minimal() +
          labs(
            # title = "Degree vs. Host availability", 
            x = "Degree centrality", y = "Host availability index")
        
        ## betweenness
        scatter_between_hai <- ggplot(df_conn, aes(x = betweenness, y = hai)) +
          geom_point() +
          geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
          theme_minimal() +
          labs(
            # title = "Betweenness vs. Host availability", 
            x = "Betweenness centrality", y = "Host availability index")
        
        ## harmonic
        scatter_harmonic_hai <- ggplot(df_conn, aes(x = harmonic, y = hai)) +
          geom_point() +
          geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
          theme_minimal() +
          labs(
            # title = "Harmonic centrality vs. Host availability", 
            x = "Harmonic centrality", y = "Host availability index")
        
        ## combine plots (adding space in between each plot)
        library(cowplot)
        scatter_networkHAI <- plot_grid(scatter_degree_hai + mytheme,
                                        scatter_between_hai + mytheme,
                                        scatter_harmonic_hai + mytheme,
                                        nrow = 1, ncol = 3)
            ## save plots
            png("figures/network/bivariate_descriptive/scatter_networkHAI.png",width=12,height=4,units="in",res=600)
            scatter_networkHAI
            dev.off()
            
        ## Spearman's rank correlation test
        cor.test(df_conn$degree, df_conn$hai)
        cor.test(df_conn$harmonic, df_conn$hai)
        cor.test(df_conn$betweenness, df_conn$hai)
        cor.test(df$degree, df$hai)
        cor.test(df$harmonic, df$hai)
        cor.test(df$betweenness, df$hai)
            
        
# Vector infestation

    ## focal house infestation status vs. degree.cat
    
        ## create labels for categories
        label_infected <- "Infested \nwith T. cruzi-infected T. infestans"
        label_infested <- "Infested with T. infestans"
        label_none <- "Non-infested"
        
        ## create cats
        df$focal.v <- NA
        df$focal.v <- factor(ifelse(df$iva > 0 & df$v.inf == "infested", label_infected,
                                    ifelse(df$v.inf == "infested", label_infested,
                                           ifelse(df$v.inf == "non-infested", label_none, df$focal.v))),
                           levels = c(label_none, label_infested, label_infected)
        )
        
        # ## for degree, recat households with no kin ties as NA
        # df$degree.nonzero <- factor(ifelse(df$degree.cat == "0", NA, as.character(df$degree.cat)),
        #                             levels = c("1-2", "3+")
        # )
        
        ## chi-square test excluding households w/ no kin connections
        table(df$degree.cat, df$focal.v)
        chisq_focal.v <- chisq.test(table(df$degree.cat, df$focal.v))
        
        # barplot
        library(ggplot2)
        bp_focal.v <- df %>%
            filter(!is.na(degree.cat) & !is.na(focal.v)) %>%
            ggplot(aes(x = degree.cat, fill = focal.v)) +
            # geom_bar(position = "fill", alpha=0.7, width = 0.5) +  # stack by proportion
            geom_bar(position = "fill", color = "black", width = 0.5) +  # stack by proportion
            scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
            # scale_fill_manual(values = c("#440154FF", "#21908CFF","#FDE725FF")) +
            scale_fill_manual(values = c("#41448733", "#41448799","#414487FF")) +
            labs(
              # title = "Proportion of houses by infestation status",
              x = "Kinship degree",
              y = "Proportion of houses",
              fill = "Infestation status"
            ) +
            theme_minimal()
        
        ## iva vs. degree.cat
        df_nozero <- df[df$iva>0,]
        length(which(!is.na(df_nozero$iva))) # n=25
        vp_iva <- ggplot(df_nozero, aes(x = degree.cat, y = iva)) +
            geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
            geom_boxplot(width = 0.1, fill = "white") + 
            scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
            labs(
                 # title = "Relative abundance of infected vector among \nhouses with at least one infected vector (n = 25)", 
                 x = "Kinship degree", y = "Abundance of infected vector",
                 fill = "Degrees \nof kinship") +
            theme_minimal()
        
        ## get legend
        legend_vp <- get_legend(
          # create some space to the left of the legend
          vp_iva + theme(legend.box.margin = margin(0, 0, 0, 12))
        )
        
        ## combine plots in a grid
        png("figures/network/bivariate_descriptive/plots_focal.vector.png", width=9,height=6,units="in",res=600)
        plot_grid(bp_focal.v + mytheme,
                  vp_iva + theme(legend.position = "none") + mytheme,
                  legend_vp,
                  nrow = 2, ncol = 2
        )
        dev.off()
        
    ## kin-house infestation status vs. degree.cat (barplot)
    
        ## create labels for categories
        label_infected <- "One or more kin-houses infested \nwith T. cruzi-infected T. infestans"
        label_infested <- "One or more kin-houses \ninfested with T. infestans"
        label_none <- "No kin-houses infested"
        
        ## create cats
        df$kin.v <- NA
        df$kin.v <- factor(ifelse(df$kin.iva.sum > 0, label_infected,
                                  ifelse(df$kin.inf.prop > 0, label_infested,
                                         ifelse(df$degree == 0, NA,
                                                ifelse(df$kin.inf.prop == 0, label_none, df$kin.v)))),
                           levels = c(label_none, label_infested, label_infected)
                           )
        
        ## for degree, recat households with no kin ties as NA
        df$degree.nonzero <- factor(ifelse(df$degree.cat == "0", NA, as.character(df$degree.cat)),
                                    levels = c("1-2", "3+")
        )
        
        ## chi-square test excluding households w/ no kin connections
        table(df$degree.nonzero, df$kin.v)
        chisq_kin.v <- chisq.test(table(df$degree.nonzero, df$kin.v))
        
        # barplot
        library(ggplot2)
        library(scales)  # for percent formatting
        bp_kin.v <- df %>% filter(!is.na(degree.nonzero) & !is.na(kin.v)) %>%
                    ggplot(aes(x = degree.nonzero, fill = kin.v)) +
                    geom_bar(position = "fill", color = "black", width = 0.5) +  # stack by proportion
                    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
                    # scale_fill_manual(values = c("#440154FF", "#21908CFF","#FDE725FF")) +
                    scale_fill_manual(values = c("#7AD1514C", "#7AD15199","#7AD151FF")) +
                    labs(
                          # title = "Proportion of houses by kin-house infestation status among connected houses",
                          x = "Kinship degree",
                          y = "Proportion of houses",
                          fill = "Infestation status"
                        ) +
                    theme_minimal() 
       
        
    ## kin-house infected vector abundance vs. degree.cat
        
        ## kin.iva.med vs. degree.cat (excluding houses w/ no connections and no infected vector)        
        df_nozero <- df[df$degree.cat != "0" & df$kin.iva.med>0,]
        table(df_nozero$degree.cat, df_nozero$kin.iva.med)
        vp_kin.iva.med <- ggplot(df_nozero, aes(x = degree.cat, y = kin.iva.med)) +
            geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
            geom_boxplot(width = 0.1, fill = "white") + 
            scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
            ylim(0,13) +
            labs(
                 # title = "Median abundance of T.cruzi-infected T. infestans \namong kin-houses (n = 19)", 
                 x = "Kinship degree", y = "Median abundance of \ninfected vector among kin-houses",
                 fill = "Degrees \nof kinship") +
            theme_minimal()
    
        # ## kin.iva.sum vs. degree.cat (excluding houses w/ no connections and no infected vector)   
        # df_nozero <- df[df$degree.cat != "0" & df$kin.iva.sum>0,]
        # table(df$kin.iva.sum)
        # vp_kin.iva.sum <- ggplot(df_nozero, aes(x = degree.cat, y = kin.iva.sum)) +
        #     geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
        #     geom_boxplot(width = 0.1, fill = "white") + 
        #     scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
        #     ylim(0,13) +
        #     labs(
        #          # title = "Cumulative abundance of T. cruzi-infected T. infestans \namong kin-houses (n = 19)", 
        #          x = "Kinship degree", y = "Cumulative abundance of \ninfected vector among kin-houses",
        #          fill = "Degrees \nof kinship") +
        #     theme_minimal()
    
        ## get bar plot legend
        legend_vp <- get_legend(
          # create some space to the left of the legend
          vp_kin.iva.med + theme(legend.box.margin = margin(0, 0, 0, 12))
        )
        
    ## combine plots and save
    png("figures/network/bivariate_descriptive/plots_kin.vector.png", width=9,height=3,units="in",res=600)
    plot_grid(bp_kin.v,
              vp_kin.iva.med + theme(legend.position="none") + mytheme,
              # vp_kin.iva.sum + theme(legend.position="none") + mytheme,
              nrow = 1, ncol = 2
    )
    dev.off()

    ###########################################################
    ### BEGIN INTERMISSION: DROP OUTLIER WHERE df$iva == 22 ###
    ###########################################################
        
    # # Drop outlier
    # df <- df[df$iva != 22,]
    # 
    # # Vector infestation
    # 
    # ## focal house infestation status vs. degree.cat
    # 
    #     ## create labels for categories
    #     label_infected <- "Infested \nwith T. cruzi-infected T. infestans"
    #     label_infested <- "Infested with T. infestans"
    #     label_none <- "Non-infested"
    #     
    #     ## create cats
    #     df$focal.v <- NA
    #     df$focal.v <- factor(ifelse(df$iva > 0 & df$v.inf == "infested", label_infected,
    #                                 ifelse(df$v.inf == "infested", label_infested,
    #                                        ifelse(df$v.inf == "non-infested", label_none, df$focal.v))),
    #                          levels = c(label_none, label_infested, label_infected)
    #     )
    #     
    #     # ## for degree, recat households with no kin ties as NA
    #     # df$degree.nonzero <- factor(ifelse(df$degree.cat == "0", NA, as.character(df$degree.cat)),
    #     #                             levels = c("1-2", "3+")
    #     # )
    #     
    #     ## chi-square test excluding households w/ no kin connections
    #     table(df$degree.cat, df$focal.v)
    #     chisq_focal.v <- chisq.test(table(df$degree.cat, df$focal.v))
    #     
    #     # barplot
    #     library(ggplot2)
    #     bp_focal.v <- df %>%
    #       filter(!is.na(degree.cat) & !is.na(focal.v)) %>%
    #       ggplot(aes(x = degree.cat, fill = focal.v)) +
    #       # geom_bar(position = "fill", alpha=0.7, width = 0.5) +  # stack by proportion
    #       geom_bar(position = "fill", color = "black", width = 0.5) +  # stack by proportion
    #       scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    #       # scale_fill_manual(values = c("#440154FF", "#21908CFF","#FDE725FF")) +
    #       scale_fill_manual(values = c("#41448733", "#41448799","#414487FF")) +
    #       labs(
    #         # title = "Proportion of houses by infestation status",
    #         x = "Kinship degree",
    #         y = "Proportion of houses",
    #         fill = "Infestation status"
    #       ) +
    #       theme_minimal()
    #     
    #     ## iva vs. degree.cat
    #     df_nozero <- df[df$iva>0,]
    #     length(which(!is.na(df_nozero$iva))) # n=25
    #     vp_iva <- ggplot(df_nozero, aes(x = degree.cat, y = iva)) +
    #       geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
    #       geom_boxplot(width = 0.1, fill = "white") + 
    #       scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
    #       labs(
    #         # title = "Relative abundance of infected vector among \nhouses with at least one infected vector (n = 25)", 
    #         x = "Kinship degree", y = "Abundance of infected vector",
    #         fill = "Degrees \nof kinship") +
    #       theme_minimal()
    #     
    #     ## get legend
    #     legend_vp <- get_legend(
    #       # create some space to the left of the legend
    #       vp_iva + theme(legend.box.margin = margin(0, 0, 0, 12))
    #     )
    #     
    #     ## kruskal wallis test
    #     test <- kruskal.test(iva ~ degree.cat, data = df)
    #     test_nozero <- kruskal.test(iva ~ degree.cat, data = df_nozero)
    #         
    #     ## combine plots in a grid
    #     png("figures/network/bivariate_descriptive/plots_focal.vector_drop22.png", width=9,height=6,units="in",res=600)
    #     plot_grid(bp_focal.v + mytheme,
    #               vp_iva + theme(legend.position = "none") + mytheme,
    #               legend_vp,
    #               nrow = 2, ncol = 2
    #     )
    #     dev.off()
    #     
    # ## kin-house infestation status vs. degree.cat (barplot)
    # 
    #     ## create labels for categories
    #     label_infected <- "One or more kin-houses infested \nwith T. cruzi-infected T. infestans"
    #     label_infested <- "One or more kin-houses \ninfested with T. infestans"
    #     label_none <- "No kin-houses infested"
    #     
    #     ## create cats
    #     df$kin.v <- NA
    #     df$kin.v <- factor(ifelse(df$kin.iva.sum > 0, label_infected,
    #                               ifelse(df$kin.inf.prop > 0, label_infested,
    #                                      ifelse(df$degree == 0, NA,
    #                                             ifelse(df$kin.inf.prop == 0, label_none, df$kin.v)))),
    #                        levels = c(label_none, label_infested, label_infected)
    #     )
    #     
    #     ## for degree, recat households with no kin ties as NA
    #     df$degree.nonzero <- factor(ifelse(df$degree.cat == "0", NA, as.character(df$degree.cat)),
    #                                 levels = c("1-2", "3+")
    #     )
    #     
    #     ## chi-square test excluding households w/ no kin connections
    #     table(df$degree.nonzero, df$kin.v)
    #     chisq_kin.v <- chisq.test(table(df$degree.nonzero, df$kin.v))
    #     
    #     # barplot
    #     library(ggplot2)
    #     library(scales)  # for percent formatting
    #     bp_kin.v <- df %>% filter(!is.na(degree.nonzero) & !is.na(kin.v)) %>%
    #       ggplot(aes(x = degree.nonzero, fill = kin.v)) +
    #       geom_bar(position = "fill", color = "black", width = 0.5) +  # stack by proportion
    #       scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    #       # scale_fill_manual(values = c("#440154FF", "#21908CFF","#FDE725FF")) +
    #       scale_fill_manual(values = c("#7AD1514C", "#7AD15199","#7AD151FF")) +
    #       labs(
    #         # title = "Proportion of houses by kin-house infestation status among connected houses",
    #         x = "Kinship degree",
    #         y = "Proportion of houses",
    #         fill = "Infestation status"
    #       ) +
    #       theme_minimal() 
    #     
    # 
    # ## kin-house infected vector abundance vs. degree.cat
    # 
    #     ## kin.iva.med vs. degree.cat (excluding houses w/ no connections and no infected vector)        
    #     df_nozero <- df[df$degree.cat != "0" & df$kin.iva.med>0,]
    #     table(df_nozero$degree.cat, df_nozero$kin.iva.med)
    #     vp_kin.iva.med <- ggplot(df_nozero, aes(x = degree.cat, y = kin.iva.med)) +
    #       geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
    #       geom_boxplot(width = 0.1, fill = "white") + 
    #       scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
    #       ylim(0,13) +
    #       labs(
    #         # title = "Median abundance of T.cruzi-infected T. infestans \namong kin-houses (n = 19)", 
    #         x = "Kinship degree", y = "Median abundance of \ninfected vector among kin-houses",
    #         fill = "Degrees \nof kinship") +
    #       theme_minimal()
    #     
    #     # ## kin.iva.sum vs. degree.cat (excluding houses w/ no connections and no infected vector)   
    #     # df_nozero <- df[df$degree.cat != "0" & df$kin.iva.sum>0,]
    #     # table(df$kin.iva.sum)
    #     # vp_kin.iva.sum <- ggplot(df_nozero, aes(x = degree.cat, y = kin.iva.sum)) +
    #     #     geom_violin(trim = TRUE, aes(fill = degree.cat), alpha = 0.7) + 
    #     #     geom_boxplot(width = 0.1, fill = "white") + 
    #     #     scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
    #     #     ylim(0,13) +
    #     #     labs(
    #     #          # title = "Cumulative abundance of T. cruzi-infected T. infestans \namong kin-houses (n = 19)", 
    #     #          x = "Kinship degree", y = "Cumulative abundance of \ninfected vector among kin-houses",
    #     #          fill = "Degrees \nof kinship") +
    #     #     theme_minimal()
    #     
    #     ## get bar plot legend
    #     legend_vp <- get_legend(
    #       # create some space to the left of the legend
    #       vp_kin.iva.med + theme(legend.box.margin = margin(0, 0, 0, 12))
    #     )
    #     
    #     ## combine plots and save
    #     png("figures/network/bivariate_descriptive/plots_kin.vector_drop22.png", width=9,height=3,units="in",res=600)
    #     plot_grid(bp_kin.v,
    #               vp_kin.iva.med + theme(legend.position="none") + mytheme,
    #               # vp_kin.iva.sum + theme(legend.position="none") + mytheme,
    #               nrow = 1, ncol = 2
    #     )
    #     dev.off()

    ###########################################################
    ### END INTERMISSION ###
    ### Droppping outlier made no qualitative difference to our results
    ###########################################################

# Serology variables
    
    ## serology vs. degree.cat
    table(df$serology, df$degree.cat)
    length(which(!is.na(df$serology)))
    vp_sero.any <- ggplot(df, aes(x = degree.cat, y = serology)) +
      geom_violin(trim = FALSE, aes(fill = degree.cat), alpha = 0.7) + 
      geom_boxplot(width = 0.1, fill = "white") + 
      scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
      labs(title = "Number of seropositive persons \nper household (n = 303)", 
           x = "Kinship degree", y = "Number of seropositive persons",
           fill = "Degrees \nof kinship") +
      theme_minimal()
        
        # kruskal-wallis test        
        kruskal.test(serology ~ degree.cat, data = df)
    
    ## serology.child vs. degree.cat
    table(df$serology.child, df$degree.cat)
    length(which(!is.na(df$serology.child)))
    vp_sero.child <- ggplot(df, aes(x = degree.cat, y = serology.child)) +
      geom_violin(trim = FALSE, aes(fill = degree.cat), alpha = 0.7) + 
      geom_boxplot(width = 0.1, fill = "white") + 
      scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
      labs(title = "Number of seropositive children \nper household (n = 226)", 
           x = "Kinship degree", y = "Number of seropositive children",
           fill = "Degrees \nof kinship") +
      theme_minimal()
    
        # kruskal-wallis test        
        kruskal.test(serology.child ~ degree.cat, data = df)
    
    ## repeat above but exclude households with no seropositive inhabitants or children
        
        ## serology vs. degree.cat
        df_nozero <- df[df$serology>0,]
        df_nozero <- df_nozero[!is.na(df_nozero$degree.cat),]
        length(which(!is.na(df_nozero$serology))) # n=214
        table(df_nozero$serology, df_nozero$degree.cat)
        vp_sero.1any <- ggplot(df_nozero, aes(x = degree.cat, y = serology)) +
          geom_violin(trim = FALSE, aes(fill = degree.cat), alpha = 0.7) + 
          geom_boxplot(width = 0.1, fill = "white") + 
          scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
          labs(title = "Number of seropositive persons among houses with at least one seropositive inhabitant (n = 214)", 
               x = "Kinship degree", y = "Number of seropositive persons",
               fill = "Degrees \nof kinship") +
          theme_minimal()
        
            # kruskal-wallis test        
            kruskal.test(serology ~ degree.cat, data = df_nozero)
            
        ## serology.child vs. degree.cat
        df_nozero <- df[df$serology.child>0,]
        df_nozero <- df_nozero[!is.na(df_nozero$degree.cat),]
        length(which(!is.na(df_nozero$serology.child))) # n=50
        table(df_nozero$serology.child, df_nozero$degree.cat)
        vp_sero.1child <- ggplot(df_nozero, aes(x = degree.cat, y = serology.child)) +
          geom_violin(trim = FALSE, aes(fill = degree.cat), alpha = 0.7) + 
          geom_boxplot(width = 0.1, fill = "white") + 
          scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
          labs(title = "Number of seropositive children among houses with at least one seropositive child (n = 50)", 
               x = "Kinship degree", y = "Number of seropositive children",
               fill = "Degrees \nof kinship") +
          theme_minimal()
    
            # kruskal-wallis test        
            kruskal.test(serology.child ~ degree.cat, data = df_nozero)
        
    # ## kin.sero.prop vs. degree.cat
    # table(df$kin.sero.prop, df$degree.cat)
    # length(which(!is.na(df$kin.sero.prop)))
    # summary(df$kin.sero.prop)
    # vp_kin.sero.prop <- ggplot(df, aes(x = degree.cat, y = kin.sero.prop)) +
    #               geom_violin(trim = FALSE, aes(fill = degree.cat), alpha = 0.7) + 
    #               geom_boxplot(width = 0.1, fill = "white") + 
    #               scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
    #               ylim(0,1) +
    #               labs(title = "Proportion of connected houses \nwith >1 seropositive person (n = 297)", 
    #                    x = "Kinship degree", y = "Proportion of connected houses",
    #                    fill = "Degrees \nof kinship") +
    #               theme_minimal()
    
    ## kin.sero.med vs. degree.cat
    df_nozero <- df[df$degree.cat != "0" & df$kin.sero.med>0,]
    df_nozero <- df_nozero[!is.na(df_nozero$degree.cat),] #filter out NAs
    table(df_nozero$degree.cat, df_nozero$kin.sero.med)
        length(which(!is.na(df_nozero$kin.sero.med))) #152
        summary(df_nozero$kin.sero.med)
    vp_kin.sero.med <- ggplot(df_nozero, aes(x = degree.cat, y = kin.sero.med)) +
                  geom_violin(trim = FALSE, aes(fill = degree.cat), alpha = 0.7) + 
                  geom_boxplot(width = 0.1, fill = "white") + 
                  scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
                  ylim(0,8) +
                  labs(title = "Median number of seropositive persons \namong kin networks (n = 152)", 
                       x = "Kinship degree", y = "Median number of seropositive persons",
                       fill = "Degrees \nof kinship") +
                  theme_minimal()
      
        ## mann-whitney/wilcox rank sum test
        wilcox.test(kin.sero.med ~ degree.cat, data = df_nozero)

    ## kin.sero.sum vs. degree.cat
    df_nozero <- df[df$degree.cat != "0" & df$kin.sero.sum>0,]
    df_nozero <- df_nozero[!is.na(df_nozero$degree.cat),] #filter out NAs
    table(df_nozero$degree.cat, df_nozero$kin.sero.sum)
        length(which(!is.na(df_nozero$kin.sero.sum))) # n = 152
        summary(df_nozero$kin.sero.sum)
    vp_kin.sero.sum <- ggplot(df_nozero, aes(x = degree.cat, y = kin.sero.sum)) +
                  geom_violin(trim = FALSE, aes(fill = degree.cat), alpha = 0.7) + 
                  geom_boxplot(width = 0.1, fill = "white") + 
                  scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
                  ylim(0,25) +
                  labs(title = "Cumulative number of seropositive persons \namong kin networks (n = 152)", 
                       x = "Kinship degree", y = "Cumulative number of seropositive persons",
                       fill = "Degrees \nof kinship") +
                  theme_minimal()
    
        ## mann-whitney/wilcox rank sum test
        wilcox.test(kin.sero.sum ~ degree.cat, data = df_nozero)
    
    ## repeat but for children
        
        # ## kin.serochild.prop vs. degree.cat
        # table(df$kin.serochild.prop, df$degree.cat)
        # length(which(!is.na(df$kin.serochild.prop)))
        # summary(df$kin.serochild.prop)
        # vp_kin.serochild.prop <- ggplot(df, aes(x = degree.cat, y = kin.serochild.prop)) +
        #               geom_violin(trim = FALSE, aes(fill = degree.cat), alpha = 0.7) + 
        #               geom_boxplot(width = 0.1, fill = "white") + 
        #               scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
        #               ylim(0,1) +
        #               labs(title = "Proportion of connected houses \nwith >1 seropositive person (n = 297)", 
        #                    x = "Kinship degree", y = "Proportion of connected houses",
        #                    fill = "Degrees \nof kinship") +
        #               theme_minimal()
        
        ## kin.serochild.med vs. degree.cat
        df_nozero <- df[df$degree.cat != "0" & df$kin.serochild.med>0,]
        df_nozero <- df_nozero[!is.na(df_nozero$degree.cat),] #filter out NAs
        table(df_nozero$degree.cat, df_nozero$kin.serochild.med)
            length(which(!is.na(df_nozero$kin.serochild.med))) # n = 52
            summary(df_nozero$kin.serochild.med)
        vp_kin.serochild.med <- ggplot(df_nozero, aes(x = degree.cat, y = kin.serochild.med)) +
          geom_violin(trim = FALSE, aes(fill = degree.cat), alpha = 0.7) + 
          geom_boxplot(width = 0.1, fill = "white") + 
          scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
          ylim(0,8) +
          labs(title = "Median number of seropositive children \namong kin networks (n = 52)", 
               x = "Kinship degree", y = "Median number of seropositive children",
               fill = "Degrees \nof kinship") +
          theme_minimal()
        
            ## mann-whitney/wilcox rank sum test
            wilcox.test(kin.serochild.med ~ degree.cat, data = df_nozero)
            
        ## kin.serochild.sum vs. degree.cat
        df_nozero <- df[df$degree.cat != "0" & df$kin.serochild.sum>0,]
        df_nozero <- df_nozero[!is.na(df_nozero$degree.cat),] #filter out NAs
        table(df_nozero$degree.cat, df_nozero$kin.serochild.sum)
            length(which(!is.na(df_nozero$kin.serochild.sum))) # n = 52
            summary(df$kin.serochild.sum)
        vp_kin.serochild.sum <- ggplot(df_nozero, aes(x = degree.cat, y = kin.serochild.sum)) +
          geom_violin(trim = FALSE, aes(fill = degree.cat), alpha = 0.7) + 
          geom_boxplot(width = 0.1, fill = "white") + 
          scale_fill_manual(values = c("0" = "#440154FF", "1-2" = "#21908CFF", "3+" = "#FDE725FF")) +
          ylim(0,25) +
          labs(title = "Cumulative number of seropositive children \namong kin networks (n = 335)", 
               x = "Kinship degree", y = "Cumulative number of seropositive children",
               fill = "Degrees \nof kinship") +
          theme_minimal()
        
            ## mann-whitney/wilcox rank sum test
            wilcox.test(kin.serochild.sum ~ degree.cat, data = df_nozero)
        
    vp_sero <- plot_grid(
                        vp_sero.child + theme(legend.position="none") + mytheme,
                        vp_sero.any + theme(legend.position="none") + mytheme,
                        legend,
                        nrow = 1, ncol = 2
    )    
    
    vp_sero.1 <- plot_grid(
                        vp_sero.1child + theme(legend.position="none") + mytheme,
                        vp_sero.1any + theme(legend.position="none") + mytheme,
                        legend,
                        nrow = 1, ncol = 2
    )    
    
    vp_kin.sero <- plot_grid(
                   vp_kin.sero.med + theme(legend.position="none") + mytheme,
                   vp_kin.sero.sum + theme(legend.position="none") + mytheme,
                   legend,
                   nrow = 1, ncol = 2
                   )
    
    vp_kin.serochild <- plot_grid(
                    vp_kin.serochild.med + theme(legend.position="none") + mytheme,
                    vp_kin.serochild.sum + theme(legend.position="none") + mytheme,
                    legend,
                    nrow = 1, ncol = 2
    )
    
# Save plots
    
    png("figures/network/bivariate_descriptive/violin_inhabitants.png",width=8,height=8,units="in",res=600)
        vp_hab
        dev.off()
    png("figures/network/bivariate_descriptive/violin_serology.png",width=8,height=5,units="in",res=600)
        vp_sero
        dev.off()
    png("figures/network/bivariate_descriptive/violin_serology_atleastone.png",width=8,height=5,units="in",res=600)
        vp_sero.1
        dev.off()
    png("figures/network/bivariate_descriptive/violin_kinsero.png",width=8,height=5,units="in",res=600)
        vp_kin.sero
        dev.off()
    png("figures/network/bivariate_descriptive/violin_kinserochild.png",width=8,height=5,units="in",res=600)
        vp_kin.serochild
        dev.off()
    png("figures/network/bivariate_descriptive/violin_svikinsvi.png",width=8,height=8,units="in",res=600)
        vp_kin.svi
        dev.off()
    png("figures/network/bivariate_descriptive/violin_kiniva.png",width=8,height=8,units="in",res=600)
        vp_kin.iva
        dev.off()
    


    
### Goal ### create a bivariate table to see how our seropositivity, sociodemographics, and ecological data change across varied levels of kin connection

# Load data
str(df_kin)

# Transform infested to factor table(df_kin$iva.bin)
df_kin$iva.bin <- factor(df_kin$iva.bin,
                         levels = c(0,1),
                         labels = c("Not infested", "Infested"))

# Review variables
table(df_kin$degree.cat)
table(df_kin$serology)
table(df_kin$serology.occ)
table(df_kin$serology.cat)
table(df_kin$serology.child)
table(df_kin$serology.child.occ)
summary(df_kin$hab_tot) # Num of total inhabitants
summary(df_kin$hab_0a14) # Num of inhabitants 0-14 yo
summary(df_kin$age) # Mean age per hhx
table(df_kin$ethnicity)
table(df_kin$mobility)
summary(df_kin$svi)
summary(df_kin$hai)
table(df_kin$v.inf)
table(df_kin$iva.bin) # Infested or not infested w/ infected vector
table(df_kin$iva.bin1) # Abundance of infected vector for houses that are infested
summary(df_kin$betweenness)
summary(df_kin$harmonic)
summary(df_kin$kin.sero.prop)
table(df_kin$kin.sero.med)
table(df_kin$kin.sero.sum)
summary
# Create table labels
library(table1)
table1::label(df_kin$degree.cat) <- "No. of households connected by kinship"
table1::label(df_kin$serology) <- "No. of seropositive persons (disc.)"
table1::label(df_kin$serology.occ) <- "Occurrence of a seropositive person"
table1::label(df_kin$serology.cat) <- "No. of seropositive persons (cat.)"
table1::label(df_kin$serology.child) <- "No. of seropositive children (disc.)"
table1::label(df_kin$serology.child.occ) <- "Occurrence of a seropositive child"
table1::label(df_kin$serology.child.cat) <- "No. of seropositive children (cat.)"
table1::label(df_kin$hab_tot) <- "No. of total inhabitants"
table1::label(df_kin$hab_0a14) <- "No. of inhabitants aged 0 to 14"
table1::label(df_kin$age) <- "Mean age per household" 
table1::label(df_kin$ethnicity) <- "Ethnicity"
table1::label(df_kin$mobility) <- "Mobility"
table1::label(df_kin$svi) <- "Social vulnerability"
table1::label(df_kin$hai) <- "Host availability index"
table1::label(df_kin$v.inf) <- "Infestation w/ vector"
table1::label(df_kin$iva.bin) <- "Infestation (w/ infected vector)"
table1::label(df_kin$iva.bin1) <- "Abundance of infected vector among infested houses"
table1::label(df_kin$svi.kin) <- "Social vulnerability of kin connected households"
table1::label(df_kin$svi.kin.mean) <- "Social vulnerability of kin connected households w/ NAs replaced by mean svi.kin"
table1::label(df_kin$svi.kin.own) <- "Social vulnerability of kin connected households w/ NAs replaced by svi of focal house"


# Create table
table_hhx_bydegree <- table1::table1(~serology + serology.occ + serology.cat + serology.child + serology.child.occ + serology.child.cat + hab_tot + hab_0a14 + df_kin$age + df_kin$ethnicity + df_kin$mobility + df_kin$svi + df_kin$svi.kin + df_kin$svi.kin.mean + df_kin$svi.kin.own + df_kin$hai + df_kin$v.inf + df_kin$iva.bin + df_kin$iva.bin1 | degree.cat, 
                                     data = df_kin,
                                     overall=c(left="Total"), 
                                     caption="Ecological and Sociodemographic Characteristics of Households Stratified by Kinship Degreeᵃ", 
                                     footnote="ᵃ The number of other houses in the study area to which a particular house is connected by kinship")

# Save table
# Convert to flextable and save
library(flextable)
library(officer)
t1flex(table_hhx_bydegree) %>% save_as_docx(path="tables/network/table_hhx_bydegree.docx",
                                            pr_section = prop_section(page_size = page_size(orient = "landscape")))

# Check correlation
cor(df_kin$svi, df_kin$svi.kin.mean, use = "complete.obs")
cor(df_kin$svi, df_kin$svi.kin.own, use = "complete.obs")


################################################################################
### MOBILITY & KINSHIP 
################################################################################

### SETUP ###

# Clean environment
rm(list=ls()) 
graphics.off()

# Set directory
setwd("/Users/~/Desktop/ChagasNetworkPublic")

# Load data
load('data/CleanData.RData')

# get df 
df <- df_kin

# Prep dataframe and variables

## subset binary variables
df_bin <- subset(df, select = c(degree.cat, betweenness, harmonic, serology.occ, serology.child.occ, ethnicity, mobility, v.inf, iva.inf))

## add degree as binary variable
df_bin$degree.bin <- factor(ifelse(df_bin$degree.cat == "0", "0", ">0"),
                            levels = c("0", ">0"))

# for serology.child.occ, set observations with 'No child' to NA, and transform to binary
str(df_bin$serology.child.occ)
df_bin$serology.child.occ <- ifelse(df_bin$serology.child.occ == "No child", NA, 
                                    ifelse(df_bin$serology.child.occ == "Occurrence", 1, 0))

## transform remaining factor variables to binary 0/1 (numeric)
library(dplyr)
str(df_bin)
df_bin <- df_bin %>% mutate(serology.occ = case_when(serology.occ == "Occurrence" ~ 1,
                                                     serology.occ == "No occurrence" ~ 0))
df_bin <- df_bin %>% mutate(ethnicity = case_when(ethnicity == "Qom" ~ 1,
                                                  ethnicity == "Creole" ~ 0))
df_bin <- df_bin %>% mutate(mobility = case_when(mobility == "mover" ~ 1,
                                                 mobility == "non-mover" ~ 0))
df_bin <- df_bin %>% mutate(v.inf = case_when(v.inf == "infested" ~ 1,
                                              v.inf == "non-infested" ~ 0))

### ANALYSIS ###

df_bin$connected <- as.numeric(df_bin$degree.bin)
str(df_bin)
prop.table(df_bin$mobility, df_bin$connected)

unique(df_kin$kin.inf.prop)

df_kin$kin.v.inf <- ifelse(is.na(df_kin$kin.inf.prop), NA, 
                           ifelse(df_kin$kin.inf.prop > 0, 1, 0))
View(df_kin[,c("kin.v.inf", "kin.inf.prop")])

chisq.test(table(df_kin$mobility, df_kin$kin.v.inf))  
chisq.test(table(df_kin$mobility, df_kin$kin.iva.sum))  

