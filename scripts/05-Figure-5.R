##
##      __________  _____
##     / ____/ __ \/ ___/
##    / /   / / / /\__ \ 
##   / /___/ /_/ /___/ / 
##   \____/_____//____/  
##
##  Creative Data Solutions
##  Vanderbilt University
##  https://cds.vanderbilt.edu
##
## Authors: Drs. Shristi Shrestha & JP Cartailler
## Date Created: 2022-08-10
## 

# libraries
library(TCseq)
library(dplyr)
library(viridis)
library(ggplot2)
library(ggpubr)
library(svglite)
library(tidyr)
library(ggsci)
library(ggthemes)


# Figure 5C----

# This input data contains regulons and its target gene counts obtained from pyscenic analysis on Adult, Old Adult, young Adult and Infants age group subset

input_data <- "data/Combined_common_TFtargets_Adult_OldAdults_YoungAdults_Infants_08032021.csv"
data = read.csv(input_data, row.names="TF", sep=",", stringsAsFactors = FALSE) # 

# Column 5:8 description for following: ratios against infanttarget as baseline
# Infants = infanttarget / infanttarget
# YoungAdults = youngadulttarget / infanttarget
# Adults = adulttarget / infanttarget
# OldAdult = oldtarget / infanttarget

data = data[,5:8] 
data <- as.matrix(data) # 1059 TFs


#Filter out regulons that contained only one motif from motif annotation database used in pyscenic analysis 
# "motifs-v9-nr.hgnc-m0.001-o0.0.tbl" motif annotation reference was used from https://resources.aertslab.org/cistarget/motif2tf/ 
# 
# list of regulons with one motif to remove:
tf_to_remove_file <- "data/worksheets/1motifTF_387.csv"
tf_to_remove = read.csv(tf_to_remove_file, sep=",", stringsAsFactors = FALSE) # 
tf_to_remove <- tf_to_remove$gene_name

data_clean = data[!row.names(data) %in% tf_to_remove,]
dim(data_clean) # 764 TFs

# Filter based on GRN networks, age-group specific

# this following input worksheet per age group contains top 15 % regulons ranked by their pyscenic "importance" metric(measure of motif-TF association strength) 
# Cells are only included from healthy beta cells

tf_adult = read.csv("02_regulons__adj_mat_Cytoscape_thres15pct_1motifTF_removed__beta_healthy_adults.100X.csv", sep=",", stringsAsFactors = FALSE) # 
tf_adult <- unique(tf_adult$TF) # 255 TFs

tf_infants = read.csv("02_regulons__adj_mat_Cytoscape_thres15pct_1motifTF_removed__beta_healthy_infants.100X.csv", sep=",", stringsAsFactors = FALSE) # 
tf_infants <- unique(tf_infants$TF) # 278 TFs

tf_old = read.csv("02_regulons__adj_mat_Cytoscape_thres15pct_1motifTF_removed__beta_healthy_old.100X.csv", sep=",", stringsAsFactors = FALSE) # 
tf_old <- unique(tf_old$TF) # 288 TFs

tf__infants_adult_old <- unique(c(tf_infants, tf_adult, tf_old)) # 499 unique TFs

#subset only unique ones
data_clean2 = data_clean[row.names(data_clean) %in% tf__infants_adult_old,]
data <- data_clean2 # 379 TFs


# Analysis
# The temporal patterns are analyzed using timeclust(tca, algo = "cm", k = 6, standardize = TRUE)
# 
# algo - # Two types of clustering algorithms are included in the package: 
# hard clustering (hierachical, pam, kmeans) and soft clustering (fuzzy cmeans):
# c("pam", "km", "hc", "cm")
# 
# k - 
# 
# standardize - instead of absolute value of different time series, one might only focus on the change patterns and
# expect time series with similar pattern to be cluster in same group. In this case, "standardize"
# parameter gives an option to perform z-score transformation on the data to be clustered, which
# reduces the noises introduced by the difference in the absolute values.

algos <- c("pam", "km", "hc", "cm")
algos <- c("cm")

#loop through the different algorithms
for(param_algo in algos){
    
    for(param_k in c(6)){
        
        #param_k = 6
        #param_algo = "cm"
        tca <- timeclust(data_clean, algo = param_algo, k = param_k, standardize = TRUE) # expects data matrix to be stricly numbers, so no extra annotations
        
        col = viridis(20, direction = -1) # generate viridis color palette
        p <- timeclustplot(tca, categories="Age Group", value="Relative FC of Target Genes", cols = 2, membership.color=col, 
                           axis.text.size=10, legend.text.size=10, legend.title.size=10, title.size=12, axis.title.size=10)
        
        # Fix ggplot objects
        # 1. 45 angle title on x-axis labels
        # 2. # of TFs per cluster in title
        
        plots = list()
        for (i in 1:length(p)){
            num_tf_in_cluster = table(tca@cluster)[[i]]
            p[[i]] = p[[i]]  + theme(axis.text.x = element_text(angle = 45, hjust=1)) + ggtitle(paste0("Cluster ", i, " (", num_tf_in_cluster, " TFs)"))
            plots[[i]] <- p[[i]]
            
            # save individual plot
            filename <- paste0(output_dir, "param_algo-", param_algo, "__param_k-", param_k, ".", i , ".svg")
            ggsave(filename, plot = p[[i]], scale=2, width=2.25, height=2)
        }
        #plots[[1]]
        figure <- ggarrange(plotlist = plots, ncol=2, nrow = 3)
        
        plot_title = paste0("Clustering used: ", param_algo, ", k: ", param_k)
        figure <- annotate_figure(figure,
                                  top = text_grob(plot_title, face = "bold", size = 12),
                                  bottom = text_grob(paste0("Input data: ", input_data, "\n Hard clustering: hierachical (hc), pam (pam), kmeans (km)\n Soft clustering: fuzzy cmeans (cm)"), hjust = 1, x = 1, size = 9)
        )
        #figure
        filename <- paste0(output_dir, "param_algo-", param_algo, "__param_k-", param_k, ".png")
        ggsave(filename, plot = figure, scale=2, width=3, height=5)
        
        filename <- paste0(output_dir, "param_algo-", param_algo, "__param_k-", param_k, ".svg")
        ggsave(filename, plot = figure, scale=2, width=3, height=5)
        
        # Summary table
        summary <- as.data.frame(table(tca@cluster))
        colnames(summary) <- c("cluster", "num_TF")
        #summary
        write.csv(summary, paste0(output_dir, "param_algo-", param_algo, "__param_k-", param_k, "_SUMMARY.csv"), row.names = FALSE)
        
        # Detailed TF/cluster assignments
        # 
        detail <- as.data.frame(tca@cluster)
        colnames(detail) <- c("cluster")
        write.csv(detail, paste0(output_dir, "param_algo-", param_algo, "__param_k-", param_k, "_DETAILS.csv"), row.names = TRUE)
        
        
    }
}

plots_backup <- plots

param_algo = "cm"
param_k = 6
tca <- timeclust(data, algo = param_algo, k = param_k, standardize = TRUE) # expects data matrix to be stricly numbers, so no extra annotations

# Sandbox - panels ----
#
regression_data <- data.frame(tca@data)
regression_data$tf <- rownames(regression_data)
regression_data$cluster <- tca@cluster

data_long <- gather(regression_data, age_group, measurement, Infants:OldAdult, factor_key=TRUE)

p <- ggplot(data_long, aes(age_group, measurement, group=tf)) + geom_line()
p + facet_grid(rows = vars(cluster))



figure <- ggplot(data_long, aes(age_group, measurement,color=as.factor(cluster), group=as.factor(cluster))) + 
    geom_point() + 
    geom_smooth(method = "loess", level=0.95) +
    scale_fill_tron() +
    scale_color_tron() +
    theme_bw()

filename <- paste0(output_dir, "combined__param_algo-", param_algo, "__param_k-", param_k, ".png")
ggsave(filename, plot = figure, scale=2, width=3, height=3)

figure <- ggplot(data_long, aes(x=age_group, y=measurement, group=as.factor(cluster), colour=as.factor(cluster), fill=as.factor(cluster))) + 
    #geom_smooth(span = 1, method="gam", aes(color=as.factor(cluster))) +
    geom_point() + 
    geom_smooth(method="loess",  level = 0.99, aes(color=as.factor(cluster))) +
    scale_colour_tableau('Tableau 10') + 
    scale_fill_tableau('Tableau 10') + 
    theme_bw()

filename <- paste0(output_dir, "combined2__param_algo-", param_algo, "__param_k-", param_k, ".png")
ggsave(filename, plot = figure, scale=2, width=3, height=3)

#
p

# Sandbox - aggregate  ----
# 
agg_data <- data.frame(tca@data)
agg_data$tf <- rownames(agg_data)
agg_data$cluster <- tca@cluster

## mean
mean_Infants <- aggregate(agg_data$Infants, by=list(agg_data$cluster), FUN=mean)
mean_YoungAdults <- aggregate(agg_data$YoungAdults, by=list(agg_data$cluster), FUN=mean)
mean_Adults <- aggregate(agg_data$Adults, by=list(agg_data$cluster), FUN=mean)
mean_OldAdult <- aggregate(agg_data$OldAdult, by=list(agg_data$cluster), FUN=mean)

## sd
sd_Infants <- aggregate(agg_data$Infants, by=list(agg_data$cluster), FUN=sd)
sd_YoungAdults <- aggregate(agg_data$YoungAdults, by=list(agg_data$cluster), FUN=sd)
sd_Adults <- aggregate(agg_data$Adults, by=list(agg_data$cluster), FUN=sd)
sd_OldAdult <- aggregate(agg_data$OldAdult, by=list(agg_data$cluster), FUN=sd)

## max
max_Infants <- aggregate(agg_data$Infants, by=list(agg_data$cluster), FUN=max)
max_YoungAdults <- aggregate(agg_data$YoungAdults, by=list(agg_data$cluster), FUN=max)
max_Adults <- aggregate(agg_data$Adults, by=list(agg_data$cluster), FUN=max)
max_OldAdult <- aggregate(agg_data$OldAdult, by=list(agg_data$cluster), FUN=max)

## min
min_Infants <- aggregate(agg_data$Infants, by=list(agg_data$cluster), FUN=min)
min_YoungAdults <- aggregate(agg_data$YoungAdults, by=list(agg_data$cluster), FUN=min)
min_Adults <- aggregate(agg_data$Adults, by=list(agg_data$cluster), FUN=min)
min_OldAdult <- aggregate(agg_data$OldAdult, by=list(agg_data$cluster), FUN=min)

combined <- mean_Infants
colnames(combined) <- c('cluster', 'mean_Infants')
#combined$cluster = as.numeric(as.character(combined$cluster))
combined$mean_YoungAdults <- mean_YoungAdults$x
combined$mean_Adults <- mean_Adults$x
combined$mean_OldAdult <- mean_OldAdult$x
combined$sd_Infants <- sd_Infants$x
combined$sd_YoungAdults <- sd_YoungAdults$x
combined$sd_Adults <- sd_Adults$x
combined$sd_OldAdult <- sd_OldAdult$x
combined$max_Infants <- max_Infants$x
combined$max_YoungAdults <- max_YoungAdults$x
combined$max_Adults <- max_Adults$x
combined$max_OldAdult <- max_OldAdult$x
combined$min_Infants <- min_Infants$x
combined$min_YoungAdults <- min_YoungAdults$x
combined$min_Adults <- min_Adults$x
combined$min_OldAdult <- min_OldAdult$x

data_long <- gather(combined, age_group, measurement, mean_Infants:mean_OldAdult, factor_key=TRUE)

write.csv(data_long, paste0(output_dir, "data_long.csv"), row.names = FALSE)

# now you have to manually re-arrange the data, moving measurements to col2, then only preserving the mean/min/max for the 
# appropriate age groups, removing the uneeded ones, then collapsing them into single mean/min/max cols
data_long2 <- read.csv(paste0(output_dir, "data_long__MANUALLY_FIXED.csv"), stringsAsFactors = TRUE)
data_long2$age_group <- factor(data_long2$age_group, levels = c("mean_Infants", "mean_YoungAdults", "mean_Adults", "mean_OldAdult"))

# min/max, adjusted
ggplot(data_long2, aes(x=age_group, y=measurement, ymin=(measurement-((measurement-min)/5)), ymax=(((max-measurement)/5)+measurement), group=cluster, fill=as.factor(cluster))) + 
    geom_line(aes(color=as.factor(cluster))) +
    geom_ribbon(alpha=0.1) + 
    scale_fill_tron() + 
    theme_bw()

# min/max
ggplot(data_long2, aes(x=age_group, y=measurement, ymin=min, ymax=max, group=cluster, fill=as.factor(cluster))) + 
    geom_line(aes(color=as.factor(cluster))) +
    geom_ribbon(alpha=0.1) + 
    scale_fill_tron() + 
    theme_bw()

# sd
ggplot(data_long2, aes(x=age_group, y=measurement, ymin=(measurement-sd), ymax=(measurement+sd), group=cluster, fill=as.factor(cluster))) + 
    geom_line(size=1.5, aes(color=as.factor(cluster))) +
    geom_ribbon(alpha=0.1) + 
    scale_colour_tableau('Tableau 10') + 
    scale_fill_tableau('Tableau 10') + 
    theme_bw()

filename <- paste0(output_dir, "combined3__param_algo-", param_algo, "__param_k-", param_k, ".png")
ggsave(filename, scale=2, width=3, height=3)

#
# sd
ggplot(data_long2, aes(x=age_group, y=measurement, ymin=(measurement-sd), ymax=(measurement+sd), group=cluster, fill=as.factor(cluster))) + 
    #geom_line(size=1.5, aes(color=as.factor(cluster))) +
    geom_ribbon(alpha=0.4) + 
    scale_colour_tableau('Tableau 10') + 
    scale_fill_tableau('Tableau 10') + 
    theme_bw()

filename <- paste0(output_dir, "combined4__param_algo-", param_algo, "__param_k-", param_k, ".png")
ggsave(filename, scale=2, width=3, height=3)

# sd
ggplot(data_long2, aes(x=age_group, y=measurement, ymin=(measurement-sd), ymax=(measurement+sd), group=cluster, fill=as.factor(cluster))) + 
    geom_smooth(stat = "smooth", formula = y ~ poly(x, 2), aes(color=as.factor(cluster))) +
    scale_colour_tableau('Tableau 10') + 
    scale_fill_tableau('Tableau 10') + 
    theme_bw()

#




# Sandbox - Fake facets ----
# 
library(grid)

plots2 = list()
for (i in 1:length(plots)){
    num_tf_in_cluster = table(tca@cluster)[[i]]
    if(i != length(plots)){
        plots[[i]] = plots[[i]]  + theme(axis.title.x=element_blank(),
                                         axis.text.x=element_blank(),
                                         axis.ticks.x=element_blank())
    }
    plots[[i]] <- plots[[i]] + ggtitle('')
    plots2[[i]] <- plots[[i]]
}


grid.newpage()
g <- grid.draw(rbind(ggplotGrob(plots2[[1]]), 
                     ggplotGrob(plots2[[2]]), 
                     ggplotGrob(plots2[[3]]), 
                     ggplotGrob(plots2[[4]]), 
                     ggplotGrob(plots2[[5]]), 
                     ggplotGrob(plots2[[6]]), 
                     size = "last"))


