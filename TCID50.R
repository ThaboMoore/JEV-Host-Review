# Load packages
library(dplyr)
library(tidyr)
library(scales)
library(pracma)

# Set directory to where you want plot .pngs to appear
setwd("C:/Users/s5248233/OneDrive - Griffith University/JEV Review/Analysis/Code and figures")

##############################################
### PART 1: DATA CLEANUP AND CONSOLIDATION ###
##############################################

### UNIT: TCID50

### DATA TO DESCRIBE VIRAEMIC PROFILES 

# Load and clean viraemia data
data_raw <- tibble(read.csv(file=file.choose(), header=T)) # load csv
data <- subset(data_raw, is.na(data_raw$omit)) # remove entries marked to omit from analysis
data$group <- with(data, paste(host_common_name)) # add 'group' column - this is a leftover from a previous grouping, but could be used to great new group at anytime. 

# Separate entries in 'days-post-infection' form (data_DPI) from those in summary form (data_sum)
data_sum <- subset(data, data$summarized == 1)
data_DPI <- subset(data, is.na(data$summarized))

# Convert DPI dataset to summarised format for each experimental replicate (ie each individual animal for which we have daily viraemia values)
data_DPI_summary <- data_DPI %>% group_by(ind_ID) %>% # 
  summarise(host_figure_name = unique(host_figure_name),
            host_common_name = unique(host_common_name),
            figure_name = unique(figure_name), 
            scientific = unique(scientific), 
            ind_ID = unique(ind_ID),
            peak_titre = max(titre), # determine each replicate's (ie ind_ID) peak titre across DPI
            duration = sum(titre>0), # determine each replicate's (ie ind_ID) duration of viraemia across DPI
            unit = unique(unit),
            group = unique(group)) %>%
  mutate(sample_size = 1, # for each replicate (ie ind_ID) add a sample size of 1 to be summed later
         n_viraemic = 1)  # for each replicate (ie ind_ID) add a 'number of viraemic' of 1 to be summed later

# Consolidate DPI summary dataset by all unique species
data_DPI_grouped <- data_DPI_summary %>% group_by(group, unit) %>% 
  summarise(host_figure_name = unique(host_figure_name),
            host_common_name = unique(host_common_name),
            figure_name = unique(figure_name), 
            scientific = unique(scientific), 
            avg_peak_titre = mean(peak_titre), # avg. peak titre across replicates
            max_peak_titre = max(peak_titre), # max. peak titre across replicates
            avg_duration = mean(duration), # avg. duration across replicates
            unit = unique(unit),
            group = unique(group),
            sample_size = sum(sample_size), # sample size for that species
            n_viraemic = sum(n_viraemic))   # number viraemic for that species

# Consolidate pre-summarized dataset by all unique species
data_sum_grouped <- data_sum %>%
  group_by(group, unit) %>%
  mutate(flag = ifelse(any(titre > 0 & n_viraemic > 0) & titre == 0 & n_viraemic == 0, FALSE, TRUE)) %>% # create a flag for the non-viraemic species - the weighted mean cannot be calculated for species with values of 0. 
  summarise(host_figure_name = unique(host_figure_name),
            host_common_name = unique(host_common_name),
            figure_name = unique(figure_name), 
            scientific = unique(scientific), 
            avg_peak_titre = ifelse(sum(n_viraemic) == 0, 0, weighted.mean(ifelse(flag, titre, NA), ifelse(flag, n_viraemic, NA), na.rm = TRUE)),
            max_peak_titre = max(max_titre, na.rm=TRUE),
            avg_duration = weighted.mean(duration, total_innoc, na.rm = TRUE),
            unit = unique(unit),
            group = unique(group),
            sample_size = sum(total_innoc),
            n_viraemic = sum(n_viraemic))

# Merge pre-summarized and DPI datasets that have now been consolidated by species
data_final_stack <- rbind(data_sum_grouped, data_DPI_grouped) 

# Final consolidation of the newly merged dataset into single entries for each species
data_final<- data_final_stack %>% group_by(group, unit) %>%
  summarise(host_figure_name = unique(host_figure_name),
            host_common_name = unique(host_common_name),
            figure_name = unique(figure_name), 
            scientific = unique(scientific), 
            avg_peak_titre = weighted.mean(avg_peak_titre, sample_size),
            max_peak_titre = max(max_peak_titre, na.rm=TRUE),
            avg_duration = weighted.mean(avg_duration, sample_size, na.rm = TRUE),
            unit = unique(unit),
            sample_size = sum(sample_size),
            n_viraemic = sum(n_viraemic), 
            prop_viraemic = round(sum(n_viraemic)/sum(sample_size), 2),
            group = unique(group)) 

#######################################################################
### PART 2: QUADRATIC APPROXIMATION, FIG 1 PLOT and AUC CALCULATION ###
#######################################################################

# Create empty data frame to store AUC output from upcoming for-loop
df_AUC_TCID50 <- tibble(group = character (),
                      host=character(), 
                      unit = character(),
                      AUC_avg_titre = numeric())

# Subset data to only include rows where unit is 'TCID50'
data_filtered <- subset(data_final, unit == "TCID50")

# split final dataset into a list of sub-dataframes, one df for each host species
by_host <- split(data_filtered, data_filtered$group)

# Initiate for-loop to cycle through each host species one-by-one and create each viraemia curve
for (h in 1: length(by_host)){
  df <- by_host[[h]] # assign the species sub-dataframe to 'df' for easy calling
  # Set up JPEG plotting interface, wider margines, and empty plot with axes for the species
  jpeg(paste(df$host_figure_name[1],'_TCID50.jpeg', sep=''), units='in', width = 5, height = 4, res=1200, pointsize=6) 
  par(mar = c(5,5,5,2))
  plot(NA, NA, ylim=c(0,max(data_final$max_peak_titre))*1.05, xlim=c(0,8), type='o', pch=16, 
       main=ifelse(df$unit[1]=='TCID50', paste(df$host_figure_name[1], '- TCID50'), df$host_figure_name),
       xlab = 'Mean Duration of Viraemia (days)', ylab = paste('Titre (',df$unit[1],'/mL)',sep=''), 
       cex.main = 1.75, cex.lab = 1.5, cex.axis =1.25)
  color_store <- c() # vector to save colors for figure legend
  
  # Initiate nested for-loop to cycle through each genotype for a given species one-by-one
  for (i in 1:nrow(df)){
    color <- '#0A54A8'
    axis <- seq(0, 8, by = 0.0001) # sets x-values as which to evaluate the quadratic from day 0-8 
    
    # Generate avg. peak titre curve
    x_avg <- c(0,(df$avg_duration[i]+1)/2, df$avg_duration[i]+1) # x-vals for the quadratic model
    y_avg <- c(0, df$avg_peak_titre[i], 0) # y-vals for the quadratic model
    mod_avg <- lm(y_avg ~ poly(x_avg, 2)) # fit the quadratic model
    pred_avg <- predict(mod_avg, newdata = data.frame(x_avg=axis)) # generate fitted values from day 0-8
    pred_avg_plot <- subset(pred_avg, pred_avg >= 0) # keep only positive values for plotting 
    x_avg_plot <- axis[1:length(pred_avg_plot)] # keep only the x-values for the positive y-values
    
    # calculate AUC
    if (df$unit[i]=='TCID50') {# only calculates AUC for TCID50-based curves
      AUC_avg <- trapz(x_avg_plot, pred_avg_plot)
      df_AUC_TCID50 <- df_AUC_TCID50 %>% add_row(group = df$group[i],
                                             host = df$host_figure_name[i], 
                                             unit = df$unit[i],
                                             AUC_avg_titre = AUC_avg)
    }
    # plot polygon for current species
    polygon(x = c(x_avg_plot, rev(x_avg_plot)), 
            y = c(pred_avg_plot, rep(0, length(x_avg_plot))), 
            col = alpha(color, 0.25), border=NA)
    points(x_avg_plot, pred_avg_plot, type='l', lwd=2, col = color) # add lower bound line
    points(x_avg, y_avg, type = 'p', pch=16, col = color,lwd=0.5, cex =1.5) # add lower bound points
    # Add a label
    # x-coordinate of the point
    x_point <- (df$avg_duration[i]+1)/2
    
    # Add a label to the right of the point
    text(x = x_point + 1.5, y = df$avg_peak_titre[i], labels = "Mean peak titre", adj = 0, cex = 1.2, col = color)
    arrows(x0 = (df$avg_duration[i]+1)/2, y0 = df$avg_peak_titre[i], x1 = (df$avg_duration[i]+1)/2 + 1, y1 = df$avg_peak_titre[i], col = color, angle = 30, code = 2, length = 0.1)
    points(x_max_plot, pred_max_plot, type='l', lty=3 , col = color) # add upper bound line
    
    # add samples sizes in upper left corner
    text(x=-0.2, y=max(data_final$max_peak_titre)*(1.05-0.065*(i-1)), adj = 0, 
         labels = paste("sample size: ", df$sample_size[i]), col = color, cex=1.5)
  }
  
  dev.off() # close plotting window
}

###############################################################
### PART 3: SCALE AUC BY PROPORTION VIRAEMIC AND PLOT FIG 2 ###
###############################################################

# Calculate AUC weighted by p_viraemic

df_AUC_scaled_TCID50 <- left_join(df_AUC_TCID50, data_final, by = c('group', 'unit')) # join AUC and p_viraemic data
df_AUC_scaled_TCID50<- df_AUC_scaled_TCID50%>%
  mutate(df_AUC_scaled_TCID50, weighted_AUC = AUC_avg_titre*prop_viraemic) %>% # calculate weighted AUC column
  mutate(df_AUC_scaled_TCID50, figure_name = paste(figure_name, sep = ' - ')) %>%
  arrange(desc(weighted_AUC)) # sort descending

cols <- c('#0A54A8','#F24384') 

######## PROTOTYPE
colfunc<-colorRampPalette(c("#005AB5"))
pal<-colfunc(1000)
cols_avg <- pal[with(df_AUC_scaled_TCID50, round(weighted_AUC/max(weighted_AUC)*1000,0))]
jpeg('TCID50.jpeg', units='in', width = 5, height = 1.4, res=1200, pointsize=6) 
par(mar=c(5,20,3,15))
plot(NA, NA, ylim=c(0, nrow(df_AUC_scaled_TCID50)+1), xlim=c(0, 15), 
     yaxs = 'i', main='',
     xlab = expression('(AUC of viraemia curve)*(proportion viraemic)'),
     yaxt = 'n', ylab='', cex.main = 1.75, cex.lab = 1.2, cex.axis =1.25)
index <- nrow(df_AUC_scaled_TCID50)
for (i in 1:index){
  lines(x = c(-0.4,df_AUC_scaled_TCID50$weighted_AUC[i]),
        y = c(index+(1-i), index+(1-i)), col=alpha(cols_avg[i], 0.5), lwd=3)}
with(df_AUC_scaled_TCID50,{
  point_colors <- ifelse(weighted_AUC == 0, "grey", cols_avg)
  points(x = weighted_AUC, y=index:1, type='p', pch=16, cex = 2, col = point_colors)
})

mtext(side=2, text=df_AUC_scaled_TCID50$scientific, at=index:1, line=8, las=1, cex=0.8, adj=0)
mtext(side=2, text=df_AUC_scaled_TCID50$figure_name, at=index:1, line=19, las=1, cex=1.2, adj=0)
mtext(side=4, text=df_AUC_scaled_TCID50$sample_size, at=index:1, line=3.1, las=1, cex=1.2, adj=0.5)  
mtext(side=4, text=df_AUC_scaled_TCID50$prop_viraemic, at=index:1, line=8.5, las=1, cex=1.2, adj=0.5)  
mtext(side=4, text = 'TCID50', line=-5.0, padj=-1.6, las=1, cex=1.7, font=2)

dev.off()









