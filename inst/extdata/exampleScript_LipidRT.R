
################################################################################
# LipidRT 1.0.0:
################################################################################

library(LipidRT)

#==============================================================================#
# Input
#==============================================================================#
# Example data files can be downloaded from
# https://.................................. ????????

setwd("yourWorkingDirectory")


# The input data corresponds to a peaklist table containing lipid annotations,
# which can be provided in .txt or .csv format. This table may originate from:
#
# - LipidMS, such as msbatch$features or msobject[[x]]$annotation$annotatedpeaklist
#
# - MS-DIAL
#
# - Other software tools, as long as the table includes the following required columns:
#     · mz: mass-to-charge ratio
#     · RTmin or RTsec: retention time (in minutes or seconds)
#     · lipidname: lipid annotation for each feature
#     · adduct: ion adduct type
#
# Alternatively, the input can be a data frame already loaded into the R environment,
# or an rtdata object resulting from previous filtering or prediction steps.

#==============================================================================#
# Read the input data
#==============================================================================#
rtdata <- readrtdata(data = "annotatedpeaklist.csv",
                     sep.csv = ",")
#==============================================================================#
# Filtering
#==============================================================================#

# Once you have an annotated peaklist, you can filter the data to retain only
# those features that are correctly annotated based on lipid retention time.
rtdata <- filterannotations(rtdata,
                            isomeric_plasmalogens = FALSE,
                            filedbs = NULL,
                            sep = ",",
                            adducts = NULL,
                            dmz = 5,
                            ppm = TRUE,
                            distC_filter = 5,
                            minSpecies_filter = 5,
                            distribution = "gamma",
                            rtscorethr = 0.5,
                            rttol_id = 5,
                            rttol_name = 1,
                            filter_plot = TRUE,
                            verbose = TRUE)

# Print annotations filtering plots to pdf
pdf("Filtering.pdf", width = 14, height = 7) # You can change the size
for(p in 1:length(rtdata$filtering$plots)){
  print(rtdata$filtering$plots[[p]])
}
dev.off()

#==============================================================================#
# Prediction
#==============================================================================#

# Once the annotated peaklist has been filtered, you can proceed to predict new
# lipid annotations based on retention time for features that lacked sufficient
# MS/MS information for confident identification.

# You can also pass the annotated peaklist directly to the prediction function
# without prior filtering. In this case, the filtering will be performed
# internally using the default parameters before the prediction step.
rtdata <- rtpredictannotations(rtdata,
                               distribution = "gamma",
                               subclasses = NULL,
                               filedbs = NULL,
                               sep = ",",
                               adducts = NULL,
                               dmz = 5,
                               ppm = TRUE,
                               rtscorethr = 0.5,
                               distC_predict = 5,
                               minSpecies_predict = 5,
                               rt_thr = 30,
                               rttol_id = 5,
                               rttol_name = 1,
                               predict_plot = TRUE,
                               verbose = TRUE)

# Print RT predictions plots to pdf
pdf("Predictions.pdf", width = 16, height = 7) # You can change the size
for(p in 1:length(rtdata$prediction$plots)){
  print(rtdata$prediction$plots[[p]])
}
dev.off()



