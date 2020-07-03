library(dplyr)
library(reshape2)
library(ggplot2)
library(growthcurver)
library(purrr)

df <- read.csv("Data.csv", header = TRUE, dec = ".")

gc_out <- SummarizeGrowthByPlate(df)

#Basic R plotting function
gc_out <- SummarizeGrowthByPlate(df, plot_fit = TRUE,
                                 plot_file = "gc_plots.pdf")

#Plotting using ggplot2

summG <- function(x) {SummarizeGrowth(df$time,x)}
lapply(df[2:ncol(df)], summG)

models.all <- lapply(df[2:ncol(df)], function(x) SummarizeGrowth(df$time, x))

df.predicted.plate <- data.frame(time = df$time)
for (i in names(df[2:ncol(df)]))
  {df.predicted.plate[[i]] <- predict(models.all[[i]]$model)}

melt1 <- melt(df, id.vars = "time", variable.name = "sample", value.name = "od")
melt2 <- melt(df.predicted.plate, id.vars = "time", variable.name = "sample", value.name = "pred.od")
df.final <- cbind(melt1, pred.od=melt2[,3])

#ggplot(df.final, aes(x=time, y=od)) + geom_point(aes(), alpha=0.5) + geom_line(aes(y=pred.od), color="red") + facet_wrap(~sample, ncol = 12) + theme_bw()

pdf("curve.pdf", height = 11.69, width = 16.53)
ggplot(df.final, aes(x=time, y=od)) + geom_point(aes(), alpha=0.5) + geom_line(aes(y=pred.od), color="red") + facet_wrap(~sample, ncol = 12) + theme_classic() +
    theme( axis.text = element_text( size = 14 ),
           axis.text.x = element_text( size = 20 ),
           axis.title = element_text( size = 16, face = "bold" ),
           legend.position="none",
           # The new stuff
           strip.text = element_text(size = 20))
dev.off()


# As in the simple example, load the package and the data.
library(growthcurver)
d <- read.csv("Data.csv", header = TRUE, dec = ".")

# Let's create an output data frame to store the results in.
# We'll create it so that it is the right size (it's faster this way!),
# but leave it empty.
num_analyses <- length(names(d)) - 1
d_gc <- data.frame(sample = character(num_analyses),
                   k = numeric(num_analyses),
                   n0  = numeric(num_analyses),
                   r = numeric(num_analyses),
                   t_mid = numeric(num_analyses),
                   t_gen = numeric(num_analyses),
                   auc_l = numeric(num_analyses),
                   auc_e = numeric(num_analyses),
                   sigma = numeric(num_analyses),
                   stringsAsFactors = FALSE)

# Truncate or trim the input data to observations occuring in the first 20 hours.
# Remember that the times in these sample data are reported in hours. To use
# minutes (or to trim at a different time), change the next line of code.
# For example, if you still would like to trim at 20 hours, but your time data
# are reported in minutes use: trim_at_time <- 20 * 60
trim_at_time <- 96

# Now, loop through all of the columns in the data frame. For each column,
# run Growthcurver, save the most useful metrics in the output data frame,
# and make a plot of all the growth curve data and their best fits.

# First, create a plot for each of the wells in the 96-well plate.
# Uncomment the next line to save the plots from your 96-well plate to a
# pdf file in the working directory.
pdf("test.pdf", height = 8.5, width = 11)
par(mfcol = c(8,12))
par(mar = c(0.25,0.25,0.25,0.25))
y_lim_max <- max(d[,setdiff(names(d), "time")]) - min(d[,setdiff(names(d), "time")])

n <- 1    # keeps track of the current row in the output data frame
for (col_name in names(d)) {

  # Don't process the column called "time".
  # It contains time and not absorbance data.
  if (col_name != "time") {

    # Create a temporary data frame that contains just the time and current col
    d_loop <- d[, c("time", col_name)]

    # Do the background correction.
    # Background correction option 1: subtract the minimum value in a column
    #                                 from all measurements in that column
        min_value <- min(d_loop[, col_name])
    d_loop[, col_name] <- d_loop[, col_name] - min_value
    # Background correction option 2: subtract the mean value of blank wells
    #                                 over the course the experiment
    #                                 (Replace B2, D8, G11 with the column
    #                                  names of your media-only wells)
    #d$blank <- apply(d[, c("B2", "D8", "G11")], 1, mean)
    #d$A1 <- d$A1 - d$blank

    # Now, call Growthcurver to calculate the metrics using SummarizeGrowth
    gc_fit <- SummarizeGrowth(data_t = d_loop[, "time"],
                              data_n = d_loop[, col_name],
                              t_trim = trim_at_time,
                              bg_correct = "min")

    # Now, add the metrics from this column to the next row (n) in the
    # output data frame, and increment the row counter (n)
    d_gc$sample[n] <- col_name
    d_gc[n, 2:9] <- c(gc_fit$vals$k,
                      gc_fit$vals$n0,
                      gc_fit$vals$r,
                      gc_fit$vals$t_mid,
                      gc_fit$vals$t_gen,
                      gc_fit$vals$auc_l,
                      gc_fit$vals$auc_e,
                      gc_fit$vals$sigma)
    n <- n + 1

    # Finally, plot the raw data and the fitted curve
    # Here, I'll just print some of the data points to keep the file size smaller
    n_obs <- length(gc_fit$data$t)
    idx_to_plot <- 1:20 / 20 * n_obs
    plot(gc_fit$data$t[idx_to_plot], gc_fit$data$N[idx_to_plot],
         pch = 20,
         xlim = c(0, trim_at_time),
         ylim = c(0, y_lim_max),
         cex = 0.6, xaxt = "n", yaxt = "n")
     text(x = trim_at_time / 4, y = y_lim_max, labels = col_name, pos = 1)
     lines(gc_fit$data$t, predict(gc_fit$model), col = "red")
  }
}

dev.off()


# Check if Growthcurver provided any notes in a plate of growthcurves returned
# from SummarizeGrowthByPlate
gc_out %>% filter(note != "")

# Check if Growthcurver provided any notes in a single growthcurve returned
# from SummarizeGrowth
gc_fit$vals$note

# Load dplyr and the sample output data
library(dplyr)
gc_out <- as_data_frame(gc_out)

# Plot a histogram of the sigma values in order to check for outliers
hist(gc_out$sigma, main = "Histogram of sigma values", xlab = "sigma")

# Load dplyr, ggplot2, and the sample data
library(dplyr)
library(ggplot2)
pca_gc_out <- as_data_frame(gc_out)

# Prepare the gc_out data for the PCA
rownames(pca_gc_out) <- pca_gc_out$sample
## Warning: Setting row names on a tibble is deprecated.
# Do the PCA
pca.res <- prcomp(pca_gc_out %>% select(k:sigma), center=TRUE, scale=TRUE)

# Plot the results
as_data_frame(list(PC1=pca.res$x[,1],
                   PC2=pca.res$x[,2],
                   samples = rownames(pca.res$x))) %>%
  ggplot(aes(x=PC1,y=PC2, label=samples)) +
  geom_text(size = 3)
