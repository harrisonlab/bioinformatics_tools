# Absorbance

library(reshape2)
library(dplyr)
library(ggplot2)

## No correction

# Rows with samples
rawdata <- read.csv("raw.csv")
# Sample, Time, OD650 columns
reshaped <- melt(rawdata, id=c("Sample", "Well"), variable.name="Well", value.name="OD600")
# Write csv file
#annotated <- inner_join(dd, platemap, by="Well")
write.csv(reshaped, "data-annotated.csv")
dd <- read.csv("data-annotated.csv")
# Group samples by time
#grouped <- group_by(dd, Sample, Time)
# Calculate statistics for each group
#stats <- summarise(grouped, N=length(OD600), Average=mean(OD600), StDev=sd(OD600))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(data) {
    n <- length(data)
    error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits in one line
stats <- dd %>%
              group_by(Sample, Time) %>%
              summarise(N=length(OD600),
                        Average=mean(OD600),
                        CI95=conf_int95(OD600)) %>%
     
# Plot results. Average line
ggplot(data=stats, aes(x=Time, y=Average, color=Sample)) + geom_line() + labs(x="Time (Hours)", y="Absorbance at 650 nm")
# Plot results. Average line separated
ggplot(data=stats, aes(x=Time, y=Average, color=Sample)) + geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=Sample),color=NA, alpha=0.3) + 
       geom_line() +
       scale_y_log10() +
       facet_grid(Sample ~ .) +
       labs(x="Time (Hours)", y="Absorbance at 650 nm")
# Plot results.
ggplot(data=stats, aes(x=Time, y=Average, color=Sample)) +
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=Sample),
color=NA, alpha=0.3) + 
geom_line() +
labs(x="Time (Hours)", y="Absorbance at 650 nm") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

# Plot results using log scale
ggplot(data=stats, aes(x=Time, y=Average, color=Sample)) +
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=Sample),
color=NA, alpha=0.3) + 
scale_y_log10() +
geom_line() +
labs(x="Time (Hours)", y="Absorbance at 650 nm") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

## Corrected with blank results

# Rows with samples
A <- read.csv("Corrected.csv")
# Sample, Time, OD650 columns
reshaped <- melt(A, id=c("Sample"), variable.name="Time", value.name="OD650")
# Write csv file
write.csv(reshaped, "Corrected2.csv")
B <- read.csv("Corrected2.csv")
# Group samples by time
grouped <- group_by(B, Sample,Time)
# Calculate statistics for each group
stats <- summarise(grouped, N=length(OD600), Average=mean(OD600), StDev=sd(OD600))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(data) {
    n <- length(data)
    error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits
stats <- summarise(grouped, N=length(OD600), Average=mean(OD600), CI95=conf_int95(OD600))
# Plot results 
ggplot(data=stats, aes(x=Time, y=Average, color=Sample)) +
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=Sample),
color=NA, alpha=0.3) + 
geom_line() +
scale_fill_manual(values=c('#228B22','#E69F00')) +
scale_color_manual(values=c('#228B22','#E69F00')) +
labs(x="Time (Hours)", y="Absorbance at 650 nm") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)
# Plot results using log scale
ggplot(data=stats, aes(x=Time, y=Average, color=Sample)) +
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=Sample),
color=NA, alpha=0.3) + 
scale_y_log10() +
geom_line() +
scale_fill_manual(values=c('#228B22','#E69F00')) +
scale_color_manual(values=c('#228B22','#E69F00')) +
labs(x="Time (Hours)", y="Absorbance at 650 nm") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

## Remove blanks

# Rows with samples
A <- read.csv("raw_noblank.csv")
# Sample, Time, OD650 columns
reshaped <- melt(A, id=c("Sample"), variable.name="Time", value.name="OD650")
# Write csv file
write.csv(reshaped, "raw2.csv")
B <- read.csv("raw2.csv")
# Group samples by time
grouped <- group_by(B, Sample,Time)
# Calculate statistics for each group
stats <- summarise(grouped, N=length(OD650), Average=mean(OD650), StDev=sd(OD650))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(data) {
    n <- length(data)
    error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits
stats <- summarise(grouped, N=length(OD650), Average=mean(OD650), CI95=conf_int95(OD650))
# Plot results 
ggplot(data=stats, aes(x=Time, y=Average, color=Sample)) +
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=Sample),
color=NA, alpha=0.3) + 
geom_line() +
scale_fill_manual(values=c('#228B22','#E69F00')) +
scale_color_manual(values=c('#228B22','#E69F00')) +
labs(x="Time (Hours)", y="Absorbance at 650 nm") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)
# Plot results using log scale
ggplot(data=stats, aes(x=Time, y=Average, color=Sample)) +
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=Sample),
color=NA, alpha=0.3) + 
scale_y_log10() +
geom_line() +
scale_fill_manual(values=c('#228B22','#E69F00')) +
scale_color_manual(values=c('#228B22','#E69F00')) +
labs(x="Time (Hours)", y="Absorbance at 650 nm") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)

# Fluorescence

# Rows with samples
A <- read.csv("Data.csv")
# Sample, Time, OD650 columns
reshaped <- melt(A, id=c("Time"), variable.name="Time", value.name="OD650")
# Write csv file
write.csv(reshaped, "raw2.csv")
B <- read.csv("raw2.csv")
# Group samples by time
grouped <- group_by(B, Sample,Time)
# Calculate statistics for each group
stats <- summarise(grouped, N=length(OD650), Average=mean(OD650), StDev=sd(OD650))
# Create a function that calculates 95% confidence intervals for the given data vector using a t-distribution
conf_int95 <- function(data) {
    n <- length(data)
    error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
    return(error)
}
# Create summary for each group containing sample size, average OD600, and 95% confidence limits
stats <- summarise(grouped, N=length(OD650), Average=mean(OD650), CI95=conf_int95(OD650))
# Plot results 
ggplot(data=stats, aes(x=Time, y=Average, color=Sample)) +
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=Sample),
color=NA, alpha=0.3) + 
geom_line() +
scale_fill_manual(values=c('#228B22','#E69F00')) +
scale_color_manual(values=c('#228B22','#E69F00')) +
labs(x="Time (Hours)", y="Absorbance at 470-15/515-20") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)
# Plot results using log scale
ggplot(data=stats, aes(x=Time, y=Average, color=Sample)) +
geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, fill=Sample),
color=NA, alpha=0.3) + 
scale_y_log10() +
geom_line() +
scale_fill_manual(values=c('#228B22','#E69F00')) +
scale_color_manual(values=c('#228B22','#E69F00')) +
labs(x="Time (Hours)", y="Absorbance at 470-15/515-20") +
theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey90"),panel.ontop = F)