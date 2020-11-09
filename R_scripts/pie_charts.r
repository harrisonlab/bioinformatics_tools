# Load library
library(ggplot2)

# All CAZY separated by categories
df <- data.frame(Category = c("AA","CBM","CE","GH","GT","PL"),value = c(108,82,143,314,110,32))
head(df)

# Create a basic bar
pie = ggplot(df, aes(x="", y=value, fill=Category)) + geom_bar(width = 1, size = 1, color = "white", stat = "identity")
# Convert to pie (polar coordinates) and add labels
rlabel = 1.3
pie = pie + coord_polar(theta = "y") + geom_text(aes(x=rlabel, label = paste0(round(100*value/sum(value)), "%")), size = 8, position = position_stack(vjust = 0.5))
# Add color scale (hex colors)
pie = pie + scale_fill_manual(values=c("#55DDE0", "#33658A", "#2F4858", "#F6AE2D", "#F26419", "#999999"))
#pie = pie + scale_fill_manual(values=c("#003f5c", "#444e86", "#955196", "#dd5182", "#ff6e54", "#ffa600"))
#pie = pie + scale_fill_manual(values=c("#2f4858", "#006573", "#00826f", "#3e9a4d", "#9aa915", "#ffa600"))


# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL)
# Tidy up the theme
pie = pie + theme_classic() + theme(axis.line = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    plot.title = element_text(hjust = 0.5, color = "#666666"))


#All secreted CAZY separated by categories
df <- data.frame(Category = c("AA","CBM","CE","GH","GT","PL"),value = c(54,51,49,141,4,30))
head(df)

# Create a basic bar
pie = ggplot(df, aes(x="", y=value, fill=Category)) + geom_bar(width = 1, size = 1, color = "white", stat = "identity")
# Convert to pie (polar coordinates) and add labels
rlabel = 1.3
pie = pie + coord_polar(theta = "y") + geom_text(aes(x=rlabel, label = paste0(round(100*value/sum(value)), "%")), size = 6, position = position_stack(vjust = 0.5))
# Add color scale (hex colors)
pie = pie + scale_fill_manual(values=c("#55DDE0", "#33658A", "#2F4858", "#F6AE2D", "#F26419", "#999999"))
#pie = pie + scale_fill_manual(values=c("#003f5c", "#444e86", "#955196", "#dd5182", "#ff6e54", "#ffa600"))
#pie = pie + scale_fill_manual(values=c("#2f4858", "#006573", "#00826f", "#3e9a4d", "#9aa915", "#ffa600"))


# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL)
# Tidy up the theme
pie = pie + theme_classic() + theme(axis.line = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    plot.title = element_text(hjust = 0.5, color = "#666666"))
