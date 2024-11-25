################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("ggplot2")
# install.packages(required_pckgs, dependencies = TRUE, Ncpus = detectCores() - 1)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Load in the data
times <- read.csv(file = "Data/Times.csv", header = FALSE)

################################################################################
## Create a data frame for the data
sample_size <- c(250, 500, 1000, 2000, 4000)
d <- c(5, 10, 15)
methods <- c("One-step - Graphical", "Two-step - Graphical", "Three-step - Graphical",
             "One-step - Saturated", "Two-step - Saturated", "Three-step - Saturated")

## This is the order we want to present the methods in for aestetic purposes
methods_order <- c("One-step - Graphical", "One-step - Saturated", 
                   "Two-step - Graphical", "Two-step - Saturated", 
                   "Three-step - Graphical", "Three-step - Saturated")

n_sample_size <- length(sample_size)
n_d <- length(d)
n_methods <- length(methods)

times_table <- data.frame(time = as.numeric(unname(do.call(c, c(times)))),
                          sample_size = rep(sample_size, each = n_d*n_methods),
                          Dimension = rep(rep(d, each = n_methods), times = n_sample_size),
                          Model = rep(methods, times = n_d*n_sample_size),
                          line_no = rep(1:(n_d*n_methods), times = n_sample_size))

## Make factors
times_table$Dimension <- factor(times_table$Dimension, levels = d)
times_table$Model <- factor(times_table$Model, levels = methods_order)
times_table$line_no <- factor(times_table$line_no, levels = 1:(n_d*n_methods))

## Create the ggplot
## guides command sets the order of the legends and manipulates them to be level-ish
pdf("Images/Simulation_Study/Time_Graph.pdf", height = 10, width = 10)
ggplot(data = times_table, aes(x = sample_size, y = log(time), group = line_no, colour = Model)) +
  geom_line(aes(linetype = Dimension), linewidth = 1) +
  geom_point(size = 3) +
  scale_linetype_manual(values = c(1,2,3)) +
  scale_x_continuous(limits = c(250, 4000), breaks = sample_size) +
  labs(x = "Number of Excesses", y = "log(Time) (seconds)") +
  theme(legend.position = "top",
        legend.justification = "centre",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.box = "horizontal") +
  guides(colour = guide_legend(order = 1),
         linetype = guide_legend(order = 2,
                                 theme = theme(legend.margin=margin(15, 0, 0, 0))))
dev.off()
################################################################################

