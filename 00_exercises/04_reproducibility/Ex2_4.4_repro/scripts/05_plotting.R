# this code generates the visualizations in the manuscript
# --------------------------------------------------------------

# create figure:
# build the matrix from the results
tot.mat <- cbind(100*scenarios,apply(beta.hat,2,mean))

# set column names
colnames(tot.mat) <- c("me.exp","me.conf","estimate")

# plotting: X and Y axis represent % of variance due to measurement error in the two traits
FIGURE <- ggplot(tot.mat, aes(me.exp, me.conf)) +

  # heatmap with color scale
  geom_tile(color="white",aes(fill = estimate)) +
  
  # numerical values in each cell
  geom_text(aes(label = round(estimate, 2))) +
  
  # custom scale centered on the reference value (true value)
  scale_fill_gradient2(low="#D55E00",mid="white",high = "#56B4E9", midpoint=ref) +

  # axis labels
  labs(x=paste("% of total variance of HbA1c due to measurement error"),
       y=paste("% of total variance of BMI due to measurement error")) +
  
  # equal aspect ratio to ensure we get square cells
  coord_equal()+
  
  # ticks only at the scenario levels
  scale_y_continuous(breaks=unique(tot.mat[,1]))+
  scale_x_continuous(breaks=unique(tot.mat[,1]))+

  # custom theme
  theme(panel.background = element_rect(fill='white', colour='grey'),
        plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.title=element_text(size=8),
        axis.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10))

# saving the plot using the here package
ggsave(
  filename = here("output", "Figure_STRATOS.tif"),
  plot = FIGURE,
  width = 6,
  height = 6,
  dpi = 300
)