bgr1 <- theme(# Plot panel
              panel.grid.major=element_blank(),       
              panel.grid.minor=element_blank(),
              panel.background=element_rect(colour="gray22", 
                                            size=1, fill=NA),
              # Main title
              plot.title=element_text(face="bold", colour="black", 
                                      size=10, hjust=0.5,
                                      margin=margin(b = 15, unit = "pt")),
              # X-axis title 
              axis.title.x=element_text(face="bold", size=20,
                                        margin=margin(t=15, unit = "pt")),
              # Y-axis title
              axis.title.y=element_text(face="bold", size=20,
                                        margin=margin(r=15, unit = "pt")),
              # All text within theme function
              text=element_text(colour="black"),
              # Axis tick labels
              axis.text.x = element_text(size=15, angle=360, colour="black"),
              axis.text.y = element_text(size=15, angle=360, colour="black"),
              # Axis ticks
              axis.ticks = element_line(size=0.7, colour="gray22"),
              axis.ticks.length=unit(.28, "cm"),
              aspect.ratio=1,
              # Legend text
              legend.text=element_text(size=15),
              legend.title=element_text(size=20, face="bold")
)

# Larger size of labels
bgr2 <- bgr1 + 
  theme(
    # X-axis title 
    axis.title.x=element_text(face="bold", size=30,
                              margin=margin(t=15, unit = "pt")),
    # Y-axis title
    axis.title.y=element_text(face="bold", size=30,
                              margin=margin(r=15, unit = "pt")),
    # Axis tick labels
    axis.text.x = element_text(size=25, angle=360, colour="black"),
    axis.text.y = element_text(size=25, angle=360, colour="black"),
    legend.text=element_text(size=20),
    legend.title=element_text(size=25, face="bold")
  )

# Lighter panel.background=element_rect()
bgr2.1 <- bgr2 +
  theme(
    panel.background=element_rect(colour="gray44", 
                                  size=1, fill=NA)
  )
  
bgr3 <- bgr1 + 
  theme(
    # X-axis title 
    axis.title.x=element_text(face="bold", size=40,
                              margin=margin(t=15, unit = "pt")),
    # Y-axis title
    axis.title.y=element_text(face="bold", size=40,
                              margin=margin(r=15, unit = "pt")),
    # Axis tick labels
    axis.text.x = element_text(size=35, angle=360, colour="black"),
    axis.text.y = element_text(size=35, angle=360, colour="black")
    )

bgr4 <- theme(panel.grid.major=element_blank(),       
              panel.grid.minor=element_blank(),
              # Modify border lines of plot
              panel.background=element_rect(colour="gray22", 
                                            size=1, fill=NA),
              # Main title
              plot.title=element_text(face="bold", colour="black", 
                                      size=18, hjust=0.5,
                                      margin=margin(b = 15, unit = "pt")),
              # X-axis title 
              axis.title.x=element_text(face="bold", size=30,
                                        margin=margin(t=15, unit = "pt")),
              # Y-axis title
              axis.title.y=element_text(face="bold", size=30,
                                        margin=margin(r=15, unit = "pt")),
              # All text within theme function
              text=element_text(colour="black"),
              # Axis tick labels
              axis.text.x = element_text(size=25, angle=360, colour="black"),
              axis.text.y = element_text(size=25, angle=360, colour="black"),
              # Axis ticks
              axis.ticks = element_line(size=1, colour="gray22"),
              axis.ticks.length=unit(.20, "cm")
              #,
              #aspect.ratio=1
)

# For the GO-like plots 
bgr5 <- theme(panel.grid.major.y=element_line(colour="gray"),
              panel.grid.minor.y=element_line(colour="gray"),
              panel.grid.major.x=element_line(colour="gray"),
              panel.background=element_rect(colour="gray22", 
                                            size=1, fill="ghostwhite"),
              plot.title=element_text(size=8, colour="black"),
              axis.title.y=element_blank(),
              axis.text.y=element_text(size=10, colour="black", face="bold"),
              axis.title.x=element_text(size=10, colour="black", face="bold"),
              axis.text.x=element_text(size=10, colour="black", face="bold"),
              legend.text=element_text(size=10, face="bold"),
              legend.title=element_text(size=10, face="bold")
)