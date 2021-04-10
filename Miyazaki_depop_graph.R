# Script for cleaning/plotting Miyazaki depopulation
    
library(tidyverse)
    
# Clean data from Miyazaki_depop_data.csv
# Change date to yyyy/mm/dd
clean_data <- function(infile = "./Miyazaki_depop_data.csv", 
                       outfile = "./Miyazaki_depop_data_edited.csv") {
    
    df_orig <- readr::read_csv(infile) %>%
        tidyr::separate(date, c("mo", "day", "other"), "[^0-9]") %>% 
        dplyr::select(-c(other)) 
    
    year <- as.character(rep(2010, dim(df_orig)[1]))
    df <- df_orig %>% 
        dplyr::mutate(year) %>% 
        tidyr::unite(col = "date", c("year", "mo", "day"), sep = "/")
    
    readr::write_csv(df, outfile)
}

# Plot Miyazaki depopulation graph
# Save to file name
epi_depop_plot <- function(infile = "./Miyazaki_depop_data_edited.csv") {
    # Setting system language to English for plotting purpose
    Sys.setlocale("LC_TIME", "C")  
    df_plt <- readr::read_csv(infile)
    plt <- ggplot(df_plt, aes(x = date)) +
        geom_line(aes(y = reported_infection), size = 1.0) +
        geom_line(aes(y = depop_start), color = "blue", size = 1.0) +
        geom_line(aes(y = depop_completed), color = "green", 
                  linetype = "dotted", size = 1.0) +
        scale_x_date(date_breaks = "1 week", date_labels ="%m/%d") +
        xlab("Date") + ylab("Number of Farms") +
        theme_bw() +
        theme(axis.text.x = element_text(size = 24, angle = 30, vjust = .5),
              axis.text.y = element_text(size = 24),
              axis.title = element_text(size = 24, face = "bold"),
              text = element_text(size = 24),
              legend.text = element_text(size = 24),
              legend.title = element_blank(),
              legend.position = c(0.9, 0.75),
              legend.key.height = unit(2, "lines"),
              plot.margin = unit(c(5, 10, 1, 1), "mm"))
    ggsave(filename = "./figs/epidemic_depopulation_curve.png", plot = plt,
           dpi = 300, width = 10.0, height = 8.0)
}


clean_data(infile = "./Miyazaki_depop_data.csv", 
           outfile = "./Miyazaki_depop_data_edited.csv")
epi_depop_plot(infile = "./Miyazaki_depop_data_edited.csv")
