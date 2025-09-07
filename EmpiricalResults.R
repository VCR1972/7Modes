########################################
####### Code generating 7 modes
####### using data from the HMD 
####### Authors: Paola Vazquez-Castillo
#######          Vladimir Canudas-Romo
####### Last update: 07 September 2025
########################################

rm(list = ls())

# Load required libraries
library(ggplot2)
library(tidyverse)
library(viridis)
library(HMDHFDplus)
library(ggrepel)

# Please enter your Human Mortality Database username and password below. No other changes to the program are required.

username.hmd = ""
password.hmd = ""

# Please note that this program will take a few minutes to run.

# Country definitions
Names <- c("AUS","AUT","BGR","BLR","CAN","CHE","CHL","CZE",
           "DEUTE","DEUTW","DEUTNP","DNK","GBR_NIR","GBR_NP",
           "GBR_SCO","GBRTENW","ESP","EST","FIN","FRATNP","HUN",
           "IRL","ISL","ISR","ITA","JPN","LTU","LUX","LVA","NLD",
           "NOR","NZL_NM","NZL_NP","POL","PRT","RUS","SVK","SVN",
           "SWE","TWN","UKR","USA")

NamesEnglish <- c("Australia","Austria","Bulgaria","Belarus",
                  "Canada","Switzerland","Chile","Czechia",
                  "Germany, East","Germany, West","Germany",
                  "Denmark","Northern Ireland","United Kingdom",
                  "Scotland","England and Wales","Spain",
                  "Estonia","Finland","France","Hungary",
                  "Ireland","Iceland","Israel","Italy","Japan",
                  "Lithuania","Luxembourg","Latvia",
                  "The Netherlands","Norway",
                  "New Zealand – non-Māori","New Zealand",
                  "Poland","Portugal","Russia","Slovakia",
                  "Slovenia","Sweden","Taiwan","Ukraine",
                  "United States of America")

# Set working directory
setwd("")

# Select year for distribution analysis
# Change this to any year you want to see the distributions (dx, fx and Cal-d)
selected_year <- 2022  

# Countries with cohort data (indices)
# CHE (6), DNK (12), GBRSCO(15), GBRENW(16), FIN(19)
# FRANP(20), ISL(23), ITA(25), NLD(30), NOR(31), SWE(39)

# Select countries to process 
# (Main analysis: CHE, DNK, Scotland, Netherlands)
countries_to_process <- c(12, 30, 15, 6)  
# for all cohort countries:
#countries_to_process <- c(6, 12, 15, 16, 19, 20, 23, 25, 30, 31, 39) 

# Creating empty elements for storing selected year data
# Store dx data for all countries
pd_year_data <- list()  
# Store dx data for processed (from the defined vector) countries
pd_year_processed <- list() 

# Store fx data for all countries
pdpe_year_data <- list() 
# Store fx data for processed countries
pdpe_year_processed <- list()  

# Store Cald data for all countries 
cald_year_data <- list()  
# Store Cald data for processed countries
cald_year_processed <- list()  

# Mode calculation function
Mode <- function(dx, Y){
  dx <- cbind(dx, seq(0, 110))
  
  MaxV <- c() 
  MaxA <- c()
  Mode <- c()
  
  for (y in 1:Y){
    MaxV[y] <- max(dx[10:110, y])
    A <- c()
    A <- ifelse(MaxV[y] == dx[, y], dx[, (Y+1)], 0)
    MaxA[y] <- max(A)
    # Using Kannisto's mode calibration to get decimal point
    Mode[y] <- max(A) + 
      ((dx[(MaxA[y]+1), y] - dx[(MaxA[y]), y]) /
         (2*dx[(MaxA[y]+1), y] - dx[(MaxA[y]), y] - 
            dx[(MaxA[y]+2), y]))
  }
  
  return(Mode)
}

# CAL function for cross-average length of life
CALfunction <- function(Mx){
  lx <- c(1)
  
  for (x in 1:111){
    px <- c()
    
    for (z in 1:x){
      px <- c(px, Mx[z, 111-x+z])
    }
    
    lx <- c(lx, prod(px))
  }
  return(lx[-112])
}

#####################################
# FUNCTIONS FOR TIMING MODE ANALYSIS
#####################################

# Function to find the year when the period mode shifts from infant to adult
THEyear <- function(M){
  dx <- matrix(M$dx, 111)
  N <- c()
  
  for (x in 1:(dim(dx)[2])){
    N <- c(N, min(which(dx[,x] == max(dx[,x]))))
  } 
  
  Year <- c(range(M$Year)[1]:range(M$Year)[2])
  Y <- min(Year[((N-1) > 1)])
  
  return(Y)
}

# Function to find the period modal age level at shifting
THElevel <- function(M){
  dx <- matrix(M$dx, 111)
  N <- c()
  
  for (x in 1:(dim(dx)[2])){
    N <- c(N, min(which(dx[,x] == max(dx[,x]))))
  } 
  
  Year <- c(range(M$Year)[1]:range(M$Year)[2])
  Y <- min(Year[((N-1) > 1)])
  
  L <- M[M$Year == Y,]$dx
  
  A <- min(which(L == max(L))) - 1
  
  # Using Kannisto's mode calibration to get decimal point
  Mode <- A + ((L[(A+1)] - L[A]) / (2*L[(A+1)] - L[(A)] - L[(A+2)]))
  
  return(Mode)
}

# Function to find the year when M-dagger shifts from infant to adult
THEyear_fx <- function(M){
  dx <- matrix(M$dx, 111)
  ex <- matrix(M$ex, 111)
  dx <- dx * ex / 100000
  
  N <- c()
  
  for (x in 1:(dim(dx)[2])){
    N <- c(N, min(which(dx[,x] == max(dx[,x]))))
  } 
  
  Year <- c(range(M$Year)[1]:range(M$Year)[2])
  Y <- min(Year[((N-1) > 1)])
  
  return(Y)
}

# Function to find the M-dagger level at shifting
THElevel_fx <- function(M){
  dx <- matrix(M$dx, 111)
  ex <- matrix(M$ex, 111)
  dx <- (dx * ex) / 100000
  
  N <- c()
  
  for (x in 1:(dim(dx)[2])){
    N <- c(N, min(which(dx[,x] == max(dx[,x]))))
  } 
  
  Year <- c(range(M$Year)[1]:range(M$Year)[2])
  Y <- min(Year[((N-1) > 1)])
  
  L <- (M[M$Year == Y,]$dx) * (M[M$Year == Y,]$ex)
  
  A <- min(which(L == max(L))) - 1
  
  # Using Kannisto's mode calibration to get decimal point
  Mode <- A + ((L[(A+1)] - L[A]) / (2*L[(A+1)] - L[(A)] - L[(A+2)]))
  
  return(Mode)
}

# Function to process a single country
process_country <- function(t) {
  cat("Processing country:", Names[t], "(" , NamesEnglish[t], ")\n")
  
  # Loading data
  # Period Life Table (both sex)
  P <- readHMDweb(Names[t], "bltper_1x1",username.hmd, password.hmd)
  # Cohort life table (both sex)
  C <- readHMDweb(Names[t], "bltcoh_1x1",username.hmd, password.hmd)
  # Actual death counts (both sex)
  A <- readHMDweb(Names[t], "Deaths_1x1",username.hmd, password.hmd)
  
  # Create matrices for period data
  Py <- matrix(P$Year, 111)
  Pd <- matrix(P$dx, 111)/100000
  Pp <- matrix(1-P$qx, 111)
  Pe <- matrix(P$ex, 111)
  
  # Extract year vectors
  PYear <- Py[1,]
  PY <- length(PYear)
  
  # Extract data for selected year and store
  if(selected_year %in% PYear) {
    year_index <- which(PYear == selected_year)
    pd_year_country <- Pd[, year_index]
    pdpe_year_country <- (Pd * Pe)[, year_index]
    
    # Store Pd (period dx) data for processed countries
    pd_year_processed[[Names[t]]] <<- data.frame(
      Age = 0:110,
      Value = pd_year_country,
      Country = NamesEnglish[t],
      CountryCode = Names[t],
      Type = "Processed"
    )
    
    # Store fx (Pd*Pe) data for processed countries
    pdpe_year_processed[[Names[t]]] <<- data.frame(
      Age = 0:110,
      Value = pdpe_year_country,
      Country = NamesEnglish[t],
      CountryCode = Names[t],
      Type = "Processed"
    )
    
    cat("  Stored", selected_year, "Pd and Pd*Pe data for", NamesEnglish[t], "\n")
  }
  
  # Create matrices for cohort data
  Cy <- matrix(C$Year, 111)
  Cd <- matrix(C$dx, 111)/100000
  Ce <- as.numeric(as.character(C$ex))
  Ce[is.na(Ce)] <- 0
  Ce <- matrix(Ce, 111)
  
  # Create matrices for actual death counts
  Ay <- matrix(A$Year, 111)
  AY <- dim(Ay)[2]
  Ad <- matrix(A$Total, 111)/t(matrix(rep(colSums(matrix(A$Total, 111)), 111), AY))
  
  # Extract year vectors
  CYear <- Cy[1,]
  CY <- length(CYear)
  AYear <- Ay[1,]
  AY <- length(AYear)
  
  # Calculate actual death counts for cohorts
  AYearc <- min(AYear):(max(AYear)-90)
  AYc <- length(AYearc)
  Adc <- matrix(0, 111, AYc)
  AAd <- matrix(A$Total, 111)
  AAd <- cbind(AAd, matrix(0, 111, 21))
  
  for (x in 1:AYc){
    for (z in 1:110){
      Adc[z,x] <- AAd[z, x+z-1]
    }
  }
  Adc <- Adc/t(matrix(rep(colSums(Adc), 111), AYc))
  
  # Calculate CAL data
  CALlx <- c()
  CALY <- (min(PYear)+111):max(PYear)
  
  for (x in CALY){
    y <- x - (min(PYear))
    Px <- Pp[, (y-110):y]
    CALlx <- cbind(CALlx, CALfunction(Px))
  }
  
  CALd <- rbind(CALlx[-111,] - CALlx[-1,], 0)
  
  # Cubic spline smoothing of CAL data
  Cald <- matrix(0, 111, dim(CALd)[2])
  for (z in 1:length(CALY)){
    Cald[, z] <- predict(smooth.spline(0:110, CALd[, z], df=12))$y
  }
  
  # extract CAL daata for selected year and store
  if(selected_year %in% CALY) {
    cal_year_index <- which(CALY == selected_year)
    cald_year_country <- Cald[, cal_year_index]
    
    # Store Cald data for processed countries
    cald_year_processed[[Names[t]]] <<- data.frame(
      Age = 0:110,
      Value = cald_year_country,
      Country = NamesEnglish[t],
      CountryCode = Names[t],
      Type = "Processed"
    )
    
    cat("  Stored", selected_year, "Cald data for", NamesEnglish[t], "\n")
  }
  
  # Calculate all modes
  Mp <- Mode(Pd, PY)                    # Mode period
  Mc <- Mode(Cd, CY)                    # Mode cohort
  Ma <- Mode(Ad, AY)                    # Mode death counts period
  Mac <- Mode(Adc, AYc)                 # Mode death counts cohort
  Md <- Mode(Pd*Pe, PY)                 # Mode deaths*le period
  Mdc <- Mode(Cd*Ce, CY)                # Mode cohort deaths*le
  Mcal <- Mode(Cald, length(CALY))      # Mode CAL
  
  # Set up ggplot2 visualization
  colors <- viridis::viridis(7)
  labels <- c("Mp"=expression(M[p]), 
              "Mc"=expression(M[c]),
              "Md"=expression(M[p]^"\u2020"),
              "Mdc"=expression(M[c]^"\u2020"),
              "Mcal"=expression(M^"CAL"),
              "Ma"=expression(M[p]^"A"),
              "Mac"=expression(M[c]^"A"))
  
  colors2 <- c("Mp"=colors[1],
               "Mc"=colors[3],
               "Md"=colors[5],
               "Mdc"=colors[6],
               "Mcal"=colors[4],
               "Ma"=colors[7],
               "Mac"=colors[2])
  
  # Create data frames for plotting
  period_data <- data.frame(Years=PYear, Mp=Mp, Ma=Ma, Md=Md)
  cohort_data <- data.frame(Years=CYear, Mc=Mc, Mac=Mac[-length(Mac)], Mdc=Mdc)
  cal_data <- data.frame(Years=CALY, Mcal=Mcal)
  
  # Join data frames
  one <- left_join(period_data, cohort_data, by="Years")
  wide <- left_join(one, cal_data, by="Years")
  
  # Add padding for early years
  aux <- data.frame(Years=1835:(PYear[1]-1), Mp=NA, Ma=NA, Md=NA, Mc=NA, Mac=NA, Mdc=NA, Mcal=NA)
  wide2 <- rbind(aux, wide)
  
  # Convert to long format
  long <- gather(wide2, key="Mode", value="Value", -Years)
  long$Mode <- factor(long$Mode, levels = c("Mp", "Mc", "Md", "Mdc", "Mcal", "Ma", "Mac")) 
  
  # Create the plot
  plot_country <- ggplot(data=long,
                         aes(Years, y=Value)) +
    geom_point(aes(col=Mode, shape=Mode)) +
    geom_smooth(aes(col=Mode), se=F) +
    scale_color_viridis_d(labels=labels) +
    scale_shape_manual(values=c(seq(15,18),3,4,8),
                       labels=labels) +
    scale_color_manual(values = colors2,
                       labels=labels) +
    ylim(c(60,90)) +
    xlab("Year") +
    ylab("Modal Age") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.background = element_blank(),
          legend.key = element_blank()) +
    guides(color=guide_legend(nrow=1)) +
    ggtitle(NamesEnglish[t]) +
    labs(color="Modes:", shape="Modes:") +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          axis.line=element_line(size=0.5)) +
    theme(
      text = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      panel.grid.major = element_line(linewidth = 0.8),
      panel.grid.minor = element_line(linewidth = 0.6),
      legend.key.size = unit(1,"cm")
    ) +
    guides(color=guide_legend(nrow = 1, override.aes = list(size = 4)),
           shape=guide_legend(nrow=1, override.aes = list(size = 4)))
  
  # Display plot
  print(plot_country)
  
  # Save plot
  # filename <- paste0(Names[t], "_plot4.pdf")
  # ggsave(filename, plot = plot_country, device = "pdf", width = 9, height = 8)
  return(plot_country)
}

# FUNCTION TO COLLECT ALL SELECTED YEAR DATA FOR ALL COUNTRIES
collect_all_year_data <- function() {
  cat("Collecting", selected_year, "data for all countries...\n")

  for(i in 1:length(Names)) {
    # Skip if already processed
    if(Names[i] %in% names(pd_year_processed)) {
      cat("Skipping", Names[i], "(already processed)\n")
      next
    }
    
    cat("Collecting", selected_year, "data for:", Names[i], "\n")
    
    tryCatch({
      # Get period life table data
      P <- readHMDweb(Names[i], "bltper_1x1", username.hmd, password.hmd)
      
      # Create matrices for period data
      Py <- matrix(P$Year, 111)
      Pd <- matrix(P$dx, 111)/100000
      Pp <- matrix(1-P$qx, 111)
      Pe <- matrix(P$ex, 111)
      
      # Extract year vector
      PYear <- Py[1,]
      
      # Check if selected year data exists
      if(selected_year %in% PYear) {
        year_index <- which(PYear == selected_year)
        pd_year_country <- Pd[, year_index]
        pdpe_year_country <- (Pd * Pe)[, year_index]
        
        # Store Pd data
        pd_year_data[[Names[i]]] <<- data.frame(
          Age = 0:110,
          Value = pd_year_country,
          Country = NamesEnglish[i],
          CountryCode = Names[i],
          Type = "Background"
        )
        
        # Store Pd*Pe data
        pdpe_year_data[[Names[i]]] <<- data.frame(
          Age = 0:110,
          Value = pdpe_year_country,
          Country = NamesEnglish[i],
          CountryCode = Names[i],
          Type = "Background"
        )
        
        # Calculate CAL data
        CALlx <- c()
        CALY <- (min(PYear)+111):max(PYear)
        
        for (x in CALY){
          y <- x - (min(PYear))
          Px <- Pp[, (y-110):y]
          CALlx <- cbind(CALlx, CALfunction(Px))
        }
        
        CALd <- rbind(CALlx[-111,] - CALlx[-1,], 0)
        
        # Smooth CAL data
        Cald <- matrix(0, 111, dim(CALd)[2])
        for (z in 1:length(CALY)){
          Cald[, z] <- predict(smooth.spline(0:110, CALd[, z], df=12))$y
        }
        
        # Store Cald data if selected year exists
        if(selected_year %in% CALY) {
          cal_year_index <- which(CALY == selected_year)
          cald_year_country <- Cald[, cal_year_index]
          
          cald_year_data[[Names[i]]] <<- data.frame(
            Age = 0:110,
            Value = cald_year_country,
            Country = NamesEnglish[i],
            CountryCode = Names[i],
            Type = "Background"
          )
        }
        
        cat("  Successfully collected", selected_year, "data for", NamesEnglish[i], "\n")
      } else {
        cat("  No", selected_year, "data available for", NamesEnglish[i], "\n")
      }
      
    }, error = function(e) {
      cat("  Error collecting data for", Names[i], ":", e$message, "\n")
    })
  }
}


###########################################
# FUNCTION TO CREATE DISTRIBUTION PLOTS
###########################################

create_distribution_plot <- function(background_data_list, processed_data_list, title, filename, ylabel) {
  cat("Creating", title, "...\n")
  
  # Convert lists to data frames
  background_df <- data.frame()
  processed_df <- data.frame()
  
  # Process background data
  if(length(background_data_list) > 0) {
    background_df <- do.call(rbind, background_data_list)
    background_df$Type <- "Background"
  }
  
  # Process processed countries data
  if(length(processed_data_list) > 0) {
    processed_df <- do.call(rbind, processed_data_list)
    processed_df$Type <- "Processed"
  }
  
  # Combine all data
  all_data <- rbind(background_df, processed_df)
  
  if(nrow(all_data) == 0) {
    cat("No", selected_year, "data available for", title, "\n")
    return(NULL)
  }
  
  # Create separate data frames for plotting
  # Background data is the non-selected countries
  background_data <- all_data[all_data$Type == "Background", ]
  processed_data <- all_data[all_data$Type == "Processed", ]
  
  # Create the plot
  plot <- ggplot() +
    # Gray lines for all countries (background)
    {if(nrow(background_data) > 0) 
      geom_line(data = background_data, 
                aes(x = Age, y = Value, group = CountryCode), 
                color = "gray70", 
                alpha = 0.6, 
                size = 0.5)
    } +
    # Colored lines for processed countries
    {if(nrow(processed_data) > 0)
      geom_line(data = processed_data, 
                aes(x = Age, y = Value, color = Country), 
                size = 1)
    } +
    scale_color_viridis_d() +
    labs(
      title = title,
      x = "Age",
      y = ylabel
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    ) +
    guides(color = guide_legend(override.aes = list(size = 1.5)))
  
  # Display the plot
  print(plot)
  
  # Save the plot
  #ggsave(filename, plot = plot, device = "pdf", width = 12, height = 8)
  
  #cat(title, "saved as:", filename, "\n")
  
  return(plot)
}

# Main loop to process all selected countries
cat("Starting analysis for", length(countries_to_process), "countries...\n\n")

plots_list <- list()
successful_countries <- c()

for (i in seq_along(countries_to_process)) {
  t <- countries_to_process[i]
  
  result <- process_country(t)
  
  if (!is.null(result)) {
    plots_list[[Names[t]]] <- result
    successful_countries <- c(successful_countries, t)
  }
}

cat("Analysis completed successfully for", length(successful_countries), "countries:\n")
cat(paste(NamesEnglish[successful_countries], collapse = ", "), "\n")

# COLLECT ALL SELECTED YEAR DATA FOR ALL COUNTRIES
collect_all_year_data()

# CREATE ALL THREE DISTRIBUTION PLOTS
cat("Creating Distribution Plots for", selected_year, "...\n")

# Plot 1: dx distribution
pd_plot <- create_distribution_plot(
  pd_year_data, 
  pd_year_processed, 
  paste("Period Distribution of Deaths (dx) by Age in", selected_year),
  paste0("dx_", selected_year, "_Distribution.pdf"),
  "dx"
)

# Plot 2: fx distribution  
pdpe_plot <- create_distribution_plot(
  pdpe_year_data, 
  pdpe_year_processed, 
  paste("Distribution of Years Lost (fx) by Age in", selected_year),
  paste0("fx_", selected_year, "_Distribution.pdf"),
  "fx"
)

# Plot 3: Cald distribution
cald_plot <- create_distribution_plot(
  cald_year_data, 
  cald_year_processed, 
  paste("Cross-sectional Average Length of Life (CAL) Distribution by Age in", selected_year),
  paste0("CAL_", selected_year, "_Distribution.pdf"),
  "CAL dx"
)

###########################################
#  TIMING MODE ANALYSIS FOR ALL COUNTRIES
###########################################

# Define regions for timing analysis
eastern_europe <- c("BLR", "BGR", "CZE", "HUN", "POL", "RUS", 
                    "SVK", "UKR", "EST","LVA", "LTU","HRV","SVN")

northern_europe <- c("DNK", "FIN", "ISL", "IRL", 
                     "NOR", "SWE", "GBR_NIR","GBR_NP","GBR_SCO","GBRTENW", "NLD")

southern_central_europe <- c("AUT", "FRA", "DEUTNP", 
                             "GRC", "ITA", "LUX", "PRT",  
                             "ESP", "CHE", "FRATNP", "DEUTE", "DEUTW")

# Color and shape mappings for regions
color_mapping <- c(
  "Eastern Europe" = "#D55E00",  # Muted Orange-Red
  "Northern Europe" = "#0072B2",  # Muted Blue (colorblind-friendly)
  "Southern & Central Europe" = "#009E73",  # Muted Green (colorblind-friendly)
  "Non-European" = "#333333"  # Dark Gray-Black
)

shape_mapping <- c(
  "Eastern Europe" = 16,  # Solid Circle
  "Northern Europe" = 17,  # Solid Triangle
  "Southern & Central Europe" = 15,  # Solid Square
  "Non-European" = 18  # Diamond shape for non-European
)

# Initialize vectors for timing analysis
YY <- c()      # Years for regular mode
YL <- c()      # Modal age levels for regular mode
YY_fx <- c()   # Years for fx
YL_fx <- c()   # Modal age levels for fx

# Process all countries for timing analysis
for(i in 1:length(Names)){
  cat("Processing timing analysis for:", Names[i], "\n")
  
  tryCatch({
    # Get period life table data
    E <- readHMDweb(Names[i], "bltper_1x1", username.hmd, password.hmd)
    E <- E[E$Year > 1925,]
    
    # Handle potential data format issues
    if(i == 3){  # Bulgaria
      E$ex <- as.numeric(as.character(E$ex))
      E$dx <- as.numeric(as.character(E$dx))
    }
    
    # Calculate timing for regular mode
    TY <- THEyear(E)
    YYY <- ifelse(TY > min(E$Year), TY, 9999)
    YY <- c(YY, YYY)
    YL <- c(YL, THElevel(E))
    
    # Calculate timing for fx mode (dx*ex)
    TY_fx <- THEyear_fx(E)
    YYY_fx <- ifelse(TY_fx > min(E$Year), TY_fx, 9999)
    YY_fx <- c(YY_fx, YYY_fx)
    YL_fx <- c(YL_fx, THElevel_fx(E))
    
  }, error = function(e) {
    cat("  Error processing", Names[i], ":", e$message, "\n")
    YY <<- c(YY, 9999)
    YL <<- c(YL, NA)
    YY_fx <<- c(YY_fx, 9999)
    YL_fx <<- c(YL_fx, NA)
  })
}

# Create data frames for timing plots
df <- data.frame(Year = YY,
                 ModalAge = YL,
                 Country = Names)

df_fx <- data.frame(Year = YY_fx,
                    ModalAge = YL_fx,
                    Country = Names)

# Add region information to both data frames
df <- df %>% 
  mutate(Region = case_when(Country %in% eastern_europe ~ 
                              "Eastern Europe",
                            Country %in% northern_europe ~
                              "Northern Europe",
                            Country %in% southern_central_europe ~
                              "Southern & Central Europe",
                            TRUE ~ "Non-European"))

df_fx <- df_fx %>% 
  mutate(Region = case_when(
    Country %in% eastern_europe ~ "Eastern Europe",
    Country %in% northern_europe ~ "Northern Europe",
    Country %in% southern_central_europe ~ "Southern & Central Europe",
    TRUE ~ "Non-European"))

df$Region <- as.factor(df$Region)
df_fx$Region <- as.factor(df_fx$Region)

# Create Figure 5 - Period Mode Timing
F5 <- df %>% 
  filter(Year < 9999) %>% 
  ggplot(aes(x = Year, y = ModalAge, color = Region, shape = Region)) +
  geom_point(size = 3.5) +  
  geom_text_repel(aes(label = Country), size = 4, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) + 
  labs(x = "Year", y = "Modal Age", title = "Timing of the change in dominance to old age modal age in the period life table distribution of deaths") +
  xlim(1925, 1975) +
  ylim(75, 82) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.key.size = unit(1.5, "cm"),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 16, hjust = 0.5))

print(F5)
#ggsave("FigureC1.pdf", plot = F5, device = "pdf", width = 9, height = 8)

# Create Figure 6 (dx*ex)
F6 <- df_fx %>% 
  filter(Year != Inf,
         is.na(ModalAge) == FALSE,
         Year < 9999) %>% 
  ggplot(aes(x = Year, y = ModalAge, color = Region, shape = Region)) +
  geom_point(size = 3.5) +  
  geom_text_repel(aes(label = Country), size = 4, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = color_mapping) +
  scale_shape_manual(values = shape_mapping) + 
  labs(x = "Year", y = "Modal Age", title = "Figure C2. Timing of the change in dominance to old age modal age in the number of years lost distribution") +
  xlim(1986, 2023) +
  ylim(60, 85) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.key.size = unit(1.5, "cm"),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 16, hjust = 0.5))

print(F6)
#ggsave("FigureC2.pdf", plot = F6, device = "pdf", width = 9, height = 8)

######### Timing table (Table 3) ######### 
timing_table <- data.frame(
  Country = NamesEnglish,
  Year_Transition_fx = ifelse(YY_fx == Inf, "-", as.character(YY_fx)),
  Mode_Transition_fx = ifelse(is.na(YL_fx), "-", round(YL_fx,1)),
  Year_Transition_dx = ifelse(YY == 9999, "-", as.character(YY)),
  Mode_Transition_dx = ifelse(YY == 9999, "-", round(YL,1)))

View(timing_table)

cat("Analysis completed! :)")
#################################################
