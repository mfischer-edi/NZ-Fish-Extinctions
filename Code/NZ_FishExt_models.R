# Fish bioturbation extinction modelling
# April 2026
# Written by Mara Fischer, based on code written by Williams et al. 2025
# Contact: m.fischer2@exeter.ac.uk

# Libraries ----
library(tidyverse)
library(patchwork)
library(mgcv)
library(visreg)
library(ggh4x)
library(emmeans)
library(broom)
library(tidytext)
library(gt)

# Starting dataframe ----

data <- read.csv("NZ24_fish_data.csv")

# Raw BRUVS data from Hauraki Gulf
# one row per species per deployment

# Key columns:
## OpCode: Unique identifier of deployment
## MaxN: Relative abundance (the maximum number of individuals of each species present in a single video frame)
## binomial: scientific species name
## Site: Survey site
## Deployment: replicate deployment number per site
## CPA.status: protection status (sites are either inside a cable protection area and therefore protected or outside of it and therefore unprotected)
## bioturbation.mode: Sediment reworking mode assigned after Fischer et al. (2025)'s classification system
## bioturbation.mode.score: categorical score assigned to bioturbation mode, reflecting the increasing influence on sediment turnover
## bioturbation.frequency: frequency of sediment reworking
## bioturbation.frequency.score: categorical score assigned to bioturbation frequency, reflecting increasing frequency of interaction with the sediment
## biomass.g: biomass of species at MaxN, calculated using allometric length-weight relationships

# tidy
data <- data %>% 
  dplyr::select(OpCode, Site, Deployment, CPA.status, binomial, common.name, bioturbation.mode,
                bioturbation.mode.score, bioturbation.frequency, bioturbation.frequency.score,
                MaxN, biomass.g) %>% 
  rename(Abundance = MaxN, Biomass = biomass.g, Prot = CPA.status,
         Mode = bioturbation.mode, Mi = bioturbation.mode.score,
         Freq = bioturbation.frequency, Fi = bioturbation.frequency.score) %>% 
  mutate(ID = match(OpCode, unique(OpCode))) %>% 
  mutate(across(c(OpCode, Site, Deployment, Prot, binomial, ID, Mode, Freq), as.factor))

# Remove non-bioturbators
data <- data %>% 
  na.omit(bioturbation.mode) %>% 
  mutate(binomial = fct_drop(binomial))

unique(data$binomial)

# Data prep ----

# aggregate across all deployments from the same site

## Calculate mean and total abundance and biomass
prot_mean <- data %>%
  group_by(Prot, binomial, Mode) %>%
  summarise(Mi = first(Mi), Fi = first(Fi), 
            Bi = mean(Biomass), 
            Ai = mean(Abundance), 
            Btot = sum(Biomass), 
            Atot = sum(Abundance),
            Bind = sum(Btot)/sum(Atot), # mean individual biomass per species per protection (to be used in BP calcs)
            .groups = "drop")

# what's species richness per site?
richness <- data %>% 
  group_by(Site, Prot) %>% 
  summarise(nsp = length(unique(binomial)))

# min richness: 3 for both
# max richness: 10 for inside, 8 for outside

## Catch value ----

# explore Sea Around Us Catch database which includes landed values
# good indicator of value as the financial gains from catch when it is landed is primary motivator for fishers to go fishing
# rather than using export value

SAU <- read.csv("SAU_catch_value_data.csv", header = T)

(unique(SAU$scientific_name)) # 217
unique(SAU$year) # 1950-2019

# use just 2019 as the most recent
sau19 <- SAU %>% 
  filter(year == 2019)

length(unique(sau19$scientific_name)) # 190
unique(sau19$scientific_name)

sau_spp <- sau19 %>% 
  filter(scientific_name %in% data$binomial)

unique(sau_spp$scientific_name) # 11

# missing species: Dipturus innominatus, Tetronarce nobiliana, Pseudocaranx georgianus,
# Caesioperca lepidoptera, Bodianus unimaculatus

# Pseudocaranx georgianus is under P. dentex, so include that
# use Rajidae value for Dipturus innominatus as they are exploited
# for all other species that are not commercially exploited, add a nominal baseline value of $0.1 USD. 

sau_spp <- sau19 %>% 
  filter(scientific_name %in% data$binomial | scientific_name == "Pseudocaranx dentex" | scientific_name == "Rajidae")

unique(sau_spp$scientific_name) # 13

# aggregate
spp_value <- sau_spp %>% 
  group_by(scientific_name, common_name) %>% 
  summarise(tonnes = mean(tonnes, na.rm = T),
            landed_value = mean(landed_value, na.rm = T)) %>% 
  mutate(tonne_value = landed_value/tonnes, # calculate unit value (per tonne)
         kg_value = tonne_value/1000) # convert to value per kg

# add to aggregated bioturbation dataframe

# first make sure al names match
spp_value <- spp_value %>% 
  mutate(binomial = case_when(scientific_name == "Pseudocaranx dentex" ~ "Pseudocaranx georgianus",
                              scientific_name == "Rajidae" ~ "Dipturus innominatus", 
                              .default = scientific_name))

prot_mean <- prot_mean %>% 
  left_join(spp_value, by = c("binomial")) %>% 
  dplyr::select(- tonnes, - landed_value, - tonne_value, 
                - scientific_name, - common_name)

# Add very low value for species not in SAU database - these are not commercially exploited,
# but setting their value to 0.01
prot_mean$kg_value[is.na(prot_mean$kg_value)] <- 0.1


# Run simulations ----

## First we define the extinction scenarios.
## The code then runs through each scenario and simulates the sequential loss of species
## after species (ordered by the different vulnerabilities defined below). Each time as species
## is lost, a revised BPc is calculated, and if functional compensation is allowed, 
## a species from the same bioturbation mode is randomly chosen to compensate for lost biomass.


# --- 1. DEFINE SCENARIO MATRIX ---
# This creates the 15 valid combinations
scenarios <- expand.grid(
  Zone = c("inside", "outside"),
  Driver = c("random", "size", "rarity", "value"),
  Comp = c("none", "functional"),
  stringsAsFactors = FALSE
)

all_results <- list()
all_contributions <- list()
all_counts <- list()  

# --- 2. THE SCENARIO RUNNER ---
for (i in 1:nrow(scenarios)) {
  
  # Get settings for this specific run
  current_zone   <- scenarios$Zone[i]
  current_driver <- scenarios$Driver[i]
  current_comp   <- scenarios$Comp[i]
  
  message(paste("Running:", current_zone, "|", current_driver, "|", current_comp))
  
  # Create a temporary list to hold steps for THIS scenario only (for contributions)
  scenario_step_list <- list()
  
  # STEP A: Prepare sim_data for this zone
  # (Assumes 'prot_mean' has both 'inside' and 'outside' data)
  temp_sim_data <- prot_mean %>%
    filter(Prot == current_zone, Ai > 0) 
  
  # STEP B: Update Extinction Probabilities based on Driver
  if (current_driver == "size") {
    # Body size (largest first)
    temp_sim_data$ExtinctProb <- temp_sim_data$Bind / sum(temp_sim_data$Bind)
  } else if (current_driver == "rarity") {
    # Rarity: Inverse of abundance (1/Ai)
    temp_sim_data$ExtinctProb <- (1/temp_sim_data$Ai) / sum(1/temp_sim_data$Ai)
  } else if (current_driver == "value") {
    # Commercial value (highest value first)
    temp_sim_data$ExtinctProb <- temp_sim_data$kg_value / sum(temp_sim_data$kg_value)
  } else {
    # Random (Flat probability)
    temp_sim_data$ExtinctProb <- 1/nrow(temp_sim_data)
  }
  
  # --- 3. THE MASTER SIMULATION LOOP ---
  nsims <- 1000
  nsp <- nrow(temp_sim_data)
  
  output <- expand.grid(Simulation = 1:nsims, Nsp = nsp:1, Value = NA)
  
  # Create a blank slate to collect sums
  rolling_sum <- temp_sim_data %>% 
    dplyr::select(Species = binomial) %>%
    expand_grid(Nsp = nsp:1) %>%
    mutate(SumBPi = 0)
  
  # Create a slate to collect the count of active species
  active_count_sum <- data.frame(Nsp = nsp:1, TotalActiveCount = 0)
  
  for (sim_count in 1:nsims) {
    sim_data <- temp_sim_data
    sim_data$AiSim <- sim_data$Ai
    sim_data$EPSim <- sim_data$ExtinctProb
    
    for (n in nsp:1) {
      # 1. Calculate Function
      current_bpc <- sum(sim_data$AiSim * (sqrt(sim_data$Bind) * sim_data$Mi * sim_data$Fi), na.rm = T)
      output[output$Simulation == sim_count & output$Nsp == n, "Value"] <- current_bpc
      
      # 2. Pick species to go extinct
      if (sum(sim_data$EPSim, na.rm = T) == 0) break
      extinct_idx <- which(cumsum(sim_data$EPSim) >= runif(1))[1]
      lost_biomass <- sim_data$AiSim[extinct_idx] * sim_data$Bind[extinct_idx]
      extinct_mode  <- sim_data$Mode[extinct_idx]
      extinct_name  <- sim_data$binomial[extinct_idx]
      
      # 3. Compensation
      if (current_comp == "functional") {
        # Pick one random survivor in the same bioturbation mode
        candidates <- which(sim_data$Mode == extinct_mode & sim_data$AiSim > 0 & sim_data$binomial != extinct_name)
        
        if (length(candidates) > 0) {
          # Randomly select one index from candidates
          comp_idx <- if(length(candidates) > 1) sample(candidates, 1) else candidates
          # Compensate biomass (Round to nearest whole fish)
          sim_data$AiSim[comp_idx] <- sim_data$AiSim[comp_idx] + round(lost_biomass / sim_data$Bind[comp_idx])
        }
        sim_data$AiSim[extinct_idx] <- 0
        
      } else { 
        # NONE
        sim_data$AiSim[extinct_idx] <- 0
      }
      
      # Record current contribution
      current_BPi <- sim_data$AiSim * (sqrt(sim_data$Bind) * sim_data$Mi * sim_data$Fi)
      # Update the specific rows in rolling_sum where Nsp == n
      rolling_sum$SumBPi[rolling_sum$Nsp == n] <- rolling_sum$SumBPi[rolling_sum$Nsp == n] + current_BPi
      
      
      # 4. Cleanup for next step
      sim_data$EPSim[extinct_idx] <- 0
      if (sum(sim_data$EPSim, na.rm = T) > 0) sim_data$EPSim <- sim_data$EPSim / sum(sim_data$EPSim, na.rm = T)
    }
  }
  
  # Summarise contributions Scenario 
  all_contributions[[i]] <- rolling_sum %>%
    mutate(MeanBPi = SumBPi / nsims,
           Driver = current_driver, 
           Comp = current_comp, 
           Zone = current_zone) %>%
    dplyr::select(-SumBPi)
  
  # Tag results and save to list
  all_results[[i]] <- output %>% 
    mutate(Zone = current_zone, Driver = current_driver, Comp = current_comp)
  
}

# Combine all scenarios into one dataframe
final_df <- bind_rows(all_results)

str(final_df)

final_df$Comp <- factor(final_df$Comp, levels = c("none", "functional"))
final_df$Driver <- factor(final_df$Driver, levels = c("random", "size", "rarity", "value"))

final_df <- final_df %>% 
  mutate(Prot = if_else(Zone == "inside", "protected", "unprotected"))

## Sort species contribution data

# Combine the list into one dataframe
heatmap_df <- bind_rows(all_contributions)

count_df <- bind_rows(all_counts)

# Calculate max BPc and tipping points per scenario ----

# --- 1. CALCULATE MEAN BPc PER STEP FOR EACH SCENARIO ---
scenario_means <- final_df %>%
  group_by(Zone, Driver, Comp, Nsp) %>%
  summarise(Mean_Bpc = mean(Value, na.rm = TRUE), .groups = "drop")

# --- 2. EXTRACT THE STARTING VALUE (Max Species Richness) ---
# This isolates the BPc before any extinctions happen
starting_values <- scenario_means %>%
  group_by(Zone, Driver, Comp) %>%
  filter(Nsp == max(Nsp)) %>%
  select(Zone, Driver, Comp, Starting_Bpc = Mean_Bpc)

# --- 3. CALCULATE % CHANGE AND EXTRACT THE PEAK ---
bpc_spikes <- scenario_means %>%
  left_join(starting_values, by = c("Zone", "Driver", "Comp")) %>%
  mutate(
    # Calculate % change relative to the starting community
    Pct_Change = ((Mean_Bpc - Starting_Bpc) / Starting_Bpc) * 100
  ) %>%
  group_by(Zone, Driver, Comp) %>%
  summarise(
    Starting_Bpc = first(Starting_Bpc),
    # Find the maximum mean BPc hit during the simulation
    Max_Bpc_Reached = max(Mean_Bpc),
    # Extract the absolute maximum % increase
    Max_Pct_Increase = max(Pct_Change),
    # Identify at exactly how many species remaining the peak occurred
    Nsp_At_Peak = Nsp[which.max(Mean_Bpc)],
    .groups = "drop"
  ) %>%
  # Filter for where an actual spike occurred (> 0% increase)
  # and sort by the largest spike
  arrange(desc(Max_Pct_Increase))

write.csv(bpc_spikes, "Outputs/BPc_compensation_perc_change.csv")


#  GAM ----

# one model with a three-way interaction between richness, protection and scenario

# Log BPc for modelling
final_df$Value_log <- log(final_df$Value)

# create the interaction between protection and scenario in the data first:

final_df <- final_df %>% 
  mutate(Scenario = interaction(Prot, Driver, Comp))

# Run GAM

gam1 <- bam(Value_log ~ s(Nsp, by = Scenario, bs = "tp", k = 5) +
              #s(Nsp, bs = "tp") + 
              Scenario,
            data = final_df, method = "fREML")
summary(gam1)
plot(gam1) 

gam.check(gam1)
visibly::plot_gam_check(gam1)

gam1_sum <- tidy(gam1)

# Compare models
# Run different combinations and do backwards stepwise selection to find best model
# using maximum restricted likelihood (ML)

GAM_1.0 <- bam(Value_log ~ s(Nsp, by = Scenario, bs = "tp", k = 5) + s(Nsp, bs = "tp", k = 5) + Scenario,
               data = final_df, method = "ML") 
GAM_1.0_summary <- summary(GAM_1.0)

GAM_1.1 <- bam(Value_log ~ s(Nsp, by = Scenario, bs = "tp", k = 5) + Scenario,
               data = final_df, method = "ML")

GAM_1.1_summary <- summary(GAM_1.1)

GAM_1.2 <- bam(Value_log ~ s(Nsp, by = Scenario, bs = "tp", k = 5) + s(Nsp, bs = "tp", k = 5),
               data = final_df, method = "ML") 

GAM_1.2_summary <- summary(GAM_1.2)

GAM_1.3 <- bam(Value_log ~ s(Nsp, bs = "tp", k = 5) + Scenario,
               data = final_df, method = "ML")

GAM_1.3_summary <- summary(GAM_1.3)

GAM_1.4 <- bam(Value_log ~ s(Nsp, by = Scenario, bs = "tp", k = 5),
               data = final_df, method = "ML")

GAM_1.4_summary <- summary(GAM_1.4)

GAM_1.5 <- bam(Value_log ~ s(Nsp, bs = "tp", k = 5),
               data = final_df, method = "ML")

GAM_1.5_summary <- summary(GAM_1.5)

GAM_1.6 <- bam(Value_log ~ Scenario,
               data = final_df, method = "ML")

GAM_1.6_summary <- summary(GAM_1.6)

GAM_1.7 <- gam(Value_log ~ 1,
               data = final_df, method = "ML")

GAM_1.7_summary <- summary(GAM_1.7)

model_summaries <- list(GAM_1.0_summary,
                        GAM_1.1_summary, 
                        GAM_1.2_summary,
                        GAM_1.3_summary, 
                        GAM_1.4_summary, 
                        GAM_1.5_summary,
                        GAM_1.6_summary, 
                        GAM_1.7_summary)

AIC_table <- AIC(GAM_1.0, GAM_1.1, GAM_1.2, GAM_1.3,
                 GAM_1.4, GAM_1.5, GAM_1.6, GAM_1.7) %>%
  mutate(
    R2 = sapply(model_summaries, function(x) x$r.sq),
    Deviance_explained = sapply(model_summaries, function(x) x$dev.expl)
  ) %>%
  arrange(AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

head(AIC_table, 8)

# GAM 1.1 is best (with interaction term plus scenario)

write.csv(AIC_table, "Outputs/GAM_model_selection.csv")

# Re-fit fnal model with f-reml
gam1 <- bam(Value_log ~ s(Nsp, by = Scenario, bs = "tp", k = 5) + Scenario,
            data = final_df, method = "fREML")

summary(gam1)

gam.check(gam1)
visibly::plot_gam_check(gam1)

gam1_sum <- tidy(gam1)

write.csv(gam1_sum, "Outputs/GAM_all_summary.csv")

# Plot GAM curves

megaDUMDAT<-data.frame(Scenario= factor(c(rep('protected.random.none', 200), 
                                          rep('protected.random.functional', 200),
                                          rep('protected.size.none', 200), 
                                          rep('protected.size.functional', 200),
                                          rep('protected.rarity.none', 200), 
                                          rep('protected.rarity.functional', 200),
                                          rep('protected.value.none', 200), 
                                          rep('protected.value.functional', 200),
                                          rep('unprotected.random.none', 200), 
                                          rep('unprotected.random.functional', 200),
                                          rep('unprotected.size.none', 200), 
                                          rep('unprotected.size.functional', 200),
                                          rep('unprotected.rarity.none', 200), 
                                          rep('unprotected.rarity.functional', 200),
                                          rep('unprotected.value.none', 200), 
                                          rep('unprotected.value.functional', 200))
),
Nsp=c(rep(seq(from=1, to=14, length=200))
))

## Predict and bind dummy data and prediction 

mega_pred <- predict(gam1, newdata = megaDUMDAT, se=T)
PREDmega<-cbind(megaDUMDAT,mega_pred)

str(PREDmega)
unique(PREDmega$Scenario)

PREDmega$Scenario <- factor(PREDmega$Scenario, levels = c("protected.random.none", "unprotected.random.none", 
                                                          "protected.random.functional", "unprotected.random.functional",
                                                          "protected.size.none", "unprotected.size.none", 
                                                          "protected.size.functional", "unprotected.size.functional",  
                                                          "protected.rarity.none", "unprotected.rarity.none",
                                                          "protected.rarity.functional", "unprotected.rarity.functional",
                                                          "protected.value.none", "unprotected.value.none",
                                                          "protected.value.functional", "unprotected.value.functional"))

PREDmega <- PREDmega %>% 
  separate(Scenario, into = c("Prot", "Driver", "Comp"), sep = "\\.", remove = F)

PREDmega$Driver <- factor(PREDmega$Driver, levels = c("random", "size", "rarity", "value"))

# add other variables

PREDmega <- PREDmega %>% 
  mutate(Comp_label = if_else(Comp == "none", "No compensation", "Functional compensation"))

PREDmega$Comp_label <- factor(PREDmega$Comp_label, levels = c("No compensation",
                                                              "Functional compensation"))

unique(PREDmega$Comp_label)

PREDmega$RowLabel <- "Extinctions ordered by:"


## Plotting 

## logged and facet_wrap

(GAM_mega_facet <- ggplot(PREDmega,aes(x = Nsp, y = (fit), group = Scenario)) +
    geom_line(colour = "blue")+
    geom_ribbon(fill = "blue", aes(ymin= (fit) + 1.96*((se.fit)),
                                   ymax= (fit) - 1.96*((se.fit))),alpha=.2)+
    theme_bw()+
    scale_y_continuous(name ="log(BPc)")+
    #scale_x_continuous(limits = c(0, 13), name = "Species richness")+
    facet_wrap(~Scenario, ncol = 4) +
    #labs(title = "Unprotected") +
    theme(legend.position="none",
          axis.text.x  = element_text(size=9),
          axis.text.y  = element_text(size=9),
          strip.background = element_rect(fill = "white"), # Clean label boxes
          #strip.text = element_text(face = "bold", size = 12), 
    )
)

# Pairwise comparisons ----

# with emmeans

# Compare the trends (curves) of the scenarios
# This calculates whether the derivatives/slopes of the curves differ pairwise

curve_comps <- emtrends(gam1, pairwise ~ Scenario, var = "Nsp")

# View the pairwise tests (look for specific pairs we are interested in, e.g., protected vs unprotected)
summary(curve_comps$contrasts)

summary(curve_comps)

pairs <- as.data.frame(curve_comps$contrasts)

slopes <- as.data.frame(curve_comps$emtrends)

# Format into clean table

# 1. Convert the emtrends objects to normal dataframes
trends_df <- as.data.frame(curve_comps$emtrends)
contrasts_df <- as.data.frame(curve_comps$contrasts)

# 2. Separate the "Scenario" string into its components 
# (e.g., "protected.size.functional" becomes Status = protected, Driver = size, Comp = functional)
trends_clean <- trends_df %>%
  separate(Scenario, into = c("Status", "Driver", "Comp"), sep = "\\.") %>%
  mutate(Combined_Scenario = paste(tools::toTitleCase(Driver), tools::toTitleCase(Comp), sep = "."))

# 3. Pivot the table so Protected and Unprotected sit on the same row and merge trend and SE with ±
table_data <- trends_clean %>%
  mutate(Trend_SE = sprintf("%.3f ± %.3f", Nsp.trend, SE)) %>%
  select(Combined_Scenario, Status, Trend_SE) %>%
  pivot_wider(
    names_from = Status, 
    values_from = Trend_SE
  )

# 4. Isolate the specific protected vs unprotected contrasts we care about
# (Filtering giant contrast matrix down to just the 8 main comparisons)
just_main_contrasts <- contrasts_df %>%
  filter(stringr::str_detect(contrast, "protected\\.(.*) - unprotected\\.\\1")) %>%
  separate(contrast, into = c("Prot", "Unprot"), sep = " - ") %>%
  separate(Prot, into = c("trash", "Driver", "Comp"), sep = "\\.") %>%
  mutate(Combined_Scenario = paste(tools::toTitleCase(Driver), tools::toTitleCase(Comp), sep = ".")) %>%
  select(Combined_Scenario, estimate, t.ratio, p.value)

# 5. Merge them together into one master data frame
table_df <- left_join(table_data, just_main_contrasts, by = "Combined_Scenario")

# Create the final table
final_table <- table_df %>%
  gt() %>%
  
  # Group the rows into "No Compensation" vs "Functional Compensation"
  tab_row_group(
    label = "Functional compensation",
    rows = stringr::str_detect(Combined_Scenario, "\\.Functional")
  ) %>%
  tab_row_group(
    label = "No compensation",
    rows = stringr::str_detect(Combined_Scenario, "\\.None")
  ) %>%
  
  # Clean up the row names by removing the ".None" or ".Functional" suffixes
  text_transform(
    locations = cells_body(columns = Combined_Scenario),
    fn = function(x) stringr::str_remove(x, "\\..*")
  ) %>%
  
  # Rename the columns
  cols_label(
    Combined_Scenario = "Extinction driver",
    protected = "Slope ± SE",
    unprotected = "Slope ± SE",
    estimate = "Contrast estimate",
    t.ratio = "t-ratio",
    p.value = "p-value"
  ) %>%
  
  # Create spanner headers to group Protected vs Unprotected
  tab_spanner(label = "Protected area", columns = c(protected)) %>%
  tab_spanner(label = "Unprotected area", columns = c(unprotected)) %>%
  
  # Format numbers to 3 decimal places
  fmt_number(columns = c(estimate, t.ratio), decimals = 3) %>%
  
  # Format the p-values to 3 decimal places
  fmt_number(columns = p.value, decimals = 3) %>%
  
  # Automatically convert 0.000 into "< 0.001"
  text_transform(
    locations = cells_body(columns = p.value),
    fn = function(x) {
      ifelse(x == "0.000" | as.numeric(x) < 0.001, "< 0.001", x)
    }
  ) %>%
  
  # Stylistic formatting: Bold the significant p-values automatically
  # (Accounts for both numeric values less than 0.05 and the new "< 0.001" string)
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p.value, 
      rows = p.value == "< 0.001" | as.numeric(p.value) < 0.05
    )
  ) %>%
  
  # Apply classic academic horizontal table lines (no vertical lines)
  opt_table_lines(extent = "default") %>%
  tab_options(
    row_group.font.weight = "bold",
    table.border.top.style = "solid",
    table.border.top.color = "black",
    table.border.bottom.style = "solid",
    table.border.bottom.color = "black",
    table_body.border.bottom.style = "solid",
    table_body.border.bottom.color = "black"
  )

# View the finished table in your RStudio viewer
final_table

gtsave(final_table, filename = "Outputs/GAM_scenario_slopes_pairs.docx",  format = "docx")

## Plot GAM SE ----

label_df2 <- PREDmega %>%
  group_by(RowLabel, Driver, Comp_label, Prot) %>%
  summarize(.groups = "drop") %>% # This gets 1 unique row per panel
  arrange(Driver, Comp_label, Prot) %>% # Standard order: across then down
  mutate(
    panel_letter = paste0(letters[1:n()]),
    Nsp = 2,   
    se.fit = 0.025 
  )

(se_plots <- ggplot(PREDmega, aes(x = Nsp, y = se.fit)) +
    geom_line(colour = "gold", 
              size = 1) +
    scale_x_continuous(expand = c(0, 0.01), breaks = c(3, 6, 9, 12, 15)) +
    
    facet_nested(rows = vars(RowLabel, Driver), 
                 cols = vars(Comp_label, Prot)) +
    geom_text(data = label_df2, 
              aes(x = Nsp, y = se.fit, label = panel_letter),
              #hjust = - 1,  # Adjust these to move text away from the edge
              vjust = -0.8, 
              size = 5,
              fontface = "bold",
              inherit.aes = FALSE) +
    theme_bw() + 
    theme(
      legend.position = "none",
      #legend.position = "right",
      strip.background = element_rect(fill = "white"), # Clean label boxes
      strip.text = element_text(face = "bold", size = 12), 
      title = element_text(face = "bold", size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Species Richness", y = "log(BPc) s.e.", 
         #title = "Functional Compensation"
    )
)

ggsave(se_plots, filename = "outputs/GAM_se_plot.png", width = 8, height = 6, dpi = 600)


# Plot extinction curves ----

# prepare extinctions output data

final_df <- final_df %>% 
  mutate(Comp_label = if_else(Comp == "none", "No compensation", "Functional compensation"),
         Prot = if_else(Zone == "inside", "protected", "unprotected"))

final_df$Comp_label <- factor(final_df$Comp_label, levels = c("No compensation",
                                                              "Functional compensation"))
final_df$Driver <- factor(final_df$Driver, levels = c("random", "size", "rarity", "value"))

final_df$RowLabel <- "Extinctions ordered by:"

# Create label_df by borrowing the structure of your existing data
label_df <- final_df %>%
  group_by(RowLabel, Driver, Comp_label, Prot) %>%
  summarize(.groups = "drop") %>% # This gets 1 unique row per panel
  arrange(Driver, Comp_label, Prot) %>% # Standard order: across then down
  mutate(
    panel_letter = paste0(letters[1:n()]),
    Nsp = 13,   
    Value = 4 
  )

(ext_plots <- ggplot(final_df, aes(x = Nsp, y = log(Value))) +
    geom_point(colour = "grey80", alpha = 0.1, size = 0.5) +
    stat_density2d(aes(fill = after_stat(level), alpha = after_stat(level)),
                   size = 3, 
                   bins = 20, 
                   geom = 'polygon', 
                   h = c(3, 0.3)
    ) +
    # Add GAM lines from mega GAM with all scenarios
    geom_line(data = PREDmega, aes(x = Nsp, y = fit), 
              inherit.aes = F, colour = "gold", 
              #linetype = "dashed", 
              size = 0.5) +
    scale_x_continuous(expand = c(0, 0.01), breaks = c(3, 6, 9, 12, 15)) +
    facet_nested(rows = vars(RowLabel, Driver), 
                 cols = vars(Comp_label, Prot)) +
    scale_fill_viridis_c(option = "magma", 
                         limits = c(0, 0.4),
                         oob = scales::squish) +
    geom_text(data = label_df, 
              aes(x = Nsp, y = Value, label = panel_letter),
              #hjust = - 1,  # Adjust these to move text away from the edge
              vjust = -0.8, 
              size = 5,
              fontface = "bold",
              inherit.aes = FALSE) +
    theme_bw() + 
    theme(
      legend.position = "none",
      #legend.position = "right",
      strip.background = element_rect(fill = "white"), # Clean label boxes
      strip.text = element_text(face = "bold", size = 12), 
      title = element_text(face = "bold", size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Species Richness", y = "log(BPc)")
)

#ggsave(ext_plots, filename = "outputs/extinctions/extinction_curves_final.png",
 #      width = 7, height = 6, dpi = 300)

# With GAM lines
ggsave(ext_plots, filename = "outputs/extinctions/extinction_curves_all_GAMs.png",
       width = 7, height = 6, dpi = 300)

# Taxonomy plots ----

# Prepare heatmap dataframe:

# 1. Calculate Relative Contribution (%) 
# Group by Scenario and Nsp so the percentages in each vertical column add to 100
heatmap_df <- heatmap_df %>%
  group_by(Zone, Driver, Comp, Nsp) %>%
  mutate(RelCont = (MeanBPi / sum(MeanBPi)) * 100) %>%
  ungroup()

# 2. Apply "Capping" logic
# This prevents the dominant species (>20%) from washing out the smaller ones
heatmap_df <- heatmap_df %>%
  mutate(HeatValue = ifelse(RelCont > 20, 20, RelCont))

# 3. Fix factors for consistent plotting order
heatmap_df$Comp <- factor(heatmap_df$Comp, levels = c("none", "functional"))
heatmap_df$Driver <- factor(heatmap_df$Driver, levels = c("random", "size", "rarity", "value"))

# 4. Clean up NAs and the >20% cap
heatmap_df <- heatmap_df %>%
  mutate(
    # If RelCont is NaN or NA (because sum was 0), make it 0
    RelCont = ifelse(is.nan(RelCont) | is.na(RelCont), 0, RelCont),
    # RE-APPLY CAP: Ensure it's numeric and explicitly capped
    HeatValue = pmin(RelCont, 20) 
  )

# 5. Calculate the rank for each species WITHIN each scenario
heatmap_df <- heatmap_df %>%
  group_by(Zone, Driver, Comp) %>%
  # Find the contribution at the start (Max Nsp) to determine the order
  mutate(StartCont = RelCont[Nsp == max(Nsp)][match(Species, Species[Nsp == max(Nsp)])]) %>%
  # Rank them (1 = highest contributor)
  mutate(Rank = dense_rank(desc(StartCont))) %>%
  ungroup()

# 6. Fix the Cap and NAs
heatmap_df <- heatmap_df %>%
  mutate(
    RelCont = replace_na(RelCont, 0),
    HeatValue = pmin(RelCont, 20)
  )

# 7. Define specific vlines 
# Make sure the names match Zone/Comp/Driver columns exactly
vline_df <- data.frame(
  Prot = rep(c("protected", "unprotected"), 2),
  Comp = rep(c("none", "functional"), each = 2),
  intercept = c(10, 8) # Example values
)

# 8. Add common names for plotting

common <- data %>% 
  select(common.name, binomial) %>% 
  distinct(common.name, binomial)
  
heatmap_df <- heatmap_df %>% 
  left_join(common, by = c("Species" = "binomial"))

head(heatmap_df)

# 9. The Plotting Function
# This creates one out of the 16 panels
make_mini_plot <- function(z, d, c, plot_tag) {
  
  # Filter data for this specific panel
  sub_data <- heatmap_df %>% 
    filter(Zone == z, Driver == d, Comp == c) %>%
    mutate(Prot = if_else(Zone == "inside", "protected", "unprotected"),
           common.name = reorder(common.name, StartCont)) # Order by starting contribution
  
  # Get the specific Prot string for This plot to find the vline
  current_prot <- if_else(z == "inside", "protected", "unprotected")
  
  # Look up the line using 'Prot'
  this_vline <- vline_df %>% 
    filter(Prot == current_prot, Comp == c) %>% 
    pull(intercept)
  
  ggplot(sub_data, aes(x = Nsp, y = common.name, fill = HeatValue)) +
    geom_tile() +
    geom_vline(aes(xintercept=3), col = "red", linetype = "dashed") +
    geom_vline(xintercept = this_vline, linetype = "dashed", color = "green") +
    scale_fill_gradient(low = "#eff3ff", 
                        high = "#08519c",
                        limits = c(0, 20), 
                        breaks = c(0, 5, 10, 15, 20),
                        labels = c("0", "5", "10", "15", ">20%"),
                        name =  "Contribution \nto BPc (%)") +
    annotate("text", x = Inf, y = -Inf, label = plot_tag, 
             hjust = 1.5, vjust = -1, 
             fontface = "bold", size = 4) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 8),
      plot.tag = element_text(face = "bold", size = 10)
    ) +
    labs(x = "Species richness")
}

# 10: Create the list of 16 plots
# This loop goes through every combination of variables
plot_list <- list()
counter <- 1

for (d in levels(heatmap_df$Driver)) {
  for (c in levels(heatmap_df$Comp)) {
    for (z in unique(heatmap_df$Zone)) {
      
      # We pass letters[counter] to give each plot a unique letter
      current_tag <- paste0(letters[counter])
      
      plot_list[[counter]] <- make_mini_plot(z, d, c, current_tag)
      
      counter <- counter + 1
    }
  }
}

# 11: Assemble with patchwork ---
# This creates the 4x4 grid (4 columns = 2 Comp x 2 Zone)
final_plot <- wrap_plots(plot_list, ncol = 4, guides = "collect", axes = "collect") 

final_plot

ggsave(final_plot, filename = "outputs/all_taxonomy_plots.png",
       width = 8, height = 6, dpi = 300)
