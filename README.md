# Fish bioturbation extinction modelling

This repository contains code and data necessary to replicate data analysis and figures in the manuscript *"Commercial extraction of large, valuable fish disproportionately reduces sediment bioturbation"* by Mara Fischer, Thomas Williams, Martin Solan, Benjamin Harris, Tara Williams, Ceri Lewis, Julie Hawkins, and Callum Roberts.

*Abstract:* Marine protected areas are widely used to conserve biodiversity, but whether they can maintain ecosystem functioning as species are lost remains unclear. Here, we investigate how commercial extraction drives fish species loss and, in turn, affects sediment bioturbation, a key mediator of benthic biogeochemical cycling. Using trait and abundance data from marine benthic fish communities sampled with baited remote underwater video systems (BRUVS), we parameterise trait-based extinction simulations and assess how different patterns of species loss affect fish community-level bioturbation. Although bioturbation generally declined with species loss across both random and non-random extinction scenarios, the magnitude and trajectory of functional loss depended heavily on extinction order and the capacity of surviving species to compensate. Specifically, the loss of large-bodied and commercially valuable species produced the greatest reductions in bioturbation, highlighting their disproportionate functional importance. Seabed protection moderated the functional consequences of size- and value-driven species loss when functional compensation was permitted but showed weaker and more variable effects in its absence. These findings identify a previously overlooked pathway through which overfishing alters benthic ecosystem functioning in soft-sediment habitats, and highlight the value of incorporating functional traits, such as maximum body size, into fisheries management to maintain key bioturbators in fished seascapes.

## Navigation

*NZ_FishExt_Models.R* - script to prepare data for modelling, run simulation models, and create outputs, including figures in the manuscript.

### Data files:

*NZ24_fish_data.csv* - Raw data with fish relative abundance (MaxN) and biomass, collected in November and December 2024 in the Hauraki Gulf, New Zealand using baited remote underwater video systems (BRUVS):

Key columns:
- OpCode: Unique identifier of deployment
- MaxN: Relative abundance (the maximum number of individuals of each species present in a single video frame)
- binomial: scientific species name
- Site: Survey site
- Deployment: replicate deployment number per site
- CPA.status: protection status (sites are either inside a cable protection area and therefore protected or outside of it and therefore unprotected)
- bioturbation.mode: Sediment reworking mode assigned after Fischer et al. (2025)'s classification system
- bioturbation.mode.score: categorical score assigned to bioturbation mode, reflecting the increasing influence on sediment turnover
- bioturbation.frequency: frequency of sediment reworking
- bioturbation.frequency.score: categorical score assigned to bioturbation frequency, reflecting increasing frequency of interaction with the sediment
- biomass.g: biomass of species at MaxN, calculated using allometric length-weight relationships

*SAU_catch_value_data.csv* - Raw landed commercial value data from the Sea Around Us Catch database

- [Outputs](https://github.com/mfischer-edi/Irish-Sea-bioturbation/tree/main/outputs): Contains outputs (model summaries and figures) produced by running the script.

## Requirements

### Software

R version 4.4.2 or higher

### Packages

`tidyverse, patchwork, mgcv, visreg, visibly, emmeans, ggh4x, gt, broom, tidytext`

*Contact: Mara Fischer, m.fischer2@exeter.ac.uk*

----
This work is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](https://creativecommons.org/licenses/by-nc-nd/4.0/).

