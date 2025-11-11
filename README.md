# REEFMOD.7.3_GBR

A coral individual-based model that simulates coral metapopulations across 2,300 km of Australia's Great Barrier Reef (GBR). ReefMod-GBR simulates coral trajectories for >3,800 individual reefs over the recent past (hindcast 2008-2024) and projects coral trajectories into the future (forecast 2025-2100) in response to temporally- and spatially-explicit forcings of water quality, cyclones, heat stress (mass coral bleaching) and the simulated population dynamics of the coral-eating crown-of-thorns starfish (CoTS). Heat stress projections are available for 5 emission scenarios of the Shared Socioeconomic Pathway (SSP) framework, derived from an ensemble of 10 coupled atmosphere-ocean general circulation models (AOGCMS) of the CMIP6 database. Following [version 7.0](https://github.com/ymbozec/REEFMOD.7.0_GBR), the model simulates mechanisms of coral adaptation to heat stress (adaptation through selection).

_This version is still under development (and testing). Please report any bug or issue._

Compared to v7.2, this version had the following implementations (see details in [MAIN_REEFMOD_GBR](https://github.com/ymbozec/REEFMOD.7.3_GBR/blob/main/MAIN_REEFMOD_GBR.m)):
- (Mar 2025): size-specific natural mortality (whole-colony) for adult corals;
- (Apr 2025): extensive revision of the CoTS control module, including the generation of reef priority lists, the conversion of CoTS density into CoTS per manta tow, the calculation of culled CoTS when culling down to the ecological threshold (ET), several corrections around ET, the tracking of visited vs culled reefs;
- (Oct 2025): revision of coral and CoTS connectivity for more flexible selection of the dispersal model (GBR4, GBR1, GBRLUP).
- (Nov 2025): automated formatting of model outputs with option to save outputs either seasonally (following native time step of 6 months) or annually (aggregated over 2 time steps).


## Extra requirements

Running the code requires:
- inclusion into the folder /data/Climatology of the CMIP6 projections of heat stress. These projections are in the folder /Climatology/Future/CMIP6 of the repository [REEFMOD.7.0_GBR](https://github.com/ymbozec/REEFMOD.7.0_GBR)


## Instructions

The code is written in MATLAB (2023b or earlier versions).
To execute the model:
1. Download all the necessary scripts and folders.
2. Add them to your current MATLAB path.
3. In the Command Window, type:
    > run('MAIN_REEFMOD_GBR.m')

This will start the simulation. The current settings run a projection of one climate change scenario (ie, one CMIP-6 climate model under a specific scenario of carbon emission SSP - as specified by the user) for the period 2008-2100, with CoTS control as the only management intervention (counterfactual simulation). See REEFMOD.6.8_GBR/MAIN_REEFMOD_GBR.m for using CMIP-5 models as input.

The number of repeat simulations can be set with 'NB_SIMULATIONS' (currently set to 20). Simulations are then executed sequentially, each identified by the iterator "simul" (eg, from 1 to 20), which sets set a specific seed for the MATLAB random number generator, ensuring reproducibility of the results. Each simulation is stochastic, incorporating several randomised components, including the timing of future heat stress within each decade, the selection of a specific scenario of future cyclones, the initialisation of coral cover and Crown-of-Thorns Starfish (CoTS) density, the magnitude of coral mortality events, the forcing scheme of water quality. Because the runtime of one complete simulation (ie, from year 2008 to year 2100) is about 2 hours, the use of HPC resources is recommended. Shorter simulations can be obtained by setting a lower number of 6-month time steps ('NB_TIME_STEPS').

## Contact
Yves-Marie Bozec, The University of Queensland (y.bozec@uq.edu.au)


## Citation
Bozec, Y-M, AA Adam, B Arellano-Nava, AK Cresswell, V Haller-Bull, T Iwanaga, L Lachs, SA Matthews, JK McWhorter, KRN Anthony, SA Condie, PR Halloran, JC Ortiz, C Riginos, and PJ Mumby (2025) A rapidly closing window for coral persistence under global warming. Nature Communications (16), 9704. https://www.nature.com/articles/s41467-025-65015-4

## Earlier model versions for the GBR
(* modified version available from first author)

_v5.6_: Bozec Y-M, PJ Mumby (2019) ReefMod-GBR. Model description and simulation results. Appendix of the technical report provided to the Australian Government by the Reef Restoration and Adaptation Program, 39 pp. https://static1.squarespace.com/static/5990d4bec534a54aef8ec205/t/61d665f9b537d40b9a99ec4e/1641440843276/Bozec_Mumby_2019_RRAP_Appendix_T6+Modelling+Methods+and+Findings_26April_FINAL3.pdf

_v6.3-6.4_: Bozec, Y-M, K Hock, RA Mason, ME Baird, C Castro-Sanguino, SA Condie, M Puotinen, A Thompson, and PJ Mumby (2022) Cumulative impacts across Australia’s Great Barrier Reef: A mechanistic evaluation. Ecological Monographs 92(1), e01494. https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecm.1494

_v6.4*_: Mason RA, Bozec Y-M, Mumby PJ (2023) Setting sustainable limits on anchoring to improve the resilience of coral reefs. Marine Pollution Bulletin 189, 114721. https://www.sciencedirect.com/science/article/pii/S0025326X23001522

_v6.4_: Mason RA, Bozec Y-M, Mumby PJ (2023) Demographic resilience may sustain significant coral populations in a 2° C-warmer world. Global Change Biology 29:4152–4160. https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16741

_v6.4_: Mason RA, Bozec Y-M, Mumby PJ (2025) Coral bleaching and mortality overestimated in projections based on Degree Heating Months. Nature Geoscience 18(2), 120-123. https://www.nature.com/articles/s41561-024-01635-7

_v6.5*_: Castro-Sanguino C, Y-M Bozec, SA Condie, CS Fletcher, K Hock, C Roelfsema, DA Westcott, PJ Mumby (2023) Control efforts of crown‐of‐thorns starfish outbreaks to limit future coral decline across the Great Barrier Reef. Ecosphere 14:e4580. https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.4580

_v6.8*_: Skinner, C, Y-M Bozec, SA Matthews, DH Williamson, R Beeden, and PJ Mumby (2024) Advancing projections of crown-of-thorns starfish to support management interventions. Science of the Total Environment 950:175282. https://www.sciencedirect.com/science/article/pii/S0048969724054329

_v6.8*_: Skinner C, Bozec Y-M, Fletcher CS, Mumby PJ (2025) Maximising the benefits of local management for coral reefs amidst near-term environmental change. Journal of Environmental Management 392, 126627. https://www.sciencedirect.com/science/article/pii/S0301479725026039

_v7.0_: Bozec, Y-M, AA Adam, B Arellano-Nava, AK Cresswell, V Haller-Bull, T Iwanaga, L Lachs, SA Matthews, JK McWhorter, KRN Anthony, SA Condie, PR Halloran, JC Ortiz, C Riginos, and PJ Mumby (2025) A rapidly closing window for coral persistence under global warming. Nature Communications (16), 9704. https://www.nature.com/articles/s41467-025-65015-4
