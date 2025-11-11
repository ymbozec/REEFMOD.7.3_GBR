% REEFMOD-GBR connectivity settings
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 09/2019
% Edited by Suki Leung 06/2025
%
% YM (10/2025)
% Multiple dispersal models available for coral/CoTS larval connectivity 
% (1) GBR4 for both coral (Acropora) and coTS for 7 years (2008-09, 2010-11 to 2012-13, 2014-15 to 2016-17), derived from eReefs
% at 4 km resolution (Hock et al. 2017, 2019) (used until v.7.0, including CF 2024)
% (2) GBR1 for coral (Acropora) with 3 spawning events over 7 years (2015-16 to 2021-22) added in Oct 2024 (used for CF 2025)
% GBR1 is also available for CoTS but not considered here because would require new calibration of larval supply/recruitment
% (3) R05 (1 km resolution) for CoTS (early version of the one presented in Choukroun et al 2025?) -> disregarded due to
% unclear methodology, survival model, release dates.
% (4) GBRLUP (1 km resolution) for CoTS for 5 years (2018 to 2022) (Choukroun et al 2025)
% Note: only GBR4 connectivity is synchronised with WQ forcing (coral reproduction and CoTS larval survival)

% YM (11/2025)
% Modified in such a way that the user can select the dispersal model of their choice for corals and CoTS separately.
% Some connectivity data (GBR1, GBRLUP) were available for multiple spawning events within a spawning season 
% (ie, connectivity year) so were averaged per spawning season to get a standardised data structure. Connectivity
% matrices for any given year yr are then stored in two separate multi-structure arrays:
% - CONNECT_CORAL(yr).LINK_STRENGTH
% - CONNECT_COTS(yr).LINK_STRENGTH
% Information on the chosen dispersal model and sequence of connectivity years are stored in META.connectivity.
% Note changing connectivity might require modifying BH_alpha, BH_beta and minimal larval retention (min_selfseed).
% The CF2025 simulations (v.7.2) have run GBR1 for corals and GBR4 for CoTS without any change in those parameters.
% __________________________________________________________________________

% Before selecting the connectivity model, create the diagonal for larval retention
self_CORAL = sparse(diag(META.coral_min_selfseed * ones(1,META.nb_reefs)));
self_COTS = sparse(diag(META.COTS_min_selfseed * ones(1,META.nb_reefs,1)));

%% GBR4 - Coral (Acropora) and COTS connectivity (Hock et al. 2017, 2019)
load('GBR_CONNECT_7years.mat'); % Connectivity matrices for Acropora and CoTS (7 spawning seasons for each)
% We've got connectivity matrices for summers 2008-2009, 2010-11, 2011-12, 2012-13, 2014-15, 2015-16, 2016-2017
years = {'2008-09',' 2010-11', '2011-12', '2012-13', '2014-15', '2015-16', '2016-17'}; % available years

% We've got 2 missing years in the time series - fill the gap by selecting at random within the available water years
% Selection will be unique to a given run
id1 = randi(7); % randomly select one of the 7 matrices for 2013-14
id2 = randi(7); % randomly select one of the 7 matrices for 2017-18
% 
% % CORAL ------------------
% % Start with 2010-11 (2008-09 is a weirdo but can be selected at random)
% CONNECT_CORAL(1).LINK_STRENGTH = self_CORAL + GBR_CONNECT(2).ACROPORA(META.reef_ID,META.reef_ID) ; % 2010-11
% CONNECT_CORAL(2).LINK_STRENGTH = self_CORAL + GBR_CONNECT(3).ACROPORA(META.reef_ID,META.reef_ID) ; % 2011-12
% CONNECT_CORAL(3).LINK_STRENGTH = self_CORAL + GBR_CONNECT(4).ACROPORA(META.reef_ID,META.reef_ID) ; % 2012-13
% CONNECT_CORAL(4).LINK_STRENGTH = self_CORAL + GBR_CONNECT(id1).ACROPORA(META.reef_ID,META.reef_ID) ; % 2013-14 is missing
% CONNECT_CORAL(5).LINK_STRENGTH = self_CORAL + GBR_CONNECT(5).ACROPORA(META.reef_ID,META.reef_ID) ; % 2014-15
% CONNECT_CORAL(6).LINK_STRENGTH = self_CORAL + GBR_CONNECT(6).ACROPORA(META.reef_ID,META.reef_ID) ; % 2015-16
% CONNECT_CORAL(7).LINK_STRENGTH = self_CORAL + GBR_CONNECT(7).ACROPORA(META.reef_ID,META.reef_ID) ; % 2016-17
% CONNECT_CORAL(8).LINK_STRENGTH = self_CORAL + GBR_CONNECT(id2).ACROPORA(META.reef_ID,META.reef_ID) ; % 2017-18 is missing
% 
% % If not simulating all 3806 reefs but only a subset, we need to determine in-degree from source reefs outside of the
% % selected region to infer larval supply from the outside. This returns a vector of the summed in-degree probabilitites
% % of 'external' larval supply (ie, originating outside of the simulated region) for each reef included in the simulation.
% % We don't do this for CoTS because their distribution is too patchy for reliable extrapolation (might revisit this)
% if isempty(META.outside_reef_ID)==0
% 
%     CONNECT_CORAL(1).LINK_STRENGTH_ext = sum(GBR_CONNECT(2).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2010-11
%     CONNECT_CORAL(2).LINK_STRENGTH_ext = sum(GBR_CONNECT(3).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2011-12
%     CONNECT_CORAL(3).LINK_STRENGTH_ext = sum(GBR_CONNECT(4).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2012-13
%     CONNECT_CORAL(4).LINK_STRENGTH_ext = sum(GBR_CONNECT(id1).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2013-14 is missing
%     CONNECT_CORAL(5).LINK_STRENGTH_ext = sum(GBR_CONNECT(5).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2014-15
%     CONNECT_CORAL(6).LINK_STRENGTH_ext = sum(GBR_CONNECT(6).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2015-16
%     CONNECT_CORAL(7).LINK_STRENGTH_ext = sum(GBR_CONNECT(7).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2016-17
%     CONNECT_CORAL(8).LINK_STRENGTH_ext = sum(GBR_CONNECT(id2).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2017-18 is missing
% end
% 
% % Capture the selected connectivity model for coral
% META.connectivity.CORAL_dispersal_model = 'GBR4';
% % Capture the sequence of connectivity years at each time step
% my_sequence = repelem([2, 3, 4, id1, 5, 6, 7, id2], 2); % will repeat each element to get same value for two consectutive time steps (=1 year)
% my_sequence_extended = [7, 7, id2, id2, repmat(my_sequence, 1, 20)]; % 2010-11 = 1st year available so cycle back for 2008-09 and 2009-10
% META.connectivity.CORAL_sequence = [years(my_sequence_extended(1:META.nb_time_steps)); num2cell(my_sequence_extended(1:META.nb_time_steps))];
% 
% COTS ------------------
if META.doing_COTS == 1 && META.doing_COTS_connectivity == 1

    CONNECT_COTS(1).LINK_STRENGTH = self_COTS + GBR_CONNECT(2).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT_COTS(2).LINK_STRENGTH = self_COTS + GBR_CONNECT(3).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT_COTS(3).LINK_STRENGTH = self_COTS + GBR_CONNECT(4).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT_COTS(4).LINK_STRENGTH = self_COTS + GBR_CONNECT(id1).COTS(META.reef_ID,META.reef_ID) ; % 2013-14 is missing
    CONNECT_COTS(5).LINK_STRENGTH = self_COTS + GBR_CONNECT(5).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT_COTS(6).LINK_STRENGTH = self_COTS + GBR_CONNECT(6).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT_COTS(7).LINK_STRENGTH = self_COTS + GBR_CONNECT(7).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT_COTS(8).LINK_STRENGTH = self_COTS + GBR_CONNECT(id2).COTS(META.reef_ID,META.reef_ID) ; % 2017-18 is missing

end

% Capture the selected connectivity model for COTS
META.connectivity.COTS_dispersal_model = 'GBR4';
% Capture the sequence of connectivity years at each time step
my_sequence = repelem([2, 3, 4, id1, 5, 6, 7, id2], 2); % will repeat each element to get same value for two consectutive time steps (=1 year)
my_sequence_extended = [7, 7, id2, id2, repmat(my_sequence, 1, 20)]; % 2010-11 = 1st year available so cycle back for 2008-09 and 2009-10
META.connectivity.COTS_sequence = [years(my_sequence_extended(1:META.nb_time_steps)); num2cell(my_sequence_extended(1:META.nb_time_steps))];


%% GBR1: using larval connectivity simulated by eReefs GBR1 (1 km resolution) particle release and dispersal
load('GBR1_CONNECT_NEW.mat'); % Connectivity matrices for coral averaged over 7 spawning seasons (2015-16 to 2021-22).
% CORAL_release_date and COTS_release_date give the date of particle release for each (yr,event)
% Note: as of 30/10/2024, spawning is disconnected from WQ forcing on coral reproduction
years = {'2015-16',' 2016-17', '2017-18', '2018-19', '2019-20', '2020-21', '2021-22'};

% CORAL ------------------
for yr = 1:7
    CONNECT_CORAL(yr).LINK_STRENGTH = self_CORAL + GBR_CONNECT(yr).CORAL(META.reef_ID,META.reef_ID) ;
    if isempty(META.outside_reef_ID)==0
        CONNECT_CORAL(yr).LINK_STRENGTH_ext = sum(GBR_CONNECT(yr).CORAL(META.outside_reef_ID,META.reef_ID),1) ;
    end
end

% Capture the selected connectivity model for coral
META.connectivity.CORAL_dispersal_model = 'GBR1';
% Capture the sequence of connectivity years at each time step
my_sequence = repelem([1:7], 2); % will repeat each element to get same value for two consectutive time steps (=1 year)
my_sequence_extended = [my_sequence, repmat(my_sequence, 1, 20)]; % 2015-16 = 1st year available so cycle back for 7 years (2008-09 to 2013-14)
META.connectivity.CORAL_sequence = [years(my_sequence_extended(1:META.nb_time_steps)); num2cell(my_sequence_extended(1:META.nb_time_steps))];

% COTS ------------------
% Not included in 'GBR1_CONNECT_NEW.mat' (but available on request) as using GBR1 for CoTS would require a new parameterisation of
% BH_alpha and BH_beta as done for GBRLUP by suki

% if META.doing_COTS == 1 && META.doing_COTS_connectivity == 1
% 
%     for yr = 1:7
%         CONNECT_COTS(yr).LINK_STRENGTH = self_COTS + GBR_CONNECT(yr).COTS(META.reef_ID,META.reef_ID) ;
%     end
% end
% 
% % Capture the selected connectivity model for COTS
% META.connectivity.COTS_dispersal_model = 'GBR1';
% % Capture the sequence of connectivity years at each time step
% my_sequence = repelem([1:7], 2); % will repeat each element to get same value for two consectutive time steps (=1 year)
% my_sequence_extended = [my_sequence, repmat(my_sequence, 1, 20)]; % 2015-16 = 1st year available so cycle back for 7 years (2008-09 to 2013-14)
% META.connectivity.COTS_sequence = [years(my_sequence_extended(1:META.nb_time_steps)); num2cell(my_sequence_extended(1:META.nb_time_steps))];


%% GBRLUP CoTS connectivity matrices (1km resolution) - from Suki (JUN 2025)
% % Only available for CoTS - requires changing BH_alpha, BH_beta and minimal larval retention in settings_COTS and f_multiple_reef
% load('GBRLUP_CONNECT_NEW.mat');
% years = {'2018-19', '2019-20', '2020-21', '2021-22', '2022-23'};
% 
% if META.doing_COTS == 1 && META.doing_COTS_connectivity == 1
% 
%     self_COTS = 0; % decided not to add self retention values to R05 matrices
%     CONNECT_COTS(1).LINK_STRENGTH = self_COTS + GBR_CONNECT(1).COTS(META.reef_ID,META.reef_ID) ;
%     CONNECT_COTS(2).LINK_STRENGTH = self_COTS + GBR_CONNECT(2).COTS(META.reef_ID,META.reef_ID) ;
%     CONNECT_COTS(3).LINK_STRENGTH = self_COTS + GBR_CONNECT(3).COTS(META.reef_ID,META.reef_ID) ;
%     CONNECT_COTS(4).LINK_STRENGTH = self_COTS + GBR_CONNECT(4).COTS(META.reef_ID,META.reef_ID) ; 
%     CONNECT_COTS(5).LINK_STRENGTH = self_COTS + GBR_CONNECT(5).COTS(META.reef_ID,META.reef_ID) ;
% end
% 
% % Capture the selected connectivity model for COTS
% META.connectivity.COTS_dispersal_model = 'GBRLUP';
% % Capture the sequence of connectivity years at each time step
% my_sequence = repelem([1:5 1:5], 2); % will repeat each element to get same value for two consectutive time steps (=1 year)
% my_sequence_extended = [my_sequence, repmat(my_sequence, 1, 20)]; % 2018-19 = 1st year available so cycle back for 10 years (2008-09 to 2017-18)
% META.connectivity.COTS_sequence = [years(my_sequence_extended(1:META.nb_time_steps)); num2cell(my_sequence_extended(1:META.nb_time_steps))];
