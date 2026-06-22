% -----------------------------------------------------------------------------------
% Yves-Marie Bozec, UQ-MSEL (y.bozec@uq.edu.au). Created 05/2026
%
% Parameteric emulator of Degree Heating Weeks (DHW) under Marine Cloud Brightening.
% Predicts DHW for every reef and every year from (1) the counterfactual DHW projected 
% by a given climate model under a given SSP, (2) shelf position, (3) latitude and 
% (4) number of years since start of simulated cooling under a given MCB scenario.
%
% The emulator combines hurdle models built using the R package glmmTMB on MCB projections
% performed by Jessica Benthuysen (AIMS).
%
% The MCB projections where obtained by applying the night-time temperature anomaly (SSTa)
% predicted on every reef from eReefs GBR4 simulations (2015/16, 2016/17, 2019/20) conducted
% by Luke Harrison (SCU) over the entire GBR for an increase of 0.2 or 0.3 albedo assumed
% representative of the MCB effect. The resulting SSTa on each reef was applied from 1 December
% for a maximum of 150 days to the SST projections of 11 CMIP6 climate models
% under SSP1-2.6, SSP2-4.5, SSP3-7.0 and SSP5-8.5, and the maximum annual DHWs calculated.
% Three scenarios of MCB duration were developed, whereby SSTa is applied each year for 50 days, 
% 100 days or 150 days.
% 
% Hurdle models were performed on the DHW difference between MCB (DHW_MCB) and the baseline (DHW_CF)
% projections (11*4 DHW baselines) per shelf position/SSP/albedo/MCB duration.
% Only the scenario of 0.2 albedo and 100 days of MCB is currently available for each shelf position
% and 3 SSPs (SSP1-2.6, SSP2-4.5, SSP3-7.0). The function returns stochastic predictions for each reef
% based on the DHW in the baseline, the year since MCB started, and the latitude of the reef,
% with random noise consistent with the modelled variance.
%_____________________________________________________________________________________

function out = f_predict_hurdle_MCB(DHW_CF, Year, LAT, model)

% Year: year since MCB started in simulations (2015)
% DHW_CF: baseline DHW predicted for that year on a reef 
% LAT: latitude of the reef centroid (as in ReefMod)
% Returns 'out' struct with:
% p_occ - probability of occurrence predicted for each reef/year
% mu_pos - expected magnitude given occurence of cooling for each reef/year
% Y - realised difference in DHW relative to DHW_CF (stochastic) for each reef/year
% DHW_MCB - realised DHW value under MCB, ie, back-transformed response as DHW_CF - Y

%% 1. Extract and apply scalings
% Predictors were scaled to stabilise the model and facilitate convergence.
sc = model.scalings; % retrieve the parameters used for the scaling (different for inshore/midshelf/outershelf)
Year_sc  = (Year  - sc.year_center)  ./ sc.year_scale;
DHW_CF_sc   = (DHW_CF - sc.dhw_center) ./ sc.dhw_scale;
DHW_CF2     = DHW_CF.^2;
DHW_CF2_sc  = (DHW_CF2 - sc.dhw2_center) ./ sc.dhw2_scale;
LAT_sc   = (LAT - sc.lat_center) ./ sc.lat_scale;

% Each hurdle model was designed with same structure, on the response Y = DHW_MAX_CF - DHW_MAX_INT  (Y >= 0)
% Part 1: P(Y > 0)  (binomial, logit link)
% (Y > 0) ~ DHW_CF_sc + Year_sc + LAT_sc,
% Part 2: Y | Y > 0 (Gamma, log link)
% Y ~ DHW_CF_sc + DHW_CF2_sc + Year_sc + DHW_CF_sc:Year_sc + DHW_CF2_sc:Year_sc+ LAT_sc,

%% 2. OCCURRENCE (logit)
b_occ   = model.beta_occ.value; % coefficients
names_o = model.beta_occ.term; % predictor names

eta_occ = zeros(size(DHW_CF));

for i = 1:length(b_occ)
    term = names_o{i};
    if strcmp(term, '(Intercept)');   x = 1;
    elseif strcmp(term, 'DHW_CF_sc'); x = DHW_CF_sc;
    elseif strcmp(term, 'Year_sc');   x = Year_sc;
    elseif strcmp(term, 'LAT_sc');    x = LAT_sc;
    end
    eta_occ = eta_occ + b_occ(i) .* x;
end

% logistic transformation
p_occ = 1 ./ (1 + exp(-eta_occ));

%% 3. MAGNITUDE (log link)
b_pos   = model.beta_pos.value; % coefficients
names_p = model.beta_pos.term; % predictor names

eta_pos = zeros(size(DHW_CF));

for i = 1:length(b_pos)
    term = names_p{i};
    if strcmp(term, '(Intercept)');      x = 1;
    elseif strcmp(term, 'DHW_CF_sc');    x = DHW_CF_sc;
    elseif strcmp(term, 'DHW_CF2_sc');   x = DHW_CF2_sc;
    elseif strcmp(term, 'Year_sc');      x = Year_sc;
    elseif strcmp(term, 'LAT_sc');       x = LAT_sc;
    elseif strcmp(term, 'DHW_CF_sc:Year_sc');    x = DHW_CF_sc .* Year_sc;
    elseif strcmp(term, 'DHW_CF2_sc:Year_sc');   x = DHW_CF2_sc .* Year_sc;
    end
    eta_pos = eta_pos + b_pos(i) .* x;
end

% log link → positive mean
mu_pos = exp(eta_pos);

%% 4. Combine hurdle (with random noise)
% Simulate occurence
occ = rand(size(p_occ)) < p_occ; % Bernoulli for occurrence

% Gamma simulation
shape = 1 / (model.sigma_pos^2);
scale = mu_pos * model.sigma_pos^2;
Y_pos =  gamrnd(shape, scale);

Y = occ .* Y_pos;

%% 5. Outputs
out.p_occ  = p_occ;
out.mu_pos = mu_pos;
out.Y  = Y;
out.DHW_INT = DHW_CF - Y;
out.DHW_INT(out.DHW_INT<0) = 0; % DHW cannot be negative (physical constraints)
