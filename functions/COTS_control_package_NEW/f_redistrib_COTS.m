% Initially written by K. Hock (Jun 2019) then re-worked by C. Castro-Sanguino (Aug 2021)
% Refined by YM (Sep 2025)

function [COTS_sites_density, stored_betarndN]=f_redistrib_COTS(current_COTS_forControl, nb_cntrl_sites)

% function that redistributes reef-level COTS population to individual
% control sites using beta infalted random numbers at the moment
mu = 0.3324884;
sigma = 0.2749365;
nu = 2.717073;
tau = 0.02195122;

nbreefs = length(nb_cntrl_sites);

COTS_sites_density=cell([nbreefs,1]);
stored_betarndN=cell([nbreefs,1]);

for i=1:length(nb_cntrl_sites)

    if nb_cntrl_sites(i)==1 % if only one site
        COTS_sites_density(i,1)={current_COTS_forControl(i,:)};

    else
        check_zero=0;
        while check_zero==0 % to avoid all random numbers beign zero and then not redistributing anything
            tot=1;
            Nsites = nb_cntrl_sites(i);
            betarndN=nan(1,Nsites);
            for s=1:(Nsites-1)
                r=betainfrnd(1, mu, sigma, nu, tau);%%Using random-intercept model
                betarndN(1,s)=tot*r;
                tot=tot-betarndN(1,s);
                if tot<0%if for some reason this becomes negative, break
                    betarndN(1,s:(Nsites-1))=0;
                    break;
                end
            end
            betarndN(1,Nsites) = 1-nansum(betarndN(1,:));
            
            if sum(betarndN)>0%repeat if all random numbers zero
                check_zero=1;
            end
        end
        
        % YM (03/2025): might generate negative props (tiny ones), just force to 0 here
        betarndN(betarndN<0)=0;

        %now use the random numbers to redistribute reef-level COTS population; all
        %size classes are redistributed usign the same number; would be interesting to know whether outbreaks have different size class distribution from background populations
        %this_reef_sites=zeros(META.COTS_control_sites(i,2),META.COTS_maximum_age);
        this_reef_sites=zeros(nb_cntrl_sites(i), size(current_COTS_forControl,2));
        
        for j=1:nb_cntrl_sites(i)
            
            %this_reef_cots_4redistrib=current_COTS(idx,:)*Newlist(i,2);%note that this redistributes the same size class distrubtion to all sites and in all circumstances
            this_reef_cots_4redistrib = current_COTS_forControl(i,:)*nb_cntrl_sites(i);%note that this redistributes the same size class distrubtion to all sites and in all circumstances
            this_reef_sites(j,:) = this_reef_cots_4redistrib.*betarndN(j);%used to be current_COTS(i,:) but this just divided the reef-level COTS to all sites
        end

        COTS_sites_density(i,1)={this_reef_sites};
        stored_betarndN(i,1)={betarndN};
    end
end