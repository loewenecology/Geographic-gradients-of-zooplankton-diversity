# Geographic-gradients-of-zooplankton-diversity
Code and data for calculating zooplankton biodiversity metrics and replicating analysis of elevational and latitudinal diversity gradients by generalized linear mixed effect/multilevel models.

Sitebysp.csv contains zooplankton occurrence data (site rows and species columns).

Sitebygeo.csv contains site geographic coordinate data (site rows and variable columns), including: ID (lake ID); 
Lake (lake name); 
Elevation (masl); 
Latitude (N); 
Longitude (E); 
Area (km2); 
Depth (m); 
Sampling.events (estimated/presumed number of sampling events at each site); and 
Sampling.years (estimated/presumed number of years sampled at each site).

Spbytrait.csv contains species trait data (species rows and trait/taxonomic columns), including: Rownames (taxa names); 
(Sp (taxa name); 
Sp (taxa name); 
Taxa (taxa name with space); 
Level (taxonomic level); 
Order (taxonomic order); 
Length (taxa length, mm); 
Length.ref (taxa length reference, see Table S1 in manuscript); 
Length.taxa (length representative taxa); 
Feed.ord (feeding guild ordinal representation, 1 = substrate-grazing, 2 = seston-filtering, 3 = stationary suspension-feeding with occasional grasping, and 4 = raptorial-feeding); and 
Full (full taxa name).

coefficients.bio.metrics.csv contains data for plotting coefficients for biodiversity by elevation/latitude models, including: Metric (biodiversity metric); 
Group (grouping); 
Parameter (Coefficient); 
Class (elevational/latitudinal zone); 
Estimate (metric estimate); 
Est.Error (metric estimate error); 
Q2.5 (2.5% quantile); and
Q97.5 (97.5 % quantile).

coefficients.clim.csv contains data for plotting coefficients for biodiversity by climate modelsi, ncluding: Metric (biodiversity metric); 
Climate (climate variable); 
Group (grouping); 
Parameter (Coefficient); 
Class (elevational/latitudinal zone); 
Estimate (metric estimate); 
Est.Error (metric estimate error); 
Q2.5 (2.5% quantile); and 
Q97.5 (97.5 % quantile).

Taxalist.rds contains taxonomic ranks of species (rds file).

Bio.metrics.pub.csv contains processed biodiversity data (site rows and biodiversity metric columns) for use with statistical code.

Variables include:
ID (lake ID); 
Lake (lake name);
Latitude (N);
Elevation (masl);
lat_class (site latitudinal zone);
elev_class (site elevational zone);
Area (km2);
Depth (m);
Sampling.events (estimated/presumed number of sampling events at each site);
Sampling.years (estimated/presumed number of years sampled at each site);
Taxo.alpha (species richness);
Glob.nbsp (number of taxa excluding juveniles);
Glob.FRic (functional richness standardized between zero and one);
Glob.FDis (functional dispersion);
CWM.Bl (community mean body length, mm);
CWM.Feed (community mean feeding guilds);
Glob.PD (phylogenetic richness);
Glob.MPD (phylogenetic mean pairwise distance);
Glob.MPD.beta (phylogenetic mean pairwise distance divided by 100 to fit bounds of zero and one);
pd.ntaxa.lat (number of taxa used in calculations for standardized effect sizes of phylogenetic richness in latitudinal zones);
pd.obs.lat (observed value used in calculations for standardized effect sizes of phylogenetic richness in latitudinal zones);
pd.rand.mean.lat (mean of null communities used in calculations for standardized effect sizes of phylogenetic richness in latitudinal zones);
pd.rand.sd.lat (standard deviation of null communities used in calculations for standardized effect sizes of phylogenetic richness in latitudinal zones);
pd.obs.rank.lat (rank of observed value used in calculations for standardized effect sizes of phylogenetic richness in latitudinal zones);
pd.obs.z.lat (standardized effect sizes of phylogenetic richness in latitudinal zones);
pd.obs.p.lat (P-values for observed value used in calculations for standardized effect sizes of phylogenetic richness in latitudinal zones);
pd.runs.lat (number of randomizations used in calculations for standardized effect sizes of phylogenetic richness in latitudinal zones);
pd.ntaxa.ele (number of taxa used in calculations for standardized effect sizes of phylogenetic richness in elevational zones);
pd.obs.ele (observed value used in calculations for standardized effect sizes of phylogenetic richness in elevational zones);
pd.rand.mean.ele (mean of null communities used in calculations for standardized effect sizes of phylogenetic richness in elevational zones);
pd.rand.sd.ele (standard deviation of null communities used in calculations for standardized effect sizes of phylogenetic richness in elevational zones);
pd.obs.rank.ele (rank of observed value used in calculations for standardized effect sizes of phylogenetic richness in elevational zones);
pd.obs.z.ele (standardized effect sizes of phylogenetic richness in elevational zones);
pd.obs.p.ele (P-values for observed value used in calculations for standardized effect sizes of phylogenetic richness in elevational zones);
pd.runs.ele (number of randomizations used in calculations for standardized effect sizes of phylogenetic richness in elevational zones);
mpd.ntaxa.lat (number of taxa used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in latitudinal zones);
mpd.obs.lat (observed value used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in latitudinal zones);
mpd.rand.mean.lat (mean of null communities used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in latitudinal zones);
mpd.rand.sd.lat (standard deviation of null communities used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in latitudinal zones);
mpd.obs.rank.lat (rank of observed value used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in latitudinal zones);
mpd.obs.z.lat (standardized effect sizes of phylogenetic mean pairwise distance in latitudinal zones);
mpd.obs.p.lat (P-values for observed value used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in latitudinal zones);
mpd.runs.lat (number of randomizations used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in latitudinal zones);
mpd.ntaxa.ele (number of taxa used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in elevational zones);
mpd.obs.ele (observed value used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in elevational zones);
mpd.rand.mean.ele (mean of null communities used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in elevational zones);
mpd.rand.sd.ele (standard deviation of null communities used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in elevational zones);
mpd.obs.rank.ele (rank of observed value used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in elevational zones);
mpd.obs.z.ele (standardized effect sizes of phylogenetic mean pairwise distance in elevational zones);
mpd.obs.p.ele (P-values for observed value used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in elevational zones);
mpd.runs.ele (number of randomizations used in calculations for standardized effect sizes of phylogenetic mean pairwise distance in elevational zones);
FRic.ntaxa.lat (number of taxa used in calculations for standardized effect sizes of functional richness in latitudinal zones);
FRic.obs.lat (observed value used in calculations for standardized effect sizes of functional richness in latitudinal zones);
FRic.rand.mean.lat (mean of null communities used in calculations for standardized effect sizes of functional richness in latitudinal zones);
FRic.rand.sd.lat (standard deviation of null communities used in calculations for standardized effect sizes of functional richness in latitudinal zones);
FRic.obs.rank.lat (rank of observed value used in calculations for standardized effect sizes of functional richness in latitudinal zones);
FRic.obs.z.lat (standardized effect sizes of functional richness in latitudinal zones);
FRic.obs.p.lat (P-values for observed value used in calculations for standardized effect sizes of functional richness in latitudinal zones);
FRic.runs.lat (number of randomizations used in calculations for standardized effect sizes of functional richness in latitudinal zones);
FRic.ntaxa.ele (number of taxa used in calculations for standardized effect sizes of functional richness in elevational zones);
FRic.obs.ele (observed value used in calculations for standardized effect sizes of functional richness in elevational zones);
FRic.rand.mean.ele (mean of null communities used in calculations for standardized effect sizes of functional richness in elevational zones);
FRic.rand.sd.ele (standard deviation of null communities used in calculations for standardized effect sizes of functional richness in elevational zones);
FRic.obs.rank.ele (rank of observed value used in calculations for standardized effect sizes of functional richness in elevational zones);
FRic.obs.z.ele (standardized effect sizes of functional richness in elevational zones);
FRic.obs.p.ele (P-values for observed value used in calculations for standardized effect sizes of functional richness in elevational zones); and
FRic.runs.ele (number of randomizations used in calculations for standardized effect sizes of functional richness in elevational zones);
FDis.ntaxa.lat (number of taxa used in calculations for standardized effect sizes of functional dispersion in latitudinal zones);
FDis.obs.lat (observed value used in calculations for standardized effect sizes of functional dispersion in latitudinal zones);
FDis.rand.mean.lat (mean of null communities used in calculations for standardized effect sizes of functional dispersion in latitudinal zones);
FDis.rand.sd.lat (standard deviation of null communities used in calculations for standardized effect sizes of functional dispersion in latitudinal zones);
FDis.obs.rank.lat (rank of observed value used in calculations for standardized effect sizes of functional dispersion in latitudinal zones);
FDis.obs.z.lat (standardized effect sizes of functional dispersion in latitudinal zones);
FDis.obs.p.lat (P-values for observed value used in calculations for standardized effect sizes of functional dispersion in latitudinal zones);
FDis.runs.lat (number of randomizations used in calculations for standardized effect sizes of functional dispersion in latitudinal zones);
FDis.ntaxa.ele (number of taxa used in calculations for standardized effect sizes of functional dispersion in elevational zones);
FDis.obs.ele (observed value used in calculations for standardized effect sizes of functional dispersion in elevational zones);
FDis.rand.mean.ele (mean of null communities used in calculations for standardized effect sizes of functional dispersion in elevational zones);
FDis.rand.sd.ele (standard deviation of null communities used in calculations for standardized effect sizes of functional dispersion in elevational zones);
FDis.obs.rank.ele (rank of observed value used in calculations for standardized effect sizes of functional dispersion in elevational zones);
FDis.obs.z.ele (standardized effect sizes of functional dispersion in elevational zones);
FDis.obs.p.ele (P-values for observed value used in calculations for standardized effect sizes of functional dispersion in elevational zones);
FDis.runs.ele (number of randomizations used in calculations for standardized effect sizes of functional dispersion in elevational zones);
pd.obs.p.lat.fix (P-values for observed value used in calculations for phylogenetic richness SES in latitudinal zones replacing values of 1);
pd.obs.p.lat.qnorm (probit-transformed quantile P-values for phylogenetic richness SES in latitudinal zones);
pd.obs.p.ele.fix (P-values for observed value used in calculations for phylogenetic richness SES in elevational zones replacing values of 1);
pd.obs.p.ele.qnorm (probit-transformed quantile P-values for phylogenetic richness SES in elevational zones);
mpd.obs.p.lat.fix (P-values for observed value used in calculations for phylogenetic mean pairwise distance SES in latitudinal zones replacing values of 1);
mpd.obs.p.lat.qnorm (probit-transformed quantile P-values for phylogenetic mean pairwise distance SES in latitudinal zones);
mpd.obs.p.ele.fix (P-values for observed value used in calculations for phylogenetic mean pairwise distance SES in elevational zones replacing values of 1);
mpd.obs.p.ele.qnorm (probit-transformed quantile P-values for phylogenetic mean pairwise distance SES in elevational zones);
FRic.obs.p.lat.fix (P-values for observed value used in calculations for functional richness SES in latitudinal zones replacing values of 1);
FRic.obs.p.lat.qnorm (probit-transformed quantile P-values for functional richness SES in latitudinal zones);
FRic.obs.p.ele.fix (P-values for observed value used in calculations for functional richness SES in elevational zones replacing values of 1);
FRic.obs.p.ele.qnorm (probit-transformed quantile P-values for functional richness SES in elevational zones);
FDis.obs.p.lat.fix (P-values for observed value used in calculations for functional dispersion SES in latitudinal zones replacing values of 1);
FDis.obs.p.lat.qnorm (probit-transformed quantile P-values for functional dispersion SES in latitudinal zones);
FDis.obs.p.ele.fix (P-values for observed value used in calculations for functional dispersion SES in elevational zones replacing values of 1);
FDis.obs.p.ele.qnorm (probit-transformed quantile P-values for functional dispersion SES in elevational zones);
pd.ntaxa.glob (number of taxa used in calculations for standardized effect sizes of phylogenetic richness with global pools);
pd.obs.glob (observed value used in calculations for standardized effect sizes of phylogenetic richness with global pools);
pd.rand.mean.glob (mean of null communities used in calculations for standardized effect sizes of phylogenetic richness with global pools);
pd.rand.sd.glob (standard deviation of null communities used in calculations for standardized effect sizes of phylogenetic richness with global pools);
pd.obs.rank.glob (rank of observed value used in calculations for standardized effect sizes of phylogenetic richness with global pools);
pd.obs.z.glob (standardized effect sizes of phylogenetic richness with global pools);
pd.obs.p.glob (P-values for observed value used in calculations for standardized effect sizes of phylogenetic richness with global pools);
pd.runs.glob (number of randomizations used in calculations for standardized effect sizes of phylogenetic richness with global pools);
mpd.ntaxa.glob (number of taxa used in calculations for standardized effect sizes of phylogenetic mean pairwise distance with global pools);
mpd.obs.glob (observed value used in calculations for standardized effect sizes of phylogenetic mean pairwise distance with global pools);
mpd.rand.mean.glob (mean of null communities used in calculations for standardized effect sizes of phylogenetic mean pairwise distance with global pools);
mpd.rand.sd.glob (standard deviation of null communities used in calculations for standardized effect sizes of phylogenetic mean pairwise distance with global pools);
mpd.obs.rank.glob (rank of observed value used in calculations for standardized effect sizes of phylogenetic mean pairwise distance with global pools);
mpd.obs.z.glob (standardized effect sizes of phylogenetic mean pairwise distance with global pools);
mpd.obs.p.glob (P-values for observed value used in calculations for standardized effect sizes of phylogenetic mean pairwise distance with global pools);
mpd.runs.glob (number of randomizations used in calculations for standardized effect sizes of phylogenetic mean pairwise distance with global pools);
FRic.ntaxa.glob (number of taxa used in calculations for standardized effect sizes of functional richness with global pools);
FRic.obs.glob (observed value used in calculations for standardized effect sizes of functional richness with global pools);
FRic.rand.mean.glob (mean of null communities used in calculations for standardized effect sizes of functional richness with global pools);
FRic.rand.sd.glob (standard deviation of null communities used in calculations for standardized effect sizes of functional richness with global pools);
FRic.obs.rank.glob (rank of observed value used in calculations for standardized effect sizes of functional richness with global pools);
FRic.obs.z.glob (standardized effect sizes of functional richness with global pools);
FRic.obs.p.glob (P-values for observed value used in calculations for standardized effect sizes of functional richness with global pools);
FRic.runs.glob (number of randomizations used in calculations for standardized effect sizes of functional richness with global pools);
FDis.ntaxa.glob (number of taxa used in calculations for standardized effect sizes of functional dispersion with global pools);
FDis.obs.glob (observed value used in calculations for standardized effect sizes of functional dispersion with global pools);
FDis.rand.mean.glob (mean of null communities used in calculations for standardized effect sizes of functional dispersion with global pools);
FDis.rand.sd.glob (standard deviation of null communities used in calculations for standardized effect sizes of functional dispersion with global pools);
FDis.obs.rank.glob (rank of observed value used in calculations for standardized effect sizes of functional dispersion with global pools);
FDis.obs.z.glob (standardized effect sizes of functional dispersion with global pools);
FDis.obs.p.glob (P-values for observed value used in calculations for standardized effect sizes of functional dispersion with global pools);
FDis.runs.glob (number of randomizations used in calculations for standardized effect sizes of functional dispersion with global pools);
pd.obs.p.glob.fix (P-values for observed value used in calculations for phylogenetic richness SES with global pools replacing values of 1);
pd.obs.p.glob.qnorm (probit-transformed quantile P-values for phylogenetic richness SES with global pools);
mpd.obs.p.glob.fix (P-values for observed value used in calculations for phylogenetic mean pairwise distance SES with global pools replacing values of 1);
mpd.obs.p.glob.qnorm (probit-transformed quantile P-values for phylogenetic mean pairwise distance SES with global pools);
FRic.obs.p.glob.fix (P-values for observed value used in calculations for functional richness SES with global pools replacing values of 1);
FRic.obs.p.glob.qnorm (probit-transformed quantile P-values for functional richness SES with global pools);
FDis.obs.p.glob.fix (P-values for observed value used in calculations for functional dispersion SES with global pools replacing values of 1); and
FDis.obs.p.glob.qnorm (probit-transformed quantile P-values for functional dispersion SES with global pools).

ClimateNA_input_elev_1964-2015Y.csv contains annual climate data across the entire study period (1964-2015) at each sampling location obtained using ClimateNA v6.40 (Wang et al. 2016).

MAT refers to mean annual temperature (°C) and TD refers to temperature difference between mean warmest month temperature and mean coldest month temperature, or continentality (°C). 
Year (year); ID1 (lake ID); ID2 (lake region); Latitude (N); Longitude (E); and Elevation (masl). See Wang et al. (2016) for details of other variables.

REFERENCES

Wang, T., Hamann, A., Spittlehouse, D. & Carroll, C. (2016). Locally downscaled and spatially customizable climate data for historical and future periods for North America. PLoS One, 11, e0156720.
