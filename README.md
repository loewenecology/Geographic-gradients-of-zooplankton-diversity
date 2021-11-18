# Geographic-gradients-of-zooplankton-diversity
Code and data for calculating zooplankton biodiversity metrics and replicating analysis of elevational and latitudinal diversity gradients by generalized linear mixed effect/multilevel models.

Sitebysp.csv contains zooplankton occurrence data (site rows and species columns);
Sitebygeo.csv contains site geographic coordinate data (site rows and variable columns);
Spbytrait.csv contains species trait data (species rows and trait/taxonomic columns); and
Taxalist.RData contains taxonomic ranks of species (RData file).

Sitebysp.csv, Sitebygeo.csv, Spbytrait.csv, and Taxalist.RData are permanently archived at (url).

Bio.metrics.pub.csv contains processed biodiversity data (site rows and biodiversity metric columns) for use with statistical code.

Variables include:
Latitude;
Elevation;
lat_class (site latitudinal zone);
elev_class (site elevational zone);
Taxo.alpha (species richness);
Glob.nbsp (number of taxa excluding juveniles);
Glob.FRic (functional richness standardized between zero and one);
Glob.FDis (functional dispersion);
CWM.Bl (community mean body length);
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
FDis.obs.p.ele (P-values for observed value used in calculations for standardized effect sizes of functional dispersion in elevational zones); and
FDis.runs.ele (number of randomizations used in calculations for standardized effect sizes of functional dispersion in elevational zones).


ClimateNA_input_elev_1964-2015Y.csv contains annual climate data across the entire study period (1964-2015) at each sampling location obtained using ClimateNA v6.40 (Wang et al. 2016).

Wang, T., Hamann, A., Spittlehouse, D. & Carroll, C. (2016). Locally downscaled and spatially customizable climate data for historical and future periods for North America. PLoS One, 11, e0156720.

MAT refers to mean annual temperature (°C) and
TD refers to temperature difference between mean warmest month temperature and mean coldest month temperature, or continentality (°C).
