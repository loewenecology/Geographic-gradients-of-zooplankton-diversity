# Geographic-gradients-of-zooplankton-diversity
Code for calculating zooplankton biodiversity metrics and replicating analysis of elevational and latitudinal diversity gradients by generalized linear mixed effect/multilevel models.

Sitebysp.csv contains zooplankton occurrence data (site rows and species columns);
Sitebygeo.csv contains site geographic coordinate data (site rows and variable columns);
Spbytrait.csv contains species trait data (species rows and trait/taxonomic columns);
Taxalist.RData contains taxonomic ranks of species (RData file); and
Sitebysp.csv, Sitebygeo.csv, Spbytrait.csv, and Taxalist.RData are permanently archived at (url).

Bio.metrics.pub.csv contains processed biodiversity data (site rows and biodiversity metric columns) for use with statistical code. 
Variables include:
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
FRic.runs.ele (number of randomizations used in calculations for standardized effect sizes of functional richness in elevational zones).
