#........................................................................................................#
#..................Biodiversity patterns diverge along geographic temperature gradients..................#
#........................................................................................................#

#........................................................................................................#
#..................Author of the script: Charlie Loewen..................................................#
#..................Date latest modifications: 2022-01-05.................................................#
#........................................................................................................#

##########################################################################################################
################## For calculating biodiversity metrics and evaluating elevational and latitudinal
################## relationships by Bayesian generalized linear mixed effect/multilevel models
##########################################################################################################

##################
# Load R packages
##################

library(tidyverse)
library(classInt)
library(FD)
library(betapart)
library(picante)
library(taxize)
library(ggplot2)
library(grid)
library(gtable)
library(ggmap)
library(ggsn)
library(ggridges)
library(brms)
library(ggeffects)

##################
# Load data in R environment
##################

taxo.alpha = read.table("Sitebysp.csv", header = T, sep = ",", row.names = 1) #matrix of sites (rows) by species (columns)
taxo.beta <- taxo.alpha[1:119] #excludes UNID cyclopoids, calanaoids, and cladocera

sites = read.table("Sitebygeo.csv", header = T, sep = ",") #matrix of sites (rows) by geographic constraints
rownames(sites) <- sites$ID

traits = read.table("Spbytrait.csv", header = T, sep = ",", row.names = 1) #matrix of species (rows) by traits (columns)
taxa.list = readRDS("Taxalist.RData") #taxonomic rankings

##################
# Mapping natural breaks
##################

# Find natural breaks in geographic distributions
lat.classes <- classIntervals(sites$Latitude, style = "fisher")
elev.classes <- classIntervals(sites$Elevation, style = "fisher")

clean.env.lat <- sites %>%
  mutate(lat_class = cut(Latitude, lat.classes$brks, include.lowest = T))
clean.env.elev <- sites %>%
  mutate(elev_class = cut(Elevation, elev.classes$brks, include.lowest = T))

# Assign latitudinal zones
clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "[36.6,39.3]", "36.6-39.3", as.character(lat_class)))

clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "(39.3,41.2]", "39.3-41.2", as.character(lat_class)))

clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "(41.2,43.2]", "41.2-43.2", as.character(lat_class)))

clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "(43.2,44.6]", "43.2-44.6", as.character(lat_class)))

clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "(44.6,46.3]", "44.6-46.3", as.character(lat_class)))

clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "(46.3,48]", "46.3-48.0", as.character(lat_class)))

clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "(48,50]", "48.0-50.0", as.character(lat_class)))

clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "(50,52.1]", "50.0-52.1", as.character(lat_class)))

clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "(52.1,54.2]", "52.1-54.2", as.character(lat_class)))

clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "(54.2,58.1]", "54.2-58.1", as.character(lat_class)))

clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "(58.1,62.1]", "58.1-62.1", as.character(lat_class)))

clean.env.lat <- clean.env.lat %>%
  mutate(lat_class = ifelse(as.character(lat_class) == "(62.1,66.2]", "62.1-66.2", as.character(lat_class)))

clean.env.lat$lat_class <- as.factor(clean.env.lat$lat_class)
clean.env.lat$lat_class <- factor(clean.env.lat$lat_class, levels = c("36.6-39.3",
                                                                      "39.3-41.2",
                                                                      "41.2-43.2",
                                                                      "43.2-44.6",
                                                                      "44.6-46.3",
                                                                      "46.3-48.0",
                                                                      "48.0-50.0",
                                                                      "50.0-52.1",
                                                                      "52.1-54.2",
                                                                      "54.2-58.1",
                                                                      "58.1-62.1",
                                                                      "62.1-66.2"))
rownames(clean.env.lat) <- clean.env.lat$ID

# Set colour palete
a60505_6e068f<-c("#a60505", "#c5600a", "#e2d110", "#a2ec29", "#67eb4d", "#70ec95", "#50ebbe", "#30d4ea", "#1476e5", "#0e1ac9", "#420aac", "#6e068f")

# Plot latitudinal zones
(p1.sites.ele <- ggplot() +
    geom_vline(xintercept = elev.classes$brks, color = "grey", size = 0.25) +
    geom_hline(yintercept = lat.classes$brks, color = "grey", size = 0.25) +
    geom_point(data = clean.env.lat, aes(Elevation, Latitude, colour = lat_class), alpha = 0.3) +
    
    scale_x_continuous("Elevation (masl)", limits = c(0, 3750), breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    scale_y_continuous("Latitude (°N)", limits = c(36, 66.17637), breaks = c(36, 42, 48, 54, 60, 66)) +
    
    scale_colour_manual(name = "Latitude", values=a60505_6e068f) +
    
    theme(strip.background = element_blank(),
          strip.text.y = element_text(colour = "black", size = 9),
          strip.text.x = element_blank(),
          
          panel.border = element_rect(colour = "black", fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", fill = "white", size = 0.75),
          panel.spacing = unit(0.12, "cm"),
          
          axis.text.x = element_text(colour = "black", size = 8),
          axis.text.y = element_text(colour = "black", size = 8),
          axis.ticks = element_line(colour = "black", size = 0.4),
          axis.title.y = element_text(colour = "black", size = 10),
          axis.title.x = element_text(colour = "black", size = 10),
          
          legend.key = element_rect("white"),
          legend.key.size = unit(0.4, "cm"),
          legend.justification = c("right", "top"),
          legend.position = c(0.99, 0.985),
          legend.direction = "horizontal",
          legend.text = element_text(size = 8, colour = "black"),
          legend.title = element_text(size = 10, colour = "black"),
          legend.box.background = element_rect(colour = "black", size = 0.8),
          legend.spacing.y = unit(0, 'cm'),
          legend.margin = margin(1, 2, 0, 1)))

# Assign elevational zones
clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "[0,258]", "0-258", as.character(elev_class)))

clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "(258,599]", "258-599", as.character(elev_class)))

clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "(599,856]", "599-856", as.character(elev_class)))

clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "(856,1.1e+03]", "856-1100", as.character(elev_class)))

clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "(1.1e+03,1.32e+03]", "1100-1320", as.character(elev_class)))

clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "(1.32e+03,1.5e+03]", "1320-1500", as.character(elev_class)))

clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "(1.5e+03,1.67e+03]", "1500-1670", as.character(elev_class)))

clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "(1.67e+03,1.87e+03]", "1670-1870", as.character(elev_class)))

clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "(1.87e+03,2.09e+03]", "1870-2090", as.character(elev_class)))

clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "(2.09e+03,2.33e+03]", "2090-2330", as.character(elev_class)))

clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "(2.33e+03,2.78e+03]", "2330-2780", as.character(elev_class)))

clean.env.elev <- clean.env.elev %>%
  mutate(elev_class = ifelse(as.character(elev_class) == "(2.78e+03,3.74e+03]", "2780-3740", as.character(elev_class)))

clean.env.elev$elev_class <- as.factor(clean.env.elev$elev_class)
clean.env.elev$elev_class <- factor(clean.env.elev$elev_class, levels = c("0-258", 
                                                                          "258-599", 
                                                                          "599-856", 
                                                                          "856-1100", 
                                                                          "1100-1320", 
                                                                          "1320-1500", 
                                                                          "1500-1670",
                                                                          "1670-1870",
                                                                          "1870-2090",
                                                                          "2090-2330",
                                                                          "2330-2780",
                                                                          "2780-3740"), ordered = TRUE)
rownames(clean.env.elev) <- clean.env.elev$ID

# Set colour palete
a60505_6e068f<-c("#a60505", "#c5600a", "#e2d110", "#a2ec29", "#67eb4d", "#70ec95", "#50ebbe", "#30d4ea", "#1476e5", "#0e1ac9", "#420aac", "#6e068f")

# Plot elevational zones 
(p1.sites.lat <- ggplot() +
   geom_hline(yintercept = elev.classes$brks, color = "grey", size = 0.25) +
   geom_vline(xintercept = lat.classes$brks, color = "grey", size = 0.25) +
   geom_point(data = clean.env.elev, aes(Latitude, Elevation, colour = elev_class), alpha = 0.3) +
   
   scale_y_continuous("Elevation (masl)", limits = c(0, 3737.901), breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
   scale_x_continuous("Latitude (°N)", limits = c(36, 66.17637), breaks = c(36, 42, 48, 54, 60, 66)) +
   
   scale_colour_manual(name = "Latitude", values = a60505_6e068f) +
   
   theme(strip.background = element_blank(),
         strip.text.y = element_text(colour = "black", size = 9),
         strip.text.x = element_blank(),
         
         panel.border = element_rect(colour = "black", fill = NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black", fill = "white", size = 0.75),
         panel.spacing = unit(0.12, "cm"),
         
         axis.text.x = element_text(colour = "black", size = 8),
         axis.text.y = element_text(colour = "black", size = 8),
         axis.ticks = element_line(colour = "black", size = 0.4),
         axis.title.y = element_text(colour = "black", size = 10),
         axis.title.x = element_text(colour = "black", size = 10),
         
         legend.key = element_rect("white"),
         legend.key.size = unit(0.4, "cm"),
         legend.justification = c("right", "top"),
         legend.position = c(0.99, 0.985),
         legend.direction = "horizontal",
         legend.text = element_text(size = 8, colour = "black"),
         legend.title = element_text(size = 10, colour = "black"),
         legend.box.background = element_rect(colour = "black", size = 0.8),
         legend.spacing.y = unit(0, 'cm'),
         legend.margin = margin(1, 2, 0, 1)))

# Combining latitudinal and elevational zones (Figure 1b in manuscript)
# requires some minor revisions (e.g. adjusting legend) in a graphics editor
p1.sites.ele <- ggplotGrob(p1.sites.ele)
p1.sites.lat<- ggplotGrob(p1.sites.lat)

g1 <- rbind(p1.sites.lat, p1.sites.ele)

g1$widths <- unit.pmax(p1.sites.lat$widths, p1.sites.ele$widths)

grid.newpage()
grid.draw(g1)

# Acquire and plot Basemaps (Figure 1a in manuscript)
# requires some minor revisions (e.g. overlaying plots) in a graphics editor
basemap <- get_stamenmap( bbox = c(left = -139.75116, bottom = 35.60120, right = -112.68687, top = 66.57637), zoom = 4, maptype = "terrain-background")

clean.env.elev <- clean.env.elev[order(clean.env.elev$elev_class), ]

(samplingmap <- ggmap(basemap) +
    geom_point(data = clean.env.elev, aes(Longitude, Latitude, color = elev_class), shape = 16, size = 3) +
    scale_x_continuous(breaks = c(-130, -120),labels = c("-130", "-120")) +
    scale_y_continuous(breaks = (c(36, 42, 48, 54, 60, 66)),
                       labels = c("36", "42", "48", "54", "60", "66")) +
    scalebar(x.min = attr(basemap, "bb")[[2]],
             y.min = attr(basemap, "bb")[[1]],
             x.max = attr(basemap, "bb")[[4]],
             y.max = attr(basemap, "bb")[[3]],
             dist = 500, anchor = c(x = -138.75116, y = 37.5),
             transform = T, location = "bottomleft", st.size = 2.5, st.dist = 0.032, dist_unit = "km") +
    theme(plot.title = element_text(colour = "orange"), 
          panel.border = element_rect(colour = "black", fill = NA)) +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f)) #requires some minor revisions in a graphics editor

basemap2 <- get_stamenmap( bbox = c(left = -168, bottom = 25, right = -85.5, top = 72), zoom = 4, maptype = "toner-background")

(samplingmap2 <- ggmap(basemap2) + 
    theme(plot.title = element_text(colour = "orange"), 
          panel.border = element_rect(colour = "black", fill = NA)))

samplingmap2 + annotate("rect", xmin = -139.75116, xmax = -112.68687, ymin = 35.60120, ymax = 66.57637,
                        alpha = .2, color = "red", size = 2) #requires some minor revisions in a graphics editor

##################
# Mapping species latitudinal and elevational distributions
##################

# Prepare data for plotting site and species distributions
sp.distributions <- taxo.alpha
sp.distributions$SITE <- rownames(sp.distributions)

colnames(sp.distributions)[which(names(sp.distributions) == "Acanthocyclops.capillatus")] <- "Acanthocyclops capillatus"
colnames(sp.distributions)[which(names(sp.distributions) == "Acanthocyclops.vernalis")] <- "Acanthocyclops vernalis"
colnames(sp.distributions)[which(names(sp.distributions) == "Acanthodiaptomus.denticornis")] <- "Acanthodiaptomus denticornis"
colnames(sp.distributions)[which(names(sp.distributions) == "Acroperus.harpae")] <- "Acroperus harpae"
colnames(sp.distributions)[which(names(sp.distributions) == "Aglaodiaptomus.forbesi")] <- "Aglaodiaptomus forbesi"
colnames(sp.distributions)[which(names(sp.distributions) == "Aglaodiaptomus.leptopus")] <- "Aglaodiaptomus leptopus"
colnames(sp.distributions)[which(names(sp.distributions) == "Aglaodiaptomus.lintoni")] <- "Aglaodiaptomus lintoni"
colnames(sp.distributions)[which(names(sp.distributions) == "Alona")] <- "Alona sp."
colnames(sp.distributions)[which(names(sp.distributions) == "Alona.affinis")] <- "Alona affinis"
colnames(sp.distributions)[which(names(sp.distributions) == "Alona.circumfimbriata")] <- "Alona circumfimbriata"
colnames(sp.distributions)[which(names(sp.distributions) == "Alona.costata")] <- "Alona costata"
colnames(sp.distributions)[which(names(sp.distributions) == "Alona.guttata")] <- "Alona guttata"
colnames(sp.distributions)[which(names(sp.distributions) == "Alona.intermedia")] <- "Alona intermedia"
colnames(sp.distributions)[which(names(sp.distributions) == "Alona.quadrangularis")] <- "Alona quadrangularis"
colnames(sp.distributions)[which(names(sp.distributions) == "Alonella.excisa")] <- "Alonella excisa"
colnames(sp.distributions)[which(names(sp.distributions) == "Alonella.nana")] <- "Alonella nana"
colnames(sp.distributions)[which(names(sp.distributions) == "Arctodiaptomus.arapahoensis")] <- "Arctodiaptomus arapahoensis"
colnames(sp.distributions)[which(names(sp.distributions) == "Bosmina.coregoni")] <- "Bosmina coregoni"
colnames(sp.distributions)[which(names(sp.distributions) == "Bosmina.hagmanni")] <- "Bosmina hagmanni"
colnames(sp.distributions)[which(names(sp.distributions) == "Bosmina.longirostris")] <- "Bosmina longirostris"
colnames(sp.distributions)[which(names(sp.distributions) == "Camptocercus.rectirostris")] <- "Camptocercus rectirostris"
colnames(sp.distributions)[which(names(sp.distributions) == "Ceriodaphnia")] <- "Ceriodaphnia sp."
colnames(sp.distributions)[which(names(sp.distributions) == "Ceriodaphnia.acanthina")] <- "Ceriodaphnia acanthina"
colnames(sp.distributions)[which(names(sp.distributions) == "Ceriodaphnia.lacustris")] <- "Ceriodaphnia lacustris"
colnames(sp.distributions)[which(names(sp.distributions) == "Ceriodaphnia.pulchella")] <- "Ceriodaphnia pulchella"
colnames(sp.distributions)[which(names(sp.distributions) == "Ceriodaphnia.quadrangula")] <- "Ceriodaphnia quadrangula"
colnames(sp.distributions)[which(names(sp.distributions) == "Ceriodaphnia.reticulata")] <- "Ceriodaphnia reticulata"
colnames(sp.distributions)[which(names(sp.distributions) == "Chydorus")] <- "Chydorus sp."
colnames(sp.distributions)[which(names(sp.distributions) == "Chydorus.ovalis")] <- "Chydorus ovalis"
colnames(sp.distributions)[which(names(sp.distributions) == "Chydorus.sphaericus")] <- "Chydorus sphaericus"
colnames(sp.distributions)[which(names(sp.distributions) == "Coronatella.rectangula")] <- "Coronatella rectangula"
colnames(sp.distributions)[which(names(sp.distributions) == "Cyclops")] <- "Cyclops sp."
colnames(sp.distributions)[which(names(sp.distributions) == "Cyclops.scutifer")] <- "Cyclops scutifer"
colnames(sp.distributions)[which(names(sp.distributions) == "Daphnia")] <- "Daphnia sp."
colnames(sp.distributions)[which(names(sp.distributions) == "Daphnia.ambigua")] <- "Daphnia ambigua"
colnames(sp.distributions)[which(names(sp.distributions) == "Daphnia.catawba")] <- "Daphnia catawba"
colnames(sp.distributions)[which(names(sp.distributions) == "Daphnia.dentifera")] <- "Daphnia dentifera"
colnames(sp.distributions)[which(names(sp.distributions) == "Daphnia.galeata")] <- "Daphnia galeata"
colnames(sp.distributions)[which(names(sp.distributions) == "Daphnia.longiremis")] <- "Daphnia longiremis"
colnames(sp.distributions)[which(names(sp.distributions) == "Daphnia.middendorffiana")] <- "Daphnia middendorffiana"
colnames(sp.distributions)[which(names(sp.distributions) == "Daphnia.pulex")] <- "Daphnia pulex"
colnames(sp.distributions)[which(names(sp.distributions) == "Daphnia.schoedleri")] <- "Daphnia schoedleri"
colnames(sp.distributions)[which(names(sp.distributions) == "Daphnia.similis")] <- "Daphnia similis"
colnames(sp.distributions)[which(names(sp.distributions) == "Diacyclops.navus")] <- "Diacyclops navus"
colnames(sp.distributions)[which(names(sp.distributions) == "Diacyclops.thomasi")] <- "Diacyclops thomasi"
colnames(sp.distributions)[which(names(sp.distributions) == "Diaphanosoma.birgei")] <- "Diaphanosoma birgei"
colnames(sp.distributions)[which(names(sp.distributions) == "Diaptomus")] <- "Diaptomus sp."
colnames(sp.distributions)[which(names(sp.distributions) == "Ectocyclops")] <- "Ectocyclops sp."
colnames(sp.distributions)[which(names(sp.distributions) == "Epischura.nevadensis")] <- "Epischura nevadensis"
colnames(sp.distributions)[which(names(sp.distributions) == "Eucyclops.agilis")] <- "Eucyclops agilis"
colnames(sp.distributions)[which(names(sp.distributions) == "Eucyclops.elegans")] <- "Eucyclops elegans"
colnames(sp.distributions)[which(names(sp.distributions) == "Eurycercus.lamellatus")] <- "Eurycercus lamellatus"
colnames(sp.distributions)[which(names(sp.distributions) == "Graptoleberis.testudinaria")] <- "Graptoleberis testudinaria"
colnames(sp.distributions)[which(names(sp.distributions) == "Hesperodiaptomus.arcticus")] <- "Hesperodiaptomus arcticus"
colnames(sp.distributions)[which(names(sp.distributions) == "Hesperodiaptomus.eiseni")] <- "Hesperodiaptomus eiseni"
colnames(sp.distributions)[which(names(sp.distributions) == "Hesperodiaptomus.franciscanus")] <- "Hesperodiaptomus franciscanus"
colnames(sp.distributions)[which(names(sp.distributions) == "Hesperodiaptomus.kenai")] <- "Hesperodiaptomus kenai"
colnames(sp.distributions)[which(names(sp.distributions) == "Hesperodiaptomus.shoshone")] <- "Hesperodiaptomus shoshone"
colnames(sp.distributions)[which(names(sp.distributions) == "Heterocope.septentrionalis")] <- "Heterocope septentrionalis"
colnames(sp.distributions)[which(names(sp.distributions) == "Holopedium.gibberum")] <- "Holopedium gibberum"
colnames(sp.distributions)[which(names(sp.distributions) == "Kurzia.latissima")] <- "Kurzia latissima"
colnames(sp.distributions)[which(names(sp.distributions) == "Latona.setifera")] <- "Latona setifera"
colnames(sp.distributions)[which(names(sp.distributions) == "Leptodiaptomus.angustilobus")] <- "Leptodiaptomus angustilobus"
colnames(sp.distributions)[which(names(sp.distributions) == "Leptodiaptomus.ashlandi")] <- "Leptodiaptomus ashlandi"
colnames(sp.distributions)[which(names(sp.distributions) == "Leptodiaptomus.connexus")] <- "Leptodiaptomus connexus"
colnames(sp.distributions)[which(names(sp.distributions) == "Leptodiaptomus.cuauhtemoci")] <- "Leptodiaptomus cuauhtemoci"
colnames(sp.distributions)[which(names(sp.distributions) == "Leptodiaptomus.novamexicanus")] <- "Leptodiaptomus novamexicanus"
colnames(sp.distributions)[which(names(sp.distributions) == "Leptodiaptomus.nudus")] <- "Leptodiaptomus nudus"
colnames(sp.distributions)[which(names(sp.distributions) == "Leptodiaptomus.sicilis")] <- "Leptodiaptomus sicilis"
colnames(sp.distributions)[which(names(sp.distributions) == "Leptodiaptomus.siciloides")] <- "Leptodiaptomus siciloides"
colnames(sp.distributions)[which(names(sp.distributions) == "Leptodiaptomus.signicauda")] <- "Leptodiaptomus signicauda"
colnames(sp.distributions)[which(names(sp.distributions) == "Leptodiaptomus.tyrrelli")] <- "Leptodiaptomus tyrrelli"
colnames(sp.distributions)[which(names(sp.distributions) == "Leptodora.kindtii")] <- "Leptodora kindtii"
colnames(sp.distributions)[which(names(sp.distributions) == "Leydigia.leydigi")] <- "Leydigia leydigi"
colnames(sp.distributions)[which(names(sp.distributions) == "Macrocyclops")] <- "Macrocyclops sp."
colnames(sp.distributions)[which(names(sp.distributions) == "Macrocyclops.albidus")] <- "Macrocyclops albidus"
colnames(sp.distributions)[which(names(sp.distributions) == "Macrocyclops.fuscus")] <- "Macrocyclops fuscus"
colnames(sp.distributions)[which(names(sp.distributions) == "Macrothrix.hirsuticornis")] <- "Macrothrix hirsuticornis"
colnames(sp.distributions)[which(names(sp.distributions) == "Macrothrix.laticornis")] <- "Macrothrix laticornis"
colnames(sp.distributions)[which(names(sp.distributions) == "Microcyclops.varicans")] <- "Microcyclops varicans"
colnames(sp.distributions)[which(names(sp.distributions) == "Onychodiaptomus.sanguineus")] <- "Onychodiaptomus sanguineus"
colnames(sp.distributions)[which(names(sp.distributions) == "Orthocyclops.modestus")] <- "Orthocyclops modestus"
colnames(sp.distributions)[which(names(sp.distributions) == "Paracyclops.poppei")] <- "Paracyclops poppei"
colnames(sp.distributions)[which(names(sp.distributions) == "Picripleuroxus.denticulatus")] <- "Picripleuroxus denticulatus"
colnames(sp.distributions)[which(names(sp.distributions) == "Pleuroxus")] <- "Pleuroxus sp."
colnames(sp.distributions)[which(names(sp.distributions) == "Polyphemus.pediculus")] <- "Polyphemus pediculus"
colnames(sp.distributions)[which(names(sp.distributions) == "Scapholeberis.kingi")] <- "Scapholeberis kingi"
colnames(sp.distributions)[which(names(sp.distributions) == "Senecella.calanoides")] <- "Senecella calanoides"
colnames(sp.distributions)[which(names(sp.distributions) == "Sida.crystallina")] <- "Sida crystallina"
colnames(sp.distributions)[which(names(sp.distributions) == "Simocephalus.serrulatus")] <- "Simocephalus serrulatus"
colnames(sp.distributions)[which(names(sp.distributions) == "Simocephalus.vetulus")] <- "Simocephalus vetulus"
colnames(sp.distributions)[which(names(sp.distributions) == "Skistodiaptomus.oregonensis")] <- "Skistodiaptomus oregonensis"
colnames(sp.distributions)[which(names(sp.distributions) == "Streblocerus.serricaudatus")] <- "Streblocerus serricaudatus"
colnames(sp.distributions)[which(names(sp.distributions) == "Tropocyclops.prasinus")] <- "Tropocyclops prasinus"

sp.distributions <- merge(sp.distributions, sites["Elevation"], by = 0)
sp.distributions <- sp.distributions[-1]
rownames(sp.distributions) <- sp.distributions$SITE

sp.distributions <- merge(sp.distributions, sites["Latitude"], by = 0)
sp.distributions <- sp.distributions[-1]
rownames(sp.distributions) <- sp.distributions$SITE

sp.distributions <- merge(sp.distributions, sites["Longitude"], by = 0)
sp.distributions <- sp.distributions[-1]
rownames(sp.distributions) <- sp.distributions$SITE

sp.distributions.working <- sp.distributions[1:119]
sp.distributions.elev <- as.data.frame(apply(sp.distributions.working, 2, function(x) sp.distributions[124] * x))
colnames(sp.distributions.elev) <- colnames(sp.distributions.working)

speciesbyelev <- gather(sp.distributions.elev, Taxa, Elevation)

sp.distributions.working <- sp.distributions[1:119]
sp.distributions.lat <- as.data.frame(apply(sp.distributions.working, 2, function(x) sp.distributions[125] * x))
colnames(sp.distributions.lat) <- colnames(sp.distributions.working)

speciesbylat <- gather(sp.distributions.lat, Taxa, Latitude)

combinbedelevlat <- merge(speciesbyelev, speciesbylat["Latitude"], by = 0)
rownames(combinbedelevlat) <- combinbedelevlat$Row.names
combinbedelevlat <- combinbedelevlat[-1]

finaleelvlat <- gather(combinbedelevlat, Geo, Value, Elevation, Latitude)

finaleelvlat[finaleelvlat == 0] <- NA

# Plot ridge lines plots (Figure 2 in manuscript)
# requires some minor revisions (e.g. combining plots) in a graphics editor
finaleelvlat$Taxaa <- factor(finaleelvlat$Taxa, levels = c("Acanthocyclops capillatus",
                                                           "Acanthocyclops vernalis",
                                                           "Acanthodiaptomus denticornis",
                                                           "Acroperus harpae",
                                                           "Aglaodiaptomus forbesi",
                                                           "Aglaodiaptomus leptopus",
                                                           "Aglaodiaptomus lintoni",
                                                           "Alona sp.",
                                                           "Alona affinis",
                                                           "Alona circumfimbriata",
                                                           "Alona costata",
                                                           "Alona guttata",
                                                           "Alona intermedia",
                                                           "Alona quadrangularis",
                                                           "Alonella excisa",
                                                           "Alonella nana",
                                                           "Arctodiaptomus arapahoensis",
                                                           "Bosmina coregoni",
                                                           "Bosmina hagmanni",
                                                           "Bosmina longirostris",
                                                           "Camptocercus rectirostris",
                                                           "Ceriodaphnia sp.",
                                                           "Ceriodaphnia acanthina",
                                                           "Ceriodaphnia lacustris",
                                                           "Ceriodaphnia pulchella",
                                                           "Ceriodaphnia quadrangula",
                                                           "Ceriodaphnia reticulata",
                                                           "Chydorus sp.",
                                                           "Chydorus ovalis",
                                                           "Chydorus sphaericus",
                                                           "Coronatella rectangula",
                                                           "Cyclops sp.",
                                                           "Cyclops scutifer",
                                                           "Daphnia sp.",
                                                           "Daphnia ambigua",
                                                           "Daphnia catawba",
                                                           "Daphnia dentifera",
                                                           "Daphnia galeata",
                                                           "Daphnia longiremis",
                                                           "Daphnia middendorffiana",
                                                           "Daphnia pulex",
                                                           "Daphnia schoedleri",
                                                           "Daphnia similis",
                                                           "Diacyclops navus",
                                                           "Diacyclops thomasi",
                                                           "Diaphanosoma birgei",
                                                           "Diaptomus sp."))

finaleelvlat$Taxab <- factor(finaleelvlat$Taxa, levels = c("Ectocyclops sp.",
                                                           "Epischura nevadensis",
                                                           "Eucyclops agilis",
                                                           "Eucyclops elegans",
                                                           "Eurycercus lamellatus",
                                                           "Graptoleberis testudinaria",
                                                           "Hesperodiaptomus arcticus",
                                                           "Hesperodiaptomus eiseni",
                                                           "Hesperodiaptomus franciscanus",
                                                           "Hesperodiaptomus kenai",
                                                           "Hesperodiaptomus shoshone",
                                                           "Heterocope septentrionalis",
                                                           "Holopedium gibberum",
                                                           "Kurzia latissima",
                                                           "Latona setifera",
                                                           "Leptodiaptomus angustilobus",
                                                           "Leptodiaptomus ashlandi",
                                                           "Leptodiaptomus connexus",
                                                           "Leptodiaptomus cuauhtemoci",
                                                           "Leptodiaptomus novamexicanus",
                                                           "Leptodiaptomus nudus",
                                                           "Leptodiaptomus sicilis",
                                                           "Leptodiaptomus siciloides",
                                                           "Leptodiaptomus signicauda",
                                                           "Leptodiaptomus tyrrelli",
                                                           "Leptodora kindtii",
                                                           "Leydigia leydigi",
                                                           "Macrocyclops sp.",
                                                           "Macrocyclops albidus",
                                                           "Macrocyclops fuscus",
                                                           "Macrothrix hirsuticornis",
                                                           "Macrothrix laticornis",
                                                           "Microcyclops varicans",
                                                           "Onychodiaptomus sanguineus",
                                                           "Orthocyclops modestus",
                                                           "Paracyclops poppei",
                                                           "Picripleuroxus denticulatus",
                                                           "Pleuroxus sp.",
                                                           "Polyphemus pediculus",
                                                           "Scapholeberis kingi",
                                                           "Senecella calanoides",
                                                           "Sida crystallina",
                                                           "Simocephalus serrulatus",
                                                           "Simocephalus vetulus",
                                                           "Skistodiaptomus oregonensis",
                                                           "Streblocerus serricaudatus",
                                                           "Tropocyclops prasinus"))

cols <- c("#a60505", "#c5600a", "#e2d110", "#a2ec29", "#67eb4d", "#70ec95", "#50ebbe", "#30d4ea", "#1476e5", "#0e1ac9", "#420aac", "#6e068f")

(tax.ridge.elev.a <- finaleelvlat %>%
    filter(Geo == "Elevation") %>%
    filter(Taxaa != "NA") %>%
    ggplot(aes(x = Value, y = fct_rev(Taxa), fill = ..x..)) +
    geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2, scale = 2, rel_min_height = 0.01) +
    scale_fill_binned(palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
                      breaks = c(258, 599, 856, 1100, 1320, 1500, 1670, 1870, 2090, 2330, 2780)) +
    scale_x_continuous("Elevational range", breaks = c(0, 750, 1500, 2250, 3000, 3750), limit=c(0, 3750)) +
    scale_y_discrete("Taxa") +
    coord_cartesian(clip = "on", xlim = c(0, 3750)) +
    theme_bw() +
    theme(
    panel.grid.major = element_line(size = 0.5), 
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
    panel.border = element_blank(), 
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    axis.title.x = element_text(size = "10", face = "plain", colour = "black"),
    axis.text.x = element_text(size = "9", face = "plain", colour = "black"),
    axis.title.y = element_text(size = "10", face = "plain", colour = "black"),
    axis.text.y = element_text(vjust = 0, size = "9", face = "plain", colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines"),
    strip.background = element_rect(fill = NA),
    strip.text = element_blank(),
    strip.text.x = element_text(size = 12, colour = "black")))

(tax.ridge.elev.b <- finaleelvlat %>%
    filter(Geo == "Elevation") %>%
    filter(Taxab != "NA") %>%
    ggplot(aes(x = Value, y = fct_rev(Taxa), fill = ..x..)) +
    geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2, scale = 2, rel_min_height = 0.01) +
    scale_fill_binned(palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
                      breaks = c(258, 599, 856, 1100, 1320, 1500, 1670, 1870, 2090, 2330, 2780)) +
    scale_x_continuous("Elevational range", breaks = c(0, 750, 1500, 2250, 3000, 3750), limit=c(0, 3750)) +
    scale_y_discrete("Taxa") +
    coord_cartesian(clip = "on", xlim = c(0, 3750)) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(size = 0.5), 
      panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
      panel.border = element_blank(), 
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1),"cm"),
      axis.title.x = element_text(size = "10", face = "plain", colour = "black"),
      axis.text.x = element_text(size = "9", face = "plain", colour = "black"),
      axis.title.y = element_text(size = "10", face = "plain", colour = "black"),
      axis.text.y = element_text(vjust = 0, size = "9", face = "plain", colour = "black"),
      axis.ticks = element_line(colour = "black"),
      legend.title = element_blank(),
      legend.position = "none",
      panel.spacing = unit(0.5, "lines"),
      strip.background = element_rect(fill = NA),
      strip.text = element_blank(),
      strip.text.x = element_text(size = 12, colour = "black")))

(tax.ridge.lat.a <- finaleelvlat %>%
  filter(Geo == "Latitude") %>%
  filter(Taxaa != "NA") %>%
  ggplot(aes(x = Value, y = fct_rev(Taxa), fill = ..x..)) +
  geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2, scale = 2, rel_min_height = 0.01) +
  scale_fill_binned(palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
                    breaks = c(39.3, 41.2, 43.2, 44.6, 46.3, 48.0, 50.0, 52.1, 54.2, 58.1, 62.1)) +
  scale_x_continuous("Latitudinal range", breaks = c(36, 42, 48, 54, 60, 66), limit = c(36, 66)) +
  scale_y_discrete("Taxa") +
  coord_cartesian(clip = "on", xlim = c(36, 66)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(size = 0.5), 
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
    panel.border = element_blank(), 
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    axis.title.x = element_text(size = "10", face = "plain", colour = "black"),
    axis.text.x = element_text(size = "9", face = "plain", colour = "black"),
    axis.title.y = element_text(size = "10", face = "plain", colour = "black"),
    axis.text.y = element_text(vjust = 0, size = "9", face = "plain", colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines"),
    strip.background = element_rect(fill = NA),
    strip.text = element_blank(),
    strip.text.x = element_text(size = 12, colour = "black")))

(tax.ridge.lat.b <- finaleelvlat %>%
    filter(Geo == "Latitude") %>%
    filter(Taxab!="NA") %>%
    ggplot(aes(x = Value, y = fct_rev(Taxa), fill = ..x..)) +
    geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2, scale = 2, rel_min_height = 0.01) +
    scale_fill_binned(palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
                      breaks = c(39.3, 41.2, 43.2, 44.6, 46.3, 48.0, 50.0, 52.1, 54.2, 58.1, 62.1)) +
    scale_x_continuous("Latitudinal range", breaks = c(36, 42, 48, 54, 60, 66), limit = c(36, 66)) +
    scale_y_discrete("Taxa") +
    coord_cartesian(clip = "on", xlim = c(36, 66)) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(size = 0.5), 
      panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
      panel.border = element_blank(), 
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      axis.title.x = element_text(size = "10", face = "plain", colour = "black"),
      axis.text.x = element_text(size = "9", face = "plain", colour = "black"),
      axis.title.y = element_text(size = "10", face = "plain", colour = "black"),
      axis.text.y = element_text(vjust = 0, size = "9", face = "plain", colour = "black"),
      axis.ticks = element_line(colour = "black"),
      legend.title = element_blank(),
      legend.position = "none",
      panel.spacing = unit(0.5, "lines"),
      strip.background = element_rect(fill = NA),
      strip.text = element_blank(),
      strip.text.x = element_text(size = 12, colour = "black")))

sp.distributions$constant <- c(1)

(sample.ridge.elev <- sp.distributions %>%
  ggplot(aes(x = Elevation, y = constant, fill = ..x..)) +
  geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2, scale = 2) +
  scale_fill_binned(palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
                    breaks = c(258, 599, 856, 1100, 1320, 1500, 1670, 1870, 2090, 2330, 2780)) +
  scale_x_continuous("Elevational range", breaks=c(0, 750, 1500, 2250, 3000, 3750), limit=c(0, 3750)) +
  coord_cartesian(clip = "on", xlim = c(0, 3750)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(size = 0.5), 
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
    panel.border = element_blank(), 
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    axis.title.x = element_text(size = "10", face = "plain", colour = "black"),
    axis.text.x = element_text(size = "9", face = "plain", colour = "black"),
    axis.title.y = element_text(size = "10", face = "plain", colour = "black"),
    axis.text.y = element_text(vjust = 0, size = "9", face = "plain", colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines"),
    strip.background = element_rect(fill = NA),
    strip.text = element_blank(),
    strip.text.x = element_text(size = 12, colour = "black")))

(sample.ridge.lat <- sp.distributions %>%
  ggplot(aes(x = Latitude, y = constant, fill = ..x..)) +
  geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2, scale = 2) +
  scale_fill_binned(palette = ggplot2:::binned_pal(scales::manual_pal(values = cols)),
                    breaks = c(39.3, 41.2, 43.2, 44.6, 46.3, 48.0, 50.0, 52.1, 54.2, 58.1, 62.1)) +
  scale_x_continuous("Latitudinal range", breaks = c(36, 42, 48, 54, 60, 66), limit = c(36, 66)) +
  coord_cartesian(clip = "on", xlim = c(36,66)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(size = 0.5), 
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
    panel.border = element_blank(), 
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    axis.title.x = element_text(size = "10", face = "plain", colour = "black"),
    axis.text.x = element_text(size = "9", face = "plain", colour = "black"),
    axis.title.y = element_text(size = "10", face = "plain", colour = "black"),
    axis.text.y = element_text(vjust = 0, size = "9", face = "plain", colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines"),
    strip.background = element_rect(fill = NA),
    strip.text = element_blank(),
    strip.text.x = element_text(size = 12, colour = "black")))

                                            
##################
# Metrics
##################

# Prepare species data for calculating biodiversity metrics
Bio.metrics <- data.frame(sites$ID)
Bio.metrics <- rename(Bio.metrics, c("ID" = "sites.ID"))
rownames(Bio.metrics) <- Bio.metrics$ID

Bio.metrics <- merge(Bio.metrics, clean.env.lat["Lake"], by = 0)
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

Bio.metrics <- merge(Bio.metrics, clean.env.lat["Latitude"], by = 0)
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

Bio.metrics <- merge(Bio.metrics, clean.env.elev["Elevation"], by = 0)
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

Bio.metrics <- merge(Bio.metrics, clean.env.lat["lat_class"], by = 0)
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

Bio.metrics <- merge(Bio.metrics, clean.env.elev["elev_class"], by = 0)
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

# Calculate taxonomic richness
taxo.alpha.sum = as.data.frame(rowSums(taxo.alpha))
Bio.metrics <- merge(Bio.metrics, taxo.alpha.sum, by = 0)
Bio.metrics <- rename(Bio.metrics, c("Taxo.alpha" = "rowSums(taxo.alpha)"))
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

# Calculate functional richness and dispersion
t <- c("Length", "Feed.ord")

sel.trait <- traits[, t] #numeric feeding groups
sel.trait <- as.matrix(sel.trait)
rownames(sel.trait) <- traits$Sp

sel.trait.ord <- as.data.frame(sel.trait)
sel.trait.ord$Feed.ord <- ordered(sel.trait.ord$Feed.ord) #define feeding groups as ordered variable

func.alpha.factor <- dbFD(sel.trait.ord, taxo.beta, ord = "metric", corr = "sqrt", stand.FRic = TRUE) #calculate metrics

Bio.metrics <- merge(Bio.metrics, func.alpha.factor["nbsp"], by = 0) #collect number of species used to claculate each metric
Bio.metrics <- rename(Bio.metrics, c("Glob.nbsp" = "nbsp"))
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

Bio.metrics <- merge(Bio.metrics, func.alpha.factor["FRic"], by = 0) #collect raw/global functional richness (standardized between 0 and 1)
Bio.metrics <- rename(Bio.metrics, c("Glob.FRic" = "FRic"))
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

Bio.metrics <- merge(Bio.metrics, func.alpha.factor["FDis"], by = 0) #collected raw/global functional dispersion
Bio.metrics <- rename(Bio.metrics, c("Glob.FDis" = "FDis"))
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

Bio.metrics$Glob.FDis[Bio.metrics$Glob.FDis == 0] <- NA #replace values of zero with NA for functional dispersion

# Calculate functional CWM body length
func.alpha.cwm <- functcomp(sel.trait, as.matrix(taxo.beta))

Bio.metrics <- merge(Bio.metrics, func.alpha.cwm["Length"], by = 0) #collect community weighted mean body lengths
Bio.metrics <- rename(Bio.metrics, c("CWM.Bl" = "Length"))
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

# Calculate phylogenetic metrics
taxa.list.tree <- class2tree(taxa.list, varstep = TRUE, check = TRUE)
tree.distmat.matrix <- as.matrix(taxa.list.tree$distmat)

# Plot phylogenetic tree as fan plot
# requires some minor revisions (e.g. adjusting colours) in a graphics editor
plot(taxa.list.tree, type = "fan", no.margin = TRUE)

# Calculate phylogenetic richness
phylo.alpha.PD <- pd(taxo.beta,taxa.list.tree$phylo,include.root = FALSE)

Bio.metrics <- merge(Bio.metrics, phylo.alpha.PD["PD"], by = 0) #collect phylogenetic richness
Bio.metrics <- rename(Bio.metrics, c("Glob.PD" = "PD"))
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

# Calculate phylogenetic mean pairwise distance
phylo.alpha.MPD <- mpd(taxo.beta,tree.distmat.matrix)
names(phylo.alpha.MPD) <- rownames(taxo.beta)
phylo.alpha.MPD <- as.data.frame(phylo.alpha.MPD)

Bio.metrics <- merge(Bio.metrics, phylo.alpha.MPD, by = 0) #collected raw/global phylogenetic mean pairwise distances
Bio.metrics <- rename(Bio.metrics, c("Glob.MPD" = "phylo.alpha.MPD"))
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

Bio.metrics$Glob.MPD.beta <- Bio.metrics$Glob.MPD / 100 #bound data between zero and one

# Calculate standardized phylogenetic richness (SES)

# Iterate function defining null communities by latitudinal zones
set.seed(1) #set seed for reproducibility
for(i in unique(clean.env.lat$lat_class)) {
  
  clean.env.lat.class = subset(clean.env.lat, lat_class == i)
  sp.class <- taxo.beta[(rownames(taxo.beta) %in% rownames(clean.env.lat.class)), ]
  sp.class <- sp.class[, colSums(sp.class) != 0]
  
  pruned.tree.phylo <- prune.sample(sp.class, taxa.list.tree$phylo)
  
  newname <- i
  assign(newname, sp.class)
  
  PD.ses.lat = ses.pd(sp.class, pruned.tree.phylo, include.root = FALSE, null.model = "independentswap", 
                      runs = 1000, iterations = 10000)

  newname.phylo.alpha.PD.ses.lat <- paste("ses.PD.lat", i, sep = "")
  assign(newname.phylo.alpha.PD.ses.lat, PD.ses.lat)

}  

# Iterate function defining null communities by elevational zones
set.seed(1) #set seed for reproducibility
for(i in unique(clean.env.elev$elev_class)) {
  
  clean.env.elev.class = subset(clean.env.elev, elev_class == i)
  sp.class <- taxo.beta[(rownames(taxo.beta) %in% rownames(clean.env.elev.class)), ]
  sp.class <- sp.class[, colSums(sp.class) != 0]
  
  pruned.tree.phylo <- prune.sample(sp.class, taxa.list.tree$phylo)

  newname <- i
  assign(newname, sp.class)
  
  PD.ses.ele = ses.pd(sp.class, pruned.tree.phylo, include.root = FALSE, null.model = "independentswap", 
                      runs = 1000, iterations = 10000)

  newname.phylo.alpha.PD.ses.ele <- paste("ses.PD.ele", i, sep = "")
  assign(newname.phylo.alpha.PD.ses.ele, PD.ses.ele)

}  

# Merge and collect results
PD.ses.lat.comb <- do.call("rbind", list(`ses.PD.lat36.6-39.3`,
                                         `ses.PD.lat39.3-41.2`,
                                         `ses.PD.lat41.2-43.2`,
                                         `ses.PD.lat43.2-44.6`,
                                         `ses.PD.lat44.6-46.3`,
                                         `ses.PD.lat46.3-48.0`,
                                         `ses.PD.lat48.0-50.0`,
                                         `ses.PD.lat50.0-52.1`,
                                         `ses.PD.lat52.1-54.2`,
                                         `ses.PD.lat54.2-58.1`,
                                         `ses.PD.lat58.1-62.1`,
                                         `ses.PD.lat62.1-66.2`))

PD.ses.lat.comb <- rename(PD.ses.lat.comb, c("pd.ntaxa.lat" = "ntaxa"))
PD.ses.lat.comb <- rename(PD.ses.lat.comb, c("pd.obs.lat" = "pd.obs"))
PD.ses.lat.comb <- rename(PD.ses.lat.comb, c("pd.rand.mean.lat" = "pd.rand.mean"))
PD.ses.lat.comb <- rename(PD.ses.lat.comb, c("pd.rand.sd.lat" = "pd.rand.sd"))
PD.ses.lat.comb <- rename(PD.ses.lat.comb, c("pd.obs.rank.lat" = "pd.obs.rank"))
PD.ses.lat.comb <- rename(PD.ses.lat.comb, c("pd.obs.z.lat" = "pd.obs.z"))
PD.ses.lat.comb <- rename(PD.ses.lat.comb, c("pd.obs.p.lat" = "pd.obs.p"))
PD.ses.lat.comb <- rename(PD.ses.lat.comb, c("pd.runs.lat" = "runs"))

Bio.metrics <- merge(Bio.metrics, PD.ses.lat.comb, by = "row.names")
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

PD.ses.ele.comb <- do.call("rbind", list(`ses.PD.ele0-258`,
                                         `ses.PD.ele258-599`,
                                         `ses.PD.ele599-856`,
                                         `ses.PD.ele856-1100`,
                                         `ses.PD.ele1100-1320`,
                                         `ses.PD.ele1320-1500`,
                                         `ses.PD.ele1500-1670`,
                                         `ses.PD.ele1670-1870`,
                                         `ses.PD.ele1870-2090`,
                                         `ses.PD.ele2090-2330`,
                                         `ses.PD.ele2330-2780`,
                                         `ses.PD.ele2780-3740`))

PD.ses.ele.comb <- rename(PD.ses.ele.comb, c("pd.ntaxa.ele" = "ntaxa"))
PD.ses.ele.comb <- rename(PD.ses.ele.comb, c("pd.obs.ele" = "pd.obs"))
PD.ses.ele.comb <- rename(PD.ses.ele.comb, c("pd.rand.mean.ele" = "pd.rand.mean"))
PD.ses.ele.comb <- rename(PD.ses.ele.comb, c("pd.rand.sd.ele" = "pd.rand.sd"))
PD.ses.ele.comb <- rename(PD.ses.ele.comb, c("pd.obs.rank.ele" = "pd.obs.rank"))
PD.ses.ele.comb <- rename(PD.ses.ele.comb, c("pd.obs.z.ele" = "pd.obs.z"))
PD.ses.ele.comb <- rename(PD.ses.ele.comb, c("pd.obs.p.ele" = "pd.obs.p"))
PD.ses.ele.comb <- rename(PD.ses.ele.comb, c("pd.runs.ele" = "runs"))

Bio.metrics <- merge(Bio.metrics, PD.ses.ele.comb, by = "row.names")
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

# Calculate standardized mean pairwise distance (SES)

# Iterate function defining null communities by latitudinal zones
set.seed(1) #set seed for reproducibility
for(i in unique(clean.env.lat$lat_class)) {
  
  clean.env.lat.class = subset(clean.env.lat, lat_class == i)
  sp.class <- taxo.beta[(rownames(taxo.beta) %in% rownames(clean.env.lat.class)), ]
  sp.class <- sp.class[, colSums(sp.class) != 0]
  
  taxa.list.tree.distmat <- tree.distmat.matrix[(rownames(tree.distmat.matrix) %in% colnames(sp.class)), ]
  taxa.list.tree.distmat <- taxa.list.tree.distmat[, (colnames(taxa.list.tree.distmat) %in% colnames(sp.class))]
  
  newname <- i
  assign(newname, sp.class)
  
  MPD.ses.lat = ses.mpd(sp.class, taxa.list.tree.distmat, null.model = "independentswap", 
                      runs = 1000, iterations = 10000)
  
  newname.phylo.alpha.MPD.ses.lat <- paste("ses.MPD.lat", i, sep = "")
  assign(newname.phylo.alpha.MPD.ses.lat, MPD.ses.lat)
  
}  

# Iterate function defining null communities by elevational zones
set.seed(1) #set seed for reproducibility
for(i in unique(clean.env.elev$elev_class)) {
  
  clean.env.elev.class = subset(clean.env.elev, elev_class == i)
  sp.class <- taxo.beta[(rownames(taxo.beta) %in% rownames(clean.env.elev.class)), ]
  sp.class <- sp.class[, colSums(sp.class) != 0]
  
  taxa.list.tree.distmat <- tree.distmat.matrix[(rownames(tree.distmat.matrix) %in% colnames(sp.class)), ]
  taxa.list.tree.distmat <- taxa.list.tree.distmat[, (colnames(taxa.list.tree.distmat) %in% colnames(sp.class))]
  
  newname <- i
  assign(newname, sp.class)
  
  MPD.ses.ele = ses.mpd(sp.class, taxa.list.tree.distmat, null.model = "independentswap", 
                      runs = 1000, iterations = 10000)
  
  newname.phylo.alpha.MPD.ses.ele <- paste("ses.MPD.ele", i, sep = "")
  assign(newname.phylo.alpha.MPD.ses.ele, MPD.ses.ele)
  
}

# Merge and collect results
MPD.ses.lat.comb <- do.call("rbind", list(`ses.MPD.lat36.6-39.3`,
                                          `ses.MPD.lat39.3-41.2`,
                                          `ses.MPD.lat41.2-43.2`,
                                          `ses.MPD.lat43.2-44.6`,
                                          `ses.MPD.lat44.6-46.3`,
                                          `ses.MPD.lat46.3-48.0`,
                                          `ses.MPD.lat48.0-50.0`,
                                          `ses.MPD.lat50.0-52.1`,
                                          `ses.MPD.lat52.1-54.2`,
                                          `ses.MPD.lat54.2-58.1`,
                                          `ses.MPD.lat58.1-62.1`,
                                          `ses.MPD.lat62.1-66.2`))

MPD.ses.lat.comb <- rename(MPD.ses.lat.comb, c("mpd.ntaxa.lat" = "ntaxa"))
MPD.ses.lat.comb <- rename(MPD.ses.lat.comb, c("mpd.obs.lat" = "mpd.obs"))
MPD.ses.lat.comb <- rename(MPD.ses.lat.comb, c("mpd.rand.mean.lat" = "mpd.rand.mean"))
MPD.ses.lat.comb <- rename(MPD.ses.lat.comb, c("mpd.rand.sd.lat" = "mpd.rand.sd"))
MPD.ses.lat.comb <- rename(MPD.ses.lat.comb, c("mpd.obs.rank.lat" = "mpd.obs.rank"))
MPD.ses.lat.comb <- rename(MPD.ses.lat.comb, c("mpd.obs.z.lat" = "mpd.obs.z"))
MPD.ses.lat.comb <- rename(MPD.ses.lat.comb, c("mpd.obs.p.lat" = "mpd.obs.p"))
MPD.ses.lat.comb <- rename(MPD.ses.lat.comb, c("mpd.runs.lat" = "runs"))

Bio.metrics <- merge(Bio.metrics, MPD.ses.lat.comb, by = "row.names")
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

MPD.ses.ele.comb <- do.call("rbind", list(`ses.MPD.ele0-258`,
                                          `ses.MPD.ele258-599`,
                                          `ses.MPD.ele599-856`,
                                          `ses.MPD.ele856-1100`,
                                          `ses.MPD.ele1100-1320`,
                                          `ses.MPD.ele1320-1500`,
                                          `ses.MPD.ele1500-1670`,
                                          `ses.MPD.ele1670-1870`,
                                          `ses.MPD.ele1870-2090`,
                                          `ses.MPD.ele2090-2330`,
                                          `ses.MPD.ele2330-2780`,
                                          `ses.MPD.ele2780-3740`))

MPD.ses.ele.comb <- rename(MPD.ses.ele.comb, c("mpd.ntaxa.ele" = "ntaxa"))
MPD.ses.ele.comb <- rename(MPD.ses.ele.comb, c("mpd.obs.ele" = "mpd.obs"))
MPD.ses.ele.comb <- rename(MPD.ses.ele.comb, c("mpd.rand.mean.ele" = "mpd.rand.mean"))
MPD.ses.ele.comb <- rename(MPD.ses.ele.comb, c("mpd.rand.sd.ele" = "mpd.rand.sd"))
MPD.ses.ele.comb <- rename(MPD.ses.ele.comb, c("mpd.obs.rank.ele" = "mpd.obs.rank"))
MPD.ses.ele.comb <- rename(MPD.ses.ele.comb, c("mpd.obs.z.ele" = "mpd.obs.z"))
MPD.ses.ele.comb <- rename(MPD.ses.ele.comb, c("mpd.obs.p.ele" = "mpd.obs.p"))
MPD.ses.ele.comb <- rename(MPD.ses.ele.comb, c("mpd.runs.ele" = "runs"))

Bio.metrics <- merge(Bio.metrics, MPD.ses.ele.comb, by = "row.names")
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

# Calculate standardized functional richness (SES)

# Create function for estimating standardized effect sizes
ses.FRic <-
  function (traits, species, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"), 
            abundance.weighted = FALSE, runs = 1000, iterations = 10000)
  {
    FRic.obs <- dbFD(traits, species, stand.FRic = FALSE, ord = "metric", corr = "sqrt")$FRic
    null.model <- match.arg(null.model)
    FRic.rand <- switch(null.model,
                        independentswap = t(replicate(runs, dbFD(traits, randomizeMatrix(species, null.model = "independentswap", iterations), stand.FRic = FALSE, ord = "metric", corr = "sqrt")$FRic)), )
    FRic.rand.mean <- apply(X = FRic.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
    FRic.rand.sd <- apply(X = FRic.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    FRic.obs.z <- (FRic.obs - FRic.rand.mean) / FRic.rand.sd
    FRic.obs.rank <- apply(X = rbind(FRic.obs, FRic.rand), MARGIN = 2,
                           FUN = rank)[1, ]
    FRic.obs.rank <- ifelse(is.na(FRic.rand.mean), NA, FRic.obs.rank)
    data.frame(ntaxa = specnumber(species),FRic.obs, FRic.rand.mean, FRic.rand.sd, FRic.obs.rank,
               FRic.obs.z, FRic.obs.p = FRic.obs.rank / (runs + 1),runs = runs, row.names = row.names(species))
  }

# Iterate function defining null communities by latitudinal zones
set.seed(1) #set seed for reproducibility
for(i in unique(clean.env.lat$lat_class)) {
  
  clean.env.lat.class = subset(clean.env.lat, lat_class == i)
  sp.class <- taxo.beta[(rownames(taxo.beta) %in% rownames(clean.env.lat.class)), ]
  sp.class <- sp.class[, colSums(sp.class) != 0]
  
  trait.class.factor <- sel.trait.ord[(rownames(sel.trait.ord) %in% colnames(sp.class)), ]

  newname <- i
  assign(newname, sp.class)
  
  FRic.ses.lat = ses.FRic(trait.class.factor, sp.class, null.model = "independentswap",
                          abundance.weighted = FALSE, runs = 1000, iterations = 10000)
  
  newname.dbFD.FRic.ses.lat <- paste("ses.FRic.lat", i, sep = "")
  assign(newname.dbFD.FRic.ses.lat, FRic.ses.lat)
}

# Iterate function defining null communities by elevational zones
set.seed(1) #set seed for reproducibility
for(i in unique(clean.env.elev$elev_class)) {
  
  clean.env.elev.class = subset(clean.env.elev, elev_class == i)
  sp.class <- taxo.beta[(rownames(taxo.beta) %in% rownames(clean.env.elev.class)), ]
  sp.class <- sp.class[, colSums(sp.class) != 0]
  
  trait.class.factor <- sel.trait.ord[(rownames(sel.trait.ord) %in% colnames(sp.class)), ]

  newname <- i
  assign(newname, sp.class)
  
  FRic.ses.ele = ses.FRic(trait.class.factor, sp.class, null.model = "independentswap",
                          abundance.weighted = FALSE, runs = 1000, iterations = 10000)
  
  newname.dbFD.FRic.ses.ele <- paste("ses.FRic.ele", i, sep = "")
  assign(newname.dbFD.FRic.ses.ele, FRic.ses.ele)
}

# Merge and collect results
FRic.ses.lat.comb <- do.call("rbind", list(`ses.FRic.lat36.6-39.3`,
                                           `ses.FRic.lat39.3-41.2`,
                                           `ses.FRic.lat41.2-43.2`,
                                           `ses.FRic.lat43.2-44.6`,
                                           `ses.FRic.lat44.6-46.3`,
                                           `ses.FRic.lat46.3-48.0`,
                                           `ses.FRic.lat48.0-50.0`,
                                           `ses.FRic.lat50.0-52.1`,
                                           `ses.FRic.lat52.1-54.2`,
                                           `ses.FRic.lat54.2-58.1`,
                                           `ses.FRic.lat58.1-62.1`,
                                           `ses.FRic.lat62.1-66.2`))

FRic.ses.lat.comb <- rename(FRic.ses.lat.comb, c("FRic.ntaxa.lat" = "ntaxa"))
FRic.ses.lat.comb <- rename(FRic.ses.lat.comb, c("FRic.obs.lat" = "FRic.obs"))
FRic.ses.lat.comb <- rename(FRic.ses.lat.comb, c("FRic.rand.mean.lat" = "FRic.rand.mean"))
FRic.ses.lat.comb <- rename(FRic.ses.lat.comb, c("FRic.rand.sd.lat" = "FRic.rand.sd"))
FRic.ses.lat.comb <- rename(FRic.ses.lat.comb, c("FRic.obs.rank.lat" = "FRic.obs.rank"))
FRic.ses.lat.comb <- rename(FRic.ses.lat.comb, c("FRic.obs.z.lat" = "FRic.obs.z"))
FRic.ses.lat.comb <- rename(FRic.ses.lat.comb, c("FRic.obs.p.lat" = "FRic.obs.p"))
FRic.ses.lat.comb <- rename(FRic.ses.lat.comb, c("FRic.runs.lat" = "runs"))

Bio.metrics <- merge(Bio.metrics, FRic.ses.lat.comb, by = "row.names")
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

FRic.ses.ele.comb <- do.call("rbind", list(`ses.FRic.ele0-258`,
                                           `ses.FRic.ele258-599`,
                                           `ses.FRic.ele599-856`,
                                           `ses.FRic.ele856-1100`,
                                           `ses.FRic.ele1100-1320`,
                                           `ses.FRic.ele1320-1500`,
                                           `ses.FRic.ele1500-1670`,
                                           `ses.FRic.ele1670-1870`,
                                           `ses.FRic.ele1870-2090`,
                                           `ses.FRic.ele2090-2330`,
                                           `ses.FRic.ele2330-2780`,
                                           `ses.FRic.ele2780-3740`))

FRic.ses.ele.comb <- rename(FRic.ses.ele.comb, c("FRic.ntaxa.ele" = "ntaxa"))
FRic.ses.ele.comb <- rename(FRic.ses.ele.comb, c("FRic.obs.ele" = "FRic.obs"))
FRic.ses.ele.comb <- rename(FRic.ses.ele.comb, c("FRic.rand.mean.ele" = "FRic.rand.mean"))
FRic.ses.ele.comb <- rename(FRic.ses.ele.comb, c("FRic.rand.sd.ele" = "FRic.rand.sd"))
FRic.ses.ele.comb <- rename(FRic.ses.ele.comb, c("FRic.obs.rank.ele" = "FRic.obs.rank"))
FRic.ses.ele.comb <- rename(FRic.ses.ele.comb, c("FRic.obs.z.ele" = "FRic.obs.z"))
FRic.ses.ele.comb <- rename(FRic.ses.ele.comb, c("FRic.obs.p.ele" = "FRic.obs.p"))
FRic.ses.ele.comb <- rename(FRic.ses.ele.comb, c("FRic.runs.ele" = "runs"))

Bio.metrics <- merge(Bio.metrics, FRic.ses.ele.comb, by = "row.names")
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

# Calculate standardized functional dispersion (SES)

# Create function for estimating standardized effect sizes
ses.FDis <-
  function (traits, species, null.model = c("taxa.labels", "richness", "frequency", "sample.pool",
                                          "phylogeny.pool", "independentswap", "trialswap"),
            abundance.weighted = FALSE, runs = 1000, iterations = 10000)
  {
    FDis.obs <- fdisp(gowdis(traits, ord = "metric"), as.matrix(species))$FDis
    null.model <- match.arg(null.model)
    FDis.rand <- switch(null.model,
                        independentswap = t(replicate(runs, fdisp(gowdis(traits, ord = "metric"), randomizeMatrix(species, null.model = "independentswap", iterations))$FDis)), )
    FDis.rand.mean <- apply(X = FDis.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
    FDis.rand.sd <- apply(X = FDis.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    FDis.obs.z <- (FDis.obs - FDis.rand.mean) / FDis.rand.sd
    FDis.obs.rank <- apply(X = rbind(FDis.obs, FDis.rand), MARGIN = 2,
                           FUN = rank)[1, ]
    FDis.obs.rank <- ifelse(is.na(FDis.rand.mean), NA, FDis.obs.rank)
    data.frame(ntaxa = specnumber(species), FDis.obs, FDis.rand.mean, FDis.rand.sd, FDis.obs.rank,
               FDis.obs.z, FDis.obs.p = FDis.obs.rank / (runs + 1),runs = runs, row.names = row.names(species))
  }

# Iterate function defining null communities by latitudinal zones
set.seed(1) #set seed for reproducibility
for(i in unique(clean.env.lat$lat_class)) {
  
  clean.env.lat.class=subset(clean.env.lat, lat_class == i)
  sp.class <- taxo.beta[(rownames(taxo.beta) %in% rownames(clean.env.lat.class)), ]
  sp.class <- sp.class[, colSums(sp.class) != 0]
  
  trait.class.factor <- sel.trait.ord[(rownames(sel.trait.ord) %in% colnames(sp.class)), ]
  
  newname <- i
  assign(newname, sp.class)
  
  FDis.ses.lat = ses.FDis(trait.class.factor, sp.class, null.model = "independentswap",
                          abundance.weighted = FALSE, runs = 1000, iterations = 10000)
  
  newname.dbFD.FDis.ses.lat <- paste("ses.FDis.lat", i, sep = "")
  assign(newname.dbFD.FDis.ses.lat, FDis.ses.lat)
}

# Iterate function defining null communities by elevational zones
set.seed(1) #set seed for reproducibility
for(i in unique(clean.env.elev$elev_class)) {
  
  clean.env.elev.class = subset(clean.env.elev, elev_class == i)
  sp.class <- taxo.beta[(rownames(taxo.beta) %in% rownames(clean.env.elev.class)), ]
  sp.class <- sp.class[, colSums(sp.class) != 0]
  
  trait.class.factor <- sel.trait.ord[(rownames(sel.trait.ord) %in% colnames(sp.class)), ]
  
  newname <- i
  assign(newname, sp.class)
  
  FDis.ses.ele = ses.FDis(trait.class.factor, sp.class, null.model = "independentswap",
                          abundance.weighted = FALSE, runs = 1000, iterations = 10000)
  
  newname.dbFD.FDis.ses.ele <- paste("ses.FDis.ele", i, sep = "")
  assign(newname.dbFD.FDis.ses.ele, FDis.ses.ele)
}

# Merge and collect results
FDis.ses.lat.comb <- do.call("rbind", list(`ses.FDis.lat36.6-39.3`,
                                           `ses.FDis.lat39.3-41.2`,
                                           `ses.FDis.lat41.2-43.2`,
                                           `ses.FDis.lat43.2-44.6`,
                                           `ses.FDis.lat44.6-46.3`,
                                           `ses.FDis.lat46.3-48.0`,
                                           `ses.FDis.lat48.0-50.0`,
                                           `ses.FDis.lat50.0-52.1`,
                                           `ses.FDis.lat52.1-54.2`,
                                           `ses.FDis.lat54.2-58.1`,
                                           `ses.FDis.lat58.1-62.1`,
                                           `ses.FDis.lat62.1-66.2`))

FDis.ses.lat.comb <- rename(FDis.ses.lat.comb, c("FDis.ntaxa.lat" = "ntaxa"))
FDis.ses.lat.comb <- rename(FDis.ses.lat.comb, c("FDis.obs.lat" = "FDis.obs"))
FDis.ses.lat.comb <- rename(FDis.ses.lat.comb, c("FDis.rand.mean.lat" = "FDis.rand.mean"))
FDis.ses.lat.comb <- rename(FDis.ses.lat.comb, c("FDis.rand.sd.lat" = "FDis.rand.sd"))
FDis.ses.lat.comb <- rename(FDis.ses.lat.comb, c("FDis.obs.rank.lat" = "FDis.obs.rank"))
FDis.ses.lat.comb <- rename(FDis.ses.lat.comb, c("FDis.obs.z.lat" = "FDis.obs.z"))
FDis.ses.lat.comb <- rename(FDis.ses.lat.comb, c("FDis.obs.p.lat" = "FDis.obs.p"))
FDis.ses.lat.comb <- rename(FDis.ses.lat.comb, c("FDis.runs.lat" = "runs"))

Bio.metrics <- merge(Bio.metrics, FDis.ses.lat.comb, by = "row.names")
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

FDis.ses.ele.comb <- do.call("rbind", list(`ses.FDis.ele0-258`,
                                           `ses.FDis.ele258-599`,
                                           `ses.FDis.ele599-856`,
                                           `ses.FDis.ele856-1100`,
                                           `ses.FDis.ele1100-1320`,
                                           `ses.FDis.ele1320-1500`,
                                           `ses.FDis.ele1500-1670`,
                                           `ses.FDis.ele1670-1870`,
                                           `ses.FDis.ele1870-2090`,
                                           `ses.FDis.ele2090-2330`,
                                           `ses.FDis.ele2330-2780`,
                                           `ses.FDis.ele2780-3740`))

FDis.ses.ele.comb <- rename(FDis.ses.ele.comb, c("FDis.ntaxa.ele" = "ntaxa"))
FDis.ses.ele.comb <- rename(FDis.ses.ele.comb, c("FDis.obs.ele" = "FDis.obs"))
FDis.ses.ele.comb <- rename(FDis.ses.ele.comb, c("FDis.rand.mean.ele" = "FDis.rand.mean"))
FDis.ses.ele.comb <- rename(FDis.ses.ele.comb, c("FDis.rand.sd.ele" = "FDis.rand.sd"))
FDis.ses.ele.comb <- rename(FDis.ses.ele.comb, c("FDis.obs.rank.ele" = "FDis.obs.rank"))
FDis.ses.ele.comb <- rename(FDis.ses.ele.comb, c("FDis.obs.z.ele" = "FDis.obs.z"))
FDis.ses.ele.comb <- rename(FDis.ses.ele.comb, c("FDis.obs.p.ele" = "FDis.obs.p"))
FDis.ses.ele.comb <- rename(FDis.ses.ele.comb, c("FDis.runs.ele" = "runs"))

Bio.metrics <- merge(Bio.metrics, FDis.ses.ele.comb, by = "row.names")
Bio.metrics <- Bio.metrics[-1]
rownames(Bio.metrics) <- Bio.metrics$ID

# Write biodiversity metrics to a csv file for later use
write.table(Bio.metrics, file = "Bio.metrics.pub.csv", sep = ",", row.names = TRUE, col.names = TRUE)

                                            
##################
# Bayesian generalized linear multilevel models
##################

# Load previously calculated biodiversity metrics (as compiled in previous scripts)
Bio.metrics = read.table("Bio.metrics.pub.csv", header = T, sep = ",")
rownames(Bio.metrics) <- Bio.metrics$ID

# Load climate data output from ClimateNA (Wang et al. 2016)
climate = read.table("ClimateNA_input_elev_1964-2015Y.csv", header = T, sep = ",")

# Organize climate data
Climate.MAT <- dcast(climate, ID1 ~ Year, value.var = "MAT")
rownames(Climate.MAT) <- Climate.MAT$ID1
Climate.MAT$x1964_2015 <- as.numeric(apply(Climate.MAT[, 2:53], 1, mean))
Climate.MAT <- Climate.MAT[(rownames(Climate.MAT) %in% rownames(Bio.metrics)), ]
Climate.MAT.x1964_2015 <- as.data.frame(Climate.MAT$x1964_2015, row.names = row.names(Climate.MAT))
Climate.MAT.x1964_2015 <- rename(Climate.MAT.x1964_2015, c("MAT" = "Climate.MAT$x1964_2015"))
Bio.metrics <- merge(Bio.metrics, Climate.MAT.x1964_2015, by = "row.names")
rownames(Bio.metrics) <- Bio.metrics$Row.names
Bio.metrics <- Bio.metrics[, -1]

Climate.TD <- dcast(climate, ID1 ~ Year, value.var = "TD")
rownames(Climate.TD) <- Climate.TD$ID1
Climate.TD$x1964_2015 <- as.numeric(apply(Climate.TD[, 2:53], 1, mean))
Climate.TD <- Climate.TD[(rownames(Climate.TD) %in% rownames(Bio.metrics)), ]
Climate.TD.x1964_2015 <- as.data.frame(Climate.TD$x1964_2015, row.names = row.names(Climate.TD))
Climate.TD.x1964_2015 <- rename(Climate.TD.x1964_2015, c("TD" = "Climate.TD$x1964_2015"))
Bio.metrics <- merge(Bio.metrics, Climate.TD.x1964_2015, by = "row.names")
rownames(Bio.metrics) <- Bio.metrics$Row.names
Bio.metrics <- Bio.metrics[, -1]

# Calculate z-scores for elevation, latitude, and climate variables
Bio.metrics$Elevation.scale <- scale(Bio.metrics$Elevation, center = TRUE, scale = TRUE)
Bio.metrics$Latitude.scale <- scale(Bio.metrics$Latitude, center = TRUE, scale = TRUE)
Bio.metrics$MAT.scale <- scale(Bio.metrics$MAT, center = TRUE, scale = TRUE)
Bio.metrics$TD.scale <- scale(Bio.metrics$TD, center = TRUE, scale = TRUE)

# Set order of latitudinal and elevational zones (useful for plotting later on)
Bio.metrics$lat_class <- factor(Bio.metrics$lat_class, levels = c("36.6-39.3",
                                                                  "39.3-41.2",
                                                                  "41.2-43.2",
                                                                  "43.2-44.6",
                                                                  "44.6-46.3",
                                                                  "46.3-48.0",
                                                                  "48.0-50.0",
                                                                  "50.0-52.1",
                                                                  "52.1-54.2",
                                                                  "54.2-58.1",
                                                                  "58.1-62.1",
                                                                  "62.1-66.2"), ordered = TRUE)

Bio.metrics$elev_class <- factor(Bio.metrics$elev_class, levels = c("0-258", 
                                                                    "258-599", 
                                                                    "599-856", 
                                                                    "856-1100", 
                                                                    "1100-1320", 
                                                                    "1320-1500", 
                                                                    "1500-1670",
                                                                    "1670-1870",
                                                                    "1870-2090",
                                                                    "2090-2330",
                                                                    "2330-2780",
                                                                    "2780-3740"), ordered = TRUE)
# Collect variables for modelling
lmdata <- Bio.metrics %>%
  melt(id.vars = c("ID", "Lake", "Elevation", "Elevation.scale", "Latitude", "Latitude.scale", "lat_class", "elev_class",
                   "MAT", "MAT.scale", "TD", "TD.scale"),
       measure.vars = c("Taxo.alpha",
                        "Glob.PD",
                        "Glob.MPD",
                        "Glob.MPD.beta",
                        "Glob.FRic",
                        "Glob.FDis",
                        "CWM.Bl",
                        "pd.obs.z.lat",
                        "pd.obs.z.ele",
                        "mpd.obs.z.lat",
                        "mpd.obs.z.ele",
                        "FRic.obs.z.lat",
                        "FRic.obs.z.ele",                        
                        "FDis.obs.z.lat",
                        "FDis.obs.z.ele"))

# Define theme for plotting
workingtheme <- theme(strip.background = element_blank(),
                      strip.text.y = element_text(colour = "black", size = 9),
                      strip.text.x = element_blank(),
                    
                      panel.border = element_rect(colour = "black", fill = NA),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", fill = "white", size = 0.75),
                      panel.spacing = unit(0.12, "cm"),
                      
                      axis.text.x = element_text(colour = "black", size = 8),
                      axis.text.y = element_text(colour = "black", size = 8),
                      axis.ticks = element_line(colour = "black", size = 0.4),
                      axis.title.y = element_text(colour = "black", size = 10),
                      axis.title.x = element_text(colour = "black", size = 10),
                      
                      legend.key = element_rect("white"),
                      legend.key.size = unit(0.4, "cm"),
                      legend.justification = c("right", "top"),
                      legend.position = c(0.99, 0.985),
                      legend.direction = "vertical",
                      legend.text = element_text(size = 10, colour = "black"),
                      legend.title = element_text(size = 12, colour = "black"),
                      legend.box.background = element_rect(colour = "black", size = 0.8),
                      legend.spacing.y = unit(0, 'cm'),
                      legend.margin = margin(1, 2, 0, 1))

# Set colour palate
a60505_6e068f<-c("#a60505", "#c5600a", "#e2d110", "#a2ec29", "#67eb4d", "#70ec95", "#50ebbe", "#30d4ea", "#1476e5", "#0e1ac9", "#420aac", "#6e068f")

# Set number of cores that can be run in parallel
options(mc.cores = parallel::detectCores())


## Model species richness

# Check distribution
dist <- subset(lmdata, variable == "Taxo.alpha")
hist(dist$value)

# SR ~ simple elevation fixed effect model
bayes.alpha.elev.1fe <- brm(formula = value ~ Elevation.scale,
                            data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                            prior = c(set_prior("normal(0,5)", class = "b")),
                            warmup = 1000, iter = 3000, chains = 4,
                            control = list(adapt_delta = 0.99),
                            seed = 1)

waic(bayes.alpha.elev.1fe, pointwise = TRUE)
loo(bayes.alpha.elev.1fe, pointwise = TRUE)

# SR ~ simple elevation multilevel lat effect model
bayes.alpha.elev.1re <- brm(formula = value ~ Elevation.scale + (1 + Elevation.scale | lat_class),
                            data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                            prior = c(set_prior("normal(0,5)", class = "b")),
                            warmup = 1000, iter = 3000, chains = 4,
                            control = list(adapt_delta = 0.99),
                            seed = 1)

waic(bayes.alpha.elev.1re, pointwise = TRUE)
loo(bayes.alpha.elev.1re, pointwise = TRUE)

# SR ~ poly elevation fixed effect model
bayes.alpha.elev.2fe <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE),
                            data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                            prior = c(set_prior("normal(0,5)", class = "b")),
                            warmup = 1000, iter = 3000, chains = 4,
                            control = list(adapt_delta = 0.99),
                            seed = 1)

waic(bayes.alpha.elev.2fe, pointwise = TRUE)
loo(bayes.alpha.elev.2fe, pointwise = TRUE)
summary(bayes.alpha.elev.2fe)

# SR ~ poly elevation fixed effect plotting
pra <- ggpredict(bayes.alpha.elev.2fe, terms = c("Elevation.scale [all]"), ppd = FALSE)

pra$x <- pra$x * sd(lmdata$Elevation) + mean(lmdata$Elevation)

(p1.alpha.ele.2fe <- ggplot() +
    geom_ribbon(data = pra, aes(ymin=conf.low, ymax=conf.high, x=x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = Elevation, y = value, colour = "red")) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Species richness", breaks = c(0, 5, 10, 15, 20, 25)) +
    scale_x_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    coord_cartesian(ylim = c(0, 26), xlim = c(-3.319242, 3750)) +
    workingtheme +
    theme(legend.position = "none") +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.alpha.ele.2fe, file = "p1.alpha.ele_post.linpred.2fe.rdata")

# SR ~ poly elevation multilevel lat effect model
bayes.alpha.elev.2re <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE) + (1 + poly(Elevation.scale, 2, raw = FALSE) | lat_class),
                            data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                            prior = c(set_prior("normal(0,5)", class = "b")),
                            warmup = 1000, iter = 3000, chains = 4,
                            control = list(adapt_delta = 0.99),
                            seed = 1)

waic(bayes.alpha.elev.2re, pointwise = TRUE)
loo(bayes.alpha.elev.2re, pointwise = TRUE)
summary(bayes.alpha.elev.2re)
prior_summary(bayes.alpha.elev.2re)
coef(bayes.alpha.elev.2re)
fixef(bayes.alpha.elev.2re)
bayes_R2(bayes.alpha.elev.2re)
control_params(bayes.alpha.elev.2re)

# SR ~ poly elevation multilevel lat effect plotting
pp_check(bayes.alpha.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.alpha.elev.2re, terms = c("Elevation.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.alpha.elev.2re, terms = c("Elevation.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 
prb$x <- prb$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 

(p1.alpha.ele <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = Elevation, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour=group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Species richness", breaks = c(0, 5, 10, 15, 20, 25)) +
    scale_x_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    coord_cartesian(ylim = c(0, 26), xlim = c(-3.319242, 3750)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.alpha.ele, file = "p1.alpha.ele_post.linpred.re2.rdata")

# SR ~ simple latitude fixed effect model
bayes.alpha.lat.1fe <- brm(formula = value ~ Latitude.scale,
                           data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                           prior = c(set_prior("normal(0,5)", class = "b")),
                           warmup = 1000, iter = 3000, chains = 4,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.alpha.lat.1fe, pointwise = TRUE)
loo(bayes.alpha.lat.1fe, pointwise = TRUE)

# SR ~ simple latitude multilevel elev effect model
bayes.alpha.lat.1re <- brm(formula = value ~ Latitude.scale + (1 + Latitude.scale | elev_class),
                           data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                           prior = c(set_prior("normal(0,5)", class = "b")),
                           warmup = 1000, iter = 3000, chains = 4,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.alpha.lat.1re, pointwise = TRUE)
loo(bayes.alpha.lat.1re, pointwise = TRUE)

# SR ~ poly latitude fixed effect model
bayes.alpha.lat.2fe <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE),
                           data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                           prior = c(set_prior("normal(0,5)", class = "b")),
                           warmup = 1000, iter = 3000, chains = 4,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.alpha.lat.2fe, pointwise = TRUE)
loo(bayes.alpha.lat.2fe, pointwise = TRUE)
summary(bayes.alpha.lat.2fe)

# SR ~ poly latitude fixed effect plotting
pra <- ggpredict(bayes.alpha.lat.2fe, terms = c("Latitude.scale [all]"), ppd = FALSE)

pra$x <- pra$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 

(p1.alpha.lat.2fe <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = Latitude, y = value, colour = "red")) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Species richness", breaks = c(0, 5, 10, 15, 20, 25)) +
    scale_x_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    coord_cartesian(ylim = c(0, 26), xlim = c(36, 66)) +
    workingtheme +
    theme(legend.position = "none") +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.alpha.lat.2fe, file = "p1.alpha.lat_post.linpred.2fe.rdata")

# SR ~ poly latitude multilevel elev effect model
bayes.alpha.lat.2re <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE) + (1 + poly(Latitude.scale, 2, raw = FALSE) | elev_class),
                           data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                           prior = c(set_prior("normal(0,5)", class = "b")),
                           warmup = 1000, iter = 3000, chains = 4,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.alpha.lat.2re, pointwise = TRUE)
loo(bayes.alpha.lat.2re, pointwise = TRUE)
summary(bayes.alpha.lat.2re)
prior_summary(bayes.alpha.lat.2re)
coef(bayes.alpha.lat.2re)
fixef(bayes.alpha.lat.2re)
bayes_R2(bayes.alpha.lat.2re)
control_params(bayes.alpha.lat.2re)

# SR ~ poly latitude multilevel elev effect plotting
pp_check(bayes.alpha.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.alpha.lat.2re, terms = c("Latitude.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.alpha.lat.2re, terms = c("Latitude.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 
prb$x <- prb$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 

(p1.alpha.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = Latitude, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Species richness", breaks = c(0, 5, 10, 15, 20, 25)) +
    scale_x_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    coord_cartesian(ylim = c(0, 26), xlim = c(36, 66)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values=a60505_6e068f))

save(p1.alpha.lat, file = "p1.alpha.lat_post.linpred.re2.rdata")


## Model phylogenetic richness

# Check distribution
dist <- subset(lmdata, variable == "Glob.PD")
hist(test$value)

# PD ~ simple elevation fixed effect model
bayes.pd.elev.1fe <- brm(formula = value ~ Elevation.scale,
                         data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                         prior = c(set_prior("normal(0,5)", class = "b")),
                         warmup = 1000, iter = 3000, chains = 4,
                         control = list(adapt_delta = 0.99),
                         seed = 1)

waic(bayes.pd.elev.1fe, pointwise = TRUE)
loo(bayes.pd.elev.1fe, pointwise = TRUE)

# PD ~ simple elevation multilevel lat effect model
bayes.pd.elev.1re <- brm(formula = value ~ Elevation.scale + (1 + Elevation.scale | lat_class),
                         data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                         prior = c(set_prior("normal(0,5)", class = "b")),
                         warmup = 1000, iter = 3000, chains = 4,
                         control = list(adapt_delta = 0.99),
                         seed = 1)

waic(bayes.pd.elev.1re, pointwise = TRUE)
loo(bayes.pd.elev.1re, pointwise = TRUE)

# PD ~ poly elevation fixed effect model
bayes.pd.elev.2fe <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE),
                         data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                         prior = c(set_prior("normal(0,5)", class = "b")),
                         warmup = 1000, iter = 3000, chains = 4,
                         control = list(adapt_delta = 0.99),
                         seed = 1)

waic(bayes.pd.elev.2fe, pointwise = TRUE)
loo(bayes.pd.elev.2fe, pointwise = TRUE)

# PD ~ poly elevation multilevel lat effect model
bayes.pd.elev.2re <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE) + (1 + poly(Elevation.scale, 2, raw = FALSE) | lat_class),
                         data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                         prior = c(set_prior("normal(0,5)", class = "b")),
                         warmup = 1000, iter = 3000, chains = 4,
                         control = list(adapt_delta = 0.99),
                         seed = 1)

waic(bayes.pd.elev.2re, pointwise = TRUE)
loo(bayes.pd.elev.2re, pointwise = TRUE)
summary(bayes.pd.elev.2re)
prior_summary(bayes.pd.elev.2re)
coef(bayes.pd.elev.2re)
fixef(bayes.pd.elev.2re)
bayes_R2(bayes.pd.elev.2re)
control_params(bayes.pd.elev.2re)

# PD ~ poly elevation multilevel lat effect plotting
pp_check(bayes.pd.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.elev.2re, terms = c("Elevation.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.elev.2re, terms = c("Elevation.scale [all]", "lat_class"), type = "re", ppd = FALSE)


raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 
prb$x <- prb$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 

(p1.pd.ele <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.PD"), alpha = 0.1, aes(x = Elevation, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Phylogenetic richness", breaks = c(0, 160, 320, 480, 640)) +
    scale_x_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    coord_cartesian(ylim = c(0, 647.1955), xlim = c(-3.319242, 3750)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.pd.ele, file = "p1.pd.ele_post.linpred.re2.rdata")

#  PD ~ simple latitude fixed effect model
bayes.pd.lat.1fe <- brm(formula = value ~ Latitude.scale,
                        data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                        prior = c(set_prior("normal(0,5)", class = "b")),
                        warmup = 1000, iter = 3000, chains = 4,
                        control = list(adapt_delta = 0.99),
                        seed = 1)

waic(bayes.pd.lat.1fe, pointwise = TRUE)
loo(bayes.pd.lat.1fe, pointwise = TRUE)

# PD ~ simple latitude multilevel elev effect model
bayes.pd.lat.1re <- brm(formula = value ~ Latitude.scale + (1 + Latitude.scale | elev_class),
                        data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                        prior = c(set_prior("normal(0,5)", class = "b")),
                        warmup = 1000, iter = 3000, chains = 4,
                        control = list(adapt_delta = 0.99),
                        seed = 1)

waic(bayes.pd.lat.1re, pointwise = TRUE)
loo(bayes.pd.lat.1re, pointwise = TRUE)

# PD ~ poly latitude fixed effect model
bayes.pd.lat.2fe <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE),
                        data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                        prior = c(set_prior("normal(0,5)", class = "b")),
                        warmup = 1000, iter = 3000, chains = 4,
                        control = list(adapt_delta = 0.99),
                        seed = 1)

waic(bayes.pd.lat.2fe, pointwise = TRUE)
loo(bayes.pd.lat.2fe, pointwise = TRUE)

# PD ~ poly latitude multilevel elev effect model
bayes.pd.lat.2re <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE) + (1 + poly(Latitude.scale, 2, raw = FALSE) | elev_class),
                        data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                        prior = c(set_prior("normal(0,5)", class = "b")),
                        warmup = 1000, iter = 3000, chains = 4,
                        control = list(adapt_delta = 0.99),
                        seed = 1)

waic(bayes.pd.lat.2re, pointwise = TRUE)
loo(bayes.pd.lat.2re, pointwise = TRUE)
summary(bayes.pd.lat.2re)
prior_summary(bayes.pd.lat.2re)
coef(bayes.pd.lat.2re)
fixef(bayes.pd.lat.2re)
bayes_R2(bayes.pd.lat.2re)
control_params(bayes.pd.lat.2re)

# PD ~ poly latitude multilevel elev effect plotting
pp_check(bayes.pd.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.lat.2re, terms = c("Latitude.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.lat.2re, terms = c("Latitude.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 
prb$x <- prb$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 

(p1.pd.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.PD"), alpha = 0.1, aes(x = Latitude, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Phylogenetic richness", breaks = c(0, 160, 320, 480, 640)) +
    scale_x_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    coord_cartesian(ylim = c(0, 647.1955), xlim = c(36, 66)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.pd.lat, file = "p1.pd.lat_post.linpred.re2.rdata")


## Model raw/global phylogenetic mean pairwise distance (bounded by 0 and 1)

# Check distribution
dist <- subset(lmdata, variable == "Glob.MPD.beta")
hist(dist$value)

# MPD ~ simple elevation fixed effect model
bayes.mpd.elev.1fe <- brm(formula = value ~ Elevation.scale,
                          data = subset(lmdata, variable == "Glob.MPD.beta"), family = zero_one_inflated_beta(link = "logit"),
                          prior = c(set_prior("normal(0,5)", class = "b")),
                          warmup = 1000, iter = 3000, chains = 4,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.mpd.elev.1fe, pointwise = TRUE)
loo(bayes.mpd.elev.1fe, pointwise = TRUE)

# MPD ~ simple elevation multilevel lat effect model
bayes.mpd.elev.1re <- brm(formula = value ~ Elevation.scale + (1 + Elevation.scale | lat_class),
                          data = subset(lmdata, variable == "Glob.MPD.beta"), family = zero_one_inflated_beta(link = "logit"),
                          prior = c(set_prior("normal(0,5)", class = "b")),
                          warmup = 1000, iter = 3000, chains = 4,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.mpd.elev.1re, pointwise = TRUE)
loo(bayes.mpd.elev.1re, pointwise = TRUE)

# MPD ~ poly elevation fixed effect model
bayes.mpd.elev.2fe <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE),
                          data = subset(lmdata, variable == "Glob.MPD.beta"), family = zero_one_inflated_beta(link = "logit"),
                          prior = c(set_prior("normal(0,5)", class = "b")),
                          warmup = 1000, iter = 3000, chains = 4,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.mpd.elev.2fe, pointwise = TRUE)
loo(bayes.mpd.elev.2fe, pointwise = TRUE)

# MPD ~ poly elevation multilevel lat effect model
bayes.mpd.elev.2re <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE) + (1 + poly(Elevation.scale, 2, raw = FALSE) | lat_class),
                          data = subset(lmdata, variable == "Glob.MPD.beta"), family = zero_one_inflated_beta(link = "logit"),
                          prior = c(set_prior("normal(0,5)", class = "b")),
                          warmup = 1000, iter = 3000, chains = 4,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.mpd.elev.2re, pointwise = TRUE)
loo(bayes.mpd.elev.2re, pointwise = TRUE)
summary(bayes.mpd.elev.2re)
prior_summary(bayes.mpd.elev.2re)
coef(bayes.mpd.elev.2re)
fixef(bayes.mpd.elev.2re)
bayes_R2(bayes.mpd.elev.2re)
control_params(bayes.mpd.elev.2re)

# MPD ~ poly elevation multilevel lat effect plotting
pp_check(bayes.mpd.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.mpd.elev.2re, terms = c("Elevation.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.mpd.elev.2re, terms = c("Elevation.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 
prb$x <- prb$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 

(p1.mpd.ele <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.MPD.beta"), alpha = 0.1, aes(x = Elevation, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Phylogenetic mean pairwise distance", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_x_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    coord_cartesian(ylim = c(0, 1), xlim = c(-3.319242, 3750)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.mpd.ele, file = "p1.mpd.ele_post.linpred.re2.rdata")

# MPD ~ simple latitude fixed effect model
bayes.mpd.lat.1fe <- brm(formula = value ~ Latitude.scale,
                         data = subset(lmdata, variable == "Glob.MPD.beta"), family = zero_one_inflated_beta(link = "logit"),
                         prior = c(set_prior("normal(0,5)", class = "b")),
                         warmup = 1000, iter = 3000, chains = 4,
                         control = list(adapt_delta = 0.99),
                         seed=1)

waic(bayes.mpd.lat.1fe, pointwise = TRUE)
loo(bayes.mpd.lat.1fe, pointwise = TRUE)

# MPD ~ simple latitude multilevel elev effect model
bayes.mpd.lat.1re <- brm(formula = value ~ Latitude.scale + (1 + Latitude.scale | elev_class),
                         data = subset(lmdata, variable == "Glob.MPD.beta"), family = zero_one_inflated_beta(link = "logit"),
                         prior = c(set_prior("normal(0,5)", class = "b")),
                         warmup = 1000, iter = 3000, chains = 4,
                         control = list(adapt_delta = 0.99),
                         seed = 1)

waic(bayes.mpd.lat.1re, pointwise = TRUE)
loo(bayes.mpd.lat.1re, pointwise = TRUE)

# MPD ~ poly latitude fixed effect model
bayes.mpd.lat.2fe <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE),
                         data = subset(lmdata, variable == "Glob.MPD.beta"), family = zero_one_inflated_beta(link = "logit"),
                         prior = c(set_prior("normal(0,5)", class = "b")),
                         warmup = 1000, iter = 3000, chains = 4,
                         control = list(adapt_delta = 0.99),
                         seed = 1)

waic(bayes.mpd.lat.2fe, pointwise = TRUE)
loo(bayes.mpd.lat.2fe, pointwise = TRUE)

# MPD ~ poly latitude multilevel elev effect model
bayes.mpd.lat.2re <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE) + (1 + poly(Latitude.scale, 2, raw = FALSE) | elev_class),
                         data = subset(lmdata, variable == "Glob.MPD.beta"), family = zero_one_inflated_beta(link = "logit"),
                         prior = c(set_prior("normal(0,5)", class = "b")),
                         warmup = 1000, iter = 3000, chains = 4,
                         control = list(adapt_delta = 0.99),
                         seed = 1)

waic(bayes.mpd.lat.2re, pointwise = TRUE)
loo(bayes.mpd.lat.2re, pointwise = TRUE)
summary(bayes.mpd.lat.2re)
prior_summary(bayes.mpd.lat.2re)
coef(bayes.mpd.lat.2re)
fixef(bayes.mpd.lat.2re)
bayes_R2(bayes.mpd.lat.2re)
control_params(bayes.mpd.lat.2re)

# MPD ~ poly latitude multilevel elev effect plotting
pp_check(bayes.mpd.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.mpd.lat.2re, terms = c("Latitude.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.mpd.lat.2re, terms = c("Latitude.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Latitude) + mean(lmdata$Latitude)
prb$x <- prb$x * sd(lmdata$Latitude) + mean(lmdata$Latitude)

(p1.mpd.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin=conf.low, ymax=conf.high, x=x), alpha=0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.MPD.beta"), alpha = 0.1, aes(x = Latitude, y = value, colour=elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour=group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Phylogenetic mean pairwise distance", breaks=c(0,0.2,0.4,0.6,0.8,1)) +
    scale_x_continuous("Latitude (°N)", breaks=c(36,42,48,54,60,66)) +
    coord_cartesian(ylim = c(0, 1), xlim = c(36, 66)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.mpd.lat, file = "p1.mpd.lat_post.linpred.re2.rdata")


## Model functional richness

# Check distribution
dist <- subset(lmdata, variable == "Glob.FRic")
hist(dist$value)

# FRic ~ simple elevation fixed effect model
mix <- mixture(Beta, Beta)
bayes.FRic.elev.1fe <- brm(formula = value ~ Elevation.scale,
                           data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                           prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                     prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                           warmup = 1000, iter = 3000, chains = 4, inits = 0,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.FRic.elev.1fe, pointwise = TRUE)
loo(bayes.FRic.elev.1fe, pointwise = TRUE)

# FRic ~ simple elevation multilevel lat effect model
mix <- mixture(Beta, Beta)
bayes.FRic.elev.1re <- brm(formula = value ~ Elevation.scale + (1 + Elevation.scale | lat_class),
                           data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                           prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                     prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                           warmup = 1000, iter = 3000, chains = 4, inits = 0,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.FRic.elev.1re, pointwise = TRUE)
loo(bayes.FRic.elev.1re, pointwise = TRUE)

# FRic ~ poly elevation fixed effect model
mix <- mixture(Beta, Beta)
bayes.FRic.elev.2fe <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE),
                           data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                           prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                     prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                           warmup = 1000, iter = 3000, chains = 4, inits = 0,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.FRic.elev.2fe, pointwise = TRUE)
loo(bayes.FRic.elev.2fe, pointwise = TRUE)

# FRic ~ poly elevation multilevel lat effect model
mix <- mixture(Beta, Beta)
bayes.FRic.elev.2re <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE) + (1 + poly(Elevation.scale, 2, raw = FALSE) | lat_class),
                           data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                           prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                     prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                           warmup = 1000, iter = 9000, chains = 4, inits = 0,
                           control = list(adapt_delta = 0.99, max_treedepth = 15),
                           seed = 1)

waic(bayes.FRic.elev.2re, pointwise = TRUE)
loo(bayes.FRic.elev.2re, pointwise = TRUE)
summary(bayes.FRic.elev.2re)
prior_summary(bayes.FRic.elev.2re)
coef(bayes.FRic.elev.2re)
fixef(bayes.FRic.elev.2re)
bayes_R2(bayes.FRic.elev.2re)
control_params(bayes.FRic.elev.2re)

# FRic ~ poly elevation multilevel lat effect plotting
pp_check(bayes.FRic.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.elev.2re, terms = c("Elevation.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.elev.2re, terms = c("Elevation.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 
prb$x <- prb$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 

(p1.FRic.ele <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.FRic"), alpha = 0.1, aes(x = Elevation, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Functional richness", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_x_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    coord_cartesian(ylim = c(0, 1), xlim = c(-3.319242, 3750)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.FRic.ele, file = "p1.FRic.ele_post.linpred.re2.rdata")

# FRic ~ simple latitude fixed effect model
mix <- mixture(Beta, Beta)
bayes.FRic.lat.1fe <- brm(formula = value ~ Latitude.scale,
                          data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                          prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                    prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                          warmup = 1000, iter = 3000, chains = 4, inits = 0,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.FRic.lat.1fe, pointwise = TRUE)
loo(bayes.FRic.lat.1fe, pointwise = TRUE) 

# FRic ~ simple latitude multilevel elev effect model
mix <- mixture(Beta, Beta)
bayes.FRic.lat.1re <- brm(formula = value ~ Latitude.scale + (1 + Latitude.scale | elev_class),
                          data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                          prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                    prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                          warmup = 1000, iter = 3000, chains = 4, inits = 0,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.FRic.lat.1re, pointwise = TRUE)
loo(bayes.FRic.lat.1re, pointwise = TRUE)

# FRic ~ poly latitude fixed effect model
mix <- mixture(Beta, Beta)
bayes.FRic.lat.2fe <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE),
                          data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                          prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                    prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                          warmup = 1000, iter = 3000, chains = 4, inits = 0,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.FRic.lat.2fe, pointwise = TRUE)
loo(bayes.FRic.lat.2fe, pointwise = TRUE)

# FRic ~ poly latitude multilevel elev effect model
mix <- mixture(Beta, Beta)
bayes.FRic.lat.2re <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE) + (1 + poly(Latitude.scale, 2, raw = FALSE) | elev_class),
                          data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                          prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                    prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                          warmup = 1000, iter = 9000, chains = 4, inits = 0,
                          control = list(adapt_delta = 0.99, max_treedepth = 15),
                          seed = 1)

waic(bayes.FRic.lat.2re, pointwise = TRUE)
loo(bayes.FRic.lat.2re, pointwise = TRUE)
summary(bayes.FRic.lat.2re)
prior_summary(bayes.FRic.lat.2re)
coef(bayes.FRic.lat.2re)
fixef(bayes.FRic.lat.2re)
bayes_R2(bayes.FRic.lat.2re)
control_params(bayes.FRic.lat.2re)

# FRic ~ poly latitude multilevel elev effect plotting
pp_check(bayes.FRic.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.lat.2re, terms = c("Latitude.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.lat.2re, terms = c("Latitude.scale [all]", "elev_class"), type = "re", ppd =FALSE)
 
raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 
prb$x <- prb$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 

(p1.FRic.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.FRic"), alpha = 0.1, aes(x = Latitude, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Functional richness", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_x_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    coord_cartesian(ylim = c(0, 1), xlim = c(36, 66)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.FRic.lat, file = "p1.FRic.lat_post.linpred.re2.rdata")


## Model global/raw functional dispersion

# Check distribution
dist <- subset(lmdata, variable == "Glob.FDis")
hist(dist$value)

# FDis ~ simple elevation fixed effect model
bayes.FDis.elev.1fe <- brm(formula = value ~ Elevation.scale,
                           data = subset(lmdata, variable == "Glob.FDis"), family = Gamma(link = "log"),
                           prior = c(set_prior("normal(0,5)", class = "b")),
                           warmup = 1000, iter = 3000, chains = 4,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.FDis.elev.1fe, pointwise = TRUE)
loo(bayes.FDis.elev.1fe, pointwise = TRUE)

# FDis ~ simple elevation multilevel lat effect model
bayes.FDis.elev.1re <- brm(formula = value ~ Elevation.scale + (1 + Elevation.scale | lat_class),
                           data = subset(lmdata, variable == "Glob.FDis"), family = Gamma(link = "log"),
                           prior = c(set_prior("normal(0,5)", class = "b")),
                           warmup = 1000, iter = 3000, chains = 4,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.FDis.elev.1re, pointwise = TRUE)
loo(bayes.FDis.elev.1re, pointwise = TRUE)

# FDis ~ poly elevation fixed effect model
bayes.FDis.elev.2fe <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE),
                           data = subset(lmdata, variable == "Glob.FDis"), family = Gamma(link = "log"),
                           prior = c(set_prior("normal(0,5)", class = "b")),
                           warmup = 1000, iter = 3000, chains = 4,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.FDis.elev.2fe, pointwise = TRUE)
loo(bayes.FDis.elev.2fe, pointwise = TRUE)

# FDis ~ poly elevation multilevel lat effect model
bayes.FDis.elev.2re <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE) + (1 + poly(Elevation.scale, 2, raw = FALSE) | lat_class),
                           data = subset(lmdata, variable == "Glob.FDis"), family = Gamma(link = "log"),
                           prior = c(set_prior("normal(0,5)", class = "b")),
                           warmup = 4000, iter = 12000, chains = 4,
                           control = list(adapt_delta = 0.99, max_treedepth = 15),
                           seed = 1)

waic(bayes.FDis.elev.2re, pointwise = TRUE)
loo(bayes.FDis.elev.2re, pointwise = TRUE)
summary(bayes.FDis.elev.2re)
prior_summary(bayes.FDis.elev.2re)
coef(bayes.FDis.elev.2re)
fixef(bayes.FDis.elev.2re)
bayes_R2(bayes.FDis.elev.2re)
control_params(bayes.FDis.elev.2re)

# FDis ~ poly elevation multilevel lat effect plotting
pp_check(bayes.FDis.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FDis.elev.2re, terms = c("Elevation.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FDis.elev.2re, terms = c("Elevation.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 
prb$x <- prb$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 

(p1.FDis.ele <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.FDis"), alpha = 0.1, aes(x = Elevation, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Functional dispersion", breaks = c(0, 0.06, 0.12, 0.18, 0.24, 0.3)) +
    scale_x_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    coord_cartesian(ylim = c(0, 0.3), xlim = c(-3.319242, 3750)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.FDis.ele, file = "p1.FDis.ele_post.linpred.re2.rdata")

# FDis ~ simple latitude fixed effect model
bayes.FDis.lat.1fe <- brm(formula = value ~ Latitude.scale,
                          data = subset(lmdata, variable == "Glob.FDis"), family = Gamma(link = "log"),
                          prior = c(set_prior("normal(0,5)", class = "b")),
                          warmup = 1000, iter = 3000, chains = 4,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.FDis.lat.1fe, pointwise = TRUE)
loo(bayes.FDis.lat.1fe, pointwise = TRUE)

# FDis ~ simple latitude multilevel elev effect model
bayes.FDis.lat.1re <- brm(formula = value ~ Latitude.scale + (1 + Latitude.scale | elev_class),
                          data = subset(lmdata, variable == "Glob.FDis"), family = Gamma(link = "log"),
                          prior = c(set_prior("normal(0,5)", class = "b")),
                          warmup = 1000, iter = 3000, chains = 4,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.FDis.lat.1re, pointwise = TRUE)
loo(bayes.FDis.lat.1re, pointwise = TRUE)

# FDis ~ poly latitude fixed effect model
bayes.FDis.lat.2fe <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE),
                          data = subset(lmdata, variable == "Glob.FDis"), family = Gamma(link = "log"),
                          prior = c(set_prior("normal(0,5)", class = "b")),
                          warmup = 1000, iter = 3000, chains = 4,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.FDis.lat.2fe, pointwise = TRUE)
loo(bayes.FDis.lat.2fe, pointwise = TRUE)

# FDis ~ poly latitude multilevel elev effect model
bayes.FDis.lat.2re <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE) + (1 + poly(Latitude.scale, 2, raw = FALSE) | elev_class),
                          data = subset(lmdata, variable == "Glob.FDis"), family = Gamma(link = "log"),
                          prior = c(set_prior("normal(0,5)", class = "b")),
                          warmup = 1000, iter = 3000, chains = 4,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.FDis.lat.2re, pointwise = TRUE)
loo(bayes.FDis.lat.2re, pointwise = TRUE)
summary(bayes.FDis.lat.2re)
prior_summary(bayes.FDis.lat.2re)
coef(bayes.FDis.lat.2re)
fixef(bayes.FDis.lat.2re)
bayes_R2(bayes.FDis.lat.2re)
control_params(bayes.FDis.lat.2re)

# FDis ~ poly latitude multilevel elev effect plotting
pp_check(bayes.FDis.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FDis.lat.2re, terms = c("Latitude.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FDis.lat.2re, terms = c("Latitude.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 
prb$x <- prb$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 

(p1.FDis.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.FDis"), alpha = 0.1, aes(x = Latitude, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Functional dispersion", breaks = c(0, 0.06, 0.12, 0.18, 0.24, 0.3)) +
    scale_x_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    coord_cartesian(ylim = c(0, 0.3), xlim = c(36, 66)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.FDis.lat, file = "p1.FDis.lat_post.linpred.re2.rdata")


## Model community weighted mean body lengths

# Check distribution
dist <- subset(lmdata, variable == "CWM.Bl")
hist(dist$value)

# CWM.Bl ~ simple elevation fixed effect model
bayes.CWM.Bl.elev.1fe <- brm(formula = value ~ Elevation.scale,
                             data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                             prior = c(set_prior("normal(0,5)", class = "b")),
                             warmup = 1000, iter = 3000, chains = 4,
                             control = list(adapt_delta = 0.99),
                             seed = 1)

waic(bayes.CWM.Bl.elev.1fe, pointwise = TRUE)
loo(bayes.CWM.Bl.elev.1fe, pointwise = TRUE)

# CWM.Bl ~ simple elevation multilevel lat effect model
bayes.CWM.Bl.elev.1re <- brm(formula = value ~ Elevation.scale + (1 + Elevation.scale | lat_class),
                             data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                             prior = c(set_prior("normal(0,5)", class = "b")),
                             warmup = 1000, iter = 3000, chains = 4,
                             control = list(adapt_delta = 0.99),
                             seed = 1)

waic(bayes.CWM.Bl.elev.1re, pointwise = TRUE)
loo(bayes.CWM.Bl.elev.1re, pointwise = TRUE)

# CWM.Bl ~ poly elevation fixed effect model
bayes.CWM.Bl.elev.2fe <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE),
                             data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                             prior = c(set_prior("normal(0,5)", class = "b")),
                             warmup = 1000, iter = 3000, chains = 4,
                             control = list(adapt_delta = 0.99),
                             seed = 1)

waic(bayes.CWM.Bl.elev.2fe, pointwise = TRUE)
loo(bayes.CWM.Bl.elev.2fe, pointwise = TRUE)

# CWM.Bl ~ poly elevation multilevel lat effect model
bayes.CWM.Bl.elev.2re <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE) + (1 + poly(Elevation.scale, 2, raw = FALSE) | lat_class),
                             data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                             prior = c(set_prior("normal(0,5)", class = "b")),
                             warmup = 1000, iter = 3000, chains = 4,
                             control = list(adapt_delta = 0.99),
                             seed = 1)

waic(bayes.CWM.Bl.elev.2re, pointwise = TRUE)
loo(bayes.CWM.Bl.elev.2re, pointwise = TRUE)
summary(bayes.CWM.Bl.elev.2re)
prior_summary(bayes.CWM.Bl.elev.2re)
coef(bayes.CWM.Bl.elev.2re)
fixef(bayes.CWM.Bl.elev.2re)
bayes_R2(bayes.CWM.Bl.elev.2re)
control_params(bayes.CWM.Bl.elev.2re)

# CWM.Bl ~ poly elevation multilevel lat effect plotting
pp_check(bayes.CWM.Bl.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.CWM.Bl.elev.2re, terms = c("Elevation.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.CWM.Bl.elev.2re, terms = c("Elevation.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 
prb$x <- prb$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 

(p1.CWM.Bl.ele <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "CWM.Bl"), alpha = 0.1, aes(x = Elevation, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Community mean body length", breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
    scale_x_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    coord_cartesian(ylim = c(0, 3.5), xlim = c(-3.319242, 3750)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.CWM.Bl.ele, file = "p1.CWM.Bl.ele_post.linpred.re2.rdata")

# CWM.Bl ~ simple latitude fixed effect model
bayes.CWM.Bl.lat.1fe <- brm(formula = value ~ Latitude.scale,
                            data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                            prior = c(set_prior("normal(0,5)", class = "b")),
                            warmup = 1000, iter = 3000, chains = 4,
                            control = list(adapt_delta = 0.99),
                            seed = 1)

waic(bayes.CWM.Bl.lat.1fe, pointwise = TRUE)
loo(bayes.CWM.Bl.lat.1fe, pointwise = TRUE)

# CWM.Bl ~ simple latitude multilevel elev effect model
bayes.CWM.Bl.lat.1re <- brm(formula = value ~ Latitude.scale + (1 + Latitude.scale | elev_class),
                            data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                            prior = c(set_prior("normal(0,5)", class = "b")),
                            warmup = 1000, iter = 3000, chains = 4,
                            control = list(adapt_delta = 0.99),
                            seed = 1)

waic(bayes.CWM.Bl.lat.1re, pointwise = TRUE)
loo(bayes.CWM.Bl.lat.1re, pointwise = TRUE)

# CWM.Bl ~ poly latitude fixed effect model
bayes.CWM.Bl.lat.2fe <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE),
                            data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                            prior = c(set_prior("normal(0,5)", class = "b")),
                            warmup = 1000, iter = 3000, chains = 4,
                            control = list(adapt_delta = 0.99),
                            seed = 1)

waic(bayes.CWM.Bl.lat.2fe, pointwise = TRUE)
loo(bayes.CWM.Bl.lat.2fe, pointwise = TRUE)

# CWM.Bl ~ poly latitude multilevel elev effect model
bayes.CWM.Bl.lat.2re <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE) + (1 + poly(Latitude.scale, 2, raw = FALSE) | elev_class),
                            data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                            prior = c(set_prior("normal(0,5)", class = "b")),
                            warmup = 1000, iter = 3000, chains = 4,
                            control = list(adapt_delta = 0.99),
                            seed = 1)

waic(bayes.CWM.Bl.lat.2re, pointwise = TRUE)
loo(bayes.CWM.Bl.lat.2re, pointwise = TRUE)
summary(bayes.CWM.Bl.lat.2re)
prior_summary(bayes.CWM.Bl.lat.2re)
coef(bayes.CWM.Bl.lat.2re)
fixef(bayes.CWM.Bl.lat.2re)
bayes_R2(bayes.CWM.Bl.lat.2re)
control_params(bayes.CWM.Bl.lat.2re)

# CWM.Bl ~ poly latitude multilevel elev effect plotting
pp_check(bayes.CWM.Bl.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.CWM.Bl.lat.2re, terms = c("Latitude.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.CWM.Bl.lat.2re, terms = c("Latitude.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 
prb$x <- prb$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 

(p1.CWM.Bl.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "CWM.Bl"), alpha = 0.1, aes(x = Latitude, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Community mean body length", breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
    scale_x_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    coord_cartesian(ylim = c(0, 3.5), xlim = c(36, 66)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.CWM.Bl.lat, file = "p1.CWM.Bl.lat_post.linpred.re2.rdata")


## Model standardized phylogenetic richness

# Check distribution
dist <- subset(lmdata, variable == "pd.obs.z.lat")
hist(dist$value)

dist <- subset(lmdata, variable == "pd.obs.z.ele")
hist(dist$value)

# PD.obs.z ~ simple elevation fixed effect model
bayes.pd.obs.z.lat.elev.1fe <- brm(formula = value ~ Elevation.scale,
                                   data = subset(lmdata, variable == "pd.obs.z.lat"),family = skew_normal(),
                                   prior = c(set_prior("normal(0,5)", class = "b")),
                                   warmup = 1000, iter = 3000, chains = 4,
                                   control = list(adapt_delta = 0.99),
                                   seed = 1)

waic(bayes.pd.obs.z.lat.elev.1fe, pointwise = TRUE)
loo(bayes.pd.obs.z.lat.elev.1fe, pointwise = TRUE)

# PD.obs.z ~ simple elevation multilevel lat effect model
bayes.pd.obs.z.lat.elev.1re <- brm(formula = value ~ Elevation.scale + (1 + Elevation.scale | lat_class),
                                   data = subset(lmdata, variable == "pd.obs.z.lat"),family = skew_normal(),
                                   prior = c(set_prior("normal(0,5)", class = "b")),
                                   warmup = 1000, iter = 3000, chains = 4,
                                   control = list(adapt_delta = 0.99),
                                   seed = 1)

waic(bayes.pd.obs.z.lat.elev.1re, pointwise = TRUE)
loo(bayes.pd.obs.z.lat.elev.1re, pointwise = TRUE)

# PD.obs.z ~ poly elevation fixed effect model
bayes.pd.obs.z.lat.elev.2fe <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE),
                                   data = subset(lmdata, variable == "pd.obs.z.lat"),family = skew_normal(),
                                   prior = c(set_prior("normal(0,5)", class = "b")),
                                   warmup = 1000, iter = 3000, chains = 4,
                                   control = list(adapt_delta = 0.99),
                                   seed = 1)

waic(bayes.pd.obs.z.lat.elev.2fe, pointwise = TRUE)
loo(bayes.pd.obs.z.lat.elev.2fe, pointwise = TRUE)

# PD.obs.z ~ poly elevation multilevel lat effect model
bayes.pd.obs.z.lat.elev.2re <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE) + (1 + poly(Elevation.scale, 2, raw = FALSE) | lat_class),
                                   data = subset(lmdata, variable == "pd.obs.z.lat"),family = skew_normal(),
                                   prior = c(set_prior("normal(0,5)", class = "b")),
                                   warmup = 1000, iter = 3000, chains = 4,
                                   control = list(adapt_delta = 0.99),
                                   seed = 1)

waic(bayes.pd.obs.z.lat.elev.2re, pointwise = TRUE)
loo(bayes.pd.obs.z.lat.elev.2re, pointwise = TRUE)
summary(bayes.pd.obs.z.lat.elev.2re)
prior_summary(bayes.pd.obs.z.lat.elev.2re)
coef(bayes.pd.obs.z.lat.elev.2re)
fixef(bayes.pd.obs.z.lat.elev.2re)
bayes_R2(bayes.pd.obs.z.lat.elev.2re)
control_params(bayes.pd.obs.z.lat.elev.2re)

# PD.obs.z ~ poly elevation multilevel lat effect plotting
pp_check(bayes.pd.obs.z.lat.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.obs.z.lat.elev.2re, terms = c("Elevation.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.obs.z.lat.elev.2re, terms = c("Elevation.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 
prb$x <- prb$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 

(p1.pd.obs.z.lat.elev <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "pd.obs.z.lat"), alpha = 0.1, aes(x = Elevation, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic richness", breaks = c(-5, -3, -1, 1, 3)) +
    scale_x_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    coord_cartesian(ylim = c(-4.75, 3.361625), xlim = c(-3.319242, 3750)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.pd.obs.z.lat.elev, file = "p1.pd.obs.z.lat.elev_post.linpred.re2.rdata")

# PD.obs.z ~ simple latitude fixed effect model
bayes.pd.obs.z.ele.lat.1fe <- brm(formula = value ~ Latitude.scale,
                                  data = subset(lmdata, variable == "pd.obs.z.ele"),family = skew_normal(),
                                  prior = c(set_prior("normal(0,5)", class = "b")),
                                  warmup = 1000, iter = 3000, chains = 4,
                                  control = list(adapt_delta = 0.99),
                                  seed = 1)

waic(bayes.pd.obs.z.ele.lat.1fe, pointwise = TRUE)
loo(bayes.pd.obs.z.ele.lat.1fe, pointwise = TRUE)

# PD.obs.z ~ simple latitude multilevel elev effect model
bayes.pd.obs.z.ele.lat.1re <- brm(formula = value ~ Latitude.scale + (1 + Latitude.scale | elev_class),
                                  data = subset(lmdata, variable == "pd.obs.z.ele"),family = skew_normal(),
                                  prior = c(set_prior("normal(0,5)", class = "b")),
                                  warmup = 1000, iter = 3000, chains = 4,
                                  control = list(adapt_delta = 0.99),
                                  seed = 1)

waic(bayes.pd.obs.z.ele.lat.1re, pointwise = TRUE)
loo(bayes.pd.obs.z.ele.lat.1re, pointwise = TRUE)

# PD.obs.z ~ poly latitude fixed effect model
bayes.pd.obs.z.ele.lat.2fe <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE),
                                  data = subset(lmdata, variable == "pd.obs.z.ele"),family = skew_normal(),
                                  prior = c(set_prior("normal(0,5)", class = "b")),
                                  warmup = 1000, iter = 3000, chains = 4,
                                  control = list(adapt_delta = 0.99),
                                  seed = 1)

waic(bayes.pd.obs.z.ele.lat.2fe, pointwise = TRUE)
loo(bayes.pd.obs.z.ele.lat.2fe, pointwise = TRUE)

# PD.obs.z ~ poly latitude multilevel elev effect model
bayes.pd.obs.z.ele.lat.2re <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE) + (1 + poly(Latitude.scale, 2, raw = FALSE) | elev_class),
                                  data = subset(lmdata, variable == "pd.obs.z.ele"),family = skew_normal(),
                                  prior = c(set_prior("normal(0,5)", class = "b")),
                                  warmup = 1000, iter = 3000, chains = 4,
                                  control = list(adapt_delta = 0.99),
                                  seed = 1)

waic(bayes.pd.obs.z.ele.lat.2re, pointwise = TRUE)
loo(bayes.pd.obs.z.ele.lat.2re, pointwise = TRUE)
summary(bayes.pd.obs.z.ele.lat.2re)
prior_summary(bayes.pd.obs.z.ele.lat.2re)
coef(bayes.pd.obs.z.ele.lat.2re)
fixef(bayes.pd.obs.z.ele.lat.2re)
bayes_R2(bayes.pd.obs.z.ele.lat.2re)
control_params(bayes.pd.obs.z.ele.lat.2re)

# PD.obs.z ~ poly latitude multilevel elev effect plotting
pp_check(bayes.pd.obs.z.ele.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.obs.z.ele.lat.2re, terms = c("Latitude.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.obs.z.ele.lat.2re, terms = c("Latitude.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 
prb$x <- prb$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 

(p1.pd.obs.z.ele.lat <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "pd.obs.z.ele"), alpha = 0.1, aes(x = Latitude, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic richness", breaks = c(-5, -3, -1, 1, 3)) +
    scale_x_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    coord_cartesian(ylim = c(-4.75, 3.361625), xlim = c(36, 66)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.pd.obs.z.ele.lat, file = "p1.pd.obs.z.ele.lat_post.linpred.re2.rdata")


## Model standardized phylogenetic MPD

# Check distribution
dist <- subset(lmdata, variable == "mpd.obs.z.lat")
hist(dist$value)

dist <- subset(lmdata, variable == "mpd.obs.z.ele")
hist(dist$value)

# MPD.obs.z ~ simple elevation fixed effect model
bayes.mpd.obs.z.lat.elev.1fe <- brm(formula = value ~ Elevation.scale,
                                    data = subset(lmdata, variable == "mpd.obs.z.lat"),family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.mpd.obs.z.lat.elev.1fe, pointwise = TRUE)
loo(bayes.mpd.obs.z.lat.elev.1fe, pointwise = TRUE)

# MPD.obs.z ~ simple elevation multilevel lat effect model
bayes.mpd.obs.z.lat.elev.1re <- brm(formula = value ~ Elevation.scale + (1 + Elevation.scale | lat_class),
                                    data = subset(lmdata, variable == "mpd.obs.z.lat"),family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.mpd.obs.z.lat.elev.1re, pointwise = TRUE)
loo(bayes.mpd.obs.z.lat.elev.1re, pointwise = TRUE)

# MPD.obs.z ~ poly elevation fixed effect model
bayes.mpd.obs.z.lat.elev.2fe <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE),
                                    data = subset(lmdata, variable == "mpd.obs.z.lat"),family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.mpd.obs.z.lat.elev.2fe, pointwise = TRUE)
loo(bayes.mpd.obs.z.lat.elev.2fe, pointwise = TRUE)

# MPD.obs.z ~ poly elevation multilevel lat effect model
bayes.mpd.obs.z.lat.elev.2re <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE) + (1 + poly(Elevation.scale, 2, raw = FALSE) | lat_class),
                                    data = subset(lmdata, variable == "mpd.obs.z.lat"),family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.mpd.obs.z.lat.elev.2re, pointwise = TRUE)
loo(bayes.mpd.obs.z.lat.elev.2re, pointwise = TRUE)
summary(bayes.mpd.obs.z.lat.elev.2re)
prior_summary(bayes.mpd.obs.z.lat.elev.2re)
coef(bayes.mpd.obs.z.lat.elev.2re)
fixef(bayes.mpd.obs.z.lat.elev.2re)
bayes_R2(bayes.mpd.obs.z.lat.elev.2re)
control_params(bayes.mpd.obs.z.lat.elev.2re)

# MPD.obs.z ~ poly elevation multilevel lat effect plotting
pp_check(bayes.mpd.obs.z.lat.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.mpd.obs.z.lat.elev.2re, terms = c("Elevation.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.mpd.obs.z.lat.elev.2re, terms = c("Elevation.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 
prb$x <- prb$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 

(p1.mpd.obs.z.lat.elev <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "mpd.obs.z.lat"), alpha = 0.1, aes(x = Elevation, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic MPD", breaks = c(-5, -2.5, 0, 2.5)) +
    scale_x_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    coord_cartesian(ylim = c(-5.255725, 2.434389), xlim = c(-3.319242, 3750)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))
 
save(p1.mpd.obs.z.lat.elev, file = "p1.mpd.obs.z.lat.elev_post.linpred.re2.rdata")

# MPD.obs.z ~ simple latitude fixed effect model
bayes.mpd.obs.z.ele.lat.1fe <- brm(formula = value ~ Latitude.scale,
                                   data = subset(lmdata, variable == "mpd.obs.z.ele"),family = skew_normal(),
                                   prior = c(set_prior("normal(0,5)", class = "b")),
                                   warmup = 1000, iter = 3000, chains = 4,
                                   control = list(adapt_delta = 0.99),
                                   seed = 1)

waic(bayes.mpd.obs.z.ele.lat.1fe, pointwise = TRUE)
loo(bayes.mpd.obs.z.ele.lat.1fe, pointwise = TRUE)

# MPD.obs.z ~ simple latitude multilevel elev effect model
bayes.mpd.obs.z.ele.lat.1re <- brm(formula = value ~ Latitude.scale + (1 + Latitude.scale | elev_class),
                                   data = subset(lmdata, variable == "mpd.obs.z.ele"),family = skew_normal(),
                                   prior = c(set_prior("normal(0,5)", class = "b")),
                                   warmup = 1000, iter = 3000, chains = 4,
                                   control = list(adapt_delta = 0.99),
                                   seed = 1)

waic(bayes.mpd.obs.z.ele.lat.1re, pointwise = TRUE)
loo(bayes.mpd.obs.z.ele.lat.1re, pointwise = TRUE)

# MPD.obs.z ~ poly latitude fixed effect model
bayes.mpd.obs.z.ele.lat.2fe <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE),
                                   data = subset(lmdata, variable == "mpd.obs.z.ele"),family = skew_normal(),
                                   prior = c(set_prior("normal(0,5)", class = "b")),
                                   warmup = 1000, iter = 3000, chains = 4,
                                   control = list(adapt_delta = 0.99),
                                   seed = 1)

waic(bayes.mpd.obs.z.ele.lat.2fe, pointwise = TRUE)
loo(bayes.mpd.obs.z.ele.lat.2fe, pointwise = TRUE)

# MPD.obs.z ~ poly latitude multilevel elev effect model
bayes.mpd.obs.z.ele.lat.2re <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE) + (1 + poly(Latitude.scale, 2, raw = FALSE) | elev_class),
                                   data = subset(lmdata, variable == "mpd.obs.z.ele"),family = skew_normal(),
                                   prior = c(set_prior("normal(0,5)", class = "b")),
                                   warmup = 1000, iter = 3000, chains = 4,
                                   control = list(adapt_delta = 0.99),
                                   seed = 1)

waic(bayes.mpd.obs.z.ele.lat.2re, pointwise = TRUE)
loo(bayes.mpd.obs.z.ele.lat.2re, pointwise = TRUE)
summary(bayes.mpd.obs.z.ele.lat.2re)
prior_summary(bayes.mpd.obs.z.ele.lat.2re)
coef(bayes.mpd.obs.z.ele.lat.2re)
fixef(bayes.mpd.obs.z.ele.lat.2re)
bayes_R2(bayes.mpd.obs.z.ele.lat.2re)
control_params(bayes.mpd.obs.z.ele.lat.2re)

# MPD.obs.z ~ poly latitude multilevel elev effect plotting
pp_check(bayes.mpd.obs.z.ele.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.mpd.obs.z.ele.lat.2re, terms = c("Latitude.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.mpd.obs.z.ele.lat.2re, terms = c("Latitude.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 
prb$x <- prb$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 

(p1.mpd.obs.z.ele.lat <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "mpd.obs.z.ele"), alpha = 0.1, aes(x = Latitude, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic MPD", breaks = c(-5, -2.5, 0, 2.5)) +
    scale_x_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    coord_cartesian(ylim = c(-5.255725, 2.434389), xlim = c(36, 66)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.mpd.obs.z.ele.lat, file = "p1.mpd.obs.z.ele.lat_post.linpred.re2.rdata")


## Model standardized functional richness

# Check distribution
dist <- subset(lmdata, variable == "FRic.obs.z.lat")
hist(dist$value)

dist <- subset(lmdata, variable == "FRic.obs.z.ele")
hist(dist$value)

# FRic.obs.z ~ simple elevation fixed effect model
bayes.FRic.obs.z.lat.elev.1fe <- brm(formula = value ~ Elevation.scale,
                                     data = subset(lmdata, variable == "FRic.obs.z.lat"), family = skew_normal(),
                                     prior = c(set_prior("normal(0,5)", class = "b")),
                                     warmup = 1000, iter = 3000, chains = 4,
                                     control = list(adapt_delta = 0.99),
                                     seed = 1)

waic(bayes.FRic.obs.z.lat.elev.1fe, pointwise = TRUE)
loo(bayes.FRic.obs.z.lat.elev.1fe, pointwise = TRUE)

# FRic.obs.z ~ simple elevation multilevel lat effect model
bayes.FRic.obs.z.lat.elev.1re <- brm(formula = value ~ Elevation.scale + (1 + Elevation.scale | lat_class),
                                     data = subset(lmdata, variable == "FRic.obs.z.lat"), family = skew_normal(),
                                     prior = c(set_prior("normal(0,5)", class = "b")),
                                     warmup = 1000, iter = 3000, chains = 4,
                                     control = list(adapt_delta = 0.99),
                                     seed = 1)

waic(bayes.FRic.obs.z.lat.elev.1re, pointwise = TRUE)
loo(bayes.FRic.obs.z.lat.elev.1re, pointwise = TRUE)

# FRic.obs.z ~ poly elevation fixed effect model
bayes.FRic.obs.z.lat.elev.2fe <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE),
                                     data = subset(lmdata, variable == "FRic.obs.z.lat"), family = skew_normal(),
                                     prior = c(set_prior("normal(0,5)", class = "b")),
                                     warmup = 1000, iter = 3000, chains = 4,
                                     control = list(adapt_delta = 0.99),
                                     seed = 1)

waic(bayes.FRic.obs.z.lat.elev.2fe, pointwise = TRUE)
loo(bayes.FRic.obs.z.lat.elev.2fe, pointwise = TRUE)

# FRic.obs.z ~ poly elevation multilevel lat effect model
bayes.FRic.obs.z.lat.elev.2re <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE) + (1 + poly(Elevation.scale, 2, raw = FALSE) | lat_class),
                                     data = subset(lmdata, variable == "FRic.obs.z.lat"), family = skew_normal(),
                                     prior = c(set_prior("normal(0,5)", class = "b")),
                                     warmup = 1000, iter = 3000, chains = 4,
                                     control = list(adapt_delta = 0.99),
                                     seed = 1)

waic(bayes.FRic.obs.z.lat.elev.2re, pointwise = TRUE)
loo(bayes.FRic.obs.z.lat.elev.2re, pointwise = TRUE)
summary(bayes.FRic.obs.z.lat.elev.2re)
prior_summary(bayes.FRic.obs.z.lat.elev.2re)
coef(bayes.FRic.obs.z.lat.elev.2re)
fixef(bayes.FRic.obs.z.lat.elev.2re)
bayes_R2(bayes.FRic.obs.z.lat.elev.2re)
control_params(bayes.FRic.obs.z.lat.elev.2re)

# FRic.obs.z ~ poly elevation multilevel lat effect plotting
pp_check(bayes.FRic.obs.z.lat.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.obs.z.lat.elev.2re, terms = c("Elevation.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.obs.z.lat.elev.2re, terms = c("Elevation.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 
prb$x <- prb$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 

(p1.FRic.obs.z.lat.elev <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FRic.obs.z.lat"), alpha = 0.1, aes(x = Elevation, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional richness", breaks = c(-4, -2, 0, 2, 4)) +
    scale_x_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    coord_cartesian(ylim = c(-4.067616, 4), xlim = c(-3.319242, 3750)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values=a60505_6e068f))

save(p1.FRic.obs.z.lat.elev, file = "p1.FRic.obs.z.lat.elev_post.linpred.re2.rdata")

# FRic.obs.z ~ simple latitude fixed effect model
bayes.FRic.obs.z.ele.lat.1fe <- brm(formula = value ~ Latitude.scale,
                                    data = subset(lmdata, variable == "FRic.obs.z.ele"), family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.FRic.obs.z.ele.lat.1fe, pointwise = TRUE)
loo(bayes.FRic.obs.z.ele.lat.1fe, pointwise = TRUE)

# FRic.obs.z ~ simple latitude multilevel elev effect model
bayes.FRic.obs.z.ele.lat.1re <- brm(formula = value ~ Latitude.scale + (1 + Latitude.scale | elev_class),
                                    data = subset(lmdata, variable == "FRic.obs.z.ele"), family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.FRic.obs.z.ele.lat.1re, pointwise = TRUE)
loo(bayes.FRic.obs.z.ele.lat.1re, pointwise = TRUE)

# FRic.obs.z ~ poly latitude fixed effect model
bayes.FRic.obs.z.ele.lat.2fe <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE),
                                    data = subset(lmdata, variable == "FRic.obs.z.ele"), family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.FRic.obs.z.ele.lat.2fe, pointwise = TRUE)
loo(bayes.FRic.obs.z.ele.lat.2fe, pointwise = TRUE)

# FRic.obs.z ~ poly latitude multilevel elev effect model
bayes.FRic.obs.z.ele.lat.2re <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE) + (1 + poly(Latitude.scale, 2, raw = FALSE) | elev_class),
                                    data = subset(lmdata, variable == "FRic.obs.z.ele"), family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.FRic.obs.z.ele.lat.2re, pointwise = TRUE)
loo(bayes.FRic.obs.z.ele.lat.2re, pointwise = TRUE)
summary(bayes.FRic.obs.z.ele.lat.2re)
prior_summary(bayes.FRic.obs.z.ele.lat.2re)
coef(bayes.FRic.obs.z.ele.lat.2re)
fixef(bayes.FRic.obs.z.ele.lat.2re)
bayes_R2(bayes.FRic.obs.z.ele.lat.2re)
control_params(bayes.FRic.obs.z.ele.lat.2re)

# FRic.obs.z ~ poly latitude multilevel elev effect plotting
pp_check(bayes.FRic.obs.z.ele.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.obs.z.ele.lat.2re, terms = c("Latitude.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.obs.z.ele.lat.2re, terms = c("Latitude.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 
prb$x <- prb$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 

(p1.FRic.obs.z.ele.lat <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FRic.obs.z.ele"), alpha = 0.1, aes(x = Latitude, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional richness", breaks = c(-4, -2, 0, 2, 4)) +
    scale_x_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    coord_cartesian(ylim = c(-4.067616, 4), xlim = c(36, 66)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.FRic.obs.z.ele.lat, file = "p1.FRic.obs.z.ele.lat_post.linpred.re2.rdata")


## Model standardized functional dispersion

# Check distribution
dist <- subset(lmdata, variable == "FDis.obs.z.lat")
hist(dist$value)

dist <- subset(lmdata, variable == "FDis.obs.z.ele")
hist(dist$value)

# FDis.obs.z ~ simple elevation fixed effect model
bayes.FDis.obs.z.lat.elev.1fe <- brm(formula = value ~ Elevation.scale,
                                     data = subset(lmdata, variable == "FDis.obs.z.lat"), family = skew_normal(),
                                     prior = c(set_prior("normal(0,5)", class = "b")),
                                     warmup = 1000, iter = 3000, chains = 4,
                                     control = list(adapt_delta = 0.99),
                                     seed = 1)

waic(bayes.FDis.obs.z.lat.elev.1fe, pointwise = TRUE)
loo(bayes.FDis.obs.z.lat.elev.1fe, pointwise = TRUE)

# FDis.obs.z ~ simple elevation multilevel lat effect model
bayes.FDis.obs.z.lat.elev.1re <- brm(formula = value ~ Elevation.scale + (1 + Elevation.scale | lat_class),
                                     data = subset(lmdata, variable == "FDis.obs.z.lat"), family = skew_normal(),
                                     prior = c(set_prior("normal(0,5)", class = "b")),
                                     warmup = 1000, iter = 3000, chains = 4,
                                     control = list(adapt_delta = 0.99),
                                     seed = 1)

waic(bayes.FDis.obs.z.lat.elev.1re, pointwise = TRUE)
loo(bayes.FDis.obs.z.lat.elev.1re, pointwise = TRUE)

# FDis.obs.z ~ poly elevation fixed effect model
bayes.FDis.obs.z.lat.elev.2fe <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE),
                                     data = subset(lmdata, variable == "FDis.obs.z.lat"), family = skew_normal(),
                                     prior = c(set_prior("normal(0,5)", class = "b")),
                                     warmup = 1000, iter = 3000, chains = 4,
                                     control = list(adapt_delta = 0.99),
                                     seed = 1)

waic(bayes.FDis.obs.z.lat.elev.2fe, pointwise = TRUE)
loo(bayes.FDis.obs.z.lat.elev.2fe, pointwise = TRUE)

# FDis.obs.z ~ poly elevation multilevel lat effect model
bayes.FDis.obs.z.lat.elev.2re <- brm(formula = value ~ poly(Elevation.scale, 2, raw = FALSE) + (1 + poly(Elevation.scale, 2, raw = FALSE) | lat_class),
                                     data = subset(lmdata, variable == "FDis.obs.z.lat"), family = skew_normal(),
                                     prior = c(set_prior("normal(0,5)", class = "b")),
                                     warmup = 1000, iter = 3000, chains = 4,
                                     control = list(adapt_delta = 0.99),
                                     seed = 1)

waic(bayes.FDis.obs.z.lat.elev.2re, pointwise = TRUE)
loo(bayes.FDis.obs.z.lat.elev.2re, pointwise = TRUE)
summary(bayes.FDis.obs.z.lat.elev.2re)
prior_summary(bayes.FDis.obs.z.lat.elev.2re)
coef(bayes.FDis.obs.z.lat.elev.2re)
fixef(bayes.FDis.obs.z.lat.elev.2re)
bayes_R2(bayes.FDis.obs.z.lat.elev.2re)
control_params(bayes.FDis.obs.z.lat.elev.2re)

# FDis.obs.z ~ poly elevation multilevel lat effect plotting
pp_check(bayes.FDis.obs.z.lat.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FDis.obs.z.lat.elev.2re, terms = c("Elevation.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FDis.obs.z.lat.elev.2re, terms = c("Elevation.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 
prb$x <- prb$x * sd(lmdata$Elevation) + mean(lmdata$Elevation) 

(p1.FDis.obs.z.lat.elev <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FDis.obs.z.lat"), alpha = 0.1, aes(x = Elevation, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional dispersion", breaks = c(-3, -1.5, 0, 1.5, 3)) +
    scale_x_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    coord_cartesian(ylim = c(-3.079376, 2.948407), xlim = c(-3.319242, 3750)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.FDis.obs.z.lat.elev, file = "p1.FDis.obs.z.lat.elev_post.linpred.re2.rdata")

# FDis.obs.z ~ simple latitude fixed effect model
bayes.FDis.obs.z.ele.lat.1fe <- brm(formula = value ~ Latitude.scale,
                                    data = subset(lmdata, variable == "FDis.obs.z.ele"), family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.FDis.obs.z.ele.lat.1fe, pointwise = TRUE)
loo(bayes.FDis.obs.z.ele.lat.1fe, pointwise = TRUE)

# FDis.obs.z ~ simple latitude multilevel elev effect model
bayes.FDis.obs.z.ele.lat.1re <- brm(formula = value ~ Latitude.scale + (1 + Latitude.scale | elev_class),
                                    data = subset(lmdata, variable == "FDis.obs.z.ele"), family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.FDis.obs.z.ele.lat.1re, pointwise = TRUE)
loo(bayes.FDis.obs.z.ele.lat.1re, pointwise = TRUE)

# FDis.obs.z ~ poly latitude fixed effect model
bayes.FDis.obs.z.ele.lat.2fe <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE),
                                    data = subset(lmdata, variable == "FDis.obs.z.ele"), family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.FDis.obs.z.ele.lat.2fe, pointwise = TRUE)
loo(bayes.FDis.obs.z.ele.lat.2fe, pointwise = TRUE)

# FDis.obs.z ~ poly latitude multilevel elev effect model
bayes.FDis.obs.z.ele.lat.2re <- brm(formula = value ~ poly(Latitude.scale, 2, raw = FALSE) + (1 + poly(Latitude.scale, 2, raw = FALSE) | elev_class),
                                    data = subset(lmdata, variable == "FDis.obs.z.ele"), family = skew_normal(),
                                    prior = c(set_prior("normal(0,5)", class = "b")),
                                    warmup = 1000, iter = 3000, chains = 4,
                                    control = list(adapt_delta = 0.99),
                                    seed = 1)

waic(bayes.FDis.obs.z.ele.lat.2re, pointwise = TRUE)
loo(bayes.FDis.obs.z.ele.lat.2re, pointwise = TRUE)
summary(bayes.FDis.obs.z.ele.lat.2re)
prior_summary(bayes.FDis.obs.z.ele.lat.2re)
coef(bayes.FDis.obs.z.ele.lat.2re)
fixef(bayes.FDis.obs.z.ele.lat.2re)
bayes_R2(bayes.FDis.obs.z.ele.lat.2re)
control_params(bayes.FDis.obs.z.ele.lat.2re)

# FDis.obs.z ~ poly latitude multilevel elev effect plotting
pp_check(bayes.FDis.obs.z.ele.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FDis.obs.z.ele.lat.2re, terms = c("Latitude.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FDis.obs.z.ele.lat.2re, terms = c("Latitude.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 
prb$x <- prb$x * sd(lmdata$Latitude) + mean(lmdata$Latitude) 

(p1.FDis.obs.z.ele.lat <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FDis.obs.z.ele"), alpha = 0.1, aes(x = Latitude, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional dispersion", breaks = c(-3, -1.5, 0, 1.5, 3)) +
    scale_x_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    coord_cartesian(ylim = c(-3.079376, 2.948407), xlim = c(36, 66)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.FDis.obs.z.ele.lat, file = "p1.FDis.obs.z.ele.lat_post.linpred.re2.rdata")

# Plotting biodiversity by elevation and latitude (Figure 3 in manuscript)
# requires some minor revisions in a graphics editor
load("p1.alpha.ele_post.linpred.re2.rdata")
load("p1.alpha.lat_post.linpred.re2.rdata")
load("p1.pd.ele_post.linpred.re2.rdata")
load("p1.pd.lat_post.linpred.re2.rdata")
load("p1.mpd.obs.z.lat.elev_post.linpred.re2.rdata")
load("p1.mpd.obs.z.ele.lat_post.linpred.re2.rdata")
load("p1.FRic.ele_post.linpred.re2.rdata")
load("p1.FRic.lat_post.linpred.re2.rdata")
load("p1.FDis.obs.z.lat.elev_post.linpred.re2.rdata")
load("p1.FDis.obs.z.ele.lat_post.linpred.re2.rdata")
load("p1.CWM.Bl.ele_post.linpred.re2.rdata")
load("p1.CWM.Bl.lat_post.linpred.re2.rdata")
load("p1.pd.obs.z.ele.lat_post.linpred.re2.rdata")
load("p1.pd.obs.z.lat.elev_post.linpred.re2.rdata")
load("p1.FRic.obs.z.ele.lat_post.linpred.re2.rdata")
load("p1.FRic.obs.z.lat.elev_post.linpred.re2.rdata")

p1.alpha.ele <- ggplotGrob(p1.alpha.ele)
p1.alpha.lat <- ggplotGrob(p1.alpha.lat)
p1.pd.ele <- ggplotGrob(p1.pd.ele)
p1.pd.lat <- ggplotGrob(p1.pd.lat)
p1.mpd.obs.z.lat.elev <- ggplotGrob(p1.mpd.obs.z.lat.elev)
p1.mpd.obs.z.ele.lat <- ggplotGrob(p1.mpd.obs.z.ele.lat)
p1.FRic.ele <- ggplotGrob(p1.FRic.ele)
p1.FRic.lat <- ggplotGrob(p1.FRic.lat)
p1.FDis.obs.z.lat.elev <- ggplotGrob(p1.FDis.obs.z.lat.elev)
p1.FDis.obs.z.ele.lat <- ggplotGrob(p1.FDis.obs.z.ele.lat)
p1.CWM.Bl.ele <- ggplotGrob(p1.CWM.Bl.ele)
p1.CWM.Bl.lat <- ggplotGrob(p1.CWM.Bl.lat)
p1.pd.obs.z.lat.elev <- ggplotGrob(p1.pd.obs.z.lat.elev)
p1.pd.obs.z.ele.lat <- ggplotGrob(p1.pd.obs.z.ele.lat)
p1.FRic.obs.z.lat.elev <- ggplotGrob(p1.FRic.obs.z.lat.elev)
p1.FRic.obs.z.ele.lat <- ggplotGrob(p1.FRic.obs.z.ele.lat)

g2.c1 <- rbind(p1.alpha.ele,
               p1.FRic.ele,
               p1.pd.ele,
               p1.FDis.obs.z.lat.elev)

g2.c1$widths <- unit.pmax(p1.alpha.ele$widths, 
                          p1.FRic.ele$widths,
                          p1.pd.ele$widths,
                          p1.FDis.obs.z.lat.elev$widths)

g2.c2 <- rbind(p1.CWM.Bl.ele,
               p1.FRic.obs.z.lat.elev,
               p1.pd.obs.z.lat.elev,
               p1.mpd.obs.z.lat.elev)

g2.c2$widths <- unit.pmax(p1.CWM.Bl.ele$widths, 
                          p1.FRic.obs.z.lat.elev$widths,
                          p1.pd.obs.z.lat.elev$widths,
                          p1.mpd.obs.z.lat.elev$widths)

g2.c3 <- rbind(p1.alpha.lat,
               p1.FRic.lat,
               p1.pd.lat,
               p1.FDis.obs.z.ele.lat)

g2.c3$widths <- unit.pmax(p1.alpha.lat$widths, 
                          p1.FRic.lat$widths,
                          p1.pd.lat$widths,
                          p1.FDis.obs.z.ele.lat$widths)

g2.c4 <- rbind(p1.CWM.Bl.lat,
               p1.FRic.obs.z.ele.lat,
               p1.pd.obs.z.ele.lat,
               p1.mpd.obs.z.ele.lat)

g2.c4$widths <- unit.pmax(p1.CWM.Bl.lat$widths, 
                          p1.FRic.obs.z.ele.lat$widths,
                          p1.pd.obs.z.ele.lat$widths,
                          p1.mpd.obs.z.ele.lat$widths)

g2 <- cbind(g2.c1, g2.c2, g2.c3, g2.c4)
grid.newpage()
grid.draw(g2)


## Model mean annual temperatures relationships to species richness, elevation, and latitude

# Check distribution
hist(lmdata$MAT)
hist(lmdata$Elevation.scale)
hist(lmdata$Latitude.scale)

# Elevation ~ poly MAT multilevel lat effect model
bayes.ELEV.MAT.lat.2re <- brm(formula = Elevation.scale ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | lat_class),
                              data = subset(lmdata, variable == "Taxo.alpha"), family = gaussian(),
                              prior = c(set_prior("normal(0,5)", class = "b")),
                              warmup = 1000, iter = 3000, chains = 4,
                              control = list(adapt_delta = 0.99, max_treedepth = 15),
                              seed = 1)

waic(bayes.ELEV.MAT.lat.2re, pointwise = TRUE)
loo(bayes.ELEV.MAT.lat.2re, pointwise = TRUE)
summary(bayes.ELEV.MAT.lat.2re)
prior_summary(bayes.ELEV.MAT.lat.2re)
coef(bayes.ELEV.MAT.lat.2re)
fixef(bayes.ELEV.MAT.lat.2re)
bayes_R2(bayes.ELEV.MAT.lat.2re)
control_params(bayes.ELEV.MAT.lat.2re)

# Elevation ~ poly MAT multilevel lat effect plotting
pp_check(bayes.ELEV.MAT.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.ELEV.MAT.lat.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.ELEV.MAT.lat.2re, terms = c("MAT.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

pra$predicted <- pra$predicted * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)
prb$predicted <- prb$predicted * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)

pra$conf.low <- pra$conf.low * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)
prb$conf.low <- prb$conf.low * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)

pra$conf.high <- pra$conf.high * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)
prb$conf.high <- prb$conf.high * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)

(p1.ELEV.MAT.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = MAT, y = Elevation, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(-3.319242, 3750), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.ELEV.MAT.lat, file = "p1.ELEV.MAT.lat_post.linpred.re2.rdata")

# Latitude ~ poly MAT multilevel elev effect model
bayes.LAT.MAT.elev.2re <- brm(formula = Latitude.scale ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | elev_class),
                              data = subset(lmdata, variable == "Taxo.alpha"), family = gaussian(),
                              prior = c(set_prior("normal(0,5)", class = "b")),
                              warmup = 1000, iter = 3000, chains = 4,
                              control = list(adapt_delta = 0.99, max_treedepth = 15),
                              seed = 1)

waic(bayes.LAT.MAT.elev.2re, pointwise = TRUE)
loo(bayes.LAT.MAT.elev.2re, pointwise = TRUE)
summary(bayes.LAT.MAT.elev.2re)
prior_summary(bayes.LAT.MAT.elev.2re)
coef(bayes.LAT.MAT.elev.2re)
fixef(bayes.LAT.MAT.elev.2re)
bayes_R2(bayes.LAT.MAT.elev.2re)
control_params(bayes.LAT.MAT.elev.2re)

# Latitude ~ poly MAT multilevel elev effect plotting
pp_check(bayes.LAT.MAT.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.LAT.MAT.elev.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.LAT.MAT.elev.2re, terms = c("MAT.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

pra$predicted <- pra$predicted * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)
prb$predicted <- prb$predicted * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)

pra$conf.low <- pra$conf.low * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)
prb$conf.low <- prb$conf.low * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)

pra$conf.high <- pra$conf.high * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)
prb$conf.high <- prb$conf.high * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)

(p1.LAT.MAT.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = MAT, y = Latitude, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(36, 66.2), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.LAT.MAT.elev, file = "p1.LAT.MAT.elev_post.linpred.re2.rdata")

# SR ~ poly MAT fixed effect model
bayes.alpha.MAT.2fe <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE)),
                           data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                           prior = c(set_prior("normal(0, 5)", class = "b")),
                           warmup = 1000, iter = 3000, chains = 4,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.alpha.MAT.2fe, pointwise = TRUE)
loo(bayes.alpha.MAT.2fe, pointwise = TRUE)
summary(bayes.alpha.MAT.2fe)

# SR ~ poly MAT fixed effect plotting
pra <- ggpredict(bayes.alpha.MAT.2fe, terms = c("MAT.scale [all]"), ppd = FALSE)

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.alpha.MAT.2fe <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = MAT, y = value, colour = "red")) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Species richness", breaks = c(0, 5, 10, 15, 20, 25)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(0, 26), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    theme(legend.position="none") +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.alpha.MAT.2fe, file = "p1.alpha.MAT_post.linpred.2fe.rdata")

# SR ~ poly MAT multilevel lat effect model
bayes.alpha.MAT.lat.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | lat_class),
                               data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                               prior = c(set_prior("normal(0,5)", class = "b")),
                               warmup = 1000, iter = 3000, chains = 4,
                               control = list(adapt_delta = 0.99),
                               seed = 1)

waic(bayes.alpha.MAT.lat.2re, pointwise = TRUE)
loo(bayes.alpha.MAT.lat.2re, pointwise = TRUE)
summary(bayes.alpha.MAT.lat.2re)
prior_summary(bayes.alpha.MAT.lat.2re)
coef(bayes.alpha.MAT.lat.2re)
fixef(bayes.alpha.MAT.lat.2re)
bayes_R2(bayes.alpha.MAT.lat.2re)
control_params(bayes.alpha.MAT.lat.2re)

# SR ~ poly MAT multilevel lat effect plotting
pp_check(bayes.alpha.MAT.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.alpha.MAT.lat.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.alpha.MAT.lat.2re, terms = c("MAT.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.alpha.MAT.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = MAT, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Species richness", breaks = c(0, 5, 10, 15, 20, 25)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(0, 26), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.alpha.MAT.lat, file = "p1.alpha.MAT.lat_post.linpred.re2.rdata")

# SR ~ poly MAT multilevel elev effect model
bayes.alpha.MAT.elev.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | elev_class),
                                data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                                prior = c(set_prior("normal(0,5)", class = "b")),
                                warmup = 1000, iter = 3000, chains = 4,
                                control = list(adapt_delta = 0.99),
                                seed = 1)

waic(bayes.alpha.MAT.elev.2re, pointwise = TRUE)
loo(bayes.alpha.MAT.elev.2re, pointwise = TRUE)
summary(bayes.alpha.MAT.elev.2re)
prior_summary(bayes.alpha.MAT.elev.2re)
coef(bayes.alpha.MAT.elev.2re)
fixef(bayes.alpha.MAT.elev.2re)
bayes_R2(bayes.alpha.MAT.elev.2re)
control_params(bayes.alpha.MAT.elev.2re)

# SR ~ poly MAT multilevel elev effect plotting
pp_check(bayes.alpha.MAT.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.alpha.MAT.elev.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.alpha.MAT.elev.2re, terms = c("MAT.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.alpha.MAT.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = MAT, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Species richness", breaks = c(0, 5, 10, 15, 20, 25)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(0, 26), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.alpha.MAT.elev, file = "p1.alpha.MAT.elev_post.linpred.re2.rdata")

                 
## Model mean temperature difference relationships to species richness, elevation, and latitude

# Check distribution
hist(lmdata$TD)
hist(lmdata$Elevation.scale)
hist(lmdata$Latitude.scale)

# Elevation ~ poly TD multilevel lat effect model
bayes.ELEV.TD.lat.2re <- brm(formula = Elevation.scale ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | lat_class),
                             data = subset(lmdata, variable == "Taxo.alpha"), family = gaussian(),
                             prior = c(set_prior("normal(0,5)", class = "b")),
                             warmup = 1000, iter = 3000, chains = 4,
                             control = list(adapt_delta = 0.99, max_treedepth = 15),
                             seed = 1)

waic(bayes.ELEV.TD.lat.2re, pointwise = TRUE)
loo(bayes.ELEV.TD.lat.2re, pointwise = TRUE)
summary(bayes.ELEV.TD.lat.2re)
prior_summary(bayes.ELEV.TD.lat.2re)
coef(bayes.ELEV.TD.lat.2re)
fixef(bayes.ELEV.TD.lat.2re)
bayes_R2(bayes.ELEV.TD.lat.2re)
control_params(bayes.ELEV.TD.lat.2re)

# Elevation ~ poly TD multilevel lat effect plotting
pp_check(bayes.ELEV.TD.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.ELEV.TD.lat.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.ELEV.TD.lat.2re, terms = c("TD.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

pra$predicted <- pra$predicted * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)
prb$predicted <- prb$predicted * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)

pra$conf.low <- pra$conf.low * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)
prb$conf.low <- prb$conf.low * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)

pra$conf.high <- pra$conf.high * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)
prb$conf.high <- prb$conf.high * sd(subset(lmdata, variable == "Taxo.alpha")$Elevation) + mean(subset(lmdata, variable == "Taxo.alpha")$Elevation)

(p1.ELEV.TD.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = TD, y = Elevation, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Elevation (masl)", breaks = c(0, 750, 1500, 2250, 3000, 3750)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(-3.319242, 3750), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.ELEV.TD.lat, file = "p1.ELEV.TD.lat_post.linpred.re2.rdata")

# Latitude ~ poly TD multilevel elev effect model
bayes.LAT.TD.elev.2re <- brm(formula = Latitude.scale ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | elev_class),
                             data = subset(lmdata, variable == "Taxo.alpha"), family = gaussian(),
                             prior = c(set_prior("normal(0,5)", class = "b")),
                             warmup = 1000, iter = 3000, chains = 4,
                             control = list(adapt_delta = 0.99, max_treedepth = 15),
                             seed = 1)

waic(bayes.LAT.TD.elev.2re, pointwise = TRUE)
loo(bayes.LAT.TD.elev.2re, pointwise = TRUE)
summary(bayes.LAT.TD.elev.2re)
prior_summary(bayes.LAT.TD.elev.2re)
coef(bayes.LAT.TD.elev.2re)
fixef(bayes.LAT.TD.elev.2re)
bayes_R2(bayes.LAT.TD.elev.2re)
control_params(bayes.LAT.TD.elev.2re)

# Latitude ~ poly TD multilevel elev effect plotting
pp_check(bayes.LAT.TD.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.LAT.TD.elev.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.LAT.TD.elev.2re, terms = c("TD.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

pra$predicted <- pra$predicted * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)
prb$predicted <- prb$predicted * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)

pra$conf.low <- pra$conf.low * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)
prb$conf.low <- prb$conf.low * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)

pra$conf.high <- pra$conf.high * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)
prb$conf.high <- prb$conf.high * sd(subset(lmdata, variable == "Taxo.alpha")$Latitude) + mean(subset(lmdata, variable == "Taxo.alpha")$Latitude)

(p1.LAT.TD.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = TD, y = Latitude, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Latitude (°N)", breaks = c(36, 42, 48, 54, 60, 66)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(36, 66.2), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.LAT.TD.elev, file = "p1.LAT.TD.elev_post.linpred.re2.rdata")

# SR ~ poly TD fixed effect model
bayes.alpha.TD.2fe <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE)),
                          data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                          prior = c(set_prior("normal(0, 5)", class = "b")),
                          warmup = 1000, iter = 3000, chains = 4,
                          control = list(adapt_delta = 0.99),
                          seed = 1)

waic(bayes.alpha.TD.2fe, pointwise = TRUE)
loo(bayes.alpha.TD.2fe, pointwise = TRUE)
summary(bayes.alpha.TD.2fe)

# SR ~ poly TD fixed effect plotting
pra <- ggpredict(bayes.alpha.TD.2fe, terms = c("TD.scale [all]"), ppd = FALSE)

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.alpha.TD.2fe <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = TD, y = value, colour = "red")) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Species richness", breaks = c(0, 5, 10, 15, 20, 25)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(0, 26), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    theme(legend.position="none") +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.alpha.TD.2fe, file = "p1.alpha.TD_post.linpred.2fe.rdata")

# SR ~ poly TD multilevel lat effect model
bayes.alpha.TD.lat.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | lat_class),
                              data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                              prior = c(set_prior("normal(0,5)", class = "b")),
                              warmup = 1000, iter = 3000, chains = 4,
                              control = list(adapt_delta = 0.99, max_treedepth = 15),
                              seed = 1)

waic(bayes.alpha.TD.lat.2re, pointwise = TRUE)
loo(bayes.alpha.TD.lat.2re, pointwise = TRUE)
summary(bayes.alpha.TD.lat.2re)
prior_summary(bayes.alpha.TD.lat.2re)
coef(bayes.alpha.TD.lat.2re)
fixef(bayes.alpha.TD.lat.2re)
bayes_R2(bayes.alpha.TD.lat.2re)
control_params(bayes.alpha.TD.lat.2re)

# SR ~ poly TD multilevel lat effect plotting
pp_check(bayes.alpha.TD.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.alpha.TD.lat.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.alpha.TD.lat.2re, terms = c("TD.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.alpha.TD.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = TD, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Species richness", breaks = c(0, 5, 10, 15, 20, 25)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(0, 26), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.alpha.TD.lat, file = "p1.alpha.TD.lat_post.linpred.re2.rdata")

# SR ~ poly TD multilevel elev effect model
bayes.alpha.TD.elev.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | elev_class),
                               data = subset(lmdata, variable == "Taxo.alpha"), family = negbinomial(),
                               prior = c(set_prior("normal(0,5)", class = "b")),
                               warmup = 1000, iter = 3000, chains = 4,
                               control = list(adapt_delta = 0.99),
                               seed = 1)

waic(bayes.alpha.TD.elev.2re, pointwise = TRUE)
loo(bayes.alpha.TD.elev.2re, pointwise = TRUE)
summary(bayes.alpha.TD.elev.2re)
prior_summary(bayes.alpha.TD.elev.2re)
coef(bayes.alpha.TD.elev.2re)
fixef(bayes.alpha.TD.elev.2re)
bayes_R2(bayes.alpha.TD.elev.2re)
control_params(bayes.alpha.TD.elev.2re)

# SP ~ poly TD multilevel elev effect plotting
pp_check(bayes.alpha.TD.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.alpha.TD.elev.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.alpha.TD.elev.2re, terms = c("TD.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.alpha.TD.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Taxo.alpha"), alpha = 0.1, aes(x = TD, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Species richness", breaks = c(0, 5, 10, 15, 20, 25)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(0, 26), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.alpha.TD.elev, file = "p1.alpha.TD.elev_post.linpred.re2.rdata")

# Plotting climate by elevation, latitude, and species richness (Figure 4 in manuscript)
# requires some minor revisions in a graphics editor
load("p1.ELEV.MAT.lat_post.linpred.re2.rdata")
load("p1.LAT.MAT.elev_post.linpred.re2.rdata")
load("p1.ELEV.TD.lat_post.linpred.re2.rdata")
load("p1.LAT.TD.elev_post.linpred.re2.rdata")
load("p1.alpha.MAT.elev_post.linpred.re2.rdata")
load("p1.alpha.MAT.lat_post.linpred.re2.rdata")
load("p1.alpha.TD.elev_post.linpred.re2.rdata")
load("p1.alpha.TD.lat_post.linpred.re2.rdata")

p1.ELEV.MAT.lat <- ggplotGrob(p1.ELEV.MAT.lat)
p1.LAT.MAT.elev <- ggplotGrob(p1.LAT.MAT.elev)
p1.ELEV.TD.lat <- ggplotGrob(p1.ELEV.TD.lat)
p1.LAT.TD.elev <- ggplotGrob(p1.LAT.TD.elev)
p1.alpha.MAT.elev <- ggplotGrob(p1.alpha.MAT.elev)
p1.alpha.MAT.lat <- ggplotGrob(p1.alpha.MAT.lat)
p1.alpha.TD.elev <- ggplotGrob(p1.alpha.TD.elev)
p1.alpha.TD.lat <- ggplotGrob(p1.alpha.TD.lat)

g3.c1 <- rbind(p1.ELEV.MAT.lat,
               p1.alpha.MAT.lat,
               p1.ELEV.TD.lat,
               p1.alpha.TD.lat)

g3.c1$widths <- unit.pmax(p1.ELEV.MAT.lat$widths,
                          p1.alpha.MAT.lat$widths, 
                          p1.ELEV.TD.lat$widths,
                          p1.alpha.TD.lat$widths)

g3.c2 <- rbind(p1.LAT.MAT.elev,
               p1.alpha.MAT.elev,
               p1.LAT.TD.elev,
               p1.alpha.TD.elev)

g3.c2$widths <- unit.pmax(p1.LAT.MAT.elev$widths,
                          p1.alpha.MAT.elev$widths, 
                          p1.LAT.TD.elev$widths,
                          p1.alpha.TD.elev$widths)

g3 <- cbind(g3.c1, g3.c2)
grid.newpage()
grid.draw(g3)

# Plotting species richness by polynomial elevation, latitude, MAT, and TD fixed effect plots (SI in manuscript)
# requires some minor revisions in a graphics editor
load("p1.alpha.ele_post.linpred.2fe.rdata")
load("p1.alpha.lat_post.linpred.2fe.rdata")
load("p1.alpha.MAT_post.linpred.2fe.rdata")
load("p1.alpha.TD_post.linpred.2fe.rdata")

p1.alpha.ele.2fe <- ggplotGrob(p1.alpha.ele.2fe)
p1.alpha.lat.2fe <- ggplotGrob(p1.alpha.lat.2fe)
p1.alpha.MAT.2fe <- ggplotGrob(p1.alpha.MAT.2fe)
p1.alpha.TD.2fe <- ggplotGrob(p1.alpha.TD.2fe)

g4.c1 <- rbind(p1.alpha.ele.2fe,
               p1.alpha.MAT.2fe)

g4.c1$widths <- unit.pmax(p1.alpha.ele.2fe$widths,
                          p1.alpha.MAT.2fe$widths)

g4.c2 <- rbind(p1.alpha.lat.2fe,
             p1.alpha.TD.2fe)

g4.c2$widths <- unit.pmax(p1.alpha.lat.2fe$widths,
                          p1.alpha.TD.2fe$widths)

g4 <- cbind(g4.c1, g4.c2)
grid.newpage()
grid.draw(g4)


##################
# Plotting coefficients from Bayesian generalized linear multilevel models
##################

# Load data for coefficient plots (obtained from prior analyses)
coef = read.table("coefficients.bio.metrics.csv", header = T, sep = ",") #use for geographic gradient analyses
coef = read.table("coefficients.clim.csv", header = T, sep = ",") #use for temperature gradient analyses

# Prepare data for plotting
coef.elev <- subset(coef, Group == "Elevation")
coef.lat <- subset(coef, Group == "Latitude")

coef.elev$Class <- as.ordered(coef.elev$Class)

coef.elev$Class <- factor(coef.elev$Class, levels = c("36.6-39.3",
                                                      "39.3-41.2",
                                                      "41.2-43.2",
                                                      "43.2-44.6",
                                                      "44.6-46.3",
                                                      "46.3-48.0",
                                                      "48.0-50.0",
                                                      "50.0-52.1",
                                                      "52.1-54.2",
                                                      "54.2-58.1",
                                                      "58.1-62.1",
                                                      "62.1-66.2"))

coef.lat$Class <- as.ordered(coef.lat$Class)

coef.lat$Class <- factor(coef.lat$Class, levels = c("0-258",
                                                    "258-599",
                                                    "599-856",
                                                    "856-1100",
                                                    "1100-1320",
                                                    "1320-1500",
                                                    "1500-1670",
                                                    "1670-1870",
                                                    "1870-2090",
                                                    "2090-2330",
                                                    "2330-2780",
                                                    "2780-3740"))


## Plot biodiversity model intercepts
(coef.int.sr.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Species richness"), 
                           aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.sr.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Species richness"), 
                            aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.cwm.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Community mean body length"), 
                            aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Community mean body length") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.cwm.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Community mean body length"), 
                             aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Community mean body length") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.fricmu1.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Functional richness mu1"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu1") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.fricmu1.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Functional richness mu1"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu1") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.fricmu2.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Functional richness mu2"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu2") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values=a60505_6e068f))

(coef.int.fricmu2.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Functional richness mu2"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu2") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.fric.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Standardized functional richness"), 
                                   aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.fric.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Standardized functional richness"), 
                                    aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.pd.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Phylogenetic richness"), 
                           aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Phylogenetic richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.pd.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Phylogenetic richness"), 
                            aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Phylogenetic richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.pd.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Standardized phylogenetic richness"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.pd.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Standardized phylogenetic richness"), 
                                  aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.fdis.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Standardized functional dispersion"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional dispersion") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.fdis.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Standardized functional dispersion"), 
                                  aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional dispersion") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.mpd.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Standardized phylogenetic MPD"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic MPD") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.mpd.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Standardized phylogenetic MPD"), 
                                  aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic MPD") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

# Combine plots of biodiversity model intercepts (SI in manuscript)
# requires some minor revisions in a graphics editor
coef.int.sr.lat <- ggplotGrob(coef.int.sr.lat)
coef.int.sr.elev <- ggplotGrob(coef.int.sr.elev)
coef.int.cwm.lat <- ggplotGrob(coef.int.cwm.lat)
coef.int.cwm.elev <- ggplotGrob(coef.int.cwm.elev)
coef.int.fricmu1.lat <- ggplotGrob(coef.int.fricmu1.lat)
coef.int.fricmu1.elev <- ggplotGrob(coef.int.fricmu1.elev)
coef.int.fricmu2.lat <- ggplotGrob(coef.int.fricmu2.lat)
coef.int.fricmu2.elev <- ggplotGrob(coef.int.fricmu2.elev)
coef.int.fric.obs.z.lat <- ggplotGrob(coef.int.fric.obs.z.lat)
coef.int.fric.obs.z.elev <- ggplotGrob(coef.int.fric.obs.z.elev)
coef.int.pd.lat <- ggplotGrob(coef.int.pd.lat)
coef.int.pd.elev <- ggplotGrob(coef.int.pd.elev)
coef.int.pd.obs.z.lat <- ggplotGrob(coef.int.pd.obs.z.lat)
coef.int.pd.obs.z.elev <- ggplotGrob(coef.int.pd.obs.z.elev)
coef.int.fdis.obs.z.lat <- ggplotGrob(coef.int.fdis.obs.z.lat)
coef.int.fdis.obs.z.elev <- ggplotGrob(coef.int.fdis.obs.z.elev)
coef.int.mpd.obs.z.lat <- ggplotGrob(coef.int.mpd.obs.z.lat)
coef.int.mpd.obs.z.elev <- ggplotGrob(coef.int.mpd.obs.z.elev)

g5.c1 <- rbind(coef.int.sr.lat,
               coef.int.fricmu1.lat,
               coef.int.fricmu2.lat,
               coef.int.pd.lat,
               coef.int.fdis.obs.z.lat)

g5.c1$widths <- unit.pmax(coef.int.sr.lat$widths, 
                          coef.int.fricmu1.lat$widths,
                          coef.int.fricmu2.lat$widths,
                          coef.int.pd.lat$widths,
                          coef.int.fdis.obs.z.lat$widths)

g5.c2 <- rbind(coef.int.cwm.lat,
               coef.int.fric.obs.z.lat,
               coef.int.fricmu2.lat,
               coef.int.pd.obs.z.lat,
               coef.int.mpd.obs.z.lat)

g5.c2$widths <- unit.pmax(coef.int.cwm.lat$widths, 
                          coef.int.fric.obs.z.lat$widths,
                          coef.int.fricmu2.lat$widths,
                          coef.int.pd.obs.z.lat$widths,
                          coef.int.mpd.obs.z.lat$widths)

g5.c3 <- rbind(coef.int.sr.elev,
               coef.int.fricmu1.elev,
               coef.int.fricmu2.elev,
               coef.int.pd.elev,
               coef.int.fdis.obs.z.elev)

g5.c3$widths <- unit.pmax(coef.int.sr.elev$widths, 
                          coef.int.fricmu1.elev$widths,
                          coef.int.fricmu2.elev$widths,
                          coef.int.pd.elev$widths,
                          coef.int.fdis.obs.z.elev$widths)

g5.c4 <- rbind(coef.int.cwm.elev,
               coef.int.fric.obs.z.elev,
               coef.int.fricmu2.elev,
               coef.int.pd.obs.z.elev,
               coef.int.mpd.obs.z.elev)

g5.c4$widths <- unit.pmax(coef.int.cwm.elev$widths, 
                          coef.int.fric.obs.z.elev$widths,
                          coef.int.fricmu2.elev$widths,
                          coef.int.pd.obs.z.elev$widths,
                          coef.int.mpd.obs.z.elev$widths)

g5<-cbind(g5.c1, g5.c2, g5.c3, g5.c4)
grid.newpage()
grid.draw(g5)


## Plot biodiversity model B1 slopes
(coef.B1.sr.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Species richness"), 
                          aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.sr.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Species richness"), 
                           aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.cwm.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Community mean body length"), 
                           aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Community mean body length") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.cwm.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Community mean body length"), 
                            aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Community mean body length") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.fricmu1.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Functional richness mu1"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu1") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.fricmu1.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Functional richness mu1"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu1") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.fricmu2.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Functional richness mu2"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu2") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.fricmu2.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Functional richness mu2"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu2") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.fric.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Standardized functional richness"), 
                                  aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept=0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.fric.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Standardized functional richness"), 
                                   aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.pd.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Phylogenetic richness"), 
                          aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Phylogenetic richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.pd.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Phylogenetic richness"), 
                           aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Phylogenetic richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.pd.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Standardized phylogenetic richness"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.pd.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Standardized phylogenetic richness"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.fdis.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Standardized functional dispersion"),
                                  aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional dispersion") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.fdis.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Standardized functional dispersion"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional dispersion") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.mpd.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Standardized phylogenetic MPD"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic MPD") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.mpd.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Standardized phylogenetic MPD"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic MPD") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

# Combine plots of biodiversity model B1 slopes (SI in manuscript)
# requires some minor revisions in a graphics editor
coef.B1.sr.lat <- ggplotGrob(coef.B1.sr.lat)
coef.B1.sr.elev <- ggplotGrob(coef.B1.sr.elev)
coef.B1.cwm.lat <- ggplotGrob(coef.B1.cwm.lat)
coef.B1.cwm.elev <- ggplotGrob(coef.B1.cwm.elev)
coef.B1.fricmu1.lat <- ggplotGrob(coef.B1.fricmu1.lat)
coef.B1.fricmu1.elev <- ggplotGrob(coef.B1.fricmu1.elev)
coef.B1.fricmu2.lat <- ggplotGrob(coef.B1.fricmu2.lat)
coef.B1.fricmu2.elev <- ggplotGrob(coef.B1.fricmu2.elev)
coef.B1.fric.obs.z.lat <- ggplotGrob(coef.B1.fric.obs.z.lat)
coef.B1.fric.obs.z.elev <- ggplotGrob(coef.B1.fric.obs.z.elev)
coef.B1.pd.lat <- ggplotGrob(coef.B1.pd.lat)
coef.B1.pd.elev <- ggplotGrob(coef.B1.pd.elev)
coef.B1.pd.obs.z.lat <- ggplotGrob(coef.B1.pd.obs.z.lat)
coef.B1.pd.obs.z.elev <- ggplotGrob(coef.B1.pd.obs.z.elev)
coef.B1.fdis.obs.z.lat <- ggplotGrob(coef.B1.fdis.obs.z.lat)
coef.B1.fdis.obs.z.elev <- ggplotGrob(coef.B1.fdis.obs.z.elev)
coef.B1.mpd.obs.z.lat <- ggplotGrob(coef.B1.mpd.obs.z.lat)
coef.B1.mpd.obs.z.elev <- ggplotGrob(coef.B1.mpd.obs.z.elev)

g6.c1 <- rbind(coef.B1.sr.lat,
               coef.B1.fricmu1.lat,
               coef.B1.fricmu2.lat,
               coef.B1.pd.lat,
               coef.B1.fdis.obs.z.lat)

g6.c1$widths <- unit.pmax(coef.B1.sr.lat$widths, 
                          coef.B1.fricmu1.lat$widths,
                          coef.B1.fricmu2.lat$widths,
                          coef.B1.pd.lat$widths,
                          coef.B1.fdis.obs.z.lat$widths)

g6.c2 <- rbind(coef.B1.cwm.lat,
               coef.B1.fric.obs.z.lat,
               coef.B1.fricmu2.lat,
               coef.B1.pd.obs.z.lat,
               coef.B1.mpd.obs.z.lat)

g6.c2$widths <- unit.pmax(coef.B1.cwm.lat$widths, 
                          coef.B1.fric.obs.z.lat$widths,
                          coef.B1.fricmu2.lat$widths,
                          coef.B1.pd.obs.z.lat$widths,
                          coef.B1.mpd.obs.z.lat$widths)

g6.c3 <- rbind(coef.B1.sr.elev,
               coef.B1.fricmu1.elev,
               coef.B1.fricmu2.elev,
               coef.B1.pd.elev,
               coef.B1.fdis.obs.z.elev)

g6.c3$widths <- unit.pmax(coef.B1.sr.elev$widths, 
                          coef.B1.fricmu1.elev$widths,
                          coef.B1.fricmu2.elev$widths,
                          coef.B1.pd.elev$widths,
                          coef.B1.fdis.obs.z.elev$widths)

g6.c4 <- rbind(coef.B1.cwm.elev,
               coef.B1.fric.obs.z.elev,
               coef.B1.fricmu2.elev,
               coef.B1.pd.obs.z.elev,
               coef.B1.mpd.obs.z.elev)

g6.c4$widths <- unit.pmax(coef.B1.cwm.elev$widths, 
                          coef.B1.fric.obs.z.elev$widths,
                          coef.B1.fricmu2.elev$widths,
                          coef.B1.pd.obs.z.elev$widths,
                          coef.B1.mpd.obs.z.elev$widths)

g6 <- cbind(g6.c1, g6.c2, g6.c3, g6.c4)
grid.newpage()
grid.draw(g6)


## Plot biodiversity model B2 curvatures
(coef.B2.sr.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Species richness"), 
                          aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.sr.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Species richness"), 
                           aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.cwm.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Community mean body length"), 
                           aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Community mean body length") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.cwm.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Community mean body length"), 
                            aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Community mean body length") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.fricmu1.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Functional richness mu1"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu1") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.fricmu1.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Functional richness mu1"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu1") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.fricmu2.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Functional richness mu2"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu2") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.fricmu2.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Functional richness mu2"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Functional richness mu2") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.fric.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Standardized functional richness"), 
                                  aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.fric.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Standardized functional richness"), 
                                   aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.pd.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Phylogenetic richness"), 
                          aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Phylogenetic richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.pd.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Phylogenetic richness"), 
                           aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Phylogenetic richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.pd.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Standardized phylogenetic richness"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.pd.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Standardized phylogenetic richness"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.fdis.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Standardized functional dispersion"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional dispersion") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.fdis.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Standardized functional dispersion"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized functional dispersion") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.mpd.obs.z.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Standardized phylogenetic MPD"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic MPD") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.mpd.obs.z.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Standardized phylogenetic MPD"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Standardized phylogenetic MPD") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

# Combine plots of biodiversity model B2 curvatures (SI in manuscript)
# requires some minor revisions in a graphics editor
coef.B2.sr.lat <- ggplotGrob(coef.B2.sr.lat)
coef.B2.sr.elev <- ggplotGrob(coef.B2.sr.elev)
coef.B2.cwm.lat <- ggplotGrob(coef.B2.cwm.lat)
coef.B2.cwm.elev <- ggplotGrob(coef.B2.cwm.elev)
coef.B2.fricmu1.lat <- ggplotGrob(coef.B2.fricmu1.lat)
coef.B2.fricmu1.elev <- ggplotGrob(coef.B2.fricmu1.elev)
coef.B2.fricmu2.lat <- ggplotGrob(coef.B2.fricmu2.lat)
coef.B2.fricmu2.elev <- ggplotGrob(coef.B2.fricmu2.elev)
coef.B2.fric.obs.z.lat <- ggplotGrob(coef.B2.fric.obs.z.lat)
coef.B2.fric.obs.z.elev <- ggplotGrob(coef.B2.fric.obs.z.elev)
coef.B2.pd.lat <- ggplotGrob(coef.B2.pd.lat)
coef.B2.pd.elev <- ggplotGrob(coef.B2.pd.elev)
coef.B2.pd.obs.z.lat <- ggplotGrob(coef.B2.pd.obs.z.lat)
coef.B2.pd.obs.z.elev <- ggplotGrob(coef.B2.pd.obs.z.elev)
coef.B2.fdis.obs.z.lat <- ggplotGrob(coef.B2.fdis.obs.z.lat)
coef.B2.fdis.obs.z.elev <- ggplotGrob(coef.B2.fdis.obs.z.elev)
coef.B2.mpd.obs.z.lat <- ggplotGrob(coef.B2.mpd.obs.z.lat)
coef.B2.mpd.obs.z.elev <- ggplotGrob(coef.B2.mpd.obs.z.elev)

g7.c1 <- rbind(coef.B2.sr.lat,
               coef.B2.fricmu1.lat,
               coef.B2.fricmu2.lat,
               coef.B2.pd.lat,
               coef.B2.fdis.obs.z.lat)

g7.c1$widths <- unit.pmax(coef.B2.sr.lat$widths, 
                          coef.B2.fricmu1.lat$widths,
                          coef.B2.fricmu2.lat$widths,
                          coef.B2.pd.lat$widths,
                          coef.B2.fdis.obs.z.lat$widths)

g7.c2 <- rbind(coef.B2.cwm.lat,
               coef.B2.fric.obs.z.lat,
               coef.B2.fricmu2.lat,
               coef.B2.pd.obs.z.lat,
               coef.B2.mpd.obs.z.lat)

g7.c2$widths <- unit.pmax(coef.B2.cwm.lat$widths, 
                          coef.B2.fric.obs.z.lat$widths,
                          coef.B2.fricmu2.lat$widths,
                          coef.B2.pd.obs.z.lat$widths,
                          coef.B2.mpd.obs.z.lat$widths)

g7.c3 <- rbind(coef.B2.sr.elev,
               coef.B2.fricmu1.elev,
               coef.B2.fricmu2.elev,
               coef.B2.pd.elev,
               coef.B2.fdis.obs.z.elev)

g7.c3$widths <- unit.pmax(coef.B2.sr.elev$widths, 
                          coef.B2.fricmu1.elev$widths,
                          coef.B2.fricmu2.elev$widths,
                          coef.B2.pd.elev$widths,
                          coef.B2.fdis.obs.z.elev$widths)

g7.c4 <- rbind(coef.B2.cwm.elev,
               coef.B2.fric.obs.z.elev,
               coef.B2.fricmu2.elev,
               coef.B2.pd.obs.z.elev,
               coef.B2.mpd.obs.z.elev)

g7.c4$widths <- unit.pmax(coef.B2.cwm.elev$widths, 
                          coef.B2.fric.obs.z.elev$widths,
                          coef.B2.fricmu2.elev$widths,
                          coef.B2.pd.obs.z.elev$widths,
                          coef.B2.mpd.obs.z.elev$widths)

g7 <- cbind(g7.c1, g7.c2, g7.c3, g7.c4)
grid.newpage()
grid.draw(g7)


## Plot climate model intercepts
(coef.int.sr.mat.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Species richness" & Climate == "Mean temperature"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.sr.mat.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Species richness" & Climate == "Mean temperature"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.sr.td.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Species richness" & Climate == "Temperature difference"), 
                              aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.sr.td.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Species richness" & Climate == "Temperature difference"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.elev.mat.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Elevation" & Climate == "Mean temperature"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Elevation") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.lat.mat.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Latitude" & Climate == "Mean temperature"), 
                                 aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Latitude") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.elev.td.lat <- ggplot(data = subset(coef.elev, Parameter == "Intercept" & Metric == "Elevation" & Climate == "Temperature difference"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Elevation") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.int.lat.td.elev <- ggplot(data = subset(coef.lat, Parameter == "Intercept" & Metric == "Latitude" & Climate == "Temperature difference"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Latitude") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

# Combine plots of climate model intercepts (SI in manuscript)
# requires some minor revisions in a graphics editor
coef.int.sr.mat.lat <- ggplotGrob(coef.int.sr.mat.lat)
coef.int.sr.mat.elev <- ggplotGrob(coef.int.sr.mat.elev)
coef.int.sr.td.lat <- ggplotGrob(coef.int.sr.td.lat)
coef.int.sr.td.elev <- ggplotGrob(coef.int.sr.td.elev)
coef.int.elev.mat.lat <- ggplotGrob(coef.int.elev.mat.lat)
coef.int.lat.mat.elev <- ggplotGrob(coef.int.lat.mat.elev)
coef.int.elev.td.lat <- ggplotGrob(coef.int.elev.td.lat)
coef.int.lat.td.elev <- ggplotGrob(coef.int.lat.td.elev)

g8.c1 <- rbind(coef.int.elev.mat.lat,
               coef.int.sr.mat.lat,
               coef.int.elev.td.lat,
               coef.int.sr.td.lat)

g8.c1$widths <- unit.pmax(coef.int.elev.mat.lat$widths,
                          coef.int.sr.mat.lat$widths,
                          coef.int.elev.td.lat$widths,
                          coef.int.sr.td.lat$widths)

g8.c2 <- rbind(coef.int.lat.mat.elev,
               coef.int.sr.mat.elev,
               coef.int.lat.td.elev,
               coef.int.sr.td.elev)

g8.c2$widths <- unit.pmax(coef.int.lat.mat.elev$widths,
                          coef.int.sr.mat.elev$widths, 
                          coef.int.lat.td.elev$widths,
                          coef.int.sr.td.elev$widths)

g8 <- cbind(g8.c1, g8.c2)
grid.newpage()
grid.draw(g8)


## Plot climate model B1 slopes
(coef.B1.sr.mat.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Species richness" & Climate == "Mean temperature"),
                              aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.sr.mat.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Species richness" & Climate == "Mean temperature"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.sr.td.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Species richness" & Climate == "Temperature difference"), 
                             aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.sr.td.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Species richness" & Climate == "Temperature difference"), 
                              aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.elev.mat.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Elevation" & Climate == "Mean temperature"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Elevation") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.lat.mat.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Latitude" & Climate == "Mean temperature"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Latitude") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.elev.td.lat <- ggplot(data = subset(coef.elev, Parameter == "B1" & Metric == "Elevation" & Climate == "Temperature difference"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Elevation") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B1.lat.td.elev <- ggplot(data = subset(coef.lat, Parameter == "B1" & Metric == "Latitude" & Climate == "Temperature difference"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Latitude") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

# Combine plots of climate model B1 slopes (SI in manuscript)
# requires some minor revisions in a graphics editor
coef.B1.sr.mat.lat <- ggplotGrob(coef.B1.sr.mat.lat)
coef.B1.sr.mat.elev <- ggplotGrob(coef.B1.sr.mat.elev)
coef.B1.sr.td.lat <- ggplotGrob(coef.B1.sr.td.lat)
coef.B1.sr.td.elev <- ggplotGrob(coef.B1.sr.td.elev)
coef.B1.elev.mat.lat <- ggplotGrob(coef.B1.elev.mat.lat)
coef.B1.lat.mat.elev <- ggplotGrob(coef.B1.lat.mat.elev)
coef.B1.elev.td.lat <- ggplotGrob(coef.B1.elev.td.lat)
coef.B1.lat.td.elev <- ggplotGrob(coef.B1.lat.td.elev)

g9.c1 <- rbind(coef.B1.elev.mat.lat,
               coef.B1.sr.mat.lat,
               coef.B1.elev.td.lat,
               coef.B1.sr.td.lat)

g9.c1$widths <- unit.pmax(coef.B1.elev.mat.lat$widths,
                          coef.B1.sr.mat.lat$widths, 
                          coef.B1.elev.td.lat$widths,
                          coef.B1.sr.td.lat$widths)

g9.c2 <- rbind(coef.B1.lat.mat.elev,
               coef.B1.sr.mat.elev,
               coef.B1.lat.td.elev,
               coef.B1.sr.td.elev)

g9.c2$widths <- unit.pmax(coef.B1.lat.mat.elev$widths,
                          coef.B1.sr.mat.elev$widths,
                          coef.B1.lat.td.elev$widths,
                          coef.B1.sr.td.elev$widths)

g9 <- cbind(g9.c1, g9.c2)
grid.newpage()
grid.draw(g9)

                 
## Plot climate model B2 curvatures
(coef.B2.sr.mat.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Species richness" & Climate == "Mean temperature"),
                              aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.sr.mat.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Species richness" & Climate == "Mean temperature"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.sr.td.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Species richness" & Climate == "Temperature difference"), 
                             aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.sr.td.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Species richness" & Climate == "Temperature difference"), 
                              aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Species richness") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.elev.mat.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Elevation" & Climate == "Mean temperature"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Elevation") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.lat.mat.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Latitude" & Climate == "Mean temperature"), 
                                aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Latitude") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.elev.td.lat <- ggplot(data = subset(coef.elev, Parameter == "B2" & Metric == "Elevation" & Climate == "Temperature difference"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Elevation") +
    scale_x_discrete("Latitudinal zone (°N)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

(coef.B2.lat.td.elev <- ggplot(data = subset(coef.lat, Parameter == "B2" & Metric == "Latitude" & Climate == "Temperature difference"), 
                               aes(x = Class, y = Estimate, colour = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
    geom_point(alpha = 0.5, size = 6) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .2, position = position_dodge(.9)) +
    scale_y_continuous("Latitude") +
    scale_x_discrete("Elevational zone (masl)") +
    workingtheme +
    theme(axis.text.x = element_text(angle = 67.5, hjust = 1)) +
    scale_colour_manual(name = "Class", values = a60505_6e068f))

# Combine plots of climate model B2 curvatures (SI in manuscript)
# requires some minor revisions in a graphics editor
coef.B2.sr.mat.lat <- ggplotGrob(coef.B2.sr.mat.lat)
coef.B2.sr.mat.elev <- ggplotGrob(coef.B2.sr.mat.elev)
coef.B2.sr.td.lat <- ggplotGrob(coef.B2.sr.td.lat)
coef.B2.sr.td.elev <- ggplotGrob(coef.B2.sr.td.elev)
coef.B2.elev.mat.lat <- ggplotGrob(coef.B2.elev.mat.lat)
coef.B2.lat.mat.elev <- ggplotGrob(coef.B2.lat.mat.elev)
coef.B2.elev.td.lat <- ggplotGrob(coef.B2.elev.td.lat)
coef.B2.lat.td.elev <- ggplotGrob(coef.B2.lat.td.elev)

g10.c1 <- rbind(coef.B2.elev.mat.lat,
                coef.B2.sr.mat.lat,
                coef.B2.elev.td.lat,
                coef.B2.sr.td.lat)

g10.c1$widths <- unit.pmax(coef.B2.elev.mat.lat$widths,
                          coef.B2.sr.mat.lat$widths,
                          coef.B2.elev.td.lat$widths,
                          coef.B2.sr.td.lat$widths)

g10.c2 <- rbind(coef.B2.lat.mat.elev,
                coef.B2.sr.mat.elev,
                coef.B2.lat.td.elev,
                coef.B2.sr.td.elev)

g10.c2$widths <- unit.pmax(coef.B2.lat.mat.elev$widths,
                          coef.B2.sr.mat.elev$widths,
                          coef.B2.lat.td.elev$widths,
                          coef.B2.sr.td.elev$widths)

g10 <- cbind(g10.c1, g10.c2)
grid.newpage()
grid.draw(g10)


## Model mean annual temperature relationships to other diversity metrics

# PD ~ poly MAT multilevel lat effect model
bayes.pd.MAT.lat.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | lat_class),
                            data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                            prior = c(set_prior("normal(0,5)", class = "b")),
                            warmup = 1000, iter = 3000, chains = 4,
                            control = list(adapt_delta = 0.99),
                            seed = 1)

waic(bayes.pd.MAT.lat.2re, pointwise = TRUE)
loo(bayes.pd.MAT.lat.2re, pointwise = TRUE)
summary(bayes.pd.MAT.lat.2re)
prior_summary(bayes.pd.MAT.lat.2re)
coef(bayes.pd.MAT.lat.2re)
fixef(bayes.pd.MAT.lat.2re)
bayes_R2(bayes.pd.MAT.lat.2re)
control_params(bayes.pd.MAT.lat.2re)

# PD ~ poly MAT multilevel lat effect plotting
pp_check(bayes.pd.MAT.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.MAT.lat.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.MAT.lat.2re, terms = c("MAT.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.pd.MAT.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.PD"), alpha = 0.1, aes(x = MAT, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Phylogenetic richness", breaks = c(0, 160, 320, 480, 640)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(0, 647.1955), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.pd.MAT.lat, file = "p1.pd.MAT.lat_post.linpred.re2.rdata")

# PD ~ poly MAT multilevel elev effect model
bayes.pd.MAT.elev.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | elev_class),
                             data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                             prior = c(set_prior("normal(0,5)", class = "b")),
                             warmup = 1000, iter = 3000, chains = 4,
                             control = list(adapt_delta = 0.99),
                             seed = 1)

waic(bayes.pd.MAT.elev.2re, pointwise = TRUE)
loo(bayes.pd.MAT.elev.2re, pointwise = TRUE)
summary(bayes.pd.MAT.elev.2re)
prior_summary(bayes.pd.MAT.elev.2re)
coef(bayes.pd.MAT.elev.2re)
fixef(bayes.pd.MAT.elev.2re)
bayes_R2(bayes.pd.MAT.elev.2re)
control_params(bayes.pd.MAT.elev.2re)

# PD ~ poly MAT multilevel elev effect plotting
pp_check(bayes.pd.MAT.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.MAT.elev.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.MAT.elev.2re, terms = c("MAT.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.pd.MAT.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.PD"), alpha = 0.1, aes(x = MAT, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Phylogenetic richness", breaks = c(0, 160, 320, 480, 640)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(0, 647.1955), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.pd.MAT.elev, file = "p1.pd.MAT.elev_post.linpred.re2.rdata")

# FRic ~ poly MAT multilevel lat effect model
mix <- mixture(Beta, Beta)
bayes.FRic.MAT.lat.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | lat_class),
                              data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                              prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                        prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                              warmup = 1000, iter = 9000, chains = 4, inits = 0,
                              control = list(adapt_delta = 0.99, max_treedepth = 15),
                              seed = 1)

waic(bayes.FRic.MAT.lat.2re, pointwise = TRUE)
loo(bayes.FRic.MAT.lat.2re, pointwise = TRUE)
summary(bayes.FRic.MAT.lat.2re)
prior_summary(bayes.FRic.MAT.lat.2re)
coef(bayes.FRic.MAT.lat.2re)
fixef(bayes.FRic.MAT.lat.2re)
bayes_R2(bayes.FRic.MAT.lat.2re)
control_params(bayes.FRic.MAT.lat.2re)

# FRic ~ poly MAT multilevel lat effect plotting
pp_check(bayes.FRic.MAT.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.MAT.lat.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.MAT.lat.2re, terms = c("MAT.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.FRic.MAT.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.FRic"), alpha = 0.1, aes(x = MAT, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Functional richness", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(0, 1), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.FRic.MAT.lat, file = "p1.FRic.MAT.lat_post.linpred.re2.rdata")

# FRic ~ poly MAT multilevel elev effect model
mix <- mixture(Beta, Beta)
bayes.FRic.MAT.elev.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | elev_class),
                               data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                               prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                         prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                               warmup = 1000, iter = 9000, chains = 4, inits = 0,
                               control = list(adapt_delta = 0.99, max_treedepth = 15),
                               seed = 1)

waic(bayes.FRic.MAT.elev.2re, pointwise = TRUE)
loo(bayes.FRic.MAT.elev.2re, pointwise = TRUE)
summary(bayes.FRic.MAT.elev.2re)
prior_summary(bayes.FRic.MAT.elev.2re)
coef(bayes.FRic.MAT.elev.2re)
fixef(bayes.FRic.MAT.elev.2re)
bayes_R2(bayes.FRic.MAT.elev.2re)
control_params(bayes.FRic.MAT.elev.2re)

# FRic ~ poly MAT multilevel elev effect plotting
pp_check(bayes.FRic.MAT.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.MAT.elev.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.MAT.elev.2re, terms = c("MAT.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.FRic.MAT.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.FRic"), alpha = 0.1, aes(x = MAT, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Functional richness", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(0, 1), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.FRic.MAT.elev, file = "p1.FRic.MAT.elev_post.linpred.re2.rdata")

# CWM.Bl ~ poly MAT multilevel lat effect model
bayes.CWM.Bl.MAT.lat.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | lat_class),
                                data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                                prior = c(set_prior("normal(0,5)", class = "b")),
                                warmup = 1000, iter = 3000, chains = 4,
                                control = list(adapt_delta = 0.99),
                                seed = 1)

waic(bayes.CWM.Bl.MAT.lat.2re, pointwise = TRUE)
loo(bayes.CWM.Bl.MAT.lat.2re, pointwise = TRUE)
summary(bayes.CWM.Bl.MAT.lat.2re)
prior_summary(bayes.CWM.Bl.MAT.lat.2re)
coef(bayes.CWM.Bl.MAT.lat.2re)
fixef(bayes.CWM.Bl.MAT.lat.2re)
bayes_R2(bayes.CWM.Bl.MAT.lat.2re)
control_params(bayes.CWM.Bl.MAT.lat.2re)

# CWM.Bl ~ poly MAT multilevel lat effect plotting
pp_check(bayes.CWM.Bl.MAT.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.CWM.Bl.MAT.lat.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.CWM.Bl.MAT.lat.2re, terms = c("MAT.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.CWM.Bl.MAT.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "CWM.Bl"), alpha = 0.1, aes(x = MAT, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Community mean body length", breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(0, 3.5), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.CWM.Bl.MAT.lat, file = "p1.CWM.Bl.MAT.lat_post.linpred.re2.rdata")

# CWM.Bl ~ poly MAT multilevel elev effect model
bayes.CWM.Bl.MAT.elev.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | elev_class),
                                 data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                                 prior = c(set_prior("normal(0,5)", class = "b")),
                                 warmup = 1000, iter = 3000, chains = 4,
                                 control = list(adapt_delta = 0.99),
                                 seed = 1)

waic(bayes.CWM.Bl.MAT.elev.2re, pointwise = TRUE)
loo(bayes.CWM.Bl.MAT.elev.2re, pointwise = TRUE)
summary(bayes.CWM.Bl.MAT.elev.2re)
prior_summary(bayes.CWM.Bl.MAT.elev.2re)
coef(bayes.CWM.Bl.MAT.elev.2re)
fixef(bayes.CWM.Bl.MAT.elev.2re)
bayes_R2(bayes.CWM.Bl.MAT.elev.2re)
control_params(bayes.CWM.Bl.MAT.elev.2re)

# CWM.Bl ~ poly MAT multilevel elev effect plotting
pp_check(bayes.CWM.Bl.MAT.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.CWM.Bl.MAT.elev.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.CWM.Bl.MAT.elev.2re, terms = c("MAT.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.CWM.Bl.MAT.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "CWM.Bl"), alpha = 0.1, aes(x = MAT, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Community mean body length", breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(0, 3.5), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.CWM.Bl.MAT.elev, file = "p1.CWM.Bl.MAT.elev_post.linpred.re2.rdata")

# PD.obs.z ~ poly MAT multilevel lat effect model
bayes.pd.obs.z.lat.MAT.lat.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | lat_class),
                                      data = subset(lmdata, variable == "pd.obs.z.lat"),family = skew_normal(),
                                      prior = c(set_prior("normal(0,5)", class = "b")),
                                      warmup = 1000, iter = 3000, chains = 4,
                                      control = list(adapt_delta = 0.99),
                                      seed = 1)

waic(bayes.pd.obs.z.lat.MAT.lat.2re, pointwise = TRUE)
loo(bayes.pd.obs.z.lat.MAT.lat.2re, pointwise = TRUE)
summary(bayes.pd.obs.z.lat.MAT.lat.2re)
prior_summary(bayes.pd.obs.z.lat.MAT.lat.2re)
coef(bayes.pd.obs.z.lat.MAT.lat.2re)
fixef(bayes.pd.obs.z.lat.MAT.lat.2re)
bayes_R2(bayes.pd.obs.z.lat.MAT.lat.2re)
control_params(bayes.pd.obs.z.lat.MAT.lat.2re)

# PD.obs.z ~ poly MAT multilevel lat effect plotting
pp_check(bayes.pd.obs.z.lat.MAT.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.obs.z.lat.MAT.lat.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.obs.z.lat.MAT.lat.2re, terms = c("MAT.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.pd.obs.z.lat.MAT.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "pd.obs.z.lat"), alpha = 0.1, aes(x = MAT, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic richness", breaks = c(-5, -3, -1, 1, 3)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(-4.75, 3.361625), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.pd.obs.z.lat.MAT.lat, file = "p1.pd.obs.z.lat.MAT.lat_post.linpred.re2.rdata")

# PD.obs.z ~ poly MAT multilevel elev effect model
bayes.pd.obs.z.ele.MAT.elev.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | elev_class),
                                       data = subset(lmdata, variable == "pd.obs.z.ele"),family = skew_normal(),
                                       prior = c(set_prior("normal(0,5)", class = "b")),
                                       warmup = 1000, iter = 3000, chains = 4,
                                       control = list(adapt_delta = 0.99),
                                       seed = 1)

waic(bayes.pd.obs.z.ele.MAT.elev.2re, pointwise = TRUE)
loo(bayes.pd.obs.z.ele.MAT.elev.2re, pointwise = TRUE)
summary(bayes.pd.obs.z.ele.MAT.elev.2re)
prior_summary(bayes.pd.obs.z.ele.MAT.elev.2re)
coef(bayes.pd.obs.z.ele.MAT.elev.2re)
fixef(bayes.pd.obs.z.ele.MAT.elev.2re)
bayes_R2(bayes.pd.obs.z.ele.MAT.elev.2re)
control_params(bayes.pd.obs.z.ele.MAT.elev.2re)

# PD.obs.z ~ poly MAT multilevel elev effect plotting
pp_check(bayes.pd.obs.z.ele.MAT.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.obs.z.ele.MAT.elev.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.obs.z.ele.MAT.elev.2re, terms = c("MAT.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.pd.obs.z.ele.MAT.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "pd.obs.z.ele"), alpha = 0.1, aes(x = MAT, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic richness", breaks = c(-5, -3, -1, 1, 3)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(-4.75, 3.361625), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.pd.obs.z.ele.MAT.elev, file = "p1.pd.obs.z.ele.MAT.elev_post.linpred.re2.rdata")

# MPD.obs.z ~ poly MAT multilevel lat effect model
bayes.mpd.obs.z.lat.MAT.lat.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | lat_class),
                                       data = subset(lmdata, variable == "mpd.obs.z.lat"),family = skew_normal(),
                                       prior = c(set_prior("normal(0,5)", class = "b")),
                                       warmup = 1000, iter = 3000, chains = 4,
                                       control = list(adapt_delta = 0.99),
                                       seed = 1)

waic(bayes.mpd.obs.z.lat.MAT.lat.2re, pointwise = TRUE)
loo(bayes.mpd.obs.z.lat.MAT.lat.2re, pointwise = TRUE)
summary(bayes.mpd.obs.z.lat.MAT.lat.2re)
prior_summary(bayes.mpd.obs.z.lat.MAT.lat.2re)
coef(bayes.mpd.obs.z.lat.MAT.lat.2re)
fixef(bayes.mpd.obs.z.lat.MAT.lat.2re)
bayes_R2(bayes.mpd.obs.z.lat.MAT.lat.2re)
control_params(bayes.mpd.obs.z.lat.MAT.lat.2re)

# MPD.obs.z ~ poly MAT multilevel lat effect plotting
pp_check(bayes.mpd.obs.z.lat.MAT.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.mpd.obs.z.lat.MAT.lat.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.mpd.obs.z.lat.MAT.lat.2re, terms = c("MAT.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.mpd.obs.z.lat.MAT.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "mpd.obs.z.lat"), alpha = 0.1, aes(x = MAT, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic MPD", breaks = c(-5, -2.5, 0, 2.5)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(-5.255725, 2.434389), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.mpd.obs.z.lat.MAT.lat, file = "p1.mpd.obs.z.lat.MAT.lat_post.linpred.re2.rdata")

# MPD.obs.z ~ poly MAT multilevel elev effect model
bayes.mpd.obs.z.ele.MAT.elev.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | elev_class),
                                        data = subset(lmdata, variable == "mpd.obs.z.ele"),family = skew_normal(),
                                        prior = c(set_prior("normal(0,5)", class = "b")),
                                        warmup = 1000, iter = 3000, chains = 4,
                                        control = list(adapt_delta = 0.99),
                                        seed = 1)

waic(bayes.mpd.obs.z.ele.MAT.elev.2re, pointwise = TRUE)
loo(bayes.mpd.obs.z.ele.MAT.elev.2re, pointwise = TRUE)
summary(bayes.mpd.obs.z.ele.MAT.elev.2re)
prior_summary(bayes.mpd.obs.z.ele.MAT.elev.2re)
coef(bayes.mpd.obs.z.ele.MAT.elev.2re)
fixef(bayes.mpd.obs.z.ele.MAT.elev.2re)
bayes_R2(bayes.mpd.obs.z.ele.MAT.elev.2re)
control_params(bayes.mpd.obs.z.ele.MAT.elev.2re)

# MPD.obs.z ~ poly MAT multilevel elev effect plotting
pp_check(bayes.mpd.obs.z.ele.MAT.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.mpd.obs.z.ele.MAT.elev.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.mpd.obs.z.ele.MAT.elev.2re, terms = c("MAT.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.mpd.obs.z.ele.MAT.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "mpd.obs.z.ele"), alpha = 0.1, aes(x = MAT, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic MPD", breaks = c(-5, -2.5, 0, 2.5)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(-5.255725, 2.434389), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.mpd.obs.z.ele.MAT.elev, file = "p1.mpd.obs.z.ele.MAT.elev_post.linpred.re2.rdata")

# FRic.obs.z ~ poly MAT multilevel lat effect model
bayes.FRic.obs.z.lat.MAT.lat.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | lat_class),
                                        data = subset(lmdata, variable == "FRic.obs.z.lat"), family = skew_normal(),
                                        prior = c(set_prior("normal(0,5)", class = "b")),
                                        warmup = 1000, iter = 3000, chains = 4,
                                        control = list(adapt_delta = 0.99),
                                        seed = 1)

waic(bayes.FRic.obs.z.lat.MAT.lat.2re, pointwise = TRUE)
loo(bayes.FRic.obs.z.lat.MAT.lat.2re, pointwise = TRUE)
summary(bayes.FRic.obs.z.lat.MAT.lat.2re)
prior_summary(bayes.FRic.obs.z.lat.MAT.lat.2re)
coef(bayes.FRic.obs.z.lat.MAT.lat.2re)
fixef(bayes.FRic.obs.z.lat.MAT.lat.2re)
bayes_R2(bayes.FRic.obs.z.lat.MAT.lat.2re)
control_params(bayes.FRic.obs.z.lat.MAT.lat.2re)

# FRic.obs.z ~ poly MAT multilevel lat effect plotting
pp_check(bayes.FRic.obs.z.lat.MAT.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.obs.z.lat.MAT.lat.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.obs.z.lat.MAT.lat.2re, terms = c("MAT.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.FRic.obs.z.lat.MAT.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FRic.obs.z.lat"), alpha = 0.1, aes(x = MAT, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional richness", breaks = c(-4, -2, 0, 2, 4)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(-4.067616, 4), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.FRic.obs.z.lat.MAT.lat, file = "p1.FRic.obs.z.lat.MAT.lat_post.linpred.re2.rdata")

# FRic.obs.z ~ poly MAT multilevel elev effect model
bayes.FRic.obs.z.ele.MAT.elev.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | elev_class),
                                         data = subset(lmdata, variable == "FRic.obs.z.ele"), family = skew_normal(),
                                         prior = c(set_prior("normal(0,5)", class = "b")),
                                         warmup = 1000, iter = 3000, chains = 4,
                                         control = list(adapt_delta = 0.99),
                                         seed = 1)

waic(bayes.FRic.obs.z.ele.MAT.elev.2re, pointwise = TRUE)
loo(bayes.FRic.obs.z.ele.MAT.elev.2re, pointwise = TRUE)
summary(bayes.FRic.obs.z.ele.MAT.elev.2re)
prior_summary(bayes.FRic.obs.z.ele.MAT.elev.2re)
coef(bayes.FRic.obs.z.ele.MAT.elev.2re)
fixef(bayes.FRic.obs.z.ele.MAT.elev.2re)
bayes_R2(bayes.FRic.obs.z.ele.MAT.elev.2re)
control_params(bayes.FRic.obs.z.ele.MAT.elev.2re)

# FRic.obs.z ~ poly MAT multilevel elev effect plotting
pp_check(bayes.FRic.obs.z.ele.MAT.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.obs.z.ele.MAT.elev.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.obs.z.ele.MAT.elev.2re, terms = c("MAT.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.FRic.obs.z.ele.MAT.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FRic.obs.z.ele"), alpha = 0.1, aes(x = MAT, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional richness", breaks = c(-4, -2, 0, 2, 4)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(-4.067616, 4), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.FRic.obs.z.ele.MAT.elev, file = "p1.FRic.obs.z.ele.MAT.elev_post.linpred.re2.rdata")

# FDis.obs.z ~ poly MAT multilevel lat effect model
bayes.FDis.obs.z.lat.MAT.lat.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | lat_class),
                                        data = subset(lmdata, variable == "FDis.obs.z.lat"), family = skew_normal(),
                                        prior = c(set_prior("normal(0,5)", class = "b")),
                                        warmup = 1000, iter = 3000, chains = 4,
                                        control = list(adapt_delta = 0.99),
                                        seed = 1)

waic(bayes.FDis.obs.z.lat.MAT.lat.2re, pointwise = TRUE)
loo(bayes.FDis.obs.z.lat.MAT.lat.2re, pointwise = TRUE)
summary(bayes.FDis.obs.z.lat.MAT.lat.2re)
prior_summary(bayes.FDis.obs.z.lat.MAT.lat.2re)
coef(bayes.FDis.obs.z.lat.MAT.lat.2re)
fixef(bayes.FDis.obs.z.lat.MAT.lat.2re)
bayes_R2(bayes.FDis.obs.z.lat.MAT.lat.2re)
control_params(bayes.FDis.obs.z.lat.MAT.lat.2re)

# FDis.obs.z ~ poly MAT multilevel lat effect plotting
pp_check(bayes.FDis.obs.z.lat.MAT.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FDis.obs.z.lat.MAT.lat.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FDis.obs.z.lat.MAT.lat.2re, terms = c("MAT.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.FDis.obs.z.lat.MAT.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FDis.obs.z.lat"), alpha = 0.1, aes(x = MAT, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional dispersion", breaks = c(-3, -1.5, 0, 1.5, 3)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(-3.079376, 2.948407), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.FDis.obs.z.lat.MAT.lat, file = "p1.FDis.obs.z.lat.MAT.lat_post.linpred.re2.rdata")

# FDis.obs.z ~ poly MAT multilevel elev effect model
bayes.FDis.obs.z.ele.MAT.elev.2re <- brm(formula = value ~ poly(MAT.scale, 2, raw = FALSE) + (1 + poly(MAT.scale, 2, raw = FALSE) | elev_class),
                                         data = subset(lmdata, variable == "FDis.obs.z.ele"), family = skew_normal(),
                                         prior = c(set_prior("normal(0,5)", class = "b")),
                                         warmup = 1000, iter = 3000, chains = 4,
                                         control = list(adapt_delta = 0.99),
                                         seed = 1)

waic(bayes.FDis.obs.z.ele.MAT.elev.2re, pointwise = TRUE)
loo(bayes.FDis.obs.z.ele.MAT.elev.2re, pointwise = TRUE)
summary(bayes.FDis.obs.z.ele.MAT.elev.2re)
prior_summary(bayes.FDis.obs.z.ele.MAT.elev.2re)
coef(bayes.FDis.obs.z.ele.MAT.elev.2re)
fixef(bayes.FDis.obs.z.ele.MAT.elev.2re)
bayes_R2(bayes.FDis.obs.z.ele.MAT.elev.2re)
control_params(bayes.FDis.obs.z.ele.MAT.elev.2re)

# FDis.obs.z ~ poly MAT multilevel elev effect plotting
pp_check(bayes.FDis.obs.z.ele.MAT.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FDis.obs.z.ele.MAT.elev.2re, terms = c("MAT.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FDis.obs.z.ele.MAT.elev.2re, terms = c("MAT.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$MAT) + mean(lmdata$MAT) 
prb$x <- prb$x * sd(lmdata$MAT) + mean(lmdata$MAT) 

(p1.FDis.obs.z.ele.MAT.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FDis.obs.z.ele"), alpha = 0.1, aes(x = MAT, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional dispersion", breaks = c(-3, -1.5, 0, 1.5, 3)) +
    scale_x_continuous("Mean annual temperature (°C)", breaks = c(-8, -4, 0, 4, 8, 12)) +
    coord_cartesian(ylim = c(-3.079376, 2.948407), xlim = c(-7.859615, 11.5)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.FDis.obs.z.ele.MAT.elev, file = "p1.FDis.obs.z.ele.MAT.elev_post.linpred.re2.rdata")


## Model mean temperature difference relationships to other diversity metrics

# PD ~ poly TD multilevel lat effect model
bayes.pd.TD.lat.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | lat_class),
                           data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                           prior = c(set_prior("normal(0,5)", class = "b")),
                           warmup = 1000, iter = 3000, chains = 4,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.pd.TD.lat.2re, pointwise = TRUE)
loo(bayes.pd.TD.lat.2re, pointwise = TRUE)
summary(bayes.pd.TD.lat.2re)
prior_summary(bayes.pd.TD.lat.2re)
coef(bayes.pd.TD.lat.2re)
fixef(bayes.pd.TD.lat.2re)
bayes_R2(bayes.pd.TD.lat.2re)
control_params(bayes.pd.TD.lat.2re)

# PD ~ poly TD multilevel lat effect plotting
pp_check(bayes.pd.TD.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.TD.lat.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.TD.lat.2re, terms = c("TD.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.pd.TD.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.PD"), alpha = 0.1, aes(x = TD, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Phylogenetic richness", breaks = c(0, 160, 320, 480, 640)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(0, 647.1955), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.pd.TD.lat, file = "p1.pd.TD.lat_post.linpred.re2.rdata")

# PD ~ poly TD multilevel elev effect model
bayes.pd.TD.elev.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | elev_class),
                            data = subset(lmdata, variable == "Glob.PD"), family = Gamma(link = "log"),
                            prior = c(set_prior("normal(0,5)", class = "b")),
                            warmup = 1000, iter = 3000, chains = 4,
                            control = list(adapt_delta = 0.99),
                            seed = 1)

waic(bayes.pd.TD.elev.2re, pointwise = TRUE)
loo(bayes.pd.TD.elev.2re, pointwise = TRUE)
summary(bayes.pd.TD.elev.2re)
prior_summary(bayes.pd.TD.elev.2re)
coef(bayes.pd.TD.elev.2re)
fixef(bayes.pd.TD.elev.2re)
bayes_R2(bayes.pd.TD.elev.2re)
control_params(bayes.pd.TD.elev.2re)

# PD ~ poly TD multilevel elev effect plotting
pp_check(bayes.pd.TD.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.TD.elev.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.TD.elev.2re, terms = c("TD.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.pd.TD.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.PD"), alpha = 0.1, aes(x = TD, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Phylogenetic richness", breaks = c(0, 160, 320, 480, 640)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(0, 647.1955), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.pd.TD.elev, file = "p1.pd.TD.elev_post.linpred.re2.rdata")

# FRic ~ poly TD multilevel lat effect model
mix <- mixture(Beta, Beta)
bayes.FRic.TD.lat.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | lat_class),
                             data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                             prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                       prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                             warmup = 1000, iter = 9000, chains = 4, inits = 0,
                             control = list(adapt_delta = 0.99, max_treedepth = 15),
                             seed = 1)

waic(bayes.FRic.TD.lat.2re, pointwise = TRUE)
loo(bayes.FRic.TD.lat.2re, pointwise = TRUE)
summary(bayes.FRic.TD.lat.2re)
prior_summary(bayes.FRic.TD.lat.2re)
coef(bayes.FRic.TD.lat.2re)
fixef(bayes.FRic.TD.lat.2re)
bayes_R2(bayes.FRic.TD.lat.2re)
control_params(bayes.FRic.TD.lat.2re)

# FRic ~ poly TD multilevel lat effect plotting
pp_check(bayes.FRic.TD.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.TD.lat.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.TD.lat.2re, terms = c("TD.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.FRic.TD.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.FRic"), alpha = 0.1, aes(x = TD, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Functional richness", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(0, 1), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.FRic.TD.lat, file = "p1.FRic.TD.lat_post.linpred.re2.rdata")

# FRic ~ poly TD multilevel elev effect model
mix <- mixture(Beta, Beta)
bayes.FRic.TD.elev.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | elev_class),
                              data = subset(lmdata, variable == "Glob.FRic"), family = mix,
                              prior = c(prior("normal(0,5)", class = "b", dpar = mu1), prior("normal(0,5)", class = "b", dpar = mu2), 
                                        prior("student_t(3,-0.5,0.2)", class = "Intercept", dpar = mu1), prior("student_t(3,0.5,0.2)", class = "Intercept", dpar = mu2)),
                              warmup = 1000, iter = 9000, chains = 4, inits = 0,
                              control = list(adapt_delta = 0.99, max_treedepth = 15),
                              seed = 1)

waic(bayes.FRic.TD.elev.2re, pointwise = TRUE)
loo(bayes.FRic.TD.elev.2re, pointwise = TRUE)
summary(bayes.FRic.TD.elev.2re)
prior_summary(bayes.FRic.TD.elev.2re)
coef(bayes.FRic.TD.elev.2re)
fixef(bayes.FRic.TD.elev.2re)
bayes_R2(bayes.FRic.TD.elev.2re)
control_params(bayes.FRic.TD.elev.2re)

# FRic ~ poly TD multilevel elev effect plotting
pp_check(bayes.FRic.TD.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.TD.elev.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.TD.elev.2re, terms = c("TD.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.FRic.TD.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "Glob.FRic"), alpha = 0.1, aes(x = TD, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Functional richness", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(0, 1), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.FRic.TD.elev, file = "p1.FRic.TD.elev_post.linpred.re2.rdata")

# CWM.Bl ~ poly TD multilevel lat effect model
bayes.CWM.Bl.TD.lat.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | lat_class),
                           data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                           prior = c(set_prior("normal(0,5)", class = "b")),
                           warmup = 1000, iter = 3000, chains = 4,
                           control = list(adapt_delta = 0.99),
                           seed = 1)

waic(bayes.CWM.Bl.TD.lat.2re, pointwise = TRUE)
loo(bayes.CWM.Bl.TD.lat.2re, pointwise = TRUE)
summary(bayes.CWM.Bl.TD.lat.2re)
prior_summary(bayes.CWM.Bl.TD.lat.2re)
coef(bayes.CWM.Bl.TD.lat.2re)
fixef(bayes.CWM.Bl.TD.lat.2re)
bayes_R2(bayes.CWM.Bl.TD.lat.2re)
control_params(bayes.CWM.Bl.TD.lat.2re)

# CWM.Bl ~ poly TD multilevel lat effect plotting
pp_check(bayes.CWM.Bl.TD.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.CWM.Bl.TD.lat.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.CWM.Bl.TD.lat.2re, terms = c("TD.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.CWM.Bl.TD.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "CWM.Bl"), alpha = 0.1, aes(x = TD, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Community mean body length", breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(0, 3.5), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.CWM.Bl.TD.lat, file = "p1.CWM.Bl.TD.lat_post.linpred.re2.rdata")

# CWM.Bl ~ poly TD multilevel elev effect model
bayes.CWM.Bl.TD.elev.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | elev_class),
                            data = subset(lmdata, variable == "CWM.Bl"), family = Gamma(link = "log"),
                            prior = c(set_prior("normal(0,5)", class = "b")),
                            warmup = 1000, iter = 3000, chains = 4,
                            control = list(adapt_delta = 0.99),
                            seed = 1)

waic(bayes.CWM.Bl.TD.elev.2re, pointwise = TRUE)
loo(bayes.CWM.Bl.TD.elev.2re, pointwise = TRUE)
summary(bayes.CWM.Bl.TD.elev.2re)
prior_summary(bayes.CWM.Bl.TD.elev.2re)
coef(bayes.CWM.Bl.TD.elev.2re)
fixef(bayes.CWM.Bl.TD.elev.2re)
bayes_R2(bayes.CWM.Bl.TD.elev.2re)
control_params(bayes.CWM.Bl.TD.elev.2re)

# CWM.Bl ~ poly TD multilevel elev effect plotting
pp_check(bayes.CWM.Bl.TD.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.CWM.Bl.TD.elev.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.CWM.Bl.TD.elev.2re, terms = c("TD.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.CWM.Bl.TD.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "CWM.Bl"), alpha = 0.1, aes(x = TD, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Community mean body length", breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(0, 3.5), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.CWM.Bl.TD.elev, file = "p1.CWM.Bl.TD.elev_post.linpred.re2.rdata")

# PD.obs.z ~ poly TD multilevel lat effect model
bayes.pd.obs.z.lat.TD.lat.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | lat_class),
                               data = subset(lmdata, variable == "pd.obs.z.lat"),family = skew_normal(),
                               prior = c(set_prior("normal(0,5)", class = "b")),
                               warmup = 1000, iter = 3000, chains = 4,
                               control = list(adapt_delta = 0.99),
                               seed = 1)

waic(bayes.pd.obs.z.lat.TD.lat.2re, pointwise = TRUE)
loo(bayes.pd.obs.z.lat.TD.lat.2re, pointwise = TRUE)
summary(bayes.pd.obs.z.lat.TD.lat.2re)
prior_summary(bayes.pd.obs.z.lat.TD.lat.2re)
coef(bayes.pd.obs.z.lat.TD.lat.2re)
fixef(bayes.pd.obs.z.lat.TD.lat.2re)
bayes_R2(bayes.pd.obs.z.lat.TD.lat.2re)
control_params(bayes.pd.obs.z.lat.TD.lat.2re)

# PD.obs.z ~ poly TD multilevel lat effect plotting
pp_check(bayes.pd.obs.z.lat.TD.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.obs.z.lat.TD.lat.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.obs.z.lat.TD.lat.2re, terms = c("TD.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.pd.obs.z.lat.TD.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "pd.obs.z.lat"), alpha = 0.1, aes(x = TD, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic richness", breaks = c(-5, -3, -1, 1, 3)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(-4.75, 3.361625), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.pd.obs.z.lat.TD.lat, file = "p1.pd.obs.z.lat.TD.lat_post.linpred.re2.rdata")

# PD.obs.z ~ poly TD multilevel elev effect model
bayes.pd.obs.z.ele.TD.elev.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | elev_class),
                                data = subset(lmdata, variable == "pd.obs.z.ele"),family = skew_normal(),
                                prior = c(set_prior("normal(0,5)", class = "b")),
                                warmup = 1000, iter = 3000, chains = 4,
                                control = list(adapt_delta = 0.99),
                                seed = 1)

waic(bayes.pd.obs.z.ele.TD.elev.2re, pointwise = TRUE)
loo(bayes.pd.obs.z.ele.TD.elev.2re, pointwise = TRUE)
summary(bayes.pd.obs.z.ele.TD.elev.2re)
prior_summary(bayes.pd.obs.z.ele.TD.elev.2re)
coef(bayes.pd.obs.z.ele.TD.elev.2re)
fixef(bayes.pd.obs.z.ele.TD.elev.2re)
bayes_R2(bayes.pd.obs.z.ele.TD.elev.2re)
control_params(bayes.pd.obs.z.ele.TD.elev.2re)

# PD.obs.z ~ poly TD multilevel elev effect plotting
pp_check(bayes.pd.obs.z.ele.TD.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.pd.obs.z.ele.TD.elev.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.pd.obs.z.ele.TD.elev.2re, terms = c("TD.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.pd.obs.z.ele.TD.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "pd.obs.z.ele"), alpha = 0.1, aes(x = TD, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic richness", breaks = c(-5, -3, -1, 1, 3)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(-4.75, 3.361625), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.pd.obs.z.ele.TD.elev, file = "p1.pd.obs.z.ele.TD.elev_post.linpred.re2.rdata")

# MPD.obs.z ~ poly TD multilevel lat effect model
bayes.mpd.obs.z.lat.TD.lat.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | lat_class),
                                      data = subset(lmdata, variable == "mpd.obs.z.lat"),family = skew_normal(),
                                      prior = c(set_prior("normal(0,5)", class = "b")),
                                      warmup = 1000, iter = 3000, chains = 4,
                                      control = list(adapt_delta = 0.99),
                                      seed = 1)

waic(bayes.mpd.obs.z.lat.TD.lat.2re, pointwise = TRUE)
loo(bayes.mpd.obs.z.lat.TD.lat.2re, pointwise = TRUE)
summary(bayes.mpd.obs.z.lat.TD.lat.2re)
prior_summary(bayes.mpd.obs.z.lat.TD.lat.2re)
coef(bayes.mpd.obs.z.lat.TD.lat.2re)
fixef(bayes.mpd.obs.z.lat.TD.lat.2re)
bayes_R2(bayes.mpd.obs.z.lat.TD.lat.2re)
control_params(bayes.mpd.obs.z.lat.TD.lat.2re)

# MPD.obs.z ~ poly TD multilevel lat effect plotting
pp_check(bayes.mpd.obs.z.lat.TD.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.mpd.obs.z.lat.TD.lat.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.mpd.obs.z.lat.TD.lat.2re, terms = c("TD.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.mpd.obs.z.lat.TD.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "mpd.obs.z.lat"), alpha = 0.1, aes(x = TD, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic MPD", breaks = c(-5, -2.5, 0, 2.5)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(-5.255725, 2.434389), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.mpd.obs.z.lat.TD.lat, file = "p1.mpd.obs.z.lat.TD.lat_post.linpred.re2.rdata")

# MPD.obs.z ~ poly TD multilevel elev effect model
bayes.mpd.obs.z.ele.TD.elev.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | elev_class),
                                       data = subset(lmdata, variable == "mpd.obs.z.ele"),family = skew_normal(),
                                       prior = c(set_prior("normal(0,5)", class = "b")),
                                       warmup = 1000, iter = 3000, chains = 4,
                                       control = list(adapt_delta = 0.99),
                                       seed = 1)

waic(bayes.mpd.obs.z.ele.TD.elev.2re, pointwise = TRUE)
loo(bayes.mpd.obs.z.ele.TD.elev.2re, pointwise = TRUE)
summary(bayes.mpd.obs.z.ele.TD.elev.2re)
prior_summary(bayes.mpd.obs.z.ele.TD.elev.2re)
coef(bayes.mpd.obs.z.ele.TD.elev.2re)
fixef(bayes.mpd.obs.z.ele.TD.elev.2re)
bayes_R2(bayes.mpd.obs.z.ele.TD.elev.2re)
control_params(bayes.mpd.obs.z.ele.TD.elev.2re)

# MPD.obs.z ~ poly TD multilevel elev effect plotting
pp_check(bayes.mpd.obs.z.ele.TD.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.mpd.obs.z.ele.TD.elev.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.mpd.obs.z.ele.TD.elev.2re, terms = c("TD.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.mpd.obs.z.ele.TD.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "mpd.obs.z.ele"), alpha = 0.1, aes(x = TD, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized phylogenetic MPD", breaks = c(-5, -2.5, 0, 2.5)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(-5.255725, 2.434389), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.mpd.obs.z.ele.TD.elev, file = "p1.mpd.obs.z.ele.TD.elev_post.linpred.re2.rdata")

# FRic.obs.z ~ poly TD multilevel lat effect model
bayes.FRic.obs.z.lat.TD.lat.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | lat_class),
                                       data = subset(lmdata, variable == "FRic.obs.z.lat"), family = skew_normal(),
                                       prior = c(set_prior("normal(0,5)", class = "b")),
                                       warmup = 1000, iter = 3000, chains = 4,
                                       control = list(adapt_delta = 0.99),
                                       seed = 1)

waic(bayes.FRic.obs.z.lat.TD.lat.2re, pointwise = TRUE)
loo(bayes.FRic.obs.z.lat.TD.lat.2re, pointwise = TRUE)
summary(bayes.FRic.obs.z.lat.TD.lat.2re)
prior_summary(bayes.FRic.obs.z.lat.TD.lat.2re)
coef(bayes.FRic.obs.z.lat.TD.lat.2re)
fixef(bayes.FRic.obs.z.lat.TD.lat.2re)
bayes_R2(bayes.FRic.obs.z.lat.TD.lat.2re)
control_params(bayes.FRic.obs.z.lat.TD.lat.2re)

# FRic.obs.z ~ poly TD multilevel lat effect plotting
pp_check(bayes.FRic.obs.z.lat.TD.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.obs.z.lat.TD.lat.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.obs.z.lat.TD.lat.2re, terms = c("TD.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.FRic.obs.z.lat.TD.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FRic.obs.z.lat"), alpha = 0.1, aes(x = TD, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional richness", breaks = c(-4, -2, 0, 2, 4)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(-4.067616, 4), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.FRic.obs.z.lat.TD.lat, file = "p1.FRic.obs.z.lat.TD.lat_post.linpred.re2.rdata")

# FRic.obs.z ~ poly TD multilevel elev effect model
bayes.FRic.obs.z.ele.TD.elev.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | elev_class),
                                        data = subset(lmdata, variable == "FRic.obs.z.ele"), family = skew_normal(),
                                        prior = c(set_prior("normal(0,5)", class = "b")),
                                        warmup = 1000, iter = 3000, chains = 4,
                                        control = list(adapt_delta = 0.99),
                                        seed = 1)

waic(bayes.FRic.obs.z.ele.TD.elev.2re, pointwise = TRUE)
loo(bayes.FRic.obs.z.ele.TD.elev.2re, pointwise = TRUE)
summary(bayes.FRic.obs.z.ele.TD.elev.2re)
prior_summary(bayes.FRic.obs.z.ele.TD.elev.2re)
coef(bayes.FRic.obs.z.ele.TD.elev.2re)
fixef(bayes.FRic.obs.z.ele.TD.elev.2re)
bayes_R2(bayes.FRic.obs.z.ele.TD.elev.2re)
control_params(bayes.FRic.obs.z.ele.TD.elev.2re)

# FRic.obs.z ~ poly TD multilevel elev effect plotting
pp_check(bayes.FRic.obs.z.ele.TD.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FRic.obs.z.ele.TD.elev.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FRic.obs.z.ele.TD.elev.2re, terms = c("TD.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.FRic.obs.z.ele.TD.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FRic.obs.z.ele"), alpha = 0.1, aes(x = TD, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional richness", breaks = c(-4, -2, 0, 2, 4)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(-4.067616, 4), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.FRic.obs.z.ele.TD.elev, file = "p1.FRic.obs.z.ele.TD.elev_post.linpred.re2.rdata")

# FDis.obs.z ~ poly TD multilevel lat effect model
bayes.FDis.obs.z.lat.TD.lat.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | lat_class),
                                       data = subset(lmdata, variable == "FDis.obs.z.lat"), family = skew_normal(),
                                       prior = c(set_prior("normal(0,5)", class = "b")),
                                       warmup = 1000, iter = 3000, chains = 4,
                                       control = list(adapt_delta = 0.99),
                                       seed = 1)

waic(bayes.FDis.obs.z.lat.TD.lat.2re, pointwise = TRUE)
loo(bayes.FDis.obs.z.lat.TD.lat.2re, pointwise = TRUE)
summary(bayes.FDis.obs.z.lat.TD.lat.2re)
prior_summary(bayes.FDis.obs.z.lat.TD.lat.2re)
coef(bayes.FDis.obs.z.lat.TD.lat.2re)
fixef(bayes.FDis.obs.z.lat.TD.lat.2re)
bayes_R2(bayes.FDis.obs.z.lat.TD.lat.2re)
control_params(bayes.FDis.obs.z.lat.TD.lat.2re)

# FDis.obs.z ~ poly TD multilevel lat effect plotting
pp_check(bayes.FDis.obs.z.lat.TD.lat.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FDis.obs.z.lat.TD.lat.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FDis.obs.z.lat.TD.lat.2re, terms = c("TD.scale [all]", "lat_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.FDis.obs.z.lat.TD.lat <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FDis.obs.z.lat"), alpha = 0.1, aes(x = TD, y = value, colour = lat_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional dispersion", breaks = c(-3, -1.5, 0, 1.5, 3)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(-3.079376, 2.948407), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Latitude", values = a60505_6e068f))

save(p1.FDis.obs.z.lat.TD.lat, file = "p1.FDis.obs.z.lat.TD.lat_post.linpred.re2.rdata")

# FDis.obs.z ~ poly TD multilevel elev effect model
bayes.FDis.obs.z.ele.TD.elev.2re <- brm(formula = value ~ poly(TD.scale, 2, raw = FALSE) + (1 + poly(TD.scale, 2, raw = FALSE) | elev_class),
                                        data = subset(lmdata, variable == "FDis.obs.z.ele"), family = skew_normal(),
                                        prior = c(set_prior("normal(0,5)", class = "b")),
                                        warmup = 1000, iter = 3000, chains = 4,
                                        control = list(adapt_delta = 0.99),
                                        seed = 1)

waic(bayes.FDis.obs.z.ele.TD.elev.2re, pointwise = TRUE)
loo(bayes.FDis.obs.z.ele.TD.elev.2re, pointwise = TRUE)
summary(bayes.FDis.obs.z.ele.TD.elev.2re)
prior_summary(bayes.FDis.obs.z.ele.TD.elev.2re)
coef(bayes.FDis.obs.z.ele.TD.elev.2re)
fixef(bayes.FDis.obs.z.ele.TD.elev.2re)
bayes_R2(bayes.FDis.obs.z.ele.TD.elev.2re)
control_params(bayes.FDis.obs.z.ele.TD.elev.2re)

# FDis.obs.z ~ poly TD multilevel elev effect plotting
pp_check(bayes.FDis.obs.z.ele.TD.elev.2re, plotfun = "dens_overlay", nsamples = 100)

pra <- ggpredict(bayes.FDis.obs.z.ele.TD.elev.2re, terms = c("TD.scale [all]"), type = "re", ppd = FALSE)
prb <- ggpredict(bayes.FDis.obs.z.ele.TD.elev.2re, terms = c("TD.scale [all]", "elev_class"), type = "re", ppd = FALSE)

raw_data <- attr(prb, "rawdata", exact = TRUE)
!is.null(raw_data)

ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE), 3))
for (i in names(ranges)) {
  remove <- prb$group == i & prb$x < ranges[[i]][1]
  prb$x[remove] <- NA
  remove <- prb$group == i & prb$x > ranges[[i]][2]
  prb$x[remove] <- NA
}

pra$x <- pra$x * sd(lmdata$TD) + mean(lmdata$TD) 
prb$x <- prb$x * sd(lmdata$TD) + mean(lmdata$TD) 

(p1.FDis.obs.z.ele.TD.elev <- ggplot() +
    geom_ribbon(data = pra, aes(ymin = conf.low, ymax = conf.high, x = x), alpha = 0.075) +
    geom_point(data = subset(lmdata, variable == "FDis.obs.z.ele"), alpha = 0.1, aes(x = TD, y = value, colour = elev_class)) +
    geom_line(data = prb, aes(x, predicted, colour = group)) +
    geom_line(data = pra, aes(x, predicted), size = 1, linetype = "longdash") +
    scale_y_continuous("Standardized functional dispersion", breaks = c(-3, -1.5, 0, 1.5, 3)) +
    scale_x_continuous("Temperature difference (°C)", breaks = c(12, 22, 32, 42)) +
    coord_cartesian(ylim = c(-3.079376, 2.948407), xlim = c(11.95769, 42.71731)) +
    workingtheme +
    scale_colour_manual(name = "Elevation", values = a60505_6e068f))

save(p1.FDis.obs.z.ele.TD.elev, file = "p1.FDis.obs.z.ele.TD.elev_post.linpred.re2.rdata")

# Plotting mean annual temperature relationships to other diversity metrics (SI in manuscript)
# requires some minor revisions in a graphics editor
load("p1.alpha.MAT.elev_post.linpred.re2.rdata")
load("p1.alpha.MAT.lat_post.linpred.re2.rdata")
load("p1.pd.MAT.elev_post.linpred.re2.rdata")
load("p1.pd.MAT.lat_post.linpred.re2.rdata")
load("p1.mpd.obs.z.ele.MAT.elev_post.linpred.re2.rdata")
load("p1.mpd.obs.z.lat.MAT.lat_post.linpred.re2.rdata")
load("p1.FRic.MAT.elev_post.linpred.re2.rdata")
load("p1.FRic.MAT.lat_post.linpred.re2.rdata")
load("p1.FDis.obs.z.ele.MAT.elev_post.linpred.re2.rdata")
load("p1.FDis.obs.z.lat.MAT.lat_post.linpred.re2.rdata")
load("p1.CWM.Bl.MAT.elev_post.linpred.re2.rdata")
load("p1.CWM.Bl.MAT.lat_post.linpred.re2.rdata")
load("p1.pd.obs.z.ele.MAT.elev_post.linpred.re2.rdata")
load("p1.pd.obs.z.lat.MAT.lat_post.linpred.re2.rdata")
load("p1.FRic.obs.z.ele.MAT.elev_post.linpred.re2.rdata")
load("p1.FRic.obs.z.lat.MAT.lat_post.linpred.re2.rdata")

p1.alpha.MAT.elev <- ggplotGrob(p1.alpha.MAT.elev)
p1.alpha.MAT.lat <- ggplotGrob(p1.alpha.MAT.lat)
p1.pd.MAT.elev <- ggplotGrob(p1.pd.MAT.elev)
p1.pd.MAT.lat <- ggplotGrob(p1.pd.MAT.lat)
p1.mpd.obs.z.ele.MAT.elev <- ggplotGrob(p1.mpd.obs.z.ele.MAT.elev)
p1.mpd.obs.z.lat.MAT.lat <- ggplotGrob(p1.mpd.obs.z.lat.MAT.lat)
p1.FRic.MAT.elev <- ggplotGrob(p1.FRic.MAT.elev)
p1.FRic.MAT.lat <- ggplotGrob(p1.FRic.MAT.lat)
p1.FDis.obs.z.ele.MAT.elev <- ggplotGrob(p1.FDis.obs.z.ele.MAT.elev)
p1.FDis.obs.z.lat.MAT.lat <- ggplotGrob(p1.FDis.obs.z.lat.MAT.lat)
p1.CWM.Bl.MAT.elev <- ggplotGrob(p1.CWM.Bl.MAT.elev)
p1.CWM.Bl.MAT.lat <- ggplotGrob(p1.CWM.Bl.MAT.lat)
p1.pd.obs.z.ele.MAT.elev <- ggplotGrob(p1.pd.obs.z.ele.MAT.elev)
p1.pd.obs.z.lat.MAT.lat <- ggplotGrob(p1.pd.obs.z.lat.MAT.lat)
p1.FRic.obs.z.ele.MAT.elev <- ggplotGrob(p1.FRic.obs.z.ele.MAT.elev)
p1.FRic.obs.z.lat.MAT.lat <- ggplotGrob(p1.FRic.obs.z.lat.MAT.lat)

g11.c1 <- rbind(p1.alpha.MAT.elev,
                p1.FRic.MAT.elev,
                p1.pd.MAT.elev,
                p1.FDis.obs.z.ele.MAT.elev)

g11.c1$widths <- unit.pmax(p1.alpha.MAT.elev$widths,
                           p1.FRic.MAT.elev$widths,
                           p1.pd.MAT.elev$widths,
                           p1.FDis.obs.z.ele.MAT.elev$widths)

g11.c2 <- rbind(p1.CWM.Bl.MAT.elev,
                p1.FRic.obs.z.ele.MAT.elev,
                p1.pd.obs.z.ele.MAT.elev,
                p1.mpd.obs.z.ele.MAT.elev)

g11.c2$widths <- unit.pmax(p1.CWM.Bl.MAT.elev$widths,
                           p1.FRic.obs.z.ele.MAT.elev$widths,
                           p1.pd.obs.z.ele.MAT.elev$widths,
                           p1.mpd.obs.z.ele.MAT.elev$widths)

g11.c3 <- rbind(p1.alpha.MAT.lat,
                p1.FRic.MAT.lat,
                p1.pd.MAT.lat,
                p1.FDis.obs.z.lat.MAT.lat)

g11.c3$widths <- unit.pmax(p1.alpha.MAT.lat$widths,
                           p1.FRic.MAT.lat$widths,
                           p1.pd.MAT.lat$widths,
                           p1.FDis.obs.z.lat.MAT.lat$widths)

g11.c4 <- rbind(p1.CWM.Bl.MAT.lat,
                p1.FRic.obs.z.lat.MAT.lat,
                p1.pd.obs.z.lat.MAT.lat,
                p1.mpd.obs.z.lat.MAT.lat)

g11.c4$widths <- unit.pmax(p1.CWM.Bl.MAT.lat$widths,
                           p1.FRic.obs.z.lat.MAT.lat$widths,
                           p1.pd.obs.z.lat.MAT.lat$widths,
                           p1.mpd.obs.z.lat.MAT.lat$widths)

g11 <- cbind(g11.c1, g11.c2, g11.c3, g11.c4)
grid.newpage()
grid.draw(g11)

# Plotting mean temperature difference relationships to other diversity metrics (SI in manuscript)
# requires some minor revisions in a graphics editor
load("p1.alpha.TD.elev_post.linpred.re2.rdata")
load("p1.alpha.TD.lat_post.linpred.re2.rdata")
load("p1.pd.TD.elev_post.linpred.re2.rdata")
load("p1.pd.TD.lat_post.linpred.re2.rdata")
load("p1.mpd.obs.z.ele.TD.elev_post.linpred.re2.rdata")
load("p1.mpd.obs.z.lat.TD.lat_post.linpred.re2.rdata")
load("p1.FRic.TD.elev_post.linpred.re2.rdata")
load("p1.FRic.TD.lat_post.linpred.re2.rdata")
load("p1.FDis.obs.z.ele.TD.elev_post.linpred.re2.rdata")
load("p1.FDis.obs.z.lat.TD.lat_post.linpred.re2.rdata")
load("p1.CWM.Bl.TD.elev_post.linpred.re2.rdata")
load("p1.CWM.Bl.TD.lat_post.linpred.re2.rdata")
load("p1.pd.obs.z.ele.TD.elev_post.linpred.re2.rdata")
load("p1.pd.obs.z.lat.TD.lat_post.linpred.re2.rdata")
load("p1.FRic.obs.z.ele.TD.elev_post.linpred.re2.rdata")
load("p1.FRic.obs.z.lat.TD.lat_post.linpred.re2.rdata")

p1.alpha.TD.elev <- ggplotGrob(p1.alpha.TD.elev)
p1.alpha.TD.lat <- ggplotGrob(p1.alpha.TD.lat)
p1.pd.TD.elev <- ggplotGrob(p1.pd.TD.elev)
p1.pd.TD.lat <- ggplotGrob(p1.pd.TD.lat)
p1.mpd.obs.z.ele.TD.elev <- ggplotGrob(p1.mpd.obs.z.ele.TD.elev)
p1.mpd.obs.z.lat.TD.lat <- ggplotGrob(p1.mpd.obs.z.lat.TD.lat)
p1.FRic.TD.elev <- ggplotGrob(p1.FRic.TD.elev)
p1.FRic.TD.lat <- ggplotGrob(p1.FRic.TD.lat)
p1.FDis.obs.z.ele.TD.elev <- ggplotGrob(p1.FDis.obs.z.ele.TD.elev)
p1.FDis.obs.z.lat.TD.lat <- ggplotGrob(p1.FDis.obs.z.lat.TD.lat)
p1.CWM.Bl.TD.elev <- ggplotGrob(p1.CWM.Bl.TD.elev)
p1.CWM.Bl.TD.lat <- ggplotGrob(p1.CWM.Bl.TD.lat)
p1.pd.obs.z.ele.TD.elev <- ggplotGrob(p1.pd.obs.z.ele.TD.elev)
p1.pd.obs.z.lat.TD.lat <- ggplotGrob(p1.pd.obs.z.lat.TD.lat)
p1.FRic.obs.z.ele.TD.elev <- ggplotGrob(p1.FRic.obs.z.ele.TD.elev)
p1.FRic.obs.z.lat.TD.lat <- ggplotGrob(p1.FRic.obs.z.lat.TD.lat)

g12.c1 <- rbind(p1.alpha.TD.elev,
                p1.FRic.TD.elev,
                p1.pd.TD.elev,
                p1.FDis.obs.z.ele.TD.elev)

g12.c1$widths <- unit.pmax(p1.alpha.TD.elev$widths,
                           p1.FRic.TD.elev$widths,
                           p1.pd.TD.elev$widths,
                           p1.FDis.obs.z.ele.TD.elev$widths)

g12.c2 <- rbind(p1.CWM.Bl.TD.elev,
                p1.FRic.obs.z.ele.TD.elev,
                p1.pd.obs.z.ele.TD.elev,
                p1.mpd.obs.z.ele.TD.elev)

g12.c2$widths <- unit.pmax(p1.CWM.Bl.TD.elev$widths,
                           p1.FRic.obs.z.ele.TD.elev$widths,
                           p1.pd.obs.z.ele.TD.elev$widths,
                           p1.mpd.obs.z.ele.TD.elev$widths)

g12.c3 <- rbind(p1.alpha.TD.lat,
                p1.FRic.TD.lat,
                p1.pd.TD.lat,
                p1.FDis.obs.z.lat.TD.lat)

g12.c3$widths <- unit.pmax(p1.alpha.TD.lat$widths,
                           p1.FRic.TD.lat$widths,
                           p1.pd.TD.lat$widths,
                           p1.FDis.obs.z.lat.TD.lat$widths)

g12.c4 <- rbind(p1.CWM.Bl.TD.lat,
                p1.FRic.obs.z.lat.TD.lat,
                p1.pd.obs.z.lat.TD.lat,
                p1.mpd.obs.z.lat.TD.lat)

g12.c4$widths <- unit.pmax(p1.CWM.Bl.TD.lat$widths,
                           p1.FRic.obs.z.lat.TD.lat$widths,
                           p1.pd.obs.z.lat.TD.lat$widths,
                           p1.mpd.obs.z.lat.TD.lat$widths)

g12 <- cbind(g12.c1, g12.c2, g12.c3, g12.c4)
grid.newpage()
grid.draw(g12)

###################################################################################################################
