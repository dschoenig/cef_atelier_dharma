#####################################################################
# 0. PRÉPARATION (AVANT L'ATELIER) ##################################
#####################################################################

# Avant d'installer les paquets, il faut s'assurer de disposer de la
# version la plus récente de R (https://cran.r-project.org/)

install.packages(c("DHARMa",
                   "data.table",
                   "lme4",
                   "Matrix",
                   "TMB",
                   "glmmTMB",
                   "mgcv",
                   "mgcViz",
                   "gratia",
                   "marginaleffects",
                   "ggplot2",
                   "colorspace",
                   "sf"))


#####################################################################
# 1. CONFIGURATION ##################################################
#####################################################################


# Charger les paquets

library(DHARMa)
library(data.table)
library(lme4)
library(glmmTMB)
library(mgcv)
library(marginaleffects)
library(ggplot2)
library(colorspace)
library(sf)
options(show.signif.stars = FALSE)


# Définir un thème pour ggplot2

base.size <- 11
plot_theme <-
  theme_light(base_size = base.size) +
  theme(plot.title = element_text(hjust = 0,
                                  face = "bold",
                                  margin = margin(l = 0, b = base.size/3, t = base.size/3)),
        plot.tag = element_text(face = "bold"),
        axis.line.x = element_line(color = "black",
                                   linewidth = rel(0.5)),
        axis.line.y = element_line(color = "black",
                                   linewidth = rel(0.5)),
        axis.title.x = element_text(margin = margin(t = base.size/2)),
        axis.title.y = element_text(margin = margin(r = base.size/2)),
        legend.position = "right",
        legend.justification = "top",
        legend.key.size = unit(base.size, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(3, 3, 3, 3),
        strip.text = element_text(size = rel(0.8),
                                  hjust = 0.5,
                                  color = "black",
                                  margin = margin(base.size/2,
                                                  base.size/2,
                                                  base.size/2,
                                                  base.size/2)),
        strip.background = element_rect(fill = "gray90", colour = NA))
theme_set(plot_theme)


# Spécifier le répositoire de travail. Dans ce cas on utilise le dossier
# qui contient les données (remplacer DOSSIER par le chemin d'accès au
# répertoire de données)

setwd("DOSSIER")



#####################################################################
# Exemple 1: Croissance de'l Épinette de Norvège dans les Alpes #####
#####################################################################


# Importer les données

gutten <- readRDS("gutten.rds")


# Aperçu

gutten


# Graphique exploratoire

ggplot(gutten) +
  geom_point(aes(x = age.base, y = volume, colour = quality),
             alpha = 0.5) +
  facet_wrap(vars(quality)) +
  labs(y = "Volume (m³/1000)",
       x = "Age (years)",
       colour = "Site quality")



#####################################################################
# Exemple 2: Effet de l'ozone sur les semis de l'épinette de Sitka ##
#####################################################################


# Importer les données

sitka <- readRDS("sitka.rds")


# Aperçu

sitka


# Graphique exploratoire

ggplot(sitka) +
  geom_line(data = sitka,
            aes(x = day, y = size,
                group = tree.id, colour = treatment),
            alpha = 0.5) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  labs(x = "Time (days)",
       y = "Size (cm²m)",
       colour = "Treatment")



#####################################################################
# Exemple 3: Distribution des lichens en Suède ######################
#####################################################################


# Importer les données

lichen <- readRDS("lichen.rds")
swe <- st_read("adm_swe.gpkg")


# Aperçu

lichen


# Modifier le thème pour ggplot2

plot_theme <- 
  plot_theme +
  theme(panel.grid.major = element_line(linewidth = rel(1)))
theme_set(plot_theme)


# Graphique exploratoire

lichen.prop.agg <-
  lichen[,
         .(occurrence = mean(occurrence)),
         by = .(ip, species, east.agg, north.agg)]

ggplot(lichen.prop.agg) +
  geom_raster(aes(x = east.agg, y = north.agg, fill = occurrence)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(ip), cols = vars(species)) +
  labs(x = NULL, y = NULL, fill = "Occurence probability")
