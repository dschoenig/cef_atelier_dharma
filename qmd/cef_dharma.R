## -----------------------------------------------------------------------------
#| eval: false
# install.packages(c("DHARMa",
#                    "data.table",
#                    "lme4",
#                    "Matrix",
#                    "TMB",
#                    "glmmTMB",
#                    "mgcv",
#                    "mgcViz",
#                    "gratia",
#                    "marginaleffects",
#                    "ggplot2",
#                    "colorspace",
#                    "sf"))


## -----------------------------------------------------------------------------
#| output: false
#| cache: false
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


## -----------------------------------------------------------------------------
#| include: false
#| cache: false
base.family <- "IBMPlexSansCondensed"
base.size <- 11
plot_theme <-
  theme_light(base_size = base.size) +
  theme(
        plot.title = element_text(hjust = 0,
                                  face = "bold",
                                  margin = margin(l = 0, b = base.size/3, t = base.size/3)),
        plot.tag = element_text(face = "bold"),
        axis.line.x = element_line(color = "black",
                                   linewidth = rel(0.5)),
        axis.line.y = element_line(color = "black",
                                   linewidth = rel(0.5)),
        axis.title.x = element_text(margin = margin(t = base.size/2)),
        axis.title.y = element_text(margin = margin(r = base.size/2)),
        # axis.text.y = element_text(color = "black", size = rel(1)),
        # axis.text.x = element_text(color = "black"),
        # legend.title = element_text(margin = margin(b = base.size)),
        legend.position = "right",
        legend.justification = "top",
        legend.key.size = unit(base.size, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # panel.spacing.x = unit(base.size, "pt"),
        # panel.spacing.y = unit(base.size/2, "pt"),
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


## -----------------------------------------------------------------------------
#| eval: false
# base.size <- 11
# plot_theme <-
#   theme_light(base_size = base.size) +
#   theme(plot.title = element_text(hjust = 0,
#                                   face = "bold",
#                                   margin = margin(l = 0, b = base.size/3, t = base.size/3)),
#         plot.tag = element_text(face = "bold"),
#         axis.line.x = element_line(color = "black",
#                                    linewidth = rel(0.5)),
#         axis.line.y = element_line(color = "black",
#                                    linewidth = rel(0.5)),
#         axis.title.x = element_text(margin = margin(t = base.size/2)),
#         axis.title.y = element_text(margin = margin(r = base.size/2)),
#         legend.position = "right",
#         legend.justification = "top",
#         legend.key.size = unit(base.size, "pt"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         plot.margin = margin(3, 3, 3, 3),
#         strip.text = element_text(size = rel(0.8),
#                                   hjust = 0.5,
#                                   color = "black",
#                                   margin = margin(base.size/2,
#                                                   base.size/2,
#                                                   base.size/2,
#                                                   base.size/2)),
#         strip.background = element_rect(fill = "gray90", colour = NA))
# 
# theme_set(plot_theme)


## -----------------------------------------------------------------------------
#| include: false
gutten <- readRDS("../data/gutten.rds")


## -----------------------------------------------------------------------------
#| eval: false
# gutten <- readRDS("gutten.rds")


## -----------------------------------------------------------------------------
gutten


## -----------------------------------------------------------------------------
ggplot(gutten) +
  geom_point(aes(x = age.base, y = volume, colour = quality),
             alpha = 0.5) +
  facet_wrap(vars(quality)) +
  labs(y = "Volume (m³/1000)",
       x = "Age (years)",
       colour = "Site quality")


## -----------------------------------------------------------------------------
mod <- glm(volume ~ age.base * quality,
           family = gaussian(link = "identity"),
           data = gutten)

summary(mod)


## -----------------------------------------------------------------------------
age.seq <- seq(min(gutten$age.base), max(gutten$age.base), length.out = 100)
quality.u <- unique(gutten$quality)

pred.grid <- datagrid(age.base = age.seq,
                      quality = quality.u,
                      model = mod)
mod.pred <- predictions(mod, newdata = pred.grid)

mod.pred


## -----------------------------------------------------------------------------
ggplot(mod.pred) +
  geom_point(data = gutten,
             aes(x = age.base, y = volume, colour = quality),
                 alpha = 0.2) +
  geom_line(aes(x = age.base, y = estimate, colour = quality)) +
  geom_ribbon(aes(x = age.base, ymin = conf.low, ymax = conf.high,
                  group = quality),
              alpha = 0.2) +
  facet_wrap(vars(quality)) +
  labs(y = "Volume (m³/1000)",
       x = "Age (years)",
       colour = "Site quality")


## -----------------------------------------------------------------------------
#| layout-ncol: 2
#| fig-width: 4
#| fig-height: 4
mod.res <- residuals(mod)
mod.fit <- fitted(mod)

ggplot() +
  geom_hline(yintercept = 0, colour = 2) +
  geom_point(aes(y = mod.res, x = mod.fit), alpha = 0.3) +
  labs(x = "Fitted", y = "Residual")

ggplot() +
  geom_hline(yintercept = 0, colour = 2) +
  geom_point(aes(y = mod.res, x = rank(mod.fit)), alpha = 0.3) +
  labs(x = "Fitted (rank)", y = "Residual")


## -----------------------------------------------------------------------------
#| warning: false
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres)


## -----------------------------------------------------------------------------
#| eval: false
# ?simulateResiduals()


## -----------------------------------------------------------------------------
testUniformity(mod.qres)
testDispersion(mod.qres)
testOutliers(mod.qres)


## -----------------------------------------------------------------------------
testQuantiles(mod.qres)


## -----------------------------------------------------------------------------
mod <- glm(volume ~ age.base * quality,
           family = Gamma(link = "log"),
           data = gutten)

summary(mod)


## -----------------------------------------------------------------------------
age.seq <- seq(min(gutten$age.base), max(gutten$age.base), length.out = 100)
quality.u <- unique(gutten$quality)

pred.grid <- datagrid(age.base = age.seq,
                      quality = quality.u,
                      model = mod)
mod.pred <- predictions(mod, newdata = pred.grid)

ggplot(mod.pred) +
  geom_point(data = gutten,
             aes(x = age.base, y = volume, colour = quality),
                 alpha = 0.2) +
  geom_line(aes(x = age.base, y = estimate, colour = quality)) +
  geom_ribbon(aes(x = age.base, ymin = conf.low, ymax = conf.high,
                  group = quality),
              alpha = 0.2) +
  facet_wrap(vars(quality)) +
  labs(y = "Volume (m³/1000)",
       x = "Age (years)",
       colour = "Site quality")


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres)


## -----------------------------------------------------------------------------
#| layout-nrow: 2
plotResiduals(mod.qres, gutten$quality)
plotResiduals(mod.qres, gutten$age.base)


## -----------------------------------------------------------------------------
mod <- glm(volume ~ poly(age.base, 2) * quality,
           family = Gamma(link = "log"),
           data = gutten)

summary(mod)


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres)


## -----------------------------------------------------------------------------
plotResiduals(mod.qres, gutten$age.base)


## -----------------------------------------------------------------------------
mod <- glm(volume ~ poly(age.base, 5) * quality,
           family = Gamma(link = "log"),
           data = gutten)

summary(mod)


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres)


## -----------------------------------------------------------------------------
#| layout-nrow: 2
plotResiduals(mod.qres, gutten$quality)
plotResiduals(mod.qres, gutten$age.base)


## -----------------------------------------------------------------------------
age.seq <- seq(min(gutten$age.base), max(gutten$age.base), length.out = 100)
quality.u <- unique(gutten$quality)

pred.grid <- datagrid(age.base = age.seq,
                      quality = quality.u,
                      model = mod)
mod.pred <- predictions(mod, newdata = pred.grid)

ggplot(mod.pred) +
  geom_point(data = gutten,
             aes(x = age.base, y = volume, colour = quality),
                 alpha = 0.3) +
  geom_line(aes(x = age.base, y = estimate, colour = quality)) +
  geom_ribbon(aes(x = age.base, ymin = conf.low, ymax = conf.high,
                  group = quality),
              alpha = 0.2) +
  facet_wrap(vars(quality)) +
  labs(y = "Volume (m³/1000)",
       x = "Age (years)",
       colour = "Site quality")


## -----------------------------------------------------------------------------

mod <- glmmTMB(volume ~ poly(age.base, 5) * quality + (1 | site + tree.id),
              family = Gamma(link = "log"),
              data = gutten)

summary(mod)


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres)


## -----------------------------------------------------------------------------
age.seq <- seq(min(gutten$age.base), max(gutten$age.base), length.out = 100)
quality.u <- unique(gutten$quality)

pred.grid <- datagrid(age.base = age.seq,
                      quality = quality.u,
                      model = mod)
mod.pred <- predictions(mod, newdata = pred.grid, re.form = NA)

ggplot(mod.pred) +
  geom_point(data = gutten,
             aes(x = age.base, y = volume, colour = quality),
                 alpha = 0.3) +
  geom_line(aes(x = age.base, y = estimate, colour = quality)) +
  geom_ribbon(aes(x = age.base, ymin = conf.low, ymax = conf.high,
                  group = quality),
              alpha = 0.2) +
  facet_wrap(vars(quality)) +
  labs(y = "Volume (m³/1000)",
       x = "Age (years)",
       colour = "Site quality")


## -----------------------------------------------------------------------------
mod <- glmmTMB(volume ~ poly(age.base, 5) * quality + (1 | site + tree.id),
              family = tweedie(link = "log"),
              data = gutten)

summary(mod)
AIC(mod)


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres)


## -----------------------------------------------------------------------------
age.seq <- seq(min(gutten$age.base), max(gutten$age.base), length.out = 100)
quality.u <- unique(gutten$quality)

pred.grid <- datagrid(age.base = age.seq,
                      quality = quality.u,
                      model = mod)
mod.pred <- predictions(mod, newdata = pred.grid, re.form = NA)

ggplot(mod.pred) +
  geom_point(data = gutten,
             aes(x = age.base, y = volume, colour = quality),
                 alpha = 0.3) +
  geom_line(aes(x = age.base, y = estimate, colour = quality)) +
  geom_ribbon(aes(x = age.base, ymin = conf.low, ymax = conf.high,
                  group = quality),
              alpha = 0.2) +
  facet_wrap(vars(quality)) +
  labs(y = "Volume (m³/1000)",
       x = "Age (years)",
       colour = "Site quality")


## -----------------------------------------------------------------------------
mod <- glmmTMB(volume ~ poly(age.base, 5) * quality + (1 | site + tree.id),
               dispformula = ~ age.base * quality,
               family = Gamma(link = "log"),
               data = gutten)

summary(mod)


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres)


## -----------------------------------------------------------------------------
age.seq <- seq(min(gutten$age.base), max(gutten$age.base), length.out = 100)
quality.u <- unique(gutten$quality)

pred.grid <- datagrid(age.base = age.seq,
                      quality = quality.u,
                      model = mod)
mod.pred <- predictions(mod, newdata = pred.grid, re.form = NA)

ggplot(mod.pred) +
  geom_point(data = gutten,
             aes(x = age.base, y = volume, colour = quality),
                 alpha = 0.3) +
  geom_line(aes(x = age.base, y = estimate, colour = quality)) +
  geom_ribbon(aes(x = age.base, ymin = conf.low, ymax = conf.high,
                  group = quality),
              alpha = 0.2) +
  facet_wrap(vars(quality)) +
  labs(y = "Volume (m³/1000)",
       x = "Age (years)",
       colour = "Site quality")


## -----------------------------------------------------------------------------
#| include: false
sitka <- readRDS("../data/sitka.rds")


## -----------------------------------------------------------------------------
#| eval: false
# sitka <- readRDS("sitka.rds")


## -----------------------------------------------------------------------------
sitka


## -----------------------------------------------------------------------------
ggplot(sitka) +
  geom_line(data = sitka,
            aes(x = day, y = size,
                group = tree.id, colour = treatment),
            alpha = 0.5) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  labs(x = "Time (days)",
       y = "Size (cm²m)",
       colour = "Treatment")


## -----------------------------------------------------------------------------
mod <- glmmTMB(size ~ day * treatment + (1 | tree.id),
               data = sitka,
               family = Gamma(link = "log"))

summary(mod)


## -----------------------------------------------------------------------------
day.seq <- seq(min(sitka$day), max(sitka$day), length.out = 100)
treatment.u <- unique(sitka$treatment)

pred.grid <- datagrid(day = day.seq,
                      treatment = treatment.u,
                      model = mod)
mod.pred <- predictions(mod, newdata = pred.grid, re.form = NA)

ggplot(mod.pred) +
  geom_line(data = sitka, aes(x = day, y = size, group = tree.id),
            alpha = 0.2) +
  geom_ribbon(aes(x = day, ymin = conf.low, ymax = conf.high,
                  fill = treatment),
              alpha = 0.2) +
  geom_line(aes(x = day, y = estimate, colour = treatment)) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  scale_fill_brewer(type = "qual", palette = "Set1",
                    guide = "none") +
  labs(x = "Time (days)",
       y = "Size (cm²m)",
       colour = "Treatment")



## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres)


## -----------------------------------------------------------------------------
#| warning: false
#| fig-width: 10
#| fig-height: 5.85
mod.qres.time <- recalculateResiduals(mod.qres, group = sitka$day)

testTemporalAutocorrelation(residuals(mod.qres.time), unique(sitka$day))


## -----------------------------------------------------------------------------
mod <- gam(size ~ s(day) + treatment + s(tree.id, bs = "re"),
           data = sitka,
           family = Gamma(link = "log"),
           method = "REML")

summary(mod)


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres)


## -----------------------------------------------------------------------------
day.seq <- seq(min(sitka$day), max(sitka$day), length.out = 100)
treatment.u <- unique(sitka$treatment)

pred.grid <- datagrid(day = day.seq,
                      treatment = treatment.u,
                      model = mod)
mod.pred <- predictions(mod,
                        newdata = pred.grid,
                        exclude = "s(tree.id)")

ggplot(mod.pred) +
  geom_line(data = sitka, aes(x = day, y = size, group = tree.id),
            alpha = 0.2) +
  geom_ribbon(aes(x = day, ymin = conf.low, ymax = conf.high,
                  fill = treatment),
              alpha = 0.2) +
  geom_line(aes(x = day, y = estimate, colour = treatment)) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  scale_fill_brewer(type = "qual", palette = "Set1",
                    guide = "none") +
  labs(x = "Time (days)",
       y = "Size (cm²m)",
       colour = "Treatment")


## -----------------------------------------------------------------------------
#| warning: false
#| fig-width: 10
#| fig-height: 5.85
mod.qres.time <- recalculateResiduals(mod.qres, group = sitka$day)

testTemporalAutocorrelation(residuals(mod.qres.time), unique(sitka$day))


## -----------------------------------------------------------------------------
mod <- gam(size ~
             s(day) +
             s(day, treatment, bs = "sz") +
             s(tree.id, bs = "re"),
           data = sitka,
           family = Gamma(link = "log"),
           method = "REML")

summary(mod)


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres)


## -----------------------------------------------------------------------------
day.seq <- seq(min(sitka$day), max(sitka$day), length.out = 100)
treatment.u <- unique(sitka$treatment)

pred.grid <- datagrid(day = day.seq,
                      treatment = treatment.u,
                      model = mod)
mod.pred <- predictions(mod,
                        newdata = pred.grid,
                        exclude = "s(tree.id)")

ggplot(mod.pred) +
  geom_line(data = sitka, aes(x = day, y = size, group = tree.id),
            alpha = 0.2) +
  geom_ribbon(aes(x = day, ymin = conf.low, ymax = conf.high,
                  fill = treatment),
              alpha = 0.2) +
  geom_line(aes(x = day, y = estimate, colour = treatment)) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  scale_fill_brewer(type = "qual", palette = "Set1",
                    guide = "none") +
  labs(x = "Time (days)",
       y = "Size (cm²m)",
       colour = "Treatment")


## -----------------------------------------------------------------------------
#| warning: false
#| fig-width: 10
#| fig-height: 5.85
mod.qres.time <- recalculateResiduals(mod.qres, group = sitka$day)

testTemporalAutocorrelation(residuals(mod.qres.time), unique(sitka$day))


## -----------------------------------------------------------------------------
gratia::draw(mod)


## -----------------------------------------------------------------------------
#| include: false
lichen <- readRDS("../data/lichen.rds")
swe <- st_read("../data/adm_swe.gpkg")


## -----------------------------------------------------------------------------
#| eval: false
# lichen <- readRDS("lichen.rds")
# swe <- st_read("adm_swe.gpkg")


## -----------------------------------------------------------------------------
lichen


## -----------------------------------------------------------------------------
lichen.usnea1 <- lichen[species == "Usnea" & ip == 1]


## -----------------------------------------------------------------------------
plot_theme <- 
  plot_theme +
  theme(panel.grid.major = element_line(linewidth = rel(1)))

theme_set(plot_theme)


## -----------------------------------------------------------------------------
#| fig-width: 8
#| fig-height: 10
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


## -----------------------------------------------------------------------------

mod <-
  glm(occurrence ~ temp * rain,
      data = lichen.usnea1,
      family = binomial(link = "logit"))

summary(mod)


## -----------------------------------------------------------------------------
mod.pred <- avg_predictions(mod, by = c("east.agg", "north.agg"))

ggplot(mod.pred) +
  geom_raster(aes(x = east.agg, y = north.agg, fill = estimate)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL, fill = "Estimated occcurence\nprobability")


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres, quantreg = TRUE)


## -----------------------------------------------------------------------------
#| warning: false
#| fig-width: 10
#| fig-height: 5.85
mod.res.agg <- recalculateResiduals(mod.qres, lichen.usnea1$rast.agg.id)

plot(mod.res.agg, quantreg = TRUE)


## -----------------------------------------------------------------------------
testSpatialAutocorrelation(mod.qres, lichen.usnea1$east, lichen.usnea1$north)


## -----------------------------------------------------------------------------
lichen.usnea1[, qres := residuals(mod.qres)]


lichen.usnea1[,
              .(qres = median(qres)),
              by = c("east.agg", "north.agg")] |>
ggplot() +
  geom_raster(aes(x = east.agg, y = north.agg, fill = qres)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_continuous_diverging("Blue-Red", mid = 0.5) +
  labs(x = NULL, y = NULL, fill = "Quantile residual (median)")



## -----------------------------------------------------------------------------
#| layout-ncol: 2
#| fig-width: 4
#| fig-height: 4
plotResiduals(mod.qres, lichen.usnea1$temp, quantreg = TRUE)
plotResiduals(mod.qres, lichen.usnea1$rain, quantreg = TRUE)


## -----------------------------------------------------------------------------
#| layout-ncol: 2
#| fig-width: 4
#| fig-height: 4
plotResiduals(mod.qres, lichen.usnea1$ndep, quantreg = TRUE)
plotResiduals(mod.qres, lichen.usnea1$mat, quantreg = TRUE)


## -----------------------------------------------------------------------------
#| layout-ncol: 2
#| fig-width: 4
#| fig-height: 4
plotResiduals(mod.qres, lichen.usnea1$dbh, quantreg = TRUE)
plotResiduals(mod.qres, lichen.usnea1$crl, quantreg = TRUE)


## -----------------------------------------------------------------------------
#| layout-ncol: 2
#| fig-width: 4
#| fig-height: 4
plotResiduals(mod.qres, lichen.usnea1$bas, quantreg = TRUE)
plotResiduals(mod.qres, lichen.usnea1$age, quantreg = TRUE)


## -----------------------------------------------------------------------------

mod <-
  glm(occurrence ~
        temp * rain +
        mat + ndep +
        dbh + crl +
        bas + age,
      data = lichen.usnea1,
      family = binomial(link = "logit"))

summary(mod)


## -----------------------------------------------------------------------------
mod.pred <- avg_predictions(mod, by = c("east.agg", "north.agg"))

ggplot(mod.pred) +
  geom_raster(aes(x = east.agg, y = north.agg, fill = estimate)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL, fill = "Estimated occcurence\nprobability")


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres, quantreg = TRUE)


## -----------------------------------------------------------------------------
#| warning: false
#| fig-width: 10
#| fig-height: 5.85
mod.res.agg <- recalculateResiduals(mod.qres, lichen.usnea1$rast.agg.id)

plot(mod.res.agg, quantreg = TRUE)


## -----------------------------------------------------------------------------
testSpatialAutocorrelation(mod.qres, lichen.usnea1$east, lichen.usnea1$north)


## -----------------------------------------------------------------------------
lichen.usnea1[, qres := residuals(mod.qres)]

lichen.usnea1[,
              .(qres = median(qres)),
              by = c("east.agg", "north.agg")] |>
ggplot() +
  geom_raster(aes(x = east.agg, y = north.agg, fill = qres)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_continuous_diverging("Blue-Red", mid = 0.5) +
  labs(x = NULL, y = NULL, fill = "Quantile residual")


## -----------------------------------------------------------------------------

mod <-
  gam(occurrence ~
        s(east, north, k = 60) +
        temp * rain +
        mat + ndep +
        dbh + crl +
        bas + age,
      data = lichen.usnea1,
      family = binomial(link = "logit"),
      method = "REML")

summary(mod)


## -----------------------------------------------------------------------------
mod.pred <- avg_predictions(mod, by = c("east.agg", "north.agg"))

ggplot(mod.pred) +
  geom_raster(aes(x = east.agg, y = north.agg, fill = estimate)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL, fill = "Estimated occcurence\nprobability")


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres, quantreg = TRUE)


## -----------------------------------------------------------------------------
#| warning: false
#| fig-width: 10
#| fig-height: 5.85
mod.res.agg <- recalculateResiduals(mod.qres, lichen.usnea1$rast.agg.id)

plot(mod.res.agg, quantreg = TRUE)


## -----------------------------------------------------------------------------
testSpatialAutocorrelation(mod.qres, lichen.usnea1$east, lichen.usnea1$north)


## -----------------------------------------------------------------------------
lichen.usnea1[, qres := residuals(mod.qres)]

lichen.usnea1[,
              .(qres = median(qres)),
              by = c("east.agg", "north.agg")] |>
ggplot() +
  geom_raster(aes(x = east.agg, y = north.agg, fill = qres)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_continuous_diverging("Blue-Red", mid = 0.5) +
  labs(x = NULL, y = NULL, fill = "Quantile residual")


## -----------------------------------------------------------------------------

lichen.usnea <- lichen[species == "Usnea"]

mod <-
  gam(occurrence ~
        s(east, north, by = ip, k = 60) +
        temp * rain +
        mat + ndep +
        dbh + crl +
        bas + age,
      data = lichen.usnea,
      family = binomial(link = "logit"),
      method = "REML")

summary(mod)


## -----------------------------------------------------------------------------
mod.pred <- avg_predictions(mod, by = c("east.agg", "north.agg", "ip"))

ggplot(mod.pred) +
  geom_raster(aes(x = east.agg, y = north.agg, fill = estimate)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_viridis_c() +
  facet_wrap(vars(ip)) +
  labs(x = NULL, y = NULL, fill = "Estimated occcurence\nprobability")


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres, quantreg = TRUE)


## -----------------------------------------------------------------------------
#| warning: false
#| fig-width: 10
#| fig-height: 5.85
mod.res.agg <- recalculateResiduals(mod.qres, lichen.usnea$rast.agg.id)

plot(mod.res.agg, quantreg = TRUE)


## -----------------------------------------------------------------------------
lichen.usnea[, qres := residuals(mod.qres)]

lichen.usnea[,
             .(qres = median(qres)),
             by = c("ip", "east.agg", "north.agg")] |>
ggplot() +
  geom_raster(aes(x = east.agg, y = north.agg, fill = qres)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_continuous_diverging("Blue-Red", mid = 0.5) +
  facet_wrap(vars(ip)) +
  labs(x = NULL, y = NULL, fill = "Quantile residual")


## -----------------------------------------------------------------------------

mod <-
  gam(occurrence ~
        s(east, north, by = ip, k = 60) +
        ti(temp, k = 5) +
        ti(rain, k = 5) +
        ti(temp, rain, k = c(5, 5)) +
        s(mat) +
        s(ndep) +
        s(dbh) +
        s(crl) +
        s(bas) +
        s(age),
      data = lichen.usnea,
      family = binomial(link = "logit"),
      method = "REML",
      select = TRUE,
      optimizer = "efs")

summary(mod)


## -----------------------------------------------------------------------------
mod.pred <- avg_predictions(mod, by = c("east.agg", "north.agg", "ip"))

ggplot(mod.pred) +
  geom_raster(aes(x = east.agg, y = north.agg, fill = estimate)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_viridis_c() +
  facet_wrap(vars(ip)) +
  labs(x = NULL, y = NULL, fill = "Estimated occcurence\nprobability")


## -----------------------------------------------------------------------------
#| fig-width: 10
#| fig-height: 5.85
mod.qres <- simulateResiduals(mod)

plot(mod.qres, quantreg = TRUE)


## -----------------------------------------------------------------------------
#| warning: false
#| fig-width: 10
#| fig-height: 5.85
mod.res.agg <- recalculateResiduals(mod.qres, lichen.usnea$rast.agg.id)

plot(mod.res.agg, quantreg = TRUE)


## -----------------------------------------------------------------------------
lichen.usnea[, qres := residuals(mod.qres)]

lichen.usnea[,
             .(qres = median(qres)),
             by = c("ip", "east.agg", "north.agg")] |>
ggplot() +
  geom_raster(aes(x = east.agg, y = north.agg, fill = qres)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_continuous_diverging("Blue-Red", mid = 0.5) +
  facet_wrap(vars(ip)) +
  labs(x = NULL, y = NULL, fill = "Quantile residual")


## -----------------------------------------------------------------------------
gratia::draw(mod, rug = FALSE, select = 1:2)


## -----------------------------------------------------------------------------
gratia::draw(mod, rug = FALSE, select = 3:5)


## -----------------------------------------------------------------------------
#| fig-height: 3
gratia::draw(mod, rug = FALSE, select = 6:7)


## -----------------------------------------------------------------------------
#| fig-height: 3
gratia::draw(mod, rug = FALSE, select = 8:9)


## -----------------------------------------------------------------------------
#| fig-height: 3
gratia::draw(mod, rug = FALSE, select = 10:11)


## -----------------------------------------------------------------------------
sessionInfo()


## -----------------------------------------------------------------------------
#| include: false
knitr::purl(input = "cef_dharma.qmd", output = "cef_dharma.R")

