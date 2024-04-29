library(DHARMa)
library(data.table)
library(lme4)
library(glmmTMB)
library(mgcv)
library(marginaleffects)
library(ggplot2)
library(colorspace)
library(sf)
library(patchwork)




# DATA OVERVIEWS -------------------------------------------------------

base.family <- "IBMPlexSansCondensed"
base.size <- 12
plot_theme <-
  theme_light(base_size = base.size, base_family = base.family) +
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


gutten <- readRDS("../data/gutten.rds")

png("../qmd/figures/gutten.png", width = 8, height = 5, units = "in", res = 150)
ggplot(gutten) +
  geom_point(aes(x = age.base, y = volume, colour = quality),
             alpha = 0.5) +
  facet_wrap(vars(quality)) +
  labs(y = "Volume (m³/1000)",
       x = "Age (years)",
       colour = "Site quality")
dev.off()


sitka <- readRDS("../data/sitka.rds")

png("../qmd/figures/sitka.png", width = 8, height = 5, units = "in", res = 150)
ggplot(sitka) +
  geom_line(data = sitka,
            aes(x = day, y = size,
                group = tree.id, colour = treatment),
            alpha = 0.5) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  labs(x = "Time (days)",
       y = "Size (cm²m)",
       colour = "Treatment")
dev.off()


lichen <- readRDS("../data/lichen.rds")
swe <- st_read("../data/adm_swe.gpkg")

png("../qmd/figures/lichen.png", width = 12, height = 5, units = "in", res = 150)
lichen.agg <-
  lichen[,
       .(occurrence = mean(occurrence)),
       by = .(ip, species, east.agg, north.agg)]
p1 <-
ggplot(lichen.agg[ip == 1]) +
  geom_raster(aes(x = east.agg, y = north.agg, fill = occurrence)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_viridis_c() +
  facet_wrap(vars(species), nrow = 1) +
  labs(x = NULL, y = NULL, fill = "Occurence\nprobability") +
  labs(title = "First inventory period")
p2 <-
ggplot(lichen.agg[ip == 2]) +
  geom_raster(aes(x = east.agg, y = north.agg, fill = occurrence)) +
  geom_sf(data = swe, fill = NA, colour = "black") +
  scale_fill_viridis_c() +
  facet_wrap(vars(species), nrow = 1) +
  labs(x = NULL, y = NULL, fill = "Occurence\nprobability") +
  labs(title = "Second inventory period") +
  theme(axis.text.y = element_blank())
p1 + p2 + plot_layout(guides = "collect") & 
  theme(panel.grid.major = element_line(linewidth = rel(1)),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
dev.off()



# RESIDUALS ------------------------------------------------------------


base.family <- "IBMPlexSansCondensed"
base.size <- 12
plot_theme <-
  theme_light(base_size = base.size, base_family = base.family) +
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



# Normal distribution, first focal point

set.seed(123)
n <- 25
a <- 10
b <- 0.4
s <- 10
x <- sample(1:100, n)
mu <- a + b * x
y <- rnorm(n, mean = mu, sd = s)
dat <- data.frame(x = x, y = y)
mod <- lm(y ~ x, data = dat)
mod.c <- coef(mod)
mod.s <- sigma(mod)
dat$mu <- fitted(mod)
dat$r <- residuals(mod)
qres <- simulateResiduals(mod, n = 500)
rq <- residuals(qres)
x.foc <- 57
y.foc <- y[which(x == x.foc)]
mu.foc <- fitted(mod)[which(x == x.foc)]
r.foc <- residuals(mod)[which(x == x.foc)]
seg.dat <- data.frame(x = x, mu = fitted(mod), y = y)
y.d <- seq(mu.foc-100, mu.foc+100, length.out = 501)
y.d <- c(y.d, y.d[1])
x.d <- dnorm(y.d, mean = mu.foc, sd = mod.s)
d.dat <- data.frame(x.d = x.d, y.d = y.d, r.d = y.d-mu.foc)
y.foc.d <- seq(mu.foc-100, y.foc, length.out = 501)
x.foc.d <- dnorm(y.foc.d, mean = mu.foc, sd = mod.s)
y.foc.d <- c(y.foc.d, y.foc.d[length(y.foc.d)], 0)
x.foc.d <- c(x.foc.d, x.foc.d[1], x.foc.d[1])
d.foc.dat <- data.frame(x.foc.d = x.foc.d, y.foc.d = y.foc.d)
pal <- c("#e41a1c", "#377eb8", "#984ea3")
lab <-
  data.frame(x = c(x.foc+2, x.foc+2, x.foc-2),
             y = c(y.foc+2, mu.foc-2, (y.foc+mu.foc)/2),
             lab1 = c("y", "µ", "r"),
             lab2 = c("observée", "prédite", "résidu"),
             hjust = c(0, 0, 1),
             col = as.character(1:3))
col.txt <- c("black", pal[2], pal[1])
names(col.txt) <- lab$col
size.pt.small <- 2 
size.pt.large <- 3
stroke.foc <- 0.2
lw.res <- 0.75

plot_res <- function(x, file) {
  png(file, width = 5.5, height = 4, units = "in", res = 300)
  print(x)
  dev.off()
}

plot_d <- function(x, file) {
  png(file, width = 2, height = 4, units = "in", res = 300)
  print(x)
  dev.off()
}

(ggplot(dat) +
   geom_point(aes(x = x, y = y),
              size = size.pt.small, colour = "grey5") +
   geom_abline(intercept = mod.c[1], slope = mod.c[2],
               colour = pal[2]) +
   scale_colour_manual(values = col.txt, guide = "none") +
   coord_cartesian(ylim = c(0, 65)) +
   labs(y = "Réponse", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.lm.1.png")




(ggplot(dat) +
   geom_point(aes(x = x, y = y),
              size = size.pt.small, colour = "grey50") +
   geom_abline(intercept = mod.c[1], slope = mod.c[2],
               colour = pal[2]) +
   geom_point(x = x.foc, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = x.foc, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   scale_colour_manual(values = col.txt, guide = "none") +
   coord_cartesian(ylim = c(0, 65)) +
   labs(y = "Réponse", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.lm.2.png")



(ggplot(dat) +
   geom_point(aes(x = x, y = y),
              size = size.pt.small, colour = "grey50") +
   geom_abline(intercept = mod.c[1], slope = mod.c[2],
               colour = pal[2]) +
   geom_vline(xintercept = x.foc, linetype = "dashed", colour = "grey5") +
   geom_point(x = x.foc, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = x.foc, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = x.foc, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   scale_colour_manual(values = col.txt, guide = "none") +
   coord_cartesian(ylim = c(0, 65)) +
   labs(y = "Réponse", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.lm.3.png")
 
 
(ggplot(dat) +
   geom_point(aes(x = x, y = y),
              size = size.pt.small, colour = "grey50") +
   geom_abline(intercept = mod.c[1], slope = mod.c[2],
               colour = pal[2]) +
   geom_vline(xintercept = x.foc, linetype = "dashed", colour = "grey5") +
   geom_segment(x = x.foc, xend = x.foc, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = x.foc, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = x.foc, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = x.foc, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_text(data = lab, aes(x = x, y = y, label = lab2, hjust = hjust, colour = col),
             family = base.family, fontface = "italic", size = 5) +
   scale_colour_manual(values = col.txt, guide = "none") +
   coord_cartesian(ylim = c(0, 65)) +
   labs(y = "Réponse", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.lm.4.png")


(ggplot(dat) +
   geom_point(aes(x = x, y = y),
              size = size.pt.small, colour = "grey50") +
   geom_abline(intercept = mod.c[1], slope = mod.c[2],
               colour = pal[2]) +
   geom_vline(xintercept = x.foc, linetype = "dashed", colour = "grey5") +
   geom_segment(x = x.foc, xend = x.foc, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = x.foc, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = x.foc, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = x.foc, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_text(data = lab, aes(x = x, y = y, label = lab1, hjust = hjust, colour = col),
             family = base.family, fontface = "italic", size = 5) +
   scale_colour_manual(values = col.txt, guide = "none") +
   coord_cartesian(ylim = c(0, 65)) +
   labs(y = "Réponse", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.lm.5.png")


(ggplot(dat) +
   geom_point(aes(x = x, y = y),
              size = size.pt.small, colour = "grey50") +
   geom_abline(intercept = mod.c[1], slope = mod.c[2],
               colour = pal[2]) +
   geom_vline(xintercept = x.foc, linetype = "dashed", colour = "grey5") +
   geom_segment(aes(x = x, xend = x, y = mu, yend = y),
                colour = pal[1], linewidth = lw.res, alpha = 0.5) +
   geom_segment(x = x.foc, xend = x.foc, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = x.foc, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = x.foc, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = x.foc, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_text(data = lab, aes(x = x, y = y, label = lab1, hjust = hjust, colour = col),
             family = base.family, fontface = "italic", size = 5) +
   scale_colour_manual(values = col.txt, guide = "none") +
   coord_cartesian(ylim = c(0, 65)) +
   labs(y = "Réponse", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.lm.6.png")


(ggplot(dat) +
  geom_point(aes(x = x, y = r),
             size = size.pt.small, colour = "grey50") +
  geom_hline(yintercept = 0, colour = pal[2]) +
  geom_vline(xintercept = x.foc, linetype = "dashed", colour = "grey5") +
  geom_segment(aes(x = x, xend = x, y = 0, yend = r),
               colour = pal[1], linewidth = lw.res, alpha = 0.5) +
  geom_segment(x = x.foc, xend = x.foc, y = 0, yend = r.foc,
               colour = pal[1], linewidth = lw.res) +
  geom_point(x = x.foc, y = r.foc, fill = 1,
             shape = 21, size = size.pt.large, stroke = stroke.foc) +
  geom_point(x = x.foc, y = r.foc, colour = "white",
             shape = 16, size = size.pt.large/2) +
  scale_colour_manual(values = col.txt, guide = "none") +
  coord_cartesian(ylim = c(0, 65)-mu.foc) +
  labs(y = "Résidu", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.lm.7.png")

(ggplot(dat) +
  geom_point(aes(x = x, y = r),
             size = size.pt.small, colour = "grey5") +
  coord_cartesian(ylim = c(0, 65)-mu.foc) +
  labs(y = "Résidu", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.lm.8.png")

(ggplot(dat) +
  geom_point(aes(x = x, y = rq),
             size = size.pt.small, colour = pal[3]) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Résidu quantile", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.lm.9.png")

 
(ggplot(dat) +
   geom_polygon(data = d.dat,
                aes(x = x.d, y = y.d), fill = "grey75", colour = NA) +
   geom_path(data = d.dat[-nrow(d.dat),],
             aes(x = x.d, y = y.d),
             colour = "grey35", linewidth = 0.3) +
   geom_vline(xintercept = -0.00375, linetype = "dashed", colour = "grey5") +
   geom_segment(x = -0.00375, xend = -0.00375, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = -0.00375, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = -0.00375, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = -0.00375, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   scale_x_continuous(breaks = c(0, 0.05)) +
   coord_cartesian(xlim = c(-0.01, 0.05), ylim = c(0, 65)) +
   theme(axis.text.x = element_text(color = "white")) +
   labs(x = "Densité", y = "Réponse")) |>
plot_d("../qmd/figures/res.d.1.png")
 
 
 

(ggplot(dat) +
   geom_polygon(data = d.dat,
                aes(x = x.d, y = y.d), fill = "grey75", colour = NA) +
   geom_polygon(data = d.foc.dat,
                aes(x = x.foc.d, y = y.foc.d),
                fill = "white", colour = NA) +
   geom_polygon(data = d.foc.dat,
                aes(x = x.foc.d, y = y.foc.d),
                fill = pal[3], colour = NA, alpha = 0.8) +
   geom_path(data = d.dat[-nrow(d.dat),],
             aes(x = x.d, y = y.d),
             colour = "grey35", linewidth = 0.3) +
   geom_vline(xintercept = -0.00375, linetype = "dashed", colour = "grey5") +
   geom_segment(x = -0.00375, xend = -0.00375, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = -0.00375, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = -0.00375, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = -0.00375, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   scale_x_continuous(breaks = c(0, 0.05)) +
   coord_cartesian(xlim = c(-0.01, 0.05), ylim = c(0, 65)) +
   theme(axis.text.x = element_text(color = "white")) +
   labs(x = "Densité", y = "Réponse")) |>
plot_d("../qmd/figures/res.d.2.png")


(ggplot(dat) +
   # geom_polygon(data = d.dat,
   #              aes(x = x.d, y = y.d), fill = "grey75", colour = NA) +
   geom_polygon(data = d.foc.dat,
                aes(x = x.foc.d, y = y.foc.d),
                fill = "white", colour = NA) +
   geom_polygon(data = d.foc.dat,
                aes(x = x.foc.d, y = y.foc.d),
                fill = pal[3], colour = NA, alpha = 0.8) +
   geom_path(data = d.dat[-nrow(d.dat),],
             aes(x = x.d, y = y.d),
             colour = "grey35", linewidth = 0.3) +
   geom_vline(xintercept = -0.00375, linetype = "dashed", colour = "grey5") +
   geom_segment(x = -0.00375, xend = -0.00375, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = -0.00375, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = -0.00375, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = -0.00375, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_text(data = data.frame(x = x.foc.d[length(x.foc.d)-2],
                               y = y.foc+(sign(y.foc) * 2),
                               lab = "résidu\nquantile"),
             aes(x = x, y = y, label = lab, vjust = vjust), colour = pal[3],
             family = base.family, fontface = "italic",
             hjust = 0, vjust = 0, size = 4) +
   scale_x_continuous(breaks = c(0, 0.05)) +
   coord_cartesian(xlim = c(-0.01, 0.05), ylim = c(0, 65)) +
   theme(axis.text.x = element_text(color = "white")) +
   labs(x = "Densité", y = "Réponse")) |>
plot_d("../qmd/figures/res.d.3.png")


(ggplot(dat) +
   geom_polygon(data = d.dat,
                aes(x = x.d, y = r.d), fill = "grey75", colour = NA) +
   geom_path(data = d.dat[-nrow(d.dat),],
             aes(x = x.d, y = r.d),
             colour = "grey35", linewidth = 0.3) +
   geom_vline(xintercept = -0.00375, linetype = "dashed", colour = "grey5") +
   geom_segment(x = -0.00375, xend = -0.00375, y = 0, yend = r.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = -0.00375, y = r.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = -0.00375, y = r.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = -0.00375, y = 0, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   scale_x_continuous(breaks = c(0, 0.05)) +
   coord_cartesian(xlim = c(-0.01, 0.05), ylim = c(0, 65)-mu.foc) +
   theme(axis.text.x = element_text(color = "white")) +
   labs(x = "Densité", y = "Résidu")) |>
plot_d("../qmd/figures/res.d.4.png")



# Normal distribution, second focal point

set.seed(123)
n <- 25
a <- 10
b <- 0.4
s <- 10
x <- sample(1:100, n)
mu <- a + b * x
y <- rnorm(n, mean = mu, sd = s)
dat <- data.frame(x = x, y = y)
mod <- lm(y ~ x, data = dat)
mod.c <- coef(mod)
mod.s <- sigma(mod)
dat$mu <- fitted(mod)
dat$r <- residuals(mod)
x.foc <- 32
y.foc <- y[which(x == x.foc)]
mu.foc <- fitted(mod)[which(x == x.foc)]
r.foc <- residuals(mod)[which(x == x.foc)]
seg.dat <- data.frame(x = x, mu = fitted(mod), y = y)
y.d <- seq(mu.foc-100, mu.foc+100, length.out = 501)
y.d <- c(y.d, y.d[1])
x.d <- dnorm(y.d, mean = mu.foc, sd = mod.s)
d.dat <- data.frame(x.d = x.d, y.d = y.d, r.d = y.d-mu.foc)
y.foc.d <- seq(mu.foc-100, y.foc, length.out = 501)
x.foc.d <- dnorm(y.foc.d, mean = mu.foc, sd = mod.s)
y.foc.d <- c(y.foc.d, y.foc.d[length(y.foc.d)], 0)
x.foc.d <- c(x.foc.d, x.foc.d[1], x.foc.d[1])
d.foc.dat <- data.frame(x.foc.d = x.foc.d, y.foc.d = y.foc.d)
pal <- c("#e41a1c", "#377eb8", "#984ea3")
lab <-
  data.frame(x = c(x.foc+2, x.foc+2, x.foc-2),
             y = c(y.foc-2, mu.foc+4, (y.foc+mu.foc)/2),
             lab1 = c("y", "µ", "r"),
             lab2 = c("observée", "prédite", "résidu"),
             hjust = c(0, 0, 1),
             col = as.character(1:3))
col.txt <- c("black", pal[2], pal[1])
names(col.txt) <- lab$col
size.pt.small <- 2 
size.pt.large <- 3
stroke.foc <- 0.2
lw.res <- 0.75


(ggplot(dat) +
   geom_point(aes(x = x, y = y),
              size = size.pt.small, colour = "grey50") +
   geom_abline(intercept = mod.c[1], slope = mod.c[2],
               colour = pal[2]) +
   geom_vline(xintercept = x.foc, linetype = "dashed", colour = "grey5") +
   geom_segment(x = x.foc, xend = x.foc, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = x.foc, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = x.foc, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = x.foc, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_text(data = lab, aes(x = x, y = y, label = lab1, hjust = hjust, colour = col),
             family = base.family, fontface = "italic", size = 5) +
   scale_colour_manual(values = col.txt, guide = "none") +
   coord_cartesian(ylim = c(0, 65)) +
   labs(y = "Réponse", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.lm.10.png")


(ggplot(dat) +
   geom_polygon(data = d.dat,
                aes(x = x.d, y = y.d), fill = "grey75", colour = NA) +
   geom_path(data = d.dat[-nrow(d.dat),],
             aes(x = x.d, y = y.d),
             colour = "grey35", linewidth = 0.3) +
   geom_vline(xintercept = -0.00375, linetype = "dashed", colour = "grey5") +
   geom_segment(x = -0.00375, xend = -0.00375, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = -0.00375, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = -0.00375, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = -0.00375, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   scale_x_continuous(breaks = c(0, 0.05)) +
   coord_cartesian(xlim = c(-0.01, 0.05), ylim = c(0, 65)) +
   theme(axis.text.x = element_text(color = "white")) +
   labs(x = "Densité", y = "Réponse")) |>
plot_d("../qmd/figures/res.d.5.png")

(ggplot(dat) +
   geom_polygon(data = d.foc.dat,
                aes(x = x.foc.d, y = y.foc.d),
                fill = "white", colour = NA) +
   geom_polygon(data = d.foc.dat,
                aes(x = x.foc.d, y = y.foc.d),
                fill = pal[3], colour = NA, alpha = 0.8) +
   geom_path(data = d.dat[-nrow(d.dat),],
             aes(x = x.d, y = y.d),
             colour = "grey35", linewidth = 0.3) +
   geom_vline(xintercept = -0.00375, linetype = "dashed", colour = "grey5") +
   geom_segment(x = -0.00375, xend = -0.00375, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = -0.00375, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = -0.00375, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = -0.00375, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_text(data = data.frame(x = x.foc.d[length(x.foc.d)-2],
                               y = y.foc-2,
                               lab = "résidu\nquantile"),
             aes(x = x, y = y, label = lab, vjust = vjust), colour = pal[3],
             family = base.family, fontface = "italic",
             hjust = 0, vjust = 1, size = 4) +
   scale_x_continuous(breaks = c(0, 0.05)) +
   coord_cartesian(xlim = c(-0.01, 0.05), ylim = c(0, 65)) +
   theme(axis.text.x = element_text(color = "white")) +
   labs(x = "Densité", y = "Réponse")) |>
plot_d("../qmd/figures/res.d.6.png")


# Gamma distribution, first focal point

set.seed(123)
n <- 25
a <- 0
b <- 0.04
s <- 10
x <- sample(1:100, n)
mu <- exp(a + b * x)
shape <- 4
scale <- mu/shape
y <- rgamma(n, scale = scale, shape = shape)
dat <- data.frame(x = x, y = y)
mod <- glm(y ~ x, family = Gamma(link = "log"), data = dat)
mod.c <- coef(mod)
mod.s <- sigma(mod)
pred.dat <- data.frame(x = 1:100)
pred.dat$yhat <- predict(mod, pred.dat, type = "response")
dat$mu <- fitted(mod)
dat$r <- residuals(mod, type = "working")
qres <- simulateResiduals(mod, n = 500)
rq <- residuals(qres)
x.foc <- 90
y.foc <- y[which(x == x.foc)]
mu.foc <- fitted(mod)[which(x == x.foc)]
r.foc <- residuals(mod)[which(x == x.foc)]
seg.dat <- data.frame(x = x, mu = fitted(mod), y = y)
y.d <- seq(mu.foc-100, mu.foc+100, length.out = 501)
y.d <- c(y.d, y.d[1])
x.d <- dgamma(y.d, shape = shape, scale = mu.foc/shape)
d.dat <- data.frame(x.d = x.d, y.d = y.d, r.d = y.d-mu.foc)
y.foc.d <- seq(mu.foc-100, y.foc, length.out = 501)
x.foc.d <- dgamma(y.foc.d, shape = shape, scale = mu.foc/shape)
y.foc.d <- c(y.foc.d, y.foc.d[length(y.foc.d)], 0)
x.foc.d <- c(x.foc.d, x.foc.d[1], x.foc.d[1])
d.foc.dat <- data.frame(x.foc.d = x.foc.d, y.foc.d = y.foc.d)
pal <- c("#e41a1c", "#377eb8", "#984ea3")
lab <-
  data.frame(x = c(x.foc+2, x.foc+2, x.foc-2),
             y = c(y.foc+2, mu.foc-2, (y.foc+mu.foc)/2),
             lab1 = c("y", "µ", "r"),
             lab2 = c("observée", "prédite", "résidu"),
             hjust = c(0, 0, 1),
             col = as.character(1:3))
col.txt <- c("black", pal[2], pal[1])
names(col.txt) <- lab$col
size.pt.small <- 2 
size.pt.large <- 3
stroke.foc <- 0.2
lw.res <- 0.75

(ggplot(dat) +
   geom_point(aes(x = x, y = y),
              size = size.pt.small, colour = "grey5") +
   geom_line(data = pred.dat, aes(x = x, y = yhat),
             colour = pal[2]) +
   coord_cartesian(ylim = c(0, 70)) +
   labs(y = "Réponse", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.glm.1.png")

(ggplot(dat) +
   geom_point(aes(x = x, y = y),
              size = size.pt.small, colour = "grey50") +
   geom_line(data = pred.dat, aes(x = x, y = yhat),
             colour = pal[2]) +
   geom_vline(xintercept = x.foc, linetype = "dashed", colour = "grey5") +
   geom_segment(x = x.foc, xend = x.foc, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = x.foc, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = x.foc, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = x.foc, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_text(data = lab, aes(x = x, y = y, label = lab1, hjust = hjust, colour = col),
             family = base.family, fontface = "italic", size = 5) +
   scale_colour_manual(values = col.txt, guide = "none") +
   coord_cartesian(ylim = c(0, 70)) +
   labs(y = "Réponse", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.glm.2.png")


(ggplot(dat) +
   geom_polygon(data = d.dat,
                aes(x = x.d, y = y.d), fill = "grey75", colour = NA) +
   geom_path(data = d.dat[-nrow(d.dat),],
             aes(x = x.d, y = y.d),
             colour = "grey35", linewidth = 0.3) +
   geom_vline(xintercept = -0.00375, linetype = "dashed", colour = "grey5") +
   geom_segment(x = -0.00375, xend = -0.00375, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = -0.00375, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = -0.00375, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = -0.00375, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   scale_x_continuous(breaks = c(0, 0.05)) +
   coord_cartesian(xlim = c(-0.01, 0.05), ylim = c(0, 70)) +
   theme(axis.text.x = element_text(color = "white")) +
   labs(x = "Densité", y = "Réponse")) |>
plot_d("../qmd/figures/res.d.7.png")


(ggplot(dat) +
   geom_polygon(data = d.dat,
                aes(x = x.d, y = y.d), fill = "grey75", colour = NA) +
   geom_polygon(data = d.foc.dat,
                aes(x = x.foc.d, y = y.foc.d),
                fill = "white", colour = NA) +
   geom_polygon(data = d.foc.dat,
                aes(x = x.foc.d, y = y.foc.d),
                fill = pal[3], colour = NA, alpha = 0.8) +
   geom_path(data = d.dat[-nrow(d.dat),],
             aes(x = x.d, y = y.d),
             colour = "grey35", linewidth = 0.3) +
   geom_vline(xintercept = -0.00375, linetype = "dashed", colour = "grey5") +
   geom_segment(x = -0.00375, xend = -0.00375, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = -0.00375, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = -0.00375, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = -0.00375, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   scale_x_continuous(breaks = c(0, 0.05)) +
   coord_cartesian(xlim = c(-0.01, 0.05), ylim = c(0, 70)) +
   theme(axis.text.x = element_text(color = "white")) +
   labs(x = "Densité", y = "Réponse")) |>
plot_d("../qmd/figures/res.d.8.png")


# Gamma distribution, second focal point

set.seed(123)
n <- 25
a <- 0
b <- 0.04
s <- 10
x <- sample(1:100, n)
mu <- exp(a + b * x)
shape <- 4
scale <- mu/shape
y <- rgamma(n, scale = scale, shape = shape)
dat <- data.frame(x = x, y = y)
mod <- glm(y ~ x, family = Gamma(link = "log"), data = dat)
mod.c <- coef(mod)
mod.s <- sigma(mod)
pred.dat <- data.frame(x = 1:100)
pred.dat$yhat <- predict(mod, pred.dat, type = "response")
dat$mu <- fitted(mod)
dat$r <- residuals(mod, type = "response")
qres <- simulateResiduals(mod, n = 500)
rq <- residuals(qres)
x.foc <- 69
y.foc <- y[which(x == x.foc)]
mu.foc <- fitted(mod)[which(x == x.foc)]
r.foc <- residuals(mod)[which(x == x.foc)]
seg.dat <- data.frame(x = x, mu = fitted(mod), y = y)
y.d <- seq(mu.foc-100, mu.foc+100, length.out = 501)
y.d <- c(y.d, y.d[1])
x.d <- dgamma(y.d, shape = shape, scale = mu.foc/shape)
d.dat <- data.frame(x.d = x.d, y.d = y.d, r.d = y.d-mu.foc)
y.foc.d <- seq(mu.foc-100, y.foc, length.out = 501)
x.foc.d <- dgamma(y.foc.d, shape = shape, scale = mu.foc/shape)
y.foc.d <- c(y.foc.d, y.foc.d[length(y.foc.d)], 0)
x.foc.d <- c(x.foc.d, x.foc.d[1], x.foc.d[1])
d.foc.dat <- data.frame(x.foc.d = x.foc.d, y.foc.d = y.foc.d)
pal <- c("#e41a1c", "#377eb8", "#984ea3")
lab <-
  data.frame(x = c(x.foc+2, x.foc+2, x.foc-2),
             y = c(y.foc+2, mu.foc-2, (y.foc+mu.foc)/2),
             lab1 = c("y", "µ", "r"),
             lab2 = c("observée", "prédite", "résidu"),
             hjust = c(0, 0, 1),
             col = as.character(1:3))
col.txt <- c("black", pal[2], pal[1])
names(col.txt) <- lab$col
size.pt.small <- 2 
size.pt.large <- 3
stroke.foc <- 0.2
lw.res <- 0.75

(ggplot(dat) +
   geom_point(aes(x = x, y = y),
              size = size.pt.small, colour = "grey50") +
   geom_line(data = pred.dat, aes(x = x, y = yhat),
             colour = pal[2]) +
   geom_vline(xintercept = x.foc, linetype = "dashed", colour = "grey5") +
   geom_segment(x = x.foc, xend = x.foc, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = x.foc, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = x.foc, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = x.foc, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_text(data = lab, aes(x = x, y = y, label = lab1, hjust = hjust, colour = col),
             family = base.family, fontface = "italic", size = 5) +
   scale_colour_manual(values = col.txt, guide = "none") +
   coord_cartesian(ylim = c(0, 70)) +
   labs(y = "Réponse", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.glm.3.png")

(ggplot(dat) +
   geom_polygon(data = d.dat,
                aes(x = x.d, y = y.d), fill = "grey75", colour = NA) +
   geom_path(data = d.dat[-nrow(d.dat),],
             aes(x = x.d, y = y.d),
             colour = "grey35", linewidth = 0.3) +
   geom_vline(xintercept = -0.005, linetype = "dashed", colour = "grey5") +
   geom_segment(x = -0.005, xend = -0.005, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = -0.005, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = -0.005, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = -0.005, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   scale_x_continuous(breaks = c(0, 0.05)) +
   coord_cartesian(xlim = c(-0.01, 0.07), ylim = c(0, 70)) +
   theme(axis.text.x = element_text(color = "white")) +
   labs(x = "Densité", y = "Réponse")) |>
plot_d("../qmd/figures/res.d.9.png")


(ggplot(dat) +
   geom_polygon(data = d.dat,
                aes(x = x.d, y = y.d), fill = "grey75", colour = NA) +
   geom_polygon(data = d.foc.dat,
                aes(x = x.foc.d, y = y.foc.d),
                fill = "white", colour = NA) +
   geom_polygon(data = d.foc.dat,
                aes(x = x.foc.d, y = y.foc.d),
                fill = pal[3], colour = NA, alpha = 0.8) +
   geom_path(data = d.dat[-nrow(d.dat),],
             aes(x = x.d, y = y.d),
             colour = "grey35", linewidth = 0.3) +
   geom_vline(xintercept = -0.005, linetype = "dashed", colour = "grey5") +
   geom_segment(x = -0.005, xend = -0.005, y = mu.foc, yend = y.foc,
                colour = pal[1], linewidth = lw.res) +
   geom_point(x = -0.005, y = y.foc, fill = 1,
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   geom_point(x = -0.005, y = y.foc, colour = "white",
              shape = 16, size = size.pt.large/2) +
   geom_point(x = -0.005, y = mu.foc, fill = pal[2],
              shape = 21, size = size.pt.large, stroke = stroke.foc) +
   scale_x_continuous(breaks = c(0, 0.05)) +
   coord_cartesian(xlim = c(-0.01, 0.07), ylim = c(0, 70)) +
   theme(axis.text.x = element_text(color = "white")) +
   labs(x = "Densité", y = "Réponse")) |>
plot_d("../qmd/figures/res.d.10.png")



(ggplot(dat) +
  geom_point(aes(x = x, y = r),
             size = size.pt.small, colour = "grey5") +
  labs(y = "Résidu", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.glm.4.png")

(ggplot(dat) +
  geom_point(aes(x = x, y = rq),
             size = size.pt.small, colour = pal[3]) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Résidu quantile", x = "Variable explicative")) |>
plot_res("../qmd/figures/res.glm.5.png")


plot(qres)


# DHARMa plots

set.seed(1)
dat <- createData(family = poisson())
mod <- glmmTMB(observedResponse ~ x + Environment1 + (1|group) , data = dat, family = poisson)
mod.qres <- simulateResiduals(mod)

png("../qmd/figures/dharma.1.png", height = 5, width = 8, units = "in", res = 300)
plot(mod.qres)
dev.off()

png("../qmd/figures/dharma.2.png", height = 5, width = 5, units = "in", res = 300)
testUniformity(mod.qres)
dev.off()

png("../qmd/figures/dharma.3.png", height = 5, width = 5, units = "in", res = 300)
plotResiduals(mod.qres)
dev.off()


set.seed(1)
dat <- createData(sampleSize = 100, family = poisson(), overdispersion = 1.5)
mod <- glmmTMB(observedResponse ~ Environment1 , data = dat, family = poisson)
mod.qres <- simulateResiduals(mod)
png("../qmd/figures/dharma.4.png", height = 5, width = 5, units = "in", res = 300)
testUniformity(mod.qres)
dev.off()

set.seed(1)
dat <- createData(sampleSize = 100, intercept = 0,
                  fixedEffects = 2, overdispersion = 0,
                  family = poisson(),
                  roundPoissonVariance = 0.001, randomEffectVariance = 0)
mod <- glmmTMB(observedResponse ~ Environment1 , data = dat, family = poisson)
mod.qres <- simulateResiduals(mod)
testUniformity(mod.qres)
png("../qmd/figures/dharma.5.png", height = 5, width = 5, units = "in", res = 300)
testUniformity(mod.qres)
dev.off()

set.seed(3)
dat <- createData(sampleSize = 500, intercept = 0,
                  fixedEffects = 4, overdispersion = 0,
                  family = poisson(),
                  roundPoissonVariance = 0.05, randomEffectVariance = 0)
mod <- glmmTMB(observedResponse ~ x + Environment1 + (1|group), data = dat, family = poisson)
mod.qres <- simulateResiduals(mod)
png("../qmd/figures/dharma.6.png", height = 5, width = 5, units = "in", res = 300)
plotResiduals(mod.qres)
dev.off()

