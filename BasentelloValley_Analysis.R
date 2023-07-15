################################################################################
#
# Multi-isotope analysis
# Basentello Valley estate:
# Vagnari vicus & San Felice Villa
#
# Angela Trentacoste
# May 2023
#
# Analysis script for the paper:
# Isotopic insights into livestock production in Roman Italy: 
# diet, seasonality, and mobility on an imperial estate
# Angela Trentacoste, Michael Mackinnon, Christopher Day, Petrus Le Roux, 
# Mike Buckley, Myles McCallum, and Maureen Carroll
#
# If you use or adapt this script please cite the paper.
#
################################################################################
# Remember to set working directory
# load packages
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(ggrepel)
library(rstatix)
library(car)
library(rKIN)
library(ggalt)
################################################################################
# >>> BONE COLLAGEN -----------------------------------------------------------
################################################################################
# assign input file
input.collagen <- c("BVestate_collagen_SI.csv")
# Load isotope data
Coldata <- read.csv(input.collagen)

##### Quality indicators -------------------------------------------------------
# . Collagen yield -----
ggplot(Coldata, aes(x="Collagen yield", y=CollagenYield, 
                    label = Sample, color = Site)) +
  geom_jitter() + 
  geom_hline(yintercept = 0.005, color = "red") + 
  geom_hline(yintercept = 0.01, color = "red", 
             linetype = "dashed")
# . C:N -----
# Guiry and Szpak 2021 = 3.15–3.45/5 (conservative), 3.9 (liberal)
# van Klinken 1999: C:N = 3.1-3.65
ggplot(Coldata, aes(x="C/N", y=C.N, label = Sample, color = Site)) +
  geom_jitter() + 
  geom_hline(yintercept = 2.9, color = "red") + 
  geom_hline(yintercept = 3.6, color = "red") +
  geom_hline(yintercept = 3.15, color = "red", linetype="dashed") + 
  geom_hline(yintercept = 3.45, color = "red", linetype="dashed") 
# . Percent C and N -----
ggplot(Coldata, aes(perC, perN, color = Site, label = Sample)) +
  geom_point() + 
  scale_x_continuous(breaks = seq(from = 15, to = 50, by = 1)) +
  scale_y_continuous(breaks = seq(from = 5, to = 20, by = 1)) +
  geom_vline(xintercept = 30) + geom_vline(xintercept = 46) + 
  geom_hline(yintercept = 10) + geom_hline(yintercept = 17) +
  geom_hline(yintercept = 16, 
             linetype = 2)
# Zoom in
ggplot(Coldata, aes(perC, perN, color = Site, label = Sample)) +
  geom_point() + 
  scale_x_continuous(breaks = seq(from = 15, to = 50, by = 1)) +
  scale_y_continuous(breaks = seq(from = 5, to = 20, by = 1)) +
  geom_vline(xintercept = 30) + geom_vline(xintercept = 46) + 
  geom_hline(yintercept = 10) + geom_hline(yintercept = 17) +
  geom_hline(yintercept = 16, 
             linetype = 2) + 
  xlim(43, 50) + 
  ylim(16, 18) +
  geom_text_repel(size = 3)
# Sample VV10 was re-run (run 08) and produced normal %C and %N 
# and comparable stable isotope values (see below) and was therefore included.

##### Replicates and duplicates ------------------------------------------------
# Create df for replicates
Reps <- Coldata %>% 
  group_by(Sample) %>% 
  filter(n()>1) 
# Calculate maximum range in d15N and d13C values for replicates
Reps.range <- Reps %>% 
  group_by(Site, Sample) %>%
  summarise(
    C.range = max(d13C)-min(d13C),
    N.range = max(d15N)-min(d15N))
# check range of values
ggplot(Reps.range, aes(C.range, y = N.range, 
                 color = Site,
                 label=Sample)) +
  geom_point() +
  geom_text_repel(size = 3)
# check individual values
ggplot(Reps, aes(d13C, d15N, color = Site)) +
  geom_point() + 
  facet_wrap(vars(Sample))
# Remove third reps that are unacceptable
Coldata.QC <- Coldata %>%
  subset(!(Sample == "SF56" & d15N < 5))
# Remove samples that are unacceptable
exclude <- c("SF65") # d13C difference >0.3, not enough collagen to repeat again
Coldata.QC <- subset(Coldata.QC, !(Sample %in% exclude))
# Average reps
Coldata.QC <- Coldata.QC %>%
  group_by(Sample, Site, Taxon, TaxGroup, TaxGroupComp, MajorGroup) %>% 
  summarise(d15N.m = mean(d15N),
            d13C.m = mean(d13C),
            C.N.m = mean(C.N)) 
Coldata.QC <- Coldata.QC %>% 
  rename(d15N = d15N.m,
         d13C = d13C.m,
         C.N = C.N.m)

#####  C:N vs SI values --------------------------------------------------------
# d13C vs C:N
ggscatter(Coldata.QC, 
          x = "C.N", y = "d13C", 
          cor.coeff.args = list(method = "pearson", color = "red"),
          add = "reg.line",
          cor.coef = TRUE) +
  facet_wrap(Site~.)
# d15N vs C:N
ggscatter(Coldata.QC, 
          x = "C.N", y = "d15N", 
          cor.coeff.args = list(method = "pearson", color = "red"),
          add = "reg.line",
          cor.coef = TRUE) +
  facet_wrap(Site~.)

# write dataset
S2 <- Coldata.QC %>% 
  select(Sample:d13C) %>%
  mutate(d15N = round(d15N, 2),
         d13C = round(d13C, 2))
write.csv(S2, "Supplement02_Table_CollagenResults.csv")

##### Summary statistics -------------------------------------------------------
# Exclude sample VV31 - unidentified rib, likely caprine
Coldata.QC <- subset(Coldata.QC, !(Sample %in% c("VV31")))

# Set taxon order --------------------------------------------------------------
#Set preferred order of taxa and taxa groups
TaxOrder <- c("Cattle", "Sheep/goat", "Sheep", "Goat", 
              "Pig", "Equid", "Red deer", "Dog")
TaxGroupOrder <- c("Bos", "Ovis/Capra", "Sus", "Equus", "Cervus elaphus", 
                   "Cervus elaphus?", "Canis")
compgroups <- c("Cattle", "Sheep/goat", "Pig", "Equid", "Deer", "Canid")
Coldata.QC$Taxon = factor(Coldata.QC$Taxon, levels = TaxOrder)
Coldata.QC$TaxGroup = factor(Coldata.QC$TaxGroup, levels = TaxGroupOrder)
Coldata.QC$TaxGroupComp = factor(Coldata.QC$TaxGroupComp, levels = compgroups)

# . By site -----
SiteSumStats <- Coldata.QC %>%
  select(Site, TaxGroupComp, d15N, d13C) %>%
  pivot_longer(d15N:d13C, 
               names_to = "Isotope", values_to = "Value") %>%
  group_by(Site, TaxGroupComp, Isotope) %>%
  summarise(
    n = n(),
    min = min(Value),
    max = max(Value),
    range = max-min,
    mean = mean(Value),
    median = median(Value),
    sd = sd(Value)) 
relocate_order <- c("Site", "TaxGroupCom", "n","d15N", "d13C")
SiteSumStats <- SiteSumStats %>%
  pivot_wider(names_from = Isotope,
               names_glue = "{Isotope}_{.value}",
               values_from = c(min, max, range, mean, median, sd)) %>% 
  relocate(starts_with(relocate_order)) %>% 
  mutate_at(4:15, round, 2)
# Table S1.5a Summary statistics for bone collagen isotopic values - by site
write.csv(SiteSumStats, "Table_S1_5a_Collagen_SumStats_bySite.csv")
# . All grouped -----
GroupedSumStats <- Coldata.QC %>%
  select(TaxGroupComp, d15N, d13C) %>%
  pivot_longer(d15N:d13C, 
               names_to = "Isotope", values_to = "Value") %>%
  group_by(TaxGroupComp, Isotope) %>%
  summarise(
    n = n(),
    min = min(Value),
    max = max(Value),
    range = max-min,
    mean = mean(Value),
    median = median(Value),
    sd = sd(Value)) 
GroupedSumStats  <- GroupedSumStats  %>%
  pivot_wider(names_from = Isotope,
              names_glue = "{Isotope}_{.value}",
              values_from = c(min, max, range, mean, median, sd)) %>% 
  relocate(starts_with(relocate_order)) %>% 
  mutate_at(3:14, round, 2) 
# Table S1.5b Summary statistics for bone collagen isotopic values - grouped
write.csv(GroupedSumStats, "Table_S1_5b_Collagen_SumStats_Grouped.csv")

##### Comparative statistics -------------------------------------------------------
# Table S1.6 Results of multi-variate statistical tests 
# on isotope values from faunal collagen 
# Select main taxa and isotope columns
maindom <- c("Bos", "Ovis/Capra","Sus")
# . All grouped -----
Stat.data <- Coldata.QC %>% 
  subset(TaxGroup %in% maindom) %>%
  ungroup() %>%
  select(TaxGroup, d13C, d15N)
# Pivot data
Coldatalong <- Stat.data %>%  
  pivot_longer(d13C:d15N, names_to = "Isotope", values_to = "SIvalue")
# Seperate N and C SI data
Ndata <- subset(Coldatalong, Isotope == "d15N")
Cdata <- subset(Coldatalong, Isotope == "d13C")
# .. Extreme outliers
ggboxplot(Cdata, x = "TaxGroup", y = "SIvalue")
ggboxplot(Ndata, x = "TaxGroup", y = "SIvalue")
Coldatalong %>% 
  group_by(TaxGroup, Isotope) %>%
  identify_outliers(SIvalue)
# .. Shapiro-Wilk’s test
Coldatalong %>%
  group_by(TaxGroup, Isotope) %>%
  shapiro_test(SIvalue)
# .. Levene test
leveneTest(SIvalue ~ TaxGroup, data = Cdata)
leveneTest(SIvalue ~ TaxGroup, data = Ndata)
# ... ANOVA 
# We want to know if there is any significant difference between the average weights of plants in the 3 experimental conditions.
Cdata %>% anova_test(SIvalue ~ TaxGroup)
# ... Kruskal Wallis
Ndata %>% kruskal_test(SIvalue ~ TaxGroup)
Ndata %>% kruskal_effsize(SIvalue ~ TaxGroup)
# . By site -----
Stat.data <- Coldata.QC %>% 
  subset(TaxGroup %in% maindom) %>%
  ungroup() %>%
  select(Site, TaxGroup, d13C, d15N)
# Pivot data
Coldatalong <- Stat.data %>%  
  pivot_longer(d13C:d15N, names_to = "Isotope", values_to = "SIvalue")
# Separate N and C SI data
Ndata <- subset(Coldatalong, Isotope == "d15N")
Cdata <- subset(Coldatalong, Isotope == "d13C")
# .. Extreme outliers
Coldatalong %>% 
  group_by(Site, TaxGroup, Isotope) %>%
  identify_outliers(SIvalue)
# .. Shapiro-Wilk’s test
Coldatalong %>%
  group_by(Site, TaxGroup, Isotope) %>%
  shapiro_test(SIvalue)
## .. Levene test
Bos <- subset(Ndata, TaxGroup == "Bos")
Ovis.Capra <- subset(Ndata, TaxGroup == "Ovis/Capra")
Sus <- subset(Ndata, TaxGroup == "Sus")
leveneTest(SIvalue ~ Site, data = Bos)
leveneTest(SIvalue ~ Site, data = Ovis.Capra)
leveneTest(SIvalue ~ Site, data = Sus)
Bos <- subset(Cdata, TaxGroup == "Bos")
Ovis.Capra <- subset(Cdata, TaxGroup == "Ovis/Capra")
Sus <- subset(Cdata, TaxGroup == "Sus")
leveneTest(SIvalue ~ Site, data = Bos)
leveneTest(SIvalue ~ Site, data = Ovis.Capra)
leveneTest(SIvalue ~ Site, data = Sus)
# .. t-tests
Coldatalong %>%
  group_by(Isotope, TaxGroup) %>%
  t_test(SIvalue ~ Site, var.equal = TRUE) %>%
  add_significance()
# .. Wilcoxon rank-sum test / Mann-Whitney U-test
# ... compare site
Coldatalong %>% 
  group_by(Isotope, TaxGroup) %>%
  wilcox_test(SIvalue ~ Site) %>%
  add_significance()
# ... compare axa
Coldatalong %>% 
  group_by(Isotope) %>%
  wilcox_test(SIvalue ~ TaxGroup) %>%
  add_significance()
# effect size
Coldatalong %>% 
  group_by(Isotope) %>%
  wilcox_effsize(SIvalue ~ TaxGroup)

##### Plots --------------------------------------------------------------------
# Load comparative data
# Human data from Semchuk, L. (2016). A stable isotope investigation of diet 
# at Vagnari. MA thesis, School of Graduate Studies, McMaster University.
# http://hdl.handle.net/11375/20498
input.comp.human <- c("Comparative_collagenSI_VagnariHuman.csv")
Human <- read.csv(input.comp.human)
#Define groups
maindom <- c("Bos", "Ovis/Capra","Sus")

##### . Figure 4 ----------------------------------------------------------------
# δ13C and δ15N values from animal bone collagen 
# compared to values from humans from the Vagnari cemetery (gray crosses)
mypal <- c("palevioletred4", "skyblue3", "turquoise4", "goldenrod",  
           "springgreen2", "coral2")
# calculate summary statistics for major groups
SumstatsPlotMajor <- Coldata.QC %>%
  select(MajorGroup, d15N, d13C) %>%
  group_by(MajorGroup) %>%
  summarise(
    n = n(),
    Nmin = min(d15N),
    Nmax = max(d15N),
    Nrange = Nmax-Nmin,
    Nmean = mean(d15N),
    Nmedian = median(d15N),
    Nsd = sd(d15N),
    Cmin = min(d13C),
    Cmax = max(d13C),
    Crange = Cmax-Cmin,
    Cmean = mean(d13C),
    Cmedian = median(d13C),
    Csd = sd(d13C))
HerbComp <- SumstatsPlotMajor[1,] 
HerbComp[1,1] <- "Omnivores"
# text annotation
Plottext <- data.frame(
  MajorGroup = c("Herbivores", "Omnivores"),
  text = c("", "herbivores"))
# draw plot
ggplot() + 
    geom_errorbarh(data = HerbComp, 
                   aes(y = Nmean, xmax = Cmax, xmin = Cmin, 
                       height = 0.15), 
                   linetype = "longdash", colour = "skyblue3", size = 0.3) +
  geom_errorbar(data = HerbComp,
                aes(x = Cmean, ymax = Nmax, ymin = Nmin,
                    width = 0.15), 
                linetype = "longdash", colour = "skyblue3", size = 0.3) +
  geom_point(data = Coldata.QC, 
             aes(x=d13C, y=d15N, 
                 shape=Site, fill = TaxGroup, colour = TaxGroup), 
             size=1.5) +
  geom_point(data=Human, aes(x=d13C, y=d15N), size=1.5, 
             shape=3, colour="gray50", alpha=0.55) +
  scale_colour_manual(values = mypal) +
  scale_shape_manual(values = c(16,1)) +
  facet_wrap(.~MajorGroup) +
  geom_text(data=Plottext, aes(x=-20.3, y=13.3,label=text), colour ="skyblue3", 
            size = 2) +
  geom_segment(data = data.frame(x = -21.6, y = 13.28, xend = -20.9, 
                                 yend = 13.28, MajorGroup="Omnivores"), 
               aes(x = x, y = y, xend = xend, yend = yend), linetype = "dashed", 
               size = 0.3, colour = "skyblue3") +
  theme_bw(base_size = 7) +
  theme(
    strip.background = element_rect(fill = "white", linetype = "blank")) +
  labs(
    y=expression(paste(delta^{15}, "N (\u2030)")), 
    x = expression(paste(delta^{13}, "C (\u2030)")),
    colour = "Taxon", fill = "Taxon")
ggsave("Fig04_BVestate_Collagen.pdf", height = 3.6, width = 5.5, units = "in")

#### . Figure 5 ------------------------------------------------------------------
# Distribution of isotope values from animal bone collagen. 
# significance stars in published version added manually
Nhist <- ggplot(data = subset(Coldata.QC, TaxGroup %in% maindom), 
                aes(x=d15N, colour=Site, fill=Site, alpha=Site)) + 
  geom_histogram(bins = 16, position="stack") + 
  facet_grid(TaxGroup~.) +
  scale_fill_manual(values = c("turquoise4", "turquoise4")) +
  scale_colour_manual(values = c("turquoise4", "turquoise4")) +
  scale_alpha_manual(values=c(0.5, 0.8)) +
  xlab(expression(paste(delta^{15}, "N (\u2030)"))) +
  ylim(0,8) +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    legend.position="top",
    panel.border = element_rect(colour = "black", fill = NA))
Chist <- ggplot(data = subset(Coldata.QC, TaxGroup %in% maindom), 
                aes(x=d13C, colour=Site, fill=Site, alpha=Site)) + 
  geom_histogram(bins = 6, position="stack") + 
  facet_grid(TaxGroup~.) +
  scale_fill_manual(values = c("turquoise4", "turquoise4")) +
  scale_colour_manual(values = c("turquoise4", "turquoise4")) +
  scale_alpha_manual(values=c(0.5, 0.8)) +
  ylim(0,8) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    legend.position="top",
    panel.border = element_rect(colour = "black", fill = NA))
ggsave(Nhist, file="Fig05_N_BVestate_Histogram.pdf", height = 3.5, width = 2, units = "in")
ggsave(Chist, file="Fig05_C_BVestate_Histogram.pdf", height = 3.5, width = 1.08, units = "in")
# N and C histograms then combined in Adobe Illustrator to create final figure.

##### . Figure 7  ------------------------------------------------------------------
# Isotope values from animal bone collagen compared to other Roman and Etruscan sites. 
input.comp.fauna <- c("Comparative_collagenSI_ItalyFauna.csv")
CompCol <- read.csv(input.comp.fauna)
compgroups <- c("Cattle", "Sheep/goat", "Pig", "Equid", "Deer", "Canid")
area.order <- c("Basentello estate",
                "Pompeii/Herculaneum",
                "Velia",
                "Portus",
                "Isola Sacra",
                "Etruscan C Italy")
CompPlot1 <- subset(CompCol, TaxGroupComp %in% compgroups)
CompPlotV <- subset(Coldata.QC, TaxGroupComp %in% compgroups)
CompPlot1$TaxGroupComp = factor(CompPlot1$TaxGroupComp, levels = compgroups)
CompPlotV$TaxGroupComp = factor(CompPlotV$TaxGroupComp, levels = compgroups)
CompPlot1$Area = factor(CompPlot1$Area, levels = area.order)
CompPlotV$Area <- "Basentello estate"
CompPlotV$Area = factor(CompPlotV$Area, levels = area.order)
ggplot() + 
  geom_point(data= CompPlotV, aes(x=d13C, y=d15N, shape=Area, colour = Area), 
             size=1.5, alpha = 0.8) +
  geom_point(data=CompPlot1, aes(x=d13C, y=d15N, shape=Area, colour = Area), 
             size=2, alpha=0.7) +
  facet_wrap(.~TaxGroupComp) +
  scale_colour_manual(values = c("turquoise4", "palevioletred4", "skyblue3", 
                                 "goldenrod",  "springgreen2", "coral2", 
                                 "olivedrab4"), 
                      limits = area.order) +
  scale_shape_manual(values = c(16, 1, 2, 0, 5, 6), limits = area.order) +
  scale_y_continuous(breaks = seq(2, 12, by = 2 )) +
  theme_bw(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "white", linetype = "blank"),
    panel.grid.minor.x = element_blank() 
        ) +
  labs(
    y=expression(paste(delta^{15}, "N (\u2030)")), 
    x = expression(paste(delta^{13}, "C (\u2030)")))
ggsave("Fig07_Comparative_Collagen_SIvalues.pdf", width=7, height = 5.5, units = "in")

##### rKIN  ------------------------------------------------------------------
comp.species <- c("Bos", "Ovis/Capra", "Cervus elaphus", "Sus")
ColdataMain <- subset(Coldata.QC, TaxGroup %in% comp.species)
df <- data.frame (C  = ColdataMain$d13C,
                  N = ColdataMain$d15N,
                  Group = ColdataMain$TaxGroup,
                  Site = ColdataMain$Site)
kin.50 <- estKIN(data=df, x="C", y="N", group="Group",
                 levels=50, smallSamp = TRUE)
kin.75 <- estKIN(data=df, x="C", y="N", group="Group",
                 levels=75, smallSamp = TRUE)
kin.95 <- estKIN(data=df, x="C", y="N", group="Group",
                 levels=95, smallSamp = TRUE)
# extract the area of each polygon at each contour interval
kin.50.area = getArea(kin.50)
kin.75.area = getArea(kin.75)
kin.95.area = getArea(kin.95)
# Table S1.7 Isotopic niche size estimated based on 
# kernel utilization density (KUD) 
write.csv(kin.50.area, "BVestate_50_area.csv")
write.csv(kin.75.area , "BVestate_75_area.csv")
write.csv(kin.95.area, "BVestate_95_area.csv")
# determine polygon overlap for all polygons at each contour interval
kin.50.olp = calcOverlap(kin.50)
kin.75.olp = calcOverlap(kin.75)
kin.95.olp = calcOverlap(kin.95)
# Table S1.8 Isotopic niche overlap estimated based on 
# kernel utilization density (KUD) 
write.csv(kin.50.olp, "BVestate_50_olp.csv")
write.csv(kin.75.olp, "BVestate_75_olp.csv")
write.csv(kin.95.olp, "BVestate_95_olp.csv")

##### . Figure 6 ---------------------------------------------------------------
# Niche size and overlap estimated using kernel utilization density (KUD) 
# at different contour intervals (50%, 75%, and 95%). 
# plots saved separately and manually combined
plot.50 = plotKIN(kin.50, scaler = 1, title = "Kernel Estimates",
                  xlab = expression({delta}^13*C~ ('\u2030')),
                  ylab = expression({delta}^15*N~ ('\u2030')))
plot.75 = plotKIN(kin.75, scaler = 1, title = "Kernel Estimates",
                  xlab = expression({delta}^13*C~ ('\u2030')),
                  ylab = expression({delta}^15*N~ ('\u2030')))
plot.95 = plotKIN(kin.95, scaler = 1, title = "Kernel Estimates",
                  xlab = expression({delta}^13*C~ ('\u2030')),
                  ylab = expression({delta}^15*N~ ('\u2030')))

################################################################################
# >>> TOOTH ENAMEL -------------------------------------------------------------
################################################################################
##### FTIR: Enamel indices  ----------------------------------------------------
input.ftir <- c("BVestate_FTIR.csv")
FTIR <- read.csv(input.ftir)
FTIR.indices <- FTIR %>%
  select(Tooth, Peak, Abs_intensity) %>%
  pivot_wider(names_from = Peak, names_glue = "Peak_{Peak}", 
              values_from = Abs_intensity) %>%
  group_by(Tooth) %>%
  summarise(C.P = (Peak_1415/Peak_1035),
            IRSF = (Peak_565 + Peak_605) / Peak_590,
            C.C = Peak_1455 / Peak_1415,
            BPI = Peak_1415/Peak_605,
            API = Peak_1540/Peak_605) %>%
  mutate(across(C.P:API, round, 2))
# Table S1.4 ATR-FTIR enamel indices (peak height ratios) 
# following France, Sugiyama et al. (2020)
write.csv(FTIR.indices,"Table_S1_4_FTIR_EnamelIndices.csv")

##### Summary stats and plots --------------------------------------------------
# assign input files
input.enamel <- c("BVestate_enamel_SI.csv")
input.enamel.comp <- c("Comparative_enamel_VagnariCemetery.csv")
Results <- read.csv(input.enamel)
Human <- read.csv(input.enamel.comp)

# select San Felice teeth and set order
SanFelice <- subset(Results, Site == "SF")
curveorder <- c("SF_04","SF_03",  "SF_05", "SF_07", "SF_01", "SF_06", "SF_08")
SanFelice$Order <- factor(SanFelice$Tooth, curveorder)  
Results$Tooth <- as.character(Results$Tooth)
toothorder <- c("1", "3", "4", "5", "6", "7", "8", 
                "9", "10", "11", "12", "13", "16", "18", "19")
Results$Order <- factor(Results$Tooth, toothorder)  
ResultsLong <- Results %>%
  pivot_longer(cols = c("C","O"), names_to = "Type", values_to = "Value")
# summary statistics
SumStats <- Results %>%
  group_by(Site, Tooth, Taxon) %>%
  summarise(
    n = n(),
    Cmin = min(C),
    Cmax = max(C),
    Cmean = mean(C),
    Omin = min(O),
    Omax = max(O),
    Omean = mean(O)
  )
SumStats <- SumStats %>%
  mutate(Crange = Cmax - Cmin,
         Orange = Omax - Omin)
# Table S1.9 Summary statistics for intra-tooth carbon and 
# oxygen isotopic values from San Felice
write.csv(subset(SumStats, Site == "SF"), "Tab_S1_9_SummaryStats_SF_enamel.csv")

##### . Figure 8 -----
# δ18Ocarb  and δ13Ccarb  values from San Felice and  Vagnari vicus 
# compared to values from Vagnari humans and weather stations
input.GNIP <- c("GNIPdata.csv")
GNIP <- read.csv(input.GNIP, header = TRUE)
GNIP.sum.loc <- GNIP %>%
  group_by(location) %>%
  summarise(n = n(),
            GNIPmin = min(d18O.avg, na.rm = TRUE),
            GNIPmax = max(d18O.avg, na.rm = TRUE),
            GNIPmean = mean(d18O.avg, na.rm = TRUE))
colnames(GNIP.sum.loc) <- c("location", "n", "Omin",  "Omax",  "Omean", "Order")
GNIP.sum.loc$Order <- GNIP.sum.loc$location
SumStats$Order <- paste(SumStats$Taxon, SumStats$Tooth,  sep = "-")
SumStats <- SumStats %>%
  bind_rows(GNIP.sum.loc)
SumStats[14:16,1] <- "WS"
O <- ggplot(data=SumStats) +
  geom_pointrange(aes(x=Order, y=Omean, ymin = Omin, ymax = Omax, shape = Site, 
                      fill = Taxon, colour = Taxon), size = 0.4, stroke = 0.55) +
  geom_jitter(data = Human, aes(y=d18O.VPDB, x=".Human"), alpha = 0.24, 
              shape = 16, size=0.75, width=0.3) +
  geom_boxplot(data = Human, aes(y=d18O.VPDB, x=".Human"), notch = TRUE, 
               alpha = 0, lwd=0.35, width = 0.85, colour = "gray20") +
  scale_shape_manual(values = c(4, 0, 16)) +
  scale_fill_manual(values = c("cadetblue4", "sienna3")) +
  scale_colour_manual(values = c("cadetblue4", "sienna3")) +
  annotate("text", x = "..Inland", y = -1, 
           label = "Weather\n stations",
           colour = "grey",
           hjust = "middle",
           size = 2) +
  theme_bw(base_size = 7) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank()) +
  labs(y=expression(paste(delta^{18}, "O (\u2030)"))) + 
  scale_y_continuous(limits = c(-12.5,1.5), breaks = seq(from=-12, to=2 ,by=2))
O
C <- ggplot(data=subset(SumStats, !(Order %in%  c("..Inland", 
                                                  "..Coastal", "..Piloni")))) +
  geom_pointrange(aes(x=Order, y=Cmean,ymin = Cmin, ymax = Cmax, 
                      shape = Site, fill = Taxon, colour = Taxon), 
                  size = 0.4, stroke = 0.55) +
  scale_shape_manual(values = c(4, 0, 16)) +
  scale_fill_manual(values = c("cadetblue4", "sienna3")) +
  scale_colour_manual(values = c("cadetblue4", "sienna3")) +
  theme_bw(base_size = 7) +
  theme(legend.position="right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank()) +
  labs(y=expression(paste(delta^{13}, "C (\u2030)"))) +
  scale_y_continuous(limits = c(-19.5, -5.5), 
                     breaks = seq(from = -19, to=-5, by=2))
C
grid.arrange(O, C, nrow = 1, ncol = 2)
g <- arrangeGrob(O, C, nrow = 1, ncol = 2)
ggsave(g, file="Fig08_O_and_C_Enamelranges.pdf", 
       height = 3.5, width = 5.35, units = "in")

# O and C isotope curves -------------------------------------------------------
# C and O together
carb.plot <- ggplot(data=SanFelice, aes(x=ERJ))+
  geom_point(aes(y=O), colour = "cadetblue4", size = 0.7) +
  geom_xspline(aes(y=O),colour = "cadetblue4", size = .32) +
  geom_point(aes(y=C+6), shape=4, colour = "sienna4", 
             size = 0.75, stroke = 0.4) +
  geom_xspline(aes(y=C+6), colour = "sienna4", size = 0.32) +
  scale_x_reverse() +
  facet_wrap(Order ~ ToothTax, nrow = 2) +
  labs(x="distance ERJ (mm)") +
  scale_y_continuous(name=expression(paste(delta^{18}, "O (\u2030)")), 
                     sec.axis = sec_axis(~ . - 6, name=expression(paste(delta^{13}, "C (\u2030)")),
                                         breaks = seq(-18, 6, 2))) +
  theme_bw(base_size=7) +
  theme(strip.background =element_rect(fill="white"),
        axis.ticks.y.right = element_line(color = "sienna4"),
        axis.text.y.right = element_text(color = "sienna4"),
        axis.title.y.right = element_text(colour = "sienna4"),
        axis.ticks.y.left = element_line(color = "cadetblue4"),
        axis.text.y.left = element_text(color = "cadetblue4"),
        axis.title.y.left = element_text(colour = "cadetblue4"))
carb.plot
ggsave(carb.plot, file = "Fig09a_OandC_SI_curves.pdf", 
       height = 3.5, width = 5.35, units = "in")

### Seasonality ----------------------------------------------------------------
# Estimated using Chazin et al (2019)

# Chazin, H., Deb, S., Falk, J., and Srinivasan, A. (2019) 
# New Statistical Approaches to Intra-individual Isotopic Analysis 
# and Modelling of Birth Seasonality in Studies of Herd Animals*, 
# Archaeometry, 61, 478– 493. https://doi.org/10.1111/arcm.12432

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Requires script file in Chazin et al (2019) Data S5. Supporting information.
# Load functions and then implement the SCEM and Cosine method 
# on the San Felice data
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!

SanFelice <- read.csv(input.enamel, header=T)
# transform data to format needed by functions
SanFelice <- SanFelice %>% 
  subset(Site == "SF") %>%
  select(Tooth, O, ERJ, Taxon) %>%
  rename(ID = Tooth, distance = ERJ, oxygen = O)
# implement methods 
SF = split(SanFelice,f = SanFelice$ID)
results = SCEM(SF,bandwidth = -0.33)
cosine = makeFits(SF)  
rownames(cosine) = 1:nrow(cosine)
results$results$Cosine = cosine$birth
dif = abs(results$results$Cosine - results$results$Season)
results$results$Difference = apply(cbind(dif,(1-dif)),1,min)
results$results$R = cosine$Pearson
rrr.a = results$results[order(results$results$Cluster),]
# adjust maxillary values
Seasonality <- results$results
Seasonality <- Seasonality %>% rename(SCEM = Season)
Seasonality$adj.SCEM <- Seasonality$SCEM
Seasonality$adj.Cosine <- Seasonality$Cosine
max.teeth <- c("SF_03", "SF_04", "SF_05", "SF_07")
sheep <- c("SF_01", "SF_04", "SF_06", "SF_08")
Seasonality <- Seasonality %>%
  mutate(
    adj.SCEM = case_when(ID %in% max.teeth ~ adj.SCEM-0.073, 
                         TRUE ~ as.numeric(adj.SCEM)),
    adj.Cosine = case_when(ID %in% max.teeth ~ adj.Cosine-0.073, 
                           TRUE ~ as.numeric(adj.Cosine)),
    Element = case_when(ID %in% max.teeth ~ "Maxillary", TRUE ~ "Mandibular"),
    Taxon = case_when(ID %in% sheep ~ "Ovis", TRUE ~ "Capra"),
    Fit = case_when(R > 0.9 ~ "good", TRUE ~ "poor"))
# Table S1.10 Birth seasonality models  for caprine teeth from San Felice
write.csv(Seasonality, "Tab_S1_10_Seasonality_SF_enamel.csv")

# . Figure 9 -----
Seasonlong <- Seasonality %>%
  pivot_longer(adj.SCEM:adj.Cosine, names_to = "est.type", values_to = "sea.est")
toothorder <- c("SF_03", 
                "SF_04",
                "SF_07",
                "SF_08",
                "SF_05",
                "SF_01",
                "SF_06")
Seasonlong$ID <- factor(Seasonlong$ID, toothorder)
Seasonlong$Order <- factor(Seasonlong$ID, toothorder)
Seasonlong$est.type2 <- factor(Seasonlong$est.type, labels = c("cosine", "SCEM"))
ggplot(Seasonlong, aes(x=Order, y=sea.est)) +
  geom_vline(xintercept = c(0, 0.25, 0.5, 0.75), 
             colour = "gray70", size = 0.2, linetype="dashed") +
  geom_hline(yintercept = c(0, 1), colour = "black", size = 0.2) +
  geom_point(data=subset(Seasonlong, 
                         Element == "Maxillary" & est.type == "adj.Cosine"), 
             aes(x=sea.est, y=1, shape = est.type, colour = Fit, fill=Fit),
             size = 1.25) +
  geom_point(data=subset(Seasonlong, 
                         Element == "Mandibular" & est.type == "adj.Cosine"), 
             aes(x=sea.est, y=1, shape = est.type, colour = Fit), 
             size = 1.25, fill = "white") +
  geom_point(data=subset(Seasonlong, 
                         Element == "Maxillary" & est.type == "adj.SCEM"), 
             aes(x=sea.est, y=1, shape = est.type),
             size = 1.25,  colour = "sienna3", fill="sienna3") +
  geom_point(data=subset(Seasonlong, 
                         Element == "Mandibular" & est.type == "adj.SCEM"), 
             aes(x=sea.est, y=1, shape = est.type), 
             size = 1.25, fill = "white", colour = "sienna3") +
  scale_colour_manual(values = c("sienna3","gray50"), name = "Cosine fit") +
  scale_fill_manual(values = c("sienna3","gray50"), name = "Cosine fit") +
  scale_shape_manual(values = c(21, 24), guide = "none") +
  geom_text(data=Seasonlong, aes(x=sea.est, y=1.25, label=ID), 
            color="gray10", size=1.75) + 
  coord_polar(direction = -1, start = -1.57) + 
  scale_x_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75)) +
  scale_y_continuous(limits = c(0, 1.25), breaks = c(0, 1)) +
  facet_grid(.~est.type2) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "left",
    legend.key = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    panel.grid  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
# Figure 09 - seasonality estimates
ggsave("Fig09b_Seasonplot.pdf", units = "in", width = 4.25, height = 1.8)

###  Strontium -----------------------------------------------------------------
### .Figure 10  -----
# Figure 10 - 87Sr/86Sr values from the Basentello Valley estate
# load files: 
# Basentello Valley Sr data and comparative data from the Vagnari cemetery
input.Sr <- c("BVestate_enamel_Strontium.csv")
input.Sr.comp <- c("Comparative_enamel_VagnariCemetery.csv")
# load files
srResults <- read.csv(input.Sr)
srResultsComp <- read.csv(input.Sr.comp)
# set groups and order
srResults$Group <- "study"
srResultsComp$Group <- "comparative"
srResults <- srResults %>%
  bind_rows(srResultsComp)
srResults$x.axis <- NA
srResults <- srResults %>%
  mutate(x.axis = case_when(Group == "study" ~ ToothTaxon,
                            Group == "comparative" ~ Material,
                            TRUE ~ as.character(x.axis)))
order <- c(
  "SF01_Ovis",
  "SF06_Ovis",
  "SF08_Ovis",
  "SF04_Ovis",
  "SF05_Capra",
  "SF07_Capra",
  "SF03_Capra",
  "VV18_Ovis",
  "VV11_Ovis",
  "VV19_Ovis",
  "VV16_Ovis",
  "VV12_Ovis",
  "VV10_Capra",
  "VV09_Capra",
  "VV13_Capra",
  "Fauna",
  "Snail shell",
  "Human",
  "Soil")
srResults$ToothTaxon <- factor(x = srResults$ToothTaxon, 
                             levels = order)
srResults$x.axis<- factor(x = srResults$x.axis, 
                               levels = order)
srResults <- srResults %>%
  mutate(SampleLoc = Sr.Sample) %>%
  mutate(SampleLoc = dplyr::recode(SampleLoc, A = "bottom", B = "mid", C = "top"))

srResults$SampleLoc <- factor(x = srResults$SampleLoc, 
                              levels = c("top", "mid", "bottom"))
# Summary stats
Stats <- subset(srResults, x.axis == "Human" | x.axis == "Fauna")
boxplot.stats(Stats$Sr)
boxplot(Stats$Sr)
quantile(Stats$Sr)
Stats %>%
  summarise(mean = mean(Sr),
            sd = sd(Sr),
            median = median(Sr))
# df with 1.5 IQR
localRange <- data.frame(
  from=as.numeric(c('0.7083175')), 
  to=as.numeric(c('0.708828')),
  range=c('1.5IQR'))
srPlot <- ggplot() +
  geom_rect(data = localRange, aes(xmin = -Inf, xmax = Inf, 
                                   ymin = from, ymax = to), alpha = 0.1) +
  geom_point(data = srResults, 
             aes(x.axis, Sr, fill= SampleLoc, 
                 shape = SampleLoc),  stroke = 0.5, size=1.5) +
  scale_colour_manual(na.translate = F, name = "Comparative materials", 
                      values = c("springgreen3", "turquoise4", "palevioletred4")) +
  scale_shape_manual(na.translate = F, values = c(21, 24, 22), 
                     name = "Intra-tooth sample") +
  scale_fill_manual(na.translate = F, 
                    values = c("coral2", "goldenrod", "turquoise3"), 
                    name = "Intra-tooth sample") +
  geom_jitter(data = subset(srResults, Group == "comparative"), 
              aes(x.axis, Sr, colour = x.axis, size = .8), 
              alpha = 0.5, width=.3, size=1.1, shape = 19) +
  geom_boxplot(data = subset(srResults, Group == "comparative"), 
               aes(x.axis, Sr),  notch = TRUE, outlier.alpha = 0, 
               colour = "gray20", width = .9, alpha = 0, lwd = 0.25) +
  labs(x = "", y=expression({}^87*"Sr/"^86*"Sr")) +
  scale_y_continuous(breaks=seq(0.707,0.710,0.0005)) +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  guides(color = "none")
srPlot
ggsave(file="Figure10_intratooth_Sr.pdf", height = 90, width = 136, units = "mm")