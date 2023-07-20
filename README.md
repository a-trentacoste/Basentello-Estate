# Basentello-Estate
## Multi-isotope analysis of fauna from the Roman imperial estate in the Basentello Valley (Puglia, Italy)
by Angela Trentacoste
July 2023

Analysis script for the paper: <br>
Isotopic insights into livestock production in Roman Italy: diet, seasonality, and mobility on an imperial estate <br>
Angela Trentacoste, Michael Mackinnon, Christopher Day, Petrus Le Roux, Mike Buckley, Myles McCallum, and Maureen Carroll <br>
under review

<b>If you use or adapt these files please cite the paper and the data.</b>

## Abstract
Agriculture is the most important intersection between farming communities and the natural world,  
with major implications for land exploitation and labour organisation. In Italy, at the heart of the Roman Empire, 
understanding of agriculture remains heavily dependant on ancient sources, which are unable to provide a regional or 
diachronic view of practices across the socio-economic spectrum. In order to provide insight into agricultural economies 
in Roman Italy and their social and environmental implications, this article reconstructs agro-pastoral strategies at an
imperial estate in southern Italy through a multi-isotope investigation of livestock bone collagen and tooth enamel. 
Analysis of carbon, nitrogen, oxygen and strontium isotopes (δ13C, δ15N, δ18O, 86Sr/87Sr) are combined to evaluate animal management 
and mobility at Vagnari vicus and the villa of San Felice in the Basentello Valley. Results reveal taxon-specific herding strategies with 
the potential for significant inputs from legume forage/fodder and/or natural environments. Caprine herding did not appear to include 
long-distance transhumance. This analysis moves past previous text-based generalisations to provide a new and nuanced perspective on 
animal production in rural southern Italy and its economic and environmental implications.
Keywords: Roman economy, Imperial estates, animal husbandry, transhumance, agriculture, Puglia

## About
This repository contains the files needed for analysis and visualisation of the date presented in the paper, including quality control for bone collagen isotope values and 
enamel preservation (using ATR-FTIR), as well as a statistical workflow for interrogating stable isotope data. 
The repository includes data and analysis of carbon, nitrogen, oxygen and strontium isotopes (δ13C, δ15N, δ18O, 86Sr/87Sr) from two sites on the Basentello Imperial estate 
(Vagnari vicus and San Felice villa) and comparative data: <list>
- Isotopic results from collagen (Semchuk 2016) and tooth enamel (Emery et al. 2018) from humans in the Vagnari cemetery
- Comparative strontium isotope values from soil and fauna relevant to the Vagnari cemetery (Emery et al. 2018)
- Faunal collagen from Pompeii and Herculaneum (Soncin et al. 2021), Velia (Craig et al. 2009), Portus (O'Connell et al. 2019), Isola Sacra (Prowse 2001), and Etruscan central Italy (Trentacoste et al. 2020). </list>
- Data from the Global Network of Isotopes in Precipitation (GNIP) on precipitation amount and its oxygen isotope composition. </list>

Collagen isotope results were assessed for normality using a Shapiro–Wilk test and for 
homogeneity using Levene’s test. On the basis of these assessments, ANOVA and Kruskal-Wallis multivariate analyses were applied as appropriate. 
Differences between groups were tested using a Mann-Whitney U test. Intra-estate differences in particular taxa were assessed using a Students t-test. 

Isotopic niche space and overlap was estimated using the kernel utilization density method in the rKIN package at 50%, 75% and 95% contours (Eckrich et al. 2020).

Sequences of oxygen isotope values from tooth enamel carbonates (δ18Ocarb) were used to estimate birth seasonality at Vagnari vicus
based on the cosine model and non-parametric splitting–coalescence–estimation method (SCEM) presented in Chazin et al. (2019),
using the R script provided in the publication: Note that this requires the script file in Chazin et al (2019) Data S5. Supporting information.
These methods estimate birth seasonality using the position of maximum δ18Ocarb value scaled relative to the length of the tooth (Balasse et al. 2012).

## Contents
BasentelloValley_Analysis.R - annotated R script<br>
BVestate_collagen_SI.csv - collagen (C, N) isotope results and sample details for Vagnari and San Felice <br>
BVestate_enamel_SI.csv - tooth enamel carbonate istope results (O, C, Sr) and sample details for Vagnari and San Felice <br>
BVestate_enamel_Strontium.csv - strontium isotope ratios from tooth enamel from Vagnari and San Felice <br>
BVestate_FTIR.csv - FTIR peak heights from Vagnari and San Felice tooth enamel<br>
Comparative_collagenSI_ItalyFauna.csv - comparative C and N isotope values from bone collagen from other Etruscan and Roman fauna in Italy<br>
Comparative_collagenSI_VagnariHuman.csv - comparative C and N isotope values from bone collagen from humans from the Vagnari cemetery<br>
Comparative_enamel_VagnariCemetery.csv - comparative O, C, and Sr isotope values from tooth enamel from humans from the Vagnari cemetery<br>




## References
Balasse M, Obein G, Ughetto-Monfrin J, Mainland I (2012) Investigating seasonality and season of birth in past herds: a reference set of sheep enamel stable oxygen isotope ratios Archaeometry 54:349–368 doi: https://doi.org/10.1111/j.1475-4754.2011.00624.x

Chazin H, Deb S, Falk J, Srinivasan A (2019) New Statistical Approaches to Intra-individual Isotopic Analysis and Modelling of Birth Seasonality in Studies of Herd Animals Archaeometry 61:478-493 doi: https://doi.org/10.1111/arcm.12432

Craig OE et al. (2009) Stable isotopic evidence for diet at the Imperial Roman coastal site of Velia (1st and 2nd Centuries AD) in Southern Italy American Journal of Physical Anthropology 139:572-583 doi: https://doi.org/10.1002/ajpa.21021

Eckrich CA, Albeke SE, Flaherty EA, Bowyer RT, Ben-David M (2020) rKIN: Kernel-based method for estimating isotopic niche size and overlap Journal of Animal Ecology 89:757-771 doi: https://doi.org/10.1111/1365-2656.13159

Emery MV, Stark RJ, Murchie TJ, Elford S, Schwarcz HP, Prowse TL (2018) Mapping the origins of Imperial Roman workers (1st–4th century CE) at Vagnari, Southern Italy, using 87Sr/86Sr and δ18O variability American Journal of Physical Anthropology 166:837–850 doi: https://doi.org/10.1002/ajpa.23473

IAEA/WMO (2019) Global Network of Isotopes in Precipitation. The GNIP Database. https://nucleus.iaea.org/wiser

O'Connell TC et al. (2019) Living and dying at the Portus Romae Antiquity 93:719–734 doi: https://doi.org/10.15184/aqy.2019.64

Prowse T (2001) Isotopic and dental evidence for diet from the necropolis of Isola Sacra (1st-3rd centuries AD), Italy. PhD dissertation, McMaster University

Semchuk L (2016) A stable isotope investigation of diet at Vagnari. MA theis, School of Graduate Studies, McMaster University. http://hdl.handle.net/11375/20498

Soncin S et al. (2021) High-resolution dietary reconstruction of victims of the 79 CE Vesuvius eruption at Herculaneum by compound-specific isotope analysis Science Advances 7:eabg5791 doi: https://doi.org/10.1126/sciadv.abg5791

Trentacoste A, Lightfoot E, Le Roux P, Buckley M, Kansa SW, Esposito C, Gleba M (2020) Heading for the hills? A multi-isotope study of sheep management in first-millennium BC Italy Journal of Archaeological Science: Reports 29:102036 doi: https://doi.org/10.1016/j.jasrep.2019.102036
