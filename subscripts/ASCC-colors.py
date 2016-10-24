# Brendan Reardon
# Van Allen Laboratory
# Dana-Farber Cancer Institute, Harvard Medical School
# The Broad Institute of MIT & Harvard
# 10 April 2016

# Genomic Evolution After Chemoradiotherapy in Anal Squamous Cell Carcinoma
# Colors Reference

## ===================================================
# We set the colors which will be referenced
# For the most part, we use The Tableau Public-published workbook: http://tableaufriction.blogspot.com/2012/11/finally-you-can-use-tableau-data-colors.html
# Some clinical annotations utilize color brewer 2: http://colorbrewer2.org/

## ===================================================
# Figure 1: Primary and Recurrent Co-Mut Plot. Tableau
colors_fig1_comut = [(230,230,230), (44,160,44), (31,119,180),(148,103,189),(214,39,40),
                    (255,127,14),(188,189,34),(227,119,194),(23,190,207), (230,230,230),(255,255,255)]
for i in range(len(colors_fig1_comut)):
    r, g, b = colors_fig1_comut[i]
    colors_fig1_comut[i] = (r / 255., g / 255., b / 255.)
    
cmap_fig1_comut = colors.ListedColormap(colors_fig1_comut)
bounds_fig1_comut = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
norm_fig1_comut = colors.BoundaryNorm(bounds_fig1_comut, cmap_fig1_comut.N) 

# Figure 1: Clinical annotation, response to chemotherapy. 
colors_tableauResp = [(94,60,153), (178,171,210), (255,255,255)]
for i in range(len(colors_tableauResp)):
    r, g, b = colors_tableauResp[i]
    colors_tableauResp[i] = (r / 255., g / 255., b / 255.)

cmapResp = colors.ListedColormap(colors_tableauResp)
boundsResp = [0,1,2,3]
normResp = colors.BoundaryNorm(boundsResp, cmapResp.N)

# Figure 1: Clinical annotation, recurrence distance. 
colors_fig1_distance = [(250,159,181),(197,27,138),(255,255,255)]
for i in range(len(colors_fig1_distance)):
    r, g, b = colors_fig1_distance[i]
    colors_fig1_distance[i] = (r / 256., g / 256., b / 256.)

cmapDistance = colors.ListedColormap(colors_fig1_distance)
boundsDistance = [0,1,2,3]
normDistance = colors.BoundaryNorm(boundsDistance, cmapDistance.N)  

# Figure 1: Clinical annotation, recurrence time, Color Brewer 2 
colors_fig1_time = [(173,221,142),(49,163,84),(255,255,255)]
#colors_fig1_time = [(253,187,132),(227,74,51),(255,255,255)]
for i in range(len(colors_fig1_time)):
    r, g, b = colors_fig1_time[i]
    colors_fig1_time[i] = (r / 256., g / 256., b / 256.)

cmapTime = colors.ListedColormap(colors_fig1_time)
boundsTime = [0,1,2,3]
normTime = colors.BoundaryNorm(boundsTime, cmapTime.N)  

# Figure 1: Clinical annotation, HPV status. Color brewer 2
colors_tableauHPV = [(255,221,113),(219,161,58), (199,199,199), (255, 255, 255)]
for i in range(len(colors_tableauHPV)):
    r, g, b = colors_tableauHPV[i]
    colors_tableauHPV[i] = (r / 256., g / 256., b / 256.)

cmapHPV = colors.ListedColormap(colors_tableauHPV)
boundsHPV = [0,1,2,3,4]
normHPV = colors.BoundaryNorm(boundsHPV, cmapHPV.N)

# Figure 1: Genomic annotation, 3q status. Tableau
colors_tableau3q = [(5,113,176),(114,158,206),(230,230,230),(237,102,93),(202,0,32),(255, 255, 255)]
for i in range(len(colors_tableau3q)):
    r, g, b = colors_tableau3q[i]
    colors_tableau3q[i] = (r / 255., g / 255., b / 255.)

cmap3q = colors.ListedColormap(colors_tableau3q)
bounds3q = [-2,-1,0,1,2,3,4]
norm3q = colors.BoundaryNorm(bounds3q, cmap3q.N)

# Figure 1: Genomic annotation, mutations per Mb. Tableau
colors_bar = [(140,86,75),(196,156,148)]
for i in range(len(colors_bar)):
    r, g, b = colors_bar[i]
    colors_bar[i] = (r / 255., g / 255., b / 255.)

## ===================================================
# Figure 2: Phylogeny
colors_fig2 = [(1, 0, 214), (127, 127, 127), (165, 0, 0)] 
for i in range(len(colors_fig2)):
    r, g, b = colors_fig2[i]
    colors_fig2[i] = (r / 256., g / 256., b / 256.)

## ===================================================
# Figure 3: Mutational and Neoantigen Burden. Tableau 
#colors_fig3 = [(0, 107, 164), (255, 128, 14), (216, 37, 38), (48, 147, 67)] 
colors_fig3 = [(31, 119, 180), (255, 127, 14), (255, 193, 86), (44, 160, 44)]
for i in range(len(colors_fig3)):
    r, g, b = colors_fig3[i]
    colors_fig3[i] = (r / 256., g / 256., b / 256.)

## ===================================================
# Supplementary Figure 4: Response and No Response. Color Brewer 2 
colors_cohort= [(161,217,155),(49,163,84),(229,245,224)]
for i in range(len(colors_cohort)):
    r, g, b = colors_cohort[i]
    colors_cohort[i] = (r / 256., g / 256., b / 256.)

cmap_cohort = colors.ListedColormap(colors_cohort)
bounds_cohort = [0,1,2,3]
norm_cohort = colors.BoundaryNorm(bounds_cohort, cmap_cohort.N)