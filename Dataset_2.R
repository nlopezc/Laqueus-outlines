#Outline analysis for Dataset 2

library(Momocs)
library(geomorph)
library(Morpho)
library(ggplot2)
library(dplyr)

#import jpgs
out.lf <- list.files('~/Extant_Laqueus_OutlineJPEGS', full.names = T)
laqueusoutlines.xy <-import_jpg(out.lf)
#factors
splaqueus <- as.factor(c(rep("erythraeus",50), rep("vancouveriensis",49)))
laqueus_localities <- as.factor(c(rep("Catalina",25), rep("Monterey",25), 
                                  rep("JuandeFuca",2), rep("Unknown",9), 
                                  rep("PortAlice",16), rep("PortEtches",2),
                                  rep("St3450WA",2),rep("Shelikof",18))) #or modify depending on file order
#create out object
LaqueusOut <-Out(laqueusoutlines.xy, fac = splaqueus) %>%
  coo_sample(800)
stack(LaqueusOut)
#add landmarks
LaqueusOutlines <- LaqueusOut
for(i in seq_along(1:12)) {
  LaqueusOutlines <- def_ldk_angle(LaqueusOutlines, angle = 0.52*i)
}
stack(LaqueusOutlines)

#panel showing all outlines
panel(LaqueusOutlines, fac = splaqueus, cols =c("#31a354", "#c2e699"))

#centroid size
centroidsize <- coo_centsize(LaqueusOutlines)

#GPA-- Generalized Procrustes Analysis
LaqueusOutlines.gpa <- fgProcrustes(LaqueusOutlines, tol = 1e-100) %>%
  coo_slide(ldk=3) #defines umbo as first point
stack(LaqueusOutlines.gpa)

## Partial Procrustes Analysis--i.e. no scaling
LaqueusOutlines.ppa <- part_procrustes_out(LaqueusOutlines) %>%
  coo_slide(ldk=3)
stack(LaqueusOutlines.ppa)
coo_centsize(LaqueusOutlines.ppa) #to verify size was not normalized

#EFA & PCA
  #GPA dataset
LaqueusOutlines.gpa.ef <- efourier(LaqueusOutlines.gpa, norm = FALSE) #norm = F because we already did alignment
LaqueusOutlines.gpa.ef$coe
LaqueusOutlines.gpa.pca <- PCA(LaqueusOutlines.gpa.ef)
##PC1-8 = 95.37% variance
plot_PCA(LaqueusOutlines.gpa.pca, f = splaqueus, morphospace_position = "xy", axesnames = F)
plot_PCA(LaqueusOutlines.gpa.pca, axes = c(1,3), f = splaqueus, morphospace_position = "xy", axesnames = F)

  #Partial Procrustes dataset
LaqueusOutlines.ppa.ef <- efourier(LaqueusOutlines.ppa, norm = FALSE)
LaqueusOutlines.ppa.ef$coe
LaqueusOutlines.ppa.pca <- PCA(LaqueusOutlines.ppa.ef)
plot_PCA(LaqueusOutlines.ppa.pca, f = splaqueus, morphospace_position = "xy", axesnames = F)
plot_PCA(LaqueusOutlines.ppa.pca, axes = c(1,3), f = splaqueus, morphospace_position = "range", axesnames = F)

###

### mean shapes GPA
meanshape <- mshapes(LaqueusOutlines.gpa.ef, fac = splaqueus)
meanshape$Coe

meanshapeery <- coo_plot(meanshape$shp$erythraeus, col = "#cfd2d3")
meanshapevan <- coo_plot(meanshape$shp$vancouveriensis, col = "#cfd2d3")

PC1min <- coo_plot(LaqueusOutlines.gpa$coo$erythraeus_DAVSJCLab_C31, col = "#cfd2d3") #DAV:SJCLab C31
PC1max <- coo_plot(LaqueusOutlines.gpa$coo$van_USNM_PAL_770910, col = "#cfd2d3") #USNM PAL 770910
PC2min <- coo_plot(LaqueusOutlines.gpa$coo$van_USNM_PAL_770907, col = "#cfd2d3") #USNM PAL 770907
PC2max <- coo_plot(LaqueusOutlines.gpa$coo$erythraeus_DAVSJCLab_C31, col = "#cfd2d3") #DAV:SJCLab C31
PC3min <- coo_plot(LaqueusOutlines.gpa$coo$van_USNM_PAL_770881, col = "#cfd2d3") #USNM PAL 770881
PC3max <- coo_plot(LaqueusOutlines.gpa$coo$van_USNM_PAL_770892, col = "#cfd2d3") #USNM PAL 770892

#Fig. 10
tpsPC1 <- tps_arr(LaqueusOutlines.gpa$coo$erythraeus_DAVSJCLab_C31,LaqueusOutlines.gpa$coo$van_USNM_PAL_770910)
tpsPC2 <- tps_arr(LaqueusOutlines.gpa$coo$van_USNM_PAL_770907,LaqueusOutlines.gpa$coo$erythraeus_DAVSJCLab_C31)
tpsPC3 <- tps_arr(LaqueusOutlines.gpa$coo$van_USNM_PAL_770881,LaqueusOutlines.gpa$coo$van_USNM_PAL_770892)

PCcontributions <- PCcontrib(LaqueusOutlines.gpa.pca, nax = 1:3)

#shapes from Partial Procrustes
PC1minp <- coo_plot(LaqueusOutlines.gpa$coo$erythraeus_DAVSJCLab_LC2.50, col = "#cfd2d3") #DAV:SJCLab LC2.50
PC1maxp <- coo_plot(LaqueusOutlines.gpa$coo$van_USNM_PAL_770918, col = "#cfd2d3") #USNM PAL 770918
PC2minp <- coo_plot(LaqueusOutlines.gpa$coo$erythraeus_DAVSJCLab_C31, col = "#cfd2d3") #DAV:SJCLab C31
PC2maxp <- coo_plot(LaqueusOutlines.gpa$coo$van_USNM_PAL_770873, col = "#cfd2d3") #USNM PAL 770873
PC3minp <- coo_plot(LaqueusOutlines.gpa$coo$van_USNM_PAL_770883, col = "#cfd2d3") #USNM PAL 770883
PC3maxp <- coo_plot(LaqueusOutlines.gpa$coo$van_USNM_PAL_770907, col = "#cfd2d3") #USNM PAL 770907

# PCA GRAPHS
#GPA
LaqueusOutlines.gpa.pca$x
OutLdk.pca.df <- data.frame(LaqueusOutlines.gpa.pca$x)
OutLdk.pca.df.LOC <- data.frame(LaqueusOutlines.gpa.pca$x)
OutLdk.pca.df.LOC$loc <- factor(laqueus_localities)
OutLdk.pca.df.LOC$sp <- factor(splaqueus,
                               labels = c("erythraeus","vancouveriensis"))
  
OutLdk.pca.df$Sp <- factor(splaqueus,
                            labels = c("erythraeus","vancouveriensis"))

OutLdk.pca.df.size <- data.frame(LaqueusOutlines.gpa.pca$x)
OutLdk.pca.df.size$size <- centroidsize
OutLdk.pca.df.size$sp <- factor(splaqueus,
                                labels = c("erythraeus","vancouveriensis"))
LaqueusOutlines.gpa.pca$eig
outldk.pcval <- as.matrix(LaqueusOutlines.gpa.pca$eig)
outldk.pcval1 <- round(outldk.pcval*100, digits=2)
#library(plyr)
hullsPCA_1_2 <- ddply(OutLdk.pca.df, .(Sp), function(df) df[chull(df$PC1, df$PC2),])
hullsPCA_1_3 <- ddply(OutLdk.pca.df, .(Sp), function(df) df[chull(df$PC1, df$PC3),])

#PCA--Fig. 9A
OutLdkPC1_2 <- ggplot(data = OutLdk.pca.df, aes(x = PC1, y = PC2, color = Sp)) +
  geom_point(size = 3, alpha = 0.9) +
  #geom_text(aes(x = PC1, y = PC2, label = rownames(OutLdk.pca.df))) +
  labs(x = paste("PC1 (",outldk.pcval1[1,1], "% )"),
       y = paste("PC2 (",outldk.pcval1[2,1], "% )")) +
  scale_colour_manual(name="Sp",
                      labels= c("erythraeus","vancouveriensis"),
                      values=c("#31a354", "#c2e699"))+
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 3, alpha = 0.9, shape = 21, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  geom_polygon(data = hullsPCA_1_2, aes(x=PC1, y=PC2, fill=Sp), alpha=0.19) +
  scale_fill_manual(values=c("#31a354", "#c2e699"))

#PCA--Fig. 9B
OutLdkPC1_3 <- ggplot(data = OutLdk.pca.df, aes(x = PC1, y = PC3, color = Sp)) +
  geom_point(size = 3, alpha = 0.9) +
  #geom_text(aes(x = PC1, y = PC3, label = rownames(OutLdk.pca.df))) +
  labs(x = paste("PC1 (",outldk.pcval1[1,1], "% )"),
       y = paste("PC3 (",outldk.pcval1[3,1], "% )")) +
  scale_colour_manual(name="Sp",
                      labels= c("erythraeus","vancouveriensis"),
                      values=c("#31a354", "#c2e699"))+
  theme_classic() + 
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 3, alpha = 0.9, shape = 21, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  geom_polygon(data = hullsPCA_1_3, aes(x=PC1, y=PC3, fill=Sp), alpha=0.19) +
  scale_fill_manual(values=c("#31a354", "#c2e699"))

#Size-coded PCA--Fig.11A
OutLdkPC1_2_size <- ggplot(data = OutLdk.pca.df.size, aes(x = PC1, y = PC2, fill = size, shape = sp)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_shape_manual(values = c(21,22))+
  viridis::scale_fill_viridis()+
  #geom_text(aes(x = PC1, y = PC2, label = rownames(OutLdk.pca.df))) +
  labs(x = paste("PC1 (",outldk.pcval1[1,1], "% )"),
       y = paste("PC2 (",outldk.pcval1[2,1], "% )")) +
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  theme(legend.position="none")+
  theme(text = element_text(size = 13))

#PCA graph with locality info--Fig. 12A
PC1_2localities <- ggplot(data = OutLdk.pca.df.LOC, aes(x = PC1, y = PC2, fill = loc, shape = sp)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_shape_manual(values = c(21,22))+
  #geom_text(aes(x = PC1, y = PC2, label = rownames(OutLdk.pca.df.LOC))) +
  labs(x = paste("PC1 (",outldk.pcval1[1,1], "% )"),
       y = paste("PC2 (",outldk.pcval1[2,1], "% )")) +
  scale_fill_manual(name="loc",
                      values= c("#c51b8a","#c6e8d2", "#fa9fb5", "#41b6c4", "#d394ff", "#6c32b7", "#fec44f","gray50"))+
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 3, alpha = 0.9, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())

##

#PPA
LaqueusOutlines.ppa.pca$x
LaqOut.ppa.pca.df <- data.frame(LaqueusOutlines.ppa.pca$x)

LaqOut.ppa.pca.loc <- data.frame(LaqueusOutlines.ppa.pca$x)
LaqOut.ppa.pca.loc$loc <- factor(laqueus_localities)
LaqOut.ppa.pca.loc$sp <- factor(splaqueus,
                                labels = c("erythraeus","vancouveriensis"))

LaqOut.ppa.pca.df$Sp <- factor(splaqueus,
                           labels = c("erythraeus","vancouveriensis"))

LaqOut.ppa.pca.size <- data.frame(LaqueusOutlines.ppa.pca$x)
LaqOut.ppa.pca.size$size <- centroidsize
LaqOut.ppa.pca.size$sp <- factor(splaqueus,
                                 labels = c("erythraeus","vancouveriensis"))

LaqueusOutlines.ppa.pca$eig
outldkppa.pcval <- as.matrix(LaqueusOutlines.ppa.pca$eig)
outldkppa.pcval1 <- round(outldkppa.pcval*100, digits=2)

hullsPCA_1_2p <- ddply(LaqOut.ppa.pca.df, .(Sp), function(df) df[chull(df$PC1, df$PC2),])
hullsPCA_1_3p <- ddply(LaqOut.ppa.pca.df, .(Sp), function(df) df[chull(df$PC1, df$PC3),])

#PCA--Fig.9C
OutLdkPC1_2p <- ggplot(data = LaqOut.ppa.pca.df, aes(x = PC1, y = PC2, color = Sp)) +
  geom_point(size = 3, alpha = 0.9) +
  #geom_text(aes(x = PC1, y = PC2, label = rownames(LaqOut.ppa.pca.df))) +
  labs(x = paste("PC1 (",outldkppa.pcval1[1,1], "% )"),
       y = paste("PC2 (",outldkppa.pcval1[2,1], "% )")) +
  scale_colour_manual(name="Sp",
                      labels= c("erythraeus","vancouveriensis"),
                      values=c("#31a354", "#c2e699"))+
  geom_point(size = 3, alpha = 0.9, shape = 21, colour = "black") +
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  geom_polygon(data = hullsPCA_1_2p, aes(x=PC1, y=PC2, fill=Sp), alpha=0.19) +
  scale_fill_manual(values=c("#31a354", "#c2e699"))

#PCA--Fig.9D
OutLdkPC1_3p <- ggplot(data = LaqOut.ppa.pca.df, aes(x = PC1, y = PC3, color = Sp)) +
  geom_point(size = 3, alpha = 0.9) +
  #geom_text(aes(x = PC1, y = PC3, label = rownames(LaqOut.ppa.pca.df))) +
  labs(x = paste("PC1 (",outldkppa.pcval1[1,1], "% )"),
       y = paste("PC3 (",outldkppa.pcval1[3,1], "% )")) +
  scale_colour_manual(name="Sp",
                      labels= c("erythraeus","vancouveriensis"),
                      values=c("#31a354", "#c2e699"))+
  geom_point(size = 3, alpha = 0.9, shape = 21, colour = "black") +
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  geom_polygon(data = hullsPCA_1_3p, aes(x=PC1, y=PC3, fill=Sp), alpha=0.19) +
  scale_fill_manual(values=c("#31a354", "#c2e699"))

#Size-coded PCA--Fig. 11B
ppa.pca_size <- ggplot(data = LaqOut.ppa.pca.size, aes(x = PC1, y = PC2, fill = size, shape = sp)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_shape_manual(values = c(21,22))+
  viridis::scale_fill_viridis()+
  #geom_text(aes(x = PC1, y = PC2, label = rownames(OutLdk.pca.df))) +
  labs(x = paste("PC1 (",outldkppa.pcval1[1,1], "% )"),
       y = paste("PC2 (",outldkppa.pcval1[2,1], "% )")) +
  theme_classic() + 
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  #geom_point(size = 3, alpha = 0.9, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))

#PCA graph on with locality info--Fig. 12B
LaqOut.ppa.pca.loc
PC1_2localities_pPr <- ggplot(data = LaqOut.ppa.pca.loc, aes(x = PC1, y = PC2, fill = loc, shape = sp)) +
  geom_point(size = 3, alpha = 0.9) +
  #geom_text(aes(x = PC1, y = PC2, label = rownames(OutLdk.pca.df.LOC))) +
  scale_shape_manual(values = c(21,22))+
  labs(x = paste("PC1 (",outldkppa.pcval1[1,1], "% )"),
       y = paste("PC2 (",outldkppa.pcval1[2,1], "% )")) +
  scale_fill_manual(name="loc",
                      values= c("#c51b8a","#c6e8d2", "#fa9fb5", "#41b6c4", "#d394ff", "#6c32b7", "#fec44f","gray50"))+
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 3, alpha = 0.9, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  theme(legend.text= element_text(face="italic"))
## ## ##

fig.9 <- cowplot::plot_grid(OutLdkPC1_2,
                                 OutLdkPC1_3,
                                 OutLdkPC1_2p,
                                 OutLdkPC1_3p,
                                 ncol = 2,
                                 align = "v")

fig.11 <- cowplot::plot_grid(OutLdkPC1_2_size,
                             ppa.pca_size,
                             ncol=2,
                             align = "hv")

fig.12 <- cowplot::plot_grid(PC1_2localities,
                                    PC1_2localities_pPr,
                                    ncol = 2,
                                    align = "hv")
## ## ##


#LDA

#GPA dataset
  # LDA on EF coeffs 
set.seed(123)
LaqueusOutLDA.efc <- LDA(LaqueusOutlines.gpa.ef, f=splaqueus, retain = 0.99)
classification_metrics(LaqueusOutLDA.efc)
LaqueusOutLDA.efc$mod.pred$x
predicted.efc.lda <- as.data.frame(LaqueusOutLDA.efc$CV.fac)

  #LDA on PC coeffs--verification
set.seed(123)
LDA_pcscores_gpa <- LDA(LaqueusOutlines.gpa.pca, f=splaqueus, retain = 64)
classification_metrics(LDA_pcscores_gpa) 
predicted.pcscores.gpa <- as.data.frame(LDA_pcscores_gpa$CV.fac)

#Partial Procrustes dataset
  #LDA on EF coeffs
set.seed(123)
LaqueusOut.ppa.lda <- LDA(LaqueusOutlines.ppa.ef, f=splaqueus, retain = 0.99)
classification_metrics(LaqueusOut.ppa.lda)
LaqueusOut.ppa.lda$mod.pred$x
predicted.efc.ppa.lda <- as.data.frame(LaqueusOut.ppa.lda$CV.fac)

  #LDA on PC coeffs--verification
set.seed(123)
LDA_pcscores_ppa <- LDA(LaqueusOutlines.ppa.pca, f=splaqueus, retain = 64)
classification_metrics(LDA_pcscores_ppa)
predicted.pca.ppa.lda <- as.data.frame(LDA_pcscores_ppa$CV.fac)

# Density plot LD1 for full Procrustes dataset
ldagen <- as.data.frame(LaqueusOutLDA.efc$mod.pred$x) %>% mutate(SampleID = row.names(.))
ldagen$Species <- splaqueus

LDA_gproc <- ggplot(ldagen, aes(LD1, fill = Species)) +
  geom_density(alpha = 0.8) +
  geom_rug(aes(color = Species)) +
  theme_bw() +
  scale_colour_manual(name="Sp",
                      labels= c("erythraeus","vancouveriensis"),
                      values=c("#31a354", "#c2e699"))+
  scale_fill_manual(values=c("#31a354", "#c2e699"))+
  theme_classic() + 
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  theme(text = element_text(size = 15))+
  theme(legend.position="none")

classtable.gpa <- as.table(LaqueusOutLDA.efc$CV.tab)
classification.gpa <- matrix(c(40,10,22,27), ncol =2, byrow = T) #based on LaqueusOutLDA.efc$CV.tab
colnames(classification.gpa)<- c("erythraeus","vancouveriensis")
rownames(classification.gpa)<- c("erythraeus","vancouveriensis") 
class.gpa <- as.data.frame(classification.gpa)
class.gpa$species <- c("erythraeus","vancouveriensis")

table.gpa.lda <- ggpubr::ggtexttable(class.gpa, rows = NULL)

# Density plot LD1 for partial Procrustes dataset
ldappa <- as.data.frame(LaqueusOut.ppa.lda$mod.pred$x) %>% mutate(SampleID = row.names(.))
ldappa$Species <- splaqueus

LDA_pproc <- ggplot(ldappa, aes(LD1, fill = Species)) +
  geom_density(alpha = 0.8) +
  geom_rug(aes(color = Species)) +
  theme_bw() +
  scale_colour_manual(name="Sp",
                      labels= c("erythraeus","vancouveriensis"),
                      values=c("#31a354", "#c2e699"))+
  scale_fill_manual(values=c("#31a354", "#c2e699"))+
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  theme(text = element_text(size = 15))+
  theme(legend.position="none")

LaqueusOut.ppa.lda$CV.tab
classification.ppa <- matrix(c(38,12,17,32), ncol =2, byrow = T) #based on LaqueusOut.ppa.lda$CV.tab
colnames(classification.ppa)<- c("erythraeus","vancouveriensis")
rownames(classification.ppa)<- c("erythraeus","vancouveriensis") 
class.ppa <- as.data.frame(classification.ppa)
class.ppa$species <- c("erythraeus","vancouveriensis")

table.ppa.lda <- ggpubr::ggtexttable(class.ppa, rows = NULL)

## ## ##

fig.13 <- cowplot::plot_grid(LDA_gproc,
                              #table.gpa.lda,
                              LDA_pproc,
                              #table.ppa.lda,
                              ncol = 2,
                              align = "vh")

#Supplementary LDA plots (on PC scores)

lda.pca.gpa <- as.data.frame(LDA_pcscores_gpa$mod.pred$x) %>% mutate(SampleID = row.names(.))
lda.pca.gpa$Species <- splaqueus

LDA_pcascores.gpa <- ggplot(lda.pca.gpa, aes(LD1, fill = Species)) +
  geom_density(alpha = 0.8) +
  geom_rug(aes(color = Species)) +
  theme_bw() +
  scale_colour_manual(name="Sp",
                      labels= c("erythraeus","vancouveriensis"),
                      values=c("#31a354", "#c2e699"))+
  scale_fill_manual(values=c("#31a354", "#c2e699"))+
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  theme(text = element_text(size = 15))+
  theme(legend.position="none")

lda.pca.ppa <- as.data.frame(LDA_pcscores_ppa$mod.pred$x) %>% mutate(SampleID = row.names(.))
lda.pca.ppa$Species <- splaqueus

LDA_pcascores.ppa <- ggplot(lda.pca.ppa, aes(LD1, fill = Species)) +
  geom_density(alpha = 0.8) +
  geom_rug(aes(color = Species)) +
  theme_bw() +
  scale_colour_manual(name="Sp",
                      labels= c("erythraeus","vancouveriensis"),
                      values=c("#31a354", "#c2e699"))+
  scale_fill_manual(values=c("#31a354", "#c2e699"))+
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  theme(text = element_text(size = 15))+
  theme(legend.position="none")

fig.supplementary1 <- cowplot::plot_grid(LDA_pcascores.gpa,
                             LDA_pcascores.ppa,
                             ncol = 2,
                             align = "hv") #same as Fig.13
## ## ##

#Stats-- same using Momocs and base stats
    #Momocs
laquesmanovatest <- MANOVA(LaqueusOutlines.gpa.ef, fac = splaqueus) # *
laquesmanovatestPPA <- MANOVA(LaqueusOutlines.ppa.ef, fac = splaqueus) # **
    #Base
LaqueusOutEFC.mx <- as.matrix(LaqueusOutlines.gpa.ef$coe)
LaqueusOutEFCp.mx <- as.matrix(LaqueusOutlines.ppa.ef$coe)
PCscores.gpa.mx <- as.matrix(LaqueusOutlines.gpa.pca$x)
PCscores.ppa.mx <- as.matrix(LaqueusOutlines.ppa.pca$x)

laqueusmanovaStats <- manova(LaqueusOutEFC.mx ~ splaqueus)
laqueusmanovaStats_p <- manova(LaqueusOutEFCp.mx ~ splaqueus)
#laqueusmanovaStats_pca_gpa <- manova(PCscores.gpa.mx ~ splaqueus) #on PC scores to verify
#laqueusmanovaStats_pca_ppa <- manova(PCscores.ppa.mx ~ splaqueus)

summary(laqueusmanovaStats, test= "Hotelling-Lawley") # *
summary(laqueusmanovaStats_p, test= "Hotelling-Lawley") # **
#summary(laqueusmanovaStats_pca_gpa, test= "Hotelling-Lawley")
#summary(laqueusmanovaStats_pca_ppa, test= "Hotelling-Lawley")

eta_sq_gpa <- heplots::etasq(laqueusmanovaStats, anova=FALSE, partial=TRUE) # *
eta_sq_ppa <- heplots::etasq(laqueusmanovaStats_p, anova=FALSE, partial=TRUE) # **
#eta_sq_pca_gpa <- heplots::etasq(laqueusmanovaStats_pca_gpa, anova=FALSE, partial=TRUE) #on PC scores to verify
#eta_sq_pca_ppa <- heplots::etasq(laqueusmanovaStats_pca_ppa, anova=FALSE, partial=TRUE)

## ## ##
library(randomForest)
#ExtFoss.lf <- list.files('~/Folder_containing_JPEGS_from_extant_and_fossil_specimens', full.names = T) #the following code 
                                                      #expects that the fossil specimens are the last ones on the list
ExtFoss.xy <-import_jpg(ExtFoss.lf)
laqueusspextfoss <- as.factor(c(rep("erythraeus",50), rep("vancouveriensis",65)))

extfoss.out <-Out(ExtFoss.xy, fac = laqueusspextfoss) %>%
  coo_sample(800)
stack(extfoss.out)
#add landmarks
laqueusextfoss.out <- extfoss.out
for(i in seq_along(1:12)) {
  laqueusextfoss.out <- def_ldk_angle(laqueusextfoss.out, angle = 0.52*i)
}
stack(laqueusextfoss.out)

##Full Procrustes, EFA, RF -- SHAPE
ext_fossils.gpa <- fgProcrustes(laqueusextfoss.out, tol = 1e-100) %>%
  coo_slide(ldk=3) 
stack(ext_fossils.gpa)

ext_fossils.gpa.efa <- efourier(ext_fossils.gpa, norm = FALSE) #16 harmonics
ext_fossils.gpa.efa$coe

ext_fossils.gpa.efa.df <- as.data.frame(ext_fossils.gpa.efa$coe)
ext_fossils.gpa.efa.df$species <- factor(laqueusspextfoss,
                                         labels = c("erythraeus","vancouveriensis"))

train1 <- as.data.frame(ext_fossils.gpa.efa.df[1:99,])
test1 <- as.data.frame(ext_fossils.gpa.efa.df[100:115,])

set.seed(1001)
rf1 <- randomForest(species ~., data = train1, ntree=5000)
fossilpred1_prob <- predict(rf1, test1, type = "prob") %>% as.data.frame()
fossilpred1_sp <- predict(rf1, test1) %>% as.data.frame() 

##Full Procrustes, EFA + centroid sizes, RF -- SHAPE + SIZE

ext_fossils.gpa.efa_cs.df <- ext_fossils.gpa.efa.df
ext_fossils.gpa.efa_cs.df$size <-as.matrix(coo_centsize(laqueusextfoss.out))

train2 <- as.data.frame(ext_fossils.gpa.efa_cs.df[1:99,])
test2 <- as.data.frame(ext_fossils.gpa.efa_cs.df[100:115,])

set.seed(1001)
rf2 <- randomForest(species ~., data = train2, ntree=5000)
fossilpred2_prob <- predict(rf2, test2, type = "prob") %>% as.data.frame()
fossilpred2_sp <- predict(rf2, test2) %>% as.data.frame()
