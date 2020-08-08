#Outline analysis for Dataset 2

library(Momocs)
library(geomorph)
library(Morpho)
library(ggplot2)
library(dplyr)
library(vegan)

#import jpgs
out.lf <- list.files('~/Extant_Laqueus_OutlineJPEGS', full.names = T)
laqueusoutlines.xy <-import_jpg(out.lf)

splaqueus <- as.factor(c(rep("erythraeus",50), rep("vancouveriensis",49)))
laqueus_localities <- as.factor(c(rep("Catalina",25), rep("Monterey",25), 
                                  rep("JuandeFuca",2), rep("Unknown",1), 
                                  rep("PortAlice",1), rep("PortEtches",2),
                                  rep("Unknown",4), rep("St3450WA",2), 
                                  rep("PortAlice",15), rep("Unknown",4),
                                  rep("Shelikof",18))) #or modify depending on file order
#create out object
LaqueusOut <-Out(laqueusoutlines.xy, fac = splaqueus) %>%
  coo_sample(800)
stack(LaqueusOut) # or pile()

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

#Full Procrustes
##defining umbo as first point and sliding coords
LaqueusOutlines.gpa <- fgProcrustes(LaqueusOutlines, tol = 1e-100) %>%
  coo_slide(ldk=3) #defines first point as the umbo
stack(LaqueusOutlines.gpa)

#Partial Procrustes Analysis- no scaling
LaqueusOutlines.ppa <- part_procrustes_out(LaqueusOutlines)
stack(LaqueusOutlines.ppa)
coo_centsize(LaqueusOutlines.ppa) #centroid sizes same after partial Procrustes

#EFA & PCA
  #GPA dataset
LaqueusOutlines.gpa.ef <- efourier(LaqueusOutlines.gpa, norm = FALSE) #norm = F because we already did GPA
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
plot_PCA(LaqueusOutlines.ppa.pca, axes = c(1,3), f = splaqueus, morphospace_position = "xy", axesnames = F)

# mean shapes GPA
meanshape <- mshapes(LaqueusOutlines.gpa.ef, fac = splaqueus)
meanshape$Coe

meanshapeery <- coo_plot(meanshape$shp$erythraeus, col = "#cfd2d3")
meanshapevan <- coo_plot(meanshape$shp$vancouveriensis, col = "#cfd2d3")
# shape changes, Fig. 10
PC1min <- coo_plot(LaqueusOutlines.gpa$coo$C31, col = "#cfd2d3") #DAV:SJCLab C31
PC1max <- coo_plot(LaqueusOutlines.gpa$coo$LvanJ, col = "#cfd2d3") #USNM PAL 770910
PC2min <- coo_plot(LaqueusOutlines.gpa$coo$LvanG, col = "#cfd2d3") #USNM PAL 770907
PC2max <- coo_plot(LaqueusOutlines.gpa$coo$C31, col = "#cfd2d3") #DAV:SJCLab C31
PC3min <- coo_plot(LaqueusOutlines.gpa$coo$Lvan18, col = "#cfd2d3") #USNM PAL 770881
PC3max <- coo_plot(LaqueusOutlines.gpa$coo$Lvan2862K, col = "#cfd2d3") #USNM PAL 770892

tpsPC1 <- tps_arr(LaqueusOutlines.gpa$coo$C31,LaqueusOutlines.gpa$coo$LvanJ)
tpsPC2 <- tps_arr(LaqueusOutlines.gpa$coo$LvanG,LaqueusOutlines.gpa$coo$C31)
tpsPC3 <- tps_arr(LaqueusOutlines.gpa$coo$Lvan18,LaqueusOutlines.gpa$coo$Lvan2862K)

PCcontributions <- PCcontrib(LaqueusOutlines.gpa.pca, nax = 1:3)

#shapes from Partial Procrustes
PC1minp <- coo_plot(LaqueusOutlines.gpa$coo$LC2.50, col = "#cfd2d3") #DAV:SJCLab LC2.50
PC1maxp <- coo_plot(LaqueusOutlines.gpa$coo$LvanS, col = "#cfd2d3") # USNM PAL 770918
PC2minp <- coo_plot(LaqueusOutlines.gpa$coo$C33, col = "#cfd2d3") #DAV:SJCLab C33
PC2maxp <- coo_plot(LaqueusOutlines.gpa$coo$C31, col = "#cfd2d3") #DAV:SJCLab C31
PC3minp <- coo_plot(LaqueusOutlines.gpa$coo$C31, col = "#cfd2d3") #DAV:SJCLab C31
PC3maxp <- coo_plot(LaqueusOutlines.gpa$coo$LvanJ, col = "#cfd2d3") #USNM PAL 770910

# *PCA GRAPHS*

#data prep for ggplot PCA--Full Procrustes
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

hullsPCA_1_2 <- plyr::ddply(OutLdk.pca.df, .(Sp), function(df) df[chull(df$PC1, df$PC2),])
hullsPCA_1_3 <- plyr::ddply(OutLdk.pca.df, .(Sp), function(df) df[chull(df$PC1, df$PC3),])

#PCA--Fig. 9A
OutLdkPC1_2 <- ggplot(data = OutLdk.pca.df, aes(x = PC1, y = PC2, color = Sp)) +
  geom_point(size = 3, alpha = 0.9) +
  #geom_text(aes(x = PC1, y = PC2, label = rownames(OutLdk.pca.df))) +
  labs(x = paste("PC1 (",outldk.pcval1[1,1], "% )"),
       y = paste("PC2 (",outldk.pcval1[2,1], "% )")) +
  scale_colour_manual(name="Sp",
                      labels= c("erythraeus","vancouveriensis"),
                      values=c("#31a354", "#c2e699"))+
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 3, alpha = 0.9, shape = 21, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  theme(legend.text= element_text(face="italic"))+
  geom_polygon(data = hullsPCA_1_2, aes(x=PC1, y=PC2, fill=Sp), alpha=0.19) +
  scale_fill_manual(values=c("#31a354", "#c2e699"))

#PCA--Fig.9B
  OutLdkPC1_3<- ggplot(data = OutLdk.pca.df, aes(x = PC1, y = PC3, color = Sp)) +
  geom_point(size = 3, alpha = 0.9) +
  #geom_text(aes(x = PC1, y = PC3, label = rownames(OutLdk.pca.df))) +
  labs(x = paste("PC1 (",outldk.pcval1[1,1], "% )"),
       y = paste("PC3 (",outldk.pcval1[3,1], "% )")) +
  scale_colour_manual(name="Sp",
                      labels= c("erythraeus","vancouveriensis"),
                      values=c("#31a354", "#c2e699"))+
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 3, alpha = 0.9, shape = 21, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  theme(legend.text= element_text(face="italic"))+
  geom_polygon(data = hullsPCA_1_3, aes(x=PC1, y=PC3, fill=Sp), alpha=0.19) +
  scale_fill_manual(values=c("#31a354", "#c2e699"))

#Size-coded PCA--Fig. 11A
OutLdkPC1_2_size <- ggplot(data = OutLdk.pca.df.size, aes(x = PC1, y = PC2, fill = size, shape = sp)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_shape_manual(values = c(21,22))+
  viridis::scale_fill_viridis()+
  #geom_text(aes(x = PC1, y = PC2, label = rownames(OutLdk.pca.df))) +
  labs(x = paste("PC1 (",outldk.pcval1[1,1], "% )"),
       y = paste("PC2 (",outldk.pcval1[2,1], "% )")) +
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  #geom_point(size = 3, alpha = 0.9, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))


#PCA graph on with locality info--Fig. 12A
PC1_2localities <- ggplot(data = OutLdk.pca.df.LOC, aes(x = PC1, y = PC2, fill = loc, shape = sp)) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_shape_manual(values = c(21,22))+
  #geom_text(aes(x = PC1, y = PC2, label = rownames(OutLdk.pca.df.LOC))) +
  labs(x = paste("PC1 (",outldk.pcval1[1,1], "% )"),
       y = paste("PC2 (",outldk.pcval1[2,1], "% )")) +
  scale_fill_manual(name="loc",
                      values= c("#c51b8a","#c6e8d2", "#fa9fb5", "#41b6c4", "#d394ff", "#6c32b7", "#fec44f","gray50"))+
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 3.5, alpha = 0.9, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  theme(legend.text= element_text(face="italic"))
  #geom_polygon(data = hullsPCA_1_2, aes(x=PC1, y=PC2, fill=Sp), alpha=0.19) +
  #scale_fill_manual(values=c("#31a354", "#c2e699"))


###

#data prep for ggplot PCA--Partial Procrustes dataset
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

hullsPCA_1_2p <- plyr::ddply(LaqOut.ppa.pca.df, .(Sp), function(df) df[chull(df$PC1, df$PC2),])
hullsPCA_1_3p <- plyr::ddply(LaqOut.ppa.pca.df, .(Sp), function(df) df[chull(df$PC1, df$PC3),])

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
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  theme(legend.text= element_text(face="italic"))+
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
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  theme(legend.text= element_text(face="italic"))+
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
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  #geom_point(size = 3, alpha = 0.9, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))

#PCA graph on with locality info--Fig. 12B
LaqOut.ppa.pca.loc
PC1_2localities_pPr <- ggplot(data = LaqOut.ppa.pca.loc, aes(x = PC1, y = PC2, fill = loc, shape = sp)) +
  geom_point(size = 3.5, alpha = 0.9) +
  #geom_text(aes(x = PC1, y = PC2, label = rownames(OutLdk.pca.df.LOC))) +
  scale_shape_manual(values = c(21,22))+
  labs(x = paste("PC1 (",outldkppa.pcval1[1,1], "% )"),
       y = paste("PC2 (",outldkppa.pcval1[2,1], "% )")) +
  scale_fill_manual(name="loc",
                      values= c("#c51b8a","#c6e8d2", "#fa9fb5", "#41b6c4", "#d394ff", "#6c32b7", "#fec44f","gray50"))+
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 3.5, alpha = 0.9, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  theme(legend.text= element_text(face="italic"))
#geom_polygon(data = hullsPCA_1_2, aes(x=PC1, y=PC2, fill=Sp), alpha=0.19) +
#scale_fill_manual(values=c("#31a354", "#c2e699"))

fig.9 <- cowplot::plot_grid(OutLdkPC1_2,
                                 OutLdkPC1_3,
                                 OutLdkPC1_2p,
                                 OutLdkPC1_3p,
                           ncol = 2,
                           align = "v")

fig.12 <- cowplot::plot_grid(PC1_2localities,
                                    PC1_2localities_pPr,
                                    ncol = 2,
                                    align = "hv")

fig.11 <- cowplot::plot_grid(OutLdkPC1_2_size,
                               ppa.pca_size,
                               ncol=2,
                               align = "hv")


##  ##  ##

#LDA on elliptical Fourier coefficients

# GPA dataset
  # LDA [MOMOCS]
LaqueusOutLDA.efc <- LDA(LaqueusOutlines.gpa.ef, f=splaqueus, retain = 0.99)
classification_metrics(LaqueusOutLDA.efc)
LaqueusOutLDA.efc$mod.pred$x

#Density plot LD1 from full Procrustes dataset
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
  theme_classic() + #background of graph
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  theme(text = element_text(size = 15))+
  theme(legend.position="none")

classtable.gpa <- as.table(LaqueusOutLDA.efc$CV.tab)
classification.gpa <- matrix(c(40,10,22,27), ncol =2, byrow = T)
colnames(classification.gpa)<- c("erythraeus","vancouveriensis")
rownames(classification.gpa)<- c("erythraeus","vancouveriensis") 
class.gpa <- as.data.frame(classification.gpa)
class.gpa$species <- c("erythraeus","vancouveriensis")

table.gpa.lda <- ggpubr::ggtexttable(class.gpa, rows = NULL)

# Partial Procrustes dataset
  #LDA [MOMOCS]
LaqueusOut.ppa.lda <- LDA(LaqueusOutlines.ppa.ef, f=splaqueus, retain = 0.99)
classification_metrics(LaqueusOut.ppa.lda)
LaqueusOut.ppa.lda$mod.pred$x

#Density plot LD1 from partial Procrustes dataset
tmp <- as.data.frame(LaqueusOut.ppa.lda$mod.pred$x) %>% mutate(SampleID = row.names(.))
tmp$Species <- splaqueus

LaqueusOutLDA.efc$mod.pred$x

LDA_pproc <- ggplot(tmp, aes(LD1, fill = Species)) +
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
classification.ppa <- matrix(c(38,12,15,34), ncol =2, byrow = T)
colnames(classification.ppa)<- c("erythraeus","vancouveriensis")
rownames(classification.ppa)<- c("erythraeus","vancouveriensis") 
class.ppa <- as.data.frame(classification.ppa)
class.ppa$species <- c("erythraeus","vancouveriensis")

table.ppa.lda <- ggpubr::ggtexttable(class.ppa, rows = NULL)


fig.13 <- cowplot::plot_grid(LDA_gproc,
                              table.gpa.lda,
                              LDA_pproc,
                              table.ppa.lda,
                              ncol = 2,
                              align = "v")

## ## ##

#Stats-- same using Momocs and base stats
	#Momocs
laquesmanovatest <- MANOVA(LaqueusOutlines.gpa.ef, fac = splaqueus)
laquesmanovatestPPA <- MANOVA(LaqueusOutlines.ppa.ef, fac = splaqueus)
	#Base
LaqueusOutEFC.mx <- as.matrix(LaqueusOutlines.gpa.ef$coe)
LaqueusOutEFCp.mx <- as.matrix(LaqueusOutlines.ppa.ef$coe)
laqueusmanovaStats <- manova(LaqueusOutEFC.mx ~ splaqueus)
laqueusmanovaStats_p <- manova(LaqueusOutEFCp.mx ~ splaqueus)
summary(laqueusmanovaStats, test= "Hotelling-Lawley")
summary(laqueusmanovaStats_p, test= "Hotelling-Lawley")

eta_sq_gpa <- heplots::etasq(laqueusmanovaStats, anova=FALSE, partial=TRUE)
eta_sq_ppa <- heplots::etasq(laqueusmanovaStats_p, anova=FALSE, partial=TRUE)

## ## ##

## FOSSILS

outfossil.lf <- list.files('~/Documents/Laqueus/OUTLINES/Fossils/Fossil_JPEGS', full.names = T)
fossils.xy <-import_jpg(outfossil.lf)

laqueusfossils <- as.factor(c(rep("vancouveriensis",16)))

laqueusf.out <-Out(fossils.xy, fac = laqueusfossils) %>%
  coo_sample(800)
stack(laqueusf.out) # or pile
#add landmarks
laqueusfossils.out <- laqueusf.out
for(i in seq_along(1:12)) {
  laqueusfossils.out <- def_ldk_angle(laqueusfossils.out, angle = 0.52*i)
}
stack(laqueusfossils.out) # or pile
#panel showing all outlines
panel(laqueusfossils.out, dim = c(8,2), fac = laqueusfossils, cols ="#b0b0b0")
size_1 <- coo_centsize(laqueusfossils.out)

#Full Procrustes
laqueusfossils.al <- fgProcrustes(laqueusfossils.out, tol = 1e-100) %>%
  coo_slide(ldk=3) #defines first point as the umbo
stack(laqueusfossils.al)
coo_centsize(laqueusfossils.al) %>% as.data.frame()

#Partial Procrustes
fossils.ppa <- part_procrustes_out(laqueusfossils.out)
stack(fossils.ppa)
coo_centsize(fossils.ppa) %>% as.data.frame()

#EFA
Laqueusfossils.ef <- efourier(laqueusfossils.al, norm = FALSE) #GPA
Laqueusfossils.ef$coe
fossils.ppa.ef <- efourier(fossils.ppa, norm = FALSE) #PPA
fossils.ppa.ef$coe

fossilsef.df <- as.data.frame(Laqueusfossils.ef$coe)
fossilsef.df$size <- as.matrix(size_1)
fossilsef.df$species <- factor(laqueusfossils,
                               labels = "vancouveriensis")

extanttraining <- efourier(LaqueusOutlines.gpa, nb.h = 11, norm = FALSE)
extanttraining$coe
extantoutlines.df <- as.data.frame(extanttraining$coe)
sizeextant <- coo_centsize(LaqueusOutlines)
extantoutlines.df$size <- as.matrix(sizeextant)
extantoutlines.df$sp <- factor(splaqueus,
                               labels = c("erythraeus","vancouveriensis"))

library(randomForest)
set.seed(1001)

random.forest <- randomForest(sp ~., data = extantoutlines.df, ntree=5000)
importance(random.forest)
fossilpred_shapesize <- predict(random.forest, fossilsef.df) %>% as.data.frame()

## complete analysis with both extant and fossils--Shape only

ExtFoss.lf <- list.files('~/Documents/Laqueus/OUTLINES/Randomforests', full.names = T)
ExtFoss.xy <-import_jpg(ExtFoss.lf)

laqueusspextfoss <- as.factor(c(rep("L.erythraeus",50), rep("vancouveriensis",65)))

extfoss.out <-Out(ExtFoss.xy, fac = laqueusspextfoss) %>%
  coo_sample(800)
stack(extfoss.out) # or pile()
#add landmarks
laqueusextfoss.out <- extfoss.out
for(i in seq_along(1:12)) {
  laqueusextfoss.out <- def_ldk_angle(laqueusextfoss.out, angle = 0.52*i)
}
stack(laqueusextfoss.out) # or pile()

laqueusspextfoss.al <- fgProcrustes(laqueusextfoss.out, tol = 1e-100) %>%
  coo_slide(ldk=3) #defines first point as the umbo
stack(laqueusspextfoss.al)

laqueusspextfoss.efa <- efourier(laqueusspextfoss.al, norm = FALSE)
laqueusspextfoss.efa$coe

extfoss.df <- as.data.frame(laqueusspextfoss.efa$coe)
extfoss.df$species <- factor(laqueusspextfoss,
                             labels = c("l.erythraeus","vancouveriensis"))

train <- as.data.frame(extfoss.df[1:99,])
test <- as.data.frame(extfoss.df[100:115,])

set.seed(123)
FE.random.forest <- randomForest(species ~., data = train, ntree=5000)
importance(FE.random.forest)
fossilpred_prob <- predict(FE.random.forest, test, type = "prob") %>% as.data.frame()
fossilpred_shape <- predict(FE.random.forest, test) %>% as.data.frame()
