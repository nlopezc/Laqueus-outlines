#Outline analysis for CT-imaged dataset (Dataset 1)

library(Momocs)
library(geomorph)
library(Morpho)
library(ggplot2)
library(dplyr)

ctout.lf <- list.files('~/Dataset1_OutlineJPEGS', full.names = T)
CToutlines.xy <-import_jpg(ctout.lf)
#factors
laqueussp <- as.factor(c(rep("blanfordi",1), rep("erythraeus",16), rep("quadratus",2),
                         rep("rubellus",9), rep("vancouveriensis",12))) #all species
laqueussp_eqrv <- as.factor(c(rep("erythraeus",16), rep("quadratus",2),
                         rep("rubellus",9), rep("vancouveriensis",12))) #without L. blanfordi
laqueussp_erv <- as.factor(c(rep("erythraeus",16), rep("rubellus",9), rep("vancouveriensis",12))) #without L. blanfordi and L. quadratus

#Out object
LaqueusCTOut <-Out(CToutlines.xy, fac = laqueussp) %>%
  coo_sample(800)
stack(LaqueusCTOut) # or pile()
#add landmarks
LaqueusCTOutlines <- LaqueusCTOut
for(i in seq_along(1:12)) {
  LaqueusCTOutlines <- def_ldk_angle(LaqueusCTOutlines, angle = 0.52*i)
}
stack(LaqueusCTOutlines)

#panel showing all outlines
panel(LaqueusCTOutlines, dim = c(8,5), fac = laqueussp, cols  = c("#feb24c", "#31a354", "#dd1c77", "#9e9ac8", "#c2e699"))

#Centroid size
CSoutlines <- coo_centsize(LaqueusCTOutlines)

#GPA- Full Procrustes
LaqueusCTOutlines.al <- fgProcrustes(LaqueusCTOutlines, tol = 1e-100) %>%
  coo_slide(ldk=3) #defines first point at the umbo 
stack(LaqueusCTOutlines.al) #to visualize superimposition

#EFA
LaqueusCTOutlines.ef <- efourier(LaqueusCTOutlines.al, norm = FALSE) #norm = F because we already did GPA
#PCA
LaqueusCTOutlines.pca <- PCA(LaqueusCTOutlines.ef)
#Standard Momocs plots
plot_PCA(LaqueusCTOutlines.pca, f = laqueussp, morphospace_position = "xy", axesnames = F)
plot_PCA(LaqueusCTOutlines.pca, axes = c(1,3), f = laqueussp, morphospace_position = "xy", axesnames = F)

#data prep for ggplot
LaqueusCTOutlines.pca$x
OutCTLdk.pca.df <- data.frame(LaqueusCTOutlines.pca$x)
OutCTLdk.pca.df$Sp <- factor(laqueussp,
                             labels = c("blanfordi"," erythraeus","quadratus","rubellus",'vancouveriensis'))
LaqueusCTOutlines.pca$eig
outCTldk.pcval <- as.matrix(LaqueusCTOutlines.pca$eig)
outCTldk.pcval1 <- round(outCTldk.pcval*100, digits=2)
#convex hulls
library(plyr)
hullsPC1_2 <- ddply(OutCTLdk.pca.df, .(Sp), function(df) df[chull(df$PC1, df$PC2),])
hullsPC1_3 <- ddply(OutCTLdk.pca.df, .(Sp), function(df) df[chull(df$PC1, df$PC3),])

##Plots
#PC1-2
OutCTLdkPC1_2 <- ggplot(data = OutCTLdk.pca.df, aes(x = PC1, y = PC2, color = Sp))+
  geom_point(size = 4, alpha = 0.95) +
  #geom_text(aes(x = PC1, y = PC2, label = rownames(OutCTLdk.pca.df))) +
  labs(x = paste("PC1 (",outCTldk.pcval1[1,1], "% )"),
       y = paste("PC2 (",outCTldk.pcval1[2,1], "% )")) +
  scale_colour_manual(name="",
                      labels= c("L. blanfordi","L. erythraeus","L. quadratus","L. rubellus","L. vancouveriensis", "T. coreanica", "T. transversa"),
                      values=c("#feb24c", "#31a354", "#dd1c77", "#9e9ac8", "#c2e699"))+
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 4, alpha = 0.95, shape = 21, colour = "black") +
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  geom_polygon(data = hullsPC1_2, aes(x=PC1, y=PC2, fill=Sp), alpha=0.15) +
  scale_fill_manual(values=c("#deebf7", "#31a354", "#2c7fb8", "#9e9ac8", "#c2e699"))

#PC1-3
OutCTLdkPC1_3 <- ggplot(data = OutCTLdk.pca.df, aes(x = PC1, y = PC3, color = Sp)) +
  geom_point(size = 4, alpha = 0.95) +
  #geom_text(aes(x = PC1, y = PC3, label = rownames(OutCTLdk.pca.df))) +
  labs(x = paste("PC1 (",outCTldk.pcval1[1,1], "% )"),
       y = paste("PC3 (",outCTldk.pcval1[3,1], "% )")) +
  scale_colour_manual(name="",
                      labels= c("L. blanfordi","L. erythraeus","L. quadratus","L. rubellus","L. vancouveriensis", "T. coreanica", "T. transversa"),
                      values=c("#feb24c", "#31a354", "#dd1c77", "#9e9ac8", "#c2e699"))+
  geom_point(size = 4, alpha = 0.95, shape = 21, colour = "black") +
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  theme(legend.position="none")+
  theme(text = element_text(size = 13))+
  theme(legend.title=element_blank())+
  geom_polygon(data = hullsPC1_3, aes(x=PC1, y=PC3, fill=Sp), alpha=0.15) +
  scale_fill_manual(values=c("#deebf7", "#31a354", "#2c7fb8", "#9e9ac8", "#c2e699"))


Fig_5ab <- cowplot::plot_grid(OutCTLdkPC1_2,
                           OutCTLdkPC1_3,
                           ncol = 2,
                           align = "hv")
#Shape change
CTPC1min <- coo_plot(LaqueusCTOutlines.al$coo$erythraeus_DAVSJCLab_0007, col = "#cfd2d3")
CTPC1max <- coo_plot(LaqueusCTOutlines.al$coo$rubellus_USNM_PAL_716075, col = "#cfd2d3")
CTPC2min <- coo_plot(LaqueusCTOutlines.al$coo$van_USNM_PAL_716058, col = "#cfd2d3")
CTPC2max <- coo_plot(LaqueusCTOutlines.al$coo$quadratus_USNM_PAL_716076, col = "#cfd2d3")
CTPC3min <- coo_plot(LaqueusCTOutlines.al$coo$van_USNM_PAL_716065, col = "#cfd2d3")
CTPC3max <- coo_plot(LaqueusCTOutlines.al$coo$van_USNM_PAL_716058, col = "#cfd2d3")
CTblanfordi <- coo_plot(LaqueusCTOutlines.al$coo$blanfordi_USNM_PAL_716078, col = "#cfd2d3")

##Thin-plate spline deformation vector fields
CTtpsPC1 <- tps_arr(LaqueusCTOutlines.al$coo$erythraeus_DAVSJCLab_0007,LaqueusCTOutlines.al$coo$rubellus_USNM_PAL_716075) #Fig.5C
CTtpsPC2 <- tps_arr(LaqueusCTOutlines.al$coo$van_USNM_PAL_716058,LaqueusCTOutlines.al$coo$quadratus_USNM_PAL_716076) #Fig.5D
CTtpsPC3 <- tps_arr(LaqueusCTOutlines.al$coo$van_USNM_PAL_716065,LaqueusCTOutlines.al$coo$van_USNM_PAL_716058) #Fig.5E
##Mean shapes--Fig.6
CTmeanshape <- mshapes(LaqueusCTOutlines.ef, fac = laqueussp)
CTmeanshapeery <- coo_plot(CTmeanshape$shp$erythraeus)
CTmeanshapequad <- coo_plot(CTmeanshape$shp$quadratus)
CTmeanshaperub <- coo_plot(CTmeanshape$shp$rubellus)
CTmeanshapevan <- coo_plot(CTmeanshape$shp$vancouveriensis)

## CVA prep
Momocs::export(LaqueusCTOutlines.ef) ## file written: coefficients.txt

ef.coefs <- read.table("coefficients.txt", header = T)
ef.coefs.df <- as.data.frame(ef.coefs) 
rownames(ef.coefs.df) <- ef.coefs.df[,1]
ef.coefs.df <- subset(ef.coefs.df, select = -c(name,value))
ef.coefs.df_no_blanf <-ef.coefs.df[-c(1),] #removing blanfordi

#CVA of EF coefficients--L. erythraeus, L. quadratus, L. rubellus, L. vancouveriensis
ct_out_CVA <- CVA(ef.coefs.df_no_blanf, groups = laqueussp_eqrv, rounds = 100000)
ct_out_classify <- classify(ct_out_CVA) #Fig.7C

ct_out_CVA$CVscores #individual canonical variate scores
ct_out_CVA$Var #variance explained by canonical variates
ct_out_CVA.df <- as.data.frame(ct_out_CVA$CVscores)
ct_out_CVA.df$Species <- factor(laqueussp_eqrv,
                                     labels = c("Laqueus erythraeus", "Laqueus quadratus", "Laqueus rubellus", "Laqueus vancouveriensis"))
ct_out_CVA.var <- (ct_out_CVA$Var)[,2]
ct_out_CVA.var <- round(ct_out_CVA.var, digits = 2)

hullsCVA1_2 <- ddply(ct_out_CVA.df, .(Species), function(df) df[chull(df$V1, df$V2),])

## CVA plot--L. erythraeus, L. quadratus, L. rubellus, L. vancouveriensis
ct_out_CVA.gg <- ggplot(data = ct_out_CVA.df, aes(x = V1, y = V2, fill = Species)) +
  geom_point(size = 4, alpha = 0.95)+
  labs(x = paste("CV1 (", ct_out_CVA.var[1], "% )"),
       y = paste("CV2 (", ct_out_CVA.var[2], "% )")) +
  scale_colour_manual(name="",
                      labels= c("L. erythraeus","L. quadratus","L. rubellus","L. vancouveriensis"),
                      values=c("#31a354", "#dd1c77", "#9e9ac8", "#c2e699"))+
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 4, alpha = 0.95, shape = 21, colour = "black") +
  theme(text = element_text(size = 13)) +
  theme(legend.position="none")+
  geom_polygon(data = hullsCVA1_2, aes(x=V1, y=V2, fill=Species, color = Species), alpha=0.5) +
  scale_fill_manual(values=c("#31a354", "#dd1c77", "#9e9ac8", "#c2e699"))

#CVA without L. quadratus
ef.coefs.df_no_quadratus <-ef.coefs.df[-c(1,18,19),]
ct_out_CVA_erv <- CVA(ef.coefs.df_no_quadratus, groups = laqueussp_erv, rounds = 100000)
ct_out_CVA_erv$CVscores #individual canonical variate scores
ct_out_CVA_erv$Var #variance explained by canonical variates
ct_out_CVA_erv.df <- as.data.frame(ct_out_CVA_erv$CVscores)
ct_out_CVA_erv.df$Species <- factor(laqueussp_erv,
                                labels = c("Laqueus erythraeus", "Laqueus rubellus", "Laqueus vancouveriensis"))
ct_out_CVA_erv.var <- (ct_out_CVA_erv$Var)[,2]
ct_out_CVA_erv.var <- round(ct_out_CVA_erv.var, digits = 2)

hullsCVA1_2_erv <- ddply(ct_out_CVA_erv.df, .(Species), function(df) df[chull(df$V1, df$V2),])

## CVA plot--L. erythraeus, L. rubellus, L. vancouveriensis
ct_out_CVA_erv.gg <- ggplot(data = ct_out_CVA_erv.df, aes(x = V1, y = V2, fill = Species)) +
  geom_point(size = 4, alpha = 0.95)+
  labs(x = paste("CV1 (", ct_out_CVA_erv.var[1], "% )"),
       y = paste("CV2 (", ct_out_CVA_erv.var[2], "% )")) +
  scale_colour_manual(name="",
                      labels= c("L. erythraeus","L. rubellus","L. vancouveriensis"),
                      values=c("#31a354", "#9e9ac8", "#c2e699"))+
  theme_classic() + 
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 4, alpha = 0.95, shape = 21, colour = "black") +
  theme(text = element_text(size = 13)) +
  theme(legend.position="none")+
  geom_polygon(data = hullsCVA1_2_erv, aes(x=V1, y=V2, fill=Species, color = Species), alpha=0.5) +
  scale_fill_manual(values=c("#31a354", "#9e9ac8", "#c2e699"))

#Fig.7 CVA of CT data
Fig_7ab <- cowplot::plot_grid(ct_out_CVA.gg,
                              ct_out_CVA_erv.gg,
                              ncol = 2,
                              align = "hv")

####

#CVA of PC scores (for comparison vs EF coefficients)

OutCTLdk.PCscores.df <- data.frame(LaqueusCTOutlines.pca$x)
dataset1.PCscores_EQRV <- OutCTLdk.PCscores.df[-1,]
dataset1.CVA_pcscores_EQRV <- CVA(dataset1.PCscores_EQRV, groups = laqueussp_eqrv, rounds = 100000) # without blanfordi 
classify(dataset1.CVA_pcscores_EQRV) # same results as Fig.7C

dataset1.PCscores_ERV <- dataset1.PCscores_EQRV[-17,]
dataset1.PCscores_ERV <- dataset1.PCscores_ERV[-17,]
dataset1.CVA_pcscores_ERV <- CVA(dataset1.PCscores_ERV, groups = laqueussp_erv, rounds = 100000) #without blanfordi and quadratus
classify(dataset1.CVA_pcscores_ERV)

#graph CVA from PC scores, L. erythraeus, quadratus, rubellus, vancouveriensis
dataset1.CVA_pcscores_EQRV$CVscores
dataset1.CVA_pcscores.df <- as.data.frame(dataset1.CVA_pcscores_EQRV$CVscores)
dataset1.CVA_pcscores.df$Species <- factor(laqueussp_eqrv,
                                labels = c("Laqueus erythraeus", "Laqueus quadratus", "Laqueus rubellus", "Laqueus vancouveriensis"))
ct_out_CVA.var1 <- (dataset1.CVA_pcscores_EQRV$Var)[,2]
ct_out_CVA.var1 <- round(ct_out_CVA.var1, digits = 2)

hullsCVA1_2_pcscores <- ddply(dataset1.CVA_pcscores.df, .(Species), function(df) df[chull(df$V1, df$V2),])

ct_out_CVA.pca.graph <- ggplot(data = dataset1.CVA_pcscores.df, aes(x = V1, y = V2, fill = Species)) +
  geom_point(size = 4, alpha = 0.95)+
  labs(x = paste("CV1 (", ct_out_CVA.var1[1], "% )"),
       y = paste("CV2 (", ct_out_CVA.var1[2], "% )")) +
  scale_colour_manual(name="",
                      labels= c("L. erythraeus","L. quadratus","L. rubellus","L. vancouveriensis"),
                      values=c("#31a354", "#dd1c77", "#9e9ac8", "#c2e699"))+
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 4, alpha = 0.95, shape = 21, colour = "black") +
  theme(text = element_text(size = 13)) +
  theme(legend.position="none")+
  geom_polygon(data = hullsCVA1_2_pcscores, aes(x=V1, y=V2, fill=Species, color = Species), alpha=0.5) +
  scale_fill_manual(values=c("#31a354", "#dd1c77", "#9e9ac8", "#c2e699"))

#graph CVA from PC scores, L. erythraeus, rubellus, vancouveriensis
dataset1.CVA_pcscores_ERV.df <- as.data.frame(dataset1.CVA_pcscores_ERV$CVscores)
dataset1.CVA_pcscores_ERV.df$Species <- factor(laqueussp_erv,
                                    labels = c("Laqueus erythraeus", "Laqueus rubellus", "Laqueus vancouveriensis"))
ct_out_CVA_erv.var1 <- (dataset1.CVA_pcscores_ERV.df$Var)[,2]
ct_out_CVA_erv.var1 <- round(ct_out_CVA_erv.var, digits = 2)

hullsCVA1_2_erv1 <- ddply(dataset1.CVA_pcscores_ERV.df, .(Species), function(df) df[chull(df$V1, df$V2),])

## CVA plot--L. erythraeus, L. rubellus, L. vancouveriensis
ct_out_CVA.pca.graph2 <- ggplot(data = dataset1.CVA_pcscores_ERV.df, aes(x = V1, y = V2, fill = Species)) +
  geom_point(size = 4, alpha = 0.95)+
  labs(x = paste("CV1 (", ct_out_CVA_erv.var1[1], "% )"),
       y = paste("CV2 (", ct_out_CVA_erv.var1[2], "% )")) +
  scale_colour_manual(name="",
                      labels= c("L. erythraeus","L. rubellus","L. vancouveriensis"),
                      values=c("#31a354", "#9e9ac8", "#c2e699"))+
  theme_classic() + 
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 4, alpha = 0.95, shape = 21, colour = "black") +
  theme(text = element_text(size = 13)) +
  theme(legend.position="none")+
  geom_polygon(data = hullsCVA1_2_erv1, aes(x=V1, y=V2, fill=Species, color = Species), alpha=0.5) +
  scale_fill_manual(values=c("#31a354", "#9e9ac8", "#c2e699"))

Fig_CVAonPCscores <- cowplot::plot_grid(ct_out_CVA.pca.graph,
                              ct_out_CVA.pca.graph2,
                              ncol = 2,
                              align = "hv")


## 3D geometric morphometric analysis of long loops

#Import landmark data
LaqueusLoop40 <- read.morphologika("~/Dataset1_landmarkcoords_raw.txt")
#define.sliders(LaqueusLoop40[,,1], 138, surfsliders=NULL, write.file=TRUE) 
curvesldks <- as.matrix(read.csv("~/curveslide6.csv", header=TRUE))

#GPA
laqueusBEQRV.gpa <- gpagen(LaqueusLoop40,curves=curvesldks, ProcD=FALSE) ## ProcD=False to minimize bending energy
plot(laqueusBEQRV.gpa)

# Morphological Integration test
Outlinecoords.al <- read.morphologika('~/Dataset1_outlinecoords_aligned.txt') 
dim(Outlinecoords.al)
dim(laqueusBEQRV.gpa$coords)

intTest <- integration.test(laqueusBEQRV.gpa$coords, Outlinecoords.al, iter = 1000) 
summary(intTest)
plot(intTest, label = laqueussp, warpgrids = T) # basically plot(intTest$XScores[,1], intTest$YScores[,1])

PLS.df <- as.data.frame(intTest$XScores[,1])
colnames(PLS.df)[1] <- "Xscores"
PLS.df$Yscores <- intTest$YScores[,1]
PLS.df$Species <- laqueussp
#Fig.8
LaqueusPLS <- ggplot(data = PLS.df, aes(x = Xscores, y = Yscores, color = Species)) +
  geom_point(size = 5, alpha = 0.85) +
  #geom_text(aes(x = PC1, y = PC2, label = rownames(LaqueusLoopsOut.pca.df))) +
  labs(x = "PLS1 block 1 (long-looped brachidia)",
       y = "PLS1 block 2 (outlines)") +
  scale_colour_manual(name="",
                      labels= c("L. blanfordi","L. erythraeus","L. quadratus","L. rubellus","L. vancouveriensis", "T. coreanica", "T. transversa"),
                      values=c("#feb24c", "#31a354", "#dd1c77", "#9e9ac8", "#c2e699"))+
  theme_classic() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.2, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.2, linetype = "solid", colour = "black"))+
  geom_point(size = 5, alpha = 0.95, shape = 21, colour = "black") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none")
