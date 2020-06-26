# Comparison of morphometric methods with 
# introduction to morphometric analyses 
# and the comparison of methods on coral reef fish
# By M. Quitzau, R. Frelat, V. Bonhomme
# Last update: 26th June 2020

# For more information, see https://rfrelat.github.io/CoralFishes.html

# 1. Load and explore the dataset ---------------
library("Momocs")

load("CoralFishes.Rdata")

table(traits$family)

# Set the color for each family
palF<-c("#dc322f", "#497dcb", "#85995e")
# Attribute a color per image
colF <- palF[traits$family]

#Represent the outline of all the images
panel(outlines, col=colF)

# 2. Traditional Morphometric -------------------
#remove the species and family information
traitsonly <- traits[, 3:15]
names(traitsonly)

# Run a PCA on the traits
pcaTM   <- PCA(traitsonly)
# Visualize the morphospace
plot(pcaTM, col=colF, cex=1.3)

par(mfrow=c(3,1), mar=c(2,5,3,1), las=1)
barplot(sort(pcaTM$rotation[,1]), 
        horiz=TRUE, main="PC1")
barplot(sort(pcaTM$rotation[,2]), 
        horiz=TRUE, main="PC2")
barplot(sort(pcaTM$rotation[,3]), 
        horiz=TRUE, main="PC3")

# 3. Landmarks Analysis -------------------------
# Here we plot the landmarks of Scarus chameleon
# 47th image of our dataset
ldk_plot(landmarks[47], pch=16, cex=1.2, centroid=FALSE)
ldk_labels(landmarks[47], cex = 1.1)

ldkPro <- fgProcrustes(landmarks)

par(mfrow=c(1,2))
stack(landmarks, title="before")
stack(ldkPro, title="after")

# Compute PCA for the landmarks coordinates  
pcaLA <- PCA(ldkPro)

#adjust sign of PC2
pcaLA$x[,2]<- -pcaLA$x[,2]
pcaLA$rotation[,2]<- -pcaLA$rotation[,2]

# Visualize the morphological space defined by PC1 and PC2
plot(pcaLA, col=colF, cex=1.3)

par(mar=c(1,1,0,0), cex=1.5)
PCldkContrib(pcaLA, ldkPro, nax = 1:3, liney=1, linex=1)

# 4. Outline Analysis ---------------------------
# Here we plot the outline of Scarus chameleon
# 47th image of our dataset
coo_plot(outlines[47], xy.axis=FALSE)

# Add the landmarks position
points(outlines[47][outlines$ldk[[47]],], pch="+", 
       col="red", cex=2)
# and their label
text(outlines[47][outlines$ldk[[47]],], labels = 1:5, 
     col="red", pos=c(1,1,4,3,3), cex=1.5)

outPro <- fgProcrustes(outlines)

par(mfrow=c(1,2))
stack(outlines, title="before")
stack(outPro, title="after")

#Determine the number of harmonics
calibrate_harmonicpower_efourier(outPro, nb.h=100)

#Run EFT
fish_efa <- efourier(outPro, norm = FALSE, nb.h = 15)

#Run a PCA
pcaOA <- PCA(fish_efa)

#adjust sign of PC1
pcaOA$x[,1]<- -pcaOA$x[,1]
pcaOA$rotation[,1]<- -pcaOA$rotation[,1]
#adjust sign of PC3
pcaOA$x[,3]<- -pcaOA$x[,3]
pcaOA$rotation[,3]<- -pcaOA$rotation[,3]

#Visualize the morphospace
plot(pcaOA, col=colF, cex = 1.5, 
     xlim=range(pcaOA$x[,1]),
     ylim=range(pcaOA$x[,2]))

PCcontrib(pcaOA, nax = 1:3,gap = 0.8)

# 5. Morphospace comparison ---------------------
res <- cbind(pcaTM$x[, 1:3], 
             pcaLA$x[, 1:3], 
             pcaOA$x[, 1:3])
colnames(res)<- paste(rep(c("TM", "LA", "OA"), each=3), 
                      rep(c("PC1", "PC2", "PC3"), 3), sep="\n")

# Correlation Matrix
corrplot::corrplot.mixed(cor(res), upper="number", lower="ellipse")

# Differentiation of taxonomic groups
par(mfrow = c(3,2))
#TM
plot(pcaTM, fac=traits$family, eigen = FALSE)
plot(pcaTM, fac=traits$family, eigen = FALSE, yax = 3)

#LA
plot(pcaLA, fac=traits$family, 
     morphospace = FALSE, eigen = FALSE)
plot(pcaLA, fac=traits$family, eigen = FALSE,
     morphospace = FALSE, yax = 3)

#OA
plot(pcaOA, fac=traits$family, 
     morphospace = FALSE, eigen = FALSE)
plot(pcaOA, fac=traits$family, eigen = FALSE,
     morphospace = FALSE, yax = 3)

# Quantification of differentiation
#Transform family into categorical numbers
fam <- as.numeric(traits$family)

#Compute the silhouette for each morphospace
sTM <- cluster::silhouette(fam, dist(pcaTM$x[, 1:3]))
sLA <- cluster::silhouette(fam, dist(pcaLA$x[, 1:3]))
sOA <- cluster::silhouette(fam, dist(pcaOA$x[, 1:3]))

#Average silhouette value for TM    
summary(sTM)$avg.width
#Average silhouette value for LA    
summary(sLA)$avg.width
#Average silhouette value for OA    
summary(sOA)$avg.width

