#' ---
#' title: "Class5: Data exploration and visualization in R "
#' author: "Barry Grant"
#' output: github_document
#' ---



# Class5 Data visualization
x <- rnorm(1000)

# some summary stats
mean(x)
sd(x)

summary(x)
boxplot(x)

hist(x)

# Section 2 scaterplots
# lets read our input file first

baby <- read.table("bimm143_05_rstats/weight_chart.txt",
                   header = TRUE)

plot(baby$Age, baby$Weight,typ="o", cex=1.5, pch=15, 
     col="blue", lwd=2, xlab="Baby age (months)",
     ylab="Baby weight (kg)")

##
#read.table(", 
feat <- read.table("bimm143_05_rstats/feature_counts.txt", 
           sep="\t", header=TRUE)

mouse <- read.delim("bimm143_05_rstats/feature_counts.txt")

par(mar=c(5,11,2,2) )
barplot(mouse$Count, names.arg = mouse$Feature, 
        horiz = TRUE, las=1)




