#! /usr/bin/env Rscript
d<-scan("stdin", quiet=TRUE)
cat(min(d), max(d), median(d), mean(d), sd(d), sep="\n")
png("fitness_histograma.png")
h<-hist(d,
	breaks = 30,
	main = "Histograma del mejor fitness de cada ejecución para CH04",
	xlab = "Mejor fitness",
	ylab = "Frecuencia de ejecuciones")

xfit<-seq(min(d)-3000,max(d)+3000,length=40) 
yfit<-dnorm(xfit,mean=mean(d),sd=sd(d)) 
yfit <- yfit*diff(h$mids[1:2])*length(d) 
lines(xfit, yfit, col="blue", lwd=2)
