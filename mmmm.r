#! /usr/bin/env Rscript
d<-scan("stdin", quiet=TRUE)
cat(min(d), max(d), median(d), mean(d), sd(d), sep="\n")
png("demanda_histograma.png")
hist(d,
	breaks = 100,
	main = "Histograma de la demanda de personas de los contenedores",
	xlab = "Cantidad de personas",
	ylab = "Frecuencia")
