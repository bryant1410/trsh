#!/usr/bin/ruby

require 'json'

def url_get_ruta origen, destino
	"https://maps.googleapis.com/maps/api/directions/json?origin=#{origen}&destination=#{destino}"
end

i = 0
puntos = []

ARGF.each do |line|
    puntos[i] = line.strip
    #puts "#{i} => #{puntos[i]}"
    i += 1
end

for i in 0..(puntos.length-1)
	for j in 0..(puntos.length-1)
		if i != j
			puts url_get_ruta(puntos[i], puntos[j])
		end
	end
end
