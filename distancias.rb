#!/usr/bin/ruby

require 'json'
require 'optparse'
require 'ostruct'

def url_get_ruta origen, destino
	"https://maps.googleapis.com/maps/api/directions/json?origin=#{origen}&destination=#{destino}"
end

def get_ruta origen, destino, proxy
	JSON.parse `curl -sS #{proxy} "#{url_get_ruta(origen, destino)}"`
end

def print_duracion puntos, i, j, options
	ruta = get_ruta puntos[i], puntos[j], options.proxy
	if ruta.nil?
		$stderr.puts "Respuesta vacía."
		exit 1
	else
		if ruta['status'] == 'OK'
			duracion = ruta['routes'][0]['legs'][0]['duration']['value']
			puts "#{i} #{j} #{duracion}"
		else
			$stderr.puts ruta
			exit 1
		end
	end
end

options = OpenStruct.new
options.proxy = ''
options.i_0 = 0
options.j_0 = 0
options.origin_destination = false

OptionParser.new do |opts|
	opts.banner = "Usage: distancias"
	opts.separator  ""
	opts.separator  "Options"

	opts.on('-u', '--use-proxy', 'Use default proxy') do
		#options.proxy = '--proxy "httpproxy.fing.edu.uy:3128"'
		#options.proxy = ''
		options.proxy = '--proxy localhost:3128'
	end

	opts.on('-p', '--proxy PROXY', 'Use proxy') do |proxy|
		options.proxy = "--proxy #{proxy}"
	end

	opts.on('-c', '--continue I,J', Array, 'Continue from certain point') do |indices|
		options.i_0 = indices[0].to_i
		options.j_0 = indices[1].to_i
	end

	opts.on('-o', '--origin-destination', 'Origin-Destination mode') do
		options.origin_destination = true
	end
end.parse!

i = 0
puntos = []

ARGF.each do |line|
	puntos[i] = line.strip
    #puts "#{i} => #{puntos[i]}"
	i += 1
end

if options.origin_destination
	for i in 2..(puntos.length-1)
		print_duracion puntos, 0, i, options
	end
	for i in 2..(puntos.length-1)
		print_duracion puntos, i, 1, options
	end
else
	for i in options.i_0..(puntos.length-1)
		for j in options.j_0..(puntos.length-1)
			if i != j
				print_duracion puntos, i, j, options
			end
		end
		options.j_0 = 0
	end
end
