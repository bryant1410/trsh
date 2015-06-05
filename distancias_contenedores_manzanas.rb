#!/usr/bin/ruby

require 'csv'
require 'ostruct'

$rad_per_deg = Math::PI/180  # PI / 180
$rkm = 6371                  # Earth radius in kilometers
$rm = $rkm * 1000            # Radius in meters

$cant_minimos = 4

# https://stackoverflow.com/questions/12966638/rails-how-to-calculate-the-distance-of-two-gps-coordinates-without-having-to-u
def distancia a, b
  dlon_rad = (b[1]-a[1]) * $rad_per_deg  # Delta, converted to rad
  dlat_rad = (b[0]-a[0]) * $rad_per_deg

  lat1_rad, _lon1_rad = a.map! {|i| i * $rad_per_deg }
  lat2_rad, _lon2_rad = b.map! {|i| i * $rad_per_deg }

  a = Math.sin(dlat_rad/2)**2 + Math.cos(lat1_rad) * Math.cos(lat2_rad) * Math.sin(dlon_rad/2)**2
  c = 2 * Math.asin(Math.sqrt(a))

  $rm * c # Delta in meters
end

def distancias_a contenedor, manzanas
	Array.new(manzanas.length) { |i| distancia([contenedor.lat, contenedor.long], [manzanas[i].lat, manzanas[i].long]) }
end

def index_of_min array
	array.each_with_index.min[1]
end

def suma array
	array.inject 0, :+
end

manzanas = []

cabezal = true

CSV.foreach(ARGV[0]) do |fila|
	if cabezal
		cabezal = false
	else
		manzana = OpenStruct.new
		manzana.poblacion = fila[16].to_i
		manzana.lat = fila[25].to_f
		manzana.long = fila[26].to_f

		manzanas << manzana
	end
end

#contenedores = []

cabezal = false

i = 0

CSV.foreach(ARGV[1]) do |fila|
	if cabezal
		cabezal = false
	else
		#contenedor = OpenStruct.new
		#contenedor.lat = fila[0].to_f
		#contenedor.long = fila[1].to_f

		lat = fila[0].to_f
		long = fila[1].to_f

		manzanas.each_index { |j|
			distancia = distancia([lat, long], [manzanas[j].lat, manzanas[j].long])
			puts "#{i} #{j} #{distancia}"
		}

		#contenedores << contenedor
		i += 1
	end
end
