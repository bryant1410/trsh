#!/usr/bin/ruby

require 'csv'
require 'ostruct'
require 'pqueue'

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

class Asignacion
	attr_accessor :i_manzana
	attr_accessor :j_contenedor
	attr_accessor :valor
end

def demanda contenedores, manzanas
	asignacion_m = Array.new(manzanas.length) { |i| [] }
	asignacion_c = Array.new(contenedores.length) { |i| [] }

	contenedores.each_index { |j|
		elegidos = PQueue.new() { |a, b| a.valor > b.valor }

		manzanas.each_index { |i|
			if manzanas[i].poblacion > 0
				asignacion = Asignacion.new
				asignacion.i_manzana = i
				asignacion.j_contenedor = j
				asignacion.valor = distancia [contenedores[j].lat, contenedores[j].long], [manzanas[i].lat, manzanas[i].long]

				if elegidos.size == $cant_minimos
					if elegidos.top.valor > asignacion.valor
						elegidos.pop
						elegidos.push asignacion
					end
				else
					elegidos.push asignacion
				end
			end
		}

		elegidos.each_pop do |e|
			e.valor = 1/(e.valor*e.valor)
			asignacion_m[e.i_manzana] << e
			asignacion_c[e.j_contenedor] << e
		end
	}

	manzanas.each_index { |i|
		suma_total = suma(asignacion_m[i].map { |e| e.valor })

		if suma_total > 0
			asignacion_m[i].map! { |e|
				e.valor = manzanas[i].poblacion*e.valor/suma_total
				e
			}
		end
	}

	asignacion_c.map { |es| suma es.map { |e| e.valor } }
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

contenedores = []

cabezal = false

CSV.foreach(ARGV[1]) do |fila|
	if cabezal
		cabezal = false
	else
		contenedor = OpenStruct.new
		contenedor.lat = fila[0].to_f
		contenedor.long = fila[1].to_f

		contenedores << contenedor
	end
end

puts demanda(contenedores, manzanas)
