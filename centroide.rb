#!/usr/bin/ruby

require 'csv'
require 'ostruct'

def centroide v
	n = v.length

	v << v[0] # El último lo pongo al final así lo puedo acceder como el n-ésimo.

	a = 0
	for i in 0..(n-1)
		a += v[i].x*v[i+1].y - v[i+1].x*v[i].y
	end
	a *= 0.5

	c = OpenStruct.new

	c.x = 0
	for i in 0..(n-1)
		c.x += (v[i].x + v[i+1].x)*(v[i].x*v[i+1].y - v[i+1].x*v[i].y)
	end
	c.x /= 6*a

	c.y = 0
	for i in 0..(n-1)
		c.y += (v[i].y + v[i+1].y)*(v[i].x*v[i+1].y - v[i+1].x*v[i].y)
	end
	c.y /= 6*a

	c
end

i = 0
shapes = []

shapes[0] = []

CSV.foreach(ARGV[0]) do |row|

	unless row[0].eql? "shapeid" # ignoro cabezal
		shapeid = row[0].to_i
		x = row[1].to_f
		y = row[2].to_f

		if shapeid > i # cambio de figura
			c = centroide(shapes[i])
			puts "#{c.x},#{c.y}"

			i += 1
			shapes[i] = []
		end

		vertice = OpenStruct.new
		vertice.x = x
		vertice.y = y

		shapes[i] << vertice
	end
end

c = centroide(shapes[i])
puts "#{c.x},#{c.y}"
