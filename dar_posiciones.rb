#!/usr/bin/ruby

puntos = []

File.open(ARGV[0]).each do |fila|
	puntos << fila
end

STDIN.each_line do |line|
	puts puntos[line.to_i]
end
