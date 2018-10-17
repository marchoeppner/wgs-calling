#!/bin/env ruby

IO.readlines(ARGV.shift).each do |line|
	next if line.match(/^@.*/)
	elements = line.strip.split("\t")
	puts elements[0] + ":" + elements[1] + "-" + elements[2]
end
