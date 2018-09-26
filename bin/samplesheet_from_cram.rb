directory = ARGV.shift.to_s
files = Dir[directory + "/*.cram"]

puts "IndividualID;sampleID;Cram;Crai"

files.each do |file|
	path = File.realpath(file)

	individual = path.split("/")[-3]
	sample = path.split("/")[-2]

	puts "#{individual};#{sample};#{path};#{path}.bai"	

end
