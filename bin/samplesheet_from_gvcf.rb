directory = ARGV.shift.to_s
files = Dir[directory + "/*.g.vcf.gz"]

puts "IndividualID;sampleID;vcf;tbi"

files.each do |file|
	path = File.realpath(file)

	individual = path.split("/")[-3]
	sample = path.split("/")[-2]

	puts "#{individual};#{sample};#{path};#{path}.tbi"	

end
