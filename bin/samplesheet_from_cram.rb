this_dir = Dir.pwd

root_directory = ARGV.shift

folders = Dir["#{root_directory}/Indi*"]

puts "Individual;Sample;Cram;Bai"

folders.each do |folder|

	sample_folders = Dir[folder + "/Sample*"]

	sample_folders.each do |sample_folder|

		cram_file = Dir[sample_folder + "/*.cram"].shift

		individual = folder.split("/")[-1]
		sample = sample_folder.split("/")[-1]

		puts "#{individual};#{sample};#{this_dir}/#{cram_file};#{this_dir}/#{cram_file}.bai"
	end
end
