require 'rubygems'
require 'bio'
require 'stringio'
require 'fileutils'
require 'tempfile'
require 'optparse'
require 'date'

$stderr.sync = true
$stdout.sync = true

HLA_A_LENGTH=787
MIN_HLA_BC_LENGTH=787
MAX_HLA_BC_LENGTH=796
EXON2_LENGTH=270
EXON3_LENGTH=276

# Read sequence type and segments from stdin.
options = {}
OptionParser.new do |opts|
	opts.banner = "Usage: ruby hla-easy.rb [options] [A|B|C]"
	opts.on("-t", "--threshold [INTEGER]", Integer, 
		"Maximum mismatches to accept") do |t|
		options[:threshold] = t;
	end
end.parse!

letter = ARGV[0]
threshold = options[:threshold] 
if threshold < 0
	threshold = nil
end
$log = IO.new(7, "w")
$details = IO.new(8, "w")

# Validate arguments.
if ! ['A', 'B', 'C'].include? letter
	$stderr.puts("Bad locus: expected A, B, or C, got #{letter}")
	exit(1)
end

AMBIG = { 
'a' => ['a'],
't' => ['t'],
'c' => ['c'],
'g' => ['g'],
'r' => ['a', 'g'], 
'y' => ['c', 't'], 
'k' => ['g', 't'], 
'm' => ['a', 'c'], 
's' => ['c', 'g'],
'w' => ['a', 't']
}
REV_AMBIG = AMBIG.invert

NUC2BIN = {
	'a' => 0b0001,
	'c' => 0b0010,
	'g' => 0b0100,
	't' => 0b1000,
	'm' => 0b0011,
	'r' => 0b0101,
	'w' => 0b1001,
	's' => 0b0110,
	'y' => 0b1010,
	'k' => 0b1100,
	'v' => 0b1110,
	'h' => 0b1101,
	'd' => 0b1011,
	'b' => 0b0111,
	'n' => 0b1111
}
BIN2NUC = NUC2BIN.invert

# Check if a sequence is the right length and has no invalid characters.
#def validate_seq(letter, seq, name)
def check_length(letter, seq, name)
	err = false
	if name.downcase.end_with?("short")
		if letter == "A" 
			err = (seq.length >= HLA_A_LENGTH)
		elsif name.include? "exon2"
			err = (seq.length >= EXON2_LENGTH)
		elsif name.include? "exon3" 
			err = (seq.length >= EXON3_LENGTH)
		else 
			err = (seq.length >= MAX_HLA_BC_LENGTH)
		end

	elsif letter == "A"
		err = (seq.length != HLA_A_LENGTH)
	elsif name.include? "exon2"
		err = (seq.length != EXON2_LENGTH)
	elsif name.include? "exon3"
		err = (seq.length != EXON3_LENGTH)
	else
		err = ((seq.length > MAX_HLA_BC_LENGTH) or 
		       (seq.length < MIN_HLA_BC_LENGTH))
	end

	if err
		err = "Sequence %s is the wrong length (%d bp). Check the locus %s." % [name, seq.length, letter]
		$log.puts err
		return false
	end
	return true
end

def check_bases(seq, name)
	if(!(seq =~ /^[atgcrykmswnbdhv]+$/i))
		$log.puts "Sequence %s has invalid characters." % name
		return false
	end
	return true
end

def calc_padding(std, seq)
	best = 10E10
	pad = std.length - seq.length
	left_pad = 0
	0.upto(pad) do |i|
		pseq = nuc2bin("n"*i) + seq + nuc2bin("n"*(pad-i))
		mismatches = std_match(std, pseq)
		if mismatches < best
			best = mismatches
			left_pad = i
		end
	end
	return [left_pad, pad-left_pad]
end

def pad_short(letter, seq, name, hla_stds)
	if name.include? "exon2"
		std = hla_stds[0][1][0,EXON2_LENGTH]
	elsif name.include? "exon3"
		std = hla_stds[0][1][EXON2_LENGTH,EXON3_LENGTH]
	else
		has_intron = true
		std = hla_stds[0][1]
	end

	if has_intron
		left_pad, _ = calc_padding(std[0,EXON2_LENGTH], seq[0,EXON2_LENGTH/2])
		_, right_pad = calc_padding(std[-EXON3_LENGTH,EXON3_LENGTH], 
			seq[-EXON3_LENGTH/2,EXON3_LENGTH/2])
	else
		left_pad, right_pad = calc_padding(std, seq)
	end
	
	return nuc2bin("n"*left_pad) + seq + nuc2bin("n"*right_pad)
end


#Gives the number of mismatches, where ambiguous alleles match their non-ambiguous counterparts.
def std_match(std, seq)
	mismatches = 0
	delta = 0
	0.upto(std.size - 1) do |i|
		if seq[i] & std[i] == 0
			mismatches += 1 
		end
	end
	return mismatches
end

def nuc2bin(seq)
	binseq = []
	0.upto(seq.size - 1) do |i|
		binseq << NUC2BIN[seq[i,1]]
	end
	return binseq
end

def bin2nuc(seq)
	nseq = ""
	0.upto(seq.length - 1) do |i|
		nseq = nseq + BIN2NUC[seq[i]]
	end
	return nseq
end

# Find all standards which have less than 5 mismatches to the query sequence.
def get_matching_stds(seq, hla_stds)
	matching_stds = []
	hla_stds.each do |std|
		allele, std_seq = std
		mismatches = std_match(std_seq, seq)
		if mismatches < 5
			matching_stds.push([allele, std_seq, mismatches])
		end
	end
	return matching_stds.sort_by{|s| s[2]}
end

# Form all pairwise combinations of matching standards.
def combine_stds(matching_stds, seq, threshold)
	alleles_hash = Hash.new
	len = matching_stds[0][1].size

	min = 9999
	if threshold == nil
		tmp_threshold = 0
	else
		tmp_threshold = threshold
	end
	combos = Hash.new

	# Go through each pair of standards, building the ambiguous version of every
	# pair. 
	matching_stds.each_with_index do |std_a, i_a|
		next if std_a[2] > [min, tmp_threshold].max
		matching_stds.each_with_index do |std_b, i_b|
			# Only look at each pair once.
			if i_a < i_b
				break
			end
			next if std_b[2] > [min, tmp_threshold].max

			std = []

			mismatches = 0
			0.upto(len-1) do |i|
				if std_b[1][i] == std_a[1][i]
					std << std_b[1][i]
				else
					std << (std_b[1][i] | std_a[1][i])
				end
				
				if (std[i] ^ seq[i]) & 15 != 0
					mismatches += 1
				end
				break if(mismatches > [min, tmp_threshold].max)
			end

			if(mismatches <= [min, tmp_threshold].max)
				if mismatches < min
					min = mismatches
				end
				combos[mismatches] = Hash.new if combos[mismatches] == nil
				combos[mismatches][std] = [] if combos[mismatches][std] == nil
				combos[mismatches][std] << [std_a[0], std_b[0]].sort
			end
		end
	end
	
	result = []
	combos.each do |c|
		cur_combos = []
		c[1].each do |std, allele_list|
			cur_combos << [std, allele_list]
		end
		result << [c[0], cur_combos]
	end

	result = result.sort_by{|c| c[0]}
	return result
end

# Read in HLA frequencies.
column = case letter
	when "A" then 0
	when "B" then 2
	else 4
end
	
hla_freqs = Hash.new
File.open("hla_frequencies.csv") do |file|
	file.each_line do |line|
		tmp = line.split(',').map{|a| a.strip }[column,2]
		tmp = tmp.map {|a| [a[0,2], a[2,2]] }
		hla_freqs[tmp] = 0 if hla_freqs[tmp] == nil
		hla_freqs[tmp] += 1
	end
end

# Read in HLA standards.
hla_stds = []
File.open("hla_#{letter.downcase}_std_reduced.csv") do |file|
	file.each_line do |line|
		tmp = line.strip.split(",")
		allele = tmp[0]
		seq = nuc2bin((tmp[1]+tmp[2]).downcase)
		hla_stds.push([tmp[0], seq])
	end
end

now = DateTime.now.strftime("%+")
mtime = ""
File.open('hla_nuc.fasta.mtime', 'r') do |file|
	mtime = file.readline.strip
end
puts "Run commencing #{now}. Allele definitions last updated #{mtime}."
puts "ENUM,ALLELES_CLEAN,ALLELES,AMBIGUOUS,HOMOZYGOUS,MISMATCH_COUNT,MISMATCHES,EXON2,INTRON,EXON3"

$details.puts "Run commencing #{now}. Allele definitions last updated #{mtime}."
$details.puts "ALLELE,MISMATCHES,EXON2,INTRON,EXON3"

unmatched = [[], []]
nseqs = 0
npats = 0

fasta = Bio::FlatFile.auto($stdin)
fasta.each do |entry|
	samp = entry.definition

	next if !check_length(letter, entry.seq, samp)
	next if !check_bases(entry.seq, samp)

	# Check if the sequence is an exon2 or exon3. If so, try to match it with an
	# existing other exon.
	is_exon = false
	matched = false
	exon2 = ""
	intron = ""
	exon3 = ""
	2.upto(3) do |exon|
		if (samp.include?("exon%d" % exon))
			is_exon = true
			samp = samp.split("_")[0]
			unmatched[3-exon].each do |other|
				if (other.definition.start_with?(samp))
					matched = true
					intron = ""
					if (exon == 2)
						exon2 = entry.naseq
						exon3 = other.naseq
					else
						exon2 = other.naseq
						exon3 = entry.naseq
					end
					unmatched[3-exon].delete(other)
					break
				end
			end
			if not matched
				unmatched[exon % 2] << entry
				break
			end
		end
	end

	# If it was an exon2 or 3 but didn't have a pair, keep going.
	next if is_exon and not matched

	if is_exon
		exon2_bin = pad_short(letter, nuc2bin(exon2.downcase), "exon2", hla_stds)
		exon3_bin = pad_short(letter, nuc2bin(exon3.downcase), "exon3", hla_stds)
		exon2 = bin2nuc(exon2_bin)
		exon3 = bin2nuc(exon3_bin)
		seq = exon2_bin + exon3_bin
	else 
		seq = pad_short(letter, nuc2bin(entry.naseq.downcase), "", hla_stds)
		exon2 = bin2nuc(seq[0,EXON2_LENGTH])
		intron = bin2nuc(seq[EXON2_LENGTH .. -EXON3_LENGTH-1])
		exon3 = bin2nuc(seq[-EXON3_LENGTH, EXON3_LENGTH])
		seq = seq[0,EXON2_LENGTH] + seq[-EXON3_LENGTH,EXON3_LENGTH]
	end

	# Find all standards which roughly match the query sequence.
	matching_stds = get_matching_stds(seq, hla_stds)
	if matching_stds.size == 0
		$log.write("Sequence #{samp} did not match any known alleles. ")
		$log.write("Please check the locus and orientation.\n")
		next
	end
		
	# Now, combine all the stds (pick up that can citizen!)
	all_combos = combine_stds(matching_stds, seq, threshold)

	# Write out individual mismatches for each possibility below the threshold.
	if threshold != nil
		all_combos.each_with_index do |combos, idx|
			if combos[0] > threshold
				if idx == 0
					$log.write("No matches found below specified threshold. ")
					$log.write("Please check the locus, orientation, and/or increase ")
					$log.write("number of mismatches.\n")
				end
				break
			end
	
			combos[1].each do |cons|	#really, we should just pick one, right?
				cons[1].each do |pair|
					$details.write("#{pair.join(" - ")},")
					misstrings = []
					0.upto(cons[0].size - 1) do |i|
						base = BIN2NUC[seq[i]].upcase
						if(cons[0][i] ^ seq[i]) != 0
							correct_base = BIN2NUC[cons[0][i]].upcase
							if letter == 'A' and i > 270
								dex = i+242
							else
								dex = i+1
							end
							misstrings << "#{dex}:#{base}->#{correct_base}"
						end
					end
					$details.write("#{misstrings.join(';')},")
					$details.write("#{exon2.upcase},#{intron.upcase},#{exon3.upcase}\n")
				end
			end
		end
	end

	best_matches = all_combos[0][1]
	mismatch_count = all_combos[0][0]

	mishash = Hash.new
	best_matches.each do |cons|	#really, we should just pick one, right?
		0.upto(cons[0].size - 1) do |i|
			base = BIN2NUC[seq[i]].upcase
			if(cons[0][i] ^ seq[i]) != 0
				correct_base = BIN2NUC[cons[0][i]]
				if letter == 'A' and i > 270
					dex = i+241
				else
					dex = i+1
				end
				mishash[i] = [] if mishash[i] == nil
				if not mishash[i].include? cons[0][i]
					mishash[i] << cons[0][i]
				end
			end
		end
	end

	mislist = []
	mishash.each do |m|
		if letter == 'A' and m[0] > 270
			dex = m[0]+241
		else
			dex = m[0]+1
		end
		base = BIN2NUC[seq[m[0]]].upcase
		correct_bases = ""
		m[1].each do |correct_bin|
			if correct_bases == ""
				correct_bases = BIN2NUC[correct_bin].upcase
			else
				correct_bases += "/" + BIN2NUC[correct_bin].upcase
			end
		end
		mislist.push("#{dex.to_s}:#{base}->#{correct_bases}")
	end
	
	mislist = mislist.sort_by{|b| b.split(":")[0].to_i}
	mismatches = mislist.join(";")

	# Clean the alleles
	
	# Column offset in frequencies file.
	fcnt = case letter
		when 'A' then 0 
		when 'B' then 2
		else 4
	end
	
	clean_allele = ''
	alleles = []
	ambig = '0'

	best_matches.each do |a|
		alleles += a[1]
	end

	# Strip leading "A:" , "B:", or "C:" from each allele.
	collection = alleles.map do |a|
		[a[0].gsub(/[^\d:]/, '').split(':'), a[1].gsub(/[^\d:]/, '').split(':')]
	end
	
	if(collection.map{|e| [e[0][0], e[1][0]]}.uniq.size != 1)
		ambig = '1'
		collection_ambig = collection.map{|e| [e[0][0 .. 1], e[1][0 .. 1]]}.uniq 

		collection_ambig.each do |a|
			freq = hla_freqs[a]
			if freq == nil
				freq = 0
			end
			a.push(freq)
		end
						
		# Try to find the allele occuring the maximum number of times. If it's a tie,
		# just pick the alphabetically first one. 
		max_allele = collection_ambig.max do |a,b| 
			if(a[2] != b[2]) #Go by frequency
				a[2] <=> b[2]
			elsif(b[0][0].to_i != a[0][0].to_i) #Then lowest first allele
				b[0][0].to_i <=> a[0][0].to_i
			elsif(b[0][1].to_i != a[0][1].to_i) 
				b[0][1].to_i <=> a[0][1].to_i
			elsif(b[1][0].to_i != a[1][0].to_i) #Then lowest second allele
				b[1][0].to_i <=> a[1][0].to_i
			else 
				b[1][1].to_i <=> a[1][1].to_i
			end
		end
	
		a1 = max_allele[0][0]
		a2 = max_allele[1][0]
		alleles.delete_if {|a| !(a[0] =~ /^#{letter}\*#{a1}:([^\s])+/)}
		alleles.delete_if {|a| !(a[1] =~ /^#{letter}\*#{a2}:([^\s])+/)}
	end
	
	#non ambiguous now, do the easy way
	
	collection = alleles.map do |a|
		[a[0].strip.split(':'), a[1].strip.split(':')]
	end

	4.downto(1) do |i|
		if (collection.map{|a| a[0][0, i]}.uniq.size == 1)
			clean_allele = collection[0][0][0,i].join(':').gsub(/[A-Z]$/, '') + " - "
			break
		end
	end

	4.downto(1) do |i|
		if (collection.map{|a| a[1][0, i]}.uniq.size == 1)
			clean_allele += collection[0][1][0,i].join(':').gsub(/[A-Z]$/, '')
			break
		end
	end

	#OK, now we must find homozygousity.	IE: Cw*0722 - Cw*0722
	homozygous = '0'
	# Lets say if we detect two of the same mixtures, its heterozygous
	best_matches.each do |a|
		a[1].each do |allele|
			if(allele[0] == allele[1])
				homozygous = '1'
			end
		end
	end
			
	alleles_all = []
	best_matches.each do |_, alleles|
		alleles.each do |a|
			alleles_all << "%s - %s" % a
		end
	end
	alleles_all = alleles_all.sort.join(";")
	if alleles_all.size > 3900  
		alleles_all = alleles_all[0 .. 3920].gsub(/;[^;]+$/, ';...TRUNCATED')
	end

	npats += 1
	if is_exon
		nseqs += 2
	else
		nseqs += 1
	end

	row = [samp, clean_allele, alleles_all, ambig, homozygous,
		mismatch_count, mismatches, exon2.upcase, intron.upcase, exon3.upcase]
	puts row.join(",")
end

unmatched[0].each do |exon2|
	$log.puts "No matching exon 3 found for #{exon2.definition}."
end

unmatched[1].each do |exon3|
	$log.puts "No matching exon 2 found for #{exon3.definition}."
end

puts "Done: #{nseqs} sequences from #{npats} patients successfully processed"
$details.write("Done\n")
$log.close()
$details.close()
