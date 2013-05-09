require 'rubygems'
require 'bio'
require 'stringio'
require 'fileutils'
require 'Tempfile'

$stderr.sync = true

HLA_A_LENGTH=787
MIN_HLA_BC_LENGTH=787
MAX_HLA_BC_LENGTH=796
EXON2_LENGTH=270
EXON3_LENGTH=276

# Read sequence type and segments from stdin.
letter = ARGV[0]
fasta_file = ARGV[1]
outfile = Tempfile.new(csv, '.')

# Validate arguments.
if ! ['A', 'B', 'C'].include? letter
	abort("Usage: ruby hla.rb [A|B|C] [fasta_file]")
end

AMBIG = { 
'A' => ['A'],
'T' => ['T'],
'C' => ['C'],
'G' => ['G'],
'R' => ['A', 'G'], 
'Y' => ['C', 'T'], 
'K' => ['G', 'T'], 
'M' => ['A', 'C'], 
'S' => ['C', 'G'],
'W' => ['A', 'T']
}

REV_AMBIG = AMBIG.invert

def all_strings(l, k)
  if k == 1
    return l
  else
    prev = all_strings(l, k-1)
    subsets = []
    l.each do |i| 
      prev.each do |p| 
        subsets << i + p 
      end 
    end 
    return subsets
  end 
end

TUPLE_LEN = 3

b570101 = 'GCTCCCACTCCATGAGGTATTTCTACACCGCCATGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACCCAGTTCGTGAGGTTCGACAGCGACGCCGCGAGTCCGAGGATGGCGCCCCGGGCGCCATGGATAGAGCAGGAGGGGCCGGAGTATTGGGACGGGGAGACACGGAACATGAAGGCCTCCGCGCAGACTTACCGAGAGAACCTGCGGATCGCGCTCCGCTACTACAACCAGAGCGAGGCCGGGTCTCACATCATCCAGGTGATGTATGGCTGCGACGTGGGGCCGGACGGGCGCCTCCTCCGCGGGCATGACCAGTCCGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGAGCTCCTGGACCGCGGCGGACACGGCGGCTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTGTGGCGGAGCAGCTGAGAGCCTACCTGGAGGGCCTGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGG'
b570102 = 'GCTCCCACTCCATGAGGTATTTCTACACCGCCATGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACCCAGTTCGTGAGGTTCGACAGCGACGCCGCGAGTCCGAGGATGGCGCCCCGGGCGCCATGGATAGAGCAGGAGGGGCCGGAGTATTGGGACGGGGAGACACGGAACATGAAGGCCTCCGCGCAGACTTACCGAGAGAACCTGCGGATCGCGCTCCGCTACTACAACCAGAGCGAGGCCGGGTCTCACATCATCCAGGTGATGTATGGCTGCGACGTGGGGCCGGACGGGCGCCTCCTCCGCGGGCATGACCAGTCTGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGAGCTCCTGGACCGCGGCGGACACGGCGGCTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTGTGGCGGAGCAGCTGAGAGCCTACCTGGAGGGCCTGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGG'
b570103 = 'GCTCCCACTCCATGAGGTATTTCTACACCGCCATGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACCCAGTTCGTGAGGTTCGACAGCGACGCCGCGAGTCCGAGGATGGCGCCCCGGGCGCCATGGATAGAGCAGGAGGGGCCGGAGTATTGGGACGGGGAGACACGGAACATGAAGGCCTCCGCGCAGACTTACCGAGAGAACCTGCGGATCGCGCTCCGCTACTACAACCAGAGCGAGGCCGGGTCTCACATCATCCAGGTGATGTATGGCTGCGACGTGGGGCCGGACGGGCGCCTCCTCCGCGGGCATGACCAGTCCGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGAGCTCCTGGACCGCGGCGGACACGGCGGCTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTGTGGCGGAGCAGCTGAGAGCCTACCTGGAGGGCCTGTGTGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGG'

#Gives the number of mismatches, where ambiguous alleles match their non-ambiguous counterparts.
def std_match(std, seq)
		mismatches = 0
		0.upto(std.size - 1) do |i|
				if(!AMBIG[seq[i,1]].include?(std[i,1]))
						mismatches += 1 
				end
		end
		return mismatches
end

# Check if a sequence is the right length and has no invalid characters.
def validate_seq(letter, seq, name)
	if name.end_with?("short")
		return true
	end

	if letter == "A"
		if seq.length != HLA_A_LENGTH
			err = "Sequence %s is the wrong length (%d bp)." % [name, seq.length]
			STDERR.puts err
			return false
		end
	else
		if name.end_with?("exon2") 
			if seq.length != EXON2_LENGTH
				err = "Sequence %s is the wrong length (%d bp)." % [name, seq.length]
				STDERR.puts err
				return false
			end
		elsif name.end_with?("exon3") 
			if seq.length != EXON3_LENGTH
				err = "Sequence %s is the wrong length (%d bp)." % [name, seq.length]
				STDERR.puts err
				return false
			end
		else 
			if seq.length < MIN_HLA_BC_LENGTH or seq.length > MAX_HLA_BC_LENGTH
				err = "Sequence %s is the wrong length (%d bp)." % [name, seq.length]
				STDERR.puts err
				return false
			end
		end
	end
	
	if(!(seq =~ /^[atgcrykmswnbdhv]+$/i))
		STDERR.puts "Sequence %d has invalid characters." % name
		return false
	end
	
	return true
end

# Find all standards which have less than 5 mismatches to the query sequence.
def get_matching_stds(seq, hla_stds)
	matching_stds = []
	hla_stds.each do |std|
		allele, std_seq = std
		if (std_match(std_seq, seq) < 5) 
			matching_stds.push([allele, std_seq])
		end
	end
	return matching_stds.sort_by{|s| s[1]}
end

# Form all pairwise combinations of matching standards.
def combine_stds(matching_stds, ambig_tuples)
	alleles_hash = Hash.new

	# Go through each pair of standards, building the ambiguous version of every
	# pair. 
	matching_stds.each_with_index do |std_a, i_a|
		prev_b="N"*std_a[1].size
		prev_combo="N"*std_a[1].size
		matching_stds.each_with_index do |std_b, i_b|
			
			# Only look at each pair once.
			if i_a < i_b:
				break
			end
			std = ''
	
			0.upto((std_a[1].size / TUPLE_LEN) - 1) do |i| 
				tuple_b = std_b[1][i*TUPLE_LEN, TUPLE_LEN]
				if tuple_b == prev_b[i*TUPLE_LEN, TUPLE_LEN]
					std += prev_combo[i*TUPLE_LEN, TUPLE_LEN]
				else
					tuple_a = std_a[1][i*TUPLE_LEN, TUPLE_LEN]
					if tuple_a == tuple_b
						std += tuple_b
					else
						pair = [tuple_a, tuple_b].sort
						std += ambig_tuples[i][pair]
					end
				end
			end 

			rest = -(std_a[1].size.modulo(TUPLE_LEN))
			-rest.upto(-1) do |i|
				pair = [std_a[1][i,1], std_b[1][i,1]].sort
				if pair[0] == pair[1]
					std += pair[0]
				else
					std += REV_AMBIG[pair]
				end
			end
	
			alleles_hash[std] = [] if(alleles_hash[std] == nil)
			alleles_hash[std].push("#{std_b[0]} - #{std_a[0]}")

			prev_combo = std
			prev_b = std_b[1]
		end
	end

	hla_consensus = []
	alleles_hash.each do |std, allele_list|
		hla_consensus.push([std, allele_list])	
	end

	return hla_consensus
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
		hla_stds.push([tmp[0], tmp[1]+tmp[2]])
	end
end

ambig_tuples = []
0.upto((hla_stds[0][1].size / TUPLE_LEN)-1) do |i|
	ambig_tuples[i] = Hash.new
	done = []
	hla_stds.each do |std|
		tuple = std[1][i*TUPLE_LEN,TUPLE_LEN]
		if not done.include?(tuple)
			done.each do |other|
				pair = [tuple, other].sort
				ambig_tuples[i][pair] = ""
				0.upto(TUPLE_LEN-1) do |j|
					if tuple[j,1] == other[j,1]:
						ambig_tuples[i][pair] += tuple[j,1]
					else
						ambig_tuples[i][pair] += REV_AMBIG[[tuple[j,1], other[j,1]].sort]
					end
				end
			end
			done << tuple
		end
	end
end

unmatched = [[], []]
fasta = Bio::FlatFile.auto(fasta_file)
fasta.each do |entry|
	samp = entry.definition

	puts samp
	next if !validate_seq(letter, entry.seq, samp)

	# Check if the sequence is an exon2 or exon3. If so, try to match it with an
	# existing other exon.
	seq = nil
	is_exon = false
	2.upto(3) do |exon|
		if (samp.end_with?("exon%d" % exon))
			is_exon = true
			samp = samp.split("_")[0]
			unmatched[3-exon].each do |other|
				if (other.definition.start_with?(samp))
					if (exon == 2)
						fullseq = entry.seq + "," + other.seq
						seq = entry.seq + other.seq
					else
						fullseq = entry.seq + "," + other.seq
						seq = other.seq + entry.seq
					end
					unmatched[3-exon].delete(other)
					break
				end
			end
			if (seq == nil)
				unmatched[exon % 2] << entry
				break
			end
		end
	end

	# If it was an exon2 or 3 but didn't have a pair, keep going.
	next if (is_exon and seq == nil)

	# It wasn't an exon2 or 3.
	if seq == nil
		fullseq = seq
		seq = entry.seq
		samp = entry.definition
	end

	if seq.length > EXON2_LENGTH+EXON3_LENGTH
		seq = seq[0 .. 270] + seq[512 .. 787]
	end

	# Find all standards which roughly match the query sequence.
	matching_stds = get_matching_stds(seq, hla_stds)
	puts matching_stds.size
		
	# Now, combine all the stds (pick up that can citizen!)
	puts "Combining standards."
	t = Time.now
	hla_consensus = combine_stds(matching_stds, ambig_tuples)
	puts hla_consensus.size
	puts "Took %d seconds." % (Time.now - t)
	
	puts "Resolving ambiguities."
	min = [9999, []]
	hla_consensus.each do |cons|
		if(cons[0] == seq) #Look for exact matches first
			min[0] = 0
			min[1].push(cons)
		end
	end
	
	if(min[0] != 0) #Otherwise....
	
		# This time we don't match to ambiguities? 
		hla_consensus.each do |cons|
			mismatches = 0
			0.upto(cons[0].size - 1) do |dex|
				if(cons[0][dex, 1] != seq[dex, 1])
					mismatches += 1
				end
				break if(mismatches > min[0])
			end
			
			if(mismatches == min[0])
				min[1].push(cons)
			elsif(mismatches < min[0])
				min[0] = mismatches
				min[1] = [cons]
			end
		end
	end

	puts min[0]
	puts min[1]
		
	#Get list of mismatches
	mislist = []
		
	min[1].each do |cons|	#really, we should just pick one, right?
		0.upto(cons[0].size - 1) do |dex|
			if(cons[0][dex,1] != seq[dex,1])
				mislist.push("#{((dex + ((letter == 'A' and dex > 270) ? 241 : 0)) + 1).to_s}:#{seq[dex]}")
				#mislist.push("#{dex}:#{seq[dex,1]}")
			end
		end
	end
	
	mislist.uniq!
	mismatch_count = min[0]
	mismatches = mislist.join(";")
	best_matches = min[1]
	
	#Clean the alleles
	
	puts "Cleaning alleles"
	
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
					
	#Need a new way to resolve ambiguities right here.	Must take into account BOTH
	#alleles to find the most common 
	
	# Strip leading "A:" , "B:", or "C:" from each allele.
	collection = alleles.map do |a|
		tmp = a.split('-')
		[tmp[0].gsub(/[^\d:]/, '').split(':'), tmp[1].gsub(/[^\d:]/, '').split(':')]
	end
	
	puts "Resolving ambiguity"
	
	if(collection.map{|e| [e[0][0], e[1][0]]}.uniq.size != 1)
		ambig = '1'
		collection_ambig = collection.map{|e| [e[0][0 .. 1], e[1][0 .. 1]]}.uniq 

		collection_ambig.each do |a|
			freq = hla_freqs[a]
			if freq == nil:
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
		alleles.delete_if {|a| !(a =~ /^#{letter}\*#{a1}:([^\s])+\s-\s#{letter}\*#{a2}:/)} 
	end
	
	#non ambiguous now, do the easy way
	
	collection = alleles.map do |a|
		tmp = a.split('-')
		tmp[0].gsub!(/[^\d:]/, '') 
		tmp[1].gsub!(/[^\d:]/, '') 
		[tmp[0].strip.split(':'), tmp[1].strip.split(':')]
	end

	4.downto(1) do |i|
		if (collection.map{|a| a[0][0, i]}.uniq.size == 1)
			clean_allele = letter + "\*" + collection[0][0][0,i].join(':') + " - "
			break
		end
	end

	4.downto(1) do |i|
		if (collection.map{|a| a[1][0, i]}.uniq.size == 1)
			clean_allele += letter + "\*" + collection[0][1][0,i].join(':')
			break
		end
	end

	#OK, now we must find homozygousity.	IE: Cw*0722 - Cw*0722
	puts "Finding homozygousity"
	homozygous = '0'
	# Lets say if we detect two of the same mixtures, its heterozygous
	best_matches.each do |a|
		a[1].each do |allele|
			tmp = allele.split(' - ')
			if(tmp[0] == tmp[1])
				homozygous = '1'
			end
		end
	end
			
	# For some reason it's important to know if the allele is B*57:01.
	if letter == 'B'
		#Find B*5701
		puts "Finding B*5701's"
		dist_from_b5701 = [0,0,0]
		[b570101, b570102, b570103].each_with_index do |cons_b5701, dex|
			0.upto(cons_b5701.size - 1) do |i|
				# If either position is ambiguous, but the ambiguities overlap completely
				# (one is a subset of the other), don't increment the distance.
				if(!(AMBIG[cons_b5701[i,1]].all? {|a| AMBIG[seq[i,1]].include?(a)} or
					   AMBIG[seq[i,1]].all? {|a| AMBIG[cons_b5701[i,1]].include?(a)}))
					dist_from_b5701[dex] += 1
				end
			end
		end
			
		dist_from_b5701 = dist_from_b5701.min
				
		b5701 = '0'
		best_matches.each do |a|
			a[1].each do |allele|
				if(allele.include?('B*57:01'))
					b5701 = '1'
				end
			end
		end
	end
	
	alleles_all = best_matches.map {|a| a[1].join(';')}.join(';')
	if alleles_all.size > 3900  
		alleles_all = alleles_all[0 .. 3920].gsub(/;[^;]+$/, ';...TRUNCATED')
	end

	if letter == "B"
		row = [samp, clean_allele, alleles_all, ambig, homozygous, mismatch_count,
			mismatches, fullseq]
	else
#		row = [samp, clean_allele, alleles_all, ambig, homozygous, mismatch_count,
#			mismatches, fullseq]
		puts clean_allele
		puts alleles_all
	end
#	puts row.join(",")
end
