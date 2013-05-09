#!/usr/bin/ruby
#
# This program downloads the newest version of hla_nuc.fasta and then 
# processes it, creating three standard files for a, b, and c.  Any sequence
# that doesn't seem to align properly will be rejected.  
# The alignment algorithm I'm using is the dirtiest possible.

require 'net/ftp'

#scores the sequence according to how many characters don't match.  An
#alignment of 0 means a perfect match.  We optimize it by assuming anything under 20 is probably
#a match.
def score(seq, align)
	maxscore = align.size
	maxseq = -1;

	0.upto(seq.size - align.size) do |i|
		score = align.size
		0.upto(align.size - 1) do |j|
			if(align[j] == seq[i + j])
				score -= 1
			end
		end
		
		if(score < maxscore)
			maxscore = score
			maxseq = i
			if(maxscore < 20)
				return [maxscore, seq[maxseq, align.size], maxseq]
			end
		end
	end
	return [maxscore, seq[maxseq, align.size], maxseq]
end


#Exon sequences to compare against(used for scoring)
a_exon2_align='GCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGACG'
a_exon3_align='GTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCCGCGGGTACCGGCAGGACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGATCACCAAGCGCAAGTGGGAGGCGGTCCATGCGGCGGAGCAGCGGAGAGTCTACCTGGAGGGCCGGTGCGTGGACGGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGG'
b_exon2_align='GCTCCCACTCCATGAGGTATTTCTACACCTCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCTCAGTGGGCTACGTGGACGACACCCAGTTCGTGAGGTTCGACAGCGACGCCGCGAGTCCGAGAGAGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCGGAACACACAGATCTACAAGGCCCAGGCACAGACTGACCGAGAGAGCCTGCGGAACCTGCGCGGCTACTACAACCAGAGCGAGGCCG'
b_exon3_align='GGTCTCACACCCTCCAGAGCATGTACGGCTGCGACGTGGGGCCGGACGGGCGCCTCCTCCGCGGGCATGACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCCTGGACCGCCGCGGACACGGCGGCTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTGAGGCGGAGCAGCGGAGAGCCTACCTGGAGGGCGAGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGACAAGCTGGAGCGCGCTG'
c_exon2_align='GCTCCCACTCCATGAAGTATTTCTTCACATCCGTGTCCCGGCCTGGCCGCGGAGAGCCCCGCTTCATCTCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGTCCGAGAGGGGAGCCGCGGGCGCCGTGGGTGGAGCAGGAGGGGCCGGAGTATTGGGACCGGGAGACACAGAAGTACAAGCGCCAGGCACAGACTGACCGAGTGAGCCTGCGGAACCTGCGCGGCTACTACAACCAGAGCGAGGCCG'
c_exon3_align='GGTCTCACACCCTCCAGTGGATGTGTGGCTGCGACCTGGGGCCCGACGGGCGCCTCCTCCGCGGGTATGACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCCTGGACCGCCGCGGACACCGCGGCTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTGAGGCGGAGCAGCGGAGAGCCTACCTGGAGGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCGCGG'

filename = 'hla_nuc.fasta'

#first check to see if the file has been updated.

#read the old modification time
mtime_old = ''
begin
	File.open('hla_nuc.fasta.mtime', 'r') do |file|
		mtime_old = file.readline.strip
	end
rescue 
#Just in case the file doesn't exist yet
end

#log in to the server and check for a newer file
ftp = Net::FTP.new('ftp.ebi.ac.uk')
ftp.login
ftp.chdir('pub/databases/imgt/mhc/hla/')
mtime = ftp.mtime('hla_nuc.fasta').to_s.strip
if(mtime != mtime_old) #download a new file
	puts "New hla_nuc.fasta found!  Downloading!"
	ftp.getbinaryfile('hla_nuc.fasta')
	puts "Finished downloading"
	File.open('hla_nuc.fasta.mtime', 'w') do |file|
		file.puts mtime
	end
end
ftp.close

puts "Parsing " + filename;

hla_a = []
hla_b = []
hla_c = []

diff_reject = 32

fasta = []
enu=[]
File.open(filename) do |file|
  file.each_line do |line|
    if(line =~ /^>/)
      fasta.push(enu)
      enu = [line.strip, '']
    else
      enu[1] += line.strip
    end
  end
  fasta.push(enu)
end

fasta.delete_if{|e| e== []}

fasta.each do |entry| #for each fasta sequence
	#title = entry.definition[entry.definition.index(' ') + 1 .. entry.definition.size]
  title = entry[0].split(' ')[1]
	type = title[0, 1]
	data = entry[1]

	data.tr!("\n\t\r ", '') #get rid of whitespace
	
	if(type == 'A')
		exon2 = score(data, a_exon2_align)
		exon3 = score(data, a_exon3_align)
		if(exon2[0] <= diff_reject and exon3[0] <= diff_reject)
            puts "Approving " + title + ":  " + exon2[0].to_s + " " + exon3[0].to_s
			hla_a.push( [title, exon2[1], exon3[1]] )
		else
			puts "***Rejecting " + title + ":  " + exon2[0].to_s + " " + exon3[0].to_s
		end
	elsif(type == 'B')
		exon2 = score(data, b_exon2_align)
		exon3 = score(data, b_exon3_align)
		if(exon2[0] <= diff_reject and exon3[0] <= diff_reject)
            puts "Approving " + title + ":  " + exon2[0].to_s + " " + exon3[0].to_s
			intron = data[exon2[2] + exon2[1].length .. exon3[2]-1]
			hla_b.push( [title, exon2[1], exon3[1]] )
		else
			puts "***Rejecting " + title + ":  " + exon2[0].to_s + " " + exon3[0].to_s
		end
	elsif(type == 'C')
		exon2 = score(data, c_exon2_align)
		exon3 = score(data, c_exon3_align)
		if(exon2[0] <= diff_reject and exon3[0] <= diff_reject)
            puts "Approving " + title + ":  " + exon2[0].to_s + " " + exon3[0].to_s
			intron = data[exon2[2] + exon2[1].length .. exon3[2]-1]
			hla_c.push( [title, exon2[1], exon3[1]] )
		else
			puts "***Rejecting " + title + ":  " + exon2[0].to_s + " " + exon3[0].to_s
		end
	end
end

#file.close

#Lets sort, just to make things easier for our eyes
hla_a.sort!
hla_b.sort!
hla_c.sort!

File.open('hla_a_std.csv', 'w') do |file|
	hla_a.each do |entry|
		file.puts entry[0] + "," + entry[1] + "," + entry[2]
	end
end
	
File.open('hla_b_std.csv', 'w') do |file|
	hla_b.each do |entry|
		file.puts entry[0] + "," + entry[1] + "," + entry[2] 
	end
end

File.open('hla_c_std.csv', 'w') do |file|
	hla_c.each do |entry|
		file.puts entry[0] + "," + entry[1] + "," + entry[2]
	end
end
