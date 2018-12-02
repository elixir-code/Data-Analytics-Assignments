import pandas as pd
import numpy as np

from matplotlib import pyplot as plt

''' Read DNA sequence strings into an array from fasta file '''
def readDNASequences(filename):

	dna_sequences = []
	with open(filename, 'r') as f:	
		dna_sequence = ''
		for line in f:
			if line.startswith('>'):
				dna_sequences.append(dna_sequence)
				dna_sequence = ''
			else:
				dna_sequence += line.rstrip('\n')
	dna_sequences.pop(0)

	return dna_sequences

''' Counts the number of ORFs and their lengths for each DNA sequence '''
def countORFLengths(dna_sequences):

	''' Processing Open Read Frames (ORFs) in dna sequence '''
	dna_seq_summaries = []

	for index, dna_sequence in enumerate(dna_sequences):
		dna_seq_summary = {'n_genes':0, 'lengths_genes':[], 'avg_length':[]}

		i = 0
		while i<len(dna_sequence)-2:
			if dna_sequence[i:i+3] == 'ATG':
				j = i+3
	
				while j<len(dna_sequence)-2:
					if dna_sequence[j:j+3] in ['TAA', 'TAG', 'TGA']:
						dna_seq_summary['n_genes'] += 1
						dna_seq_summary['lengths_genes'].append(j-i+3)
						break
					j += 1
				i = j+3
			else:	
				i += 1

		dna_seq_summary['avg_length'] = np.mean(dna_seq_summary['lengths_genes'])
		
		dna_seq_summaries.append(dna_seq_summary)

	return dna_seq_summaries

''' count the percentages of A, G, C, T in each dna sequence '''
def countAGCT(dna_sequences):

	dna_sequences_agct = []

	for dna_sequence in dna_sequences:
		dna_sequence_agct = [0, 0, 0, 0, 0]
		for char in dna_sequence:
			if char == 'A':
				dna_sequence_agct[0] += 1
			elif char == 'G':
				dna_sequence_agct[1] += 1
			elif char == 'C':
				dna_sequence_agct[2] += 1
			elif char == 'T':
				dna_sequence_agct[3] += 1

			else:
				dna_sequence_agct[4] += 1

		length_dna_seq = sum(dna_sequence_agct)
		dna_sequence_agct[0] = (dna_sequence_agct[0]/length_dna_seq)*100
		dna_sequence_agct[1] = (dna_sequence_agct[1]/length_dna_seq)*100
		dna_sequence_agct[2] = (dna_sequence_agct[2]/length_dna_seq)*100
		dna_sequence_agct[3] = (dna_sequence_agct[3]/length_dna_seq)*100
		dna_sequence_agct[4] = (dna_sequence_agct[4]/length_dna_seq)*100

		dna_sequences_agct.append(dna_sequence_agct)

	return dna_sequences_agct

''' count number of occurenence of certain patterns in a gene sequqnce '''
def countOccurences(dna_sequence, pattern):

	count = 0

	for i in range(len(dna_sequence)-len(pattern)+1):
		if dna_sequence[i:i+len(pattern)] == pattern:
			count += 1

	return count

''' Count number of occurrences of a pattern for each dna sequence '''
def countDNASeqOccurrences(dna_sequences, pattern):

	return [countOccurences(dna_sequence, pattern) for dna_sequence in dna_sequences]


if __name__ == '__main__':

	dna_sequences = readDNASequences('DNASequences.fasta')

	# Table 0
	length_dna_sequences = [len(dna_sequence) for dna_sequence in dna_sequences]
	
	# Table 1
	dna_sequence_orfs = countORFLengths(dna_sequences)

	# Table 2
	dna_agct_percents = np.array(countAGCT(dna_sequences), dtype=np.float16)

	# Table 3
	cag_n_occurences = np.array(countDNASeqOccurrences(dna_sequences, 'CAG'), dtype='int')
	ttaggg_n_occurences = np.array(countDNASeqOccurrences(dna_sequences, 'TTAGGG'), dtype='int')
	ac_n_occurences = np.array(countDNASeqOccurrences(dna_sequences, 'AC'), dtype='int')

	""" Export all data as excel """
	writer = pd.ExcelWriter('output.xlsx')

	lengths_df = pd.DataFrame(length_dna_sequences, index=range(1,len(dna_sequences)+1), columns=['Length of DNA sequence'])
	lengths_df.to_excel(writer, 'DNA Sequence Length')

	dna_orf_count_length = pd.DataFrame(dna_sequence_orfs, index=range(1,len(dna_sequences)+1))
	dna_orf_count_length.to_excel(writer, 'ORFs')

	dna_agct_df = pd.DataFrame(dna_agct_percents, index=range(1,len(dna_sequences)+1), columns=['Percentage A', 'Percentage G', 'Percentage C', 'Percentage T', 'Percentage Others'])
	dna_agct_df.to_excel(writer, 'AGCT Percentages')

	dna_pattern_counts_df = pd.DataFrame({'Count CAG': cag_n_occurences, 'Count TTAGGG': ttaggg_n_occurences, 'Count AC': ac_n_occurences}, index=range(1,len(dna_sequences)+1))
	dna_pattern_counts_df.to_excel(writer, 'Count Patterns')

	writer.save()


	""" All plots with ACGT percentages """

	''' Plot charts for all data '''
	dna_agct_means = np.mean(dna_agct_percents, axis=0)[:4]
	dna_agct_labels = 'Adenine %', 'Guanine %', 'Cytosine %', 'Thymine %'

	fig, ax1 = plt.subplots()

	# ax1.set_title('AGCT Percentages - Pie Chart')
	ax1.pie(dna_agct_means, labels=dna_agct_labels, autopct='%1.1f%%', shadow=True, startangle=90)
	ax1.axis('equal')

	plt.show()

	''' Box plot for each AGCT percentages data '''
	dna_agct_labels = 'Adenine %', 'Guanine %', 'Cytosine %', 'Thymine %'
	
	fig, ax1 = plt.subplots()

	ax1.boxplot(dna_agct_percents[:,:4], labels=dna_agct_labels)

	plt.show()

	''' Histogram for AGCT percentages '''
	fig, axs = plt.subplots(2, 2, sharey=True, tight_layout=True)

	axs[0,0].set_xlabel("Percentage of 'A'")
	axs[0,0].set_ylabel("Count of DNA sequences")
	axs[0,0].hist(dna_agct_percents[:, 0])

	axs[0,1].set_xlabel("Percentage of 'G'")
	axs[0,1].set_ylabel("Count of DNA sequences")
	axs[0,1].hist(dna_agct_percents[:, 1])

	axs[1,0].set_xlabel("Percentage of 'C'")
	axs[1,0].set_ylabel("Count of DNA sequences")
	axs[1,0].hist(dna_agct_percents[:, 2])

	axs[1,1].set_xlabel("Percentage of 'T'")
	axs[1,1].set_ylabel("Count of DNA sequences")
	axs[1,1].hist(dna_agct_percents[:, 3])

	plt.show()

	''' Violin plot for AGCT percentages '''

	fig, axs = plt.subplots(2, 2, sharey=True, tight_layout=True)

	axs[0,0].set_title('Adenine %')
	axs[0,0].violinplot(dna_agct_percents[:, 0], showmeans=True, showextrema=True, showmedians=True)

	axs[0,1].set_title('Guanine %')
	axs[0,1].violinplot(dna_agct_percents[:, 1], showmeans=True, showextrema=True, showmedians=True)

	axs[1,0].set_title('Cytosine %')
	axs[1,0].violinplot(dna_agct_percents[:, 2], showmeans=True, showextrema=True, showmedians=True)

	axs[1,1].set_title('Thymine %')
	axs[1,1].violinplot(dna_agct_percents[:, 3], showmeans=True, showextrema=True, showmedians=True)

	plt.show()

	""" All plots with count of patterns AC, CAG, TTAGGG """

	fig, ax = plt.subplots()

	n_dna_sequence_samples = 10
	dna_sequence_samples = np.random.choice(np.arange(len(dna_sequences)), n_dna_sequence_samples, replace=False)

	bar_width = 0.3
	opacity = 0.4

	rects_ac = ax.bar(np.arange(n_dna_sequence_samples)-bar_width, ac_n_occurences[dna_sequence_samples], bar_width, alpha=opacity, color='b', label="pattern 'AC'")
	rects_cag = ax.bar(np.arange(n_dna_sequence_samples), cag_n_occurences[dna_sequence_samples], bar_width, alpha=opacity, color='r', label="pattern 'CAG'")
	rects_ttaggg = ax.bar(np.arange(n_dna_sequence_samples)+bar_width, ttaggg_n_occurences[dna_sequence_samples], bar_width, alpha=1.0, color='g', label="pattern 'TTAGGG'")

	ax.set_xlabel('DNA Sequences')
	ax.set_ylabel('Pattern Occurrence counts')

	ax.set_xticks(np.arange(n_dna_sequence_samples))
	ax.set_xticklabels(np.arange(n_dna_sequence_samples)+1)

	ax.legend()

	fig.tight_layout()
	
	plt.show()


	''' Plot histogram for each pattern 'AC', 'CAG', 'TTAGGG' '''
	fig, axs = plt.subplots(1, 3, sharey=True, tight_layout=True)

	#axs[0].set_title("Number of occurrences of 'AC'")
	axs[0].set_xlabel("Number of occurrences of 'AC'")
	axs[0].set_ylabel("Frequency of DNA sequences")
	axs[0].hist(ac_n_occurences)

	#axs[1].set_title("Number of occurrences of 'CAG'")
	axs[1].set_xlabel("Number of occurrences of 'CAG'")
	axs[1].set_ylabel("Frequency of DNA sequences")
	axs[1].hist(cag_n_occurences)

	#axs[2].set_title("Number of occurrences of 'TTAGGG'")
	axs[2].set_xlabel("Number of occurrences of 'TTAGGG'")
	axs[2].set_ylabel("Frequency of DNA sequences")
	axs[2].hist(ttaggg_n_occurences)

	plt.show()

	''' Box plot for each each pattern 'AC', 'CAG', 'TTAGGG' '''
	dna_seq_pattern_labels = "Pattern 'AC'", "Pattern 'CAG'", "Pattern 'TTTAGG'"

	fig, ax1 = plt.subplots()

	# ax1.set_title("Number of occurrences of patterns 'AC', 'CAG', 'TTAGGG'")
	ax1.set_ylabel('Number of occurrences of pattern in a DNA sequence')
	ax1.boxplot([ac_n_occurences, cag_n_occurences, ttaggg_n_occurences], labels=dna_seq_pattern_labels)

	plt.show()

	''' Violin plot for each each pattern 'AC', 'CAG', 'TTAGGG' '''
	fig, axs = plt.subplots(1, 3, sharey=True, tight_layout=True)

	axs[0].set_title("Np. of occurrences of 'AG'")
	axs[0].violinplot(ac_n_occurences, showmeans=True, showextrema=True, showmedians=True)

	axs[1].set_title("No. of occurrences of 'CAG'")
	axs[1].violinplot(cag_n_occurences, showmeans=True, showextrema=True, showmedians=True)

	axs[2].set_title("No. of ccurrences of 'TTAGGG'")
	axs[2].violinplot(ttaggg_n_occurences, showmeans=True, showextrema=True, showmedians=True)

	plt.show()

	''' Scatter plot for each each pattern 'AC', 'CAG', 'TTAGGG' '''
	fig, ax1 = plt.subplots()

	# ax1.set_title("Number of occurrences of patterns 'AC' and 'CAG'")
	ax1.set_xlabel("Number of occurrences of pattern 'AC'")
	ax1.set_ylabel("Number of occurrences of pattern 'CAG'")

	ax1.scatter(ac_n_occurences, cag_n_occurences)
	plt.show()

	""" Draw plots for some selected dna gene sequence ORF length """
	choosen_sequence_index = np.random.randint(len(dna_sequences))
	choosen_sequence_orf_lens = dna_sequence_orfs[choosen_sequence_index]['lengths_genes']

	''' Draw box plot for the orf lengths for choosen dna sequence'''
	fig, ax1 = plt.subplots()

	ax1.set_title("Length of ORFs for the {}th dna sequence".format(choosen_sequence_index+1))
	ax1.boxplot(choosen_sequence_orf_lens)

	plt.show()

	''' Draw histogram for the orf lengths for choosen dna sequence'''
	fig, ax1 = plt.subplots()

	ax1.set_title("Length of ORFs for the {}th dna sequence".format(choosen_sequence_index+1))
	ax1.hist(choosen_sequence_orf_lens)

	plt.show()

	''' Draw violin plot for the orf lengths for choosen dna sequence'''
	fig, ax1 = plt.subplots()

	ax1.set_title("Length of ORFs for the {}th dna sequence".format(choosen_sequence_index+1))
	ax1.violinplot(choosen_sequence_orf_lens, showmeans=True, showextrema=True, showmedians=True)

	plt.show()

	""" Draw plots length of each dna sequence """

	''' Draw box plot for the length of each dna sequence '''
	fig, ax1 = plt.subplots()

	ax1.set_ylabel("Length of DNA sequence")
	ax1.boxplot(length_dna_sequences)

	plt.show()

	''' Draw histogram for the length of each dna sequence '''
	fig, ax1 = plt.subplots()

	ax1.set_xlabel("Length of DNA sequences")
	ax1.set_ylabel("Frequency of DNA sequences")
	ax1.hist(length_dna_sequences)

	plt.show()

	''' Draw histogram for the length of each dna sequence '''
	fig, ax1 = plt.subplots()

	ax1.violinplot(length_dna_sequences, showmeans=True, showextrema=True, showmedians=True)

	plt.show()


