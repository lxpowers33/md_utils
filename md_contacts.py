import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("talk", font_scale=1.4)
sns.set_style("whitegrid")
import statistics

from scipy.stats import ttest_ind_from_stats

def plot_compare():
	x = []
	y = [] 
	xerr = []
	yerr = []
	i = j = 0
	with open('data.tsv','r') as tsvin:
		file = csv.reader(tsvin, delimiter='\t')
		for row in file:
			zero_var = ((float(row[2]) == 0) & (float(row[5]) == 0))
			if (not zero_var):
				try: 
					test_result = ttest_ind_from_stats(float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]), float(row[6]), equal_var=False)
				except ZeroDivisionError:
					print('Oops! ERROR') 

			i = i + 1
			to_small = (float(row[1]) < 5) & (float(row[4]) < 5)
			if ((test_result[1] < 0.05) & ((not to_small) & (not zero_var))):
				if ('AP8' in row[0]):
					j = j + 1
				x.append(float(row[1]))
				y.append(float(row[4]))
				xerr.append(float(row[2]))
				yerr.append(float(row[5]))

	print(i)
	print(len(x))
	print(j)
	#plt.scatter(x, y, alpha=0.5)
	plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o')
	plt.xlabel('contact % (G103E)')
	plt.ylabel('contact % (+AP8 G103E)')
	plt.ylim([0, 100])
	plt.xlim([0, 100])
	plt.show()


def run_filter(input_file='data.tsv', save_name='data_filtered_p10.tsv'):
	out = filter_compare(input_file)
	with open(save_name, 'w') as csvout:
		csvout = csv.writer(csvout, delimiter='\t')
		for row in out:
			csvout.writerow(row)
	return out

def filter_compare(input_file, p_value = 0.10):
	filtered = []
	with open(input_file,'r') as tsvin:
		tsvin = csv.reader(tsvin, delimiter='\t')
		for row in tsvin:
			zero_var = ((float(row[2]) == 0) & (float(row[5]) == 0))
			if (not zero_var):
				test_result = ttest_ind_from_stats(float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]), float(row[6]), equal_var=False)
			to_small = (float(row[1]) < 5) & (float(row[4]) < 5)
			if ((test_result[1] < p_value) & ((not to_small) & (not zero_var))):
				#write row out
				if (float(row[1]) < float(row[4])):
					if (float(row[1]) == 0):
						ratio = -1000
					else:
						ratio = -float(row[4])/float(row[1])
				else:
					if (float(row[4]) == 0):
						ratio = 1000
					else:
						ratio = float(row[1])/float(row[4])

				item = [row[0], row[1][0:6], row[2][0:6], row[3][0:6], row[4][0:6], row[5][0:6], row[6][0:6], ratio]
				filtered.append(item)

	out = sorted(filtered, key=lambda x: abs(x[7]), reverse=True)
	return out
	
					

def run_compare(conditions,name):
	out = compare_conditions(conditions)
	with open(name, 'w') as csvfile:
		csvout = csv.writer(csvfile, delimiter='\t')
		for row in out:
			print(row)
			csvout.writerow(row)
	return out

def run_compare_test():
	'''
	Test the computation of average and stardard deviation for 2 conditions
	'''
	condition1 = ['testing/test1_contacts.tsv','testing/test2_contacts.tsv']
	condition2 = ['testing/test1_contacts.tsv','testing/test3_contacts.tsv']
	out = run_compare([condition1, condition2], 'testing/dataout.tsv')
	test = [['v--A--A', 100.0, 0.0, 2, 50.0, 70.71067811865476, 2],
	['v--C--C', 37.5, 17.67766952966369, 2, 37.5, 17.67766952966369, 2],
	['v--B--B', 25.0, 35.35533905932738, 2, 50.0, 0.0, 2]]
	for i in range(0,3): 
		if (out[i] != test[i]):
			print('Test Failed')
			return
	print('Test Passed')


def compare_conditions(condition_files):
	'''
	Input list of list of strings, each list is strings for file names of condition
		[['1/file1.txt', '1/file2.txt'], ['2/file1.txt', '2/file2.txt']]
	Output: list of lists 
		each list contains [interaction key, freq in 1, std in 1, freq in 2, std in 2, ...]
	note the statistics module gives a correct standard deviation (1/N-1 for the averaging)
	'''
	conditions = []
	all_keys = []
	for files in condition_files:
		data = frequencies_over_reps(files)
		conditions.append(data)
		all_keys.append(data.keys())

	unique_keys = set().union(*all_keys)

	compare = []
	for key in unique_keys:
		keyitem = [key]
		for i in range(0, len(condition_files)):
			freq = conditions[i].get(key, [0, 0, 0, 0, 0])
			nobs = len(condition_files[i])
			freq = fill_blanks(freq, nobs)
			keyitem.append(statistics.mean(freq))
			keyitem.append(statistics.stdev(freq))
			keyitem.append(nobs)
		compare.append(keyitem)

	#sort values by the mean frequency of the first condition
	out = sorted(compare, key=lambda x: x[1], reverse=True)

	return out

def fill_blanks(partial_list, expected_length):
	for i in range(0, expected_length-len(partial_list)):
		partial_list.append(0)
	return partial_list

def Union(lst1, lst2):
    final_list = list(set(lst1) | set(lst2))
    return final_list

def frequencies_over_reps(files):
	"""
	output dictionairy of contacts
	key is contact key type--atom1--atom2
	value is list of frequencies for each replicate
	"""
	totals = {}
	for file in files:
		file = open(file, "r")
		out = frequencies_from_contacts(file)
		#aggregate with previous
		for key in out:
			if (key in totals): 
				totals[key].append(out[key])
			else: 
				totals[key] = [out[key]]

	return totals

def frequencies_from_contacts(f):
	"""
	file = open("ternary_AP8_MK6_G103E_v1_1_contacts.tsv", "r")
	compute_frequencies(file)

	Parameters
		f: an open file
	Returns
		dictionairy mapping contact to frequency
	"""
	#read f to the list format
	all_contacts, nframes = parse_contacts(f)
	f.close()
	counts = {}
	for item in all_contacts:
		key = item[1]+'--'+item[2]+'--'+item[3]
		if (key in counts): 
			counts[key] = counts[key] + 1
		else: 
			counts[key] = 1

	frequency = {k: v*100.0/nframes for k, v in counts.items()}
	return frequency



def parse_contacts(input_lines, itypes=None):
    """
    Read a contact-file (tab-separated file with columns: frame, i-type, atomid1, atomid2[, atomid3[, atomid4]] and
    return it as a list of lists with frames converted to ints. The total number of frames is also returned.
    Example
    -------
        parse_contacts([
            "# total_frames:2\n",
            "0  hbbb    A:ALA:1:N   A:THR:10:O\n",
            "0  vdw     A:ALA:1:CB  B:CYS:3:H\n",
            "1  vdw     A:ALA:1:N   A:THR:10:C\n"
        ])
        # returns:
        # ([
        #        [0, "hbbb", "A:ALA:1:N", "A:THR:10:O"],
        #        [0, "vdw", "A:ALA:1:CB", "B:CYS:3:H"],
        #        [1, "vdw", "A:ALA:1:N", "A:THR:10:C"]
        #  ], 2)
    Parameters
    ----------
    input_lines: iterable
        Iterator of over a set of strings. Can be a file-handle
    itypes: set of str | None
        Interactions to include in the output
    Returns
    -------
    (list of list, int)
        The list of interactions and the total number of frames
    """
    ret = []
    total_frames = 0
    for line in input_lines:
        line = line.strip()
        if "total_frames" in line:
            tokens = line.split(" ")
            total_frames = int(tokens[1][tokens[1].find(":")+1:])

        if len(line) == 0 or line[0] == "#":
            continue

        tokens = line.split("\t")
        tokens[0] = int(tokens[0])

        if itypes is None or tokens[1] in itypes:
            ret.append(tokens)

    return ret, total_frames

if __name__== "__main__":
	#run_compare_test()
	#plot_compare()
	run_filter()
