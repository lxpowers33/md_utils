import csv

def run_compare(condition1, condition2):
	out = compare_conditions(condition1, condition2)
	with open('testing/dataout.tsv', 'w') as csvfile:
		csvout = csv.writer(csvfile, delimiter='\t')
		for row in out:
			print(row)
			csvout.writerow(row)

def run_compare_test():
	condition1 = ['testing/test1_contacts.tsv','testing/test2_contacts.tsv']
	condition2 = ['testing/test1_contacts.tsv','testing/test3_contacts.tsv']
	run_compare(condition1, condition2)

def compare_conditions(files1, files2):
	freq1 = frequencies_over_reps(files1)
	freq2 = frequencies_over_reps(files2)
	#make an intersection of the two lists
	keys = Union(freq1.keys(),freq2.keys())
	compare = []
	for key in keys:
		compare.append([key, freq1.get(key, 0), freq2.get(key, 0)])
	out = sorted(compare, key=lambda x: x[1], reverse=True)
	return out

def Union(lst1, lst2):
    final_list = list(set(lst1) | set(lst2))
    return final_list

def frequencies_over_reps(files):
	"""
	output dictionairy of contacts, frequency, standard deviation
	"""
	totals = {}
	for file in files:
		file = open(file, "r")
		out = frequencies_from_contacts(file)
		#aggregate with previous
		for key in out:
			if (key in totals): 
				totals[key] = totals[key] + out[key]
			else: 
				totals[key] = out[key]

	frequency = {k: v/len(files)*100.0 for k, v in totals.items()}
	return frequency

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

	frequency = {k: v/1.0/nframes for k, v in counts.items()}
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
	run_compare_test()
