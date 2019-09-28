#!/usr/bin/env python3
#%%
import sys, argparse, csv, re
import numpy as np
import pandas as pd
from datetime import date
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

base_pairs = [('A','U'),('G','C'),('G','U')]

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Create a 1D topology map for polymer connectivity.')
	contact_file = parser.add_mutually_exclusive_group(required=True)
	contact_file.add_argument('-fr3d','--contacts_path', help='Path to fr3d contacts file from RiboVision.')
	contact_file.add_argument('-jv','--jalview_path', help='Path to Jalview annotation file.')
	contact_file.add_argument('-bpseq','--bpseq_path', help='Path to bpseq file from XRNA.')
	parser.add_argument('-o','--output_path', help='Path to output file.', type=str)
	parser.add_argument('-cm','--color_map', help='Colormap to use.', type=str, default="viridis")
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def read_csv(csv_path):
	'''	Reads fr3d csv; ignores the first line; returns list of base-pairing positions.
	'''
	cannonical_tuples=[]
	with open(csv_path, encoding="utf-16") as csv_file:
		csv_reader = csv.reader(csv_file, delimiter='\t')
		line_count = 0
		for row in csv_reader:
			if line_count == 0:
				line_count += 1
				pass
			else:
				line_count += 1
				for cannonical_bp in base_pairs:
					cannonical_tuples.append((int(row[0].split(':')[1]),int(row[2].split(':')[1])))
	return cannonical_tuples

def read_bpseq(bpseq_path):
	'''	Reads a bpseq contacts file for XRNA; returns list of base-pairing positions.
	'''
	cannonical_tuples=[]
	with open(bpseq_path, encoding="utf8") as bp_file:
		for row in bp_file:
			if int(row.split()[0]) == 0 or int(row.split()[2]) == 0:
				continue
			elif len(cannonical_tuples) < 1:
				cannonical_tuples.append((int(row.split()[0]),int(row.split()[2])))
			elif (row.split()[0],row.split()[2]) in cannonical_tuples or (row.split()[2],row.split()[0]) in cannonical_tuples:
				pass
			else:
				cannonical_tuples.append((int(row.split()[0]),int(row.split()[2])))
	return cannonical_tuples

def read_jv(jalview_path):
	'''Reads a jalview annotation file for helices; returns list of base-pairing positions.
		Assumes just one annotation line.
	'''
	with open(jalview_path, encoding="utf8") as jv_file:
		for line in jv_file:
			if re.match('^NO_GRAPH	Secondary structure', line):
				anno_resi = [w.replace(',[000000]', '') for w in line.split('\t')[2].split('|')]

	return find_parens(anno_resi)

def find_parens(s):
	toret = {}
	pstack = []
		
	for i, c in enumerate(s):
		if c == '(':
			pstack.append(i)
		elif c == ')':
			if len(pstack) == 0:
				raise IndexError("No matching closing parens at: " + str(i))
			toret[pstack.pop()] = i
	if len(pstack) > 0:
		raise IndexError("No matching opening parens at: " + str(pstack.pop()))
	tupled_toret = []
	for k,v in toret.items():
		tupled_toret.append((k,v))
	return tupled_toret
	
def create_helices(cann_tups):
	'''
	'''
	ordered_basepairs = sorted(list(set(cann_tups)))
	helix_num=0
	helix_basepairs={}
	#print(ordered_basepairs)
	for bp_num in range(0,len(ordered_basepairs)):
		if bp_num > 0:		#Think how to make that to allow for up to 3 unpaired residues...
			if (ordered_basepairs[bp_num][0]-1 == ordered_basepairs[bp_num-1][0]) or (ordered_basepairs[bp_num][1]-1 == ordered_basepairs[bp_num-1][1]):
				helix_basepairs[helix_num].append(ordered_basepairs[bp_num])
			else:
				helix_num+=1
				helix_basepairs[helix_num]=[]
		else:		#Possible bug here if the first element is not part of the first helix
			helix_basepairs[helix_num]=[]
			helix_basepairs[helix_num].append(ordered_basepairs[bp_num])
	return helix_basepairs, helix_num

def main(commandline_args):
	comm_args = create_and_parse_argument_options(commandline_args)
	if comm_args.contacts_path:
		cann_tups = read_csv(comm_args.contacts_path)
	if comm_args.bpseq_path:
		cann_tups = read_bpseq(comm_args.bpseq_path)
	if comm_args.jalview_path:
		cann_tups = read_jv(comm_args.jalview_path)
	helix_basepairs, helix_num = create_helices(cann_tups)
	
	#Get the maximum values for the residue numbers
	maxresi = max([max(max(cann_tups,key=lambda item:item[0])),max(max(cann_tups,key=lambda item:item[1]))])
	minresi = min([min(min(cann_tups,key=lambda item:item[0])),min(min(cann_tups,key=lambda item:item[1]))])

	#Plotting
	cmap = plt.cm.get_cmap(comm_args.color_map)
	hexes = cmap(np.linspace(0, 1, helix_num+1))
	plt.style.use('default')
	fig, ax = plt.subplots(1,1)
	ax.get_yaxis().set_visible(False)
	for helices in sorted(helix_basepairs.keys()):
		for basepair in helix_basepairs[helices]:
			pac = mpatches.Arc([sum(basepair)/2,0], abs(basepair[0]-basepair[1]), abs(basepair[0]-basepair[1]), angle=180, color=hexes[helices])
			ax.add_patch(pac)
	ax.axis([minresi-(maxresi*0.01),maxresi*1.01,minresi,maxresi/1.5])
	fig.canvas.draw()

#%%

#Execute script within a Jupyter notebook
main(['-jv','/home/ppenev/Dropbox-Gatech/Programs/topology_1d/test_data/ES39_helices_v4','-cm', 'tab10'])
#main(['-fr3d','/home/ppenev/Dropbox-Gatech/Programs/topology_1d/test_data/ECOLI_5S.csv', '-cm', 'rainbow'])
#main(['-bpseq','/home/ppenev/Dropbox-Gatech/Programs/topology_1d/test_data/PYRFU.bpseq','-cm', 'plasma'])

#%%
if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))

#%%
