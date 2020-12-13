#! /usr/bin/python3
#! /usr/bin/python3

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import sys
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from forcefield import VdwParamset

NACCESS_BIN = '/home/anabel/Documents/BIOPH2020/soft/NACCESS/naccess'

parser = argparse.ArgumentParser(
    prog='structure_setup',
    description='basic structure setup'
)

parser.add_argument(
    '--naccess',
    action='store',
    dest='naccess_bin',
    default=os.path.dirname(os.path.abspath(__file__)) + '/soft/NACCESS/naccess',
    help='Vdw parameters'
)
parser.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default=os.path.dirname(os.path.abspath(__file__)) + '/data/vdwprm',
    help='Vdw parameters'
)

parser.add_argument('pdb_file',help='Input PDB', type=open)
parser.add_argument('pdbqt_file',help='Input PDBQT', type=open)

args = parser.parse_args()

print("PDB File:", args.pdb_file.name)
print("PDBQT File:", args.pdbqt_file.name)
print("Parameters File:", args.vdwprm_file)

# loading VdW parameters
ff_params = VdwParamset(args.vdwprm_file)

parser = PDBParser(PERMISSIVE=1)
print('Parsing PDB', args.pdb_file.name)

# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file.name)

# We will use the xtra attribute in Bio.PDB.Atom to hold the new data
# Getting Charges and Atom type from PDBQT
print('Parsing PDBQT', args.pdbqt_file.name)
params=[{}]

#Fix aton numbers when they do not start in 1
i = 1
for at in st.get_atoms():
    at.serial_number = i
    i += 1

for line in args.pdbqt_file:
    line = line.rstrip()
    #Skip TER records from PDBQT
    if line.find('TER') != -1:
        continue
    params.append({'charge': line[69:76], 'type': line[77:].replace(' ','')})

total_charge = 0.
for at in st.get_atoms():
    at.xtra['atom_type'] = params[at.serial_number]['type']
    at.xtra['charge'] = float(params[at.serial_number]['charge'])
    at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]
    total_charge += at.xtra['charge']
print('Total Charge: {:8.2f}'.format(total_charge))

# Calculating surfaces
# Srf goes to .xtra['EXP_NACCESS'] field
srf = NACCESS_atomic(st[0], naccess_binary=args.naccess_bin)



###Define functions that I will need for exercise 3:
def MH_diel(r):
    '''Mehler-Solmajer dielectric'''
    return 86.9525 / (1 - 7.7839 * math.exp(-0.3153 * r)) - 8.5525

def elec_int(at1, at2, r):
    return 332.16 * at1.xtra['charge'] * at2.xtra['charge'] / MH_diel(r) / r

def vdw_int(at1, at2, r):
    '''Vdw interaction energy between two atoms'''
    eps12 = math.sqrt(at1.xtra['vdw'].eps * at2.xtra['vdw'].eps)
    sig12_2 = at1.xtra['vdw'].sig * at2.xtra['vdw'].sig
    return 4 * eps12 * (sig12_2*6/r12 - sig12_23/r*6)

###I have transformed those functions of github gelpi repositorie to two functions that valorate these three at the same time.
def electrostatics(atom_1,atom_2): 
    diel_list = []
	r = atom_2-atom_1 
	dielectric = [1,80,86.9525/(1-7.7839*exp(-0.3153*r))-8.5525]
	for x in dielectric:
		electro= 332.16*((atom_1.xtra["charge"]*atom_2.xtra["charge"])/(x*r))
		diel_list.append(electro)
	return diel_list


def vdw(atom_1,atom_2):
	r = atom_1-atom_2
	epsilon = sqrt(atom_1.xtra["vdw"].eps * atom_2.xtra["vdw"].eps)
	sigma = sqrt(atom_1.xtra["vdw"].sig * atom_2.xtra["vdw"].sig) 
	vdw = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)
	return vdw


van_der_waals = 0
total_electro= [0,0,0]


for residue in st.get_residues():
	total_res= [0,0,0]
	van_der_waals = 0
	for atom_1 in residue.get_atoms():
		for atom_2 in st.get_atoms():
			if atom_2.get_parent() != residue:
				electro1= electrostatics(atom_1,atom_2)
				vdw__ = vdw(atom_1,atom_2)
				van_der_waals += vdw__
				res_vdw += vdw__
				for x in range(3): #because we want distance > 2 and range is -1!
					total_electro[x] += electro[x]/2
					total_res[x] += electro[x]/2

print("Total electrostatic of the molecule is :", total_electrostatic)
print("Total van der waals of the molecule is :", van_der_waals)




###Functions needed on exercise 4:
solvation = 0

'''def calc_solvation(st, res):
    solv = 0
    for at in res.get_atoms():
        if 'EXP_NACCESS' not in at.xtra:
            continue
        s = float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf
        solv += s
        if at.id in ala_atoms:
            solv_ala += s
    return solv, solv_ala'''


def solvation_energy(atom):
    sigma1 = atom.xtra['vdw'].fsrf
    asa_approx = atom.xtra['EXP_NACCESS']
    return sigma1*float(asa_approx) #= solvation

solv_energy = 0

for res in st.get_residues():
    resolv = 0
    for at in res.get_atoms():
        if at.element != "H":
            solv = solvation_energy(at)
            solv_energy += solv
            resolv += solv
    num = num + 1

print("Total solvation is :",solv_energy)
