import os, sys, glob
from pathlib import Path
import Galaxy
from Galaxy.core import Template
from libfr_site2 import SiteTemplate
from Galaxy.utils.libtm import TM_result

import multiprocessing as mp
#random.seed(50)
n_proc = 96

chembl_home = '/home/j2ho/DB/ChEMBL/chembl_34'
RCSB_PATH='http://www.rcsb.org/pdb/files/%s.pdb'
PDBNR_LISTFILE = f'{chembl_home}/PDBNR_ALL.list'
# make dictionary of PDBID_CHAINID that has nonmetal ligand deposited in PDB
sitedb = []
# with open('/home/j2ho/GalaxyPipe/db/site_casp16/ligands.240510','r') as f:
#     lines = f.readlines()
# for ln in lines:
#     pdbid_chain = ln.split()[0]
#     ligands = ln.strip().split()[2].split(',')
#     for lig in ligands: 
#         ligtype = lig.split('_')[1]
#         if ligtype == '0': 
#          # 1:lipids, 2:small non-biological, 3: ions(metals), 4: glycans (af3latest)
#         if not pdbid_chain in sitedb:
#             sitedb.append(pdbid_chain)


def read_targets(file_path):
    with open(file_path, 'r') as file:
        targets = file.readlines()
    targets = [target.strip() for target in targets]
    return targets

targets = read_targets(PDBNR_LISTFILE)

with open(f'{chembl_home}/filtered/list_csvs/sitepdb.list','r') as f: 
    for ln in f.readlines(): 
        sitedb.append(ln.strip()) 

# check if in sitedb
def check_holo(pdblist):
    holo_list = []
    hasholo = False
    with open(pdblist, 'r') as f:
        lines = f.readlines()
    for ln in lines:
        if not ln.startswith("#"):
            pdbid = ln.split()[0]
            chainids = ln.split()[3].split('=')[0].split('/')
            for chainid in chainids: 
                pdbid_chain = '%s_%s'%(pdbid, chainid)
                if pdbid_chain in sitedb:
                    hasholo = True
                    holo_list.append(pdbid_chain)
    return hasholo, holo_list


def get_and_align(uniid, holopdb):
    templ_home = f'{chembl_home}/PDBNR/{uniid}'
    templ_home_metal = f'{chembl_home}/PDBNR/{uniid}/pdb_with_metal'
    os.makedirs(templ_home_metal, exist_ok=True) 
    afpdb = f'{templ_home}/afmodel.pdb'
    if not os.path.exists(afpdb):
        os.system(f'wget https://alphafold.ebi.ac.uk/files/AF-{uniid}-F1-model_v4.pdb -O {afpdb}') 
    if Path(afpdb).stat().st_size == 0:
        return False
    else:
        templ = SiteTemplate(holopdb)
        templ.pdb_fn = Galaxy.core.FilePath('%s/%s.pdb'%(uniid, holopdb))
        templ.write(templ_home_metal)
        templ.pdb_fn = Galaxy.core.FilePath('%s.pdb'%(holopdb))
        # sys.exit()
        if not os.path.exists('%s/%s.pdb'%(templ_home_metal,holopdb)): 
            print ("WARNING (TEMPL NOT FOUND):", templ.pdb_fn, 'is not in PDB', flush=True)
        else:
            tm = Galaxy.utils.TM_align(templ.pdb_fn, afpdb)
            templ.tm = tm
            if (len(tm.tr) == 0):
                print ("TM ERROR: no rotation matrix. UNI: %s, TEMP: %s"%(uniid, holopdb),flush=True)
            else: 
                templ.rewrite(templ_home_metal)
        return True

# targets = []
# with open('/home/j2ho/DB/ChEMBL/chembl_34/filtered/list_csvs/all_active_count.list','r') as f: 
#     for line in f:
#         targets.append(line.strip().split()[0]) # line.strip().split()[0] = target name


# targets = []
# with open('/home/j2ho/projects/vs_benchmarks/karma.err.log','r') as f:
#     for ln in f.readlines(): 
#         ln = ln.strip().split()
#         uniid = ln[-1]
#         targets.append(uniid)

def main(uniid):
    haspdb = False
    hasholo = False
    pdblist = f'{chembl_home}/filtered/{uniid}/pdbid.list'
    # if not os.path.exists(pdblist): 
    #     os.system(f'rm -r {target}')
    if os.path.exists(pdblist):
        haspdb = True
        hasholo, holo_list = check_holo(pdblist)
        for holo in holo_list:
            inafdb = get_and_align(uniid, holo)
    else:
        haspdb = False
    
    if hasholo == True:
        if inafdb == False: 
            print ('WARNING: %s has PDB, but not in AFDB'%uniid, flush=True)
        else: 
            print (f'{uniid} {holo_list}', flush=True) 

if __name__=='__main__':
    # holo_list = ['4KS7_A', '4KS8_A']
    # get_and_align('P54149','1FVG_A') 
    # get_and_align('Q9NQU5','4KS8_A')
    # for i, target in enumerate(targets):
    #     print (i, target)
    #     main(target)
    pool = mp.Pool(n_proc)
    pool.map(main,targets)
    
