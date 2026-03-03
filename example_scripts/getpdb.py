import sys, os, glob 
import requests
import multiprocessing as mp
#random.seed(50)
n_proc = 96

def getpdblist(uniid):
    pdblist = []
    url = 'https://rest.uniprot.org/uniprot/%s'%uniid
    html = requests.get(url).json()
    for key in html:
        if key == 'uniProtKBCrossReferences':
            structures = html[key]
            for elem in structures:
                if elem['database'] == 'PDB':
                    pdblist.append([elem['id'],elem['properties'][0]['value'],elem['properties'][1]['value'],elem['properties'][2]['value']])
    return pdblist
   
withpdb = [] 
with open('/home/j2ho/DB/uni_pdb_map/uni_to_pdb_chain.csv','r') as f: 
    for ln in f.readlines():
        uniid = ln.strip().split(',')[0]
        withpdb.append(uniid) 

def main(target): 
    pdblist = getpdblist(target)
    wrt = ['# PDBID method resolution chains=residues\n']
    if len(pdblist) > 0:
        print (target)
    if os.path.exists(f'{targetdir}pdbid.list'): 
        return True 
    for pdb in pdblist: 
        pdbid = pdb[0]
        method = pdb[1]
        resolution = pdb[2]
        if resolution != '-':
            resolution = (pdb[2][:-1].strip())
        chains = pdb[3]
        wrt.append('%s  %s  %s  %s\n'%(pdbid, method, resolution, chains))
    with open(f'{targetdir}pdbid.list','wt') as f: 
        f.writelines(wrt)
    return True

if __name__=='__main__':
    pool = mp.Pool(n_proc)
    pool.map(main,targets) 
