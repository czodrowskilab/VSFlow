from rdkit.Chem import AllChem as Chem
from molvs.tautomer import TautomerCanonicalizer
from molvs.standardize import Standardizer


def do_standard_mp(mol, n, ntauts):
    try:
        mol_sta = Standardizer().charge_parent(Standardizer().fragment_parent(mol), skip_standardize=True)
        mol_can = TautomerCanonicalizer(max_tautomers=ntauts).canonicalize(mol_sta)
    except:
        mol_sta = mol
        mol_can = mol
    return (n, mol_sta, mol_can)


def do_standard(mols, ntauts):
    for n in mols:
        try:
            mols[n]["mol_sta"] = Standardizer().charge_parent(Standardizer().fragment_parent(mols[n]["mol"]), skip_standardize=True)
            mols[n]["mol_can"] = TautomerCanonicalizer(max_tautomers=ntauts).canonicalize(mols[n]["mol_sta"])
        except:
            mols[n]["mol_sta"] = mols[n]["mol"]
            mols[n]["mol_can"] = mols[n]["mol"]


def gen_confs_mp(mol, num, nconfs, seed, nthreads):
    params = Chem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    #params.pruneRmsThresh = threshold
    params.numThreads = nthreads
    params.randomSeed = seed
    mol_H = Chem.AddHs(mol)
    Chem.EmbedMultipleConfs(mol_H, numConfs=nconfs, params=params)
    try:
        Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=nthreads, maxIters=2000)
    except:
        pass
    print(num)
    return (mol_H, num, Chem.MolToSmiles(mol))


def gen_confs(mols, nconfs, seed, key, nthreads):
    params = Chem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    # params.pruneRmsThresh = threshold
    params.numThreads = nthreads
    params.randomSeed = seed
    counter = 0
    for i in mols:
        mol_H = Chem.AddHs(mols[i][key])
        Chem.EmbedMultipleConfs(mol_H, numConfs=nconfs, params=params)
        try:
            Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=nthreads, maxIters=2000)
        except:
            pass
        mols[i]["confs"] = mol_H
        mols[i]["pattern"] = Chem.MolToSmiles(mols[i][key])
        counter += 1
        print(counter)