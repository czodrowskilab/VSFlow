from rdkit.Chem import AllChem as Chem
from molvs.tautomer import TautomerCanonicalizer
from molvs.standardize import Standardizer


def do_standard(mol, dict, ntauts):
    mol_dict = {}
    mol2d = Chem.RemoveHs(mol)
    Chem.Compute2DCoords(mol2d)
    try:
        mol_sta = Standardizer().charge_parent(Standardizer().fragment_parent(mol2d), skip_standardize=True)
        mol_can = TautomerCanonicalizer(max_tautomers=ntauts).canonicalize(mol_sta)
        # mol_dict["mol_sta"] = mol_sta
        # mol_dict["mol_can"] = mol_can
        # mol_dict["props"] = dict
        # mol_dict["mol"] = mol
        mol_dict["mol_sta"] = mol_sta
        mol_dict["mol_can"] = mol_can
        mol_dict["props"] = dict
        mol_dict["mol"] = mol2d
        mol_dict["pattern"] = Chem.MolToSmiles(mol)
        try:
            if mol.GetConformer().Is3D():
                mol_dict["pattern"] = Chem.MolToSmiles(mol)
                mol_dict["confs"] = Chem.AddHs(mol, addCoords=True)
        except ValueError:
            pass
        # return (mol_can, dict)
        return mol_dict
    except:
        # return (mol, dict)
        mol_dict["mol_sta"] = mol2d
        mol_dict["mol_can"] = mol2d
        mol_dict["props"] = dict
        mol_dict["mol"] = mol2d
        #mol_dict["pattern"] = Chem.MolToSmiles(mol)
        try:
            if mol.GetConformer().Is3D():
                mol_dict["pattern"] = Chem.MolToSmiles(mol2d)
                mol_dict["confs"] = Chem.AddHs(mol, addCoords=True)
        except ValueError:
            pass
        # mol_dict["mol_sta"] = mol
        # mol_dict["mol_can"] = mol
        # mol_dict["props"] = dict
        # mol_dict["mol"] = mol
        return mol_dict
        # return ([mol], dict)


def gen_confs_mp(mol, num, nconfs, seed, nthreads):
    params = Chem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    #params.pruneRmsThresh = threshold
    params.numThreads = nthreads
    params.randomSeed = seed
    mol_H = Chem.AddHs(mol)
    Chem.EmbedMultipleConfs(mol_H, numConfs=nconfs, params=params)
    #Chem.EmbedMultipleConfs(mol_H, numConfs=nconfs, randomSeed=seed, ETversion=2, numThreads=0)
    try:
        Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=nthreads, maxIters=2000)
    except:
        pass
    #mol_3D = Chem.RemoveHs(mol_H)
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
        # Chem.EmbedMultipleConfs(mol_H, numConfs=args.nconfs, randomSeed=seed, ETversion=2, numThreads=0)
        try:
            Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=nthreads, maxIters=2000)
        except:
            pass
        # mol_3D = Chem.RemoveHs(mol_H)
        mols[i]["confs"] = mol_H
        mols[i]["pattern"] = Chem.MolToSmiles(mols[i][key])
        counter += 1
        print(counter)