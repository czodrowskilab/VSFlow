from molvs.standardize import Standardizer
from molvs.tautomer import TautomerCanonicalizer
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs, Torsions


def do_standard_mp(mol: Chem.rdchem.Mol, n: int = 0, ntauts: int = 100) -> tuple:
    """
    Function to standardize and canonicalize and RDKit Mol.
    :param mol: molecule to be standardized
    :type mol: rdkit.Chem.rdchem.Mol
    :param n: molecule number, default = 0
    :type n: int
    :param ntauts: maximum number of tautomers to be enumerated in canonicalization step, default = 100
    :return: tuple containing the molecule number, the standardized molecule and the canonical tautomer (n, Mol, Mol)
    :rtype: tuple
    """
    try:
        mol_sta = Standardizer().charge_parent(Standardizer().fragment_parent(mol), skip_standardize=True)
        mol_can = TautomerCanonicalizer(max_tautomers=ntauts).canonicalize(mol_sta)
        if Chem.MolToSmiles(mol_sta) == Chem.MolToSmiles(mol_can):
            mol_can = mol_sta
    except:
        mol_sta = mol
        mol_can = mol
    return (n, mol_sta, mol_can)


def do_standard(mols, ntauts):
    """
    Adds a standardized molecule and canonicalized tautomer to the database
    :param mols: dictionary containing rdkit.Chem.rdchem.Mol objects as follows: {0: {"mol": rdkit.Chem.rdchem.Mol}, 1: {"mol": rdkit.Chem.rdchem.Mol} ...}
    :type mols: dict
    :param ntauts: maximum number of tautomers to be enumerated in canonicalization step, default = 100
    :type ntauts: int
    """
    for n in mols:
        try:
            mol_sta = Standardizer().charge_parent(Standardizer().fragment_parent(mols[n]["mol"]),
                                                   skip_standardize=True)
            mol_can = TautomerCanonicalizer(max_tautomers=ntauts).canonicalize(mol_sta)
            mols[n]["mol"] = mol_sta
            if Chem.MolToSmiles(mol_sta) == Chem.MolToSmiles(mol_can):
                mols[n]["mol_can"] = mol_sta
            else:
                mols[n]["mol_can"] = mol_can
        except:
            mols[n]["mol_can"] = mols[n]["mol"]


def standardize_mp(mol, n):
    mol_sta = Standardizer().charge_parent(Standardizer().fragment_parent(mol), skip_standardize=True)
    return (n, mol_sta)


def standardize(mols):
    for n in mols:
        mol_sta = Standardizer().charge_parent(Standardizer().fragment_parent(mols[n]["mol"]),
                                               skip_standardize=True)
        mols[n]["mol"] = mol_sta


def canonicalize_mp(mol, n, ntauts):
    mol_can = TautomerCanonicalizer(max_tautomers=ntauts).canonicalize(mol)
    if Chem.MolToSmiles(mol) == Chem.MolToSmiles(mol_can):
        mol_can = mol
    return (n, mol_can)


def canonicalize(mols, ntauts):
    for n in mols:
        mol_can = TautomerCanonicalizer(max_tautomers=ntauts).canonicalize(mols[n]["mol"])
        if Chem.MolToSmiles(mols[n]["mol"]) == Chem.MolToSmiles(mol_can):
            mols[n]["mol_can"] = mols[n]["mol"]
        else:
            mols[n]["mol_can"] = mol_can


def fp_morgan_std_mp(mol, mol_can, i, radius, features, chiral, nBits):
    fp_mol = Chem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, useFeatures=features,
                                                useChirality=chiral)
    if id(mol) == id(mol_can):
        fp_can = fp_mol
    else:
        fp_can = Chem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, useFeatures=features,
                                                    useChirality=chiral)
    return (i, fp_mol, fp_can)


def fp_morgan_std(mols, radius, features, chiral, nBits):
    for i in mols:
        fp_mol = Chem.GetMorganFingerprintAsBitVect(mols[i]["mol"], radius, nBits=nBits, useFeatures=features,
                                                    useChirality=chiral)
        mols[i]["fp"] = fp_mol
        if id(mols[i]["mol"]) == id(mols[i]["mol_can"]):
            mols[i]["fp_can"] = fp_mol
        else:
            mols[i]["fp_can"] = Chem.GetMorganFingerprintAsBitVect(mols[i]["mol_can"], radius, nBits=nBits,
                                                                   useFeatures=features,
                                                                   useChirality=chiral)


def fp_rdkit_std_mp(mol, mol_can, i, nBits):
    fp_mol = Chem.RDKFingerprint(mol, fpSize=nBits)
    if id(mol) == id(mol_can):
        fp_can = fp_mol
    else:
        fp_can = Chem.RDKFingerprint(mol_can, fpSize=nBits)
    return (i, fp_mol, fp_can)


def fp_rdkit_std(mols, nBits):
    for i in mols:
        fp_mol = Chem.RDKFingerprint(mols[i]["mol"], fpSize=nBits)
        mols[i]["fp"] = fp_mol
        if id(mols[i]["mol"]) == id(mols[i]["mol_can"]):
            mols[i]["fp_can"] = fp_mol
        else:
            mols[i]["fp_can"] = Chem.RDKFingerprint(mols[i]["mol_can"], fpSize=nBits)


def fp_tt_std_mp(mol, mol_can, i, nBits, chiral):
    fp_mol = Torsions.GetHashedTopologicalTorsionFingerprint(mol, nBits=nBits, includeChirality=chiral)
    if id(mol) == id(mol_can):
        fp_can = fp_mol
    else:
        fp_can = Torsions.GetHashedTopologicalTorsionFingerprint(mol_can, nBits=nBits, includeChirality=chiral)
    return (i, fp_mol, fp_can)


def fp_tt_std(mols, nBits, chiral):
    for i in mols:
        fp_mol = Torsions.GetHashedTopologicalTorsionFingerprint(mols[i]["mol"], nBits=nBits, includeChirality=chiral)
        mols[i]["fp"] = fp_mol
        if id(mols[i]["mol"]) == id(mols[i]["mol_can"]):
            mols[i]["fp_can"] = fp_mol
        else:
            mols[i]["fp_can"] = Torsions.GetHashedTopologicalTorsionFingerprint(mols[i]["mol_can"], nBits=nBits,
                                                                                includeChirality=chiral)


def fp_ap_std_mp(mol, mol_can, i, nBits, chiral):
    fp_mol = Pairs.GetHashedAtomPairFingerprint(mol, nBits=nBits, includeChirality=chiral)
    if id(mol) == id(mol_can):
        fp_can = fp_mol
    else:
        fp_can = Pairs.GetHashedAtomPairFingerprint(mol_can, nBits=nBits, includeChirality=chiral)
    return (i, fp_mol, fp_can)


def fp_ap_std(mols, nBits, chiral):
    for i in mols:
        fp_mol = Pairs.GetHashedAtomPairFingerprint(mols[i]["mol"], nBits=nBits, includeChirality=chiral)
        mols[i]["fp"] = fp_mol
        if id(mols[i]["mol"]) == id(mols[i]["mol_can"]):
            mols[i]["fp_can"] = fp_mol
        else:
            mols[i]["fp_can"] = Pairs.GetHashedAtomPairFingerprint(mols[i]["mol_can"], nBits=nBits,
                                                                   includeChirality=chiral)


def fp_maccs_std_mp(mol, mol_can, i):
    fp_mol = MACCSkeys.GenMACCSKeys(mol)
    if id(mol) == id(mol_can):
        fp_can = fp_mol
    else:
        fp_can = MACCSkeys.GenMACCSKeys(mol_can)
    return (i, fp_mol, fp_can)


def fp_maccs_std(mols):
    for i in mols:
        fp_mol = MACCSkeys.GenMACCSKeys(mols[i]["mol"])
        mols[i]["fp"] = fp_mol
        if id(mols[i]["mol"]) == id(mols[i]["mol_can"]):
            mols[i]["fp_can"] = fp_mol
        else:
            mols[i]["fp_can"] = MACCSkeys.GenMACCSKeys(mols[i]["mol_can"])


def gen_confs_mp(mol, num, nconfs, seed, threshold, nthreads):
    params = Chem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.pruneRmsThresh = threshold
    params.numThreads = nthreads
    params.randomSeed = seed
    mol_H = Chem.AddHs(mol)
    Chem.EmbedMultipleConfs(mol_H, numConfs=nconfs, params=params)
    try:
        Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=nthreads, maxIters=2000)
    except:
        pass
    return (mol_H, num, Chem.MolToSmiles(mol))


def gen_confs(mols, nconfs, seed, threshold, key, nthreads):
    params = Chem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.pruneRmsThresh = threshold
    params.numThreads = nthreads
    params.randomSeed = seed
    for i in mols:
        mol_H = Chem.AddHs(mols[i][key])
        Chem.EmbedMultipleConfs(mol_H, numConfs=nconfs, params=params)
        try:
            Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=nthreads, maxIters=2000)
        except:
            pass
        mols[i]["confs"] = mol_H
        mols[i]["pattern"] = Chem.MolToSmiles(mols[i][key])


def fp_morgan(mols, radius, features, chiral, nBits):
    for i in mols:
        fp_mol = Chem.GetMorganFingerprintAsBitVect(mols[i]["mol"], radius, nBits=nBits, useFeatures=features,
                                                    useChirality=chiral)
        mols[i]["fp"] = fp_mol
