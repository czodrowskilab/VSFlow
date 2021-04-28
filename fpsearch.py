from rdkit.Chem import AllChem as Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit import DataStructs
from rdkit.Chem import Descriptors
from itertools import groupby
import copy

sim_dict = {"tan": DataStructs.TanimotoSimilarity, "dice": DataStructs.DiceSimilarity, "cos": DataStructs.CosineSimilarity,
            "sok": DataStructs.SokalSimilarity, "russ": DataStructs.RusselSimilarity, "kulc": DataStructs.KulczynskiSimilarity,
            "mcco": DataStructs.McConnaugheySimilarity, "tver": DataStructs.TverskySimilarity}

sim_dict_name = {"tan": "TanimotoSimilarity", "dice": "DiceSimilarity", "cos": "CosineSimilarity",
                 "sok": "SokalSimilarity", "russ": "RusselSimilarity", "kulc": "KulczynskiSimilarity",
                 "mcco": "McConnaugheySimilarity", "tver": "TverskySimilarity"}


def filter_res(mol, filter_dict):
    filter_func = {"mw": Descriptors.MolWt, "logp": Descriptors.MolLogP, "tpsa": Descriptors.TPSA,
                   "hdon": Descriptors.NumHDonors, "hacc": Descriptors.NumHAcceptors,
                   "rotb": Descriptors.NumRotatableBonds, "narom": Descriptors.NumAromaticRings,
                   "nhet": Descriptors.NumAromaticHeterocycles}
    filt_mol = mol
    for prop in filter_dict:
        if filter_func[prop](mol) <= filter_dict[prop]:
            pass
        else:
            filt_mol = None
            break
    return filt_mol


def fp_rdkit(mols, key, nBits):
    for i in mols:
        fp = Chem.RDKFingerprint(mols[i][key], fpSize=nBits)
        mols[i]["fp"] = fp
        #mols[i]["props"]["Fingerprint"] = "RDKit"


def fp_morgan(mols, key, radius, nBits, features, chiral):
    for i in mols:
        fp = Chem.GetMorganFingerprintAsBitVect(mols[i][key], radius, nBits=nBits, useFeatures=features,
                                                useChirality=chiral)
        mols[i]["fp"] = fp


def fp_atompairs(mols, key, nBits, chiral):
    for i in mols:
        fp = Pairs.GetHashedAtomPairFingerprint(mols[i][key], nBits=nBits, includeChirality=chiral)
        mols[i]["fp"] = fp


def fp_maccs(mols, key):
    for i in mols:
        fp = MACCSkeys.GenMACCSKeys(mols[i][key])
        mols[i]["fp"] = fp


def fp_torsion(mols, key, nBits, chiral):
    for i in mols:
        fp = Torsions.GetHashedTopologicalTorsionFingerprint(mols[i][key], nBits=nBits, includeChirality=chiral)
        mols[i]["fp"] = fp


def sim(mols, query, key, cutoff, similarity, filter_dict, name):
    results = {}
    counter = 0
    simis = []
    for j in query:
        for i in mols:
    # for i in mols:
    #     for j in query:
            sim = sim_dict[similarity](mols[i]["fp"], query[j]["fp"])
            if sim >= cutoff:
                filt_mol = filter_res(mols[i][key], filter_dict)
                if filt_mol:
                    simis.append((sim, i ,j))
    print(len(simis))
    grouped_simis = [sorted(list(group), reverse=True) for k, group in groupby(sorted(simis, key=lambda entry: entry[2]), lambda x: x[2])]
    print("---")
    for entry in grouped_simis:
        for element in entry:
            props = copy.deepcopy(mols[element[1]]["props"])
            props["QuerySmiles"] = query[element[2]]["pattern"]
            props["Fingerprint"] = name
            props[sim_dict_name[similarity]] = round(element[0], 5)
            #mols[element[1]]["props"]["QuerySmiles"] = query[element[2]]["pattern"]
            #mols[element[1]]["props"]["Fingerprint"] = name
            #mols[element[1]]["props"][sim_dict_name[similarity]] = round(element[0], 5)
            results[counter] = {"mol": mols[element[1]][key], "props": props, "q_num": element[2]}
            counter += 1
                    # props = copy.deepcopy(mols[i]["props"])
                    # props["QuerySmiles"] = query[j]["pattern"]
                    # props["Fingerprint"] = name
                    # props[sim_dict_name[similarity]] = round(sim, 5)
                    # # mols[i]["props"]["QuerySmiles"] = query[j]["pattern"]
                    # # mols[i]["props"]["Fingerprint"] = name
                    # # mols[i]["props"][sim_dict_name[similarity]] = round(sim, 5)
                    # results[counter] = {"mol": mols[i][key], "props": props, "q_num": j}
                    # counter += 1
    return results


def sim_tver(mols, query, key, cutoff, similarity, filter_dict, name, tva, tvb):
    results = {}
    counter = 0
    simis = []
    for j in query:
        for i in mols:
    # for i in mols:
    #     for j in query:
            sim = sim_dict[similarity](mols[i]["fp"], query[j]["fp"], tva, tvb)
            if sim >= cutoff:
                filt_mol = filter_res(mols[i][key], filter_dict)
                if filt_mol:
                    simis.append((sim, i, j))
    print(len(simis))
    grouped_simis = [sorted(list(group), reverse=True) for k, group in
                     groupby(sorted(simis, key=lambda entry: entry[2]), lambda x: x[2])]
    print("---")
    for entry in grouped_simis:
        for element in entry:
            props = copy.deepcopy(mols[element[1]]["props"])
            props["QuerySmiles"] = query[element[2]]["pattern"]
            props["Fingerprint"] = name
            props[sim_dict_name[similarity]] = round(element[0], 5)
            # mols[element[1]]["props"]["QuerySmiles"] = query[element[2]]["pattern"]
            # mols[element[1]]["props"]["Fingerprint"] = name
            # mols[element[1]]["props"][sim_dict_name[similarity]] = round(element[0], 5)
            results[counter] = {"mol": mols[element[1]][key], "props": props, "q_num": element[2]}
            counter += 1
                    # props = copy.deepcopy(mols[i]["props"])
                    # props["QuerySmiles"] = query[j]["pattern"]
                    # props["Fingerprint"] = name
                    # props[sim_dict_name[similarity]] = round(sim, 5)
                    # # mols[i]["props"]["QuerySmiles"] = query[j]["pattern"]
                    # # mols[i]["props"]["Fingerprint"] = name
                    # # mols[i]["props"][sim_dict_name[similarity]] = round(sim, 5)
                    # results[counter] = {"mol": mols[i][key], "props": props, "q_num": j}
                    # counter += 1
    return results


def sim_top(mols, query, key, top_hits, similarity, filter_dict, name):
    results = {}
    simis = []
    for i in mols:
        filt_mol = filter_res(mols[i][key], filter_dict)
        if filt_mol:
            for j in query:
                sim = sim_dict[similarity](mols[i]["fp"], query[j]["fp"])
                simis.append((sim, i, j))
    print(len(simis))
    grouped_simis = [list(group) for k, group in groupby(sorted(simis, key=lambda entry: entry[2]), lambda x: x[2])]
    print(len(grouped_simis))
    counter = 0
    for entry in grouped_simis:
        nearest_hits = sorted(entry, reverse=True)[:top_hits]
        print(nearest_hits)
        for element in nearest_hits:
            props = copy.deepcopy(mols[element[1]]["props"])
            props["QuerySmiles"] = query[element[2]]["pattern"]
            props["Fingerprint"] = name
            props[sim_dict_name[similarity]] = round(element[0], 5)
            #mols[element[1]]["props"]["QuerySmiles"] = query[element[2]]["pattern"]
            #mols[element[1]]["props"]["Fingerprint"] = name
            #mols[element[1]]["props"][sim_dict_name[similarity]] = round(element[0], 5)
            results[counter] = {"mol": mols[element[1]][key], "props": props, "q_num": element[2]}
            counter += 1
    print(results)
    return results


def sim_top_tver(mols, query, key, top_hits, similarity, filter_dict, name, tva, tvb):
    results = {}
    simis = []
    for i in mols:
        filt_mol = filter_res(mols[i][key], filter_dict)
        if filt_mol:
            for j in query:
                sim = sim_dict[similarity](mols[i]["fp"], query[j]["fp"], tva, tvb)
                simis.append((sim, i, j))
    grouped_simis = [list(group) for k, group in groupby(sorted(simis, key=lambda entry: entry[2]), lambda x: x[2])]
    counter = 0
    for entry in grouped_simis:
        nearest_hits = sorted(entry, reverse=True)[:top_hits]
        for element in nearest_hits:
            props = copy.deepcopy(mols[element[1]]["props"])
            props["QuerySmiles"] = query[element[2]]["pattern"]
            props["Fingerprint"] = name
            props[sim_dict_name[similarity]] = round(element[0], 5)
            # mols[element[1]]["props"]["QuerySmiles"] = query[element[2]]["pattern"]
            # mols[element[1]]["props"]["Fingerprint"] = name
            # mols[element[1]]["props"][sim_dict_name[similarity]] = round(element[0], 5)
            results[counter] = {"mol": mols[element[1]][key], "props": props, "q_num": element[2]}
            counter += 1
    return results


def set_fp_mp(fps, mols):
    for entry in fps:
        mols[entry[0]]["fp"] = entry[1]


def fp_rdkit_mp(mol, i, nBits):
    fp = Chem.RDKFingerprint(mol, fpSize=nBits)
    return (i, fp)


def fp_morgan_mp(mol, i, radius, features, chiral, nBits):
    fp = Chem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, useFeatures=features,
                                            useChirality=chiral)
    return (i, fp)


def fp_atompairs_mp(mol, i, nBits, chiral):
    fp = Pairs.GetHashedAtomPairFingerprint(mol, nBits=nBits, includeChirality=chiral)
    return (i, fp)


def fp_torsion_mp(mol, i, nBits, chiral):
    fp = Torsions.GetHashedTopologicalTorsionFingerprint(mol, nBits=nBits, includeChirality=chiral)
    return (i, fp)


def fp_maccs_mp(mol, i):
    fp = MACCSkeys.GenMACCSKeys(mol)
    return (i, fp)
