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


def fp_rdkit_taut(query, nBits):
    for i in query:
        for j in range(len(query[i]["tauts"])):
            fp = Chem.RDKFingerprint(query[i]["tauts"][j], fpSize=nBits)
            query[i][f"fp{j}"] = fp


def fp_morgan(mols, key, radius, nBits, features, chiral):
    for i in mols:
        fp = Chem.GetMorganFingerprintAsBitVect(mols[i][key], radius, nBits=nBits, useFeatures=features,
                                                useChirality=chiral)
        mols[i]["fp"] = fp


def fp_morgan_taut(query, radius, nBits, features, chiral):
    for i in query:
        for j in range(len(query[i]["tauts"])):
            fp = Chem.GetMorganFingerprintAsBitVect(query[i]["tauts"][j], radius, nBits=nBits, useFeatures=features,
                                                    useChirality=chiral)
            query[i][f"fp{j}"] = fp


def fp_atompairs(mols, key, nBits, chiral):
    for i in mols:
        fp = Pairs.GetHashedAtomPairFingerprint(mols[i][key], nBits=nBits, includeChirality=chiral)
        mols[i]["fp"] = fp


def fp_atompairs_taut(query, nBits, chiral):
    for i in query:
        for j in range(len(query[i]["tauts"])):
            fp = Pairs.GetHashedAtomPairFingerprint(query[i]["tauts"][j], nBits=nBits, includeChirality=chiral)
            query[i][f"fp{j}"] = fp


def fp_maccs(mols, key):
    for i in mols:
        fp = MACCSkeys.GenMACCSKeys(mols[i][key])
        mols[i]["fp"] = fp


def fp_maccs_taut(query):
    for i in query:
        for j in range(len(query[i]["tauts"])):
            fp = MACCSkeys.GenMACCSKeys(query[i]["tauts"][j])
            query[i][f"fp{j}"] = fp


def fp_torsion(mols, key, nBits, chiral):
    for i in mols:
        fp = Torsions.GetHashedTopologicalTorsionFingerprint(mols[i][key], nBits=nBits, includeChirality=chiral)
        mols[i]["fp"] = fp


def fp_torsion_taut(query, nBits, chiral):
    for i in query:
        for j in range(len(query[i]["tauts"])):
            fp = Torsions.GetHashedTopologicalTorsionFingerprint(query[i]["tauts"][j], nBits=nBits, includeChirality=chiral)
            query[i][f"fp{j}"] = fp


def sim(mols, query, key, fp_key, cutoff, similarity, filter_dict, name, mode):
    results = {}
    counter = 0
    simis = []
    if mode == "all_tauts":
        for j in query:
            for i in mols:
                int_simis = []
                filt_mol = filter_res(mols[i][key], filter_dict)
                if filt_mol:
                    for k in range(len(query[j]["tauts"])):
                        sim = sim_dict[similarity](mols[i][fp_key], query[j][f"fp{k}"])
                        if sim >= cutoff:
                            int_simis.append(sim)
                if int_simis:
                    simis.append((max(int_simis), i, j))
    else:
        for j in query:
            for i in mols:
                sim = sim_dict[similarity](mols[i][fp_key], query[j]["fp"])
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


def sim_tver(mols, query, key, fp_key, cutoff, similarity, filter_dict, name, tva, tvb, mode):
    results = {}
    counter = 0
    simis = []
    if mode == "all_tauts":
        for j in query:
            for i in mols:
                int_simis = []
                filt_mol = filter_res(mols[i][key], filter_dict)
                if filt_mol:
                    for k in range(len(query[j]["tauts"])):
                        sim = sim_dict[similarity](mols[i][fp_key], query[j][f"fp{k}"])
                        if sim >= cutoff:
                            int_simis.append(sim)
                if int_simis:
                    simis.append((max(int_simis), i, j))
    else:
        for j in query:
            for i in mols:
                sim = sim_dict[similarity](mols[i][fp_key], query[j]["fp"], tva, tvb)
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


def sim_top(mols, query, key, fp_key, top_hits, similarity, filter_dict, name, mode):
    results = {}
    simis = []
    if mode == "all_tauts":
        for j in query:
            for i in mols:
                int_simis = []
                filt_mol = filter_res(mols[i][key], filter_dict)
                if filt_mol:
                    for k in range(len(query[j]["tauts"])):
                        sim = sim_dict[similarity](mols[i][fp_key], query[j][f"fp{k}"])
                        int_simis.append(sim)
                if int_simis:
                    print(len(int_simis))
                    simis.append((max(int_simis), i, j))
    else:
        for i in mols:
            filt_mol = filter_res(mols[i][key], filter_dict)
            if filt_mol:
                for j in query:
                    sim = sim_dict[similarity](mols[i][fp_key], query[j]["fp"])
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


def sim_top_tver(mols, query, key, fp_key, top_hits, similarity, filter_dict, name, tva, tvb, mode):
    results = {}
    simis = []
    if mode == "all_tauts":
        for j in query:
            for i in mols:
                int_simis = []
                filt_mol = filter_res(mols[i][key], filter_dict)
                if filt_mol:
                    for k in range(len(query[j]["tauts"])):
                        sim = sim_dict[similarity](mols[i][fp_key], query[j][f"fp{k}"])
                        int_simis.append(sim)
                if int_simis:
                    simis.append((max(int_simis), i, j))
    else:
        for i in mols:
            filt_mol = filter_res(mols[i][key], filter_dict)
            if filt_mol:
                for j in query:
                    sim = sim_dict[similarity](mols[i][fp_key], query[j]["fp"], tva, tvb)
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


def set_fp_taut_mp(fps, mols):
    for entry in fps:
        mols[entry[0]][f"fp{entry[2]}"] = entry[1]


def fp_rdkit_mp(mol, i, nBits):
    fp = Chem.RDKFingerprint(mol, fpSize=nBits)
    return (i, fp)


def fp_rdkit_taut_mp(taut, i, k, nBits):
    fp = Chem.RDKFingerprint(taut, fpSize=nBits)
    return (i, fp, k)


def fp_morgan_mp(mol, i, radius, features, chiral, nBits):
    fp = Chem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, useFeatures=features,
                                            useChirality=chiral)
    return (i, fp)


def fp_morgan_taut_mp(taut, i, k, radius, features, chiral, nBits):
    fp = Chem.GetMorganFingerprintAsBitVect(taut, radius, nBits=nBits, useFeatures=features,
                                            useChirality=chiral)
    return (i, fp, k)


def fp_atompairs_mp(mol, i, nBits, chiral):
    fp = Pairs.GetHashedAtomPairFingerprint(mol, nBits=nBits, includeChirality=chiral)
    return (i, fp)


def fp_atompairs_taut_mp(taut, i, k, nBits, chiral):
    fp = Pairs.GetHashedAtomPairFingerprint(taut, nBits=nBits, includeChirality=chiral)
    return (i, fp, k)


def fp_torsion_mp(mol, i, nBits, chiral):
    fp = Torsions.GetHashedTopologicalTorsionFingerprint(mol, nBits=nBits, includeChirality=chiral)
    return (i, fp)


def fp_torsion_taut_mp(taut, i, k, nBits, chiral):
    fp = Torsions.GetHashedTopologicalTorsionFingerprint(taut, nBits=nBits, includeChirality=chiral)
    return (i, fp, k)


def fp_maccs_mp(mol, i):
    fp = MACCSkeys.GenMACCSKeys(mol)
    return (i, fp)


def fp_maccs_taut_mp(taut, i, k):
    fp = MACCSkeys.GenMACCSKeys(taut)
    return (i, fp, k)
