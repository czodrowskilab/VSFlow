from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors


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


def set_attrs_mp(pool_results, mols, key, query, results):
    counter = 0
    for entry in pool_results:
        if entry:
            mols[entry[0]]["props"]["QuerySmiles"] = query[entry[1]]["pattern"]
            results[counter] = {"mol": mols[entry[0]][key], "props": mols[entry[0]]["props"], "match": entry[2], "q_num": entry[1]}
            counter += 1


def substruct_mult(mol, i, query_mol, j, filter_dict):
    match = mol.GetSubstructMatches(query_mol)
    if match:
        filt_mol = filter_res(mol, filter_dict)
        if filt_mol:
            return (i, j, list(match))


def substruct_mult_fm(mol, i, query_mol, j, filter_dict):
    match = mol.GetSubstructMatches(query_mol)
    if match:
        if mol.GetNumHeavyAtoms() == len(tuple(l for k in match for l in k)):
            filt_mol = filter_res(mol, filter_dict)
            if filt_mol:
                return (i, j, list(match))


def substruct_mult_fm_nost(mol, i, query_mol, j, filter_dict):
    match = mol.GetSubstructMatches(query_mol)
    if match:
        frags = Chem.GetMolFrags(mol, asMols=True)
        frag_max = max(frags, key=lambda m: m.GetNumAtoms())
        if frag_max.GetNumHeavyAtoms() == len(tuple(l for k in match for l in k)):
            filt_mol = filter_res(mol, filter_dict)
            if filt_mol:
                return (i, j, list(match))


def sss_fm_mp(query, mols, key, filter_dict, results, pool):
    #pool = mp.Pool(processes=np)
    argslist = [(mols[i][key], i, query[j]["mol"], j, filter_dict) for i in mols for j in query]
    pool_results = pool.starmap(substruct_mult_fm, argslist)
    set_attrs_mp(pool_results, mols, key, query, results)
    #pool.close()


def sss_fm_nost_mp(query, mols, key, filter_dict, results, pool):
    #pool = mp.Pool(processes=np)
    argslist = [(mols[i][key], i, query[j]["mol"], j, filter_dict) for i in mols for j in query]
    pool_results = pool.starmap(substruct_mult_fm_nost, argslist)
    set_attrs_mp(pool_results, mols, key, query, results)
    #pool.close()


def sss_fm_taut_mp(query, mols, key, filter_dict, results, pool):
    #pool = mp.Pool(processes=np)
    argslist = [(mols[i][key], i, taut, j, filter_dict) for i in mols for j in query for taut in query[j]["tauts"]]
    pool_results = pool.starmap(substruct_mult_fm, argslist)
    set_attrs_mp(pool_results, mols, key, query, results)
    #pool.close()


def sss_mp(query, mols, key, filter_dict, results, pool):
    #pool = mp.Pool(processes=np)
    argslist = [(mols[i][key], i, query[j]["mol"], j, filter_dict) for i in mols for j in query]
    pool_results = pool.starmap(substruct_mult, argslist)
    set_attrs_mp(pool_results, mols, key, query, results)
    #pool.close()


def sss_mp_taut(query, mols, key, filter_dict, results, pool):
    #pool = mp.Pool(processes=np)
    argslist = [(mols[i][key], i, taut, j, filter_dict) for i in mols for j in query for taut in query[j]["tauts"]]
    pool_results = pool.starmap(substruct_mult, argslist)
    set_attrs_mp(pool_results, mols, key, query, results)
    #pool.close()


def sss(query, mols, key, filter_dict, results):
    counter = 0
    for i in mols:
        mol = mols[i][key]
        for j in query:
            match = mol.GetSubstructMatches(query[j]["mol"])
            if match:
                filt_mol = filter_res(mol, filter_dict)
                if filt_mol:
                    #set_attrs(mol, query[j][1], j, match)
                    mols[i]["props"]["QuerySmiles"] = query[j]["pattern"]
                    #mols[i]["props"]["Similarity"] = match
                    #results[i] = (mol, mols[i]["props"])
                    results[counter] = {"mol": mol, "props": mols[i]["props"], "match": list(match), "q_num": j}
                    counter += 1


def sss_taut(query, mols, key, filter_dict, results):
    counter = 0
    for i in mols:
        mol = mols[i][key]
        for j in query:
            for taut in query[j]["tauts"]:
                match = mol.GetSubstructMatches(taut)
                if match:
                    filt_mol = filter_res(mol, filter_dict)
                    if filt_mol:
                        mols[i]["props"]["QuerySmiles"] = query[j]["pattern"]
                        #set_attrs(mol, query[j][1], j, match)
                        results[counter] = {"mol": mol, "props": mols[i]["props"], "match": list(match), "q_num": j}
                        counter += 1
                        break


def sss_fm(query, mols, key, filter_dict, results):
    counter = 0
    for i in mols:
        mol = mols[i][key]
        for j in query:
            match = mol.GetSubstructMatches(query[j]["mol"])
            if match:
                if mol.GetNumHeavyAtoms() == len(tuple(l for k in match for l in k)):
                    filt_mol = filter_res(mol, filter_dict)
                    if filt_mol:
                        mols[i]["props"]["QuerySmiles"] = query[j]["pattern"]
                        #set_attrs(mol, query[j][1], j, match)
                        results[counter] = {"mol": mol, "props": mols[i]["props"], "match": list(match), "q_num": j}


def sss_fm_nost(query, mols, key, filter_dict, results):
    counter = 0
    for i in mols:
        mol = mols[i][key]
        for j in query:
            match = mol.GetSubstructMatches(query[j]["mol"])
            if match:
                frags = Chem.GetMolFrags(mol, asMols=True)
                frag_max = max(frags, key=lambda m: m.GetNumAtoms())
                if frag_max.GetNumHeavyAtoms() == len(tuple(l for k in match for l in k)):
                    filt_mol = filter_res(mol, filter_dict)
                    if filt_mol:
                        mols[i]["props"]["QuerySmiles"] = query[j]["pattern"]
                        #set_attrs(mol, query[j][1], j, match)
                        results[counter] = {"mol": mol, "props": mols[i]["props"], "match": list(match), "q_num": j}


def sss_fm_taut(query, mols, key, filter_dict, results):
    counter = 0
    for i in mols:
        mol = mols[i][key]
        for j in query:
            for taut in query[j]["tauts"]:
                match = mol.GetSubstructMatches(taut)
                if match:
                    if mol.GetNumHeavyAtoms() == len(tuple(l for k in match for l in k)):
                        filt_mol = filter_res(mol, filter_dict)
                        if filt_mol:
                            mols[i]["props"]["QuerySmiles"] = query[j]["pattern"]
                            #set_attrs(mol, query[j][1], j, match)
                            results[counter] = {"mol": mol, "props": mols[i]["props"], "match": list(match), "q_num": j}
                            counter += 1
                            break
