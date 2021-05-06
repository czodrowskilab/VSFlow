import copy

from rdkit.Chem import AllChem as Chem

import utils


def set_attrs_mp(pool_results, mols, key, query, results):
    counter = 0
    sort_res = sorted([entry for entry in pool_results if entry], key=lambda entry: entry[1])
    for entry in sort_res:
        props = copy.deepcopy(mols[entry[0]]["props"])
        props["QuerySmiles"] = query[entry[1]]["pattern"]
        results[counter] = {"mol": mols[entry[0]][key], "props": props, "match": entry[2], "q_num": entry[1],
                            "num": entry[0]}
        counter += 1


def substruct_mult(mol, i, query_mol, j, filter_dict):
    match = mol.GetSubstructMatches(query_mol)
    if match:
        filt_mol = utils.filter_res(mol, filter_dict)
        if filt_mol:
            return (i, j, list(match))


def substruct_mult_fm(mol, i, query_mol, j, filter_dict):
    match = mol.GetSubstructMatch(query_mol)
    if match:
        if mol.GetNumHeavyAtoms() == len(match):
            filt_mol = utils.filter_res(mol, filter_dict)
            if filt_mol:
                return (i, j, list(match))


def substruct_mult_fm_nost(mol, i, query_mol, j, filter_dict):
    match = mol.GetSubstructMatch(query_mol)
    if match:
        frags = Chem.GetMolFrags(mol, asMols=True)
        frag_max = max(frags, key=lambda m: m.GetNumAtoms())
        if frag_max.GetNumHeavyAtoms() == len(match):
            filt_mol = utils.filter_res(mol, filter_dict)
            if filt_mol:
                return (i, j, list(match))


def sss_fm_mp(query, mols, key, filter_dict, results, pool):
    argslist = [(mols[i][key], i, query[j]["mol"], j, filter_dict) for i in mols for j in query]
    pool_results = pool.starmap(substruct_mult_fm, argslist)
    set_attrs_mp(pool_results, mols, key, query, results)


def sss_fm_nost_mp(query, mols, key, filter_dict, results, pool):
    argslist = [(mols[i][key], i, query[j]["mol"], j, filter_dict) for i in mols for j in query]
    pool_results = pool.starmap(substruct_mult_fm_nost, argslist)
    set_attrs_mp(pool_results, mols, key, query, results)


def sss_fm_taut_mp(query, mols, key, filter_dict, results, pool):
    argslist = [(mols[i][key], i, taut, j, filter_dict) for i in mols for j in query for taut in query[j]["tauts"]]
    pool_results = pool.starmap(substruct_mult_fm, argslist)
    set_attrs_mp(pool_results, mols, key, query, results)


def sss_mp(query, mols, key, filter_dict, results, pool):
    argslist = [(mols[i][key], i, query[j]["mol"], j, filter_dict) for i in mols for j in query]
    pool_results = pool.starmap(substruct_mult, argslist)
    set_attrs_mp(pool_results, mols, key, query, results)


def sss_mp_taut(query, mols, key, filter_dict, results, pool):
    argslist = [(mols[i][key], i, taut, j, filter_dict) for i in mols for j in query for taut in query[j]["tauts"]]
    pool_results = pool.starmap(substruct_mult, argslist)
    set_attrs_mp(pool_results, mols, key, query, results)


def sss(query, mols, key, filter_dict, results):
    counter = 0
    for j in query:
        for i in mols:
            mol = mols[i][key]
            match = mol.GetSubstructMatches(query[j]["mol"])
            if match:
                filt_mol = utils.filter_res(mol, filter_dict)
                if filt_mol:
                    props = copy.deepcopy(mols[i]["props"])
                    props["QuerySmiles"] = query[j]["pattern"]
                    results[counter] = {"mol": mol, "props": props, "match": list(match), "q_num": j, "num": i}
                    counter += 1


def sss_taut(query, mols, key, filter_dict, results):
    counter = 0
    for i in mols:
        mol = mols[i][key]
        for j in query:
            for taut in query[j]["tauts"]:
                match = mol.GetSubstructMatches(taut)
                if match:
                    filt_mol = utils.filter_res(mol, filter_dict)
                    if filt_mol:
                        props = copy.deepcopy(mols[i]["props"])
                        props["QuerySmiles"] = query[j]["pattern"]
                        results[counter] = {"mol": mol, "props": props, "match": list(match), "q_num": j, "num": i}
                        counter += 1
                        break


def sss_fm(query, mols, key, filter_dict, results):
    counter = 0
    for i in mols:
        mol = mols[i][key]
        for j in query:
            match = mol.GetSubstructMatch(query[j]["mol"])
            if match:
                if mol.GetNumHeavyAtoms() == len(match):
                    filt_mol = utils.filter_res(mol, filter_dict)
                    if filt_mol:
                        props = copy.deepcopy(mols[i]["props"])
                        props["QuerySmiles"] = query[j]["pattern"]
                        results[counter] = {"mol": mol, "props": props, "match": list(match), "q_num": j, "num": i}
                        counter += 1


def sss_fm_nost(query, mols, key, filter_dict, results):
    counter = 0
    for i in mols:
        mol = mols[i][key]
        for j in query:
            match = mol.GetSubstructMatch(query[j]["mol"])
            if match:
                frags = Chem.GetMolFrags(mol, asMols=True)
                frag_max = max(frags, key=lambda m: m.GetNumAtoms())
                if frag_max.GetNumHeavyAtoms() == len(match):
                    filt_mol = utils.filter_res(mol, filter_dict)
                    if filt_mol:
                        props = copy.deepcopy(mols[i]["props"])
                        props["QuerySmiles"] = query[j]["pattern"]
                        results[counter] = {"mol": mol, "props": props, "match": list(match), "q_num": j, "num": i}
                        counter += 1


def sss_fm_taut(query, mols, key, filter_dict, results):
    counter = 0
    for i in mols:
        mol = mols[i][key]
        for j in query:
            for taut in query[j]["tauts"]:
                match = mol.GetSubstructMatch(taut)
                if match:
                    if mol.GetNumHeavyAtoms() == len(match):
                        filt_mol = utils.filter_res(mol, filter_dict)
                        if filt_mol:
                            props = copy.deepcopy(mols[i]["props"])
                            props["QuerySmiles"] = query[j]["pattern"]
                            results[counter] = {"mol": mol, "props": props, "match": list(match), "q_num": j, "num": i}
                            counter += 1
                            break
