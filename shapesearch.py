from rdkit import DataStructs
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate

factory = Gobbi_Pharm2D.factory


### generate query conformations (and 3D pharmacophore fps) for 2D input
def gen_query_conf_pfp(query: dict, num_confs: int, seed: int, max_confs, nthreads):
    factory = Gobbi_Pharm2D.factory
    params = Chem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    # params.pruneRmsThresh = -1.0
    params.numThreads = nthreads
    params.randomSeed = seed
    for i in query:
        mol_H = Chem.AddHs(query[i]["mol"], addCoords=True)
        mol_conf = Chem.Mol(mol_H)
        try:
            if mol_H.GetConformer().Is3D():
                pass
            else:
                mol_conf.RemoveAllConformers()
                print("Generating 3D conformer(s) for query molecule(s)")
                Chem.EmbedMultipleConfs(mol_H, numConfs=num_confs, params=params)
                #Chem.EmbedMultipleConfs(mol_H, numConfs=num_confs, randomSeed=seed, ETversion=2, numThreads=0)
                try:
                    opts = Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=nthreads, maxIters=2000)
                    energies_list = [(opts[j][1], j) for j in range(len(opts))]
                    if len(opts) >= max_confs:
                        sort_energies = sorted(energies_list)[:max_confs]
                    else:
                        sort_energies = sorted(energies_list)[:len(opts)]
                    for entry in sort_energies:
                        mol_conf.AddConformer(mol_H.GetConformer(entry[1]), assignId=True)
                except:
                    if mol_H.GetNumConformers() >= max_confs:
                        for k in range(max_confs):
                            mol_conf.AddConformer(mol_H.GetConformer(k), assignId=True)
                    else:
                        mol_conf = mol_H
        except ValueError:
            mol_conf.RemoveAllConformers()
            print("Generating 3D conformer(s) for query molecule(s)")
            Chem.EmbedMultipleConfs(mol_H, numConfs=num_confs, params=params)
            try:
                opts = Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=nthreads, maxIters=2000)
                energies_list = [(opts[j][1], j) for j in range(len(opts))]
                if len(opts) >= max_confs:
                    sort_energies = sorted(energies_list)[:max_confs]
                else:
                    sort_energies = sorted(energies_list)[:len(opts)]
                for entry in sort_energies:
                    mol_conf.AddConformer(mol_H.GetConformer(entry[1]), assignId=True)
            except:
                if mol_H.GetNumConformers() >= max_confs:
                    for k in range(max_confs):
                        mol_conf.AddConformer(mol_H.GetConformer(k), assignId=True)
                else:
                    mol_conf = mol_H
        query[i]["param"] = Chem.MMFFGetMoleculeProperties(mol_conf)
        fp_mol = []
        for k in range(mol_conf.GetNumConformers()):
            fp = Generate.Gen2DFingerprint(mol_conf, factory, dMat=Chem.Get3DDistanceMatrix(mol_conf, confId=k))
            fp_mol.append(fp)
        query[i]["fp_shape"] = fp_mol
        query[i]["confs"] = mol_conf


def gen_query_conf_pfp_mp(mol, i, num_confs, seed, max_confs, nthreads):
    factory = Gobbi_Pharm2D.factory
    params = Chem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    #params.pruneRmsThresh = -1.0
    params.numThreads = nthreads
    params.randomSeed = seed
    mol_H = Chem.AddHs(mol)
    mol_conf = Chem.Mol(mol_H)
    mol_conf.RemoveAllConformers()
    Chem.EmbedMultipleConfs(mol_H, numConfs=num_confs, params=params)
    try:
        opts = Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=nthreads, maxIters=2000)
        energies_list = [(opts[j][1], j) for j in range(len(opts))]
        if len(opts) >= max_confs:
            sort_energies = sorted(energies_list)[:max_confs]
        else:
            sort_energies = sorted(energies_list)[:len(opts)]
        for entry in sort_energies:
            mol_conf.AddConformer(mol_H.GetConformer(entry[1]), assignId=True)
    except:
        if mol_H.GetNumConformers() >= max_confs:
            for k in range(max_confs):
                mol_conf.AddConformer(mol_H.GetConformer(k), assignId=True)
        else:
            mol_conf = mol_H
    fp_mol = []
    for k in range(mol_conf.GetNumConformers()):
        fp = Generate.Gen2DFingerprint(mol_conf, factory, dMat=Chem.Get3DDistanceMatrix(mol_conf, confId=k))
        fp_mol.append(fp)
    #fp_mol = Generate.Gen2DFingerprint(mol_H, factory, dMat=Chem.Get3DDistanceMatrix(mol_H, confId=0))
    return (i, mol_conf, fp_mol)


def gen_query_pfp_mp(mol, i):
    mol_H = Chem.AddHs(mol, addCoords=True)
    fp_mol = Generate.Gen2DFingerprint(mol_H, factory, dMat=Chem.Get3DDistanceMatrix(mol_H, confId=0))
    fp_mol = [fp_mol]
    return (i, mol_H, fp_mol)


### perform shape screening depending on selected parameters


def shape_pfp_tani(mols, query, nthreads):
    counter = 0
    factory = Gobbi_Pharm2D.factory
    score = []
    for j in mols:
        db_mol_params = Chem.MMFFGetMoleculeProperties(mols[j]["confs"])
        for i in query:
            for confid in range(query[i]["confs"].GetNumConformers()):
                algns = rdMolAlign.GetO3AForProbeConfs(mols[j]["confs"], query[i]["confs"],
                                                       prbPyMMFFMolProperties=db_mol_params,
                                                       refPyMMFFMolProperties=query[i]["param"], refCid=confid,
                                                       numThreads=nthreads)
                shape_simis = []
                for k in range(len(algns)):
                    algns[k].Align()
                    shape_simis.append((1 - Chem.ShapeTanimotoDist(mols[j]["confs"], query[i]["confs"], confId1=k,
                                                                   confId2=confid), k))
                max_shape_sim = max(shape_simis)[0]
                fp_db = Generate.Gen2DFingerprint(mols[j]["confs"], factory,
                                                  dMat=Chem.Get3DDistanceMatrix(mols[j]["confs"],
                                                                                confId=max(shape_simis)[1]))
                pfp_sim = DataStructs.TanimotoSimilarity(query[i]["fp_shape"][confid], fp_db)
                combo = max_shape_sim + pfp_sim
                score.append((combo, max_shape_sim, pfp_sim, i, j, max(shape_simis)[1], mols[j]["confs"]))
                print(counter)
                counter += 1
    return score


def shape_tani_mp(db_mol, i, query_mol, j, query_fp, confid, nthreads):
    fp_q = query_fp[confid]
    algns = rdMolAlign.GetO3AForProbeConfs(db_mol, query_mol,
                                           prbPyMMFFMolProperties=Chem.MMFFGetMoleculeProperties(db_mol),
                                           refPyMMFFMolProperties=Chem.MMFFGetMoleculeProperties(query_mol),
                                           refCid=confid, numThreads=nthreads)
    shape_simis = []
    for k in range(len(algns)):
        algns[k].Align()
        shape_simis.append((1 - Chem.ShapeTanimotoDist(db_mol, query_mol, confId1=k, confId2=confid), k))
    max_shape_sim = max(shape_simis)[0]
    fp_db = Generate.Gen2DFingerprint(db_mol, factory,
                                          dMat=Chem.Get3DDistanceMatrix(db_mol, confId=max(shape_simis)[1]))
    pfp_sim = DataStructs.TanimotoSimilarity(fp_q, fp_db)
    combo = max_shape_sim + pfp_sim
    return (combo, max_shape_sim, pfp_sim, j, i, max(shape_simis)[1], db_mol, confid)

