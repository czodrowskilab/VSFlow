from rdkit import DataStructs
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate


def gen_query_conf_pfp(query, num_confs, seed):
    factory = Gobbi_Pharm2D.factory
    params = Chem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    # params.pruneRmsThresh = -1.0
    params.numThreads = 0
    params.randomSeed = seed
    for i in query:
        if query[i]["mol"].GetNumConformers() == 0:
            mol_H = Chem.AddHs(query[i]["mol"])
            Chem.EmbedMultipleConfs(mol_H, numConfs=num_confs, params=params)
            #Chem.EmbedMultipleConfs(mol_H, numConfs=num_confs, randomSeed=seed, ETversion=2, numThreads=0)
            try:
                Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=0, maxIters=2000)
            except:
                pass
            mol_3D = Chem.RemoveHs(mol_H)
            query[i]["confs"] = mol_3D
            query[i]["param"] = Chem.MMFFGetMoleculeProperties(mol_3D)
            fp_mol = Generate.Gen2DFingerprint(mol_3D, factory, dMat=Chem.Get3DDistanceMatrix(mol_3D, confId=0))
            query[i]["fp_shape"] = fp_mol
        else:
            mol_3D = query[i]["mol"]
            query[i]["confs"] = mol_3D
            query[i]["param"] = Chem.MMFFGetMoleculeProperties(mol_3D)
            fp_mol = Generate.Gen2DFingerprint(mol_3D, factory, dMat=Chem.Get3DDistanceMatrix(mol_3D, confId=0))
            query[i]["fp_shape"] = fp_mol



def shape_pfp_tani(mols, query):
    counter = 0
    factory = Gobbi_Pharm2D.factory
    score = []
    for j in mols:
        db_mol_params = Chem.MMFFGetMoleculeProperties(mols[j]["confs"])
        for i in query:
            algns = rdMolAlign.GetO3AForProbeConfs(mols[j]["confs"], query[i]["confs"],
                                                   prbPyMMFFMolProperties=db_mol_params,
                                                   refPyMMFFMolProperties=query[i]["param"], refCid=0, numThreads=0)
            shape_simis = []
            pharm_simis = []
            for k in range(len(algns)):
                algns[k].Align()
                shape_simis.append(
                    (1 - Chem.ShapeTanimotoDist(mols[j]["confs"], query[i]["confs"], confId1=k, confId2=0), k))
                fp_db = Generate.Gen2DFingerprint(mols[j]["confs"], factory,
                                                  dMat=Chem.Get3DDistanceMatrix(mols[j]["confs"], confId=k))
                pharm_simis.append((DataStructs.TanimotoSimilarity(query[i]["fp_shape"], fp_db), k))
            max_shape_simis = max(shape_simis)[0]
            max_pharm_simis = max(pharm_simis)[0]
            combo = max(shape_simis)[0] + max(pharm_simis)[0]
            score.append((combo, max_shape_simis, max_pharm_simis, i, j, max(shape_simis)[1], mols[j]["confs"]))
            print(counter)
            counter += 1
    print(max(score))
    aligns = sorted(score, reverse=True)
    return aligns


def shape_pfp_tani_mp(db_mol, i, query_mol, query_fp, j):
    factory = Gobbi_Pharm2D.factory
    algns = rdMolAlign.GetO3AForProbeConfs(db_mol, query_mol,
                                           prbPyMMFFMolProperties=Chem.MMFFGetMoleculeProperties(db_mol),
                                           refPyMMFFMolProperties=Chem.MMFFGetMoleculeProperties(query_mol),
                                           refCid=0, numThreads=0)
    shape_simis = []
    pharm_simis = []
    for k in range(len(algns)):
        algns[k].Align()
        shape_simis.append((1 - Chem.ShapeTanimotoDist(db_mol, query_mol, confId1=k, confId2=0), k))
        fp_db = Generate.Gen2DFingerprint(db_mol, factory,
                                          dMat=Chem.Get3DDistanceMatrix(db_mol, confId=k))
        pharm_simis.append((DataStructs.TanimotoSimilarity(query_fp, fp_db), k))
    max_shape_simis = max(shape_simis)[0]
    max_pharm_simis = max(pharm_simis)[0]
    combo = max(shape_simis)[0] + max(pharm_simis)[0]
    return (combo, max_shape_simis, max_pharm_simis, j, i, max(shape_simis)[1], db_mol)


def gen_query_conf_mp(mol, i, num_confs, seed):
    factory = Gobbi_Pharm2D.factory
    if mol.GetNumConformers() == 0:
        params = Chem.ETKDGv3()
        params.useSmallRingTorsions = True
        params.useMacrocycleTorsions = True
        #params.pruneRmsThresh = -1.0
        params.numThreads = 0
        params.randomSeed = seed
        mol_H = Chem.AddHs(mol)
        Chem.EmbedMultipleConfs(mol_H, numConfs=num_confs, params=params)
        # Chem.EmbedMultipleConfs(mol_H, numConfs=num_confs, randomSeed=seed, ETversion=2, numThreads=0)
        try:
            Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=0, maxIters=2000)
        except:
            pass
        mol_3D = Chem.RemoveHs(mol_H)
        fp_mol = Generate.Gen2DFingerprint(mol_3D, factory, dMat=Chem.Get3DDistanceMatrix(mol_3D, confId=0))
    else:
        print("No confs")
        mol_3D = mol
        fp_mol = Generate.Gen2DFingerprint(mol_3D, factory, dMat=Chem.Get3DDistanceMatrix(mol_3D, confId=0))
    #query[i]["confs"] = mol_3D
    # query[i]["param"] = Chem.MMFFGetMoleculeProperties(mol_3D)
    # mol_params = Chem.MMFFGetMoleculeProperties(mol_3D)
    # print(type(mol_params))
    # fp_mol = Generate.Gen2DFingerprint(mol_3D, factory, dMat=Chem.Get3DDistanceMatrix(mol_3D, confId=0))
    return (i, mol_3D, fp_mol)#mol_params)#, fp_mol
