from rdkit import DataStructs
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate, SigFactory
from rdkit.Chem import ChemicalFeatures
import os

script_path = os.path.dirname(os.path.abspath(__file__))

def init_sig(fdefName):
    featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    #sigFactory = SigFactory.SigFactory(featFactory, minPointCount=2, maxPointCount=3, trianglePruneBins=True, includeBondOrder=True, useCounts=True, shortestPathsOnly=False)
    sigFactory = SigFactory.SigFactory(featFactory, minPointCount=2, maxPointCount=3, trianglePruneBins=False)
    sigFactory.SetBins([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 100)])
    #sigFactory.SetBins([(2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 100)])
    #sigFactory.SetBins([(0, 2), (2, 5), (5, 8)])
    sigFactory.Init()
    return sigFactory


shape_dist = {"tan": Chem.ShapeTanimotoDist, "protr": Chem.ShapeProtrudeDist}
mol_align = {"mmff": Chem.GetO3AForProbeConfs, "crippen": Chem.GetCrippenO3AForProbeConfs}
align_contribs = {"mmff": Chem.MMFFGetMoleculeProperties, "crippen": rdMolDescriptors._CalcCrippenContribs}
sim_measures = {"tan": DataStructs.TanimotoSimilarity, "dice": DataStructs.DiceSimilarity,
                "cos": DataStructs.CosineSimilarity,
                "sok": DataStructs.SokalSimilarity, "russ": DataStructs.RusselSimilarity,
                "kulc": DataStructs.KulczynskiSimilarity,
                "mcco": DataStructs.McConnaugheySimilarity}
feat_factory = {"gobbi": Gobbi_Pharm2D.factory, "basic": init_sig(f'{script_path}/resources/BaseFeatures.fdef'),
                "minimal": init_sig(f'{script_path}/resources/MinimalFeatures.fdef')}


def sim(fp1, fp2, simi, tver_a=0.5, tver_b=0.5):

    if simi == "tver":
        pfp_sim = DataStructs.TverskySimilarity(fp1, fp2, tver_a, tver_b)
    else:
        pfp_sim = sim_measures[simi](fp1, fp2)
    return pfp_sim


### generate query conformations (and 3D pharmacophore fps) for 2D input
def gen_query_conf_pfp(query: dict, num_confs: int, seed: int, max_confs, nthreads, align_method, pharm_def):
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
                print(f"Generating 3D conformer(s) for query molecule {i}")
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
        query[i]["param"] = align_contribs[align_method](mol_conf)
        fp_mol = []
        for k in range(mol_conf.GetNumConformers()):
            print("start")
            print(k)
            fp = Generate.Gen2DFingerprint(mol_conf, feat_factory[pharm_def], dMat=Chem.Get3DDistanceMatrix(mol_conf, confId=k))
            fp_mol.append(fp)
        print("end")
        query[i]["fp_shape"] = fp_mol
        query[i]["confs"] = mol_conf


def gen_query_conf_pfp_mp(mol, i, num_confs, seed, max_confs, nthreads):
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
    # fp_mol = []
    # for k in range(mol_conf.GetNumConformers()):
    #     fp = Generate.Gen2DFingerprint(mol_conf, feat_factory[pharm_def], dMat=Chem.Get3DDistanceMatrix(mol_conf, confId=k))
    #     fp_mol.append(fp)
    #fp_mol = Generate.Gen2DFingerprint(mol_H, factory, dMat=Chem.Get3DDistanceMatrix(mol_H, confId=0))
    return (i, mol_conf)


def gen_query_pfp_mp(mol, i):
    mol_H = Chem.AddHs(mol, addCoords=True)
    # print("start1")
    # fp_mol = Generate.Gen2DFingerprint(mol_H, feat_factory[pharm_def], dMat=Chem.Get3DDistanceMatrix(mol_H, confId=-1))
    # print("end1")
    # fp_mol = [fp_mol]
    return (i, mol_H)


### perform shape screening depending on selected parameters


def shape_search(mols, query, nthreads, align_method, dist, fp_simi, pharm_def, tva, tvb):
    counter = 0
    score = []
    for j in mols:
        if "confs" in mols[j]:
            try:
                db_mol_params = align_contribs[align_method](mols[j]["confs"])
            except:
                db_mol_params = None
                pass
            for i in query:
                for confid in range(query[i]["confs"].GetNumConformers()):
                    # algns = Chem.GetO3AForProbeConfs(mols[j]["confs"], query[i]["confs"],
                    #                                        prbPyMMFFMolProperties=db_mol_params,
                    #                                        refPyMMFFMolProperties=query[i]["param"], refCid=confid,
                    #                                        numThreads=nthreads)
                    try:
                        algns = mol_align[align_method](mols[j]["confs"], query[i]["confs"], nthreads,
                                                        db_mol_params,
                                                        query[i]["param"],
                                                        refCid=confid)
                    except:
                        algns = []
                    if algns:
                        shape_simis = []
                        if dist == "tver":
                            for k in range(len(algns)):
                                algns[k].Align()
                                shape_simis.append(
                                    (Chem.ShapeTverskyIndex(mols[j]["confs"], query[i]["confs"], tva, tvb, confId1=k, confId2=confid), k))
                        else:
                            for k in range(len(algns)):
                                algns[k].Align()
                                shape_simis.append((1 - shape_dist[dist](mols[j]["confs"], query[i]["confs"], confId1=k,
                                                                               confId2=confid), k))
                        max_shape_sim = max(shape_simis)[0]
                        fp_db = Generate.Gen2DFingerprint(mols[j]["confs"], feat_factory[pharm_def],
                                                          dMat=Chem.Get3DDistanceMatrix(mols[j]["confs"],
                                                                                        confId=max(shape_simis)[1]))
                        pfp_sim = sim(query[i]["fp_shape"][confid], fp_db, fp_simi, tva, tvb)
                        combo = (max_shape_sim + pfp_sim) / 2
                        score.append((combo, max_shape_sim, pfp_sim, i, j, max(shape_simis)[1], Chem.Mol(mols[j]["confs"]), confid, mols[j]["pattern"], query[i]["pattern"]))
                        print(counter)
                        counter += 1
    return score


def shape_mp(db_mol, i, db_smi, query_mol, j, q_smi, confid, nthreads, dist, fp_simi, tva, tvb, align_method, pharm_def):
    #fp_q = query_fp[confid]
    fp_q = Generate.Gen2DFingerprint(query_mol, feat_factory[pharm_def],
                              dMat=Chem.Get3DDistanceMatrix(query_mol, confId=confid))
    try:
        algns = mol_align[align_method](db_mol, query_mol, nthreads,
                                               align_contribs[align_method](db_mol),
                                               align_contribs[align_method](query_mol),
                                               refCid=confid)
    except:
        algns = []
    if algns:
        shape_simis = []
        if dist == "tver":
            for k in range(len(algns)):
                algns[k].Align()
                shape_simis.append((Chem.ShapeTverskyIndex(db_mol, query_mol, tva, tvb, confId1=k, confId2=confid), k))
        else:
            for k in range(len(algns)):
                algns[k].Align()
                shape_simis.append((1 - shape_dist[dist](db_mol, query_mol, confId1=k, confId2=confid), k))

        max_shape_sim = max(shape_simis)[0]
        # mol = Chem.Mol(db_mol)
        # mol.RemoveAllConformers()
        # mol.AddConformer(db_mol.GetConformer(max(shape_simis)[1]))
        # fp_db = Generate.Gen2DFingerprint(mol, feat_factory[pharm_def],
        #                                       dMat=Chem.Get3DDistanceMatrix(mol, confId=-1))
        fp_db = Generate.Gen2DFingerprint(db_mol, feat_factory[pharm_def],
                                              dMat=Chem.Get3DDistanceMatrix(db_mol, confId=max(shape_simis)[1]))
        pfp_sim = sim(fp_q, fp_db, fp_simi, tva, tvb)
        combo = (max_shape_sim + pfp_sim) / 2
        print(i)
        return (combo, max_shape_sim, pfp_sim, j, i, max(shape_simis)[1], Chem.Mol(db_mol), confid, db_smi, q_smi)
