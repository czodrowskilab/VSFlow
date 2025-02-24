import gzip
import ssl
from itertools import groupby
from urllib.request import urlopen

from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from xlrd import open_workbook
from rdkit.Chem.MolStandardize import rdMolStandardize


def conv_smi(entry):
    return {"mol": Chem.MolFromSmiles(entry[1]), "props": {"chembl_id": entry[0].decode()}}


def req_chembl(nproc, pool):
    gcontext = ssl.SSLContext()
    r = urlopen("https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_28/chembl_28_chemreps.txt.gz",
                context=gcontext)
    fd = gzip.GzipFile(fileobj=r, mode="rb")
    byte_cont = fd.readlines()
    content = {}
    for i in range(1, len(byte_cont)):
        content[i] = byte_cont[i].split(b'\t')[:2]
    del byte_cont
    if nproc:
        sub = {}
        preps = pool.map(conv_smi, [content[n] for n in content])
        for i in range(len(preps)):
            sub[i] = preps[i]
        del preps
    else:
        sub = {}
        for n in content:
            mol = Chem.MolFromSmiles(content[n][1])
            if mol:
                sub[n] = {"mol": mol, "props": {"chembl_id": content[n][0].decode()}}
    del content
    return sub


def req_pdb():
    r = urlopen("http://ligand-expo.rcsb.org/dictionaries/Components-smiles-stereo-oe.smi")
    c = r.readlines()
    cont = {}
    for i in range(len(c)):
        cont[i] = c[i].strip(b"\n").split(b"\t")
    sub = {}
    for n in cont:
        try:
            mol = Chem.MolFromSmiles(cont[n][0])
            if mol:
                sub[n] = {"mol": mol,
                          "props": {"ligand_id": cont[n][1].decode(), "name": cont[n][2].decode()}}
        except:
            continue
    return sub


def query_standardize(mol):
    try:
        mol_sta = rdMolStandardize.ChargeParent(mol)
        return mol_sta
    except:
        return mol


def query_canonicalize(mol, ntauts=1):
    try:
        mol_sta = rdMolStandardize.ChargeParent(mol)
        mol_can = rdMolStandardize.TautomerParent(mol_sta)
        return mol_can
    except:
        return mol


def query_enumerate(mol, ntauts):
    try:
        mol_sta = rdMolStandardize.ChargeParent(mol)
        mol_tauts = [mol for mol in rdMolStandardize.TautomerEnumerator().Enumerate(mol_sta)]
        return (mol_sta, mol_tauts)
    except:
        return (mol, [mol])


def read_smiles(smiles, mode, ntauts):
    sub = {}
    if mode == "std":
        for i in range(len(smiles)):
            mol = Chem.MolFromSmiles(smiles[i])
            if mol:
                mol_std = query_standardize(mol)
                sub[i] = {"mol": mol_std, "pattern": smiles[i]}
    elif mode == "can_taut":
        for i in range(len(smiles)):
            mol = Chem.MolFromSmiles(smiles[i])
            if mol:
                mol_std = query_canonicalize(mol, ntauts)
                sub[i] = {"mol": mol_std, "pattern": smiles[i]}
    elif mode == "all_tauts":
        for i in range(len(smiles)):
            mol = Chem.MolFromSmiles(smiles[i])
            if mol:
                mol_tauts = query_enumerate(mol, ntauts)
                sub[i] = {"mol": mol_tauts[0], "pattern": smiles[i], "tauts": mol_tauts[1]}
    else:
        for i in range(len(smiles)):
            mol = Chem.MolFromSmiles(smiles[i])
            if mol:
                frags = Chem.GetMolFrags(mol, asMols=True)
                frag_max = max(frags, key=lambda m: m.GetNumAtoms())
                # sub[i] = (frag_max, smiles[i])
                sub[i] = {"mol": frag_max, "pattern": smiles[i]}
    return sub


def read_smarts(smarts):
    sub = {}
    for i in range(len(smarts)):
        mol = Chem.MolFromSmarts(smarts[i])
        if mol:
            sub[i] = {"mol": mol, "pattern": smarts[i], "tauts": [mol]}
    return sub


def read_smarts_file(file):
    sub = {}
    with open(file, "r") as sma_file:
        content = sma_file.readlines()
    for i in range(len(content)):
        sma = content[i].strip("\n")
        if sma != "":
            mol = Chem.MolFromSmarts(sma)
            if mol:
                sub[i] = {"mol": mol, "pattern": sma, "tauts": [mol]}
    return sub


def read_sd(infile, mode, ntauts, gz=False):
    sub = {}
    if gz:
        with gzip.open(infile, mode="rt") as inf:
            content = inf.readlines()
    else:
        with open(infile, "r") as sd_file:
            content = sd_file.readlines()
    try:
        sd_blocks = [list(group) for k, group in groupby(content, lambda x: x == "$$$$\n") if not k]
    except ValueError:
        return sub
    del content
    if mode == "std":
        for i in range(len(sd_blocks)):
            mol_block_list = sd_blocks[i][:sd_blocks[i].index("M  END\n") + 1]
            mol_block = ''.join([elem for elem in mol_block_list])
            mol = Chem.MolFromMolBlock(mol_block)
            if mol:
                mol_std = query_standardize(mol)
                sub[i] = {"mol": mol_std, "pattern": Chem.MolToSmiles(mol)}
    elif mode == "can_taut":
        for i in range(len(sd_blocks)):
            mol_block_list = sd_blocks[i][:sd_blocks[i].index("M  END\n") + 1]
            mol_block = ''.join([elem for elem in mol_block_list])
            mol = Chem.MolFromMolBlock(mol_block)
            if mol:
                mol_std = query_canonicalize(mol, ntauts)
                sub[i] = {"mol": mol_std, "pattern": Chem.MolToSmiles(mol)}
    elif mode == "all_tauts":
        for i in range(len(sd_blocks)):
            mol_block_list = sd_blocks[i][:sd_blocks[i].index("M  END\n") + 1]
            mol_block = ''.join([elem for elem in mol_block_list])
            mol = Chem.MolFromMolBlock(mol_block)
            if mol:
                mol_std = query_enumerate(mol, ntauts)
                sub[i] = {"mol": mol_std[0], "pattern": Chem.MolToSmiles(mol), "tauts": mol_std[1]}
    else:
        for i in range(len(sd_blocks)):
            mol_block_list = sd_blocks[i][:sd_blocks[i].index("M  END\n") + 1]
            mol_block = ''.join([elem for elem in mol_block_list])
            mol = Chem.MolFromMolBlock(mol_block)
            if mol:
                sub[i] = {"mol": mol, "pattern": Chem.MolToSmiles(mol)}
    return sub


def read_csv(filename, smiles_column, delimiter, mode="std", ntauts=100, header=None, db=False, gz=False):
    sub = {}
    mol_func = Chem.MolFromSmiles
    if gz:
        with gzip.open(filename, mode="rt") as inf:
            content = inf.readlines()
    else:
        with open(filename, "r") as file:
            content = file.readlines()
    if delimiter is None:
        test_del = [";", ",", "\t"]
        if len(content) > 1:
            for sep in test_del:
                if len(content[0].strip("\n").split(sep)) == len(content[1].strip("\n").split(sep)):
                    if len(content[0].strip("\n").split(sep)) > 1:
                        delimiter = sep
        if delimiter is None:
            counter = 0
            for line in content:
                if line != "":
                    mol = mol_func(line.strip("\n"))
                    if mol:
                        break
                    mol = Chem.MolFromInchi(line.strip("\n"))
                    if mol:
                        mol_func = Chem.MolFromInchi
                        break
                    counter += 1
                if counter >= 100:
                    return sub
            counter = 0
            if mode == "std":
                for line in content:
                    mol = mol_func(line.strip("\n"))
                    if mol:
                        props = {}
                        try:
                            name = mol.GetProp("_Name")
                            props["Name"] = name
                        except KeyError:
                            pass
                        mol_sta = query_standardize(mol)
                        sub[counter] = {"mol": mol_sta, "pattern": line.strip("\n"), "props": props}
                        counter += 1
            elif mode == "can_taut":
                for line in content:
                    mol = mol_func(line.strip("\n"))
                    if mol:
                        props = {}
                        try:
                            name = mol.GetProp("_Name")
                            props["Name"] = name
                        except KeyError:
                            pass
                        mol_can = query_canonicalize(mol, ntauts)
                        sub[counter] = {"mol": mol_can, "pattern": line.strip("\n"), "props": props}
                        counter += 1
            elif mode == "all_tauts":
                for line in content:
                    mol = mol_func(line.strip("\n"))
                    if mol:
                        props = {}
                        try:
                            name = mol.GetProp("_Name")
                            props["Name"] = name
                        except KeyError:
                            pass
                        mol_tauts = query_enumerate(mol, ntauts)
                        sub[counter] = {"mol": mol_tauts[0], "pattern": line.strip("\n"), "tauts": mol_tauts[1], "props": props}
                        counter += 1
            else:
                for line in content:
                    mol = mol_func(line.strip("\n"))
                    if mol:
                        props = {}
                        try:
                            name = mol.GetProp("_Name")
                            props["Name"] = name
                        except KeyError:
                            pass
                        sub[counter] = {"mol": mol, "pattern": line.strip("\n"), "props": props}
                        counter += 1
            return sub
    proc_content = {}
    for i in range(len(content)):
        proc_content[i] = content[i].strip("\n").split(delimiter)
    del content
    if smiles_column is None:
        smi_pos = None
        counter = 0
        for n in proc_content:
            for i in range(len(proc_content[n])):
                if proc_content[n][i] != "":
                    mol = Chem.MolFromSmiles(proc_content[n][i])
                    if mol:
                        smi_pos = i
                        break
                    mol = Chem.MolFromInchi(proc_content[n][i])
                    if mol:
                        mol_func = Chem.MolFromInchi
                        smi_pos = i
                        break
            if smi_pos is not None:
                break
            counter += 1
            if counter >= 100:
                return sub
    else:
        smi_pos = None
        try:
            smi_pos = int(smiles_column)
            counter = 0
            for n in proc_content:
                if proc_content[n][smi_pos] != "":
                    mol = Chem.MolFromSmiles(proc_content[n][smi_pos])
                    if mol:
                        break
                    mol = Chem.MolFromInchi(proc_content[n][smi_pos])
                    if mol:
                        mol_func = Chem.MolFromInchi
                        break
                    counter += 1
                    if counter >= 100:
                        return sub
        except ValueError:
            counter = 0
            for n in proc_content:
                if smiles_column in proc_content[n]:
                    smi_pos = proc_content[n].index(smiles_column)
                    break
                counter += 1
                if counter >= 100:
                    return sub
            counter = 0
            for n in proc_content:
                if proc_content[n][smi_pos] != "":
                    mol = Chem.MolFromSmiles(proc_content[n][smi_pos])
                    if mol:
                        break
                    mol = Chem.MolFromInchi(proc_content[n][smi_pos])
                    if mol:
                        mol_func = Chem.MolFromInchi
                        break
                    counter += 1
                    if counter >= 100:
                        return sub
        if smi_pos is None:
            return sub
    if db:
        if header:
            titles = proc_content[header - 1]
            for n in proc_content:
                try:
                    if proc_content[n][smi_pos] != "":
                        mol = mol_func(proc_content[n][smi_pos])
                        if mol:
                            props = {}
                            for k in range(len(titles)):
                                props[titles[k]] = proc_content[n][k]
                            sub[n] = {"mol": mol, "pattern": proc_content[n][smi_pos], "props": props}
                except IndexError:
                    continue
        else:
            for n in proc_content:
                try:
                    if proc_content[n][smi_pos] != "":
                        mol = mol_func(proc_content[n][smi_pos])
                        if mol:
                            props = {}
                            for i in range(len(proc_content[n])):
                                # props = {name: proc_content[n][smi_pos]}

                                props[i] = proc_content[n][i]
                            sub[n] = {"mol": mol, "pattern": proc_content[n][smi_pos], "props": props}
                except IndexError:
                    continue
    else:
        if smi_pos is not None:
            if mode == "std":
                for n in proc_content:
                    try:
                        if proc_content[n][smi_pos] != "":
                            mol = mol_func(proc_content[n][smi_pos])
                            if mol:
                                props = {}
                                try:
                                    name = mol.GetProp("_Name")
                                    props["Name"] = name
                                except KeyError:
                                    pass
                                mol_std = query_standardize(mol)
                                sub[n] = {"mol": mol_std, "pattern": proc_content[n][smi_pos], "props": props}
                    except IndexError:
                        continue
            elif mode == "can_taut":
                for n in proc_content:
                    try:
                        if proc_content[n][smi_pos] != "":
                            mol = mol_func(proc_content[n][smi_pos])
                            if mol:
                                props = {}
                                try:
                                    name = mol.GetProp("_Name")
                                    props["Name"] = name
                                except KeyError:
                                    pass
                                mol_std = query_canonicalize(mol, ntauts)
                                sub[n] = {"mol": mol_std, "pattern": proc_content[n][smi_pos], "props": props}
                    except IndexError:
                        continue
            elif mode == "all_tauts":
                for n in proc_content:
                    try:
                        if proc_content[n][smi_pos] != "":
                            mol = mol_func(proc_content[n][smi_pos])
                            if mol:
                                props = {}
                                try:
                                    name = mol.GetProp("_Name")
                                    props["Name"] = name
                                except KeyError:
                                    pass
                                mol_std = query_enumerate(mol, ntauts)
                                sub[n] = {"mol": mol_std[0], "pattern": proc_content[n][smi_pos], "tauts": mol_std[1], "props": props}
                    except IndexError:
                        continue
            else:
                for n in proc_content:
                    try:
                        if proc_content[n][smi_pos] != "":
                            mol = mol_func(proc_content[n][smi_pos])
                            if mol:
                                props = {}
                                try:
                                    name = mol.GetProp("_Name")
                                    props["Name"] = name
                                except KeyError:
                                    pass
                                sub[n] = {"mol": mol, "pattern": proc_content[n][smi_pos], "props": props}
                    except IndexError:
                        continue
    return sub


def read_excel(filename, smiles_column, mode="std", ntauts=100, header=None, db=False):
    sub = {}
    wb = open_workbook(filename)
    mol_func = Chem.MolFromSmiles
    if smiles_column is None:
        smi_pos = None
        for sheet in wb.sheets():
            counter = 0
            for i in range(sheet.nrows):
                for j in range(sheet.ncols):
                    smi = sheet.row(i)[j].value
                    if smi != "":
                        mol = Chem.MolFromSmiles(sheet.row(i)[j].value)
                        if mol:
                            smi_pos = j
                            break
                        mol = Chem.MolFromInchi(sheet.row(i)[j].value)
                        if mol:
                            mol_func = Chem.MolFromInchi
                            smi_pos = j
                            break
                if smi_pos is not None:
                    break
                counter += 1
                if counter >= 100:
                    return sub
    else:
        smi_pos = None
        try:
            smi_pos = int(smiles_column)
            for sheet in wb.sheets():
                counter = 0
                for i in range(sheet.nrows):
                    smi = sheet.row(i)[smiles_column].value
                    if smi != "":
                        mol = Chem.MolFromSmiles(smi)
                        if mol:
                            break
                        mol = Chem.MolFromInchi(sheet.row(i)[smiles_column].value)
                        if mol:
                            mol_func = Chem.MolFromInchi
                            break
                if smi_pos is not None:
                    break
                counter += 1
                if counter >= 100:
                    return sub
        except ValueError:
            for sheet in wb.sheets():
                counter = 0
                for i in range(sheet.nrows):
                    for j in range(sheet.ncols):
                        if sheet.row(i)[j].value == smiles_column:
                            smi_pos = j
                            break
                    if smi_pos is not None:
                        break
                    counter += 1
                    if counter >= 100:
                        return sub
            for sheet in wb.sheets():
                counter = 0
                for i in range(sheet.nrows):
                    for j in range(sheet.ncols):
                        if sheet.row(i)[smi_pos] != "":
                            mol = Chem.MolFromSmiles(sheet.row(i)[j].value)
                            if mol:
                                break
                            mol = Chem.MolFromInchi(sheet.row(i)[j].value)
                            if mol:
                                mol_func = Chem.MolFromInchi
                                break
                    counter += 1
                    if counter >= 100:
                        return sub
    if db:
        titles = []
        if header:
            for sheet in wb.sheets():
                for k in range(sheet.ncols):
                    titles.append(sheet.row(header - 1)[k].value)
        if smi_pos is not None:
            if header:
                for sheet in wb.sheets():
                    for i in range(sheet.nrows):
                        try:
                            mol = mol_func(sheet.row(i)[smi_pos].value)
                            if mol:
                                props = {}
                                for k in range(sheet.ncols):
                                    props[titles[k]] = sheet.row(i)[k].value
                                sub[i] = {"mol": mol, "pattern": sheet.row(i)[smi_pos].value, "props": props}
                        except IndexError:
                            continue
            else:
                for sheet in wb.sheets():
                    for i in range(sheet.nrows):
                        try:
                            mol = mol_func(sheet.row(i)[smi_pos].value)
                            if mol:
                                props = {}
                                for k in range(len(sheet.row(i))):
                                    props[k] = sheet.row(i)[k].value
                                sub[i] = {"mol": mol, "pattern": sheet.row(i)[smi_pos].value, "props": props}
                        except IndexError:
                            continue
    else:
        if smi_pos is not None:
            if mode == "std":
                for sheet in wb.sheets():
                    for i in range(sheet.nrows):
                        try:
                            smi = sheet.row(i)[smi_pos].value
                            if smi != "":
                                mol = mol_func(smi)
                                if mol:
                                    mol_std = query_standardize(mol)
                                    sub[i] = {"mol": mol_std, "pattern": sheet.row(i)[smi_pos].value}
                        except IndexError:
                            continue
            elif mode == "can_taut":
                for sheet in wb.sheets():
                    for i in range(sheet.nrows):
                        try:
                            smi = sheet.row(i)[smi_pos].value
                            if smi != "":
                                mol = mol_func(smi)
                                if mol:
                                    mol_std = query_canonicalize(mol, ntauts)
                                    sub[i] = {"mol": mol_std, "pattern": sheet.row(i)[smi_pos].value}
                        except IndexError:
                            continue
            elif mode == "all_tauts":
                for sheet in wb.sheets():
                    for i in range(sheet.nrows):
                        try:
                            smi = sheet.row(i)[smi_pos].value
                            if smi != "":
                                mol = mol_func(smi)
                                if mol:
                                    mol_std = query_enumerate(mol, ntauts)
                                    sub[i] = {"mol": mol_std[0], "pattern": sheet.row(i)[smi_pos].value,
                                              "tauts": mol_std[1]}
                        except IndexError:
                            continue
            else:
                for sheet in wb.sheets():
                    for i in range(sheet.nrows):
                        try:
                            smi = sheet.row(i)[smi_pos].value
                            if smi != "":
                                mol = mol_func(smi)
                                if mol:
                                    sub[i] = {"mol": mol, "pattern": sheet.row(i)[smi_pos].value}
                        except IndexError:
                            continue
    return sub


def read_file(filename, smiles_column, delimiter, mode, ntauts, smarts=False):
    if filename.endswith(".sdf"):
        sub = read_sd(filename, mode, ntauts)
    elif filename.endswith(".sdf.gz"):
        sub = read_sd(filename, mode, ntauts, gz=True)
    elif filename.endswith(".csv") or filename.endswith(".smi") or filename.endswith(".ich") or filename.endswith(
            ".tsv") or filename.endswith(".txt"):
        sub = read_csv(filename, smiles_column, delimiter, mode, ntauts)
    elif filename.endswith(".csv.gz") or filename.endswith(".tsv.gz") or filename.endswith(".txt.gz"):
        sub = read_csv(filename, smiles_column, delimiter, mode, ntauts, gz=True)
    elif filename.endswith(".xlsx"):
        sub = read_excel(filename, smiles_column, mode, ntauts)
    elif filename.endswith(".sma"):
        if smarts:
            sub = read_smarts_file(filename)
        else:
            sub = {}
    else:
        sub = {}
    return sub


def read_tags(name, tags):
    props = {}
    if name:
        props["Title"] = name
    for line in tags:
        if line.startswith(">  <") and not line.strip("\n").endswith(">"):
            line_strip = line.strip("\n").strip(">  <")
            key = line_strip[:line_strip.index(">")]
            value = tags[tags.index(line) + 1].strip("\n")
            props[key] = value
        elif line.startswith(">  <") and line.strip("\n").endswith(">"):
            key = line.strip("\n").strip("> <")
            value = tags[tags.index(line) + 1].strip("\n")
            props[key] = value
        elif line.startswith("> <") and not line.strip("\n").endswith(">"):
            line_strip = line.strip("\n").strip("> <")
            key = line_strip[:line_strip.index(">")]
            value = tags[tags.index(line) + 1].strip("\n")
            props[key] = value
        elif line.startswith("> <") and line.strip("\n").endswith(">"):
            key = line.strip("\n").strip("> <")
            value = tags[tags.index(line) + 1].strip("\n")
            props[key] = value
    return props


def read_db_from_sd(infile, gz=False):
    sub = {}
    failed = []
    if gz:
        with gzip.open(infile, mode="rt") as inf:
            content = inf.readlines()
    else:
        with open(infile, "r") as sd_file:
            content = sd_file.readlines()
    try:
        sd_blocks = [list(group) for k, group in groupby(content, lambda x: x == "$$$$\n") if not k]
    except ValueError:
        return sub, failed
    del content
    for i in range(len(sd_blocks)):
        mol_block_list = sd_blocks[i][:sd_blocks[i].index("M  END\n") + 1]
        mol_block = ''.join([elem for elem in mol_block_list])
        mol = Chem.MolFromMolBlock(mol_block)
        if mol:
            name = mol.GetProp("_Name")
            tags = sd_blocks[i][sd_blocks[i].index("M  END\n") + 1:]
            props = read_tags(name, tags)
            sub[i] = {"mol": mol, "props": props}
        else:
            failed.append(i)
    return sub, failed


def read_mol_block(block):
    mol_block_list = block[:block.index("M  END\n") + 1]
    mol_block = ''.join([elem for elem in mol_block_list])
    mol = Chem.MolFromMolBlock(mol_block)
    if mol:
        name = mol.GetProp("_Name")
        tags = block[block.index("M  END\n") + 1:]
        props = read_tags(name, tags)
        return {"mol": mol, "props": props}


def read_prepare_mol_block(block):
    mol_block_list = block[:block.index("M  END\n") + 1]
    mol_block = ''.join([elem for elem in mol_block_list])
    mol = Chem.MolFromMolBlock(mol_block)
    if mol:
        name = mol.GetProp("_Name")
        mol2d = Chem.RemoveHs(mol)
        Chem.Compute2DCoords(mol2d)
        tags = block[block.index("M  END\n") + 1:]
        props = read_tags(name, tags)
        if mol.GetConformer().Is3D():
            mol3d = Chem.AddHs(mol, addCoords=True)
            return {"mol": mol2d, "props": props, "confs": mol3d, "pattern": Chem.MolToSmiles(mol2d)}
        else:
            return {"mol": mol2d, "props": props}


def read_3d_mol_block(block):
    mol_block_list = block[:block.index("M  END\n") + 1]
    mol_block = ''.join([elem for elem in mol_block_list])
    mol = Chem.MolFromMolBlock(mol_block)
    if mol:
        if mol.GetConformer().Is3D():
            name = mol.GetProp("_Name")
            tags = block[block.index("M  END\n") + 1:]
            props = read_tags(name, tags)
            mol3d = Chem.AddHs(mol, addCoords=True)
            return {"props": props, "confs": mol3d, "pattern": Chem.MolToSmiles(mol)}


def read_sd_mp(infile, pool, mode="read", gz=False):
    sub = {}
    failed = []
    if gz:
        with gzip.open(infile, mode="rt") as inf:
            content = inf.readlines()
    else:
        with open(infile, "r") as sd_file:
            content = sd_file.readlines()
    try:
        sd_blocks = [list(group) for k, group in groupby(content, lambda x: x == "$$$$\n") if not k]
    except ValueError:
        return sub, failed
    del content
    # pool = mp.Pool(processes=nproc)
    if mode == "read":
        mol_dict = pool.map(read_mol_block, [block for block in sd_blocks])
    elif mode == "3d":
        mol_dict = pool.map(read_3d_mol_block, [block for block in sd_blocks])
    else:
        mol_dict = pool.map(read_prepare_mol_block, [block for block in sd_blocks])
    for i in range(len(mol_dict)):
        if mol_dict[i]:
            sub[i] = mol_dict[i]
        else:
            failed.append(i)
    # pool.close()
    return sub, failed


def read_db_from_sd_3d(infile, gz=False):
    sub = {}
    failed = []
    if gz:
        with gzip.open(infile, mode="rt") as inf:
            content = inf.readlines()
    else:
        with open(infile, "r") as sd_file:
            content = sd_file.readlines()
    try:
        sd_blocks = [list(group) for k, group in groupby(content, lambda x: x == "$$$$\n") if not k]
    except ValueError:
        return sub, failed
    del content
    for i in range(len(sd_blocks)):
        mol_block_list = sd_blocks[i][:sd_blocks[i].index("M  END\n") + 1]
        mol_block = ''.join([elem for elem in mol_block_list])
        mol = Chem.MolFromMolBlock(mol_block)
        if mol:
            if mol.GetConformer().Is3D():
                name = mol.GetProp("_Name")
                mol = Chem.AddHs(mol, addCoords=True)
                tags = sd_blocks[i][sd_blocks[i].index("M  END\n") + 1:]
                props = read_tags(name, tags)
                sub[i] = {"confs": mol, "props": props, "pattern": Chem.MolToSmiles(mol)}
        else:
            failed.append(i)
    return sub, failed


def read_prepare_db_from_sd(infile, gz=False):
    sub = {}
    failed = []
    if gz:
        with gzip.open(infile, mode="rt") as inf:
            content = inf.readlines()
    else:
        with open(infile, "r") as sd_file:
            content = sd_file.readlines()
    try:
        sd_blocks = [list(group) for k, group in groupby(content, lambda x: x == "$$$$\n") if not k]
    except ValueError:
        return sub, failed
    del content
    for i in range(len(sd_blocks)):
        mol_block_list = sd_blocks[i][:sd_blocks[i].index("M  END\n") + 1]
        mol_block = ''.join([elem for elem in mol_block_list])
        mol = Chem.MolFromMolBlock(mol_block)
        if mol:
            name = mol.GetProp("_Name")
            mol2d = Chem.RemoveHs(mol)
            Chem.Compute2DCoords(mol2d)
            tags = sd_blocks[i][sd_blocks[i].index("M  END\n") + 1:]
            props = read_tags(name, tags)
            if mol.GetConformer().Is3D():
                mol3d = Chem.AddHs(mol, addCoords=True)
                sub[i] = {"mol": mol2d, "props": props, "confs": mol3d, "pattern": Chem.MolToSmiles(mol2d)}
            else:
                sub[i] = {"mol": mol2d, "props": props}
        else:
            failed.append(i)
    return sub, failed
