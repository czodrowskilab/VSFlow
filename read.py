from molvs.tautomer import TautomerCanonicalizer, TautomerEnumerator
from molvs.standardize import Standardizer
from rdkit.Chem import AllChem as Chem
from itertools import groupby
from xlrd import open_workbook


def query_standardize(mol):
    try:
        mol_sta = Standardizer().charge_parent(Standardizer().fragment_parent(mol), skip_standardize=True)
        return mol_sta
    except:
        return mol


def query_canonicalize(mol, ntauts):
    try:
        mol_sta = Standardizer().charge_parent(Standardizer().fragment_parent(mol), skip_standardize=True)
        mol_can = TautomerCanonicalizer(max_tautomers=ntauts).canonicalize(mol_sta)
        return mol_can
    except:
        return mol


def query_enumerate(mol, ntauts):
    try:
        mol_sta = Standardizer().charge_parent(Standardizer().fragment_parent(mol), skip_standardize=True)
        mol_tauts = TautomerEnumerator(max_tautomers=ntauts).enumerate(mol_sta)
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
                #sub[i] = (frag_max, smiles[i])
                sub[i] = {"mol": frag_max, "pattern": smiles[i]}
    return sub


def read_smarts(smarts):
    sub = {}
    for i in range(len(smarts)):
        mol = Chem.MolFromSmarts(smarts[i])
        if mol:
            sub[i] = {"mol": mol, "pattern": smarts[i], "tauts": [mol]}
    return sub


def read_smiles_std(smiles):
    sub = {}
    for i in range(len(smiles)):
        mol = Chem.MolFromSmiles(smiles[i])
        if mol:
            mol_std = query_standardize(mol)
            sub[i] = {"mol": mol_std, "pattern": smiles[i]}
    return sub


def read_sd(infile, mode, ntauts):
    sub = {}
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


def read_sd_3d(infile):
    sub = {}
    with open(infile, "r") as sd_file:
        content = sd_file.readlines()
    try:
        sd_blocks = [list(group) for k, group in groupby(content, lambda x: x == "$$$$\n") if not k]
    except ValueError:
        return sub
    del content
    for i in range(len(sd_blocks)):
        mol_block_list = sd_blocks[i][:sd_blocks[i].index("M  END\n") + 1]
        mol_block = ''.join([elem for elem in mol_block_list])
        mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
        if mol:
            sub[i] = {"mol": mol, "pattern": Chem.MolToSmiles(mol)}
    return sub




def read_csv(filename, smiles_column, delimiter, mode, ntauts):
    sub = {}
    with open(filename, "r") as file:
        content = file.readlines()
    if delimiter is None:
        test_del = [";", ",", "\t"]
        for sep in test_del:
            if len(content[0].strip("\n").split(sep)) == len(content[1].strip("\n").split(sep)):
                if len(content[0].strip("\n").split(sep)) > 1:
                    delimiter = sep
        if delimiter is None:
            return sub
            #return "CSV Error: Delimiter could not be guessed"
    proc_content = {}
    for i in range(len(content)):
        proc_content[i] = content[i].strip("\n").split(delimiter)
    del content
    mol_func = Chem.MolFromSmiles
    if smiles_column is None:
        smi_pos = None
        counter = 0
        for n in proc_content:
            for i in range(len(proc_content[n])):
                mol = Chem.MolFromSmiles(proc_content[n][i])
                if mol:
                    #mol_func = Chem.MolFromSmiles
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
                #return "CSV Error: Mol column not found"
    else:
        smi_pos = None
        try:
            smi_pos = int(smiles_column)
            counter = 0
            for n in proc_content:
                mol = Chem.MolFromSmiles(proc_content[n][smi_pos])
                if mol:
                    #mol_func = Chem.MolFromSmiles
                    break
                mol = Chem.MolFromInchi(proc_content[n][smi_pos])
                if mol:
                    mol_func = Chem.MolFromInchi
                    break
                counter += 1
                if counter >= 100:
                    return sub
                    #return "CSV Error: No valid SMILES or InChI found in mol column"
        except ValueError:
            counter = 0
            for n in proc_content:
                if smiles_column in proc_content[n]:
                    smi_pos = proc_content[n].index(smiles_column)
                    break
                counter += 1
                if counter >= 100:
                    return sub
                    #return "CSV Error: No valid SMILES or InChI found in mol column"
            counter = 0
            for n in proc_content:
                mol = Chem.MolFromSmiles(proc_content[n][smi_pos])
                if mol:
                    #mol_func = Chem.MolFromSmiles
                    break
                mol = Chem.MolFromInchi(proc_content[n][smi_pos])
                if mol:
                    mol_func = Chem.MolFromInchi
                    break
                counter += 1
                if counter >= 100:
                    return sub
                    #return "CSV Error: No valid SMILES or InChI found in mol column"
        if smi_pos is None:
            return sub
            #return "CSV Error: Mol Column not found"

    if smi_pos is not None:
        if mode == "std":
            for n in proc_content:
                try:
                    mol = mol_func(proc_content[n][smi_pos])
                    if mol:
                        mol_std = query_standardize(mol)
                        sub[n] = {"mol": mol_std, "pattern": proc_content[n][smi_pos]}
                except IndexError:
                    continue
        elif mode == "can_taut":
            for n in proc_content:
                try:
                    mol = mol_func(proc_content[n][smi_pos])
                    if mol:
                        mol_std = query_canonicalize(mol, ntauts)
                        sub[n] = {"mol": mol_std, "pattern": proc_content[n][smi_pos]}
                except IndexError:
                    continue
        elif mode == "all_tauts":
            for n in proc_content:
                try:
                    mol = mol_func(proc_content[n][smi_pos])
                    if mol:
                        mol_std = query_enumerate(mol, ntauts)
                        sub[n] = {"mol": mol_std[0], "pattern": proc_content[n][smi_pos], "tauts": mol_std[1]}
                except IndexError:
                    continue
        else:
            for n in proc_content:
                try:
                    mol = mol_func(proc_content[n][smi_pos])
                    if mol:
                        sub[n] = {"mol": mol, "pattern": proc_content[n][smi_pos]}
                except IndexError:
                    continue
    return sub


def read_excel(filename, smiles_column, mode, ntauts):
    sub = {}
    wb = open_workbook(filename)
    mol_func = Chem.MolFromSmiles
    if smiles_column is None:
        smi_pos = None
        for sheet in wb.sheets():
            counter = 0
            for i in range(sheet.nrows):
                for j in range(sheet.ncols):
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
                    #return "Excel Error: Mol Column not found"
    else:
        smi_pos = None
        try:
            smi_pos = int(smiles_column)
            for sheet in wb.sheets():
                counter = 0
                for i in range(sheet.nrows):
                    mol = Chem.MolFromSmiles(sheet.row(i)[smiles_column].value)
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
                    #return "Excel Error: No valid SMILES or InChI found in mol column"
        except ValueError:
            for sheet in wb.sheets():
                counter = 0
                for i in range(sheet.nrows):
                    for j in range(sheet.ncols):
                        if sheet.row(i)[j].value == smiles_column:
                            smi_pos = j
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
                        #return "Excel Error: No valid SMILES or InChI found in mol column"
    if smi_pos is not None:
        if mode == "std":
            for sheet in wb.sheets():
                for i in range(sheet.nrows):
                    try:
                        mol = mol_func(sheet.row(i)[smi_pos].value)
                        if mol:
                            mol_std = query_standardize(mol)
                            sub[i] = {"mol": mol_std, "pattern": sheet.row(i)[smi_pos].value}
                    except IndexError:
                        continue
        elif mode == "can_taut":
            for sheet in wb.sheets():
                for i in range(sheet.nrows):
                    try:
                        mol = mol_func(sheet.row(i)[smi_pos].value)
                        if mol:
                            mol_std = query_canonicalize(mol, ntauts)
                            sub[i] = {"mol": mol_std, "pattern": sheet.row(i)[smi_pos].value}
                    except IndexError:
                        continue
        elif mode == "all_tauts":
            for sheet in wb.sheets():
                for i in range(sheet.nrows):
                    try:
                        mol = mol_func(sheet.row(i)[smi_pos].value)
                        if mol:
                            mol_std = query_enumerate(mol, ntauts)
                            sub[i] = {"mol": mol_std[0], "pattern": sheet.row(i)[smi_pos].value, "tauts": mol_std[1]}
                    except IndexError:
                        continue
        else:
            for sheet in wb.sheets():
                for i in range(sheet.nrows):
                    try:
                        mol = mol_func(sheet.row(i)[smi_pos].value)
                        if mol:
                            sub[i] = {"mol": mol, "pattern": sheet.row(i)[smi_pos].value}
                    except IndexError:
                        continue
    return sub


def read_file(filename, smiles_column, delimiter, mode, ntauts):
    if filename.endswith(".sdf"):# or file_format == "sdf":
        sub = read_sd(filename, mode, ntauts)
    elif filename.endswith(".csv"):# or file_format == "csv":
        sub = read_csv(filename, smiles_column, delimiter, mode, ntauts)
    elif filename.endswith(".xlsx"):# or file_format == "xlsx":
        sub = read_excel(filename, smiles_column, mode, ntauts)
    else:
        sub = {}
    return sub


def read_db_from_sd(infile):
    sub = {}
    failed = []
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
        if mol.GetConformer().Is3D():
            mol3d = Chem.AddHs(mol, addCoords=True)
            return {"mol": mol2d, "props": props, "confs": mol3d, "pattern": Chem.MolToSmiles(mol2d)}
        else:
            return {"mol": mol2d, "props": props}


def read_sd_mp(infile, pool, mode="read"):
    sub = {}
    failed = []
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
    else:
        mol_dict = pool.map(read_prepare_mol_block, [block for block in sd_blocks])
    for i in range(len(mol_dict)):
        if mol_dict[i]:
            sub[i] = mol_dict[i]
        else:
            failed.append(i)
    # pool.close()
    return sub, failed


def read_db_from_sd_3d(infile):
    sub = {}
    failed = []
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
            mol = Chem.AddHs(mol, addCoords=True)
            tags = sd_blocks[i][sd_blocks[i].index("M  END\n") + 1:]
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
            sub[i] = {"confs": mol, "props": props, "pattern": Chem.MolToSmiles(mol)}
        else:
            failed.append(i)
    return sub, failed


def read_prepare_db_from_sd(infile):
    sub = {}
    failed = []
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
            if mol.GetConformer().Is3D():
                mol3d = Chem.AddHs(mol, addCoords=True)
                sub[i] = {"mol": mol2d, "props": props, "confs": mol3d, "pattern": Chem.MolToSmiles(mol2d)}
            else:
                sub[i] = {"mol": mol2d, "props": props}
        else:
            failed.append(i)
    return sub, failed
