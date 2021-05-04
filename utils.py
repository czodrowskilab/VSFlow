from rdkit.Chem import Descriptors


def check_filter(filter_list, parser):
    choices = ["mw", "logp",
               "nhdo",
               "nhac", "nrob",
               "naro", "nhet",
               "tpsa"]
    filter_dict = {}
    for prop in choices:
        for entry in filter_list:
            if prop in entry:
                if entry.startswith(f"{prop}_"):
                    try:
                        x = float(entry.split("_")[1])
                        filter_dict[prop] = x
                    except ValueError:
                        parser.exit(status=1, message=f"Filter {entry} not supported. Please check the documentation "
                                                      f"which property filters are supported!")
                else:
                    parser.exit(status=1, message=f"Filter {entry} not supported. Please check the documentation which "
                                                  f"property filters are supported!")
    return filter_dict


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


def calc_props(mol, props):
    props["MW (g/mol)"] = str(round(Descriptors.MolWt(mol), 2))
    props["cLogP"] = str(round(Descriptors.MolLogP(mol), 2))
    props["TPSA (A\u00b2)"] = str(round(Descriptors.TPSA(mol), 2))
    props["HDon"] = str(Descriptors.NumHDonors(mol))
    props["HAcc"] = str(Descriptors.NumHAcceptors(mol))
    props["RotBonds"] = str(Descriptors.NumRotatableBonds(mol))
    props["AromRings"] = str(Descriptors.NumAromaticRings(mol))
    props["HetAromRings"] = str(Descriptors.NumAromaticHeterocycles(mol))
