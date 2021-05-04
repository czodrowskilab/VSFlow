import csv
import os

import xlsxwriter
from rdkit.Chem import AllChem as Chem


def write_sdf(mol, props, output):
    block = Chem.MolToMolBlock(mol).rstrip("\n").split("\n")
    try:
        block[0] = props["Title"]
    except KeyError:
        pass
    block[1] = "   VSFlow 1.0 (RDKit 2D)"
    for tag in props.items():
        block.append(f">  <{tag[0]}>")
        block.append(f"{tag[1]}")
        block.append("")
    for line in block:
        output.write(f"{line}\n")
    output.write("$$$$\n")


def write_sdf_conformer(mol, props, confId, output):
    block = Chem.MolToMolBlock(mol, confId=confId).rstrip("\n").split("\n")
    try:
        block[0] = props["Title"]
    except KeyError:
        pass
    block[1] = "   VSFlow 1.0 (RDKit 3D)"
    for tag in props.items():
        block.append(f">  <{tag[0]}>")
        block.append(f"{tag[1]}")
        block.append("")
    for line in block:
        output.write(f"{line}\n")
    output.write("$$$$\n")


def write_excel(lines, sorteddict, out_file):
    workbook = xlsxwriter.Workbook(out_file)
    worksheet = workbook.add_worksheet()
    for i in range(len(list(sorteddict.keys()))):
        worksheet.write(0, i, list(sorteddict.keys())[i])
    for j in range(len(lines)):
        for k in range(len(lines[j])):
            try:
                worksheet.write(j + 1, k, float(lines[j][k]))
            except ValueError:
                worksheet.write(j + 1, k, lines[j][k])
    workbook.close()


def write_csv(lines, sorteddict, out_file):
    with open(out_file, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=";")
        csvwriter.writerow(list(sorteddict.keys()))
        for line in lines:
            csvwriter.writerow(line)


def sort_props(props, sorteddict):
    for prop in props:
        if prop not in sorteddict:
            try:
                sorteddict[prop] = max(sorteddict.values()) + 1
            except ValueError:
                sorteddict[prop] = 0


def prepare_lines(mol, props, sorteddict, lines):
    line = [""] * len(sorteddict)
    line[0] = Chem.MolToSmiles(mol)
    for prop in props:
        line[sorteddict[prop]] = props[prop]
    lines.append(line)


def gen_csv_xls_mult(query, results, output):
    out_prefix = output.rsplit(".", maxsplit=1)[0]
    out_ext = output.rsplit(".", maxsplit=1)[1]
    if out_ext == "csv":
        writer = write_csv
    else:
        writer = write_excel
    counter = 1
    for m in query:
        out_file = f"{out_prefix}_{counter}.{out_ext}"
        lines = []
        sorteddict = {"Smiles": 0}
        for n in results:
            if results[n]["q_num"] == m:
                sort_props(results[n]["props"], sorteddict)
                prepare_lines(results[n]["mol"], results[n]["props"], sorteddict, lines)
        if lines:
            writer(lines, sorteddict, out_file)
            counter += 1


def gen_sdf_mult(query, results, output):
    out_prefix = output.rsplit(".sdf", maxsplit=1)[0]
    counter = 1
    for m in query:
        out_file = f"{out_prefix}_{counter}.sdf"
        query_count = 0
        with open(out_file, "w") as output:
            for n in results:
                if results[n]["q_num"] == m:
                    query_count += 1
                    write_sdf(results[n]["mol"], results[n]["props"], output)
        if query_count > 0:
            counter += 1
        else:
            os.remove(f"{out_prefix}_{counter}.sdf")


def gen_csv_xls(results, output):
    if output.endswith(".csv"):
        writer = write_csv
    else:
        writer = write_excel
    lines = []
    sorteddict = {"Smiles": 0}
    for n in results:
        sort_props(results[n]["props"], sorteddict)
        prepare_lines(results[n]["mol"], results[n]["props"], sorteddict, lines)
    if lines:
        writer(lines, sorteddict, output)


def gen_sdf(results, output):
    out_file = output.rsplit(".sdf", maxsplit=1)[0]
    with open(f"{out_file}.sdf", "w") as out:
        for n in results:
            write_sdf(results[n]["mol"], results[n]["props"], out)
