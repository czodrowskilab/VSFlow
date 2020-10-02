import argparse
import csv
import math
import multiprocessing as mp
import os
import time
from subprocess import run, PIPE, Popen

import pandas as pd
import requests
import xlsxwriter
from bs4 import BeautifulSoup
from fpdf import FPDF
from pdfrw import PdfReader
from pdfrw import PdfWriter
from pymol import cmd
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw, MACCSkeys
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Draw import SimilarityMaps

RDLogger.logger().setLevel(RDLogger.CRITICAL)


## Read path and identity tags of integrated databases in DATABASES.csv

try:
    DATABASES = {}
    IDENTITY = {}
    with open("DATABASES.csv", "r") as file:
        content = file.readlines()
    for line in csv.reader(content):
        if "name" and "path" and "identity" in line:
            continue
        else:
            DATABASES[line[0]] = line[1]
            try:
                IDENTITY[line[0]] = line[2]
            except IndexError:
                IDENTITY[line[0]] = ""
    default_db = list(DATABASES.keys())[0]
except FileNotFoundError:
    DATABASES = {}
    IDENTITY = {}
    default_db = ""

parser = argparse.ArgumentParser(description="VSFlow")
subparsers = parser.add_subparsers(title="mode", help="specify mode of vsflow")

# Visualize a database as pdf

visualize = subparsers.add_parser("visualize", description="visualize structures in a sdf file as pdf")
visualize.add_argument("-in", "--input", help="specify input file", required=True)
visualize.add_argument("-out", "--output", help="specify name of output pdf file", default="visualization.pdf")
visualize.add_argument("-np", "--mpi_np", default=4, type=int,
                       help="Specifies the number of processors n when the application is run in"
                            " MPI mode. [Default = 4]")


def return_mols(test_data):
    results = []
    for mol in test_data:
        if mol:
            props = mol.GetPropsAsDict()
            results.append([mol, props])
    return results


def gen_grid_mols_2(i, pool_results):
    Chem.Compute2DCoords(pool_results[i][0])
    grid_mol = Draw.MolsToGridImage([pool_results[i][0]], molsPerRow=1, subImgSize=(600, 600))
    filename = f"grid_mol_{i}.jpeg"
    grid_mol.save(filename)
    return filename


def gen_grid_mols(i, pool_results):
    Chem.Compute2DCoords(pool_results[i][0])
    grid_mol = Draw.MolsToGridImage([pool_results[i][0]], highlightAtomLists=[pool_results[i][1]], molsPerRow=1,
                                    subImgSize=(600, 600))
    filename = f"grid_mol_{i}.jpeg"
    grid_mol.save(filename)
    return filename


def export_page_id(counter, image_list, results):
    page_mols = image_list[counter:counter + 6]
    page_props = results[counter:counter + 6]
    img_place_y_even = 10
    txt_place_y_even = 16
    num_place_y_even = 98
    img_place_y_odd = 10
    txt_place_y_odd = 16
    num_place_y_odd = 98
    pdf = FPDF()
    pdf.add_page()
    pdf.add_font("DejaVu", "", "DejaVuSansMono.ttf", uni=True)
    pdf.set_font("DejaVu", "", 11)
    for i in range(len(page_mols)):
        div = str(i / 2).rsplit(".")[1]
        if div == str(0):
            pdf.image(f"grid_mol_{counter + i}.jpeg", 10, img_place_y_even, 90)
            pdf.rect(10, img_place_y_even, 90, 90)
            pdf.text(12, num_place_y_even, str(counter + i + 1))
            pdf.text(12, txt_place_y_even, page_props[i][2])
            img_place_y_even += 94
            txt_place_y_even += 94
            num_place_y_even += 94
        else:
            pdf.image(f"grid_mol_{counter + i}.jpeg", 110, img_place_y_odd, 90)
            pdf.rect(110, img_place_y_odd, 90, 90)
            pdf.text(112, num_place_y_odd, str(counter + i + 1))
            pdf.text(112, txt_place_y_odd, page_props[i][2])
            img_place_y_odd += 94
            txt_place_y_odd += 94
            num_place_y_odd += 94
    out = f"page_{counter}.pdf"  # counter + mol_count
    pdf.output(out)
    return out


def export_page(counter, image_list, results, prop_dict):
    page_mols = image_list[counter:counter + 3]
    page_props = results[counter:counter + 3]
    img_place_y = 10
    txt_place_y = 6
    num_place_y = 98
    pdf = FPDF()
    pdf.add_page()
    pdf.add_font("DejaVu", "", "DejaVuSansMono.ttf", uni=True)
    pdf.set_font("DejaVu", "", 11)
    for i in range(len(page_mols)):
        pdf.image(page_mols[i], 10, img_place_y, 90)
        pdf.rect(10, img_place_y, 190, 90)
        pdf.dashed_line(100, img_place_y, 100, img_place_y + 90)
        pdf.text(12, num_place_y, str(counter + i + 1))
        img_place_y += 94
        num_place_y += 94
        txt_space_y = 10
        for tag in page_props[i][prop_dict].items():
            string = f"{tag[0]}: {tag[1]}"
            str_fact = 0
            str_space = 0
            if len(string) < 40:
                pdf.text(103, txt_place_y + txt_space_y, string)
                txt_space_y += 6
            else:
                if len(string) <= 80:
                    while str_fact < len(string) - 40:
                        pdf.text(103, txt_place_y + txt_space_y + str_space, string[str_fact:str_fact + 40])
                        str_space += 5
                        str_fact += 40
                    if str_fact >= len(string) - 40:
                        pdf.text(103, txt_place_y + txt_space_y + str_space, string[str_fact:len(string)])
                        str_space += 5
                    txt_space_y = txt_space_y + str_space + 1
                else:
                    write_string = string[:77] + "..."
                    while str_fact < len(write_string) - 40:
                        pdf.text(103, txt_place_y + txt_space_y + str_space, write_string[str_fact:str_fact + 40])
                        str_space += 5
                        str_fact += 40
                    if str_fact >= len(write_string) - 40:
                        pdf.text(103, txt_place_y + txt_space_y + str_space, write_string[str_fact:len(write_string)])
                        str_space += 5
                    txt_space_y = txt_space_y + str_space + 1
        txt_place_y += 94
    out = f"page_{counter}.pdf"
    pdf.output(out)
    return out


def write_pdf(pages, out_file):
    pdf_writer = PdfWriter()
    for page in pages:
        pdf_reader = PdfReader(page)
        pdf_writer.addPage(pdf_reader.getPage(0))
    with open(out_file, "wb") as file:
        pdf_writer.write(file)


def gen_counter(grids, split):
    x = math.ceil(len(grids)/split)
    counter = []
    for i in range(x):
        counter.append(i * 3)
    return counter


def pdf_out(pages, last_idx, out_file, mol_factor):
    if len(pages) <= 400:
        write_pdf(pages, out_file)
    else:
        mol_count = 0
        out_prefix = out_file.rsplit(".")[0]
        while mol_count < len(pages) - 400:
            file_pages = pages[
                         mol_count:mol_count + 400]  # oder Seitenzahl als output, dann muss last_idx nicht übergeben werden
            outfile = f"{out_prefix}_{mol_count * mol_factor + 1}_{mol_count * mol_factor + 1200}.pdf"
            write_pdf(file_pages, outfile)
            mol_count += 400
        if mol_count >= len(pages) - 400:
            file_pages = pages[mol_count:len(pages)]
            outfile = f"{out_prefix}_{mol_count * mol_factor + 1}_{last_idx}.pdf"
            write_pdf(file_pages, outfile)
    for file in os.listdir("."):
        if file.startswith("grid_mol_") or file.startswith("page_"):
            os.remove(file)


def gen_pdf(args):
    start_time = time.time()
    print(f"Start: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    data = Chem.SDMolSupplier(args.input, sanitize=False, strictParsing=False)
    mols = return_mols(data)
    print(f"Finished reading database: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    pool = mp.Pool(processes=args.mpi_np)
    grids = pool.starmap(gen_grid_mols_2, [(i, mols) for i in range(len(mols))])
    print(f"Finished generating grid mols: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    counter = gen_counter(grids, 3)
    pages = pool.starmap(export_page, [(j, grids, mols, 1) for j in counter])
    print(f"Finished generating pdf pages: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    pdf_out(pages, len(grids), args.output, 3)
    print(f"Finished generating pdf file(s): {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    pool.close()
    end_time = time.time()
    duration = round(end_time - start_time)
    print(f"Finished in {duration} seconds")


visualize.set_defaults(func=gen_pdf)

## Substructure search

substructure = subparsers.add_parser("substructure", description="perform a substructure search")
group = substructure.add_mutually_exclusive_group(required=True)
group.add_argument("-in", "--input", help="input file")
group.add_argument("-smi", "--smiles", help="specify smiles for substructure search", action="append")
group.add_argument("-sma", "--smarts", help="specify smarts for substructure search", action="append")
substructure.add_argument("-out", "--output",
                          help="specify name of output file. Supported formats are sdf, csv and xlsx."
                               "If no file extension is provided, all possible files are generated."
                          , default="vsflow_substructure.sdf")
substructure.add_argument("-np", "--mpi_np", default=4, type=int,
                          help="Specifies the number of processors n when the application is run in"
                               " MPI mode. [Default = 4]")
substructure.add_argument("-pdf", "--PDF", help="generate a pdf file for substructure matches", action="store_true")
substructure.add_argument("-db", "--database", help="select database or provide path to database", default=default_db)
substructure.add_argument("-id", "--identity",
                          help="specify the identity tag for molecules in database, only required if non-implemented database is used")
substructure.add_argument("-props", "--properties",
                          help="specifies if calculated molecular properties are written to the output files",
                          action="store_true")
substructure.add_argument("-mf", "--multfile", help="generate separate output files for every query molecule",
                          action="store_true")
substructure.add_argument("-fi", "--input_format",
                          help="Specify file typ if no file extension is present in input file name")
substructure.add_argument("-col", "--smiles_column",
                          help="Specify name of smiles column in csv file",
                          default="smiles")
substructure.add_argument("-del", "--delimiter", help="Specify delimiter of csv file", default=";")
substructure.add_argument("-head", "--header", help="Specify row of file to be used as column names", type=int,
                          default=1)
substructure.add_argument("-fm", "--fullmatch", help="when specified, only full matches are returned",
                          action="store_true")
substructure.add_argument("-filt", "--filter", help="specify property to filter screening results", action="append")


def read_smiles(smiles):
    sub = []
    for i in range(len(smiles)):
        mol = Chem.MolFromSmiles(smiles[i])
        frags = Chem.GetMolFrags(mol, asMols=True)
        frag_max = max(frags, key=lambda m: m.GetNumAtoms())
        sub.append((frag_max, i))
    return sub


def read_smarts(smarts):
    sub = []
    for i in range(len(smarts)):
        mol = Chem.MolFromSmarts(smarts[i])
        sub.append((mol, i))
    return sub


def read_file(filename, file_format, smiles_column, delimiter, header):
    sub = []
    if filename.endswith(".sdf") or file_format == "sdf":
        mols = Chem.SDMolSupplier(filename)
        for i in range(len(mols)):
            frags = Chem.GetMolFrags(mols[i], asMols=True)
            frag_max = max(frags, key=lambda m: m.GetNumAtoms())
            sub.append((frag_max, i))
    elif filename.endswith(".csv") or file_format == "csv":
        csv_df = pd.read_csv(filename, header=header - 1, delimiter=delimiter)
        for i, rn in csv_df.iterrows():
            mol = Chem.MolFromSmiles(rn[smiles_column])
            frags = Chem.GetMolFrags(mol, asMols=True)
            frag_max = max(frags, key=lambda m: m.GetNumAtoms())
            sub.append((frag_max, i))
        # eigener Code
        # content = []
        # with open(filename, "r") as file:
        #     for line in file:
        #         int_cont = line.strip().split(delimiter)
        #         content.append(int_cont)
        # col_names = content.pop(header - 1)
        # if smiles_column in col_names:
        #     index = col_names.index(smiles_column)
        # else:
        #     substructure.error(message="Please specify (correct) name of smiles column (-col)")
        # sub = []
        # for i in range(len(content)):
        #     mol = Chem.MolFromSmiles(content[i][index])
        #     sub.append((mol, i))
    elif filename.endswith(".xls") or filename.endswith(".xlsx") or file_format == "xls":
        xls = pd.read_excel(filename, header=header - 1)
        for i, rn in xls.iterrows():
            mol = Chem.MolFromSmiles(rn[smiles_column])
            frags = Chem.GetMolFrags(mol, asMols=True)
            frag_max = max(frags, key=lambda m: m.GetNumAtoms())
            sub.append((frag_max, i))

    else:
        if file_format is None:
            substructure.error(message="Please specify input file format (-fi) [sdf, csv and xlsx are supported]")
    return sub


def substruct(mol, name, query, query_num, fullmatch):
    match = tuple(j for k in mol.GetSubstructMatches(query) for j in k)
    if match:
        if fullmatch:
            frags = Chem.GetMolFrags(mol, asMols=True)
            frag_max = max(frags, key=lambda m: m.GetNumAtoms())
            if frag_max.GetNumAtoms() == len(match):
                query_smi = Chem.MolToSmiles(query)
                return [mol, match, name, query_smi, query_num]
        else:
            query_smi = Chem.MolToSmiles(query)
            return [mol, match, name, query_num, query_smi]




def read_input(smiles, smarts, infile, input_format, smiles_column, delimiter, header):
    if smiles:
        query = read_smiles(smiles)
    elif smarts:
        query = read_smarts(smarts)
    else:
        query = read_file(infile, input_format, smiles_column, delimiter, header)
    return query


def set_props_sub(properties, sub_results, identity):
    results = []
    if properties:
        for mol in sub_results:
            if mol:
                mol[0].SetProp(identity, mol[2])
                mol[0].SetProp("MW (g/mol)", str(round(Descriptors.MolWt(mol[0]), 2)))
                mol[0].SetProp("cLogP", str(round(Descriptors.MolLogP(mol[0]), 2)))
                mol[0].SetProp("TPSA (A\u00b2)", str(round(Descriptors.TPSA(mol[0]), 2)))
                mol[0].SetProp("HDon", str(Descriptors.NumHDonors(mol[0])))
                mol[0].SetProp("HAcc", str(Descriptors.NumHAcceptors(mol[0])))
                mol[0].SetProp("RotBonds", str(Descriptors.NumRotatableBonds(mol[0])))
                mol[0].SetProp("AromRings", str(Descriptors.NumAromaticRings(mol[0])))
                mol[0].SetProp("HetAromRings", str(Descriptors.NumAromaticHeterocycles(mol[0])))
                mol[0].SetProp("_Name", mol[2])
                mol[0].SetProp("QuerySmiles", mol[4])
                mol[0].SetProp("MatchAtoms", f"{mol[1]}")
                mol.append(mol[0].GetPropsAsDict())
                results.append(mol)
    else:
        for mol in sub_results:
            if mol:
                mol[0].SetProp(identity, mol[2])
                mol[0].SetProp("_Name", mol[2])
                mol[0].SetProp("QuerySmiles", mol[4])
                mol[0].SetProp("MatchAtoms", f"{mol[1]}")
                mol.append(mol[0].GetPropsAsDict())
                results.append(mol)
    return results


def sdf_write(pool_results, output):
    writer = Chem.SDWriter(output)
    for res in pool_results:
        writer.write(res[0])
    writer.close()


def sdf_write_self(pool_results, output):
    mol_blocks = []
    for entry in pool_results:
        block = Chem.MolToMolBlock(entry[0]).rstrip("\n").split("\n")
        for line in block:
            if line.startswith("     RDKit          2D"):
                new_line = "   VSFlow   version1.0"
                block[block.index("     RDKit          2D")] = new_line
        for tag in entry[5].items():
            block.append(f">  <{tag[0]}>")
            block.append(f"{tag[1]}")
            block.append("")
        mol_blocks.append(block)
    with open(output, "w") as writefile:
        for block in mol_blocks:
            for line in block:
                writefile.write(f"{line}\n")
            writefile.write("$$$$\n")


def write_csv(pool_results, output):
    header = []
    values = []
    for element in pool_results[0][5]:
        header.append(element)
    header.append("Smiles")
    for entry in pool_results:
        entry[5]["Smiles"] = Chem.MolToSmiles(entry[0])
        values.append(entry[5].values())
    with open(output, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=";")
        csvwriter.writerow(header)
        for val in values:
            csvwriter.writerow(val)


def write_xls(pool_results, output):
    header = []
    values = []
    for element in pool_results[0][5]:
        header.append(element)
    header.append("Smiles")
    for entry in pool_results:
        entry[5]["Smiles"] = Chem.MolToSmiles(entry[0])
        int_vals = []
        for val in entry[5].values():
            int_vals.append(val)
        values.append(int_vals)
    workbook = xlsxwriter.Workbook(output)
    worksheet = workbook.add_worksheet()
    for i in range(len(header)):
        worksheet.write(0, i, header[i])
    for j in range(len(values)):
        for i in range(len(values[j])):
            worksheet.write(j + 1, i, values[j][i])
    workbook.close()


def sub_pdf(pool_sub, grids, results, out_name, properties):
    if properties:
        counter = gen_counter(grids, 3)
        pages = pool_sub.starmap(export_page, [(j, grids, results, 5) for j in counter])
        mol_factor = 3
    else:
        counter = gen_counter(grids, 6)
        pages = pool_sub.starmap(export_page_id, [(j, grids, results) for j in counter])
        mol_factor = 6
    pdf_out(pages, len(grids), f"{out_name}.pdf", mol_factor)


def read_db(database):
    if database in DATABASES:
        data = Chem.SDMolSupplier(DATABASES[database])
    else:
        try:
            data = Chem.SDMolSupplier(database)
        except FileNotFoundError:
            substructure.error(message="File not found. Please make sure you specified the correct path!")
    return data


def read_id(database, identity):
    if database in IDENTITY:
        if IDENTITY[database]:
            ident = IDENTITY[database]
        else:
            ident = "_Name"
    else:
        if identity:
            ident = identity
        else:
            ident = "_Name"
    return ident


def prop_filt(filtered, filter_dict, results, prop_func, key):
    int_filt = []
    if filtered:
        for entry in filtered[-1]:
            if prop_func(entry[0]) <= filter_dict[key]:
                int_filt.append(entry)
    else:
        for entry in results:
            if prop_func(entry[0]) <= filter_dict[key]:
                int_filt.append(entry)
    return int_filt


def check_filter(filter_list):
    choices = ["mw", "logp",
               "hdon",
               "hacc", "rotb",
               "narom", "nhet",
               "tpsa"]
    filter_dict = {}
    for prop in choices:
        for entry in filter_list:
            if prop in entry:
                if entry.startswith(prop):
                    try:
                        x = float(entry.split(prop)[1])
                        filter_dict[prop] = x
                    except ValueError:
                        parser.error(message=f"Filter {entry} not supported")
    return filter_dict


def filter_res(filter_dict, results):
    filtered = []
    if "mw" in filter_dict:
        int_filt = prop_filt(filtered, filter_dict, results, Descriptors.MolWt, "mw")
        filtered.append(int_filt)
    if "logp" in filter_dict:
        int_filt = prop_filt(filtered, filter_dict, results, Descriptors.MolLogP, "logp")
        filtered.append(int_filt)
    if "hdon" in filter_dict:
        int_filt = prop_filt(filtered, filter_dict, results, Descriptors.NumHDonors, "hdon")
        filtered.append(int_filt)
    if "hacc" in filter_dict:
        int_filt = prop_filt(filtered, filter_dict, results, Descriptors.NumHAcceptors, "hacc")
        filtered.append(int_filt)
    if "rotb" in filter_dict:
        int_filt = prop_filt(filtered, filter_dict, results, Descriptors.NumRotatableBonds, "rotb")
        filtered.append(int_filt)
    if "narom" in filter_dict:
        int_filt = prop_filt(filtered, filter_dict, results, Descriptors.NumAromaticRings, "narom")
        filtered.append(int_filt)
    if "nhet" in filter_dict:
        int_filt = prop_filt(filtered, filter_dict, results, Descriptors.NumAromaticHeterocycles, "nhet")
        filtered.append(int_filt)
    if "tpsa" in filter_dict:
        int_filt = prop_filt(filtered, filter_dict, results, Descriptors.TPSA, "tpsa")
        filtered.append(int_filt)
    return filtered[-1]


def output_files(output, results):
    if output.endswith(".sdf"):
        sdf_write_self(results, output)
    elif output.endswith(".csv"):
        write_csv(results, output)
    elif output.endswith(".xlsx"):
        write_xls(results, output)
    else:
        sdf_write_self(results, f"{output}.sdf")
        write_csv(results, f"{output}.csv")
        write_xls(results, f"{output}.xlsx")


def sub_search(args):
    start_time = time.time()
    print(f"Start: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    # check input
    if args.filter:
        filter_dict = check_filter(args.filter)
    # read input files or input strings
    data = read_db(args.database)
    identity = read_id(args.database, args.identity)
    query = read_input(args.smiles, args.smarts, args.input, args.input_format, args.smiles_column, args.delimiter,
                       args.header)
    print(f"Finished reading query molecules: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    data_names = []
    for mol in data:
        if mol:
            name = mol.GetProp(identity)
            data_names.append((mol, name))
    argslist = [(mol[0], mol[1], query_mol[0], query_mol[1], args.fullmatch) for mol in data_names for query_mol in
                query]
    print(f"Finished reading database: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    # perform substructure search
    pool_sub = mp.Pool(processes=args.mpi_np)
    sub_results = pool_sub.starmap(substruct, argslist)
    print(f"Finished substructure search: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    # sort and filter substructure results
    results = sorted([entry for entry in sub_results if entry], key=lambda entry: entry[3])
    if args.filter:
        results = filter_res(filter_dict, results)
    # set properties for output files
    results = set_props_sub(args.properties, results, identity)
    # generate output files
    if args.fullmatch:
        print(f"{len(results)} full matches found")
    else:
        print(f"{len(results)} substructure matches found")
    if results:
        print(results)
        if args.output.endswith(".sdf") or args.output.endswith(".csv") or args.output.endswith(".xlsx"):
            out_prefix = args.output.rsplit(".", maxsplit=1)[0]
            out_ext = args.output.rsplit(".", maxsplit=1)[1]
        else:
            out_prefix = args.output
            out_ext = ""
        if args.multfile:
            for i in range(len(query)):
                output_mols = []
                for entry in results:
                    if entry[3] == i:
                        output_mols.append(entry)
                if output_mols:
                    out_name = f"{out_prefix}_query_{i + 1}.{out_ext}"
                    out_name_pdf = f"{out_prefix}_query_{i + 1}"
                    output_files(out_name, output_mols)
                    if args.PDF:
                        grids = pool_sub.starmap(gen_grid_mols, [(i, output_mols) for i in range(len(output_mols))])
                        sub_pdf(pool_sub, grids, output_mols, out_name_pdf, args.properties)
        else:
            output_files(args.output, results)
            if args.PDF:
                grids = pool_sub.starmap(gen_grid_mols, [(i, results) for i in range(len(results))])
                sub_pdf(pool_sub, grids, results, out_prefix, args.properties)
        print(f"Finished generating output file(s): {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    else:
        print(f"No substructure matches found for {len(query)} query molecules")
    pool_sub.close()
    end_time = time.time()
    duration = round(end_time - start_time)
    print(f"Finished in {duration} seconds")


substructure.set_defaults(func=sub_search)

## Fingerprint similarity search

fp_sim = subparsers.add_parser("fpsim", description="molecular similarity search using fingerprints")
group_fp = fp_sim.add_mutually_exclusive_group(required=True)
group_fp.add_argument("-in", "--input", help="input file")
group_fp.add_argument("-smi", "--smiles", help="specify smiles on command line", action="append")
group_fp.add_argument("-sma", "--smarts", help="specify smiles on command line", action="append")
fp_sim.add_argument("-np", "--mpi_np", default=4, type=int,
                    help="Specifies the number of processors n when the application is run in"
                         " MPI mode. [Default = 4]")
fp_sim.add_argument("-out", "--output", help="specify name of output file", default="vsflow_fingerprint.sdf")
fp_sim.add_argument("-pdf", "--PDF", help="generate a pdf file for substructure matches", action="store_true")
fp_sim.add_argument("-db", "--database", help="select database", default=default_db)
fp_sim.add_argument("-id", "--identity",
                    help="specify the identity tag for molecules in database, only required if non-implemented database is used")
fp_sim.add_argument("-props", "--properties",
                    help="specifies if calculated molecular properties are written to the output files",
                    action="store_true")
fp_sim.add_argument("-mf", "--multfile", help="generate separate output files for every query molecule",
                    action="store_true")
fp_sim.add_argument("-fi", "--input_format", help="Specify file typ if no file extension is present in input file name")
fp_sim.add_argument("-col", "--smiles_column", help="Specify name of smiles column in csv file", default="smiles")
fp_sim.add_argument("-del", "--delimiter", help="Specify delimiter of csv file", default=";")
fp_sim.add_argument("-head", "--header", help="Specify row of file to be used as column names", type=int, default=1)
fp_sim.add_argument("-fp", "--fingerprint", help="specify fingerprint to be used", choices=["rdkit", "ecfp", "fcfp", "ap", "tt", "maccs"],
                    default="fcfp")
fp_sim.add_argument("-sim", "--similarity", help="specify fingerprint similarity metric to be used",
                    choices=["tan", "dice", "cos", "sok", "russ", "kulc", "mcco", "tver"], default="tan")
fp_sim.add_argument("-fpr", "--fpradius", help="radius of circular fingerprints ecfp and fcfp", type=int,
                    default=3)
fp_sim.add_argument("-nbits", "--NBITS", help="number of bits used to generate ecfp and fcfp fingerprints", type=int,
                    default=4096)
fp_sim.add_argument("-top", "--top_hits", type=int, default=10,
                    help="Maximum number of molecules with highest similarity to keep. [Default = 10]")
fp_sim.add_argument("-map", "--simmap", help="generates similarity maps for fingerprints in pdf file", action="store_true")
fp_sim.add_argument("-cut", "--cutoff", help="specify cutoff value for similarity coefficient", type=float)
fp_sim.add_argument("-filt", "--filter", help="specify property to filter screening results", action="append")


def query_fp_rdkit(query_mol, query_num, nbits):
    q_fp = Chem.RDKFingerprint(query_mol, fpSize=nbits)
    return (query_mol, query_num, q_fp)


def database_fp_rdkit(mol, name, nbits):
    mol_fp = Chem.RDKFingerprint(mol, fpSize=nbits)
    return (mol, mol_fp, name)


def query_fp_ecfp(query_mol, query_num, nbits, radius):
    q_fp = Chem.GetMorganFingerprintAsBitVect(query_mol, radius, nBits=nbits)
    return (query_mol, query_num, q_fp)


def database_fp_ecfp(mol, name, nbits, radius):
    mol_fp = Chem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    return (mol, mol_fp, name)


def query_fp_fcfp(query_mol, query_num, nbits, radius):
    q_fp = Chem.GetMorganFingerprintAsBitVect(query_mol, radius, nBits=nbits, useFeatures=True)
    return (query_mol, query_num, q_fp)


def database_fp_fcfp(mol, name, nbits, radius):
    mol_fp = Chem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits, useFeatures=True)
    return (mol, mol_fp, name)


def query_fp_ap(query_mol, query_num):
    q_fp = Pairs.GetAtomPairFingerprintAsBitVect(query_mol)
    return (query_mol, query_num, q_fp)


def database_fp_ap(mol, name):
    mol_fp = Pairs.GetAtomPairFingerprintAsBitVect(mol)
    return (mol, mol_fp, name)


def query_fp_tt(query_mol, query_num, path_length):
    q_fp = Torsions.GetTopologicalTorsionFingerprint(query_mol, targetSize=path_length)
    return (query_mol, query_num, q_fp)


def database_fp_tt(mol, name, path_length):
    mol_fp = Torsions.GetTopologicalTorsionFingerprint(mol, targetSize=path_length)
    return (mol, mol_fp, name)

def query_fp_maccs(query_mol, query_num):
    q_fp = MACCSkeys.GenMACCSKeys(query_mol)
    return (query_mol, query_num, q_fp)


def database_fp_maccs(mol, name):
    mol_fp = MACCSkeys.GenMACCSKeys(mol)
    return (mol, mol_fp, name)


def sim_tan(mol, mol_fp, name, query_fp, query_num):
    coef = DataStructs.TanimotoSimilarity(query_fp, mol_fp)
    return [mol, coef, name, query_num]


def sim_dice(mol, mol_fp, name, query_fp, query_num):
    coef = DataStructs.DiceSimilarity(query_fp, mol_fp)
    return [mol, coef, name, query_num]


def sim_tver(mol, mol_fp, name, query_fp, query_num):
    coef = DataStructs.TverskySimilarity(query_fp, mol_fp)
    return [mol, coef, name, query_num]


def sim_cos(mol, mol_fp, name, query_fp, query_num):
    coef = DataStructs.CosineSimilarity(query_fp, mol_fp)
    return [mol, coef, name, query_num]


def sim_sok(mol, mol_fp, name, query_fp, query_num):
    coef = DataStructs.SokalSimilarity(query_fp, mol_fp)
    return [mol, coef, name, query_num]


def sim_russ(mol, mol_fp, name, query_fp, query_num):
    coef = DataStructs.RusselSimilarity(query_fp, mol_fp)
    return [mol, coef, name, query_num]


def sim_kulc(mol, mol_fp, name, query_fp, query_num):
    coef = DataStructs.KulczynskiSimilarity(query_fp, mol_fp)
    return [mol, coef, name, query_num]


def sim_mcco(mol, mol_fp, name, query_fp, query_num):
    coef = DataStructs.McConnaugheySimilarity(query_fp, mol_fp)
    return [mol, coef, name, query_num]


def fp_sim_metric(mol, mol_fp, name, query_fp, query_num, metric):
    coef = metric(query_fp, mol_fp)
    return [mol, coef, name, query_num]


def set_props_fp(properties, sub_results, identity, similarity, fingerprint):
    results = []
    if properties:
        for mol in sub_results:
            if mol:
                mol[0].SetProp(identity, mol[2])
                mol[0].SetProp("MW (g/mol)", str(round(Descriptors.MolWt(mol[0]), 2)))
                mol[0].SetProp("cLogP", str(round(Descriptors.MolLogP(mol[0]), 2)))
                mol[0].SetProp("TPSA (A\u00b2)", str(round(Descriptors.TPSA(mol[0]), 2)))
                mol[0].SetProp("HDon", str(Descriptors.NumHDonors(mol[0])))
                mol[0].SetProp("HAcc", str(Descriptors.NumHAcceptors(mol[0])))
                mol[0].SetProp("RotBonds", str(Descriptors.NumRotatableBonds(mol[0])))
                mol[0].SetProp("AromRings", str(Descriptors.NumAromaticRings(mol[0])))
                mol[0].SetProp("HetAromRings", str(Descriptors.NumAromaticHeterocycles(mol[0])))
                mol[0].SetProp("_Name", mol[2])
                mol[0].SetProp("QuerySmiles", mol[4])
                mol[0].SetProp("Fingerprint", fingerprint)
                mol[0].SetProp("Similarity", f"{round(mol[1], 5)} ({similarity})")
                mol.append(mol[0].GetPropsAsDict())
                results.append(mol)
    else:
        for mol in sub_results:
            if mol:
                mol[0].SetProp(identity, mol[2])
                mol[0].SetProp("_Name", mol[2])
                mol[0].SetProp("QuerySmiles", mol[4])
                mol[0].SetProp("Fingerprint", fingerprint)
                mol[0].SetProp("Similarity", f"{round(mol[1], 5)} ({similarity})")
                mol.append(mol[0].GetPropsAsDict())
                results.append(mol)
    return results


def sim_map(results, fp_func, metric):
    grids = []
    for i in range(len(results)):
        fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(Chem.MolFromSmiles(results[i][4]), results[i][0],
                                                                       fp_func, metric=metric)
        fig.set_figwidth(3.255)
        fig.set_figheight(3.255)
        filename = f"grid_mol_{i}.jpeg"
        fig.savefig(filename, bbox_inches="tight")
        grids.append(filename)
    return grids


def fp_gen_pdf(fp_pool, results, metric, simmap, fingerprint, fpradius, nbits, features, output, properties):
    if simmap:  # generate similarity maps for supported fingerprints
        if fingerprint == "ecfp" or fingerprint == "fcfp":
            fp_func = lambda m, idx: SimilarityMaps.GetMorganFingerprint(m, atomId=idx, radius=fpradius,
                                                                         fpType='bv', nBits=nbits,
                                                                         useFeatures=features)
            grids = sim_map(results, fp_func, metric)
        elif fingerprint == "rdkit":
            fp_func = lambda m, idx: SimilarityMaps.GetRDKFingerprint(m, atomId=idx, fpType="bv", nBits=nbits)
            grids = sim_map(results, fp_func, metric)
        elif fingerprint == "ap":
            fp_func = lambda m, idx: SimilarityMaps.GetAPFingerprint(m, atomId=idx, fpType="bv", nBits=nbits)
            grids = sim_map(results, fp_func, metric)
        elif fingerprint == "tt":
            fp_func = lambda m, idx: SimilarityMaps.GetTTFingerprint(m, atomId=idx, fpType="bv", nBits=nbits)
            grids = sim_map(results, fp_func, metric)
        else:
            grids = fp_pool.starmap(gen_grid_mols_2, [(i, results) for i in range(len(results))])
    else:
        grids = fp_pool.starmap(gen_grid_mols_2, [(i, results) for i in range(len(results))])
    sub_pdf(fp_pool, grids, results, output, properties)


def fp_search(args):
    start_time = time.time()
    print(f"Start: {time.strftime('%d/%m/%Y, %H:%M:%S', time.localtime())}")
    # check input arguments
    if args.cutoff:
        if args.cutoff > 1:
            fp_sim.error(message="Similarity cutoff value must be between 0 and 1")
    if args.filter:
        filter_dict = check_filter(args.filter)
    # read input files or input strings
    data = read_db(args.database)
    identity = read_id(args.database, args.identity)
    query = read_input(args.smiles, args.smarts, args.input, args.input_format, args.smiles_column, args.delimiter,
                       args.header)
    print(f"Finished reading query molecules: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    data_names = []
    for mol in data:
        if mol:
            name = mol.GetProp(identity)
            data_names.append((mol, name))
    print(f"Finished reading database: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    # calculate fingerprints and similarities
    fp_pool = mp.Pool(processes=args.mpi_np)
    if args.fingerprint == "rdkit":
        fp_name = f"RDKitFingerprint {args.NBITS} bit"
        fp_query = fp_pool.starmap(query_fp_rdkit, [(query_mol[0], query_mol[1], args.NBITS) for query_mol in query])
        data_fp = fp_pool.starmap(database_fp_rdkit, [(mol[0], mol[1], args.NBITS) for mol in data_names])
    elif args.fingerprint == "ecfp":
        fp_name = f"ECFP{args.fpradius*2} {args.NBITS} bit"
        features = False
        fp_query = fp_pool.starmap(query_fp_ecfp, [(query_mol[0], query_mol[1], args.NBITS, args.fpradius) for query_mol in query])
        data_fp = fp_pool.starmap(database_fp_ecfp, [(mol[0], mol[1], args.NBITS, args.fpradius) for mol in data_names])
    elif args.fingerprint == "fcfp":
        fp_name = f"FCFP{args.fpradius*2} {args.NBITS} bit"
        features = True
        fp_query = fp_pool.starmap(query_fp_fcfp, [(query_mol[0], query_mol[1], args.NBITS, args.fpradius) for query_mol in query])
        data_fp = fp_pool.starmap(database_fp_fcfp, [(mol[0], mol[1], args.NBITS, args.fpradius) for mol in data_names])
    elif args.fingerprint == "ap":
        fp_name = "AtomPairs"
        fp_query = fp_pool.starmap(query_fp_ap, [(query_mol[0], query_mol[1]) for query_mol in query])
        data_fp = fp_pool.starmap(database_fp_ap, [(mol[0], mol[1]) for mol in data_names])
    elif args.fingerprint == "tt":
        fp_name = "TopTorsion"
        fp_query = fp_pool.starmap(query_fp_tt, [(query_mol[0], query_mol[1], args.pathlength) for query_mol in query])
        data_fp = fp_pool.starmap(database_fp_tt, [(mol[0], mol[1], args.pathlength) for mol in data_names])
    elif args.fingerprint == "maccs":
        fp_name = "MACCS"
        fp_query = fp_pool.starmap(query_fp_maccs, [(query_mol[0], query_mol[1]) for query_mol in query])
        data_fp = fp_pool.starmap(database_fp_maccs, [(mol[0], mol[1]) for mol in data_names])
    argslist = [(mol[0], mol[1], mol[2], fp[2], fp[1]) for mol in data_fp for fp in fp_query]
    if args.similarity == "tan":
        similarity = "Tanimoto"
        metric = DataStructs.TanimotoSimilarity
        fp_results = fp_pool.starmap(sim_tan, argslist)
    elif args.similarity == "dice":
        similarity = "Dice"
        fp_results = fp_pool.starmap(sim_dice, argslist)
    elif args.similarity == "cos":
        similarity = "Cosine"
        fp_results = fp_pool.starmap(sim_cos, argslist)
    elif args.similarity == "sok":
        similarity = "Sokal"
        fp_results = fp_pool.starmap(sim_sok, argslist)
    elif args.similarity == "russ":
        similarity = "Russel"
        fp_results = fp_pool.starmap(sim_russ, argslist)
    elif args.similarity == "kulc":
        similarity = "Kulczynski"
        fp_results = fp_pool.starmap(sim_kulc, argslist)
    elif args.similarity == "mcco":
        similarity = "McConnaughey"
        fp_results = fp_pool.starmap(sim_mcco, argslist)
    elif args.similarity == "tver":
        similarity = "Tversky"
        fp_results = fp_pool.starmap(sim_tver, argslist)
    print(f"Finished fingerprint search: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    # sort and filter results
    sorted_results = []
    for i in range(len(query)):
        query_results = []
        for entry in fp_results:
            if entry[3] == i:
                query_results.append(entry)
        if args.cutoff:
            sort_results = sorted(query_results, key=lambda entry: entry[1], reverse=True)
            for entry in sort_results:
                if entry[1] >= args.cutoff:
                    q_smi = Chem.MolToSmiles(query[i][0])
                    entry.append(q_smi)
                    sorted_results.append(entry)
        else:
            sort_results = sorted(query_results, key=lambda entry: entry[1], reverse=True)[:args.top_hits]
            for element in sort_results:
                q_smi = Chem.MolToSmiles(query[i][0])
                element.append(q_smi)
                sorted_results.append(element)
    if args.filter:
        sorted_results = filter_res(filter_dict, sorted_results)
    # set properties for output files
    results = set_props_fp(args.properties, sorted_results, identity, similarity, fp_name)
    # generate output files
    if results:
        if args.output.endswith(".sdf") or args.output.endswith(".csv") or args.output.endswith(".xlsx"):
            out_prefix = args.output.rsplit(".", maxsplit=1)[0]
            out_ext = args.output.rsplit(".", maxsplit=1)[1]
        else:
            out_prefix = args.output
            out_ext = ""
        if args.multfile:
            for i in range(len(query)):
                output_mols = []
                for entry in results:
                    if entry[3] == i:
                        output_mols.append(entry)
                if output_mols:
                    out_name = f"{out_prefix}_query_{i + 1}.{out_ext}"
                    out_name_pdf = f"{out_prefix}_query_{i + 1}"
                    print(out_name)
                    output_files(out_name, output_mols)
                    if args.PDF:
                        fp_gen_pdf(fp_pool, output_mols, metric, args.simmap, args.fingerprint, args.fpradius, args.NBITS,
                                   features, out_name_pdf, args.properties)
                else:
                    print(f"No matches found for query molecule {i}")
        else:
            output_files(args.output, results)
            if args.PDF:
                fp_gen_pdf(fp_pool, results, metric, args.simmap, args.fingerprint, args.fpradius, args.NBITS, features,
                           out_prefix, args.properties)
        print(f"Finished generating output files: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    else:
        print("No similarity matches found for query molecules with the specified filter criteria")
    fp_pool.close()
    end_time = time.time()
    duration = round(end_time - start_time)
    print(f"Finished in {duration} seconds")


fp_sim.set_defaults(func=fp_search)


## generate 3D coordinates with RDKit

rd3d = subparsers.add_parser("rd3d", description="generate 3D coordinates for input molecules")
group_rd3d = rd3d.add_mutually_exclusive_group(required=True)
group_rd3d.add_argument("-in", "--input", help="input file")
group_rd3d.add_argument("-smi", "--smiles", help="specify smiles on command line", action="append")
group_rd3d.add_argument("-sma", "--smarts", help="specify smiles on command line", action="append")
rd3d.add_argument("-out", "--output", help="specify name of output file", default="structures_3D.sdf")
rd3d.add_argument("-col", "--smiles_column", help="Specify name of smiles column in csv file", default="smiles")
rd3d.add_argument("-del", "--delimiter", help="Specify delimiter of csv file", default=";")
rd3d.add_argument("-head", "--header", help="Specify row of file to be used as column names", type=int, default=1)
rd3d.add_argument("-fi", "--input_format", help="Specify file typ if no file extension is present in input file name")
rd3d.add_argument("-np", "--mpi_np", default=4, type=int,
                    help="Specifies the number of processors n when the application is run in"
                         " MPI mode. [Default = 4]")
rd3d.add_argument("-confs", "--conformers", help="specify number of conformers to be generated for each molecule",
                  type=int)
rd3d.add_argument("-ff", "--ffield", help="specify forcefield to be used for energy minimization", choices=["MMFF94",
                                                                                                            "MMFF94s",
                                                                                                            "UFF"],
                  default="MMFF94")


def gen_coords_mmff(mol, num, ffield):
    mol_H = Chem.AddHs(mol)
    Chem.EmbedMolecule(mol_H)
    Chem.MMFFOptimizeMolecule(mol_H, mmffVariant=ffield)
    mol_3D = Chem.RemoveHs(mol_H)
    return (mol_3D, num)


def gen_coords_uff(mol, num):
    mol_H = Chem.AddHs(mol)
    Chem.EmbedMolecule(mol_H)
    Chem.UFFOptimizeMolecule(mol_H)
    mol_3D = Chem.RemoveHs(mol_H)
    return (mol_3D, num)


def gen_confs_mmff(mol, num, confs, ffield):
    mol_H = Chem.AddHs(mol)
    Chem.EmbedMultipleConfs(mol_H, numConfs=confs)
    Chem.MMFFOptimizeMoleculeConfs(mol_H, mmffVariant=ffield)
    mol_3D = Chem.RemoveHs(mol_H)
    return (mol_3D, num)


def gen_confs_uff(mol, num, confs):
    mol_H = Chem.AddHs(mol)
    Chem.EmbedMultipleConfs(mol_H, numConfs=confs)
    Chem.UFFOptimizeMoleculeConfs(mol_H)
    mol_3D = Chem.RemoveHs(mol_H)
    return (mol_3D, num)


def sdf_write_3d(pool_results, output):
    mol_blocks = []
    for entry in pool_results:
        block = Chem.MolToMolBlock(entry[0]).rstrip("\n").split("\n")
        for line in block:
            if line.startswith("     RDKit          2D"):
                new_line = "   VSFlow   version1.0"
                block[block.index("     RDKit          2D")] = new_line
        for tag in entry[0].GetPropsAsDict().items():
            block.append(f">  <{tag[0]}>")
            block.append(f"{tag[1]}")
            block.append("")
        mol_blocks.append(block)
    with open(output, "w") as writefile:
        for block in mol_blocks:
            for line in block:
                writefile.write(f"{line}\n")
            writefile.write("$$$$\n")


def sdf_writer_conf(mol_list, out_file):
    writer = Chem.SDWriter(out_file)
    for entry in mol_list:
        confs = entry[0].GetConformers()
        for i in range(len(confs)):
            writer.write(entry[0], confId=i)
    writer.close()


def gen_3d(args):
    mols = read_input(args.smiles, args.smarts, args.input, args.input_format, args.smiles_column, args.delimiter, args.header)
    print(Chem.MolToMolBlock(mols[0][0]))
    pool_3d = mp.Pool(processes=args.mpi_np)
    if args.conformers:
        if args.ffield == "MMFF94" or args.ffield == "MMFF94s":
            mols_3d = pool_3d.starmap(gen_confs_mmff, [(mol[0], mol[1], args.conformers, args.ffield) for mol in mols])
        else:
            mols_3d = pool_3d.starmap(gen_confs_uff, [(mol[0], mol[1], args.conformers) for mol in mols])
        results_3d = sorted(mols_3d, key=lambda entry: entry[1])
        sdf_writer_conf(results_3d, args.output)
    else:
        if args.ffield == "MMFF94" or args.ffield == "MMFF94s":
            mols_3d = pool_3d.starmap(gen_coords_mmff, [(mol[0], mol[1], args.ffield) for mol in mols])
        else:
            mols_3d = pool_3d.starmap(gen_coords_uff, [(mol[0], mol[1]) for mol in mols])
        results_3d = sorted(mols_3d, key=lambda entry: entry[1])
        sdf_write_3d(results_3d, args.output)
    pool_3d.close()


rd3d.set_defaults(func=gen_3d)

## PDB database search  ## Beta version: under development

pdb = subparsers.add_parser("pdblig", description="download ligands from a pdb entry")
group_pdb = pdb.add_mutually_exclusive_group(required=True)
group_pdb.add_argument("-in", "--input", help="specify input text file containing comma-separated pdb codes")
group_pdb.add_argument("-id", "--pdbid", help="specifiy pdb code on command line", action="append")
group_pdb.add_argument("-s", "--search", help="specify text to search for in pdb database")
pdb.add_argument("-out", "--output", help="specifiy name of output file")
pdb.add_argument("-ex", "--exclude", help="if set, common small molecule impurties (e.g. GOL, ACT etc.) are not considered as ligands,"
                                          " ions are never considered as ligands.", action="store_true")


def read_pdb(infile, pdbid):
    if infile:
        with open(infile) as file:
            content = file.read()
            pdb_ids = content.split(",")
    if pdbid:
        pdb_ids = pdbid
    return pdb_ids


def req_lig(pdb_code, exclude):
    r = requests.get(f"https://www.rcsb.org/structure/{pdb_code}")
    if exclude:
        non_ligs = ['PG4', '1PE', 'PEG', 'PGE', 'KHZ', 'DMS', 'ACT', 'ZN"', 'SO4', 'GOL', 'CA"', 'EDO', 'OXL', 'MG"', 'K">',
                'NA"', 'CL"', 'BR"', 'I">']
    else:
        non_ligs = []
    entry_ligs = []
    if r.status_code == 200:
        counter = r.text.find("ligand_row_")
        ligs = []
        while counter != -1:
            lig_name = r.text[counter + 11:counter + 14]
            ligs.append(lig_name)
            counter = r.text.find("ligand_row_", counter + 1)
        for lig in ligs:
            if not lig in non_ligs:
                lig_r = requests.get(f"https://www.rcsb.org/ligand/{lig}")
                if lig_r.status_code == 200:
                    doc = BeautifulSoup(lig_r.text, "html5lib")
                    smi = doc.select_one("#chemicalIsomeric td").text
                    ident = doc.select_one("#chemicalIdentifiers td").text
                    name = doc.select_one("#chemicalName td").text
                    mw = doc.select_one("#chemicalMolecularWeight td").text
                    inchi = doc.select_one("#chemicalInChI td").text
                    mol = Chem.MolFromSmiles(smi)
                    mol.SetProp("IUPAC", ident)
                    mol.SetProp("_Name", name)
                    mol.SetProp("ID", lig)
                    mol.SetProp("PDB", pdb_code)
                    mol.SetProp("Smiles", smi)
                    mol.SetProp("MW", mw)
                    mol.SetProp("InChI", inchi)
                    entry_ligs.append(mol)
    else:
        print(f"Invalid pdb code. No entry found for pdb code {pdb_code}.")
    return entry_ligs


def search_pdb(args):
    start_time = time.time()
    print(f"Start: {time.strftime('%d/%m/%Y, %H:%M:%S', time.localtime())}")
    r_rest = requests.get(
        "http://www.rcsb.org/pdb/rest/customReport.xml?pdbids=*&customReportColumns=structureId,ligandId,ligandSmiles,"
        "InChI,ligandName,uniprotAcc,uniprotRecommendedName,uniprotAlternativeNames&service=wsfile&format=csv")
    print(f"Finished reading database: {time.strftime('%d/%m/%Y, %H:%M:%S', time.localtime())}")
    content = r_rest.text.split("\n")
    head = content.pop(0)
    pop = content.pop(-1)
    result_mols = []
    non_ligs = ['PO4', 'PG4', '1PE', 'PEG', 'PGE', 'KHZ', 'DMS', 'ACT', 'ZN', 'SO4', 'GOL', 'CA', 'EDO', 'OXL', 'MG',
                'K', 'NA', 'CL', 'BR', 'I', 'PG6']
    ids = set()
    for line in csv.reader(content):
        for entry in line:
            if args.search.lower() in entry.lower():
                ids.add(line[0])
    results = sorted(ids)
    print(f"Finished searching database: {time.strftime('%d/%m/%Y, %H:%M:%S', time.localtime())}")
    files = []
    for entry in results:
        cmd.fetch(entry, type="pdb")
        print(entry)
    # for entry in results:
    #     cmd.save(f"{entry}.pdb", selection=f"{entry}")
    print(f"Finished saving pdb files: {time.strftime('%d/%m/%Y, %H:%M:%S', time.localtime())}")
                # if line[2] in non_ligs:
                #     continue
                # else:
                #     print(line)
                    # mol = Chem.MolFromInchi(line[4])
                    # mol.SetProp("_Name", line[2])
                    # mol.SetProp("ID", line[2])
                    # mol.SetProp("Smiles", line[3])
                    # mol.SetProp("InChI", line[4])
                    # mol.SetProp("PDB_ID", line[0])
                    # mol.SetProp("ChemicalName", line[5])
                    # mol.SetProp("UniprotID", line[6])
                    # result_mols.append((mol, line[2], line[3], line[4], line[5], line[0], line[6]))



def get_ligs(args):
    pdb_list = read_pdb(args.input, args.pdbid)
    results_ligs = []
    for code in pdb_list:
        code_ligs = req_lig(code, args.exclude)
        results_ligs.append(code_ligs)
    print(results_ligs)


pdb.set_defaults(func=search_pdb)


# Omega

omega = subparsers.add_parser("omega", description="generate 3D conformers for input molecules")
group_2 = omega.add_mutually_exclusive_group(required=True)
group_2.add_argument("-in", "--input", help="input file")
group_2.add_argument("-smi", "--smiles", help="specify smiles for substructure search")
omega.add_argument("-flag", "--FLAG", help="classic, macrocycle, rocs, pose or dense", choices=["classic","macrocycle","rocs","pose","dense"], default="classic")
omega.add_argument("-out", "--output", help="name of output file", default="omega_results")
omega.add_argument("-maxconfs", default="50", help="Maximum number of conformations to be saved. [Default = 50]")
omega.add_argument("-mpi_np", default="6", help="Specifies the number of processors n when the application is run in"
                                                 " MPI mode. [Default = 6]")
omega.add_argument("-fo", "--output_format", help="specify output file format", choices=["sdf", "oeb", "oeb.gz"], default="sdf")


def run_omega(args):
    if args.output == "omega_results":
        output = f"{args.output}_{args.FLAG}.{args.output_format}"
    elif args.output.endswith(".sdf") or args.output.endswith(".oeb") or args.output.endswith(".oeb.gz"):
        output = args.output
    else:
        output = f"{args.output}.{args.output_format}"
    if args.smiles:
        smi = args.smiles.split(",")
        l = []
        for smiles in smi:
            x = smiles+"\n"
            l.append(x)
        # noinspection PyArgumentList
        p = Popen(["oeomega", args.FLAG, "-out", output, "-maxconfs", args.maxconfs, "-in", ".txt", "-mpi_np", args.mpi_np], stdin=PIPE, text=True)
        f = p.stdin
        f.writelines(l)
        f.close()
        p.communicate()
        p.wait()
    else:
        run(["oeomega", args.FLAG, "-in", args.input, "-out", output, "-maxconfs", args.maxconfs, "-mpi_np", args.mpi_np])


omega.set_defaults(func=run_omega)


## ROCS

rocs = subparsers.add_parser("rocs", description="perform ROCS screening")
group_3 = rocs.add_mutually_exclusive_group(required=True)
group_3.add_argument("-in", "--input", help="input file")
group_3.add_argument("-smi", "--smiles", help="specify smiles for substructure search")
rocs.add_argument("-db", "--database", help="select database", choices=["comas_3d"], default="comas_3d")
rocs.add_argument("-out", "--output", help="name of output file without file extension", default="rocs")
rocs.add_argument("-top", "--top_hits", help="number of top hits to keep", default="100")
rocs.add_argument("-fo", "--output_format", help="specify output file format", choices=["sdf", "oeb", "oeb.gz"], default="sdf")
rocs.add_argument("-py", "--pymol", help="generate pymol session file for rocs results", choices=["ob_1", "ob_split"])
rocs.add_argument("-report", "--REPORT", help="visualizes rocs hits as pdf file", action="store_true")
rocs.add_argument("-mpi_np", default="6", help="Specifies the number of processors n when the application is run in"
                                                 " MPI mode. [Default = 6]")


def rocs_hits(wd_list, output):
    hit_list = []
    for filename in wd_list:
        if f"{output}_hits" in filename:
            hit_list.append(filename)
    return hit_list


def run_rocs(args):
    run(["rocs", "-query", args.input, "-dbase", DATABASES[args.database], "-prefix", args.output, "-besthits", args.top_hits, "-oformat", args.output_format, "-mpi_np", args.mpi_np])
    if args.pymol:
        wd_list = os.listdir(".")
        hits = rocs_hits(wd_list, args.output)
        for file in hits:
            if file.endswith(".sdf") or file.endswith(".mol2") or file.endswith(".pdb"):
                py_object = file.split(".")[0]
                cmd.load(filename=file)
                if args.pymol == "ob_1":
                    cmd.save(f"{py_object}.pse")
                else:
                    cmd.split_states(object=py_object)
                    cmd.save(f"{py_object}_split_states.pse")
    if args.REPORT:
        wd_list = os.listdir(".")
        hits = rocs_hits(wd_list, args.output)
        for file in hits:
            out_file = file.split(".")[0]
            run(["rocs_report", "-in", file, "-out", f"{out_file}.pdf"])


rocs.set_defaults(func=run_rocs)


## Omega ROCS combination


omrocs = subparsers.add_parser("omrocs", description="run rocs startinf from 2d structure query")
group_4 = omrocs.add_mutually_exclusive_group(required=True)
group_4.add_argument("-in", "--input", help="input file")
group_4.add_argument("-smi", "--smiles", help="specify smiles on command line")
# omrocs.add_argument("-flag", "--FLAG", help="classic, macrocycle, rocs, pose or dense", choices=["classic","macrocycle","rocs","pose","dense"], default="classic")
omrocs.add_argument("-maxconfs", default="50", help="Maximum number of conformations to be saved. [Default = 50]")
omrocs.add_argument("-db", "--database", help="select database", choices=["comas_3d"], default="comas_3d")
omrocs.add_argument("-out", "--output", help="name of output file without file extension", default="omega_rocs")
omrocs.add_argument("-top", "--top_hits", help="number of top hits to keep", default="100")
omrocs.add_argument("-fo", "--output_format", help="specify output file format", choices=["sdf", "oeb", "oeb.gz"], default="sdf")
omrocs.add_argument("-mpi_np", default="6", help="Specifies the number of processors n when the application is run in"
                                                 " MPI mode. [Default = 6]")


def run_omega_rocs(args):
    if args.smiles:
        smi = args.smiles.split(",")
        l = []
        for smiles in smi:
            x = smiles+"\n"
            l.append(x)
        # noinspection PyArgumentList
        p = Popen(["oeomega", "rocs", "-out", "query.sdf", "-maxconfs", args.maxconfs, "-in", ".txt", "-mpi_np", args.mpi_np], stdin=PIPE, text=True)
        f = p.stdin
        f.writelines(l)
        f.close()
        p.communicate()  # stdout, _ = p.communicate()
        p.wait()
    else:
        run(["oeomega", "rocs", "-in", args.input, "-out", "query.sdf", "-maxconfs", args.maxconfs, "-mpi_np", args.mpi_np])

    run(["rocs", "-query", "query.sdf", "-dbase", DATABASES[args.database], "-prefix", args.output, "-besthits",
         args.top_hits, "-oformat", args.output_format, "-mpi_np", args.mpi_np])


omrocs.set_defaults(func=run_omega_rocs)


def main():
    args = parser.parse_args()
    if "func" in args:
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
