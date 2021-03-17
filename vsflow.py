import argparse
import csv
import math
import multiprocessing as mp
import os
import time
from subprocess import run, PIPE, Popen
from itertools import groupby
import pickle
import random

import pandas as pd
import requests
import xlsxwriter
#from bs4 import BeautifulSoup
from fpdf import FPDF, set_global
#from pdfrw import PdfReader
#from pdfrw import PdfWriter
from pymol import cmd
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw, MACCSkeys
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Draw import SimilarityMaps
from molvs.tautomer import TautomerCanonicalizer, TautomerEnumerator
from molvs.standardize import Standardizer
import matplotlib as mpl
from rdkit.Chem import rdMolAlign
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from xlrd import open_workbook
mpl.rc('figure', max_open_warning=0)
RDLogger.logger().setLevel(RDLogger.CRITICAL)
import visualize
import sss
import read
import write_output
import fpsearch
import  shapesearch
from rdkit.Chem.Draw import rdMolDraw2D


# set path and global variables
script_path = os.path.dirname(os.path.abspath(__file__))
home = os.path.expanduser("~")
#database_path = f"{home}/.vsflow/DATABASES.csv"
#database_global = f"{script_path}/DATABASES.csv"
ttf_path = f"{script_path}/resources/DejaVuSansMono.ttf"
set_global("FPDF_CACHE_MODE", 2)
set_global("FPDF_CACHE_DIR", script_path)

## Read path of integrated databases or generate generate/update config files
if os.path.exists(f"{home}/.vsflow"):
    config = pickle.load(open(f"{home}/.vsflow/.config", "rb"))
    db_config = pickle.load(open(f"{home}/.vsflow/.db_config", "rb"))
    db_default = pickle.load(open(f"{home}/.vsflow/.db_default", "rb"))
else:
    os.mkdir(f"{home}/.vsflow")
    config = {"global_db": f"{script_path}/DATABASES", "local_db": f"{home}/VSFlow_Databases"}
    db_config = {}
    db_default = ""
print(db_config)
for key in config:
    try:
        for file in os.listdir(config[key]):
            if file.endswith(".vsdb"):
                if file.rsplit(".vsdb", maxsplit=1)[0] in db_config:
                    if db_config[file.rsplit(".vsdb", maxsplit=1)[0]][0] == time.ctime(os.path.getmtime(f"{config[key]}/{file}")):
                        continue
                print("Updating databases")
                db = pickle.load(open(f"{config[key]}/{file}", "rb"))
                db_info = db["config"]
                n_mols = db_info[2]
                standardized = db_info[0]
                n_confs = db_info[1]
                num_seed = db_info[3]
                db_config[file.rsplit(".vsdb", maxsplit=1)[0]] = [time.ctime(os.path.getmtime(f"{config[key]}/{file}")),
                                                                  standardized,
                                                                  n_confs,
                                                                  n_mols,
                                                                  num_seed]
                del db
    except FileNotFoundError:
        continue
to_remove = []
for db_info in db_config:
    try:
        if f"{db_info}.vsdb" not in os.listdir(config["global_db"]):
            if f"{db_info}.vsdb" not in os.listdir(config["local_db"]):
                #db_config.pop(db_info)
                to_remove.append(db_info)
    except FileNotFoundError:
        try:
            if f"{db_info}.vsdb" not in os.listdir(config["global_db"]):
                #db_config.pop(db_info)
                to_remove.append(db_info)
        except FileNotFoundError:
            try:
                if f"{db_info}.vsdb" not in os.listdir(config["local_db"]):
                    #db_config.pop(db_info)
                    to_remove.append(db_info)
            except FileNotFoundError:
                pass
for entry in to_remove:
    db_config.pop(entry)
    if entry == db_default:
        db_default = ""
pickle.dump(config, open(f"{home}/.vsflow/.config", "wb"))
pickle.dump(db_config, open(f"{home}/.vsflow/.db_config", "wb"))
pickle.dump(db_default, open(f"{home}/.vsflow/.db_default", "wb"))
print(config)
print(db_config)
print(db_default)
parser = argparse.ArgumentParser(description="Virtual Screening Workflow")
print('''\
**************************

 VV        VV  SSSSSSS             VSFlow   
  VV      VV  SSS    SS       Virtual Screening
   VV    VV    SSSS               Workflow
    VV  VV       SSSS         
     VVVV     SS    SSS       
      VV       SSSSSSS           

**************************
''')
subparsers = parser.add_subparsers(title="mode", help="specify mode of vsflow")

substructure = subparsers.add_parser("substructure", description="perform a substructure search")
group = substructure.add_mutually_exclusive_group(required=True)
group.add_argument("-in", "--input", help="input file")
group.add_argument("-smi", "--smiles", help="specify smiles for substructure search", action="append")
group.add_argument("-sma", "--smarts", help="specify smarts for substructure search", action="append")
substructure.add_argument("-out", "--output",
                          help="specify name of output file. Supported formats are sdf, csv and xlsx."
                               "If no file extension is provided, all possible files are generated."
                          , default="vsflow_substructure.sdf")
substructure.add_argument("-np", "--mpi_np", type=int,
                          help="Specify the number of processors n to run the application in"
                               " MPI mode.")
substructure.add_argument("-pdf", "--PDF", help="generate a pdf file for substructure matches", action="store_true")
substructure.add_argument("-db", "--database", help="select database or provide path to database", default=db_default)
substructure.add_argument("-props", "--properties",
                          help="specifies if calculated molecular properties are written to the output files",
                          action="store_true")
substructure.add_argument("-mf", "--multfile", help="generate separate output files for every query molecule",
                          action="store_true")
substructure.add_argument("-fi", "--input_format",
                          help="Specify file typ if no file extension is present in input file name")
substructure.add_argument("-col", "--smiles_column",
                          help="Specify name of smiles column in csv file")
substructure.add_argument("-del", "--delimiter", help="Specify delimiter of csv file")
# substructure.add_argument("-head", "--header", help="Specify row of file to be used as column names", type=int,
#                           default=1)
substructure.add_argument("-fm", "--fullmatch", help="when specified, only full matches are returned",
                          action="store_true")
substructure.add_argument("-filt", "--filter", help="specify property to filter screening results", action="append")
substructure.add_argument("-nt", "--ntauts", help="maximum number of tautomers of query molecules to be enumerated", type=int, default=100)
substructure.add_argument("-m", "--mode", help="choose a mode for substructure search", choices=["std", "all_tauts",
                                                                                                 "can_taut", "no_std"],
                          default="no_std")


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


def calc_props(mol, props):
    props["MW (g/mol)"] = str(round(Descriptors.MolWt(mol), 2))
    props["cLogP"] = str(round(Descriptors.MolLogP(mol), 2))
    props["TPSA (A\u00b2)"] = str(round(Descriptors.TPSA(mol), 2))
    props["HDon"] = str(Descriptors.NumHDonors(mol))
    props["HAcc"] = str(Descriptors.NumHAcceptors(mol))
    props["RotBonds"] = str(Descriptors.NumRotatableBonds(mol))
    props["AromRings"] = str(Descriptors.NumAromaticRings(mol))
    props["HetAromRings"] = str(Descriptors.NumAromaticHeterocycles(mol))


def read_database(args, pool):
    mols = {}
    if args.database in db_config:
        try:
            mols = pickle.load(open(f"{config['local_db']}/{args.database}.vsdb", "rb"))
        except FileNotFoundError:
            try:
                mols = pickle.load(open(f"{config['global_db']}/{args.database}.vsdb", "rb"))
            except FileNotFoundError:
                substructure.error(
                    message=f"{args.database} not found. Please make sure you specified the correct shortcut")
    else:
        if os.path.exists(args.database):
            if args.database.endswith(".vsdb"):
                try:
                    mols = pickle.load(open(args.database, "rb"))
                except:
                    substructure.error(message=f"{args.output} could not be opened. Please make sure the file has the correct "
                                               f"format")

            else:
                if args.mpi_np:
                    mols, failed = read.read_sd_mp(args.database, pool)
                else:
                    mols, failed = read.read_db_from_sd(args.database)
                if failed:
                    print(f"{len(failed)} of {len(mols) + len(failed)} molecules in {args.database} could not be processed")
                if not mols:
                    substructure.error(message="No molecules could be read from SD file. Please make sure it has the right "
                                               "format")
        else:
            substructure.error(message=f"File {args.database} not found. Please make sure you specified the correct path")
    return mols


def read_input(args):
    if args.smarts:
        query = read.read_smarts(args.smarts)
        if not query:
            substructure.error(message="No valid molecule(s) could be generated from the provided SMARTS.")
    elif args.smiles:
        query = read.read_smiles(args.smiles, args.mode, args.ntauts)
        if not query:
            substructure.error(message="No valid molecule(s) could be generated from the provided SMILES.")
    else:
        if os.path.exists(args.input):
            query = read.read_file(args.input, args.input_format, args.smiles_column, args.delimiter, args.mode, args.ntauts)
            if not query:
                if args.input.endswith(".sdf") or args.input_format == "sdf":
                    substructure.error(message="No valid molecules could be read from SD file.")
                elif args.input.endswith(".csv") or args.input_format == "csv":
                    substructure.error(message="No valid molecules could be read from input file. Please check/specify "
                                               "name of SMILES/InChI containing column (--mol_column) or check/specify the"
                                               "separator (--delimiter)")
                elif args.input.endswith(".xlsx") or args.input_format == "xlsx":
                    substructure.error(message="No valid molecules could be read from input file. Please check/specify "
                                               "name of SMILES/InChI containing column (--mol_column))")
                else:
                    substructure.error(message="File format not recognized. Please specify the file format (--file_format)")
        else:
            query = {}
            substructure.error(message=f"File {args.input} not found. Please make sure you specified the correct path")
    return query


def substruct(args):
    start_time = time.time()
    print(f"Start: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    if args.filter:
        filter_dict = check_filter(args.filter)
    else:
        filter_dict = {}
    pool = None
    if args.mpi_np:
        if 1 < args.mpi_np <= mp.cpu_count():
            pool = mp.Pool(processes=args.mpi_np)
        else:
            pool = mp.Pool(processes=int(mp.cpu_count()/2))
    ### Load database from vsdb file or from SD file
    print(f"Loading database {args.database} ...")
    mols = {}
    if args.database in db_config:
        try:
            mols = pickle.load(open(f"{config['local_db']}/{args.database}.vsdb", "rb"))
        except FileNotFoundError:
            try:
                mols = pickle.load(open(f"{config['global_db']}/{args.database}.vsdb", "rb"))
            except FileNotFoundError:
                substructure.error(message=f"{args.database} not found. Please make sure you specified the correct shortcut")
    else:
        if os.path.exists(args.database):
            if args.database.endswith(".vsdb"):
                try:
                    mols = pickle.load(open(args.database, "rb"))
                except:
                    substructure.error(message=f"{args.output} could not be opened. Please make sure the file has the correct "
                                               f"format")

            else:
                if args.mpi_np:
                    mols, failed = read.read_sd_mp(args.database, pool)
                else:
                    mols, failed = read.read_db_from_sd(args.database)
                if failed:
                    print(f"{len(failed)} of {len(mols) + len(failed)} molecules in {args.database} could not be processed")
                if not mols:
                    substructure.error(message="No molecules could be read from SD file. Please make sure it has the right "
                                               "format")
        else:
            substructure.error(message=f"File {args.database} not found. Please make sure you specified the correct path")
    db_desc = mols.pop("config")
    ### Read input query
    print("Reading input molecules ...")
    if args.smarts:
        query = read.read_smarts(args.smarts)
        if not query:
            substructure.error(message="No valid molecules could be generated from the provided SMARTS.")
    elif args.smiles:
        query = read.read_smiles(args.smiles, args.mode, args.ntauts)
        if not query:
            substructure.error(message="No valid molecules could be generated from the provided SMILES.")
    else:
        if os.path.exists(args.input):
            query = read.read_file(args.input, args.input_format, args.smiles_column, args.delimiter, args.mode, args.ntauts)
            if not query:
                if args.input.endswith(".sdf") or args.input_format == "sdf":
                    substructure.error(message="No valid molecules could be read from SD file.")
                elif args.input.endswith(".csv") or args.input_format == "csv":
                    substructure.error(message="No valid molecules could be read from input file. Please check/specify "
                                               "name of SMILES/InChI containing column (--mol_column) or check/specify the"
                                               "separator (--delimiter)")
                elif args.input.endswith(".xlsx") or args.input_format == "xlsx":
                    substructure.error(message="No valid molecules could be read from input file. Please check/specify "
                                               "name of SMILES/InChI containing column (--mol_column))")
                else:
                    substructure.error(message="File format not recognized. Please specify the file format (--file_format)")
        else:
            query = {}
            substructure.error(message=f"File {args.input} not found. Please make sure you specified the correct path")
    ### perform substructure search depending on selected mode and parameters
    print("Performing substructure search ...")
    sub_time = time.time()
    results = {}
    if db_desc[0] == "yes":
        if args.mode == "std":
            key = "mol_sta"
        elif args.mode == "can_taut":
            key = "mol_can"
        elif args.mode == "all_tauts":
            key = "mol_sta"
        else:
            key = "mol"
    else:
        key = "mol"
    if args.fullmatch:
        if args.mpi_np:
            if args.mode == "std" or args.mode == "mol_can":
                sss.sss_fm_mp(query, mols, key, filter_dict, results, pool)
            elif args.mode == "no_std":
                sss.sss_fm_nost_mp(query, mols, key, filter_dict, results, pool)
            else:
                sss.sss_fm_taut_mp(query, mols, key, filter_dict, results, pool)
        else:
            if args.mode == "std" or args.mode == "mol_can":
                sss.sss_fm(query, mols, key, filter_dict, results)
            elif args.mode == "no_std":
                sss.sss_fm_nost(query, mols, key, filter_dict, results)
            else:
                sss.sss_fm_taut(query, mols, key, filter_dict, results)
    else:
        if args.mpi_np:
            if args.mode == "std" or args.mode == "can_taut" or args.mode == "no_std":
                sss.sss_mp(query, mols, key, filter_dict, results, pool)
            else:
                sss.sss_mp_taut(query, mols, key, filter_dict, results, pool)
        else:
            if args.mode == "std" or args.mode == "can_taut" or args.mode == "no_std":
                sss.sss(query, mols, key, filter_dict, results)
            else:
                sss.sss_taut(query, mols, key, filter_dict, results)
    sub_time_2 = time.time()
    sub_dur = sub_time_2 - sub_time
    print(sub_dur)
    print("Finished substructure search")
    del mols
    ### calculate properties if desired
    if args.properties:
        for i in results:
            calc_props(results[i]["mol"], results[i]["props"])
    ### write results to output file(s)
    print(f"{len(results)} matches found")
    print("Generating output file(s) ...")
    if args.multfile:
        if results:
            if args.output.endswith(".csv"):
                write_output.gen_csv_xls_mult(query, results, args.output)
            elif args.output.endswith(".xlsx") or args.output.endswith(".xls"):
                write_output.gen_csv_xls_mult(query, results, args.output)
            else:
                write_output.gen_sdf_mult(query, results, args.output)
            if args.PDF:
                print("Generating PDF file(s) ...")
                if args.output.endswith(".sdf") or args.output.endswith(".csv") or args.output.endswith(".xlsx") or args.output.endswith(".xls"):
                    out_file = args.output.rsplit(".", maxsplit=1)[0]
                else:
                    out_file = args.output
                visualize.gen_pdf_mf(query, results, out_file, ttf_path)
    else:
        if results:
            if args.output.endswith(".csv") or args.output.endswith(".xlsx") or args.output.endswith(".xls"):
                write_output.gen_csv_xls(results, args.output)
            else:
                write_output.gen_sdf(results, args.output)
            if args.PDF:
                print("Generating PDF file(s) ...")
                if args.output.endswith(".sdf") or args.output.endswith(".csv") or args.output.endswith(
                        ".xlsx") or args.output.endswith(".xls"):
                    out_file = f"{args.output.rsplit('.', maxsplit=1)[0]}.pdf"
                else:
                    out_file = f"{args.output}.pdf"
                visualize.gen_pdf(query, results, out_file, ttf_path)
    end_time = time.time()
    print(f"Finished: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    duration = round(end_time - start_time, 5)
    if args.mpi_np:
        pool.close()
    print(f"Finished in {duration} seconds")


substructure.set_defaults(func=substruct)


## Fingerprint similarity search

fp_sim = subparsers.add_parser("fpsim", description="molecular similarity search using fingerprints")
group_fp = fp_sim.add_mutually_exclusive_group(required=True)
group_fp.add_argument("-in", "--input", help="specify path of input file [sdf, csv, xlsx]", metavar="")
group_fp.add_argument("-smi", "--smiles", help="specify SMILES string on command line in double quotes",
                      action="append", metavar="")
group_fp.add_argument("-sma", "--smarts", help="specify SMARTS string on command line in double quotes",
                      action="append", metavar="")
fp_sim.add_argument("-db", "--database", help="specify the database file [sdf or vsdb] or specify the shortcut for an "
                                              "integrated database", default=db_default, metavar="")
fp_sim.add_argument("-out", "--output", help="specify name of output file", default="vsflow_fingerprint.sdf", metavar="")
fp_sim.add_argument("-m", "--mode", help="choose a mode for substructure search [std, all_tauts, can_taut, no_std]",
                    choices=["std", "all_tauts", "can_taut", "no_std"], default="std", metavar="")
fp_sim.add_argument("-np", "--mpi_np", type=int, help="Specify the number of processors used when the application is "
                                                      "run in multiprocessing mode.", metavar="")
fp_sim.add_argument("-pdf", "--PDF", help="generate a pdf file for substructure matches", action="store_true")
fp_sim.add_argument("-props", "--properties",
                    help="specifies if calculated molecular properties are written to the output files",
                    action="store_true")
fp_sim.add_argument("-mf", "--multfile", help="generate separate output files for every query molecule",
                    action="store_true")
fp_sim.add_argument("-fi", "--input_format", help="Specify the file typ if no file extension is present in input file "
                                                  "name [sdf, csv, xlsx]", metavar="")
fp_sim.add_argument("-col", "--smiles_column", help="Specify name of smiles column in csv file", default="smiles", metavar="")
fp_sim.add_argument("-del", "--delimiter", help="Specify delimiter of csv file", default=";", metavar="")
fp_sim.add_argument("-head", "--header", help="Specify row of file to be used as column names", type=int, default=1, metavar="")
fp_sim.add_argument("-fp", "--fingerprint", help="specify fingerprint to be used", choices=["rdkit", "ecfp", "fcfp", "ap", "tt", "maccs"],
                    default="fcfp", metavar="")
fp_sim.add_argument("-sim", "--similarity", help="specify fingerprint similarity metric to be used",
                    choices=["tan", "dice", "cos", "sok", "russ", "kulc", "mcco", "tver"], default="tan", metavar="")
fp_sim.add_argument("-fpr", "--radius", help="radius of circular fingerprints ecfp and fcfp", type=int,
                    default=3, metavar="")
fp_sim.add_argument("-nbits", "--NBITS", help="number of bits used to generate ecfp and fcfp fingerprints", type=int,
                    default=4096, metavar="")
fp_sim.add_argument("-top", "--top_hits", type=int, default=10,
                    help="Maximum number of molecules with highest similarity to keep. [Default = 10]", metavar="")
fp_sim.add_argument("-map", "--simmap", help="generates similarity maps for fingerprints in pdf file", action="store_true")
fp_sim.add_argument("-cut", "--cutoff", help="specify cutoff value for similarity coefficient", type=float, metavar="")
fp_sim.add_argument("-filt", "--filter", help="specify property to filter screening results", action="append", metavar="")
fp_sim.add_argument("-nt", "--ntauts", help="maximum number of tautomers to be enumerated", type=int, default=100, metavar="")
fp_sim.add_argument("-c", "--chiral", help="specify if chirality should be considered", action="store_true")
fp_sim.add_argument("-tva", "--tver_alpha", help="specify alpha parameter for Tversky similarity", default=0.5, type=float, metavar="")
fp_sim.add_argument("--tver_beta", help="specify beta parameter for Tversky similarity", default=0.5, type=float, metavar="")


def fingerprint(args):
    start_time = time.time()
    print(f"Start: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    if args.filter:
        filter_dict = check_filter(args.filter)
    else:
        filter_dict = {}
    pool = mp.Pool(processes=args.mpi_np)
    print(f"Loading database {args.database} ...")
    sub_time = time.time()
    mols = read_database(args, pool)
    db_desc = mols.pop("config")
    sub_time_2 = time.time()
    sub_dur = sub_time_2 - sub_time
    print(sub_dur)
    print("Reading query input ...")
    query = read_input(args)
    if db_desc[0] == "yes":
        if args.mode == "std":
            key = "mol_sta"
        elif args.mode == "can_taut":
            key = "mol_can"
        elif args.mode == "all_tauts":
            key = "mol_sta"
        else:
            key = "mol"
    else:
        key = "mol"
    print("Calculating fingerprints ...")
    sub_time = time.time()
    features = False
    if args.fingerprint == "fcfp" or args.fingerprint == "ecfp":
        if args.fingerprint == "fcfp":
            name = f"FCFP{args.radius * 2}-like Morgan {args.NBITS} bits"
            features = True
        else:
            name = f"ECFP{args.radius * 2}-like Morgan {args.NBITS} bits"
        if args.mpi_np:
            argslist = [(mols[i][key], i, args.radius, features, args.chiral, args.NBITS) for i in mols]
            fps = pool.starmap(fpsearch.fp_morgan_mp, argslist)
            fpsearch.set_fp_mp(fps, mols)
            argslist = [(query[j]["mol"], j, args.radius, features, args.chiral, args.NBITS) for j in query]
            fps = pool.starmap(fpsearch.fp_morgan_mp, argslist)
            fpsearch.set_fp_mp(fps, query)
            del argslist
            del fps
        else:
            fpsearch.fp_morgan(mols, key, args.radius, args.NBITS, features, args.chiral)
            fpsearch.fp_morgan(query, "mol", args.radius, args.NBITS, features, args.chiral)
    elif args.fingerprint == "rdkit":
        name = f"RDKit {args.NBITS} bits"
        if args.mpi_np:
            argslist = [(mols[i][key], i, args.NBITS) for i in mols]
            fps = pool.starmap(fpsearch.fp_rdkit_mp, argslist)
            fpsearch.set_fp_mp(fps, mols)
            argslist = [(query[j]["mol"], j, args.NBITS) for j in query]
            fps = pool.starmap(fpsearch.fp_rdkit_mp, argslist)
            fpsearch.set_fp_mp(fps, query)
            del argslist
            del fps
        else:
            fpsearch.fp_rdkit(mols, key, args.NBITS)
            fpsearch.fp_rdkit(query, "mol", args.NBITS)
    elif args.fingerprint == "ap":
        name = f"AtomPairs {args.NBITS} bits"
        if args.mpi_np:
            argslist = [(mols[i][key], i, args.NBITS, args.chiral) for i in mols]
            fps = pool.starmap(fpsearch.fp_atompairs_mp, argslist)
            fpsearch.set_fp_mp(fps, mols)
            argslist = [(query[j]["mol"], j, args.NBITS, args.chiral) for j in query]
            fps = pool.starmap(fpsearch.fp_atompairs_mp, argslist)
            fpsearch.set_fp_mp(fps, query)
            del argslist
            del fps
        else:
            fpsearch.fp_atompairs(mols, key, args.NBITS, args.chiral)
            fpsearch.fp_atompairs(query, "mol", args.NBITS, args.chiral)
    elif args.fingerprint == "tt":
        name = f"TopologicalTorsion {args.NBITS} bits"
        if args.mpi_np:
            argslist = [(mols[i][key], i, args.NBITS, args.chiral) for i in mols]
            fps = pool.starmap(fpsearch.fp_torsion_mp, argslist)
            fpsearch.set_fp_mp(fps, mols)
            argslist = [(query[j]["mol"], j, args.NBITS, args.chiral) for j in query]
            fps = pool.starmap(fpsearch.fp_torsion_mp, argslist)
            fpsearch.set_fp_mp(fps, query)
            del argslist
            del fps
        else:
            fpsearch.fp_torsion(mols, key, args.NBITS, args.chiral)
            fpsearch.fp_torsion(query, "mol", args.NBITS, args.chiral)
    else:
        name = "MACCS"
        if args.mpi_np:
            argslist = [(mols[i][key], i) for i in mols]
            fps = pool.starmap(fpsearch.fp_maccs_mp, argslist)
            fpsearch.set_fp_mp(fps, mols)
            argslist = [(query[j]["mol"], j) for j in query]
            fps = pool.starmap(fpsearch.fp_maccs_mp, argslist)
            fpsearch.set_fp_mp(fps, query)
            del argslist
            del fps
        else:
            fpsearch.fp_maccs(mols, key)
            fpsearch.fp_maccs(query, "mol")
    print("Calculating similarities ...")
    if args.cutoff:
        if args.similarity == "tver":
            results = fpsearch.sim_tver(mols, query, key, args.cutoff, args.similarity, filter_dict, name,
                                        args.tver_alpha, args.tver_beta)
        else:
            results = fpsearch.sim(mols, query, key, args.cutoff, args.similarity, filter_dict, name)
    else:
        if args.similarity == "tver":
            results = fpsearch.sim_top_tver(mols, query, key, args.top_hits, args.similarity, filter_dict, name,
                                            args.tver_alpha, args.tver_beta)
            print(len(results))
        else:
            results = fpsearch.sim_top(mols, query, key, args.top_hits, args.similarity, filter_dict, name)
    sub_time_2 = time.time()
    sub_dur = sub_time_2 - sub_time
    print(sub_dur)
    del mols
    print(f"{len(results)} matches found")
    if args.properties:
        for i in results:
            calc_props(results[i]["mol"], results[i]["props"])
    print("Generating output file(s) ...")
    if args.multfile:
        if results:
            if args.output.endswith(".csv"):
                write_output.gen_csv_xls_mult(query, results, args.output)
            elif args.output.endswith(".xlsx") or args.output.endswith(".xls"):
                write_output.gen_csv_xls_mult(query, results, args.output)
            else:
                write_output.gen_sdf_mult(query, results, args.output)
            if args.PDF:
                print("Generating PDF file ...")
                if args.output.endswith(".sdf") or args.output.endswith(".csv") or args.output.endswith(".xlsx") or args.output.endswith(".xls"):
                    out_file = args.output.rsplit(".", maxsplit=1)[0]
                else:
                    out_file = args.output
                if args.simmap:
                    print(f"Calculating similarity maps for {len(results)} matches ...")
                    visualize.fp_maps(results, query, args.fingerprint, args.radius, args.NBITS, features,
                                      args.similarity, out_file, ttf_path, args.multfile)
                else:
                    visualize.gen_pdf_mf(query, results, out_file, ttf_path)
    else:
        if results:
            if args.output.endswith(".csv") or args.output.endswith(".xlsx") or args.output.endswith(".xls"):
                write_output.gen_csv_xls(results, args.output)
            else:
                write_output.gen_sdf(results, args.output)
            if args.PDF:
                print("Generating PDF file(s) ...")
                if args.output.endswith(".sdf") or args.output.endswith(".csv") or args.output.endswith(
                        ".xlsx") or args.output.endswith(".xls"):
                    out_file = f"{args.output.rsplit('.', maxsplit=1)[0]}.pdf"
                else:
                    out_file = f"{args.output}.pdf"
                if args.simmap:
                    if args.fingerprint == "MACCS":
                        print("For MACCS keys no similarity map can be computed")
                        visualize.gen_pdf(query, results, out_file, ttf_path)
                    elif args.similarity == "tver":
                        print("For Tversky Similarity no similarity map can be computed")
                        visualize.gen_pdf(query, results, out_file, ttf_path)
                    else:
                        print(f"Calculating similarity maps for {len(results)} matches ...")
                        visualize.fp_maps(results, query, args.fingerprint, args.radius, args.NBITS, features, args.similarity, out_file, ttf_path, args.multfile)
                else:
                    visualize.gen_pdf(query, results, out_file, ttf_path)
    end_time = time.time()
    print(f"Finished: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    duration = round(end_time - start_time, 5)

    print(f"Finished in {duration} seconds")

    pool.close()


fp_sim.set_defaults(func=fingerprint)


shape_sim = subparsers.add_parser("shape")
shape_group = shape_sim.add_mutually_exclusive_group(required=True)
shape_group.add_argument("-in", "--input")
shape_group.add_argument("-smi", "--smiles", action="append")
shape_sim.add_argument("-out", "--output", default="shape.sdf")
shape_sim.add_argument("-db", "--database", default=db_default)
shape_sim.add_argument("-np", "--mpi_np", type=int)
shape_sim.add_argument("-py", "--pymol", action="store_true")
shape_sim.add_argument("-top", "--top_hits", default=10, type=int)
shape_sim.add_argument("-am", "--align_method", choices=["mmff", "crippen"], default="mmff")
shape_sim.add_argument("-pdf", "--PDF", action="store_true")



def shape(args):
    mols = {}
    if args.database in db_config:
        try:
            mols = pickle.load(open(f"{config['local_db']}/{args.database}.vsdb", "rb"))
        except FileNotFoundError:
            try:
                mols = pickle.load(open(f"{config['global_db']}/{args.database}.vsdb", "rb"))
            except FileNotFoundError:
                parser.error(
                    message=f"{args.database} not found. Please make sure you specified the correct shortcut")
    else:
        if os.path.exists(args.database):
            if args.database.endswith(".vsdb"):
                try:
                    mols = pickle.load(open(args.database, "rb"))
                except:
                    parser.error(message=f"{args.database} could not be opened. Please make sure the file has the correct "
                                               f"format")
            else:
                parser.error(message="Database must be in format .vsdb. Use mode preparedb to prepare a database for"
                                     " shape similarity screening.")
        else:
            parser.error(message=f"{args.database} could not be opened. Please make sure you specified the correct path")

    db_desc = mols.pop("config")
    if db_desc[1] == 0:
        parser.error(message=f"Database {args.database} does not contain conformers. Use mode preparedb to generate "
                             f"conformers and prepare a database for shape similarity screening")
    seed = db_desc[3]
    num_confs = db_desc[1]

    if args.smiles:

        query = read.read_smiles_std(args.smiles)
    else:
        query = read.read_sd_3d(args.input)

        print(query)
    # perform shape screening with specified parameters
    if args.mpi_np:
        pool_shape = mp.Pool(processes=args.mpi_np)
        query_confs = pool_shape.starmap(shapesearch.gen_query_conf_mp, [(query[i]["mol"], i, num_confs, seed) for i in query])
        print(query_confs)
        for entry in query_confs:
            query[entry[0]]["confs"] = entry[1]
            query[entry[0]]["fp_shape"] = entry[2]
        del query_confs
        print(query)
        aligns = pool_shape.starmap(shapesearch.shape_pfp_tani_mp, [(mols[i]["confs"], i, query[j]["confs"], query[j]["fp_shape"], j)
                                                for i in mols for j in query])
        print(aligns)
        print(max(aligns))
        aligns = sorted(aligns, reverse=True)
        pool_shape.close()
    else:
        shapesearch.gen_query_conf_pfp(query, num_confs, seed)
        aligns = shapesearch.shape_pfp_tani(mols, query)
    # prepare and write results to output files
    grouped = [res[:args.top_hits] for res in (list(group) for k, group in groupby(aligns, lambda x: x[3]))]
    results = {}
    counter = 0
    for entry in grouped:
        for feat in entry:
            mols[feat[4]]["props"]["Combo_Score"] = feat[0]
            mols[feat[4]]["props"]["Shape_Similarity"] = feat[1]
            mols[feat[4]]["props"]["3D_FP_Similarity"] = feat[2]
            mols[feat[4]]["props"]["QuerySmiles"] = query[feat[3]]["pattern"]
            # results[counter] = {"mol": mols[feat[4]]["confs"], "props": mols[feat[4]]["props"], "top_conf": feat[5],
            #                     "q_num": feat[3]}
            results[counter] = {"mol": feat[6], "props": mols[feat[4]]["props"], "top_conf": feat[5],
                                "q_num": feat[3]}
            counter += 1
    print(results)
    out_file = args.output.rsplit(".sdf", maxsplit=1)[0]
    for j in query:
        with open(f"{out_file}_{j + 1}.sdf", "w") as out:
            for i in results:
                if results[i]["q_num"] == j:
                    write_output.write_sdf_conformer(results[i]["mol"], results[i]["props"], results[i]["top_conf"], out)
        with open(f"{out_file}_{j + 1}_query.sdf", "w") as out_q:
            write_output.write_sdf_conformer(query[j]["confs"], {"Smiles": query[j]["pattern"]}, 0, out_q)
    if args.pymol:
        #path = os.path.dirname(os.path.abspath(args.output))
        for j in query:
            #print(f"{out_file}_{j + 1}".split("/")[])
            #pref = f"{out_file}_{j + 1}_query".split("/")[-1]
            visualize.export_pymol(f"{out_file}_{j + 1}_query.sdf", f"{out_file}_{j + 1}.sdf")
    if args.PDF:
        visualize.gen_pdf_shape(query, results, out_file, ttf_path)




shape_sim.set_defaults(func=shape)




canon = subparsers.add_parser("preparedb")
# canon_group = canon.add_mutually_exclusive_group()
canon_group = canon.add_mutually_exclusive_group(required=True)
canon.add_argument("-in", "--input", required=True)
canon.add_argument("-np", "--mpi", type=int)
canon_group.add_argument("-out", "--output")
canon_group.add_argument("-int", "--integrate", help="specify shortcut for database")
canon.add_argument("-intg", "--int_global", help="Stores database by default within the folder of the script", action="store_true")
canon.add_argument("-nt", "--ntauts", help="maximum number of tautomers to be enumerated during standardization process", type=int, default=100)
canon.add_argument("-st", "--standardize", help="standardizes molecules, removes salts and associated charges", action="store_true")
canon.add_argument("-c", "--conformers", help="generates multiple 3D conformers, required for mode shape", action="store_true")
canon.add_argument("-nc", "--nconfs", help="number of conformers generated", type=int, default=20)
canon.add_argument("-t", "--threshold", help="Retain only the conformations out of nconfs that are at least "
                                             "this far apart from each other (RMSD calculated on heavy atoms)", type=float)
# canon_group.add_argument("-can", "--canonicalize", help="standardizes molecules, removes salts and associated charges and canonicalizes tautomers", action="store_true")
#canon_group.add_argument("-ats", "--all_tauts", help="generate all possible tautomers for molecule up to a maximum specified in --ntauts", action="store_true")




def do_standard(mol, dict, ntauts):
    mol_dict = {}
    try:
        mol_sta = Standardizer().charge_parent(Standardizer().fragment_parent(mol), skip_standardize=True)
        mol_can = TautomerCanonicalizer(max_tautomers=ntauts).canonicalize(mol_sta)
        # mol_dict["mol_sta"] = mol_sta
        # mol_dict["mol_can"] = mol_can
        # mol_dict["props"] = dict
        # mol_dict["mol"] = mol
        mol_dict["mol_sta"] = Chem.MolFromSmiles(Chem.MolToSmiles(mol_sta))
        mol_dict["mol_can"] = Chem.MolFromSmiles(Chem.MolToSmiles(mol_can))
        mol_dict["props"] = dict
        mol_dict["mol"] = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        # return (mol_can, dict)
        return mol_dict
    except:
        # return (mol, dict)
        mol_dict["mol_sta"] = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        mol_dict["mol_can"] = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        mol_dict["props"] = dict
        mol_dict["mol"] = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        # mol_dict["mol_sta"] = mol
        # mol_dict["mol_can"] = mol
        # mol_dict["props"] = dict
        # mol_dict["mol"] = mol
        return mol_dict
        # return ([mol], dict)


def gen_confs_mmff(mol, num, nconfs, seed):
    params = Chem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    #params.pruneRmsThresh = threshold
    params.numThreads = 0
    params.randomSeed = seed
    mol_H = Chem.AddHs(mol)
    Chem.EmbedMultipleConfs(mol_H, numConfs=nconfs, params=params)
    #Chem.EmbedMultipleConfs(mol_H, numConfs=nconfs, randomSeed=seed, ETversion=2, numThreads=0)
    try:
        Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=0, maxIters=2000)
    except:
        pass
    mol_3D = Chem.RemoveHs(mol_H)
    print(num)
    return (mol_3D, num)



def canon_mol(args):
    start_time = time.time()
    print(f"Start: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")

    if args.integrate:
        db_name = args.integrate
        if args.int_global:
            db_path = config["global_db"]
            try:
                if os.path.exists(db_path):
                    with open(f"{db_path}/.test", "w") as test_file:
                        test_file.write("")
                else:
                    os.mkdir(db_path)
            except PermissionError:
                canon.error(message="You do not have the permission to integrate a database globally. Please contact "
                                    "your system administrator or re-run with sudo permissions !")
        else:
            db_path = config["local_db"]
            try:
                if os.path.exists(db_path):
                    pass
                else:
                    os.mkdir(db_path)
            except FileNotFoundError:
                canon.error(message="Path not valid. Please make sure you specified a correct path !")
            except PermissionError:
                canon.error(message="Permission denied")
            except OSError:
                canon.error(message="Permission denied")
        if args.integrate in db_config:
            choice = input(
                f"A database with name {db_name} is already integrated in VSFlow. Press 'o' to override the database, "
                f"press 'n' to choose a new name for the database or press 'c' to cancel and exit.")
            while choice != "o" and choice != "n" and choice != "c":
                choice = input("Press 'o' to override the database, press 'n' to choose a new name for the "
                                "database or press 'c' to cancel and exit.")
            if choice == "o":
                try:
                    test_load = pickle.load(open(f"{config['global_db']}/{args.integrate}", "rb"))
                    del test_load
                    with open(f"{config['global_db']}/.test", "w") as test_file:
                        test_file.write("")
                except FileNotFoundError:
                    pass
                except PermissionError:
                    sec_choice = input(f"You do not have the permission to change the database {args.integrate}. Please contact "
                          f"your system administrator or run again with sudo. Press 'c' to cancel or press 'n' to "
                          f"enter a different name")
                    while sec_choice != "n" and sec_choice != "c":
                        sec_choice = input(f"Press 'c' to cancel or press 'n' to enter a different name")
                    if sec_choice == "n":
                        db_name = input("Please enter different name:")
                        while db_name == args.integrate:
                            db_name = input("Please enter different name:")
                    else:
                        exit()
            elif choice == "n":
                db_name = input("Please enter different name:")
                while db_name == args.integrate:
                    db_name = input("Please enter different name:")
            else:
                exit()
        out_path = f"{db_path}/{db_name}.vsdb"
    else:
        if args.output.endswith(".vsdb"):
            out_path = args.output
        else:
            out_path = f"{args.output}.vsdb"
    standardized = "no"
    conformers = 0
    seed = None
    if args.standardize:
        standardized = "yes"
        if args.mpi:
            can_pool = mp.Pool(processes=args.mpi)
            data, failed = read.read_sd_mp(args.input, can_pool)
            print(f"Finished reading sdf file: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
            if failed:
                print(f"{len(failed)} molecules out of {len(data)} could not be processed")
            data_can = can_pool.starmap(do_standard, [(data[n]["mol"], data[n]["props"], args.ntauts) for n in data])
            mols = {}
            for i in range(len(data_can)):
                mols[i] = data_can[i]
            can_pool.close()
        else:
            data, failed = read.read_db_from_sd(args.input)
            print(f"Finished reading sdf file: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
            if failed:
                print(f"{len(failed)} molecules out of {len(data)} could not be processed")

            mols = {}
            for n in data:
                std = do_standard(data[n]["mol"], data[n]["props"], args.ntauts)
                mols[n] = std
        print(f"Finished standardizing molecules: {time.strftime('%m/%d/%Y, %H:%M:%S', time.localtime())}")
    else:
        if args.mpi:
            can_pool = mp.Pool(processes=args.mpi)
            mols, failed = read.read_sd_mp(args.input, can_pool)
            can_pool.close()
        else:
            mols, failed = read.read_db_from_sd(args.input)
        if failed:
            print(f"{len(failed)} molecules out of {len(mols) + len(failed)} could not be processed")
    if args.conformers:
        conformers = args.nconfs
        for i in mols:
            try:
                mol = mols[i]["mol_sta"]
                key = "mol_sta"
            except KeyError:
                key = "mol"
            break
        # if args.threshold:
        #     threshold = args.threshold
        # else:
        #     threshold = -1.0
        print(key)
        seed = random.randint(0, 10000)
        if args.mpi:
            print("mpi")
            can_pool = mp.Pool(processes=args.mpi)
            confs = can_pool.starmap(gen_confs_mmff, [(mols[i][key], i, args.nconfs, seed) for i in mols])
            for entry in confs:
                mols[entry[1]]["confs"] = entry[0]
            can_pool.close()
        else:
            params = Chem.ETKDGv3()
            params.useSmallRingTorsions = True
            params.useMacrocycleTorsions = True
            #params.pruneRmsThresh = threshold
            params.numThreads = 0
            params.randomSeed = seed
            counter = 0
            for i in mols:
                mol_H = Chem.AddHs(mols[i][key])
                Chem.EmbedMultipleConfs(mol_H, numConfs=100, params=params)
                #Chem.EmbedMultipleConfs(mol_H, numConfs=args.nconfs, randomSeed=seed, ETversion=2, numThreads=0)
                try:
                    Chem.MMFFOptimizeMoleculeConfs(mol_H, numThreads=0, maxIters=2000)
                except:
                    pass
                #mol_3D = Chem.RemoveHs(mol_H)
                mols[i]["confs"] = mol_H
                counter += 1
                print(counter)
    mols["config"] = [standardized, conformers, len(mols), seed]

    #pickle.dump(mols, open(f"{db_path}/{db_name}.vsdb", "wb"))
    pickle.dump(mols, open(out_path, "wb"))
    if args.integrate:
        db_config[db_name] = [time.ctime(os.path.getmtime(out_path)),
                              standardized,
                              conformers,
                              len(mols),
                              seed]
        pickle.dump(db_config, open(f"{home}/.vsflow/.db_config", "wb"))
        print(f"{args.input} was integrated as database {db_name} in VSFlow. You can now search the database calling -db "
              f"{db_name}.")

    end_time = time.time()
    duration = round(end_time - start_time)
    print(f"Finished in {duration} seconds")



canon.set_defaults(func=canon_mol)




show_db = subparsers.add_parser("managedb")
show_db.add_argument("-d", "--default", help="specify name of database to be set as default")
show_db.add_argument("-s", "--show", help="Show currently integrated databases in VSFlow", action="store_true")
show_db.add_argument("-rm", "--remove", help="specify name of database to be removed")


def get_db(args):
    db_shortcut = ["shortcut"]
    db_create = ["created"]
    db_standard = ["standardized"]
    db_conformers = ["conformers"]
    db_length = ["number of cpds"]
    for db in db_config:
        db_shortcut.append(db)
        db_create.append(db_config[db][0])
        db_standard.append(db_config[db][1])
        db_conformers.append(str(db_config[db][2]))
        db_length.append(db_config[db][3])
    default_db = pickle.load(open(f"{home}/.vsflow/.db_default", "rb"))
    if args.default:
        if args.default in db_shortcut:
            pickle.dump(args.default, open(f"{home}/.vsflow/.db_default", "wb"))
            print(f"'{args.default}' has been set as default database. This database is now used if the -db flag is not"
                  f" specified.")
        else:
            print(f"Database '{args.default}' is not integrated in VSFlow. Use preparedb -h to see how a database can "
                  f"be integrated.")
    if args.show:
        for i in range(len(db_shortcut)):
            if i == 0:
                print('\033[1m')
                print("DATABASES\n")
                print(f"{db_shortcut[i]}" + " "*(max([len(string) for string in db_shortcut]) + 5 - len(db_shortcut[i])) +
                      f"{db_create[i]}" + " "*(max([len(string) for string in db_create]) + 5 - len(db_create[i])) +
                      #f"{db_source[i]}" + " " * (max([len(string) for string in db_source]) + 5 - len(db_source[i])) +
                      f"{db_standard[i]}" + " " * (max([len(string) for string in db_standard]) + 5 - len(db_standard[i])) +
                      f"{db_conformers[i]}" + " " * (max([len(string) for string in db_conformers]) + 5 - len(db_conformers[i])) +
                      #f"{db_all_tauts[i]}" + " " * (max([len(string) for string in db_all_tauts]) + 5 - len(db_all_tauts[i])) +
                      f"{db_length[i]}")
                print('\033[0m')
            else:
                print(f"{db_shortcut[i]}" + " "*(max([len(string) for string in db_shortcut]) + 5 - len(db_shortcut[i])) +
                      f"{db_create[i]}" + " "*(max([len(string) for string in db_create]) + 5 - len(db_create[i])) +
                      #f"{db_source[i]}" + " " * (max([len(string) for string in db_source]) + 5 - len(db_source[i])) +
                      f"{db_standard[i]}" + " " * (max([len(string) for string in db_standard]) + 5 - len(db_standard[i])) +
                      f"{db_conformers[i]}" + " " * (max([len(string) for string in db_conformers]) + 5 - len(db_conformers[i])) +
                      #f"{db_all_tauts[i]}" + " " * (max([len(string) for string in db_all_tauts]) + 5 - len(db_all_tauts[i])) +
                      f"{db_length[i]}")

        print("\n")
        print('\033[1m')
        print(f"Default database: {default_db}")
        print('\033[0m')
        print("You can set (or change) a default database by calling managedb --default {shortcut}")
    if args.remove:
        if args.remove in db_config:
            choice = input(f"Are you sure you want to remove {args.remove} ? [y/n]")
            while choice != "y" and choice != "n":
                choice = input(f"Are you sure you want to remove {args.remove} ? [y/n]")
            if choice == "n":
                pass
            else:
                try:
                    os.remove(f"{db_config['local_db']}/{args.remove}")
                    db_config.pop(args.remove)
                    pickle.dump(db_config, open(f"{home}/.vsflow/.db_config", "wb"))
                except FileNotFoundError:
                    try:
                        os.remove(f"{db_config['global_db']}/{args.remove}")
                        db_config.pop(args.remove)
                        pickle.dump(db_config, open(f"{home}/.vsflow/.db_config", "wb"))
                        print(f"{args.remove} was successfully removed !")
                    except PermissionError:
                        print(f"You do not have the permission to remove {args.remove}. Please contact your system"
                              f"administrator!")
        else:
            print(f"No database with name {args.remove} integrated in VSFlow !")



show_db.set_defaults(func=get_db)


def main():
    args = parser.parse_args()
    if "func" in args:
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

