import os
from fpdf import FPDF
from pdfrw import PdfReader
from pdfrw import PdfWriter
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import SimilarityMaps
from pymol import cmd
import fpsearch


def add_colours_to_map(els, cols, col_num, COLS):
    for el in els:
        if el not in cols:
            cols[el] = []
        if COLS[col_num] not in cols[el]:
            cols[el].append(COLS[col_num])


def write_props(pdf, props_list, txt_place_y, txt_space_y, size):
    pdf.set_font("DejaVu", "", size)
    for tag in props_list:
        string = f"{tag[0]}: {tag[1]}"
        str_fact = 0
        str_space = 0
        one_line = int((size * -10 + 190) / 2)
        two_lines = size * -10 + 190
        space = size / 2
        if len(string) < one_line:
            pdf.text(103, txt_place_y + txt_space_y, string)
            txt_space_y += space
        else:
            if len(string) <= two_lines:
                while str_fact < len(string) - one_line:
                    pdf.text(103, txt_place_y + txt_space_y + str_space, string[str_fact:str_fact + one_line])
                    str_space += space
                    str_fact += one_line
                if str_fact >= len(string) - one_line:
                    pdf.text(103, txt_place_y + txt_space_y + str_space, string[str_fact:len(string)])
                    str_space += space
                txt_space_y = txt_space_y + str_space
            else:
                write_string = string[:two_lines - 3] + "..."
                while str_fact < len(write_string) - one_line:
                    pdf.text(103, txt_place_y + txt_space_y + str_space, write_string[str_fact:str_fact + one_line])
                    str_space += space
                    str_fact += one_line
                if str_fact >= len(write_string) - one_line:
                    pdf.text(103, txt_place_y + txt_space_y + str_space, write_string[str_fact:len(write_string)])
                    str_space += space
                txt_space_y = txt_space_y + str_space


def gen_pdf_pages(mol_keys, results, ttf_path):
    sort_pages = [mol_keys[i * 3:(i + 1) * 3] for i in range((len(mol_keys) + 3 - 1) // 3)]
    counter = 1
    pages = []
    for page in sort_pages:
        pdf = FPDF()
        pdf.add_font("DejaVu", "", ttf_path, uni=True)
        pdf.add_page()
        img_place_y = 10
        txt_place_y = 15
        num_place_y = 98
        for mol_n in page:
            pdf.set_font("DejaVu", "", 10)
            pdf.image(f"{mol_n}.png", 10, img_place_y, 90)
            os.remove(f"{mol_n}.png")
            pdf.rect(10, img_place_y, 190, 90)
            pdf.dashed_line(100, img_place_y, 100, img_place_y + 90)
            pdf.text(12, num_place_y, str(counter))
            counter += 1
            img_place_y += 94
            num_place_y += 94
            txt_space_y = 0
            line_counter = 0
            props_list = list(results[mol_n]["props"].items())
            for tag in props_list:
                entry_length = len(str(tag[0])) + len(str(tag[1])) + 2
                if entry_length <= 45:
                    line_counter += 1
                else:
                    line_counter += 2
            if line_counter <= 17:
                size = 10
                write_props(pdf, props_list, txt_place_y, txt_space_y, size)
            elif 17 < line_counter <= 19:
                size = 9
                write_props(pdf, props_list, txt_place_y, txt_space_y, size)
            elif 19 < line_counter <= 21:
                size = 8
                write_props(pdf, props_list, txt_place_y, txt_space_y, size)
            elif 21 < line_counter <= 23:
                size = 7
                write_props(pdf, props_list, txt_place_y, txt_space_y, size)
            elif 23 < line_counter <= 25:
                size = 6
                write_props(pdf, props_list, txt_place_y, txt_space_y, size)
            elif 25 < line_counter <= 27:
                size = 5
                write_props(pdf, props_list, txt_place_y, txt_space_y, size)
            else:
                size = 5
                write_props_list = props_list[:27]
                write_props(pdf, write_props_list, txt_place_y, txt_space_y, size)
            txt_place_y += 94
        out = f"{counter}.pdf"
        pdf.output(out)
        pages.append(f"{counter}.pdf")
    return pages


def export_pdf(pages, out_file):
    pdf_writer = PdfWriter()
    for page in pages:
        pdf_reader = PdfReader(page)
        pdf_writer.addPage(pdf_reader.getPage(0))
        os.remove(page)
    with open(out_file, "wb") as file:
        pdf_writer.write(file)


def gen_pdf_mf(query, results, out_file, ttf_path):
    COLOR = [(0.8, 0.0, 0.0)]
    file_counter = 1
    for m in query:
        mol_keys = []
        for n in results:
            if results[n]["q_num"] == m:
                Chem.Compute2DCoords(results[n]["mol"])
                acols = {}
                bcols = {}
                h_rads = {}
                h_lw_mult = {}
                alist = []
                blist = []
                try:
                    for match in results[n]["match"]:
                        alist.extend(match)
                    for bond in query[m]["mol"].GetBonds():
                        aid1 = alist[bond.GetBeginAtomIdx()]
                        aid2 = alist[bond.GetEndAtomIdx()]
                        blist.append(results[n]["mol"].GetBondBetweenAtoms(aid1, aid2).GetIdx())
                    add_bonds = []
                    for ha1 in alist:
                        for ha2 in alist:
                            if ha1 > ha2:
                                b = results[n]["mol"].GetBondBetweenAtoms(ha1, ha2)
                                if b:
                                    add_bonds.append(b.GetIdx())
                    if len(add_bonds) % len(blist) == 0:
                        for bo2 in add_bonds:
                            if bo2 not in blist:
                                blist.append(bo2)
                    add_colours_to_map(alist, acols, 0, COLOR)
                    add_colours_to_map(blist, bcols, 0, COLOR)
                except:
                    pass
                d = rdMolDraw2D.MolDraw2DCairo(600, 600)
                d.DrawMoleculeWithHighlights(results[n]["mol"], "", acols, bcols, h_rads, h_lw_mult, -1)
                d.FinishDrawing()
                d.WriteDrawingText(f"{n}.png")
                mol_keys.append(n)
        pages = gen_pdf_pages(mol_keys, results, ttf_path)
        export_pdf(pages, f"{out_file}_{file_counter}.pdf")
        file_counter += 1


def gen_pdf(query, results, out_file, ttf_path):
    for i in results:
        COLS = [(0.8, 0.0, 0.0), (0.0, 0.8, 0.0),
                (0.0, 0.0, 0.8), (1.0, 0.55, 1.0)]
        acols = {}
        bcols = {}
        h_rads = {}
        h_lw_mult = {}
        Chem.Compute2DCoords(results[i]["mol"])
        alist = []
        blist = []
        try:
            for match in results[i]["match"]:
                alist.extend(match)
            for bond in query[results[i]["q_num"]]["mol"].GetBonds():
                aid1 = alist[bond.GetBeginAtomIdx()]
                aid2 = alist[bond.GetEndAtomIdx()]
                blist.append(results[i]["mol"].GetBondBetweenAtoms(aid1, aid2).GetIdx())
            add_bonds = []
            for ha1 in alist:
                for ha2 in alist:
                    if ha1 > ha2:
                        b = results[i]["mol"].GetBondBetweenAtoms(ha1, ha2)
                        if b:
                            add_bonds.append(b.GetIdx())
            if len(add_bonds) % len(blist) == 0:
                for bo2 in add_bonds:
                    if bo2 not in blist:
                        blist.append(bo2)
            add_colours_to_map(alist, acols, 0, COLS)
            add_colours_to_map(blist, bcols, 0, COLS)
        except:
            pass
        d = rdMolDraw2D.MolDraw2DCairo(600, 600)
        d.DrawMoleculeWithHighlights(results[i]["mol"], "", acols, bcols, h_rads, h_lw_mult, -1)
        d.FinishDrawing()
        d.WriteDrawingText(f"{i}.png")
    pages = gen_pdf_pages(list(results.keys()), results, ttf_path)
    export_pdf(pages, out_file)


def gen_pdf_shape(query, results, out_file, ttf_path):
    for m in query:
        mol_keys = []
        for n in results:
            if results[n]["q_num"] == m:
                Chem.Compute2DCoords(results[n]["mol"])
                acols = {}
                bcols = {}
                h_rads = {}
                h_lw_mult = {}
                d = rdMolDraw2D.MolDraw2DCairo(600, 600)
                d.DrawMoleculeWithHighlights(results[n]["mol"], "", acols, bcols, h_rads, h_lw_mult, -1)
                d.FinishDrawing()
                d.WriteDrawingText(f"{n}.png")
                mol_keys.append(n)
        pages = gen_pdf_pages(mol_keys, results, ttf_path)
        export_pdf(pages, f"{out_file}_{m + 1}.pdf")


def sim_map(results, query, fp_func, metric, out_file, ttf_path):
    for i in results:
        fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(query[results[i]["q_num"]]["mol"], results[i]["mol"],
                                                                       fp_func, metric=fpsearch.sim_dict[metric])
        fig.set_figwidth(3.255)
        fig.set_figheight(3.255)
        filename = f"{i}.png"
        fig.savefig(filename, bbox_inches="tight")
    pages = gen_pdf_pages(list(results.keys()), results, ttf_path)
    export_pdf(pages, out_file)


def sim_map_mf(results, query, fp_func, metric, out_file, ttf_path):
    file_counter = 1
    for m in query:
        mol_keys = []
        for i in results:
            if results[i]["q_num"] == m:
                fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(query[results[i]["q_num"]]["mol"],
                                                                               results[i]["mol"],
                                                                               fp_func,
                                                                               metric=fpsearch.sim_dict[metric])
                fig.set_figwidth(3.255)
                fig.set_figheight(3.255)
                filename = f"{i}.png"
                fig.savefig(filename, bbox_inches="tight")
                mol_keys.append(i)
        pages = gen_pdf_pages(mol_keys, results, ttf_path)
        export_pdf(pages, f"{out_file}_{file_counter}.pdf")
        file_counter += 1


def fp_maps(results, query, fingerprint, fpradius, nbits, features, metric, out_file, ttf_path, multfile):
    if fingerprint == "ecfp" or fingerprint == "fcfp":
        fp_func = lambda m, idx: SimilarityMaps.GetMorganFingerprint(m, atomId=idx, radius=fpradius,
                                                                     fpType='bv', nBits=nbits,
                                                                     useFeatures=features)
    elif fingerprint == "rdkit":
        fp_func = lambda m, idx: SimilarityMaps.GetRDKFingerprint(m, atomId=idx, fpType="bv", nBits=nbits)
    elif fingerprint == "ap":
        fp_func = lambda m, idx: SimilarityMaps.GetAPFingerprint(m, atomId=idx, fpType="bv", nBits=nbits)
    else:
        #fingerprint == "tt":
        fp_func = lambda m, idx: SimilarityMaps.GetTTFingerprint(m, atomId=idx, fpType="bv", nBits=nbits)
    if multfile:
        sim_map_mf(results, query, fp_func, metric, out_file, ttf_path)
    else:
        sim_map(results, query, fp_func, metric, out_file, ttf_path)


def export_pymol(file1, file2):
    py_object1 = file1.rsplit(".sdf", maxsplit=1)[0]
    py_object2 = file2.rsplit(".sdf", maxsplit=1)[0]
    pref1 = py_object1.split("/")[-1]
    pref2 = py_object2.split("/")[-1]
    cmd.load(filename=file1)
    cmd.set_name(pref1, "query_conf")
    cmd.load(filename=file2)
    cmd.split_states(object="query_conf")
    cmd.split_states(object=pref2)
    cmd.delete("query_conf")
    cmd.delete(pref2)
    cmd.save(f"{py_object2}.pse")


