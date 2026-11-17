#!/usr/bin/env python3
#
# VocProfile - 01_report_compounds.py
#
# Compound table with structures and the overview of classes: molecular weight, logP
# and rotatable bonds from the SMILES, one drawing per compound, the table as Excel
# and HTML, the stacked bar of compounds per class and subclass, and two structure
# grids (Fig. 1, Fig. 2a-b, Table S6).
#
# in   COMPOUND_FILE
# out  results/01_report_compounds/compound_report/ (table)
#      figures/01_report_compounds/compound_report/ (structures, grids, bar, HTML)
#      results/01_report_compounds/01_qc.txt
#
# Usage:
#   python scripts/01_report_compounds.py <config.txt>
#
# Notes: docs/NOTES.md, step 01
#
import base64
import math
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from openpyxl import Workbook
from openpyxl.drawing.image import Image
from PIL import Image as PILImage, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import MolDrawOptions, MolToImage

from lib_config import load_config, input_path, output_dir, die
from lib_run import start_step, finish_step, save_figure, note, notes, use_font

cfg, step, args = load_config("01")
start_step(cfg, 'compounds')
LEVELS = cfg["significance_levels"]
TREATMENT = (
    'compounds without a structure (no SMILES) stay in the table and are not drawn',
)
TABLE_PX = step["table_structure_px"]                       # structure images of the table (Excel, HTML)
FIG1_PX = step["fig1_structure_px"]  # structure images drawn in all_compounds_structures.png
FIT_LABELS = step["fig1_fit_labels"]
SUBMITTED_LABEL_SIZE = 7.5           # label size of the submitted figure; it fixes the layout
GRID_SCALE = step["grid_scale"]
GRID_NAME_FACTOR = step["grid_name_factor"]


def stored_layouts():
    """Molecules with stored 2D coordinates, by compound name, from STRUCTURE_LAYOUT_FILE."""
    if not cfg["settings"].flat["STRUCTURE_LAYOUT_FILE"]:
        note("no STRUCTURE_LAYOUT_FILE: every structure is laid out from its SMILES")
        return {}
    path = cfg["settings"].path("STRUCTURE_LAYOUT_FILE")
    if not os.path.isfile(path):
        die(f"STRUCTURE_LAYOUT_FILE: file not found: {path}")
    layouts = {}
    for molecule in Chem.SDMolSupplier(path, removeHs=False):
        if molecule is not None:
            layouts[molecule.GetProp("_Name")] = molecule
    return layouts


def add_drawings(frame):
    """Column 'Drawing': the stored layout of a compound, or the molecule built from its SMILES."""
    shown_as = cfg["labels"]["all_display_names"]
    by_shown_name = {shown_as.get(name, name): molecule for name, molecule in LAYOUTS.items()}
    frame['Drawing'] = [LAYOUTS.get(name, by_shown_name.get(name, molecule))
                        for name, molecule in zip(frame['Compound'], frame['Molecule'])]
    from_smiles = [name for name, molecule in zip(frame['Compound'], frame['Molecule'])
                   if molecule is not None and name not in LAYOUTS and name not in by_shown_name]
    if LAYOUTS and not any("stored layout" in line for line in notes()):
        note(f"structures drawn from the stored layout: {int(frame['Molecule'].notna().sum()) - len(from_smiles)}; "
             f"laid out from SMILES: {from_smiles or 'none'}")


def scaled_options(factor):
    """Drawing options with line width and font limits, which are in pixels, multiplied by a factor."""
    options = MolDrawOptions()
    options.bondLineWidth = options.bondLineWidth * factor
    options.minFontSize = int(options.minFontSize * factor)
    options.maxFontSize = int(options.maxFontSize * factor)
    return options


def fig1_structure(row):
    """Structure for Fig. 1: the table image, or the same drawing rendered at a higher resolution."""
    if FIG1_PX == TABLE_PX:
        return PILImage.open(row['Image_File'])
    factor = FIG1_PX / TABLE_PX
    options = scaled_options(factor)
    if step["fig1_atom_label_min_px"]:
        options.minFontSize = int(step["fig1_atom_label_min_px"] * factor)
    return MolToImage(row['Drawing'], size=(FIG1_PX, FIG1_PX), options=options)


def two_lines(name, width_of, limit):
    """Name on one line if it fits, otherwise broken once after a space, hyphen or comma."""
    if width_of(name) <= limit:
        return name, True
    best = None
    for i, char in enumerate(name[:-1]):
        if char not in ' -,':
            continue
        first = name[:i] if char == ' ' else name[:i + 1]
        candidate = first + '\n' + name[i + 1:]
        width = width_of(candidate)
        if best is None or width < best[1]:
            best = (candidate, width)
    if best is None:
        return name, False
    return best[0], best[1] <= limit


def fit_labels(fig, axes, labels):
    """Enlarge the labels without moving anything else; return the names that do not fit their cell."""
    renderer = fig.canvas.get_renderer()
    limit = (axes[1].get_window_extent(renderer).x0 - axes[0].get_window_extent(renderer).x0) * 0.98
    too_wide = []
    for name, text in labels:
        text.set_in_layout(False)
        text.set_wrap(False)
        text.set_fontsize(step["structure_label_size"])
        text.set_y(step["structure_label_offset"])
        text.set_linespacing(1.0)

        def width_of(candidate):
            text.set_text(candidate)
            return text.get_window_extent(renderer).width

        wrapped, fits = two_lines(name, width_of, limit)
        # a name that is still too wide on two lines is drawn smaller until it fits
        while not fits and text.get_fontsize() > step["structure_label_min_size"]:
            text.set_fontsize(text.get_fontsize() - 0.5)
            wrapped, fits = two_lines(name, width_of, limit)
        text.set_text(wrapped)
        if text.get_fontsize() < step["structure_label_size"]:
            print(f"Label drawn at {text.get_fontsize():g} pt to fit its cell: {name}")
        if not fits:
            too_wide.append(name)
    return too_wide

file_path = input_path(cfg, step["input"])
output_folder = output_dir(cfg, "01_report_compounds", step["output_folder"])
img_folder = os.path.join(output_folder,"compound_structures")
os.makedirs(img_folder, exist_ok=True)

df = pd.read_excel(file_path) if file_path.lower().endswith('.xlsx') else pd.read_csv(file_path, sep='\t')

LAYOUTS = stored_layouts()
df['Molecule'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(str(x)) if pd.notna(x) else None)
add_drawings(df)

df['Molecular_Weight'] = df['Molecule'].apply(lambda x: Descriptors.MolWt(x) if x else None)
df['XLogP'] = df['Molecule'].apply(lambda x: Descriptors.MolLogP(x) if x else None)
df['Num_Rotatable_Bonds'] = df['Molecule'].apply(lambda x: Descriptors.NumRotatableBonds(x) if x else None)

df['Image_File'] = ""
for idx,row in df.iterrows():
    mol = row['Drawing']
    if mol:
        img = MolToImage(mol,size=(TABLE_PX,TABLE_PX))
        filename = re.sub(r'[^\w\-]', '_', str(row['Compound'])) + ".png"
        img_path = os.path.join(img_folder,filename)
        img.save(img_path)
        df.at[idx,'Image_File'] = img_path

excel_path = os.path.join(output_folder,"complete_compound_table_with_images.xlsx")
wb = Workbook()
ws = wb.active
ws.title = "Volatile Compounds"

columns = ["Compound","SMILES","IUPAC_Name","Main_Class","Subclass",
           "Molecular_Weight","XLogP","Num_Rotatable_Bonds","Structure"]
ws.append(columns)

for idx,row in df.iterrows():
    ws.append([
        row['Compound'],
        row['SMILES'],
        row.get('IUPAC_Name',''),
        row['Main_Class'],
        row['Subclass'],
        row['Molecular_Weight'],
        row['XLogP'],
        row['Num_Rotatable_Bonds'],
        ""
    ])

for idx,row in enumerate(df.itertuples(), start=2):
    if row.Image_File:
        img = Image(row.Image_File)
        img.width = 40
        img.height = 40
        ws.add_image(img,f"I{idx}")

wb.save(excel_path)
print(f"Excel report with images saved: {excel_path}")

main_classes = df['Main_Class'].unique()
Subclasss = df['Subclass'].unique()
colors = sns.color_palette('tab20',n_colors=len(Subclasss))

plt.figure(figsize=(12,6))
bottom_vals = [0]*len(main_classes)
for i, sg in enumerate(Subclasss):
    counts = [len(df[(df['Main_Class']==mc) & (df['Subclass']==sg)]) for mc in main_classes]
    plt.bar(main_classes, counts, bottom=bottom_vals, color=colors[i], edgecolor='black', label=sg)
    for j, val in enumerate(counts):
        if val>0:
            plt.text(j, bottom_vals[j]+val/2, str(val), ha='center', va='center', fontsize=9, color='white', fontweight='bold')
    bottom_vals = [sum(x) for x in zip(bottom_vals, counts)]

plt.ylabel("Number of Compounds")
plt.title("Compounds per Main Class and Subclass")
plt.xticks(rotation=0, ha='center')
plt.legend(title="Subclass", loc='upper left', bbox_to_anchor=(0.01,0.99))
plt.tight_layout()
save_figure(os.path.join(output_folder,"stacked_bar_classes_with_counts.png"))
plt.close()
print("Stacked bar plot saved.")

html_path = os.path.join(output_folder,"complete_compound_table.html")
with open(html_path,"w",encoding="utf-8") as f:
    f.write("<html><head><meta charset='utf-8'><title>Volatile Compounds</title></head><body>")
    f.write("<h1>Volatile Compounds Table</h1>")
    f.write("<table border='1' style='border-collapse: collapse; text-align:center;'>")
    f.write("<tr>")
    for col in columns[:-1]:
        f.write(f"<th>{col}</th>")
    f.write("<th>Structure</th></tr>")

    for idx,row in df.iterrows():
        f.write("<tr>")
        for col in columns[:-1]:
            val = row.get(col,'')
            if isinstance(val,float):
                f.write(f"<td>{val:.2f}</td>")
            else:
                f.write(f"<td>{val}</td>")
        if row['Image_File']:
            with open(row['Image_File'],'rb') as img_file:
                b64_img = base64.b64encode(img_file.read()).decode('utf-8')
            f.write(f"<td><img src='data:image/png;base64,{b64_img}' width='80' height='80'/></td>")
        else:
            f.write("<td></td>")
        f.write("</tr>")

    f.write("</table></body></html>")
print(f"HTML report saved: {html_path}")

plt.rcParams['font.family'] = use_font(step["font_family"])
plt.rcParams['font.size'] = 8

cell_width = 150
cell_height = 200
cols_per_row = 6
rows = (len(df) + cols_per_row - 1) // cols_per_row

fig_width = cols_per_row * 2.5
fig_height = rows * 3.0
fig, axes = plt.subplots(rows, cols_per_row, figsize=(fig_width, fig_height))

if rows == 1 and cols_per_row == 1:
    axes = [axes]
elif rows == 1:
    axes = list(axes)
elif cols_per_row == 1:
    axes = list(axes)
else:
    axes = axes.flatten()

labels = []
for idx, (i, row) in enumerate(df.iterrows()):
    ax = axes[idx]
    ax.axis('off')

    if row['Image_File'] and os.path.exists(row['Image_File']):
        mol_img = fig1_structure(row)
        ax.imshow(mol_img)

        compound_name = str(row['Compound'])

        max_chars = 20
        if len(compound_name) > max_chars:
            if '-' in compound_name:
                parts = compound_name.split('-')
                lines = []
                current = parts[0]
                for part in parts[1:]:
                    if len(current + '-' + part) <= max_chars:
                        current += '-' + part
                    else:
                        lines.append(current)
                        current = part
                lines.append(current)
                compound_name = '\n'.join(lines)
            elif len(compound_name) > max_chars * 2:
                compound_name = compound_name[:max_chars*2-3] + '...'

        label = ax.text(0.5, -0.12, compound_name,
                ha='center', va='top',
                fontsize=SUBMITTED_LABEL_SIZE if FIT_LABELS else step["structure_label_size"],
                transform=ax.transAxes,
                family=plt.rcParams['font.family'][0],
                wrap=True)
        labels.append((str(row['Compound']), label))

for idx in range(len(df), len(axes)):
    axes[idx].axis('off')

plt.tight_layout(pad=1.5)
canvas_path = os.path.join(output_folder,"all_compounds_structures.png")
if FIT_LABELS:
    # layout and canvas are those of the submitted figure; only the labels are redrawn larger
    canvas_box = fig.get_tightbbox(fig.canvas.get_renderer()).padded(plt.rcParams['savefig.pad_inches'])
    too_wide = fit_labels(fig, axes, labels)
    if too_wide:
        print(f"Labels wider than their cell: {too_wide}")
    save_figure(canvas_path, dpi=step["dpi_structures"], bbox_inches=canvas_box, facecolor='white')
else:
    save_figure(canvas_path, dpi=step["dpi_structures"], bbox_inches='tight', facecolor='white')
plt.close()
print(f"All compounds structure image saved: {canvas_path}")



def create_a4_grid_image(df, output_path, cols=5):
    """Grid of all structures with names on an A4 canvas (2480 x 3508 px times GRID_SCALE)."""
    mols = df['Drawing'].tolist()
    legends = df['Compound'].tolist()
    num_mols = len(mols)

    a4_width_px = 2480 * GRID_SCALE
    a4_height_px = 3508 * GRID_SCALE

    mols_per_row = cols
    num_rows = math.ceil(num_mols / mols_per_row)
    margin = 80 * GRID_SCALE

    cell_width = (a4_width_px - (mols_per_row + 1) * margin) // mols_per_row
    cell_height = (a4_height_px - (num_rows + 1) * margin) // num_rows

    canvas = PILImage.new("RGB", (a4_width_px, a4_height_px), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    try:
        font = ImageFont.truetype(step["grid_font_file"], round(45 * GRID_SCALE * GRID_NAME_FACTOR))
    except IOError:
        font = ImageFont.load_default()
        note(f"structure grid: font file {step['grid_font_file']} not found, default PIL font used")

    mol_img_width = int(cell_width * 1.2)
    mol_img_height = int(cell_height * 1)
    mol_img_size = (mol_img_width, mol_img_height)

    x_pos_start, y_pos_start = margin, margin

    for idx in range(num_mols):
        mol = mols[idx]
        legend = legends[idx]

        row_idx = idx // mols_per_row
        col_idx = idx % mols_per_row

        x_pos = x_pos_start + col_idx * (cell_width + margin)
        y_pos = y_pos_start + row_idx * (cell_height + margin)

        if mol:
            try:
                img_mol = MolToImage(mol, size=mol_img_size, options=scaled_options(GRID_SCALE))
                canvas.paste(img_mol, (x_pos, y_pos + (mol_img_height - img_mol.height) // 2))
            except Exception as e:
                print(f"Failed to render molecule {legend}: {e}")

            words = str(legend).split()
            lines = []
            current_line = ""
            for word in words:
                if draw.textlength(current_line + " " + word, font=font) < mol_img_width:
                    current_line += " " + word if current_line else word
                else:
                    lines.append(current_line)
                    current_line = word
            lines.append(current_line)
            lines = [line for line in lines if line]

            text_y = y_pos + mol_img_height + 15 * GRID_SCALE
            for line in lines:
                text_width, _ = draw.textlength(line, font=font), font.getmetrics()[0]
                text_x = x_pos + (mol_img_width - text_width) / 2
                draw.text((text_x, text_y), line, fill=(0, 0, 0), font=font)
                text_y += draw.textlength(line, font=font) + 10 * GRID_SCALE

    canvas.save(output_path, dpi=(step["dpi_a4_grid"] * GRID_SCALE, step["dpi_a4_grid"] * GRID_SCALE))
    print(f"A4 grid image saved to: {output_path}")


if __name__ == "__main__":
    try:
        df = pd.read_excel(file_path) if file_path.lower().endswith('.xlsx') else pd.read_csv(file_path, sep='\t')
        if 'SMILES' not in df.columns:
            print("\nError: 'SMILES' column not found. Cannot process compounds.")
            exit()

        df['Molecule'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(str(x)) if pd.notna(x) else None)
        df.dropna(subset=['Molecule'], inplace=True)
        add_drawings(df)

        if df.empty:
            print("\nError: No valid compounds found after processing SMILES strings.")
            exit()

        output_image_path = os.path.join(output_folder, "all_compounds_structures_A4_high_quality.png")
        create_a4_grid_image(df, output_image_path, cols=4)

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}. Please check the path and filename.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

finish_step(cfg, TREATMENT)
