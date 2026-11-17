#!/usr/bin/env python3
#
# VocProfile - lib_run.py
#
# What every step does besides its own computation: check the input files, write
# figures and tables with the wording of the label state, and leave a QC file.
#
# in   the settings of lib_config
# out  results/<step>/<id>_qc.txt; figures moved to figures/<step>/
#
# Usage:
#   imported by the step scripts
#
# Notes: docs/NOTES.md, what every step writes
#
import os
import platform
import re
import shutil
import sys

import numpy as np
import pandas as pd

from lib_config import INPUT_KEYS, die, version

FIGURE_TYPES = (".png", ".html")
_state = {"cfg": None, "scope": "", "inputs": [], "notes": [], "later": {}}


def start_step(cfg, *inputs):
    """Check the named input files and remember what was read for the QC file."""
    _state.update(cfg=cfg, scope="", inputs=[], notes=[], data_names=False, later={})
    for kind in inputs:
        _state["inputs"].extend(check_input(cfg, kind))


def set_scope(grouping):
    """Grouping of the input set in work ('clusters' or 'k_groups'); decides some figure wording."""
    _state["scope"] = grouping or ""


def note(text):
    """Add a line to the QC file of the step."""
    _state["notes"].append(str(text))


def notes():
    """The lines noted so far."""
    return list(_state["notes"])


def _duplicates(values):
    series = pd.Series(list(values))
    return sorted(series[series.duplicated()].astype(str).unique())


def check_input(cfg, kind):
    """Checks on one input file. Returns the lines that describe it in the QC file."""
    settings = cfg["settings"]
    path = settings.path(INPUT_KEYS[kind])
    if not os.path.isfile(path):
        die(f"{INPUT_KEYS[kind]}: file not found: {path}")
    lines = [f"{INPUT_KEYS[kind]}: {os.path.basename(path)}"]
    accession, replicate, code = cfg["accession_column"], cfg["replicate_column"], cfg["nondetect_code"]
    if kind == "replicates":
        book = pd.read_excel(path, sheet_name=None, na_values=code, engine="openpyxl")
        reference = None
        for dataset in settings.datasets:
            sheet = settings.sheet(dataset)
            if sheet not in book:
                die(f"sheet '{sheet}' (SHEET_{dataset.upper()}) not found in {path}")
            frame = book[sheet]
            for column in (accession, replicate):
                if column not in frame.columns:
                    die(f"column '{column}' not found in sheet '{sheet}' of {path}")
            values = frame.drop(columns=[accession, replicate])
            text_columns = [str(c) for c in values.columns if not pd.api.types.is_numeric_dtype(values[c])
                            and values[c].notna().any()]
            if text_columns:
                die(f"non-numeric values in sheet '{sheet}' of {path}: columns {text_columns[:5]}")
            names = sorted(frame[accession].astype(str).unique())
            if reference is None:
                reference = names
            elif names != reference:
                differing = sorted(set(names) ^ set(reference))
                die(f"accession names differ between sheets of {path}: {differing[:10]}")
            counts = frame.groupby(accession).size()
            lines.append(f"  {dataset}: {len(frame)} runs, {len(names)} accessions, {values.shape[1]} compounds; "
                         f"replicates per accession {counts.min()} to {counts.max()}; "
                         f"runs without detection {int(values.isna().sum().sum())} of {values.size} values")
    elif kind in ("means", "means_clusters"):
        book = pd.read_excel(path, sheet_name=None, header=None, engine="openpyxl")
        datasets = settings.datasets if kind == "means" else [settings.flat["MAIN_DATASET"]]
        rows = len(cfg["header_rows"])
        for dataset in datasets:
            sheet = settings.sheet(dataset)
            if sheet not in book:
                die(f"sheet '{sheet}' (SHEET_{dataset.upper()}) not found in {path}")
            frame = book[sheet]
            names = frame.iloc[rows:, 0].astype(str)
            if _duplicates(names):
                die(f"duplicated accession names in sheet '{sheet}' of {path}: {_duplicates(names)[:10]}")
            values = frame.iloc[rows:, 2:].apply(pd.to_numeric, errors="coerce")
            lines.append(f"  {dataset}: {len(names)} accessions, {values.shape[1]} compounds, {rows} header rows; "
                         f"groups {sorted(frame.iloc[rows:, 1].astype(str).unique())}; "
                         f"accession means without detection {int(values.isna().sum().sum())} of {values.size}")
    elif kind == "blup_clusters":
        frame = pd.read_csv(path, header=cfg["header_rows"], index_col=[0, 1])
        names = frame.index.get_level_values(0).astype(str)
        if _duplicates(names):
            die(f"duplicated accession names in {path}: {_duplicates(names)[:10]}")
        lines.append(f"  {len(frame)} accessions, {frame.shape[1]} compounds; "
                     f"clusters {frame.index.get_level_values(1).value_counts().to_dict()}")
    elif kind == "compounds":
        frame = pd.read_excel(path) if path.lower().endswith(".xlsx") else pd.read_csv(path, sep="\t")
        for column in ("Compound", "SMILES", "Main_Class", "Subclass"):
            if column not in frame.columns:
                die(f"column '{column}' not found in {path}")
        without = frame.loc[frame["SMILES"].isna(), "Compound"].astype(str).tolist()
        lines.append(f"  {len(frame)} compounds; without a structure (no SMILES, not drawn): {without or 'none'}")
    return lines


def check_same_accessions(cfg, first, second):
    """Stop if two input files do not hold the same accession names."""
    def names(kind):
        path = cfg["settings"].path(INPUT_KEYS[kind])
        if kind == "replicates":
            sheet = cfg["settings"].sheet(cfg["settings"].flat["MAIN_DATASET"])
            return set(pd.read_excel(path, sheet_name=sheet, engine="openpyxl")[cfg["accession_column"]].astype(str))
        if kind == "blup_clusters":
            return set(pd.read_csv(path, header=cfg["header_rows"], index_col=[0, 1]).index.get_level_values(0).astype(str))
        sheet = cfg["settings"].sheet(cfg["settings"].flat["MAIN_DATASET"])
        frame = pd.read_excel(path, sheet_name=sheet, header=None, engine="openpyxl")
        return set(frame.iloc[len(cfg["header_rows"]):, 0].astype(str))
    one, two = names(first), names(second)
    if one != two:
        die(f"accession names differ between {INPUT_KEYS[first]} and {INPUT_KEYS[second]}: {sorted(one ^ two)[:10]}")
    note(f"accession names equal in {INPUT_KEYS[first]} and {INPUT_KEYS[second]} ({len(one)})")


def _figure_terms():
    labels = _state["cfg"]["labels"]
    if not labels["revised_terms"]:
        return []
    return [(old, new) for old, new, scope in labels["figure_terms"]
            if scope == "figures" or (scope == "cluster_runs" and _state["scope"] == "clusters")]


def _reword(text, pairs):
    for old, new in pairs:
        text = re.sub(r"(?<![A-Za-z])" + re.escape(old) + r"(?![A-Za-z])", new, text)
    return text


def _display_pairs():
    """(name in the data, name shown) for the label state; longest first."""
    names = _state["cfg"]["labels"]["compound_display_names"]
    return sorted(names.items(), key=lambda pair: -len(pair[0]))


def _replace_names(text, pairs):
    for old, new in pairs:
        text = text.replace(old, new)
    return text


def shown(value):
    """A compound name, or a text that holds one, as shown to a reader in the label state."""
    return _replace_names(value, _display_pairs()) if isinstance(value, str) else value


def later(text, replacement):
    """Text for a figure: the submitted wording now, its replacement when the figure is saved.
    The layout is thus computed with the same texts in both label states."""
    if replacement != text and text:
        _state["later"][text] = replacement
    return text


def deferred(function):
    """A wording function whose result is applied when the figure is saved (see later)."""
    return lambda text: later(text, function(text))


def minimum_dpi(dpi):
    """Resolution of a panel of the paper: the step value, raised to MINIMUM_PANEL_DPI of the label state."""
    lowest = _state["cfg"]["settings"].label("MINIMUM_PANEL_DPI")
    return max(dpi, lowest) if lowest else dpi


def _final_text(text, pairs, names):
    """A figure text in the label state: wording set aside with later(), word pairs, names shown."""
    for old, new in sorted(_state["later"].items(), key=lambda pair: -len(pair[0])):
        text = re.sub(r"(?<![A-Za-z0-9])" + re.escape(old) + r"(?![A-Za-z0-9])", lambda match: new, text)
    return _replace_names(_reword(text, pairs), names)


def _reword_figure(figure, pairs, names):
    """Put the texts of a figure in the terms of the label state without touching the layout."""
    from matplotlib.text import Text
    from matplotlib.ticker import FixedFormatter, FixedLocator

    def new_text(old):
        return _final_text(old, pairs, names)
    for axes in figure.axes:
        for axis in (axes.xaxis, axes.yaxis):
            for minor in (False, True):
                labels = [label.get_text() for label in axis.get_ticklabels(minor=minor)]
                changed = [new_text(label) for label in labels]
                if changed != labels:
                    locations = axis.get_ticklocs(minor=minor)
                    setter = (axis.set_minor_locator, axis.set_minor_formatter) if minor else \
                        (axis.set_major_locator, axis.set_major_formatter)
                    setter[0](FixedLocator(locations))
                    setter[1](FixedFormatter(changed))
    for item in figure.findobj(Text):
        old = item.get_text()
        new = new_text(old)
        if new != old:
            item.set_text(new)


def save_figure(path, panel=False, after_layout=None, **options):
    """plt.savefig, after the text of the figure is put in the terms of the label state.

    With LABELS = revised the canvas is the one of the figure before its texts are replaced,
    so both label states give the same size and the same position of the axes. panel=True
    marks a panel of the paper, saved at MINIMUM_PANEL_DPI or more. after_layout(figure, box) is
    called in that state once the canvas (box, in inches, or None) is fixed."""
    import matplotlib.pyplot as plt
    from matplotlib import cbook
    figure = plt.gcf()
    pairs, names = _figure_terms(), _display_pairs()
    if revised():
        box = getattr(figure, "canvas_box", None)
        dpi = options.get("dpi", figure.dpi)
        if box is not None:
            options.pop("pad_inches", None)
            options["bbox_inches"] = box
        elif options.get("bbox_inches") == "tight":
            pad = options.pop("pad_inches", plt.rcParams["savefig.pad_inches"])
            with cbook._setattr_cm(figure, dpi=dpi):
                figure.canvas.draw()
                box = figure.get_tightbbox(figure.canvas.get_renderer(),
                                           bbox_extra_artists=options.get("bbox_extra_artists"))
            box = box.padded(pad)
            figure.canvas_box = box
            options["bbox_inches"] = box
        figure.canvas.draw()
        _reword_figure(figure, pairs, names)
        if after_layout is not None:
            after_layout(figure, box)
    if panel and "dpi" in options:
        options["dpi"] = minimum_dpi(options["dpi"])
    plt.savefig(path, **options)


def revised():
    return _state["cfg"]["settings"].state == "revised"


def save_legend(path, handles, labels, title, dpi, **options):
    """A legend alone, as its own image, for figures assembled from several panels."""
    import matplotlib.pyplot as plt
    figure = plt.figure(figsize=(2, 2))
    if revised():
        labels = [_final_text(label, _figure_terms(), _display_pairs()) for label in labels]
    legend = figure.legend(handles=handles, labels=labels, title=title,
                           loc="center", frameon=False, **options)
    figure.canvas.draw()
    box = legend.get_window_extent().transformed(figure.dpi_scale_trans.inverted())
    figure.savefig(path, dpi=minimum_dpi(dpi), bbox_inches=box.padded(0.05), facecolor="white")
    plt.close(figure)


def keep_data_names():
    """Tables of this step are read back by later steps: compound names stay those of the data."""
    _state["data_names"] = True


def write_table(frame, target, data_names=False, **options):
    """DataFrame.to_excel, with headings and compound names in the terms of the label state."""
    labels = _state["cfg"]["labels"]
    if labels["revised_terms"]:
        pairs = sorted(labels["terms"].items(), key=lambda pair: -len(pair[0]))
        frame = frame.rename(columns=lambda name: _reword(name, pairs) if isinstance(name, str) else name)
    names = _display_pairs()
    if names and not data_names and not _state.get("data_names"):
        frame = frame.rename(columns=shown, index=shown)
        text_columns = [position for position, kind in enumerate(frame.dtypes) if kind == object]
        if text_columns:
            frame = frame.copy()
            for position in text_columns:
                frame.isetitem(position, frame.iloc[:, position].map(shown))
    frame.to_excel(target, **options)


def read_own_table(path, **options):
    """Read a table this step wrote with write_table, with the compound names of the data."""
    frame = pd.read_excel(path, **options)
    back = {new: old for old, new in _display_pairs()}
    return frame.rename(columns=lambda name: back.get(name, name), index=lambda name: back.get(name, name))


def show_names_in_text_files(folder, keep=()):
    """Replace the compound names of the data by the names shown, in the CSV files of a folder."""
    names = _display_pairs()
    if not names:
        return
    for name in sorted(os.listdir(folder)):
        if name.lower().endswith(".csv") and name not in keep:
            path = os.path.join(folder, name)
            with open(path, encoding="utf-8", newline="") as handle:
                text = handle.read()
            changed = _replace_names(text, names)
            if changed != text:
                with open(path, "w", encoding="utf-8", newline="") as handle:
                    handle.write(changed)


def use_font(family):
    """Set the font family of the figures; fall back to DejaVu Sans when it is not installed."""
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    installed = {font.name for font in font_manager.fontManager.ttflist}
    used = family if family in installed else "DejaVu Sans"
    plt.rcParams["font.family"] = used
    note(f"font family asked: {family}; used: {used}")
    return used


def finish_step(cfg, treatment=()):
    """Move the figures of the step to figures/<step>/ and write the QC file."""
    settings, step_id = cfg["settings"], cfg["step_id"]
    script = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    results = os.path.join(settings.work_dir, "results", script)
    figures = os.path.join(settings.work_dir, "figures", script)
    os.makedirs(results, exist_ok=True)
    moved = 0
    for folder, _, files in os.walk(results):
        for name in files:
            if name.lower().endswith(FIGURE_TYPES):
                target = os.path.join(figures, os.path.relpath(folder, results))
                os.makedirs(target, exist_ok=True)
                shutil.move(os.path.join(folder, name), os.path.join(target, name))
                moved += 1
    for folder, _, _ in os.walk(results, topdown=False):
        if folder != results and not os.listdir(folder):
            os.rmdir(folder)
    import matplotlib
    import scipy
    import sklearn
    lines = [f"pipeline version: {version()}", f"step {step_id}: {script}.py", f"config: {settings.config_path}",
             f"LABELS = {settings.state}", "", "inputs"]
    lines += [f"  {line}" for line in _state["inputs"]] or ["  none"]
    lines += ["", "treatment of the data"] + [f"  {line}" for line in treatment]
    lines += ["", "notes"] + ([f"  {line}" for line in _state["notes"]] or ["  none"])
    lines += ["", "seeds", f"  SEED = {settings.flat['SEED']} (k-means, t-SNE, random forest, odor index steps)",
              f"  JITTER_SEED = {settings.flat['JITTER_SEED'] or 'empty (unseeded)'} (jitter of points in box plots)"]
    lines += ["", "settings"] + [f"  {key} = {value}" for key, value in settings.flat.items()]
    lines += ["", "versions", f"  Python {platform.python_version()}", f"  numpy {np.__version__}", f"  pandas {pd.__version__}",
              f"  scipy {scipy.__version__}", f"  scikit-learn {sklearn.__version__}", f"  matplotlib {matplotlib.__version__}"]
    for module in ("seaborn", "statsmodels", "openpyxl", "rdkit", "PIL", "adjustText"):
        if module in sys.modules:
            lines.append(f"  {module} {getattr(sys.modules[module], '__version__', 'unknown')}")
    lines += ["", f"figures moved to figures/{script}: {moved}"]
    with open(os.path.join(results, f"{step_id}_qc.txt"), "w", encoding="utf-8", newline="\n") as handle:
        handle.write(_replace_names("\n".join(lines), _display_pairs()) + "\n")
