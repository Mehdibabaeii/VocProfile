#!/usr/bin/env python3
#
# VocProfile - lib_config.py
#
# Reads config.txt and the tables it names, checks them, and hands every step its
# settings. Nothing is computed here.
#
# in   config.txt (KEY = VALUE), the tables under config/
# out  none
#
# Usage:
#   python scripts/lib_config.py <config.txt> [--set KEY=VALUE ...]     check and print the run plan
#
# Notes: docs/NOTES.md, configuration
#
import csv
import os
import re
import sys

CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE = os.path.join(CODE_DIR, "config.txt")
STATES = ("submitted", "revised")
INPUT_KEYS = {"compounds": "COMPOUND_FILE", "replicates": "REPLICATE_FILE", "means": "MEANS_FILE",
              "means_clusters": "MEANS_CLUSTER_FILE", "blup_clusters": "BLUP_CLUSTER_FILE"}


def die(message):
    sys.exit(f"ERROR: {message}")


def version():
    with open(os.path.join(CODE_DIR, "VERSION"), encoding="utf-8") as handle:
        return handle.read().strip()


def read_keys(path):
    """KEY = VALUE lines of a config file, in order. '#' starts a comment."""
    if not os.path.isfile(path):
        die(f"configuration file not found: {path}")
    values = {}
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            line = line.split("#", 1)[0].strip() if not line.lstrip().startswith("#") else ""
            if not line:
                continue
            if "=" not in line:
                die(f"{path}, line {number}: expected KEY = VALUE")
            key, value = (part.strip() for part in line.split("=", 1))
            if key in values:
                die(f"{path}: key {key} is set twice")
            values[key] = value
    return values


def typed(value):
    """Text of a setting as a number, a truth value, None (empty) or text."""
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    if text.lower() in ("true", "yes", "on"):
        return True
    if text.lower() in ("false", "no", "off"):
        return False
    for convert in (int, float):
        try:
            return convert(text)
        except ValueError:
            pass
    return text


def as_list(value):
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def read_table(path):
    """Rows of a tab-separated table as dictionaries; lines starting with '#' are skipped."""
    with open(path, encoding="utf-8", newline="") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


class Settings:
    """The settings of one run: the keys of config.txt, the tables and the label state."""

    def __init__(self, config_path, overrides=None):
        self.config_path = os.path.abspath(config_path)
        template = read_keys(TEMPLATE)
        flat = read_keys(self.config_path)
        self.overrides = dict(overrides or {})
        for key in self.overrides:
            if key not in template:
                die(f"--set {key}: config.txt declares no such setting")
        flat.update(self.overrides)
        unknown = [key for key in flat if key not in template]
        if unknown:
            die(f"unknown key(s) in {self.config_path}: {', '.join(unknown)}")
        missing = [key for key in template if key not in flat]
        if missing:
            die(f"key(s) missing from {self.config_path}: {', '.join(missing)}")
        self.flat = flat
        for key in ("WORK_DIR", "INPUT_DIR"):
            if not flat[key]:
                die(f"{key} is empty")
        if flat["LABELS"] not in STATES:
            die(f"LABELS = {flat['LABELS']}: use 'revised' or 'submitted'")
        self.state = flat["LABELS"]
        self.work_dir = os.path.abspath(self.expand(flat["WORK_DIR"]))
        self.tables = {}
        for key in flat:
            if key.endswith("_TABLE"):
                path = self.path(key)
                if not os.path.isfile(path):
                    die(f"{key}: file not found: {path}")
                self.tables[key] = read_table(path)
        self.datasets = as_list(flat["DATASETS"])
        for name in self.datasets:
            for prefix in ("SHEET_", "LABEL_", "SHORT_"):
                if prefix + name.upper() not in flat:
                    die(f"dataset '{name}' of DATASETS has no {prefix}{name.upper()}")
        for key in ("MAIN_DATASET", "ABUNDANCE_DATASET"):
            if flat[key] not in self.datasets:
                die(f"{key} = {flat[key]} is not one of DATASETS")

    def expand(self, value):
        return value.replace("CODE_DIR", CODE_DIR, 1) if value.startswith("CODE_DIR") else value

    def path(self, key):
        """Absolute path of a file setting. Input files are taken from INPUT_DIR unless absolute."""
        value = self.expand(self.flat[key])
        value = value.replace("WORK_DIR", self.work_dir, 1) if value.startswith("WORK_DIR") else value
        if os.path.isabs(value):
            return value
        if key.endswith("_TABLE"):
            return os.path.join(os.path.dirname(self.config_path), value)
        return os.path.join(os.path.abspath(self.expand(self.flat["INPUT_DIR"])), value)

    def check_input_files(self):
        """Input files must exist at the start, except those a step of this run writes under WORK_DIR."""
        for key in INPUT_KEYS.values():
            path = self.path(key)
            if not path.startswith(self.work_dir + os.sep) and not os.path.isfile(path):
                die(f"{key}: file not found: {path}")

    def value(self, key):
        return typed(self.flat[key])

    def label(self, setting):
        """Value of a row of LABELS_TABLE for the label state of the run."""
        for row in self.tables["LABELS_TABLE"]:
            if row["setting"] == setting:
                return typed(row[self.state])
        die(f"LABELS_TABLE has no row '{setting}'")

    def tables_label(self, setting, state):
        """Value of a row of LABELS_TABLE for a named label state."""
        for row in self.tables["LABELS_TABLE"]:
            if row["setting"] == setting:
                return typed(row[state])
        die(f"LABELS_TABLE has no row '{setting}'")

    def sheet(self, dataset):
        return self.flat["SHEET_" + dataset.upper()]

    def dataset_label(self, dataset):
        return self.flat["LABEL_" + dataset.upper()] or dataset


def parse_arguments(argv):
    """<config.txt> [--set KEY=VALUE ...] [--input-set NAME]"""
    if len(argv) < 2 or argv[1].startswith("-"):
        die(f"usage: python {os.path.basename(argv[0])} <config.txt> [--set KEY=VALUE ...] [--input-set NAME]")
    overrides, input_set, position = {}, None, 2
    while position < len(argv):
        option = argv[position]
        if option == "--set" and position + 1 < len(argv) and "=" in argv[position + 1]:
            key, value = argv[position + 1].split("=", 1)
            overrides[key.strip()] = value.strip()
        elif option == "--input-set" and position + 1 < len(argv):
            input_set = argv[position + 1]
        else:
            die(f"unknown argument: {option}")
        position += 2
    return argv[1], overrides, input_set


class Arguments:
    def __init__(self, input_set):
        self.input_set = input_set


def figure_style(settings, name):
    rows = [row for row in settings.tables["FIGURE_STYLES_TABLE"] if row["style"] == name]
    if not rows:
        die(f"FIGURE_STYLES_TABLE has no style '{name}'")
    return {row["setting"]: typed(row[settings.state]) for row in rows}


def colours(settings, group):
    return {row["key"]: row["colour"] for row in settings.tables["COLOUR_TABLE"] if row["group"] == group}


def clusters(settings):
    """Rows of CLUSTER_TABLE that carry a name, in table order (largest cluster first)."""
    return [row for row in settings.tables["CLUSTER_TABLE"] if row["name"]]


def input_sets(settings, step_id):
    sets = {}
    cluster_colours = {row["name"]: row["colour"] for row in clusters(settings)}
    for row in settings.tables["INPUT_SETS_TABLE"]:
        if row["step"] != step_id:
            continue
        datasets = as_list(row["datasets"])
        entry = {
            "input": row["input"],
            "sheets": [settings.dataset_label(name) for name in datasets],
            "output_label": row["output_label"],
            "output_folder": row["output_folder"],
            "compounds": row["compounds"],
            "cluster_legend": typed(row["cluster_legend"]),
            "layout": row["layout"],
            "transform": row["transform"],
            "grouping": row["grouping"],
            "group_axis_label": row["group_axis_label_" + settings.state],
            "compound_legend_title": row["compound_legend_title_" + settings.state],
            "group_colours": cluster_colours if row["grouping"] == "clusters" else colours(settings, "k_group"),
        }
        for style in ("heatmap", "correlation"):
            if row[style + "_style"]:
                entry[style] = figure_style(settings, row[style + "_style"])
        sets[row["set"]] = entry
    return sets


def label_block(settings):
    wording = settings.tables["LABEL_WORDING_TABLE"]
    names = {int(row["id"]): row["name"] for row in clusters(settings)}
    return {
        "revised_terms": settings.state == "revised",
        "revised_terms_in_titles": settings.state == "revised",
        "terms": {row["old"]: row["new"] for row in wording if row["scope"] == "terms"},
        "figure_terms": [(row["old"], row["new"], row["scope"]) for row in wording
                         if row["scope"] in ("figures", "cluster_runs")],
        "all_display_names": {row["data_name"]: row["display_name"] for row in settings.tables["DISPLAY_NAMES_TABLE"]},
        "compound_display_names": ({row["data_name"]: row["display_name"]
                                    for row in settings.tables["DISPLAY_NAMES_TABLE"]}
                                   if settings.state == "revised" else {}),
        "fig2c_y_axis": settings.label("DETECTION_Y_AXIS"),
        "cluster_names": names if settings.label("CLUSTER_NAMES_ON_FIGURES") else {},
        "dendrogram_legend_title": settings.label("DENDROGRAM_LEGEND_TITLE"),
        "scatter_legend_title": settings.label("SCATTER_LEGEND_TITLE"),
        "revised_odor_wording": settings.state == "revised",
        "odor_index_wording": [[row["old"], row["new"]] for row in wording if row["scope"] == "odor_index"],
        "fig_s5_significance": settings.label("ODOR_RANKING_SIGNIFICANCE"),
    }


def step_block(settings, step_id):
    """Settings of one step, under the names the step scripts use."""
    value, label = settings.value, settings.label
    all_labels = [settings.dataset_label(name) for name in settings.datasets]
    main_label = settings.dataset_label(settings.flat["MAIN_DATASET"])
    blocks = {
        "01": lambda: {"input": "compounds", "output_folder": "compound_report",
                       "dpi_structures": label("STRUCTURE_FIGURE_DPI"), "structure_label_size": label("STRUCTURE_LABEL_SIZE"),
                       "fig1_fit_labels": label("STRUCTURE_LABELS_FITTED"), "structure_label_min_size": value("STRUCTURE_LABEL_MIN_SIZE"),
                       "structure_label_offset": value("STRUCTURE_LABEL_OFFSET"), "fig1_structure_px": label("STRUCTURE_RENDER_PX"),
                       "fig1_atom_label_min_px": label("STRUCTURE_ATOM_LABEL_MIN_PX"), "dpi_a4_grid": value("STRUCTURE_GRID_DPI"),
                       "table_structure_px": value("TABLE_STRUCTURE_PX"), "font_family": value("FONT_FAMILY"),
                       "grid_font_file": value("STRUCTURE_GRID_FONT_FILE"),
                       "grid_scale": label("STRUCTURE_GRID_SCALE"), "grid_name_factor": label("STRUCTURE_GRID_NAME_FACTOR")},
        "02": lambda: {"input": "means", "sheets": all_labels, "output_folder": "statistics", "dpi": value("DPI_STATISTICS")},
        "03": lambda: {"input": "replicates", "sheets": all_labels, "output_file": "dataset_comparison.xlsx",
                       "short_names": {settings.dataset_label(name): settings.flat["SHORT_" + name.upper()]
                                       for name in settings.datasets}},
        "04": lambda: {"input": "replicates", "sheets": all_labels, "output_file": "replicate_statistics.xlsx"},
        "05": lambda: {"input": "means", "sheets": all_labels, "output_file": "accession_means_filled.xlsx"},
        "06": lambda: {"input": "replicates", "sheets": all_labels, "output_folder": "blup",
                       "cluster_remap": cluster_remap(settings), "dpi": value("DPI_CLUSTER_FIGURE"),
                       "n_init": value("KMEANS_N_INIT"), "clip_limit": value("CLIP_LIMIT"),
                       "variance_floor": value("ACCESSION_VARIANCE_FLOOR"), "nondetect_fraction": value("NONDETECT_FRACTION")},
        "07": lambda: {"input": "replicates", "sheets": all_labels, "output_folder": "cluster_number",
                       "k_min": value("K_MIN"), "k_max": value("K_MAX"), "dpi": value("DPI_CLUSTER_NUMBER"),
                       "n_init": value("KMEANS_N_INIT"), "variance_floor": value("ACCESSION_VARIANCE_FLOOR"),
                       "nondetect_fraction": value("NONDETECT_FRACTION")},
        "08": lambda: {"input": "means", "blup_workbook": f"blup_{main_label}/blup_results_{main_label}.xlsx",
                       "blup_sheet": "2_BLUP_Estimates", "cluster_names": None,
                       "expected_sizes": expected_sizes(settings)},
        "09": lambda: {"dpi": value("DPI_CLASS_SUMMARIES")},
        "10": lambda: {"input": "means", "sheets": all_labels, "output_folder": "detection", "dpi": value("DPI_DETECTION")},
        "11": lambda: {"input": "means", "sheets": all_labels, "output_folder": "distributions", "dpi": value("DPI_DISTRIBUTIONS")},
        "12": lambda: {"input": "means", "sheets": all_labels, "output_folder": "distributions_filled",
                       "dpi": value("DPI_DISTRIBUTIONS")},
        "13": lambda: {"dpi": value("DPI_GROUP_PROFILES")},
        "14": lambda: {"dpi": value("DPI_GROUP_DISTRIBUTIONS"), "heatmap_file": "compound_heatmap.png",
                       "clip_limit": value("CLIP_LIMIT"), "tukey_alpha": value("TUKEY_ALPHA"),
                       "unit_labels": {transform: {"abundance": label(f"UNIT_ABUNDANCE_{transform.upper()}"),
                                                   "concentration": label(f"UNIT_CONCENTRATION_{transform.upper()}")}
                                       for transform in ("log1p", "filled", "none")}},
        "15": lambda: {"input": "means_clusters", "sheet": main_label, "output_folder": "correlation_exports",
                       "use_display_names": True, "clip_limit": value("CLIP_LIMIT")},
        "16": lambda: {"input": "replicates", "sheet": main_label, "labels": "blup_clusters",
                       "dpi": value("DPI_ODOR_INDEX"), "dpi_high": value("DPI_ODOR_INDEX_HIGH"),
                       "font_family": value("FONT_FAMILY"), "clip_limit": value("CLIP_LIMIT")},
        "17": lambda: {"input": "replicates", "sheet": main_label, "labels": "blup_clusters",
                       "output_folder": "radar", "show_legend": value("RADAR_LEGEND"), "dpi": value("DPI_RADAR"),
                       "font_family": value("FONT_FAMILY"), "clip_limit": value("CLIP_LIMIT")},
    }
    block = blocks[step_id]() if step_id in blocks else {}
    sets = input_sets(settings, step_id)
    if sets:
        block["input_sets"] = sets
    return block


def cluster_remap(settings):
    pairs = as_list(settings.flat["CLUSTER_REMAP"])
    return {int(pair.split(":")[0]): int(pair.split(":")[1]) for pair in pairs} or None


def expected_sizes(settings):
    sizes = as_list(settings.flat["CLUSTER_SIZES_EXPECTED"])
    names = [row["name"] for row in clusters(settings)]
    return {name: int(size) for name, size in zip(names, sizes)} or None


def load_config(step_id):
    """Settings for a step script: (cfg, step, args). The script is called as <script> <config.txt> [...]."""
    config_path, overrides, chosen_set = parse_arguments(sys.argv)
    settings = Settings(config_path, overrides)
    table = settings.tables
    cfg = {
        "settings": settings,
        "step_id": step_id,
        "seed": settings.value("SEED"),
        "random": {"jitter_seed": settings.value("JITTER_SEED")},
        "n_clusters": settings.value("N_CLUSTERS"),
        "clusters": [row["name"] for row in clusters(settings)],
        "cluster_rows": table["CLUSTER_TABLE"],
        "markers": [row["compound"] for row in table["MARKER_TABLE"]],
        "abundance_sheet": settings.dataset_label(settings.flat["ABUNDANCE_DATASET"]),
        "main_sheet": settings.dataset_label(settings.flat["MAIN_DATASET"]),
        "sheet_of": {settings.dataset_label(name): settings.sheet(name) for name in settings.datasets},
        "accession_column": settings.flat["ACCESSION_COLUMN"],
        "replicate_column": settings.flat["REPLICATE_COLUMN"],
        "nondetect_code": settings.flat["NONDETECT_CODE"],
        "header_rows": list(range(settings.value("MEANS_HEADER_ROWS"))),
        "main_class_colours": colours(settings, "main_class"),
        "subclass_order": {name: [row["subclass"] for row in table["SUBCLASS_ORDER_TABLE"] if row["order"] == name]
                           for name in ("distributions", "profiles")},
        "significance_levels": [float(level) for level in as_list(settings.flat["SIGNIFICANCE_LEVELS"])],
        "cluster_colours": {row["name"]: row["colour"] for row in clusters(settings)},
        "cluster_colours_by_id": {int(row["id"]): row["colour"] for row in table["CLUSTER_TABLE"]},
        "cluster_descriptors": {row["name"]: row["descriptor"] for row in clusters(settings)},
        "odor_thresholds": {row["compound"]: {"ODT": float(row["threshold"]), "odor": row["odor"], "category": row["category"]}
                            for row in table["ODOR_THRESHOLD_TABLE"]},
        "labels": label_block(settings),
        "revision": {key: settings.value(key) for key in ("N_RESAMPLE", "N_SEEDS", "N_MONTE_CARLO")},
    }
    step = step_block(settings, step_id)
    return cfg, step, Arguments(chosen_set)


def input_path(cfg, key):
    """Absolute path of an input file; stops if it does not exist."""
    settings = cfg["settings"]
    if key not in INPUT_KEYS:
        die(f"unknown input '{key}'")
    path = settings.path(INPUT_KEYS[key])
    if not os.path.isfile(path):
        die(f"{INPUT_KEYS[key]}: file not found: {path}")
    return path


def table_path(cfg, key):
    return cfg["settings"].path(key)


def require_sheets(path, sheets):
    """Stop if the workbook lacks one of the sheets."""
    from openpyxl import load_workbook
    workbook = load_workbook(path, read_only=True)
    present = workbook.sheetnames
    workbook.close()
    missing = [sheet for sheet in sheets if sheet not in present]
    if missing:
        die(f"sheet(s) {missing} not found in {path} (available: {present})")


def output_dir(cfg, *parts):
    """Create and return a folder under WORK_DIR/results."""
    path = os.path.join(cfg["settings"].work_dir, "results", *parts)
    os.makedirs(path, exist_ok=True)
    return path


def selected_sets(step, args):
    """Input sets of a step, restricted to --input-set if given."""
    sets = step["input_sets"]
    if args.input_set is None:
        return sets
    if args.input_set not in sets:
        die(f"unknown input set '{args.input_set}' (available: {list(sets)})")
    return {args.input_set: sets[args.input_set]}


def term_function(cfg, switch="revised_terms"):
    """Function that rewrites the old terms in a label or table heading when LABELS = revised."""
    labels = cfg["labels"]
    if not labels.get(switch):
        return lambda text: text
    pairs = sorted(labels["terms"].items(), key=lambda pair: -len(pair[0]))

    def rewrite(text):
        if not isinstance(text, str):
            return text
        for old, new in pairs:
            text = re.sub(r"(?<![A-Za-z])" + re.escape(old) + r"(?![A-Za-z])", new, text)
        return text
    return rewrite


def display_names(cfg):
    """Compound names as shown on figures and in exported tables: {name in the data: name shown}."""
    return dict(cfg["labels"]["compound_display_names"])


def print_plan(settings):
    """Lines printed by run.sh after the banner."""
    flat = settings.flat
    print(f"  config     {settings.config_path}")
    print(f"  reads      {os.path.abspath(settings.expand(flat['INPUT_DIR']))}")
    print(f"  reads      {os.path.join(CODE_DIR, 'config')}")
    for name in ("results", "figures", "logs"):
        print(f"  writes     {os.path.join(settings.work_dir, name)}")
    for key, value in settings.overrides.items():
        print(f"  --set      {key} = {value}")
    print(f"  LABELS     {settings.state}")


if __name__ == "__main__":
    wanted = None
    if "--print" in sys.argv:
        position = sys.argv.index("--print")
        wanted = sys.argv[position + 1]
        del sys.argv[position:position + 2]
    path, given, _ = parse_arguments(sys.argv)
    run_settings = Settings(path, given)
    if wanted == "WORK_DIR":
        print(run_settings.work_dir)
    elif wanted:
        print(run_settings.flat[wanted])
    else:
        run_settings.check_input_files()
        print_plan(run_settings)
