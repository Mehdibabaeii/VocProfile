#!/usr/bin/env python3
#
# VocProfile - compare_outputs.py
#
# Compares the tables of two run folders cell by cell. Numbers are equal when they agree
# to a number of significant digits (default 6); text must be identical. Figures and QC
# files are not compared. Tables under a folder named with --machine-dependent are
# compared and reported in a block of their own; a difference there does not fail the
# comparison.
#
# in   two folders of results; optionally a name map with the columns old, new
#      (relative paths) when the first folder uses other file names
# out  a summary on the terminal; exit status 1 if a table differs or is missing
#
# Usage:
#   python tools/compare_outputs.py <folder A> <folder B> [--map names.tsv] [--digits 6]
#                                   [--machine-dependent <folder> ...]
#
# Notes: docs/NOTES.md, comparing two runs
#
import csv
import math
import os
import re
import sys

from openpyxl import load_workbook

TABLES = (".xlsx", ".csv", ".tsv")
SKIPPED_FOLDERS = ("logs",)


def list_files(folder):
    found = {}
    for path, _, names in os.walk(folder):
        relative = os.path.relpath(path, folder)
        if relative.split(os.sep)[0] in SKIPPED_FOLDERS:
            continue
        for name in names:
            found[os.path.normpath(os.path.join(relative, name)).replace(os.sep, "/")] = os.path.join(path, name)
    return found


def as_number(value):
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def same(first, second, digits):
    one, two = as_number(first), as_number(second)
    if one is not None and two is not None:
        if one == two or (math.isnan(one) and math.isnan(two)):
            return True
        return abs(one - two) <= 0.5 * 10 ** (-digits) * max(abs(one), abs(two))
    empty = (None, "")
    return first == second or (first in empty and second in empty)


def read_sheets(path):
    """{sheet: rows}; a text table is one sheet."""
    if path.lower().endswith(".xlsx"):
        book = load_workbook(path, read_only=True, data_only=True)
        sheets = {sheet.title: [list(row) for row in sheet.iter_rows(values_only=True)] for sheet in book.worksheets}
        book.close()
        return sheets
    delimiter = "\t" if path.lower().endswith(".tsv") else ","
    with open(path, newline="", encoding="utf-8") as handle:
        return {"table": [row for row in csv.reader(handle, delimiter=delimiter)]}


def compare_tables(first, second, digits):
    """List of differences, at most five per sheet."""
    one, two = read_sheets(first), read_sheets(second)
    if list(one) != list(two):
        return [f"sheets differ: {list(one)} and {list(two)}"]
    problems = []
    for sheet in one:
        rows_a, rows_b = one[sheet], two[sheet]
        if len(rows_a) != len(rows_b):
            problems.append(f"{sheet}: {len(rows_a)} and {len(rows_b)} rows")
            continue
        shown = 0
        for number, (row_a, row_b) in enumerate(zip(rows_a, rows_b), start=1):
            width = max(len(row_a), len(row_b))
            row_a = list(row_a) + [None] * (width - len(row_a))
            row_b = list(row_b) + [None] * (width - len(row_b))
            for column, (a, b) in enumerate(zip(row_a, row_b), start=1):
                if not same(a, b, digits) and shown < 5:
                    problems.append(f"{sheet}, row {number}, column {column}: {a!r} and {b!r}")
                    shown += 1
    return problems


def renamer(name_map):
    """Function from an old relative path to the present one. {word} in a row of the map
    stands for one folder or file name part and is carried over to the new path."""
    exact = {old: new for old, new in name_map.items() if "{" not in old}
    patterns = []
    for old, new in name_map.items():
        if "{" in old:
            seen = set()

            def group(match):
                word = match.group(1)
                if word in seen:
                    return f"(?P={word})"
                seen.add(word)
                return f"(?P<{word}>[^/]+?)"
            expression = re.sub(r"\\\{(\w+)\\\}", group, re.escape(old))
            patterns.append((re.compile("^" + expression + "$"), new))

    def rename(name):
        if name in exact:
            return exact[name]
        for pattern, new in patterns:
            match = pattern.match(name)
            if match:
                return new.format(**match.groupdict())
        return name
    return rename


def main(argv):
    if len(argv) < 3:
        sys.exit("usage: python tools/compare_outputs.py <folder A> <folder B> [--map names.tsv] [--digits 6] "
                 "[--machine-dependent <folder> ...]")
    folder_a, folder_b = argv[1], argv[2]
    digits, name_map, loose = 6, {}, []
    rename = lambda name: name
    options = argv[3:]
    while options:
        option = options.pop(0)
        if option == "--digits" and options:
            digits = int(options.pop(0))
        elif option == "--machine-dependent" and options:
            loose.append(options.pop(0).strip("/") + "/")
        elif option == "--map" and options:
            with open(options.pop(0), encoding="utf-8") as handle:
                for row in csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"):
                    name_map[row["old"]] = row["new"]
            rename = renamer(name_map)
        else:
            sys.exit(f"ERROR: unknown argument: {option}")
    for folder in (folder_a, folder_b):
        if not os.path.isdir(folder):
            sys.exit(f"ERROR: folder not found: {folder}")
    files_a = {rename(name): path for name, path in list_files(folder_a).items()}
    files_b = list_files(folder_b)
    equal, different, missing = 0, [], []
    loose_equal, loose_different = 0, []
    for name in sorted(set(files_a) | set(files_b)):
        if not name.lower().endswith(TABLES):
            continue
        if name not in files_a or name not in files_b:
            missing.append((name, "A" if name not in files_a else "B"))
            continue
        problems = compare_tables(files_a[name], files_b[name], digits)
        if name.startswith(tuple(loose)) and loose:
            if problems:
                loose_different.append(name)
            else:
                loose_equal += 1
        elif problems:
            different.append((name, problems))
        else:
            equal += 1
    for name, side in missing:
        print(f"only in {'B' if side == 'A' else 'A'}: {name}")
    for name, problems in different:
        print(f"DIFFERENT  {name}")
        for problem in problems:
            print(f"    {problem}")
    missing_tables = [name for name, _ in missing]
    print(f"tables compared: {equal + len(different)}, equal to {digits} significant digits: {equal}, "
          f"different: {len(different)}, tables present on one side only: {len(missing_tables)}")
    if loose:
        print(f"machine-dependent tables ({', '.join(folder.rstrip('/') for folder in loose)}), not part of the result: "
              f"{loose_equal} equal, {len(loose_different)} different")
        for name in loose_different:
            print(f"    differs on this machine: {name}")
    if different or missing_tables:
        print("RESULT: the two runs differ")
        return 1
    print("RESULT: all tables are equal")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
