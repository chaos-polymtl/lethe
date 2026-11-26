#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Generate the PDF report of a Lethe validation run.

validate_lethe.sh, and the add_to_report function of validation_functions.sh
called by the validate.sh script of every case, list the content of the report
in report_contents.txt, in the run folder. Every line of that file is one of:

  computer <name of the computer on which the run is made>
  case <path of the example> <number of MPI processes> <solver>
  file <case folder> <file name>

where the file lines give, in order, the files of a case folder of the run
that are shown for the case listed before them. The date of the run is the
modification time of current_git_hash.txt, which validate_lethe.sh writes when
the run starts.

The report contains a cover page, a table of contents grouped by solver and,
for every case, a title page followed by its files. The title page lists the
files of the case, and those that are missing or cannot be read. The pages of
the PDF files are scaled onto pages of the report without being rasterized, so
that the plots remain vector graphics, under a header that names the case.
Text files are shown verbatim. The contents link to the title pages, which
link to the files, and the bookmarks of the PDF follow the same structure.

Usage:
  python3 generate_validation_report.py --run <run folder> [--output <pdf>]

Requires matplotlib and pypdf (python3 -m pip install --user pypdf).
"""

import argparse
import io
import itertools
import os
import sys
import textwrap
import time
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

try:
    from pypdf import PdfReader, PdfWriter, Transformation
    from pypdf.annotations import Link
except ImportError:
    sys.exit("Error: the Python package pypdf is required to generate the "
             "validation report. Install it with: "
             "python3 -m pip install --user pypdf")

#############################################################################
# Constants
#############################################################################

# Files of the run folder
CONTENTS_NAME = "report_contents.txt"
HASH_FILE = "current_git_hash.txt"
REPORT_NAME = "report.pdf"

# Solver shown for a case whose solver is not given in the contents
UNKNOWN_SOLVER = "unknown solver"

# Every page is a US Letter page in landscape. Distances are in points
# (1/72 inch) and, unless stated otherwise, measured from the top left corner.
POINTS_PER_INCH = 72.0
PAGE_WIDTH = 792.0
PAGE_HEIGHT = 612.0

# Colours
TEXT_COLOR = "#1b1f24"
MUTED_COLOR = "#5f6b76"
ACCENT_COLOR = "#1f5f99"
RULE_COLOR = "#c3cbd3"
ALERT_COLOR = "#b3261e"
PAPER_COLOR = "white"

SANS_FONT = "DejaVu Sans"
MONOSPACE_FONT = "DejaVu Sans Mono"
# Width of a character of the monospace font, as a fraction of its size
MONOSPACE_CHARACTER_WIDTH = 0.602

# Margins of the pages
SIDE_MARGIN = 54.0
BOTTOM_MARGIN = 54.0
HEADER_BASELINE = 30.0
# Distance between the bottom of a page and the baseline of its footer
FOOTER_BASELINE = 30.0

# Cover and title pages
TITLE_LEFT = 81.0
TITLE_RIGHT = PAGE_WIDTH - TITLE_LEFT
TITLE_TOP = 150.0
DETAILS_LABEL_WIDTH = 118.0
DETAILS_LEADING = 19.0
# Baseline of the first file listed on a title page, below two rows of details
FILES_FIRST_BASELINE = TITLE_TOP + 108.0 + 2.0 * DETAILS_LEADING + 8.0
FILE_LEADING = 15.0
FILE_NUMBER_WIDTH = 24.0

# Table of contents
CONTENTS_COLUMNS = 2
CONTENTS_TITLE_BASELINE = 84.0
CONTENTS_FIRST_ROW = 116.0
CONTENTS_COLUMN_GAP = 40.0
CONTENTS_HEADING_GAP = 10.0
CONTENTS_HEADING_HEIGHT = 20.0
CONTENTS_CASE_HEIGHT = 15.5
CONTENTS_CASE_INDENT = 14.0
# Marker of the cases with missing output
PROBLEM_MARKER = "\u25cf"

# Text files
TEXT_SIZE = 9.0
TEXT_LEADING = 12.5
TEXT_FIRST_LINE = 70.0
TEXT_COLUMNS = int((PAGE_WIDTH - 2.0 * SIDE_MARGIN) /
                   (MONOSPACE_CHARACTER_WIDTH * TEXT_SIZE))
TEXT_LINES_PER_PAGE = int((PAGE_HEIGHT - TEXT_FIRST_LINE - BOTTOM_MARGIN) /
                          TEXT_LEADING) + 1

# Box in which the pages of the PDF files are scaled, between the header and
# the footer
PLOT_SIDE_MARGIN = 36.0
PLOT_TOP_MARGIN = 44.0
PLOT_BOTTOM_MARGIN = 40.0

matplotlib.rcParams.update({
    # Embed the fonts as TrueType, so that the text can be selected
    "pdf.fonttype": 42,
    "font.family": SANS_FONT,
})

#############################################################################
# Content of the report
#############################################################################


@dataclass
class ReportFile:
    """A file of a case folder shown in the report."""
    folder: str
    name: str
    # Reader of a PDF file
    reader: object = None
    # Lines of every page of a text file
    text_pages: list = field(default_factory=list)
    # Why the file is not in the report; empty if it is
    problem: str = ""
    # Index of its first page in the report
    first_page: int = 0

    @property
    def label(self):
        return f"{self.folder}/{self.name}"

    @property
    def n_pages(self):
        if self.problem:
            return 0
        if self.reader is not None:
            return len(self.reader.pages)
        return len(self.text_pages)


@dataclass
class Case:
    """A validation case and the files shown for it."""
    path: str
    n_proc: str
    solver: str
    files: list = field(default_factory=list)
    # Index of its title page in the report
    first_page: int = 0

    @property
    def name(self):
        return PurePosixPath(self.path).name

    @property
    def folders(self):
        """Case folders of the files of the case, in order."""
        return list(dict.fromkeys(f.folder for f in self.files))

    @property
    def has_problem(self):
        return not self.files or any(f.problem for f in self.files)


def read_contents(path):
    """Read the cases of the report and the computer of the run."""
    computer = "unknown"
    cases = []
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for line_number, line in enumerate(lines, start=1):
        words = line.split()
        if not words or words[0].startswith("#"):
            continue
        keyword, arguments = words[0], words[1:]
        if keyword == "computer" and arguments:
            computer = " ".join(arguments)
        elif keyword == "case" and len(arguments) in (2, 3):
            solver = arguments[2] if len(arguments) == 3 else UNKNOWN_SOLVER
            cases.append(Case(arguments[0], arguments[1], solver))
        elif keyword == "file" and len(arguments) == 2:
            # Files listed before any case come from a validate.sh script
            # that was run by itself; they form a case named after their
            # folder.
            if not cases:
                cases.append(Case(arguments[0], "unknown", UNKNOWN_SOLVER))
            cases[-1].files.append(ReportFile(arguments[0], arguments[1]))
        else:
            print(f"Warning: line {line_number} of {path} is not "
                  f"understood and is ignored: {line}", file=sys.stderr)
    return computer, cases


def load_file(report_file, run_path):
    """Open a file of the report, or record why it cannot be shown."""
    path = run_path / report_file.folder / report_file.name
    if not path.is_file():
        report_file.problem = "not produced"
        return
    if path.suffix.lower() == ".pdf":
        try:
            reader = PdfReader(path)
            n_pages = len(reader.pages)
        # A broken plot, for example the output of a post-processing script
        # that crashed, must not prevent the rest of the report.
        except Exception as error:  # pylint: disable=broad-except
            report_file.problem = f"unreadable PDF ({error})"
            return
        if n_pages == 0:
            report_file.problem = "PDF without any page"
            return
        report_file.reader = reader
        return

    lines = []
    text = path.read_text(encoding="utf-8", errors="replace").expandtabs()
    for line in text.splitlines():
        lines.extend(textwrap.wrap(line, TEXT_COLUMNS, drop_whitespace=False)
                     or [""])
    if not "".join(lines).strip():
        report_file.problem = "empty file"
        return
    report_file.text_pages = [lines[i:i + TEXT_LINES_PER_PAGE] for i in
                              range(0, len(lines), TEXT_LINES_PER_PAGE)]


def group_by_solver(cases):
    """Group the consecutive cases that use the same solver."""
    return [(solver, list(group))
            for solver, group in itertools.groupby(cases,
                                                   key=lambda c: c.solver)]

#############################################################################
# Layout of the table of contents and of the title pages
#############################################################################


@dataclass
class ContentsRow:
    """A row of the table of contents: a solver heading or a case."""
    page: int
    column: int
    # Top of the row
    top: float
    heading: str = ""
    case: Case = None


def paginate_contents(cases):
    """Place the rows of the table of contents on its pages."""
    rows = []
    bottom = PAGE_HEIGHT - BOTTOM_MARGIN
    page = 0
    column = 0
    top = CONTENTS_FIRST_ROW

    def next_column():
        nonlocal page, column, top
        column += 1
        if column == CONTENTS_COLUMNS:
            column = 0
            page += 1
        top = CONTENTS_FIRST_ROW

    for solver, group in group_by_solver(cases):
        gap = CONTENTS_HEADING_GAP if top > CONTENTS_FIRST_ROW else 0.0
        # A heading is never left alone at the bottom of a column
        if (top + gap + CONTENTS_HEADING_HEIGHT + CONTENTS_CASE_HEIGHT >
                bottom):
            next_column()
            gap = 0.0
        top += gap
        rows.append(ContentsRow(page, column, top, heading=solver))
        top += CONTENTS_HEADING_HEIGHT
        for case in group:
            if top + CONTENTS_CASE_HEIGHT > bottom:
                next_column()
                rows.append(ContentsRow(page, column, top,
                                        heading=f"{solver} (continued)"))
                top += CONTENTS_HEADING_HEIGHT
            rows.append(ContentsRow(page, column, top, case=case))
            top += CONTENTS_CASE_HEIGHT
    return rows, page + 1


def contents_column_extent(column):
    """Horizontal extent of a column of the table of contents."""
    width = ((PAGE_WIDTH - 2.0 * SIDE_MARGIN -
              (CONTENTS_COLUMNS - 1) * CONTENTS_COLUMN_GAP) /
             CONTENTS_COLUMNS)
    left = SIDE_MARGIN + column * (width + CONTENTS_COLUMN_GAP)
    return left, left + width


def title_page_file_rows(case):
    """Baselines of the files listed on the title page of a case.

    Return the baselines of the listed files and, if they do not all fit on
    the page, the baseline of a last row that says how many are left out.
    """
    n_rows = int((PAGE_HEIGHT - BOTTOM_MARGIN - FILES_FIRST_BASELINE) /
                 FILE_LEADING) + 1
    baselines = [FILES_FIRST_BASELINE + i * FILE_LEADING
                 for i in range(min(n_rows, len(case.files)))]
    if len(case.files) > n_rows:
        return baselines[:-1], baselines[-1]
    return baselines, None


def link_rectangle(left, right, top, bottom):
    """Rectangle of a link, in PDF coordinates, from distances to the top."""
    return (left, PAGE_HEIGHT - bottom, right, PAGE_HEIGHT - top)

#############################################################################
# Drawing of the generated pages
#############################################################################


class Canvas:
    """A page drawn with matplotlib, addressed in points from its top left."""

    def __init__(self):
        self.figure = Figure(figsize=(PAGE_WIDTH / POINTS_PER_INCH,
                                      PAGE_HEIGHT / POINTS_PER_INCH),
                             facecolor=PAPER_COLOR)

    def text(self, x, y, string, size, color=TEXT_COLOR, weight="normal",
             font=SANS_FONT, align="left", masked=False):
        """Write a string whose baseline is at the distance y from the top."""
        options = {}
        # A masked string hides the dotted leaders drawn under it
        if masked:
            options["bbox"] = dict(boxstyle="square,pad=0.2",
                                   facecolor=PAPER_COLOR, edgecolor="none")
        self.figure.text(x / PAGE_WIDTH, 1.0 - y / PAGE_HEIGHT, string,
                         fontsize=size, color=color, fontweight=weight,
                         family=font, ha=align, va="baseline", **options)

    def rule(self, x_start, x_end, y, color=RULE_COLOR, width=0.6,
             dotted=False):
        """Draw a horizontal line at the distance y from the top."""
        style = (0, (0.6, 2.4)) if dotted else "-"
        self.figure.add_artist(Line2D(
            [x_start / PAGE_WIDTH, x_end / PAGE_WIDTH],
            [1.0 - y / PAGE_HEIGHT] * 2, transform=self.figure.transFigure,
            color=color, linewidth=width, linestyle=style,
            solid_capstyle="butt", dash_capstyle="round", zorder=1))

    def page_number(self, index, n_pages):
        self.text(PAGE_WIDTH - SIDE_MARGIN, PAGE_HEIGHT - FOOTER_BASELINE,
                  f"{index + 1} / {n_pages}", 8, MUTED_COLOR, align="right")

    def details(self, y, items):
        """Write rows of labels and values, and return the y below them."""
        for label, value, font in items:
            self.text(TITLE_LEFT, y, label, 9.5, MUTED_COLOR)
            self.text(TITLE_LEFT + DETAILS_LABEL_WIDTH, y, value,
                      9.5 if font == MONOSPACE_FONT else 10.5, font=font)
            y += DETAILS_LEADING
        return y


def draw_cover(canvas, git_hash, started, computer, cases, n_pages):
    canvas.text(TITLE_LEFT, TITLE_TOP, "Lethe", 14, ACCENT_COLOR,
                weight="bold")
    canvas.text(TITLE_LEFT, TITLE_TOP + 40, "Validation report", 32)
    canvas.rule(TITLE_LEFT, TITLE_RIGHT, TITLE_TOP + 64, ACCENT_COLOR,
                width=1.2)
    n_problems = sum(case.has_problem for case in cases)
    n_cases = f"{len(cases)}"
    if n_problems:
        n_cases += f", {n_problems} with missing output"
    canvas.details(TITLE_TOP + 98, [
        ("Started", started, SANS_FONT),
        ("Git hash", git_hash, MONOSPACE_FONT),
        ("Computer", computer, SANS_FONT),
        ("Cases", n_cases, SANS_FONT),
        ("Pages", f"{n_pages}", SANS_FONT),
    ])
    canvas.text(TITLE_LEFT, PAGE_HEIGHT - BOTTOM_MARGIN - 13,
                "Every case starts with a title page; the contents and the "
                "bookmarks link to them.", 8.5, MUTED_COLOR)
    canvas.text(TITLE_LEFT, PAGE_HEIGHT - BOTTOM_MARGIN,
                "Generated by contrib/validation/generate_validation_report.py"
                f" on {time.strftime('%Y-%m-%d %H:%M %Z')}.", 8.5,
                MUTED_COLOR)


def draw_contents_page(canvas, rows, page, n_pages):
    canvas.text(SIDE_MARGIN, CONTENTS_TITLE_BASELINE, "Contents", 22)
    for row in rows:
        if row.page != page:
            continue
        left, right = contents_column_extent(row.column)
        if row.heading:
            baseline = row.top + CONTENTS_HEADING_HEIGHT * 0.7
            canvas.text(left, baseline, row.heading, 10.5, ACCENT_COLOR,
                        weight="bold")
            continue
        case = row.case
        color = ALERT_COLOR if case.has_problem else TEXT_COLOR
        baseline = row.top + CONTENTS_CASE_HEIGHT * 0.72
        canvas.rule(left + CONTENTS_CASE_INDENT, right, baseline - 0.5,
                    MUTED_COLOR, width=0.8, dotted=True)
        if case.has_problem:
            canvas.text(left + 2.0, baseline - 0.5, PROBLEM_MARKER, 7,
                        ALERT_COLOR)
        canvas.text(left + CONTENTS_CASE_INDENT, baseline, case.name, 10,
                    color, masked=True)
        canvas.text(right, baseline, f"{case.first_page + 1}", 10, color,
                    align="right", masked=True)
    if any(row.case is not None and row.case.has_problem for row in rows):
        canvas.text(SIDE_MARGIN, PAGE_HEIGHT - FOOTER_BASELINE,
                    f"{PROBLEM_MARKER}  Missing output: a file of the case "
                    "was not produced or cannot be read; see its title page.",
                    8, ALERT_COLOR)
    canvas.page_number(1 + page, n_pages)


def draw_title_page(canvas, case, index, n_cases, n_pages):
    canvas.text(TITLE_LEFT, TITLE_TOP, case.solver, 12, ACCENT_COLOR,
                weight="bold")
    canvas.text(TITLE_LEFT, TITLE_TOP + 38, case.name, 28)
    canvas.text(TITLE_LEFT, TITLE_TOP + 62, case.path, 10, MUTED_COLOR,
                font=MONOSPACE_FONT)
    canvas.rule(TITLE_LEFT, TITLE_RIGHT, TITLE_TOP + 80, ACCENT_COLOR,
                width=1.2)
    folders = case.folders
    canvas.details(TITLE_TOP + 108, [
        ("MPI processes", case.n_proc, SANS_FONT),
        ("Output folders" if len(folders) > 1 else "Output folder",
         ", ".join(folders) if folders else "unknown", MONOSPACE_FONT),
    ])

    canvas.text(TITLE_LEFT, FILES_FIRST_BASELINE,
                "Files" if len(case.files) > 1 else "File", 9.5, MUTED_COLOR)
    if not case.files:
        canvas.text(TITLE_LEFT + DETAILS_LABEL_WIDTH, FILES_FIRST_BASELINE,
                    "Nothing was added to the report for this case; see its "
                    "logs.", 10.5, ALERT_COLOR, weight="bold")
    baselines, overflow_baseline = title_page_file_rows(case)
    number_right = TITLE_LEFT + DETAILS_LABEL_WIDTH + FILE_NUMBER_WIDTH
    name_left = number_right + 10.0
    for report_file, baseline in zip(case.files, baselines):
        if report_file.problem:
            canvas.text(number_right, baseline, "–", 9.5, ALERT_COLOR,
                        align="right")
            canvas.text(name_left, baseline,
                        f"{report_file.label}: {report_file.problem}", 9.5,
                        ALERT_COLOR, font=MONOSPACE_FONT)
        else:
            canvas.text(number_right, baseline,
                        f"{report_file.first_page + 1}", 9.5, MUTED_COLOR,
                        align="right")
            canvas.text(name_left, baseline, report_file.label, 9.5,
                        font=MONOSPACE_FONT)
    if overflow_baseline is not None:
        canvas.text(name_left, overflow_baseline,
                    f"and {len(case.files) - len(baselines)} more files", 9.5,
                    MUTED_COLOR)

    canvas.text(SIDE_MARGIN, PAGE_HEIGHT - FOOTER_BASELINE,
                f"Case {index + 1} of {n_cases}", 8, MUTED_COLOR)
    canvas.page_number(case.first_page, n_pages)


def draw_file_header(canvas, case, report_file, part, n_pages):
    """Draw the header and the footer of a page that shows a file."""
    canvas.text(SIDE_MARGIN, HEADER_BASELINE, case.name, 9, weight="bold")
    canvas.text(PAGE_WIDTH - SIDE_MARGIN, HEADER_BASELINE, case.solver, 9,
                ACCENT_COLOR, align="right")
    canvas.rule(SIDE_MARGIN, PAGE_WIDTH - SIDE_MARGIN, HEADER_BASELINE + 7)
    label = report_file.label
    if report_file.n_pages > 1:
        label += f"  (page {part + 1} of {report_file.n_pages})"
    canvas.text(SIDE_MARGIN, PAGE_HEIGHT - FOOTER_BASELINE, label, 8,
                MUTED_COLOR, font=MONOSPACE_FONT)
    canvas.page_number(report_file.first_page + part, n_pages)


def draw_text_page(canvas, case, report_file, part, n_pages):
    draw_file_header(canvas, case, report_file, part, n_pages)
    y = TEXT_FIRST_LINE
    for line in report_file.text_pages[part]:
        canvas.text(SIDE_MARGIN, y, line, TEXT_SIZE, font=MONOSPACE_FONT)
        y += TEXT_LEADING

#############################################################################
# Assembly of the report
#############################################################################


def place_on_page(page, source):
    """Scale a page of a PDF file onto a page of the report and centre it."""
    if source.rotation:
        source.transfer_rotation_to_content()
    box = source.mediabox
    width = float(box.width)
    height = float(box.height)
    box_width = PAGE_WIDTH - 2.0 * PLOT_SIDE_MARGIN
    box_height = PAGE_HEIGHT - PLOT_TOP_MARGIN - PLOT_BOTTOM_MARGIN
    scale = min(box_width / width, box_height / height)
    x = PLOT_SIDE_MARGIN + 0.5 * (box_width - scale * width)
    y = PLOT_BOTTOM_MARGIN + 0.5 * (box_height - scale * height)
    page.merge_transformed_page(
        source, Transformation().translate(-float(box.left),
                                           -float(box.bottom))
        .scale(scale, scale).translate(x, y))


def add_link(writer, page, rectangle, target_page):
    writer.add_annotation(page_number=page,
                          annotation=Link(rect=rectangle,
                                          target_page_index=target_page,
                                          border=[0, 0, 0]))


def generate_report(run_path, output_path):
    """Generate the report of a run and return its number of pages and cases.
    """
    computer, cases = read_contents(run_path / CONTENTS_NAME)
    git_hash = "unknown"
    started = "unknown"
    hash_path = run_path / HASH_FILE
    if hash_path.is_file():
        git_hash = hash_path.read_text().strip()
        started = time.strftime("%Y-%m-%d %H:%M %Z",
                                time.localtime(hash_path.stat().st_mtime))

    for case in cases:
        for report_file in case.files:
            load_file(report_file, run_path)

    # Number the pages: cover, contents, then every case
    rows, n_contents_pages = paginate_contents(cases)
    page_index = 1 + n_contents_pages
    for case in cases:
        case.first_page = page_index
        page_index += 1
        for report_file in case.files:
            report_file.first_page = page_index
            page_index += report_file.n_pages
    n_pages = page_index

    # Draw every page of the report, in order, and keep the page of a PDF
    # file to place on each of them
    generated = io.BytesIO()
    sources = []
    with PdfPages(generated) as pdf:

        def add(draw, *arguments, source=None):
            canvas = Canvas()
            draw(canvas, *arguments)
            pdf.savefig(canvas.figure)
            sources.append(source)

        add(draw_cover, git_hash, started, computer, cases, n_pages)
        for page in range(n_contents_pages):
            add(draw_contents_page, rows, page, n_pages)
        for index, case in enumerate(cases):
            add(draw_title_page, case, index, len(cases), n_pages)
            for report_file in case.files:
                for part in range(report_file.n_pages):
                    if report_file.reader is None:
                        add(draw_text_page, case, report_file, part, n_pages)
                    else:
                        add(draw_file_header, case, report_file, part,
                            n_pages, source=report_file.reader.pages[part])

    writer = PdfWriter()
    for page, source in zip(PdfReader(generated).pages, sources):
        page = writer.add_page(page)
        if source is not None:
            place_on_page(page, source)

    # Links of the contents to the title pages, and of the title pages to the
    # files
    for row in rows:
        if row.case is not None:
            left, right = contents_column_extent(row.column)
            add_link(writer, 1 + row.page,
                     link_rectangle(left, right, row.top,
                                    row.top + CONTENTS_CASE_HEIGHT),
                     row.case.first_page)
    for case in cases:
        baselines, _ = title_page_file_rows(case)
        for report_file, baseline in zip(case.files, baselines):
            if report_file.n_pages:
                add_link(writer, case.first_page,
                         link_rectangle(TITLE_LEFT + DETAILS_LABEL_WIDTH,
                                        TITLE_RIGHT,
                                        baseline - 0.75 * FILE_LEADING,
                                        baseline + 0.25 * FILE_LEADING),
                         report_file.first_page)

    # Bookmarks: solver, case and, for a case with several folders, folder
    writer.add_outline_item(title="Contents", page_number=1)
    for solver, group in group_by_solver(cases):
        parent = writer.add_outline_item(title=solver,
                                         page_number=group[0].first_page)
        for case in group:
            item = writer.add_outline_item(title=case.name,
                                           page_number=case.first_page,
                                           parent=parent, is_open=False)
            if len(case.folders) < 2:
                continue
            for folder in case.folders:
                pages = [f.first_page for f in case.files
                         if f.folder == folder and f.n_pages]
                if pages:
                    writer.add_outline_item(title=folder, page_number=pages[0],
                                            parent=item)
    writer.page_mode = "/UseOutlines"
    writer.add_metadata({
        "/Title": "Lethe validation report",
        "/Subject": f"Git hash {git_hash}",
        "/Creator": "contrib/validation/generate_validation_report.py",
    })

    # Write the report next to its final location and then move it, so that
    # an interrupted generation never leaves a broken report.
    temporary_path = output_path.with_name(output_path.name + ".tmp")
    with open(temporary_path, "wb") as stream:
        writer.write(stream)
    os.replace(temporary_path, output_path)
    return n_pages, len(cases)

#############################################################################
# Command line
#############################################################################


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate the PDF report of a Lethe validation run from "
                    f"the list of files in {CONTENTS_NAME}, written in the "
                    "run folder by validate_lethe.sh and add_to_report.")
    parser.add_argument(
        "-r", "--run", type=Path, required=True,
        help="Folder of the validation run.")
    parser.add_argument(
        "-o", "--output", type=Path,
        help=f"Path of the report. Defaults to {REPORT_NAME} in the run "
             "folder.")
    return parser.parse_args()


def main():
    args = parse_arguments()
    contents_path = args.run / CONTENTS_NAME
    if not contents_path.is_file():
        sys.exit(f"Error: {contents_path} does not exist. It is written by "
                 "validate_lethe.sh while the cases run.")
    output_path = (args.output or args.run / REPORT_NAME).resolve()
    n_pages, n_cases = generate_report(args.run.resolve(), output_path)
    print(f"Report written to {output_path} ({n_cases} cases, {n_pages} "
          "pages)")


if __name__ == "__main__":
    main()
