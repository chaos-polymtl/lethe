#!/usr/bin/env python3
"""
Combine dated release-note markdown files and extract MAJOR entries.

Expects a folder of files named like:  yyyy_mm_dd_name.md

Produces two files:
  1. combined_release_notes.md
     The full text of every release note, most recent date first, with
     same-day notes merged under one date heading.

  2. major_changes.md
     Every bullet point that starts with the MAJOR keyword (e.g.
     "- MAJOR This PR ..."), grouped under its release date, most
     recent first.

Usage:
    python generate_release_summary.py -i /path/to/release_notes_folder
    python generate_release_summary.py -i /path/to/folder -o /path/to/output_folder
"""

import argparse
import re
from collections import OrderedDict
from datetime import date
from itertools import groupby
from pathlib import Path


# Matches "..._2026_06_16_Blais_2.md" -> groups: 2026, 06, 16, Blais_2
FILENAME_DATE_RE = re.compile(r"(\d{4})_(\d{2})_(\d{2})_(.+)\.md$", re.IGNORECASE)

# Matches a bullet line that starts with the MAJOR keyword, e.g.:
#   "- MAJOR This PR adds ..."
# This is case-insensitive and discards the white space before and in between the bullet and the keyword.
MAJOR_BULLET_RE = re.compile(r"^\s*[-*]\s*MAJOR\b", re.IGNORECASE)

# Matches a top-level "## [Master] - 2026/06/16" style date heading.
TOP_HEADING_RE = re.compile(r"^##\s")

# Matches a "### Added" style section heading.
SECTION_HEADING_RE = re.compile(r"^###\s+(.*)")


 
def parse_release_file(path: Path):
    """Return {date, name, text} for a release note file, or None if the
    filename doesn't match the expected yyyy_mm_dd_name.md pattern."""
    match = FILENAME_DATE_RE.search(path.name)
    if not match:
        return None
    year, month, day, name = match.groups()
    try:
        file_date = date(int(year), int(month), int(day))
    except ValueError:
        return None
    text = path.read_text(encoding="utf-8")
    return {"date": file_date, "name": name, "text": text, "path": path}
 
 
def collect_release_notes(folder: Path):
    """Find, parse, and date-sort (newest first) every release note in folder."""
    notes = []
    skipped = []
    for path in sorted(folder.glob("*.md")):
        parsed = parse_release_file(path)
        if parsed is None:
            skipped.append(path.name)
            continue
        notes.append(parsed)
 
    if skipped:
        print("Skipped (filename doesn't match yyyy_mm_dd_name.md):")
        for name in skipped:
            print(f"  - {name}")
 
    notes.sort(key=lambda n: n["date"], reverse=True)
    return notes
 
 
def group_by_date(notes):
    """Group already date-sorted notes into (date, [notes...]) pairs,
    preserving the overall newest-first order."""
    return [(file_date, list(group)) for file_date, group in groupby(notes, key=lambda n: n["date"])]
 
 
def strip_blank_edges(lines):
    """Drop leading/trailing blank lines from a list of lines."""
    start, end = 0, len(lines)
    while start < end and lines[start].strip() == "":
        start += 1
    while end > start and lines[end - 1].strip() == "":
        end -= 1
    return lines[start:end]
 
 
def parse_sections(text: str):
    """Split a release note's body into an ordered list of (title, content_lines),
    skipping the top-level "## ..." date heading."""
    lines = text.splitlines()
    idx = 0
    while idx < len(lines) and lines[idx].strip() == "":
        idx += 1
    if idx < len(lines) and TOP_HEADING_RE.match(lines[idx]):
        idx += 1
 
    sections = []
    current_title = None
    current_content = []
    for line in lines[idx:]:
        match = SECTION_HEADING_RE.match(line)
        if match:
            if current_title is not None:
                sections.append((current_title, strip_blank_edges(current_content)))
            current_title = match.group(1).strip()
            current_content = []
        else:
            current_content.append(line)
    if current_title is not None:
        sections.append((current_title, strip_blank_edges(current_content)))
    return sections
 
 
def merge_sections(notes_for_date):
    """Merge the sections of several same-date notes into one ordered
    list of (title, merged_content_lines), combining matching titles."""
    merged = OrderedDict()
    for note in notes_for_date:
        for title, content in parse_sections(note["text"]):
            merged.setdefault(title, []).extend(content)
    return list(merged.items())
 
 
def format_merged_note(file_date, sections):
    lines = [f"## [Master] - {file_date.strftime('%Y/%m/%d')}", ""]
    for title, content in sections:
        lines.append(f"### {title}")
        lines.append("")
        lines.extend(content)
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"
 
 
def extract_major_bullets(text: str):
    """Return every bullet line in text that starts with the MAJOR keyword."""
    return [line.strip() for line in text.splitlines() if MAJOR_BULLET_RE.match(line)]
 
 
def write_combined_notes(notes, output_path: Path):
    parts = []
    for file_date, notes_for_date in group_by_date(notes):
        sections = merge_sections(notes_for_date)
        parts.append(format_merged_note(file_date, sections).rstrip())
    output_path.write_text("\n\n---\n\n".join(parts) + "\n", encoding="utf-8")
 
 
def write_major_changes(notes, output_path: Path):
    lines = ["# Major Changes", ""]
    found_any = False
    for file_date, notes_for_date in group_by_date(notes):
        bullets = []
        for note in notes_for_date:
            bullets.extend(extract_major_bullets(note["text"]))
        if not bullets:
            continue
        found_any = True
        lines.append(f"## {file_date.strftime('%Y/%m/%d')}")
        lines.append("")
        lines.extend(bullets)
        lines.append("")
    if not found_any:
        lines.append("No MAJOR entries found.")
    output_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")



def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input-folder", help="Folder containing the release note .md files")
    parser.add_argument(
        "-o", "--output-folder", default=None,
        help="Where to write the two output files (default: same as input folder)",
    )
    args = parser.parse_args()

    input_folder = Path(args.input_folder).expanduser().resolve()
    if not input_folder.is_dir():
        raise SystemExit(f"Not a folder: {input_folder}")

    output_folder = Path(args.output_folder).expanduser().resolve() if args.output_folder else input_folder
    output_folder.mkdir(parents=True, exist_ok=True)

    notes = collect_release_notes(input_folder)
    if not notes:
        raise SystemExit("No matching release note files found (expected yyyy_mm_dd_name.md).")

    combined_path = output_folder / "combined_release_notes.md"
    major_path = output_folder / "major_changes.md"

    write_combined_notes(notes, combined_path)
    write_major_changes(notes, major_path)

    print(f"Processed {len(notes)} release note file(s), newest to oldest:")
    for note in notes:
        print(f"  - {note['date'].isoformat()}  {note['path'].name}")
    print(f"\nCombined notes written to: {combined_path}")
    print(f"Major changes written to:  {major_path}")


if __name__ == "__main__":
    main()