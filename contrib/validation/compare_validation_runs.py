#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

"""
Compare a Lethe validation run with the two previous runs and report the result.

Every case folder of a validation run (see validate_lethe.sh) contains data
files (*.dat, *.txt), plots (*.pdf) and logs (log*). For every case, this
script:
  - compares the data files with those of the two most recent previous runs
    that contain the same case, first byte by byte and then value by value;
  - checks that the case produced at least one data file and one plot;
  - compares the wall time recorded by `time` at the end of the logs.

The previous runs are the folders next to the current run that contain a
current_git_hash.txt file, ordered by the modification time of that file.

The result is written as a Markdown report in the run folder and printed on
the standard output. The diffs of the data files that changed since the
previous run are written next to the report, in comparison_diffs/<case>.diff.
If the environment variable LETHE_VALIDATION_SLACK_WEBHOOK contains the URL of
a Slack incoming webhook, a summary of every case, with a legend of the
symbols it uses, is also posted to Slack. It is followed by one message per
case that changed or has a problem, which gives an excerpt of the diff of
every changed file and, for a case with a problem, the end of its section of
the Slurm log (slurm_log.txt in the run folder).

Usage:
  python3 compare_validation_runs.py --current <run folder> [options]

Run with --help for the list of options.
"""

import argparse
import difflib
import fnmatch
import json
import math
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

#############################################################################
# Constants
#############################################################################

# File written by validate_lethe.sh at the start of every run
HASH_FILE = "current_git_hash.txt"

# Files of a case folder that hold the validation data, plots and logs
DATA_PATTERNS = ("*.dat", "*.txt")
PLOT_PATTERN = "*.pdf"
LOG_PATTERN = "log*"

REPORT_NAME = "comparison_report.md"
SLACK_WEBHOOK_VARIABLE = "LETHE_VALIDATION_SLACK_WEBHOOK"
GITHUB_COMPARE_URL = "https://github.com/chaos-polymtl/lethe/compare/{}...{}"
GIT_HASH_REGEX = re.compile(r"^[0-9a-f]{7,40}$")
SHORT_HASH_LENGTH = 7

# Values smaller than this fraction of the largest value of a file are
# compared in absolute terms, so that round-off around zero does not produce
# huge relative differences.
RELATIVE_FLOOR = 1e-12

# Wall-time changes are only flagged when the case ran for at least this long
# in the current or in the previous run, since short timings are mostly noise.
MIN_TIMED_SECONDS = 60.0

# Only the end of each log is read to find the wall time; some logs are
# hundreds of MB.
LOG_TAIL_BYTES = 4096

# Wall time printed by the bash `time` keyword, e.g. "real	82m43.277s"
WALL_TIME_REGEX = re.compile(r"^real\s+(?:(\d+)h)?(\d+)m([\d.]+)s\s*$",
                             re.MULTILINE)

MAX_COMMITS_LISTED = 15
GIT_TIMEOUT_SECONDS = 60
SLACK_TIMEOUT_SECONDS = 30

# Folder, next to the report, where the diffs of the changed cases are written
DIFF_FOLDER_NAME = "comparison_diffs"

# Output of the Slurm job, copied into the run folder by the job script
SLURM_LOG_NAME = "slurm_log.txt"

# Lines written by validate_lethe.sh before and after every case
PROCESSING_REGEX = re.compile(r"^---> Processing (\S+) using")
FINISHED_REGEX = re.compile(r"^\s*Finished processing")

# Size of the end of the Slurm log shown for a case with a problem
MAX_LOG_EXCERPT_LINES = 30
MAX_LOG_LINE_LENGTH = 300

# Slack recommends messages of at most 4000 characters, and incoming webhooks
# accept about one message per second
MAX_SLACK_DIFF_LINES = 20
MAX_SLACK_MESSAGE_LENGTH = 4000
SLACK_SECONDS_BETWEEN_MESSAGES = 1.0

# Status of one data file compared with the same file of a reference run
IDENTICAL = "identical"
SAME_VALUES = "same values"
CHANGED = "changed"
STRUCTURE_CHANGED = "structure changed"
MISSING = "missing"
NEW = "new"

# Symbol of the cases whose wall time changed significantly
WALL_TIME_SYMBOL = "⏱️"

#############################################################################
# Data structures
#############################################################################


@dataclass(frozen=True)
class Verdict:
    """How a verdict is shown in the reports."""
    symbol: str
    # Label used in the summary and on Slack
    label: str
    # Description used in the headings of the details of the report
    description: str
    # Explanation of the verdict given in the legend of the reports
    meaning: str


# Verdict for a case, from the most to the least severe
VERDICTS = {
    "problem": Verdict(
        "❌", "with problems", "problem",
        "the case is missing from this run, produced no data file or no "
        "plot, or one of its data files is missing or has a different "
        "structure than in the previous run"),
    "changed": Verdict(
        "🔶", "changed", "changed in this run",
        "the results differ from the previous run"),
    "reverted": Verdict(
        "⚠️", "reverted",
        "different from the previous run but identical to the one before it "
        "(non-deterministic results or reverted change?)",
        "the results differ from the previous run but are identical to the "
        "run before it (non-deterministic results or a reverted change?)"),
    "stabilized": Verdict(
        "ℹ️", "changed in the previous run",
        "changed in the previous run, unchanged since",
        "the results are the same as in the previous run, but that run "
        "differed from the one before it (the change was already reported "
        "last time)"),
    "new": Verdict(
        "🆕", "new", "new case, no previous run to compare with",
        "no previous run contains this case, so there is nothing to compare "
        "with"),
    "unchanged": Verdict(
        "✅", "unchanged", "unchanged",
        "the results are the same as in the previous runs"),
}
SEVERITY = list(VERDICTS)

# Verdicts of the cases that are detailed in a message of their own
DETAILED_VERDICTS = ("problem", "changed", "reverted")


@dataclass
class Run:
    """A validation run: one folder of the archive."""
    path: Path
    git_hash: str
    start_time: float
    cases: dict = field(default_factory=dict)

    @property
    def name(self):
        return self.path.name

    @property
    def short_hash(self):
        return self.git_hash[:SHORT_HASH_LENGTH] if self.git_hash else "unknown"


@dataclass
class FileComparison:
    """Comparison of one data file with the same file of a reference run."""
    name: str
    status: str
    max_relative_difference: float = 0.0
    line: int = 0
    detail: str = ""
    # Paths of the file in the reference and in the current run, if it exists
    reference_path: Path = None
    current_path: Path = None

    def describe(self):
        if self.status == CHANGED:
            difference = format_relative_difference(
                self.max_relative_difference)
            return (f"{self.name}: changed, max. relative difference "
                    f"{difference} (line {self.line})")
        if self.detail:
            return f"{self.name}: {self.status}, {self.detail}"
        return f"{self.name}: {self.status}"


@dataclass
class ReferenceComparison:
    """Comparison of all the data files of a case with a reference run."""
    run: Run
    files: list

    @property
    def outcome(self):
        """Return "same", "changed" or "problem"."""
        statuses = {f.status for f in self.files}
        if statuses & {STRUCTURE_CHANGED, MISSING}:
            return "problem"
        if statuses & {CHANGED, NEW}:
            return "changed"
        return "same"

    @property
    def max_relative_difference(self):
        return max((f.max_relative_difference for f in self.files
                    if f.status == CHANGED), default=0.0)

    def summary(self):
        if not self.files:
            return "no data"
        outcome = self.outcome
        if outcome == "same":
            if all(f.status == IDENTICAL for f in self.files):
                return "identical"
            return "same values"
        if outcome == "changed":
            if self.max_relative_difference > 0.0:
                difference = format_relative_difference(
                    self.max_relative_difference)
                return f"changed (max. rel. diff. {difference})"
            return "changed (new files)"
        return ", ".join(sorted({f.status for f in self.files
                                 if f.status in (STRUCTURE_CHANGED, MISSING)}))

    def differences(self):
        return [f for f in self.files if f.status != IDENTICAL]

    def changes(self):
        """Return the files whose values differ from the reference."""
        return [f for f in self.files
                if f.status not in (IDENTICAL, SAME_VALUES)]


@dataclass
class CaseResult:
    """Everything that is reported for one case."""
    name: str
    problems: list
    comparisons: list
    wall_time: float = None
    previous_wall_time: float = None
    wall_time_flagged: bool = False
    verdict: str = "unchanged"
    # Diff of every changed data file with the previous run, as a dict
    # file name -> lines
    diffs: dict = field(default_factory=dict)
    # End of the section of the Slurm log of a case with a problem
    log_excerpt: list = field(default_factory=list)

#############################################################################
# Discovery of the runs and cases
#############################################################################


def read_git_hash(run_path):
    """Return the git hash stored in a run folder, or "" if it is unusable."""
    try:
        git_hash = (run_path / HASH_FILE).read_text().strip()
    except OSError:
        return ""
    return git_hash if GIT_HASH_REGEX.match(git_hash) else ""


def find_cases(run_path):
    """Return the case folders of a run, as a dict name -> path.

    A case folder is a sub-folder that contains a log, since every validate.sh
    writes the output of the simulation into one.
    """
    cases = {}
    for path in sorted(run_path.iterdir()):
        if path.is_dir() and any(p.is_file() for p in path.glob(LOG_PATTERN)):
            cases[path.name] = path
    return cases


def load_run(run_path):
    """Load a run folder. The start time is the time of its hash file."""
    hash_file = run_path / HASH_FILE
    start_time = hash_file.stat().st_mtime if hash_file.exists() else 0.0
    return Run(path=run_path, git_hash=read_git_hash(run_path),
               start_time=start_time, cases=find_cases(run_path))


def find_previous_runs(current, archive, pattern):
    """Return the runs of the archive that started before the current run.

    The runs are sorted from the most recent to the oldest. Only the folders
    whose name matches the glob pattern are considered.
    """
    runs = []
    for path in archive.iterdir():
        if (path.resolve() == current.path.resolve() or
                not path.is_dir() or
                not fnmatch.fnmatch(path.name, pattern) or
                not (path / HASH_FILE).is_file()):
            continue
        run = load_run(path)
        if run.start_time < current.start_time:
            runs.append(run)
    return sorted(runs, key=lambda run: run.start_time, reverse=True)

#############################################################################
# Comparison of the data files
#############################################################################


def tokenize(content):
    """Split a text file into tokens.

    Returns a list of (line number, token, value) where value is the token as
    a float, or None if the token is not a number.
    """
    tokens = []
    for line_number, line in enumerate(content.splitlines(), start=1):
        for token in line.split():
            try:
                value = float(token)
            except ValueError:
                value = None
            tokens.append((line_number, token, value))
    return tokens


def relative_difference(a, b, floor):
    """Return |a-b| / max(|a|, |b|, floor), handling non-finite values."""
    if a == b or (math.isnan(a) and math.isnan(b)):
        return 0.0
    if not (math.isfinite(a) and math.isfinite(b)):
        return math.inf
    return abs(a - b) / max(abs(a), abs(b), floor)


def format_relative_difference(difference):
    """Return a relative difference as a percentage, e.g. "1.09%"."""
    if not math.isfinite(difference):
        return "non-finite"
    return f"{100.0 * difference:.3g}%"


def compare_data_file(name, current_path, reference_path):
    """Compare a data file of the current run with the one of a reference."""
    current_bytes = current_path.read_bytes()
    reference_bytes = reference_path.read_bytes()
    if current_bytes == reference_bytes:
        return FileComparison(name, IDENTICAL)

    current_tokens = tokenize(current_bytes.decode("utf-8", errors="replace"))
    reference_tokens = tokenize(
        reference_bytes.decode("utf-8", errors="replace"))
    if len(current_tokens) != len(reference_tokens):
        return FileComparison(
            name, STRUCTURE_CHANGED,
            detail=f"{len(current_tokens)} entries instead of "
                   f"{len(reference_tokens)}")

    magnitudes = [abs(value) for _, _, value in current_tokens + reference_tokens
                  if value is not None and math.isfinite(value)]
    floor = RELATIVE_FLOOR * max(magnitudes, default=0.0)

    result = FileComparison(name, SAME_VALUES)
    for (line, token, value), (_, reference_token, reference_value) in zip(
            current_tokens, reference_tokens):
        if value is None or reference_value is None:
            if token != reference_token:
                return FileComparison(
                    name, STRUCTURE_CHANGED,
                    detail=f"line {line}: '{token}' instead of "
                           f"'{reference_token}'")
            continue
        difference = relative_difference(value, reference_value, floor)
        if difference > result.max_relative_difference:
            result.max_relative_difference = difference
            result.line = line

    return result


def data_files(case_path):
    """Return the data files of a case folder, as a dict name -> path."""
    files = {}
    for pattern in DATA_PATTERNS:
        for path in case_path.glob(pattern):
            if path.is_file():
                files[path.name] = path
    return files


def compare_case_with_reference(case_path, reference, case_name, rtol):
    """Compare all the data files of a case with the ones of a reference."""
    current_files = data_files(case_path)
    reference_files = data_files(reference.cases[case_name])

    files = []
    for name in sorted(current_files.keys() | reference_files.keys()):
        if name not in reference_files:
            comparison = FileComparison(name, NEW)
        elif name not in current_files:
            comparison = FileComparison(name, MISSING)
        else:
            comparison = compare_data_file(name, current_files[name],
                                           reference_files[name])
            if (comparison.status == SAME_VALUES and
                    comparison.max_relative_difference > rtol):
                comparison.status = CHANGED
        comparison.reference_path = reference_files.get(name)
        comparison.current_path = current_files.get(name)
        files.append(comparison)
    return ReferenceComparison(run=reference, files=files)

#############################################################################
# Diffs and Slurm log
#############################################################################


def diff_lines(reference_path, current_path, reference_label, current_label):
    """Return the unified diff, without context, of two versions of a file.

    When both versions have the same number of lines, as when only values
    changed, every changed line gets a hunk of its own, so that each old line
    is directly followed by the new one. Otherwise, the usual unified diff is
    returned. In both cases, the diff can be applied with patch.
    """
    reference = reference_path.read_text(encoding="utf-8",
                                         errors="replace").splitlines()
    current = current_path.read_text(encoding="utf-8",
                                     errors="replace").splitlines()
    if len(reference) != len(current):
        return list(difflib.unified_diff(reference, current, reference_label,
                                         current_label, n=0, lineterm=""))

    lines = []
    for number, (old, new) in enumerate(zip(reference, current), start=1):
        if old != new:
            lines += [f"@@ -{number} +{number} @@", f"-{old}", f"+{new}"]
    if not lines:
        return []
    return [f"--- {reference_label}", f"+++ {current_label}"] + lines


def case_diffs(case, current):
    """Return the diffs of the changed data files of a case.

    The files are compared with the previous run that contains the case. The
    result is a dict file name -> lines of the diff.
    """
    if not case.comparisons:
        return {}
    reference = case.comparisons[0].run
    diffs = {}
    for comparison in case.comparisons[0].files:
        if comparison.status not in (CHANGED, STRUCTURE_CHANGED):
            continue
        lines = diff_lines(
            comparison.reference_path, comparison.current_path,
            f"{reference.name}/{case.name}/{comparison.name}",
            f"{current.name}/{case.name}/{comparison.name}")
        if lines:
            diffs[comparison.name] = lines
    return diffs


def diff_excerpt(lines, max_lines):
    """Return the beginning of a diff and the number of lines left out.

    The diff is cut between two hunks, so that an old line is never shown
    without the new one, unless its first hunk alone is too long.
    """
    if len(lines) <= max_lines:
        return lines, 0
    hunk_starts = [index for index, line in enumerate(lines)
                   if line.startswith("@@ ")]
    fitting = [start for start in hunk_starts[1:] if start <= max_lines]
    end = fitting[-1] if fitting else max_lines
    return lines[:end], len(lines) - end


def read_slurm_log_sections(run_path):
    """Return the lines of the Slurm log written for every case of a run.

    The result is a dict case name -> lines, where the case name is the last
    component of the path printed by validate_lethe.sh. Only the text after
    the last carriage return of a line is kept, as a terminal shows it, so
    that the progress animation of validate.sh takes a single line.
    Returns an empty dict if the run has no Slurm log.
    """
    # The bytes are decoded explicitly, since reading in text mode would turn
    # the carriage returns into line breaks
    try:
        content = (run_path / SLURM_LOG_NAME).read_bytes().decode(
            "utf-8", errors="replace")
    except OSError:
        return {}

    sections = {}
    section = None
    for line in content.split("\n"):
        frames = [frame for frame in line.split("\r") if frame.strip()]
        line = frames[-1].rstrip() if frames else ""
        match = PROCESSING_REGEX.match(line)
        if match:
            section = sections[Path(match.group(1)).name] = []
        if section is None:
            continue
        section.append(line)
        if FINISHED_REGEX.match(line):
            section = None
    return sections


def log_excerpt(sections, case_name):
    """Return the end of the section of the Slurm log of a case.

    A case written in a folder with another name than its example, such as
    3d-rectangular-hopper-periodic, gets the section of the example whose
    name is the longest prefix of its own.
    """
    names = [name for name in sections
             if case_name == name or case_name.startswith(name + "-")]
    if not names:
        return []
    lines = sections[max(names, key=len)][-MAX_LOG_EXCERPT_LINES:]
    return [line if len(line) <= MAX_LOG_LINE_LENGTH
            else line[:MAX_LOG_LINE_LENGTH] + "…" for line in lines]


def write_diffs(cases, folder):
    """Write the diffs of every case into folder/<case>.diff.

    The diffs of a previous comparison are removed first. Returns a dict
    case name -> path of its diff.
    """
    if folder.is_dir():
        for old_diff in folder.glob("*.diff"):
            old_diff.unlink()
    paths = {}
    for case in cases:
        if not case.diffs:
            continue
        folder.mkdir(exist_ok=True)
        path = folder / f"{case.name}.diff"
        path.write_text("".join("\n".join(lines) + "\n"
                                for lines in case.diffs.values()),
                        encoding="utf-8")
        paths[case.name] = path
    return paths

#############################################################################
# Wall time
#############################################################################


def read_wall_time(case_path):
    """Return the total wall time in seconds recorded in the logs of a case.

    Returns None if no log contains a wall time, for instance when the
    simulation was killed.
    """
    total = None
    for log in sorted(case_path.glob(LOG_PATTERN)):
        if not log.is_file():
            continue
        with open(log, "rb") as stream:
            stream.seek(0, os.SEEK_END)
            stream.seek(max(0, stream.tell() - LOG_TAIL_BYTES))
            tail = stream.read().decode("utf-8", errors="replace")
        for match in WALL_TIME_REGEX.finditer(tail):
            hours, minutes, seconds = match.groups()
            total = ((total or 0.0) + 3600.0 * int(hours or 0) +
                     60.0 * int(minutes) + float(seconds))
    return total


def format_duration(seconds):
    if seconds is None:
        return "n/a"
    if seconds < 60.0:
        return f"{seconds:.1f} s"
    minutes, seconds = divmod(round(seconds), 60)
    if minutes < 60:
        return f"{minutes}m {seconds:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes:02d}m"


def wall_time_change(case):
    """Return the relative change of wall time, or None if unknown."""
    if case.wall_time is None or not case.previous_wall_time:
        return None
    return case.wall_time / case.previous_wall_time - 1.0

#############################################################################
# Comparison of the cases
#############################################################################


def find_problems(case_path):
    """Return the reasons why a case of the current run is incomplete."""
    problems = []
    if not data_files(case_path):
        problems.append("no data file")
    if not any(p.is_file() for p in case_path.glob(PLOT_PATTERN)):
        problems.append("no plot")
    return problems


def decide_verdict(case):
    if case.problems:
        return "problem"
    if not case.comparisons:
        return "new"
    latest = case.comparisons[0].outcome
    before = (case.comparisons[1].outcome if len(case.comparisons) > 1
              else None)
    if latest == "problem":
        return "problem"
    if latest == "same":
        return "unchanged" if before in (None, "same") else "stabilized"
    return "reverted" if before == "same" else "changed"


def compare_cases(current, previous_runs, rtol, time_threshold):
    """Compare every case of the current run with the previous runs."""
    names = set(current.cases)
    if previous_runs:
        names |= set(previous_runs[0].cases)
    log_sections = read_slurm_log_sections(current.path)

    results = []
    for name in sorted(names):
        references = [run for run in previous_runs if name in run.cases][:2]
        if name not in current.cases:
            case = CaseResult(name=name,
                              problems=["case missing from this run"],
                              comparisons=[], verdict="problem")
        else:
            case_path = current.cases[name]
            case = CaseResult(
                name=name,
                problems=find_problems(case_path),
                comparisons=[compare_case_with_reference(case_path, reference,
                                                         name, rtol)
                             for reference in references],
                wall_time=read_wall_time(case_path))
            if references:
                case.previous_wall_time = read_wall_time(
                    references[0].cases[name])
            change = wall_time_change(case)
            case.wall_time_flagged = (
                change is not None and abs(change) > time_threshold and
                max(case.wall_time, case.previous_wall_time) >
                MIN_TIMED_SECONDS)
            case.verdict = decide_verdict(case)

        if case.verdict in DETAILED_VERDICTS:
            case.diffs = case_diffs(case, current)
        if case.verdict == "problem":
            case.log_excerpt = log_excerpt(log_sections, name)
        results.append(case)
    return results

#############################################################################
# Git information
#############################################################################


def compare_url(old_run, new_run):
    if not (old_run.git_hash and new_run.git_hash):
        return None
    return GITHUB_COMPARE_URL.format(old_run.short_hash, new_run.short_hash)


def list_commits(repository, old_run, new_run):
    """Return the subjects of the commits between two runs, or None."""
    if not (repository and old_run.git_hash and new_run.git_hash):
        return None
    try:
        process = subprocess.run(
            ["git", "-C", str(repository), "log", "--oneline", "--no-decorate",
             f"{old_run.git_hash}..{new_run.git_hash}"],
            capture_output=True, text=True, check=True,
            timeout=GIT_TIMEOUT_SECONDS)
    except (OSError, subprocess.SubprocessError):
        return None
    return process.stdout.splitlines()

#############################################################################
# Reports
#############################################################################


def count_verdicts(cases):
    counts = {verdict: 0 for verdict in SEVERITY}
    for case in cases:
        counts[case.verdict] += 1
    return counts


def summary_line(cases):
    counts = count_verdicts(cases)
    parts = [f"{counts[verdict]} {VERDICTS[verdict].label}"
             for verdict in reversed(SEVERITY) if counts[verdict]]
    n_flagged = sum(case.wall_time_flagged for case in cases)
    if n_flagged:
        parts.append(f"{n_flagged} wall-time change"
                     f"{'s' if n_flagged > 1 else ''}")
    return " · ".join(parts)


def reference_cell(case, index, default_reference):
    if index >= len(case.comparisons):
        return "—"
    comparison = case.comparisons[index]
    text = comparison.summary()
    if default_reference is None or comparison.run.name != default_reference.name:
        text += f" (vs {comparison.run.name})"
    return text


def wall_time_text(case):
    """Return the previous and current wall times of a case and their change.

    Only the current wall time is given if the previous one is unknown.
    """
    text = format_duration(case.wall_time)
    if case.previous_wall_time is not None:
        text = f"{format_duration(case.previous_wall_time)} → {text}"
    change = wall_time_change(case)
    if change is not None:
        text += f" ({change:+.1%})"
    if case.wall_time_flagged:
        text += f" {WALL_TIME_SYMBOL}"
    return text


def legend_lines(cases, time_threshold):
    """Return the explanation of every symbol used for these cases."""
    used = {case.verdict for case in cases}
    lines = [f"{VERDICTS[verdict].symbol} {VERDICTS[verdict].label}: "
             f"{VERDICTS[verdict].meaning}"
             for verdict in SEVERITY if verdict in used]
    if any(case.wall_time_flagged for case in cases):
        lines.append(f"{WALL_TIME_SYMBOL} wall time changed by more than "
                     f"{100 * time_threshold:g}% since the previous run (only "
                     f"checked for cases that ran for more than "
                     f"{format_duration(MIN_TIMED_SECONDS)})")
    return lines


def case_details(case, diff_path):
    """Return the lines that explain why a case is not unchanged.

    diff_path is the path of the diff of the case, or None if it has none.
    """
    lines = [f"- {problem}" for problem in case.problems]
    for comparison in case.comparisons:
        differences = comparison.differences()
        if not differences:
            continue
        lines.append(f"- vs {comparison.run.name} "
                     f"({comparison.run.short_hash}):")
        lines.extend(f"  - {f.describe()}" for f in differences)
    if diff_path:
        lines.append(f"- full diff: `{DIFF_FOLDER_NAME}/{diff_path.name}`")
    if case.wall_time_flagged:
        lines.append(f"- wall time {format_duration(case.previous_wall_time)}"
                     f" → {format_duration(case.wall_time)}")
    if case.log_excerpt:
        lines += [f"- end of the {SLURM_LOG_NAME} section of this case:", "",
                  "```"] + case.log_excerpt + ["```"]
    return lines


def build_markdown(current, previous_runs, cases, commits, repository,
                   time_threshold, diff_paths):
    previous = previous_runs[0] if previous_runs else None
    before = previous_runs[1] if len(previous_runs) > 1 else None

    lines = [f"# Lethe validation comparison: {current.name}", "",
             "| | Run | Git hash |", "|---|---|---|",
             f"| Current | {current.name} | {current.short_hash} |"]
    if previous:
        lines.append(f"| Previous | {previous.name} | {previous.short_hash} |")
    if before:
        lines.append(f"| Before previous | {before.name} | "
                     f"{before.short_hash} |")
    lines.append("")
    if not previous:
        lines.append("No previous run was found, so there is nothing to "
                     "compare with.")

    for old, new, label in ((previous, current, "previous run"),
                            (before, previous, "run before the previous one")):
        url = compare_url(old, new) if old else None
        if url:
            lines.append(f"Changes since the {label}: {url}  ")
    lines += ["", f"**Summary:** {summary_line(cases)}", "", "## Cases", "",
              "| | Case | vs previous | vs before previous | "
              "Wall time (previous → current) |",
              "|---|---|---|---|---|"]
    for case in cases:
        lines.append(f"| {VERDICTS[case.verdict].symbol} | {case.name} | "
                     f"{reference_cell(case, 0, previous)} | "
                     f"{reference_cell(case, 1, before)} | "
                     f"{wall_time_text(case)} |")
    lines += ["", "**Legend:**", ""]
    lines += [f"- {line}" for line in legend_lines(cases, time_threshold)]

    flagged = [case for case in cases
               if case.verdict != "unchanged" or case.wall_time_flagged]
    if flagged:
        lines += ["", "## Details"]
        for case in sorted(flagged, key=lambda c: SEVERITY.index(c.verdict)):
            verdict = VERDICTS[case.verdict]
            timer = f" {WALL_TIME_SYMBOL}" if case.wall_time_flagged else ""
            lines += ["", f"### {verdict.symbol}{timer} {case.name}: "
                          f"{verdict.description}", ""]
            lines += case_details(case, diff_paths.get(case.name))

    if commits is not None:
        lines += ["", f"## Commits since the previous run ({len(commits)})", ""]
        lines += [f"- {commit}" for commit in commits[:MAX_COMMITS_LISTED]]
        if len(commits) > MAX_COMMITS_LISTED:
            lines.append(f"- ... and {len(commits) - MAX_COMMITS_LISTED} more")
    elif repository and previous:
        lines += ["", f"The commits since the previous run could not be listed "
                      f"from {repository}: one of the hashes is unknown there "
                      f"(was a run made from another branch?)."]
    lines.append("")
    return "\n".join(lines)


def slack_reasons(case):
    """Return the short explanations of a verdict shown on Slack."""
    if case.verdict == "problem":
        reasons = list(case.problems)
        if case.comparisons:
            reasons.extend(f.describe() for f in case.comparisons[0].files
                           if f.status in (STRUCTURE_CHANGED, MISSING))
        return reasons
    if case.verdict in ("changed", "reverted"):
        reasons = [f.describe() for f in case.comparisons[0].changes()]
        if case.verdict == "reverted":
            reasons.append(f"identical to {case.comparisons[1].run.name}")
        return reasons
    if case.verdict == "stabilized":
        return [f"{case.comparisons[1].summary()} vs "
                f"{case.comparisons[1].run.name}"]
    return []


def escape_slack(text):
    """Escape the characters that Slack interprets as markup."""
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def code_block(lines):
    """Return lines as a Slack code block."""
    return "```\n" + "\n".join(escape_slack(line) for line in lines) + "\n```"


def build_slack_message(current, previous_runs, cases, commits,
                        time_threshold):
    previous = previous_runs[0] if previous_runs else None
    title = f"*Lethe validation* `{current.name}` (`{current.short_hash}`)"
    if previous:
        title += f" vs `{previous.name}` (`{previous.short_hash}`)"
    lines = [title, summary_line(cases)]

    for case in sorted(cases, key=lambda c: SEVERITY.index(c.verdict)):
        verdict = VERDICTS[case.verdict]
        reasons = slack_reasons(case)
        if case.wall_time is not None or case.previous_wall_time is not None:
            reasons.append(f"wall time {wall_time_text(case)}")
        lines.append(f"{verdict.symbol} *{case.name}*: {verdict.label}"
                     + (f" — {escape_slack('; '.join(reasons))}"
                        if reasons else ""))

    url = compare_url(previous, current) if previous else None
    if url:
        n_commits = f"{len(commits)} commits" if commits is not None else "diff"
        lines.append(f"<{url}|Changes since the previous run ({n_commits})>")
    lines.append("_Legend:_")
    lines += legend_lines(cases, time_threshold)
    return "\n".join(lines)


def has_slack_details(case):
    """Return whether a case gets a detailed Slack message of its own.

    Only the cases that changed or have a problem do, and only if there is
    more to show than the summary: a diff or the end of the Slurm log.
    """
    return (case.verdict in DETAILED_VERDICTS and
            bool(case.diffs or case.log_excerpt))


def build_slack_detail_message(case, current, diff_path):
    """Return the Slack message that details why a case changed or failed.

    The message states what was compared, lists every changed data file and
    gives the beginning of its diff and, for a case with a problem, the end of
    its section of the Slurm log. Excerpts that would make the message longer
    than MAX_SLACK_MESSAGE_LENGTH are left out, but every file is listed.
    diff_path is the path of the full diff of the case, or None.
    """
    verdict = VERDICTS[case.verdict]
    head = [f"{verdict.symbol} *{case.name}*: {verdict.description}"]
    files = []
    tail = []
    if case.comparisons:
        reference = case.comparisons[0].run
        head.append(f"Data files of `{current.name}` (`{current.short_hash}`)"
                    f" compared with `{reference.name}` "
                    f"(`{reference.short_hash}`).")
        if case.verdict == "reverted":
            before = case.comparisons[1].run
            head.append(f"They are identical to `{before.name}` "
                        f"(`{before.short_hash}`).")
        files = case.comparisons[0].changes()
        url = compare_url(reference, current)
        if url:
            tail.append(f"<{url}|Changes since {reference.name}>")
    head += [f"• {escape_slack(problem)}" for problem in case.problems]
    if diff_path:
        tail.insert(0, f"Full diff: `{diff_path}`")
    omitted_note = "_Excerpts left out to keep this message short._"

    # The text that is always shown. The excerpts are added while they fit.
    length = len("\n".join(head + [f"• {escape_slack(f.describe())}"
                                   for f in files] + tail + [omitted_note]))
    omitted = False

    log_block = []
    if case.log_excerpt:
        title = f"End of the {SLURM_LOG_NAME} section of this case:"
        excerpt = case.log_excerpt
        while excerpt:
            block = [title, code_block(excerpt)]
            if length + len("\n".join(block)) + 1 <= MAX_SLACK_MESSAGE_LENGTH:
                log_block = block
                length += len("\n".join(block)) + 1
                break
            excerpt = excerpt[1:]
            omitted = True

    lines = list(head)
    for comparison in files:
        lines.append(f"• {escape_slack(comparison.describe())}")
        if comparison.name not in case.diffs:
            continue
        excerpt, n_left_out = diff_excerpt(case.diffs[comparison.name],
                                           MAX_SLACK_DIFF_LINES)
        if n_left_out:
            excerpt = excerpt + [f"… {n_left_out} more lines"]
        block = code_block(excerpt)
        if length + len(block) + 1 <= MAX_SLACK_MESSAGE_LENGTH:
            lines.append(block)
            length += len(block) + 1
        else:
            omitted = True
    lines += log_block
    if omitted:
        lines.append(omitted_note)
    lines += tail
    return "\n".join(lines)


def post_to_slack(webhook_url, text):
    """Post a message to a Slack incoming webhook.

    Returns whether the message was posted. Never raises.
    """
    request = urllib.request.Request(
        webhook_url, data=json.dumps({"text": text}).encode("utf-8"),
        headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request,
                                    timeout=SLACK_TIMEOUT_SECONDS) as response:
            response.read()
    except (urllib.error.URLError, OSError, ValueError) as error:
        print(f"Warning: a message could not be posted to Slack: {error}",
              file=sys.stderr)
        return False
    return True


def post_messages_to_slack(webhook_url, messages):
    """Post the summary and then the detailed messages to Slack.

    The detailed messages are not posted if the summary could not be.
    """
    n_posted = 0
    for index, message in enumerate(messages):
        if index > 0:
            time.sleep(SLACK_SECONDS_BETWEEN_MESSAGES)
        if post_to_slack(webhook_url, message):
            n_posted += 1
        elif index == 0:
            break
    print(f"{n_posted} of {len(messages)} messages posted to Slack.")

#############################################################################
# Main
#############################################################################


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Compare a Lethe validation run with the two previous "
                    "runs and report the differences.",
        epilog=f"The diffs of the changed data files are written next to the "
               f"report, in {DIFF_FOLDER_NAME}/. If the environment variable "
               f"{SLACK_WEBHOOK_VARIABLE} contains the URL of a Slack incoming "
               "webhook, a summary is posted to Slack, followed by one message "
               "per case that changed or has a problem.")
    parser.add_argument(
        "-c", "--current", type=Path, required=True,
        help="Folder of the run to check (the output folder of "
             "validate_lethe.sh).")
    parser.add_argument(
        "-a", "--archive", type=Path,
        help="Folder that contains the previous runs. Defaults to the parent "
             "of the current run.")
    parser.add_argument(
        "-p", "--previous", type=Path, nargs="+",
        help="Previous runs to compare with, from the most recent to the "
             "oldest, instead of searching the archive.")
    parser.add_argument(
        "-m", "--match", default="*",
        help="Only consider the archive folders whose name matches this glob "
             "pattern, e.g. '*_Lucille'. Default: '*'.")
    parser.add_argument(
        "-r", "--repo", type=Path,
        help="Lethe git repository, used to list the commits between runs.")
    parser.add_argument(
        "-o", "--report", type=Path,
        help=f"Path of the Markdown report. Defaults to "
             f"<current>/{REPORT_NAME}.")
    parser.add_argument(
        "--rtol", type=float, default=0.0,
        help="Data files whose maximal relative difference is at most this "
             "value are considered unchanged. Default: 0, since validation "
             "runs are bitwise reproducible.")
    parser.add_argument(
        "--time-threshold", type=float, default=0.05,
        help="Relative change of wall time above which a case is flagged. "
             "Default: 0.05.")
    parser.add_argument(
        "--no-slack", action="store_true",
        help="Do not post the summary to Slack.")
    return parser.parse_args()


def main():
    args = parse_arguments()
    # The report contains emojis, which must not make the script fail when the
    # output (e.g. the Slurm log) does not use UTF-8
    sys.stdout.reconfigure(errors="replace")
    if not args.current.is_dir():
        print(f"Error: '{args.current}' is not a folder.", file=sys.stderr)
        return 1

    current = load_run(args.current)
    if args.previous:
        previous_runs = []
        for path in args.previous:
            if not path.is_dir():
                print(f"Error: '{path}' is not a folder.", file=sys.stderr)
                return 1
            previous_runs.append(load_run(path))
    else:
        archive = args.archive or args.current.resolve().parent
        previous_runs = find_previous_runs(current, archive, args.match)

    cases = compare_cases(current, previous_runs, args.rtol,
                          args.time_threshold)
    commits = (list_commits(args.repo, previous_runs[0], current)
               if previous_runs else None)

    report_path = (args.report or current.path / REPORT_NAME).resolve()
    try:
        diff_paths = write_diffs(cases, report_path.parent / DIFF_FOLDER_NAME)
        report = build_markdown(current, previous_runs, cases, commits,
                                args.repo, args.time_threshold, diff_paths)
        report_path.write_text(report, encoding="utf-8")
    except OSError as error:
        print(f"Error: the report could not be written: {error}",
              file=sys.stderr)
        return 1
    print(report)
    print(f"Report written to {report_path}")

    webhook_url = os.environ.get(SLACK_WEBHOOK_VARIABLE)
    if args.no_slack:
        pass
    elif webhook_url:
        messages = [build_slack_message(current, previous_runs, cases,
                                        commits, args.time_threshold)]
        messages += [build_slack_detail_message(case, current,
                                                diff_paths.get(case.name))
                     for case in sorted(cases,
                                        key=lambda c: SEVERITY.index(c.verdict))
                     if has_slack_details(case)]
        post_messages_to_slack(webhook_url, messages)
    else:
        print(f"{SLACK_WEBHOOK_VARIABLE} is not set: nothing posted to Slack.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
