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
the standard output. If the environment variable
LETHE_VALIDATION_SLACK_WEBHOOK contains the URL of a Slack incoming webhook,
a summary of every case, with a legend of the symbols it uses, is also posted
to Slack.

Usage:
  python3 compare_validation_runs.py --current <run folder> [options]

Run with --help for the list of options.
"""

import argparse
import fnmatch
import json
import math
import os
import re
import subprocess
import sys
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

    def describe(self):
        if self.status == CHANGED:
            return (f"{self.name}: changed, max. relative difference "
                    f"{self.max_relative_difference:.1e} (line {self.line})")
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
                return (f"changed (max. rel. diff. "
                        f"{self.max_relative_difference:.1e})")
            return "changed (new files)"
        return ", ".join(sorted({f.status for f in self.files
                                 if f.status in (STRUCTURE_CHANGED, MISSING)}))

    def differences(self):
        return [f for f in self.files if f.status != IDENTICAL]


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
            files.append(FileComparison(name, NEW))
        elif name not in current_files:
            files.append(FileComparison(name, MISSING))
        else:
            comparison = compare_data_file(name, current_files[name],
                                           reference_files[name])
            if (comparison.status == SAME_VALUES and
                    comparison.max_relative_difference > rtol):
                comparison.status = CHANGED
            files.append(comparison)
    return ReferenceComparison(run=reference, files=files)

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

    results = []
    for name in sorted(names):
        references = [run for run in previous_runs if name in run.cases][:2]
        if name not in current.cases:
            results.append(CaseResult(name=name,
                                      problems=["case missing from this run"],
                                      comparisons=[], verdict="problem"))
            continue

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
            max(case.wall_time, case.previous_wall_time) > MIN_TIMED_SECONDS)
        case.verdict = decide_verdict(case)
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


def case_details(case):
    """Return the lines that explain why a case is not unchanged."""
    lines = [f"- {problem}" for problem in case.problems]
    for comparison in case.comparisons:
        differences = comparison.differences()
        if not differences:
            continue
        lines.append(f"- vs {comparison.run.name} "
                     f"({comparison.run.short_hash}):")
        lines.extend(f"  - {f.describe()}" for f in differences)
    if case.wall_time_flagged:
        lines.append(f"- wall time {format_duration(case.previous_wall_time)}"
                     f" → {format_duration(case.wall_time)}")
    return lines


def build_markdown(current, previous_runs, cases, commits, repository,
                   time_threshold):
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
            lines += case_details(case)

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
        reasons = [difference_summary(case.comparisons[0])]
        if case.verdict == "reverted":
            reasons.append(f"identical to {case.comparisons[1].run.name}")
        return reasons
    if case.verdict == "stabilized":
        return [f"{case.comparisons[1].summary()} vs "
                f"{case.comparisons[1].run.name}"]
    return []


def difference_summary(comparison):
    if comparison.max_relative_difference > 0.0:
        return f"max. rel. diff. {comparison.max_relative_difference:.1e}"
    return "new data files"


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
                     + (f" — {'; '.join(reasons)}" if reasons else ""))

    url = compare_url(previous, current) if previous else None
    if url:
        n_commits = f"{len(commits)} commits" if commits is not None else "diff"
        lines.append(f"<{url}|Changes since the previous run ({n_commits})>")
    lines.append("_Legend:_")
    lines += legend_lines(cases, time_threshold)
    return "\n".join(lines)


def post_to_slack(webhook_url, text):
    """Post a message to a Slack incoming webhook. Never raises."""
    request = urllib.request.Request(
        webhook_url, data=json.dumps({"text": text}).encode("utf-8"),
        headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request,
                                    timeout=SLACK_TIMEOUT_SECONDS) as response:
            response.read()
    except (urllib.error.URLError, OSError, ValueError) as error:
        print(f"Warning: the summary could not be posted to Slack: {error}",
              file=sys.stderr)
        return
    print("The summary was posted to Slack.")

#############################################################################
# Main
#############################################################################


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Compare a Lethe validation run with the two previous "
                    "runs and report the differences.",
        epilog=f"If the environment variable {SLACK_WEBHOOK_VARIABLE} "
               "contains the URL of a Slack incoming webhook, a summary is "
               "posted to Slack.")
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

    report = build_markdown(current, previous_runs, cases, commits,
                            args.repo, args.time_threshold)
    report_path = args.report or current.path / REPORT_NAME
    try:
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
        post_to_slack(webhook_url,
                      build_slack_message(current, previous_runs, cases,
                                          commits, args.time_threshold))
    else:
        print(f"{SLACK_WEBHOOK_VARIABLE} is not set: nothing posted to Slack.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
