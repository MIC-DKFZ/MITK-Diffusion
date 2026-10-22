#!/usr/bin/env python3
"""Search or create issues in mic/personal/group6/neuro/mitk-diffusion on git.dkfz.de for the issue
skills (bug-report, feature-request, improvement-request), which create an issue only when the user
asks for it explicitly. The token comes from GITLAB_TOKEN in the environment or the mitk-diffusion
repo's .env and is never printed.

  create_gitlab_issue.py search <words...>
  create_gitlab_issue.py labels
  create_gitlab_issue.py create gitlab_issues/new-<branch>.md [--dry-run]
"""
import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://git.dkfz.de/api/v4"
PROJECT = "mic/personal/group6/neuro/mitk-diffusion"
REPO_ROOT = Path(__file__).resolve().parents[2]
LABEL_LINE = re.compile(r"^/label\s+(.+)$")
LABEL = re.compile(r'~"([^"]+)"|~(\S+)')


def token():
    value = os.environ.get("GITLAB_TOKEN", "")
    env = REPO_ROOT / ".env"
    if not value and env.is_file():
        # Last line wins, as when a shell sources the file.
        for line in env.read_text().splitlines():
            if line.startswith("GITLAB_TOKEN="):
                value = line.split("=", 1)[1].strip().strip("\"'")
    if not value:
        sys.exit("GITLAB_TOKEN is not set in the environment or in .env")
    return value


def call(method, path, data=None):
    url = f"{API}/projects/{urllib.parse.quote(PROJECT, safe='')}{path}"
    body = urllib.parse.urlencode(data).encode() if data else None
    request = urllib.request.Request(url, data=body, method=method,
                                     headers={"PRIVATE-TOKEN": token()})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return json.load(response)
    except urllib.error.HTTPError as err:
        sys.exit(f"GitLab {method} {path.split('?')[0]}: HTTP {err.code} "
                 f"{err.read().decode(errors='replace')[:300]}")


def parse(path):
    """Title from the '# ' first line; /label quick-action lines become the label list."""
    lines = path.read_text().splitlines()
    if not lines or not lines[0].startswith("# "):
        sys.exit(f"{path}: the first line must be the title as '# <title>'")
    labels, body = [], []
    for line in lines[1:]:
        match = LABEL_LINE.match(line.strip())
        if match:
            labels += [quoted or bare for quoted, bare in LABEL.findall(match.group(1))]
        else:
            body.append(line)
    return lines[0][2:].strip(), labels, "\n".join(body).strip() + "\n"


def search(args):
    query = urllib.parse.urlencode({"search": " ".join(args.words), "per_page": 20})
    for issue in call("GET", f"/issues?{query}"):
        print(f"#{issue['iid']} {issue['state']} {issue['title']}  {issue['web_url']}")


def known_labels():
    return {label["name"] for label in
            call("GET", "/labels?per_page=100&include_ancestor_groups=true")}


def labels(args):
    for name in sorted(known_labels()):
        print(name)


def create(args):
    path = Path(args.file)
    # Only a file still waiting for its number may be created; anything else names an issue
    # that already exists, and creating it again would be a duplicate.
    if not path.name.startswith("new-"):
        sys.exit(f"{path}: not a new-<branch>.md file, so it names an existing issue; not creating")
    title, labels, description = parse(path)
    known = known_labels()
    missing = [label for label in labels if label not in known]
    if missing:
        # The API creates a label it does not know, so a typo would become a new project label.
        sys.exit(f"labels not in {PROJECT}: {missing}; known: {sorted(known)}")
    query = urllib.parse.urlencode({"search": title, "in": "title", "state": "opened"})
    same = [issue for issue in call("GET", f"/issues?{query}") if issue["title"] == title]
    if same:
        sys.exit(f"an open issue has this title: #{same[0]['iid']} {same[0]['web_url']}")
    if args.dry_run:
        print(f"dry run, nothing created: {PROJECT}: {title!r} labels={labels} "
              f"description={len(description)} chars")
        return
    issue = call("POST", "/issues",
                 {"title": title, "description": description, "labels": ",".join(labels)})
    target = path.with_name(f"{issue['iid']}-{path.name[len('new-'):]}")
    path.rename(target)
    print(f"created #{issue['iid']} {issue['web_url']} labels={issue['labels']} file={target}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    find = commands.add_parser("search", help="list issues whose title or description match")
    find.add_argument("words", nargs="+")
    find.set_defaults(run=search)
    commands.add_parser("labels", help="list the labels the project has").set_defaults(run=labels)
    new = commands.add_parser("create", help="create the issue held in a gitlab_issues/ file")
    new.add_argument("file")
    new.add_argument("--dry-run", action="store_true", help="check everything, create nothing")
    new.set_defaults(run=create)
    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
