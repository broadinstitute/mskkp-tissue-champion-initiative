#!/usr/bin/env python3
"""
Regenerate docs/mskkp_schema_docs.html's embedded data block from the YAML
schema files in schemas/.

The HTML page renders entirely client-side from a JSON blob embedded in a
<script type="application/json" id="embedded-schema-data"> tag, so it can be
opened directly (file://, double-click) without a local server. This script
is what keeps that blob in sync with the YAML sources — run it after editing
any of mskkp_rnaseq_schema.yaml, tissue_profiles.yaml, champion_criteria.yaml,
or champions_roster.yaml.

Usage:
    python3 docs/build_docs.py
"""

import json
import re
import sys
from pathlib import Path

import yaml

DOCS_DIR = Path(__file__).resolve().parent
REPO_ROOT = DOCS_DIR.parent
SCHEMA_DIR = REPO_ROOT / "schemas"
HTML_PATH = DOCS_DIR / "mskkp_schema_docs.html"

EMBED_RE = re.compile(
    r'(<script type="application/json" id="embedded-schema-data">)(.*?)(</script>)',
    re.DOTALL,
)


def load_yaml(name):
    """Load a YAML file from schemas/ by filename.

    Args:
        name (str): Filename within SCHEMA_DIR, e.g. "mskkp_rnaseq_schema.yaml".

    Returns:
        dict: The parsed YAML document.
    """
    with open(SCHEMA_DIR / name) as f:
        return yaml.safe_load(f)


def build_embedded_data():
    """Load all four source YAML files into the shape init() expects.

    Returns:
        dict: Keys ``schema``, ``tissueData``, ``criteriaData``,
        ``rosterData`` — matching the destructuring in the page's
        ``init()`` function.
    """
    return {
        "schema": load_yaml("mskkp_rnaseq_schema.yaml"),
        "tissueData": load_yaml("tissue_profiles.yaml"),
        "criteriaData": load_yaml("champion_criteria.yaml"),
        "rosterData": load_yaml("champions_roster.yaml"),
    }


def to_embeddable_json(data):
    """Serialize data for safe inclusion inside an HTML <script> tag.

    Args:
        data (dict): The object to serialize.

    Returns:
        str: Compact JSON with "</" escaped to "<\\/" so a literal
        "</script>" inside any string value can't prematurely close the
        surrounding <script> tag.
    """
    return json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def main():
    """Rebuild the embedded JSON block in docs/mskkp_schema_docs.html in place.

    Returns:
        None.

    Raises:
        SystemExit: If the HTML file is missing the
        ``#embedded-schema-data`` script tag to replace.
    """
    html = HTML_PATH.read_text()
    if not EMBED_RE.search(html):
        sys.exit(f"Couldn't find #embedded-schema-data script tag in {HTML_PATH}")

    payload = to_embeddable_json(build_embedded_data())
    new_html = EMBED_RE.sub(lambda m: m.group(1) + payload + m.group(3), html, count=1)
    HTML_PATH.write_text(new_html)
    print(f"Regenerated {HTML_PATH.relative_to(REPO_ROOT)} ({len(payload):,} bytes of embedded data).")


if __name__ == "__main__":
    main()
