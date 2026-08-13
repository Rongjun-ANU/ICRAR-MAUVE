#!/usr/bin/env python3
"""Resolve the review DOI seed list against Crossref and write an auditable registry."""

from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ASSET_DIR = Path(__file__).resolve().parent
SEED_PATH = ASSET_DIR / "doi_seeds.tsv"
USER_AGENT = (
    "MAUVE-analytic-environment-review/20260813 "
    "(bibliographic metadata verification; contact via institutional project)"
)


def first(value, default=""):
    if isinstance(value, list):
        return value[0] if value else default
    return value if value is not None else default


def date_year(message):
    for field in ("published-print", "published-online", "published", "issued"):
        parts = message.get(field, {}).get("date-parts", [])
        if parts and parts[0]:
            return parts[0][0]
    return ""


def page_value(message):
    return message.get("page") or message.get("article-number") or ""


def author_value(message):
    names = []
    for author in message.get("author", []):
        literal = author.get("literal")
        if literal:
            names.append(literal)
            continue
        family = author.get("family", "")
        given = author.get("given", "")
        names.append(", ".join(part for part in (family, given) if part))
    return "; ".join(names)


def resolve(seed):
    doi = seed["doi"].strip().lower()
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error = ""
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                payload = json.load(response)
            message = payload["message"]
            return {
                **seed,
                "status": "resolved",
                "http_status": 200,
                "title": first(message.get("title")),
                "authors": author_value(message),
                "year": date_year(message),
                "journal": first(message.get("container-title")),
                "volume": message.get("volume", ""),
                "issue": message.get("issue", ""),
                "page": page_value(message),
                "type": message.get("type", ""),
                "publisher": message.get("publisher", ""),
                "crossref_url": message.get("URL", "https://doi.org/" + doi),
                "doi_url": "https://doi.org/" + doi,
                "error": "",
            }
        except urllib.error.HTTPError as exc:
            last_error = f"HTTP {exc.code}"
            if exc.code == 404:
                break
        except Exception as exc:  # network and malformed-response audit trail
            last_error = f"{type(exc).__name__}: {exc}"
        time.sleep(0.75 * (attempt + 1))
    return {
        **seed,
        "status": "failed",
        "http_status": "",
        "title": "",
        "authors": "",
        "year": "",
        "journal": "",
        "volume": "",
        "issue": "",
        "page": "",
        "type": "",
        "publisher": "",
        "crossref_url": "",
        "doi_url": "https://doi.org/" + doi,
        "error": last_error,
    }


def main():
    with SEED_PATH.open(newline="", encoding="utf-8") as handle:
        seeds = list(csv.DictReader(handle, delimiter="\t"))

    records = []
    with ThreadPoolExecutor(max_workers=5) as executor:
        future_map = {executor.submit(resolve, seed): seed for seed in seeds}
        for future in as_completed(future_map):
            records.append(future.result())
    order = {seed["key"]: index for index, seed in enumerate(seeds)}
    records.sort(key=lambda record: order[record["key"]])

    payload = {
        "checked_on": "2026-08-13",
        "service": "Crossref REST API",
        "n_seeded_dois": len(records),
        "n_resolved": sum(record["status"] == "resolved" for record in records),
        "n_failed": sum(record["status"] != "resolved" for record in records),
        "checks": records,
    }
    (ASSET_DIR / "doi_resolution_audit.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    fields = [
        "key", "family", "year", "title", "authors", "journal", "volume",
        "issue", "page", "doi", "doi_url", "type", "publisher",
        "inclusion_role", "status", "error",
    ]
    with (ASSET_DIR / "source_registry.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)

    print(json.dumps({key: payload[key] for key in ("n_seeded_dois", "n_resolved", "n_failed")}))
    for record in records:
        if record["status"] != "resolved":
            print(f"FAILED\t{record['key']}\t{record['doi']}\t{record['error']}")


if __name__ == "__main__":
    main()
