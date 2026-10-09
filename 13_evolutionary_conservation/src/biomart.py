"""Minimal Ensembl BioMart client (GET only; the server rejects POST).

BioMart rejects queries whose attributes span more than one attribute page, so
callers fetch one attribute group per query and merge locally on stable IDs.
"""
import io
import time

import pandas as pd
import requests

URL = "https://www.ensembl.org/biomart/martservice"


def query(dataset, attrs, filters=None, timeout=900, retries=5):
    flt = "".join(
        f'<Filter name="{k}" value="{",".join(v) if isinstance(v, (list, tuple)) else v}"/>'
        for k, v in (filters or {}).items()
    )
    att = "".join(f'<Attribute name="{a}"/>' for a in attrs)
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
        '<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" '
        'count="" datasetConfigVersion="0.6">'
        f'<Dataset name="{dataset}" interface="default">{flt}{att}</Dataset></Query>'
    )
    last = None
    for attempt in range(retries):
        try:
            r = requests.get(URL, params={"query": xml}, timeout=timeout)
            r.raise_for_status()
            if r.text.lstrip().startswith("Query ERROR"):
                raise RuntimeError(r.text[:300])
            if r.text.lstrip()[:15].lower().startswith(("<html", "<!doctype")):
                # Ensembl returns a 200 'Service unavailable' page under load
                raise ConnectionError("BioMart returned an HTML error page")
            return pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)
        except RuntimeError:
            raise  # malformed query: retrying will not help
        except Exception as e:  # network / 5xx
            last = e
            time.sleep(30 * (attempt + 1))
    raise last


def release():
    r = requests.get(URL, params={"type": "registry"}, timeout=60)
    import re
    m = re.search(r'database="ensembl_mart_(\d+)"', r.text)
    return m.group(1) if m else "unknown"


def attribute_names(dataset):
    """Set of valid attribute internal names for a BioMart dataset."""
    r = requests.get(URL, params={"type": "attributes", "dataset": dataset}, timeout=120)
    r.raise_for_status()
    return {line.split("\t")[0] for line in r.text.splitlines() if line}
