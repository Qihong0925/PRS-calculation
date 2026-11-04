#!/usr/bin/env python3
"""Look up GRCh37 CanRisk sites in Ensembl: rsID at the GRCh37 position and its GRCh38 mapping."""
import json
import sys
import urllib.request

SITES = [  # chr, GRCh37 pos, ref, effect, our GRCh38 pos
    ("10", 38523626, "C", "A", 38234698),
    ("10", 71335574, "C", "T", 69575818),
    ("17", 29168077, "G", "T", 30841059),
    ("6",  87803819, "T", "C", 87094101),
    ("5",  58241712, "C", "T", 58945885),
    ("4", 106069013, "G", "T", 105147856),
    ("22", 45319953, "G", "A", 44924073),
]


def get(url):
    sep = "&" if "?" in url else "?"
    with urllib.request.urlopen(f"{url}{sep}content-type=application/json", timeout=30) as r:
        return json.load(r)


for chrom, pos, ref, eff, our38 in SITES:
    try:
        v37 = get(f"https://grch37.rest.ensembl.org/overlap/region/human/{chrom}:{pos}-{pos}?feature=variation")
    except Exception as e:
        print(f"{chrom}:{pos}\tERROR {e}", file=sys.stderr)
        continue
    hits = [v for v in v37 if v["start"] == pos and v["id"].startswith("rs")]
    if not hits:
        print(f"{chrom}:{pos} {ref}>{eff}\tno rsID at GRCh37 position")
        continue
    for v in hits:
        d = get(f"https://rest.ensembl.org/variation/human/{v['id']}")
        maps = [f"{m['seq_region_name']}:{m['start']}({m['allele_string']})"
                for m in d.get("mappings", []) if m.get("assembly_name") == "GRCh38"]
        print(f"{chrom}:{pos} {ref}>{eff}\t{v['id']}\tGRCh37 alleles={v.get('alleles')}\t"
              f"GRCh38={';'.join(maps)}\tours=chr{chrom}:{our38}\tMAF={d.get('MAF')} minor={d.get('minor_allele')}")
