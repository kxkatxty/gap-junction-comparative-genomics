# Standalone CLIs

Python scripts that are not `pipeline.*` modules. Run from the repo root:

```bash
python3 tools/discover_innexins.py --family innexin --help
python3 tools/extract_exon_structures.py project/metadata/species_config.csv
python3 tools/plot_exon_maps.py
```

They import each other (e.g. discovery uses the GFF downloader), so keep this folder on `PYTHONPATH` (pytest already does).
