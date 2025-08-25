# Dependency Inventory

Erzeugt: 2025-08-24T14:32:39.226539+00:00  
Quelle JSON: `tools/inventory/dependency_graph.json`

## Pakete

### electrofit (electrofit)

Version: 0.1.0

Dependencies (deklariert):

- app
- typer
- tqdm
- tomli_w
- tomli>=2.0.1; python_version < "3.11"
- openmol
Scripts:

- electrofit -> electrofit.cli.app:main
- electrofit-plot-charges -> electrofit.cli.plot_charges:main
- electrofit-write-symmetry -> electrofit.cli.write_symmetry:main

### electrofit-analysis (electrofit-analysis)

Version: 0.1.0

Scripts:

- electrofit-analysis -> electrofit_analysis.cli:main

## Paket-Kanten (interne Cross-Package Imports)

- electrofit_analysis -> electrofit: 14

## Externe Top-Level Imports (Häufigkeit)

- os: 47
- logging: 45
- __future__: 45
- pathlib: 36
- typing: 27
- argparse: 25
- matplotlib: 21
- numpy: 19
- sys: 14
- json: 11
- shutil: 9
- re: 9
- seaborn: 9
- MDAnalysis: 8
- pandas: 7
- dataclasses: 7
- rdkit: 6
- subprocess: 5
- time: 5
- contextlib: 4
- fnmatch: 4
- traceback: 3
- tomli: 3
- tomllib: 3
- mdtraj: 3
- networkx: 3
- atexit: 2
- signal: 2
- threading: 2
- copy: 2
- filecmp: 2
- hashlib: 2
- glob: 2
- uuid: 2
- concurrent: 2
- multiprocessing: 2
- math: 2
- itertools: 2
- scipy: 2
- types: 2
- mpl_toolkits: 2
- psutil: 1
- runpy: 1
- io: 1
- tomli_w: 1
- faulthandler: 1
- gc: 1
- random: 1
- warnings: 1
- openbabel: 1
- openmol: 1
- errno: 1
- tqdm: 1
- mmap: 1

## Modulzählung

Anzahl Module: 94

## Hinweise

- Dynamische / optionale Imports sind nicht speziell markiert.
- Nutzung (Symbole) nicht aufgelöst; es handelt sich nur um Modulreferenzen.
- Bei Refactor erneut Skript ausführen und Änderungen committen.