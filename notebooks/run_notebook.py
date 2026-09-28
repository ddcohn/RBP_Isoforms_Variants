import sys

import nbformat
from nbclient import NotebookClient

# usage: run_notebook.py <notebook.ipynb> [...]
# Executes each notebook in place. Run from notebooks/ so relative data
# paths inside the notebooks (e.g. "clinvar_spliceai_scores.tsv") resolve.
for path in sys.argv[1:]:
    nb = nbformat.read(path, as_version=4)
    NotebookClient(nb, timeout=1800, kernel_name="python3").execute()
    nbformat.write(nb, path)
    print(f"executed {path}")
