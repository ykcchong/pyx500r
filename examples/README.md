# pyx500r examples

Runnable applications built on the `pyx500r` library.

| Example | What it is |
|---------|------------|
| [`server/`](server/) | Generic **FastAPI/Uvicorn** JSON API over `.wiff2` acquisitions and `.qsession` results — the backend reference for a custom TypeScript GUI. See [`../docs/GUI_INTEGRATION.md`](../docs/GUI_INTEGRATION.md). |
| [`library_browser/`](library_browser/) | A self-contained **GUI app** to browse a SCIEX **LibraryView** `.sqlite` database: search compounds, view their settings, and render their (often multiple) reference spectra as interactive plots. |

Both require the FastAPI extras on top of the base install:

```bash
pip install -e ".[numba]"
pip install -r examples/<example>/requirements.txt
```

Then follow the per-example README to launch `uvicorn`.
