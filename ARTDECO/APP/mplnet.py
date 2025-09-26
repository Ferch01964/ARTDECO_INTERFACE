# mplnet_downloader.py  (versión revisada)
from __future__ import annotations
import zipfile
from datetime import datetime, timedelta, date
from pathlib import Path
from typing import Union, List

import requests


BASE_URL = "https://mplnet.gsfc.nasa.gov/download"


def download_mplnet_range(
    target_date: Union[str, date],
    *,
    site: str = "Barcelona",
    product: str = "AER",
    var: str = "aerosol_top",
    version: str = "V3",
    level: str = "L1",
    download_dir: Union[str, Path] = Path("INPUTS/resources"),
    unzip: bool = True,
    keep_zip: bool = True,
) -> Union[Path, List[Path]]:
    """
    Descarga el ZIP MPLNET para `target_date` y el día anterior.
    Si `unzip=True` devuelve *una lista* con las rutas completas de los
    archivos extraídos.  Si `unzip=False` devuelve la ruta del ZIP.

    Parameters
    ----------
    target_date : str | date | datetime
        Fecha objetivo (YYYYMMDD).
    download_dir : str | Path
        Carpeta base para guardar el ZIP y los archivos extraídos.
    keep_zip : bool
        Si `False`, elimina el ZIP tras la extracción.

    Returns
    -------
    Path | list[Path]
        • `unzip=True`  →  lista de Path (archivos dentro del ZIP).  
        • `unzip=False` →  Path al ZIP descargado.
    """
    # ------------------------------------------------------------------
    # 1. Validación / normalización de la fecha
    # ------------------------------------------------------------------
    if isinstance(target_date, str):
        if len(target_date) != 8 or not target_date.isdigit():
            raise ValueError("La fecha debe tener formato 'YYYYMMDD'.")
        target_date = datetime.strptime(target_date, "%Y%m%d").date()
    elif isinstance(target_date, datetime):
        target_date = target_date.date()
    elif not isinstance(target_date, date):
        raise TypeError("target_date debe ser str, datetime o date.")

    start_date = target_date - timedelta(days=1)

    # ------------------------------------------------------------------
    # 2. Parámetros de la query
    # ------------------------------------------------------------------
    params = {
        "app": "download_tool",
        "site": site,
        "year":  start_date.year,
        "month": start_date.month,
        "day":   start_date.day,
        "eyear":  target_date.year,
        "emonth": target_date.month,
        "eday":   target_date.day,
        "version": version,
        "level": level,
        "product": product,
        "no_qa": "false",
        "no_flags": "false",
        "var": var,
    }

    # ------------------------------------------------------------------
    # 3. Descarga del ZIP
    # ------------------------------------------------------------------
    download_dir = Path(download_dir)
    download_dir.mkdir(parents=True, exist_ok=True)

    print(f"↓  Solicitando datos MPLNET ({start_date} → {target_date})…")
    response = requests.get(BASE_URL, params=params, stream=True, timeout=120)
    response.raise_for_status()

    # Nombre del ZIP (header o fallback)
    cd = response.headers.get("Content-Disposition", "")
    filename = (
        cd.split("filename=")[-1].strip('"')
        if "filename=" in cd
        else f"mplnet_{site}_{start_date:%Y%m%d}-{target_date:%Y%m%d}.zip"
    )
    zip_path = download_dir / filename

    # Guardar a disco
    with open(zip_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"✓ ZIP descargado: {zip_path}")

    # ------------------------------------------------------------------
    # 4. Descompresión (si procede) y retorno de *paths*
    # ------------------------------------------------------------------
    if not unzip:
        return zip_path

    extract_dir = zip_path.with_suffix("")          # Ej.: xxx.zip → xxx
    print(f"⇢  Extrayendo en {extract_dir} …")

    with zipfile.ZipFile(zip_path) as zf:
        members = [
            m for m in zf.namelist() if not m.endswith("/")  # sólo archivos
        ]
        zf.extractall(extract_dir)

    print("✓ Extracción completada.")

    if not keep_zip:
        zip_path.unlink(missing_ok=True)

    # Construir lista de rutas completas de los archivos extraídos
    extracted_files = [extract_dir / m for m in members]
    return extracted_files