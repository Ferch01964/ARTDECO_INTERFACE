# main_kdis_disort_artdeco.py
"""
Función que genera los ficheros main.py / kdis.py / disort.py para ARTDECO
en los rangos SW-LW.

Autor original (MATLAB): Cristina Gil Díaz
Traducción a Python: 2025-06-24

Parámetros
----------
base_path : str | pathlib.Path
    Ruta al directorio raíz que contiene la carpeta INPUTS.
min_wvl, max_wvl : float
    Longitudes de onda mínima y máxima (micras).
date : datetime.datetime
    Fecha/hora que formará parte del nombre de la campaña.
bool_cld, bool_aer : bool
    Indican si la escena atmosférica contiene nubes y/o aerosoles.
bool_sw, bool_lw : bool
    Indican si se generan ficheros para simulación en SW y/o LW.

Esta función no devuelve nada; crea directorios y escribe ficheros.
"""

from __future__ import annotations
import shutil
from datetime import datetime
from pathlib import Path


def _read_template(path: Path) -> list[str]:
    """Lee un fichero de plantilla y devuelve la lista de líneas."""
    return path.read_text(encoding="utf-8").splitlines()


def _write_file(path: Path, lines: list[str]) -> None:
    """Escribe todas las líneas en path (añade \n)."""
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _replace_line(
    lines: list[str], starts_with: str, new_line: str, fallback_idx: int | None = None
) -> None:
    """
    Sustituye la línea que comienza con `starts_with` por `new_line`.
    Si no se encuentra, usa fallback_idx (mismo comportamiento que el MATLAB original
    donde se empleaban índices absolutos).
    """
    for i, line in enumerate(lines):
        if line.strip().startswith(starts_with):
            lines[i] = new_line
            return
    if fallback_idx is not None and 0 <= fallback_idx < len(lines):
        lines[fallback_idx] = new_line


def main_kdis_disort_artdeco(
    base_path: str | Path,
    min_wvl: float,
    max_wvl: float,
    date: datetime,
    bool_cld: bool,
    bool_aer: bool,
    bool_sw: bool,
    bool_lw: bool,
) -> None:
    base_path = Path(base_path)

    # 1. Tipo de escena atmosférica
    if bool_cld and bool_aer:
        str_type = "GAC"
    elif bool_cld and not bool_aer:
        str_type = "CLD"
    elif not bool_cld and bool_aer:
        str_type = "AER"
    else:
        str_type = "GAS"

    # 2. Plantilla para main.py
    main_template_path = base_path  / "templates" / "main.py"
    main_template_lines = _read_template(main_template_path)

    # 3. Bucle para SW y LW (evitamos repetir código)
    for band, band_flag in (("SW", bool_sw), ("LW", bool_lw)):
        if not band_flag:
            continue

        # 3.1. Directorio de destino (…/str_type/SW|LW/SW_yyyymmddhh/config)
        campaign = f"{band}_{date:%Y%m%d%H}"
        dest_dir = (
            base_path / str_type / band / campaign / "config"
        )
        dest_dir.mkdir(parents=True, exist_ok=True)

        # 3.2. Copiar disort.py
        shutil.copy2(
            base_path / "templates" / "disort.py",
            dest_dir / "disort.py",
        )

        # 3.3. Personalizar main.py
        lines = main_template_lines.copy()  # partimos siempre de la plantilla original

        _replace_line(
            lines,
            "With_particles",
            f'With_particles = {bool(bool_cld or bool_aer)}',
            fallback_idx=6,  # línea 7 en MATLAB (0-based aquí)
        )
        _replace_line(
            lines,
            "Prefix_of_input_files",
            f'Prefix_of_input_files = "{str_type}/{band}/{campaign}"',
            fallback_idx=38,  # línea 39 en MATLAB
        )
        _replace_line(
            lines,
            "Output_directory",
            f'Output_directory = "{str_type}/{band}/{campaign}"',
            fallback_idx=50,  # línea 51 en MATLAB
        )

        _write_file(dest_dir / "main.py", lines)

    # 4. Plantilla y personalización de kdis.py
    #    Se escribe en el último dest_dir calculado (igual que en el MATLAB original).
    #    Si sólo se ha creado SW, ése será el destino; si sólo LW, idem.
    #    Si se crean ambos, el fichero termina en LW_…/config.
    kdis_template_path = base_path / "templates" / "kdis.py"
    kdis_lines = _read_template(kdis_template_path)

    k_dis_model = '"gamesw_1"' if bool_sw else '"gamelw_1"'
    _replace_line(
        kdis_lines,
        "K_distribution_parameterization",
        f"K_distribution_parameterization = {k_dis_model}",
        fallback_idx=9,  # línea 10 MATLAB
    )
    _replace_line(
        kdis_lines,
        "Wavelength_start",
        f"Wavelength_start = {min_wvl:.3f}",
        fallback_idx=29,
    )
    _replace_line(
        kdis_lines,
        "Wavelength_end",
        f"Wavelength_end = {max_wvl:.3f}",
        fallback_idx=30,
    )

    _write_file(dest_dir / "kdis.py", kdis_lines)