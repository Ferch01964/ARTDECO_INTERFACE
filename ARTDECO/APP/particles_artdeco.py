# particles_artdeco.py  –– versión corregida
from __future__ import annotations

import os
from pathlib import Path
from datetime import datetime
from typing import Sequence, Tuple, Union, List

import numpy as np
import pandas as pd


# Tipos “equivalentes” a los cell‑array de MATLAB
GeoProps            = np.ndarray                # [cbh, cth, phase, delta]
CldMatrix           = np.ndarray                # m_ext | m_asy | m_ssa  (heights×λ)
AeronetProf         = np.ndarray                # (Nz, 3)  [v_weight, v_aod_550, v_ext_550]
AeronetInterp       = np.ndarray                # (4, M)   λ | AOD | SSA | ASY
CloudProperties     = Tuple[GeoProps, CldMatrix, CldMatrix, CldMatrix]
AerosolProperties   = Tuple[AeronetProf, AeronetInterp]


# ──────────────────────────────────────────────────────────────────────────────
# UTILIDADES
# ──────────────────────────────────────────────────────────────────────────────
def _read_column(path: Path) -> np.ndarray:
    """Lee un fichero de una sola columna de floats (espaciado libre)."""
    # ‣ FIX 1 : evitar el TypeError → usar keyword dtype o, más simple, np.loadtxt
    return np.loadtxt(path, dtype=float)


def _read_template(path: Path) -> List[str]:
    """Lee un archivo de plantilla y devuelve sus líneas (incluye saltos NL)."""
    return path.read_text(encoding="utf-8").splitlines()


def _write_text(path: Path, lines: Sequence[str]) -> None:
    """Escribe un archivo línea a línea añadiendo '\n'."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


# ──────────────────────────────────────────────────────────────────────────────
# FUNCIÓN PRINCIPAL
# ──────────────────────────────────────────────────────────────────────────────
def particles_artdeco(
    base_path: str | Path,
    date: datetime,
    bool_cld: bool,
    bool_aer: bool,
    bool_sw: bool,
    bool_lw: bool,
    prop_cld: Union[CloudProperties, None],
    prop_aer: Union[AerosolProperties, None],
) -> None:
    """
    Crea los ficheros *particles.py* y *opt_*.dat para ARTDECO (SW / LW).
    """
    base_path = Path(base_path)

    # ------------------------------------------------------------------ #
    # 1. Tipo de escena (GAC / CLD / AER / GAS)
    # ------------------------------------------------------------------ #
    if bool_cld and bool_aer:
        str_type = "GAC"
    elif bool_cld and not bool_aer:
        str_type = "CLD"
    elif not bool_cld and bool_aer:
        str_type = "AER"
    else:
        str_type = "GAS"

    # ------------------------------------------------------------------ #
    # 2. Leer ficheros de niveles y longitudes de onda
    # ------------------------------------------------------------------ #
    z_artdeco_sw = _read_column(base_path / "resources/z_p_ARTDECO_SW.txt")
    z_artdeco_lw = _read_column(base_path / "resources/z_p_ARTDECO_LW.txt")
    λ_sw         = _read_column(base_path / "resources/lambda_ARTDECO_SW.txt")
    λ_lw         = _read_column(base_path / "resources/lambda_ARTDECO_LW.txt")

    # ------------------------------------------------------------------ #
    # 3. Cargar plantilla de particles.py
    # ------------------------------------------------------------------ #
    template_particles = _read_template(base_path / "templates/particles.py")

    # ------------------------------------------------------------------ #
    # 4. Desempacar propiedades
    # ------------------------------------------------------------------ #
    geo = m_ext = m_asy = m_ssa = None
    aer_prof = aer_interp = None

    if bool_cld and prop_cld:
        geo, m_ext, m_asy, m_ssa = prop_cld
        # ‣ FIX 2 : asegurar dtype=float para evitar “setting an array element…”
        geo = np.asarray(geo, dtype=float).ravel()
        cbh, cth = float(geo[0]), float(geo[1])
        m_ext = np.asarray(m_ext, dtype=float)
        m_asy = np.asarray(m_asy, dtype=float)
        m_ssa = np.asarray(m_ssa, dtype=float)

    if bool_aer and prop_aer:
        aer_prof, aer_interp = prop_aer
        aer_prof   = np.asarray(aer_prof,   dtype=float)
        aer_interp = np.asarray(aer_interp, dtype=float)

    # ------------------------------------------------------------------ #
    # 5. Preparativos comunes
    # ------------------------------------------------------------------ #
    fmt_date       = date.strftime("%Y%m%d%H")
    levels_artdeco = np.arange(0, 21)              # 0–20 km (resolución 1 km)

    def _add_line(lines: list[str], idx: int, new_line: str) -> None:
        """Reemplaza la línea *idx* (base=1) por `new_line`."""
        lines[idx - 1] = new_line

    # ------------------------------------------------------------------ #
    # 6. Cálculo COD 550 nm por capa de nube (si procede)
    # ------------------------------------------------------------------ #
    v_cod_550: np.ndarray | None = None            # se crea sólo si hay nube
    layers_cld: np.ndarray | None = None

    if bool_cld:
        v_cod_550 = np.zeros(21, dtype=float)
        layers_cld = np.flatnonzero(m_ext[:, 1] != 0).astype(int)

        if layers_cld.size:
            idx_start = layers_cld[0]
            idx_end   = layers_cld[-1]

            # Grosor de la primera y última capa (en km → m)
            dz_start = abs(levels_artdeco[idx_start + 1] - cbh) * 1e3
            dz_end   =       (cth - levels_artdeco[idx_end]) * 1e3

            v_cod_550[idx_start] = m_ext[idx_start, 1] * dz_start
            v_cod_550[idx_end]   = m_ext[idx_end,   1] * dz_end

            # Capas intermedias
            for z in layers_cld[1:-1]:
                dz = (levels_artdeco[z + 1] - levels_artdeco[z]) * 1e3
                v_cod_550[z] = m_ext[z, 1] * dz

    # ------------------------------------------------------------------ #
    # 7. Capas de aerosol con AOD ≠ 0
    # ------------------------------------------------------------------ #
    if bool_aer:
        layers_aer_sw = z_artdeco_sw[np.flatnonzero(aer_prof[:, 1] != 0)]

    # ================================================================== #
    #                       GENERACIÓN DE ARCHIVOS
    # ================================================================== #
    tmp_sw: list[str] | None = None
    tmp_lw: list[str] | None = None

    # ------------- PARTICLES SW -------------
    if bool_sw:
        tmp_sw     = template_particles.copy()
        insert_row = 15                             # tabla comienza en línea 15

        if bool_aer:
            for idx, z in enumerate(layers_aer_sw, start=1):
                _add_line(
                    tmp_sw,
                    insert_row + idx,
                    f'   ["{fmt_date}_{idx}_aer_sw", False, True, True, '
                    f'{aer_prof[idx - 1, 1]:.5f}, "layer", {z:.3f}],',
                )
            insert_row += len(layers_aer_sw)

        if bool_cld:
            for idx, z in enumerate(layers_cld, start=1):
                _add_line(
                    tmp_sw,
                    insert_row + idx,
                    f'   ["{fmt_date}_{idx}_cld_sw", False, True, True, '
                    f'{v_cod_550[z]:.5f}, "layer", {z}],',
                )
            insert_row += len(layers_cld)

        _add_line(tmp_sw, insert_row + 1, "            ]")

    # ------------- PARTICLES LW -------------
    if bool_lw:
        tmp_lw     = template_particles.copy()
        insert_row = 15

        if bool_aer:
            # primera capa (0‑1 km)
            aod_first = aer_prof[z_artdeco_sw < 1, 1].sum()
            _add_line(
                tmp_lw,
                insert_row + 1,
                f'   ["aerosol_lw", False, True, True, {aod_first:.5f}, "layer", 0],',
            )
            insert_row += 1

            # resto de capas
            for idx, z in enumerate(layers_aer_sw[11:], start=12):
                _add_line(
                    tmp_lw,
                    insert_row + idx - 11,
                    f'   ["aerosol_lw", False, True, True, {aer_prof[idx - 1, 1]:.5f}, '
                    f'"layer", {z}],',
                )
            insert_row += max(0, len(layers_aer_sw) - 11)

        if bool_cld:
            for idx, z in enumerate(layers_cld, start=1):
                _add_line(
                    tmp_lw,
                    insert_row + idx,
                    f'   ["{fmt_date}_{idx}_cld_lw", False, True, True, '
                    f'{v_cod_550[z]:.5f}, "layer", {z}],',
                )
            insert_row += len(layers_cld)

        _add_line(tmp_lw, insert_row + 1, "            ]")

    # ------------------------------------------------------------------ #
    # 8. Escribir archivos particles.py
    # ------------------------------------------------------------------ #
    if bool_sw and tmp_sw:
        dir_sw = base_path / "" / str_type / "SW" / f"SW_{fmt_date}" / "config"
        _write_text(dir_sw / "particles.py", tmp_sw)

    if bool_lw and tmp_lw:
        dir_lw = base_path / "INPUTS" / str_type / "LW" / f"LW_{fmt_date}" / "config"
        _write_text(dir_lw / "particles.py", tmp_lw)

    # ------------------------------------------------------------------ #
    # 9. Plantilla común de archivos ópticos (*opt_X.dat*)
    # ------------------------------------------------------------------ #
    opt_template_lines = _read_template(base_path / "opt_X.dat")

    def _write_opt(
        directory: Path,
        file_name: str,
        wavelengths: np.ndarray,
        ext: np.ndarray,
        ssa_: np.ndarray,
        asy_: np.ndarray,
    ) -> None:
        lines      = opt_template_lines.copy()
        lines[10]  = str(len(wavelengths))          # línea 11 en MATLAB → idx 10
        print("Longitud de lines:", len(lines))
        print("Número de longitudes de onda:", len(wavelengths))
        for i, wav in enumerate(wavelengths, start=1):
            lines[10 + i] = (
                f"    {wav:.6e}    190     {ext[i-1]:.6e}     "
                f"{ssa_[i-1]:.6e}    {asy_[i-1]:.6e}"
            )
        _write_text(directory / file_name, lines)

    # ------------------------------------------------------------------ #
    # 10. Archivos ópticos SW
    # ------------------------------------------------------------------ #
    if bool_sw:
        dir_sw = base_path / "INPUTS" / str_type / "SW" / f"SW_{fmt_date}" / "config"

        # nube
        if bool_cld:
            for idx, z in enumerate(layers_cld, start=1):
                _write_opt(
                    dir_sw,
                    f"opt_{fmt_date}_{idx}_cld_sw.dat",
                    λ_sw,
                    m_ext[z],
                    m_ssa[z],
                    m_asy[z],
                )

        # aerosol
        if bool_aer:
            for idx, z in enumerate(layers_aer_sw, start=1):
                dz = z_artdeco_sw[idx] - z_artdeco_sw[idx - 1 if idx > 1 else 0]
                ext_layer = aer_interp[1] / (dz * 1e3)        # [m⁻¹]
                _write_opt(
                    dir_sw,
                    f"opt_{fmt_date}_{idx}_aer_sw.dat",
                    aer_interp[0],
                    ext_layer,
                    aer_interp[2],
                    aer_interp[3],
                )

    # ------------------------------------------------------------------ #
    # 11. Archivos ópticos LW  (sólo nube)
    # ------------------------------------------------------------------ #
    if bool_lw and bool_cld:
        dir_lw = base_path / "INPUTS" / str_type / "LW" / f"LW_{fmt_date}" / "config"

        for idx, z in enumerate(layers_cld, start=1):
            # λ_lw + canal 0.55 µm (posición λ_sw[3])
            wavs = np.insert(λ_lw, 0, λ_sw[3])
            ext  = np.insert(m_ext[z, 7:], 0, m_ext[z, 3])
            ssa_ = np.insert(m_ssa[z, 7:], 0, m_ssa[z, 3])
            asy_ = np.insert(m_asy[z, 7:], 0, m_asy[z, 3])

            _write_opt(
                dir_lw,
                f"opt_{fmt_date}_cld_case_lw.dat",
                wavs,
                ext,
                ssa_,
                asy_,
            )