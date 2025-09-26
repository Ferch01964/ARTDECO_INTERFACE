# prop_rad_aer_artdeco.py
import os
import glob
from datetime import datetime, timedelta
from mplnet import download_mplnet_range
import numpy as np
import pandas as pd
from netCDF4 import Dataset


def prop_rad_aer_artdeco(path: str, date: datetime):
    """
    Conversión depurada de la función MATLAB `prop_rad_aer_artdeco`.

    Parameters
    ----------
    path : str
        Ruta con la carpeta INPUTS/resources y los ficheros auxiliares.
    date : datetime
        Instante de la escena.

    Returns
    -------
    data_aeronet_prof : ndarray (Nz, 3)
        v_weight | v_aod_550 | v_ext_550
    data_aeronet_interpn : ndarray (4, M)
        λ | AOD | SSA | ASY (interpolados a λ_ARTDECO)
    """
    # ------------------------------------------------------------------ #
    # 1. Perfil vertical y longitudes de onda de ARTDECO
    # ------------------------------------------------------------------ #
    z_artdeco_sw = pd.read_csv(
        os.path.join(path, "resources", "z_p_ARTDECO_SW.txt"),
        sep=r'\s+', header=None
    ).iloc[:, 0].to_numpy(dtype=float)  # km

    lambda_artdeco_sw = pd.read_csv(
        os.path.join(path, "resources", "lambda_ARTDECO_SW.txt"),
        sep=r'\s+', header=None
    ).iloc[:, 0].to_numpy(dtype=float)  # µm
    M = lambda_artdeco_sw.size

    # ------------------------------------------------------------------ #
    # 2. Lectura de fichero AERONET (Level‑15, inversión)
    # ------------------------------------------------------------------ #
    aeronet_files = glob.glob(os.path.join(
        path,  "resources", "*AERONET_INVERSION_L15*"))
    if not aeronet_files:
        raise FileNotFoundError("No se encontró *AERONET_INVERSION_L15* en INPUTS/resources")
    aeronet_file = aeronet_files[0]

    usecols = [
        "Date(dd:mm:yyyy)", "Time(hh:mm:ss)",
        "Surface_Albedo[440m]", "Surface_Albedo[675m]",
        "Surface_Albedo[870m]", "Surface_Albedo[1020m]",
        "Single_Scattering_Albedo[440nm]", "Single_Scattering_Albedo[675nm]",
        "Single_Scattering_Albedo[870nm]", "Single_Scattering_Albedo[1020nm]",
        "Asymmetry_Factor-Total[440nm]", "Asymmetry_Factor-Total[675nm]",
        "Asymmetry_Factor-Total[870nm]", "Asymmetry_Factor-Total[1020nm]",
    ]
    df_aero = pd.read_csv(aeronet_file, skiprows=6, usecols=usecols)

    # -- fechas a datetime
    df_aero["Date(dd:mm:yyyy)"] = pd.to_datetime(
        df_aero["Date(dd:mm:yyyy)"], format="%d:%m:%Y", errors="coerce")
    df_aero["Time(hh:mm:ss)"] = pd.to_datetime(
        df_aero["Time(hh:mm:ss)"], format="%H:%M:%S", errors="coerce").dt.time
    df_aero["DatetimeAERONET"] = [
        datetime.combine(d, t) if pd.notnull(d) and pd.notnull(t) else np.nan
        for d, t in zip(df_aero["Date(dd:mm:yyyy)"], df_aero["Time(hh:mm:ss)"])
    ]

    # -- estilo MATLAB datenum
    epoch = datetime(1970, 1, 1)
    datenum = lambda dt: (dt - epoch).total_seconds() / 86400.0 if dt is not np.nan else np.nan
    data_times_num = df_aero["DatetimeAERONET"].map(datenum).to_numpy(dtype=float)

    # -- helper para columnas float
    def col(name): return df_aero[name].to_numpy(dtype=float)

    aod = np.vstack([col(c) for c in (
        "Surface_Albedo[440m]", "Surface_Albedo[675m]",
        "Surface_Albedo[870m]", "Surface_Albedo[1020m]")])
    ssa = np.vstack([col(c) for c in (
        "Single_Scattering_Albedo[440nm]", "Single_Scattering_Albedo[675nm]",
        "Single_Scattering_Albedo[870nm]", "Single_Scattering_Albedo[1020nm]")])
    asy = np.vstack([col(c) for c in (
        "Asymmetry_Factor-Total[440nm]", "Asymmetry_Factor-Total[675nm]",
        "Asymmetry_Factor-Total[870nm]", "Asymmetry_Factor-Total[1020nm]")])

    data_aeronet = np.vstack([data_times_num, aod, ssa, asy])        # (13, N)
    lambda_aeronet = np.array([0.440, 0.675, 0.870, 1.020])

    # ------------------------------------------------------------------ #
    # 3. Índice temporal AERONET más cercano a `date`
    # ------------------------------------------------------------------ #
    target_num = datenum(date)
    idx_time, n_tries = [], 1
    while not idx_time and n_tries < 10:
        if 6 < date.hour < 18:                                       # horario diurno
            diffs = np.abs(data_aeronet[0] - target_num)
            idx_time = np.where(diffs == np.nanmin(diffs))[0]
        else:                                                        # ±n días
            lo = datenum(date - timedelta(days=n_tries))
            hi = datenum(date + timedelta(days=n_tries))
            idx_time = np.where((data_aeronet[0] >= lo) & (data_aeronet[0] <= hi))[0]
        n_tries += 1

    # ------------------------------------------------------------------ #
    # 4. Interpolación espectral AERONET → λ_ARTDECO
    # ------------------------------------------------------------------ #
    data_aeronet_interpn = np.full((4, M), np.nan, dtype=float)
    data_aeronet_interpn[0] = lambda_artdeco_sw

    if idx_time:                                                    # hay datos válidos
        for j, lam in enumerate(lambda_artdeco_sw):
            k = np.where(np.isclose(lambda_aeronet, lam))[0]
            if k.size:
                k = k[0]
                data_aeronet_interpn[1, j] = np.nanmean(data_aeronet[1 + k, idx_time])      # AOD
                data_aeronet_interpn[2, j] = np.nanmean(data_aeronet[1 + 4 + k, idx_time])  # SSA
                data_aeronet_interpn[3, j] = np.nanmean(data_aeronet[1 + 8 + k, idx_time])  # ASY

        # --- completar AOD con ley de Ångström ---
        ind_ok = np.flatnonzero(~np.isnan(data_aeronet_interpn[1]))
        ind_nan = np.flatnonzero(np.isnan(data_aeronet_interpn[1]))
        if ind_ok.size >= 2:
            for j in ind_nan:
                dif = np.abs(data_aeronet_interpn[0, ind_ok] - data_aeronet_interpn[0, j])
                k1, k2 = ind_ok[np.argsort(dif)[:2]]
                a1, a2 = data_aeronet_interpn[1, [k1, k2]]
                l1, l2 = data_aeronet_interpn[0, [k1, k2]]
                if a1 > 0 and a2 > 0 and l1 != l2:
                    ae = -np.log(a2 / a1) / np.log(l2 / l1)
                    data_aeronet_interpn[1, j] = a1 * (l1 / data_aeronet_interpn[0, j]) ** ae

        # --- completar SSA y ASY por interpolación lineal ---
        from numpy import interp
        lam_ok = data_aeronet_interpn[0, ind_ok]
        for row in (2, 3):
            val_ok = data_aeronet_interpn[row, ind_ok]
            data_aeronet_interpn[row, ind_nan] = interp(
                data_aeronet_interpn[0, ind_nan], lam_ok, val_ok, left=np.nan, right=np.nan
            )

    # ------------------------------------------------------------------ #
    # 5. Lectura del producto MPLNET (extinction)
    # ------------------------------------------------------------------ #
    nc1_files_mplnet = download_mplnet_range(date)
    print("Archivos extraídos:")
    for p in nc1_files_mplnet:
        print(" •", p)
    print(f"MPLNET path files: {nc1_files_mplnet[0]}")
    nc_files = glob.glob(os.path.join(path,  "resources", "*prove*"))
    if not nc_files:
        raise FileNotFoundError("No se encontró ningún netCDF de MPLNET.")
    nc_file = nc_files[0]

    with Dataset(nc_file) as nc:
        extinction = nc.variables["extinction"][:]      # MaskedArray
        altitude   = nc.variables["altitude"][:]

    # -- ajustar orden (height, time)
    if extinction.shape[0] == altitude.shape[-1]:       # llegada (time,height)
        extinction = extinction.T

    altitude = np.ma.filled(altitude, np.nan).astype(float)
    altitude = altitude[:, 0] if altitude.ndim == 2 else altitude
    altitude = altitude - altitude[0]

    # -- serie temporal 1 min del día
    start_day = datetime(date.year, date.month, date.day)
    time_day_num = np.array([datenum(start_day + timedelta(minutes=i)) for i in range(1440)])

    # ------------------------------------------------------------------ #
    # 6. Promedio de extinción en un rango ±30 min·k
    # ------------------------------------------------------------------ #
    idx_time2, bool_ext = np.array([], int), False
    for k_loop in range(5):                         # k = 0,1,2,3,4  (0–120 min)
        lo = datenum(date - timedelta(minutes=30 * k_loop))
        hi = datenum(date + timedelta(minutes=30 * k_loop))
        idx_time2 = np.where((time_day_num >= lo) & (time_day_num <= hi))[0]
        if idx_time2.size:
            ext_mean = np.ma.mean(extinction[:, idx_time2], axis=1)
            if np.any(ext_mean <= 0):
                bool_ext = True
                break

    # ------------------------------------------------------------------ #
    # 7. Perfil medio de extinción y altitud top SW
    # ------------------------------------------------------------------ #
    extinction_prof, ind_top_sw = None, 12
    if bool_ext:
        ext_mean = np.ma.mean(extinction[:, idx_time2], axis=1)
        # ► conversión segura a vector 1‑D de floats
        extinction_prof = np.ma.filled(ext_mean, np.nan).astype(float).ravel()

        aux = np.where((extinction_prof <= 0) | np.isnan(extinction_prof))[0]
        if aux.size:
            ind_top_sw = np.argmin(np.abs(z_artdeco_sw - altitude[aux[0]]))

    # ------------------------------------------------------------------ #
    # 8. Construcción del perfil 550 nm (v_weight, v_ext_550, v_aod_550)
    # ------------------------------------------------------------------ #
    Nz = z_artdeco_sw.size
    v_weight = np.zeros(Nz)
    v_ext_550 = np.zeros(Nz)
    v_aod_550 = np.zeros(Nz)
    ind_base = 0

    if extinction_prof is not None:                    # con MPLNET
        # capa base
        v_ext_550[ind_base] = float(extinction_prof[0]) * 1e-3 * z_artdeco_sw[0] / 0.0749
        v_aod_550[ind_base] = v_ext_550[ind_base] * z_artdeco_sw[0] * 1e3

        # 0‑0.06 km (índices 1‑6, paso 0.01 km)
        for z in range(1, 7):
            dz = z_artdeco_sw[z] - z_artdeco_sw[z - 1]
            v_ext_550[z] = float(extinction_prof[0]) * 1e-3 * dz / 0.0749
            v_aod_550[z] = v_ext_550[z] * dz * 1e3

        # resto hasta ind_top_sw
        for z in range(7, ind_top_sw + 1):
            mask = (altitude >= z_artdeco_sw[z - 1]) & (altitude < z_artdeco_sw[z])
            val  = np.nanmean(extinction_prof[mask]) if mask.any() else np.nan
            v_ext_550[z] = val * 1e-3
            dz = z_artdeco_sw[z] - z_artdeco_sw[z - 1]
            v_aod_550[z] = v_ext_550[z] * dz * 1e3

        total_aod = np.nansum(v_aod_550[: ind_top_sw + 1])
        if total_aod:
            v_weight[: ind_top_sw + 1] = v_aod_550[: ind_top_sw + 1] / total_aod
    else:                                             # sin MPLNET: perfil plano
        aod_550_equiv = data_aeronet_interpn[1, 2] if M > 2 else 0.0
        v_ext_550[:] = aod_550_equiv / z_artdeco_sw[-1] * 1e-3
        v_aod_550[:] = aod_550_equiv / z_artdeco_sw[-1]
        v_weight[:]  = v_aod_550 / max(np.nansum(v_aod_550), 1e-12)

    # -- remplaza NaN por 0
    for arr in (v_weight, v_ext_550, v_aod_550):
        arr[np.isnan(arr)] = 0.0

    data_aeronet_prof = np.column_stack([v_weight, v_aod_550, v_ext_550])
    return data_aeronet_prof, data_aeronet_interpn
