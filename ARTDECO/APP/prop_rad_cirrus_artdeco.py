import os
import numpy as np
import pandas as pd
import sympy as sp
from scipy.interpolate import interp1d

def prop_rad_cirrus_artdeco(path, prop_geo, v_h_lidar, v_temp_cld, v_ext_cld):
    """
    Función que crea las propiedades ópticas de dispersión: ext, asy y SSA con la parametrización de Baran 2018 
    en el espectro SW/LW, tal como en el código MATLAB original.
    
    Entradas:
      - path: Ruta principal donde se encuentra la carpeta INPUTS.
      - prop_geo: Tuple (cbh, cth) con la altura base y tope del cirrus [km].
      - v_h_lidar: Perfil vertical de alturas del MPL [km].
      - v_temp_cld: Perfil vertical de temperatura en la nube [K].
      - v_ext_cld: Perfil vertical de coeficientes de extinción [m^-1].
    
    Salidas:
      - m_ext_cld: Matriz de coeficientes de extinción [m^-1] (filas: niveles, columnas: longitudes).
      - m_asy_cld: Matriz de factores de asimetría.
      - m_ssa_cld: Matriz de albedo de dispersión simple.
      - m_sc_cld: Matriz de dispersión [m^-1].
      - m_abs_cld: Matriz de absorción [m^-1].
      - v_iwc_cld: Vector de contenido de agua en forma de hielo [kg/m^3].
    """

    # --- Bloque 1: Lectura de coeficientes de dispersión ---
    file_sca_coef = os.path.join(path,  "resources", "sca_coef_Baran.txt")
    df_sca = pd.read_csv(file_sca_coef, delimiter=';', header=None, engine='python', dtype=str)
    df_sca = df_sca.apply(lambda col: col.map(lambda x: pd.to_numeric(x.replace(',', ''), errors='coerce') 
                                              if isinstance(x, str) else x))
    # Se usan las primeras 6 columnas
    coef_sc = df_sca.iloc[:, :6].to_numpy()

    # --- Bloque 2: Lectura de longitudes de onda de "lambda_Baran.txt" ---
    file_lambda_baran = os.path.join(path,  "resources", "lambda_Baran.txt")
    df_lambda = pd.read_csv(file_lambda_baran, delimiter=';', header=None, engine='python', dtype=float)
    
    # Concatenar todas las columnas (29 filas × 5 columnas = 145 elementos)
    lambda_baran_all = np.concatenate([df_lambda[col].values for col in df_lambda.columns])
    lambda_baran = np.sort(lambda_baran_all)[::-1]

    # Verificar la dimensión
    assert len(lambda_baran) == coef_sc.shape[0], (
        f"El vector de longitudes de onda tiene {len(lambda_baran)} elementos, "
        f"pero coef_sc tiene {coef_sc.shape[0]} filas. Revisa el formato de lambda_Baran.txt."
    )

    # --- Bloque 3: Lectura de coeficientes de absorción ---
    file_abs_coef = os.path.join(path,  "resources", "abs_coef_Baran.txt")
    df_abs = pd.read_csv(file_abs_coef, delimiter=';', header=None, engine='python', dtype=str)
    df_abs = df_abs.apply(lambda col: col.map(lambda x: pd.to_numeric(x.replace('/)', '').replace(',', ''), errors='coerce') 
                                              if isinstance(x, str) else x))
    coef_abs = df_abs.iloc[:, :6].to_numpy()

    # --- Bloque 4: Lectura de coeficientes de asimetría ---
    file_asym_coef = os.path.join(path,  "resources", "asym_coef_Baran.txt")
    df_asym = pd.read_csv(file_asym_coef, delimiter=';', header=None, engine='python', dtype=str)
    df_asym = df_asym.apply(lambda col: col.map(lambda x: pd.to_numeric(x.replace(',', ''), errors='coerce') 
                                                if isinstance(x, str) else x))
    coef_g = df_asym.iloc[:, :3].to_numpy()

    # --- Bloque 5: Lectura de longitudes de onda ARTDECO SW ---
    file_lambda_artdeco_sw = os.path.join(path,  "resources", "lambda_ARTDECO_SW.txt")
    df_lambda_sw = pd.read_csv(file_lambda_artdeco_sw, delimiter=';', header=None, engine='python', dtype=float)
    lambda_artdeco_sw = df_lambda_sw.iloc[:, 0].to_numpy()

    # --- Bloque 6: Lectura de longitudes de onda ARTDECO LW ---
    file_lambda_artdeco_lw = os.path.join(path,  "resources", "lambda_ARTDECO_LW.txt")
    df_lambda_lw = pd.read_csv(file_lambda_artdeco_lw, delimiter=';', header=None, engine='python', dtype=float)
    lambda_artdeco_lw = df_lambda_lw.iloc[:, 0].to_numpy()

    # --- Bloque 7: Lectura de niveles ARTDECO LW ---
    file_z_artdeco = os.path.join(path,  "resources", "z_p_ARTDECO_LW.txt")
    df_z = pd.read_csv(file_z_artdeco, delimiter=r"\s+", header=None, engine='python', dtype=float)
    z_artdeco_lw = df_z.iloc[:, 0].to_numpy()

    # --- Bloque 8: Definición de funciones simbólicas con sympy ---
    t, a, coef_1, coef_2, coef_3, coef_4, coef_5, coef_6 = sp.symbols('t a coef_1 coef_2 coef_3 coef_4 coef_5 coef_6', real=True)
    f_liwc_expr = -0.5*(coef_3 + coef_6*t)/coef_5 + 0.5*sp.sqrt(((coef_3+coef_6*t)/coef_5)**2 + 4*(sp.log(a,10) - coef_1 - coef_2*t - coef_4*t*t)/coef_5)
    f_liwc_func = sp.lambdify((t, a, coef_1, coef_2, coef_3, coef_4, coef_5, coef_6), f_liwc_expr, modules=["numpy"])

    liwc_sym = sp.symbols('liwc', real=True)
    f3_coef_expr = coef_1 + coef_2*t + coef_3*liwc_sym
    f3_coef_func = sp.lambdify((t, liwc_sym, coef_1, coef_2, coef_3), f3_coef_expr, modules=["numpy"])

    fl6_coef_expr = coef_1 + coef_2*t + coef_3*liwc_sym + coef_4*t*t + coef_5*liwc_sym**2 + coef_6*t*liwc_sym
    fl6_coef_func = sp.lambdify((t, liwc_sym, coef_1, coef_2, coef_3, coef_4, coef_5, coef_6), fl6_coef_expr, modules=["numpy"])

    # --- Bloque 9: Conversión de altitudes a ARTDECO y cálculo de perfiles ---
    cbh = prop_geo[0]
    cth = prop_geo[1]
    cbh_mpl = np.argmin(np.abs(v_h_lidar - cbh))
    cth_mpl = np.argmin(np.abs(v_h_lidar - cth))
    v_h_lidar_cloud = v_h_lidar[cbh_mpl:cth_mpl+1]

    indices_base = np.where(z_artdeco_lw <= np.min(v_h_lidar_cloud))[0]
    base_artdeco = indices_base[-1] if len(indices_base) > 0 else 0
    indices_top = np.where(z_artdeco_lw >= np.max(v_h_lidar_cloud))[0]
    top_artdeco = indices_top[0] if len(indices_top) > 0 else len(z_artdeco_lw)-1

    v_temp_cld_cloud = v_temp_cld[cbh_mpl:cth_mpl+1]
    v_ext_cld_cloud = v_ext_cld[cbh_mpl:cth_mpl+1]

    v_temp_cld_artdeco = np.full(21, np.nan)
    v_ext_cld_artdeco = np.full(21, np.nan)
    for z_p in range(base_artdeco, top_artdeco):
        mask = (v_h_lidar_cloud >= z_artdeco_lw[z_p]) & (v_h_lidar_cloud < z_artdeco_lw[z_p+1])
        if np.any(mask):
            v_temp_cld_artdeco[z_p] = np.nanmean(v_temp_cld_cloud[mask])
            v_ext_cld_artdeco[z_p] = np.nanmean(v_ext_cld_cloud[mask])

    # --- Bloque 10: Cálculo de IWC ---
    l_val = 0.532
    lambda_index = np.argmin(np.abs(lambda_baran - l_val))
    coef_sc_row = coef_sc[lambda_index, :]
    f_liwc_values = f_liwc_func(v_temp_cld_artdeco, v_ext_cld_artdeco,
                                coef_sc_row[0], coef_sc_row[1], coef_sc_row[2],
                                coef_sc_row[3], coef_sc_row[4], coef_sc_row[5])
    iwc_sol_artdeco = 10 ** (np.array(f_liwc_values))
    v_iwc_cld = iwc_sol_artdeco

    # --- Bloque 11: Cálculo de las matrices espectrales ---
    lambda_game = np.concatenate([lambda_artdeco_sw, lambda_artdeco_lw])
    m_ext_cld = np.zeros((21, len(lambda_game)))
    m_asy_cld = np.zeros((21, len(lambda_game)))
    m_ssa_cld = np.zeros((21, len(lambda_game)))
    m_sc_cld = np.zeros((21, len(lambda_game)))
    m_abs_cld = np.zeros((21, len(lambda_game)))
    
    for j in range(21):
        if not np.isnan(v_temp_cld_artdeco[j]):
            log_iwc = np.log10(iwc_sol_artdeco[j]) if iwc_sol_artdeco[j] > 0 else -np.inf
            v_sc = 10 ** (fl6_coef_func(v_temp_cld_artdeco[j], log_iwc,
                                        coef_sc[:, 0], coef_sc[:, 1], coef_sc[:, 2],
                                        coef_sc[:, 3], coef_sc[:, 4], coef_sc[:, 5]))
            v_abs = 10 ** (fl6_coef_func(v_temp_cld_artdeco[j], log_iwc,
                                         coef_abs[:, 0], coef_abs[:, 1], coef_abs[:, 2],
                                         coef_abs[:, 3], coef_abs[:, 4], coef_abs[:, 5]))
            v_asy = f3_coef_func(v_temp_cld_artdeco[j],
                                 np.log10(iwc_sol_artdeco[j]) if iwc_sol_artdeco[j] > 0 else -np.inf,
                                 coef_g[:, 0], coef_g[:, 1], coef_g[:, 2])
            v_ssa = v_sc / (v_sc + v_abs)
            v_ext = v_sc + v_abs

            lambda_baran_asc = lambda_baran[::-1]
            v_ext_asc = v_ext[::-1]
            v_asy_asc = v_asy[::-1]
            v_ssa_asc = v_ssa[::-1]
            v_sc_asc = v_sc[::-1]
            v_abs_asc = v_abs[::-1]

            f_ext = interp1d(lambda_baran_asc, v_ext_asc, kind='nearest', bounds_error=False, fill_value="extrapolate")
            f_asy = interp1d(lambda_baran_asc, v_asy_asc, kind='nearest', bounds_error=False, fill_value="extrapolate")
            f_ssa = interp1d(lambda_baran_asc, v_ssa_asc, kind='nearest', bounds_error=False, fill_value="extrapolate")
            f_sc = interp1d(lambda_baran_asc, v_sc_asc, kind='nearest', bounds_error=False, fill_value="extrapolate")
            f_abs = interp1d(lambda_baran_asc, v_abs_asc, kind='nearest', bounds_error=False, fill_value="extrapolate")
            
            m_ext_cld[j, :] = np.real(f_ext(lambda_game))
            m_asy_cld[j, :] = np.real(f_asy(lambda_game))
            m_ssa_cld[j, :] = np.real(f_ssa(lambda_game))
            m_sc_cld[j, :] = np.real(f_sc(lambda_game))
            m_abs_cld[j, :] = np.real(f_abs(lambda_game))
    
    return m_ext_cld, m_asy_cld, m_ssa_cld, m_sc_cld, m_abs_cld, v_iwc_cld
