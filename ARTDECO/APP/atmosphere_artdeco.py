import os
import numpy as np
import datetime
from datetime import timedelta
import glob
import pandas as pd
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from radiosoun import radiosoun

def atmosphere_artdeco(path, date, bool_cld, bool_aer, bool_sw, bool_lw, lat_place, lon_place):
    """
    Crea los archivos atmosféricos para ARTDECO (SW y LW) a partir de insumos y perfiles termodinámicos.
    
    Parámetros:
      - path (str): Ruta base donde se encuentran los archivos de entrada y se generarán los archivos.
      - date (datetime.datetime): Fecha y hora de la simulación.
      - bool_cld (bool): Indica si la escena atmosférica tiene nubes.
      - bool_aer (bool): Indica si la escena atmosférica tiene aerosoles.
      - bool_sw (bool): Indica si se crea el archivo para simulación SW.
      - bool_lw (bool): Indica si se crea el archivo para simulación LW.
      - lat_place (float): Latitud del lugar de la simulación.
      - lon_place (float): Longitud del lugar de la simulación.
    
    Se mantiene la misma estructura de directorios que en MATLAB y se realizan los siguientes pasos:
      1. Se determina el tipo de atmósfera a usar según la presencia de nubes y aerosoles.
      2. Se lee y modifica la plantilla de archivo (atmosphere.py) para cada tipo (SW/LW).
      3. Se leen los perfiles verticales (alturas) de los archivos correspondientes.
      4. Se leen los datos de concentración de ozono (si existen) desde un archivo NetCDF.
      5. Se leen el perfil estándar de atmósfera y los datos de radiosondeo (RS) a partir de archivos.
      6. Se calcula el perfil vertical de contenido de vapor (DH en ppmv) y se interpola
         (usando interpolación de vecino más cercano) para los perfiles SW y LW.
      7. Se corrigen valores faltantes y se reordena el perfil en orden descendente de altitud.
      8. Finalmente, se escribe el archivo final de salida con un encabezado modificado.
    """
    # 1. Determinar el tipo de escena atmosférica según nubes y aerosoles
    if bool_cld and bool_aer:
        str_type = 'GAC'
    elif bool_cld and not bool_aer:
        str_type = 'CLD'
    elif not bool_cld and bool_aer:
        str_type = 'AER'
    else:
        str_type = 'GAS'
    
    # 2. Lectura de la plantilla "atmosphere.py" que se modificará según el caso
    template_path = os.path.join(path,  'templates','atmosphere.py')
    with open(template_path, 'r') as f:
        template_lines = f.readlines()
    
    def update_template_line(sim_type):
        """
        Actualiza la línea 10 (índice 9) de la plantilla para asignar el nombre del perfil atmosférico.
        sim_type: 'sw' o 'lw'
        """
        new_line = f'Atmospheric_profile = "{date.strftime("%Y%m%d%H")}_{sim_type}_ppmv"\n'
        updated_lines = template_lines.copy()
        if len(updated_lines) >= 10:
            updated_lines[9] = new_line  # Se reemplaza la línea 10
        else:
            updated_lines.append(new_line)
        return updated_lines
    
    # 3. Creación de los archivos de configuración SW y LW según los booleanos
    if bool_sw:
        dir_sw = os.path.join(path,  str_type, 'SW', f'SW_{date.strftime("%Y%m%d%H")}', 'config')
        os.makedirs(dir_sw, exist_ok=True)
        updated_lines_sw = update_template_line('sw')
        sw_file_path = os.path.join(dir_sw, 'atmosphere.py')
        with open(sw_file_path, 'w') as f:
            f.writelines(updated_lines_sw)
    
    if bool_lw:
        dir_lw = os.path.join(path,  str_type, 'LW', f'LW_{date.strftime("%Y%m%d%H")}', 'config')
        os.makedirs(dir_lw, exist_ok=True)
        updated_lines_lw = update_template_line('lw')
        lw_file_path = os.path.join(dir_lw, 'atmosphere.py')
        with open(lw_file_path, 'w') as f:
            f.writelines(updated_lines_lw)
    
    # 4. Lectura de los niveles verticales (alturas) en km desde archivos de texto
    if bool_sw:
        sw_levels_file = os.path.join(path, 'resources', 'z_p_ARTDECO_SW.txt')
        z_artdeco_sw = np.loadtxt(sw_levels_file, usecols=0)
    if bool_lw:
        lw_levels_file = os.path.join(path,  'resources','z_p_ARTDECO_LW.txt')
        z_artdeco_lw = np.loadtxt(lw_levels_file, usecols=0)
    
    # 5. Lectura del perfil vertical de ozono desde archivo NetCDF (CAMS_O3.nc)
    o3_file_path = os.path.join(path, 'resources', 'CAMS_O3.nc')
    bool_o3 = False
    if os.path.exists(o3_file_path):
        nc = Dataset(o3_file_path, 'r')
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
        level = nc.variables['level'][:]
        aux_hr = nc.variables['time'][:]
        base_time = datetime.datetime(1900, 1, 1)
        # Conversión de horas a objetos datetime
        time_vals = np.array([base_time + timedelta(hours=float(hr)) for hr in aux_hr])
        # Búsqueda del punto más cercano en longitud y latitud
        ind_lon = np.argmin(np.abs(lon - lon_place))
        ind_lat = np.argmin(np.abs(lat - lat_place))
        # Búsqueda del índice temporal que coincida con la fecha (a nivel de año, mes, día y hora)
        time_match = [i for i, t in enumerate(time_vals) 
                      if (t.year == date.year and t.month == date.month and t.day == date.day and t.hour == date.hour)]
        if time_match:
            ind_time = time_match[0]
            data_o3 = nc.variables['go3'][ind_lon, ind_lat, :, ind_time]
            data_o3 = data_o3 * 28.97 / 48  # Conversión de unidades
            bool_o3 = True
        nc.close()
    
    # 6. Lectura del perfil estándar de atmósfera (us67)
    thermo_files = glob.glob(os.path.join(path,  'resources','thermo_profil_us67*'))
    if thermo_files:
        df_us67 = pd.read_csv(thermo_files[0], delimiter='\t', skiprows=1, header=None,
                              names=['z', 'p', 'T', 'h2o', 'o3'])
        data_us62 = df_us67[['z', 'p', 'T', 'h2o', 'o3']].to_numpy()
    else:
        data_us62 = None
    
    # 7. Lectura de los datos de radiosondeo (RS) mediante la función auxiliar reading_rs
    path_rs = os.path.join(path,  'resources','Radiosoundings')
    data_rs = reading_rs(date, path_rs)
    # Se ajusta la altitud para que el nivel más bajo sea 0
    data_rs[:, 0] = data_rs[:, 0] - data_rs[0, 0]
    
    # 8. Interpolación del perfil de ozono (si está disponible)
    if bool_o3:
        z_cams = []
        # Para niveles altos: usar los 7 primeros niveles de 'level'
        for l in level[:7]:
            idx = np.argmin(np.abs(data_us62[:, 1] - l))
            z_cams.append(data_us62[idx, 0])
        # Para niveles bajos: usar altitudes de RS para el resto de niveles (convertidas a km)
        for l in level[7:]:
            idx = np.argmin(np.abs(data_rs[:, 1] - l))
            z_cams.append(data_rs[idx, 0] * 1e-3)
        z_cams = np.array(z_cams)
        z_unique = np.unique(z_cams)
        data_o3_unique = []
        for z_val in z_unique:
            indices = np.where(z_cams == z_val)[0]
            data_o3_unique.append(np.nanmean(data_o3[indices]))
        z_cams_unique = z_unique
        data_o3_unique = np.array(data_o3_unique)
    
    # 9. Cálculo del perfil vertical DH (ppmv) a partir de datos de RS
    # Esquema: es -> ep -> w -> dh
    es = 6.108 * np.exp(17.269 * (data_rs[:, 2] - 273.15) / ((data_rs[:, 2] - 273.15) + 237.19))
    ep = (data_rs[:, 3] / 100.0) * es
    w = 621.5 * ep / (data_rs[:, 1] - ep)
    dh = w * 1e-3 * 18 / 28.97
    dh[dh < 0] = np.nan  # Se reemplazan valores negativos por NaN
    
    # 10. Interpolación de perfiles termodinámicos para SW y LW (usando vecino más cercano)
    def interpolate_profile(z_target, x, y):
        """
        Interpola los valores de y para los puntos z_target a partir de la grilla x usando el método 'nearest'.
        """
        interp_func = interp1d(x, y, kind='nearest', fill_value="extrapolate")
        return interp_func(z_target)
    
    # Perfil SW
    if bool_sw:
        # Se obtienen altitudes únicas (de RS) en km
        rs_alt_km, unique_indices = np.unique(data_rs[:, 0] * 1e-3, return_index=True)
        p_sw = interpolate_profile(z_artdeco_sw, data_rs[unique_indices, 0] * 1e-3, data_rs[unique_indices, 1])
        T_sw = interpolate_profile(z_artdeco_sw, data_rs[unique_indices, 0] * 1e-3, data_rs[unique_indices, 2])
        dh_sw = interpolate_profile(z_artdeco_sw, data_rs[unique_indices, 0] * 1e-3, dh[unique_indices])
        if bool_o3:
            o3_sw = interpolate_profile(z_artdeco_sw, z_cams_unique, data_o3_unique)
            profil_thermo_sw = np.column_stack((z_artdeco_sw, p_sw, T_sw, dh_sw, o3_sw))
        else:
            o3_sw = np.full_like(z_artdeco_sw, np.nan)
            profil_thermo_sw = np.column_stack((z_artdeco_sw, p_sw, T_sw, dh_sw, o3_sw))
    
    # Perfil LW
    if bool_lw:
        rs_alt_km, unique_indices = np.unique(data_rs[:, 0] * 1e-3, return_index=True)
        p_lw = interpolate_profile(z_artdeco_lw, data_rs[unique_indices, 0] * 1e-3, data_rs[unique_indices, 1])
        T_lw = interpolate_profile(z_artdeco_lw, data_rs[unique_indices, 0] * 1e-3, data_rs[unique_indices, 2])
        dh_lw = interpolate_profile(z_artdeco_lw, data_rs[unique_indices, 0] * 1e-3, dh[unique_indices])
        if bool_o3:
            o3_lw = interpolate_profile(z_artdeco_lw, z_cams_unique, data_o3_unique)
            profil_thermo_lw = np.column_stack((z_artdeco_lw, p_lw, T_lw, dh_lw, o3_lw))
        else:
            o3_lw = np.full_like(z_artdeco_lw, np.nan)
            profil_thermo_lw = np.column_stack((z_artdeco_lw, p_lw, T_lw, dh_lw, o3_lw))
    
    # 11. Corrección de valores faltantes y continuidad en el perfil
    def fix_profile(profile, data_std):
        """
        Rellena valores NaN en presión y temperatura usando el perfil estándar, 
        realiza 'forward fill' y 'backward fill' para DH y O3, y corrige duplicados en presión.
        """
        # Reemplazar NaN en presión (columna 1) y temperatura (columna 2)
        for col in [1, 2]:
            nan_idx = np.where(np.isnan(profile[:, col]))[0]
            if nan_idx.size > 0:
                interp_std = interp1d(data_std[:, 0], data_std[:, col], kind='nearest', fill_value="extrapolate")
                profile[nan_idx, col] = interp_std(profile[nan_idx, 0])
        # Ajuste de DH (columna 3)
        for i in range(len(profile)):
            if np.isnan(profile[i, 3]) or profile[i, 3] == 0:
                if i > 0:
                    profile[i, 3] = profile[i-1, 3]
        for i in range(len(profile)-1, -1, -1):
            if np.isnan(profile[i, 3]) or profile[i, 3] == 0:
                if i < len(profile)-1:
                    profile[i, 3] = profile[i+1, 3]
        # Ajuste de O3 (columna 4)
        if bool_o3:
            for i in range(len(profile)):
                if np.isnan(profile[i, 4]) or profile[i, 4] == 0:
                    if i > 0:
                        profile[i, 4] = profile[i-1, 4]
            for i in range(len(profile)-1, -1, -1):
                if np.isnan(profile[i, 4]) or profile[i, 4] == 0:
                    if i < len(profile)-1:
                        profile[i, 4] = profile[i+1, 4]
        else:
            interp_std = interp1d(data_std[:, 0], data_std[:, 4], kind='nearest', fill_value="extrapolate")
            nan_idx = np.where(np.isnan(profile[:, 4]))[0]
            profile[nan_idx, 4] = interp_std(profile[nan_idx, 0])
        # Ajuste de presiones duplicadas: se añade un incremento muy pequeño
        for i in range(len(profile)):
            dup_idx = np.where(np.isclose(profile[i, 1], profile[:, 1]))[0]
            if dup_idx.size > 1:
                for j, idx in enumerate(dup_idx):
                    profile[idx, 1] += (dup_idx.size - j) * 0.01
        return profile
    
    if bool_sw:
        profil_thermo_sw = fix_profile(profil_thermo_sw, data_us62)
        # Ordenar de mayor a menor altitud
        profil_thermo_sw = profil_thermo_sw[np.argsort(-profil_thermo_sw[:, 0])]
    if bool_lw:
        profil_thermo_lw = fix_profile(profil_thermo_lw, data_us62)
        profil_thermo_lw = profil_thermo_lw[np.argsort(-profil_thermo_lw[:, 0])]
    
    # 12. Escritura de los archivos finales con encabezado modificado
    header_file = os.path.join(path,  'resources','atm_X_ppmv.dat')
    with open(header_file, 'r') as f:
        header_lines = f.readlines()
    
    # Para SW: se asigna 31 muestras
    if bool_sw:
        header_lines_sw = header_lines.copy()
        if len(header_lines_sw) >= 4:
            header_lines_sw[3] = "31\n"
        sw_output_path = os.path.join(path,  str_type, 'SW',
                                      f'SW_{date.strftime("%Y%m%d%H")}', 'config',
                                      f'atm_{date.strftime("%Y%m%d%H")}_sw_ppmv.dat')
        with open(sw_output_path, 'w') as f:
            f.writelines(header_lines_sw)
            for row in profil_thermo_sw:
                line = "\t".join([f"{val:1.8e}" for val in row]) + "\n"
                f.write(line)
    
    # Para LW: se asigna 40 muestras
    if bool_lw:
        header_lines_lw = header_lines.copy()
        if len(header_lines_lw) >= 4:
            header_lines_lw[3] = "40\n"
        lw_output_path = os.path.join(path,  str_type, 'LW',
                                      f'LW_{date.strftime("%Y%m%d%H")}', 'config',
                                      f'atm_{date.strftime("%Y%m%d%H")}_lw_ppmv.dat')
        with open(lw_output_path, 'w') as f:
            f.writelines(header_lines_lw)
            for row in profil_thermo_lw:
                line = "\t".join([f"{val:1.8e}" for val in row]) + "\n"
                f.write(line)

# --- Funciones auxiliares ---

def reading_rs(date, path):
    """
    Lee los datos de radiosondeo correspondientes a la fecha.
    
    Se utiliza la función ConvertDate2String para determinar el archivo adecuado.
    Retorna una matriz NumPy con columnas:
      [Altitud (m), Presión (hPa), Temperatura (K), Humedad Relativa (%)]
    """
   
    month_str, day_str, hour_str, adjusted_date = ConvertDate2String(date)
    year_str = str(adjusted_date.year)
    pattern = os.path.join(path, f"DORSBCNT0001_{hour_str}{year_str}{month_str}{day_str}*")
    
        
    print("--------------------")
    print("Este es date que se envía: ",date)
    print("Este es el path que se envía: ",path)
    print("Este es el month que se envía: ",month_str)
    print("Este es el year_str que se envía: ",year_str)
    print("Este es el pattern que se envía: ",pattern)          
    print("--------------------")
    
    rs_files= radiosoun(date)
    # rs_file = rs_files[0]
    # Se asume un archivo delimitado por tabulaciones, saltando la primera línea
    df_rs = pd.read_csv(rs_files, delimiter='\t', skiprows=1, header=None)
    # Se extraen las columnas 2, 3, 4 y 5 (índices 1, 2, 3, 4) y se convierte la temperatura a Kelvin
    data_rs = df_rs.iloc[:, [1, 2, 3, 4]].to_numpy()
    data_rs[:, 2] = data_rs[:, 2] + 273.15
    return data_rs

def ConvertDate2String(date):
    """
    Convierte un objeto datetime a cadenas de dos dígitos para mes y día, y determina la hora a usar:
      - Si la hora es antes de las 8 a.m., se usa '00' (radiosondeo 00 UTC).
      - Si está entre 8 y 20, se usa '12' (radiosondeo 12 UTC).
      - Si es después de 20, se utiliza el radiosondeo del día siguiente a las 00 UTC.
    
    Retorna:
      - month_str (str): Mes en dos dígitos.
      - day_str (str): Día en dos dígitos.
      - hour_str (str): Hora ('00' o '12').
      - adjusted_date (datetime): Fecha ajustada en caso de que la hora sea >20.
    """
    adjusted_date = date
    if date.hour < 8:
        hour_str = '00'
    elif 8 <= date.hour <= 20:
        hour_str = '12'
    else:
        adjusted_date = date + timedelta(days=1)
        hour_str = '00'
    
    day_str = f"{adjusted_date.day:02d}"
    month_str = f"{adjusted_date.month:02d}"
    return month_str, day_str, hour_str, adjusted_date
