import os
import glob
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from datetime import datetime, timedelta

r'''
Variables que usaré para validar la función...las copio directamente de Main_inputs_ARTDECO.m


path = 'C:\Users\ferna\Desktop\ARTDECO - FERNANDO'
min_wvl = 0.350; % Minimum wavelenth in the spectral range (microns)
max_wvl = 0.408; % Maximum wavelenth in the spectral range (microns)
v_date = datenum(2025,2,3,12,0,0);
v_date = datetime(2025, 2, 3, 12, 0, 0)
bool_cld = False
bool_aer = True
bool_sw = True
bool_lw = False
lat_place = 41.38
lon_place = 2.11
'''

def surface_artdeco(path, min_wvl, max_wvl, date, bool_cld, bool_aer, bool_sw, bool_lw, lat_place, lon_place):
    """
    Función que genera los archivos de superficie para ARTDECO en simulaciones SW y LW.
    
    Parámetros:
        path (str): Ruta base de trabajo.
        min_wvl (float): Longitud de onda mínima (en micras) del rango espectral.
        max_wvl (float): Longitud de onda máxima (en micras) del rango espectral.
        date (datetime): Fecha y hora de interés (objeto datetime).
            date --> fecha y hora de simulación en matlab se carga con 
                    v_date = datenum(2025,2,3,12,0,0);  -  % All the dates are scanned - 
                    if hour(v_date(k)) == 0
                        bool_sw = false;
                    end

        bool_cld (bool): Indica si la escena atmosférica tiene nubes.
        bool_aer (bool): Indica si la escena atmosférica tiene aerosoles.
        bool_sw (bool): Si se desea generar el archivo para simulación en SW.
        bool_lw (bool): Si se desea generar el archivo para simulación en LW.
        lat_place (float): Latitud del lugar de la escena.
        lon_place (float): Longitud del lugar de la escena.
    """
    
    #############################################
    # 1. Determinación del tipo de escena (str_type)
    #############################################
    if bool_cld and bool_aer:
        str_type = 'GAC'
    elif bool_cld and not bool_aer:
        str_type = 'CLD'
    elif not bool_cld and bool_aer:
        str_type = 'AER'
    else:
        str_type = 'GAS'
    
    #############################################
    # 2. Lectura del archivo plantilla "surface.py"
    #############################################
    template_path = os.path.join(path,  "templates","surface.py")
    with open(template_path, "r") as file:
        # Se lee todo el contenido y se almacena en una lista de líneas
        template_lines = file.readlines()

    
    #############################################
    # 3. Lectura de datos NetCDF para cálculo de SST y albedo
    #############################################
    #intentamos recuperar los archivos de internet, falta el codigo para recuperar de esta pagina 
    #pagina web: https://urs.earthdata.nasa.gov/users/new
    #user: fercho
    #pw: ElJivarito1993!


    # Buscamos el primer archivo que contenga 'NOAA20' y extensión .nc en el directorio 'path'
    nc_files = glob.glob(os.path.join(path,"resources","NOAA20.nc"))
    if not nc_files:
        raise FileNotFoundError("No se encontró ningún archivo NOAA20.nc en la ruta indicada.")
    
    nc_file = nc_files[0]
    
    # Abrir el archivo NetCDF en modo lectura
    with Dataset(nc_file, 'r') as nc:
        # Imprimir las variables y grupos disponibles para depuración
        #print("Variables en el grupo raíz:", nc.variables.keys())
        #print("Grupos en el archivo:", nc.groups.keys())
        #print("Información de la variable 'Conditions':", nc.variables['Conditions'])
        
        # Acceder al grupo 'Time_and_Position'
        tp_group = nc.groups['Time_and_Position']
        # Leer la longitud y latitud
        lon_noaa = tp_group.variables['instrument_fov_longitude'][:]
        lat_noaa = tp_group.variables['instrument_fov_latitude'][:]
        # Leer el vector de tiempo (se asume que representa días desde 1970-01-01)
        aux_time = nc.variables['time'][:]
        # Convertir cada valor en aux_time a un objeto datetime
        time_noaa = np.array([datetime(1970, 1, 1) + timedelta(days=float(t)) for t in aux_time])
        
        # Acceder al grupo 'TOA_and_Surface_Fluxes' para la emisividad
        toa_group = nc.groups['TOA_and_Surface_Fluxes']
        emissivity = toa_group.variables['longwave_surface_emissivity'][:]
        
        # Acceder al grupo 'Full_Footprint_Area' para la temperatura superficial
        ffa_group = nc.groups['Full_Footprint_Area']
        sst_noaa = ffa_group.variables['surface_skin_temperature'][:]
        print("lon: ",lon_noaa,"lat: ", lat_noaa,"aux_time: ", aux_time,"toa: ", toa_group,"time_noaa: ", time_noaa,"emi: ", emissivity,"ffa: ", ffa_group,"sst: ", sst_noaa)
   

    # Selección de datos cercanos a la ubicación y fecha de interés
    # Se define una máscara para los datos que estén dentro de ±0.5° de la latitud y longitud deseadas
    mask = (
        (lon_noaa < (lon_place + 0.5)) & (lon_noaa > (lon_place - 0.5)) &
        (lat_noaa < (lat_place + 0.5)) & (lat_noaa > (lat_place - 0.5))
    )
    # Se filtran los tiempos para que coincidan en año, mes y día con la fecha objetivo
    time_years = np.array([t.year for t in time_noaa])
    time_months = np.array([t.month for t in time_noaa])
    time_days = np.array([t.day for t in time_noaa])
    mask &= (time_years == date.year) & (time_months == date.month) & (time_days == date.day)
    
    indices = np.where(mask)[0]

    if len(indices) > 1:
        # Si hay más de un registro, se selecciona el que tenga menor diferencia temporal con 'date'
        diffs = np.array([abs((time_noaa[i] - date).total_seconds()) for i in indices])
        idx_selected = indices[np.argmin(diffs)]
        sst = sst_noaa[idx_selected]
        albedo = 1 - emissivity[idx_selected]
    elif len(indices) == 1:
        idx_selected = indices[0]
        sst = sst_noaa[idx_selected]
        albedo = 1 - emissivity[idx_selected]
    else:
        # Si no se encontró coincidencia exacta en día, se busca registros que coincidan en mes y hora
        mask_alt = (
            (lon_noaa < (lon_place + 0.5)) & (lon_noaa > (lon_place - 0.5)) &
            (lat_noaa < (lat_place + 0.5)) & (lat_noaa > (lat_place - 0.5))
        )
        time_hours = np.array([t.hour for t in time_noaa])
        mask_alt &= (time_months == date.month) & (time_hours == date.hour)
        indices_alt = np.where(mask_alt)[0]
        if len(indices_alt) == 0:
            raise ValueError("No se encontraron registros alternativos para la fecha y ubicación especificadas.")
        # Se calcula el promedio (ignorando posibles NaN) de los valores encontrados
        sst = np.nanmean(sst_noaa[indices_alt])
        albedo = 1 - np.nanmean(emissivity[indices_alt])

        print("------------------------------")
        print(sst,albedo)
        print("------------------------------")

    #############################################
    # 4. Creación de archivos de configuración para SW
    #############################################
    if bool_sw:
        # Definición de la ruta de destino para archivos SW
        sw_dir = os.path.join(path,  str_type, "SW", "SW_" + date.strftime("%Y%m%d%H"), "config")
        os.makedirs(sw_dir, exist_ok=True)  # Crear directorio si no existe
        
        # Copia de la plantilla y modificación de líneas específicas:
        # Línea para Surface_nature (índice 18 en Python, equivalente a línea 19 en MATLAB)
        template_lines[18] = 'Surface_nature = "sw_' + date.strftime("%Y%m%d%H") + '"\n'
        # Línea para Surface_temperature (índice 26 en Python)
        template_lines[26] = 'Surface_temperature = ' + format(sst, ".4f") + '\n'
        
        # Escritura del archivo modificado en la ruta destino
        sw_file_path = os.path.join(sw_dir, "surface.py")
        with open(sw_file_path, "w") as f:
            f.writelines(template_lines)
  
    #############################################
    # 5. Creación de archivos de configuración para LW
    #############################################
    if bool_lw:
        # Definición de la ruta de destino para archivos LW
        lw_dir = os.path.join(path,  str_type, "LW", "LW_" + date.strftime("%Y%m%d%H"), "config")
        os.makedirs(lw_dir, exist_ok=True)
        
        # Modificación de la plantilla para LW:
        template_lines[18] = 'Surface_nature = "cste"\n'
        # Línea para Surface_albedo (índice 21 en Python, equivalente a línea 22 en MATLAB)
        template_lines[21] = 'Surface_albedo = ' + format(albedo, ".4f") + '\n'
        template_lines[26] = 'Surface_temperature = ' + format(sst, ".4f") + '\n'
        
        # Escritura del archivo modificado en la ruta destino
        lw_file_path = os.path.join(lw_dir, "surface.py")
        with open(lw_file_path, "w") as f:
            f.writelines(template_lines)
     
    #############################################
    # 6. Creación de archivo con perfil atmosférico (AERONET) para SW
    #############################################
    if bool_sw:
        # Buscar el archivo que contenga 'AERONET_INVERSION_L15' en su nombre
        aeronet_files = glob.glob(os.path.join(path, "resources","*AERONET_INVERSION_L15*"))
        print(path)
        if not aeronet_files:
            raise FileNotFoundError("No se encontró ningún archivo AERONET_INVERSION_L15 en la ruta indicada.")
        aeronet_file = aeronet_files[0]
         
        # Lectura del archivo CSV usando pandas.
        # Se asume que las primeras 6 filas son de encabezado, archivo BARCELONA_AERONET_INVER....se hace skip de las primeras filas de encabezado
        usecols = ["Date(dd:mm:yyyy)", "Time(hh:mm:ss)", "Surface_Albedo[440m]", 
                   "Surface_Albedo[675m]", "Surface_Albedo[870m]", "Surface_Albedo[1020m]"]
        df = pd.read_csv(aeronet_file, delimiter=",", skiprows=6, usecols=usecols)

        # Convertir las columnas de fecha y hora a objetos datetime.
        # Se asume que el formato de fecha es "dd:MM:yyyy" y de hora "HH:mm:ss"
        df["Dateddmmyyyy"] = pd.to_datetime(df["Date(dd:mm:yyyy)"], format="%d:%m:%Y", errors='coerce')
        df["Timehhmmss"] = pd.to_datetime(df["Time(hh:mm:ss)"], format="%H:%M:%S", errors='coerce').dt.time
        # Combinar fecha y hora en una única columna datetime
        df["datetime"] = df.apply(lambda row: datetime.combine(row["Dateddmmyyyy"].date(), row["Timehhmmss"]), axis=1)
        
        # Encontrar el registro cuya fecha-hora sea el más cercano a la fecha de interés
        df["diff"] = df["datetime"].apply(lambda t: abs((t - date).total_seconds()))
        idx = df["diff"].idxmin()
        selected_row = df.loc[idx]
        
        # Definición de los valores de longitudes de onda AERONET (en micras)
        lambda_aeronet = np.array([440, 675, 870, 1020]) * 1e-3  # conversión de nm a micras
        num_albedo = len(lambda_aeronet)
        # Se comprueba si es necesario añadir un valor extra para cubrir el rango espectral
        bool_min = False
        bool_max = False
        if min_wvl < lambda_aeronet[0]:
            bool_min = True
            num_albedo += 1
        if max_wvl > lambda_aeronet[-1]:
            bool_max = True
            num_albedo += 1
        
        # Definir la ruta de salida para el archivo de albedo SW
        alb_file_name = "surfalb_sw_" + date.strftime("%Y%m%d%H") + ".dat"
        alb_file_path = os.path.join(path,  str_type, "SW", "SW_" + date.strftime("%Y%m%d%H"), "config", alb_file_name)
        os.makedirs(os.path.dirname(alb_file_path), exist_ok=True)
        
        # Apertura y escritura del archivo de albedo
        with open(alb_file_path, "w") as f:
            # Escribir cabecera
            f.write("#\n")
            f.write("# Surface albedo :\n")
            f.write("#\n")
            f.write("# Number of sampled values :\n")
            f.write(f"{num_albedo}\n")
            f.write("# wavelength        albedo\n")
            f.write("#  (microns)\n")
            
            # Si es necesario, escribir el primer valor (para min_wvl)
            if bool_min:
                # Se aplica una corrección para fijar la resolución (ejemplo: truncar a 1 decimal)
                wavelength_min = np.fix(min_wvl * 10) / 10
                # Se utiliza el albedo correspondiente al primer canal (440 nm) como referencia
                alb_min = np.nanmean([selected_row["Surface_Albedo[440m]"]])
                f.write(f"{wavelength_min:1.4f} \t {alb_min:1.4f}\n")
            
            # Escribir para cada longitud de onda definida en AERONET
            # Se asume el orden de columnas: 440, 675, 870 y 1020 nm
            albedo_values = [selected_row["Surface_Albedo[440m]"],
                             selected_row["Surface_Albedo[675m]"],
                             selected_row["Surface_Albedo[870m]"],
                             selected_row["Surface_Albedo[1020m]"]]
            for lam, alb in zip(lambda_aeronet, albedo_values):
                f.write(f"{lam:1.4f} \t {np.nanmean([alb]):1.4f}\n")
            
            # Si es necesario, escribir el último valor (para max_wvl)
            if bool_max:
                alb_max = np.nanmean([selected_row["Surface_Albedo[1020m]"]])
                f.write(f"{max_wvl:1.4f} \t {alb_max:1.4f}\n")
r'''
'''