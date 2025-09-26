

import os
import math
from datetime import datetime

def geometrics_artdeco(path, date, bool_cld, bool_aer, bool_sw, bool_lw, lat_place, lon_place):
    """
    Función que crea los archivos geometrics para ARTDECO en SW/LW.
    
    Parámetros:
        path (str): Directorio base.
        date (datetime): Fecha y hora de la simulación (se asume un objeto datetime).
        bool_cld (bool): Indica si la escena atmosférica tiene nubes.
        bool_aer (bool): Indica si la escena atmosférica tiene aerosoles.
        bool_sw (bool): Indica si se debe crear el archivo para simulación en SW.
        bool_lw (bool): Indica si se debe crear el archivo para simulación en LW.
        lat_place (float): Latitud del lugar de la escena.
        lon_place (float): Longitud del lugar de la escena.
    """

    # ----------------------------------------------------------
    # Selección del tipo de escena atmosférica
    # ----------------------------------------------------------
    if bool_cld and bool_aer:
        str_type = 'GAC'
    elif bool_cld and not bool_aer:
        str_type = 'CLD'
    elif not bool_cld and bool_aer:
        str_type = 'AER'
    else:
        str_type = 'GAS'
    
    # ----------------------------------------------------------
    # Lectura del archivo template "geometrics.py"
    # ----------------------------------------------------------
    template_file = os.path.join(path,  'templates','geometrics.py')
    with open(template_file, 'r') as file:
        tmp = file.readlines()
    
    # ----------------------------------------------------------
    # Cálculo del ángulo cenital solar (SZA)
    # Se verifica si la hora UTC es 00. En ese caso se asigna 90, 
    # de lo contrario se calcula usando la función anidada.
    # Aunque en Matlab se guarda la variable 'sza', luego se escribe en 
    # el template utilizando la llamada a sza_computation con lat y lon.
    # ----------------------------------------------------------
    if date.hour == 0:
        sza = 90
    else:
        sza = sza_computation(date, lat_place, lon_place)
    
    # Modificación de la línea 33 del template (en Python, índice 32)
    # Se usa el valor de sza_computation con latitud y longitud y se formatea a 3 decimales.
    tmp[32] = 'Solar_zenith_angle_list = [%.3f]\n' % sza_computation(date, lat_place, lon_place)

    # ----------------------------------------------------------
    # Creación y escritura del archivo para SW (cortas)
    # ----------------------------------------------------------
    if bool_sw:
        # Construir el nombre del directorio para SW
        dir_name = os.path.join(path,  str_type, 'SW', 'SW_' + date.strftime('%Y%m%d%H'), 'config')
        # Verificar si el directorio existe; si no, crearlo
        if not os.path.exists(dir_name):
            os.makedirs(dir_name, exist_ok=True)
        # Ruta completa del archivo a escribir
        sw_file_path = os.path.join(dir_name, 'geometrics.py')
        # Escritura del archivo con el contenido modificado del template
        with open(sw_file_path, 'w') as f:
            for line in tmp:
                f.write(line)
    
    # ----------------------------------------------------------
    # Creación y escritura del archivo para LW (largas)
    # ----------------------------------------------------------
    if bool_lw:
        # Construir el nombre del directorio para LW
        dir_name = os.path.join(path,  str_type, 'LW', 'LW_' + date.strftime('%Y%m%d%H'), 'config')
        # Verificar si el directorio existe; si no, crearlo
        if not os.path.exists(dir_name):
            os.makedirs(dir_name, exist_ok=True)
        # Ruta completa del archivo a escribir
        lw_file_path = os.path.join(dir_name, 'geometrics.py')
        # Escritura del archivo con el contenido modificado del template
        with open(lw_file_path, 'w') as f:
            for line in tmp:
                f.write(line)

    
    # ----------------------------------------------------------
    # Función anidada para el cálculo del SZA
    # Se implementa exactamente la misma lógica que en Matlab.
    # ----------------------------------------------------------
def sza_computation(date, lat_place, lon_place):
        # Variables
        lat = lat_place
        lon = lon_place
        # Obtener el día del año
        jday = date.timetuple().tm_yday
        # Hora UTC
        utc_time = date.hour
        # Constantes para el cálculo
        b1 = 0.006918
        b2 = 0.399912
        b3 = 0.070257
        b4 = 0.006758
        b5 = 0.000907
        b6 = 0.002697
        b7 = 0.001480
        a1 = 0.000075
        a2 = 0.001868
        a3 = 0.032077
        a4 = 0.014615
        a5 = 0.040849
        
        # Desfase horario en función de la longitud
        # Se asume que cada grado equivale a 4 minutos; se convierte a horas.
        hlong = lon * 4 / 60.0  
        
        # Conversión de la latitud a radianes
        latr = lat * math.pi / 180.0
        # Ángulo para calcular la declinación solar
        tet = 2.0 * math.pi * jday / 365.0
        # Cálculo de la declinación solar
        dec_day = (b1 - b2 * math.cos(tet) + b3 * math.sin(tet) -
                   b4 * math.cos(2.0 * tet) + b5 * math.sin(2.0 * tet) -
                   b6 * math.cos(3.0 * tet) + b7 * math.sin(3.0 * tet))
        # Ecuación del tiempo
        time_eq = (a1 + a2 * math.cos(tet) - a3 * math.sin(tet) -
                   a4 * math.cos(2.0 * tet) - a5 * math.sin(2.0 * tet)) * 12.0 * 60.0 / math.pi
        # Tiempo solar verdadero (TSV)
        tsv = utc_time + hlong + time_eq / 60.0
        # Hora angular en radianes
        ahr = (tsv - 12.0) * 15.0 * math.pi / 180.0
        # Cálculo de la elevación solar
        h = math.asin(math.sin(latr) * math.sin(dec_day) + math.cos(latr) * math.cos(dec_day) * math.cos(ahr))
        # Cálculo del ángulo cenital solar
        teta_S = (math.pi / 2.0 - max(h, 0.0)) + 0.00000001
        # Mu0 se calcula pero no se utiliza en la salida
        mu0 = math.cos(teta_S)
        # Conversión del ángulo a grados
        sza = teta_S * 180.0 / math.pi
        return sza
