import requests
from requests.auth import HTTPBasicAuth
from datetime import datetime, date as date_cls
import os

def radiosoun(date_input: str | date_cls | datetime | None = None):

    # Credenciales y URL base
    username = 'fernando.pincay'
    password = 'Mimadre196455!'
    base_url = 'https://drive.tsc.upc.edu/index.php/login'
    file_to_return = "maybe you should check the download of radiosoundings"
    if date_input is None:                       # modo interactivo
        date_input = input("Fecha (YYYY-MM-DD): ").strip()
    elif isinstance(date_input, (date_cls, datetime)):
        date_input = date_input.strftime("%Y-%m-%d")
    else:                                        # viene como str
        try:
            datetime.strptime(date_input, "%Y-%m-%d")
        except ValueError as exc:
            raise ValueError(
                f"'{date_input}' no cumple formato YYYY-MM-DD"
            ) from exc
    
    try:
        # Convertir la fecha ingresada
        date_obj = datetime.strptime(date_input, "%Y-%m-%d")
        year = date_obj.year
        
        # Determinar el formato del nombre del archivo
        if year < 2013:
            file_name = date_obj.strftime("%y%m%d")  # Formato yymmdd para años antes de 2013
        else:
            file_name = date_obj.strftime("%Y%m%d")  # Formato yyyymmdd para años desde 2013
        #https://drive.tsc.upc.edu/index.php/apps/files/?dir=/Lidar/upclidar_Data/Radiosoundings/2025&fileid=91892648
        # Construir la ruta del archivo
        file_path = f'/Lidar/upclidar_Data/Radiosoundings/{year}/DORSBCNT0001_00{file_name}'
        download_url = f'https://drive.tsc.upc.edu/remote.php/webdav{file_path}'
        print(f"{download_url}")
        # Descargar el archivo
        response = requests.get(download_url, auth=HTTPBasicAuth(username, password))
        if response.status_code == 200:
            file_path = os.path.join("INPUTS/resources/Radiosoundings", f"DORSBCNT0001_00{file_name}")

            with open(file_path, 'wb') as file:
                file.write(response.content)
            print(f"Archivo {file_name} descargado exitosamente.")
        else:
            print(f"No se pudo descargar el archivo. Código de estado: {response.status_code}")
    except ValueError:
        print("Formato de fecha inválido. Por favor, use el formato YYYY-MM-DD.")

    return file_path