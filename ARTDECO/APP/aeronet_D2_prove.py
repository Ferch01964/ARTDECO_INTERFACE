import requests

def download_with_referer():
    # Solicitar rango de fechas al usuario
    print("Ingrese el período deseado para la descarga de datos:")
    start_date = input("Fecha de inicio (YYYYMMDD): ").strip()
    end_date = input("Fecha de fin (YYYYMMDD): ").strip()

    # Validar las fechas ingresadas
    if not (len(start_date) == 8 and len(end_date) == 8 and start_date.isdigit() and end_date.isdigit()):
        print("Las fechas deben estar en formato YYYYMMDD.")
        return

    # Ubicación fija
    location = "Barcelona"

    # URL del archivo
    base_url = "https://aeronet.gsfc.nasa.gov/zip_files_v3/inv"
    file_name = f"{start_date}_{end_date}_{location}.zip"
    url = f"{base_url}/{file_name}"

    # Imprimir URL para verificar
    print(f"Generando URL: {url}")

    # Cabeceras incluyendo el Referer
    headers = {
        "Referer": "https://aeronet.gsfc.nasa.gov/cgi-bin/webtool_inv_v3",
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36"
    }

    # Realizar la solicitud GET
    print(f"Descargando datos desde: {url}...")
    response = requests.get(url, headers=headers)

    # Verificar si la solicitud fue exitosa
    if response.status_code == 200:
        # Guardar el archivo con el nombre correspondiente
        with open(file_name, "wb") as file:
            file.write(response.content)
        print(f"Archivo descargado correctamente y guardado como {file_name}.")
    elif response.status_code == 404:
        print(f"Error 404: No se encontró el archivo para el período {start_date} a {end_date}.")
    else:
        print(f"Error al descargar el archivo: {response.status_code}")
        print("Verifique que las fechas sean correctas.")

# Ejecutar la función
if __name__ == "__main__":
    download_with_referer()
