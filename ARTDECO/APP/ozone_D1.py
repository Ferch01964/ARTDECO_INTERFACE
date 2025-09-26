import cdsapi

# Solicitar al usuario el año y el mes deseados
print("Por favor, ingrese el año y el mes para la descarga de datos.")
year = input("Año (formato YYYY, ej. 2023): ").strip()
month = input("Mes (formato MM, ej. 08): ").strip()

print (year, month)

# Validar entradas
if not (year.isdigit() and len(year) == 4):
    print("El año ingresado no es válido. Debe tener el formato YYYY.")
    exit()

if not (month.isdigit() and len(month) == 2 and 1 <= int(month) <= 12):
    print("El mes ingresado no es válido. Debe tener el formato MM (01-12).")
    exit()

# Crear cliente de CDSAPI
client = cdsapi.Client()

# Configurar parámetros básicos
dataset = "cams-europe-air-quality-reanalyses"
request = {
    "variable": ["ozone"],
    "model": ["monarch"],
    "level": [
        "0",
        "50",
        "100",
        "250",
        "500",
        "750",
        "1000",
        "2000",
        "3000",
        "5000"
    ],
    "type": ["interim_reanalysis"],
    "year": [year],
    "month": [month]
}

# Intentar realizar la solicitud
try:
    print(f"Solicitando datos para el año {year} y el mes {month}...")
    target = 'download.grib'

    client.retrieve(dataset, request, target)
    print("¡Descarga completada!")
except Exception as e:
    print("Error al realizar la solicitud:")
    print(e)
