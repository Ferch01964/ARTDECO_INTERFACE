import os
import numpy as np  # Asegurarse de importar numpy, ya que se usa en el código
import tkinter as tk
from datetime import datetime, timedelta
from tkinter import ttk
from tkcalendar import Calendar
from PIL import Image, ImageTk 
from datetime import date

# Importar los módulos propios (asegúrate de que estas rutas o el paquete estén correctamente configurados)
from atmosphere_artdeco import atmosphere_artdeco
from main_kdis_disort_artdeco import main_kdis_disort_artdeco
from surface_artdeco import surface_artdeco
from geometrics_artdeco import geometrics_artdeco
from prop_rad_cirrus_artdeco import prop_rad_cirrus_artdeco
from prop_rad_aer_artdeco import prop_rad_aer_artdeco
from particles_artdeco import particles_artdeco

# GUI Logic
def show_selections():
    # Step 1: Get checkbox states
    bool_cld = bool_cld_check.get()
    bool_aer = bool_aer_check.get()
    bool_sw = bool_sw_check.get()
    bool_lw = bool_lw_check.get()

    # Step 2: Convert calendar date
    selected_date = calendar.get_date()
    v_date = datetime.strptime(selected_date, "%m/%d/%y").replace(hour=12)

    # Step 3: Update GUI with result
    # result = f"Checkboxes: {[bool_cld, bool_aer, bool_sw, bool_lw]}\nSelected Date: {v_date}"
    # result_label.config(text=result)

    # Ejecución de funciones

    # La función main_kdis_disort_artdeco funciona correctamente
    main_kdis_disort_artdeco(path, min_wvl, max_wvl, v_date, bool_cld, bool_aer, bool_sw, bool_lw)

    # atmosphere_artdeco "necesita retoques"; 
    # review implementation 
    #  puede requerir ajustes en los parámetros o en la forma de ejecución.
    atmosphere_artdeco(path, v_date, bool_cld, bool_aer, bool_sw, bool_lw, lat_place, lon_place)

    # La función surface_artdeco funciona correctamente
    surface_artdeco(path, min_wvl, max_wvl, v_date, bool_cld, bool_aer, bool_sw, bool_lw, lat_place, lon_place)

    # Se ejecuta la función geometrics_artdeco
    geometrics_artdeco(path, v_date, bool_cld, bool_aer, bool_sw, bool_lw, lat_place, lon_place)

    # Ejemplo de variables para la nube:
    # v_cbh y v_cth representan las alturas base y tope del cirrus [km]
    v_cbh = [1.0]   # Ejemplo: 1 km
    v_cth = [2.0]   # Ejemplo: 2 km
    k = 0  # Índice de la escena (en este ejemplo, el primero)

    # v_h_lidar: perfil vertical de alturas (ejemplo de 0 a 20 km con 100 puntos)
    y_nrb = np.linspace(0, 20, 100)
    # t_par_z: perfil de temperatura en la nube (por ejemplo, de 250K a 270K)
    t_par_z = [np.linspace(250, 270, 100)]
    # a_par_z: perfil de coeficientes de extinción (por ejemplo, de 0.01 a 0.05 m^-1)
    a_par_z = [np.linspace(0.01, 0.05, 100)]

    # Inicialización de las variables de salida de prop_rad_cirrus_artdeco
    m_ext_cld, m_asy_cld, m_ssa_cld, m_sc_cld, m_abs_cld, v_iwc_cld = (None,)*6

    if bool_cld:
        m_ext_cld, m_asy_cld, m_ssa_cld, m_sc_cld, m_abs_cld, v_iwc_cld = prop_rad_cirrus_artdeco(
            path,
            (v_cbh[k], v_cth[k]),      # prop_geo: (cbh, cth)
            y_nrb,                     # Perfil de alturas del MPL
            t_par_z[0] + 273.15,       # Temperatura en Kelvin (se suma 273.15 para conversión, si es necesario)
            a_par_z[0] * 1e-3          # Extinción en m^-1 (multiplicada por 1e-3)
        )
    if bool_aer:
        # Llamada a la función de aerosoles
        data_aeronet_prof, data_aeronet_interpn = prop_rad_aer_artdeco(path, v_date)

    # Construcción de la "cell array" de propiedades de la nube
    # (en MATLAB: prop_cld{1,1} = prop_geo; prop_cld{1,2} = m_ext_cld; etc.)
    prop_cld = [None]*4
    prop_cld[0] = m_ext_cld
    prop_cld[1] = m_asy_cld
    prop_cld[2] = m_ssa_cld
    prop_cld[3] = m_sc_cld

    # Construcción de la "cell array" de aerosoles
    # (en MATLAB: prop_aer{1,1} = data_aeronet_prof; prop_aer{1,2} = data_aeronet_interpn)
    prop_aer = [None]*2
    prop_aer[0] = data_aeronet_prof
    prop_aer[1] = data_aeronet_interpn

    particles_artdeco(path,v_date, bool_cld, bool_aer, bool_sw, bool_lw, prop_cld, prop_aer)

def resize_bg(event):
    global bg_image, bg_photo

    new_width = event.width
    new_height = event.height
    x_center = new_width // 2

    # Resize and update background
    resized_bg = bg_image.resize((new_width, new_height), Image.Resampling.LANCZOS)
    bg_photo = ImageTk.PhotoImage(resized_bg)
    canvas.itemconfig(bg_image_id, image=bg_photo)

    # Reposition logo (top-right with margin)
    canvas.coords(logo_image_id, new_width - 10, 10)

    # Stack checkboxes
    start_y = 30
    spacing = 25
    for i, cb_id in enumerate(checkbox_ids):
        canvas.coords(cb_id, x_center, start_y + i * spacing)

    # Reposition calendar under checkboxes
    total_cb_height = start_y + len(checkbox_ids) * spacing
    canvas.coords(calendar_id, x_center, total_cb_height + 20)

    # Use actual calendar bounding box to place submit button
    bbox = canvas.bbox(calendar_id)
    if bbox:
        _, _, _, calendar_bottom = bbox
        canvas.coords(submit_id, x_center, calendar_bottom + 20)

# Variables
path = os.path.abspath("INPUTS")
min_wvl = 0.350
max_wvl = 0.408
lat_place = 41.38
lon_place = 2.11
BACKGROUND_PATH = "006_Tibidabo.jpg"
LOGO_PATH = "EETAC.png"

# GUI Setup
root = tk.Tk()
root.title("ARTDECO")
root.minsize(500, 500)

# Canvas base
canvas = tk.Canvas(root, width=500, height=500)
canvas.pack(fill="both", expand=True)

# Background image
bg_image = Image.open(BACKGROUND_PATH).resize((500, 500), Image.Resampling.LANCZOS)
logo_image = Image.open(LOGO_PATH).resize((150, 150), Image.Resampling.LANCZOS)

bg_photo = ImageTk.PhotoImage(bg_image)
logo_photo = ImageTk.PhotoImage(logo_image)
# canvas.create_image(0, 0, anchor='nw', image=bg_photo)
bg_image_id = canvas.create_image(0, 0, anchor='nw', image=bg_photo)
logo_image_id = canvas.create_image(0, 0, anchor='ne', image=logo_photo)

root.bg_photo = bg_photo  # prevent garbage collection
root.logo_photo = logo_photo

# Checkboxes
bool_cld_check = tk.BooleanVar()
bool_aer_check = tk.BooleanVar()
bool_sw_check = tk.BooleanVar()
bool_lw_check = tk.BooleanVar()

checkboxes = [
    tk.Checkbutton(root, text="CLD", variable=bool_cld_check),
    tk.Checkbutton(root, text="AER", variable=bool_aer_check),
    tk.Checkbutton(root, text="SW", variable=bool_sw_check),
    tk.Checkbutton(root, text="LW", variable=bool_lw_check)
]

checkbox_ids = []
for i, cb in enumerate(checkboxes):
    checkbox_ids.append(canvas.create_window(150, 20 + i * 25, anchor='n', window=cb))  # Centered at x=150

# Calendar widget
calendar = Calendar(root, selectmode='day',
                    year=date.today().year,
                    month=date.today().month,
                    day=date.today().day)
calendar_id = canvas.create_window(165, 140, anchor='n', window=calendar)

# # Result label
# result_label = tk.Label(root, text="", wraplength=250)
# canvas.create_window(20, 310, anchor='nw', window=result_label)

# Submit button
submit_btn = tk.Button(root, text="Submit", command=show_selections)
submit_id = canvas.create_window(150, 400, anchor='n', window=submit_btn)

canvas.bind("<Configure>", resize_bg)

# Start app
root.mainloop()