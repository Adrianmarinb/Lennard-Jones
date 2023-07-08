# -*- coding: utf-8 -*-
# ================================================================================
# ANIMACION SISTEMA SOLAR
#
# Genera una animación a partir de un fichero de datos con las posiciones
# de los planetas en diferentes instantes de tiempo.
# 
# El fichero debe estructurarse de la siguiente forma:
# 
#   t_1
#   x1_1, y1_1
#   x2_1, y2_1
#   x3_1, y3_1
#   (...)
#   xN_1, yN_1
#   
#   t_2
#   x1_2, y1_2
#   x2_2, y2_2
#   x3_2, y3_2
#   (...)
#   xN_2, yN_2
#
#   t_3
#   x1_3, y1_3
#   x2_3, y2_3
#   x3_3, y3_3
#   (...)
#   xN_3, yN_3
#   
#   (...)
#
# donde xi_j es la componente x del planeta i-ésimo en el instante de
# tiempo j-ésimo, e yi_j lo mismo en la componente y. El programa asume que
# el nº de planetas es siempre el mismo.
# ¡OJO! Los datos están separados por comas.
# 
# Si solo se especifica un instante de tiempo, se genera una imagen en pdf
# en lugar de una animación
#
# Se puede configurar la animación cambiando el valor de las variables
# de la sección "Parámetros"
#
# ================================================================================

# Importa los módulos necesarios
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle, Rectangle
import numpy as np
from matplotlib import cm

# Parámetros
# ========================================
file_in = "lennard-jones-output.txt"  # Nombre del fichero de datos
file_out = "particles"  # Nombre del fichero de salida (sin extensión)

# Límites de los ejes X e Y
x_min = -1
x_max = 5
y_min = -1
y_max = 6

interval = 20  # Tiempo entre fotogramas en milisegundos
show_trail = False  # Muestra la "estela" del planeta
trail_width = 1  # Ancho de la estela
save_to_file = False  # False: muestra la animación por pantalla,
# True: la guarda en un fichero
dpi = 150  # Calidad del vídeo de salida (dots per inch)

# Radio del planeta, en las mismas unidades que la posición
# Puede ser un número (el radio de todos los planetas) o una lista con
# el radio de cada uno
planet_radius = 0.1
# planet_radius = [0.5, 0.7, 1.1]

# Lectura del fichero de datos
# ========================================
# Lee el fichero a una cadena de texto
with open(file_in, "r") as f:
    data_str = f.read()

# Inicializa la lista con los datos de cada fotograma.
# frames_data[j] contiene los datos del fotograma j-ésimo
frames_data = list()

# Itera sobre los bloques de texto separados por líneas vacías
# (cada bloque corresponde a un instante de tiempo)
for frame_data_str in data_str.split("\n\n"):
    # Inicializa la lista con la posición de cada planeta
    frame_data = list()

    # Divide las líneas del bloque en dos partes:
    # - La primera línea es el tiempo (t)
    # - Las siguientes líneas son las posiciones de los planetas
    lines = frame_data_str.strip().split("\n")
    t_line = lines[0].strip()  # Get the first line (time line)
    if t_line:
        t = float(t_line)  # Pasamos a float
        frame_data.append(t)  # Añadimos el tiempo

        # Iteramos sobre el resto de líneas 
        for planet_pos_str in lines[1:]:
            planet_pos = np.fromstring(planet_pos_str, sep=",")
            # Si la línea no está vacía
            if planet_pos.size > 0:
                frame_data.append(np.fromstring(planet_pos_str, sep=","))

        # Add the data of this frame to the list
        frames_data.append(frame_data)

# El número de planetas es el número de líneas en cada bloque - 1
# Lo calculamos del primer bloque
nplanets = len(frames_data[0]) - 1

# Generamos los colores que queramos visualizar
colors = cm.viridis(np.linspace(0, 1, nplanets))

# Creación de la animación/gráfico
# ========================================
# Crea los objetos figure y axis
fig, ax = plt.subplots(figsize=(6, 6.54))

# Define el rango de los ejes
ax.axis("equal")  # Misma escala para ejes X e Y
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

# Definimos x_label y y_label
plt.xlabel('x(σ)', fontsize = 15)
plt.ylabel('y(σ)', fontsize = 15)

ax.tick_params(axis="both", which="major", labelsize=15, pad=10)

box = Rectangle((0, 0), 4, 4, edgecolor='green', linewidth=4, fill=False)
ax.add_patch(box)

# Si solo se ha dado un radio para todos los planetas, conviértelo a una
# lista con todos los elementos iguales
if not hasattr(planet_radius, "__iter__"):
    planet_radius = planet_radius * np.ones(nplanets)
# En caso contrario, comprueba que el nº de radios coincide con el
# nº de planetas y devuelve error en caso contrario
else:
    if not nplanets == len(planet_radius):
        raise ValueError(
            "El número de radios especificados no coincide con el número "
            "de planetas")

# Representa el primer fotograma
# Pinta un punto en la posición de cada paneta y guarda el objeto asociado
# al punto en una lista
planet_points = list()
planet_trails = list()
for planet_pos, radius, color in zip(frames_data[0][1:], planet_radius, colors):
    x, y = planet_pos
    planet_point = Circle((x, y), radius)
    planet_point.set_facecolor(color)
    ax.add_artist(planet_point)
    planet_points.append(planet_point)

    # Inicializa las estelas (si especificado en los parámetros)
    if show_trail:
        planet_trail, = ax.plot(
            x, y, "-", linewidth=trail_width,
            color=planet_points[-1].get_facecolor())
        planet_trails.append(planet_trail)

# Creamos el texto que muestra el tiempo
time_text = ax.text(0.07, 0.90, '', transform=ax.transAxes, fontsize=15)

# Función que actualiza la posición de los planetas en la animación 
def update(frame, frames_data, planet_points, planet_trails, show_trail):
    # Valor del tiempo para el frame
    t = frames_data[frame][0]

    # Actualiza la posición del correspondiente a cada planeta
    for j_planet, planet_pos in enumerate(frames_data[frame][1:]):
        x, y = planet_pos
        planet_points[j_planet].center = (x, y)

        if show_trail:
            xs_old, ys_old = planet_trails[j_planet].get_data()
            xs_new = np.append(xs_old, x)
            ys_new = np.append(ys_old, y)

            planet_trails[j_planet].set_data(xs_new, ys_new)

    # Actualizamos el valor del tiempo
    time_text.set_text('Tiempo: {:.3f} s'.format(t))

    return planet_points + planet_trails + [time_text]

def init_anim(): 
    if show_trail:
        for j_planet in range(nplanets):
            planet_trails[j_planet].set_data(list(), list())

    # Actualizamos el valor del tiempo
    t = frames_data[0][0]
    time_text.set_text('Tiempo: {:.3f} s'.format(t))

    return planet_points + planet_trails + [time_text]

# Calcula el nº de frames
nframes = len(frames_data)

# Si hay más de un instante de tiempo, genera la animación
if nframes > 1:
    # Info sobre FuncAnimation: https://matplotlib.org/stable/api/animation_api.html
    animation = FuncAnimation(
        fig, update, init_func=init_anim,
        fargs=(frames_data, planet_points, planet_trails, show_trail),
        frames=len(frames_data), blit=True, interval=interval)

    # Muestra por pantalla o guarda según parámetros
    if save_to_file:
        animation.save("{}.mp4".format(file_out), dpi=dpi)
    else:
        plt.show()
# En caso contrario, muestra o guarda una imagen
else:
    # Muestra por pantalla o guarda según parámetros
    if save_to_file:
        fig.savefig("{}.pdf".format(file_out))
    else:
        plt.show()

plt.subplots_adjust(
top=0.88,
bottom=0.13,
left=0.145,
right=0.91,
hspace=0.2,
wspace=0.2
)
