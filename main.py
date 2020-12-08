import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import box, Point, LineString
from shapely import affinity
import math


# This is making sure that the picture will fit
PICTURE_SIDE = 0.4
margin = PICTURE_SIDE * np.sqrt(2) / 2
low = np.ceil(margin * 100)
high = np.floor((1-margin) * 100)
CENTER_OF_PICTURE = (np.random.randint(low=low, high=high, size=2)) / 100


def rotate_about_origin(x, y, origin, deg):
    x, y = np.asarray(x), np.asarray(y)
    ox, oy = origin

    xx = ox + (x-ox) * np.cos(np.deg2rad(deg)) + (y-oy) * np.sin(np.deg2rad(deg))
    yy = oy - (x-ox) * np.sin(np.deg2rad(deg)) + (y-oy) * np.cos(np.deg2rad(deg))
    return xx, yy


coords = np.random.rand(100, 2)
ids = list(range(1, len(coords)+1))

x, y = zip(*coords)
catalogue = pd.DataFrame({'ID': ids, 'X': x, 'Y': y, 'grid': None})
catalogue['grid'] = catalogue['grid'].astype(object)

GRID_SIZE = 0.2
GRID_N = 10
SEARCH_RADIUS = GRID_SIZE


def make_pattern(star, close_stars):

    if not close_stars:
        return np.array([])

    grid = np.zeros((GRID_N, GRID_N))
    closest_star = min(close_stars, key=lambda x: np.sqrt((x[0] - star[0]) ** 2 + (x[1] - star[1]) ** 2))

    angle_to_rotate = math.atan2(closest_star[1] - star[1], closest_star[0] - star[0])

    unzipped_coords = list(zip(*close_stars))
    rotated_x, rotated_y = rotate_about_origin(
        unzipped_coords[0],
        unzipped_coords[1],
        (star[0], star[1]),
        np.rad2deg(angle_to_rotate)
    )

    x_linspace = np.linspace(star[0] - GRID_SIZE, star[0] + GRID_SIZE, num=GRID_N, endpoint=True)
    y_linspace = np.linspace(star[1] - GRID_SIZE, star[1] + GRID_SIZE, num=GRID_N, endpoint=True)[::-1]

    for close_star in list(zip(rotated_x, rotated_y)):
        x_idx = (np.abs(x_linspace - close_star[0])).argmin()
        y_idx = (np.abs(y_linspace - close_star[1])).argmin()
        grid[y_idx, x_idx] = 1

    return np.flatnonzero(grid)


"""""""""""""""""""""""""""""""""""""""""""""""""""
 Make grids/patterns for all stars in the catalogue
"""""""""""""""""""""""""""""""""""""""""""""""""""
for index, row in catalogue.iterrows():

    # This thing below includes current star
    close_stars_df = catalogue.loc[np.sqrt((catalogue['X']-row['X']) ** 2 + (catalogue['Y']-row['Y']) ** 2) <= SEARCH_RADIUS]

    # I'm pretty sure this is making a copy and not modifying the dataframe? I hope so.
    close_stars = list(zip(close_stars_df['X'].tolist(), close_stars_df['Y'].tolist()))

    # Get rid of the actual star we're considering
    close_stars.remove((row['X'], row['Y']))

    pattern = make_pattern((row['X'], row['Y']), close_stars)
    catalogue.loc[index, 'grid'] = [[pattern]]

    if not index % 10:
        print(f"Processed star nr: {index}")


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Draw stars and make a picture frame, then collect stars inside
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
fig, (ax, ax2) = plt.subplots(1, 2,figsize=(16, 8))
ax.set_aspect('equal', adjustable='box')
ax2.set_aspect('equal', adjustable='box')
fig.suptitle('Star tracker', fontsize=20)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

ax.scatter(x, y, marker='*')
for i, id in enumerate(ids):
    ax.annotate(id, (x[i], y[i]))

# Make a box/frame
b = box(
    CENTER_OF_PICTURE[0] - PICTURE_SIDE / 2,    # minx
    CENTER_OF_PICTURE[1] - PICTURE_SIDE / 2,    # miny
    CENTER_OF_PICTURE[0] + PICTURE_SIDE / 2,    # maxx
    CENTER_OF_PICTURE[1] + PICTURE_SIDE / 2     # maxy
)


picture_angle = np.random.randint(0, 90)
rotated_box = affinity.rotate(b, picture_angle, 'center')
box_representative_point = rotated_box.representative_point().coords[:][0]

ax.plot(*rotated_box.exterior.xy)
ax.annotate('picture', box_representative_point)

box_x, box_y = rotated_box.exterior.coords.xy
rotated_box_coords = list(zip(box_x, box_y))
top_line = LineString([rotated_box_coords[1], rotated_box_coords[2]])
ax.plot(*top_line.xy)
# top_line_representative_point = top_line.representative_point().coords[:][0]
# plt.annotate('top', top_line_representative_point)

# Collect points inside
points_in_picture = []

for p in coords:
    if rotated_box.contains(Point(p[0], p[1])):
        points_in_picture.append(p)

# Draw stars in the picture as the tracker sees them, just for show
#fig2, ax2 = plt.subplots(1, figsize=(6, 6))
# plt.xlim(0, 1)
# plt.ylim(0, 1)

cx, cy = zip(*points_in_picture)
rx, ry = rotate_about_origin(cx, cy, CENTER_OF_PICTURE, picture_angle)
ax2.scatter(rx, ry, marker='*')

"""""""""""""""""""""""""""""""""""""""""""""""""""""
 Search through the catalogue for all stars in picture
"""""""""""""""""""""""""""""""""""""""""""""""""""""

candidates = []
for star in points_in_picture:

    most_plausible_hits = []
    close_stars = [p for p in points_in_picture if (p != star).all() and np.sqrt((p[0]-star[0]) ** 2 + (p[1]-star[1]) ** 2) <= SEARCH_RADIUS]
    pattern = make_pattern(star, close_stars)
    max_intersect = 2   # Min amount of close stars we want to confirm

    # If we found less than 2 neighbours, it immediately disqualifies this star
    if len(pattern) < max_intersect:
        candidates.append([])
        continue

    # Search through the database and find possible candidates
    # (I know itterrows is lame but let's chill out, there's only like a 100 entries here)
    for index, row in catalogue.iterrows():
        row_patt = row['grid'][0]

        if len(row_patt) == len(pattern) and (row_patt == pattern).all():

            print('jackpot for', star, "match:", row['ID'])
            most_plausible_hits = [row['ID'], 'O']
            break

        intersect = len(np.intersect1d(row_patt, pattern))
        if intersect > max_intersect:
            max_intersect = intersect
            most_plausible_hits = [row['ID']]

        elif intersect == max_intersect:
            most_plausible_hits.append(row['ID'])

    candidates.append(most_plausible_hits)

"""""""""""""""""
 Display and label
"""""""""""""""""
all_cases = 0
true_cases = 0

for i, cand in enumerate(candidates):
    annotation = cand[0] if cand else ''

    if len(cand) == 1:                                             # We've got one most likely candidate, color it green
        ax2.annotate(annotation, (rx[i], ry[i]), color='green')
        true_id = catalogue.loc[(catalogue['X'] == cx[i]) & (catalogue['Y'] == cy[i])]['ID'].values[0]
        if true_id == annotation:
            true_cases += 1
        all_cases += 1
    elif 'O' in cand:                                              # We found a jackpot for this one, color it blue
        ax2.annotate(cand[0], (rx[i], ry[i]), color='blue')
        true_cases += 1
        all_cases += 1
    else:
        ax2.annotate(cand, (rx[i], ry[i]))                         # We have few, equally likely candidates for this
        all_cases += 1

if all_cases > 0:
    print("accuracy: ", true_cases/all_cases, " %")

else:
    print("duppa")

plt.show()
