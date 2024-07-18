import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt

def generate_random_coordinates(ra_center, dec_center, min_separation_arcmin, num_points, seed=None):

    if seed is not None:
        np.random.seed(seed)

    # Convert minimum separation to degrees
    min_separation_deg = min_separation_arcmin / 60.0

    # Convert the radius to degrees
    max_radius_deg = 1 / 60.0

    # Convert the center coordinates to a SkyCoord object
    center_coord = SkyCoord(ra=ra_center, dec=dec_center, unit=(u.deg, u.deg), frame='icrs')

    # List to hold the generated coordinates
    coordinates = []

    # Generate the coordinates ensuring sufficient separation and within the radius
    while len(coordinates) < num_points:
        # Generate random points within a circle of radius 1 arcmin
        r = max_radius_deg * np.sqrt(np.random.uniform(0, 1))
        theta = np.random.uniform(0, 2 * np.pi)
        delta_ra_deg = r * np.cos(theta)
        delta_dec_deg = r * np.sin(theta)

        # Create the new coordinate by applying the offsets
        new_coord = SkyCoord(ra=center_coord.ra.deg + delta_ra_deg,
                             dec=center_coord.dec.deg + delta_dec_deg,
                             unit=u.deg,
                             frame='icrs')

        # Check the separation with all existing points
        sufficient_separation = True
        for coord in coordinates:
            sep = new_coord.separation(coord).deg
            if sep < min_separation_deg:
                sufficient_separation = False
                break

        # If the point has sufficient separation, add it to the list
        if sufficient_separation:
            coordinates.append(new_coord)

    return coordinates

def plot_coordinates(coordinates, ra_center_deg, dec_center_deg):
    # Extract RA and Dec from SkyCoord objects
    ra_degrees = [coord.ra.deg for coord in coordinates]
    dec_degrees = [coord.dec.deg for coord in coordinates]

    # Create the plot
    plt.figure(figsize=(8, 8))
        # Plot the circle with 1 arcmin radius
    circle = plt.Circle((ra_center_deg, dec_center_deg), 1 / 60.0, color='green', fill=False, label='1 arcmin radius')
    plt.gca().add_patch(circle)
    plt.scatter(ra_degrees, dec_degrees, color='blue', marker='o', label='Random Points')
    plt.scatter(ra_center_deg, dec_center_deg, color='red', marker='x', label='Center')



    plt.xlabel('RA (degrees)')
    plt.ylabel('Dec (degrees)')
    plt.title('Random Coordinates with Minimum Separation')
    plt.legend()
    plt.grid(True)
    plt.gca().invert_xaxis()  # RA increases to the left
    plt.tight_layout()
    plt.savefig('coords.pdf')
    # plt.show()

# Example usage
ra_center = 322.49304  # RA in degrees
dec_center = 12.16700  # Dec in degrees
min_separation_arcmin = 0.1  # Minimum separation in arcminutes
num_points = 6
seed = 42

coordinates = generate_random_coordinates(ra_center, dec_center, min_separation_arcmin, num_points, seed)

# Print the generated coordinates
for coord in coordinates:
    print(coord.to_string('hmsdms'))

# Plot the coordinates with the center and the circle
plot_coordinates(coordinates, ra_center, dec_center)
