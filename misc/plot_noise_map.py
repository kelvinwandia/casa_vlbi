import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j1 # bessel func of order 1


def calculate_attenuation(diameters,offset,wavelength,pb_model):
    """
    Calculate the attenuation for a Gaussian and an Airy Disk modelled
    using the Bessel's function of the first kind

    offset: in radians
    """
    attenuation_factors = []

    offset = offset * np.pi /(180*60) # convert to radians

    P = 1.0 # Peak response, assuming normalized to 1

    if len(pb_model) != len(diameters):
        raise ValueError("Length of pb_models must match length of diameters")
    
    first_bessel_zero = 3.8317  # First zero of the Bessel function j1

    for diameter,pb_model in zip(diameters,pb_model):
        if pb_model == 'G':
            factor_gauss = 4 * np.log(2) * diameter**2 * offset**2
            attenuation = P * np.exp(-factor_gauss / wavelength**2)
        elif pb_model == 'B':
            factor_bessel = (np.pi / wavelength)*diameter*np.sin(offset)
            if factor_bessel == 0:
                # Handle phasecenter where offset is zero
                attenuation = P 
            elif factor_bessel > first_bessel_zero:
                # Handle the zeros of the bessel function -- goes negative after the first dark; which breaks the code
                # the Bessel function will produce multiple darks giving wromg attenuation values for
                # different darks e.g, the positive darks will give positive values of attenuation
                # set all values outside the first dark to zero to avoid this
                attenuation = 0.0000001
            else:
                attenuation = P*(2*j1(factor_bessel)/factor_bessel)
        else:
            print("Attenuations due to primary beams not calculated")
        attenuation_factors.append(attenuation)
    # print(attenuation_factors)
    attenuation_factors = np.sqrt(attenuation_factors) # sqrt is important
    # for attenuation in attenuation_factors:
    #     print(f"{attenuation:.6f}")
    return np.array(attenuation_factors)

def calculate_sefd_array_scaled(sefd_list, attenuation_factors):

    """
    Takes in a list of antenna SEFDs and calculate the SEFD of the entire array and
    then scales the arrays SEFD using the primary beam
    """

    N = len(sefd_list)
    sum_inv_sefd_pairs = 0.0


    for i in range(N-1):
        for j in range(i+1,N):
            sefd_i_scaled =  sefd_list[i] / attenuation_factors[i]
            sefd_j_scaled =  sefd_list[j] / attenuation_factors[j]

            sum_inv_sefd_pairs = sum_inv_sefd_pairs + (sefd_i_scaled*sefd_j_scaled) ** -1

    sefd_array = sum_inv_sefd_pairs ** -0.5
    # print(sum_inv_sefd_pairs**-0.5)
    # print(f"The SEFD of the array is : {sefd_array:.2f} Jy")

    return sefd_array


def format_thermal_noise(thermal_noise):
    """
    Automatically format the thermal noise based on its magnitude.
    
    Returns:
    - Formatted thermal noise and the appropriate unit ('mJy' or 'µJy')
    """
    max_noise = np.max(thermal_noise)
    print(max_noise)
    if max_noise < 1e-6:
        # Convert to micro Jy
        thermal_noise *= 1e6
        unit = 'µJy'
    elif max_noise < 1e-3:
        # Convert to micro Jy
        thermal_noise *= 1e3
        unit = 'mJy'
    else:
        unit = 'Jy'

    return thermal_noise, unit

def calculate_primary_beam(diameter, wavelength):
    return 1.22 * wavelength / diameter

def calculate_thermal_noise(sefd_list,efficiency, bandwidth,obs_time,num_pol,attenuation_factors):

    """ 
    Calculate the thermal noise at the phasecenter in a naturally weighted image given
    the efficiency of the array, the bandwidth, the total observing time and the number
    of polarisations
    """

    sefd_array = calculate_sefd_array_scaled(sefd_list,attenuation_factors)
    thermal_noise = (1/efficiency)*(sefd_array)*(1/np.sqrt(2*bandwidth*obs_time*num_pol))

    thermal_noise = thermal_noise*1e6

    # formatted_noise, unit = format_thermal_noise(thermal_noise)

    # print(f"The rms noise for a naturally weighted image is {thermal_noise:.3f} µJy ")

    return thermal_noise


def create_noise_map(sefd_list,efficiency, bandwidth,obs_time,num_pol,diameters,offset,wavelength,pb_model, grid_size):

    # Create a 2D grid of radial offsets from the phase center
    max_offset_rad = offset * np.pi / (180 * 60)  # Convert arcminutes to radians
    x = np.linspace(-max_offset_rad, max_offset_rad, grid_size)
    y = np.linspace(-max_offset_rad, max_offset_rad, grid_size)
    X, Y = np.meshgrid(x, y)
    radial_offsets = np.sqrt(X**2 + Y**2)  # Radial distance from the center
    
    # Initialize an empty noise map
    noise_map = np.zeros_like(radial_offsets)
    
    # Calculate noise at each point in the grid
    for i in range(grid_size):
        for j in range(grid_size):
            offset_rad = radial_offsets[i, j]
            ## You have called the function calculate attenuation here
            noise_map[i, j] =calculate_thermal_noise(sefd_list,efficiency, bandwidth,obs_time,num_pol,
                            calculate_attenuation(diameters,offset_rad,wavelength,pb_model))
            
    
    # Calculate noise at the phase center
    # The scaling factor will be 1 at the phase_centre
    ones_list = [1.0 for _ in range(len(sefd_list))]
    noise_center =calculate_thermal_noise(sefd_list,efficiency, bandwidth,obs_time,num_pol,ones_list)
    # print(noise_center)
    noise_map = noise_map*1e6
    print(f'Min Noise: {noise_map.min()}, Max Noise: {noise_map.max()}')

    return X, Y, noise_map, noise_center


def main():

    sefd_list = [40,560,19,310,700,300,740,3,10]  
    diameters = [67,25,76,25,32,32,32,213,70.2]  # Diameters of the antennas in meters
    antennas = ['JB','WB','EF','ON','MC','TR','NT','AR','GB']
    pb_model = ['G','G','G','G','G','G','G','B','B']
    efficiency= 0.7
    bandwidth = 128e6
    ## obs_time in list format allows unequal obs_times
    obs_time = [3.6,3.6,3.6,3.6,3.6,3.6,3.6,3.6,3.6] # in minutes or hours 
    obs_time = np.mean(obs_time)*3600 # change to hours
    num_pol = 2
    wavelength = 0.18  # Observing wavelength in meters (e.g., 21 cm for 1.4 GHz)
    grid_size = 100  # Size of the grid for the noise map
    offset = 3.1  # Maximum radial offset in arcminutes




    attenuation_factors = calculate_attenuation(diameters,offset,wavelength,pb_model)

    calculate_thermal_noise(sefd_list,efficiency,bandwidth,obs_time,num_pol,attenuation_factors)


    X, Y, noise_map,noise_center =create_noise_map(sefd_list,efficiency, bandwidth,obs_time,num_pol,diameters,offset,wavelength,pb_model, grid_size)


    # Define the colorbar limits based on the noise center
    noise_center_value = noise_center * 1e6  # Convert noise center to µJy
    vmin = min(noise_map.min(), noise_center_value)
    vmax = max(noise_map.max(), noise_center_value)
    print(vmin)

    # Plotting the noise map
    plt.figure(figsize=(10, 8))
    cmap = plt.get_cmap('inferno')
    cmap.set_under('gray')  # Set color for values below the minimum

    # Display the noise map with adjusted colorbar limits
    cmap = plt.imshow(noise_map, extent=[X.min() * 180 / np.pi * 60, X.max() * 180 / np.pi * 60, Y.min() * 180 / np.pi * 60, Y.max() * 180 / np.pi * 60],
                    origin='lower', cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
    plt.colorbar(label='Thermal Noise (µJy)')






    # Calculate the primary beam for the largest antenna (smallest beamwidth)
    largest_antenna_diameter = max(diameters)
    primary_beam_rad = calculate_primary_beam(largest_antenna_diameter, wavelength)
    primary_beam_arcmin = primary_beam_rad * 180 / np.pi * 60  # Convert primary beam from radians to arcminutes

    # Overlay a circle representing the primary beam of the largest antenna
    circle = plt.Circle((0, 0), primary_beam_arcmin / 2, color='white', fill=False, linestyle='--', label=f'Arecibo PB ({largest_antenna_diameter}m)')
    plt.gca().add_artist(circle)

    plt.xlabel('Offset (arcminutes)')
    plt.ylabel('Offset (arcminutes)')
    plt.legend()
    plt.show()
    plt.savefig('rms_noise_map.pdf',dpi=300)


if __name__=="__main__":
    main()