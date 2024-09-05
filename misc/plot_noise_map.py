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
    print(attenuation_factors)
    attenuation_factors = np.sqrt(attenuation_factors) # sqrt is important
    for attenuation in attenuation_factors:
        print(f"{attenuation:.6f}")
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
    print(sum_inv_sefd_pairs**-0.5)
    print(f"The SEFD of the array is : {sefd_array} Jy")

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

    print(f"The rms noise for a naturally weighted image is {thermal_noise:.3f} µJy ")

    return thermal_noise


def make_noise_map(sefd_list,efficiency,bandwidth,obs_time,num_pol,attenuation_factors):
    """
    Make an rms noise map
    """

    attenuation_factors = calculate_attenuation(diameters,offset,wavelength,pb_model)

    calculate_thermal_noise(sefd_list,efficiency,bandwidth,obs_time,num_pol,attenuation_factors)


# Example values (adjust these as per your data)
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
offset = 3.0  # Maximum radial offset in arcminutes



attenuation_factors = calculate_attenuation(diameters,offset,wavelength,pb_model)

calculate_thermal_noise(sefd_list,efficiency,bandwidth,obs_time,num_pol,attenuation_factors)

