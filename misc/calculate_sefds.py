import numpy as np
import matplotlib.pyplot as plt


def calculate_sefd_array(sefd_list):

    """
    Takes in a list of antenna SEFDs and calculate the SEFD of the entire array
    """

    N = len(sefd_list)
    sum_inv_sefd_pairs = 0.0

    for i in range(N-1):
        for j in range(i+1,N):
            sum_inv_sefd_pairs = sum_inv_sefd_pairs + (sefd_list[i]*sefd_list[j]) ** -1
    
    sefd_array = sum_inv_sefd_pairs ** -0.5
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

def calculate_thermal_noise(sefd_list, efficiency, bandwidth,obs_time,num_pol):

    """
    Calculate the thermal noise at the phasecenter in a naturally weighted image given
    the efficiency of the array, the bandwidth, the total observing time and the number
    of polarisations
    """

    sefd_array = calculate_sefd_array(sefd_list)
    thermal_noise = (1/efficiency)*(sefd_array)*(1/np.sqrt(2*bandwidth*obs_time*num_pol))

    formatted_noise, unit = format_thermal_noise(thermal_noise)

    # Check if formatted_noise is an array and handle accordingly
    if isinstance(formatted_noise, np.ndarray):
            print(f"The rms noise for a naturally weighted image is {formatted_noise[0]:.3f} {unit}")
    else:
        print(f"The rms noise for a naturally weighted image is {formatted_noise:.3f} {unit}")


    return thermal_noise, unit

# sefd_list = [40,560,19,310,700,300,740,3,10] 
# antennas = [JB,WB,EF,ON,MC,TR,NT,AR,GB]


sefd_list = [40,560,19,310,700,300,740,10] 
efficiency= 0.7
bandwidth = 31.25*1e3
## obs_time in list format allows unequal obs_times
obs_time = [216,216,216,216,216,216,216,216,216] # in minutes or hours 
obs_time = np.mean(obs_time)*60 # change to seconds
num_pol = 2
calculate_thermal_noise(sefd_list, efficiency, bandwidth, obs_time, num_pol)


#### RSG12
# sefd_list = [35,420,20,490,740,350,250,360,300,330,67]  
# efficiency= 0.7
# bandwidth = 8e6
# ## obs_time in list format allows unequal obs_times
# obs_time = [20,20,20,20,20,20,20,20,20,20,20] # in minutes or hours 
# obs_time = np.mean(obs_time)*60 # change to hours
# num_pol = 2
# calculate_thermal_noise(sefd_list, efficiency, bandwidth, obs_time, num_pol)
