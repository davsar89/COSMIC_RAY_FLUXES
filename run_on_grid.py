import csv
import subprocess
import numpy as np
import matplotlib.pyplot as plt

def run_executable_with_parameters(ID_nb, energy_thres, altitude):
    # Run the executable with parameters and capture the output
    command = ['./electron_fluxes', str(ID_nb), str(energy_thres), str(altitude)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    # Check for errors
    if process.returncode != 0:
        print("Error:", stderr.decode())
        return None

    # Decode and split the output lines
    output_lines = stdout.decode().split('\n')

    # Parse the output lines to extract required information
    energy_threshold = None
    altitude = None
    energy_integrated_flux = None

    for line in output_lines:
        if 'Energy Threshold:' in line:
            energy_threshold = float(line.split(':')[1].split()[0])
        elif 'Altitude:' in line:
            altitude = float(line.split(':')[1].split()[0])
        elif 'Energy integrated flux :' in line:
            energy_integrated_flux = float(line.split(':')[1].split()[0])

    return energy_threshold, altitude, energy_integrated_flux

###########################################################################################
if __name__ == "__main__":

    output_file = 'results.csv'

    ID_nb = 31 # (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)

    ID_list = []

    eners_thresholds_MeV = [0.3, 1] # MeV
    altitudes_km = np.arange(10, 20 + 0.2, 0.2)

    energy_threshold_list = []
    altitude_list = []
    energy_integrated_flux_list = []

    for energy_thres_MeV in eners_thresholds_MeV:
        for altitude_km in altitudes_km:

            #energy_thres_MeV = 0.3  # Example value for the first parameter
            #altitude_km = 10.0  # Example value for the second parameter
            _, _, energy_integrated_flux = run_executable_with_parameters(ID_nb, energy_thres_MeV, altitude_km)

            if energy_integrated_flux is not None:
                print("Energy Integrated Flux:", energy_integrated_flux, "cm^-2 s^-1")

            ID_list.append(ID_nb)
            energy_threshold_list.append(energy_thres_MeV)
            altitude_list.append(altitude_km)
            energy_integrated_flux_list.append(energy_integrated_flux)

    
    energy_threshold_list = np.round(np.array(energy_threshold_list),2)
    altitude_list = np.round(np.array(altitude_list),2)
    energy_integrated_flux_list = np.round(np.array(energy_integrated_flux_list),4)

    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Particle ID (0:neut, 1-28:ions, 29-30:mu+-, 31:e-, 32:e+, 33:ph)', 'Energy Threshold (MeV)', 'Altitude (kilometers)', 'Energy Integrated Flux (cm^-2 s^-1)'])
        for ii in range(len(energy_threshold_list)):
            writer.writerow([ID_list[ii],energy_threshold_list[ii], altitude_list[ii], energy_integrated_flux_list[ii]])


    # Plot energy_integrated_flux versus altitude
    print(len(altitude_list))
    nn = round(len(altitude_list)/2)
    plt.plot(altitude_list[0:nn], energy_integrated_flux_list[0:nn], label='Min energy thres = 0.3 MeV')
    plt.xlabel('Altitude (kilometers)')
    plt.ylabel('Angular and Energy Integrated Flux (cm^-2 s^-1)')
    plt.title('Electrons')
    plt.grid(True)
    plt.plot(altitude_list[nn+1:], energy_integrated_flux_list[nn+1:], label='Min energy thres = 1.0 MeV')
    plt.legend()
    plt.savefig('flux_vs_altitude_electron.png')
    plt.show()
    
