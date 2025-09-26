import os, shutil
import fileinput as f
from ArtdecoConfigErrors import *
from ArtdecoTools import *

class ArtdecoConfig:
    ## Initializes a configuration with default parameters.
    def __init__(self):

        self.Lib_directory = os.getcwd() + "/lib/"
        self.Input_directory = os.getcwd() + "/input/"
        self.Output_directory = os.getcwd() + "/out/"
        self.Output_subdirectory = "default_out/"
        self.Output_format = "ascii"
        self.Prefix_of_input_files = "default_prefix"

        self.Verbose = False
        self.Warnings = False
        self.Log_file = "default_log.txt"

        self.Mode = "mono"
        self.Wavelength_unit = "microm"
        self.Wavelength_list = [0.550]
        self.Depolarization_value = -1
        self.Instrumental_filters_list = []
        self.K_distributions_parameterization = "default_kdis"
        self.Gas_list = []
        self.Gas_continuum_list = []
        self.Wavelength_start = 0.2
        self.Wavelength_end = 50.0

        self.RTE_solver = "none"
        self.Compute_radiance = True
        self.Legendre_and_truncation = True
        self.Rayleigh_scattering = True
        self.Stokes_components = 1
        self.Thermal_only = False
        self.Thermal_in_situ = False
        self.Surface_parameters = ["lambert", "cste", 0.0]
        self.Surface_temperature = -1
        self.Wind_speed = 5.0
        self.Pigment_concentration = 0.0
        self.Ocean_salinity = 34.3
        self.Shadowing = False
        self.Incident_spectrum = "kurudz_medium"
        self.FBeam = []
        self.Solar_configuration_mode = "angle"
        self.Solar_configuration_list = [[60.0]]
        self.View_zenith_angle_list = [15.0]
        self.View_azimuth_angle_list = [0.0]
    
        self.With_particles = False
        self.Particle_option = []
        self.Truncation_method = "none"
        self.Theta_min = 14.0
        self.Theta_max = 15.0
        self.Theta_cut = 0.0
        self.Fit_all = False
        self.Number_of_betal = -1
        self.TMS_correction = False

        self.With_atmosphere = False
        self.Atmosphere_in_ppmv = False
        self.Atmospheric_profile = "noatm"
        self.Atmospheric_profile_original = "noatm"
        self.Uniform_gas_list = []

        self.Print_downward_radiance = False
        self.Print_betal = False
        self.Print_recomposed_phase_matrix = False

        self.keywords = []

    ## Keywords update.
    def AddKeyword(self, keyword):
        if not(keyword in self.keywords):
            self.keywords.append(keyword)

    def UpdateKeywords(self):
        if self.Verbose:
            self.AddKeyword("verbose")
        if not(self.Warnings):
            self.AddKeyword("nowarn")
        if not(self.Rayleigh_scattering):
            self.AddKeyword("no_rayleigh")
        if self.RTE_solver == "none":
            if self.Legendre_and_truncation:
                self.AddKeyword("betal_only")
            else:
                self.AddKeyword("opt_only")
        if self.Print_betal:
            self.AddKeyword("print_betal")
        if self.Print_recomposed_phase_matrix:
            self.AddKeyword("print_recomp")
        if not(self.Compute_radiance):
            self.AddKeyword("flux_only")
        if self.Thermal_only:
            self.AddKeyword("thermal_only")
        if self.Print_downward_radiance:
            self.AddKeyword("print_down_rad")
        if self.Atmosphere_in_ppmv:
            self.AddKeyword("atm_in_ppmv")
        
    ## RTE solver configuration and specific parameters files.
    def RadiativeModelConfig(self):
        if self.RTE_solver == "disort":
            import disort
            self.Compute_radiance, self.Print_downward_radiance, \
                self.Thermal_only, self.Thermal_in_situ = \
                disort.BasicParameters()
            self.N_stream, self.Convergence_criterion = \
                disort.AdvancedParameters()

        elif self.RTE_solver == "doad":
            import doad
            self.Print_downward_radiance = doad.BasicParameters()
            self.N_stream, self.Adding_epsilon_accuracy, \
                self.Maximum_fourier_terms, self.Minimum_fourier_terms = \
                doad.AdvancedParameters()

        elif self.RTE_solver == "mcrad1d":
            import mcrad1d
            self.N_photons, self.Maximum_interactions, \
                self.Variance_reduction_methods, self.Epsilon_ddis, \
                self.Mc_vrm_n_firstcp, self.Mc_vrm_n_sccp, self.Mc_vrm_n_lecp, \
                self.Critical_weight_splitting, \
                self.Maximum_number_secondary_photons, \
                self.Critical_weight_russian_roulette, \
                self.Minimum_probability_photon \
                = mcrad1d.AdvancedParameters()

    def WriteRadiativeModelConfigFile(self):
        if self.RTE_solver == "disort":
            file_out = open(self.Input_directory + "/" + self.Prefix_of_input_files + \
                                "_disort_spec.dat", "w")
            file_out.write(str(int(self.N_stream)) + "\n#\n")
            file_out.write(str(self.Convergence_criterion) + "\n")
            file_out.close()

        elif self.RTE_solver == "doad":
            file_out = open(self.Input_directory + self.Prefix_of_input_files + \
                                "_doad_spec.dat", "w")
            file_out.write(str(self.Adding_epsilon_accuracy) + "\n#\n")
            file_out.write(str(int(self.N_stream)) + "\n#\n")
            file_out.write(str(int(self.Maximum_fourier_terms)) + "\n#\n")
            file_out.write(str(int(self.Minimum_fourier_terms)) + "\n")
            file_out.close()

        elif self.RTE_solver == "mcrad1d":
            file_out = open(self.Input_directory + self.Prefix_of_input_files + \
                                "_mcrad1d_spec.dat", "w")
            file_out.write(str(int(self.N_photons)) + "\n#\n")
            file_out.write(str(self.Maximum_interactions) + "\n#\n")
            file_out.write("0" * (not(self.Variance_reduction_methods)) + \
                               "1" * self.Variance_reduction_methods + "\n#\n")
            file_out.write(str(self.Epsilon_ddis) + "\n#\n")
            file_out.write(str(self.Mc_vrm_n_firstcp) + " " + str(self.Mc_vrm_n_sccp) + \
                               " " + str(self.Mc_vrm_n_lecp) + "\n#\n")
            file_out.write(str(self.Critical_weight_splitting) + "\n#\n")
            file_out.write(str(self.Maximum_number_secondary_photons) + "\n#\n")
            file_out.write(str(self.Critical_weight_russian_roulette) + "\n#\n")
            file_out.write(str(self.Minimum_probability_photon) + "\n")
            file_out.close()

    def SymlinkToRadiativeModelConfigFile(self):
        if self.RTE_solver == "disort":
            if os.path.lexists(self.Input_directory + "/od_spec.dat"):
                os.remove(self.Input_directory + "/od_spec.dat")
            os.symlink(self.Input_directory + "/" + self.Prefix_of_input_files + \
                           "_disort_spec.dat", self.Input_directory + "/od_spec.dat")

        elif self.RTE_solver == "doad":
            if os.path.lexists(self.Input_directory + "doad_spec.dat"):
                os.remove(self.Input_directory + "doad_spec.dat")
            os.symlink(self.Input_directory + self.Prefix_of_input_files + \
                           "_doad_spec.dat", self.Input_directory + "doad_spec.dat")

        elif self.RTE_solver == "mcrad1d":
            if os.path.lexists(self.Input_directory + "mcrad1d_spec.dat"):
                os.remove(self.Input_directory + "mcrad1d_spec.dat")
            os.symlink(self.Input_directory + self.Prefix_of_input_files + \
                           "_mcrad1d_spec.dat", self.Input_directory + \
                           "mcrad1d_spec.dat")
        
    ## Mode (monochromatic or k-distribution).
    def ModeConfig(self):
        if self.Mode == "mono":
            import mono
            self.Wavelength_unit, self.Wavelength_list, \
                self.Depolarization_value = mono.BasicParameters()
        elif self.Mode == "kdis":
            import kdis
            self.K_distributions_parameterization, self.Gas_list, \
                self.Gas_continuum_list, self.Instrumental_filters_list, \
                self.Wavelength_unit, self.Wavelength_start, self.Wavelength_end, \
                self.Depolarization_value = kdis.BasicParameters()

    ## In kdis mode, copies the definition file of kdis, and writes 
    ## the continuum flags.
    def WriteKdisFiles(self):
        if self.Mode == "kdis":
            source_file = "./lib/kdis/reference/kdis_" + \
                self.K_distributions_parameterization + "_def.dat"
            target_file = "./lib/kdis/kdis_" + \
                self.K_distributions_parameterization + "_def.dat"
            stream_in = open(source_file, "r")
            lines_in = stream_in.readlines()
            stream_in.close()
            stream_out = open(target_file, "w")
            for line in lines_in:
                items_list = line.split(" ")
                if "continuum_flag" in items_list:
                    value = "0.000"
                    for gc in self.Gas_continuum_list:
                        if gc in items_list:
                            value = "1.000"
                    items_list = [value if item == "continuum_flag" else item for item in items_list]
                if "continuum_flag\n" in items_list:
                    value = "0.000\n"
                    for gc in self.Gas_continuum_list:
                        if gc in items_list:
                            value = "1.000\n"
                    items_list = [value if item == "continuum_flag\n" else item for item in items_list]
                if "continuum_flag\t" in items_list:
                    value = "0.000\t"
                    for gc in self.Gas_continuum_list:
                        if gc in items_list:
                            value = "1.000\t"
                    items_list = [value if item == "continuum_flag\t" else item for item in items_list]
                stream_out.write(" ".join(items_list))
            stream_out.close()

            for gas in self.Gas_list:
                source_file = "./lib/kdis/reference/kdis_" + \
                    self.K_distributions_parameterization + "_" + gas + ".dat"
                target_file =  "./lib/kdis/kdis_" + \
                    self.K_distributions_parameterization + "_" + gas + ".dat"
                shutil.copyfile(source_file, target_file)
     
    ## Atmosphere options.
    def AtmosphereConfig(self):
        import atmosphere
        self.Atmosphere_in_ppmv, self.Atmospheric_profile, \
            self.Uniform_gas_list = atmosphere.BasicParameters()

    ## Atmosphere profile.
    def WriteAtmosphericProfile(self):
        ## Takes the atmospheric profile file specified in atmosphere.py
        ## and creates a new one with only the gas in Gas_list of kdis.py
        ## that are not in Uniform_gas_list.
        if self.With_atmosphere:
            file_in = open("./lib/atm/reference/atm_" + self.Atmospheric_profile + ".dat", "r")
            if self.Atmosphere_in_ppmv:
                file_out = open("./lib/atm/atm_current_atm_ppmv.dat", "w")
            else:
                file_out = open("./lib/atm/atm_current_atm.dat", "w")
            
            lines_in = file_in.readlines()
            file_in.close()
            index = 0
            current_line = lines_in[index]
            while current_line[0] == "#":
                index += 1
                current_line = lines_in[index]
            Nz = int(current_line)
            index += 1
            current_line = lines_in[index]
            while current_line[0] == "#":
                index += 1
                current_line = lines_in[index]
            Ng = int(current_line)
            index += 1
            current_line = lines_in[index]
            while current_line[0] == "#":
                index += 1
                current_line = lines_in[index]
            
            col_list = []
            atm_gas_list = []
            index_gas = 0
            uniform_gas_list = [item[0] for item in self.Uniform_gas_list]
            while current_line[0] != "#":
                gas = current_line.strip()
                if (gas.lower() in self.Gas_list) \
                        and (not(gas.lower() in uniform_gas_list)):
                    col_list.append(index_gas)
                    atm_gas_list.append(gas)
                index_gas += 1
                index += 1
                current_line = lines_in[index]

            file_out.write("#\n")
            file_out.write(str(Nz) + "\n#\n")
            file_out.write(str(len(col_list)) + "\n#\n")
            for gas in atm_gas_list:
                file_out.write(gas + "\n")
            if self.Atmosphere_in_ppmv:
                information_line = "# HGT [km] PRE [mb] TEM [K] "
                for gas in atm_gas_list:
                    information_line += gas + " [ppmv] "
                information_line += "\n"
            else:
                information_line = "# HGT [km] PRE [mb] TEM [K] air [cm-3] "
                for gas in atm_gas_list:
                    information_line += gas + " [cm-3] "
                information_line += "\n"
            file_out.write(information_line)
            index += 1
            for iz in range(Nz):
                current_line = lines_in[index]
                values_list = clean_list(current_line.split(None))
                if self.Atmosphere_in_ppmv:
                    file_out.write(values_list[0] + " " + values_list[1] + " " + \
                                       values_list[2] + " ")
                    for ig in range(len(atm_gas_list)):
                        file_out.write(values_list[col_list[ig] + 3] + " ")
                    file_out.write("\n")
                else:
                    file_out.write(values_list[0] + " " + values_list[1] + " " + \
                                       values_list[2] + " " + values_list[3] + " ")
                    for ig in range(len(atm_gas_list)):
                        file_out.write(values_list[col_list[ig] + 4] + " ")
                    file_out.write("\n")
                index += 1
                      
            file_out.close()
            self.Atmospheric_profile_original = self.Atmospheric_profile
            self.Atmospheric_profile = "current_atm"
            if self.Atmosphere_in_ppmv:
                self.Atmospheric_profile += "_ppmv"
        else:
            file_in = "./lib/atm/reference/atm_noatm.dat"
            file_out = "./lib/atm/atm_noatm.dat"
            shutil.copyfile(file_in, file_out)
        
                
    def WriteVdistFiles(self):
        import os, shutil
        if not(os.path.isdir(self.Input_directory + "/vdist/")):
            os.mkdir(self.Input_directory + "/vdist/")
        for particle in self.Particle_option:
            if particle[5] == "user":
                shutil.copyfile(self.Lib_directory + "/vdist/vdist_" + \
                                    particle[0] + ".dat", \
                                    self.Input_directory + \
                                    "/vdist/vdist_" + particle[0] + ".dat")
        
                
    def WriteReferenceInputFiles(self):
        import os, shutil
        if not(os.path.lexists(self.Input_directory + "/artdeco_constants.dat")):
            shutil.copyfile(self.Lib_directory + \
                                "/constants/artdeco_constants.dat", \
                                self.Input_directory + "/artdeco_constants.dat")

    def WriteOceanConfigFile(self):
        if self.Surface_parameters[1] == "ocean":
            file_out = open(self.Input_directory + self.Prefix_of_input_files + \
                                "_ocean_spec.dat", "w")
            file_out.write(str(self.Wind_speed) + "\n#\n" + \
                               str(self.Ocean_salinity) + "\n#\n" + \
                               str(self.Pigment_concentration) + "\n#\n" + \
                               ("yes" * self.Shadowing + \
                                    "no" * (not(self.Shadowing))) + "\n#\n")
            file_out.close()
            if os.path.lexists(self.Input_directory + "ocean_spec.dat"):
                os.remove(self.Input_directory + "ocean_spec.dat")
            os.symlink(self.Input_directory + self.Prefix_of_input_files + \
                           "_ocean_spec.dat", \
                           self.Input_directory + "ocean_spec.dat")
            

    def WriteUniformGasesConcentrationFile(self):
        
        ## Writes the file with the constant uniform distribution 
        ## of concentration.
        file_out = open(self.Input_directory + self.Prefix_of_input_files + \
                            "_uniform_gases_conc.dat", "w")
        file_out.write(str(len(self.Uniform_gas_list)) + "\n#")
        for gas in self.Uniform_gas_list:
            file_out.write("\n" + gas[0] + " " + str(gas[1]))
        file_out.write("\n#\n")
        file_out.close()
        if os.path.lexists(self.Input_directory + "uniform_gases_conc.dat"):
            os.remove(self.Input_directory + "uniform_gases_conc.dat")
        os.symlink(self.Input_directory + self.Prefix_of_input_files + \
                       "_uniform_gases_conc.dat", \
                        self.Input_directory + "uniform_gases_conc.dat")

    ## Particle options.
    def ParticleConfig(self):
        try:
            import particles
        except ImportError:
            mesg = "particles.py : The program can not find particles.py.\n"
        except:
            mesg = "particles.py : Error with at least one argument.\n"
            try:
                import traceback
                tback = traceback.format_exc()
                tlines = tback.split("\n")
                mesg += "The line with an error is :\n   "
                mesg += (tlines[-3]).strip()
                mesg += "\n"
                mesg += (tlines[-2]).strip()
                mesg += "\n"
            except:
                print("no traceback module")
            RaiseError(self.Log_file, mesg)
        self.Particle_option, self.Truncation_method, \
            self.Number_of_betal, self.TMS_correction = particles.BasicParameters()
        ## Phase matrix truncation options.
        if self.Truncation_method == "potter":
            self.Theta_min, self.Theta_max = particles.PotterParameters()
        elif self.Truncation_method == "dfit":
            self.Theta_cut, self.Fit_all = particles.DfitParameters()

    ## Truncation options.
    def WriteTruncationMethodConfigFile(self):
        if self.Truncation_method == "potter":
            file_out = open(self.Input_directory + \
                                self.Prefix_of_input_files + \
                            "_potter_spec.dat", "w")
            file_out.write(str(self.Theta_min) + "\n#\n" + \
                               str(self.Theta_max) + "\n")
            file_out.close()
        elif self.Truncation_method == "dfit":
            file_out = open(self.Input_directory + \
                                self.Prefix_of_input_files + \
                            "_dfit_spec.dat", "w")
            file_out.write(str(self.Theta_cut) + "\n#\n" + \
                               "1" * self.Fit_all + \
                               "0" * (not(self.Fit_all)) + "\n")
            file_out.close()

    def SymlinkToTruncationMethodConfigFile(self):
        if self.Truncation_method == "potter":
            if os.path.lexists(self.Input_directory + "potter_spec.dat"):
                os.remove(self.Input_directory + "potter_spec.dat")
            os.symlink(self.Input_directory + \
                           self.Prefix_of_input_files + \
                           "_potter_spec.dat", \
                           self.Input_directory + "potter_spec.dat")
        elif self.Truncation_method == "dfit":
            if os.path.lexists(self.Input_directory + "dfit_spec.dat"):
                os.remove(self.Input_directory + "dfit_spec.dat")
            os.symlink(self.Input_directory + \
                           self.Prefix_of_input_files + \
                           "_dfit_spec.dat", \
                           self.Input_directory + "dfit_spec.dat")
            
            
    ## Surface configuration.
    def SurfaceConfig(self):
        import surface
        self.Surface_parameters, \
            self.Surface_temperature = surface.BasicParameters()
        if self.Surface_parameters[1] == "ocean":
            self.Wind_speed, self.Ocean_salinity, self.Pigment_concentration, \
                self.Shadowing = surface.OceanParameters()

    ## Geometric configuration
    def GeometricsConfig(self):
        import geometrics
        self.Incident_spectrum, self.FBeam, \
            self.Solar_configuration_mode, self.Solar_configuration_list, \
            self.View_zenith_angle_list, self.View_azimuth_angle_list = \
            geometrics.BasicParameters(self.Mode)


    ## Writes the main configuration file.
    def WriteMainConfigurationFile(self):
        file_out = open(self.Input_directory + self.Prefix_of_input_files + \
                            "_artdeco_in.dat", "w")
        ## Writes keywords.
        if len(self.keywords) == 0:
            file_out.write("none")
        else:
            for kw in self.keywords:
                file_out.write(kw + " ")
        file_out.write("\n#\n")

        ## kdis Mode, Gas list, filters and spectral domain.
        if self.Mode == "kdis":
            file_out.write(self.Mode + " " + self.K_distributions_parameterization + "\n#\n")
            if len(self.Gas_list) == 0:
                file_out.write("none")
            else:
                file_out.write(self.Gas_list[0])
                for i in range(len(self.Gas_list) - 1):
                    file_out.write("," + self.Gas_list[i + 1])
            file_out.write("\n#\n")
            if len(self.Instrumental_filters_list) == 0:
                file_out.write("none\n")
            else:
                for instrumental_filter in self.Instrumental_filters_list:
                    file_out.write(instrumental_filter + "\n")
            file_out.write("#\n")
            file_out.write(str(self.Wavelength_start) + " " + str(self.Wavelength_end))
            file_out.write("\n#\n")
        ## Mono mode.
        elif self.Mode == "mono":
            file_out.write(self.Mode + "\n#\n")
            file_out.write("none\n#\nnone\n#\n")
            for wavelength in self.Wavelength_list:
                file_out.write(str(wavelength) + " ")
            file_out.write("\n#\n")

        ## Depolarization.
        file_out.write(str(self.Depolarization_value) + " ")
        file_out.write("\n#\n")
        
        ## Particles.
        if self.With_particles:
            for particle in self.Particle_option:
                file_out.write(particle[0] + " " + \
                                   ("yes" * particle[2] + \
                                        "no" * (not(particle[2]))) + " " + \
                                   ("yes" * particle[3] + \
                                        "no" * (not(particle[3]))) + " " + \
                                   ("yes" * particle[1] + \
                                        "no" * (not(particle[1]))) + " " + \
                                   str(particle[4]) + " " + \
                                   particle[5])
                if len(particle) > 6:
                    file_out.write(" " + str(particle[6]))
                if len(particle) > 7:
                    file_out.write(" " + str(particle[7]))
                file_out.write("\n")
            file_out.write("#\n" + self.Truncation_method + "\n#\n" + \
                               str(self.Number_of_betal) + "\n#\n")
        else:
            file_out.write("none\n#\nnone\n#\n-1\n#\n")

        ## If any, RTE solver, and related options.
        #if self.RTE_solver != "none":
        file_out.write(self.RTE_solver + "\n#\n")
        file_out.write("yes" * self.TMS_correction + \
                           "no" * (not(self.TMS_correction)) + "\n#\n")
        file_out.write("yes" * self.Thermal_in_situ + \
                           "no" * (not(self.Thermal_in_situ)) + "\n#\n")
        file_out.write(str(self.Stokes_components) + "\n#\n")
        if self.Atmospheric_profile == "noatm":
            file_out.write(self.Atmospheric_profile + "\n#\n")
        elif self.Atmosphere_in_ppmv:
            file_out.write("current_atm_ppmv\n#\n")
        else:
            file_out.write("current_atm\n#\n")
        
        for i in range(len(self.Surface_parameters)):
            file_out.write(str(self.Surface_parameters[i]) + " ")
        file_out.write("\n#\n")
        file_out.write(str(self.Surface_temperature) + "\n#\n")
        
        file_out.write(self.Incident_spectrum + "\n#\n")
        if len(self.FBeam) == 0:
            file_out.write("no\n#\n")
        else:
            file_out.write("yes " + str(self.FBeam[0]) + "\n#\n")

        if len(self.Solar_configuration_list[0]) == 1:
            for i in range(len(self.Solar_configuration_list)):
                file_out.write(str(self.Solar_configuration_list[i][0]) + " ")
            file_out.write("\n#\nnone\n#\n")
        else:
            file_out.write("-1\n#\n")
            for i in range(len(self.Solar_configuration_list)):
                file_out.write(str(self.Solar_configuration_list[i][0]) + " " + \
                                   str(self.Solar_configuration_list[i][1]) + " " + \
                                   str(self.Solar_configuration_list[i][2]) + " " + \
                                   str(self.Solar_configuration_list[i][3]) + "\n#\n")
                
        for i in range(len(self.View_zenith_angle_list)):
            file_out.write(str(self.View_zenith_angle_list[i]) + " ")
        file_out.write("\n#\n")
        for i in range(len(self.View_azimuth_angle_list)):
            file_out.write(str(self.View_azimuth_angle_list[i]) + " ")
        file_out.write("\n")
            
        file_out.close()

    ## Writes specific format outputs.
    def WriteOutput(self, format):
        if format == "netcdf" or format == "hdf5" or format == "all":
            import numpy, glob
            output_dir = "./out/" + self.Output_subdirectory + "/"
            
            #########################
            ## Read radiances.
            for radiance_file in glob.glob(output_dir + "*_radiance.dat"):
                direction = (radiance_file.split("/")[-1]).split("_")[0]

                ## Reads ascii file and erase comment lines.
                file_in = open(radiance_file, "r")
                raw_lines_in = file_in.readlines()
                file_in.close()
                lines_in = clean_comment_lines(raw_lines_in)
                del raw_lines_in

                ## Number of Stokes parameters.
                Ns = int(clean_line(lines_in[2])[0])
                
                ## dimensions.
                Nw = int(clean_line(lines_in[3])[0])
                Nsza = int(clean_line(lines_in[1])[0])
                Nvza = int(clean_line(lines_in[1])[1])
                Nvaa = int(clean_line(lines_in[1])[2])

                ## dimension variables.
                Wv = numpy.zeros([Nw], dtype = "f8")
                Sza = numpy.zeros([Nsza], dtype = "f8")
                Vza = numpy.zeros([Nvza], dtype = "f8")
                Vaa = numpy.zeros([Nvaa], dtype = "f8")

                ## Upward radiance array.
                up_rad = numpy.zeros([Ns, Nw, Nsza, Nvza, Nvaa], dtype = "f8")

                for isza in range(Nsza):
                    Sza[isza] = float(clean_line(lines_in[isza * Nvza * Nvaa + 6])[0])
                for ivza in range(Nvza):
                    Vza[ivza] = float(clean_line(lines_in[ivza * Nvaa + 6])[1])
                for ivaa in range(Nvaa):
                    Vaa[ivaa] = float(clean_line(lines_in[ivaa + 6])[2])
                for iw in range(Nw):
                    Wv[iw] = float(clean_line(lines_in[iw * (Nsza * Nvza * Nvaa + 2) + 4])[0])           
                                
                for iw in range(Nw):
                    for isza in range(Nsza):
                        for ivza in range(Nvza):
                            for ivaa in range(Nvaa):
                                for iparam in range(Ns):
                                    up_rad[iparam, iw, isza, ivza, ivaa] = \
                                        float(clean_line(lines_in[iw * (Nsza * Nvza * Nvaa + 2) + \
                                                            isza * Nvza * Nvaa + \
                                                            ivza * Nvaa + \
                                                            ivaa + 6])[4 + iparam])
                           
                if format == "netcdf" or format == "all":
                    import netCDF4, time
                    ThisDataset = netCDF4.Dataset(output_dir + direction + "_radiance.nc", "w", format = "NETCDF4")

                    ThisDataset.created = "Created " + time.ctime(time.time())
                    ThisDataset.Source = "ARTDECO V1.1"
                    ThisDataset.URL = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisDataset.case_name = self.Output_subdirectory
                    ThisDataset.mode = self.Mode
                    if self.Mode == "kdis":
                        ThisDataset.k_distribution = self.K_distributions_parameterization
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisDataset.gas_list = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisDataset.gas_continuum_list = gas_continuum_list_str
                    ThisDataset.radiative_transfer_model = self.RTE_solver
                    ThisDataset.stokes_components = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisDataset.depolarization_value = str(float(depol_value))
                    ThisDataset.atmosphere = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisDataset.particles = particles_string
                    ThisDataset.phase_function_truncation = self.Truncation_method
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisDataset.Fbeam = str(float(self.FBeam[0]))
                        else:
                            ThisDataset.Fbeam = "1.0"
                    if self.Mode == "kdis":
                        ThisDataset.incident_spectrum = self.Incident_spectrum
                        

                    Nw_dim = ThisDataset.createDimension("Nw", Nw)
                    Nsza_dim = ThisDataset.createDimension("Nsza", Nsza)
                    Nvza_dim = ThisDataset.createDimension("Nvza", Nvza)
                    Nvaa_dim = ThisDataset.createDimension("Nvaa", Nvaa)

                    Wav_var = ThisDataset.createVariable("Wv", "f8", ("Nw",))
                    Sza_var = ThisDataset.createVariable("Sza", "f8", ("Nsza",))
                    Vza_var = ThisDataset.createVariable("Vza", "f8", ("Nvza",))
                    Vaa_var = ThisDataset.createVariable("Vaa", "f8", ("Nvaa",))

                    Wav_var.long_name = "Wavelengths"
                    Wav_var.units = "micrometers"
                    Sza_var.long_name = "Solar zenith angles"
                    Sza_var.units = "degrees"
                    Vza_var.long_name = "Viewing zenith angles"
                    Vza_var.units = "degrees"
                    Vaa_var.long_name = "Viewing azimuth angles"
                    Vaa_var.units = "degrees"

                    Wav_var[:] = Wv
                    Sza_var[:] = Sza
                    Vza_var[:] = Vza
                    Vaa_var[:] = Vaa
                    
                    I_var = ThisDataset.createVariable("I", "f8", ("Nw","Nsza","Nvza","Nvaa",))
                    I_var.long_name = "Upward radiance (first Stokes parameter) at TOA"
                    I_var.units = "W/m2/sr"
                    I_var[:] = up_rad[0]

                    if Ns > 1:
                        Q_var = ThisDataset.createVariable("Q", "f8", ("Nw","Nsza","Nvza","Nvaa",))
                        Q_var.long_name = "Upward second Stokes parameter at TOA"
                        Q_var.units = "W/m2/sr"
                        Q_var[:] = up_rad[1]
                        U_var = ThisDataset.createVariable("U", "f8", ("Nw","Nsza","Nvza","Nvaa",))
                        U_var.long_name = "Upward third Stokes parameter at TOA"
                        U_var.units = "W/m2/sr"
                        U_var[:] = up_rad[2]
                        if Ns == 4:
                            V_var = ThisDataset.createVariable("V", "f8", ("Nw","Nsza","Nvza","Nvaa",))
                            U_var.long_name = "Upward fourth Stokes parameter at TOA"
                            V_var.units = "W/m2/sr"
                            V_var[:] = up_rad[3]
                    
                    ThisDataset.close()

                if format == "hdf5" or format == "all":
                    import h5py
                    ThisFile = h5py.File(output_dir + direction + "_radiance.hdf5", "w")
                    ThisFile.attrs['created'] = "Created " + time.ctime(time.time())
                    ThisFile.attrs['Source'] = "ARTDECO V1.1"
                    ThisFile.attrs['URL'] = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisFile.attrs['case_name'] = self.Output_subdirectory
                    ThisFile.attrs['mode'] = self.Mode
                    if self.Mode == "kdis":
                        ThisFile.attrs['k_distribution'] = self.K_distributions_parameterization
                        ThisFile.attrs['wv_min'] = self.Wavelength_start
                        ThisFile.attrs['wv_max'] = self.Wavelength_end
                        filter_list = ""
                        if len(self.Instrumental_filters_list) > 0:
                            filter_list += self.Instrumental_filters_list[0]
                            for iff in range(len(self.Instrumental_filters_list) - 1):
                                filter_list += " " + self.Instrumental_filters_list[iff + 1]
                        else:
                            filter_list = "none"
                        ThisFile.attrs['instrumental_filter'] = filter_list
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisFile.attrs['gas_list'] = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisFile.attrs['gas_continuum_list'] = gas_continuum_list_str
                    ThisFile.attrs['radiative_transfer_model'] = self.RTE_solver
                    ThisFile.attrs['stokes_components'] = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisFile.attrs['depolarization_value'] = str(float(depol_value))
                    ThisFile.attrs['atmosphere'] = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisFile.attrs['phase_function_truncation'] = self.Truncation_method
                    ThisFile.attrs['particles'] = particles_string
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisFile.attrs['Fbeam'] = str(float(self.FBeam[0]))
                        else:
                            ThisFile.attrs['Fbeam'] = "1.0"
                    if self.Mode == "kdis":
                        ThisFile.attrs['incident_spectrum'] = self.Incident_spectrum
                    ThisFile.attrs['Nw'] = Nw
                    ThisFile.attrs['Nsza'] = Nsza
                    ThisFile.attrs['Nvza'] = Nvza
                    ThisFile.attrs['Nvaa'] = Nvaa

                    Wv_dset = ThisFile.create_dataset("Wv", data = Wv, dtype = "f8")
                    Wv_dset.attrs['long_name'] = "Wavelengths"
                    Wv_dset.attrs['units'] = "degrees"
                    Wv_dset.attrs['dimensions'] = "(Nw,)"

                    Sza_dset = ThisFile.create_dataset("Sza", data = Sza, dtype = "f8")
                    Sza_dset.attrs['long_name'] = "Solar zenith angles"
                    Sza_dset.attrs['units'] = "degrees"
                    Sza_dset.attrs['dimensions'] = "(Nsza,)"

                    Vza_dset = ThisFile.create_dataset("Vza", data = Vza, dtype = "f8")
                    Vza_dset.attrs['long_name'] = "Viewing zenith angles"
                    Vza_dset.attrs['units'] = "degrees"
                    Vza_dset.attrs['dimensions'] = "(Nvza,)"

                    Vaa_dset = ThisFile.create_dataset("Vaa", data = Vaa, dtype = "f8")
                    Vaa_dset.attrs['long_name'] = "Viewing azimuthal angles"
                    Vaa_dset.attrs['units'] = "degrees"
                    Vaa_dset.attrs['dimensions'] = "(Nvaa,)"

                    I_np = up_rad[0]
                    I_dset = ThisFile.create_dataset("I", data = I_np, dtype = "f8")
                    del I_np
                    I_dset.attrs['long_name'] = "Upward radiance (first Stokes parameter) at TOA"
                    I_dset.attrs['units'] = "W/m2/sr"
                    I_dset.attrs['dimensions'] = "(Nw, Nsza, Nvza, Nvaa,)"
                    
                    if self.Stokes_components > 1:

                        Q_np = up_rad[0]
                        Q_dset = ThisFile.create_dataset("Q", data = Q_np, dtype = "f8")
                        del Q_np
                        Q_dset.attrs['long_name'] = "Upward second Stokes parameter at TOA"
                        Q_dset.attrs['units'] = "W/m2/sr"
                        Q_dset.attrs['dimensions'] = "(Nw, Nsza, Nvza, Nvaa,)"
                        del Q_dset
                        
                        U_np = up_rad[0]
                        U_dset = ThisFile.create_dataset("U", data = U_np, dtype = "f8")
                        del U_np
                        U_dset.attrs['long_name'] = "Upward third Stokes parameter at TOA"
                        U_dset.attrs['units'] = "W/m2/sr"
                        U_dset.attrs['dimensions'] = "(Nw, Nsza, Nvza, Nvaa,)"
                        del U_dset
                        
                        if self.Stokes_components == 4:
                            V_np = up_rad[0]
                            V_dset = ThisFile.create_dataset("V", data = V_np, dtype = "f8")
                            del V_np
                            V_dset.attrs['long_name'] = "Upward fourth Stokes parameter at TOA"
                            V_dset.attrs['units'] = "W/m2/sr"
                            V_dset.attrs['dimensions'] = "(Nw, Nsza, Nvza, Nvaa,)"
                            del V_dset
                            
                    ThisFile.close()
                    del Wv_dset, Sza_dset, Vza_dset, Vaa_dset, I_dset

            #########################
            ## Reads Betal_files.
            for betal_file in glob.glob(output_dir + "Betal_*.dat"):
                particle_name = (betal_file.split("_")[-1]).split(".da")[0]

                ## Reads ascii file and erase comment lines.
                file_in = open(betal_file, "r")
                raw_lines_in = file_in.readlines()
                file_in.close()
                lines_in = clean_comment_lines(raw_lines_in)
                del raw_lines_in

                
                ## dimensions.
                Nw = int(clean_line(lines_in[1])[0])
                Nbetal = int(clean_line(lines_in[2])[1])

                ## dimension variables.
                Wv = numpy.zeros([Nw], dtype = "f8")
                for iw in range(Nw):
                    Wv[iw] = float(clean_line(lines_in[iw + 2])[0])
                l = numpy.zeros([Nbetal], dtype = "i2")
                for ib in range(Nbetal):
                    l[ib] = ib
                
                ## some variables.
                TruncCoeff = numpy.zeros([Nw], dtype = "f8")
                Cext = numpy.zeros([Nw], dtype = "f8")
                SSA = numpy.zeros([Nw], dtype = "f8")
                g = numpy.zeros([Nw], dtype = "f8")
                alpha1 = numpy.zeros([Nw, Nbetal], dtype = "f8")
                alpha2 = numpy.zeros([Nw, Nbetal], dtype = "f8")
                alpha3 = numpy.zeros([Nw, Nbetal], dtype = "f8")
                alpha4 = numpy.zeros([Nw, Nbetal], dtype = "f8")
                beta1 = numpy.zeros([Nw, Nbetal], dtype = "f8")
                beta2 = numpy.zeros([Nw, Nbetal], dtype = "f8")
                
                ## reads the variables.
                for iw in range(Nw):
                    TruncCoeff[iw] = float(clean_line(lines_in[iw + 2])[2])
                    Cext[iw] = float(clean_line(lines_in[iw + 2])[3])
                    SSA[iw] = float(clean_line(lines_in[iw + 2])[4])
                    g[iw] = float(clean_line(lines_in[iw + 2])[5])
                    
                    for ib in range(Nbetal):
                        alpha1[iw, ib] = float(clean_line(lines_in[iw * Nbetal + ib + Nw + 2])[1])
                        alpha2[iw, ib] = float(clean_line(lines_in[iw * Nbetal + ib + Nw + 2])[2])
                        alpha3[iw, ib] = float(clean_line(lines_in[iw * Nbetal + ib + Nw + 2])[3])
                        alpha4[iw, ib] = float(clean_line(lines_in[iw * Nbetal + ib + Nw + 2])[4])
                        beta1[iw, ib] = float(clean_line(lines_in[iw * Nbetal + ib + Nw + 2])[5])
                        beta2[iw, ib] = float(clean_line(lines_in[iw * Nbetal + ib + Nw + 2])[6])
                           
                if format == "netcdf" or format == "all":
                    import netCDF4, time
                    ThisDataset = netCDF4.Dataset(output_dir + "Betal_" + \
                                                      particle_name + ".nc", "w", 
                                                  format = "NETCDF4")

                    ThisDataset.created = "Created " + time.ctime(time.time())
                    ThisDataset.Source = "ARTDECO V1.1"
                    ThisDataset.URL = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisDataset.case_name = self.Output_subdirectory
                    ThisDataset.mode = self.Mode
                    if self.Mode == "kdis":
                        ThisDataset.k_distribution = self.K_distributions_parameterization
                        filter_list = ""
                        if len(self.Instrumental_filters_list) > 0:
                            filter_list += self.Instrumental_filters_list[0]
                            for iff in range(len(self.Instrumental_filters_list) - 1):
                                filter_list += " " + self.Instrumental_filters_list[iff + 1]
                        else:
                            filter_list = "none"
                        ThisDataset.instrumental_filter = filter_list
                    ThisDataset.particle = particle_name
                    ThisDataset.phase_function_truncation = self.Truncation_method                      

                    Nw_dim = ThisDataset.createDimension("Nw", Nw)
                    Nbetal_dim = ThisDataset.createDimension("Nbetal", Nbetal)

                    Wv_var = ThisDataset.createVariable("Wv", "f8", ("Nw",))
                    Betal_var = ThisDataset.createVariable("l", "f8", ("Nbetal",))

                    Wv_var.long_name = "Wavelengths"
                    Wv_var.units = "micrometers"
                    Betal_var.long_name = "index of Legendre expansion coefficients"

                    Wv_var[:] = Wv
                    Betal_var[:] = l

                    Trunc_var = ThisDataset.createVariable("TruncCoeff", "f8", ("Nw",))
                    Trunc_var.long_name = "Truncation coefficient"
                    Trunc_var.units = ""
                    Trunc_var[:] = TruncCoeff

                    Cext_var = ThisDataset.createVariable("Cext", "f8", ("Nw",))
                    Cext_var.long_name = "Extinction coefficient"
                    Cext_var.units = ""
                    Cext_var[:] = Cext

                    SSA_var = ThisDataset.createVariable("SSA", "f8", ("Nw",))
                    SSA_var.long_name = "Single scattering albedo"
                    SSA_var.units = ""
                    SSA_var[:] = SSA

                    g_var = ThisDataset.createVariable("g", "f8", ("Nw",))
                    g_var.long_name = "Asymmetry coefficient"
                    g_var.unit = ""
                    g_var[:] = g

                    alpha1_var = ThisDataset.createVariable("alpha1", "f8", ("Nw","Nbetal",))
                    alpha1_var.units = ""
                    alpha1_var[:] = alpha1
                    
                    alpha2_var = ThisDataset.createVariable("alpha2", "f8", ("Nw","Nbetal",))
                    alpha2_var.units = ""
                    alpha2_var[:] = alpha2
                
                    alpha3_var = ThisDataset.createVariable("alpha3", "f8", ("Nw","Nbetal",))
                    alpha3_var.units = ""
                    alpha3_var[:] = alpha3
                    
                    alpha4_var = ThisDataset.createVariable("alpha4", "f8", ("Nw","Nbetal",))
                    alpha4_var.units = ""
                    alpha4_var[:] = alpha4
                    
                    beta1_var = ThisDataset.createVariable("beta1", "f8", ("Nw","Nbetal",))
                    beta1_var.units = ""
                    beta1_var[:] = beta1
            
                    beta2_var = ThisDataset.createVariable("beta2", "f8", ("Nw","Nbetal",))
                    beta2_var.units = ""
                    beta2_var[:] = beta2
                    
                    ThisDataset.close()

                if format == "hdf5" or format == "all":
                    import h5py
                    ThisFile = h5py.File(output_dir + "Betal_" + \
                                             particle_name + ".hdf5", "w")
                    ThisFile.attrs['created'] = "Created " + time.ctime(time.time())
                    ThisFile.attrs['Source'] = "ARTDECO V1.1"
                    ThisFile.attrs['URL'] = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisFile.attrs['case_name'] = self.Output_subdirectory
                    ThisFile.attrs['mode'] = self.Mode
                    if self.Mode == "kdis":
                        ThisFile.attrs['k_distribution'] = self.K_distributions_parameterization
                        ThisFile.attrs['wv_min'] = self.Wavelength_start
                        ThisFile.attrs['wv_max'] = self.Wavelength_end
                        filter_list = ""
                        if len(self.Instrumental_filters_list) > 0:
                            filter_list += self.Instrumental_filters_list[0]
                            for iff in range(len(self.Instrumental_filters_list) - 1):
                                filter_list += " " + self.Instrumental_filters_list[iff + 1]
                        else:
                            filter_list = "none"
                        ThisFile.attrs['instrumental_filter'] = filter_list
                    ThisFile.attrs['phase_function_truncation'] = self.Truncation_method
                    ThisFile.attrs['particles'] = particle_name

                    ThisFile.attrs['Nw'] = Nw
                    ThisFile.attrs['Nbetal'] = Nbetal

                    Wv_dset = ThisFile.create_dataset("Wv", data = Wv, dtype = "f8")
                    Wv_dset.attrs['long_name'] = "Wavelengths"
                    Wv_dset.attrs['units'] = "micrometers"
                    Wv_dset.attrs['dimensions'] = "(Nw,)"

                    Trunc_dset = ThisFile.create_dataset("TruncCoeff", data=TruncCoeff, dtype = "f8")
                    Trunc_dset.attrs['dimensions'] = "(Nw,)"
                    Trunc_dset.attrs['long_name'] = "Truncation coefficient"
                    Trunc_dset.attrs['units'] = ""

                    Cext_dset = ThisFile.create_dataset("Cext", data=Cext, dtype = "f8")
                    Cext_dset.attrs['dimensions'] = "(Nw,)"
                    Cext_dset.attrs['long_name'] = "Extinction coefficient"
                    Cext_dset.attrs['units'] = ""

                    SSA_dset = ThisFile.create_dataset("SSA", data=SSA, dtype = "f8")
                    SSA_dset.attrs['dimensions'] = "(Nw,)"
                    SSA_dset.attrs['long_name'] = "Single scattering albedo"
                    SSA_dset.attrs['units'] = ""

                    g_dset = ThisFile.create_dataset("g", data=g, dtype = "f8")
                    g_dset.attrs['dimensions'] = "(Nw,)"
                    g_dset.attrs['long_name'] = "Asymmetry coefficient"
                    g_dset.attrs['units'] = ""

                    alpha1_dset = ThisFile.create_dataset("alpha1", data=alpha1, dtype = "f8")
                    alpha1_dset.attrs['dimensions'] = "(Nw, Nbetal)"
                    alpha1_dset.attrs['long_name'] = ""
                    alpha1_dset.attrs['units'] = ""

                    alpha2_dset = ThisFile.create_dataset("alpha2", data=alpha2, dtype = "f8")
                    alpha2_dset.attrs['dimensions'] = "(Nw, Nbetal)"
                    alpha2_dset.attrs['long_name'] = ""
                    alpha2_dset.attrs['units'] = ""

                    alpha3_dset = ThisFile.create_dataset("alpha3", data=alpha3, dtype = "f8")
                    alpha3_dset.attrs['dimensions'] = "(Nw, Nbetal)"
                    alpha3_dset.attrs['long_name'] = ""
                    alpha3_dset.attrs['units'] = ""

                    alpha4_dset = ThisFile.create_dataset("alpha4", data=alpha4, dtype = "f8")
                    alpha4_dset.attrs['dimensions'] = "(Nw, Nbetal)"
                    alpha4_dset.attrs['long_name'] = ""
                    alpha4_dset.attrs['units'] = ""

                    beta1_dset = ThisFile.create_dataset("beta1", data=beta1, dtype = "f8")
                    beta1_dset.attrs['dimensions'] = "(Nw, Nbetal)"
                    beta1_dset.attrs['long_name'] = ""
                    beta1_dset.attrs['units'] = ""

                    beta2_dset = ThisFile.create_dataset("beta2", data=beta2, dtype = "f8")
                    beta2_dset.attrs['dimensions'] = "(Nw, Nbetal)"
                    beta2_dset.attrs['long_name'] = ""
                    beta2_dset.attrs['units'] = ""

                    ThisFile.close()
                    del Wv_dset, Trunc_dset, Cext_dset, SSA_dset, g_dset, alpha1_dset
                    del alpha2_dset, alpha3_dset, alpha4_dset, beta1_dset, beta2_dset

            #########################
            ## Reads Recomposed phase matrices files.
            for opt_file in glob.glob(output_dir + "Recomposed_*.dat"):
                particle_name = (opt_file.split("_")[-1]).split(".da")[0]

                ## Reads ascii file and erase comment lines.
                file_in = open(opt_file, "r")
                raw_lines_in = file_in.readlines()
                file_in.close()
                lines_in = clean_comment_lines(raw_lines_in)
                del raw_lines_in
                
                ## dimensions.
                Nw = int(clean_line(lines_in[1])[0])
                Nangle = int(clean_line(lines_in[2])[2])

                ## dimension variables.
                Wv = numpy.zeros([Nw], dtype = "f8")
                for iw in range(Nw):
                    Wv[iw] = float(clean_line(lines_in[iw + 2])[0])
                angle = numpy.zeros([Nangle], dtype = "i2")
                for ib in range(Nangle):
                    angle[ib] = ib
                
                ## some variables.
                TruncCoeff = numpy.zeros([Nw], dtype = "f8")
                Cext = numpy.zeros([Nw], dtype = "f8")
                sumF11 = numpy.zeros([Nw], dtype = "f8")
                SSA = numpy.zeros([Nw], dtype = "f8")
                g = numpy.zeros([Nw], dtype = "f8")
                u = numpy.zeros([Nw, Nangle], dtype = "f8")
                F11 = numpy.zeros([Nw, Nangle], dtype = "f8")
                F22 = numpy.zeros([Nw, Nangle], dtype = "f8")
                F33 = numpy.zeros([Nw, Nangle], dtype = "f8")
                F44 = numpy.zeros([Nw, Nangle], dtype = "f8")
                F21 = numpy.zeros([Nw, Nangle], dtype = "f8")
                F34 = numpy.zeros([Nw, Nangle], dtype = "f8")
                
                ## reads the variables.
                for iw in range(Nw):
                    TruncCoeff[iw] = float(clean_line(lines_in[iw + 2])[3])
                    sumF11[iw] = float(clean_line(lines_in[iw + 2])[4])
                    Cext[iw] = float(clean_line(lines_in[iw + 2])[5])
                    SSA[iw] = float(clean_line(lines_in[iw + 2])[6])
                    g[iw] = float(clean_line(lines_in[iw + 2])[7])
                    
                    for ib in range(Nangle):
                        u[iw, ib] = float(clean_line(lines_in[iw * Nangle + ib + Nw + 2])[0])
                        F11[iw, ib] = float(clean_line(lines_in[iw * Nangle + ib + Nw + 2])[1])
                        F22[iw, ib] = float(clean_line(lines_in[iw * Nangle + ib + Nw + 2])[2])
                        F33[iw, ib] = float(clean_line(lines_in[iw * Nangle + ib + Nw + 2])[3])
                        F44[iw, ib] = float(clean_line(lines_in[iw * Nangle + ib + Nw + 2])[4])
                        F21[iw, ib] = float(clean_line(lines_in[iw * Nangle + ib + Nw + 2])[5])
                        F34[iw, ib] = float(clean_line(lines_in[iw * Nangle + ib + Nw + 2])[6])
                           
                if format == "netcdf" or format == "all":
                    import netCDF4, time
                    ThisDataset = netCDF4.Dataset(output_dir + "Recomposed_phase_matrix_" + \
                                                      particle_name + ".nc", "w", 
                                                  format = "NETCDF4")

                    ThisDataset.created = "Created " + time.ctime(time.time())
                    ThisDataset.Source = "ARTDECO V1.1"
                    ThisDataset.URL = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisDataset.case_name = self.Output_subdirectory
                    ThisDataset.mode = self.Mode
                    if self.Mode == "kdis":
                        ThisDataset.k_distribution = self.K_distributions_parameterization
                        filter_list = ""
                        if len(self.Instrumental_filters_list) > 0:
                            filter_list += self.Instrumental_filters_list[0]
                            for iff in range(len(self.Instrumental_filters_list) - 1):
                                filter_list += " " + self.Instrumental_filters_list[iff + 1]
                        else:
                            filter_list = "none"
                        ThisDataset.instrumental_filter = filter_list
                    ThisDataset.particle = particle_name
                    ThisDataset.phase_function_truncation = self.Truncation_method                      

                    Nw_dim = ThisDataset.createDimension("Nw", Nw)
                    Nangle_dim = ThisDataset.createDimension("Nangle", Nangle)

                    Wv_var = ThisDataset.createVariable("Wv", "f8", ("Nw",))
                    Angle_var = ThisDataset.createVariable("angle_index", "f8", ("Nangle",))

                    Wv_var.long_name = "Wavelengths"
                    Wv_var.units = "micrometers"
                    Angle_var.long_name = "index of scattering angle"

                    Wv_var[:] = Wv
                    Angle_var[:] = angle

                    Trunc_var = ThisDataset.createVariable("TruncCoeff", "f8", ("Nw",))
                    Trunc_var.long_name = "Truncation coefficient"
                    Trunc_var.units = ""
                    Trunc_var[:] = TruncCoeff

                    sumF11_var = ThisDataset.createVariable("sumF11", "f8", ("Nw",))
                    sumF11_var.long_name = "Sum of F11 entries of the matrix"
                    sumF11_var.units = ""
                    sumF11_var[:] = sumF11

                    Cext_var = ThisDataset.createVariable("Cext", "f8", ("Nw",))
                    Cext_var.long_name = "Extinction coefficient"
                    Cext_var.units = ""
                    Cext_var[:] = Cext

                    SSA_var = ThisDataset.createVariable("SSA", "f8", ("Nw",))
                    SSA_var.long_name = "Single scattering albedo"
                    SSA_var.units = ""
                    SSA_var[:] = SSA

                    g_var = ThisDataset.createVariable("g", "f8", ("Nw",))
                    g_var.long_name = "Asymmetry coefficient"
                    g_var.unit = ""
                    g_var[:] = g
                    
                    u_var = ThisDataset.createVariable("u", "f8", ("Nw","Nangle",))
                    u_var.units = ""
                    u_var[:] = u
                    
                    F11_var = ThisDataset.createVariable("F11", "f8", ("Nw","Nangle",))
                    F11_var.units = ""
                    F11_var[:] = F11
                    
                    F22_var = ThisDataset.createVariable("F22", "f8", ("Nw","Nangle",))
                    F22_var.units = ""
                    F22_var[:] = F22
                
                    F33_var = ThisDataset.createVariable("F33", "f8", ("Nw","Nangle",))
                    F33_var.units = ""
                    F33_var[:] = F33
                    
                    F44_var = ThisDataset.createVariable("F44", "f8", ("Nw","Nangle",))
                    F44_var.units = ""
                    F44_var[:] = F44
                    
                    F21_var = ThisDataset.createVariable("F21", "f8", ("Nw","Nangle",))
                    F21_var.units = ""
                    F21_var[:] = F21
            
                    F34_var = ThisDataset.createVariable("F34", "f8", ("Nw","Nangle",))
                    F34_var.units = ""
                    F34_var[:] = F34
                    
                    ThisDataset.close()

                if format == "hdf5" or format == "all":
                    import h5py
                    ThisFile = h5py.File(output_dir + "Recomposed_phase_matrix_" + \
                                                      particle_name + ".hdf5", "w")
                    ThisFile.attrs['created'] = "Created " + time.ctime(time.time())
                    ThisFile.attrs['Source'] = "ARTDECO V1.1"
                    ThisFile.attrs['URL'] = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisFile.attrs['case_name'] = self.Output_subdirectory
                    ThisFile.attrs['mode'] = self.Mode
                    if self.Mode == "kdis":
                        ThisFile.attrs['k_distribution'] = self.K_distributions_parameterization
                        ThisFile.attrs['wv_min'] = self.Wavelength_start
                        ThisFile.attrs['wv_max'] = self.Wavelength_end
                        filter_list = ""
                        if len(self.Instrumental_filters_list) > 0:
                            filter_list += self.Instrumental_filters_list[0]
                            for iff in range(len(self.Instrumental_filters_list) - 1):
                                filter_list += " " + self.Instrumental_filters_list[iff + 1]
                        else:
                            filter_list = "none"
                        ThisFile.attrs['instrumental_filter'] = filter_list
                    ThisFile.attrs['phase_function_truncation'] = self.Truncation_method
                    ThisFile.attrs['particles'] = particle_name

                    ThisFile.attrs['Nw'] = Nw
                    ThisFile.attrs['Nangle'] = Nangle

                    Wv_dset = ThisFile.create_dataset("Wv", data = Wv, dtype = "f8")
                    Wv_dset.attrs['long_name'] = "Wavelengths"
                    Wv_dset.attrs['units'] = "micrometers"
                    Wv_dset.attrs['dimensions'] = "(Nw,)"

                    Trunc_dset = ThisFile.create_dataset("TruncCoeff", data=TruncCoeff, dtype = "f8")
                    Trunc_dset.attrs['dimensions'] = "(Nw,)"
                    Trunc_dset.attrs['long_name'] = "Truncation coefficient"
                    Trunc_dset.attrs['units'] = ""

                    sumF11_dset = ThisFile.create_dataset("sumF11", data=sumF11, dtype = "f8")
                    sumF11_dset.attrs['dimensions'] = "(Nw,)"
                    sumF11_dset.attrs['long_name'] = "Sum of F11 entries of the matrix"
                    sumF11_dset.attrs['units'] = ""

                    Cext_dset = ThisFile.create_dataset("Cext", data=Cext, dtype = "f8")
                    Cext_dset.attrs['dimensions'] = "(Nw,)"
                    Cext_dset.attrs['long_name'] = "Extinction coefficient"
                    Cext_dset.attrs['units'] = ""

                    SSA_dset = ThisFile.create_dataset("SSA", data=SSA, dtype = "f8")
                    SSA_dset.attrs['dimensions'] = "(Nw,)"
                    SSA_dset.attrs['long_name'] = "Single scattering albedo"
                    SSA_dset.attrs['units'] = ""

                    g_dset = ThisFile.create_dataset("g", data=g, dtype = "f8")
                    g_dset.attrs['dimensions'] = "(Nw,)"
                    g_dset.attrs['long_name'] = "Asymmetry coefficient"
                    g_dset.attrs['units'] = ""

                    u_dset = ThisFile.create_dataset("u", data=u, dtype = "f8")
                    u_dset.attrs['dimensions'] = "(Nw, Nangle)"
                    u_dset.attrs['long_name'] = ""
                    u_dset.attrs['units'] = ""

                    F11_dset = ThisFile.create_dataset("F11", data=F11, dtype = "f8")
                    F11_dset.attrs['dimensions'] = "(Nw, Nangle)"
                    F11_dset.attrs['long_name'] = ""
                    F11_dset.attrs['units'] = ""

                    F22_dset = ThisFile.create_dataset("F22", data=F22, dtype = "f8")
                    F22_dset.attrs['dimensions'] = "(Nw, Nangle)"
                    F22_dset.attrs['long_name'] = ""
                    F22_dset.attrs['units'] = ""

                    F33_dset = ThisFile.create_dataset("F33", data=F33, dtype = "f8")
                    F33_dset.attrs['dimensions'] = "(Nw,  Nangle)"
                    F33_dset.attrs['long_name'] = ""
                    F33_dset.attrs['units'] = ""

                    F44_dset = ThisFile.create_dataset("F44", data=F44, dtype = "f8")
                    F44_dset.attrs['dimensions'] = "(Nw,  Nangle)"
                    F44_dset.attrs['long_name'] = ""
                    F44_dset.attrs['units'] = ""

                    F21_dset = ThisFile.create_dataset("F21", data=F21, dtype = "f8")
                    F21_dset.attrs['dimensions'] = "(Nw,  Nangle)"
                    F21_dset.attrs['long_name'] = ""
                    F21_dset.attrs['units'] = ""

                    F34_dset = ThisFile.create_dataset("F34", data=F34, dtype = "f8")
                    F34_dset.attrs['dimensions'] = "(Nw,  Nangle"
                    F34_dset.attrs['long_name'] = ""
                    F34_dset.attrs['units'] = ""

                    ThisFile.close()
                    del Wv_dset, Trunc_dset, Cext_dset, SSA_dset, g_dset, F11_dset
                    del F22_dset, F33_dset, F44_dset, F21_dset, F34_dset
                    del u_dset, sumF11
            
            ## Read upward radiances for instrument filters.
            for filter_radiance_file in glob.glob(output_dir + "Upward_radiance_*.dat"):
                filter_name = (filter_radiance_file.split("_")[-1]).split(".da")[0]

                ## Reads ascii file and erase comment lines.
                file_in = open(filter_radiance_file, "r")
                raw_lines_in = file_in.readlines()
                file_in.close()
                lines_in = clean_comment_lines(raw_lines_in)
                del raw_lines_in

                ## Number of Stokes parameters.
                Ns = int(clean_line(lines_in[2])[0])
                
                ## dimensions.
                Nsza = int(clean_line(lines_in[1])[0])
                Nvza = int(clean_line(lines_in[1])[1])
                Nvaa = int(clean_line(lines_in[1])[2])

                ## dimension variables.
                Sza = numpy.zeros([Nsza], dtype = "f8")
                Vza = numpy.zeros([Nvza], dtype = "f8")
                Vaa = numpy.zeros([Nvaa], dtype = "f8")

                ## Upward radiance array.
                up_rad = numpy.zeros([Ns, Nsza, Nvza, Nvaa], dtype = "f8")

                for isza in range(Nsza):
                    Sza[isza] = float(clean_line(lines_in[isza * Nvza * Nvaa + 4])[0])
                for ivza in range(Nvza):
                    Vza[ivza] = float(clean_line(lines_in[ivza * Nvaa + 4])[1])
                for ivaa in range(Nvaa):
                    Vaa[ivaa] = float(clean_line(lines_in[ivaa + 4])[2])         
                                
                for isza in range(Nsza):
                    for ivza in range(Nvza):
                        for ivaa in range(Nvaa):
                            for iparam in range(Ns):
                                up_rad[iparam, isza, ivza, ivaa] = \
                                    float(clean_line(lines_in[isza * Nvza * Nvaa + \
                                                                  ivza * Nvaa + \
                                                                  ivaa + 4])[4 + iparam])
                           
                if format == "netcdf" or format == "all":
                    import netCDF4, time
                    ThisDataset = netCDF4.Dataset(output_dir + "Upward_radiance_" + \
                                                      filter_name + ".nc", "w", 
                                                  format = "NETCDF4")

                    ThisDataset.created = "Created " + time.ctime(time.time())
                    ThisDataset.Source = "ARTDECO V1.1"
                    ThisDataset.URL = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisDataset.case_name = self.Output_subdirectory
                    ThisDataset.mode = self.Mode
                    if self.Mode == "kdis":
                        ThisDataset.k_distribution = self.K_distributions_parameterization
                        ThisDataset.instrumental_filter = filter_name
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisDataset.gas_list = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisDataset.gas_continuum_list = gas_continuum_list_str
                    ThisDataset.radiative_transfer_model = self.RTE_solver
                    ThisDataset.stokes_components = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisDataset.depolarization_value = str(float(depol_value))
                    ThisDataset.atmosphere = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisDataset.phase_function_truncation = self.Truncation_method
                    ThisDataset.particles = particles_string
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisDataset.Fbeam = str(float(self.FBeam[0]))
                        else:
                            ThisDataset.Fbeam = "1.0"
                    if self.Mode == "kdis":
                        ThisDataset.incident_spectrum = self.Incident_spectrum                        

                    Nsza_dim = ThisDataset.createDimension("Nsza", Nsza)
                    Nvza_dim = ThisDataset.createDimension("Nvza", Nvza)
                    Nvaa_dim = ThisDataset.createDimension("Nvaa", Nvaa)

                    Sza_var = ThisDataset.createVariable("Sza", "f8", ("Nsza",))
                    Vza_var = ThisDataset.createVariable("Vza", "f8", ("Nvza",))
                    Vaa_var = ThisDataset.createVariable("Vaa", "f8", ("Nvaa",))

                    Sza_var.long_name = "Solar zenith angles"
                    Sza_var.units = "degrees"
                    Vza_var.long_name = "Viewing zenith angles"
                    Vza_var.units = "degrees"
                    Vaa_var.long_name = "Viewing azimuth angles"
                    Vaa_var.units = "degrees"

                    Sza_var[:] = Sza
                    Vza_var[:] = Vza
                    Vaa_var[:] = Vaa
                    
                    I_var = ThisDataset.createVariable("I", "f8", ("Nsza","Nvza","Nvaa",))
                    I_var.units = "W/m2/sr"
                    I_var[:] = up_rad[0]

                    if Ns > 1:
                        Q_var = ThisDataset.createVariable("Q", "f8", ("Nsza","Nvza","Nvaa",))
                        Q_var.units = "W/m2/sr"
                        Q_var[:] = up_rad[1]
                        U_var = ThisDataset.createVariable("U", "f8", ("Nsza","Nvza","Nvaa",))
                        U_var.units = "W/m2/sr"
                        U_var[:] = up_rad[2]
                        if Ns == 4:
                            V_var = ThisDataset.createVariable("V", "f8", ("Nsza","Nvza","Nvaa",))
                            V_var.units = "W/m2/sr"
                            V_var[:] = up_rad[3]
                    
                    ThisDataset.close()

                if format == "hdf5" or format == "all":
                    import h5py
                    ThisFile = h5py.File(output_dir + "Upward_radiance_" + \
                                             filter_name + ".hdf5", "w")
                    ThisFile.attrs['created'] = "Created " + time.ctime(time.time())
                    ThisFile.attrs['Source'] = "ARTDECO V1.1"
                    ThisFile.attrs['URL'] = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisFile.attrs['case_name'] = self.Output_subdirectory
                    ThisFile.attrs['mode'] = self.Mode
                    if self.Mode == "kdis":
                        ThisFile.attrs['k_distribution'] = self.K_distributions_parameterization
                        ThisFile.attrs['wv_min'] = self.Wavelength_start
                        ThisFile.attrs['wv_max'] = self.Wavelength_end
                        filter_list = ""
                        if len(self.Instrumental_filters_list) > 0:
                            filter_list += self.Instrumental_filters_list[0]
                            for iff in range(len(self.Instrumental_filters_list) - 1):
                                filter_list += " " + self.Instrumental_filters_list[iff + 1]
                        else:
                            filter_list = "none"
                        ThisFile.attrs['instrumental_filter'] = filter_list
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisFile.attrs['gas_list'] = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisFile.attrs['gas_continuum_list'] = gas_continuum_list_str
                    ThisFile.attrs['radiative_transfer_model'] = self.RTE_solver
                    ThisFile.attrs['stokes_components'] = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisFile.attrs['depolarization_value'] = str(float(depol_value))
                    ThisFile.attrs['atmosphere'] = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisFile.attrs['phase_function_truncation'] = self.Truncation_method
                    ThisFile.attrs['particles'] = particles_string
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisFile.attrs['Fbeam'] = str(float(self.FBeam[0]))
                        else:
                            ThisFile.attrs['Fbeam'] = "1.0"
                    if self.Mode == "kdis":
                        ThisFile.attrs['incident_spectrum'] = self.Incident_spectrum
                    ThisFile.attrs['Nsza'] = Nsza
                    ThisFile.attrs['Nvza'] = Nvza
                    ThisFile.attrs['Nvaa'] = Nvaa

                    Sza_dset = ThisFile.create_dataset("Sza", data = Sza, dtype = "f8")
                    Sza_dset.attrs['long_name'] = "Solar zenith angles"
                    Sza_dset.attrs['units'] = "degrees"
                    Sza_dset.attrs['dimensions'] = "(Nsza,)"

                    Vza_dset = ThisFile.create_dataset("Vza", data = Vza, dtype = "f8")
                    Vza_dset.attrs['long_name'] = "Viewing zenith angles"
                    Vza_dset.attrs['units'] = "degrees"
                    Vza_dset.attrs['dimensions'] = "(Nvza,)"

                    Vaa_dset = ThisFile.create_dataset("Vaa", data = Vaa, dtype = "f8")
                    Vaa_dset.attrs['long_name'] = "Viewing azimuthal angles"
                    Vaa_dset.attrs['units'] = "degrees"
                    Vaa_dset.attrs['dimensions'] = "(Nvaa,)"

                    I_np = up_rad[0]
                    I_dset = ThisFile.create_dataset("I", data = I_np, dtype = "f8")
                    del I_np
                    I_dset.attrs['long_name'] = "Upward radiance (first Stokes parameter) at TOA"
                    I_dset.attrs['units'] = "W/m2/sr"
                    I_dset.attrs['dimensions'] = "(Nsza, Nvza, Nvaa,)"
                    
                    if self.Stokes_components > 1:

                        Q_np = up_rad[0]
                        Q_dset = ThisFile.create_dataset("Q", data = Q_np, dtype = "f8")
                        del Q_np
                        Q_dset.attrs['long_name'] = "Upward second Stokes parameter at TOA"
                        Q_dset.attrs['units'] = "W/m2/sr"
                        Q_dset.attrs['dimensions'] = "(Nsza, Nvza, Nvaa,)"
                        del Q_dset
                        
                        U_np = up_rad[0]
                        U_dset = ThisFile.create_dataset("U", data = U_np, dtype = "f8")
                        del U_np
                        U_dset.attrs['long_name'] = "Upward third Stokes parameter at TOA"
                        U_dset.attrs['units'] = "W/m2/sr"
                        U_dset.attrs['dimensions'] = "(Nsza, Nvza, Nvaa,)"
                        del U_dset
                        
                        if self.Stokes_components == 4:
                            V_np = up_rad[0]
                            V_dset = ThisFile.create_dataset("V", data = V_np, dtype = "f8")
                            del V_np
                            V_dset.attrs['long_name'] = "Upward fourth Stokes parameter at TOA"
                            V_dset.attrs['units'] = "W/m2/sr"
                            V_dset.attrs['dimensions'] = "(Nsza, Nvza, Nvaa,)"
                            del V_dset
                            
                    ThisFile.close()
                    del Sza_dset, Vza_dset, Vaa_dset, I_dset
            
            ## Read incoming radiation flux for kdis.
            for solrad_file in glob.glob(output_dir + "Incoming_radiation_flux.dat"):

                ## Reads ascii file and erase comment lines.
                file_in = open(solrad_file, "r")
                raw_lines_in = file_in.readlines()
                file_in.close()
                lines_in = clean_comment_lines(raw_lines_in)
                del raw_lines_in
                
                ## dimensions.
                Nw = int(clean_line(lines_in[0])[0])
                Nsza = int(clean_line(lines_in[1])[0])

                ## dimension variables.
                Wv = numpy.zeros([Nw], dtype = "f8")
                Sza = numpy.zeros([Nsza], dtype = "f8")

                ## Incoming flux.
                F0 = numpy.zeros([Nw], dtype = "f8")
                ## Varsol, lon, lat, day of year and hour (UT).
                varsol = numpy.zeros([Nsza], dtype = "f8")
                lon = numpy.zeros([Nsza], dtype = "f8")
                lat = numpy.zeros([Nsza], dtype = "f8")
                day = numpy.zeros([Nsza], dtype = "f8")
                hour = numpy.zeros([Nsza], dtype = "f8")

                for iw in range(Nw):
                    Wv[iw] = float(clean_line(lines_in[iw + 2])[0])
                    F0[iw] = float(clean_line(lines_in[iw + 2])[1])
                for isza in range(Nsza):
                    Sza[isza] = float(clean_line(lines_in[isza +Nw + 2])[0]) 
                    varsol[isza] = float(clean_line(lines_in[isza +Nw + 2])[1]) 
                    lon[isza] = float(clean_line(lines_in[isza +Nw + 2])[2]) 
                    lat[isza] = float(clean_line(lines_in[isza +Nw + 2])[3]) 
                    day[isza] = float(clean_line(lines_in[isza +Nw + 2])[4]) 
                    hour[isza] = float(clean_line(lines_in[isza +Nw + 2])[5])      
                           
                if format == "netcdf" or format == "all":
                    import netCDF4, time
                    ThisDataset = netCDF4.Dataset(output_dir + "Incoming_radiation_flux.nc", "w", 
                                                  format = "NETCDF4")

                    ThisDataset.created = "Created " + time.ctime(time.time())
                    ThisDataset.Source = "ARTDECO V1.1"
                    ThisDataset.URL = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisDataset.case_name = self.Output_subdirectory
                    ThisDataset.mode = self.Mode
                    if self.Mode == "kdis":
                        ThisDataset.k_distribution = self.K_distributions_parameterization
                        filter_list = ""
                        if len(self.Instrumental_filters_list) > 0:
                            filter_list += self.Instrumental_filters_list[0]
                            for iff in range(len(self.Instrumental_filters_list) - 1):
                                filter_list += " " + self.Instrumental_filters_list[iff + 1]
                        else:
                            filter_list = "none"
                        ThisDataset.instrumental_filter = filter_list
                        ThisDataset.incident_spectrum = self.Incident_spectrum                        

                    Nw_dim = ThisDataset.createDimension("Nw", Nw)
                    Nsza_dim = ThisDataset.createDimension("Nsza", Nsza)

                    Wv_var = ThisDataset.createVariable("Wv", "f8", ("Nw",)) 
                    Sza_var = ThisDataset.createVariable("Sza", "f8", ("Nsza",))

                    Wv_var.long_name = "Wavelengths"
                    Wv_var.units = "micrometers"
                    Sza_var.long_name = "Solar zenith angles"
                    Sza_var.units = "degrees"

                    Wv_var[:] = Wv
                    Sza_var[:] = Sza
                    
                    F0_var = ThisDataset.createVariable("F0", "f8", ("Nw",))
                    F0_var.long_name = "Integrated incoming flux over the band"
                    F0_var.units = "W/m2"
                    F0_var[:] = F0

                    varsol_var = ThisDataset.createVariable("varsol", "f8", ("Nsza",))
                    varsol_var.long_name = "Corrective coefficient for the incoming radiation flux, taking into account the position and dates."
                    varsol_var.units = ""
                    varsol_var[:] = varsol

                    lon_var = ThisDataset.createVariable("lon", "f8", ("Nsza",))
                    lon_var.long_name = "Longitude"
                    lon_var.units = "degrees_west"
                    lon_var[:] = lon

                    lat_var = ThisDataset.createVariable("lat", "f8", ("Nsza",))
                    lat_var.long_name = "Latitude"
                    lat_var.units = "degrees_north"
                    lat_var[:] = lat

                    day_var = ThisDataset.createVariable("day", "f8", ("Nsza",))
                    day_var.long_name = "Day of the year"
                    day_var.units = ""
                    day_var[:] = day

                    hour_var = ThisDataset.createVariable("hour", "f8", ("Nsza",))
                    hour_var.long_name = "Hour of the day (UT)"
                    hour_var.units = ""
                    hour_var[:] = hour
                    
                    ThisDataset.close()

                if format == "hdf5" or format == "all":
                    import h5py
                    ThisFile = h5py.File(output_dir + "Incoming_radiation_flux.hdf5", "w")
                    ThisFile.attrs['created'] = "Created " + time.ctime(time.time())
                    ThisFile.attrs['Source'] = "ARTDECO V1.1"
                    ThisFile.attrs['URL'] = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisFile.attrs['case_name'] = self.Output_subdirectory
                    ThisFile.attrs['mode'] = self.Mode
                    if self.Mode == "kdis":
                        ThisFile.attrs['k_distribution'] = self.K_distributions_parameterization
                        ThisFile.attrs['wv_min'] = self.Wavelength_start
                        ThisFile.attrs['wv_max'] = self.Wavelength_end
                        filter_list = ""
                        if len(self.Instrumental_filters_list) > 0:
                            filter_list += self.Instrumental_filters_list[0]
                            for iff in range(len(self.Instrumental_filters_list) - 1):
                                filter_list += " " + self.Instrumental_filters_list[iff + 1]
                        else:
                            filter_list = "none"
                        ThisFile.attrs['instrumental_filter'] = filter_list
                    if self.Mode == "kdis":
                        ThisFile.attrs['incident_spectrum'] = self.Incident_spectrum

                    ThisFile.attrs['Nw'] = Nw
                    ThisFile.attrs['Nsza'] = Nsza

                    Wv_dset = ThisFile.create_dataset("Wv", data = Wv, dtype = "f8")
                    Wv_dset.attrs['long_name'] = "Wavelengths"
                    Wv_dset.attrs['units'] = "micrometers"
                    Wv_dset.attrs['dimensions'] = "(Nw,)"

                    Sza_dset = ThisFile.create_dataset("Sza", data = Sza, dtype = "f8")
                    Sza_dset.attrs['long_name'] = "Solar zenith angles"
                    Sza_dset.attrs['units'] = "degrees"
                    Sza_dset.attrs['dimensions'] = "(Nsza,)"

                    F0_dset = ThisFile.create_dataset("F0", data = F0, dtype = "f8")
                    F0_dset.attrs['long_name'] = "Integrated incoming flux over the band"
                    F0_dset.attrs['units'] = "W/m2"
                    F0_dset.attrs['dimensions'] = "(Nw,)"

                    varsol_dset = ThisFile.create_dataset("varsol", data = varsol, dtype = "f8")
                    varsol_dset.attrs['long_name'] = "Corrective coefficient for the incoming radiation flux, taking into account the position and dates"
                    varsol_dset.attrs['units'] = ""
                    varsol_dset.attrs['dimensions'] = "(Nsza,)"

                    lon_dset = ThisFile.create_dataset("lon", data = lon, dtype = "f8")
                    lon_dset.attrs['long_name'] = "Longitude"
                    lon_dset.attrs['units'] = "degrees_west"
                    lon_dset.attrs['dimensions'] = "(Nsza,)"

                    lat_dset = ThisFile.create_dataset("lat", data = lat, dtype = "f8")
                    lat_dset.attrs['long_name'] = "Latitude"
                    lat_dset.attrs['units'] = "degrees_north"
                    lat_dset.attrs['dimensions'] = "(Nsza,)"

                    day_dset = ThisFile.create_dataset("day", data = day, dtype = "f8")
                    day_dset.attrs['long_name'] = "day of the year"
                    day_dset.attrs['units'] = ""
                    day_dset.attrs['dimensions'] = "(Nsza,)"

                    hour_dset = ThisFile.create_dataset("hour", data = hour, dtype = "f8")
                    hour_dset.attrs['long_name'] = "hour of the day (UT)"
                    hour_dset.attrs['units'] = ""
                    hour_dset.attrs['dimensions'] = "(Nsza,)"

                    ThisFile.close()
                    del Sza_dset, Wv_dset, F0_dset, lon_dset, lat_dset, day_dset, hour_dset

            ## Read integrated flux.
            for flux_file in glob.glob(output_dir + "Integrated_flux*.dat"):
                filter_name = "none"
                if flux_file.split("_")[-1] != "flux.dat":
                    filter_name = (flux_file.split("_")[-1]).split(".da")[0]
                filename = flux_file.split(".da")[0]

                ## Reads ascii file and erase comment lines.
                file_in = open(flux_file, "r")
                raw_lines_in = file_in.readlines()
                file_in.close()
                lines_in = clean_comment_lines(raw_lines_in)
                del raw_lines_in
                
                ## dimensions.
                Nsza = int(clean_line(lines_in[0])[0])
                Nalt = int(clean_line(lines_in[1])[0])

                ## dimension variables.
                Sza = numpy.zeros([Nsza], dtype = "f8")
                Alt = numpy.zeros([Nalt], dtype = "f8")

                ## Pressure field
                Pres = numpy.zeros([Nalt], dtype = "f8")
                flux_dw = numpy.zeros([Nsza, Nalt], dtype = "f8")
                flux_up = numpy.zeros([Nsza, Nalt], dtype = "f8")
                flux_dir = numpy.zeros([Nsza, Nalt], dtype = "f8")
                
                for isza in range(Nsza):
                    Sza[isza] = float(clean_line(lines_in[isza * (Nalt + 1) + 2])[0])
                for ia in range(Nalt):
                    Alt[ia] = float(clean_line(lines_in[ia + 3])[0])
                    Pres[ia] = float(clean_line(lines_in[ia + 3])[1])
                    for isza in range(Nsza):
                        flux_dw[isza, ia] = float(clean_line(lines_in[isza * (Nalt + 1) + ia + 3])[2])
                        flux_up[isza, ia] = float(clean_line(lines_in[isza * (Nalt + 1) + ia + 3])[3])
                        flux_dir[isza, ia] = float(clean_line(lines_in[isza * (Nalt + 1) + ia + 3])[4])
                           
                if format == "netcdf" or format == "all":
                    import netCDF4, time

                    ThisDataset = netCDF4.Dataset(filename + ".nc", "w", 
                                                  format = "NETCDF4")
                    
                    ThisDataset.created = "Created " + time.ctime(time.time())
                    ThisDataset.Source = "ARTDECO V1.1"
                    ThisDataset.URL = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisDataset.case_name = self.Output_subdirectory
                    ThisDataset.mode = self.Mode
                    if self.Mode == "kdis":
                        ThisDataset.k_distribution = self.K_distributions_parameterization
                        ThisDataset.instrumental_filter = filter_name
                        if filter_name == "none":
                            ThisDataset.wv_min = self.Wavelength_start
                            ThisDataset.wv_max = self.Wavelength_end
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisDataset.gas_list = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisDataset.gas_continuum_list = gas_continuum_list_str
                    ThisDataset.radiative_transfer_model = self.RTE_solver
                    ThisDataset.stokes_components = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisDataset.depolarization_value = str(float(depol_value))
                    ThisDataset.atmosphere = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisDataset.phase_function_truncation = self.Truncation_method
                    ThisDataset.particles = particles_string
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisDataset.Fbeam = str(float(self.FBeam[0]))
                        else:
                            ThisDataset.Fbeam = "1.0"
                    if self.Mode == "kdis":
                        ThisDataset.incident_spectrum = self.Incident_spectrum                      

                    Nsza_dim = ThisDataset.createDimension("Nsza", Nsza)
                    Nalt_dim = ThisDataset.createDimension("Nalt", Nalt)

                    Alt_var = ThisDataset.createVariable("Altitude", "f8", ("Nalt",)) 
                    Sza_var = ThisDataset.createVariable("Sza", "f8", ("Nsza",))

                    Alt_var.long_name = "Altitude of the atmospheric layer interface"
                    Alt_var.units = "kilometers"
                    Sza_var.long_name = "Solar zenith angles"
                    Sza_var.units = "degrees"

                    Alt_var[:] = Alt
                    Sza_var[:] = Sza
                    
                    Pres_var = ThisDataset.createVariable("Pres", "f8", ("Nalt",))
                    Pres_var.long_name = "Air pressure at the given altitude"
                    Pres_var.units = "hPa"
                    Pres_var[:] = Pres

                    Flux_dw_var = ThisDataset.createVariable("flux_dw", "f8", ("Nsza", "Nalt",))
                    Flux_dw_var.long_name = "Total downward integrated flux"
                    Flux_dw_var.units = "W/m2"
                    Flux_dw_var[:] = flux_dw

                    Flux_dir_var = ThisDataset.createVariable("flux_dir", "f8", ("Nsza", "Nalt",))
                    Flux_dir_var.long_name = "Directional downward integrated flux"
                    Flux_dir_var.units = "W/m2"
                    Flux_dir_var[:] = flux_dir

                    Flux_up_var = ThisDataset.createVariable("flux_up", "f8", ("Nsza", "Nalt",))
                    Flux_up_var.long_name = "Total upward integrated flux"
                    Flux_up_var.units = "W/m2"
                    Flux_up_var[:] = flux_up
                    
                    ThisDataset.close()

                if format == "hdf5" or format == "all":
                    import h5py
                    ThisFile = h5py.File(filename + ".hdf5", "w")
                    ThisFile.attrs['created'] = "Created " + time.ctime(time.time())
                    ThisFile.attrs['Source'] = "ARTDECO V1.1"
                    ThisFile.attrs['URL'] = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisFile.attrs['case_name'] = self.Output_subdirectory
                    ThisFile.attrs['mode'] = self.Mode
                    if self.Mode == "kdis":
                        ThisFile.attrs['k_distribution'] = self.K_distributions_parameterization
                        ThisFile.attrs['wv_min'] = self.Wavelength_start
                        ThisFile.attrs['wv_max'] = self.Wavelength_end
                        filter_list = ""
                        if len(self.Instrumental_filters_list) > 0:
                            filter_list += self.Instrumental_filters_list[0]
                            for iff in range(len(self.Instrumental_filters_list) - 1):
                                filter_list += " " + self.Instrumental_filters_list[iff + 1]
                        else:
                            filter_list = "none"
                        ThisFile.attrs['instrumental_filter'] = filter_list
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisFile.attrs['gas_list'] = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisFile.attrs['gas_continuum_list'] = gas_continuum_list_str
                    ThisFile.attrs['radiative_transfer_model'] = self.RTE_solver
                    ThisFile.attrs['stokes_components'] = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisFile.attrs['depolarization_value'] = str(float(depol_value))
                    ThisFile.attrs['atmosphere'] = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisFile.attrs['phase_function_truncation'] = self.Truncation_method
                    ThisFile.attrs['particles'] = particles_string
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisFile.attrs['Fbeam'] = str(float(self.FBeam[0]))
                        else:
                            ThisFile.attrs['Fbeam'] = "1.0"
                    if self.Mode == "kdis":
                        ThisFile.attrs['incident_spectrum'] = self.Incident_spectrum  
                    ThisFile.attrs['Nalt'] = Nalt
                    ThisFile.attrs['Nsza'] = Nsza
                    

                    Alt_dset = ThisFile.create_dataset("Altitude", data = Alt, dtype = "f8")
                    Sza_dset = ThisFile.create_dataset("Sza", data = Sza, dtype = "f8")

                    Alt_dset.attrs['long_name'] = "Altitude of the atmospheric layer interface"
                    Alt_dset.attrs['units'] = "kilometers"
                    Alt_dset.attrs['dimensions'] = "(Nalt,)"
                    Sza_dset.attrs['long_name'] = "Solar zenith angles"
                    Sza_dset.attrs['units'] = "degrees"
                    Sza_dset.attrs['dimensions'] = "(Nsza,)"
                    
                    Pres_dset = ThisFile.create_dataset("Pres", data = Pres, dtype = "f8")
                    Pres_dset.attrs['long_name'] = "Air pressure at the given altitude"
                    Pres_dset.attrs['units'] = "hPa"
                    Pres_dset.attrs['dimensions'] = "(Nalt,)"

                    Flux_dw_dset = ThisFile.create_dataset("flux_dw", data=flux_dw, dtype = "f8") 
                    Flux_dw_dset.attrs['dimensions'] = "(Nsza, Nalt,)"
                    Flux_dw_dset.attrs['long_name'] = "Total downward integrated flux"
                    Flux_dw_dset.attrs['units'] = "W/m2"

                    Flux_dir_dset = ThisFile.create_dataset("flux_dir", data=flux_dir, dtype = "f8")
                    Flux_dir_dset.attrs['dimensions'] = "(Nsza, Nalt,)"
                    Flux_dir_dset.attrs['long_name'] = "Directional downward integrated flux"
                    Flux_dir_dset.attrs['units'] = "W/m2"

                    Flux_up_dset = ThisFile.create_dataset("flux_up", data=flux_up, dtype = "f8")
                    Flux_up_dset.attrs['dimensions'] = "(Nsza, Nalt,)"
                    Flux_up_dset.attrs['long_name'] = "Total upward integrated flux"
                    Flux_up_dset.attrs['units'] = "W/m2"       

                    ThisFile.close()
                    del Sza_dset, Alt_dset, Pres_dset, Flux_dw_dset, Flux_up_dset, Flux_dir_dset

            ## Read kdis bands flux.
            for flux_file in glob.glob(output_dir + "Kdis_bands_flux.dat"):
                filename = flux_file.split(".da")[0]

                ## Reads ascii file and erase comment lines.
                file_in = open(flux_file, "r")
                raw_lines_in = file_in.readlines()
                file_in.close()
                lines_in = clean_comment_lines(raw_lines_in)
                del raw_lines_in
                
                ## dimensions.
                Nw = int(clean_line(lines_in[0])[0])
                Nsza = int(clean_line(lines_in[1])[0])
                Nalt = int(clean_line(lines_in[2])[0])

                ## dimension variables.
                WV = numpy.zeros([Nw], dtype = "f8")
                Sza = numpy.zeros([Nsza], dtype = "f8")
                Alt = numpy.zeros([Nalt], dtype = "f8")

                ## Pressure field
                Pres = numpy.zeros([Nalt], dtype = "f8")
                flux_dw = numpy.zeros([Nw, Nsza, Nalt], dtype = "f8")
                flux_up = numpy.zeros([Nw, Nsza, Nalt], dtype = "f8")
                flux_dir = numpy.zeros([Nw, Nsza, Nalt], dtype = "f8")
                
                for iw in range(Nw):
                    Wv[iw] = float(clean_line(lines_in[iw * (Nsza * (Nalt + 1) + 1) + 3])[0])
                for isza in range(Nsza):
                    Sza[isza] = float(clean_line(lines_in[isza * (Nalt + 1) + 4])[0])
                for ia in range(Nalt):
                    Alt[ia] = float(clean_line(lines_in[ia + 5])[0])
                    Pres[ia] = float(clean_line(lines_in[ia + 5])[1])
                    for isza in range(Nsza):
                        for iw in range(Nw):
                            index = iw * (Nsza * (Nalt + 1) + 1) + isza * (Nalt + 1) + ia + 5
                            flux_dw[iw, isza, ia] = float(clean_line(lines_in[index])[2])
                            flux_up[iw, isza, ia] = float(clean_line(lines_in[index])[3])
                            flux_dir[iw, isza, ia] = float(clean_line(lines_in[index])[4])
                           
                if format == "netcdf" or format == "all":
                    import netCDF4, time

                    ThisDataset = netCDF4.Dataset(filename + ".nc", "w", 
                                                  format = "NETCDF4")
                    
                    ThisDataset.created = "Created " + time.ctime(time.time())
                    ThisDataset.Source = "ARTDECO V1.1"
                    ThisDataset.URL = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisDataset.case_name = self.Output_subdirectory
                    ThisDataset.mode = self.Mode
                    if self.Mode == "kdis":
                        ThisDataset.k_distribution = self.K_distributions_parameterization
                        filter_list = ""
                        if len(self.Instrumental_filters_list) > 0:
                            filter_list += self.Instrumental_filters_list[0]
                            for iff in range(len(self.Instrumental_filters_list) - 1):
                                filter_list += " " + self.Instrumental_filters_list[iff + 1]
                        else:
                            filter_list = "none"
                        ThisDataset.instrumental_filter = filter_list
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisDataset.gas_list = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisDataset.gas_continuum_list = gas_continuum_list_str
                    ThisDataset.radiative_transfer_model = self.RTE_solver
                    ThisDataset.stokes_components = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisDataset.depolarization_value = str(float(depol_value))
                    ThisDataset.atmosphere = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisDataset.phase_function_truncation = self.Truncation_method
                    ThisDataset.particles = particles_string
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisDataset.Fbeam = str(float(self.FBeam[0]))
                        else:
                            ThisDataset.Fbeam = "1.0"
                    if self.Mode == "kdis":
                        ThisDataset.incident_spectrum = self.Incident_spectrum                      

                    Nw_dim = ThisDataset.createDimension("Nw", Nw)
                    Nsza_dim = ThisDataset.createDimension("Nsza", Nsza)
                    Nalt_dim = ThisDataset.createDimension("Nalt", Nalt)

                    Wv_var = ThisDataset.createVariable("Wv", "f8", ("Nw",))
                    Alt_var = ThisDataset.createVariable("Altitude", "f8", ("Nalt",)) 
                    Sza_var = ThisDataset.createVariable("Sza", "f8", ("Nsza",))

                    Wv_var.long_name = "Wavelengths"
                    Wv_var.units = "micrometers"
                    Alt_var.long_name = "Altitude of the atmospheric layer interface"
                    Alt_var.units = "kilometers"
                    Sza_var.long_name = "Solar zenith angles"
                    Sza_var.units = "degrees"

                    Wv_var[:] = Wv
                    Alt_var[:] = Alt
                    Sza_var[:] = Sza
                    
                    Pres_var = ThisDataset.createVariable("Pres", "f8", ("Nalt",))
                    Pres_var.long_name = "Air pressure at the given altitude"
                    Pres_var.units = "hPa"
                    Pres_var[:] = Pres

                    Flux_dw_var = ThisDataset.createVariable("flux_dw", "f8", ("Nw", "Nsza", "Nalt",))
                    Flux_dw_var.long_name = "Total downward flux"
                    Flux_dw_var.units = "W/m2"
                    Flux_dw_var[:] = flux_dw

                    Flux_dir_var = ThisDataset.createVariable("flux_dir", "f8", ("Nw", "Nsza", "Nalt",))
                    Flux_dir_var.long_name = "Directional downward flux"
                    Flux_dir_var.units = "W/m2"
                    Flux_dir_var[:] = flux_dir

                    Flux_up_var = ThisDataset.createVariable("flux_up", "f8", ("Nw", "Nsza", "Nalt",))
                    Flux_up_var.long_name = "Total upward flux"
                    Flux_up_var.units = "W/m2"
                    Flux_up_var[:] = flux_up
                    
                    ThisDataset.close()

                if format == "hdf5" or format == "all":
                    import h5py
                    ThisFile = h5py.File(filename + ".hdf5", "w")
                    ThisFile.attrs['created'] = "Created " + time.ctime(time.time())
                    ThisFile.attrs['Source'] = "ARTDECO V1.1"
                    ThisFile.attrs['URL'] = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisFile.attrs['case_name'] = self.Output_subdirectory
                    ThisFile.attrs['mode'] = self.Mode
                    if self.Mode == "kdis":
                        ThisFile.attrs['k_distribution'] = self.K_distributions_parameterization
                        ThisFile.attrs['wv_min'] = self.Wavelength_start
                        ThisFile.attrs['wv_max'] = self.Wavelength_end
                        filter_list = ""
                        if len(self.Instrumental_filters_list) > 0:
                            filter_list += self.Instrumental_filters_list[0]
                            for iff in range(len(self.Instrumental_filters_list) - 1):
                                filter_list += " " + self.Instrumental_filters_list[iff + 1]
                        else:
                            filter_list = "none"
                        ThisFile.attrs['instrumental_filter'] = filter_list
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisFile.attrs['gas_list'] = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisFile.attrs['gas_continuum_list'] = gas_continuum_list_str
                    ThisFile.attrs['radiative_transfer_model'] = self.RTE_solver
                    ThisFile.attrs['stokes_components'] = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisFile.attrs['depolarization_value'] = str(float(depol_value))
                    ThisFile.attrs['atmosphere'] = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisFile.attrs['phase_function_truncation'] = self.Truncation_method
                    ThisFile.attrs['particles'] = particles_string
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisFile.attrs['Fbeam'] = str(float(self.FBeam[0]))
                        else:
                            ThisFile.attrs['Fbeam'] = "1.0"
                    if self.Mode == "kdis":
                        ThisFile.attrs['incident_spectrum'] = self.Incident_spectrum  
                    ThisFile.attrs['Nalt'] = Nalt
                    ThisFile.attrs['Nsza'] = Nsza
                    ThisFile.attrs['Nw'] = Nw

                    Wv_dset = ThisFile.create_dataset("Wv", data = Wv, dtype = "f8")
                    Alt_dset = ThisFile.create_dataset("Altitude", data = Alt, dtype = "f8")
                    Sza_dset = ThisFile.create_dataset("Sza", data = Sza, dtype = "f8")

                    Wv_dset.attrs['long_name'] = "Wavelengths"
                    Wv_dset.attrs['units'] = "micrometers"
                    Wv_dset.attrs['dimensions'] = "(Nw,)"
                    Alt_dset.attrs['long_name'] = "Altitude of the atmospheric layer interface"
                    Alt_dset.attrs['units'] = "kilometers"
                    Alt_dset.attrs['dimensions'] = "(Nalt,)"
                    Sza_dset.attrs['long_name'] = "Solar zenith angles"
                    Sza_dset.attrs['units'] = "degrees"
                    Sza_dset.attrs['dimensions'] = "(Nsza,)"
                    
                    Pres_dset = ThisFile.create_dataset("Pres", data = Pres, dtype = "f8")
                    Pres_dset.attrs['long_name'] = "Air pressure at the given altitude"
                    Pres_dset.attrs['units'] = "hPa"
                    Pres_dset.attrs['dimensions'] = "(Nalt,)"

                    Flux_dw_dset = ThisFile.create_dataset("flux_dw", data=flux_dw, dtype = "f8") 
                    Flux_dw_dset.attrs['dimensions'] = "(Nw, Nsza, Nalt,)"
                    Flux_dw_dset.attrs['long_name'] = "Total downward flux"
                    Flux_dw_dset.attrs['units'] = "W/m2"

                    Flux_dir_dset = ThisFile.create_dataset("flux_dir", data=flux_dir, dtype = "f8")
                    Flux_dir_dset.attrs['dimensions'] = "(Nw, Nsza, Nalt,)"
                    Flux_dir_dset.attrs['long_name'] = "Directional downward flux"
                    Flux_dir_dset.attrs['units'] = "W/m2"

                    Flux_up_dset = ThisFile.create_dataset("flux_up", data=flux_up, dtype = "f8")
                    Flux_up_dset.attrs['dimensions'] = "(Nw, Nsza, Nalt,)"
                    Flux_up_dset.attrs['long_name'] = "Total upward flux"
                    Flux_up_dset.attrs['units'] = "W/m2"

                    ThisFile.close()
                    del Sza_dset, Alt_dset, Wv_dset, Pres_dset, Flux_dw_dset, Flux_up_dset, Flux_dir_dset

            ## Read warming rates.
            for in_file in glob.glob(output_dir + "Warming_rate.dat"):
                filename = in_file.split(".da")[0]

                ## Reads ascii file and erase comment lines.
                file_in = open(in_file, "r")
                raw_lines_in = file_in.readlines()
                file_in.close()
                lines_in = clean_comment_lines(raw_lines_in)
                del raw_lines_in
                
                ## dimensions.
                Nsza = int(clean_line(lines_in[0])[0])
                Nalt = int(clean_line(lines_in[1])[0])

                ## dimension variables.
                Sza = numpy.zeros([Nsza], dtype = "f8")

                ## variables.
                up = numpy.zeros([Nalt], dtype = "f8")
                low = numpy.zeros([Nalt], dtype = "f8")
                rate = numpy.zeros([Nsza, Nalt], dtype = "f8")
                
                for isza in range(Nsza):
                    Sza[isza] = float(clean_line(lines_in[isza * (Nalt + 1) + 2])[0])
                for ia in range(Nalt):
                    up[ia] = float(clean_line(lines_in[ia + 3])[1])
                    low[ia] = float(clean_line(lines_in[ia + 3])[0])
                    for isza in range(Nsza):
                        rate[isza, ia] = float(clean_line(lines_in[isza * (Nalt + 1) + ia + 3])[2])
                           
                if format == "netcdf" or format == "all":
                    import netCDF4, time

                    ThisDataset = netCDF4.Dataset(filename + ".nc", "w", 
                                                  format = "NETCDF4")
                    
                    ThisDataset.created = "Created " + time.ctime(time.time())
                    ThisDataset.Source = "ARTDECO V1.1"
                    ThisDataset.URL = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisDataset.case_name = self.Output_subdirectory
                    ThisDataset.mode = self.Mode
                    if self.Mode == "kdis":
                        ThisDataset.k_distribution = self.K_distributions_parameterization
                        ThisDataset.wv_min = self.Wavelength_start
                        ThisDataset.wv_max = self.Wavelength_end
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisDataset.gas_list = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisDataset.gas_continuum_list = gas_continuum_list_str
                    ThisDataset.radiative_transfer_model = self.RTE_solver
                    ThisDataset.stokes_components = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisDataset.depolarization_value = str(float(depol_value))
                    ThisDataset.atmosphere = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisDataset.phase_function_truncation = self.Truncation_method
                    ThisDataset.particles = particles_string
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisDataset.Fbeam = str(float(self.FBeam[0]))
                        else:
                            ThisDataset.Fbeam = "1.0"
                    if self.Mode == "kdis":
                        ThisDataset.incident_spectrum = self.Incident_spectrum                      

                    Nsza_dim = ThisDataset.createDimension("Nsza", Nsza)
                    Nalt_dim = ThisDataset.createDimension("Nalt", Nalt)

                    Sza_var = ThisDataset.createVariable("Sza", "f8", ("Nsza",))
                    Sza_var.long_name = "Solar zenith angles"
                    Sza_var.units = "degrees"
                    Sza_var[:] = Sza

                    up_var = ThisDataset.createVariable("Upper_bound", "f8", ("Nalt",))
                    up_var.long_name = "Altitude of the atmospheric layer upper interface"
                    up_var.units = "kilometers"
                    up_var[:] = up

                    low_var = ThisDataset.createVariable("Lower_bound", "f8", ("Nalt",))
                    low_var.long_name = "Altitude of the atmospheric layer lower interface"
                    low_var.units = "kilometers"
                    low_var[:] = low

                    rate_var = ThisDataset.createVariable("Warming_rate", "f8", ("Nsza", "Nalt",))
                    rate_var.long_name = "Warming rate"
                    rate_var.units = "K/day"
                    rate_var[:] = rate
                    
                    ThisDataset.close()

                if format == "hdf5" or format == "all":
                    import h5py
                    ThisFile = h5py.File(filename + ".hdf5", "w")
                    ThisFile.attrs['created'] = "Created " + time.ctime(time.time())
                    ThisFile.attrs['Source'] = "ARTDECO V1.1"
                    ThisFile.attrs['URL'] = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisFile.attrs['case_name'] = self.Output_subdirectory
                    ThisFile.attrs['mode'] = self.Mode
                    if self.Mode == "kdis":
                        ThisFile.attrs['k_distribution'] = self.K_distributions_parameterization
                        ThisFile.attrs['wv_min'] = self.Wavelength_start
                        ThisFile.attrs['wv_max'] = self.Wavelength_end
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisFile.attrs['gas_list'] = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisFile.attrs['gas_continuum_list'] = gas_continuum_list_str
                    ThisFile.attrs['radiative_transfer_model'] = self.RTE_solver
                    ThisFile.attrs['stokes_components'] = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisFile.attrs['depolarization_value'] = str(float(depol_value))
                    ThisFile.attrs['atmosphere'] = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisFile.attrs['phase_function_truncation'] = self.Truncation_method
                    ThisFile.attrs['particles'] = particles_string
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisFile.attrs['Fbeam'] = str(float(self.FBeam[0]))
                        else:
                            ThisFile.attrs['Fbeam'] = "1.0"
                    if self.Mode == "kdis":
                        ThisFile.attrs['incident_spectrum'] = self.Incident_spectrum 
                    ThisFile.attrs['Nsza'] = Nsza
                    ThisFile.attrs['Nalt'] = Nalt

                    Sza_dset = ThisFile.create_dataset("Sza", data = Sza, dtype = "f8")
                    Sza_dset.attrs['long_name'] = "Solar zenith angles"
                    Sza_dset.attrs['dimensions'] = "(Sza,)"
                    Sza_dset.attrs['units'] = "degrees"

                    up_dset = ThisFile.create_dataset("Upper_bound", data = up, dtype = "f8")
                    up_dset.attrs['long_name'] = "Altitude of the atmospheric layer upper interface"
                    up_dset.attrs['dimensions'] = "(Nalt,)"
                    up_dset.attrs['units'] = "kilometers"

                    low_dset = ThisFile.create_dataset("Lower_bound", data = low, dtype = "f8")
                    low_dset.attrs['long_name'] = "Altitude of the atmospheric layer lower interface"
                    low_dset.attrs['dimensions'] = "(Nalt,)"
                    low_dset.attrs['units'] = "kilometers"

                    rate_dset = ThisFile.create_dataset("Net_radiative_heating_rate", data = rate, dtype = "f8")
                    rate_dset.attrs['dimensions'] = "(Sza, Nalt,)"
                    rate_dset.attrs['long_name'] = "Warming rate"
                    rate_dset.attrs['units'] = "K/day"
                    ThisFile.close()
                    del Sza_dset, up_dset, low_dset, rate_dset

            ## Read net radiative heating rates.
            for in_file in glob.glob(output_dir + "Net_radiative_heating_rate.dat"):
                filename = in_file.split(".da")[0]

                ## Reads ascii file and erase comment lines.
                file_in = open(in_file, "r")
                raw_lines_in = file_in.readlines()
                file_in.close()
                lines_in = clean_comment_lines(raw_lines_in)
                del raw_lines_in
                
                ## dimensions.
                Nsza = int(clean_line(lines_in[0])[0])
                Nalt = int(clean_line(lines_in[1])[0])

                ## dimension variables.
                Sza = numpy.zeros([Nsza], dtype = "f8")

                ## variables.
                up = numpy.zeros([Nalt], dtype = "f8")
                low = numpy.zeros([Nalt], dtype = "f8")
                rate = numpy.zeros([Nsza, Nalt], dtype = "f8")
                
                for isza in range(Nsza):
                    Sza[isza] = float(clean_line(lines_in[isza * (Nalt + 1) + 2])[0])
                for ia in range(Nalt):
                    up[ia] = float(clean_line(lines_in[ia + 3])[1])
                    low[ia] = float(clean_line(lines_in[ia + 3])[0])
                    for isza in range(Nsza):
                        rate[isza, ia] = float(clean_line(lines_in[isza * (Nalt + 1) + ia + 3])[2])
                           
                if format == "netcdf" or format == "all":
                    import netCDF4, time

                    ThisDataset = netCDF4.Dataset(filename + ".nc", "w", 
                                                  format = "NETCDF4")
                    
                    ThisDataset.created = "Created " + time.ctime(time.time())
                    ThisDataset.Source = "ARTDECO V1.1"
                    ThisDataset.URL = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisDataset.case_name = self.Output_subdirectory
                    ThisDataset.mode = self.Mode
                    if self.Mode == "kdis":
                        ThisDataset.k_distribution = self.K_distributions_parameterization
                        ThisDataset.wv_min = self.Wavelength_start
                        ThisDataset.wv_max = self.Wavelength_end
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisDataset.gas_list = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisDataset.gas_continuum_list = gas_continuum_list_str
                    ThisDataset.radiative_transfer_model = self.RTE_solver
                    ThisDataset.stokes_components = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisDataset.depolarization_value = str(float(depol_value))
                    ThisDataset.atmosphere = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisDataset.phase_function_truncation = self.Truncation_method
                    ThisDataset.particles = particles_string
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisDataset.Fbeam = str(float(self.FBeam[0]))
                        else:
                            ThisDataset.Fbeam = "1.0"
                    if self.Mode == "kdis":
                        ThisDataset.incident_spectrum = self.Incident_spectrum                      

                    Nsza_dim = ThisDataset.createDimension("Nsza", Nsza)
                    Nalt_dim = ThisDataset.createDimension("Nalt", Nalt)

                    Sza_var = ThisDataset.createVariable("Sza", "f8", ("Nsza",))
                    Sza_var.long_name = "Solar zenith angles"
                    Sza_var.units = "degrees"
                    Sza_var[:] = Sza

                    up_var = ThisDataset.createVariable("Upper_bound", "f8", ("Nalt",))
                    up_var.long_name = "Altitude of the atmospheric layer upper interface"
                    up_var.units = "kilometers"
                    up_var[:] = up

                    low_var = ThisDataset.createVariable("Lower_bound", "f8", ("Nalt",))
                    low_var.long_name = "Altitude of the atmospheric layer lower interface"
                    low_var.units = "kilometers"
                    low_var[:] = low

                    rate_var = ThisDataset.createVariable("Net_radiative_heating_rate", "f8", ("Nsza", "Nalt",))
                    rate_var.long_name = "Net radiative heating rate"
                    rate_var.units = "W/m2/m"
                    rate_var[:] = rate
                    
                    ThisDataset.close()

                if format == "hdf5" or format == "all":
                    import h5py
                    ThisFile = h5py.File(filename + ".hdf5", "w")
                    ThisFile.attrs['created'] = "Created " + time.ctime(time.time())
                    ThisFile.attrs['Source'] = "ARTDECO V1.1"
                    ThisFile.attrs['URL'] = "http://www.icare.univ-lille1.fr/projects/artdeco"
                    ThisFile.attrs['case_name'] = self.Output_subdirectory
                    ThisFile.attrs['mode'] = self.Mode
                    if self.Mode == "kdis":
                        ThisFile.attrs['k_distribution'] = self.K_distributions_parameterization
                        ThisFile.attrs['wv_min'] = self.Wavelength_start
                        ThisFile.attrs['wv_max'] = self.Wavelength_end
                        gas_list_str = ""
                        if len(self.Gas_list) > 0:
                            gas_list_str += self.Gas_list[0]
                            for ig in range(len(self.Gas_list) - 1):
                                gas_list_str += " " + self.Gas_list[ig + 1]
                        ThisFile.attrs['gas_list'] = gas_list_str
                        gas_continuum_list_str = ""
                        if len(self.Gas_continuum_list) > 0:
                            gas_continuum_list_str += self.Gas_continuum_list[0]
                            for ig in range(len(self.Gas_continuum_list) - 1):
                                gas_continuum_list_str += " " + self.Gas_continuum_list[ig + 1]
                        ThisFile.attrs['gas_continuum_list'] = gas_continuum_list_str
                    ThisFile.attrs['radiative_transfer_model'] = self.RTE_solver
                    ThisFile.attrs['stokes_components'] = str(int(self.Stokes_components))
                    if self.Depolarization_value == -1:
                        depol_value = 0.0279
                    else:
                        depol_value = self.Depolarization_value
                    ThisFile.attrs['depolarization_value'] = str(float(depol_value))
                    ThisFile.attrs['atmosphere'] = "none" * (not(self.With_atmosphere)) \
                        + self.Atmospheric_profile_original * self.With_atmosphere
                    particles_string = ""
                    if self.With_particles:
                        particles_string += self.Particle_option[0][0]
                        for ip in range(len(self.Particle_option) - 1):
                            particles_string += " " + self.Particle_option[ip + 1][0]
                    else:
                        particles_string = "none" * (not(self.With_particles))
                    ThisFile.attrs['phase_function_truncation'] = self.Truncation_method
                    ThisFile.attrs['particles'] = particles_string
                    if self.Mode == "mono":
                        if len(self.FBeam) == 1:
                            ThisFile.attrs['Fbeam'] = str(float(self.FBeam[0]))
                        else:
                            ThisFile.attrs['Fbeam'] = "1.0"
                    if self.Mode == "kdis":
                        ThisFile.attrs['incident_spectrum'] = self.Incident_spectrum 
                    ThisFile.attrs['Nsza'] = Nsza
                    ThisFile.attrs['Nalt'] = Nalt

                    Sza_dset = ThisFile.create_dataset("Sza", data = Sza, dtype = "f8")
                    Sza_dset.attrs['long_name'] = "Solar zenith angles"
                    Sza_dset.attrs['dimensions'] = "(Sza,)"
                    Sza_dset.attrs['units'] = "degrees"

                    up_dset = ThisFile.create_dataset("Upper_bound", data = up, dtype = "f8")
                    up_dset.attrs['long_name'] = "Altitude of the atmospheric layer upper interface"
                    up_dset.attrs['dimensions'] = "(Nalt,)"
                    up_dset.attrs['units'] = "kilometers"

                    low_dset = ThisFile.create_dataset("Lower_bound", data = low, dtype = "f8")
                    low_dset.attrs['long_name'] = "Altitude of the atmospheric layer lower interface"
                    low_dset.attrs['dimensions'] = "(Nalt,)"
                    low_dset.attrs['units'] = "kilometers"

                    rate_dset = ThisFile.create_dataset("Net_radiative_heating_rate", data = rate, dtype = "f8")
                    rate_dset.attrs['dimensions'] = "(Sza, Nalt,)"
                    rate_dset.attrs['long_name'] = "Net radiative heating rate"
                    rate_dset.attrs['units'] = "W/m2/m"
                    ThisFile.close()
                    del Sza_dset, up_dset, low_dset, rate_dset
            
            

    ## Checks for configuration errors.
    def ChecksConfigurationErrors(self):
        ## module netcdf and hdf5 if needed.
        if self.Output_format == "all" \
                or self.Output_format == "netcdf" \
                or self.Output_format == "netCDF" \
                or self.Output_format == "NetCDF":
            try:
                import netCDF4
            except:
                mesg = "main.py : netCDF4 python module can not be loaded."
                mesg += " For this reason, \""
                mesg += self.Output_format
                mesg += "\" output format cannot be chosen in main.py.\n"
                RaiseError(self.Log_file, mesg)

        if self.Output_format == "all" \
                or self.Output_format == "hdf5" \
                or self.Output_format == "HDF5":
            try:
                import h5py
            except:
                mesg = "main.py : h5py python module can not be loaded."
                mesg += " For this reason, \""
                mesg += self.Output_format
                mesg += "\" output format cannot be chosen in main.py.\n"
                RaiseError(self.Log_file, mesg)

        if self.Output_format != "ascii" \
                and self.Output_format != "netcdf" \
                and self.Output_format != "netCDF" \
                and self.Output_format != "NetCDF" \
                and self.Output_format != "hdf5" \
                and self.Output_format != "HDF5" \
                and self.Output_format != "all":
            mesg = "main.py : Output_format \""
            mesg += self.Output_format
            mesg += "\" is not recognized. It should be \"ascii\","
            mesg += "\"netcdf\", \"hdf5\" or \"all\"."
            RaiseError(self.Log_file, mesg)
        
        ## Arguments not valid.
        if self.RTE_solver != "disort" \
                and self.RTE_solver != "doad" \
                and self.RTE_solver != "mcrad1d" \
                and self.RTE_solver != "sinsca" \
                and self.RTE_solver != "none":
            mesg = "main.py : the specified RT model (\""
            mesg += self.RTE_solver + "\") is not recognized. Please check "
            mesg += "the flag RTE_solver in " + "main.py"
            RaiseError(self.Log_file, mesg)

        if self.Mode != "mono" and self.Mode != "kdis":
            mesg = "main.py : the specified mode (\""
            mesg += self.Mode + "\") is not recognized. It should be"
            mesg += " either \"kdis\" or \"mono\"."
            RaiseError(self.Log_file, mesg)

        if self.Mode == "mono":
            if len(self.Wavelength_list) == 0:
                mesg = "mono.py : At least one wavelength must be specified "
                mesg += "in Wavelength_list."
                RaiseError(self.Log_file, mesg)
            for iw in range(len(self.Wavelength_list)):
                if self.Wavelength_list[iw] < 0.:
                    mesg = "mono.py : Wavelengths specified in Wavelength_list "
                    mesg += "must be positive."
                    RaiseError(self.Log_file, mesg)

        if self.Truncation_method != "none" \
                and self.Truncation_method != "dm" \
                and self.Truncation_method != "dfit" \
                and self.Truncation_method != "potter":
            mesg = "particles.py : the specified truncation method (\""
            mesg += self.Truncation_method + "\") is not recognized. It"
            mesg += " should be either \"dm\", \"dfit\", \"potter\" or \"none\"."
            RaiseError(self.Log_file, mesg)
        
        if self.Wavelength_unit != "cm-1" and self.Wavelength_unit != "cm^-1" \
                and self.Wavelength_unit != "microm" \
                and self.Wavelength_unit != "micrometers" \
                and self.Wavelength_unit != "mum":
            mesg = self.Mode + " : the specified unit (\""
            mesg += self.Wavelength_unit + "\") for the wavelengths "
            mesg += "is not valid."
            RaiseError(self.Log_file, mesg)

        if self.Depolarization_value != -1 \
                and not(isinstance(self.Depolarization_value, float)) \
                and not(isinstance(self.Depolarization_value, int)):
            mesg = self.Mode = " : Depolarization_value is not a valid value."
            RaiseError(self.Log_file, mesg)            

        ## Kdis parameters.
        if self.Mode == "kdis":
            input_file = "./lib/kdis/reference/kdis_"
            input_file += self.K_distributions_parameterization
            input_file += "_def.dat"
            if not(os.path.lexists(input_file)):
                mesg = "kdis.py : The specified K-distribution"
                mesg += " parameterization is \""
                mesg += self.K_distributions_parameterization
                mesg += "\", but the file \"" + input_file
                mesg += "\" is not a valid file."
                RaiseError(self.Log_file, mesg)

            for gas in self.Gas_list:
                input_file = "./lib/kdis/reference/kdis_"
                input_file += self.K_distributions_parameterization
                input_file += "_" + gas + ".dat"
                if not(os.path.lexists(input_file)):
                    mesg = "kdis.py : The specified K-distribution"
                    mesg += " parameterization is \""
                    mesg += self.K_distributions_parameterization
                    mesg += "\", but the file \"" + input_file
                    mesg += "\" giving parameterization for gas \""
                    mesg += gas + "\" is not a valid file. "
                    mesg += "Please check the Gas_list in kdis.py"
                    RaiseError(self.Log_file, mesg)

            for gc in self.Gas_continuum_list:
                if not(gc in self.Gas_list):
                    mesg = "kdis.py : The gas \"" + str(gc) + "\" is in the "
                    mesg += "Gas_continuum_list, but not in the Gas_list."
                    RaiseError(self.Log_file, mesg)

        if self.RTE_solver == "none" and self.Legendre_and_truncation \
                and self.Number_of_betal == -1:
            mesg = "main.py | particles.py : If no RT solver is select"
            mesg += "ed (RTE_solver = \"none\") and Legendre_and_truncation = "
            mesg += "True, you must set With_particles = True and a float val"
            mesg += "ue for Number_of_betal (e.g not \"optimized\" or -1)."
            RaiseError(self.Log_file, mesg)

        ## disort parameters.
        if self.RTE_solver == "disort":
            if self.Truncation_method == "potter" \
                    and self.TMS_correction:
                mesg = "main.py | disort.py : It is not possible to "
                mesg += "use in the same calculation:\n     - disort RTE "
                mesg += "solver,\n     - Potter truncation,\n     - TMS "
                mesg += "correction.\n"
                RaiseError(self.Log_file, mesg)

            if not(self.Compute_radiance) and self.Print_downward_radiance:
                mesg = "disort.py : Mismatch between two arguments : "
                mesg += "Compute_radiance and Print_downward_radiance."
                RaiseError(self.Log_file, mesg)

            if self.Convergence_criterion > 0.001:
                mesg = "disort.py : The value for the convergence criterion ("
                mesg += str(Convergence_criterion) + ") might be too high. "
                mesg += "It should be between 0.000 and 0.001."
                RaiseError(self.Log_file, mesg)
            
        ## Monochromatic mode parameters.
        if self.Thermal_in_situ and self.Mode == "mono":
            mesg = "main.py | disort.py : Thermal local emission "
            mesg += "(keyword \"Thermal_emission\" in disort.py) can "
            mesg += "not be accounted in monochromatic mode, but only in "
            mesg += "k-distribution mode."
            RaiseError(self.Log_file, mesg)

        ## TMS correction parameters.
        if self.RTE_solver != "none":
            if self.TMS_correction and self.Number_of_betal != -1:
                mesg = "particles.py : When using TMS_correction, Number"
                mesg += "_of_betal must be equal to -1 (optimized value)."
                RaiseError(self.Log_file, mesg)
        if self.TMS_correction and self.Print_downward_radiance:
            mesg = self.RTE_solver + ".py | particles.py : Downward radiances can "
            mesg += "not be computed when using TMS_correction."
            RaiseError(self.Log_file, mesg)
                
        ## mcrad1d parameters.
        if self.RTE_solver == "mcrad1d" and self.TMS_correction:
            mesg = "main.py | particles.py : No intensity TMS_corr"
            mesg += "ection is performed when using \"mcrad1d\" RTE solver."
            RaiseError(self.Log_file, mesg)

        if self.RTE_solver == "mcrad1d" and self.Mode == "kdis":
            mesg = "main.py : K-distribution mode is not implemented "
            mesg += "with the mcrad1d RTE solver."
            RaiseError(self.Log_file, mesg)

        ## addoub parameters.
        if self.RTE_solver == "addoub" and self.TMS_correction:
            mesg = "main.py | particles.py : No intensity TMS_corr"
            mesg += "ection is performed when using \"addoub\" RTE solver."
            RaiseError(self.Log_file, mesg)
            
        ## Polarization parameters (Stokes vector components).
        if self.Stokes_components != 1 \
                and self.Stokes_components != 3 \
                and self.Stokes_components != 4:
            mesg = "main.py : The number of Stokes parameters (\"Stoke"
            mesg += "s_components\" = " + self.Stokes_components + " is not r"
            mesg += "ecognized. It should be either 1, 3,  or 4."
            RaiseError(self.Log_file, mesg)

        if self.Stokes_components != 1 and self.RTE_solver == "disort":
            mesg = "main.py : When calling disort RTE solver, the numb"
            mesg += "er of Stokes parameters should be equal to 1 ("
            mesg += str(self.Stokes_components) + " specified)."
            RaiseError(self.Log_file, mesg)

        ## Atmosphere parameters.
        if self.With_atmosphere:
            if self.Atmosphere_in_ppmv \
                    and self.Atmospheric_profile.split("_")[-1] != "ppmv":
                mesg = "atmosphere.py : The atmosphere definition file in "
                mesg += "ppmv (as specified in atmosphere.py) must contain "
                mesg += "the suffix _ppmv (currently \""
                mesg += self.Atmospheric_profile + "\")." 
                RaiseError(self.Log_file, mesg)

            input_file = "./lib/atm/reference/atm_" + self.Atmospheric_profile + ".dat"
            if not(os.path.lexists(input_file)):
                mesg = "atmosphere.py : The file \"" + input_file + "\" neces"
                mesg += "sary for the specified atmospheric profile \""
                mesg += self.Atmospheric_profile + "\" is not a valid file."
                RaiseError(self.Log_file, mesg)          
            
        ## Particles parameters.
        if self.With_particles:
         
            if not(isinstance(self.Number_of_betal, float)) \
                    and not(isinstance(self.Number_of_betal, int)):
                mesg = "particles.py : Number_of_betal (\""
                mesg += self.Number_of_betal + "\") is not valid."
                RaiseError(self.Log_file, mesg)

            if self.RTE_solver == "none" and self.Number_of_betal == -1:
                mesg = "main.py | particles.py: When using RTE_solver "
                mesg += "= \"none\" in main.py, Number_of_betal must be "
                mesg += "specified in particles.py (not -1, nor \"optimized\")."
                RaiseError(self.Log_file, mesg)

            if (self.RTE_solver == "mcrad1d" or self.RTE_solver == "sinsca") \
                    and self.Truncation_method == "none":
                if self.Number_of_betal == -1:
                    mesg = "main.py | particles.py : When using \"" 
                    mesg += self.RTE_solver
                    mesg += "\" RTE solver, Number_of_betal must be specified "
                    mesg += "in particles.py (not -1)."
                    RaiseError(self.Log_file, mesg)

            if self.Number_of_betal != -1 and  self.Number_of_betal < 2:
                mesg = "particles.py : Number_of_betal must be -1 or greater than 2."
                RaiseError(self.Log_file, mesg)

            if len(self.Particle_option) == 0:
                mesg = "main.py | particles.py : The particle list is "
                mesg += "empty in particles.py. If no particle are "
                mesg += "considered in the calculation, the keyword "
                mesg += "\"With_particles\" must be set to \"False\" in "
                mesg += "main.py."
                RaiseError(self.Log_file, mesg)

            for i1 in range(len(self.Particle_option)):
                for i2 in range(i1):
                    if self.Particle_option[i1][0] == self.Particle_option[i2][0]:
                        mesg = "particles.py : Problem with particle list. The "
                        mesg += "particle \"" + self.Particle_option[i1][0]
                        mesg += "\" appears more that once."
                        RaiseError(self.Log_file, mesg)
            
            hg_flag = False
            for particle in self.Particle_option:
                if particle[2]:
                    hg_flag = True
            if self.Stokes_components > 1 and hg_flag:
                mesg = "main.py | particles.py : Henyey-Greenstein app"
                mesg += "roximation is activated for the particle "
                mesg += particle[0] + " in particles.py, but is not allowed w"
                mesg += "ith polarized radiative transfer (Stokes_component"
                mesg += "s = " + self.Stokes_components + " in " + "main.py"
                RaiseError(self.Log_file, mesg)

            for particle in self.Particle_option:
                if particle[1] and particle[2]:
                    mesg = "particles.py : \"use_betal\" and \"interp\" can"
                    mesg += "not be used at the same time (concerned "
                    mesg += "particle type : \"" + particle[0] + "\")."
                    RaiseError(self.Log_file, mesg)
                if particle[1] and self.TMS_correction:
                    mesg = "particles.py : \"user_betal\" and \"TMS_correct"
                    mesg += "ion\" can not be used at the same time (concer"
                    mesg += "ned particle type : \"" + particle[0] + "\")."
                    RaiseError(self.Log_file, mesg)
                if particle[1] and particle[3]:
                    mesg = "particles.py : \"use_betal\" and \"H-G\" can"
                    mesg += "not be used at the same time (concerned "
                    mesg += "particle type : \"" + particle[0] + "\")."
                    RaiseError(self.Log_file, mesg)
                if particle[5] == "layer" or particle[5] == "sh":
                    if len(particle) != 7:
                        mesg = "particles.py : For the particle type \""
                        mesg += particle[0] + "\" the specified "
                        mesg += "parameterization for the vertical "
                        mesg += "distribution is \"" + particle[5]
                        mesg += "\", hence the number of parameters "
                        mesg += "should be equal to 1 ("
                        mesg += str(len(particle) - 6) + " specified"
                        mesg += ")."
                        RaiseError(self.Log_file, mesg)
                elif particle[5] == "gauss":
                    if len(particle) != 8:
                        mesg = "particles.py : For the particle type \""
                        mesg += particle[0] + "\" the specified "
                        mesg += "parameterization for the vertical "
                        mesg += "distribution is"
                        mesg += "\"gauss\", hence the "
                        mesg += "number of parameters "
                        mesg += "should be equal to 2 ("
                        mesg += str(len(particle) - 6) + " specified"
                        mesg += ")."
                        RaiseError(self.Log_file, mesg)
                elif particle[5] == "user":
                    if len(particle) != 6:
                        mesg = "particles.py : For the particle type \""
                        mesg += particle[0] + "\" the specified "
                        mesg += "parameterization for the vertical "
                        mesg += "distribution "
                        mesg += "is \"" + particle[5] + "\", hence "
                        mesg += "the number of parameters "
                        mesg += "should be equal to 0 ("
                        mesg += str(len(particle) - 6) + " specified"
                        mesg += ") and the distribution is read in "
                        mesg += "the file ./lib/vdist/vdist_"
                        mesg += particle[0] + ".dat."
                        RaiseError(self.Log_file, mesg)
                    else:
                        if not(os.path.lexists(self.Lib_directory + \
                                                   "/vdist/vdist_" + \
                                                   particle[0] + ".dat")):
                            mesg = "particles.py : The file "
                            mesg += self.Lib_directory
                            mesg += "/vdist/vdist_" + particle[0] + ".dat,"
                            mesg += " necessary to specify user-defined "
                            mesg += "vertical distribution of particle \""
                            mesg += particle[0] + "\", is not a valid file."
                            RaiseError(self.Log_file, mesg)
                else:
                    mesg = "particles.py : The specified type of vertical "
                    mesg += "distribution parameterization (\"" + particle[5]
                    mesg += "\") is not recognized for the particle "
                    mesg += "type \"" + particle[0] + "\"."
                    RaiseError(self.Log_file, mesg)

        ## Instrumental filter parameters.
        for filt in self.Instrumental_filters_list:
            if not(os.path.lexists("./lib/filter/filter_" + filt + ".dat")):
                mesg = "kdis.py : The file ./lib/filter/filter_"
                mesg += str(filt) + ".dat, necessary for the instrumental filter"
                mesg += " \"" + str(filt) + "\" (as specified in kdis.py) is not"
                mesg += " a valid file."
                RaiseError(self.Log_file, mesg)
       
        if self.RTE_solver != "none":

            ## Fbeam parameters.
            if len(self.FBeam) != 0 and self.Mode == "kdis":
                mesg = "main.py | geometrics.py : the use of FBeam (as "
                mesg += "specified in geometrics.py) is not compatible"
                mesg += " with the k-distribution mode.\n             "
                mesg += "In kdis mode, \"Fbeam_user = "
                mesg += "False\" and the incident radiation is co"
                mesg += "mputed from the spectrum specified in ge"
                mesg += "ometrics.py (currently \"" + self.Incident_spectrum + "\")."
                RaiseError(self.Log_file, mesg)

            ## Surface parameters
            if self.Surface_parameters[0] == "lambert":
                if self.Surface_parameters[1] != "cste":
                    input_file = "./lib/surface/surfalb_"
                    input_file += self.Surface_parameters[1]
                    input_file += ".dat"
                    if not(os.path.lexists(input_file)):
                        mesg = "surface.py : The specified surface type is"
                        mesg += "\"" + self.Surface_parameters[1] + "\""
                        mesg += ", but the file \"" + input_file
                        mesg +=  "\" is not a valid file."
                        RaiseError(self.Log_file, mesg)
                if self.Surface_parameters[1] == "ocean":
                    mesg = "surface.py : The specified surface type is "
                    mesg += "\"ocean\" but the surface nature is "
                    mesg += "\"lambert\". It should be \"brdf\"."
                    RaiseError(self.Log_file, mesg)
            elif self.Surface_parameters[0] == "brdf":
                if self.Surface_parameters[1] != "ocean": #MODIF Mathieu BRDF
                    input_file = "./lib/surface/surfbrdf_"
                    input_file += self.Surface_parameters[1]
                    input_file += ".dat"
                    if not(os.path.lexists(input_file)):
                        mesg = "surface.py : The specified surface type is"
                        mesg += "\"" + self.Surface_parameters[1] + "\""
                        mesg += ", but the file \"" + input_file
                        mesg += "\" is not a valid file."
                        RaiseError(self.Log_file, mesg) #MODIF MATHIEU BRDF
            else:
                mesg = "surface.py : The surface family (\""
                mesg += self.Surface_parameters[0] + "\") is not recognized."
                mesg += " It should be either \"brdf\" or \"lambert\"."
                RaiseError(self.Log_file, mesg)

            ## geometrics parameters.
            input_file = "./lib/solrad/solrad_" + self.Incident_spectrum +\
                ".dat"
            if not(os.path.lexists(input_file)):
                mesg = "geometrics.py : The specified incident spectrum is \""
                mesg += self.Incident_spectrum + "\" but the file \""
                mesg += input_file + "\" is not a valid file."
                RaiseError(self.Log_file, mesg)

            if self.Solar_configuration_mode == "angle":
                for item in self.Solar_configuration_list:
                    if not(isinstance(item[0], float)) \
                            and not(isinstance(item[0], int)):
                        mesg = "geometrics.py : One item in \"Solar_zenith_an"
                        mesg += "gle_list\" is not a single float (\""
                        mesg += str(item[0]) + "\")."
                        RaiseError(self.Log_file, mesg)
            elif self.Solar_configuration_mode == "position":
                for item in self.Solar_configuration_list:
                    if len(item) != 4:
                        mesg = "geometrics.py : One geographic position (\""
                        mesg += str(item) + "\") is not valid. The template"
                        mesg += " is [longitude, latitude, day of the year,"
                        mesg += " time (decimal, UT)]."
                        RaiseError(self.Log_file, mesg)
                    else:
                        for value in item:
                            if not(isinstance(value, float)) \
                                    and not(isinstance(value, int)):
                                mesg = "geometrics.py : One geographic positi"
                                mesg += "on (\"" + str(item) + "\") is not valid."
                                mesg += "The template"
                                mesg += " is [longitude, latitude, day of the year,"
                                mesg += " time (decimal, UT)], each item being a "
                                mesg += "float or integer."
                                RaiseError(self.Log_file, mesg)                                
            else:
                mesg = "geometrics.py : The specified solar configuration mode"
                mesg += " (\"" + self.Solar_configuration_mode + "\") is not "
                mesg += "not recognized."
                RaiseError(self.Log_file, mesg)
