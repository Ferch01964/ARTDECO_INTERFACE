##############################################################################
##                                                                          ##
## No user changes should be necessary under this line.                     ##
##                                                                          ##
##############################################################################
import os, os.path, sys, glob, shutil
sys.path.insert(0, "./src/python")
sys.path.insert(0, "./input")

if len(sys.argv) != 3:
    mesg = "\nUsage :\n      - python Artdeco.py usage : displays"
    mesg += " this message.\n      - python Artdeco.py"
    mesg += " new <configuration_directory> : Creates a new "
    mesg += "configuration directory <configuration_directory> and copies "
    mesg += "reference configuration files in it.\n      - python Artdeco.py"
    mesg += " help-configuration <configuration_directory> : "
    mesg += "Indicates the additional configuration " 
    mesg += "files (depending on the configuration in ./config/<configuration"
    mesg += "_directory>.main.py) that must be checked before running the"
    mesg += " code.\n      - python Artdeco.py check-configuration "
    mesg += "<configuration_directory> : Checks for any configuratio"
    mesg += "n errors.\n      - python Artdeco.py configure "
    mesg += "<configuration_directory> : Checks for any "
    mesg += "configuration errors and edits the text "
    mesg += "input files.\n      - python Artdeco.py execute "
    mesg += "<configuration_directory> : Check"
    mesg += "s for any configuration errors, edits"
    mesg += " the text input files and launches the program.\n"
    print(mesg)

else:
    config_dir = sys.argv[2] + "/"
    if sys.argv[1] == "new":
        if os.path.lexists(config_dir):
            mesg = "The directory " + config_dir + \
                " cannot be created because it already exists."
            print(mesg)
        else:
            os.mkdir(config_dir)
            for f in glob.glob("./config/reference/*.py"):
                name = f.split("/")[-1]
                shutil.copyfile(f, config_dir + name)

    else:
        if not(os.path.isdir(config_dir)):
            mesg = "\nThe directory " + config_dir + " is not a valid one.\n"
            mesg += "\nUsage :\n      - python Artdeco.py usage : displays"
            mesg += " this message.\n      - python Artdeco.py"
            mesg += " new <configuration_directory> : Creates a new "
            mesg += "configuration directory <configuration_directory> and copies "
            mesg += "reference configuration files in it.\n      - python Artdeco.py"
            mesg += " help-configuration <configuration_directory> : "
            mesg += "Indicates the additional configuration " 
            mesg += "files (depending on the configuration in ./config/<configuration"
            mesg += "_directory>.main.py) that must be checked before running the"
            mesg += " code.\n      - python Artdeco.py check-configuration "
            mesg += "<configuration_directory> : Checks for any configuratio"
            mesg += "n errors.\n      - python Artdeco.py configure "
            mesg += "<configuration_directory> : Checks for any "
            mesg += "configuration errors and edits the text "
            mesg += "input files.\n      - python Artdeco.py execute "
            mesg += "<configuration_directory> : Check"
            mesg += "s for any configuration errors, edits"
            mesg += " the text input files and launches the program.\n"
            print(mesg)

        elif not(os.path.isfile(config_dir + "main.py")):
            mesg = "\nThe file " + config_dir + "main.py is not a valid file.\n"
            mesg += "\nUsage :\n      - python Artdeco.py usage : displays"
            mesg += " this message.\n     - python Artdeco.py"
            mesg += " new <configuration_directory> : Creates a new "
            mesg += "configuration directory <configuration_directory> and copies "
            mesg += "reference configuration files in it.\n       - python Artdeco.py"
            mesg += " help-configuration <configuration_directory> : "
            mesg += "Indicates the additional configuration " 
            mesg += "files (depending on the configuration in ./config/<configuration"
            mesg += "_directory>.main.py) that must be checked before running the"
            mesg += " code.\n      - python Artdeco.py check-configuration "
            mesg += "<configuration_directory> : Checks for any configuratio"
            mesg += "n errors.\n      - python Artdeco.py configure "
            mesg += "<configuration_directory> : Checks for any "
            mesg += "configuration errors and edits the text "
            mesg += "input files.\n      - python Artdeco.py execute "
            mesg += "<configuration_directory> : Check"
            mesg += "s for any configuration errors, edits"
            mesg += " the text input files and launches the program.\n"
            print(mesg)
        
        else:
            for f in glob.glob("./input/*.py"):
                os.remove(f)
            shutil.copyfile(config_dir + "main.py", "./input/main.py")
            import main
            Mode, With_particles, With_atmosphere, RTE_solver, \
                Legendre_and_truncation, Rayleigh_scattering, \
                Stokes_components, Prefix_of_input_files, Verbose, Warnings, \
                Log_file, Output_directory, Output_format, Print_betal, \
                Print_recomposed_phase_matrix = main.BasicParameters()

            if sys.argv[1] == "help-configuration":
                mesg = "\nBefore running the program, you should check and modify some "
                mesg += "additional configuration files in " + config_dir + " :\n"
                mesg += "    - "
                mesg += Mode + ".py\n"
                if RTE_solver != "none":
                    if RTE_solver != "sinsca":
                        mesg += "    - " + RTE_solver + ".py\n"
                    mesg += "    - surface.py\n    - geometrics.py\n"
                if With_particles:
                    mesg += "    - particles.py\n"
                if With_atmosphere:
                    mesg += "    - atmosphere.py\n"
                print(mesg)       
    
            elif sys.argv[1] == "check-configuration" \
                    or sys.argv[1] == "execute" \
                    or sys.argv[1] == "configure":

                ## Copies needed configuration files in ./input.
                config_list = [Mode + ".py"]
                if RTE_solver != "none":
                    config_list.append("surface.py")
                    config_list.append("geometrics.py")
                    if RTE_solver != "sinsca":
                        config_list.append(RTE_solver + ".py")
                if With_atmosphere:
                    config_list.append("atmosphere.py")
                if With_particles:
                    config_list.append("particles.py")
                is_valid_file = True
                for item in config_list:
                    if not(os.path.isfile(config_dir + item)):
                        mesg = "\nThe file " + config_dir + item + " is not a valid file.\n"
                        is_valid_file = False
                        print(mesg)
                    else:
                        shutil.copyfile(config_dir + item, "./input/" + item)

                if is_valid_file:
                    ## Initializes this configuration.
                    from ArtdecoClass import *
                    ThisConfig = ArtdecoConfig()

                    ThisConfig.Verbose = Verbose
                    ThisConfig.Warnings = Warnings

                    ## Updates the directories and file names.
                    ThisConfig.Input_directory = os.getcwd() + "/input/"
                    ThisConfig.Prefix_of_input_files = Prefix_of_input_files
                    ThisConfig.Output_subdirectory = Output_directory
                    ThisConfig.Log_file = Log_file

                    ## Output options.
                    ThisConfig.Output_format = Output_format
                    ThisConfig.Print_recomposed_phase_matrix = Print_recomposed_phase_matrix
                    ThisConfig.Print_betal = Print_betal

                    ## RTE_solver specific files.
                    ThisConfig.RTE_solver = RTE_solver
                    ThisConfig.RadiativeModelConfig()

                    ## Polarization model (1, 3 or 4).
                    ThisConfig.Stokes_components = Stokes_components

                    ## Mode.
                    ThisConfig.Mode = Mode
                    ThisConfig.ModeConfig()

                    ## Atmosphere.
                    ThisConfig.With_atmosphere = With_atmosphere
                    if ThisConfig.With_atmosphere:
                        ThisConfig.AtmosphereConfig()

                        ## Particles and phase matrix truncation.
                    ThisConfig.Legendre_and_truncation = Legendre_and_truncation
                    ThisConfig.With_particles = With_particles
                    if ThisConfig.With_particles:
                        ThisConfig.ParticleConfig()

                    ## If a RTE solver is used, reads the configuration for surface 
                    ## and geometrics.
                    if ThisConfig.RTE_solver != "none":
                        ThisConfig.SurfaceConfig()
                        ThisConfig.GeometricsConfig()

                    ## Updates the keywords with the parameters given by the user.
                    ThisConfig.Rayleigh_scattering = Rayleigh_scattering
                    ThisConfig.UpdateKeywords()

                    ## Checks for configuration errors.
                    ThisConfig.ChecksConfigurationErrors()

                    if sys.argv[1] == "configure" \
                            or sys.argv[1] == "execute":

                        ## Writes the main configuration file.
                        ThisConfig.WriteMainConfigurationFile()

                        ## Writes specific configuration files for RTE solvers.
                        ThisConfig.WriteRadiativeModelConfigFile()
                        ThisConfig.SymlinkToRadiativeModelConfigFile()

                        ## If necessary, writes specific configuration files for
                        ## the phase matrix truncation method.
                        ThisConfig.WriteTruncationMethodConfigFile()
                        ThisConfig.SymlinkToTruncationMethodConfigFile()

                        ## If necessary, writes definition file for the 
                        ## Kdis parameterization.
                        ThisConfig.WriteKdisFiles()

                        ## If necessary, writes definition file for the
                        ## ocean parameterization.
                        ThisConfig.WriteOceanConfigFile()

                        ## If necessary, writes definition files for the
                        ## atmosphere profile.
                        ThisConfig.WriteAtmosphericProfile()
                        ThisConfig.WriteUniformGasesConcentrationFile()

                        ## If necessary, copies the reference input files.
                        ThisConfig.WriteReferenceInputFiles()

                        ## If necessary, copies the vertical distribution files.
                        ThisConfig.WriteVdistFiles()
    
                        if sys.argv[1] == "execute":
                            import subprocess
                            os.environ["PATH_ARTDECO"] = os.getcwd() + "/"
                            log_output = open(ThisConfig.Log_file, "w")

                            ##Cleans the output directory.
                            if os.path.isdir("./out/" + ThisConfig.Output_subdirectory):
                                shutil.rmtree("./out/" + ThisConfig.Output_subdirectory)
                            os.mkdir("./out/" + ThisConfig.Output_subdirectory)

                            ## Checks if artdeco is already compiled.
                            if not(os.path.lexists("./src/artdeco")):
                                mesg = "Artdeco.py : the program artdeco is not compiled in ./src directory."
                                print(mesg)
                                
                            else:
                                ## Calls the main program ARTDECO.
                                subprocess.call(["src/artdeco", "input/" + \
                                                     ThisConfig.Prefix_of_input_files + \
                                                     "_artdeco_in.dat", \
                                                     ThisConfig.Output_subdirectory], \
                                                    stdout = log_output, stderr = log_output)
                                
                                ## Writes output file in the specified format.
                                ThisConfig.WriteOutput(ThisConfig.Output_format)

                                ## Tags the outpur files name with ARTDECO_*
                                for f in glob.glob("./out/" + ThisConfig.Output_subdirectory + "/*.dat"):
                                    f_name = f.split("/")[-1]
                                    os.rename(f, "./out/" + ThisConfig.Output_subdirectory + "/ARTDECO_" + f_name)
                                for f in glob.glob("./out/" + ThisConfig.Output_subdirectory + "/*.nc"):
                                    f_name = f.split("/")[-1]
                                    os.rename(f, "./out/" + ThisConfig.Output_subdirectory + "/ARTDECO_" + f_name)
                                for f in glob.glob("./out/" + ThisConfig.Output_subdirectory + "/*.hdf5"):
                                    f_name = f.split("/")[-1]
                                    os.rename(f, "./out/" + ThisConfig.Output_subdirectory + "/ARTDECO_" + f_name)
                                
                                ## Copies the log file in the output directory.
                                output_dir = "./out/" + ThisConfig.Output_subdirectory
                                output_dir += "/log/"
                                os.mkdir(output_dir)
                                shutil.copyfile(Log_file, output_dir + Log_file)

                                ## Copies configuration files in the output directory.
                                output_dir = "./out/" + ThisConfig.Output_subdirectory
                                output_dir += "/config/"
                                os.mkdir(output_dir)
                                shutil.copyfile(config_dir + "main.py", output_dir + "main.py")
                                
                                filename = ThisConfig.Mode + ".py"
                                shutil.copyfile(config_dir + filename, output_dir + filename)
                                if ThisConfig.RTE_solver != "none":
                                    if ThisConfig.RTE_solver != "sinsca":
                                        filename = ThisConfig.RTE_solver + ".py"
                                        shutil.copyfile(config_dir + filename, output_dir + filename)
                                    shutil.copyfile(config_dir + "surface.py", output_dir + "surface.py")
                                    shutil.copyfile(config_dir + "geometrics.py", output_dir + "geometrics.py")
                                if ThisConfig.With_particles:
                                    shutil.copyfile(config_dir + "particles.py", output_dir + "particles.py")
                                if ThisConfig.With_atmosphere:
                                    shutil.copyfile(config_dir + "atmosphere.py", output_dir + "atmosphere.py")

            else:
                mesg = "\nUsage :\n      - python Artdeco.py usage : displays"
                mesg += " this message.\n      - python Artdeco.py"
                mesg += " new <configuration_directory> : Creates a new "
                mesg += "configuration directory <configuration_directory> and copies "
                mesg += "reference configuration files in it.\n      - python Artdeco.py"
                mesg += " help-configuration : Indicates the additional configuration " 
                mesg += "files (depending on the configuration in Artdeco.py"
                mesg += ") that must be checked before running the code.\n      - python "
                mesg += "Artdeco.py check-configuration : Checks for any configuratio"
                mesg += "n errors.\n      - python Artdeco.py configure : Checks for any "
                mesg += "configuration errors and edits the text "
                mesg += "input files.\n      - python Artdeco.py execute : Check"
                mesg += "s for any configuration errors, edits"
                mesg += " the text input files and launches the program.\n"
                print(mesg)
