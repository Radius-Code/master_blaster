import os
import sys
from functools import partial
from master_blaster import program

class program2():

    def __init__(self):

        self.program = program()
        self.pending_functions = ["launch_blast", "parsing", "file_list",
                                     "next_selected_file", "format_file", "store_queries_in_file",
                                     "homolog_filter", "comparation", "add_sequence_to_list"]
        self.functions_with_args = {
            "format_file": self.program.format_file,
            "store_queries_in_file": self.program.store_queries_in_file,
            "homolog_filter": self.program.homolog_filter,
        }

        self.program.blast = True

    def restart_elements(self):

        if self.pending_functions == self.program.run_functions:
            self.program.run_functions = []
            self.program.last_run_function = None



    def run_program(self):
      
        try:
            self.restart_elements()

            for function in self.pending_functions:
                if function not in self.program.run_functions:
                    if function in self.functions_with_args:
                        argument = self.functions_with_args[function]
                        partial_func = partial(getattr(self.program, function), argument)
                        partial_func()
                        self.program.run_functions.append(function)
                        self.program.last_run_function = function
                    else:
                        full_func = getattr(self.program, function)
                        full_func()
                        self.program.run_functions.append(function)
                        self.program.last_run_function = function
                else:
                    pass

           

        except KeyboardInterrupt:
            self.program.save_data()
            print("Data Saved")
            sys.exit(0)

        except Exception as e:
            print(e)

    def ejecution(self):


        while self.program.blast:

            self.program.load_pickle_data()

            if len(self.program.sequences) != self.program.counter:

                self.program.inputs()

                self.program.end(self.program.blast)

                if not self.program.blast:
                    print("Program executed correctly")
                    os.remove("counter.pickle")
                    break


                # Ejecutar el program
                self.run_program()
                

                self.program.counters()

                self.program.save_pickle_data()

            else:

                self.program.inputs()

                self.delete_final_files()

               
    def delete_previous_files(self):

        patrons = ["qblast_blastn_", "output_", "counter"]

        found_files = False

        for file_name in os.listdir(os.getcwd()):
            for patron in patrons:
                if file_name.startswith(patron):
                    found_files = True
                    break

        if found_files:
            delete_files = input("¿Do you want to delete previous blast files? (Yes/No): ")
            if delete_files == "Yes":
                for file_name in os.listdir(os.getcwd()):
                    for patron in patrons:
                        if file_name.startswith(patron):
                            os.remove(file_name)
                            print(f"The file {file_name} has been deleted.")
                           
        else:
            pass



    def delete_final_files(self):
            

        self.xml_files = [f"qblast_blastn_{n}.xml" for n in range(len(self.program.sequences))]  
        self.csv_files = [f"output_{o}.csv" for o in range(len(self.program.sequences))]
        
        for archivo1, archivo2 in zip(self.xml_files, self.csv_files):

            if os.path.isfile(archivo1) or os.path.isfile(archivo2) or os.path.isfile("counter.pickle"):
                self.delete()
                
            elif not (os.path.isfile(archivo1) or os.path.isfile(archivo2)) and os.path.isfile("counter.pickle"):
                continue
                


    def delete(self):

        
        # Generate lists of filenames

        delete_files = input("Do you want to delete current xml, csv and pickle files?: Yes / No: ")

        if delete_files == "Yes":
            # Remove "counter.pickle" file
            try:
                os.remove("counter.pickle")
                print("counter.pickle removed successfully.")
            except Exception as e:
                print(f"Error while removing counter.pickle: {e}")

            # Remove CSV files
            for csv_file in self.csv_files:
                try:
                    os.remove(csv_file)
                    print(f"File '{csv_file}' removed successfully.")
                except Exception as e:
                    print(f"Error while removing file '{csv_file}': {e}")

            # Remove XML files
            for xml_file in self.xml_files:
                try:
                    os.remove(xml_file)
                    print(f"File '{xml_file}' removed successfully.")
                except Exception as e:
                    print(f"Error while removing file '{xml_file}': {e}")

        else:
            sys.exit(0)



        relaunch = input("Do you want to relaunch the script?: Yes / No: ")
        
        if relaunch == "Yes":
            self.program.blast = True
            self.ejecution()

        else:
            sys.exit(0)

    

program2 = program2()


program2.delete_previous_files()

program2.ejecution()


#mover {secuencia proteina}_final_report.txt a carpeta saves homologues
