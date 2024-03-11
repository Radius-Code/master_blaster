import pandas as pd
import numpy as np
from io import open
import os
import subprocess
from Bio.Blast import NCBIWWW
import xml.etree.ElementTree as ET
import csv
import time
import pickle
import sys


class program():

    def __init__ (self):


        self.n = 0
        self.o = 0
        self.counter = 0
        self.th = 0

        self.data = []

        self.actual_directory = os.getcwd()
        self.files_xml = [file for file in os.listdir(self.actual_directory) if file.endswith('.xml')]  # Obtener una lista de files XML en el directorio actual
        self.first_time_blast = True
        self.blast = True

        self.next_file = None
        
        self.entrez_query = ""
        self.sequence = ""

        self.csvs = []

        
        self.sequences = []

        self.new_line = ""

        self.run_functions = []
        self.last_run_functions = None 
        

        
    def counters(self):

        if not self.first_time_blast:

            self.counter+=1

            self.o+=1

            self.n+=1



    def load_pickle_data(self):

        try:

            if os.path.isfile('counter.pickle'):
                with open('counter.pickle', 'rb') as file:
        
                    data = pickle.load(file)


                
            
                    self.counter, self.o, self.n, self.sequences, self.blast, self.first_time_blast, self.th, self.sequence, self.entrez_query, self.run_functions, self.last_run_functions = data
                    #print(self.counter, self.o, self.n, self.sequences, self.blast, self.first_time_blast, self.th, self.sequence, self.entrez_query, self.run_functions, self.last_run_functions)
                    
                file.close()

            else:
                pass
            
        except Exception as e:
            print(e)






    def inputs(self):



        try:

            if self.first_time_blast:
                
                self.sequence = input("Introduce a sequence to run: ").replace(".1", "")
                self.entrez_query = input("Introduce the code / name of Organism to scan (ej: taxid:'number' e.g.: 4686 OR 'name [Organism]' e.g.: Asparagus officinalis [Organism]): ")

                self.sequences.append(self.sequence)
                self.first_time_blast = False
                
        

            else:  


                if len(self.sequences) != 1:
                    for i in range(len(self.sequences)):
                        next_sequence = self.sequences[self.counter] # Append the next sequence to the sequences list
                        self.sequence = next_sequence
                        break

                print("Organism: ", self.entrez_query)
                print("Sequence to blast: ", self.sequence)
                print("\n")
                
            if self.th == 0:

                self.th = input("Introduce a threshold of Identity and Coverage Values (default E Value is 10^-20): ")
                self.th = int(self.th)

            else:

                self.th = int(self.th)
       

        except Exception as e:

            print("No more sequences to analyze.", "\n", "End of script")
            




    def launch_blast(self):


        try:

            
            print("\n")
            print("Blasting file", f"qblast_blastn_{self.o}.xml", "...")

            result_handle = NCBIWWW.qblast("blastp", "refseq_protein", self.sequence, entrez_query=self.entrez_query, format_type="XML")
            result_handle

            print("Blast ended, saving file...")

            save_file = open(f"qblast_blastn_{self.o}.xml", "w")

            save_file.write(result_handle.read())

            save_file.close()

            result_handle.close()


            print("File saved")

             

        except Exception as e:
                    
            print("End of script, " + str(e))
            sys.exit(0)
            
            

        


    def parsing(self):


        print("\n")
        print("Parsing XML", f"qblast_blastn_{self.n}.xml", "...")


        # Line with the parse command
        tree = ET.parse(f"qblast_blastn_{self.n}.xml")
        root = tree.getroot()

        file = f'output_{self.n}.csv'
        # Open a CSV file to write the extracted data
        with open(file, 'w', newline='') as csvfile:
            fieldnames = ['query', 'subject', 'identity', 'align_length', 'mismatches', 'gap_opens', 'q.start', 'q.end', 's.start', 's.end', 'e_value', 'bit_score', '%positives']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            # Find Hit_accession elements and extract information
            for hit in root.findall('.//Hit'):
                hit_accession = hit.find('Hit_accession').text
                for hsp in hit.findall('.//Hit_hsps/Hsp'):
                    datas = {
                        'query': hit.find('.//Hit_def').text,
                        'subject': hit_accession,
                        'identity': hsp.find('Hsp_identity').text,
                        'align_length': hsp.find('Hsp_align-len').text,
                        'mismatches': hsp.find('Hsp_mismatch').text if hsp.find('Hsp_mismatch') is not None else None,
                        'gap_opens': hsp.find('Hsp_gaps').text,
                        'q.start': hsp.find('Hsp_query-from').text,
                        'q.end': hsp.find('Hsp_query-to').text,
                        's.start': hsp.find('Hsp_hit-from').text,
                        's.end': hsp.find('Hsp_hit-to').text,
                        'e_value': hsp.find('Hsp_evalue').text,
                        'bit_score': hsp.find('Hsp_bit-score').text,
                        '%positives': hsp.find('Hsp_positive').text
                    }
                    writer.writerow(datas)
                
            print("Parsed completed")
           
    


    def file_list(self):

        # Get a list of CSV files in the current directory
        self.csvs = [x for x in os.listdir(self.actual_directory) if os.path.isfile(os.path.join(self.actual_directory, x)) and x.endswith(".csv")]

        # Sort list of CSV files by modification date
        self.csvs.sort(key=lambda x: os.stat(os.path.join(self.actual_directory, x)).st_mtime)


       

    def next_selected_file(self):

        #print(self.csvs)
        if len(self.csvs) == 1:
            for file in self.csvs:
                self.next_file = file
                break

        
        else:
            for file in self.csvs:
                ruta_file = os.path.join(self.actual_directory, file)
                if self.next_file is None or os.stat(ruta_file).st_mtime > os.stat(os.path.join(self.actual_directory, self.next_file)).st_mtime:
                    self.next_file = file




    def format_file(self, data):


        try:

            df = pd.read_csv(self.next_file, index_col=False)   # Open file with no header


            if not df.empty and "align_length" in df.columns:
                first_align_length = df["align_length"].iloc[0]
                
                if first_align_length != 0:
                    df['coverage'] = df['align_length'].div(first_align_length) * 100 # Get the coverage of columns as a percentage

                else:
                    print("The first value of 'align_length' is zero, cannot divide by zero.")
            else:
                print("DataFrame is empty or 'align_length' column doesn't exist.")

            coverage = df['coverage']                           
            df = df.drop(columns=['coverage'])                
            df.insert(loc=3, column='coverage', value=coverage) # Add coverage column with previous calculated values


            df['coverage'] = df['coverage'].round(2) # Round Identity and Coverage columns to two decimal places
            df['identity'] = df['identity'].round(2)

            df.to_csv(self.next_file, index=False) 

            return df

        except:
            print("There is no more homologous sequences")


    def store_queries_in_file(self, data): #Aplicamos filtro que elimina la o las query de la tabla y las almacena en un documento final o lista (Se modifica la tabla).


        df = pd.read_csv(self.next_file)

        rows = df[df['coverage'] == 100.00 ] # Add rows with 100 value in a new value

        queries = pd.Series(dtype="string")

        queries = rows['subject'].astype(str).drop_duplicates() # Drop duplicates and add the result to a value

        #print(queries)

        # Add to the end of lines
        with open("final_report.txt", 'a') as file:
        # Verify if fil to write exists
            if os.path.exists("final_report.txt"):
                # If file exists, read it and add all content to a variable
                with open("final_report.txt", 'r') as file:
                    content = file.read()
                    
                # Divide file in queries and add queries to list
                existent_queries = set(content.strip().split("\n"))
                queries_to_write = [query for query in queries if query not in existent_queries] # Delete all empty queries

                # Append to file the queries
                with open("final_report.txt", 'a') as file:
                    # Escribir las querys que no están en el file
                    for query in queries_to_write:
                        file.write(query + "\n")
                    

        df.drop(rows.index, inplace = True) # delete the query/s from the table

        df.to_csv(self.next_file, index=False) # Save file

        return df



    def homolog_filter(self, data):


        global homologous_sequences



        df = pd.read_csv(self.next_file)

        df_mask = (df.identity >= self.th) & (df.coverage >= self.th) & (df.e_value <= 10**-20) #creo una "máscara" para los valores de identidad y covertura mayores o iguales a 50 y con un e value menor o igual a 10^-20

        homologous_sequences = df[df_mask] # Add a mask to the homomologous_sequences (get the ones above that threshold)
       

        with open("last_homologous.txt", 'w') as file: # Introduce them into a file
            for element in homologous_sequences['subject']:
                print(element, file = file)

        return homologous_sequences 



    def comparation(self):

        # Read all the lines from the first file (final report) and store them in a list
        with open("final_report.txt", "r") as file1:
            lines_file1 = file1.readlines()

        # Read all lines from the second file (last homologous)
        with open("last_homologous.txt", "r") as file2:
            lines_file2 = file2.readlines()

                
        # Add the last homologous lines to the final report file
        with open("final_report.txt", "a") as file1:
            for line in lines_file2:
                if line not in lines_file1:
                    self.new_line = file1.write(line)
        
        # Empty last_homologoues file to posterior use of it
        with open("last_homologous.txt", "w") as file2:
            file2.write("")


        

    def add_sequence_to_list(self):
        try:
            with open("final_report.txt", "r") as file:
                for line in file:
                    element = line.strip("\n")
                    if element not in self.sequences:
                        self.sequences.append(element)
                        
            print("\n")
            print("List of sequences to analyze: ", self.sequences)

        except Exception as e:
            print(e)

    
    def save_pickle_data(self):
        
        
        data_to_save = [self.counter, self.o, self.n, self.sequences, self.blast, self.first_time_blast, self.th, self.sequence, self.entrez_query, self.run_functions, self.last_run_functions]
       

        with open("counter.pickle", "wb") as pickle_file:
            pickle.dump(data_to_save, pickle_file)



    def end(self, blast):


        if len(self.sequences)-1 < self.counter:
            self.blast = False
        else:
            self.blast


    