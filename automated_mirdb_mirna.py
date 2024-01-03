from selenium import webdriver
from selenium.webdriver.support.ui import Select
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from bs4 import BeautifulSoup
from collections import OrderedDict
import pandas as pd
from Bio import SeqIO
import argparse
import logging
import sqlite3


class MirdbSearch():
    
    # Constructor:
    def __init__(self, fasta, species, cutoff):
        
        self.url = 'http://www.mirdb.org/custom.html'
        self.fasta = self._open_fasta(fasta)
        self.species = species
        self.submission = 'miRNA Sequence'
        self.cutoff = cutoff
    
    # Methods:
    # _open_fasta() method: the file is automatically closed after operations are
    # performed. Use '-> list' to specify that we are interested in returning
    # a list. 
    def _open_fasta(self, fasta) -> list:
        # Use try and except to handle the IOError exception that occurs when the
        # open() function does not open the file. 
        try:
            with open(fasta) as handle:
                # SeqIO.parse(filename, format) returns a SeqRecord object:
                return list(SeqIO.parse(handle, 'fasta'))
        except IOError:
            logging.exception("Could not read file.")
            raise IOError("Could not read file.")


class Crawler():
    
    # Constructor:
    def __init__(self, visible):

        # 'driver' (selenium package) allows us to interact with the browser and websites.
        self.driver = self._create_driver(visible)

    # Methods:
    # _create_driver() method: initialises the driver. The FireFox browser will open
    # and if the visible variable is True, the window will be shown, while if False,
    # the browser will open in headless mode. By default, visible is False:
    def _create_driver(self, visible):
        options = Options()
        if args.visible:
            return webdriver.Firefox(options=options)
        else:
            options.add_argument('--headless')
            return webdriver.Firefox(options=options)
    
    # select_element() method: uses different functions to locate and select an element
    # and finally select the value that we want in the element:
    def select_element(self, name, value):
        select = Select(self.driver.find_element('name', name))
        select.select_by_visible_text(value)
    
    # enter_sequence() method: locates an element and use send_keys() to input text
    # into that element:
    def enter_sequence(self, sequence):
        self.driver.find_element('name','customSub').send_keys(sequence)
    
    # continue_to_results() method: locate an element and click on it.  
    def continue_to_results(self):
        try:
            self.driver.find_element('xpath', '/html/body/table[2]/tbody/tr/td[3]/form/table/tbody/tr[5]/td/input[1]').click()
        except:
            pass

    # wait_results() method: includes a timeout to show an element that allows to 
    # access the results. If this element appears, click on it, if not, the TimeoutError
    # exception is thrown.
    def wait_results(self):
        timeout = 30
        result = EC.presence_of_element_located((By.XPATH, '/html/body/form/input[2]'))
        try:
            WebDriverWait(self.driver, timeout).until(result)
            self.driver.find_element('xpath', '/html/body/form/input[2]').click()
        except TimeoutError:
            logging.exception("TimeoutError")

class Scraper():
    
    # Constructor:
    def __init__(self):
        
        self.soup = None
    
    # Methods:
    # parse() method: formats the HTML code.
    def parse(self, page):
        self.soup = BeautifulSoup(page, 'html.parser')
    
    # get_above_cutoff() method: saves row indexes of table results that have a
    # target score value greater than or equal to the cutoff.
    def get_above_cutoff(self, cutoff):
        passed_cutoff = []
        rows = self.soup.find('table', id='table1').find('tbody').find_all('tr')
        for i in range(1, len(rows)):
            try:
                cells = rows[i].find_all('td')
                if int(cells[2].text) >= cutoff:
                    passed_cutoff.append(i)
            except AttributeError:
                logging.exception("AttributeError")
                pass
        return passed_cutoff
    
    # get_score() method: keeps the score of the results.
    def get_score(self):
        # usually the score is in cell #7, but sometimes there is an extra row for miRNA previous name
        table = self.soup.find_all('td')
        try:
            if table[7].text.isdigit():  # so check if cell #7 text is a digit
                score = int(table[7].text)
            else:                       # if it isn't then score should be in cell #9
                score = int(table[9].text)
        except AttributeError:
            logging.exception("AttributeError")
        return score
    
    # get_number_of_seeds() method: returns the number of seeds (find seeds by text colour).
    def get_number_of_seeds(self):
        seeds = self.soup.find_all('font', {'color': '#0000FF'}) 
        try:
            number_of_seeds = len(seeds)
        except AttributeError:
            number_of_seeds = None
        return number_of_seeds

    # get_gene_symbol() method: keeps Gene Symbol of the results (cell #19):
    def get_gene_symbol(self):
        table = self.soup.find_all('td')
        try:
            table[19].text.isupper()  # Check that all letters are capitalized
            gene_symbol = table[19].text
        except AttributeError:
            logging.exception("AttributeError")
        return gene_symbol
    
    # get_genbank_accession() method: searches for <a> (HTML link tag) and saves its
    # text (the link tag corresponding to GenBank Accession is the second one).
    def get_genbank_accession(self):
        links = scraper.soup.find_all('a', href=True)
        try:
            genbank_accession = links[2].font.text
        except AttributeError:
            genbank_accession = None
        return genbank_accession

class Target():
    
    # Constructor:
    def __init__(self):
        
        self.sequence = None
        self.mirna_name = None
        self.score = None
        self.number_of_seeds = None
        self.gene_symbol = None
        self.genbank_accession = None

class Database():
    
    # Constructor:
    def __init__(self, database: str):
        
        # __connection() method: SQLite allows us to create a database for internal 
        # data storage. Use sqlite3.connect() to generate a database and establish 
        # a connection for accessing it. 
        self.__connection = sqlite3.connect(database)
        self.__create_tables()

    # Connection property store the connection to the database (returns a sqlite3.Connection).
    # In the constructor we are initializing __connection and here we return the value 
    # of __connection:
    @property
    def connection(self) -> sqlite3.Connection:
        return self.__connection
    
    # Methods:
    # Private __execute_query() method: executes the query (returns a list):
    def __execute_query(self, *args) -> list:
        try:
            # Initialize the database cursor (connection between the database and the query):
            cursor = self.connection.cursor()
            # Execute the query provided in the *args variable.
            cursor.execute(*args)
            self.connection.commit()
            # Extract the query results:
            result = cursor.fetchall()
        # Finally is always executed after a try. Close the cursor:
        finally:
            cursor.close()
        return result

    # Private __create_tables() method: creates the database table.
    def __create_tables(self):
        self.__execute_query("""
        CREATE TABLE IF NOT EXISTS prueba_final_2 (
        sequence TEXT,
        mirna TEXT,
        score INTEGER,
        seeds INTEGER,
        GeneSymbol TEXT,
        GenBankAccession TEXT 
    );""")

    # insert_target() method: inserts the data into the table.
    def insert_target(self, target: Target):
        self.__execute_query("""INSERT INTO prueba_final_2 (sequence, mirna, score, seeds, GeneSymbol, GenBankAccession)
                                    VALUES (?,?,?,?,?,?)""",
                            [target.sequence,target.mirna_name,target.score,target.number_of_seeds,target.gene_symbol,target.genbank_accession])
    
    # export_to_csv() method: exports all data in a .csv file.
    def export_to_csv(self, output_name: str):
        table = pd.read_sql_query("SELECT * FROM prueba_final_2", self.connection)
        table.to_csv(output_name, index=None, header=True)


if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='Automated miRDB custom microRNA target prediction search: '
                                                '\nThis script uses a webdriver to access the miRDB website and search'
                                                ' for mRNA targets prediction with user-provided miRNA.')
    
    # Positional arguments:
    parser.add_argument('inp', type=str, help='Input FASTA file with miRNA sequences')
    parser.add_argument('out', type=str, help='Name for output file')
    parser.add_argument('sp', type=str, choices=["Human","Rat","Mouse","Chicken", "Dog"], help='Species')
    
    # Optional arguments (default values can be changed here):
    parser.add_argument('-c', '--cutoff', type=int, help='Score cut-off <int> (default: 95)', default=95)
    parser.add_argument('-v', '--visible', action='store_true', help='Shows browser window during the process'
                                                                      '(default: False)', default=True)

    args = parser.parse_args()

    logging.basicConfig(filename="mirdb_error.log",format='%(asctime)s - %(message)s', level=logging.INFO)

    search = MirdbSearch(args.inp, args.sp, args.cutoff)
    crawler = Crawler(args.visible)
    scraper = Scraper()
    database = Database(f'{args.out}.db')

    try:
        for n, sequence in enumerate(search.fasta):
            if 17 <= len(sequence.seq) <= 30:  # Length restriction in miRDB
                print(f'\rSearching targets for sequence: {n+1}/{len(search.fasta)} - {sequence.id}', end='', flush=True)
                crawler.driver.get(search.url) # Load the miRDB website
                crawler.select_element('searchSpecies', search.species)
                crawler.select_element('subChoice', search.submission)
                crawler.enter_sequence(sequence.seq)
                crawler.continue_to_results()
                crawler.wait_results()
                html = crawler.driver.page_source # Extract the HTML code
                scraper.parse(html)
                passed_cutoff_rows = scraper.get_above_cutoff(search.cutoff)

                for row in passed_cutoff_rows:
                    details = crawler.driver.find_elements('name', '.submit')  # The first is the "Return" button, the others are "Target Details"
                    try:
                        details[row].click()
                    except IndexError:
                        logging.exception(f"{sequence.id}")
                        continue
                    html = crawler.driver.page_source
                    scraper.parse(html)
                    target = Target()
                    target.sequence = sequence.seq._data
                    target.mirna_name = sequence.id
                    target.score = scraper.get_score()
                    target.number_of_seeds = scraper.get_number_of_seeds()
                    target.gene_symbol = scraper.get_gene_symbol()
                    target.genbank_accession = scraper.get_genbank_accession()
                    database.insert_target(target)
                    crawler.driver.back() # Go to previous page in the browser
            else:
                print(f'\nFailed to search {sequence.id}. Sequence length out of range ({len(sequence)} nt).')
                logging.error(f"{sequence.id} length ({len(sequence)} nt) out of range (17 - 30 nt)")
    except:
        logging.exception("Error")
        print(f'\nAn exception has occurred while searching: {sequence.id}')
    finally:
        database.export_to_csv(f'{args.out}.csv')
        print(f'\nResults saved to {args.out}.csv')
        crawler.driver.close() # Close the browser
        database.connection.close # Close the connection to the database
