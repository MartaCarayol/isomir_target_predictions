
# Utilizamos selenium para automatizar la busqueda en el navegador
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

# Para ejecutar el script en la consola de Spyder, usamos chdir() para
# cambiar el directorio de trabajo, estableciendo como directorio donde se
# encuentre el archivo FASTA. Para ello, tenemos que escribir en la consola:

# import os
# os.chdir('C:/Users/marta/Desktop/masterBioinformaticaBioestadistica/TFM')

# Clase MirdbSearch:
class MirdbSearch():
    
    # Constructor:
    def __init__(self, fasta, species, cutoff):
        
        # Atributos:
        self.url = 'http://www.mirdb.org/custom.html'
        self.fasta = self._open_fasta(fasta)
        self.species = species
        self.submission = 'miRNA Sequence'
        self.cutoff = cutoff
    
    # Metodos (funciones):
    # Metodo _open_fasta(): con 'with open('archivo') as nombre_variable', el 
    # archivo se cierra automaticamente tras realizar las operaciones con el.
    # Utilizamos '-> list' para especificar que estamos interesados en que el 
    # metodo devuelva una lista (type hint), pero puede no serlo y no 
    # obtendremos un error.
    def _open_fasta(self, fasta) -> list:
        # Utilizamos try and except para controlar la excepcion IOError que se
        # produce cuando la funcion open() no abre el archivo.
        try:
            with open(fasta) as handle:
                # SeqIO.parse(nombre_archivo, formato) devuelve un objeto 
                # SeqRecord (se utiliza para contener, entre otros, una secuencia
                # (objeto Seq), un identificador (id), un nombre (name), etc.). 
                return list(SeqIO.parse(handle, 'fasta'))
        except IOError:
            logging.exception("Could not read file.")
            raise IOError("Could not read file.")


# Clase Crawler:
class Crawler():
    
    # Constructor:
    def __init__(self, visible):

        # Atributos:
        # El 'driver' es el objeto fundamental en selenium que nos permitira
        # interactuar con el navegador y los sitios web. 
        self.driver = self._create_driver(visible)

    # Metodos (funciones):
    # Metodo _create_driver(): se utiliza para inicializar el 'driver' en su 
    # version mas simple. En concreto, se abrira el navegador FireFox y si la
    # variable 'visible' es True, veremos la ventana, mientras que en caso de 
    # ser False, se abrira el navegador en modo headless (se oculta la ventana). 
    # Por defecto, visible es False:
    def _create_driver(self, visible):
        options = Options()
        if args.visible:
            return webdriver.Firefox(options=options)
        else:
            options.add_argument('--headless')
            return webdriver.Firefox(options=options)
    
    # Metodo select_element(): primero utilizamos find_element() para localizar
    # un elemento, despues seleccionamos dicho elemento utilizando Select() y 
    # finalmente seleccionamos el valor que queremos que aparezca en el elemento
    # con select_by_visible_text():
    def select_element(self, name, value):
        select = Select(self.driver.find_element('name', name))
        select.select_by_visible_text(value)
    
    # Metodo enter_sequence(): localizamos el elemento y utilizamos send_keys()
    # para enviar el texto a dicho elemento.
    def enter_sequence(self, sequence):
        self.driver.find_element('name','customSub').send_keys(sequence)
    
    # Metodo continue_to_results(): localizamos el elemento y utilizamos click()
    # para hacer click en el. 
    def continue_to_results(self):
        try:
            self.driver.find_element('xpath', '/html/body/table[2]/tbody/tr/td[3]/form/table/tbody/tr[5]/td/input[1]').click()
        except:
            pass

    # Metodo wait_results(): hay un timeout para que aparezca el elemento que nos
    # permite acceder a los resultados. Si aparece dicho elemento, hacemos click,
    # sino, se lanza la excepcion TimeoutError.
    def wait_results(self):
        timeout = 30
        result = EC.presence_of_element_located((By.XPATH, '/html/body/form/input[2]'))
        try:
            WebDriverWait(self.driver, timeout).until(result)
            self.driver.find_element('xpath', '/html/body/form/input[2]').click()
        except TimeoutError:
            logging.exception("TimeoutError")

# Clase Scraper:
class Scraper():
    
    # Constructor:
    def __init__(self):
        
        # Atributos:
        self.soup = None
    
    # Metodo parse(): formatea el  codigo html.
    def parse(self, page):
        self.soup = BeautifulSoup(page, 'html.parser')
    
    # Metodo get_above_cutoff(): guarda los indices de las filas de la tabla 
    # resultado que tienen un valor de target score superior o igual al cutoff.
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
    
    # Metodo get_score(): obtenemos el score de los resultados.
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
    
    # Metodo get_number_of_seeds(): devuelve el numero de seeds buscandolas por
    # el color.
    def get_number_of_seeds(self):
        seeds = self.soup.find_all('font', {'color': '#0000FF'})  # find seeds by text color
        try:
            number_of_seeds = len(seeds)
        except AttributeError:
            number_of_seeds = None
        return number_of_seeds

    # Metodo get_gene_symbol(): obtenemos el Gene Symbol de los resultados.
    def get_gene_symbol(self):
        # Corresponde con la posicion 19 (viendo el codigo fuente de la pagina)
        table = self.soup.find_all('td')
        try:
            table[19].text.isupper()  # Comprobamos que todas las letras son mayusculas
            gene_symbol = table[19].text
        except AttributeError:
            logging.exception("AttributeError")
        return gene_symbol
    
    # Metodo get_genbank_accession(): busca las <a> (etiqueta link de html) y guarda
    # su texto (en concreto, para cada resultado, la etiqueta link correspondiente a
    # GenBank Accession es la [2]).
    def get_genbank_accession(self):
        links = scraper.soup.find_all('a', href=True)
        try:
            genbank_accession = links[2].font.text
        except AttributeError:
            genbank_accession = None
        return genbank_accession

# Clase Target:
class Target():
    
    # Constructor:
    def __init__(self):
        
        # Atributos:
        self.sequence = None
        self.mirna_name = None
        self.score = None
        self.number_of_seeds = None
        self.gene_symbol = None
        self.genbank_accession = None

# Clase Database:
class Database():
    
    # Constructor:
    def __init__(self, database: str):
        
        # Atributos:
        # SQLite nos permite crear una base de datos para el almacenamiento
        # interno de datos. Utilizamos sqlite3.connect() para generar una base
        # de datos y crear una conexion para poder acceder a ella. 
        self.__connection = sqlite3.connect(database)
        self.__create_tables()

    # Se crea la propiedad connection para guardar la conexion a la base de
    # datos. Indicamos que devuelve un tipo sqlite3.Connection. 
    # En el constructor estamos inicializando __connection y aqui devolvemos
    # el valor de __connection para poder usarlo como connection.
    @property
    def connection(self) -> sqlite3.Connection:
        return self.__connection
    
    # Metodo privado __execute_query: ejecuta la consulta. Indicamos que devuelve una 
    # lista:
    def __execute_query(self, *args) -> list:
        try:
            # Inicializamos el cursor de la base de datos (conexion entre la 
            # base de datos y la consulta).
            cursor = self.connection.cursor()
            # Ejecutamos la consulta que nos llega en la variable *args.
            cursor.execute(*args)
            self.connection.commit()
            # Extraemos los resultados de la consulta.
            result = cursor.fetchall()
        # Finally se ejecuta siempre despues de un try o un except.
        finally:
            # Cerramos el cursor (obligatorio).
            cursor.close()
        return result

    # Metodo privado __create_tables: para crear la tabla de la base de datos.
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

    # Metodo insert_target: para insertar los datos en las columnas de la tabla.
    def insert_target(self, target: Target):
        self.__execute_query("""INSERT INTO prueba_final_2 (sequence, mirna, score, seeds, GeneSymbol, GenBankAccession)
                                    VALUES (?,?,?,?,?,?)""",
                            [target.sequence,target.mirna_name,target.score,target.number_of_seeds,target.gene_symbol,target.genbank_accession])
    
    # Metodo export_to_csv: exporta todos los datos de la tabla en un archivo csv.
    def export_to_csv(self, output_name: str):
        table = pd.read_sql_query("SELECT * FROM prueba_final_2", self.connection)
        table.to_csv(output_name, index=None, header=True)


if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='Automated miRDB custom microRNA target prediction search: '
                                                '\nThis script uses a webdriver to access the miRDB website and search'
                                                ' for mRNA targets prediction with user-provided miRNA.')
    
    # Argumentos posicionales:
    parser.add_argument('inp', type=str, help='Input FASTA file with miRNA sequences')
    parser.add_argument('out', type=str, help='Name for output file')
    parser.add_argument('sp', type=str, choices=["Human","Rat","Mouse","Chicken", "Dog"], help='Species')
    
    # Argumentos opcionales (los valores por defecto hay que cambiarlos aqui):
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
            if 17 <= len(sequence.seq) <= 30:  # Restriccion longitud miRDB
                print(f'\rSearching targets for sequence: {n+1}/{len(search.fasta)} - {sequence.id}', end='', flush=True)
                crawler.driver.get(search.url) # Cargamos la pagina web miRDB
                crawler.select_element('searchSpecies', search.species)
                crawler.select_element('subChoice', search.submission)
                crawler.enter_sequence(sequence.seq)
                crawler.continue_to_results()
                crawler.wait_results()
                html = crawler.driver.page_source # Extraemos el codigo html de la pagina
                scraper.parse(html)
                passed_cutoff_rows = scraper.get_above_cutoff(search.cutoff)

                for row in passed_cutoff_rows:
                    details = crawler.driver.find_elements('name', '.submit')  # the first is the "Return" button, the others are "Target Details"
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
                    crawler.driver.back() # Va a la pagina anterior en el navegador
            else:
                print(f'\nFailed to search {sequence.id}. Sequence length out of range ({len(sequence)} nt).')
                logging.error(f"{sequence.id} length ({len(sequence)} nt) out of range (17 - 30 nt)")
    except:
        logging.exception("Error")
        print(f'\nAn exception has occurred while searching: {sequence.id}')
    finally:
        database.export_to_csv(f'{args.out}.csv')
        print(f'\nResults saved to {args.out}.csv')
        crawler.driver.close() # Cierra el navegador
        database.connection.close # Cierra la conexion con la base de datos
