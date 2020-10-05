import string
import xlrd

class Taxon:
    """
    A class that stores all the data concerning a taxon in the target Excel file.
    """

    def __init__(self, name, charset, nchars, LAD, FAD):
        self.name = name
        self.charset = charset
        self.nchars = nchars
        self.LAD = LAD
        self.FAD = FAD

def load_data(loc, sheet_index):
    """
    Loads a spreadsheet at the given index which contains all the coded data.
    """

    wb = xlrd.open_workbook(loc)
    sheet = wb.sheet_by_index(sheet_index)
    return sheet

def get_taxon_for_row(sheet, row_index):
    """
    When given a row on a spreadsheet, returns a Taxon instance corresponding to that row of the spreadsheet.
    """

    name = sheet.cell_value(row_index, 0) # Get and clean taxon name (many characters are not allowed in the names of taxa).
    name = clean_taxon_name(name)

    charset = "" # Get and clean taxon characters so that the data can be read correctly by MrBayes.
    col_index = 1
    while col_index + 3 <= sheet.ncols:
        cell_value = sheet.cell_value(row_index, col_index)
        cleaned_cell_value = clean_cell_value(cell_value)
        charset = charset + cleaned_cell_value
        col_index += 1 # Final incrementation means that col index is now for LAD column.

    nchars = sheet.ncols - 3 # -1 for the 0th col (names), -2 for FAD, LAD cols.

    LAD = str(int(sheet.cell_value(row_index, col_index)))
    FAD = str(int(sheet.cell_value(row_index, col_index + 1)))

    row_taxon = Taxon(name, charset, nchars, LAD, FAD) # Since col_index started from 1, also serves to count the number of characters for the taxon.
    return row_taxon

def get_taxon_from_rows(sheet, row_start, row_end):
    """
    Returns an array of taxa from a spreadsheet from the specified rows.
    """

    taxa_arr = []
    for i in range(row_start, row_end):
        taxa_arr.append(get_taxon_for_row(sheet, i))
    return taxa_arr

def generate_nexus_file(sheet, row_start, row_end):
    """
    When given a range of rows, generates a nexus file of taxa data corresponding to those rows. Includes MrBayes code at the end.
    """

    taxa_arr = get_taxon_from_rows(sheet, row_start, row_end)
    ntaxa = len(taxa_arr)
    nchars = taxa_arr[0].nchars
    matrix = generate_matrix(taxa_arr)
    mrbayes = generate_mrbayes(taxa_arr)

    f = open("generated_file.nex", "w+")
    f.write("#NEXUS\n\nBEGIN DATA;\n")
    f.write(f"    DIMENSIONS NTAX={ntaxa} NCHAR={nchars};\n    FORMAT Datatype=Standard Symbols=\"0123456\" Missing=? Gap=-;\n")
    f.write(matrix)
    f.write(mrbayes)
    f.close
    return

def clean_cell_value(cell_value):
    """
    Ensures that the cell data is formatted correctly so that it can be added to the respective taxon data and read by MrBayes.
    """

    coded_chars = ["?", "-"]

    if cell_value in coded_chars or isinstance(cell_value, float): # If the cell only contains one coded character, we can just return it as a string
        if isinstance(cell_value, float):
            return str(int(cell_value))
        elif isinstance(cell_value, str):
            return cell_value
        else:
            print(cell_value)
            print("Error cleaning this cell character")
    
    else: # If we've gotten here it should be because the cell has more than one coded character.
        char_lst = cell_value.split(",")

        instance_counter = {"?": 0, " ?": 0} # This will handle and remove any instances of ? in the cell if it contains more than one character.
        for char in char_lst:
            if char == "?" or char == " ?":
                instance_counter[char] += 1
        for key in instance_counter.keys():
            while instance_counter[key] > 0:
                char_lst.remove(key)
                instance_counter[key] -= 1

        if len(char_lst) == 1: # If removing ? means that there is now only one character left, then just return that character.
            return char_lst[0]
        else:
            middle_chars = "".join(char_lst)
            return f"({middle_chars})"
    
def clean_taxon_name(dirty_name):
    """
    If there are any characters that aren't allowed in the taxon name then replace them with "_".
    """

    str_arr = list(dirty_name)
    whtspc_chars = string.whitespace
    disallowed_chars_counter = {"(": 0, ")": 0 , "?": 0 , "=": 0 , "+": 0, "&": 0, ";": 0, ",": 0}
    for i, char in enumerate(str_arr):
        if char in whtspc_chars:
            str_arr[i] = "_"
        elif char in disallowed_chars_counter.keys():
            disallowed_chars_counter[char] += 1
    for key in disallowed_chars_counter.keys():
        while disallowed_chars_counter[key] > 0:
            str_arr.remove(key)
            disallowed_chars_counter[key] -= 1
    return "".join(str_arr)

def generate_matrix(taxa_arr): # generates a matrix when given an array of taxa
    """
    Returns a string that can be interpreted as a matrix of the taxa data by MrBayes.
    """

    largest_name_len = 0
    for taxon in taxa_arr:
        if len(taxon.name) > largest_name_len:
            largest_name_len = len(taxon.name)
    matrix = "MATRIX\n"
    for taxon in taxa_arr: # Ensures that the matrix string looks nice when read by human eye as well.
        name_diff = largest_name_len - len(taxon.name)
        spaces = "    "
        for _ in range(name_diff):
            spaces = spaces + " "
        matrix = matrix + f"    {taxon.name}{spaces}{taxon.charset}\n"
    matrix = matrix + "    ;\nEND;\n\n"
    return matrix

def generate_mrbayes(taxa_arr):
    """
    Generates a string that can be read by MrBayes as commands for the analysis.
    """

    calibration_data = generate_calibration_data(taxa_arr)

    mrbayes = "BEGIN MrBayes;\n"
    relaxed_clock = "    [relaxed clock model]\n    prset clockvarpr = igr;\n    prset igrvarpr = exp(10);\n\n"
    tip_dating = f"    [tip dating]\n    {calibration_data}"
    mcmc = generate_mcmc(1000000)
    exe = "    mcmc;\n    sumt;\n    sump;\n"

    mrbayes = mrbayes + relaxed_clock + tip_dating + mcmc + exe + "END;"
    return mrbayes

def generate_calibration_data(taxa_arr):
    """
    Generates a string that can be read by MrBayes as calibration data for tip dating.
    """
    calibration_data = "calibrate\n"
    for taxon in taxa_arr:
        calibration_data = calibration_data + f"        {taxon.name} = unif({taxon.LAD}, {taxon.FAD})\n"
    calibration_data = calibration_data + "    ;\n    prset nodeagepr = calibrated;\n\n"
    return calibration_data

def generate_mcmc(ngen):
    """
    Generates the MCMC commands for the analysis of the data.
    """
    
    mcmc_line1 = "    [mcmc settings]\n"
    mcmc_line2 = f"    mcmcp ngen = {ngen} samplefr = {0.05 * ngen} printfr = {0.05 * ngen} diagnfr = {0.125 * ngen};\n"
    mcmc_line3 = "    mcmcp filename = \"analysis\";\n\n"
    return mcmc_line1 + mcmc_line2 + mcmc_line3