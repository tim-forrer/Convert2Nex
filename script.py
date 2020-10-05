import convert2nex

# Edit these lines of code in order to convert the desired spreadsheet into a .nex file.
file_loc = "/home/tim/Desktop/Plectambonitoidea (2)/data.xlsx"
sheet_index = 1
row_start_index = 3
row_end_index = 118

# Nothing beyond here needs to be changed.
sheet = convert2nex.load_data(file_loc, sheet_index)
convert2nex.generate_nexus_file(sheet, row_start_index, row_end_index + 1)