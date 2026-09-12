def range_to_table(workbook_path):
    from openpyxl import load_workbook
    from openpyxl.styles import Font
    from openpyxl.worksheet.table import Table, TableStyleInfo
    from openpyxl.styles import Alignment

    #^################################################################################
    #^ Grab the table number of all the tables in all the worksheets:
    #^################################################################################
    workbook      = load_workbook(workbook_path)
    table_numbers = list()
    for sheet in workbook.worksheets:
        if len(sheet.tables) != 0:
            for table_name in sheet.tables:
                table_number = int(table_name[5:]) # "table_name" variable will store values like "Table1", "Table2" & so on. I am grabbing only the integer part of it
                table_numbers.append(table_number)
    table_numbers = sorted(table_numbers)

    #^################################################################################
    #^ Add table to those sheets which don't have one:
    #^################################################################################
    if len(table_numbers) == 0: # if the workbook had no previous tables, then ctr should be equal to 0
        ctr = 0
    if len(table_numbers) != 0:
        ctr = table_numbers[-1]
    for sheet in workbook.worksheets:
        ctr += 1
        if len(sheet.tables) == 0: # if a sheet has no tables
            print("converting the following worksheet to table:", sheet.title)
            table_range          = sheet.dimensions
            table                = Table(displayName = "Table" + str(ctr), ref = table_range)
            style                = TableStyleInfo(name = 'TableStyleLight8', showFirstColumn = False, showLastColumn = False, showRowStripes = True, showColumnStripes = True)
            table.tableStyleInfo = style
            sheet.add_table(table)

        #^######################################################################################
        #^ Adjust the width of each column of the current sheet's table to best readable fit:
        #^######################################################################################
        for column_cells in sheet.columns:
            list_widths   = [len(str(cell.value)) for cell in column_cells] 
            max_width     = max(list_widths)
            best_width    = (max_width + 3) * 1.2          # to get the best width for the column, add a bit of padding and adjust the width slightly
            column_letter = column_cells[0].column_letter  # grab the column letter of the current column e.g., 'A', 'B', ...
            sheet.column_dimensions[column_letter].width = best_width

        # left align each cell of each row (except colheader row) in each column
        for row_index, row in enumerate(sheet.iter_rows(min_col = 2, max_col = len(list(sheet.columns))), 1): #"min_col" and "max_col" denote column numbers (1 = A, 2 = B and so on) to iterate over
            if row_index != 1:
                for cell in row:
                    cell.alignment = Alignment(horizontal = "left")

    #^################################################################################
    #^ Change the font color of the topmost row to white:
    #^################################################################################
    for sheet in workbook.worksheets:
        for cell in sheet[1]:                  # "sheet[1]" corresponds to 1st row of the worksheet
            cell.font = Font(color = "FFFFFF") # change the font color of the topmost row to white
    workbook.save(workbook_path)