################################################################################
# Read all sheets of excel workbooks; output is a tibble
# Source: https://stackoverflow.com/questions/8188415/save-excel-spreadsheet-as-csv-with-r
# Dependencies:
## library(readxl)
################################################################################
################################################################################
readExcelWb <- function(filepath=".xlsx",
                        colnames=FALSE) { # TRUE of vector of names
  sheets.nme <- readxl::excel_sheets(filepath)
  sheets.lst <- lapply(sheets.nme, function(x){
    #as.data.frame(readxl::read_excel(filepath, sheet=x))
    df <- readxl::read_excel(filepath, sheet=x,
                             col_names=colnames)
  
    return(data.frame(df, stringsAsFactors=FALSE))
    #return(df)
    
  })
  names(sheets.lst) <- sheets.nme
  return(sheets.lst)
}
