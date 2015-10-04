library(xlsx) # use this here as reading takes longer
# for bigger files

merge1 <- read.xlsx("DP_PG_PKP2vsCtr.xlsx", 1)

wb <- loadWorkbook("DP_PG_PKP2vsCtr.xlsx")

fo1 <- Fill(foregroundColor="lightgreen")   # create fill object # 1
cs1 <- CellStyle(wb, fill=fo1)        # create cell style # 1
fo2 <- Fill(foregroundColor="indianred1") # create fill object # 2
cs2 <- CellStyle(wb, fill=fo2)  

sheets <- getSheets(wb)               # get all sheets
sheet <- sheets[["oneFC>=2"]]          # get specific sheet

rows <- getRows(sheet, rowIndex=2:(nrow(merge1)+1)) # get rows

cells <- getCells(rows, colIndex = 3:8)         # get cells

values <- lapply(cells, getCellValue)

highlightgreen <- NULL
for (i in names(values)) {
  x <- as.numeric(values[i])
  if (x>2 & !is.na(x)) {
    highlightgreen <- c(highlightgreen, i)
  }    
}

# find cells meeting conditional criteria < 2
highlightred <- NULL
for (i in names(values)) {
  x <- as.numeric(values[i])
  if (x<(-2) & !is.na(x)) {
    highlightred <- c(highlightred, i)
  }    
}

lapply(names(cells[highlightgreen]),
       function(ii)setCellStyle(cells[[ii]],cs1))

lapply(names(cells[highlightred]),
       function(ii)setCellStyle(cells[[ii]],cs2))

saveWorkbook(wb, "DPPGPKP2color.xlsx")

rm(list=ls())