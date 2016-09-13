# Export data frame to LaTeX table 

export_table = function (filename='', data, aligns='') {
    # If aligns is not set, by default just center every column
    if (aligns == '') {
        aligns = paste0('|', paste(rep('c', ncol(data)), collapse='|'), '|')
    }
    
    # Generate latex header
    f = file(filename, open='w')
    writeLines(c('\\begin{table}[ht]', '\t\\centering', 
                 paste0('\t\\begin{tabular}{', aligns, '}'),
                 '\t\t\\hline\\hline'), f)
    
    # Generate table header
    writeLines(paste0('\t\t', paste(names(data), collapse=' & '), ' \\\\'), f)
    writeLines('\t\t\\hline\\hline', f)
    invisible(apply(data, 1, function (row) {
        writeLines(paste0('\t\t', paste(row, collapse=' & '), ' \\\\'), f)
    }))
    writeLines(c('\t\t\\hline', '\t\\end{tabular}', '\\end{table}'), f)
    close(f)
}