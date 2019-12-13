def print_columns(data, out_file=None):
    col_width = [max([len(data[j][i]) for j in range(len(data))]) for i in range(len(data[0]))]
    first = True
    last_row = col_width
    for row in data:
        #if not first and cat != row[0]:
        #    print("| " + "".join(" "*l + " | " for l in col_width))
        #    cat = row[0]
        s = "| " + "".join(
                (word.ljust(col_width[i]) if last_row[:i+1] != row[:i+1] \
                                          else " "*col_width[i]) + " | " for i, word in enumerate(row)
              )
        if not out_file is None:
            out_file.write(s)
            out_file.write("\n")
        print(s)
        if first:
            s2 = "-" * (sum(col_width) + len(col_width)*3 + 1)
            print(s2)
            if not out_file is None:
                out_file.write(s2)
                out_file.write("\n")
            first = False
        last_row = row