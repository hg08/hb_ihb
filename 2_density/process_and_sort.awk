BEGIN { OFS = " " }
NR == 1 { print; next }
{
    # Modify the first field
    $1 = (($1 + 8.0) * 10000) % 312808 / 10000;
    # Collect all lines and keys
    lines[NR] = $0;
    keys[NR] = $1;
}
END {
    # Print header
    print lines[1];
    # Sorting algorithm (Bubble sort for simplicity)
    for (i = 2; i <= NR; i++) {
        for (j = 2; j <= NR - i + 1; j++) {
            if (keys[j] > keys[j + 1]) {
                # Swap keys
                tmp_key = keys[j];
                keys[j] = keys[j + 1];
                keys[j + 1] = tmp_key;
                # Swap lines
                tmp_line = lines[j];
                lines[j] = lines[j + 1];
                lines[j + 1] = tmp_line;
            }
        }
    }
    # Print sorted lines
    for (i = 2; i <= NR; i++) {
        print lines[i];
    }
}

