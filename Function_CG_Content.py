def cg_content(Seq):
    c = Seq.count('C')  # get the number of C in the sequence
    g = Seq.count('G')  # get the number of G in the sequence
    cg_content = round((c + g) / len(Seq) * 100, 1)  # Rounding to first decimal
    return cg_content
