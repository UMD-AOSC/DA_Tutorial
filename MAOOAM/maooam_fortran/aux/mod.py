with open("IC_atm99_ocn66.txt") as f:
    a = f.read().split()
with open("IC_atm99_ocn66mod.txt", "w") as f:
    for c in a:
        f.write(c)
        f.write("\n")
