import pandas as pd

bbp_eos = pd.read_csv("/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/bbp.csv")

# rho,P,nb,A,Z,Gamma
rho = bbp_eos["rho"]
P = bbp_eos["P"]
nb = bbp_eos["nb"]
