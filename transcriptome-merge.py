c0=pd.read_table('c0.txt',header=None)
c1=pd.read_table('c1.txt',header=None)
l1=pd.read_table('l1.txt',header=None)
l2=pd.read_table('l2.txt',header=None)
e1=pd.read_table('e1.txt',header=None)
e2=pd.read_table('e2.txt',header=None)


for i in range(5377):
    for k in range(5377):
        if data.iat[i,0]==h1.iat[k,0]:
            data.iat[i,7]=h1.iat[k,6]
            continue
