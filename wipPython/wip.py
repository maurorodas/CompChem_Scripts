class inputData:
    def __init__(self, totalProts):
        self.totalProts = totalProts
        self.pdbcodes = []
        self.chains = []

    def add_pdb(self, code):
        self.pdbcodes.append(code)
    
    def add_chain(self, chain):
        self.chains.append(chain)
    
data = inputData( 5 )    
data.add_pdb('1FM6')
data.add_pdb('3CS8')
data.add_chain('A')
data.add_chain('B')

for x in range(len(data.pdbcodes)):
    print(data.pdbcodes[x], data.chains[x], sep=' ')
