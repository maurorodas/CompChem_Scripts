import re
class inputData:
    def __init__(self, totalProts, ligandName, spacing, method, adtPath):
        self.totalProts = totalProts
        self.pdbcodes = []
        self.chains = []
        self.ligandName = ligandName
        self.boxCenter = []
        self.boxSize = []
        self.spacing = spacing
        self.method = method
        self.adtPath = adtPath
        
    def add_pdb(self, code):
        self.pdbcodes.append(code)
    
    def add_chain(self, chain):
        self.chains.append(chain)
        
    def add_boxCenter(self, pos):
        self.boxCenter.append(pos)
        
    def add_boxSize(self, size):
        self.boxSize.append(size)
    
inputFile = open("parameters.inp", "r")

for line in inputFile:
    line = line.strip()
    if re.match("proteinNumber", line):
        print(line.split()[1])
    
# data = inputData( 5 )    
# data.add_pdb('1FM6')
# data.add_pdb('3CS8')
# data.add_chain('A')
# data.add_chain('B')
# data.add_pdb('3DZY')
# data.add_chain('D')

# for x in range(len(data.pdbcodes)):
#     print(data.pdbcodes[x], data.chains[x], sep=' ')
