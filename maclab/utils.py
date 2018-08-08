import re, operator

digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))

def readfile(filepath):
    if os.path.exists(filepath):
        output = subprocess.run(["cat", filepath],stdout=subprocess.PIPE).stdout.decode('utf-8')
        print(output)
    else:
        print('The filepath does not exist')

def which(x, op, val):
    mappings = {'>': operator.lt, '>=': operator.le,
                '==': operator.eq, '<': operator.gt,
               '<=': operator.ge} 
    return([i for i, x in enumerate(x) if mappings[op](x, val) ])

def asnumber(number, dtype = 'float32'):
    formtype = '{0:.%sf}' %dtype[-2:]
    print(formtype.format(number))
    
def get(x):
    """Equiavlent of R get? Get name of variable as string"""
    return([ k for k,v in locals().items() if v is x][0])