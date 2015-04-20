import multiprocessing

def f(x):
    created = multiprocessing.Process()
    current = multiprocessing.current_process()
    print x, 'running:', current._identity[0] # current.name,
    #print x, 'created:', created.name, created._identity
    return x * x

var = range(3)#[1,2,3,5,6,6,6,6,6,5,4,36,,2,1]

p = multiprocessing.Pool(2)
print p.map(f, var)
