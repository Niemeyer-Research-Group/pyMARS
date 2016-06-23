
import argparse
class test(object):
    def __init__(self, var1, var2):
        self.var1 = var1
        self.var2 = var2
    variable=2


def convert(a, b='none'):
    value=a
    other_file=b
    print value
    print other_file



def func( var='none', *arg, **argv):
        print var
        print arg


        if var =='none':
            print argv
            if 'thermo' in argv:
                thermo_file = argv['thermo']
            if 'transport' in argv:
                transport_file = argv['transport']
            if 'species' in argv:
                species = argv['species']
            else:
                species = []

input1='gri30.cti'
class arguments:
    fileinput=input1
args=argparse.Namespace()
a={'plot': 'false'}
args.plot = 'True'
print args
