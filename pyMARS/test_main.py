

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



def func(file_name, var='none', **sys_args):
    print file_name
    print var
    if 'plot' in sys_args:
        print 'plotting'
