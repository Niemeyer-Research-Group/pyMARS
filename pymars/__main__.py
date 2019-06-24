"""Main module"""
import sys

from .pymars import pymars


def main(args=None):
    if args is None:
        args = sys.argv[1:]
        pymars(args)


if __name__ == '__main__':
    sys.exit(main())
