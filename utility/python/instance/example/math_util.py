#!/usr/bin/python

import sys, codecs
import time
import numpy as np

def funx2(a):
    data = [[1,2],[3,4],[5,6]]
    x = np.array(data)
    print(x)
    print("call funx2 numpy2")
    return a*2


def add_func(a, b):
    return (a + b)


def sub_func(a, b):
    return (a - b)


def mainfunc(a, b):
    starttime = time.time()
    funx2(a)
    print("success!")
    result_1 = add_func(a, b)
    result_2 = sub_func(a, b)
    print("a+b= ", result_1, " a-b= ", result_2)
    endtime = time.time()
    print("spend time: ", endtime - starttime)


def main():
    mainfunc(2, 6)


if __name__ == '__main__':
    main()