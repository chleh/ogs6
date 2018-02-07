#!/usr/bin/python

import inspect

class test:
    def __init__(self):
        print("test.__init__")
        # print(inspect.getsource(self.test))
        ls = []
        for i, l in enumerate(inspect.getsourcelines(self.test)[0][1:]):
            if i == 0:
                indent = len(l) - len(l.lstrip())
            ls.append(l[indent:])
        self.Script = "".join(ls)

    def test(self):
        print(x)

class test2(test):
    # def __init__(self):
    #     print("test2.__init__")
    #     # super(test2, self).__init__()

    def test(self):
        print(y)

t = test2()

x = { 2: 3, "4": False, "5\n5": 1e-16 }
for i in range(50):
    x[i] = i+1
# import pprint
# print(pprint.pformat(x).split("\n"))

x = {0: 1,
 1: 2,
 2: 3,
 3: 4,
 4: 5,
 5: 6,
 6: 7,
 7: 8,
 8: 9,
 9: 10,
 10: 11,
 11: 12,
 12: 13,
 13: 14,
 14: 15,
 15: 16,
 16: 17,
 17: 18,
 18: 19,
 19: 20,
 20: 21,
 21: 22,
 22: 23,
 23: 24,
 24: 25,
 25: 26,
 26: 27,
 27: 28,
 28: 29,
 29: 30,
 30: 31,
 31: 32,
 32: 33,
 33: 34,
 34: 35,
 35: 36,
 36: 37,
 37: 38,
 38: 39,
 39: 40,
 40: 41,
 41: 42,
 42: 43,
 43: 44,
 44: 45,
 45: 46,
 46: 47,
 47: 48,
 48: 49,
 49: 50,
 '4': False,
 '5\n5': 1e-16}
