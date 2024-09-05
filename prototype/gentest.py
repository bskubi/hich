def gen():
    i = yield 5

g = gen()
next(g)