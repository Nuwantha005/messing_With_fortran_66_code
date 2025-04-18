import numpy as np
import f1b1
print(f1b1.fib.__doc__)

a = np.zeros(8, 'd')
f1b1.fib(a)
print(a)