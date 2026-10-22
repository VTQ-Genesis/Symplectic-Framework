from quantum_anticodes import *
from galois import GF, FieldArray
import numpy as np

_BINARY_FIELD = GF(2)

gx = _BINARY_FIELD([[1,1,1,1], [1,1,0,0]])
gz = _BINARY_FIELD([[1,1,0,0], [0,0,1,1],[1,0,1,0]])

c = choi(gx, gz)

print(c)