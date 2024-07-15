from phcpy.solver import random_trinomials
from phcpy.solver import solve

def hola():
	f = random_trinomials()
	sols = solve(f)
	for s in sols:
		print(s)