### import the library (should accessible through PYTHON_PATH)
import sac_de_billes
import time

L = [1, 1, 1]
seed = 0
desiredRPhi = [[0.01, 0.55]]
tabPhases = [0]
mindist = 0

tic = time.time()

sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, seed, desiredRPhi, tabPhases, mindist)
sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.WP, sac_de_billes.NameShape.Tore, L, seed, desiredRPhi, tabPhases, mindist)

print(time.time()-tic)
