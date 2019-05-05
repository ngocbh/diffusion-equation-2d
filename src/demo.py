import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

def file_to_matrix(path):
	matrix = []
	with open(path) as f:
		for line in f:
			submatrix = [float(x) for x in line.split(" ")]
			matrix.append(submatrix)
	return matrix

def draw_graph(matrix, m = 20, n = 20):
	# Create figure
	fig = plt.figure(figsize=[12,8])
	ax = fig.add_subplot(111)
	div = make_axes_locatable(ax)
	cax = div.append_axes('right', '5%', '5%')
	ax.cla()
	cax.cla()
	ax.set_xlabel("$x$ Position", labelpad=15, size=18)
	ax.set_ylabel("$y$ Position", labelpad=15, size=18)
	ax.set_aspect('equal')
	
	xgrid = np.linspace(0, 1, m)
	ygrid = np.linspace(0, 1, n)
	X, Y = np.meshgrid(xgrid, ygrid)
	Z = matrix
	sol = ax.contourf(X, Y, Z, cmap= "plasma")
	clb = fig.colorbar(sol, cax=cax, spacing='proportional')
	clb.ax.set_ylabel("Density", labelpad=15, size=18)
	plt.show()
	
def main():
	matrix = []
	for i in range(20):
		submatrix = []
		for j in range(20):
			submatrix.append(j)
		matrix.append(submatrix)
	draw_graph(matrix)
	
if __name__ == '__main__':
	main()