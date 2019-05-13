import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation
import sys

def file_to_matrix(path = 'out.txt', m = 20):
	list_of_matrix = []
	matrix = []
	i = 0
	with open(path) as f:
		for line in f:
			i += 1
			submatrix = [float(x) for x in line.split()]
			matrix.append(submatrix)
			if i % m == 0:
				list_of_matrix.append(matrix)
				matrix = []
		list_of_matrix.append(matrix)
	return list_of_matrix

def make_video(list_of_matrix, m = 20, n = 20, name = 'demo.mp4'):
	fig = plt.figure(figsize=[12,8])
	ax = fig.add_subplot(111)
	div = make_axes_locatable(ax)
	cax = div.append_axes('right', '5%', '5%')
	nFrames = len(list_of_matrix) - 1
	def draw_graph(i, list_of_matrix, m, n):
		xgrid = np.linspace(0, 1, m)
		ygrid = np.linspace(0, 1, n)
		X, Y = np.meshgrid(xgrid, ygrid)
		Z = list_of_matrix[i]
		sol = ax.contourf(X, Y, Z, 1000, cmap= "plasma")
		if i < 51:
			ax.set_title("n = {}".format(i))
		elif i != nFrames - 1:
			ax.set_title("n = {}".format(20*(i-50)))
		else:
			ax.set_title("n = {}".format(list_of_matrix[i+1][0][0]))
		ax.set_xlabel("$l$", labelpad=15, size=18)
		ax.set_ylabel("$m$", labelpad=15, size=18)
		ax.set_aspect('equal')
		ticks = np.linspace(25.0, 80.0, num=11)
		clb = fig.colorbar(sol, cax=cax, ticks = ticks, spacing='proportional')
		clb.ax.set_ylabel("Concentration", labelpad=15, size=18)
		if i == 0:
			plt.savefig('Before100')
		if i == nFrames - 1:
			plt.savefig('After100')
			sys.exit(0)
	ani = animation.FuncAnimation(fig, draw_graph, fargs=(list_of_matrix,m,n), blit=False,
                              repeat=False, save_count=0, frames=nFrames)
	ani.save(name, writer='ffmpeg', dpi=180, fps=5)

def main():
	#After
	m = 100
	n = 100
	fig = plt.figure()
	list_of_matrix = file_to_matrix('out.txt',m)
	print(len(list_of_matrix))
	make_video(list_of_matrix, m, n, 'demo100_3.mp4')
if __name__ == '__main__':
	main()