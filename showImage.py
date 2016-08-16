from astropy.io import fits
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import numpy as np
def main(): 


	#fileName = sys.argv[1]
	fileName = 'horseshoe_test/HorseShoe_large_ADU.fits'
	fileName = 'sie_test/img_mod_0.fits'
	fileName = 'sie_test/img_src_0_2.fits'
	fileName = 'sie_test/sie_test.fits'
	fileName = 'sie_test/img_res_0.fits'
	hdulist = fits.open(fileName)
	data = hdulist[0].data
	print data.shape
	print len(data)
	# reverse the data; 
	n = len(data)

	n = data.shape[0]
	m = data.shape[1]
	new_data = np.zeros(data.shape)

	for i in range(n): 	
		new_data[i] = data[n-1-i]
		
	## cut off the image; 
	cut_ratio = 0/3
	new_n = int(n*(1-2*cut_ratio))
	new_m = int(m*(1-2*cut_ratio))
	cut_data = np.zeros((new_n,new_m))
	for i in range(new_n):
		for j in range(new_m):
			cut_data[i][j] = new_data[i+cut_ratio*n][j+cut_ratio*m]-10

	# img1: [-1.4, 0.0]
	# img2: [-0.5, 1.0]
	ax = plt.gca()
	im = ax.imshow(cut_data,clim=(0, 3))
	plt.axis("off")
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)

	plt.colorbar(im, cax=cax)
	plt.show()




if __name__=="__main__":  
	main()