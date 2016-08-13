from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import numpy as np
def main(): 

	hdulist = fits.open('horseshoe_test/HorseShoe_large.fits')
	data = hdulist[0].data
	

	# reverse the data; 
	n = len(data)
	new_data = np.zeros(data.shape)

	for i in range(n): 
		
		new_data[i] = data[n-1-i]
		
		
	# img1: [-1.4, 0.0]
	# img2: [-0.5, 1.0]
	ax = plt.gca()
	im = ax.imshow(new_data, clim=(-0.02, 0.16))
	plt.axis("off")
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)

	plt.colorbar(im, cax=cax)
	plt.show()




if __name__=="__main__":  
	main()