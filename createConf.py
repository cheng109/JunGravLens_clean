import numpy as np 

def createConf(confName): 
    header= """
## Setting primary parameters: 
omega     0.28600
lambda	  0.70000
weos	  -1.000000
hubble	  0.696000
srcZ	2.379
lenZ	0.4457


## Setting secondary parameters: 
srcRes	0.1
imgRes	0.04
potRes	0.1

## Files
imageFileName	horseshoe_test/HorseShoe_large.fits
regionFileName	horseshoe_test/mask_horseshoe_large_4.reg
varFileName	  jun_var.fits
psfFileName		horseshoe_test/HorseShoe_psf.fits
criticalName	horseshoe_test/critical.reg
causticName	horseshoe_test/caustic.reg


## Source plane size
srcX	100
srcY	100
"""

    f = open(confName, 'w')
    f.write(header)
    xRange = np.arange(-150, 150, 3)
    yRange = np.arange(-150, 150, 3)
    critRange = np.arange(0.1, 5.5, 0.1)

    PTMASS_critFrom = 0.0
    PTMASS_critTo   = 5.5
    PTMASS_critInc  = 0.1

    for x in xRange:
        for y in yRange:
            f.write("PTMASS "+ str(x)       + " " + str(x)       + " " + "0.1 "
                             + str(y)       + " " + str(y)       + " " + "0.1 "
                             + str(PTMASS_critFrom) + " " + str(PTMASS_critTo) + " " + str(PTMASS_critInc)
                             + "\n")


    f.close()
    print len(xRange)*len(yRange)*len(critRange)


def main():

    confName="horseshoe_test/conf_test.txt"

    createConf(confName)





if __name__=="__main__": 
    main()