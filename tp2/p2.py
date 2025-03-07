from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

mic_idxs = [1, 2, 3]

path = Path(f"tp2/data/")

tif_files = [path / f"Micrographie_{i}.tif" for i in mic_idxs]

images = [plt.imread(f) for f in tif_files]

# fft of each image

for im in images[:]:

    plt.imshow(im)
    plt.show()

    fft = np.fft.fft2(im)
    fft = np.fft.fftshift(fft)

    # Apply a Gaussian filter in Fourier space
    sigma = 2  # Adjust sigma to control filtering strength
    filtered_fft = gaussian_filter(np.abs(fft), sigma=sigma)


    # Display results
    plt.figure(figsize=(10,5))
    plt.subplot(1,2,1)
    plt.title("Original Diffraction Pattern")
    plt.imshow(np.log(1 + np.abs(fft)), cmap='gray')

    plt.subplot(1,2,2)
    plt.title("Filtered Diffraction Pattern")
    plt.imshow(np.log(1 + np.abs(filtered_fft)), cmap='gray')
    plt.show()

    ### TODO : Link Pattern to what is requested in the TP

    """

    p1 = (137, 371)
    p2 = (888, 657)
    range_around_p = 1

    # only keep fft pixels with a square around p1 and p2
    mask = np.zeros_like(fft)
    # mask[p1[0]-range_around_p:p1[0]+range_around_p, p1[1]-range_around_p:p1[1]+range_around_p] = 1
    # mask[p2[0]-range_around_p:p2[0]+range_around_p, p2[1]-range_around_p:p2[1]+range_around_p] = 1

    mask[p1[0], p1[1]] = 1
    mask[p2[0], p2[1]] = 1
    
    masked = fft * mask

    #plt.subplot(121)
    #plt.imshow(normalized)
    #plt.imshow(mask)
    #plt.subplot(122)
    #plt.imshow(masked)
    #plt.show()


    ifft = np.fft.ifft2(np.fft.ifftshift(masked))

    plt.imshow(np.abs(ifft), cmap='gray')
    plt.show()

    """