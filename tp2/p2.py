from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from findpeaks import findpeaks

mic_idxs = [1, 2, 3]

path = Path(f"tp2/data/")

tif_files = [path / f"Micrographie_{i}.tif" for i in mic_idxs]

images = [plt.imread(f) for f in tif_files]

# fft of each image

for im in images[0:1]:
    fft = np.fft.fft2(im)
    fft = np.fft.fftshift(fft)

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
    # plt.subplot(121)
    # plt.imshow(normalized)
    # plt.subplot(122)
    # plt.imshow(masked)
    # plt.show()


    ifft = np.fft.ifft2(np.fft.ifftshift(masked))

    plt.imshow(np.abs(ifft), cmap='gray')
    plt.show()