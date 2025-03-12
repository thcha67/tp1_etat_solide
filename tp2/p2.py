from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max

mic_idxs = [2, 3, 4]

path = Path(f"tp2/data/")

tif_files = [path / f"Micrographie_{i}.tif" for i in mic_idxs]

images = [plt.imread(f) for f in tif_files]

# fft of each image

for im in images[:]:

    plt.imshow(im)
    plt.show()

    fft = np.fft.fft2(im)
    fft = np.fft.fftshift(fft)

    # Cela converti un pixel dans la transfo de fourrier en mm-1. La distance entre 2 point est désormais en mm-1. 
    spatial_res = 42.468 * 10**(-12)
    freq_res = 1/(1024 * spatial_res)    
    freq = np.fft.fftfreq(1024, d=spatial_res)
    freq = np.fft.fftshift(freq)

    fft_magnitude = np.abs(fft)
    sigma = 4
    filtered_fft = gaussian_filter(fft_magnitude, sigma=sigma)
    fft_log = np.log(filtered_fft + 1)

    
    COMPARISON_FFT = False
    if COMPARISON_FFT:
        # Display results
        plt.figure(figsize=(10,5))
        plt.subplot(1,2,1)
        plt.title("Original Diffraction Pattern")
        plt.imshow(np.log(1 + np.abs(fft)), cmap='gray')

        plt.subplot(1,2,2)
        plt.title("Filtered Diffraction Pattern")
        plt.imshow(np.log(1 + np.abs(filtered_fft)), cmap='gray')
        plt.show()

    # Find local maxima in the FFT image (spots)
    # Use a threshold to avoid detecting noise
    threshold = np.percentile(filtered_fft, 99.2)  # Can adjust this as needed
    spots = peak_local_max(filtered_fft, min_distance=20, threshold_abs=threshold)

    # Display the identified spots
    plt.figure(figsize=(6, 6))
    plt.imshow(fft_log, cmap='gray')
    plt.scatter(spots[:, 1], spots[:, 0], color='red', s=20)  # Plot spots
    # Annotate each point with its coordinates
    for i, (y_val, x_val) in enumerate(spots):
        print(f"Noeud[{i}] = ({freq[x_val]:.2e} mm-1,{freq[y_val]:.2e} mm-1)")
        plt.text(x_val, y_val - 15, f"({freq[x_val]:.2e},{freq[y_val]:.2e})", fontsize=12 , color='red', ha = "center")
        pass
    plt.title("Reciprocal Lattice")
    plt.xlabel("Pixels [-]", fontsize = 14)
    plt.ylabel("Pixels [-]", fontsize = 14)
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


    [228.47, 228.47, 228.39, 228.29, 228.47, 228.47]

    [136.71, 137.72, 135.62, 136.71, 136.71, 135.62]

    """