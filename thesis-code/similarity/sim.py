import Image
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def ssd(r1, r2):
    rd = r1-r2
    return np.sum(rd*rd)/np.size(rd)


def ncc(r1, r2):
    mr1 = np.mean(r1)
    mr2 = np.mean(r2)
    sr1 = np.std(r1)
    sr2 = np.std(r2)
    return np.sum((r1 - mr1)*(r2-mr2))/(sr1*sr2)/np.size(r1)


def ssd_test():
    img = Image.open("/Prime/Thesis-Data/SimilarityMetric/UNC_train_Case05_T1_sagittal_slice_0265.png").convert("L")
    arr = 1.0*np.array(img)
    (nx, ny) = arr.shape

    values = []
    for j in range(1, nx/2, 1):
        r1 = np.squeeze(arr[:, (nx/2-j):])
        r2 = np.squeeze(arr[:, 0:(nx/2+j)])
        values.append(ssd(r1, r2))

    values.extend(values[-2::-1])
    plt.plot(np.arange(-255,254),values)
    plt.ylabel('SSD')
    plt.xlabel('Displacement')

def ssd_test2():
    img = Image.open("/Prime/Thesis-Data/SimilarityMetric/UNC_train_Case05_T1_sagittal_slice_0265.png").convert("L")
    arr = 1.0*np.array(img)
    arr2 = -1.0*arr-2
    (nx, ny) = arr.shape

    values = []
    for j in range(1, nx/2, 1):
        r1 = np.squeeze(arr[:, (nx/2-j):])
        r2 = np.squeeze(arr2[:, 0:(nx/2+j)])
        values.append(ssd(r1, r2))
    values.extend(values[-2::-1])
    plt.plot(np.arange(-255,254),values, 'g')
    plt.ylabel('SSD')
    plt.xlabel('Displacement')


def ncc_test():
    img = Image.open("/Prime/Thesis-Data/SimilarityMetric/UNC_train_Case05_T1_sagittal_slice_0265.png").convert("L")
    arr = 1.0*np.array(img)
    (nx, ny) = arr.shape

    values = []
    for j in range(1, nx/2, 1):
        r1 = np.squeeze(arr[:, (nx/2-j):])
        r2 = np.squeeze(arr[:, 0:(nx/2+j)])
        values.append(-np.abs(ncc(r1,r2)))
    values.extend(values[-2::-1])
    plt.plot(np.arange(-255,254), values, 'bo')
    plt.ylabel('NCC')
    plt.xlabel('Displacement')


def ncc_test2():
    img = Image.open("/Prime/Thesis-Data/SimilarityMetric/UNC_train_Case05_T1_sagittal_slice_0265.png").convert("L")
    arr = 1.0*np.array(img)
    arr2 = -1.0*arr +254
    (nx, ny) = arr.shape

    imgx = Image.fromarray(arr2.astype('uint8'))
    imgx.save("/Prime/Thesis-Data/SimilarityMetric/UNC_train_Case05_T1_sagittal_slice_0265_inv.png")
    values = []
    for j in range(1, nx/2, 1):
        r1 = np.flatten(arr[:, (nx/2-j):])
        r2 = np.flatten(arr2[:, 0:(nx/2+j)])
        values.append(-np.abs(ncc(r1,r2)))
    values.extend(values[-2::-1])
    plt.plot(np.arange(-255,254),values, 'g+')
    plt.ylabel('NCC')
    plt.xlabel('Displacement')


def ssd_test_mm():
    img1 = Image.open("/Prime/Thesis-Data/SimilarityMetric/UNC_train_Case05_T1_sagittal_slice_0265.png").convert("L")
    img2 = Image.open("/Prime/Thesis-Data/SimilarityMetric/UNC_train_Case05_T2_sagittal_slice_0265.png").convert("L")
    arr1 = 1.0*np.array(img1)
    arr2 = 1.0*np.array(img2)

    (nx, ny) = arr1.shape

    values = []
    for j in range(1, nx/2, 1):
        r1 = np.squeeze(arr1[:, (nx/2-j):])
        r2 = np.squeeze(arr2[:, 0:(nx/2+j)])
        values.append(ssd(r1, r2))

    values.extend(values[-2::-1])
    plt.plot(np.arange(-255,254),values, 'b')
    plt.ylabel('SSD')
    plt.xlabel('Displacement')


def ncc_test_mm():
    img1 = Image.open("/Prime/Thesis-Data/SimilarityMetric/UNC_train_Case05_T1_sagittal_slice_0265.png").convert("L")
    img2 = Image.open("/Prime/Thesis-Data/SimilarityMetric/UNC_train_Case05_T2_sagittal_slice_0265.png").convert("L")
    arr1 = 1.0*np.array(img1)
    arr2 = 1.0*np.array(img2)

    (nx, ny) = arr1.shape

    values = []
    for j in range(1, nx/2, 1):
        r1 = np.squeeze(arr1[:, (nx/2-j):])
        r2 = np.squeeze(arr2[:, 0:(nx/2+j)])
        values.append(-ncc(r1, r2))

    values.extend(values[-2::-1])
    plt.plot(np.arange(-255,254),values, 'g+')
    plt.ylabel('NCC')
    plt.xlabel('Displacement')


def histo2d(r1, r2):
    nBin = 8
    h1 = r1*nBin/255.0
    h2 = r2*nBin/255.0
    plt.hist2d(h1, h2, bins=nBin+1, norm=LogNorm())
    plt.colorbar()
    plt.show()


def compute_MI(data):
    hx = np.sum(data,axis=0)
    hy = np.sum(data,axis=1)
    hh = np.sum(data)
    (nx,ny) = data.shape
    mi = 0
    for j in range(0,nx):
        for k in range(0,ny):
            if hx[j] > 0 and hy[k] > 0:
                px = data[j,k]/hx[j]
                py = data[j,k]/hy[k]
                pxy = data[j,k]/hh
                if (px>0 and py> 0 and pxy > 0):
                    mi += pxy * np.log(pxy/px*py)
                    if np.isnan(mi):
                        print "NaN", mi
    return mi


def mi_test():
    img = Image.open("/Prime/Thesis-Data/SimilarityMetric/UNC_train_Case05_T1_sagittal_slice_0265.png").convert("L")
    arr = 1.0*np.array(img)
    (nx, ny) = arr.shape
    values = []
    j = 0
    r1 = arr[:, (nx/2-j):].flatten()
    r2 = arr[:, 0:(nx/2+j)].flatten()
    #histo2d(r1,r2)
    nBins = 8
    for j in range(1, nx/2, 8):
        print j
        r1 = arr[:, (nx/2-j):].flatten()
        r2 = arr[:, 0:(nx/2+j)].flatten()
        (H, xedges, yedges) = np.histogram2d(r1, r2, bins=nBins)
        mi = compute_MI(H)
        print mi
        values.append(mi)
        if j == 1:
            plt.figure()
            plt.hist2d(r1, r2, bins=nBins, norm=LogNorm())
            plt.colorbar()
            plt.savefig("/Prime/Thesis-Data/SimilarityMetric/mitest1.pdf")
        if nx/2-j<8:
            plt.figure()
            plt.hist2d(r1, r2, bins=nBins, norm=LogNorm())
            plt.colorbar()
            plt.savefig("/Prime/Thesis-Data/SimilarityMetric/mitest2.pdf")
    values.extend(values[-2::-1])
    plt.plot(values, 'g+')
    plt.ylabel('MI')
    plt.xlabel('Displacement')


def main():
    # ssd_test()
    # ssd_test2()
    # plt.savefig("/Prime/Thesis-Data/SimilarityMetric/ssdtest.pdf")
    # plt.figure()
    # ncc_test()
    mi_test()
    plt.savefig("/Prime/Thesis-Data/SimilarityMetric/mitest.pdf")
    # ssd_test_mm()
    # plt.savefig("/Prime/Thesis-Data/SimilarityMetric/ssdtest_mm.pdf")
    # plt.figure()
    # ncc_test_mm()
    # plt.savefig("/Prime/Thesis-Data/SimilarityMetric/ncctest_mm.pdf")

if __name__ == "__main__":
    main()