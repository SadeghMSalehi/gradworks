import Image
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, PowerNorm


def two_scale(t, s1, s2, fname = ""):
    plt.ticklabel_format(sytle='sci')
    fig = plt.figure(figsize=(8,4))

    ax1 = plt.gca()
    ax1.plot(t, s1, 'b-')
    ax1.set_xlabel('Translation')
    # Make the y-axis label and tick labels match the line color.
    ax1.set_ylabel('MI', color='b')
    ax1.ticklabel_format(style='sci')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')

    ax2 = ax1.twinx()
    ax2.plot(t, s2, 'r.')
    ax2.set_ylabel('NMI', color='r')
    ax2.ticklabel_format(style='sci')

    for tl in ax2.get_yticklabels():
        tl.set_color('r')

    if fname == "":
        plt.show()
    else:
        plt.savefig(fname)


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
                    hh += -pxy*np.log(pxy)
                    if np.isnan(mi) and np.isnan(hh):
                        print "NaN", mi
    return (mi, mi/hh*1e5)


def mi_test():
    img = Image.open("/Prime/Thesis-Data/SimilarityMetric/UNC_train_Case05_T1_sagittal_slice_0265.png").convert("L")
    arr = 1.0*np.array(img)
    (nx, ny) = arr.shape

    t_values = []
    mi_values = []
    nmi_values = []

    j = 0
    r1 = arr[:, (nx/2-j):].flatten()
    r2 = arr[:, 0:(nx/2+j)].flatten()
    #histo2d(r1,r2)
    cm = plt.cm.get_cmap('RdYlBu_r')

    jvalues = range(1, nx/2, 8)
    nBins = 64
    for j in jvalues:
        print j
        r1 = arr[:, (nx/2-j):].flatten()
        r2 = arr[:, 0:(nx/2+j)].flatten()
        (H, xedges, yedges) = np.histogram2d(r1, r2, bins=nBins)
        mi, nmi = compute_MI(H)
        print mi, nmi
        t_values.append(j)
        mi_values.append(mi)
        nmi_values.append(nmi)
        if j == 1:
            np.savetxt('/Prime/Thesis-Data/SimilarityMetric/mitest1.txt', np.array((r1,r2)))
            plt.figure()
            n, xedges, yedges, patches = plt.hist2d(r1, r2, bins=nBins, norm=LogNorm())
            plt.colorbar()
            plt.savefig("/Prime/Thesis-Data/SimilarityMetric/mitest1.pdf")
        if nx/2-j<8:
            np.savetxt('/Prime/Thesis-Data/SimilarityMetric/mitest2.txt', np.array((r1,r2)))
            plt.figure()
            plt.hist2d(r1, r2, bins=nBins, norm=LogNorm())
            plt.colorbar()
            plt.savefig("/Prime/Thesis-Data/SimilarityMetric/mitest2.pdf")


    #jvalues.extend(list(range(nx/2,nx,8)))

    jvalues.extend(list(range(nx/2,nx-8,8)))
    mi_values.extend(mi_values[-2::-1])
    nmi_values.extend(nmi_values[-2::-1])

    print len(jvalues), len(mi_values), len(nmi_values)

    two_scale(jvalues, mi_values, nmi_values, "/Prime/Thesis-Data/SimilarityMetric/miplot.pdf")
    #
    # plt.figure()
    # plt.plot(mi_values, 'g')
    # plt.plot(nmi_values, 'b')
    # plt.ylabel('MI')
    # plt.xlabel('Displacement')
    # plt.savefig("/Prime/Thesis-Data/SimilarityMetric/miplot.pdf")


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