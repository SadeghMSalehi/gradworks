__author__ = 'joohwi'

import random,csv,glob

root = "/Prime/Thesis-Data/Rodent-Thickness/RPV_Thickness/RPV1+3/RPV3_AIE_vs_RPV3_Control_Right/Statistical/Permtest"

files = glob.glob("/Prime/Thesis-Data/Rodent-Thickness/RPV_Thickness/RPV1+3/RPV3_AIE_vs_RPV3_Control_Right/Statistical/Permtest/*_regions.txt")

for f in files:
    with open(f, 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')
        data = []
        total_area = 0
        nregions = 0
        for row in spamreader:
            area = float(row[0])
            minp = float(row[1])
            maxp = float(row[2])
            avgp = float(row[3])
            total_area += area
            nregions += 1
            data.append((area, minp, maxp, avgp))
        sumnp = 0.0
        sumxp = 0.0
        sumgp = 0.0
        for ar,np,xp,gp in data:
            sumnp += np/ar
            sumxp += xp/ar
            sumgp += gp/ar

        print sumnp,sumxp,sumgp,nregions,total_area

