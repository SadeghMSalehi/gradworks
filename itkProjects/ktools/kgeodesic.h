//
//  kgeodesic.h
//  ktools
//
//  Created by Joowhi Lee on 8/20/15.
//
//

#ifndef __ktools__kgeodesic__
#define __ktools__kgeodesic__

#include <stdio.h>

#include "vGraph.h"
#include <vtkPolygon.h>

#include "piOptions.h"


void processGeodesicOptions(pi::Options& opts);
void processGeodesicCommands(pi::Options &opts, pi::StringVector &args);

#endif /* defined(__ktools__kgeodesic__) */
