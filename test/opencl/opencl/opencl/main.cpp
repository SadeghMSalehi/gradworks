//
//  main.cpp
//  opencl
//
//  Created by Joohwi Lee on 4/11/13.
//  Copyright (c) 2013 UNC. All rights reserved.
//

#include <iostream>
#include <OpenCL/OpenCL.h>
#include "squares.cl.h"


static void print_device_info(cl_device_id device) {
    char name[128];
    char vendor[128];

    clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(name), name, NULL);
    clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(name), name, NULL);
    fprintf(stdout, "%s : %s\n", vendor, name);
}

int main(int argc, const char * argv[])
{
    dispatch_queue_t queue = gcl_create_dispatch_queue(CL_DEVICE_TYPE_GPU, NULL);
    if (queue == NULL) {
        std::cout << "OpenCL is not available" << std::endl;
        return 0;
    }
    cl_device_id gpu = gcl_get_device_id_with_dispatch_queue(queue);
    print_device_info(gpu);

    dispatch_release(queue);

    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}

