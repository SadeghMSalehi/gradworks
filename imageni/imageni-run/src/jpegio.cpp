#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jpeglib.h"
#include "jpegio.h"

#include "Image.h"

void init_jpeg_imageinfo(jpeg_imageinfo_t &jinfo) {
    memset(&jinfo, 0, sizeof(jinfo));
}

void destroy_jpeg_imageinfo(jpeg_imageinfo_t &jinfo) {
    if (jinfo.imageptr) {
        free(jinfo.imageptr);
    }
}

int load2raw(const char* fname, jpeg_imageinfo_t &jinfo) {
    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_error_mgr jerr;
    jpeg_decompress_struct cinfo;
    jpeg_create_decompress(&cinfo);

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(fname, "rb");

    if (!infile) {
        printf("Error opening jpeg file %s\n!", fname);
        return 0;
    }

    unsigned long location = 0;
    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);

    /* setup decompression process and source, then read JPEG header */

    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);

    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Uncomment the following to output image information, if needed. */

    jinfo.width = cinfo.image_width;
    jinfo.height = cinfo.image_height;
    jinfo.colors = cinfo.num_components;

    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    size_t szBuffer = cinfo.output_width * cinfo.output_height * cinfo.num_components;
    /* allocate memory to hold the uncompressed image */
    unsigned char* imagePtr = (unsigned char*) malloc(szBuffer);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        int nStride = cinfo.image_width * cinfo.num_components;
        memcpy(imagePtr + location, row_pointer[0], nStride);
        location += nStride;
    }
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);

    free(row_pointer[0]);
    fclose(infile);
    /* yup, we succeeded! */
    jinfo.imageptr = imagePtr;

    printf("success in loading %s\n", fname);
    return 1;
}

bool raw2jpeg(jpeg_imageinfo_t& jinfo, unsigned char **outbuffer, int *outlen) {
    struct jpeg_compress_struct cinfo = { 0 };
    struct jpeg_error_mgr jerr;
    JSAMPROW row_ptr[1];
    int row_stride;

    *outbuffer = NULL;
    *outlen = 0;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    unsigned long int outbufferLen = 0;
    jpeg_mem_dest(&cinfo, outbuffer, &outbufferLen);

    cinfo.image_width = jinfo.width;
    cinfo.image_height = jinfo.height;
    cinfo.input_components = jinfo.colors;
    if (jinfo.colors == 3) {
        cinfo.in_color_space = JCS_RGB;
    } else if (jinfo.colors == 1) {
        cinfo.in_color_space = JCS_GRAYSCALE;
    }

    jpeg_set_defaults(&cinfo);
    jpeg_start_compress(&cinfo, TRUE);
    row_stride = jinfo.width * cinfo.input_components;

    while (cinfo.next_scanline < cinfo.image_height) {
        row_ptr[0] = &jinfo.imageptr[cinfo.next_scanline * row_stride];
        jpeg_write_scanlines(&cinfo, row_ptr, 1);
    }
    jpeg_finish_compress(&cinfo);

    *outlen = outbufferLen;
    jpeg_destroy_compress(&cinfo);

    return true;
}

bool rgb2gray(jpeg_imageinfo_t& jinfo, jpeg_imageinfo_t& outinfo) {
    outinfo.width = jinfo.width;
    outinfo.height = jinfo.height;
    outinfo.colors = 1;
    outinfo.color_space = JCS_GRAYSCALE;
    outinfo.imageptr = (unsigned char*) malloc(jinfo.width * jinfo.height);
    for (int i = 0; i < jinfo.width * jinfo.height; i++) {
        if (jinfo.colors == 3) {
            int r = jinfo.imageptr[3 * i];
            int g = jinfo.imageptr[3 * i + 1];
            int b = jinfo.imageptr[3 * i + 2];
            int y = (int) ((r + g + b) / 3.0 + .5);
            outinfo.imageptr[i] = y;
        } else {
            outinfo.imageptr[i] = jinfo.imageptr[i];
        }
    }
    return true;
}

void uchar2int(unsigned char* in, int* out, int sz) {
    for (int i = 0; i < sz; i++) {
        out[i] = in[i];
    }
}

void int2uchar(int* in, unsigned char* out, int sz) {
    for (int i = 0; i < sz; i++) {
        out[i] = (unsigned char) in[i];
    }
}

void save2file(const char* fname, jpeg_imageinfo_t& img) {
    unsigned char* jpegBuff;
    int jpegBuffLen;
    raw2jpeg(img, &jpegBuff, &jpegBuffLen);
    printf("saving %s ...\n", fname);
    FILE *fout = fopen(fname, "wb");
    if (fout) {
        int nwrite = jpegBuffLen;
        while (nwrite > 0) {
            nwrite -= fwrite(jpegBuff, 1, nwrite, fout);
        }
        fclose(fout);
    }
}

bool loadImage(const char* fname, IntImage& imgIn) {
    jpeg_imageinfo_t rgb, gray;
    init_jpeg_imageinfo(rgb);
    if (load2raw(fname, rgb) == 0) {
        return false;
    }
    rgb2gray(rgb, gray);
    imgIn.createCopyOf<unsigned char>(gray.imageptr, gray.width, gray.height);
    destroy_jpeg_imageinfo(rgb);
    destroy_jpeg_imageinfo(gray);
    return true;
}

void saveImage(const char* fname, IntImage& imgIn) {
    CImage<unsigned char> tmp;
    tmp.createCopyOf(imgIn);
    jpeg_imageinfo_t info;
    info.width = tmp.getWidth();
    info.height = tmp.getHeight();
    info.imageptr = tmp.getData();
    info.colors = 1;
    info.color_space = JCS_GRAYSCALE;
    save2file(fname, info);
}
