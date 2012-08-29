/*
 * jpegio.h
 *
 *  Created on: Jun 27, 2012
 *      Author: joohwile
 */

#ifndef JPEGIO_H_
#define JPEGIO_H_

#include "jpeglib.h"
#include "Image.h"

typedef struct jpeg_imageinfo {
	int width;
	int height;
	int colors;
	int color_space;
	unsigned char* imageptr;
} jpeg_imageinfo_t;

void init_jpeg_imageinfo(jpeg_imageinfo_t &jinfo);
void destroy_jpeg_imageinfo(jpeg_imageinfo_t &jinfo);

int load2raw(const char* fname, jpeg_imageinfo_t &jinfo);
bool raw2jpeg(jpeg_imageinfo_t& jinfo, unsigned char **outbuffer, int *outlen);
bool rgb2gray(jpeg_imageinfo_t& jinfo, jpeg_imageinfo_t& out);
void save2file(const char* fname, jpeg_imageinfo_t &img);
bool loadImage(const char* fname, IntImage& imgIn);

template <typename T>
void writeImage(const char* fname, CImage<T>& img) {
	CImage<unsigned char> grayImg;
	grayImg.createCopyOf(img);

	jpeg_imageinfo_t imgInfo;
	init_jpeg_imageinfo(imgInfo);

	imgInfo.width = img.getWidth();
	imgInfo.height = img.getHeight();
	imgInfo.colors = 1;
	imgInfo.color_space = JCS_GRAYSCALE;
	imgInfo.imageptr = grayImg.getData();
	save2file(fname, imgInfo);
}

#endif /* JPEGIO_H_ */
