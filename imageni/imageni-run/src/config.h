/*
 * android-extension.h
 *
 *  Created on: Jul 19, 2012
 *      Author: joohwile
 */

#ifndef CONFIG_H
#define CONFIG_H

//#define ANDROID

#ifdef ANDROID
#include <android/log.h>
#include <android/bitmap.h>
#define LOGI(...) ((void)__android_log_print(ANDROID_LOG_INFO, "FACEAR", __VA_ARGS__))
#define LOGD(...) ((void)__android_log_print(ANDROID_LOG_DEBUG, "FACEAR", __VA_ARGS__))
#else
#ifdef DEBUG
#define LOGI(...) (printf(__VA_ARGS__))
#define LOGD(...) (printf(__VA_ARGS__))
#else
#define LOGI(...)
#define LOGD(...)
#endif
#endif

#ifdef WIN32
//#include <Windows.h>
#endif

#endif /* ANDROID_EXTENSION_H_ */
