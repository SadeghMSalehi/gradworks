/*
 * jnihelper.h
 *
 *  Created on: Aug 10, 2012
 *      Author: joohwile
 */

#ifndef JNIHELPER_H_
#define JNIHELPER_H_

#include <jni.h>

#define lockArray(var,type) env->Get##type##ArrayElements(j##var, NULL)
#define releaseArray(var,type) env->Release##type##ArrayElements(j##var, p##var, 0)

#define safeFindClassMacro(clsVar, clsName) env->FindClass(clsName); \
	if (clsVar == NULL) {\
		cout << "Can't find class " << #clsName << endl;\
		return NULL;\
	}

#define safeFindMethodMacro(methodVar, clsVar, methodName, methodSig) \
	env->GetMethodID(clsVar, methodName, methodSig); \
	if (methodVar == NULL) { \
		cout << "Can't find the method of " << methodName << methodSig << endl;\
		jmethodID method = env->GetMethodID(clsVar, "getName", "()Ljava/lang/String;");\
		if (method == NULL) { \
			cout << "Can't fetch className" << endl; \
		}\
		jstring name = (jstring) env->CallObjectMethod(clsVar, method);\
		jboolean copy;\
		const char* clsName = env->GetStringUTFChars(name, &copy);\
		cout << "Can't find constructor for " << clsName << "." << methodName\
			<< methodSig << endl;\
		return NULL;\
	}

#define arrayCopyToJNIMacro(var, tname, jname, src) \
	jname* p##var = env->Get##tname##ArrayElements(j##var, NULL); \
	memcpy(p##var, src, sizeof(jname) * env->GetArrayLength(j##var)); \
	env->Release##tname##ArrayElements(j##var, p##var, 0);

#define arrayCopyFromJNIMacro(var, tname, jname, dst) \
	jname* p##var = env->Get##tname##ArrayElements(j##var, NULL); \
	memcpy((void*) dst, p##var, sizeof(jname) * env->GetArrayLength(j##var)); \
	env->Release##tname##ArrayElements(j##var, p##var, 0);

inline jstring jni_translate_string(JNIEnv *jni, const char string[]) {
	// translate c string to java string
	return jni->NewStringUTF(string);
}

inline const char* jni_translate_string(JNIEnv *jni, jstring string) {
	// copy flag
	jboolean copy;

	// translate java string to c string
	return jni->GetStringUTFChars(string, &copy);
}

inline const char* jni_class_name(JNIEnv *jni, jobject object) {
	// get class name
	jclass cls = jni->GetObjectClass(object);
	jmethodID method = jni->GetMethodID(cls, "getName", "()Ljava/lang/String;");
	jstring name = (jstring) jni->CallObjectMethod(cls, method);

	// translate to c string
	return jni_translate_string(jni, name);
}

#endif /* JNIHELPER_H_ */
