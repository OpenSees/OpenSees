/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2007-04-06 03:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/java/OpenSeesEvaluator.h,v $
// 
// Written: fmk 
// Created: 07/04
//

#include <jni.h>
/* Header for class OpenSeesEvaluator */

#ifndef _Included_OpenSeesEvaluator
#define _Included_OpenSeesEvaluator
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     OpenSeesEvaluator
 * Method:    openSeesInit
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_OpenSeesEvaluator_openSeesInit
  (JNIEnv *, jobject);

/*
 * Class:     OpenSeesEvaluator
 * Method:    openSeesEval
 * Signature: (Ljava/lang/String;I)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_OpenSeesEvaluator_openSeesEval
  (JNIEnv *, jobject, jstring, jint);

/*
 * Class:     OpenSeesEvaluator
 * Method:    openSeesQuit
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_OpenSeesEvaluator_openSeesQuit
  (JNIEnv *, jobject);

#ifdef __cplusplus
}
#endif
#endif
