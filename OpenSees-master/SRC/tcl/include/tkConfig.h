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
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2001-08-18 00:48:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/include/tkConfig.h,v $
                                                                        
                                                                        
/*
 * tkConfig.h --
 *
 *	This file is included by all of the Tk C files.  It contains
 *	information that may be configuration-dependent, such as
 *	#includes for system include files and a few other things.
 *
 * Copyright 1991 Regents of the University of California
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any purpose and without
 * fee is hereby granted, provided that this copyright
 * notice appears in all copies.  The University of California
 * makes no representations about the suitability of this
 * software for any purpose.  It is provided "as is" without
 * express or implied warranty.
 *
 * $Header: /usr/local/cvs/OpenSees/SRC/tcl/include/tkConfig.h,v 1.1 2001-08-18 00:48:55 fmk Exp $ SPRITE (Berkeley)
 */

#ifndef _TKCONFIG
#define _TKCONFIG

/*
 * Macro to use instead of "void" for arguments that must have
 * type "void *" in ANSI C;  maps them to type "char *" in
 * non-ANSI systems.  This macro may be used in some of the include
 * files below, which is why it is defined here.
 */

#ifndef VOID
#   ifdef __STDC__
#       define VOID void
#   else
#       define VOID char
#   endif
#endif

#include <stdio.h>
#include <ctype.h>
#include <fcntl.h>
#include <math.h>
#include <pwd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/time.h>
#ifndef _TCL
#   include <tcl.h>
#endif
#include <X11/Xlib.h>
#include <X11/cursorfont.h>
#include <X11/keysym.h>
#include <X11/Xatom.h>
#include <X11/Xproto.h>
#include <X11/Xresource.h>
#include <X11/Xutil.h>

/*
 * Macro to use instead of "void" for arguments that must have
 * type "void *" in ANSI C;  maps them to type "char *" in
 * non-ANSI systems.
 */

#ifndef VOID
#   ifdef __STDC__
#       define VOID void
#   else
#       define VOID char
#   endif
#endif

/*
 * Not all systems declare the errno variable in errno.h. so this
 * file does it explicitly.
 */

extern int errno;

/*
 * Define OPEN_MAX if it isn't already defined for this system.
 */

#ifndef OPEN_MAX
#   define OPEN_MAX 256
#endif

/*
 * The following macro defines the type of the mask arguments to
 * select:
 */

#if (defined(sun) && !defined(sprite)) || defined(linux)
#   define SELECT_MASK fd_set
#else
#   if defined(_IBMR2)
#	define SELECT_MASK void
#   else
#	define SELECT_MASK int
#   endif
#endif

/*
 * Declarations for various library procedures that aren't declared
 * in a header file.
 */

extern int		close _ANSI_ARGS_((int fd));
//extern int		gettimeofday _ANSI_ARGS_((struct timeval *tp,
//			    struct timezone *tzp));
extern uid_t		getuid _ANSI_ARGS_((void));
#if !(defined(_CRAY) || defined(sparc) || defined(_IBMR2))
extern int		open _ANSI_ARGS_((CONST char *path, int flags, ...));
#endif
//extern void		panic _ANSI_ARGS_(VARARGS);
extern int		read _ANSI_ARGS_((int fd, char *buf, int numBytes));
//extern int		select _ANSI_ARGS_((size_t nfds, SELECT_MASK *readfds,
//			    SELECT_MASK *writefds, SELECT_MASK *exceptfds,
//			    const struct timeval *timeout));

#endif /* _TKCONFIG */
