/* Mike Gertz - 2-Aug-98 */

#ifndef SNOPTFILEWRAPPER
#define SNOPTFILEWRAPPER

#ifdef __cplusplus
extern "C" {
#endif

void
snoptopenappend_(integer *iunit, char *name,
		 integer *inform, ftnlen name_len);

void
snoptfilewrapper_(char *name__, integer *ispec, integer *inform__,
		  char *cw, integer *lencw, integer *iw,
		  integer *leniw, doublereal *rw, integer *lenrw,
		  ftnlen name_len, ftnlen cw_len);

void snoptclose_(integer *iunit);

void snoptopenread_(integer *iunit, char *name, integer *inform,
		    ftnlen name_len);


#ifdef __cplusplus
}
#endif

#endif
