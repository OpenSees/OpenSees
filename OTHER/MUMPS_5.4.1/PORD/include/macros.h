/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: macros.h
/
/ author        J"urgen Schulze, University of Paderborn
/ created       99jan24
/
/ This file contains some useful macros
/
******************************************************************************/


/* Some compilers (VC++ for instance) define a min and a max in the stdlib */
#ifdef min
# undef min
#endif
#ifdef max
# undef max
#endif
#define  min(a,b)  ((a) < (b) ? (a) : (b))
#define  max(a,b)  ((a) > (b) ? (a) : (b))

#define  mymalloc(ptr, nr, type) \
           if (!(ptr = (type*)malloc((max(nr,1)) * sizeof(type)))) \
            { printf("malloc failed on line %d of file %s (nr=%d)\n", \
                     __LINE__, __FILE__, nr); \
              exit(ERR); \
            }

#define myrealloc(ptr, nr, type) \
           if (!(ptr = (type*)realloc(ptr, (nr) * sizeof(type)))) \
            { printf("realloc failed on line %d of file %s (nr=%d)\n", \
                     __LINE__, __FILE__, nr); \
              exit(ERR); \
            }

#define myrandom(range) \
           rand() % (range);

#define swap(a, b, tmp) \
           { (tmp) = (a); (a) = (b); (b) = (tmp); }

#define seed() \
           srand((PORD_INT)time(0) % 10000);

#define bit(var, d) \
           ((var) & (1 << (d)))

#define negbit(var, d) \
           ((var) ^ (1 << (d)))

#define waitkey() \
           { char _s[MAX_LINE_LEN]; printf("\n<RETURN>"); gets(_s); }

#define pord_resettimer(var) \
           var = 0;

#define pord_starttimer(var) \
           var -= ((FLOAT)clock()/CLOCKS_PER_SEC);

#define pord_stoptimer(var) \
           var += ((FLOAT)clock()/CLOCKS_PER_SEC);

#define quit() \
           exit(ERR);

#ifdef PARIX
#undef pord_starttimer(var)
#ifdef __EPX
#define pord_starttimer(var) \
           var -= ((FLOAT)TimeNow()/CLOCK_TICK);
#else
#define pord_starttimer(var) \
           var -= ((FLOAT)TimeNowHigh()/CLK_TCK_HIGH);
#endif
#undef pord_stoptimer(var)
#ifdef __EPX
#define pord_stoptimer(var) \
           var += ((FLOAT)TimeNow()/CLOCK_TICK);
#else
#define pord_stoptimer(var) \
           var += ((FLOAT)TimeNowHigh()/CLK_TCK_HIGH);
#endif
#undef quit()
#define quit() \
           exit(ERR);
#endif
