// Define 32 bit signed and unsigned integers.
// Change these definitions, if necessary, on 64 bit computers


#include <stdlib.h>
#include <localMath.h>

uint32 _lrotl (uint32 x, int r) {
  return (x << r) | (x >> (sizeof(x)*8-r));
}


char *utoa(unsigned value, char *digits, int base)
{
    char *s, *p;

    s = "0123456789abcdefghijklmnopqrstuvwxyz"; /* don't care if s is in
                                                 * read-only memory
                                                 */
    if (base == 0)
        base = 10;
    if (digits == NULL || base < 2 || base > 36)
        return NULL;
    if (value < (unsigned) base) {
        digits[0] = s[value];
        digits[1] = '\0';
    } else {
        for (p = utoa(value / ((unsigned)base), digits, base);
             *p;
             p++);
        utoa( value % ((unsigned)base), p, base);
    }
    return digits;
}   

char *itoa(int value, char *digits, int base)
{
    char *d;
    unsigned u; /* assume unsigned is big enough to hold all the
                 * unsigned values -x could possibly be -- don't
                 * know how well this assumption holds on the
                 * DeathStation 9000, so beware of nasal demons
                 */

    d = digits;
    if (base == 0)
        base = 10;
    if (digits == NULL || base < 2 || base > 36)
        return NULL;
    if (value < 0) {
        *d++ = '-';
        u = -value;
    } else
        u = value;
    utoa(u, d, base);
    return digits;
}


