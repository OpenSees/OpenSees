#include <unistd.h>

int check_mem_leak(char *where)
{ 
    void *addr;
    addr = sbrk(0);
    printf("\tsbrk(0) %s: addr = %u\n", where, addr);
}
