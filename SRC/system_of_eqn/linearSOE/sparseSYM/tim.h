/*
 * File:  tim.h
 * ============
 *
 */

#ifndef TIM_H
#define TIM_H

/* Timer Codes */

#define SYMBOLIC 0
#define FACTOR 1
#define SOLVE 2
#define END 3


/*
 * Function:  timer(tcode)
 * =======================
 * This routine keeps track of the CPU statistics.
 *
 * INPUT
 * -----
 * tcode => timer codes : name of processes to be timed, including 'begin',
 *			  'end', 'back', which are the control input 
 */

void timer(int tcode);


#endif
