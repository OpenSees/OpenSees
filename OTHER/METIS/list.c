/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * list.c
 *
 * This file contains routines that manipulate link lists
 *
 * 8/31/94
 * George
 *
 * $Id: list.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"


/*************************************************************************
* This function adds a node stored in ListNodeType variable to a list. 
* It returns the list itself.
**************************************************************************/
ListNodeType *addthisnode(ListNodeType *hptr, ListNodeType *nptr)
{
  nptr->link = hptr;

  return nptr;
}


/*************************************************************************
* This function deletes a node from an adjacency list
**************************************************************************/
ListNodeType *delthisnode(ListNodeType *adj, int node, ListNodeType **retptr)
{
  ListNodeType *ptr;
  ListNodeType *prevptr;

  if (adj == NULL) 
    errexit("\nTrying to delete from an empty list!");

  if (adj->link == NULL) {
    if (adj->pos == node) {
      *retptr = adj;
      return NULL;
    }
    else 
      errexit("\nElement not found in list!");
  }

  if (adj->pos == node) {
    ptr = adj;
    adj = adj->link;
    *retptr = ptr;
    return adj;
  }

  prevptr = adj;
  for (ptr = adj->link; ptr != NULL; ptr = ptr->link) {
    if (ptr->pos == node)
      break;
    prevptr = ptr;
  }
  if (ptr != NULL) {
    prevptr->link = ptr->link;
    *retptr = ptr;
    return adj;
  }
  else 
    errexit("\nElement not found in list!");
}





