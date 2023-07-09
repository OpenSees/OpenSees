#ifndef CONTAINER_H
#define CONTAINER_H

#include <OPS_Globals.h>

template <class T> class LIST_CONTAINER
{
// friend class LIST<T>;
// private:
  public:
  LIST_CONTAINER()
	// Empty initialization. Used only for sentinel 
	// nodes in circular lists.
    { 
      previous = this; 
      next = this;
      data = 0;
      sort_field = 0;
    }
  LIST_CONTAINER(T *_data, float *_sort_field = 0)
    { 
      previous = 0;
      next = 0;
      data = _data; 
      sort_field = _sort_field;
    }
  ~LIST_CONTAINER()
    {
      Unlink();
      delete data;
    }
    
  void LinkBefore(LIST_CONTAINER<T> *x)
	// Object places itself in circular linked list
	// before object x.
    {
      next = x;
      previous = x->previous;
      x->previous = this;
      previous->next = this;
    }
  void LinkAfter(LIST_CONTAINER<T> *t)
	// Object places itself in circular linked list
	// after object x.
    {
      previous = t;
      next = t->next;
      t->next = this;
      next->previous = this;
    }
  void Unlink()
	// Object removes itself from doubly linked list
    {
      previous->next = next;
      next->previous = previous;
    }
  LIST_CONTAINER *previous;
  LIST_CONTAINER *next;
  float *sort_field;
  T *data;
};


#define FOR_EACH(_E,_L) for(_E = (_L).SetToHead();(_E!=0);_E = (_L).Forward())
#define FOR_EACH_REVERSE(_E,_L) for(_E = (_L).SetToTail();(_E!=0);_E = (_L).Backward())
#define FOR_EACH_PAIR(_E1,_E2,_L,_S) \
  _E1 = _L.SetToHead(); \
  _E2 = _L.GetNext(); \
  for (; _E2 != 0; _E1 = _L.Forward(), _E2 = _L.GetNext()) \
{_S} \
_E2 = _L.SetToHead(); \
{_S}


#define ORDER_INCREASING 0
#define ORDER_DECREASING 1
#define NOT_SORTED 2
				       
// class: LIST
// description: Doubly linked list with sentinels
template <class T> class LIST
{
public:
  LIST(int _order = NOT_SORTED)
    {
      sentinel = new LIST_CONTAINER<T>;
      current = sentinel;
      order = _order;
      numEntries = 0;
    }
  ~LIST()
	// Remove all list_containers from the list.
    {
      Clear();
      delete sentinel;
    }
    
  int Count()
    {
      return numEntries;
    }
  void Clear()
    // Removes every element from the list. Note that these elements are not
    // deleted, simply removed from the list. Only the list containers that
    // contain these elements are deleted.
    {
      LIST_CONTAINER<T> *ptr, *temp_ptr;
      for (ptr = sentinel->next; ptr != sentinel;)
	{
	  temp_ptr = ptr->next;
	  delete ptr;
	  ptr = temp_ptr;
	}
      current = sentinel;
      numEntries = 0;
    }
  void KILL()
    // Goes through every element of the list and deletes it.
    {
      LIST_CONTAINER<T> *ptr, *temp_ptr;
      for (ptr = sentinel->next; ptr != sentinel;)
	{
	  temp_ptr = ptr->next;
	  delete ptr->data;
	  delete ptr;
	  ptr = temp_ptr;
	}
      sentinel->next = sentinel;
      sentinel->previous = sentinel;
      current = sentinel;
      numEntries = 0;
    }
  void AddHead(T &_Element)
    // Add _Element to the head of the list.
    {
      LIST_CONTAINER<T> *newContainer = new LIST_CONTAINER<T>(&_Element);
      newContainer->LinkAfter(sentinel);
      numEntries++;
    }
  void AddTail(T &_Element)
    // Add _Element to the tail of the list.
    {
      LIST_CONTAINER<T> *newContainer = new LIST_CONTAINER<T>(&_Element);
      newContainer->LinkBefore(sentinel);
      numEntries++;
    }
  T *GetHead()
    // Returns the element at the head of the list.
    {
      if (sentinel->next == sentinel)
	return 0;
      return sentinel->next->data;
    }
  T *GetTail()
    // Returns the element at the tail of the list.
    {
      if (sentinel->previous == sentinel)
	return 0;
      return sentinel->previous->data;
    }
  T *RemoveHead()
    // Removes the head of the list from the list.
    {
      if (sentinel->next == sentinel)
	{
	  opserr << "List is empty, cannot remove head" << endln;
	  return 0;
	}
      if (current == sentinel->next)
	{
	  current = sentinel->next->next;
	}
      T *data = sentinel->next->data;
      delete sentinel->next;
      numEntries--;
      return data;
    }
  T *RemoveTail()
    // Removes the tail of the list from the list.
    {
      if (sentinel->previous == sentinel)
	{
	  opserr << "List is empty, cannot remove tail" << endln;
	  return 0;
	}
      if (current == sentinel->previous)
	{
	  current = sentinel->previous->previous;
	}
      T *data = sentinel->previous->data;
      delete sentinel->previous;
      numEntries--;
      return data;
    }
  int AddSorted(T &_Element, float *_sort_field)
    // Adds an element in sorted order into the list. The order of sorting
    // is determined at the time of creation of the list. _sort_field is
    // a pointer into the field that will determine the position of the 
    // element in the sorted list. (A Pointer to the sort field is necessary
    // as opposed to just a value for the purpose of being able to resort
    // the list when the values change.)
    {
      if (order == NOT_SORTED)
	{
	  opserr << "Calling AddSorted on an unsorted list" << endln;
	  return -1;
	}
      LIST_CONTAINER<T> *ptr;
      LIST_CONTAINER<T> *newContainer = 
	new LIST_CONTAINER<T>(&_Element,_sort_field);
      
      numEntries++;

      if (order == ORDER_INCREASING)
	{
	  for (ptr = sentinel->next; ptr != sentinel; ptr = ptr->next)
	    {
	      if (*_sort_field <= *(ptr->sort_field))
		{
		  newContainer->LinkBefore(ptr);
		  return 0;
		}
	    }
	  newContainer->LinkBefore(ptr);
	  return 0;	  
	} else {
	  for (ptr = sentinel->next; ptr != sentinel; ptr = ptr->next)
	    {
	      if (*_sort_field >= *(ptr->sort_field))
		{
		  newContainer->LinkBefore(ptr);
		  return 0;
		}
	    }
	  newContainer->LinkBefore(ptr);
	  return 0;	  
	}
    }
  void Resort()
    // Resorts the list.
    {
      if (order == NOT_SORTED)
	{
	  opserr << "Calling Resort on an unsorted list" << endln;
	  return;
	}
      LIST_CONTAINER<T> *ptr1, *ptr2, *oldNext;

				// Insertion sort
      for (ptr1 = sentinel->next; ptr1 != sentinel;)
	{

	  
	  for (ptr2 = ptr1->previous; ptr2 != sentinel; ptr2 = ptr2->previous)
	    {
				// The following is the standard C (x ? y : z) construction
				// If order is ORDER_INCREASING, then we check for >=
				// else we check for <=
	      if ((order == ORDER_INCREASING) 
		  ? 
		  (*(ptr1->sort_field) >= *(ptr2->sort_field)) 
		  :
		  (*(ptr1->sort_field) <= *(ptr2->sort_field)) )
		{
		  oldNext = ptr1->next;
		  ptr1->Unlink();
		  ptr1->LinkAfter(ptr2);
		  ptr1 = oldNext;
		  break;
		}
	    }
	  if (ptr2 == sentinel)
	    {
	      oldNext = ptr1->next;
	      ptr1->Unlink();
	      ptr1->LinkAfter(ptr2);
	      ptr1 = oldNext;
	    }
	}
		  
      current = sentinel;
    }
	  
  void Push(T *data)
    // Treats the list as a stack and adds "data" as the first element
    // of the list.
    {
      AddHead(*data);
    }
  T *Pop()
    // Treats the list as a stack and removes and returns the first element
    // of the list.
    {
      if (numEntries == 0)
	{
	  opserr << "ERROR: Tried to Pop an empty list." << endln;
	  return 0;
	}
      if (current == sentinel->next)
	{
	  current = sentinel->next->next;
	}
      numEntries--;
      T *data = sentinel->next->data;
      delete sentinel->next;
      return data;
    }
  T *Peek()
    // Treats the list as a stack and returns the first element of the list.
    {
      if (numEntries == 0)
	return 0;
      else 
	return sentinel->next->data;
    }
  
  // The following are iterative functions on the list. They serve to iterate
  // through the list and to perform operations on the list members.
  T *SetToHead()
    // Sets the "current" pointer to the beginning of the list and returns
    // the data element associated with the beginning of the list.
    {
      current = sentinel->next;
      if (current == sentinel)
	return 0;
      return current->data;
    }
  T *SetToTail()
    // Sets the "current" pointer to the end of the list and returns
    // the data element associated with the end of the list.
    {
      current = sentinel->previous;
      if (current == sentinel)
	return 0;
      return current->data;
    }
  T *Forward()
    // Advances the "current" pointer by one element through the list and
    // returns the data associated with this new element. 
    {
      current = current->next;
      if (current==sentinel)
	return 0;
      return current->data;
    }
  T *Backward()
    // Similar to "Forward()" but in the opposite direction.
    {
      current = current->previous;
      if (current==sentinel)
	return 0;
      return current->data;
    }
  T *GetCurrent()
    // Returns the data element associated with the container to which 
    // "current" is pointing.
    {
      if (current == sentinel)
	return 0;
      return current->data;
    }
  T *GetNext()
    // Returns the data element associated with the container that is one
    // element closer to the tail of the list than the container to which
    // "current" is pointing.
    {
      if (current->next == sentinel)
	return 0;
      return current->next->data;
    }
  T *RemoveCurrent()
    // Removes the element to which "current" is pointing from the list.
    {
      if (current == sentinel)
	{
	  opserr << "Cannot remove sentinel" << endln;
	  return 0;
	}
      LIST_CONTAINER<T> *newCurrent = current->next;
      T *data = current->data;
      prevcurrent = current->previous;
      delete current;
      current = prevcurrent;
      numEntries--;
      return data;
    }
    void SetMark() 
    {
       currentmark = current;
    }
    void ReSetMark()
    {
       current = currentmark;
    }

private:
  LIST_CONTAINER<T> *sentinel;
  LIST_CONTAINER<T> *current;
  LIST_CONTAINER<T> *prevcurrent;
  LIST_CONTAINER<T> *currentmark;

  int order;
  int numEntries;
};

#endif

