#ifndef _PLSTRING_H
#define _PLSTRING_H

#include <iostream.h>

// A first class string type
class String {
public:
   // Constructors and destructor
   String(const char s[] = ""); // Constructs a deep copy of a C-string.
                                //  ASSUME: s is a valid C-string.
   String(const String& s);     // Constructs a deep copy of s.
   ~String();                   // Deallocates String memory.          

   // Assignment operators
   String& operator= (const String& rhs);   // Assigns a deep copy of rhs. 
   String& operator+= (const String& rhs);  // Adds a deep copy of rhs on the
                                            //   end of this string.

   char& operator[](int i);          // The element at subscript i. 
                                     // ASSUME: i < length of the string
   char operator[](int i) const;     // The element at subscript i.
                                     // ASSUME: i < length of the string

   int length() const;               // Number of string characters.
   const char* charString( ) const;  // C-String equivalent value.

   // Comparison operators
   friend bool operator== (const String& s, const String& t);
   friend bool operator!= (const String& s, const String& t);
   friend bool operator<  (const String& s, const String& t);
   friend bool operator<= (const String& s, const String& t);
   friend bool operator>  (const String& s, const String& t);
   friend bool operator>= (const String& s, const String& t);

   friend ostream& operator<<(ostream& out, const String& s);
               // Writes the C-string equivalent to out.
   friend istream& operator>> (istream& in, String & s);   
               // Reads at most 999 characters up to the next newline from 
               //   from in. The newline is extracted but not assigned.   
   friend String operator+(const String& s, const String& t);  
               // A deep copy of s with a deep copy of t appended to the end.

private:
   char* info_;
};

#endif

