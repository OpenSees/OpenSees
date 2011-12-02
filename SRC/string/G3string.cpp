#include <string.h>
#include <G3string.h>

String::String(const char s[])
{
   info_ = new char[strlen(s) + 1];      // Leave room for the '\0'. 
   strcpy(info_, s);                     // Copy from s to info.
}

String::String(const String& s)
{
   info_ = new char[strlen(s.info_) + 1]; // Allocate memory for the copy.
   strcpy(info_, s.info_);                // Copy the argument's characters.
}

String::~String( )
{
   delete [] info_;    // Deallocate the array.
}

String& String::operator= (const String& rhs) 
{
   if (this != &rhs) {                       // Cover the case of s = s.
      delete [] info_;                       //   Deallocate the old buffer.
      info_ = new char[strlen(rhs.info_) +1];// Allocate memory for a new one.
      strcpy(info_, rhs.info_);              // Copy the characters from the
   }                                         //  right side to the left
   return *this;                             // Return this object.
}     

String& String::operator+= (const String& s)
{
   char* temp = new char[strlen(info_) + strlen(s.info_) + 1]; // Create a new 
                           // array to hold the two string arrays.
   strcpy(temp,info_);     // Copy the characters from this array into temp. 
   strcat(temp,s.info_);   // Then copy the characters of s.info_ into temp. 

   delete [] info_;        // Replace the old value for info_ by temp.
   info_ = temp;
   return *this;
}

String operator+ (const String& s, const String& t)
{
   String temp(s);    // temp is a copy of s.
   temp += t;         // Add t to the end of temp.
   return temp;      
}

char& String::operator[](int i)
{
   return info_[i];
} 

char String::operator[](int i) const
{
   return info_[i];
} 

const char* String::charString() const       
{  
   return info_;
}

ostream& operator<< (ostream& out, const String& s)
{
   out << s.charString();
   return out;
}

istream& operator>> (istream& in, String& s)   
{
   char buffer[1000];             // Buffer to store the stream characters
   in.getline(buffer,1000,'\n');  // Remove up to 999 characters from in, 
                                  //   up to and including first occurrence                                
                                  //   of '\n'.  Store all but '\n' in the  
                                  //   buffer; terminate the buffer with '\0'.   
   s = String(buffer);            // Create a new String from the buffer and
                                  //   assign it to s.
   return in;
}

int String::length() const
{
   return strlen(info_); 
}

bool operator== (const String& s, const String& t)
{
   return strcmp(s.info_,t.info_) == 0;
}

bool operator!= (const String& s, const String& t)
{
   return strcmp(s.info_,t.info_) != 0;
}

bool operator< (const String& s, const String& t)
{
   return strcmp(s.info_,t.info_) < 0;
}

bool operator<= (const String& s, const String& t)
{
   return strcmp(s.info_,t.info_) <= 0;
}

bool operator> (const String& s, const String& t)
{
   return strcmp(s.info_,t.info_) > 0;
}

bool operator>= (const String& s, const String& t)
{
   return strcmp(s.info_,t.info_) >= 0;
}


