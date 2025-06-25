// The MIT License (MIT)
//
// Copyright (c) 2019 Luigi Pertoldi
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// ============================================================================
//
// Very simple progress bar for c++ loops with internal running variable
//
// Author: Luigi Pertoldi
// Created: 3 Dec 2016
//
#ifndef __PROGRESSBAR_HPP
#define __PROGRESSBAR_HPP

#include <iostream>
#include <ostream>
#include <string>
#include <stdexcept>
#include <codecvt>
#include <locale>

class ProgressBar {

public:
  // default destructor
  ~ProgressBar() = default;

  // delete everything else
  ProgressBar(ProgressBar const &)            = delete;
  ProgressBar &operator=(ProgressBar const &) = delete;
  ProgressBar(ProgressBar &&)                 = delete;
  ProgressBar &operator=(ProgressBar &&)      = delete;

  // default constructor, must call set_niter later
  inline ProgressBar();
  inline ProgressBar(int n, bool showbar = true, std::wostream &out = std::wcerr);

  // reset bar to use it again
  inline void reset();

  // set number of loop iterations
  inline void set_niter(int iter);

  // chose your style
  inline void
  set_done_char(const std::wstring &sym)
  {
    done_char = sym;
  }

  inline void
  set_todo_char(const std::wstring &sym)
  {
    todo_char = sym;
  }

  inline void
  set_opening_char(const std::wstring &sym)
  {
    opening_char = sym;
  }

  inline void
  set_closing_char(const std::wstring &sym)
  {
    closing_char = sym;
  }

  // to show only the percentage
  inline void
  show_bar(bool flag = true)
  {
    do_show_bar = flag;
  }

  // set the output stream
  inline void
  set_output_stream(const std::wostream &stream)
  {
    output.rdbuf(stream.rdbuf());
  }

  // main update function
  inline int update(std::string);

private:
  int msg_width = 0;
  int bar_width = 50;

  int progress;
  int n_cycles;
  int last_perc;
  bool do_show_bar;
  bool update_is_called;

  std::wstring done_char;
  std::wstring todo_char;
  std::wstring opening_char;
  std::wstring closing_char;

  std::wostream &output;
};

inline ProgressBar::ProgressBar()
    : progress(0), n_cycles(0), last_perc(0), do_show_bar(true),
      update_is_called(false), done_char(L"#"), todo_char(L" "),
      opening_char(L"["), closing_char(L"]"), output(std::wcerr)
{
}

inline ProgressBar::ProgressBar(int n, bool showbar, std::wostream &out)
    : progress(0), n_cycles(n), last_perc(0), do_show_bar(showbar),
      update_is_called(false), done_char(L"#"), todo_char(L" "),
      opening_char(L"["), closing_char(L"]"), output(out)
{
}

inline void
ProgressBar::reset()
{
  progress = 0, update_is_called = false;
  last_perc = 0;
  return;
}

inline void
ProgressBar::set_niter(int niter)
{
  if (niter <= 0)
    throw std::invalid_argument(
        "ProgressBar::set_niter: number of iterations null or negative");
  n_cycles = niter;
  return;
}

inline int
ProgressBar::update(std::string message_ = "")
{
  std::wstring message = std::wstring_convert<std::codecvt_utf8<wchar_t>>().from_bytes(message_);

  if (n_cycles == 0) {
    std::cerr <<  "ProgressBar::update: number of cycles not set";
    return -1;
  }

  for (int _ = 0; _ < msg_width; _++)
    output << '\b';

  if (!update_is_called) {
    if (do_show_bar == true) {
      output << opening_char;
      for (int _ = 0; _ < bar_width; _++)
        output << todo_char;
      output << closing_char << " 0%";
    } else
      output << "0%";
  }
  update_is_called = true;

  int perc = 0;

  // compute percentage, if did not change, do nothing and return
  perc = progress * 100. / (n_cycles - 1);
  if (perc < last_perc)
    return 1;

  // update percentage each unit
  if (perc == last_perc + 1) {
    // erase the correct  number of characters
    if (perc <= 10)
      output << "\b\b" << perc << '%';
    else if (perc > 10 && perc < 100)
      output << "\b\b\b" << perc << '%';
    else if (perc == 100)
      output << "\b\b\b" << perc << '%';
  }
  if (do_show_bar == true) {
    // update bar every ten units
    if (perc % 2 == 0) {
      // erase closing bracket
      output << std::wstring(closing_char.size(), '\b');
      // erase trailing percentage characters
      if (perc < 10)
        output << "\b\b\b";
      else if (perc >= 10 && perc < 100)
        output << "\b\b\b\b";
      else if (perc == 100)
        output << "\b\b\b\b\b";

      // erase 'todo_char'
      for (int j = 0; j < bar_width - (perc - 1) / 2; ++j) {
        output << std::wstring(todo_char.size(), '\b');
      }

      // add one additional 'done_char'
      if (perc == 0)
        output << todo_char;
      else
        output << done_char;

      // refill with 'todo_char'
      for (int j = 0; j < bar_width - (perc - 1) / 2 - 1; ++j)
        output << todo_char;

      // readd trailing percentage characters
      output << closing_char << ' ' << perc << '%';
      
    }
  }
  last_perc = perc;
  ++progress;

  msg_width = message.length();

  if (msg_width > 0) {
    output << L" -- " << message;
    // correct for " -- " chars
    msg_width += 4;
  }

  if (perc == 100)
    output << "\n";

  output << std::flush;

  return 1;
}

#endif
