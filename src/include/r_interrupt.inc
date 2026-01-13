// Copyright 2018  Marc Hofmann
//
// This file is part of the 'lmSubsets' R extension.
//
// 'lmSubsets' is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// 'lmSubsets' is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with 'lmSubsets'.  If not, see <http://www.gnu.org/licenses/>.



#define MCS_CHECK_INTERRUPT()  r_interrupt_check()



const int r_interrupt_period_ = 8192;



int r_interrupt_tick_;
bool r_interrupt_flag_;



void r_interrupt_setup();
bool r_interrupt_check();
bool r_interrupt_flag();
void r_interrupt_check_aux(void *dummy);



inline void
r_interrupt_setup()
{
  r_interrupt_tick_ = r_interrupt_period_;
  r_interrupt_flag_ = false;
}



// See also:  https://stackoverflow.com/questions/40563522/r-how-to-write-interruptible-c-function-and-recover-partial-results
inline bool
r_interrupt_check()
{
  --r_interrupt_tick_;

  if (r_interrupt_tick_ == 0) {
    r_interrupt_tick_ = r_interrupt_period_;
    r_interrupt_flag_ = (R_ToplevelExec(r_interrupt_check_aux, NULL) == FALSE);
  }

  return r_interrupt_flag_;
}



inline bool
r_interrupt_flag()
{
  return r_interrupt_flag_;
}



inline void
r_interrupt_check_aux(void*)
{
  R_CheckUserInterrupt();
}
