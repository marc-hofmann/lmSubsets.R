// Copyright 2018  Marc Hofmann
//
// This file is part of the 'mcs' library (see
// <https://github.com/marc-hofmann/mcs/>).
//
// 'mcs' is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// 'mcs' is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with 'mcs'.  If not, see <http://www.gnu.org/licenses/>.



#ifndef MCS_DETAIL_INTERRUPT_HH
#define MCS_DETAIL_INTERRUPT_HH



#ifndef MCS_CHECK_INTERRUPT
#define MCS_INTERRUPT_FLAG()   ((void) 0)
#define MCS_INTERRUPT_BREAK()  ((void) 0)
#else
#define MCS_INTERRUPT_FLAG()   MCS_CHECK_INTERRUPT()
#define MCS_INTERRUPT_BREAK()  if (MCS_INTERRUPT_FLAG())  break
#endif



#endif
