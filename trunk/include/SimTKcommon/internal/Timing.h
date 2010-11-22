#ifndef SimTK_SimTKCOMMON_TIMING_H_
#define SimTK_SimTKCOMMON_TIMING_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**@file This file ensures that we have access to the Posix time functions
clock_getttime() and nanosleep(), and also provides some convenient methods
for use in common timing situations. **/

/**@defgroup TimingFunctions Timing Functions
 * @ingroup GlobalFunctions
 *
 * These functions provide a convenient way to do timings, either for real
 * time use or performance measurement. Both elapsed timings and CPU timings
 * are supported, with the latter on a per-process or per-thread basis.
 * Times are returned as a double precision floating point number of seconds,
 * which is usually the most convenient form. Alternatives are available that
 * return timings as a 64 bit integer count of nanosecond ticks, providing
 * the highest resolution for very short measurements.
 *
 * Note that you can also use the Posix functions clock_gettime() and 
 * nanosleep() on any SimTK-supported platform.
 *
 * On Linux or Mac systems, use of any of these timing functions may require 
 * linking with the librt realtime library (-lrt).
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/Constants.h"

// This header is needed on Mac and Linux for the Posix time functions
// and the timespec struct.
#include <ctime>

#if defined(_MSC_VER)
    /* On Windows, the timespec struct is not defined. However, note that the 
     * timespec struct is also defined in the pthread.h header on Windows, so
     * the guard symbols must match here to avoid a duplicate.
     */
    #ifndef HAVE_STRUCT_TIMESPEC
    #define HAVE_STRUCT_TIMESPEC 1
    struct timespec {
            long tv_sec;  // TODO: this should be time_t but must fix in pthreads too
            long tv_nsec;
    };
    #endif /* HAVE_STRUCT_TIMESPEC */
#endif

#if defined(_MSC_VER) || defined(__APPLE__)
    /* On Windows and OSX, these Posix time functions are missing.
     */
    typedef long clockid_t;

    /* These constants are the clock ids we support. All the varieties
     * of CLOCK_MONOTONIC are high resolution with no NTP adjustments.
     * I measured the resolutions on a single Windows 7 machine; hopefully
     * they are typical (resolution here means how often they are updated):
     *   - MONOTONIC (counter):    0.001ms      1us
     *   - REALTIME (time of day):    1ms    1000us
     *   - CPUTIME (either):         20ms   20000us
     * These are slightly conservative resolutions so you should be able
     * to achieve them in practice.
     */
    #define CLOCK_REALTIME              1   /* time of day clock, from 1/1/1970 */ 
    #define CLOCK_MONOTONIC             2   /* counter from last boot time */
    #define CLOCK_MONOTONIC_HR          3   /* "high resolution" (same) */
    #define CLOCK_MONOTONIC_RAW         4   /* "not subject to NTP adjustments" (same) */
    #define CLOCK_THREAD_CPUTIME_ID     5   /* lifetime cpu time (kernel+user) of the current thread */ 
    #define CLOCK_PROCESS_CPUTIME_ID    6   /* cumulative cpu time of all threads of this process */

    /* Returns zero if it succeeds (or if tp==NULL); otherwise EINVAL. 
     * On a Linux system, this requires including <time.h> (or <ctime>)
     * and linking with -lrt to get the realtime library.
     */
    SimTK_SimTKCOMMON_EXPORT int clock_gettime(clockid_t clock_id, struct timespec *tp); 

    /* Posix nanosleep() sleeps the indicated number of nanoseconds and returns 0, or 
     * if it is interrupted early it returns how much time was left in rem and returns
     * EINTR. Ours is not interruptable so will always succeed and return rem==0. It is
     * OK if rem is NULL, but req==NULL or req<0 returns EINVAL. A time of req==0 is
     * allowed and our interpretation is that the thread relinquishes its time slice
     * to another ready-to-run thread if there is one, otherwise returns immediately.
     * This implementation rounds the desired sleep time to the nearest millisecond.
     * On a Linux system, this requires including <time.h> (or <ctime>).
     */
    SimTK_SimTKCOMMON_EXPORT int nanosleep(const struct timespec* req, struct timespec* rem);
#endif


namespace SimTK {

/** @defgroup TimeConversions Timespec/Nanosecond/Second Conversions
    @ingroup TimingFunctions

These inline functions provide a fast and convenient way for doing arithmetic
with the ugle Posix timespec struct. Use them to convert the timespec to
a long long integer number of nanoseconds, do arithmetic in that form, and
then convert back. Negative times are handled correctly (they come up as
the result of subtraction and comparisons). 

We usually prefer to deal with times as a double precision floating point
number of seconds and functions are provided for converting between 
nanoseconds and seconds in this format. 

@par Cautions:
    - For long intervals the precision of time in seconds will necessarily be 
      less than the precision of the nanosecond count, since IEEE double 
      precision has a 53 bit mantissa, while 63 bits are available for the count.
    - A signed long long integer containing a count of nanosecond ticks
      is limited to time intervals of about +/- 292 years, which is 
      substantially less than can be contained in a timespec, but is plenty 
      for interval timing. 
**/

/**@{**/
/** Convert a time stored in a timespec struct to the equivalent number
of nanoseconds (as a signed quantity). @see nsToTimespec() **/    
inline long long timespecToNs(const timespec& ts)
{   return (long long)ts.tv_sec*1000000000LL + (long long)ts.tv_nsec; }

/** Given a signed number of nanoseconds, convert that into seconds and 
leftover nanoseconds in a timespec struct. @see timespecToNs() **/
inline void nsToTimespec(const long long& ns, timespec& ts) {
    ts.tv_sec  = (long)(ns / 1000000000LL); // signed
    if (ns >= 0) ts.tv_nsec =  (long)(  ns  % 1000000000LL);
    else         ts.tv_nsec = -(long)((-ns) % 1000000000LL);
}

/** Given a count of nanosecond ticks as a signed 64 bit integer, return
the same time interval as a double precision floating point number of
seconds. See @ref TimeConversions for cautions. @see secToNs() **/
inline double nsToSec(const long long& ns) 
{   return (double)(ns*SimTK_NS_TO_S); }

/** Given a signed time interval as a double precision floating point number of
seconds, return the same time interval as a count of nanosecond ticks in a 
signed 64 bit integer. See @ref TimeConversions for cautions. @see nsToSec() **/
inline long long secToNs(const double& s) 
{   return (long long)(s*SimTK_S_TO_NS); }
/**@}**/

/** @defgroup CPUTimers Measuring CPU Time
    @ingroup TimingFunctions

These functions provide measurement of CPU time consumed by a process or
by individual threads. Time included both kernel and user time together,
and is reported as a double precision floating point value in seconds.
CPU timers typically have a very coarse resolution, likely to be in the 
10-50ms range depending on the particulars of your system. That means you 
won't get repeatable results unless you measure substantial amounts of 
CPU time; don't expect to get meaningful information measuring CPU times
of less than a second or so. **/
/**@{**/

/** Return the cumulative CPU time in seconds (both kernel and user time) that 
has been used so far by any of the threads in the currently executing process.

@return CPU time used since this process was created, in seconds, as
a double precision floating point number.
@see threadCpuTime() **/
inline double cpuTime() {
    timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return (double)(timespecToNs(ts)*SimTK_NS_TO_S);
}

/** Return the total CPU time in seconds (both kernel and user time) that 
has been used so far by the currently executing thread.

@return CPU time used since this thread was created, in seconds, as
a double precision floating point number.
@see cpuTime()
**/
inline double threadCpuTime() {
    timespec ts;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return (double)(timespecToNs(ts)*SimTK_NS_TO_S);
}
/**@}**/

/** @defgroup ElapsedTime High-Resolution Elapsed Time Measurement
    @ingroup TimingFunctions

These functions provide access to the system's high resolution interval
timer, which measures elapsed time from some arbitrary starting point
(typically since the system was last booted). It is expected that this
timer provides very precise measurement of short time intervals, but
cannot be depended upon to measure very long periods without drifting.
That is, it may not be synchronized to the system time of day.

Generally it is most convenient to measure intervals as a floating
point number of seconds, however this provides a variable amount of
precision as the absolute value of the timer increases. For maximum
precision, you can obtain the timer value as an integer number of
nanoseconds instead. You can improve accuracy by subtacting the integer 
counts first before converting to seconds. The actual resolution is 
system-dependent, but it should be able to accurately measure elapsed 
times of 1ms or less, substantially less on some systems. **/
/**@{**/

/** Return current time on the high-resolution interval timer in
nanoseconds, as a 64-bit integer count. Generally it is more convenient
to use realTime() which reports the interval time in
seconds instead, but the nanosecond count is best for maximum
accuracy.

@return Elapsed nanoseconds since some arbitrary time, as a 64 bit 
integer count.

@see realTime(), nsToSec()
**/
inline long long realTimeInNs() {
    timespec ts;
    #ifdef CLOCK_MONOTONIC_RAW
        clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    #else
        clock_gettime(CLOCK_MONOTONIC, &ts);
    #endif
    return timespecToNs(ts);
}

/** Return current time on the high-resolution interval timer in
seconds. For maximum precision, you can improve repeatability and 
accuracy somewhat by obtaining the interval times as integer counts 
using realTimeInNs(). 

@return Elapsed seconds since some arbitrary time, as a double 
precision floating point number.

@see realTimeInNs()
**/ 
inline double realTime() {
    return nsToSec(realTimeInNs());
}
/**@}**/

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_TIMING_H_
