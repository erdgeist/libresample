/**********************************************************************

  CPP single class implementation of libresample's high quality mode only
  by Dirk Engling <erdgeist@erdgeist.org>.

  All contributions by erdgeist@erdgeist.org are under beerware license.

  Effectively this means that the original two licenses below fully apply.

 **********************************************************************/
/**********************************************************************

  Real-time library interface by Dominic Mazzoni

  Based on resample-1.7:

  http://www-ccrma.stanford.edu/~jos/resample/

  Dual-licensed as LGPL and BSD; see README.md and LICENSE* files.

  License and warranty:

All of the files in this package are Copyright 2003 by Dominic Mazzoni
dominic@minorninth.com.  This library was based heavily on Resample-1.7,
Copyright 1994-2002 by Julius O. Smith III jos@ccrma.stanford.edu, all rights
reserved.

Permission to use and copy is granted subject to the terms of the "GNU Lesser
General Public License" (LGPL) as published by the Free Software Foundation;
either version 2.1 of the License, or any later version. In addition, Julius O.
Smith III requests that a copy of any modified files be sent by email to
jos@ccrma.stanford.edu so that he may incorporate them into the CCRMA version.

This library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
details.

Permission to use and copy is also granted subject to the terms of the BSD
license found in LICENSE-BSD.txt.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 **********************************************************************/

#pragma once

#include <inttypes.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

class Resampler {

public:
    Resampler() = delete;
    Resampler(double minFactor, double maxFactor) : _minFactor(minFactor), _maxFactor(maxFactor)
    {
        /* Just exit if we get invalid factors */
        if (_minFactor <= 0.0 || _maxFactor <= 0.0 || _maxFactor < _minFactor)
            throw std::invalid_argument("Factors out of range");

        _LpScl = 1.0;

        double Rolloff = 0.90;
        double frq = 0.5 * Rolloff;
        double Beta = 6.0;
        double IBeta = 1.0 / Izero(Beta);

        /* Calculate ideal lowpass filter impulse response coefficients: */
        _Imp[0] = 2.0 * frq;
        for (size_t i = 1; i < _NWING; i++) {
            double temp = M_PI * (double)i / (double)_Npc;
            _Imp[i] = sin(2.0 * temp * frq) / temp; /* Analog sinc function, cutoff = frq */
        }

        /*
         * Calculate and Apply Kaiser window to ideal lowpass filter.
         * Note: last window value is IBeta which is NOT zero.
         * You're supposed to really truncate the window here, not ramp
         * it to zero. This helps reduce the first sidelobe.
         */
        double inm1 = 1.0 / ((double)(_NWING - 1));
        for (size_t i = 1; i < _NWING; i++) {
            double temp = (double)i * inm1;
            double temp1 = 1.0 - temp * temp;
            temp1 = (temp1<0? 0: temp1); /* make sure it's not negative since
                                            we're taking the square root - this
                                            happens on Pentium 4's due to tiny
                                            roundoff errors */
            _Imp[i] *= Izero(Beta*sqrt(temp1)) * IBeta;
        }

        /* Storing deltas in ImpD makes linear interpolation
           of the filter coefficients faster */
        for (size_t i = 0; i < _NWING - 1; i++)
            _ImpD[i] = _Imp[i + 1] - _Imp[i];

        /* Last coeff. not interpolated */
        _ImpD[_NWING - 1] = - _Imp[_NWING - 1];

        /* Calc reach of LP filter wing (plus some creeping room) */
        size_t Xoff_min = 18.0 * std::max(1.0, 1.0 / _minFactor) + 10.0;
        size_t Xoff_max = 18.0 * std::max(1.0, 1.0 / _maxFactor) + 10.0;
        _Xoff = std::max(Xoff_min, Xoff_max);

        /* Make the inBuffer size at least 4096, but larger if necessary
           in order to store the minimum reach of the LP filter and then some.
           Then allocate the buffer an extra Xoff larger so that
           we can zero-pad up to Xoff zeros at the end when we reach the
           end of the input samples. */
        _XSize = std::max<size_t>(2 * _Xoff + 10, 4096);
        _X = new float[_XSize + _Xoff];
        _Xp = _Xoff;
        _Xread = _Xoff;

        /* Need Xoff zeros at begining of X buffer */
        for (size_t i = 0; i < _Xoff; i++)
            _X[i] = 0;

        /* Make the outBuffer long enough to hold the entire processed
           output of one inBuffer */
        _YSize = (size_t)(((double)_XSize) * maxFactor + 2.0);
        _Y = new float[_YSize];
        _Yp = 0;

        _Time = (double)_Xoff; /* Current-time pointer for converter */
    }

    ~Resampler()
    {
        delete[] _X;
        delete[] _Y;
    }

    size_t process(double       factor,
                   const short *inBuffer,
                   size_t       inBufferLen,
                   size_t      &inBufferUsed, /* output param */
                   short       *outBuffer,
                   size_t       outBufferLen)
    {
        /* Initialize inBufferUsed and outSampleCount to 0 */
        inBufferUsed = 0;
        size_t outSampleCount = 0;

        if (factor < _minFactor || factor > _maxFactor)
            throw std::invalid_argument("factor out of range");

        /* Start by copying any samples still in the Y buffer to the output
           buffer */
        if (_Yp && (outBufferLen - outSampleCount) > 0) {
            size_t len = std::min<size_t>(outBufferLen - outSampleCount, _Yp);
            for (size_t i = 0; i < len; i++)
                outBuffer[outSampleCount + i] = (short)_Y[i];
            outSampleCount += len;
            for (size_t i = 0; i < _Yp - len; i++)
                _Y[i] = _Y[i + len];
            _Yp -= len;
        }

        /* If there are still output samples left, return now - we need
           the full output buffer available to us... */
        if (_Yp)
            return outSampleCount;

        /* Account for increased filter gain when using factors less than 1 */
        float LpScl = _LpScl;
        if (factor < 1)
            LpScl = LpScl * factor;

        for (;;) {
            /* Copy as many samples as we can from the input buffer into X */
            size_t len = _XSize - _Xread;

            if (len >= (inBufferLen - inBufferUsed))
                len = (inBufferLen - inBufferUsed);

            for (size_t i = 0; i < len; i++)
                _X[_Xread + i] = (float)inBuffer[inBufferUsed + i];

            inBufferUsed += len;
            _Xread += len;

            if (_Xread <= 2 * _Xoff)
                break;

            /* Resample stuff in input buffer */
            size_t Nout, Nx = _Xread - 2 * _Xoff;
            if (factor >= 1)     /* SrcUp() is faster if we can use it */
                Nout = lrsSrcUp(factor, Nx, LpScl);
            else
                Nout = lrsSrcUD(factor, Nx, LpScl);

            _Time -= Nx;         /* Move converter Nx samples back in time */
            _Xp += Nx;           /* Advance by number of samples processed */

            /* Calc time accumulation in Time */
            size_t Ncreep = (size_t)(_Time) - _Xoff;
            if (Ncreep) {
                _Time -= Ncreep;  /* Remove time accumulation */
                _Xp += Ncreep;    /* and add it to read pointer */
            }

            /* Copy part of input signal that must be re-used */
            size_t Nreuse = _Xread - (_Xp - _Xoff);
            for (size_t i = 0; i < Nreuse; i++)
                _X[i] = _X[i + (_Xp - _Xoff)];

            _Xread = Nreuse;  /* Pos in input buff to read new data into */
            _Xp = _Xoff;

            /* Check to see if output buff overflowed (shouldn't happen!) */
            if (Nout > _YSize)
                throw std::runtime_error("libresample: Output array overflow!");

            _Yp = Nout;

            /* Copy as many samples as possible to the output buffer */
            if (_Yp && (outBufferLen - outSampleCount) > 0) {
                len = std::min(outBufferLen - outSampleCount, _Yp);
                for (size_t i = 0; i < len; i++)
                    outBuffer[outSampleCount + i] = (short)_Y[i];
                outSampleCount += len;
                for (size_t i = 0; i < _Yp - len; i++)
                    _Y[i] = _Y[i + len];
                _Yp -= len;
            }

            /* If there are still output samples left, return now,
               since we need the full output buffer available */
            if (_Yp)
                break;
        }

        return outSampleCount;
    }

private:
    enum {
      _Npc   = 4096,
      _NWING = _Npc*17 /* # of filter coeffs in right wing */
    };

    double Izero(double x)
    {
        const double IzeroEPSILON = 1E-21; /* Max error acceptable in Izero */

        double sum = 1.0;
        double u = 1.0;
        size_t n = 1;
        do {
            double temp = ( x / 2.0 ) / (double)n;
            temp *= temp;
            u *= temp;
            sum += u;
            n += 1;
        } while (u >= IzeroEPSILON*sum);
        return sum;
    }

    float lrsFilterUp(float *Xp,    /* Current sample */
                      double Ph,    /* Phase */
                      int Inc)      /* increment (1 for right wing or -1 for left) */
    {
        Ph *= _Npc;                 /* Npc is number of values per 1/delta in impulse response */

        float v = 0.0; /* The output value */
        float *Hp = &_Imp[(int)Ph];
        float *Hdp = &_ImpD[(int)Ph];
        float *End = &_Imp[_NWING];
        double a = Ph - std::floor(Ph); /* fractional part of Phase */

        if (Inc == 1)           /* If doing right wing...              */
        {                       /* ...drop extra coeff, so when Ph is  */
            End--;              /*    0.5, we don't do too many mult's */
            if (Ph == 0)        /* If the phase is zero...             */
            {                   /* ...then we've already skipped the   */
                Hp += _Npc;     /*    first sample, so we must also    */
                Hdp += _Npc;    /*    skip ahead in _Imp[] and _ImpD[] */
            }
        }

        while (Hp < End) {
            float t = *Hp;      /* Get filter coeff */
            t += (*Hdp)*a;      /* t is now interp'd filter coeff */
            Hdp += _Npc;        /* Filter coeff differences step */
            t *= *Xp;           /* Mult coeff by input sample */
            v += t;             /* The filter output */
            Hp += _Npc;         /* Filter coeff step */
            Xp += Inc;          /* Input signal step. NO CHECK ON BOUNDS */
        }

        return v;
    }

    float lrsFilterUD(float *Xp,  /* Current sample */
                      double Ph,  /* Phase */
                      int Inc,    /* increment (1 for right wing or -1 for left) */
                      double dhb) /* filter sampling period */
    {

        float v = 0.0; /* The output value */
        double Ho = Ph * dhb;
        float *End = &_Imp[_NWING];

        if (Inc == 1)           /* If doing right wing...              */
        {                       /* ...drop extra coeff, so when Ph is  */
            End--;              /*    0.5, we don't do too many mult's */
            if (Ph == 0)        /* If the phase is zero...             */
                Ho += dhb;      /* ...then we've already skipped the   */
        }                       /*    first sample, so we must also    */
        /*    skip ahead in _Imp[] and _ImpD[] */

        float *Hp;
        while ((Hp = &_Imp[(int)Ho]) < End) {
            float t = *Hp;                  /* Get IR sample */
            float *Hdp = &_ImpD[(int)Ho];   /* get interp bits from diff table*/
            float a = Ho - floor(Ho);       /* a is logically between 0 and 1 */
            t += (*Hdp)*a;                  /* t is now interp'd filter coeff */
            t *= *Xp;                       /* Mult coeff by input sample */
            v += t;                         /* The filter output */
            Ho += dhb;                      /* IR step */
            Xp += Inc;                      /* Input signal step. NO CHECK ON BOUNDS */
        }

        return v;
    }

    /* Sampling rate up-conversion only subroutine;
     * Slightly faster than down-conversion;
     */
    size_t lrsSrcUp(double factor, size_t Nx, float LpScl)
    {
        double endTime = _Time + Nx;            /* When Time reaches EndTime, return to user */
        double dt = 1.0 / factor;               /* Output sampling period, Step through input signal */ 

        float *Y = _Y;
        while (_Time < endTime)
        {
            double LeftPhase = _Time - floor(_Time);
            double RightPhase = 1.0 - LeftPhase;

            float *Xp = &_X[(size_t)_Time];        /* Ptr to current input sample */
            /* Perform left-wing inner product */
            float v = lrsFilterUp(Xp, LeftPhase, -1);
            /* Perform right-wing inner product */
            v += lrsFilterUp(Xp + 1, RightPhase, 1);

            v *= LpScl;                         /* Normalize for unity filter gain */
            *Y++ = v;                           /* Deposit output */
            _Time += dt;                        /* Move to next sample by time increment */
        }

        return Y - _Y;                          /* Return the number of output samples */
    }

    /* Sampling rate conversion subroutine */
    size_t lrsSrcUD(double factor, size_t Nx, float LpScl)
    {

        double endTime = _Time + Nx;            /* When Time reaches EndTime, return to user */
        double dh = std::min<float>(_Npc, factor * _Npc); /* Filter sampling period, Step through filter impulse response */
        double dt = 1.0 / factor;               /* Output sampling period, Step through input signal */ 

        float *Y = _Y;
        while (_Time < endTime)
        {
            double LeftPhase = _Time - std::floor(_Time);
            double RightPhase = 1.0 - LeftPhase;

            float *Xp = &_X[(size_t)_Time];        /* Ptr to current input sample */
            /* Perform left-wing inner product */
            float v = lrsFilterUD(Xp, LeftPhase, -1, dh);
            /* Perform right-wing inner product */
            v += lrsFilterUD(Xp + 1, RightPhase, 1, dh);

            v *= LpScl;                         /* Normalize for unity filter gain */
            *Y++ = v;                           /* Deposit output */
            _Time += dt;                        /* Move to next sample by time increment */
        }

        return Y - _Y;                          /* Return the number of output samples */
    }

private:
    float   _Imp[_NWING];
    float   _ImpD[_NWING];
    float   _LpScl;
    double  _minFactor;
    double  _maxFactor;
    size_t  _XSize;
    float  *_X = nullptr;
    float  *_Y = nullptr;
    size_t  _Xp; /* Current "now"-sample pointer for input */
    size_t  _Xread; /* Position to put new samples */
    size_t  _Xoff;
    size_t  _YSize;
    size_t  _Yp;
    double  _Time;
};
