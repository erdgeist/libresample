#include "resampler.h"
#include <sys/uio.h>
#include <unistd.h>

/*
int main() {
    short inbuffer[44100];
    short outbuffer[8000];
    double factor = ((double)sizeof(outbuffer))/((double)sizeof(inbuffer));
    Resampler resample(factor, 1.0 / factor);

    while (1) {
        size_t in = read(0, inbuffer, sizeof(inbuffer));

        if (in <= 0)
            return 0;

        size_t inbufused;
        size_t out = resample.process(factor, inbuffer, in / sizeof(short), &inbufused, outbuffer, sizeof(outbuffer) / sizeof(short));
        write(1, outbuffer, out * sizeof(short));
        if (inbufused != sizeof(inbuffer))
            memmove(inbuffer, inbuffer + inbufused, (sizeof(inbuffer) / sizeof(short) ) - inbufused);
    }
    return 0;
}
*/


int main() {
    double inrate = 8000;
    double outrate = 48000;

    short outbuffer[960];
    short inbuffer[240];
    double factor = (outrate / inrate);
    Resampler resample( std::min(factor, 1.0 / factor ), std::max(factor, 1.0 / factor) );

    while (1) {
        size_t in = read(0, inbuffer, sizeof(inbuffer));

        if (in <= 0)
            return 0;

        size_t inbufused;
        size_t out = resample.process(factor, inbuffer, in / sizeof(short), &inbufused, outbuffer, sizeof(outbuffer) / sizeof(short));

        fprintf( stderr, "yielded: %zd, used: %zd\n", out, inbufused);

        write(1, outbuffer, out * sizeof(short));
        if (inbufused != sizeof(inbuffer))
            memmove(inbuffer, inbuffer + inbufused, (sizeof(inbuffer) / sizeof(short) ) - inbufused);
    }
    return 0;
}
