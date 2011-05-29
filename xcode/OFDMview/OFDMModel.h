// OFDMModel definition
// OFDMview

#import <Cocoa/Cocoa.h>
#include <fftw3.h>

@interface OFDMModel : NSObject {
	double *samples;
	int nsamples;
	fftw_complex *symbols;
	int nsymbols;
}

- (int)loadData :(char*)filename;
- (void)adjust : (int)frq;
- (int)getSampleCount;
- (void)runFFT : (int)startsamp;
- (void)constellationIter :(id)sender :(int)carrier;

@end
